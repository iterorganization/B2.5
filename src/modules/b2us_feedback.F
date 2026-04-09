      module b2us_feedback
      use b2mod_types
      use b2mod_indirect
      use b2mod_geometry
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2plot
     & , only : b2_fnmti
      use b2mod_diag
      use b2mod_time
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
      use b2mod_ad
     & , only : nsdmax
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_user_namelist
     & , only : nomp, omp, icsepomp, nimp, imp, icsepimp,
     &          npfr_cvs, pfr_cvs
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_wneutrals
     & , only : wldna, wldra, wldpa, wldnm, wldrm, wldpm,
     &          ewldrp_res, ewldt_res
      use eirmod_eirbra
     & , only : fnmti_B2
#else
      use eirdiag          ! IGNORE
#endif
#else
      use b2mod_eirdiag
#endif
      use b2mod_dimensions
      implicit none
      integer, parameter :: nvac=4, nvacreg=4, nMxFbReg = 1000

*   ..variables
      real (kind=R8), save ::
     * cum_volrec, dt_prev, ddtim, she_rad_tot, she_eir_tot
      logical, save :: feedback_namelist_used
      logical, save :: feedback_read, control_read
      logical, save :: feedback_allocated = .false.
      real (kind=R8) :: charge_frac(DEF_NSD)

      integer, save :: icped = 0, b2sral_style = 2, ncallfeedback = 0
      integer, save :: geometryType

      integer, allocatable, save :: fbRegP(:,:) ! (nFb,2) pointing for each feedback to the first index in the fbReg list, and its number of CVs or faces
      integer, allocatable, save :: fbReg(:)  !         listing for each feedback the cells or faces
      integer, save :: nfbreg = 0, nfb = 0
*   ..namelists
! csc new namelist merging all feedback schemes
      real (kind=R8), save ::
     & fb_target(DEF_NATM) = 0.0_R8, fb_time(DEF_NATM) = 0.0_R8,
     & saved_fb_actuator(DEF_NATM) = 0.0_R8,
     & saved_fb_prev(DEF_NATM) = 0.0_R8,
     & fb_alpha(DEF_NATM) = 0.001_R8, fb_beta(DEF_NATM) = 1.0_R8,
     & fb_const(DEF_NATM) = 0.0_R8, fb_puff_min(DEF_NATM) = 0.0_R8,
     & fb_puff_max(DEF_NATM) = 0.0_R8, fb_overshoot(DEF_NATM) = 0.0_R8,
     & fb_current(DEF_NATM) = 0.0_R8, fb_rescale(DEF_NATM) = 1.0_R8,
     & fb_current_prev(DEF_NATM) = 0.0_R8

      integer, save :: fb_reg_par(DEF_NATM,2) = 0, !contains face label to define list of faces for feedback
     & fb_species(DEF_NATM) = 0, fb_ib(DEF_NATM) = -1-DEF_NSTS,
     & fb_type(DEF_NATM) = 0, fb_rescale_option(DEF_NATM) = 0,
     & fb_actuator(DEF_NATM) = 0, fb_istra(DEF_NATM) = 0,
     & fb_reg(nMxFbReg) = 0, fb_regp(DEF_NATM,2) = 0 !these two arrays have same meaning as fbreg and fbregp above but are needed for conversion and namelist
      logical, save :: fb_type_inverse(DEF_NATM) = .false.

! csc old switches
      integer, save :: fhiycore_kinetic_energy
      integer, save :: coreregno, coreregn2, solregno,
     &                 pfrregno1, pfrregno2,
     &                 isfeedback, nesepm_istra, iyped
      real (kind=R8), save ::
     &  fnaycore, fnaycore_alpha,
     1  fheycore, fheycore_alpha,
     1  fhiycore, fhiycore_alpha,
     2  fchycore, fchycore_alpha,
     3  nesepm, nesepm_sol, nepedm_sol, nesepm_pfr,
     3  nesepm_alpha, nesepm_beta, nesepm_gamma,
     3  nesepm_minpuff, nesepm_maxpuff, nesepm_overshoot,
     4  volrec_sol, volrec_alpha, volrec_beta, volrec_overshoot,
     4  ndes, ndes_sol, private_flux_puff,
     9  save_core_flux_while_overshot, save_gas_puff_while_overshot,
     9  save_pfr1_puff_while_overshot, save_pfr2_puff_while_overshot
      data
     1  fnaycore/-1.0e30_R8/, fnaycore_alpha/1.0e-1_R8/,
     1  fheycore/0.0_R8/, fheycore_alpha/1.0e-1_R8/,
     2  fhiycore/0.0_R8/, fhiycore_alpha/1.0e-1_R8/,
     2  fhiycore_kinetic_energy/0/,
     2  fchycore/-1.0e30_R8/, fchycore_alpha/1.0e-1_R8/,
     3  nesepm/0.0_R8/, nesepm_sol/0.0_R8/, nepedm_sol/0.0_R8/,
     3  nesepm_pfr/0.0_R8/, nesepm_alpha/1.0e-3_R8/,
     3  nesepm_beta/0.0_R8/, nesepm_gamma/0.99_R8/,
     3  nesepm_minpuff/0.0_R8/, nesepm_maxpuff/0.0_R8/,
     3  nesepm_overshoot/0.0_R8/,
     4  volrec_sol/0.0_R8/, volrec_alpha/1.0e-3_R8/,
     4  volrec_beta/1.0_R8/, volrec_overshoot/0.0_R8/,
     5  ndes/0.0_R8/, ndes_sol/0.0_R8/,
     7  private_flux_puff/0.0_R8/,
     9  save_core_flux_while_overshot/0.0_R8/, !WG_RM not used anymore with new structure
     9  save_gas_puff_while_overshot/0.0_R8/,  !WG_RM not used anymore with new structure
     9  save_pfr1_puff_while_overshot/0.0_R8/, !WG_RM not used anymore with new structure
     9  save_pfr2_puff_while_overshot/0.0_R8/, !WG_RM not used anymore with new structure
     i  coreregno/1/, coreregn2/0/, solregno/6/,
     i  pfrregno1/4/, pfrregno2/5/, isfeedback/0/, nesepm_istra/-1/

      namelist /feedback/ saved_fb_actuator, saved_fb_prev
! new namelist for specifying more complicated feedback options
      namelist /feedback_control/
! csc new namelist merging all feedback schemes
     & nfb,
     & fb_target, fb_type, fb_rescale_option,
     & fb_actuator, fb_alpha, fb_beta,
     & fb_const, fb_time,
     & fb_ib, fb_species, fb_type_inverse,
     & fb_puff_min, fb_puff_max,
     & fb_overshoot, fb_reg, fb_regp,
     & fb_reg_par, fb_istra
      data
     & feedback_namelist_used/.false./, feedback_read/.false./,
     & control_read/.false./

      character*72 message

      contains

      subroutine read_b2us_feedback(nread,edition,found)
      implicit none
      integer nread
      character*4 edition
      logical found
      character*256 filename
      logical file_ok

      if (edition.eq.'    ') then
        filename = 'b2.feedback_save.parameters'
      else
        filename = 'b2.feedback_save.parameters.'//edition
      endif
      call find_file(filename,file_ok)
      if(file_ok.and.(.not.feedback_read.or.edition.ne.'    ')) then
        open(nread,file=trim(filename))
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,feedback,err=999)
#else
        read(nread,feedback)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
        if (edition.eq.'    ') feedback_read = .true.
      else
        if (.not.file_ok) write(6,*) 'No ',trim(filename)
      endif
      found = maxval(saved_fb_actuator).gt.0.0_R8 .or.
     &        maxval(saved_fb_prev).gt.0.0_R8
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading FEEDBACK namelist!')
#endif
      end subroutine read_b2us_feedback

      subroutine write_b2us_feedback(nwrite,edition)
      implicit none
      integer nwrite
      character*4 edition
      integer is
      character*26 hlp_frm
      character*256 filename

c      feedback_namelist_used = feedback_namelist_used .or.
c     .            maxval(saved_fb_actuator).gt.0.0_R8 .or.
c     .            maxval(saved_fb_prev).gt.0.0_R8

      if (edition.eq.'    ') then
        filename = 'b2.feedback_save.parameters'
        feedback_read = .false.
      else
        filename = 'b2.feedback_save.parameters.'//edition
      endif
cxpb We do the write explicitly to avoid array overflows when the problem
cxpb size is changed in b2mod_dimensions
      if(feedback_namelist_used) then
        open(nwrite,file=trim(filename))
        write(nwrite,'(a9)') '&FEEDBACK'
        write(hlp_frm,'(a,i2,a)')
     >   '(a,1p,',nFb,'(1x,1e15.8,'',''))'
        write(nwrite,hlp_frm) 'SAVED_FB_ACTUATOR=',
     >   (saved_fb_actuator(is),is=1,nFb)
        write(nwrite,hlp_frm) 'SAVED_FB_PREV=',
     >   (saved_fb_prev(is),is=1,nFb)
        write(nwrite,'(a1)') '/'
        close(nwrite)
        write(*,*) 'Saved ',trim(filename)
      endif
      return

      end subroutine write_b2us_feedback

      subroutine read_b2us_feedback_control(nread)
      implicit none
      integer nread
      character*256 filename
      logical file_ok

      if (control_read) return
      filename='b2.feedback_control.parameters'
      call find_file(filename,file_ok)
      if(file_ok) then
        open(nread,file=filename)
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,feedback_control,err=999)
#else
        read(nread,feedback_control)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
        control_read = .true.
      else
        write(6,*) 'No ',trim(filename)
      endif
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading FEEDBACK_CONTROL namelist!')
#endif
      end subroutine read_b2us_feedback_control

      subroutine write_b2us_feedback_control(nwrite)
      implicit none
      integer nwrite
      open(nwrite,file='b2.feedback_control.parameters')
      write(nwrite,feedback_control)
      close(nwrite)
      control_read = .false.
      end subroutine write_b2us_feedback_control

      subroutine init_feedback(nFc, ns, switch, mpg, geo)
      implicit none
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      integer :: nFc, ns
      integer ifb, ifbb, istrai, icount
      integer ifbreg, is, ib, is0, is00, is_start, is_end, iatm
      logical consistent, done

      call subini('init_feedback')
      call read_b2us_feedback_control(99)
      call read_b2us_feedback(99,'    ',feedback_namelist_used)
      call xertst(icsepomp.gt.0,
     & 'Invalid icsepomp value, check rzomp in b2.user.parameters')
      if (any(fb_type.eq.9)) call xertst(icsepimp.gt.0,
     & 'Invalid icsepimp value, check rzimp in b2.user.parameters')
      call ipgetr('b2mndr_dtim', ddtim)
      call xertst(ddtim.gt.0.0_R8, 'faulty parameter dtim')
      dt_prev = ddtim
      call ipgeti ('b2sral_style', b2sral_style)
*  ..first read data for other feedback schemes given directly in b2mn.dat
*    This data will be included in the new structure
      ifb = nFb
      call ipgetr ('b2stbc_fnaycore', fnaycore)
      call ipgetr ('b2stbc_fnaycore_alpha', fnaycore_alpha)
      call xertst (0.0_R8.le.fnaycore_alpha,
     1  'faulty internal parameter fnaycore_alpha')
      call ipgetr ('b2stbc_fheycore', fheycore)
      call xertst (0.0_R8.le.fheycore,
     1  'faulty internal parameter fheycore')
      call ipgetr ('b2stbc_fheycore_alpha', fheycore_alpha)
      call xertst (0.0_R8.le.fheycore_alpha,
     1  'faulty internal parameter fheycore_alpha')
      call ipgetr ('b2stbc_fhiycore', fhiycore)
      call xertst (0.0_R8.le.fhiycore,
     1  'faulty internal parameter fhiycore')
      call ipgetr ('b2stbc_fhiycore_alpha', fhiycore_alpha)
      call xertst (0.0_R8.le.fhiycore_alpha,
     1  'faulty internal parameter fhiycore_alpha')
      call ipgeti ('b2tsbc_fhiycore_kinetic_energy',
     1 fhiycore_kinetic_energy)
      call ipgetr ('b2stbc_fchycore', fchycore)
      call ipgetr ('b2stbc_fchycore_alpha', fchycore_alpha)
      call xertst (0.0_R8.le.fchycore_alpha,
     1  'faulty internal parameter fchycore_alpha')
      call ipgetr ('b2stbc_nesepm', nesepm)
      call xertst (0.0_R8.le.nesepm,
     1  'faulty internal parameter nesepm')
      call ipgetr ('b2stbc_nesepm_sol', nesepm_sol)
      call xertst (0.0_R8.le.nesepm_sol,
     1  'faulty internal parameter nesepm_sol')
      call ipgetr ('b2stbc_nepedm_sol', nepedm_sol)
      call xertst (0.0_R8.le.nepedm_sol,
     1  'faulty internal parameter nepedm_sol')
      call ipgetr ('b2stbc_nesepm_pfr', nesepm_pfr)
      call xertst (0.0_R8.le.nesepm_pfr,
     1  'faulty internal parameter nesepm_pfr')
      call ipgetr ('b2stbc_volrec_sol', volrec_sol)
      call xertst (0.0_R8.le.volrec_sol,
     1  'faulty internal parameter volrec_sol')
      if (volrec_sol.gt.0.0_R8) call xertst(switch%use_eirene.ne.0,
     1 'b2stbc_volrec_sol needs implementation of cum_volrec for'//
     2 ' fluid neutrals')
      call ipgetr ('b2stbc_nesepm_alpha', nesepm_alpha)
      call xertst (0.0_R8.le.nesepm_alpha,
     1  'faulty internal parameter nesepm_alpha')
      call ipgetr ('b2stbc_nesepm_beta', nesepm_beta)
      call ipgetr ('b2stbc_nesepm_gamma', nesepm_gamma) !WG_RM not used anymore with new structure
      call xertst (0.0_R8.le.nesepm_gamma,
     1  'faulty internal parameter nesepm_gamma')
      call ipgetr ('b2stbc_nesepm_minpuff', nesepm_minpuff)
      call xertst (0.0_R8.le.nesepm_minpuff,
     1  'faulty internal parameter nesepm_minpuff')
      call ipgetr ('b2stbc_nesepm_maxpuff', nesepm_maxpuff)
      call xertst (0.0_R8.le.nesepm_maxpuff,
     1  'faulty internal parameter nesepm_maxpuff')
      call ipgetr ('b2stbc_nesepm_overshoot', nesepm_overshoot)
      call xertst (1.0_R8.le.nesepm_overshoot.or.
     1 0.0_R8.eq.nesepm_overshoot,
     1 'faulty internal parameter nesepm_overshoot')
      call ipgetr ('b2stbc_volrec_alpha', volrec_alpha)
      call xertst (0.0_R8.le.volrec_alpha,
     1  'faulty internal parameter volrec_alpha')
      call ipgetr ('b2stbc_volrec_beta', volrec_beta)
      call ipgetr ('b2stbc_volrec_overshoot', volrec_overshoot)
      call xertst (1.0_R8.le.volrec_overshoot.or.
     1 0.0_R8.eq.volrec_overshoot,
     1 'faulty internal parameter volrec_overshoot')
      if (mpg%nnreg(0).le.2) then
        coreregno = 0
      else
        coreregno = 1
      end if
      call ipgeti ('b2stbc_coreregno', coreregno)
      call ipgeti ('b2stbc_coreregn2', coreregn2)
      call ipgeti ('b2stbc_solregno', solregno)
      call ipgetr ('b2stbc_ndes_sol', ndes_sol)
      call xertst (0.0_R8.le.ndes_sol,
     . 'faulty internal parameter ndes_sol')
#ifdef WG_RM
      if(fnaycore.ne.-1.0e30_R8.or.fheycore.gt.0.0_R8.or.
     &   fhiycore.gt.0.0_R8.or.fchycore.ne.-1.0e30_R8.or.
     2   nesepm.gt.0.0_R8) then
        call xertst (cbirso.le.coreregno.and.
     1         coreregno.lt.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregno')
        if (mpg%nncut.ge.2) then
         call xerrab('b2stbc_nesepm needs revision for DN')
         call xertst (cbirso.le.coreregn2.and.
     1   coreregn2.lt.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregn2')
        endif
      endif
#endif
      if(nesepm_sol.gt.0.0_R8.or.volrec_sol.gt.0.0_R8.or.
     1   nepedm_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8) then
#ifdef WG_RM
        call xertst (cbirno.le.solregno.and.
     1   solregno.le.cbirno+cbnrno-1,
     2   'faulty internal parameter solregno')
#endif
        call xertst(
     .   (ndes_sol.gt.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .    nesepm_sol.eq.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .   (ndes_sol.eq.0.0_R8 .and. volrec_sol.gt.0.0_R8 .and.
     .    nesepm_sol.eq.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .   (ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .    nesepm_sol.gt.0.0_R8 .and. nepedm_sol.eq.0.0_R8) .or.
     .   (ndes_sol.eq.0.0_R8 .and. volrec_sol.eq.0.0_R8 .and.
     .    nesepm_sol.eq.0.0_R8 .and. nepedm_sol.gt.0.0_R8),
     .   'Simultaneous SOL feedback schemes not allowed!')
      endif
      iyped = icsepomp/2 - 2 !This is needed to default to same vale as ST code
      icped = icsepomp/2
      if(nepedm_sol.gt.0.0_R8) then
        call ipgeti ('b2stbc_iyped', iyped)
        icped = iyped+2
        call xertst (icped.ge.1.and.icped.le.nomp, 'faulty input iyped')
        call xertst (mpg%cvOnClosedSurface(omp(icped)),
     .   'location of top of pedestal outside core region!')
      endif
      call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
      call ipgetr ('b2stbc_ndes', ndes)
      call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
      call ipgetr ('b2stbc_private_flux_puff', private_flux_puff)
      call xertst (0.0_R8.le.private_flux_puff,
     1 'faulty internal parameter private_flux_puff')
      call ipgeti ('b2stbc_pfrregno1', pfrregno1)
      call ipgeti ('b2stbc_pfrregno2', pfrregno2)
      call ipgeti ('b2stbc_isfeedback', isfeedback)
      call xertst (0.le.isfeedback.and.isfeedback.lt.ns,
     1 'faulty internal parameter isfeedback')
      if (ndes.gt.0.0_R8.or.private_flux_puff.gt.0.0_R8.or.
     &    nesepm_pfr.gt.0.0_R8) then
#ifdef WG_RM
        call xertst (cbirso.le.pfrregno1.and.
     1   pfrregno1.le.cbirso+cbnrso-1,
     2   'faulty internal parameter pfrregno1')
        call xertst (cbirso.le.pfrregno2.and.
     1   pfrregno2.le.cbirso+cbnrso-1,
     2   'faulty internal parameter pfrregno2')
        if (pfrregno1.gt.cbirso) then
          ixl1 = max(0, nint(cbrbrk(pfrregno1-1)*mpg%nx))
        else
          ixl1 = 0
        endif
        ixr1 = min(mpg%nx, nint(cbrbrk(pfrregno1)*mpg%nx))-1
        if (pfrregno2.gt.cbirso) then
          ixl2 = max(0, nint(cbrbrk(pfrregno2-1)*mpg%nx))
        else
          ixl2 = 0
        endif
        ixr2 = min(mpg%nx, nint(cbrbrk(pfrregno2)*mpg%nx))-1
        call xertst (pfrregno1.eq.pfrregno2 .or.
     &   rightix(ixr1,-1).eq.ixl2 .or. rightix(ixr2,-1).eq.ixl1,
     .   'pfrregno1 and pfrregno2 are not connected!')
#endif
        call xertst(
     .   (ndes.gt.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and.
     .    nesepm_pfr.eq.0.0_R8) .or.
     .   (ndes.eq.0.0_R8 .and. private_flux_puff.gt.0.0_R8 .and.
     .    nesepm_pfr.eq.0.0_R8) .or.
     .   (ndes.eq.0.0_R8 .and. private_flux_puff.eq.0.0_R8 .and.
     .    nesepm_pfr.gt.0.0_R8),
     .   'Simultaneous PFR feedback schemes not allowed!')
      endif
      call ipgeti ('eirene_nesepm_istra', nesepm_istra)
      if(nesepm_istra.gt.0) then
        call xertst (nesepm_istra.le.nstrat,
     1    'faulty internal parameter nesepm_istra')
        call xertst (crcstra(nesepm_istra).eq.'C',
     1    'Stratum nesepm_istra is not declared as a puff stratum!')
      else
        do istrai = 1, nstrat
          if(nesepm_istra.le.0.and.CRCSTRA(istrai).EQ.'C') then
            nesepm_istra=istrai
          endif
        enddo
        if ((ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &   private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *   volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &   nepedm_sol.gt.0.0_R8).and.switch%use_eirene.eq.1) then
          icount = 0
          do istrai = 1, nstrat
            if (crcstra(istrai).eq.'C' .and.
     .          species_start(istrai).le.isfeedback .and.
     .          species_end(istrai).ge.isfeedback) icount = icount+1
          end do
          if (icount.eq.0) then
            message =
     .       'No puff stratum found for feedback scheme!'
          else if (icount.gt.1) then
            message =
     .       'Multiple definitions of puff strata, use nesepm_istra'
          end if
          call xertst (icount.eq.1, trim(message))
        endif
      endif
      if (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &   private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *   volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &   nepedm_sol.gt.0.0_R8) then
        call xertst (nesepm_istra.gt.0,
     &    'No "C" stratum for feedback scheme specified !')
      endif

! now add the feedbacks turned on with switches into the new format
! the idea should be that any switch overrides what is in b2.feedback_control.parameters
      if (fnaycore.ne.-1.0e30_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_fnaycore not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.10) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 10
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = fnaycore
        fb_ib(ifbb) = coreregno
        fb_actuator(ifbb) = 3
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = fnaycore_alpha
        fb_istra(ifbb) = nesepm_istra
      endif
      if (fheycore.ne.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_fheycore not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.11) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 11
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = fheycore
        fb_ib(ifbb) = coreregno
        fb_actuator(ifbb) = 4
        fb_alpha(ifbb) = fheycore_alpha
      endif
      if (fhiycore.ne.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_fchycore not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.12) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 12
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = fhiycore
        fb_ib(ifbb) = coreregno
        fb_actuator(ifbb) = 5
        fb_alpha(ifbb) = fhiycore_alpha
      endif
      if (fchycore.ne.-1.0e30_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_fchycore not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.13) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 13
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = fchycore
        fb_ib(ifbb) = coreregno
        fb_actuator(ifbb) = 6
        fb_alpha(ifbb) = fchycore_alpha
      endif
      if (ndes.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_ndes not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.14) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 14
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = ndes
        fb_ib(ifbb) = pfrregno1
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_istra(ifbb) = nesepm_istra
      endif
      if (nesepm.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_nesepm not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 3
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = nesepm
        fb_ib(ifbb) = coreregno
        fb_actuator(ifbb) = 3
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = nesepm_alpha
        fb_overshoot(ifbb) = nesepm_overshoot
      endif
      if (nesepm_sol.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_nesepm_sol not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 3
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = nesepm_sol
        fb_ib(ifbb) = solregno
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = nesepm_alpha
        fb_beta(ifbb) = nesepm_beta
        fb_puff_min(ifbb) = nesepm_minpuff
        fb_puff_max(ifbb) = nesepm_maxpuff
        fb_overshoot(ifbb) = nesepm_overshoot
        fb_istra(ifbb) = nesepm_istra
      endif
      if (ndes_sol.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_ndes_sol not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.14) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 14
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = ndes_sol
        fb_ib(ifbb) = solregno
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = nesepm_alpha
        fb_puff_min(ifbb) = nesepm_minpuff
        fb_puff_max(ifbb) = nesepm_maxpuff
        fb_overshoot(ifbb) = nesepm_overshoot
        fb_istra(ifbb) = nesepm_istra
      endif
      if (nesepm_pfr.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_nesepm_pfr not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 3
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = nesepm_pfr
        fb_ib(ifbb) = pfrregno1
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = nesepm_alpha
        fb_puff_min(ifbb) = nesepm_minpuff
        fb_puff_max(ifbb) = nesepm_maxpuff
        fb_overshoot(ifbb) = nesepm_overshoot
        fb_istra(ifbb) = nesepm_istra
      endif
      if (private_flux_puff.gt.0.0_R8) then
        call xerrab('private_flux_puff not converted for US code')
        ! not sure what to do with this..
        write(*,*) 'WARNING! Feedback b2stbc_private_flux_puff '//
     &  'not fully backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.13) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 13
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = private_flux_puff
        fb_ib(ifbb) = pfrregno1
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = nesepm_alpha
        fb_puff_min(ifbb) = nesepm_minpuff
        fb_puff_max(ifbb) = nesepm_maxpuff
        fb_overshoot(ifbb) = nesepm_overshoot
        fb_istra(ifbb) = nesepm_istra
      endif
      if (volrec_sol.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_volrec_sol '//
     &  'not fully backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.15) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 15
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = volrec_sol
        fb_ib(ifbb) = solregno
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_alpha(ifbb) = volrec_alpha
        fb_beta(ifbb) = volrec_beta
        fb_puff_min(ifbb) = nesepm_minpuff
        fb_puff_max(ifbb) = nesepm_maxpuff
        fb_overshoot(ifbb) = volrec_overshoot
        fb_istra(ifbb) = nesepm_istra
      endif
      if (nepedm_sol.gt.0.0_R8) then
        write(*,*) 'WARNING! Feedback b2stbc_nepedm_sol not fully '//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.16) ifbb=ifb !to override what is present
        enddo
        if (ifbb.eq.0) then
          nFb = nFb + 1
          ifbb = nFb
        endif
        call xertst(nFb.le.DEF_NATM,
     &   'init_feedback: increase size of DEF_NATM')
        fb_type(ifbb) = 16
        fb_rescale_option(ifbb) = 1
        fb_target(ifbb) = nepedm_sol
        fb_ib(ifbb) = solregno
        fb_actuator(ifbb) = 1
        fb_species(ifbb) = isfeedback
        fb_istra(ifbb) = nesepm_istra
      endif

      call xertst(nFb.le.DEF_NATM,
     & 'b2us_feedback: increase size of DEF_NATM')

      do ifb = 1, nFb
        if(fb_actuator(ifb).ne.1.or.switch%use_eirene.eq.0) cycle
        if(fb_istra(ifb).gt.0) then
          call xertst (fb_istra(ifb).le.nstrat,
     1     'faulty internal parameter fb_istra')
          call xertst (crcstra(fb_istra(ifb)).eq.'C',
     1     'Stratum fb_istra is not declared as a puff stratum!')
        else
          do istrai = 1, nstrat
            if(fb_istra(ifb).le.0.and.CRCSTRA(istrai).EQ.'C') then
              fb_istra(ifb)=istrai
            endif
          enddo
          icount = 0
          do istrai = 1, nstrat
            if (crcstra(istrai).eq.'C' .and.
     .          species_start(istrai).le.fb_species(ifb) .and.
     .          species_end(istrai).ge.fb_species(ifb))
     .          icount = icount+1
          end do
          if (icount.eq.0) then
            message =
     .       'No puff stratum found for feedback scheme!'
          else if (icount.gt.1) then
            message =
     .       'Multiple definitions of puff strata, use nesepm_istra'
          end if
          call xertst (icount.eq.1, trim(message))
        endif
        call xertst (fb_istra(ifb).gt.0,
     &    'No "C" stratum for feedback scheme specified !')
      enddo
! sanity checks of namelist and converted input
      do iFb = 1, nFb
        call xertst(0.le.fb_species(iFb).and.fb_species(iFb).lt.ns,
     .   'faulty parameter fb_species!')
        call xertst(0.le.fb_type(iFb), 'faulty parameter fb_type!')
        call xertst(0.le.fb_rescale_option(iFb),
     .   'faulty parameter fb_rescale_option!')
        call xertst(0.le.fb_actuator(iFb),
     .   'faulty parameter fb_actuator!')
        call xertst( (fb_actuator(iFb).ne.8 .and.
     .   ((0.lt.fb_ib(iFb).and.fb_ib(iFb).le.nbc).or.
     .    fb_ib(iFb).eq.(-1-DEF_NSTS))) .or.
     .   (fb_actuator(iFb).eq.8.and.
     .    ((abs(fb_ib(iFb)).le.nsts.and.switch%use_eirene.ne.0).or.
     .     (0.le.fb_ib(iFb).and.fb_ib(iFb).le.nbc.and.
     .      switch%use_eirene.eq.0))),
     .   'faulty parameter fb_ib!')
        call xertst(0.le.fb_istra(iFb).and.fb_istra(iFb).le.nstrai,
     .   'faulty parameter fb_istra!')
        call xertst(0.0_R8.le.fb_puff_min(iFb).and.
     .    fb_puff_min(iFb).le.fb_puff_max(iFb),
     .   'faulty parameters fb_puff_min/max!')
        call xertst(0.0_R8.le.fb_overshoot(iFb),
     .   'faulty parameter fb_overshoot!')
        call xertst(0.le.fb_regp(iFb,1).and.0.le.fb_regp(iFb,2),
     .   'faulty parameters fb_regp!')
        call xertst(0.le.fb_reg_par(iFb,1).and.0.le.fb_reg_par(iFb,2),
     .   'faulty parameters fb_reg_par!')
      end do

      consistent = .true.
      do ifb = 1, nFb
        if (((fb_ib(ifb).le.0.and.switch%use_eirene.eq.0) .or.   ! not using EIRENE, but no fb_ib
     &       (fb_istra(ifb).le.0 .and. switch%use_eirene.ne.0    ! using EIRENE, but no puff stratum
     &        .and.fb_actuator(ifb).eq.1)) .and.
     &      .not.fb_actuator(ifb).eq.2 .and.                     ! feedback not 2 => either puff stratum or fb_ib needed
     &      .not.fb_actuator(ifb).eq.8 )                         ! feedback = 8 values already tested above
     &    call xerrab
     .     ('Init_feedback: fb_ib <= 0 and actuator != (2 or 8)')
        call xertst(0.le.fb_species(ifb).and.fb_species(ifb).le.ns-1,
     &    'Init_feedback: fb_species>ns!')
        iatm = b2espcr(fb_species(ifb))
        is0 = eb2spcr(iatm)
        do ifbb = ifb+1, nFb
          is00 = eb2spcr(b2espcr(fb_species(ifbb)))
          call xertst(.not.(is0.eq.is00),'init_feedback: cannot '//
     &     'define multiple feedbacks on same atom sequence')
          if ( (fb_type(ifb).eq.fb_type(ifbb)) .and. .not.
     &         ((fb_type(ifb).eq.1).or.(fb_type(ifb).eq.4).or.
     &          (fb_type(ifb).eq.5).or.(fb_type(ifb).eq.7).or.
     &          (fb_type(ifb).eq.8).or.(fb_type(ifb).eq.17).or.
     &          (fb_type(ifb).eq.22).or.(fb_type(ifb).eq.27)) )
     &          call xerrab('init_feedback:'//
     &      ' set of feedbacks cannot be applied multiple times')
        enddo
        is_start = is0+1
        is_end = is0+nfluids(iatm)
        if (fb_actuator(ifb).ne.2) then
          consistent = consistent .and.
     &     ((fb_target(ifb) .eq. 0.0_R8 .or.
     &       fb_type(ifb) .eq. 0 .or.
     &       fb_rescale_option(ifb) .eq. 0 .or.
     &       fb_actuator(ifb) .eq. 0) .or.
! csc For fluid neutrals (or fluid neutral strata, TBD) change directly the BC
!     For kinetic neutrals (or kinetic strata, TBD), change directly userfluxparm, no issue with fluid BC
!     Can be arbitrarily expanded to allow applying feedback to any already existent BC
     &      (fb_actuator(ifb) .eq. 1 .and. !neutral flux (or gas puff)
     &       ( bccon(is0,fb_ib(ifb)) .eq. 8 .or.
     &        switch%use_eirene.ne.0) ) .or.
     &      (fb_actuator(ifb) .eq. 3 .and. !charged particle flux ONLY for fb_type 6, to be improved?
     &       maxval(bccon(is_start:is_end,
     .                    fb_ib(ifb))).eq.8 .and.
     &       minval(bccon(is_start:is_end,
     .                    fb_ib(ifb))).eq.8 .and.
     &       fb_type(ifb).eq.6) .or.
     &      (fb_actuator(ifb) .eq. 3 .and. !charged particle flux, more general
     &       (bccon(fb_species(ifb),fb_ib(ifb)) .eq. 8 .or.
     &        bccon(fb_species(ifb),fb_ib(ifb)) .eq. 13 ) .and.
     &       .not.fb_type(ifb).eq.6) .or.
     &      (fb_actuator(ifb) .eq. 4 .and. !electron heat flux
     &       bcene(fb_ib(ifb)) .eq. 8) .or.
     &      (fb_actuator(ifb) .eq. 5 .and. !ion heat flux
     &       bceni(fb_ib(ifb)) .eq. 8) .or.
     &      (fb_actuator(ifb) .eq. 6 .and. !current
     &       bcpot(fb_ib(ifb)) .eq. 8) .or.
     &      (fb_actuator(ifb) .eq. 7 .and. !density
     &       (bccon(fb_species(ifb),fb_ib(ifb)) .eq. 1 .or.
     &        bccon(fb_species(ifb),fb_ib(ifb)) .eq. 28)) .or.
     &       fb_actuator(ifb) .eq. 8)
        end if
      end do
      consistent = consistent .or. b2sral_style.eq.1 .or.
     & switch%b2stbc_boundary_namelist.lt.1
      do ifb = 1, nFb
        if (fb_actuator(ifb).eq.1) cycle
        consistent = consistent .and.
     &    .not.(fb_rescale_option(ifb).eq.4 .or.
     &          fb_rescale_option(ifb).eq.5)
      end do
      do ifb = 1, nFb
        if (fb_actuator(ifb).ne.2) cycle
        do is = is_start, is_end
          do ib = 1, nbc
            consistent = consistent .and.
     &       (bccon(is,ib).eq.0 .or. bccon(is,ib).eq.2 .or.
     &        bccon(is,ib).eq.3 .or. bccon(is,ib).eq.4 .or.
     &        bccon(is,ib).eq.9 .or. bccon(is,ib).eq.10.or.
     &        bccon(is,ib).eq.14.or. bccon(is,ib).eq.15.or.
     &       (bccon(is,ib).eq.5 .and.conpar(is,ib,1).eq.0.0_R8).or.
     &       (bccon(is,ib).eq.8 .and.conpar(is,ib,1).eq.0.0_R8))
          end do
        end do
      end do
      call xertst (consistent,
     . 'feedback scheme not consistent with B.C.''s!')

      if (.not.feedback_allocated) then
        allocate(fbRegP(nFb,2),fbReg(nMxFbReg))
        feedback_allocated = .true.
      endif

!   ..Calculate and fill list of CVs and faces needed for each type of feedback
      geometryType = geometryId ( mpg, geo, 2 )
      ifbreg = 1
      if  (maxval(fb_type(1:nFb)).gt.0) then
        do ifb = 1, nFb
          ! first check if fbregP has already been filled
          done = .false.
          if (fb_RegP(ifb,2).gt.0) then
            fbRegP(ifb,1) = ifbreg
            fbRegP(ifb,2) = fb_RegP(ifb,2)
            fbReg(ifbreg:ifbreg+fbRegP(ifb,2)-1) =
     &        fb_Reg(fb_RegP(ifb,1):fb_RegP(ifb,1)+fb_RegP(ifb,2)-1)
            ifbreg = ifbreg + fb_RegP(ifb,2)
            done = .true.
          endif
        ! include default to outer target for certain feedback schemes
          if (fb_type(ifb).ge.23 .and. fb_type(ifb).le.26 .and.
     &        fb_reg_par(ifb,1).eq.0) then
            write(*,*) 'No fb_reg_par is found. '//
     &       'Setting fb_reg_par to 2 (presumably the outer target)'
            fb_reg_par(ifb,1) = 2
          end if
          if (fb_type(ifb).ge.23 .and. fb_type(ifb).le.26 .and.
     &        fb_reg_par(ifb,2).eq.0) fb_reg_par(ifb,2) = 2
        ! now look for feedback types requiring faces
          if ((fb_type(ifb).eq.6  .or. fb_type(ifb).eq.10 .or.
     &         fb_type(ifb).eq.11 .or. fb_type(ifb).eq.12 .or.
     &         fb_type(ifb).eq.13 .or.
     &        (fb_type(ifb).ge.23 .and. fb_type(ifb).le.26))
     &        .and. .not.done) then
            call xertst(fb_reg_par(ifb,1).gt.0,
     &        'init_feedback: At lest one segment must be specified'//
     &        ' in fb_reg_par for fb_type = 6, 7, 10 to 13, or 15')
            ! check if same region/boundary has already been done
            do ifbb = 1, ifb-1
              if (fb_reg_par(ifb,1).eq.fb_reg_par(ifbb,1) .and.
     &            fb_reg_par(ifb,2).eq.fb_reg_par(ifbb,2)) then
                fbRegP(ifb,:) = fbRegP(ifbb,:)
                done = .true.
              endif
            end do
            if (done) cycle
            fbRegP(ifb,1) = ifbreg
            call do_faces(nFc, mpg, ifb, ifbreg)
            fbRegP(ifb,2) = ifbreg - fbRegP(ifb,1)
          endif
        ! now look for feedback types requiring volumes
          if ((fb_type(ifb).eq.1 .or. fb_type(ifb).eq.2 .or.
     &        fb_type(ifb).eq.4  .or. fb_type(ifb).eq.5 .or.
     &        fb_type(ifb).eq.8) .and. .not.done) then
            call xertst(fb_reg_par(ifb,1).gt.0,
     &        'init_feedback: At lest one segment must be specified'//
     &        ' in fb_reg_par for fb_type = 1, 2, 5, or 8')
            ! check if same region/boundary has already been done
            do ifbb = 1, ifb-1
              if (fb_reg_par(ifb,1).eq.fb_reg_par(ifbb,1)
     &         .and.fb_reg_par(ifb,2).eq.fb_reg_par(ifbb,2)) then
                fbRegP(ifb,:) = fbRegP(ifbb,:)
                done = .true.
              endif
            end do
            if (done) cycle
            fbRegP(ifb,1) = ifbreg
            call do_volumes(nFc, mpg, ifb, ifbreg)
            fbRegP(ifb,2) = ifbreg - fbRegP(ifb,1)
          endif
        end do
      endif
      nfbreg = ifbreg - 1

      call subend()
      return
      end subroutine init_feedback

      subroutine compute_feedback(nCv, nFc, ns, ismain,
     &  switch, geo, mpg, pl, dv, rt, sr, srw, psnc, psnl, main_call)
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent (in) :: nCv, nFc, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2Source), intent (in) :: sr
      type (B2SourceWork), intent (in) :: srw
      type (B2PlasmaSnapshot), intent (inout) :: psnc, psnl
      logical main_call
*   ..local variables
      integer :: ifb, iFc, iFc1, iFc2, iCv, ic1, ic2, is, is_start,
     & is_end, iss, ibc, istra, i, j, k, l, ic, iFs, find_loc,
     & iatm, imol
      real (kind=R8) :: sna_pump, sna_puff, sna_core
      integer, save :: ncall = 0
      real (kind=R8) :: sum_e, total_flux, vol_feedback, sum_s,
     & t0, t1, zaf
      real (kind=R8), external :: intface_s

      call subini ('compute_feedback')
      if (ncall.eq.0) then
*   ..test nCv, nFc, ns
        call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
        call xertst (1.le.ns, 'faulty argument ns')
        call xertst (0.le.ismain.and.ismain.lt.ns,
     .   'invalid main plasma species index ismain')
        call xertst (.not.is_neutral(ismain),
     .   'invalid main plasma species ismain; must not be neutral')
      endif

      if(maxval(fb_target(1:nFb)).gt.0.0_R8 .or.
     &   maxval(fb_rescale_option(1:nFb)).gt.0) then
! we have a request for feedback control
! the idea here is to break the actions into 3 phases
! 1/ calculate something based on the present state of the plasma depending on
!    "fb_type"
! 2/ compare this with the reference value and calculate a change
! 3/ apply this in some way based on "fb_actuator"

        do ifb = 1, nFb
          if( fb_type(ifb) .gt. 0 .and.
     &        fb_rescale_option(ifb) .gt. 0 .and.
     &        fb_actuator(ifb) .gt.0) then

! phase 1:
            select case(fb_type(ifb))

            case (0) ! no action

            case (1) ! local density

              if(ncall.eq.0) then
                write(*,*)
     &           'compute_feedback_type: using ',
     &           'local_density for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              fb_current(ifb) = 0.0_R8
              vol_feedback = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  ic1 = fbRegP(ifb,1)
                  ic2 = ic1 + fbRegP(ifb,2) - 1
                  do iCv = ic1, ic2
                    fb_current(ifb) = fb_current(ifb) +
     &                 pl%na(fbReg(iCv),is) * geo%cvVol(fbReg(iCv))
                    if(is.eq.is_start) then
                      vol_feedback = vol_feedback +
     &                               geo%cvVol(fbReg(iCv))
                    endif
                  enddo
                  is_end=is
                endif
              enddo
              fb_current(ifb) = fb_current(ifb) /
     &                 vol_feedback

            case (2) ! local electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'local_electron_density for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              fb_current(ifb) = 0.0_R8
              vol_feedback = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                 if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                   if(is_start.eq.-1) is_start=is
                   is_end=is
                 endif
              enddo
              ic1 = fbRegP(ifb,1)
              ic2 = ic1 + fbRegP(ifb,2) - 1
              do iCv = ic1, ic2
                fb_current(ifb) = fb_current(ifb) +
     &            dv%ne(fbReg(iCv)) * geo%cvVol(fbReg(iCv))
                vol_feedback = vol_feedback +
     &                         geo%cvVol(fbReg(iCv))
              enddo
              fb_current(ifb) = fb_current(ifb) /
     &                 vol_feedback

            case (3) ! nesepm electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'nesepm_electron_density for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              fb_current(ifb) = (geo%cvVol(omp(icsepomp-1))*
     &          dv%ne(omp(icsepomp-1))+geo%cvVol(omp(icsepomp))*
     &          dv%ne(omp(icsepomp)))/(geo%cvVol(omp(icsepomp-1))+
     &          geo%cvVol(omp(icsepomp)))
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (4) ! total particle content for that sequence

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'total particle content for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  ic1 = fbRegP(ifb,1)
                  ic2 = ic1 + fbRegP(ifb,2) - 1
                  do iCv = ic1, ic2
                    if(mpg%cvReg(fbReg(iCv)).ne.0 .and.
     &                 fbReg(iCv).le.mpg%nCi) then
                      fb_current(ifb) = fb_current(ifb) +
     &                  pl%na(fbReg(iCv),is) * geo%cvVol(fbReg(iCv))
                    endif
                  enddo
                  is_end=is
                endif
              enddo

            case (5) ! total ion content for that sequence (not including neutrals)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'total ion content for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  ic1 = fbRegP(ifb,1)
                  ic2 = ic1 + fbRegP(ifb,2) - 1
                  do iCv = ic1, ic2
                    if(mpg%cvReg(fbReg(iCv)).ne.0 .and.
     &                 fbReg(iCv).le.mpg%nCi) then
                      fb_current(ifb) = fb_current(ifb) +
     &                  pl%na(fbReg(iCv),is) * geo%cvVol(fbReg(iCv))
                    endif
                  enddo
                  is_end=is
                endif
              enddo

            case (6) ! neutral particle flux through the core boundary (??)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'core boundary flux for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = fb_const(ifb)
              write(*,*) 'WARNING: fb_type = 6 needs to be'//
     &         ' adapted for hybrid methods!'

              call xertst(switch%use_eirene.ne.0,
     &         'fb_type=6 not yet adapted for pure fluid cases')
              do k=1,ncns !{
                l=nlim+lcns(k)
                do j=1,nnatmi !{
                  if (latmscl(j).eq.b2espcr(fb_species(ifb))) then !{
                    fb_current(ifb)=fb_current(ifb)
     .                +wldna(l,j,0)-wldra(l,j,0)-wldpa(l,j,0)
                    do i=1,nnmoli !{
                      fb_current(ifb)=fb_current(ifb)+mlcmp(j,i)*
     .                 (wldnm(l,i,0)-wldrm(l,i,0)-wldpm(l,i,0))
                    end do !}
                  end if !}
                end do !}
              end do !}
              is_start = eb2spcr(b2espcr(fb_species(ifb)))
              is_end = eb2spcr(b2espcr(fb_species(ifb)))+
     &          nfluids(b2espcr(fb_species(ifb)))
              charge_frac = 0.0_R8
              sum_s = 0.0_R8
              do is = is_start, is_end
                iss = is-is_start+1
                if (is_neutral(is).and.switch%use_eirene.ne.0) cycle
                ic1 = fbRegP(ifb,1)
                ic2 = ic1 + fbRegP(ifb,2) - 1
                do iCv = ic1, ic2
                  if (mpg%cvOnClosedSurface(iCv)) then
                    sum_s = sum_s + pl%na(iCv,is) * geo%cvVol(iCv)
                    charge_frac(iss) = charge_frac(iss) +
     .                pl%na(iCv,is) * geo%cvVol(iCv)
                  end if
                end do
              end do
              charge_frac = charge_frac / sum_s

            case (7) ! relative average concentration of this species at the separatrix

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'average separatrix concentration for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = 0.0_R8
              sum_e = 0.0_R8
              is_start = -1
              is_end = -1
              do is = 0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is).and.switch%use_eirene.ne.0) cycle
                  iFc1 = mpg%fsFcP(mpg%iFssep,1)
                  iFc2 = iFc1 + mpg%fsFcP(mpg%iFssep,2) - 1
! csc Here the idea is to follow the separatrix along its flux surface, then use the two CVs
!     belonging to each face of the flux surface
                  do iFc = iFc1, iFc2
                    ic1 = mpg%fcCv(mpg%fsFc(iFc),1)
                    ic2 = mpg%fcCv(mpg%fsFc(iFc),2)
                    if(mpg%cvOnClosedSurface(ic1) .or.
     &                 mpg%cvOnClosedSurface(ic2)) then
                      fb_current(ifb) = fb_current(ifb) +
     &                  (pl%na(ic1,is) + pl%na(ic2,is)) *
     &                  geo%fcS(mpg%fsFc(iFc))
                      if (is.eq.is_start+1) then
                        sum_e = sum_e + (dv%ne(ic1) +
     &                    dv%ne(ic2)) * geo%fcS(mpg%fsFc(iFc))
                      end if
                    end if
                  end do
                  if (mpg%iFssep2.gt.0) then
                    iFc1 = mpg%fsFcP(mpg%iFssep2,1)
                    iFc2 = iFc1 + mpg%fsFcP(mpg%iFssep2,2) - 1
                    do iFc = iFc1, iFc2
                      ic1 = mpg%fcCv(mpg%fsFc(iFc),1)
                      ic2 = mpg%fcCv(mpg%fsFc(iFc),2)
                      if(mpg%cvOnClosedSurface(ic1) .or.
     &                   mpg%cvOnClosedSurface(ic2)) then
                        fb_current(ifb) = fb_current(ifb) +
     &                    (pl%na(ic1,is) + pl%na(ic2,is)) *
     &                    geo%fcS(mpg%fsFc(iFc))
                        if (is.eq.is_start+1) then
                          sum_e = sum_e + (dv%ne(ic1) +
     &                     dv%ne(ic2)) * geo%fcS(mpg%fsFc(iFc))
                        end if
                      end if
                    end do
                  end if
                  is_end = is
                end if
              end do
              fb_current(ifb) = fb_current(ifb)/sum_e

            case (8) ! relative average concentration of this species on desired domain

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using average ',
     &            'concentration for sequence ',
     &            b2espcr(fb_species(ifb))
              end if
              fb_current(ifb) = 0.0_R8
              sum_e = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is).and.switch%use_eirene.ne.0) cycle
                  ic1 = fbRegP(ifb,1)
                  ic2 = ic1 + fbRegP(ifb,2) - 1
                  do iCv = ic1, ic2
                    fb_current(ifb) = fb_current(ifb) +
     &                pl%na(fbReg(iCv),is) * geo%cvVol(fbReg(iCv))
                    if (is.eq.is_start+1) then
                      sum_e = sum_e +
     &                  dv%ne(fbReg(iCv)) * geo%cvVol(fbReg(iCv))
                    end if
                  enddo
                  is_end=is
                endif
              enddo
              fb_current(ifb) = fb_current(ifb)/sum_e

            case (9) ! inner midplane electron density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'inner_midplance_electron_density for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              fb_current(ifb) = 0.5_R8 *
     *          (dv%ne(imp(icsepimp-1))+dv%ne(imp(icsepimp)))
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (10) ! (core) boundary radial particle flux

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'radial particle flux for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = 0.0_R8
              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fna(fbReg(iFc),1,fb_species(ifb))
              end do
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (11) ! (core) boundary radial electron heat flux

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &            'compute_feedback_type: using ',
     &            'electron heat flux'
              endif
              fb_current(ifb) = 0.0_R8
              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fhe(fbReg(iFc),1)
              end do
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (12) ! (core) boundary radial ion heat flux

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &            'compute_feedback_type: using ',
     &            'ion heat flux '
              endif
              fb_current(ifb) = 0.0_R8
              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fhi(fbReg(iFc),1)
              enddo
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (13) ! (core) boundary radial current

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &            'compute_feedback_type: using ',
     &            'current '
              endif
              fb_current(ifb) = 0.0_R8
              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fch(fbReg(iFc),1)
              enddo
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (14) ! total particle content

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'total particle content for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              is_start=-1
              is_end=-1
              fb_current(ifb) = 0.0_R8
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  do iCv = 1, mpg%nCi
                    fb_current(ifb) = fb_current(ifb) +
     &               pl%na(iCv,is)*geo%cvVol(iCv)
                  enddo
                  is_end=is
                endif
              enddo

            case (15) ! volumetric recombination

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &            'compute_feedback_type: using ',
     &            'volumetric recombination '
                if (cum_volrec .lt. 1.0e4_R8) then
                  cum_volrec = 1.0e4_R8
                write(*,'(a,a)')
     &            'Warning! compute_feedback_type: setting ',
     &            'cum_volrec=1e4 for first timestep'
                endif
              endif
              ! in principle nothing to be done here...
              fb_current(ifb) = cum_volrec
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (16) ! nepedm

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &            'compute_feedback_type: using ',
     &            'pedestal density '
              endif
              fb_current(ifb) = dv%ne(omp(icped))
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

            case (17) !  density of ions at outer midplane

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'density of ions at outer midplane for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  fb_current(ifb) = fb_current(ifb) + 0.5_R8 *
     &              (pl%na(omp(icsepomp-1),is)+pl%na(omp(icsepomp),is))
                  is_end=is
                endif
              enddo

            case (18) !  integrated electron cooling

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'integrated electron cooling'
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = she_rad_tot + she_eir_tot
              if (fb_current(ifb).eq.0.0_R8)
     &          fb_current(ifb) = fb_target(ifb)

            case (20) !  neutral pressure in PFR

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'average neutral pressure for sequence ',
     &           b2espcr(fb_species(ifb))
                write(*,*) pfr_cvs(1:npfr_cvs)
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  is_end=is
                endif
              enddo

#ifndef SOLPS4_3
              do k = 1,npfr_cvs
                iCv = pfr_cvs(k)
                if (switch%use_eirene.ne.0) then
                  do iatm=1,nnatmi
                    if (latmscl(iatm).eq.fb_species(ifb)) then
                      t0 = dab2(iCv,iatm,1)*tab2(iCv,iatm,1)
                    end if
                  end do
                  do imol=1,nnmoli
                    if (lmolscl(imol).eq.fb_species(ifb)) then
                      t0 = t0 + dmb2(iCv,imol,1)*tmb2(iCv,imol,1)
                    end if
                  end do
                else
                  if (nint(zn(fb_species(ifb))).eq.1) then
                    t0 = pl%na(iCv,fb_species(ifb))*pl%tn(iCv)
                  else
                    t0 = pl%na(iCv,fb_species(ifb))*pl%ti(iCv)
                  end if
                end if
                fb_current(ifb) = fb_current(ifb) + t0*geo%cvVol(iCv)
              end do
#endif
              fb_current(ifb) = fb_current(ifb)/
     &         sum(geo%cvVol(pfr_cvs(1:npfr_cvs)))

            case (21) !  total radiation losses

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'total radiation losses'
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = sum(srw%rqbrm)+sum(srw%rqrad)
#ifndef SOLPS4_3
              if (switch%use_eirene.ne.0)
     >         fb_current(ifb) = fb_current(ifb) -
     &          (sum(emolrad(:,:,0))+sum(eneutrad(:,:,0))
     &          +sum(eionrad(:,:,0)))
#endif

            case (22) ! average density of this species at the separatrix

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using ',
     &            'average separatrix density for sequence ',
     &            b2espcr(fb_species(ifb))
              endif
              fb_current(ifb) = 0.0_R8
              sum_e = 0.0_R8
              is_start = -1
              is_end = -1
              do is = 0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is).and.switch%use_eirene.ne.0) cycle
! csc Here the idea is to follow the separatrix along its flux surface, then use the two CVs
!     belonging to each face of the flux surface
                  do ic = 1, 2
                    if (ic.eq.1) then
                      iFs = mpg%iFssep
                    elseif (ic.eq.2.and.
     &               geometryType.eq.GEOMETRY_CDN) then
                      iFs = mpg%iFssep2
                    else
                      cycle
                    end if
                    do k = mpg%fsFcP(iFs,1),
     &                     mpg%fsFcP(iFs,1)+mpg%fsFcP(iFs,2)-1
                      iFc = mpg%fsFc(k)
                      ic1 = mpg%fcCv(iFc,1)
                      ic2 = mpg%fcCv(iFc,2)
                      if(mpg%cvOnClosedSurface(ic1) .or.
     &                   mpg%cvOnClosedSurface(ic2)) then
                        fb_current(ifb) = fb_current(ifb) +
     &                    (pl%na(ic1,is)*geo%fcHc(iFc,2) + pl%na(ic2,is)
     &                    *geo%fcHc(iFc,1)) /
     &                    sum(geo%fcHc(iFc,:))*geo%fcS(iFc)
                        if (is.eq.is_start+1) then
                          sum_e = sum_e + geo%fcS(iFc)
                        end if
                      end if
                    enddo
                  enddo
                  is_end = is
                endif
              enddo
              fb_current(ifb) = fb_current(ifb)/sum_e

            case (23) !  Total particle heat flux at outer target

              if(ncall.eq.0) then
                if (fb_reg_par(ifb,1).eq.fb_reg_par(ifb,2)) then
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             ' particle heat flux density along boundary ',
     &              fb_reg_par(ifb,1), ' for sequence ',
     &              b2espcr(fb_species(ifb))
                else
                  write(*,'(a,a,i3,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             ' particle heat flux density along boundaries ',
     &              fb_reg_par(ifb,1),' to ',fb_reg_par(ifb,2),
     &             ' for sequence ',b2espcr(fb_species(ifb))
                end if
                fb_target(ifb) = 1.0_R8/fb_target(ifb)
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = 0.0_R8
              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                ic = mpg%fcBc(fbReg(iFc),1)
                t1 = mpg%bcFcOr(ic)
                t0 = dv%fht(fbReg(iFc),0)+dv%fht(fbReg(iFc),1)
                t0 = t1*t0/geo%fcS(fbReg(iFc))
                if (switch%use_eirene.ne.0) then
                  t0 = t0 + ewldt_res(b2_fnmti(fbReg(iFc)))
     .                    - ewldrp_res(b2_fnmti(fbReg(iFc)))
                endif
                fb_current(ifb) = max(fb_current(ifb),t0)
              enddo
              if (fb_current(ifb).ne.0.0_R8) then
                fb_current(ifb) = 1.0_R8/abs(fb_current(ifb))
              else
                fb_current(ifb) = fb_target(ifb)
              endif

            case (24) ! peak electron temperature at specific boundary

              if(ncall.eq.0) then
                if (fb_reg_par(ifb,1).eq.fb_reg_par(ifb,2)) then
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak electron temperature along boundary ',
     &              fb_reg_par(ifb,1),' for sequence ',
     &              b2espcr(fb_species(ifb))
                else
                  write(*,'(a,a,i3,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak electron temperature along boundaries ',
     &              fb_reg_par(ifb,1),' to ',fb_reg_par(ifb,2),
     &             ' for sequence ',b2espcr(fb_species(ifb))
                end if
                fb_target(ifb) = 1.0_R8/fb_target(ifb)
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = 0.0_R8

              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                iCv = min(mpg%fcCv(fbReg(iFc),1),mpg%fcCv(fbReg(iFc),2))
                fb_current(ifb) = max(fb_current(ifb),pl%te(iCv)/qe)
              enddo
              fb_current(ifb) = 1.0_R8/fb_current(ifb)

            case (25) ! peak electron density at specific boundary

              if(ncall.eq.0) then
                if (fb_reg_par(ifb,1).eq.fb_reg_par(ifb,2)) then
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak electron density along boundary ',
     &              fb_reg_par(ifb,1),' for sequence ',
     &              b2espcr(fb_species(ifb))
                else
                  write(*,'(a,a,i3,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak electron density along boundaries ',
     &              fb_reg_par(ifb,1),' to ',fb_reg_par(ifb,2),
     &             ' for sequence ',b2espcr(fb_species(ifb))
                endif
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = 0.0_R8

              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do iFc = iFc1, iFc2
                iCv = min(mpg%fcCv(fbReg(iFc),1),mpg%fcCv(fbReg(iFc),2))
                fb_current(ifb) = max(fb_current(ifb),dv%ne(iCv))
              enddo

            case (26) ! peak parallel saturation current

              if(ncall.eq.0) then
                if (fb_reg_par(ifb,1).eq.fb_reg_par(ifb,2)) then
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak parallel saturation current along boundary ',
     &             fb_reg_par(ifb,1),' for sequence ',
     &             b2espcr(fb_species(ifb))
                else
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_type: using ',
     &             'peak parallel saturation current along boundaries ',
     &             fb_reg_par(ifb,1),' to ',fb_reg_par(ifb,2),
     &             ' for sequence ',b2espcr(fb_species(ifb))
                end if
                if (fb_type_inverse(ifb))
     &           fb_target(ifb) = 1.0_R8/fb_target(ifb)
              endif

              fb_current(ifb) = 0.0_R8
              is_start=-1
              is_end=-1
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  if(is_neutral(is)) cycle
                  is_end=is
                endif
              enddo
              fb_current(ifb) = 0.0_R8

              iFc1 = fbRegP(ifb,1)
              iFc2 = iFc1 + fbRegP(ifb,2) - 1
              do k = iFc1, iFc2
                iFc = fbReg(k)
                ic = mpg%fcBc(iFc,1)
                t1 = mpg%bcFcOr(ic)
                t0 = 0.0_R8
                do is = 0, ns-1
                  if (is_neutral(is)) cycle
                  if (abs(geo%fcPbs(iFc)/geo%fcS(iFc))
     .                  .lt.geo%Qalfmin) cycle
                  zaf = intface_s(iFc,nCv,nFc,mpg%fcCv,geo%fcVol,
     .             rt%rza(1,is))
                  t0 = t0 + t1*(dv%fna(iFc,0,is)+
     .                          dv%fna(iFc,1,is))*
     .                 zaf*qe/geo%fcS(iFc)/
     .                (geo%fcQalf(iFc,0)*geo%fcBb(iFc,0)/
     .                 geo%fcBb(iFc,3))
                enddo
                fb_current(ifb) = max(fb_current(ifb),t0)
              enddo
              if (fb_type_inverse(ifb))
     &         fb_current(ifb) = 1.0_R8/fb_current(ifb)

            case (27) ! total particle content in SOL+PFR for given
                      ! isonuclear sequence (including neutrals)
                      ! test ions from EIRENE are not included for now

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_type: using ',
     &           'total particle content in SOL for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              is_start=-1
              is_end=-1
              fb_current(ifb) = 0.0_R8
              do is=0, ns-1
                if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                  if(is_start.eq.-1) is_start=is
                  do iCv = 1, mpg%nCi
                    if(mpg%cvOnClosedSurface(iCv)) cycle
                    if(mpg%cvReg(iCv).eq.0) cycle
                    fb_current(ifb) = fb_current(ifb) +
     &               pl%na(iCv,is)*geo%cvVol(iCv)
                  enddo
                  is_end=is
                endif
              enddo
#ifndef SOLPS4_3
              if (switch%use_eirene.ne.0) then
                do iatm=1,nnatmi
                  if (latmscl(iatm).eq.fb_species(ifb)) then
                    do iCv = 1, mpg%nCi
                      if(mpg%cvOnClosedSurface(iCv)) cycle
                      if(mpg%cvReg(iCv).eq.0) cycle
                      fb_current(ifb) = fb_current(ifb) +
     &                   dab2(iCv,iatm,1)*geo%cvVol(iCv)
                      do imol=1,nnmoli
                        if (mlcmp(iatm,imol).gt.0) then
                          fb_current(ifb) = fb_current(ifb) +
     &                       dmb2(iCv,imol,1)*mlcmp(iatm,imol)
     &                      *geo%cvVol(iCv)
                        endif
                      enddo
                    enddo
                  endif
                enddo
              endif
#endif

            case default
               write(*,*) 'fb_type not coded ',
     &              fb_type(ifb)
               call xerrab('fb_type not coded')

            end select

! phase 2:

            select case(fb_rescale_option(ifb))

            case (0) ! no action

            case (1) ! rescale slowed by fb_alpha

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_option: using ',
     &            'slowed_rescale for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              fb_rescale(ifb) =
     &             (1.0_R8 + fb_alpha(ifb) *
     &             fb_target(ifb) /
     &             fb_current(ifb)) /
     &             (1.0_R8 + fb_alpha(ifb))

            case (2) ! pure rescale

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_option: using ',
     &           'rescale for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              fb_rescale(ifb) = fb_target(ifb) /
     &             fb_current(ifb)

            case (3) ! rescale slowed by tanh_log

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_option: using ',
     &           'tanh_log_rescale for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

! 2**(tanh(log(x)/b)*log(a)/log(2))
              fb_rescale(ifb) = 2.0_R8 **
     &             (tanh(log(fb_target(ifb) /
     &             fb_current(ifb))/
     &             fb_beta(ifb))*
     &             log(fb_alpha(ifb))/log(2.0_R8))

            case (4) ! rescale done according to SOLPS4 formula:
                     ! N = A + B*(time-T)

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_option: using ',
     &           'SOLPS4 linear increase for sequence ',
     &           b2espcr(fb_species(ifb))
                if(fb_const(ifb).le.0.0_R8)
     &             fb_const(ifb) = fb_target(ifb)
                if(saved_fb_prev(ifb).le.0.0_R8)
     &             saved_fb_prev(ifb) = fb_current(ifb)
              endif

              fb_target(ifb) =
     .          fb_const(ifb) +
     .          fb_beta(ifb)*
     .          (tim-fb_time(ifb))
              fb_rescale(ifb) =
     .          fb_alpha(ifb) *
     .          ((fb_target(ifb) -
     .          fb_current(ifb))/ddtim +
     .          (saved_fb_prev(ifb) -
     .          fb_current(ifb))/dt_prev)

            case (5) ! rescale done according to SOLPS4 formula:
                     ! Nt = A*exp((time-T)*B)

              if(ncall.eq.0) then
                write(*,*)
     &           'compute_feedback_option: using ',
     &           'SOLPS4 exponential decay for sequence ',
     &           b2espcr(fb_species(ifb))
                if(fb_const(ifb).le.0.0_R8)
     &             fb_const(ifb) = fb_target(ifb)
                if(saved_fb_prev(ifb).le.0.0_R8)
     &             saved_fb_prev(ifb) = fb_current(ifb)
              endif

              fb_target(ifb) =
     .           fb_const(ifb) *
     .           exp(fb_beta(ifb)*
     .           (tim-fb_time(ifb)))
              fb_rescale(ifb) =
     .           fb_alpha(ifb) *
     .           ((fb_target(ifb) -
     .           fb_current(ifb))/ddtim +
     .           (saved_fb_prev(ifb) -
     .           fb_current(ifb))/dt_prev)

            case (6) ! rescale slowed by fb_alpha (SOLPS4 style)

              fb_rescale(ifb) =
     &           fb_alpha(ifb) * fb_current(ifb)

            case (7) ! rescale slowed by tanh_log

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_option: using ',
     &           'tanh_log_rescale with plateau for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

! 2**(tanh(log(x)/b)*log(a)/log(2))
              t0 = fb_target(ifb)*(1.0_R8-fb_const(ifb))
              t1 = fb_target(ifb)*(1.0_R8+fb_const(ifb))
              if(fb_current(ifb) < t0) then
                fb_rescale(ifb) = 2.0_R8 **
     &             (tanh(log(t0 /
     &             fb_current(ifb))/
     &             fb_beta(ifb))*
     &             log(fb_alpha(ifb))/log(2.0_R8))
              elseif(fb_current(ifb) > t1) then
                fb_rescale(ifb) = 2.0_R8 **
     &             (tanh(log(t1 /
     &             fb_current(ifb))/
     &             fb_beta(ifb))*
     &             log(fb_alpha(ifb))/log(2.0_R8))
              else
                fb_rescale(ifb) = 1.0_R8
              endif

            case (8) ! rescaling done according to SOLPS4 formula:
                     ! Feedback target changes as:
                     ! S --> S + D, where
                     ! D = A((N_tar - N_cur)/dt + (N_prev - N_cur))/dt_prev
                     ! if actuator is not number of particles it is
                     ! assumed to be linearly proportional to the total number
                     ! of particles in the edge (including neutrals)

              if(ncall.eq.0) then
                write(*,'(5a,i3)')
     &           'compute_feedback_option: using ',
     &           'b2stbc_fb_pre_option: using ',
     &           'SOLPS4 formula and assuming linear relation between ',
     &           'target and total number of particles ',
     &           'in the edge for sequence ',
     &           b2espcr(fb_species(ifb))
                if ((fb_type(ifb).eq.4).or.(fb_type(ifb).eq.5).or.
     &              (fb_type(ifb).eq.14).or.(fb_type(ifb).eq.27)) then
                  write(*,'(2a,i3,a,i3,2a)')
     .             'compute_feedback_option: Warning! ',
     .             'feedback for sequence ',b2espcr(fb_species(ifb)),
     .             ' using fb_type ',fb_type(ifb),
     .             '. No linearisation needed. Falling back to ',
     .             'fb_type 4 with fb_const = 0.0'
                endif
              endif

              if ((fb_type(ifb).ne.4).and.(fb_type(ifb).ne.5).and.
     &            (fb_type(ifb).ne.14).and.(fb_type(ifb).ne.27)) then
                t0 = 0.0_R8
                is_start=-1
                is_end=-1
                do is = 0, ns-1
                  if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                    if(is_start.eq.-1) is_start=is
                    do iCv = 1, mpg%nCi
                      if(mpg%cvOnClosedSurface(iCv)) cycle
                      if(mpg%cvReg(iCv).eq.0) cycle
                      t0 = t0 + pl%na(iCv,is)*geo%cvVol(iCv)
                    enddo
                    is_end=is
                  endif
                enddo
#ifndef SOLPS4_3
                if (switch%use_eirene.ne.0) then
                  do iatm=1,nnatmi
                    if (latmscl(iatm).eq.fb_species(ifb)) then
                      do iCv = 1, mpg%nCi
                        if(mpg%cvOnClosedSurface(iCv)) cycle
                        if(mpg%cvReg(iCv).eq.0) cycle
                        t0 = t0 + dab2(iCv,iatm,1)*geo%cvVol(iCv)
                        do imol = 1, nnmoli
                          if (mlcmp(iatm,imol).gt.0) then
                            t0 = t0 + dmb2(iCv,imol,1)*mlcmp(iatm,imol)*
     &                                geo%cvVol(iCv)
                          endif
                        enddo
                      enddo
                    endif
                  enddo
                endif
#endif
                if (ncall.eq.0.and.saved_fb_prev(ifb).le.0.0_R8)
     &           saved_fb_prev(ifb) = t0
                if (fb_current(ifb).ne.0.0_R8) then
                  fb_rescale(ifb) = fb_alpha(ifb)*
     .             ((fb_target(ifb)/fb_current(ifb) - 1.0_R8)*t0/ddtim +
     .              (saved_fb_prev(ifb) - t0)/dt_prev)
                 else
                   fb_rescale(ifb) = 1.0_R8
                 end if
                 fb_current_prev(ifb) = t0
              else
                if (ncall.eq.0.and.saved_fb_prev(ifb).le.0.0_R8)
     &           saved_fb_prev(ifb) = fb_current(ifb)
                fb_rescale(ifb) = fb_alpha(ifb)*
     .           ((fb_target(ifb) - fb_current(ifb))/ddtim +
     .            (saved_fb_prev(ifb) - fb_current(ifb))/dt_prev)
              endif

            case (9) ! puff depending from pump rescaled by tanh_log with plateau

              if(ncall.eq.0) then
                write(*,'(a,a,a,i3)')
     &           'compute_feedback_option: using ',
     &           'puff depending from pump ',
     &           'tanh_log_rescale with plateau for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

! 2**(tanh(log(x)/b)*log(a)/log(2))
              t0 = fb_target(ifb)*(1.0_R8-fb_const(ifb))
              t1 = fb_target(ifb)*(1.0_R8+fb_const(ifb))
              if(fb_current(ifb) < t0) then
                fb_rescale(ifb) = 2.0_R8 **
     &           (tanh(log(t0 / fb_current(ifb))/ fb_beta(ifb))*
     &            log(fb_alpha(ifb))/log(2.0_R8))
              elseif(fb_current(ifb) > t1) then
                fb_rescale(ifb) = 2.0_R8 **
     &           (tanh(log(t1 / fb_current(ifb))/ fb_beta(ifb))*
     &            log(fb_alpha(ifb))/log(2.0_R8))
              else
                fb_rescale(ifb) = 1.0_R8
              endif

              if (switch%use_eirene.eq.0) then
                call calc_puff_pump_b25(switch,
     &           mpg, nCv, dv, pl, sr,
     &           is_start, is_end, sna_pump, sna_puff, sna_core)
                if (sna_pump < -0.0_R8) then
                  fb_rescale(ifb) = max(0.0_R8,
     &              (abs(sna_pump) - sna_core)*fb_rescale(ifb))
                endif
#ifdef B25_EIRENE
              else
                call calc_puff_pump(switch, ns, mpg, nCv, dv,
     &             b2espcr(fb_species(ifb)),
     &             sna_pump,sna_puff, sna_core)
                  sna_puff = 0.0_R8
                  do istra=1,mpg%nRc
                    if (crcstra(istra).ne.'C') cycle
                    if (istra.eq.fb_istra(ifb)) cycle
                    if (species_start(fb_istra(ifb)).ne.
     &                  species_start(istra)) cycle
                    sna_puff = sna_puff + userfluxparm(istra,1)
                  enddo
                  write(*,'(a,i6,es17.6)') 'Total puff',ifb,sna_puff
                  write(*,'(a,i6,es17.6)') 'Total pump',ifb,sna_pump
                  write(*,'(a,i6,es17.6)') 'Total core',ifb,sna_core
                  fb_rescale(ifb) = max(0.0_R8,
     &              (sna_pump - sna_core)*fb_rescale(ifb) - sna_puff)
#endif
              endif

            case default
              write(*,*) 'fb_rescale_option not coded ',
     &          fb_rescale_option(ifb)

            end select

            write(*,'(a,1x,1p,i3,6g15.6)')
     &        'feedback_rescale ', b2espcr(fb_species(ifb)),
     &        fb_current(ifb),
     &        fb_target(ifb),
     &        fb_rescale(ifb),
     &        saved_fb_prev(ifb),
     &        ddtim, dt_prev

! phase 3:

            select case(fb_actuator(ifb))

            case (0) ! no action

            case (1) ! neutral particle flux => gas puff

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_actuator: using ',
     &           'gas_puff for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              if(ncall.eq.0 .and. saved_fb_actuator(ifb).eq.0.0_R8) then
                if (switch%use_eirene.eq.0) then
                  if (conpar(is_start, fb_ib(ifb), 1).gt.0.0_R8) then
                    saved_fb_actuator(ifb) =
     &               conpar(is_start, fb_ib(ifb), 1)
                    write(message,'(a,i3,a,i2,a)')
     .               'Initial gas puff strength set by CONPAR(',
     .                is_start,',',fb_ib(ifb),',1) !'
                  else
                    saved_fb_actuator(ifb) = 1.0e18_R8
                    write(*,*) 'Warning!  ','The gas puff was ',
     .               'started with an initial rate of 1e18!'
                  endif
                else if (userfluxparm(fb_istra(ifb),1).gt.0.0_R8) then
                  saved_fb_actuator(ifb) =
     &             userfluxparm(fb_istra(ifb),1)
                  write(message,'(a,i3,a,i2,a)')
     .             'Initial gas puff strength set by USERFLUXPARM(',
     .              fb_istra(ifb),',1)  !'
                else
                  saved_fb_actuator(ifb) = 1.0e18_R8
                  write(*,*) 'Warning!  ','The gas puff was ',
     .             'started with an initial rate of 1e18!'
                endif
              endif

              if (fb_rescale_option(ifb).eq.4 .or.
     &            fb_rescale_option(ifb).eq.5 .or.
     &            fb_rescale_option(ifb).eq.8) then
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) + fb_rescale(ifb)
              elseif (fb_rescale_option(ifb).eq.6 .or.
     &                fb_rescale_option(ifb).eq.9) then
                saved_fb_actuator(ifb) = fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) * fb_rescale(ifb)
              end if
              if(fb_puff_min(ifb).ge.0.0_R8)
     >         saved_fb_actuator(ifb) =
     &          max(fb_puff_min(ifb),saved_fb_actuator(ifb))
              if(fb_puff_max(ifb).gt.0.0_R8)
     >         saved_fb_actuator(ifb) =
     &          min(fb_puff_max(ifb),saved_fb_actuator(ifb))

              if(fb_overshoot(ifb).gt.1.0_R8 .and.
     &           fb_current(ifb).gt.fb_target(ifb)*fb_overshoot(ifb))
     >         then
! csc Here the idea is to change conpar if fluid neutrals are used, otherwise userfluxparm
!     This way it is easier to generalize the schemes when hybrid neutrals are used (to be done!)
                if (switch%use_eirene.eq.0) then
                  conpar(is_start,fb_ib(ifb),1) = 0.0_R8
                else
                  userfluxparm(fb_istra(ifb),1) = 0.0_R8
                endif
              else
                if (switch%use_eirene.eq.0) then
                  conpar(is_start,fb_ib(ifb),1) = saved_fb_actuator(ifb)
                else
                  userfluxparm(fb_istra(ifb),1) = saved_fb_actuator(ifb)
                endif
              endif

            case (2) ! rescale na

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_actuator: using ',
     &           'density rescale for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              do is = is_start, is_end
                pl%na(:,is) = pl%na(:,is) * fb_rescale(ifb)
                dv%kinrgy(:,is) =
     &           dv%kinrgy(:,is) * fb_rescale(ifb)
                dv%ne=dv%ne+(fb_rescale(ifb)-1.0_R8)*
     &           pl%na(:,is)*rt%rza(:,is)
                dv%fna(:,0:1,is) =
     &           dv%fna(:,0:1,is)*fb_rescale(ifb)
                dv%fna_mdf(:,0:1,is) =
     &           dv%fna_mdf(:,0:1,is)*fb_rescale(ifb)
                dv%fna_fcor(:,0:1,is) =
     &           dv%fna_fcor(:,0:1,is)*fb_rescale(ifb)
                dv%fna_nodrift(:,0:1,is) =
     &           dv%fna_nodrift(:,0:1,is)*fb_rescale(ifb)
                dv%fna_he(:,0:1,is) =
     &           dv%fna_he(:,0:1,is)*fb_rescale(ifb)
                dv%fnaPSch(:,0:1,is) =
     &           dv%fnaPSch(:,0:1,is)*fb_rescale(ifb)
                dv%fna_eir(:,0:1,is) =
     &           dv%fna_eir(:,0:1,is)*fb_rescale(ifb)
                dv%fna_32(:,0:1,is) =
     &           dv%fna_32(:,0:1,is)*fb_rescale(ifb)
                dv%fna_52(:,0:1,is) =
     &           dv%fna_52(:,0:1,is)*fb_rescale(ifb)
              end do
*   ..initialise previous state and fluxes
              psnl%na = pl%na
              psnc%na = pl%na
              psnl%ne = dv%ne
              psnc%ne = dv%ne
              psnl%fna = dv%fna
              psnc%fna = dv%fna
              psnl%kinrgy = dv%kinrgy
              psnc%kinrgy = dv%kinrgy
*   ..compute initial ni
              call b2xpni (nCv, ns, pl%na, dv%ni)
              call b2xpnn (nCv, ns, pl%na, dv%nn)
*   ..compute ni0
              psnl%ni = dv%ni
              psnc%ni = dv%ni

            case (3) ! charged species particle flux

              if(ncall.eq.0) then
               write(*,'(a,a,i3)')
     &          'compute_feedback_actuator: using ',
     &          'boundary flux condition for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              is_start = eb2spcr(b2espcr(fb_species(ifb)))+1
              is_end = eb2spcr(b2espcr(fb_species(ifb))) +
     &          nfluids(b2espcr(fb_species(ifb)))

              if (feedback_namelist_used .and.
     &            ncall .eq. 0 .and.
     &            saved_fb_actuator(ifb) .eq. 0.0_R8)
     &         saved_fb_actuator(ifb) =
     &          conpar(is_end, fb_ib(ifb), 1)

              if (fb_type(ifb).eq.6) then !csc, this is for backward compatibility, can probably be removed
                saved_fb_actuator(ifb) =
     &           (1.0_R8 - fb_alpha(ifb)) *
     &           saved_fb_actuator(ifb) +
     &           fb_rescale(ifb)
              else
                if (fb_rescale_option(ifb).lt.4) then
                  saved_fb_actuator(ifb) =
     &             saved_fb_actuator(ifb) * fb_rescale(ifb)
                else
                  saved_fb_actuator(ifb) =
     &             fb_rescale(ifb)
                endif
              end if

              if (fb_puff_min(ifb) .ge. 0.0_R8)
     &         saved_fb_actuator(ifb) =
     &          max(fb_puff_min(ifb),saved_fb_actuator(ifb))
              if (fb_puff_max(ifb) .gt. 0.0_R8)
     &         saved_fb_actuator(ifb) =
     &          min(fb_puff_max(ifb),saved_fb_actuator(ifb))

              if (fb_overshoot(ifb) .gt. 1.0_R8 .and.
     &            fb_current(ifb) .gt.
     &            fb_target(ifb) * fb_overshoot(ifb)) then
                total_flux = 0.0_R8
              else
                total_flux = saved_fb_actuator(ifb)
              endif

              if (fb_type(ifb).eq.6) then
                do is = is_start, is_end
                  iss = is-is_start+1
                  conpar(is, fb_ib(ifb), 1) =
     &             total_flux * charge_frac(iss)
                end do
              else
                conpar(fb_species(ifb), fb_ib(ifb), 1) = total_flux
              endif

            case (4) ! electron heat flux

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &           'compute_feedback_actuator: using ',
     &           'electron heat flux boundary condition'
              endif

              if (feedback_namelist_used .and.
     &            ncall .eq. 0 .and.
     &            saved_fb_actuator(ifb) .eq. 0.0_R8)
     &         saved_fb_actuator(ifb) = enepar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &           fb_rescale(ifb)
              end if

              enepar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (5) ! ion heat flux

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &           'compute_feedback_actuator: using ',
     &           'ion heat flux boundary condition'
              endif

              if (feedback_namelist_used .and.
     &            ncall .eq. 0 .and.
     &            saved_fb_actuator(ifb) .eq. 0.0_R8)
     &         saved_fb_actuator(ifb) = enipar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &           fb_rescale(ifb)
              end if

              enipar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (6) ! current

              if(ncall.eq.0) then
                write(*,'(a,a)')
     &           'compute_feedback_actuator: using ',
     &           'electric current boundary condition'
              endif

              if (feedback_namelist_used .and.
     &            ncall .eq. 0 .and.
     &            saved_fb_actuator(ifb) .eq. 0.0_R8)
     &         saved_fb_actuator(ifb) = potpar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &           fb_rescale(ifb)
              end if

              potpar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (7) ! density

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_actuator: using ',
     &           'density condition for sequence ',
     &            b2espcr(fb_species(ifb))
              endif

              if (ncall .eq. 0 .and.
     &            saved_fb_actuator(ifb) .eq. 0.0_R8)
     &         saved_fb_actuator(ifb) =
     &          conpar(fb_species(ifb), fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &           saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &           fb_rescale(ifb)
              end if

              conpar(fb_species(ifb),fb_ib(ifb), 1) =
     &         saved_fb_actuator(ifb)

            case (8) ! pump

              if(ncall.eq.0) then
                if (fb_ib(ifb).eq.0) then
                  write(*,'(a,a)')
     &             'compute_feedback_actuator: using ',
     &             'gas pump for all isonuclear sequences'
                elseif (fb_ib(ifb).lt.0) then
                  write(*,'(a,a,i3,a,i3)')
     &             'compute_feedback_actuator: using ',
     &             'gas pump for isonuclear sequence ',
     &              b2espcr(fb_species(ifb)),
     &             ' at non-default standard surface number',
     &              abs(fb_ib(ifb))
                else
                  write(*,'(a,a,a,i3)')
     &             'compute_feedback_actuator: using ',
     &             'gas pump for all isonuclear sequences ',
     &             'at non-default standard surface number',
     &              fb_ib(ifb)
                endif
              endif

              if(ncall.eq.0 .and. saved_fb_actuator(ifb).eq.0.0_R8) then
                saved_fb_actuator(ifb) = fb_puff_max(ifb)
                write(*,*) 'Warning! The pump was ',
     &           'started with an initial rate of ',
     &            saved_fb_actuator(ifb)
              endif

              if (fb_rescale_option(ifb).lt.4 .or.
     &            fb_rescale_option(ifb).eq.7) then
                if (switch%use_eirene.eq.0) then
                  call calc_puff_pump_b25(switch,
     &             mpg, nCv, dv, pl, sr,
     &             is_start, is_end, sna_pump, sna_puff, sna_core)
                  if (sna_pump < 0.0_R8) then
                    saved_fb_actuator(ifb) =
     &               abs((sna_puff+sna_core)/sna_pump)
     &               *saved_fb_actuator(ifb)/fb_rescale(ifb)
                  endif
#ifdef B25_EIRENE
                else
                  call calc_puff_pump(switch, ns, mpg, nCv, dv,
     &             fb_species(ifb), sna_pump, sna_puff, sna_core)
                  if (mod(ncall,eirene_mod).eq.0 .and.
     &                ncall.gt.eirene_mod) then
                    t0 = 1.0_R8/fb_rescale(ifb) * (sna_puff+sna_core)
     &               / sna_pump
                    saved_fb_actuator(ifb) = saved_fb_actuator(ifb)*t0
                  endif
#endif
                endif
              else
                saved_fb_actuator(ifb) =
     &           fb_rescale(ifb)
              end if

    !           if(fb_puff_min(ifb).le.0.0_R8)
    !  &
              saved_fb_actuator(ifb) =
     &         max(fb_puff_min(ifb),saved_fb_actuator(ifb))
    !           if(fb_puff_max(ifb).le.0.0_R8)
    !  &
              saved_fb_actuator(ifb) =
     &         min(fb_puff_max(ifb),saved_fb_actuator(ifb))

              if(fb_overshoot(ifb).gt.1.0_R8 .and.
     &           fb_current(ifb).gt.fb_target(ifb)*fb_overshoot(ifb))
     >         then
! csc Here the idea is to change conpar if fluid neutrals are used, otherwise userfluxparm
!     This way it is easier to generalize the schemes when hybrid neutrals are used (to be done!)
              else
                if (switch%use_eirene.eq.0) then
                  if (switch%recycle_afn.eq.0 .or.
     &                nint(zn(fb_species(ifb))).ne.1) then
                    do istra = 1, nbc
                      if (bccon(is_start,istra).eq.10 .and.
     &                   (istra.eq.fb_ib(ifb).or.fb_ib(ifb).eq.0)) then
                        do ibc = 1, mpg%bcCvP(istra,2)
                          conpar(is_start,istra,1) =
     &                     -saved_fb_actuator(ifb)
                        enddo
                      endif
                    enddo
                  else
                    do istra = 1, nstrai
                      if (((1.0_R8-recyc(is_start,istra)).ge.1.0e-7_R8)
     &                  .and. (species_start(istra) .eq. is_start).and.
     &                   (istra.eq.fb_ib(ifb).or.fb_ib(ifb).eq.0)) then
                        b2recyc(is_start,istra) =
     &                   (1.0_R8-saved_fb_actuator(ifb))
                        recyc(is_start,istra) =
     &                   b2recyc(is_start,istra)
                      endif
                    enddo
                  endif
#ifdef B25_EIRENE
                elseif (mod(ncall,eirene_mod).eq.0 .and.
     &                  ncall.gt.0) then
                  call eirene_recyct_update(b2espcr(fb_species(ifb)),
     &             saved_fb_actuator(ifb),fb_ib(ifb))
#endif
                endif
              endif

            case default
              write(*,*) 'fb_actuator not coded ',
     &          fb_type(ifb)

            end select

            if (fb_actuator(ifb).eq.1 .or.
     &          fb_actuator(ifb).eq.3 .or.
     &          fb_actuator(ifb).eq.8)
     >        write(*,'(a,1p,i3,g15.6)') 'feedback_actuator ',
     &         b2espcr(fb_species(ifb)), saved_fb_actuator(ifb)

          endif

          feedback_namelist_used = feedback_namelist_used .or.
     &     fb_actuator(ifb).eq.1 .or. fb_actuator(ifb).eq.3 .or.
     &     fb_actuator(ifb).eq.8

          if (fb_current_prev(ifb).gt.0.0_R8) then
            saved_fb_prev(ifb) = fb_current_prev(ifb)
          else
            saved_fb_prev(ifb) = fb_current(ifb)
          endif

        enddo !ifb

        dt_prev = ddtim
      endif


      ncall = ncall + 1
      call subend ()
      return
      end subroutine compute_feedback

      subroutine dealloc_feedback
      implicit none
      if (feedback_allocated) then
        deallocate(fbRegP,fbReg)
        feedback_allocated = .false.
      endif
      end subroutine dealloc_feedback

      subroutine do_faces(nFc,mpg,ifb,ifbreg)
      implicit none
      type (mapping), intent (in) :: mpg
      integer, intent (in) :: ifb, nFc
      integer, intent (inout) :: ifbreg
      integer iFc

      do iFc = 1, nFc
        if ((mpg%fcReg(iFc).eq.fb_reg_par(ifb,1)) .or.
     &      (fb_reg_par(ifb,2).gt.0.and.
     &      (mpg%fcReg(iFc).eq.fb_reg_par(ifb,2)))) then
          fbReg(ifbreg) = iFc
          ifbreg = ifbreg + 1
         endif
      end do

      return
      end subroutine do_faces

      subroutine do_volumes(nCv,mpg,ifb,ifbreg)
      implicit none
      type (mapping), intent (in) :: mpg
      integer, intent (in) :: ifb, nCv
      integer, intent (inout) :: ifbreg
      integer iCv

      do iCv = 1, nCv
        if ((mpg%cvReg(iCv).eq.fb_reg_par(ifb,1)) .or.
     &      (fb_reg_par(ifb,2).gt.0.and.
     &      (mpg%cvReg(iCv).eq.fb_reg_par(ifb,2)))) then
          fbReg(ifbreg) = iCv
          ifbreg = ifbreg + 1
         endif
      end do

      return
      end subroutine do_volumes

***********************************************************************

      subroutine conv_fbparms (m)
      use b2mod_feedback
      implicit none
      type (mapping), intent (inout) :: m
      integer :: is, ix, iy, iCv, ifb, iss, ifbreg
      logical :: done

      nFb = nspecies !number of feedback imposed
      call xertst(nFb.le.DEF_NATM,
     & 'b2us_feedback: increase size of DEF_NATM')
      fb_type(1:nfb) = na_feedback_choice(0:nspecies-1)
      fb_rescale_option(1:nfb) = na_feedback_option(0:nspecies-1)
      fb_actuator(1:nfb) = na_feedback_actuator(0:nspecies-1)
      fb_target(1:nfb) = na_feedback_target(0:nspecies-1)
      fb_alpha(1:nfb) = na_feedback_alpha(0:nspecies-1)
      fb_beta(1:nfb) = na_feedback_beta(0:nspecies-1)
      fb_const(1:nfb) = na_feedback_const(0:nspecies-1)
      fb_time(1:nfb) = na_feedback_time(0:nspecies-1)
      fb_ib(1:nfb) = na_feedback_ib(0:nspecies-1)
      do ifb = 1, nfb
        if (fb_ib(ifb).eq.-1.and.fb_actuator(ifb).ne.8)
     >   fb_ib(ifb) = -1-DEF_NSTS
      end do
      fb_puff_min(1:nfb) = na_feedback_puff_min(0:nspecies-1)
      fb_puff_max(1:nfb) = na_feedback_puff_max(0:nspecies-1)
      fb_overshoot(1:nfb) = na_feedback_overshoot(0:nspecies-1)
      ifbreg = 1
      do is = 0, nspecies-1
        fb_species(is+1) = eb2spcr(is+1) !new format specifies feedbacks on fluid species
        if (na_feedback_ix1(is).gt.-2.or.na_feedback_ix2(is).gt.-2.or.
     &      na_feedback_iy1(is).gt.-2.or.na_feedback_iy2(is).gt.-2) then
          done = .false.
          ! check if same region has already been done
          do iss = 0, is-1
            if (na_feedback_ix1(is).eq.na_feedback_ix1(iss) .and.
     &          na_feedback_ix2(is).eq.na_feedback_ix2(iss) .and.
     &          na_feedback_iy1(is).eq.na_feedback_iy1(iss) .and.
     &          na_feedback_iy2(is).eq.na_feedback_iy2(iss)) then
              fb_RegP(is+1,:) = fb_RegP(iss+1,:)
              done = .true.
            endif
          end do
          if (done) cycle
          fb_RegP(is+1,1) = ifbreg
          do ix = na_feedback_ix1(is),na_feedback_ix2(is)
            do iy = na_feedback_iy1(is),na_feedback_iy2(is)
              call xertst(ifbreg.le.nMxFbReg,
     &          'b2us_feedback: increase size of nMxFbReg')
              iCv = m%imapCv(ix,iy)
              if (iCv.eq.0) cycle
              fb_Reg(ifbreg) = iCv
              ifbreg = ifbreg + 1
            enddo
          enddo
          fb_RegP(is+1,2) = ifbreg - fb_RegP(is+1,1)
        endif
      end do
      nfbreg = ifbreg - 1

      return
      end subroutine conv_fbparms

***********************************************************************

      subroutine write_fbparms
      use b2mod_version
      implicit none
      integer :: ifb, iout
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

#ifndef LEGACYCOMP
      open(newunit=iout,file='b2us.feedback_control.parameters')
#else
      iout = newunit()
      open(unit=iout,file='b2us.feedback_control.parameters')
#endif
      call cfverw(iout,newversion)

      write (iout,'(1x,a)') '&feedback_control'
      write (iout,'(1x,a,i0)') 'NFB     = ',nFb

      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'FB_TYPE =  ',(fb_type(ifb),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T24,10(i5,a1)))') 'FB_RESCALE_OPTION =  ',
     .       (fb_rescale_option(ifb),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T17,10(i5,a1)))')
     .       'FB_ACTUATOR =  ',(fb_actuator(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_TARGET=',(fb_target(ifb),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'FB_IB =  ',(fb_ib(ifb),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'FB_SPECIES =  ',(fb_species(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_ALPHA=',(fb_alpha(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_BETA=',(fb_beta(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_PUFF_MIN=',(fb_puff_min(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_PUFF_MAX=',(fb_puff_max(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T16,10(es10.2,a1)))')
     .       'FB_OVERSHOOT=',(fb_overshoot(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_TIME=',(fb_time(ifb),',',ifb=1,nFb)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'FB_CONST=',(fb_const(ifb),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T24,10(i5,a1)))')
     .       'FB_REG_PAR(1,1) =  ',(fb_reg_par(ifb,1),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T24,10(i5,a1)))')
     .       'FB_REG_PAR(1,2) =  ',(fb_reg_par(ifb,2),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T20,10(i5,a1)))')
     .       'FB_REGP(1,1) =  ',(fb_regp(ifb,1),',',ifb=1,nFb)
      write(iout,'(1x,a,(1x,T20,10(i5,a1)))')
     .       'FB_REGP(1,2) =  ',(fb_regp(ifb,2),',',ifb=1,nFb)
      if (nFbReg.gt.0)
     & write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     &       'FB_REG =  ',(fb_reg(ifb),',',ifb=1,nFbReg)

      write (iout,'(1x,a)') '/'
      close (iout)

      return
      end subroutine write_fbparms

***********************************************************************

      subroutine write_fbsave
      use b2mod_version
      implicit none
      integer :: is, iout
      character*26 :: hlp_frm
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

#ifndef LEGACYCOMP
      open(newunit=iout,file='b2us.feedback_save.parameters')
#else
      iout = newunit()
      open(unit=iout,file='b2us.feedback_save.parameters')
#endif
      call cfverw(iout,newversion)

      write (iout,'(1x,a)') '&feedback'
      write(hlp_frm,'(a,i2,a)') '(1x,a,(1x,',nFb,'(es22.15,a1)))'
      if (maxval(saved_fb_actuator).gt.0.0_R8)
     > write(iout,hlp_frm) 'SAVED_FB_ACTUATOR = ',
     .  (saved_fb_actuator(is),',',is=1,nFb)
      if (maxval(saved_fb_prev).gt.0.0_R8)
     > write(iout,hlp_frm) 'SAVED_FB_PREV = ',
     .  (saved_fb_prev(is),',',is=1,nFb)

      write (iout,'(1x,a)') '/'
      close (iout)

      return
      end subroutine write_fbsave

      subroutine calc_puff_pump_b25(switch, mpg, nCv, dv, pl, sr,
     & is_start, is_end, sna_pump, sna_puff, sna_core)
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      implicit none
      real (kind=R8) :: sna_pump, sna_puff, sna_core
      integer, intent(in) :: nCv, is_start, is_end
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (B2Plasma), intent (inout) :: pl
      type (B2Source), intent (in) :: sr
      type (switches), intent (in) :: switch
      integer :: i, is, iFc, iCv, iCv1, ib, istra, find_loc

      sna_core = 0.0_R8
      do is = is_start, is_end
        do iCv1 = mpg%nCi+1, nCv
          iFc = mpg%cvFc(mpg%cvFcP(iCv1,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if(mpg%cvOnClosedSurface(iCv)) then
            if (switch%mdf_fnb.ne.0) then
              sna_core = sna_core + dv%fna_mdf(iFc,1,is)
            else
              sna_core = sna_core + dv%fna(iFc,1,is)
            end if
          endif
        enddo
      enddo

      sna_puff = 0.0_R8
      is = is_start
      do ib = 1, nbc
        if (BCCON(IS,IB).eq.8) then
          sna_puff = sna_puff+CONPAR(is,ib,1)
        endif
      enddo

      sna_pump = 0.0_R8

      if (switch%recycle_afn.eq.0.or.nint(zn(is_start)).ne.1) then
        do istra = 1, nbc
          if (bccon(is_start,istra).eq.10) then
            do ib = 1, mpg%bcCvP(istra,2)
              do is = is_start, is_end
                iCv1 = mpg%bcCv(mpg%bcCvP(istra,1)+ib-1,1)
                sna_pump = sna_pump + sr%sna(iCv1,0,is)
     &           +pl%na(iCv1,is)*sr%sna(iCv1,1,is)
              enddo
            enddo
          endif
        enddo

      else
        do istra = 1, nstrai
          if (((1.0_R8-recyc(is_start,istra)).ge.1.0e-7_R8)
     &     .and. (species_start(istra) .eq. is_start)) then
            do is=is_start,is_end
              ib = find_loc(BCSTART, NBC, RCSTART(istra))
              write(*,*) "istra, ib index",RCSTART(istra), ib
              do i = 1, mpg%bcCvP(ib,2)
                iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+i-1,1)
                sna_pump = sna_pump + sr%sna(iCv1,0,is)
     &           +pl%na(iCv1,is)*sr%sna(iCv1,1,is)
              enddo
            enddo
          endif
        enddo
      endif

      write(*,*) 'pump ', sna_pump,
     & ' puff+core ', sna_core+sna_puff

      return
      end subroutine calc_puff_pump_b25

#ifdef B25_EIRENE
      subroutine calc_puff_pump(switch, ns, mpg, nCv, dv,
     & ispecies, sna_pump, sna_puff, sna_core)

      use b2mod_neutrals_namelist
      use eirmod_wneutrals
     , , only : eirpump
      use b2mod_switches
      implicit none

      integer, intent(in) :: nCv, ns, ispecies
      real (kind=R8) :: sna_pump, sna_puff, sna_core
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer :: is, iFc, iCv, iCv1, k, im

      sna_pump = 0.0_R8
      sna_pump = sna_pump + eirpump(ispecies)
      do im = 1, nnmoli
        if (mlcmp(ispecies,im).gt.0) then
          sna_pump = sna_pump +
     &     mlcmp(ispecies,im) * eirpump(nnatmi+im)
        endif
      end do

      sna_puff = 0.0_R8
      do is = 1, nstrat
        if(latmscl(ispecies).eq.0) cycle
        if(crcstra(is).eq.'C') then
          if(gpfc(ispecies,is).ne.0.0_R8) sna_puff =
     &     sna_puff + userfluxparm(is,1)*gpfc(ispecies,is) !*recyctmp
        end if
      end do

      sna_core = 0.0_R8
      do is = 0, ns
        k = b2espcr(is)
        if (is_neutral(is) .or. k.ne.ispecies) cycle
        do iCv1 = mpg%nCi+1, nCv
          iFc = mpg%cvFc(mpg%cvFcP(iCv1,1))
          iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
          if(mpg%cvOnClosedSurface(iCv)) then
            if (switch%mdf_fnb.ne.0) then
              sna_core = sna_core + dv%fna_mdf(iFc,1,is)
            else
              sna_core = sna_core + dv%fna(iFc,1,is)
            end if
          endif
        enddo
      enddo

      write(*,*) 'pump ', sna_pump,
     & ' puff+core ', sna_core+sna_puff

      return
      end subroutine calc_puff_pump

      subroutine eirene_recyct_update(ispecies, factor, ists)
#ifndef SOLPS4_3
      use eirmod_clgin             ! IGNORE
#else
      use clgin                    ! IGNORE
#endif
      use b2mod_neutrals_namelist
      implicit none
      integer, intent(in) :: ispecies         ! Index of B2.5 species (or isonuclear sequence)
      integer, intent(in) :: ists             ! number of non-default standard surface
      real (kind=R8), intent(in) :: factor    ! factor to multiply the absorption coefficient
                                              ! note that recyct is (1.0 - absorption coefficient)
      logical :: set_for_all_sorts, is_albedo_changed

      real (kind=R8) :: t0
      integer ja, l, jm, k, ilim

      is_albedo_changed = .false.
      if (ists >= 0) set_for_all_sorts = .true.
      if (ists < 0) set_for_all_sorts = .false.
      ilim = nnlimi + abs(ists)

      write (*,*) "Change absorption for all", set_for_all_sorts

      do l = 1, nnlimi + nsts  ! loop over Eirene surfaces
        if (ists.ne.0 .and. l.ne.ilim) cycle
        write (*,*) "Scan ilim", ilim, ists
        do ja = 1,nnatmi
          write (*,*) "ja", ja,latmscl(ja),ispecies
          if (latmscl(ja).eq.ispecies .or. set_for_all_sorts) then       ! Eirene atom corresponds to B2.5 isonuclear sequence
            if (recyct(ja,l) > 0.0_R8 .and. recyct(ja,l) < 1.0_R8) then  ! this condition is required to not disturb those surfaces where albedo is exactly 1.0
              t0 = 1.0_R8 - factor
              if (t0 > 0.0_R8 .and. t0 < 1.0_R8) recyct(ja,l) = t0       ! again, not to make albedo =1.0 for those surfaces where it was set less than 1.0
    !           if (t0 > 0.0_R8 .and. t0 < 1.0_R8)
    !  &         write(*,*) "surf:", l, ja, recyct(ja,l)
              is_albedo_changed = .true.
            endif
            do jm = 1,nnmoli                                             ! check whether this atom constitutes a part of the molecule and update the albedo for such molecule in the same manner
              if (mlcmp(ja,jm) .gt. 0) then
                if (recyct(jm+nnatmi,l) > 0.0_R8 .and.
     &              recyct(jm+nnatmi,l) < 1.0_R8) then
                  t0 = 1.0_R8 - factor
                  if (t0 > 0.0_R8 .and. t0 < 1.0_R8)
     >             recyct(jm+nnatmi,l) = t0
    !               if (t0 > 0.0_R8 .and. t0 < 1.0_R8)
    !  &             write(*,*) "surf:",  l, jm, recyct(jm+nnatmi,l)
                  is_albedo_changed = .true.
                endif
              endif
            enddo
          endif
        enddo
      enddo

      call xertst(is_albedo_changed,
     &  "feedback: albedo has not been changed anywhere!"//
     &   " Check the FB_IB!")
      return
      end subroutine eirene_recyct_update
#endif

      end module b2us_feedback
