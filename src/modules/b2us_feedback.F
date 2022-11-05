      module b2us_feedback
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
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
     & , only : nomp, omp, icsepomp, ifssep, nimp, imp, icsepimp
#ifdef B25_EIRENE
      use eirmod_extrab25
     & , only : wldna, wldra, wldpa, wldnm, wldrm, wldpm
#endif
      implicit none
      integer, parameter :: nvac=4, nvacreg=4, nMxFbReg = 1000
#     include "DIMENSIONS.F"

*   ..variables
      real (kind=R8), save ::
     * cum_volrec, dt_prev, ddtim
      logical, save :: feedback_namelist_used
      logical, save :: feedback_read, control_read
      logical, save :: feedback_allocated=.false.
      real (kind=R8) :: charge_frac(DEF_NSD)

      integer, save :: icped = 0, b2sral_style = 2, ncallfeedback = 0

      integer, allocatable, save :: fbRegP(:,:) ! (nFb,2) pointing for each feedback to the first index in the fbReg list, and its number of CVs or faces
      integer, allocatable, save :: fbReg(:)  !         listing for each feedback the cells or faces
      integer, save :: nfbreg = 0, nfb = 0
*   ..namelists
! csc new namelist merging all feedback schemes
      real (kind=R8), save ::
     & fb_target(DEF_NATM) = 0.0_R8, fb_time(DEF_NATM) = 0.0_R8,
     & saved_fb_actuator(DEF_NATM) = 0.0_R8,
     & fb_alpha(DEF_NATM) = 0.001_R8, fb_beta(DEF_NATM) = 1.0_R8,
     & fb_const(DEF_NATM) = 0.0_R8, fb_puff_min(DEF_NATM) = 0.0_R8,
     & fb_puff_max(DEF_NATM) = 0.0_R8, fb_overshoot(DEF_NATM) = 0.0_R8,
     & fb_current(DEF_NATM) = 0.0_R8, fb_prev(DEF_NATM) = 0.0_R8,
     & fb_rescale(DEF_NATM) = 1.0_R8

      integer, save :: fb_reg_par(DEF_NATM,2) = 0, !contains face label to define list of faces for feedback
     & fb_species(DEF_NATM) = 0, fb_ib(DEF_NATM) = -1,
     & fb_type(DEF_NATM) = 0, fb_rescale_option(DEF_NATM) = 0,
     & fb_actuator(DEF_NATM) = 0, fb_istra(DEF_NATM) = 0,
     & fb_reg(nMxFbReg) = 0, fb_regp(DEF_NATM,2) = 0 !these two arrays have same meaning as fbreg and fbregp above but are needed for conversion and namelist

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

      namelist /feedback/ saved_fb_actuator
! new namelist for specifying more complicated feedback options
      namelist /feedback_control/
! csc new namelist merging all feedback schemes
     & nfb,
     & fb_target, fb_type, fb_rescale_option,
     & fb_actuator, fb_alpha, fb_beta,
     & fb_const, fb_time, 
     & fb_ib, fb_species,
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
      found = maxval(saved_fb_actuator).gt. 0.0_R8
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
c     .            saved_cbsna_sol, maxval(saved_na_feedback_actuator))
c     .            .gt.0.0_R8

      if (edition.eq.'    ') then
        filename = 'b2.feedback_save.parameters'
        feedback_read = .false.
      else
        filename = 'b2.feedback_save.parameters.'//edition
      endif
cxpb We do the write explicitly to avoid array overflows when the problem
cxpb size is changed in DIMENSIONS.F
      if(feedback_namelist_used) then
        open(nwrite,file=trim(filename))
        write(nwrite,'(a9)') '&FEEDBACK'
        write(hlp_frm,'(a,i2,a)')
     >   '(a,1p,',nFb,'(1x,1e15.8,'',''))'
        write(nwrite,hlp_frm) 'SAVED_FB_ACTUATOR=',
     >   (saved_fb_actuator(is),is=1,nFb)
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

      subroutine init_feedback(nCv, nFc, nVx, ns, natmi, switch, geo,
     & mpg)
      implicit none
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      integer :: nCv, nFc, nVx, ns, natmi
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
#ifdef WG_RM
        call xertst (iyped.gt.-2 .and. iyped.le.mpg%ny 
     &   .and. iyped.le.nomp, 'faulty input iyped')
#endif
        call xertst (mpg%cvOnClosedSurface(omp(icped)),
     .   'location of top of pedestal outside core region!')
      endif
      call ipgeti ('b2stbc_ncallfeedback', ncallfeedback)
      call ipgetr ('b2stbc_ndes', ndes)
      call xertst (0.0_R8.le.ndes, 'faulty internal parameter ndes')
#ifdef WG_RM
      call ipgetr ('b2stbc_cbsnafac', cbsnafac)
      call xertst (0.0_R8.lt.cbsnafac,
     1  'faulty internal parameter cbsnafac')
#endif
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
        if (ndes.gt.0.0_R8.or.nesepm_pfr.gt.0.0_R8.or.
     &   private_flux_puff.gt.0.0_R8.or.nesepm_sol.gt.0.0_R8.or.
     *   volrec_sol.gt.0.0_R8.or.ndes_sol.gt.0.0_R8.or.
     &   nepedm_sol.gt.0.0_R8) then
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
        write(*,*) 'WARNING! Feedback b2stbc_fnaycore not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.10) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_fheycore not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.11) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_fchycore not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.12) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_fchycore not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.13) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_ndes not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.14) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_nesepm not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_nesepm_sol not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_ndes_sol not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.14) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_nesepm_pfr not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.3) ifbb=ifb !to override what is present..
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
          if (fb_type(ifb).eq.13) ifbb=ifb !to override what is present..
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
          if (fb_type(ifb).eq.15) ifbb=ifb !to override what is present..
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
        write(*,*) 'WARNING! Feedback b2stbc_nepedm_sol not fully'//
     &  'backward compatible in WG'
        ifbb=0
        do ifb=1,nFb
          if (fb_type(ifb).eq.16) ifbb=ifb !to override what is present..
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
     & 'b2mod_feedback: increase size of DEF_NATM')

      do ifb = 1, nFb
        if (fb_actuator(ifb).ne.1) cycle
        if(fb_istra(ifb).gt.0) then
          call xertst (fb_istra(ifb).le.nstrat,
     1      'faulty internal parameter fb_istra')
          call xertst (crcstra(fb_istra(ifb)).eq.'C',
     1      'Stratum fb_istra is not declared as a puff stratum!')
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

      consistent = .true.
      do ifb = 1, nFb
        if (fb_ib(ifb).le.0 .and. .not.fb_actuator(ifb).eq.2)
     &    call xerrab('Init_feedback: fb_ib<=0 and actuator=/2')
        call xertst(0.le.fb_species(ifb).and.fb_species(ifb).le.ns-1,
     &    'Init_feedback: fb_species>ns!')
        iatm = b2espcr(fb_species(ifb))
        is0 = eb2spcr(iatm)
        do ifbb = ifb+1, nFb
          is00 = eb2spcr(b2espcr(fb_species(ifbb)))
          call xertst(.not.(is0.eq.is00),'init_feedback: cannot '//
     &     'define multiple feedbacks on same atom sequence')
          if ( (fb_type(ifb).eq.fb_type(ifbb)).and. .not.
     &         ((fb_type(ifb).eq.1).or.(fb_type(ifb).eq.4).or.
     &          (fb_type(ifb).eq.5).or.(fb_type(ifb).eq.7).or.
     &          (fb_type(ifb).eq.8)) ) call xerrab('init_feedback:'//
     &      ' set of feedbacks cannot be applied multiple times')
        enddo
        is_start = is0+1
        is_end = is0+nfluids(iatm)
        consistent = consistent .and.
     &   ((fb_target(ifb) .eq. 0.0_R8 .or.
     &     fb_type(ifb) .eq. 0 .or.
     &     fb_rescale_option(ifb) .eq. 0 .or.
     &     fb_actuator(ifb) .eq. 0) .or.
! csc For fluid neutrals (or fluid neutral strata, TBD) change directly the BC
!     For kinetic neutrals (or kinetic strata, TBD), change directly userfluxparm, no issue with fluid BC 
!     Can be arbitrarily expanded to allow applying feedback to any already existent BC
     &    (fb_actuator(ifb) .eq. 1 .and. !neutral flux (or gas puff)
     &     ( bccon(is0,fb_ib(ifb)) .eq. 8 .or.
     &      switch%use_eirene.ne.0) ) .or.
     &    fb_actuator(ifb) .eq. 2 .or.
     &    (fb_actuator(ifb) .eq. 3 .and. !charged particle flux ONLY for fb_type 6, to be improved?
     &     maxval(bccon(is_start:is_end,
     .                  fb_ib(ifb))).eq.8 .and.
     &     minval(bccon(is_start:is_end,
     .                  fb_ib(ifb))).eq.8 .and.
     &     fb_type(ifb).eq.6) .or.
     &    (fb_actuator(ifb) .eq. 3 .and. !charged particle flux, more general
     &     (bccon(fb_species(ifb),fb_ib(ifb)) .eq. 8 .or. 
     &      bccon(fb_species(ifb),fb_ib(ifb)) .eq. 13 ) .and.
     &     .not.fb_type(ifb).eq.6) .or.
     &    (fb_actuator(ifb) .eq. 4 .and. !electron heat flux
     &     bcene(fb_ib(ifb)) .eq. 8) .or.
     &    (fb_actuator(ifb) .eq. 5 .and. !ion heat flux
     &     bceni(fb_ib(ifb)) .eq. 8) .or.
     &    (fb_actuator(ifb) .eq. 6 .and. !current
     &     bcpot(fb_ib(ifb)) .eq. 8) .or.
     &    (fb_actuator(ifb) .eq. 7 .and. !density
     &     (bccon(fb_species(ifb),fb_ib(ifb)) .eq. 1 .or.
     &      bccon(fb_species(ifb),fb_ib(ifb)) .eq. 28)))
      end do
      consistent = consistent .or. lbndusr .or. b2sral_style.eq.1 .or.
     & switch%b2stbc_boundary_namelist.lt.1
      do ifb = 1, nFb
        if (fb_ib(ifb).ge.0) cycle
        consistent = consistent .and.
     &    .not.(fb_actuator(ifb).eq.1 .or.
     &          fb_actuator(ifb).eq.3) .and.
     &    .not.(fb_rescale_option(ifb).eq.4)
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
            ifbreg = ifbreg +   fb_RegP(ifb,2)
            done = .true.
          endif
        ! now look for feedback types requiring faces
          if ((fb_type(ifb).eq.6 .or. fb_type(ifb).eq.10 .or.  
     &        fb_type(ifb).eq.11 .or. fb_type(ifb).eq.12 .or.
     &        fb_type(ifb).eq.13)
     &        .and. .not.done) then
            call xertst(fb_reg_par(ifb,1).gt.0,
     &        'init_feedback: At lest one segment must be specified'//
     &        ' in fb_reg_par for fb_type=  6, 7, 10-13 and 15')
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
            call do_faces(nFc, mpg, ifb, ifbreg)
            fbRegP(ifb,2) = ifbreg - fbRegP(ifb,1)
          endif
        ! now look for feedback types requiring volumes
          if ((fb_type(ifb).eq.1 .or. fb_type(ifb).eq.2 .or.  
     &        fb_type(ifb).eq.4  .or. fb_type(ifb).eq.5 .or.
     &        fb_type(ifb).eq.8) .and. .not.done) then
            call xertst(fb_reg_par(ifb,1).gt.0,
     &        'init_feedback: At lest one segment must be specified'//
     &        ' in fb_reg_par for fb_type=  1, 2, 5 and 8')
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

      subroutine compute_feedback(nCv, nFc, nVx, ns, ismain,
     &  switch, geo, mpg, pl, dv, rt, srw, psnc, psnl, main_call,
     &  diagno)
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, diagno
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2PlasmaSnapshot), intent (inout) :: psnc, psnl
      logical main_call
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
*   ..local variables
      integer :: ifb, ifc, ifc1, ifc2, icv, ic1, ic2, is, is_start,
     & is_end, iss
#ifdef B25_EIRENE
      integer :: i, j, k, l
#endif
      integer, save :: ncall = 0
      real(kind=r8) :: sum_e, total_flux, vol_feedback, sum_s

      call subini ('compute_feedback')
      if (ncall.eq.0) then
*   ..test nCv, nFc, ns
        call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
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
     &       fb_rescale_option(ifb) .gt. 0 .and.
     &       fb_actuator(ifb) .gt.0) then

            if (main_call) then
              fb_prev(ifb) = saved_fb_actuator(ifb)
            end if

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
     &                 pl%na(fbReg(icv),is) * geo%cvVol(fbReg(icv))
                    if(is.eq.is_start) then
                      vol_feedback = vol_feedback +
     &                               geo%cvVol(fbReg(icv))
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
     &            dv%ne(fbReg(icv)) * geo%cvVol(fbReg(icv))
                vol_feedback = vol_feedback +
     &                         geo%cvVol(fbReg(icv))
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

              fb_current(ifb) = 0.5_R8 *
     *          (dv%ne(omp(icsepomp-1))+dv%ne(omp(icsepomp)))
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
                    if(mpg%cvReg(fbReg(icv)).ne.0 .and.
     &                 fbReg(icv).le.mpg%nCi) then
                      fb_current(ifb) = fb_current(ifb) +
     &                  pl%na(fbReg(icv),is) * geo%cvVol(fbReg(icv))
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
                    if(mpg%cvReg(fbReg(icv)).ne.0 .and.
     &                 fbReg(icv).le.mpg%nCi) then
                      fb_current(ifb) = fb_current(ifb) +
     &                  pl%na(fbReg(icv),is) * geo%cvVol(fbReg(icv))
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

#ifdef B25_EIRENE
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
#else
              call xerrab('fb_type=6 not yet adapted for'//
     &          ' pure fluid cases')
#endif
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
                  if (mpg%cvOnClosedSurface(icv)) then
                    sum_s = sum_s + pl%na(icv,is) * geo%cvVol(icv)
                    charge_frac(iss) = charge_frac(iss) +
     .                pl%na(icv,is) * geo%cvVol(icv)
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
                  ifc1 = mpg%fsFcP(ifssep,1)
                  ifc2 = ifc1 + mpg%fsFcP(ifssep,2) - 1
! csc Here the idea is to follow the separatrix along its flux surface, then use the two CVs
!     belonging to each face of the flux surface
                  do iFc = ifc1, ifc2
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
                  enddo
                  is_end = is
                endif
              enddo
              fb_current(ifb) = fb_current(ifb)/sum_e

            case (8) ! relative average concentration of this species on desired domain

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &            'compute_feedback_type: using average ',
     &            'concentration for sequence ', 
     &            b2espcr(fb_species(ifb))   
              endif
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
              ifc1 = fbRegP(ifb,1)
              ifc2 = ifc1 + fbRegP(ifb,2) - 1
              do iFc = ifc1, ifc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fna(fbReg(ifc),1,fb_species(ifb))
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
#ifdef COMPLICATED
              call xerrab('COMPLICATED not yet converted for WG')
! the 'complicated' calculation downwards was correct assuming we look anly at core boundary
! while now we generalize for any boundary/surface 
              teave1=0.0_R8
              neave1=0.0_R8
              teave2=0.0_R8
              neave2=0.0_R8
#endif
              fb_current(ifb) = 0.0_R8
              ifc1 = fbRegP(ifb,1)
              ifc2 = ifc1 + fbRegP(ifb,2) - 1
              do iFc = ifc1, ifc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fhe(fbReg(ifc),1)
#ifdef COMPLICATED
            ic1 = mpg%bcCvP(ib,1)
            ic2 = ic1 + mpg%bcCvP(ib,2) -1
            do icv = ic1, ic2
              teave1=teave1+pl%te(mpg%bcCv(icv,1))*
     .                      dv%ne(mpg%bcCv(icv,1))*
     .                     geo%cvVol(mpg%bcCv(icv,1))
              neave1=neave1+dv%ne(mpg%bcCv(icv,1))*
     .                     geo%cvVol(mpg%bcCv(icv,1))
              teave2=teave2+pl%te(mpg%bcCv(icv,2))*
     .                      dv%ne(mpg%bcCv(icv,2))*
     .                     geo%cvVol(mpg%bcCv(icv,2))
              neave2=neave2+dv%ne(mpg%bcCv(icv,2))*
     .                     geo%cvVol(mpg%bcCv(icv,2))
            enddo
#endif
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
#ifdef COMPLICATED
              call xerrab('COMPLICATED not yet converted for WG')
! the 'complicated' calculation downwards was correct assuming we look anly at core boundary
! while now we generalize for any boundary/surface              
              tiave1=0.0_R8
              niave1=0.0_R8
              tiave2=0.0_R8
              niave2=0.0_R8
#endif
              fb_current(ifb) = 0.0_R8
              ifc1 = fbRegP(ifb,1)
              ifc2 = ifc1 + fbRegP(ifb,2) - 1
              do iFc = ifc1, ifc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fhi(fbReg(ifc),1)
#ifdef COMPLICATED
            ic1 = mpg%bcCvP(ib,1)
            ic2 = ic1 + mpg%bcCvP(ib,2) -1
            do icv = ic1, ic2
              tiave1=tiave1+pl%ti(mpg%bcCv(icv,1))*
     .                      dv%ni(mpg%bcCv(icv,1),0)*
     .                     geo%cvVol(mpg%bcCv(icv,1))
              niave1=niave1+dv%ni(mpg%bcCv(icv,1),0)*
     .                     geo%cvVol(mpg%bcCv(icv,1))
              tiave2=tiave2+pl%ti(mpg%bcCv(icv,2))*
     .                      dv%ni(mpg%bcCv(icv,2)),0*
     .                     geo%cvVol(mpg%bcCv(icv,2))
              niave2=niave2+dv%ni(mpg%bcCv(icv,2),0)*
     .                     geo%cvVol(mpg%bcCv(icv,2))
            enddo
#endif
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
              ifc1 = fbRegP(ifb,1)
              ifc2 = ifc1 + fbRegP(ifb,2) - 1
              do iFc = ifc1, ifc2
                fb_current(ifb)=fb_current(ifb)+
     &           dv%fch(fbReg(ifc),1)
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
                  do icv = 1, nCv
                   fb_current(ifb) = fb_current(ifb) +
     &               pl%na(icv,is)*geo%cvVol(icv)
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

            case default
               write(*,*) 'fb_type not coded ',
     &              fb_type(ifb)
               call xerrab('fb_type not coded')

            end select

! phase 2:

            select case(fb_rescale_option(ifb))

            case (0) ! no action

            case (1) ! rescale slowed by na_feedback_alpha

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
                if(fb_const(ifb).lt.0.0_R8) then
                  fb_const(ifb) = 0.0_R8
                  is_start=-1
                  is_end=-1
                  do is=0, ns-1
                    if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                      if(is_start.eq.-1) is_start=is
                      ic1 = fbRegP(ifb,1)
                      ic2 = ic1 + fbRegP(ifb,2) - 1
                      do iCv = ic1, ic2
                        if(mpg%cvReg(fbReg(icv)).ne.0 .and.
     &                     fbReg(icv).le.mpg%nCi) then
                          fb_const(ifb) = 
     &                      fb_const(ifb) + 
     &                      pl%na(fbReg(icv),is) * geo%cvVol(fbReg(icv))
                        endif
                      enddo
                      is_end=is
                    endif
                  enddo
                endif
                  fb_prev(ifb) = fb_current(ifb)
              endif

              fb_target(ifb) =
     .          fb_const(ifb) +
     .          fb_beta(ifb)*
     .          (tim-fb_time(ifb))
              fb_rescale =
     .          fb_alpha(ifb) *
     .          ((fb_target(ifb) -
     .          fb_current(ifb))/ddtim +
     .          (fb_prev(ifb) -
     .          fb_current(ifb))/dt_prev)
              saved_fb_actuator(ifb) = fb_current(ifb)

            case (5) ! rescale done according to SOLPS4 formula:
                        ! Nt = A*exp((time-T)*B)

              if(ncall.eq.0) then
                write(*,*)
     &           'compute_feedback_option: using ',
     &           'SOLPS4 exponential decay for sequence ',
     &           b2espcr(fb_species(ifb))
                if(fb_const(ifb).lt.0.0_R8) then
                  fb_const(ifb) = 0.0_R8
                  is_start=-1
                  is_end=-1
                  do is=0, ns-1
                    if(b2espcr(is).eq.b2espcr(fb_species(ifb))) then
                      if(is_start.eq.-1) is_start=is
                      do iCv = 1, mpg%nCi
                        if(mpg%cvReg(iCv).ne.0) then
                          fb_const(ifb) =
     &                      fb_const(ifb) +
     &                      pl%na(iCv,is) * geo%cvVol(iCv)
                        endif
                      enddo
                      is_end=is
                    endif
                  enddo
                endif
                fb_prev(ifb) = fb_current(ifb)
              endif

              fb_target(ifb) =
     .           fb_const(ifb) *
     .           exp(fb_beta(ifb)*
     .           (tim-fb_time(ifb)))
              fb_rescale(ifb) =
     .           fb_alpha(ifb) *
     .           ((fb_target(ifb) -
     .           fb_current(ifb))/ddtim +
     .           (fb_prev(ifb) -
     .           fb_current(ifb))/dt_prev)
              saved_fb_actuator(ifb) = fb_current(ifb)

            case (6) ! rescale slowed by fb_alpha (SOLPS4 style)

                  fb_rescale(ifb) =
     &               fb_alpha(ifb) * fb_current(ifb)

            case default
              write(*,*) 'fb_rescale_option not coded ',
     &          fb_rescale_option(ifb)

            end select

            write(*,'(a,1x,1p,i3,6g15.6)')
     &        'feedback_rescale ', b2espcr(fb_species(ifb)), 
     &        fb_current(ifb),
     &        fb_target(ifb),
     &        fb_rescale(ifb),
     &        fb_prev(ifb),
     &        ddtim, dt_prev

! phase 3:

            select case(fb_actuator(ifb))

            case (0) ! no action

            case (1) ! neutral particle flux => gas puff

              if(ncall.eq.0) then
                write(*,'(a,a,i3)')
     &           'compute_feedback_actuator: using ',
     &           'gas_puff for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              if(ncall.eq.0 .and. saved_fb_actuator(ifb).eq.0.0_R8) then
                if (conpar(is_start, fb_ib(ifb), 1).gt.0.0_R8) then
                  saved_fb_actuator(ifb) =
     &              conpar(is_start, fb_ib(ifb), 1)
                  write(message,'(a,i3,a,i2,a)')
     .             'Initial gas puff strength set by CONPAR(',
     .             is_start,',',fb_ib(ifb),',1) !'
                elseif (userfluxparm(fb_istra(ifb),1).gt.0.0_R8) then
                  saved_fb_actuator(ifb) =
     &              userfluxparm(fb_istra(ifb),1)
                  write(message,'(a,i3,a,i2,a)')
     .             'Initial gas puff strength set by USERFLUXPARM(',
     .             fb_istra(ifb),',1)  !'
                else
                  saved_fb_actuator(ifb) = 1.0e18_R8
                    write(*,*) 'Warning!  ','The gas puff was ',
     .               'started with an initial rate of 1e18!'
                endif
              endif

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &            saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &            fb_rescale(ifb)
              end if
              if(fb_puff_min(ifb) .ge. 0.0_R8)
     &          saved_fb_actuator(ifb) =
     &            max(fb_puff_min(ifb),
     &            saved_fb_actuator(ifb))
              if(fb_puff_max(ifb) .gt. 0.0_R8)
     &          saved_fb_actuator(ifb) =
     &            min(fb_puff_max(ifb),
     &            saved_fb_actuator(ifb))

              if(fb_overshoot(ifb) .gt. 1.0_R8 .and.
     &           fb_current(ifb) .gt.
     &           fb_target(ifb) *
     &           fb_overshoot(ifb)) then
! csc Here the idea is to change conpar if fluid neutrals are used, otherwise userfluxparm
!     This way it is easier to generalize the schemes when hybrid neutrals are used (to be done!)
                if (switch%use_eirene.eq.0) then
                  conpar(is_start, fb_ib(ifb), 1) = 0.0_R8
                else
                  userfluxparm(fb_istra(ifb),1) = 0.0_R8
                endif
              else
                if (switch%use_eirene.eq.0) then
                  conpar(is_start, fb_ib(ifb), 1) = 
     &              saved_fb_actuator(ifb)
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
     &           'compute_feedback_actuator: using ',
     &           'boundary flux condition for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              is_start = eb2spcr(b2espcr(fb_species(ifb)))+1
              is_end = eb2spcr(b2espcr(fb_species(ifb))) +
     &          nfluids(b2espcr(fb_species(ifb)))

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_fb_actuator(ifb) .eq. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             conpar(is_end, fb_ib(ifb), 1)

              if (fb_type(ifb).eq.6) then !csc, this is for backward compatibility, can probably be removed
                saved_fb_actuator(ifb) =
     &             (1.0_R8 - fb_alpha(ifb)) *
     &             saved_fb_actuator(ifb) +
     &             fb_rescale(ifb)
              else
                if (fb_rescale_option(ifb).lt.4) then
                  saved_fb_actuator(ifb) =
     &              saved_fb_actuator(ifb) * fb_rescale(ifb)
                else
                  saved_fb_actuator(ifb) =
     &              fb_rescale(ifb)
                endif
              end if

              if(fb_puff_min(ifb) .ge. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             max(fb_puff_min(ifb),
     &             saved_fb_actuator(ifb))
              if(fb_puff_max(ifb) .gt. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             min(fb_puff_max(ifb),
     &             saved_fb_actuator(ifb))

              if(fb_overshoot(ifb) .gt. 1.0_R8 .and.
     &             fb_current(ifb) .gt.
     &             fb_target(ifb) *
     &             fb_overshoot(ifb)) then
                 total_flux = 0.0_R8
              else
                 total_flux = saved_fb_actuator(ifb)
              endif

              if (fb_type(ifb).eq.6) then
                do is = is_start, is_end
                  iss = is-is_start+1
                  conpar(is, fb_ib(ifb), 1) =
     &              total_flux * charge_frac(iss)
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

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_fb_actuator(ifb) .eq. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             enepar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &            saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &            fb_rescale(ifb)
              end if

              enepar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (5) ! ion heat flux

              if(ncall.eq.0) then
                     write(*,'(a,a)')
     &           'compute_feedback_actuator: using ',
     &           'ion heat flux boundary condition'
              endif

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_fb_actuator(ifb) .eq. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             enipar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &            saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &            fb_rescale(ifb)
              end if

              enipar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (6) ! current

              if(ncall.eq.0) then
                     write(*,'(a,a)')
     &           'compute_feedback_actuator: using ',
     &           'electric current boundary condition'
              endif

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_fb_actuator(ifb) .eq. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             potpar(fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &            saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &            fb_rescale(ifb)
              end if

              potpar(fb_ib(ifb), 1) = saved_fb_actuator(ifb)

            case (7) ! density

              if(ncall.eq.0) then
                     write(*,'(a,a,i3)')
     &           'compute_feedback_actuator: using ',
     &           'density condition for sequence ',
     &           b2espcr(fb_species(ifb))
              endif

              if(feedback_namelist_used .and.
     &             ncall .eq. 0 .and.
     &             saved_fb_actuator(ifb) .eq. 0.0_R8)
     &             saved_fb_actuator(ifb) =
     &             conpar(fb_species(ifb),fb_ib(ifb), 1)

              if (fb_rescale_option(ifb).lt.4) then
                saved_fb_actuator(ifb) =
     &            saved_fb_actuator(ifb) * fb_rescale(ifb)
              else
                saved_fb_actuator(ifb) =
     &            fb_rescale(ifb)
              end if

              conpar(fb_species(ifb),fb_ib(ifb), 1) = 
     &          saved_fb_actuator(ifb)

            case default
              write(*,*) 'fb_actuator not coded ',
     &          fb_type(ifb)

            end select

            if (fb_actuator(ifb).eq.1 .or.
     &          fb_actuator(ifb).eq.3)
     >        write(*,'(a,1p,i3,g15.6)') 'feedback_actuator ',
     &         b2espcr(fb_species(ifb)), saved_fb_actuator(ifb)

          endif

          feedback_namelist_used = feedback_namelist_used .or.
     &     fb_actuator(ifb).eq.1 .or. fb_actuator(ifb).eq.3

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
      integer ifc

      do ifc = 1, nFc
        if ((mpg%fcReg(ifc).eq.fb_reg_par(ifb,1)) .or.
     &      (fb_reg_par(ifb,2).gt.0.and.
     &      (mpg%fcReg(ifc).eq.fb_reg_par(ifb,2)))) then
          fbReg(ifbreg) = ifc
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
      integer icv

      do icv = 1, nCv
        if ((mpg%cvReg(icv).eq.fb_reg_par(ifb,1)) .or.
     &      (fb_reg_par(ifb,2).gt.0.and.
     &      (mpg%cvReg(icv).eq.fb_reg_par(ifb,2)))) then
          fbReg(ifbreg) = icv
          ifbreg = ifbreg + 1
         endif
      end do

      return
      end subroutine do_volumes

***********************************************************************

      subroutine conv_fbparms (nx, ny, ns, m)
      use b2mod_feedback
      implicit none
      integer, intent(in) :: nx, ny, ns
      type (mapping), intent (inout) :: m
      integer :: is, ix, iy, icv, iss, ifbreg
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
      fb_puff_min(1:nfb) = na_feedback_puff_min(0:nspecies-1)
      fb_puff_max(1:nfb) = na_feedback_puff_max(0:nspecies-1)
      fb_overshoot(1:nfb) = na_feedback_overshoot(0:nspecies-1)
      ifbreg = 1
      do is = 0,nspecies-1
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
              icv = m%imapCv(ix,iy)
              if (icv.eq.0) cycle
              fb_Reg(ifbreg) = icv
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

      subroutine write_fbparms (ns)
      use b2mod_version
      implicit none
      integer, intent(in) :: ns
      integer :: ifb, iout
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
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

      end module b2us_feedback
