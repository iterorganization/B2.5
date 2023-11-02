      module b2mod_usrspb
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_geo
      use b2mod_numerics_namelist
      use b2mod_user_namelist
      use b2mod_b2plot_wall_loading
      use b2mod_dimensions
      implicit none

c
c  version : 31.08.2023 by Ilya Senichenkov
c
c----------------------------------------------------------------------
c*** This tracing file is alternative to user.trc - it contains different traced quantities and
c*** is more straightforvard in programming and processing

!    Q_pneut            - neutral pressure in PFR in Pa - mean (the same as in user.trc)
!    Q_pneut_trg(:)     - neutral pressure in PFR in Pa on two targets (the same as in user.trc)
!    Q_Te_min_LCFS      - minimal Te (in eV) on the outermost closed flux tube - of interest in XPR regime
!    Q_Te_max(:)        - peak Te on each target, in eV
!    Q_det_wdth(:)      - on each target, distance from the strike point to the point where Te < Te_det_threshold (2.0 eV by default), or 0.0 if Te_sep > Te_det_threshold.
!    Area_pneut_trg(:)  - auxiliary variable which means the area (in m^2) of those surfaces (on each of two targets) where the neutral pressure is computed
!    Q_I_tot_trg(:)     - on each target, the net electric current (in A, with sign) collected by target
!    Q_jsat_max_trg(:)  - on each target, peak saturation current density (in A/m^2, with sign)
!    Q_Isat_trg(:)      - on each target, the integral of saturation current density along the target (in A, with sign)
!    Q_qmax_trg(:)      - on each target, the peak energy flux density with account of volumetric radiation
!                         i.e. is the same as computed by `echo "1101 wlld" | b2plot`, AND DIFFERS FROM WHAT IS WRITTEN TO user.trc

c----------------------------------------------------------------------
c
c*** lv           :  version number to appear in the file
c*** N_Q_global:  number of output quantities
c*** N_Q_target:   number of output per each target
c
      integer lv, N_Q_global, N_Q_target  !iyv 11.04.22
      parameter (lv=20230831, N_Q_global=6, N_Q_target=6)  ! IYS 22.03.2022 iyv 11.04.22

      integer, save :: N_Q_total
      integer ltrc
      parameter (ltrc=ntrc_uss)

c----------------------------------------------------------------------
c*** Output data (for each species):
c***
c----------------------------------------------------------------------

      integer, save :: lun, mintopcut
      real (kind=R8), save :: Q_pneut, Q_Te_min_LCFS, 
     &    Q_rad_tot, Q_rad_SOL
      real (kind=R8), allocatable, save ::
     &    Q_pneut_trg(:),
     &    Q_Te_max(:), Q_det_wdth(:), Area_pneut_trg(:),
     &    Q_I_tot_trg(:), Q_jsat_trg(:,:), Q_Isat_trg(:),
     &    Q_qmax_trg(:), Q_jsat_max_trg(:)

      integer, allocatable, save :: ix_pneut_beg(:), ix_pneut_end(:)  ! range of cells to average neutral pressure over
      character*12, allocatable, save :: headers(:)
      character*9, save :: header_tmpl(N_Q_target)
c*** Column headers:
      data header_tmpl /
     1    '   Te_max', ' det_wdth', '    q_max',
     2    '  j_s_max','  I_s_tot', '    I_tot'/

      character, save :: hlp_frm*256
      

      public ::     calc_pneut_init, calc_pneut

      contains

      subroutine b2usstrc(nx,ny,ns,BoRiS)
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_plasma
      use b2mod_rates
      use b2mod_sources
      use b2mod_indirect
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2mod_dimensions

c----------------------------------------------------------------------
c*** debug_flags usage:
c***  90  >0 : debugging printout on
c***  90, 91 : 'non-default surfaces' (90) to (91)
c***  92, 93 : strata (92) to (93); 0 means the sum over all strata
c***  94, 95 : more detailed output for strata (94) to (95)
c----------------------------------------------------------------------

      implicit none

      integer nx, ny, ns

      integer it  ! counter for targets
      integer ix, iy, is, k, l, j, ixt, ixt_fl, iyp                 !iyv 11.04.22

      logical first

      real (kind=R8) :: BoRiS

      real (kind=R8) :: u, t0

c----------------------------------------------------------------------
c*** Taken from the common blocks:
c***
c*** lcbs: boundary segments representing the core boundary in B2.5
c*** lcns: 'non-default surfaces' for the core boundary (Eirene)
c*** ltns: 'non-default surfaces' for the targets (Eirene)
c*** lsns: list of Eirene surfaces comprising the recycling sources
c*** msns: msns(1,k): index of first surface in lsns list for Eirene
c***                  recycling source k
c***       msns(2,k): number of surfaces for Eirene recycling source k
c*** mlcmp: composition of molecules.
c***        Molecule of kind j contains mlcmp(i,j) atoms of kind i
c*** gpfc: composition of gas puff sources.
c***       gpfc(i,j)=(fraction of atoms of kind i in stratum number j)
c*** sclstra: scaling factors for species and strata
c*** lmolscl: pointer to the species with which these molecules scale
c----------------------------------------------------------------------
*   ..initialisation
      data first / .true. /

c----------------------------------------------------------------------
      save first
c----------------------------------------------------------------------

      if(first_time_step) return

      if(first) then
       first=.false.
       from_b2plot = .false.
      end if

! Calculate the peak power loads (with account of radiation)
      irdparm=1
      call b2ptrgl(BoRiS)

! neutral pressure in PFR
! Q_pneut - mean, Q_pneut_trg(1) - inner PFR, Q_pneut_trg(2) - outer PFR
      call calc_pneut(ns, ny)

      l=0
      u=-1.0_R8
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        it=itrg(j)
        if (it.gt.0) then !{
          ixt   =ltrgpos(it)+l ! narrow cells for temperatures
          ixt_fl=ltrgpos(it)+1 ! for currents                             !} iyv 11.04.22

! Maximal electron temperature along target
          Q_Te_max(it) = maxval(te(ixt,:)) / eV

c------------------------------
! So called width of detached target
! It is measured along the target from the strike point
! (lseppos(it) points to first cell in the PFR at given target)
! up to the point where Te rises above the prescribed threshold
          if (te(topix(ixt,lseppos(it)),topiy(ixt,lseppos(it)))
     &                    / eV .le. Te_det_threshold) then

            Q_det_wdth(it) = 0.5_R8 *
     &        hy(topix(ixt,lseppos(it)),topiy(ixt,lseppos(it)))
            iy = topiy(ixt,topiy(ixt,lseppos(it)))

            do while (iy .lt. ny .and. iy.gt.lseppos(it)+1 .and.
     &           te(ixt,iy) / eV .le. Te_det_threshold)
              Q_det_wdth(it) = Q_det_wdth(it) +
     &          0.5_R8 * (hy(bottomix(ixt,iy),bottomiy(ixt,iy)) +
     &                    hy(ixt,iy))
              iy = topiy(ixt,iy)
            enddo

! Linear interpolation
! Only perform it the temperature is increasing radially outwards
            if (te(ixt,iy).gt.te(bottomix(ixt,iy),bottomiy(ixt,iy)))
     .       then
              Q_det_wdth(it) = Q_det_wdth(it) + 0.5_R8 *
     &          (hy(bottomix(ixt,iy),bottomiy(ixt,iy))+hy(ixt,iy))
     &        / (te(ixt,iy)-te(bottomix(ixt,iy),bottomiy(ixt,iy)))
     &        * (Te_det_threshold * eV -
     &           te(bottomix(ixt,iy),bottomiy(ixt,iy)))
            else if (iy.eq.ny) then
! We have reached the edge of the grid and the temperature is decreasing
! So we only add the width of the cell, since we do not know if or where
! the profile will rise again
              Q_det_wdth(it) = Q_det_wdth(it) +
     &         0.5_R8 * hy(bottomix(ixt,ny),bottomiy(ixt,ny))
            endif

          else
! Width of detached target is zero if strike point temperature is
! above the threshold
            Q_det_wdth(it) = 0.0_R8
          endif
c------------------------------
! Total poloidal electric current flowing between plasma and targets
          Q_I_tot_trg(it) = sum(fch(ixt_fl,0:ny-1,0))

c------------------------------
! Saturation current density to targets and net saturation current
          Q_jsat_trg(-1:ny,it) = 0.0_R8

          do is = 0, ns-1
            Q_jsat_trg(-1:ny,it) = Q_jsat_trg(-1:ny,it) +
     &        u * fna(ixt_fl,-1:ny,0,is) * rza(ixt_fl,-1:ny,is)
!  factor 'u' is used here to make saturation current positive on Eastern and Western targets
          enddo

          Q_Isat_trg(it) = sum(Q_jsat_trg(-1:ny,it)) * qe

          Q_jsat_trg(-1:ny,it) = Q_jsat_trg(-1:ny,it) /
     &             gs(ixt_fl,-1:ny,0) * qe

c------------------------------
! Maximal saturation current density to targets
          Q_jsat_max_trg(it) = maxval(Q_jsat_trg(:,it))

c------------------------------
! Maximal energy loads on targets
          Q_qmax_trg(it) = wpk_tt(it)
        end if !}
        l=1-l
        u=-u
      enddo

c------------------------------
! Minimal Te on last closed flux surface
!        iy_te_min = mintopcut-1
!        ix_te_min = nmdpl  ! ix_te_min = nx/2
      Q_Te_min_LCFS = te(nmdpl,mintopcut-1)

      do ix=-1,nx
        if(on_closed_surface(ix,mintopcut-1)) then
          if(te(ix,mintopcut-1) .lt. Q_Te_min_LCFS) then
            Q_Te_min_LCFS = te(ix,mintopcut-1)
!                ix_te_min = ix
          endif
        endif
      enddo

      Q_Te_min_LCFS = Q_Te_min_LCFS / eV
      
c------------------------------
! Total radiated power and radiated power outside the separatrix
      Q_rad_tot =  0.0_R8
      Q_rad_SOL =  0.0_R8
       
      do ix=0,nx-1
        do iy=0,ny-1
          t0 = 0.0_R8
          do is=0,ns-1
            t0 = t0 + rqbrm(ix,iy,is) + rqrad(ix,iy,is)
          enddo
#ifdef B25_EIRENE
          t0 = t0 - sum(eneutrad(ix+1,iy+1,1:nnatmi,0)) - 
     &        sum(emolrad(ix+1,iy+1,1:nnmoli,0)) -
     &           sum(eionrad(ix+1,iy+1,1:nnioni,0))              
#endif
          
          Q_rad_tot = Q_rad_tot + t0  
          if (.not. ON_CLOSED_SURFACE(ix,iy)) Q_rad_SOL = Q_rad_SOL + t0
        enddo
      enddo


c*** Write the data in a separate file for plotting...
c
      write(lun,hlp_frm) tim,
     &     Q_pneut,
     &    (Q_pneut_trg(it), it=1,2),
     &     Q_Te_min_LCFS, Q_rad_tot, Q_rad_SOL,
     &    (Q_Te_max(it), it=1,ntrgts),
     &    (Q_det_wdth(it), it=1,ntrgts),
     &    (Q_qmax_trg(it), it=1,ntrgts),         !iyv 11.04.22
     &    (Q_jsat_max_trg(it), it=1,ntrgts),
     &    (Q_Isat_trg(it), it=1,ntrgts),
     &    (Q_I_tot_trg(it), it=1,ntrgts)

      return
      end subroutine b2usstrc

c======================================================================
c======================================================================
c======================================================================
c======================================================================

      subroutine b2usstri(nx,ny,nxtl,nxtr,jsep)
      use b2mod_b2plot
     & , only : nx_b2plot => nx, ny_b2plot => ny, jsep_b2plot => jsep,
     &          nxtl_b2plot => nxtl, nxtr_b2plot => nxtr
      implicit none
      integer nx, ny, nxtl, nxtr, jsep
      integer i, j, nv, N_ttl_read, N_per_trg_read, Ntrgts_read
      real (kind=R8), save :: stim = 0.0_R8
      character f*(lenfnm), hlp_c
      logical :: ex

c----------------------------------------------------------------------
c*** Initialise the file for tracing data
c*** Only called with filedata(ltrc)=.true.
c----------------------------------------------------------------------

      if(first_time_step) then
        call alloc_b2plot_wall_loading(nlim,nsgmx)
        from_b2plot = .false.
        nx_b2plot = nx
        ny_b2plot = ny
        nxtr_b2plot = nxtr
        nxtl_b2plot = nxtl
        jsep_b2plot = jsep
      endif
c*** Input the DG data and prepare the geometry (this routine works once)
      call b2pwinp

      N_Q_total = N_Q_global + ntrgts * N_Q_target

      if (.not. allocated(Q_Te_max)) then
        allocate(
!     &   Q_pneut_trg(2), Area_pneut_trg(ntrgts),
     &   Q_Te_max(ntrgts), Q_det_wdth(ntrgts),
     &   Q_I_tot_trg(ntrgts), Q_jsat_trg(-1:ny,ntrgts),
     &   Q_jsat_max_trg(ntrgts), Q_Isat_trg(ntrgts),
     &   Q_qmax_trg(ntrgts))
!        allocate(
!     &   ix_pneut_beg(2), ix_pneut_end(2))
        allocate(headers(N_Q_total))
      endif

! The allocating of Q_pneut_trg, Area_pneut_trg, ix_pneut_beg, ix_pneut_end
! is moved into calc_pneut_init

      call calc_pneut_init(nx,ny,nxtl,nxtr)
      
      mintopcut = ny
      if(nncut > 0) mintopcut = minval(topcut(1:nncut))

      headers(1) = ' p_neut_PFR'
      headers(2) = 'p_neut_PFRi'
      headers(3) = 'p_neut_PFRo'
      headers(4) = 'Te_min_LCFS'
      headers(5) = '  Q_rad_tot'
      headers(6) = '  Q_rad_SOL'


      if (ntrgts .eq. 2) then
        do i=1,N_Q_target
          headers(N_Q_global+(i-1)*ntrgts+1) = header_tmpl(i)//'_i '  ! inner target
          headers(N_Q_global+(i-1)*ntrgts+2) = header_tmpl(i)//'_o '  ! outer target
        enddo
      elseif (ntrgts .eq. 4) then
        do i=1,N_Q_target
          headers(N_Q_global+(i-1)*ntrgts+1) = header_tmpl(i)//'_i '  ! inner target
          headers(N_Q_global+(i-1)*ntrgts+2) = header_tmpl(i)//'_iu'  ! outer target
          headers(N_Q_global+(i-1)*ntrgts+3) = header_tmpl(i)//'_ou'  ! inner upper target
          headers(N_Q_global+(i-1)*ntrgts+4) = header_tmpl(i)//'_o '  ! outer upper target
        enddo
      endif

      call ipgetr ('b2mndr_stim', stim)

      if(.not.filedata(ntrc_uss)) return

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
!      nc = 0
      lun = ntrc(ltrc)
      f = fn_ntrc(ltrc)

      call find_file(f,ex)

      open(lun, file=f, err=90)
      rewind(lun)

      if(ex) read(lun, * , end=91,err=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)

c----------------------------------------------------------------------
      if(ex .and. stim.lt.0.0_R8) then

        read(lun, *, err=40)
     &      hlp_c, nv, N_ttl_read, Ntrgts_read, N_per_trg_read
        write(*,*) hlp_c, nv, N_ttl_read, Ntrgts_read, N_per_trg_read
        read(lun, *) hlp_c
        go to 42
   40   read(lun,*,err=92) hlp_c, i, j

        lv_trc(ltrc)=lv

   42   continue

! Check version
        if(nv.ne.lv) then
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2mod_usrspb: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
          filedata(ltrc)=.false.
          return
        end if

! Check total number of columns
        if(N_ttl_read .ne. N_Q_total) then
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2mod_usrspb: current number of columns ',N_Q_total,
     ,     ' differs from that in the file ',trim(f),' : ',N_ttl_read
          write(*,*) 'Writing into file ',trim(f),
     ,   ' cannot be continued - remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if

! Check total number of targets
        if(N_per_trg_read .ne. N_Q_target) then
          write(*,'(a,i2,3a,i2)')
     ,   '!!! b2mod_usrspb: current number of quantities per targets ',
     ,      N_Q_target,' differs from that in the file ',trim(f),' : ',
     ,      N_per_trg_read
          write(*,*) 'Writing into file ',trim(f),
     ,   ' cannot be continued - remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if

! Check the number of quantities per target
        if(Ntrgts_read .ne. ntrgts) then
          write(*,'(a,i2,3a,i2)')
     ,   '!!! b2mod_usrspb: the number of targets ',
     ,      ntrgts,' differs from that in the file ',trim(f),' : ',
     ,      Ntrgts_read
          write(*,*) 'Writing into file ',trim(f),
     ,   ' cannot be continued - remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if

   50   read(lun,'(a1)',end=55,err=55) hlp_c
        go to 50
   55   backspace lun

c----------------------------------------------------------------------
      else ! file 'user_SPb.trc' does not yet exist

        if(.not.ex)
     >   write(*,*) 'File user_SPb.trc not present. Will be created...'

        write(lun, 700) lv, N_Q_total, ntrgts, N_Q_target, trim(label)

        write(hlp_frm,'(a,i3,a)')
     &   '(a1,6x,a4,3x,',N_Q_total,'(2x,a19))'

        write(lun, hlp_frm) '#','time',
     &   (headers(i), i=1, N_Q_total)

      end if

      write(hlp_frm,'(a,i3,a)')
     &  '(1p,e16.7,',N_Q_total,'(3x,e16.7,2x))'
      return

c----------------------------------------------------------------------
   90 write(*,*) 'b2mod_usrspb: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2mod_usrspb: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return

c----------------------------------------------------------------------
  700 format('#   ',4i8/'#  ',a)
c----------------------------------------------------------------------
      end subroutine b2usstri

c======================================================================
c======================================================================
c======================================================================
c======================================================================


      subroutine b2usstrcc(iline)
      implicit none
      integer, intent(out) :: iline
      real (kind=R8) :: tdummy
      character f*(lenfnm), hlp_c
      logical :: ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------

      if(.not.filedata(ltrc)) return

      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)

      open(lun,file=f,err=190)
      rewind(lun)

      iline = 0
      read(lun,'(a1)',end=155) hlp_c

      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',end=155) hlp_c
      end do

      backspace lun

      read(lun,*,end=155) tdummy

      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,end=155) tdummy
      end do

  155 continue
!      close(lun)
      return
c----------------------------------------------------------------------
  190 write(*,*) 'b2mod_usrspb: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
c----------------------------------------------------------------------
      end subroutine b2usstrcc
c
      subroutine b2usstrs
c----------------------------------------------------------------------
c*** - to be called BEFORE each timestep
c----------------------------------------------------------------------
c
      return
      end subroutine b2usstrs

c======================================================================
c======================================================================
c======================================================================
c======================================================================

      subroutine b2usstrf
      implicit none
c----------------------------------------------------------------------
c*** - finalization step
c----------------------------------------------------------------------
      logical file_open

      lun=ntrc(ltrc)
      inquire(unit=lun,opened=file_open)
      if (file_open) close(lun)

      if (allocated(Q_Te_max)) then
        deallocate(
     &   Q_Te_max, Q_det_wdth, Q_I_tot_trg,
     &   Q_jsat_trg, Q_Isat_trg, Q_qmax_trg)
        deallocate(headers)
       endif

      if (allocated(Q_pneut_trg)) then
        deallocate ( Q_pneut_trg)
        deallocate (ix_pneut_beg, ix_pneut_end)
      endif
      
      return
      end subroutine b2usstrf

c======================================================================
c======================================================================
c======================================================================
c======================================================================


      subroutine calc_pneut_init(nx,ny,nxtl,nxtr)
      implicit none
      
      integer nx, ny, nxtl, nxtr
      
! Subroutine returns ix_pneut_beg(1:2) and ix_pneut_end(1:2):
! poloidal cell number ranges over which neutral pressure will be averaged

! This subroutine is called either from b2usstri or b2stbc_fb_pre, so it is worth checking whether necessary arrays are allocated
      if (.not. allocated(Q_pneut_trg)) then
! More safety here since Q_pneut_trg may be allocated from b2stbc_fb
        allocate(Q_pneut_trg(2), Area_pneut_trg(2))
        allocate(ix_pneut_beg(2), ix_pneut_end(2))
      endif
  
! Western divertor
      if(lpfrs_pmp.ne.0) then !  iyv 11.04.22
c***  lpfrs_pmp - location of the pump: 0 - no pump at all,
c***              1 - in lower PFR, 2 - on the outer wall,
c***              3 - in upper PFR, 4 - on the inner wall

        ix_pneut_beg(1) = lpfrb_i !{
        if (ix_pneut_beg(1) .le. 0) then !{
          if (lpfrs_pmp .eq. 1) then !{
            ix_pneut_beg(1) = leftcut(1) - 1 + ix_pneut_beg(1)
          elseif (lpfrs_pmp .eq. 2) then !}{
            ix_pneut_beg(1) = nx - 1 + ix_pneut_beg(1)
          elseif(lpfrs_pmp .eq.3 .or. lpfrs_pmp .eq. 4) then !}{
            ix_pneut_beg(1) = nxtl - 1 + ix_pneut_beg(1)
          end if !}
        end if !}
        ix_pneut_end(1) = lpfrt_i
        if (ix_pneut_end(1) .le. 0) then !{
          if (lpfrs_pmp .eq. 1) then !{
            ix_pneut_end(1) = leftcut(1) - 1 + ix_pneut_end(1)
          elseif (lpfrs_pmp .eq. 2) then !}{
            ix_pneut_end(1) = nx - 1 + ix_pneut_end(1)
          elseif (lpfrs_pmp .eq. 3) then !}{
            ix_pneut_end(1) = leftcut(2) - ix_pneut_end(1)
          elseif(lpfrs_pmp.eq.4) then !}{
            ix_pneut_end(1) = nxtl - 1 + ix_pneut_end(1)
          end if !}
        end if !}

! Eastern divertor

        ix_pneut_beg(2) = lpfrt_o
        if (ix_pneut_beg(2) .le. 0) then !{
          if (lpfrs_pmp .eq. 1) then !{
            ix_pneut_beg(2) = rightcut(1) - ix_pneut_beg(2)
          elseif(lpfrs_pmp .eq. 2) then !}{
            ix_pneut_beg(2) = nx - 1 + ix_pneut_beg(2)
          elseif(lpfrs_pmp .eq. 3) then !}{
            ix_pneut_beg(2) = rightcut(2) - 1 + ix_pneut_beg(2)
          elseif(lpfrs_pmp.eq.4) then !}{
            ix_pneut_beg(2) = nxtr - 1 - ix_pneut_beg(2)
          end if !}
        end if !}
        ix_pneut_end(2) = lpfrb_o
        if (ix_pneut_end(2) .le. 0) then !{
          if(lpfrs_pmp .eq. 1 .or. lpfrs_pmp .eq. 2) then !{
            ix_pneut_end(2) = nx - 1 + ix_pneut_end(2)
          else if(lpfrs_pmp .eq. 3) then !}{
            ix_pneut_end(2) = nxtr + 1 - ix_pneut_end(2)
          else if(lpfrs_pmp .eq. 4) then !}{
            ix_pneut_end(2) = nxtl - 1 + ix_pneut_end(2)
          end if !}
        end if !}
      else !}{
        ix_pneut_beg(1) = 1
        ix_pneut_end(1) = 0
        ix_pneut_beg(2) = 1
        ix_pneut_end(2) = 0
      end if !}
      
      return
      end subroutine calc_pneut_init
c======================================================================

      subroutine calc_pneut(ns,ny)
      use b2mod_driver
     ,    , only : use_eirene
      use b2mod_plasma
     ,    , only : na, ti
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      
      integer ix,iy,iyp,it, ns,ny, k, is

! neutral pressure in PFR
! Q_pneut - mean, Q_pneut_trg(1) - inner PFR, Q_pneut_trg(2) - outer PFR

      iy=0
      iyp=1
      if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.4) then
        iy=ny
        iyp=ny
      end if

      do it = 1, 2
        Area_pneut_trg(it) = 0.0_R8
        Q_pneut_trg(it) = 0.0_R8

        do ix = ix_pneut_beg(it), ix_pneut_end(it)
          Area_pneut_trg(it) = Area_pneut_trg(it) + gs(ix,iy,1)

#ifdef B25_EIRENE
          if(use_eirene.ne.0) then
            do k=1,nnatmi
              Q_pneut_trg(it) = Q_pneut_trg(it) +
     &         dab2(ix+1,iyp,k,1)*tab2(ix+1,iyp,k,1) * gs(ix,iy,1)
            enddo
            do k=1,nnmoli
              Q_pneut_trg(it) = Q_pneut_trg(it) +
     &         dmb2(ix+1,iyp,k,1)*tmb2(ix+1,iyp,k,1) * gs(ix,iy,1)
            enddo
          else
            do is=0,ns-1
              if (is_neutral(is)) Q_pneut_trg(it) = Q_pneut_trg(it) +
     &          na(ix,iy,is)*ti(ix,iy) * gs(ix,iy,1)
            enddo
          endif
#else
          do is=0,ns-1
            if (is_neutral(is)) Q_pneut_trg(it) = Q_pneut_trg(it) +
     &        na(ix,iy,is)*ti(ix,iy) * gs(ix,iy,1)
          enddo
#endif
        enddo

        if(Area_pneut_trg(it).gt.0.0_R8)
     &   Q_pneut_trg(it) = Q_pneut_trg(it) / Area_pneut_trg(it)
      end do

      if(Area_pneut_trg(1).gt.0.0_R8 .and. Area_pneut_trg(2).gt.0.0_R8)
     > Q_pneut = ( Q_pneut_trg(1) * Area_pneut_trg(1) +
     &    Q_pneut_trg(2) * Area_pneut_trg(2) ) /
     &   (Area_pneut_trg(1) + Area_pneut_trg(2))

      return
      end subroutine calc_pneut

      end module b2mod_usrspb

