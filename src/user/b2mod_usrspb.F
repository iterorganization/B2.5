      module b2mod_usrspb
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_geo
      use b2mod_numerics_namelist
      use b2mod_user_namelist
      use b2mod_b2plot_wall_loading
      implicit none
#include <DIMENSIONS.F>

c
c  version : 06.21.2022 by Ilya Senichenkov
c
c----------------------------------------------------------------------
c*** Computes some quantities which are of interest for the SPb group

c----------------------------------------------------------------------
c
c*** lv           :  version number to appear in the file
c*** N_Q_global:  number of output quantities
c*** N_Q_target:   number of output per each target
c
      integer lv, N_Q_global, N_Q_target  !iyv 11.04.22
      parameter (lv=20221006, N_Q_global=4, N_Q_target=6)  ! IYS 22.03.2022 iyv 11.04.22
      
      integer, save :: N_Q_total
      integer ltrc
      parameter (ltrc=ntrc_uss)
      integer it  ! counter for targets


c----------------------------------------------------------------------
c*** Output data (for each species):
c***
c----------------------------------------------------------------------

      integer, save :: nc, lun, mintopcut
      real (kind=R8), allocatable, save :: 
     &          Q_pneut, Q_Te_min_LCFS, Q_pneut_trg(:),
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

      contains

      subroutine b2usstrc(nx,ny,ns,BoRis)
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
      use eirmod_clgin
     ,   , only : iliin
#endif

#include "DIMENSIONS.F"

c----------------------------------------------------------------------
c*** debug_flags usage:
c***  90  >0 : debugging printout on
c***  90, 91 : 'non-default surfaces' (90) to (91)
c***  92, 93 : strata (92) to (93); 0 means the sum over all strata
c***  94, 95 : more detailed output for strata (94) to (95)
c----------------------------------------------------------------------

      implicit none

      integer nx, ny, ns
      real (kind=R8) :: dtim

      integer idir, ix, iy, is, k, l, j, ixt, ixt_fl, iyp       !iyv 11.04.22

      logical first
      integer, save :: mdf_fnb = 0, style_b2mndt = 1                ! IYS 12.01.2018

      real (kind=R8) :: hlp_pls(nbcd,ns)
      integer use_recyceir, use_eirene
      real (kind=R8) :: recyctmp
      real (kind=R8) :: BoRiS


      real (kind=R8) :: u, hlpdt, hlp, hlp_rec, hlp_sput,
     &                        bsign, hlp_nw, hlp_pump
      logical
     &     hlp_pump_counted(nnatmi+nnmoli+nnioni+nnplsi,nnstsi+nnlimi)
       logical from_b2plot
     
      real (kind=R8) :: Te_det_threshold
      parameter (Te_det_threshold = 2.0_R8)

c----------------------------------------------------------------------
c
      character hlp_prt*26
c
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
      data use_recyceir /1/, use_eirene /0/
c
      data first / .true. /
      data from_b2plot / .false. /
      
c----------------------------------------------------------------------
      save !i, j, k, l, m, n, ic, is, ix, iy,
!     1 js, ks, ls, kks, 
     2 use_recyceir, use_eirene, first
c----------------------------------------------------------------------

      if(first_time_step) return
      
      if(first) then
       call ipgeti ('b2mndr_eirene', use_eirene)   
       first=.false.
      end if

! Calculate the peak power loads (with account of radiation)      
      irdparm=1
      call b2ptrgl(BoRiS, from_b2plot, use_eirene)

! neutral pressure in PFR 
! Q_pneut - mean, Q_pneut_trg(1) - inner PFR, Q_pneut_trg(2) - outer PFR
      do it = 1, 2
        iy=0
        iyp=1
        if(lpfrs_pmp.eq.2 .or. lpfrs_pmp.eq.4) then
          iy = ny
          iyp=ny
        end if  
        Area_pneut_trg(it) = 0.0_R8
        Q_pneut_trg(it) = 0.0_R8
      
          do ix = ix_pneut_beg(it), ix_pneut_end(it)
            Area_pneut_trg(it) = Area_pneut_trg(it) + gs(ix,iy,1)

            do k=1,nnatmi
              Q_pneut_trg(it) = Q_pneut_trg(it) + 
     &        dab2(ix+1,iyp,k,1)*tab2(ix+1,iyp,k,1) * gs(ix,iy,1)
            enddo
            do k=1,nnmoli
              Q_pneut_trg(it) = Q_pneut_trg(it) + 
     &        dmb2(ix+1,iyp,k,1)*tmb2(ix+1,iyp,k,1) * gs(ix,iy,1)
            enddo
          enddo
        
          Q_pneut_trg(it) = Q_pneut_trg(it) / Area_pneut_trg(it) 
      end do    

          Q_pneut = ( Q_pneut_trg(1) * Area_pneut_trg(1) + 
     &    Q_pneut_trg(2) * Area_pneut_trg(2) ) / 
     &     (Area_pneut_trg(1) + Area_pneut_trg(2))  
          
!          u = -1.0_R8    ! for Western targets x-direction is opposite to target normal
!          u = 1.0_R8     ! for Eastern targets x-direction and normal direction coincide

      l=0
      u=-1.0_R8
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        it=itrg(j) 
        if (it.gt.0) then !{
          ixt   =ltrgpos(it)+l ! narrow cells for temperatures
          ixt_fl=ltrgpos(it)+1 ! for currents                             !} iyv 11.04.22 

! Maximal electron temperature along target
          Q_Te_max(it) = maxval(te(ixt,:)) / qe

c------------------------------
! So called width of detached target
! It is measured along the target from the strike point 
! (nsep(it) should point to first cell in the SOL at given target)
! up to the point where Te rises above the prescribed threshold
          if (te(ixt,lseppos(it)) / qe .le. Te_det_threshold) then 

            Q_det_wdth(it) = 0.5_R8 * hy(ixt,lseppos(it))
            iy = lseppos(it) + 1

            do while (iy .lt. ny .and. 
     &           te(ixt,iy) / qe .le. Te_det_threshold)
              Q_det_wdth(it) = Q_det_wdth(it) + 
     &          0.5_R8 * (hy(ixt,iy-1) + hy(ixt,iy)) 

              iy = iy + 1                        
            enddo
 
          else
! Width of detached target is zero if strike point temperature is
! above the threshold
            Q_det_wdth(it) = 0.0_R8
          endif
c------------------------------
! Net electric current to targets
          Q_I_tot_trg(it) = sum(fch(ixt_fl,0:ny-1,0))

c------------------------------
! Saturation current density to targets and net saturation current
          Q_jsat_trg(-1:ny,it) = 0.0_R8

          do is=0, ns
            Q_jsat_trg(-1:ny,it) = Q_jsat_trg(-1:ny,it) +
     &        u * fna(ixt_fl,-1:ny,0,is) * rza(ixt_fl,-1:ny,is)
!  factor 'u' is used here to make saturation current positive on Eastern and Western targets
          enddo

          Q_Isat_trg(it) = sum(Q_jsat_trg(-1:ny,it)) * qe  

          Q_jsat_trg(-1:ny,it) = Q_jsat_trg(-1:ny,it) / 
     &             hy(ixt,-1:ny) / hz(ixt,-1:ny) * qe

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
              Q_Te_min_LCFS    = te(ix,mintopcut-1)
!                ix_te_min = ix
            endif
          endif
      enddo

      Q_Te_min_LCFS = Q_Te_min_LCFS / qe

c*** Write the data in a separate file for plotting...
c
c?        if(nc.gt.np_tmp .or. nc.le.0) then
c?          write(*,*) 'b2mod_usrspb: unexpected error. nc =',nc
c?          filedata(ltrc)=.false.
c?        else
          write(lun,hlp_frm) tim,
     &    Q_pneut,  
     &    (Q_pneut_trg(it), it=1,2), 
     &    Q_Te_min_LCFS,
     &    (Q_Te_max(it), it=1,ntrgts),
     &    (Q_det_wdth(it), it=1,ntrgts),
     &    (Q_qmax_trg(it), it=1,ntrgts),         !iyv 11.04.22
     &    (Q_jsat_max_trg(it), it=1,ntrgts),
     &    (Q_Isat_trg(it), it=1,ntrgts),
     &    (Q_I_tot_trg(it), it=1,ntrgts)
c?        end if

      return
      end subroutine b2usstrc
      
c======================================================================
c======================================================================
c======================================================================
c======================================================================

      subroutine b2usstri(nx,ny,ns,nxtl,nxtr)
      
      implicit none
      integer nx, ny, ns, nxtl, nxtr
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
      endif
      
c*** Input the DG data and prepare the geometry (this routine works once)
      call b2pwinp
      
      N_Q_total = N_Q_global + ntrgts * N_Q_target     

      if (.not. allocated(Q_pneut)) then
        allocate(
     &   Q_pneut, Q_pneut_trg(2), Area_pneut_trg(ntrgts),
     &   Q_Te_max(ntrgts), Q_det_wdth(ntrgts),
     &   Q_I_tot_trg(ntrgts), Q_jsat_trg(-1:ny,ntrgts),
     &   Q_jsat_max_trg(ntrgts), Q_Isat_trg(ntrgts), 
     &   Q_qmax_trg(ntrgts), Q_Te_min_LCFS)
        allocate(
     &   ix_pneut_beg(2), ix_pneut_end(2))
        allocate(headers(N_Q_total))
      endif
      
      mintopcut = ny   
      if(nncut > 0) mintopcut = minval(topcut(1:nncut))

      headers(1) = ' p_neut_PFR'
      headers(2) = 'p_neut_PFRi'
      headers(3) = 'p_neut_PFRo'
      headers(4) = 'Te_min_LCFS'
          
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
         ix_pneut_beg(1) =1
         ix_pneut_end(1) =0
         ix_pneut_beg(2) =1
         ix_pneut_end(2) =0
      end if !}
    
      if (ntrgts .eq. 2) then
        do i=1,N_Q_target
       headers(N_Q_global+(i-1)*ntrgts+1) = header_tmpl(i) // '_i '  ! inner target
       headers(N_Q_global+(i-1)*ntrgts+2) = header_tmpl(i) // '_o '  ! outer target
        enddo        
      elseif (ntrgts .eq. 4) then
        do i=1,N_Q_target
       headers(N_Q_global+(i-1)*ntrgts+1) = header_tmpl(i) // '_i '  ! inner target
       headers(N_Q_global+(i-1)*ntrgts+2) = header_tmpl(i) // '_iu'  ! outer target
       headers(N_Q_global+(i-1)*ntrgts+3) = header_tmpl(i) // '_ou'  ! inner upper target
       headers(N_Q_global+(i-1)*ntrgts+4) = header_tmpl(i) // '_o '  ! outer upper target
        enddo
      endif

      call ipgetr ('b2mndr_stim', stim)
      
!      write(*,*) 'filedata, ntrc = ', filedata(ntrc_uss), ntrc_uss 

      if(.not.filedata(ntrc_uss)) return

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      nc = 0
      lun = ntrc(ltrc)
      f = fn_ntrc(ltrc)
      
!      write(*,*) 'lun = ', lun
!      write(*,*) 'f = ', f
      
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
        
c?        if(i.lt.960727) then
c?          nc=16
c?        else
c?          go to 92
c?        end if
        
   42   continue
   
! Check version   
        if(nv.ne.lv) then
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2mod_usrspb: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
! FIXME
! These lines also should be revised, since they might lead to run-time error
! or, in some other cases, prevent it (depending on what is in input)
            filedata(ltrc)=.false.
            return
! FIXME
        end if
        
! Check total number of columns   
        if(N_ttl_read .ne. N_Q_total) then
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2mod_usrspb: current version ',lv,
     ,     ' differs from tha current number of columns ',N_Q_total,
     ,     ' differs from that in the file ',trim(f),' : ',N_ttl_read
c?          if(nc.gt.np_tmp) then
c?            write(*,*) 'No ',trim(f),' tracing file - ',
c?     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            return
c?          end if
c?          write(*,*) ' You should remove the file!'
        end if

! Check total number of targets   
        if(N_per_trg_read .ne. N_Q_target) then
          write(*,'(a,i2,3a,i2)')
     ,   '!!! b2mod_usrspb: current number of quantities per targets ',
     ,      N_Q_target,' differs from that in the file ',trim(f),' : ',
     ,      N_per_trg_read
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                              'remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if

! Check the number of quantities per target   
        if(Ntrgts_read .ne. ntrgts) then
          write(*,'(a,i2,3a,i2)')
     ,   '!!! b2mod_usrspb: the number of targets ',
     ,      ntrgts,' differs from that in the file ',trim(f),' : ',
     ,      Ntrgts_read
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                              'remove the existing file first!'
          filedata(ltrc)=.false.
          return
        end if
        
!!! FIXME
! This block is still leads to run-time error...
   50   read(lun,fmt='(a1)',end=55,err=55) hlp_c
        go to 50
   55   backspace lun
!!! FIXME

c----------------------------------------------------------------------
      else ! file 'user_SPb.trc' does not exist
      
        write(*,*) 'File user_SPb.trc does not exist!'

c?        nc = np_tmp

        write(lun, 700) lv, N_Q_total, ntrgts, N_Q_target, trim(label)
! What is 'label' and how is it defined?        
        
        write(hlp_frm,'(a,i3,a)')
     &   '(a1,6x,a4,3x,',N_Q_total,'(2x,a19))'
c       print *,'b2blns: nspecies= ',nspecies                       !###
c       print *,'species_list :',(species_list(j),':',j=1,nspecies) !###

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

      if (allocated(Q_pneut))
     &  deallocate(
     &   Q_pneut, Q_Te_min_LCFS, Q_pneut_trg, Area_pneut_trg,
     &   Q_Te_max, Q_det_wdth, Q_I_tot_trg,
     &   Q_jsat_trg, Q_Isat_trg, Q_qmax_trg)
        deallocate( 
     &   ix_pneut_beg, ix_pneut_end)
        deallocate(headers)     

      return
      end subroutine b2usstrf
      
c======================================================================
c======================================================================
c======================================================================
c======================================================================


      end module b2mod_usrspb

