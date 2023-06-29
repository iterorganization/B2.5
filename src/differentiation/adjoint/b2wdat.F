      subroutine b2wdat (nCv, nFc, nVx, ns, ismain, ismain0,
     .                   mpg, geo, state, state_ext, BoRiS)     !srv 05.10.01
      use b2mod_types
      use b2mod_constants
      use b2mod_neutrals_namelist_diff
      use b2mod_boundary_namelist_diff
      use b2mod_elements_diff
      use b2mod_b2cmgs
      use b2mod_b2cmpa_diff
      use b2mod_b2cmpb_diff
      use b2mod_b2cmpt_diff
      use b2mod_b2cmrc_diff
#ifdef B25_EIRENE
      use eirmod_parmmod
     *   , only : natm, nmol, nstra,  npls
      use eirmod_extrab25
#endif
      use b2us_geo_diff
      use b2us_map_diff
      use b2us_plasma_diff
      use b2mod_user_namelist_diff
     , , only : nomp, omp, icsepomp, nimp, imp, icsepimp
#ifdef WG_DEBUG
      use b2us_debug
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      real (kind=R8) :: BoRiS
      type(mapping)    :: mpg
      type(geometry)   :: geo
      type(B2State)    :: state
      type(B2StateExt) :: state_ext
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WDAT is an output routine that writes out B2 arrays into ASCII
*     data files.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer iCv, iCv1, iCv2, iFc, ibc, iReg
      integer ncall, iout, is, ic, ib, iout_iout_ASTRA             !srv 05.10.01 03.04.13 19.04.13
      integer nx, ny                                               !srv 12.07.10
      integer use_astra
      integer boundary_namelist
      integer, save :: ishift=0                                    !srv 05.11.09
      integer, save :: use_mms = 0                                 !srv 26.02.18

      integer iRFc, nRFc, icell
#ifdef WG_DEBUG
      integer nmdpl, nmdpl_i, n_closed_surface, ix, iy,
     ,        ilast_closed_surface, ifirst_open_surface
      logical first
      logical found
#endif
      integer, allocatable :: regFaces(:)
      real (kind=R8), save :: shift = 0.0_R8                       !srv
      real (kind=R8), allocatable :: fun_st(:,:), clength_r(:),
     *  clength_s(:)
!   clength_r(i) is the length from the cell center of 1th cell 
!                to the cell center of i-th cell in "radial direction"
!   clength_s(i) is the length along a boundary from the face center of 1st face
!                to the face center of i-th face

#ifdef WG_TODO
      real (kind=R8) ::                                            !srv
     *  Total_current(-1:ny), B0out, Total_flux(-1:ny),            !srv 31.08.00 19.02.02
     *  Total_abscurrent(-1:ny), Total_absflux(-1:ny),
     *  Total_energy_flux_i(-1:ny), Total_energy_flux_e(-1:ny),    !srv 23.11.00
#endif
      real (kind=R8) ::                                            !srv
     *  wrk2(nCv)                                                  !srv 23.11.00
#ifdef WG_TODO
      real (kind=R8) ::                                            !srv
     *  Total_energy_flux(-1:ny), power,                           !srv 23.11.00
     *  Total_energy_flux_mdf(-1:ny),
     *  Total_energy_flux_i_mdf(-1:ny),
     *  Total_energy_flux_e_mdf(-1:ny),
!     *  kintmp,
     *  sfuadia(-1:ny,0:ns-1), sfvadia(-1:ny,0:ns-1),
     *  sfvaecrb(-1:ny,0:ns-1), sfna(-1:ny,0:ns-1),
     &  sfna_mdf(-1:ny,0:ns-1), sfn(-1:ny), sfn_mdf(-1:ny)
      real (kind=R8) ::                                            !srv
     &  pz(-1:nx,-1:ny), rz(-1:nx,-1:ny),
     &  pr(-1:nx,-1:ny), ro(-1:nx,-1:ny), cst(-1:nx,-1:ny)         !srv 30.01.15
      real (kind=R8) :: bsp(-1:ny), csp(-1:ny), dsp(-1:ny), dy,
     &     cr(4*(ny+1)+1), dseval, wdim(8*ns)                      !srv 03.04.13
#endif
      character chis*3, chib*2                                     !srv 27.04.01
      character*32 hlp_frm                                         !srv 03.04.13
      character(7), save :: folder                                 !srv 26.09.17
#ifdef B25_EIRENE
      real (kind=R8), save :: hlp_neu_peak
      character chns*3
      data hlp_neu_peak / 5.0_R8 /                                 !iyv 15.10.09
#endif
*   ..procedures
      external ipgeti
*   ..initialisation
      save ncall, iout, use_astra, boundary_namelist
      data ncall/0/, iout/0/, use_astra/0/, boundary_namelist/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wdat')
*   ..initialise files and common

*   ..set internal parameters on first (and only) call
      if (ncall.eq.0) then
       call ipgeti('b2wdat_iout', iout)                           !srv 05.10.01
       call ipgeti('b2wdat_ishift', ishift)                       !srv 05.10.01
       call ipgeti('b2mndr_use_mms',use_mms)                      !srv 26.02.18
       call ipgeti('b2mndr_astra', use_astra)
       call ipgeti('b2stbc_boundary_namelist',boundary_namelist)
*    ..try to create a file in optional output folder
       open(70,file='output/__1__st.dat',err=10)
*    ..folder exists: future output will be done in folder
       folder='output/'
       write(70,*) 'Test file to check on output folder'
       close(70)
       go to 20
*    ..folder does not exist: future output will be done in current folder
   10  folder=''
   20  continue
      endif

      ! At the moment, some output is only supported for grids that were
      ! originally structured
      if (mpg%isClassicalGrid.eq.1) then

*    ..set nx, ny
       nx = mpg%nx
       ny = mpg%ny

*    ..allocate and initialize some arrays for structured format
       allocate(fun_st(-1:nx,-1:ny))
       fun_st = 0.0_R8

      endif

*   ..initialisation
      if (use_astra.ne.0) then
        iout_iout_ASTRA = 1
      else
        iout_iout_ASTRA = 0
      endif

#ifdef WG_DEBUG
      if (mpg%isClassicalGrid.eq.1) then

*      ..find the number of the last closed surface in CORE
         ilast_closed_surface = -2
         n_closed_surface = 0
         ix = nx/2    !  nmdpl; may not work for DN?
         do iy = -1, ny
          if(mpg%imapCv(ix,iy) /= 0) then
           if(mpg%cvOnClosedSurface(mpg%imapCv(ix,iy))) then
            n_closed_surface = n_closed_surface + 1
            ilast_closed_surface = iy
           endif
          endif
         enddo
         if (ilast_closed_surface.eq.-2) 
     &    call xerrab('Problem identifying last closed surface')
*      ..find the number of the first open surface in SOL
         ifirst_open_surface = ilast_closed_surface
         found = .false.
         do iy = ilast_closed_surface + 1, ny
          if(mpg%imapCv(ix,iy) /= 0 .and. .not.found) then
           if(.not.mpg%cvOnClosedSurface(mpg%imapCv(ix,iy))) then
            ifirst_open_surface = iy
            found = .true.
           endif
          endif
         enddo
         if (ifirst_open_surface.eq.ilast_closed_surface) 
     &    call xerrab('Problem identifying first open surface')
      
*      ..find the poloidal number of midplane cell as a number of cell with min abs(Bz)
*        and inner midplane as cell with max abs(Bz)
         iy=ilast_closed_surface
         first=.true.
         do ix = -1, nx
          if(mpg%imapCv(ix,iy) /= 0) then
           if(mpg%cvOnClosedSurface(mpg%imapCv(ix,iy))) then
            if(first) then
                nmdpl = ix
                nmdpl_i = ix
                first = .false.
            endif         
            if(abs(geo%cvBb(mpg%imapCv(ix,iy),2)) .lt. 
     &         abs(geo%cvBb(mpg%imapCv(nmdpl,iy),2))) nmdpl = ix
            if(abs(geo%cvBb(mpg%imapCv(ix,iy),2)) .gt. 
     &         abs(geo%cvBb(mpg%imapCv(nmdpl_i,iy),2))) nmdpl_i = ix
           endif
          endif
         enddo
      endif
#endif

      if (iout.eq.1) then                                          !srv begin if for iout
*    .. iout = 1 for outputting plasma state files
#ifdef B25_EIRENE
       if (iout.ne.4) then
        open(70,file=trim(folder)//'eirdiag.dat')                   !srv 26.09.17
        write(70,"(a/6i5)") 'nlimps natm nmol npls nstra nstsid',
     &                       nlimps, natm, nmol, npls, nstra, nstsid
        write(70,"(a/6i5)") 'natmid nmolid ', natmid, nmolid
        write(70,"(a/6i5)") 'nnatmi nnmoli ', nnatmi, nnmoli
        close(70)
        do is = 1,nnatmi                                            !srv 04.06.09 16.02.12{
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,0,dab2(1,is,1),'dab2'//chns)
         call my_out_us(70,nCv,0,tab2(1,is,1),'tab2'//chns)
         call my_out_us(70,nCv,0,rfluxa(1,is,1),'rfluxa'//chns)
         call my_out_us(70,nCv,0,refluxa(1,is,1),
     .                       'refluxa'//chns)
         call my_out_us(70,nCv,0,pfluxa(1,is,1),'pfluxa'//chns)
         call my_out_us(70,nCv,0,pefluxa(1,is,1),
     .                       'pefluxa'//chns)
        enddo
        do is = 1,nnmoli                                            !srv 16.02.12
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,0,dmb2(1,is,1),'dmb2'//chns) !srv 29.04.10
         call my_out_us(70,nCv,0,tmb2(1,is,1),'tmb2'//chns) !srv 29.04.10
         call my_out_us(70,nCv,0,rfluxm(1,is,1),'rfluxm'//chns)
         call my_out_us(70,nCv,0,refluxm(1,is,1),'refluxm'//chns)
         call my_out_us(70,nCv,0,pfluxm(1,is,1),'pfluxm'//chns)
         call my_out_us(70,nCv,0,pefluxm(1,is,1),'pefluxm'//chns)
        enddo                                                       !srv 04.06.09 }
       endif
#endif

       wrk2=0.0_R8
       do is = 0,ns-1
         write (chis,'(i3.3)') is
         call my_out_us(70,nCv,0,state%pl%na(1,is),'na'//chis)           !srv
         call my_out_us(70,nCv,0,state%pl%ua(1,is),'ua'//chis)           !srv
       end do

       call my_out_us(70,nCv,0,state%dv%ne,'ne')                               !srv
       wrk2=state%pl%te/ev
       call my_out_us(70,nCv,0,wrk2,'te_eV')
       wrk2=state%pl%ti/ev
       call my_out_us(70,nCv,0,wrk2,'ti_eV')
       call my_out_us(70,nCv,0,state%pl%po,'po')                                !srv

        wrk2=0.0_R8
        do ib = 1, nbc
            write (chib,'(i2.2)') ib
            is=1
            write (chis,'(i3.3)') is
            do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
              iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
              iCv2 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,2)                   ! number of the corresponding domain cell
              iFc  = mpg%cvFc(mpg%cvFcP(iCv1,1))                         ! number of the guard cell face
              wrk2(iCv1)=state%pl%na(iCv1,is)
            enddo          
            call my_out_us(70,nCv,0,wrk2,'na_'//chis//'_bnd_'//chib)
        enddo
       
        do ib = 1, nbc
          write (chib,'(i2.2)') ib
          wrk2=0.0_R8
            do ibc = 1, mpg%bcCvP(ib,2)                                  ! loop over number of cells in the boundary
              iCv1 = mpg%bcCv(mpg%bcCvP(ib,1)+ibc-1,1)                   ! number of the guard cell
              wrk2(iCv1)=state%dv%ne(iCv1)
            enddo          
          call my_out_us(70,nCv,0,wrk2,'ne_bnd_'//chib)
        enddo
      
      endif              !srv   endif for iout=1   
          
      if (iout.eq.2) then                                               !srv begin for iout=2
*    .. iout = 2 for outputting profiles of ne, te and ti at the plates and midplane
        if (mpg%isClassicalGrid.eq.1) then
          !! Identification of target boundaries currently based on (face) region
          !! numbers as used in structured code.
          !! To be generalized to unstructured grids; e.g. loop over all face regions?

          !! INNER TARGET
          open(70,file=trim(folder)//'inner.dat')

          ! identify face region number of inner target
          iReg = 0
          if (maxval(mpg%cvReg).gt.1) iReg = 1
          if (iReg.eq.0) call xerrab('Problem identifying inner target')

          ! count number of faces in the region
          nRFc = count(mpg%fcReg==iReg)
          if (nRFc.eq.0) call xerrab('No faces on inner target')
          allocate (regFaces(nRFc))
          allocate (clength_s(nRFc))

          ! sort faces into continuous surface
          call sort_faces_nodiff(nRFc,regFaces,nFc,mpg%fcReg,iReg,mpg)

!          write (*,*) 'Inner target faces, nRFc = ', nRFc, 
!     &     ', indices = ', regFaces(1:nRFC)

          ! compute distance along inner target
          call calc_dist_f_nodiff(mpg,geo,regFaces,nRFc,0,clength_s)

         ! write properties along inner target
          do iRFc = 1, nRFc
           iFc  = regFaces(iRFc)
           iCv1 = max(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))        ! number of guard cell
!         write (*,*) 'iFc, iCv1 = ', iFc, iCv1
           if (iCv1.le.mpg%nCi) 
     &       call xerrab('Problem identifying guard cell')
           write(70,"(6es14.6)")  clength_s(iRFc)- shift,
     &           state%dv%ne(iCv1), state%pl%te(iCv1)/ev,
     &           state%pl%ti(iCv1)/ev
          enddo          
          deallocate (clength_s)
          deallocate (regFaces)
          close(70)

          !! OUTER TARGET
          open(70,file=trim(folder)//'outer.dat')

          ! identify face region number of outer target
          iReg = 0
          if (maxval(mpg%cvReg).eq.2) iReg = 2 ! limiter case
          if (maxval(mpg%cvReg).gt.2) iReg = 4 ! SN case
          if (maxval(mpg%cvReg).gt.5) iReg = 8 ! DN case
          if (iReg.eq.0) call xerrab('Problem identifying outer target')

          ! count number of faces in the region
          nRFc = count(mpg%fcReg==iReg)
          if (nRFc.eq.0) call xerrab('No faces on outer target')
          allocate (regFaces(nRFc))
          allocate (clength_s(nRFc))

          ! sort faces into continuous surface
          call sort_faces_nodiff(nRFc,regFaces,nFc,mpg%fcReg,iReg,mpg)

!          write (*,*) 'Outer target faces, nRFc = ', nRFc, 
!     &     ', indices = ', regFaces(1:nRFC)

          ! compute distance along outer target
          call calc_dist_f_nodiff(mpg,geo,regFaces,nRFc,0,clength_s)

          ! write properties along outer target
          do iRFc = 1, nRFc
           iFc  = regFaces(iRFc)
           iCv1 = max(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))        ! number of guard cell
!         write (*,*) 'iFc, iCv1 = ', iFc, iCv1
           if (iCv1.le.mpg%nCi) 
     &       call xerrab('Problem identifying guard cell')
           write(70,"(6es14.6)")  clength_s(iRFc)- shift,
     &           state%dv%ne(iCv1), state%pl%te(iCv1)/ev,
     &           state%pl%ti(iCv1)/ev
          enddo          
          deallocate (clength_s)
          deallocate (regFaces)
          close(70)

          !! OUTER MIDPLANE
          open(70,file=trim(folder)//'midpln.dat')
          allocate (clength_r(1:nomp))
          call calc_dist_nodiff(mpg,geo,omp,nomp,icsepomp,clength_r)

          icell = 1
          do icell = 1, nomp
            iCv = omp(icell)
            write(70,"(4es14.6)") clength_r(icell)-shift,
     &            state%dv%ne(iCv), state%pl%te(iCv)/ev, 
     &            state%pl%ti(iCv)/ev
          enddo
          deallocate (clength_r)
          close(70)

          !! INNER MIDPLANE
          open(70,file=trim(folder)//'midpln_i.dat')
          allocate (clength_r(1:nimp))
          call calc_dist_nodiff(mpg,geo,imp,nimp,icsepimp,clength_r)

          icell = 1
          do icell = 1, nimp
            iCv = imp(icell)
            write(70,"(4es14.6)") clength_r(icell)-shift,
     &            state%dv%ne(iCv), state%pl%te(iCv)/ev, 
     &            state%pl%ti(iCv)/ev
          enddo
          deallocate (clength_r)
          close(70)

        endif   
      endif                                                             !srv   end for iout=2   
      
      if(iout.eq.3 .or. iout.eq.4 .or. iout.eq.1) then                  !srv begin for iout=3
*    .. iout = 3 for outputting geometry files

! output of cell-centered values
       call my_out_us(70,nCv,0,geo%cvBb(1,0),'Cv.B_pol')
       call my_out_us(70,nCv,0,geo%cvBb(1,2),'Cv.B_tor')
       call my_out_us(70,nCv,0,geo%cvBb(1,3),'Cv.B_full')
       call my_out_us(70,nCv,0,geo%cvEb(1,0),'Cv.b_R_cyl')
       call my_out_us(70,nCv,0,geo%cvEb(1,1),'Cv.b_Z_cyl')
       call my_out_us(70,nCv,0,geo%cvEb(1,2),'Cv.b_phi_cyl')
       call my_out_us(70,nCv,0,geo%cvHx,'Cv.lenght_pol')
       call my_out_us(70,nCv,0,geo%cvHz,'Cv.hz')
       call my_out_us(70,nCv,0,geo%cvVol,'Cv.vol')
       call my_out_us(70,nCv,0,geo%cvSz,'Cv.Sz')
       call my_out_us(70,nCv,0,geo%cvOnedBsq,'Cv.OnedBsq')
       call my_out_us(70,nCv,0,geo%cvBzB,'Cv.BzB')
       call my_out_us(70,nCv,0,geo%cvX,'Cv.cRc')                       ! R-coordinate of the (ix,iy) cell centre
       call my_out_us(70,nCv,0,geo%cvY,'Cv.cZc')                       ! Z-coordinate of the (ix,iy) cell centre 
       call my_outi_us(70,nCv,0,mpg%cvReg,'Cv.reg')                    !lkw 20.09.2022
         open(70,file=trim(folder)//'czc.dat')                         !srv 26.09.17
         write(hlp_frm,"(a1,i4,a9)") '(',10,'es24.16))'                !srv 07.08.15
         write(70,hlp_frm) (geo%cvY(ic),ic=1,nCv)
         close(70)
         open(70,file=trim(folder)//'crc.dat')                         !srv 26.09.17
         write(hlp_frm,"(a1,i4,a9)") '(',10,'es24.16))'                !srv 07.08.15
         write(70,hlp_frm) (geo%cvX(ic),ic=1,nCv)
         close(70)
       
#ifdef WG_DEBUG
       if (mpg%isClassicalGrid.eq.1) then
         open(70,file=trim(folder)//'info.dat')
         write(70,"(' Number of outer midplane cell =',i4)") nmdpl    !srv 24.03.06 }
         write(70,"(' Number of inner midplane cell =',i4)") nmdpl_i  !srv 24.03.06 }
         write (70,"(' Number of the last closed surface in CORE ='"//
     &           ",i4)") ilast_closed_surface
         write (70,"(' Number of closed surfaces in CORE =',i4)")
     &           n_closed_surface
         write (70,"(' Number of the first open surface in SOL =',i4)")
     &           ifirst_open_surface
         close(70)
        do iy = -1, ny_st
         do ix = -1, nx_st
           if(imapFcx(ix,iy) /= 0) then                        ! poloidal faces
             fun_st(ix,iy) = 2.0_R8*geo%fcHc(imapFcx(ix,iy),2) ! get hx in cell centre
           endif
         enddo
        enddo
*     ..restore the geometry of guard cells (should be redone for DN case)
        ix=-1
        fun_st(ix,:)=fun_st(ix+1,:)/1.0e+2_R8 ! for west 
        iy=-1
        fun_st(:,iy)=fun_st(:,iy+1)           ! for south
        iy=ny
        fun_st(:,iy)=fun_st(:,iy-1)           ! for north
        call my_out(70, nx, ny, fun_st, 'hx_st')
       
        do iy = -1, ny_st
         do ix = -1, nx_st
           if(imapFcy(ix,iy) /= 0) then                        ! radial faces
             fun_st(ix,iy) = 2.0_R8*geo%fcHc(imapFcy(ix,iy),2) ! get hy  in cell centre 
           endif
         enddo
        enddo
*     ..restore the geometry of guard cells (should be redone for DN case)
        iy=-1
        fun_st(:,iy)=fun_st(:,iy+1)/1.0e+2_R8   ! for south
        ix=-1
        fun_st(ix,:)=fun_st(ix+1,:)             ! for west 
        ix=nx
        fun_st(ix,:)=fun_st(ix-1,:)             ! for east 
        call my_out(70, nx, ny, fun_st, 'hy_st')
        
        do iy = -1, ny_st
         do ix = -1, nx_st
           if(imapCv(ix,iy) /= 0) then
             fun_st(ix,iy) = fun_st(ix,iy)*geo%cvQgam(imapCv(ix,iy),0) ! get hy1=hy*qz(ix,iy,1)  in cell centre qz(ix,iy,1) is cos
           endif
         enddo
        enddo
        call my_out(70, nx, ny, fun_st, 'hy1_st')
       endif 
#endif
      
! output of values corresponding to cell faces
       call my_out_us(70,nFc,1,geo%fcHt,'Fc.ht')
       call my_out_us(70,nFc,1,geo%fcPbs,'Fc.pbs')
       call my_out_us(70,nFc,1,geo%fcBb(1,0),'Fc.B_pol')
       call my_out_us(70,nFc,1,geo%fcBb(1,2),'Fc.B_tor')
       call my_out_us(70,nFc,1,geo%fcBb(1,3),'Fc.B_full')
       call my_out_us(70,nFc,0,geo%fcEb(1,0),'Fc.b_R_cyl')
       call my_out_us(70,nFc,0,geo%fcEb(1,1),'Fc.b_Z_cyl')
       call my_out_us(70,nFc,0,geo%fcEb(1,2),'Fc.b_phi_cyl')
       call my_out_us(70,nFc,1,geo%fcHz,'Fc.hz')
       call my_out_us(70,nFc,1,geo%fcS,'Fc.S')
       call my_out_us(70,nFc,1,geo%fcPbsHz,'Fc.pbshz')
       call my_out_us(70,nFc,1,geo%fcQgam(1,0),'Fc.cos_gamma')
       call my_out_us(70,nFc,1,geo%fcQgam(1,1),'Fc.sin_gamma')
       call my_out_us(70,nFc,1,geo%fcQbet(1,0),'Fc.cos_beta')
       call my_out_us(70,nFc,1,geo%fcQbet(1,1),'Fc.sin_beta')
       call my_out_us(70,nFc,1,geo%fcQalf(1,0),'Fc.cos_alpha')
       call my_out_us(70,nFc,1,geo%fcQalf(1,1),'Fc.sin_alpha')

! output of values corresponding to vertices
       call my_out_us(70,nVx,2,geo%vxX,'Vx.cR')
       call my_out_us(70,nVx,2,geo%vxY,'Vx.cZ')
       call my_out_us(70,nVx,2,geo%vxHz,'Vx.hz')
       call my_out_us(70,nVx,2,geo%vxBb(1,0),'Vx.B_pol')
       call my_out_us(70,nVx,2,geo%vxBb(1,2),'Vx.B_tor')
       call my_out_us(70,nVx,2,geo%vxBb(1,3),'Vx.B_full')
       call my_out_us(70,nVx,2,geo%vxFPsi,'Vx.FPsi')
       call my_out_us(70,nVx,2,geo%vxFfbz,'Vx.FfBz')
       call my_out_us(70,nVx,2,geo%vxOnedBsq,'Vx.OnedBsq')
       call my_out_us(70,nVx,2,geo%vxBzB,'Vx.BzB')
      endif                                                             !srv   end for iout=3   

      if(use_mms.eq.1 .and. iout.eq.10) then                                 !srv 26.02.18 {
!WG_TODO        call calc_err(nx, ny, ns)
      endif                                                             !srv 26.02.18 }

      if (mpg%isClassicalGrid.eq.1) then

*    ..deallocate arrays for structured format
       if (allocated(fun_st)) deallocate(fun_st)

      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wdat

      end subroutine b2wdat
