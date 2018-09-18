*-----------------------------------------------------------------------
*.specification

      subroutine b2srst (nx, ny, ns, guard_flows, na, ua, te, ti, po,   !srv 01.07.08
     &  ne, ni, sna, smo, she, shi, sch, sne)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     1 , only : b2srst_sna0to1, b2srst_smo0to3, b2srst_she0to3,
     2          b2srst_shi0to3, balance_netcdf
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows                                   !srv 01.07.08
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
*   ..input/output arguments
      real (kind=R8) ::
     *  sna(-1:nx,-1:ny,0:1,0:ns-1), smo(-1:nx,-1:ny,0:3,0:ns-1),
     &  she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3),
     &  sch(-1:nx,-1:ny,0:3), sne(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRST modifies the source coefficients to provide numerical
*     stabilisation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, k, use_eirene
      integer, save :: iout = 0, !iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
!     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,            !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2wdat = 0
!     &                 iout_b2npmo =0, iout_b2npp7 = 0, iout_b2nppo = 0
      real (kind=R8) ::
     *  snast(-1:nx,-1:ny,0:1,0:ns-1), smost(-1:nx,-1:ny,0:3,0:ns-1), 
     &  shest(-1:nx,-1:ny,0:3), shist(-1:nx,-1:ny,0:3),
     &  schst(-1:nx,-1:ny,0:3), snest(-1:nx,-1:ny,0:1),
     &  wrk0(-1:nx,-1:ny)
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *   rf0, rf1, rf2, rf3, t0
      save ncall, rf0, rf1, rf2, rf3, use_eirene
*   ..procedures
      intrinsic min, max
      external subini, subend, xertst, ipgetr, b2xvsg
*   ..initialisation
      data ncall/0/, rf0/1.0_R8/, rf1/1.0_R8/, rf2/1.0_R8/, rf3/1.0_R8/,
     * use_eirene/0/
*     (rf0/1/2/3 for numerical stabilisation)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srst')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2srst_rf0', rf0)
       call ipgetr ('b2srst_rf1', rf1)
       call ipgetr ('b2srst_rf2', rf2)
       call ipgetr ('b2srst_rf3', rf3)
       call xertst (0.0_R8.le.rf0.and.0.0_R8.le.rf1.and.
     &              0.0_R8.le.rf2.and.0.0_R8.le.rf3,
     &              'faulty internal parameter rf0/1/2/3')
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('b2srst_iout', iout)                            !srv 30.06.08
!       call ipgeti ('b2npc9_iout', iout_b2npc9)                    !srv 11.09.09 {
!       call ipgeti ('b2npco_iout', iout_b2npco)
!       call ipgeti ('b2npc7_iout', iout_b2npc7)                    !srv 02.04.18
!       call ipgeti ('b2npc11_iout', iout_b2npc11)                  !srv 02.04.18
!       call ipgeti ('b2npmo_iout', iout_b2npmo)
!       call ipgeti ('b2npp7_iout', iout_b2npp7)
!       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                     !srv 11.09.09 }
       call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
      endif
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of na, ni, ne, te, ti
      if(ncall.lt.3) then
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (n2, ne, 1, 'ne', '.gt.')
        call b2xvsg (n2, te, 1, 'te', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.') 
      endif

*  .. save sources
      snast = sna
      smost = smo
      shest = she
      shist = shi
      schst = sch
      snest = sne

*  ..stabilise
*   ..modify sna
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
#ifndef NO_OPENMP_B2SRST
!$OMP PARALLEL
!$OMP& default(none)
!$OMP& private(is, iy, ix, t0)
!$OMP& shared(nx,ns,ny,guard_flows,is_neutral,use_eirene)
!$OMP& shared(leftix,rightix,bottomiy,topiy)
!$OMP& shared(na,rf0,rf1,rf2,rf3,ua,te,ti,po,ne,sne,snest)
! these are the result arrays:
!$OMP& shared(sna,smo,she,shi,sch)
!$OMP& shared(snast,smost,shest,shist,schst)
!$OMP DO
#endif
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
           if(guard_flows.eq.2 .and.
     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle 
         t0 = sna(ix,iy,0,is)+sna(ix,iy,1,is)*na(ix,iy,is)
         sna(ix,iy,0,is) =
     &     max(sna(ix,iy,0,is),(1.0_R8+rf0)*t0,-rf0*t0)
         sna(ix,iy,1,is) = (1.0_R8/na(ix,iy,is))*
     &     min(sna(ix,iy,1,is)*na(ix,iy,is),-rf0*t0,(1.0_R8+rf0)*t0)
        enddo
       enddo
       if(guard_flows.eq.2 .or.                                      !srv 06.07.08 {
     &    (is_neutral(is) .and. use_eirene.ne.0)) then               !srv 01.07.08 {
*     The corner cell value is defined by interpolation in b2npco and b2npc9
*   ..There is some problem for 9-point stencil case.
        do iy=-1, ny, ny+1                                           !srv 27.06.08 {
         do ix = -1, nx
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 ) then
             if(sna(ix,iy,1,is).eq.0.0_R8) sna(ix,iy,1,is) = -1.0_R8
           endif
         enddo
        enddo                                                        !srv 27.06.08 }
       endif                                                         !srv 01.07.08 }
      enddo
#ifndef NO_OPENMP_B2SRST
!$OMP END DO NOWAIT
#endif
*   ..modify smo
#ifndef NO_OPENMP_B2SRST
!$OMP DO
#endif
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
*      ..make smo(,,1,) negative
         smo(ix,iy,0,is) = smo(ix,iy,0,is)+
     &     max(0.0_R8,smo(ix,iy,1,is))*ua(ix,iy,is)
         smo(ix,iy,1,is) = min(smo(ix,iy,1,is),0.0_R8)
*      ..make smo(,,3,) negative
         smo(ix,iy,2,is) = smo(ix,iy,2,is)+
     &     max(0.0_R8,smo(ix,iy,3,is))*ua(ix,iy,is)
         smo(ix,iy,3,is) = min(smo(ix,iy,3,is),0.0_R8)
        enddo
       enddo
      enddo
#ifndef NO_OPENMP_B2SRST
!$OMP END DO NOWAIT
#endif
*   ..modify she, shi, sch
#ifndef NO_OPENMP_B2SRST
!$OMP DO
#endif
      do iy = -1, ny
       do ix = -1, nx
*     ..modify she(,,0:1)
        t0 = she(ix,iy,0)+she(ix,iy,1)*te(ix,iy)
        she(ix,iy,0) = max(she(ix,iy,0),(1.0_R8+rf1)*t0,-rf1*t0)
        she(ix,iy,1) = (1.0_R8/te(ix,iy))*
     &    min(she(ix,iy,1)*te(ix,iy),-rf1*t0,(1.0_R8+rf1)*t0)
*     ..modify she(,,2:3)
        t0 = she(ix,iy,2)+she(ix,iy,3)*te(ix,iy)
        she(ix,iy,2) = max(she(ix,iy,2),(1.0_R8+rf1)*t0,-rf1*t0)
        she(ix,iy,3) = (1.0_R8/te(ix,iy))*
     &    min(she(ix,iy,3)*te(ix,iy),-rf1*t0,(1.0_R8+rf1)*t0)
*     ..modify shi(,,0:1)
        t0 = shi(ix,iy,0)+shi(ix,iy,1)*ti(ix,iy)
        shi(ix,iy,0) = max(shi(ix,iy,0),(1.0_R8+rf2)*t0,-rf2*t0)
        shi(ix,iy,1) = (1.0_R8/ti(ix,iy))*
     &    min(shi(ix,iy,1)*ti(ix,iy),-rf2*t0,(1.0_R8+rf2)*t0)
*     ..modify shi(,,2:3)
        t0 = shi(ix,iy,2)+shi(ix,iy,3)*ti(ix,iy)
        shi(ix,iy,2) = max(shi(ix,iy,2),(1.0_R8+rf2)*t0,-rf2*t0)
        shi(ix,iy,3) = (1.0_R8/ti(ix,iy))*
     &    min(shi(ix,iy,3)*ti(ix,iy),-rf2*t0,(1.0_R8+rf2)*t0)
*     ..make sch(,,1) negative
        sch(ix,iy,0) = sch(ix,iy,0)+max(0.0_R8,sch(ix,iy,1))*po(ix,iy)
        sch(ix,iy,1) = min(sch(ix,iy,1),0.0_R8)
*     ..make sch(,,3) negative
        sch(ix,iy,2) = sch(ix,iy,2)+max(0.0_R8,sch(ix,iy,3))*po(ix,iy)
        sch(ix,iy,3) = min(sch(ix,iy,3),0.0_R8)
*     ..modify sne
        t0 = sne(ix,iy,0)+sne(ix,iy,1)*ne(ix,iy)
        sne(ix,iy,0) = max(sne(ix,iy,0),(1.0_R8+rf3)*t0,-rf3*t0)
        sne(ix,iy,1) = (1.0_R8/ne(ix,iy))*
     &   min(sne(ix,iy,1)*ne(ix,iy),-rf3*t0,(1.0_R8+rf3)*t0)
       enddo
      enddo
#ifndef NO_OPENMP_B2SRST
!$OMP END DO
#endif

#ifndef NO_OPENMP_B2SRST
!$OMP WORKSHARE
#endif
*   ..calculate the source coefficients of numerical stabilisation.
      snast = sna - snast
      smost = smo - smost
      shest = she - shest
      shist = shi - shist
      schst = sch - schst
      snest = sne - snest
#ifndef NO_OPENMP_B2SRST
!$OMP END WORKSHARE
!$OMP END PARALLEL
#endif
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif

cdjm Jan2017 Store for balance
      if (balance_netcdf.ne.0) then
        b2srst_sna0to1 = snast
        b2srst_smo0to3 = smost
        b2srst_she0to3 = shest
        b2srst_shi0to3 = shist
      endif

!      if (iout_b2npmo.eq.1) then                           !srv 11.09.09 {
!       do is = 0, ns-1
!         write (chns,'(i3.3)') is
!         wrk0=smost(:,:,0,is)+              !srv 07.10.09
!     &        smost(:,:,1,is)*ua(:,:,is)+
!     &        smost(:,:,2,is)*na(:,:,is)*mp*am(is)+
!     &        smost(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
!         call my_out(70,nx,ny,wrk0,'b2srst_smost'//chns)
!       enddo
!      endif
!      if (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1 .or.
!          iout_b2npc7.eq.1 .or. iout_b2npc11.eq.1) then
!       do is = 0, ns-1
!         write (chns,'(i3.3)') is
!         wrk0=snast(:,:,0,is)+              !srv 07.10.09
!     &        snast(:,:,1,is)*na(:,:,is)
!         call my_out(70,nx,ny,wrk0,'b2srst_snast'//chns)
!       enddo
!      endif
      if (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1 .or. iout_b2wdat.eq.4)
     &          then
        wrk0=shest(:,:,0)+shest(:,:,1)*te+  !srv 07.10.09
     &       shest(:,:,2)*ne + shest(:,:,3)*te*ne
        call my_out(70,nx,ny,wrk0,'b2srst_shest')
        wrk0=shist(:,:,0)+shist(:,:,1)*ti+
     &       shist(:,:,2)*ni(:,:,0)+shist(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrk0,'b2srst_shist')
      endif
!      if (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1) then
!        wrk0=schst(:,:,0)+schst(:,:,1)*po+  !srv 07.10.09
!     &       schst(:,:,2)*ne+schst(:,:,3)*po*ne(:,:)
!        call my_out(70,nx,ny,wrk0,'b2srst_schst')
!      endif

      if (iout.eq.1) then                       !srv 11.09.09 {
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          do k=0,3
            write (chk,'(i1)') k
            call my_out(70,nx,ny,smost(-1,-1,k,is),
     &                       'b2srst_smost'//chk//chns)
          enddo
          do k=0,1
            write (chk,'(i1)') k
            call my_out(70,nx,ny,snast(-1,-1,k,is),
     &                       'b2srst_snast'//chk//chns)
          enddo
        enddo
        do k=0,3
          write (chk,'(i1)') k
          call my_out(70,nx,ny,shest(-1,-1,k),'b2srst_shest'//chk)
          call my_out(70,nx,ny,shist(-1,-1,k),'b2srst_shist'//chk)
          call my_out(70,nx,ny,schst(-1,-1,k),'b2srst_schst'//chk)
        enddo
      endif                                                 !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srst

      end subroutine b2srst
