*-----------------------------------------------------------------------
*.specification

      subroutine b2srst (nCv, ns, switch, na, ua, te, ti, tn, po,       !srv 01.07.08
     &  ne, ni, nn, kt, zt, sr)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_plasma
      use b2mod_balance !djm Jan2017
     1 , only : b2srst_sna0to1, b2srst_smo0to3, b2srst_she0to3,
     2          b2srst_shi0to3, balance_netcdf
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2srst, b2srst_kt_eps, b2srst_zt_eps
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns                                  !srv 01.07.08
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1),
     &  te(nCv), ti(nCv), tn(nCv), po(nCv),
     &  ne(nCv), ni(nCv,0:1), nn(nCv), kt(nCv), zt(nCv)
*   ..input/output arguments
      type (B2Source), intent (inout) :: sr
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
      integer iCv, is, k
      real (kind=R8) ::
     *  snast(nCv,0:1,0:ns-1), smost(nCv,0:3,0:ns-1),
     &  shest(nCv,0:3), shist(nCv,0:3),
     &  schst(nCv,0:3), snest(nCv,0:1),
     &  sktst(nCv,0:3), sztst(nCv,0:3), wrk0(nCv)
      character chns*3, chk*1                                         !srv 11.09.09 }
      real (kind=R8) :: t0
*   ..procedures
      intrinsic min, max
      external xertst, ipgetr, b2xvsg
*   ..initialisation
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srst')
*   ..set internal parameters on first call
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of na, ni, ne, te, ti
      if(ncall_b2srst.lt.3) then
        call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
        call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (nCv, nn, 1, 'nn', '.gt.')
        call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
        call b2xvsg (nCv, te, 1, 'te', '.gt.')
        call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
        call b2xvsg (nCv, kt, 1, 'kt', '.ge.')
        call b2xvsg (nCv, zt, 1, 'zt', '.ge.')
      endif

*  .. save sources
      snast = sr%sna
      smost = sr%smo
      shest = sr%she
      shist = sr%shi
      schst = sr%sch
      snest = sr%sne
      sktst = sr%skt
      sztst = sr%szt

*  ..stabilise
*   ..modify sna
      do is = 0, ns-1
       do iCv = 1, nCv
!WG_TODO           if(guard_flows.eq.2 .and.
!WG_TODO     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!WG_TODO     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle 
         t0 = sr%sna(iCv,0,is)+sr%sna(iCv,1,is)*na(iCv,is)
         sr%sna(iCv,0,is) =
     &     max(sr%sna(iCv,0,is),(1.0_R8+switch%b2srst_rf0)*t0,
     &     -switch%b2srst_rf0*t0)
         sr%sna(iCv,1,is) = (1.0_R8/na(iCv,is))*
     &     min(sr%sna(iCv,1,is)*na(iCv,is),-switch%b2srst_rf0*t0,
     &     (1.0_R8+switch%b2srst_rf0)*t0)
       enddo
!WG_TODO       if(guard_flows.eq.2 .or.                                      !srv 06.07.08 {
!WG_TODO     &    (is_neutral(is) .and. use_eirene.ne.0)) then               !srv 01.07.08 {
!WG_TODO*     The corner cell value is defined by interpolation in b2npco
!WG_TODO*   ..There is some problem for 9-point stencil case.
!WG_TODO        do iy=-1, ny, ny+1                                           !srv 27.06.08 {
!WG_TODO          do ix = -1, nx
!WG_TODO            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 ) then
!WG_TODO              if(sna(ix,iy,1,is).eq.0.0_R8) sna(ix,iy,1,is) = -1.0_R8
!WG_TODO            endif
!WG_TODO          enddo
!WG_TODO        enddo                                                        !srv 27.06.08 }
!WG_TODO       endif                                                         !srv 01.07.08 }
      enddo
*   ..modify smo
      do is = 0, ns-1
       do iCv = 1, nCv
*      ..make smo(,,1,) negative
         sr%smo(iCv,0,is) = sr%smo(iCv,0,is)+
     &     max(0.0_R8,sr%smo(iCv,1,is))*ua(iCv,is)
         sr%smo(iCv,1,is) = min(sr%smo(iCv,1,is),0.0_R8)
*      ..make smo(,,3,) negative
         sr%smo(iCv,2,is) = sr%smo(iCv,2,is)+
     &     max(0.0_R8,sr%smo(iCv,3,is))*ua(iCv,is)
         sr%smo(iCv,3,is) = min(sr%smo(iCv,3,is),0.0_R8)
       enddo
      enddo
*   ..modify she, shi, sch, skt, szt
      do iCv = 1, nCv
*     ..modify she(,0:1)
        t0 = sr%she(iCv,0)+sr%she(iCv,1)*te(iCv)
        sr%she(iCv,0) = max(sr%she(iCv,0),(1.0_R8+switch%b2srst_rf1)*t0,
     &                  -switch%b2srst_rf1*t0)
        sr%she(iCv,1) = (1.0_R8/te(iCv))*
     &    min(sr%she(iCv,1)*te(iCv),-switch%b2srst_rf1*t0,
     &    (1.0_R8+switch%b2srst_rf1)*t0)
*     ..modify she(,2:3)
        t0 = sr%she(iCv,2)+sr%she(iCv,3)*te(iCv)
        sr%she(iCv,2) = max(sr%she(iCv,2),(1.0_R8+switch%b2srst_rf1)*t0,
     &                  -switch%b2srst_rf1*t0)
        sr%she(iCv,3) = (1.0_R8/te(iCv))*
     &    min(sr%she(iCv,3)*te(iCv),-switch%b2srst_rf1*t0,
     &    (1.0_R8+switch%b2srst_rf1)*t0)
*     ..modify shi(,0:1)
        t0 = sr%shi(iCv,0)+sr%shi(iCv,1)*ti(iCv)
        sr%shi(iCv,0) = max(sr%shi(iCv,0),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%shi(iCv,1) = (1.0_R8/ti(iCv))*
     &    min(sr%shi(iCv,1)*ti(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify shi(,2:3)
        t0 = sr%shi(iCv,2)+sr%shi(iCv,3)*ti(iCv)
        sr%shi(iCv,2) = max(sr%shi(iCv,2),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%shi(iCv,3) = (1.0_R8/ti(iCv))*
     &    min(sr%shi(iCv,3)*ti(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify shn(,0:1)
        t0 = sr%shn(iCv,0)+sr%shn(iCv,1)*tn(iCv)
        sr%shn(iCv,0) = max(sr%shn(iCv,0),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%shn(iCv,1) = (1.0_R8/tn(iCv))*
     &    min(sr%shn(iCv,1)*tn(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify shn(,2:3)
        t0 = sr%shn(iCv,2)+sr%shn(iCv,3)*tn(iCv)
        sr%shn(iCv,2) = max(sr%shn(iCv,2),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%shn(iCv,3) = (1.0_R8/tn(iCv))*
     &    min(sr%shn(iCv,3)*tn(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify skt(,0:1)
        t0 = sr%skt(iCv,0)+sr%skt(iCv,1)*kt(iCv)
        sr%skt(iCv,0) = max(sr%skt(iCv,0),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%skt(iCv,1) = (1.0_R8/(kt(iCv)+b2srst_kt_eps))*
     &    min(sr%skt(iCv,1)*kt(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify skt(,2:3)
        t0 = sr%skt(iCv,2)+sr%skt(iCv,3)*kt(iCv)
        sr%skt(iCv,2) = max(sr%skt(iCv,2),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%skt(iCv,3) = (1.0_R8/(kt(iCv)+b2srst_kt_eps))*
     &    min(sr%skt(iCv,3)*kt(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify szt(,0:1)
        t0 = sr%szt(iCv,0)+sr%szt(iCv,1)*zt(iCv)
        sr%szt(iCv,0) = max(sr%szt(iCv,0),(1.0_R8+switch%b2srst_rf2)*t0,
     &                  -switch%b2srst_rf2*t0)
        sr%szt(iCv,1) = (1.0_R8/(zt(iCv)+b2srst_zt_eps))*
     &    min(sr%szt(iCv,1)*zt(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)
*     ..modify szt(,2:3)
        t0 = sr%szt(iCv,2)+sr%szt(iCv,3)*zt(iCv)
        sr%szt(iCv,2) = max(sr%szt(iCv,2),(1.0_R8+switch%b2srst_rf2)*t0,
     &                 -switch%b2srst_rf2*t0)
        sr%szt(iCv,3) = (1.0_R8/(zt(iCv)+b2srst_zt_eps))*
     &    min(sr%szt(iCv,3)*zt(iCv),-switch%b2srst_rf2*t0,
     &    (1.0_R8+switch%b2srst_rf2)*t0)


*     ..make sch(,1) negative
        sr%sch(iCv,0) = sr%sch(iCv,0)+max(0.0_R8,sr%sch(iCv,1))*po(iCv)
        sr%sch(iCv,1) = min(sr%sch(iCv,1),0.0_R8)
*     ..make sch(,3) negative
        sr%sch(iCv,2) = sr%sch(iCv,2)+max(0.0_R8,sr%sch(iCv,3))*po(iCv)
        sr%sch(iCv,3) = min(sr%sch(iCv,3),0.0_R8)
*     ..modify sne
        t0 = sr%sne(iCv,0)+sr%sne(iCv,1)*ne(iCv)
        sr%sne(iCv,0) = max(sr%sne(iCv,0),(1.0_R8+switch%b2srst_rf3)*t0,
     &                  -switch%b2srst_rf3*t0)
        sr%sne(iCv,1) = (1.0_R8/ne(iCv))*
     &   min(sr%sne(iCv,1)*ne(iCv),-switch%b2srst_rf3*t0,
     &   (1.0_R8+switch%b2srst_rf3)*t0)
      enddo

*   ..calculate the source coefficients of numerical stabilisation.
      snast = sr%sna - snast
      smost = sr%smo - smost
      shest = sr%she - shest
      shist = sr%shi - shist
      schst = sr%sch - schst
      snest = sr%sne - snest
      sktst = sr%skt - sktst
      sztst = sr%szt - sztst

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
!         wrk0=smost(:,0,is)+                               !srv 07.10.09
!     &        smost(:,1,is)*ua(:,is)+
!     &        smost(:,2,is)*na(:,is)*mp*am(is)+
!     &        smost(:,3,is)*ua(:,is)*na(:,is)*mp*am(is)
!         call my_out_us(70,nCv,0,wrk0,'b2srst_smost'//chns)
!       enddo
!      endif
!      if (iout_b2npco.eq.1) then
!       do is = 0, ns-1
!         write (chns,'(i3.3)') is
!         wrk0=snast(:,0,is)+              !srv 07.10.09
!     &        snast(:,1,is)*na(:,is)
!         call my_out_us(70,nCv,0,wrk0,'b2srst_snast'//chns)
!       enddo
!      endif
      if (switch%b2npht_iout.eq.1 .or. 
     &                              switch%iout_b2wdat.eq.4)   then
        wrk0=shest(:,0)+shest(:,1)*te+  !srv 07.10.09
     &       shest(:,2)*ne + shest(:,3)*te*ne
        call my_out_us(70,nCv,0,wrk0,'b2srst_shest')
        wrk0=shist(:,0)+shist(:,1)*ti+
     &       shist(:,2)*ni(:,0)+shist(:,3)*ti*ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2srst_shist')
      endif
!      if (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1) then
!        wrk0=schst(:,0)+schst(:,1)*po+  !srv 07.10.09
!     &       schst(:,2)*ne+schst(:,3)*po*ne(:)
!        call my_out_us(70,nCv,0,wrk0,'b2srst_schst')
!      endif

      if (switch%b2srst_iout.eq.1) then                       !srv 11.09.09 {
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          do k=0,3
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,smost(1,k,is),
     &                       'b2srst_smost'//chk//chns)
          enddo
          do k=0,1
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,snast(1,k,is),
     &                       'b2srst_snast'//chk//chns)
          enddo
        enddo
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,shest(1,k),'b2srst_shest'//chk)
          call my_out_us(70,nCv,0,shist(1,k),'b2srst_shist'//chk)
          call my_out_us(70,nCv,0,schst(1,k),'b2srst_schst'//chk)
        enddo
      endif                                                 !srv 11.09.09 }

* ..return
      ncall_b2srst = ncall_b2srst+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srst

      end subroutine b2srst
