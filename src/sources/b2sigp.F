*-----------------------------------------------------------------------
*.specification

      subroutine b2sigp (nCv, nFc, nVx, isb, geo, mpg, rzb, nb, ti,
     &  ne, te, po, smbgp, wrk0, wrk1, wrk2)
      use b2mod_types
      use b2mod_constants                                         !srv 23.03.10
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance                                           !djm Jan2017
!WG_TODO     & , only : b2sigp_smogp0to3, b2sigp_pstat0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rzb(nCv), 
     &  nb(nCv), ti(nCv),
     &  ne(nCv), te(nCv), po(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbgp(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     * wrk0(nCv), wrk1(nCv), wrk2(nCv),
     * wrk3(nCv), wrk4(nCv),                      !djm Jan2017
     * wrk5(nCv), wrk6(nCv)                       !srv 28.03.17

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIGP computes the pressure gradient term on the right hand
*     of the parallel momentum balance equation for one species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, k
      integer, save :: style = 2                                     !srv 13.01.17
      integer, save :: iout = 0, iout_b2npmo = 0, iout_b2wdat = 0
      integer, save :: pot_eq = 1, no_current = 0  !srv 17.06.02 11.09.09 23.03.10
      character chns*3 !, chk*1                                       !srv 17.06.02
      real (kind=R8) :: phm0 = 1.0_R8                                 !srv 26.06.18
      real (kind=R8) ::
     *   t0, wrkf(nFc), wrkv(nVx)
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg
      intrinsic abs
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sigp')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sigp_iout', iout)                             !srv 17.06.02
       call ipgetr ('b2sigp_phm0', phm0)                             !srv 26.06.18
       call ipgeti ('b2wdat_iout', iout_b2wdat)                      !iys 17.04.17
       call ipgeti ('b2sigp_style', style)                           !srv 02.07.08
       call ipgeti ('b2npmo_iout', iout_b2npmo)                      !srv 11.09.09
       call ipgeti ('b2news_poteq', pot_eq)                          !srv 23.03.10
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
      endif
*   ..test nCv, nFv
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, ti, ne, te, rzb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, rzb, 1, 'rzb', '.ge.')
      endif

      if (phm0.ne.0.0_R8) then                                     !srv 26.06.18
* ..compute pressure gradient term
*   ..compute smbgp
      if(style.eq.0) then                                            !srv 02.07.08
        call xerrab('b2sigp_style=0 not available for WG')
#ifdef WG_TODO
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           wrk0(ix,iy) = 0.0_R8
           wrk1(ix,iy) = 0.0_R8
           wrk2(ix,iy) = rzb(ix,iy)*nb(ix,iy)/ne(ix,iy)
          else
           wrk0(ix,iy) = nb(ix,iy)*ti(ix,iy)-
     &      nb(leftix(ix,iy),leftiy(ix,iy))*
     &      ti(leftix(ix,iy),leftiy(ix,iy))
           wrk1(ix,iy) = ne(ix,iy)*te(ix,iy)-
     &      ne(leftix(ix,iy),leftiy(ix,iy))*
     &      te(leftix(ix,iy),leftiy(ix,iy))
           wrk2(ix,iy) = ((rzb(ix,iy)*nb(ix,iy)/ne(ix,iy))*
     &                  vol(leftix(ix,iy),leftiy(ix,iy))+
     &                  (rzb(leftix(ix,iy),leftiy(ix,iy))*
     &                   nb(leftix(ix,iy),leftiy(ix,iy))/
     &                   ne(leftix(ix,iy),leftiy(ix,iy)))*vol(ix,iy))/
     &                 (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
          endif
          wrk3(ix,iy) = nb(ix,iy)*ti(ix,iy)                          !djm Jan2017
          wrk4(ix,iy) = ne(ix,iy)*te(ix,iy)*wrk2(ix,iy)
         enddo
        enddo
        call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
        call b2xxgs (nx, ny, 0.0_R8, wrk1, 1)
        call b2xxgs (nx, ny, 0.0_R8, wrk2, 1)
*   ..  compute smbgp
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
              wrk5(ix,iy) = -hz(ix,iy)*
     &         (vol(ix,iy)/hx(ix,iy)/2.0_R8)*                         !srv 16.06.09
     &          bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                   !srv 11.01.13
     &         (wrk0(rightix(ix,iy),rightiy(ix,iy)))/
     &         (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
              wrk6(ix,iy) = -hz(ix,iy)*
     &         (vol(ix,iy)/hx(ix,iy)/2.0_R8)*                         !srv 16.06.09
     &          bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                   !srv 11.01.13
     &          wrk1(rightix(ix,iy),rightiy(ix,iy))*
     &          wrk2(rightix(ix,iy),rightiy(ix,iy))/
     &         (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
              smbgp(ix,iy,0) = wrk5(ix,iy) + wrk6(ix,iy)
              smbgp(ix,iy,1) = 0.0_R8
              smbgp(ix,iy,2) = 0.0_R8
              smbgp(ix,iy,3) = 0.0_R8
              if (balance_netcdf.ne.0) then !djm Jan2017
!     djm Mar2017 NB - THE FOLLOWING DOES NOT GIVE EXACT BALANCE! FOR EXACT PRESSURE BALANCE IN
!                      THIS REGION, USE B2SIGP_STYLE = 1
                if (ncall.lt.3) then
                 write(*,'(a,i4,a,i4)')
     &       'warning, tot. press. balance will not be exact at ix=',
     &            ix,',iy=',iy
                 write(*,*) 'Use b2sigp_style.eq.1 or 2'
                end if
                b2sigp_pstat0to3(ix,iy,0,isb) = hz(ix,iy)*
     &           (vol(ix,iy)/hx(ix,iy))*
     &           bb(ix,iy,0)/bb(ix,iy,3)*
     &           (wrk3(ix,iy)+wrk4(ix,iy))
                b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
              endif
            elseif(rightix(ix,iy).eq.nx+1) then
              wrk5(ix,iy) = -hz(ix,iy)*
     &         (vol(ix,iy)/hx(ix,iy)/2.0_R8)*                         !srv 16.06.09
     &         bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                    !srv 11.01.13
     &         wrk0(ix,iy)/
     &         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
              wrk6(ix,iy) = -hz(ix,iy)*
     &         (vol(ix,iy)/hx(ix,iy)/2.0_R8)*                         !srv 16.06.09
     &         bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                    !srv 11.01.13
     &         wrk1(ix,iy)*wrk2(ix,iy)/
     &         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
              smbgp(ix,iy,0) = wrk5(ix,iy) + wrk6(ix,iy)
              smbgp(ix,iy,1) = 0.0_R8
              smbgp(ix,iy,2) = 0.0_R8
              smbgp(ix,iy,3) = 0.0_R8
              if (balance_netcdf.ne.0) then !djm Jan2017
                b2sigp_pstat0to3(ix,iy,0,isb) = hz(ix,iy)*
     &           (vol(ix,iy)/hx(ix,iy))*
     &           bb(ix,iy,0)/bb(ix,iy,3)*
     &           ((wrk3(leftix(ix,iy),leftiy(ix,iy))+
     &             wrk4(leftix(ix,iy),leftiy(ix,iy)))*vol(ix,iy)+
     &            (wrk3(ix,iy)+
     &             wrk4(ix,iy))*vol(leftix(ix,iy),leftiy(ix,iy)))/
     &           (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
                b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
              endif
            else
              wrk5(ix,iy) = -hz(ix,iy)*(vol(ix,iy)/hx(ix,iy))*        !srv 16.06.09
     &         bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                    !srv 11.01.13
     &        (wrk0(ix,iy)/
     &         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &         wrk0(rightix(ix,iy),rightiy(ix,iy))/
     &         (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy))))
              wrk6(ix,iy) = -hz(ix,iy)*(vol(ix,iy)/hx(ix,iy))*        !srv 16.06.09
     &         bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)*                    !srv 11.01.13
     &        (wrk1(ix,iy)*wrk2(ix,iy)/
     &         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &         wrk1(rightix(ix,iy),rightiy(ix,iy))*
     &          wrk2(rightix(ix,iy),rightiy(ix,iy))/
     &         (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy))))
              smbgp(ix,iy,0) = wrk5(ix,iy) + wrk6(ix,iy)
              smbgp(ix,iy,1) = 0.0_R8
              smbgp(ix,iy,2) = 0.0_R8
              smbgp(ix,iy,3) = 0.0_R8
              if (balance_netcdf.ne.0) then !djm Jan2017
                b2sigp_pstat0to3(ix,iy,0,isb) = hz(ix,iy)*
     &           (vol(ix,iy)/hx(ix,iy))*
     &           bb(ix,iy,0)/bb(ix,iy,3)*
     &           ((wrk3(leftix(ix,iy),leftiy(ix,iy))+
     &             wrk4(leftix(ix,iy),leftiy(ix,iy)))*vol(ix,iy)+
     &            (wrk3(ix,iy)+
     &             wrk4(ix,iy))*vol(leftix(ix,iy),leftiy(ix,iy)))/
     &           (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
                b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
                b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
              endif
            endif
          enddo
        enddo
#endif
      elseif(style.eq.1) then                                         !srv 02.07.08 {
*   ..  compute gradients on faces and interpolate to centers
        wrk2 = nb*ti
        wrk3 = ne*te
        call grad_p(nCv, nFc, nVx, 0, geo, mpg, wrk2, wrkv, wrkf)
        call intcell(nFc,nCv, mpg, mpg%intcellP, wrkf, wrk0)
        call grad_p(nCv, nFc, nVx, 0, geo, mpg, wrk3, wrkv, wrkf)
        call intcell(nFc,nCv, mpg, mpg%intcellP, wrkf, wrk1)
* ..compute pressure gradient term
        smbgp = 0.0e0_R8
        wrk5  = 0.0e0_R8
        wrk6  = 0.0e0_R8
*   ..  compute smbgp
        do iCv = 1, mpg%nCi
         t0=geo%cvHz(iCv)*geo%cvVol(iCv)*geo%cvBb(iCv,0)/geo%cvBb(iCv,3)           !srv 18.10.02 10.06.05
         wrk5(iCv) = -t0*wrk0(iCv)
         wrk6(iCv) = -t0*rzb(iCv)*(nb(iCv)/ne(iCv))*wrk1(iCv)
         smbgp(iCv,0) = wrk5(iCv) + wrk6(iCv)
        enddo
cdjm compute b2sigp_pstat0to3 on left cell face
#ifdef WG_TODO
        if (balance_netcdf.ne.0) then
          b2sigp_pstat0to3 = 0.0_R8
          do iy = -1, ny
           do ix = -1, nx
            if(leftix(ix,iy).eq.-2.or.
     &         bottomiy(ix,iy).eq.-2.or.topiy(ix,iy).eq.ny+1) then
              b2sigp_pstat0to3(ix,iy,0,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
            else
              t0=hz(ix,iy)*vol(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)
              b2sigp_pstat0to3(ix,iy,0,isb) =
     &         t0*((wrk2(leftix(ix,iy),leftiy(ix,iy))+
     &              rzb(leftix(ix,iy),leftiy(ix,iy))*
     &              (nb(leftix(ix,iy),leftiy(ix,iy))/
     &              ne(leftix(ix,iy),leftiy(ix,iy)))*
     &              wrk3(leftix(ix,iy),leftiy(ix,iy)))*hx(ix,iy)+
     &             (wrk2(ix,iy)+rzb(ix,iy)*(nb(ix,iy)/ne(ix,iy))*
     &              wrk3(ix,iy))*hx(leftix(ix,iy),leftiy(ix,iy)))/
     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))/hx(ix,iy)
              b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
            endif
           enddo
          enddo
        endif
#endif
      else                                                            !srv 23.03.10 {
        call xertst(pot_eq.eq.1,
     &   'b2sigp_style.eq.2 requires pot_eq.eq.1!')
* ..compute pressure gradient term
        smbgp = 0.0e0_R8
        wrk5  = 0.0e0_R8
        wrk6  = 0.0e0_R8
*   ..  compute gradients on faces and interpolate to centers
        call grad_p(nCv, nFc, nVx, 0, geo, mpg, nb*ti, wrkv, wrkf)
        call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf,wrk0)
        call grad_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, wrkf)
        call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf,wrk1)
*   ..  compute smbgp
        do iCv = 1, mpg%nCi
         t0=geo%cvHz(iCv)*geo%cvVol(iCv)*geo%cvBb(iCv,0)/geo%cvBb(iCv,3)
         wrk5(iCv) = -t0*wrk0(iCv)
         wrk6(iCv) = -t0*qe*rzb(iCv)*nb(iCv)*wrk1(iCv)
         smbgp(iCv,0) = wrk5(iCv) + wrk6(iCv)
        enddo                                                         !srv 23.03.10 }
cdjm compute b2sigp_pstat0to3 on left cell face
#ifdef WG_TODO
        if (balance_netcdf.ne.0) then
          b2sigp_pstat0to3 = 0.0_R8
          do iy = -1, ny
           do ix = -1, nx
            if(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1 .or.
     &         bottomiy(ix,iy).eq.-2.or.topiy(ix,iy).eq.ny+1) then
              b2sigp_pstat0to3(ix,iy,0,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
            else
              t0=hz(ix,iy)*vol(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)
              b2sigp_pstat0to3(ix,iy,0,isb) =
     &         -t0*wrk0(ix,iy)
     &         -t0*qe*wrk1(ix,iy)*
     &         (rzb(ix,iy)*nb(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+
     &          rzb(leftix(ix,iy),leftiy(ix,iy))*
     &          nb(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/
     &         (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
              b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
              b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
            endif
           enddo
          enddo
        endif
#endif
      endif                                                           !srv 02.07.08 }

cxpb apply restriction to pressure gradient for minority ions
#ifdef WG_TODO
      do iy = -1, ny                                                  !lk 23.06.17{  
       do ix = -1, nx
        if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1 .and.
     &     nb(ix,iy)/ne(ix,iy) .lt. 0.1e0_R8) then
          smbgp(ix,iy,0) = smbgp(ix,iy,0)*min(1.0e0_R8, 
     &    (2.0e0_R8*hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy))
     &                       +hx(rightix(ix,iy),rightiy(ix,iy)))*
     &     nb(ix,iy)*ti(ix,iy)/
     &    ((hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*
     &     nb(rightix(ix,iy),rightiy(ix,iy))
     &    *ti(rightix(ix,iy),rightiy(ix,iy))+
     &    (hx(ix,iy)+hx(rightix(ix,iy),rightiy(ix,iy)))*
     &     nb(leftix(ix,iy),leftiy(ix,iy))
     &    *ti(leftix(ix,iy),leftiy(ix,iy))))     
        endif   
       enddo
      enddo                                                           !lk 23.06.17}
#endif
      else                                                            !srv 26.06.18 {
        smbgp = 0.0_R8
#ifdef WG_TODO
        b2sigp_pstat0to3 = 0.0_R8
#endif
        wrk5 = 0.0_R8
        wrk6 = 0.0_R8       
      endif                                                           !srv 26.06.18 }
      
      if (iout.ne.0 .or. iout_b2wdat.eq.4) then                       !srv 17.06.02 {
       write (chns,'(i3.3)') isb
!      do k=0,3
!        write (chk,'(i1)') k
!        call my_out_us(70,nCv,0,smbgp(-1,-1,k),'b2sigp_smogp'//chk//chns)
!      enddo
       call my_out_us(70,nCv,0,smbgp(1,0),'b2sigp_smogp'//chns)
       if(style.eq.0 .or. style.eq.1 ) then
          call my_out_us(70,nCv,0,wrk5,'b2sigp_smogpi'//chns)
          call my_out_us(70,nCv,0,wrk6,'b2sigp_smogpe'//chns)
       elseif(style.eq.2) then
          call my_out_us(70,nCv,0,wrk5,'b2sigp_smogpi'//chns)
          call my_out_us(70,nCv,0,wrk6,'b2sigp_smogpo'//chns)
       endif
      endif                                                           !srv 17.06.02 }

cdjm Jan2017 Store linearised source for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2sigp_smogp0to3(-1:nx,-1:ny,0:3,isb)=smbgp
!WG_TODO      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sigp

      end subroutine b2sigp
