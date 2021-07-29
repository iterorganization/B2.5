*-----------------------------------------------------------------------
*.specification

      subroutine b2siav (nCv, nFc, nVx, ns, isb, switch, geo, mpg, mb,
     & ub, te, ti, nb, ni, ne2, rz2, lnlam, hcix_c, cvsbhz_cl, smbvh,
     & smbvv, kt_neo, wrk0, wrk1, wrk2)                                 !srv 23.07.99 04.10.99 17.01.01 09.01.01 31.07.02 01.07.05
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance                                                 !djm Jan2017
!WG_TODO     & , only : b2siav_smovh0to3, b2siav_smovv0to3, balance_netcdf
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2siav
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     & mb, lnlam(nCv),                                                  !srv 01.07.09
     & ub(nCv),
     & te(nCv), ti(nCv), nb(nCv), ni(nCv), ne2(nCv),
     & hcix_c(nCv), cvsbhz_cl(nFc,0:1), rz2(nCv,0:ns-1)                 !srv 31.07.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbvh(nCv,0:3), smbvv(nCv,0:3),
     &   kt_neo(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(nCv), wrk1(nCv), wrk2(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAV computes a linearised representation of the additional
*     viscosity.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ic, iCv, iFc, iFt, k
      integer, save :: qstyle = 0, lluciani = 3                         !lk 03.12.09
      integer, save :: iout = 0                                         !srv 17.06.02
      integer, save :: nsmooth = 2
      character chns*3, chk*1                                           !srv 17.06.02
      real (kind=R8) ::
     *     t0, ft                                                       !srv 16.10.17
      real (kind=R8), save ::
     &     addvis=0.0_R8, addvis1=1.0_R8, cqip1 = 0.5e0_R8,             !srv 28.03.06 !srv 19.10.08
     &     cvsa_mltpl = 1.0e0_R8                                        !iyv 06.12.13
      real (kind=R8) ::
     &     wrks(nCv), wrkf(nFc), wrkv(nVx), wght(nFc,2),
     *     tauia(nCv,0:ns-1), qip0(nCv),                                !srv 23.07.99
     *     qip1(nCv), qipgen(nCv),                                      !srv 23.07.99 19.01.05
     &     tt0(mpg%nFt), tt1(mpg%nFt), bbaverage2(mpg%nFt),             !lk 23.05.06 05.07.07 {
     &     tauia1(nCv), alpha1(nCv), nu1(nCv), s2,                      !lk !iyv 02.12.13
     &     gti(nCv,0:1)
      integer mintopcut                                                 !lk 23.05.06 05.07.07 }
*   ..procedures
      intrinsic abs, sqrt
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2siav')
*   ..set internal parameters on first call
      if (ncall_b2siav.eq.0) then
       call ipgetr ('b2siav_addvis', addvis)                            !srv 23.07.99
       call ipgetr ('b2siav_addvis1', addvis1)                          !srv 04.10.99
       call ipgetr ('b2siav_cqip1', cqip1)                              !srv 28.03.06
       call ipgeti ('b2siav_iout', iout)                                !lk 23.05.06 05.07.07
       call ipgeti ('b2siav_style_qip', qstyle)                         !lk 23.05.06 05.07.07
       call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)	                !iyv 06.12.13
      endif
      ft = 4.0e0_R8/3.0e0_R8                                            !srv 16.10.17
      s2 = sqrt(2.0_R8)                                                 !iyv 02.12.13
      smbvh = 0.0_R8
      smbvv = 0.0_R8

      if(nint(zn(isb)).eq.1 .and. nint(zamax(isb)).eq.1) then           !srv 26.10.08

        if (addvis.ne.0.0_R8) then                                      !srv 19.01.05 28.03.06 {

          call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)             !srv 01.07.09

*       ..compute radial gradients of ti in cell centers
          call gradc(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)

* ..calculate the parallel heat flux of ions in the center of the cell

          tt0 = 0.0e0_R8
          tt1 = 0.0e0_R8
          bbaverage2 = 0.0e0_R8
          do iFt = 1, mpg%nFt 
           if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
            do ic = mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
             iCv = mpg%ftCv(ic)
             tt0(iFt) = tt0(iFt) + geo%cvBb(iCv,3)*geo%cvBb(iCv,3)*
     .         geo%cvVol(iCv)
             tt1(iFt) = tt1(iFt) + geo%cvVol(iCv)
            enddo
            bbaverage2(iFt) = tt0(iFt)/tt1(iFt)
           endif
          enddo

          qip0 = 0.0e0_R8

          if (qstyle .eq. 0) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -2.5e0_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*bbaverage2(iFt))*(ni(iCv)*ti(iCv)/qe)*
     &           gti(iCv,1)
              enddo
             endif
            enddo

            !smoothing: interpolate to faces and back to cell centers
            wght = 1.0_R8
            do k = 1, nsmooth
             call intface(nCv,nFc,mpg%fcCv,wght,qip0,wrkf)
             call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf,qip0)
            enddo

          endif

          if (qstyle.eq.1) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -hcix_c(iCv)*gti(iCv,0)*geo%cvBb(iCv,3)/
     &           geo%cvBb(iCv,0) - 2.5e0_R8*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*geo%cvBb(iCv,3))*(ni(iCv)*ti(iCv)/qe)*
     &           gti(iCv,1)
              enddo
             endif
            enddo

          endif

          tauia1 = 0.0e0_R8
          nu1    = 0.0e0_R8
          alpha1 = 0.0e0_R8
          kt_neo = 0.0e0_R8                                             !iyv 02.12.13
          do iFt = 1, mpg%nFt 
            if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
                iCv = mpg%ftCv(ic)
                if (lluciani.eq.3) then
                  tauia1(iCv) = cvsa_mltpl*tauia(iCv,isb)/              !lk 20.11.07 { iyv 06.12.13
     &             (1.0e0_R8+15.12e16_R8/geo%ftConn(iFt)*
     &             (ti(iCv)/ev)**2/ne2(iCv))
     &             /(1.0e0_R8+15.12e16_R8*geo%ftEps(iFt)**1.5e0_R8
     &             /geo%ftConn(iFt)*(ti(iCv)/ev)**2/ne2(iCv))
                elseif (lluciani.eq.1) then
                  tauia1(iCv) = tauia(iCv,isb)/
     1             (1.0_R8+7.5e16_R8/geo%ftConn(iFt)*
     2             (ti(iCv)/ev)**2/ne2(iCv))
                else
                  tauia1(iCv) = tauia(iCv,isb)
                endif
                nu1(iCv) = geo%ftConn(iFt)/(tauia(iCv,isb)
     &           *sqrt(ti(iCv)/(mp*am(isb)))
     &           *2.0e0_R8*pi*sqrt(geo%ftEps(iFt)**3))
                alpha1(iCv) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(iCv))+
     &           2.7e0_R8*nu1(iCv)**2*geo%ftEps(iFt)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu1(iCv))+
     &           nu1(iCv)**2*geo%ftEps(iFt)**3)-1.0e0_R8)*0.533e0_R8   !lk 20.11.07 }
                kt_neo(iCv) =
     &           ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(iCv)*s2)               !iyv 02.12.13
     &           +2.7e0_R8*(nu1(iCv)*s2)**2*geo%ftEps(iFt)**3)/
     &           (1.0e0_R8+0.7e0_R8*sqrt(nu1(iCv)*s2)+
     &           (nu1(iCv)*s2)**2*geo%ftEps(iFt)**3))
              enddo
            endif
          enddo

*       ..compute qip1
          qip1=cqip1*ni*ti*sqrt(te+ti)/sqrt(am(isb)*mp)

*       ..compute qipgen
          qipgen=qip1*datan(qip0/qip1)*sqrt(geo%cvBb(:,3))

*       ..gradient qipgen in cell centers   
          !wdk modified computation compared to structured code to 
          !wdk facilitate implementation on unstructured grids
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, qipgen, wrkv, wrk0)
          wrk0 = wrk0*geo%cvBb(:,0)/geo%cvBb(:,3)**3*tauia1*alpha1
          
*       ..finally calculate the r.h.s. contribution to the parallel balance
*         corresponding to parallel viscosity connected with heat flow
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk1)
          smbvh(:,0) = addvis*geo%cvVol*geo%cvBb(:,0)*
     &      sqrt(geo%cvBb(:,3))*geo%cvHz*wrk1

*       ..zero out in guard cells
          smbvh(mpg%nCi+1:nCv,0) = 0.0_R8

        endif                                                           !lk 23.05.06 05.07.07

        if (addvis1.ne.0.0_R8) then                                     !srv 19.01.05 30.01.10

*       ..prepare wrk arrays
          wrk2 = sqrt(geo%cvBb(:,3))                                    ! wrk2 is used below twice
          wrk1 = log(geo%cvHz*wrk2)

*       ..poloidal differences of wrk1 on cell faces
          wrkv = log(geo%vxHz*sqrt(geo%vxBb(:,3)))
          call diff_p(nCv, nFc, nVx, 1, geo, mpg, wrk1, wrkv, wrkf)

*       ..compute flow and its divergence
          wrkf = wrkf*ft*cvsbhz_cl(:,0)/sqrt(geo%fcBb(:,3))
          call div0(nCv, nFc, mpg, wrkf, wrk0)

*       ..add to r.h.s. source
          do iCv = 1, mpg%nCi
            t0 = wrk0(iCv)*wrk2(iCv)                                    !srv 01.07.05 }
            if (t0.gt.0.0e0_R8) then
              smbvv(iCv,0) = t0*ub(iCv)
            else
              smbvv(iCv,1) = t0
            endif
          enddo                                                         !srv 19.01.05 }

        endif
      endif                                                             !srv 26.10.08

      if (iout.ne.0) then                                               !srv 17.06.02 {
        call my_out_us(70,nCv,0,alpha1,'b2siav_alpha1')
        call my_out_us(70,nCv,0,kt_neo,'b2siav_kt_neo')
        call my_out_us(70,nCv,0,tauia1,'b2siav_tauia1')
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,smbvv(1,k),'b2siav_smovv'//chk//chns)
          call my_out_us(70,nCv,0,smbvh(1,k),'b2siav_smovh'//chk//chns)
        enddo
        wrks=smbvv(:,0)+smbvv(:,1)*ub+
     &       smbvv(:,2)*nb*mb+smbvv(:,3)*ub*nb*mb
        call my_out_us(70,nCv,0,wrks,'b2siav_smovv'//chns)
        wrks=smbvh(:,0)+smbvh(:,1)*ub+
     &       smbvh(:,2)*nb*mb+smbvh(:,3)*ub*nb*mb
        call my_out_us(70,nCv,0,wrks,'b2siav_smovh'//chns)
        call my_out_us(70,nFc,1,cvsbhz_cl(1,0),'b2siav_cvsbhzx'//chns)
        call my_out_us(70,nFc,1,cvsbhz_cl(1,1),'b2siav_cvsbhzy'//chns)
      endif                                                             !srv 17.06.02 }

cdjm Jan2017 Store linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2siav_smovh0to3(-1:nx,-1:ny,0:3,isb)=smbvh
!WG_TODO        b2siav_smovv0to3(-1:nx,-1:ny,0:3,isb)=smbvv
!WG_TODO      endif

* ..return
      ncall_b2siav = ncall_b2siav+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2siav

      end subroutine b2siav
