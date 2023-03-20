      subroutine b2tvsq (nCv, nFc, nVx, ns, switch, geo, mpg, fac_vis,
     &       ne2, rz2, lnlam, te, ti, ni, hcix_c, fchvisq)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2ttia
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   fac_vis(nFc), ne2(nCv), lnlam(nCv), te(nCv), ti(nCv), ni(nCv),
     &   rz2(nCv,0:ns-1), hcix_c(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchvisq(nFc,0:1)                                               !srv 11.07.99
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ic, iCv, iFc, iFt, is, k
      integer, save :: nsmooth = 2
      integer, save :: ncall = 0
      real (kind=R8) ::
     *    fac_vism, wrkc(nCv), wrkf(nFc), wrkv(nVx), wght(nFc,2),       !srv 11.07.99
     *    qip0(nCv), qip1(nCv), qipgen(nCv), gqipgen(nFc),              !srv 20.12.99
     *    tauia1(nCv), alpha1(nCv), nu1(nCv), tauia(nCv,0:ns-1),        !srv 28.03.99 11.07.99
     *    gti(nCv,0:1), gOnedBsq(nFc,0:1)
*   ..procedures
      intrinsic abs
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvsq')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      
      fchvisq   = 0.0e0_R8                                              !srv 11.07.99
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. switch%fhe_vis_q.ne.0.0e0_R8.and.     !srv 16.10.17
     &   switch%no_current.eq.0) then

*     ..compute gradients of ti in cell centers
        call gradc(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

        do is = 0,ns-1                                                  !srv 12.05.11
         if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then            !srv 12.05.11

          call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)             !srv 01.07.09

*       ..calculate the parallel heat flux of ion at the center of cell !lk 23.05.06 05.07.07 {

          qip0 = 0.0e0_R8

          if (switch%b2siav_qstyle .eq. 0) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -2.5e0_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*geo%ftBbav2(iFt))*
     &           (ni(iCv)*ti(iCv)/qe)*gti(iCv,1)
              enddo
             endif
            enddo
            !to ensure that subsequent steps lead to well-defined value on separatrix:
            !copy value of qip0 across separatrix
            !wdk todo: rationalize this when separatrix-structure is introduced?
            do iFc = 1, nFc
              if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
              endif
              if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                 mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
              endif
            enddo

            !smoothing: interpolate to faces and back to cell centers
            wght = 1.0_R8
            do k = 1, nsmooth
             call intface(nCv,nFc,mpg%fcCv,wght,qip0,wrkf)
             call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf,qip0)
             !copy value of qip0 across separatrix
             !wdk todo: rationalize this when separatrix-structure is introduced?
             do iFc = 1, nFc
               if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &             .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                 qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
               endif
               if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                  mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                 qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
               endif
             enddo
            enddo

          endif

          if (switch%b2siav_qstyle.eq.1) then

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

            !to ensure that subsequent steps lead to well-defined value on separatrix:
            !copy value of qip0 across separatrix
            !wdk todo: rationalize this when separatrix-structure is introduced?
            do iFc = 1, nFc
              if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
              endif
              if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                 mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
              endif
            enddo

          endif

*       ..compute qip1
          qip1 = switch%cqip1*ni*ti*sqrt(te+ti)/sqrt(am(is)*mp)

*       ..compute qipgen and its poloidal gradient on faces
          qipgen = qip1*datan(qip0/qip1)*sqrt(geo%cvBb(:,3))                                                                  
          call grad_p(nCv, nFc, nVx, 0, geo, mpg, qipgen, wrkv, gqipgen)

          !set gqipgen to zero on faces not in touch with core
          !wdk todo: rationalize this when separatrix-structure is introduced?
          do iFc = 1, nFc
            if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &          .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
              gqipgen(iFc) = 0.0_R8
            endif
          enddo

*       ..compute tau, nu, alpha, and interpolate to faces
          tauia1 = 0.0e0_R8
          nu1    = 0.0e0_R8
          alpha1 = 0.0e0_R8
          do iFt = 1, mpg%nFt 
            if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
                iCv = mpg%ftCv(ic)
                if (switch%lluciani.eq.3) then
                  tauia1(iCv) = switch%cvsa_mltpl*tauia(iCv,is)/        !lk 20.11.07 { iyv 06.12.13
     &             (1.0e0_R8+15.12e16_R8/geo%ftConn(iFt)*
     &             (ti(iCv)/ev)**2/ne2(iCv))
     &             /(1.0e0_R8+15.12e16_R8*geo%ftEps(iFt)**1.5e0_R8
     &             /geo%ftConn(iFt)*(ti(iCv)/ev)**2/ne2(iCv))
                elseif (switch%lluciani.eq.1) then
                  tauia1(iCv) = tauia(iCv,is)/
     1             (1.0_R8+7.5e16_R8/geo%ftConn(iFt)*
     2             (ti(iCv)/ev)**2/ne2(iCv))
                else
                  tauia1(iCv) = tauia(iCv,is)
                endif
                nu1(iCv) = geo%ftConn(iFt)/(tauia(iCv,is)
     &           *sqrt(ti(iCv)/(mp*am(is)))
     &           *2.0e0_R8*pi*sqrt(geo%ftEps(iFt)**3))
                alpha1(iCv) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(iCv))+
     &           2.7e0_R8*nu1(iCv)**2*geo%ftEps(iFt)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu1(iCv))+
     &           nu1(iCv)**2*geo%ftEps(iFt)**3)-1.0e0_R8)*0.533e0_R8    !lk 20.11.07 }
              enddo
            endif
          enddo

          !to ensure that subsequent steps lead to well-defined value on separatrix:
          !copy values of tauia1, nu1, alpha1 across separatrix
          !wdk todo: rationalize this when separatrix-structure is introduced?
          do iFc = 1, nFc
            if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &          .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
              tauia1(mpg%fcCv(iFc,2)) = tauia1(mpg%fcCv(iFc,1))
              alpha1(mpg%fcCv(iFc,2)) = alpha1(mpg%fcCv(iFc,1))
            endif
            if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &               mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
              tauia1(mpg%fcCv(iFc,1)) = tauia1(mpg%fcCv(iFc,2))
              alpha1(mpg%fcCv(iFc,1)) = alpha1(mpg%fcCv(iFc,2))
            endif
          enddo

          wrkc = tauia1*alpha1
          call intface(nCv,nFc,mpg%fcCv,wght,wrkc,wrkf)
 
*       ..add prefactor from magnetic field
          wrkf = wrkf*geo%fcBb(:,2)*geo%fcBb(:,0)/(geo%fcBb(:,3)*
     &      sqrt(geo%fcBb(:,3)))

*       ..calculate the x-component of the parallel viscosity current
*         connected with the heat flux 
          fchvisq(:,0) = fchvisq(:,0) - 0.24e0_R8*switch%fhe_vis_q*
     &       fac_vis*wrkf*gOnedBsq(:,1)*gqipgen*geo%fcS*geo%fcQalf(:,0)
          fchvisq(:,1) = fchvisq(:,1) + 0.24e0_R8*switch%fhe_vis_q*
     &       fac_vis*wrkf*gOnedBsq(:,0)*gqipgen*geo%fcS*geo%fcQalf(:,1)

         endif
        enddo
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvsq

      end subroutine b2tvsq
