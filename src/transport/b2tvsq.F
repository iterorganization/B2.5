      subroutine b2tvsq (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &   fac_vis, ne2, rz2, lnlam, te, ti, na, ni, rza, vsaf_hAdp_albe,
     &   vsaf_hBdp_al, hcix_c, hci_a, hci_al_ast, fchvisq_a, fchvisq)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2ttia
#endif
      use b2mod_ad
     & , only : nsmooth_tvsq
      use b2mod_subsys
      use b2mod_zhfrtf
     & , only : nal, Z_to_m1_ast, Ia
      use b2mod_b2mrel
     & , only : amnucl, s2nucl
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   fac_vis(nFc), ne2(nCv), lnlam(nCv), te(nCv), ti(nCv),
     &   na(nCv,0:ns-1), ni(nCv), rza(nCv,0:ns-1), rz2(nCv,0:ns-1),
     &   hcix_c(nCv), hci_a(nCv,0:ns-1), hci_al_ast(nCv,0:nspecies-1),
     &   vsaf_hAdp_albe(nFc,0:1,0:nspecies-1,0:nspecies-1),
     &   vsaf_hBdp_al(nFc,0:1,0:nspecies-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchvisq_a(nFc,0:1,0:ns-1), fchvisq(nFc,0:1)                    !srv 11.07.99
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSQ computes the [poloidal and radial components of the heat
*     viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ic, iCv, iFc, iFt, is, inucl, k
      integer, save :: ncall = 0
      real (kind=R8) ::
     *    fac_vism, wrkc(nCv), wrkf(nFc), wrkv(nVx), wght(nFc,2),       !srv 11.07.99
     *    qip0(nCv), qip1(nCv), qipgen(nCv), gqipgen(nFc),              !srv 20.12.99
     *    gqipgen_al(nFc,0:nspecies-1),
     *    tauia1(nCv), alpha1(nCv), nu1(nCv), tauia(nCv,0:ns-1),        !srv 28.03.99 11.07.99
     *    gti(nCv,0:1), gOnedBsq(nFc,0:1), naf(nFc), Iaf(nFc)
*   ..procedures
      external xertst
      external b2xvsg
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvsq')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      
      fchvisq   = 0.0_R8                                                !srv 11.07.99
      fchvisq_a = 0.0_R8
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. switch%fhe_vis_q.ne.0.0_R8.and.       !srv 16.10.17
     &   switch%no_current.eq.0) then

       wght = 1.0_R8

*    ..compute gradients of ti in cell centers
       call gradc(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)

*    ..compute gradients of OnedBsq on cell faces
       call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &   geo%vxOnedBsq, gOnedBsq)

       if (switch%Zhdanov_closure.eq.0 .or.
     &     switch%Zhdanov_test.eq.1 .or.
     &     switch%Zhdanov_vish.eq.0) then                               !som 20.10.21

        do is = 0,ns-1                                                  !srv 12.05.11
         if(is.eq.ismain) then                                          !srv 12.05.11

          call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)             !srv 01.07.09

*       ..calculate the parallel heat flux of ion at the center of cell !lk 23.05.06 05.07.07 {

          qip0 = 0.0_R8

          if (switch%b2siav_qstyle .eq. 0) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -2.5_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
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
            do k = 1, nsmooth_tvsq
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
     &           geo%cvBb(iCv,0) - 2.5_R8*geo%cvBb(iCv,2)/
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
          tauia1 = 0.0_R8
          nu1    = 0.0_R8
          alpha1 = 0.0_R8
          do iFt = 1, mpg%nFt 
            if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
                iCv = mpg%ftCv(ic)
                if (switch%lluciani.eq.3) then
                  tauia1(iCv) = switch%cvsa_mltpl*tauia(iCv,is)/        !lk 20.11.07 { iyv 06.12.13
     &             (1.0_R8+15.12e16_R8/geo%ftConn(iFt)*
     &             (ti(iCv)/ev)**2/ne2(iCv))
     &             /(1.0_R8+15.12e16_R8*geo%ftEps(iFt)**1.5_R8
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
     &           *2.0_R8*pi*sqrt(geo%ftEps(iFt)**3))
                alpha1(iCv) = ((-0.17_R8+1.05_R8*sqrt(nu1(iCv))+
     &           2.7_R8*nu1(iCv)**2*geo%ftEps(iFt)**3)/(1.0_R8+
     &           0.7_R8*sqrt(nu1(iCv))+
     &           nu1(iCv)**2*geo%ftEps(iFt)**3)-1.0_R8)*0.533_R8        !lk 20.11.07 }
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
          fchvisq_a(:,0,is) = -0.24_R8*switch%fhe_vis_q*
     &       fac_vis*wrkf*gOnedBsq(:,1)*gqipgen*geo%fcS*geo%fcQalf(:,0)
          fchvisq_a(:,1,is) =  0.24_R8*switch%fhe_vis_q*
     &       fac_vis*wrkf*gOnedBsq(:,0)*gqipgen*geo%fcS*geo%fcQalf(:,1)

         endif

*      ..accumulate individual contributions into fchvisq
         fchvisq = fchvisq + fchvisq_a(:,:,is)

        enddo

       else ! Zhdanov closure

        do is = 0,ns-1                                                  !srv 12.05.11
         if (.not.is_neutral(is)) then

*       ..calculate the parallel heat flux of ion at the center of cell !lk 23.05.06 05.07.07 {

          qip0 = 0.0_R8

          if (switch%b2siav_qstyle .eq. 0) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -2.5_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*geo%ftBbav2(iFt))*
     &           (na(iCv,is)*ti(iCv)/(qe*rza(iCv,is)))*gti(iCv,1)
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
            do k = 1, nsmooth_tvsq
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
               qip0(iCv) = -hci_a(iCv,is)*gti(iCv,0)*geo%cvBb(iCv,3)/
     &           geo%cvBb(iCv,0) - 2.5_R8*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*geo%cvBb(iCv,3))*
     &           (na(iCv,is)*ti(iCv)/(qe*rza(iCv,is)))*gti(iCv,1)
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
          qip1 = switch%cqip1*na(:,is)*ti*sqrt(te+ti)/sqrt(am(is)*mp)

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

*       ..calculate the parallel charge state averaged heat flux of all species
          do inucl= 0, nspecies-1

            qip0 = 0.0e0_R8

            if (switch%b2siav_qstyle.eq.0) then

              do iFt = 1, mpg%nFt 
               if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1))))
     &           then
                do ic=mpg%ftCvP(iFt,1),
     &                mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
                 iCv = mpg%ftCv(ic)
                 qip0(iCv) = -2.5e0_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
     &             (geo%cvBb(iCv,0)*geo%ftBbav2(iFt))*
     &             (nal(iCv,inucl)*ti(iCv)/qe*Z_to_m1_ast(iCv,inucl))*
     &             gti(iCv,1)
                enddo
               endif
              enddo

              !to ensure that subsequent steps lead to well-defined value on separatrix:
              !copy value of qip0 across separatrix
              !wdk todo: rationalize this when separatrix-structure is introduced?
              do iFc = 1, nFc
                if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &              .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                  qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
                endif
                if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                   mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                  qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
                endif
              enddo

              !smoothing: interpolate to faces and back to cell centers
              do k = 1, nsmooth_tvsq
               call intface(nCv,nFc,mpg%fcCv,wght,qip0,wrkf)
               call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf,qip0)
               !copy value of qip0 across separatrix
               !wdk todo: rationalize this when separatrix-structure is introduced?
               do iFc = 1, nFc
                 if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &               .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                   qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
                 endif
                 if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                    mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                   qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
                 endif
               enddo
              enddo

            endif

            if (switch%b2siav_qstyle.eq.1)  then

              do iFt = 1, mpg%nFt 
               if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1))))
     &           then
                do ic=mpg%ftCvP(iFt,1),
     &                mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
                 iCv = mpg%ftCv(ic)
                 qip0(iCv) = -hci_al_ast(iCv,inucl)*gti(iCv,0)*
     &             geo%cvBb(iCv,3)/geo%cvBb(iCv,0) - 2.5e0_R8*
     &             geo%cvBb(iCv,2)/(geo%cvBb(iCv,0)*geo%cvBb(iCv,3))*
     &             (nal(iCv,inucl)*ti(iCv)/qe*Z_to_m1_ast(iCv,inucl))*
     &            gti(iCv,1)
                enddo
               endif
              enddo

              !to ensure that subsequent steps lead to well-defined value on separatrix:
              !copy value of qip0 across separatrix
              !wdk todo: rationalize this when separatrix-structure is introduced?
              do iFc = 1, nFc
                if (     mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &              .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                  qip0(mpg%fcCv(iFc,2)) = qip0(mpg%fcCv(iFc,1))
                endif
                if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &                   mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                  qip0(mpg%fcCv(iFc,1)) = qip0(mpg%fcCv(iFc,2))
                endif
              enddo

            endif

*         ..compute qip1
            qip1 = switch%cqip1*nal(:,inucl)*ti*sqrt(te+ti)/
     &        sqrt(amnucl(inucl)*mp)

*         ..compute qipgen and its poloidal gradient on faces
            qipgen = qip1*datan(qip0/qip1)*
     &        sqrt(geo%cvBb(:,3))
            call grad_p(nCv, nFc, nVx, 0, geo, mpg, qipgen, wrkv,
     &        gqipgen_al(:,inucl))

            !set gqipgen_al to zero on faces not in touch with core
            !wdk todo: rationalize this when separatrix-structure is introduced?
            do iFc = 1, nFc
              if (.not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            .not.mpg%cvOnClosedSurface(mpg%fcCv(iFc,2))) then
                gqipgen_al(iFc, inucl) = 0.0_R8
              endif
            enddo

          enddo

*       ..interpolate pre-factors to faces
          wrkc = na(:,is)/nal(:,s2nucl(is)-1)
          call intface(nCv,nFc,mpg%fcCv,wght,wrkc,naf)
          call intface(nCv,nFc,mpg%fcCv,wght,Ia(:,is),Iaf)

          wrkf = 0.0_R8
          do inucl = 0, nspecies-1
            wrkf = wrkf + vsaf_hAdp_albe(:,0,s2nucl(is)-1,inucl)*
     &        gqipgen_al(:,inucl)
          enddo
          if (switch%zhcscorr.gt.0 .and. switch%zhcscorr_vq.gt.0) then
            wrkf = wrkf + vsaf_hBdp_al(:,0,s2nucl(is)-1)*
     &       (gqipgen/Iaf - gqipgen_al(:,s2nucl(is)-1))
          endif
          wrkf = wrkf*naf*geo%fcBb(:,2)*sqrt(geo%fcBb(:,3))/
     &      geo%fcBb(:,0)/3.0_R8

*       ..calculate the components of the parallel viscosity current
*         connected with the heat flux 
          fchvisq_a(:,0,is) = -0.4_R8*wrkf*gOnedBsq(:,1)*geo%fcS*
     &      geo%fcQalf(:,0)
          fchvisq_a(:,1,is) =  0.4_R8*wrkf*gOnedBsq(:,0)*geo%fcS*
     &      geo%fcQalf(:,1)

         endif

*      ..accumulate individual contributions into fchvisq
         fchvisq = fchvisq + fchvisq_a(:,:,is)

        enddo
       endif
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvsq

      end subroutine b2tvsq
