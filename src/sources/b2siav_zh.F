*-----------------------------------------------------------------------
*.specification

      subroutine b2siav_zh (nCv, nFc, nVx, ns, isb, switch, geo, mpg,
     & mb, rzb, te, ti, na, ua, ni, ne2, rz2, lnlam, hci_a, hci_al_ast,
     & cvsbhz_cl, cvsbhz_drho, vsbf_drho, vsaf_hAdp_albe, vsaf_hBdp_al,
     & vsaf_uAdp_albe, vsaf_uBdp_al,
     & smbvh, smbvv, smbvi, rho_rel, avm_u,
     & wrk0, wrk1, wrk2)                                         !srv 23.07.99 04.10.99 17.01.01 09.01.01 31.07.02 01.07.05
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_balance                                                 !djm Jan2017
     & , only : b2siav_smovh0to3, b2siav_smovv0to3, balance_netcdf
      use b2mod_ad
     & , only : ncall_b2siav
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2ttia
#endif
      use b2mod_zhfrtf
     & , only : nal, Z_to_m1_ast, Ia
      use b2mod_b2mrel
     & , only : amnucl, s2nucl, nucl2s, znnucl_int
      use b2mod_b2zhco
     & , only : zhcscorr, zhcscorr_vq
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     & mb, lnlam(nCv),                                                  !srv 01.07.09
     & na(nCv,0:ns-1), ua(nCv,0:ns-1),
     & te(nCv), ti(nCv), rzb(nCv), ni(nCv), ne2(nCv),
     & hci_a(nCv,0:ns-1),
     & hci_al_ast(nCv,0:nspecies-1),
     & cvsbhz_cl(nFc,0:1),
     & cvsbhz_drho(nFc,0:1),
     & vsbf_drho(nFc,0:1),
     & rz2(nCv,0:ns-1),                                                 !srv 31.07.02
     & vsaf_hAdp_albe(nFc,0:1,0:nspecies-1,0:nspecies-1),               !som 08.10.21
     & vsaf_hBdp_al(nFc,0:1,0:nspecies-1),                              !som 08.10.21
     & vsaf_uAdp_albe(nFc,0:1,0:nspecies-1,0:nspecies-1),               !som 08.10.21
     & vsaf_uBdp_al(nFc,0:1,0:nspecies-1),                              !som 08.10.21
     & rho_rel(nCv), avm_u(nCv)                                         !som 13.08.21
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbvh(nCv,0:3), smbvv(nCv,0:3), smbvi(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(nCv), wrk1(nCv), wrk2(nCv), wrkns(nFc,0:nspecies-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAV_ZH computes a linearised representation of the additional
*     viscosity if the Zhdanov closure is applied.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ic, iCv, iFc, iFt, k, inucl, iz
      integer, save :: nsmooth = 2
      character chns*3, chk*1                                           !srv 17.06.02
      real (kind=R8) ::
     *     t0, ft, ft_h                                                 !srv 16.10.17
      real (kind=R8) ::
     &     wrks(nCv), wrkf(nFc), wrkv(nVx), wght(nFc,2),
     *     wrkf0(nFc), wrkf1(nFc), wrkf2(nFc), wrkf3(nFc),
     *     tauia(nCv,0:ns-1), qip0(nCv),                                !srv 23.07.99
     *     qip1(nCv), qipgen(nCv), qipgen_al(nCv,0:nspecies-1),         !srv 23.07.99 19.01.05
     &     gti(nCv,0:1), u_res_b(nCv)
*   ..procedures
      intrinsic sqrt
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2siav')

      ft = 4.0e0_R8/3.0e0_R8                                            !srv 16.10.17
      ft_h = 8.0_R8/15.0_R8
      smbvh = 0.0_R8
      smbvv = 0.0_R8
      smbvi = 0.0_R8

      if(.not.is_neutral(isb)) then                                     !srv 26.10.08

        if (switch%addvis.ne.0.0_R8) then !srv 19.01.05 28.03.06 {

          ! init
          qip1   = 0.0e0_R8
          qipgen = 0.0e0_R8
          wrk0   = 0.0e0_R8

          call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)             !srv 01.07.09

*       ..compute radial gradients of ti in cell centers
          call gradc(nCv, nFc, nVx, 0, geo, mpg, ti, wrkv, gti)

*       ..calculate the parallel heat flux of species isb in the center of the cell
          qip0   = 0.0e0_R8

          if (switch%b2siav_qstyle.eq.0.and.minval(rzb).gt.0.0_R8) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -2.5e0_R8*geo%cvBb(iCv,3)*geo%cvBb(iCv,2)/
     &           (geo%cvBb(iCv,0)*geo%ftBbav2(iFt))*
     &           (na(iCv,isb)*ti(iCv)/(qe*rzb(iCv)))*gti(iCv,1)
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

          if (switch%b2siav_qstyle.eq.1.and.minval(rzb).gt.0.0_R8) then

            do iFt = 1, mpg%nFt 
             if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFt,1)))) then
              do ic=mpg%ftCvP(iFt,1),mpg%ftCvP(iFt,1)+mpg%ftCvP(iFt,2)-1
               iCv = mpg%ftCv(ic)
               qip0(iCv) = -hci_a(iCv,isb)*gti(iCv,0)*
     &           geo%cvBb(iCv,3)/geo%cvBb(iCv,0) - 2.5e0_R8*
     &           geo%cvBb(iCv,2)/(geo%cvBb(iCv,0)*geo%cvBb(iCv,3))*
     &           (na(iCv,isb)*ti(iCv)/(qe*rzb(iCv)))*gti(iCv,1)
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
          qip1 = switch%cqip1*na(:,isb)*ti*sqrt(te+ti)/sqrt(am(isb)*mp)

*       ..compute qipgen
          qipgen = qip1*datan(qip0/qip1)*sqrt(geo%cvBb(:,3))

*       ..calculate the parallel charge state averaged heat flux of all species
          do inucl= 0, nspecies-1

            qip0 = 0.0e0_R8

            if (switch%b2siav_qstyle.eq.0.and.
     &          minval(rzb).gt.0.0_R8) then

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
              wght = 1.0_R8
              do k = 1, nsmooth
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

            if (switch%b2siav_qstyle.eq.1.and.
     &          minval(rzb).gt.0.0_R8)  then

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

*         ..compute qipgen
            qipgen_al(:,inucl) = qip1*datan(qip0/qip1)*
     &        sqrt(geo%cvBb(:,3))

          enddo

*       ..main contribution

*       ..gradient qipgen_al in cell centers   
          !wdk modified computation compared to structured code to 
          !wdk facilitate implementation on unstructured grids
          do inucl = 0, nspecies-1
            call grad_p(nCv, nFc, nVx, 0, geo, mpg, qipgen_al(:,inucl),
     &        wrkv, wrkns(:,inucl))
          enddo

          wrkf0 = 0.0_R8
          do inucl = 0, nspecies-1
            wrkf0 = wrkf0 + wrkns(:,inucl)*
     &        vsaf_hAdp_albe(:,0,s2nucl(isb)-1,inucl)
          enddo

          wrk0 = na(:,isb)/nal(:,s2nucl(isb)-1)
          call intface(nCv,nFc,mpg%fcCv,wght,wrk0,wrkf1)

          wrkf0 = wrkf0*wrkf1/(geo%fcBb(:,0)*geo%fcBb(:,3))

*       ..finally calculate the r.h.s. contribution
          call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf0, wrk0)
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk1)

          !zero out spurious contributions outside of separatrix
          !wdk todo: rationalize this when separatrix-structure is introduced?
          do iCv = 1, nCv
            if (.not.mpg%cvOnClosedSurface(iCv)) then
              wrk1(iCv) = 0.0_R8
            endif
          enddo

          !calculate the actual source
          smbvh(:,0) = switch%addvis*ft_h*geo%cvVol*geo%cvBb(:,0)*
     &      sqrt(geo%cvBb(:,3))*geo%cvHz*wrk1

*       ..zero out in guard cells
          smbvh(mpg%nCi+1:nCv,0) = 0.0_R8

*       ..charge state correction
          if (zhcscorr.gt.0 .and. zhcscorr_vq.gt.0) then

*         ..gradient qipgen in cell centers   
            call grad_p(nCv, nFc, nVx, 0, geo, mpg, qipgen, wrkv, wrkf0)

*         ..calculate the difference from parallel heat flux of ion
            wrk0 = na(:,isb)/nal(:,s2nucl(isb)-1)
            call intface(nCv,nFc,mpg%fcCv,wght,wrk0,wrkf1)

            wrkf0 = wrkf1/(geo%fcBb(:,0)*geo%fcBb(:,3))*
     &        vsaf_hBdp_al(:,0,s2nucl(isb)-1)*
     &        (wrkf0/Ia(:,isb) - wrkns(:,s2nucl(isb)-1))

*         ..calculate the r.h.s. contribution to the parallel balance
*           corresponding to parallel viscosity connected with heat flow
            call intcell(nFc, nCv, mpg, mpg%intcellP, wrkf0, wrk0)
            call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk1)

            !zero out spurious contributions outside of separatrix
            !wdk todo: rationalize this when separatrix-structure is introduced?
            do iCv = 1, nCv
              if (.not.mpg%cvOnClosedSurface(iCv)) then
                wrk1(iCv) = 0.0_R8
              endif
            enddo

            !calculate the actual source
            smbvh(:,0) = smbvh(:,0) + 
     &        switch%addvis*ft_h*geo%cvVol*geo%cvBb(:,0)*
     &        sqrt(geo%cvBb(:,3))*geo%cvHz*wrk1

*         ..zero out in guard cells
            smbvh(mpg%nCi+1:nCv,0) = 0.0_R8

          endif

        endif                                                           !lk 23.05.06 05.07.07

        if (switch%addvis1.ne.0.0_R8) then                              !srv 19.01.05 30.01.10

*       ..prepare wrk arrays
          wrk2 = sqrt(geo%cvBb(:,3))                                    ! wrk2 is used below twice
          wrk1 = log(geo%cvHz*wrk2)

*       ..poloidal differences of wrk1 on cell faces
          wrkv = log(geo%vxHz*sqrt(geo%vxBb(:,3))+1.0e-30_R8)
          call diff_p(nCv, nFc, nVx, 1, geo, mpg, wrk1, wrkv, wrkf)

*       ..compute flow and its divergence
          if (switch%zh_vis_type.eq.2) then
            wrkf = wrkf*cvsbhz_cl(:,0)/sqrt(geo%fcBb(:,3))
            call div0(nCv, nFc, mpg, wrkf, wrk0)
          else
            wrkf = wrkf*cvsbhz_drho(:,0)/sqrt(geo%fcBb(:,3))
            call div0(nCv, nFc, mpg, wrkf, wrk0)
          endif

*       ..add to r.h.s. source
          if (switch%zh_vis_type.eq.0.or.switch%zh_vis_type.eq.2) then
            wrk1 = ft*wrk0*wrk2
          elseif (switch%zh_vis_type.eq.1) then
            wrk1 = ft*wrk0*wrk2*rho_rel
          else
            call xerrab ('wrong value of zh_vis_type')
          endif
          do iCv = 1, mpg%nCi
            t0 = wrk1(iCv)
            if (t0.gt.0.0e0_R8) then
              smbvv(iCv,0) = t0*ua(iCv,isb)
            else
              smbvv(iCv,1) = t0
            endif
          enddo                                                         !srv 19.01.05 }


*       ..viscosity source with contribution from other ion velocities  !som 13.08.21

*       ..calculate part of the mass-averaged velocity which depends on the other
*         ion velocities
          if (switch%zh_vis_type.eq.1) then

            u_res_b = (avm_u - rho_rel*ua(:,isb))*wrk2

*         ..gradient u_res in cell centers   
            call grad_p(nCv, nFc, nVx, 0, geo, mpg, u_res_b, wrkv,
     &        wrkf0)

*         ..calculate the difference from parallel heat flux of ion
            wrkf0 = wrkf0*vsbf_drho(:,0)/(geo%fcBb(:,0)*geo%fcBb(:,3))

*         ..calculate the r.h.s. contribution to the parallel balance
*           corresponding to parallel viscosity connected with heat flow
            call intcell(nFc,nCv,mpg, mpg%intcellP, wrkf0,wrk0)
            call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk1)

            !calculate the actual source
            smbvi(:,0) = ft*wrk1*wrk2*geo%cvVol*geo%cvBb(:,0)*geo%cvHz

*         ..zero out in guard cells
            smbvi(mpg%nCi+1:nCv,0) = 0.0_R8

          elseif (switch%zh_vis_type.eq.2) then

            wrkf = 0.0_R8
            do inucl=0,nspecies-1

              if(inucl.ne.s2nucl(isb)-1) then

                wrkf0=0.0_R8
                do iz= 1, znnucl_int(inucl)

                  u_res_b = ua(:,nucl2s(inucl,iz))*wrk2
                  call grad_p(nCv, nFc, nVx, 0, geo, mpg, u_res_b,
     &              wrkv, wrkf1)

                  wrk1 = na(:,nucl2s(inucl,iz))*ti
                  call intface(nCv,nFc,mpg%fcCv,wght,wrk1,wrkf2)

                  wrkf0 = wrkf0 + wrkf1*wrkf2

                enddo

                wrkf = wrkf +
     &            wrkf0*vsaf_uAdp_albe(:,0,s2nucl(isb)-1,inucl)

              endif
            enddo

            wrkf0=0.0_R8
            do iz = 1, znnucl_int(s2nucl(isb)-1)
              !wdk not clear if Zh-closure applies to bundling;
              !wdk if yes: need to update if-condition here
              if(iz.ne.nint(maxval(rzb)).or.
     &           iz.ne.nint(minval(rzb))) then

                u_res_b = ua(:,nucl2s(s2nucl(isb)-1,iz))*wrk2
                call grad_p(nCv, nFc, nVx, 0, geo, mpg, u_res_b,
     &            wrkv, wrkf1)

                wrk1 = na(:,nucl2s(s2nucl(isb)-1,iz))*ti
                call intface(nCv,nFc,mpg%fcCv,wght,wrk1,wrkf2)

                wrkf0 = wrkf0 + wrkf1*wrkf2

              endif
            enddo

            wrkf = wrkf + wrkf0*
     &        (vsaf_uAdp_albe(:,0,s2nucl(isb)-1,s2nucl(isb)-1) - 
     &         zhcscorr*vsaf_uBdp_al(:,0,s2nucl(isb)-1))

            !wdk common multipliers
            wrk1 = na(:,isb)/nal(:,s2nucl(isb)-1)
            call intface(nCv,nFc,mpg%fcCv,wght,wrk1,wrkf1)

            wrkf = wrkf*wrkf1/(geo%fcBb(:,0)*geo%fcBb(:,3))

*         ..calculate the r.h.s. contribution to the parallel balance
*           corresponding to parallel viscosity connected with heat flow
            call intcell(nFc,nCv,mpg, mpg%intcellP, wrkf,wrk0)
            call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk1)

            !calculate the actual source
            smbvi(:,0) = ft*wrk1*wrk2*geo%cvVol*geo%cvBb(:,0)*geo%cvHz

*         ..zero out in guard cells
            smbvi(mpg%nCi+1:nCv,0) = 0.0_R8

          endif

        endif

      endif                                                             !srv 26.10.08

      if (switch%b2siav_iout.ne.0) then                                 !srv 17.06.02 {
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,smbvv(1,k),'b2siav_smovv'//chk//chns)
          call my_out_us(70,nCv,0,smbvh(1,k),'b2siav_smovh'//chk//chns)
        enddo
        wrks=smbvv(:,0)+smbvv(:,1)*ua(:,isb)+
     &       smbvv(:,2)*na(:,isb)*mb+smbvv(:,3)*ua(:,isb)*na(:,isb)*mb
        call my_out_us(70,nCv,0,wrks,'b2siav_smovv'//chns)
        wrks=smbvh(:,0)+smbvh(:,1)*ua(:,isb)+
     &       smbvh(:,2)*na(:,isb)*mb+smbvh(:,3)*ua(:,isb)*na(:,isb)*mb
        call my_out_us(70,nCv,0,wrks,'b2siav_smovh'//chns)
        call my_out_us(70,nFc,1,cvsbhz_drho(1,0),
     &       'b2siav_cvsbhz_drho_th'//chns)
        call my_out_us(70,nFc,1,cvsbhz_drho(1,1),
     &       'b2siav_cvsbhz_drho_r'//chns)
      endif                                                             !srv 17.06.02 }

cdjm Jan2017 Store linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2siav_smovh0to3(1:nCv,0:3,isb)=smbvh
        b2siav_smovv0to3(1:nCv,0:3,isb)=smbvv
      endif

* ..return
      ncall_b2siav = ncall_b2siav+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2siav_zh

      end subroutine b2siav_zh
