*-----------------------------------------------------------------------
*.specification

      subroutine b2tvspa (nx, ny, ns, ismain, bb,                       !srv 20.09.06 { 16.10.17 19.02.21
     &  na, ua, ti, cvsa_cl, cvsa_drho, fac_vis, 
     &  cvsa_uAdp_albe, cvsa_uBdp_al,                                   !som 26.04.24
     &  fchvispar_a, fchvispar)                                         !srv 22.06.08 som 18.08.21
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_geo
     , , only : hy1
      use b2mod_zhfrtf
     , , only: avm_u, nal,                                              !som 17.08.21 30.04.24
     &         b2mod_zhfrtf_avm_u, b2mod_zhfrtf_prep                    !som 18.08.21 30.04.24
      use b2mod_b2mrel
     & , only : s2nucl, nucl2s, znnucl_int
      use b2mod_subsys
      use b2mod_b2zhco
     , , only: zh_vis_type
      use b2mod_rates                                                   !som 07.05.24
     & , only : rza
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                        !srv 19.02.21
      real (kind=R8) ::
     &  na(-1:nx,-1:ny,0:ns-1),
     &  ua(-1:nx,-1:ny,0:ns-1),
     &  ti(-1:nx,-1:ny),
     &  cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1),                                !srv 16.10.17
     &  cvsa_drho(-1:nx,-1:ny,0:1,0:ns-1),                              !som 18.08.21
     1  bb(-1:nx,-1:ny,0:3), fac_vis(-1:nx,-1:ny),                       !srv 22.06.08
     &  cvsa_uAdp_albe(-1:nx,-1:ny,0:1,0:nspecies-1,0:nspecies-1),      !som 26.04.24
     &  cvsa_uBdp_al(-1:nx,-1:ny,0:1,0:nspecies-1)                      !som 26.04.24
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchvispar_a(-1:nx,-1:ny,0:1,0:ns-1),
     &  fchvispar(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSPA computes the X and Y components of the parallel viscosity
*     current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, ix, iy, ib, is, inucl, iz
      integer, save :: no_current = 0
      real (kind=R8) ::
     &  vis_par
      real (kind=R8) ::
     *   dvpar(-1:nx,-1:ny),
     *   wrk0(-1:nx,-1:ny),
     *   t0
      save ncall, vis_par
      integer, save :: Zhdanov_closure = 0, Zhdanov_test = 0             !som 16.08.21
      integer, save :: Zhdanov_visu = 1                                  !som 16.08.21
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, vis_par/1.0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspa')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_vis_par', vis_par)
       call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)          !som 16.08.21
       call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)                !som 16.08.21
       call ipgeti ('b2tral_Zhdanov_visu', Zhdanov_visu)                !som 16.08.21
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of cvsa_cl
       do is = 0, ns-1
        call b2xvff (nx, ny, cvsa_cl(-1,-1,0,is), 'cvsa_cl')            !srv 16.10.17
       enddo
       call ipgeti ('b2tfhe_no_current', no_current)
      endif

      if (Zhdanov_closure.ne.0 .and. Zhdanov_test.ne.1 .and.
     &    Zhdanov_visu.ne.0) then
*   ..re-calculate mass-averaged velocity
        call b2mod_zhfrtf_avm_u(nx, ny, ns, na, ua)                     !som 16.08.21
*   ..re-calculate mass-averaged velocity
        call b2mod_zhfrtf_prep(nx,ny,ns,na,ua)                          !som 30.04.24
      endif

c Drift terms
      fchvispar =0.0e0_R8                                               !srv 09.07.99
      fchvispar_a=0.0e0_R8
      dvpar=0.0e0_R8                                                    !srv 24.06.99
      if (vis_par.ne.0.0e0_R8.and.maxval(fac_vis).ne.0.0_R8 .and.
     &   no_current.eq.0) then
       do is = 0, ns-1                                                  !srv 24.06.99
        if (Zhdanov_closure.eq.0 .or. Zhdanov_test.eq.1 .or.
     &      Zhdanov_visu.eq.0) then                                     !som 16.08.21
!          if (nint(zn(is)).eq.1 .and. nint(zamax(is)).eq.1 .and.
!     &        .not.is_neutral(is)) then                                !srv 24.06.99
         if(is.eq.ismain .and. .not.is_neutral(is)) then                !srv 19.02.21
          do iy =-1, ny                                                 !srv 18.02.00
           do ix = 0, nx
            if(leftix(ix,iy).ne.-2) then                                !srv 26.11.02 23.07.05 {
!           it has been considered that cvsa_cl was multiplied by bx**2
             dvpar(ix,iy)=                                              !srv 24.06.99
     &        (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &         bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &         bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8
     &         *cvsa_cl(ix,iy,0,is)*
     &         (ua(ix,iy,is)*sqrt(bb(ix,iy,3))                          !srv 16.10.17
     &         -ua(leftix(ix,iy),leftiy(ix,iy),is)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))                 !srv 18.02.00 23.07.05 }
!              dvpar(ix,iy)=                                            !srv 24.06.99
!     &            bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)            !srv
!     &            *cvsa_cl(ix,iy,0,is)*                                !srv 18.02.00
!     &        (ua(rightix(ix,iy),rightiy(ix,iy),is)                    !srv
!     &         *sqrt((bb(rightix(ix,iy),rightiy(ix,iy),3)              !srv 18.02.00
!     &             +bb(ix,iy,3))*0.5e0_R8)                             !srv 18.02.00
!     &         -ua(ix,iy,is)*sqrt((bb(leftix(ix,iy),leftiy(ix,iy),3)
!     &             +bb(ix,iy,3))*0.5e0_R8))                            !srv 18.02.00
            endif                                                       !srv 26.11.02
           enddo                                                        !srv
          enddo                                                         !srv
         else
          dvpar=0.0_R8
         endif
        else
         if (.not.is_neutral(is)) then                                  !som 09.11.21
          if (zh_vis_type.eq.0) then                                    !som 13.01.22
           do iy = -1, ny
            do ix = 0, nx
             if(leftix(ix,iy).ne.-2) then
!           it has been considered that cvsa_cl was multiplied by bx**2
              dvpar(ix,iy) =
     &        (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &         bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &              bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     &         cvsa_drho(ix,iy,0,is)*
     &        (ua(ix,iy,is)*sqrt(bb(ix,iy,3))
     &        -ua(leftix(ix,iy),leftiy(ix,iy),is)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))
             endif
            enddo
           enddo
          else if (zh_vis_type.eq.1) then
           do iy = -1, ny
            do ix = 0, nx
             if(leftix(ix,iy).ne.-2) then
!           it has been considered that cvsa_cl was multiplied by bx**2
              dvpar(ix,iy) =
     &        (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &         bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &              bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     &         cvsa_drho(ix,iy,0,is)*
     &        (avm_u(ix,iy)*sqrt(bb(ix,iy,3))
     &        -avm_u(leftix(ix,iy),leftiy(ix,iy))*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))
             endif
            enddo
           enddo
          else if (zh_vis_type.eq.2) then
           do iy = -1, ny
            do ix = 0, nx
             if(leftix(ix,iy).ne.-2) then
!           it has been considered that cvsa_cl was multiplied by bx**2
              dvpar(ix,iy) =
     &        (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &         bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &              bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     &         cvsa_cl(ix,iy,0,is)*
     &        (ua(ix,iy,is)*sqrt(bb(ix,iy,3))
     &        -ua(leftix(ix,iy),leftiy(ix,iy),is)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))
     
              do inucl=0,nspecies-1                                     !som 30.04.24
               if(inucl.ne.s2nucl(is)-1) then
                t0=0.0e0_R8
                do iz= 1, znnucl_int(inucl)
                 t0 = t0 + 0.5e0_R8*(na(ix,iy,nucl2s(inucl,iz))*
     &           ti(ix,iy)+
     &           na(leftix(ix,iy),leftiy(ix,iy),nucl2s(inucl,iz))*
     &           ti(leftix(ix,iy),leftiy(ix,iy)))*
     &           (ua(ix,iy,nucl2s(inucl,iz))*sqrt(bb(ix,iy,3))
     &           -ua(leftix(ix,iy),leftiy(ix,iy),nucl2s(inucl,iz))*
     &           sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))
                enddo
                
                dvpar(ix,iy) = dvpar(ix,iy) +
     &          (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &           bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &           sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &                bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     &           (na(ix,iy,is)/nal(ix,iy,s2nucl(is)-1)+
     &           na(leftix(ix,iy),leftiy(ix,iy),is)/
     &           nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(is)-1))*
     &           0.5e0_R8*cvsa_uAdp_albe(ix,iy,0,s2nucl(is)-1,inucl)*
     &           t0            
               endif
              enddo
              
              t0=0.0e0_R8                                               !som 07.05.24
              do iz= 1, znnucl_int(s2nucl(is)-1)
               if(iz.ne.nint(rza(ix,iy,is))) then
                t0 = t0 + 0.5e0_R8*(na(ix,iy,nucl2s(s2nucl(is)-1,iz))*
     &          ti(ix,iy)+
     &          na(leftix(ix,iy),leftiy(ix,iy),
     &          nucl2s(s2nucl(is)-1,iz))*
     &          ti(leftix(ix,iy),leftiy(ix,iy)))*
     &          (ua(ix,iy,nucl2s(s2nucl(is)-1,iz))*sqrt(bb(ix,iy,3)) -
     &          ua(leftix(ix,iy),leftiy(ix,iy),
     &          nucl2s(s2nucl(is)-1,iz))*
     &          sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))
               endif
              enddo
              
              dvpar(ix,iy) = dvpar(ix,iy) +
     &        (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     &         bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &         sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     &              bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     &         (na(ix,iy,is)/nal(ix,iy,s2nucl(is)-1)+
     &         na(leftix(ix,iy),leftiy(ix,iy),is)/
     &         nal(leftix(ix,iy),leftiy(ix,iy),s2nucl(is)-1))*
     &         0.5e0_R8*
     &         (cvsa_uAdp_albe(ix,iy,0,s2nucl(is)-1,s2nucl(is)-1)
     &         -cvsa_uBdp_al(ix,iy,0,s2nucl(is)-1))*
     &         t0
         
             endif
            enddo
           enddo
          endif
         else
          dvpar=0.0_R8
         endif
        endif

*       X-COMPONENT OF VISCOSITY PARALLEL CURRENT                       !srv
         do iy = 0, ny-1                                                !srv
          do ix = 0, nx                                                 !srv 01.07.99
           if(leftix(ix,iy).ne.-2) then                                 !srv 26.11.02
            fchvispar_a(ix,iy,0,is) =                                   !srv 16.10.17
     &        -vis_par*fac_vis(ix,iy)*                        !signmf*  !srv 25.04.00 05.02.05 22.06.08 11.01.13
     &         2.0e0_R8*dvpar(ix,iy)                                    !srv 23.07.05
     &        /(hy1(leftix(ix,iy),leftiy(ix,iy))+hy1(ix,iy))            !srv 18.02.00
     &        *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv
     &        /3.0e0_R8                                                 !srv 18.02.00
!     &        (dvpar(ix,iy)+dvpar(leftix(ix,iy),leftiy(ix,iy)))        !srv
!     &        /(hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy))             !srv 18.02.00
!     &        *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))   !srv
!     &        /3.0e0_R8                                                !srv 18.02.00
           endif                                                        !srv 26.11.02
          enddo
         enddo
*        Y-COMPONENT OF VISCOSITY PARALLEL CURRENT                      !srv
         wrk0 = 0.0e0_R8                                                !srv 23.07.05 {
         do iy = -1, ny
          do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then
             wrk0(ix,iy) = (dvpar(ix,iy)+
     &         dvpar(rightix(ix,iy),rightiy(ix,iy)))*0.5e0_R8           !srv 19.07.06
           endif
          enddo
         enddo                                                          !srv 23.07.05 }
         do iy = 0, ny                                                  !srv 01.07.99
          do ix = 0, nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
            fchvispar_a(ix,iy,1,is)=                                    !srv 16.10.17
     &        vis_par*fac_vis(ix,iy)*                       ! signmf*   !srv 25.04.00 21.04.00 22.06.08 11.01.13
     &       (wrk0(ix,iy)+wrk0(bottomix(ix,iy),bottomiy(ix,iy)))        !srv 23.07.05
     &       /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))         !srv 18.02.00
     &       *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv
     &       /3.0e0_R8                                                  !srv 18.02.00
           endif                                                        !srv 26.11.02
          enddo
         enddo
         fchvispar(:,:,:)=fchvispar(:,:,:)+fchvispar_a(:,:,:,is)
       enddo

      endif                                                             !srv

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchvispar_a(ix,iy,0,:)=0.0e0_R8
            fchvispar(ix,iy,0)=0.0e0_R8
          end do
        end do
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchvispar_a(ix,iy,1,:)=0.0e0_R8
            fchvispar(ix,iy,1)=0.0e0_R8
          end do
        end do
      end do
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspa

      end subroutine b2tvspa                                            !srv 20.09.06 }
