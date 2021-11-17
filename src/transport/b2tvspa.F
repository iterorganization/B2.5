*-----------------------------------------------------------------------
*.specification

      subroutine b2tvspa (nx, ny, ns, ismain, bb, na, ua, cvsa_cl,      !srv 20.09.06 { 16.10.17 19.02.21
     & cvsa_drho, fac_vis, fchvispar_a, fchvispar)                       !srv 22.06.08 som 18.08.21
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_geo
     , , only : hy1
      use b2mod_zhfrtf ,
     & only: avm_u, rho_A_rel,                                          !som 17.08.21
     & b2mod_zhfrtf_avm_u                                               !som 18.08.21
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain                                        !srv 19.02.21
      real (kind=R8) ::
     & na(-1:nx,-1:ny,0:ns-1),
     & ua(-1:nx,-1:ny,0:ns-1),
     & cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 16.10.17
     & cvsa_drho(-1:nx,-1:ny,0:1,0:ns-1),                               !som 18.08.21
     1 bb(-1:nx,-1:ny,0:3), fac_vis(-1:nx,-1:ny)                        !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchvispar_a(-1:nx,-1:ny,0:1,0:ns-1), fchvispar(-1:nx,-1:ny,0:1)
*   ..common blocks
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
      integer ncall, ix, iy, ib, is
      integer, save :: no_current = 0
      real (kind=R8) ::
     &  vis_par
      real (kind=R8) ::
     *   dvpar(-1:nx,-1:ny),
     *   wrk0(-1:nx,-1:ny)
      save ncall, vis_par
      integer, save :: Zhdanov_closure = 0, Zhdanov_test = 0           !som 16.08.21
      integer, save :: Zhdanov_visu = 1                                !som 16.08.21
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
       call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)        !som 16.08.21
       call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)              !som 16.08.21
       call ipgeti ('b2tral_Zhdanov_visu', Zhdanov_visu)              !som 16.08.21
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of cvsa_cl
        call b2xvff (nx, ny, cvsa_cl, 'cvsa_cl')                        !srv 16.10.17
        call ipgeti ('b2tfhe_no_current', no_current)
      endif
      
      if ((Zhdanov_closure.ne.0) .and. (Zhdanov_test.ne.1).and.
     & Zhdanov_visu.ne.0) then
*   ..re-calculate mass-average velocity  
        call b2mod_zhfrtf_avm_u(nx, ny, ns, na, ua)                     !som 16.08.21    
      endif

c Drift terms
      fchvispar =0.0e0_R8                                               !srv 09.07.99
      fchvispar_a=0.0e0_R8
      dvpar=0.0e0_R8                                                    !srv 24.06.99
      if (vis_par.ne.0.0e0_R8.and.maxval(fac_vis).ne.0.0_R8 .and.   
     &   no_current.eq.0) then
       do is = 0, ns-1                                                  !srv 24.06.99
        if ((Zhdanov_closure.eq.0) .or. (Zhdanov_test.eq.1).or.
     & Zhdanov_visu.eq.0) then                                          !som 16.08.21
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
        else                                                            !som 16.08.21      
         if (.not.is_neutral(is)) then                                  !som 09.11.21                                  
          do iy = -1, ny                                           
           do ix = 0, nx
            if(leftix(ix,iy).ne.-2) then                           
!           it has been considered that cvsa_cl was multiplied by bx**2
             dvpar(ix,iy)=                           
     & (bb(ix,iy,2)*sqrt(bb(ix,iy,3))/bb(ix,iy,0)+
     & bb(leftix(ix,iy),leftiy(ix,iy),2)*
     & sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))/
     & bb(leftix(ix,iy),leftiy(ix,iy),0))*0.5e0_R8*
     & cvsa_drho(ix,iy,0,is)*
     & (avm_u(ix,iy)*sqrt(bb(ix,iy,3))                    
     & -avm_u(leftix(ix,iy),leftiy(ix,iy))*
     & sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3)))         
            endif                                                   
           enddo                                                     
          enddo 
         else
          dvpar=0.0_R8
         endif                                                                                                                
        endif
        
*       X-COMPONENT OF VISCOSITY PARALLEL CURRENT                       !srv
         do iy = 0, ny-1                                                !srv
          do ix = 0, nx                                                 !srv 01.07.99
           if(leftix(ix,iy).ne.-2) then                                 !srv 26.11.02
            fchvispar_a(ix,iy,0,is)=                                    !srv 16.10.17
     &        -vis_par*fac_vis(ix,iy)*                        !signmf*  !srv 25.04.00 05.02.05 22.06.08 11.01.13
     &         2.0e0_R8*dvpar(ix,iy)                                    !srv 23.07.05
     &        /(hy1(leftix(ix,iy),leftiy(ix,iy))+hy1(ix,iy))            !srv 18.02.00
     &        *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv
     &        /3.0e0_R8                                                 !srv 18.02.00
!     &        (dvpar(ix,iy)+dvpar(leftix(ix,iy),leftiy(ix,iy)))         !srv
!     &        /(hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy))              !srv 18.02.00
!     &        *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv
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
         fchvispar(:,:,:) = fchvispar(:,:,:) + fchvispar_a(:,:,:,is)
       enddo

      endif                                                             !srv
      
! it will be soon removed {
       do is = 0, ns-1
        if (is.eq.ismain) then
         fchvispar_a(:,:,:,is) = fchvispar(:,:,:) 
        else
         fchvispar_a(:,:,:,is) = 0.0_R8
        endif
       enddo
! }      

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
