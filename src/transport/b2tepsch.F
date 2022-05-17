      subroutine b2tepsch(nx, ny, ns, hx, hy, qs, facdrift, cdde, te,    !srv 19.02.21 !som 03.11.21
     & ne, cone, fhePSch)
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_geo
     &   , only : hzl,hzb
      use b2mod_constants
     &  ,only : qe
      use b2mod_geo_corner
     &  , only : bzbl,bzbv
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &  facdrift(-1:nx,-1:ny), cdde(-1:nx,-1:ny,0:1,0:1),                !som 03.11.21
     &  te(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qs(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  cone(-1:nx,-1:ny,0:1,0:1,0:2),                                   !srv 02.01.07
     &  fhePSch(-1:nx,-1:ny,0:1,0:1)                                     !srv 13.10.06
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ETPSCH computes the X and Y components of the Pfirsch-Schlueter
*     heat flux.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, is, ixl, iyl, ixb, iyb
      integer, save :: ncall = 0
      integer, save :: PSch_style
      real (kind=R8), save :: PSch=1.0e0_R8                              !srv 17.01.07
      real (kind=R8) ::
     &  facdriftm, nete2(-1:nx,-1:ny), dnete2(-1:nx,-1:ny,0:1,0:1),
     &  wght(-1:nx,-1:ny)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tepsch')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
       if (redef_gmtry.eq.0) then
         PSch_style = 1
       else
         PSch_style = 2
      end if
       call ipgeti ('b2tfnb_PSch_style', PSch_style)
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

      fhePSch = 0.0e0_R8                                                 !srv 13.10.06
      cone    = 0.0e0_R8                                                 !srv 02.01.07
      facdriftm = maxval(facdrift)
      if(facdriftm.ne.0.0e0_R8) then

*   ..precompute differences
        wght = 1.0_R8
        nete2=ne*te*te
        call diff(nx, ny, hx, hy, wght, qs, nete2, dnete2)

*   ..compute P.Sch heat x-flux                                          !srv 13.10.06 {
!       do is = 0, ns-1                                                  !srv 12.05.11
!        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then             !srv 12.05.11
!       if(is.eq.ismain) then                                            !srv 19.02.21
!        print *, 'style=0'
        if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
         do iy = 0, ny-1
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
            fhePSch(ix,iy,0,0:1)=0.0e0_R8
           else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            fhePSch(ix,iy,0,0)=-0.25e0_R8*PSch*(
     &       cdde(topix(ix,iy),topiy(ix,iy),1,1)*                        !som 03.11.21
     &       dnete2(topix(ix,iy),topiy(ix,iy),1,1)+
     &       cdde(topix(ixl,iyl),topiy(ixl,iyl),1,1)*                    !som 03.11.21
     &       dnete2(topix(ixl,iyl),topiy(ixl,iyl),1,1)+
     &       cdde(ixl,iyl,1,1)*                                          !som 03.11.21
     &       dnete2(ixl,iyl,1,1)+
     &       cdde(ix,iy,1,1)*                                            !som 03.11.21
     &       dnete2(ix,iy,1,1)
     &                               )
cwdk computation of radial PSch flux on x-faces
cwdk somewhat messy interpolation of geometric factors, for exact
cwdk cancellation of some terms in poloidal and radial flux components
            fhePSch(ix,iy,0,1)=-0.25e0_R8*PSch*(
     &       cdde(topix(ix,iy),topiy(ix,iy),1,1)*
     &       dnete2(topix(ix,iy),topiy(ix,iy),1,0)*
     &       qs(topix(ix,iy),topiy(ix,iy),1)*
     &      (hy(ix,iy) + hy(topix(ix,iy),topiy(ix,iy)))/
     &      (hx(ix,iy) + hx(topix(ix,iy),topiy(ix,iy)))+
     &       cdde(topix(ixl,iyl),topiy(ixl,iyl),1,1)*
     &       dnete2(topix(ixl,iyl),topiy(ixl,iyl),1,0)*
     &       qs(topix(ixl,iyl),topiy(ixl,iyl),1)*
     &      (hy(ixl,iyl) + hy(topix(ixl,iyl),topiy(ixl,iyl)))/
     &      (hx(ixl,iyl) + hx(topix(ixl,iyl),topiy(ixl,iyl)))+
     &       cdde(ixl,iyl,1,1)*
     &       dnete2(ixl,iyl,1,0)*qs(ixl,iyl,1)*
     &      (hy(ixl,iyl) + hy(bottomix(ixl,iyl),bottomiy(ixl,iyl)))/
     &      (hx(ixl,iyl) + hx(bottomix(ixl,iyl),bottomiy(ixl,iyl)))+
     &       cdde(ix,iy,1,1)*
     &       dnete2(ix,iy,1,0)*qs(ix,iy,1)*
     &      (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &      (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))
     &                               )
            cone(ix,iy,0,0,1) = -0.5e0_R8*PSch*cdde(ix,iy,1,1)*          ! 0.25*2  !srv 02.01.07
     *                                             ne(ix,iy)*te(ix,iy)
            cone(ix,iy,0,1,1) = -cone(ix,iy,0,0,1)
            cone(ix,iy,0,0,2) = -0.5e0_R8*PSch*cdde(ix,iy,1,1)*          ! 0.25*2  !srv 02.01.07
     &                      ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                      te(bottomix(ix,iy),bottomiy(ix,iy))
            cone(ix,iy,0,1,2) = -cone(ix,iy,0,0,2)
           endif
          enddo
          !print *,fhePSch(:,iy,0,0)
         enddo
        else
!          print *,'style=2'
          do iy = -1, ny
             do ix = -1, nx
               if(topiy(ix,iy).eq.ny+1 .or. bottomiy(ix,iy).eq.-2
     &           .or. leftix(ix,iy).eq.-2 ) then
                fhePSch(ix,iy,0,0:1)= 0.0e0_R8
                cone(ix,iy,0,:,:) = 0.0e0_R8
              else
                fhePSch(ix,iy,0,0)=-PSch*
     &          dnete2(ix,iy,1,1)*
     &          bzbl(ix,iy)*hzl(ix,iy)/qe
                cone(ix,iy,0,0,1) = -0.5e0_R8*PSch*cdde(ix,iy,1,1)*      ! 0.25*2  !srv 02.01.07
     *                                         ne(ix,iy)*te(ix,iy)
                cone(ix,iy,0,1,1) = -cone(ix,iy,0,0,1)
                cone(ix,iy,0,0,2) = -0.5e0_R8*PSch*cdde(ix,iy,1,1)*      ! 0.25*2  !srv 02.01.07
     &                      ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                      te(bottomix(ix,iy),bottomiy(ix,iy))
                cone(ix,iy,0,1,2) = -cone(ix,iy,0,0,2)

      !           print *, iy,ix,topiy(ix,iy),topix(ix,iy)
      !           print *, nete2c(topix(ix,iy),topiy(ix,iy)),
      !&          nete2c(ix,iy) ,te(ix,iy)
              endif
            enddo
            !print *, fhePSch(:,iy,0,0)
          enddo
        endif

*   ..compute P.Sch y-flux
         iy = -1
        if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
         fhePSch(:,iy,1,0:1)=0.0e0_R8
         do iy = 0, ny
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhePSch(ix,iy,1,0:1) = 0.0e0_R8
            cone(ix,iy,1,0:1,:) = 0.0e0_R8
           else
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            fhePSch(ix,iy,1,1) = 0.25e0_R8*PSch*(
     &       cdde(rightix(ix,iy),rightiy(ix,iy),0,0)*                    !som 03.11.21
     &       dnete2(rightix(ix,iy),rightiy(ix,iy),0,0)+
     &       cdde(ix,iy,0,0)*                                            !som 03.11.21
     &       dnete2(ix,iy,0,0)+
     &       cdde(bottomix(ix,iy),bottomiy(ix,iy),0,0)*                  !som 03.11.21
     &       dnete2(bottomix(ix,iy),bottomiy(ix,iy),0,0)+
     &       cdde(rightix(ixb,iyb),rightiy(ixb,iyb),0,0)*                !som 03.11.21
     &       dnete2(rightix(ixb,iyb),rightiy(ixb,iyb),0,0)
     &                               )
            cone(ix,iy,1,1,1) = -0.5e0_R8*PSch*cdde(ix,iy,0,0)*          ! 0.25*2 !srv 02.01.07 {
     &                                      ne(ix,iy)*te(ix,iy)
            cone(ix,iy,1,1,2) = -0.5e0_R8*PSch*cdde(ix,iy,0,0)*          ! 0.25*2
     &                        ne(leftix(ix,iy),leftiy(ix,iy))*
     &                        te(leftix(ix,iy),leftiy(ix,iy))            !srv 02.01.07 }
           endif
          enddo
         enddo
        else
          do iy = -1, ny
           do ix = -1, nx
            if(bottomiy(ix,iy).eq.-2
     &       .or. leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
             fhePSch(ix,iy,1,0:1)= 0.0e0_R8
             cone(ix,iy,1,:,:) = 0.0e0_R8
            else
              fhePSch(ix,iy,1,1) = PSch*
     &        dnete2(ix,iy,0,0)*
     &        bzbv(ix,iy)*hzb(ix,iy)/qe
              cone(ix,iy,1,1,1) = -0.5e0_R8*PSch*cdde(ix,iy,0,0)*     ! 0.25*2 !som 03.11.21 !srv 02.01.07 {
     &                                           ne(ix,iy)*te(ix,iy)
              cone(ix,iy,1,1,2) = -0.5e0_R8*PSch*cdde(ix,iy,0,0)*     ! 0.25*2 !som 03.11.21
     &                        ne(leftix(ix,iy),leftiy(ix,iy))*
     &                        te(leftix(ix,iy),leftiy(ix,iy))            !srv 02.01.07 }
            endif
           enddo
          enddo
        endif

         do iy = -1, ny
          do ix = -1, nx                                                 !srv 26.11.02 {
           fhePSch(ix,iy,:,:) =
     &       2.5e0_R8*fhePSch(ix,iy,:,:)*facdrift(ix,iy)
           cone(ix,iy,:,:,1:2) = cone(ix,iy,:,:,1:2)*facdrift(ix,iy)     !srv 02.01.07
!           cone(ix,iy,:,2) = cone(ix,iy,:,2)*facdrift(ix,iy)            !srv 02.01.07
          enddo
         enddo
!       endif
!      enddo

      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tepsch
      end subroutine b2tepsch
