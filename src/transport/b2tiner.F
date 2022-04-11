*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner (nx, ny, ns, ismain, vol, gs, bb,               !srv 20.09.06 {
     &                    na, ua, facdrift, fchinert_a, fchinert)        !srv 22.06.08
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_geo
     &   , only : hyl, hxb                                               ! IYS 04.03.2019
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     & na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1 bb(-1:nx,-1:ny,0:3), facdrift(-1:nx,-1:ny)                        !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     & fchinert_a(-1:nx,-1:ny,0:1,0:ns-1), fchinert(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINER computes the X and Y components of the inertial current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, ix, iy, ib, is
      integer, save :: fch_inert_style                                   !srv 19.04.13
      integer, save :: no_current = 0
      integer, save :: inert_all_ions = 1
      integer, save :: guard_flows
      real (kind=R8) ::
     &  fch_inert
      save ncall, fch_inert
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_inert/1.0e0_R8/, guard_flows/2/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tiner')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tiner_inert', fch_inert)                          !srv 07.07.99
       call ipgeti ('b2tfhe_no_current', no_current)
       call ipgeti ('b2tiner_inert_all_ions', inert_all_ions)
       call ipgeti ('b2news_guard_flows', guard_flows)
       if (redef_gmtry.eq.0) then
          fch_inert_style = 1
       else
          fch_inert_style = 2
       end if
       call ipgeti ('b2tiner_fch_inert_style', fch_inert_style)          !srv 07.07.99
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

c Drift terms
      fchinert   = 0.0e0_R8
      fchinert_a = 0.0e0_R8
      if(fch_inert.ne.0.0e0_R8.and.maxval(facdrift).ne.0.0_R8.and.
     .   no_current.eq.0)then                                            !srv 07.07.99

        if(fch_inert_style.eq.0) then

         do is = 0, ns-1                                                 !srv 07.07.99
          if (.not.is_neutral(is)) then                                  !srv 07.07.99
*         X-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny-1                                                  !srv 07.07.99
             do ix=0,nx                                                  !srv 07.07.99
               if(leftix(ix,iy).ne.-2) then                              !srv 26.11.02
                fchinert_a(ix,iy,0,is)=                                  !srv 07.07.99
     &         +0.5e0_R8*am(is)*mp*0.5e0_R8*fch_inert*      ! signmf*    !srv 23.07.05 11.01.13
     &          facdrift(ix,iy)*                                         !srv 22.06.08
     &          (ua(ix,iy,is)**2+ua(leftix(ix,iy),leftiy(ix,iy),is)**2)  !srv 07.07.99 05.02.05 23.07.05
     &         *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8 !srv 07.07.99
     &         *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8   !srv 07.07.99
     &         *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))            !srv 07.07.99
     &         /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))          !srv 07.07.99
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
*         Y-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny                                                    !srv 07.07.99
             do ix=0,nx-1                                                !srv 07.07.99
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then !srv 26.11.02
               fchinert_a(ix,iy,1,is)=                                   !srv 07.07.99
     &          -0.5e0_R8*am(is)*mp*fch_inert*facdrift(ix,iy)  ! signmf* !srv 07.07.99 21.04.00 22.06.08 11.01.13
     &          *(                                                       !srv 07.07.99
     &           na(ix,iy,is)+na(bottomix(ix,iy),bottomiy(ix,iy),is)     !srv 07.07.99
     &           )/2.0e0_R8                                              !srv 07.07.99
     &          *0.5e0_R8*(ua(ix,iy,is)**2 +                             !srv 07.07.99 18.02.00 23.07.05
     &                     ua(bottomix(ix,iy),bottomiy(ix,iy),is)**2)    !srv 07.07.99 23.07.05
     &          *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          /2.0e0_R8                                                !srv 07.07.99
     &          *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 07.07.99
     &          /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
          endif                                                          !srv 07.07.99
          fchinert(:,:,:) = fchinert(:,:,:) + fchinert_a(:,:,:,is)
         enddo                                                           !srv 07.07.99


        elseif(fch_inert_style.eq.1) then

         do is = 0, ns-1                                                 !srv 07.07.99
          if (.not.is_neutral(is)) then                                  !srv 07.07.99
*         X-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny-1                                                  !srv 07.07.99
             do ix=0,nx                                                  !srv 07.07.99
               if(leftix(ix,iy).ne.-2) then                              !srv 26.11.02
                fchinert_a(ix,iy,0,is)=                                  !srv 07.07.99
     &         +0.5e0_R8*am(is)*mp*0.5e0_R8*fch_inert*      ! signmf*    !srv 23.07.05 11.01.13
     &          facdrift(ix,iy)*                                         !srv 22.06.08
     &          (ua(ix,iy,is)**2*na(ix,iy,is)+                           !IYS 01.02.2013
     &           ua(leftix(ix,iy),leftiy(ix,iy),is)**2*                  !IYS 01.02.2013
     &          na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
     &         *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8   !srv 07.07.99
     &         *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))            !srv 07.07.99
     &         /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))          !srv 07.07.99
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
*         Y-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny                                                    !srv 07.07.99
             do ix=0,nx-1                                                !srv 07.07.99
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then !srv 26.11.02
                fchinert_a(ix,iy,1,is)=                                  !srv 07.07.99
     &          -0.5e0_R8*am(is)*mp*fch_inert*facdrift(ix,iy) ! signmf*  !srv 07.07.99 21.04.00 22.06.08 11.01.13
     &          *0.5e0_R8*(na(ix,iy,is)*ua(ix,iy,is)**2 +                !IYS 01.02.2013
     &          na(bottomix(ix,iy),bottomiy(ix,iy),is)*                  !IYS 01.02.2013
     &          ua(bottomix(ix,iy),bottomiy(ix,iy),is)**2)               !IYS 01.02.2013 !sv 07.07.99
     &          *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          /2.0e0_R8                                                !srv 07.07.99
     &          *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 07.07.99
     &          /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
          endif                                                          !srv 07.07.99
          fchinert(:,:,:) = fchinert(:,:,:) + fchinert_a(:,:,:,is)
         enddo                                                           !srv 07.07.99

        elseif(fch_inert_style.eq.2) then

         do is = 0, ns-1                                                 !srv 07.07.99
          if (.not.is_neutral(is)) then                                  !srv 07.07.99
*         X-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny-1                                                  !srv 07.07.99
             do ix=0,nx                                                  !srv 07.07.99
               if(leftix(ix,iy).ne.-2) then                              !srv 26.11.02
                fchinert_a(ix,iy,0,is)=                                  !srv 07.07.99
     &         +0.5e0_R8*am(is)*mp*0.5e0_R8*fch_inert*      ! signmf*    !srv 23.07.05 11.01.13
     &          facdrift(ix,iy)*                                         !srv 22.06.08
     &         (ua(ix,iy,is)**2 * na(ix,iy,is) +                         ! IYS 04.03.2019
     &          ua(leftix(ix,iy),leftiy(ix,iy),is)**2 *
     &          na(leftix(ix,iy),leftiy(ix,iy),is)) *
     &          gs(ix,iy,0) * wbbl(ix,iy,2) / hyl(ix,iy)
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
*         Y-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny                                                    !srv 07.07.99
             do ix=0,nx-1                                                !srv 07.07.99
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then !srv 26.11.02
                fchinert_a(ix,iy,1,is)=                                  !srv 07.07.99
     &          -0.5e0_R8*am(is)*mp*fch_inert*facdrift(ix,iy)  ! signmf* !srv 07.07.99 21.04.00 22.06.08 11.01.13
     &          *0.5e0_R8*( na(ix,iy,is) * ua(ix,iy,is)**2 +             ! IYS 04.03.2019
     &           na(bottomix(ix,iy),bottomiy(ix,iy),is) *
     &           ua(bottomix(ix,iy),bottomiy(ix,iy),is)**2) *
     &           gs(ix,iy,1) * wbbv(ix,iy,2) / hxb(ix,iy)
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 07.07.99
               endif                                                     !srv 26.11.02
             enddo                                                       !srv 07.07.99
           enddo                                                         !srv 07.07.99
          endif                                                          !srv 07.07.99
          fchinert(:,:,:) = fchinert(:,:,:) + fchinert_a(:,:,:,is)
         enddo                                                           !srv 07.07.99
        else
          CALL XERRAB ('B2TINER -- INCORRECT VALUE OF fch_inert_style')
        endif
      endif                                                              !srv 07.07.99

      if (inert_all_ions.eq.0) then
       do is = 0, ns-1
        if (is.eq.ismain) then
         fchinert_a(:,:,:,is) = fchinert(:,:,:) 
        else
         fchinert_a(:,:,:,is) = 0.0_R8
        endif
       enddo
      endif

      if(guard_flows.eq.0) then                                          !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fchinert_a(ix,iy,1,:)=0.0_R8
              fchinert(ix,iy,1)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                fchinert_a(bc_list_x(ix,ib),bc_list_y(ix,ib),0,:) =
     .                                                    0.0_R8
                fchinert(bc_list_x(ix,ib),bc_list_y(ix,ib),0) =
     .                                                    0.0_R8
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              fchinert_a(ix,-1,0,:) = 0.0_R8
              fchinert(ix,-1,0) = 0.0_R8
            endif
          enddo
        endif
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
       DO IB=1,NNISO
         DO IX=NXISO1(IB), NXISO2(IB)
           DO IY=NYISO1(IB), NYISO2(IB)
             fchinert_a(ix,iy,0,:)=0.0e0_R8
             fchinert(ix,iy,0)=0.0e0_R8
           end do
         end do
         DO IX=NXISO1(IB), NXISO2(IB)
           DO IY=NYISO1(IB), NYISO2(IB)
             fchinert_a(ix,iy,1,:)=0.0e0_R8
             fchinert(ix,iy,1)=0.0e0_R8
           end do
         end do
       end do
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner

      end subroutine b2tiner                                            !srv 20.09.06 }
