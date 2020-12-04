*-----------------------------------------------------------------------
*.specification
      subroutine b2tipsch (nx, ny, ns, facdrift, !vol, hx, hy,       !srv 02.04.18 {
     1  na, ti, cddi, coni, fhiPSch)                                    !srv 13.10.06
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_geo_corner
     &  , only : bzbl
      use b2mod_geo
!     &   , only : hzl,hzb
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  facdrift(-1:nx,-1:ny), ! vol(-1:nx,-1:ny),
!     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),          !srv 13.10.06
     1  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny),
     &  cddi(-1:nx,-1:ny,0:1,0:ns-1),
     &  nati2c(-1:nx,-1:ny,0:ns-1),nati2(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fhiPSch(-1:nx,-1:ny,0:1), coni(-1:nx,-1:ny,0:1,0:2),           !srv 13.10.06
     &   weight(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TIPSCH computes the Pfirsch-Schlueter all atom heat flux between
*     neighboring cells, and the coefficients that specify its linearization.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the Pfirsch-Schlueter atom heat flux, fhi,
*     and the coefficients coni that linearize the response of the
*     atom heat flux to changes in the atom temperature.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib, ixl, iyl, ixb, iyb
      integer, save :: iout = 0
      real (kind=R8), save :: PSch=1.0e0_R8                             !srv 17.01.07
      integer, save :: PSch_style
      real (kind=R8) ::
     *   facdriftm                                                      !srv 25.01.02
      save ncall
*   ..procedures
      external xertst
      external b2xvsg, b2xvff
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tipsch')
*   ..set internal parameters on first call                             !srv 17.05.99
      if (ncall.eq.0) then                                              !srv
       call ipgeti ('b2tipsch_iout', iout)                              !srv 16.09.02
       call ipgetr ('b2tfhi_PSch', PSch)                                !srv 17.01.07
       if (redef_gmtry.eq.0) then
         PSch_style = 1
       else
         PSch_style = 2
       end if
       call ipgeti ('b2tfnb_PSch_style', PSch_style)
      endif                                                             !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ti
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

      facdriftm = maxval(facdrift)
      fhiPSch = 0.0e0_R8                                                !srv 13.10.06
      coni    = 0.0e0_R8                                                !srv 02.01.07
      if(facdriftm.ne.0.0e0_R8) then
*   ..compute P.Sch heat x-flux                                         !srv 13.10.06 {
!      iy = -1
!       do ix = -1, nx
!         if(leftix(ix,iy).eq.-2) then
!         fhiPSch(ix,iy,0)=0.0e0_R8
!         else
!          do is = 0, ns-1
!           if(za(is) .ne. 0.0e0_R8) then
!          fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.5e0_R8*(
!     &    cddi(ix,topiy(ix,iy),1,is)*
!     &    (na(ix,topiy(ix,iy),is)*ti(ix,topiy(ix,iy))**2-
!     &                        na(ix,iy,is)*ti(ix,iy)**2)+
!     &    cddi(leftix(ix,iy),topiy(ix,iy),1,is)*
!     &     (na(leftix(ix,iy),topiy(ix,iy),is)*
!     &      ti(leftix(ix,iy),topiy(ix,iy))**2
!     &                 -na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2)
!     &                               )
!           endif
!          enddo
!         endif
!        enddo
       if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           fhiPSch(ix,iy,0) = 0.0e0_R8
           coni(ix,iy,0,:) = 0.0e0_R8
          else
           do is = 0, ns-1
            if(.not.is_neutral(is)) then
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.25e0_R8*PSch*(
     &      cddi(topix(ix,iy),topiy(ix,iy),1,is)*
     &      (na(topix(ix,iy),topiy(ix,iy),is)*
     &      ti(topix(ix,iy),topiy(ix,iy))**2-na(ix,iy,is)*ti(ix,iy)**2)+
     &      cddi(topix(ixl,iyl),topiy(ixl,iyl),1,is)*
     &      (na(topix(ixl,iyl),topiy(ixl,iyl),is)*
     &      ti(topix(ixl,iyl),topiy(ixl,iyl))**2
     &                -na(ixl,iyl,is)*ti(ixl,iyl)**2)+
     &      cddi(ixl,iyl,1,is)*
     &      (na(ixl,iyl,is)*ti(ixl,iyl)**2
     &                     -na(bottomix(ixl,iyl),bottomiy(ixl,iyl),is)*
     &                      ti(bottomix(ixl,iyl),bottomiy(ixl,iyl))**2)+
     &      cddi(ix,iy,1,is)*
     &      (na(ix,iy,is)*ti(ix,iy)**2
     &             -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &              ti(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            coni(ix,iy,0,1) = coni(ix,iy,0,1) + 0.5e0_R8*PSch*  ! 0.25*2  !srv 02.01.07
     &               cddi(ix,iy,1,is)*na(ix,iy,is)*ti(ix,iy)
            coni(ix,iy,0,2) = coni(ix,iy,0,2) + 0.5e0_R8*PSch*  ! 0.25*2  !srv 02.01.07
     &      cddi(ix,iy,1,is)*na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
            endif
           enddo
          endif
         enddo
        enddo
       else
!        print *,'b2tipsch style=2'
        weight = 1.0e0_R8
        nati2 = 0.0
        do is = 0, ns-1
          if(.not.is_neutral(is)) then
!            print *, is
            nati2(:,:,is)=na(:,:,is)*ti(:,:)*ti(:,:)
          endif
        enddo
        call intcor1(nx,ny,ns-1,weight,nati2,nati2c)
!        print *, 'intcore done'
        do iy = -1, ny
         do ix = -1, nx
          if(topiy(ix,iy).eq.ny+1 .or. bottomiy(ix,iy).eq.-2 
     &       .or. leftix(ix,iy).eq.-2) then
            fhiPSch(ix,iy,0) = 0.0e0_R8
            coni(ix,iy,0,:) = 0.0e0_R8
          else
            do is = 0, ns-1
              if(.not.is_neutral(is)) then
                fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 
     &          PSch*(nati2c(topix(ix,iy),topiy(ix,iy),is)-
     &       nati2c(ix,iy,is))/rza(ix,iy,is)/qe*bzbl(ix,iy)*hzl(ix,iy)
                coni(ix,iy,0,1) = coni(ix,iy,0,1) + 0.5e0_R8*PSch*  ! 0.25*2  !srv 02.01.07
     &               cddi(ix,iy,1,is)*na(ix,iy,is)*ti(ix,iy)
                coni(ix,iy,0,2) = coni(ix,iy,0,2) + 0.5e0_R8*PSch*  ! 0.25*2  !srv 02.01.07
     &      cddi(ix,iy,1,is)*na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
              endif
            enddo
          endif
          enddo
        enddo
!        print *, 'psch x done'
       endif

       do is = 0, ns-1
         if(.not.is_neutral(is)) then
         endif
       enddo
!      iy = ny
!       do ix = -1, nx
!         if(leftix(ix,iy).eq.-2) then
!         fhiPSch(ix,iy,0)=0.0e0_R8
!         else
!          do is = 0, ns-1
!          fhiPSch(ix,iy,0) = fhiPSch(ix,iy,0) + 0.5e0_R8*(
!     &     cddi(leftix(ix,iy),iy,1,is)*
!     &    (na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2
!     &                         -na(leftix(ix,iy),bottomiy(ix,iy),is)*
!     &                            ti(leftix(ix,iy),bottomiy(ix,iy))**2)+
!     &     cddi(ix,iy,1,is)*
!     &    (na(ix,iy,is)*ti(ix,iy)**2
!     &             -na(ix,bottomiy(ix,iy),is)*ti(ix,bottomiy(ix,iy))**2)
!     &                               )
!          enddo
!         endif
!        enddo
       if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
*   ..compute P.Sch y-flux
        iy = -1
        fhiPSch(:,iy,1)=0.0e0_R8
        do iy = 0, ny
         do ix = -1, nx
          do is = 0, ns-1
  !          if(leftix(ix,iy).eq.-2) then
  !          fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.5e0_R8*(
  !     &     cddi(rightix(ix,iy),iy,0,is)*
  !     &     (na(rightix(ix,iy),iy,is)*ti(rightix(ix,iy),iy)**2
  !     &                                     -na(ix,iy,is)*ti(ix,iy)**2)+
  !     &     cddi(rightix(ix,iy),bottomiy(ix,iy),0,is)*
  !     &     (na(rightix(ix,iy),bottomiy(ix,iy),is)*
  !     &             ti(rightix(ix,iy),bottomiy(ix,iy))**2
  !     &           -na(ix,bottomiy(ix,iy),is)*ti(ix,bottomiy(ix,iy))**2)
  !     &                               )
  !          elseif(rightix(ix,iy).eq.nx+1) then
  !          fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.5e0_R8*(
  !     &     cddi(ix,iy,0,is)*(na(ix,iy,is)*ti(ix,iy)**2
  !     &                -na(leftix(ix,iy),iy,is)*ti(leftix(ix,iy),iy)**2)+
  !     &     cddi(ix,bottomiy(ix,iy),0,is)*
  !     &     (na(ix,bottomiy(ix,iy),is)*ti((ix),bottomiy(ix,iy))**2
  !     &                         -na(leftix(ix,iy),bottomiy(ix,iy),is)*
  !     &                            ti(leftix(ix,iy),bottomiy(ix,iy))**2)
  !     &                               )
  !          else
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhiPSch(ix,iy,1) = 0.0e0_R8
            coni(ix,iy,1,:)  = 0.0e0_R8
           else
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 0.25e0_R8*PSch*(
     &      cddi(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &      (na(rightix(ix,iy),rightiy(ix,iy),is)*
     &       ti(rightix(ix,iy),rightiy(ix,iy))**2
     &                                     -na(ix,iy,is)*ti(ix,iy)**2)+
     &      cddi(ix,iy,0,is)*(na(ix,iy,is)*ti(ix,iy)**2
     &                -na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                 ti(leftix(ix,iy),leftiy(ix,iy))**2)+
     &      cddi(bottomix(ix,iy),bottomiy(ix,iy),0,is)*
     &      (na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &       ti(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                         -na(leftix(ixb,iyb),leftiy(ixb,iyb),is)*
     &                          ti(leftix(ixb,iyb),leftiy(ixb,iyb))**2)+
     &      cddi(rightix(ixb,iyb),rightiy(ixb,iyb),0,is)*
     &      (na(rightix(ixb,iyb),rightiy(ixb,iyb),is)*
     &       ti(rightix(ixb,iyb),rightiy(ixb,iyb))**2
     &             -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &              ti(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            coni(ix,iy,1,1) = coni(ix,iy,1,1)+0.5e0_R8*PSch*     ! 0.25*2 !srv 02.01.07 {
     &                      cddi(ix,iy,0,is)*na(ix,iy,is)*ti(ix,iy)
            coni(ix,iy,1,2) = coni(ix,iy,1,2)+0.5e0_R8*PSch*     ! 0.25*2
     &       cddi(ix,iy,0,is)*na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                        ti(leftix(ix,iy),leftiy(ix,iy))           !srv 02.01.07 }
           endif
          enddo
         enddo
        enddo
       else
!        print *,'psch y start'
        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2 
     &       .or. leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhiPSch(ix,iy,1) = 0.0e0_R8
            coni(ix,iy,1,:)  = 0.0e0_R8
          else
            do is = 0, ns-1
              if(.not.is_neutral(is)) then
                fhiPSch(ix,iy,1) = fhiPSch(ix,iy,1) - 
     &          PSch*(nati2c(rightix(ix,iy),rightiy(ix,iy),is)-
     &       nati2c(ix,iy,is))/rza(ix,iy,is)/qe*bzbl(ix,iy)*hzl(ix,iy)
                coni(ix,iy,1,1) = coni(ix,iy,1,1)+0.5e0_R8*PSch*     ! 0.25*2 !srv 02.01.07 {
     &                      cddi(ix,iy,0,is)*na(ix,iy,is)*ti(ix,iy)
                coni(ix,iy,1,2) = coni(ix,iy,1,2)+0.5e0_R8*PSch*     ! 0.25*2
     &       cddi(ix,iy,0,is)*na(leftix(ix,iy),leftiy(ix,iy),is)*
     &                        ti(leftix(ix,iy),leftiy(ix,iy))           !srv 02.01.07 }
              endif
            enddo
          endif
          enddo
        enddo
!        print *, 'psch y done'
       endif

       do iy = -1, ny
        do ix = 0, nx
         fhiPSch(ix,iy,:) = 2.5e0_R8*fhiPSch(ix,iy,:)*facdrift(ix,iy)
         coni(ix,iy,:,1:2) = coni(ix,iy,:,1:2)*facdrift(ix,iy)          !srv 02.01.07
        enddo
       enddo

c  ..set null flux in guard cells  ?
!      fhiPSch(-1, :,:)=0.0e0_R8    ! MODIFY for Double Case
!      fhiPSch(nx, :,1)=0.0e0_R8
!      fhiPSch( :,-1,:)=0.0e0_R8
!      fhiPSch( :,ny,0)=0.0e0_R8
      endif

csw 18feb2014 remove fluxes between guard cells
      do ix=-1,nx
        do iy=-1,ny
          if(iffboundx(ix,iy).ne.0) then
            fhiPSch(ix,iy,0)=0.0_R8
            coni(ix,iy,0,:)=0.0_R8
          endif
          if(iffboundy(ix,iy).ne.0) then
            fhiPSch(ix,iy,1)=0.0_R8
            coni(ix,iy,1,:)=0.0_R8
          endif
        enddo
      enddo

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
       DO IX=NXISO1(IB), NXISO2(IB)
        DO IY=NYISO1(IB), NYISO2(IB)
         fhipSch(ix,iy,0) = 0.0e0_R8
         coni(ix,iy,0,:)  = 0.0e0_R8                                     !srv 02.01.07
        enddo
       end do
       DO IX=NXISO1(IB), NXISO2(IB)
        DO IY=NYISO1(IB), NYISO2(IB)
          fhipSch(ix,iy,1) = 0.0e0_R8
          coni(ix,iy,1,:)  = 0.0e0_R8                                    !srv 02.01.07
        enddo
       enddo
      enddo
      if(iout.ne.0) then                                                !srv 18.06.08 {
       call my_out(70,nx,ny,coni(-1,-1,0,0),'b2tipsch_conix0')
       call my_out(70,nx,ny,coni(-1,-1,0,1),'b2tipsch_conix1')
       call my_out(70,nx,ny,coni(-1,-1,0,2),'b2tipsch_conix2')
       call my_out(70,nx,ny,coni(-1,-1,1,0),'b2tipsch_coniy0')
       call my_out(70,nx,ny,coni(-1,-1,1,1),'b2tipsch_coniy1')
       call my_out(70,nx,ny,coni(-1,-1,1,2),'b2tipsch_coniy2')
       call my_out(70,nx,ny,fhiPSch(-1,-1,0),'b2tipsch_fhiPSchx')
       call my_out(70,nx,ny,fhiPSch(-1,-1,1),'b2tipsch_fhiPSchy')
      endif                                                             !srv 18.06.08 }
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tipsch

      end subroutine b2tipsch                                           !srv 02.04.18 }
