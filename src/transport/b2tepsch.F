      subroutine b2tepsch(nx, ny, ns, facdrift,
     &                    cddi, te, ne, cone, fhePSch)
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_geo
     &   , only : hzl,hzb
      use b2mod_constants
     &  ,only : qe
      use b2mod_geo_corner
     &  , only : bzbl,bzbv
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &  facdrift(-1:nx,-1:ny), cddi(-1:nx,-1:ny,0:1,0:ns-1),
     &  te(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  cone(-1:nx,-1:ny,0:1,0:2),                                      !srv 02.01.07
     &  fhePSch(-1:nx,-1:ny,0:1)                                        !srv 13.10.06
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ETPSCH computes the X and Y components of the Pfirsch-Schlueter
*     heat flux
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
     &  facdriftm, nete2(-1:nx,-1:ny), nete2c(-1:nx,-1:ny),
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
        call intcor(nx,ny,wght,nete2,nete2c)

*   ..compute P.Sch heat x-flux                                          !srv 13.10.06 {
        do is = 0,ns-1                                                   !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!        print *, 'style=0'
        if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
         do iy = 0, ny-1
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
            fhePSch(ix,iy,0)=0.0e0_R8
           else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            fhePSch(ix,iy,0)=-0.25e0_R8*PSch*(
     &       cddi(topix(ix,iy),topiy(ix,iy),1,is)*
     &       (ne(topix(ix,iy),topiy(ix,iy))*
     &       te(topix(ix,iy),topiy(ix,iy))**2-ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(topix(ixl,iyl),topiy(ixl,iyl),1,is)*
     &       (ne(topix(ixl,iyl),topiy(ixl,iyl))*
     &        te(topix(ixl,iyl),topiy(ixl,iyl))**2
     &                  -ne(ixl,iyl)*te(ixl,iyl)**2)+
     &       cddi(ixl,iyl,1,is)*
     &       (ne(ixl,iyl)*te(ixl,iyl)**2
     &                     -ne(bottomix(ixl,iyl),bottomiy(ixl,iyl))*
     &                      te(bottomix(ixl,iyl),bottomiy(ixl,iyl))**2)+
     &       cddi(ix,iy,1,is)*
     &       (ne(ix,iy)*te(ix,iy)**2-ne(bottomix(ix,iy),bottomiy(ix,iy))
     &                          *te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            cone(ix,iy,0,1) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*           ! 0.25*2  !srv 02.01.07
     *                                             ne(ix,iy)*te(ix,iy)
            cone(ix,iy,0,2) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*           ! 0.25*2  !srv 02.01.07
     &                      ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                      te(bottomix(ix,iy),bottomiy(ix,iy))
           endif
          enddo
          !print *,fhePSch(:,iy,0)
         enddo
        else
!          print *,'style=2'
          do iy = -1, ny
             do ix = -1, nx
               if(topiy(ix,iy).eq.ny+1 .or. bottomiy(ix,iy).eq.-2 
     &           .or. leftix(ix,iy).eq.-2 ) then
                fhePSch(ix,iy,0)= 0.0e0_R8
                cone(ix,iy,0,:) = 0.0e0_R8
              else
                fhePSch(ix,iy,0) = -PSch*
     &          (nete2c(topix(ix,iy),topiy(ix,iy))-nete2c(ix,iy))*
     &          bzbl(ix,iy)*hzl(ix,iy)/qe
                cone(ix,iy,0,1) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*       ! 0.25*2  !srv 02.01.07
     *                                         ne(ix,iy)*te(ix,iy)
                cone(ix,iy,0,2) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*       ! 0.25*2  !srv 02.01.07
     &                      ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                      te(bottomix(ix,iy),bottomiy(ix,iy))

      !           print *, iy,ix,topiy(ix,iy),topix(ix,iy)
      !           print *, nete2c(topix(ix,iy),topiy(ix,iy)),
      !&          nete2c(ix,iy) ,te(ix,iy)
              endif
            enddo
            !print *, fhePSch(:,iy,0)
          enddo
        endif

*   ..compute P.Sch y-flux
         iy = -1
        if((PSch_style.eq.0).or.(PSch_style.eq.1)) then
         fhePSch(:,iy,1)=0.0e0_R8
         do iy = 0, ny
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhePSch(ix,iy,1) = 0.0e0_R8
            cone(ix,iy,1,:) = 0.0e0_R8
           else
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            fhePSch(ix,iy,1) = 0.25e0_R8*PSch*(
     &       cddi(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &       (ne(rightix(ix,iy),rightiy(ix,iy))*
     &       te(rightix(ix,iy),rightiy(ix,iy))**2
     &                                         -ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(ix,iy,0,is)*
     &       (ne(ix,iy)*te(ix,iy)**2
     &                   -ne(leftix(ix,iy),leftiy(ix,iy))*
     &                    te(leftix(ix,iy),leftiy(ix,iy))**2)+
     &       cddi(bottomix(ix,iy),bottomiy(ix,iy),0,is)*
     &       (ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &       te(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                         -ne(leftix(ixb,iyb),leftiy(ixb,iyb))*
     &                          te(leftix(ixb,iyb),leftiy(ixb,iyb))**2)+
     &       cddi(rightix(ixb,iyb),rightiy(ixb,iyb),0,is)*
     &       (ne(rightix(ixb,iyb),rightiy(ixb,iyb))*
     &                       te(rightix(ixb,iyb),rightiy(ixb,iyb))**2
     &                      -ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                       te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            cone(ix,iy,1,1) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*ne(ix,iy)* ! 0.25*2 !srv 02.01.07 {
     &                                                te(ix,iy)
            cone(ix,iy,1,2) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*           ! 0.25*2
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
             fhePSch(ix,iy,1)= 0.0e0_R8
             cone(ix,iy,1,:) = 0.0e0_R8
            else
              fhePSch(ix,iy,1) = PSch*                                  
     &        (nete2c(rightix(ix,iy),rightiy(ix,iy))-nete2c(ix,iy))*
     &        bzbv(ix,iy)*hzb(ix,iy)/qe
              cone(ix,iy,1,1) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*         ! 0.25*2 !srv 02.01.07 {
     &                                           ne(ix,iy)*te(ix,iy)
              cone(ix,iy,1,2) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*         ! 0.25*2
     &                        ne(leftix(ix,iy),leftiy(ix,iy))*
     &                        te(leftix(ix,iy),leftiy(ix,iy))            !srv 02.01.07 }
            endif
           enddo
          enddo
        endif

         do iy = -1, ny
          do ix = -1, nx                                                 !srv 26.11.02 {
           fhePSch(ix,iy,:) =
     &       2.5e0_R8*fhePSch(ix,iy,:)*facdrift(ix,iy)
           cone(ix,iy,:,1:2) = cone(ix,iy,:,1:2)*facdrift(ix,iy)         !srv 02.01.07
!           cone(ix,iy,:,2) = cone(ix,iy,:,2)*facdrift(ix,iy)            !srv 02.01.07
          enddo
         enddo
        endif
       enddo

      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tepsch
      end subroutine b2tepsch
