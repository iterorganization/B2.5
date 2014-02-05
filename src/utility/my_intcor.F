      subroutine my_intcor1(nx,ny,nz,vol,centre,corner)
      use b2mod_types
      implicit none
      integer nx,ny,nz
      real (kind=R8) ::
     * vol(-1:nx,-1:ny),centre(-1:nx,-1:ny,0:nz),
     1 corner(-1:nx+1,-1:ny+1,0:nz)
      integer is
      do is=0,nz
        call my_intcor(nx,ny,vol,centre(-1,-1,is),corner(-1,-1,is))
      enddo
      return
      end
c
      subroutine my_intcor(nx,ny,vol,centre,corner)
      use b2mod_types
      use b2mod_indirect
      implicit none

      integer nx,ny
      real (kind=R8) ::
     * vol(-1:nx,-1:ny),centre(-1:nx,-1:ny),corner(-1:nx+1,-1:ny+1)

      integer ix,iy,ix1,ix2,iy1,iy2

cxpb  This routine is coded to interpolate a function (such as B) which is
cxpb  DISCONTINUOUS across the cut in the case of diverted island geometries

      corner=0

      do iy=-1,ny
       do ix=-1,nx
        if(bottomiy(ix,iy).eq.-2) then
         if(leftix(ix,iy).eq.-2) then
          corner(ix,iy)=centre(ix,iy)                                    !xpb
         elseif(periodic_bc.ne.1 .or. leftix(ix,iy).eq.ix-1 .or.
     &         nncut.eq.0) then                                          !xpb
          corner(ix,iy)=(
     1        centre(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     2        centre(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))/
     3        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         else                                                            !xpb
          corner(ix,iy)=centre(ix,iy)                                    !xpb
      write(*,*) 'my_intcor: special treatment for cell ix,iy = ',ix,iy
         endif
        endif                                                            !xpb
       enddo
       if(bottomiy(nx,iy).eq.-2) corner(nx+1,iy)=centre(nx,iy)
      enddo

      do iy=-1,ny
        do ix=-1,nx
          if(bottomiy(ix,iy).ne.-2) then
            if(leftix(ix,iy).eq.-2) then
              corner(ix,iy)=(
     1         centre(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+
     2         centre(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     3         (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
            elseif((periodic_bc.ne.1 .or. leftix(ix,iy).eq.ix-1          !xpb
     &         .or. region(ix,iy,0).eq.5 .or. nncut.eq.0) .and.
     &        leftix(ix,iy).ne.-2 .and.
     &        bottomiy(leftix(ix,iy),leftiy(ix,iy)).ne.-2)  then         !xpb
              corner(ix,iy)=(
     1           centre(ix,iy)*
     2           vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     3               bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     4           centre(leftix(ix,iy),leftiy(ix,iy))*
     5           vol(bottomix(ix,iy),bottomiy(ix,iy))+
     6           centre(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     7                  bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     8           vol(ix,iy)+
     9           centre(bottomix(ix,iy),bottomiy(ix,iy))*
     1           vol(leftix(ix,iy),leftiy(ix,iy)))/
     1           (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     3                bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     4           vol(bottomix(ix,iy),bottomiy(ix,iy))+
     5           vol(ix,iy)+
     6           vol(leftix(ix,iy),leftiy(ix,iy)))
            elseif((periodic_bc.ne.1 .or. leftix(ix,iy).eq.ix-1          !xpb
     &         .or. region(ix,iy,0).eq.5 .or. nncut.eq.0) .and.
     &        leftix(ix,iy).ne.-2 .and.
     &        leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.-2)  then       !xpb
              corner(ix,iy)=(
     1           centre(ix,iy)*
     2           vol(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     3               leftiy(bottomix(ix,iy),bottomiy(ix,iy)))+
     4           centre(leftix(ix,iy),leftiy(ix,iy))*
     5           vol(bottomix(ix,iy),bottomiy(ix,iy))+
     6           centre(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     7                  leftiy(bottomix(ix,iy),bottomiy(ix,iy)))*
     8           vol(ix,iy)+
     9           centre(bottomix(ix,iy),bottomiy(ix,iy))*
     1           vol(leftix(ix,iy),leftiy(ix,iy)))/
     1           (vol(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     3                leftiy(bottomix(ix,iy),bottomiy(ix,iy)))+
     4           vol(bottomix(ix,iy),bottomiy(ix,iy))+
     5           vol(ix,iy)+
     6           vol(leftix(ix,iy),leftiy(ix,iy)))
            else                                                         !xpb
              corner(ix,iy)=(                                            !xpb
     1         centre(ix,iy)*vol(bottomix(ix,iy),bottomiy(ix,iy))+       !xpb
     2         centre(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/      !xpb
     4         (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))         !xpb
      write(*,*) 'my_intcor: special treatment for cell ix,iy = ',ix,iy
            endif
          endif                                                          !xpb
        enddo
        if(bottomiy(nx,iy).ne.-2) corner(nx+1,iy)=(
     1        centre(nx,iy)*vol(bottomix(nx,iy),bottomiy(nx,iy))+
     2        centre(bottomix(nx,iy),bottomiy(nx,iy))*vol(nx,iy))/
     4        (vol(bottomix(nx,iy),bottomiy(nx,iy))+vol(nx,iy))
      enddo

      do iy=-1,ny
       do ix=-1,nx
         if(leftix(ix,iy).eq.-2 .and. topiy(ix,iy).eq.ny+1) then
           corner(ix,iy+1)=centre(ix,iy)
         elseif (topiy(ix,iy).eq.ny+1) then
           corner(ix,iy+1)=(
     1        centre(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     3        centre(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))/
     5        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         endif
       enddo
       if(topiy(nx,iy).eq.ny+1) corner(nx+1,iy+1)=centre(nx,iy)
      enddo

c special treatment of standard X-points
      do iy = 0, ny-1                                                    !xpb
        do ix = 0, nx-1                                                  !xpb
          if(leftix(ix,iy).ne.-2) then                                    !srv 28.07.10
           if(leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.                !xpb !srv 28.07.10
     .       bottomix(leftix(ix,iy),leftiy(ix,iy))) then                 !xpb
            ix1=ix                                                       !xpb
            ix2=rightix(topix(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     1                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))),  !xpb
     2                  topiy(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     3                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))))  !xpb
cxpb Check that this is a standard X-point                               !xpb
            if(ix2.ge.0 .and. ix2.le.nx-1 .and. periodic_bc.ne.1) then   !xpb
              corner(ix1,iy)=                                            !xpb
     1          (centre(ix1,iy)*vol(ix2,iy)+                             !xpb
     2           centre(bottomix(ix1,iy),bottomiy(ix1,iy))*              !xpb
     2           vol(bottomix(ix2,iy),bottomiy(ix2,iy))+                 !xpb
     3           centre(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),       !xpb
     3                  leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))*      !xpb
     3           vol(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),          !xpb
     3               leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))+         !xpb
     4           centre(leftix(ix1,iy),leftiy(ix1,iy))*                  !xpb
     4           vol(leftix(ix2,iy),leftiy(ix2,iy))+                     !xpb
     5           centre(ix2,iy)*vol(ix1,iy)+                             !xpb
     6           centre(bottomix(ix2,iy),bottomiy(ix2,iy))*              !xpb
     6           vol(bottomix(ix1,iy),bottomiy(ix2,iy))+                 !xpb
     7           centre(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),       !xpb
     7                  leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))*      !xpb
     7           vol(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),          !xpb
     7               leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))+         !xpb
     8           centre(leftix(ix2,iy),leftiy(ix2,iy))*                  !xpb
     8           vol(leftix(ix1,iy),leftiy(ix1,iy)))                     !xpb
              corner(ix1,iy)=corner(ix1,iy)/                             !xpb
     1          (vol(ix2,iy)+vol(bottomix(ix2,iy),bottomiy(ix2,iy))+     !xpb
     2           vol(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),          !xpb
     3               leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))+         !xpb
     4           vol(leftix(ix2,iy),leftiy(ix2,iy))+                     !xpb
     5           vol(ix1,iy)+vol(bottomix(ix1,iy),bottomiy(ix1,iy))+     !xpb
     6           vol(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),          !xpb
     7               leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))+         !xpb
     8           vol(leftix(ix1,iy),leftiy(ix1,iy)))                     !xpb
              corner(ix2,iy)=corner(ix1,iy)                              !xpb
            elseif(ix2.ge.0 .and. ix2.le.nx-1 .and. periodic_bc.eq.1)    !xpb
     >       then                                                        !xpb
cxpb  Must treat the two sides separately when B is not a periodic variable
      write(*,*) 'my_intcor: special treatment for cell ix,iy = ',ix,iy
              corner(ix1,iy)=                                            !xpb
     1         (centre(bottomix(leftix(ix1,iy),leftiy(ix1,iy)),          !xpb
     1                 bottomiy(leftix(ix1,iy),leftiy(ix1,iy)))*         !xpb
     1          vol(bottomix(ix1,iy),bottomiy(ix1,iy))+                  !xpb
     2          centre(leftix(ix1,iy),leftiy(ix1,iy))*vol(ix1,iy)+       !xpb
     3          centre(ix1,iy)*vol(leftix(ix1,iy),leftiy(ix1,iy))+       !xpb
     4          centre(bottomix(ix1,iy),bottomiy(ix1,iy))*               !xpb
     4          vol(bottomix(leftix(ix1,iy),leftiy(ix1,iy)),             !xpb
     4              bottomiy(leftix(ix1,iy),leftiy(ix1,iy))))/           !xpb
     1         (vol(bottomix(ix1,iy),bottomiy(ix1,iy))+                  !xpb
     2          vol(ix1,iy)+vol(leftix(ix1,iy),leftiy(ix1,iy))+          !xpb
     4          vol(bottomix(leftix(ix1,iy),leftiy(ix1,iy)),             !xpb
     4              bottomiy(leftix(ix1,iy),leftiy(ix1,iy))))            !xpb
              corner(ix2,iy)=                                            !xpb
     1         (centre(bottomix(leftix(ix2,iy),leftiy(ix2,iy)),          !xpb
     1                 bottomiy(leftix(ix2,iy),leftiy(ix2,iy)))*         !xpb
     1          vol(bottomix(ix2,iy),bottomiy(ix2,iy))+                  !xpb
     2          centre(leftix(ix2,iy),leftiy(ix2,iy))*vol(ix2,iy)+       !xpb
     3          centre(ix2,iy)*vol(leftix(ix2,iy),leftiy(ix2,iy))+       !xpb
     4          centre(bottomix(ix2,iy),bottomiy(ix2,iy))*               !xpb
     4          vol(bottomix(leftix(ix2,iy),leftiy(ix2,iy)),             !xpb
     4              bottomiy(leftix(ix2,iy),leftiy(ix2,iy))))/           !xpb
     1         (vol(bottomix(ix2,iy),bottomiy(ix2,iy))+                  !xpb
     2          vol(ix2,iy)+vol(leftix(ix2,iy),leftiy(ix2,iy))+          !xpb
     4          vol(bottomix(leftix(ix2,iy),leftiy(ix2,iy)),             !xpb
     4              bottomiy(leftix(ix2,iy),leftiy(ix2,iy))))            !xpb
            endif                                                        !xpb
           endif                                                         !srv 28.07.10
		endif                                                          !xpb
        enddo                                                            !xpb
      enddo                                                              !xpb

c special treatment of diverted island X-points                          !xpb
      do iy = -1, ny-1                                                   !xpb
        do ix = -1, nx-1                                                 !xpb
          if(topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1) then  !xpb
            ix1=topix(rightix(ix,iy),rightiy(ix,iy))                     !xpb
            iy1=topiy(rightix(ix,iy),rightiy(ix,iy))                     !xpb
            ix2=topix(ix,iy)                                             !xpb
            corner(ix1,iy1)=(centre(ix1,iy1)*vol(ix2,iy)+                !xpb
     2        centre(bottomix(ix1,iy1),bottomiy(ix1,iy1))*               !xpb
     2        vol(leftix(ix2,iy),leftiy(ix2,iy))+                        !xpb
     3        centre(leftix(ix1,iy),leftiy(ix1,iy))*                     !xpb
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+                      !xpb
     4        centre(ix2,iy)*vol(ix1,iy1)+                               !xpb
     5        centre(leftix(ix2,iy),leftiy(ix2,iy))*                     !xpb
     5        vol(bottomix(ix1,iy1),bottomiy(ix1,iy1))+                  !xpb
     6        centre(leftix(ix1,iy1),leftiy(ix1,iy1))*                   !xpb
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))/                       !xpb
     1       (vol(ix2,iy)+vol(leftix(ix2,iy),leftiy(ix2,iy))+            !xpb
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+vol(ix1,iy1)+         !xpb
     5        vol(bottomix(ix1,iy1),bottomiy(ix1,iy1))+                  !xpb
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))                        !xpb
cxpb No corner(ix2,iy1) is assigned as this would map into a dead cell   !xpb 
          endif                                                          !xpb
        enddo                                                            !xpb
      enddo                                                              !xpb

c special treatment of limiter tangency X-points
      ix = -1
      do iy = 0, ny
        if(bottomiy(ix,iy).eq.-2) then
          ix1=rightix(ix,iy)
          iy1=rightiy(ix,iy)
          ix2=rightix(
     1     topix(leftix(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     2                  bottomiy(rightix(ix,iy),rightiy(ix,iy))),
     3           leftiy(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     4                  bottomiy(rightix(ix,iy),rightiy(ix,iy)))),
     1     topiy(leftix(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     2                  bottomiy(rightix(ix,iy),rightiy(ix,iy))),
     3            leftiy(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     4                   bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
          iy2=rightiy(
     1     topix(leftix(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     2                  bottomiy(rightix(ix,iy),rightiy(ix,iy))),
     3           leftiy(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     4                  bottomiy(rightix(ix,iy),rightiy(ix,iy)))),
     1     topiy(leftix(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     2                  bottomiy(rightix(ix,iy),rightiy(ix,iy))),
     3            leftiy(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     4                   bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
          corner(ix1,iy1)=(centre(ix1,iy1)*vol(ix2,iy)+
     2        centre(topix(ix1,iy1),topiy(ix1,iy1))*
     2        vol(leftix(ix2,iy),leftiy(ix2,iy))+
     3        centre(leftix(ix1,iy),leftiy(ix1,iy))*
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+
     4        centre(ix2,iy)*vol(ix1,iy1)+
     5        centre(leftix(ix2,iy),leftiy(ix2,iy))*
     5        vol(topix(ix1,iy1),topiy(ix1,iy1))+
     6        centre(leftix(ix1,iy1),leftiy(ix1,iy1))*
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))/
     1       (vol(ix2,iy)+vol(leftix(ix2,iy),leftiy(ix2,iy))+
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+vol(ix1,iy1)+
     5        vol(topix(ix1,iy1),topiy(ix1,iy1))+
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))
          corner(ix2,iy2)=corner(ix1,iy1)
        endif
      enddo

      return
      end
