      subroutine intcor1(nx,ny,nz,vol,centre,corner)
      use b2mod_types
      implicit none
      integer nx,ny,nz
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny,0:nz),
     1 corner(-1:nx,-1:ny,0:nz)
      integer is
      do is=0,nz
        call intcor(nx,ny,vol,centre(-1,-1,is),corner(-1,-1,is))
      enddo
      return
      end
c
      subroutine intcor(nx,ny,vol,centre,corner)
      use b2mod_types
      use b2mod_indirect
      implicit none

      integer nx,ny
      real (kind=R8) ::
     *   vol(-1:nx,-1:ny),centre(-1:nx,-1:ny),corner(-1:nx,-1:ny)
      
      integer ix,iy,ix1,ix2,iy1,iy2

      corner=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2.and.bottomiy(ix,iy).ne.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)/vol(ix,iy) + 
     2       centre(leftix(ix,iy),leftiy(ix,iy))/
     3          vol(leftix(ix,iy),leftiy(ix,iy)) + 
     4       centre(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     5              bottomiy(leftix(ix,iy),leftiy(ix,iy)))/
     6          vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     7              bottomiy(leftix(ix,iy),leftiy(ix,iy))) + 
     8       centre(bottomix(ix,iy),bottomiy(ix,iy))/
     9          vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     1       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     2              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     3        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     4        1/vol(ix,iy)+
     5        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          elseif(leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).ne.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)/vol(ix,iy) + 
     2       centre(bottomix(ix,iy),bottomiy(ix,iy))/
     3          vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     4       (1/vol(bottomix(ix,iy),bottomiy(ix,iy)) + 1/vol(ix,iy))
          elseif(leftix(ix,iy).ne.-2.and.bottomiy(ix,iy).eq.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)/vol(ix,iy) +
     2       centre(leftix(ix,iy),leftiy(ix,iy))/
     3          vol(leftix(ix,iy),leftiy(ix,iy)))/
     4       (1/vol(ix,iy)+1/vol(leftix(ix,iy),leftiy(ix,iy)))
          elseif(leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).eq.-2) then
            corner(ix,iy)=centre(ix,iy)
          endif
        enddo
      enddo

c special treatment of standard X-points
      do iy = 0, ny-1                                                    !xpb
        do ix = 0, nx-1                                                  !xpb
!          if(leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2 .and.
!     .       leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.                 !xpb
!     .       bottomix(leftix(ix,iy),leftiy(ix,iy))) then                 !xpb
          if(leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2) then       !srv 19.04.13
           if(leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.                !xpb
     .       bottomix(leftix(ix,iy),leftiy(ix,iy))) then                 !xpb
            ix1=ix                                                       !xpb
            ix2=rightix(topix(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     1                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))),  !xpb
     2                  topiy(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     3                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))))  !xpb
cxpb Check that this is a standard X-point                               !xpb
            if(ix2.ge.0 .and. ix2.le.nx) then                            !xpb
              corner(ix1,iy)=                                            !xpb
     1          (centre(ix1,iy)/vol(ix1,iy) + 
     2           centre(bottomix(ix1,iy),bottomiy(ix1,iy))/              !xpb          
     2              vol(bottomix(ix1,iy),bottomiy(ix1,iy)) +             !xpb
     3           centre(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),       !xpb
     3                  leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))/      !xpb
     3              vol(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),       !xpb
     3                  leftiy(bottomix(ix1,iy),bottomiy(ix1,iy))) +     !xpb
     4           centre(leftix(ix1,iy),leftiy(ix1,iy))/                  !xpb
     4              vol(leftix(ix1,iy),leftiy(ix1,iy)) +                 !xpb
     5           centre(ix2,iy)/vol(ix2,iy) +                            !xpb
     6           centre(bottomix(ix2,iy),bottomiy(ix2,iy))/              !xpb          
     6              vol(bottomix(ix2,iy),bottomiy(ix2,iy)) +             !xpb
     7           centre(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),       !xpb
     7                  leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))/      !xpb
     7              vol(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),       !xpb
     7                  leftiy(bottomix(ix2,iy),bottomiy(ix2,iy))) +     !xpb
     8           centre(leftix(ix2,iy),leftiy(ix2,iy))/
     8              vol(leftix(ix2,iy),leftiy(ix2,iy)))
              corner(ix1,iy)=corner(ix1,iy)/
     1          (1/vol(ix2,iy)+1/vol(bottomix(ix2,iy),bottomiy(ix2,iy))+ !xpb
     2           1/vol(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),        !xpb
     3                 leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))+       !xpb
     4           1/vol(leftix(ix2,iy),leftiy(ix2,iy))+                   !xpb
     5           1/vol(ix1,iy)+1/vol(bottomix(ix1,iy),bottomiy(ix1,iy))+ !xpb
     6           1/vol(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),        !xpb
     7                 leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))+       !xpb
     8           1/vol(leftix(ix1,iy),leftiy(ix1,iy)))                   !xpb
              corner(ix2,iy)=corner(ix1,iy)                              !xpb
            endif                                                        !xpb
          endif                                                          !xpb
	   endif
        enddo                                                            !xpb
      enddo                                                              !xpb


#ifdef WG_TEST
      corner=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2.and.bottomiy(ix,iy).ne.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)*
     2        vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     3            bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     4       centre(leftix(ix,iy),leftiy(ix,iy))*
     5        vol(bottomix(ix,iy),bottomiy(ix,iy))+
     6       centre(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     7            bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     8        vol(ix,iy)+
     9       centre(bottomix(ix,iy),bottomiy(ix,iy))*
     1        vol(leftix(ix,iy),leftiy(ix,iy)))/
     1       (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     3          bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     4        vol(bottomix(ix,iy),bottomiy(ix,iy))+
     5        vol(ix,iy)+
     6        vol(leftix(ix,iy),leftiy(ix,iy)))
          elseif(leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).ne.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)*
     5        vol(bottomix(ix,iy),bottomiy(ix,iy))+
     9       centre(bottomix(ix,iy),bottomiy(ix,iy))*vol(ix,iy))/
     1       (vol(bottomix(ix,iy),bottomiy(ix,iy))+
     5        vol(ix,iy))
          elseif(leftix(ix,iy).ne.-2.and.bottomiy(ix,iy).eq.-2) then
            corner(ix,iy)=(
     1       centre(ix,iy)*vol(leftix(ix,iy),leftiy(ix,iy))+
     4       centre(leftix(ix,iy),leftiy(ix,iy))*vol(ix,iy))/
     1       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
          elseif(leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).eq.-2) then
            corner(ix,iy)=centre(ix,iy)
          endif
        enddo
      enddo

c special treatment of standard X-points
      do iy = 0, ny-1                                                    !xpb
        do ix = 0, nx-1                                                  !xpb
!          if(leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2 .and.
!     .       leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.                 !xpb
!     .       bottomix(leftix(ix,iy),leftiy(ix,iy))) then                 !xpb
          if(leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2) then       !srv 19.04.13
           if(leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.                !xpb
     .       bottomix(leftix(ix,iy),leftiy(ix,iy))) then                 !xpb
            ix1=ix                                                       !xpb
            ix2=rightix(topix(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     1                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))),  !xpb
     2                  topiy(leftix(bottomix(ix,iy),bottomiy(ix,iy)),   !xpb
     3                        leftiy(bottomix(ix,iy),bottomiy(ix,iy))))  !xpb
cxpb Check that this is a standard X-point                               !xpb
            if(ix2.ge.0 .and. ix2.le.nx) then                            !xpb
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
              corner(ix1,iy)=corner(ix1,iy)/
     1          (vol(ix2,iy)+vol(bottomix(ix2,iy),bottomiy(ix2,iy))+     !xpb
     2           vol(leftix(bottomix(ix2,iy),bottomiy(ix2,iy)),          !xpb
     3               leftiy(bottomix(ix2,iy),bottomiy(ix2,iy)))+         !xpb
     4           vol(leftix(ix2,iy),leftiy(ix2,iy))+                     !xpb
     5           vol(ix1,iy)+vol(bottomix(ix1,iy),bottomiy(ix1,iy))+     !xpb
     6           vol(leftix(bottomix(ix1,iy),bottomiy(ix1,iy)),          !xpb
     7               leftiy(bottomix(ix1,iy),bottomiy(ix1,iy)))+         !xpb
     8           vol(leftix(ix1,iy),leftiy(ix1,iy)))                     !xpb
              corner(ix2,iy)=corner(ix1,iy)                              !xpb
            endif                                                        !xpb
          endif                                                          !xpb
	   endif
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
cxpb No corner(ix2,iy2) is assigned as this would map into a dead cell   !xpb 
          endif                                                          !xpb
        enddo                                                            !xpb
      enddo                                                              !xpb

c special treatment of limiter tangency X-points                         !xpb
      do iy = 0, ny-1                                                    !xpb
        do ix = -1, nx-1                                                 !xpb
          if(bottomiy(bottomix(ix,iy),bottomiy(ix,iy)).eq.iy 
     &     .and. ny.gt.1) then                                           !xpb
            ix1=rightix(ix,iy)                                           !xpb
            iy1=rightiy(ix,iy)                                           !xpb
            ix2=bottomix(ix,iy)                                          !xpb
            iy2=bottomiy(ix,iy)
            corner(ix1,iy1)=(centre(ix1,iy1)*vol(ix2,iy)+                !xpb
     2        centre(topix(ix1,iy1),topiy(ix1,iy1))*                     !xpb
     2        vol(leftix(ix2,iy),leftiy(ix2,iy))+                        !xpb
     3        centre(leftix(ix1,iy),leftiy(ix1,iy))*                     !xpb
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+                      !xpb
     4        centre(ix2,iy)*vol(ix1,iy1)+                               !xpb
     5        centre(leftix(ix2,iy),leftiy(ix2,iy))*                     !xpb
     5        vol(topix(ix1,iy1),topiy(ix1,iy1))+                        !xpb
     6        centre(leftix(ix1,iy1),leftiy(ix1,iy1))*                   !xpb
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))/                       !xpb
     1       (vol(ix2,iy)+vol(leftix(ix2,iy),leftiy(ix2,iy))+            !xpb
     3        vol(leftix(ix1,iy1),leftiy(ix1,iy1))+vol(ix1,iy1)+         !xpb
     5        vol(topix(ix1,iy1),topiy(ix1,iy1))+                        !xpb
     6        vol(leftix(ix1,iy),leftiy(ix1,iy)))                        !xpb
            corner(ix2,iy2)=corner(ix1,iy1)
          endif                                                          !xpb
        enddo                                                            !xpb
      enddo                                                              !xpb
#endif
      return
      end
