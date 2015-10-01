      subroutine map(nx,ny,vol,inp,outp)
      use b2mod_types
      use b2mod_indirect
      use b2mod_boundary_namelist
c map to north-east vertices
      integer nx, ny
      real (kind=R8) :: inp(-1:nx,-1:ny),outp(-1:nx,-1:ny),
     1 vol(-1:nx,-1:ny)
      integer ix, iy, nc, ix1, ix2
                    
c standard case
      do iy=-1,ny
        do ix=-1,nx
          if(rightix(ix,iy).ne.nx+1 .and. topiy(ix,iy).ne.ny+1) then
            outp(ix,iy)=
     1       (inp(ix,iy)*vol(topix(rightix(ix,iy),rightiy(ix,iy)),
     1                        topiy(rightix(ix,iy),rightiy(ix,iy)))+
     2        inp(topix(ix,iy),topiy(ix,iy))*
     2        vol(rightix(ix,iy),rightiy(ix,iy))+
     3        inp(rightix(ix,iy),rightiy(ix,iy))*
     3        vol(topix(ix,iy),topiy(ix,iy))+
     4        inp(topix(rightix(ix,iy),rightiy(ix,iy)),
     4            topiy(rightix(ix,iy),rightiy(ix,iy)))*vol(ix,iy))/
     1       (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy))+
     3        vol(rightix(ix,iy),rightiy(ix,iy))+
     4        vol(topix(rightix(ix,iy),rightiy(ix,iy)),
     4            topiy(rightix(ix,iy),rightiy(ix,iy))))
           elseif(rightix(ix,iy).eq.nx+1 .and. topiy(ix,iy).ne.ny+1) 
     >      then
c eastern edges
             outp(ix,iy)=
     1        (inp(ix,iy)*vol(topix(ix,iy),topiy(ix,iy))+
     2         inp(topix(ix,iy),topiy(ix,iy))*vol(ix,iy))/
     1        (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
           elseif(rightix(ix,iy).ne.nx+1 .and. topiy(ix,iy).eq.ny+1) 
     >      then
c northern edges
             outp(ix,iy)=
     1        (inp(ix,iy)*vol(rightix(ix,iy),rightiy(ix,iy))+
     2         inp(rightix(ix,iy),rightiy(ix,iy))*vol(ix,iy))/
     1        (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
           else
c north-east corners
             outp(ix,iy)=inp(ix,iy)
           endif
         enddo
      enddo

c isolating cuts      
      do nc=1,nniso
        ix=nxiso1(nc)
        do iy=nyiso1(nc),nyiso2(nc)-1
          outp(ix,iy)=
     1   (inp(ix,iy)*vol(topix(ix,iy),topiy(ix,iy))+
     1    inp(topix(ix,iy),topiy(ix,iy))*vol(ix,iy))/
     2           (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
        enddo
        iy=nyiso2(nc)
        outp(ix,iy)=inp(ix,iy)
        ix=nxiso2(nc)
        do iy=nyiso1(nc),nyiso2(nc)-1
          outp(ix,iy)=
     1   (inp(ix,iy)*vol(topix(ix,iy),topiy(ix,iy))+
     1    inp(topix(ix,iy),topiy(ix,iy))*vol(ix,iy))/
     2   (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))
        enddo
        iy=nyiso2(nc)
        outp(ix,iy)=inp(ix,iy)
      enddo

c special treatment of standard X-points
      do iy=0,ny-1
        do ix=0,nx-1
          if(rightix(ix,iy).ne.nx+1 .and.
     .       rightix(topix(ix,iy),topiy(ix,iy)).ne.
     .       topix(rightix(ix,iy),rightiy(ix,iy))) then
            ix1=ix
            ix2=leftix(bottomix(rightix(topix(ix1,iy),topiy(ix1,iy)),
     .                          rightiy(topix(ix1,iy),topiy(ix1,iy))),
     .                 bottomiy(rightix(topix(ix1,iy),topiy(ix1,iy)),
     .                          rightiy(topix(ix1,iy),topiy(ix1,iy))))
c check that this is a standard X-point
            if(ix2.gt.-1 .and. ix2.lt.nx) then
              outp(ix1,iy)=
     1         (inp(ix1,iy)*vol(ix2,iy)+
     2          inp(topix(ix1,iy),topiy(ix1,iy))*
     2          vol(topix(ix2,iy),topiy(ix2,iy))+
     3          inp(rightix(topix(ix1,iy),topiy(ix1,iy)),
     3              rightiy(topix(ix1,iy),topiy(ix1,iy)))*
     3          vol(rightix(topix(ix2,iy),topiy(ix2,iy)),
     3              rightiy(topix(ix2,iy),topiy(ix2,iy)))+
     4          inp(rightix(ix1,iy),rightiy(ix1,iy))*
     4          vol(rightix(ix2,iy),rightiy(ix2,iy))+
     5          inp(ix2,iy)*vol(ix1,iy)+
     6          inp(topix(ix2,iy),topiy(ix2,iy))*
     6          vol(topix(ix1,iy),topiy(ix1,iy))+
     7          inp(rightix(topix(ix2,iy),topiy(ix2,iy)),
     7              rightiy(topix(ix2,iy),topiy(ix2,iy)))*
     7          vol(rightix(topix(ix1,iy),topiy(ix1,iy)),
     7              rightiy(topix(ix1,iy),topiy(ix1,iy)))+
     8          inp(rightix(ix2,iy),rightiy(ix2,iy))*
     8          vol(rightix(ix1,iy),rightiy(ix1,iy)))
              outp(ix1,iy)=outp(ix1,iy)/
     1         (vol(ix1,iy)+vol(topix(ix1,iy),topiy(ix1,iy))+
     3          vol(rightix(topix(ix1,iy),topiy(ix1,iy)),
     3              rightiy(topix(ix1,iy),topiy(ix1,iy)))+
     4          vol(rightix(ix1,iy),rightiy(ix1,iy))+
     5          vol(ix2,iy)+vol(topix(ix2,iy),topiy(ix2,iy))+
     7          vol(rightix(topix(ix2,iy),topiy(ix2,iy)),
     7              rightiy(topix(ix2,iy),topiy(ix2,iy)))+
     8          vol(rightix(ix2,iy),rightiy(ix2,iy)))
              outp(ix2,iy)=outp(ix1,iy)
            endif
          endif
        enddo
      enddo

c special treatment of diverted island X-points
      do iy=-1,ny-1
        do ix=-1,nx-1
          if(topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1) then
            ix1=ix
            ix2=bottomix(leftix(topix(rightix(ix1,iy),rightiy(ix1,iy)),
     .                          topiy(rightix(ix1,iy),rightiy(ix1,iy))),
     .                   leftiy(topix(rightix(ix1,iy),rightiy(ix1,iy)),
     .                          topiy(rightix(ix1,iy),rightiy(ix1,iy))))
            outp(ix1,iy)=(inp(ix1,iy)*vol(topix(ix2,iy),topiy(ix2,iy))+
     2       inp(rightix(ix1,iy),rightiy(ix1,iy))*vol(ix2,iy)+
     3       inp(topix(rightix(ix1,iy),rightiy(ix1,iy)),
     3           topiy(rightix(ix1,iy),rightiy(ix1,iy)))*
     3       vol(rightix(ix2,iy),rightiy(ix2,iy))+
     4       inp(ix2,iy)*vol(rightix(ix1,iy),rightiy(ix1,iy))+
     5       inp(topix(ix2,iy),topiy(ix2,iy))*vol(ix1,iy)+
     6       inp(rightix(ix2,iy),rightiy(ix2,iy))*
     6       vol(topix(rightix(ix1,iy),rightiy(ix1,iy)),
     6           topiy(rightix(ix1,iy),rightiy(ix1,iy))))/
     1      (vol(ix1,iy)+vol(rightix(ix1,iy),rightiy(ix1,iy))+
     3       vol(topix(rightix(ix1,iy),rightiy(ix1,iy)),
     3           topiy(rightix(ix1,iy),rightiy(ix1,iy)))+
     4       vol(ix2,iy)+vol(topix(ix2,iy),topiy(ix2,iy))+
     6       vol(rightix(ix2,iy),rightiy(ix2,iy)))
            outp(ix2,iy)=outp(ix1,iy)
          endif
        enddo
      enddo

c special treatment of limiter tangency X-points
      do iy=-1,ny-1
        do ix=-1,nx-1
          if(bottomiy(ix,iy).eq.-2) cycle
          if(bottomiy(bottomix(ix,iy),bottomiy(ix,iy)).eq.iy 
     &     .and. ny.gt.1) then
            ix1=ix
            ix2=topix(leftix(bottomix(rightix(ix1,iy),rightiy(ix1,iy)),
     .                       bottomiy(rightix(ix1,iy),rightiy(ix1,iy))),
     .                leftiy(bottomix(rightix(ix1,iy),rightiy(ix1,iy)),
     .                       bottomiy(rightix(ix1,iy),rightiy(ix1,iy))))
            outp(ix1,iy)=
     1      (inp(ix1,iy)*vol(bottomix(ix2,iy),bottomiy(ix2,iy))+
     2       inp(rightix(ix1,iy),rightiy(ix1,iy))*vol(ix2,iy)+
     3       inp(bottomix(rightix(ix1,iy),rightiy(ix1,iy)),
     3           bottomiy(rightix(ix1,iy),rightiy(ix1,iy)))*
     3       vol(rightix(ix2,iy),rightiy(ix2,iy))+
     4       inp(ix2,iy)*vol(rightix(ix1,iy),rightiy(ix1,iy))+
     5       inp(bottomix(ix2,iy),bottomiy(ix2,iy))*vol(ix1,iy)+
     6       inp(rightix(ix2,iy),rightiy(ix2,iy))*
     6       vol(bottomix(rightix(ix1,iy),rightiy(ix1,iy)),
     6           bottomiy(rightix(ix1,iy),rightiy(ix1,iy))))/
     1      (vol(ix1,iy)+vol(rightix(ix1,iy),rightiy(ix1,iy))+
     3       vol(bottomix(rightix(ix1,iy),rightiy(ix1,iy)),
     3           bottomiy(rightix(ix1,iy),rightiy(ix1,iy)))+
     4       vol(ix2,iy)+vol(bottomix(ix2,iy),bottomiy(ix2,iy))+
     6       vol(rightix(ix2,iy),rightiy(ix2,iy)))
            outp(ix2,iy)=outp(ix1,iy)
          endif
        enddo
      enddo

      return
      end
