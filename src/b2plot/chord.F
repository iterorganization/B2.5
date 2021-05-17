c
      function chord_value(nx,ny,
     1 tt,crx,cry,bb,
     2 xp,yp,zp,
     3 dxp,dyp,dzp,vector)

      use b2mod_types
      use b2mod_constants
      implicit none
      
      real (kind=R8) ::
     * chord_value
      integer nx,ny
      real (kind=R8) ::
     * tt(-1:nx,-1:ny),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp,yp,zp,dxp,dyp,dzp
      logical vector

      real (kind=R8) ::
     * x,y,z,dl,dx,dy,dz,dbr,dbz,dbphi,dbx,dby,dbs,dotprod,
     2 phi,rt,zt
      integer ix,iy
      logical found

      real (kind=R8) ::
     * rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
      
      real (kind=R8) :: pit
      intrinsic sqrt,sin,cos,abs,atan2
      external chord_find, subini, subend, xertst

      pit(ix,iy)=sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3)

      call subini ('chord_value')
      x=xp
      y=yp
      z=zp
      dx=dxp
      dy=dyp
      dz=dzp

      chord_value=0.0_R8

      call chord_find(nx,ny,crx,cry,x,y,z,found,ix,iy)

      if(.not.found) then
        if(.not.vector.and.a0.gt.0.0_R8.and.b0.gt.0.0_R8) then
          zt=z+zshift
          rt=sqrt(x**2+y**2)+rshift
          if(((rt-r0)/a0)**2+((zt-z0)/b0)**2.le.1.0_R8) then
            chord_value=tt0
          endif
        endif
        call subend ()
        return
      endif

      if(vector) then
c in r,z,phi
        dbr=((crx(ix,iy,1)+crx(ix,iy,3))-
     .       (crx(ix,iy,2)+crx(ix,iy,0)))/2.0_R8
        dbz=((cry(ix,iy,1)+cry(ix,iy,3))-
     .       (cry(ix,iy,2)+cry(ix,iy,0)))/2.0_R8
        dbs=sqrt(dbr**2+dbz**2)
        dbr=dbr/dbs*pit(ix,iy)
        dbz=dbz/dbs*pit(ix,iy)
        dbphi=1.0*sqrt(1.0-pit(ix,iy)**2)
c in x,y,z
        phi=atan2(y,x)
        dbz=dbz
        dbx=dbr*cos(phi)+dbphi*cos(phi+pi/2.0_R8)
        dby=dbr*sin(phi)+dbphi*sin(phi+pi/2.0_R8)
c assert
        dotprod=dbx*dbx+dby*dby+dbz*dbz
        if(abs(dotprod-1.0).gt.1e-10_R8) then
          write(*,*) 'Assertion failed, norm of unit vector = ',
     1     dotprod
        endif
c end-assertion
        dl=sqrt(dx**2+dy**2+dz**2)
        call xertst(dl.ne.0., 'dl in chord.F == 0')
        dx=dx/dl
        dy=dy/dl
        dz=dz/dl
        dotprod=dbx*dx+dby*dy+dbz*dz
c       write(*,*) 'ix,iy = ',ix,iy
c       write(*,*) 'dl = ',dx,dy,dz
c       write(*,*) 'bl = ',dbr,dbz,dbphi
c       write(*,*) 'phi = ',phi
c       write(*,*) 'bl = ',dbx,dby,dbz
c       write(*,*) 'dl.bl = ',dotprod
        chord_value=tt(ix,iy)*dotprod
      else
        chord_value=tt(ix,iy)
      endif
      call subend ()
      return
      end function chord_value
c
      function chord_integral(nx,ny,tt,crx,cry,bb,npoints,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector)
      
      use b2mod_types
      implicit none
      
      real (kind=R8) ::
     * chord_integral 
      integer nx,ny,npoints
      real (kind=R8) ::
     * tt(-1:nx,-1:ny),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      logical vector

      real (kind=R8) ::
     * chord_value
      integer i
      real (kind=R8) ::
     * f,x,y,z,dl,dx,dy,dz

      intrinsic sqrt
      external subini, subend

      call subini('chord_integral')
      chord_integral=0.0

      dx=xp2-xp1
      dy=yp2-yp1
      dz=zp2-zp1
      dl=sqrt(dx**2+dy**2+dz**2)
      if(dl.eq.0.0) return
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dl=dl/real(npoints)

      do i=1,npoints
        f=real(i-0.5)/real(npoints)
        x=xp1*(1.0-f)+xp2*f
        y=yp1*(1.0-f)+yp2*f
        z=zp1*(1.0-f)+zp2*f
        chord_integral=chord_integral+dl*
     1   chord_value(nx,ny,tt,crx,cry,bb,x,y,z,dx,dy,dz,vector)
      enddo

!DIAG      call chord_diagnostics      ! diag

      call subend()
      return
      end function chord_integral
c
      subroutine chord_compile(npoints,nx,ny,crx,cry,bb,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector,
     2 indchd,indix,indiy,indscl)

      use b2mod_types
      implicit none

      integer npoints,nx,ny
      real (kind=R8) ::
     * tt(-1:nx,-1:ny),score(-1:nx,-1:ny),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      integer indchd,indix(*),indiy(*)
      real (kind=R8) ::
     * indscl(*)
      logical vector

      real (kind=R8) ::
     * chord_value
      integer i,ix,iy
      real (kind=R8) ::
     * f,x,y,z,dl,dxp,dyp,dzp
      logical found
      external subini, subend

      call subini('chord_compile')
      do iy=-1,ny
        do ix=-1,nx
          tt(ix,iy)=1.0_R8
          score(ix,iy)=0.0_R8
        enddo
      enddo

      dxp=xp2-xp1
      dyp=yp2-yp1
      dzp=zp2-zp1
      dl=sqrt(dxp**2+dyp**2+dzp**2)
      dxp=dxp/dl
      dyp=dyp/dl
      dzp=dzp/dl
      dl=dl/real(npoints)

      do i=1,npoints
        f=real(i-0.5)/real(npoints)
        x=xp1*(1.0-f)+xp2*f
        y=yp1*(1.0-f)+yp2*f
        z=zp1*(1.0-f)+zp2*f
        call chord_find(nx,ny,crx,cry,x,y,z,found,ix,iy)
        if(found) then
          score(ix,iy)=score(ix,iy)+dl*
     1     chord_value(nx,ny,tt,crx,cry,bb,
     2                 x,y,z,dxp,dyp,dzp,vector)
        endif
      enddo

      i=0
      do iy=-1,ny
        do ix=-1,nx
          if(score(ix,iy).ne.0) then
            i=i+1
            indix(i)=ix
            indiy(i)=iy
            indscl(i)=score(ix,iy)
          endif
        enddo
      enddo
      indchd=i

      call subend()
      return
      end subroutine chord_compile
c
      subroutine chord_find(nx,ny,
     1 crx,cry,xp,yp,zp,found,ixpos,iypos)

      use b2mod_types
      use b2mod_indirect
      implicit none

!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag

      integer nx,ny
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3),cry(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *   xp,yp,zp
      logical found
      integer ixpos,iypos

      logical inside_rectangle
      integer i,ix,iy,icnt1,icnt2
      real (kind=R8) ::
     *   x,y,z,r,rmin,rmax,zmin,zmax

      logical first
      integer ir,iz,irpos,izpos
      real (kind=R8) ::
     *   rpos,zpos
      integer lngcov,lngind
      parameter(lngcov=200,lngind=500 000)
      integer covered(lngcov,lngcov),
     1 xcov(lngind),ycov(lngind),icov(lngind),
     2 ipos,jpos

#ifdef EASY
      integer ir1,ir2,iz1,iz2
#else
      integer ip,irmin,irmax,izmin,izmax
      integer next(0:3),
     1 xmin(lngcov),xmax(lngcov),ymin,ymax
      real (kind=R8) ::
     *   r1,r2,z1,z2
      data next/1,3,0,2/
#endif

      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0

      external subini, subend, xerrab
      save covered,rmin,rmax,zmin,zmax,first,xcov,ycov,icov,ix,iy

      data first/.true./

      irpos(r)=(r-rmin)/(rmax-rmin)*lngcov*0.999999_R8+1
      izpos(z)=(z-zmin)/(zmax-zmin)*lngcov*0.999999_R8+1
      rpos(ir)=rmin+real(ir-1)*(rmax-rmin)/(real(lngcov)*0.999999_R8)
      zpos(iz)=zmin+real(iz-1)*(zmax-zmin)/(real(lngcov)*0.999999_R8)

      call subini ('chord_find')
      if(first) then
        first=.false.
        rmin=crx(0,0,0)
        rmax=crx(0,0,0)
        zmin=cry(0,0,0)
        zmax=cry(0,0,0)
        do i=0,3
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &           bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &           region(ix,iy,0).ne.0) then
                rmin=min(rmin,crx(ix,iy,i))
                rmax=max(rmax,crx(ix,iy,i))
                zmin=min(zmin,cry(ix,iy,i))
                zmax=max(zmax,cry(ix,iy,i))
              endif
            enddo
          enddo
        enddo
        write(*,*) 'rmin,rmax,zmin,zmax ',rmin,rmax,zmin,zmax

        do iz=1,lngcov
          do ir=1,lngcov
            covered(ir,iz)=0
          enddo
        enddo
        ipos=0
        icnt2=0

        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
#ifdef EASY
              ir1=
     1         irpos(min(crx(ix,iy,0),crx(ix,iy,1),
     2                   crx(ix,iy,2),crx(ix,iy,3)))
              ir2=
     1         irpos(max(crx(ix,iy,0),crx(ix,iy,1),
     2                   crx(ix,iy,2),crx(ix,iy,3)))
              iz1=
     1         izpos(min(cry(ix,iy,0),cry(ix,iy,1),
     2                   cry(ix,iy,2),cry(ix,iy,3)))
              iz2=
     1         izpos(max(cry(ix,iy,0),cry(ix,iy,1),
     2                   cry(ix,iy,2),cry(ix,iy,3)))
              do iz=iz1,iz2
                do ir=ir1,ir2
                  ipos=ipos+1
                  if(ipos.gt.lngind) then
                    write(*,*) ix,iy,ipos
                    call xerrab ('increase size of {xyi}cov')
                  endif
                  xcov(ipos)=ix
                  ycov(ipos)=iy
                  icov(ipos)=0
                  if(covered(ir,iz).eq.0) then
                    covered(ir,iz)=ipos
                  else
                    jpos=covered(ir,iz)
                    icnt1=1
                    do while(icov(jpos).ne.0)
                      jpos=icov(jpos)
                      icnt1=icnt1+1
                    enddo
                    icov(jpos)=ipos
                    icnt1=icnt1+1
                    icnt2=max(icnt1,icnt2)
                  endif
                enddo
              enddo
#else
              do ip=1,lngcov
                xmin(ip)=lngcov+1
                xmax(ip)=0
              enddo
              ymin=lngcov+1
              ymax=0
              do ip=0,3
                r1=crx(ix,iy,ip)
                z1=cry(ix,iy,ip)
                r2=crx(ix,iy,next(ip))
                z2=cry(ix,iy,next(ip))
                ir=irpos(r1)
                iz=izpos(z1)
c assert ir,iz is inside
                ymin=min(ymin,iz)
                ymax=max(ymax,iz)
                xmin(iz)=min(xmin(iz),ir)
                xmax(iz)=max(xmax(iz),ir)
c look for intersections for integer ir
c       write(*,*) 'r1,r2 ',r1,r2
                irmin=irpos(min(r1,r2))
                irmax=irpos(max(r1,r2))
c       write(*,*) 'irmin,irmax ',irmin,irmax
                do ir=irmin+1,irmax
                  r=rpos(ir)
                  z=z1+(z2-z1)/(r2-r1)*(r-r1)
                  iz=izpos(z)
c       write(*,*) 'ir,r,z,iz ',ir,r,z,iz
c assert ir-1,iz and ir,iz covered
                  ymin=min(ymin,iz)
                  ymax=max(ymax,iz)
                  xmin(iz)=min(xmin(iz),ir-1)
                  xmax(iz)=max(xmax(iz),ir)
                enddo
c look for intersections for integer iz
c       write(*,*) 'z1,z2 ',z1,z2
                izmin=izpos(min(z1,z2))
                izmax=izpos(max(z1,z2))
c       write(*,*) 'izmin,izmax ',izmin,izmax
                do iz=izmin+1,izmax
                  z=zpos(iz)
                  r=r1+(r2-r1)/(z2-z1)*(z-z1)
                  ir=irpos(r)
c         write(*,*) 'iz,z,r,ir ',iz,z,r,ir
c assert ir,iz-1 and ir,iz covered
                  ymin=min(ymin,iz-1)
                  ymax=max(ymax,iz)
                  xmin(iz-1)=min(xmin(iz-1),ir)
                  xmax(iz-1)=max(xmax(iz-1),ir)
                  xmin(iz)=min(xmin(iz),ir)
                  xmax(iz)=max(xmax(iz),ir)
                enddo
c       write(*,*) 'ymin,ymax ',ymin,ymax
c       write(*,*) 'xmin ',xmin
c       write(*,*) 'xmax ',xmax
              enddo
              do iz=ymin,ymax
                do ir=xmin(iz),xmax(iz)
                  ipos=ipos+1
                  if(ipos.gt.lngind) then
                    write(*,*) 'ix,iy ',ix,iy
                    call xerrab ('increase lngind')
                  endif
                  xcov(ipos)=ix
                  ycov(ipos)=iy
                  icov(ipos)=0
                  icnt2=1
                  if(covered(ir,iz).eq.0) then
                    covered(ir,iz)=ipos
                  else
                    jpos=covered(ir,iz)
                    icnt2=icnt2+1
   20               if(icov(jpos).ne.0) then
                      jpos=icov(jpos)
                      icnt2=icnt2+1
                      goto 20
                    endif
                    icov(jpos)=ipos
                  endif
                enddo
              enddo
#endif
            endif
          enddo
        enddo

        icnt1=0
        do iz=1,lngcov
          do ir=1,lngcov
            if(covered(ir,iz).ne.0) icnt1=icnt1+1
          enddo
        enddo
        write(*,*) 100.0*icnt1/(lngcov*lngcov),' % covered'
        write(*,*) ipos,' positions used'
        write(*,*) icnt2,' is the maximum chain length'
        ix=0
        iy=ny-1
!DIAG   rectangles=0      ! diag
!DIAG   triangles=0      ! diag
!DIAG   chain=0      ! diag
      endif

      found=.false.

      x=xp
      y=yp
      z=zp+zshift
      r=sqrt(x**2+y**2)+rshift
      if(r.lt.rmin.or.r.gt.rmax.or.z.lt.zmin.or.z.gt.zmax) then
        call subend ()
        return
      end if
      jpos=covered(irpos(r),izpos(z))
      if(jpos.ne.0) then
        if(inside_rectangle(r,z,
     1   crx(ix,iy,1),cry(ix,iy,1),
     2   crx(ix,iy,3),cry(ix,iy,3),
     3   crx(ix,iy,2),cry(ix,iy,2),
     4   crx(ix,iy,0),cry(ix,iy,0))) goto 2100
      endif
c next try the chain
      do while(jpos.ne.0)
!DIAG   chain=chain+1      ! diag
        ix=xcov(jpos)
        iy=ycov(jpos)
        jpos=icov(jpos)
        if(inside_rectangle(r,z,
     1   crx(ix,iy,1),cry(ix,iy,1),
     2   crx(ix,iy,3),cry(ix,iy,3),
     3   crx(ix,iy,2),cry(ix,iy,2),
     4   crx(ix,iy,0),cry(ix,iy,0))) goto 2100
      enddo
      call subend ()
      return

 2100 continue
      ixpos=ix
      iypos=iy
      found=.true.
      call subend ()
      return
      end subroutine chord_find
c
      function inside_rectangle(r,z,x1,y1,x2,y2,x3,y3,x4,y4)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *  r,z,x1,y1,x2,y2,x3,y3,x4,y4
      logical inside_rectangle,inside_triangle
!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag
!DIAG      rectangles=rectangles+1      ! diag
      inside_rectangle=.true.
c       rmin=min(x1,x2,x3,x4))
c       rmax=max(x1,x2,x3,x4))
c       zmin=min(y1,y2,y3,y4))
c       zmax=max(y1,y2,y3,y4))
c       if(r.lt.rmin.or.r.gt.rmax.or.z.lt.zmin.or.z.gt.zmax) then
c         inside_rectangle=.false.
c         return
c       endif
      if(inside_triangle(r,z,x1,y1,x2,y2,x3,y3)) return
      if(inside_triangle(r,z,x3,y3,x4,y4,x1,y1)) return
      inside_rectangle=.false.
      return
      end function inside_rectangle
c
      function inside_triangle(r,z,x1,y1,x2,y2,x3,y3)
      use b2mod_types
      implicit none
      logical inside_triangle
      real (kind=R8) ::
     *  r,z,x1,y1,x2,y2,x3,y3,s1,s2,s3
!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag
!DIAG      triangles=triangles+1      ! diag
      s1=(x2-x1)*(z-y1)-(y2-y1)*(r-x1)
      s2=(x3-x2)*(z-y2)-(y3-y2)*(r-x2)
      s3=(x1-x3)*(z-y3)-(y1-y3)*(r-x3)
      inside_triangle=((s1.le.0.and.s2.le.0.and.s3.le.0).or.
     1 (s1.ge.0.and.s2.ge.0.and.s3.ge.0))
      return
      end function inside_triangle

      block data bd_chord_shift
      use b2mod_types
      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
      data rshift, zshift, r0, z0, a0, b0, tt0
     1    /0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8, 0.0_R8/
      end block data bd_chord_shift

!DIAG      subroutine chord_diagnostics      ! diag
!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag
!DIAG      write(*,*) chain,rectangles,triangles      ! diag
!DIAG      end subroutine chord_diagnostics      ! diag

