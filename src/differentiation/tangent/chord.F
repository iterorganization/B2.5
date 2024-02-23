c
      function chord_value(nCv,
     1 tt,
     2 xp,yp,zp,
     3 dxp,dyp,dzp,vector)

      use b2mod_types
      use b2us_data_diffv
      use b2mod_constants
      use b2mod_subsys
      implicit none
      
      real (kind=R8) ::
     * chord_value
      integer nCv
      real (kind=R8) ::
     * tt(nCv)
      real (kind=R8) ::
     * xp,yp,zp,dxp,dyp,dzp
      logical vector

      real (kind=R8) ::
     * x,y,z,dl,dx,dy,dz,dbr,dbz,dbphi,dbx,dby,dbs,dotprod,
     2 phi,rt,zt
      integer iCv
      logical found

      real (kind=R8) ::
     * rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
      
      intrinsic sqrt,sin,cos,abs,atan2
      external chord_find, xertst

      call subini ('chord_value')
      x=xp
      y=yp
      z=zp
      dx=dxp
      dy=dyp
      dz=dzp

      chord_value=0.0_R8

      call chord_find(mpg,geo,x,y,z,found,iCv)

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
#ifdef WG_TODO
c in r,z,phi
        dbr=((crx(ix,iy,1)+crx(ix,iy,3))-
     .       (crx(ix,iy,2)+crx(ix,iy,0)))/2.0_R8
        dbz=((cry(ix,iy,1)+cry(ix,iy,3))-
     .       (cry(ix,iy,2)+cry(ix,iy,0)))/2.0_R8
#endif
        dbs=sqrt(dbr**2+dbz**2)
        dbr=dbr/dbs*pitch(geo,iCv)
        dbz=dbz/dbs*pitch(geo,iCv)
        dbphi=1.0*sqrt(1.0-pitch(geo,iCv)**2)
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
c       write(*,*) 'iCv = ',iCv
c       write(*,*) 'dl = ',dx,dy,dz
c       write(*,*) 'bl = ',dbr,dbz,dbphi
c       write(*,*) 'phi = ',phi
c       write(*,*) 'bl = ',dbx,dby,dbz
c       write(*,*) 'dl.bl = ',dotprod
        chord_value=tt(iCv)*dotprod
      else
        chord_value=tt(iCv)
      endif
      call subend ()
      return
      end function chord_value
c
      function chord_integral(nCv,tt,npoints,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector)
      
      use b2mod_types
      use b2mod_subsys
      implicit none
      
      real (kind=R8) ::
     * chord_integral 
      integer nCv,npoints
      real (kind=R8) ::
     * tt(nCv)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      logical vector

      real (kind=R8) ::
     * chord_value
      integer i
      real (kind=R8) ::
     * f,x,y,z,dl,dx,dy,dz

      intrinsic sqrt

      call subini('chord_integral')
      chord_integral=0.0_R8

      dx=xp2-xp1
      dy=yp2-yp1
      dz=zp2-zp1
      dl=sqrt(dx**2+dy**2+dz**2)
      if(dl.eq.0.0_R8) return
      dx=dx/dl
      dy=dy/dl
      dz=dz/dl
      dl=dl/real(npoints)

      do i=1,npoints
        f=(real(i,R8)-0.5_R8)/real(npoints,R8)
        x=xp1*(1.0_R8-f)+xp2*f
        y=yp1*(1.0_R8-f)+yp2*f
        z=zp1*(1.0_R8-f)+zp2*f
        chord_integral=chord_integral+dl*
     1   chord_value(nCv,tt,x,y,z,dx,dy,dz,vector)
      enddo

!DIAG      call chord_diagnostics      ! diag

      call subend()
      return
      end function chord_integral
c
      subroutine chord_compile(npoints,mpg,geo,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector,
     2 indchd,indCv,indscl)

      use b2mod_types
      use b2us_map_diffv
      use b2us_geo_diffv
      use b2mod_subsys
      implicit none

      integer npoints
      type(mapping) :: mpg
      type(geometry) :: geo
      real (kind=R8) ::
     * tt(mpg%nCv),score(mpg%nCv)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      integer indchd,indCv(*)
      real (kind=R8) ::
     * indscl(*)
      logical vector

      real (kind=R8) ::
     * chord_value
      integer i,iCv
      real (kind=R8) ::
     * f,x,y,z,dl,dxp,dyp,dzp
      logical found

      call subini('chord_compile')
      tt=1.0_R8
      score=0.0_R8

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
        call chord_find(mpg,geo,x,y,z,found,iCv)
        if(found) then
          score(iCv)=score(iCv)+dl*
     1     chord_value(mpg%nCv,tt,x,y,z,dxp,dyp,dzp,vector)
        endif
      enddo

      i=0
      do iCv=1,mpg%nCv
        if(score(iCv).ne.0.0_R8) then
          i=i+1
          indCv(i)=iCv
          indscl(i)=score(iCv)
        endif
      enddo
      indchd=i

      call subend()
      return
      end subroutine chord_compile
c
      subroutine chord_find(mpg,geo,xp,yp,zp,found,cvpos)

      use b2mod_types
      use b2us_map_diffv
      use b2us_geo_diffv
      use b2mod_subsys
      implicit none

!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag

      type(mapping) :: mpg
      type(geometry) :: geo
      real (kind=R8) ::
     *   xp,yp,zp
      logical found
      integer cvpos

#ifdef WG_TODO
      logical inside_rectangle
#endif
      integer iCv,icnt1,icnt2
      real (kind=R8) ::
     *   x,y,z,r,rmin,rmax,zmin,zmax

      logical first
      integer ir,iz
      integer lngcov,lngind
      parameter(lngcov=200,lngind=500 000)
      integer covered(lngcov,lngcov),
     1 cvcov(lngind),icov(lngind),
     2 ipos,jpos

#ifdef EASY
      integer ir1,ir2,iz1,iz2
#else
      integer ip
#ifdef WG_TODO
      integer irmin,irmax,izmin,izmax
      integer next(0:3)
      real (kind=R8) ::
     *   r1,r2,z1,z2
      data next/1,3,0,2/
#endif
      integer xmin(lngcov),xmax(lngcov),ymin,ymax
#endif

      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0

      external xerrab
      save covered,rmin,rmax,zmin,zmax,first,cvcov,icov,iCv

      data first/.true./

      call subini ('chord_find')
      if(first) then
        first=.false.
        rmin=minval(geo%vxX)
        rmax=maxval(geo%vxX)
        zmin=minval(geo%vxY)
        zmax=maxval(geo%vxY)
        write(*,*) 'rmin,rmax,zmin,zmax ',rmin,rmax,zmin,zmax

        do iz=1,lngcov
          do ir=1,lngcov
            covered(ir,iz)=0
          enddo
        enddo
        ipos=0
        icnt2=0

        do iCv=1,mpg%nCi
#ifdef EASY
#ifdef WG_TODO
          ir1= irpos(min(crx(ix,iy,0),crx(ix,iy,1),
     2                   crx(ix,iy,2),crx(ix,iy,3)))
          ir2= irpos(max(crx(ix,iy,0),crx(ix,iy,1),
     2                   crx(ix,iy,2),crx(ix,iy,3)))
          iz1= izpos(min(cry(ix,iy,0),cry(ix,iy,1),
     2                   cry(ix,iy,2),cry(ix,iy,3)))
          iz2= izpos(max(cry(ix,iy,0),cry(ix,iy,1),
     2                   cry(ix,iy,2),cry(ix,iy,3)))
#endif
          do iz=iz1,iz2
            do ir=ir1,ir2
              ipos=ipos+1
              if(ipos.gt.lngind) then
                write(*,*) iCv,ipos
                call xerrab ('increase size of {xyi}cov')
              endif
              cvcov(ipos)=iCv
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
#ifdef WG_TODO
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
#endif
          do iz=ymin,ymax
            do ir=xmin(iz),xmax(iz)
              ipos=ipos+1
              if(ipos.gt.lngind) then
                write(*,*) 'iCv ',iCv
                call xerrab ('increase lngind')
              endif
              cvcov(ipos)=iCv
              icov(ipos)=0
              icnt2=1
              if(covered(ir,iz).eq.0) then
                covered(ir,iz)=ipos
              else
                jpos=covered(ir,iz)
                icnt2=icnt2+1
   20           if(icov(jpos).ne.0) then
                  jpos=icov(jpos)
                  icnt2=icnt2+1
                  goto 20
                endif
                icov(jpos)=ipos
              endif
            enddo
          enddo
#endif
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
        iCv=1
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
#ifdef WG_TODO
        if(inside_rectangle(r,z,
     1   crx(ix,iy,1),cry(ix,iy,1),
     2   crx(ix,iy,3),cry(ix,iy,3),
     3   crx(ix,iy,2),cry(ix,iy,2),
     4   crx(ix,iy,0),cry(ix,iy,0))) goto 2100
#endif
      endif
c next try the chain
      do while(jpos.ne.0)
!DIAG   chain=chain+1      ! diag
        iCv=cvcov(jpos)
        jpos=icov(jpos)
#ifdef WG_TODO
        if(inside_rectangle(r,z,
     1   crx(ix,iy,1),cry(ix,iy,1),
     2   crx(ix,iy,3),cry(ix,iy,3),
     3   crx(ix,iy,2),cry(ix,iy,2),
     4   crx(ix,iy,0),cry(ix,iy,0))) goto 2100
#endif
      enddo
      call subend ()
      return

 2100 continue
      cvpos=iCv
      found=.true.
      call subend ()
      return

      contains

      integer function irpos(r)
      implicit none
      real(kind=R8), intent(in) :: r
      irpos = int((r-rmin)/(rmax-rmin)*lngcov*0.999999_R8+1)
      return
      end function irpos

      integer function izpos(z)
      implicit none
      real(kind=R8), intent(in) :: z
      izpos = int((z-zmin)/(zmax-zmin)*lngcov*0.999999_R8+1)
      return
      end function izpos

      function rpos(ir)
      implicit none
      integer, intent(in) :: ir
      real(kind=R8) :: rpos
      rpos = rmin+
     .       real(ir-1,R8)*(rmax-rmin)/(real(lngcov,R8)*0.999999_R8)
      return
      end function rpos

      function zpos(iz)
      implicit none
      integer, intent(in) :: iz
      real(kind=R8) :: zpos
      zpos = zmin+
     .       real(iz-1,R8)*(zmax-zmin)/(real(lngcov,R8)*0.999999_R8)
      return
      end function zpos

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

