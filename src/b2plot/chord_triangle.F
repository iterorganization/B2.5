c
      function chord_tria_value(ntrii,nrknot,nx,ny,tt,
     1 xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     2 xp,yp,zp,dxp,dyp,dzp,vector)

      use b2mod_types
      use b2mod_constants
      implicit none
      
      real (kind=R8) ::
     * chord_tria_value
      integer nx,ny,ntrii,nrknot
      integer necke(3,ntrii),ixtri(ntrii),iytri(ntrii)
      real (kind=R8) ::
     * tt(1:ntrii),
     1 xtrian(nrknot),
     2 ytrian(nrknot),
     3 crx(-1:nx,-1:ny,0:3),
     4 cry(-1:nx,-1:ny,0:3),
     5 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp,yp,zp,dxp,dyp,dzp
      logical vector

      real (kind=R8) ::
     * x,y,z,dl,dx,dy,dz,dbr,dbz,dbphi,dbx,dby,dbs,dotprod,
     2 phi,rt,zt
      integer it
      logical found

      real (kind=R8) ::
     * rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0
      
      integer ix, iy
      integer ix_e2b
      real (kind=R8) :: pit
      intrinsic sqrt,sin,cos,abs,atan2
      external chord_find, subini, subend, xertst, ix_e2b

      pit(ix,iy)=sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3)

      call subini ('chord_tria_value')
      x=xp
      y=yp
      z=zp
      dx=dxp
      dy=dyp
      dz=dzp

      chord_tria_value=0.0_R8

      call chord_tria_find(ntrii,nrknot,necke,xtrian,ytrian,
     1 x,y,z,found,it)

      if(.not.found) then
        if(.not.vector.and.a0.gt.0.0_R8.and.b0.gt.0.0_R8) then
          zt=z+zshift
          rt=sqrt(x**2+y**2)+rshift
          if(((rt-r0)/a0)**2+((zt-z0)/b0)**2.le.1.0_R8) then
            chord_tria_value=tt0
          endif
        endif
        call subend ()
        return
      endif

      if(vector) then
        if (ixtri(it).ne.-1 .and. iytri(it).ne.-1) then
c in r,z,phi
          dbr=((crx(ix_e2b(ixtri(it)),iytri(it)-1,1)+
     .          crx(ix_e2b(ixtri(it)),iytri(it)-1,3))-
     .         (crx(ix_e2b(ixtri(it)),iytri(it)-1,2)+
     .          crx(ix_e2b(ixtri(it)),iytri(it)-1,0)))/2.0_R8
          dbz=((cry(ix_e2b(ixtri(it)),iytri(it)-1,1)+
     .          cry(ix_e2b(ixtri(it)),iytri(it)-1,3))-
     .         (cry(ix_e2b(ixtri(it)),iytri(it)-1,2)+
     .          cry(ix_e2b(ixtri(it)),iytri(it)-1,0)))/2.0_R8
          dbs=sqrt(dbr**2+dbz**2)
          dbr=dbr/dbs*pit(ix_e2b(ixtri(it)),iytri(it)-1)
          dbz=dbz/dbs*pit(ix_e2b(ixtri(it)),iytri(it)-1)
          dbphi=1.0*sqrt(1.0-pit(ix_e2b(ixtri(it)),iytri(it)-1)**2)
c in x,y,z
          phi=atan2(y,x)
          dbz=dbz
          dbx=dbr*cos(phi)+dbphi*cos(phi+pi/2.0_R8)
          dby=dbr*sin(phi)+dbphi*sin(phi+pi/2.0_R8)
c assert
          dotprod=dbx*dbx+dby*dby+dbz*dbz
          if(abs(dotprod-1.0).gt.1.0e-10_R8) then
            write(*,*) 'Assertion failed, norm of unit vector = ',
     1       dotprod
          endif
c end-assertion
          dl=sqrt(dx**2+dy**2+dz**2)
          call xertst(dl.ne.0., 'dl in chord.F == 0')
          dx=dx/dl
          dy=dy/dl
          dz=dz/dl
          dotprod=dbx*dx+dby*dy+dbz*dz
c         write(*,*) 'it,ix,iy = ',it,ix_e2b(ixtri(it)),iytri(it)-1
c         write(*,*) 'dl = ',dx,dy,dz
c         write(*,*) 'bl = ',dbr,dbz,dbphi
c         write(*,*) 'phi = ',phi
c         write(*,*) 'bl = ',dbx,dby,dbz
c         write(*,*) 'dl.bl = ',dotprod
          chord_tria_value=tt(it)*dotprod
        endif
      else
        chord_tria_value=tt(it)
      endif
      call subend ()
      return
      end function chord_tria_value
c
      function chord_tria_integral(ntrii,nrknot,nx,ny,tt,
     2 xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,npoints,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector)
      
      use b2mod_types
      implicit none
      
      real (kind=R8) ::
     * chord_tria_integral 
      integer nx,ny,ntrii,nrknot,npoints
      integer necke(3,ntrii),ixtri(ntrii),iytri(ntrii)
      real (kind=R8) ::
     * tt(1:ntrii),
     1 xtrian(nrknot),
     2 ytrian(nrknot),
     3 crx(-1:nx,-1:ny,0:3),
     4 cry(-1:nx,-1:ny,0:3),
     5 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      logical vector

      integer i
      real (kind=R8) ::
     * f,x,y,z,dl,dx,dy,dz

      real (kind=R8) ::
     * chord_tria_value
      external chord_tria_value
      intrinsic sqrt
      external subini, subend

      call subini('chord_tria_integral')
      chord_tria_integral=0.0

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
        chord_tria_integral=chord_tria_integral+dl*
     1   chord_tria_value(ntrii,nrknot,nx,ny,tt,
     2                    xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     2                    x,y,z,dx,dy,dz,vector)
      enddo

!DIAG      call chord_diagnostics      ! diag

      call subend()
      return
      end function chord_tria_integral
c
      subroutine chord_tria_compile(npoints,ntrii,nrknot,nx,ny,
     * xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     1 xp1,yp1,zp1,xp2,yp2,zp2,vector,
     2 indchd,indit,indscl)

      use b2mod_types
      implicit none

      integer npoints,ntrii,nrknot,nx,ny
      integer necke(3,ntrii),ixtri(ntrii),iytri(ntrii)
      real (kind=R8) ::
     * tt(ntrii),score(ntrii),
     * xtrian(nrknot),
     * ytrian(nrknot),
     1 crx(-1:nx,-1:ny,0:3),
     2 cry(-1:nx,-1:ny,0:3),
     3 bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     * xp1,yp1,zp1,xp2,yp2,zp2
      integer indchd,indit(*)
      real (kind=R8) ::
     * indscl(*)
      logical vector

      real (kind=R8) ::
     * chord_tria_value
      integer i,it
      real (kind=R8) ::
     * f,x,y,z,dl,dxp,dyp,dzp
      logical found
      external subini, subend

      call subini('chord_tria_compile')
      do it=1,ntrii
        tt(it)=1.0_R8
        score(it)=0.0_R8
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
        f=real(i-0.5_R8)/real(npoints)
        x=xp1*(1.0_R8-f)+xp2*f
        y=yp1*(1.0_R8-f)+yp2*f
        z=zp1*(1.0_R8-f)+zp2*f
        call chord_tria_find(ntrii,nrknot,necke,xtrian,ytrian,
     &   x,y,z,found,it)
        if(found) then
          score(it)=score(it)+dl*
     1     chord_tria_value(ntrii,nrknot,nx,ny,tt,
     2      xtrian,ytrian,necke,ixtri,iytri,crx,cry,bb,
     2      x,y,z,dxp,dyp,dzp,vector)
        endif
      enddo

      i=0
      do it=1,ntrii
        if(score(it).ne.0.0_R8) then
            i=i+1
            indit(i)=it
            indscl(i)=score(it)
        endif
      enddo
      indchd=i

      call subend()
      return
      end subroutine chord_tria_compile
c
      subroutine chord_tria_find(ntrii,nrknot,necke,xtrian,ytrian,
     1 xp,yp,zp,found,itpos)

      use b2mod_types
      use b2mod_indirect
      implicit none

!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag

      integer ntrii,nrknot
      integer necke(3,nrknot)
      real (kind=R8) ::
     *   xtrian(ntrii), ytrian(ntrii)
      real (kind=R8) ::
     *   xp,yp,zp
      logical found
      integer itpos

      logical inside_triangle
      integer i,icnt1,icnt2
      real (kind=R8) ::
     *   x,y,z,r,rmin,rmax,zmin,zmax

      logical first
      integer ir,iz,irpos,izpos
      real (kind=R8) ::
     *   rpos,zpos
      integer lngcov,lngind
      parameter(lngcov=200,lngind=500 000)
      integer covered(lngcov,lngcov),
     1 tcov(lngind),icov(lngind),
     2 ipos,jpos

#ifdef EASY
      integer ir1,ir2,iz1,iz2
#else
      integer ip,irmin,irmax,izmin,izmax
      integer next(1:3),
     1 xmin(lngcov),xmax(lngcov),ymin,ymax
      real (kind=R8) ::
     *   r1,r2,z1,z2
      data next/2,3,1/
#endif

      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0
      common /chord_shift/ rshift, zshift, r0, z0, a0, b0, tt0

      external subini, subend, xerrab
      save covered,rmin,rmax,zmin,zmax,first,tcov,icov

      data first/.true./

      irpos(r)=int((r-rmin)/(rmax-rmin)*lngcov*0.999999_R8+1)
      izpos(z)=int((z-zmin)/(zmax-zmin)*lngcov*0.999999_R8+1)
      rpos(ir)=rmin+real(ir-1)*(rmax-rmin)/(real(lngcov)*0.999999_R8)
      zpos(iz)=zmin+real(iz-1)*(zmax-zmin)/(real(lngcov)*0.999999_R8)

      call subini ('chord_find')
      if(first) then
        first=.false.
        rmin=xtrian(1)
        rmax=xtrian(1)
        zmin=ytrian(1)
        zmax=ytrian(1)
        do i=2,nrknot
          rmin=min(rmin,xtrian(i))
          rmax=max(rmax,xtrian(i))
          zmin=min(zmin,ytrian(i))
          zmax=max(zmax,ytrian(i))
        enddo
        write(*,*) 'rmin,rmax,zmin,zmax ',rmin,rmax,zmin,zmax

        do iz=1,lngcov
          do ir=1,lngcov
            covered(ir,iz)=0
          enddo
        enddo
        ipos=0
        icnt2=0

        do i=1,ntrii
#ifdef EASY
          ir1=irpos(min(xtrian(necke(1,i)),xtrian(necke(2,i)),
     .                  xtrian(necke(3,i))))
          ir2=irpos(max(xtrian(necke(1,i)),xtrian(necke(2,i)),
     .                  xtrian(necke(3,i))))
          iz1=izpos(min(ytrian(necke(1,i)),ytrian(necke(2,i)),
     .                  ytrian(necke(3,i))))
          iz2=izpos(max(ytrian(necke(1,i)),ytrian(necke(2,i)),
     .                  ytrian(necke(3,i))))
          do iz=iz1,iz2
            do ir=ir1,ir2
              ipos=ipos+1
              if(ipos.gt.lngind) then
                write(*,*) i,ipos
                call xerrab ('increase size of {ti}cov')
              endif
              tcov(ipos)=i
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
          do ip=1,3
            r1=xtrian(necke(ip,i))
            z1=ytrian(necke(ip,i))
            r2=xtrian(necke(next(ip),i))
            z2=ytrian(necke(next(ip),i))
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
                write(*,*) 'i ',i
                call xerrab ('increase lngind')
              endif
              tcov(ipos)=i
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
        i=1
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
        if(inside_triangle(r,z,
     1   xtrian(necke(1,i)),ytrian(necke(1,i)),
     2   xtrian(necke(2,i)),ytrian(necke(2,i)),
     3   xtrian(necke(3,i)),ytrian(necke(3,i)))) goto 2100
      endif
c next try the chain
      do while(jpos.ne.0)
!DIAG   chain=chain+1      ! diag
        i=tcov(jpos)
        jpos=icov(jpos)
        if(inside_triangle(r,z,
     1   xtrian(necke(1,i)),ytrian(necke(1,i)),
     2   xtrian(necke(2,i)),ytrian(necke(2,i)),
     3   xtrian(necke(3,i)),ytrian(necke(3,i)))) goto 2100
      enddo
      call subend ()
      return

 2100 continue
      itpos=i
      found=.true.
      call subend ()
      return
      end subroutine chord_tria_find
c
!DIAG      subroutine chord_diagnostics      ! diag
!DIAG      integer rectangles, triangles, chain      ! diag
!DIAG      common /chord_diag/ rectangles, triangles, chain      ! diag
!DIAG      write(*,*) chain,rectangles,triangles      ! diag
!DIAG      end subroutine chord_diagnostics      ! diag

