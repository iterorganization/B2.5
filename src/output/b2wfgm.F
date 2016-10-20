*-----------------------------------------------------------------------
*.specification

      subroutine b2wfgm (nprt, idout, nx, ny)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_indirect
      use b2mod_b2cmgs
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny
      character idout*(*)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFGM writes to a formatted file a summary of the b2 geometry
*     and magnetic field.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to a formatted file reduced information on
*     the b2 geometry and magnetic field. The information includes
*     cell volumes, cell diameters, and cell surface areas, all reduced
*     to a function of one coordinate by summing over the other
*     coordinate.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     idout - character*(*), input.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     crx, cry, fpsi, ffbz - real*8 array, input.
*     bb, vol, hx, hy, qz, qc, gs - real*8 array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, i
      real (kind=R8) ::
     *   xc0, xc1, yc0, yc1, rvx, rvy, rhx, rhy, rsx, rsy,
     &  tmin, tmax, qmin, qmax
*   ..procedures
      real (kind=R8) ::
     *   smin, smax, b2ssum
      logical strmas
      real (kind=R8) ::
     *   crxc, cryc
      external subini, subend, xertst, smin, smax, b2ssum, strmas
      external b2xvgm
      intrinsic max
      crxc(ix,iy) = (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+
     &  crx(ix,iy,3))/4.0_R8
      cryc(ix,iy) = (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+
     &  cry(ix,iy,3))/4.0_R8

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfgm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nprt, nx, ny
      call xertst (1.le.nprt, 'faulty argument nprt')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test geometry
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, gs, pbs)

* ..write the metric data
*   ..write heading
      write (nprt,101) nx, ny
*   ..write pgmm segment?
      if (strmas(idout,'pgmm')) then
*    ..write range of crx, cry, bb
       write (nprt,'(/2x,a6)') '+pgmm'
       tmin = smin(n2*4,crx,1)
       tmax = smax(n2*4,crx,1)
       write (nprt,'(/2x,a,1p,2e11.3)') 'Range of crx: ', tmin, tmax
       tmin = smin(n2*4,cry,1)
       tmax = smax(n2*4,cry,1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of cry: ', tmin, tmax
       tmin = smin(n2*4,fpsi,1)
       tmax = smax(n2*4,fpsi,1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of fpsi:', tmin, tmax
       tmin = smin(n2*4,ffbz,1)
       tmax = smax(n2*4,ffbz,1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of ffbz:', tmin, tmax
       tmin = smin(n2,bb(-1,-1,0),1)
       tmax = smax(n2,bb(-1,-1,0),1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of bx:  ', tmin, tmax
       tmin = smin(n2,bb(-1,-1,1),1)
       tmax = smax(n2,bb(-1,-1,1),1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of by:  ', tmin, tmax
       tmin = smin(n2,bb(-1,-1,2),1)
       tmax = smax(n2,bb(-1,-1,2),1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of bz:  ', tmin, tmax
       tmin = smin(n2,bb(-1,-1,3),1)
       tmax = smax(n2,bb(-1,-1,3),1)
       write (nprt,'(2x,a,1p,2e11.3)') 'Range of bb:  ', tmin, tmax
      endif
*   ..write pgcc segment?
      if (strmas(idout,'pgcc')) then
*    ..write 5*5 corners
       write (nprt,'(/2x,a6)') '+pgcc'
       write (nprt,102)
       do iy = 0, ny, max(ny/4,1)
        do ix = 0, nx, max(nx/4,1)
c---DPC problematic
         if(leftix(ix,iy).ne.-2) then
          if(bottomiy(ix,iy).ne.-2) then
           if(bottomiy(leftix(ix,iy),leftiy(ix,iy)).ne.-2) then
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        crx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)),3), 
     &        cry(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)),3),
     &        crx(bottomix(ix,iy),bottomiy(ix,iy),2),
     1        cry(bottomix(ix,iy),bottomiy(ix,iy),2),
     &        crx(leftix(ix,iy),leftiy(ix,iy),1),
     1        cry(leftix(ix,iy),leftiy(ix,iy),1),
     &        crx(ix,iy,0), cry(ix,iy,0)
           else
            if(leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.-2) then
             write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        crx(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     &            leftiy(bottomix(ix,iy),bottomiy(ix,iy)),3), 
     &        cry(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     &            leftiy(bottomix(ix,iy),bottomiy(ix,iy)),3),
     &        crx(bottomix(ix,iy),bottomiy(ix,iy),2),
     1        cry(bottomix(ix,iy),bottomiy(ix,iy),2),
     &        crx(leftix(ix,iy),leftiy(ix,iy),1),
     1        cry(leftix(ix,iy),leftiy(ix,iy),1),
     &        crx(ix,iy,0), cry(ix,iy,0)
            else
             write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        -999.0_R8,-999.0_R8,
     &        crx(bottomix(ix,iy),bottomiy(ix,iy),2),
     1        cry(bottomix(ix,iy),bottomiy(ix,iy),2),
     &        crx(leftix(ix,iy),leftiy(ix,iy),1),
     1        cry(leftix(ix,iy),leftiy(ix,iy),1),
     &        crx(ix,iy,0), cry(ix,iy,0)
            endif
           endif
          else
           if(bottomiy(leftix(ix,iy),leftiy(ix,iy)).ne.-2) then
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        crx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)),3), 
     &        cry(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)),3),
     &        -999.0_R8,-999.0_R8,
     &        crx(leftix(ix,iy),leftiy(ix,iy),1),
     1        cry(leftix(ix,iy),leftiy(ix,iy),1),
     &        crx(ix,iy,0), cry(ix,iy,0)
           else
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        -999.0_R8,-999.0_R8,-999.0_R8,-999.0_R8,
     &        crx(leftix(ix,iy),leftiy(ix,iy),1),
     1        cry(leftix(ix,iy),leftiy(ix,iy),1),
     &        crx(ix,iy,0), cry(ix,iy,0)
           endif
          endif
         else
          if(bottomiy(ix,iy).ne.-2) then
           if(leftix(bottomix(ix,iy),bottomiy(ix,iy)).ne.-2) then
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        crx(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     &            leftiy(bottomix(ix,iy),bottomiy(ix,iy)),3), 
     &        cry(leftix(bottomix(ix,iy),bottomiy(ix,iy)),
     &            leftiy(bottomix(ix,iy),bottomiy(ix,iy)),3),
     &        crx(bottomix(ix,iy),bottomiy(ix,iy),2),
     1        cry(bottomix(ix,iy),bottomiy(ix,iy),2),
     &        -999.0_R8,-999.0_R8,
     &        crx(ix,iy,0), cry(ix,iy,0)
           else
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        -999.0_R8,-999.0_R8,
     &        crx(bottomix(ix,iy),bottomiy(ix,iy),2),
     1        cry(bottomix(ix,iy),bottomiy(ix,iy),2),
     &        -999.0_R8,-999.0_R8,
     &        crx(ix,iy,0), cry(ix,iy,0)
           endif
          else
            write (nprt,'(1p,2x,i4,1x,i4,4(2x,2e11.3))') ix, iy,
     &        -999.0_R8,-999.0_R8,-999.0_R8,-999.0_R8,
     &        -999.0_R8,-999.0_R8,
     &        crx(ix,iy,0), cry(ix,iy,0)
          endif
         endif
        enddo
       enddo
      endif
*   ..write pgnl segment?
      if (strmas(idout,'pgnl')) then
       write (nprt,'(/2x,a6)') '+pgnl'
       if (nlreg.eq.0) then
        write (nprt,'(/2x,a)') 'There are no singular regions.'
       else
        write (nprt,'(/2x,a)')
     &    'There are singular regions; the parameters are:'
        write (nprt,'(4x,a)') ' xlo xhi ylo yhi loc'
        do i = 0, nlreg-1
          write (nprt,'(4x,5i4)')
     &      nlxlo(i), nlxhi(i), nlylo(i), nlyhi(i), nlloc(i)
        enddo
       endif
      endif
*   ..write pgmx segment?
      if (strmas(idout,'pgmx')) then
*    ..write functions of x
       write (nprt,'(/2x,a6)') '+pgmx'
       write (nprt,103)
       xc0 = -hx(-1,-1)
       xc1 = -hx(-1,ny)
       do ix = -1, nx
*     ..update xc0, xc1
        xc0 = xc0+hx(ix,-1)/2.0_R8
        xc1 = xc1+hx(ix,ny)/2.0_R8
*     ..compute 'vol(ix)'
        rvx = b2ssum(ny+2,vol(ix,-1),nx+2)
*     ..compute 'hy(ix)'
        rhy = b2ssum(ny+2,hy(ix,-1),nx+2)
*     ..compute 'qmin(ix)'
        qmin = smin(ny+2,qz(ix,-1,0),nx+2)
*     ..compute 'qmax(ix)'
        qmax = smax(ny+2,qz(ix,-1,0),nx+2)
*     ..compute 'sx(ix)'
        if (leftix(ix,-1).eq.-2) then
         rsx = b2ssum(ny+2,gs(rightix(ix,-1),rightiy(ix,-1),0),nx+2)
        else if (rightix(ix,-1).eq.nx+1) then
         rsx = b2ssum(ny+2,gs(ix,-1,0),nx+2)
        else
         rsx = (b2ssum(ny+2,gs(ix,-1,0),nx+2)+
     &          b2ssum(ny+2,gs(rightix(ix,-1),
     &                         rightiy(ix,-1),0),nx+2))/2.0_R8
        endif
*     ..write line
        write (nprt,'(2x,i4,1x,1p,7e11.3)')
     &    ix, xc0, xc1, rvx, rhy, rsx, qmin, qmax
*     ..update xc0, xc1
        xc0 = xc0+hx(ix,-1)/2.0_R8
        xc1 = xc1+hx(ix,ny)/2.0_R8
       enddo
      endif
*   ..write pgmy segment?
      if (strmas(idout,'pgmy')) then
*    ..write functions of y
       write (nprt,'(/2x,a6)') '+pgmy'
       write (nprt,104)
       yc0 = -hy(-1,-1)
       yc1 = -hy(nx,-1)
       do iy = -1, ny
*     ..update yc0, yc1
        yc0 = yc0+hy(-1,iy)/2.0_R8
        yc1 = yc1+hy(nx,iy)/2.0_R8
*     ..compute 'vol(iy)'
        rvy = b2ssum(nx+2,vol(-1,iy),1)
*     ..compute 'hx(iy)'
        rhx = b2ssum(nx+2,hx(-1,iy),1)
*     ..compute 'qmin(iy)'
        qmin = smin(nx+2,qz(-1,iy,0),1)
*     ..compute 'qmax(iy)'
        qmax = smax(nx+2,qz(-1,iy,0),1)
*     ..compute 'sy(iy)'
        if (bottomiy(-1,iy).eq.-2) then
         rsy = b2ssum(nx+2,gs(topix(-1,iy),topiy(-1,iy),1),1)
        else if (topiy(-1,iy).eq.ny+1) then
         rsy = b2ssum(nx+2,gs(-1,iy,1),1)
        else
         rsy = (b2ssum(nx+2,gs(-1,iy,1),1)+
     &          b2ssum(nx+2,gs(topix(-1,iy),topiy(-1,iy),1),1))/2.0_R8
        endif
*     ..write line
        write (nprt,'(2x,i4,1x,1p,7e11.3)')
     &    iy, yc0, yc1, rvy, rhx, rsy, qmin, qmax
*     ..update yc0, yc1
        yc0 = yc0+hy(-1,iy)/2.0_R8
        yc1 = yc1+hy(nx,iy)/2.0_R8
       enddo
      endif
*   ..write pgso segment?
      if (strmas(idout,'pgso')) then
       write (nprt,'(/2x,a6)') '+pgso'
*    ..write geometry on south boundary
       write (nprt,111)
       do ix = -1, nx
        write (nprt,'(2x,i4,1x,1p,2e11.3,2x,4e11.3,2(2x,2e11.3))')
     &    ix, crxc(ix,-1), cryc(ix,-1), bb(ix,-1,0), bb(ix,-1,1),
     &    bb(ix,-1,2), bb(ix,-1,3), hx(ix,-1), gs(ix,0,1),
     &    qz(ix,-1,0), qz(ix,-1,1)
       enddo
      endif
*   ..write pgno segment?
      if (strmas(idout,'pgno')) then
       write (nprt,'(/2x,a6)') '+pgno'
*    ..write geometry on north boundary
       write (nprt,112)
       do ix = -1, nx
        write (nprt,'(2x,i4,1x,1p,2e11.3,2x,4e11.3,2(2x,2e11.3))')
     &    ix, crxc(ix,ny), cryc(ix,ny), bb(ix,ny,0), bb(ix,ny,1),
     &    bb(ix,ny,2), bb(ix,ny,3), hx(ix,ny), gs(ix,ny,1),
     &    qz(ix,ny,0), qz(ix,ny,1)
       enddo
      endif
*   ..write pgwe segment?
      if (strmas(idout,'pgwe')) then
       write (nprt,'(/2x,a6)') '+pgwe'
*    ..write geometry on west boundary
       write (nprt,113)
       do iy = -1, ny
        write (nprt,'(2x,i4,1x,1p,2e11.3,2x,4e11.3,2(2x,2e11.3))')
     &    iy, crxc(-1,iy), cryc(-1,iy), bb(-1,iy,0), bb(-1,iy,1),
     &    bb(-1,iy,2), bb(-1,iy,3), hy(-1,iy), gs(0,iy,0),
     &    qz(-1,iy,0), qz(-1,iy,1)
       enddo
      endif
*   ..write pgea segment?
      if (strmas(idout,'pgea')) then
       write (nprt,'(/2x,a6)') '+pgea'
*    ..write geometry on east boundary
       write (nprt,114)
       do iy = -1, ny
        write (nprt,'(2x,i4,1x,1p,2e11.3,2x,4e11.3,2(2x,2e11.3))')
     &    iy, crxc(nx,iy), cryc(nx,iy), bb(nx,iy,0), bb(nx,iy,1),
     &    bb(nx,iy,2), bb(nx,iy,3), hy(nx,iy), gs(nx,iy,0),
     &    qz(nx,iy,0), qz(nx,iy,1)
       enddo
      endif
*   ..write blank line
      write (nprt,'()')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  101 format (/'  The grid contains nx*ny interior cells, ',
     &  'where nx =',i4,' and ny =',i4,'.'/
     &  '  With boundary cells included the size is (nx+2)*(ny+2).')
  102 format (/'  Locations of 5*5 corners:'/
     &  t5,'ix   iy',t17,'crx(--)    cry(--)',6x,'crx(+-)    cry(+-)',
     &  6x,'crx(-+)    cry(-+)',6x,'crx(++)    cry(++)')
  103 format (/'  The grid metric parameters reduced to a function',
     &  ' of ix only are:'/
     &  t5,'ix',t12,'xc0',t23,'xc1',t33,'vol(ix)',t44,'hy(ix)',
     &  t55,'sx(ix)',t65,'qmin(ix)',t76,'qmax(ix)')
  104 format (/'  The grid metric parameters reduced to a function',
     &  ' of iy only are:'/
     &  t5,'iy',t12,'yc0',t23,'yc1',t33,'vol(iy)',t44,'hx(iy)',
     &  t55,'sy(iy)',t65,'qmin(iy)',t76,'qmax(iy)')
  111 format (/'  Geometry on the south boundary:'/
     &  t4,'ix',t12,'crx',t23,'cry',t37,'bx',t48,'by',t59,'bz',
     &  t70,'bb',t83,'hx',t94,'sy',t106,'qz0',t117,'qz1')
  112 format (/'  Geometry on the north boundary:'/
     &  t4,'ix',t12,'crx',t23,'cry',t37,'bx',t48,'by',t59,'bz',
     &  t70,'bb',t83,'hx',t94,'sy',t106,'qz0',t117,'qz1')
  113 format (/'  Geometry on the west boundary:'/
     &  t4,'iy',t12,'crx',t23,'cry',t37,'bx',t48,'by',t59,'bz',
     &  t70,'bb',t83,'hy',t94,'sx',t106,'qz0',t117,'qz1')
  114 format (/'  Geometry on the east boundary:'/
     &  t4,'iy',t12,'crx',t23,'cry',t37,'bx',t48,'by',t59,'bz',
     &  t70,'bb',t83,'hy',t94,'sx',t106,'qz0',t117,'qz1')

*-----------------------------------------------------------------------
*.end b2wfgm

      end subroutine b2wfgm
