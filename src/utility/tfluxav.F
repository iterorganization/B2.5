      function tfluxav(quant,bb,hx,vol,nx,ny,iy)

*---------------------------------------------------------------------
*.documentation
*
*   1. purpose
*
*      FLUXAV computes the flux surface average of a quantity
*
*   2. parameters
*
*      quant - array(-1:nx,-1:ny) quantity to be averaged
*      bb    - array(-1:nx,-1:ny,0:3) magnetic field at the center
*              of the (ix,iy) cell
*      hx    - array(-1:nx,-1:ny) x-diameter of the (ix,iy) cell
*      vol   - array(-1:nx,-1:ny) volume of the (ix,iy) cell
*      nx,ny - dimensions
*      iy    - y-label of the flux surface
*
*   3. description
*
*      This function computes and returns the flux surface average
*      of the quantity quant on the flux surface with the y-coordinate
*      iy+1. It assumes that flux surfaces be parameterized by ix=-1,nx
*      for a constant iy and be closed, i.e., inside the core or the
*      the pedestal region. The integration is performed through the
*      Trapezian rule.
*----------------------------------------------------------------------
*.declarations
      use b2mod_types
      use b2mod_indirect

      real (kind=R8) precision, tfluxav
      parameter (precision=1e-8)
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *   quant(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny),
     *   hx(-1:nx,-1:ny)
*   ..local variables
      real (kind=R8) :: 
     *   qsum, qterm, wsum, bface, hxface, hxpface, hxnface
      integer nx, ny, iy, ix
*   ..procedures
      external xertst

*.flux surface averaging
      call xertst (-1.le.iy.and.iy.le.ny, 'faulty argument iy')
      qsum=0
      wsum=0
      do ix=-1,nx
       if (mod(region(ix,iy,0),4).eq.1) then
        if(iy.eq.ny) then
*   ..extrapolation on boundary
          qterm=quant(ix,ny)-vol(ix,ny)*(-quant(ix,ny)+
     &          quant(bottomix(ix,ny),bottomiy(ix,ny)))/
     &          (vol(ix,ny)+vol(bottomix(ix,ny),bottomiy(ix,ny)))
          bface=bb(ix,ny,0)-vol(ix,ny)*(-bb(ix,ny,0)+
     &          bb(bottomix(ix,ny),bottomiy(ix,ny),0))/
     &          (vol(ix,ny)+vol(bottomix(ix,ny),bottomiy(ix,ny)))
          hxface=hx(ix,ny)-vol(ix,ny)*(-hx(ix,ny)+
     &          hx(bottomix(ix,ny),bottomiy(ix,ny)))/
     &          (vol(ix,ny)+vol(bottomix(ix,ny),bottomiy(ix,ny)))
          hxpface=hx(leftix(ix,ny),ny)-vol(leftix(ix,ny),ny)*
     &          (-hx(leftix(ix,ny),ny)+
     &          hx(bottomix(leftix(ix,ny),ny),
     &          bottomiy(leftix(ix,ny),ny)))/
     &          (vol(leftix(ix,ny),ny)+
     &          vol(bottomix(leftix(ix,ny),ny),
     &          bottomiy(leftix(ix,ny),ny)))
          hxnface=hx(rightix(ix,ny),ny)-vol(rightix(ix,ny),ny)*
     &          (-hx(rightix(ix,ny),ny)+
     &          hx(bottomix(rightix(ix,ny),ny),
     &          bottomiy(rightix(ix,ny),ny)))/
     &          (vol(rightix(ix,ny),ny)+
     &          vol(bottomix(rightix(ix,ny),ny),
     &          bottomiy(rightix(ix,ny),ny)))
        else
          qterm=(vol(ix,iy)*quant(topix(ix,iy),topiy(ix,iy))+
     &           vol(topix(ix,iy),topiy(ix,iy))*quant(ix,iy))/
     &          (vol(topix(ix,iy),topiy(ix,iy))+vol(ix,iy))
          bface=(vol(ix,iy)*bb(topix(ix,iy),topiy(ix,iy),0)+
     &           vol(topix(ix,iy),topiy(ix,iy))*bb(ix,iy,0))/
     &          (vol(topix(ix,iy),topiy(ix,iy))+vol(ix,iy))
          hxface=(vol(ix,iy)*hx(topix(ix,iy),topiy(ix,iy))+
     &           vol(topix(ix,iy),topiy(ix,iy))*hx(ix,iy))/
     &          (vol(topix(ix,iy),topiy(ix,iy))+vol(ix,iy))
          hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &            hx(topix(leftix(ix,iy),leftiy(ix,iy)),
     &            topiy(leftix(ix,iy),leftiy(ix,iy)))+
     &            vol(topix(leftix(ix,iy),leftiy(ix,iy)),
     &            topiy(leftix(ix,iy),leftiy(ix,iy)))*
     &            hx(leftix(ix,iy),leftiy(ix,iy)))/
     &            (vol(topix(leftix(ix,iy),leftiy(ix,iy)),
     &            topiy(leftix(ix,iy),leftiy(ix,iy)))+
     &            vol(leftix(ix,iy),leftiy(ix,iy)))
          hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &            hx(topix(rightix(ix,iy),rightiy(ix,iy)),
     &            topiy(rightix(ix,iy),rightiy(ix,iy)))+
     &            vol(topix(rightix(ix,iy),rightiy(ix,iy)),
     &            topiy(rightix(ix,iy),rightiy(ix,iy)))*
     &            hx(rightix(ix,iy),rightiy(ix,iy)))/
     &            (vol(topix(rightix(ix,iy),rightiy(ix,iy)),
     &            topiy(rightix(ix,iy),rightiy(ix,iy)))+
     &            vol(rightix(ix,iy),rightiy(ix,iy)))
        endif
        qsum=qsum+0.5_R8*qterm/bface*((hxpface+hxnface)/2._R8+
     &                                hxface)
        wsum=wsum+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/bface
       endif
      enddo
      if(abs(wsum).lt.precision) then		! psauter
        tfluxav=0.0_R8
      else
        tfluxav=qsum/wsum
      endif

*   ..return
      return
      end
