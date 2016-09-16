      function sfluxav(quant,bb,hx,vol,nx,ny,iy,nlambda,bmax)

*---------------------------------------------------------------------
*.documentation
*
*   1. purpose
*
*      SFLUXAV computes a specific flux surface average of a quantity
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
*      This function computes and returns a specific flux surface average
*      of the quantity quant on the flux surface with the y-coordinate
*      iy. It assumes that flux surfaces be parameterized by ix=-1,nx
*      for a constant iy and be closed, i.e., inside the core or the
*      the pedestal region. The integration is performed through the
*      Trapezian rule.
*----------------------------------------------------------------------
*.declarations
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) precision, sfluxav
      parameter (precision=1e-8)
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *   quant(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny),
     *   hx(-1:nx,-1:ny), nlambda, bmax, quantn(-1:nx,-1:ny)
*   ..local variables
      real (kind=R8) :: 
     *   qsum, qterm, wsum, bface, hxface, hxpface, hxnface
      integer iy, ix
*   ..procedures
      external xertst

*.flux surface averaging
      call xertst (-1.le.iy.and.iy.le.ny, 'faulty argument iy')
      qsum=0
      wsum=0
*   ..working variable quantn
      quantn=quant
      do ix=-1,nx
       if (mod(region(ix,iy,0),4).eq.1) then
        quantn(ix,iy)=sqrt(1._R8-nlambda*bb(ix,iy,3)/bmax)
        if (1._R8-nlambda*bb(topix(ix,iy),topiy(ix,iy),3)/bmax.lt.0) 
     &   then
          quantn(topix(ix,iy),topiy(ix,iy))=0._R8
         else
          quantn(topix(ix,iy),topiy(ix,iy))=real(sqrt(1._R8-nlambda*
     &           bb(topix(ix,iy),topiy(ix,iy),3)/bmax), R8)
        endif
        if (1._R8-nlambda*bb(bottomix(ix,iy),bottomiy(ix,iy),3)/bmax
     &   .lt.0) then
          quantn(bottomix(ix,iy),bottomiy(ix,iy))=0._R8
         else
          quantn(bottomix(ix,iy),bottomiy(ix,iy))=real(
     &          sqrt(1._R8-nlambda*bb(bottomix(ix,iy),
     &          bottomiy(ix,iy),3)/bmax),R8)
        endif

        if(iy.eq.-1) then
*   ..extrapolation on boundary
          qterm=quantn(ix,-1)-vol(ix,-1)*(-quantn(ix,-1)+
     &          quantn(topix(ix,-1),topiy(ix,-1)))/
     &          (vol(ix,-1)+vol(topix(ix,-1),topiy(ix,-1)))
          bface=bb(ix,-1,0)-vol(ix,-1)*(-bb(ix,-1,0)+
     &          bb(topix(ix,-1),topiy(ix,-1),0))/
     &          (vol(ix,-1)+vol(topix(ix,-1),topiy(ix,-1)))
          hxface=hx(ix,-1)-vol(ix,-1)*(-hx(ix,-1)+
     &          hx(topix(ix,-1),topiy(ix,-1)))/
     &          (vol(ix,-1)+vol(topix(ix,-1),topiy(ix,-1)))
          hxpface=hx(leftix(ix,-1),-1)-vol(leftix(ix,-1),-1)*
     &          (-hx(leftix(ix,-1),-1)+
     &          hx(topix(leftix(ix,-1),-1),
     &          topiy(leftix(ix,-1),-1)))/
     &          (vol(leftix(ix,-1),-1)+
     &          vol(topix(leftix(ix,-1),-1),
     &          topiy(leftix(ix,-1),-1)))
          hxnface=hx(rightix(ix,-1),-1)-vol(rightix(ix,-1),-1)*
     &          (-hx(rightix(ix,-1),-1)+
     &          hx(topix(rightix(ix,-1),-1),
     &          topiy(rightix(ix,-1),-1)))/
     &          (vol(rightix(ix,-1),-1)+
     &          vol(topix(rightix(ix,-1),-1),
     &          topiy(rightix(ix,-1),-1)))
        else
          qterm=(vol(ix,iy)*quantn(bottomix(ix,iy),bottomiy(ix,iy))+
     &           vol(bottomix(ix,iy),bottomiy(ix,iy))*quantn(ix,iy))/
     &          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
          bface=(vol(ix,iy)*bb(bottomix(ix,iy),bottomiy(ix,iy),0)+
     &           vol(bottomix(ix,iy),bottomiy(ix,iy))*bb(ix,iy,0))/
     &          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
          hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &           vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &          (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
          hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &            hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &            vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &            hx(leftix(ix,iy),leftiy(ix,iy)))/
     &            (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &            vol(leftix(ix,iy),leftiy(ix,iy)))
          hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &            hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &            vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &            hx(rightix(ix,iy),rightiy(ix,iy)))/
     &            (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &            vol(rightix(ix,iy),rightiy(ix,iy)))
        endif
        qsum=qsum+0.5_R8*qterm/bface*((hxpface+hxnface)/2._R8+
     &                                hxface)
        wsum=wsum+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/bface
       endif
      enddo
      if(abs(wsum).lt.precision) then		! psauter
        sfluxav=0.0_R8
      else
        sfluxav=qsum/wsum
      endif

*   ..return
      return
      end function sfluxav
