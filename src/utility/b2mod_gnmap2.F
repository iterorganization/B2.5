      module b2mod_gnmap2
      use b2mod_types
      implicit none

      integer nmax, nxmap, nymap
      save nxmap, nymap
      parameter (nmax=300)
      real (kind=R8) ::
     *   xmap(0:nmax), ymap(0:nmax)
      save xmap, ymap
      data nxmap/-1/, nymap/-1/

      contains

      subroutine gnmap2 (xinp, yinp, xotp, yotp)
      implicit none
      real (kind=R4) :: xinp, yinp, xotp, yotp
*     ------------------------------------------------------------------
*     GNMAP2 (together with the subroutine GNDFM2) is a service
*     routine for use with the CONPACK routines in the NCARG graphics
*     package. It defines a rectangular map for CONPACK. The map is
*     specified by a call to the subroutine gndfm2.
*     ------------------------------------------------------------------
      integer kx, ky
      intrinsic nint, min, max
      external xertst
*     ------------------------------------------------------------------
*!!   Note: nmax must be set to a reasonable upper bound for the number
*     of mesh points in any coordinate direction.
*     ------------------------------------------------------------------
*     This routine defines a piecewise linear mapping xinp->xotp,
*     yinp->yotp.
*     If gndfm2 is never called, then nxmap and nymap will both have
*     the value -1, and gnmap2 will provide the identity mapping:
*     xotp=xinp and yotp=yinp.
*     If gndfm2 has been called, then the arrays xmap and ymap contain
*     data for linear interpolation. There are nxmap intervals in the
*     x-direction and nymap intervals in the y-direction (the number
*     of data points is nxmap+1 and nymap+1). In that case the input
*     arguments should normally be in the range 0.le.xinp.le.nxmap
*     and 0.le.yinp.le.nymap; if they are out of range then linear
*     extrapolation is used.
*     ------------------------------------------------------------------
      if (nxmap.eq.-1) then
       xotp = xinp
      else if (nxmap.eq.0) then
       xotp = real(xmap(0),R4)
      else
       kx = max(0,min(nxmap-1,nint(xinp-0.5_R8)))
       xotp = (kx+1-xinp)*real(xmap(kx),R4)+
     &        (xinp-kx)*real(xmap(kx+1),R4)
      endif
      if (nymap.eq.-1) then
       yotp = yinp
      else if (nymap.eq.0) then
       yotp = real(ymap(0),R4)
      else
       ky = max(0,min(nymap-1,nint(yinp-0.5_R8)))
       yotp = (ky+1-yinp)*real(ymap(ky),R4)+
     &        (yinp-ky)*real(ymap(ky+1),R4)
      endif
      return
      end subroutine gnmap2
*     ------------------------------------------------------------------
      subroutine gndfm2 (xmp, nxmp, ymp, nymp)
      implicit none
      integer nxmp, nymp
      real (kind=R8) ::
     *   xmp(0:nxmp), ymp(0:nymp)
      integer kx, ky
      call xertst (-1.le.nxmp.and.-1.le.nymp,
     &  'gndfm2--faulty argument nxmp, nymp')
      call xertst (nxmp.le.nmax.and.nymp.le.nmax,
     &  'gndfm2--insufficient internal memory')
      nxmap = nxmp
      do kx = 0, nxmp
       xmap(kx) = xmp(kx)
      enddo
      nymap = nymp
      do ky = 0, nymp
       ymap(ky) = ymp(ky)
      enddo
      return
      end subroutine gndfm2
*     ------------------------------------------------------------------
      end module b2mod_gnmap2
