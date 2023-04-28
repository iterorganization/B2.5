      module b2mod_gnmap3
      use b2mod_types
      implicit none

      integer, save :: nxmap, nymap
      integer nmax2d
      parameter (nmax2d=20000)
      real (kind=R8) ::
     *   xmap1(0:nmax2d-1), ymap1(0:nmax2d-1)
      save xmap1, ymap1
      data nxmap/-1/, nymap/-1/

      contains

      subroutine gnmap3 (xinp, yinp, xotp, yotp)
      implicit none
      real (kind=R4) :: xinp, yinp, xotp, yotp
*     ------------------------------------------------------------------
*     GNMAP3 (together with the subroutine GNDFM3) is a service
*     routine for use with the CONPACK routines in the NCARG graphics
*     package. It defines a curvilinear, logically rectangular map for
*     CONPACK. The map is specified by a call to the subroutine gndfm3.
*     ------------------------------------------------------------------
      integer kx, ky
      intrinsic nint, min, max
      external xertst
*     ------------------------------------------------------------------
*!!   Note: nmax2d must be set to a reasonable upper bound for the
*     total number of mesh points.
*     ------------------------------------------------------------------
*     This routine defnes a piecewise bilinear mapping (xinp,yinp)->
*     (xotp,yotp).
*     If gndfm3 is never called, then nxmap and nymap will both have
*     the value -1, and gnmap3 will provide the identity mapping:
*     xotp=xinp and yotp=yinp.
*     If gndfm3 has been called, then the arrays xmap and ymap contain
*     data for linear interpolation. There are nxmap intervals in the
*     x-direction and nymap intervals in the y-direction (the number
*     of data points is (nxmap+1)*(nymap+1)). In that case the input
*     arguments should normally be in the range 0.le.xinp.le.nxmap
*     and 0.le.yinp.le.nymap; if they are out of range then linear
*     extrapolation is used.
*     ------------------------------------------------------------------
      call xertst (nxmap.ne.0.and.nymap.ne.0,
     &  'gnmap3--special case nxmap.eq.0 or nymap.eq.0 not allowed')
      if (nxmap.eq.-1.or.nymap.eq.-1) then
       xotp = xinp
       yotp = yinp
      else
       kx = max(0,min(nxmap-1,nint(xinp-0.5_R8)))
       ky = max(0,min(nymap-1,nint(yinp-0.5_R8)))
       xotp = (kx+1-xinp)*(ky+1-yinp)*real(xmap1(ind1d(kx,ky)),R4)+
     &   (xinp-kx)*(ky+1-yinp)*real(xmap1(ind1d(kx+1,ky)),R4)+
     &   (kx+1-xinp)*(yinp-ky)*real(xmap1(ind1d(kx,ky+1)),R4)+
     &   (xinp-kx)*(yinp-ky)*real(xmap1(ind1d(kx+1,ky+1)),R4)
       yotp = (kx+1-xinp)*(ky+1-yinp)*real(ymap1(ind1d(kx,ky)),R4)+
     &   (xinp-kx)*(ky+1-yinp)*real(ymap1(ind1d(kx+1,ky)),R4)+
     &   (kx+1-xinp)*(yinp-ky)*real(ymap1(ind1d(kx,ky+1)),R4)+
     &   (xinp-kx)*(yinp-ky)*real(ymap1(ind1d(kx+1,ky+1)),R4)
      endif
      return
      end subroutine gnmap3
*     ------------------------------------------------------------------
      subroutine gndfm3 (nxmpd, nympd, nxmp, nymp, xmp, ymp)
      implicit none
      integer nxmpd, nympd, nxmp, nymp
      real (kind=R8) ::
     *   xmp(0:nxmpd,0:nympd), ymp(0:nxmpd,0:nympd)
      integer kx, ky
      call xertst (1.le.nxmp.and.1.le.nymp,
     &  'gndfm3--faulty argument nxmp, nymp')
      call xertst ((nxmp+1)*(nymp+1).le.nmax2d,
     &  'gndfm3--insufficient internal memory')
      nxmap = nxmp
      nymap = nymp
      do ky = nympd-nymp, nympd
       do kx = nxmpd-nxmp, nxmpd
        xmap1(ind1d(kx,ky)) = xmp(kx,ky)
        ymap1(ind1d(kx,ky)) = ymp(kx,ky)
       enddo
      enddo
      return
      end subroutine gndfm3
*     ------------------------------------------------------------------
      integer function ind1d(kx,ky)
      implicit none
      integer, intent(in) :: kx, ky
      ind1d = kx+ky*(nxmap+1)
      return
      end function ind1d

      end module b2mod_gnmap3
