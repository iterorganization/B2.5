*-----------------------------------------------------------------------
*.specification

      subroutine gradc(nx, ny, hx, hy1, qz, fun, dfun)
      use b2mod_types
      implicit none
*   ..input arguments
      integer nx,ny
      real (kind=R8) ::
     *   hx(-1:nx,-1:ny), hy1(-1:nx,-1:ny),
     &   qz(-1:nx,-1:ny,0:1), fun(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     *   dfun(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*     GRADC computes gradients of a cell centered quantity fun in cell
*     centers.
*     dfun(:,:,0) will hold the gradient in the poloidal direction
*     dfun(:,:,1) will hold the gradient in the radial direction
*
*     The routine is a simple interface to GRADC_P and GRADC_R, that 
*     compute the actual poloidal and radial gradients
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables

*   ..procedures
      external gradc_p, gradc_r

*-----------------------------------------------------------------------
*.computation      

*   ..compute poloidal and radial gradients
      call gradc_p(nx, ny, hx, fun, dfun(-1,-1,0))
      call gradc_r(nx, ny, hx, hy1, qz, fun, dfun(-1,-1,1))

*   ..return
      return

*-----------------------------------------------------------------------
*.end gradc

      end subroutine gradc
