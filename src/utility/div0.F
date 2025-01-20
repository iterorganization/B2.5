*-----------------------------------------------------------------------
*.specification

      subroutine div0(nCv, nFc, mpg, flo, dflo)
      use b2mod_types
      use b2us_map
      implicit none
*   ..input arguments
      integer, intent (in) :: nCv, nFc
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   flo(nFc)
*   ..output arguments
      real (kind=R8) ::
     *   dflo(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     DIV0 computes the divergence of a flow field flo which has only a 
*     single flo component (poloidal or radial)
*
*-----------------------------------------------------------------------
*.declarations

      integer iFc

*-----------------------------------------------------------------------
*.computation      

*   ..initialize to zero
      dflo = 0.0_R8

*   ..compute divergence
      do iFc = 1, nFc
        dflo(mpg%fcCv(iFc,1)) = 
     &    dflo(mpg%fcCv(iFc,1)) + flo(iFc)
        dflo(mpg%fcCv(iFc,2)) = 
     &    dflo(mpg%fcCv(iFc,2)) - flo(iFc)
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end div

      end subroutine div0
