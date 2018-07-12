*-----------------------------------------------------------------------
*.specification

      subroutine div(nCv, nFc, mpg, flo, dflo)
      use b2mod_types
      use b2us_map
      implicit none
*   ..input arguments
      integer nCv, nFc
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   flo(nFc,0:1)
*   ..output arguments
      real (kind=R8) ::
     *   dflo(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     DIV computes the divergence of a flow field flo.
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
     &    dflo(mpg%fcCv(iFc,1)) + flo(iFc,0) + flo(iFc,1)
        dflo(mpg%fcCv(iFc,2)) = 
     &    dflo(mpg%fcCv(iFc,2)) - flo(iFc,0) - flo(iFc,1)
      enddo

*   ..return
      return

*-----------------------------------------------------------------------
*.end div

      end subroutine div
