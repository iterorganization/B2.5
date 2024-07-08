*-----------------------------------------------------------------------
*.specification

      subroutine b2trht (nnucl, iCv, pit, zetap,
     & zeffne, ti, Ainv, hci0_al)
      use b2mod_types
      use b2mod_subsys
      use b2mod_constants     
      use b2mod_zhfrtf
     , , only: nal
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl, iCv
      real (kind=R8) ::
     &  pit, zetap, zeffne, ti
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1)
*   ..output arguments 
      real (kind=R8) :: 
     &  hci0_al(0:nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRHT computes the ion heat conductivity coefficient.
*
*  3. description 
*
*     This routine computes:
*     -hci0_al: the ion heat conductivity coefficient 
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: j

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trht')
*   ..main
*   ..sums using inverted matrix 
      do j = 0, nnucl-1
        hci0_al(j)=-2.5_R8*pit**2*(zetap/zeffne**2)*ti/mp*
     &   nal(iCv,j)*sum(Ainv(j,0:nnucl-1)*
     &   nal(iCv,0:nnucl-1))
      enddo

      call subend ()
      return

      end subroutine b2trht


