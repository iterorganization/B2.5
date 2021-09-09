*-----------------------------------------------------------------------
*.specification

      subroutine b2trht_an(nnucl,ix,iy,pit,zetap,
     & zeffne,ti,c_hT_an,nz2n,rmd2,
     & hci0_al)
      use b2mod_types
      use b2mod_subsys
      use b2mod_constants     
      use b2mod_zhfrtf, only: nal
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl,ix,iy
      real (kind=R8) ::
     & pit,zetap,zeffne,ti
      real (kind=R8) :: c_hT_an(0:nnucl-1),nz2n(0:nnucl-1),
     & rmd2(0:nnucl-1)

*   ..output arguments 
      real (kind=R8) :: 
     & hci0_al(0:nnucl-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2trtr_an computes the ion heat conductivity coefficient 
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
      call subini ('b2trht_an')
*   ..main
*   ..sums using analitical coefficient
      do j = 0, nnucl-1
        hci0_al(j)=pit**2*(zetap/zeffne**2)*ti/mp*c_hT_an(j)/rmd2(j)*
     & nal(ix,iy,j)**2/nz2n(j)**2
      enddo

      call subend ()
      return


      end subroutine b2trht_an


