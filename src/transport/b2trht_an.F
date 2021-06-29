*-----------------------------------------------------------------------
*.specification

      subroutine b2trht_an(nnucl,ix,iy,pit,zetap,
     & zeffne,ti,c_hTA_an,nz2n,rmd2,
     & hci0)
      use b2mod_types
      use b2mod_subsys
      use b2mod_constants     
      use b2mod_zhfrtf, only: nal
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl,ix,iy
      real (kind=R8) ::
     & pit,zetap,zeffne,ti
      real (kind=R8) :: c_hTA_an(0:nnucl-1),nz2n(0:nnucl-1),
     & rmd2(0:nnucl-1)

*   ..output arguments 
      real (kind=R8) :: 
     & hci0
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
*     -hci0: the ion heat conductivity coefficient 
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      real (kind=R8) :: temp

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trht_an')
*   ..main
*   ..sums using analitical coefficient
      temp=sum(c_hTA_an(0:nnucl-1)/rmd2(0:nnucl-1)*
     & nal(ix,iy,0:nnucl-1)**2/nz2n(0:nnucl-1)**2)

      hci0=pit**2*(zetap/zeffne**2)*ti/mp*temp

      call subend ()
      return


      end subroutine b2trht_an


