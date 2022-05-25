      subroutine b2xgbs (g,m)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none

      type(mapping) :: m
      type(geometry):: g
      

*     ------------------------------------------------------------------
*     B2XGBS computes the product (magnetic field pitch)*(surface area)
*     at cell faces fcPbs. 
*     fcPbs thus contains the parallel contact area.
*     ------------------------------------------------------------------


*   ..preliminaries
      call subini('b2xgbs')
*   ..test input


!      stuctured code of the calculation
!*   ..compute pbs(,,0:1)
!      do iy = -1, ny
!       do ix = -1, nx
!        if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
!        if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) 
!     &    pbs(ix,iy,0) = gs(ix,iy,0)*qac(ix,iy,0)*abs(wbbl(ix,iy,0))/
!     &      wbbl(ix,iy,3)
!        if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) 
!     &    pbs(ix,iy,1) = gs(ix,iy,1)*qac(ix,iy,1)*abs(wbbv(ix,iy,0))/
!     &      wbbv(ix,iy,3)
!       enddo
!      enddo


      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xgbs
