      subroutine b2tbfl(nx, ny, lcore, flux, tbfl, div)              !srv 12.02.16 {
      use b2mod_types
      use b2mod_boundary_namelist
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      logical lcore
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: flux(-1:nx,-1:ny,0:1), tbfl(1:nbc), div
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TBFL computes the total flux through the boundaries
*            if lcore is true then all boundaries are included
*            if lcore is false then Core boundary is not included
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: ib, ix, iy
      integer, save :: ncall = 0
!$OMP THREADPRIVATE(ncall)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tbfl')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      div = 0.0_R8
      do ib=1, nbc
        tbfl(ib)=0.0e0_R8
        IF(BCCHAR(IB).EQ.'S') THEN
          DO IX=1, BC_LIST_SIZE(IB)
            if(lcore.or.
     &       mod(region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0),4).ne.1)
     >        tbfl(ib) = tbfl(ib) +
     &         flux(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &              TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1)
          end do
          div = div - tbfl(ib)
        ELSEIF(BCCHAR(IB).EQ.'N') THEN
          DO IX=1, BC_LIST_SIZE(IB)
            if(lcore.or.
     &         region(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),0).ne.5)
     >       tbfl(ib) = tbfl(ib) +
     &        flux(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1)
          end do
          div = div + tbfl(ib)
        ELSEIF(BCCHAR(IB).EQ.'W') THEN
          DO IY=1, BC_LIST_SIZE(IB)
            tbfl(ib) = tbfl(ib) +
     &        flux(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &             rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0)
          end do
          div = div - tbfl(ib)
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          DO IY=1, BC_LIST_SIZE(IB)
            tbfl(ib) = tbfl(ib) +
     &        flux(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0)
          end do
          div = div + tbfl(ib)
        ENDIF
      enddo
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tbfl
      end subroutine b2tbfl                                !srv 12.02.16 }
