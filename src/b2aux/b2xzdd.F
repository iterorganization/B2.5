      subroutine b2xzdd (nx, ny, nd, src)
      use b2mod_types
      use b2mod_indirect
      use b2mod_boundary_namelist
      implicit none
      integer nx, ny, nd
      real (kind=R8) ::
     *   src(-1:nx, -1:ny, 0:nd)
*     ------------------------------------------------------------------
*      B2XZDD zeroes out the array elements corresponding to dead cells
*      and core boundary cells.
*     ------------------------------------------------------------------
      integer ib, id, ix, iy
      integer ncall, zero_dead_and_core
      save ncall, zero_dead_and_core
*   ..procedures
      external ipgeti
*   ..initialisation
      data ncall/0/, zero_dead_and_core/1/
*     ------------------------------------------------------------------

      if (ncall.eq.0) then
        call ipgeti('b2xzdd_zero_dead_and_core', zero_dead_and_core)
      endif
      if (zero_dead_and_core.ge.1) then
        do iy = -1, ny
          do ix = -1, nx
            if(region(ix,iy,0).eq.0) then
              do id = 0, nd
                src(ix,iy,id) = 0.0_R8
              enddo
            endif
          enddo
        enddo
      endif
      if (zero_dead_and_core.ge.2) then
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if (core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                do id = 0, nd
                  src(bc_list_x(ix,ib),bc_list_y(ix,ib),id) = 0.0_R8
                enddo
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              do id = 0, nd
                src(ix,-1,id) = 0.0_R8
              enddo
            end if
          end do
        endif
      endif
      ncall=ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xzdd
