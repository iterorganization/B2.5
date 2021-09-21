*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner_a (nx, ny, ns, vol, gs, bb, na, ua,             !srv 20.09.06 { 08.09.21 {
     &                    facdrift, fchinert_a, fchinert)                !srv 22.06.08
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     & na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1 bb(-1:nx,-1:ny,0:3), facdrift(-1:nx,-1:ny)                        !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &     fchinert_a(-1:nx,-1:ny,0:1,0:ns-1), fchinert(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINER_A computes the X and Y components of the inertial current 
*     for each species and calculates the X and Y components of 
*     the total inertial current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer ncall, ix, iy, ib, is
      integer, save :: fch_inert_style                                   !srv 19.04.13
      integer, save :: no_current = 0
      integer, save :: guard_flows
      real (kind=R8) ::
     &  fch_inert
      save ncall, fch_inert
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_inert/1.0e0_R8/, guard_flows/2/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tiner_a')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tiner_inert', fch_inert)                          !srv 07.07.99
       call ipgeti ('b2tfhe_no_current', no_current)
       call ipgeti ('b2news_guard_flows', guard_flows)
       if (redef_gmtry.eq.0) then
         fch_inert_style = 1
       else
         fch_inert_style = 2
       end if
       call ipgeti ('b2tiner_fch_inert_style', fch_inert_style)          !srv 07.07.99
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

c Drift terms
      fchinert =0.0e0_R8
      fchinert_a =0.0e0_R8
      if(fch_inert.ne.0.0e0_R8.and.maxval(facdrift).ne.0.0_R8.and.
     .   no_current.eq.0)then                                            !srv 07.07.99
       do is = 0, ns-1
        if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
         call b2tiner (nx, ny, ns, vol, gs, bb, na, ua,
     &                    facdrift, fchinert_a(-1,-1,0,is))
        else
         fchinert_a(:,:,:,is) = 0.0_R8
        endif
        fchinert(:,:,:) = fchinert(:,:,:) + fchinert_a(:,:,:,is)
       enddo
      endif                                                              !srv 07.07.99

      if(guard_flows.eq.0) then                                          !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fchinert(ix,iy,1)=0.0_R8
              fchinert_a(ix,iy,1,:)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                fchinert(bc_list_x(ix,ib),bc_list_y(ix,ib),0)=0.0_R8
                fchinert_a(bc_list_x(ix,ib),bc_list_y(ix,ib),0,:)=0.0_R8
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              fchinert(ix,-1,0) = 0.0_R8
              fchinert_a(ix,-1,0,:) = 0.0_R8
            endif
          enddo
        endif
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchinert(ix,iy,0)=0.0e0_R8
            fchinert_a(ix,iy,0,:)=0.0e0_R8
          end do
        end do
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchinert(ix,iy,1)=0.0e0_R8
            fchinert_a(ix,iy,1,:)=0.0e0_R8
          end do
        end do
      end do

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner_a

      end subroutine b2tiner_a                                            !srv 20.09.06 } 08.09.21 }
