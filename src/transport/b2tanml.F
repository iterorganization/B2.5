*-----------------------------------------------------------------------
*.specification

      subroutine b2tanml (nx, ny, csig_an, po, fchanml)                 !srv 20.09.06 {
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
      integer nx, ny
      real (kind=R8) ::
     & po(-1:nx,-1:ny),
     1 csig_an(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     & fchanml(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TANML computes the x,y-components of anomalous current
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, ix, iy, ib
      integer, save :: no_current = 0                          !srv 27.02.13
      integer, save :: pot_eq
      integer, save :: guard_flows
      real (kind=R8) :: fch_anomalous
      save ncall, fch_anomalous
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_anomalous/1.0e0_R8/, pot_eq/1/, guard_flows/2/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tanml')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)
       if (fch_anomalous.eq.1.0_R8)
     >  call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgeti ('b2tfhe_no_current', no_current)                !srv 27.02.13
       call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2news_guard_flows', guard_flows)
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

      if(no_current.ne.0) then                                !srv 27.02.13
        fchanml = 0.0e0_R8
      else                                                    !srv 27.02.13
        if(fch_anomalous.ne.0.0e0_R8) then
          do iy = -1, ny
            do ix = -1, nx
              if(leftix(ix,iy).eq.-2)then
                fchanml(ix,iy,0)   = 0.0e0_R8
              else
                fchanml(ix,iy,0) = - csig_an(ix,iy,0)*
     &           (po(ix,iy) - po(leftix(ix,iy),leftiy(ix,iy)))
              endif
            enddo
          enddo

          do ix = -1, nx
            fchanml(ix,-1,1) = 0.0e0_R8
          enddo
          do iy = 0, ny
            do ix = -1, nx
              if(bottomiy(ix,iy).eq.-2)then
                fchanml(ix,iy,1) = 0.0e0_R8
              else
                fchanml(ix,iy,1) = - csig_an(ix,iy,1)*
     &           (po(ix,iy) - po(bottomix(ix,iy),bottomiy(ix,iy)))
              endif
            enddo
          enddo
        else
          fchanml = 0.0e0_R8
        endif
      endif                                                   !srv 27.02.13

      if(guard_flows.eq.0) then                                         !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fchanml(ix,iy,1)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                fchanml(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              fchanml(ix,-1,0) = 0.0_R8
            endif
          enddo
        endif
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchanml(ix,iy,0)=0.0e0_R8
          end do
        end do
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchanml(ix,iy,1)=0.0e0_R8
          end do
        end do
      end do

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tanml

      end subroutine b2tanml                                                            !srv 20.09.06 }
