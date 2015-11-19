*-----------------------------------------------------------------------
*.specification

      subroutine b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po,       !srv 20.09.06 {
     &                    fchanml)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), po(-1:nx,-1:ny), 
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 csig_an(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &    fchanml(-1:nx,-1:ny,0:1)
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
      integer ncall, n2, ix, iy, ib, is
      integer, save :: no_current=0                          !srv 27.02.13
      integer, save :: pot_eq
      real (kind=R8) :: fch_anomalous
      save ncall, fch_anomalous
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_anomalous/1.0e0_R8/, pot_eq/1/
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
       call ipgeti ('b2tfhe_no_current',no_current)                 !srv 27.02.13
       call ipgeti ('b2news_poteq', pot_eq)                         !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
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
              fchanml(ix,iy,1) = - csig_an(ix,iy,1)*
     &           (po(ix,iy) - po(bottomix(ix,iy),bottomiy(ix,iy)))
	    enddo
	  enddo
        else
          fchanml = 0.0e0_R8
        endif
      endif                                                   !srv 27.02.13

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
	DO IX=NXISO1(IB), NXISO2(IB)+1
	  DO IY=NYISO1(IB), NYISO2(IB)
	    fchanml(ix,iy,0)=0.0e0_R8
	  end do
	end do
	DO IX=NXISO1(IB), NXISO2(IB)
	  DO IY=NYISO1(IB), NYISO2(IB)+1
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

      end                                                             !srv 20.09.06 }
