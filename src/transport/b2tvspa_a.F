*-----------------------------------------------------------------------
*.specification

      subroutine b2tvspa_a (nx, ny, ns, bb, ua,                          !srv 20.09.06 { 16.10.17 08.09.21 {
     &           cvsa_cl, fac_vis, fchvispar_a, fchvispar)               !srv 22.06.08
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
     & ua(-1:nx,-1:ny,0:ns-1),
     & cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 16.10.17
     1 bb(-1:nx,-1:ny,0:3), fac_vis(-1:nx,-1:ny)                        !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchvispar_a(-1:nx,-1:ny,0:1,0:ns-1), fchvispar(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSPA_A computes the X and Y components of the parallel viscosity
*     current for each species and calculates the X and Y components of 
*     the total parallel viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer ncall, ix, iy, ib, is
      integer, save :: no_current = 0
      real (kind=R8) ::
     &  vis_par
      save ncall, vis_par
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, vis_par/1.0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspa_a')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_vis_par', vis_par)
       call ipgeti ('b2tfhe_no_current', no_current)
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of cvsa_cl
        call b2xvff (nx, ny, cvsa_cl, 'cvsa_cl')                        !srv 16.10.17
      endif

c Drift terms
      fchvispar =0.0e0_R8                                               !srv 09.07.99
      fchvispar_a=0.0e0_R8
      if (vis_par.ne.0.0e0_R8.and.maxval(fac_vis).ne.0.0_R8.and.
     .   no_current.eq.0) then
       do is = 0, ns-1
        if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
         call b2tvspa (nx, ny, ns, bb, ua, cvsa_cl, fac_vis,
     &                 fchvispar_a(-1,-1,0,is))
        else
         fchvispar_a(:,:,:,is)=0.0_R8
        endif
        fchvispar(:,:,:) = fchvispar(:,:,:) + fchvispar_a(:,:,:,is)
       enddo
      endif                                                             !srv

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchvispar(ix,iy,0)=0.0e0_R8
            fchvispar_a(ix,iy,0,:)=0.0_R8
          end do
        end do
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fchvispar(ix,iy,1)=0.0e0_R8
            fchvispar_a(ix,iy,1,:)=0.0_R8
          end do
        end do
      end do

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspa_a

      end subroutine b2tvspa_a                                            !srv 20.09.06 } 08.09.21 }
