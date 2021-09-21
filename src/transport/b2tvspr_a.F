      subroutine b2tvspr_a(nx, ny, ns, fac_vis, vol,                      !srv 10.10.17 08.09.21 {
     &   po, ti, na, avisper, bvisper, cvisper, fchvisper_a, fchvisper)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo_corner
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &   fac_vis(-1:nx,-1:ny), vol(-1:nx,-1:ny),
     &   po(-1:nx,-1:ny), ti(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),            !srv 27.01.00
     &   cvisper(-1:nx,-1:ny,0:1), fchvisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   fchvisper(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSPR_A computes the X and Y components of the perpendicular viscosity
*     current for each species and calculates the X and Y components of 
*     the total perpendicular viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
      integer, save :: ncall = 0, iout = 0                              !srv 06.12.17
      integer, save :: no_current = 0
      real (kind=R8), save :: vis_per = 1.0e0_R8                        !srv 22.11.99
      real (kind=R8) :: fac_vism
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvspr_a')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                          !srv 22.11.99
       call ipgeti ('b2tfhe_no_current', no_current)
       call ipgeti ('b2tvspr_iout', iout)                               !srv 06.12.17
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..calculate the coefficients in difference approximation of viscosity perpendicular current
      call b2txvspr(nx, ny, ns, fac_vis, vol,                           !srv 10.10.17
     &              avisper, bvisper, cvisper)
* ..calculate the contribution of viscosity perpendicular current       !srv 22.11.99 
      fchvisper=0.0e0_R8                                                !srv 22.11.99
      fchvisper_a=0.0e0_R8
        
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_per.ne.0.0e0_R8 .and.
     .   no_current.eq.0) then
       do is = 0, ns-1
        if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
          call b2tvspr (nx, ny, ns, fac_vis, vol, po, ti, na,           !srv 10.10.17
     &     avisper, bvisper, cvisper, fchvisper_a(-1,-1,0,is))
        else
          fchvisper_a(:,:,:,is)=0.0_R8
        endif
        fchvisper(:,:,:) = fchvisper(:,:,:) + fchvisper_a(:,:,:,is)
       enddo
      endif                                                             !srv 22.11.99

      if(iout.eq.1) then
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tvspr_a_fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tvspr_a_fchvispery')
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvspr_a

      end subroutine b2tvspr_a                                           !srv 10.10.17 08.09.21 }

