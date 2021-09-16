      subroutine b2tvsq_a (nx, ny, ns, fac_vis,                         !srv 08.09.21 {
     &           ne2, lnlam, te, ti, ni, chci, fchvisq_a, fchvisq)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_geo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &   fac_vis(-1:nx,-1:ny),
     &   ne2(-1:nx,-1:ny), lnlam(-1:nx,-1:ny), te(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), ni(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchvisq_a(-1:nx,-1:ny,0:1,0:ns-1), fchvisq(-1:nx,-1:ny,0:1)    !srv 11.07.99
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TVSQ_A computes the X and Y components of the heat viscosity
*     current for each species and calculates the X and Y components of 
*     the total heat viscosity current.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, is
      integer, save :: ncall = 0
      integer, save :: no_current = 0
      real (kind=R8), save :: vis_q = 1.0e0_R8                          !srv 11.07.99
      real (kind=R8) :: fac_vism                                        !srv 16.10.17
*   ..procedures
      intrinsic abs
      external xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tvsq_a')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfhe_vis_q', vis_q)                              !srv 11.07.99
       call ipgeti ('b2tfhe_no_current', no_current)
      end if
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

      fchvisq   = 0.0e0_R8                                              !srv 11.07.99
      fchvisq_a = 0.0e0_R8
      fac_vism = maxval(fac_vis)
      if(fac_vism.ne.0.0_R8 .and. vis_q.ne.0.0e0_R8 .and.
     .   no_current.eq.0) then                                          !srv 16.10.17
       do is = 0, ns-1
        if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
         call b2tvsq (nx, ny, ns, fac_vis,
     &           ne2, lnlam, te, ti, ni, chci, fchvisq_a(-1,-1,0,is))
        else
         fchvisq_a(:,:,:,is)=0.0_R8
        endif
         fchvisq(:,:,:) = fchvisq(:,:,:) + fchvisq_a(:,:,:,is)
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tvsq_a

      end subroutine b2tvsq_a                                           !srv 08.09.21 }
