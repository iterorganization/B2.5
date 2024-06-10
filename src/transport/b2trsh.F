*-----------------------------------------------------------------------
*.specification

      subroutine b2trsh (nx, ny, vol, hy1, po, te, ne, csig_stoch,       !srv 24.08.21
     &                   chce_stoch)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     & hy1(-1:nx,-1:ny), vol(-1:nx,-1:ny), po(-1:nx,-1:ny),
     1 te(-1:nx,-1:ny), ne(-1:nx,-1:ny), csig_stoch(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  chce_stoch(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRSH computes the X and Y components of electron thermal conductivity
*     associated with stochastization of the magnetic field.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer n2, ncall, ix, iy, bottomdelta1, ixref
      integer jxi, jxa, jsep
      real (kind=R8), save :: fch_stochastic = 0.0_R8,                   !lk 31.01.08 }
     &                        alpha_stoch = 1.0_R8                       !srv 17.12.13
      real (kind=R8) :: sig_stoch, delta_stoch, width
      save ncall, jsep, ixref, sig_stoch, delta_stoch, bottomdelta1
*   ..procedures
      external ipgetr, ipgeti, xertst, xerrab
      external b2xvsg, b2xvff, b2xvfx
      external get_jsep
*   ..initialization
      data ncall/0/, delta_stoch/0.0_R8/, sig_stoch/1.0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trsh')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign
       call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tstch_delta', delta_stoch)
        call ipgetr ('b2tstch_sig', sig_stoch)
        call xertst (0.0_R8.le.sig_stoch,
     &   'faulty internal parameter sig_stoch')
        call ipgetr ('b2tfhe_stochastic', fch_stochastic)                !lk 31.01.08
        call xertst (0.0_R8.lt.delta_stoch.or.fch_stochastic.eq.0.0_R8,
     &   'faulty internal parameter delta_stoch')
        call ipgetr ('b2trno_con_e_stochastic', alpha_stoch)             !lk 31.01.08
        call get_jsep(nx,ny,jxi,jxa,jsep)
        ixref=jxa
        call ipgeti ('b2tqna_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument b2tqna_ixref')

        bottomdelta1 = jsep+1
        width = 0.0_R8
        do iy = jsep,0,-1
          width = width + hy1(ixref,iy)
          if(width.gt.delta_stoch) then
            bottomdelta1 = iy+1
            exit
          endif
        enddo
        bottomdelta1 = bottomdelta1 + 1         !srv 17.12.13   plus 1 specially
!       write(*,"(a,3i5)") ' b2tstch: ixref, bottomdelta1, topcut(1)',
!     &                               ixref, bottomdelta1, topcut(1)
      endif

      chce_stoch = 0.0_R8
      do iy = bottomdelta1, jsep+1
        do ix = -1, nx
         if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then
          chce_stoch(ix,iy,1)=alpha_stoch*1.6_R8*
     &      csig_stoch(ix,iy,1)/qe/qe*
     &     (te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
         endif
        enddo
      enddo

!      call my_out(70,nx,ny,chce_stoch(-1,-1,1),'b2tstch_chcey_stoch')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2trsh

      end subroutine b2trsh                                            !srv 24.08.21 }
