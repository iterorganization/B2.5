*-----------------------------------------------------------------------
*.specification

      subroutine b2tstch (nx, ny, vol, po, te, ne,                      !srv 07.07.14
     *           csig_stoch, fchstoch)                                  !lk 31.01.08 {
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_geo
     , , only : hy1
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), po(-1:nx,-1:ny),
     1 te(-1:nx,-1:ny), ne(-1:nx,-1:ny), csig_stoch(-1:nx,-1:ny,0:1)    !srv 24.08.21
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     & fchstoch(-1:nx,-1:ny,0:1)                                        !srv 24.08.21
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TSTCH computes the Y component of the electron current
*     associated with stochastization of the magnetic field.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, ix, iy, bottomdelta1, ixref
      integer jxi, jxa, jsep
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8                 !lk 31.01.08 }
      real (kind=R8) :: sig_stoch, delta_stoch, width
      save ncall, jsep, ixref, sig_stoch, delta_stoch, bottomdelta1
*   ..procedures
      external ipgetr, ipgeti, xertst, xerrab
      external b2xvsg, b2xvff, b2xvfx
      external get_jsep
*   ..initialization
      data ncall/0/, delta_stoch/0.0e0_R8/, sig_stoch/1.0e0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tstch')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2tstch_delta', delta_stoch)
        call ipgetr ('b2tstch_sig', sig_stoch)
        call xertst (0.0_R8.le.sig_stoch,
     &   'faulty internal parameter sig_stoch')
        call ipgetr ('b2tfhe_stochastic', fch_stochastic)                !lk 31.01.08
        call xertst (0.0_R8.lt.delta_stoch.or.fch_stochastic.eq.0.0_R8,
     &   'faulty internal parameter delta_stoch')
        call get_jsep(nx,ny,jxi,jxa,jsep)
        ixref=jxa
        call ipgeti ('b2tqna_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument b2tqna_ixref')

        bottomdelta1 = jsep+1
        width = 0.0e0_R8
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

!! FIXME: Consider adding XY current and fourth dimension
      fchstoch   = 0.0e0_R8
      if (fch_stochastic .ne. 0.0e0_R8 ) then                           !lk 31.01.08 {
        call b2trsc (nx, ny, vol, po, te, ne, csig_stoch)               !srv 24.08.21
        do iy = bottomdelta1, jsep+1
          do ix = -1, nx
           if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then
            fchstoch(ix,iy,1) = -csig_stoch(ix,iy,1)*                   !srv 24.08.21
     &       ((po(ix,iy) - po(bottomix(ix,iy),bottomiy(ix,iy)))-
     &       0.5e0_R8/qe*
     &       (te(ix,iy) - te(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -(te(ix,iy) + te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8/
     &        qe
     &      *(log(ne(ix,iy))-
     &        log(ne(bottomix(ix,iy),bottomiy(ix,iy)))))
           endif
          enddo
        enddo
      endif                                                             !lk 31.01.08 }

!      call my_out(70,nx,ny,fchstoch(-1,-1,1),'b2tstch_fchstochy')      !srv 24.08.21

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tstch

      end subroutine b2tstch                                            !lk 31.01.08 }
