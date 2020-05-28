*-----------------------------------------------------------------------
*.specification

      subroutine b2tfcc (nx, ny, ns, na, te, ti, cdpa, floe, cone,
     &  floi, coni)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  cdpa(-1:nx,-1:ny,0:1,0:ns-1)
*   ..input/output arguments
      real (kind=R8) ::
     *   floe(-1:nx,-1:ny,0:1), cone(-1:nx,-1:ny,0:1,0:2),           !srv 02.01.07
     &  floi(-1:nx,-1:ny,0:1), coni(-1:nx,-1:ny,0:1,0:2)             !srv 02.01.07
*   ..common blocks
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is
      real (kind=R8) ::
     *   xfac
      save ncall, xfac
*   ..procedures
      external subini, subend, xertst, ipgetr
      external b2xvsg, b2xvff
*   ..initialisation
      data ncall/0/, xfac/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfcc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tfcc_xfac', xfac)
       call xertst (0.le.xfac, 'faulty internal parameter xfac')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of cdpa
       do is = 0, ns-1
        call b2xvff (nx, ny, cdpa(-1,-1,0,is), 'cdpa')
       enddo
*    ..test edge values of floe, cone, floi, coni
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0), 'cone0')                !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,2), 'coni2')                !srv 02.01.07 }
*    ..test sign of cdpa, cone, coni
       call b2xvsg (n2*2*ns, cdpa, 1, 'cdpa', '.ge.')
       call b2xvsg (n2*2, cone, 1, 'cone', '.ge.')
       call b2xvsg (n2*2, coni, 1, 'coni', '.ge.')
      endif

* ..correct floe, cone, floi, coni
*   ..correct x-components
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2) then
          floe(ix,iy,0) = floe(ix,iy,0)-
     &     (xfac*2.5_R8/2.0_R8)*cdpa(ix,iy,0,is)*
     &     ((rza(ix,iy,is)+1.0_R8)*rza(ix,iy,is)*na(ix,iy,is)-
     &      (rza(leftix(ix,iy),leftiy(ix,iy),is)+1.0_R8)*
     &       rza(leftix(ix,iy),leftiy(ix,iy),is)*
     &       na(leftix(ix,iy),leftiy(ix,iy),is))*
     &      (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
          cone(ix,iy,0,0) = cone(ix,iy,0,0)+                         !srv 02.01.07
     &     (xfac*2.5_R8/4.0_R8)*cdpa(ix,iy,0,is)*
     &     ((rza(ix,iy,is)+1.0_R8)*rza(ix,iy,is)*na(ix,iy,is)+
     &      (rza(leftix(ix,iy),leftiy(ix,iy),is)+1.0_R8)*
     &       rza(leftix(ix,iy),leftiy(ix,iy),is)*
     &       na(leftix(ix,iy),leftiy(ix,iy),is))*
     &     (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
          floi(ix,iy,0) = floi(ix,iy,0)-
     &     (xfac*2.5_R8/2.0_R8)*cdpa(ix,iy,0,is)*
     &     ((rza(ix,iy,is)+1.0_R8)*na(ix,iy,is)-
     &      (rza(leftix(ix,iy),leftiy(ix,iy),is)+1.0_R8)*
     &      na(leftix(ix,iy),leftiy(ix,iy),is))*
     &     (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))
          coni(ix,iy,0,0) = coni(ix,iy,0,0)+                         !srv 02.01.07
     &     (xfac*2.5_R8/4.0_R8)*cdpa(ix,iy,0,is)*
     &     ((rza(leftix(ix,iy),leftiy(ix,iy),is)+1.0_R8)*
     &       na(leftix(ix,iy),leftiy(ix,iy),is)+
     &      (rza(ix,iy,is)+1.0_R8)*na(ix,iy,is))*
     &     (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))
         endif
        enddo
       enddo
      enddo
*   ..correct y-components
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).ne.-2) then
          floe(ix,iy,1) = floe(ix,iy,1)-
     &     (xfac*2.5_R8/2.0_R8)*cdpa(ix,iy,1,is)*
     &     ((rza(ix,iy,is)+1.0_R8)*rza(ix,iy,is)*na(ix,iy,is)-
     &      (rza(bottomix(ix,iy),bottomiy(ix,iy),is)+1.0_R8)*
     &      rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &      na(bottomix(ix,iy),bottomiy(ix,iy),is))*
     &     (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))
          cone(ix,iy,1,0) = cone(ix,iy,1,0)+                         !srv 02.01.07
     &     (xfac*2.5_R8/4.0_R8)*cdpa(ix,iy,1,is)*
     &     ((rza(bottomix(ix,iy),bottomiy(ix,iy),is)+1.0_R8)*
     &      rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &      na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &      (rza(ix,iy,is)+1.0_R8)*rza(ix,iy,is)*na(ix,iy,is))*
     &     (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))
          floi(ix,iy,1) = floi(ix,iy,1)-
     &     (xfac*2.5_R8/2.0_R8)*cdpa(ix,iy,1,is)*
     &     ((rza(ix,iy,is)+1.0_R8)*na(ix,iy,is)-
     &      (rza(bottomix(ix,iy),bottomiy(ix,iy),is)+1.0_R8)*
     &      na(bottomix(ix,iy),bottomiy(ix,iy),is))*
     &     (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))
          coni(ix,iy,1,0) = coni(ix,iy,1,0)+                         !srv 02.01.07
     &     (xfac*2.5_R8/4.0_R8)*cdpa(ix,iy,1,is)*
     &     ((rza(bottomix(ix,iy),bottomiy(ix,iy),is)+1.0_R8)*
     &      na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &      (rza(ix,iy,is)+1.0_R8)*na(ix,iy,is))*
     &     (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))
         endif
        enddo
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfcc

      end
