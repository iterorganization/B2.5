*-----------------------------------------------------------------------
*.specification

      subroutine b2upht (nx, ny, ns, nregionv, po_solve, ua_solve,
     &  solvereg, solvireg, solvtreg, solvpreg, solvmreg,                !srv 22.05.18
     &  rxf, pcm0, cortt, corte, corti,
     &  pccm0, na, ua, te, ti, po, ne, ni, wrk0, vol)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_numerics_namelist
     , , only : corr_core_dt, add_te_corr_to_po                ! IYS 19.12.2017, 21.09.2018
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv
      real (kind=R8) ::
     *  rxf, pcm0, cortt(-1:nx,-1:ny), corte(-1:nx,-1:ny),
     &  corti(-1:nx,-1:ny), pccm0(-1:nx,-1:ny), vol(-1:nx,-1:ny)
      logical po_solve, ua_solve,
     &  solvereg(0:nregionv), solvireg(0:nregionv),
     &  solvtreg(0:nregionv),                                            !srv 22.05.18
     &  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), wrk0(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPHT adjusts the state of the plasma towards the solution
*     of the heat equations.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, stylec, pot_eq, ix, iy, is, ireg 
     &  , nrings, mintopcut                ! IYS 27.03.2019
      real (kind=R8) ::
     *  ctt, dte, dti, t0, t1, t2
      real (kind=R8) :: 
     *  dte1(-1:nx,-1:ny), dti1(-1:nx,-1:ny),                             !lk 26.10.17
     *  tt0(-1:ny), tt1(-1:ny), tt2(-1:ny), tt3(-1:ny)                    !lk 26.10.17
*   ..procedures

      save ncall, stylec, pot_eq, nrings, mintopcut                ! IYS 27.03.2019
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *  med
      external subini, subend, xertst, ipgeti, b2xvsg
      med(t0,t1,t2) = min(max(t0,t1),max(t1,t2),max(t2,t0))
*   ..initialisation
      data ncall/0/, stylec/0/, pot_eq/1/, nrings/0/
*!!   (stylec is a numerical switch, needs experiments)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upht')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2upht_stylec', stylec)
       call ipgeti ('nrings_for_no_speedUP_averaging', nrings)     ! IYS 27.03.2019
       call xertst (0.le.stylec, 'faulty internal parameter stylec')
       call ipgeti ('b2news_poteq', pot_eq)
       mintopcut = minval(topcut(1:nncut))     ! IYS 27.03.2019
      endif
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8, 
     & 'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')
*   ..test sign of na, ni, vol, te, ti, ne
      if(ncall.lt.3) then
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (n2, vol, 1, 'vol', '.gt.')
        call b2xvsg (n2, te, 1, 'te', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.')
        call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif 

* ..apply the correction
*   ..adjust te, ti, po
      tt0 = 0.0e0_R8
      tt1 = 0.0e0_R8
      tt2 = 0.0e0_R8
      tt3 = 0.0e0_R8
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (stylec.eq.0) then
         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then  !srv 22.05.18
          dte = rxf*max(-te(ix,iy),(corte(ix,iy)+cortt(ix,iy))/2.0_R8)
          dti = rxf*max(-ti(ix,iy),(corti(ix,iy)+cortt(ix,iy))/2.0_R8)
         else
          dte = rxf*max(-te(ix,iy),corte(ix,iy))
          dti = rxf*max(-ti(ix,iy),corti(ix,iy))
         endif
        else if (stylec.eq.1) then
         if (solvereg(ireg).and.solvireg(ireg).and.solvtreg(ireg)) then  !srv 22.05.18
          ctt = cortt(ix,iy)-med(corte(ix,iy),cortt(ix,iy),corti(ix,iy))
          dte = rxf*max(-te(ix,iy),corte(ix,iy)+ctt)
          dti = rxf*max(-ti(ix,iy),corti(ix,iy)+ctt)
         else
          dte = rxf*max(-te(ix,iy),corte(ix,iy))
          dti = rxf*max(-ti(ix,iy),corti(ix,iy))
         endif
        else if (stylec.eq.2) then
*        !! This for testing purposes only
         dte = rxf*max(-te(ix,iy),corte(ix,iy))
         dti = 0.0_R8
        else if (stylec.eq.3) then
*        !! This for testing purposes only
         dte = 0.0_R8
         dti = rxf*max(-ti(ix,iy),corti(ix,iy))
        else if (stylec.eq.4) then
*        !! This for testing purposes only
         dte = rxf*max(-te(ix,iy),corte(ix,iy))
         dti = rxf*max(-ti(ix,iy),corti(ix,iy))
        endif
        dte1(ix,iy) = dte
        dti1(ix,iy) = dti
*   ..smooth the correction for inner flux surfaces      !lk 26.10.2017{
        if(on_closed_surface(ix,iy) .and. 
     *    iy .lt. mintopcut - nrings .and. bottomiy(ix,iy).ne.-2) then  ! IYS 27.03.2019
          tt0(iy) = tt0(iy) + vol(ix,iy)*ni(ix,iy,1)
          tt1(iy) = tt1(iy) + vol(ix,iy)*ni(ix,iy,1)*dti1(ix,iy)
          tt2(iy) = tt2(iy) + vol(ix,iy)*ne(ix,iy)
          tt3(iy) = tt3(iy) + vol(ix,iy)*ne(ix,iy)*dte1(ix,iy)
        endif
       enddo
      enddo

      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(on_closed_surface(ix,iy) .and. 
     *    iy .lt. mintopcut - nrings .and. bottomiy(ix,iy).ne.-2) then  ! IYS 27.03.2019
          dti1(ix,iy) = tt1(iy)/tt0(iy) +
     *         corr_core_dt*(dti1(ix,iy)-tt1(iy)/tt0(iy))    ! IYS 19.12.2017
          dte1(ix,iy) = tt3(iy)/tt2(iy) +
     *         corr_core_dt*(dte1(ix,iy)-tt3(iy)/tt2(iy))    ! IYS 19.12.2017
        endif
        dte = dte1(ix,iy)
        dti = dti1(ix,iy)                                       !srv 25.07.17
        if (leftix(ix,iy).ne.-2) 
     >   wrk0(ix,iy) = ne(ix,iy)*dte+ni(ix,iy,0)*dti-
     &                 ne(leftix(ix,iy),leftiy(ix,iy))*dte-
     &                 ni(leftix(ix,iy),leftiy(ix,iy),0)*dti
        if (solvereg(ireg)) te(ix,iy) = max(1.0e-5_R8*ev,te(ix,iy)+dte)   !xpb
        if ((solvpreg(ireg).and.po_solve).and.pot_eq.eq.1.and.
     &      leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &      bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1
     &                            .and. add_te_corr_to_po(ireg)) then
                            po(ix,iy) = po(ix,iy)+dte/qe
        else if (pot_eq.eq.2) then
          po(ix,iy) = 3.1_R8*te(ix,iy)/qe
        end if
        if (solvireg(ireg)) ti(ix,iy) = max(1.0e-5_R8*ev,ti(ix,iy)+dti)   !xpb
       enddo
      enddo
*   ..adjust ua
      if (ua_solve) then
       do iy = -1, ny
        do ix = -1, nx
         ireg = region(ix,iy,0)
         if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          t0 = pcm0*pccm0(ix,iy)*vol(ix,iy)*
     &     (wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &      wrk0(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
         elseif(rightix(ix,iy).eq.nx+1) then
          t0 = pcm0*pccm0(ix,iy)*vol(ix,iy)*
     &     wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         elseif(leftix(ix,iy).eq.-2) then
          t0 = pcm0*pccm0(ix,iy)*vol(ix,iy)*
     &     wrk0(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy)))
         endif
         do is = 0, ns-1
          if (solvmreg(is,ireg)) ua(ix,iy,is) = ua(ix,iy,is)-t0
         enddo
        enddo
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht

      end subroutine b2upht
