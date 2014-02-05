*-----------------------------------------------------------------------
*.specification

      subroutine b2upht_ (nx, ny, ns, rxf, pcm0, vol, cortt, corte, 
     &  corti, pccm0, na, ua, te, ti, po, ne, ni, wrk0)                 !srv 01.07.05
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
	use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   rxf, pcm0, cortt(-1:nx,-1:ny), corte(-1:nx,-1:ny),
     &  corti(-1:nx,-1:ny), pccm0(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), vol(-1:nx,-1:ny)
*   ..input/output arguments
      real (kind=R8) ::
     *   ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny)
*   ..output arguments
      real (kind=R8) ::
     &   wrk0(-1:nx,-1:ny)
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
      integer ncall, n2, stylec, ix, iy, is
	integer, save :: pot_eq = 1                             !srv 13.07.05
      real (kind=R8) ::
     *   ctt, dte, dti, t0, t1, t2
      save ncall, stylec
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *   med
      external subini, subend, xertst, ipgeti
      external b2xvfx
      med(t0,t1,t2) = min(max(t0,t1),max(t1,t2),max(t2,t0))
*   ..initialization
      data ncall/0/, stylec/0/
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
       call ipgeti ('b2news_poteq', pot_eq)                           !srv 13.07.05
       call xertst (0.le.stylec, 'faulty internal parameter stylec')
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test rxf, pcm0
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
      call xertst (0.le.pcm0, 'faulty argument pcm0')
*   ..test sign of na, ni, vol, te, ti, ne
      if (ncall.lt.3) then
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
        call b2xvsg (n2, vol, 1, 'vol', '.gt.')
        call b2xvsg (n2, te, 1, 'te', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.')
        call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

* ..apply the correction
*   ..adjust te, ti, po
      do iy = -1, ny
       do ix = -1, nx
        if (stylec.eq.0) then
         dte = rxf*max(-te(ix,iy),(corte(ix,iy)+cortt(ix,iy))/2)
         dti = rxf*max(-ti(ix,iy),(corti(ix,iy)+cortt(ix,iy))/2)
        else if (stylec.eq.1) then
         ctt = cortt(ix,iy)-
     &     med(corte(ix,iy),cortt(ix,iy),corti(ix,iy))
         dte = rxf*max(-te(ix,iy),corte(ix,iy)+ctt)
         dti = rxf*max(-ti(ix,iy),corti(ix,iy)+ctt)
        else if (stylec.eq.2) then
*        !! This for testing purposes only
         dte = rxf*max(-te(ix,iy),corte(ix,iy))
         dti = 0
        else if (stylec.eq.3) then
*        !! This for testing purposes only
         dte = 0
         dti = rxf*max(-ti(ix,iy),corti(ix,iy))
        else if (stylec.eq.4) then
*        !! This for testing purposes only
         dte = rxf*max(-te(ix,iy),corte(ix,iy))
         dti = rxf*max(-ti(ix,iy),corti(ix,iy))
        else if (stylec.eq.5) then                                    !srv 02.01.07
*       
         dte = corte(ix,iy)
         dti = corti(ix,iy)
        endif
        if (leftix(ix,iy).ne.-2) 
     >   wrk0(ix,iy) = ne(ix,iy)*dte+ni(ix,iy,0)*dti-
     &                 ne(leftix(ix,iy),leftiy(ix,iy))*dte-
     &                 ni(leftix(ix,iy),leftiy(ix,iy),0)*dti
        te(ix,iy) = te(ix,iy)+dte
        ti(ix,iy) = ti(ix,iy)+dti
	  if(stylec.ne.5) then                                          !srv 02.01.07 {
	   if(pot_eq .ne. 0) po(ix,iy) = po(ix,iy)+dte/qe               !srv 17.07.05
	  endif                                                         !srv 02.01.07 }
       enddo
      enddo
*   ..adjust ua
      do iy = -1, ny
       do ix = -1, nx
        if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 17.07.05 {
         t0 = pcm0*pccm0(ix,iy)*vol(ix,iy)*
     &    (wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &     wrk0(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &       vol(rightix(ix,iy),rightiy(ix,iy))))
        elseif(rightix(ix,iy).eq.nx+1) then
         t0 = pcm0*pccm0(ix,iy)*vol(ix,iy)*
     &    wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
        elseif(leftix(ix,iy).eq.-2) then
         t0 = 0.0_R8
        endif
        do is = 0, ns-1
         ua(ix,iy,is) = ua(ix,iy,is)-t0
        enddo
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht_

      end
