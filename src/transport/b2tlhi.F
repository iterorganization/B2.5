      subroutine b2tlhi (nx, ny, mp, cflmi, pbsx, qc, ni, nirm, ti,
     &  chvimx)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  mp, cflmi, 
     *  pbsx(-1:nx,-1:ny), qc(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  nirm(-1:nx,-1:ny), ti(-1:nx,-1:ny), chvimx(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TLHI computes the all atom heat transport flux limit velocity.
*     ------------------------------------------------------------------
*   ..local variables
      integer ncall, ix, iy, area_fix
      real (kind=R8) ::
     *   t0
      save ncall, area_fix
*   ..procedures
      intrinsic max, sqrt
      external subini, subend, xertst, ipgeti
      external b2xvfx
*   ..initialisation
      data ncall/0/, area_fix/3/
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlhi')
*   ..test mp, cflmi
      call xertst (0.0_R8.lt.mp, 'faulty argument mp')
      call xertst (0.0_R8.le.cflmi, 'faulty argument cflmi')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
      endif

*   ..compute flux limit velocity
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         chvimx(ix,iy) = 0.0_R8
        else
         t0 = max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy))
         chvimx(ix,iy) = cflmi*abs(pbsx(ix,iy))*sqrt(t0/mp)*            !srv 11.01.13
     &    (nirm(leftix(ix,iy),leftiy(ix,iy))
     &     +nirm(ix,iy))/(ni(leftix(ix,iy),leftiy(ix,iy),0)
     &                                 +ni(ix,iy,0))
         if (area_fix.eq.2) chvimx(ix,iy) = chvimx(ix,iy)*qc(ix,iy)
        endif
       enddo
      enddo
      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlhi
