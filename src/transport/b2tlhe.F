      subroutine b2tlhe (nx, ny, me, cflme, bb, gs, pbsx, qc, ne, te, 
     &  chvemx)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  me, cflme, 
     &  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2),
     &  pbsx(-1:nx,-1:ny), qc(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), chvemx(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TLHE computes the electron heat transport flux limit velocity.
*     ------------------------------------------------------------------
*   ..local variables
      integer ncall, n2, ix, iy, area_fix, model
      real (kind=R8) ::
     *   t0, t1
      save ncall, area_fix, model
*   ..procedures
      real (kind=R8) :: pit
      intrinsic max, sqrt
      external subini, subend, xertst, b2xvsg, ipgeti
      pit(ix,iy) = abs(bb(ix,iy,0))/bb(ix,iy,3)
*   ..initialisation
      data ncall/0/, area_fix/3/, model/1/
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlhe')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx, ny')
*   ..test me, cflme
      call xertst (0.0_R8.lt.me, 'faulty argument me')
      call xertst (0.0_R8.le.cflme, 'faulty argument cflme')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tqce_model', model)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif

*   ..compute flux limit velocity
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
          chvemx(ix,iy) = 0.0_R8
        else
cwdk perhaps the if-test should be based on a 'mode' switch, instead of b2tqce_model, 
c    because now both options are interchangeable
          if (model.eq.1.or.model.eq.3) then  ! SOLPS5.0 treatments
            t0 = max(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
            chvemx(ix,iy) = cflme*sqrt(t0/me)*abs(pbsx(ix,iy))        !srv 11.01.13
            if (area_fix.eq.2) chvemx(ix,iy) = chvemx(ix,iy)*qc(ix,iy)
          else if (model.eq.2) then  ! SOLPS4.0 treatment
            t0 = te(ix,iy)
            t1 = te(leftix(ix,iy),leftiy(ix,iy))
cwdk 10.03.16 Divide by ne*te at cell faces, so chvemx contains the same information 
c           as case b2tqce_model.eq.1
c           Removes need for different treatment later on in the code (b2tfhe), and 
c           ensures consistency with b2trcl_conductive_limit .ne. 0
	    chvemx(ix,iy) = cflme * gs(ix,iy,0) *
     .	     (ne(ix,iy)*pit(ix,iy)*sqrt(t0/me)*t0 +
     .	      ne(leftix(ix,iy),leftiy(ix,iy))*
     .        pit(leftix(ix,iy),leftiy(ix,iy))*sqrt(t1/me)*t1)/
     .       ((ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*max(t0,t1))  
            if (area_fix.ge.2) chvemx(ix,iy) = chvemx(ix,iy)*qc(ix,iy)
          endif
        endif
       enddo
      enddo
      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlhe
