      subroutine b2tlhi (nx, ny, mp, cfl, gs, pbsx, qc, ni, nirm,
     &  ti, chvimx)
      use b2mod_types
      use b2mod_indirect
      use b2mod_external
      use b2mod_geo
     , , only : pitx
      use b2mod_subsys
      implicit none
      integer nx, ny
      real (kind=R8) ::
     *  mp, cfl(-1:nx,-1:ny),
     *  gs(-1:nx,-1:ny,0:2),
     *  pbsx(-1:nx,-1:ny), qc(-1:nx,-1:ny,0:1), ni(-1:nx,-1:ny,0:1),
     *  nirm(-1:nx,-1:ny), ti(-1:nx,-1:ny), chvimx(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TLHI computes the all atom heat transport flux limit velocity.
*     ------------------------------------------------------------------
*   ..local variables
      integer ncall, n2, ix, iy, area_fix, model
      real (kind=R8) ::
     *   t0, t1
      save ncall, area_fix, model
*   ..procedures
      intrinsic max, sqrt
      external xertst, b2xvfx, ipgeti
*   ..initialisation
      data ncall/0/, area_fix/3/, model/1/
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlhi')
*   ..set n2
      n2=(nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx, ny')
*   ..test mp, cflmi
      call xertst (0.0_R8.lt.mp, 'faulty argument mp')
      call b2xvsg (n2, cfl, 1, 'cflmi', '.ge.')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tqca_model', model)
      endif

*   ..compute flux limit velocity
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         chvimx(ix,iy) = 0.0_R8
        else
cwdk perhaps the if-test should be based on a 'mode' switch, instead of b2tqca_model,
c    because now both options are interchangeable
          if (model.eq.1) then  ! SOLPS5.0 treatments
            t0 = max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy))
            chvimx(ix,iy) = cfl(ix,iy)*abs(pbsx(ix,iy))*sqrt(t0/mp)*       !srv 11.01.13
     &        (nirm(leftix(ix,iy),leftiy(ix,iy))
     &        +nirm(ix,iy))/(ni(leftix(ix,iy),leftiy(ix,iy),0)+
     &                    ni_ext(leftix(ix,iy),leftiy(ix,iy),0)+
     &                    ni(ix,iy,0)+ni_ext(ix,iy,0))
            if (area_fix.eq.2) chvimx(ix,iy) = chvimx(ix,iy)*qc(ix,iy,0)
          else if (model.eq.2) then  ! SOLPS4.0 treatment
            t0 = ti(ix,iy)
            t1 = ti(leftix(ix,iy),leftiy(ix,iy))
cwdk 10.03.16 Divide by ni*ti at cell faces, so chvimx contains the same information
c           as case b2tqca_model.eq.1
c           Removes need for different treatment later on in the code (b2tfhi), and
c           ensures consistency with b2trcl_conductive_limit .ne. 0
            chvimx(ix,iy) = cfl(ix,iy)*gs(ix,iy,0)*
     &       (nirm(ix,iy)*pitx(ix,iy)*sqrt(t0/mp)*t0 +
     &        nirm(leftix(ix,iy),leftiy(ix,iy))*
     &        pitx(leftix(ix,iy),leftiy(ix,iy))*sqrt(t1/mp)*t1)/
     &       ((ni(leftix(ix,iy),leftiy(ix,iy),0)+
     &         ni_ext(leftix(ix,iy),leftiy(ix,iy),0)+
     &         ni(ix,iy,0)+ni_ext(ix,iy,0))*max(t0,t1))
            if (area_fix.ge.2) chvimx(ix,iy) = chvimx(ix,iy)*qc(ix,iy,0)
          endif
        endif
       enddo
      enddo
      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlhi
