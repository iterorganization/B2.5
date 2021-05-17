      subroutine b2tlmv (nx, ny, isb, cflmv, pbsx, qc, nb, tb, ub,
     &   cvsbx, cvsbxhz, fl)                                         !srv 09.01.01
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      implicit none
      integer nx, ny, isb
      real (kind=R8) ::
     *  cflmv, pbsx(-1:nx,-1:ny), qc(-1:nx,-1:ny),
     &  nb(-1:nx,-1:ny), tb(-1:nx,-1:ny),
     &  ub(-1:nx,-1:ny), cvsbx(-1:nx,-1:ny), cvsbxhz(-1:nx,-1:ny),   !srv 09.01.01
     &  fl(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2TLMV applies a flux limit to the viscosity.
*     ------------------------------------------------------------------
      integer ncall, n2, ix, iy, area_fix
      integer, save :: style = 1, iout = 0                           !srv 30.05.08
      character charns*3                                             !srv 12.04.04
      save ncall, area_fix
      real (kind=R8) ::
     *   qcl, qfl, tfl
      intrinsic abs
      external subini, subend, xertst, ipgeti
      external  b2xvsg
*   ..initialisation
      data ncall/0/, area_fix/3/
*     ------------------------------------------------------------------
      call subini ('b2tlmv')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test cflmv
      call xertst (0.0_R8.le.cflmv, 'faulty argument cflmv')
*   ..read switches
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tlmv_style', style)                          !srv 30.05.08
        call ipgeti ('b2tlmv_iout', iout)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, tb
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
       call b2xvsg (n2, tb, 1, 'tb', '.gt.')
      endif
*   ..main computation
      if(style.eq.0) then                                            !srv 30.05.08
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2) then
          qcl = -cvsbx(ix,iy)*
     &          (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
          qfl = cflmv*abs(pbsx(ix,iy))*                              !srv 11.01.13
     &     (nb(leftix(ix,iy),leftiy(ix,iy))*
     &      tb(leftix(ix,iy),leftiy(ix,iy))+
     &      nb(ix,iy)*tb(ix,iy))/2.0_R8
          if (area_fix.eq.2) qfl = qfl*qc(ix,iy)
          if (qfl.gt.0.0_R8) then
            tfl = 1.0_R8 + abs(qcl/qfl)
          else
            tfl = 1.0_R8
          endif
          cvsbx(ix,iy) = cvsbx(ix,iy)/tfl
          cvsbxhz(ix,iy) = cvsbxhz(ix,iy)/tfl                        !srv 02.07.08
          fl(ix,iy) = tfl
         else
          fl(ix,iy) = 1.0_R8
         endif
        enddo
       enddo
csw 18sep2014 TESTING only (SOLPS4 discr), i.e. switch off mom.flux
c             limit
c      else if(style.eq.2) then
c       fl(ix,iy)=1.0d0
      else                                                           !srv 30.05.08
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2) then
          if(.not.(.not.is_neutral(isb) .and.
     &              on_closed_surface(ix,iy))) then                  !srv 12.04.08
           qcl = -cvsbx(ix,iy)*
     &          (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
           qfl = cflmv*abs(pbsx(ix,iy))*                             !srv 11.01.13
     &          (nb(leftix(ix,iy),leftiy(ix,iy))*
     &           tb(leftix(ix,iy),leftiy(ix,iy))+
     &           nb(ix,iy)*tb(ix,iy))/2.0_R8
           if (area_fix.eq.2) qfl = qfl*qc(ix,iy)
           if (qfl.gt.0.0_R8) then
             tfl = 1.0_R8 + abs(qcl/qfl)
           else
             tfl = 1.0_R8
           endif
           cvsbx(ix,iy) = cvsbx(ix,iy)/tfl
           cvsbxhz(ix,iy) = cvsbxhz(ix,iy)/tfl
           fl(ix,iy) = tfl
          else
           fl(ix,iy) = 1.0_R8
          endif                                                      !srv 12.04.08
         else
          fl(ix,iy) = 1.0_R8
         endif
        enddo
       enddo
      endif                                                          !srv 30.05.08

      if(iout.ne.0) then                                             !srv 12.04.04 {
       write (charns,'(i3.3)') isb
       call my_out(70,nx,ny,fl,'cflmv-1'//charns)
      endif                                                          !srv 12.04.04 }

      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlmv
