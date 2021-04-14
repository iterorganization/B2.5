*-----------------------------------------------------------------------
*.specification

      subroutine b2xehx (nx, ny, te, po, ne, ehx)                        !srv 13.01.17
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     &   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehx(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XEHX computes the modified electrical field on x - cell faces.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: mode_ehx = 0                                      !srv 06.04.15
      integer, save :: pot_eq = 1, no_current= 0                         !srv 29.01.20
      real (kind=R8) ::
     *   t0
      save ncall
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialization
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xehx')
*   ..set internal parameters on first call                              !srv
      if (ncall.eq.0) then                                               !srv
       call ipgeti ('b2tfhe_mode_ehx', mode_ehx)                         !srv 06.04.15
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 29.01.20 {
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')       !srv 29.01.20 }
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif

      if(pot_eq.ne.1 .and. no_current.ne.0) then                         !srv 29.01.20 {
        ehx = 0.0_R8
      else

*   ..compute ehx
        do iy = -1, ny
          do ix = -1, nx                                                 !srv 26.11.02 {
            if(leftix(ix,iy).eq.-2) then                                 !srv 26.11.02
              ehx(ix,iy)  = 0.0e0_R8
            else                                                         !srv 26.11.02 }
              if(mode_ehx.eq.0) then                                     !srv 06.04.15 {
                ehx(ix,iy) =
     &         -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &          (ne(ix,iy)*te(ix,iy)                                     !srv 98
     &          -ne(leftix(ix,iy),leftiy(ix,iy))*                        !srv
     &           te(leftix(ix,iy),leftiy(ix,iy)))*                       !srv
     &          (2.0_R8/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))  !srv
              else if(mode_ehx.eq.1) then
                t0= (ne(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy)+
     &               ne(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
                ehx(ix,iy) =
     &             -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &              (ne(ix,iy)*te(ix,iy)
     &              -ne(leftix(ix,iy),leftiy(ix,iy))*
     &               te(leftix(ix,iy),leftiy(ix,iy)))/(qe*t0)
              else
                ehx(ix,iy) =
     &           -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &            (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))+            !srv
     &            (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))*0.5e0_R8*  !srv
     &            (log(ne(ix,iy))-
     &             log(ne(leftix(ix,iy),leftiy(ix,iy)))))/qe             !srv .98
              endif                                                      !srv 06.04.15 }
            endif                                                        !srv 26.11.02 02.04.15
          enddo
        enddo
      endif                                                              !srv 29.01.20 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehx

      end subroutine b2xehx
