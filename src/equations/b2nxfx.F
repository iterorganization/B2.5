*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfx (nx, ny, qe, vol, hx, bb, ne, te, po, ehxp,
     &     wrk0, wrk1)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   qe, bb(-1:nx,-1:ny,0:3), ne(-1:nx,-1:ny), te(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny),vol(-1:nx,-1:ny),hx(-1:nx,-1:ny),
     &  wrk0(-1:nx,-1:ny),wrk1(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   ehxp(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFX computes an expression that occurs in the electron-atom
*     thermal force.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine evaluates an expression that is symbolically
*     written as  qe*ne*ehatx*pbsx*hx , where ehat is the modified
*     electric field and the other factors are familiar. Some
*     interpolation from cell centers to x-faces is involved.
*     This expression occurs in the electron-atom thermal force.
*     pbsx is written as (vol/hx)*abs(bb(,,0))/bb(,,3) to replace it
*     with its cell-centred value.
*
*
*     ehxp - (-1:nx,-1:ny) real array, output.
*     For (ix,iy) in (0:nx-1,-1:ny), ehxp(ix,iy) is defined by:
*       ehxp(ix,iy) = pbs(ix,iy,0)*
*         (ne("ix+1/2",iy)*te("ix+1/2",iy)-
*          ne("ix-1/2",iy)*te("ix-1/2",iy)-
*          qe*ne(ix,iy)*(po("ix+1/2",iy)-po("ix-1/2",iy)))
*     For ix.eq.-1 or ix.eq.nx, ehxp(ix,-1:ny) will hold 0.
*     ehxp is proportional to the modified electric field, and occurs
*     in the electron-atom friction force.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy
      integer, save :: style = 1, iout = 0                            !srv 26.10.08
      integer, save :: pot_eq = 1, no_current=0                       !srv 16.07.12
      save ncall
*   ..procedures
      external subini, subend, xertst, ipgeti
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxfx')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti('b2nxfx_style', style)                            !srv 26.10.08
       call ipgeti('b2nxfx_iout', iout)
       call ipgeti('b2news_poteq', pot_eq)                               !srv 16.07.12
       call xertst(0.le.pot_eq.and.pot_eq.le.2,'faulty argument pot_eq') !srv 16.07.12 27.02.13
       call ipgeti ('b2tfhe_no_current', no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      if(pot_eq.eq.0 .or. pot_eq.eq.2) then                           !srv 16.07.12 27.02.13
       ehxp = 0.0_R8
      elseif(pot_eq.eq.1) then
       if(style.eq.0) then                                             !srv 26.10.08
*   ..compute wrk0 and wrk1
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk0(ix,iy) = 0.0_R8
          wrk1(ix,iy) = 0.0_R8
         else
          wrk0(ix,iy) = ne(ix,iy)*te(ix,iy)-
     &         ne(leftix(ix,iy),leftiy(ix,iy))
     &        *te(leftix(ix,iy),leftiy(ix,iy))
          wrk1(ix,iy) = po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy))
         endif
        enddo
       enddo
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
       call b2xxgs (nx, ny, 0.0_R8, wrk1, 1)
*   ..compute ehxp
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
           if(rightix(ix,iy).eq.nx+1) then
             ehxp(ix,iy) = 0.0_R8
           else
             ehxp(ix,iy) = 
     &        (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3)*        !srv 11.01.13
     &         vol(ix,iy)*
     &        (wrk0(rightix(ix,iy),rightiy(ix,iy))-qe*ne(ix,iy)*
     &         wrk1(rightix(ix,iy),rightiy(ix,iy)))/
     &        (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy)))
           endif
         elseif(rightix(ix,iy).eq.nx+1) then
          ehxp(ix,iy) = 
     &     (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3)*           !srv 11.01.13
     &      vol(ix,iy)*
     &     (wrk0(ix,iy)-qe*ne(ix,iy)*wrk1(ix,iy))/
     &     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         else
          ehxp(ix,iy) = 
     &     (vol(ix,iy)/hx(ix,iy))*bb(ix,iy,0)/bb(ix,iy,3)*           !srv 11.01.13
     &      vol(ix,iy)*
     &     ((wrk0(ix,iy)-qe*ne(ix,iy)*wrk1(ix,iy))/
     &      (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+
     &      (wrk0(rightix(ix,iy),rightiy(ix,iy))-qe*ne(ix,iy)*
     &       wrk1(rightix(ix,iy),rightiy(ix,iy)))/
     &      (vol(ix,iy)+vol(rightix(ix,iy),rightiy(ix,iy))))
         endif
        enddo
       enddo
      else ! style                                     !srv 26.10.08 {
*   ..compute wrk0 and wrk1
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk0(ix,iy) = 0.0e0_R8
          wrk1(ix,iy) = 0.0e0_R8
         else
          wrk0(ix,iy) = 
     &        (ne(ix,iy)*te(ix,iy)
     &        -ne(leftix(ix,iy),leftiy(ix,iy))*
     &         te(leftix(ix,iy),leftiy(ix,iy)))
     &       /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          wrk1(ix,iy) = 
     &        (po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))
     &       /(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
         endif
        enddo
       enddo
*   ..compute ehxp
       ehxp = 0.0e0_R8
       do iy = 0, ny-1
        do ix = 0, nx-1
         if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
          ehxp(ix,iy) = 0.0e0_R8
         else
          ehxp(ix,iy) = ((wrk0(ix,iy)+wrk0(leftix(ix,iy),leftiy(ix,iy)))
     &                  *0.5e0_R8-
     &                  (wrk1(ix,iy)+wrk1(leftix(ix,iy),leftiy(ix,iy)))
     &                  *0.5e0_R8*qe*ne(ix,iy))
     &                  *bb(ix,iy,0)/bb(ix,iy,3)*vol(ix,iy)
         endif
        enddo
       enddo
      endif                                               !srv 26.10.08 }
      else                                                !srv 16.07.12

      endif                                               !srv 16.07.12

      if (iout.eq.1) then                                 !srv 11.09.09
        call my_out(70,nx,ny,ehxp,'b2nxfx_ehxp')
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfx

      end subroutine b2nxfx
