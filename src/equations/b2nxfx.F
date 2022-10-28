*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfx (nCv, nFc, nVx, switch, geo, mpg, qe, ne, te,
     &     po, ehxp, wrk0, wrk1)
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   qe, ne(nCv), te(nCv), po(nCv), wrk0(nFc), wrk1(nFc)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   ehxp(nCv)

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
*     ehxp - (nFc) real array, output.
*     For (iFc) in (1:nFc), ehxp(iFc) is defined by:
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
      real (kind=R8) :: wrkc0(nCv), wrkc1(nCv), wrkv(nVx)
*   ..procedures
      external xertst
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2nxfx')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      if(switch%pot_eq.eq.0 .or. switch%pot_eq.eq.2) then               !srv 16.07.12 27.02.13
       ehxp = 0.0_R8
      elseif(switch%pot_eq.eq.1) then
      if(switch%b2nxfx_style.eq.0) then                                 !srv 26.10.08
       call xerrab ('b2nxfx -- style 0 not supported for WG')
#ifdef WG_TODO
*   ..compute wrk0 and wrk1
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrk0)
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, po   , wrkv, wrk1)
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
#endif
      else ! style                                                   !srv 26.10.08 {
*   ..compute wrk0 and wrk1
       call gradc_p(nCv, nFc, nVx, 0, geo, mpg, ne*te, wrkv, wrkc0)
       call gradc_p(nCv, nFc, nVx, 0, geo, mpg, po   , wrkv, wrkc1)
*   ..compute ehxp
       ehxp = (wrkc0-wrkc1*qe*ne)*geo%cvBb(:,0)/geo%cvBb(:,3)*geo%cvVol
*   ..zero out guard cells
       ehxp(mpg%nCi+1:nCv) = 0.0_R8
      endif                                                             !srv 26.10.08 }
      else                                                              !srv 16.07.12
        call xerrab ('wrong value of pot_eq')
      endif                                                             !srv 16.07.12

      if (switch%b2nxfx_iout.eq.1) then                                 !srv 11.09.09
        call my_out_us(70,nFc,1,ehxp,'b2nxfx_ehxp')
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfx

      end subroutine b2nxfx
