*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfc (nCv, nFc, nVx, isb, geo, mpg, rob, ub, flub,
     &  cvsb, flcb, cvcb, wrkf)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rob(nCv), ub(nCv),
     &  flub(nFc,0:1), cvsb(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flcb(nFc,0:1), cvcb(nFc,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrkf(nFc), flub0(nFc,0:1)                          !srv 02.07.08

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFC computes the transport terms in the parallel momentum
*     correction equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, iFc, ncallout, ib
      integer, save :: iout = 0                                     !srv 17.06.02
      character charns*3                                            !srv 17.06.02
      save ncall, style, ncallout
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr, b2scopy
      external b2xvff, b2xvsg, hbout, intfaceh
*   ..initialisation
      data ncall/0/, style/1/, ncallout/-1/
*!!   (style is a numerical switch; needs experiments)

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     style - integer.
*     style specifies the precise form of interpolation used in the
*     computation of flcb and cvcb.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2nxfc_style', style)
       call xertst (0.le.style.and.style.lt.2,
     &   'faulty internal parameter style')
       call ipgeti ('b2news_ncallout', ncallout)
       call ipgeti ('b2nxfc_iout', iout)                             !srv 17.06.02
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       if (ncall.eq.ncallout) then
        write(*,*) 'ncall',ncall
!WG_TODO        call hbout('b2nxfc_xvsg ',nx,ny,2,rob,'rob',ub,'ub',cvsb,
!WG_TODO     &  'cvsb',cvsb,'ua',cvsb,'na',cvsb,'ua')
       endif
*    ..test sign of rob
       call b2xvsg (nCv, rob, 1, 'rob', '.gt.')
      endif

* ..compute convection factors
      flub0 = flub
      if (style.eq.0) then
*     ..compute wrkf
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ub*rob,wrkf)
*     ..compute flcbx, cvcbx
        flub0(:,0) = flub0(:,0) + geo%fcPbshz*wrkf
      elseif (style.eq.1) then
*     ..compute wrkf
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ub,wrkf)
*     ..compute flcbx, cvcbx        
        do iFc = 1, nFc
          flub0(iFc,0) = flub0(iFc,0) + geo%fcPbs(iFc)*wrkf(iFc)*
     &      (rob(mpg%fcCv(iFc,1))+rob(mpg%fcCv(iFc,2)))*0.5e0_R8
        enddo
      else
        call xerrab ('Only style 0 or 1 in b2nxfc are included')
      endif
* ..compute linearization
!WG_TODO: change the mode "0" to an input to select different
!WG_TODO: discretization schemes
      call calccoef(nCv,nFc,nVx,0,geo,mpg,flub0,cvsb,flcb,cvcb)

      if (iout.ne.0) then                                              !sv 17.06.02 {
        write (charns,'(i3.3)') isb    
        call my_out_us(70,nFc,flub(1,0),'b2nxfc_flubx'//charns)
        call my_out_us(70,nFc,flub(1,1),'b2nxfc_fluby'//charns)
        call my_out_us(70,nFc,flcb(1,0),'b2nxfc_flcbx'//charns)
        call my_out_us(70,nFc,flcb(1,1),'b2nxfc_flcby'//charns)
        call my_out_us(70,nFc,cvcb(1,0),'b2nxfc_cvcbx'//charns)
        call my_out_us(70,nFc,cvcb(1,1),'b2nxfc_cvcby'//charns)
        call my_out_us(70,nCv,rob,'b2nxfc_rob'//charns)
      endif                                                            !sv 17.06.02 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfc

      end subroutine b2nxfc
