*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfc (nCv, nFc, nVx, isb, switch, geo, mpg, rob,
     &  ub, flub, cvsb, flcb, cvcb, wrkf)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2nxfc
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) ::
     *  rob(nCv), ub(nCv),
     &  flub(nFc,0:1), cvsb(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   flcb(nFc,0:1), cvcb(nFc,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8), intent(inout) ::
     *   wrkf(nFc)                                                  !srv 02.07.08

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
      integer iFc
      character charns*3                                            !srv 17.06.02
*   ..procedures
      external xertst, ipgeti, ipgetr
      external b2xvsg, hbout, intfaceh
*   ..initialisation
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
*   ..test nCv, nFc, nVx
      call xertst (0.lt.nCv.and.0.lt.nFc.and.0.lt.nVx,
     . 'faulty argument nCv, nFc, nVx')
*   ..extensive tests on first few calls
      if (ncall_b2nxfc.lt.3) then
       if (ncall_b2nxfc.eq.switch%b2news_ncallout) then
        write(*,*) 'ncall',ncall_b2nxfc
!WG_TODO        call hbout('b2nxfc_xvsg ',nx,ny,2,rob,'rob',ub,'ub',cvsb,
!WG_TODO     &  'cvsb',cvsb,'ua',cvsb,'na',cvsb,'ua')
       endif
*    ..test sign of rob
       call b2xvsg (nCv, rob, 1, 'rob', '.gt.')
      endif

*   ..compute flcbx, cvcbx
      if (switch%b2nxfc_style.eq.0) then
*     ..compute wrkf
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ub*rob,wrkf)
*     ..compute flcbx, cvcbx
        flcb(:,0) = flub(:,0) + geo%fcPbshz*wrkf
        cvcb(:,0) = cvsb(:,0)
      elseif (switch%b2nxfc_style.eq.1) then
*     ..compute wrkf
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ub,wrkf)
*     ..compute flcbx, cvcbx
        do iFc = 1, nFc
          flcb(iFc,0) = flub(iFc,0) + geo%fcPbshz(iFc)*wrkf(iFc)*
     &      (rob(mpg%fcCv(iFc,1))+rob(mpg%fcCv(iFc,2)))*0.5_R8
        enddo
        cvcb(:,0) = cvsb(:,0)
      endif
*   ..compute flcby, cvcby
      flcb(:,1) = flub(:,1)
      cvcb(:,1) = cvsb(:,1)

      if (switch%b2nxfc_iout.ne.0) then                                !sv 17.06.02 {
        if (in_parallel()) then
#ifdef _OPENMP
!$OMP MASTER
#endif
          write(*,*) 'B2NXFC OpenMP warning: no file output in ',
     &     'parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
        else
          write (charns,'(i3.3)') isb
          call my_out_us(70,nFc,1,flub(1,0),'b2nxfc_flub_th'//charns)
          call my_out_us(70,nFc,1,flub(1,1),'b2nxfc_flub_r'//charns)
          call my_out_us(70,nFc,1,flcb(1,0),'b2nxfc_flcb_th'//charns)
          call my_out_us(70,nFc,1,flcb(1,1),'b2nxfc_flcb_r'//charns)
          call my_out_us(70,nFc,1,cvcb(1,0),'b2nxfc_cvcb_th'//charns)
          call my_out_us(70,nFc,1,cvcb(1,1),'b2nxfc_cvcb_r'//charns)
          call my_out_us(70,nCv,0,rob,'b2nxfc_rob'//charns)
        endif
      endif                                                            !sv 17.06.02 }
* ..return
      ncall_b2nxfc = ncall_b2nxfc+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfc

      end subroutine b2nxfc
