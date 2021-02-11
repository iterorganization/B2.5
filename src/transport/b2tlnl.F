*-----------------------------------------------------------------------
*.specification

      subroutine b2tlnl (nCv, switch, icase, te, ti, ne, lnlam)                 !srv 01.07.09 {
      use b2mod_types
      use b2mod_constants
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tlnl
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      integer nCv, icase                                                !srv 20.09.11
      real (kind=R8) ::
     &   te(nCv), ti(nCv), ne(nCv)                                      !srv 20.09.11
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   lnlam(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TLNL computes the Coulomb logarithm according to Braginskii or 
*     Wesson formulas.
*
*  3. description (see also routine b2cdca)
*
*  4. references
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      real (kind=R8) :: lambda
      integer iCv
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, sfill
      external b2xvsg
*   ..initialization

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lnlam is the log of the plasma parameter.
*
*     Description of some local procedures:
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tlnl')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
*   ..extensive tests on first few calls
      if (ncall_b2tlnl.lt.3) then
*    ..test state
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      endif
      lambda = switch%b2trcl_lambda
* ..compute lnlam
      if (lambda.lt.0.0_R8) then
	select case(icase)
	 case(0)      ! Braginsky formula
        do iCv = 1, nCv
	    if(te(iCv)/ev .le. 50.0_R8) then
           lnlam(iCv)=max(-lambda,
     &     23.4_R8 - 1.15_R8*log10(ne(iCv)/1.0e6_R8) +
     &                            3.45_R8*log10(te(iCv)/ev))
	    else
           lnlam(iCv)=max(-lambda,
     &             25.3_R8 - 1.15_R8*log10(ne(iCv)/1.0e6_R8) + 
     &                            2.30_R8*log10(te(iCv)/ev))
	    endif
	    if(lnlam(iCv).le.0.0_R8) then
	     open(70,file='info_b2tlnl.dat')
	     write(70,"(a)")
     &     ' ix   iy   ne(iCv)      te(iCv)  lnlam(iCv)'
	     write(70,"(i5,3es13.4)")
     &       iCv,ne(iCv),te(iCv),lnlam(iCv)
	     close(70)
	     write(*,"(a)")
     &     ' ix   iy   ne(iCv)      te(iCv)  lnlam(iCv)'
	     write(*,"(i5,3es13.4)")
     &       iCv,ne(iCv),te(iCv),lnlam(iCv)
	     call xertst (lnlam(iCv).gt.0.0_R8,
     &                  'Negative Coulomb logarithm')     !srv 12.01.16
	    endif
        enddo
	 case(1)     ! Wesson formula for ion-ion collisions
        do iCv = 1, nCv
           lnlam(iCv)=max(-lambda,17.3_R8-
     &       log(ne(iCv)/1.0e20_R8)/2.0_R8+
     &       1.5_R8*log(ti(iCv)/ev/1000.0_R8))
        enddo
	 case(2)    ! Wesson formula for electron-electron collisions
        do iCv = 1, nCv
           lnlam(iCv)=max(-lambda,14.9_R8-
     &	   log(ne(iCv)/1.0e20_R8)/2.0_R8+
     &	   log(te(iCv)/ev/1000.0_R8))
        enddo
	 case(3)   !  Wesson formula for electron-ion collisions
        do iCv = 1, nCv
           lnlam(iCv)=max(-lambda,15.2-
     &	  log(ne(iCv)/1e20_R8)/2.0_R8+
     &	  log(te(iCv)/ev/1000.0_R8))
        enddo
	 case default
	  write(*,*) ' b2tlnl: incorrect value icase =',icase
	  call xerrab ('Incorrect icase value in b2tlnl')
	 end select
      else
        lnlam = lambda
      endif
* ..return
      ncall_b2tlnl = ncall_b2tlnl+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tlnl

      end                                                            !srv 01.07.09 }
