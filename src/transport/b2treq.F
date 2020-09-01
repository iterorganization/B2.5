*-----------------------------------------------------------------------
*.specification

      subroutine b2treq (nCv, switch, vol, te, ti, ne, ni, ne2m, lnlam,
     &   st_ext, ceqp)   !xpb !srv 01.07.09
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_plasma
      use b2mod_switches
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv
      type (switches), intent (in) :: switch
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  vol(nCv), ni(nCv,0:1), 
     &  te(nCv), ti(nCv),
     &  ne(nCv), ne2m(nCv), lnlam(nCv)        !xpb !srv 01.07.09
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  ceqp(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TREQ computes the cell-integrated energy equipartition
*     multiplier.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the electron-atom energy equipartition
*     multiplier ceqp, defined so that on each cell the signed
*     equipartition energy transfer from atoms to electrons is
*     (ceqp*ne*ni/te**1.5)*(te-ti) .
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv
      real (kind=R8) ::
     *  cfeqp
      parameter (cfeqp=4.8e-15_R8)
      save ncall
*     lnlam is the log of the plasma parameter.
*     cfeqp is a numerical coefficient in the energy transfer rate.
*   ..procedures
      intrinsic sqrt, log, max
      external subini, subend, xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2treq')
*   ..test nCv
      call xertst (0.le.nCv, 'faulty argument nCv')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test vol, state
       call b2xvsg (nCv, vol, 1, 'vol', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, ne2m, 1, 'ne2m', '.gt.')
      endif

      call b2tlnl (nCv, switch, switch%icase_ei, te, ti, ne, lnlam)                     !srv 20.09.11
*   ..compute ceqp
      do iCv = 1, nCv
        ceqp(iCv) = switch%b2treq_phm0*
     &    (cfeqp*lnlam(iCv)*ev*sqrt(ev))*vol(iCv)*                               !srv 01.07.09
     &    (ne2m(iCv)/(ni(iCv,0)+st_ext%ni(iCv,0)))
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2treq

      end subroutine b2treq

