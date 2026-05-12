*-----------------------------------------------------------------------
*.specification

      subroutine b2upht_mol (nCv, nFc, nVx, switch, geo, mpg,
     &  ua_solve, rxf, cortm, pccm0, nm0, um0, tm0)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_fluid_molecules
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxf, cortm(nCv), pccm0(nCv)
      logical :: ua_solve
*   ..input/output arguments
      real (kind=R8) ::
     *  nm0(nCv), um0(nCv), tm0(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPHT_MOL adjusts the state of the fluid molecules towards the solution
*     of the heat equations.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, ireg
      integer, save :: ncall_b2upht_mol = 0
      real (kind=R8) ::
     *  dtm(nCv), t0,
     *  wrk(nCv), wrk0(nCv), wrkv(nVx)
      real (kind=R8) :: min_tm, max_tm, new_tm
*   ..procedures
      external xertst, b2xvsg, xerrab
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upht_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     & 'faulty argument rxf')
      call xertst (0.0_R8.le.switch%b2npht_pcm0, 'faulty argument pcm0')
*   ..test sign of vol, tm
      if(ncall_b2upht_mol.lt.3) then
        call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
        call b2xvsg (nCv, tm0, 1, 'tm', '.gt.')
      endif

* ..apply the correction
*   ..adjust tm
      min_tm = 0.02_R8*ev !! TODO make user-defined inputs
      max_tm = 10000.0_R8*ev

      dtm = rxf*cortm
      do iCv = 1, nCv
        new_tm = tm0(iCv) + dtm(iCv)
        new_tm = min(max(new_tm,min_tm),max_tm)
        tm0(iCv) = new_tm
      enddo

*   ..adjust ua
      wrk0 = nm0*dtm
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk)
      if (ua_solve) then
       do iCv = 1, mpg%nCi
        ireg = mpg%cvReg(iCv)
        t0 = switch%b2npht_pcm0*pccm0(iCv)*wrk(iCv)
         um0(iCv) = um0(iCv)-t0
       enddo
      endif

* ..return
      ncall_b2upht_mol = ncall_b2upht_mol+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht_mol

      end subroutine b2upht_mol
