*-----------------------------------------------------------------------
*.specification

      subroutine b2upht_mol_3D (nCv, nFc, nVx, switch, geo, mpg,
     &  ua_solve, rxf, cortm, pccm0_pol, pccm0_rad,
     &  nm0, um0_pol, um0_rad, tm0)
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
     *  rxf, cortm(nCv), pccm0_pol(nCv), pccm0_rad(nCv)
      logical :: ua_solve
*   ..input/output arguments
      real (kind=R8) ::
     *  nm0(nCv), um0_pol(nCv), um0_rad(nCv), tm0(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPHT_MOL_3D adjusts the state of the fluid molecules towards the solution
*     of the heat equations.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, ireg
      integer, save :: ncall_b2upht_mol = 0
      real (kind=R8) ::
     *  dtm(nCv), t0,
     *  wrk(nCv,0:1), wrk0(nCv), wrkv(nVx,0:1)
      real (kind=R8) :: min_tm, max_tm, new_tm
*   ..procedures
      external xertst, b2xvsg, xerrab
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upht_mol_3D')
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
      min_tm = 0.02_R8*ev !! FM_TODO make user-defined inputs
      max_tm = 10000.0_R8*ev

      dtm = rxf*cortm
      do iCv = 1, nCv
        new_tm = tm0(iCv) + dtm(iCv)
        new_tm = min(max(new_tm,min_tm),max_tm)
        tm0(iCv) = new_tm
      enddo

*   ..adjust ua
      wrk0 = nm0*dtm
      call gradc(nCv, nFc, nVx, 0, geo, mpg, wrk0, wrkv, wrk)
      if (ua_solve) then
       do iCv = 1, mpg%nCi
        ireg = mpg%cvReg(iCv)
        t0 = switch%b2npht_pcm0*pccm0_pol(iCv)*wrk(iCv,0)
        um0_pol(iCv) = um0_pol(iCv)-t0
        t0 = switch%b2npht_pcm0*pccm0_pol(iCv)*wrk(iCv,1)
        um0_rad(iCv) = um0_rad(iCv)-t0
       enddo
      endif

* ..return
      ncall_b2upht_mol = ncall_b2upht_mol+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upht_mol

      end subroutine b2upht_mol_3D
