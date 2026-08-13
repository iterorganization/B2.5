*-----------------------------------------------------------------------
*.specification

      subroutine b2smgp_perp (nCv, nFc, nVx, switch, geo, mpg,
     &  nm, tm, smgp_rad,smgp_dia)
      use b2mod_types
      use b2mod_constants                                                !srv 23.03.10
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, nFc, nVx
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) :: nm(nCv), tm(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   smgp_rad(nCv), smgp_dia(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SMGP_PERP computes the perpendicular pressure gradient term 
*     on the right hand of the perpendicular momentum balance equation
*     for fluid molecules.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc
      character chns*3
      real (kind=R8) ::
     *   t0, wrkf(nFc), wrkv(nVx), nmf(nFc), tmf(nFc),
     *   wrk0(nCv), wrk1(nCv), wrk(nCv)
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2smgp_perp')
*   ..test nCv, nFv
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

!*   .. check output writing
!      if ((switch%b2sigp_iout.ne.0.or.switch%iout_b2wdat.eq.4)) then                       !srv 17.06.02 {
!        if (in_parallel() .and. ncall_b2sigp == 0) then
!#ifdef _OPENMP
!$OMP MASTER
!#endif
!           if (ncall_b2sigp==0) then
!             write(*,*) 'b2sigp OpenMP warning: no file output ',
!     &                  'in parallel mode'
!           end if
!#ifdef _OPENMP
!$OMP END MASTER
!#endif
!        end if
!      end if

* ..interpolate some variables to cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,nm,nmf)
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,tm,tmf)

* ..compute pressure gradient term
*   ..compute smbgp
*   ..  compute gradients in cell centers
          wrk = nm*tm
          call gradc_r(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkv, wrk0)
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkv, wrk1)
* ..compute pressure gradient term
          smgp_rad = 0.0_R8
          smgp_dia = 0.0_R8
*   ..  compute smbgp
          do iCv = 1, mpg%nCi
           t0=geo%cvVol(iCv)
           smgp_rad(iCv) = -t0*wrk0(iCv)
           smgp_dia(iCv) = -t0*wrk1(iCv) !! TODO_FM poloidal vs diamagnetic !!!
          enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2smgp_perp
      end subroutine b2smgp_perp
