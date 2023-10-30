      subroutine b2tlhe (nCv, nFc, me, cflme, switch, geo, mpg, ne, te,
     &  chvemx)
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tlhe
      use b2mod_subsys
      implicit none
      integer nCv, nFc
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  me, cflme, ne(nCv), te(nCv), chvemx(nFc)
*     ------------------------------------------------------------------
*     B2TLHE computes the electron heat transport flux limit velocity.
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc
      real (kind=R8) ::
     *   t0, t1
*   ..procedures
      intrinsic max, sqrt
      external xertst, b2xvsg

*   ..initialisation
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlhe')
*   ..test nCv, nFc
      call xertst (0.le.nCv .and. 0.le.nFc, 'faulty argument nCv, nFc')
*   ..test me, cflme
      call xertst (0.0_R8.lt.me, 'faulty argument me')
      call xertst (0.0_R8.le.cflme, 'faulty argument cflme')
*   ..extensive tests on first few calls
      if (ncall_b2tlhe.lt.3) then
*    ..test sign of ne, te
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
      endif

*   ..compute flux limit velocity
      if (switch%b2tqce_model.eq.1.or.switch%b2tqce_model.eq.3) then  ! SOLPS5.0 treatments
cwdk perhaps the if-test should be based on a 'mode' switch, instead of b2tqce_model,
c    because now both options are interchangeable
        do iFc = 1, nFc
          t0 = max(te(mpg%fcCv(iFc,1)),te(mpg%fcCv(iFc,2)))
          chvemx(iFc) = cflme*sqrt(t0/me)*abs(geo%fcPbs(iFc))            !srv 11.01.13
        enddo
      else if (switch%b2tqce_model.eq.2) then  ! SOLPS4.0 treatment
        do iFc = 1, nFc
          t0 = te(mpg%fcCv(iFc,2))
          t1 = te(mpg%fcCv(iFc,1))
cwdk 10.03.16 Divide by ne*te at cell faces, so chvemx contains the same information
c           as case b2tqce_model.eq.1
c           Removes need for different treatment later on in the code (b2tfhe), and
c           ensures consistency with b2trcl_conductive_limit .ne. 0
          chvemx(iFc) = cflme * geo%fcS(iFc) * abs(geo%fcQalf(iFc,0)) *
     .       (ne(mpg%fcCv(iFc,2))*pit(mpg%fcCv(iFc,2))*sqrt(t0/me)*t0 +
     .        ne(mpg%fcCv(iFc,1))*pit(mpg%fcCv(iFc,1))*sqrt(t1/me)*t1)/
     .       ((ne(mpg%fcCv(iFc,1))+ne(mpg%fcCv(iFc,2)))*max(t0,t1))
        enddo
      endif
      ncall_b2tlhe = ncall_b2tlhe+1
      call subend ()
      return
*     ------------------------------------------------------------------
      contains

      function pit(iCv)
      implicit none
      real (kind=R8) :: pit
      integer, intent(in) :: iCv

      pit = abs(geo%cvBb(iCv,0))/geo%cvBb(iCv,3)

      return
      end function pit

      end subroutine b2tlhe
