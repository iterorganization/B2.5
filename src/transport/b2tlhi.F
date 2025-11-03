      subroutine b2tlhi (nCv, nFc, mp, cflmi, switch, geo, mpg,
     &  ni, nirm, ni_ext, ti, chvimx)
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      integer, intent (in) :: nCv, nFc
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  mp, cflmi, ni(nCv,0:1), nirm(nCv), ni_ext(nCv,0:1), ti(nCv),
     *  chvimx(nFc)
*     ------------------------------------------------------------------
*     B2TLHI computes the all atom heat transport flux limit velocity.
*     ------------------------------------------------------------------
*   ..local variables
      integer iFc
      real (kind=R8) ::
     *   t0, t1, cflmi_loc(nFc), nif(nFc)
*   ..procedures
      intrinsic max, sqrt
      external xertst
*   ..initialisation
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2tlhi')
*   ..test mp, cflmi
      call xertst (0.0_R8.lt.mp, 'faulty argument mp')
      call xertst (0.0_R8.le.cflmi, 'faulty argument cflmi')

      ! Init local flux limit to constant value
      cflmi_loc = cflmi

      ! Enhanced flux limit for low-density regions
      if (switch%b2tlhi_far_sol .eq. 1) then
        ! Attempt enhanced flux limit for low-density regions
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,ni(:,0),nif)
        do iFc = 1, nFc
          if (nif(iFc) .le. switch%b2tlhi_ni_min) then
           ! Set at min value
           cflmi_loc(iFc) = switch%b2tlhi_cflmi_min
          elseif (nif(iFc) .lt. switch%b2tlhi_ni_max ) then
           ! Interpolate for smooth transition
           cflmi_loc(iFc) =
     .      (cflmi*(nif(iFc) - switch%b2tlhi_ni_min) +
     .       switch%b2tlhi_cflmi_min*(switch%b2tlhi_ni_max - nif(iFc)))/
     .      (switch%b2tlhi_ni_max - switch%b2tlhi_ni_min)
          endif
        enddo
      endif

*   ..compute flux limit velocity
      if (switch%b2tqca_model.eq.1) then  ! SOLPS5.0 treatments
        do iFc = 1, nFc
cwdk perhaps the if-test should be based on a 'mode' switch, instead of b2tqca_model,
c    because now both options are interchangeable
          t0 = max(ti(mpg%fcCv(iFc,1)),ti(mpg%fcCv(iFc,2)))
          chvimx(iFc) = cflmi_loc(iFc)*abs(geo%fcPbs(iFc))*sqrt(t0/mp)*!srv 11.01.13
     &      (nirm(mpg%fcCv(iFc,1))+nirm(mpg%fcCv(iFc,2)))/
     &      (ni(mpg%fcCv(iFc,1),0)+ni_ext(mpg%fcCv(iFc,1),0)+
     &       ni(mpg%fcCv(iFc,2),0)+ni_ext(mpg%fcCv(iFc,2),0))
        enddo
      else if (switch%b2tqca_model.eq.2) then  ! SOLPS4.0 treatment
       do iFc = 1, nFc
        t0 = ti(mpg%fcCv(iFc,2))
        t1 = ti(mpg%fcCv(iFc,1))
cwdk 10.03.16 Divide by ni*ti at cell faces, so chvimx contains the same information
c           as case b2tqca_model.eq.1
c           Removes need for different treatment later on in the code (b2tfhi), and
c           ensures consistency with b2trcl_conductive_limit .ne. 0
        chvimx(iFc) = cflmi_loc(iFc)*
     &    geo%fcS(iFc)*abs(geo%fcQalf(iFc,0))*
     &   (nirm(mpg%fcCv(iFc,2))*pit(mpg%fcCv(iFc,2))*sqrt(t0/mp)*t0 +
     &    nirm(mpg%fcCv(iFc,1))*pit(mpg%fcCv(iFc,1))*sqrt(t1/mp)*t1)/
     &   ((ni(mpg%fcCv(iFc,1),0)+ni_ext(mpg%fcCv(iFc,1),0)+
     &     ni(mpg%fcCv(iFc,2),0)+ni_ext(mpg%fcCv(iFc,2),0))*max(t0,t1))
       enddo
      endif
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

      end subroutine b2tlhi
