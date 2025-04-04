*-----------------------------------------------------------------------
*.specification
c**************************************************************************************
c**********************  NEW KU Leuven module for ECH plasmas  ************************
c**************************************************************************************
! ap  01/09/23
! main author: Anthony Piras
      module b2mod_ech
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_math
      use b2mod_ad
     & , only : ncall_calcECH

      implicit none
*   ..declarations allocatable variable
      real (kind=R8), allocatable, save ::
     &  chi(:), she_ech(:,:)
      integer, allocatable, save :: list_ci_ech(:)

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ECH is intended to simulate power depositions in Electron
*     Cyclotron Heated plasmas. To turn on the model, the flag
*     'b2sral_use_ECH' has to be set to '1'. This module computes an
*     additional heating source due to EC heating injected either at the
*     cyclotron frequency or its harmonics.
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      contains

c alloc
      subroutine alloc_b2mod_ech(mpg, nCv)
      implicit none
      type(mapping), intent(in):: mpg
      integer, intent(in) :: nCv

      allocate(chi(mpg%nCi))
      allocate(she_ech(nCv,0:3))
      allocate(list_ci_ech(mpg%nCi))

      chi(:) = 0.0_R8
      she_ech(:,:) = 0.0_R8
      list_ci_ech = 0

      return
      end subroutine alloc_b2mod_ech
c end alloc

c ECH power computation
      subroutine CalcECH(nCv,geo,mpg,dv,switch)
c     This routine calculates the ECH power deposition
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_ad
     & , only : ncall_calcECH
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c     Output variables
c      real (kind=R8), intent(out):: she_ech(nCv,0:3)
      character chns*3


* ..preliminaries
*   ..subprogram start-up calls
      call subini ('CalcECH')
c      call sfill (nCv*4, 0.0_R8, she0_ech, 1) not used anymore

      if(ncall_calcECH.lt.1) then
        write (0,*) 'ECH module used.'
        write (*,*) 'ECH module used.'
      endif

      select case (switch%ech_coupling)

* -- COUPLING = 0 -- Fixed power distributed as a Dirac delta in volumes intersected at given radial position.

      case(0)
        call compute_ech_0(nCv,geo,mpg,switch)


* -- COUPLING = 1 -- Fixed fraction of launched power is coupled with a Gaussian-like distribution.

      case(1)
        call compute_ech_1(nCv,geo,mpg,switch)


* -- COUPLING = 2 -- Budden analysis ! TODO
      case(2)
        call compute_ech_2(nCv,geo,mpg,dv,switch)


* -- COUPLING = 3 -- Ray-tracing. ! TODO

      case(3)
        call xerrab ('b2mod_ech -- ECH_COUPLING = 3 not coded yet')

      end select

      if (switch%ech_out.eq.1.or.switch%iout_b2wdat.eq.4) then
c        call my_out_us(70,nCv,0,chi,'b2mod_ECH_shapefun')
c        call my_outi_us(70,nCv,0,maskf,'b2mod_ECH_maskf')
c        do cntr = 0,3
c          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,she_ech(:,0),'b2mod_ECH_she_ech_0')
c        enddo
      endif

      ncall_calcECH = ncall_calcECH + 1

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end CalcECH

      end subroutine CalcECH


c Compute ECH - Dirac's delta
      subroutine compute_ech_0(nCv,geo,mpg,switch)
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_ad
     & , only : ncall_calcECH
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c Output variable
c      real (kind=R8), intent(out) :: she_ech(nCv,0:3)
c Internal variables
      real (kind=R8) :: segment(2,2), segmdist(mpg%nCi),
     &                  totaldist, crction, summation
      integer :: list_cv(nCv), nncv, nnci,
     &                 list_fc(mpg%nFc), nnfc
      logical :: ifail
      integer :: iCv, cntr

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('compute_ech_0')
      call sfill (nCv*4, 0.0_R8, she_ech, 1)

* ..precomputation
      if(ncall_calcECH.lt.1) then
        write (0,*) 'Dirac delta coupling scheme used.'
        write (*,*) 'Dirac delta coupling scheme used.'
        summation = 0.0_R8
        totaldist = 0.0_R8
        crction = switch%geom_match_dist !ap correction distance to avoid numerical error when defining segment for cv_intersection
        call sfill (mpg%nCi, 0.0_R8, chi, 1)
        segment = reshape((/ switch%ech_Rdep, minval(geo%vxY)-crction,
     &      switch%ech_Rdep, maxval(geo%vxY)+crction /), shape(segment))

        call calc_face_dist(geo, mpg, nCv, segment, list_cv, nncv,
     &   list_ci_ech, nnci, list_fc, nnfc, segmdist, totaldist, ifail)

        chi(:) = segmdist/totaldist
        summation = sum(chi(:)*geo%cvVol(:)) ! normalisation condition
        chi = chi/summation
      endif

* ..computation
      do cntr = 1, size(list_ci_ech)
        if(list_ci_ech(cntr).ne.0) then
          iCv = list_ci_ech(cntr)
          she_ech(iCv,0) = chi(iCv)*switch%ech_absfract*
     &                     switch%ech_power*geo%cvVol(iCv)
        endif
      enddo

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end compute_ech_0

      end subroutine compute_ech_0


c Compute ECH - Gaussian-like
      subroutine compute_ech_1(nCv,geo,mpg,switch)
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_ad
     & , only : ncall_calcECH
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c Output variable
c      real (kind=R8), intent(out) :: she_ech(nCv,0:3)
c Internal variables
      real (kind=R8):: sx, sy, summation
      integer :: iCv

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('compute_ech_1')
      call sfill (nCv*4, 0.0_R8, she_ech, 1)
      summation = 0.0_R8

* ..precomputation

      if(ncall_calcECH.lt.1) then
        write (0,*) 'Gaussian-like coupling scheme used.'
        write (*,*) 'Gaussian-like coupling scheme used.'
        sx = switch%ech_gaussian_width_x ! widthech (X direction) in TOMATOR 1D
        sy = switch%ech_gaussian_width_y ! widthech (Y direction)
        call sfill (mpg%nCi, 0.0_R8, chi, 1)
        if (switch%ech_gaussian_type.eq.1) then ! gaussian symmetric in poloidal direction
          write (0,*) 'Poloidally symmetric gaussian used. ' //
     &                'Power spread over the whole domain in '
          write (0,*) 'poloidal direction. Should be used for '
     &                // 'symmetric cases only!'
          write (*,*) 'Poloidally symmetric gaussian used. ' //
     &                'Power spread over the whole domain in '
          write (*,*)  'poloidal direction. Should be used for '
     &                // 'symmetric cases only!'
          do iCv = 1, mpg%nCi
            chi(iCv) = geo%cvX(iCv)*expu(-
     &            ((geo%cvX(iCv)-switch%ech_Rdep)/sx)**2)+
     &            1.0E-7*geo%cvX(iCv)*expu(- ! 1.0E-7 is the background value (from TOMATOR 1D)
     &            ((geo%cvX(iCv)-switch%ech_Rdep)/switch%a0)**2)
          enddo
        elseif (switch%ech_gaussian_type.eq.2) then ! standard binormal distribution
          write (0,*) 'Standard binormal distribution used. ' //
     &                'No toroidal effects considered.'
          write (*,*) 'Standard binormal distribution used. ' //
     &                'No toroidal effects considered.'
          do iCv = 1, mpg%nCi
            chi(iCv) = 1.0_R8/(2*pi*sx*sy)*exp(-1.0_R8/2.0_R8*
     &                 ((geo%cvX(iCv)-switch%ech_Rdep)**2/sx**2 +
     &                  (geo%cvY(iCv)-switch%ech_Zdep)**2/sy**2))
          enddo
        else ! 1/R*binormal to account for toroidal effects (see A. Piras PhD)
          write (0,*) 'Standard binormal distribution used. ' //
     &                'Toroidal effects considered.'
          write (*,*) 'Standard binormal distribution used. ' //
     &                'Toroidal effects considered.'
          do iCv = 1, mpg%nCi
            chi(iCv) = 1.0_R8/geo%cvX(iCv)*(1.0_R8/(2*pi*sx*sy)*
     &                   exp(-1.0_R8/2.0_R8*((geo%cvX(iCv)-
     &                   switch%ech_Rdep)**2/sx**2 + (geo%cvY(iCv)
     &                   -switch%ech_Zdep)**2/sy**2)))
          enddo
        endif

        summation = sum(chi(:)*geo%cvVol(:))
        chi = chi/summation
      endif

      do iCv = 1, mpg%nCi
        she_ech(iCv,0) = switch%ech_power*chi(iCv)*
     &                   switch%ech_absfract*geo%cvVol(iCv)
      enddo

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end compute_ech_1

      end subroutine compute_ech_1


c Compute ECH - Budden analysis !TODO
      subroutine compute_ech_2(nCv,geo,mpg,dv,switch)
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_ad
     & , only : ncall_calcECH
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c Output variable
c      real (kind=R8), intent(out) :: she_ech(nCv,0:3)
c Internal variables
      real (kind=R8):: segment(2,2), crction
      real (kind=R8) :: wce, wceV(nCv), wpe(nCv),
     &                  wRcut(nCv), wLcut(nCv), wUHR(nCv), alphaCE(nCv),
     &                  alphaUHR(nCv)
      real (kind=R8) :: XI(mpg%nFc), YI(mpg%nFc),
     &                  neRcut(nCv), neLcut(nCv), neUHR(nCv),
     &                  R_UHR(nCv), BUHR(nCv)
      real (kind=R8) :: gradne(nCv), gradB(nCv), funv(mpg%nVx),
     &                  Lne(nCv), LBphi(nCv)
c Other variables
      real (kind=R8) :: omegace(nCv), R_res
c Internal variables
      integer :: list_cv(nCv), nncv, list_ci(mpg%nCi), nnci,
     &           list_fc(mpg%nFc), nnfc
      integer, allocatable :: list_cutoff(:)
      integer :: iCv, iCv2, iFc, cntr
      logical :: ifail

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('compute_ech_2')
      call sfill (nCv*4, 0.0_R8, she_ech, 1)

* ..precomputation
      crction = switch%geom_match_dist !ap correction distance to avoid numerical error when defining segment for cv_intersection

c-----------------Locate resonance-----------------
        call find_resonance(nCv, geo, mpg, switch, omegace, wce, R_res)

        segment = reshape((/ R_res, minval(geo%cvY)-crction,
     &       R_res, maxval(geo%cvY)+crction /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv, list_ci,
     &   list_fc, nncv, nnci, nnfc, nCv, XI, YI, 1, ifail)

        call sfill (nCv, wce, wceV, 1)

        do iCv = 1, nncv
          iCv2 = list_cv(iCv)
          wpe(iCv2) = sqrt(dv%ne(iCv2)*qe**2/me/eps0)
        enddo
        alphaCE = wpe**2/wce**2

c-----------------Locate R-cutoff-----------------

        wRcut = 0.5_R8*omegace + 0.5_R8*sqrt(4.0_R8*dv%ne*
     &          qe*qe/me/eps0 + omegace**2) ! frequency for UH resonance

        do iCv = 1, nCv
          if (wRcut(iCv).le.wceV(iCv)) then
            cntr = cntr + 1
          endif
        enddo

        allocate(list_cutoff(cntr))

        iCv2 = 0
        do iCv = 1, nCv
          if (wRcut(iCv).le.wceV(iCv)) then
            iCv2 = iCv2 + 1
            list_cutoff(iCv2) = iCv
          endif
        enddo
        
        neRcut(list_cutoff) = dv%ne(list_cutoff)
        deallocate (list_cutoff)
        cntr = 0

c-----------------Locate L-cutoff-----------------

        wLcut = -0.5_R8*omegace + 0.5_R8*sqrt(4.0_R8*dv%ne*
     &          qe*qe/me/eps0 + omegace**2) ! frequency for UH resonance

        do iCv = 1, nCv
          if (wLcut(iCv).le.wceV(iCv)) then
            cntr = cntr + 1
          endif
        enddo

        allocate(list_cutoff(cntr))

        iCv2 = 0
        do iCv = 1, nCv
          if (wLcut(iCv).le.wceV(iCv)) then
            iCv2 = iCv2 + 1
            list_cutoff(iCv2) = iCv
          endif
        enddo
        
        neLcut(list_cutoff) = dv%ne(list_cutoff)
        deallocate (list_cutoff)
        cntr = 0

c-----------------Locate UHR-layer-----------------

        wUHR = sqrt(dv%ne*qe*qe/me/eps0 + omegace**2) ! frequency for UH resonance

        do iCv = 1, nCv
          if (wUHR(iCv).le.wceV(iCv)) then
            cntr = cntr + 1
          endif
        enddo

        allocate(list_cutoff(cntr))

        iCv2 = 0
        do iCv = 1, nCv
          if (wUHR(iCv).le.wceV(iCv)) then
            iCv2 = iCv2 + 1
            list_cutoff(iCv2) = iCv
          endif
        enddo

        neUHR(list_cutoff) = dv%ne(list_cutoff)
        BUHR(list_cutoff) = geo%cvBb(list_cutoff,3)
        R_UHR(list_cutoff) = geo%cvX(list_cutoff) ! UHR radial location
        deallocate (list_cutoff)
        cntr = 0

c-----------------Budden Analysis @ UHR-----------------

c     Lne = abs(nr.ne[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (nr.ne[id - 1] - nr.ne[id + 1])); ! density scale length [m]
        call gradc_r(nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     &               neUHR, funv, gradne)

c     LBr = abs(Br[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (Br[id - 1] - Br[id + 1])); ! magnetic scal length [m]
        call gradc_r(nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     &               BUHR, funv, gradB)

        Lne = neUHR/gradne
        LBphi = BUHR/gradB
        alphaUHR = sqrt(neUHR*qe**2/me/eps0) / omegace ! wpe/wce @ UHR
#ifdef WG_TODO
        Buddenparam = omegace*Lne/c * alphaUHR/sqrt(alphaUHR**2 +
     &                2.0_R8*(Lne/LBphi))*((sqrt(1.0 +
     &		      alphaUHR**2)-1.0_R8)/(alphaUHR**2+(Lne/LBr)*
     &                sqrt(1.0 + alphaUHR**2))**0.5_R8 ! Budden parameter

        T = expu(-2.0_R8 * Buddenparam) ! Tunneling fraction
        R = (1.0_R8 - Tunnel)**2        ! Reflected part
        C = 1.0_R8 - Tunnel - Reflect   ! XB converted fraction

c-----------------O-absorption-----------------

      DensO = res_freq**2*me*eps0/(qe**2)
      if (maxval(dv%ne).le.DensO) then
        PabsO = (1.0-exp(-0.25*pi*R/100.0/c/c/c*Teabsce * 11604.45*kb / me * wpe * wpe / wce * sqrt(1.0 - alphaCE))); ! 1 pass aborption for O-wave
        PabsO = PabsO / (1.0 - (1.0 - muw) * (1.0 - PabsO)); ! multi-pass absorption with losses to the wall (muw)
      endif

c-----------------X-absorption-----------------

     PabsX = (1.0 - exp(-0.5 * pi * R / 100.0 / c / c / c * Teabsce * 11604.45 * kb / me * wpe * wpe / wce * pow(2.0 - alphaCE, 2.5) * (1.0 +
		   sqrt(alphaCE)))); ! 1 pass aborption for X-wave below Rcutoff
    DensRcut = 3.13e9;
    Densuh = 6.26e9;
    if (nemax < DensRcut) {
        PabsX = PabsX * abs(DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE - 1.0e-5) / (DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE +
		1.0e-5); ! adaptation for X-wave nearing Rcutoff
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX)); ! multi-pass absorption with losses to the wall (muw)
    }
    else if (nemax < Densuh) {
        PabsX = (1.0-Reflect)*PabsX/(1.0 - (1.0 - muw) * (1.0 - PabsX)); ! multi-pass absorption of the tunneling wave with losses to the wall (muw)
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX)); ! multi-pass absorption with losses to the wall (muw)
    }
    ! Not reflected part tunnels and reaches ECR,, reflecting on wall and coupling at second pass
    else if (nemax < DensLcut) {
        PabsX = (1.0 - PRdep2) * Tunnel * PabsX + (1.0 - PRdep4) * Tunnel * (1.0 - PabsX) * (1.0 - muw) * PabsX;
    }
    ! X-wave tunneling through R-cutoff, reflecting on wall and coupling at second pass
    ! third pass not possible because of evanescent region, SX-B conversion possible
    else {
        PabsX = 0.0;
    } ! All X-wave is converted to B-wave

*-----------------B-absorption-----------------

    double PabsB = 0.0; ! EBW coupling near ECR
    id = 0;
    while (Br[id] > Bce)
        ++id;
    if ((abs(Br[id] - Bce)) > (abs(Br[id + 1] - Bce))) {
        ++id;
    }
    !no shift to LFS of ECR implemented for the moment, doppler shift to be added ????
     RdepB = Rdep3; ! Hardcoded position of B-coupling

    if (nemax < Densuh) {
        PabsB = 0.0;
    }
    ! No B-wave present

    else if (nemax < DensO) {
        PabsB = PRdep2 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE)); ! percentage of converted X-wave to B-wave coupling
    }
    ! B-wave appears when UHR appears
    else {
        PabsB = PRdep3 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE));
        PabsB += PRdep3 * (1.0 - PRdep4) * (1.0 - PRdep4) * Tunnel;
    }
    ! double conversion to B-wave present: 1st converted at UHR, 2nd from tunneling X-wave reflected at Lcutoff and converted to B-wave at UHR

*-----------------UH-resonance-----------------

    double PabsU = 0.0; ! coupling at the UHR

    if (nemax < Densuh) {
        PabsU = 0.0;
    }
    ! No UHR present
    else if (nemax < DensO) {
        PabsU = PRdep4 * (Tunnel); ! PRdep4 = percentage of UHR coupling of tunneled X-wave
    } else {
        PabsU = PRdep4 * (Tunnel); ! 1st resonance at UHR
        PabsU += PRdep4 * (1.0 - PRdep4) * (Tunnel); ! 2nd resonance at UHR from tunneled X-wave reflected at Lcutoff
    }

    ! #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { ! line integrated density
        eprof1[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Rdepce) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0)); ! normaalverdeling met hoogte 1 voor O/X-wave deposition at ECR
        eprof2[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - RdepB) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));  ! B-wave deposition at ECR+1.5cm (should be moving from 0 to + 2.0 ???)
        eprof3[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Ruh) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));    ! UH Resonance
    }

    for (int id = 0; id < NMESHP - 1; ++id) { ! line integrated density
        inteprof1 += 2.0 * b * pi * 0.5 * (eprof1[id] + eprof1[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof2 += 2.0 * b * pi * 0.5 * (eprof2[id] + eprof2[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof3 += 2.0 * b * pi * 0.5 * (eprof3[id] + eprof3[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    }

    for (int id = 0; id < NMESHP; ++id) { ! line integrated density
        eprof1[id] = eprof1[id] / inteprof1;
        eprof2[id] = eprof2[id] / inteprof2;
        eprof3[id] = eprof3[id] / inteprof3;
    }

    for (int id = 0; id < NMESHP; ++id) {
        Pabs1[id] = (PRdep1 * PabsO * eprof1[id]) * switch%ech_power ! Power (kW) of O-wave absorption at ECR
        Pabs2[id] = (PRdep2 * PabsX * eprof1[id]) * switch%ech_power ! Power (kW) of X-wave absorption at ECR (inluding tunneled)
        Pabs3[id] = (PRdep2 * PabsB * eprof2[id]) * switch%ech_power ! Power (kW) of B-wave absorption near ECR (converted from X-wave => PRdep2)
        Pabs4[id] = (PRdep2 * PabsU * eprof3[id]) * (switch%ech_power * 1.0e3); ! Power (kW) of UHR absorption at UHR
    }
}

        do iCv = 1, mpg%nCi
          she0_ech(iCv,0) = switch%ech_power*chi(iCv)*
     &                      switch%ech_absfract*geo%cvVol(iCv)
        enddo
#endif

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end compute_ech_2

      end subroutine compute_ech_2


c Find resonance layer
      subroutine find_resonance(nCv,geo,mpg,switch,omegace,
     &           res_freq,R_res)
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out) :: res_freq, omegace(nCv), R_res
c     Internal variables
      real (kind=R8):: B_res

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('find_resonance')
      call sfill (nCv, 0.0_R8, omegace, 1)

      call calc_elec_cycl_freq(nCv,geo,omegace)

      res_freq = 2*pi*switch%ech_RFfreq/switch%ech_harmonics
      B_res = me*res_freq/qe

      R_res = geo%cvBb(1,3)*geo%cvX(1)/B_res !ap field assumed mostly toroidal

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end find_resonance

      end subroutine find_resonance

c Electron Cyclotron Frequency computation
      subroutine calc_elec_cycl_freq(nCv,geo,omegace)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      type (geometry), intent (in) :: geo
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: omegace(nCv)

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('calc_elec_cycl_freq')
      call sfill (nCv, 0.0_R8, omegace, 1)

      omegace(:) = qe*geo%cvBb(:,3)/me

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end calc_elec_cycl_freq

      end subroutine calc_elec_cycl_freq

c dealloc
      subroutine dealloc_b2mod_ech
      implicit none

      if (allocated(she_ech)) deallocate(she_ech)
      if (allocated(chi)) deallocate(chi)
      if (allocated(list_ci_ech)) deallocate(list_ci_ech)

      return
      end subroutine dealloc_b2mod_ech
c End dealloc

*-----------------------------------------------------------------------
*.end b2mod_ech

      end module b2mod_ech
