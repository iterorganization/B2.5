*-----------------------------------------------------------------------
*.specification
c**************************************************************************************
c**********************  NEW KU Leuven module for ECH plasmas  ************************
c**************************************************************************************
! mb  01/09/23
! main author: Anthony Piras
      module b2mod_ech
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2us_plasma

      implicit none
*   ..declarations allocatable variable
      real (kind=R8), allocatable, save ::
     &  sECH(:,:)

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

c Alloc
      subroutine alloc_b2mod_ech(nCv)
      implicit none
      integer, intent(in) :: nCv

      allocate(sECH(nCv,0:3))

      sECH(:,:) = 0.0_R8

      return
      end subroutine alloc_b2mod_ech
c End alloc

c ECH power computation
      subroutine CalcECH(nCv,geo,mpg,dv,switch,she0_ech)
c     This routine calculates the ECH power deposition
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_math
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: she0_ech(nCv,0:3)
c     Other variables
      real (kind=R8) :: omegace(nCv), Rdep_vect(nCv)
      integer :: maskf(nCv)
c     Internal variables
      real (kind = R8) :: chi(nCv), sx, sy, intchi, normlz, h,
     &                    summation, res_freq
      real (kind = R8) :: segment(2,2), segmdist(nCv), totaldist
      real (kind = R8) :: fcCnX1, fcCnX2, fcCnY1, fcCnY2
      integer :: list_cv(nCv), nncv, list_ci(mpg%nCi), nnci,
     &           list_fc(mpg%nFc), nnfc, face1, face2
      integer :: iCv, iCv2, iFc, cntr
      character chns*3
      logical:: ifail

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('CalcECH')
      call sfill (nCv*4, 0.0_R8, she0_ech, 1)
      call ifill (nCv, 0, maskf, 1)
      call sfill (nCv, switch%ech_Rdep, Rdep_vect, 1)
      call sfill (nCv, 0.0_R8, chi, 1)
      call sfill (4, 0.0_R8, segment, 1)
      call sfill (nCv, 0.0_R8, segmdist, 1)
      summation = 0.0_R8
      totaldist = 0.0_R8
      fcCnX1 = 0.0_R8
      fcCnX2 = 0.0_R8
      fcCnY1 = 0.0_R8
      fcCnY2 = 0.0_R8
      cntr = 0

      select case (switch%ech_coupling)


* -- COUPLING = 0 -- Fixed power distributed as a Dirac delta in volumes intersected at given radial position.

      case(0)
        segment = reshape((/ switch%ech_Rdep, -0.1_R8,
     &       switch%ech_Rdep, maxval(geo%cvY)+0.1_R8 /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv, list_ci,
     &   list_fc, nncv, nnci, nnfc, nCv, 1, ifail)
        do iFc = 1, nnfc-1
            iCv = list_ci(iFc)
            face1 = list_fc(iFc)
            face2 = list_fc(iFc+1)
            fcCnX1 = (geo%vxX(mpg%fcVx(face1,1))+
     &                geo%vxX(mpg%fcVx(face1,2)))/2
            fcCnY1 = (geo%vxY(mpg%fcVx(face1,1))+
     &                geo%vxY(mpg%fcVx(face1,2)))/2
            fcCnX2 = (geo%vxX(mpg%fcVx(face2,1))+
     &                geo%vxX(mpg%fcVx(face2,2)))/2
            fcCnY2 = (geo%vxY(mpg%fcVx(face2,1))+
     &                geo%vxY(mpg%fcVx(face2,2)))/2
            segmdist(iCv) = sqrt((fcCnX1-fcCnX2)**2+(fcCnY1-fcCnY2)**2)
            totaldist = totaldist + segmdist(iCv)
        enddo
        chi(:) = segmdist/totaldist ! (distance of 2 face centers (vxX(fcVx(ifc,1)+vxX(fcVx(ifc,2))/2 for X, same for Y
        summation = sum(chi(:)*geo%cvVol(:)) ! normalisation condition
        chi = chi/summation
        do cntr = 1, size(list_ci)
          if(list_ci(cntr).ne.0) then
            iCv = list_ci(cntr)
            she0_ech(iCv,0) = chi(iCv)*switch%ech_absfract*
     &                        switch%ech_power*geo%cvVol(iCv)
          endif
        enddo

* -- COUPLING = 1 -- Fixed fraction of launched power is coupled with a Gaussian-like distribution.

      case(1)
        sx = switch%ech_gaussian_width ! widthech in TOMATOR 1D
        sy = sx ! assumed equal for now for 2D
        if (switch%ech_gaussian_symmetry.eq.1) then ! symmetry on poloidal direction
          do iCv = 1, mpg%nCi
            chi(iCv) = geo%cvX(iCv)*expu(-
     &            ((geo%cvX(iCv)-switch%ech_Rdep)/sx)**2)+
     &            1.0E-7*geo%cvX(iCv)*expu(-
     &            ((geo%cvX(iCv)-switch%ech_Rdep)/switch%a0)**2)
          enddo

          summation = sum(chi(:)*geo%cvVol(:))
          chi = chi/summation

        else
          chi = 1.0_R8/(2*pi*sx*sy)*exp(-1.0_R8/2.0_R8*
     &            ((geo%cvX-switch%ech_Rdep)**2/sx**2 +
     &             (geo%cvY-switch%ech_Zdep)**2/sy**2))

          do iCv = 1, mpg%nCi
            summation = summation + chi(iCv)*geo%cvVol(iCv) ! normalisation condition
          enddo
          chi = chi/summation ! normalisation condition

        endif

        do iCv = 1, mpg%nCi
          she0_ech(iCv,0) = switch%ech_power*chi(iCv)*
     &                      switch%ech_absfract*geo%cvVol(iCv)
        enddo


* -- COUPLING = 2 -- TOMAS - deposition on ECR (Budden analysis). !!TODO

      case(2)
#ifdef WG_TODO
        call find_resonance(nCv,geo,mpg,switch,omegace,res_freq,
     &       maskf_res)

        call sfill (nCv, res_freq, res_freqV, 1)

        omega_plasma = sqrt(dv%ne*qe**2/me/eps0)*maskf
        alphaCE = omega_plasma**2/res_freq**2

*-----------------Locate R-cutoff-----------------

        wRcut = 0.5_R8*geo%cvBb(:,3)*qe/me + 0.5_R8*sqrt(4.0*dv%ne*
     &          qe**2/me/eps0 + (geo%cvBb(:,3)*qe/me)**2) ! frequency for UH resonance

        call maskfun(mpg%nCi, wRcut, res_freqV,
     &              1.0E-3_R8, maskf_Rcut) ! 1E-3 is the tolerance

    double DensRcut = nr.ne[id]; !!!!!!!!!!!!!!!!!!!!!!!! TODO

*-----------------Locate L-cutoff-----------------

        wLcut = -0.5_R8*geo%cvBb(:,3)*qe/me + 0.5_R8*sqrt(4.0*dv%ne*
     &          qe**2/me/eps0 + (geo%cvBb(:,3)*qe/me)**2) ! frequency for UH resonance

        call maskfun(mpg%nCi, wLcut, res_freqV,
     &              1.0E-3_R8, maskf_Lcut) ! 1E-3 is the tolerance

    double DensLcut = nr.ne[id]; !!!!!!!!!!!!!!!!!!!!!!!! TODO

*-----------------Locate UHR-layer-----------------

        wuh = sqrt(dv%ne*qe**2/me/eps0 + (geo%cvBb(:,3)*qe/me)**2) ! frequency for UH resonance

        call maskfun(mpg%nCi, wuh, res_freqV,
     &              1.0E-3_R8, maskf_uh) ! 1E-3 is the tolerance

    Ruh = aR[id]; !!!!!!!!!!!!!!!!!!!!!!!! TODO
    double Densuh = nr.ne[id]; !!!!!!!!!!!!!!!!!!!!!!!! TODO

*-----------------Budden Analysis @ UHR-----------------

    double Lne = abs(nr.ne[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (nr.ne[id - 1] - nr.ne[id + 1])); // density scale length [m]
    double LBr = abs(Br[id] * abs(aR[id + 1] - aR[id - 1]) / 100.0 / (Br[id - 1] - Br[id + 1])); // magnetic scal length [m]
    double alphaUHR = sqrt(nr.ne[id] * 1.0e6 * qe * qe / me / eps0) / (Br[id] * qe / me); // wpe/wce @ UHR
    double Buddenparam = (Br[id]*qe/me)*Lne/c * alphaUHR / sqrt(alphaUHR * alphaUHR + 2.0 * (Lne / LBr)) * pow((sqrt(1.0 + alphaUHR * alphaUHR) - 1.0)/
			 (alphaUHR * alphaUHR + (Lne / LBr) * sqrt(1.0 + alphaUHR * alphaUHR)), 0.5); // Budden parameter
    double Tunnel = exp(-2.0 * Buddenparam); // Tunneling fraction
    double Reflect = pow((1.0 - Tunnel), 2.0); // Reflected part
    double Convert = 1.0 - Tunnel - Reflect; // XB converted fraction

*-----------------O-absorption-----------------

      DensO = res_freq**2*me*eps0/(qe**2)
      if (maxval(dv%ne).le.DensO) then
        PabsO = (1.0-exp(-0.25*pi*R/100.0/c/c/c*Teabsce * 11604.45*kb / me * wpe * wpe / wce * sqrt(1.0 - alphaCE))); // 1 pass aborption for O-wave
        PabsO = PabsO / (1.0 - (1.0 - muw) * (1.0 - PabsO)); // multi-pass absorption with losses to the wall (muw)
      endif

*-----------------X-absorption-----------------

    double PabsX = (1.0 - exp(-0.5 * pi * R / 100.0 / c / c / c * Teabsce * 11604.45 * kb / me * wpe * wpe / wce * pow(2.0 - alphaCE, 2.5) * (1.0 +
		   sqrt(alphaCE)))); // 1 pass aborption for X-wave below Rcutoff
    DensRcut = 3.13e9;
    Densuh = 6.26e9;
    if (nemax < DensRcut) {
        PabsX = PabsX * abs(DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE - 1.0e-5) / (DensRcut * 1.0e6 * me / eps0 / Bce / Bce - alphaCE +
		1.0e-5); // adaptation for X-wave nearing Rcutoff
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX)); // multi-pass absorption with losses to the wall (muw)
    }
    else if (nemax < Densuh) {
        PabsX = (1.0-Reflect)*PabsX/(1.0 - (1.0 - muw) * (1.0 - PabsX)); // multi-pass absorption of the tunneling wave with losses to the wall (muw)
        PabsX = PabsX / (1.0 - (1.0 - muw) * (1.0 - PabsX)); // multi-pass absorption with losses to the wall (muw)
    }
    // Not reflected part tunnels and reaches ECR,, reflecting on wall and coupling at second pass
    else if (nemax < DensLcut) {
        PabsX = (1.0 - PRdep2) * Tunnel * PabsX + (1.0 - PRdep4) * Tunnel * (1.0 - PabsX) * (1.0 - muw) * PabsX;
    }
    // X-wave tunneling through R-cutoff, reflecting on wall and coupling at second pass
    // third pass not possible because of evanescent region, SX-B conversion possible
    else {
        PabsX = 0.0;
    } // All X-wave is converted to B-wave

*-----------------B-absorption-----------------

    double PabsB = 0.0; // EBW coupling near ECR
    id = 0;
    while (Br[id] > Bce)
        ++id;
    if ((abs(Br[id] - Bce)) > (abs(Br[id + 1] - Bce))) {
        ++id;
    }
    //no shift to LFS of ECR implemented for the moment, doppler shift to be added ????
    double RdepB = Rdep3; // Hardcoded position of B-coupling

    if (nemax < Densuh) {
        PabsB = 0.0;
    }
    // No B-wave present

    else if (nemax < DensO) {
        PabsB = PRdep2 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE)); // percentage of converted X-wave to B-wave coupling
    }
    // B-wave appears when UHR appears
    else {
        PabsB = PRdep3 * Convert * (1.0 - exp(-0.5 * pi * R / 100.0 / c * wce * alphaCE));
        PabsB += PRdep3 * (1.0 - PRdep4) * (1.0 - PRdep4) * Tunnel;
    }
    // double conversion to B-wave present: 1st converted at UHR, 2nd from tunneling X-wave reflected at Lcutoff and converted to B-wave at UHR

*-----------------UH-resonance-----------------

    double PabsU = 0.0; // coupling at the UHR

    if (nemax < Densuh) {
        PabsU = 0.0;
    }
    // No UHR present
    else if (nemax < DensO) {
        PabsU = PRdep4 * (Tunnel); // PRdep4 = percentage of UHR coupling of tunneled X-wave
    } else {
        PabsU = PRdep4 * (Tunnel); // 1st resonance at UHR
        PabsU += PRdep4 * (1.0 - PRdep4) * (Tunnel); // 2nd resonance at UHR from tunneled X-wave reflected at Lcutoff
    }

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof1[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Rdepce) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0)); // normaalverdeling met hoogte 1 voor O/X-wave deposition at ECR
        eprof2[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - RdepB) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));  // B-wave deposition at ECR+1.5cm (should be moving from 0 to + 2.0 ???)
        eprof3[id] = 0.01 * aR[id] * pow(1.0 + pow((aR[id] - Ruh) / (pow(2.0 * a / 100.0, 0.5) * 0.5), 2.0), -(0.5 + 1.0));    // UH Resonance
    }

    for (int id = 0; id < NMESHP - 1; ++id) { // line integrated density
        inteprof1 += 2.0 * b * pi * 0.5 * (eprof1[id] + eprof1[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof2 += 2.0 * b * pi * 0.5 * (eprof2[id] + eprof2[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
        inteprof3 += 2.0 * b * pi * 0.5 * (eprof3[id] + eprof3[id + 1]) * (pow(aR[id + 1], 2.0) - pow(aR[id], 2.0));
    }

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) { // line integrated density
        eprof1[id] = eprof1[id] / inteprof1;
        eprof2[id] = eprof2[id] / inteprof2;
        eprof3[id] = eprof3[id] / inteprof3;
    }

    // #pragma omp parallel for
    for (int id = 0; id < NMESHP; ++id) {
        Pabs1[id] = (PRdep1 * PabsO * eprof1[id]) * (Prf * 1.0e3); // Power (kW) of O-wave absorption at ECR
        Pabs2[id] = (PRdep2 * PabsX * eprof1[id]) * (Prf * 1.0e3); // Power (kW) of X-wave absorption at ECR (inluding tunneled)
        Pabs3[id] = (PRdep2 * PabsB * eprof2[id]) * (Prf * 1.0e3); // Power (kW) of B-wave absorption near ECR (converted from X-wave => PRdep2)
        Pabs4[id] = (PRdep2 * PabsU * eprof3[id]) * (Prf * 1.0e3); // Power (kW) of UHR absorption at UHR
        PRFe_array[id] = 6.24e18 * Pabs1[id] + 6.24e18 * Pabs2[id] + 6.24e18 * Pabs3[id] + 6.24e18 * Pabs4[id];
        PRFHi_array[id] = 0.0;
        PRFH2i_array[id] = 0.0;
        PRFH3i_array[id] = 0.0;
        PRFHeII_array[id] = 0.0;
        PRFHeIII_array[id] = 0.0;
    }
}
#endif
* -- COUPLING = 3 -- Ray-tracing with TORAY. !!TODO

      case(3)
        call xerrab ('b2mod_ech -- ECH_COUPLING = 3 not coded yet')

      end select

      sECH = she0_ech
      if (switch%ech_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nCv,0,chi,'b2mod_ECH_shapefun')
        call my_outi_us(70,nCv,0,maskf,'b2mod_ECH_maskf')
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,sECH(:,cntr),'b2mod_ECH_she_'//chns)
        enddo
      endif

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end CalcECH

      end subroutine CalcECH

c Find resonance layer
      subroutine find_resonance(nCv,geo,mpg,switch,omegace,
     &           res_freq,maskf)
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
      real (kind=R8):: omegace(nCv), res_freq, B_res, B_resV(nCv)
c     Output variables
      integer :: maskf(nCv)

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('find_resonance')
      call sfill (nCv, 0.0_R8, omegace, 1)
      call ifill (nCv, 0, maskf, 1)

      call elcyclfr(nCv,geo,omegace)

      res_freq = 2*pi*switch%ech_RFfreq/switch%ech_harmonics
      B_res = me*res_freq/qe
      call sfill (nCv, B_res, B_resV, 1)

      call maskfun(mpg%nCi, B_resV, geo%cvBb(:,3),
     &             1.0E-3_R8, maskf) ! 1E-3 is the tolerance

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end elcyclfr

      end subroutine find_resonance

c Electron Cyclotron Frequency computation
      subroutine elcyclfr(nCv,geo,omegace)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      type (geometry), intent (in) :: geo
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: omegace(nCv)

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('elcyclfr')
      call sfill (nCv, 0.0_R8, omegace, 1)

      omegace(:) = qe*geo%cvBb(:,3)/me

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end elcyclfr

      end subroutine elcyclfr

c Dealloc
      subroutine dealloc_b2mod_ech
      implicit none
      if (allocated(sECH)) then
        deallocate(sECH)
      endif

      return
      end subroutine dealloc_b2mod_ech
c End dealloc

*-----------------------------------------------------------------------
*.end b2mod_ech

      end module b2mod_ech
