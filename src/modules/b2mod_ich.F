*-----------------------------------------------------------------------
*.specification
c**************************************************************************************
c**********************  NEW KU Leuven module for ECH plasmas  ************************
c**************************************************************************************
! mb  16/10/23
! main author: Elise Chambard
      module b2mod_ich
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2us_plasma

      implicit none
*   ..declarations allocatable variable
      real (kind=R8), allocatable, save ::
     &  sICH(:,:)

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ICH is intended to simulate power depositions in Ion
*     Cyclotron Heated plasmas. To turn on the model, the flag
*     'b2sral_use_ICH' has to be set to '1'. This module computes an
*     additional heating source due to IC heating injected either at the
*     cyclotron frequency or its harmonics.
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      contains

c Alloc
      subroutine alloc_b2mod_ich(nCv)
      implicit none
      integer, intent(in) :: nCv

      allocate(sICH(nCv,0:3))

      sICH(:,:) = 0.0_R8

      return
      end subroutine alloc_b2mod_ich
c End alloc

c ICH power computation
      subroutine CalcICH(nCv,ns,geo,mpg,pl,dv,switch,shi0_ich,she0_ich)
c     This routine calculates the ICH power deposition
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_math
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv,ns
c     Output variables
      real (kind=R8), intent(out):: shi0_ich(nCv,0:3), she0_ich(nCv,0:3)
c     Other variables
      real (kind=R8) :: omegaci(nCv), Rdep_vect(nCv)
      integer :: maskf(nCv)
c     Internal variables
      real (kind = R8) :: chi(nCv), sx, sigma_y, intchi, normlz, h,
     &                    summation, res_freq
      real (kind = R8) :: segment(2,2), segmdist(nCv), totaldist
      real (kind = R8) :: fcCnX1, fcCnX2, fcCnY1, fcCnY2
      integer :: list_cv(nCv), nncv, list_ci(mpg%nCi), nnci,
     &           list_fc(mpg%nFc), nnfc, face1, face2
      integer :: iCv, iCv2, iFc, cntr, is
      character chns*3
      logical:: ifail
c     Variables for coupling 4
      real (kind = R8) :: P_rest, P_norm, chi_tot(nCv), NF,
     &                    eps_perp(mpg%nCi), sign_eps_perp, R_lhr(nCv),
     &                    lhr_background, test_lhr(mpg%nCi)
      real (kind = R8) :: omega2_in, omega2_cycl_e(nCv), 
     &                    omega2_cycl_Hi(nCv),omega2_plas_e(nCv),
     &                    omega2_plas_Hi(nCv), Rp_min,
     &                    omega2_cycl_Hei(nCv), omega2_plas_Hei(nCv),
     &                    omega2_cycl_Heii(nCv), omega2_plas_Heii(nCv),
     &                    real_ich_Rdep, XI(mpg%nFc), YI(mpg%nFc)
      integer :: iAnt, il_HFS, id_lhr, N_lhr
      integer :: list_cv_Rp_min(nCv), nncv_Rp_min, 
     &           list_ci_Rp_min(mpg%nCi), nnci_Rp_min,
     &           list_fc_Rp_min(mpg%nFc), nnfc_Rp_min, i_R_lhr(nCv)
      integer :: list_cv_R_ant(nCv), nncv_R_ant, list_ci_R_ant(mpg%nCi),
     &           nnci_R_ant,
     &           list_fc_R_ant(mpg%nFc), nnfc_R_ant
      integer :: list_cv_Z(nCv), nncv_Z, list_ci_Z(mpg%nCi), nnci_Z,
     &           list_fc_Z(mpg%nFc), nnfc_Z, i_R_Rp
      integer :: id, id2, i_R_Ant, i_Z_Ant, i_Cell_Ant,
     &           i_Z_Rp, i_Cell_Rp


* ..preliminaries
*   ..subprogram start-up calls
      call subini ('CalcICH')
      call sfill (nCv*4, 0.0_R8, shi0_ich, 1)
      call sfill (nCv*4, 0.0_R8, she0_ich, 1)
      call ifill (nCv, 0, maskf, 1)
      call sfill (nCv, switch%ich_Rdep, Rdep_vect, 1)
      call sfill (nCv, 0.0_R8, chi, 1)
      call sfill (4, 0.0_R8, segment, 1)
      call sfill (nCv, 0.0_R8, segmdist, 1)
c      inteprof = 0.0_R8
c      i_R_Rp = -120 !why this value?
      summation = 0.0_R8
      totaldist = 0.0_R8
      fcCnX1 = 0.0_R8
      fcCnX2 = 0.0_R8
      fcCnY1 = 0.0_R8
      fcCnY2 = 0.0_R8
      cntr = 0
c      normlz = 0.0_R8

      select case (switch%ich_coupling)


* -- COUPLING = 0 -- Fixed power distributed as a Dirac delta in volumes intersected at given radial position.

      case(0)
        segment = reshape((/ switch%ich_Rdep, -0.1_R8,
     &       switch%ich_Rdep, maxval(geo%cvY)+0.1_R8 /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv, list_ci,
     &   list_fc, nncv, nnci, nnfc, nCv, XI, YI, 1, ifail)

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
        summation = sum(chi(:)*geo%cvVol(:))
        chi = chi/summation
        do cntr = 1, size(list_ci)
          if(list_ci(cntr).ne.0) then
            iCv = list_ci(cntr)
            shi0_ich(iCv,0) = chi(iCv)*switch%ich_absfract*
     &                        switch%ich_power*geo%cvVol(iCv)
          endif
        enddo


* -- COUPLING = 1 -- Fixed fraction of launched power is coupled/Gaussian injection.

      case(1)
        sx = switch%ich_gaussian_width_x ! widthich in T1D
        sigma_y = switch%ich_gaussian_width_y ! used to be defined as "sy" but this is also a function in geo so doesn't work
        if (switch%ich_gaussian_symmetry.eq.1) then
          do iCv = 1, mpg%nCi
            chi(iCv) = geo%cvX(iCv)*
     &            expu(-1.0_R8*
     &            ((geo%cvX(iCv)-switch%ich_Rdep)/sx)**2)+
     &            1.0E-7*geo%cvX(iCv)*expu(-1.0_R8*
     &            ((geo%cvX(iCv)-switch%ich_Rdep)/switch%a0)**2) ! symmetry on poloidal direction
          enddo

c     The mask function chi was changed on 28/04 again to have the exact same profile as in TOMATOR: no -1/2 but -1 in exp

c          segment = reshape((/
c     &       minval(geo%cVx), maxval(geo%cvY)/2.0_R8,
c     &       maxval(geo%cVx), maxval(geo%cvY)/2.0_R8 /),
c     &              shape(segment))

c          call cv_intersections (geo, mpg, segment, list_cv, list_ci,
c     &   list_fc, nncv, nnci, nnfc, nCv, 1, ifail)

          summation = sum(chi(:)*geo%cvVol(:))

c          do cntr = 1, size(list_ci)
c            if(list_ci(cntr).ne.0) then
c              iCv = list_ci(cntr)
c              summation = summation + chi(:)*geo%cvVol(:)
c            endif
c          enddo

          chi = chi/summation

        else
          chi = 1.0_R8/(2*pi*sx*sigma_y)*exp(-1.0_R8/2.0_R8*
     &            ((geo%cvX-switch%ich_Rdep)**2/sx**2 +
     &             (geo%cvY-switch%ich_Zdep)**2/sigma_y**2)) ! ich_Zdep to be defined

          do iCv = 1, mpg%nCi
            summation = summation + chi(iCv)*geo%cvVol(iCv) ! normalisation condition to ensure grid independence
          enddo
          chi = chi/summation ! normalisation condition to ensure grid independence

        endif

        do iCv = 1, mpg%nCi
          shi0_ich(iCv,0) = switch%ich_power*chi(iCv)*
     &                      (1-switch%ich_fracpne)*
     &                      switch%ich_absfract*geo%cvVol(iCv)
          she0_ich(iCv,0) = switch%ich_power*chi(iCv)*
     &                      switch%ich_fracpne*
     &                      switch%ich_absfract*geo%cvVol(iCv)
        enddo

* -- COUPLING = 2 -- lhr coupling !!TODO

      case(2)

      P_rest = 1.0_R8
      P_norm = 0.0_R8
      N_lhr = 0
      NF = 0.0_R8
      iAnt = 0
      sx = switch%ich_gaussian_width_x

      omega2_in = (2*pi*switch%ich_RFfreq)**2
      Rp_min = switch%R0 - switch%limiter_width * switch%a0 ! this is lHFS in T1D

      segment = reshape((/ Rp_min, -0.1_R8,
     &       Rp_min, maxval(geo%cvY)+0.1_R8 /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv_Rp_min,
     &   list_ci_Rp_min, list_fc_Rp_min, nncv_Rp_min, nnci_Rp_min,
     &   nnfc_Rp_min, nCv, XI, YI, 1, ifail)

      segment = reshape((/ switch%ich_Rdep, -0.1_R8,
     &       switch%ich_Rdep, maxval(geo%cvY)+0.1_R8 /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv_R_ant,
     &   list_ci_R_ant, list_fc_R_ant, nncv_R_ant, nnci_R_ant,
     &   nnfc_R_ant, nCv, XI, YI, 1, ifail)


      segment = reshape((/ minval(geo%cVx) -0.1_R8, 
     &          maxval(geo%cvY)/2.0_R8, maxval(geo%cVx)+0.1_R8,
     &          maxval(geo%cvY)/2.0_R8 /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv_Z, list_ci_Z,
     &   list_fc_Z, nncv_Z, nnci_Z, nnfc_Z, nCv, XI, YI, 1, ifail)


      call sfill (nCv, 0.0_R8, chi, 1)
      call sfill (nCv, 0.0_R8, chi_tot, 1)
      call sfill (nCv, 0.0_R8, R_lhr, 1)

      do iCv=1, nCv
        i_R_lhr(iCv) = 0
      enddo

c      do while (geo%cvX(iAnt) < switch%ich_Rdep) ! can be left out because we have i_Ant as a vector of indexes from cv_intersections
c          iAnt = iAnt+1
c      enddo
c      do while (geo%cvX(iAnt) > Rp_min) ! can be left out because we have il_HFS as a vector of indexes from cv_intersections
c          il_HFS = il_HFS-1
c      enddo

      omega2_cycl_e = (geo%cvBb(:,3)*qe/me)**2
      omega2_plas_e = (dv%ne(:)*qe**2/me/eps0)

      omega2_cycl_Hi = (geo%cvBb(:,3)*qe/(am(1)*mp))**2
      omega2_plas_Hi = (pl%na(:,1)*qe**2/(am(1)*mp)/eps0)
      omega2_cycl_Hei = (geo%cvBb(:,3)*qe/(am(3)*mp))**2
      omega2_plas_Hei = (pl%na(:,3)*qe**2/(am(3)*mp)/eps0)
      omega2_cycl_Heii = (geo%cvBb(:,3)*(2*qe)/(am(4)*mp))**2
      omega2_plas_Heii = (pl%na(:,4)*(2*qe)**2/(am(4)*mp)/eps0)


      do iCv = 1,mpg%nCi
        eps_perp(iCv) = 1.0_R8 -omega2_plas_e(iCv)/
     &                 (omega2_in-omega2_cycl_e(iCv))
     &                 -omega2_plas_Hi(iCv)/
     &                 (omega2_in-omega2_cycl_Hi(iCv))
     &                 -omega2_plas_Hei(iCv)/
     &                 (omega2_in-omega2_cycl_Hei(iCv))
     &                 -omega2_plas_Heii(iCv)/
     &                 (omega2_in-omega2_cycl_Heii(iCv))
      enddo


      do id = 1,nnci_Z
          do id2 = 1,nnci_R_ant
            if (list_ci_Z(id) == list_ci_R_ant(id2)) then
              i_R_Ant = id ! index in the radial direction
              i_Z_Ant = id2
              i_Cell_Ant = list_ci_Z(id)
            endif
          enddo
      enddo


      do id = 1,nnci_Z
          do id2 = 1,nnci_Rp_min
            if (list_ci_Z(id) == list_ci_Rp_min(id2)) then
              i_R_Rp = id
              i_Z_Rp = id2 ! normally this should be equal to i_Z_Ant
              i_Cell_Rp = list_ci_Z(id)
            endif
          enddo
      enddo

      sign_eps_perp = sign(1.0_R8,eps_perp(i_Cell_Ant))

      real_ich_Rdep = geo%cvX(i_Cell_Ant)

      do id = i_R_Ant-1,i_R_Rp,-1
          iCv = list_ci_Z(id)
          if (sign_eps_perp /= sign(1.0_R8,eps_perp(iCv))) then
            N_lhr = N_lhr+1
            R_lhr(iCv)=P_rest*switch%ich_absfract
            i_R_lhr(iCv) = 1
            P_rest = P_rest*(1-switch%ich_absfract)
          endif
          sign_eps_perp = sign(1.0_R8,eps_perp(iCv))
      enddo

      do id = i_R_Rp,i_R_Ant
          iCv = list_ci_Z(id)
          if (i_R_lhr(iCv) /= 0) then
            R_lhr(iCv) = R_lhr(iCv)+P_rest*switch%ich_absfract
            P_rest = P_rest*(1-switch%ich_absfract)
            P_norm = P_norm+R_lhr(iCv)
          endif
      enddo

      do id = i_R_Rp,i_R_Ant
          iCv = list_ci_Z(id)
          if (i_R_lhr(iCv) /= 0) then
            R_lhr(iCv) = R_lhr(iCv)/P_norm
          endif
      enddo

      call sfill (mpg%nCi, 0.0_R8, test_lhr, 1)

      if (N_lhr == 0) then
          test_lhr(1) = 1.0_R8
          N_lhr = 1
          i_R_lhr(i_Cell_Ant) = 1
          R_lhr(i_Cell_Ant) = 1.0_R8
      endif

      lhr_background = 1.0E-8_R8
      do id = i_R_Ant,i_R_Rp,-1
          id_lhr = list_ci_Z(id)
          if (i_R_lhr(id_lhr) /= 0) then
            if (switch%ich_gaussian_symmetry.eq.1) then
              do iCv = 1,mpg%nCi
                chi(iCv) = geo%cvX(iCv)*expu(-1.0_R8*((geo%cvX(iCv)-
     &                     geo%cvX(id_lhr))/sx)**2)+
     &                     lhr_background*geo%cvX(iCv)*expu(-1.0_R8*
     &                     ((geo%cvX(iCv)-geo%cvX(id_lhr))/
     &                     switch%a0)**2)
              enddo
              NF = 0.0_R8
              NF = sum(chi(:)*geo%cvVol(:))
              do iCv=1,mpg%nCi
                chi(iCv) = chi(iCv)/NF*R_lhr(id_lhr)
                chi_tot(iCv) = chi_tot(iCv)+chi(iCv)
              enddo
            endif
          endif
      enddo

      do iCv = 1, mpg%nCi
          shi0_ich(iCv,0) = switch%ich_power*chi_tot(iCv)*
     &                      (1-switch%ich_fracpne)*
     &                      geo%cvVol(iCv)
          she0_ich(iCv,0) = switch%ich_power*chi_tot(iCv)*
     &                      switch%ich_fracpne*
     &                      geo%cvVol(iCv)
      enddo


* -- COUPLING = 3 -- 2D Gauss coupling !!TODO

      case(3)

        sx = switch%ich_gaussian_width_x
        sigma_y = switch%ich_gaussian_width_y

        do iCv = 1,mpg%nCi
          chi(iCv) = 1.0_R8/(2.0_R8*pi*sx*sigma_y)*expu(-1.0_R8/2.0_R8*
     &          (((geo%cvX(iCv)-switch%ich_Rdep)/sx)**2+
     &          ((geo%cvY(iCv)-switch%ich_Zdep)/sigma_y)**2))
        enddo

        NF = sum(chi(:)*geo%cvVol(:))

        chi = chi/NF

        do iCv = 1, mpg%nCi
          shi0_ich(iCv,0) = switch%ich_power*chi(iCv)*
     &                      (1-switch%ich_fracpne)*
     &                      switch%ich_absfract*geo%cvVol(iCv)
          she0_ich(iCv,0) = switch%ich_power*chi(iCv)*
     &                      switch%ich_fracpne*
     &                      switch%ich_absfract*geo%cvVol(iCv)
        enddo

      end select

      sICH = shi0_ich
      if (switch%ich_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nCv,0,chi,'b2mod_ICH_shapefun')
        call my_out_us(70,nCv,0,eps_perp(:),'b2mod_ICH_eps_perp')
        call my_out_us(70,nCv,0,test_lhr(:),'b2mod_ICH_test_lhr')
        call my_outi_us(70,nCv,0,maskf,'b2mod_ECH_maskf')
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,sICH(:,cntr),'b2mod_ICH_shi_'//chns)
          call my_out_us(70,nCv,0,she0_ich(:,cntr),
     &         'b2mod_ICH_she_'//chns)
        enddo
      endif

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end CalcICH

      end subroutine CalcICH

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

      res_freq = 2*pi*switch%ich_RFfreq/switch%ich_harmonics
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
      subroutine dealloc_b2mod_ich
      implicit none
      if (allocated(sICH)) then
        deallocate(sICH)
      endif

      return
      end subroutine dealloc_b2mod_ich
c End dealloc

*-----------------------------------------------------------------------
*.end b2mod_ich

      end 

