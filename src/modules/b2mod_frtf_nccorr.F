      module b2mod_frtf_nccorr
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2zhco
     & , only : c_R_TA_an_save, c_hTA_an_fl_save
      use b2mod_b2mrel
      use b2mod_subsys
      use b2mod_zhfrtf
      use b2us_geo
      implicit none
      real (kind=R8), allocatable ::
     &  alpha_HS_style(:,:), g_HS_style(:),
     &  corr_tfia(:,:), corr_fria(:,:)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_FRTF_NCCORR stores a parameter of relative impurity presence,
*     main ion collisionality (neoclassical style, see the manual) and
*     correction coefficients,
*     that can be used to make impurity transport inside the separatrix
*     closer to neoclassical.
*     It contains subroutines to calculate correction coefficients and
*     the friction force (Zhdanov closure style, see b2mod_zhfrtf_fr) with
*     neoclassical correction applied.
*
*  3. description
*
*     This module contains:
*     - alpha_HS_style: ratio of the squared charge weighted density of
*                       impurity ions to the main ion density
*     - g_HS_style    : collisionality of main ions multiplied by \eps^{3/2}
*                       and divided by a component containing impurity presence
*     - corr_fria     : a correction coefficient that is designed
*                       only for friction of chosen medium-to-high
*                       charged impurity with main ions
*     - corr_tfia     : a correction coefficient that is designed
*                       only for thermal force of chosen medium-to-high
*                       charged impurity with main ions
*-----------------------------------------------------------------------
      contains
!.......................................................................

      subroutine alloc_b2mod_frtf_nccorr(nCv, ns)
      implicit none
      integer, intent(in) :: nCv, ns

      allocate(alpha_HS_style(1:nCv,0:ns-1))
      allocate(g_HS_style(1:nCv))
      allocate(corr_tfia(1:nCv,0:ns-1))
      allocate(corr_fria(1:nCv,0:ns-1))

      return
      end subroutine alloc_b2mod_frtf_nccorr

      subroutine dealloc_b2mod_frtf_nccorr
      implicit none

      if (.not.allocated(alpha_HS_style)) return
      deallocate(alpha_HS_style)
      deallocate(g_HS_style)
      deallocate(corr_tfia)
      deallocate(corr_fria)

      return
      end subroutine dealloc_b2mod_frtf_nccorr

!.......................................................................

      subroutine b2mod_frtf_nccorr_prep(nCv,nFt,mpg,geo,pl,co,rt,zeff)
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none

      integer :: nCv, nFt
      type (mapping) :: mpg
      type (geometry) :: geo
      type (B2Plasma) :: pl
      type (B2Coeff) :: co
      type (B2Rates) :: rt
      real (kind=R8) :: zeff(1:nCv)
*   ..local variables
      integer :: iCv, is

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_FRTF_NCCORR_PREP stores a parameter of relative impurity presence
*     and main ion collisionality (neoclassical style, see the manual).
*
*  3. description
*
*     This module contains:
*     - alpha_HS_style: ratio of the squared charge weighted density of
*                       impurity ions to the main ion density
*     - g_HS_style    : collisionality of main ions multiplied by \eps^{3/2}
*                       and divided by a component containing impurity presence
*
*-----------------------------------------------------------------------

      do iCv = 1, mpg%nCi
        do is = 0, pl%ns-1
          alpha_HS_style(iCv,is)=pl%na(iCv,is)*
     &     rt%rz2(iCv,is)/pl%na(iCv,pl%ismain)
        enddo
        if (mpg%cvFt(iCv).gt.0) g_HS_style(iCv)=co%nu2(iCv)*
     &   sqrt(geo%ftEps(mpg%cvFt(iCv))**3)/zeff(iCv)
       enddo

      return
      end subroutine b2mod_frtf_nccorr_prep

      subroutine b2mod_frtf_nccorr_coefs(nCv,nFt,mpg,geo,
     &  pl,co,rt,switch,chTAancv,cRTAancv,z2n_cv,zeff)
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_switches
      implicit none
      integer, intent(in) :: nCv, nFt
      type (mapping) :: mpg
      type (geometry) :: geo
      type (B2Plasma) :: pl
      type (B2Coeff) :: co
      type (B2Rates) :: rt
      type (switches) :: switch
      real (kind=R8) ::
     &  zeff(1:nCv),
     &  chTAancv(1:nCv,0:nspecies-1),
     &  cRTAancv(1:nCv,0:nspecies-1,0:nspecies-1),
     &  z2n_cv(1:nCv,0:nspecies-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_FRTF_NCCORR_COEFS calculates correction coefficients.
*
*  3. description
*
*     This routine computes:
*     - corr_fria: a correction coefficient that is designed
*                  only for friction of chosen medium-to-high
*                  charged impurity with main ions
*     - corr_tfia: a correction coefficient that is designed
*                  only for thermal force of chosen medium-to-high
*                  charged impurity with main ions
*
*-----------------------------------------------------------------------
*  ..local variabless
      character chns*3
      integer, save :: iout = 0
      integer :: corr_spec
      integer :: iCv, iFt, is
      real (kind=R8) :: temp_corr(1:nCv,0:2,0:pl%ns-1),
     &  corr_fria_D(1:nCv),corr_tfia_D(1:nCv),
     &  k2_mean(1:nCv),k2_tmp(1:nFt,0:1)
*-----------------------------------------------------------------------

      corr_spec = switch%Zhdanov_nc_spec - 1
      call ipgeti ('b2sifrtf_iout', iout)
      call b2mod_frtf_nccorr_prep(nCv,nFt,mpg,geo,pl,co,rt,zeff)
      corr_fria_D=0.0_R8
      corr_tfia_D=0.0_R8
      corr_fria=1.0_R8
      corr_tfia=1.0_R8
      temp_corr=0.0_R8
      k2_mean=0.0_R8
      k2_tmp=0.0_R8
      do iCv = 1, mpg%nCi
        if (mpg%cvOnClosedSurface(iCv)) then
          iFt = mpg%cvFt(iCv)
          if (iFt.gt.0) then  
            k2_tmp(iFt,0)=k2_tmp(iFt,0)+co%k2(iCv)*geo%cvVol(iCv)
            k2_tmp(iFt,1)=k2_tmp(iFt,1)+geo%cvVol(iCv)
          endif
        endif
      enddo
      do iCv = 1, mpg%nCi
        if (mpg%cvOnClosedSurface(iCv)) then
          iFt = mpg%cvFt(iCv)
          if (iFt.gt.0) then  
            k2_mean(iCv)=k2_tmp(iFt,0)/k2_tmp(iFt,1)
          endif 
        endif
      enddo

      do iCv = 1, mpg%nCi
* coefficients are calculated only inside the separatrix and for plasma with low main ion collisionality
        if (mpg%cvOnClosedSurface(iCv) .and.
     &     (g_HS_style(iCv).le.5.0_R8)) then !this is a temporary decision
          iFt = mpg%cvFt(iCv)
          do is = 0, pl%ns-1
c                if (switch%Zhdanov_closure.eq.1) then

             if (rt%rza(iCv,is).ge.switch%Zhdanov_nc_min.and.iFt.gt.0)
     >        then
* coefficients are calculated only for chosen ions of the chosen impurity
              corr_fria(iCv,is)=
     &           1.0_R8+(-0.36_R8*alpha_HS_style(iCv,is)+
     &                    0.44_R8*alpha_HS_style(iCv,is)**(0.77_R8)+
     &                    0.01_R8)/(1.0_R8+g_HS_style(iCv)**2)

              temp_corr(iCv,0,is)=
     &         (chTAancv(iCv,s2nucl(pl%ismain)-1)*sqrt(2.0_R8))/
     &         (cRTAancv(iCv,s2nucl(pl%ismain)-1,s2nucl(is)-1)*2.5_R8)
              temp_corr(iCv,1,is)=
     &         (2.4_R8*sqrt(geo%ftEps(iFt)**3))/(1.0_R8+(0.36_R8+0.4_R8*
     &          alpha_HS_style(iCv,is))*g_HS_style(iCv)**2)
              temp_corr(iCv,2,is)=
     &         (1.6_R8*(geo%ftEps(iFt)**3)*g_HS_style(iCv)**2)/
     &         (1.0_R8+sqrt(geo%ftEps(iFt)**3)*g_HS_style(iCv)**2)
              corr_tfia(iCv,is)=
     &         (temp_corr(iCv,0,is)*temp_corr(iCv,1,is)+
     &          temp_corr(iCv,2,is))/k2_mean(iCv)
            endif

* the correction for main ions is calculated to make total forces zero
            if ((s2nucl(is)-1).eq.corr_spec) then

              corr_fria_D(iCv)=corr_fria_D(iCv)+
     &          corr_fria(iCv,is)*rt%rz2(iCv,is)*pl%na(iCv,is)/
     &          z2n_cv(iCv,s2nucl(is)-1)
              corr_tfia_D(iCv)=corr_tfia_D(iCv)+
     &          corr_tfia(iCv,is)*rt%rz2(iCv,is)*pl%na(iCv,is)/
     &          z2n_cv(iCv,s2nucl(is)-1)
            endif
          enddo
          corr_fria(iCv,pl%ismain)=corr_fria_D(iCv)
          corr_tfia(iCv,pl%ismain)=corr_tfia_D(iCv)

        endif
      enddo

      if (iout.ne.0) then !  .or. iout_b2wdat.eq.4) then
        do is = 0, pl%ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,temp_corr(:,0,is),
     &     'b2mod_frtf_nccorr_coefs_temp0_'//chns)
          call my_out_us(70,nCv,0,temp_corr(:,1,is),
     &     'b2mod_frtf_nccorr_coefs_temp1_'//chns)
          call my_out_us(70,nCv,0,temp_corr(:,2,is),
     &     'b2mod_frtf_nccorr_coefs_temp2_'//chns)
        enddo
        call my_out_us(70,nCv,0,k2_mean,
     &     'b2mod_frtf_nccorr_coefs_k2mean_'//chns)
      endif

      return
      end subroutine b2mod_frtf_nccorr_coefs

      end module b2mod_frtf_nccorr
