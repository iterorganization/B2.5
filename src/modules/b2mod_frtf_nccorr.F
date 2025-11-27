      module b2mod_frtf_nccorr
      use b2mod_types
      use b2mod_geo
     & , only : vol
      use b2mod_b2cmpa
      use b2mod_b2zhco
     & , only : c_R_TA_an_save, c_hTA_an_fl_save
      use b2mod_b2mrel
      use b2mod_subsys
      use b2mod_indirect
     , , only : on_closed_surface
      use b2mod_transport
     , , only : nu2, eps1, k2
      use b2mod_zhfrtf
      implicit none
      real (kind=R8), allocatable ::
     &  alpha_HS_style(:,:,:), g_HS_style(:,:),
     &  corr_tfia(:,:,:), corr_fria(:,:,:)
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

      subroutine alloc_b2mod_frtf_nccorr(nx, ny, ns)
      implicit none
      integer, intent(in) :: nx, ny, ns

      allocate(alpha_HS_style(-1:nx,-1:ny,0:ns-1))
      allocate(g_HS_style(-1:nx,-1:ny))
      allocate(corr_tfia(-1:nx,-1:ny,0:ns-1))
      allocate(corr_fria(-1:nx,-1:ny,0:ns-1))

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

      subroutine b2mod_frtf_nccorr_prep(nx,ny,ns,na,ismain,rz2,
     &    nu2,eps1,zeff)
      use b2mod_indirect
      implicit none

      integer :: nx, ny, ns, ismain
      real (kind=R8) :: na(-1:nx,-1:ny,0:ns-1),
     &  rz2(-1:nx,-1:ny,0:ns-1),
     &  nu2(-1:nx,-1:ny), eps1(-1:ny),zeff(-1:nx,-1:ny)
*   ..local variables
      integer :: ix, iy, is

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

      do iy = 0, ny-1
        do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            do is = 0, ns-1
              alpha_HS_style(ix,iy,is)=na(ix,iy,is)*
     &         rz2(ix,iy,is)/na(ix,iy,ismain)
            enddo
            g_HS_style(ix,iy)=nu2(ix,iy)*
     &       sqrt(eps1(iy)**3)/zeff(ix,iy)
          endif
        enddo
      enddo

      return
      end subroutine b2mod_frtf_nccorr_prep

      subroutine b2mod_frtf_nccorr_coefs(na,nx,ny,ns,ismain,rza,
     &   rz2,corr_Z_min,corr_spec,chTAanxy,cRTAanxy,z2n_xy,zeff)
      use b2mod_indirect
      implicit none
      integer :: nx, ny, ns, ismain, corr_Z_min,
     &    corr_spec
      real (kind=R8) ::
     &  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1),
     &  zeff(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),
     &  chTAanxy(-1:nx,-1:ny,0:nspecies-1),
     &  cRTAanxy(-1:nx,-1:ny,0:nspecies-1,0:nspecies-1),
     &  z2n_xy(-1:nx,-1:ny,0:nspecies-1)
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
      integer :: ix, iy, is
      real (kind=R8) :: temp_corr(-1:nx,-1:ny,0:2,0:ns-1),
     &  corr_fria_D(-1:nx,-1:ny),corr_tfia_D(-1:nx,-1:ny),
     &  k2_mean(-1:nx,-1:ny),k2_tmp(-1:ny,0:1)
*-----------------------------------------------------------------------

      call ipgeti ('b2sifrtf_iout', iout)
      call b2mod_frtf_nccorr_prep(nx,ny,ns,na,ismain,rz2,
     &  nu2,eps1,zeff)
      corr_fria_D(:,:)=0.0e0_R8
      corr_tfia_D(:,:)=0.0e0_R8
      corr_fria(:,:,:)=1.0e0_R8
      corr_tfia(:,:,:)=1.0e0_R8
      temp_corr(:,:,:,:)=0.0e0_R8
      k2_mean(:,:)=0.0e0_R8
      k2_tmp(:,:)=0.0e0_R8
      do iy=-1,ny
        do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            if (on_closed_surface(ix,iy)) then
              k2_tmp(iy,0)=k2_tmp(iy,0)+k2(ix,iy)*vol(ix,iy)
              k2_tmp(iy,1)=k2_tmp(iy,1)+vol(ix,iy)
            endif
          endif
        enddo
        do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            if (on_closed_surface(ix,iy)) then
              k2_mean(ix,iy)=k2_tmp(iy,0)/k2_tmp(iy,1)
            endif
          endif
        enddo
      enddo

      do iy = 0, ny-1
        do ix = 0, nx-1
* coefficients are calculated only inside the separatrix and for plasma with low main ion collisionality
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            if (on_closed_surface(ix,iy) .and.
     &         (g_HS_style(ix,iy).le.5.0e0_R8)) then  !this is a temporary decision
              do is = 0, ns-1
c                if (Zhdanov_closure.eq.1) then

                if (rza(ix,iy,is).ge.corr_Z_min) then
* coefficients are calculated only for chosen ions of the chosen impurity
                  corr_fria(ix,iy,is)=
     &             1.0_R8+(-0.36_R8*alpha_HS_style(ix,iy,is)+
     &                      0.44_R8*alpha_HS_style(ix,iy,is)**(0.77_R8)+
     &                      0.01_R8)/(1+g_HS_style(ix,iy)**2)

                  temp_corr(ix,iy,0,is)=
     &           (chTAanxy(ix,iy,s2nucl(ismain)-1)*sqrt(2.0_R8))/
     &           (cRTAanxy(ix,iy,s2nucl(ismain)-1,s2nucl(is)-1)*2.5_R8)
                  temp_corr(ix,iy,1,is)=
     &             (2.4_R8*sqrt(eps1(iy)**3))/(1.0_R8+(0.36_R8+0.4_R8*
     &              alpha_HS_style(ix,iy,is))*g_HS_style(ix,iy)**2)
                  temp_corr(ix,iy,2,is)=
     &             (1.6_R8*(eps1(iy)**3)*g_HS_style(ix,iy)**2)/
     &             (1.0_R8+sqrt(eps1(iy)**3)*g_HS_style(ix,iy)**2)
                  corr_tfia(ix,iy,is)=
     &             (temp_corr(ix,iy,0,is)*temp_corr(ix,iy,1,is)+
     &              temp_corr(ix,iy,2,is))/k2_mean(ix,iy)
                endif

* the correction for main ions is calculated to make total forces zero
                if ((s2nucl(is)-1).eq.corr_spec) then

                  corr_fria_D(ix,iy)=corr_fria_D(ix,iy)+
     &             corr_fria(ix,iy,is)*rz2(ix,iy,is)*na(ix,iy,is)/
     &             z2n_xy(ix,iy,s2nucl(is)-1)
                  corr_tfia_D(ix,iy)=corr_tfia_D(ix,iy)+
     &             corr_tfia(ix,iy,is)*rz2(ix,iy,is)*na(ix,iy,is)/
     &             z2n_xy(ix,iy,s2nucl(is)-1)
                endif
              enddo
              corr_fria(ix,iy,ismain)=corr_fria_D(ix,iy)
              corr_tfia(ix,iy,ismain)=corr_tfia_D(ix,iy)

            endif
          endif
        enddo
      enddo

      if (iout.ne.0) then !  .or. iout_b2wdat.eq.4) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out(70,nx,ny,temp_corr(:,:,0,is),
     &     'b2mod_frtf_nccorr_coefs_temp0_'//chns)
          call my_out(70,nx,ny,temp_corr(:,:,1,is),
     &     'b2mod_frtf_nccorr_coefs_temp1_'//chns)
          call my_out(70,nx,ny,temp_corr(:,:,2,is),
     &     'b2mod_frtf_nccorr_coefs_temp2_'//chns)
        enddo
        call my_out(70,nx,ny,k2_mean(:,:),
     &     'b2mod_frtf_nccorr_coefs_k2mean_'//chns)
      endif

      return
      end subroutine b2mod_frtf_nccorr_coefs

      end module b2mod_frtf_nccorr
