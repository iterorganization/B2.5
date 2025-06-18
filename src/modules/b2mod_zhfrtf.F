      module b2mod_zhfrtf
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2zhco
      use b2mod_b2mrel
      use b2mod_subsys
      implicit none
      integer, allocatable, save ::
     &  is_i(:)
      real (kind=R8), allocatable, save ::
     &  z2n_xy(:,:), nal(:,:), Ia(:,:), av_ualpha(:,:),
     &  gt_Ac(:,:), gtalc(:,:), avm_u(:), rho_A_rel(:,:),
     &  gavm_uc(:), Z_to_m1_ast(:,:)
      real (kind=R8), allocatable, save ::
     &  w_out(:,:), hTdp_out(:,:), hwdp_out(:,:),
     &  rTdpgamma_out(:,:), rwdpgamma_out(:,:), piu_out(:,:)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF stores temperature gradients, average velocities and
*     other variables, which are used in the Zhdanov closure calculation.
*     It contains subroutines to calculate temperature gradients etc...,
*     as well as the friction and thermal forces.
*
*  3. description
*
*     This module contains:
*     z2n_xy: squared charge weighted average density of the isonuclear sequence
*     nal: average density of the isonuclear sequence
*     Ia: ratio of the squared charge weighted density of an ion to
*         the squared charge weighted average density of its isonuclear sequence
*     av_ualpha: Ia weighted average flow velocity of the isonuclear sequence
*     avm_u: mass-weighted average flow velocity for all ions
*     gtalc: effective temperature gradient of the isonuclear sequence
*           (cell-centred quantity)
*     gt_Ac: temperature gradient of each ion charge state
*           (cell-centred quantity)
*
*-----------------------------------------------------------------------
      contains
!...................................................................................
      subroutine alloc_b2mod_zhfrtf(nCv, ns)
      implicit none
      integer :: nCv, ns, is

      allocate(is_i(0:ns-1))
*   ..if ion sequence
      do is=0,ns-1
        if (.not.is_neutral(is)) then
          is_i(is)=1
        else
          is_i(is)=0
        endif
      enddo

      allocate(z2n_xy(nCv,0:nspecies-1))
      allocate(nal(nCv,0:nspecies-1))
      allocate(Ia(nCv,0:ns-1))
      allocate(av_ualpha(nCv,0:nspecies-1))
      allocate(avm_u(nCv))
      allocate(gtalc(nCv,0:nspecies-1))
      allocate(gt_Ac(nCv,0:ns-1))
      allocate(rho_A_rel(nCv,0:ns-1))
      allocate(gavm_uc(nCv))
      allocate(Z_to_m1_ast(nCv,0:nspecies-1))

      return
      end subroutine alloc_b2mod_zhfrtf

      subroutine dealloc_b2mod_zhfrtf
      implicit none
      if (.not.allocated(is_i)) return
      deallocate(is_i)
      deallocate(z2n_xy)
      deallocate(nal)
      deallocate(Ia)
      deallocate(av_ualpha)
      deallocate(avm_u)
      deallocate(gtalc)
      deallocate(gt_Ac)
      deallocate(rho_A_rel)
      deallocate(gavm_uc)
      deallocate(Z_to_m1_ast)
      return
      end subroutine dealloc_b2mod_zhfrtf

!...................................................................................
      subroutine alloc_b2mod_zhfrtf_df(nCv, ns)
      implicit none
      integer :: nCv, ns

      allocate(w_out(nCv,0:ns-1))
      allocate(hTdp_out(nCv,0:ns-1))
      allocate(hwdp_out(nCv,0:ns-1))
      allocate(rTdpgamma_out(nCv,0:ns-1))
      allocate(rwdpgamma_out(nCv,0:ns-1))
      allocate(piu_out(nCv,0:ns-1))

      return
      end subroutine alloc_b2mod_zhfrtf_df

      subroutine dealloc_b2mod_zhfrtf_df
      implicit none
      if (.not.allocated(w_out)) return
      deallocate(w_out)
      deallocate(hTdp_out)
      deallocate(hwdp_out)
      deallocate(rTdpgamma_out)
      deallocate(rwdpgamma_out)
      deallocate(piu_out)
      return
      end subroutine dealloc_b2mod_zhfrtf_df

!...................................................................................
      subroutine b2mod_zhfrtf_prep(nCv,ns,na,ua)
      implicit none
      integer :: nCv, ns
      real (kind=R8) :: na(nCv,0:ns-1),
     &  ua(nCv,0:ns-1)
*   ..local variables
      integer :: inucl, iz
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF_PREP calculates averaged densities, the Ia ratio and
*     the Ia weighted average velocity for each isonuclear sequence.
*
*  3. description
*
*     This routine computes:
*     z2n_xy: squared charge weighted average density
*             of the isonuclear sequence
*     nal: average density of the isonuclear sequence
*     Ia: ratio of the squared charge weighted density to
*         the squared charge weighted average density
*         of the isonuclear sequence
*     av_ualpha: Ia weighted average flow velocityi
*         of the isonuclear sequence
*
*-----------------------------------------------------------------------

*   ..calculate the average Z**2*n and numerical density of the isonuclear sequence
      z2n_xy = 0.0_R8
      nal = 0.0_R8
      do inucl = 0, nspecies-1
        do iz = 1, znnucl_int(inucl)
          z2n_xy(:,inucl) = z2n_xy(:,inucl) +
     &      iz**2*na(:,nucl2s(inucl,iz))
          nal(:,inucl) = nal(:,inucl) +
     &      na(:,nucl2s(inucl,iz))
        enddo
      enddo

*   ..Ia ratio and Ia weighted average velocity of the one type of ions
      Ia = 0.0_R8
      av_ualpha = 0.0_R8
      Z_to_m1_ast = 0.0_R8
      do inucl = 0, nspecies-1
        do iz= 1, znnucl_int(inucl)
          Ia(:,nucl2s(inucl,iz)) =
     &      iz**2*na(:,nucl2s(inucl,iz))/z2n_xy(:,inucl)
          av_ualpha(:,inucl) = av_ualpha(:,inucl) +
     &      Ia(:,nucl2s(inucl,iz))*ua(:,nucl2s(inucl,iz))
          Z_to_m1_ast(:,inucl) = Z_to_m1_ast(:,inucl) +
     &      na(:,nucl2s(inucl,iz))/(iz*nal(:,inucl))
        enddo
      enddo

      return
      end subroutine b2mod_zhfrtf_prep

!...................................................................................
      subroutine b2mod_zhfrtf_avm_u(nCv, ns, na, ua)
      use b2mod_constants
      implicit none
      integer nCv, ns

*.declarations
*   ..input arguments
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1)

*   ..local variables
      integer is
      real (kind=R8) :: rho(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF_AVM_U calculates the mass-weighted average flow
*     velocity of all ions.
*
*  3. description
*
*     This routine computes:
*     rho_A_rel: relative mass density of the ion "a" to
*                the total ion mass density
*     avm_u: mass-weighted average flow velocity of all ions
*
*-----------------------------------------------------------------------
*.computation

      rho = 0.0_R8
      avm_u = 0.0_R8
      do is = 0, ns-1
        rho_A_rel(:,is) = is_i(is)*am(is)*na(:,is)
        rho   = rho   + rho_A_rel(:,is)
        avm_u = avm_u + rho_A_rel(:,is)*ua(:,is)
      end do
      do is = 0, ns-1
        rho_A_rel(:,is) = rho_A_rel(:,is)/rho
      end do
      avm_u = avm_u / rho

      return
      end subroutine b2mod_zhfrtf_avm_u

!...................................................................................
      subroutine b2mod_zhfrtf_tgta(nCv, nFc, nVx, ns, switch, geo, mpg,
     &  na, ua, ti)
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2mod_switches
      implicit none

*.declarations
*   ..input arguments
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     *  na(nCv,0:ns-1), ua(nCv,0:ns-1), ti(nCv)

*   ..local variables
      integer is, iz, inucl
      real (kind=R8) :: rho(nCv), w_A(nCv,0:ns-1),
     &   t_A(nCv,0:ns-1), wrkv(nVx)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF_TGTA calculates temperature gradients and mass-weighted
*     average flow velocities for all ions.
*
*  3. description
*
*     This routine computes:
*     rho_A_rel: relative mass density of the ion "a" to the total ion mass density
*     avm_u: mass-weighted average flow velocity of all ions
*     gtalc: effective temperature gradient of the isonuclear sequence
*           (cell-centred quantity)
*     gt_Ac: temperature gradient of each ion charge state
*           (cell-centred quantity)
*
*-----------------------------------------------------------------------
*.computation

*   ..relative mass density and mass-weighted velocity
      rho = 0.0_R8
      avm_u = 0.0_R8
      do is = 0, ns-1
        rho_A_rel(:,is) = is_i(is)*am(is)*na(:,is)
        rho   = rho   + rho_A_rel(:,is)
        avm_u = avm_u + rho_A_rel(:,is)*ua(:,is)
      end do
      do is = 0, ns-1
        rho_A_rel(:,is) = rho_A_rel(:,is)/rho
      end do
      avm_u = avm_u / rho

*   ..temperature of each ion charge state
      do is = 0, ns-1
        w_A(:,is) = ua(:,is) - avm_u
        t_A(:,is) = ti + switch%zhd_corr*switch%zhd_corrT*
     &    (am(is)*mp*w_A(:,is)**2)/3.0_R8
      end do

*   ..temperature gradient of each ion charge state  in the cell center
      gt_Ac = 0.0_R8
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, t_A(:,is), wrkv,
     .      gt_Ac(:,is))
          gt_Ac(:,is) = gt_Ac(:,is)*geo%cvBb(:,0)/geo%cvBb(:,3)
        end if
      end do

*   ..effective temperature gradient of each isonuclear sequence in the cell center
      gtalc = 0.0_R8
      do inucl = 0, nspecies-1
        do iz = 1, znnucl_int(inucl)
          gtalc(:,inucl)= gtalc(:,inucl) +
     &     na(:,nucl2s(inucl,iz))/nal(:,inucl)*gt_Ac(:,nucl2s(inucl,iz))
        end do
      end do

*   ..average velocity gradient in the cell center
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, avm_u, wrkv, gavm_uc)
      gavm_uc = gavm_uc*geo%cvBb(:,0)/geo%cvBb(:,3)

      return
      end subroutine b2mod_zhfrtf_tgta

!...................................................................................
      subroutine b2mod_zhfrtf_fr(iCv,is,coef,smfr_p2,smfr_p1)
      implicit none
      integer :: iCv, is, inucl
      real (kind=R8) :: coef, smfr_p1, smfr_p2
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF_FR calculates the friction force.
*
*  3. description
*
*     This routine computes:
*     -smfr_p1: smbfrial(iCv,1) part of the friction force
*     -smfr_p2: smbfrial(iCv,0) part of the friction force
*
*-----------------------------------------------------------------------
*  ..local variabless
      real (kind=R8) :: temp

*   ..the first part of the friction force
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(iCv,inucl)*
     &   rmu(inucl,s2nucl(is)-1)*c_R_w(iCv,inucl,s2nucl(is)-1)
      enddo
      smfr_p1=-coef*temp

*   ..the second part of the friction force
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(iCv,inucl)*
     &   rmu(inucl,s2nucl(is)-1)*c_R_w(iCv,inucl,s2nucl(is)-1)*
     &   av_ualpha(iCv,inucl)
      enddo
      smfr_p2=coef*temp

      return
      end subroutine b2mod_zhfrtf_fr

!...................................................................................
      subroutine b2mod_zhfrtf_tf(iCv,is,zhcscorr,nais,smtf,smtf_nofl)
      implicit none
      integer :: iCv, is, zhcscorr
      real (kind=R8) :: smtf, smtf_nofl,nais
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ZHFRTF_TF calculates the thermal force.
*
*  3. description
*
*     This routine computes:
*     -smtf: the thermal force
*     -smtf_nofl: the thermal force without flux limiting
*
*-----------------------------------------------------------------------

*   ..the thermal force before charge state correction
      smtf=-Ia(iCv,is)*sum(nal(iCv,0:nspecies-1)*
     & c_R_TA(iCv,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(iCv,0:nspecies-1))

*   ..the thermal force before charge state correction without flux limiting
      smtf_nofl=-Ia(iCv,is)*sum(nal(iCv,0:nspecies-1)*
     & c_R_TA_nofl(iCv,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(iCv,0:nspecies-1))

*   ..apply the charge state correction
      if (zhcscorr.eq.1) then
        smtf=smtf-c_R_TB(iCv,s2nucl(is)-1)*
     &    nal(iCv,s2nucl(is)-1)*
     &   (nais/nal(iCv,s2nucl(is)-1)*gt_Ac(iCv,is)-
     &    Ia(iCv,is)*gtalc(iCv,s2nucl(is)-1))

        smtf_nofl=smtf_nofl-c_R_TB_nofl(iCv,s2nucl(is)-1)*
     &    nal(iCv,s2nucl(is)-1)*
     &   (nais/nal(iCv,s2nucl(is)-1)*gt_Ac(iCv,is)-
     &    Ia(iCv,is)*gtalc(iCv,s2nucl(is)-1))
      endif

      return
      end subroutine b2mod_zhfrtf_tf

      end module b2mod_zhfrtf
