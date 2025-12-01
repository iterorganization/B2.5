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
     &  z2n_xy(:,:,:), nal(:,:,:), Ia(:,:,:), av_ualpha(:,:,:),
     &  gt_Ac(:,:,:), gtalc(:,:,:), avm_u(:,:), rho_A_rel(:,:,:),
     &  gavm_uc(:,:), Z_to_m1_ast(:,:,:)
      real (kind=R8), allocatable, save ::
     &  w_out(:,:,:), hTdp_out(:,:,:), hwdp_out(:,:,:),
     &  rTdpgamma_out(:,:,:), rwdpgamma_out(:,:,:), piu_out(:,:,:)
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
      subroutine alloc_b2mod_zhfrtf(nx, ny, ns)
      implicit none
      integer :: nx, ny, ns, is

      allocate(is_i(0:ns-1))
*   ..if ion sequence
      do is=0,ns-1
        if (.not.is_neutral(is)) then
          is_i(is)=1
        else
          is_i(is)=0
        endif
      enddo

      allocate(z2n_xy(-1:nx,-1:ny,0:nspecies-1))
      allocate(nal(-1:nx,-1:ny,0:nspecies-1))
      allocate(Ia(-1:nx,-1:ny,0:ns-1))
      allocate(av_ualpha(-1:nx,-1:ny,0:nspecies-1))
      allocate(avm_u(-1:nx,-1:ny))
      allocate(gtalc(-1:nx,-1:ny,0:nspecies-1))
      allocate(gt_Ac(-1:nx,-1:ny,0:ns-1))
      allocate(rho_A_rel(-1:nx,-1:ny,0:ns-1))
      allocate(gavm_uc(-1:nx,-1:ny))
      allocate(Z_to_m1_ast(-1:nx,-1:ny,0:nspecies-1))

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
      subroutine alloc_b2mod_zhfrtf_df(nx, ny, ns)
      implicit none
      integer :: nx, ny, ns

      allocate(w_out(-1:nx,-1:ny,0:ns-1))
      allocate(hTdp_out(-1:nx,-1:ny,0:ns-1))
      allocate(hwdp_out(-1:nx,-1:ny,0:ns-1))
      allocate(rTdpgamma_out(-1:nx,-1:ny,0:ns-1))
      allocate(rwdpgamma_out(-1:nx,-1:ny,0:ns-1))
      allocate(piu_out(-1:nx,-1:ny,0:ns-1))

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
      subroutine b2mod_zhfrtf_prep(nx,ny,ns,na,ua)
      implicit none
      integer :: nx, ny, ns
      real (kind=R8) :: na(-1:nx,-1:ny,0:ns-1),
     &  ua(-1:nx,-1:ny,0:ns-1)
*   ..local variables
      integer :: ix, iy, inucl, iz
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

      do iy = -1, ny
        do ix = -1, nx
*   ..calculate the average Z**2*n and numerical density of the isonuclear sequence
          do inucl = 0, nspecies-1
            z2n_xy(ix,iy,inucl) = 0.0_R8
            nal(ix,iy,inucl) = 0.0_R8
            do iz= 1, znnucl_int(inucl)
              z2n_xy(ix,iy,inucl)=z2n_xy(ix,iy,inucl)+
     &         iz**2*na(ix,iy,nucl2s(inucl,iz))
              nal(ix,iy,inucl)=nal(ix,iy,inucl)+
     &         na(ix,iy,nucl2s(inucl,iz))
            enddo
          enddo
*   ..Ia ratio and Ia weighted average velocity of the one type of ions
          Ia(ix,iy,:)=0.0_R8
          do inucl = 0, nspecies-1
            av_ualpha(ix,iy,inucl) = 0.0_R8
            Z_to_m1_ast(ix,iy,inucl) = 0.0_R8
            do iz= 1, znnucl_int(inucl)
              Ia(ix,iy,nucl2s(inucl,iz))=
     &         iz**2*na(ix,iy,nucl2s(inucl,iz))/z2n_xy(ix,iy,inucl)
              av_ualpha(ix,iy,inucl)=av_ualpha(ix,iy,inucl)+
     &         Ia(ix,iy,nucl2s(inucl,iz))*ua(ix,iy,nucl2s(inucl,iz))
              Z_to_m1_ast(ix,iy,inucl)=Z_to_m1_ast(ix,iy,inucl)+
     &         na(ix,iy,nucl2s(inucl,iz))/(iz*nal(ix,iy,inucl))
            enddo
          enddo
        enddo
      enddo

      return
      end subroutine b2mod_zhfrtf_prep

!...................................................................................
      subroutine b2mod_zhfrtf_avm_u(nx, ny, ns, na, ua)
      use b2mod_constants
      use b2mod_indirect
      implicit none
      integer nx, ny, ns

*.declarations
*   ..input arguments
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1)

*   ..local variables
      integer ix, iy
      real (kind=R8) :: rho_A(0:ns-1), rho

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
      do iy = -1, ny
       do ix = -1, nx
         rho_A(0:ns-1)=is_i(0:ns-1)*am(0:ns-1)*na(ix,iy,0:ns-1)
         rho=sum(rho_A)
         rho_A_rel(ix,iy,0:ns-1)=rho_A(0:ns-1)/rho
         avm_u(ix,iy)=sum(rho_A(0:ns-1)*ua(ix,iy,0:ns-1))/rho
       enddo
      enddo

      return
      end subroutine b2mod_zhfrtf_avm_u

!...................................................................................
      subroutine b2mod_zhfrtf_tgta(nx, ny, ns, na, ua, ti, hx, bb)
      use b2mod_constants
      use b2mod_indirect
      implicit none

*.declarations
*   ..input arguments
      integer nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny),
     *  hx(-1:nx,-1:ny),
     *  bb(-1:nx,-1:ny,0:3)

*   ..local variables
      integer ix, iy, is, iz, inucl
      real (kind=R8) :: rho_A(0:ns-1), rho, w_A(0:ns-1),
     &   t_A(-1:nx,-1:ny,0:ns-1),
     &  gt_A(-1:nx,-1:ny,0:ns-1), gavm_u(-1:nx,-1:ny)

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

*   ..temperature of each ion charge state
      do iy = -1, ny
       do ix = -1, nx
         rho_A(0:ns-1)=is_i(0:ns-1)*am(0:ns-1)*na(ix,iy,0:ns-1)
         rho=sum(rho_A)
         rho_A_rel(ix,iy,0:ns-1)=rho_A(0:ns-1)/rho
         avm_u(ix,iy)=sum(rho_A(0:ns-1)*ua(ix,iy,0:ns-1))/rho
         w_A=ua(ix,iy,0:ns-1)-avm_u(ix,iy)
         t_A(ix,iy,0:ns-1)=ti(ix,iy)+
     &    zhd_corr*zhd_corrT*(am(0:ns-1)*mp*w_A(0:ns-1)**2)/3.0_R8
       enddo
      enddo

*   ..temperature gradient of each ion charge state
      gt_A = 0.0_R8
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
           gt_A(ix,iy,is) =
     &     (t_A(ix,iy,is)-t_A(leftix(ix,iy),leftiy(ix,iy),is))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          endif
         enddo
        enddo
       endif
      enddo

*   ..temperature gradient of each ion charge state in the cell center
      gt_Ac = 0.0_R8
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           gt_Ac(ix,iy,is) = 0.5e0_R8*bb(ix,iy,0)/bb(ix,iy,3)*
     &     (gt_A(ix,iy,is)+gt_A(rightix(ix,iy),rightiy(ix,iy),is))
          endif
         enddo
        enddo
       endif
      enddo

*   ..effective temperature gradient of each isonuclear sequence in the cell center
      gtalc = 0.0_R8
      do iy = 0, ny-1
        do ix = 0, nx-1
          if (leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
            do inucl = 0, nspecies-1
              do iz = 1, znnucl_int(inucl)
                gtalc(ix,iy,inucl)=gtalc(ix,iy,inucl)+
     &           na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl)*
     &           gt_Ac(ix,iy,nucl2s(inucl,iz))
              enddo
            enddo
          endif
        enddo
        do ix = -1, nx
          if (leftix(ix,iy).eq.-2.and.rightix(ix,iy).ne.nx+1) then
            do inucl = 0, nspecies-1
              gtalc(ix,iy,inucl)=
     &         gtalc(rightix(ix,iy),rightiy(ix,iy),inucl)
            enddo
          else if (rightix(ix,iy).eq.nx+1.and.leftix(ix,iy).ne.-2) then
            do inucl = 0, nspecies-1
              gtalc(ix,iy,inucl)=
     &         gtalc(leftix(ix,iy),leftiy(ix,iy),inucl)
            enddo
          endif
        enddo
      enddo
      do ix = -1, nx
        do inucl = 0, nspecies-1
          gtalc(ix,-1,inucl)=
     &     gtalc(topix(ix,-1),topiy(ix,-1),inucl)
          gtalc(ix,ny,inucl)=
     &     gtalc(bottomix(ix,ny),bottomiy(ix,ny),inucl)
        enddo
      enddo

*   ..average velocity gradient
      gavm_u = 0.0_R8
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
           gavm_u(ix,iy) =
     &     (avm_u(ix,iy)-avm_u(leftix(ix,iy),leftiy(ix,iy)))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          endif
        enddo
      enddo

*   ..average velocity gradient in the cell center
      gavm_uc = 0.0_R8
      do iy = 0, ny-1
        do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            gavm_uc(ix,iy) = 0.5e0_R8*bb(ix,iy,0)/bb(ix,iy,3)*
     &      (gavm_u(ix,iy)+gavm_u(rightix(ix,iy),rightiy(ix,iy)))
          endif
        enddo
        do ix = -1, nx
          if (leftix(ix,iy).eq.-2.and.rightix(ix,iy).ne.nx+1) then
            gavm_uc(ix,iy) = gavm_uc(rightix(ix,iy),rightiy(ix,iy))
          else if (rightix(ix,iy).eq.nx+1.and.leftix(ix,iy).ne.-2) then
            gavm_uc(ix,iy) = gavm_uc(leftix(ix,iy),leftiy(ix,iy))
          endif
        enddo
      enddo
      do ix = -1, nx
        gavm_uc(ix,-1) = gavm_uc(topix(ix,-1),topiy(ix,-1))
        gavm_uc(ix,ny) = gavm_uc(bottomix(ix,ny),bottomiy(ix,ny))
      enddo

      return
      end subroutine b2mod_zhfrtf_tgta

!...................................................................................
      subroutine b2mod_zhfrtf_fr(ix,iy,is,coef,smfr_p2,smfr_p1,
     &  corr_fria)
      implicit none
      integer :: ix, iy, is, inucl
      real (kind=R8) :: coef, smfr_p1, smfr_p2, corr_fria
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
*     -smfr_p1: smbfrial(ix,iy,1) part of the friction force
*     -smfr_p2: smbfrial(ix,iy,0) part of the friction force
*
*-----------------------------------------------------------------------
*  ..local variabless
      real (kind=R8) :: temp

*   ..the first part of the friction force
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(ix,iy,inucl)*
     &   rmu(inucl,s2nucl(is)-1)*c_R_w(ix,iy,inucl,s2nucl(is)-1)
      enddo
      smfr_p1=-coef*temp*corr_fria

*   ..the second part of the friction force
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(ix,iy,inucl)*
     &   rmu(inucl,s2nucl(is)-1)*c_R_w(ix,iy,inucl,s2nucl(is)-1)*
     &   av_ualpha(ix,iy,inucl)
      enddo
      smfr_p2=coef*temp*corr_fria

      return
      end subroutine b2mod_zhfrtf_fr

!...................................................................................
      subroutine b2mod_zhfrtf_tf(ix,iy,is,nais,smtf,smtf_nofl,corr_tfia)
      implicit none
      integer :: ix, iy, is
      real (kind=R8) :: smtf, smtf_nofl,nais, corr_tfia
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
      smtf=-Ia(ix,iy,is)*sum(nal(ix,iy,0:nspecies-1)*
     & c_R_TA(ix,iy,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(ix,iy,0:nspecies-1))

*   ..the thermal force before charge state correction without flux limiting
      smtf_nofl=-Ia(ix,iy,is)*sum(nal(ix,iy,0:nspecies-1)*
     & c_R_TA_nofl(ix,iy,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(ix,iy,0:nspecies-1))

*   ..apply the charge state correction
      if (zhcscorr.eq.1) then
        smtf=smtf-c_R_TB(ix,iy,s2nucl(is)-1)*
     &    nal(ix,iy,s2nucl(is)-1)*
     &   (nais/nal(ix,iy,s2nucl(is)-1)*gt_Ac(ix,iy,is)-
     &    Ia(ix,iy,is)*gtalc(ix,iy,s2nucl(is)-1))

        smtf_nofl=smtf_nofl-c_R_TB_nofl(ix,iy,s2nucl(is)-1)*
     &    nal(ix,iy,s2nucl(is)-1)*
     &   (nais/nal(ix,iy,s2nucl(is)-1)*gt_Ac(ix,iy,is)-
     &    Ia(ix,iy,is)*gtalc(ix,iy,s2nucl(is)-1))
      endif

      smtf = smtf*corr_tfia
      smtf_nofl = smtf_nofl*corr_tfia
      return
      end subroutine b2mod_zhfrtf_tf

      end module b2mod_zhfrtf
