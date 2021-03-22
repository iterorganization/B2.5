      module b2mod_zhfrtf
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2zhco
      use b2mod_b2mrel 
      use b2mod_subsys
      implicit none
      integer, allocatable ::
     & is_i(:)
      real (kind=R8), allocatable ::
     & z2n_xy(:,:,:),nal(:,:,:),Ia(:,:,:),av_ualpha(:,:,:),
     & gt_Ac(:,:,:),gtalc(:,:,:),avm_u(:,:)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_zhfrtf store tempetrture gradients, average velocities and other variebles, which are used in the Zhdanov closure calculation 
*     b2mod_zhfrtf store subroutines to calculate temperature gradients and etc, the friction and thermal forces 
*
*  3. description 
*
*     This routine computes:
*     z2n_xy: squared charge weighted average density of the isonuclear ion
*     nal: average density of the isonuclear ion
*     Ia: squared charge weighted density to squared charge weighted average density of the isonuclear ion ratio
*     av_ualpha: Ia weighted average flow velocity of the isonuclear ion
*     avm_u: mass weighted average flow velocity of all ion
*     gtalc: effective temperature gradient of the isonuclear ion in the cell center
*     gt_Ac: temperature gradient of the ion each charge state in the cell center
*
*-----------------------------------------------------------------------
      contains
!...................................................................................
      subroutine b2mod_zhfrtf_init(nx, ny, ns)
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

      end subroutine b2mod_zhfrtf_init

!...................................................................................
      subroutine b2mod_zhfrtf_prep(nx,ny,ns,na,ua)
      implicit none
      integer :: nx, ny, ns
      real (kind=R8) :: na(-1:nx,-1:ny,0:ns-1),
     & ua(-1:nx,-1:ny,0:ns-1)
*   ..local vars
      integer :: ix,iy,is,inucl,iz
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_zhfrtf_prep calculates averaged densities, Ia ratio and Ia weighted average velocity of the isonuclear ion
*
*  3. description 
*
*     This routine computes:
*     z2n_xy: squared charge weighted average density of the isonuclear ion
*     nal: average density of the isonuclear ion
*     Ia: squared charge weighted density to squared charge weighted average density of the isonuclear ion ratio
*     av_ualpha: Ia weighted average flow velocity of the isonuclear ion
*
*-----------------------------------------------------------------------

      call subini ('b2mod_zhfrtf_prep')
      do iy = -1, ny
        do ix = -1, nx
*   ..calculate the average Z**2*n and numerical density of the isonuclear ion
          do inucl = 0, nspecies-1
            z2n_xy(ix,iy,inucl) = 0
            nal(ix,iy,inucl) = 0
            do iz= 1, znnucl_int(inucl)
              z2n_xy(ix,iy,inucl)=z2n_xy(ix,iy,inucl)+
     & iz**2*na(ix, iy, nucl2s(inucl,iz))
              nal(ix,iy,inucl)=nal(ix,iy,inucl)+
     & na(ix, iy, nucl2s(inucl,iz))
            enddo
          enddo
*   ..Ia ratio and Ia weighted average velocity of the isonuclear ion
          Ia(ix, iy, :)=0
          do inucl = 0, nspecies-1
            av_ualpha(ix,iy,inucl) = 0
            do iz= 1, znnucl_int(inucl)
              Ia(ix, iy, nucl2s(inucl,iz))=
     & iz**2*na(ix, iy, nucl2s(inucl,iz))/z2n_xy(ix,iy,inucl)
              av_ualpha(ix,iy,inucl)=av_ualpha(ix,iy,inucl)+
     & Ia(ix, iy, nucl2s(inucl,iz))*ua(ix, iy, nucl2s(inucl,iz))
            enddo
          enddo
        enddo
      enddo

      call subend ()
      end subroutine b2mod_zhfrtf_prep

!...................................................................................
      subroutine b2mod_zhfrtf_tgta(nx, ny, ns, na, ua, ti, hx, bb)
      use b2mod_constants
      use b2mod_indirect
      implicit none
      integer nx, ny, ns

*   ..input arguments 
      real (kind=R8) :: 
     * na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny),
     * hx(-1:nx,-1:ny), 
     * bb(-1:nx,-1:ny,0:3)

*.declarations
*   ..local variables
      integer ix, iy, is,iz,inucl
      real (kind=R8) :: rho_A(0:ns-1), rho, w_A(0:ns-1),
     & t_A(-1:nx,-1:ny,0:ns-1),gti(-1:nx,-1:ny),
     & gt_A(-1:nx,-1:ny,0:ns-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_zhfrtf_tgta calculates temperature gradients and mass weighted average flow velocity of all ion
*
*  3. description 
*
*     This routine computes:
*     avm_u: mass weighted average flow velocity of all ion
*     gtalc: effective temperature gradient of the isonuclear ion in the cell center
*     gt_Ac: temperature gradient of the ion each charge state in the cell center
*
*-----------------------------------------------------------------------
*.computation
      call subini ('b2mod_zhfrtf_tgta')
*   ..temperature of the ion each charge state
      do iy = -1, ny
       do ix = -1, nx
         rho_A(0:ns-1)=is_i(0:ns-1)*am(0:ns-1)*na(ix,iy,0:ns-1)
         rho=sum(rho_A)
         avm_u(ix,iy)=sum(rho_A(0:ns-1)*ua(ix,iy,0:ns-1))/rho
         w_A=ua(ix,iy,0:ns-1)-avm_u(ix,iy)
         t_A(ix,iy,0:ns-1)=ti(ix,iy)+
     & zhd_corr*(am(0:ns-1)*mp*w_A(0:ns-1)**2)/3.0_R8
       enddo
      enddo

*   ..temperature gradient of the ion each charge state
      gt_A = 0.0_R8
      do is=0,ns-1
       if (.not.is_neutral(is)) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
           gt_A(ix,iy,is) = 
     & (t_A(ix,iy,is)-t_A(leftix(ix,iy),leftiy(ix,iy),is))/
     & (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
          endif
         enddo
        enddo
       endif
      enddo

*   ..temperature gradient of the ion each charge state in the cell center 
      gt_Ac = 0.0_R8
      do is=0,ns-1
       if (.not.is_neutral(is)) then
        do iy = 0, ny-1
         do ix = 0, nx-1
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           gt_Ac(ix,iy,is) = bb(ix,iy,0)/bb(ix,iy,3)* 
     & (gt_A(ix,iy,is)+gt_A(rightix(ix,iy),rightiy(ix,iy),is))*0.5e0_R8
          endif
         enddo
        enddo
       endif
      enddo

*   ..effective temperature gradient of each isonuclear ion in the cell center 
      gtalc = 0.0_R8
      do iy = 0, ny-1
        do ix = 0, nx-1
          do inucl = 0, nspecies-1
            do iz= 1, znnucl_int(inucl)
              gtalc(ix,iy,inucl)=gtalc(ix,iy,inucl)+
     & na(ix, iy, nucl2s(inucl,iz))/nal(ix, iy, inucl)*
     & gt_Ac(ix,iy,nucl2s(inucl,iz))
            enddo
          enddo
        enddo
      enddo
      call subend ()
      end subroutine b2mod_zhfrtf_tgta

!...................................................................................
      subroutine b2mod_zhfrtf_fr(ix,iy,is,coef,smfr_p2,
     & smfr_p1)
      implicit none
      integer :: ix, iy, is, inucl
      real (kind=R8) :: coef, smfr_p1, smfr_p2
*  ..local vars
      real (kind=R8) :: temp
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_zhfrtf_fr calculates the friction force 
*
*  3. description 
*
*     This routine computes:
*     -smfr_p1: smbfrial(ix,iy,1) part of the friction force
*     -smfr_p2: smbfrial(ix,iy,0) part of the friction force
*
*-----------------------------------------------------------------------

      call subini ('b2mod_zhfrtf_fr')
*   ..the first part of the friction force 
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(ix,iy,inucl)*
     & rmu(s2nucl(is)-1,inucl)*c_R_w(ix,iy,s2nucl(is)-1,inucl)
      enddo
      smfr_p1=-coef*temp

*   ..the second part of the friction force 
      temp=0
      do inucl = 0, nspecies-1
        temp=temp+z2n_xy(ix,iy,inucl)*
     & rmu(s2nucl(is)-1,inucl)*c_R_w(ix,iy,s2nucl(is)-1,inucl)*
     & av_ualpha(ix,iy,inucl)
      enddo
      smfr_p2=coef*temp

      call subend ()
      end subroutine b2mod_zhfrtf_fr 

!...................................................................................
      subroutine b2mod_zhfrtf_tf(ix,iy,is,nais,smtf,smtf_nofl)
      implicit none
      integer :: ix, iy, is
      real (kind=R8) :: smtf, smtf_nofl,nais
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2mod_zhfrtf_tf calculates the thermal force 
*
*  3. description 
*
*     This routine computes:
*     -smtf: the thermal force
*     -smtf_nofl: the thermal force without flux limiting
*
*-----------------------------------------------------------------------

      call subini ('b2mod_zhfrtf_tf')
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
     & nal(ix,iy,s2nucl(is)-1)*
     & (nais/nal(ix,iy,s2nucl(is)-1)*gt_Ac(ix,iy,is)-
     & Ia(ix,iy,is)*gtalc(ix,iy,s2nucl(is)-1))

         smtf_nofl=smtf_nofl-c_R_TB_nofl(ix,iy,s2nucl(is)-1)*
     & nal(ix,iy,s2nucl(is)-1)*
     & (nais/nal(ix,iy,s2nucl(is)-1)*gt_Ac(ix,iy,is)-
     & Ia(ix,iy,is)*gtalc(ix,iy,s2nucl(is)-1))
       endif

      call subend ()
      end subroutine b2mod_zhfrtf_tf
      end module b2mod_zhfrtf
