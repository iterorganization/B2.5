      module b2mod_zhfrtf
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2zhco
      use b2mod_b2mrel 
!    & only: nucl2s, znnucl_int, rmu

      integer, allocatable ::
     & is_i(:)
      real (kind=R8), allocatable ::
     & z2n_xy(:,:,:),nal(:,:,:),Ia(:,:,:),av_ualpha(:,:,:),
     & gtalc(:,:,:)
      contains
!----------------------------------------------------------
      subroutine b2mod_zhfrtf_init(nx, ny, ns)
      integer :: nx, ny, ns, is
      
      allocate(is_i(0:ns-1))
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
      allocate(gtalc(-1:nx,-1:ny,0:nspecies-1))

      end subroutine b2mod_zhfrtf_init

!----------------------------------------------------------
      subroutine b2mod_zhfrtf_prep(nx,ny,ns,na,ua)
      integer :: nx, ny, ns
      real (kind=R8) :: na(-1:nx,-1:ny,0:ns-1),
     & ua(-1:nx,-1:ny,0:ns-1)
*   ..local vars
      integer :: ix,iy,is,inucl,iz

      do iy = -1, ny
        do ix = -1, nx
*   ..calculate the average Z**2*n
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
!           print *, 'z2n:'
!           print *,  z2n_xy(ix, iy,:)
!           print *, 'Ia:'
!           print *,  Ia(ix, iy,:)
!           print *, 'av_ualpha:'
!           print *,  av_ualpha(ix, iy,:)
        enddo
      enddo

      end subroutine b2mod_zhfrtf_prep

!----------------------------------------------------------
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
      real (kind=R8) :: rho_A(0:ns-1), rho, avm_u, w_A(0:ns-1),
     & t_A(-1:nx,-1:ny,0:ns-1),gti(-1:nx,-1:ny),
     & gt_A(-1:nx,-1:ny,0:ns-1), gt_Ac(-1:nx,-1:ny,0:ns-1)

*-----------------------------------------------------------------------
*.computation
!     is_i(0:ns-1)=.not.is_neutral(0:ns-1)
!     print *, is_neutral
!     print *, is_i
*   ..temperature of each species
      do iy = -1, ny
       do ix = -1, nx
         rho_A(0:ns-1)=is_i(0:ns-1)*am(0:ns-1)*na(ix,iy,0:ns-1)
!        print *,rho_A
         rho=sum(rho_A)
!        print *,rho
         avm_u=sum(rho_A(0:ns-1)*ua(ix,iy,0:ns-1))/rho
!        print *,avm_u
         w_A=ua(ix,iy,0:ns-1)-avm_u
!        print *,w_A
         t_A(ix,iy,0:ns-1)=ti(ix,iy)+
     & zhd_corr*(am(0:ns-1)*mp*w_A(0:ns-1)**2)/3.0_R8
!        print *,ti(ix,iy)
!        print *,tal
!        print *,zhd_corr
       enddo
      enddo

*   ..temperature gradient
!     do iy = -1, ny
!      do ix = -1, nx
!       if(leftix(ix,iy).eq.-2) then
!        gti(ix,iy) = 0.0_R8
!       else
!        gti(ix,iy) = (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))/
!    &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
!       endif
!      enddo
!     enddo
!     print *, 'gti:'
!     do, iy=-1,ny
!       print *, (gti(ix,iy), ix=-1,nx )
!     enddo

*   ..temperature gradient of each species
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
!     print *, 'gt_A:'
!     do is=0,ns-1
!      print *, 'is:',is
!      do, iy=-1,ny
!       print *, (gt_A(ix,iy,is), ix=-1,nx )
!      enddo
!     enddo

*   ..temperature gradient of each species in the cell center 
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
!      print *, 'gt_Ac:'
!      do is=0,ns-1
!       print *, 'is:',is
!       do, iy=-1,ny
!        print *, (gt_Ac(ix,iy,is), ix=-1,nx )
!       enddo
!      enddo

*   ..effective temperature gradient of each isonuclear type in the cell center 
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
!     print *, 'gtalc:'
!     do inucl=0,nspecies-1
!      print *, 'inucl:',inucl
!      do, iy=-1,ny
!       print *, (gtalc(ix,iy,inucl), ix=-1,nx )
!      enddo
!     enddo
      end subroutine b2mod_zhfrtf_tgta

*------------------------------------------------------------------------------
      subroutine b2mod_zhfrtf_fr(ix,iy,is,coef,smfr_p2,
     & smfr_p1)
      integer :: ix, iy, is, inucl
      real (kind=R8) :: coef, smfr_p1, smfr_p2
*  ..local vars
      real (kind=R8) :: temp

      temp=0
      do inucl = 0, nspecies-1
!     if (s2nucl(is)-1.ne.inucl) then
        temp=temp+z2n_xy(ix,iy,inucl)*
     & rmu(s2nucl(is)-1,inucl)*c_R_w(ix,iy,s2nucl(is)-1,inucl)
!     endif
      enddo
!     temp=temp+z2n_xy(ix,iy,s2nucl(is)-1)*
!    & rmu(s2nucl(is)-1,s2nucl(is)-1)
      smfr_p1=-coef*temp

      temp=0
      do inucl = 0, nspecies-1
!     if (s2nucl(is)-1.ne.inucl) then
        temp=temp+z2n_xy(ix,iy,inucl)*
     & rmu(s2nucl(is)-1,inucl)*c_R_w(ix,iy,s2nucl(is)-1,inucl)*
     & av_ualpha(ix,iy,inucl)
!     endif
      enddo
!     temp=temp+z2n_xy(ix,iy,s2nucl(is)-1)*
!    & rmu(s2nucl(is)-1,s2nucl(is)-1)*
!    & av_ualpha(ix,iy,s2nucl(is)-1)
      smfr_p2=coef*temp

      end subroutine b2mod_zhfrtf_fr 

      subroutine b2mod_zhfrtf_tf(ix,iy,is,smtf,smtf_nofl)
      integer :: ix, iy, is
      real (kind=R8) :: smtf, smtf_nofl

        smtf=-Ia(ix,iy,is)*sum(nal(ix,iy,0:nspecies-1)*
     & c_R_TA(ix,iy,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(ix,iy,0:nspecies-1))

        smtf_nofl=-Ia(ix,iy,is)*sum(nal(ix,iy,0:nspecies-1)*
     & c_R_TA_nofl(ix,iy,s2nucl(is)-1,0:nspecies-1)*
     & gtalc(ix,iy,0:nspecies-1))
      end subroutine b2mod_zhfrtf_tf
      end module b2mod_zhfrtf
