*-----------------------------------------------------------------------
*.specification

      subroutine b2trzh (nx, ny, ns, znmax, na, ua, ti, lnlam, hx, bb,
     & vol,gs,qz,mode,gamma, chci0_nofl,chci0,chci0_luc)  !som 05.02.2021  
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants     
      use b2mod_subsys
      use b2mod_b2mrel
      use b2mod_b2zhco
      use b2mod_zhfrtf
      use b2mod_transport , only: f_luc_ki
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, znmax, mode
      real (kind=R8) :: gamma                             !som 20.03.21
      real (kind=R8) ::
     &  na(-1:nx,-1:ny,0:ns-1), ti(-1:nx,-1:ny), ua(-1:nx,-1:ny,0:ns-1),
     &  lnlam(-1:nx,-1:ny), hx(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  qz(-1:nx,-1:ny,0:1)

*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  chci0(-1:nx,-1:ny,0:1), chci0_nofl(-1:nx,-1:ny,0:1),
     &  chci0_luc(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRZH computes the transport coefficient using 21N-moment approach, which  
*     is discussed in [V. M. Zhdanov. Transport Processes in Multicomponent Plasma,
*     vol. 44. 10 2002.]
*
*  3. description 
*
*     This routine computes:
*     - c_R_TA: the thermal force transport coefficient for isonuclear ions
*     - c_R_TA_nofl: the thermal force transport coefficient for isonuclear ions without flux limiting
*     - c_R_TB: the thermal force transport coefficient correction for each charge state
*     - c_R_TB_nofl: the thermal force transport coefficient correction for each charge state without flux limiting
*     - c_R_wA: the friction force transport coefficient for isonuclear ions
*     - c_R_wB: the friction force transport coefficient correction for each charge state
*     - c_R_w: the total friction force transport coefficient 
*     - c_hwA: the velocity dependent ion heat flux tansport coefficient for isonuclear ions
*     - temperature gradients  
*     - mass-related variables  
*
*-----------------------------------------------------------------------
*.declarations
      integer ncall
      integer ix, iy, is, inucl, iz,i, j 
      real (kind=R8) :: temp
      real (kind=R8) :: zn_temp, am_temp, amti(0:nspecies-1),
     & zeffne, nz2n(0:nspecies-1),
     & G2(0:nspecies-1,0:nspecies-1), G8(0:nspecies-1,0:nspecies-1),
     & Ainv(0:2*nspecies-1,0:2*nspecies-1),
     & Ainv_mdf(0:nspecies-1,0:2*nspecies-1),
     & Ainv_fl(0:2*nspecies-1,0:2*nspecies-1),
     & Ainv_fl_h(0:2*nspecies-1,0:2*nspecies-1),
     & Ainv_mdf_fl(0:nspecies-1,0:2*nspecies-1),
     & c_rwA(0:nspecies-1,0:nspecies-1),
     & c_R_TA_an_nofl(0:nspecies-1,0:nspecies-1),
     & c_R_TA_an(0:nspecies-1,0:nspecies-1), c_rTA_an(0:nspecies-1),
     & c_hTA_an(0:nspecies-1), c_hTA_an_fl(0:nspecies-1),
     & t0, t1, czetap,zetap(-1:nx,-1:ny),
     & cF_limh(0:nspecies-1),
     & cF_lim(0:nspecies-1),cF_limB(0:nspecies-1),
     & hci0_nofl(-1:nx,-1:ny,0:1),hci0(-1:nx,-1:ny,0:1),
     & c_hw(0:nspecies-1,0:nspecies-1),c_rw(0:nspecies-1,0:nspecies-1),
     & c_hwB(0:nspecies-1),c_rwB(0:nspecies-1),
     & c_hTB(0:nspecies-1),c_rTB(0:nspecies-1),
     & kappa_hTA(0:nspecies-1,0:nspecies-1),kappa_hTB(0:nspecies-1),
     & kappa_rTA(0:nspecies-1,0:nspecies-1),kappa_rTB(0:nspecies-1)

      save ncall
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trzh')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
*   ..allocate
        call init_mrel(nspecies,znmax,ns)
        call init_zhco(nspecies, nx, ny)
        call b2mod_zhfrtf_init(nx, ny, ns)
        if (zh_dist_fun.eq.1) then
          call b2mod_zhfrtf_init_df(nx, ny, ns)
        endif

*   ..calculate the atom number and the atomic mass of each type of ions
        inucl=1
        znnucl(0)=zn(0)
        amnucl(0)=am(0)         
        zn_temp=zn(0)
        am_temp=am(0)
        do is = 0, ns-1
          if (zn(is).ne.zn_temp .and. am(is).ne.am_temp) then
           inucl=inucl+1 
           znnucl(inucl-1)=zn(is)
           amnucl(inucl-1)=am(is)  
           zn_temp=zn(is)
           am_temp=am(is)
          endif 
*   ..correspondence matrix
          nucl2s(inucl-1,nint(zamax(is)))=is 
          s2nucl(is)=inucl
          znnucl_int(inucl-1)=nint(znnucl(inucl-1))
        enddo        
*   ..calculate the N-coefficients and all mass related variables 
        call b2trcn (nspecies, amnucl, N2, N5, N6, N8, N9, N10, 
     & N11, N12, mu, rmu, rmd2, mum, mdm)
      endif
       
*   ..main part  
*   ..calculate the average Z**2*n and density of each type
      call b2mod_zhfrtf_prep(nx,ny,ns,na,ua)

*   ..calculate temperature gradients  
      call b2mod_zhfrtf_tgta(nx, ny, ns, na, ua, ti, hx, bb)

*   ..compute zetap 
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam

*   ..calculate transport coefficietns  
      do iy = -1, ny
        do ix = -1, nx
          amti=amnucl*mp/ti(ix,iy)
          zeffne=sum(z2n_xy(ix,iy,0:nspecies-1)) 
*    ..normalized average Z**2*n
          nz2n=z2n_xy(ix,iy,0:nspecies-1)/zeffne
*    ..calculate the correction for each charge state
          call b2trcc(nspecies,mum,mdm,nz2n,
     & c_R_TB_nofl(ix,iy,0:nspecies-1),
     & c_R_wB(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_hwB,c_rwB,c_hTB,c_rTB)

          if (zh_imp_analyt.eq.0) then
*    ..calculate the inversed matrix for heat flux/r-moment
            call b2trca(nspecies, amti, N5, N6, N9,  
     & N10, N11, N12, nz2n, Ainv)

*    ..re-normalize part of the inversed matrix  
            do j = 0, 2*nspecies-1
              Ainv_mdf(0:nspecies-1,j)=
     & amti(0:nspecies-1)*Ainv(nspecies:2*nspecies-1,j)
            enddo

*    ..Flux limiting coefficient	
            call b2trfl(ix,iy,Ainv,zeffne,amti,zetap(ix,iy),gamma,
     & cF_lim,cF_limh,cF_limB)
            do i = 0, nspecies-1
              Ainv_fl(i,0:nspecies-1)=
     & cF_lim(i)*Ainv(i,0:nspecies-1)
              Ainv_mdf_fl(i,0:nspecies-1)=
     & cF_lim(i)*Ainv_mdf(i,0:nspecies-1)
              Ainv_fl_h(i,0:nspecies-1)=
     & cF_limh(i)*Ainv(i,0:nspecies-1)
            enddo
*    ..G-objects for the thermal and friction forces	 
            do j = 0, nspecies-1
              G2(:,j)=nz2n(:)*nz2n(j)*N2(:,j)
              G8(:,j)=nz2n(:)*nz2n(j)*N8(:,j)
            enddo

*    ..calculate the classical ion conductivity without flux limiting 
            call b2trht(nspecies,ix,iy,bb(ix,iy,0)/bb(ix,iy,3),
     & zetap(ix,iy),zeffne,ti(ix,iy), Ainv, 
     & hci0_nofl(ix,iy,0))
*    ..calculate the classical ion conductivity with flux limiting 
            call b2trht(nspecies,ix,iy,bb(ix,iy,0)/bb(ix,iy,3),
     & zetap(ix,iy),zeffne,ti(ix,iy), Ainv_fl_h, 
     & hci0(ix,iy,0))
*    ..calculate the velocity parts heat flux and r transport coefficient 
            call b2trhw(nspecies, G2, G8, mum, rmu,
     & z2n_xy(ix,iy,0:nspecies-1), Ainv, Ainv_mdf, 
     & c_hwA(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_rwA(0:nspecies-1,0:nspecies-1))
*    ..calculate the thermal force transport coefficient without flux limiting 
            call b2trtf(nspecies, G2, G8, mum, mdm,
     & Ainv, Ainv_mdf, c_R_TA_nofl(ix,iy,0:nspecies-1,0:nspecies-1))
*    ..calculate the thermal force transport coefficient with flux limiting 
            call b2trtf(nspecies, G2, G8, mum, mdm,
     & Ainv_fl, Ainv_mdf_fl, c_R_TA(ix,iy,0:nspecies-1,0:nspecies-1))

*    ..for the distribution function 
            if (zh_dist_fun.eq.1) then
              do i = 0, nspecies-1
                kappa_hTA(i,0:nspecies-1)=
     &   2.5_R8*(zetap(ix,iy)/(zeffne**2*mp))*
     &   Ainv(i,0:nspecies-1)*nal(ix,iy,0:nspecies-1)
                kappa_hTB(i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_hTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)

                kappa_rTA(i,0:nspecies-1)=
     &   2.5_R8*(zetap(ix,iy)/(zeffne**2*mp))*
     &   Ainv_mdf(i,0:nspecies-1)*nal(ix,iy,0:nspecies-1)
                kappa_rTB(i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_rTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
              enddo
            endif

          else
*    ..the analytical thermal force transport coefficient  
            call b2trhrtf_an(nspecies, mum, mdm,nz2n,
     & c_hwA(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_rwA(0:nspecies-1,0:nspecies-1),  
     & c_hTA_an,c_rTA_an,c_R_TA_an_nofl)
*    ..calculate the classical ion conductivity without flux limiting 
            call b2trht_an(nspecies,ix,iy,bb(ix,iy,0)/bb(ix,iy,3),
     & zetap(ix,iy),zeffne,ti(ix,iy),c_hTA_an,nz2n,rmd2,  
     & hci0_nofl(ix,iy,0))
*    ..the thermal force transport coefficient without flux limiting
            do i = 0, nspecies-1
              temp=sum(c_R_TA_an_nofl(i,0:nspecies-1)*
     & (nz2n(0:nspecies-1)/nz2n(i)))
              c_R_TA_nofl(ix,iy,i,0:nspecies-1)=
     & -c_R_TA_an_nofl(0:nspecies-1,i)*(nz2n(i)/nz2n(0:nspecies-1))
              c_R_TA_nofl(ix,iy,i,i)=c_R_TA_nofl(ix,iy,i,i)+temp
            enddo

*    ..Flux limiting coefficient (analytical)	
            call b2trfl_an(ix,iy,c_hTA_an,
     & amnucl,amti,zetap(ix,iy),gamma,
     & cF_lim,cF_limh,cF_limB)
            do i = 0, nspecies-1
              c_R_TA_an(i,0:nspecies-1)=
     & cF_lim(i)*c_R_TA_an_nofl(i,0:nspecies-1)
              c_hTA_an_fl(i)=
     & cF_limh(i)*c_hTA_an(i)
            enddo
*    ..calculate the classical ion conductivity with flux limiting 
            call b2trht_an(nspecies,ix,iy,bb(ix,iy,0)/bb(ix,iy,3),
     & zetap(ix,iy),zeffne,ti(ix,iy),c_hTA_an_fl,nz2n,rmd2,  
     & hci0(ix,iy,0))
*    ..the thermal force transport coefficient with flux limiting
            do i = 0, nspecies-1
              temp=sum(c_R_TA_an(i,0:nspecies-1)*
     & (nz2n(0:nspecies-1)/nz2n(i)))
              c_R_TA(ix,iy,i,0:nspecies-1)=
     & -c_R_TA_an(0:nspecies-1,i)*(nz2n(i)/nz2n(0:nspecies-1))
              c_R_TA(ix,iy,i,i)=c_R_TA(ix,iy,i,i)+temp
            enddo

*    ..for the distribution function 
            if (zh_dist_fun.eq.1) then
              kappa_hTA=0.0_R8
              kappa_rTA=0.0_R8
              do i = 0, nspecies-1
                kappa_hTA(i,i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_hTA_an(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
                kappa_hTB(i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_hTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)

                kappa_rTA(i,i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_rTA_an(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
                kappa_rTB(i)=
     &   -(zetap(ix,iy)/(zeffne**2*mp))*
     & (c_rTB(i)/rmd2(i)*nal(ix,iy,i)/nz2n(i)**2)
              enddo
            endif

          endif

*    ..calculate the thermal force correction for each charge state
          c_R_TB(ix,iy,0:nspecies-1)=cF_limB(0:nspecies-1)*
     & c_R_TB_nofl(ix,iy,0:nspecies-1)

*    ..calculate the friction force transport coefficient 
          call b2trfr(nspecies, mum, rmu,
     & z2n_xy(ix,iy,0:nspecies-1), 
     & c_hwA(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_rwA(0:nspecies-1,0:nspecies-1),
     & c_R_wB(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_R_wA(ix,iy,0:nspecies-1,0:nspecies-1),
     & c_R_w(ix,iy,0:nspecies-1,0:nspecies-1))

          if (zh_dist_fun.eq.1) then
*   ..apply charge state correction for the heat flux
            c_hw(0:nspecies-1,0:nspecies-1)=
     & c_hwA(ix,iy,0:nspecies-1,0:nspecies-1)
            do inucl = 0, nspecies-1
              c_hw(inucl,inucl)=c_hw(inucl,inucl)+zhcscorr*
     & c_hwB(inucl)-sum(c_hwA(ix,iy,0:nspecies-1,inucl))
            enddo

*   ..apply charge state correction for the r-moment
            c_rw(0:nspecies-1,0:nspecies-1)=
     & c_rwA(0:nspecies-1,0:nspecies-1)
            do inucl = 0, nspecies-1
              c_rw(inucl,inucl)=c_rw(inucl,inucl)+zhcscorr*
     & c_rwB(inucl)-sum(c_rwA(0:nspecies-1,inucl))
            enddo

*    ..moments of the distribution function 
            do inucl = 0, nspecies-1
              do iz= 1, znnucl_int(inucl)
                w_out(ix,iy,nucl2s(inucl,iz)) =
     & ua(ix,iy,nucl2s(inucl,iz))-avm_u(ix,iy)

                hTdp_out(ix,iy,nucl2s(inucl,iz)) =
     & sum(kappa_hTA(inucl,0:nspecies-1)*gtalc(ix,iy,0:nspecies-1))+
     & zhcscorr*kappa_hTB(inucl)*
     & (na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl)/
     & Ia(ix,iy,nucl2s(inucl,iz))*gt_Ac(ix,iy,nucl2s(inucl,iz))-
     & gtalc(ix,iy,inucl))

                hwdp_out(ix,iy,nucl2s(inucl,iz)) =
     & sum(c_hw(0:nspecies-1,inucl)*
     & (ua(ix,iy,nucl2s(inucl,iz))-av_ualpha(ix,iy,0:nspecies-1)))

                rTdpgamma_out(ix,iy,nucl2s(inucl,iz)) =
     & sum(kappa_rTA(inucl,0:nspecies-1)*gtalc(ix,iy,0:nspecies-1))+
     & zhcscorr*kappa_rTB(inucl)*
     & (na(ix,iy,nucl2s(inucl,iz))/nal(ix,iy,inucl)/
     & Ia(ix,iy,nucl2s(inucl,iz))*gt_Ac(ix,iy,nucl2s(inucl,iz))-
     & gtalc(ix,iy,inucl))

                rwdpgamma_out(ix,iy,nucl2s(inucl,iz)) =
     & sum(c_rw(0:nspecies-1,inucl)*
     & (ua(ix,iy,nucl2s(inucl,iz))-av_ualpha(ix,iy,0:nspecies-1)))
              enddo
            enddo
          endif

        enddo
      enddo
      
      call b2txcx(nx,ny,mode,vol,gs,hx,qz,hci0(-1,-1,0),
     & chci0(-1,-1,0))
      call b2txcx(nx,ny,mode,vol,gs,hx,qz,hci0_nofl(-1,-1,0),
     & chci0_nofl(-1,-1,0))
        do iy=-1,ny
          do ix=0,nx
            if (on_closed_surface(ix,iy)) then
              chci0_luc(ix,iy,0)=chci0(ix,iy,0)*f_luc_ki(ix,iy)
            else
              chci0_luc(ix,iy,0)=chci0(ix,iy,0)
            endif
          enddo
        enddo

      ncall = ncall+1
      call subend ()
      return


      end subroutine b2trzh
