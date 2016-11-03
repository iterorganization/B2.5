*-----------------------------------------------------------------------
*.specification

      subroutine b2tfch (nx, ny, ns, guard_flows, 
     1  facdrift, fac_vis, BoRiS,
     1  vol, gs, hx, hy, qz, bb, te, po, ne,                    	 !xpb
     1  ti, na, fna, uadia, vedia, wadia, vaecrb, veecrb, ua,
     2  fne_32, fne_52, csig, csigin, chce, chve,               	 !xpb
     &  calf, chvemx, cvsa, fch, conc, fchdia,
     &  fch_p, fchvispar, fchvisq, fchin, fch_32, fch_52,                !srv 25.11.08
     1  avisper, bvisper, cvisper, fchvisper, csig_cl, calf_cl)
      use b2mod_types
      use b2mod_rates
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport                             	 !xpb
      use b2mod_geo_corner                                      	 !xpb
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows                           	 !xpb
      real (kind=R8) ::
     * te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     & vol(-1:nx,-1:ny), fna(-1:nx,-1:ny,0:1,0:ns-1),
     & gs(-1:nx,-1:ny,0:2), qz(-1:nx,-1:ny,0:1),                	 !xpb
     & uadia(-1:nx,-1:ny,0:1,0:ns-1), 
     & ua(-1:nx,-1:ny,0:ns-1),
     & wadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1), 	 !xpb
     & fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     1 csig(-1:nx,-1:ny,0:1), 
     & csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),                   	 !xpb
     & chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     & calf(-1:nx,-1:ny,0:1), chvemx(-1:nx,-1:ny),
     & cvsa(-1:nx,-1:ny,0:1,0:ns-1),na(-1:nx,-1:ny,0:ns-1),
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 bb(-1:nx,-1:ny,0:3), ti(-1:nx,-1:ny), 
     2 facdrift(-1:nx,-1:ny), BoRiS,
     3 fac_vis(-1:nx,-1:ny),
     4 csig_cl(-1:nx,-1:ny,0:1), calf_cl(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  fch(-1:nx,-1:ny,0:1), 
     &  conc(-1:nx,-1:ny,0:1),
     *  vedia(-1:nx,-1:ny,0:1), veecrb(-1:nx,-1:ny,0:1),
     2  fchdia(-1:nx,-1:ny,0:1), 
     3  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),              !xpb srv 10.06.08
     4  fch_32(-1:nx,-1:ny,0:1), fch_52(-1:nx,-1:ny,0:1),
     5  avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),
     6  cvisper(-1:nx,-1:ny,0:1), fchvisper(-1:nx,-1:ny,0:1),
     7  fchvisq(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1)                 !srv 25.11.08
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFCH computes the electric current between neighboring cells,
*     and the coefficients that specify the linearised response to
*     changes in electric potential.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, is, is0, 
     *   no_current, conduction_only, no_hybr, pot_eq
      integer, save :: iout = 0                                         !srv 18.06.08
      integer, save :: drift_style = 1
      real (kind=R8) ::
     *   ehx, ehy, flo, dif, floe0, flomx, cone0, t0,
     &   tex, tex0, ehx_po, ehy_po,
     1   flo_vedia, prefactor, floe0_nd, cone0_nd, t0_nd, tex0_nd
      real (kind=R8) ::                                                 !xpb
     * fch_ion_neutral, fch_anomalous, vis_par,                         !xpb
     * dvpar(-1:nx,-1:ny), vis_per,                                     !xpb
     * va(-1:nx,-1:ny), dva(-1:nx,-1:ny),
     * fnadspl(-1:nx,-1:ny,0:1,0:ns-1)
      character*1 charnc
#ifdef CLASSICAL_CURRENT
      real (kind=R8) ::
     *   gampvd(-1:nx,-1:ny), gamrvd(-1:nx,-1:ny),
     *   gampvdc(-1:nx,-1:ny)
#endif
      save ncall, fch_ion_neutral, fch_anomalous, vis_par, pot_eq,
     *   no_current, conduction_only, no_hybr, vis_per
*   ..procedures
      intrinsic max, abs, sqrt, log, mod                        	!xpb
      real (kind=R8) ::
     *   hybr, hy1
      external subini, subend, xertst, ipgetr, ipgeti
      external b2xvsg, b2xvff, b2xvfx, intcor
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
*   ..initialisation
      data ncall/0/, fch_ion_neutral/0.0_R8/, fch_anomalous/1.0_R8/,  	!xpb
     * vis_par/0.0_R8/, no_current/0/, conduction_only/0/, 
     1 vis_per/0.0_R8/, no_hybr/0/, pot_eq/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfch')
      if (ncall.eq.0) then                                      	!xpb
       call ipgetr ('b2tfhe_neutral', fch_ion_neutral)          	!xpb
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)          	!xpb
       if (fch_anomalous.eq.1.0_R8)
     >  call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
       call ipgetr ('b2tfhe_vis_par', vis_par)                  	!xpb
       call ipgetr ('b2tfhe_vis_per', vis_per)                  	!xpb
       call xertst (fch_ion_neutral.ne.0.0_R8 .or. vis_par.ne.0.0_R8 	!xpb
     &           .or. fch_anomalous.ne.0.0_R8 
     &           .or. maxval(facdrift).ne.0.0_R8
     &           .or. vis_per.ne.0.0_R8,
     & 'faulty argument fch_ion_neutral/fch_anomalous/facdrift/vis_p.r') !xpb
       call xertst (0.0_R8.le.vis_par .and. 1.0_R8.ge.vis_par,          !xpb
     &  'faulty argument vis_par')                                      !xpb
       call xertst (0.0_R8.le.vis_per .and. 1.0_R8.ge.vis_per,          !xpb
     &  'faulty argument vis_per')                                      !xpb
       call xertst (0.0_R8.le.fch_ion_neutral .and.                     !xpb
     &  1.0_R8.ge.fch_ion_neutral, 'faulty argument fch_ion_neutral')   !xpb
       call ipgeti ('b2tfhe_no_current', no_current)                    !xpb
       call ipgeti ('b2news_poteq', pot_eq)                             !srv 27.02.13
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_conduction_only', conduction_only)          !xpb
       call ipgeti ('b2tfhe_no_hybr', no_hybr)                          !xpb
       call ipgeti ('b2tfch_iout', iout)                                !srv 18.06.08
       call ipgeti ('b2tfnb_drift_style', drift_style)
      endif                                                             !xpb
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fne
       call b2xvff (nx, ny, fne_32, 'fne_32')
       call b2xvff (nx, ny, fne_52, 'fne_52')
*    ..test sign of ne, te, csig, csig_an, calf_an, csigin, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
       call b2xvsg (n2*2, csig_an, 1, 'csig_an', '.ge.')
       call b2xvsg (n2*2, calf_an, 1, 'calf_an', '.ge.')
       call b2xvsg (n2*2*ns*ns, csigin, 1, 'csigin', '.ge.')
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
*    ..test edge values of csig, calf, csig_an, calf_an, csigin, chce, chve, chvemx
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, calf, 'calf')
       call b2xvff (nx, ny, csig_an, 'csig_an')
       call b2xvff (nx, ny, calf_an, 'calf_an')
       do is = 0, ns-1
         do is0 = 0, ns-1
           call b2xvff (nx, ny, csigin(-1,-1,0,is,is0), 'csigin')
         enddo
       enddo
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvfx (nx, ny, chvemx, 'chvemx')
      endif

c Drift terms
      fchdia=0.0_R8
      if(maxval(facdrift).ne.0.0_R8) then
*       X-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT             !srv
        do iy = -1, ny                                                   !xpb
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
cxpb  Must check for the case of the island X-point
              if(topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then   	 !xpb
cxpb  Standard case
                fchdia(ix,iy,0)=fchdia(ix,iy,0)                          !srv
     &           +(ne(ix,iy)*(te(ix,iy)+ti(ix,iy))*bb(ix,iy,2)+
     &             ne(leftix(ix,iy),leftiy(ix,iy))                       !srv
     &           *(te(leftix(ix,iy),leftiy(ix,iy))                       !srv
     &            +ti(leftix(ix,iy),leftiy(ix,iy)))                      !srv
     &            *bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8             !xpb
     &           *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))          !srv
     &           /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))        !srv
     &           *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))  !srv
              else                                                       !xpb
                fchdia(ix,iy,0)=fchdia(ix,iy,0)                          !srv
     &           +(ne(ix,iy)*(te(ix,iy)+ti(ix,iy))*bb(ix,iy,2)+
     &             ne(leftix(ix,iy),leftiy(ix,iy))                       !srv
     &           *(te(leftix(ix,iy),leftiy(ix,iy))                       !srv
     &            +ti(leftix(ix,iy),leftiy(ix,iy)))                      !srv
     &            *bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8             !xpb
     &           *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))          !srv
     &           /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))        !srv
     &           *((vol(topix(ix,iy),topiy(ix,iy))/wbbl(ix,iy,3)**2+     !xpb
     &              vol(ix,iy)/wbbr(topix(ix,iy),topiy(ix,iy),3)**2)/    !xpb
     &             (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))           !xpb
     &           -OnedBsqc(ix,iy))                                       !srv
              endif                                                      !xpb
            endif
          enddo
        enddo

*       Y-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT             !srv
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
cxpb  must check for the diverted island case where B is not periodic
              if(periodic_bc.ne.1 .or. rightix(ix,iy).eq.ix+1 .or.
     &           region(ix,iy,0).eq.5 .or. nncut.eq.0) then
cxpb  standard case
                fchdia(ix,iy,1)=fchdia(ix,iy,1)                          !srv
     &           -(ne(ix,iy)*(te(ix,iy)+ti(ix,iy))*bb(ix,iy,2)+
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))                   !srv
     &           *(te(bottomix(ix,iy),bottomiy(ix,iy))                   !srv
     &            +ti(bottomix(ix,iy),bottomiy(ix,iy)))                  !srv
     &            *bb(bottomix(ix,iy),bottomiy(ix,iy),2))/2.0_R8         !srv
     &           *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))      !srv
     &           /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))    !srv
     &     *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))    !srv
              else
cxpb  geometrical cut to the right with jump in B
      if(ncall.eq.0)
     > write(*,*) 'b2tfch: special treatment for cell ix,iy = ',ix,iy
                fchdia(ix,iy,1)=fchdia(ix,iy,1)                          !srv
     &           -(ne(ix,iy)*(te(ix,iy)+ti(ix,iy))*bb(ix,iy,2)+
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))                   !srv
     &           *(te(bottomix(ix,iy),bottomiy(ix,iy))                   !srv
     &            +ti(bottomix(ix,iy),bottomiy(ix,iy)))                  !srv
     &            *bb(bottomix(ix,iy),bottomiy(ix,iy),2))/2.0_R8         !srv
     &           *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))      !srv
     &           /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))    !srv
     &           *(1.0_R8/wbbv(ix,iy,3)**2-OnedBsqc(ix,iy))              !srv
              endif
            endif
          enddo
        enddo

#ifdef CLASSICAL_CURRENT
crfs inertia terms
        gampvd=0.0_R8
        gamrvd=0.0_R8
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            do iy = -1, ny
              do ix = -1, nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1.and.
     &             bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1.and.
     &             region(ix,iy,0).ne.0)
     >            gampvd(ix,iy)=gampvd(ix,iy)+am(is)*mp*
     1            0.5_R8*(fna(ix,iy,0,is)*uadia(ix,iy,0,is)+
     2            fna(rightix(ix,iy),rightiy(ix,iy),0,is)
     2            *uadia(rightix(ix,iy),rightiy(ix,iy),0,is))
              enddo
            enddo
            do ix = -1, nx
              do iy = -1, ny
                if(rightix(ix,iy).ne.nx+1 .and. topiy(ix,iy).ne.ny+1
     &       .and. leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2
     &       .and. region(ix,iy,0).ne.0) then
cxpb  Must check for case of the island X-point
                  if(topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then     	 !xpb
cxpb  Standard case
                    gamrvd(ix,iy)=gamrvd(ix,iy)+am(is)*mp*
     1               0.5_R8*(fna(ix,iy,1,is)+
     2               fna(topix(ix,iy),topiy(ix,iy),1,is))*
     1               0.5_R8*(uadia(ix,iy,0,is)+
     2               uadia(rightix(ix,iy),rightiy(ix,iy),0,is))
                  else                                                   !xpb
cxpb  Island X-point cells
                    gamrvd(ix,iy)=gamrvd(ix,iy)+am(is)*mp*
     1               fna(ix,iy,1,is)*                                    !xpb
     2               0.5_R8*(uadia(ix,iy,0,is)+
     2               uadia(rightix(ix,iy),rightiy(ix,iy),0,is))
                  endif
                endif                                           	 !xpb
              enddo
            enddo
          endif
        enddo
        call intcor(nx,ny,vol,gampvd,gampvdc)
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
              fchdia(ix,iy,1) = fchdia(ix,iy,1)+
     1         0.5_R8/wbbv(ix,iy,3)*
     2         ((gampvdc(rightix(ix,iy),rightiy(ix,iy))
     3          -gampvdc(ix,iy))/
     6         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
     2         (gamrvd(ix,iy)
     3          -gamrvd(bottomix(ix,iy),bottomiy(ix,iy)))/
     8         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
            endif
          enddo
        enddo
       	! /* CLASSICAL_CURRENT */
#endif

        fchdia(:,:,0)=fchdia(:,:,0)*facdrift(:,:)                        !xpb
        fchdia(:,:,1)=fchdia(:,:,1)*facdrift(:,:)                        !xpb
      endif                                                              !xpb

      fchvispar=0.0_R8                                                   !xpb
      dvpar=0.0_R8                                                       !srv 24.06.99 
      if (vis_par.ne.0.0_R8) then                                        !srv
        do is = 0, ns-1                                                  !srv 24.06.99
          if (.not.is_neutral(is)) then                                  !srv 24.06.99
            do iy= -1, ny                                       !xpb     !srv 18.02.00
              do ix= -1, nx                                              !srv
               if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &            bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &            region(ix,iy,0).ne.0) then
                 dvpar(ix,iy)=dvpar(ix,iy)                               !srv 24.06.99
     &             +wbbl(ix,iy,3)**2/sqrt(wbbl(ix,iy,3))/wbbl(ix,iy,0)   !srv
     &             *cvsa(ix,iy,0,is)*                                    !xpb
     &             ((ua(ix,iy,is)+bb(ix,iy,3)/bb(ix,iy,0)*               !srv 11.01.13
     &              (wadia(ix,iy,0,is)+vaecrb(ix,iy,0,is)))
     &              *sqrt(bb(ix,iy,3))             
     &             -(ua(leftix(ix,iy),leftiy(ix,iy),is)       
     &              +bb(leftix(ix,iy),leftiy(ix,iy),3)/
     &               bb(leftix(ix,iy),leftiy(ix,iy),0)*
     &              (wadia(leftix(ix,iy),leftiy(ix,iy),0,is)+
     &               vaecrb(leftix(ix,iy),leftiy(ix,iy),0,is)))*
     &               sqrt(wbbl(ix,iy,3)))
               endif
              enddo                                                      !srv
            enddo                                                        !srv
          endif                                                          !srv 24.06.99
        enddo                                                            !srv 24.06.99
        do iy = -1, ny
          do ix = -1, nx
            if(bottomiy(ix,iy).eq.-2)
     >       dvpar(ix,iy)=dvpar(topix(ix,iy),topiy(ix,iy))       !xpb    !srv
            if(topiy(ix,iy).eq.ny+1)
     >       dvpar(ix,iy)=dvpar(bottomix(ix,iy),bottomiy(ix,iy)) !xpb    !srv
          enddo
        enddo
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2) 
     >       dvpar(ix,iy)=dvpar(rightix(ix,iy),rightiy(ix,iy))   !xpb    !srv
            if(rightix(ix,iy).eq.nx+1)
     >       dvpar(ix,iy)=dvpar(leftix(ix,iy),leftiy(ix,iy))     !xpb    !srv
          enddo
        enddo

*       X-COMPONENT OF VISCOSITY PARALLEL CURRENT                        !srv
        do iy = -1, ny                                                   !srv
          do ix = -1, nx                                !xpb     !srv 01.07.99
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
cxpb  Must check for the island X-point
              if(topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then   	 !xpb
cxpb  Standard case
                fchvispar(ix,iy,0)= -vis_par*fac_vis(ix,iy)*             !srv 25.04.00
     &           (dvpar(ix,iy)+dvpar(leftix(ix,iy),leftiy(ix,iy)))       !srv
     &           /(hy1(leftix(ix,iy),leftiy(ix,iy))+hy1(ix,iy))
     &           *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))  !srv
     &           /3.0_R8                                         !srv 18.02.00
              else                                                       !xpb
cxpb  Island X-point
                fchvispar(ix,iy,0)= -vis_par*fac_vis(ix,iy)*             !srv 25.04.00
     &           (dvpar(ix,iy)+dvpar(leftix(ix,iy),leftiy(ix,iy)))       !srv
     &           /(hy1(leftix(ix,iy),leftiy(ix,iy))+hy1(ix,iy))
     &           *((vol(topix(ix,iy),topiy(ix,iy))/wbbl(ix,iy,3)**2+     !xpb
     &              vol(ix,iy)/wbbr(topix(ix,iy),topiy(ix,iy),3)**2)/    !xpb
     &             (vol(ix,iy)+vol(topix(ix,iy),topiy(ix,iy)))           !xpb
     &           -OnedBsqc(ix,iy))                                       !xpb
     &           /3.0_R8                                          !srv 18.02.00
              endif
            endif                                               	 !xpb
          enddo         
        enddo 

*       Y-COMPONENT OF VISCOSITY PARALLEL CURRENT                        !srv
        do iy = -1, ny                                  !xpb     !srv 01.07.99
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
cxpb  must check for the diverted island case where B is not periodic
              if(periodic_bc.ne.1 .or. rightix(ix,iy).eq.ix+1 .or.
     &           region(ix,iy,0).eq.5 .or. nncut.eq.0) then
cxpb  standard case
                fchvispar(ix,iy,1)= vis_par*fac_vis(ix,iy)*      !srv 25.04.00
     &           abs(wbbv(ix,iy,2)/wbbv(ix,iy,3))*
     &           (dvpar(ix,iy)+dvpar(bottomix(ix,iy),bottomiy(ix,iy)))   !srv
     &           /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv
     &           /3.0_R8                                         !srv 18.02.00
              else
cxpb  geometrical cut to the right with jump in B
                fchvispar(ix,iy,1)= vis_par*fac_vis(ix,iy)*      !srv 25.04.00
     &           abs(wbbv(ix,iy,2)/wbbv(ix,iy,3))*
     &           (dvpar(ix,iy)+dvpar(bottomix(ix,iy),bottomiy(ix,iy)))   !srv
     &           /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
     &           *(1.0_R8/wbbv(ix,iy,3)**2-OnedBsqc(ix,iy))      !xpb 
     &           /3.0_R8                                         !srv 18.02.00
              endif
            endif
          enddo
        enddo 
      endif                                                              !srv

c perpendicular viscosity driven current

      fchvisper=0.0_R8                                          	 !xpb
      avisper=0.0_R8
      cvisper=0.0_R8
      if (vis_per.ne.0.0_R8) then                                        !srv
        do iy= -1, ny
          do ix= -1, nx
            if(bottomiy(ix,iy).ne.-2) then
              cvisper(ix,iy,1)=2.0_R8/wbbv(ix,iy,3)/
     1         (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
              avisper(ix,iy,1)=cvisper(ix,iy,1)
            endif
          enddo
        enddo
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            bvisper=0.0_R8
            va=0.0_R8
            dva=0.0_R8
            do iy= -1, ny
              do ix= -1, nx
                bvisper(ix,iy,1)=
     1           rza(ix,iy,is)*na(ix,iy,is)/ne(ix,iy)* 
     2           0.25_R8*vsa0(ix,iy,is)*vol(ix,iy)/hy1(ix,iy)**2
              enddo
            enddo
            do iy= -1, ny
              do ix= -1, nx
                if(bottomiy(ix,iy).ne.-2) then
                  va(ix,iy)=-cvisper(ix,iy,1)*
     1             (po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     &             -facdrift(ix,iy)*cvisper(ix,iy,1)*
     &               (na(ix,iy,is)*ti(ix,iy)-
     1                na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     2                  ti(bottomix(ix,iy),bottomiy(ix,iy)))*
     3                  4.0_R8/
     .                  ((rza(ix,iy,is)+
     +                    rza(bottomix(ix,iy),bottomiy(ix,iy),is))*qe)/
     4                  (na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     5                   na(ix,iy,is))
                endif
              enddo
            enddo
            do iy= -1, ny
              do ix= -1, nx
                if(topiy(ix,iy).ne.ny+1) then
                  dva(ix,iy)=-bvisper(ix,iy,1)*
     1             (va(topix(ix,iy),topiy(ix,iy))-va(ix,iy))
                endif
              enddo
            enddo
            do iy= -1, ny
              do ix= -1, nx
                if(bottomiy(ix,iy).ne.-2.and.
     2             topiy(ix,iy).ne.ny+1.and.
     5             leftix(ix,iy).ne.-2.and.
     6             rightix(ix,iy).ne.nx+1) then
                  if(
     1             bottomiy(bottomix(ix,iy),bottomiy(ix,iy)).ne.-2.and.
     3             bottomiy(bottomix(ix,iy),bottomiy(ix,iy)).ne.iy.and.
     4             topiy(topix(ix,iy),topiy(ix,iy)).ne.iy)
     >             fchvisper(ix,iy,1)=fchvisper(ix,iy,1)-
     1              avisper(ix,iy,1)*vis_per*fac_vis(ix,iy)*
     2              (dva(ix,iy)-dva(bottomix(ix,iy),bottomiy(ix,iy)))
                endif
              enddo
            enddo
          endif
        enddo
c       write (charnc,'(i1)') ncall
c       call my_out(70,nx,ny,avisper(:,:,1),'avispery'//charnc)
c       call my_out(70,nx,ny,bvisper(:,:,1),'bvispery'//charnc)
c       call my_out(70,nx,ny,cvisper(:,:,1),'cvispery'//charnc)
c       call my_out(70,nx,ny,fchvisper(:,:,1),'fchvispery'//charnc)
c       call xerrab (' ')
      endif

* ..calculate the parallel viscosity current connected with the heat flux 
      fchvisq = 0.0_R8                     ! Temporary                   !srv 25.11.08

      fchin=0.0_R8                                                       !xpb
      if (fch_ion_neutral .ne. 0.0_R8) then                              !srv 01.07.99
        do is0 = 0, ns-1
          if (is_neutral(is0)) then
*       X-COMPONENT OF ION-NEUTRAL CURRENT                               !srv
            call intcor(nx,ny,vol,fna(-1,-1,1,is0),fnadspl(-1,-1,1,is0)) !srv 28.03.99
            do is = 0, ns-1
              if (.not.is_neutral(is)) then
                do iy = -1, ny
                  do ix = -1, nx                                         !srv
                if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1
     &        .and. bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1
     &        .and. region(ix,iy,0).ne.0) then
             fchin(ix,iy,0)=fchin(ix,iy,0)+csigin(ix,iy,0,is,is0)*((
     &       -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))                !srv
     &       -(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))+                !srv
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))*0.5_R8*        !srv
     &        (log(na(ix,iy,is))-
     &         log(na(leftix(ix,iy),leftiy(ix,iy),is)))                  !srv
     &        )/(qe*rza(ix,iy,is)))                                      !srv
     1       *((bb(ix,iy,2)/bb(ix,iy,3))**2+                             !srv 28.03.99
     2         (bb(leftix(ix,iy),leftiy(ix,iy),2)/                       !srv
     3          bb(leftix(ix,iy),leftiy(ix,iy),3))**2)*0.5_R8 !  *signmf !srv 11.01.13
     &        +(bb(ix,iy,2)*gs(ix,iy,2)/vol(ix,iy)                       !xpb
     &                  /na(ix,iy,is0)+                                  !srv 28.03.99
     &          bb(leftix(ix,iy),leftiy(ix,iy),2)*
     &          gs(leftix(ix,iy),leftiy(ix,iy),2)                        !srv 28.03.99
     &         /vol(leftix(ix,iy),leftiy(ix,iy))/                        !srv 28.03.99
     &          na(leftix(ix,iy),leftiy(ix,iy),is0)
     &         )*0.5_R8*fnadspl(ix,iy,1,is0))*
     1	       fch_ion_neutral*fac_vis(ix,iy)
                endif
                  enddo                                                  !srv
                enddo                                                    !srv
              endif
            enddo
*       Y-COMPONENT OF ION-NEUTRAL CURRENT                               !srv
            call intcor(nx,ny,vol,fna(-1,-1,0,is0),fnadspl(-1,-1,0,is0)) !srv 28.03.99
            do is = 0, ns-1
              if (.not.is_neutral(is)) then
                do iy = -1, ny                                           !srv
                  do ix= -1, nx                                          !srv
                if (bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1
     &        .and. leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1
     &        .and. region(ix,iy,0).ne.0)
     >           fchin(ix,iy,1)=fchin(ix,iy,1)+csigin(ix,iy,1,is,is0)*(  !srv
     &       -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))            !srv
     &       -(                                                          !srv
     &         ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy))+            !srv 28.03.99
     &        (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))*0.5_R8*    !srv
     &        (log(na(ix,iy,is))-                                        !srv
     &         log(na(bottomix(ix,iy),bottomiy(ix,iy),is))))/
     &        (qe*rza(ix,iy,is))                                         !srv
     &       -(bb(ix,iy,3)*gs(ix,iy,2)/vol(ix,iy)                        !xpb
     &        /na(ix,iy,is0)+                                            !srv Vpern=Vxn
     &         bb(bottomix(ix,iy),bottomiy(ix,iy),3)*                    !srv 28.03.99
     &         gs(bottomix(ix,iy),bottomiy(ix,iy),2)                     !srv 28.03.99
     &        /vol(bottomix(ix,iy),bottomiy(ix,iy))/                     !srv 28.03.99
     &         na(bottomix(ix,iy),bottomiy(ix,iy),is0)                   !srv
     &        )*0.5_R8*fnadspl(ix,iy,0,is0))*                            !srv
     1	       fch_ion_neutral*fac_vis(ix,iy)
                  enddo                                                  !srv
                enddo                                                    !srv
              endif
            enddo
          endif
        enddo
      endif
        
* ..compute fluxes and linearisation
*  ..do not forget equivalent code in B2TFHE !!!!!!!!!!!!!!!!!!!!!
*   ..compute fchx, concx
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
            conc(ix,iy,0) = 0.0_R8
            fch(ix,iy,0) = 0.0_R8
            fch_p(ix,iy,0) = 0.0_R8                                  !srv 10.06.08
            fch_32(ix,iy,0) = 0.0_R8
            fch_52(ix,iy,0) = 0.0_R8
          else
!xpb        ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
!xpb     &         (ne(ix,iy)*te(ix,iy)
!xpb     &         -ne(leftix(ix,iy),leftiy(ix,iy))*
!xpb     &         te(leftix(ix,iy),leftiy(ix,iy)))*
!xpb     &         (2/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))
            if (no_current.eq.1) then
              ehx = 0.0_R8
              ehx_po = 0.0_R8
            else
              ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+     	 !xpb
     &         (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))+      	 !xpb
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))*0.5_R8*    	 !xpb
     &         (log(ne(ix,iy))-log(ne(leftix(ix,iy),leftiy(ix,iy))))  	 !xpb
     &         )/qe                                             	 !xpb
              ehx_po = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))   	 !xpb
            endif
            tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
cxpb Changed the solution to (3/2+BoRiS) for benchmarking !!!!!!!!!!!!!!!!!!!
            if (conduction_only.eq.1) then
              floe0_nd = 0.0_R8
              flo_vedia=0
            else
              floe0_nd = (1.5_R8+BoRiS)*fne_32(ix,iy,0)+
     1         2.5_R8*fne_52(ix,iy,0)+
     &         chve(ix,iy,0)*
     1         (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))/2.0_R8
              prefactor=
     1         (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))    	 !xpb
     &         /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))     	 !xpb
     &         *(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))/2.0_R8
              if (drift_style.eq.0) then
                flo_vedia=0.5_R8*prefactor*(vedia(ix,iy,0) + 
     1            vedia(leftix(ix,iy),leftiy(ix,iy),1))
              else
                flo_vedia=prefactor*vedia(ix,iy,0)
              endif
            endif
            floe0=floe0_nd+2.5_R8*flo_vedia
            if (no_hybr.eq.0) then
              cone0 = hybr(floe0,chce(ix,iy,0))-abs(floe0/2.0_R8)
              cone0_nd = hybr(floe0_nd,chce(ix,iy,0))-abs(floe0_nd/2)
            else
              cone0 = chce(ix,iy,0)
              cone0_nd = chce(ix,iy,0)
            endif
cwdk 10.03.16 unified treatment of chvemx in b2tlhe.F. Test for b2tqce_model no longer required
            flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &        ne(ix,iy))/2.0_R8
            if (flomx.gt.0.0_R8) then
              t0 = cone0*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
              t0_nd = cone0_nd*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            else
              t0 = 0.0_R8
              t0_nd = 0.0_R8
            endif
            tex0 = tex/(1.0_R8+t0)
            tex0_nd = tex/(1.0_R8+t0_nd)
            fch(ix,iy,0) = csig_cl(ix,iy,0)*ehx+csig_an(ix,iy,0)*ehx_po
     &       -calf(ix,iy,0)*tex0 	                                 !xpb
     &       +fchdia(ix,iy,0)+fchin(ix,iy,0)+fchvispar(ix,iy,0) 	 !xpb
     &       +fchvisper(ix,iy,0)
            fch_32(ix,iy,0) = csig_cl(ix,iy,0)*ehx
     1       +fchin(ix,iy,0)+fchvispar(ix,iy,0)
     &       +fchvisper(ix,iy,0)
            fch_52(ix,iy,0) = csig_an(ix,iy,0)*ehx_po
     &                       -calf(ix,iy,0)*tex0_nd
            if (no_current.ne.1) then
              fch_p(ix,iy,0) = csig_cl(ix,iy,0)*ehx
     &                        -calf(ix,iy,0)*tex0                        !srv 10.06.08
            endif
            conc(ix,iy,0) = csig(ix,iy,0)
            if(fch_ion_neutral.ne.0.0_R8) then                  	 !xpb
              do is = 0, ns-1
                if (.not.is_neutral(is)) then
                  do is0 = 0, ns-1
                    if (is_neutral(is0)) conc(ix,iy,0) = conc(ix,iy,0) + 
     &               csigin(ix,iy,0,is,is0)
                  enddo
                endif
              enddo
            endif                                               	 !xpb
          endif
        enddo
      enddo
*   ..compute fchy, concy
      do iy = -1, ny
        do ix = -1, nx
          if (bottomiy(ix,iy).eq.-2) then
            conc(ix,iy,1) = 0.0_R8
            fch(ix,iy,1) = 0.0_R8
            fch_p(ix,iy,1) = 0.0_R8                                  !srv 10.06.08
            fch_32(ix,iy,1) = 0.0_R8
            fch_52(ix,iy,1) = 0.0_R8
          else
            if (no_current.eq.1) then
              ehy = 0.0_R8
              ehy_po = 0.0_R8
            else
              ehy = -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &         (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy))+
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))*0.5_R8*
     &         (log(ne(ix,iy))-log(ne(bottomix(ix,iy),bottomiy(ix,iy))))
     &         )/qe    
              ehy_po = -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))  !xpb
            endif
            fch(ix,iy,1)= csig_cl(ix,iy,1)*ehy+csig_an(ix,iy,1)*ehy_po
     .       -calf(ix,iy,1)*
     .       (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))             !xpb
     .       +fchdia(ix,iy,1)+fchin(ix,iy,1)+fchvispar(ix,iy,1)          !xpb
     &       +fchvisper(ix,iy,1)
            fch_p(ix,iy,1)= 0.0_R8                                       !srv 10.06.08
            fch_32(ix,iy,1)= csig_cl(ix,iy,1)*ehy
     .       +fchin(ix,iy,1)+fchvispar(ix,iy,1)                          !xpb
     &       +fchvisper(ix,iy,1)
            fch_52(ix,iy,1)= csig_an(ix,iy,1)*ehy_po-calf(ix,iy,1)*
     .       (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy))) 
            conc(ix,iy,1) = csig(ix,iy,1)
            if(fch_ion_neutral .ne. 0.0_R8) then                	 !xpb
              do is = 0, ns-1
                if (.not.is_neutral(is)) then
                  do is0 = 0, ns-1
                    if (is_neutral(is0)) conc(ix,iy,1) = conc(ix,iy,1) 
     .               + csigin(ix,iy,1,is,is0)                     	 !xpb
                  enddo
                endif
              enddo
            endif
          endif                                                 	 !xpb
        enddo
      enddo

      if(guard_flows.eq.0) then                                 	 !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              conc(ix,iy,1)=0.0_R8
              fch(ix,iy,1)=0.0_R8
              fch_32(ix,iy,1)=0.0_R8
              fch_52(ix,iy,1)=0.0_R8
              fchdia(ix,iy,1)=0.0_R8
              fchin(ix,iy,1)=0.0_R8
              fchvispar(ix,iy,1)=0.0_R8
              fchvisper(ix,iy,1)=0.0_R8
              fchvisq(ix,iy,1)=0.0_R8
            endif
          enddo        
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(mod(region(ix,-1,0),4).eq.1) then
            conc(ix,-1,0) = 0.0_R8
            fch(ix,-1,0) = 0.0_R8
            fch_32(ix,-1,0) = 0.0_R8
            fch_52(ix,-1,0) = 0.0_R8
            fch_p(ix,-1,0) = 0.0_R8                                  !srv 10.06.08
            fchdia(ix,-1,0) = 0.0_R8
            fchin(ix,-1,0) = 0.0_R8
            fchvispar(ix,-1,0) = 0.0_R8
            fchvisper(ix,-1,0) = 0.0_R8
            fchvisq(ix,-1,0) = 0.0_R8
          endif
        enddo
      endif                                                     	 !xpb

cxpb turn off current terms if requested
      if (no_current.eq.1) then
        conc=0.0_R8
        fch=0.0_R8
        fch_p=0.0_R8
        fch_32=0.0_R8
        fch_52=0.0_R8
        fchdia=0.0_R8
        fchin=0.0_R8
        fchvispar=0.0_R8
        fchvisper=0.0_R8
        fchvisq=0.0_R8
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            conc(ix,iy,0)=0.0_R8
            fch(ix,iy,0)=0.0_R8
            fch_32(ix,iy,0)=0.0_R8
            fch_52(ix,iy,0)=0.0_R8
            fchdia(ix,iy,0)=0.0_R8
            fch_p(ix,iy,:) =0.0_R8
            fchvispar(ix,iy,0)=0.0_R8
            fchvisper(ix,iy,0)=0.0_R8
            fchvisq(ix,iy,0)=0.0_R8
            fchin(ix,iy,0)=0.0_R8
            conc(ix,iy,1)=0.0_R8
            fch(ix,iy,1)=0.0_R8
            fch_32(ix,iy,1)=0.0_R8
            fch_52(ix,iy,1)=0.0_R8
            fchdia(ix,iy,1)=0.0_R8
            fchin(ix,iy,1)=0.0_R8
            fchvispar(ix,iy,1)=0.0_R8
            fchvisper(ix,iy,1)=0.0_R8
            fchvisq(ix,iy,1)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                             !srv 18.06.08 {
       call my_out(70,nx,ny,fch(-1,-1,0),'b2tfch_fchx')
       call my_out(70,nx,ny,fch(-1,-1,1),'b2tfch_fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0),'b2tfch_chinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'b2tfch_chiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'b2tfch_chvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'b2tfch_chvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0),'b2tfch_chvisqx')
!       call my_out(70,nx,ny,fchvisq(-1,-1,1),'b2tfch_chvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tfch_chvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tfch_chvispery')
!       call my_out(70,nx,ny,fchinert(-1,-1,0),'b2tfch_chinertx')
!       call my_out(70,nx,ny,fchinert(-1,-1,1),'b2tfch_chinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0),'b2tfch_chdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1),'b2tfch_chdiay')
!       call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfch_chanmlx')
!       call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfch_chanmly')
       call my_out(70,nx,ny,conc(-1,-1,0),'b2tfch_concx')
       call my_out(70,nx,ny,conc(-1,-1,1),'b2tfch_concy')
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfch_ch_px')
       call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfch_ch_py')
       call my_out(70,nx,ny,avisper(-1,-1,0),'b2tfch_avisperx')
       call my_out(70,nx,ny,avisper(-1,-1,1),'b2tfch_avispery')
       call my_out(70,nx,ny,bvisper(-1,-1,0),'b2tfch_bvisperx')
       call my_out(70,nx,ny,bvisper(-1,-1,1),'b2tfch_bvispery')
       call my_out(70,nx,ny,cvisper(-1,-1,0),'b2tfch_cvisperx')
       call my_out(70,nx,ny,cvisper(-1,-1,1),'b2tfch_cvispery')
      endif                                                          !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

cxpb  This version of the routine uses a (1.5+BoRiS) factor for the
cxpb  poloidal convective heat flux instead of the usual 3/2 or 5/2 
cxpb  so as to solve the same equation as in BoRiS for benchmark purposes. 
cxpb  The default is BoRiS = 0.0.
cxpb  For comparison, set BoRiS = 1.0.
cxpb  Routine rewritten to explicitly separate 3/2 and 5/2 terms.
cxpb  The contributions from fne_32 and fne_52 are NO DRIFTS
cxpb  contributions, i.e. diamagnetic and ExB terms are treated
cxpb  separately. By the same token, the computed fch_32 and fch_52
cxpb  terms are also NO DRIFTS currents, i.e. without contribution
cxpb  from diamagnetic or ExB terms.

*-----------------------------------------------------------------------
*.end b2tfch

      end subroutine b2tfch
