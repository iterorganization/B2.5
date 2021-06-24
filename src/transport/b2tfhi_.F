*-----------------------------------------------------------------------
*.specification
csw 2013 added fna_he, cvsa, vadia into argument list
      subroutine b2tfhi_ (nx, ny, ns, facdrift,
     1  vol, hx, hy, qz, bb, ti, ni, na, fni, fna_he, ua,
     1  chci, cddi, chvi, chvimx, cvsa, vadia,                          !srv 17.05.99 27.12.01
     &  fchvispar, fchvisper, fchvisq, fchinert, fchanml,               !lk 30.05.11
     &  fhi, fhi_mdf, fhiPSch, floi, coni,                              !srv 13.10.06
     &  balance)                                                        !djm Jan2017
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_balance                                                 !djm Jan2017
     & , only : fhi_32, fhi_52, fhi_cond, fhi_dia, fhi_ecrb,
     &          fhi_strange, fhi_pschused, fhi_inert, fhi_vispar,
     &          fhi_visper, fhi_visq, fhi_anml, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     * ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     & vol(-1:nx,-1:ny), chci(-1:nx,-1:ny,0:1),
     & chvi(-1:nx,-1:ny,0:1),
     & cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     & cddi(-1:nx,-1:ny,0:1,0:ns-1), qz(-1:nx,-1:ny,0:1),               !srv 13.10.06
     & chvimx(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     & facdrift(-1:nx,-1:ny),
     & vadia(-1:nx,-1:ny,0:1,0:ns-1),                                   !srv 17.05.99 27.12.01 25.01.02
     & fchvispar(-1:nx,-1:ny,0:1), fchvisq(-1:nx,-1:ny,0:1),            !lk 30.05.11
     & fchinert(-1:nx,-1:ny,0:1), fchanml(-1:nx,-1:ny,0:1),             !lk 30.05.11
csw 30apr2014 added fna_he, ia
     . fchvisper(-1:nx,-1:ny,0:1),
     . fna_he(-1:nx,-1:ny,0:1,0:ns-1),
     . ua(-1:nx,-1:ny,0:ns-1)
      logical balance                                                   !djm Jan2017
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   fhi(-1:nx,-1:ny,0:1),
     &   fhi_mdf(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),            !srv 13.10.06
     1   floi(-1:nx,-1:ny,0:1), coni(-1:nx,-1:ny,0:1,0:2)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHI computes the all atom heat flux between neighboring cells,
*     and the coefficients that specify its linearization.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the atom heat flux, fhi, and the
*     coefficients floi and coni that linearize the response of the
*     atom heat flux to changes in the atom temperature.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ib
      integer, save :: iout = 0, iout_b2wdat = 0,                       !srv 16.09.02
     &                 iout_b2nph9 = 0, iout_b2npht = 0
      integer, save :: lim_flux = 0, mdf_fhi = 0                        !srv 04.02.03 13.10.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0             !srv 06.07.06
      integer, save :: conductive_limit = 1                             !srv 16.06.08 17.01.12
      integer, save :: b2tfhi_hybr2=0, b2_upwind=0, flo53=0
      integer, save :: b2kin=0, b2vis=0
      real (kind=R8), save :: PSch=1.0e0_R8                             !srv 17.01.07
      real (kind=R8), save :: xcur = 0.0_R8                             !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                             !srv 07.07.08
      real (kind=R8), save :: vis_per = 1.0_R8, vis_q = 1.0_R8          !srv 02.10.17
      real (kind=R8) ::
     *   flo, dif, floi0, flomx, coni0, t0, t1, dti, dti0, vdia_par,    !srv 17.05.99
     &   floi0_mdf, coni0_mdf, t0_mdf, t1_mdf                           !srv 13.10.06 17.06.09
cpax
      real (kind=R8) ::
     *   flidia(-1:nx,-1:ny,0:1), facdriftm                             !srv 25.01.02
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)          !srv 16.09.02
      real (kind=R8) :: upfac, hybfac, scur                             !djm Jan2017
      real (kind=R8) :: cutlo
      parameter (cutlo = 1.0e-60_R8)
      save ncall, vdia_par                                              !srv 17.05.99
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   hybr, hy1, upwind
      external subini, subend, xertst
      external b2xvsg, b2xvff
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 07.05.08
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
*   ..initialization
      data ncall/0/, vdia_par/1.0e0_R8/                                 !srv 17.05.99

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhi_')
*   ..set internal parameters on first call                             !srv 17.05.99
      if (ncall.eq.0) then                                              !srv
       call ipgetr ('b2tfhi_vdia_par', vdia_par)                        !srv 17.05.99
       call ipgeti ('b2tfhi_iout', iout)                                !srv 16.09.02
       call ipgeti ('b2tfhi_lim_flux', lim_flux)                        !srv 04.02.03
       call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)                 !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)                 !srv 06.07.06
       call ipgeti ('b2tfhi_mdf', mdf_fhi)                              !srv 13.10.06
       call ipgetr ('b2tfhi_PSch', PSch)                                !srv 17.01.07
       call ipgetr ('b2tfnb_xcur', xcur)                                !xpb 27.09.17
       call ipgetr ('b2tfnb_ycur', ycur)                                !srv 17.07.08
       call ipgetr ('b2tfnb_vis_per', vis_per)                          !srv 02.10.17
       call ipgetr ('b2tfnb_vis_q', vis_q)                              !srv 02.10.17
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)
       call ipgeti ('b2trcl_conductive_limit', conductive_limit)        !srv 16.06.08 17.01.12
       call ipgeti ('b2tfhi_upwind', b2_upwind)
       call ipgeti ('b2tfhi_hybr2', b2tfhi_hybr2)
       call ipgeti ('b2tfhi_flo53', flo53)
       call ipgeti ('b2tfhi_b2kin', b2kin)
       call ipgeti ('b2tfhi_b2vis', b2vis)
      endif                                                             !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fni
       call b2xvff (nx, ny, fni, 'fni')
*    ..test sign of ni, ti, chci
       call b2xvsg (n2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
*    ..test edge values of chci, chvi
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
      endif

      call b2tipsch (nx, ny, ns, facdrift, !vol, hx, hy,                 !srv 02.04.18
     1  na, ti, cddi, coni, fhiPSch)                                    !srv 13.10.06

      facdriftm = maxval(facdrift)
      flidia = 0.0_R8                                                   !srv 25.01.02 {
      if(facdriftm.ne.0.0e0_R8) then                                    !srv 18.01.07
       do is = 0, ns-1
        if ( .not.is_neutral(is) ) then
* ..compute flidiax
         do iy = -1, ny
         do ix = 0, nx
          if(leftix(ix,iy).ne.-2)then                                   !srv 26.11.02
           flidia(ix,iy,0) = flidia(ix,iy,0)
     &       +(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))             !srv 17.05.99 11.03.00 25.01.02 was 5/2
     &       /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))               !srv 17.05.99 11.03.00
     &       *(na(leftix(ix,iy),leftiy(ix,iy),is)                       !srv 17.05.99 11.03.00
!     &       *bb(leftix(ix,iy),leftiy(ix,iy),2)                        !srv 17.05.99 11.03.00 05.02.05
!     &       /bb(leftix(ix,iy),leftiy(ix,iy),3)                        !srv 17.05.99 11.03.00 05.02.05
     &        +na(ix,iy,is)
!     &        *bb(ix,iy,2)/bb(ix,iy,3)                                 !srv 05.02.05
     &                    )*0.5e0_R8                                    !srv 17.05.99 11.03.00
     &       *vadia(ix,iy,0,is)*2.5e0_R8*vdia_par                       !srv 17.05.99 11.03.00 21.04.00 13.09.02 18.01.07
!     &          +(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))         !srv 25.01.02 was 5/2
!     &           /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
!     &           *(bb(leftix(ix,iy),leftiy(ix,iy),2)
!     &            /bb(leftix(ix,iy),leftiy(ix,iy),3)
!     &             +bb(ix,iy,2)/bb(ix,iy,3))*0.5e0_R8
!     &     *(na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))/2.0e0_R8
!     &     *vadia(ix,iy,0,is)
          endif                                                         !srv 26.11.02
         enddo
         enddo
* ..compute flidiay
         do iy = 0, ny
         do ix = -1, nx
          flidia(ix,iy,1) = flidia(ix,iy,1)
     &        +(vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))        !srv 25.01.02 was 5/2
     &        /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &   *(na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))/2.0e0_R8
     &   *vadia(ix,iy,1,is)*2.5e0_R8*vdia_par                           !srv 13.09.02 18.01.07
         enddo
         enddo
        endif
       enddo                                                            !srv 25.01.02 }
      endif                                                             !srv 18.01.07

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &   .or. (iout_b2wdat.eq.4 .and. balance)) then                                    !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny)
     *   ,w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                             !srv 16.09.02 }

      fhi_mdf = 0.0_R8                                                  !srv 17.06.09
* ..compute fluxes and linearization
*   ..compute fhix, floix, conix
      do iy = -1, ny
       do ix = -1, nx                                                   !srv 26.11.02
        if(leftix(ix,iy).eq.-2)then                                     !srv 26.11.02
         floi(ix,iy,0)     = 0.0e0_R8
         coni(ix,iy,0,0)   = 0.0e0_R8                                   !srv 02.01.07
         fhi(ix,iy,0)      = 0.0e0_R8
         fhi_mdf(ix,iy,0)  = 0.0e0_R8                                   !srv 13.10.06
         if (balance.and.(balance_netcdf.ne.0)) then
           fhi_32(ix,iy,0) = 0.0_R8                               !djm Jan2017
           fhi_52(ix,iy,0) = 0.0_R8
           fhi_cond(ix,iy,0) = 0.0_R8
           fhi_dia(ix,iy,0) = 0.0_R8
           fhi_ecrb(ix,iy,0) = 0.0_R8
           fhi_strange(ix,iy,0) = 0.0_R8
           fhi_pschused(ix,iy,0) = 0.0_R8
           fhi_inert(ix,iy,0) = 0.0_R8
           fhi_vispar(ix,iy,0) = 0.0_R8
           fhi_visper(ix,iy,0) = 0.0_R8
           fhi_visq(ix,iy,0) = 0.0_R8
           fhi_anml(ix,iy,0) = 0.0_R8
         endif
        else
         dti = ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
         scur = 1.5e0_R8/qe*(fchinert(ix,iy,0)+fchvispar(ix,iy,0)
     &                      +fchanml(ix,iy,0)+fchvisper(ix,iy,0)
     &                      +fchvisq(ix,iy,0))*xcur
         floi0 = 1.5e0_R8*fni(ix,iy,0)+chvi(ix,iy,0)*
     &       (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8
     &           +flidia(ix,iy,0) ! *signmf                             !srv 25.01.02 18.01.07 11.01.13
     &           +scur
         floi0_mdf = 1.5e0_R8*fni(ix,iy,0)+chvi(ix,iy,0)*
     &    (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8
     &    +scur

         if(lim_flux.eq.1 .and. conductive_limit.eq.0) then             !srv 17.01.12

csw 17feb2012
          if(no_Ptncr_x_hi.eq.0) then
            coni0 = hybr(floi0,chci(ix,iy,0))-abs(floi0/2.e0)
          else
            coni0 = chci(ix,iy,0)-abs(floi0/2.e0)
          endif

cwdk 10.03.16 unified treatment of chvimx in b2tlhi.F. Test for b2tqca_model no longer required
          flomx = chvimx(ix,iy)*(ni(leftix(ix,iy),leftiy(ix,iy),0)
     &                          +ni(ix,iy,0))/2.d0
          if(flomx.gt.0.0_R8) then
            t0 =coni0*abs(dti)/(flomx*
     &              max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
          else
            t0 = 0.0_R8
          endif

          dti0 = dti/(1.0_R8+t0)
          t1 = t0/(1.0_R8+t0)

csw 13mar2013 modified hybrid SOLPS4 style (overwrite cone0)
          if(b2tfhi_hybr2 .ne. 0) then

csw 30apr2014 for TESTING only
           if(flo53 .ne. 0) then
            floi0=floi0*5.d0/3.d0
           endif

           coni0=dim(chci(ix,iy,0)/(1.d0+t0), abs(floi0)*
     .      max(0.0_R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .      max(cutlo,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
          endif


csw 27nov2012  added check upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
           fhi(ix,iy,0) = upwind(floi0,ti(leftix(ix,iy),leftiy(ix,iy))
     &                                  ,ti(ix,iy))
     &                    -coni0*dti
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.eq.0)) then !djm Jan2017
             if(floi0.ge.0.0_R8) then
               upfac=1.0_R8
             else
               upfac=0.0_R8
             endif
             fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_52(ix,iy,0) = 0.0_R8
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &        (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &        (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &        (1.0_R8-upfac)*ti(ix,iy))
             fhi_dia(ix,iy,0) = flidia(ix,iy,0)*
     &                          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*ti(ix,iy))
             fhi_cond(ix,iy,0) = -coni0*dti
             fhi_ecrb(ix,iy,0) = 0.0_R8
             fhi_pschused(ix,iy,0) = 0.0_R8
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
           endif
          else
           fhi(ix,iy,0) = floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                                  +ti(ix,iy))/2.0_R8
     &                    -coni0*dti0
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.eq.0)) then !djm Jan2017
             fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_52(ix,iy,0) = 0.0_R8
             fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &        (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_dia(ix,iy,0) = flidia(ix,iy,0)*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_cond(ix,iy,0) = -coni0*dti
             fhi_ecrb(ix,iy,0) = 0.0_R8
             fhi_pschused(ix,iy,0) = 0.0_R8
             fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
           endif
          endif


csw 2014 add kinetic term (ONLY FOR TESTING w/ SOLPS4 discretization)
          if(b2kin .ne. 0) then
           do is=0,ns-1
             fhi(ix,iy,0) = fhi(ix,iy,0) + fna_he(ix,iy,0,is) *
     .          0.5_R8 * am(is)*mp*
     .              (
     .           ua(leftix(ix,iy), leftiy(ix,iy),is)
     .              )**2
           enddo
          endif

csw 2014 add viscous term (assume: already flux limited cvsa, if any)
c     (ONLY FOR TESTING w/ SOLPS4 discretization)
          if(b2vis .ne. 0) then
           do is=0,ns-1
             fhi(ix,iy,0) = fhi(ix,iy,0) -
     .          0.5_R8*(cvsa(ix,iy,0,is)
     .          +cvsa(leftix(ix,iy),leftiy(ix,iy),0,is))/2.0_R8*
     .         (ua(ix,iy,is) + ua(leftix(ix,iy),leftiy(ix,iy),is))*
     .         (ua(ix,iy,is) - ua(leftix(ix,iy),leftiy(ix,iy),is))
           enddo
          endif



          if(mdf_fhi.eq.0) then

csw 30apr2014 for SOLPS4 style
           if(b2tfhi_hybr2 .ne.0) then
            floi(ix,iy,0) = floi0
            coni(ix,iy,0,0) = coni0
           else
            floi(ix,iy,0) = floi0-
     &      sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &      min(0.0e0_R8,t1)**2
            coni(ix,iy,0,0) = abs(floi0/2.e0)+coni0/(1+t0)              !srv 02.01.07
           endif

          else

csw 17feb2012
           if(no_Ptncr_x_hi.eq.0) then
            coni0_mdf =
     .        hybr(floi0_mdf,chci(ix,iy,0))-abs(floi0_mdf/2.0)          !srv 17.06.09 {
           else
            coni0_mdf = chci(ix,iy,0)-abs(floi0_mdf/2.0)
           endif

csw 13mar2013 SOLPS4 style (overwrite coni0_mdf)
           if(b2tfhi_hybr2 .ne. 0) then
            coni0_mdf=dim(chci(ix,iy,0), abs(floi0_mdf)*
     .       max(0._R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .       max(1.e-30_R8,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
           endif

           t0_mdf =coni0_mdf*abs(dti)/(flomx*
     &                  max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
           t1_mdf = t0_mdf/(1.0_R8+t0_mdf)                              !srv 17.06.09 {

csw 27nov2012 added check upwind (SOLPS4 style)
           if(b2_upwind .ne. 0) then
             fhi_mdf(ix,iy,0) = upwind(floi0_mdf,
     &           ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy))
     &       -coni0_mdf*dti0 - fhiPSch(ix,iy,0)                         !srv 17.06.09 }
             floi(ix,iy,0) = floi0_mdf-
     &        sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &        min(0.0e0_R8,t1_mdf)**2
             coni(ix,iy,0,0) = coni0_mdf/(1.0_R8+t0_mdf)                !srv 02.01.07
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhi.ne.0)) then !djm Jan2017
               if(floi0_mdf.ge.0.0_R8) then
                 upfac=1.0_R8
               else
                 upfac=0.0_R8
               endif
               fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &                           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*ti(ix,iy))
               fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &           (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &           (1.0_R8-upfac)*ti(ix,iy))
               fhi_52(ix,iy,0) = 0.0_R8
               fhi_cond(ix,iy,0) = -coni0_mdf*dti0
               fhi_dia(ix,iy,0) = 0.0_R8
               fhi_ecrb(ix,iy,0) = 0.0_R8
               fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
               fhi_inert(ix,iy,0) =
     &           1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_vispar(ix,iy,0) =
     &           1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_visper(ix,iy,0) =
     &           1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_visq(ix,iy,0) =
     &           1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_anml(ix,iy,0) =
     &           1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
             endif
           else
             fhi_mdf(ix,iy,0) = floi0_mdf*
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
     &       -abs(floi0_mdf/2.0)*dti-coni0_mdf*dti0 - fhiPSch(ix,iy,0)  !srv 17.06.09 }
             floi(ix,iy,0) = floi0_mdf-
     &        sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &        min(0.0e0_R8,t1_mdf)**2
             coni(ix,iy,0,0)=
     &        abs(floi0_mdf/2)+coni0_mdf/(1.0_R8+t0_mdf)                !srv 02.01.07
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhi.ne.0)) then !djm Jan2017
               !djm Jan2017 share the hybrid part according to relative contributions
               if (floi0_mdf.ge.0.0_R8) then
                 hybfac=1.0_R8
               else
                 hybfac=-1.0_R8
               endif
               fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &           hybfac*1.5e0_R8*fni(ix,iy,0)/2.0_R8*dti
               fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &          (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &          (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &           hybfac*chvi(ix,iy,0)*
     &          (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8/
     &           2.0_R8*dti
               fhi_52(ix,iy,0) = 0.0_R8
               fhi_cond(ix,iy,0) = -coni0_mdf*dti0
               fhi_dia(ix,iy,0) = 0.0_R8
               fhi_ecrb(ix,iy,0) = 0.0_R8
               fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
               fhi_inert(ix,iy,0) =
     &              1.5e0_R8/qe*fchinert(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_vispar(ix,iy,0) =
     &              1.5e0_R8/qe*fchvispar(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_visper(ix,iy,0) =
     &              1.5e0_R8/qe*fchvisper(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_visq(ix,iy,0) =
     &              1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_anml(ix,iy,0) =
     &              1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             endif
           endif
          endif

         else  !lim_flux  check                                         !srv 04.02.03 {

          if(no_Ptncr_x_hi.eq.0) then                                   !srv 06.07.06 {
           coni0_mdf = hybr(floi0_mdf,chci(ix,iy,0))
           coni0     = coni0_mdf
          else
           coni0     = chci(ix,iy,0)
           coni0_mdf = chci(ix,iy,0)
          endif                                                         !srv 06.07.06 }

csw 13mar2013 SOLPS4 style (overwrite coni0/coni0_mdf)
          if(b2tfhi_hybr2 .ne. 0) then
            coni0=dim(chci(ix,iy,0), abs(floi0)*
     .       max(0._R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .       max(1.e-30_R8,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
            coni0_mdf=dim(chci(ix,iy,0), abs(floi0_mdf)*
     .       max(0._R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .       max(1.e-30_R8,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
          endif

          dti0 = dti

          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                    .or. (iout_b2wdat.eq.4 .and. balance)) then   !srv 16.09.02 {
           if(mdf_fhi.eq.0) then                                        !srv 13.10.06 {
            w1(ix,iy)=floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                        +ti(ix,iy))/2.0_R8
            w2(ix,iy)= -coni0*dti0
            w6(ix,iy)=  coni0
           else
            w1(ix,iy)=floi0_mdf*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                        +ti(ix,iy))/2.0_R8
            w2(ix,iy)= -coni0_mdf*dti0
            w6(ix,iy)=  coni0_mdf
           endif
           w5(ix,iy) = flidia(ix,iy,0)*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
          endif


csw 27nov2012 added check upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
            fhi(ix,iy,0) = upwind(floi0,
     &          ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy))
     &           - coni0*dti0
            fhi_mdf(ix,iy,0) = upwind(floi0_mdf,                        !srv 13.10.06 {
     &           ti(leftix(ix,iy),leftiy(ix,iy)) ,ti(ix,iy))
     &       -coni0_mdf*dti0 - fhiPSch(ix,iy,0)                         !srv 13.10.06 }
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.eq.0)) then !djm Jan2017
              if(floi0.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &         (1.0_R8-upfac)*ti(ix,iy))
              fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &         (1.0_R8-upfac)*ti(ix,iy))
              fhi_dia(ix,iy,0) = flidia(ix,iy,0)*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &         (1.0_R8-upfac)*ti(ix,iy))
              fhi_52(ix,iy,0) = 0.0_R8
              fhi_cond(ix,iy,0) = -coni0*dti0
              fhi_ecrb(ix,iy,0) = 0.0_R8
              fhi_pschused(ix,iy,0) = 0.0_R8
              fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.ne.0)) then !djm Jan2017
              if(floi0_mdf.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &         (1.0_R8-upfac)*ti(ix,iy))
              fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &         (1.0_R8-upfac)*ti(ix,iy))
              fhi_dia(ix,iy,0) = 0.0_R8
              fhi_52(ix,iy,0) = 0.0_R8
              fhi_cond(ix,iy,0) = -coni0_mdf*dti0
              fhi_ecrb(ix,iy,0) = 0.0_R8
              fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
              fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &                     xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
              fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
            endif
          else
            fhi(ix,iy,0) = floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                            +ti(ix,iy))/2.0_R8
     &      - coni0*dti0
            fhi_mdf(ix,iy,0) = floi0_mdf*                               !srv 13.10.06 {
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
     &       -coni0_mdf*dti0 - fhiPSch(ix,iy,0)                         !srv 13.10.06 }
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.eq.0)) then !djm Jan2017
              fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_52(ix,iy,0) = 0.0_R8
              fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &               (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_dia(ix,iy,0) = flidia(ix,iy,0)*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_cond(ix,iy,0) = -coni0*dti0
              fhi_ecrb(ix,iy,0) = 0.0_R8
              fhi_pschused(ix,iy,0) = 0.0_R8
              fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.ne.0)) then !djm Jan2017
              fhi_32(ix,iy,0) = 1.5e0_R8*fni(ix,iy,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_52(ix,iy,0) = 0.0_R8
              fhi_strange(ix,iy,0) = chvi(ix,iy,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_dia(ix,iy,0) = 0.0_R8
              fhi_cond(ix,iy,0) = -coni0_mdf*dti0
              fhi_ecrb(ix,iy,0) = 0.0_R8
              fhi_pschused(ix,iy,0) = -fhiPSch(ix,iy,0)
              fhi_inert(ix,iy,0) = 1.5e0_R8/qe*fchinert(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_vispar(ix,iy,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visper(ix,iy,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visq(ix,iy,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_anml(ix,iy,0) = 1.5e0_R8/qe*fchanml(ix,iy,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
            endif
          endif

          if(mdf_fhi.eq.0) then                                         !srv 13.10.06 {
           floi(ix,iy,0)   = floi0
           coni(ix,iy,0,0) = coni0                                      !srv 02.01.07
          else
           floi(ix,iy,0)   = floi0_mdf
           coni(ix,iy,0,0) = coni0_mdf                                  !srv 02.01.07
          endif                                                         !srv 13.10.06 }
         endif                                                          !srv 04.02.03 }
        endif                                                           !srv 26.11.02
       enddo
      enddo

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &   .or. (iout_b2wdat.eq.4 .and. balance)) then                    !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhi__qi_32x')
       call my_out(70,nx,ny,w2,'b2tfhi__qi_ki_gTx')
!       w3=w1+w2
!       call my_out(70,nx,ny,w3,'qix3')
       call my_out(70,nx,ny,w5,'b2tfhi__qidiax')
       call my_out(70,nx,ny,w6,'b2tfhi__chcix_Pat')
       call my_out(70,nx,ny,fhiPSch(-1,-1,0),'b2tfhi__fhiPSchx')
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                             !srv 16.09.02 }

*   ..compute fhiy, floiy, coniy
      do ix = -1, nx
       floi(ix,-1,1)    = 0
       coni(ix,-1,1,0)  = 0                                             !srv 02.01.07
       fhi(ix,-1,1)     = 0.0e0_R8
       fhi_mdf(ix,-1,1) = 0.0e0_R8                                      !srv 13.10.06
       if (balance.and.(balance_netcdf.ne.0)) then
         fhi_32(ix,-1,1) = 0.0_R8                                       !djm Jan2017
         fhi_52(ix,-1,1) = 0.0_R8
         fhi_cond(ix,-1,1) = 0.0_R8
         fhi_dia(ix,-1,1) = 0.0_R8
         fhi_ecrb(ix,-1,1) = 0.0_R8
         fhi_strange(ix,-1,1) = 0.0_R8
         fhi_pschused(ix,-1,1) = 0.0_R8
         fhi_inert(ix,-1,1) = 0.0_R8
         fhi_vispar(ix,-1,1) = 0.0_R8
         fhi_visper(ix,-1,1) = 0.0_R8
         fhi_visq(ix,-1,1) = 0.0_R8
         fhi_anml(ix,-1,1) = 0.0_R8
       endif
      enddo
      do iy = 0, ny
       do ix = -1, nx
        scur = 1.5e0_R8/qe*(fchinert(ix,iy,1)+fchvispar(ix,iy,1)        !lk 30.05.11
     &                     +fchanml(ix,iy,1)+vis_per*fchvisper(ix,iy,1) !srv 02.10.17
     &                     +vis_q*fchvisq(ix,iy,1))*ycur                !lk 30.05.11 !srv 02.10.17
        floi0 = 1.5e0_R8*fni(ix,iy,1)+                                  !srv 11.05.99
     &    chvi(ix,iy,1)*(ni(bottomix(ix,iy),bottomiy(ix,iy),0)
     &                   +ni(ix,iy,0))/2.0_R8
     &        +flidia(ix,iy,1)                                          !srv 25.01.02 18.01.07
     &    +scur
        floi0_mdf = 1.5e0_R8*fni(ix,iy,1)+                              !srv 11.05.99
     &    chvi(ix,iy,1)*(ni(bottomix(ix,iy),bottomiy(ix,iy),0)
     &                   +ni(ix,iy,0))/2.0_R8
     &   +scur                                                          !lk 30.05.11
        if(no_Ptncr_y_hi.eq.0) then                                     !srv 06.07.06 {
         coni0_mdf = hybr(floi0_mdf,chci(ix,iy,1))
         coni0     = coni0_mdf
        else
         coni0     = chci(ix,iy,1)
         coni0_mdf = chci(ix,iy,1)
        endif                                                           !srv 06.07.06 }

csw 13mar2013 SOLPS4 hybrid (overwrite coni0/coni0_mdf)
        if(b2tfhi_hybr2 .ne. 0) then

csw 30apr2014 only for TESTING
         if(flo53 .ne .0) then
          floi0=floi0*5.0d0/3.0d0
         endif

         coni0=dim(chci(ix,iy,1), abs(floi0)*
     .   max(0._R8,min(ti(ix,iy),ti(bottomix(ix,iy), bottomiy(ix,iy))))/
     .   max(1.e-30_R8,ti(ix,iy)+ti(bottomix(ix,iy), bottomiy(ix,iy))))
         coni0_mdf=dim(chci(ix,iy,1), abs(floi0_mdf)*
     .   max(0._R8,min(ti(ix,iy),ti(bottomix(ix,iy), bottomiy(ix,iy))))/
     .   max(1.e-30_R8,ti(ix,iy)+ti(bottomix(ix,iy), bottomiy(ix,iy))))
        endif


        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &          .or. (iout_b2wdat.eq.4 .and. balance)) then             !srv 16.09.02 {
         if(mdf_fhi.eq.0) then                                          !srv 13.10.06 {
          w1(ix,iy)=1.5e0_R8*fni(ix,iy,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w2(ix,iy)=flidia(ix,iy,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w3(ix,iy)=-coni0*
     1    (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=0.0e0_R8
          w6(ix,iy)=  coni0
         else
          w1(ix,iy)=1.5e0_R8*fni(ix,iy,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w2(ix,iy)=0.0e0_R8
          w3(ix,iy)=-coni0_mdf*
     1     (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=0.0e0_R8
          w6(ix,iy)=  coni0_mdf
         endif
         w5(ix,iy) = flidia(ix,iy,1)*
     &                 (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2

        endif                                                           !srv 16.09.02 }

csw 13mar2013 added check upwind (SOLPS4)
        if(b2_upwind .ne. 0) then
         fhi(ix,iy,1) = upwind(floi0,
     1    ti(bottomix(ix,iy),bottomiy(ix,iy)),ti(ix,iy))
     &   -coni0*(ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
         fhi_mdf(ix,iy,1) = upwind(floi0_mdf,                           !srv 13.10.06 {
     1    ti(bottomix(ix,iy),bottomiy(ix,iy)) , ti(ix,iy))
     &   -coni0_mdf*(ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &    -fhiPSch(ix,iy,1)                                             !srv 13.10.06 }
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.eq.0)) then !djm Jan2017
           if(floi0.ge.0.0_R8) then
             upfac=1.0_R8
           else
             upfac=0.0_R8
           endif
           fhi_32(ix,iy,1) = 1.5e0_R8*fni(ix,iy,1)*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_strange(ix,iy,1) = chvi(ix,iy,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &      (1.0_R8-upfac)*ti(ix,iy))
           fhi_dia(ix,iy,1) = flidia(ix,iy,1)*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_inert(ix,iy,1) = 1.5e0_R8/qe*fchinert(ix,iy,1)*ycur*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_vispar(ix,iy,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1)*
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_visper(ix,iy,1) =
     &                      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1)*     !srv 02.10.17
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                      (1.0_R8-upfac)*ti(ix,iy))
           fhi_visq(ix,iy,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1)*      !srv 02.10.17
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_anml(ix,iy,1) = 1.5e0_R8/qe*fchanml(ix,iy,1)*ycur*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_52(ix,iy,1) = 0.0_R8
           fhi_cond(ix,iy,1) = -coni0*
     &                   (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1) = 0.0_R8
           fhi_pschused(ix,iy,1) = 0.0_R8
         endif
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.ne.0)) then !djm Jan2017
           if(floi0_mdf.ge.0.0_R8) then
             upfac=1.0_R8
           else
             upfac=0.0_R8
           endif
           fhi_32(ix,iy,1) = 1.5e0_R8*fni(ix,iy,1)*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_strange(ix,iy,1) = chvi(ix,iy,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &      (1.0_R8-upfac)*ti(ix,iy))
           fhi_inert(ix,iy,1) = 1.5e0_R8/qe*fchinert(ix,iy,1)*ycur*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_vispar(ix,iy,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1)*
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_visper(ix,iy,1) =
     &                      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1)*     !srv 02.10.17
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                      (1.0_R8-upfac)*ti(ix,iy))
           fhi_visq(ix,iy,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1)*      !srv 02.10.17
     &                  ycur*(upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_anml(ix,iy,1) = 1.5e0_R8/qe*fchanml(ix,iy,1)*ycur*
     &                       (upfac*ti(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*ti(ix,iy))
           fhi_dia(ix,iy,1) = 0.0_R8
           fhi_52(ix,iy,1) = 0.0_R8
           fhi_cond(ix,iy,1) = -coni0_mdf*
     &                   (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1) = 0.0_R8
           fhi_pschused(ix,iy,1) = -fhiPSch(ix,iy,1)
         endif
        else
         fhi(ix,iy,1) = floi0*
     &         (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8 -
     &         coni0*
     &         (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
         fhi_mdf(ix,iy,1) = floi0_mdf*                                  !srv 13.10.06 {
     &         (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8 -
     &         coni0_mdf*
     &         (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &    -fhiPSch(ix,iy,1)                                             !srv 13.10.06 }
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.eq.0)) then !djm Jan2017
           fhi_32(ix,iy,1) = 1.5e0_R8*fni(ix,iy,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_strange(ix,iy,1) = chvi(ix,iy,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_dia(ix,iy,1) = flidia(ix,iy,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_inert(ix,iy,1) = 1.5e0_R8/qe*fchinert(ix,iy,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_vispar(ix,iy,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1)*
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visper(ix,iy,1) =
     &      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1)*ycur*                !srv 02.10.17
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visq(ix,iy,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1)*      !srv 02.10.17
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_anml(ix,iy,1) = 1.5e0_R8/qe*fchanml(ix,iy,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_52(ix,iy,1) = 0.0_R8
           fhi_cond(ix,iy,1) = -coni0*
     &      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1) = 0.0_R8
           fhi_pschused(ix,iy,1) = 0.0_R8
         endif
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.ne.0)) then !djm Jan2017
           fhi_32(ix,iy,1) = 1.5e0_R8*fni(ix,iy,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_strange(ix,iy,1) = chvi(ix,iy,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_dia(ix,iy,1) = 0.0_R8
           fhi_inert(ix,iy,1) = 1.5e0_R8/qe*fchinert(ix,iy,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_vispar(ix,iy,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1)*
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visper(ix,iy,1) =
     &      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1)*ycur*                !srv 02.10.17
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visq(ix,iy,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1)*      !srv 02.10.17
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_anml(ix,iy,1) = 1.5e0_R8/qe*fchanml(ix,iy,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_52(ix,iy,1) = 0.0_R8
           fhi_cond(ix,iy,1) = -coni0_mdf*
     &      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1) = 0.0_R8
           fhi_pschused(ix,iy,1) = -fhiPSch(ix,iy,1)
         endif
        endif


csw 2014 add kinetic term (ONLY FOR TESTING w/ SOLPS4 discretization)
          if(b2kin .ne. 0) then
           do is=0,ns-1
             fhi(ix,iy,1) = fhi(ix,iy,1) + fna_he(ix,iy,1,is) *
     .          0.5_R8 * am(is)*mp*
     .              (
     .           ua(ix,iy,is)
     .              )**2
           enddo
          endif

csw 2014 add viscous term   (ONLY FOR TESTING w/ SOLPS4 discretization)
          if(b2vis .ne. 0) then
           do is=0,ns-1
             fhi(ix,iy,1) = fhi(ix,iy,1) -
     .          0.5_R8*cvsa(ix,iy,1,is)*
     .         (ua(ix,iy,is) - ua(bottomix(ix,iy),bottomiy(ix,iy),is))*
     .         (ua(ix,iy,is) + ua(bottomix(ix,iy),bottomiy(ix,iy),is))
           enddo
          endif

         if(mdf_fhi.eq.0) then                                          !srv 13.10.06 {
          floi(ix,iy,1)   = floi0
          coni(ix,iy,1,0) = coni0                                       !srv 02.01.07
         else
          floi(ix,iy,1)   = floi0_mdf
          coni(ix,iy,1,0) = coni0_mdf                                   !srv 02.01.07
         endif                                                          !srv 13.10.06 }

       enddo
      enddo

csw 18feb2014 remove fluxes between guard cells
        do ix=-1,nx
         do iy=-1,ny
           if(iffboundx(ix,iy).ne.0) then
             fhi(ix,iy,0)=0.0_R8
             fhi_mdf(ix,iy,0)=0.0_R8
             fhiPSch(ix,iy,0)=0.0_R8
             floi(ix,iy,0)=0.0_R8
             coni(ix,iy,0,:)=0.0_R8
             if (balance.and.(balance_netcdf.ne.0)) then
               fhi_32(ix,iy,0) = 0.0_R8 !djm Jan2017
               fhi_52(ix,iy,0) = 0.0_R8
               fhi_cond(ix,iy,0) = 0.0_R8
               fhi_dia(ix,iy,0) = 0.0_R8
               fhi_ecrb(ix,iy,0) = 0.0_R8
               fhi_strange(ix,iy,0) = 0.0_R8
               fhi_pschused(ix,iy,0) = 0.0_R8
               fhi_inert(ix,iy,0) = 0.0_R8
               fhi_vispar(ix,iy,0) = 0.0_R8
               fhi_visper(ix,iy,0) = 0.0_R8
               fhi_visq(ix,iy,0) = 0.0_R8
               fhi_anml(ix,iy,0) = 0.0_R8
             endif
           endif
           if(iffboundy(ix,iy).ne.0) then
             fhi(ix,iy,1)=0.0_R8
             fhi_mdf(ix,iy,1)=0.0_R8
             fhiPSch(ix,iy,1)=0.0_R8
             floi(ix,iy,1)=0.0_R8
             coni(ix,iy,1,:)=0.0_R8
             if (balance.and.(balance_netcdf.ne.0)) then
               fhi_32(ix,iy,1) = 0.0_R8 !djm Jan2017
               fhi_52(ix,iy,1) = 0.0_R8
               fhi_cond(ix,iy,1) = 0.0_R8
               fhi_dia(ix,iy,1) = 0.0_R8
               fhi_ecrb(ix,iy,1) = 0.0_R8
               fhi_strange(ix,iy,1) = 0.0_R8
               fhi_pschused(ix,iy,1) = 0.0_R8
               fhi_inert(ix,iy,1) = 0.0_R8
               fhi_vispar(ix,iy,1) = 0.0_R8
               fhi_visper(ix,iy,1) = 0.0_R8
               fhi_visq(ix,iy,1) = 0.0_R8
               fhi_anml(ix,iy,1) = 0.0_R8
             endif
           endif
         enddo
        enddo


      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &           .or. (iout_b2wdat.eq.4 .and. balance)) then            !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhi__qi_32y')
!       call my_out(70,nx,ny,w2,'qiy2')
       call my_out(70,nx,ny,w3,'b2tfhi__qi_ki_gTy')
!       call my_out(70,nx,ny,w4,'qiy4')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'qiy5')
       call my_out(70,nx,ny,w5,'b2tfhi__qidiay')
       call my_out(70,nx,ny,w6,'b2tfhi__chciy_Pat')                     !srv 16.09.02 }
       call my_out(70,nx,ny,fhiPSch(-1,-1,1),'b2tfhi__fhiPSchy')
       call my_out(70,nx,ny,fhi(-1,-1,0),'b2tfhi__fhix')
       call my_out(70,nx,ny,fhi(-1,-1,1),'b2tfhi__fhiy')
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &           .or. (iout_b2wdat.eq.4 .and. balance)) then
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
       DO IX=NXISO1(IB), NXISO2(IB)
        DO IY=NYISO1(IB), NYISO2(IB)
         fhi(ix,iy,0)   =0.0e0_R8
         floi(ix,iy,0)  =0.0e0_R8
         coni(ix,iy,0,:)=0.0e0_R8                                       !srv 02.01.07
         if (balance.and.(balance_netcdf.ne.0)) then
           fhi_32(ix,iy,0) = 0.0_R8 !djm Jan2017
           fhi_52(ix,iy,0) = 0.0_R8
           fhi_cond(ix,iy,0) = 0.0_R8
           fhi_dia(ix,iy,0) = 0.0_R8
           fhi_ecrb(ix,iy,0) = 0.0_R8
           fhi_strange(ix,iy,0) = 0.0_R8
           fhi_pschused(ix,iy,0) = 0.0_R8
           fhi_inert(ix,iy,0) = 0.0_R8
           fhi_vispar(ix,iy,0) = 0.0_R8
           fhi_visper(ix,iy,0) = 0.0_R8
           fhi_visq(ix,iy,0) = 0.0_R8
           fhi_anml(ix,iy,0) = 0.0_R8
         endif
        enddo
       end do
       DO IX=NXISO1(IB), NXISO2(IB)
        DO IY=NYISO1(IB), NYISO2(IB)
          fhi(ix,iy,1)   =0.0e0_R8
          floi(ix,iy,1)  =0.0e0_R8
          coni(ix,iy,1,:)=0.0e0_R8                                      !srv 02.01.07
          if (balance.and.(balance_netcdf.ne.0)) then
            fhi_32(ix,iy,1) = 0.0_R8 !djm Jan2017
            fhi_52(ix,iy,1) = 0.0_R8
            fhi_cond(ix,iy,1) = 0.0_R8
            fhi_dia(ix,iy,1) = 0.0_R8
            fhi_ecrb(ix,iy,1) = 0.0_R8
            fhi_strange(ix,iy,1) = 0.0_R8
            fhi_pschused(ix,iy,1) = 0.0_R8
            fhi_inert(ix,iy,1) = 0.0_R8
            fhi_vispar(ix,iy,1) = 0.0_R8
            fhi_visper(ix,iy,1) = 0.0_R8
            fhi_visq(ix,iy,1) = 0.0_R8
            fhi_anml(ix,iy,1) = 0.0_R8
          endif
        enddo
       enddo
      enddo
      if(iout.ne.0) then                                                !srv 18.06.08 {
       call my_out(70,nx,ny,floi(-1,-1,0),'b2tfhi__floix')
       call my_out(70,nx,ny,floi(-1,-1,1),'b2tfhi__floiy')
       call my_out(70,nx,ny,coni(-1,-1,0,0),'b2tfhi__conix0')
       call my_out(70,nx,ny,coni(-1,-1,0,1),'b2tfhi__conix1')
       call my_out(70,nx,ny,coni(-1,-1,0,2),'b2tfhi__conix2')
       call my_out(70,nx,ny,coni(-1,-1,1,0),'b2tfhi__coniy0')
       call my_out(70,nx,ny,coni(-1,-1,1,1),'b2tfhi__coniy1')
       call my_out(70,nx,ny,coni(-1,-1,1,2),'b2tfhi__coniy2')
       call my_out(70,nx,ny,fhi(-1,-1,0),'b2tfhi__fhix')
       call my_out(70,nx,ny,fhi(-1,-1,1),'b2tfhi__fhiy')
       call my_out(70,nx,ny,fhi_mdf(-1,-1,0),'b2tfhi__fhi_mdfx')
       call my_out(70,nx,ny,fhi_mdf(-1,-1,1),'b2tfhi__fhi_mdfy')
       call my_out(70,nx,ny,fhiPSch(-1,-1,0),'b2tfhi__fhiPSchx')
       call my_out(70,nx,ny,fhiPSch(-1,-1,1),'b2tfhi__fhiPSchy')
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'b2tfhi__fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'b2tfhi__fchvispary')
       call my_out(70,nx,ny,fchvisq(-1,-1,0),'b2tfhi__fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1),'b2tfhi__fchvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tfhi__fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tfhi__fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0),'b2tfhi__fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1),'b2tfhi__fchinerty')
       call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfhi__fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfhi__fchanmly')
       call my_out(70,nx,ny,fni(-1,-1,0),'b2tfhi__fnix')
       call my_out(70,nx,ny,fni(-1,-1,1),'b2tfhi__fniy')
       call my_out(70,nx,ny,ni(-1,-1,0),'b2tfhi__ni0')
       call my_out(70,nx,ny,ni(-1,-1,1),'b2tfhi__ni1')
      endif                                                             !srv 18.06.08 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Reconsider factor 2.5d0 versus 1.5d0 in transport and linearization.

*-----------------------------------------------------------------------
*.end b2tfhi_

      end subroutine b2tfhi_
