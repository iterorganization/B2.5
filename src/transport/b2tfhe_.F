*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhe_ (nx, ny, ns, ismain, facdrift,                  !srv 23.11.10
     1 vol, gs, hx, hy, qz, bb, te, po, ne, ne2, lnlam, pz,              !srv 11.07.99 18.01.02 01.07.09
     1 ti, ni, na, fna, vadia, ua, fne, csig, csigin, chce, chve, chci,  !srv 17.05.99 11.07.99
     & calf, chvemx, cvsa, cdna, cdpa, cddi, fch, fhe, fhe_mdf, fhePSch, !srv 27.01.00 13.10.06
     & conc, avisper, bvisper, cvisper, floe, cone, fchdia, fchin,       !srv 27.01.00 30.09.02
     & fch_p, fchvispar, fchvisper, fchvisq, fchinert, fchanml)          !srv 07.07.99 11.07.99 22.11.99 24.03.07
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport                                      !srv
      use b2mod_geo_corner                                         !srv
      use b2mod_rates                                                 !srv 22.06.08
      use b2mod_b2cmfs                                                !srv 22.06.08
      use b2mod_b2cmpa                                                !srv 22.06.08
      use b2mod_work, only : iffboundx, iffboundy
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     & vol(-1:nx,-1:ny), fna(-1:nx,-1:ny,0:1,0:ns-1), 
     & gs(-1:nx,-1:ny,0:2), ne2(-1:nx,-1:ny), pz(-1:nx,-1:ny),     !srv 11.07.99 18.01.02
     & vadia(-1:nx,-1:ny,0:1,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),      !srv 17.05.99 25.01.02
     & fne(-1:nx,-1:ny,0:1), csig(-1:nx,-1:ny,0:1),
     & csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1), chci(-1:nx,-1:ny,0:1), !srv 11.07.99 24.06.08
     & chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     & calf(-1:nx,-1:ny,0:1), chvemx(-1:nx,-1:ny),
     & cvsa(-1:nx,-1:ny,0:1,0:ns-1),na(-1:nx,-1:ny,0:ns-1),
     & cdna(-1:nx,-1:ny,0:1,0:ns-1), cdpa(-1:nx,-1:ny,0:1,0:ns-1), !srv 27.01.00
     & cddi(-1:nx,-1:ny,0:1,0:ns-1), lnlam(-1:nx,-1:ny),           !srv 13.10.06 01.07.09
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     1 bb(-1:nx,-1:ny,0:3), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny),
     2 facdrift(-1:nx,-1:ny)
      real (kind=R8) :: cutlo
      parameter (cutlo = 1.0e-60_R8)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   fch(-1:nx,-1:ny,0:1), fhe(-1:nx,-1:ny,0:1),
     &   fhe_mdf(-1:nx,-1:ny,0:1), fhePSch(-1:nx,-1:ny,0:1),      !srv 13.10.06
     &   fch_p(-1:nx,-1:ny,0:1),fchvispar(-1:nx,-1:ny,0:1),       !srv
     &   fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),      !srv 22.11.99 30.09.02
     &   conc(-1:nx,-1:ny,0:1), fchvisq(-1:nx,-1:ny,0:1),         !srv 11.07.99
     1   floe(-1:nx,-1:ny,0:1), cone(-1:nx,-1:ny,0:1,0:2),        !srv 02.01.07
     2   fchdia(-1:nx,-1:ny,0:1), fchinert(-1:nx,-1:ny,0:1),      !srv 07.07.99
     3   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),      !srv 27.01.00
     $   cvisper(-1:nx,-1:ny,0:1), fchanml(-1:nx,-1:ny,0:1)       !srv 24.03.07
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFHE_ computes the electric current and the electron heat flux
*     between neighboring cells, and the coefficients that specify the
*     linearized response to changes in electric potential and electron
*     temperature.
*
*
*
*  8. further comments
*
*     In the evaluation of the thermo-electric term in the electron
*     heat flux, the expression harm(te(),te()) cannot be replaced
*     by the arithmetic average (this is an empirical result).
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, n2, ix, iy, ib, is, ixl, iyl, ixb, iyb, is0
      integer, save :: iout = 0,                                      !srv 16.09.02
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0
      integer, save :: lim_flux = 1, mdf_fhe = 0                      !srv 04.02.03 13.10.06
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0           !srv 06.07.06
      integer, save :: pot_eq = 1                                     !srv 26.07.06
      integer, save :: mode_ehx = 0                                   !srv 21.10.09
      integer, save :: qstyle = 0                                     !lk 23.05.06 05.07.07 {
      integer, save :: conductive limit = 0                           !srv 16.06.08 17.01.12
      integer mintopcut                                               !lk  
      integer, save :: no_current=0                          !srv 27.02.13
csw 10oct2012
      integer, save :: b2_upwind =0     
      integer, save :: b2tfhe_hybr2=0, flo53=0

      character(2) :: charns, charns0
      real (kind=R8), save :: alfTeEh = 1.0_R8                     !srv 28.04.08 30.05.08
      real (kind=R8), save :: fch_pTe = 0.0_R8                        !srv 30.12.09
      real (kind=R8), save :: cqip1 = 0.5e0_R8                        !lk   
      real (kind=R8) :: tt0(-1:ny), tt1(-1:ny), bbaverage2(-1:ny)     !lk
     &                  , conn(-1:ny), tt2(-1:ny), tt3(-1:ny)
     &                  , eps1(-1:ny)
      real (kind=R8) :: qip0a(-1:nx,-1:ny), qip0b(-1:nx,-1:ny)
     &                  , tauia1(-1:nx,-1:ny), alpha1(-1:nx,-1:ny)
     &                  , nu1(-1:nx,-1:ny)                                  !lk 23.05.06 05.07.07 }
      real (kind=R8) :: csig_stoch(-1:nx,-1:ny),                     !lk 31.01.08 {
     &                  fchstoch(-1:nx,-1:ny)       
      real (kind=R8), save :: fch_stochastic = 0.0e0_R8               !lk 31.01.08 }
      real (kind=R8), save :: PSch=1.0e0_R8                           !srv 17.01.07
      real (kind=R8) ::
     *   ehx, ehy, flo, dif, floe0, flomx, cone0, alfe0, t0, t1,
     &  tex, tex0, dia_cur, vis_par, vis_q, vdia_par, !srv 17.05.99 11.07.99
     &  vis_per, floe0_mdf, cone0_mdf, alfe0_mdf, t0_mdf, tex0_mdf    !srv 22.11.99 13.10.06 17.06.09
      real (kind=R8) ::
     *    fch_ion_neutral, fch_anomalous       !srv 30.09.02
     *   ,fnadspl(-1:nx,-1:ny,0:1,0:ns-1), tauia(-1:nx,-1:ny,0:ns-1)  !srv 28.03.99 11.07.99
     *   ,fch_inert                                                   !srv 07.07.99 
     *   ,va(-1:nx,-1:ny,0:1)                                         !srv 22.11.99
     *   ,dva(-1:nx,-1:ny), fledia(-1:nx,-1:ny,0:1)                   !srv 22.11.99 25.01.02
!      real (kind=R8) ::                                                          !srv
!     *   pra(-1:nx,-1:ny), prac(-1:nx,-1:ny), nac(-1:nx,-1:ny)        !srv 30.01.00
!      real (kind=R8) ::
!     *   bbc(-1:nx,-1:ny,0:3)
!      real (kind=R8) ::
!     *   gampvd(-1:nx,-1:ny),gamrvd(-1:nx,-1:ny)
!      real (kind=R8) ::
!     *   gampvdc(-1:nx,-1:ny)                                        !srv
      real (kind=R8) ::
     *   dvpar(-1:nx,-1:ny)  !,dvdia(-1:nx,-1:ny)                     !srv
     *  , wrk0(-1:nx,-1:ny), qip0(-1:nx,-1:ny)                        !srv 11.07.99
     *  , qip1(-1:nx,-1:ny), qipgen(-1:nx,-1:ny)                      !srv 20.12.99
      real (kind=R8) :: facdriftm
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)               !srv 16.09.02
!      real (kind=R8) ::
!     *   dvparc(-1:nx,-1:ny),dvdiac(-1:nx,-1:ny)

      save ncall, fch_ion_neutral, fch_anomalous, dia_cur, fch_inert, !srv 07.07.99
     &     vis_par, vis_per, vis_q, vdia_par                          !srv 17.05.99 11.07.99 22.11.99
      integer, save :: model = 1
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   harm, hybr, hy1, upwind
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                         !srv 11.07.99
csw aug2014 fix hy1 as suggested by xpb (multiply qz instead divide)
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                              !srv 07.05.08
      harm(t0,t1) = 2*t0*t1/(t0+t1)
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
csw 2013 uwpind scheme (for modified hybrid, SOLPS4 style)
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
*   ..initialization
      data ncall/0/, fch_ion_neutral/1.0e0_R8/, fch_inert/1.0e0_R8/      !srv 07.07.99 22.03.02    
     *        , fch_anomalous/1.0e0_R8/                                  !srv 22.03.02
     *        , dia_cur/1.0e0_R8/, vis_par/1.0e0_R8/, vdia_par/1.0e0_R8/ !srv 17.05.99 22.03.02
     *        , vis_q/1.0e0_R8/, vis_per/1.0e0_R8/                       !srv 11.07.99 22.11.99 22.03.02

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhe_')
*   ..set internal parameters on first call                              !srv
      if (ncall.eq.0) then                                               !srv
       call ipgetr ('b2tfhe_neutral', fch_ion_neutral)                   !srv
       call ipgetr ('b2tfhe_anomalous', fch_anomalous)                   !srv
       call ipgetr ('b2tfhe_dia_cur', dia_cur)                           !srv
       call ipgetr ('b2tfhe_vis_par', vis_par)                           !srv
       call ipgetr ('b2tfhe_vis_per', vis_per)                           !srv 22.11.99
       call ipgetr ('b2tfhe_vis_q', vis_q)                               !srv 11.07.99
       call ipgetr ('b2tfhe_vdia_par', vdia_par)                         !srv 17.05.99
       call ipgetr ('b2tfhe_inert', fch_inert)                           !srv 07.07.99
       call ipgeti ('b2tfhe_lim_flux', lim_flux)                         !srv 04.02.03
       call ipgetr ('b2tfhe_stochastic', fch_stochastic)                 !lk 31.01.08
       if(fch_ion_neutral .eq. 0.0e0_R8 .and.                            !srv
     *      fch_anomalous .eq. 0.0e0_R8) then                            !srv
         stop 'fch_ion_neutral .eq. 0.0 .and. fch_anomalous .eq. 0.0'    !srv
       endif                                                             !srv
       call ipgeti ('b2tfhe_iout', iout)                                 !srv 16.09.02
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                  !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                  !srv 06.07.06
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
       call ipgeti ('b2tfhe_mdf', mdf_fhe)                               !srv 13.10.06
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
       call ipgeti ('b2siav_style_qip', qstyle)                          !lk 23.05.06 05.07.07
       call ipgetr ('b2tfhe_alfTeEh', alfTeEh)                           !srv 28.04.08
       call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                           !srv 30.12.09
       if(fch_pTe .eq. 1.0_R8) alfTeEh = 0.0_R8                          !srv 16.04.10
       call xertst (.not.(0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe),
     &        ' Do not set 0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe')
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2tfhe_mode_ehx',mode_ehx)                          !srv 21.10.09
       call ipgeti ('b2trcl_conductive limit', conductive limit)         !srv 16.06.08 17.01.12
       call ipgeti ('b2tfhe_no_current',no_current)
       call ipgeti ('b2tqce_model', model)
       call ipgeti ('b2tfhe_upwind',b2_upwind)
       call ipgeti ('b2tfhe_hybr2',b2tfhe_hybr2)
       call ipgeti ('b2tfhe_flo53',flo53)
       endif                                                              !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fne
       call b2xvff (nx, ny, fne, 'fne')
*    ..test sign of ne, te, csig, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, csig, 1, 'csig', '.ge.')
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
*    ..test edge values of csig, chce, chve, chvemx, calf
       call b2xvff (nx, ny, csig, 'csig')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvfx (nx, ny, chvemx, 'chvemx')
       call b2xvff (nx, ny, calf, 'calf')
      endif

      facdriftm = maxval(facdrift)
*   ..calculate the contribution of ion-neutral current                 !srv 11.02.04 {
      fchin     = 0.0e0_R8                                             !srv 09.07.99
      if(pot_eq.ne.0) then                                             !srv 26.07.06
       if (fch_ion_neutral .ne. 0.0e0_R8) then                             !srv 01.07.99
        call b2tinnt (nx, ny, ns, ismain, vol, gs, hx, hy, qz, bb,  
     1                   te, po, ti, na, fna, csigin, fchin)            !srv 11.05.09
       endif                                                            !srv 11.02.04 }
      endif                                                              !srv 26.07.06
c Drift terms
      fchdia    =0.0e0_R8
      fchvispar =0.0e0_R8                                             !srv 09.07.99
      fchvisper =0.0e0_R8                                             !srv 22.11.99
      fchvisq   =0.0e0_R8                                             !srv 11.07.99
      fchinert  =0.0e0_R8                                             !srv 09.07.99
      avisper  = 0.0e0_R8; bvisper = 0.0e0_R8; cvisper = 0.0e0_R8     !srv 27.01.00 22.03.02 
      if(facdriftm.ne.0) then
        if (dia_cur.ne.0.0e0_R8) then                                  !srv
*       X-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT             !srv

!         do iy=-1,ny                                                                                              !srv
         do iy=0,ny-1                                                                                              !srv
           do ix=0,nx                                                                                                  !srv
            if(leftix(ix,iy).ne.-2) then                                 !srv 26.11.02
             fchdia(ix,iy,0)=fchdia(ix,iy,0)+                                                  !srv
     &                       dia_cur* ! signmf*                          !srv 25.04.00 05.02.05 11.01.13
     &      (pz(ix,iy)*bb(ix,iy,2)+                                      !srv 31.03.00 18.01.02
     &        pz(leftix(ix,iy),leftiy(ix,iy))                            !srv 31.03.00
     &       *bb(leftix(ix,iy),leftiy(ix,iy),2))*0.5e0_R8                !srv 31.03.00
     &      *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))               !srv
     &      /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))             !srv
     &      *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))         !srv
!     &       (ne(ix,iy)+ne(leftix(ix,iy),leftiy(ix,iy)))/2.0e0_R8      !srv
!     &      *(te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))                              !srv
!     &       +ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/2.0e0_R8      !srv
!     &      *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0e0_R8  !srv
!     &      *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))             !srv
!     &      /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))           !srv
!     &      *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))       !srv
            endif                                                        !srv 26.11.02
           enddo                                                                                                               !srv
         enddo                                                                                                         !srv 

*       Y-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT                   !srv

         do iy=0,ny                                                                                                    !srv
           do ix=0,nx-1                                                                                            !srv
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
!           do ix=-1,nx                                                                                            !srv
!             fchdia(ix,iy,1)=fchdia(ix,iy,1) - dia_cur*signmf*            !srv 25.04.00 21.04.00 21.03.02
             fchdia(ix,iy,1)= - dia_cur* ! signmf*                          !srv 25.04.00 21.04.00 21.03.02 11.01.13
!     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))*pz(ix,iy)*bb(ix,iy,2)+ !srv 31.03.00 18.01.02 { interpolaion for bottom
!     &       vol(ix,iy)*pz(bottomix(ix,iy),bottomiy(ix,iy))                        !srv 31.03.00
!     &       *bb(bottomix(ix,iy),bottomiy(ix,iy),2))/            !srv 31.03.00
!     &      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))            !srv 18.01.02 }
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))           !srv
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))         !srv
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))    !srv
     &       (pz(ix,iy)*bb(ix,iy,2)+                                     !srv 31.03.00 18.01.02 {
     &        pz(bottomix(ix,iy),bottomiy(ix,iy))                        !srv 31.03.00
     &       *bb(bottomix(ix,iy),bottomiy(ix,iy),2))*0.5e0_R8            !srv 31.03.00 18.01.02 }
     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))           !srv
     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))         !srv
     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))     !srv

!     &      -(ne(ix,iy)+ne(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8  !srv
!     &      *(te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy))                       !srv
!     &       +ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8  !srv
!     &      *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))         !srv
!     ^       /2.0e0_R8                                                 !srv
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))           !srv
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))         !srv
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))    !srv
            endif                                                        !srv 26.11.02
           enddo                                                                                                                   !srv
         enddo                                                                                                             !srv

        endif
*  ..calculate  X,Y-COMPONENTS OF VISCOSITY PARALLEL CURRENT
        if(vis_par .ne. 0.0e0_R8) 
     &   call b2tvspa (nx, ny, ns, vol, hx, hy, qz, bb, ua, 
     &       cvsa, facdrift, fchvispar) !srv 20.09.06
 
* ..calculate the contribution of viscosity perpendicular current        !srv 22.11.99 
        va=0.0e0_R8
         dva=0.0e0_R8

        if(vis_per.ne.0.0e0_R8) then                                        !srv 22.11.99

!      analogous part for avisper and va takes place in b2siav          !srv 17.01.01

         is = ismain                                         !srv 23.11.10

!         avisper(:,-1,1)=0.0e0_R8                          !srv 06.05.02 {
!         do iy=0,ny    
!          do ix=-1,nx
!               avisper(ix,iy,1)=
!!     *           (vol(ix,iy)/(hy(ix,iy)**2*bb(ix,iy,3))+
!!     *          vol(bottomix(ix,iy),bottomiy(ix,iy))/
!!     *          (hy(bottomix(ix,iy),bottomiy(ix,iy))**2*
!!     *           bb(bottomix(ix,iy),bottomiy(ix,iy),3))
!!     *         )*0.5e0_R8
!     *      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))*4.0e0_R8
!     *         /(hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))**2
!     *         /(bb(ix,iy,3)+bb(bottomix(ix,iy),bottomiy(ix,iy),3))
!          enddo
!         enddo                                                !srv 06.05.02 }

         avisper(:,-1,1) = 0.0e0_R8                                  !srv 06.05.02
         cvisper(:,-1,1) = 0.0e0_R8
         do iy=0,ny
          do ix=-1,nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then   !srv 26.11.02
             cvisper(ix,iy,1)=
     &         4.0e0_R8/(bb(ix,iy,3)+
     1         bb(bottomix(ix,iy),bottomiy(ix,iy),3))/
     &         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
             avisper(ix,iy,1) = cvisper(ix,iy,1)                      !srv 06.05.02
           endif                                                      !srv 26.11.02
          enddo
         enddo

!         do is = 0,ns-1
!          if(za(is).gt.0.0e0_R8) then

         do iy=-1,ny
          do ix=-1,nx
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then   !srv 26.11.02
!            bvisper(ix,iy,1) = 0.25e0_R8*vsa0(ix,iy,is)/hy(ix,iy)  !srv 06.05.02
            bvisper(ix,iy,1) = 0.25e0_R8*vsa0(ix,iy,is)            !srv 06.05.02
     &                                    *vol(ix,iy)/hy1(ix,iy)**2   !srv 06.05.02
!            pra(ix,iy) = na(ix,iy,is)*ti(ix,iy) 
           endif                                                      !srv 26.11.02
          enddo
         enddo

!          call intcor(nx,ny,vol,pra,prac)
!            call intcor(nx,ny,vol,na(-1,-1,is),nac)

         do iy=0,ny
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then   !srv 26.11.02
              va(ix,iy,0)=-cvisper(ix,iy,1)
     &        *(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
!     &        -cdna(ix,iy,1,is)*2.0e0_R8*
!     &        (nac(ix,iy)-nac(rightix(ix,iy),rightiy(ix,iy)))/
!     &        (hx(ix,iy)*na(ix,iy,is)+
!     &        hx(bottomix(ix,iy),bottomiy(ix,iy))*
!     &        na(bottomix(ix,iy),bottomiy(ix,iy),is))
!     &        -cdpa(ix,iy,1,is)*2.0e0_R8*
!     &        (prac(ix,iy)-prac(rightix(ix,iy),rightiy(ix,iy)))/
!     &        (hx(ix,iy)*na(ix,iy,is)+
!     &        hx(bottomix(ix,iy),bottomiy(ix,iy))*
!     &        na(bottomix(ix,iy),bottomiy(ix,iy),is))
     &         -cvisper(ix,iy,1)*(na(ix,iy,is)*ti(ix,iy)
     &         -na(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &         ti(bottomix(ix,iy),bottomiy(ix,iy)))*
     &         (2.0e0_R8/(rza(ix,iy,is)*qe))
     &         /(na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))
           endif                                                      !srv 26.11.02
          enddo
         enddo

         do iy=0,ny-1
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then   !srv 26.11.02
            dva(ix,iy)=bvisper(ix,iy,1)
     *       *(va(topix(ix,iy),topiy(ix,iy),0)-va(ix,iy,0))
           endif                                                      !srv 26.11.02
          enddo
         enddo

*          Y-COMPONENT OF VISCOSITY PERPENDICULAR CURRENT                 !srv 22.11.99
         fchvisper(:,0,1) = 0.0e0_R8      ! because of boundary condition !srv 28.01.00
         do iy=1,ny-1    
          do ix=0,nx-1
           if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then   !srv 26.11.02
            fchvisper(ix,iy,1)=fchvisper(ix,iy,1)
     *       -avisper(ix,iy,1)*vis_per
     *       *(dva(ix,iy)-dva(bottomix(ix,iy),bottomiy(ix,iy)))
           endif                                                      !srv 26.11.02
          enddo
         enddo
         fchvisper(:,ny,1) = 0.0e0_R8     ! because of boundary condition !srv 28.01.00
*          X-COMPONENT OF VISCOSITY PERPENDICULAR CURRENT                 !srv 22.11.99
         fchvisper(:,:,0)=fchvisper(:,:,0)+0.0e0_R8
!          endif
!         enddo
        endif                                                            !srv 22.11.99
*
* ..calculate the parallel heat flux of ion at the center of cell        !lk 23.05.06 05.07.07 {

        if(vis_q .ne. 0.0e0_R8) then
         do is = 0,ns-1                                                    !srv 12.05.11
          if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!           is = ismain    ! 1                                  
           qip0  = 0.0e0_R8
           qip1  = 0.0e0_R8
           qipgen= 0.0e0_R8
           qip0a = 0.0e0_R8
           qip0b = 0.0e0_R8
           wrk0 = 0.0e0_R8

           call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)                    !srv 11.07.99
           
           tt0 = 0.0e0_R8                                                 
           tt1 = 0.0e0_R8
           tt2 = 0.0e0_R8                                                 
           tt3 = 0.0e0_R8
           conn = 0.0e0_R8
           bbaverage2 = 0.0e0_R8
           eps1 = 0.0e0_R8

           bbaverage2 = 0.0e0_R8                                         
           mintopcut = minval(topcut(1:nncut))                           
        
           do iy = 0,mintopcut    ! for inner surfaces and first open surface                                                                                           
            do ix=0,nx                                                  
             if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then 
              tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
              tt1(iy) = tt1(iy) +vol(ix,iy)
              conn(iy) = conn(iy)+hx(ix,iy)*abs(bb(ix,iy,3)/bb(ix,iy,0)) !srv 11.01.13
              tt2(iy) = tt2(iy) +hx(ix,iy)
              tt3(iy) = tt3(iy) +vol(ix,iy)/hy1(ix,iy)
             endif
            enddo
           enddo

           do iy = 0,mintopcut                                                                                  
            bbaverage2(iy) = tt0(iy)/tt1(iy)
            eps1(iy) = tt2(iy)**2/tt3(iy)
           enddo                                                        
         
           if (qstyle .eq. 0) then

            do iy = 0,mintopcut-1    ! for inner surfaces                                                                           
             do ix = 0,nx                                                  
              if(mod(region(bottomix(ix,iy),
     .                      bottomiy(ix,iy),0),4).eq.1) then 
               qip0a(ix,iy) = (-2.5e0_R8)*(bb(ix,iy,3)*bb(ix,iy,2)
     &          /bb(ix,iy,0)/bbaverage2(iy))
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*
     &          (ti(topix(ix,iy),topiy(ix,iy))
     &          -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &          0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy)))    
              endif
             enddo
            enddo
       
            do iy = 0,mintopcut-1                                                                                                   
             do ix = 0,nx                                                  
              if(mod(region(bottomix(ix,iy),
     .                      bottomiy(ix,iy),0),4).eq.1) then                         
               qip0b(ix,iy) = 0.25e0_R8*(qip0a(ix,iy)*2.0e0_R8
     &          +qip0a(leftix(ix,iy),leftiy(ix,iy))
     &          +qip0a(rightix(ix,iy),rightiy(ix,iy)))
              endif
             enddo
            enddo
                  
            do iy = 0,mintopcut-1                                                                                                   
             do ix = 0,nx                                                  
              if(mod(region(bottomix(ix,iy),
     .                      bottomiy(ix,iy),0),4).eq.1) then                         
               qip0(ix,iy) = 0.25e0_R8*(qip0b(ix,iy)*2.0e0_R8
     &          +qip0b(leftix(ix,iy),leftiy(ix,iy))
     &          +qip0b(rightix(ix,iy),rightiy(ix,iy)))
              endif
             enddo
            enddo
       
           endif

           if (qstyle.eq.1) then

            qip0a = 0.0e0_R8
            do iy=0,ny-1                                                 
             do ix=-1,nx                                                 
              if(leftix(ix,iy).ne.-2) then                               
               qip0a(ix,iy)=-chci(ix,iy,0)
     &          *(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
              endif                                                      
             enddo                                                       
            enddo 
         
            do iy=0,ny-1                                                 
             do ix=-1,nx-1                                                
              if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then                         
               qip0(ix,iy)=0.5e0_R8*(qip0a(ix,iy)+
     &          qip0a(rightix(ix,iy),rightiy(ix,iy)))
     &          *hx(ix,iy)*bb(ix,iy,3)/bb(ix,iy,0)/vol(ix,iy)
     &          -2.5e0_R8*bb(ix,iy,2)
     &          /bb(ix,iy,0)/bb(ix,iy,3)
     &          *(ni(ix,iy)*ti(ix,iy)/qe)*
     &          (ti(topix(ix,iy),topiy(ix,iy))
     &          -ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &          (hy1(ix,iy)+0.5e0_R8*hy1(topix(ix,iy),topiy(ix,iy))+
     &          0.5e0_R8*hy1(bottomix(ix,iy),bottomiy(ix,iy))) 
              endif                                                      
             enddo                                                       
            enddo 
                                                                  
           endif
        
           do iy=0,ny-1                                                 
            do ix=0,nx-1                                                
             qip1(ix,iy)=cqip1*ni(ix,iy)*ti(ix,iy)*
     &        sqrt(te(ix,iy)+ti(ix,iy))/sqrt(am(is)*mp)
            enddo                                                       
           enddo                                                        

           do iy=0,ny-1                                                 
            do ix=0,nx-1                                                
             wrk0(ix,iy)=qip1(ix,iy)*datan(qip0(ix,iy)/qip1(ix,iy))
     &        *sqrt(bb(ix,iy,3))
            enddo                                                       
           enddo                                                          

           do iy = 0,mintopcut                                                                                                   
            do ix = 0,nx                                                  
             if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then                         
              qipgen(ix,iy)=(wrk0(ix,iy)
     &         +wrk0(leftix(ix,iy),leftiy(ix,iy)))*0.5e0_R8
             endif                                                      
            enddo
           enddo                                                                                                                  

           tauia1 = 0.0e0_R8
           nu1 = 0.0e0_R8
           alpha1 = 0.0e0_R8

           do iy = 0,mintopcut                                                                                            
            do ix=0,nx
             if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then                         
              tauia1(ix,iy) = tauia(ix,iy,is)/                             !lk 20.11.07 {
     &         (1.0e0_R8+15.12e16_R8/conn(iy)*
     &         (ti(ix,iy)/ev)**2/ne2(ix,iy))
     &         /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     &         /conn(iy)*(ti(ix,iy)/ev)**2/ne2(ix,iy))
              nu1(ix,iy) = conn(iy)/(tauia(ix,iy,is)
     &         *sqrt(ti(ix,iy)/(mp*am(is)))
     &         *2.0e0_R8*pi*eps1(iy)**1.5e0_R8) 
              alpha1(ix,iy) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu1(ix,iy))+
     &         2.7e0_R8*nu1(ix,iy)**2*eps1(iy)**3)/(1.0e0_R8+
     &         0.7e0_R8*sqrt(nu1(ix,iy))+nu1(ix,iy)**2*eps1(iy)**3)
     &         -1.0e0_R8)
     &         *0.533e0_R8                                                   !lk 20.11.07 }
             endif    
            enddo
           enddo                                                                                                            !lk 23.05.06 05.07.07 }

* ..calculate the x-component of the parallel viscosity current
*                                  connected with the heat flux 

           do iy = 0,mintopcut-1                                                                                 !lk 23.05.06 05.07.07 {                                                                                             
            do ix = 0,nx                                                  
             if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then                         
              fchvisq(ix,iy,0)=-0.24e0_R8*vis_q*(                ! signmf*                !srv 11.01.13
     &         vol(ix,iy)/hx(ix,iy)**2/hy1(ix,iy)
     &         *tauia1(ix,iy)
     &         *alpha1(ix,iy)
     &         *bb(ix,iy,2)/bb(ix,iy,3)
     &         *bb(ix,iy,0)/sqrt(bb(ix,iy,3))
     &         +vol(leftix(ix,iy),leftiy(ix,iy))
     &         /hx(leftix(ix,iy),leftiy(ix,iy))**2
     &         /hy1(leftix(ix,iy),leftiy(ix,iy))
     &         *tauia1(leftix(ix,iy),leftiy(ix,iy))
     &         *alpha1(leftix(ix,iy),leftiy(ix,iy))
     &         *bb(leftix(ix,iy),leftiy(ix,iy),2)
     &         /bb(leftix(ix,iy),leftiy(ix,iy),3)
     &         *bb(leftix(ix,iy),leftiy(ix,iy),0)
     &         /sqrt(bb(leftix(ix,iy),leftiy(ix,iy),3))
     &         )*0.5e0_R8
     &         *(wrk0(ix,iy)-wrk0(leftix(ix,iy),leftiy(ix,iy)))
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))
             endif                                                          !lk 23.05.06 05.07.07 }
            enddo
           enddo         

* ..calculate the y-component of the parallel viscosity current
*                                  connected with the heat flux 
           do iy = 0,mintopcut                                                                 !lk 23.05.06 05.07.07 {                                                                                               
            do ix = 0,nx                                                  
             if(mod(region(bottomix(ix,iy),
     .                     bottomiy(ix,iy),0),4).eq.1) then                         
              fchvisq(ix,iy,1)=0.24e0_R8*vis_q*(             ! signmf*      !srv 11.01.13                    
     &         vol(ix,iy)/hx(ix,iy)**2/hy1(ix,iy)
     &         *tauia1(ix,iy)
     &         *alpha1(ix,iy)
     &         *bb(ix,iy,2)/bb(ix,iy,3)
     &         *bb(ix,iy,0)/sqrt(bb(ix,iy,3))
     &         +vol(bottomix(ix,iy),bottomiy(ix,iy))
     &         /hx(bottomix(ix,iy),bottomiy(ix,iy))**2
     &         /hy1(bottomix(ix,iy),bottomiy(ix,iy))
     &         *tauia1(bottomix(ix,iy),bottomiy(ix,iy))
     &         *alpha1(bottomix(ix,iy),bottomiy(ix,iy))
     &         *bb(bottomix(ix,iy),bottomiy(ix,iy),2)
     &         /bb(bottomix(ix,iy),bottomiy(ix,iy),3)
     &         *bb(bottomix(ix,iy),bottomiy(ix,iy),0)
     &         /sqrt(bb(bottomix(ix,iy),bottomiy(ix,iy),3))
     &                      )*0.5e0_R8
     &         *(qipgen(rightix(ix,iy),rightiy(ix,iy))
     &         +qipgen(rightix(bottomix(ix,iy),bottomiy(ix,iy)),
     &             rightiy(bottomix(ix,iy),bottomiy(ix,iy)))
     &         -qipgen(ix,iy)-qipgen(bottomix(ix,iy),bottomiy(ix,iy))
     &            )*0.5e0_R8
     &         *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))
     &         -OnedBsqc(ix,iy))
             endif                                                         !lk 23.05.06 05.07.07 }
            enddo
           enddo   
          endif
         enddo
        endif                                                            
!*   ..calculate the contribution of ion-neutral current
*     INERTIA TERMS                                                      !srv 07.07.99
        call b2tiner (nx, ny, ns, vol, gs, hx, hy, bb, na, ua,    !srv 22.11.12
     &                     facdrift,fchinert)
        do iy = -1, ny
         do ix = -1, nx
          fchdia(ix,iy,:) = fchdia(ix,iy,:)*facdrift(ix,iy)
!          fchvispar(ix,iy,:) = fchvispar(ix,iy,:)*facdrift(ix,iy)        !srv          25.04.00 { 22.11.12 multyplied in b2tvspa.F
          fchvisper(ix,iy,:) = fchvisper(ix,iy,:) !*facdrift(ix,iy)      !srv 22.11.99 07.02.03
          fchvisq(ix,iy,:)   = fchvisq(ix,iy,:)*facdrift(ix,iy)          !srv 11.07.99
!          fchin(ix,iy,:)     = fchin(ix,iy,:) !*facdrift(ix,iy)          !srv 09.07.99 07.02.03
!          fchinert(ix,iy,:)  = fchinert(ix,iy,:)*facdrift(ix,iy)         !srv 09.07.99 25.04.00 } 22.11.12 multyplied in b2tiner.F
         enddo
        enddo
      endif !facdriftm!

      if(fch_anomalous.ne.0.0e0_R8 ) then                              !srv 23.03.07 {
       call b2tanml (nx, ny, ns, vol, hx, hy, csig_an, po, fchanml)
      endif                                                            !srv 23.03.07 }

      if(fch_stochastic.ne.0.0e0_R8 ) then                             !lk 31.01.08 { 
       call b2tstch (nx, ny, vol, hy, bb, po, te, ne, csig_stoch, 
     &      fchstoch)
      else
       fchstoch = 0.0e0_R8
      endif                                                            !lk 31.01.08 }
       
      fledia = 0.0_R8                                                  !srv 25.01.02 {
      if(facdriftm.ne.0.0e0_R8) then                            !srv 18.01.07
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!         is = ismain   ! 1    ! for is = 1 (za = 1)                        !srv 12.05.11
* ..compute flediax
         do iy = -1, ny
          do ix = 0, nx
           if(leftix(ix,iy).ne.-2) then                                 !srv 26.11.02
            fledia(ix,iy,0) = 
     &       -(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))        !srv 25.01.02 was 5/2
     &       /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))          !srv 17.05.99 11.03.00
     &       *(te(leftix(ix,iy),leftiy(ix,iy))                     !srv 17.05.99 11.03.00
     &       /ti(leftix(ix,iy),leftiy(ix,iy))                      !srv 17.05.99 11.03.00
     &       +te(ix,iy)/ti(ix,iy))*0.5e0_R8                       !srv 17.05.99 11.03.00 25.01.02 13.09.02
     &       *vadia(ix,iy,0,is)*vdia_par ! for is = 1 (za = 1)     !srv 17.05.99 11.03.00 21.04.00 25.01.02 18.01.07
     1       *(bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 13.09.02 {
     2       bb(leftix(ix,iy),leftiy(ix,iy),2)/
     3       bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0e0_R8                !srv 13.09.02 }
           endif                                                        !srv 26.11.02
          enddo
         enddo
* ..compute flediay
         do iy = 0, ny
          do ix = -1, nx
           fledia(ix,iy,1) = 
     &      -(vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))    !srv 17.05.99 11.03.00 25.01.02 was 5/2
     &      /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))     !srv 17.05.99 11.03.00
     &      *(te(bottomix(ix,iy),bottomiy(ix,iy))                !srv 17.05.99 11.03.00 13.09.02
     &      /ti(bottomix(ix,iy),bottomiy(ix,iy))                 !srv 17.05.99 11.03.00
     &      +te(ix,iy)/ti(ix,iy))*0.5e0_R8                       !srv 17.05.99 11.03.00 25.01.02 13.09.02
     &      *vadia(ix,iy,1,is)*vdia_par                          !srv 17.05.99 11.03.00 18.01.07
     1      *(bb(ix,iy,2)/bb(ix,iy,3)+                            !srv 13.09.02 {
     2      bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     3      bb(bottomix(ix,iy),bottomiy(ix,iy),3))/2.0e0_R8      !srv 13.09.02 }
          enddo
         enddo                                                        !srv 25.01.02 }
        endif
       enddo
      endif                                                           !srv 18.01.07

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny),
     *    w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                !srv 16.09.02 }

      fhePSch = 0.0e0_R8                                              !srv 13.10.06
      cone    = 0.0e0_R8                                              !srv 02.01.07
      if(facdriftm.ne.0.0e0_R8) then
*   ..compute P.Sch heat x-flux                                 !srv 13.10.06 {
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!        is = ismain !  1                                                 !srv 12.05.11  
!        iy = -1
!        do ix = -1, nx
!         if(leftix(ix,iy).eq.-2) then
!          fhePSch(ix,iy,0)=0.0e0_R8
!         else
!          fhePSch(ix,iy,0)=-0.5e0_R8*(
!     &     cddi(ix,topiy(ix,iy),1,is)*
!     &     (ne(ix,topiy(ix,iy))*te(ix,topiy(ix,iy))**2-
!     &                        ne(ix,iy)*te(ix,iy)**2)+
!     &     cddi(leftix(ix,iy),topiy(ix,iy),1,is)*
!     &      (ne(leftix(ix,iy),topiy(ix,iy))*
!     &      te(leftix(ix,iy),topiy(ix,iy))**2
!     &                    -ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2)
!     &                               )
!          endif
!         enddo

         do iy = 0, ny-1
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
            fhePSch(ix,iy,0)=0.0e0_R8
           else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            fhePSch(ix,iy,0)=-0.25e0_R8*PSch*(
     &       cddi(topix(ix,iy),topiy(ix,iy),1,is)*
     &       (ne(topix(ix,iy),topiy(ix,iy))*
     &       te(topix(ix,iy),topiy(ix,iy))**2-ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(topix(ixl,iyl),topiy(ixl,iyl),1,is)*
     &       (ne(topix(ixl,iyl),topiy(ixl,iyl))*
     &        te(topix(ixl,iyl),topiy(ixl,iyl))**2
     &                  -ne(ixl,iyl)*te(ixl,iyl)**2)+
     &       cddi(ixl,iyl,1,is)*
     &       (ne(ixl,iyl)*te(ixl,iyl)**2
     &                     -ne(bottomix(ixl,iyl),bottomiy(ixl,iyl))*
     &                      te(bottomix(ixl,iyl),bottomiy(ixl,iyl))**2)+
     &       cddi(ix,iy,1,is)*
     &       (ne(ix,iy)*te(ix,iy)**2-ne(bottomix(ix,iy),bottomiy(ix,iy))
     &                          *te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            cone(ix,iy,0,1) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*ne(ix,iy)*  ! 0.25*2  !srv 02.01.07
     *                                                te(ix,iy)
            cone(ix,iy,0,2) = -0.5e0_R8*PSch*cddi(ix,iy,1,is)*            ! 0.25*2  !srv 02.01.07
     &                      ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                      te(bottomix(ix,iy),bottomiy(ix,iy))
           endif
          enddo
         enddo

!         iy = ny
!          do ix = -1, nx
!           if(leftix(ix,iy).eq.-2) then
!            fhePSch(ix,iy,0)=0.0e0_R8
!           else
!            fhePSch(ix,iy,0)=-0.5e0_R8*(
!     &       cddi(leftix(ix,iy),iy,1,is)*
!     &       (ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2
!     &                         -ne(leftix(ix,iy),bottomiy(ix,iy))*
!     &                            te(leftix(ix,iy),bottomiy(ix,iy))**2)+
!     &       cddi(ix,iy,1,is)*
!     &       (ne(ix,iy)*te(ix,iy)**2
!     &                -ne(ix,bottomiy(ix,iy))*te(ix,bottomiy(ix,iy))**2)
!     &                               )
!           endif
!          enddo

*   ..compute P.Sch y-flux
         iy = -1
         fhePSch(:,iy,1)=0.0e0_R8
         do iy = 0, ny
          do ix = -1, nx
!           if(leftix(ix,iy).eq.-2) then
!            fhePSch(ix,iy,1) = 0.5e0_R8*(
!     &       cddi(rightix(ix,iy),iy,0,is)*
!     &        (ne(rightix(ix,iy),iy)*te(rightix(ix,iy),iy)**2
!     &                                        -ne(ix,iy)*te(ix,iy)**2)+
!     &       cddi(rightix(ix,iy),bottomiy(ix,iy),0,is)*
!     &       (ne(rightix(ix,iy),bottomiy(ix,iy))*
!     &             te(rightix(ix,iy),bottomiy(ix,iy))**2
!     &                -ne(ix,bottomiy(ix,iy))*te(ix,bottomiy(ix,iy))**2)
!     &                               )
!           elseif(rightix(ix,iy).eq.nx+1) then
!            fhePSch(ix,iy,1) = 0.5e0_R8*(
!     &       cddi(ix,iy,0,is)*
!     &       (ne(ix,iy)*te(ix,iy)**2
!     &                   -ne(leftix(ix,iy),iy)*te(leftix(ix,iy),iy)**2)+
!     &       cddi(ix,bottomiy(ix,iy),0,is)*
!     &       (ne(ix,bottomiy(ix,iy))*te((ix),bottomiy(ix,iy))**2
!     &                         -ne(leftix(ix,iy),bottomiy(ix,iy))*
!     &                            te(leftix(ix,iy),bottomiy(ix,iy))**2)
!     &                               )
!           else
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fhePSch(ix,iy,1) = 0.0e0_R8
           else
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            fhePSch(ix,iy,1) = 0.25e0_R8*PSch*(
     &       cddi(rightix(ix,iy),rightiy(ix,iy),0,is)*
     &       (ne(rightix(ix,iy),rightiy(ix,iy))*
     &       te(rightix(ix,iy),rightiy(ix,iy))**2
     &                                         -ne(ix,iy)*te(ix,iy)**2)+
     &       cddi(ix,iy,0,is)*
     &       (ne(ix,iy)*te(ix,iy)**2
     &                   -ne(leftix(ix,iy),leftiy(ix,iy))*
     &                    te(leftix(ix,iy),leftiy(ix,iy))**2)+
     &       cddi(bottomix(ix,iy),bottomiy(ix,iy),0,is)*
     &       (ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &       te(bottomix(ix,iy),bottomiy(ix,iy))**2
     &                         -ne(leftix(ixb,iyb),leftiy(ixb,iyb))*
     &                          te(leftix(ixb,iyb),leftiy(ixb,iyb))**2)+
     &       cddi(rightix(ixb,iyb),rightiy(ixb,iyb),0,is)*
     &       (ne(rightix(ixb,iyb),rightiy(ixb,iyb))*
     &                       te(rightix(ixb,iyb),rightiy(ixb,iyb))**2
     &                      -ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &                       te(bottomix(ix,iy),bottomiy(ix,iy))**2)
     &                               )
            cone(ix,iy,1,1) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*ne(ix,iy)*   ! 0.25*2 !srv 02.01.07 {
     &                                                te(ix,iy)
            cone(ix,iy,1,2) = -0.5e0_R8*PSch*cddi(ix,iy,0,is)*             ! 0.25*2
     &                        ne(leftix(ix,iy),leftiy(ix,iy))*
     &                        te(leftix(ix,iy),leftiy(ix,iy))      !srv 02.01.07 }
           endif
          enddo
         enddo

         do iy = -1, ny
          do ix = -1, nx                                                  !srv 26.11.02 {
           fhePSch(ix,iy,:)=2.5e0_R8*fhePSch(ix,iy,:)*facdrift(ix,iy)
           cone(ix,iy,:,1:2) = cone(ix,iy,:,1:2)*facdrift(ix,iy)                           !srv 02.01.07
!           cone(ix,iy,:,2) = cone(ix,iy,:,2)*facdrift(ix,iy)                           !srv 02.01.07
          enddo
         enddo
        endif
       enddo

c  ..set null flux in guard cells  ?
!      fhePSch(-1, :,:)=0.0e0_R8    ! MODIFY for Double Case 
!      fhePSch(nx, :,1)=0.0e0_R8
!      fhePSch( :,-1,:)=0.0e0_R8
!      fhePSch( :,ny,0)=0.0e0_R8
      endif

      fhe_mdf = 0.0e0_R8                                              !srv 13.10.06
 
* ..compute fluxes and linearization
*   ..compute fhex, floex, conex, fchx, concx
      do iy = -1, ny
       do ix = -1, nx                                                  !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2)then                                    !srv 26.11.02
         floe(ix,iy,0)    = 0.0e0_R8
         cone(ix,iy,0,0)  = 0.0e0_R8                                   !srv 02.01.07
         fhe(ix,iy,0)     = 0.0e0_R8
         fhe_mdf(ix,iy,0) = 0.0e0_R8                                   !srv 13.10.06 
         conc(ix,iy,0)    = 0.0e0_R8
         fch(ix,iy,0)     = 0.0e0_R8
         fch_p(ix,iy,0)   = 0.0e0_R8
        else                                                           !srv 26.11.02 }
         if(mode_ehx.eq.0) then
          ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &      (ne(ix,iy)*te(ix,iy)                                      !srv .98
     &      -ne(leftix(ix,iy),leftiy(ix,iy))*                         !srv
     &      te(leftix(ix,iy),leftiy(ix,iy)))*                         !srv
     &      (2/qe)/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))        !srv
         else
          ehx = -(po(ix,iy)-po(leftix(ix,iy),leftiy(ix,iy)))+
     &    (                                                             !srv .98
     &     te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))+                   !srv
     &     (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))*0.5e0_R8*        !srv
     &     (log(ne(ix,iy))-log(ne(leftix(ix,iy),leftiy(ix,iy))))        !srv 
     &    )/qe                                                          !srv .98 
         endif
         tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
         fch_p(ix,iy,0) = (csig(ix,iy,0)*ehx-calf(ix,iy,0)*tex)        !srv 30.12.09 }
         if(pot_eq.ne.1 .and. no_current.ne.0) then                    !srv 27.02.13
          fch_p = 0.0e0_R8
          ehx = 0.0_R8                                                 !srv 27.02.13
         endif                                                         !srv 30.12.09 }
         floe0 = 1.5e0_R8*fne(ix,iy,0)+
     &    (chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*                   !srv 13.09.02 18.01.07 11.01.13  *signmf
     &    (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &    -0.71e0_R8*fch_p(ix,iy,0)/qe*fch_pTe                        !srv 30.12.09
         floe0_mdf = 1.5e0_R8*fne(ix,iy,0)+                            !srv 13.06.06 {
     &    chve(ix,iy,0)*
     &    (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.06.06 }
     &    -0.71e0_R8*fch_p(ix,iy,0)/qe*fch_pTe                        !srv 30.12.09


         if(lim_flux.eq.1 .and. conductive limit.eq.0) then           !srv 04.02.03 17.01.12

          if(no_Ptncr_x_he.eq.0) then
           cone0 = hybr(floe0,chce(ix,iy,0))-abs(floe0/2.0)
          else
           cone0 = chce(ix,iy,0)-abs(floe0/2.0)
          endif

          if (chce(ix,iy,0).eq.0) then
           alfe0 = 0
          else
           alfe0 = calf(ix,iy,0)*cone0/chce(ix,iy,0)
          endif

csw 05oct2012 check for solps4 model
          if(model.eq.1.or.model.eq.3) then   ! SOLPS5.0 treatments
           flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &       ne(ix,iy))/2
           if(flomx.gt.0.0_R8) then
             t0 = cone0*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
           else
             t0=0.0_R8
           endif
          else
csw 2014 changed flux limiting expression (no t0 cutoff for neg. flomx)
c        i.e.for  model .eq. 2 (SOLPS4)
           flomx = chvemx(ix,iy)
           t0 = cone0*abs(tex/(flomx+1.d-60))
          endif 

          tex0 = tex/(1.d0+t0)

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0)
          if(b2tfhe_hybr2 .ne. 0) then

csw 30apr2014 flo53 ONLY FOR TESTING (SOLPS4)
           if(flo53 .ne. 0) then
            floe0=floe0*5.d0/3.d0
           endif

           cone0=dim(chce(ix,iy,0)/(1.d0+t0), abs(floe0)*
     .      max(0.0d0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .      max(cutlo,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
          endif


csw 27nov2012 add check for upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
             fhe(ix,iy,0) =
     &        upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)), te(ix,iy))
     &         -cone0*(te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))) 
     &         + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
          else
             fhe(ix,iy,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.d0
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
          endif

          if(mdf_fhe.eq.0) then                                         !srv 17.06.09 {

csw 30apr2014
           if(b2tfhe_hybr2 .ne.0) then
            floe(ix,iy,0) = floe0
            cone(ix,iy,0,0) = cone0
           else
            floe(ix,iy,0) = floe0
            cone(ix,iy,0,0) = abs(floe0/2.0)+cone0/(1+t0)                  !srv 02.01.07
           endif


          else

csw 17feb2012
           if(no_Ptncr_x_he.eq.0) then
            cone0_mdf=hybr(floe0_mdf,chce(ix,iy,0))-abs(floe0_mdf/2.0)
           else
            cone0_mdf = chce(ix,iy,0)-abs(floe0_mdf/2.0)
           endif

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0_mdf)
           if(b2tfhe_hybr2 .ne. 0) then
            cone0_mdf=dim(chce(ix,iy,0), abs(floe0_mdf)*
     .        max(0.0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .        max(1.e-30,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
           endif

           if (chce(ix,iy,0).eq.0) then
            alfe0_mdf = 0
           else
            alfe0_mdf = calf(ix,iy,0)*cone0_mdf/chce(ix,iy,0)
           endif

           t0_mdf = cone0_mdf*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
           tex0_mdf = tex/(1.0_R8+t0_mdf)

csw 27nov2012 added check for upwind (SOLPS4 style)
           if(b2_upwind .ne. 0) then
             fhe_mdf(ix,iy,0) =
     &       upwind(floe0_mdf,te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &       -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &       harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
           else
             fhe_mdf(ix,iy,0) =
     &       floe0_mdf*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.d0
     &       -abs(floe0_mdf/2.0)*tex
     &       -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &       harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
           endif
            floe(ix,iy,0) = floe0_mdf
            cone(ix,iy,0,0) = cone0_mdf/(1.0_R8+t0_mdf) !srv 02.01.07
          endif                                                          !srv 17.06.09 }
 
          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then                !srv 16.09.02 {
           w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
           w2(ix,iy)=(chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*          !srv 13.09.02 11.01.13 *signmf
     &     (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
           w3(ix,iy)=-abs(floe0/2)*tex-cone0*tex0
           if(alfTeEh.eq.1.0_R8) 
     &      w4(ix,iy)=
     *       alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
           if(fch_pTe.eq.1.0_R8)
     &      w4(ix,iy)=-0.71e0_R8*fch_p(ix,iy,0)/qe*
     *       (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
           endif                                                         !srv 16.09.02 }
           fch_p(ix,iy,0) = (csig(ix,iy,0)*ehx-calf(ix,iy,0)*tex0)       !srv 30.12.09 }
           if(pot_eq.ne.1 .and. no_current.ne.0) then                    !srv 27.02.13
            fch_p = 0.0e0_R8
            ehx = 0.0_R8                                                 !srv 27.02.13
           endif                                                         !srv 30.12.09 }

         else !lim_flux                                                          !srv 04.02.03 {
!           if(mdf_fhe.eq.0) then                                       !srv 13.10.06 {
          if(no_Ptncr_x_he.eq.0) then                                 !srv 06.07.06 {
            cone0     = hybr(floe0,chce(ix,iy,0))
            cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0))

csw 13mar2013 overwrite cone0/cone0_mdf (SOLPS4 style)
            if(b2tfhe_hybr2 .ne. 0) then
             cone0=dim(chce(ix,iy,0), abs(floe0)*
     .        max(0.0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .        max(1.e-30,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
             cone0_mdf=dim(chce(ix,iy,0), abs(floe0_mdf)*
     .        max(0.0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .        max(1.e-30,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
            endif

          else
             cone0     = chce(ix,iy,0)
             cone0_mdf = chce(ix,iy,0)
          endif                                                        !srv 06.07.06 }

!          if (chce(ix,iy,0).eq.0) then                   !iyv 12.05.08 !srv 03.10.08 {
!           alfe0 = 0.0_R8
!          else
!           alfe0 = calf(ix,iy,0)*cone0_mdf/chce(ix,iy,0)
!          endif
!          alfe0 = 0.0e0_R8 ! calf(ix,iy,0)               !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08  }

          alfe0 = calf(ix,iy,0)                          !srv 16.06.07
          tex0  = tex

csw 13mar2013 added check upwind (SOLPS4 style)
          if(b2_upwind .eq. 0) then
            fhe(ix,iy,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0) =                                          !srv 13.06.06 {
     &         floe0_mdf*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
     &         -cone0_mdf*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0)                                         !srv 13.06.06 }
          else
            fhe(ix,iy,0) =
     &         upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)) , te(ix,iy))
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0) =                                          !srv 13.06.06 {
     &         upwind(floe0_mdf,
     &                        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0)                                         !srv 13.06.06 }
          endif


          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then !srv 16.09.02 {
           if(mdf_fhe.eq.0) then                                        !srv 13.10.06 {
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
            w2(ix,iy)=(chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*          !srv 13.09.02 11.01.13 *signmf
     &      (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
            w3(ix,iy)= -cone0*tex0
            if(alfTeEh.eq.1.0_R8) 
     &       w4(ix,iy)=
     *        alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            if(fch_pTe.eq.1.0_R8)
     &       w4(ix,iy)=-0.71e0_R8*fch_p(ix,iy,0)/qe*
     *        (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
            w6(ix,iy) = cone0
           else
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
            w2(ix,iy)= chve(ix,iy,0)*                                    !srv 13.09.02
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &       *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2
            w3(ix,iy)= -cone0_mdf*tex0
            w4(ix,iy)=
     *       alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            w6(ix,iy) = cone0_mdf
           endif
          endif                                                !srv 16.09.02 }
          if(mdf_fhe.eq.0) then                                        !srv 13.10.06 {
           floe(ix,iy,0) = floe0
           cone(ix,iy,0,0) = cone0                          !srv 02.01.07
          else
           floe(ix,iy,0) = floe0_mdf
           cone(ix,iy,0,0) = cone0_mdf                          !srv 02.01.07
          endif                                                        !srv 13.10.06 }
         endif                                                   !srv 04.02.03 }

         fch(ix,iy,0) = fch_p(ix,iy,0)                            !srv 03.02.03
     &     +fchin(ix,iy,0)                                      !srv
     &     +fchdia(ix,iy,0)                                     !srv 22.10.98 21.04.00 05.02.05
     &     +fchvispar(ix,iy,0)                                  !srv 26.02.99
     &     +fchvisq(ix,iy,0)                                    !srv 11.07.99
     &     +fchinert(ix,iy,0)                                   !srv 07.07.99 21.04.00 05.02.05
     &     +fchanml(ix,iy,0)                                    !srv 24.03.07
!     &     +fchvisper(ix,iy,0)                                 !srv 22.11.99
!     &  /(bb(ix,iy,0)/bb(ix,iy,2)                              !srv 11.10.02 {
!     &   +bb(leftix(ix,iy),leftiy(ix,iy),0)/
!     &    bb(leftix(ix,iy),leftiy(ix,iy),2))                   !srv
!     &  /(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))         !srv
!     &  *(hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0       !srv 11.10.02 } 
         conc(ix,iy,0) = csig(ix,iy,0)
         if(fch_anomalous .ne. 0.0e0_R8 ) then                          !srv
          conc(ix,iy,0) =conc(ix,iy,0)  + csig_an(ix,iy,0)               !srv
         endif                                                            !srv
         if(fch_ion_neutral .ne. 0.0e0_R8 ) then                        !srv
          do is = 0, ns-1                                            !srv 24.06.08 {
           if (.not.is_neutral(is)) then
            do is0 = 0, ns-1
             conc(ix,iy,0) = conc(ix,iy,0) + csigin(ix,iy,0,is,is0) !*facdrift    !srv 07.02.03
     &        *(wbbl(ix,iy,2)/wbbl(ix,iy,3))**2                    !srv 15.02.11
!     1        *((bb(ix,iy,2)/bb(ix,iy,3))**2+                          !srv 25.02.03 {
!     2        (bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3        bb(leftix(ix,iy),leftiy(ix,iy),3))**2)*0.5e0_R8          !srv 25.02.03 }
            enddo
           endif
          enddo                                                      !srv 24.06.08 }
         endif                                                            !srv
        endif                                                         !srv 26.11.02
       enddo
      enddo

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then  !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32x')
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_strx')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTx')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehx')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'qex5')
       call my_out(70,nx,ny,w6,'b2tfhe__chcex_Pat')

       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;

      endif                                  !srv 16.09.02 }

*   ..compute fhey, floey, coney, fchy, concy
      fch_p(:,:,1) = 0.0e0_R8
      do ix = -1, nx
       floe(ix,-1,1)    = 0
       cone(ix,-1,1,0)    = 0                                     !srv 02.01.07
       fhe(ix,-1,1)     = 0
       fhe_mdf(ix,-1,1) = 0.0e0_R8                                !srv 13.10.06
       conc(ix,-1,1)    = 0
       fch(ix,-1,1)     = 0
      enddo
      do iy = 0, ny
       do ix = -1, nx
        if(pot_eq.eq.1 .and. no_current.eq.0) then                !srv 27.02.13
         ehy = -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &    (ne(ix,iy)*te(ix,iy)-ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &    te(bottomix(ix,iy),bottomiy(ix,iy)))*
     &    (2/qe)/(ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))
        else
         ehy=0.0_R8
        endif                                                        !srv 27.02.13

        floe0 = 1.5e0_R8*fne(ix,iy,1)+                                 !srv 11.05.99
!* in the next line there is the coeffisient +4 because                 !srv 13.08.02
!* we need this term with +5/2=-3/2+4                                   !srv 13.08.02 19.08.02
!     &   +4.0e0_R8*fledia(ix,iy,1)*vdia_par*signmf                  !srv 13.09.02
     &   (chve(ix,iy,1)+2.5e0_R8*fledia(ix,iy,1))*                  !srv 13.09.02 14.05.03 18.01.07 11.01.13 *signmf
     &   (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8   !srv 13.09.02
        floe0_mdf = 1.5e0_R8*fne(ix,iy,1)+                             !srv 11.05.99
     &   chve(ix,iy,1)*                                             !srv 13.09.02 14.05.03
     &   (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8   !srv 13.09.02
!        if(mdf_fhe.eq.0) then                                         !srv 13.10.06 {
        if(no_Ptncr_y_he.eq.0) then                                  !srv 06.07.06 {
         if(mdf_fhe.eq.0) then                                       !srv 13.10.06 {
          cone0     = hybr(floe0,chce(ix,iy,1))
         else
          if(facdriftm.eq.0.0e0_R8) then
           cone0    = hybr(floe0,chce(ix,iy,1))
          else
           cone0    = chce(ix,iy,1)
          endif
         endif
         cone0_mdf = hybr(floe0_mdf,chce(ix,iy,1))
        else
          cone0     = chce(ix,iy,1)
          cone0_mdf = chce(ix,iy,1)
        endif 

csw 13mar2013 overwrite cone0 (SOLPS4 style)
        if(b2tfhe_hybr2 .ne. 0) then

csw 30apr2014 flo53 ONLY FOR TESTING (SOLPS4)
         if(flo53 .ne. 0) then
          floe0=floe0*5.d0/3.d0
         endif

         cone0=dim(chce(ix,iy,1), abs(floe0)*
     .   max(0.0,min(te(ix,iy),te(bottomix(ix,iy),bottomiy(ix,iy))))/
     .   max(cutlo,te(ix,iy)+te(bottomix(ix,iy), bottomiy(ix,iy))))
         cone0_mdf=dim(chce(ix,iy,1), abs(floe0_mdf)*
     .   max(0.0,min(te(ix,iy),te(bottomix(ix,iy),bottomiy(ix,iy))))/
     .   max(cutlo,te(ix,iy)+te(bottomix(ix,iy), bottomiy(ix,iy))))
        endif

csw 13mar2013 added check upwind (SOLPS4 style)
        if(b2_upwind .ne. 0) then
         fhe(ix,iy,1) = upwind(floe0,
     &         te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))
     &        -cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &        +alfTeEh*calf(ix,iy,1)*                                !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy  !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
         fhe_mdf(ix,iy,1) = upwind(floe0_mdf,                                 !srv 13.10.06 {
     &         te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))
     &       -cone0_mdf*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhePSch(ix,iy,1)                                         !srv 13.10.06 }
     &        +alfTeEh*calf(ix,iy,1)*                                !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy  !iyv 12.05.08 !srv 03.10.08
        else
         fhe(ix,iy,1) = floe0*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2
     &         -cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &        +alfTeEh*calf(ix,iy,1)*                                !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy  !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
         fhe_mdf(ix,iy,1) = floe0_mdf*                                 !srv 13.10.06 {
     &        (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2
     &      -cone0_mdf*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhePSch(ix,iy,1)                                         !srv 13.10.06 }
     &        +alfTeEh*calf(ix,iy,1)*                                !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy  !iyv 12.05.08 !srv 03.10.08
        endif


        if(mdf_fhe.eq.0) then                                         !srv 13.10.06 {
         floe(ix,iy,1) = floe0
         cone(ix,iy,1,0) = cone0                           !srv 02.01.07            !srv 23.06.99
        else
         floe(ix,iy,1) = floe0_mdf
         cone(ix,iy,1,0) = cone0_mdf                           !srv 02.01.07          !srv 23.06.99
        endif                                                          !srv 13.10.06 }
        fch(ix,iy,1) =                                           !srv                             
     1      +fchin(ix,iy,1)                                      !srv
     2      +fchdia(ix,iy,1)
     3      +fchvispar(ix,iy,1)                                  !srv 26.02.99
     4      +fchvisq(ix,iy,1)                                    !srv 11.07.99
     5      +fchinert(ix,iy,1)                                   !srv 07.07.99
     6      +fchvisper(ix,iy,1)                                  !srv 22.11.99
     7      +fchanml(ix,iy,1)                                    !srv 24.03.07
     8      +fchstoch(ix,iy)                                     !lk 31.01.08 
        if(fch_anomalous.ne.0.0e0_R8) then  
          conc(ix,iy,1) = csig_an(ix,iy,1)                       !srv
        endif                                                    !srv
        if(fch_ion_neutral .ne. 0.0e0_R8 ) then                  !srv
         do is = 0, ns-1                                             !srv 24.06.08 {
          if (.not.is_neutral(is)) then
           do is0 = 0, ns-1
            conc(ix,iy,1) = conc(ix,iy,1) + csigin(ix,iy,1,is,is0) !*facdrift !srv 27.03.99 07.02.03 31.01.08
           enddo
          endif
         enddo                                                       !srv 24.06.08 }
        endif                                                    !srv
        if(fch_stochastic.ne.0.0e0_R8 ) then                     !lk 31.01.08 { 
          conc(ix,iy,1) = conc(ix,iy,1) + csig_stoch(ix,iy)
        endif                                                    !lk 31.01.08 }
         
       enddo
      enddo

csw 18feb2014 remove fluxes between guard cells 
      do ix=-1,nx
       do iy=-1,ny
         if(iffboundx(ix,iy).ne.0) then
           fhe(ix,iy,0)=0.0_R8
           fhe_mdf(ix,iy,0)=0.0_R8
           fhePSch(ix,iy,0)=0.0_R8
           floe(ix,iy,0)=0.0_R8
           cone(ix,iy,0,:)=0.0_R8
         endif
          if(iffboundy(ix,iy).ne.0) then
           fhe(ix,iy,1)=0.0_R8
           fhe_mdf(ix,iy,1)=0.0_R8
           fhePSch(ix,iy,1)=0.0_R8
           floe(ix,iy,1)=0.0_R8
           cone(ix,iy,1,:)=0.0_R8
         endif
       enddo
      enddo


      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then  !srv 16.09.02 {
       if(mdf_fhe.eq.0) then                                         !srv 13.10.06 {
        w1(ix,iy)=1.5e0_R8*fne(ix,iy,1)*
     1   (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2
        w2(ix,iy)=2.5e0_R8*fledia(ix,iy,1)*
     &    (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1   *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2
        w3(ix,iy)=-cone0*
     1    (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
        w4(ix,iy)=calf(ix,iy,1)*
     &     harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
        w6(ix,iy) = cone0
       else
        w1(ix,iy)=1.5e0_R8*fne(ix,iy,1)*
     1   (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2
        w2(ix,iy)=0.0e0_R8
        w3(ix,iy)=-cone0_mdf*
     1    (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
        w4(ix,iy)=calf(ix,iy,1)*
     &     harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
        w6(ix,iy) = cone0_mdf
       endif
      endif                                                !srv 16.09.02 }



      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then  !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32y')        !srv 16.09.02 {
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_stry')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTy')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehy')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'qey5')
       call my_out(70,nx,ny,w6,'b2tfhe__chcey_Pat')
      endif                                   !srv 16.09.02 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then  !srv 16.09.02 {
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif                                   !srv 16.09.02 }


*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
      DO IX=NXISO1(IB), NXISO2(IB)+1
        DO IY=NYISO1(IB), NYISO2(IB)
          fhe(ix,iy,0)    =0.0e0_R8
          floe(ix,iy,0)   =0.0e0_R8
          cone(ix,iy,0,:) =0.0e0_R8                          !srv 02.01.07
          fch(ix,iy,0)    =0.0e0_R8
          conc(ix,iy,0)   =0.0e0_R8
          enddo
        enddo
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)+1
            fhe(ix,iy,1)   =0.0e0_R8
            floe(ix,iy,1)  =0.0e0_R8
            cone(ix,iy,1,:)=0.0e0_R8                          !srv 02.01.07
            fch(ix,iy,1)   =0.0e0_R8
            conc(ix,iy,1)  =0.0e0_R8
          end do
        enddo
      enddo

      if(iout.ne.0) then                                             !srv 18.06.08 {
       call my_out(70,nx,ny,fch(-1,-1,0),'fchx')
       call my_out(70,nx,ny,fch(-1,-1,1),'fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0),'fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'fchvispary')
       call my_out(70,nx,ny,fchvisq(-1,-1,0),'fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1),'fchvisqy')
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0),'fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1),'fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0),'fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1),'fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0),'fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1),'fchanmly')
       call my_out(70,nx,ny,fchstoch,'fchstoch')                 !lk 31.01.08
       call my_out(70,nx,ny,fne(-1,-1,0),'fnex')
       call my_out(70,nx,ny,fne(-1,-1,1),'fney')
       call my_out(70,nx,ny,conc(-1,-1,0),'concx')
       call my_out(70,nx,ny,conc(-1,-1,1),'concy')
       call my_out(70,nx,ny,floe(-1,-1,0),'floex')
       call my_out(70,nx,ny,floe(-1,-1,1),'floey')
       call my_out(70,nx,ny,cone(-1,-1,0,0),'conex0')
       call my_out(70,nx,ny,cone(-1,-1,0,1),'conex1')
       call my_out(70,nx,ny,cone(-1,-1,0,2),'conex2')
       call my_out(70,nx,ny,cone(-1,-1,1,0),'coney0')
       call my_out(70,nx,ny,cone(-1,-1,1,1),'coney1')
       call my_out(70,nx,ny,cone(-1,-1,1,2),'coney2')
       call my_out(70,nx,ny,avisper(-1,-1,0),'avisperx')
       call my_out(70,nx,ny,avisper(-1,-1,1),'avispery')
       call my_out(70,nx,ny,bvisper(-1,-1,0),'bvisperx')
       call my_out(70,nx,ny,bvisper(-1,-1,1),'bvispery')
       call my_out(70,nx,ny,cvisper(-1,-1,0),'cvisperx')
       call my_out(70,nx,ny,cvisper(-1,-1,1),'cvispery')
       call my_out(70,nx,ny,fch_p(-1,-1,0),'fch_px')
       call my_out(70,nx,ny,fch_p(-1,-1,1),'fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0),'csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1),'csig_any')
       do is = 0, ns-1
        write (charns,'(i2.2)') is 
        call my_out(70,nx,ny,na(-1,-1,is),'na'//charns)
        do is0=0,ns-1
         write (charns0,'(i2.2)') is0 
         call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'csiginx'//charns//
     &        '_'//charns0)
         call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'csiginy'//charns//
     &        '_'//charns0)
        enddo
       enddo
       call my_out(70,nx,ny,fhe(-1,-1,0),'fhex')
       call my_out(70,nx,ny,fhe(-1,-1,1),'fhey')
       call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'fhe_mdfx')
       call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0),'fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1),'fhePSchy')
      endif                                                          !srv 18.06.08 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then                !srv 11.09.09 {
       call my_out(70,nx,ny,fhe(-1,-1,0),'b2tfhe_fhe_no_mdfx')
       call my_out(70,nx,ny,fhe(-1,-1,1),'b2tfhe_fhe_no_mdfy')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'b2tfhe_fhe_mdfx')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'b2tfhe_fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0),'b2tfhe_fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1),'b2tfhe_fhePSchy')
      endif

      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
!       call my_out(70,nx,ny,fch(-1,-1,0),'b2tfhe_fchx')
!       call my_out(70,nx,ny,fch(-1,-1,1),'b2tfhe_fchy')
       call my_out(70,nx,ny,fchin(-1,-1,0),'b2tfhe_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'b2tfhe_fchiny')
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'b2tfhe_fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'b2tfhe_fchvispary')
!       call my_out(70,nx,ny,fchvisq(-1,-1,0),'b2tfhe_fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1),'b2tfhe_fchvisqy')
!       call my_out(70,nx,ny,fchvisper(-1,-1,0),'b2tfhe_fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'b2tfhe_fchvispery')
       call my_out(70,nx,ny,fchinert(-1,-1,0),'b2tfhe_fchinertx')
       call my_out(70,nx,ny,fchinert(-1,-1,1),'b2tfhe_fchinerty')
       call my_out(70,nx,ny,fchdia(-1,-1,0),'b2tfhe_fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1),'b2tfhe_fchdiay')
       call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfhe_fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfhe_fchanmly')
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfhe_fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfhe_fch_py')
       call my_out(70,nx,ny,csig_an(-1,-1,0),'b2tfhe_csig_anx')
       call my_out(70,nx,ny,csig_an(-1,-1,1),'b2tfhe_csig_any')
       call my_out(70,nx,ny,fchstoch,'b2tfhe_fchstoch')                 !lk 31.01.08
      endif                                                          !srv 11.09.09 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfhe_

      end
