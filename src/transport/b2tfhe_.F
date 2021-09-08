*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhe_ (nx, ny, ns, facdrift, fac_ExB,                 !srv 23.11.10
     1 vol, hx, te, po, ne, na, vedia, fne,                              !srv 11.07.99
     2 chce, chve, calf, chvemx, cddi, fch_p, fchstoch,                  !srv 17.05.99 11.07.99
     & floe, cone, fhe, fhe_mdf, fhePSch, balance)                       !djm Jan2017
      use b2mod_types
      use b2mod_math
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_geo
     & , only : hy1
      use b2mod_balance                                                  !djm Jan2017
     & , only : fhe_32, fhe_52, fhe_thermj, fhe_cond, fhe_dia, fhe_ecrb,
     &          fhe_strange, fhe_pschused, balance_netcdf
      use b2mod_transport_fun                                            !srv 13.01.17
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &   facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),
     &   vol(-1:nx,-1:ny),
     1   hx(-1:nx,-1:ny),
     &   te(-1:nx,-1:ny), po(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &   na(-1:nx,-1:ny,0:ns-1),
     &   vedia(-1:nx,-1:ny,0:1), fne(-1:nx,-1:ny,0:1),
     &   chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     &   calf(-1:nx,-1:ny,0:1), chvemx(-1:nx,-1:ny),
     &   cddi(-1:nx,-1:ny,0:1,0:ns-1),                                   !srv 13.10.06
     &   fch_p(-1:nx,-1:ny,0:1), fchstoch(-1:nx,-1:ny,0:1)               !lk 31.01.08 { 24.08.21
      logical balance                                                    !djm Jan2017
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   floe(-1:nx,-1:ny,0:1), cone(-1:nx,-1:ny,0:1,0:2),
     *   fhe(-1:nx,-1:ny,0:1), fhe_mdf(-1:nx,-1:ny,0:1),
     &   fhePSch(-1:nx,-1:ny,0:1)                                        !srv 13.10.06
*   ..local variables
      real (kind=R8) :: upfac, hybfac !djm Jan2017
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

*   ..local variables
      integer ncall, n2, ix, iy, ib, is
      integer, save :: iout = 0, iout_b2wdat = 0,                        !srv 16.09.02
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0
      integer, save :: lim_flux = 0, mdf_fhe = 0                         !srv 04.02.03 13.10.06
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0              !srv 06.07.06
      integer, save :: pot_eq = 1                                        !srv 26.07.06
      integer, save :: conductive_limit = 1                              !srv 16.06.08 17.01.12
      integer, save :: style = 2                                         !srv 13.01.17
      integer, save :: no_current= 0                                     !srv 27.02.13
csw 10oct2012
      integer, save :: b2_upwind = 0
      integer, save :: b2tfhe_hybr2 = 0, flo53 = 0

      integer, save :: style_b2mndt = 1
      character(3) :: charns
      real (kind=R8), save :: alfTeEh = 0.0_R8                           !srv 28.04.08 30.05.08
      real (kind=R8), save :: fch_pTe = 1.0_R8                           !srv 30.12.09
      real (kind=R8), save :: PSch=1.0e0_R8                              !srv 17.01.07
      real (kind=R8) ::
     *  ehx, ehy, floe0, flomx, cone0, alfe0, t0,
     &  tex, tex0, floe0_mdf, cone0_mdf, alfe0_mdf, t0_mdf, tex0_mdf     !srv 22.11.99 13.10.06 17.06.09
      real (kind=R8) ::
     *    fledia(-1:nx,-1:ny,0:1), vdia_par,                             !srv 22.11.99 25.01.02
     &    ehx_(-1:nx,-1:ny), ehy_(-1:nx,-1:ny),                          !srv 13.01.17
     &    zeff(-1:nx,-1:ny), c071(-1:nx,-1:ny)                           !srv 13.01.17
      real (kind=R8) :: facdriftm, fac_ExBm                              !srv 23.12.19
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)           !srv 16.09.02

      save ncall, vdia_par
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization
      data ncall/0/, vdia_par/1.0e0_R8/                                  !srv 17.05.99 22.03.02
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhe_')
*   ..set internal parameters on first call                              !srv
      if (ncall.eq.0) then                                               !srv
       call ipgetr ('b2tfhe_vdia_par', vdia_par)                         !srv 17.05.99
       call ipgeti ('b2tfhe_lim_flux', lim_flux)                         !srv 04.02.03
       call ipgeti ('b2tfhe_iout', iout)                                 !srv 16.09.02
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                  !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                  !srv 06.07.06
       call ipgeti ('b2news_poteq', pot_eq)                              !srv 26.07.06
       call ipgeti ('b2tfhe_mdf', mdf_fhe)                               !srv 13.10.06
       call ipgetr ('b2tfhe_PSch', PSch)                                 !srv 17.01.07
       call ipgeti ('b2mndt_style', style_b2mndt)
       call ipgetr ('b2tfhe_alfTeEh', alfTeEh)                           !srv 28.04.08
       call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                           !srv 30.12.09
       if(fch_pTe .eq. 1.0_R8) alfTeEh = 0.0_R8                          !srv 16.04.10
       call xertst (.not.(0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe),
     &        ' Do not set 0.0_R8.ne.alfTeEh.and.0.0_R8.ne.fch_pTe')
       if (style_b2mndt.ge.1) then
         call xertst (.not.(alfTeEh.eq.1.0_R8.and.fch_pTe.eq.0.0_R8),
     &    ' Cannot set alfTeEh to 1.0 and fch_pTe to 0.0 when using'//
     &    ' the SOLPS5.2 physics model!')
       end if
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2wdat_iout', iout_b2wdat)                          !iys 17.04.17
       call ipgeti ('b2trcl_conductive_limit', conductive_limit)         !srv 16.06.08 17.01.12
       call ipgeti ('b2tfhe_no_current',no_current)
       call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &  'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
       call ipgeti ('b2tfhe_upwind', b2_upwind)
       call ipgeti ('b2tfhe_hybr2', b2tfhe_hybr2)
       call ipgeti ('b2tfhe_flo53', flo53)
       call ipgeti ('b2sigp_style', style)                               !srv 13.01.17
       endif                                                             !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of fne
       call b2xvff (nx, ny, fne, 'fne')
*    ..test sign of ne, te, chce
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
*    ..test edge values of chce, chve, chvemx, calf, calf_an
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvfx (nx, ny, chvemx, 'chvemx')
       call b2xvff (nx, ny, calf, 'calf')
       call b2xvff (nx, ny, calf_an, 'calf_an')
      endif

      facdriftm = maxval(facdrift)
      fac_ExBm  = maxval(fac_ExB)                                        !srv 23.12.19

      fledia = 0.0_R8                                                    !srv 25.01.02 {
      if(facdriftm.ne.0.0e0_R8) then                                     !srv 18.01.07
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
!         is = ismain   ! 1    ! for is = 1 (za = 1)                     !srv 12.05.11
* ..compute flediax
         do iy = -1, ny
          do ix = 0, nx
           if(leftix(ix,iy).ne.-2) then                                  !srv 26.11.02
            fledia(ix,iy,0) =
     &        (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))              !srv 25.01.02 was 5/2
     &       /(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))                !srv 17.05.99 11.03.00
     &       *vedia(ix,iy,0)*vdia_par                                    !srv 17.05.99 11.03.00 21.04.00 25.01.02 18.01.07
!     1       *(bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 13.09.02 {
!     2       bb(leftix(ix,iy),leftiy(ix,iy),2)/
!     3       bb(leftix(ix,iy),leftiy(ix,iy),3))/2.0e0_R8                 !srv 13.09.02 }
           endif                                                         !srv 26.11.02
          enddo
         enddo
* ..compute flediay
         do iy = 0, ny
          do ix = -1, nx
           fledia(ix,iy,1) =
     &       (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))           !srv 17.05.99 11.03.00 25.01.02 was 5/2
     &      /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))           !srv 17.05.99 11.03.00
     &      *vedia(ix,iy,1)*vdia_par                                     !srv 17.05.99 11.03.00 18.01.07
!     1      *(bb(ix,iy,2)/bb(ix,iy,3)+                                  !srv 13.09.02 {
!     2      bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
!     3      bb(bottomix(ix,iy),bottomiy(ix,iy),3))/2.0e0_R8             !srv 13.09.02 }
          enddo
         enddo                                                           !srv 25.01.02 }
        endif
       enddo
      endif                                                              !srv 18.01.07

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                         .or. iout_b2wdat.eq.4) then               !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny),
     *    w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                              !srv 16.09.02 }

      call b2tepsch(nx, ny, ns, facdrift, cddi, te, ne,
     .  cone, fhePSch)                                                   !srv 16.10.17

      fhe_mdf = 0.0e0_R8                                                 !srv 13.10.06

* ..compute parallel current
      call b2xehx (nx, ny, hx, te, po, ne, ehx_)                         !srv 13.01.17
      if(style.eq.2) then
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        do iy = -1, ny
          do ix = -1, nx
            c071(ix,iy) = zeff(ix,iy)*fce2n(zeff(ix,iy))                 !srv 21.07.21
          end do
        end do
      else
        c071 = 0.71e0_R8                                                 !srv 13.01.17
      endif

* ..compute fluxes and linearization
*   ..compute fhex, floex, conex
      do iy = -1, ny
       do ix = -1, nx                                                    !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then                                     !srv 26.11.02
         floe(ix,iy,0)    = 0.0e0_R8
         cone(ix,iy,0,0)  = 0.0e0_R8                                     !srv 02.01.07
         fhe(ix,iy,0)     = 0.0e0_R8
         fhe_mdf(ix,iy,0) = 0.0e0_R8                                     !srv 13.10.06
         if (balance.and.(balance_netcdf.ne.0)) then
           fhe_32(ix,iy,0) = 0.0_R8                                      !djm Jan2017
           fhe_52(ix,iy,0) = 0.0_R8
           fhe_thermj(ix,iy,0) = 0.0_R8
           fhe_cond(ix,iy,0) = 0.0_R8
           fhe_dia(ix,iy,0) = 0.0_R8
           fhe_ecrb(ix,iy,0) = 0.0_R8
           fhe_strange(ix,iy,0) = 0.0_R8
           fhe_pschused(ix,iy,0) = 0.0_R8
         endif
        else                                                             !srv 26.11.02 }
         ehx = ehx_(ix,iy)                                               !srv 13.01.17
         tex = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
         if(pot_eq.ne.1 .and. no_current.ne.0) then                      !srv 27.02.13
          ehx = 0.0_R8                                                   !srv 27.02.13
         endif                                                           !srv 30.12.09 }
         floe0 = 1.5e0_R8*fne(ix,iy,0)+
     &    (chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*                      !srv 13.09.02 18.01.07 11.01.13  *signmf
     &    (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8           !srv 13.09.02
     &    -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe                         !srv 30.12.09 13.01.17
         floe0_mdf = 1.5e0_R8*fne(ix,iy,0)+                              !srv 13.06.06 {
     &    chve(ix,iy,0)*
     &    (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8           !srv 13.06.06 }
     &    -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe                         !srv 30.12.09 13.01.17

        if(lim_flux.eq.1 .and. conductive_limit.eq.0) then               !srv 04.02.03 17.01.12

          if(no_Ptncr_x_he.eq.0) then
           cone0 = hybr(floe0,chce(ix,iy,0))-abs(floe0/2.0_R8)
          else
           cone0 = chce(ix,iy,0)-abs(floe0/2.0_R8)
          endif

          if (chce(ix,iy,0).eq.0) then
           alfe0 = 0.0_R8
          else
           alfe0 = calf(ix,iy,0)*cone0/chce(ix,iy,0)
          endif

cwdk 10.03.16 unified treatment of chvemx in b2tlhe.F. Test for b2tqce_model no longer required
          flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &      ne(ix,iy))/2
          if(flomx.gt.0.0_R8) then
            t0 = cone0*abs(tex)/(flomx*max(
     &        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
          else
            t0=0.0_R8
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
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.eq.0)) then !djm Jan2017
               if(floe0.ge.0.0_R8) then
                 upfac=1.0_R8
               else
                 upfac=0.0_R8
               endif
               fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &                           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*te(ix,iy))
               fhe_52(ix,iy,0) = 0.0_R8
               fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_dia(ix,iy,0) = 2.5e0_R8*fledia(ix,iy,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_thermj(ix,iy,0) =
     &              -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe*
     &              (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &               (1.0_R8-upfac)*te(ix,iy))+
     &              alfTeEh*alfe0*
     &              harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0) = 0.0_R8
               fhe_cond(ix,iy,0) = -cone0*
     &                      (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
               fhe_pschused(ix,iy,0) = 0.0_R8
             endif
          else
             fhe(ix,iy,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.eq.0)) then !djm Jan2017
               fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_52(ix,iy,0) = 0.0_R8
               fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_dia(ix,iy,0) = 2.5e0_R8*fledia(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_thermj(ix,iy,0) = -c071(ix,iy)*fch_p(ix,iy,0)/qe*
     &            fch_pTe*
     &            (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &            alfTeEh*alfe0*
     &            harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0) = 0.0_R8
               fhe_cond(ix,iy,0) = -cone0*tex0
               fhe_pschused(ix,iy,0) = 0.0_R8
             endif
          endif

          if(mdf_fhe.eq.0) then                                         !srv 17.06.09 {

csw 30apr2014
            if(b2tfhe_hybr2 .ne.0) then
              floe(ix,iy,0) = floe0
              cone(ix,iy,0,0) = cone0
            else
              floe(ix,iy,0) = floe0
              cone(ix,iy,0,0) = abs(floe0/2.0_R8)+cone0/(1.0_R8+t0)     !srv 02.01.07
            endif

          else

csw 17feb2012
            if(no_Ptncr_x_he.eq.0) then
              cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0))-
     &                    abs(floe0_mdf/2.0_R8)
            else
              cone0_mdf = chce(ix,iy,0)-abs(floe0_mdf/2.0_R8)
            endif

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0_mdf)
            if(b2tfhe_hybr2 .ne. 0) then
              cone0_mdf=dim(chce(ix,iy,0), abs(floe0_mdf)*
     .         max(0._R8,
     .         min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(1.e-30_R8,te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))))
            endif

            if (chce(ix,iy,0).eq.0) then
              alfe0_mdf = 0.0_R8
            else
              alfe0_mdf = calf(ix,iy,0)*cone0_mdf/chce(ix,iy,0)
            endif

            t0_mdf = cone0_mdf*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            tex0_mdf = tex/(1.0_R8+t0_mdf)

csw 27nov2012 added check for upwind (SOLPS4 style)
            if(b2_upwind .ne. 0) then
              fhe_mdf(ix,iy,0) =
     &         upwind(floe0_mdf,
     &                te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.ne.0)) then !djm Jan2017
                if(floe0_mdf.ge.0.0_R8) then
                  upfac=1.0_R8
                else
                  upfac=0.0_R8
                endif
                fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*te(ix,iy))
                fhe_52(ix,iy,0) = 0.0_R8
                fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
                fhe_dia(ix,iy,0) = 0.0_R8
                fhe_thermj(ix,iy,0) =
     &             -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))+
     &             alfTeEh*alfe0_mdf*
     &             harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0) = 0.0_R8
                fhe_cond(ix,iy,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0) = 0.0_R8
              endif
            else
              fhe_mdf(ix,iy,0) = floe0_mdf*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -abs(floe0_mdf/2.0_R8)*tex
     &         -cone0_mdf*tex0_mdf + alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhe.ne.0)) then !djm Jan2017
               !djm Jan2017 share the hybrid part according to relative contributions
                if (floe0_mdf.ge.0.0_R8) then
                  hybfac=1.0_R8
                else
                  hybfac=-1.0_R8
                endif
                fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &            hybfac*1.5e0_R8*fne(ix,iy,0)/2.0_R8*tex
                fhe_52(ix,iy,0) = 0.0_R8
                fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &           hybfac*chve(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8/
     &           2.0_R8*tex
                fhe_dia(ix,iy,0) = 0.0_R8
                fhe_thermj(ix,iy,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &           hybfac*c071(ix,iy)*fch_p(ix,iy,0)/qe*
     &           fch_pTe/2.0_R8*tex+
     &           alfTeEh*alfe0_mdf*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0) = 0.0_R8
                fhe_cond(ix,iy,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0) = 0.0_R8
              endif
            endif
            floe(ix,iy,0) = floe0_mdf
            cone(ix,iy,0,0) = cone0_mdf/(1.0_R8+t0_mdf)                  !srv 02.01.07
          endif                                                          !srv 17.06.09 }

          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                            .or. iout_b2wdat.eq.4) then            !srv 16.09.02 {
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)=(chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*          !srv 13.09.02 11.01.13 *signmf
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)=-abs(floe0/2)*tex-cone0*tex0
            if(alfTeEh.eq.1.0_R8)
     &       w4(ix,iy)=
     *        alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            if(fch_pTe.eq.1.0_R8)
     &       w4(ix,iy)=-c071(ix,iy)*fch_p(ix,iy,0)/qe*                   !srv 13.01.17
     *        (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w5(ix,iy) = 2.5e0_R8*fledia(ix,iy,0)*                        !srv 13.09.02 11.01.13 *signmf
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*       !srv 13.09.02
     &       (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
          endif                                                          !srv 16.09.02 }
          if(pot_eq.ne.1 .and. no_current.ne.0) then                     !srv 27.02.13
            ehx = 0.0_R8                                                 !srv 27.02.13
          endif                                                          !srv 30.12.09 }

        else !lim_flux                                                   !srv 04.02.03 {

          if(no_Ptncr_x_he.eq.0) then                                    !srv 06.07.06 {
            cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0))
            cone0     = cone0_mdf
csw 13mar2013 overwrite cone0/cone0_mdf (SOLPS4 style)
            if(b2tfhe_hybr2 .ne. 0) then
              cone0=dim(chce(ix,iy,0), abs(floe0)*
     .         max(0._R8,
     .         min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(1.e-30_R8,te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))))
              cone0_mdf=dim(chce(ix,iy,0), abs(floe0_mdf)*
     .         max(0._R8,
     .         min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(1.e-30_r8,te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))))
            endif

          else
            cone0     = chce(ix,iy,0)
            cone0_mdf = chce(ix,iy,0)
          endif                                                          !srv 06.07.06 }

!          if (chce(ix,iy,0).eq.0) then                    !iyv 12.05.08 !srv 03.10.08 {
!           alfe0 = 0.0_R8
!          else
!           alfe0 = calf(ix,iy,0)*cone0_mdf/chce(ix,iy,0)
!          endif
!          alfe0 = 0.0e0_R8 ! calf(ix,iy,0)  !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08 }

          alfe0 = calf(ix,iy,0)                                          !srv 16.06.07
          tex0  = tex

csw 13mar2013 added check upwind (SOLPS4 style)
          if(b2_upwind .eq. 0) then
            fhe(ix,iy,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0) = floe0_mdf*                                !srv 13.06.06 {
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0_mdf*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0)                                          !srv 13.06.06 }
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
              fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &               (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              fhe_52(ix,iy,0) = 0.0_R8
              fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &         (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              if (mdf_fhe.eq.0) then
                fhe_dia(ix,iy,0) = 2.5e0_R8*fledia(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
                fhe_cond(ix,iy,0) = -cone0*tex0
                fhe_pschused(ix,iy,0) = 0.0_R8
              else
                fhe_dia(ix,iy,0) = 0.0_R8
                fhe_cond(ix,iy,0) = -cone0_mdf*tex0
                fhe_pschused(ix,iy,0) = -fhePSch(ix,iy,0)
              endif
              fhe_thermj(ix,iy,0) = -c071(ix,iy)*fch_p(ix,iy,0)/qe*
     &          fch_pTe*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
     &          /2.0_R8 + alfTeEh*alfe0*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0) = 0.0_R8
            endif
          else
            fhe(ix,iy,0) =
     &         upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)) , te(ix,iy))
     &         -cone0*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0) =                                           !srv 13.06.06 {
     &         upwind(floe0_mdf,
     &                        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0 + alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0)                                          !srv 13.06.06 }
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhe.eq.0)) then !djm Jan2017
              if(floe0.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0) = 0.0_R8
              fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0) = 2.5e0_R8*fledia(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_thermj(ix,iy,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0) = 0.0_R8
              fhe_cond(ix,iy,0) = -cone0*tex0
              fhe_pschused(ix,iy,0) = 0.0_R8
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhe.ne.0)) then !djm Jan2017
              if(floe0_mdf.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0) = 1.5e0_R8*fne(ix,iy,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0) = 0.0_R8
              fhe_strange(ix,iy,0) = chve(ix,iy,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0) = 0.0_R8
              fhe_thermj(ix,iy,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0) = 0.0_R8
              fhe_cond(ix,iy,0) = -cone0_mdf*tex0
              fhe_pschused(ix,iy,0) = -fhePSch(ix,iy,0)
            endif
          endif


          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                             .or. iout_b2wdat.eq.4) then           !srv 16.09.02 {
           if(mdf_fhe.eq.0) then                                         !srv 13.10.06 {
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)=(chve(ix,iy,0)+2.5e0_R8*fledia(ix,iy,0))*          !srv 13.09.02 11.01.13 *signmf
     &      (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0*tex0
            w6(ix,iy) = cone0
           else
            w1(ix,iy)=1.5e0_R8*fne(ix,iy,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)= chve(ix,iy,0)*                                    !srv 13.09.02
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &       *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0_mdf*tex0
            w6(ix,iy) = cone0_mdf
           endif
           if(alfTeEh.eq.1.0_R8)
     &       w4(ix,iy)=
     *        alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
           if(fch_pTe.eq.1.0_R8)
     &       w4(ix,iy)=-c071(ix,iy)*fch_p(ix,iy,0)/qe*                   !srv 13.01.17
     *        (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
           w5(ix,iy) = 2.5e0_R8*fledia(ix,iy,0)*                         !srv 13.09.02 11.01.13 *signmf
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
          endif                                                          !srv 16.09.02 }
          if(mdf_fhe.eq.0) then                                          !srv 13.10.06 {
           floe(ix,iy,0) = floe0
           cone(ix,iy,0,0) = cone0                                       !srv 02.01.07
          else
           floe(ix,iy,0) = floe0_mdf
           cone(ix,iy,0,0) = cone0_mdf                                   !srv 02.01.07
          endif                                                          !srv 13.10.06 }
         endif                                                           !srv 04.02.03 }

        endif                                                            !srv 26.11.02
       enddo
      enddo

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                      .or. iout_b2wdat.eq.4) then                  !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32x')
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_strx')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTx')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehx')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'b2tfhe__qex5')
       call my_out(70,nx,ny,w5,'b2tfhe__qediax')
       call my_out(70,nx,ny,w6,'b2tfhe__chcex_Pat')

       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;

      endif                                                              !srv 16.09.02 }

*   ..compute fhey, floey, coney
      call b2xehy (nx, ny, te, po, ne, ehy_)                             !srv 13.01.17
      do ix = -1, nx
       floe(ix,-1,1)    = 0.0_R8
       cone(ix,-1,1,0)  = 0.0_R8                                         !srv 02.01.07
       fhe(ix,-1,1)     = 0.0_R8
       fhe_mdf(ix,-1,1) = 0.0_R8                                         !srv 13.10.06
       if (balance.and.(balance_netcdf.ne.0)) then
         fhe_32(ix,-1,1) = 0.0_R8                                        !djm Jan2017
         fhe_52(ix,-1,1) = 0.0_R8
         fhe_thermj(ix,-1,1) = 0.0_R8
         fhe_cond(ix,-1,1) = 0.0_R8
         fhe_dia(ix,-1,1) = 0.0_R8
         fhe_ecrb(ix,-1,1) = 0.0_R8
         fhe_strange(ix,-1,1) = 0.0_R8
         fhe_pschused(ix,-1,1) = 0.0_R8
       endif
      enddo
      do iy = 0, ny
       do ix = -1, nx
        if(pot_eq.eq.1 .and. no_current.eq.0) then                       !srv 27.02.13
         ehy = ehy_(ix,iy)                                               !srv 13.01.17
        else
         ehy = 0.0_R8
        endif                                                            !srv 27.02.13
        floe0 = 1.5e0_R8*(fne(ix,iy,1)-                                  !srv 11.05.99
     &     5.0e0_R8/3.0e0_R8*fchstoch(ix,iy,1)/qe)+                      !lk 21.02.08 17.12.13 24.08.21
!* in the next line there is the coeffisient +4 because                  !srv 13.08.02
!* we need this term with +5/2=-3/2+4                                    !srv 13.08.02 19.08.02
!     &   +4.0e0_R8*fledia(ix,iy,1)*vdia_par*signmf                      !srv 13.09.02
     &   (chve(ix,iy,1)+2.5e0_R8*fledia(ix,iy,1))*                       !srv 13.09.02 14.05.03 18.01.07 11.01.13 *signmf
     &   (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
        floe0_mdf = 1.5e0_R8*(fne(ix,iy,1)-                              !srv 11.05.99
     &    5.0e0_R8/3.0e0_R8*fchstoch(ix,iy,1)/qe)+                       !lk 21.02.08 17.12.13 24.08.21
     &   chve(ix,iy,1)*                                                  !srv 13.09.02 14.05.03
     &   (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
        if(no_Ptncr_y_he.eq.0) then                                      !srv 06.07.06 {
          cone0_mdf = hybr(floe0_mdf,chce(ix,iy,1))
          cone0     = cone0_mdf
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
     .   max(0._R8,min(te(ix,iy),te(bottomix(ix,iy),bottomiy(ix,iy))))/
     .   max(cutlo,te(ix,iy)+te(bottomix(ix,iy), bottomiy(ix,iy))))
         cone0_mdf=dim(chce(ix,iy,1), abs(floe0_mdf)*
     .   max(0._R8,min(te(ix,iy),te(bottomix(ix,iy),bottomiy(ix,iy))))/
     .   max(cutlo,te(ix,iy)+te(bottomix(ix,iy), bottomiy(ix,iy))))
        endif

        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                        .or. iout_b2wdat.eq.4) then                !srv 16.09.02 {
         if(mdf_fhe.eq.0) then                                           !srv 13.10.06 {
          w1(ix,iy)=1.5e0_R8*fne(ix,iy,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=2.5e0_R8*fledia(ix,iy,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w3(ix,iy)=-cone0*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0
         else
          w1(ix,iy)=1.5e0_R8*fne(ix,iy,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=0.0e0_R8
          w3(ix,iy)=-cone0_mdf*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0_mdf
         endif
         w5(ix,iy)=2.5e0_R8*fledia(ix,iy,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
        endif                                                            !srv 16.09.02 }

csw 13mar2013 added check upwind (SOLPS4 style)
        if(b2_upwind .ne. 0) then
         fhe(ix,iy,1) = upwind(floe0,
     &         te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))
     &        -cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &        +alfTeEh*calf(ix,iy,1)*                                    !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy   !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
         fhe_mdf(ix,iy,1) = upwind(floe0_mdf,                            !srv 13.10.06 {
     &         te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))
     &       -cone0_mdf*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhePSch(ix,iy,1)                                            !srv 13.10.06 }
     &        +alfTeEh*calf(ix,iy,1)*                                    !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy   !iyv 12.05.08 !srv 03.10.08
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhe.eq.0)) then !djm Jan2017
           if(floe0.ge.0.0_R8) then
             upfac=1.0_R8
           else
             upfac=0.0_R8
           endif
           fhe_32(ix,iy,1) = 1.5e0_R8*fne(ix,iy,1)*
     &                       (upfac*te(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*te(ix,iy))
           fhe_52(ix,iy,1) = 0.0_R8
           fhe_strange(ix,iy,1) = chve(ix,iy,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (upfac*te(bottomix(ix,iy),bottomiy(ix,iy))+
     &          (1.0_R8-upfac)*te(ix,iy))
           fhe_dia(ix,iy,1) = 2.5e0_R8*fledia(ix,iy,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (upfac*te(bottomix(ix,iy),bottomiy(ix,iy))+
     &          (1.0_R8-upfac)*te(ix,iy))
           fhe_cond(ix,iy,1) = -cone0*
     &                 (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
           fhe_pschused(ix,iy,1) = 0.0_R8
           fhe_thermj(ix,iy,1) = alfTeEh*calf(ix,iy,1)*
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
           fhe_ecrb(ix,iy,1) = 0.0_R8
         endif
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhe.ne.0)) then
           if(floe0_mdf.ge.0.0_R8) then
             upfac=1.0_R8
           else
             upfac=0.0_R8
           endif
           fhe_32(ix,iy,1) = 1.5e0_R8*fne(ix,iy,1)*
     &                       (upfac*te(bottomix(ix,iy),bottomiy(ix,iy))+
     &                        (1.0_R8-upfac)*te(ix,iy))
           fhe_52(ix,iy,1) = 0.0_R8
           fhe_strange(ix,iy,1) = chve(ix,iy,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (upfac*te(bottomix(ix,iy),bottomiy(ix,iy))+
     &          (1.0_R8-upfac)*te(ix,iy))
           fhe_dia(ix,iy,1) = 0.0_R8
           fhe_cond(ix,iy,1) = -cone0_mdf*
     &                 (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
           fhe_pschused(ix,iy,1) = -fhePSch(ix,iy,1)
           fhe_thermj(ix,iy,1) = alfTeEh*calf(ix,iy,1)*
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
           fhe_ecrb(ix,iy,1) = 0.0_R8
         endif
        else
         fhe(ix,iy,1) = floe0*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &        +alfTeEh*calf(ix,iy,1)*                                    !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy   !srv 03.06.03 !iyv 12.05.08 !srv 03.10.08
         fhe_mdf(ix,iy,1) = floe0_mdf*                                   !srv 13.10.06 {
     &        (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
     &      -cone0_mdf*(te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
     &      -fhePSch(ix,iy,1)                                            !srv 13.10.06 }
     &        +alfTeEh*calf(ix,iy,1)*                                    !iyv 12.05.08 !srv 03.10.08
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy   !iyv 12.05.08 !srv 03.10.08
         if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
           fhe_32(ix,iy,1) = 1.5e0_R8*fne(ix,iy,1)*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           fhe_52(ix,iy,1) = 0.0_R8
           fhe_strange(ix,iy,1) = chve(ix,iy,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           if (mdf_fhe.eq.0) then
             fhe_cond(ix,iy,1) = -cone0*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1) = 2.5e0_R8*fledia(ix,iy,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          else
             fhe_cond(ix,iy,1) = -cone0_mdf*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1) = 0.0_R8
           endif
           fhe_pschused(ix,iy,1) = -fhePSch(ix,iy,1)
           fhe_thermj(ix,iy,1) = alfTeEh*calf(ix,iy,1)*
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
           fhe_ecrb(ix,iy,1) = 0.0_R8
         endif
        endif

        if(mdf_fhe.eq.0) then                                            !srv 13.10.06 {
         floe(ix,iy,1)   = floe0
         cone(ix,iy,1,0) = cone0                           !srv 02.01.07 !srv 23.06.99
        else
         floe(ix,iy,1)   = floe0_mdf
         cone(ix,iy,1,0) = cone0_mdf                       !srv 02.01.07 !srv 23.06.99
        endif                                                            !srv 13.10.06 }

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
           if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
             fhe_32(ix,iy,0) = 0.0_R8
             fhe_52(ix,iy,0) = 0.0_R8
             fhe_thermj(ix,iy,0) = 0.0_R8
             fhe_cond(ix,iy,0) = 0.0_R8
             fhe_dia(ix,iy,0) = 0.0_R8
             fhe_ecrb(ix,iy,0) = 0.0_R8
             fhe_strange(ix,iy,0) = 0.0_R8
             fhe_pschused(ix,iy,0) = 0.0_R8
           endif
         endif
          if(iffboundy(ix,iy).ne.0) then
           fhe(ix,iy,1)=0.0_R8
           fhe_mdf(ix,iy,1)=0.0_R8
           fhePSch(ix,iy,1)=0.0_R8
           floe(ix,iy,1)=0.0_R8
           cone(ix,iy,1,:)=0.0_R8
           if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
             fhe_32(ix,iy,1) = 0.0_R8
             fhe_52(ix,iy,1) = 0.0_R8
             fhe_thermj(ix,iy,1) = 0.0_R8
             fhe_cond(ix,iy,1) = 0.0_R8
             fhe_dia(ix,iy,1) = 0.0_R8
             fhe_ecrb(ix,iy,1) = 0.0_R8
             fhe_strange(ix,iy,1) = 0.0_R8
             fhe_pschused(ix,iy,1) = 0.0_R8
           endif
         endif
       enddo
      enddo

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then   !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32y')                    !srv 16.09.02 {
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_stry')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTy')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehy')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'b2tfhe__qey5')
       call my_out(70,nx,ny,w5,'b2tfhe__qediay')
       call my_out(70,nx,ny,w6,'b2tfhe__chcey_Pat')
      endif                                                         !srv 16.09.02 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then   !srv 16.09.02 {
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif                                                         !srv 16.09.02 }


*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fhe(ix,iy,0)    = 0.0e0_R8
            floe(ix,iy,0)   = 0.0e0_R8
            cone(ix,iy,0,:) = 0.0e0_R8                              !srv 02.01.07
            if (balance.and.(balance_netcdf.ne.0)) then
              fhe_32(ix,iy,0) = 0.0_R8                              !djm Jan2017
              fhe_52(ix,iy,0) = 0.0_R8
              fhe_thermj(ix,iy,0) = 0.0_R8
              fhe_cond(ix,iy,0) = 0.0_R8
              fhe_dia(ix,iy,0) = 0.0_R8
              fhe_ecrb(ix,iy,0) = 0.0_R8
              fhe_strange(ix,iy,0) = 0.0_R8
              fhe_pschused(ix,iy,0) = 0.0_R8
            endif
          enddo
        enddo
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fhe(ix,iy,1)   = 0.0e0_R8
            floe(ix,iy,1)  = 0.0e0_R8
            cone(ix,iy,1,:)= 0.0e0_R8                               !srv 02.01.07
            if (balance.and.(balance_netcdf.ne.0)) then
              fhe_32(ix,iy,1) = 0.0_R8                              !djm Jan2017
              fhe_52(ix,iy,1) = 0.0_R8
              fhe_thermj(ix,iy,1) = 0.0_R8
              fhe_cond(ix,iy,1) = 0.0_R8
              fhe_dia(ix,iy,1) = 0.0_R8
              fhe_ecrb(ix,iy,1) = 0.0_R8
              fhe_strange(ix,iy,1) = 0.0_R8
              fhe_pschused(ix,iy,1) = 0.0_R8
            endif
          end do
        enddo
      enddo

      if(iout.ne.0) then                                             !srv 18.06.08 {
       call my_out(70,nx,ny,fchstoch(-1,-1,1),'b2tfhe__fchstochy')   !lk 31.01.08 24.08.21
       call my_out(70,nx,ny,fne(-1,-1,0),'b2tfhe__fnex')
       call my_out(70,nx,ny,fne(-1,-1,1),'b2tfhe__fney')
       call my_out(70,nx,ny,floe(-1,-1,0),'b2tfhe__floex')
       call my_out(70,nx,ny,floe(-1,-1,1),'b2tfhe__floey')
       call my_out(70,nx,ny,cone(-1,-1,0,0),'b2tfhe__conex0')
       call my_out(70,nx,ny,cone(-1,-1,0,1),'b2tfhe__conex1')
       call my_out(70,nx,ny,cone(-1,-1,0,2),'b2tfhe__conex2')
       call my_out(70,nx,ny,cone(-1,-1,1,0),'b2tfhe__coney0')
       call my_out(70,nx,ny,cone(-1,-1,1,1),'b2tfhe__coney1')
       call my_out(70,nx,ny,cone(-1,-1,1,2),'b2tfhe__coney2')
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfhe__fch_px')
       call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfhe__fch_py')
       call my_out(70,nx,ny,chve(-1,-1,0),'b2tfhe__chvex')
       call my_out(70,nx,ny,chve(-1,-1,1),'b2tfhe__chvey')
       call my_out(70,nx,ny,calf(-1,-1,0),'b2tfhe__calfx')
       call my_out(70,nx,ny,calf(-1,-1,1),'b2tfhe__calfy')
       call my_out(70,nx,ny,chce(-1,-1,0),'b2tfhe__chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'b2tfhe__chcey')
       call my_out(70,nx,ny,te,'b2tfhe__te')
       call my_out(70,nx,ny,ne,'b2tfhe__ne')
       call my_out(70,nx,ny,po,'b2tfhe__po')
       do is = 0, ns-1
        write (charns,'(i3.3)') is
        call my_out(70,nx,ny,na(-1,-1,is),'b2tfhe__na'//charns)
       enddo
       call my_out(70,nx,ny,fhe(-1,-1,0),'b2tfhe__fhex')
       call my_out(70,nx,ny,fhe(-1,-1,1),'b2tfhe__fhey')
       call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'b2tfhe__fhe_mdfx')
       call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'b2tfhe__fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0),'b2tfhe__fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1),'b2tfhe__fhePSchy')
       call my_out(70,nx,ny,ehx_,'b2tfhe__ehx_')                     !srv 13.01.17
       call my_out(70,nx,ny,ehy_,'b2tfhe__ehy_')                     !srv 13.01.17
      endif                                                          !srv 18.06.08 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                          .or. iout_b2wdat.eq.4) then          !srv 11.09.09 {
       call my_out(70,nx,ny,fhe(-1,-1,0),'b2tfhe__fhe_no_mdfx')
       call my_out(70,nx,ny,fhe(-1,-1,1),'b2tfhe__fhe_no_mdfy')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'b2tfhe__fhe_mdfx')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'b2tfhe__fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0),'b2tfhe__fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1),'b2tfhe__fhePSchy')
      endif

      if((iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0)
     &                          .and. .not. balance) then
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfhe__for_po_fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfhe__fch_py')
       call my_out(70,nx,ny,fchstoch(-1,-1,1),
     &            'b2tfhe__for_po_fchstochy')                        !lk 31.01.08  24.08.21
      endif                                                          !srv 11.09.09 }

      if(iout_b2wdat.eq.4 .and. balance) then
       call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfhe__fch_px')
!      call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfhe__fch_py')
       call my_out(70,nx,ny,fchstoch(-1,-1,1),'b2tfhe__fchstochy')   !lk 31.01.08  24.08.21
      endif                                                          !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfhe_

      end subroutine b2tfhe_
