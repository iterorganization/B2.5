*-----------------------------------------------------------------------
*.specification

      subroutine b2tfhe_ (nCv, nFc, nVx, ns, switch, geo, mpg, 
     &  pl, dv, co, rt, st_ext, balance)                                 !djm Jan2017 !srv 29.01.20
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance                                                  !djm Jan2017
!WG_TODO     & , only : fhe_32, fhe_52, fhe_thermj, fhe_cond, fhe_dia, fhe_ecrb,
!WG_TODO     &          fhe_strange, fhe_pschused, balance_netcdf
      use b2mod_transport_fun                                            !srv 13.01.17
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tfhe, ncall_b2xehy, ncall_b2xehx, b2tfhe_cutlo
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl 
      type (B2Rates), intent (in) :: rt 
      type (B2StateExt), intent (in) :: st_ext
      logical balance                                                   !djm Jan2017
*   ..output arguments (unspecified on entry)
*   ..input/output arguments (unspecified on entry)
      type (B2Coeff), intent (inout) :: co
      type (B2Derivatives), intent (inout) :: dv
*   ..local variables
      real (kind=R8) :: upfac, hybfac, harm !djm Jan2017
      external harm
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
      integer iCv, is
csw 10oct2012
      character(3) :: charns, charns0
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nFc)                   !lk 31.01.08 {
      real (kind=R8) ::
     *  ehx(nFc), ehy(nFc), floe0(nFc,0:1), cone0(nFc,0:1), flomx,
     &  alfe0(nFc), floe0_mdf(nFc,0:1)                                   !srv 22.11.99 13.10.06 17.06.09
      real (kind=R8) ::
     *    fledia(nFc,0:1), c071(nFc), c071c(nCv),                        !srv 22.11.99 25.01.02
     &    zeff(nCv)                                                      !srv 13.01.17
      real (kind=R8) :: facdriftm, fac_ExBm                              !srv 23.12.19
      real (kind=R8) ::
     *    wrk(nFc,0:1), wrkv(nVx), tef(nFc), nef(nFc), tefh(nFc),
     *    dte(nFc,0:1), dpo(nFc,0:1), dpe(nFc,0:1), wght(nFc,2),
     *    fhe0_mdf(nFc,0:1), dumm0(nFc,0:1), dumm1(nFc,0:1)
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)           !srv 16.09.02
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c      integer, save :: b2_upwind = 0
c      integer, save :: b2tfhe_hybr2 = 0, flo53 = 0
c      integer, save :: lim_flux = 0                                      !srv 04.02.03 
c      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0              !srv 06.07.06
c      real (kind=R8), save :: drift_hyb = 1.0_R8
c      integer, save :: conductive_limit = 1                              !srv 16.06.08 17.01.12
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfhe_')
*   ..set internal parameters on first call                              !srv
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c       call ipgeti ('b2tfhe_upwind', b2_upwind)
c       call ipgeti ('b2tfhe_hybr2', b2tfhe_hybr2)
c       call ipgeti ('b2tfhe_flo53', flo53)
c       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                  !srv 06.07.06
c       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                  !srv 06.07.06
c       call ipgeti ('b2tfhe_lim_flux', lim_flux)                         !srv 04.02.03
c       call ipgetr ('b2tfhe_drift_hyb',drift_hyb)
c       call ipgeti ('b2trcl_conductive_limit', conductive_limit)         !srv 16.06.08 17.01.12
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2tfhe.lt.3) then
*    ..test sign of ne, te, chce
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nFc, co%chce(:,0)*geo%fcQalf(:,0),1,'chce0','.ge.')
       call b2xvsg (nFc, co%chce(:,1)*geo%fcQalf(:,1),1,'chce1','.ge.')
      endif

      wght = 1.0_R8

      facdriftm = maxval(dv%facdrift)
      fac_ExBm  = maxval(dv%fac_ExB)                                     !srv 23.12.19
      dv%fhePSch   = 0.0_R8
      fchstoch     = 0.0_R8
      csig_stoch   = 0.0_R8
      fledia       = 0.0_R8 
      if(facdriftm.ne.0.0e0_R8) then                                     !srv 18.01.07
       do is = 0,ns-1                                                    !srv 12.05.11
        if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then              !srv 12.05.11
*      ..compute fledia
         fledia(:,0) = switch%fhe_vdia_par*dv%vedia(:,0)*
     &                                     geo%fcS*geo%fcQalf(:,0)
         fledia(:,1) = switch%fhe_vdia_par*dv%vedia(:,1)*
     &                                     geo%fcS*geo%fcQalf(:,1)
        endif
       enddo
      endif                                                              !srv 18.01.07
#ifdef WG_TODO
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                         .or. iout_b2wdat.eq.4) then               !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny),
     *    w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                              !srv 16.09.02 }
#endif
      call b2tepsch(nCv, nFc, nVx, ns, switch, geo, mpg, dv%facdrift,
     &    co%cddi, pl%te, dv%ne, dv%floe, dv%cone, dv%fhePSch)         !srv 16.10.17

* ..compute parallel current
      call b2xehx (nCv, nFc, nVx, switch, geo, mpg, pl%te, 
     &                                         pl%po, dv%ne, ehx)      !srv 13.01.17 {
      call b2xehy (nCv, nFc, nVx, switch, geo, mpg, pl%te, 
     &                                         pl%po, dv%ne, ehy)      !srv 13.01.17
      if(switch%b2sigp_style.eq.2) then
        call b2xzef (nCv, ns, rt%rz2, pl%na, dv%ne, zeff, st_ext)
        do iCv = 1, nCv
          c071c(iCv) = zeff(iCv)*fce2n(zeff(iCv))
        end do
        call intface(nCv,nFc,mpg%fcCv,wght,c071c,c071)
      else
        c071 = 0.71e0_R8                                               !srv 13.01.17
      endif

*   ..interpolate ne and te to cell faces
      call intface(nCv,nFc,mpg%fcCv,wght,dv%ne,nef)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%te,tef)

*   ..compute differences of te, po, and pe
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dte)
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%po, wrkv, dpo)
      call diff(nCv, nFc, nVx, 0, geo, mpg, pl%te*dv%ne, wrkv, dpe)

* ..compute fluxes and linearization
       floe0(:,0) = 1.5e0_R8*dv%fne_he(:,0) +
     &    (co%chve(:,0) + 2.5e0_R8*fledia(:,0))*nef                    !srv 13.09.02 18.01.07 11.01.13  *signmf
     &    -c071*dv%fch_p(:,0)/qe*switch%fch_pTe                        !srv 30.12.09 13.01.17
       floe0(:,1) = 1.5e0_R8*dv%fne_he(:,1) +                           
     &    (co%chve(:,1) + 2.5e0_R8*fledia(:,1))*nef -                  !srv 11.05.99
     &    5.0e0_R8/3.0e0_R8*fchstoch/qe

       floe0_mdf(:,0) = 1.5e0_R8*dv%fne_he(:,0) +                      !srv 13.06.06 {
     &    co%chve(:,0)*nef - c071*dv%fch_p(:,0)/qe*
     &                                      switch%fch_pTe             !srv 30.12.09 13.01.17
       floe0_mdf(:,1) = 1.5e0_R8*dv%fne_he(:,1) +
     &    co%chve(:,1)*nef - 5.0e0_R8/3.0e0_R8*fchstoch/qe             !srv 11.05.99

*   ..apply discretization scheme
      call calcflow(nCv,nFc,nVx,switch%b2tfhe_discr_meth,
     &  geo,mpg,pl%te,floe0_mdf,
     &  co%chce,fhe0_mdf,dumm0,dumm1)

*   ..compute fluxes
      alfe0 = co%calf(:,0)
      tefh  = tef !WG_TODO: harmonic average to faces

      dv%fhe(:,0) = fhe0_mdf(:,0) +  (floe0(:,0) - floe0_mdf(:,0))*tef + 
     &  switch%alfTeEh*alfe0*tefh*ehx
      dv%fhe(:,1) = fhe0_mdf(:,1) +  (floe0(:,1) - floe0_mdf(:,1))*tef + 
     &  switch%alfTeEh*co%calf(:,1)*tefh*ehy

      dv%fhe_mdf(:,0)  = fhe0_mdf(:,0) - dv%fhePSch(:,0)  + 
     &  switch%alfTeEh*alfe0*tefh*ehx
      dv%fhe_mdf(:,1)  = fhe0_mdf(:,1) - dv%fhePSch(:,1)  + 
     &  switch%alfTeEh*co%calf(:,1)*tefh*ehy

*   ..compute floe, cone
      call calccoef(nCv,nFc,nVx,switch%b2tfhe_discr_meth,
     &  geo,mpg,floe0_mdf,co%chce,
     &  floe0,cone0)
      !wdk now we always add PSch linearization; do we need to condition
      !    this on mdf /= 0?
      dv%floe        = dv%floe + floe0
      dv%cone(:,:,0) = dv%cone(:,:,0) + cone0


#ifdef WG_TODO
            if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &               (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &         (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
              if (switch%mdf_fhe.eq.0) then
                fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
                fhe_cond(ix,iy,0,0) = -cone0*tex0
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              else
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0
                fhe_pschused(ix,iy,0,0) = -fhePSch(ix,iy,0,0)
              endif
              fhe_thermj(ix,iy,0,0) = -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &          fch_pTe*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
     &          /2.0_R8 + switch%alfTeEh*alfe0*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
            endif
          else
            fhe(ix,iy,0,0) =
     &         upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)) , te(ix,iy))
     &         -cone0*tex0 + switch%alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            fhe_mdf(ix,iy,0,0) =                                         !srv 13.06.06 {
     &         upwind(floe0_mdf,
     &                        te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0 + switch%alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
     &        -fhePSch(ix,iy,0,0)                                        !srv 13.06.06 }
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (switch%mdf_fhe.eq.0)) then !djm Jan2017
              if(floe0.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           switch%alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
              fhe_cond(ix,iy,0,0) = -cone0*tex0
              fhe_pschused(ix,iy,0,0) = 0.0_R8
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (switch%mdf_fhe.ne.0)) then !djm Jan2017
              if(floe0_mdf.ge.0.0_R8) then
                upfac=1.0_R8
              else
                upfac=0.0_R8
              endif
              fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*te(ix,iy))
              fhe_52(ix,iy,0,0) = 0.0_R8
              fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))
              fhe_dia(ix,iy,0,0) = 0.0_R8
              fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &            (1.0_R8-upfac)*te(ix,iy))+
     &           switch%alfTeEh*alfe0*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              fhe_ecrb(ix,iy,0,0) = 0.0_R8
              fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0
              fhe_pschused(ix,iy,0,0) = -fhePSch(ix,iy,0,0)
            endif
          endif


          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                             .or. iout_b2wdat.eq.4) then           !srv 16.09.02 {
           if(switch%mdf_fhe.eq.0) then                                         !srv 13.10.06 {
            w1(ix,iy)=1.5e0_R8*fne_he(ix,iy,0,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)=(chve(ix,iy,0,0)+2.5e0_R8*fledia(ix,iy,0,0))*      !srv 13.09.02 11.01.13 *signmf
     &      (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8         !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0*tex0
            if(switch%alfTeEh.eq.1.0_R8)
     &       w4(ix,iy)=
     *        alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            if(fch_pTe.eq.1.0_R8)
     &       w4(ix,iy)=-c071(ix,iy)*fch_p(ix,iy,0,0)/qe*                 !srv 13.01.17
     *        (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w6(ix,iy) = cone0
           else
            w1(ix,iy)=1.5e0_R8*fne_he(ix,iy,0,0)
     &           *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w2(ix,iy)= chve(ix,iy,0,0)*                                  !srv 13.09.02
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &       *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
            w3(ix,iy)= -cone0_mdf*tex0
            w4(ix,iy)=
     *       alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
            w6(ix,iy) = cone0_mdf
           endif
           w5(ix,iy) = 2.5e0_R8*fledia(ix,iy,0,0)*                       !srv 13.09.02 11.01.13 *signmf
     &       (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8        !srv 13.09.02
     &          *(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
          endif                                                          !srv 16.09.02 }

#endif 


#ifdef WG_TODO
        if(lim_flux.eq.1 .and. conductive_limit.eq.0) then               !srv 04.02.03 17.01.12

          if(no_Ptncr_x_he.eq.0) then
           cone0 = hybr(floe0,chce(ix,iy,0,0))-abs(floe0/2.0_R8)
          else
           cone0 = chce(ix,iy,0,0)-abs(floe0/2.0_R8)
          endif

          if (chce(ix,iy,0,0).eq.0) then
           alfe0 = 0.0_R8
          else
           alfe0 = calf(ix,iy,0,0)*cone0/chce(ix,iy,0,0)
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

           cone0=dim(chce(ix,iy,0,0)/(1.d0+t0), abs(floe0)*
     .      max(0.0d0,min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .      max(b2tfhe_cutlo,te(ix,iy)+te(leftix(ix,iy), leftiy(ix,iy))))
          endif


csw 27nov2012 add check for upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
             fhe(ix,iy,0,0) =
     &        upwind(floe0,te(leftix(ix,iy),leftiy(ix,iy)), te(ix,iy))
     &         -cone0*(te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
     &         + switch%alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (switch%mdf_fhe.eq.0)) then !djm Jan2017
               if(floe0.ge.0.0_R8) then
                 upfac=1.0_R8
               else
                 upfac=0.0_R8
               endif
               fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &                           (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*te(ix,iy))
               fhe_52(ix,iy,0,0) = 0.0_R8
               fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
               fhe_thermj(ix,iy,0,0) =
     &              -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &              (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &               (1.0_R8-upfac)*te(ix,iy))+
     &              switch%alfTeEh*alfe0*
     &              harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0,0) = 0.0_R8
               fhe_cond(ix,iy,0,0) = -cone0*
     &                      (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
               fhe_pschused(ix,iy,0,0) = 0.0_R8
             endif
          else
             fhe(ix,iy,0,0) =
     &         floe0*(te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -cone0*tex0 + switch%alfTeEh*
     &         alfe0*harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (switch%mdf_fhe.eq.0)) then !djm Jan2017
               fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_52(ix,iy,0,0) = 0.0_R8
               fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_dia(ix,iy,0,0) = 2.5e0_R8*fledia(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
               fhe_thermj(ix,iy,0,0) = -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &            fch_pTe*
     &            (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &            switch%alfTeEh*alfe0*
     &            harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
               fhe_ecrb(ix,iy,0,0) = 0.0_R8
               fhe_cond(ix,iy,0,0) = -cone0*tex0
               fhe_pschused(ix,iy,0,0) = 0.0_R8
             endif
          endif

          if(switch%mdf_fhe.eq.0) then                                         !srv 17.06.09 {

csw 30apr2014
            if(b2tfhe_hybr2 .ne.0) then
              floe(ix,iy,0,0) = floe0
              cone(ix,iy,0,0,0) = cone0
            else
              floe(ix,iy,0,0) = floe0
              cone(ix,iy,0,0,0) = abs(floe0/2.0_R8)+cone0/(1.0_R8+t0)   !srv 02.01.07
            endif

          else

csw 17feb2012
            if(no_Ptncr_x_he.eq.0) then
              cone0_mdf = hybr(floe0_mdf,chce(ix,iy,0,0))-
     &                    abs(floe0_mdf/2.0_R8)
            else
              cone0_mdf = chce(ix,iy,0,0)-abs(floe0_mdf/2.0_R8)
            endif

csw 13mar2013 SOLPS4 modified f. flux limited hybrid style (overwrite cone0_mdf)
            if(b2tfhe_hybr2 .ne. 0) then
              cone0_mdf=dim(chce(ix,iy,0,0), abs(floe0_mdf)*
     .         max(0._R8,
     .         min(te(ix,iy),te(leftix(ix,iy),leftiy(ix,iy))))/
     .         max(1.e-30_R8,te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))))
            endif

            if (chce(ix,iy,0,0).eq.0) then
              alfe0_mdf = 0.0_R8
            else
              alfe0_mdf = calf(ix,iy,0,0)*cone0_mdf/chce(ix,iy,0,0)
            endif

            t0_mdf = cone0_mdf*abs(tex)/(flomx*max(
     &         te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            tex0_mdf = tex/(1.0_R8+t0_mdf)

csw 27nov2012 added check for upwind (SOLPS4 style)
            if(b2_upwind .ne. 0) then
              fhe_mdf(ix,iy,0,0) =
     &         upwind(floe0_mdf,
     &                te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))
     &         -cone0_mdf*tex0_mdf + switch%alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (switch%mdf_fhe.ne.0)) then !djm Jan2017
                if(floe0_mdf.ge.0.0_R8) then
                  upfac=1.0_R8
                else
                  upfac=0.0_R8
                endif
                fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &                          (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*te(ix,iy))
                fhe_52(ix,iy,0,0) = 0.0_R8
                fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &             (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_thermj(ix,iy,0,0) =
     &             -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &             (upfac*te(leftix(ix,iy),leftiy(ix,iy))+
     &              (1.0_R8-upfac)*te(ix,iy))+
     &             switch%alfTeEh*alfe0_mdf*
     &             harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              endif
            else
              fhe_mdf(ix,iy,0,0) = floe0_mdf*
     &         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8
     &         -abs(floe0_mdf/2.0_R8)*tex
     &         -cone0_mdf*tex0_mdf + switch%alfTeEh*alfe0_mdf*
     &          harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
              if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (switch%mdf_fhe.ne.0)) then !djm Jan2017
               !djm Jan2017 share the hybrid part according to relative contributions
                if (floe0_mdf.ge.0.0_R8) then
                  hybfac=1.0_R8
                else
                  hybfac=-1.0_R8
                endif
                fhe_32(ix,iy,0,0) = 1.5e0_R8*fne_he(ix,iy,0,0)*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &            hybfac*1.5e0_R8*fne_he(ix,iy,0,0)/2.0_R8*tex
                fhe_52(ix,iy,0,0) = 0.0_R8
                fhe_strange(ix,iy,0,0) = chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8-
     &           hybfac*chve(ix,iy,0,0)*
     &           (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*0.5_R8/
     &           2.0_R8*tex
                fhe_dia(ix,iy,0,0) = 0.0_R8
                fhe_thermj(ix,iy,0,0) =
     &           -c071(ix,iy)*fch_p(ix,iy,0,0)/qe*fch_pTe*
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8+
     &           hybfac*c071(ix,iy)*fch_p(ix,iy,0,0)/qe*
     &           fch_pTe/2.0_R8*tex+
     &           switch%alfTeEh*alfe0_mdf*
     &           harm(te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy))*ehx
                fhe_ecrb(ix,iy,0,0) = 0.0_R8
                fhe_cond(ix,iy,0,0) = -cone0_mdf*tex0_mdf
                fhe_pschused(ix,iy,0,0) = 0.0_R8
              endif
            endif
            floe(ix,iy,0,0) = floe0_mdf
            cone(ix,iy,0,0,0) = cone0_mdf/(1.0_R8+t0_mdf)                !srv 02.01.07
          endif                                                          !srv 17.06.09 }

        else !lim_flux 
#endif 

#ifdef WG_TODO
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
#endif
#ifdef WG_TODO

        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                        .or. iout_b2wdat.eq.4) then                !srv 16.09.02 {
         if(mdf_fhe.eq.0) then                                           !srv 13.10.06 {
          w1(ix,iy)=1.5e0_R8*fne_he(ix,iy,1,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=2.5e0_R8*fledia(ix,iy,1,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w3(ix,iy)=-cone0*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0
         else
          w1(ix,iy)=1.5e0_R8*fne_he(ix,iy,1,1)*
     1      (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          w2(ix,iy)=0.0e0_R8
          w3(ix,iy)=-cone0_mdf*
     1      (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=calf(ix,iy,1,1)*
     &       harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
          w6(ix,iy) = cone0_mdf
         endif
         w5(ix,iy)=2.5e0_R8*fledia(ix,iy,1,1)*
     &      (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8
     1     *(te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
        endif                                                            !srv 16.09.02 }
#endif

#ifdef WG_TODO

         if (balance.and.(balance_netcdf.ne.0)) then !djm Jan2017
           fhe_32(ix,iy,1,1) = 1.5e0_R8*fne_he(ix,iy,1,1)*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           fhe_52(ix,iy,1,1) = 0.0_R8
           fhe_strange(ix,iy,1,1) = chve(ix,iy,1,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
           if (mdf_fhe.eq.0) then
             fhe_cond(ix,iy,1,1) = -cone0*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1,1) = 2.5e0_R8*fledia(ix,iy,1,1)*
     &         (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))*0.5e0_R8*
     &         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8
          else
             fhe_cond(ix,iy,1,1) = -cone0_mdf*
     &                (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
             fhe_dia(ix,iy,1,1) = 0.0_R8
           endif
           fhe_pschused(ix,iy,1,1) = -fhePSch(ix,iy,1,1)
           fhe_thermj(ix,iy,1,1) = switch%alfTeEh*calf(ix,iy,1,1)*
     &         harm(te(bottomix(ix,iy),bottomiy(ix,iy)),te(ix,iy))*ehy
           fhe_ecrb(ix,iy,1,1) = 0.0_R8
         endif
        endif
#endif


#ifdef WG_TODO

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then        !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhe__qe_32y')                         !srv 16.09.02 {
!       call my_out(70,nx,ny,w2,'b2tfhe__qe_stry')
       call my_out(70,nx,ny,w3,'b2tfhe__qe_ke_gTy')
       call my_out(70,nx,ny,w4,'b2tfhe__qe_alphaTehy')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'b2tfhe__qey5')
       call my_out(70,nx,ny,w5,'b2tfhe__qediay')
       call my_out(70,nx,ny,w6,'b2tfhe__chcey_Pat')
      endif                                                              !srv 16.09.02 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                                .or. iout_b2wdat.eq.4) then        !srv 16.09.02 {
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif                                                              !srv 16.09.02 }
#endif


      if(switch%b2tfhe_iout.ne.0) then                                   !srv 18.06.08 {
       call my_out_us(70,nFc,1,fchstoch,'b2tfhe__fchstochy')             !lk 31.01.08
       call my_out_us(70,nFc,1,dv%fne_he(1,0),'b2tfhe__fnex')
       call my_out_us(70,nFc,1,dv%fne_he(1,1),'b2tfhe__fney')
       call my_out_us(70,nFc,1,dv%floe(1,0),'b2tfhe__floex')
       call my_out_us(70,nFc,1,dv%floe(1,1),'b2tfhe__floey')
       call my_out_us(70,nFc,1,dv%cone(1,0,0),'b2tfhe__conex0')
       call my_out_us(70,nFc,1,dv%cone(1,0,1),'b2tfhe__conex1')
       call my_out_us(70,nFc,1,dv%cone(1,0,2),'b2tfhe__conex2')
       call my_out_us(70,nFc,1,dv%cone(1,1,0),'b2tfhe__coney0')
       call my_out_us(70,nFc,1,dv%cone(1,1,1),'b2tfhe__coney1')
       call my_out_us(70,nFc,1,dv%cone(1,1,2),'b2tfhe__coney2')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfhe__fch_px')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfhe__fch_py')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2tfhe__csig_anx')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2tfhe__csig_any')
       call my_out_us(70,nFc,1,co%chve(1,0),'b2tfhe__chvex')
       call my_out_us(70,nFc,1,co%chve(1,1),'b2tfhe__chvey')
       call my_out_us(70,nFc,1,co%calf(1,0),'b2tfhe__calfx')
       call my_out_us(70,nFc,1,co%calf(1,1),'b2tfhe__calfy')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2tfhe__chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2tfhe__chciy')
       call my_out_us(70,nFc,1,co%chce(1,0),'b2tfhe__chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2tfhe__chcey')
       call my_out_us(70,nCv,0,pl%te,'b2tfhe__te')
       call my_out_us(70,nCv,0,dv%ne,'b2tfhe__ne')
       call my_out_us(70,nCv,0,pl%po,'b2tfhe__po')
       do is = 0, ns-1
        write (charns,'(i3.3)') is
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfhe__na'//charns)
       enddo
       call my_out_us(70,nFc,1,dv%fhe(1,0),'b2tfhe__fhex')
       call my_out_us(70,nFc,1,dv%fhe(1,1),'b2tfhe__fhey')
       call my_out_us(70,nFc,1,dv%fhe_mdf(1,0),'b2tfhe__fhe_mdfx')
       call my_out_us(70,nFc,1,dv%fhe_mdf(1,1),'b2tfhe__fhe_mdfy')
       call my_out_us(70,nFc,1,dv%fhePSch(1,0),'b2tfhe__fhePSchx')
       call my_out_us(70,nFc,1,dv%fhePSch(1,1),'b2tfhe__fhePSchy')
       call my_out_us(70,nFc,1,ehx,'b2tfhe__ehx')                        !srv 13.01.17
       call my_out_us(70,nFc,1,ehy,'b2tfhe__ehy')                        !srv 13.01.17
      endif                                                              !srv 18.06.08 }
#ifdef WG_TODO
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                          .or. iout_b2wdat.eq.4) then              !srv 11.09.09 {
       call my_out(70,nx,ny,fhe(-1,-1,0,0),'b2tfhe__fhe_no_mdfx')
       call my_out(70,nx,ny,fhe(-1,-1,1,1),'b2tfhe__fhe_no_mdfy')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,0,0),'b2tfhe__fhe_mdfx')
!       call my_out(70,nx,ny,fhe_mdf(-1,-1,1,1),'b2tfhe__fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0,0),'b2tfhe__fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1,1),'b2tfhe__fhePSchy')
      endif

#endif
* ..return
      ncall_b2tfhe = ncall_b2tfhe+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfhe_

      end
