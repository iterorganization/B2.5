*-----------------------------------------------------------------------
*.specification

      subroutine b2tfrn (nCv, nFc, nVx, ns, switch, geo, mpg,
     .  pl, dv, co, rt, st_ext)
      use b2mod_types
!      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_transport_fun                                            !srv 13.01.17
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
#ifndef NAGFOR
      use b2mod_b2zhco
     & , only : c_hw_save
#endif
      use b2mod_ad
     & , only : ncall_b2tfrn
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent (in) :: nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl
      type(B2Rates), intent (in) :: rt
      type(B2Coeff), intent (in) :: co
      type(B2StateExt), intent (in) :: st_ext
*   ..input/output arguments
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFRN computes the particle and energy fluxes between
*     neighbouring cells which are to be provided to Eirene.
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle fluxes
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, is, isb, meth
      character chns*3                                                   !srv 16.09.02
      real (kind=R8) ::
     *   flo_eir(nFc,0:1), scur(nFc,0:1), fle_eir(nFc,0:1),              !srv 14.07.10
     &   floe_eir(nFc,0:1), floi_eir(nFc,0:1), floi_vhx(nFc),
     &   fnef(nFc,0:1), fnec(nFc,0:1)
      real (kind=R8) ::
     &   nbf(nFc), pbf(nFc), rzf(nFc), dpb(nFc,0:1), dpbc(nCv),
     &   dnete(nFc,0:1), weight(nFc,2)
      real (kind=R8) ::
     *   wrk0(nFc), wrk1(nFc), wrk2(nFc), wrkVx(nVx),
     *   wrk00(nFc,2), wrk11(nFc,2), wrkCv(nCv),
     *   zeff(nCv), c071(nCv), c071f(nFc)                                !srv 13.01.17
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0             !srv 06.07.06
c      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0             !srv 06.07.06
c      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0             !srv 06.07.06
*   ..procedures
      external xertst, ipgeti, ipgetr
      external b2xvsg, intface
*   ..initialisation
      scur(:,:) = 0.0_R8

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfrn')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..set internal parameters on first call
c      if (ncall_b2tfrn.eq.0) then
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)                 !srv 06.07.06
c        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)                 !srv 06.07.06
c        call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)                 !srv 06.07.06
c        call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)                 !srv 06.07.06
c        call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                 !srv 06.07.06
c        call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                 !srv 06.07.06
c      endif
*   ..extensive tests on first few calls
      if (ncall_b2tfrn.lt.3) then
*    ..test sign of ne
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
      endif

      meth = switch%b2tfnb_discr_meth

      do is = 0, ns-1

*    ..compute partial pressure
       call b2xppb (nCv, rt%rza(:,is), pl%na(:,is), pl%te, pl%ti,
     &   pl%tn, is, dv%pa(:,is))

*   ..interpolate na and pa to cell faces
       weight = 1.0_R8
       call intface(nCv,nFc,mpg%fcCv,weight,pl%na(:,is),nbf)
       call intface(nCv,nFc,mpg%fcCv,weight,dv%pa(:,is),pbf)
       call intface(nCv,nFc,mpg%fcCv,weight,rt%rza(:,is),rzf)

*   ..compute differences of pb
       call diff(nCv, nFc, nVx, 0, geo, mpg, dv%pa(:,is), wrkVx, dpb)
       call gradc_p(nCv, nFc, nVx, 1, geo, mpg, dv%pa(:,is), wrkVx,
     &  dpbc)

*   ..compute wrk0
       call intface(nCv,nFc,mpg%fcCv,geo%fcVol,pl%ua(:,is),wrk0)

*   ..compute convective coefficients (poloidal and radial)
       flo_eir(:,0) = geo%fcPbs*wrk0 + co%cvla(:,0,is) +
     &  switch%fnb_drift_hyb*
     &  (switch%xvecrb*dv%vaecrb(:,0,is)+switch%xwdia*dv%wadia(:,0,is))*
     &   geo%fcS*geo%fcQalf(:,0)
       flo_eir(:,1) = co%cvla(:,1,is) +
     &  switch%fnb_drift_hyb*
     &  (switch%xvecrb*dv%vaecrb(:,1,is)+switch%xwdia*dv%wadia(:,1,is))*
     &   geo%fcS*geo%fcQalf(:,1)

*   ..add contribution from currents
cwdk   still to be checked: vis_per and vis_q also in b2tfnb?
       if (.not.is_neutral(is)) then
        scur(:,0) =
     &   (dv%fchinert_a(:,0,is) + dv%fchvispar_a(:,0,is) +
     &    switch%fnb_vis_per*dv%fchvisper_a(:,0,is) +
     &    dv%fchanml_a(:,0,is) +
     &    switch%fnb_vis_q*dv%fchvisq_a(:,0,is))*switch%b2tfnb_xcur/
     &   (rzf*qe*nbf) +
     &    dv%fchviskt_a(:,0,is)/(rzf*qe*nbf)
        scur(:,1) =
     &   (dv%fchinert_a(:,1,is) + dv%fchvispar_a(:,1,is) +
     &    switch%fnb_vis_per*dv%fchvisper_a(:,1,is) +
     &    dv%fchanml_a(:,1,is) +
     &    switch%fnb_vis_q*dv%fchvisq_a(:,1,is))*switch%b2tfnb_ycur/
     &   (rzf*qe*nbf) +
     &    dv%fchviskt_a(:,1,is)/(rzf*qe*nbf)
        flo_eir   = flo_eir + switch%fnb_drift_hyb*scur
       else
        scur = 0.0_R8 
       end if

*   ..apply discretization scheme for density-part of fna_eir
       call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%na(:,is),flo_eir,
     &  co%cdna(:,:,is),dv%fna_eir(:,:,is),wrk00,wrk11)

*   ..add non-hybrid drift contribution
       dv%fna_eir(:,0,is) = dv%fna_eir(:,0,is) +
     &  (1.0_R8-switch%fnb_drift_hyb)*
     &  ((switch%xvecrb*dv%vaecrb(:,0,is) +
     &   switch%xwdia*dv%wadia(:,0,is))*
     &   geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf
       dv%fna_eir(:,1,is) = dv%fna_eir(:,1,is) +
     &  (1.0_R8-switch%fnb_drift_hyb)*
     &  ((switch%xvecrb*dv%vaecrb(:,1,is) +
     &   switch%xwdia*dv%wadia(:,1,is))*
     &   geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf

*   ..add contribution from pressure diffusion
       dv%fna_eir(:,:,is) = dv%fna_eir(:,:,is) - co%cdpa(:,:,is)*dpb

*   ..add Rhie-Chow term
       do iFc = 1, nFc
        wrk1(iFc) = switch%b2npco_pcm1*0.5_R8*geo%fcPbs(iFc)*
     &   (dv%pcca(mpg%fcCv(iFc,1),1,is)*pl%na(mpg%fcCv(iFc,1),is) +
     &    dv%pcca(mpg%fcCv(iFc,2),1,is)*pl%na(mpg%fcCv(iFc,2),is))/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
        wrk2(iFc) = 0.5_R8*geo%fcQgam(iFc,0)*
     &   (dpbc(mpg%fcCv(iFc,1)) + dpbc(mpg%fcCv(iFc,2)))*
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))
       enddo
       dv%fna_eir(:,0,is) = dv%fna_eir(:,0,is) - wrk1*(dpb(:,0) - wrk2)

      enddo !is

c Drift terms
      if((maxval(dv%facdrift).ne.0.0_R8 .or.
     &    maxval(dv%fac_ExB).ne.0.0_R8)) then

*     ..computation of diamagnetic drift
cwdk    Note: use implementation consistent with vedia (b2tfed) and
cwdk          wadia (b2tfnb)
cwdk          (original code: cell center approach)
        wrkCv = dv%ne*pl%te
        call grad(nCv, nFc, nVx, 0, geo, mpg, wrkCv, wrkVx, dnete)
        wrk0 = geo%fcBb(:,2)/geo%fcBb(:,3)**2
cwdk    to be checked: now we use inverse volume weighting in wkr0
        call intvertex(nCv,nVx,mpg,geo%vxVol,dv%ne,wrkVx)
        call intface(nVx,nFc,mpg%fcVx,weight,wrkVx,wrk1)
        dv%wedia(:,0) =  wrk0*dnete(:,1)/(wrk1*qe)
        dv%wedia(:,1) = -wrk0*dnete(:,0)/(wrk1*qe)
      else
        dv%wedia =0.0_R8
      endif

*   ..compute wrk0
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,dv%ue,wrk0)

      if(switch%b2sigp_style.eq.2) then
        call b2xzef (nCv, ns, rt%rz2, pl%na, dv%ne, zeff, st_ext)
        do iCv = 1, nCv
          c071(iCv) = zeff(iCv)*fce2n(zeff(iCv))
        enddo
      else
        c071 = 0.71_R8                                                   !srv 13.01.17
      endif

cwdk  interpolate to cell faces!
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,c071,c071f)

*   ..velocity-dependent part of the heat flux
      floi_vhx = 0.0_R8                                                 !som 13.07.21
      if (switch%Zhdanov_closure.eq.1 .and.
     &    switch%Zhdanov_test.eq.0 .and.
     &    switch%Zhdanov_vel_heat.eq.1) then                            !som 13.07.21
        call b2tfvh (nCv, nFc, ns, geo, mpg, pl%na, pl%ua, floi_vhx)
      endif

* ..compute fluxes
cwdk  * why the species loop + check for zamax = 1 here?
cwdk    should we not compute fne_eir as sum(fna_eir,is) - fch?
cwdk    (in that case, fch including full diamagnetic current?)
cwdk  * adapted check: added zn(is).eq.1 to avoid double counting
*   ..compute convective coefficients (poloidal and radial)
      fle_eir(:,0) = geo%fcPbs*wrk0 +
     &  (switch%xvecrb*dv%veecrb(:,0) + switch%xwdia*dv%wedia(:,0))*
     &   geo%fcS*geo%fcQalf(:,0)
      fle_eir(:,1) =
     &  (switch%xvecrb*dv%veecrb(:,1) + switch%xwdia*dv%wedia(:,1))*
     &   geo%fcS*geo%fcQalf(:,1)

      dv%fne_eir = 0.0_R8
      do is = 0,ns-1
        if(.not.is_neutral(is)) then
*       ..apply discretization scheme
          !accumulate conductive piece per species(fnec)
          call calcflow(nCv,nFc,nVx,meth,geo,mpg,dv%ne,fle_eir,
     &      co%cdna(:,:,is),wrk00,fnef,fnec)
          dv%fne_eir = dv%fne_eir + fnec
        endif
      enddo
      !add convective piece (once)
      dv%fne_eir = dv%fne_eir + fnef

*   ..compute convective coefficients electron heat flux
      floe_eir(:,0) = 2.5_R8*dv%fne_eir(:,0) -
     &  c071f*dv%fch_p(:,0)/qe*switch%fch_pTe
      floe_eir(:,1) = 2.5_R8*dv%fne_eir(:,1)

*   ..apply discretization scheme
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%te,floe_eir,
     &  co%chce,dv%fhe_eir,wrk00,wrk11)

*   ..compute convective coefficients ion heat flux
      floi_eir = 0.0_R8
      do isb = 0, ns-1
        if(.not.is_neutral(isb)) then
          floi_eir = floi_eir + dv%fna_eir(:,:,isb)
        endif
      enddo
      floi_eir(:,0) = 2.5_R8*floi_eir(:,0) + floi_vhx
      floi_eir(:,1) = 2.5_R8*floi_eir(:,1)

*   ..apply discretization scheme
      call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%ti,floi_eir,
     &  co%chci,dv%fhi_eir,wrk00,wrk11)

      if(switch%b2tfrn_iout.ne.0) then                                   !srv 18.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,dv%fna_eir(1,0,is),
     .   'b2tfrn_fna_eir_th'//chns)
        call my_out_us(70,nFc,1,dv%fna_eir(1,1,is),
     .   'b2tfrn_fna_eir_r'//chns)
        call my_out_us(70,nFc,1,dv%vaecrb(1,0,is),
     &                             'b2tfrn_vecrb_th'//chns)
        call my_out_us(70,nFc,1,dv%vaecrb(1,1,is),
     &                             'b2tfrn_vecrb_r'//chns)
        call my_out_us(70,nFc,1,dv%wadia(1,0,is),
     &                             'b2tfrn_wadia_th'//chns)
        call my_out_us(70,nFc,1,dv%wadia(1,1,is),
     &                             'b2tfrn_wadia_r'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2tfrn_cdna_th'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2tfrn_cdna_r'//chns)
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2tfrn_ua'//chns)
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfrn_na'//chns)
        call my_out_us(70,nFc,1,dv%fchvispar_a(1,0,is),
     &   'b2tfrn_fchvispar_th'//chns)
        call my_out_us(70,nFc,1,dv%fchvispar_a(1,1,is),
     &   'b2tfrn_fchvispar_r'//chns)
        call my_out_us(70,nFc,1,dv%fchinert_a(1,0,is),
     &   'b2tfrn_fchinert_th'//chns)
        call my_out_us(70,nFc,1,dv%fchinert_a(1,1,is),
     &   'b2tfrn_fchinert_r'//chns)
        call my_out_us(70,nFc,1,dv%fchanml_a(1,0,is),
     &   'b2tfrn_fchanml_th'//chns)
        call my_out_us(70,nFc,1,dv%fchanml_a(1,1,is),
     &   'b2tfrn_fchanml_r'//chns)
       enddo
       call my_out_us(70,nFc,1,dv%fne_eir(1,0),'b2tfrn_fne_eir_th')
       call my_out_us(70,nFc,1,dv%fne_eir(1,1),'b2tfrn_fne_eir_r')
       call my_out_us(70,nFc,1,dv%fhe_eir(1,0),'b2tfrn_fhe_eir_th')
       call my_out_us(70,nFc,1,dv%fhe_eir(1,1),'b2tfrn_fhe_eir_r')
       call my_out_us(70,nFc,1,dv%fhi_eir(1,0),'b2tfrn_fhi_eir_th')
       call my_out_us(70,nFc,1,dv%fhi_eir(1,1),'b2tfrn_fhi_eir_r')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfrn_fch_p_th')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfrn_fch_p_r')
       call my_out_us(70,nFc,1,dv%wedia(1,0),'b2tfrn_wedia_th')
       call my_out_us(70,nFc,1,dv%wedia(1,1),'b2tfrn_wedia_r')
       call my_out_us(70,nFc,1,dv%ue,'b2tfrn_ve')
       call my_out_us(70,nCv,0,dv%ne,'b2tfrn_ne')
       call my_out_us(70,nCv,0,pl%te,'b2tfrn_te')
       call my_out_us(70,nCv,0,pl%ti,'b2tfrn_ti')
       call my_out_us(70,nCv,0,pl%po,'b2tfrn_po')
      endif                                                              !srv 18.06.08 }

* ..return
      ncall_b2tfrn = ncall_b2tfrn+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfrn

      end subroutine b2tfrn
