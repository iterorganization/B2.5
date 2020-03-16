*-----------------------------------------------------------------------
*.specification

      subroutine b2tfrn (nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt, 
     .  st_ext)
      use b2mod_types
!      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_transport_fun                                            !srv 13.01.17
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Rates), intent (in) :: rt
      type(B2Coeff), intent (in) :: co
      type(B2StateExt), intent (in) :: st_ext
      logical balance                                                   !djm Jan2017
*   ..input/output arguments
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFRN computes the particle and energy fluxes between
*     neighbouring cells which are to be provided to Eirene
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle fluxes
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, iFc, is, isb, ib
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0              !srv 06.07.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0              !srv 06.07.06
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0              !srv 06.07.06
      integer, save :: no_hybr = 0
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0, iout_b2npmo =0
      integer, save :: style = 2                                         !srv 13.01.17
      integer, save :: meth = 2
      real (kind=R8), save :: xvecrb = 1.0_R8, xwdia = 1.0_R8            !srv 27.01.12
      character chns*3                                                   !srv 16.09.02
      real (kind=R8), save :: drift_hyb = 0.0_R8
      real (kind=R8), save :: fch_pTe = 1.0_R8                           !srv 30.12.09
      real (kind=R8), save :: pcm1 = 1.0_R8                              !srv 05.11.09 30.01.10
      real (kind=R8), save :: xcur = 0.0_R8                              !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                              !srv 07.07.08
      real (kind=R8), save :: vis_per = 0.0_R8, vis_q = 1.0_R8           !srv 02.10.17
      real (kind=R8) ::
     *   flo_eir(nFc,0:1), scur(nFc,0:1), fle_eir(nFc,0:1),   !srv 14.07.10
     &   floe_eir(nFc,0:1), floi_eir(nFc,0:1)
      real (kind=R8) :: 
     &   nef(nFc), nbf(nFc), pbf(nFc), dpb(nFc,0:1), dpbc(nCv),
     &   dnete(nFc,0:1), ve(nFc), weight(nFc,2)
      real (kind=R8) ::
     *   wrk(nCv), wrk0(nFc), wrk1(nFc), wrk2(nFc), wrkVx(nVx),
     *   wrk00(nFc,2), wrk11(nFc,2),
!pb     *   zeff(nCv), c071(nCv), c071f(nCv)                                !srv 13.01.17
     *   zeff(nCv), c071(nCv), c071f(nFc)                                !srv 13.01.17
      real (kind=R8) :: t0,t1,t2
      save ncall
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr
      external b2xvsg, b2xvff, intface
      intrinsic abs, max, min, mod, sqrt
*   ..initialisation
      data ncall/0/ 

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfrn')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tfrn_iout', iout)
        call ipgeti ('b2npc9_iout', iout_b2npc9)
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                          !srv 30.12.09
        call ipgetr ('b2npco_pcm1', pcm1)
        call ipgetr ('b2tfrn_xvecrb', xvecrb)                            !srv 27.01.12
        call ipgetr ('b2tfrn_xwdia', xwdia)                              !srv 27.01.12
        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)                 !srv 06.07.06
        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)                 !srv 06.07.06
        call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)                 !srv 06.07.06
        call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)                 !srv 06.07.06
        call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                 !srv 06.07.06
        call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                 !srv 06.07.06
        call ipgetr ('b2tfnb_xcur', xcur)                                !xpb 27.09.17
        call ipgetr ('b2tfnb_ycur', ycur)                                !srv 17.07.08
        call ipgetr ('b2tfnb_vis_per', vis_per)                          !srv 02.10.17
        call ipgetr ('b2tfnb_vis_q', vis_q)                              !srv 02.10.17
        call ipgeti ('b2tfnb_no_hybr', no_hybr)
          if (no_hybr .eq.1) meth = 0
        call ipgeti ('b2sigp_style', style)                              !srv 13.01.17
        call ipgeti ('b2tfnb_discr_meth', meth)
        call ipgetr ('b2tfnb_drift_hyb',drift_hyb)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
      endif

      do is = 0, ns-1

*    ..compute partial pressure
       call b2xppb (nCv, rt%rza(1,is), pl%na(1,is), pl%te, pl%ti, 
     &   dv%pa(1,is))

*   ..interpolate na and pa to cell faces
       weight = 1.0_R8
       call intface(nCv,nFc,mpg%fcCv,weight,pl%na(:,is),nbf)
       call intface(nCv,nFc,mpg%fcCv,weight,dv%pa(:,is),pbf)

*   ..compute differences of pb
       call diff(nCv, nFc, nVx, 0, geo, mpg, dv%pa(:,is), wrkVx, dpb)
       call gradc_p(nCv, nFc, nVx, 1, geo, mpg, dv%pa(:,is), wrkVx,
     &  dpbc)

*   ..compute wrk0
       call intface(nCv,nFc,mpg%fcCv,geo%fcVol,pl%ua(:,is),wrk0)

*   ..compute convective coefficients (poloidal and radial)
       flo_eir(:,0) = geo%fcPbs*wrk0 + co%cvla(:,0,is) +
     &  drift_hyb*(xvecrb*dv%vaecrb(:,0,is) + xwdia*dv%wadia(:,0,is))*
     &   geo%fcS*geo%fcQalf(:,0)
       flo_eir(:,1) = co%cvla(:,1,is) + 
     &  drift_hyb*(xvecrb*dv%vaecrb(:,1,is) + xwdia*dv%wadia(:,1,is))*
     &   geo%fcS*geo%fcQalf(:,1)

*   ..add contribution from currents
cwdk   condition made consistent with b2tfnb => check zn(is) instead of rza(ix,iy,is)
cwdk   still to be checked: vis_per and vis_q also in b2tfnb?
       if (nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then
        scur(:,0) = 
     &   (dv%fchinert(:,0) + dv%fchvispar(:,0) + 
     &    vis_per*dv%fchvisper(:,0) + dv%fchanml(:,0) + 
     &    vis_q*dv%fchvisq(:,0))*xcur/(qe*nbf)
        scur(:,1) = 
     &   (dv%fchinert(:,1) + dv%fchvispar(:,1) + 
     &    vis_per*dv%fchvisper(:,0) + dv%fchanml(:,1) + 
     &    vis_q*dv%fchvisq(:,1))*ycur/(qe*nbf)
        flo_eir   = flo_eir + drift_hyb*scur
       end if

*   ..apply discretization scheme for density-part of fna_eir
       call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%na(:,is),flo_eir,
     &  co%cdna(:,:,is),dv%fna_eir(:,:,is),wrk00,wrk11)

*   ..add non-hybrid drift contribution
       dv%fna_eir(:,0,is) = dv%fna_eir(:,0,is) + (1.0_R8-drift_hyb)*
     &  ((xvecrb*dv%vaecrb(:,0,is) + xwdia*dv%wadia(:,0,is))*
     &   geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf
       dv%fna_eir(:,1,is) = dv%fna_eir(:,1,is) + (1.0_R8-drift_hyb)*
     &  ((xvecrb*dv%vaecrb(:,1,is) + xwdia*dv%wadia(:,1,is))*
     &   geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf

*   ..add contribution from pressure diffusion
       dv%fna_eir(:,:,is) = dv%fna_eir(:,:,is) - co%cdpa(:,:,is)*dpb

*   ..add Rhie-Chow term
       do iFc = 1, nFc
        wrk1(iFc) = pcm1*0.5_R8*geo%fcPbs(iFc)*
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
        call grad(nCv, nFc, nVx, 0, geo, mpg, dv%ne*pl%te, wrkVx, dnete)
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

      if(style.eq.2) then
        call b2xzef (nCv, ns, rt%rz2, pl%na, dv%ne, zeff, st_ext) 
        do iCv = 1, nCv
          c071(iCv) = fce2n(zeff(iCv))
        enddo
      else
        c071 = 0.71e0_R8                                                 !srv 13.01.17
      endif

cwdk  interpolate to cell faces!
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,c071,c071f)

* ..compute fluxes
cwdk  * why the species loop + check for zamax = 1 here?
cwdk    should we not compute fne_eir as sum(fna_eir,is) - fch?
cwdk    (in that case, fch including full diamagnetic current?)
cwdk  * adapted check: added zn(is).eq.1 to avoid double counting
      do is = 0,ns-1
       if(nint(zn(is)).eq.1.and.nint(zamax(is)).eq.1) then

*   ..compute convective coefficients (poloidal and radial)
        fle_eir(:,0) = geo%fcPbs*wrk0 + 
     &    (xvecrb*dv%veecrb(:,0) + xwdia*dv%wedia(:,0))*
     &     geo%fcS*geo%fcQalf(:,0)
        fle_eir(:,1) = 
     &    (xvecrb*dv%veecrb(:,1) + xwdia*dv%wedia(:,1))*
     &     geo%fcS*geo%fcQalf(:,1)

*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,dv%ne,fle_eir,
     &    co%cdna(:,:,is),dv%fne_eir,wrk00,wrk11)

*   ..compute convective coefficients electron heat flux
        floe_eir(:,0) = 2.5_R8*dv%fne_eir(:,0) - 
     &    c071f*dv%fch_p(:,0)/qe*fch_pTe
        floe_eir(:,1) = 2.5_R8*dv%fne_eir(:,1)

*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%te,floe_eir,
     &    co%chce,dv%fhe_eir,wrk00,wrk11)

*   ..compute convective coefficients ion heat flux
        floi_eir = 0.0_R8
        do isb = 0, ns-1
         if(.not.is_neutral(isb)) then
          floi_eir = floi_eir + dv%fna_eir(:,:,isb)
         endif
        enddo
        floi_eir(:,0) = 2.5_R8*floi_eir(:,0)
        floi_eir(:,1) = 2.5_R8*floi_eir(:,1)

*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%ti,floi_eir,
     &    co%chci,dv%fhi_eir,wrk00,wrk11)

       endif
      enddo

      if(iout.ne.0) then                                                 !srv 18.06.08 {
       do is = 0,ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,dv%fna_eir(1,0,is),
     .   'b2tfrn_fna_eirx'//chns)
        call my_out_us(70,nFc,1,dv%fna_eir(1,1,is),
     .   'b2tfrn_fna_eiry'//chns)
        call my_out_us(70,nFc,1,dv%vaecrb(1,0,is),'b2tfrn_vecrbx'//chns)
        call my_out_us(70,nFc,1,dv%vaecrb(1,1,is),'b2tfrn_vecrby'//chns)
        call my_out_us(70,nFc,1,dv%wadia(1,0,is),'b2tfrn_wadiax'//chns)
        call my_out_us(70,nFc,1,dv%wadia(1,1,is),'b2tfrn_wadiay'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2tfrn_cdnax'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2tfrn_cdnay'//chns)
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2tfrn_ua'//chns)
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfrn_na'//chns)
       enddo
       call my_out_us(70,nFc,1,dv%fne_eir(1,0),'b2tfrn_fne_eirx')
       call my_out_us(70,nFc,1,dv%fne_eir(1,1),'b2tfrn_fne_eiry')
       call my_out_us(70,nFc,1,dv%fhe_eir(1,0),'b2tfrn_fhe_eirx')
       call my_out_us(70,nFc,1,dv%fhe_eir(1,1),'b2tfrn_fhe_eiry')
       call my_out_us(70,nFc,1,dv%fhi_eir(1,0),'b2tfrn_fhi_eirx')
       call my_out_us(70,nFc,1,dv%fhi_eir(1,1),'b2tfrn_fhi_eiry')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfrn_fch_px')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfrn_fch_py')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfrn_fchvisparx')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfrn_fchvispary')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfrn_fchinertx')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfrn_fchinerty')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfrn_fchanmlx')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfrn_fchanmly')
       call my_out_us(70,nFc,1,dv%wedia(1,0),'b2tfrn_wediax')
       call my_out_us(70,nFc,1,dv%wedia(1,1),'b2tfrn_wediay')
       call my_out_us(70,nFc,1,dv%ue,'b2tfrn_ve')
       call my_out_us(70,nCv,0,dv%ne,'b2tfrn_ne')
       call my_out_us(70,nCv,0,pl%te,'b2tfrn_te')
       call my_out_us(70,nCv,0,pl%ti,'b2tfrn_ti')
       call my_out_us(70,nCv,0,pl%po,'b2tfrn_po')
      endif                                                              !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfrn

      end subroutine b2tfrn
