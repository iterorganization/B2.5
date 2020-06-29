*-----------------------------------------------------------------------
*.specification
      subroutine b2tfhi_ (nCv, nFc, nVx, ns, geo, mpg, pl, dv, co, rt,  !srv 13.10.06
     &  balance)                                                        !djm Jan2017
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance                                                 !djm Jan2017
!WG_TODO     & , only : fhi_32, fhi_52, fhi_cond, fhi_dia, fhi_ecrb,
!WG_TODO     &          fhi_strange, fhi_pschused, fhi_inert, fhi_vispar,
!WG_TODO     &          fhi_visper, fhi_visq, fhi_anml, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Rates), intent (in) :: rt
      logical balance                                                   !djm Jan2017
*   ..input/output arguments (unspecified on entry)
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
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
      integer ncall, iCv, iFc, is, ib, ixl, iyl, ixb, iyb
      integer, save :: iout = 0, iout_b2wdat = 0,                       !srv 16.09.02
     &                 iout_b2nph9 = 0, iout_b2npht = 0
      integer, save :: lim_flux = 0, mdf_fhi = 0                        !srv 04.02.03 13.10.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0             !srv 06.07.06
      integer, save :: conductive_limit = 1                             !srv 16.06.08 17.01.12
      integer, save :: b2tfhi_hybr2=0, b2_upwind=0, flo53=0
      integer, save :: b2kin=0, b2vis=0
      integer, save :: meth = 2, tn_style = 0
      real (kind=R8), save :: drift_hyb = 1.0_R8
      real (kind=R8), save :: PSch=1.0e0_R8                             !srv 17.01.07
      real (kind=R8), save :: xcur = 0.0_R8                             !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                             !srv 07.07.08
      real (kind=R8), save :: vis_per = 0.0_R8, vis_q = 1.0_R8          !srv 02.10.17
      real (kind=R8) ::
     *   floi0(nFc,0:1), flomx, coni0(nFc,0:1), vdia_par,     !srv 17.05.99
     &   floi0_mdf(nFc,0:1), wght(nFc,2), dumm0(nFc,0:1),
     &   dumm1(nFc,0:1), fhi0_mdf(nFc,0:1),
     *   flon0(nFc,0:1), conn0(nFc,0:1), fhn0(nFc,0:1)
cpax      
      real (kind=R8) ::
     *   flidia(nFc,0:1), facdriftm, nif(nFc), nnf(nFc),               !srv 25.01.02
     *   tif(nFc), tnf(nFc),
     *   nati(nCv), nati2(nCv), nativ(nVx), nati2v(nVx),
     *   dnati2(nFc,0:1), cddi(nFc,0:1)
#ifdef WG_TODO
      real (kind=R8), allocatable ::
     *    w1(:,:), w2(:,:), w3(:,:), w4(:,:), w5(:,:), w6(:,:)          !srv 16.09.02
#endif
      real (kind=R8) :: upfac, hybfac, scur(nFc,0:1) !djm Jan2017
      real (kind=R8) :: cutlo
      parameter (cutlo = 1.0e-60_R8)
      save ncall, vdia_par                                              !srv 17.05.99
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff
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
       call ipgetr ('b2tfhi_drift_hyb',drift_hyb)
       call ipgeti ('b2mn_tn_style', tn_style)
      endif                                                             !srv
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ni, ti, chci
       call b2xvsg (nCv, dv%ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nFc, co%chci(:,0)*geo%fcQalf(:,0),1,'chci0','.ge.')
       call b2xvsg (nFc, co%chci(:,1)*geo%fcQalf(:,1),1,'chci1','.ge.')

       if (tn_style.ne.0.and.mdf_fhi.eq.1) then
         call xerrab('mdf_fhi not yet available for',
     &               'non-default energy equation')
       endif
      endif

*   ..set dummy arrays
      dumm0 = 0.0_R8
      dumm1 = 0.0_R8


      facdriftm  = maxval(dv%facdrift)
      dv%fhiPSch = 0.0e0_R8                                                !srv 13.10.06
      dv%coni    = 0.0e0_R8                                                !srv 02.01.07
      dv%floi    = 0.0e0_R8
      if(facdriftm.ne.0.0e0_R8) then
       do is = 0, ns-1
        if(.not.is_neutral(is)) then

*   ..precompute face and vertex values of nati and nati2
         nati  = pl%na(:,is)*pl%ti
         nati2 = pl%na(:,is)*pl%ti*pl%ti
         call intvertex(nCv,nVx,mpg,geo%vxVol,nati,nativ)
         call intvertex(nCv,nVx,mpg,geo%vxVol,nati2,nati2v)
         call diff(nCv, nFc, nVx, 1, geo, mpg, nati2, nati2v, dnati2)

*   ..compute P.Sch heat flux
         dv%fhiPSch(:,0) =  2.5_R8*PSch*dv%facdrift*co%cddi(:,0,is)*
     &    dnati2(:,1)
         dv%fhiPSch(:,1) = -2.5_R8*PSch*dv%facdrift*co%cddi(:,1,is)*
     &    dnati2(:,0)
*     ..linearization coefficients for P.Sch flows
         cddi(:,0) = -5.0_R8*PSch*dv%facdrift*co%cddi(:,1,is)
         cddi(:,1) =  5.0_R8*PSch*dv%facdrift*co%cddi(:,0,is)
         call calccoef(nCv,nFc,nVx,0,geo,mpg,dumm0,cddi,floi0,coni0)
         do iFc = 1, nFc
          dv%floi(iFc,0)   = dv%floi(iFc,0) + 
     &       0.5_R8*floi0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,1)) + nati(mpg%fcCv(iFc,2))) -
     &       coni0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,2)) - nati(mpg%fcCv(iFc,1)))
          dv%floi(iFc,1)   = dv%floi(iFc,1) + 
     &       0.5_R8*floi0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,1))+nativ(mpg%fcVx(iFc,2))) -
     &       coni0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,2))-nativ(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into coni(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to dv%coni(,,1) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          dv%coni(iFc,0,0) = dv%coni(iFc,0,0) +
     &       0.5_R8*coni0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,1))+nati(mpg%fcCv(iFc,2))) -
     &       0.25_R8*floi0(iFc,0)*
     &      (nati(mpg%fcCv(iFc,2))-nati(mpg%fcCv(iFc,1)))
          dv%coni(iFc,1,0) = dv%coni(iFc,1,0) +
     &       0.5_R8*coni0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,1))+nativ(mpg%fcVx(iFc,2))) -
     &       0.25_R8*floi0(iFc,1)*
     &      (nativ(mpg%fcVx(iFc,2))-nativ(mpg%fcVx(iFc,1)))
         enddo
        endif
       enddo
      endif

      wght   = 1.0_R8
      flidia = 0.0_R8                                                   !srv 25.01.02 {
      if(facdriftm.ne.0.0e0_R8) then                                    !srv 18.01.07
       do is = 0, ns-1
        if ( .not.is_neutral(is) ) then
*      ..interpolate na(:,is) to cell faces
         call intface(nCv,nFc,mpg%fcCv,wght,pl%na(:,is),nif)
*      ..compute flidia
         flidia(:,0) = flidia(:,0) + 2.5e0_R8*vdia_par*
     &     nif*dv%vadia(:,0,is)*geo%fcS*geo%fcQalf(:,0)
         flidia(:,1) = flidia(:,1) + 2.5e0_R8*vdia_par*
     &     nif*dv%vadia(:,1,is)*geo%fcS*geo%fcQalf(:,1)
        endif
       enddo                                                            !srv 25.01.02 }
      endif                                                             !srv 18.01.07

#ifdef WG_TODO
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &   .or. iout_b2wdat.eq.4) then                                    !srv 16.09.02 {
       allocate(
     *    w1(-1:nx,-1:ny), w2(-1:nx,-1:ny), w3(-1:nx,-1:ny)
     *   ,w4(-1:nx,-1:ny), w5(-1:nx,-1:ny), w6(-1:nx,-1:ny)
     *         )
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                             !srv 16.09.02 }
#endif

*   ..interpolate ni and ti to cell faces

      call intface(nCv,nFc,mpg%fcCv,wght,pl%ti,tif)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%tn,tnf)        !! in case of default energy equation, tn will be equal to ti

      if (tn_style.eq.0) then                              !! default combined ion-neutral energy equation
        call intface(nCv,nFc,mpg%fcCv,wght,dv%ni(:,0),nif)
      else                                                 !! pure ion and separate ion-neutral energy equation
        call intface(nCv,nFc,mpg%fcCv,wght,dv%ni(:,1),nif)
        call intface(nCv,nFc,mpg%fcCv,wght,dv%ni(:,0)-dv%ni(:,1),nnf)
      endif


*   ..contributions from currents
      scur(:,0) = 1.5e0_R8/qe*(dv%fchinert(:,0) + dv%fchvispar(:,0)+
     &  dv%fchanml(:,0) + dv%fchvisper(:,0) + dv%fchvisq(:,0))*xcur
      scur(:,1) = 1.5e0_R8/qe*(dv%fchinert(:,1) + dv%fchvispar(:,1)+
     &  dv%fchanml(:,1) + vis_per*dv%fchvisper(:,1) + 
     &  vis_q*dv%fchvisq(:,1))*ycur

      if (iout.eq.2) then
        call my_out_us(70,nFc,1,co%chcn(:,0),'chcn0')
        call my_out_us(70,nFc,1,co%chci(:,0),'chci0')
        call my_out_us(70,nFc,1,co%chcn(:,1),'chcn1')
        call my_out_us(70,nFc,1,co%chci(:,1),'chci1')
      endif


*   ..convective coefficients
      if (tn_style.eq.0) then
        floi0(:,0) = 1.5e0_R8*dv%fni_he(:,0) + co%chvi(:,0)*nif + 
     &    flidia(:,0) + scur(:,0)
        floi0(:,1) = 1.5e0_R8*dv%fni_he(:,1) + co%chvi(:,1)*nif + 
     &    flidia(:,1) + scur(:,1)

        flon0(:,0) = 0.0_R8
        flon0(:,1) = 0.0_R8
      elseif (tn_style.eq.1) then
!! AFN_TODO when more contributions are considered (e.g. drifts), or multi-species cases, 
!! additions (floi0 = floi0 + ...) must be considered instead of overwrites (floi0 = ...)
        do is= 0,ns-1
          if (.not.is_neutral(is)) then
            floi0(:,0) = 1.5e0_R8*dv%fna_he(:,0,is)
            floi0(:,1) = 1.5e0_R8*dv%fna_he(:,1,is)
          endif 
          flon0(:,0) = 0.0_R8
          flon0(:,1) = 0.0_R8
        enddo
      else
        do is= 0,ns-1
          if (.not.is_neutral(is)) then
            floi0(:,0) = 1.5e0_R8*dv%fna_he(:,0,is)
            floi0(:,1) = 1.5e0_R8*dv%fna_he(:,1,is)
          else
            flon0(:,0) = 1.5e0_R8*dv%fna_he(:,0,is)
            flon0(:,1) = 1.5e0_R8*dv%fna_he(:,1,is)
          endif
        enddo
      endif


      if (tn_style.eq.0) then !! default combined ion-neutral energy equation
      
        floi0_mdf(:,0) = 1.5e0_R8*dv%fni_he(:,0) + co%chvi(:,0)*nif + 
     &    scur(:,0)
        floi0_mdf(:,1) = 1.5e0_R8*dv%fni_he(:,1) + co%chvi(:,1)*nif + 
     &    scur(:,1)

*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%ti,floi0_mdf,
     &    co%chci,fhi0_mdf,dumm0,dumm1)

*   ..compute fluxes
        dv%fhi(:,0) = fhi0_mdf(:,0) +  (floi0(:,0) - floi0_mdf(:,0))*tif
        dv%fhi(:,1) = fhi0_mdf(:,1) +  (floi0(:,1) - floi0_mdf(:,1))*tif

        dv%fhi_mdf(:,0)  = fhi0_mdf(:,0) - dv%fhiPSch(:,0)
        dv%fhi_mdf(:,1)  = fhi0_mdf(:,1) - dv%fhiPSch(:,1)

*   ..compute floi, coni
        floi0 = dv%floi
        coni0 = dv%coni(:,:,0)
        call calccoef(nCv,nFc,nVx,meth,geo,mpg,floi0_mdf,co%chci,
     &    dv%floi,dv%coni)
      !wdk now we always add PSch linearization; do we need to condition
      !    this on mdf /= 0?
        dv%floi = dv%floi + floi0
        dv%coni(:,:,0) = dv%coni(:,:,0) + coni0

      elseif (tn_style.eq.1) then
        if(mdf_fhi.eq.1) then
          call xerrab('non-default energy eqs not yet for mdf')
        endif
*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%ti,floi0,
     &    co%chci,dv%fhi,dumm0,dumm1)
*   ..compute floi, coni
        call calccoef(nCv,nFc,nVx,meth,geo,mpg,floi0,co%chci,
     &    dv%floi,dv%coni)
      else
        if(mdf_fhi.eq.1) then
          call xerrab('non-default energy eqs not yet for mdf')
        endif

*   ..apply discretization scheme
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%ti,floi0,
     &    co%chci,dv%fhi,dumm0,dumm1)
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%tn,flon0,
     &    co%chcn,dv%fhn,dumm0,dumm1)

*   ..compute floi, coni
        call calccoef(nCv,nFc,nVx,meth,geo,mpg,floi0,co%chci,
     &    dv%floi,dv%coni)
*   ..compute flon, conn
        call calccoef(nCv,nFc,nVx,meth,geo,mpg,flon0,co%chcn,
     &    dv%flon,dv%conn)

      endif



#ifdef WG_TODO
         if(lim_flux.eq.1 .and. conductive_limit.eq.0) then             !srv 17.01.12
          call xerrab ('b2tfhi_: option not yet available for WG')



csw 17feb2012
          if(no_Ptncr_x_hi.eq.0) then
            coni0 = hybr(floi0,chci(ix,iy,0,0))-abs(floi0/2.e0)
          else
            coni0 = chci(ix,iy,0,0)-abs(floi0/2.e0)
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

           coni0=dim(chci(ix,iy,0,0)/(1.d0+t0), abs(floi0)*
     .      max(0.0_R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .      max(cutlo,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
          endif


csw 27nov2012  added check upwind (SOLPS4 style)
          if(b2_upwind .ne. 0) then
           fhi(ix,iy,0,0) = upwind(floi0,ti(leftix(ix,iy),leftiy(ix,iy))
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
             fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_52(ix,iy,0,0) = 0.0_R8
             fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &        (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &        (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &        (1.0_R8-upfac)*ti(ix,iy))
             fhi_dia(ix,iy,0,0) = flidia(ix,iy,0,0)*
     &                          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                           (1.0_R8-upfac)*ti(ix,iy))
             fhi_cond(ix,iy,0,0) = -coni0*dti
             fhi_ecrb(ix,iy,0,0) = 0.0_R8
             fhi_pschused(ix,iy,0,0) = 0.0_R8
             fhi_inert(ix,iy,0,0) = 1.5e0_R8/qe*fchinert(ix,iy,0,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_vispar(ix,iy,0,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visper(ix,iy,0,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0,0)*
     &                    xcur*(upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_visq(ix,iy,0,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
             fhi_anml(ix,iy,0,0) = 1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &                         (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                          (1.0_R8-upfac)*ti(ix,iy))
           endif
          else
           fhi(ix,iy,0,0) = floi0*(ti(leftix(ix,iy),leftiy(ix,iy))
     &                                  +ti(ix,iy))/2.0_R8
     &                    -coni0*dti0
           if (balance.and.
     &         (balance_netcdf.ne.0).and.
     &         (mdf_fhi.eq.0)) then !djm Jan2017
             fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_52(ix,iy,0,0) = 0.0_R8
             fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &        (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_dia(ix,iy,0,0) = flidia(ix,iy,0,0)*
     &        (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_cond(ix,iy,0,0) = -coni0*dti
             fhi_ecrb(ix,iy,0,0) = 0.0_R8
             fhi_pschused(ix,iy,0,0) = 0.0_R8
             fhi_inert(ix,iy,0,0) = 1.5e0_R8/qe*fchinert(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_vispar(ix,iy,0,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visper(ix,iy,0,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_visq(ix,iy,0,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             fhi_anml(ix,iy,0,0) = 1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
           endif
          endif


csw 2014 add kinetic term (ONLY FOR TESTING w/ SOLPS4 discretization)
          if(b2kin .ne. 0) then
           do is=0,ns-1
             fhi(ix,iy,0,0) = fhi(ix,iy,0,0) + fna_he(ix,iy,0,0,is) *
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
             fhi(ix,iy,0,0) = fhi(ix,iy,0,0) - 
     .          0.5_R8*cvsa(ix,iy,0,0,is)*
     .         (ua(ix,iy,is) + ua(leftix(ix,iy),leftiy(ix,iy),is))*
     .         (ua(ix,iy,is) - ua(leftix(ix,iy),leftiy(ix,iy),is))
           enddo
          endif



          if(mdf_fhi.eq.0) then

csw 30apr2014 for SOLPS4 style
           if(b2tfhi_hybr2 .ne.0) then
            floi(ix,iy,0,0) = floi0
            coni(ix,iy,0,0,0) = coni0
           else
            floi(ix,iy,0,0) = floi0-
     &      sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &      min(0.0e0_R8,t1)**2
            coni(ix,iy,0,0,0) = abs(floi0/2.e0)+coni0/(1+t0)              !srv 02.01.07
           endif

          else

csw 17feb2012
           if(no_Ptncr_x_hi.eq.0) then
            coni0_mdf =
     .        hybr(floi0_mdf,chci(ix,iy,0,0))-abs(floi0_mdf/2.0)          !srv 17.06.09 {
           else
            coni0_mdf = chci(ix,iy,0,0)-abs(floi0_mdf/2.0)
           endif

csw 13mar2013 SOLPS4 style (overwrite coni0_mdf)
           if(b2tfhi_hybr2 .ne. 0) then
            coni0_mdf=dim(chci(ix,iy,0,0), abs(floi0_mdf)*
     .       max(0._R8,min(ti(ix,iy),ti(leftix(ix,iy), leftiy(ix,iy))))/
     .       max(1.e-30_R8,ti(ix,iy)+ti(leftix(ix,iy), leftiy(ix,iy))))
           endif

           t0_mdf =coni0_mdf*abs(dti)/(flomx*
     &                  max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
           dti0_mdf = dti/(1.0_R8+t0_mdf)
           t1_mdf = t0_mdf/(1.0_R8+t0_mdf)                              !srv 17.06.09 {

csw 27nov2012 added check upwind (SOLPS4 style)
           if(b2_upwind .ne. 0) then
             fhi_mdf(ix,iy,0,0) = upwind(floi0_mdf,
     &           ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)) 
     &       -coni0_mdf*dti0 - fhiPSch(ix,iy,0,0)                       !srv 17.06.09 }
             floi(ix,iy,0,0) = floi0_mdf-
     &        sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &        min(0.0e0_R8,t1_mdf)**2
             coni(ix,iy,0,0,0) = coni0_mdf/(1.0_R8+t0_mdf)              !srv 02.01.07
             if (balance.and.
     &           (balance_netcdf.ne.0).and.
     &           (mdf_fhi.ne.0)) then !djm Jan2017
               if(floi0_mdf.ge.0.0_R8) then
                 upfac=1.0_R8
               else
                 upfac=0.0_R8
               endif
               fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &                           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &                            (1.0_R8-upfac)*ti(ix,iy))
               fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &           (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &           (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &           (1.0_R8-upfac)*ti(ix,iy))
               fhi_52(ix,iy,0,0) = 0.0_R8
               fhi_cond(ix,iy,0,0) = -coni0_mdf*dti0
               fhi_dia(ix,iy,0,0) = 0.0_R8
               fhi_ecrb(ix,iy,0,0) = 0.0_R8
               fhi_pschused(ix,iy,0,0) = -fhiPSch(ix,iy,0,0)
               fhi_inert(ix,iy,0,0) =
     &           1.5e0_R8/qe*fchinert(ix,iy,0,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_vispar(ix,iy,0,0) =
     &           1.5e0_R8/qe*fchvispar(ix,iy,0,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_visper(ix,iy,0,0) =
     &           1.5e0_R8/qe*fchvisper(ix,iy,0,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_visq(ix,iy,0,0) =
     &           1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
               fhi_anml(ix,iy,0,0) =
     &           1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &          (upfac*ti(leftix(ix,iy),leftiy(ix,iy))+
     &          (1.0_R8-upfac)*ti(ix,iy))
             endif
           else
             fhi_mdf(ix,iy,0,0) = floi0_mdf*
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0 
     &       -abs(floi0_mdf/2.0)*dti-coni0_mdf*dti0 - fhiPSch(ix,iy,0,0)!srv 17.06.09 }
             floi(ix,iy,0,0) = floi0_mdf-
     &        sign(flomx,ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     &        min(0.0e0_R8,t1_mdf)**2
             coni(ix,iy,0,0,0)=
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
               fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &           hybfac*1.5e0_R8*fni_he(ix,iy,0,0)/2.0_R8*dti
               fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &          (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &          (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8-
     &           hybfac*chvi(ix,iy,0,0)*
     &          (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8/
     &           2.0_R8*dti
               fhi_52(ix,iy,0,0) = 0.0_R8
               fhi_cond(ix,iy,0,0) = -coni0_mdf*dti0
               fhi_dia(ix,iy,0,0) = 0.0_R8
               fhi_ecrb(ix,iy,0,0) = 0.0_R8
               fhi_pschused(ix,iy,0,0) = -fhiPSch(ix,iy,0,0)
               fhi_inert(ix,iy,0,0) =
     &              1.5e0_R8/qe*fchinert(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_vispar(ix,iy,0,0) =
     &              1.5e0_R8/qe*fchvispar(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_visper(ix,iy,0,0) =
     &              1.5e0_R8/qe*fchvisper(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_visq(ix,iy,0,0) =
     &              1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
               fhi_anml(ix,iy,0,0) =
     &              1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
             endif
           endif
          endif
#endif



#ifdef WG_TODO
                       

            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.eq.0)) then !djm Jan2017
              fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_52(ix,iy,0,0) = 0.0_R8
              fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &               (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_dia(ix,iy,0,0) = flidia(ix,iy,0,0)*
     &               (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0
              fhi_cond(ix,iy,0,0) = -coni0*dti0
              fhi_ecrb(ix,iy,0,0) = 0.0_R8
              fhi_pschused(ix,iy,0,0) = 0.0_R8
              fhi_inert(ix,iy,0,0) = 1.5e0_R8/qe*fchinert(ix,iy,0,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_vispar(ix,iy,0,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visper(ix,iy,0,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0,0)*
     &           xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visq(ix,iy,0,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_anml(ix,iy,0,0) = 1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
            endif
            if (balance.and.
     &          (balance_netcdf.ne.0).and.
     &          (mdf_fhi.ne.0)) then !djm Jan2017
              fhi_32(ix,iy,0,0) = 1.5e0_R8*fni_he(ix,iy,0,0)*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_52(ix,iy,0,0) = 0.0_R8
              fhi_strange(ix,iy,0,0) = chvi(ix,iy,0,0)*
     &         (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_dia(ix,iy,0,0) = 0.0_R8
              fhi_cond(ix,iy,0,0) = -coni0_mdf*dti0
              fhi_ecrb(ix,iy,0,0) = 0.0_R8
              fhi_pschused(ix,iy,0,0) = -fhiPSch(ix,iy,0,0)
              fhi_inert(ix,iy,0,0) = 1.5e0_R8/qe*fchinert(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_vispar(ix,iy,0,0) = 1.5e0_R8/qe*fchvispar(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visper(ix,iy,0,0) = 1.5e0_R8/qe*fchvisper(ix,iy,0,0)*
     &         xcur*(ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_visq(ix,iy,0,0) = 1.5e0_R8/qe*fchvisq(ix,iy,0,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
              fhi_anml(ix,iy,0,0) = 1.5e0_R8/qe*fchanml(ix,iy,0,0)*xcur*
     &         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
            endif
          endif

         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.eq.0)) then !djm Jan2017
           fhi_32(ix,iy,1,1) = 1.5e0_R8*fni_he(ix,iy,1,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_strange(ix,iy,1,1) = chvi(ix,iy,1,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_dia(ix,iy,1,1) = flidia(ix,iy,1,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_inert(ix,iy,1,1) = 1.5e0_R8/qe*fchinert(ix,iy,1,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_vispar(ix,iy,1,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1,1)*
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visper(ix,iy,1,1) =
     &      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1,1)*ycur*               !srv 02.10.17
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visq(ix,iy,1,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1,1)*  !srv 02.10.17
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_anml(ix,iy,1,1) = 1.5e0_R8/qe*fchanml(ix,iy,1,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_52(ix,iy,1,1) = 0.0_R8
           fhi_cond(ix,iy,1,1) = -coni0*
     &      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1,1) = 0.0_R8
           fhi_pschused(ix,iy,1,1) = 0.0_R8
         endif
         if (balance.and.
     &       (balance_netcdf.ne.0).and.
     &       (mdf_fhi.ne.0)) then !djm Jan2017
           fhi_32(ix,iy,1,1) = 1.5e0_R8*fni_he(ix,iy,1,1)*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_strange(ix,iy,1,1) = chvi(ix,iy,1,1)*
     &      (ni(bottomix(ix,iy),bottomiy(ix,iy),0)+ni(ix,iy,0))/2.0_R8*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_dia(ix,iy,1,1) = 0.0_R8
           fhi_inert(ix,iy,1,1) = 1.5e0_R8/qe*fchinert(ix,iy,1,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_vispar(ix,iy,1,1) = 1.5e0_R8/qe*fchvispar(ix,iy,1,1)*
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visper(ix,iy,1,1) =
     &      1.5e0_R8/qe*vis_per*fchvisper(ix,iy,1,1)*ycur*                  !srv 02.10.17
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_visq(ix,iy,1,1) = 1.5e0_R8/qe*vis_q*fchvisq(ix,iy,1,1)*   !srv 02.10.17
     &      ycur*(ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_anml(ix,iy,1,1) = 1.5e0_R8/qe*fchanml(ix,iy,1,1)*ycur*
     &      (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8
           fhi_52(ix,iy,1,1) = 0.0_R8
           fhi_cond(ix,iy,1,1) = -coni0_mdf*
     &      (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
           fhi_ecrb(ix,iy,1,1) = 0.0_R8
           fhi_pschused(ix,iy,1,1) = -fhiPSch(ix,iy,1,1)
         endif         
        endif
#endif


#ifdef WG_TODO


          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                    .or. iout_b2wdat.eq.4) then                   !srv 16.09.02 {
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
           w5(ix,iy) = flidia(ix,iy,0,0)*
     &              (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8
          endif   

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &   .or. iout_b2wdat.eq.4) then                   !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhi__qi_32x')
       call my_out(70,nx,ny,w2,'b2tfhi__qi_ki_gTx')
!       w3=w1+w2
!       call my_out(70,nx,ny,w3,'qix3')
       call my_out(70,nx,ny,w5,'b2tfhi__qidiax')
       call my_out(70,nx,ny,w6,'b2tfhi__chcix_Pat')
       call my_out(70,nx,ny,fhiPSch(-1,-1,0,0),'b2tfhi__fhiPSchx')
       w1=0.0e0_R8; w2=0.0e0_R8; w3=0.0e0_R8; w4=0.0e0_R8; w5=0.0e0_R8;
       w6=0.0e0_R8;
      endif                                                             !srv 16.09.02 }


        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &          .or. iout_b2wdat.eq.4) then                 !srv 16.09.02 {
         if(mdf_fhi.eq.0) then                                          !srv 13.10.06 {
          w1(ix,iy)=1.5e0_R8*fni_he(ix,iy,1,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w2(ix,iy)=flidia(ix,iy,1,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w3(ix,iy)=-coni0*
     1    (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=0.0e0_R8
          w6(ix,iy)=  coni0
         else
          w1(ix,iy)=1.5e0_R8*fni_he(ix,iy,1,1)*
     1   (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2
          w2(ix,iy)=0.0e0_R8
          w3(ix,iy)=-coni0_mdf*
     1     (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
          w4(ix,iy)=0.0e0_R8
          w6(ix,iy)=  coni0_mdf
         endif
         w5(ix,iy) = flidia(ix,iy,1,1)*
     &                 (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2

        endif                                                           !srv 16.09.02 }

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &           .or. iout_b2wdat.eq.4) then                   !srv 16.09.02 {
       call my_out(70,nx,ny,w1,'b2tfhi__qi_32y')
!       call my_out(70,nx,ny,w2,'qiy2')
       call my_out(70,nx,ny,w3,'b2tfhi__qi_ki_gTy')
!       call my_out(70,nx,ny,w4,'qiy4')
!       w5=w1+w2+w3+w4
!       call my_out(70,nx,ny,w5,'qiy5')
       call my_out(70,nx,ny,w5,'b2tfhi__qidiay')
       call my_out(70,nx,ny,w6,'b2tfhi__chciy_Pat')                     !srv 16.09.02 }
       call my_out(70,nx,ny,fhiPSch(-1,-1,1,1),'b2tfhi__fhiPSchy')
       call my_out(70,nx,ny,fhi(-1,-1,0,0),'b2tfhi__fhix')
       call my_out(70,nx,ny,fhi(-1,-1,1,1),'b2tfhi__fhiy')
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &           .or. iout_b2wdat.eq.4) then
       deallocate( w1, w2, w3, w4, w5, w6 )
      endif
#endif

      if(iout.ne.0) then                                                !srv 18.06.08 {
       call my_out_us(70,nFc,1,dv%floi(1,0),'b2tfhi__floix')
       call my_out_us(70,nFc,1,dv%floi(1,1),'b2tfhi__floiy')
       call my_out_us(70,nFc,1,dv%coni(1,0,0),'b2tfhi__conix0')
       call my_out_us(70,nFc,1,dv%coni(1,0,1),'b2tfhi__conix1')
       call my_out_us(70,nFc,1,dv%coni(1,0,2),'b2tfhi__conix2')
       call my_out_us(70,nFc,1,dv%coni(1,1,0),'b2tfhi__coniy0')
       call my_out_us(70,nFc,1,dv%coni(1,1,1),'b2tfhi__coniy1')
       call my_out_us(70,nFc,1,dv%coni(1,1,2),'b2tfhi__coniy2')
       call my_out_us(70,nFc,1,dv%fhi(1,0),'b2tfhi__fhix')
       call my_out_us(70,nFc,1,dv%fhi(1,1),'b2tfhi__fhiy')
       call my_out_us(70,nFc,1,dv%fhi_mdf(1,0),'b2tfhi__fhi_mdfx')
       call my_out_us(70,nFc,1,dv%fhi_mdf(1,1),'b2tfhi__fhi_mdfy')
       call my_out_us(70,nFc,1,dv%fhiPSch(1,0),'b2tfhi__fhiPSchx')
       call my_out_us(70,nFc,1,dv%fhiPSch(1,1),'b2tfhi__fhiPSchy')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfhi__fchvisparx')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfhi__fchvispary')
       call my_out_us(70,nFc,1,dv%fchvisq(1,0),'b2tfhi__fchvisqx')
       call my_out_us(70,nFc,1,dv%fchvisq(1,1),'b2tfhi__fchvisqy')
       call my_out_us(70,nFc,1,dv%fchvisper(1,0),'b2tfhi__fchvisperx')
       call my_out_us(70,nFc,1,dv%fchvisper(1,1),'b2tfhi__fchvispery')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfhi__fchinertx')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfhi__fchinerty')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfhi__fchanmlx')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfhi__fchanmly')
       call my_out_us(70,nFc,1,dv%fni_he(1,0),'b2tfhi__fnix')
       call my_out_us(70,nFc,1,dv%fni_he(1,1),'b2tfhi__fniy')
       call my_out_us(70,nCv,0,dv%ni(:,0),'b2tfhi__ni0')
       call my_out_us(70,nCv,0,dv%ni(:,1),'b2tfhi__ni1')
c neutrals
       call my_out_us(70,nFc,1,dv%flon(1,0),'b2tfhi__flonx')
       call my_out_us(70,nFc,1,dv%flon(1,1),'b2tfhi__flony')
       call my_out_us(70,nFc,1,dv%conn(1,0,0),'b2tfhi__connx0')
       call my_out_us(70,nFc,1,dv%conn(1,0,1),'b2tfhi__connx1')
       call my_out_us(70,nFc,1,dv%conn(1,0,2),'b2tfhi__connx2')
       call my_out_us(70,nFc,1,dv%conn(1,1,0),'b2tfhi__conny0')
       call my_out_us(70,nFc,1,dv%conn(1,1,1),'b2tfhi__conny1')
       call my_out_us(70,nFc,1,dv%conn(1,1,2),'b2tfhi__conny2')
       call my_out_us(70,nFc,1,dv%fhn(1,0),'b2tfhi__fhnx')
       call my_out_us(70,nFc,1,dv%fhn(1,1),'b2tfhi__fhny')
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

      end
