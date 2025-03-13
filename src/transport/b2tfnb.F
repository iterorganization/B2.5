*-----------------------------------------------------------------------
*.specification

      subroutine b2tfnb (nCv, nFc, nVx, isb, ismain,
     &  switch, geo, mpg, pl, dv, co, rt, balance)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifndef FIXED_POINT
      use b2mod_balance                                                 !djm Jan2017
     & , only : fna_pinch, fna_pll, fna_drift, fna_nanom, fna_panom,
     &          fna_ch, fna_pschused, balance_netcdf
#endif
      use b2mod_ad
     & , only : ncall_b2tfnb
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : b2tfnb_cutlo
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl
      type(B2Rates), intent (in) :: rt
      logical balance                                                   !djm Jan2017
*   ..input/output arguments
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFNB computes the particle flux of one atomic species between
*     neighbouring cells, and the coefficients that specify its
*     linearisation.
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle flux, fnb, of one atomic
*     species, and the coefficients flob and conb that linearise the
*     response of the particle flux to changes in the pressure.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv1, iCv2, iFc, meth                                      !srv 13.10.06
      character chns*3                                                   !srv 16.09.02
      real (kind=R8) :: drift_hyb
      real (kind=R8) ::
     *   flo(nFc,0:1), flo_mdf(nFc,0:1), flo_nodrift(nFc,0:1),
     *   flo2dia(nFc,0:1), flo_53(nFc,0:1),
     *   alpha, gamma, vbar, t0, t1, t2, t3, maxfluxlimit(0:1),
     *   term(nFc), flob0(nFc,0:1), conb0(nFc,0:1), scur(nFc,0:1)
      real (kind=R8) ::
     &   kmprti(nCv), kmprtif(nFc), kmprtiv(nVx),                        !srv 21.10.10
     &   kmpr(nCv), kmprv(nVx),
     &   tbfl_fnb(1:nbc), div_fnb, tbfl_fnb_mdf(1:nbc), div_fnb_mdf,     !srv 05.03.16
     &   rdiv                                                            !srv 05.03.16
      real (kind=R8) ::
     *   wrk0(nFc), wrk1(nFc), wrk2(nFc), wrk3(nFc), wrk4(nFc),
     *   wrk6(nFc), wrk8(nFc),
!    *   wrk5(nFc), wrk7(nFc,0:1),
     *   wrk(nCv), wrkf(nFc,0:1), wrkVx(nVx),
     *   zeros(nFc,0:1), dumm1(nFc,0:1), dumm2(nFc,0:1), dumm3(nFc,0:1),
     *   dumm4(nFc,0:1), zeros1(nFc,0:1)
      real (kind=R8) ::
     *   ub2dia(nFc,0:1),
     *   nbcor(nVx),                                                     !srv 01.10.99 !lk 09.06.11
     *   nbti(nCv),
     *   nbzb(nCv), nbzbc(nVx),                                          !srv 13.05.16
     *   dnb(nFc,0:1), dpb(nFc,0:1), dpbc(nCv), dpbm(nFc),
     *   dpo(nFc,0:1), dnbti(nFc,0:1),
     *   nbf(nFc),
     *   dOnedBsq(nFc,0:1),
     *   fna_mdf0(nFc,0:1), flo_factor(nFc,0:1),
     *   cddb(nFc,0:1), cdpb(nFc,0:1)
      real (kind=R8) :: weight(nFc,2)                                    !lk
      integer :: xcount
*   ..procedures
      external xertst, ipgeti, ipgetr                                    !xpb
      external b2xvsg, intface, diff
      intrinsic abs, max, sqrt, sign
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfnb')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..set internal parameters on first call
      drift_hyb = switch%fnb_drift_hyb
*   ..extensive tests on first few calls
      if (ncall_b2tfnb.lt.3) then
*    ..test sign of nb, ne, pb, cdnb, cdpb
       call b2xvsg (nCv, pl % na(:,isb), 1, 'nb', '.gt.')
       call b2xvsg (nCv, dv % ne       , 1, 'ne', '.gt.')
       call b2xvsg (nCv, dv % pa(:,isb), 1, 'pb', '.gt.')
       wrk0(:) = co%cdna(:,0,isb)*geo%fcQalf(:,0)
       call b2xvsg (nFc, wrk0, 1, 'cdna0','.ge.')
       wrk1(:) = co%cdna(:,1,isb)*geo%fcQalf(:,1)
       call b2xvsg (nFc, wrk1, 1, 'cdna1','.ge.')
       wrk0(:) = co%cdpa(:,0,isb)*geo%fcQalf(:,0)
       call b2xvsg (nFc, wrk0, 1, 'cdpa0','.ge.')
       wrk1(:) = co%cdpa(:,1,isb)*geo%fcQalf(:,1)
       call b2xvsg (nFc, wrk1, 1, 'cdna1','.ge.')
      endif

      meth = switch%b2tfnb_discr_meth
      alpha = switch%b2tfnb_alpha
      gamma = switch%b2tfnb_gamma

*   ..set dummy array
      zeros = 0.0_R8
      zeros1 = 0.0_R8

*   ..interpolate na and pa to cell faces
      weight = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,weight,pl%na(:,isb),nbf)

*   ..precompute kmpr and its differences
      kmpr = pl % na(:,isb) / dv % pa(:,isb)
      call intvertex(nCv,nVx,mpg,geo%vxVol,kmpr,kmprv)

      write(chns,'(i3.3)') isb
      if(switch%mdf_fnb .ne. 0 .and.
     &  (maxval(dv % facdrift)   .ne. 0.0_R8  .or.
     &   maxval(dv % fac_ExB )   .ne. 0.0_R8) .and.
     &   minval(rt % rza(:,isb)) .gt. 0.0_R8) then
*   ..compute kmpr*ti on the faces and vertices, and differences
        kmprti = kmpr*pl%ti
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,kmprti,kmprtif)
        call intvertex(nCv,nVx,mpg,geo%vxVol,kmprti,kmprtiv)
      else
        kmprti  = 0.0_R8
        kmprtif = 0.0_R8
        kmprtiv = 0.0_R8
      endif                                                              !srv 21.10.10 }

c Drift terms
      maxfluxlimit        = 0.0_R8
      dv%uadia(:,:,isb)   = 0.0_R8
      dv%vadia(:,:,isb)   = 0.0_R8
      dv%wadia(:,:,isb)   = 0.0_R8
      dv%vaecrb(:,:,isb)  = 0.0_R8
      dv%fnaPSch(:,:,isb) = 0.0_R8                                       !srv 13.10.06
      dv%fna_mdf(:,:,isb) = 0.0_R8                                       !srv 19.11.09
      dv%flob  = 0.0_R8
      dv%conb  = 0.0_R8                                                  !srv 13.10.06
      ub2dia   = 0.0_R8                                                  !srv 01.10.99
      if((maxval(dv%facdrift).ne.0.0_R8 .or.
     &    maxval(dv%fac_ExB) .ne.0.0_R8) .and.
     &    minval(rt%rza(:,isb)).gt.0.0_R8) then
cwdk    in the unstructured solver, we only support drifts defined
cwdk    on cell faces (drift_style = 1)
        nbti = pl%na(:,isb)*pl%ti
        nbzb = pl%na(:,isb)*rt%rza(:,isb)
*     ..precomputation of gradients on faces
cwdk    to be checked: now we use inverse volume weighting
cwdk                   (original code: simple averaging)
        call intvertex(nCv,nVx,mpg,geo%vxVol,pl%na(:,isb),nbcor)
        call intvertex(nCv,nVx,mpg,geo%vxVol,nbzb,nbzbc)
        call grad(nCv, nFc, nVx, 0, geo, mpg, pl%po, wrkVx, dpo)
        call grad(nCv, nFc, nVx, 0, geo, mpg, nbti, wrkVx, dnbti)
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &   geo%vxOnedBsq, dOnedBsq)
*     ..computation of ExB drift
        wrk0 = geo%fcBb(:,2)/geo%fcBb(:,3)**2
        dv%vaecrb(:,0,isb) = -switch%poleldr*dv%fac_ExB*wrk0*dpo(:,1)
        dv%vaecrb(:,1,isb) =  switch%radeldr*dv%fac_ExB*wrk0*dpo(:,0)
*     ..restriction ExB drift
        xcount = 0
        if (switch%ion_vlct_restrict.eq.2) then                          !iyv 06.02.14 {
          do iFc = 1, nFc
            !wdk for unstructured grids: apply ExB velocity restriction
            !    to total ExB velocity (instead of individual poloidal and radial
            !    components). Otherwise: not perpendicular to gradient of potential
            t0 = sqrt((dv%vaecrb(iFc,0,isb)*geo%fcBb(iFc,3)/
     .        geo%fcBb(iFc,2))**2 + dv%vaecrb(iFc,1,isb)**2)
            if(t0.gt.2.0_R8*co%cssb(iFc)) then
              xcount = xcount+1
              if (switch%vlct_diagno.ge.2) then
                write(6,'(a,a,i3,a,i4,a)')
     .            'Total ExB velocity restriction applied ',
     .            'to species ',isb,' on face (',iFc,'):'
                write(6,'(a,1p,1e12.5,a,1e12.5,a,1e12.5,a,1e12.5)')
     .            'Mach 2: ',2.0_R8*co%cssb(iFc),
     .            ' Vx_ExB:',dv%vaecrb(iFc,0,isb),
     .            ' Vy_ExB:',dv%vaecrb(iFc,1,isb),
     .            ' |V_ExB|:',t0
              end if
              t1 = 2.0_R8*co%cssb(iFc)/t0
              dv%vaecrb(iFc,0,isb) = dv%vaecrb(iFc,0,isb)*t1
              dv%vaecrb(iFc,1,isb) = dv%vaecrb(iFc,1,isb)*t1
            endif
          enddo
        endif
        if (switch%vlct_diagno.ge.1) then
          if (xcount.gt.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Total ExB velocity limit of Mach 2 was applied ',
     .         xcount,' times for species ',isb
          else if (xcount.eq.1) then
            write(6,'(a,a,i6,a,i3)') 'In b2tfnb:',
     .       ' Total ExB velocity limit of Mach 2 was applied ',
     .         xcount,' time for species ',isb
          end if
        end if
*     ..computation of grad-B drift
cwdk    use (na*Ti)_f/na_f for consistency with fchdia
cwdk    check treatment of rza for consistency with other drifts
        wrk = pl%na(:,isb)*pl%ti/(qe*rt%rza(:,isb))
        call intface(nCv,nFc,mpg%fcCv,weight,wrk,wrk1)
cwdk    to be checked: interpolation of geometric factors (hy1)
cwdk                   now consistent between x and y faces
cwdk                   (and consistent with treatment of other drifts)
cwdk                   (implicitly through grad-routine)
        dv%vadia(:,0,isb) =
     &    switch%poldidr*
     &    dv%facdrift*wrk1*dOnedBsq(:,1)*geo%fcBb(:,2)/nbf
        dv%vadia(:,1,isb) =
     &   -switch%raddidr*
     &    dv%facdrift*wrk1*dOnedBsq(:,0)*geo%fcBb(:,2)/nbf
*     ..computation of diamagnetic drift
cwdk    to be checked: now we use inverse volume weighting in wkr0
cwdk                   (original code: simple averaging)
        call intface(nVx,nFc,mpg%fcVx,weight,nbzbc,wrk1)
        dv%wadia(:,0,isb) = -dv%facdrift*wrk0*dnbti(:,1)/(wrk1*qe)
        dv%wadia(:,1,isb) =  dv%facdrift*wrk0*dnbti(:,0)/(wrk1*qe)
*     ..additional drift variables
        dv%uadia(:,:,isb) = dv%vaecrb(:,:,isb) + dv%vadia(:,:,isb)
        ub2dia = dv%vaecrb(:,:,isb) + 2.0_R8*dv%vadia(:,:,isb)
*     ..computation of P.Sch flows
        wrk  = 1.0_R8/(qe*rt%rza(:,isb))
        call intface(nCv,nFc,mpg%fcCv,weight,wrk,wrk1)
        wrk1 = wrk1*geo%fcBzb
        call intface(nVx,nFc,mpg%fcVx,weight,nbcor,wrk2)
        dv%fnaPSch(:,0,isb) =  switch%PSch*dv%facdrift*wrk1*dnbti(:,1)*
     &     geo%fcS*geo%fcQalf(:,0)!*nbf/wrk2
        dv%fnaPSch(:,1,isb) = -switch%PSch*dv%facdrift*wrk1*dnbti(:,0)*
     &     geo%fcS*geo%fcQalf(:,1)!*nbf/wrk2
*     ..linearization coefficients for P.Sch flows
        cddb(:,0) = -switch%PSch*dv%facdrift*co%cddi(:,1,isb)!*nbf/wrk2
        cddb(:,1) =  switch%PSch*dv%facdrift*co%cddi(:,0,isb)!*nbf/wrk2
        call calccoef(nCv,nFc,nVx,0,geo,zeros,cddb,flob0,conb0)
        do iFc = 1, nFc
          dv%flob(iFc,0)   = 0.5_R8*flob0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,1))+kmprti(mpg%fcCv(iFc,2))) -
     &       conb0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,2))-kmprti(mpg%fcCv(iFc,1)))
          dv%flob(iFc,1)   = 0.5_R8*flob0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,1))+kmprtiv(mpg%fcVx(iFc,2))) -
     &       conb0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,2))-kmprtiv(mpg%fcVx(iFc,1)))
!wdk      For now, directly add these contributions into conb(,,0),
!wdk      so they are directly added into the stencil.
!wdk      Might need to move back to dv%conb(,,3) and dedicated treatment
!wdk      in b2ux9p if we decide to move away from inverse volume weighting
!wdk      to compute vertex values of drifts!
          dv%conb(iFc,0,0) = 0.5_R8*conb0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,1))+kmprti(mpg%fcCv(iFc,2))) -
     &       0.25_R8*flob0(iFc,0)*
     &      (kmprti(mpg%fcCv(iFc,2))-kmprti(mpg%fcCv(iFc,1)))
          dv%conb(iFc,1,0) = 0.5_R8*conb0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,1))+kmprtiv(mpg%fcVx(iFc,2))) -
     &       0.25_R8*flob0(iFc,1)*
     &      (kmprtiv(mpg%fcVx(iFc,2))-kmprtiv(mpg%fcVx(iFc,1)))
        enddo
!WG_TODO        dv%conb(:,:,4) = dv%conb(:,:,3) ! Not needed anymore?
        if (isb.ne.ismain .and. .not.is_neutral(isb)) then               !lk 02.06.11 {
          do iFc = 1, nFc
            !wdk for unstructured grids: apply flux limits
            !    to total drift velocities (instead of individual poloidal and radial
            !    components). Otherwise: not perpendicular to gradients
            iCv1 = mpg%fcCv(iFc,1)
            iCv2 = mpg%fcCv(iFc,2)
            t0 = max(pl%ti(iCv1), pl%ti(iCv2))
            t1 = 0.5_R8*(pl%na(iCv1,isb) + pl%na(iCv2,isb))
            vbar = geo%fcBb(iFc,0)/geo%fcBb(iFc,3)*
     &        sqrt(8.0_R8*t0/(pi*am(isb)*mp))
            t2 = 1.0_R8 + sqrt((switch%PSch*dv%facdrift(iFc)*wrk1(iFc)*
     .        dnbti(iFc,1)*geo%fcBb(iFc,3)/geo%fcBb(iFc,2))**2 +
     .        (switch%PSch*dv%facdrift(iFc)*wrk1(iFc)*dnbti(iFc,0))**2)/
     .        (vbar*t1)
            dv%fnaPSch(iFc,0,isb) = dv%fnaPSch(iFc,0,isb)/t2
            dv%fnaPSch(iFc,1,isb) = dv%fnaPSch(iFc,1,isb)/t2
            t3 = sqrt((dv%vadia(iFc,0,isb)*geo%fcBb(iFc,3)/
     .        geo%fcBb(iFc,2))**2 + dv%vadia(iFc,1,isb)**2)
            dv%vadia(iFc,0,isb) = dv%vadia(iFc,0,isb)/(1.0_R8 + t3/vbar)
            dv%vadia(iFc,1,isb) = dv%vadia(iFc,1,isb)/(1.0_R8 + t3/vbar)
            t3 = sqrt((dv%wadia(iFc,0,isb)*geo%fcBb(iFc,3)/
     .        geo%fcBb(iFc,2))**2 + dv%wadia(iFc,1,isb)**2)
            dv%wadia(iFc,0,isb) = dv%wadia(iFc,0,isb)/(1.0_R8 + t3/vbar)
            dv%wadia(iFc,1,isb) = dv%wadia(iFc,1,isb)/(1.0_R8 + t3/vbar)
            t3 = sqrt((dv%uadia(iFc,0,isb)*geo%fcBb(iFc,3)/
     .        geo%fcBb(iFc,2))**2 + dv%uadia(iFc,1,isb)**2)
            dv%uadia(iFc,0,isb) = dv%uadia(iFc,0,isb)/(1.0_R8 + t3/vbar)
            dv%uadia(iFc,1,isb) = dv%uadia(iFc,1,isb)/(1.0_R8 + t3/vbar)
            t3 = sqrt((ub2dia(iFc,0)*geo%fcBb(iFc,3)/
     .        geo%fcBb(iFc,2))**2 + ub2dia(iFc,1)**2)
            ub2dia(iFc,0) = ub2dia(iFc,0)/(1.0_R8 + t3/vbar)
            ub2dia(iFc,1) = ub2dia(iFc,1)/(1.0_R8 + t3/vbar)
          enddo
        endif                                                            !lk 02.06.11 }
        if (isb.eq.ismain) then                                          !lkw 25.10.2022 {
          do iFc = 1, nFc
            iCv1 = mpg%fcCv(iFc,1)
            iCv2 = mpg%fcCv(iFc,2)
            t0 = max(pl%ti(iCv1), pl%ti(iCv2))
            t1 = 0.5_R8*(pl%na(iCv1,isb) + pl%na(iCv2,isb))
            vbar =   sqrt(8.0_R8*t0/(pi*am(isb)*mp))
            t2 = 1.0_R8 + sqrt((switch%PSch*dv%facdrift(iFc)*wrk1(iFc)*
     .        dnbti(iFc,1)*geo%fcBb(iFc,3)/geo%fcBb(iFc,2))**2 +
     .        (switch%PSch*dv%facdrift(iFc)*wrk1(iFc)*dnbti(iFc,0))**2)/
     .        (vbar*t1)
            dv%fnaPSch(iFc,0,isb) = dv%fnaPSch(iFc,0,isb)/t2
            dv%fnaPSch(iFc,1,isb) = dv%fnaPSch(iFc,1,isb)/t2
           enddo
         endif                                                           !lkw 25.10.2022}
      endif

*   ..compute differences of pa at cell faces and poloidal gradients
*     at cell centers
      call diff(nCv, nFc, nVx, 0, geo, mpg, dv%pa(:,isb), wrkVx, dpb)
      call gradc_p(nCv, nFc, nVx, 1, geo, mpg, dv%pa(:,isb), wrkVx,
     &  dpbc)

* ..poloidal fluxes

*   ..interpolate ua to cell faces
      if(switch%b2tfnb_style_int_vel.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,pl%ua(:,isb),wrk0)
      else
        call intface(nCv,nFc,mpg%fcCv,geo%fcHc,pl%ua(:,isb),wrk0)
      endif

*   ..compute convective coefficients (poloidal and radial)
      flo(:,0) = geo%fcPbs*wrk0 + co%cvla(:,0,isb) +
     &  (drift_hyb*dv%vaecrb(:,0,isb) + dv%vadia(:,0,isb))*
     &  geo%fcS*geo%fcQalf(:,0)
      flo(:,1) = co%cvla(:,1,isb) +
     &  (drift_hyb*dv%vaecrb(:,1,isb) + dv%vadia(:,1,isb))*
     &  geo%fcS*geo%fcQalf(:,1)
      flo_mdf(:,0) = geo%fcPbs*wrk0 + co%cvla(:,0,isb) +
     &  drift_hyb*dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)
      flo_mdf(:,1) = co%cvla(:,1,isb) +
     &  drift_hyb*dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)
      flo_53(:,0) = geo%fcPbs*wrk0 +
     &  dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)
      flo_53(:,1) =
     &  dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)
      flo2dia(:,0) = geo%fcPbshz*wrk0 + co%cvlahz(:,0,isb) +
     &  (drift_hyb*dv%vaecrb(:,0,isb) + 2.0_R8*dv%vadia(:,0,isb))*
     &  geo%fcS*geo%fcQalf(:,0)*geo%fcHz
      flo2dia(:,1) = co%cvlahz(:,1,isb) +
     &  (drift_hyb*dv%vaecrb(:,1,isb) + 2.0_R8*dv%vadia(:,1,isb))*
     &  geo%fcS*geo%fcQalf(:,1)*geo%fcHz
      flo_nodrift(:,0) = co%cvla(:,0,isb) + geo%fcPbs*wrk0
      flo_nodrift(:,1) = co%cvla(:,1,isb)

*   ..add contribution from currents
      if (nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
        scur(:,0) =
     &     (dv%fchinert(:,0) + dv%fchvispar(:,0) + dv%fchvisper(:,0) +
     &      dv%fchanml(:,0)  + dv%fchvisq(:,0))*switch%b2tfnb_xcur/
     &     (qe*nbf) + dv%fchviskt(:,0)/(qe*nbf)
        scur(:,1) =
     &     (dv%fchinert(:,1) + dv%fchvispar(:,1) +
     &      switch%fnb_vis_per*dv%fchvisper(:,1) + dv%fchanml(:,1) +
     &      switch%fnb_vis_q*dv%fchvisq(:,1))*switch%b2tfnb_ycur/
     &     (qe*nbf) + dv%fchviskt(:,1)/(qe*nbf)
        flo     = flo     + drift_hyb*scur
        flo_mdf = flo_mdf + drift_hyb*scur
        flo2dia(:,0) = flo2dia(:,0) + drift_hyb*scur(:,0)*geo%fcHz
        flo2dia(:,1) = flo2dia(:,1) + drift_hyb*scur(:,1)*geo%fcHz
      else
        scur = 0.0_R8
      end if

*   ..apply discretization scheme for density-part of fna_mdf
      call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &     flo_mdf, co%cdna(:,:,isb), fna_mdf0(:,:),
     &     dv%fna_32(:,:,isb), dv%fna_52(:,:,isb))
*   ..limit flow by sound speed                                           ! lkw 26.10.2022{
       do iFc = 1, nFc
         if (.not.is_neutral(isb).or.switch%tn_style.eq.0.or.
     &      (switch%tn_style.eq.2.and.nint(zn(isb)).ne.1)) then
           t0 = max(pl%ti(mpg%fcCv(iFc,1)), pl%ti(mpg%fcCv(iFc,2)))
         else if (switch%tn_style.eq.1) then
           t0 = max(pl%tn(mpg%fcCv(iFc,1)), pl%tn(mpg%fcCv(iFc,2)))
         else if (switch%tn_style.eq.2.and.nint(zn(isb)).eq.1) then
           t0 = max(pl%tn(mpg%fcCv(iFc,1)), pl%tn(mpg%fcCv(iFc,2)))
         end if
         t3 = fna_mdf0(iFc,0) + fna_mdf0(iFc,1)
         t1 = pl%na(mpg%fcCv(iFc,1),isb)
         if (t3.lt.0.0_R8) t1 = pl%na(mpg%fcCv(iFc,2),isb)
         t3 = abs(t3)
         vbar = sqrt(8.0_R8 * t0 / (pi * am(isb)*mp)) *
     &        geo%fcS(iFc)
         wrk4(iFc) = 1.0_R8/(1.0_R8 +
     &        (t3/(0.5_R8*vbar*t1))**gamma)**(1.0_R8/gamma)
         fna_mdf0(iFc,0)   = fna_mdf0(iFc,0)*wrk4(iFc)
         fna_mdf0(iFc,1)   = fna_mdf0(iFc,1)*wrk4(iFc)
       enddo                                                             ! lkw 26.10.2022}

!   ..keep the convective flow "correction" factor for further output
      do iFc = 1, nFc
        if (flo_mdf(iFc,0) == 0.0_R8 .and.
     &        dv%fna_32(iFc,0,isb) == 0.0_R8) then
          flo_factor(iFc,0) = 1.0_R8
        else
          call xertst (flo_mdf(iFc,0) /= 0.0_R8, 'faulty transformation'
     &     // ' of poloidal convective particle flow, is = '//chns)
          flo_factor(iFc,0) = dv%fna_32(iFc,0,isb) /
     &              (flo_mdf(iFc,0) * nbf(iFc))
        endif
        if (flo_mdf(iFc,1) == 0.0_R8 .and.
     &      dv%fna_32(iFc,1,isb) == 0.0_R8) then
          flo_factor(iFc,1) = 1.0_R8
        else
          call xertst (flo_mdf(iFc,1) /= 0.0_R8, 'faulty transformation'
     &     // ' of radial convective particle flow, is = '//chns)
          flo_factor(iFc,1) = dv%fna_32(iFc,1,isb) /
     &              (flo_mdf(iFc,1) * nbf(iFc))
        endif
      enddo

*   ..compute Rhie-Chow term
      do iFc = 1, nFc
        cdpb(iFc,0) = switch%b2npco_pcm1*0.5_R8*geo%fcPbs(iFc)*
     &   (dv%pcca(mpg%fcCv(iFc,1),1,isb)*pl%na(mpg%fcCv(iFc,1),isb) +
     &    dv%pcca(mpg%fcCv(iFc,2),1,isb)*pl%na(mpg%fcCv(iFc,2),isb))/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
        dpbm(iFc) = 0.5_R8*geo%fcQgam(iFc,0)*
     &   (dpbc(mpg%fcCv(iFc,1)) + dpbc(mpg%fcCv(iFc,2)))*
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))
      enddo
      cdpb(:,1) = 0.0_R8
      term   = -cdpb(:,0)*(dpb(:,0) - dpbm)
      fna_mdf0(:,0)      = fna_mdf0(:,0) + term
      dv%fna_32(:,0,isb) = dv%fna_32(:,0,isb) + term

*   ..compute other fluxes
      dv%fna(:,0,isb)      = fna_mdf0(:,0) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf +
     &  (flo(:,0) - flo_mdf(:,0))*nbf - co%cdpa(:,0,isb)*dpb(:,0)

      dv%fna(:,1,isb)      = fna_mdf0(:,1) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf +
     &  (flo(:,1) - flo_mdf(:,1))*nbf - co%cdpa(:,1,isb)*dpb(:,1)

      dv%fna_mdf(:,0,isb)  = fna_mdf0(:,0) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf -
     &   co%cdpa(:,0,isb)*dpb(:,0) - dv%fnaPSch(:,0,isb)
      dv%fna_mdf(:,1,isb)  = fna_mdf0(:,1) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf -
     &   co%cdpa(:,1,isb)*dpb(:,1) - dv%fnaPSch(:,1,isb)

      dv%fna_fcor(:,0,isb) = fna_mdf0(:,0)*geo%fcHz +
     &  (1.0_R8-drift_hyb)*geo%fcHz*
     &  (dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf +
     &  (flo2dia(:,0) - flo_mdf(:,0)*geo%fcHz)*nbf -
     &   co%cdpahz(:,0,isb)*dpb(:,0)
      dv%fna_fcor(:,1,isb) = fna_mdf0(:,1)*geo%fcHz +
     &  (1.0_R8-drift_hyb)*geo%fcHz*
     &  (dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf +
     &  (flo2dia(:,1) - flo_mdf(:,1)*geo%fcHz)*nbf -
     &   co%cdpahz(:,1,isb)*dpb(:,1)

      dv%fna_he(:,0,isb)   = fna_mdf0(:,0) + (1.0_R8-drift_hyb)*
     &   dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)*nbf -
     &   co%cdpa(:,0,isb)*dpb(:,0) - drift_hyb*scur(:,0)*nbf
      dv%fna_he(:,1,isb)   = fna_mdf0(:,1) + (1.0_R8-drift_hyb)*
     &   dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)*nbf -
     &   co%cdpa(:,1,isb)*dpb(:,1) - drift_hyb*scur(:,1)*nbf

      dv%fna_53(:,0,isb)   = flo_53(:,0)*nbf
      dv%fna_53(:,1,isb)   = flo_53(:,1)*nbf

      dv%fna_32(:,0,isb)   = dv%fna_32(:,0,isb) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0) + scur(:,0))*nbf +
     &  (flo(:,0) - flo_mdf(:,0))*nbf
      dv%fna_32(:,1,isb)   = dv%fna_32(:,1,isb) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1) + scur(:,1))*nbf +
     &  (flo(:,1) - flo_mdf(:,1))*nbf

      dv%fna_52(:,0,isb)   = dv%fna_52(:,0,isb) -
     &   co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_52(:,1,isb)   = dv%fna_52(:,1,isb) -
     &   co%cdpa(:,1,isb)*dpb(:,1)

*   ..apply discretization scheme for density-part of fna_nodrift
      call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &     flo_nodrift, co%cdna(:,:,isb), dv%fna_nodrift(:,:,isb),
     &     dv%fna_32nd(:,0,isb), dv%fna_52nd(:,0,isb))

*   ..compute _nodrift fluxes
      dv%fna_nodrift(:,0,isb) = dv%fna_nodrift(:,0,isb) -
     &   co%cdpa(:,0,isb)*dpb(:,0) + term
      dv%fna_nodrift(:,1,isb) = dv%fna_nodrift(:,1,isb) -
     &   co%cdpa(:,1,isb)*dpb(:,1)

      dv%fna_32nd(:,0,isb) = dv%fna_32nd(:,0,isb) + term

      dv%fna_52nd(:,0,isb) = dv%fna_52nd(:,0,isb) -
     &   co%cdpa(:,0,isb)*dpb(:,0)
      dv%fna_52nd(:,1,isb) = dv%fna_52nd(:,1,isb) -
     &   co%cdpa(:,1,isb)*dpb(:,1)

*   ..compute ExB fluxes for k-model
      dv%fna_ExB = 0.0_R8
      if (isb .eq. ismain .and. switch%solve_keps.gt.0) then
        zeros = 0.0_R8
        zeros1 = 0.0_R8
        call calcflow(nCv,nFc,nVx,meth,geo,mpg,pl%na(:,isb),zeros,
     &    co%cdna_ExB,dv%fna_ExB,zeros1,dumm1)
      endif

*   ..compute flob, conb
      call calccoef(nCv,nFc,nVx,meth,geo,
     &    flo_mdf,co%cdna(:,:,isb),flob0,conb0)
      flob0(:,0) = flob0(:,0) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0) + scur(:,0))
      flob0(:,1) = flob0(:,1) + (1.0_R8-drift_hyb)*
     &  (dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1) + scur(:,1))
      if (switch%mdf_fnb .eq. 0) then
        flob0 = flob0 + flo - flo_mdf
      endif
      do iFc = 1, nFc
        dv%flob(iFc,0)   = dv%flob(iFc,0) +
     &     0.5_R8*flob0(iFc,0)*
     &     (kmpr(mpg%fcCv(iFc,1))+kmpr(mpg%fcCv(iFc,2))) -
     &     conb0(iFc,0)*(kmpr(mpg%fcCv(iFc,2))-kmpr(mpg%fcCv(iFc,1)))
        dv%flob(iFc,1)   = dv%flob(iFc,1) +
     &     0.5_R8*flob0(iFc,1)*
     &     (kmprv(mpg%fcVx(iFc,1))+kmprv(mpg%fcVx(iFc,2))) -
     &     conb0(iFc,1)*(kmprv(mpg%fcVx(iFc,2))-kmprv(mpg%fcVx(iFc,1)))
        dv%conb(iFc,0,0) = dv%conb(iFc,0,0) +
     &     0.5_R8*conb0(iFc,0)*
     &     (kmpr(mpg%fcCv(iFc,1))+kmpr(mpg%fcCv(iFc,2))) - 0.25_R8*
     &     flob0(iFc,0)*(kmpr(mpg%fcCv(iFc,2))-kmpr(mpg%fcCv(iFc,1)))
        dv%conb(iFc,1,0) = dv%conb(iFc,1,0) +
     &     0.5_R8*conb0(iFc,1)*
     &     (kmprv(mpg%fcVx(iFc,1))+kmprv(mpg%fcVx(iFc,2))) - 0.25_R8*
     &     flob0(iFc,1)*(kmprv(mpg%fcVx(iFc,2))-kmprv(mpg%fcVx(iFc,1)))
      enddo

*   ..contribution from pressure diffusion and Rhie-Chow term
      cdpb = cdpb + co%cdpa(:,:,isb)
      zeros = 0.0_R8
      call calccoef(nCv,nFc,nVx,0,geo,zeros,cdpb,flob0,conb0)
      dv%flob = dv%flob + flob0                                     !wdk flob0 should be zero here
      dv%conb(:,:,0) = dv%conb(:,:,0) + conb0
      dv%conb(:,0,1) = conb0(:,0)*kmprtif                           !wdk is this needed? correct? probably always zero?
      dv%conb(:,1,1) = conb0(:,1)*kmprtif                           !wdk is this needed? correct? probably always zero?
      dv%conb(:,:,2) = dv%conb(:,:,1)

*   ..flux limit
      if ((alpha.gt.0.0_R8 .and. maxval(rt%rza(:,isb)).eq.0.0_R8) .and.  !srv 20.11.08
     *     switch%b2tlc0_alpha.eq.0.0_R8 ) then                          ! neutral flux limit
       call xerrab('Flux limit to be adapted for WG')
cwdk try to introduce a common flux limit here?
       do iFc = 1,nFc
        t2 = max(switch%b2tfnb_flux_limit_min_ti,
     *        pl%ti(mpg%fcCv(iFc,1)), pl%ti(mpg%fcCv(iFc,2)))
        vbar = sqrt(8.0_R8 * t2 / (pi * am(isb)*mp))
        t0 = (1.0_R8 + (abs(dv%fna(iFc,0,isb) / (0.5_R8*
     &   (pl%na(mpg%fcCv(iFc,1),isb) + pl%na(mpg%fcCv(iFc,2),isb))*
     &    alpha*vbar/4.0_R8*geo%fcS(iFc)*geo%fcQalf(iFc,0))))**gamma)**
     &   (1.0_R8/gamma)
        t1 = (1.0_R8 + (abs(dv%fna(iFc,1,isb) / (0.5_R8*
     &   (pl%na(mpg%fcCv(iFc,1),isb) + pl%na(mpg%fcCv(iFc,2),isb))*
     &    alpha*vbar/4.0_R8*geo%fcS(iFc)*geo%fcQalf(iFc,1))))**gamma)**
     &   (1.0_R8/gamma)
*     ..poloidal direction
        maxfluxlimit(0)=max(maxfluxlimit(0),t0)
        dv%fna(iFc,0,isb) = dv%fna(iFc,0,isb) / t0
        dv%fna_32(iFc,0,isb) = dv%fna_32(iFc,0,isb) / t0
        dv%fna_52(iFc,0,isb) = dv%fna_52(iFc,0,isb) / t0
        dv%fna_nodrift(iFc,0,isb) = dv%fna_nodrift(iFc,0,isb) / t0
        dv%fna_32nd(iFc,0,isb) = dv%fna_32nd(iFc,0,isb) / t0
        dv%fna_52nd(iFc,0,isb) = dv%fna_52nd(iFc,0,isb) / t0
        dv%flob(iFc,0) = dv%flob(iFc,0) / t0                             ! to check: flob(:,0), conb(:,0,:) refers to connector, not poloidal... do flux limit before computation of flob/conb?
        dv%conb(iFc,0,0:4) = dv%conb(iFc,0,0:4) / t0                     !srv 13.10.06 {
        dv%fna_he(iFc,0,isb)   = dv%fna_he(iFc,0,isb)/t0
        dv%fna_fcor(iFc,0,isb) = dv%fna_fcor(iFc,0,isb)/t0
        dv%fna_mdf(iFc,0,isb) = dv%fna_mdf(iFc,0,isb)/t0
        dv%fnaPSch(iFc,0,isb)  = dv%fnaPSch(iFc,0,isb)/t0                !srv 13.10.06 }
        co%fllim0fna(iFc,0,isb) = t0
*     ..radial direction
        maxfluxlimit(1)=max(maxfluxlimit(1),t1)
        dv%fna(iFc,1,isb) = dv%fna(iFc,1,isb) / t1
        dv%fna_32(iFc,1,isb) = dv%fna_32(iFc,1,isb) / t1
        dv%fna_52(iFc,1,isb) = dv%fna_52(iFc,1,isb) / t1
        dv%fna_nodrift(iFc,1,isb) = dv%fna_nodrift(iFc,1,isb) / t1
        dv%fna_32nd(iFc,1,isb) = dv%fna_32nd(iFc,1,isb) / t1
        dv%fna_52nd(iFc,1,isb) = dv%fna_52nd(iFc,1,isb) / t1
        dv%flob(iFc,1) = dv%flob(iFc,1) / t1                             ! to check: flob(:,1), conb(:,1,:) refers to tangential, not radial...
        dv%conb(iFc,1,0:4) = dv%conb(iFc,1,0:4) / t1                     !srv 13.10.06 {
        dv%fna_he(iFc,1,isb)   = dv%fna_he(iFc,1,isb)/t1
        dv%fna_fcor(iFc,1,isb) = dv%fna_fcor(iFc,1,isb)/t1
        dv%fna_mdf(iFc,1,isb) = dv%fna_mdf(iFc,1,isb)/t1
        dv%fnaPSch(iFc,1,isb)  = dv%fnaPSch(iFc,1,isb)/t1                !srv 13.10.06 }
        co%fllim0fna(iFc,1,isb) = t1
       enddo
      endif

      if(switch%b2npco_iout.ne.0
     &                     .or. switch%b2tfnb_iout .eq. 1 .or.
     &                               switch%iout_b2wdat.eq.4) then       !srv 11.09.09 {
cwdk to be considered: con_mdf, con_2dia not easily recoverable?
cwdk should this be conb in general? (but without the rescaling for pressure)
       call diff(nCv, nFc, nVx, 0, geo, mpg, pl%na(:,isb), wrkVx, dnb)
       wrk1 = geo%fcPbs*wrk0*nbf
       wrk2 = term
       wrk3 = dv%vaecrb(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)
     &   * (1.0_R8 + (flo_factor(:,0) - 1.0_R8) * drift_hyb) ! to account possible correction due to hybridization
!       wrk4 = -con_2dia*dnb(:,0)
!       wrk5 = -co%cdpa(:,0,isb)*dpb(:,0)
       wrk6 = term*geo%fcHz
!       wrk7 = -con_mdf*dnb(:,0)
       wrk8 = co%cvla(:,0,isb)*nbf
       if(.not. is_neutral(isb)) then
        call my_out_us(70,nFc,1,wrk3,'b2tfnb_vaecrbna_th'//chns)
        call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &                                'b2tfnb_fnbPSch_th'//chns)
!        call my_out_us(70,nFc,1,wrk7,'b2tfnb_dPat_mdf_gradnax'//chns)
!        call my_out_us(70,nFc,1,wrk4,'b2tfnb_dPat_2diagradnax'//chns)
!       else
!        call my_out_us(70,nFc,1,wrk5,'b2tfnb_dgradpb_th'//chns)
!*****
!        Note that dv%fna_52 already contains conductive terms caused by both pressure and density gradient. It is printed below
!*****
       endif
       call my_out_us(70,nFc,1,wrk1,'b2tfnb_bxuana_th'//chns)
       call my_out_us(70,nFc,1,wrk2,'b2tfnb_dpccorna_th'//chns)
       call my_out_us(70,nFc,1,dpb(1,0),'b2tfnb_dpb_th'//chns)
       call my_out_us(70,nFc,1,dpbm,'b2tfnb_dpbm'//chns)
       if (switch%iout_b2wdat .ne. 4)
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_dpccorhzna_th'//chns)       !srv 05.11.09
       if (switch%b2tfnb_iout.eq.1)
     &  call my_out_us(70,nFc,1,wrk8,'b2tfnb_cvlbna_th'//chns)           !iys 17.04.2012
!       call my_out(70,nx,ny,con_2dia,'b2tfnb_dPat_2diax'//chns)
!       call my_out(70,nx,ny,con_mdf,'b2tfnb_dPat_mdfx'//chns)
      endif

      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4 .or.
     &                                  switch%b2tfnb_iout.eq.1) then
       write(chns,'(i3.3)') isb
       wrk6 = dv%vadia(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)
       if(.not. is_neutral(isb))
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadiana_th'//chns)
       wrk6 = dv%vadia(:,0,isb)*nbf*geo%fcS*geo%fcQalf(:,0)*geo%fcHz
       if(.not. is_neutral(isb) .and. switch%iout_b2wdat.ne.4)
     &   call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianahz_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &              'b2tfnb_fnb_fcor_th'//chns)
      endif                                                              !srv 11.09.09 }

      if (alpha.gt.0.0_R8 .and. maxval(rt%rza(:,isb)).eq.0.0_R8
     > .and. ncall_b2tfnb.lt.100) then
       write(*,'(a,1p,e12.5,a,e12.5,a,2e12.5)')
     .   'b2tfnb: maxfluxlimit for zb = ',maxval(rt%rza(:,isb)),
     1   ' mb = ',am(isb),' are ',maxfluxlimit
      endif

      if(switch%b2npco_iout.ne.0
     & .or. switch%b2tfnb_iout .eq. 1 .or. switch%iout_b2wdat.eq.4) then !srv 11.09.09 {
       wrk3 = dv%vaecrb(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)
     &   * (1.0_R8 + (flo_factor(:,1) - 1.0_R8) * drift_hyb) ! to account possible correction due to hybridization
!       wrk4 = -con_2dia*dnb(:,1)
!       wrk5 = -co%cdpa(:,1,isb)*dpb(:,1)
!       wrk7 = -con_mdf*dnb(:,1)
       wrk8 = co%cvla(:,1,isb)*nbf
       write(chns,'(i3.3)') isb
       if(.not. is_neutral(isb)) then
        call my_out_us(70,nFc,1,dv%fnaPSch(:,1,isb),
     &                                'b2tfnb_fnbPSch_r'//chns)
        call my_out_us(70,nFc,1,wrk3,'b2tfnb_vaecrbna_r'//chns)
!        call my_out_us(70,nFc,1,wrk7,'b2tfnb_dPat_mdf_gradnay'//chns)
!        call my_out_us(70,nFc,1,wrk4,'b2tfnb_dPat_2diagradnay'//chns)
!       else
!        call my_out_us(70,nFc,1,wrk5,'b2tfnb_dgradpb_r'//chns)
!*****
!        Note that dv%fna_52 already contains conductive terms caused by both pressure and density gradient. It is printed below
!*****
       endif
       call my_out_us(70,nFc,1,wrk8,'b2tfnb_cvlbna_r'//chns)               !iys 17.04.2012
!       call my_out_us(70,nFc,1,con_2dia,'b2tfnb_dPat_2diay'//chns)
!       call my_out_us(70,nFc,1,con_mdf,'b2tfnb_dPat_mdfy'//chns)
      endif

      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4
     &                         .or. switch%b2tfnb_iout.eq.1) then
       write(chns,'(i3.3)') isb
       wrk6 = dv%vadia(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)
       if(.not. is_neutral(isb))
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadiana_r'//chns)
       wrk6 = dv%vadia(:,1,isb)*nbf*geo%fcS*geo%fcQalf(:,1)*geo%fcHz
       if(.not. is_neutral(isb) .and. switch%iout_b2wdat.ne.4)
     &  call my_out_us(70,nFc,1,wrk6,'b2tfnb_vadianahz_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &  'b2tfnb_fnb_fcor_r'//chns)
      endif                                                              !srv 11.09.09 }

* ..compute kinetic energy of species
      dv%kinrgy(:,isb) = 0.5_R8*am(isb)*mp*pl%ua(:,isb)**2

*   ..flux correction
*   ..we rescale the fnb fluxes to have the same integral over the boundaries
*   ..as the fnb_mdf fluxes
      if(switch%pflux_cor.ne.0 .and. switch%mdf_fnb.ne.0) then          !srv 05.03.16 {
        if(.not. is_neutral(isb)) then
          call b2tbfl(nCv, nFc, .true., mpg, dv%fna(:,:,isb),
     &      tbfl_fnb, div_fnb)
          call b2tbfl(nCv, nFc, .true., mpg, dv%fna_mdf(:,:,isb),
     &      tbfl_fnb_mdf, div_fnb_mdf)
          if(div_fnb .ne. 0.0_R8) then
            rdiv= div_fnb_mdf/div_fnb
            dv%fna(:,:,isb) = rdiv*dv%fna(:,:,isb)
          else
            if (ncall_b2tfnb.eq.0) write(*,"(a,i4)")
     >       ' Boundary flux correction not attempted for is = ',isb
          endif
          if(switch%b2tfnb_iout.ne.0) then
            write(*,*)
            write(*,*) 'Before flux correction in b2tfnb, is = ', isb
            write(*,'(a,i4,a,es18.11,a,es18.11)')
     $       ' isb = ',isb,' fnb_tot = ',div_fnb,
     $       ' fnb_mdf_tot = ',div_fnb_mdf

            call b2tbfl(nCv, nFc, .true., mpg, dv%fna(:,:,isb),
     &        tbfl_fnb, div_fnb)
            call b2tbfl(nCv, nFc, .true., mpg, dv%fna_mdf(:,:,isb),
     &        tbfl_fnb_mdf, div_fnb_mdf)
            write(*,*) 'After flux correction in b2tfnb, is = ', isb
            write(*,'(a,i4,a,es18.11,a,es18.11)')
     $       ' isb = ',isb,' fnb_tot = ',div_fnb,
     $       ' fnb_mdf_tot = ',div_fnb_mdf
            write(*,*)
!            write(*,'(a,i3,a,es15.7,a,es15.7,a,es15.7)')
!     &       ' is =',isb,' rdiv =',rdiv,' div_fnb =',div_fnb,
!     &       ' div_fnb_mdf =',div_fnb_mdf
          endif
        endif
      endif                                                             !srv 05.03.16 }
#ifndef FIXED_POINT
* ..compute variables for balance.nc
      if (balance.and.(balance_netcdf.ne.0)) then    !nh 08.01.24
        dumm1 = 0.0_R8
        dumm2 = 0.0_R8
        dumm3 = 0.0_R8
        dumm4 = 0.0_R8
        wrkf = sign(co%cvla(:,:,isb),flo_mdf)
        call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &    wrkf, dumm1, dumm2, fna_pinch(:,:,isb), dumm3)
        do iFc = 1,nFc !nh limit by sound speed
          fna_pinch(iFc,:,isb) = wrk4(iFc)*fna_pinch(iFc,:,isb)
        enddo
        fna_pinch(:,:,isb) = sign(fna_pinch(:,:,isb),co%cvla(:,:,isb))
        dumm4(:,0) = geo%fcPbs*wrk0
        wrkf = sign(dumm4,flo_mdf)
        call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &    wrkf, dumm1, dumm2, fna_pll(:,:,isb), dumm3)
        do iFc = 1,nFc !nh limit by sound speed
          fna_pll(iFc,:,isb) = wrk4(iFc)*fna_pll(iFc,:,isb)
        enddo
        fna_pll(:,:,isb) = sign(fna_pll(:,:,isb),dumm4)
        fna_pll(:,0,isb) = fna_pll(:,0,isb) + term
        dumm4(:,0) = drift_hyb*dv%vaecrb(:,0,isb)*geo%fcS*
     &    geo%fcQalf(:,0)
        dumm4(:,1) = drift_hyb*dv%vaecrb(:,1,isb)*geo%fcS*
     &    geo%fcQalf(:,1)
        wrkf = sign(dumm4,flo_mdf)
        call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &    wrkf, dumm1, dumm2, fna_drift(:,:,isb), dumm3)
        do iFc = 1,nFc !nh limit by sound speed
          fna_drift(iFc,:,isb) = wrk4(iFc)*fna_drift(iFc,:,isb)
        enddo
        fna_drift(:,:,isb) = sign(fna_drift(:,:,isb),dumm4)
        fna_drift(:,0,isb) = fna_drift(:,0,isb) + (1.0_R8-drift_hyb)*
     &    dv%vaecrb(:,0,isb)*geo%fcS*geo%fcQalf(:,0)*nbf
        fna_drift(:,1,isb) = fna_drift(:,1,isb) + (1.0_R8-drift_hyb)*
     &    dv%vaecrb(:,1,isb)*geo%fcS*geo%fcQalf(:,1)*nbf
        if (switch%mdf_fnb.eq.0) then
          fna_drift(:,0,isb) = fna_drift(:,0,isb) +
     &      dv%vadia(:,0,isb)*geo%fcS*geo%fcQalf(:,0)*nbf
          fna_drift(:,1,isb) = fna_drift(:,1,isb) +
     &      dv%vadia(:,1,isb)*geo%fcS*geo%fcQalf(:,1)*nbf
          fna_pschused(:,:,isb) = 0.0_R8
        else
          fna_pschused(:,:,isb) = -dv%fnaPSch(:,:,isb)
        endif
        wrkf = sign(drift_hyb*scur,flo_mdf)
        call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &    wrkf, dumm1, dumm2, fna_ch(:,:,isb), dumm3)
        do iFc = 1,nFc !nh limit by sound speed
          fna_ch(iFc,:,isb) = wrk4(iFc)*fna_ch(iFc,:,isb)
        enddo
        fna_ch(:,:,isb) = sign(fna_ch(:,:,isb),drift_hyb*scur)
        fna_ch(:,0,isb) = fna_ch(:,0,isb) + (1.0_R8-drift_hyb)*
     &    scur(:,0)*nbf
        fna_ch(:,1,isb) = fna_ch(:,1,isb) + (1.0_R8-drift_hyb)*
     &    scur(:,1)*nbf
        call calcflow(nCv, nFc, nVx, meth, geo, mpg, pl%na(:,isb),
     &    flo_mdf, co%cdna(:,:,isb), dumm2, dumm3, fna_nanom(:,:,isb))
        do iFc = 1,nFc !nh limit by sound speed
          fna_nanom(iFc,:,isb) = wrk4(iFc)*fna_nanom(iFc,:,isb)
        enddo
        fna_panom(:,:,isb) = -co%cdpa(:,:,isb)*dpb(:,:)
      endif
#endif

      if (switch%iout_b2wdat .eq. 4) then
       call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2tfnb_fnb_th'//chns)
       call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2tfnb_fnb_r'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &            'b2tfnb_fnb_fcor_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &            'b2tfnb_fnb_fcor_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,0,isb),
     &            'b2tfnb_fnb_he_th'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,1,isb),
     &            'b2tfnb_fnb_he_r'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &            'b2tfnb_fnbPSch_th'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,1,isb),
     &            'b2tfnb_fnbPSch_r'//chns)
       if (isb .eq. ismain) then
          call my_out_us(70,nFc,1,dv%vaecrb(1,0,isb),
     &            'b2tfnb_vbecrb_th'//chns)
          call my_out_us(70,nFc,1,dv%vaecrb(1,1,isb),
     &            'b2tfnb_vbecrb_r'//chns)
       endif
!       call my_out_us(70,nFc,1,dv%uadia(1,0,isb),'b2tfnb_ubdia_th'//chns)
!       call my_out_us(70,nFc,1,dv%uadia(1,1,isb),'b2tfnb_ubdia_r'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,0,isb),
     &            'b2tfnb_vbdia_th'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,1,isb),
     &            'b2tfnb_vbdia_r'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,0,isb),
     &            'b2tfnb_wbdia_th'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,1,isb),
     &            'b2tfnb_wbdia_r'//chns)
       call my_out_us(70,nCv,0,dv%kinrgy,'b2tfnb_kbnrgy'//chns)
! convective fluxes (drifts + currents + turbulent convection) with account of hybridization and flux limit
       call my_out_us(70,nFc,1,dv%fna_32(1,0,isb),
     &            'b2tfnb_fnb_32_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_32(1,1,isb),
     &            'b2tfnb_fnb_32_r'//chns)
! conductive fluxes (-D_n \nabla n - D_p \nabla p) with account of flux limit
       call my_out_us(70,nFc,1,dv%fna_52(1,0,isb),
     &            'b2tfnb_fnb_52_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_52(1,1,isb),
     &            'b2tfnb_fnb_52_r'//chns)
      endif

      if(switch%b2tfnb_iout.ne.0) then                                                !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2tfnb_fnb_th'//chns)
       call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2tfnb_fnb_r'//chns)
       call my_out_us(70,nFc,1,dv%fna_nodrift(1,0,isb),
     &             'b2tfnb_fnb_nodrift_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_nodrift(1,1,isb),
     &             'b2tfnb_fnb_nodrift_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,0,isb),
     &             'b2tfnb_fnb_fcor_th'//chns)
      call my_out_us(70,nFc,1,dv%fna_fcor(1,1,isb),
     &             'b2tfnb_fnb_fcor_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,0,isb),
     &            'b2tfnb_fnb_he_th'//chns)
      call my_out_us(70,nFc,1,dv%fna_he(1,1,isb),
     &            'b2tfnb_fnb_he_r'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,0,isb),
     &             'b2tfnb_fnbPSch_th'//chns)
       call my_out_us(70,nFc,1,dv%fnaPSch(1,1,isb),
     &             'b2tfnb_fnbPSch_r'//chns)
       call my_out_us(70,nFc,1,dv%fna_mdf(1,0,isb),
     &             'b2tfnb_fnb_mdf_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_mdf(1,1,isb),
     &             'b2tfnb_fnb_mdf_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_32(1,0,isb),
     &            'b2tfnb_fnb_32_th'//chns)
      call my_out_us(70,nFc,1,dv%fna_32(1,1,isb),
     &            'b2tfnb_fnb_32_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_52(1,0,isb),
     &            'b2tfnb_fnb_52_th'//chns)
      call my_out_us(70,nFc,1,dv%fna_52(1,1,isb),
     &            'b2tfnb_fnb_52_r'//chns)
      call my_out_us(70,nFc,1,dv%fna_32nd(1,0,isb),
     &             'b2tfnb_fnb_32nd_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_32nd(1,1,isb),
     &             'b2tfnb_fnb_32nd_r'//chns)
       call my_out_us(70,nFc,1,dv%fna_52nd(1,0,isb),
     &             'b2tfnb_fnb_52nd_th'//chns)
       call my_out_us(70,nFc,1,dv%fna_52nd(1,1,isb),
     &             'b2tfnb_fnb_52nd_r'//chns)
       call my_out_us(70,nFc,1,co%fllim0fna(1,0,isb),
     &             'b2tfnb_fllim0fnb_th'//chns)
       call my_out_us(70,nFc,1,co%fllim0fna(1,1,isb),
     &             'b2tfnb_fllim0fnb_r'//chns)
      call my_out_us(70,nFc,1,dv%vaecrb(1,0,isb),
     &            'b2tfnb_vbecrb_th'//chns)
      call my_out_us(70,nFc,1,dv%vaecrb(1,1,isb),
     &            'b2tfnb_vbecrb_r'//chns)
       call my_out_us(70,nFc,1,dv%uadia(1,0,isb),
     &            'b2tfnb_ubdia_th'//chns)
       call my_out_us(70,nFc,1,dv%uadia(1,1,isb),
     &            'b2tfnb_ubdia_r'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,0,isb),
     &            'b2tfnb_vbdia_th'//chns)
       call my_out_us(70,nFc,1,dv%vadia(1,1,isb),
     &            'b2tfnb_vbdia_r'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,0,isb),
     &            'b2tfnb_wbdia_th'//chns)
       call my_out_us(70,nFc,1,dv%wadia(1,1,isb),
     &            'b2tfnb_wbdia_r'//chns)
       call my_out_us(70,nCv,0,dv%kinrgy,'b2tfnb_kbnrgy'//chns)
       call my_out_us(70,nFc,1,co%cvla(1,0,isb),
     &            'b2tfnb_cvlb_th'//chns)
       call my_out_us(70,nFc,1,co%cvla(1,1,isb),
     &            'b2tfnb_cvlb_r'//chns)
      call my_out_us(70,nFc,1,co%cvlahz(1,0,isb),
     &            'b2tfnb_cvlbhz_th'//chns)
      call my_out_us(70,nFc,1,co%cvlahz(1,1,isb),
     &            'b2tfnb_cvlbhz_r'//chns)
      call my_out_us(70,nFc,1,co%cdpahz(1,0,isb),
     &            'b2tfnb_cdpbhz_th'//chns)
      call my_out_us(70,nFc,1,co%cdpahz(1,1,isb),
     &            'b2tfnb_cdpbhz_r'//chns)
       call my_out_us(70,nFc,1,co%cdpa(1,0,isb),
     &            'b2tfnb_cdpb_th'//chns)
       call my_out_us(70,nFc,1,co%cdpa(1,1,isb),
     &            'b2tfnb_cdpb_r'//chns)
       call my_out_us(70,nFc,1,co%cdna(1,0,isb),
     &            'b2tfnb_cdnb_th'//chns)
       call my_out_us(70,nFc,1,co%cdna(1,1,isb),
     &            'b2tfnb_cdnb_r'//chns)
       call my_out_us(70,nFc,1,co%cddi(1,0,isb),
     &            'b2tfnb_cddb_th'//chns)
       call my_out_us(70,nFc,1,co%cddi(1,1,isb),
     &            'b2tfnb_cddb_r'//chns)
       call my_out_us(70,nFc,1,term,'b2tfnb_term'//chns)
       call my_out_us(70,nFc,1,term*geo%fcHz,'b2tfnb_termhz'//chns)
       call my_out_us(70,nFc,1,geo%fcPbs,'b2tfnb_pbs')
       call my_out_us(70,nFc,1,geo%fcPbshz,'b2tfnb_pbshz')
       call my_out_us(70,nFc,1,ub2dia(1,0),'b2tfnb_ub2dia_th'//chns)
       call my_out_us(70,nFc,1,ub2dia(1,1),'b2tfnb_ub2dia_r'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,0),'b2tfnb_conb_th0'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,1),'b2tfnb_conb_th1'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,2),'b2tfnb_conb_th2'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,3),'b2tfnb_conb_th3'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,0,4),'b2tfnb_conb_th4'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,0),'b2tfnb_conb_r0'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,1),'b2tfnb_conb_r1'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,2),'b2tfnb_conb_r2'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,3),'b2tfnb_conb_r3'//chns)
       call my_out_us(70,nFc,1,dv%conb(1,1,4),'b2tfnb_conb_r4'//chns)
       call my_out_us(70,nFc,1,dv%flob(1,0),'b2tfnb_flob_th'//chns)
       call my_out_us(70,nFc,1,dv%flob(1,1),'b2tfnb_flob_r'//chns)
       call my_out_us(70,nFc,1,kmprtif(1),'b2tfnb_kmprtix'//chns)
       call my_out_us(70,nCv,0,pl%ua(1,isb),'b2tfnb_ua'//chns)
       call my_out_us(70,nCv,0,pl%na(1,isb),'b2tfnb_na'//chns)
       call my_out_us(70,nCv,0,dv%pa(1,isb),'b2tfnb_pa'//chns)
       call my_out_us(70,nCv,0,pl%te,'b2tfnb_te')
       call my_out_us(70,nCv,0,pl%ti,'b2tfnb_ti')
       call my_out_us(70,nCv,0,pl%po,'b2tfnb_po')
       call my_out_us(70,nCv,0,rt%rza(1,isb),'b2tfnb_zb'//chns)
       call my_out_us(70,nFc,1,dv%facdrift,'b2tfnb_facdrift'//chns)
       call my_out_us(70,nFc,1,dv%fac_ExB,'b2tfnb_fac_ExB'//chns)
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfnb_fchvispar_th')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfnb_fchvispar_r')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfnb_fchinert_th')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfnb_fchinert_r')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfnb_fchanml_th')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfnb_fchanml_r')
      endif                                                             !srv 18.06.08 }

* ..return
      ncall_b2tfnb = ncall_b2tfnb+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfnb

      end subroutine b2tfnb
