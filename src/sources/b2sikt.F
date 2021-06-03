*-----------------------------------------------------------------------
*.specification

      subroutine b2sikt (nCv, nFc, nVx, ns, ismain, switch, geo, mpg, 
     &   pl, dv, rza, co, st_ext, she0, shi0, skt0, skt_prod, skt_diss)
      use b2mod_types
!      use b2mod_indirect
!      use b2mod_geo_corner
!      use b2mod_geo
!     & , only : crx
!      use b2mod_sources
!     & , only : skt_prod, skt_diss
!      use b2mod_rates
!     & , only : rza
!      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
c      use b2mod_keps
      use b2mod_user_namelist
     , , only : nomp, omp, icsep
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sikt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Coeff), intent (in) :: co
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) :: rza(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   she0(nCv,0:3), shi0(nCv,0:3), 
     *   skt0(nCv,0:3), skt_prod(nCv), skt_diss(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIKT computes the sources in the energy equations due to the
*     k-eps model.
*
*-----------------------------------------------------------------------
*.declarations
      
*   ..local variables
      integer ic

      real (kind=R8) :: wrks(nCv), wrkv(nVx), rz(nCv), 
     &  pz(nCv), cs(nCv), t0, t1, t2, t3,
     &  wrk0(nCv), wrk1(nCv), wrk2(nCv), wrk3(nCv), 
     &  wrkg0(nCv,0:1), wrkg1(nCv,0:1),
     &  dnac(nCv,0:1), dlbc(nCv,0:1), dtec(nCv,0:1), dtic(nCv,0:1),
     &  fac_diss, fac_sheath, csfs(nCv), connfs(nCv), rhol(nCv),
     &  lpar(nCv)
      integer, save :: min_source = 0
      character chk*1
*   ..procedures
      intrinsic min, max, sqrt
      external subini, subend, xertst, ipgeti, ipgetr, sfill
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sikt')
*   ..set internal parameters on first call
      if (ncall_b2sikt.eq.0) then
       call ipgeti('b2sikt_min_source',min_source)
       if ((switch%b2sikt_local.ne.1).or.
     &     (switch%b2sikt_sheath_local.ne.1)) then
        call xertst(nomp.gt.0,
     &    'No CVs in omp list, check rzomp in b2.user.parameters')
        call xertst(icsep.gt.0,
     &    'Invalid icsep value, check rzomp in b2.user.parameters')
       endif
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2sikt.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of na, te, ti, ne, ni
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%kt, 1, 'kt', '.ge.')
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv*2, dv%ni, 1, 'ni', '.gt.')
      endif

* ..compute heat source terms
*   ..initialize sources to zero
      call sfill (nCv*4, 0.0e0_R8, she0, 1)
      call sfill (nCv*4, 0.0e0_R8, shi0, 1)
      call sfill (nCv*4, 0.0e0_R8, skt0, 1)
      call sfill (nCv, 0.0e0_R8, skt_prod, 1)
      call sfill (nCv, 0.0e0_R8, skt_diss, 1)
      call sfill (nCv*2, 0.0e0_R8, dnac, 1)
      call sfill (nCv*2, 0.0e0_R8, dlbc, 1)
      call sfill (nCv*2, 0.0e0_R8, dtec, 1)
      call sfill (nCv*2, 0.0e0_R8, dtic, 1)
      if (switch%b2sikt_model.eq.1) then
cwdk  KU Leuven k model
cwdk  source : driven by ExB energy flux
cwdk  sink   : 
cwdk   - ad hoc term, ~1/(connection length), as a proxy for enhanced
cwdk     parallel transport due to parallel current fluctuations
        call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
        call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
        cs   = sqrt(pz/rz)
        if (switch%b2sikt_sheath_local .eq. 1) then
          csfs = cs
          connfs = geo%ftConn
          rhol = sqrt(2.0_R8*pl%ti*(am(ismain)*mp))/
     &       abs(rza(:,ismain)*qe*geo%cvBb(:,3))
        else
          csfs = cs(omp(icsep))
          connfs = geo%ftConn(omp(icsep))
          rhol = sqrt(2.0_R8*pl%ti(omp(icsep))*(am(ismain)*mp))/
     &       abs(rza(omp(icsep),ismain)*qe*geo%cvBb(omp(icsep),3))
        endif
        call xertst(minval(connfs(1:mpg%nCi)).gt.0.0_R8,
     .    ' connfs in b2sikt <=0!!')
        lpar = connfs
        if (switch%b2sikt_fac_diss_core_mode.eq.1) then
          where (mod(mpg%cvReg,4) .eq. 1)
            lpar = rhol*geo%cvBb(:,2)/geo%cvBb(:,0)
          end where
        endif
*   ..compute wrk arrays
        if (switch%b2sikt_style .eq. 0) then
*      ..logarithm of magnetic field in cell center
         wrks(:) = 2.0_R8*log(geo%cvBb(:,3))
*      ..gradients na and log(bb) in poloidal and radial direction
         call gradc(nCv, nFc, nVx, 0, geo, mpg, 
     &      pl%na(:,ismain), wrkv, dnac)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, wrks, wrkv, dlbc)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dtec)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, dtic)
*      ..computations for source due to Reynolds stress
         if (switch%b2sikt_fac_vis_RS.ne.0.0_R8) then
          call intcell(nFc,nCv,mpg,mpg%intcellP,dv%vaecrb(:,0,ismain),
     &      wrk0)
          call intcell(nFc,nCv,mpg,mpg%intcellR,dv%vaecrb(:,1,ismain),
     &      wrk1)
          wrk2 = wrk0*dlbc(:,0) + wrk1*dlbc(:,1)
          wrks = -geo%cvVol*pl%na(:,ismain)*pl%kt*wrk2
          call gradc_div(nCv, nFc, nVx, 1, geo, mpg, wrk0, 
     &      dv%vaecrb(:,0,ismain), wrkg0)
          call gradc_div(nCv, nFc, nVx, 1, geo, mpg, wrk1, 
     &      dv%vaecrb(:,1,ismain), wrkg1)
          wrkg0(:,0) = wrkg0(:,0)*geo%cvBb(:,2)/geo%cvBb(:,3)
          wrkg1(:,0) = wrkg1(:,0)*geo%cvBb(:,2)/geo%cvBb(:,3)
          wrks = wrks - switch%b2sikt_fac_vis_RS*
     &     geo%cvVol*co%dna_ExB*am(ismain)*mp*pl%na(:,ismain)*
     &     (wrkg0(:,0)**2 + wrkg1(:,1)**2 +
     &      0.5_R8*(wrkg1(:,0)**2 + wrkg0(:,1)**2) +
     &      wrkg1(:,0)*wrkg0(:,1))
         else
           wrks = 0.0_R8
         endif
*      ..individual source contributions
         wrk0 = geo%cvVol*(pl%ti*co%dna_ExB*
     &      (dnac(:,0)*dlbc(:,0) + dnac(:,1)*dlbc(:,1)) +
     &       switch%b2sikt_fac_aniso*2.0_R8/3.0_R8*co%hci_ExB*
     &      (dtic(:,0)*dlbc(:,0) + dtic(:,1)*dlbc(:,1)))
         wrk1 = geo%cvVol*(pl%te*co%dna_ExB*
     &      (dnac(:,0)*dlbc(:,0) + dnac(:,1)*dlbc(:,1)) +
     &       switch%b2sikt_fac_aniso*2.0_R8/3.0_R8*co%hce_ExB*
     &      (dtec(:,0)*dlbc(:,0) + dtec(:,1)*dlbc(:,1)))
         wrk2 = dv%ni(:,1)*pl%kt*geo%cvVol*csfs/connfs
         wrk3 = co%sigx_kt*pl%kt*geo%cvVol/lpar**2
        else
         call xerrab('wrong value of b2sikt_style')
        endif

*     ..sources
        do ic = 1, mpg%nCi
         if (mod(mpg%cvReg(ic),4) .eq. 1) then
          fac_sheath = switch%b2sikt_fac_sheath_core
          fac_diss = switch%b2sikt_fac_diss_core
         else
          fac_sheath = switch%b2sikt_fac_sheath
          fac_diss = switch%b2sikt_fac_diss
         endif
         t0 = wrk0(ic)
         t1 = wrk1(ic)
         t2 = fac_sheath*wrk2(ic) + fac_diss*wrk3(ic)
         t3 = wrks(ic)
         if (min_source.eq.1 .and. (t0+t1-t2).lt.0.0_R8) then
           t0 = 0.0_R8
           t1 = 0.0_R8
           t2 = 0.0_R8
         endif
         she0(ic,0)   = -t1 + t2
         shi0(ic,0)   = -t0
         skt0(ic,0)   =  t0 + t1 - t2 + t3
         skt_prod(ic) =  t0 + t1
         skt_diss(ic) =  t2 - t3
        enddo
      else
        call xerrab('wrong value of b2sikt_model')
      endif

* ..output
      if (switch%b2sikt_iout.eq.1) then
        wrks=she0(:,0)+she0(:,1)*pl%te+
     &       she0(:,2)*dv%ne + she0(:,3)*pl%te*dv%ne
        call my_out_us(70,nCv,0,wrks,'b2sikt_she')
        wrks=shi0(:,0)+shi0(:,1)*pl%ti+
     &       shi0(:,2)*dv%ni(:,0)+shi0(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrks,'b2sikt_shi')
        wrks=skt0(:,0)+skt0(:,1)*pl%kt+
     &       skt0(:,2)*dv%ni(:,1)+skt0(:,3)*pl%kt*dv%ni(:,1)
        call my_out_us(70,nCv,0,wrks,'b2sikt_skt')
        call my_out_us(70,nCv,0,skt_prod,'b2sikt_skt_prod')
        call my_out_us(70,nCv,0,skt_diss,'b2sikt_skt_diss')
        wrks = geo%cvVol(:)*(
     &     (pl%ti(:)+pl%te(:))*co%dna_ExB(:)*
     &     (dnac(:,0)*dlbc(:,0)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_na_x')
        wrks = geo%cvVol(:)*(
     &     (pl%ti(:)+pl%te(:))*co%dna_ExB(:)*
     &     (dnac(:,1)*dlbc(:,1)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_na_y')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hce_ExB(:)*switch%b2sikt_fac_aniso*
     &     (dtec(:,0)*dlbc(:,0)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_te_x')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hce_ExB(:)*switch%b2sikt_fac_aniso*
     &     (dtec(:,1)*dlbc(:,1)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_te_y')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hci_ExB(:)*switch%b2sikt_fac_aniso*
     &     (dtic(:,0)*dlbc(:,0)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_ti_x')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hci_ExB(:)*switch%b2sikt_fac_aniso*
     &     (dtic(:,1)*dlbc(:,1)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_ti_y')
        call my_out_us(70,nCv,0,connfs,'b2sikt_connfs')
        call my_out_us(70,nCv,0,csfs,'b2sikt_csfs')
        call my_out_us(70,nCv,0,dnac(:,0),'b2sikt_dnac0')
        call my_out_us(70,nCv,0,dnac(:,1),'b2sikt_dnac1')
        call my_out_us(70,nCv,0,dtec(:,0),'b2sikt_dtec0')
        call my_out_us(70,nCv,0,dtec(:,1),'b2sikt_dtec1')
        call my_out_us(70,nCv,0,dtic(:,0),'b2sikt_dtic0')
        call my_out_us(70,nCv,0,dtic(:,1),'b2sikt_dtic1')
        call my_out_us(70,nCv,0,dlbc(:,0),'b2sikt_dlbc0')
        call my_out_us(70,nCv,0,dlbc(:,1),'b2sikt_dlbc1')
!!some temporary/experimental output computations for source due to Reynolds stress
        call intcell(nFc,nCv,mpg,mpg%intcellP,dv%vaecrb(:,0,ismain),
     &    wrk0)
        call intcell(nFc,nCv,mpg,mpg%intcellR,dv%vaecrb(:,1,ismain),
     &    wrk1)
        wrk2 = wrk0*dlbc(:,0) + wrk1*dlbc(:,1)
        wrks = -geo%cvVol*pl%na(:,ismain)*pl%kt*wrk2
        call my_out_us(70,nCv,0,wrks,'b2sikt_RS_symmetric')

        call gradc_div(nCv, nFc, nVx, 1, geo, mpg, wrk0, 
     &    dv%vaecrb(:,0,ismain), wrkg0)
        call gradc_div(nCv, nFc, nVx, 1, geo, mpg, wrk1, 
     &    dv%vaecrb(:,1,ismain), wrkg1)
        wrkg0(:,0) = wrkg0(:,0)*geo%cvBb(:,2)/geo%cvBb(:,3)
        wrkg1(:,0) = wrkg1(:,0)*geo%cvBb(:,2)/geo%cvBb(:,3)
        wrks = geo%cvVol*co%dna_ExB*am(ismain)*mp*pl%na(:,ismain)*
     &   (wrkg0(:,0)**2 + wrkg1(:,1)**2 +
     &    0.5_R8*(wrkg1(:,0)**2 + wrkg0(:,1)**2) +
     &    wrkg1(:,0)*wrkg0(:,1))
        call my_out_us(70,nCv,0,wrks,'b2sikt_RS_deviatoric')
      endif

* ..return
      ncall_b2sikt = ncall_b2sikt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sikt

      end


