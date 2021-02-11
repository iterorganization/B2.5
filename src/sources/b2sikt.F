*-----------------------------------------------------------------------
*.specification

      subroutine b2sikt (nCv, nFc, nVx, ns, ismain, geo, mpg, pl, 
     &  dv, rza, co, st_ext, she0, shi0, skt0, skt_prod, skt_diss)
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
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sikt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
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
      integer, save :: iout = 0
      integer, save :: style = 0
      integer, save :: model = 1, keps_local = 1, fac_sheath_local=1
      real (kind=R8) :: wrks(nCv), wrkv(nVx), rz(nCv), 
     &  pz(nCv), cs(nCv), t0, t1, t2, wrk0(nCv), wrk1(nCv), wrk2(nCv),
     &  dnac(nCv,0:1), dlbc(nCv,0:1), dtec(nCv,0:1), dtic(nCv,0:1),
     &  fac_diss, fac_sheath, csfs(nCv), connfs(nCv)
      real (kind=R8), save ::
     &  fac_sheath_rest=1.0_r8, fac_sheath_core=0.0_r8,
     &  fac_aniso=1.0_R8
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
!       call get_jsep(nx,ny,jxi,jxa,jsep)
       call ipgeti ('b2sikt_iout', iout)
       call ipgeti ('b2sikt_style', style)
       call ipgeti ('b2sikt_model', model)
       call ipgeti ('b2sikt_keps_local', keps_local)
       call ipgetr ('b2sikt_aniso', fac_aniso)
       call ipgetr ('b2sikt_fac_sheath', fac_sheath_rest)
       call ipgetr ('b2sikt_fac_sheath_core', fac_sheath_core)
       call ipgeti ('b2sikt_fac_sheath_local', fac_sheath_local)
       call ipgetr ('b2sikt_anisothermal', fac_aniso)
       if ((keps_local.ne.1).or.(fac_sheath_local.ne.1)) then
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
      if (model.eq.1) then
cwdk  KU Leuven k model
cwdk  source : driven by ExB energy flux
cwdk  sink   : 
cwdk   - ad hoc term, ~1/(connection length), as a proxy for enhanced
cwdk     parallel transport due to parallel current fluctuations
        call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
        call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
        cs  = sqrt(pz/rz)
        if (fac_sheath_local.eq.1) then
          csfs = cs
          connfs = geo%ftConn
        else
          csfs = cs(omp(icsep))
          connfs = geo%ftConn(omp(icsep))
        endif
        call xertst(minval(connfs(1:mpg%nCi)).gt.0.0_R8,
     .    ' connfs in b2sikt <=0!!')
*   ..compute wrk arrays
        if (style.eq.0) then
*      ..logarithm of magnetic field in cell center
         wrks(:) = 2.0_R8*log(geo%cvBb(:,3))
*      ..gradients na and log(bb) in poloidal and radial direction
         call gradc(nCv, nFc, nVx, 0, geo, mpg, 
     &      pl%na(:,ismain), wrkv, dnac)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, wrks, wrkv, dlbc)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dtec)
         call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, dtic)
*      ..individual source contributions
         wrk0 = geo%cvVol*(pl%ti*co%dna_ExB*
     &      (dnac(:,0)*dlbc(:,0) + dnac(:,1)*dlbc(:,1)) +
     &       fac_aniso*2.0_R8/3.0_R8*co%hci_ExB*
     &      (dtic(:,0)*dlbc(:,0) + dtic(:,1)*dlbc(:,1)))
         wrk1 = geo%cvVol*(pl%te*co%dna_ExB*
     &      (dnac(:,0)*dlbc(:,0) + dnac(:,1)*dlbc(:,1)) +
     &       fac_aniso*2.0_R8/3.0_R8*co%hce_ExB*
     &      (dtec(:,0)*dlbc(:,0) + dtec(:,1)*dlbc(:,1)))
         wrk2 = dv%ni(:,1)*pl%kt*geo%cvVol*csfs/connfs
        else
         call xerrab('wrong value of b2sikt_style')
        endif

*     ..sources
        do ic = 1, mpg%nCi
         if (mod(mpg%cvReg(ic),4) .eq. 1) then
          fac_sheath = fac_sheath_core
         else
          fac_sheath = fac_sheath_rest
         endif
         t0 = wrk0(ic)
         t1 = wrk1(ic)
         t2 = fac_sheath*wrk2(ic)
         she0(ic,0)   = -t1 + t2
         shi0(ic,0)   = -t0
         skt0(ic,0)   =  t0 + t1 - t2
         skt_prod(ic) =  t0 + t1
         skt_diss(ic) =  t2
        enddo
      else
        call xerrab('wrong value of b2sikt_model')
      endif

* ..output
      if (iout.eq.1) then
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
     &      2.0_R8/3.0_R8*co%hce_ExB(:)*fac_aniso*
     &     (dtec(:,0)*dlbc(:,0)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_te_x')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hce_ExB(:)*fac_aniso*
     &     (dtec(:,1)*dlbc(:,1)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_te_y')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hci_ExB(:)*fac_aniso*
     &     (dtic(:,0)*dlbc(:,0)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_ti_x')
        wrks = geo%cvVol(:)*(
     &      2.0_R8/3.0_R8*co%hci_ExB(:)*fac_aniso*
     &     (dtic(:,1)*dlbc(:,1)))
        call my_out_us(70,nCv,0,wrks,'b2sikt_prod_ti_y')
      endif

* ..return
      ncall_b2sikt = ncall_b2sikt+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sikt

      end


