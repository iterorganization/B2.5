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
      integer ncall, ic
      integer, save :: iout = 0
      integer, save :: solve_keps = 0
      integer, save :: model = 1, keps_local = 1, fac_sheath_local=1
      real (kind=R8) :: wrks(nCv), wrkv(nVx), ro(nCv), 
     &  pr(nCv), cs(nCv), t0, t1,
     &  dnac(nCv,0:1), dlbc(nCv,0:1), dtec(nCv,0:1), dtic(nCv,0:1),
     &  rhol, css, fac_diss, delta, fac_sheath, csfs, connfs
      real (kind=R8), save :: delta_om = 0.0e-1_R8,
     &  fac_theta = 1.0e0_R8, fac_gamma_diss = 0.00e0_R8, 
     &  delta_om_core = 0.0-1_R8, fac_gamma_diss_core = 0.00e0_R8,
     &  fac_sheath_rest=1.0_r8, fac_sheath_core=0.0_r8,
     &  fac_aniso=0.0_R8
      character chk*1
      save ncall
*   ..procedures
      intrinsic min, max, sqrt
      external subini, subend, xertst, ipgeti, ipgetr, sfill
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sikt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
!       call get_jsep(nx,ny,jxi,jxa,jsep)
       call ipgeti ('b2sikt_iout', iout)
       call ipgeti ('b2sikt_model', model)
       call ipgetr ('b2sikt_fac_gamma_diss', fac_gamma_diss)
       call ipgetr ('b2sikt_fac_gamma_diss_core', fac_gamma_diss_core)
       call ipgetr ('b2sikt_fac_theta', fac_theta)
       call ipgetr ('b2sikt_delta_om', delta_om)
       call ipgetr ('b2sikt_delta_om_core', delta_om_core)
       call ipgeti ('b2sikt_keps_local', keps_local)
       if (keps_local.ne.1) then
        call xertst(nomp.gt.0,
     &    'No CVs in omp list, check rzomp in b2.user.parameters')
        call xertst(icsep.gt.0,
     &    'Invalid icsep value, check rzomp in b2.user.parameters')
       endif
       call ipgeti ('b2mndr_solve_keps', solve_keps)
       call ipgetr ('b2sikt_fac_sheath', fac_sheath_rest)
       call ipgetr ('b2sikt_fac_sheath_core', fac_sheath_core)
       call ipgeti ('b2sikt_fac_sheath_local', fac_sheath_local)
       call ipgetr ('b2sikt_anisothermal', fac_aniso)
      endif
*   ..test nCv, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
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
      if (model.eq.1) then
cwdk  KU Leuven k model
cwdk  source : driven by ExB particle flux
cwdk  sink   : add hoc term, k^2 scaling, controled by delta_om
cwdk  model adaptations based on T2D simulations: 
cwdk           strong direct local dissipation of the source term,  
cwdk           with factor fac_gamma_diss (~97%)
        call b2xpro (nCv, ns, mp, am, pl%na, ro, st_ext)
        call b2xppr (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pl%tn, pr, 
     &             st_ext)
        cs  = sqrt(pr/ro)
*   ..compute wrk arrays
*     ..logarithm of magnetic field in cell center
        wrks(:) = log(geo%cvBb(:,3)**2)
*     ..gradients na and log(bb) in poloidal and radial direction
csc FIXME   if(leftix(ix,iy).eq.-2) then
c            wrk2(ix,iy) = 0.0_R8
c            wrk3(ix,iy) = 0.0_R8
c           else
        dnac = 0.0_R8
        dlbc = 0.0_R8
        dtec = 0.0_R8
        dtic = 0.0_R8
        call gradc(nCv, nFc, nVx, 0, geo, mpg, 
     &     pl%na(:,ismain), wrkv, dnac)
        call gradc(nCv, nFc, nVx, 0, geo, mpg, wrks, wrkv, dlbc)
        call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dtec)
        call gradc(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, dtic)

*     ..sources
        do ic = 1, mpg%nCi
csc FIXME  if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1 .and.
c     &       bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1) then
         if (mod(mpg%cvReg(ic),4) .eq. 1) then
          fac_diss = fac_gamma_diss_core
          delta    = delta_om_core
          fac_sheath = fac_sheath_core
         else
          fac_diss = fac_gamma_diss
          delta    = delta_om
          fac_sheath = fac_sheath_rest
         endif
         if (keps_local.eq.1) then
          rhol = am(ismain)*mp*sqrt(2.0_R8*pl%ti(ic)/(am(ismain)*mp))/
     &       abs(rza(ic,ismain)*qe*geo%cvBb(ic,3))
          css = cs(ic)
         else
          rhol = am(ismain)*mp*sqrt(2.0_R8*pl%ti(omp(icsep+1))/
     &       (am(ismain)*mp))/abs(rza(omp(icsep+1),ismain)*
     &       qe*geo%cvBb(omp(icsep+1),3))
          css = cs(omp(icsep+1))
         endif
         if (fac_sheath_local.eq.1) then
          csfs = cs(ic)
          connfs = geo%ftConn(ic)
         else
          csfs = cs(omp(icsep+1))
          connfs = geo%ftConn(omp(icsep+1))
         endif
! FIXME check here if it's OK, before it was (dnac(ic)+dnac(right(ic))*(...)*0.5
         t0 = geo%cvVol(ic)*(
     &     (pl%ti(ic)+pl%te(ic))*co%dna_ExB(ic)*
     &     (dnac(ic,0)*dlbc(ic,0) + dnac(ic,1)*dlbc(ic,1)) +
     &      co%hce_ExB(ic)*fac_aniso*
     &     (dtec(ic,0)*dlbc(ic,0) + dtec(ic,1)*dlbc(ic,1)) +
     &      co%hci_ExB(ic)*fac_aniso*
     &     (dtic(ic,0)*dlbc(ic,0) + dtic(ic,1)*dlbc(ic,1)))
         t1 = delta*geo%cvVol(ic)*dv%ni(ic,1)*pl%kt(ic)**2/
     &       (mp*am(ismain)*css*rhol) + fac_sheath*dv%ni(ic,1)*
     &       pl%kt(ic)*geo%cvVol(ic)*csfs/connfs
         she0(ic,0)   =  fac_diss*t0 + t1
         shi0(ic,0)   = -t0
         skt0(ic,0)   =  (1.0_R8 - fac_diss)*t0 - t1
         skt_prod(ic) =  (1.0_R8 - fac_diss)*t0
         skt_diss(ic) =  t1
        enddo
      else
        call xerrab('wrong value of b2sikt_model')
      endif

* ..output
      if (iout.eq.1) then
        wrks=she0(:,0)+she0(:,1)*pl%te+
     &       she0(:,2)*dv%ne + she0(:,3)*pl%te*dv%ne
        call my_out(70,nCv,0,wrks,'b2sikt_she')
        wrks=shi0(:,0)+shi0(:,1)*pl%ti+
     &       shi0(:,2)*dv%ni(:,0)+shi0(:,3)*pl%ti*dv%ni(:,0)
        call my_out(70,nCv,0,wrks,'b2sikt_shi')
        wrks=skt0(:,0)+skt0(:,1)*pl%kt+
     &       skt0(:,2)*dv%ni(:,1)+skt0(:,3)*pl%kt*dv%ni(:,1)
        call my_out(70,nCv,0,wrks,'b2sikt_skt')
        call my_out(70,nCv,0,skt_prod,'b2sikt_skt_prod')
        call my_out(70,nCv,0,skt_diss,'b2sikt_skt_diss')
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sikt

      end


