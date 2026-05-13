*-----------------------------------------------------------------------
*.specification

      subroutine b2tfcc_mol (nCv, nFc, nVx, ns, switch, geo, mpg, pl,
     &                   dv, co)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_fluid_molecules
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tfcc
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Coeff), intent (in) :: co
*   ..input/output arguments
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      real (kind=R8) ::
     *   nw(nCv), nwf(nFc),
     *   dnw(nFc,0:1), tmf(nFc), tm(nCv),
     *   wght(nFc,2), wrkf(nFc), wrkVx(nVx),
     *   flo(nFc,0:1), con(nFc,0:1), flo0(nFc,0:1), con0(nFc,0:1)
*   ..procedures
      external xertst
      external b2xvsg, diff, intface
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfcc_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

*   ..interpolate te and ti to cell faces
      wght = 1.0_R8
      tm = t_mol(:,0)
      call intface(nCv,nFc,mpg%fcCv,wght,tm,tmf)

* ..correct coefficients

      nw = n_mol(:,0)

*    ..pre-compute radial differences nw
      call diff(nCv, nFc, nVx, 0, geo, mpg, nw, wrkVx, dnw)

*    ..interpolate nw to cell faces
      call intface(nCv,nFc,mpg%fcCv,wght,nw,nwf)



*   ..correct floh, conh, for molecules

*    ..compute poloidal/radial coefficients for pressure diffusion term

      flo(:,0) = -switch%xfac*2.5_R8*co%cdp_mol(:,0)*dnw(:,0)*tmf 
      flo(:,1) = -switch%xfac*2.5_R8*co%cdp_mol(:,1)*dnw(:,1)*tmf
      con(:,0) =  switch%xfac*2.5_R8*co%cdp_mol(:,0)*nwf*tmf
      con(:,1) =  switch%xfac*2.5_R8*co%cdp_mol(:,1)*nwf*tmf


*    ..compute coefficients for stencil
      call calccoef(nCv,nFc,nVx,0,geo,flo,con,flo0,con0)

*    ..add to previous contributions
      dv%floh_mol        = dv%floh_mol        + flo0
      dv%conh_mol(:,:,0) = dv%conh_mol(:,:,0) + con0


*   ..add contributions from Rhie-Chow term to floh, conh, for molecules

      do iFc = 1, nFc
        con(iFc,0) = switch%b2npht_pcm1*1.25_R8*geo%fcPbs(iFc)*
     &    (dv%pcc_mol(mpg%fcCv(iFc,1),1)*
     &     nw(mpg%fcCv(iFc,1))*tm(mpg%fcCv(iFc,1))  +
     &     dv%pcc_mol(mpg%fcCv(iFc,2),1)*
     &     nw(mpg%fcCv(iFc,2))*tm(mpg%fcCv(iFc,2)) )/
     &    (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
      enddo

*    ..compute poloidal/radial coefficients for pressure diffusion term
       flo(:,0) = -con(:,0)*dnw(:,0)
       con(:,0) =  con(:,0)*nwf

       flo(:,1) = 0.0_R8
       con(:,1) = 0.0_R8

*    ..compute coefficients for stencil
       call calccoef(nCv,nFc,nVx,0,geo,flo,con,flo0,con0)

*    ..add to previous contributions
       dv%floh_mol        = dv%floh_mol        + flo0
       dv%conh_mol(:,:,0) = dv%conh_mol(:,:,0) + con0

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfcc_mol

      end subroutine b2tfcc_mol
