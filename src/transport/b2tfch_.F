*-----------------------------------------------------------------------
*.specification

      subroutine b2tfch_ (nCv, nFc, nVx, ns, ismain, ismain0, switch, !srv 29.01.20 {
     &  geo, mpg, pl, dv, co, rt, st_ext)
      use b2mod_types
!      use b2mod_boundary_namelist
      use b2mod_constants
!      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ad
     & , only : ncall_b2tfch
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : my_out_folder, ncall_b2xehx, ncall_b2tcpa,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2ttia
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, ismain0
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl
      type(B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
*   ..input/output arguments (unspecified on entry)
      type(B2Coeff), intent (inout) :: co
      type(B2Derivatives), intent (inout) :: dv
*   ..local variables
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFCH_ computes the currents and the coefficients that specify the
*     linearized response to changes in electric potential.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, is0
      character(3) :: charns, charns0
      real (kind=R8) :: csig_stoch(nFc), fchstoch(nFc)                   !lk 31.01.08
      real (kind=R8) ::
     &    ehx(nFc)                                    !srv 17.05.99
      real (kind=R8) :: facdriftm, fac_ExBm                              !srv 23.12.19
      real (kind=R8) ::
     *    wrk(nFc,0:1), dumm0(nFc,0:1), dumm1(nFc,0:1), zeff(nCv)
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c       real (kind=R8), save :: drift_hyb = 1.0_R8, fch_inert = 1.0_R8
c      real (kind=R8), save :: fch_stochastic = 0.0_R8                   !lk 31.01.08
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2xvfx, b2ttia                            !srv 11.07.99
*   ..initialization

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfch_')
*   ..set internal parameters on first call
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c       call ipgetr ('b2tfhe_drift_hyb',drift_hyb)
c       call ipgetr ('b2tiner_inert', fch_inert)                          !srv 07.07.99
c       call ipgetr ('b2tfhe_stochastic', fch_stochastic)                 !lk 31.01.08
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2tfch.lt.3) then
*    ..test sign of ne, te, csig, chce
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       wrk(:,0) = co%csig(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc, wrk(:,0), 1, 'csig0', '.ge.')
       wrk(:,1) = co%csig(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc, wrk(:,1), 1, 'csig1', '.ge.')
       wrk(:,0) = co%chce(:,0)*geo%fcQalf(:,0)
       call b2xvsg (nFc, wrk(:,0), 1, 'chce0', '.ge.')
       wrk(:,1) = co%chce(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc, wrk(:,1), 1, 'chce1', '.ge.')
      endif

      facdriftm = maxval(dv%facdrift)
      fac_ExBm  = maxval(dv%fac_ExB)                                     !srv 23.12.19
      fchstoch     = 0.0_R8
      csig_stoch   = 0.0_R8
!WG_TODO      avisper = 0.0_R8
!WG_TODO      bvisper = 0.0_R8
!WG_TODO      cvisper = 0.0_R8
      if (switch%no_current.ne.0) then
        dv%fch       = 0.0_R8
        dv%fchin     = 0.0_R8
        dv%fchdia    = 0.0_R8
        dv%fchvispar = 0.0_R8
        dv%fchvisper = 0.0_R8
        dv%fchvisq   = 0.0_R8
        dv%fchinert  = 0.0_R8
        dv%fch_p     = 0.0_R8
        dv%fchviskt  = 0.0_R8
        dv%fchanml   = 0.0_R8
        dv%conc      = 0.0_R8
      else
*     ..calculate the contribution of ion-neutral current                !srv 11.02.04 {
        call b2tinnt (nCv, nFc, nVx, ns, switch, geo, mpg, rt,
     &    dv%facdrift, dv%fac_ExB, ismain, ismain0, pl,
     &    dv%fna, co%csigin, dv%fchin)                     !srv 11.05.09
*     ..compute x,y-components of divergent part of diamagnetic current
        call b2tdia(nCv, nFc, nVx, ns, switch, geo, mpg, dv%pz, dv%fna,
     &    dv%uadia, dv%facdrift, dv%fchdia)
*     ..calculate x,y-components of viscosity parallel current
        call b2tvspa (nCv, nFc, nVx, ns, switch, geo, mpg, pl%ua,
     &    co%vsaf_cl, dv%fac_vis, dv%fchvispar)
*     ..calculate the contribution of viscosity perpendicular current    !srv 22.11.99
        call b2tvspr(nCv, nFc, nVx, ns, switch, geo, mpg,
     &    pl%po, pl%ti, pl%na, rt%rza, co%vsa0, dv%fac_vis,dv%fchvisper)
*     ..calculate the parallel viscosity current connected with the heat flux
        call b2tvsq (nCv, nFc, nVx, ns, switch, geo, mpg, dv%fac_vis,
     &    dv%ne2, rt%rz2, dv%lnlam, pl%te, pl%ti, dv%ni(:,1), co%hcix_c,
     &    dv%fchvisq)
*     ..calculate inertia current                                        !srv 07.07.99
        call b2tiner (nCv, nFc, nVx, ns, switch, geo, mpg, pl%na, pl%ua,
     &    dv%facdrift, dv%fchinert)
*     ..calculate the x,y-components of anomalous current
        call b2tanml (nCv, nFc, ns, switch, geo, mpg, co%csig_an,
     &    pl%po,dv%fchanml)
#ifdef WG_TODO
*     ..calculate the y-component of electron current associated with
*       stochastization of magnetic field
        call b2tstch (nx, ny, po, te, ne, csig_stoch, fchstoch)          !srv 25.01.02 {
#endif
*     ..calculate the parallel current
        call b2tcpa (nCv, nFc, nVx, ns, switch, geo, mpg, pl%te, pl%po,
     *  dv%ne, pl%na, pl%ua, rt%rza, rt%rz2, zeff, co%csig, co%calf,
     &  ehx, st_ext, dv%fch_p, dv%fch_pi_c, dv%fch_pi_f)                 !srv 13.01.17
*     ..compute current due to RS of kt-model
        call b2tvskt(nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &    pl%po, pl%na, pl%kt, dv%vaecrb, co%dna_ExB, dv%fchviskt)

*     ..compute total currents
        dv%fch = dv%fch_p + dv%fchin + dv%fchdia + dv%fchvispar +
     &    dv%fchvisq + dv%fchinert + dv%fchanml + dv%fchviskt +
     &    dv%fchvisper
        dv%fch(:,1) = dv%fch(:,1) + fchstoch

*     ..compute conc
        dumm0 = 0.0_R8
        wrk   = co%csig
        if(switch%fch_ion_neutral.ne.0.0_R8 .and.
     &    (facdriftm.gt.0.0_R8 .or. fac_ExBm.gt.0.0_R8)) then            !srv 23.12.19
         do is = 0, ns-1
          if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
            do is0 = 0, ns-1
              if (is_neutral(is0).and.nint(zn(is0)).eq.1) then
                wrk(:,0) = wrk(:,0) + co%csigin(:,0,is,is0)*
     &            (geo%fcBb(:,2)/geo%fcBb(:,3))**2
                wrk(:,1) = wrk(:,1) + co%csigin(:,1,is,is0)
              endif
            enddo
          endif
         enddo
        endif
#ifdef WG_TODO
        if(fch_stochastic.ne.0.0_R8 ) then                               !lk 31.01.08 {
          wrk(:,1) = wrk(:,1) + csig_stoch
        endif                                                            !lk 31.01.08 }
#endif
        call calccoef(nCv,nFc,nVx,0,geo,dumm0,wrk,dumm1,dv%conc)

      endif !no_current

      if(switch%b2tfch_iout.ne.0) then                                   !srv 18.06.08 {
       call my_out_us(70,nFc,1,dv%fch(1,0),'b2tfch__fch_th')
       call my_out_us(70,nFc,1,dv%fch(1,1),'b2tfch__fch_r')
       call my_out_us(70,nFc,1,dv%fchin(1,0),'b2tfch__fchin_th')
       call my_out_us(70,nFc,1,dv%fchin(1,1),'b2tfch__fchin_r')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfch__fchvispar_th')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfch__fchvispar_r')
       call my_out_us(70,nFc,1,dv%fchvisq(1,0),'b2tfch__fchvisq_th')
       call my_out_us(70,nFc,1,dv%fchvisq(1,1),'b2tfch__fchvisq_r')
       call my_out_us(70,nFc,1,dv%fchvisper(1,0),'b2tfch__fchvisper_th')
       call my_out_us(70,nFc,1,dv%fchvisper(1,1),'b2tfch__fchvisper_r')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfch__fchinert_th')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfch__fchinert_r')
       call my_out_us(70,nFc,1,dv%fchdia(1,0),'b2tfch__fchdia_th')
       call my_out_us(70,nFc,1,dv%fchdia(1,1),'b2tfch__fchdia_r')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfch__fchanml_th')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfch__fchanml_r')
       call my_out_us(70,nFc,1,fchstoch,'b2tfch__fchstoch_r')             !lk 31.01.08
       call my_out_us(70,nFc,1,dv%conc(1,0),'b2tfch__conc_th')
       call my_out_us(70,nFc,1,dv%conc(1,1),'b2tfch__conc_r')
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfch__fch_p_th')
       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfch__fch_p_r')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2tfch__csig_an_th')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2tfch__csig_an_r')
       call my_out_us(70,nFc,1,co%csig(1,0),'b2tfch__csig_th')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2tfch__csig_r')
       call my_out_us(70,nFc,1,co%calf(1,0),'b2tfch__calf_th')
       call my_out_us(70,nFc,1,co%calf(1,1),'b2tfch__calf_r')
       call my_out_us(70,nCv,0,pl%te,'b2tfch__te')
       call my_out_us(70,nCv,0,dv%ne,'b2tfch__ne')
       call my_out_us(70,nCv,0,pl%po,'b2tfch__po')
       do is = 0, ns-1
        write (charns,'(i3.3)') is
        call my_out_us(70,nCv,0,pl%na(1,is),'b2tfch__na'//charns)
        do is0=0,ns-1
         write (charns0,'(i3.3)') is0
         call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),
     &        'b2tfch__csigin_th'//charns//'_'//charns0)
         call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),
     &        'b2tfch__csigin_r'//charns//'_'//charns0)
        enddo
       enddo
      endif                                                              !srv 18.06.08 }

      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0
     &                           .or. switch%iout_b2wdat.eq.4) then
       call my_out_us(70,nFc,1,dv%fch(1,0),'b2tfch__fch_th')
       call my_out_us(70,nFc,1,dv%fch(1,1),'b2tfch__fch_r')
       call my_out_us(70,nFc,1,dv%fchin(1,0),'b2tfch__fchin_th')
       call my_out_us(70,nFc,1,dv%fchin(1,1),'b2tfch__fchin_r')
       call my_out_us(70,nFc,1,dv%fchvispar(1,0),'b2tfch__fchvispar_th')
       call my_out_us(70,nFc,1,dv%fchvispar(1,1),'b2tfch__fchvispar_r')
       call my_out_us(70,nFc,1,dv%fchvisq(1,0),'b2tfch__fchvisq_th')
       call my_out_us(70,nFc,1,dv%fchvisq(1,1),'b2tfch__fchvisq_r')
       call my_out_us(70,nFc,1,dv%fchvisper(1,0),'b2tfch__fchvisper_th')
       call my_out_us(70,nFc,1,dv%fchvisper(1,1),'b2tfch__fchvisper_r')
       call my_out_us(70,nFc,1,dv%fchinert(1,0),'b2tfch__fchinert_th')
       call my_out_us(70,nFc,1,dv%fchinert(1,1),'b2tfch__fchinert_r')
       call my_out_us(70,nFc,1,dv%fchdia(1,0),'b2tfch__fchdia_th')
       call my_out_us(70,nFc,1,dv%fchdia(1,1),'b2tfch__fchdia_r')
       call my_out_us(70,nFc,1,dv%fchanml(1,0),'b2tfch__fchanml_th')
       call my_out_us(70,nFc,1,dv%fchanml(1,1),'b2tfch__fchanml_r')
       call my_out_us(70,nFc,1,fchstoch,'b2tfch__fchstoch_r')             !lk 31.01.08
       call my_out_us(70,nFc,1,dv%fch_p(1,0),'b2tfch__fch_p_th')
!       call my_out_us(70,nFc,1,dv%fch_p(1,1),'b2tfch__fch_p_r')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2tfch__csig_an_th')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2tfch__csig_an_r')
      endif                                                              !srv 11.09.09 }
* ..return
      ncall_b2tfch = ncall_b2tfch+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfch_

      end subroutine b2tfch_                                             !srv 09.01.20 }
