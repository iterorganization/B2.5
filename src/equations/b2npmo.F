*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo (nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &  itcnt, solving, rxf, pl, dv, rt, co, sr, st_ext, ierr)
      use b2mod_types
      use b2mod_constants
!WG_TODO      use b2mod_sources
!WG_TODO     1 , only : b2npmo_smaf, b2npmo_smag, b2npmo_smav,
!WG_TODO     2          smpr, smpt, smfr, smcf
!WG_TODO: smof -> part of b2mod_sources
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     1 , only : update_balance_sources_momentum, fmo_b2nxfv,
!WG_TODO     2          balance_netcdf
      use b2mod_numerics_namelist
      use b2mod_transport_fun                                          !srv 13.01.17
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ma28_for_us
     , , only : restart_ma28_for_us
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : b2npmo_solvedum, ncall_b2npmo, ncall_b2usmo,
     &          ncall_b2urmo, ncall_b2sicf, ncall_b2nxfc, ncall_b2sigp,
     &          ncall_b2nxdu, ncall_b2ttia, ncall_b2sifr,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          my_out_folder, ncall_b2tlnl, ncall_b2sifrtf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain, itcnt                          !srv 13.01.17
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Coeff), intent (in) :: co
      type (B2Rates), intent (in) :: rt
      type (B2Source), intent (in) :: sr
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  rxf
      logical solving
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
*   ..output arguments (unspecified on entry)
      integer ierr
!WG_TODO     &  smof(nCv,0:3,0:ns-1)
!WG_TODO_TMP these sources used to be in b2mod_sources
!            => add again to B2Source?
      real (kind=R8) ::
     *  smcf(nCv,0:ns-1), smpr(nCv,0:ns-1), smpt(nCv,0:ns-1),
     *  smfr(nCv,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPMO performs one iteration of a relaxation process for the
*     solution of the b2 equations for parallel velocity.
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equations that
*     govern the parallel velocity ua.
*     An outline of the procedure follows.
*     ..compute auxiliaries
*     ..do is=0,ns-1
*      ..compute residual for species (is)
*      ..compute correction for species (is)
*      ..compute contribution to total residual
*     ..enddo
*     ..compute simultaneous correction for all species
*     ..apply correction to all species
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, is, ireg, k
      integer :: vcount(0:ns-1)
      character chns*3, chk*1                                          !srv 04.07.08
      real (kind=R8) ::
     *  t0, t1, czetap, czetae                                !srv 13.01.17
      real (kind=R8) ::
     *  smbch(nCv,0:3), smbtf(nCv,0:3),                !srv 11.09.09
     *  smbgp(nCv,0:3),
     *  smbvh(nCv,0:3), smbvv(nCv,0:3),                !srv 23.07.99 04.10.99
     *  smbcf(nCv,0:3), smban(nCv,0:3),                !srv 09.07.01 25.09.17
     1  pb(nCv), cvct(nFc,0:1),
     2  rob(nCv), flub(nFc,0:1),
     *  flubv(nFc,0:1),                                    !srv 11.09.09
     3  wrkf(nFc,0:1), wrkff(nFc,0:1), wrkfff(nFc,0:1),
     *  wrkv(nVx), aa(mpg%nCmxNv),
     4  smb(nCv,0:3), flcb(nFc,0:1),
     5  cvcb(nFc,0:1), ro(nCv), rz(nCv),
     6  pr(nCv), pz(nCv), rob1(nCv), pb1(nCv),
     8  ehxp(nCv), smt(nCv,0:3), cst,                  !srv 17.01.14
     9  flct(nFc,0:1), wrk1(nCv), wrks(nCv),
     1  smqdu(nCv), wrk0(nCv), wrk2(nCv),
     2  ctcfb(nCv), ctcft(nCv),                        !srv 06.04.07
     3  zeff(nCv), zmff(nCv),                          !srv 13.01.17 {
     4  ce1(nCv), ce2(nCv),
     5  cimp1(nCv), cimp2(nCv),
     6  zetap(nCv), zetae(nCv),
     7  gti(nCv), gte(nCv), gtn(nCv),                   !srv 13.01.17 }
     8  wrkf0(nFc), wrkf1(nFc)
      real (kind=R8) ::
     *  smbfria(nCv,0:3), smbfrea(nCv,0:3),
     *  smbtfia(nCv,0:3), smbtfea(nCv,0:3)             !srv 13.01.17 05.07.17
      logical new_matrix, main_lightest
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *  med
      external subini, subend, xertst, ipgetr, sfill, b2scopy,
     &  b2sscal, b2saxpy
      external
     &  b2xvsg, b2xvff, b2xpne, b2xprz, b2xpro, b2xppr, b2xppz,
     &  b2xppb, b2nxdu, b2nxdv, b2nxfx, b2nxfc, b2nxst,
     &  b2sigp, b2siav, b2sicf, b2nxcm, b2urmo, b2usmo, med
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo')
*   ..set internal parameters on first call
      if (ncall_b2npmo.eq.0) then
        write(*,*) 'nsmin, nsmax = ', switch%nsmin, switch%nsmax
        main_lightest = .true.
        do is = 0, ns-1
          main_lightest = main_lightest .and.
     .     nint(am(is)).ge.nint(am(ismain))
        end do
        if (.not.main_lightest .and. switch%b2sigp_style.eq.2) then
         write(0,'(2a)') 'Warning: b2sigp_style.eq.2 not appropriate ',
     .    'for plasmas where the main species is not the lightest one!'
        endif
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8, 
     . 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall_b2npmo.lt.3) then
*    ..test sign of na, te, ti
       do is = 0, ns-1
        call b2xvsg (nCv, pl%na(:,is), 1, 'na', '.gt.')
       enddo
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%tn, 1, 'tn', '.gt.')
*    ..test sign of cvsa
       do is = 0, ns-1
        call b2xvsg (nFc,co%cvsa(:,0,is)*geo%fcQalf(:,0),1,
     .    'cvsa0','.ge.')
        call b2xvsg (nFc,co%cvsa(:,1,is)*geo%fcQalf(:,1),1,
     .    'cvsa1','.ge.')
       enddo
*    ..test sign of smo, smq
       do is = 0, ns-1
        call b2xvsg (nCv, sr%smo(:,1,is), 1, 'smo 1', '.le.')
        call b2xvsg (nCv, sr%smo(:,3,is), 1, 'smo 3', '.le.')
        call b2xvsg (nCv, sr%smq(:,1,is), 1, 'smq 1', '.le.')
        call b2xvsg (nCv, sr%smq(:,3,is), 1, 'smq 3', '.le.')
       enddo
      endif
*   ..compute ne, ne2, ro, rz, pr, pz
      call b2xpni (nCv, ns, pl%na, dv%ni)
      call b2xpnn (nCv, ns, pl%na, dv%nn)
      call b2xpne (nCv, ns, rt%rza, pl%na, st_ext%ne, dv%ne)
      call b2xpne (nCv, ns, rt%rz2, pl%na, st_ext%ne2, dv%ne2)
      call b2xpro (nCv, ns, mp, am, pl%na, ro, st_ext)
      call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
      call b2xppr (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pl%tn, pr, 
     &             st_ext)
      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
      call b2xzef (nCv, ns, rt%rz2, pl%na, dv%ne, zeff, st_ext)                      !srv 13.01.17 {
      call b2xzmf (nCv, ns, ismain, rt%rz2, pl%na, zmff)
      do iCv = 1, nCv
        ce1(iCv)=fce1(zeff(iCv))
        ce2(iCv)=fce2(zeff(iCv))
        cimp1(iCv)=fce1(zmff(iCv))
        cimp2(iCv)=fce2(zmff(iCv))
      enddo                                                            !srv 13.01.17 }
*   ..compute the Coulomb logarithm
      call b2tlnl (nCv, switch, switch%icase_ii, pl%te, pl%ti, 
     &             dv%ne, dv%lnlam)
*   ..compute zetap and zetae
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      czetae = t0*sqrt(me)
      zetap = czetap*(pl%ti/ev)*sqrt(pl%ti/ev)/dv%lnlam
      zetae = czetae*(pl%te/ev)*sqrt(pl%te/ev)/dv%lnlam
*   ..compute gti, gtn, gte
      !wdk gti, gtn and gte are now computed directly in cell centers
      !wdk (instead of faces as in structured code)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, gti)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%tn, wrkv, gtn)
      call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, gte)
*   ..compute ehxp
      if(switch%b2sigp_style.ne.2) then 
        call b2nxfx (nCv, nFc, nVx, switch, geo, mpg, qe, dv%ne, 
     &     pl%te, pl%po, ehxp, wrkf0, wrkf1)
      endif                                                            !srv 13.01.17 }
*   ..compute smqdu
      call b2nxdu (nCv, ns, pl%na, pl%ua, sr%smq, smqdu, wrk0, wrk1)
      if (switch%b2npmo_iout.ne.0) then
       call my_out_us(70,nCv,0,smqdu,'b2npmo_smqdu')
      endif
*   ..initialise
      call sfill (nCv*4, 0.0_R8, smt, 1)
      call sfill (nFc*2, 0.0_R8, flct, 1)
      call sfill (nFc*2, 0.0_R8, cvct, 1)
      call sfill (nCv, 0.0_R8, dv%resmt, 1)
      call sfill (nCv, 0.0e0_R8, ctcft, 1)                              !srv 06.04.07

* ..main computation
*   ..compute individual species corrections
!      do is = 0, ns-1
      do is = switch%nsmin, switch%nsmax-1
*    ..compute rob, pb
       call b2scopy (nCv, pl%na(:,is), 1, rob, 1)
       call b2sscal (nCv, am(is)*mp, rob, 1)
       call b2xppb (nCv, rt%rza(:,is), pl%na(:,is), pl%te, pl%ti, 
     &               pl%tn, is, pb)
*    ..compute extended friction term
       if(switch%b2sigp_style.eq.2) then                               !srv 13.01.17
        call b2sifrtf (nCv, nFc, nVx, ns, is, ismain, switch, geo,
     &   mpg, pl%na, pl%ua, dv%ue, pl%te, pl%ti, dv%ni, dv%ne, dv%ne2,
     &   rt%rza, rt%rz2, pl%po,
     &   zetap, zetae, gti, gte, ce1, ce2, cimp1, cimp2, zeff, 
     &   co%f_luc_sg, co%alfx_c, co%sigx_c, st_ext,
     &   smbfrea, smbfria, smbtfea, smbtfia, smbch, smbtf)             ! smbch and smbtf will be multiplied by vol*hz below !srv 05.07.17
        do k = 0, 3
         smbch(:,k) = geo%cvVol*geo%cvHz*smbch(:,k)
         smbtf(:,k) = geo%cvVol*geo%cvHz*smbtf(:,k)
        enddo
       else
        smbfria = 0.0_R8
        smbfrea = 0.0_R8
        smbtfia = 0.0_R8
        smbtfea = 0.0_R8
        call b2sifr_ (nCv, nFc, nVx, ns, is, switch, geo, mpg,
     &   pl%na, pl%ua, pl%te, pl%ti, dv%ni, dv%ne, dv%ne2, rt%rza,
     &   rt%rz2, ehxp, dv%lnlam, co%cthe, co%cthi, st_ext,
     &   smbch, smbtf, smpr(:,is), smpt(:,is), smfr(:,is),
     &   wrk0, gti, gte)
       endif

*    ..compute pressure gradient term
       call b2sigp (nCv, nFc, nVx, is, switch, geo, mpg, rt%rza(:,is),
     &   pl%na(:,is), pl%ti, pl%tn, dv%ne, pl%te, pl%po, smbgp, wrk0, 
     &   wrk1, wrk2)
*    ..compute flub
       call b2scopy (nFc*2, dv%fna_fcor(:,:,is), 1, flub, 1)
       call b2sscal (nFc*2, am(is)*mp, flub, 1)
*    ..add convective term due to new form of ion viscosity term only
       call b2nxfv (nCv, nFc, nVx, is, switch, geo, mpg,  
     &   co%cvsahz_cl(:,:,is), co%cvsahz(:,:,is), flubv)
!WG_TODO       if (balance_netcdf.ne.0) then
!WG_TODO         fmo_b2nxfv(:,:,:,:,is)=flubv !djm Jan2017
!WG_TODO       endif
       call b2saxpy (nFc*2, 1.0_R8, flubv, 1, flub, 1)                  !srv 11.09.09
*    ..compute flcb, cvcb
       call b2nxfc (nCv, nFc, nVx, is, switch, geo, mpg, rob, 
     &  pl%ua(:,is), flub, co%cvsahz(:,:,is), flcb, cvcb, wrkf0)

*    ..compute additional viscosity term
       call b2siav (nCv, nFc, nVx, ns, is, switch, geo, mpg, am(is)*mp,
     &   pl%ua(:,is), pl%te, pl%ti, pl%na(:,is), dv%ni(:,1), dv%ne2,
     &   rt%rz2, dv%lnlam, co%hcix_c, co%cvsahz_cl(:,:,is), smbvh,
     &   smbvv, co%kt_neo, wrk0, wrk1, wrk2)
*    ..compute centrifugal force term                                   !srv 09.07.01
       call b2sicf (nCv, nFc, nVx, is, am(is)*mp, switch,  geo, mpg,    !srv 09.07.01 28.01.02 22.07.05
     *   pl%na(:,is), pl%ua(:,is), smbcf, smcf(:,is), ctcfb, wrk0)
*    ..compute correcting source due to artificial anomalous electric conductivity
       call b2sian (nCv, nFc, is, switch, geo, mpg, dv%fchanml, smban)
*    ..compute total parallel momentum source
       call b2scopy (nCv*4, sr%smo(:,:,is), 1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, sr%smq(:,:,is), 1, smb, 1)
!WG_TODO       b2npmo_smaf(:,:,is)=smbch(:,:)+smbtf(:,:)                !srv 11.09.09 {
       call b2saxpy (nCv*4, 1.0_R8, smbch, 1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbtf, 1, smb, 1)                   !srv 11.09.09 }
!WG_TODO       b2npmo_smag(:,:,is)=smbgp(:,:)
       call b2saxpy (nCv*4, 1.0_R8, smbgp, 1, smb, 1)
!WG_TODO       b2npmo_smav(:,:,is)=smbvh(:,:)+smbvv(:,:)                !srv 03.07.08 {
       call b2saxpy (nCv*4, 1.0_R8, smbvh,  1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbvv,  1, smb, 1)
       call b2saxpy (nCv*4, 1.0_R8, smbcf,  1, smb, 1)                  !srv 03.07.08 }
       call b2saxpy (nCv*4, 1.0_R8, smban,  1, smb, 1)                  !srv 25.09.17
       call b2nxdv (nCv,is,switch,mpg,rob,sr%smq(:,:,is),smqdu,smb)     !srv ???
*    ..save final source for species (is)                               !srv 10.11.02
!WG_TODO       smof(:,:,is) = smb                                               !srv 10.11.02

cdjm Jan2017 update momentum source components for balance
!WG_TODO       if (balance_netcdf.ne.0) then
!WG_TODO         call update_balance_sources_momentum(nx,ny,is,ua(-1,-1,is),rob)
!WG_TODO       endif

*    ..compute residual
      if((solving.and.any(solvemo(is,0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then                  !srv 22.05.18
        call b2urmo (nCv, nFc, nVx, switch, geo, mpg, is, pl%ua(:,is),
     .    rob, smb, flub, flubv, co%cvsahz(:,:,is), dv%resmo(:,is),
     .    dv%fmo(:,:,is), dv%fhm(:,:,is), co%cvsahz_eff(:,:,is))       !srv 04.07.08
      endif                                                            !srv 22.05.18

       write (chns,'(i3.3)') is                                        !srv 18.03.10

       if (switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat .eq. 4) then                     !srv 12.07.02 {
         call my_out_us(70,nCv,0,dv%resmo(1,is),'b2npmo_resmo'//chns)
        call my_out_us(70,nFc,1,dv%fmo(1,0,is),'b2npmo_fmox'//chns)
        call my_out_us(70,nFc,1,dv%fmo(1,1,is),'b2npmo_fmoy'//chns)
         if(.not. is_neutral(is)) then
            call my_out_us(70,nFc,1,flubv(1,0),'b2npmo_flubvx'//chns)
            wrks=smbvv(:,0)+smbvv(:,1)*pl%ua(:,is)+
     &         smbvv(:,2)*pl%na(:,is)*mp*am(is)+
     &         smbvv(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
            call my_out_us(70,nCv,0,wrks,'b2npmo_smovv'//chns)
            wrks=smbvh(:,0)+smbvh(:,1)*pl%ua(:,is)+
     &         smbvh(:,2)*pl%na(:,is)*mp*am(is)+
     &         smbvh(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
            call my_out_us(70,nCv,0,wrks,'b2npmo_smovh'//chns)
         endif
         wrks=smb(:,0)+smb(:,1)*pl%ua(:,is)+
     &       smb(:,2)*pl%na(:,is)*mp*am(is)+
     &       smb(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrks,'b2npmo_smb'//chns)                  !srv 28.02.17
         wrks=smbcf(:,0)+smbcf(:,1)*pl%ua(:,is)+
     &       smbcf(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbcf(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrks,'b2npmo_smocf'//chns)
         wrks = geo%cvVol*geo%cvHz*(smbfria(:,0)+
     &        smbfria(:,1)*pl%ua(:,is)+      !srv 05.07.17 {
     &        smbfria(:,2)*pl%na(:,is)*mp*am(is)+
     &        smbfria(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is))
	 call my_out_us(70,nCv,0,wrks,'b2npmo_smofria'//chns)

         wrks = geo%cvVol*geo%cvHz*(smbfrea(:,0)+
     &                             smbfrea(:,1)*pl%ua(:,is))
	 call my_out_us(70,nCv,0,wrks,'b2npmo_smofrea'//chns)
         wrks = geo%cvVol*geo%cvHz*smbtfia(:,0)
	 call my_out_us(70,nCv,0,wrks,'b2npmo_smotfia'//chns)
         wrks = geo%cvVol*geo%cvHz*smbtfea(:,0)
	 call my_out_us(70,nCv,0,wrks,'b2npmo_smotfea'//chns)
         call my_out_us(70,nCv,0,dv%ue,'b2npmo_ue')                        !srv 05.07.17 }
         wrks=smban(:,0)+smban(:,1)*pl%ua(:,is)+
     &       smban(:,2)*pl%na(:,is)*mp*am(is)+
     &       smban(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
         call my_out_us(70,nCv,0,wrks,'b2npmo_smoan'//chns)               !srv 25.09.17
       endif                                                           !srv 05.07.17

       if (switch%b2npmo_iout .ne. 0) then
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2npmo_ua_input'//chns)
        wrks=sr%smo(:,0,is)+sr%smo(:,1,is)*pl%ua(:,is)+
     &       sr%smo(:,2,is)*pl%na(:,is)*mp*am(is)+
     &       sr%smo(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smo'//chns)                  !srv 28.02.17
        wrks=smbch(:,0)+smbch(:,1)*pl%ua(:,is)+
     &       smbch(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbch(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smoch'//chns)
        wrks=smbtf(:,0)+smbtf(:,1)*pl%ua(:,is)+
     &       smbtf(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbtf(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smotf'//chns)
        wrks=smbgp(:,0)+smbgp(:,1)*pl%ua(:,is)+
     &       smbgp(:,2)*pl%na(:,is)*mp*am(is)+
     &       smbgp(:,3)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrks,'b2npmo_smogp'//chns)
       endif                                                           !srv 12.07.02 }

       if (solving.and.any(solvemo(is,0:mpg%nnreg(0)))) then           !srv 06.07.10 22.05.18
         b2npmo_solvedum(0:mpg%nnreg(0))=solvemo(is,0:mpg%nnreg(0))
c     ..compare matrices
         new_matrix = .false.
         do ireg = 0, mpg%nnreg(0)
          new_matrix = new_matrix.or.
     &     solvemo(is,ireg).neqv.last_solve_9(ireg)
         end do
         if (new_matrix.or..true.) call restart_ma28_for_us
*     ..compute correction
         pb1 = pb !csc to avoid aliasing in Tapenade
         rob1 = rob !csc to avoid aliasing in Tapenade
         call b2usmo (nCv, nFc, nVx, is, switch, geo, mpg, mpg%nnreg(0),
     &     b2npmo_solvedum, itcnt, switch%b2npmo_rxg, rob, rob1, pb,
     &     pb1, pl%ua(:,is), smb,
     &     flcb, cvcb, dv%resmo(:,is),  ctcfb, dv%corua(:,is),
     &     dv%pcca(:,:,is), aa, wrk1, 'b2usmo_'//chns)
         last_solve_9(0:mpg%nnreg(0)) = solvemo(is,0:mpg%nnreg(0))
csw 21feb2014 added this else branch
       else
         dv%pcca(:,:,is)  = 0.0_R8
!         dv%resmo(:,is) = 0.0_R8       ! it is useful to know resmo    !srv 22.05.18
         dv%corua(:,is) = 0.0_R8
       endif
*    ..compute contribution to total residual

csw 21feb2014 and do this only for the solved species, i.e. added next branch
      if((solving.and.any(solvemo(0:ns-1,0:mpg%nnreg(0)))
     &   .and.any(solvemt(0:mpg%nnreg(0))))) then !srv 22.05.18 resmt is non-zero when all equations are solved
!     &   .and.all(solvemt(0:mpg%nnreg(0)))).or.get_residuals.eq.1) then !srv 22.05.18
        call b2scopy (nCv*4, sr%smo(:,:,is), 1, smb, 1)
        if(switch%b2sigp_style.eq.2) then                                             !srv 13.01.17 {
*       (in this case thermal force and the friction are taken into account)
         call b2saxpy (nCv*4, 1.0_R8, smbtf, 1, smb, 1)                 !srv 08.10.14
         call b2saxpy (nCv*4, 1.0_R8, smbch, 1, smb, 1)
        else
*       (now the friction cancels out)
        endif                                                          !srv 13.01.17 }
        call b2saxpy (nCv*4, 1.0_R8, smbgp, 1, smb, 1)
        call b2saxpy (nCv*4, 1.0_R8, smbvh, 1, smb, 1)
        call b2saxpy (nCv*4, 1.0_R8, smbvv, 1, smb, 1)                  !srv 03.07.08
        call b2saxpy (nCv*4, 1.0_R8, smbcf, 1, smb, 1)                  !srv 09.07.01
        call b2saxpy (nCv*4, 1.0_R8, smban, 1, smb, 1)                  !srv 25.09.17
        call b2urmo (nCv, nFc, nVx, switch, geo, mpg,  is, pl%ua(:,is),
     .     rob, smb, flub, flubv, co%cvsahz(:,:,is), wrk0, wrkf, wrkff,
     .    wrkfff)
*    ..accumulate residual
        call b2saxpy (nCv, 1.0_R8, wrk0, 1, dv%resmt, 1)
*    ..accumulate smb
        call b2nxst (nCv, rob, smb, ro, smt)
*    ..accumulate transport coefficients
        call b2saxpy (nFc*2, 1.0_R8, flcb, 1, flct, 1)
        call b2saxpy (nFc*2, 1.0_R8, cvcb, 1, cvct, 1)
        call b2saxpy (nCv, 1.0_R8, ctcfb, 1, ctcft, 1)                  !srv 06.04.07
       endif                                                            !srv 22.05.18
      enddo

csw 21feb2014, execute simultaneous corrections always
      if (solving.and.any(solvemo(0:ns-1,0:mpg%nnreg(0)))
     &                  .and.any(solvemt(0:mpg%nnreg(0)))) then           !srv 06.07.10 22.05.18
*     ..compute simultaneous correction
        call sfill (nCv, 0.0_R8, wrk0, 1)
c     ..compare 9-pt stencil matrices
        new_matrix = .false.
        do ireg = 0, mpg%nnreg(0)
         new_matrix = new_matrix.or.
     &    solvemt(ireg).neqv.last_solve_9(ireg)
        end do
        if (new_matrix.or..true.) call restart_ma28_for_us
*     ..compute correction
        call b2usmo (nCv, nFc, nVx, ns, switch, geo, mpg, mpg%nnreg(0), 
     &    solvemt, itcnt, switch%b2npmo_rxg, ro, rz, pr, pz, wrk0,
     &    smt, flct, cvct,
     &    dv%resmt, ctcft, dv%corut, dv%pccm, aa, wrk1, 'b2usmo_total')
        last_solve_9(0:mpg%nnreg(0)) = solvemt(0:mpg%nnreg(0))
      else
        dv%pccm  = 0.0_R8
!        dv%resmt = 0.0_R8                 ! it is useful to know resmt     !srv 22.05.18
        dv%corut = 0.0_R8
      endif                                                            !srv 06.07.10

*   ..apply correction
      if (rxf.ne.0.0_R8) then
       do is = 0, ns-1
        do iCv = 1, mpg%nCi ! internal cells
         ireg = mpg%cvReg(iCv)
         if (solvemo(is,ireg)) then
          if(solvemt(ireg)) then                                   !srv 22.05.18
            pl%ua(iCv,is) = pl%ua(iCv,is)+
     &        rxf*med(dv%corua(iCv,is),dv%corut(iCv),
     &        dv%corua(iCv,is)+dv%corut(iCv))
          else                                                     !srv 22.05.18  {
            pl%ua(iCv,is) = pl%ua(iCv,is) + rxf*dv%corua(iCv,is)
          endif                                                    !srv 22.05.18  }
         endif
        enddo
        do iCv = mpg%nCi+1, nCv ! guard cells
         ireg = mpg%cvReg(iCv)
         if (solvemo(is,ireg)) then
           pl%ua(iCv,is) = pl%ua(iCv,is) + rxf*dv%corua(iCv,is)
         endif
        enddo
       enddo
      endif

      vcount = 0
      if( switch%ion_vlct_restrict .ne. 0 ) then                        !iyv 19.07.12 { !srv 17.01.14
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            do iCv = 1, nCv
              cst=sqrt(pz(iCv)/rz(iCv))                                 !iyv 03.02.14
              if( abs(pl%ua(iCv,is)).gt.
     .            switch%ion_vlct_restrict_M*cst ) then
                vcount(is) = vcount(is) + 1
                if (switch%vlct_diagno.ge.2) then
                  write(6,'(a,a,i3,a,i4,a,i4,a)')
     .             'Parallel velocity restriction applied ',
     .             'to species ',is,' in cell (',iCv,'):'
                  write(6,'(a,1p,1e12.5,a,1e12.5,a,1e12.5)')
     .             'Mach ',switch%ion_vlct_restrict_M,': ',
     .                     switch%ion_vlct_restrict_M*cst,
     .              ' Ua:',pl%ua(iCv,is)
                end if
                pl%ua(iCv,is)= sign(switch%ion_vlct_restrict_M*cst,
     .                              pl%ua(iCv,is))                      !iyv 03.02.14
              endif
            enddo
          endif
        enddo
        if (switch%vlct_diagno.ge.1) then
          do is = 0, ns-1
            if (vcount(is).gt.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_M, 
     .         ' was applied ',vcount(is),
     .         ' times for species ',is
            else if (vcount(is).eq.1) then
              write(6,'(a,a,1e12.5,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach ', switch%ion_vlct_restrict_M,
     .         ' was applied ',vcount(is),
     .         ' time for species ',is
            end if
          end do
        endif
      endif

*   ..correct pcca
      call b2nxcm (nCv, ns, dv%pccm, dv%pcca)
      if (switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat .eq. 4) then !srv 12.07.02 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is                                       !srv 18.05.02
        call my_out_us(70,nCv,0,pl%ua(1,is),'b2npmo_ua'//chns)
       enddo
      endif                                                            !srv 12.07.02 }
      if (switch%b2npmo_iout.ne.0) then
       call my_out_us(70,nCv,0,dv%pccm(1,0),'b2npmo_pccmx')            !srv 11.01.13
       call my_out_us(70,nCv,0,dv%pccm(1,1),'b2npmo_pccmy')            !srv 11.01.13
      endif                                                            !srv 12.07.02 }
*   ..set error parameter
*     (at present all errors cause an abort through xerrab)
      ierr = 0
* ..return
      ncall_b2npmo = ncall_b2npmo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npmo

      end subroutine b2npmo
