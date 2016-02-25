*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo (nx, ny, ns, ismain, itcnt, solving, rxf,      !srv 13.01.17
     &  qz, vol, hx, hy, hz, bb, pbs, pbshz, lnlam,
     &  na, te, ti, po, ue, fna_fcor,                                  !srv 09.01.01 01.07.09 05.07.17
     &  cvsahz, cvsa, cvsa_cl, cvsahz_cl, kt_neo,
     &  smo, smq, cthe, cthi, chci, ua, smof, wadia, vaecrb,
     &  resmo, resmt, corua, corut, pcca, pccm, ierr)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_sources
     1 , only : b2npmo_smaf, b2npmo_smag, b2npmo_smav,
     2          smpr, smpt, smfr, smcf
      use b2mod_balance !djm Jan2017
     1 , only : update_balance_sources_momentum, fmo_b2nxfv,
     2          balance_netcdf
      use b2mod_plasma
     & , only : fmo, fchanml                                           !srv 25.09.17
      use b2mod_numerics_namelist
      use b2mod_transport_fun                                          !srv 13.01.17
      use b2mod_b2cmpa
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, itcnt                                !srv 13.01.17
      real (kind=R8) ::
     *  rxf, 
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),             !srv 09.01.01
     &  qz(-1:nx,-1:ny,0:1),
     &  vol(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), pbs(-1:nx,-1:ny,0:1),
     &  pbshz(-1:nx,-1:ny,0:1), lnlam(-1:nx,-1:ny), ue(-1:nx,-1:ny),   !srv 01.07.09 05.07.17
     &  na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny), fna_fcor(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsahz(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  cvsa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1), 
     &  cvsahz_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1),                         !srv 09.01.01
     &  smo(-1:nx,-1:ny,0:3,0:ns-1), smq(-1:nx,-1:ny,0:3,0:ns-1),
     &  chci(-1:nx,-1:ny,0:1,0:1),                                     !srv 10.11.02
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1)
      logical solving
*   ..input/output arguments
      real (kind=R8) ::
     *  ua(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *  kt_neo(-1:nx,-1:ny),
     *  resmo(-1:nx,-1:ny,0:ns-1), resmt(-1:nx,-1:ny),
     &  corua(-1:nx,-1:ny,0:ns-1), corut(-1:nx,-1:ny),
     &  pcca(-1:nx,-1:ny,0:1,0:1,0:ns-1), pccm(-1:nx,-1:ny,0:1,0:1),
     &  smof(-1:nx,-1:ny,0:3,0:ns-1)
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
      integer ncall, n2, ix, iy, is, diagno, ireg, nsmin, nsmax, k
      integer, save :: iout = 0, b2sifr_ver = 1, iout_b2wdat = 0       !srv 04.07.08 26.10.08
      integer, save :: icase = 0                                       !srv 13.01.17
      integer, save :: ion_vlct_restrict = 0                           !srv 17.01.14
      integer, save :: vlct_diagno = 1
      integer :: vcount(0:ns-1)
      integer, save :: style = 2                                       !srv 08.10.14 13.01.17
      character chns*3, chk*1                                          !srv 04.07.08
      real (kind=R8) ::
     *  rxg, t0, t1, t2, czetap, czetae                                !srv 13.01.17
      real (kind=R8) ::
     *  smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3),                !srv 11.09.09
     *  smbgp(-1:nx,-1:ny,0:3),
     *  smbvh(-1:nx,-1:ny,0:3), smbvv(-1:nx,-1:ny,0:3),                !srv 23.07.99 04.10.99
     *  smbcf(-1:nx,-1:ny,0:3), smban(-1:nx,-1:ny,0:3),                !srv 09.07.01 25.09.17
     1  pb(-1:nx,-1:ny), cvct(-1:nx,-1:ny,0:1,0:1),
     2  rob(-1:nx,-1:ny), flub(-1:nx,-1:ny,0:1,0:1),
     *  flubv(-1:nx,-1:ny,0:1,0:1),                                    !srv 11.09.09
     3  wrkf(-1:nx,-1:ny,0:1,0:1), aa5(-1:nx,-1:ny,-2:2),
     4  smb(-1:nx,-1:ny,0:3), flcb(-1:nx,-1:ny,0:1,0:1),
     5  cvcb(-1:nx,-1:ny,0:1,0:1), ro(-1:nx,-1:ny),
     6  pr(-1:nx,-1:ny), ne2(-1:nx,-1:ny),
     7  ni(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     8  ehxp(-1:nx,-1:ny), smt(-1:nx,-1:ny,0:3), cst,                  !srv 17.01.14
     9  flct(-1:nx,-1:ny,0:1,0:1), wrk1(-1:nx,-1:ny), wrks(-1:nx,-1:ny),
     1  smqdu(-1:nx,-1:ny), wrk0(-1:nx,-1:ny), wrk2(-1:nx,-1:ny),
     2  ctcfb(-1:nx,-1:ny), ctcft(-1:nx,-1:ny),                        !srv 06.04.07
     3  zeff(-1:nx,-1:ny), zmff(-1:nx,-1:ny),                          !srv 13.01.17 {
     4  ce1(-1:nx,-1:ny), ce2(-1:nx,-1:ny),
     5  cimp1(-1:nx,-1:ny), cimp2(-1:nx,-1:ny),
     6  zetap(-1:nx,-1:ny), zetae(-1:nx,-1:ny),
     7  gti(-1:nx,-1:ny), gte(-1:nx,-1:ny)                             !srv 13.01.17 }
      real (kind=R8) ::
     *  smbfria(-1:nx,-1:ny,0:3), smbfrea(-1:nx,-1:ny,0:3),
     *  smbtfia(-1:nx,-1:ny,0:3), smbtfea(-1:nx,-1:ny,0:3)             !srv 13.01.17 05.07.17
      logical new_matrix, solvedum(0:NNREGMAX), main_lightest
      save ncall, rxg, diagno, solvedum, nsmin, nsmax
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *  med
      external subini, subend, xertst, ipgetr, sfill, b2scopy,
     &  b2sscal, b2saxpy
      external b2xvsg, b2xvff, b2xpne, b2xpro, b2xppr,
     &  b2xppb, b2nxdu, b2nxdv, b2nxfx, b2nxfm, b2nxfc, b2nxst,
     &  b2sifr, b2sigp, b2siav, b2sicf, b2nxcm, b2urmo, b2usmo
      med(t0,t1,t2) = min(max(t0,t1),max(t1,t2),max(t2,t0))
*   ..initialisation
      data ncall/0/, rxg/1.0e6_R8/, diagno/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        nsmin=0
        nsmax=ns
        call ipgetr ('b2npmo_rxg', rxg)
        call xertst (0.0_R8.lt.rxg, 'faulty internal parameter rxg')
        call ipgeti ('b2npmo_diagno', diagno)
        call ipgeti ('b2npmo_iout', iout)                              !srv 04.07.08
        call ipgeti ('b2npmo_ion_vlct_restrict', ion_vlct_restrict)    !iyv 17.01.14
        call ipgeti ('b2npmo_vlct_diagno', vlct_diagno)
        call ipgeti ('b2wdat_iout', iout_b2wdat)                       !iys 17.04.17
        call ipgeti ('b2npmo_b2sifr_', b2sifr_ver)                     !srv 04.07.08
        call ipgeti ('b2news_nsmin', nsmin)
        if(nsmin.gt.ns) then
           nsmin=ns
        elseif(nsmin.lt.0) then
           nsmin=0
        endif
        call ipgeti ('b2news_nsmax', nsmax)
        if(nsmax.gt.ns) then
           nsmax=ns
        elseif(nsmax.lt.0) then
           nsmax=0
        endif
        write(*,*) 'nsmin, nsmax = ', nsmin, nsmax
        call ipgeti ('b2sigp_style', style)                            !srv 08.10.14
        main_lightest = .true.
        do is = 0, ns-1
          main_lightest = main_lightest .and.
     .     nint(am(is)).ge.nint(am(ismain))
        end do
        if (.not.main_lightest .and. style.eq.2) then
         write(0,'(2a)') 'Warning: b2sigp_style.eq.2 not appropriate ',
     .    'for plasmas where the main species is not the lightest one!'
        endif
        call ipgeti ('b2tlnl_ii', icase)                               !srv 13.01.17
        call xertst (0.eq.icase.or.1.eq.icase,'faulty input b2tlnl_ii') 
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8, 
     . 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of na, te, ti
       do is = 0, ns-1
        call b2xvsg (n2, na(-1,-1,is), 1, 'na', '.gt.')
       enddo
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test edge values of fna_fcor
       do is = 0, ns-1
        call b2xvff (nx, ny, fna_fcor(-1,-1,0,0,is), 'fna_fcor')
       enddo
*    ..test sign of cvsa
       call b2xvsg (n2*4*ns, cvsa, 1, 'cvsa', '.ge.')
*    ..test edge values of cvsa
       do is = 0, ns-1
        call b2xvff (nx, ny, cvsa(-1,-1,0,0,is), 'cvsa')
       enddo
*    ..test sign of smo, smq
       do is = 0, ns-1
        call b2xvsg (n2, smo(-1,-1,1,is), 1, 'smo 1', '.le.')
        call b2xvsg (n2, smo(-1,-1,3,is), 1, 'smo 3', '.le.')
        call b2xvsg (n2, smq(-1,-1,1,is), 1, 'smq 1', '.le.')
        call b2xvsg (n2, smq(-1,-1,3,is), 1, 'smq 3', '.le.')
       enddo
      endif
*   ..compute ne, ne2, ro, pr
      call b2xpni (nx, ny, ns, na, ni)
      call b2xpne (nx, ny, ns, rza, na, ne_ext, ne)
      call b2xpne (nx, ny, ns, rz2, na, ne2_ext, ne2)
      call b2xpro (nx, ny, ns, mp, am, na, ro)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
      call b2xzef (nx, ny, ns, rz2, na, ne, zeff)                      !srv 13.01.17 {
      call b2xzmf (nx, ny, ns, ismain, rz2, na, zmff)
      do iy = -1, ny
        do ix = -1, nx
	  ce1(ix,iy)=fce1(zeff(ix,iy))
	  ce2(ix,iy)=fce2(zeff(ix,iy))
          cimp1(ix,iy)=fce1(zmff(ix,iy))
          cimp2(ix,iy)=fce2(zmff(ix,iy))
        enddo
      enddo                                                            !srv 13.01.17 }
*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)
*   ..compute zetap and zetae
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      czetap = t0*sqrt(mp)
      czetae = t0*sqrt(me)
      zetap = czetap*(ti/ev)*sqrt(ti/ev)/lnlam
      zetae = czetae*(te/ev)*sqrt(te/ev)/lnlam
*   ..compute gti,gte
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         gti(ix,iy) = 0.0_R8
         gte(ix,iy) = 0.0_R8
        else
         gti(ix,iy) = (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
         gte(ix,iy) = (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
     &     (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8
        endif
       enddo
      enddo
*   ..compute ehxp
      if(style.ne.2) then 
        call b2nxfx (nx, ny, qe, vol, hx, bb, ne, te, po,
     &   ehxp, wrk0, wrk1)
      endif                                                            !srv 13.01.17 }
*   ..compute smqdu
      call b2nxdu (nx, ny, ns, na, ua, smq, smqdu, wrk0, wrk1)
      if (iout.ne.0) then
       call my_out(70,nx,ny,smqdu,'b2npmo_smqdu')
      endif
*   ..initialise
      call sfill (n2*4, 0.0_R8, smt, 1)
      call sfill (n2*4, 0.0_R8, flct, 1)
      call sfill (n2*4, 0.0_R8, cvct, 1)
      call sfill (n2, 0.0_R8, resmt, 1)
      call sfill (n2, 0.0e0_R8, ctcft, 1)                              !srv 06.04.07

* ..main computation
*   ..compute individual species corrections
!      do is = 0, ns-1
      do is = nsmin, nsmax-1
*    ..compute rob, pb
       call b2scopy (n2, na(-1,-1,is), 1, rob, 1)
       call b2sscal (n2, am(is)*mp, rob, 1)
       call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute extended friction term
       if(style.eq.2) then                                             !srv 13.01.17
         call b2sifrtf (nx, ny, ns, is, ismain, hx, bb, na, ua,        !srv 05.07.17
     &    ue, te, ti, ni, ne, ne2, zetap, zetae, gti, gte, ce1, ce2,   !srv 05.07.17
     &    cimp1, cimp2, zeff, smbfrea, smbfria, smbtfea, smbtfia,      !srv 05.07.17
     &    smbch, smbtf)  ! smbch and smbtf will be multiplied by vol*hz below !srv 05.07.17
         do iy =-1,ny
           do ix =-1,nx
             smbch(ix,iy,0:3) = vol(ix,iy)*hz(ix,iy)*smbch(ix,iy,0:3)
             smbtf(ix,iy,0:3) = vol(ix,iy)*hz(ix,iy)*smbtf(ix,iy,0:3)
           enddo
         enddo
       else
        if(b2sifr_ver.eq.0) then                                       !srv 26.10.08 {
         call b2sifr (nx, ny, ns, is, vol, hx, hz, bb,                 !srv 09.01.01
     &    na, ua, te, ti, ni, ne, ne2, ehxp, lnlam, cthe, cthi,        !srv 11.09.09
     &    smbch, smbtf,
     &    smpr(-1,-1,is), smpt(-1,-1,is), smfr(-1,-1,is), 
     &    wrk0, wrk1, wrk2)                                            !srv 02.07.08 01.07.09
        else
         call b2sifr_ (nx, ny, ns, is, hx, hz, vol, bb,                !srv 09.01.01 18.10.02
     &    na, ua, te, ti,ni, ne, ne2, ehxp, lnlam, cthe, cthi,
     &    smbch, smbtf,                                                !srv 11.09.09
     &    smpr(-1,-1,is), smpt(-1,-1,is), smfr(-1,-1,is), 
     &    wrk0, wrk1, wrk2)                                            !srv 01.07.09
        endif                                                          !srv 26.10.08 }
       endif
*    ..compute pressure gradient term
       call b2sigp (nx, ny, is, rza(-1,-1,is), vol, hx, hz, bb,  
     &   na(-1,-1,is), ti, ne, te, po, smbgp, wrk0, wrk1, wrk2)        !srv 02.07.08
*    ..compute flub
       call b2scopy (n2*4, fna_fcor(-1,-1,0,0,is), 1, flub, 1)
       call b2sscal (n2*4, am(is)*mp, flub, 1)
!       call b2nxfm (nx, ny, am(is)*mp, fna_fcor(-1,-1,0,0,is), flub)
*    ..add convective term due to new form of ion viscosity term only
       call b2nxfv (nx, ny, is, hz, cvsahz_cl(-1,-1,0,0,is), flubv)    !srv 01.07.05 11.09.09 16.10.17
       if (balance_netcdf.ne.0) then
         fmo_b2nxfv(:,:,:,:,is)=flubv !djm Jan2017
       endif
       call b2saxpy (n2*4, 1.0_R8, flubv, 1, flub, 1)                  !srv 11.09.09
*    ..compute flcb, cvcb
       call b2nxfc (nx, ny, is, vol, pbshz, rob, 
     &   ua(-1,-1,is), flub, cvsahz(-1,-1,0,0,is), flcb, cvcb, wrk0)   !srv 09.01.01

*    ..compute additional viscosity term
       call b2siav (nx, ny, ns, is, am(is)*mp, hx, hy, hz, vol, 
     &   qz, bb, ua(-1,-1,is), wadia(-1,-1,0,is), vaecrb(-1,-1,0,is),  
     &   te, ti, na(-1,-1,is), po, ni(-1,-1,1), ne2, lnlam, chci,      !srv 01.07.09
     &   cvsa(-1,-1,0,0,is), cvsa_cl(-1,-1,0,0,is), 
     &   cvsahz(-1,-1,0,0,is), 
     &   cvsahz_cl(-1,-1,0,0,is), smbvh, smbvv, kt_neo, 
     &   wrk0, wrk1, wrk2)       !srv 23.07.99 04.10.99 17.01.01 09.01.01 31.07.02 01.07.05
*    ..compute centrifugal force term                                  !srv 09.07.01
       call b2sicf (nx, ny, is, am(is)*mp, bb, vol, hx, hz, 
     *     na(-1,-1,is), ua(-1,-1,is), smbcf, smcf(-1,-1,is), ctcfb, 
     *     wrk0)                                                       !srv 09.07.01 31.07.02 22.07.05} 06.04.07
*    ..compute correcting source due to artificial anomalous electric conductivity
       call b2sian (nx, ny, is, bb, vol, hx, hy, hz, qz, fchanml,      !srv 25.09.17
     &           smban)
*    ..compute total parallel momentum source
       call b2scopy (n2*4, smo(-1,-1,0,is), 1, smb, 1)
       call b2saxpy (n2*4, 1.0_R8, smq(-1,-1,0,is), 1, smb, 1)
       b2npmo_smaf(:,:,:,is)=smbch(:,:,:)+smbtf(:,:,:)                 !srv 11.09.09 {
       call b2saxpy (n2*4, 1.0_R8, smbch, 1, smb, 1)
       call b2saxpy (n2*4, 1.0_R8, smbtf, 1, smb, 1)                   !srv 11.09.09 }
       b2npmo_smag(:,:,:,is)=smbgp(:,:,:)
       call b2saxpy (n2*4, 1.0_R8, smbgp, 1, smb, 1)
       b2npmo_smav(:,:,:,is)=smbvh(:,:,:)+smbvv(:,:,:)                 !srv 03.07.08 {
       call b2saxpy (n2*4, 1.0_R8, smbvh,  1, smb, 1)
       call b2saxpy (n2*4, 1.0_R8, smbvv, 1, smb, 1)
       call b2saxpy (n2*4, 1.0_R8, smbcf,  1, smb, 1)                  !srv 03.07.08 }
       call b2saxpy (n2*4, 1.0_R8, smban,  1, smb, 1)                  !srv 25.09.17
       call b2nxdv (nx, ny, is, rob, smq(-1,-1,0,is), smqdu, smb)      !srv ???
*    ..save final source for species (is)                              !srv 10.11.02
       smof(:,:,:,is) = smb                                            !srv 10.11.02

cdjm Jan2017 update momentum source components for balance
       if (balance_netcdf.ne.0) then
         call update_balance_sources_momentum(nx,ny,is,ua(-1,-1,is),rob)
       endif

*    ..compute residual
       call b2urmo (nx, ny, is, ua(-1,-1,is), rob, smb,
     &   flub, cvsahz(-1,-1,0,0,is), resmo(-1,-1,is), wrkf)              !srv 04.07.08
       fmo(:,:,:,:,is)=wrkf

       write (chns,'(i3.3)') is                                        !srv 18.03.10

       if (iout.ne.0 .or. iout_b2wdat .eq. 4) then                     !srv 12.07.02 {
         call my_out(70,nx,ny,resmo(-1,-1,is),'b2npmo_resmo'//chns)
         call my_out(70,nx,ny,fmo(-1,-1,0,0,is),'b2npmo_fmox'//chns)
         call my_out(70,nx,ny,fmo(-1,-1,1,1,is),'b2npmo_fmoy'//chns)
         if(.not. is_neutral(is)) then
            call my_out(70,nx,ny,flubv(-1,-1,0,0),'b2npmo_flubvx'//chns)
            wrks=smbvv(:,:,0)+smbvv(:,:,1)*ua(:,:,is)+
     &         smbvv(:,:,2)*na(:,:,is)*mp*am(is)+
     &         smbvv(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
            call my_out(70,nx,ny,wrks,'b2npmo_smovv'//chns)
            wrks=smbvh(:,:,0)+smbvh(:,:,1)*ua(:,:,is)+
     &         smbvh(:,:,2)*na(:,:,is)*mp*am(is)+
     &         smbvh(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
            call my_out(70,nx,ny,wrks,'b2npmo_smovh'//chns)
         endif
         wrks=smb(:,:,0)+smb(:,:,1)*ua(:,:,is)+
     &       smb(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smb(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
         call my_out(70,nx,ny,wrks,'b2npmo_smb'//chns)                  !srv 28.02.17
         wrks=smbcf(:,:,0)+smbcf(:,:,1)*ua(:,:,is)+
     &       smbcf(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbcf(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
         call my_out(70,nx,ny,wrks,'b2npmo_smocf'//chns)
         wrks = vol*hz*(smbfria(:,:,0)+smbfria(:,:,1)*ua(:,:,is)+      !srv 05.07.17 {
     &        smbfria(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbfria(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is))
	 call my_out(70,nx,ny,wrks,'b2npmo_smofria'//chns)

         wrks = vol*hz*(smbfrea(:,:,0)+smbfrea(:,:,1)*ua(:,:,is))
	 call my_out(70,nx,ny,wrks,'b2npmo_smofrea'//chns)
         wrks = vol*hz*smbtfia(:,:,0)
	 call my_out(70,nx,ny,wrks,'b2npmo_smotfia'//chns)
         wrks = vol*hz*smbtfea(:,:,0)
	 call my_out(70,nx,ny,wrks,'b2npmo_smotfea'//chns)
         call my_out(70,nx,ny,ue,'b2npmo_ue')                          !srv 05.07.17 }
         wrks=smban(:,:,0)+smban(:,:,1)*ua(:,:,is)+
     &       smban(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smban(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
         call my_out(70,nx,ny,wrks,'b2npmo_smoan'//chns)               !srv 25.09.17
       endif                                                           !srv 05.07.17

       if (iout .ne. 0) then
        call my_out(70,nx,ny,ua(-1,-1,is),'b2npmo_ua_input'//chns)
        wrks=smo(:,:,0,is)+smo(:,:,1,is)*ua(:,:,is)+
     &       smo(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &       smo(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrks,'b2npmo_smo'//chns)                  !srv 28.02.17
        wrks=smbch(:,:,0)+smbch(:,:,1)*ua(:,:,is)+
     &       smbch(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbch(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrks,'b2npmo_smoch'//chns)
        wrks=smbtf(:,:,0)+smbtf(:,:,1)*ua(:,:,is)+
     &       smbtf(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbtf(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrks,'b2npmo_smotf'//chns)
        wrks=smbgp(:,:,0)+smbgp(:,:,1)*ua(:,:,is)+
     &       smbgp(:,:,2)*na(:,:,is)*mp*am(is)+
     &       smbgp(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrks,'b2npmo_smogp'//chns)
       endif                                                           !srv 12.07.02 }

       if (solving) then                                               !srv 06.07.10
        solvedum(0:nnreg(0))=solvemo(is,0:nnreg(0))
c    ..compare 5-pt stencil matrices
        new_matrix = .false.
        do ireg = 0, nnreg(0)
         new_matrix = new_matrix.or.
     &    solvemo(is,ireg).neqv.last_solve_5(ireg)
        end do
        if (new_matrix) call restart_ma28
*    ..compute correction
        call b2usmo (nx, ny, is, nnreg(0), solvedum, itcnt, rxg,       !srv 04.07.08
     &   vol, hx, hz, bb, rob, pb, ua(-1,-1,is),                       !srv 09.01.01
     &   smb, flcb, cvcb, resmo(-1,-1,is), ctcfb, corua(-1,-1,is),     !srv 06.04.07
     &   pcca(-1,-1,0,0,is), aa5, wrk1, 'b2usmo_'//chns)
        last_solve_5(0:nnreg(0)) = solvemo(is,0:nnreg(0))
csw 21feb2014 added this else branch
       else
         pcca(:,:,:,:,is)  = 0.0_R8
         resmo(:,:,is) = 0.0_R8
         corua(:,:,is) = 0.0_R8
       endif
*    ..compute contribution to total residual

csw 21feb2014 and do this only for the solved species, i.e. added next branch
       if (solving) then
        call b2scopy (n2*4, smo(-1,-1,0,is), 1, smb, 1)
        if(style.eq.2) then                                            !srv 13.01.17 {
*       (in this case thermal force and the friction are taken into account)
         call b2saxpy (n2*4, 1.0_R8, smbtf, 1, smb, 1)                 !srv 08.10.14
         call b2saxpy (n2*4, 1.0_R8, smbch, 1, smb, 1)
        else
*       (now the friction cancels out)
        endif                                                          !srv 13.01.17 }
        call b2saxpy (n2*4, 1.0_R8, smbgp, 1, smb, 1)
        call b2saxpy (n2*4, 1.0_R8, smbvh, 1, smb, 1)
        call b2saxpy (n2*4, 1.0_R8, smbvv, 1, smb, 1)                  !srv 03.07.08
        call b2saxpy (n2*4, 1.0_R8, smbcf, 1, smb, 1)                  !srv 09.07.01
        call b2saxpy (n2*4, 1.0_R8, smban, 1, smb, 1)                  !srv 25.09.17
        call b2urmo (nx, ny, is, ua(-1,-1,is), rob, smb,               !srv 09.01.01
     &   flub, cvsahz(-1,-1,0,0,is), wrk0, wrkf)
*    ..accumulate residual
        call b2saxpy (n2, 1.0_R8, wrk0, 1, resmt, 1)
*    ..accumulate smb
        call b2nxst (nx, ny, rob, smb, ro, smt)
*    ..accumulate transport coefficients
        call b2saxpy (n2*4, 1.0_R8, flcb, 1, flct, 1)
        call b2saxpy (n2*4, 1.0_R8, cvcb, 1, cvct, 1)
        call b2saxpy (n2, 1.0_R8, ctcfb, 1, ctcft, 1)                  !srv 06.04.07
       endif
      enddo

csw 21feb2014, execute simultaneous corrections always
c      if (solving.and.solvemt(0:nnreg(0))) then                       !srv 06.07.10
      if (solving) then
*   ..compute simultaneous correction
       call sfill (n2, 0.0_R8, wrk0, 1)
       new_matrix = .false.
       do ireg = 0, nnreg(0)
         new_matrix = new_matrix.or.
     &    solvemt(ireg).neqv.last_solve_5(ireg)
       end do
       if (new_matrix) call restart_ma28
       call b2usmo (nx, ny, ns, nnreg(0), solvemt, itcnt, rxg,         !srv 04.07.08
     &  vol, hx, hz, bb, ro, pr, wrk0, smt, flct, cvct,                !srv 09.01.01
     &  resmt, ctcft, corut, pccm, aa5, wrk1, 'b2usmo_total')          !srv 06.04.07
       last_solve_5(0:nnreg(0)) = solvemt(0:nnreg(0))
      else
       pccm  = 0.0_R8
       resmt = 0.0_R8
       corut = 0.0_R8
      endif                                                            !srv 06.07.10

*   ..apply correction
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
          ireg = region(ix,iy,0)
          if (solvemo(is,ireg)) then
            if((leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) .or.
     &       (bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) ) then  !srv 27.01.11
               ua(ix,iy,is) = ua(ix,iy,is) + corua(ix,iy,is)
            else
              ua(ix,iy,is) = ua(ix,iy,is)+
     &         rxf*med(corua(ix,iy,is),corut(ix,iy),
     &         corua(ix,iy,is)+corut(ix,iy))
            endif
          endif
        enddo
       enddo
      enddo

      vcount = 0
      if( ion_vlct_restrict .ne. 0 ) then                              !iyv 19.07.12 { !srv 17.01.14
        do iy = -1,ny
          do ix = -1,nx
            do is = 0, ns-1
              if (.not.is_neutral(is)) then
                cst=sqrt(pr(ix,iy)/ro(ix,iy))                          !iyv 03.02.14
!                cst=sqrt( (ti(ix,iy)+te(ix,iy)*rza(ix,iy,is))/
!     *                                              ( am(is)*mp ))
                if( abs(ua(ix,iy,is)) .gt. 3.0_R8*cst ) then
                  vcount(is) = vcount(is) + 1
                  if (vlct_diagno.ge.2) then
                    write(6,'(a,a,i3,a,i4,a,i4,a)')
     .               'Parallel velocity restriction applied ',
     .               'to species ',is,' in cell (',ix,',',iy,'):'
                    write(6,'(a,1p,1e12.5,a,1e12.5)')
     .               'Mach 3: ',3.0_R8*cst,' Ua:',ua(ix,iy,is)
                  end if
                  ua(ix,iy,is)= sign(3.0_R8*cst,ua(ix,iy,is))          !iyv 03.02.14
                endif
              endif
            enddo
          enddo
        enddo
        if (vlct_diagno.ge.1) then
          do is = 0, ns-1
            if (vcount(is).gt.1) then
              write(6,'(a,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach 3 was applied ',vcount(is),
     .         ' times for species ',is
            else if (vcount(is).eq.1) then
              write(6,'(a,a,i5,a,i3)')
     .         'In b2npmo, parallel velocity restriction ',
     .         'to Mach 3 was applied ',vcount(is),
     .         ' time for species ',is
            end if
          end do
        endif
      endif

*   ..correct pcca
      call b2nxcm (nx, ny, ns, pccm, pcca)
      if (iout.ne.0 .or. iout_b2wdat .eq. 4) then                      !srv 12.07.02 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is                                       !srv 18.05.02
        call my_out(70,nx,ny,ua(-1,-1,is),'b2npmo_ua'//chns)
       enddo
      endif                                                            !srv 12.07.02 }
      if (iout.ne.0) then
       call my_out(70,nx,ny,pccm(-1,-1,0,0),'b2npmo_pccmx')            !srv 11.01.13
       call my_out(70,nx,ny,pccm(-1,-1,1,1),'b2npmo_pccmy')            !srv 11.01.13
      endif                                                            !srv 12.07.02 }
*   ..set error parameter
*     (at present all errors cause an abort through xerrab)
      ierr = 0

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npmo

      end subroutine b2npmo
