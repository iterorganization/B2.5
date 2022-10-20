*-----------------------------------------------------------------------
*.specification

      subroutine b2trno (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &  switch, geo, mpg, pl, dv, rt, st_ext, co)
      use b2mod_types
!      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
!      use b2mod_neoclassical
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_transport_namelist
     & , only : parm_hce, parm_hci, parm_dna, parm_dpa, parm_vla,
     &          parm_vsa, parm_alf, parm_sig
      use b2mod_input_profile
     & , only : ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &          tr_ip_elm_count, no_div, addspec, nrr, nss, nkind_data,
     &          nkind_coeff, poloidal_scaling, transport_ip_filename,
     &          nscale
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2trno, ncall_b2tqna,
     &          ncall_transp_keps, nsdmax, cvregmax, ncall_b2tlh0
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1),
     *  ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      type (B2Coeff), intent (inout) :: co
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRNO computes the anomalous transport coefficients, including
*     geometric factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete anomalous transport
*     coefficients.
*     - cdna: particle transport proportional to grad na.
*     - cdpa: particle transport proportional to grad pa.
*     - cvla: anomalous velocity.
*     - cvsa: viscosity.
*     - chce: electron thermal conductivity.
*     - chve: electron thermal strange velocity.
*     - chci: atom thermal conductivity.
*     - chvi: atom thermal strange velocity.
*     - csig: electrical conductivity.
*     - calf: thermo-electric coefficient.
*
*     The local differential transport coefficients, defined on cells,
*     are obtained by a calls to b2tqna. To compute from these the
*     discrete transport coefficients involves multiplication by
*     geometric factors (vol, vol/hx**2, vol/hy**2), and, usually,
*     interpolation to cell faces.
*
*     We treat the anomalous electrical conductivity separately, and
*     do not add it to the classical conductivity. See treatment in
*     b2tfch.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, is, mode
      character(3) chns                                               !srv 17.06.08
      integer, save :: flux_limit_style = 2                           !srv 07.12.15
      real (kind=R8) ::
     &  csig_stoch(nFc), fchstoch(nFc),                               !lk 31.01.08
     &  chce_stoch(nFc)                                               !srv 20.09.17
      real (kind=R8) ::
     *  wrk(nFc), wrk1(nFc,2), wrkv(nVx), wrkc(nCv),
     1  hvi0(nCv,0:1), hve0(nCv,0:1),
     2  flc(nFc,0:1,0:ns-1),
     3  flv(nFc,0:1,0:ns-1),
     4  cdpa0(nFc,0:ns-1), cvsa0(nFc,0:ns-1), chcib(nFc,0:ns-1)
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external xertst, sfill, smin, smax, b2scopy
      external ipgetr, ipgeti
      external b2tqna, b2txcx, b2txcy, b2txcv, b2txsx, b2txsy
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trno')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall_b2trno.eq.0) then
      end if
*   ..extensive tests on first few calls
      if (ncall_b2trno.lt.3) then
*    ..test bb
       call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test sign of pbs
!       call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                    !srv 11.01.13
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8,
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8,
     &   'faulty argument range fcQgam(,0)')
*    ..test sign of na, te, ti
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
      endif

      mode = switch%b2tral_mode

* ..compute transport coefficients
*   ..compute differential anomalous transport coefficients
      call b2tqna (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &  switch, geo, mpg, pl, dv, rt, st_ext, co%dna0, co%dpa0, co%vla0,
     &  co%vma0, co%vsa0, co%hci0, co%hcn0, co%hcib, co%hce0, co%sig0,
     &  co%alf0, hvi0, hve0, co%dkt0, co%dzt0, co%dna_ExB, co%hce_ExB,
     &  co%hci_ExB)
*   ..incorporate geometry
      do is = 0, ns-1
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dna0(1:nCv,is), co%cdna(1:nFc,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dna0(1:nCv,is), co%cdna(1:nFc,1,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1:nCv,is), co%cdpa(1:nFc,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1:nCv,is), co%cdpa(1:nFc,1,is))
        call b2txcc (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1:nCv,is), cdpa0(1:nFc,is))
        call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   co%vla0(1:nCv,0,is), co%cvla(1:nFc,0,is))
        call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   co%vla0(1:nCv,1,is), co%cvla(1:nFc,1,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%vsa0(1:nCv,is), co%cvsa(1:nFc,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%vsa0(1:nCv,is), co%cvsa(1:nFc,1,is))
        call b2txcc (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%vsa0(1:nCv,is), cvsa0(1:nFc,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1:nCv,is), co%chcb(1:nFc,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1:nCv,is), co%chcb(1:nFc,1,is))
        call b2txcc (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1:nCv,is), chcib(1:nFc,is))
      enddo
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     & co%dna_ExB, co%cdna_ExB(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     & co%dna_ExB, co%cdna_ExB(1,1))

*   ..incorporate geometry and hz                                     !srv 09.01.01
      wrk = geo%fcS*geo%fcHz
      do iFc = 1, nFc
        wrk1(iFc,1) = geo%fcVol(iFc,1)*geo%cvHz(mpg%fcCv(iFc,1))
        wrk1(iFc,2) = geo%fcVol(iFc,2)*geo%cvHz(mpg%fcCv(iFc,2))
      enddo
      do is = 0, ns-1
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dna0(1:nCv,is), co%cdnahz(1:nFc,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dna0(1:nCv,is), co%cdnahz(1:nFc,1,is))
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dpa0(1:nCv,is), co%cdpahz(1:nFc,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dpa0(1:nCv,is), co%cdpahz(1:nFc,1,is))
       call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, wrk,
     &   co%vla0(1:nCv,0,is), co%cvlahz(1:nFc,0,is))
       call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, wrk,
     &   co%vla0(1:nCv,1,is), co%cvlahz(1:nFc,1,is))
       call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, wrk,
     &   co%vma0(1:nCv,0,is), co%cvmahz(1:nFc,0,is))
       call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, wrk,
     &   co%vma0(1:nCv,1,is), co%cvmahz(1:nFc,1,is))
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,              !srv 01.07.05 {
     &   co%vsa0(1:nCv,is), co%cvsahz(1:nFc,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%vsa0(1:nCv,is), co%cvsahz(1:nFc,1,is))                    !srv 01.07.05 }
      enddo                                                           !srv 09.01.01


*   ..Correction for advanced fluid neutral model                     !mb 09.10.17
c   poloidal pressure diffusion coefficient is only from diamagnetic contribution --> projection
      if (switch%transport_afn.eq.1) then 
       do is = 0, ns-1
         if (is_neutral(is).and.zn(is).eq.1.0_R8) then
           wrkc =  co%dpa0(:,is)*(geo%cvbb(:,2)/geo%cvbb(:,3))**2
           call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &                 wrkc, co%cdpa(1:nFc,0,is))   
           call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &                 wrkc, co%cdpahz(1:nFc,0,is))
         endif
       enddo
      endif


*   ..save the anomalous transport coefficients                       !srv 04.10.99
      co%cvsa_an = co%cvsa                                            !srv 04.10.99
      call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS,
     &   hvi0(1:nCv,0), co%chvi(1:nFc,0))
      call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS,
     &   hvi0(1:nCv,1), co%chvi(1:nFc,1))
      call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS,
     &   hve0(1:nCv,0), co%chve(1:nFc,0))
      call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS,
     &   hve0(1:nCv,1), co%chve(1:nFc,1))

      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hce0, co%chce(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hce0, co%chce(1:nFc,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%sig0, co%csig(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%sig0, co%csig(1:nFc,1))
      call b2txcx (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS,
     &  co%alf0, co%calf(1:nFc,0))
      call b2txcy (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS,
     &  co%alf0, co%calf(1:nFc,1))
      wrkc = co%dkt0*(1.0_R8 + switch%b2tfhi_fconkt)
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  wrkc, co%cdkt(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%dkt0, co%cdkt(1:nFc,1))
      wrkc = co%dzt0*(1.0_R8 + switch%b2tfhi_fconzt)
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  wrkc, co%cdzt(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%dzt0, co%cdzt(1:nFc,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hce_ExB, co%chce_ExB(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hce_ExB, co%chce_ExB(1:nFc,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hci_ExB, co%chci_ExB(1:nFc,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hci_ExB, co%chci_ExB(1:nFc,1))

*   ..set csigx = 0 and calfx = 0
      if (switch%csig_an_style.eq.1) then
        call sfill (nFc, 0.0e0_R8, co%csig(1,0), 1)
        call sfill (nFc, 0.0e0_R8, co%calf(1,0), 1)
      end if
*   ..set strange velocities chve, chvi
c      call sfill (nFc*2, 0.0_R8, co%chve, 1)
c      call sfill (nFc*2, 0.0_R8, co%chvi, 1)
*   ..add stochastic part of the electron heat conductivity
      chce_stoch = 0.0_R8                                                 !srv 20.09.17
      if(switch%fch_stochastic.ne.0.0e0_R8)then                            !lk 31.01.08 { !srv 17.12.13
      call xerrab('b2trno -- option to be adapted for WG')
#ifdef WG_TODO
        call b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne, csig_stoch,    !srv 07.07.14
     &      fchstoch)
        do iy = -1, ny                                                    !srv 20.09.17
          do ix = -1, nx
           if(bottomiy(ix,iy).ne.-2) then
            chce_stoch(ix,iy)=switch%b2trno_alpha_stoch*1.6e0_R8*
     &      csig_stoch(ix,iy)/
     &      qe/qe*(te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
           endif                                                          !srv 20.09.17
	    enddo
        enddo
        call b2saxpy (n2, 1.0_R8, chce_stoch, 1, chce(-1,-1,1,1), 1)      !srv 20.09.17 }
#endif
      endif                                                               !lk 31.01.08 } !srv 17.12.13
*   ..scale parallel coefficients if desired
      if (switch%pol_anom_scale.ne.1.0_R8) then
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            co%cdna(:,0,is) = co%cdna(:,0,is) * switch%pol_anom_scale
            co%cdpa(:,0,is) = co%cdpa(:,0,is) * switch%pol_anom_scale
            co%cvla(:,0,is) = co%cvla(:,0,is) * switch%pol_anom_scale
            co%cvsa(:,0,is) = co%cvsa(:,0,is) * switch%pol_anom_scale
            co%cdnahz(:,0,is) = co%cdnahz(:,0,is) *switch%pol_anom_scale      !wdk 16.04.14
            co%cdpahz(:,0,is) = co%cdpahz(:,0,is) *switch%pol_anom_scale
            co%cvlahz(:,0,is) = co%cvlahz(:,0,is) *switch%pol_anom_scale
            co%cvsahz(:,0,is) = co%cvsahz(:,0,is) *switch%pol_anom_scale
            co%chcb(:,0,is) = co%chcb(:,0,is) * switch%pol_anom_scale
          end if
        end do
        co%chce(:,0) = co%chce(:,0) * switch%pol_anom_scale
        co%calf(:,0) = co%calf(:,0) * switch%pol_anom_scale
        co%csig(:,0) = co%csig(:,0) * switch%pol_anom_scale
      end if
*   ..save the anomalous transport coefficients
      co%csig_an = co%csig
      co%calf_an = co%calf

      if(switch%flux_limit_to_dpa .ne. 0) then                            !srv 15.06.08 {
*   .. apply flux limit to neutrals contribution to cdpa and cdpahz
        call b2tlc0 (nCv, nFc, nVx, ns, switch, geo, mpg, pl%na, pl%te, 
     &    pl%ti, pl%tn, rt%rza, cdpa0, co%cdpa, co%cdpahz, flc)
      endif                                                               !srv 15.06.08 }

      if(switch%flux_limit_to_vsa .ne. 0) then                            !srv 15.06.08 {
*   .. apply flux limit to neutrals contribution to cvsa and cvsahz
        call b2tlv0 (nCv, nFc, nVx, ns, switch, geo, mpg, pl%na, pl%ua,
     &    pl%te, pl%tn, rt%rza, cvsa0, co%cvsa, co%cvsahz, flv)
      endif

*   ..apply flux limit to neutrals contribution to chcb
      call b2tlh0 (nCv, nFc, nVx, ns, switch, geo, mpg, pl%na, pl%ti,
     &  pl%tn, chcib, co%chcb, co%fllim0fhi)
       co%chci = 0.0_R8
       co%chcn = 0.0_R8
      do is = 0, ns-1
c the logic statements below could be simplfied but for now we go for readability
        if (switch%tn_style.eq.0) then                             !! default combined ion-neutral energy equation
          co%chci = co%chci + co%chcb(:,:,is)
        elseif (switch%tn_style.eq.1) then
         if (.not.is_neutral(is)) then                             !! pure ion energy equation, only ion contributions
          co%chci = co%chci + co%chcb(:,:,is)
         endif
        else                                                       !! separate ion-neutral energy equation
          if (.not.is_neutral(is).or.zn(is).ne.1) then
            co%chci = co%chci + co%chcb(:,:,is)
          endif
        endif
        if (is_neutral(is).and.zn(is).eq.1) then                   !! hcn for hydrogenic atoms
         co%chcn = co%chcn + co%chcb(:,:,is)
        endif
      enddo 

*   ..test sign of chce, chci
      call b2xvsg (nFc, co%chce(:,0)*geo%fcQalf(:,0), 1, 'chce', '.ge.')
      call b2xvsg (nFc, co%chce(:,1)*geo%fcQalf(:,1), 1, 'chce', '.ge.')
      call b2xvsg (nFc, co%chci(:,0)*geo%fcQalf(:,0), 1, 'chci', '.ge.')
      call b2xvsg (nFc, co%chci(:,1)*geo%fcQalf(:,1), 1, 'chci', '.ge.')
      call b2xvsg (nFc, co%chcn(:,0)*geo%fcQalf(:,0), 1, 'chcn', '.ge.')
      call b2xvsg (nFc, co%chcn(:,1)*geo%fcQalf(:,1), 1, 'chcn', '.ge.')
      call b2xvsg (nFc, co%cdkt(:,0)*geo%fcQalf(:,0), 1, 'cdkt', '.ge.')
      call b2xvsg (nFc, co%cdkt(:,1)*geo%fcQalf(:,1), 1, 'cdkt', '.ge.')

      if(switch%b2npco_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then                  !srv 19.05.10 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if (is_neutral(is) ) then
          call my_out_us(70,nFc,1,co%cdpa(1,0,is),
     &                          'b2trno_cdpa_th'//chns)
          call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2trno_cdpa_r'//chns)
        else
          call my_out_us(70,nFc,1,co%cdna(1,0,is),
     &                          'b2trno_cdna_th'//chns)
          call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2trno_cdna_r'//chns)
        endif
!        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2trno_cvlax'//chns)
!        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2trno_cvlay'//chns)
       enddo
      endif
      if(switch%b2npmo_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'b2trno_cvsa_th'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'b2trno_cvsa_r'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'b2trno_cvsahz_th'//
     &              chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'b2trno_cvsahz_r'//
     &              chns)
       enddo
      endif
      if(switch%b2npht_iout.ne.0 .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,co%chce(1,0),'b2trno_chce_th')
        call my_out_us(70,nFc,1,co%chce(1,1),'b2trno_chce_r')
        call my_out_us(70,nFc,1,co%chci(1,0),'b2trno_chci_th')
        call my_out_us(70,nFc,1,co%chci(1,1),'b2trno_chci_r')
        call my_out_us(70,nFc,1,chce_stoch,'b2trno_chcey_stoch')
      endif
      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0
     &                  .or. switch%iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,co%csig(1,0),'b2trno_csig_th')
        call my_out_us(70,nFc,1,co%csig(1,1),'b2trno_csig_r')
        call my_out_us(70,nFc,1,csig_stoch,'b2trno_csigy_stoch')
      endif                                                          !srv 19.05.10 }
      if(switch%b2trno_iout.ne.0) then                                             !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2trno_cdna_th'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2trno_cdna_r'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'b2trno_cdpa_th'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2trno_cdpa_r'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2trno_cvla_th'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2trno_cvla_r'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'b2trno_cvsa_th'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'b2trno_cvsa_r'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,0,is),
     &                           'b2trno_dnahz_th'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,1,is),
     &                           'b2trno_dnahz_r'//chns)
       call my_out_us(70,nFc,1,co%cdpahz(1,0,is),
     &                           'b2trno_cdpahz_th'//chns)
       call my_out_us(70,nFc,1,co%cdpahz(1,1,is),
     &                           'b2trno_cdpahz_r'//chns)
       call my_out_us(70,nFc,1,co%cvlahz(1,0,is),
     &                           'b2trno_cvlahz_th'//chns)
       call my_out_us(70,nFc,1,co%cvlahz(1,1,is),
     &                           'b2trno_cvlahz_r'//chns)
       call my_out_us(70,nFc,1,co%cvsahz(1,0,is),
     &                           'b2trno_cvsahz_th'//chns)
       call my_out_us(70,nFc,1,co%cvsahz(1,1,is),
     &                           'b2trno_cvsahz_r'//chns)
       call my_out_us(70,nFc,1,co%cvmahz(1,0,is),
     &                           'b2trno_cvmahz_th'//chns)
       call my_out_us(70,nFc,1,co%cvmahz(1,1,is),
     &                           'b2trno_cvmahz_r'//chns)
        if (flux_limit_style.eq.2) then                                 !srv 07.12.15 {
         call my_out_us(70,nFc,1,co%chcb(1,0,is),
     &                           'b2trno_chci_th_a'//chns)
         call my_out_us(70,nFc,1,co%chcb(1,1,is),
     &                           'b2trno_chci_r_a'//chns)
         call my_out_us(70,nFc,1,co%fllim0fhi(1,0,is),
     &                'b2trno_fllim0fhi_th'//chns)
         call my_out_us(70,nFc,1,co%fllim0fhi(1,1,is),
     &                'b2trno_fllim0fhi_r'//chns)
        endif                                                          !srv 07,12.15 }
        call my_out_us(70,nFc,1,flc(1,0,is),'b2trno_flc_th'//chns)
        call my_out_us(70,nFc,1,flc(1,1,is),'b2trno_flc_r'//chns)
        call my_out_us(70,nFc,1,flv(1,0,is),'b2trno_flv_th'//chns)
        call my_out_us(70,nFc,1,flv(1,1,is),'b2trno_flv_r'//chns)
       enddo
       call my_out_us(70,nFc,1,co%chce(1,0),'b2trno_chce_th')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2trno_chce_r')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2trno_chci_th')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2trno_chci_r')
       call my_out_us(70,nFc,1,co%chve(1,0),'b2trno_chve_th')
       call my_out_us(70,nFc,1,co%chve(1,1),'b2trno_chve_r')
       call my_out_us(70,nFc,1,co%chvi(1,0),'b2trno_chvi_th')
       call my_out_us(70,nFc,1,co%chvi(1,1),'b2trno_chvi_r')
       call my_out_us(70,nFc,1,co%csig(1,0),'b2trno_csig_th')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2trno_csig_r')
       call my_out_us(70,nFc,1,co%calf_an(1,0),'b2trno_calf_an_th')
       call my_out_us(70,nFc,1,co%calf_an(1,1),'b2trno_calf_an_r')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2trno_csig_an_th')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2trno_csig_an_r')
       call my_out_us(70,nFc,1,co%cdkt(1,0),'b2trno_cdkt_th')
       call my_out_us(70,nFc,1,co%cdkt(1,1),'b2trno_cdkt_r')
      endif                                                          !srv 17.06.08 }

* ..return
      ncall_b2trno = ncall_b2trno+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Consider to apply a flux limit to the anomalous particle transport.
*!!   Still must consider how to deal with a non-orthogonal grid.

*-----------------------------------------------------------------------
*.end b2trno

      end subroutine b2trno

