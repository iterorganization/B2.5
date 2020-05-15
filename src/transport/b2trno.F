*-----------------------------------------------------------------------
*.specification

      subroutine b2trno (nCv, nFc, nVx, ns, nscx, nscxmax, iscx,
     &  geo, mpg, pl, dv, rt, st_ext, co)
      use b2mod_types
!      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
!      use b2mod_neoclassical
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1)
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
      integer ncall, iCv, iFc, is
      integer mode
      integer, save :: iout = 0, iout_b2wdat = 0,                     !srv 17.06.08
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0,
     &                 tn_style = 0
      character(3) chns                                               !srv 17.06.08
      integer, save :: flux_limit_to_dpa = 1, csig_an_style = 1       !srv 15.06.08
      integer, save :: flux_limit_style = 2                           !srv 07.12.15
      integer, save :: user_transport = 0                             !mb  11.10.17  
      real (kind=R8), save :: fch_stochastic = 0.0_R8,
     &                        alpha_stoch = 1.0e0_R8                  !srv 17.12.13
      real (kind=R8) ::
     &  csig_stoch(nFc), fchstoch(nFc),                               !lk 31.01.08
     &  chce_stoch(nFc)                                               !srv 20.09.17
      real (kind=R8) ::
     *  wrk(nFc), wrk1(nFc,2), wrkv(nVx), wrk_afn(nCv),
     1  hvi0(nCv,0:1), hve0(nCv,0:1),
     2  flc(nFc,0:1,0:ns-1),          !srv 09.07.03 30.11.15
     3  flchz(nFc,0:1,0:ns-1),                            !srv 09.07.03 30.11.15
     4	pol_anom_scale,
     5  cdpa0(nFc,0:ns-1), cdpa0hz(nFc,0:ns-1), chcib(nFc,0:ns-1)
      save ncall, pol_anom_scale, mode
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external subini, subend, xertst, sfill, smin, smax, b2scopy
      external ipgetr, ipgeti
      external b2tqna, b2txcx, b2txcy, b2txcv, b2txsx, b2txsy
*   ..initialisation
      data ncall/0/, pol_anom_scale/1.0_R8/, mode/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trno')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tqna_user_transport', user_transport)        !mb  11.10.17
        call ipgetr ('b2trno_pol_anom_scale', pol_anom_scale)
        call ipgeti ('b2tral_mode', mode)
        call ipgeti ('b2trno_flux_limit_to_dpa', flux_limit_to_dpa)
        call ipgeti ('b2trno_csig_an_style', csig_an_style)          !srv 15.06.08
        call ipgeti ('b2trno_iout', iout)                            !srv 17.06.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 19.05.10 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 19.05.10 }
        call ipgetr ('b2tfhe_stochastic', fch_stochastic)            !lk 31.01.08
        call ipgetr ('b2trno_con_e_stochastic', alpha_stoch)         !lk 31.01.08
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
        call ipgeti ('b2mn_tn_style', tn_style)                  
      end if
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
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

* ..compute transport coefficients
*   ..compute differential anomalous transport coefficients
      call b2tqna (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, geo, mpg,
     &  pl, dv, rt, st_ext, co%dna0, co%dpa0, co%vla0, co%vma0, co%vsa0,
     &  co%hci0, co%hcn0, co%hcib, co%hce0, co%sig0, co%alf0, hvi0, 
     &  hve0)
*   ..incorporate geometry
      do is = 0, ns-1
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dna0(1,is), co%cdna(1,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dna0(1,is), co%cdna(1,1,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1,is), co%cdpa(1,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1,is), co%cdpa(1,1,is))
        call b2txcc (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%dpa0(1,is), cdpa0(1,is))
        call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   co%vla0(1,0,is), co%cvla(1,0,is))
        call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   co%vla0(1,1,is), co%cvla(1,1,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%vsa0(1,is), co%cvsa(1,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%vsa0(1,is), co%cvsa(1,1,is))
        call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1,is), co%chcb(1,0,is))
        call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1,is), co%chcb(1,1,is))
        call b2txcc (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   co%hcib(1,is), chcib(1,is))
      enddo

*   ..incorporate geometry and hz                                     !srv 09.01.01
      wrk = geo%fcS*geo%fcHz
      do iFc = 1, nFc
        wrk1(iFc,1) = geo%fcVol(iFc,1)*geo%cvHz(mpg%fcCv(iFc,1))
        wrk1(iFc,2) = geo%fcVol(iFc,2)*geo%cvHz(mpg%fcCv(iFc,2))
      enddo
      do is = 0, ns-1
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dna0(1,is), co%cdnahz(1,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dna0(1,is), co%cdnahz(1,1,is))
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dpa0(1,is), co%cdpahz(1,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dpa0(1,is), co%cdpahz(1,1,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%dpa0(1,is), cdpa0hz(1,is))
       call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, wrk, 
     &   co%vla0(1,0,is), co%cvlahz(1,0,is))
       call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, wrk, 
     &   co%vla0(1,1,is), co%cvlahz(1,1,is))
       call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, wrk, 
     &   co%vma0(1,0,is), co%cvmahz(1,0,is))
       call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, wrk, 
     &   co%vma0(1,1,is), co%cvmahz(1,1,is))
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,              !srv 01.07.05 {
     &   co%vsa0(1,is), co%cvsahz(1,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   co%vsa0(1,is), co%cvsahz(1,1,is))                            !srv 01.07.05 }
      enddo                                                           !srv 09.01.01


*   ..Correction for advanced fluid neutral model                     !mb 09.10.17
c   poloidal pressure diffusion coefficient is only from diamagnetic contribution --> projection
      if (user_transport.eq.10) then 
       do is = 0, ns-1
         if (is_neutral(is).and.zn(is).eq.1.0_R8) then
           wrk_afn =  co%dpa0(:,is)*geo%cvbb(:,2)/geo%cvbb(:,3)
           call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &                 wrk_afn, co%cdpa(1,0,is))   
           call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &                 wrk_afn, co%cdpahz(1,0,is))
         endif
       enddo
      endif


*   ..save the anomalous transport coefficients                       !srv 04.10.99
      co%cvsa_an = co%cvsa                                            !srv 04.10.99
      call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   hvi0(1,0), co%chvi(1,0))
      call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   hvi0(1,1), co%chvi(1,1))
      call b2txsx (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   hve0(1,0), co%chve(1,0))
      call b2txsy (nCv, nFc, geo, mpg, geo%fcVol, geo%fcS, 
     &   hve0(1,1), co%chve(1,1))

      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS, 
     &  co%hce0, co%chce(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%hce0,  co%chce(1,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  co%sig0, co%csig(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS, 
     &  co%sig0, co%csig(1,1))
      call b2txcx (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS, 
     &  co%alf0, co%calf(1,0))
      call b2txcy (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS, 
     &  co%alf0, co%calf(1,1))

*   .. set csigx = 0 and calfx = 0
      if (csig_an_style.eq.1) then
        call sfill (nFc, 0.0e0_R8, co%csig(1,0), 1)
        call sfill (nFc, 0.0e0_R8, co%calf(1,0), 1)
      end if
*   ..set strange velocities chve, chvi
c      call sfill (nFc*2, 0.0_R8, co%chve, 1)
c      call sfill (nFc*2, 0.0_R8, co%chvi, 1)
*   ..add stochastic part of the electron heat conductivity
      chce_stoch = 0.0_R8                                                 !srv 20.09.17
      if(fch_stochastic.ne.0.0e0_R8)then                                  !lk 31.01.08 { !srv 17.12.13
      call xerrab('b2trno -- option to be adapted for WG')
#ifdef WG_TODO
        call b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne, csig_stoch,    !srv 07.07.14
     &      fchstoch)
        do iy = -1, ny                                                    !srv 20.09.17
          do ix = -1, nx
           if(bottomiy(ix,iy).ne.-2) then
            chce_stoch(ix,iy)=alpha_stoch*1.6e0_R8*csig_stoch(ix,iy)/
     &      qe/qe*(te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0_R8
           endif                                                          !srv 20.09.17
	    enddo
        enddo
        call b2saxpy (n2, 1.0_R8, chce_stoch, 1, chce(-1,-1,1,1), 1)        !srv 20.09.17 }
#endif
      endif                                                               !lk 31.01.08 } !srv 17.12.13
*   ..scale parallel coefficients if desired
      if (pol_anom_scale.ne.1.0_R8) then
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            co%cdna(:,0,is) = co%cdna(:,0,is) * pol_anom_scale
            co%cdpa(:,0,is) = co%cdpa(:,0,is) * pol_anom_scale
            co%cvla(:,0,is) = co%cvla(:,0,is) * pol_anom_scale
            co%cvsa(:,0,is) = co%cvsa(:,0,is) * pol_anom_scale
            co%cdnahz(:,0,is) = co%cdnahz(:,0,is) * pol_anom_scale      !wdk 16.04.14
            co%cdpahz(:,0,is) = co%cdpahz(:,0,is) * pol_anom_scale
            co%cvlahz(:,0,is) = co%cvlahz(:,0,is) * pol_anom_scale
            co%cvsahz(:,0,is) = co%cvsahz(:,0,is) * pol_anom_scale
            co%chcb(:,0,is) = co%chcb(:,0,is) * pol_anom_scale
          end if
        end do
        co%chce(:,0) = co%chce(:,0) * pol_anom_scale
        co%calf(:,0) = co%calf(:,0) * pol_anom_scale
        co%csig(:,0) = co%csig(:,0) * pol_anom_scale
      end if
*   .. save the anomalous transport coefficients
      co%csig_an = co%csig
      co%calf_an = co%calf

      if(flux_limit_to_dpa .ne. 0) then                                   !srv 15.06.08 {
*   .. apply flux limit to neutrals contribution to cdpa
        wrk = 1.0_R8
        call b2tlc0 (nCv, nFc, nVx, ns, geo, mpg, wrk, pl%na, pl%te, 
     &    pl%ti, pl%tn, rt%rza, cdpa0, co%cdpa, flc) !srv 09.07.03 30.11.15
        call b2tlc0 (nCv, nFc, nVx, ns, geo, mpg, geo%fcHz, pl%na,
     &    pl%te, pl%ti, pl%tn, rt%rza, cdpa0hz, co%cdpahz, flchz) 
      endif                                                               !srv 15.06.08 }

*   ..apply flux limit to neutrals contribution to chcb
      call b2tlh0 (nCv, nFc, nVx, ns, geo, mpg, pl%na, pl%ti,
     &  pl%tn, chcib, co%chcb, co%fllim0fhi)
       co%chci = 0.0_R8
       co%chcn = 0.0_R8
      do is = 0, ns-1
c the logic statements below could be simplfied but for now we go for readability
        if (tn_style.eq.0) then                             !! default combined ion-neutral energy equation
          co%chci(:,:) = co%chci(:,:) + co%chcb(:,:,is)
        elseif (tn_style.eq.1.and..not.is_neutral(is)) then !! pure ion energy equation, only ion contributions
          co%chci(:,:) = co%chci(:,:) + co%chcb(:,:,is)
        else                                                !! separate ion-neutral energy equation
          if (.not.is_neutral(is)) then
            co%chci(:,:) = co%chci(:,:) + co%chcb(:,:,is)
          else
            co%chcn(:,:) = co%chcn(:,:) + co%chcb(:,:,is)
          endif
        endif
      enddo 

*    ..test sign of chce, chci
      call b2xvsg (nFc, co%chce(:,0)*geo%fcQalf(:,0), 1, 'chce', '.ge.')
      call b2xvsg (nFc, co%chce(:,1)*geo%fcQalf(:,1), 1, 'chce', '.ge.')
      call b2xvsg (nFc, co%chci(:,0)*geo%fcQalf(:,0), 1, 'chci', '.ge.')
      call b2xvsg (nFc, co%chci(:,1)*geo%fcQalf(:,1), 1, 'chci', '.ge.')
      call b2xvsg (nFc, co%chcn(:,0)*geo%fcQalf(:,0), 1, 'chcn', '.ge.')
      call b2xvsg (nFc, co%chcn(:,1)*geo%fcQalf(:,1), 1, 'chcn', '.ge.')

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0
     &                  .or. iout_b2wdat.eq.4) then                  !srv 19.05.10 {
       do is = 0, ns-1
	write (chns,'(i3.3)') is
        if (is_neutral(is) ) then
          call my_out_us(70,nFc,1,co%cdpa(1,0,is),'b2trno_cdpax'//chns)
          call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2trno_cdpay'//chns)
        else
          call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2trno_cdnax'//chns)
          call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2trno_cdnay'//chns)
        endif
!        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2trno_cvlax'//chns)
!        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2trno_cvlay'//chns)
       enddo
      endif
      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'b2trno_cvsax'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'b2trno_cvsay'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'b2trno_cvsahzx'//
     &              chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'b2trno_cvsahzy'//
     &              chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,co%chce(1,0),'b2trno_chcex')
        call my_out_us(70,nFc,1,co%chce(1,1),'b2trno_chcey')
        call my_out_us(70,nFc,1,co%chci(1,0),'b2trno_chcix')
        call my_out_us(70,nFc,1,co%chci(1,1),'b2trno_chciy')
        call my_out_us(70,nFc,1,chce_stoch,'b2trno_chcey_stoch')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,co%csig(1,0),'b2trno_csigx')
        call my_out_us(70,nFc,1,co%csig(1,1),'b2trno_csigy')
        call my_out_us(70,nFc,1,csig_stoch,'b2trno_csigy_stoch')
      endif                                                          !srv 19.05.10 }
      if(iout.ne.0) then                                             !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2trno_cdnax'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2trno_cdnay'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'b2trno_cdpax'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2trno_cdpay'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2trno_cvlax'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2trno_cvlay'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'b2trno_cvsax'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'b2trno_cvsay'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,0,is),'b2trno_dnahzx'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,1,is),'b2trno_dnahzy'//chns)
       call my_out_us(70,nFc,1,co%cdpahz(1,0,is),'b2trno_cdpahzx'//chns)
       call my_out_us(70,nFc,1,co%cdpahz(1,1,is),'b2trno_cdpahzy'//chns)
       call my_out_us(70,nFc,1,co%cvlahz(1,0,is),'b2trno_cvlahzx'//chns)
       call my_out_us(70,nFc,1,co%cvlahz(1,1,is),'b2trno_cvlahzy'//chns)
       call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'b2trno_cvsahzx'//chns)
       call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'b2trno_cvsahzy'//chns)
       call my_out_us(70,nFc,1,co%cvmahz(1,0,is),'b2trno_cvmahzx'//chns)
       call my_out_us(70,nFc,1,co%cvmahz(1,1,is),'b2trno_cvmahzy'//chns)
        if (flux_limit_style.eq.2) then                                 !srv 07.12.15 {
         call my_out_us(70,nFc,1,co%chcb(1,0,is),'b2trno_chcixa'//chns)
         call my_out_us(70,nFc,1,co%chcb(1,1,is),'b2trno_chciya'//chns)
         call my_out_us(70,nFc,1,co%fllim0fhi(1,0,is),
     &                'b2trno_fllim0fhix'//chns)
         call my_out_us(70,nFc,1,co%fllim0fhi(1,1,is),
     &                'b2trno_fllim0fhiy'//chns)
        endif                                                          !srv 07,12.15 }
        call my_out_us(70,nFc,1,flc(1,0,is),'b2trno_flcx'//chns)
        call my_out_us(70,nFc,1,flc(1,1,is),'b2trno_flcy'//chns)
        call my_out_us(70,nFc,1,flchz(1,0,is),'b2trno_flchzx'//chns)
        call my_out_us(70,nFc,1,flchz(1,1,is),'b2trno_flchzy'//chns)
       enddo
       call my_out_us(70,nFc,1,co%chce(1,0),'b2trno_chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2trno_chcey')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2trno_chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2trno_chciy')
       call my_out_us(70,nFc,1,co%chve(1,0),'b2trno_chvex')
       call my_out_us(70,nFc,1,co%chve(1,1),'b2trno_chvey')
       call my_out_us(70,nFc,1,co%chvi(1,0),'b2trno_chvix')
       call my_out_us(70,nFc,1,co%chvi(1,1),'b2trno_chviy')
       call my_out_us(70,nFc,1,co%csig(1,0),'b2trno_csigx')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2trno_csigy')
       call my_out_us(70,nFc,1,co%calf_an(1,0),'b2trno_calf_anx')
       call my_out_us(70,nFc,1,co%calf_an(1,1),'b2trno_calf_any')
       call my_out_us(70,nFc,1,co%csig_an(1,0),'b2trno_csig_anx')
       call my_out_us(70,nFc,1,co%csig_an(1,1),'b2trno_csig_any')
      endif                                                          !srv 17.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Consider to apply a flux limit to the anomalous particle transport.
*!!   Still must consider how to deal with a non-orthogonal grid.

*-----------------------------------------------------------------------
*.end b2trno

      end

