*-----------------------------------------------------------------------
*.specification

      subroutine b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,!srv 23.11.10
     1 switch, geo, mpg, pl, dv, rt, st_ext, co)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
!WG_TBD      use b2mod_neoclassical
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_b2cmpt
     &  , only : cfvla, cflim, cfvsa, cfalf, cfdpa, cfsig, cfdna, 
     &           cfhce, cfhci
      use b2mod_transport_namelist
     & , only : parm_hce, parm_hci, parm_dna
      use b2mod_input_profile
     & , only : ndata, no_pflux, tr_ip_new_files, region_flags, tdata,
     &          tr_ip_elm_count, no_div, addspec, alloc_input_profile,
     &          nss, nrr, nkind_coeff, nkind_data
      use b2mod_ad
     & , only : my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          b2news_solving, ncall_b2tral,
     &          ncall_b2tlmv, ncall_b2tqin, ncall_b2tqce, ncall_b2tlnl,
     &          ncall_b2tqca, ncall_b2tlh0, ncall_b2tqna, ncall_b2trno, 
     &          ncall_transp_keps, nsdmax, cvregmax, ncall_b2trcl
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1), 
     &  ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..input/output arguments
      type (B2Coeff), intent (inout) :: co
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRAL computes the transport coefficients, including geometric
*     factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete transport
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
*     - csigin: ion-neutral conductivity.
*     - calf: thermo-electric coefficient.
*
*     The routine operates as a simple interface to routines b2trcl
*     and b2trno, for the classical and the anomalous contributions,
*     respectively.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, is, is0
      character(3) chns, chns0                                          !srv 17.06.08
      real (kind=R8) ::
     1  chvi0(nFc,0:1), csig0(nFc,0:1), calf0(nFc,0:1), chci0(nFc,0:1),
     3  cvsa0(nFc,0:1,0:ns-1), cvsahz0(nFc,0:1,0:ns-1),                 !srv 09.01.01
     4  chce0(nFc,0:1), chve0(nFc,0:1), cdkt0(nFc,0:1),
     1	ft,                                                             !srv 16.10.17
     *  wrk(nFc), wrkmin
      character*72 message
*   ..procedures
      external subini, subend, xertst, ipgetr, b2saxpy
      external b2trcl, b2trno
*   ..initialisation
      intrinsic huge
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tral')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall_b2tral.eq.0) then
        b2news_solving(1:4)=switch%no_solve.eq.0
        do is=1,4
          b2news_solving(is)=b2news_solving(is).or.
     &    (switch%no_solve.lt.0.and.
     &     mod(abs(switch%no_solve),2**is).ge.2**(is-1))
        enddo
        call alloc_input_profile(nCV,ns)
      end if

* ..compute transport coefficients
*   ..compute anomalous terms
      call b2trno (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &  switch, geo, mpg, pl, dv, rt, st_ext, co)
*   ..compute classical terms
      call b2trcl (nCv, nFc, nVx, ns, nscx, iscx, ismain,         !srv 23.11.10
     1  switch, geo, mpg, pl, dv, rt, co, st_ext, cvsahz0,        !srv 09.01.01 03.06.03
     3  cvsa0, chce0, chve0, chci0, chvi0, csig0, calf0, cdkt0,   !xpb
     4  co%csigin, co%cthe, co%cthi, co%vsaf_cl, co%cvsa_cl, 
     5  co%cvsahz_cl, co%fllimvisc, co%csig_cl, co%calf_cl)
*   ..add to previous contributions=
      ft = 4.0e0_R8/3.0e0_R8                                      !srv 16.10.17 {
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        co%cvsa(:,0,is)   = ft*co%cvsa_cl(:,0,is)   + co%cvsa(:,0,is)
        co%cvsa(:,1,is)   =    co%cvsa_cl(:,1,is)   + co%cvsa(:,1,is)
        co%cvsahz(:,0,is) = ft*co%cvsahz_cl(:,0,is) + co%cvsahz(:,0,is)
        co%cvsahz(:,1,is) =    co%cvsahz_cl(:,1,is) + co%cvsahz(:,1,is)
       else
        co%cvsa(:,:,is)   =    co%cvsa_cl(:,:,is) + co%cvsa(:,:,is)
        co%cvsahz(:,:,is) =  co%cvsahz_cl(:,:,is) + co%cvsahz(:,:,is)
       end if
      end do                                                      !srv 16.10.17 }
      if (switch%fch_anomalous.ne.1.0_R8) then
        co%csig = co%csig * switch%fch_anomalous
        co%calf = co%calf * switch%fch_anomalous
        co%csig_an = co%csig_an * switch%fch_anomalous
        co%calf_an = co%calf_an * switch%fch_anomalous
      end if
      call b2saxpy (nFc*2, 1.0_R8, chce0, 1, co%chce, 1)
      call b2saxpy (nFc*2, 1.0_R8, chve0, 1, co%chve, 1)
      call b2saxpy (nFc*2, 1.0_R8, chci0, 1, co%chci, 1)
      call b2saxpy (nFc*2, 1.0_R8, chvi0, 1, co%chvi, 1)
      call b2saxpy (nFc*2, 1.0_R8, csig0, 1, co%csig, 1)
      call b2saxpy (nFc*2, 1.0_R8, calf0, 1, co%calf, 1)
      call b2saxpy (nFc*2, 1.0_R8, cdkt0, 1, co%cdkt, 1)
*   ..compute terms from diamagnetic parts                        !srv 13.10.06 {
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        !wdk we use the newly-defined mode = 3 here, for consistency
        !    with the computation of PSch-flows in b2tfnb and b2tfhe/i
        !wdk note: factor 1/(qe*rza) might require more careful 
        !wdk       treatment (check with other drifts in b2tfnb)
        !wdk also: need to use BzB defined at faces, simular to fcBb?
        call b2txcx (nCv, nFc, 3, geo, mpg, geo%fcVol, geo%fcS,
     &   geo%cvBzb/(qe*rt%rza(:,is)), co%cddi(1,0,is))
        call b2txcy (nCv, nFc, 3, geo, mpg, geo%fcVol, geo%fcS,
     &   geo%cvBzb/(qe*rt%rza(:,is)), co%cddi(1,1,is))
       else
        co%cddi(:,:,is) = 0.0_R8
       endif 
      enddo                                                           !srv 13.10.06 }
*  ..check that we have some particle diffusion
      if (ncall_b2tral.le.3) then                                            !srv 12.03.12 {
        do is = 0, ns-1
cwdk combined check for poloidal or radial diffusion
          wrk(:) = co%cdna(:,0,is)*geo%fcQalf(:,0) + 
     .             co%cdpa(:,0,is)*geo%fcQalf(:,0) + 
     .             co%cdna(:,1,is)*geo%fcQalf(:,1) + 
     .             co%cdpa(:,1,is)*geo%fcQalf(:,1)
          wrkmin = huge(1.0_R8)
          do iFc = 1, nFc
            wrkmin = min(wrkmin, wrk(iFc))
          enddo
          write(message,'(a, i3)') 
     .     'No particle diffusion for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. 
     .     switch%pol_anom_scale.eq.0.0_R8 .or. .not.b2news_solving(2),
     .     message)
*  ..check that we have some anomalous ion heat transport
          wrkmin = minval(co%hcib(:,is))
          write(message,'(a, i3)')
     .     'No ion heat diffusivity for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. .not.b2news_solving(4),
     .      message)
        enddo
*  ..check that we have some anomalous perpendicular currents
!WG_TODO
!WG_TODO        wrk(:) = co%calf_an(:,1)*geo%fcQalf(:,1) +
!WG_TODO     .           co%csig_an(:,1)*geo%fcQalf(:,1)
!WG_TODO        wrkmin = huge(1.0_R8)
!WG_TODO        do iFc = 1, nFc
!WG_TODO              wrkmin = min(wrkmin, wrk(iFc))
!WG_TODO        enddo
!WG_TODO        call xertst (wrkmin.gt.0.0_R8 .or. pot_eq.ne.1 .or. 
!WG_TODO     .   fch_anomalous.eq.0.0_R8 .or. .not.b2news_solving(3) .or. 
!WG_TODO     .   no_current.eq.1, 'No anomalous perpendicular currents')
      endif                                                           !srv 12.03.12 }

      if(switch%b2tral_iout.ne.0) then                                              !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'cdnax'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'cdnay'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'cdpax'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'cdpay'//chns)
        call my_out_us(70,nFc,1,co%cddi(1,0,is),'cddix'//chns)
        call my_out_us(70,nFc,1,co%cddi(1,1,is),'cddiy'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,0,is),'cvlax'//chns)
        call my_out_us(70,nFc,1,co%cvla(1,1,is),'cvlay'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,0,is),'cvsax'//chns)
        call my_out_us(70,nFc,1,co%cvsa(1,1,is),'cvsay'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,0,is),'cdnahzx'//chns)
        call my_out_us(70,nFc,1,co%cdnahz(1,1,is),'cdnahzy'//chns)
        call my_out_us(70,nFc,1,co%cdpahz(1,0,is),'cdpahzx'//chns)
        call my_out_us(70,nFc,1,co%cdpahz(1,1,is),'cdpahzy'//chns)
        call my_out_us(70,nFc,1,co%cvlahz(1,0,is),'cvlahzx'//chns)
        call my_out_us(70,nFc,1,co%cvlahz(1,1,is),'cvlahzy'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'cvsahzx'//chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'cvsahzy'//chns)
        call my_out_us(70,nFc,1,co%cvmahz(1,0,is),'cvmahzx'//chns)
        call my_out_us(70,nFc,1,co%cvmahz(1,1,is),'cvmahzy'//chns)
        call my_out_us(70,nFc,1,co%cvsa_cl(1,0,is),'cvsa_clx'//chns)
        call my_out_us(70,nFc,1,co%cvsa_cl(1,1,is),'cvsa_cly'//chns)
        call my_out_us(70,nFc,1,co%cvsahz_cl(1,0,is),'cvsahz_clx'//chns)
        call my_out_us(70,nFc,1,co%cvsahz_cl(1,1,is),'cvsahz_cly'//chns)
        call my_out_us(70,nFc,1,co%cthe(1,is),'cthe'//chns)
        call my_out_us(70,nFc,1,co%cthi(1,is),'cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then 
          call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),'csiginx'//
     &        chns//'_'//chns0)
          call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),'csiginy'//
     &        chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out_us(70,nFc,1,co%chce(1,0),'chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'chcey')
       call my_out_us(70,nFc,1,co%chci(1,0),'chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'chciy')
       call my_out_us(70,nFc,1,co%chve(1,0),'chvex')
       call my_out_us(70,nFc,1,co%chve(1,1),'chvey')
       call my_out_us(70,nFc,1,co%chvi(1,0),'chvix')
       call my_out_us(70,nFc,1,co%chvi(1,1),'chviy')
       call my_out_us(70,nFc,1,co%csig(1,0),'csigx')
       call my_out_us(70,nFc,1,co%csig(1,1),'csigy')
       call my_out_us(70,nFc,1,co%calf(1,0),'calfx')
       call my_out_us(70,nFc,1,co%calf(1,1),'calfy')
       call my_out_us(70,nFc,1,co%cdkt(1,0),'cdktx')
       call my_out_us(70,nFc,1,co%cdkt(1,1),'cdkty')
      endif                                                          !srv 17.06.08 }

      if(switch%b2npc9_iout.ne.0 .or. switch%b2npco_iout.ne.0) then                !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,co%cdna(1,0,is),'b2tral_cdnax'//chns)
        call my_out_us(70,nFc,1,co%cdna(1,1,is),'b2tral_cdnay'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,0,is),'b2tral_cdpax'//chns)
        call my_out_us(70,nFc,1,co%cdpa(1,1,is),'b2tral_cdpay'//chns)
        if(.not.is_neutral(is)) then
         call my_out_us(70,nFc,1,co%cddi(1,0,is),'b2tral_cddix'//chns)
         call my_out_us(70,nFc,1,co%cddi(1,1,is),'b2tral_cddiy'//chns)
        endif
!        call my_out_us(70,nFc,1,co%cvla(1,0,is),'b2tral_cvlax'//chns)
!        call my_out_us(70,nFc,1,co%cvla(1,1,is),'b2tral_cvlay'//chns)
       enddo
      endif
        
      if(switch%b2npmo_iout.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,co%cvsahz(1,0,is),'b2tral_cvsahzx'//
     &              chns)
        call my_out_us(70,nFc,1,co%cvsahz(1,1,is),'b2tral_cvsahzy'//
     &              chns)
!        call my_out_us(70,nFc,1,co%cvsa_cl(1,0,is),'b2tral_cvsa_clx'//
!     &              chns)
!        call my_out_us(70,nFc,1,co%cvsa_cl(1,1,is),'b2tral_cvsa_cly'//
!     &              chns)
!        call my_out_us(70,nFc,1,co%cvsahz_cl(1,0,is),'b2tral_cvsahz_clx'//
!     &             chns)
!        call my_out_us(70,nFc,1,co%cvsahz_cl(1,1,is),'b2tral_cvsahz_cly'//
!     &             chns)
        call my_out_us(70,nFc,1,co%cthe(1,is),'b2tral_cthe'//chns)
        call my_out_us(70,nFc,1,co%cthi(1,is),'b2tral_cthi'//chns)
       enddo
      endif
      if(switch%b2nph9_iout.ne.0 .or. switch%b2npht_iout.ne.0) then
       call my_out_us(70,nFc,1,co%chce(1,0),'b2tral_chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2tral_chcey')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2tral_chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2tral_chciy')
!       call my_out_us(70,nFc,1,co%chve(1,0),'b2tral_chvex')
!       call my_out_us(70,nFc,1,co%chve(1,1),'b2tral_chvey')
!       call my_out_us(70,nFc,1,co%chvi(1,0),'b2tral_chvix')
!       call my_out_us(70,nFc,1,co%chvi(1,1),'b2tral_chviy')
      endif
      if(switch%b2npp7_iout.ne.0 .or. switch%b2nppo_iout.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then 
          call my_out_us(70,nFc,1,co%csigin(1,0,is,is0),'b2tral_csiginx'
     &        //chns//'_'//chns0)
          call my_out_us(70,nFc,1,co%csigin(1,1,is,is0),'b2tral_csiginy'
     &        //chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out_us(70,nFc,1,co%csig(1,0),'b2tral_csigx')
       call my_out_us(70,nFc,1,co%csig(1,1),'b2tral_csigy')
      endif                                                          !srv 11.09.09 }

* ..return
      ncall_b2tral = ncall_b2tral+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tral

      end subroutine b2tral
