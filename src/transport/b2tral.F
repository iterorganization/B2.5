*-----------------------------------------------------------------------
*.specification

      subroutine b2tral (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,!srv 23.11.10
     1 geo, mpg, pl, dv, rt, st_ext, co)
      use b2mod_types
      use b2mod_constants
!WG_TBD      use b2mod_input_profile
      use b2mod_b2cmpa
!WG_TBD      use b2mod_neoclassical
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1), 
     &  ismain
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
      integer iCv, iFc, is, ncall, pot_eq, no_current, no_solve, is0
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0
      character(3) chns, chns0                                          !srv 17.06.08
      real (kind=R8) ::
     1  chvi0(nFc,0:1), csig0(nFc,0:1), calf0(nFc,0:1), chci0(nFc,0:1),
     3  cvsa0(nFc,0:1,0:ns-1), cvsahz0(nFc,0:1,0:ns-1),                 !srv 09.01.01
     4  chce0(nFc,0:1), chve0(nFc,0:1),
     1	pol_anom_scale, ft,                                             !srv 16.10.17
     *  wrk(nFc), wrkmin, fch_anomalous
      character*72 message
      logical solving(4)
      save ncall, pol_anom_scale, pot_eq, no_solve, solving, 
     .  fch_anomalous, no_current
*   ..procedures
      external subini, subend, xertst, ipgetr, b2saxpy
      external b2trcl, b2trno
*   ..initialisation
      data ncall/0/, pol_anom_scale/1.0_R8/, pot_eq/1/, no_solve/0/,
     1  fch_anomalous/1.0_R8/, no_current/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tral')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall.eq.0) then
        call ipgetr ('b2trno_pol_anom_scale', pol_anom_scale)
        call ipgetr ('b2tfhe_anomalous', fch_anomalous)
        if (fch_anomalous.eq.1.0_R8)
     >   call ipgetr ('b2tanml_anomalous', fch_anomalous) ! for backward compatibility
        call ipgeti ('b2tfhe_no_current', no_current)
        call ipgeti ('b2news_poteq', pot_eq)
        call xertst (pot_eq.eq.1.or.no_current.eq.1,
     &   'If b2news_poteq.ne.1, then b2tfhe_no_current must be 1!')
        call ipgeti ('b2news_no_solve', no_solve)
        call xertst (-15.le.no_solve.and.no_solve.le.1,
     &   'faulty input no_solve')
        solving(1:4)=no_solve.eq.0
        do is=1,4
          solving(is)=solving(is).or.
     &    (no_solve.lt.0.and.mod(abs(no_solve),2**is).ge.2**(is-1))
        enddo
        call ipgeti ('b2tral_iout', iout)                         !srv 17.06.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                  !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                  !srv 11.09.09 }
!WG_TBD        call alloc_input_profile(nx,ny,ns)
      end if

* ..compute transport coefficients
*   ..compute anomalous terms
      call b2trno (nCv, nFc, nVx, ns, nscx, nscxmax, iscx,
     &  geo, mpg, pl, dv, rt, st_ext, co)
*   ..compute classical terms
      call b2trcl (nCv, nFc, nVx, ns, nscx, iscx, ismain,         !srv 23.11.10
     1  geo, mpg, pl, dv, rt, co, st_ext, cvsahz0,                !srv 09.01.01 03.06.03
     3  cvsa0, chce0, chve0, chci0, chvi0, csig0, calf0,          !xpb
     4  co%csigin, co%cthe, co%cthi, co%cvsa_cl, co%cvsahz_cl,    !srv 09.01.01
     5  co%fllimvisc, co%csig_cl, co%calf_cl)
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
      if (fch_anomalous.ne.1.0_R8) then
        co%csig = co%csig * fch_anomalous
        co%calf = co%calf * fch_anomalous
        co%csig_an = co%csig_an * fch_anomalous
        co%calf_an = co%calf_an * fch_anomalous
      end if
      call b2saxpy (nFc*2, 1.0_R8, chce0, 1, co%chce, 1)
      call b2saxpy (nFc*2, 1.0_R8, chve0, 1, co%chve, 1)
      call b2saxpy (nFc*2, 1.0_R8, chci0, 1, co%chci, 1)
      call b2saxpy (nFc*2, 1.0_R8, chvi0, 1, co%chvi, 1)
      call b2saxpy (nFc*2, 1.0_R8, csig0, 1, co%csig, 1)
      call b2saxpy (nFc*2, 1.0_R8, calf0, 1, co%calf, 1)
*   ..compute terms from diamagnetic parts                        !srv 13.10.06 {
      co%cddi = 0.0_R8
#ifdef WG_DRIFTS
      do is = 0, ns-1
       if(.not.is_neutral(is)) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
           cddi(ix,iy,0,0,is) =
     &      +2.0e0_R8*(vol(leftix(ix,iy),leftiy(ix,iy))*
     &         bzb(leftix(ix,iy),leftiy(ix,iy))
     &         +vol(ix,iy)*bzb(ix,iy))/
     &       ((hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*
     &        (hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy)))/
     &        (rza(ix,iy,is)*qe)
           cddi(ix,iy,0,1,is) = cddi(ix,iy,0,0,is)
          endif
         enddo
        enddo

        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).ne.-2) then
           cddi(ix,iy,1,1,is) =
     &       2.0e0_R8*(vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &         bzb(bottomix(ix,iy),bottomiy(ix,iy))
     &         +vol(ix,iy)*bzb(ix,iy))/
     &       ((hx(bottomix(ix,iy),bottomiy(ix,iy))+hx(ix,iy))*
     &        (hy(bottomix(ix,iy),bottomiy(ix,iy))+hy(ix,iy)))/
     &        (rza(ix,iy,is)*qe)
          endif
         enddo
        enddo

       endif 
      enddo                                                           !srv 13.10.06 }
#endif
*  ..check that we have some particle diffusion
      if (ncall.le.3) then                                            !srv 12.03.12 {
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
          call xertst (wrkmin.gt.0.0_R8 .or. pol_anom_scale.eq.0.0_R8
     .     .or. .not.solving(2),message)
*  ..check that we have some anomalous ion heat transport
          wrkmin = minval(co%hcib(:,is))
          write(message,'(a, i3)')
     .     'No ion heat diffusivity for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. .not.solving(4), message)
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
!WG_TODO     .   fch_anomalous.eq.0.0_R8 .or. .not.solving(3) .or. 
!WG_TODO     .   no_current.eq.1, 'No anomalous perpendicular currents')
      endif                                                           !srv 12.03.12 }

      if(iout.ne.0) then                                              !srv 17.06.08 {
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
      endif                                                          !srv 17.06.08 }

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0) then                !srv 11.09.09 {
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
        
      if(iout_b2npmo.ne.0) then
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
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
       call my_out_us(70,nFc,1,co%chce(1,0),'b2tral_chcex')
       call my_out_us(70,nFc,1,co%chce(1,1),'b2tral_chcey')
       call my_out_us(70,nFc,1,co%chci(1,0),'b2tral_chcix')
       call my_out_us(70,nFc,1,co%chci(1,1),'b2tral_chciy')
!       call my_out_us(70,nFc,1,co%chve(1,0),'b2tral_chvex')
!       call my_out_us(70,nFc,1,co%chve(1,1),'b2tral_chvey')
!       call my_out_us(70,nFc,1,co%chvi(1,0),'b2tral_chvix')
!       call my_out_us(70,nFc,1,co%chvi(1,1),'b2tral_chviy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
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
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tral

      end subroutine b2tral
