*-----------------------------------------------------------------------
*.specification

      subroutine b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,       !srv 23.11.10
     1 bb, conn, vol, gs, hx, hy, hz, qz, qc, pbs, crx, cry, bzb,       !srv 09.01.01 01.07.09
     1 lnlam, fch_p, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,   !srv 03.06.03
     2 cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,                  !srv 13.10.06
     3 chvi, csig, csigin, calf, cthe, cthi,
     4 cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,      !srv 09.01.01 15.12.05
     4 fllim0fhi, fllimvisc, csig_cl, calf_cl)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_input_profile
      use b2mod_transport_namelist
      use b2mod_anomalous_transport
      use b2mod_time
      use b2mod_rates
     , , only : rza
      use b2mod_b2cmpa
      use b2mod_b2zhco ,
     & only: zhflcorr,zhcscorr,zh_imp_analyt,zh_analyt_mdf,             !som 21.03.21
     & zhflcorrh,zh_dist_fun 
      use b2mod_zhfrtf,
     & only: avm_u,w_out,hTdp_out,hwdp_out,
     & rTdpgamma_out,rwdpgamma_out
      use b2mod_neoclassical
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, nscxmax, iscx(0:nscxmax-1), ismain      !srv 23.11.10
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),              !srv 09.01.01
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1), bzb(-1:nx,-1:ny),    !srv 13.10.06
     &  lnlam(-1:nx,-1:ny),                                             !srv 01.07.09
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne2(-1:nx,-1:ny),
     &  fch_p(-1:nx,-1:ny,0:1), conn(-1:ny), qc(-1:nx,-1:ny),           !srv 10.06.08
     &  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)                      !srv 03.06.03
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cdna(-1:nx,-1:ny,0:1,0:ns-1), cdpa(-1:nx,-1:ny,0:1,0:ns-1),
     *  cddi(-1:nx,-1:ny,0:1,0:ns-1),                                   !srv 13.10.06
     &  cvla(-1:nx,-1:ny,0:1,0:ns-1), cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  cdnahz(-1:nx,-1:ny,0:1,0:ns-1), cdpahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvlahz(-1:nx,-1:ny,0:1,0:ns-1), cvsahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvsahz_cl(-1:nx,-1:ny,0:1,0:ns-1),                              !srv 09.01.01
     &  cvmahz(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 15.12.05
     &  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     &  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  csig(-1:nx,-1:ny,0:1), csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),   !xpb
     &  calf(-1:nx,-1:ny,0:1),
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     1  cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1),
     1  fllim0fhi(-1:nx,-1:ny,0:1,0:ns-1),
     2  fllimvisc(-1:nx,-1:ny,0:ns-1),
     3  csig_cl(-1:nx,-1:ny,0:1), calf_cl(-1:nx,-1:ny,0:1),
     &  chvemx(-1:nx,-1:ny), chvimx(-1:nx,-1:ny)
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
      integer ix, iy, is, is0, n2
      integer ncall, pot_eq, no_current, no_solve, style, elm_count
      integer, save :: nnucl, zmaxmax                                   !som 09.02.21
      integer, save :: transport_namelist = 0
      integer, save :: Zhdanov_closure = 0                              !som 25.01.21
      integer, save :: Zhdanov_test = 0                                 !som 16.06.21
      integer, save :: Zhdanov_cond = 1                                 !som 22.06.21
      integer, save :: icase = 0                                        !som 05.03.21
      integer, save :: mode = 1                                         !som 10.06.21
      real (kind=R8), save :: gamma = 1.0_R8                            !som 20.03.21
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,               !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, 
     &                 iout_b2npmo = 0, iout_b2wdat = 0
      character(3) chns, chns0                                          !srv 17.06.08
      real (kind=R8) ::
     1  chvi0(-1:nx,-1:ny,0:1), csig0(-1:nx,-1:ny,0:1),
     2  calf0(-1:nx,-1:ny,0:1), chci0(-1:nx,-1:ny,0:1),
     3  cvsa0(-1:nx,-1:ny,0:1,0:ns-1), cvsahz0(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     4  chce0(-1:nx,-1:ny,0:1), chve0(-1:nx,-1:ny,0:1),
     2  chci0_nofl(-1:nx,-1:ny,0:1), chci0_luc(-1:nx,-1:ny,0:1),        !som 16.06.21
     2  chci0_fl(-1:nx,-1:ny,0:1),                                      !som 16.06.21
     1  pol_anom_scale, ft,                                             !srv 16.10.17
     *  wrk(-1:nx,-1:ny), wrkmin, fch_anomalous
      character*72 message
      logical solving(4), catch_up
      save ncall, pol_anom_scale, pot_eq, no_solve, solving,
     .  fch_anomalous, no_current, style, elm_count
*   ..procedures
      external xertst, ipgetr, b2saxpy
      external b2trcl, b2trql, b2trno
*   ..initialisation
      data ncall/0/, pol_anom_scale/1.0_R8/, pot_eq/1/, no_solve/0/,
     1  fch_anomalous/1.0_R8/, no_current/0/, style/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tral')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tral_cvsa_style', style)
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
        call ipgeti ('b2wdat_iout', iout_b2wdat)                  !som 16.06.21
        call ipgeti ('b2npc9_iout', iout_b2npc9)                  !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npc7_iout', iout_b2npc7)                  !srv 02.04.18
        call ipgeti ('b2npc11_iout', iout_b2npc11)                !srv 02.04.18
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                  !srv 11.09.09 }
        call alloc_input_profile(nx,ny,ns)
        call ipgeti ('b2tqna_transport_namelist', transport_namelist)
        call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)   !som 25.01.21
        call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)         !som 16.06.21
        call ipgeti ('b2tral_Zhdanov_cond', Zhdanov_cond)         !som 22.06.21
        call ipgetr ('b2tral_zhflcorr', zhflcorr)
        call ipgetr ('b2tral_zhflcorrh', zhflcorrh)               !som 08.06.21
        call ipgeti ('b2tral_zhcscorr', zhcscorr)
        call ipgeti ('b2tlnl_ii', icase)                          !som 05.03.21
        call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
        call ipgetr ('b2tlc0_gamma', gamma)                       !som 20.03.21
        call xertst (0.lt.gamma, 'invalid input parameter gamma')      
        call ipgeti ('b2tral_zh_imp_analyt', zh_imp_analyt)       !som 21.03.21
        call ipgeti ('b2tral_zh_analyt_mdf', zh_analyt_mdf)       !som 21.03.21
        call ipgeti ('b2tral_zh_dist_fun', zh_dist_fun)           !som 10.06.21
        call ipgeti ('b2tral_zh_tf_toff', zh_tf_toff)             !som 26.06.21
        call ipgeti ('b2tral_mode', mode)                         !som 15.03.21
        if(transport_namelist.eq.1) then
          call alloc_transport_namelist(ns)
          call read_b2mod_transport_namelist
          call write_b2mod_transport_namelist
          if(transport_time_mod.gt.0.0_R8) then
            elm_count=int(tim/transport_time_mod)
          else
            elm_count=0
          endif
          catch_up=transport_time_mod.gt.0.0_R8
          do while(catch_up)
            if(mod(tim,transport_time_mod).gt.transport_time_switch)
     1       then
              write(*,*) 'b2tqna: catching up ...'
              call read_b2mod_transport_namelist
              call write_b2mod_transport_namelist
              catch_up=transport_time_mod.gt.0.0_R8
            else
              catch_up=.false.
            endif
          enddo
        endif
      else
        if(transport_namelist.eq.1.and.
     1   transport_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,transport_time_mod),
     1     transport_time_switch,int(tim/transport_time_mod),elm_count
          if(mod(tim,transport_time_mod).gt.transport_time_switch.or.
     1     int(tim/transport_time_mod).gt.elm_count) then
            elm_count=int(tim/transport_time_mod)
            call read_b2mod_transport_namelist
            call write_b2mod_transport_namelist
          endif
        endif
      end if

* ..compute transport coefficients
*   ..compute flux limit coefficients
      call b2trql (nx, ny, ns, gs, vol, hx, hy, qz, qc, pbs,
     &  na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx)
*   ..compute anomalous terms
      call b2trno (nx, ny, ns, nscx, nscxmax, iscx,
     &  bb, vol, gs, hx, hy, hz, qz, pbs, crx, cry, conn,         !srv 09.01.01
     &  na, ua, te, ti, po, ne, ni, ne2, cdna, cdpa, cvla, cvsa,
     &  cdnahz, cdpahz, cvlahz, cvmahz, cvsahz,                   !srv 09.01.01 15.12.05
     &  chce, chve, chci, chvi, csig, calf, fllim0fhi)
*   ..compute classical terms
      call b2trcl (nx, ny, ns, nscx, iscx, ismain,
     &  bb, conn, vol, gs, hx, hy, hz,                            !srv 23.11.10
     1  qz, qc, pbs, fch_p, na, ua, te, ti, po, ne, ni, ne2, lnlam,   !srv 01.07.09
     2  chvemx, chvimx, cvsahz0,                                  !srv 09.01.01 03.06.03
     3  cvsa0, chce0, chve0, chci0, chvi0, csig0, calf0, csigin,
     4  cthe, cthi, cvsa_cl, cvsahz_cl, fllimvisc,                !srv 09.01.01
     5  csig_cl, calf_cl)
*   ..compute Zhdanov closure coefficients
      if(Zhdanov_closure.eq.1) then                                !som 09.02.21
        if (ncall.eq.0) then
          call b2trnu (ns, nnucl, zmaxmax)
          if (nnucl /= nspecies) then
            stop 'Error during making isonuclear sequence!'
          end if
        endif
        call b2tlnl (nx, ny, icase, te, ti, ne, lnlam) 
        if(Zhdanov_test.eq.0.and.Zhdanov_cond.eq.1) then              !som 25.06.21
          call b2trzh (nx, ny, ns, zmaxmax, na, ua, ti, lnlam,hx,bb,
     & vol,gs,qz,mode,gamma,chci0_nofl,chci0_fl,chci0)

          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                  .or. iout_b2wdat.eq.4) then
            call my_out(70,nx,ny,chci0_nofl(-1,-1,0),
     & 'b2trzh_chcix')
            call my_out(70,nx,ny,chci0_fl(-1,-1,0),
     & 'b2trzh_fllim_chcix')
            call my_out(70,nx,ny,chci0(-1,-1,0),
     & 'b2trzh_luciani_fllim_chcix')
          endif
        else
          call b2trzh (nx, ny, ns, zmaxmax, na, ua, ti, lnlam,hx,bb,
     & vol,gs,qz,mode,gamma,chci0_nofl,chci0_fl,chci0_luc)

          if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                  .or. iout_b2wdat.eq.4) then
            call my_out(70,nx,ny,chci0_nofl(-1,-1,0),
     & 'b2trzh_chcix')
            call my_out(70,nx,ny,chci0_fl(-1,-1,0),
     & 'b2trzh_fllim_chcix')
            call my_out(70,nx,ny,chci0_luc(-1,-1,0),
     & 'b2trzh_luciani_fllim_chcix')
          endif
        endif

        if((iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     & .or. iout_b2wdat.eq.4) .and. zh_dist_fun.eq.1) then
          call my_out(70,nx,ny,avm_u(-1,-1),
     & 'b2trzh_avm_u')
          do is = 0, ns-1
            write (chns,'(i3.3)') is
            call my_out(70,nx,ny,w_out(-1,-1,is),
     & 'b2trzh_w'//chns)
            call my_out(70,nx,ny,hTdp_out(-1,-1,is),
     & 'b2trzh_hTdp'//chns)
            call my_out(70,nx,ny,hwdp_out(-1,-1,is),
     & 'b2trzh_hwdp'//chns)
            call my_out(70,nx,ny,rTdpgamma_out(-1,-1,is),
     & 'b2trzh_rTdpgamma'//chns)
            call my_out(70,nx,ny,rwdpgamma_out(-1,-1,is),
     & 'b2trzh_rwdpgamma'//chns)
          enddo
        endif
      endif     
*   ..add to previous contributions
      ft = 4.0e0_R8/3.0e0_R8                                      !srv 16.10.17
      if (style.eq.1) then
       do is = 0, ns-1
        if (.not.is_neutral(is)) then
          cvsa(:,:,0,is)   = ft*cvsa_cl(:,:,0,is)   + cvsa(:,:,0,is)
          cvsa(:,:,1,is)   =    cvsa_cl(:,:,1,is)   + cvsa(:,:,1,is)
          cvsahz(:,:,0,is) = ft*cvsahz_cl(:,:,0,is) + cvsahz(:,:,0,is)
          cvsahz(:,:,1,is) =    cvsahz_cl(:,:,1,is) + cvsahz(:,:,1,is)
        else
          cvsa(:,:,:,is)   =    cvsa_cl(:,:,:,is) + cvsa(:,:,:,is)
          cvsahz(:,:,:,is) =  cvsahz_cl(:,:,:,is) + cvsahz(:,:,:,is)
        end if
       end do
      else
       call b2saxpy(n2*2*ns, 1.0_R8,   cvsa_cl, 1, cvsa, 1)
       call b2saxpy(n2*2*ns, 1.0_R8, cvsahz_cl, 1, cvsahz, 1)
      end if
      if (fch_anomalous.ne.1.0_R8) then
        csig = csig * fch_anomalous
        calf = calf * fch_anomalous
        csig_an = csig_an * fch_anomalous
        calf_an = calf_an * fch_anomalous
      end if
      call b2saxpy (n2*2, 1.0_R8, chce0, 1, chce, 1)
      call b2saxpy (n2*2, 1.0_R8, chve0, 1, chve, 1)
      call b2saxpy (n2*2, 1.0_R8, chci0, 1, chci, 1)
      call b2saxpy (n2*2, 1.0_R8, chvi0, 1, chvi, 1)
      call b2saxpy (n2*2, 1.0_R8, csig0, 1, csig, 1)
      call b2saxpy (n2*2, 1.0_R8, calf0, 1, calf, 1)
*   ..compute terms from diamagnetic parts                        !srv 13.10.06 {
      do is = 0, ns-1
       if(is_neutral(is)) then
        cddi(:,:,0,is) = 0.0_R8
        cddi(:,:,1,is) = 0.0_R8
       else
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           cddi(ix,iy,0,is) = 0.0_R8
          else
           cddi(ix,iy,0,is) =
     &      +2.0e0_R8*(vol(leftix(ix,iy),leftiy(ix,iy))*
     &         bzb(leftix(ix,iy),leftiy(ix,iy))
     &         +vol(ix,iy)*bzb(ix,iy))/
     &       ((hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*
     &        (hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy)))/
     &        (rza(ix,iy,is)*qe)
          endif
         enddo
        enddo

        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) then
           cddi(ix,iy,1,is) = 0.0e0_R8
          else
           cddi(ix,iy,1,is) =
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

*  ..check that we have some particle diffusion
      if (ncall.le.3) then                                            !srv 12.03.12 {
        do is = 0, ns-1
          wrk(:,:) = cdna(:,:,0,is)+cdpa(:,:,0,is)
          wrkmin = huge(1.0_R8)
          do iy = -1, ny
            do ix = -1, nx
              if (leftix(ix,iy).ne.-2) then
                wrkmin = min(wrkmin, wrk(ix,iy))
              endif
            enddo
          enddo
          write(message,'(a, i3, a)')
     .     'No particle diffusion for species ',is,
     .     ' in poloidal direction'
          call xertst (wrkmin.gt.0.0_R8 .or. pol_anom_scale.eq.0.0_R8
     .     .or. .not.solving(2),message)
          wrk(:,:) = cdna(:,:,1,is)+cdpa(:,:,1,is)
          wrkmin = huge(1.0_R8)
          do iy = -1, ny
            do ix = -1, nx
              if (bottomiy(ix,iy).ne.-2) then
                wrkmin = min(wrkmin, wrk(ix,iy))
              endif
            enddo
          enddo
          write(message,'(a, i3, a)')
     .     'No particle diffusion for species ',is,
     .     ' in radial direction'
          call xertst (wrkmin.gt.0.0_R8 .or. .not.solving(2), message)
*  ..check that we have some anomalous ion heat transport
          wrkmin = minval(hcib(-1:nx,-1:ny,is))
          write(message,'(a, i3)')
     .     'No ion heat diffusivity for species ',is
          call xertst (wrkmin.gt.0.0_R8 .or. .not.solving(4), message)
        enddo
*  ..check that we have some anomalous perpendicular currents
        wrk(:,:) = calf_an(:,:,1)+csig_an(:,:,1)
        wrkmin = huge(1.0_R8)
        do iy = -1, ny
          do ix = -1, nx
            if (bottomiy(ix,iy).ne.-2) then
              wrkmin = min(wrkmin, wrk(ix,iy))
            endif
          enddo
        enddo
        call xertst (wrkmin.gt.0.0_R8 .or. pot_eq.ne.1 .or.
     .   fch_anomalous.eq.0.0_R8 .or. .not.solving(3) .or.
     .   no_current.eq.1, 'No anomalous perpendicular currents')
      endif                                                           !srv 12.03.12 }

      if(iout.ne.0) then                                              !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'cdpay'//chns)
        call my_out(70,nx,ny,cddi(-1,-1,0,is),'cddix'//chns)
        call my_out(70,nx,ny,cddi(-1,-1,1,is),'cddiy'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,0,is),'cvlax'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,1,is),'cvlay'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'cvsay'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,0,is),'cdnahzx'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,1,is),'cdnahzy'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,0,is),'cdpahzx'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,1,is),'cdpahzy'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,0,is),'cvlahzx'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,1,is),'cvlahzy'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'cvsahzy'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,0,is),'cvmahzx'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,1,is),'cvmahzy'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'cvsa_clx'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,1,is),'cvsa_cly'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,0,is),'cvsahz_clx'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,1,is),'cvsahz_cly'//chns)
        call my_out(70,nx,ny,cthe(-1,-1,is),'cthe'//chns)
        call my_out(70,nx,ny,cthi(-1,-1,is),'cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then
          call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'csiginx'//chns//
     &        '_'//chns0)
          call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'csiginy'//chns//
     &        '_'//chns0)
         endif
        enddo
       enddo
       call my_out(70,nx,ny,chce(-1,-1,0),'chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'chcey')
       call my_out(70,nx,ny,chci(-1,-1,0),'chcix')
       call my_out(70,nx,ny,chci(-1,-1,1),'chciy')
       call my_out(70,nx,ny,chve(-1,-1,0),'chvex')
       call my_out(70,nx,ny,chve(-1,-1,1),'chvey')
       call my_out(70,nx,ny,chvi(-1,-1,0),'chvix')
       call my_out(70,nx,ny,chvi(-1,-1,1),'chviy')
       call my_out(70,nx,ny,csig(-1,-1,0),'csigx')
       call my_out(70,nx,ny,csig(-1,-1,1),'csigy')
       call my_out(70,nx,ny,calf(-1,-1,0),'calfx')
       call my_out(70,nx,ny,calf(-1,-1,1),'calfy')
      endif                                                          !srv 17.06.08 }

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0 .or.
     &   iout_b2npc7.ne.0 .or. iout_b2npc11.ne.0) then               !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'b2tral_cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'b2tral_cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'b2tral_cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'b2tral_cdpay'//chns)
        if(.not.is_neutral(is)) then
         call my_out(70,nx,ny,cddi(-1,-1,0,is),'b2tral_cddix'//chns)
         call my_out(70,nx,ny,cddi(-1,-1,1,is),'b2tral_cddiy'//chns)
        endif
!        call my_out(70,nx,ny,cvla(-1,-1,0,is),'b2tral_cvlax'//chns)
!        call my_out(70,nx,ny,cvla(-1,-1,1,is),'b2tral_cvlay'//chns)
       enddo
      endif

      if(iout_b2npmo.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'b2tral_cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'b2tral_cvsahzy'//chns)
!        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'b2tral_cvsa_clx'//
!     &              chns)
!        call my_out(70,nx,ny,cvsa_cl(-1,-1,1,is),'b2tral_cvsa_cly'//
!     &              chns)
!        call my_out(70,nx,ny,cvsahz_cl(-1,-1,0,is),'b2tral_cvsahz_clx'//
!     &             chns)
!        call my_out(70,nx,ny,cvsahz_cl(-1,-1,1,is),'b2tral_cvsahz_cly'//
!     &             chns)
        call my_out(70,nx,ny,cthe(-1,-1,is),'b2tral_cthe'//chns)
        call my_out(70,nx,ny,cthi(-1,-1,is),'b2tral_cthi'//chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
       call my_out(70,nx,ny,chce(-1,-1,0),'b2tral_chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'b2tral_chcey')
       call my_out(70,nx,ny,chci(-1,-1,0),'b2tral_chcix')
       call my_out(70,nx,ny,chci(-1,-1,1),'b2tral_chciy')
!       call my_out(70,nx,ny,chve(-1,-1,0),'b2tral_chvex')
!       call my_out(70,nx,ny,chve(-1,-1,1),'b2tral_chvey')
!       call my_out(70,nx,ny,chvi(-1,-1,0),'b2tral_chvix')
!       call my_out(70,nx,ny,chvi(-1,-1,1),'b2tral_chviy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then
          call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'b2tral_csiginx'
     &        //chns//'_'//chns0)
          call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'b2tral_csiginy'
     &        //chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out(70,nx,ny,csig(-1,-1,0),'b2tral_csigx')
       call my_out(70,nx,ny,csig(-1,-1,1),'b2tral_csigy')
      endif                                                          !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tral

      end subroutine b2tral
