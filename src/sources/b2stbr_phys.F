*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr_phys (nx, ny, ns, nscx, iscx,
     &  dtim, BoRiS, vol, hx, hy, hz, gs, qz, pbs, bb, na, ua,
     &  te, ti, po, ne, ni, fna, fhe, fhi, fht, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0, sne0, tchem, tchee, tphys, tphye,
     &  thevp, thvpe, trese, tresn, trfln, trfle,
     &  neutral_sources_rescale, core_sources_rescale,
     &  sput_src, sput_chem_alpha, sput_phys_alpha, alpha,
     &  sput_chem_model, therm_evap, sput_res,
     &  reflection_on, sputter_energy_on, main_call, output,
     &  sput_frac_flag, sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1  sput_mixed_alpha, sput_mixed_beta,
     2  new_sputter_namelist, shi0_ff, f_redep,
     3  potential_at_guard_cell, recycled_neutrals_contr)
      use b2mod_types
      use b2mod_rates
      use b2mod_sputter
      use b2mod_tallies
      use b2mod_elements
      use b2mod_indirect
      use b2mod_constants
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_external
      use b2mod_diag
      use b2mod_geo
     , , only: hy1
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, iscx(0:nscxmax-1)
      integer sput_src, sput_chem_model,
     *     reflection_on, sputter_energy_on, output,
     *     sput_frac_flag, potential_at_guard_cell
      logical main_call, new_sputter_namelist
      real (kind=R8) :: dtim, BoRiS,
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  hz(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  gs(-1:nx,-1:ny,0:2), bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  fht(-1:nx,-1:ny,0:1), f_redep(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  neutral_sources_rescale, core_sources_rescale,
     &  sput_chem_alpha, sput_phys_alpha, alpha,
     &  therm_evap, sput_res,
     &  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1  sput_mixed_alpha, sput_mixed_beta,
     2  recycled_neutrals_contr
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1),
     &  shi0_ff(-1:nx,-1:ny,0:nscx-1),
     &  tchem(0:ns-1), tchee(0:ns-1), tphys(0:ns-1), tphye(0:ns-1),
     &  thevp(0:ns-1), thvpe(0:ns-1), trese(0:ns-1), tresn(0:ns-1),
     &  trfln(0:ns-1), trfle(0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBR_PHYS provides linearized expressions for the sources of
*     density, parallel momentum and heat, due to special boundary
*     conditions and the neutral gas recycling model, using the neutrals
*     namelist.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, i, ix, iy, is, ib, istra, ncall, iwall, itargsp, zbulk
!     integer k
      integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,            !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0, iout_b2npp7 = 0, iout_b2nppo = 0
      character chns*3
!     character chk*1
      character*15 hlp_frm
      integer, allocatable :: flux_x(:), flux_y(:)
      real (kind=R8), allocatable :: area(:)
      real (kind=R8) :: wrk0(-1:nx,-1:ny)
      real (kind=R8) ::
     *   amtol, t0, t1, t2, t3, b2stbr_phys_sna(0:ns-1),
     *   qpl, qex, surface_temp,
     *   Ei, Pflx, vbar, yldchm, yldres, GammaEvap, Phi_app,
     *   lnorm, deposition_error(ntrack), reflection_max,
     *   redep_max, yself_max
      logical l_chem_sput, lcore, match
      integer sput_dst_bulk
      parameter (amtol=1.0e-2_R8)
*   ..procedures
      intrinsic abs, exp
      external subini, subend, xertst
      external b2xvsg, b2xvff, smin, smax
      real (kind=R8) ::
     *   smax, smin
*   ..initialisation
      save ncall
      data ncall/0/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Note: This version of b2stbr is for the case when recycling is
*     calculated via the neutrals namelist.
*
*     For given species IS on boundary ISTRA:
*     The particle recycling coefficient is given by RECYC(is,istra)
*     The energy recycling coefficient is given by ERECYC(is,istra)
*     The recycling coefficient into the next ionized species is RCION(is,istra)
*     The chemical sputtering yield is given by chem_sput(is,istra)
*     The physical sputtering yield is computed from TRIM tables and multiplied
*     by phys_sput(is,istra)
*     Species with nuclear charge sput_src chemically sputter species
*     with index sput_dst
*
*     Description of some local variables:
*
*     amtol - real.
*     amtol specifies the mass difference (in amu) up to which two
*     species may be considered the same for the purpose of the
*     recycling model.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbr_phys')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                          !srv 22.08.08 {
       write(*,*) 'First call of b2stbr_phys'
       call ipgeti ('b2stbr_phys_iout', iout)                       !srv 11.07.05
       call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npc7_iout', iout_b2npc7)                     !srv 02.04.18
       call ipgeti ('b2npc11_iout', iout_b2npc11)                   !srv 02.04.18
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2npp7_iout', iout_b2npp7)
       call ipgeti ('b2nppo_iout', iout_b2nppo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2nphe_iout', iout_b2npht)                     !srv 11.09.09 }
      endif                                                         !srv 22.08.08 }
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test range of qz
       call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &   smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     &   'faulty argument range qz(,,0)')
       call xertst ( 0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &   smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     &   'faulty argument range qz(,,1)')
*    ..test sign of pbs
!       call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                  !srv 11.01.13
*    ..test edge values of pbs
       call b2xvff (nx, ny, pbs, 'pbs')
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif
*   ..test edge values of fluxes
      do is = 0, ns-1
       call b2xvff (nx, ny, fna(-1,-1,0,is), 'fna')
      enddo
      call b2xvff (nx, ny, fhe, 'fhe')
      call b2xvff (nx, ny, fhi, 'fhi')
      call b2xvff (nx, ny, fch, 'fch')
      call b2xvff (nx, ny, fne, 'fne')
      call b2xvff (nx, ny, fni, 'fni')

      if(new_sputter_namelist) then
         do itargsp=1, ntargsp
            if(track_chem_sput(itargsp)) then
               if(nint(zn(target_intcode(itargsp))).ne.6) then
                  write(*,*) 'b2stbr_phys: chemical sputtering turned ',
     &                 'off for tracked species ',itargsp,
     &                 ' which has nuclear charge ',
     &                 zn(target_intcode(itargsp)),
     &                 ' not 6 (i.e. is not C)'
                  track_chem_sput(itargsp)=.false.
               endif
            endif
            write(*,*) 'b2stbr_phys: chemical sputtering ',
     &           'for tracked species ',itargsp,
     &           ' (nuclear charge ',zn(target_intcode(itargsp)),
     &           ' ) is ',track_chem_sput(itargsp)
         enddo
      endif

      b2stbr_phys_sna=0.0_R8
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            b2stbr_phys_sna(is)=b2stbr_phys_sna(is)+
     1       sna0(ix,iy,0,is)+sna0(ix,iy,1,is)*na(ix,iy,is)
          enddo
        enddo
      enddo
      if((main_call.and.output.ge.1).or.output.gt.1)
     *     write(*,*) 'b2stbr_phys_sna(start): ',b2stbr_phys_sna

      if(ncall.eq.0.and.sput_dst.ge.0.and.sput_dst.lt.ns) then
        do itargsp=1,ntargsp
          deposition_error(itargsp)=0.0_R8
          lnorm=am(target_intcode(itargsp))/tgdns(itargsp)/
     &     avogr/1.0e6_R8/monolayer_l(itargsp)
          do iwall=1, nbnd
            if(deposition(iwall,itargsp).lt.0.0_R8) then
              if(output.ge.8) then
                write(*,*) 'itrack,iwall,deposition_error',
     1           itargsp,iwall,deposition(iwall,itargsp)
              endif
              deposition_error(itargsp)=deposition_error(itargsp)-
     1         deposition(iwall,itargsp)
              deposition(iwall,itargsp)=0.0_R8
            endif
            if(plate_area(iwall).gt.0.0_R8) then
              monolayer_deposition(iwall,itargsp)=
     &         deposition(iwall,itargsp)/plate_area(iwall)*lnorm
              monolayer_erosion(iwall,itargsp)=
     &         erosion(iwall,itargsp)/plate_area(iwall)*lnorm
            endif
          enddo
          if(deposition_error(itargsp).ne.0.0_R8) then
            write(*,*) 'itrack,deposition_error ',
     1       itargsp,deposition_error(itargsp)
          endif
        enddo
      endif

      reflection_max = 0.0_R8
      redep_max = 0.0_R8
      yself_max = 0.0_R8

      do istra = 1, nstrat

        sput_dst_bulk=targsp(istra,1)
        if (sput_dst_bulk.eq.-1) then
          l_chem_sput=.false.
          zbulk=0
        else
          zbulk=nint(zn(sput_dst_bulk))
          l_chem_sput=(sput_dst_bulk.eq.sput_dst.or.chemsp(istra)).and.
     &       zbulk.eq.6
        endif
        if(ncall.eq.0) then
          if(sput_frac_flag.gt.0 .and. istra.eq.1) write(*,*)
     &      'Mixed materials sputtering model is active!'
          if (crcstra(istra).eq.'S' .or. crcstra(istra).eq.'N' .or.
     &        crcstra(istra).eq.'E' .or. crcstra(istra).eq.'W') then
            if (sput_dst_bulk.eq.-1) then
              write(*,*)
     &         'No bulk material sputtering specified for stratum ',
     &          istra
            else if(sput_frac_flag.eq.0) then
              write(*,'(a,a2,a,i2,a,l1)')
     &         'Sputtering ',is_codes(sput_dst_bulk),
     &         ' for stratum ',istra,
     &         ' ; chemical sputtering is ',l_chem_sput
            else
              write(*,'(a,i2,a,a2,a,l1)')
     &         'Bulk material for stratum ',istra,
     &         ' is ',is_codes(sput_dst_bulk),
     &         ' ; chemical sputtering is ',l_chem_sput
              do i = 1, rc_list_size(istra)
                ix = rc_list_x(i,istra)
                iy = rc_list_y(i,istra)
                ib = xymap(ix,iy)
                match = .false.
                is = 1
                do while (.not.match .and.
     &           is.le.matsurf_nconstituents(ib))
                  match = zbulk.eq.matsurf_nzconstituents(ib,is)
                  is = is+1
                enddo
                if(.not.match) then
                  write(*,*) 'Could not find ', zbulk, ' in ',
     &             matsurf_nzconstituents(
     &                 ib,1:matsurf_nconstituents(ib)),
     &             ' for surface element ', ib
                endif
                call xertst (match,
     &           'Mismatch between target material declarations!')
              enddo
            endif
          endif
        endif

        b2stbr_phys_sna=0.0_R8
        do is=0,ns-1
          do iy=-1,ny
            do ix=-1,nx
              b2stbr_phys_sna(is)=b2stbr_phys_sna(is)+
     1         sna0(ix,iy,0,is)+sna0(ix,iy,1,is)*na(ix,iy,is)
            enddo
          enddo
        enddo
        if((main_call.and.output.ge.1).or.output.gt.1)
     *     write(*,*) 'stratum,b2stbr_phys_sna(start): ',
     1   istra, b2stbr_phys_sna

* ..recycling model at south boundary
        if (crcstra(istra).eq.'S') then

*   ..loop over boundary cells
          do i = 1, rc_list_size(istra)
            Phi_app = get_phi_applied(boundary_namelist,
     .       rc_list_x(i,istra),rc_list_y(i,istra))
            lcore = on_closed_surface(rc_list_x(i,istra),
     .                                rc_list_y(i,istra)) ! xpb: no sputtering in core!
*   ..compute 2-D target temperature evolution
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. i.eq.1) then
              allocate(flux_x(rc_list_size(istra)),
     .                 flux_y(rc_list_size(istra)),
     .                 area(rc_list_size(istra)))
              do ix = 1, rc_list_size(istra)
                flux_x(ix) = topix(rc_list_x(ix,istra),
     .                             rc_list_y(ix,istra))
                flux_y(ix) = topiy(rc_list_x(ix,istra),
     .                             rc_list_y(ix,istra))
                area(ix) = gs(flux_x(ix),flux_y(ix),1)
              enddo
              call prehd2D(ns,rc_list_size(istra),
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),dtim,-1,
     &         flux_x(1:rc_list_size(istra)),
     &         flux_y(1:rc_list_size(istra)),1,
     &         area(1:rc_list_size(istra)),BoRiS)
              deallocate(flux_x,flux_y,area)
            endif

            call recycle(
     &           species_start(istra), species_end(istra),
     &           rc_list_x(i,istra),rc_list_y(i,istra),
     &           topix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           topiy(rc_list_x(i,istra),rc_list_y(i,istra)),-1,
     &           topix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           topiy(rc_list_x(i,istra),rc_list_y(i,istra)),1,
     &           hy1(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           hy1(topix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &               topiy(rc_list_x(i,istra),rc_list_y(i,istra))),
     &           Phi_app,
     &           gs(topix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &              topiy(rc_list_x(i,istra),rc_list_y(i,istra)),1),
     &           lcore)
          enddo

* ..recycling model at north boundary
        else if (crcstra(istra).eq.'N') then

*   ..loop over boundary cells
          do i = 1, rc_list_size(istra)
            Phi_app = get_phi_applied(boundary_namelist,
     .       rc_list_x(i,istra),rc_list_y(i,istra))
            lcore=on_closed_surface(rc_list_x(i,istra),
     .                              rc_list_y(i,istra))
*   ..compute 2-D target temperature evolution
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. i.eq.1) then
              allocate(area(rc_list_size(istra)))
              do ix = 1, rc_list_size(istra)
                area(ix) = gs(rc_list_x(ix,istra),rc_list_y(ix,istra),1)
              enddo
              call prehd2D(ns,rc_list_size(istra),
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),dtim,+1,
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),1,
     &         area(1:rc_list_size(istra)),BoRiS)
              deallocate(area)
            endif

            call recycle(
     &           species_start(istra),species_end(istra),
     &           rc_list_x(i,istra),rc_list_y(i,istra),
     &           bottomix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           bottomiy(rc_list_x(i,istra),rc_list_y(i,istra)),+1,
     &           rc_list_x(i,istra),rc_list_y(i,istra),1,
     &           hy1(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           hy1(bottomix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &               bottomiy(rc_list_x(i,istra),rc_list_y(i,istra))),
     &           Phi_app,
     &           gs(rc_list_x(i,istra),rc_list_y(i,istra),1),
     &           lcore)
          enddo

* ..recycling model at west boundary
        else if (crcstra(istra).eq.'W') then

*   ..loop over boundary cells
          do i = 1, rc_list_size(istra)
            Phi_app = get_phi_applied(boundary_namelist,
     .       rc_list_x(i,istra),rc_list_y(i,istra))
            lcore=.false.
*   ..compute 2-D target temperature evolution
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. i.eq.1) then
              allocate(flux_x(rc_list_size(istra)),
     .                 flux_y(rc_list_size(istra)),
     .                 area(rc_list_size(istra)))
              do ix = 1, rc_list_size(istra)
                flux_x(ix) = rightix(rc_list_x(ix,istra),
     .                               rc_list_y(ix,istra))
                flux_y(ix) = rightiy(rc_list_x(ix,istra),
     .                               rc_list_y(ix,istra))
                area(ix) = gs(flux_x(ix),flux_y(ix),0)
              enddo
              call prehd2D(ns,rc_list_size(istra),
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),dtim,-1,
     &         flux_x(1:rc_list_size(istra)),
     &         flux_y(1:rc_list_size(istra)),0,
     &         area(1:rc_list_size(istra)),BoRiS)
              deallocate(flux_x,flux_y,area)
            end if

            call recycle(
     &           species_start(istra),species_end(istra),
     &           rc_list_x(i,istra),rc_list_y(i,istra),
     &           rightix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           rightiy(rc_list_x(i,istra),rc_list_y(i,istra)),-1,
     &           rightix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           rightiy(rc_list_x(i,istra),rc_list_y(i,istra)),0,
     &           hx(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           hx(rightix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &              rightiy(rc_list_x(i,istra),rc_list_y(i,istra))),
     &           Phi_app,
     &           gs(rightix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &              rightiy(rc_list_x(i,istra),rc_list_y(i,istra)),0),
     &           lcore)
          enddo

* ..recycling model at east boundary
        else if (crcstra(istra).eq.'E') then

*   ..loop over boundary cells
          do i = 1, rc_list_size(istra)
            Phi_app = get_phi_applied(boundary_namelist,
     .       rc_list_x(i,istra),rc_list_y(i,istra))
            lcore=.false.
*   ..compute 2-D target temperature evolution
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. i.eq.1) then
              allocate(area(rc_list_size(istra)))
              do ix = 1, rc_list_size(istra)
                area(ix) = gs(rc_list_x(ix,istra),rc_list_y(ix,istra),0)
              enddo
              call prehd2D(ns,rc_list_size(istra),
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),dtim,+1,
     &         rc_list_x(1:rc_list_size(istra),istra),
     &         rc_list_y(1:rc_list_size(istra),istra),0,
     &         area(1:rc_list_size(istra)),BoRiS)
              deallocate(area)
            endif

            call recycle(
     &           species_start(istra), species_end(istra),
     &           rc_list_x(i,istra),rc_list_y(i,istra),
     &           leftix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           leftiy(rc_list_x(i,istra),rc_list_y(i,istra)),+1,
     &           rc_list_x(i,istra),rc_list_y(i,istra),0,
     &           hx(rc_list_x(i,istra),rc_list_y(i,istra)),
     &           hx(leftix(rc_list_x(i,istra),rc_list_y(i,istra)),
     &              leftiy(rc_list_x(i,istra),rc_list_y(i,istra))),
     &           Phi_app,
     &           gs(rc_list_x(i,istra),rc_list_y(i,istra),0),
     &           lcore)
          enddo

        endif       ! crcstra(istra)

        b2stbr_phys_sna=0.0_R8
        do is=0,ns-1
          do iy=-1,ny
            do ix=-1,nx
              b2stbr_phys_sna(is)=b2stbr_phys_sna(is)+
     1         sna0(ix,iy,0,is)+sna0(ix,iy,1,is)*na(ix,iy,is)
            enddo
          enddo
        enddo
        if((main_call.and.output.ge.1).or.output.gt.1)
     *       write(*,*) 'stratum,b2stbr_phys_sna(end): ',
     1   istra, b2stbr_phys_sna

      enddo      ! istra

      b2stbr_phys_sna=0.0_R8
      do is=0,ns-1
        do iy=-1,ny
          do ix=-1,nx
            b2stbr_phys_sna(is)=b2stbr_phys_sna(is)+
     1       sna0(ix,iy,0,is)+sna0(ix,iy,1,is)*na(ix,iy,is)
          enddo
        enddo
      enddo
      if((main_call.and.output.ge.1).or.output.gt.1)
     *     write(*,*) 'b2stbr_phys_sna(end): ',b2stbr_phys_sna

      if(main_call.and.sput_dst.ge.0.and.sput_dst.lt.ns) then
        do itargsp=1,ntargsp
          deposition_error(itargsp)=0.0_R8
          lnorm=am(target_intcode(itargsp))/tgdns(itargsp)/
     &     avogr/1.0e6_R8/monolayer_l(itargsp)
          do iwall=1, nbnd
            if(deposition(iwall,itargsp).lt.0.0_R8) then
              if(output.ge.8) then
                write(*,*) 'itrack,iwall,deposition_error',
     1           itargsp,iwall,deposition(iwall,itargsp)
              endif
              deposition_error(itargsp)=deposition_error(itargsp)-
     1         deposition(iwall,itargsp)
              deposition(iwall,itargsp)=0.0_R8
            endif
            if(plate_area(iwall).gt.0.0_R8) then
              monolayer_deposition(iwall,itargsp)=
     &         deposition(iwall,itargsp)/plate_area(iwall)*lnorm
              monolayer_erosion(iwall,itargsp)=
     &         erosion(iwall,itargsp)/plate_area(iwall)*lnorm
            endif
          enddo
          if(deposition_error(itargsp).ne.0.0_R8) then
            write(*,*) 'itrack,deposition_error ',
     1       itargsp,deposition_error(itargsp)
          endif
        enddo
      endif

      if (iout.eq.1 .or. iout_b2npmo.eq.1) then                   !srv 11.07.05 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!         do k=0,3
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,smo0(-1,-1,k,is),
!     &                       'b2stbr_phys_smo'//chk//chns)
!         enddo
        wrk0=smo0(:,:,0,is)+
     &       smo0(:,:,1,is)*ua(:,:,is)+
     &       smo0(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &       smo0(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
        call my_out(70,nx,ny,wrk0,'b2stbr_phys_smo'//chns)
       enddo
      endif
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1 .or.
     &                    iout_b2npc7.eq.1 .or. iout_b2npc11.eq.1)) then !srv 02.04.18
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!         do k=0,1
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,sna0(-1,-1,k,is),
!     &                       'b2stbr_phys_sna'//chk//chns)
!         enddo
        wrk0=sna0(:,:,0,is)+sna0(:,:,1,is)*na(:,:,is)
       call my_out(70,nx,ny,wrk0,'b2stbr_phys_sna'//chns)
       enddo
      endif
      if (iout.eq.1 .or. (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1)) then
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,she0(-1,-1,k),'b2stbr_phys_she'//chk)
!         call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbr_phys_shi'//chk)
!         call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbr_phys_sch'//chk)
!        enddo
         wrk0=she0(:,:,0)+she0(:,:,1)*te+
     &        she0(:,:,2)*ne + she0(:,:,3)*te*ne
         call my_out(70,nx,ny,wrk0,'b2stbr_phys_she')
         wrk0=shi0(:,:,0)+shi0(:,:,1)*ti+
     &        shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0)
         call my_out(70,nx,ny,wrk0,'b2stbr_phys_shi')
      endif
      if (iout.eq.1 .or. (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)) then
         wrk0=sch0(:,:,0)+sch0(:,:,1)*po+
     &        sch0(:,:,2)*ne+sch0(:,:,3)*po*ne(:,:)
         call my_out(70,nx,ny,wrk0,'b2stbr_phys_sch')
      endif
      if(iout.eq.1) then
        open(70,file='t_files.dat')
        write (hlp_frm,'(a4,i4,a7)') '(a7,',ns,'es13.5)'               !srv 11.03.13 {
        write(70,hlp_frm) ' tchem=',tchem
        write(70,hlp_frm) ' tchee=',tchee
        write(70,hlp_frm) ' tphys=',tphys
        write(70,hlp_frm) ' tphye=',tphye
        write(70,hlp_frm) ' thevp=',thevp
        write(70,hlp_frm) ' thvpe=',thvpe
        write(70,hlp_frm) ' trese=',trese
        write(70,hlp_frm) ' tresn=',tresn
        write(70,hlp_frm) ' trfln=',trfln
        write(70,hlp_frm) ' trfle=',trfle                              !srv 11.03.13 }
        close(70)
      endif                                                            !srv 11.07.05 }

      if(reflection_max.gt.1.0_R8 .or.
     .  ((main_call.and.output.ge.1).or.output.gt.1))
     1 write(*,*) 'Maximum reflection coefficient ', reflection_max
      if (redep_max.gt.0.0_R8 .and. (main_call.or.output.gt.1))
     . write(*,*) 'Maximum redeposition coefficient ', redep_max
      if(redep_max.gt.0.0_R8 .and. (yself_max.ge.1.0_R8 .or.
     .  ((main_call.and.output.ge.1).or.output.gt.1)))
     . write(*,*) 'Maximum self-sputtering coefficient ', yself_max
* ..return
      ncall = ncall + 1
      call subend ()
      return

*-----------------------------------------------------------------------

      contains

      subroutine recycle(nsmin,nsmax,ix,iy,ixn,iyn,
     &     isign,ixf,iyf,idir,
     &     hh,hhn,Phi_app,area,
     &     lcore)
c ix,iy           position of the cell
c ixn,iyn         position of the neighbour cell
c isign           multiplies fluxes to get the right sign
c ixf,iyf         position where fluxes are evaluated
c idir            index for direction of flux (x or y)
c hh              either hx or hy1 @ ix,iy
c hhn             either hx or hy1 @ ixn,iyn
c Phi_app         applied potential
c area            area element
c lcore           is this a core element?
      implicit none
      integer nsmin, nsmax, ix, iy, ixn, iyn, isign, ixf, iyf, idir
      real (kind=R8) hh, hhn, Phi_app, area, pof
      logical lcore, found

      integer is0, is1, is2, is, js, k, iwall, itrack
      real (kind=R8) :: tef, tif, taf, zaf(0:ns-1), ptf,
     &     sput_frac(0:ntrack), sput_frac_chem(0:ntrack), yself, eself,
     &     sput_frac_chem_sum, be_frac, be_suppression_factor,
     &     ytmp, ztmp, wrk(ndepth_nml)
      real (kind=R8) mass_density, local_reflection
      character*13 hlp_frm
      external mass_density
      logical lnext
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol

      iwall=xymap(ix,iy)
cdpc-exp
      if(sput_frac_flag.eq.0.or.lcore.or.ntargsp.eq.0.or.
     .   iwall.eq.0) then
c no mixed materials
        sput_frac=0.0_R8
        sput_frac(0)=1.0_R8
        sput_frac_chem=0.0_R8
        if(l_chem_sput) then
          sput_frac_chem(0)=1.0_R8
        else
          sput_frac_chem(0)=0.0_R8
        endif
        sput_frac_chem_sum=sput_frac_chem(0)
      else  ! mixed materials
c calculate the fraction from each constituent of the mixed material
        if(output.ge.10)
     &       write(*,*) 'b2stbr_phys.recycle: ix,iy,iwall',ix,iy,iwall
c the base material ("0") is assumed to have an unnormalised weight
c sput_mixed_alpha (default 1)
        sput_frac(0)=sput_mixed_alpha
c the remaining constituents each have an unnormalised weight equal
c to the number of monolayers deposited raised to the power
c sput_mixed_beta (default 1)
        do itrack=1,ntargsp
          sput_frac(0)=sput_frac(0)+
     1     monolayer_deposition(iwall,itrack)**sput_mixed_beta
        enddo
c the remaining constituents are then normalised by the sum of weights
        do itrack=1,ntargsp
          sput_frac(itrack)=
     1      monolayer_deposition(iwall,itrack)**sput_mixed_beta/
     &      sput_frac(0)
        enddo
c the base material weight is calculated to ensure that the sum of
c normalised weights is 1
        sput_frac(0)=1.0_R8
        do itrack=1,ntargsp
          sput_frac(0)=sput_frac(0)-sput_frac(itrack)
        enddo
c now calculate the fractions for chemical sputtering
c the base chemically sputters if l_chem_sput is true
        if(l_chem_sput) then
          sput_frac_chem(0)=sput_frac(0)
        else
          sput_frac_chem(0)=0.0_R8
        endif
        sput_frac_chem_sum=sput_frac_chem(0)
c now check each tracked deposited species
        do itrack=1,ntargsp
          if(itrack.eq.track_index(sput_dst).or.
     &        track_chem_sput(itrack)) then
            sput_frac_chem(itrack)=sput_frac(itrack)*
     1         chemical_erosion_redep_fac
          else
            sput_frac_chem(itrack)=0.0_R8
          endif
          sput_frac_chem_sum=sput_frac_chem_sum+sput_frac_chem(itrack)
        enddo
c Be suppression of chemical sputtering
        if(chemical_erosion_be_fac) then
          be_frac=0.0_R8
          if(zbulk.eq.4) be_frac=be_frac+sput_frac(0)
          do itrack=1,ntargsp
            if(nint(zn(target_intcode(itrack))).eq.4)
     1        be_frac=be_frac+sput_frac(itrack)
          enddo
          call xertst(0.0_R8.le.be_frac.and.
     1     be_frac.le.1.0_R8,'be_frac out of range')
          be_suppression_factor=
     1     1.0_R8-chemical_erosion_be_fac_c/2.0_R8*
     2     (tanh((be_frac-chemical_erosion_be_fac_a)/
     3     chemical_erosion_be_fac_b)-
     4      tanh((-chemical_erosion_be_fac_a)/
     5     chemical_erosion_be_fac_b))
          call xertst(0.0_R8.le.be_suppression_factor.and.
     1      be_suppression_factor.le.1.0_R8,
     2     'be_suppression_factor out of range')
          sput_frac_chem_sum=sput_frac_chem_sum*be_suppression_factor
          sput_frac_chem(1:ntargsp)=sput_frac_chem(1:ntargsp)*
     1     be_suppression_factor
          if(output.ge.10) then
            write(*,'(a,1p,2g12.5)')
     1       'b2stbr_phys.recycle: be_frac, be_suppression_factor',
     2       be_frac, be_suppression_factor
          endif
        endif
        if(output.ge.10) then
          write(hlp_frm,'(a6,i1,a6)') '(a,1p,',ntargsp+1,'g12.5)'
          write(*,hlp_frm) 'b2stbr_phys.recycle: sput_frac ',
     &          sput_frac(0:ntargsp)
          write(*,hlp_frm) 'b2stbr_phys.recycle: sput_frac_chem',
     &         sput_frac_chem(0:ntargsp)
          write(*,*) 'b2stbr_phys.recycle: sput_frac_chem_sum',
     &         sput_frac_chem_sum
        endif
cxpb  Contributions from the bulk
        layer_relconstituents(iwall,1:matsurf_nconstituents(iwall))=
     1    matsurf_relconstituents
     .    (iwall,1:matsurf_nconstituents(iwall))*
     2   sput_frac(0)
cxpb  Contributions from the deposited species
        do itrack=1,ntargsp
           layer_nzconstituents
     .      (iwall,matsurf_nconstituents(iwall)+itrack)=
     1      nint(zn(target_intcode(itrack)))
           layer_relconstituents
     .      (iwall,matsurf_nconstituents(iwall)+itrack)=
     .      sput_frac(itrack)
        enddo
cxpb  Update layer thickness (deposited - eroded) and layer depth
cxpb  (sum of deposition), minimum depth over which material mixing occurs
        layer_thickness(iwall)=0.0_R8
        layer_depth(iwall)=0.0_R8
        do itrack=1,ntargsp
          layer_thickness(iwall)=layer_thickness(iwall)+
     &      (monolayer_deposition(iwall,itrack)-
     &       monolayer_erosion(iwall,itrack))*monolayer_l(itrack)
          layer_depth(iwall)=layer_depth(iwall)+
     &       monolayer_deposition(iwall,itrack)*monolayer_l(itrack)
        enddo
        call add_layer_thickness(ns,iwall)
        call compose_layer(iwall,sput_frac)
cxpb  Update layer mass density
cxpb  with bulk and deposited contributions
        layer_mass_density(iwall)=mass_density(iwall,2)*sput_frac(0)
        do itrack=1,ntargsp
          layer_mass_density(iwall)=layer_mass_density(iwall)+
     &     sput_frac(itrack)*(tgdns(itrack)*1.0e3_R8)
        enddo
      endif
cdpc-end
      tef = (te(ix,iy)*hhn+te(ixn,iyn)*hh)/(hhn+hh)
      tif = (ti(ix,iy)*hhn+ti(ixn,iyn)*hh)/(hhn+hh)
      if (potential_at_guard_cell.eq.1) then
!        write(6,*) 'Potential taken at guard cell (phys)'
        pof = po(ix,iy)
      else
!        write(6,*) 'Potential taken at real cell (phys)'
        pof = po(ixn,iyn)
      endif
      zaf(0:ns-1) =
     & (rza(ix,iy,0:ns-1)*hhn+rza(ixn,iyn,0:ns-1)*hh)/(hhn+hh)
      qpl = max(0.0_R8,isign*fht(ixf,iyf,idir))/area
      qex = max(0.0_R8,isign*fhi_ext(ixf,iyf,idir))/area
      do is=0,ns_ext-1
         ptf =
     &    (pt_ext(ix,iy,is)*hhn+pt_ext(ixn,iyn,is)*hh)/(hhn+hh)
         taf =
     &    (ta_ext(ix,iy,is)*hhn+ta_ext(ixn,iyn,is)*hh)/(hhn+hh)
         qex = qex + max(0.0_R8,isign*fa_ext(ixf,iyf,idir,is))
     &        *((0.5_R8*am_ext(is)*mp*ua_ext(ix,iy,is)**2+taf)*
     &          (1.0_R8-BoRiS)
     &          + ptf*ev)/area
      enddo
      if(lcore) then
         surface_temp = 0.0_R8
      else
         if(iwall.eq.0) then
            write(*,*) 'xymap returned zero for ix,iy = ',ix,iy
            call xerrab ('Invalid iwall setting!')
         endif
         if (plate_model.eq.0) then
            call heatdiff0D(iwall,qpl+qex,surface_temp)
            target_temp(iwall,1) = surface_temp
         else if (plate_model.eq.1) then
            if (main_call) then
              wrk(1:ndepth_nml) = target_temp(iwall,1:ndepth_nml)
              call heatdiff1D(iwall,ns,qpl+qex,
     .           dtim*plate_time_factor(iwall),wrk,output)
              target_temp(iwall,1:ndepth_nml) = wrk(1:ndepth_nml)
            end if
            surface_temp = target_temp(iwall,1)
         else if (plate_model.eq.2) then
           surface_temp = target_temp(iwall,1)
         endif
      endif

*    ..loop over species
      is0 = 0
      is1 = 0
      do is = nsmin, nsmax
         if (zn(is0).ne.zn(is).or.amtol.lt.abs(am(is0)-am(is))) then
            is0 = is
         endif
cdpc
         if(is_neutral(is)) is1=is
         if(zn(is1).ne.zn(is).or.amtol.lt.abs(am(is1)-am(is))) then
            write(*,*) 'b2stbr problem: is, is0, is1 ',is,is0,is1
         endif
         if(is0.ne.is1) then
            if(ncall.eq.0)
     1           write(*,*) 'b2stbr: setting is0 to is1 ',is0,is1
            is0=is1
         endif
         is2=-1
         k=0
         found=.false.
         do while (.not.found .and. k.lt.nscx)
           if (is0.eq.iscx(k)) then
             found = .true.
             is2 = k
           end if
           k = k + 1
         enddo
cdpc
*       (species (is) recycles into is0)
*       (is0 corresponds to CX species is2)
         t0 = max(0.0_R8,isign*fna(ixf,iyf,idir,is)) *
     1        neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         Ei=max(0.0_R8,tif/ev+(pof-Phi_app)*zaf(is))
                                 ! xpb: Incident energy (eV)
         Pflx = t0/area          ! xpb: Incident particle flux (part/m^2/s)
         if (is_neutral(is)) then
            vbar = sqrt(8.0_R8 * ti(ix,iy) / (pi * am(is) * mp))
            Pflx = Pflx + na(ix,iy,is)*vbar*alpha *
     1           neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         endif
         sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0)+recyc(is,istra)*t0
         if (is_neutral(is0)) then                 ! IYS 18.09.2018
           sna_stbr_neut(b2espcr(is0),lbcstra(istra)) =
     &      sna_stbr_neut(b2espcr(is0),lbcstra(istra)) +
     &      recyc(is,istra)*t0
           if (output.gt.2) write(*,*) 'Recycling: ',
     &      'Contribution to sna_stbr_neut. is0 = ',
     &       is0,' istra = ',istra,' ix = ',ix,' iy = ',iy
         end if

         itrack=track_index(is0)
cxpb Whatever is not recycled is counted as deposition
         if(itrack.gt.0.and.itrack.le.ntrack
     1        .and.iwall.gt.0.and.main_call) then
            deposition(iwall,itrack)=deposition(iwall,itrack)+
     1           dtim*plate_time_factor(iwall)*
     2           (1.0_R8-recyc(is,istra))*t0
         endif
cwdk  recycle toroidal componenent of parallel momentum, factor mrecyc
c     (toroidal piec only, because poloidal piece of parallel velocity flips sign after reflection)
         smo0(ix,iy,0,is0) = smo0(ix,iy,0,is0) + mrecyc(is,istra)*
     &        am(is)*mp*ua(ix,iy,is)*bb(ix,iy,2)**2/bb(ix,iy,3)**2*
     &        recyc(is,istra)*t0*hz(ix,iy)
         shi0(ix,iy,0) = shi0(ix,iy,0)+recycled_neutrals_contr*max(
     &        erecyc(is,istra)*(1.5_R8+BoRiS)*Ei*ev,
     &        (1.5_R8+BoRiS)*surface_temp*KBolt)*recyc(is,istra)*t0
         if (is2.ne.-1) shi0_ff(ix,iy,is2) = shi0_ff(ix,iy,is2)+
     &        recycled_neutrals_contr*max(
     &        erecyc(is,istra)*(1.5_R8+BoRiS)*Ei*ev,
     &        (1.5_R8+BoRiS)*surface_temp*KBolt)*recyc(is,istra)*t0
         ! xpb: recycled neutrals come back with at least thermal energy
         if(.not.lcore) then    ! xpb: no sputtering in core!
c chemical sputtering
            if(chem_sput(is,istra).gt.0.0_R8.and.
     &           nint(zn(is)).eq.sput_src.and.
     &           (l_chem_sput.or.sput_frac_flag.ne.0)) then
               if(zbulk.eq.6.or.sput_frac_flag.ne.0) then      ! chemical sputtering only for C
                  if(sput_chem_model.eq.1) then
                     if (Ei.gt.0.0_R8) then
                        ytmp = chem_sput(is,istra) *
     &                       exp(-(sput_chem_cutoff_alpha/Ei)**
     &                       sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
                        t1=t0 * ytmp
                        chemical_sputtering(iwall,is) = ytmp
                     else
                        t1=0.0_R8
                        chemical_sputtering(iwall,is) = 0.0_R8
                     endif
                  else if (zbulk.eq.6) then
                     call chmsput(iwall,Ei,zaf(is),zn(is),
     &                    zn(sput_dst_bulk),
     &                    am(is),am(sput_dst_bulk),surface_temp,
     &                    Pflx/neutral_sources_rescale,yldchm)
                     t1=Pflx*neutral_sources_rescale*
     &                    yldchm*chem_sput(is,istra)*area
                     chemical_sputtering(iwall,is) = yldchm
                  else
                     t1=0.0_R8
                     chemical_sputtering(iwall,is) = 0.0_R8
                  endif
               else    ! not C
                  t1=0.0_R8
                  chemical_sputtering(iwall,is) = 0.0_R8
               endif
               if(sput_frac_flag.ne.0) then
                  if (zbulk.ne.6 .and. sput_chem_model.ne.1) then
                     call chmsput(iwall,Ei,zaf(is),zn(is),
     &                    6.0_R8,am(is),12.011_R8,surface_temp,
     &                    Pflx/neutral_sources_rescale,yldchm)
                     t1=Pflx*neutral_sources_rescale*
     &                    yldchm*chem_sput(is,istra)*area
                     chemical_sputtering(iwall,is) = yldchm
                  endif
                  t1=t1*sput_frac_chem_sum
                  chemical_sputtering(iwall,is) =
     &                 chemical_sputtering(iwall,is)*sput_frac_chem_sum
               endif
               if(ncall.ne.0 .or. sputter_namelist_used) then
                  t1=(sputter_yield(ix,iy,is,1)+sput_chem_alpha*t1)/
     *                 (1+sput_chem_alpha)
               endif
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               sputter_yield(ix,iy,is,1)=t1
               tchem(is)=tchem(is)+t1
               if (sput_dst_bulk.ne.-1) then
                 sna0(ix,iy,0,sput_dst_bulk) =
     &              sna0(ix,iy,0,sput_dst_bulk) + t1
                 sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra)) =
     &            sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra))
     &            + t1
                 if (output.gt.2) write(*,*)
     &            'Chemical sputtering from bulk: ',
     &            'Contribution to sna_stbr_sput. is = ',sput_dst_bulk,
     &            ' istra = ',istra,' ix = ',ix,' iy = ',iy
               endif
cxpb  Erosion from the bulk material if chemically sputtered
cxpb  without mixed materials model
               if(sput_frac_flag.eq.0) then
                  if (sput_dst_bulk.eq.-1) then
                    itrack=0
                  else
                    itrack=track_index(sput_dst_bulk)
                  endif
                  if(itrack.gt.0.and.itrack.le.ntrack.and.
     1                 iwall.gt.0.and.main_call)
     1                 erosion(iwall,itrack) = erosion(iwall,itrack) +
     1                 dtim*plate_time_factor(iwall)*t1
               else
                  if(iwall.gt.0.and.main_call
     &                 .and.sput_frac_chem_sum.gt.0.0_R8) then
                     if (sput_dst_bulk.eq.-1) then
                       itrack=0
                     else
                       itrack=track_index(sput_dst_bulk)
                     endif
cxpb  Erosion from the bulk material if chemically sputtered
cxpb  with mixed materials model
                     if(itrack.gt.0.and.itrack.le.ntrack)
     &                    erosion(iwall,itrack) =
     &                    erosion(iwall,itrack) +
     1                    dtim*plate_time_factor(iwall)*t1*
     &                    sput_frac_chem(0)/sput_frac_chem_sum
cxpb  Removal from the deposited material if chemically sputtered
cxpb  with mixed materials model
                     do itrack=1,ntargsp
                        if(itrack.gt.0.and.itrack.le.ntrack)
     &                       deposition(iwall,itrack) =
     &                       deposition(iwall,itrack) -
     1                       dtim*plate_time_factor(iwall)*t1*
     &                       sput_frac_chem(itrack)/sput_frac_chem_sum
                     enddo
                  endif
               endif
               if(sputter_energy_on.ne.0) then
                  if(ncall .ne. 0 .or. sputter_namelist_used) then
                     t2=(sputter_yield2(ix,iy,is,1)+sput_chem_alpha*t2)/
     *                    (1+sput_chem_alpha)
                  endif
                  sputter_yield2(ix,iy,is,1)=t2
                  tchee(is)=tchee(is)+t2
                  shi0(ix,iy,0) = shi0(ix,iy,0)+t2
                  if (is2.ne.-1)
     >             shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               else
                  sputter_yield2(ix,iy,is,1)=0.0_R8
               endif
            else
               sputter_yield(ix,iy,is,1)=0.0_R8
               sputter_yield2(ix,iy,is,1)=0.0_R8
               chemical_sputtering(iwall,is) = 0.0_R8
            endif
c physical sputtering
            if(phys_sput(is,istra).gt.0.0_R8) then
cxpb We first treat the bulk material
               if(sput_dst_bulk.ne.-1) then
                  itrack=track_index(sput_dst_bulk)
                  redep_max=
     .             max(redep_max,f_redep(ixn,iyn,sput_dst_bulk))
                  yself=sputter(sput_dst_bulk+1,zaf(sput_dst_bulk+1),
     1                          sput_dst_bulk,
     .                  (pof-Phi_app)*zaf(sput_dst_bulk+1)*ev)*
     2                   phys_sput(sput_dst_bulk+1,istra)
                  if(f_redep(ixn,iyn,sput_dst_bulk).gt.0.0_R8 .and.
     1               yself.gt.0.0_R8 .and.
     2              (output.gt.2 .or. (main_call.and.output.ge.2) .or.
     3               yself.ge.1.0_R8))
     >             write(*,'(a,1x,a2,2i4,1p,2g12.5)')
     .              'Self-sputtering coefficient ',
     .               is_codes(sput_dst_bulk),
     .               ix, iy, pof-Phi_app, yself
                  yself_max=max(yself_max,yself)
                  if(sput_frac_flag.eq.0) then
                    ytmp=sputter(is,zaf(is),sput_dst_bulk,Ei*ev) *
     1               phys_sput(is,istra) *
     2               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))/
     3               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself)
! the two lines above are the reduction of sputtering due to prompt redeposition
! and then the multiplication due to self-sputtering of the redeposited particles
                  else
                    ytmp=sputter(is,zaf(is),sput_dst_bulk,Ei*ev) *
     1               phys_sput(is,istra) *
     2               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))/
     3               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself*
     4                sput_frac(0))
! we only multiply the self-sputtering cascade contribution by sput_frac(0) since the
! primary sputtered particles will be multiplied by this factor just below
                  endif
                  t1=t0 * ytmp
                  if(ncall.ne.0 .or. sputter_namelist_used) then
                     t1=(sputter_yield(ix,iy,is,2)+sput_phys_alpha*t1)/
     *                 (1.0_R8+sput_phys_alpha)
                  endif
                  if(sput_frac_flag.ne.0) then
                     ytmp=ytmp*sput_frac(0)
                     t1=t1*sput_frac(0)
                  endif
                  if(output.ge.11) then
                    write(*,'(a,3(1x,i3),1x,1p,3g15.8)')
     &               'b2stbr_phys.recycle: phys_sput,redep,self_sput ',
     &                 0,is,sput_dst_bulk,ytmp,
     &                 f_redep(ixn,iyn,sput_dst_bulk),yself
                  endif
                  physical_sputtering(iwall,is) = ytmp
                  sputter_yield(ix,iy,is,2) = t1
                  tphys(is) = tphys(is) + t1
                  sna0(ix,iy,0,sput_dst_bulk) =
     &              sna0(ix,iy,0,sput_dst_bulk) + t1
                  sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra)) =
     &             sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra))
     &             + t1
                  if (output.gt.2) write(*,*)
     &             'Physical sputtering from bulk: ',
     &             'Contribution to sna_stbr_sput. is = ',sput_dst_bulk,
     &             ' istra = ',istra,' ix = ',ix,' iy = ',iy

cxpb  Erosion from the bulk material if physically sputtered
                  if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call) then
cxpb  We count here as eroded all the sputtered material, including the promptly redeposited fraction
cxpb  We add this promptly redeposited fraction in the deposited material
                    if(sput_frac_flag.eq.0) then
                      erosion(iwall,itrack) = erosion(iwall,itrack) +
     1                 dtim*plate_time_factor(iwall)*t1/
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))*
     3                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself)
                      deposition(iwall,itrack) =
     .                 deposition(iwall,itrack) +
     1                 dtim*plate_time_factor(iwall)*t1*(1.0_R8-1.0_R8/
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst))*
     3                 (1.0_R8-f_redep(ixn,iyn,sput_dst)*yself))
                    else
                      erosion(iwall,itrack) = erosion(iwall,itrack) +
     1                 dtim*plate_time_factor(iwall)*t1/
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))*
     3                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself*
     4                  sput_frac(0))
                      deposition(iwall,itrack) =
     .                 deposition(iwall,itrack) +
     1                 dtim*plate_time_factor(iwall)*t1*(1.0_R8-1.0_R8/
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst))*
     3                 (1.0_R8-f_redep(ixn,iyn,sput_dst)*yself*
     4                  sput_frac(0)))
                    endif
                 endif
               else
                 physical_sputtering(iwall,is) = 0.0_R8
                 sputter_yield(ix,iy,is,2) = 0.0_R8
               endif
               if(sput_frac_flag.ne.0) then
                  do itrack=1,ntargsp
                    redep_max=max(redep_max,
     1               f_redep(ixn,iyn,target_intcode(itrack)))
                    yself=sputter(target_intcode(itrack)+1,
     .                            zaf(target_intcode(itrack)+1),
     1                            target_intcode(itrack),
     .                   (pof-Phi_app)*
     .                    zaf(target_intcode(itrack)+1)*ev)*
     2                   phys_sput(target_intcode(itrack)+1,istra)
                    if(f_redep(ixn,iyn,target_intcode(itrack)).gt.0.0_R8
     1                 .and. yself.gt.0.0_R8 .and.
     2                (output.gt.2 .or. (main_call.and.output.ge.2) .or.
     3                 yself.ge.1.0_R8))
     >               write(*,'(a,1x,i3,2i4,1p,2g12.5)')
     .                'Self-sputtering coefficient ',
     .                 target_intcode(itrack),
     .                 ix, iy, pof-Phi_app, yself
                    yself_max=max(yself_max,yself)
                    ytmp=
     .               sputter(is,zaf(is),target_intcode(itrack),Ei*ev)*
     1                   phys_sput(is,istra) * sput_frac(itrack) *
     2                   (1.0_R8-
     3                    f_redep(ixn,iyn,target_intcode(itrack)))/
     4                   (1.0_R8-
     5                    f_redep(ixn,iyn,target_intcode(itrack))*
     6                    yself*sput_frac(itrack))
                    if(output.ge.11) then
                      write(*,'(a,3(1x,i3),1x,1p,3g15.8)')
     &                'b2stbr_phys.recycle: phys_sput,redep,self_sput ',
     &                       itrack,is,target_intcode(itrack),ytmp,
     &                       f_redep(ixn,iyn,target_intcode(itrack)),
     &                       yself
                    endif
                    t1=t0 * ytmp
                    physical_sputtering(iwall,is) =
     &                    physical_sputtering(iwall,is) + ytmp
                    sputter_yield(ix,iy,is,2) =
     &                    sputter_yield(ix,iy,is,2) + t1
                    tphys(is) = tphys(is) + t1
                    sna0(ix,iy,0,target_intcode(itrack)) =
     &                    sna0(ix,iy,0,target_intcode(itrack)) + t1
                    sna_stbr_sput(b2espcr(target_intcode(itrack)),
     &                   lbcstra(istra)) = t1 +
     &               sna_stbr_sput(b2espcr(target_intcode(itrack)),
     &                   lbcstra(istra))
!                 write(*,*)
!     &     'Erosion. Contribution to sna_stbr_sput. is = ',
!     $    target_intcode(itrack),' istra = ',istra,' ix = ',
!     &      ix,' iy = ',iy

cxpb  Removal from the deposited material if physically sputtered
cxpb  with mixed materials model
                    if(iwall.gt.0.and.main_call)
     1                    deposition(iwall,itrack) =
     &                    deposition(iwall,itrack) -
     1                    dtim*plate_time_factor(iwall)*t1
                  enddo
               endif
               if(sputter_energy_on.ne.0) then
                  if(sput_dst_bulk.ne.-1) then
                     eself=sputter2(sput_dst_bulk+1,sput_dst_bulk,
     1                (pof-Phi_app)*zaf(sput_dst_bulk+1)*ev)*
     2                phys_sput(sput_dst_bulk+1,istra)
                     yself=sputter(sput_dst_bulk+1,zaf(sput_dst_bulk+1),
     1                          sput_dst_bulk,
     .                  (pof-Phi_app)*zaf(sput_dst_bulk+1)*ev)*
     2                   phys_sput(sput_dst_bulk+1,istra)
                     ytmp=sputter(is,zaf(is),sput_dst_bulk,Ei*ev) *
     1                 phys_sput(is,istra)
                     ztmp=sputter2(is,sput_dst_bulk,Ei*ev) *
     1                 phys_sput(is,istra) *
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))
                     if(sput_frac_flag.ne.0) then
                       ytmp=ytmp*sput_frac(0)
                       ztmp=ztmp*sput_frac(0)
                       yself=yself*sput_frac(0)
                     endif
                     t2=t0 * ztmp * Ei * ev
cxpb  We now must remove the energy lost in ionizing the promptly redeposited material
cxpb  but from the electrons!
                     t3=-t0*ev*rpi(ixn,iyn,sput_dst_bulk+1)*
     1                 ytmp*f_redep(ixn,iyn,sput_dst_bulk)/
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself)
cxpb  We must now add the energy brought back in by the self-sputtered material
                     t2=t2+t0*f_redep(ixn,iyn,sput_dst_bulk)*
     1                 ytmp*eself*(pof-Phi_app)*
     .                 zaf(sput_dst_bulk+1)*ev/
     2                 (1.0_R8-
     3                  f_redep(ixn,iyn,sput_dst_bulk)*yself)
                     physical_sputtering_energy(iwall,is) = ztmp
                     if(ncall .ne. 0 .or. sputter_namelist_used) then
                        t2=(sputter_yield2(ix,iy,is,2)+
     *                   sput_phys_alpha*t2)/(1.0_R8+sput_phys_alpha)
                     endif
                     if(output.ge.11) then
                        write(*,'(a,3(1x,i3),1p,5g15.8)')
     &             'b2stbr_phys.recycle: phys_sput2,redep,e_self_sput ',
     &                       0,is,sput_dst_bulk,ztmp,
     &                       f_redep(ixn,iyn,sput_dst_bulk),
     &                       eself,t2,t3
                     endif
                     sputter_yield2(ix,iy,is,2) = t2
                     tphye(is) = tphye(is) + t2 + t3
                     shi0(ix,iy,0) = shi0(ix,iy,0) + t2
                     she0(ix,iy,0) = she0(ix,iy,0) + t3
                     if (is2.ne.-1)
     >                shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
                  else
                     physical_sputtering_energy(iwall,is) = 0.0_R8
                     sputter_yield2(ix,iy,is,2) = 0.0_R8
                  endif
                  if(sput_frac_flag.ne.0) then
                     do itrack=1,ntargsp
                        eself=sputter2(target_intcode(itrack)+1,
     .                   target_intcode(itrack),
     1                   (pof-Phi_app)*
     .                    zaf(target_intcode(itrack)+1)*ev)*
     2                   phys_sput(target_intcode(itrack)+1,istra)
                        yself=sputter(target_intcode(itrack)+1,
     .                            zaf(target_intcode(itrack)+1),
     1                            target_intcode(itrack),
     .                            (pof-Phi_app)*
     .                            zaf(target_intcode(itrack)+1)*ev)*
     2                   phys_sput(target_intcode(itrack)+1,istra)*
     3                   sput_frac(itrack)
                        ytmp=sputter(is,zaf(is),
     1                               target_intcode(itrack),Ei*ev) *
     1                       phys_sput(is,istra) * sput_frac(itrack)
                        ztmp=sputter2(is,target_intcode(itrack),Ei*ev) *
     1                       phys_sput(is,istra) * sput_frac(itrack) *
     2                       (1.0_R8-
     3                        f_redep(ixn,iyn,target_intcode(itrack)))
                        t2=t0 * ztmp * Ei * ev
cxpb  We now must remove the energy lost in ionizing the promptly redeposited material
cxpb  but from the electrons
                        t3=-t0*ev*rpi(ixn,iyn,target_intcode(itrack)+1)*
     1                   ytmp*f_redep(ixn,iyn,target_intcode(itrack))/
     2                   (1.0_R8-
     3                    f_redep(ixn,iyn,target_intcode(itrack))*yself)
cxpb  We must now add the energy brought back in by the self-sputtered material
                        t2=t2+
     .                   t0*f_redep(ixn,iyn,target_intcode(itrack))*
     1                   ytmp*eself*(pof-Phi_app)*
     .                   zaf(target_intcode(itrack)+1)*ev/
     2                   (1.0_R8-
     3                    f_redep(ixn,iyn,target_intcode(itrack))*yself)
                        physical_sputtering_energy(iwall,is) =
     &                       physical_sputtering_energy(iwall,is) + ztmp
                        if(output.ge.11) then
                       write(*,'(a,3(1x,i3),1p,3g15.8)')
     &             'b2stbr_phys.recycle: phys_sput2,redep,e_self_sput ',
     &                       itrack,is,target_intcode(itrack),ztmp,
     &                       f_redep(ixn,iyn,target_intcode(itrack)),
     &                       eself
                        endif
                        sputter_yield2(ix,iy,is,2) =
     &                       sputter_yield2(ix,iy,is,2) + t2
                        tphye(is) = tphye(is) + t2 + t3
                        shi0(ix,iy,0) = shi0(ix,iy,0) + t2
                        she0(ix,iy,0) = she0(ix,iy,0) + t3
                        if (is2.ne.-1)
     >                   shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
                     enddo
                  endif
               else
                  physical_sputtering_energy(iwall,is) = 0.0_R8
                  sputter_yield2(ix,iy,is,2) = 0.0_R8
               endif
            else
               physical_sputtering(iwall,is) = 0.0_R8
               physical_sputtering_energy(iwall,is) = 0.0_R8
               sputter_yield(ix,iy,is,2) = 0.0_R8
               sputter_yield2(ix,iy,is,2) = 0.0_R8
            endif
c thermal evaporation     ! no fix for mixed materials yet
            itrack=track_index(is)
            if(therm_evap.gt.0.0_R8.and.iwall.gt.0.and.
     .        (is.eq.sput_dst_bulk .or.
     .        (is_neutral(is).and.(itrack.gt.0.and.itrack.le.ntrack)))
     >        ) then
               call thermev(iwall, is, surface_temp, GammaEvap, t2)
               thermal_evaporation(iwall,is) = GammaEvap
               t1=GammaEvap*area*therm_evap
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
               sna_stbr_sput(b2espcr(is),lbcstra(istra)) =   ! IYS 18.09.2018
     &          sna_stbr_sput(b2espcr(is),lbcstra(istra)) + t1
               if (output.gt.2) write(*,*)
     &          'Thermal evaporation : ',
     &          'Contribution to sna_stbr_sput. is = ',is,
     &          ' istra = ',istra,' ix = ',ix,' iy = ',iy
cxpb  Erosion from the bulk material if thermally evaporated
               if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call)
     1              erosion(iwall,itrack) = erosion(iwall,itrack) +
     1              dtim*plate_time_factor(iwall)*t1
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               thevp(is)=thevp(is)+t1
               thvpe(is)=thvpe(is)+t2
            else
               thermal_evaporation(iwall,is) = 0.0_R8
            endif
c radiation enhanced sublimation
            if(sput_res.gt.0.0_R8.and.
     &        (l_chem_sput.or.sput_frac_flag.ne.0)) then
               if(zbulk.eq.6) then
                  call ressput(Ei,zaf(is),zn(is),zn(sput_dst_bulk),
     &              am(is),am(sput_dst_bulk),surface_temp,
     &              Pflx/neutral_sources_rescale,sput_dst_bulk,yldres)
                  RES_sputtering(iwall,is) = yldres
                  t1=yldres*neutral_sources_rescale*
     &              area*sput_res
               endif
               if(sput_frac_flag.ne.0) then
                  if(zbulk.ne.6) then
                    call ressput(Ei,zaf(is),zn(is),6.0_R8,
     &               am(is),12.011_R8,surface_temp,
     &               Pflx/neutral_sources_rescale,sput_dst_bulk,yldres)
                    RES_sputtering(iwall,is) = yldres
                    t1=yldres*neutral_sources_rescale*
     &               area*sput_res
                  endif
                  t1=t1*sput_frac_chem_sum
                  RES_sputtering(iwall,is) =
     &                 RES_sputtering(iwall,is)*sput_frac_chem_sum
               endif
               if(sput_dst_bulk.ne.-1) then
                  sna0(ix,iy,0,sput_dst_bulk) =
     &              sna0(ix,iy,0,sput_dst_bulk) + t1
                  sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra)) =   ! IYS 18.09.2018
     &              sna_stbr_sput(b2espcr(sput_dst_bulk),lbcstra(istra))
     &              + t1
                  if (output.gt.2) write(*,*)
     &             'RES sputtering: ',
     &             'Contribution to sna_stbr_sput. is = ',sput_dst_bulk,
     &             ' istra = ',istra,' ix = ',ix,' iy = ',iy
               endif

cxpb  Erosion from the bulk material if RES sputtered
cxpb  without mixed materials model
               if(sput_frac_flag.eq.0) then
                  if (sput_dst_bulk.eq.-1) then
                    itrack=0
                  else
                    itrack=track_index(sput_dst_bulk)
                  endif
                  if(itrack.gt.0.and.itrack.le.ntrack.and.
     1               iwall.gt.0.and.main_call)
     1               erosion(iwall,itrack) = erosion(iwall,itrack) +
     1               dtim*plate_time_factor(iwall)*t1
               else
                  if(iwall.gt.0.and.main_call
     &                 .and.sput_frac_chem_sum.gt.0) then
                     if(sput_dst_bulk.eq.-1) then
                       itrack=0
                     else
                       itrack=track_index(sput_dst_bulk)
                     endif
cxpb  Erosion from the bulk material if RES sputtered
cxpb  with mixed materials model
                     if(itrack.gt.0.and.itrack.le.ntrack)
     &                    erosion(iwall,itrack) =
     &                    erosion(iwall,itrack) +
     1                    dtim*plate_time_factor(iwall)*t1*
     &                    sput_frac_chem(0)/sput_frac_chem_sum
cxpb  Removal from the deposited material if RES sputtered
cxpb  with mixed materials model
                     do itrack=1,ntargsp
                        if(itrack.gt.0.and.itrack.le.ntrack)
     &                    deposition(iwall,itrack) =
     &                    deposition(iwall,itrack) -
     1                    dtim*plate_time_factor(iwall)*t1*
     &                    sput_frac_chem(itrack)/sput_frac_chem_sum
                     enddo
                  endif
               endif
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               tresn(is)=tresn(is)+t1
               trese(is)=trese(is)+t2
            else
               RES_sputtering(iwall,is) = 0.0_R8
            endif
c back scattering ! mixed materials included thru layer constituents definitions
            if (recyc(is,istra).eq.0.0_R8 .and.
     &         erecyc(is,istra).eq.0.0_R8 .and. .not.is_neutral(is)
     &   .and. reflection_on.ne.0) then ! use reflection model
               local_reflection = reflect(is,Ei*ev,iwall)
               reflection_max = max(reflection_max, local_reflection)
               t1=t0 * local_reflection
               t2=t0 * reflect2(is,Ei*ev,iwall)*Ei*ev
               itrack=track_index(is0)
cxpb  Removal from the deposited material if backscattered
cxpb  because incoming non-recycling flux was originally counted
cxpb  as deposition
               if(itrack.gt.0.and.itrack.le.ntrack
     1              .and.iwall.gt.0.and.main_call) then
                  deposition(iwall,itrack)=deposition(iwall,itrack)-
     1                 dtim*plate_time_factor(iwall)*t1
               endif
               backscattering(iwall,is) = reflect(is,Ei*ev,iwall)
               backscattering_energy(iwall,is) =
     &              reflect2(is,Ei*ev,iwall)
               if(reflection_on.eq.1) then
                 sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
                 if (is_neutral(is)) then                           ! IYS 18.09.2018
                   sna_stbr_neut(b2espcr(is),lbcstra(istra)) =
     &              sna_stbr_neut(b2espcr(is),lbcstra(istra)) + t1
                   if (output.gt.2) write(*,*) 'Reflection:',
     >              ' Contribution to sna_stbr_neut. is = ',is,
     &              ' istra = ',istra,' ix = ',ix,' iy = ',iy
                 end if
               else
                 sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0) + t1    ! reflection to the neutral
                 if (is_neutral(is0)) then                          ! IYS 18.09.2018
                   sna_stbr_neut(b2espcr(is0),lbcstra(istra)) =
     &              sna_stbr_neut(b2espcr(is0),lbcstra(istra)) + t1
                   if (output.gt.2)
     >              write(*,*) 'Reflection into neutral:',
     &               'Contribution to sna_stbr_neut. is0 = ',is0,
     &               ' istra = ',istra,' ix = ',ix,' iy = ',iy
                 endif
               endif
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               trfln(is)=trfln(is)+t1
               trfle(is)=trfle(is)+t2
            else
               backscattering(iwall,is) = 0.0_R8
               backscattering_energy(iwall,is) = 0.0_R8
            endif
         endif
*    ..recycling into next ionized species
         if(is.lt.ns-1) then
           if(lnext(is,is+1)) then
              sna0(ix,iy,0,is+1) = sna0(ix,iy,0,is+1)+
     2        rcion(is,istra)*t0/neutral_sources_rescale
              if (is_neutral(is+1)) then           ! IYS 18.09.2018
                sna_stbr_neut(b2espcr(is+1),lbcstra(istra)) =
     &              sna_stbr_neut(b2espcr(is+1),lbcstra(istra)) +
     2              rcion(is,istra)*t0/neutral_sources_rescale
                if (output.gt.2) write(*,*) 'Recycling into ions: ',
     &           'Contribution to sna_stbr_neut. is+1 = ',is+1,
     &           ' istra = ',istra,' ix = ',ix,' iy = ',iy
              endif
           endif
         endif
      enddo      ! do is = 0, ns-1
*    ..recycling at core into highest ionized stage
      if (lcore) then
*    ..loop over species
         is0 = ns-1
         is1 = ns-1
         do is = nsmax, nsmin, -1
            if (zn(is0).ne.zn(is).or.
     &           amtol.lt.abs(am(is0)-am(is))) then
               is0 = is
            endif
cxpb
            if(zamax(is).eq.zn(is)) is1=is
            if(zn(is1).ne.zn(is).or.amtol.lt.abs(am(is1)-am(is))) then
              write(*,'(a,3i3)')
     .         'b2stbr problem: is, is0, is1 ',is,is0,is1
            endif
            if(is0.ne.is1) then
               if(ncall.eq.0)
     1          write(*,'(a,2i3)') 'b2stbr: setting is0 to is1 ',is0,is1
               is0=is1
            endif
cxpb
*          (species (is) recycles into is0)
            t0 = max(0.0_R8,isign*fna(ixf,iyf,idir,is)) *
     .           core_sources_rescale ! xpb: 2002.02.01
            if (is.ne.is0) then
               sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0)+
     &          recyc(is,istra)*t0
               if (is_neutral(is0)) then          ! IYS 18.09.2018
                 sna_stbr_neut(b2espcr(is0),lbcstra(istra)) =
     &            sna_stbr_neut(b2espcr(is0),lbcstra(istra)) +
     &            recyc(is,istra)*t0
                 if (output.gt.2) write(*,*) 'Core recycling: ',
     &            'Contribution to sna_stbr_neut. is0 = ',is0,
     &            ' istra = ',istra,' ix = ',ix,' iy = ',iy
               end if

cwdk  recycling of true parallel momentum
c     (at core, cell faces will be aligned with flux surfaces)
               smo0(ix,iy,0,is0) = smo0(ix,iy,0,is0)+
     .           mrecyc(is,istra)*am(is)*mp*ua(ix,iy,is)*
     .           recyc(is,istra)*t0*hz(ix,iy)
               shi0(ix,iy,0) = shi0(ix,iy,0)+erecyc(is,istra)*
     .          (1.5_R8+BoRiS)*(zaf(is)*tef+tif)*recyc(is,istra)*t0
     &          *recycled_neutrals_contr  !srv
            endif
         enddo      ! do is = nsmax, nsmin, -1
      endif

      return
      end subroutine recycle

*-----------------------------------------------------------------------
*.end b2stbr_phys

      end subroutine b2stbr_phys

