*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr (nx, ny, ns, nscx, nscxmax, iscx, ismain,  
     &  dtim, BoRiS, 
     &  vol, hx, hy, gs, qz, pbs, na, ua, kinrgy,
     &  te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni, rsa, bb,
     &  sna0, smo0, she0, shi0, sch0, sne0, wrk1, main_call)
      use b2mod_types
      use b2mod_wall
      use b2mod_diag
      use b2mod_time
      use b2mod_layer
      use b2mod_elements
      use b2mod_sputter
      use b2mod_first_flight
      use b2mod_constants
      use b2mod_indirect
      use b2mod_neutrals_namelist
      use b2mod_sources
     1 , only : b2stbr_sna,b2stbr_smo,b2stbr_she,b2stbr_shi,b2stbr_sch,
     2          b2stbr_sne
#ifdef B25_EIRENE
      use eirmod_extrab25
      use b2mod_work
     1 , only : sna0_eir_h, smo0_eir_h, she0_eir_h,       !IYS 28.03.12 {
     &          shi0_eir_h, rfluxa_h, rfluxm_h,
csw 25feb2014 fix IVS's version of source term averaging
     &          alloc_hist_arrays                         !IYS 28.03.12 } 
#endif
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#ifndef NO_CDF
      use b2mod_tallies
#endif
#ifdef ASTRA
      use b2mod_b2_to_astra
     , , only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 20.02.12 11.06.13
#endif	                       
#ifndef ALLOCATE_AND_NAMELIST
#     include "DIMENSIONS.F"
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, nscxmax, iscx(0:nscxmax-1), ismain
      real (kind=R8) :: dtim, BoRiS,
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  rsa(-1:nx,-1:ny,0:ns-1), bb(-1:nx,-1:ny,0:3)
      logical main_call
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk1(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBR provides linearized expressions for the sources of
*     density, parallel momentum and heat, due to special boundary
*     conditions and the neutral gas recycling model.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is, js, ib, ic, ibnd, itrack, k, ia, im, iih
      integer, save :: jxi, jxa, jsep
      integer, save :: ihist = 0                                     !IYS 29.03.2012
      integer, save :: no_b2stbr_phys_call=0
      integer, save :: out_tfiles = 0                                !srv 07.04.09
      integer, save :: iout = 0,        iout_b2npmo = 0,             !srv 11.09.09 {
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0,
     &                 ncall_iout = -1,                              !srv 23.03.12
     &                 nhist = 1                                     !IYS 28.03.12
      real (kind=R8), save :: recycled_neutrals_contr = 1.0_R8       !srv 12.04.08
      character chns*3, chk*1, ch3*3                                 !srv 11.09.09 } 23.03.12
      integer nx_tmp, ny_tmp, ns_tmp
      real (kind=R8) ::
     *  lnorm, tchem(0:ns-1), tchee(0:ns-1), tphys(0:ns-1),
     *  tphye(0:ns-1), thevp(0:ns-1), thvpe(0:ns-1), tresn(0:ns-1),
     *  trese(0:ns-1), trfln(0:ns-1), trfle(0:ns-1), t_frac, t0
      real (kind=R8) ::                                              !srv 29.04.10 {
     &  sna0_eir(-1:nx,-1:ny,0:1,0:ns-1),
     &  smo0_eir(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0_eir(-1:nx,-1:ny,0:3), shi0_eir(-1:nx,-1:ny,0:3),         
     &  fatm_eir(-1:nx,-1:ny,0:1,0:ns-1), datm_eir(-1:nx,-1:ny,0:ns-1) !srv 29.04.10 }
      integer ncall, sput_src, first_flight, 
     *  in_no_of_angles, in_size_of_table, in_no_of_start_points,
     1  use_eirene, output, coreregno, coreregn2, neutrals_namelist,
     2  sput_chem_model, boundary_namelist,
     3  reflection_on, sputter_energy_on, 
     4  b2wall_netcdf, potential_at_guard_cell,
     5  sput_frac_flag
      real (kind=R8) ::
     *  sput_phys, sput_frc, sput_chem_alpha, 
     1  sput_phys_alpha, in_dl, alpha, therm_evap, sput_res,
     3  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     4  sput_mixed_alpha, sput_mixed_beta, redep_alpha
      real (kind=R8) ::
     *  start_time,end_time,pfrpuff,mainpuff,corepuff,
     *  sum_start,sum_end,neutral_rescale,neutral_sources_rescale,
     *  sum_strti,sum_ndi,sum_nde,core_sources_rescale,
     *  shi0_ff(-1:nx,-1:ny,0:nscx-1),
     *  work(-1:nx,-1:ny,0:1), wrk2(-1:nx,-1:ny), wrk3(-1:nx,-1:ny),
     *  f_redep(-1:nx,-1:ny,0:ns-1)
      real (kind=R8), allocatable, save ::
     *  old_erosion(:,:), old_deposition(:,:)
      real (kind=R8) :: amtol
      parameter (amtol=1.0e-2_R8)
      logical catch_up, new_sputter_namelist, wall_namelist_used
      logical equilibrium(ntrack), over_erosion
#ifdef B25_EIRENE
      integer ih_end
      integer, save :: neutr_avg = 0
      logical, save :: full_history = .false.
#endif
#ifndef NO_CDF
#include <netcdf.inc>
      integer :: status, ncid
      integer, save :: ntstep_b2wall
#endif
*   ..procedures
      logical lnext
      external subini, subend, xertst, xerrab, sfill, ipgetr, ipgeti
      external b2xvsg, b2xvff, x05baf, smax, smin, get_jsep
      real (kind=R8) ::
     *  x05baf, smax, smin
      intrinsic abs, mod, max
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol
*   ..initialisation
      save ncall, sput_src, sput_frc, sput_phys, 
     *     sput_chem_alpha, sput_phys_alpha, first_flight, in_dl, 
     *     in_no_of_angles, in_size_of_table, in_no_of_start_points, 
     1     use_eirene, neutral_rescale, neutral_sources_rescale, output,
     2     coreregno, coreregn2, core_sources_rescale, alpha,
     3     neutrals_namelist, therm_evap, sput_res, 
     5     sput_chem_model, boundary_namelist, 
     7     reflection_on, sputter_energy_on, 
     8     b2wall_netcdf, sput_frac_flag,
     9     potential_at_guard_cell,
     1     sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     4     sput_mixed_alpha, sput_mixed_beta, redep_alpha
      data ncall/0/, sput_src/1/, sput_frc/0.0_R8/,
     1     sput_phys/0.0_R8/, sput_chem_alpha/1.0e-2_R8/, 
     *     sput_phys_alpha/1.0e-2_R8/, first_flight/0/,
     *     in_dl/0.001_R8/, in_no_of_angles/9/, alpha/0.25_R8/,
     *     in_size_of_table/200000/, 
     1     use_eirene/0/, neutral_rescale/1.0_R8/, 
     *     neutral_sources_rescale/1.0_R8/, output/0/,
     1     coreregno/1/, coreregn2/4/, core_sources_rescale/1.0_R8/,
     2     neutrals_namelist/0/, therm_evap/0.0_R8/, sput_res/0.0_R8/,
     7     sput_chem_model/0/, boundary_namelist/0/,
     8     reflection_on/1/, sputter_energy_on/1/, b2wall_netcdf/0/, 
     2     sput_frac_flag/0/, sput_chem_cutoff_alpha/1.0_R8/, 
     3     sput_chem_cutoff_beta/3.0_R8/,
     4     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/,
     5     redep_alpha/0.0_R8/, potential_at_guard_cell/1/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Note: This version of b2stbr is for the case when recycling is
*     calculated via the fluid model.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbr')
      new_sputter_namelist=.false.
      if (ncall.eq.0) then
        call get_jsep(nx,ny,jxi,jxa,jsep)
        call ipgeti ('b2stbr_output', output)
        call ipgeti ('b2stbr_out_tfiles', out_tfiles)                !srv 07.04.09
        call ipgeti ('b2stbr_iout', iout)                            !srv 11.09.09 {
        call ipgeti ('b2npc9_iout', iout_b2npc9)
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 11.09.09 }
        call ipgeti ('b2stbr_ncall_iout', ncall_iout)                !srv 23.03.12
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2stbr_eir_src_nhist', nhist)                  !IYS 28.03.12
#ifdef B25_EIRENE
        call ipgeti ('eirene_neutr_avg', neutr_avg)
        call xertst (.not.(neutr_avg.gt.1 .and. nhist.gt.1),
     .   'Cannot have both Eirene averaging schemes simultaneously!')
        if (nhist.gt.1) call alloc_hist_arrays(nx,ny,ns,nhist)
#endif
        call ipgetr ('b2stbr_recycled_neutrals_contr',               !srv 12.04.08 
     &                      recycled_neutrals_contr)
        call xertst (0.0_R8.le.recycled_neutrals_contr, 
     &      'faulty internal parameter recycled_neutrals_contr')
        call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
        call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
        call ipgeti ('b2stbr_first_flight', first_flight)
        if(first_flight.ne.0) then
           call ipgetr ('b2stbr_first_flight_dl', in_dl)
           call ipgeti ('b2stbr_first_flight_no_of_flights', 
     *          in_no_of_angles)
           call ipgeti ('b2stbr_first_flight_table_size', 
     *          in_size_of_table)
           in_no_of_start_points = 2*(nx+2)+2*max(nncut,1)*(ny+2)
           call ipgeti ('b2stbr_first_flight_no_of_start_points',
     *       in_no_of_start_points)
           call xertst (0.lt.in_no_of_start_points,
     *      'faulty input first_flight_no_of_start_points')
           start_time=x05baf()
           call alloc_b2mod_first_flight(
     *          in_no_of_angles,in_size_of_table,
     *          in_no_of_start_points,nx,ny,in_dl)
           do ix=-1,nx
             do ib=-1,ny
               if(leftix(ix,ib).eq.-2 .and. region(ix,ib,0).ne.0)
     >          call setup_b2mod_first_flight(ix,ib,nx,ny,'W',output)
               if(rightix(ix,ib).eq.nx+1 .and. region(ix,ib,0).ne.0)
     >          call setup_b2mod_first_flight(ix,ib,nx,ny,'E',output)
             enddo
           enddo
           do ib=-1,nx
             do iy=-1,ny
               if(topiy(ib,iy).eq.ny+1 .and. region(ib,iy,0).ne.0)
     >          call setup_b2mod_first_flight(ib,iy,nx,ny,'N',output)
               if(bottomiy(ib,iy).eq.-2 .and. region(ib,iy,0).ne.0)
     >          call setup_b2mod_first_flight(ib,iy,nx,ny,'S',output)
             enddo
           enddo
           end_time=x05baf()
           write(*,'(a,i8,a,i8,a,g10.3,a,g12.5,a)') 
     *          'Used ',ptr_table,' out of ',size_of_table,
     *          ' (',real(ptr_table)/real(size_of_table)*100,'%)',
     *          end_time-start_time,' cpu-seconds'
        endif
        call ipgeti ('b2stbr_sput_src', sput_src)
        call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
        call ipgeti ('b2stbr_sput_frac_flag', sput_frac_flag)
        if(sput_frac_flag.ne.0) then
          call ipgetr ('b2stbr_sput_mixed_alpha', sput_mixed_alpha)
          call xertst (sput_mixed_alpha.ge.0.0_R8,
     .     'faulty parameter sput_mixed_alpha')
          call ipgetr ('b2stbr_sput_mixed_beta', sput_mixed_beta)
          call xertst (sput_mixed_beta.gt.0.0_R8,
     .     'faulty parameter sput_mixed_beta')
        endif
        call ipgetr ('b2stbr_sput_frc', sput_frc)
        if(sput_frc.ne.0.0_R8) then
          call ipgetr ('b2stbr_sput_chem_alpha', sput_chem_alpha)
          call xertst (0.0_R8.le.sput_chem_alpha,
     1     'faulty b2stbr_sput_chem_alpha')
          call ipgetr ('b2stbr_sput_chem_cutoff_alpha', 
     1     sput_chem_cutoff_alpha)
          call xertst (0.0_R8.le.sput_chem_cutoff_alpha,
     1     'faulty b2stbr_sput_chem_cutoff_alpha')
          call ipgetr ('b2stbr_sput_chem_cutoff_beta', 
     1     sput_chem_cutoff_beta)
          call xertst (0.0_R8.le.sput_chem_cutoff_beta,
     1     'faulty b2stbr_sput_chem_cutoff_beta')
        endif
        call ipgetr ('b2stbr_sput_phys', sput_phys)
        call xertst (0.0_R8.le.sput_phys,
     1   'faulty input b2stbr_sput_phys')
        if(sput_phys.ne.0.0_R8) then
          call xertst (sput_dst.lt.ns,
     1     'faulty b2stbr_sput_dst')
          call ipgetr ('b2stbr_sput_phys_alpha', sput_phys_alpha)
          call xertst (0.0_R8.le.sput_phys_alpha,
     1     'faulty b2stbr_sput_phys_alpha')
        endif
        call ipgeti ('b2stbr_b2wall_netcdf', b2wall_netcdf)
        call ipgetr ('b2stbr_therm_evap', therm_evap)
        call xertst (0.0_R8.le.therm_evap, 'faulty b2stbr_therm_evap')
        call ipgetr ('b2stbr_sput_res', sput_res)
        call xertst (0.0_R8.le.sput_res, 'faulty b2stbr_sput_res')
        call ipgetr ('b2stbr_alpha', alpha)
        call xertst (0.0_R8.le.alpha, 'faulty b2stbr_alpha')
        call ipgeti ('b2stbr_reflection_on', reflection_on)
        call ipgeti ('b2stbr_sputter_energy_on', sputter_energy_on)
        call ipgetr ('b2mndr_rescale_neutrals', neutral_rescale)
        call ipgetr ('b2mndr_rescale_neutrals_sources',
     *       neutral_sources_rescale)
        if(use_eirene.eq.0) neutral_rescale=1.0_R8
        if(use_eirene.eq.0) neutral_sources_rescale=1.0_R8
        call ipgeti ('b2stbc_coreregno', coreregno)
        call xertst (cbirso.le.coreregno.and.
     1   coreregno.le.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregno')
        if (nncut.ge.2) then
          call ipgeti ('b2stbc_coreregn2', coreregn2)
          call xertst (cbirso.le.coreregn2.and.
     1     coreregn2.le.cbirso+cbnrso-1,
     2     'faulty internal parameter coreregn2')
        endif
        call ipgetr ('b2stbr_core_sources_rescale', 
     1   core_sources_rescale)
        call xertst (0.0_R8.le.core_sources_rescale,
     1   'faulty internal parameter core_sources_rescale')
        call ipgeti ('b2stbr_no_b2stbr_phys_call', no_b2stbr_phys_call)
        call ipgetr ('b2stbr_redep_alpha', redep_alpha)
        call xertst (0.0_R8.le.redep_alpha,
     1   'faulty internal parameter redep_alpha')
        call ipgeti ('b2stbr_potential_at_guard_cell', 
     &                potential_at_guard_cell)
        allocate(
     *   old_erosion(nwall,ntrack), old_deposition(nwall,ntrack))
      else
        if(neutrals_time_mod.gt.0.0_R8) then
          write(*,*) 'NEUTRALS: ',tim,mod(tim,neutrals_time_mod),
     1     neutrals_time_switch,int(tim/neutrals_time_mod),
     2     neutrals_count
          if(mod(tim,neutrals_time_mod).gt.neutrals_time_switch.or.
     1     int(tim/neutrals_time_mod).gt.neutrals_count) then
            neutrals_count=tim/neutrals_time_mod
            call read_neutrals_namelist(nx, ny, ns, .true.)
            new_sputter_namelist=.true.
            call xertst (nstrat.le.nstraid, 
     &       'Increase DEF_NSTRA in DIMENSIONS.F')
            call write_b2mod_neutrals_namelist
          endif
        endif
      endif
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
!       call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                    !srv 11.01.13
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
*   ..initialize sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0_R8, smo0, 1)
      call sfill (n2*2, 0.0_R8, sne0, 1)
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)
      call sfill (n2*nscx, 0.0_R8, shi0_ff, 1)

      call sfill (n2*2*ns, 0.0_R8, sna0_eir, 1)                      !srv 29.04.10 {
      call sfill (n2*4*ns, 0.0_R8, smo0_eir, 1)
      call sfill (n2*4, 0.0_R8, she0_eir, 1)
      call sfill (n2*4, 0.0_R8, shi0_eir, 1)
      fatm_eir = 0.0_R8
      datm_eir = 0.0_R8                                            !srv 29.04.10 }

      tchem=0.0_R8
      tchee=0.0_R8
      tphys=0.0_R8
      tphye=0.0_R8
      thevp=0.0_R8
      thvpe=0.0_R8
      trese=0.0_R8
      tresn=0.0_R8
      trfle=0.0_R8
      trfln=0.0_R8

      if (main_call) then
        old_erosion(1:nbnd,1:ntargsp)=erosion(1:nbnd,1:ntargsp)
        old_deposition(1:nbnd,1:ntargsp)=deposition(1:nbnd,1:ntargsp)
      endif

csw 09oct2012 SOLPS4: explicitly avoid b2stbr call here
c   --> better solution wanted here..
      if(no_b2stbr_phys_call.eq.0 ) then
       if(redep_alpha.eq.0.0_R8) then
        f_redep=0.0_R8
       else
        do is=0, ns-1
          if(is_neutral(is).and.lnext(is,is+1)) then
            do iy=-1,ny
              do ix=-1,nx
                f_redep(ix,iy,is) = 1.0_R8 /
     1           ( 1.0_R8 + redep_alpha*(qe*bb(ix,iy,3)/(am(is)*mp)/
     2           (rsa(ix,iy,is)*ne(ix,iy)))**2 )
              enddo
            enddo
            if(main_call.and.output.ge.1) then
              write(*,'(a,i3,16384(1x,f9.6))') 'is,f_redep_W ',
     1         is,f_redep(0,:,is)
              write(*,'(a,i3,16384(1x,f9.6))') 'is,f_redep_E ',
     1         is,f_redep(nx-1,:,is)
            endif
          else
            f_redep(:,:,is)=0.0_R8
          endif
        enddo
       endif
       if (neutrals_namelist.eq.1 .or. use_eirene.ne.0) then
        if(ncall.eq.0) write(*,*) 'b2stbr: using b2stbr_phys'
        call b2stbr_phys (nx, ny, ns, nscx, iscx, dtim, BoRiS,
     &   vol, hx, hy, gs, qz, pbs, na, ua, kinrgy,
     &   te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &   sna0, smo0, she0, shi0, sch0, sne0, tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &   neutral_sources_rescale, core_sources_rescale,  
     &   sput_src, sput_chem_alpha, sput_phys_alpha, alpha,
     &   sput_chem_model, therm_evap, sput_res,  
     &   reflection_on, sputter_energy_on, main_call, output, 
     &   sput_frac_flag, sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1   sput_mixed_alpha, sput_mixed_beta,
     2   new_sputter_namelist, shi0_ff, f_redep,
     3   potential_at_guard_cell, recycled_neutrals_contr)
       else
        if(ncall.eq.0) write(*,*) 'b2stbr: using b2stbr_bas'
        call b2stbr_bas (nx, ny, ns, nscx, iscx, dtim, BoRiS,
     &   vol, hx, hy, gs, qz, pbs, na, ua, kinrgy,
     &   te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &   sna0, smo0, she0, shi0, sch0, sne0, tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &   coreregno, coreregn2, output,
     &   neutral_sources_rescale, core_sources_rescale,  
     &   sput_src, sput_frc, sput_phys,
     &   sput_chem_alpha, sput_phys_alpha, alpha,
     &   sput_chem_model, therm_evap, sput_res, 
     &   reflection_on, sputter_energy_on, main_call,
     &   sput_chem_cutoff_alpha, sput_chem_cutoff_beta, shi0_ff,
     &   f_redep, potential_at_guard_cell, recycled_neutrals_contr)
       endif
      endif

      if(first_flight.ne.0) then
        start_time=x05baf()
        if (output.ge.1) then
          open(88,file='start')
          open(89,file='end')
        endif
        do is = 0, nscx-1
          if (output.ge.1) then
            do iy=-1,ny
              do ix=-1,nx
                write(88,*) ix,iy,sna0(ix,iy,0,iscx(is)),
     .                            shi0(ix,iy,0),she0(ix,iy,0)
              enddo
              write(88,'()')
            enddo
          endif
          sum_start=0.0_R8
          sum_strti=sum(shi0_ff(:,:,is))
          pfrpuff=0.0_R8
          mainpuff=0.0_R8
          corepuff=0.0_R8
          work=0.0_R8
          wrk2=0.0_R8
          wrk3=0.0_R8
          ic=0
          do ix=-1,nx
            do ib=-1,ny
              if((leftix(ix,ib).eq.-2 .or. rightix(ix,ib).eq.nx+1)
     &          .and. region(ix,ib,0).ne.0) then
                ic=ic+1
                call process_b2mod_first_flight
     *           (ic,ix,ib,nx,ny,ns,nscx,iscx,iscx(is),jsep,
     *            dtim,sna0(ix,ib,0,iscx(is)),shi0_ff(ix,ib,is),
     *            work,wrk2,wrk3,
     *            pfrpuff,mainpuff,corepuff)
                shi0(ix,ib,0)=shi0(ix,ib,0)-shi0_ff(ix,ib,is)
                sum_start=sum_start+sna0(ix,ib,0,iscx(is))
                sna0(ix,ib,0,iscx(is))=0.0_R8
                shi0_ff(ix,ib,is)=0.0_R8
              endif
            enddo
          enddo
          do ib=-1,nx
            do iy=-1,ny
              if((topiy(ib,iy).eq.ny+1 .or. bottomiy(ib,iy).eq.-2)
     &          .and. region(ib,iy,0).ne.0) then
                ic=ic+1
                call process_b2mod_first_flight
     *           (ic,ib,iy,nx,ny,ns,nscx,iscx,iscx(is),jsep,
     *            dtim,sna0(ib,iy,0,iscx(is)),shi0_ff(ib,iy,is),
     *            work,wrk2,wrk3,
     *            pfrpuff,mainpuff,corepuff)
                shi0(ib,iy,0)=shi0(ib,iy,0)-shi0_ff(ib,iy,is)
                sum_start=sum_start+sna0(ib,iy,0,iscx(is))
                sna0(ib,iy,0,iscx(is))=0.0_R8
                shi0_ff(ib,iy,is)=0.0_R8
              endif
            enddo
          enddo
          sum_end=sum(work(:,:,:))
          sum_ndi=sum(wrk2(:,:))
          sum_nde=sum(wrk3(:,:))
          if(main_call)
     >     write(*,*) 'first_flight: is = ',is,
     *      pfrpuff,mainpuff,corepuff,sum_start,sum_end,
     *      sum_start-(sum_end+pfrpuff+mainpuff+corepuff),
     *      sum_strti,sum_ndi,sum_nde
          if (sum_end.ne.0.0_R8) then
            sna0(:,:,0,iscx(is))=sna0(:,:,0,iscx(is))
     .       +work(:,:,0)*sum_start/sum_end
            sna0(:,:,0,iscx(is)+1)=sna0(:,:,0,iscx(is)+1)
     .       +work(:,:,1)*sum_start/sum_end
            she0(:,:,0)=she0(:,:,0)+wrk3(:,:)*sum_start/sum_end
          endif
          if (sum_ndi.ne.0.0_R8) 
     >     shi0(:,:,0)=shi0(:,:,0)+wrk2(:,:)*sum_strti/sum_ndi
          if (output.gt.1) then
            do iy=-1,ny
              do ix=-1,nx
                write(89,*) ix,iy,sna0(ix,iy,0,iscx(is))
     .                           ,shi0(ix,iy,0),she0(ix,iy,0)
              enddo
              write(89,'()')
            enddo
          endif
          if (main_call) then
            write(*,*) 'Total first-flight-source na:', iscx(is), 
     &       sum(sna0(:,:,0,iscx(is)))
            write(*,*) 'Total first-flight-source na:', iscx(is)+1, 
     &       sum(sna0(:,:,0,iscx(is)+1))
          endif
        enddo
        if (main_call) then
          write (*,*)
     >     'Total first-flight-source ion power=', sum(shi0(:,:,0))
          write (*,*)
     >     'Total first-flight-source electron power=', sum(she0(:,:,0))
        endif
        if (output.ge.1) then
          close(88)
          close(89)
        endif
        end_time=x05baf()
        if(main_call)
     >   write(*,*) end_time-start_time,' cpu-seconds'
      endif

      if(main_call.and.out_tfiles.ne.0) then                        ! xpb  !srv 07.04.09
        write(*,'(a,1p,(10000e12.4))') 'tphys',sum(tphys),tphys
        write(*,'(a,1p,(10000e12.4))') 'tphye',sum(tphye),tphye
        write(*,'(a,1p,(10000e12.4))') 'tchem',sum(tchem),tchem
        write(*,'(a,1p,(10000e12.4))') 'tchee',sum(tchee),tchee
        write(*,'(a,1p,(10000e12.4))') 'thevp',sum(thevp),thevp
        write(*,'(a,1p,(10000e12.4))') 'thvpe',sum(thvpe),thvpe
        write(*,'(a,1p,(10000e12.4))') 'trese',sum(trese),trese
        write(*,'(a,1p,(10000e12.4))') 'tresn',sum(tresn),tresn
        write(*,'(a,1p,(10000e12.4))') 'trfln',sum(trfln),trfln
        write(*,'(a,1p,(10000e12.4))') 'trfle',sum(trfle),trfle
      endif                                                         ! xpb

      if(use_eirene.ne.0) then
#ifdef B25_EIRENE
        call eirene_mc(nx, ny, ns, ismain, dtim, BoRiS, sna0_eir,     !srv 29.04.10 {
     &                 smo0_eir, she0_eir, shi0_eir, main_call)

csw 25feb2014  fix IVS's version of source term averaging (running average)
c fixed SMOOTH SOURCES FROM EIRENE , averaging procedure
c (25feb2014 different from SOLPS4.3 averaging procedure as implemented in
c eirene_mc)
        if(nhist .gt. 1 .and. main_call) then
          ihist = ihist + 1
          if (ihist.eq.nhist) full_history = .true.
          if (ihist.eq.nhist+1) ihist = 1
          write(*,*) 'b2stbr: EIRENE source averaging ', nhist, ihist

          sna0_eir_h(:,:,:,:,ihist) = sna0_eir
          smo0_eir_h(:,:,:,:,ihist) = smo0_eir
          she0_eir_h(:,:,:,ihist) = she0_eir
          shi0_eir_h(:,:,:,ihist) = shi0_eir
          do ia = 1,nnatmi
            rfluxa_h(:,:,ia,1,ihist) = rfluxa(:,:,ia,1)
          enddo
          do im = 1,nnmoli
            rfluxm_h(:,:,im,1,ihist) = rfluxm(:,:,im,1)
          enddo                  
          
          if (ihist.gt.1) then
            sna0_eir = 0.0_R8
            smo0_eir = 0.0_R8
            she0_eir = 0.0_R8
            shi0_eir = 0.0_R8
            rfluxa = 0.0_R8      
            rfluxm = 0.0_R8      
csw 25feb2014  fix IVS's version of source term averaging
            if (full_history) then
              ih_end = nhist
            else
              ih_end = ihist
            end if
            do iih=1,ih_end
              sna0_eir_s = sna0_eir_s + sna0_eir_h(:,:,:,:,iih)
              smo0_eir_s = smo0_eir_s + smo0_eir_h(:,:,:,:,iih)
              she0_eir_s = she0_eir_s + she0_eir_h(:,:,:,iih)
              shi0_eir_s = shi0_eir_s + shi0_eir_h(:,:,:,iih) 
              do ia = 1,nnatmi
                rfluxa_s(:,:,ia,1) = rfluxa_s(:,:,ia,1) 
     .                           + rfluxa_h(:,:,ia,1,iih)
              enddo
              do im = 1,nnmoli
                rfluxm_s(:,:,im,1) = rfluxm_s(:,:,im,1)
     .                             + rfluxm_h(:,:,im,1,iih)
              enddo                  
            enddo
csw 25feb2014  fix IVS's version of source term averaging
c            sna0_eir = sna0_eir / min(dble(nhist),dble(ncall+1))
c            smo0_eir = smo0_eir / min(dble(nhist),dble(ncall+1))
c            she0_eir = she0_eir / min(dble(nhist),dble(ncall+1))
c            shi0_eir = shi0_eir / min(dble(nhist),dble(ncall+1))
c            rfluxa = rfluxa / min(dble(nhist),dble(ncall+1))                    
c            rfluxm = rfluxm / min(dble(nhist),dble(ncall+1))                    
            sna0_eir = sna0_eir_s / dble(ih_end)
            smo0_eir = smo0_eir_s / dble(ih_end)
            she0_eir = she0_eir_s / dble(ih_end)
            shi0_eir = shi0_eir_s / dble(ih_end)
            rfluxa = rfluxa_s / dble(ih_end)                    
            rfluxm = rfluxm_s / dble(ih_end)                    

!----  end of averaging procedure -----------------
          endif
        endif                                                            !IYS 28.03.12}
        if(ncall.le.999. and. ncall.le.ncall_iout) then
         do is = 0, ns-1
          write (chns,'(i3.3)') is
          write (ch3 ,'(i3.3)') ncall
          if(.not.is_neutral(is)) then
           call my_out(70,nx,ny,sna0_eir(:,:,0,is)+sna0_eir(:,:,1,is)* 
     &                    na(:,:,is),'b2stbr_sna_eir'//chns//'_'//ch3)
          endif
         enddo
        endif
        sna0 = sna0 + sna0_eir
        smo0 = smo0 + smo0_eir
        she0 = she0 + she0_eir
        shi0 = shi0 + shi0_eir
#ifdef ASTRA
        if (use_astra.ne.0) then
          do is = 0, ns-1                                                 !srv 28.07.11 20.02.12 {
            if(is.eq.ismain) then
              call b2_to_astra_src(nx, ny, vol, sna0_eir(:,:,0,is)+
     &         sna0_eir(:,:,1,is)*na(:,:,is), svol, se)
            endif
          enddo                                                           !srv 28.07.11 20.02.12 }
          call b2_to_astra_src(nx, ny, vol, she0_eir(:,:,0)+
     &     she0_eir(:,:,1)*te(:,:)+she0_eir(:,:,2)*ne+
     &     she0_eir(:,:,3)*ne*te, svol, seh)                            !srv 11.06.13
          call b2_to_astra_src(nx, ny, vol, shi0_eir(:,:,0)+
     &     shi0_eir(:,:,1)*ti+shi0_eir(:,:,2)*ni(:,:,0)+
     &     shi0_eir(:,:,3)*ni(:,:,0)*ti, svol, sih)                     !srv 11.06.13
        endif
#endif	                       
        do ia = 1, natm
         is = eb2spcr(ia)
         do iy = -1,ny
          do ix = -1,nx
           fatm_eir(ix,iy,0,is) = pfluxa(ix+1,iy+1,ia,1)
           fatm_eir(ix,iy,1,is) = rfluxa(ix+1,iy+1,ia,1)
           datm_eir(ix,iy,is)   = dab2(ix+1,iy+1,ia,1)
          enddo
         enddo
        enddo
#else
        call xerrab ('b2stbr not compiled with the B25_EIRENE option')
#endif
      else
cwdk remove this if-test, to allow runs without Eirene if code is compiled with Eirene
        call eirene_mc_init             ! EIRDIAG = 0   !srv 03.04.13
      endif

* ..sources
      b2stbr_sna(-1:nx,-1:ny,0:ns-1)=
     1 sna0(-1:nx,-1:ny,0,0:ns-1)+
     2 sna0(-1:nx,-1:ny,1,0:ns-1)*na(-1:nx,-1:ny,0:ns-1)
      do is=0,ns-1
        b2stbr_smo(-1:nx,-1:ny,is)=
     1   smo0(-1:nx,-1:ny,0,is)+
     2   smo0(-1:nx,-1:ny,1,is)*ua(-1:nx,-1:ny,is)+
     3   smo0(-1:nx,-1:ny,2,is)*am(is)*mp*na(-1:nx,-1:ny,is)+
     4   smo0(-1:nx,-1:ny,3,is)*am(is)*mp*na(-1:nx,-1:ny,is)*
     5   ua(-1:nx,-1:ny,is)
      enddo
      b2stbr_she(-1:nx,-1:ny)=
     1 she0(-1:nx,-1:ny,0)+
     2 she0(-1:nx,-1:ny,1)*te(-1:nx,-1:ny)+
     3 she0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 she0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*te(-1:nx,-1:ny)
      b2stbr_shi(-1:nx,-1:ny)=
     1 shi0(-1:nx,-1:ny,0)+
     2 shi0(-1:nx,-1:ny,1)*ti(-1:nx,-1:ny)+
     3 shi0(-1:nx,-1:ny,2)*ni(-1:nx,-1:ny,0)+
     4 shi0(-1:nx,-1:ny,3)*ni(-1:nx,-1:ny,0)*ti(-1:nx,-1:ny)
      b2stbr_sch(-1:nx,-1:ny)=
     1 sch0(-1:nx,-1:ny,0)+
     2 sch0(-1:nx,-1:ny,1)*po(-1:nx,-1:ny)+
     3 sch0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 sch0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*po(-1:nx,-1:ny)
      b2stbr_sne(-1:nx,-1:ny)=
     1 sne0(-1:nx,-1:ny,0)+
     2 sne0(-1:nx,-1:ny,1)*ne(-1:nx,-1:ny)

      if(main_call) then
        do ibnd = 1, nbnd
          over_erosion = .false.
          equilibrium(1:ntargsp) = .false.
          t_frac = 1.0_R8
          do itrack = 1, ntargsp
            if (particle_inventory(ibnd,itrack).gt.0.0_R8 .or.
     1          itrack.eq.itrack_0(ibnd)) then
              t0 = particle_inventory(ibnd,itrack)+
     1         (deposition(ibnd,itrack)-old_deposition(ibnd,itrack))-
     2         (erosion(ibnd,itrack)-old_erosion(ibnd,itrack))
            else ! case of continuous erosion of deposited species
c-xpb  This only has meaning if we are following multiple mixed species and
c-xpb  this is not the first timestep 
              t0 = 0.0_R8
              equilibrium(itrack) = ntargsp.gt.1 .and. ncall.gt.0 .and.
     &         sput_frac_flag.ne.0
            endif
            if (t0.lt.0.0_R8) then  ! we have over-eroded, so we correct
              if (itrack.ne.itrack_0(ibnd)) then  ! over-erosion of a deposited species
                t_frac = min(t_frac, 
     1            particle_inventory(ibnd,itrack) /
     2           (particle_inventory(ibnd,itrack) - t0))
                over_erosion = .true.
              else  ! we are keeping track of the bulk and eroded it all
c-xpb  sput_frac_flag.eq.0 means we are not following mixed materials
c-xpb  plate_thickness.eq.0 means we are not doing full inventories 
                call xertst (sput_frac_flag.eq.0 .or. 
     .                       plate_thickness(ibnd).eq.0.0_R8,
     1           'Eroded the entire bulk contents of the target!')
              endif
            endif
          enddo
          if (.not.over_erosion) then
            plate_time(ibnd)=plate_time(ibnd)+
     1       dtim*plate_time_factor(ibnd)
            do itrack = 1, ntargsp
              if (.not.equilibrium(itrack)) then
                if (sput_frac_flag.ne.0) 
     >           particle_inventory(ibnd,itrack)=
     .            particle_inventory(ibnd,itrack)+
     1               (    deposition(ibnd,itrack)-
     .                old_deposition(ibnd,itrack))-
     2               (       erosion(ibnd,itrack)-
     .                   old_erosion(ibnd,itrack))
              else ! we only erode what we just deposited
                particle_inventory(ibnd,itrack)=0.0_R8
                erosion(ibnd,itrack)=old_erosion(ibnd,itrack)+
     1           deposition(ibnd,itrack)-old_deposition(ibnd,itrack)
                lnorm=am(target_intcode(itrack))/tgdns(itrack)/
     &           avogr/1.0e6_R8/monolayer_l(itrack)
                if(plate_area(ibnd).gt.0.0_R8) then
                  monolayer_erosion(ibnd,itrack)=
     &             erosion(ibnd,itrack)/plate_area(ibnd)*lnorm
                endif
              endif
            enddo
          else 
            plate_time(ibnd)=plate_time(ibnd)+
     1       dtim*plate_time_factor(ibnd)*t_frac
            do itrack = 1, ntargsp
              particle_inventory(ibnd,itrack)= max( 0.0_R8,
     .         particle_inventory(ibnd,itrack)+t_frac*
     1         (    deposition(ibnd,itrack)-
     .          old_deposition(ibnd,itrack))-
     2         (       erosion(ibnd,itrack)-
     .             old_erosion(ibnd,itrack)) )
              deposition(ibnd,itrack)=old_deposition(ibnd,itrack)+
     1          t_frac*
     2         (deposition(ibnd,itrack)-old_deposition(ibnd,itrack))
              if (.not.equilibrium(itrack)) then
                erosion(ibnd,itrack)=old_erosion(ibnd,itrack)+
     1            t_frac*
     2           (erosion(ibnd,itrack)-old_erosion(ibnd,itrack))
              else ! we only erode what we just deposited
                erosion(ibnd,itrack)=old_erosion(ibnd,itrack)+
     1           deposition(ibnd,itrack)-old_deposition(ibnd,itrack)
              endif
              lnorm=am(target_intcode(itrack))/tgdns(itrack)/
     &         avogr/1.0e6_R8/monolayer_l(itrack)
              if(plate_area(ibnd).gt.0.0_R8) then
                monolayer_deposition(ibnd,itrack)=
     &           deposition(ibnd,itrack)/plate_area(ibnd)*lnorm
                monolayer_erosion(ibnd,itrack)=
     &           erosion(ibnd,itrack)/plate_area(ibnd)*lnorm
              endif
            enddo
          endif
        enddo
      endif

      if(b2wall_netcdf.ne.0 .and. ntargsp.gt.0) then
#ifndef NO_CDF
        if(ncall.eq.0) then
          call define_cdf_b2mod_wall(nbnd,ntargsp)
          status=nf_open('b2wall.nc',NF_WRITE,ncid)
          call check_cdf_status(status)
          write(*,*) 'Writing track_species'
          call write_cdf_char_variable(ncid,-1,
     *     'track_species',track_species(1:ntargsp))
          write(*,*) 'Writing plate_area'
          call write_cdf_variable(ncid,-1,
     *     'plate_area',plate_area(1:nbnd))
          status=nf_close(ncid)
          call check_cdf_status(status)
          write(*,*) 'Closed b2wall.nc'
          ntstep_b2wall=0
        endif
        if(main_call.and.mod(ncall,b2wall_netcdf).eq.0) then
          ntstep_b2wall=ntstep_b2wall+1
          status=nf_open('b2wall.nc',NF_WRITE,ncid)
          call check_cdf_status(status)
          call write_cdf_variable(ncid,ntstep_b2wall,
     *     'plate_time',plate_time(1:nbnd))
          call write_cdf_variable(ncid,ntstep_b2wall,
     *     'erosion',erosion(1:nbnd,1:ntargsp))
          call write_cdf_variable(ncid,ntstep_b2wall,
     *     'deposition',deposition(1:nbnd,1:ntargsp))
          call write_cdf_variable(ncid,ntstep_b2wall,
     *     'monolayer_erosion',monolayer_erosion(1:nbnd,1:ntargsp))
          call write_cdf_variable(ncid,ntstep_b2wall,
     *     'monolayer_deposition',
     *     monolayer_deposition(1:nbnd,1:ntargsp))
          status=nf_close(ncid)
          call check_cdf_status(status)
        endif
#else
        call xerrab('CDF routines called when no CDF available')
#endif
      endif

      if (iout.eq.1 .or. iout_b2npmo.eq.1) then                       !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if (use_eirene.ne.0) then
         if(.not.is_neutral(is)) then
          do k=0,1 ! 3                                                  !srv 29.04.10 {
           write (chk,'(i1)') k
            call my_out(70,nx,ny,smo0_eir(-1,-1,k,is),
     *              'b2stbr_smo_eir'//chk//chns)
          enddo
          call my_out(70,nx,ny,smo0_eir(:,:,0,is)+
     &        smo0_eir(:,:,1,is)*ua(:,:,is)+
     &        smo0_eir(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smo0_eir(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2stbr_smo_eir'//chns)
         endif
        else
!	  do k=0,3
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,smo0(-1,-1,k,is),'b2stbr_smo'//chk//chns)
!	  enddo
         call my_out(70,nx,ny,b2stbr_smo(-1,-1,is),'b2stbr_smo'//chns)
        endif
!	   call my_out(70,nx,ny,ua(-1,-1,is),'b2stbr_ua'//chns)
!	   call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
       enddo                                                          !srv 29.04.10 }
      endif
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!	  call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
        if(use_eirene.ne.0) then
         if(.not.is_neutral(is)) then
          do k=0,1                                                      !srv 29.04.10
           write (chk,'(i1)') k
           call my_out(70,nx,ny,sna0_eir(-1,-1,k,is),
     *               'b2stbr_sna_eir'//chk//chns)
          enddo
          call my_out(70,nx,ny,sna0_eir(:,:,0,is)+ 
     &             sna0_eir(:,:,1,is)*na(:,:,is),'b2stbr_sna_eir'//chns)
         endif
        else
!	  do k=0,1
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,sna0(-1,-1,k,is),'b2stbr_sna'//chk//chns)
!	  enddo
         call my_out(70,nx,ny,b2stbr_sna(-1,-1,is),'b2stbr_sna'//chns)
        endif
       enddo                                                          !srv 29.04.10
#ifdef B25_EIRENE
       do ia = 1,nnatmi
        write (chns,'(i2.2)') ia
        call my_out(70,nx,ny,rfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2stbr_rfluxa_eir'//chns)
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2stbr_pfluxa_eir'//chns)
        call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,ia,1),
     &                'b2stbr_dab2_eir'//chns)
       enddo
       do im = 1,nmol
        write (chns,'(i2.2)') im
        call my_out(70,nx,ny,rfluxm(0:nx+1,0:ny+1,im,1),
     &              'b2stbr_rfluxm_eir'//chns)
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,im,1),
     &              'b2stbr_pfluxm_eir'//chns)
        call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,im,1),
     &              'b2stbr_dmb2_eir'//chns)
       enddo
#endif
      endif
      if (iout.eq.1 .or. (iout_b2npht.eq.1 .or. iout_b2nph9.eq.1)) then
       if (use_eirene.ne.0) then
        do k=0,1  ! 3
         write (chk,'(i1)') k
         call my_out(70,nx,ny,she0_eir(-1,-1,k),'b2stbr_she_eir'//chk)
         call my_out(70,nx,ny,shi0_eir(-1,-1,k),'b2stbr_shi_eir'//chk)
        enddo
        call my_out(70,nx,ny,she0_eir(:,:,0)+she0_eir(:,:,1)*ti+
     &     she0_eir(:,:,2)*ne+she0_eir(:,:,3)*ti*ne,'b2stbr_she_eir')
        call my_out(70,nx,ny,shi0_eir(:,:,0)+shi0_eir(:,:,1)*ti+
     &     shi0_eir(:,:,2)*ni(:,:,0)+shi0_eir(:,:,3)*ti*ni(:,:,0),
     &     'b2stbr_shi_eir')
       endif
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,she0(-1,-1,k),'b2stbr_she'//chk)
!	  call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbr_shi'//chk)
!	 enddo
       call my_out(70,nx,ny,b2stbr_she,'b2stbr_she')
       call my_out(70,nx,ny,b2stbr_shi,'b2stbr_shi')               !srv 29.04.10
      endif
      if (iout.eq.1 .or. (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)) then
       if (use_eirene.eq.0) then
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbr_sch'//chk)
!	 enddo
        call my_out(70,nx,ny,b2stbr_sch,'b2stbr_sch')
       endif
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return
      end subroutine b2stbr

*-----------------------------------------------------------------------
*.specification

*.end b2stbr

      subroutine b2stbr_neutr_scl(nx, ny, ns)
      implicit none
      integer nx, ny, ns
      integer ncall, use_eirene, ank_mods
      save ncall, use_eirene, ank_mods
      data ncall/0/, use_eirene/0/, ank_mods/0/
      call subini ('b2stbr_neutr_scl')
      if (ncall.eq.0) then
        call ipgeti('b2mndr_eirene',use_eirene)
        call ipgeti('eirene_ank_mods',ank_mods)
      endif
      if(use_eirene.ne.0.and.ank_mods.ne.0) then
#ifdef B25_EIRENE
c        write(*,*) 'Calling neutr_scl'
        call neutr_scl(nx, ny, ns)
c        write(*,*) 'Called neutr_scl'
#else
        call xerrab ('b2stbr not compiled with the B25_EIRENE option')
#endif
      endif
      ncall = ncall+1
      call subend ()
      return

      end subroutine b2stbr_neutr_scl

*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr_init (nx, ny, ns, BoRiS)
      use b2mod_types
      use b2mod_time
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_sputter
      use b2mod_indirect
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
#ifndef ALLOCATE_AND_NAMELIST
#     include "DIMENSIONS.F"
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nx, ny, ns
      real (kind=R8), intent(in) :: BoRiS
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBR_INIT reads the initial data necessary for the neutrals
*     recycling and wall model boundary conditions.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nx_tmp, ny_tmp, ns_tmp
      integer itargsp, is, istra
      integer use_eirene, neutrals_namelist
      integer sput_src, sput_chem_model
      real (kind=R8) :: sput_frc
      logical new_sputter_namelist, catch_up, wall_namelist_used
      character*256 filename
*   ..procedures
      intrinsic mod, nint
      external subini, subend, xertst, xerrab, ipgeti, ipgetr
*   ..initialisation
      data use_eirene/0/, neutrals_namelist/0/, 
     *  sput_src/1/, sput_chem_model/0/, sput_frc/0.0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbr_init')
      write(*,*) 'Start of b2stbr_init'
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read in needed data
      call ipgeti ('b2mndr_eirene', use_eirene)
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call ipgeti ('b2stbr_sput_src', sput_src)
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      call ipgetr ('b2stbr_sput_frc', sput_frc)
      call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
      if(sput_frc.ne.0.0_R8) then
        call xertst (0.le.sput_src.and.sput_src.lt.ns,
     1     'faulty b2stbr_sput_src')
        call xertst (sput_dst.lt.ns, 'faulty b2stbr_sput_dst')
        call xertst (sput_dst2.lt.ns, 'faulty b2stbr_sput_dst2')
        call xertst (sput_dst3.lt.ns, 'faulty b2stbr_sput_dst3')
        call xertst (0.0_R8.le.sput_frc, 'faulty b2stbr_sput_frc')
        if(sput_dst.ge.0 .and. use_eirene.eq.0 .and.
     .     neutrals_namelist.eq.0) then
          if (sput_chem_model.eq.1) then
            write(*,'(1x,a,i3,a,i3,a,g10.3)') 
     .       'Chemical sputtering from ',sput_src,
     1       ' producing ',sput_dst,' at rate ',sput_frc
          else if (sput_chem_model.eq.0) then
            write(*,'(1x,a,i3,a,i3,a,g10.3)') 
     .       'Chemical sputtering from ',sput_src,
     1       ' producing ',sput_dst,
     2       ' with empirical formula multiplier ',sput_frc
          endif
        endif
        if(sput_dst2.ge.0)
     1     write(*,*) 'Additional sputtering producing ',sput_dst2
        if(sput_dst3.ge.0)
     1     write(*,*) 'Additional sputtering producing ',sput_dst3
      endif
      call ipgeti ('b2stbr_plate_model', plate_model)
      if (plate_model.eq.0) then
        write(*,*) 
     .     'Using O-D time-independent plate temperature model'
      else if (plate_model.eq.1) then
        write(*,*) 'Using 1-D time-dependent plate temperature model'
      else if (plate_model.eq.2) then
        write(*,*) 'Using 2-D time-dependent plate temperature model'
      else
        call xerrab ('Invalid value for b2stbr_plate_model')
      endif
      call ipgeti ('b2stbr_plate_option', plate_option)
      call xertst (1.le.plate_option .and.
     &             3.ge.plate_option, 'faulty b2stbr_plate_option')
      call ipgetr ('b2stbr_plate_temp', plate_temp)
      call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
      call ipgetr ('b2stbr_plate_thick', plate_thick)
      call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')

c-xpb Reading of neutrals namelist moved here as it is needed by init_wall
      new_sputter_namelist=neutrals_namelist.eq.1.or.use_eirene.ne.0
      call read_neutrals_namelist
     1   (nx, ny, ns, new_sputter_namelist)
      call xertst (nstrat.le.nstraid, 
     &   'Increase DEF_NSTRA in DIMENSIONS.F')
      call alloc_b2mod_neutr_src_scaling(ns,nstrat)
      if(new_sputter_namelist) call write_b2mod_neutrals_namelist
      if(neutrals_time_mod.gt.0.0_R8) then
        neutrals_count=tim/neutrals_time_mod
      else
        neutrals_count=0
      endif
      catch_up=neutrals_time_mod.gt.0.0_R8
      do while(catch_up)
        if(mod(tim,neutrals_time_mod).gt.neutrals_time_switch) then
          write(*,*) 'b2stbr_init: catching up ...'
          call read_neutrals_namelist(nx, ny, ns, .true.)
          new_sputter_namelist=.true.
          call xertst (nstrat.le.nstraid, 
     &       'Increase DEF_NSTRA in DIMENSIONS.F')
          call write_b2mod_neutrals_namelist
          catch_up=neutrals_time_mod.gt.0.0_R8
        else
          catch_up=.false.
        endif
      enddo
      if (neutrals_namelist.eq.1 .and. use_eirene.eq.0) then
        do istra = 1, nstrai
          if (crcstra(istra).eq.'S' .or. crcstra(istra).eq.'N' .or.
     &        crcstra(istra).eq.'W' .or. crcstra(istra).eq.'E') then
            write(6,'(a,i2,a)') 'For recycling boundary ',istra,' : '
            do is = 0, ns-1
              if (phys_sput(is,istra).gt.0.0_R8) then
                write(6,'(a,i3,a,i3,a,g10.3)') 
     .           'Physical sputtering from species ', is,
     &           ' producing species ', sput_dst, 
     &           ' with rate multiplied by ', phys_sput(is,istra)
              else
                write(6,'(a,i3,a)') 
     &           'No physical sputtering due to species ', is,'.'
              end if
              if (chem_sput(is,istra).gt.0.0_R8 .and.
     &            targsp(istra,1).ne.-1) then
                if (sput_chem_model.eq.0) then
                  write(6,'(a,i3,a,g10.3)')  
     .             'Chemical sputtering from species ', is,
     &             ' with rate multiplied by ', chem_sput(is,istra)
                else if (sput_chem_model.eq.1) then
                  write(6,'(a,i3,a,g10.3)')  
     .             'Constant chemical sputtering from species ', is,
     &             ' with rate ', chem_sput(is,istra)
                end if
              else
                write(6,'(a,i3,a)') 
     &           'No chemical sputtering due to species ', is,'.'
              end if
            end do
          end if
        end do
        if(sput_dst2.ge.0)
     1   write(*,*) 'Additional sputtering producing ',sput_dst2
        if(sput_dst3.ge.0)
     1   write(*,*) 'Additional sputtering producing ',sput_dst3
      end if

*   ..reading wall and sputter namelist data on first call
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns)
      write(*,*) 'Calling alloc_b2mod_wall'
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
      call alloc_b2mod_wall(ns,natm,nmol,nion)
      call init_wall(nx,ny)
      call read_b2mod_wall_namelist(99,ns,'    ',wall_namelist_used)
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0, 
     .   'Inconsistent values of ndepth_nml and plate_model')
      call get_wall_data (nx,ny,ns,BoRiS,wall_namelist_used)
#ifdef ALLOCATE_AND_NAMELIST
      nx_tmp = nx
      ny_tmp = ny
      ns_tmp = ns
#else
      nx_tmp = DEF_NXD
      ny_tmp = DEF_NYD
      ns_tmp = DEF_NSD
#endif
      filename='b2.sputter_save.parameters'
      inquire(file='../b2.sputter_save.parameters',
     1 exist=sputter_namelist_used)
      if (sputter_namelist_used) then
        call read_b2mod_sputter_namelist
     .     (99,'    ',sputter_namelist_used)
      else
        sputter_yield =0.0_R8
        sputter_yield2=0.0_R8
      endif
      nx_sp = nx_tmp
      ny_sp = ny_tmp
      ns_sp = ns_tmp  

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbr_init

      end subroutine b2stbr_init

