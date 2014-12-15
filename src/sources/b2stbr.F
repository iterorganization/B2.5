*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr (nx, ny, ns, nscx, nscxmax, iscx, ismain, 
     &  dtim, BoRiS, vol, hx, hy, gs, qz, pbs, na, ua, kinrgy, 
     &  te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,            !srv 12.04.12
     &  sna0, smo0, she0, shi0, sch0, wrk1, main_call)   !dpc
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
     1 , only : b2stbr_sna,b2stbr_smo,b2stbr_she,b2stbr_shi,b2stbr_sch
#ifdef USE_EIRENE
      use b2mod_work, only : sna0_eir_h, smo0_eir_h, she0_eir_h,       !IYS 28.03.12 {
     &    shi0_eir_h, rfluxa_h, rfluxm_h, max_eir_hist_slc,
csw 25feb2014 fix IVS's version of source term avergaing
     &    alloc_hist_arrays, clear_hist_arrays                         !IYS 28.03.12 } 
#endif
      use b2mod_b2cmpa
      use b2mod_b2cmpb
#ifndef NO_CDF
      use b2mod_tallies
#endif
#ifdef ASTRA
      use b2mod_b2_to_astra, only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 20.02.12 11.06.13
#endif	                       
#ifndef ALLOCATE_AND_NAMELIST
#     include "DIMENSIONS.F"
#endif
      use eirmod_extrab25
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
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1)
      logical main_call                                                  !dpc
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3) 
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
*     density, momentum and heat, due to special boundary conditions
*     and the neutral gas recycling model.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is, ib, ic, ibnd, itrack, elm_count, k, ia, im
      integer ihist,iih                                                  !IYS 29.03.2012
      integer, save :: nx_sp, ny_sp, ns_sp
      integer, save :: no_b2stbr_phys_call=0
      integer, save :: out_tfiles = 0                                !srv 07.04.09
      integer, save :: iout = 0,        iout_b2npmo = 0,             !srv 11.09.09 {
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0,
     &                 ncall_iout = -1,                              !srv 23.03.12
     &                 nhist = 1                                     !IYS 28.03.12
      character chns*2, chk*1, ch3*3                                 !srv 11.09.09 } 23.03.12
      integer nx_tmp, ny_tmp, ns_tmp
      real (kind=R8) ::
     *  amtol, lnorm, tchem(0:ns-1), tchee(0:ns-1), tphys(0:ns-1),
     *  tphye(0:ns-1), thevp(0:ns-1), thvpe(0:ns-1), tresn(0:ns-1),
     *  trese(0:ns-1), trfln(0:ns-1), trfle(0:ns-1), t_frac, t0
      real (kind=R8) ::                                              !srv 29.04.10 {
     &  sna0_eir(-1:nx,-1:ny,0:1,0:ns-1),
     &  smo0_eir(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0_eir(-1:nx,-1:ny,0:3), shi0_eir(-1:nx,-1:ny,0:3),         
     &  fatm_eir(-1:nx,-1:ny,0:1,0:ns-1), datm_eir(-1:nx,-1:ny,0:ns-1) !srv 29.04.10 }
      parameter (amtol=1.0e-2_R8)
      integer ncall, sput_src, sput_dst, sput_phys_col, first_flight, 
     *  in_no_of_angles, in_size_of_table, in_no_of_start_points,
     1  use_eirene, output, coreregno, coreregn2, neutrals_namelist,
     2  sput_phys_model, sput_chem_model, boundary_namelist,
     3  refl_model, reflection_on, sputter_energy_on, plate_model,
     4  plate_option, b2wall_netcdf,
     5  sput_dst2, sput_dst3, sput_frac_flag
      real (kind=R8) ::
     *  sput_phys, sput_frc, sput_chem_alpha, 
     1  sput_phys_alpha, in_dl, in_eion, alpha,
     2  plate_temp, plate_thick, therm_evap, sput_res,
     3  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     4  sput_mixed_alpha, sput_mixed_beta
      real (kind=R8) ::
     *  start_time,end_time,divpuff,mainpuff,corepuff,
     *  sum_start,sum_end,neutral_rescale,neutral_sources_rescale,
     *  sum_strti,sum_ndi,wrk2(-1:nx,-1:ny), core_sources_rescale,
     *  old_erosion(nwall,ntrack), old_deposition(nwall,ntrack)
      logical catch_up, new_sputter_namelist, sputter_namelist_used
      logical equilibrium(ntrack), over_erosion
#ifndef NO_CDF
#include <netcdf.inc>
      integer :: status, ncid
      integer, save :: ntstep_b2wall
#endif
      namelist /sputter_save/ nx_sp, ny_sp, ns_sp,
     *   sputter_yield, sputter_yield2
*   ..procedures
      external subini, subend, xertst, xerrab, sfill, ipgetr, ipgeti
      external b2xvsg, b2xvff, x05baf, smax, smin
      real (kind=R8) ::
     *  x05baf, smax, smin
      intrinsic mod, max
*   ..initialisation
      save ncall, sput_src, sput_dst, sput_frc, sput_phys, in_eion,
     *     sput_chem_alpha, sput_phys_alpha, first_flight, in_dl, 
     *     in_no_of_angles, in_size_of_table, in_no_of_start_points, 
     1     use_eirene, neutral_rescale, neutral_sources_rescale, output,
     2     coreregno, coreregn2, core_sources_rescale, alpha,
     3     sput_phys_col, neutrals_namelist, sput_phys_model,
     4     plate_temp, plate_thick, therm_evap, sput_res, 
     5     sput_chem_model,
     6     boundary_namelist, refl_model, sputter_namelist_used,
     7     reflection_on, sputter_energy_on, plate_model, plate_option,
     8     b2wall_netcdf, sput_dst2, sput_dst3,
     9     sput_frac_flag, 
     1     sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     4     sput_mixed_alpha, sput_mixed_beta, elm_count
      data ncall/0/, sput_src/1/, sput_dst/-1/, sput_frc/0.0_R8/,
     1     sput_phys/0.0_R8/, sput_chem_alpha/1.0e-2_R8/, 
     *     sput_phys_alpha/1.0e-2_R8/, first_flight/0/,
     *     in_dl/0.001_R8/, in_no_of_angles/9/, alpha/0.25_R8/,
     *     in_size_of_table/200000/, in_eion/1.0_R8/,
     1     use_eirene/0/, neutral_rescale/1.0_R8/, 
     *     neutral_sources_rescale/1.0_R8/, output/0/,
     1     coreregno/1/, coreregn2/4/, core_sources_rescale/1.0_R8/,
     2     sput_phys_col/3/, neutrals_namelist/0/,
     3     sput_phys_model/1/, plate_temp/300.0_R8/, refl_model/1/,
     4     plate_thick/0.00_R8/, therm_evap/0.0_R8/, sput_res/0.0_R8/,
     7     sput_chem_model/0/, boundary_namelist/0/,
     8     sputter_namelist_used/.false./, reflection_on/1/,
     9     sputter_energy_on/1/, plate_model/0/, plate_option/3/,
     1     b2wall_netcdf/0/, sput_dst2/-1/, sput_dst3/-1/,
     2     sput_frac_flag/0/, sput_chem_cutoff_alpha/1.0_R8/, 
     3     sput_chem_cutoff_beta/3.0_R8/,
     4     sput_mixed_alpha/1.0_R8/, sput_mixed_beta/1.0_R8/

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
        write(*,*) 'Start of b2stbr with ncall==0'
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
#ifdef USE_EIRENE
csw 25feb2014   fix IVS's version of source term avergaing
        if (nhist .gt. 1) nhist = min(max_eir_hist_slc,nhist)        !IYS 28.03.12 
#endif
        call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
        call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
        call ipgeti ('b2stbr_first_flight', first_flight)
        if(first_flight.ne.0) then
           call ipgetr ('b2stbr_first_flight_dl', in_dl)
           call ipgetr ('b2stbr_first_flight_eion', in_eion)
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
     *          in_no_of_start_points,nx,ny,in_dl,in_eion)
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
        call ipgeti ('b2stbr_sput_dst', sput_dst)
        call ipgeti ('b2stbr_sput_dst2', sput_dst2)
        call ipgeti ('b2stbr_sput_dst3', sput_dst3)
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
          call xertst (0.lt.sput_src,
     1     'faulty b2stbr_sput_src')
          call xertst (sput_dst.lt.ns,
     1     'faulty b2stbr_sput_dst')
          call xertst (sput_dst2.lt.ns,
     1     'faulty b2stbr_sput_dst2')
          call xertst (sput_dst3.lt.ns,
     1     'faulty b2stbr_sput_dst3')
          call xertst (0.0_R8.le.sput_frc,
     1     'faulty b2stbr_sput_frc')
          if(sput_dst.ge.0)
     .     write(*,'(1x,a,i2,a,i2,a,g10.3)') 
     .     'Chemical sputtering from ',sput_src,
     1     ' producing ',sput_dst,' at rate ',sput_frc
          if(sput_dst2.ge.0)
     1     write(*,*) 'Additional sputtering producing ',sput_dst2
          if(sput_dst3.ge.0)
     1     write(*,*) 'Additional sputtering producing ',sput_dst3
          call ipgetr ('b2stbr_sput_chem_alpha', sput_chem_alpha)
          call xertst (0.0_R8.le.sput_chem_alpha,
     1     'faulty b2stbr_sput_chem_alpha')
          call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
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
          call ipgeti ('b2stbr_sput_phys_col', sput_phys_col)
c the TRIM sputtering files have 9 columns corresponding to
c 0 15 30 45 55 65 75 80 85 degree impact
          call xertst (0.le.sput_phys_col.and.sput_phys_col.le.9,
     1     'faulty b2stbr_sput_phys_col')
          call ipgeti ('b2stbr_sput_phys_model', sput_phys_model)
        endif
        call ipgeti ('b2stbr_b2wall_netcdf', b2wall_netcdf)
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
     &               3.ge.plate_option, 'faulty b2stbr_plate_option')
        call ipgetr ('b2stbr_plate_temp', plate_temp)
        call xertst (0.0_R8.lt.plate_temp, 'faulty b2stbr_plate_temp')
        call ipgetr ('b2stbr_plate_thick', plate_thick)
        call xertst (0.0_R8.le.plate_thick, 'faulty b2stbr_plate_thick')
        call ipgetr ('b2stbr_therm_evap', therm_evap)
        call xertst (0.0_R8.le.therm_evap, 'faulty b2stbr_therm_evap')
        call ipgetr ('b2stbr_sput_res', sput_res)
        call xertst (0.0_R8.le.sput_res, 'faulty b2stbr_sput_res')
        call ipgetr ('b2stbr_alpha', alpha)
        call xertst (0.0_R8.le.alpha, 'faulty b2stbr_alpha')
        call ipgeti ('b2stbr_refl_model', refl_model)
        call ipgeti ('b2stbr_reflection_on', reflection_on)
        call ipgeti ('b2stbr_sputter_energy_on', sputter_energy_on)

c-xpb Reading of neutrals namelist moved here as it is needed by init_wall
        new_sputter_namelist=neutrals_namelist.eq.1.or.use_eirene.ne.0
        call read_b2mod_neutrals_namelist
     1   (nx, ny, ns, new_sputter_namelist)
        call xertst (nstrai.le.nstraid, 
     &   'Increase DEF_NSTRA in DIMENSIONS.F')
        call alloc_b2mod_neutr_src_scaling(ns,nstrai)
        if(new_sputter_namelist) call write_b2mod_neutrals_namelist
        if(neutrals_time_mod.gt.0.0_R8) then
          elm_count=tim/neutrals_time_mod
        else
          elm_count=0
        endif
        catch_up=neutrals_time_mod.gt.0.0_R8
        do while(catch_up)
          if(mod(tim,neutrals_time_mod).gt.neutrals_time_switch)
     1     then
            write(*,*) 'b2stbr: catching up ...'
            call read_b2mod_neutrals_namelist(nx, ny, ns, .true.)
            new_sputter_namelist=.true.
            call xertst (nstrai.le.nstraid, 
     &       'Increase DEF_NSTRA in DIMENSIONS.F')
            call write_b2mod_neutrals_namelist
            catch_up=neutrals_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo

*   ..read sputter namelist values
        write(*,*) 'Calling alloc_b2mod_sputter'
        call alloc_b2mod_sputter(nx,ny,ns,
     1   sput_dst,sput_dst2,sput_dst3)
        write(*,*) 'Calling alloc_b2mod_wall'
        call init_wall (nx,ny,ns,BoRiS,
     *   plate_model,plate_option,plate_temp,plate_thick)
        call xertst (ndepth_nml.gt.1.or.plate_model.eq.0, 
     .   'Inconsistent values of ndepth_nml and plate_model')
        call read_sputter_data(ns,
     1   sput_phys_col,sput_phys_model,refl_model)
        call wall_data (nx,ny,ns,BoRiS,
     *   plate_model,plate_option,plate_temp,plate_thick)
#ifdef ALLOCATE_AND_NAMELIST
        nx_tmp = nx
        ny_tmp = ny
        ns_tmp = ns
#else
        nx_tmp = DEF_NXD
        ny_tmp = DEF_NYD
        ns_tmp = DEF_NSD
#endif
        inquire(file='../b2.sputter_save.parameters',
     1   exist=sputter_namelist_used)
        if(sputter_namelist_used) then
          open(99,file='../b2.sputter_save.parameters')
          read(99,sputter_save)
          close(99)
        else
          sputter_yield =0.0_R8
          sputter_yield2=0.0_R8
        endif
        nx_sp = nx_tmp
        ny_sp = ny_tmp
        ns_sp = ns_tmp

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
        call ipgeti('b2stbr_no_b2stbr_phys_call',no_b2stbr_phys_call)
      else
        if(neutrals_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,neutrals_time_mod),
     1     neutrals_time_switch,int(tim/neutrals_time_mod),elm_count
          if(mod(tim,neutrals_time_mod).gt.neutrals_time_switch.or.
     1     int(tim/neutrals_time_mod).gt.elm_count) then
            elm_count=tim/neutrals_time_mod
            call read_b2mod_neutrals_namelist(nx, ny, ns, .true.)
            new_sputter_namelist=.true.
            call xertst (nstrai.le.nstraid, 
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
!       call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                     !srv 11.01.13
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
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)

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

csw 09oct2012 SOLPS4: explicity avoid b2stbr call here
c   --> better solution wanted here..
      if(no_b2stbr_phys_call.eq.0 ) then
       if (neutrals_namelist.eq.1 .or. use_eirene.ne.0) then
        if(ncall.eq.0) write(*,*) 'b2stbr: using b2stbr_phys'
        call b2stbr_phys (nx, ny, ns, dtim, BoRiS,
     &   vol, hx, hy, gs, qz, pbs, na, ua, kinrgy,
     &   te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni, 
     &   sna0, smo0, she0, shi0, sch0, tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &   neutral_sources_rescale, core_sources_rescale,  
     &   sput_src, sput_dst, sput_chem_alpha, sput_phys_alpha, alpha,
     &   sput_chem_model, therm_evap, sput_res, plate_model,  
     &   sputter_namelist_used, reflection_on, 
     &   sputter_energy_on, main_call, output, sput_frac_flag,
     &   sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1   sput_mixed_alpha, sput_mixed_beta,
     2   new_sputter_namelist)
       else
        if(ncall.eq.0) write(*,*) 'b2stbr: using b2stbr_bas'
        call b2stbr_bas (nx, ny, ns, dtim, BoRiS,
     &   vol, hx, hy, gs, qz, pbs, na, ua, kinrgy,
     &   te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &   sna0, smo0, she0, shi0, sch0, tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &   coreregno, coreregn2, output,
     &   neutral_sources_rescale, core_sources_rescale,  
     &   sput_src, sput_dst, sput_frc, sput_phys,
     &   sput_chem_alpha, sput_phys_alpha, alpha,
     &   sput_chem_model, therm_evap, sput_res, plate_model, 
     &   sputter_namelist_used, reflection_on, 
     &   sputter_energy_on, main_call,
     &   sput_chem_cutoff_alpha, sput_chem_cutoff_beta)
       endif
      endif

      if(first_flight.ne.0) then
        start_time=x05baf()
        if (output.ge.1) then
          open(88,file='start')
          open(89,file='end')
        endif
        do is = 0, ns-1
          if (is_neutral(is) .and. zn(is).eq.1.0_R8) then
            if (output.ge.1) then
              do iy=-1,ny
                do ix=-1,nx
                  write(88,*) ix,iy,sna0(ix,iy,0,is),shi0(ix,iy,0)
                enddo
                write(88,'()')
              enddo
            endif
            sum_start=sum(sna0(:,:,0,is))
            sum_strti=sum(shi0(:,:,0))
            divpuff=0.0_R8
            mainpuff=0.0_R8
            corepuff=0.0_R8
            wrk1=0.0_R8
            wrk2=0.0_R8
            ic=0
            do ix=-1,nx
              do ib=-1,ny
                if((leftix(ix,ib).eq.-2 .or. rightix(ix,ib).eq.nx+1)
     &            .and. region(ix,ib,0).ne.0) then
                  ic=ic+1
                  call process_b2mod_first_flight
     *             (ic,ix,ib,nx,ny,ns,nscx,iscx,is,
     *              sna0(ix,ib,0,is),shi0(ix,ib,0),wrk1,wrk2,
     *              divpuff,mainpuff,corepuff)
                  sna0(ix,ib,0,is)=0.0_R8
                  shi0(ix,ib,0)=0.0_R8
                endif
              enddo
            enddo
            do ib=-1,nx
              do iy=-1,ny
                if((topiy(ib,iy).eq.ny+1 .or. bottomiy(ib,iy).eq.-2)
     &            .and. region(ib,iy,0).ne.0) then
                  ic=ic+1
                  call process_b2mod_first_flight
     *             (ic,ib,iy,nx,ny,ns,nscx,iscx,is,
     *              sna0(ib,iy,0,is),shi0(ib,iy,0),wrk1,wrk2,
     *              divpuff,mainpuff,corepuff)
                  sna0(ib,iy,0,is)=0.0_R8
                  shi0(ib,iy,0)=0.0_R8
                endif
              enddo
            enddo
            sna0(:,:,0,is)=sna0(:,:,0,is)+wrk1(:,:)
            shi0(:,:,0)=shi0(:,:,0)+wrk2(:,:)
            if (output.gt.1) then
              do iy=-1,ny
                do ix=-1,nx
                  write(89,*) ix,iy,sna0(ix,iy,0,is),shi0(ix,iy,0)
                enddo
                write(89,'()')
              enddo
            endif
            sum_end=sum(sna0(:,:,0,is))
            sum_ndi=sum(shi0(:,:,0))
            if(main_call)
     >       write(*,*) 'first_flight: is = ',is,
     *        divpuff,mainpuff,corepuff,sum_start,sum_end,
     *        sum_start-(sum_end+divpuff+mainpuff+corepuff),
     *        sum_strti,sum_ndi
            if (sum_end.ne.0.0_R8) 
     >       sna0(:,:,0,is)=sna0(:,:,0,is)*sum_start/sum_end
            if (sum_ndi.ne.0.0_R8) 
     >       shi0(:,:,0)=shi0(:,:,0)*sum_strti/sum_ndi
          endif
        enddo
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
#ifdef USE_EIRENE
        call eirene_mc(nx, ny, ns, ismain, dtim, BoRiS, sna0_eir,     !srv 29.04.10 {
     &                 smo0_eir, she0_eir, shi0_eir, main_call)
        if (nhist .gt. 1 .and. ncall .eq. 0)                          !IYS 28.03.12 {
csw 25feb2014  fix IVS's version of source term avergaing
     &     call alloc_hist_arrays(nx,ny,ns,nhist) 

csw 25feb2014  fix IVS's version of source term avergaing
c fixed SMOOTH SOURCES FROM EIRENE , averaging procedure
c (25feb2014 different from SOLPS4.3 avergaing procedure as implemented in
c eirene_mc)
        if(nhist .gt. 1 .and. main_call) then                         
          ihist = mod(ncall, nhist)+1
          write(*,*) 'b2stbr: EIRENE source averaging ', nhist, ihist
csw 25feb2014  fix IVS's version of source term avergaing
          if(ihist .eq. 1) call clear_hist_arrays

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
csw 25feb2014  fix IVS's version of source term avergaing
            do iih=1,ihist 
              sna0_eir = sna0_eir + sna0_eir_h(:,:,:,:,iih)
              smo0_eir = smo0_eir + smo0_eir_h(:,:,:,:,iih)
              she0_eir = she0_eir + she0_eir_h(:,:,:,iih)
              shi0_eir = shi0_eir + shi0_eir_h(:,:,:,iih) 
              do ia = 1,nnatmi
                rfluxa(:,:,ia,1) = rfluxa(:,:,ia,1) 
     .                           + rfluxa_h(:,:,ia,1,iih)
              enddo
              do im = 1,nnmoli
                rfluxm(:,:,im,1) = rfluxm(:,:,im,1)
     .                           + rfluxm_h(:,:,im,1,iih)
              enddo                  
            enddo
csw 25feb2014  fix IVS's version of source term avergaing
c            sna0_eir = sna0_eir / min(dble(nhist),dble(ncall+1))
c            smo0_eir = smo0_eir / min(dble(nhist),dble(ncall+1))
c            she0_eir = she0_eir / min(dble(nhist),dble(ncall+1))
c            shi0_eir = shi0_eir / min(dble(nhist),dble(ncall+1))
c            rfluxa = rfluxa / min(dble(nhist),dble(ncall+1))                    
c            rfluxm = rfluxm / min(dble(nhist),dble(ncall+1))                    
            sna0_eir = sna0_eir / dble(ihist)
            smo0_eir = smo0_eir / dble(ihist)
            she0_eir = she0_eir / dble(ihist)
            shi0_eir = shi0_eir / dble(ihist)
            rfluxa = rfluxa / dble(ihist)                    
            rfluxm = rfluxm / dble(ihist)                    

!----  end of averaging procedure -----------------
          endif
        endif                                                            !IYS 28.03.12}
        if(ncall.le.999. and. ncall.le.ncall_iout) then
         do is = 0, ns-1
          write (chns,'(i2.2)') is
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

      do is = 0, ns-1                                                 !srv 28.07.11 20.02.12 {
       if(is.eq.ismain) then
        call b2_to_astra_sna(nx, ny, vol, sna0_eir(:,:,0,is)+
     &   sna0_eir(:,:,1,is)*na(:,:,is), svol, se)
       endif
      enddo                                                           !srv 28.07.11 20.02.12 }
      call b2_to_astra_src(nx, ny, vol, she0_eir(:,:,0)+
     &   she0_eir(:,:,1)*te(:,:)+she0_eir(:,:,2)*ne+
     &   she0_eir(:,:,3)*ne*te, svol, seh)                            !srv 11.06.13
      call b2_to_astra_src(nx, ny, vol, shi0_eir(:,:,0)+
     &   shi0_eir(:,:,1)*ti+shi0_eir(:,:,2)*ni(:,:,0)+
     &   shi0_eir(:,:,3)*ni(:,:,0)*ti, svol, sih)                     !srv 11.06.13

#endif	                       
        ia = 1    ! only base atoms
        do iy = -1,ny
         do ix = -1,nx
          fatm_eir(ix,iy,0,ia-1) = pfluxa(ix+1,iy+1,ia,1)
          fatm_eir(ix,iy,1,ia-1) = rfluxa(ix+1,iy+1,ia,1)
          datm_eir(ix,iy,ia-1)   = dab2(ix+1,iy+1,ia,1)
         enddo
        enddo
#else
        call xerrab ('b2stbr not compiled with the USE_EIRENE option')
#endif
      else
#ifndef USE_EIRENE
        if(ncall.eq.0) call eirene_mc_init             ! EIRDIAG = 0   !srv 03.04.13
#endif
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

      if(b2wall_netcdf.ne.0) then
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
        write (chns,'(i2.2)') is
#ifdef USE_EIRENE
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
#else
!	  do k=0,3
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,smo0(-1,-1,k,is),'b2stbr_smo'//chk//chns)
!	  enddo
        call my_out(70,nx,ny,b2stbr_smo(-1,-1,is),'b2stbr_smo'//chns)
#endif
!	   call my_out(70,nx,ny,ua(-1,-1,is),'b2stbr_ua'//chns)
!	   call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
       enddo                                                          !srv 29.04.10 }
      endif
      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)) then
       do is = 0, ns-1
        write (chns,'(i2.2)') is
!	  call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
#ifdef USE_EIRENE
        if(.not.is_neutral(is)) then
         do k=0,1                                                      !srv 29.04.10
           write (chk,'(i1)') k
           call my_out(70,nx,ny,sna0_eir(-1,-1,k,is),
     *               'b2stbr_sna_eir'//chk//chns)
         enddo
         call my_out(70,nx,ny,sna0_eir(:,:,0,is)+ 
     &             sna0_eir(:,:,1,is)*na(:,:,is),'b2stbr_sna_eir'//chns)
        endif
#else
!	  do k=0,1
!	   write (chk,'(i1)') k
!	   call my_out(70,nx,ny,sna0(-1,-1,k,is),'b2stbr_sna'//chk//chns)
!	  enddo
        call my_out(70,nx,ny,b2stbr_sna(-1,-1,is),'b2stbr_sna'//chns)
#endif
       enddo                                                          !srv 29.04.10
#ifdef USE_EIRENE
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
#ifdef USE_EIRENE
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
#else
#endif
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,she0(-1,-1,k),'b2stbr_she'//chk)
!	  call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbr_shi'//chk)
!	 enddo
       call my_out(70,nx,ny,b2stbr_she,'b2stbr_she')
       call my_out(70,nx,ny,b2stbr_shi,'b2stbr_shi')               !srv 29.04.10
      endif
      if (iout.eq.1 .or. (iout_b2npp7.eq.1 .or. iout_b2nppo.eq.1)) then
#ifdef USE_EIRENE
#else
!	 do k=0,3
!	  write (chk,'(i1)') k
!	  call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbr_sch'//chk)
!	 enddo
       call my_out(70,nx,ny,b2stbr_sch,'b2stbr_sch')
#endif
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

* ..finalize sputter yields
      entry b2stbr_sputter_final(nx, ny, ns)
      if(maxval(sputter_yield).gt.0.0_R8.or.
     .   maxval(sputter_yield2).gt.0.0_R8) then
        open(99,file='../b2.sputter_save.parameters')
cxpb        write(99,sputter_save)
cxpb  Replaced with explicit writing of non-zero elements only
        write(99,'(a13)') '&SPUTTER_SAVE'
        write(99,'(a8,i4,a1)') 'NX_SP = ', nx_sp, ','
        write(99,'(a8,i4,a1)') 'NY_SP = ', ny_sp, ','
        write(99,'(a8,i4,a1)') 'NS_SP = ', ns_sp, ','
        do ix = -1, nx_sp
          do iy = -1, ny_sp
            do is = 0, ns_sp-1
              if ( sputter_yield(ix,iy,is,1).ne.0.0_R8 .or.
     &             sputter_yield(ix,iy,is,2).ne.0.0_R8 .or.
     &            sputter_yield2(ix,iy,is,1).ne.0.0_R8 .or.
     &            sputter_yield2(ix,iy,is,2).ne.0.0_R8) then
                do ib = 1, 2
                  write(99,'(a14,i4,a1,i4,a1,i3,a1,i1,a2,
     &                       1pe21.13e3,
     &                       a17,i4,a1,i4,a1,i3,a1,i1,a2,
     &                       1pe21.13e3,a1)')
     .             'SPUTTER_YIELD(',ix,',',iy,',',is,',',ib,')=',
     .              sputter_yield(ix,iy,is,ib),
     .           ', SPUTTER_YIELD2(',ix,',',iy,',',is,',',ib,')=',
     .              sputter_yield2(ix,iy,is,ib),','
                end do
              end if
            end do
          end do
        end do
        write(99,'(a1)') '/'
        close(99)
      endif
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
#ifdef USE_EIRENE
c        write(*,*) 'Calling neutr_scl'
        call neutr_scl(nx, ny, ns)
c        write(*,*) 'Called neutr_scl'
#else
        call xerrab ('b2stbr not compiled with the USE_EIRENE option')
#endif
      endif
      ncall = ncall+1
      call subend ()
      return

      end subroutine b2stbr_neutr_scl

