*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr (nCv, nFc, nVx, ns,
     &  nscx, nscxmax, iscx, ismain, dtim,
     &  switch, geo, mpg, st, st_ext, st_avg, main_call)
      use b2mod_types
      use b2mod_openmp
!WG_TODO      use b2mod_wall
      use b2mod_diag
      use b2mod_time
!WG_TODO      use b2mod_layer
!WG_TODO      use b2mod_elements
!WG_TODO      use b2mod_sputter
      use b2mod_first_flight
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_balance !djm Jan2017
     & , only : b2stbr_phys_sna0to1,b2stbr_phys_smo0to3,
     &          b2stbr_phys_she0to3,b2stbr_phys_shi0to3,
#ifdef WG_TODO
     &          b2stbr_bas_sna0to1,b2stbr_bas_smo0to3,
     &          b2stbr_bas_she0to3,b2stbr_bas_shi0to3,
     &          b2stbr_first_flight_sna0to1,
     &          b2stbr_first_flight_she0to3,b2stbr_first_flight_shi0to3,
#endif
     &          balance_netcdf
      use b2mod_subsys
#ifdef B25_EIRENE
      use b2mod_eirene
     ,  , only : eirene_mc
      use eirmod_wneutrals
      use b2mod_running_average
     & , only : run_av_get_sources
#endif
      use b2us_feedback
     & , only : she_eir_tot
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_numerics_namelist
#ifndef NO_CDF
      use b2mod_tallies
#endif
#ifdef ASTRA
      use b2mod_b2_to_astra
     , , only : svol, se, seh, sih, b2_to_astra_src  !srv 28.07.11 20.02.12 11.06.13
#endif
      use b2mod_dimensions
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ad
     & , only : ncall_b2stbr, old_deposition, old_erosion
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_recycle
     & , only : int0r, int0l, int3r, int3l, int6r, int6l, int2r, int2l,
     &          int5r, int5l, int1r, int1l, int4r, int4l, fna_mol
      use b2mod_boundary_namelist
     & , only : nbcd, potpar
      use b2mod_ad
     & , only : ntstep_b2wall,
     &          my_out_folder, filename_b2w,
#ifdef WG_DEBUG
     &          my_out_sf,                                              !srv 27.06.21
#endif
     &          ncall_b2stbr_phys
      use b2mod_math
     & , only : cutlo, cutll, b2mod_math_initialised, small_r4_constant
#endif
      implicit none
*   ..input arguments
      integer nCv, nFc, nVx, ns, iCv, iFc
      integer nscx, nscxmax, iscx(0:nscxmax-1), ismain
      real (kind=R8) :: dtim
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (inout) :: mpg
      type (B2State), intent (inout) :: st
      type (B2StateExt), intent (in) :: st_ext
      type (B2Average), intent(in) :: st_avg
      logical main_call
*   ..local variables
      real (kind=R8) :: snatmp(nCv,0:1,0:ns-1), !djm Mar2017
     &                  smotmp(nCv,0:3,0:ns-1),
     &                  shetmp(nCv,0:3),
     &                  shitmp(nCv,0:3)
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
      integer is, k
#ifdef WG_TODO
      integer ib, ic, ibnd, itrack
      integer, save :: jxi, jxa, jsep
#endif
      integer, save :: no_b2stbr_phys_call = 0
      integer, save :: out_tfiles = 1                                !xpb 28.09.16
      integer, save :: ncall_iout = -1
#ifdef B25_EIRENE
      integer ia, im
#ifdef ASTRA
      integer, save :: use_astra = 0
#endif
      character ch3*3
#endif
      character chns*3, chk*1                                        !srv 11.09.09 } 23.03.12
      real (kind=R8) ::
     *  tchem(0:ns-1), tchee(0:ns-1), tphys(0:ns-1),
     *  tphye(0:ns-1), thevp(0:ns-1), thvpe(0:ns-1), tresn(0:ns-1),
     *  trese(0:ns-1), trfln(0:ns-1), trfle(0:ns-1)
#ifdef WG_TODO
     *, t_frac, t0
#endif
      real (kind=R8) ::                                              !srv 29.04.10 {
     &  sna0_eir(nCv,0:1,0:ns-1),
     &  smo0_eir(nCv,0:3,0:ns-1),
     &  she0_eir(nCv,0:3), shi0_eir(nCv,0:3), shn0_eir(nCv,0:3)    
#ifdef B25_EIRENE
     &, fatm_eir(nFc,0:1,nnatmi),                                    !srv 29.04.10 }
     &  fmol_eir(nFc,0:1,nnmoli),                                    !IYS 27.01.17
     &  flxtmp(nCv),
     &  sna_puff(nspecies),  !sna_eir_tot(nspecies),
     &  sna_core(nspecies), sna_pump(nspecies), sna_sptt(nspecies),     !iyv 07.03.18 {
     &  t(nspecies), t01(nspecies), ext_sna_tot(nspecies),
     &  sna_eir(nspecies)
      real (kind=R8), allocatable, save :: t1(:)
      integer ip, j, nplsi   
#endif
      integer sput_src,
     *  in_no_of_angles, in_size_of_table, in_no_of_start_points,
#ifdef WG_TODO
     1  coreregn2,
#endif
     2  sput_chem_model, reflection_on, sputter_energy_on
      real (kind=R8) ::
#ifdef WG_TODO
     *  pfrpuff, mainpuff, corepuff,
     *  sum_start, sum_end, lnorm,
     *  sum_strti, sum_ndi, sum_nde,
     *  work(nCv,0:1), wrk2(nCv), wrk3(nCv),
#endif
     *  shi0_ff(nCv,0:nscx-1), wrk0(nCv),
     *  f_redep(nCv,0:ns-1)
#ifdef _OPENMP
      real (kind=R8) :: start_time, end_time
#else
      real (kind=R4) :: start_time, end_time
#endif
#ifdef WG_TODO
      logical equilibrium(ntrack), over_erosion
#endif
      logical new_sputter_namelist, catch_up
#ifndef NO_CDF
#include <netcdf.inc>
      integer :: status, ncid, ndims, nvars, natts, unlimid
      integer dimids(NF_MAX_VAR_DIMS), dimsize(NF_MAX_VAR_DIMS)
      character*(NF_MAX_NAME) :: dimname(NF_MAX_VAR_DIMS)
      logical ex
      intrinsic cpu_time
#endif
*   ..procedures
      external xertst, xerrab, sfill, ipgetr, ipgeti
      external b2xvsg, b2xvff, smax, smin, get_jsep
      real (kind=R8) :: smax, smin
      intrinsic mod
#if ( defined(_OPENMP) && defined(CRAY) )
      real (kind=R8) :: omp_get_wtime
      external omp_get_wtime
#endif
*   ..initialisation
      save sput_src,
     *     in_no_of_angles, in_size_of_table, in_no_of_start_points,
#ifdef WG_TODO
     1     coreregn2,
#endif
     5     sput_chem_model, reflection_on, sputter_energy_on
      data sput_src/1/,
     *     in_no_of_angles/9/,
     *     in_size_of_table/200000/,
#ifdef WG_TODO
     1     coreregn2/4/,
#endif
     7     sput_chem_model/0/,
     8     reflection_on/1/, sputter_energy_on/1/

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
      if (ncall_b2stbr.eq.0) then
!WG_TODO        call get_jsep(nx,ny,jxi,jxa,jsep)
        if (switch%use_eirene.ne.0) out_tfiles=0
#if ( defined(B25_EIRENE) && defined(ASTRA) )
        call ipgeti ('b2mndr_astra', use_astra)
#endif
        call ipgeti ('b2stbr_out_tfiles', out_tfiles)                !srv 07.04.09
        call ipgeti ('b2stbr_ncall_iout', ncall_iout)                !srv 23.03.12
        if(switch%b2stbr_first_flight.ne.0) then
          call ipgeti ('b2stbr_first_flight_no_of_flights',
     *      in_no_of_angles)
          call ipgeti ('b2stbr_first_flight_table_size',
     *      in_size_of_table)
          in_no_of_start_points = nCv - mpg%nCi
          call ipgeti ('b2stbr_first_flight_no_of_start_points',
     *      in_no_of_start_points)
          call xertst (0.lt.in_no_of_start_points,
     *     'faulty input first_flight_no_of_start_points')
#ifdef _OPENMP
          start_time = omp_get_wtime()
#else
          call cpu_time(start_time)
#endif
          call alloc_b2mod_first_flight(
     *      in_no_of_angles,in_size_of_table,
     *      in_no_of_start_points,switch%b2stbr_in_dl)
          do iCv=mpg%nCi+1,nCv
            iFc=mpg%cvFc(mpg%cvFcP(iCv,1))
            call setup_b2mod_first_flight(mpg,geo,iCv,iFc,
     *       switch%b2stbr_output)
          enddo
#ifdef _OPENMP
          end_time = omp_get_wtime()
#else
          call cpu_time(end_time)
#endif
          write(*,'(a,i8,a,i8,a,g10.3,a,g12.5,a)')
     *     'Used ',ptr_table,' out of ',size_of_table,
     *     ' (',real(ptr_table)/real(size_of_table)*100,'%)',
     *      end_time-start_time,' cpu-seconds'
        endif
        call ipgeti ('b2stbr_sput_src', sput_src)
        call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
        if(switch%b2stbr_sput_phys.ne.0.0_R8) then
#ifdef WG_TODO
          call xertst (sput_dst.lt.ns,
     1     'faulty b2stbr_sput_dst')
#endif
        endif
        call ipgeti ('b2stbr_reflection_on', reflection_on)
        call ipgeti ('b2stbr_sputter_energy_on', sputter_energy_on)
        if(switch%use_eirene.eq.0) switch%neutral_rescale=1.0_R8
        if(switch%use_eirene.eq.0) switch%neutral_sources_rescale=1.0_R8
#ifdef WG_TODO
        call xertst (cbirso.le.switch%b2stbr_coreregno.and.
     1   switch%b2stbr_coreregno.le.cbirso+cbnrso-1,
     2   'faulty internal parameter coreregno')
        if (nncut.ge.2) then
          call ipgeti ('b2stbc_coreregn2', coreregn2)
          call xertst (cbirso.le.coreregn2.and.
     1     coreregn2.le.cbirso+cbnrso-1,
     2     'faulty internal parameter coreregn2')
        endif
#endif
        call ipgeti ('b2stbr_no_b2stbr_phys_call', no_b2stbr_phys_call)
      else
        if(neutrals_time_mod.gt.0.0_R8) then
          write(*,*) 'NEUTRALS: ',tim,mod(tim,neutrals_time_mod),
     1     neutrals_time_switch,int(tim/neutrals_time_mod),
     2     neutrals_count
          catch_up=
     1    (mod(tim,neutrals_time_mod).ge.neutrals_time_switch.and.
     1     neutrals_time_switch.gt.0.0_R8).or.
     2     int(tim/neutrals_time_mod).gt.neutrals_count
        else if (neutrals_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.neutrals_time_switch
          if(catch_up) write(*,'(a,a)')
     .     'Switching to next neutrals namelist file ',
     .      trim(neutrals_filename)
        else
          catch_up=.false.
        end if
        if(catch_up) then
          if(neutrals_time_mod.gt.0.0_R8)
     >     neutrals_count=int(tim/neutrals_time_mod)
          call read_neutrals_namelist(ns, mpg, switch, .true.)
          new_sputter_namelist=.true.
          call xertst (nstrat.le.nstraid,
     &     'Increase DEF_NSTRA in b2mod_dimensions')
          call write_b2mod_neutrals_namelist
        endif
      endif
*   ..test nCv, nFc, ns
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     & 'faulty argument ismain')
*   ..extensive tests on first few calls
      if(ncall_b2stbr.lt.3) then
*    ..test sign of vol
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
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
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (nCv*ns, st%pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv*2, st%dv%ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, st%dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, st%pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, st%pl%ti, 1, 'ti', '.gt.')
      endif
*   ..initialize sources to 0
      call sfill (nCv*2*ns, 0.0_R8, st%srw%sna0, 1)
      call sfill (nCv*4*ns, 0.0_R8, st%srw%smo0, 1)
      call sfill (nCv*2, 0.0_R8, st%srw%sne0, 1)
      call sfill (nCv*4, 0.0_R8, st%srw%she0, 1)
      call sfill (nCv*4, 0.0_R8, st%srw%shi0, 1)
      call sfill (nCv*4, 0.0_R8, st%srw%shn0, 1)
      call sfill (nCv*4, 0.0_R8, st%srw%sch0, 1)
      call sfill (nCv*nscx, 0.0_R8, shi0_ff, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_eir, 1)                    !srv 29.04.10 {
      call sfill (nCv*4*ns, 0.0_R8, smo0_eir, 1)
      call sfill (nCv*4, 0.0_R8, she0_eir, 1)
      call sfill (nCv*4, 0.0_R8, shi0_eir, 1)
      call sfill (nCv*4, 0.0_R8, shn0_eir, 1)
#ifdef B25_EIRENE
      fatm_eir = 0.0_R8
#endif

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
#ifdef B25_EIRENE
      if( sum(sna_corr) .gt. 0.0_R8 .and. .not.allocated(t1)) then
        allocate(t1(1:nspecies))
        t1=0.0_R8
      endif 
#endif
      if (main_call) then
        old_erosion(1:nbnd,1:ntargsp)=erosion(1:nbnd,1:ntargsp)
        old_deposition(1:nbnd,1:ntargsp)=deposition(1:nbnd,1:ntargsp)
      endif

csw 09oct2012 SOLPS4: explicitly avoid b2stbr call here
c   --> better solution wanted here..
      if(no_b2stbr_phys_call.eq.0 ) then
       if(switch%b2stbr_redep_alpha.eq.0.0_R8) then
        f_redep=0.0_R8
       else
        call xerrab ('WG - redep_alpha .ne. 0')
#ifdef WG_TODO
        do is=0, ns-1
          if(is_neutral(is).and.lnext(is,is+1)) then
            do iy=-1,ny
              do ix=-1,nx
                f_redep(ix,iy,is) = 1.0_R8 /
     1           ( 1.0_R8 + switch%b2stbr_redep_alpha*
     2           (qe*bb(ix,iy,3)/(am(is)*mp)/
     3           (rsa(ix,iy,is)*ne(ix,iy)))**2 )
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
#endif
       endif
       if (switch%b2stbr_neutrals_namelist.eq.1
     &                          .or. switch%use_eirene.ne.0) then
        if(ncall_b2stbr.eq.0) write(*,*) 'b2stbr: using b2stbr_phys'
        if (balance_netcdf.ne.0) then
          !djm Mar2017 save sources before hand
          snatmp=st%srw%sna0
          smotmp=st%srw%smo0
          shetmp=st%srw%she0
          shitmp=st%srw%shi0
        endif
        call b2stbr_phys (nCv, nFc, nVx, ns, nscx, nscxmax, iscx,
     &   dtim, switch, geo, mpg,
     &   st%pl, st%dv, st%co, st%rt, st%rtw, st_ext, st%srw,
     &   tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle,
     &   sput_src, sput_chem_model, reflection_on, sputter_energy_on,
     &   main_call, new_sputter_namelist, shi0_ff, f_redep)
        if (balance_netcdf.ne.0) then
          !djm Mar2017 store source contributions
          b2stbr_phys_sna0to1=st%srw%sna0-snatmp
          b2stbr_phys_smo0to3=st%srw%smo0-smotmp
          b2stbr_phys_she0to3=st%srw%she0-shetmp
          b2stbr_phys_shi0to3=st%srw%shi0-shitmp
        endif
       else
#ifdef WG_TODO
        if(ncall_b2stbr.eq.0) write(*,*) 'b2stbr: using b2stbr_bas'
        if (balance_netcdf.ne.0) then
          !djm Mar2017 save sources beforehand
          snatmp=sna0
          smotmp=smo0
          shetmp=she0
          shitmp=shi0
        endif
        call b2stbr_bas (nx, ny, ns, nscx, iscx, dtim, switch%BoRiS,
     &   na, ua, kinrgy,
     &   te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &   sna0, smo0, she0, shi0, sch0, tchem, tchee, tphys, tphye,
     &   thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &   switch%b2stbr_coreregno, coreregn2, output,
     &   switch%neutral_sources_rescale,
     &   switch%b2stbr_core_sources_rescale,
     &   sput_src, switch%b2stbr_sput_frc, switch%b2stbr_sput_phys,
     &   switch%b2stbr_sput_chem_alpha, switch%sput_phys_alpha,
     &   switch%b2stbr_alpha, sput_chem_model, switch%b2stbr_therm_evap,
     &   switch%b2stbr_sput_res, reflection_on, sputter_energy_on,
     &   main_call, switch%b2stbr_sput_chem_cutoff_alpha,
     &   switch%b2stbr_sput_chem_cutoff_beta, shi0_ff,
     &   f_redep, switch%b2stbr_potential_at_guard_cell,
     &   switch%recycled_neutrals_contr)
        if (balance_netcdf.ne.0) then
          !djm Mar2017 store source contributions
          b2stbr_bas_sna0to1=sna0-snatmp
          b2stbr_bas_smo0to3=smo0-smotmp
          b2stbr_bas_she0to3=she0-shetmp
          b2stbr_bas_shi0to3=shi0-shitmp
        endif
#endif
       endif
      endif

      if(switch%b2stbr_first_flight.ne.0) then
#ifdef WG_TODO
        if (balance_netcdf.ne.0) then
          !djm Mar2017 save sources beforehand
          snatmp=sna0
          shetmp=she0
          shitmp=shi0
        endif
#ifdef _OPENMP
        start_time = omp_get_wtime()
#else
        call cpu_time(start_time)
#endif
        if (output.ge.1) then
          open(88,file='start')
          open(89,file='end')
        endif
        do is = 0, nscx-1
          if (output.ge.1) then
            do iCv=1,nCv
              write(88,*) iCv,sna0(iCv,0,iscx(is)),
     .                        shi0(iCv,0),she0(iCv,0)
            enddo
            write(88,'()')
          endif
          sum_start=0.0_R8
          sum_strti=sum(shi0_ff(:,is))
          pfrpuff=0.0_R8
          mainpuff=0.0_R8
          corepuff=0.0_R8
          work=0.0_R8
          wrk2=0.0_R8
          wrk3=0.0_R8
          ic=0
          do iCv=1,nCv
!TODO the loop must over all edge cells
            if((leftix(ix,ib).eq.-2 .or. rightix(ix,ib).eq.nx+1)
     &        .and. region(ix,ib,0).ne.0) then
              ic=ic+1
              call process_b2mod_first_flight
     *         (mpg,st,ic,iCv,ns,iscx(is),
     *          dtim,sna0(iCv,0,iscx(is)),shi0_ff(iCv,is),
     *          work,wrk2,wrk3,
     *          pfrpuff,mainpuff,corepuff)
              shi0(iCv,0)=shi0(iCv,0)-shi0_ff(iCv,is)
              sum_start=sum_start+sna0(iCv,0,iscx(is))
              sna0(iCv,0,iscx(is))=0.0_R8
              shi0_ff(iCv,is)=0.0_R8
            endif
          enddo
          sum_end=sum(work(:,:))
          sum_ndi=sum(wrk2(:))
          sum_nde=sum(wrk3(:))
          if(main_call)
     >     write(*,*) 'first_flight: is = ',is,
     *      pfrpuff,mainpuff,corepuff,sum_start,sum_end,
     *      sum_start-(sum_end+pfrpuff+mainpuff+corepuff),
     *      sum_strti,sum_ndi,sum_nde
          if (sum_end.ne.0.0_R8) then
            sna0(:,0,iscx(is))=sna0(:,0,iscx(is))
     .       +work(:,0)*sum_start/sum_end
            sna0(:,0,iscx(is)+1)=sna0(:,0,iscx(is)+1)
     .       +work(:,:,1)*sum_start/sum_end
            she0(:,0)=she0(:,0)+wrk3(:)*sum_start/sum_end
          endif
          if (sum_ndi.ne.0.0_R8)
     >     shi0(:,0)=shi0(:,0)+wrk2(:)*sum_strti/sum_ndi
          if (output.gt.1) then
            do iCv=1,nCv
              write(89,*) iCv,sna0(iCv,0,iscx(is))
     .                       ,shi0(iCv,0),she0(iCv,0)
            enddo
            write(89,'()')
          endif
          if (main_call) then
            write(*,*) 'Total first-flight-source na:', iscx(is),
     &       sum(sna0(:,0,iscx(is)))
            write(*,*) 'Total first-flight-source na:', iscx(is)+1,
     &       sum(sna0(:,0,iscx(is)+1))
          endif
        enddo
        if (main_call) then
          write (*,*)
     >     'Total first-flight-source ion power=', sum(shi0(:,0))
          write (*,*)
     >     'Total first-flight-source electron power=', sum(she0(:,0))
        endif
        if (output.ge.1) then
          close(88)
          close(89)
        endif
#ifdef _OPENMP
        end_time = omp_get_wtime()
#else
        call cpu_time(end_time)
#endif
        if(main_call)
     >   write(*,*) end_time-start_time,' cpu-seconds'
        if (balance_netcdf.ne.0) then
          !djm Mar2017 store source contributions
          b2stbr_first_flight_sna0to1=sna0-snatmp
          b2stbr_first_flight_she0to3=she0-shetmp
          b2stbr_first_flight_shi0to3=shi0-shitmp
        endif
#endif
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

      if(switch%use_eirene.ne.0) then
#ifdef B25_EIRENE
        call eirene_mc(nCv, nFc, ns, dtim,                          !srv 29.04.10 {
     &                 sna0_eir, smo0_eir, she0_eir, shi0_eir, shn0_eir,
     &                 main_call, switch, geo, mpg, st, st_ext)

!!!!!!
*   ..main part of implementation of method of effective sources (ES)
        if (sum(sna_corr) .gt. 0.0_R8 ) then                             !iyv 07.03.18 {
          ! if ( iout .ne. 0 )
           write(*,*) '******* sna_eir correction calculations *******'
          sna_core = 0.0_R8
          do is=0,st%pl%ns-1
            do iCv=mpg%nCi+1, nCv
              k = b2espcr(is)
              if (is_neutral(is) .or. (.not. do_sna_corr_core(k))) cycle
              if(mpg%cvOnClosedSurface(iCv)) then
                iFc = mpg%cvFc(mpg%cvFcP(iCv,1))  
                if (switch%mdf_fnb.ne.0) then
                  sna_core(k) = sna_core(k) + st%dv%fna_mdf(iFc,1,is)
                else
                  sna_core(k) = sna_core(k) + st%dv%fna(iFc,1,is)
                end if
              endif
            enddo
          enddo

          sna_eir= 0.0_R8
          ext_sna_tot= 0.0_R8
          print *, "max sna0_eir", maxval(abs(sna0_eir)),
     &    size(sna0_eir,1), size(sna0_eir,2),size(sna0_eir,3)
          do is=0,st%pl%ns-1
            do iCv=1,nCv
              if( is_neutral(is) ) cycle
              k = b2espcr(is)
              sna_eir(k) = sna_eir(k) +
     &       sna0_eir(iCv,0,is) + sna0_eir(iCv,1,is)*st%pl%na(iCv,is)
              ext_sna_tot(k) = ext_sna_tot(k) + st_ext%sna(iCv,is)
            enddo
          enddo

          sna_puff =0.0_R8
          do is=1,nstrat
            do ia=1,nnatmi
              if(latmscl(ia).eq.0) cycle
              if(crcstra(is).eq.'C') then
                if(gpfc(ia,is) .ne. 0.0_R8) sna_puff(latmscl(ia)) =
     &           sna_puff(latmscl(ia)) + userfluxparm(is,1)*gpfc(ia,is) !*recyctmp
               end if
            end do
          enddo
          
          sna_pump = 0.0_R8
          sna_sptt = 0.0_R8
          do ia = 1, nnatmi
            k = latmscl(ia)
            if (k.eq.0) cycle
            sna_sptt(latmscl(ia)) = sna_sptt(latmscl(ia)) + eirspta(ia)
            sna_pump(latmscl(ia)) = sna_pump(latmscl(ia)) + eirpump(ia)
            do im = 1, nnmoli
              if (mlcmp(ia,im).gt.0) then
                sna_pump(latmscl(ia)) = sna_pump(latmscl(ia)) +
     &             mlcmp(ia,im) * eirpump(nnatmi+im)
                sna_sptt(latmscl(ia)) = sna_sptt(latmscl(ia)) +
     &             mlcmp(ia,im) * eirsptm(im)
              endif
            enddo
          end do
          do ip = 1, nnplsi
            is = latmscl(lkindp(ip))
            if (is.eq.0) cycle
            sna_pump(is) = sna_pump(is) +
     &                         eirpump(nnatmi+nnmoli+nnioni+ip)
          enddo


!           if ( iout .ne. 0 ) then
            write(*,'(a5,7a16)') 'ia', 'sna_core', 'sna_eir',
     &  'sna_puff', 'sna_pump', 'sna_sptt', 'sna_ext', 't1'
            do k = 1, nspecies
              write(*, '(i5, 7es16.7)') k, sna_core(k),  sna_eir(k),
     &      sna_puff(k), sna_pump(k), sna_sptt(k), ext_sna_tot(k), t1(k)
            enddo

            write(*,'(a3,3a17)') 'ia','dG/sna_eir','t1(ia)','t(ia)'
!           endif

          do k = 1, nspecies
            if (sna_eir(k).ne.0.0_R8) then
              t01(k) =  sna_corr(k)*( sna_puff(k) + sna_core(k) +
     &            sna_sptt(k) - sna_pump(k) + ext_sna_tot(k) ) /
     &            sna_eir(k)
            else
              t01(k) = 0.0_R8
            end if
            t1(k) = t01(k)* ( 1.0_R8 - exp(-0.02) ) + !c_pumpavr
     &                     exp(-0.02)*t1(k) !c_pumpavr
 !           t1(k) = t01(k)
            t(k) = 1.0_R8 +
     &              max( min( taumax(k),t1(k) ), -taumax(k) )
            !if ( iout .ne. 0 )
             write(*,'(i3,3es17.4)') k, t01(k), t1(k), t(k)
!            sna_eir(k) = sna_eir(k)*t(k)
          enddo

          sna_eir = 0.0_R8
          do is = 0, ns-1
            k = b2espcr(is)
            if( is_neutral(is) ) cycle
!            if( eirspta(k) .eq. 0.0_R8) then  ! not for sputtered particles yet
            do iCv=1,nCv
              if (do_sna_corr_core(k) .or. 
     &                (.not. mpg%cvOnClosedSurface(iCv))) then
                sna0_eir(iCv,0,is) = sna0_eir(iCv,0,is)*t(k)
                sna0_eir(iCv,1,is) = sna0_eir(iCv,1,is)*t(k)
              end if
              sna_eir(k) = sna_eir(k) +
     &           sna0_eir(iCv,0,is) +
     &           sna0_eir(iCv,1,is)*st%pl%na(iCv,is)
            enddo
          end do

        end if ! if b2stbr_sna_corr ...                                   !iyv 07.03.18 }
!!!!!!
        if(ncall_b2stbr.le.999. and. ncall_b2stbr.le.ncall_iout) then
         do is = 0, ns-1
          write (chns,'(i3.3)') is
          write (ch3 ,'(i3.3)') ncall_b2stbr
          if(.not.is_neutral(is)) then
           wrk0=sna0_eir(:,0,is)+
     &          sna0_eir(:,1,is)*st%pl%na(:,is)
           call my_out_us(70,nCv,0,wrk0,'b2stbr_sna_eir'//chns//
     &                    '_'//ch3)
          endif
         enddo
        endif

*       compute total Eirene sources for averaging
cwdk only total Eirene source are averaged!
cwdk note: with new variables sna0_eir_tot, smo0_eir_tot,
cwdk       she0_eir_tot and shi0_eir_tot now available,
cwdk       might want to avoid recomputation of several wrk0 arrays further
cwdk       down in this routine
        do is = 0, ns-1
          st%srw%sna0_eir_tot(:,is) = sna0_eir(:,0,is) +
     &      sna0_eir(:,1,is)*st%pl%na(:,is)
          st%srw%smo0_eir_tot(:,is) = smo0_eir(:,0,is) +
     &      smo0_eir(:,1,is)*st%pl%ua(:,is) +
     &      smo0_eir(:,2,is)*am(is)*mp*st%pl%na(:,is) +
     &      smo0_eir(:,3,is)*am(is)*mp*st%pl%na(:,is)*st%pl%ua(:,is)
        enddo
        st%srw%she0_eir_tot = she0_eir(:,0) + she0_eir(:,1)*st%pl%te +
     &    she0_eir(:,2)*st%dv%ne + she0_eir(:,3)*st%dv%ne*st%pl%te
        st%srw%shi0_eir_tot = shi0_eir(:,0) + shi0_eir(:,1)*st%pl%ti +
     &    shi0_eir(:,2)*st%dv%ni(:,0) +
     &    shi0_eir(:,3)*st%dv%ni(:,0)*st%pl%ti
        st%srw%shn0_eir_tot = shn0_eir(:,0) + shn0_eir(:,1)*st%pl%tn +
     1    shn0_eir(:,2)*st%dv%nn + 
     1    shn0_eir(:,3)*st%pl%tn*st%dv%nn

*       get averaged sources
        if (switch%ird_aver > 0) then

cwdk if ird_aver > 0, all sources computed by Eirene are overwritten (!!!),
cwdk i.e. fully replaced by the averaged values from b2faveri
          sna0_eir = 0.0_R8
          smo0_eir = 0.0_R8
          she0_eir = 0.0_R8
          shi0_eir = 0.0_R8

cwdk get the averaged sources
          call run_av_get_sources(
     &     st%srw%sna0_eir_tot, st%srw%smo0_eir_tot,
     &     st%srw%she0_eir_tot, st%srw%shi0_eir_tot,
     &     st%srw%shn0_eir_tot,
     &     nCv, ns, st_avg)

cwdk note: only total Eirene source are averaged!
cwdk these are now put back into the constant part of the sources (0-component)
          do is = 0, ns-1
            sna0_eir(:,0,is) = st%srw%sna0_eir_tot(:,is)
            smo0_eir(:,0,is) = st%srw%smo0_eir_tot(:,is)
          enddo
          she0_eir(:,0) = st%srw%she0_eir_tot
          shi0_eir(:,0) = st%srw%shi0_eir_tot
          shn0_eir(:,0) = st%srw%shn0_eir_tot
        endif

        st%srw%sna0 = st%srw%sna0 + sna0_eir
        st%srw%smo0 = st%srw%smo0 + smo0_eir
        st%srw%she0 = st%srw%she0 + she0_eir
        st%srw%shi0 = st%srw%shi0 + shi0_eir
        if (switch%tn_style.eq.0) then     !! default, combined ion-neutral energy equation
          st%srw%shi0 = st%srw%shi0 + shn0_eir 
        elseif (switch%tn_style.eq.1) then !! pure ion energy equation
          !! combination of condensation reactions (leading to non-zero shn from eirene) and tn_style = 1 is not allowed.
          call xertst (SUM(shn0_eir).eq.0,
     &     'shn_eir.ne.0 is not allowed for tn_style = 1')
        else                               !! separate neutral energy equation
          st%srw%shn0 = st%srw%shn0 + shn0_eir 
        endif

#ifdef ASTRA
        if (use_astra.ne.0) then
          do is = 0, ns-1                                                !srv 28.07.11 20.02.12 {
            if(is.eq.ismain) then
              wrk0=sna0_eir(:,0,is)+
     &             sna0_eir(:,1,is)*st%pl%na(:,is)
              call b2_to_astra_src(mpg, ny, geo%cvVol, wrk0, svol, se)
            endif
          enddo                                                          !srv 28.07.11 20.02.12 }
          wrk0=she0_eir(:,:,0)+
     &     she0_eir(:,:,1)*st%pl%te(:,:)+she0_eir(:,:,2)*st%dv%ne+
     &     she0_eir(:,:,3)*st%dv%ne*st%pl%te
          call b2_to_astra_src(mpg, ny, vol, wrk0, svol, seh)             !srv 11.06.13
          wrk0=shi0_eir(:,:,0)+
     &     shi0_eir(:,:,1)*st%pl%ti+shi0_eir(:,:,2)*st%dv%ni(:,:,0)+
     &     shi0_eir(:,:,3)*st%dv%ni(:,:,0)*st%pl%ti
          call b2_to_astra_src(mpg, ny, vol, wrk0, svol, sih)             !srv 11.06.13
        endif
#endif
        do ia = 1, nnatmi
! IYS 27.01.2017
! interpolate atom fluxes to cell faces (for output only)
         flxtmp(:) = pfluxa(:,ia,1)
         call intface(nCv,nFc,mpg%fcCv,geo%fcVol,
     &                flxtmp,fatm_eir(:,0,ia))
         flxtmp(:) = rfluxa(:,ia,1)
         call intface(nCv,nFc,mpg%fcCv,geo%fcVol,
     &                flxtmp,fatm_eir(:,1,ia))
! IYS 27.01.2017
        enddo
! IYS 27.01.2017
        do im = 1, nnmoli
! interpolate molecule fluxes to cell faces (for output only)
         flxtmp(:) = pfluxm(:,im,1)
         call intface(nCv,nFc,mpg%fcCv,geo%fcVol,
     &                flxtmp,fmol_eir(:,0,im))
         flxtmp(:) = rfluxm(:,im,1)
         call intface(nCv,nFc,mpg%fcCv,geo%fcVol,
     &                flxtmp,fmol_eir(:,1,im))
        enddo
! IYS 27.01.2017

#else
        call xerrab ('b2stbr not compiled with the B25_EIRENE option')
#endif
      else
cwdk remove this if-test, to allow runs without Eirene if code is compiled with Eirene
        call eirene_mc_init(nCv,nFc)             ! EIRDIAG = 0   !srv 03.04.13
      endif

* ..sources
      st%srw%b2stbr_sna(:,0:ns-1)=
     1 st%srw%sna0(:,0,0:ns-1)+
     2 st%srw%sna0(:,1,0:ns-1)*st%pl%na(:,0:ns-1)
      do is=0,ns-1
        st%srw%b2stbr_smo(:,is)=
     1   st%srw%smo0(:,0,is)+
     2   st%srw%smo0(:,1,is)*st%pl%ua(:,is)+
     3   st%srw%smo0(:,2,is)*am(is)*mp*st%pl%na(:,is)+
     4   st%srw%smo0(:,3,is)*am(is)*mp*st%pl%na(:,is)*st%pl%ua(:,is)
      enddo
      st%srw%b2stbr_she(:)=
     1 st%srw%she0(:,0)+
     2 st%srw%she0(:,1)*st%pl%te(:)+
     3 st%srw%she0(:,2)*st%dv%ne(:)+
     4 st%srw%she0(:,3)*st%dv%ne(:)*st%pl%te(:)
      st%srw%b2stbr_shi(:)=
     1 st%srw%shi0(:,0)+
     2 st%srw%shi0(:,1)*st%pl%ti(:)+
     3 st%srw%shi0(:,2)*st%dv%ni(:,0)+
     4 st%srw%shi0(:,3)*st%dv%ni(:,0)*st%pl%ti(:)
      st%srw%b2stbr_shn(:)=
     1 st%srw%shn0(:,0)+
     2 st%srw%shn0(:,1)*st%pl%ti(:)+
     3 st%srw%shn0(:,2)*st%dv%nn(:)+
     4 st%srw%shn0(:,3)*st%dv%nn(:)*st%pl%tn(:)
      st%srw%b2stbr_sch(:)=
     1 st%srw%sch0(:,0)+
     2 st%srw%sch0(:,1)*st%pl%po(:)+
     3 st%srw%sch0(:,2)*st%dv%ne(:)+
     4 st%srw%sch0(:,3)*st%dv%ne(:)*st%pl%po(:)
      st%srw%b2stbr_skt(:)=
     1 st%srw%skt0(:,0)+
     2 st%srw%skt0(:,1)*st%pl%kt(:)+
     3 st%srw%skt0(:,2)*st%dv%ni(:,1)+
     4 st%srw%skt0(:,3)*st%dv%ni(:,1)*st%pl%kt(:)
      st%srw%b2stbr_szt(:)=
     1 st%srw%szt0(:,0)+
     2 st%srw%szt0(:,1)*st%pl%zt(:)+
     3 st%srw%szt0(:,2)*st%dv%ni(:,1)+
     4 st%srw%szt0(:,3)*st%dv%ni(:,1)*st%pl%zt(:)
      st%srw%b2stbr_sne(:)=
     1 st%srw%sne0(:,0)+
     2 st%srw%sne0(:,1)*st%dv%ne(:)

#ifdef WG_TODO
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
cxpb  This only has meaning if we are following multiple mixed species and
cxpb  this is not the first timestep
              t0 = 0.0_R8
              equilibrium(itrack) =
     &         ntargsp.gt.1 .and. ncall_b2stbr.gt.0 .and.
     &         switch%b2stbr_sput_frac_flag.ne.0
            endif
            if (t0.lt.0.0_R8) then  ! we have over-eroded, so we correct
              if (itrack.ne.itrack_0(ibnd)) then  ! over-erosion of a deposited species
                t_frac = min(t_frac,
     1            particle_inventory(ibnd,itrack) /
     2           (particle_inventory(ibnd,itrack) - t0))
                over_erosion = .true.
              else  ! we are keeping track of the bulk and eroded it all
cxpb  sput_frac_flag.eq.0 means we are not following mixed materials
cxpb  plate_thickness.eq.0 means we are not doing full inventories
                call xertst (switch%b2stbr_sput_frac_flag.eq.0 .or.
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
                if (switch%b2stbr_sput_frac_flag.ne.0)
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
#endif
      if(switch%b2stbr_b2wall_netcdf.ne.0 .and. ntargsp.gt.0) then
#ifndef NO_CDF
        if(ncall_b2stbr.eq.0) then
          filename_b2w='b2wall.nc'
          call find_file(filename_b2w,ex)
          if (.not.ex.or.switch%b2mndr_stim.ge.0.0_R8) then
            ntstep_b2wall = 0
            write(6,'(a)') trim(filename_b2w)//' will be created'
            call define_cdf_b2mod_wall(filename_b2w,ncid,nbnd,ntargsp)
          else if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
            status=nf_open(filename_b2w,or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(status)
            status=nf_inq(ncid,ndims,nvars,natts,unlimid)
            call check_cdf_status(status)
            status=nf_inq_dim(ncid,dimids(1),dimname(1),dimsize(1))
            call check_cdf_status(status)
            ntstep_b2wall=dimsize(1)
            status=nf_close(ncid)
            call check_cdf_status(status)
            write(6,'(a)') trim(filename_b2w)//' will be appended'
          else
            ntstep_b2wall=0
            write(6,'(a)') trim(filename_b2w)//' will be replaced'
            call define_cdf_b2mod_wall(filename_b2w,ncid,nbnd,ntargsp)
          end if
          write(6,'(a,i8)') trim(filename_b2w)//' : ntstep_b2wall = ',
     .                                          ntstep_b2wall
          if (ntstep_b2wall.eq.0) then
            status=nf_open(filename_b2w,or(NF_WRITE,NF_SHARE),ncid)
            call check_cdf_status(status)
            write(*,*) 'Writing track_species'
            call write_cdf_char_variable(ncid,-1,
     *       'track_species',track_species(1:ntargsp))
            write(*,*) 'Writing plate_area'
            call write_cdf_variable(ncid,-1,
     *       'plate_area',plate_area(1:nbnd))
            status=nf_close(ncid)
            call check_cdf_status(status)
            write(6,'(a)') 'Closed '//trim(filename_b2w)
          end if
        endif
        if(main_call.and.
     *     mod(ncall_b2stbr,switch%b2stbr_b2wall_netcdf).eq.0) then
          ntstep_b2wall=ntstep_b2wall+1
          status=nf_open(filename_b2w,or(NF_WRITE,NF_SHARE),ncid)
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

      if (switch%b2stbr_iout.eq.1 .or. switch%b2npmo_iout.eq.1
     &                            .or. switch%iout_b2wdat.eq.4) then   !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if (switch%use_eirene.ne.0) then
         if(.not.is_neutral(is) .and. switch%iout_b2wdat.eq.4) then
!          do k=0,1 ! 3                                                  !srv 29.04.10 {
!           write (chk,'(i1)') k
!           call my_out(70,nx,ny,smo0_eir(-1,-1,k,is),
!     *              'b2stbr_smo_eir'//chk//chns)
!          enddo
          wrk0=smo0_eir(:,0,is)+
     &         smo0_eir(:,1,is)*st%pl%ua(:,is)+
     &         smo0_eir(:,2,is)*st%pl%na(:,is)*mp*am(is)+
     &         smo0_eir(:,3,is)*st%pl%ua(:,is)*st%pl%na(:,is)*mp*am(is)
          call my_out_us(70,nCv,0,wrk0,'b2stbr_smo_eir'//chns)
         endif
        else
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,smo0(-1,-1,k,is),'b2stbr_smo'//chk//chns)
!        enddo
         call my_out_us(70,nCv,0,st%srw%b2stbr_smo(1,is),
     &                  'b2stbr_smo'//chns)
        endif
!       call my_out(70,nx,ny,ua(-1,-1,is),'b2stbr_ua'//chns)
!       call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
       enddo                                                          !srv 29.04.10 }
      endif
      if (switch%b2stbr_iout.eq.1 .or.
     &    switch%b2npco_iout.eq.1 .or. switch%iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!       call my_out(70,nx,ny,na(-1,-1,is),'b2stbr_na'//chns)
        if(switch%use_eirene.ne.0) then
         if(.not.is_neutral(is)) then
          if (switch%iout_b2wdat.ne.4) then
           do k=0,1                                                   !srv 29.04.10
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,sna0_eir(1,k,is),
     *               'b2stbr_sna_eir'//chk//chns)
           enddo
          endif
          wrk0=sna0_eir(:,0,is)+sna0_eir(:,1,is)*st%pl%na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2stbr_sna_eir'//chns)
         endif
        else
!        do k=0,1
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,sna0(-1,-1,k,is),'b2stbr_sna'//chk//chns)
!        enddo
         call my_out_us(70,nCv,0,st%srw%b2stbr_sna(1,is),
     *                  'b2stbr_sna'//chns)
        endif
       enddo                                                          !srv 29.04.10
#ifdef B25_EIRENE
       do ia = 1,nnatmi
        write (chns,'(i3.3)') ia
        wrk0(:)=rfluxa(:,ia,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_rfluxa_eir'//chns)
        wrk0(:)=tfluxa(:,ia,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_tfluxa_eir'//chns)
        call my_out_us(70,nFc,0,fatm_eir(:,1,ia),              ! IYS 27.01.2017
     &                'b2stbr_fatmy_eir'//chns)
        wrk0(:)=pfluxa(:,ia,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_pfluxa_eir'//chns)
        call my_out_us(70,nFc,0,fatm_eir(:,0,ia),              ! IYS 27.01.2017
     &                'b2stbr_fatmx_eir'//chns)
        wrk0(:)=dab2(:,ia,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_dab2_eir'//chns)
        wrk0(:)=tab2(:,ia,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_tab2_eir'//chns)
       enddo
       do im = 1,nnmoli
        write (chns,'(i3.3)') im
        wrk0(:)=rfluxm(:,im,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_rfluxm_eir'//chns)
        wrk0(:)=tfluxm(:,im,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_tfluxm_eir'//chns)
        call my_out_us(70,nFc,0,fmol_eir(:,1,im),              ! IYS 27.01.2017
     &              'b2stbr_fmoly_eir'//chns)
        wrk0(:)=pfluxm(:,im,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_pfluxm_eir'//chns)
        call my_out_us(70,nFc,0,fmol_eir(:,0,im),              ! IYS 27.01.2017
     &              'b2stbr_fmolx_eir'//chns)
        wrk0(:)=dmb2(:,im,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_dmb2_eir'//chns)
        wrk0(:)=tmb2(:,im,1)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_tmb2_eir'//chns)
       enddo
#endif
      endif
      wrk0=she0_eir(:,0)+she0_eir(:,1)*st%pl%te+
     &       she0_eir(:,2)*st%dv%ne+she0_eir(:,3)*st%pl%te*st%dv%ne
      she_eir_tot = -sum(wrk0)
      if (switch%b2stbr_iout.eq.1 .or. switch%b2npht_iout.eq.1
     &                            .or. switch%iout_b2wdat.eq.4 ) then
       if (switch%use_eirene.ne.0) then
        if (switch%iout_b2wdat.ne.4) then
         do k=0,1  ! 3
          write (chk,'(i1)') k
          call my_out_us(70,nCv,0,she0_eir(1,k),'b2stbr_she_eir'//chk)
          call my_out_us(70,nCv,0,shi0_eir(1,k),'b2stbr_shi_eir'//chk)
          call my_out_us(70,nCv,0,shn0_eir(1,k),'b2stbr_shn_eir'//chk)
         enddo
        endif
        wrk0=she0_eir(:,0)+she0_eir(:,1)*st%pl%te+
     &       she0_eir(:,2)*st%dv%ne+she0_eir(:,3)*st%pl%te*st%dv%ne
        call my_out_us(70,nCv,0,wrk0,'b2stbr_she_eir')
        wrk0=shi0_eir(:,0)+shi0_eir(:,1)*st%pl%ti+
     &       shi0_eir(:,2)*st%dv%ni(:,0)+
     &       shi0_eir(:,3)*st%pl%ti*st%dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_shi_eir')
        wrk0=shn0_eir(:,0)+shn0_eir(:,1)*st%pl%tn+
     &       shn0_eir(:,2)*(st%dv%ni(:,0)-st%dv%ni(:,1))+
     &       shn0_eir(:,3)*st%pl%tn*(st%dv%ni(:,0)-st%dv%ni(:,1))
        call my_out_us(70,nCv,0,wrk0,'b2stbr_shn_eir')

       endif
!      do k=0,3
!       write (chk,'(i1)') k
!       call my_out_us(70,nCv,0,she0(-1,-1,k),'b2stbr_she'//chk)
!       call my_out_us(70,nCv,0,shi0(-1,-1,k),'b2stbr_shi'//chk)
!      enddo
       call my_out_us(70,nCv,0,st%srw%b2stbr_she,'b2stbr_she')
       call my_out_us(70,nCv,0,st%srw%b2stbr_shi,'b2stbr_shi')        !srv 29.04.10
       call my_out_us(70,nCv,0,st%srw%b2stbr_shn,'b2stbr_shn')        !srv 29.04.10
       call my_out_us(70,nCv,0,st%srw%b2stbr_skt,'b2stbr_skt')        !srv 29.04.10
       call my_out_us(70,nCv,0,st%srw%b2stbr_szt,'b2stbr_szt')        !srv 29.04.10
      endif
      if (switch%b2stbr_iout.eq.1 .or.
     &  (switch%b2npp7_iout.eq.1 .or. switch%b2nppo_iout.eq.1)
     &                           .or. switch%iout_b2wdat.eq.4 ) then
       if (switch%use_eirene.eq.0) then
!        do k=0,3
!          write (chk,'(i1)') k
!          call my_out_us(70,nCv,0,sch0(-1,-1,k),'b2stbr_sch'//chk)
!       enddo
        call my_out_us(70,nCv,0,st%srw%b2stbr_sch,'b2stbr_sch')
       endif
      endif                                                           !srv 11.09.09 }

* ..return
      if(ncall_b2stbr.eq.0 .or. main_call) ncall_b2stbr = ncall_b2stbr+1
      call subend ()
      return
      end subroutine b2stbr

*-----------------------------------------------------------------------
*.specification

*.end b2stbr

      subroutine b2stbr_neutr_scl(nCv, nFc, ns, mpg, geo, pl, dv,
     & use_eirene, ank_mods)
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef B25_EIRENE
      use b2mod_eirene
     , , only : neutr_scl
#endif
      use b2mod_subsys
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     , , only : ncall_b2stbr_neutr_scl
      implicit none
      integer, intent(in) :: nCv, nFc, ns, use_eirene, ank_mods
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv

      call subini ('b2stbr_neutr_scl')
      if (ncall_b2stbr_neutr_scl.eq.0) then
        call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
        call xertst (1.le.ns, 'faulty argument ns')
      endif
      if(use_eirene.ne.0.and.ank_mods.ne.0) then
#ifdef B25_EIRENE
c        write(*,*) 'Calling neutr_scl'
        call neutr_scl(nCv, nFc, ns, mpg, geo, pl, dv)
c        write(*,*) 'Called neutr_scl'
#else
        call xerrab ('b2stbr not compiled with the B25_EIRENE option')
#endif
      endif
      ncall_b2stbr_neutr_scl = ncall_b2stbr_neutr_scl+1
      call subend ()
      return

      end subroutine b2stbr_neutr_scl

*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr_init (ns, switch, mpg)
      use b2mod_types
      use b2mod_time
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_sputter
      use b2mod_neutrals_namelist
      use b2mod_neutr_src_scaling
      use b2mod_recycle
      use b2mod_switches
      use b2us_map
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: ns
      type (switches), intent(in) :: switch
      type (mapping), intent(inout) :: mpg
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
#ifdef WG_TODO
      integer nCg_tmp, ns_tmp
#endif
      integer is, istra
      integer neutrals_namelist
      integer sput_src, sput_chem_model
      logical new_sputter_namelist, catch_up, wall_namelist_used
#ifdef WG_TODO
      character*256 filename
#endif
*   ..procedures
      intrinsic mod
      external xertst, xerrab, ipgeti, ipgetr
*   ..initialisation
      data neutrals_namelist/0/,
     *  sput_src/1/, sput_chem_model/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbr_init')
      write(*,*) 'Start of b2stbr_init'
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read in needed data
      call ipgeti ('b2stbr_neutrals_namelist', neutrals_namelist)
      call ipgeti ('b2stbr_sput_src', sput_src)
      call ipgeti ('b2stbr_sput_dst', sput_dst)
      call ipgeti ('b2stbr_sput_dst2', sput_dst2)
      call ipgeti ('b2stbr_sput_dst3', sput_dst3)
      call ipgeti ('b2stbr_sput_chem_model', sput_chem_model)
      if(switch%b2stbr_sput_frc.ne.0.0_R8) then
        call xertst (0.le.sput_src.and.sput_src.lt.ns,
     1     'faulty b2stbr_sput_src')
        call xertst (sput_dst.lt.ns, 'faulty b2stbr_sput_dst')
        call xertst (sput_dst2.lt.ns, 'faulty b2stbr_sput_dst2')
        call xertst (sput_dst3.lt.ns, 'faulty b2stbr_sput_dst3')
        if(sput_dst.ge.0 .and. switch%use_eirene.eq.0 .and.
     .     neutrals_namelist.eq.0) then
          if (sput_chem_model.eq.1) then
            write(*,'(1x,a,i3,a,i3,a,g10.3)')
     .       'Chemical sputtering from ',sput_src,
     1       ' producing ',sput_dst,' at rate ',switch%b2stbr_sput_frc
          else if (sput_chem_model.eq.0) then
            write(*,'(1x,a,i3,a,i3,a,g10.3)')
     .       'Chemical sputtering from ',sput_src,
     1       ' producing ',sput_dst,
     2       ' with empirical formula multiplier ',
     3       switch%b2stbr_sput_frc
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

cxpb Reading of neutrals namelist moved here as it is needed by init_wall
      new_sputter_namelist=
     &         neutrals_namelist.eq.1.or.switch%use_eirene.ne.0

      call read_neutrals_namelist(ns, mpg, switch, new_sputter_namelist)
      call xertst (nstrat.le.nstraid,
     &   'Increase DEF_NSTRA in b2mod_dimensions')
      call alloc_b2mod_neutr_src_scaling(ns,nstrat)
      if(new_sputter_namelist) call write_b2mod_neutrals_namelist
      if(neutrals_time_mod.gt.0.0_R8) then
        neutrals_count=int(tim/neutrals_time_mod)
      else
        neutrals_count=0
      endif
      if(neutrals_time_mod.gt.0.0_R8) then
        catch_up=mod(tim,neutrals_time_mod).ge.neutrals_time_switch
     &    .and.neutrals_time_switch.gt.0.0_R8
      else if (neutrals_time_switch.gt.0.0_R8) then
        catch_up=tim.ge.neutrals_time_switch
      else
        catch_up=.false.
      end if
      do while(catch_up)
        write(*,*) 'b2stbr_init: catching up ...'
        call read_neutrals_namelist(ns, mpg, switch, .true.)
        new_sputter_namelist=.true.
        call xertst (nstrat.le.nstraid,
     &       'Increase DEF_NSTRA in b2mod_dimensions')
        call write_b2mod_neutrals_namelist
        if(neutrals_time_mod.gt.0.0_R8) then
          catch_up=mod(tim,neutrals_time_mod).ge.neutrals_time_switch
     &      .and.neutrals_time_switch.gt.0.0_R8
        else if (neutrals_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.neutrals_time_switch
        else
          catch_up=.false.
        endif
      enddo
      if (switch%b2stbr_neutrals_namelist.eq.1
     &                           .and. switch%use_eirene.eq.0) then
        do istra = 1, nstrai
          if (crcstra(istra).eq.'S' .or. crcstra(istra).eq.'N' .or.
     &        crcstra(istra).eq.'W' .or. crcstra(istra).eq.'E' .or.
     &        crcstra(istra).eq.'A') then
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
#ifdef WG_TODO
*   ..reading wall and sputter namelist data on first call
      write(*,*) 'Calling alloc_b2mod_sputter'
      call alloc_b2mod_sputter(nx,ny,ns)
      write(*,*) 'Calling alloc_b2mod_wall'
      call alloc_b2mod_layer(nx,ny,nncut,nlim,nsts)
      call alloc_b2mod_wall_init(nx,ny,nncut)
#endif
      call alloc_b2mod_wall(ns)
      call init_wall(mpg)
      call read_b2mod_wall_namelist(99,'    ',wall_namelist_used)
#ifdef WG_TODO
      call xertst (ndepth_nml.gt.1.or.plate_model.eq.0,
     .   'Inconsistent values of ndepth_nml and plate_model')
      call get_wall_data (nx,ny,ns,switch%BoRiS,wall_namelist_used)
#ifdef ALLOCATE_AND_NAMELIST
      nCg_tmp = mpg%nCg
      ns_tmp = ns
#else
      nCg_tmp = 2*(DEF_NXD+2)+(DEF_NCUT+2)*(DEF_NYD+2)
      ns_tmp = DEF_NSD
#endif
      filename='b2.sputter_save.parameters'
      call find_file(filename,sputter_namelist_used)
      if (sputter_namelist_used) then
        call read_b2mod_sputter_namelist
     .     (99,'    ',sputter_namelist_used)
      endif
      nCg_sp = nCg_tmp
      ns_sp = ns_tmp
#endif

cmb  Load pre-evaluated integrals for advanced recyling mocdel based on TRIM data
      if (switch%recycle_afn.eq.1) then
        call read_recycle_fnn(ns,mpg%nCg,switch%afn_bcs_use_coarse)            !mb 12.07.17
      endif
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbr_init

      end subroutine b2stbr_init
