*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr_phys (nCv, nFc, nVx, ns, nscx, nscxmax, iscx,
     &  dtim, switch, geo, mpg, pl, dv, co, rt, rtw, st_ext, srw,
     &  tchem, tchee, tphys, tphye,
     &  thevp, thvpe, trese, tresn, trfln, trfle, boundary_namelist,
     &  neutral_sources_rescale, core_sources_rescale,
     &  sput_src, sput_chem_alpha, sput_phys_alpha, alpha,
     &  sput_chem_model, therm_evap, sput_res,
     &  reflection_on, sputter_energy_on, main_call, output,
     &  sput_frac_flag, sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1  sput_mixed_alpha, sput_mixed_beta,
     2  new_sputter_namelist, shi0_ff, f_redep,
     3  temperature_at_guard_cell, potential_at_guard_cell,
     4  recycled_neutrals_contr)
      use b2mod_types
!WG_TODO      use b2mod_sputter
      use b2mod_tallies
!WG_TODO      use b2mod_elements
      use b2mod_constants
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
!WG_TODO      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_recycle
      use b2mod_ad
     & , only : ncall_b2stbr_phys
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : my_out_folder
#ifdef WG_DEBUG
     &        , my_out_sf                                               !srv 27.06.21
#endif
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2Rates), intent (in) :: rt
      type (B2RatesWork), intent (in) :: rtw
      type (B2StateExt), intent (in) :: st_ext
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1)
      integer sput_src, sput_chem_model, boundary_namelist,
     *     reflection_on, sputter_energy_on, output,
     *     sput_frac_flag, potential_at_guard_cell,
     *     temperature_at_guard_cell
      logical main_call, new_sputter_namelist
      real (kind=R8) :: dtim,
     &  f_redep(nCv,0:ns-1)
      real (kind=R8) ::
     *  neutral_sources_rescale, core_sources_rescale,
     &  sput_chem_alpha, sput_phys_alpha, alpha,
     &  therm_evap, sput_res,
     &  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     1  sput_mixed_alpha, sput_mixed_beta,
     2  recycled_neutrals_contr
*   ..output arguments (unspecified on entry)
      type (B2SourceWork), intent (inout) :: srw
      real (kind=R8) ::
     &  shi0_ff(nCv,0:nscx-1),
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
      integer i, iRc1, iRc2, iFc, is, istra
!     integer k
      character chns*3!, chk*1                                        !srv 11.07.05
      character*15 hlp_frm
#ifdef WG_TODO
      integer ib, itargsp, iwall, zbulk
      integer sput_dst_bulk
      integer, allocatable :: flux_x(:), flux_y(:)
      real (kind=R8), allocatable :: area(:)
      real (kind=R8) :: b2stbr_phys_sna(0:ns-1)
      real (kind=R8) :: lnorm, GammaEvap, redep_max,
     *   t0, t1, t2, t3, yldchm, yldres, yself_max,
     *   reflection_max, deposition_error(ntrack)
      logical l_chem_sput, match
#endif
      real (kind=R8) :: wrk0(nCv)
      real (kind=R8) ::
     *   qin, qel, surface_temp,
     *   Einc, Pflx, vbar, Phi_app
      logical lcore
c for KUL advanced fluid neutral boundary conditions
      real (kind=R8) :: dnndy, dnndz, dnn(nFc,0:1)
c     output for analysis/debugging
      real (kind=R8) :: fmomni_out(nCv),fmomrefl_out(nCv),
     &                  fmomrec_out(nCv), feneni_out(nCv),
     &                  fenerefl_out(nCv), fenerec_out(nCv),
     &                  fnni_out(nCv)
*   ..procedures
      intrinsic abs
      external xertst
      external b2xvsg, b2xvff, smin, smax
      real (kind=R8) ::
     *   smax, smin
*   ..initialisation
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Note: This version of b2stbr is for the case when recycling is
*     calculated via the neutrals namelist.
*
*     For given species IS on boundary ISTRA:
*     The particle recycling coefficient is given by B2RECYC(is,istra)
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
      if (ncall_b2stbr_phys.eq.0) then                              !srv 22.08.08 {
       write(*,*) 'First call of b2stbr_phys'
      endif                                                         !srv 22.08.08 }
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if(ncall_b2stbr_phys.lt.3) then
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
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv*2, dv%ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, pl%tn, 1, 'tn', '.gt.')
      endif
#ifdef WG_TODO
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
        do iCv=1,nCv
            b2stbr_phys_sna(is)=b2stbr_phys_sna(is)+
     1       srw%sna0(iCv,0,is)+srw%sna0(iCv,1,is)*pl%na(iCv,is)
        enddo
      enddo
      if((main_call.and.output.ge.1).or.output.gt.1)
     *     write(*,*) 'b2stbr_phys_sna(start): ',b2stbr_phys_sna

      if(ncall_b2stbr_phys.eq.0.and.sput_dst.ge.0.and.
     &   sput_dst.lt.ns) then
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
#endif

c pre-computations for KUL fluid neutrals
      dnn = 0.0_R8      
      if (switch%recycle_afn.ne.0) then

c     initialize debugging/analysis arrays
        fmomrec_out  = 0.0_R8
        fmomrefl_out = 0.0_R8
        fmomni_out   = 0.0_R8
        fenerec_out  = 0.0_R8
        fenerefl_out = 0.0_R8
        feneni_out   = 0.0_R8
        fnni_out     = 0.0_R8
        fna_mol      = 0.0_R8  !nh 30.03.21

        call precompute_kul_quant(nCv, nFc, nVx, switch, geo, mpg,
     &   pl, dv, co, rt, dnn)
      endif !! end pre-computations KUL fluid neutrals

      do istra = 1, nstrat
#ifdef WG_TODO
        sput_dst_bulk=targsp(istra,1)
        if (sput_dst_bulk.eq.-1) then
          l_chem_sput=.false.
          zbulk=0
        else
          zbulk=nint(zn(sput_dst_bulk))
          l_chem_sput=(sput_dst_bulk.eq.sput_dst.or.chemsp(istra)).and.
     &       zbulk.eq.6
        endif
        if(ncall_b2stbr_phys.eq.0) then
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
#endif
* ..recycling model at south boundary
        if (crcstra(istra).eq.'A') then

*   ..loop over boundary cells
          do i = 1, mpg%rcCvP(istra,2)
            iRc1 = mpg%rcCv(mpg%rcCvP(istra,1)+i-1,1) ! the guard cell
            iRc2 = mpg%rcCv(mpg%rcCvP(istra,1)+i-1,2) ! the domain cell
            iFc  = mpg%cvFc(mpg%cvFcP(iRc1,1))        ! the guard face

            Phi_app = getPhiApplied(boundary_namelist,mpg,iRc1) ! the potential

            lcore = mpg%cvOnClosedSurface(iRc1) ! xpb: no sputtering in core!

*   ..compute 2-D target temperature evolution
#ifdef WG_TODO
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
     &         area(1:rc_list_size(istra)),switch%BoRiS)
              deallocate(flux_x,flux_y,area)
            endif
#endif

c for KUL recycling/reflection model
            dnndy =   mpg%rcFcOr(mpg%rcCvP(istra,1)+i-1)*dnn(iFc,1)
            dnndz =  -mpg%rcFcOr(mpg%rcCvP(istra,1)+i-1)*dnn(iFc,0)

            call recycle(istra,
     &           b2species_start(istra), b2species_end(istra),
     &           iRc1, iRc2, iFc, mpg%rcFcOr(mpg%rcCvP(istra,1)+i-1),
     &           Phi_app, lcore, dnndy, dnndz)
          enddo      ! end loop over boundary cells

        endif       ! crcstra(istra)

#ifdef WG_TODO
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
#endif
      enddo      ! istra

#ifdef WG_TODO
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
#endif
      if (switch%b2stbr_phys_iout.eq.1.or.switch%b2npmo_iout.eq.1) then                   !srv 11.07.05 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!         do k=0,3
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,smo0(-1,-1,k,is),
!     &                       'b2stbr_phys_smo'//chk//chns)
!         enddo
        wrk0=srw%smo0(:,0,is)+
     &       srw%smo0(:,1,is)*pl%ua(:,is)+
     &       srw%smo0(:,2,is)*pl%na(:,is)*mp*am(is)+
     &       srw%smo0(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_phys_smo'//chns)
       enddo
      endif
      if (switch%b2stbr_phys_iout.eq.1 .or.
     &                        switch%b2npco_iout.eq.1) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
!         do k=0,1
!          write (chk,'(i1)') k
!          call my_out(70,nx,ny,sna0(-1,-1,k,is),
!     &                       'b2stbr_phys_sna'//chk//chns)
!         enddo
        wrk0=srw%sna0(:,0,is)+srw%sna0(:,1,is)*pl%na(:,is)
        call my_out_us(70,nCv,0,wrk0,'b2stbr_phys_sna'//chns)
       enddo
      endif
      if (switch%b2stbr_phys_iout.eq.1 .or.
     &                       switch%b2npht_iout.eq.1) then
!        do k=0,3
!         write (chk,'(i1)') k
!         call my_out(70,nx,ny,she0(-1,-1,k),'b2stbr_phys_she'//chk)
!         call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbr_phys_shi'//chk)
!         call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbr_phys_sch'//chk)
!        enddo
         wrk0=srw%she0(:,0)+srw%she0(:,1)*pl%te+
     &        srw%she0(:,2)*dv%ne + srw%she0(:,3)*pl%te*dv%ne
         call my_out_us(70,nCv,0,wrk0,'b2stbr_phys_she')
         wrk0=srw%shi0(:,0)+srw%shi0(:,1)*pl%ti+
     &        srw%shi0(:,2)*dv%ni(:,0)+srw%shi0(:,3)*pl%ti*dv%ni(:,0)
         call my_out_us(70,nCv,0,wrk0,'b2stbr_phys_shi')
      endif
      if (switch%b2stbr_phys_iout.eq.1 .or.
     &   (switch%b2npp7_iout.eq.1 .or. switch%b2nppo_iout.eq.1)) then
         wrk0=srw%sch0(:,0)+srw%sch0(:,1)*pl%po+
     &        srw%sch0(:,2)*dv%ne+srw%sch0(:,3)*pl%po*dv%ne
         call my_out_us(70,nCv,0,wrk0,'b2stbr_phys_sch')
      endif
      if(switch%b2stbr_phys_iout.eq.1) then
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

      if (switch%AFN_out.ne.0) then
        call my_out_us(70,nCv,0,fmomrec_out,'AFN_fmomrec')
        call my_out_us(70,nCv,0,fmomni_out,'AFN_fmomni')
        call my_out_us(70,nCv,0,fmomrefl_out,'AFN_fmomrefl')
        call my_out_us(70,nCv,0,fenerec_out,'AFN_fenerec')
        call my_out_us(70,nCv,0,feneni_out,'AFN_feneni')
        call my_out_us(70,nCv,0,fenerefl_out,'AFN_fenerefl')
        call my_out_us(70,nCv,0,fnni_out,'AFN_fnni')
      endif

#ifdef WG_TODO
      if(reflection_max.gt.1.0_R8 .or.
     .  ((main_call.and.output.ge.1).or.output.gt.1))
     1 write(*,*) 'Maximum reflection coefficient ', reflection_max
      if (redep_max.gt.0.0_R8 .and. (main_call.or.output.gt.1))
     . write(*,*) 'Maximum redeposition coefficient ', redep_max
      if(redep_max.gt.0.0_R8 .and. (yself_max.ge.1.0_R8 .or.
     .  ((main_call.and.output.ge.1).or.output.gt.1)))
     . write(*,*) 'Maximum self-sputtering coefficient ', yself_max
#endif
* ..return
      ncall_b2stbr_phys = ncall_b2stbr_phys + 1
      call subend ()
      return

*-----------------------------------------------------------------------

      contains

      subroutine recycle(istra, nsmin,nsmax,iCv,iCn,iFc,isign,Phi_app,
     &                   lcore,dnndy, dnndz)
c iCv             number of the guard cell
c iCn             number of the neighbour cell
c iFc             number of the guard cell face
c isign           multiplies fluxes to get the right sign
c Phi_app         applied potential
c area            area element
c lcore           is this a core element?
c dnndy           the gradient of the main neutral species density along the length of the boundary face
c dnndz           the gradient of the main neutral species density across the boundary face
      implicit none
      integer, intent(in) :: istra, nsmin, nsmax, iCv, iCn, iFc
      real (kind=R8) Phi_app, pof, isign, dnndy, dnndz
      logical lcore, found

      integer is0, is1, is2, is, isi, k, iwall
      real (kind=R8) :: tef, tif, tnf, taf, zaf(0:ns-1),
     &     rpf(0:ns-1), ptf, t0
#ifdef WG_TODO
     &   , sput_frac(0:ntrack), sput_frac_chem(0:ntrack), yself, eself,
     &     sput_frac_chem_sum, be_frac, be_suppression_factor,
     &     ytmp, ztmp, wrk(ndepth_nml)
      integer js, itrack
      real (kind=R8) mass_density, local_reflection
      character*13 hlp_frm
      external mass_density
#endif
c new variables for advanced fluid neutral reflection/recycling model
      real (kind=R8) :: fnni,fmomni,nni,nnwwni,feneni,fnnrefl,fnnrec,
     &     fmomrefl,fmomrec,fenerec,fenerefl,fene_el, nnrefl,nnwwnrefl,
     &     nnrec,nnwwnrec,nnsum,nnwwnsum,snan,smon,farea,q_conv,
     &     sina, cosa

#ifdef WG_TODO
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
        layer_nrelconstituents(iwall,1:matsurf_nconstituents(iwall))=
     1    matsurf_relconstituents
     .    (iwall,1:matsurf_nconstituents(iwall))*
     2   sput_frac(0)
cxpb  Contributions from the deposited species
        do itrack=1,ntargsp
           layer_nzconstituents
     .      (iwall,matsurf_nconstituents(iwall)+itrack)=
     1      nint(zn(target_intcode(itrack)))
           layer_nrelconstituents
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
#endif

c if nsmin = b2species_start = -1, this means that the spatially hybrid method for neutrals is used
c and this strata should be ignored here
      if (nsmin.lt.0) then
        return
      endif

      if (temperature_at_guard_cell.eq.1) then
        tef = pl%te(iCv)
        tif = pl%ti(iCv)
        tnf = pl%tn(iCv)
      else
        tef = (pl%te(mpg%fcCv(iFc,1))*geo%fcHc(iFc,2)+
     &         pl%te(mpg%fcCv(iFc,2))*geo%fcHc(iFc,1))/
     &        (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
        tif = (pl%ti(mpg%fcCv(iFc,1))*geo%fcHc(iFc,2)+
     &         pl%ti(mpg%fcCv(iFc,2))*geo%fcHc(iFc,1))/
     &        (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
        tnf = (pl%tn(mpg%fcCv(iFc,1))*geo%fcHc(iFc,2)+
     &         pl%tn(mpg%fcCv(iFc,2))*geo%fcHc(iFc,1))/
     &        (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
      endif
      if (potential_at_guard_cell.eq.1) then
!        write(6,*) 'Potential taken at guard cell (phys)'
        pof = pl%po(iCv)
      else
!        write(6,*) 'Potential taken at real cell (phys)'
        pof = pl%po(iCn)
      endif
      zaf(0:ns-1) =
     & (rt%rza(mpg%fcCv(iFc,1),0:ns-1)*geo%fcHc(iFc,2)+
     &  rt%rza(mpg%fcCv(iFc,2),0:ns-1)*geo%fcHc(iFc,1))/
     & (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
      rpf(0:ns-1) =
     & (rt%rpt(mpg%fcCv(iFc,1),0:ns-1)*geo%fcHc(iFc,2)+
     &  rt%rpt(mpg%fcCv(iFc,2),0:ns-1)*geo%fcHc(iFc,1))/
     & (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
      qin = max(0.0_R8,isign*(dv%fhi(iFc,0) + dv%fhi(iFc,1) +
     &  st_ext%fhi(iFc,0) + st_ext%fhi(iFc,1)))/geo%fcS(iFc)
      do is=0,ns-1
        qin = qin+max(0.0_R8,isign*(dv%fna(iFc,0,is)+dv%fna(iFc,1,is)))*
     &    ((dv%kinrgy(iCv,is)+tif)*(1.0_R8-switch%BoRiS) + rpf(is)*ev)/
     &     geo%fcS(iFc)
      enddo
      do is=0,st_ext%ns-1
        ptf =
     &   (st_ext%pt(mpg%fcCv(iFc,1),is)*geo%fcHc(iFc,2)+
     &    st_ext%pt(mpg%fcCv(iFc,2),is)*geo%fcHc(iFc,1))/
     &   (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
        taf =
     &   (st_ext%ta(mpg%fcCv(iFc,1),is)*geo%fcHc(iFc,2)+
     &    st_ext%ta(mpg%fcCv(iFc,2),is)*geo%fcHc(iFc,1))/
     &   (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))
        qin = qin + max(0.0_R8,isign*(st_ext%fa(iFc,0,is) +
     &                            st_ext%fa(iFc,1,is)))
     &        *((0.5_R8*st_ext%am(is)*mp*st_ext%ua(iCv,is)**2+taf)*
     &          (1.0_R8-switch%BoRiS)
     &          + ptf*ev)/geo%fcS(iFc)
      enddo
      qel = max(0.0_R8,isign*(dv%fhe(iFc,0)+dv%fhe(iFc,1)+
     &  (dv%fne(iFc,0)+dv%fne(iFc,1))*tef*(1.0_R8-switch%BoRiS)))/
     &   geo%fcS(iFc)
      if(lcore) then
         surface_temp = 0.0_R8
      else
         surface_temp = 0.0_R8
#ifdef WG_TODO
         if(iwall.eq.0) then
            write(*,*) 'xymap returned zero for ix,iy = ',ix,iy
            call xerrab ('Invalid iwall setting!')
         endif
         if (plate_model.eq.0) then
            call heatdiff0D(iwall,qin,qel,surface_temp)
            target_temp(iwall,1) = surface_temp
         else if (plate_model.eq.1) then
            if (main_call) then
              wrk(1:ndepth_nml) = target_temp(iwall,1:ndepth_nml)
              call heatdiff1D(iwall,ns,qin+qel,
     .           dtim*plate_time_factor(iwall),wrk,output)
              target_temp(iwall,1:ndepth_nml) = wrk(1:ndepth_nml)
            end if
            surface_temp = target_temp(iwall,1)
         else if (plate_model.eq.2) then
           surface_temp = target_temp(iwall,1)
         endif
#endif
      endif

*    ..loop over species
cmb   Loop implicitly expects that the species are arranged per isonuclear sequence and that the neutral atom of that sequence always comes first. If not, erroneous results will be obtained.
      is0 = 0
      is1 = 0
      nnsum = 0.0_R8
      nnwwnsum = 0.0_R8
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
            if(ncall_b2stbr_phys.eq.0)
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
         t0 = max(0.0_R8,isign*(dv%fna(iFc,0,is) + dv%fna(iFc,1,is))) *
     1        neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         Einc=max(0.0_R8,tif/ev+(pof-Phi_app)*zaf(is))
                                 ! xpb: Incident energy (eV)
         Pflx = t0/geo%fcS(iFc)  ! xpb: Incident particle flux (part/m^2/s)
         if (is_neutral(is)) then
            vbar = sqrt(8.0_R8 * pl%ti(iCv) / (pi * am(is) * mp))
            Pflx = Pflx + pl%na(iCv,is)*vbar*alpha *
     1           neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         endif
         if (switch%recycle_afn.eq.0.or.zn(is0).ne.1) then !! default neutral recycling/reflection model (always for non-hydrogenic species)
           srw%sna0(iCv,0,is0) = srw%sna0(iCv,0,is0) +
     &      b2recyc(is,istra)*t0
#ifdef WG_TODO
         itrack=track_index(is0)
cxpb Whatever is not recycled is counted as deposition
         if(itrack.gt.0.and.itrack.le.ntrack
     1        .and.iwall.gt.0.and.main_call) then
            deposition(iwall,itrack)=deposition(iwall,itrack)+
     1           dtim*plate_time_factor(iwall)*
     2           (1.0_R8-b2recyc(is,istra))*t0
         endif
#endif
cwdk  recycle toroidal componenent of parallel momentum, factor mrecyc
c     (toroidal piece only, because poloidal piece of parallel velocity flips sign after reflection)
         srw%smo0(iCv,0,is0) = srw%smo0(iCv,0,is0) + mrecyc(is,istra)*
     &      am(is)*mp*pl%ua(iCv,is)*geo%cvBb(iCv,2)**2/
     &      geo%cvBb(iCv,3)**2*b2recyc(is,istra)*t0*geo%cvHz(iCv)
         if (switch%lin_shi0.eq.1) then
           srw%shi0(iCv,1) = srw%shi0(iCv,1)+recycled_neutrals_contr*
     &        max(erecyc(is,istra)*(1.5_R8+switch%BoRiS)*Einc*ev,
     &        (1.5_R8+switch%BoRiS)*surface_temp*KBolt)*
     &        b2recyc(is,istra)*t0/pl%ti(iCv)
         else
           srw%shi0(iCv,0) = srw%shi0(iCv,0)+recycled_neutrals_contr*
     &        max(erecyc(is,istra)*(1.5_R8+switch%BoRiS)*Einc*ev,
     &        (1.5_R8+switch%BoRiS)*surface_temp*KBolt)*
     &        b2recyc(is,istra)*t0
         endif
         if (is2.ne.-1) shi0_ff(iCv,is2) = shi0_ff(iCv,is2)+
     &        recycled_neutrals_contr*max(
     &        erecyc(is,istra)*(1.5_R8+switch%BoRiS)*Einc*ev,
     &        (1.5_R8+switch%BoRiS)*surface_temp*KBolt)*
     &        b2recyc(is,istra)*t0
         ! xpb: recycled neutrals come back with at least thermal energy

         else !! use advanced fluid neutral boundary conditions by N. Horsten et. al
c lvn1
c**************************************************************************************
c*****************     New KU Leuven recycling model for neutrals *********************
c**************************************************************************************

c main reference:
c   Horsten, Niels, Giovanni Samaey, and M. Baelmans.
c  "Development and assessment of 2D fluid neutral models that include atomic databases and a microscopic reflection model."
c   Nuclear Fusion 57.11 (2017): 116043.

           farea = geo%fcs(iFc) !! area of the face
           cosa = geo%fcQalf(iFc,0)
           sina = geo%fcQalf(iFc,1)


!! end of preliminary calculations
           if(erecyc(is,istra)>0.0_R8.or.
     2      mrecyc(is,istra)>0.0_R8)
     3     write(*,*)'recycle_afn incompatible with user-defined '//
     4     'momentum and energy recycling coefficients > 0'

           if (is.ne.is0) then !! ion species, recycling


cmb   net momentum sources from neutral particle recombination
c            Calculate recycled neutral particle flux (assuming thin boundary cells)
             if (.not.lcore) then

               call CalcRecycledFluxes(iCv,mpg%nCi,iFc,isign,
     &              is0,is2,is,istra,iwall,
     &              Phi_app,farea,b2recyc(is,istra),recycm(is,istra),
     &              dv%fluid_frac_hyb(iFc),pl,tif,tef,pof,
     &              pl%na(iCv,is),dv%ne(iCv),geo,mpg,t0,
     &              cosa, sina, switch%use_uy_uz_0,
     &              fnnrec,fmomrec,nnrec,nnwwnrec,fenerec)
             else
                fnnrec = 0.0_R8
                fmomrec = 0.0_R8
                nnrec = 0.0_R8
                nnwwnrec = 0.0_R8
                fenerec = 0.0_R8
             endif

             snan = fnnrec*neutral_sources_rescale
             srw%sna0(iCv,0,is0) = srw%sna0(iCv,0,is0) + snan
             smon = fmomrec*geo%cvhz(iCv)*neutral_sources_rescale
             srw%smo0(iCv,0,is0) = srw%smo0(iCv,0,is0) + smon

             nnsum = nnsum + nnrec
             nnwwnsum = nnwwnsum + nnwwnrec
cmb   net neutral energy source from recycling

             if (switch%tn_style.eq.0) then
               srw%shi0(iCv,0) = srw%shi0(iCv,0) + fenerec*
     &                              neutral_sources_rescale
             elseif (switch%tn_style.eq.1) then
               !! neutral energy boundary conditions are not used
             else
               srw%shn0(iCv,0) = srw%shn0(iCv,0) + fenerec*
     &                              neutral_sources_rescale
             endif

             if (switch%AFN_out.ne.0) then
               fmomrec_out(iCv)   = fmomrec
               fenerec_out(iCv)   = fenerec
             endif

           else !! this else belongs to "if(.not.(is.eq.is0))", hence 'is' will be neutral here

c          Identify corresponding ion species
             if (lnext(is,is+1)) then
               isi = is+1
             else
               call xerrab('Problem with species labeling')
             endif

c          A conversion is necessary because the BC implemented below were developed in a code that uses a total energy equation,
c          while B2.5 uses internal energy equations.


             q_conv =
     &             (isign*
     &             (tnf*pl%ua(iCv,is0)*geo%fcpbs(iFc)*
     &             pl%na(iCv,is0) +
     &             (dv%fhm(iFc,0,is0)+dv%fhm(iFc,1,is0)) ))*
     &             neutral_sources_rescale

             if (switch%tn_style.eq.0) then
               srw%shi0(iCv,0) = srw%shi0(iCv,0) + q_conv
             elseif (switch%tn_style.eq.1) then
               !! do nothing, no neutral contributions in the boundary conditions of the pure ion energy equation
             else
               srw%shn0(iCv,0) = srw%shn0(iCv,0) + q_conv
             endif

             if (maxw(istra).eq.0) then
               call CalcIncidentFluxesDiffusion(iCv,iCn,iFc,isign,
     &              is0,is2,isi,istra,
     &              Phi_app,farea,pl,tif,tnf,tef,pof,
     &              dv%ne(iCv),geo,rt,dnndy,dnndz,cosa,sina,
     &              switch%use_uy_uz_0,switch%use_uy_uz_1,
     &              switch%diffusion_bc_safeguard,
     &              fnni,fmomni,nni,nnwwni,feneni)
             elseif (maxw(istra).eq.1) then
               call CalcIncidentFluxesMaxwellian(iCv,iFc,is0,isign,
     &              istra,farea,pl,tnf,geo,cosa,sina,
     &              fnni,fmomni,nni,nnwwni,feneni)
             elseif (maxw(istra).eq.2) then
               call CalcIncidentFluxesKn(switch%Kn_b1,
     &              switch%Kn_b2,
     &              iCv,iCn,iFc,isign,
     &              is0,is2,isi,istra,
     &              Phi_app,farea,pl,tif,tnf,tef,pof,
     &              dv%ne(iCv),geo,rt,dnndy,dnndz,cosa,sina,
     &              switch%use_uy_uz_0,switch%use_uy_uz_1,
     &              switch%L_macro_AFN,
     &              fnni,fmomni,nni,nnwwni,feneni)
             else
               call CalcIncidentFluxesMinKn(iCv,iCn,iFc,isign,
     &              is0,is2,isi,istra,
     &              Phi_app,farea,pl,tif,tnf,tef,pof,
     &              dv%ne(iCv),geo,rt,dnndy,dnndz,cosa,sina,
     &              switch%use_uy_uz_0,switch%use_uy_uz_1,
     &              fnni,fmomni,nni,nnwwni,feneni)

             endif

             srw%sna0(iCv,0,is0) = srw%sna0(iCv,0,is0) -
     &         fnni*neutral_sources_rescale

!! save the incident fluid neutral flux
             dv%fnn_inc(iFc,is0) = fnni*isign

cmb   net momentum sources from neutral particle reflection
             if (.not.lcore) then
               if (maxw(istra).eq.0) then
                 call CalcReflectedFluxesDiffusion(iCv,iCn,mpg%nCi,
     &              iFc,isign,is0,is2,isi,istra,iwall,
     &              Phi_app,farea,b2recyc(is,istra),recycm(is,istra),
     &              dv%fluid_frac_hyb(iFc),
     &              pl,tif,tnf,
     &              tef,pof,dv%ne(iCv),
     &              geo,mpg,rt,dnndy,dnndz,
     &              cosa, sina, switch%use_uy_uz_0,switch%use_uy_uz_1,
     &              switch%diffusion_bc_safeguard,
     &              fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &              fene_el)
               elseif (maxw(istra).eq.1) then
                 call CalcReflectedFluxesMaxwellian(iCv,mpg%nCi,iFc,
     &              isign,is0,
     &              istra,iwall,farea,b2recyc(is,istra),
     &              recycm(is,istra),dv%fluid_frac_hyb(iFc),pl,tnf,
     &              geo,mpg,cosa, sina,
     &              fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &              fene_el)
               elseif (maxw(istra).eq.2) then
                 call CalcReflectedFluxesKn(switch%Kn_b1,
     &              switch%Kn_b2,
     &              iCv,iCn,mpg%nCi,
     &              iFc,isign,is0,is2,isi,istra,iwall,
     &              Phi_app,farea,b2recyc(is,istra),recycm(is,istra),
     &              dv%fluid_frac_hyb(iFc),pl,tif,tnf,
     &              tef,pof,dv%ne(iCv),
     &              geo,mpg,rt,dnndy,dnndz,
     &              cosa, sina, switch%use_uy_uz_0, switch%use_uy_uz_1,
     &              switch%L_macro_AFN, fnnrefl,fmomrefl,nnrefl,
     &              nnwwnrefl,fenerefl,fene_el)
               else
                 call CalcReflectedFluxesMinKn(iCv,iCn,mpg%nCi,
     &              iFc,isign,is0,is2,isi,istra,iwall,
     &              Phi_app,farea,b2recyc(is,istra),recycm(is,istra),
     &              dv%fluid_frac_hyb(iFc),pl,tif,tnf,
     &              tef,pof,dv%ne(iCv),
     &              geo,mpg,rt,dnndy,dnndz,
     &              cosa, sina, switch%use_uy_uz_0, switch%use_uy_uz_1,
     &              fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
               endif
             else
               fnnrefl = 0.0_R8
               fmomrefl = 0.0_R8
               nnrefl = 0.0_R8
               nnwwnrefl = 0.0_R8
               fenerefl = 0.0_R8
               fene_el = 0.0_R8
             endif

cnh  adding the reflected particle source (with the molecular part subtracted)
             srw%sna0(iCv,0,is0) = srw%sna0(iCv,0,is0) +
     1         fnnrefl*neutral_sources_rescale
             snan = fnnrefl


             if (switch%AFN_out.ne.0) then
               fmomni_out(iCv)   = fmomni
               fmomrefl_out(iCv) = fmomrefl
               feneni_out(iCv)   = feneni
               fenerefl_out(iCv) = fenerefl
               fnni_out(iCv)     = fnni
             endif

             smon = (fmomni+fmomrefl)*geo%cvhz(iCv)*
     &              neutral_sources_rescale
             srw%smo0(iCv,0,is0) = srw%smo0(iCv,0,is0) + smon

             nnsum = nnsum + nni + nnrefl
             nnwwnsum = nnwwnsum + nnwwni + nnwwnrefl

             if (switch%tn_style.eq.0) then
               srw%shi0(iCv,0) = srw%shi0(iCv,0) + (feneni + fenerefl)*
     &                                          neutral_sources_rescale
             elseif (switch%tn_style.eq.1) then
               !! neutral energy boundary conditions are not used
             else
               srw%shn0(iCv,0) = srw%shn0(iCv,0) + (feneni + fenerefl)*
     &                                          neutral_sources_rescale
             endif

c remove the Franck-Condon energy, that has been provided to the neutrals, from the electrons at the boundary
c temporary switch to evaluate the effect.
             if (switch%remove_FC_el.eq.1) then
               srw%she0(iCv,0) = srw%she0(iCv,0) + (fene_el)*
     &                                          neutral_sources_rescale
             endif

           endif
           if (zamax(is).eq.zn(is)) then   !! see eq. 62-63 of the main reference.
               if (nnsum.eq.0.0_R8) nnsum = 1.0_R8
               smon = -isign*
     &            cosa*am(is0)*mp/3.0_R8*nnwwnsum**2/
     &            nnsum*farea*abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))*
     &            geo%cvhz(iCv)
               srw%smo0(iCv,0,is0) = srw%smo0(iCv,0,is0) + smon*
     &                                    neutral_sources_rescale
             nnsum = 0.0_R8
             nnwwnsum = 0.0_R8
           endif

       endif !! end KU Leuven advanced fluid neutral models
#ifdef WG_TODO
         if(.not.lcore) then    ! xpb: no sputtering in core!
c chemical sputtering
            if(chem_sput(is,istra).gt.0.0_R8.and.
     &           nint(zn(is)).eq.sput_src.and.
     &           (l_chem_sput.or.sput_frac_flag.ne.0)) then
               if(zbulk.eq.6.or.sput_frac_flag.ne.0) then      ! chemical sputtering only for C
                  if(sput_chem_model.eq.1) then
                     if (Einc.gt.0.0_R8) then
                        ytmp = chem_sput(is,istra) *
     &                       exp(-(sput_chem_cutoff_alpha/Einc)**
     &                       sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
                        t1=t0 * ytmp
                        chemical_sputtering(iwall,is) = ytmp
                     else
                        t1=0.0_R8
                        chemical_sputtering(iwall,is) = 0.0_R8
                     endif
                  else if (zbulk.eq.6) then
                     call chmsput(iwall,Einc,zaf(is),zn(is),
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
                     call chmsput(iwall,Einc,zaf(is),zn(is),
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
               if(ncall_b2stbr_phys.ne.0 .or.
     &            sputter_namelist_used) then
                  t1=(sputter_yield(ix,iy,is,1)+sput_chem_alpha*t1)/
     *                 (1+sput_chem_alpha)
               endif
               t2=t1*(1.5_R8+switch%BoRiS)*surface_temp*KBolt
               sputter_yield(ix,iy,is,1)=t1
               tchem(is)=tchem(is)+t1
               if (sput_dst_bulk.ne.-1) then
                 sna0(ix,iy,0,sput_dst_bulk) =
     &              sna0(ix,iy,0,sput_dst_bulk) + t1
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
                  if(ncall_b2stbr_phys .ne. 0 .or.
     &               sputter_namelist_used) then
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
                    ytmp=sputter(is,zaf(is),sput_dst_bulk,Einc*ev) *
     1               phys_sput(is,istra) *
     2               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))/
     3               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself)
! the two lines above are the reduction of sputtering due to prompt redeposition
! and then the multiplication due to self-sputtering of the redeposited particles
                  else
                    ytmp=sputter(is,zaf(is),sput_dst_bulk,Einc*ev) *
     1               phys_sput(is,istra) *
     2               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))/
     3               (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk)*yself*
     4                sput_frac(0))
! we only multiply the self-sputtering cascade contribution by sput_frac(0) since the
! primary sputtered particles will be multiplied by this factor just below
                  endif
                  t1=t0 * ytmp
                  if(ncall_b2stbr_phys.ne.0 .or.
     &               sputter_namelist_used) then
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
     .               sputter(is,zaf(is),target_intcode(itrack),Einc*ev)*
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
                     ytmp=sputter(is,zaf(is),sput_dst_bulk,Einc*ev) *
     1                 phys_sput(is,istra)
                     ztmp=sputter2(is,sput_dst_bulk,Einc*ev) *
     1                 phys_sput(is,istra) *
     2                 (1.0_R8-f_redep(ixn,iyn,sput_dst_bulk))
                     if(sput_frac_flag.ne.0) then
                       ytmp=ytmp*sput_frac(0)
                       ztmp=ztmp*sput_frac(0)
                       yself=yself*sput_frac(0)
                     endif
                     t2=t0 * ztmp * Einc * ev
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
                     if(ncall_b2stbr_phys .ne. 0 .or.
     &                  sputter_namelist_used) then
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
     1                               target_intcode(itrack),Einc*ev) *
     1                       phys_sput(is,istra) * sput_frac(itrack)
                        ztmp=sputter2(is,target_intcode(itrack),Einc*ev)
     1                     * phys_sput(is,istra) * sput_frac(itrack) *
     2                       (1.0_R8-
     3                        f_redep(ixn,iyn,target_intcode(itrack)))
                        t2=t0 * ztmp * Einc * ev
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
               t2=t1*(1.5_R8+switch%BoRiS)*surface_temp*KBolt
               sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
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
                  call ressput(Einc,zaf(is),zn(is),zn(sput_dst_bulk),
     &              am(is),am(sput_dst_bulk),surface_temp,
     &              Pflx/neutral_sources_rescale,sput_dst_bulk,yldres)
                  RES_sputtering(iwall,is) = yldres
                  t1=yldres*neutral_sources_rescale*
     &              area*sput_res
               endif
               if(sput_frac_flag.ne.0) then
                  if(zbulk.ne.6) then
                    call ressput(Einc,zaf(is),zn(is),6.0_R8,
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
               if(sput_dst_bulk.ne.-1) sna0(ix,iy,0,sput_dst_bulk) =
     &              sna0(ix,iy,0,sput_dst_bulk) + t1
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
               t2=t1*(1.5_R8+switch%BoRiS)*surface_temp*KBolt
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               tresn(is)=tresn(is)+t1
               trese(is)=trese(is)+t2
            else
               RES_sputtering(iwall,is) = 0.0_R8
            endif
c back scattering ! mixed materials included thru layer constituents definitions
            if (b2recyc(is,istra).eq.0.0_R8 .and.
     &         erecyc(is,istra).eq.0.0_R8 .and. .not.is_neutral(is)
     &   .and. reflection_on.ne.0) then ! use reflection model
               local_reflection = reflect(is,Einc*ev,iwall)
               reflection_max = max(reflection_max, local_reflection)
               t1=t0 * local_reflection
               t2=t0 * reflect2(is,Einc*ev,iwall)*Einc*ev
               itrack=track_index(is0)
cxpb  Removal from the deposited material if backscattered
cxpb  because incoming non-recycling flux was originally counted
cxpb  as deposition
               if(itrack.gt.0.and.itrack.le.ntrack
     1              .and.iwall.gt.0.and.main_call) then
                  deposition(iwall,itrack)=deposition(iwall,itrack)-
     1                 dtim*plate_time_factor(iwall)*t1
               endif
               backscattering(iwall,is) = reflect(is,Einc*ev,iwall)
               backscattering_energy(iwall,is) =
     &              reflect2(is,Einc*ev,iwall)
               if(reflection_on.eq.1) then
                 sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
               else
                 sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0) + t1    ! reflection to the neutral
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
#endif
*    ..recycling into next ionized species
         if(is.lt.ns-1) then
           if(lnext(is,is+1)) then
              srw%sna0(iCv,0,is+1) = srw%sna0(iCv,0,is+1)+
     2        rcion(is,istra)*t0/neutral_sources_rescale
              if (is_neutral(is+1)) then           ! IYS 18.09.2018
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
               if(ncall_b2stbr_phys.eq.0)
     1          write(*,'(a,2i3)') 'b2stbr: setting is0 to is1 ',is0,is1
               is0=is1
            endif
cxpb
*          (species (is) recycles into is0)
            t0 = max(0.0_R8,isign*(dv%fna(iFc,0,is) + dv%fna(iFc,1,is)))
     .          *core_sources_rescale ! xpb: 2002.02.01
            if (is.ne.is0) then
               srw%sna0(iCv,0,is0) = srw%sna0(iCv,0,is0)+
     &          b2recyc(is,istra)*t0

cwdk  recycling of true parallel momentum
c     (at core, cell faces will be aligned with flux surfaces)
               srw%smo0(iCv,0,is0) = srw%smo0(iCv,0,is0)+
     .           mrecyc(is,istra)*am(is)*mp*pl%ua(iCv,is)*
     .           b2recyc(is,istra)*t0*geo%cvHz(iCv)
               srw%shi0(iCv,0) = srw%shi0(iCv,0)+erecyc(is,istra)*
     .          (1.5_R8+switch%BoRiS)*(zaf(is)*tef+tif)*
     &           b2recyc(is,istra)*t0*recycled_neutrals_contr            !srv
            endif
         enddo      ! do is = nsmax, nsmin, -1
      endif

      return
      end subroutine recycle

*-----------------------------------------------------------------------
*.end b2stbr_phys

      end subroutine b2stbr_phys
