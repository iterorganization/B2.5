c**************************************************************************************
c********************  NEW KU Leuven module for artificial sinks  *********************
c**************************************************************************************
! mb  10/11/23
! main author: Anthony Piras
 
      module b2mod_artificialsinks
*   ..other modules called 
      use b2mod_types
      use b2mod_constants

      implicit none
*   ..declarations allocatable variable
*   ..The allocatable variables hereunder have the same dimensions of the ones defined in b2us_plasma
*   ..plus an extra dimension equal to 3. The logic is as follows:
*   ..sna_losses(:,:,:,1) = sna0_limloss(:,:,:)   for limiter losses
*   ..sna_losses(:,:,:,2) = sna0_driftloss(:,:,:) for drift losses
*   ..sna_losses(:,:,:,3) = sna0_polloss(:,:,:)   for poloidal losses
*   ..and so on and so forth for the other variables
      real (kind=R8), allocatable, save ::
     &  sna_losses(:,:,:,:), ! particles losses
     &  she_losses(:,:,:),   ! electron energy losses
     &  shi_losses(:,:,:)    ! ion energy losses

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ARTIFICIALSINKS is intended to simulate artificial sinks due to losses in the poloidal direction,
*     drifts and eventual losses to limiters. These sinks are used when having SOLPS-ITER running
*     1D models (i.e. symmetry BCs applied in poloidal direction, purely radial model).
*     Subroutines contained in this module should not be used when running a full 2D simulation.
*     Verify that losses_limiters, losses_drifts and losses_poloidal switches have their default
*     value (zero), thus being turned off.
*     This losses model comes from the TOMATOR1D code (credits: Tom Wauters tom.wauters@iter.org).
*     The losses are stored in vectors sna0_XXXloss, she0_XXXloss and shi0_XXXloss, all of them of
*     the same size of st%sr%sna/she/shi.
*     Particle losses are assumed to be of the form n_i/tau, hence only sna0_XXXloss(:,1,:) is
*     non-zero. The subroutines return sna0_XXXloss(:,1,:) = -1.0/tau, then the multiplication by
*     n_i is carried out in b2sral.
*     Energy losses are assumed to be of the form n*T*gamma/tau, hence only shY0_XXXloss(:,3) is non-zero.
*     In the case of limiter losses, energy losses for electrons are assumed to be of the form
*     n_a*T_e*gamma/tau, hence only she0_XXXloss(:,1) is non-zero.
*
*
*-----------------------------------------------------------------------

      contains

c Alloc
      subroutine alloc_b2mod_artificialsinks(nCv, ns)
      implicit none
      integer, intent(in) :: nCv, ns

      allocate(sna_losses(nCv,0:1,ns,3))
      allocate(she_losses(nCv,0:3,3))
      allocate(shi_losses(nCv,0:3,3))

      sna_losses(:,:,:,:) = 0.0_R8
      she_losses(:,:,:) = 0.0_R8
      shi_losses(:,:,:) = 0.0_R8

      return
      end subroutine alloc_b2mod_artificialsinks
c End alloc

c Limiter losses computation
      subroutine limiterloss(nCv, ns, is,
     &   switch, geo, mpg, pl, dv, rt, st_ext, 
     &   sna0_limloss, shi0_limloss, shn0_limloss, she0_limloss)
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_switches
      implicit none
      
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..input/output arguments
      real (kind=R8) :: sna0_limloss(nCv,0:1,0:ns-1),
     &  shi0_limloss(nCv,0:3), she0_limloss(nCv,0:3),
     &  shn0_limloss(nCv,0:3)
*   ..workspace arguments
      integer :: is, iCv
      real (kind=R8) :: c_sound(nCv), tauion_lim(nCv,0:ns-1),
     &  Rp_min, Rp_max
*   ..local variables
      real (kind=R8) :: rz(nCv), c_sound2
      integer :: cntr, tmpry
      character chns*3

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('limiterloss')

      call sfill (nCv*2*ns, 0.0_R8, sna0_limloss, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_lim, 1)
      call sfill (nCv*4, 0.0_R8, shi0_limloss, 1)
      call sfill (nCv*4, 0.0_R8, shn0_limloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_limloss, 1)

      Rp_min = switch%R0 - switch%limiter_width * switch%a0
      Rp_max = switch%R0 + switch%limiter_width * switch%a0

      do is = 0, ns-1
        if(is_neutral(is)) cycle
          do iCv = 1, mpg%nCi
            if(geo%cvX(iCv).le.Rp_min.or.geo%cvX(iCv).ge.Rp_max) then
              c_sound2 = sqrt(rt%rza(iCv,is)*pl%te(iCv)/(am(is)*mp)) !! CHECK B2NEWS_.F FOR AUXILIARY  PPZ AND PRZ
              call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
              call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti,
     &                  dv%pz, st_ext)
              c_sound = sqrt(dv%pz/rz)

              tauion_lim(iCv,is) = 0.66_R8*2.0_R8*pi*geo%cvX(iCv)/ ! 0.66 numerical factor by Tom (see emails)
     &                        switch%nlimiters/c_sound(iCv)
              if(tauion_lim(iCv,is).ne.0.0_R8) then
                sna0_limloss(iCv,1,is) = -1.0_R8/tauion_lim(iCv,is)*
     &                                   geo%cvVol(iCv)
                tmpry = int(floor(is - rt%rza(iCv,is))) !! LOOK INTO B2STBR_PHYS (RECYCLE) TO MODIFY
                sna0_limloss(iCv,0,tmpry) = pl%na(iCv,is)/
     &                             tauion_lim(iCv,is)*geo%cvVol(iCv)
c ap SINK AT T = TI
                shi0_limloss(iCv,3) = -1.0_R8/(tauion_lim(iCv,is)/
     &                             switch%bohm_gammaD)*geo%cvVol(iCv)
c ap SOURCE AT T = T_wall (=Ta0 in TOMATOR 1D)
                if (switch%tn_style.eq.2) then
                  shn0_limloss(iCv,0) = 3.0_R8/2.0_R8*pl%na(iCv,is)/
     &                          tauion_lim(iCv,is)*geo%cvVol(iCv)*
     &                          (switch%T_wall*ev) ! Ta0, has to be generalised
                else
                  shi0_limloss(iCv,0) = 3.0_R8/2.0_R8*pl%na(iCv,is)/
     &                          tauion_lim(iCv,is)*geo%cvVol(iCv)*
     &                          (switch%T_wall*ev) ! Ta0, has to be generalised
                endif
                she0_limloss(iCv,1) = she0_limloss(iCv,1) -
     &                                3.0_R8/2.0_R8*rt%rza(iCv,is)*
     &                                pl%na(iCv,is)/(tauion_lim(iCv,is)/
     &                                switch%bohm_gammaD)*geo%cvVol(iCv)
              endif
            endif
        enddo
      enddo

      if (switch%artloss_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        do is = 0,ns-1
          if(is_neutral(is)) cycle
            write (chns,'(i3.3)') is
            call my_out_us(70,nCv,0,sna0_limloss(:,0,is),
     &             'b2mod_artificialsinks_sna_limloss_0_species_'//chns)
            call my_out_us(70,nCv,0,sna0_limloss(:,0,is),
     &             'b2mod_artificialsinks_sna_limloss_1_species_'//chns)
        enddo
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,shi0_limloss(:,cntr),
     &                   'b2mod_artificialsinks_shi_limloss_'//chns)
          call my_out_us(70,nCv,0,she0_limloss(:,cntr),
     &                   'b2mod_artificialsinks_she_limloss_'//chns)
        enddo
      endif

      sna_losses(:,:,:,1) = sna0_limloss(:,:,:)
      she_losses(:,:,1) = she0_limloss(:,:)
      shi_losses(:,:,1) = shi0_limloss(:,:)

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end limiterloss

      end subroutine limiterloss

c Drift losses computation
      subroutine driftloss(nCv, ns, is,
     &   switch, geo, mpg, pl, dv, rt,
     &   sna0_driftloss, shi0_driftloss, she0_driftloss)
      use b2mod_b2cmpa
      use b2us_geo
      use b2mod_subsys
      use b2us_plasma
      use b2us_map
      use b2mod_switches

      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      real (kind=R8) :: sna0_driftloss(nCv,0:1,0:ns-1),
     &  shi0_driftloss(nCv,0:3), she0_driftloss(nCv,0:3)
*   ..workspace arguments
      integer :: is, iCv
      real (kind=R8) :: charge(nCv,0:ns-1), tauion_drift(nCv,0:ns-1)
      integer :: cntr
      character chns*3

* ..preliminaries
*   ..subprogram start-up calls
      call subini('driftloss')

      call sfill (nCv*2*ns, 0.0_R8, sna0_driftloss, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_drift, 1)
      call sfill (nCv*4, 0.0_R8, shi0_driftloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_driftloss, 1)

      do is = 0, ns-1
        if(is_neutral(is)) cycle
          charge(:,is) = rt%rza(:,is)*qe
          tauion_drift(:,is) = switch%a0/(2.0*pl%ti/charge(:,is)/
     &                         geo%cvBb(:,3)/geo%cvX(:))
          do iCv = 1, mpg%nCi
            sna0_driftloss(iCv,1,is) = -1.0_R8/tauion_drift(iCv,is)*
     &                        geo%cvVol(iCv)
            shi0_driftloss(iCv,3) = -1.0_R8/(tauion_drift(iCv,is)/
     &                             switch%bohm_gammaD)*geo%cvVol(iCv)
            she0_driftloss(iCv,3) = she0_driftloss(iCv,3) -
     &                              1.0_R8*rt%rza(iCv,is)/
     &                              (tauion_drift(iCv,is)/
     &                              switch%bohm_gammaD)*geo%cvVol(iCv)
          enddo
      enddo

      if (switch%artloss_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        do is = 0,ns-1
          if(is_neutral(is)) cycle
            write (chns,'(i3.3)') is
            call my_out_us(70,nCv,0,sna0_driftloss(:,0,is),
     &           'b2mod_artificialsinks_sna_driftloss_0_species_'//chns)
            call my_out_us(70,nCv,0,sna0_driftloss(:,0,is),
     &           'b2mod_artificialsinks_sna_driftloss_1_species_'//chns)
        enddo
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,shi0_driftloss(:,cntr),
     &                   'b2mod_artificialsinks_shi_driftloss_'//chns)
          call my_out_us(70,nCv,0,she0_driftloss(:,cntr),
     &                   'b2mod_artificialsinks_she_driftloss_'//chns)
        enddo
      endif

      sna_losses(:,:,:,2) = sna0_driftloss(:,:,:)
      she_losses(:,:,2) = she0_driftloss(:,:)
      shi_losses(:,:,2) = shi0_driftloss(:,:)

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end driftloss

      end subroutine driftloss

c Poloidal losses computation
      subroutine poloidloss(nCv, ns, is,
     &   switch, geo, mpg, st, pl, dv, rt,
     &   sna0_polloss, shi0_polloss, she0_polloss)
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_switches
      implicit none
      
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2State), intent (in) :: st
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      real (kind=R8) :: sna0_polloss(nCv,0:1,0:ns-1),
     &  shi0_polloss(nCv,0:3), she0_polloss(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      integer :: is, iCv
      real (kind=R8) :: tauion_pol(nCv,0:ns-1)
      integer :: cntr
      character chns*3

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('poloidloss')

      call sfill (nCv*ns, 0.0_R8, tauion_pol, 1)
      call sfill (nCv*2*ns, 0.0_R8, sna0_polloss, 1)
      call sfill (nCv*4, 0.0_R8, shi0_polloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_polloss, 1)

      do is = 0, ns-1
        if(is_neutral(is)) cycle
          tauion_pol(:,is) = switch%b**2/(2.0_R8*st%co%dna_bohm(:,0,is))
          do iCv = 1, mpg%nCi
            sna0_polloss(iCv,1,is) = -1.0_R8/tauion_pol(iCv,is)*
     &                               geo%cvVol(iCv)
            shi0_polloss(iCv,3) = -1.0_R8/(tauion_pol(iCv,is)/
     &                            switch%bohm_gammaD)*geo%cvVol(iCv)
            she0_polloss(iCv,3) = she0_polloss(iCv,3) -
     &                            1.0_R8*rt%rza(iCv,is)/
     &                            (tauion_pol(iCv,is)/
     &                            switch%bohm_gammaD)*geo%cvVol(iCv)
          enddo
      enddo

      if (switch%artloss_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        do is = 0,ns-1
          if(is_neutral(is)) cycle
            write (chns,'(i3.3)') is
            call my_out_us(70,nCv,0,sna0_polloss(:,0,is),
     &           'b2mod_artificialsinks_sna_polloss_0_species_'//chns)
            call my_out_us(70,nCv,0,sna0_polloss(:,0,is),
     &           'b2mod_artificialsinks_sna_polloss_1_species_'//chns)
        enddo
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,shi0_polloss(:,cntr),
     &                   'b2mod_artificialsinks_shi_polloss_'//chns)
          call my_out_us(70,nCv,0,she0_polloss(:,cntr),
     &                   'b2mod_artificialsinks_she_polloss_'//chns)
        enddo
      endif

      sna_losses(:,:,:,3) = sna0_polloss(:,:,:)
      she_losses(:,:,3) = she0_polloss(:,:)
      shi_losses(:,:,3) = shi0_polloss(:,:)

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end poloidloss

      end subroutine poloidloss

c Dealloc
      subroutine dealloc_b2mod_artificialsinks
      implicit none
      if (allocated(sna_losses)) then
        deallocate(sna_losses)
      endif
      if (allocated(she_losses)) then
        deallocate(she_losses)
      endif
      if (allocated(shi_losses)) then
        deallocate(shi_losses)
      endif

*   ..return
      return

*-----------------------------------------------------------------------
*.end dealloc_b2mod_artificialsinks
      end subroutine dealloc_b2mod_artificialsinks

*-----------------------------------------------------------------------
*.end b2mod_artificialsinks

      end module b2mod_artificialsinks
