c**************************************************************************************
c********************  NEW KU Leuven module for artificial sinks  *********************
c**************************************************************************************
! mb  10/11/23
! main author: Anthony Piras
c This module is intended to simulate artificial sinks due to losses in the poloidal direction,
c drifts and eventual losses to limiters. These sinks are used when "tweaking" SOLPS-ITER to
c run 1D models (i.e. symmetry BCs applied in poloidal direction, purely radial model).
c Subroutines contained in this module should not be used when running a full 2D simulation.
c Verify that losses_limiters, losses_drifts and losses_poloidal switches have their default
c value (zero), thus being turned off.
c This losses model comes from the TOMATOR1D code (credits: Tom Wauters tom.wauters@iter.org).
c The losses are stored in vectors sna0_XXXloss, she0_XXXloss and shi0_XXXloss, all of them of
c the same size of st%sr%sna/she/shi.
c Particle losses are assumed to be of the form n_i/tau, hence only sna0_XXXloss(:,1,:) is
c non-zero. The subroutines return sna0_XXXloss(:,1,:) = -1.0/tau, then the multiplication by
c n_i is carried out in b2sral.
c Energy losses are assumed to be of the form n*T/tau, hence only shY0_XXXloss(:,3) is non-zero.
 
      module b2mod_artificialsinks
*   ..other modules called 
      use b2mod_types
      use b2mod_constants

      implicit none

c The allocatable variables hereunder have the same dimensions of the ones defined in b2us_plasma
c plus an extra dimension equal to 3. The logic is as follows:
c sna_losses(:,:,:,1) = sna0_limloss(:,:,:)   for limiter losses
c sna_losses(:,:,:,2) = sna0_driftloss(:,:,:) for drift losses
c sna_losses(:,:,:,3) = sna0_polloss(:,:,:)   for poloidal losses
c and so on and so forth for the other variables
      real (kind=R8), allocatable, save ::
     &  sna_losses(:,:,:,:), ! particles losses
     &  she_losses(:,:,:),   ! electron energy losses
     &  shi_losses(:,:,:)    ! ion energy losses

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
     &   switch, geo, mpg, pl, dv, rt, 
     &   sna0_limloss, shi0_limloss, she0_limloss)
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
*   ..input/output arguments
      real (kind=R8) :: sna0_limloss(nCv,0:1,0:ns-1),
     &  shi0_limloss(nCv,0:3), she0_limloss(nCv,0:3)
*   ..workspace arguments
      integer :: is, iCv
      real (kind=R8) :: c_sound, tauion_lim(nCv,0:ns-1),
     &  Rp_min, Rp_max

      call subini ('limiterloss')

      call sfill (nCv*2*ns, 0.0_R8, sna0_limloss, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_lim, 1)
      call sfill (nCv*4, 0.0_R8, shi0_limloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_limloss, 1)

      Rp_min = switch%R0 - switch%limiter_width * switch%a0;
      Rp_max = switch%R0 + switch%limiter_width * switch%a0;

      do iCv = 1, mpg%nCi
        do is = 0, ns-1
          if(is_neutral(is)) cycle
          if((geo%cvX(iCv).le.Rp_min.or.geo%cvX(iCv).ge.Rp_max)) then
            c_sound = sqrt(1.6667*rt%rza(iCv,is)*pl%te(iCv)/(am(is)*mp))

            tauion_lim(iCv,is) = 0.66_R8*2.0_R8*pi*geo%cvX(iCv)/ ! 0.66 numerical factor by Tom (see emails)
     &                        switch%nlimiters/c_sound
            if(tauion_lim(iCv,is).ne.0.0_R8) then
              sna0_limloss(iCv,1,is) = -1.0_R8/tauion_lim(iCv,is)*
     &                        geo%cvVol(iCv)
              shi0_limloss(iCv,3) = -1.0_R8/(tauion_lim(iCv,is)/
     &                             switch%bohm_gammaD)*geo%cvVol(iCv)
              she0_limloss(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                        (tauion_lim(iCv,is)/switch%bohm_gammaD)*
     &                        geo%cvVol(iCv)
            endif
          endif
        enddo
      enddo

      sna_losses(:,:,:,1) = sna0_limloss(:,:,:)
      she_losses(:,:,1) = she0_limloss(:,:)
      shi_losses(:,:,1) = shi0_limloss(:,:)

*   ..return
      call subend ()
      return
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

      call subini('driftloss')

      call sfill(nCv*2*ns, 0.0_R8, sna0_driftloss, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_drift, 1)
      call sfill (nCv*4, 0.0_R8, shi0_driftloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_driftloss, 1)

      do is = 0, ns-1
        if(.not.is_neutral(is)) then
          charge(:,is) = rt%rza(:,is)*qe
          tauion_drift(:,is) = switch%a0/(2.0*pl%ti/charge(:,is)/
     &                         geo%cvBb(:,3)/geo%cvX(:))
          do iCv = 1, mpg%nCi
            sna0_driftloss(iCv,1,is) = -1.0_R8/tauion_drift(iCv,is)*
     &                        geo%cvVol(iCv)
            shi0_driftloss(iCv,3) = -1.0_R8/(tauion_drift(iCv,is)/
     &                             switch%bohm_gammaD)*geo%cvVol(iCv)
            she0_driftloss(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                      (tauion_drift(iCv,is)/switch%bohm_gammaD)*
     &                        geo%cvVol(iCv)
          enddo
        endif
      enddo

      sna_losses(:,:,:,2) = sna0_driftloss(:,:,:)
      she_losses(:,:,2) = she0_driftloss(:,:)
      shi_losses(:,:,2) = shi0_driftloss(:,:)

*   ..return
      call subend ()
      return
      end subroutine driftloss

c Poloidal losses computation
      subroutine poloidloss(nCv, ns, is,
     &   switch, geo, mpg, pl, dv, rt,
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
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      real (kind=R8) :: dna_bohm_vert(nCv),
     &  sna0_polloss(nCv,0:1,0:ns-1),
     &  shi0_polloss(nCv,0:3), she0_polloss(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      integer :: is, iCv
      real (kind=R8) :: rholi(nCv), mass, vthi(nCv),
     &  charge(nCv,0:ns-1), vthn(nCv), mfpi(nCv), nui(nCv),
     &  nuab(nCv,0:ns,0:ns), nua(nCv,0:ns), tauion_pol(nCv,0:ns-1)
*   ..local variables for output

      call subini ('poloidloss')

      call sfill (nCv*2*ns, 0.0_R8, sna0_polloss, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_pol, 1)
      call sfill (nCv, 0.0_R8, dna_bohm_vert, 1)
      call sfill (nCv*4, 0.0_R8, shi0_polloss, 1)
      call sfill (nCv*4, 0.0_R8, she0_polloss, 1)

      call b2tnua(nCv, ns, switch, pl, dv, rt, nuab, nua)

      do is = 0, ns-1
        if(.not.is_neutral(is)) then
          mass = am(is)*mp
          vthi = sqrt(2.0_R8*pl%ti/mass)
          charge(:,is) = rt%rza(:,is)*qe

*   ..compute local larmor radius
          rholi(:) = mass*vthi/abs(charge(:,is)*geo%cvBb(:,3))
          mfpi(:) = sqrt((pl%te + pl%ti*pl%na(:,is)/dv%ne)/mass)/
     &              nua(:,is)

*   ..compute dna_bohm_vert
          dna_bohm_vert(:) = dna_bohm_vert + 1.00_R8/3.00_R8*nua(:,is)*
     &                  mfpi*(rholi+mfpi*geo%cvBb(:,0)/1.5_R8)*0.659_R8
          tauion_pol(:,is) = switch%b**2/(2.0_R8*dna_bohm_vert(:))
          do iCv = 1, mpg%nCi
            sna0_polloss(iCv,1,is) = -1.0_R8/tauion_pol(iCv,is)*
     &                               geo%cvVol(iCv)
            shi0_polloss(iCv,3) = -1.0_R8/(tauion_pol(iCv,is)/
     &                            switch%bohm_gammaD)*geo%cvVol(iCv)
            she0_polloss(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                         (tauion_pol(iCv,is)/switch%bohm_gammaD)*
     &                         geo%cvVol(iCv)
          enddo
        endif
      enddo

      sna_losses(:,:,:,3) = sna0_polloss(:,:,:)
      she_losses(:,:,3) = she0_polloss(:,:)
      shi_losses(:,:,3) = shi0_polloss(:,:)

*   ..return
      call subend ()
      return
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

      return
      end subroutine dealloc_b2mod_artificialsinks
c End dealloc

      end module b2mod_artificialsinks
