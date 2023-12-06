c**************************************************************************************
c********************  NEW KU Leuven module for artificial sinks  *********************
c**************************************************************************************
! mb  10/11/23
! main author: Anthony Piras
c   This module is intended to
      module b2mod_artificialsinks
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
c      use b2mod_wall
c      use b2mod_math
c      use b2mod_types
c      use b2mod_math
c      use b2mod_transport_namelist
c      use b2mod_input_profile
c      use b2mod_constants
c      use b2mod_time
c      use b2mod_diag
c      use b2mod_b2cmpa
c      use b2mod_b2cmpt
c      use b2mod_b2cmpb
c      use b2mod_ranges
c      use b2mod_neoclassical
c      use b2mod_geo_corner
c      use b2mod_user_namelist
c      use b2mod_switches
c      use b2us_geo
c      use b2us_map
c      use b2us_plasma

      implicit none

      real (kind=R8), allocatable, save ::
     &  sna_pol(:,:,:), sna_lim(:,:,:), sna_dri(:,:,:)

      contains

c Alloc
      subroutine alloc_b2mod_artificialsinks(nCv, ns)
      implicit none
      integer, intent(in) :: nCv, ns

      allocate(sna_pol(nCv,0:1,ns))
      allocate(sna_lim(nCv,0:1,ns))
      allocate(sna_dri(nCv,0:1,ns))

      sna_pol(:,:,:) = 0.0_R8
      sna_lim(:,:,:) = 0.0_R8
      sna_dri(:,:,:) = 0.0_R8

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

      call sfill(nCv*2*ns, 0.0_R8, sna0_limloss, 1)
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
              sna0_limloss(iCv,1,is) = -1.0_R8/tauion_lim(iCv,is)
              shi0_limloss(iCv,3) = -1.0_R8/(tauion_lim(iCv,is)/
     &                               switch%bohm_gammaD)
              she0_limloss(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                        (tauion_lim(iCv,is)/switch%bohm_gammaD)
            endif
          endif
        enddo
      enddo

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
            sna0_driftloss(iCv,1,is) = -1.0_R8/tauion_drift(iCv,is)
            shi0_driftloss(iCv,3) = -1.0_R8/(tauion_drift(iCv,is)/
     &                             switch%bohm_gammaD)
            she0_driftloss(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                      (tauion_drift(iCv,is)/switch%bohm_gammaD)
          enddo
        endif
      enddo

*   ..return
      call subend ()
      return
      end subroutine driftloss

c Poloidal losses computation
      subroutine poloidloss(nCv, ns, is,
     &   switch, geo, mpg, pl, dv, rt,
     &   sna0_pollos, shi0_pollos, she0_pollos)
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
     &  sna0_pollos(nCv,0:1,0:ns-1),
     &  shi0_pollos(nCv,0:3), she0_pollos(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      integer :: is, iCv
      real (kind=R8) :: rholi(nCv), mass, vthi(nCv),
     &  charge(nCv,0:ns-1), vthn(nCv), mfpi(nCv), nui(nCv),
     &  nuab(nCv,0:ns,0:ns), nua(nCv,0:ns), tauion_pol(nCv,0:ns-1)
*   ..local variables for output

      call subini ('poloidloss')

      call sfill (nCv*2*ns, 0.0_R8, sna0_pollos, 1)
      call sfill (nCv*ns, 0.0_R8, tauion_pol, 1)
      call sfill (nCv, 0.0_R8, dna_bohm_vert, 1)
      call sfill (nCv*4, 0.0_R8, shi0_pollos, 1)
      call sfill (nCv*4, 0.0_R8, she0_pollos, 1)

      call b2tnua(nCv, ns, switch, pl, dv, rt, nuab, nua)

      do is = 0, ns-1
        if(.not.is_neutral(is)) then
          mass = am(is)*mp
          vthi = sqrt(2.0_R8*pl%ti/mass)
          charge(:,is) = rt%rza(:,is)*qe

*   ..compute local larmor radius
          rholi(:) = mass*vthi/abs(charge(:,is)*geo%cvBb(:,3))
c          call sfill (nCv, 4.03E+5_R8, nui, 1)
          mfpi(:) = sqrt((pl%te + pl%ti*pl%na(:,is)/dv%ne)/mass)/
     &              nua(:,is)

*   ..compute dna_bohm_vert
          dna_bohm_vert(:) = dna_bohm_vert + 1.00_R8/3.00_R8*nua(:,is)*
     &                  mfpi*(rholi+mfpi*geo%cvBb(:,0)/1.5_R8)*0.659_R8
          tauion_pol(:,is) = switch%b**2/(2.0_R8*dna_bohm_vert(:))
          do iCv = 1, mpg%nCi
            sna0_pollos(iCv,1,is) = -1.0_R8/tauion_pol(iCv,is)
            shi0_pollos(iCv,3) = -1.0_R8/(tauion_pol(iCv,is)/
     &                               switch%bohm_gammaD)
            she0_pollos(iCv,3) = -1.0_R8*rt%rza(iCv,is)/
     &                        (tauion_pol(iCv,is)/switch%bohm_gammaD)
          enddo
        endif
      enddo

*   ..return
      call subend ()
      return
      end subroutine poloidloss

c Dealloc
      subroutine dealloc_b2mod_artificialsinks
      implicit none
      if (allocated(sna_pol)) then
        deallocate(sna_pol)
      endif
      if (allocated(sna_lim)) then
        deallocate(sna_lim)
      endif
      if (allocated(sna_dri)) then
        deallocate(sna_dri)
      endif

      return
      end subroutine dealloc_b2mod_artificialsinks
c End dealloc

      end module b2mod_artificialsinks
