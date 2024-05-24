*-----------------------------------------------------------------------
*.specification
c**************************************************************************************
c**********************  NEW KU Leuven module for ECH plasmas  ************************
c**************************************************************************************
! mb  01/09/23
! main author: Anthony Piras
      module b2mod_ech
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2us_plasma

      implicit none
*   ..declarations allocatable variable
      real (kind=R8), allocatable, save ::
     &  sECH(:,:)

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MOD_ECH is intended to simulate power depositions in Electron Cyclotron Heated
*     plasmas. To turn on the model, the flag 'b2sral_use_ECH' has to be set to '1'.
*     This module computes an additional heating source due to EC heating injected
*     either at the cyclotron frequency or its harmonics.
*
*-----------------------------------------------------------------------
*-----------------------------------------------------------------------

      contains

c Alloc
      subroutine alloc_b2mod_ech(nCv)
      implicit none
      integer, intent(in) :: nCv

      allocate(sECH(nCv,0:3))

      sECH(:,:) = 0.0_R8

      return
      end subroutine alloc_b2mod_ech
c End alloc

c ECH power computation
      subroutine CalcECH(nCv,geo,mpg,dv,switch,she0_ech)
c     This routine calculates the ECH power
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Derivatives), intent (in) :: dv
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: she0_ech(nCv,0:3)
c     Other variables
      real (kind=R8) :: omegace(nCv), Rdep_vect(nCv)
      integer :: maskf(nCv)
c     Internal variables
      real (kind = R8) :: eprof(nCv), sx, sy, inteprof, normlz, h,
     &                    summation
      real (kind = R8) :: segment(2,2)
      integer :: list_cv(nCv), nncv, list_ci(mpg%nCi), nnci
      integer :: iCv, iCv2, cntr
      character chns*3
      logical:: ifail

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('CalcECH')
      call sfill (nCv*4, 0.0_R8, she0_ech, 1)
      call ifill (nCv, 0, maskf, 1)
      call sfill (nCv, switch%ech_Rdep, Rdep_vect, 1)
      call sfill (nCv, 0.0_R8, eprof, 1)
      call sfill (4, 0.0_R8, segment, 1)
      inteprof = 0.0_R8
      summation = 0.0_R8
      normlz = 0.0_R8

      select case (switch%ech_coupling)

* -- COUPLING=0 -- Fixed power distributed as a Dirac's delta in volumes intersected at given radial position. !!TODO Generalise to 2D

      case(0)
        segment = reshape((/ switch%ech_Rdep, 0.0_R8,
     &       switch%ech_Rdep, maxval(geo%cvY) /), shape(segment))

        call cv_intersections (geo, mpg, segment, list_cv, list_ci,
     &   nncv, nnci, nCv, 1, ifail)
        do cntr = 1, size(list_ci)
          if(list_ci(cntr).ne.0) then
            iCv = list_ci(cntr)
            she0_ech(iCv,0) = switch%ech_absfract*switch%ech_power ! absorption fraction to be generalised with a switch lookup tables
          endif
        enddo

* -- COUPLING=1 -- Fixed fraction of launched power is coupled/Gaussian injection. !!TODO Check implementation

      case(1)
        sx = 1.8667E-2_R8 ! widthech in T1D -> to be generalised
        sy = sx ! assumed equal for now for 2D
        if (switch%ech_gaussian_symmetry.eq.1) then

          segment = reshape((/ 0.0_R8, maxval(geo%cvY)/2.0_R8,
     &       maxval(geo%cVx), maxval(geo%cvY)/2.0_R8 /), shape(segment))

          call cv_intersections (geo, mpg, segment, list_cv,
     &         list_ci, nncv, nnci, nCv, 1, ifail)

          do cntr = 1, nnci
            iCv = list_ci(cntr)
            eprof(iCv) = geo%cvX(iCv)*exp(-1.0_R8/2.0_R8*
     &            ((geo%cvX(iCv)-switch%ech_Rdep)/sx)**2) ! symmetry on poloidal direction
          enddo

          do cntr = 1, size(list_ci)-1
            if(list_ci(cntr).ne.0.and.list_ci(cntr+1).ne.0) then
              iCv = list_ci(cntr)
              iCv2 = list_ci(cntr+1)
              inteprof = inteprof+2.0_R8*switch%b*pi*0.5_R8*(eprof(iCv)+
     &                   eprof(iCv2))*(geo%cvX(iCv2)**2-geo%cvX(iCv)**2)
            endif
            if(cntr.ne.1.and.cntr.ne.(size(list_ci)-1)) then
              summation = summation + eprof(iCv)
            endif
          enddo

          h = (maxval(geo%cvX)-minval(geo%cvX))/nnci
          eprof = eprof/inteprof
          normlz = h*(eprof(list_ci(1))+2.0_R8*summation+
     &             eprof(list_ci(nnci)))/2.0_R8
          eprof = eprof/inteprof/normlz

        else
          eprof = 1.0_R8/(2*pi*sx*sy)*exp(-1.0_R8/2.0_R8*
     &            ((geo%cvX-switch%ech_Rdep)**2/sx**2 +
     &             (geo%cvY-switch%ech_Zdep)**2/sy**2)) ! ech_Zdep to be defined
        endif

        do iCv = 1, mpg%nCi
          she0_ech(iCv,0) = switch%ech_power*eprof(iCv)*
     &                      switch%ech_absfract
        enddo

* -- COUPLING=2 -- TOMAS - deposition on ECR (Budden analysis). !!TODO

      case(2)
        call xerrab ('b2mod_ech -- ECH_COUPLING = 2 not coded yet')

* -- COUPLING=3 -- Ray-tracing with TORAY. !!TODO

      case(3)
        call xerrab ('b2mod_ech -- ECH_COUPLING = 3 not coded yet')

      end select

      sECH = she0_ech
      if (switch%ech_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,sECH(:,cntr),'b2mod_ECH_she_'//chns)
        enddo
      endif

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end CalcECH

      end subroutine CalcECH

c Electron Cyclotron Frequency computation
      subroutine elcyclfr(nCv,geo,omegace)
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      type (geometry), intent (in) :: geo
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: omegace(nCv)

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('elcyclfr')
      call sfill (nCv, 0.0_R8, omegace, 1)

      omegace(:) = qe*geo%cvBb(:,3)/me

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end elcyclfr

      end subroutine elcyclfr

c Dealloc
      subroutine dealloc_b2mod_ech
      implicit none
      if (allocated(sECH)) then
        deallocate(sECH)
      endif

      return
      end subroutine dealloc_b2mod_ech
c End dealloc

*-----------------------------------------------------------------------
*.end b2mod_ech

      end module b2mod_ech
