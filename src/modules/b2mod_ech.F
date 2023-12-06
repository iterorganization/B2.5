c**************************************************************************************
c********************  NEW KU Leuven module for ECH/ECWC plasmas  *********************
c**************************************************************************************
! mb  01/09/23
! main author: Anthony Piras
c   This module is intended to simulate power depositions in Electron Cyclotron Heated
c     plasmas. To turn on the model, the flag b2stbr_ech has to be set to '1'. This
c     module computes an additional heating source due to EC heating injected either
c     at the cyclotron frequency or its harmonics
      module b2mod_ech
*   ..other modules called 
      use b2mod_types
      use b2mod_constants

      implicit none

      real (kind=R8), allocatable, save ::
     &  sECH(:,:) 

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
      subroutine CalcECH(nCv,geo,mpg,switch,she0_ech)
c     This routine calculates the ECH power
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      use b2mod_switches
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (switches), intent (in) :: switch
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: she0_ech(nCv)
c     Other variables
      real (kind=R8) :: omegace(nCv), dummy1(nCv)
      integer :: maskf(nCv)

c     Internal variables
      integer :: iCv, cntr
      character chns*3

      call subini ('CalcECH')
      call sfill (nCv, 0.0_R8, she0_ech, 1)
c      call sfill (nCv, 0.0_R8, omegace, 1)
      call ifill (nCv, 0, maskf, 1)
c      call sfill (nCv, power, she0_ech, 1) useless for now
c      call elcycl(nCv,geo,omegace)

      call sfill(nCv, switch%ech_Rdep, dummy1, 1)

      call maskfun(mpg%nCi, dummy1, geo%cvX,
     &             1.0E-2_R8, maskf) ! 1E-2 is the tolerance

c do loop maybe useless
      do iCv = 1, mpg%nCi
        she0_ech(iCv) = maskf(iCv)*switch%ech_power/sum(maskf)
      enddo

      sECH(:,0) = she0_ech
      if (switch%ech_out.eq.1.or.switch%iout_b2wdat.eq.4) then
        do cntr = 0,3
          write (chns,'(i3.3)') cntr
          call my_out_us(70,nCv,0,sECH(:,cntr),'source_ECH'//chns)
        enddo
      endif

*   ..return
      call subend ()
      return
      end subroutine CalcECH

c Electron Cyclotron Frequency computation
      subroutine elcycl(nCv,geo,omegace)
      use b2mod_b2cmpa
      use b2us_geo
      use b2mod_subsys
      implicit none
      type (geometry), intent (in) :: geo
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out):: omegace(nCv)


      call subini ('elcycl')
      call sfill (nCv, 0.0_R8, omegace, 1)

      omegace(:) = qe*geo%cvBb(:,3)/me

*   ..return
      call subend ()
      return
      end subroutine elcycl

c Dealloc
      subroutine dealloc_b2mod_ech
      implicit none
      if (allocated(sECH)) then
        deallocate(sECH)
      endif

      return
      end subroutine dealloc_b2mod_ech
c End dealloc

      end module b2mod_ech
