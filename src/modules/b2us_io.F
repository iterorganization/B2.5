      module b2us_io

      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2us_map
      use b2us_geo
      use b2us_plasma

      implicit none

      private

      public :: read_b2fgmtry,  write_b2fgmtry,
     .          read_b2fstate,  write_b2fstate,
     .          read_b2fplasma, write_b2fplasma

      contains

***********************************************************************
*     Module containing basic routines for reading/writing data files
*     in unstructured format.
*     Note: file headers are read elsewhere
***********************************************************************

      subroutine read_b2fgmtry (ninp, mpg, gm)
      implicit none
      integer, intent(in) :: ninp
      type (mapping), intent (inout) :: mpg
      type (geometry), intent (inout) :: gm

* ..read the mapping and geometry

      call read_mapping (ninp, mpg)
      call read_geometry (ninp, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv, gm)

      return
      end subroutine read_b2fgmtry

***********************************************************************

      subroutine write_b2fgmtry (nout, mpg, gm)
      implicit none
      integer :: nout
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: gm

* ..write the mapping and geometry

      call write_mapping (nout, mpg)
      call write_geometry (nout, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, gm)

      return
      end subroutine write_b2fgmtry

***********************************************************************

      subroutine read_b2fstate (ninp, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: ninp, nCv, nFc, ns
      type(B2State), intent(inout) :: st

* ..read the state

      call read_state (ninp, nCv, nFc, ns, st)

      return
      end subroutine read_b2fstate

***********************************************************************

      subroutine write_b2fstate (nout, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: nout, nCv, nFc, ns
      type(B2State), intent(in) :: st

* ..write the state

      call write_state (nout, nCv, nFc, ns, st)

      return
      end subroutine write_b2fstate

***********************************************************************

      subroutine read_b2fplasma (ninp, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: ninp, nCv, nFc, ns
      type(B2State), intent(inout) :: st

* ..read the plasma state

      call read_plasma (ninp, nCv, nFc, ns, st)

      return
      end subroutine read_b2fplasma

***********************************************************************

      subroutine write_b2fplasma (nout, nCv, nFc, ns, st)
      implicit none
      integer, intent(in) :: nout, nCv, nFc, ns
      type(B2State), intent(in) :: st

* ..write the plasma state

      call write_plasma (nout, nCv, nFc, ns, st)

      return
      end subroutine write_b2fplasma

      end module b2us_io
