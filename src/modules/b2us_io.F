      module b2us_io

      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2mod_geometry, only : GEOMETRY_GENERAL
      use b2us_map
      use b2us_geo
      use b2us_plasma

      implicit none

      private

      public :: read_b2fgmtry,  write_b2fgmtry,
     .          read_b2fstate,  write_b2fstate,
     .          read_b2fplasma, write_b2fplasma,
     .          read_b2fgmtry_header, write_b2fgmtry_header,
     .          read_b2fgmtry_topology_header

      contains

***********************************************************************
*     Module containing basic routines for reading/writing data files
*     in unstructured format.
***********************************************************************

      subroutine reset_b2fgmtry_topology_header(m)
      implicit none
      type(mapping), intent(inout) :: m

* ..Legacy .000 and .001 files contain no explicit topology metadata.
      m%hasTopologicalData = .false.
      m%topoID = 0
      m%geometryID = GEOMETRY_GENERAL
      m%nXpt = 0
      m%nOpt = 0
      m%nStr = 0
      m%nTgc = 0
      m%nDiv = 0
      m%nfcDiv = 0

      return
      end subroutine reset_b2fgmtry_topology_header

***********************************************************************

      subroutine read_b2fgmtry_topology_header(ninp, version, m)
      implicit none
      integer, intent(in) :: ninp
      character(10), intent(in) :: version
      type(mapping), intent(inout) :: m
      integer :: idum(0:5)

      call reset_b2fgmtry_topology_header(m)

* ..The complete topology header is required beginning with .002.
      if (version.lt.'03.002.002') return

      idum(0) = -1
      call cfruin_opt (ninp, 1, idum, 'topoflag')
      if (idum(0).eq.-1) call xerrab (
     &  'incomplete 03.002.002 b2fgmtry topology header; ' //
     &  'required records are topoflag, GEOMETRY_ID, and ' //
     &  'nX,nO,nS,nT,nDiv,nDivFc; regenerate the grid with ' //
     &  'current GOAT or B2.5')
      m%topoID = idum(0)
      call xertst (m%topoID.ge.0.and.m%topoID.le.4,
     &  'invalid b2fgmtry topoflag; expected 0 through 4')

      idum(0) = -huge(idum(0))
      call cfruin_opt (ninp, 1, idum, 'GEOMETRY_ID')
      if (idum(0).eq.-huge(idum(0))) call xerrab (
     &  'incomplete 03.002.002 b2fgmtry topology header; ' //
     &  'required records are topoflag, GEOMETRY_ID, and ' //
     &  'nX,nO,nS,nT,nDiv,nDivFc; regenerate the grid with ' //
     &  'current GOAT or B2.5')
      m%geometryID = idum(0)

      idum = -1
      call cfruin_opt (ninp, 6, idum,
     &  'nX,nO,nS,nT,nDiv,nDivFc')
      if (any(idum.eq.-1)) call xerrab (
     &  'incomplete 03.002.002 b2fgmtry topology header; ' //
     &  'required records are topoflag, GEOMETRY_ID, and ' //
     &  'nX,nO,nS,nT,nDiv,nDivFc; regenerate the grid with ' //
     &  'current GOAT or B2.5')
      call xertst (all(idum.ge.0),
     &  'negative count in b2fgmtry topology header')

      m%hasTopologicalData = .true.
      m%nXpt = idum(0)
      m%nOpt = idum(1)
      m%nStr = idum(2)
      m%nTgc = idum(3)
      m%nDiv = idum(4)
      m%nfcDiv = idum(5)

      return
      end subroutine read_b2fgmtry_topology_header

***********************************************************************

      subroutine write_b2fgmtry_topology_header(nout, m)
      implicit none
      integer, intent(in) :: nout
      type(mapping), intent(in) :: m
      integer :: idum(0:5)

      call xertst (m%hasTopologicalData,
     &  'cannot write .002 b2fgmtry without topology metadata')
      call xertst (m%topoID.ge.0.and.m%topoID.le.4,
     &  'invalid explicit topology flag; expected 0 through 4')
      idum(0) = m%topoID
      call cfwuin (nout, 1, idum, 'topoflag')

      idum(0) = m%geometryID
      call cfwuin (nout, 1, idum, 'GEOMETRY_ID')

      idum(0) = m%nXpt
      idum(1) = m%nOpt
      idum(2) = m%nStr
      idum(3) = m%nTgc
      idum(4) = m%nDiv
      idum(5) = m%nfcDiv
      call xertst (all(idum.ge.0),
     &  'negative count in b2fgmtry topology header')
      call cfwuin (nout, 6, idum, 'nX,nO,nS,nT,nDiv,nDivFc')

      return
      end subroutine write_b2fgmtry_topology_header

***********************************************************************

      subroutine read_b2fgmtry_header(ninp, version, m, lblgm, isymm)
      implicit none
      integer, intent(in) :: ninp
      character(10), intent(in) :: version
      type(mapping), intent(inout) :: m
      character(120), intent(out) :: lblgm
      integer, intent(out) :: isymm
      integer :: idum(0:9)

      if (version.lt.'03.002.000') then
        call xerrab('read_b2fgmtry_header: version older than ' //
     &    '03.002.000; convert the case with b2us first')
      endif

      call cfruin (ninp, 7, idum, 'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      m%nCi = idum(0)
      m%nCg = idum(1)
      m%nCv = idum(2)
      m%nFc = idum(3)
      m%nVx = idum(4)
      m%nFs = idum(5)
      m%nFt = idum(6)

      call cfruin (ninp, 5, idum,
     &  'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      m%nCmxVx = idum(0)
      m%nCmxFc = idum(1)
      m%nVmxCv = idum(2)
      m%nVmxFc = idum(3)
      m%nCmxNv = idum(4)

      call cfruin (ninp, 1, idum, 'isClassicalGrid')
      m%isClassicalGrid = idum(0)
      call cfruin_opt (ninp, 3, idum, 'nx,ny,nncut')
      if (m%isClassicalGrid.eq.1) then
        m%nx = idum(0)
        m%ny = idum(1)
        m%nncut = idum(2)
      else
        m%nx = -1
        m%ny = -1
        m%nncut = -1
      endif

      call read_b2fgmtry_topology_header(ninp, version, m)

      call cfruch (ninp, 120, lblgm, 'label')
      call cfruin (ninp, 1, idum, 'isymm')
      isymm = idum(0)

      return
      end subroutine read_b2fgmtry_header

***********************************************************************

      subroutine write_b2fgmtry_header(nout, m, lblgm, isymm)
      implicit none
      integer, intent(in) :: nout
      type(mapping), intent(in) :: m
      character(120), intent(in) :: lblgm
      integer, intent(in) :: isymm
      integer :: idum(0:9)

      call cfverw (nout, newversion)

      idum(0) = m%nCi
      idum(1) = m%nCg
      idum(2) = m%nCv
      idum(3) = m%nFc
      idum(4) = m%nVx
      idum(5) = m%nFs
      idum(6) = m%nFt
      call cfwuin (nout, 7, idum, 'nCi,nCg,nCv,nFc,nVx,nFs,nFt')

      idum(0) = m%nCmxVx
      idum(1) = m%nCmxFc
      idum(2) = m%nVmxCv
      idum(3) = m%nVmxFc
      idum(4) = m%nCmxNv
      call cfwuin (nout, 5, idum,
     &  'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')

      idum(0) = m%isClassicalGrid
      call cfwuin (nout, 1, idum, 'isClassicalGrid')
      if (m%isClassicalGrid.eq.1) then
        idum(0) = m%nx
        idum(1) = m%ny
        idum(2) = m%nncut
      else
        idum(0) = -1
        idum(1) = -1
        idum(2) = -1
      endif
      call cfwuin (nout, 3, idum, 'nx,ny,nncut')

* ..Do not emit the new records until the global version is advanced.
      if (newversion.ge.'03.002.002') then
        call write_b2fgmtry_topology_header(nout, m)
      endif

      call cfwuch (nout, 120, lblgm, 'label')
      idum(0) = isymm
      call cfwuin (nout, 1, idum, 'isymm')

      return
      end subroutine write_b2fgmtry_header

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
