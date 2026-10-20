! Geometry-only regression for a legacy structured linear grid.
!
! Build with the optional Makefile target linear_topology_test, then run in
! a scratch directory containing a structured b2fgmtry.  The driver
! chooses an interior flux surface as b2mwti_jsep, converts the geometry to
! unstructured form, and verifies the metadata-free round trip.  It does not
! require b2mn.dat, b2fpardf, or a plasma state.
program check_linear_topology_roundtrip
  use b2mod_types
  use b2mod_version
  use b2mod_geo
  use b2mod_indirect
  use b2mod_b2cmfs
  use b2mod_ipmain
  use b2mod_xerset
  use b2mod_geometry, only: GEOMETRY_LINEAR
  use b2us_map
  use b2us_geo
  use b2us_prep, only: calc_mapping, prep_geometry, finalize_legacy_topology
  use b2us_io
  implicit none

  type(mapping) :: original, reloaded
  type(geometry) :: original_geo, reloaded_geo
  character(120) :: lblgm, reloaded_label
  character(32) :: reference_text
  integer :: idum(0:9), nx, ny, reference_fs, i, reloaded_isymm
  integer, parameter :: input_unit = 51, output_unit = 61, reload_unit = 71
  external :: b2rugm, cfruin, cfruch, cfopen, xertst

  call prgini('check_linear_topology_roundtrip')
  call cfopen(input_unit, 'b2fgmtry', 'old', 'un*formatted')
  call xerset(0)
  call cfverr(input_unit, b2fgmtry_version)
  call xertst(b2fgmtry_version.lt.'03.002.000', &
       'check_linear_topology_roundtrip: expected structured input')

  call cfruin(input_unit, 2, idum, 'nx,ny')
  nx = idum(0)
  ny = idum(1)
  call alloc_b2mod_geo(nx, ny)
  call alloc_b2mod_indirect(nx, ny, nncutmax)
  call cfruch(input_unit, 120, lblgm, 'label')
  call cfruin(input_unit, 1, idum, 'isymm')
  isymm = idum(0)
  call b2rugm(input_unit, nx, ny, isymm)
  hz = vol/gs(:,:,2)

  call calc_mapping(nx, ny, original)
  call xertst(maxval(original%cvReg).eq.1, &
       'check_linear_topology_roundtrip: expected one cell region')

  ! Pick a nonempty interior surface so the legacy mapper constructs the
  ! two artificial endpoint strikes without introducing an X-point.
  reference_fs = 0
  do i = 2, original%nFs-1
    if (original%fsFcP(i,2).gt.1) then
      reference_fs = i
      exit
    endif
  enddo
  call xertst(reference_fs.gt.0, &
       'check_linear_topology_roundtrip: no interior flux surface')
  write(reference_text, '(i0)') reference_fs
  call ipsetc('b2mwti_jsep', trim(reference_text))

  call prep_geometry(nx, ny, original_geo, original)
  call finalize_legacy_topology(original, original_geo)
  call xertst(original%geometryID.eq.GEOMETRY_LINEAR, &
       'check_linear_topology_roundtrip: geometry is not linear')
  call xertst(.not.original%hasTopologicalData, &
       'check_linear_topology_roundtrip: endpoint strikes became explicit')
  call xertst(original%nXpt.eq.0.and.original%nStr.eq.2, &
       'check_linear_topology_roundtrip: expected two endpoints and no X-point')
  call xertst(allocated(original%divFcP), &
       'check_linear_topology_roundtrip: target pointers missing')
  call xertst(size(original%divFcP,1).eq.2, &
       'check_linear_topology_roundtrip: expected two targets')

  call cfopen(output_unit, 'b2fgmtry.linear_test', 'new', 'formatted')
  call write_b2fgmtry_header(output_unit, original, lblgm, isymm)
  call write_b2fgmtry(output_unit, original, original_geo)
  close(output_unit)

  call cfopen(reload_unit, 'b2fgmtry.linear_test', 'old', 'un*formatted')
  call cfverr(reload_unit, b2fgmtry_version)
  call xertst(b2fgmtry_version.eq.'03.002.001', &
       'check_linear_topology_roundtrip: expected 03.002.001 output')
  call read_b2fgmtry_header(reload_unit, b2fgmtry_version, &
       reloaded, reloaded_label, reloaded_isymm)
  call alloc_mapping(reloaded)
  call alloc_geometry(reloaded_geo, reloaded%nCv, reloaded%nFc, &
       reloaded%nVx, reloaded%nFs, reloaded%nFt, reloaded%nVmxCv)
  call read_b2fgmtry(reload_unit, reloaded, reloaded_geo)

  call xertst(reloaded%nXpt.eq.original%nXpt.and. &
       reloaded%nStr.eq.original%nStr, &
       'check_linear_topology_roundtrip: endpoint count changed')
  call xertst(reloaded%iFssep.eq.original%iFssep, &
       'check_linear_topology_roundtrip: reference surface changed')
  call xertst(allocated(reloaded%divFcP), &
       'check_linear_topology_roundtrip: reloaded targets missing')
  call xertst(all(reloaded%strVx.eq.original%strVx).and. &
       all(reloaded%strDiv.eq.original%strDiv), &
       'check_linear_topology_roundtrip: endpoint mapping changed')
  call xertst(all(shape(reloaded%divFcP).eq.shape(original%divFcP)), &
       'check_linear_topology_roundtrip: target count changed')
  call xertst(all(reloaded%divFcP.eq.original%divFcP), &
       'check_linear_topology_roundtrip: target pointers changed')
  call xertst(size(reloaded%divFc).eq.size(original%divFc), &
       'check_linear_topology_roundtrip: target face count changed')
  call xertst(all(reloaded%divFc.eq.original%divFc), &
       'check_linear_topology_roundtrip: target faces changed')

  write(*,'(a,i0,a,i0,a,i0)') 'PASS: reference surface ', reference_fs, &
       ', endpoints ', reloaded%nStr, ', targets ', size(reloaded%divFcP,1)
end program check_linear_topology_roundtrip
