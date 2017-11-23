!!  Legend:
!!     !> ................ Documentation comment (file description, function
!!                         description etc.). Also intended for doxygen
!!                         generated documentation
!!     !> @note .......... Documentation notes, intended for doxygen generated
!!                         documentation
!!     !!  ............... variables description, additional (helpful)
!!                         information etc.
!!     ! IGNORE    ....... Used to ignore this module in list dependency when
!!                         building
!!     !   ............... Commented part of code
!! -----------------------------------------------------------------------------
!! DOCUMENTATION:
!> 1. purpose
!>
!>      Note:   This script is OUTDATED! The development continued under new
!>              script b2_ual_write_b2mod.
!>              The development 'moved' to new script due to different
!>              approach of writing data to IDS using b2mod routines.
!>
!>      b2_ual_write_gsl.f90 script is used to generate b2_ual_write_gsl.exe
!>      (main program), which is a post-processor for b2.
!>      Same as b2_ual_write.F90 the script it reads the plasma grid geometry
!>      ( including Nodes and Cells grid subsets) and
!>      plasma state (electron density, electron temperature, ion temperature)
!>      and writes it to IDS database but with the use of IMAS GGD Grid
!>      Service Library routines.
!>
!>
!> 2. specification
!>
!>      Main program.
!>
!>
!> 3. description (see also routine b2cdca)
!>
!>      The complete program performs post-processing of the
!>      result of a b2 calculation.
!>      This program unit opens and closes the input/output units, and
!>      may perform some other system-dependent operations.
!>
!>      The input units are:
!>      ninp(0): formatted; provides output control parameters.
!>      ninp(1): un*formatted; provides the geometry.
!>      ninp(2): un*formatted; provides the run parameters.
!>      ninp(3): un*formatted; provides the plasma state.
!>      ninp(4): unformatted; provides the detailed plasma state.
!>      ninp(5): formatted; provides the run switches.
!>      ninp(6): un*formatted; provides the atomic data.
!>
!>      The output units are:
!>      nout(0): formatted; provides printed output.
!>
!>      (See routine b2cdca for the meaning of 'un*formatted'.)
!>
!>
!> 4. parameters (see also routine b2cdcv)
!>
!>      None.
!>
!>
!> 5. error indicators
!>
!>      In case an error condition is detected, a call is made to the
!>      routine xerrab. This causes an error message to be printed,
!>      after which the program halts.
!>
!! -----------------------------------------------------------------------------

program b2_ual_write_gsl
    use b2mod_types , B2R8 => R8
    ! use b2mod_version
    ! use b2mod_geo
    ! use b2mod_plasma
    use b2mod_rates
    ! use b2mod_residuals
    ! use b2mod_sources
    ! use b2mod_transport
    ! use b2mod_anomalous_transport
    ! use b2mod_work
    ! use b2mod_time
    ! use b2mod_ppout
    ! use b2mod_tallies
    ! use b2mod_indirect
    ! use b2mod_neutrals_namelist
    ! use b2mod_neutr_src_scaling
    ! use b2mod_b2cmfs
    ! use b2mod_constants
    ! use b2mod_grid_mapping
    ! use b2mod_ual_io_grid
    ! use b2mod_ual_io_data
    ! use b2mod_ual
    use b2mod_ual_io
    ! use b2mod_geo_corner
    ! use b2mod_b2cmrc
    ! use b2mod_b2cmgs
    ! use b2mod_b2cmpa
    ! use b2mod_b2cmpb
    ! use b2mod_b2cmpt
    ! use b2mod_b2cmwg

    use ids_schemas     ! IGNORE
                        !> These are the Fortran type definitions for the
                        !> Physics Data Model
    use ids_routines    ! IGNORE
                        !> These are the Access Layer routines + management of
                        !> IDS structures
    use ids_assert      ! IGNORE
    use ids_grid_common &       ! IGNORE
        & , IDS_COORDTYPE_R => COORDTYPE_R    &
        & , IDS_COORDTYPE_Z => COORDTYPE_Z
    use ids_string              ! IGNORE
    use ids_grid_subgrid        ! IGNORE
    use ids_grid_objectlist     ! IGNORE
    use ids_grid_examples       ! IGNORE
    use ids_grid_unstructured   ! IGNORE
    use ids_grid_structured     ! IGNORE

    implicit none

    !>--------------------------------------------------------------------------

    !>.declarations

    !>..common blocks

    !>..local variables
    integer ninp(0:6), nout(0:2), nx, ny, ns, idum(0:9), io
    real (kind=B2R8), allocatable   ::  ne1(:), te1(:), ti1(:)
    integer     ::  cf_size1(22)
    integer     ::  idx
    integer     ::  ix, ixx
    integer     ::  shot, run
    character(len=24)           ::  treename, username, device, version
    !> character(len=255)          ::  imas_connect_url
    type(ids_edge_profiles)     ::  edge_profiles
    !>  ..procedures
    external prgini, prgend, xerset, xertst, cfopen, cfruin
    !>  ..initialize input/output units
    data ninp/50, 51, 52, 53, 54, 55, 56/
    data nout/60, 61, 62/

    !>--------------------------------------------------------------------------
    !>.documentation-internal
    !>
    !>      The following common blocks have their outermost declaration in
    !>      this routine; they need not be preserved between calls.
    !>
    !>      Description of some local variables:
    !>
    !>      ninp - (0:6) integer array.
    !>      ninp specifies the input unit numbers.
    !>
    !>      nout - (0:2) integer array.
    !>      nout specifies the output unit numbers.
    !>
    !>      nx, ny - integer.
    !>      nx and ny specify the number of interior cells along the first
    !>      and the second coordinate, respectively. The total number of
    !>      cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
    !>      It will hold that 0.le.nx and 0.le.ny.
    !>
    !>      ns - integer.
    !>      ns specifies the number of atomic species in the calculation.
    !>      The species are indexed by (0:ns-1).
    !>      It will hold that 1.le.ns.
    !>
    !>--------------------------------------------------------------------------
    !>.computation

    !> Read main data
    call read_b2fgmtry_b2fstate()

    !> Read additional data
    call read_additional(ninp, nout, nx, ny, ns, ne1, te1, ti1)

    treename    = "ids"
    shot        = 16151
    run         = 1001
    username    = "penkod"
    device     = "solps-iter"
    version     = "3"

    !> Set the data to edge_profiles IDS
    call write_ids_edge_profiles(treename, shot, run, idx, username, &
                                        & device, version, ne1, te1, ti1)

    ! call read_ids(treename, shot, run, idx, username, &
    !                                     & device, version)

contains

    subroutine read_b2fgmtry_b2fstate()

        integer                 ::  i, j, k
        character(len=120)      ::  lblgm

        !>..program start_up calls
        call prgini('b2_ual_write')

        !>..open files
        !>...input files
        ! call cfopen (ninp(0),'b2_ual_write.dat','old','formatted')
        call cfopen(ninp(1),'b2fgmtry','old','un*formatted')
        ! call cfopen (ninp(2),'b2fparam','old','un*formatted')
        call cfopen(ninp(3),'b2fstate','old','un*formatted')
        ! call cfopen (ninp(4),'b2fplasma','old','unformatted')
        ! call cfopen (ninp(5),'b2mn.dat','old','formatted')
        ! call cfopen (ninp(6),'b2frates','old','formatted')
        !>...output files
        call cfopen(nout(0),'b2_ual_write.prt','new','formatted')
        call cfopen(nout(1),'b2_ual_writ2.prt','new','formatted')

        !>..mark output unit for error messages
        call xerset(0)
        !>..obtain version numbers
        call cfverr(ninp(1), b2fgmtry_version)
        ! call cfverr (ninp(2), b2fparam_version)
        call cfverr(ninp(3), b2fstate_version)
        ! call cfverr (ninp(6), b2frates_version)
        !>..obtain nx, ny, ns
        call cfruin(ninp(3), 3, idum, 'nx,ny,ns')
        nx = idum(0)
        ny = idum(1)
        ns = idum(2)
        call xertst(0.le.nx.and.0.le.ny.and.1.le.ns, &
            & 'faulty input nx, ny, ns from plasma state file')
        call cfruin (ninp(1), 2, idum(0), 'nx,ny')
        call xertst(idum(0).eq.nx.and.idum(1).eq.ny, &
            & 'faulty input nx, ny from geometry file')

        call alloc_b2mod_geo(nx,ny)
        call alloc_b2mod_plasma(nx,ny,ns)
        call alloc_b2mod_indirect(nx,ny,nncutmax)
        call alloc_b2mod_sources(nx,ny,ns)

        call cfruch (ninp(1), 120, lblgm, 'label')
        call cfruin (ninp(1), 1, idum, 'isymm')
        call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz, &
         & bb, vol, hx, hy, qz, qc, qcb, gs, pbs, &
         & wbbl, wbbr, wbbv, wbbc, cell_width, cell_height, gmap)

        call xertst(size(crx).eq.size(cry), &
            & "The number of x and y coordinates is not the same.")

        corner: do k = 0, 3, 1
            xi: do j = -1, ny, 1
                yi: do i = -1, nx, 1
                    write (nout(0),*) i, j, k, crx(i, j, k)
                end do yi
            end do xi
        end do corner

        call prgend()

        ! stop 'b2_ual_write'

    end subroutine read_b2fgmtry_b2fstate

    subroutine read_additional(ninp, nout, nx, ny, ns, ne, te, ti)

        !> Read and compute additional data (taken from b2mddr.F)

        implicit none

        !>   ..input arguments (unchanged on exit)
        integer ninp(0:6), nout(0:2), nx, ny, ns
        !>   ..output arguments (unspecified on entry)
        !     (none)
        !>   ..common blocks

        !>   ..local variables
        integer k, is, idum(0:9) !>, nscx, iscx(0:nscxmax-1), ismain
        character lblgm*120, lblcp*120, lblmn*120, lblrc*120
        character cnamip*80, cvalip*80
        !>   ..procedures
        external subini, subend, xertst, xerrab, cfruch, cfruin, cfrure
        !>   ..namelist
        character   :: exp*128, comment*128
        integer     :: shot, overwrite_shotnumber
        logical     :: timedep, snapshot, tallies, movies
        real (kind=B2R8), allocatable   ::  zamin(:), zamax(:), zn(:), am(:),   &
            &   na(:), ua(:), uadia(:), po(:), fna(:), fhe(:), fhi(:), fch(:),  &
            &   fch_32(:), fch_52(:), kinrgy(:), fch_p(:)
        real (kind=B2R8), intent(inout), allocatable    :: ne(:), te(:), ti(:)
        real (kind=B2R8)    ::  time
        integer             ::  cf_sizes(22)

        namelist /b2md_namelist/ exp, shot, time, comment, timedep, snapshot, &
            & tallies, movies, overwrite_shotnumber

        !>   ..subprogram start-up calls
        call subini ('b2mddr')
        !>   ..test ninp, nout

        call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.&
             & 1.le.ninp(3).and.1.le.ninp(5).and.1.le.ninp(6).and.&
             & 1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),&
             & 'faulty argument ninp, nout')
        !>   ..test dimensions

        call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
        call xertst (1.le.ns, 'faulty argument ns')

        call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

        !> Get array sizes from the b2fstate. We open the file again as a new
        !> unit to not disrupt the old unit.
        !> A routine should exist to do all that work but so far with no success
        !> of finding/properly use it.
        cf_sizes = read_additional_sizes(22)

        ! allocate(lblgm(cf_sizes(2)))
        allocate(zamin(cf_sizes(3)))
        allocate(zamax(cf_sizes(4)))
        allocate(zn(cf_sizes(5)))
        allocate(am(cf_sizes(6)))
        allocate(na(cf_sizes(7)))
        allocate(ne(cf_sizes(8)))
        allocate(ua(cf_sizes(9)))
        allocate(uadia(cf_sizes(10)))
        allocate(te(cf_sizes(11)))
        allocate(ti(cf_sizes(12)))
        allocate(po(cf_sizes(13)))
        allocate(fna(cf_sizes(14)))
        allocate(fhe(cf_sizes(15)))
        allocate(fhi(cf_sizes(16)))
        allocate(fch(cf_sizes(17)))
        allocate(fch_32(cf_sizes(18)))
        allocate(fch_52(cf_sizes(19)))
        !> The program gives an error saying the kinrgy and fch_p are
        !> already allocated. Where, when?
        ! allocate(kinrgy(cf_sizes(20)))
        ! allocate(fch_p(cf_sizes(22)))

        ! call cfruch(ninp(3), cf_sizes(2), lblgm, 'label')
        call cfruch(ninp(3), 120, lblgm, 'label')
        call cfrure(ninp(3), cf_sizes(3), zamin, 'zamin')
        call cfrure(ninp(3), cf_sizes(4), zamax, 'zamax')
        call cfrure(ninp(3), cf_sizes(5), zn, 'zn')
        call cfrure(ninp(3), cf_sizes(6), am, 'am')
        call cfrure(ninp(3), cf_sizes(7), na, 'na')
        call cfrure(ninp(3), cf_sizes(8), ne, 'ne')
        call cfrure(ninp(3), cf_sizes(9), ua,'ua')
        call cfrure(ninp(3), cf_sizes(10), uadia,'uadia')
        call cfrure(ninp(3), cf_sizes(11), te, 'te')
        call cfrure(ninp(3), cf_sizes(12), ti, 'ti')
        call cfrure(ninp(3), cf_sizes(13), po, 'po')
        call cfrure(ninp(3), cf_sizes(14), fna, 'fna')
        call cfrure(ninp(3), cf_sizes(15), fhe, 'fhe')
        call cfrure(ninp(3), cf_sizes(16), fhi, 'fhi')
        call cfrure(ninp(3), cf_sizes(17), fch,'fch')
        call cfrure(ninp(3), cf_sizes(18), fch_32, 'fch_32')
        call cfrure(ninp(3), cf_sizes(19), fch_52, 'fch_52')
        ! call cfrure(ninp(3), cf_sizes(20), kinrgy, 'kinrgy')
        ! call cfrure(ninp(3), 1, time, 'time')
        ! call cfrure(ninp(3), cf_sizes(22), fch_p, 'fch_p')

    end subroutine read_additional

    !> Get array sizes from the b2fstate.
    function read_additional_sizes(array_size)
        integer, intent(in)     ::  array_size
        character(len=256)      ::  rdline, rdstring
        integer                 ::  cfcount, stat
        integer    ::  cf_sizes(array_size)
        integer, dimension(:), allocatable     ::  read_additional_sizes

        allocate(read_additional_sizes(array_size))

        cfcount = 0
        open(1,file='b2fstate', status="old")
        do
            read(1,'(A)', iostat=io) rdline
            if (io/=0) exit
            if (rdline(1:3).eq."*cf") then
                cfcount = cfcount + 1
                rdstring = adjustl(rdline(20:28))
                read(rdstring,*, iostat=stat)  cf_sizes(cfcount)
            end if
        end do
        close(1)
        read_additional_sizes = cf_sizes

    end function read_additional_sizes


    subroutine write_ids_edge_profiles(treename, shot, run, idx, username, &
                                        & device, version, ne, te, ti)
        use ids_grid_unstructured   ! IGNORE

        integer ::  shot, run, idx
        integer ::  num_nodes_all, num_nodes, num_gridSubsets
        integer ::  gridSubset_index
        integer ::  numCellsX, numCellsY, num_cells, cellId
        integer ::  numEdgesX, numEdgesY, num_edges, edgeId, count, count2
        integer ::  num_ne_gridSubset, num_ne_values
        integer ::  num_te_gridSubset, num_te_values
        integer ::  num_ti_gridSubset, num_ti_values, num_ti_species, ion_specie
        integer ::  num_obj_0D, obj_0D_id
        integer ::  num_obj_2D, obj_2D_id
        integer ::  gridSubset_dim_index
        integer ::  i, j, k, icount, n
        integer ::  homogeneous_time
        integer ::  coordtype(2)
        integer, allocatable        ::  edgesNodesList(:,:), cellsNodesList(:,:)
        integer, allocatable        ::  edgeIndicesRepeat(:)
        character(len=24)           ::  treename, username, device, version
        character(len=255)          ::  imas_connect_url
        character(len=255)          ::  grid_description
        character(len=132)          ::  gridSubset_name
        real (kind=B2R8), intent(in)    :: ne(:), te(:), ti(:)
        !> TODO: get time out from b2fstate
        ! real(B2R8)    ::  time
        real(IDS_real)    ::  time
        real(IDS_real), allocatable   ::  nodesGeoList(:,:)
        real(IDS_real)    ::  scalarCells(2)
        type(ids_edge_profiles)     ::  edge_profiles
        type (ids_edge_sources)     ::  edge_sources
        type (ids_edge_transport)   ::  edge_transport
        type(ids_edge_profiles_time_slice), pointer      ::  ggd
        type(ids_generic_grid_dynamic), pointer          ::  grid
        type(ids_generic_grid_dynamic_space), pointer    ::  space
        type(ids_generic_grid_scalar), pointer    :: idsField

        !> ===  SET UP IDS ===
        write(0,*) "IDS parameters"
        write(0,*) "shot: ", shot
        write(0,*) "run:", run
        write(0,*) "username: ", username
        write(0,*) "device: ", device
        write(0,*) "version: ", version

        !> We already went through the check if size(crx)==size(cry)
        !> so we can safely assume that num_nodes == size(crx)
        num_nodes_all   = size(crx) !> Number of all available coordinates

        write(0,*) "Setting data for edge_profiles IDS"

        !> Preparing database for writing
        !> Through practice it was disclosed that there are some mandatory
        !> steps to be done in order to assure for data to be successfully
        !> written to IDS. Without going through those steps errors and failed
        !> process of writing to IDS are to be expected.
        !> This can be done using exampleSetIDSFundamentals routine
        homogeneous_time = 1
        time = 0.0_IDS_real
        call exampleSetIDSFundamentals( edge_profiles, homogeneous_time, time)

        !> Allocate ggd slice
        allocate(edge_profiles%ggd(1))
        !> Set pointers to top IDS nodes/structures
        ggd => edge_profiles%ggd(1)
        grid => ggd%grid

        !> Set IDS grid description
        grid_description = "This IDS was written by b2_ual_write_gsl"
        grid%identifier%description = grid_description

        !> Set IDS comment under ids_properties substructure
        allocate( edge_profiles%ids_properties%comment(1) )
        edge_profiles%ids_properties%comment(1) = &
            &   "This IDS was created by b2_ual_write_gsl using Grid Service    &
            &   Library (GSL)."
        !> Set IDS code name
        allocate( edge_profiles%code%name(1) )
        edge_profiles%code%name(1) = "b2_ual_write_gsl"

        !> === SET UP GRID===

        !> The 2D structured grid is in our case composed out of one 2D
        !> structured space
        !> Set definition of the coordinate system of the space
        coordtype(:) = (/ IDS_COORDTYPE_R, IDS_COORDTYPE_Z /)

        !> --- Set up grid space objects for Class 1 objects - points ---

        !> We are in 2D dimension space and the nodes are defined by coordinates
        !> in a way
        !> n1=[r1, z1], n2=[r2,z2], ..., nn=[rn, zn].

        !> Set geometry (R,Z) coordinate of each node object
        num_nodes = num_nodes_all
        allocate(nodesGeoList( num_nodes_all, 2))

        !> To get nodes coordinates to 2D array in correct order
        icount = 0
        do k = 0, 3
            do j = -1, ny
                do i = -1, nx
                    icount = icount + 1
                    nodesGeoList(icount, 1) = crx(i,j,k)
                    nodesGeoList(icount, 2) = cry(i,j,k)
                enddo
            enddo
        enddo

        !> --- (TODO) Set up connectivity array of grid space objects for   ---
        !> --- Class 2 objects - edges                                      ---
        !> TODO! Currently only placeholder edges are given

        !> Set (placeholder) list of indices for nodes defining each edge object
        allocate(edgesNodesList(1,2))
        edgesNodesList(1,1) = 0
        edgesNodesList(1,2) = 0

        !> --- Set up connectivity array of grid space objects for      ---
        !> --- Class 3 objects - 2D cells                               ---

        !> Set list of indices for nodes defining each cell object
        !>
        !> We are are in 2D dimension and each cell is defined by 4 nodes
        !> Note: cell nodes must be ordered cyclically (here anti-clockwise)
        numCellsX = nx + 2
        numCellsY = ny + 2
        num_cells = numCellsX * numCellsY
        allocate(cellsNodesList( num_cells, 4))

        cellId = 1
        do j = 1,numCellsY
            do i = 1, numCellsX
                cellsNodesList(cellId, 1) = cellId+0*numCellsX*numCellsY
                cellsNodesList(cellId, 2) = cellId+1*numCellsX*numCellsY
                cellsNodesList(cellId, 3) = cellId+3*numCellsX*numCellsY
                cellsNodesList(cellId, 4) = cellId+2*numCellsX*numCellsY
                cellId = cellId + 1
            enddo
        enddo

        !> --- Set the grid space objects and grid subsets ---
        !! For that we use GSL routine gridSetup2dSpace
        call gridSetup2dSpace(  grid, coordtype,                &
                            &   geo_0dObj   = nodesGeoList,     &
                            &   conn_1dObj  = edgesNodesList,   &
                            &   conn_2dObj  = cellsNodesList,   &
                            &   createGridSubsets = .true. )

        !> --- (Optional) Set grid subsets custom description   ---
        grid%grid_subset(1)%identifier%description = "All nodes in the domain."
        grid%grid_subset(3)%identifier%description = "All cells in the domain."

        !> (optional) Change name of the second grid subset
        grid%grid_subset(2)%identifier%name = "Edges - empty."

        !> === SET VALUES (ne, te, ti) for "Cells" grid subset ===
        !> For that we use GSL routine gridStructWriteData1d
        gridSubset_index = 3

        !> --- Set ne (electron density) ---
        allocate(ggd%electrons%density(1))
        idsField => ggd%electrons%density(1)
        call gridStructWriteData1d( grid, idsField, gridSubset_index, ne)

        !> --- Set te (electron temperature) ---
        allocate(ggd%electrons%temperature(1))
        idsField => ggd%electrons%temperature(1)
        !> convert to eV (1 J = 6.242e18 eV)
        call gridStructWriteData1d( grid, idsField, gridSubset_index,   &
            &   te*6.242e18)

        !> --- Set ti (ion temperature) ---
        allocate(ggd%ion(1))
        allocate(ggd%ion(1)%temperature(1))
        idsField => ggd%ion(1)%temperature(1)
        !> convert to eV (1 J = 6.242e18 eV)
        call gridStructWriteData1d( grid, idsField, gridSubset_index,   &
            &   ti*6.242e18)

        !> Set data to edge_profiles IDS
        write(0,*) "Writing to edge_profiles IDS"

        !> Create and modify new shot/run
        call imas_create_env(treename, shot, run, 0, 0, idx, username, &
            device, version)

        !> Or open and modify existing shot/run (might work much faster than
        !> imas_create_env)
        ! call imas_open_env('treename', shot, run, idx, username, device, version)

        !> Put data to IDS
        call ids_put_slice(idx,"edge_profiles",edge_profiles)
        call ids_put(idx,"edge_profiles",edge_profiles)

        !> Close IDS
        call ids_deallocate(edge_profiles)
        call imas_close(idx)

        write(0,*) "IDS write finished"

    end subroutine write_ids_edge_profiles

    !> subroutine for reading edge_profiles IDS
    subroutine read_ids(treename, shot, run, idx, username, &
                                        & device, version)
        !> Example for reading IDS database in Fortran90
        integer                 ::  shot, run, idx
        integer                 ::  gridSubset_index
        character(len=24)       ::  treename, username, device, version
        type(ids_edge_profiles) ::  edge_profiles
        character(len=255)      ::  imas_connect_url

        gridSubset_index = 3

        !> Open input datafile from local database
        write (0,*) "Started reading input IDS", idx, shot, run

        call imas_open_env('treename', shot, run, idx, username, device, version)
        call ids_get(idx, 'edge_profiles', edge_profiles)

        write(0,*) 'homogeneous_time = ', &
            & edge_profiles%ids_properties%homogeneous_time
        write(0,*) "Grid subset 3 name = ", &
            & edge_profiles%ggd(1)%grid%grid_subset(gridSubset_index)% &
            & identifier%name
        write(0,*) "Grid subset 3 index = ", &
            & edge_profiles%ggd(1)%grid%grid_subset(gridSubset_index)% &
            & identifier%index
        ! write(0,*) "Time = ", edge_profiles%time(1)
        call ids_deallocate(edge_profiles)
        call imas_close(idx)
        write (0,*) "Finished reading input IDS"

    end subroutine read_ids

end program b2_ual_write_gsl

!!!Local Variables:
!!! mode: f90
!!! End:
