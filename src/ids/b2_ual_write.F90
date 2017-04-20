!!------------------------------------------------------------------------------
!!.specification

program b2_ual_write
    use b2mod_types , B2R8 => R8
    use b2mod_version
    use b2mod_geo
    use b2mod_plasma
    use b2mod_rates
    use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_anomalous_transport
    use b2mod_work
    use b2mod_time
    use b2mod_ppout
    use b2mod_tallies
    use b2mod_indirect
    use b2mod_neutrals_namelist
    use b2mod_neutr_src_scaling
    use b2mod_b2cmfs
    use b2mod_constants
    use b2mod_grid_mapping
    use b2mod_ual_io_grid
    use b2mod_ual_io_data
    use ggd
    use b2mod_ual
    use b2mod_ual_io

    implicit none

    !!--------------------------------------------------------------------------
    !!.documentation

    !! 1. purpose
    !!
    !!      b2_ual_write.f90 script is used to generate b2_ual_write.exe (main
    !!      program), which is a post-processor for b2. It reads the plasma 
    !!      state writes it to IDS database
    !!
    !!
    !! 2. specification
    !!
    !!      Main program.
    !!
    !!
    !! 3. description (see also routine b2cdca)
    !!
    !!      The complete program performs post-processing of the
    !!      result of a b2 calculation.
    !!      This program unit opens and closes the input/output units, and
    !!      may perform some other system-dependent operations. 
    !!
    !!      The input units are:
    !!      ninp(0): formatted; provides output control parameters.
    !!      ninp(1): un*formatted; provides the geometry.
    !!      ninp(2): un*formatted; provides the run parameters.
    !!      ninp(3): un*formatted; provides the plasma state.
    !!      ninp(4): unformatted; provides the detailed plasma state.
    !!      ninp(5): formatted; provides the run switches.
    !!      ninp(6): un*formatted; provides the atomic data.
    !!
    !!      The output units are:
    !!      nout(0): formatted; provides printed output.
    !!
    !!      (See routine b2cdca for the meaning of 'un*formatted'.)
    !!
    !!
    !! 4. parameters (see also routine b2cdcv)
    !!
    !!      None.
    !!
    !!
    !! 5. error indicators
    !!
    !!      In case an error condition is detected, a call is made to the
    !!      routine xerrab. This causes an error message to be printed,
    !!      after which the program halts.
    !!
    !!--------------------------------------------------------------------------
    !!.declarations

    !!..common blocks

    !!..local variables
    integer ninp(0:6), nout(0:2), nx, ny, ns, idum(0:9), io
    integer                     ::  idx
    integer                     ::  ix, ixx
    integer                     ::  shot, run
    integer                     ::  test_int1, test_int2, test_int3, test_int4
    character(len=24)           ::  treename, username, machine, version
    !! character(len=255)          ::  imas_connect_url
    type(ids_edge_profiles)     ::  edge_profiles

    !!  ..procedures
    external prgini, prgend, xerset, xertst, cfopen, cfruin
    !!  ..initialize input/output units
    data ninp/50, 51, 52, 53, 54, 55, 56/
    data nout/60, 61, 62/
    
    !!--------------------------------------------------------------------------
    !!.documentation-internal
    !!
    !!      The following common blocks have their outermost declaration in
    !!      this routine; they need not be preserved between calls.
    !!
    !!      Description of some local variables:
    !!
    !!      ninp - (0:6) integer array.
    !!      ninp specifies the input unit numbers.
    !!
    !!      nout - (0:2) integer array.
    !!      nout specifies the output unit numbers.
    !!
    !!      nx, ny - integer.
    !!      nx and ny specify the number of interior cells along the first
    !!      and the second coordinate, respectively. The total number of
    !!      cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
    !!      It will hold that 0.le.nx and 0.le.ny.
    !!
    !!      ns - integer.
    !!      ns specifies the number of atomic species in the calculation.
    !!      The species are indexed by (0:ns-1).
    !!      It will hold that 1.le.ns.
    !!          
    !!--------------------------------------------------------------------------
    !!.computation

    ! call read_b2fgmtry_b2fstate()

    treename    = "ids"
    shot        = 8148
    run         = 3
    username    = "penkod"
    machine     = "solps-iter"
    version     = "3"

    call write_ids_edge_profiles(treename, shot, run, idx, username, &
                                        & machine, version)

    call read_ids(treename, shot, run, idx, username, &
                                        & machine, version)

contains 

    subroutine read_b2fgmtry_b2fstate()

        integer                 ::  i, j, k
        character(len=120)      ::  lblgm

        !!..program start_up calls
        call prgini('b2_ual_write')

        !!..open files
        !!...input files
        ! call cfopen (ninp(0),'b2_ual_write.dat','old','formatted')
        call cfopen(ninp(1),'b2fgmtry','old','un*formatted')
        ! call cfopen (ninp(2),'b2fparam','old','un*formatted')
        call cfopen(ninp(3),'b2fstate','old','un*formatted')
        ! call cfopen (ninp(4),'b2fplasma','old','unformatted')
        ! call cfopen (ninp(5),'b2mn.dat','old','formatted')
        ! call cfopen (ninp(6),'b2frates','old','formatted')
        !!...output files
        call cfopen(nout(0),'b2_ual_write.prt','new','formatted')
        call cfopen(nout(1),'b2_ual_writ2.prt','new','formatted')

        !!..mark output unit for error messages
        call xerset(0)
        !!..obtain version numbers
        call cfverr(ninp(1), b2fgmtry_version)
        ! call cfverr (ninp(2), b2fparam_version)
        call cfverr(ninp(3), b2fstate_version)
        ! call cfverr (ninp(6), b2frates_version)
        !!..obtain nx, ny, ns
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
        write (0,*) "size(crx): ", size(crx)

        write (nout(1),*) crx
        
        corner: do k = 0, 3, 1
            xi: do j = -1, ny, 1
                yi: do i = -1, nx, 1
                    write (nout(0),*) i, j, k, crx(i, j, k)
                end do yi
            end do xi
        end do corner

    ! TODO read ne, ni, te and ti 

    end subroutine read_b2fgmtry_b2fstate

    subroutine write_ids_edge_profiles(treename, shot, run, idx, username, &
                                        & machine, version)
        integer                 ::  shot, run, idx
        character(len=24)       ::  treename, username, machine, version
        type(ids_edge_profiles) ::  edge_profiles
        real(B2R8)              :: time

        !!.. Create IDS 
        write(0,*) "Writing IDS"
        write(0,*) "shot: ", shot
        write(0,*) "run:", run
        write(0,*) "username: ", username
        write(0,*) "machine: ", machine
        write(0,*) "version: ", version

        ! call imas_create(treename, shot, run, 0, 0, idx)
        call imas_create_env(treename, shot, run, 0, 0, idx, username, machine, &
            & version)

        allocate(edge_profiles%profiles_1d(1))
        allocate(edge_profiles%ggd(1))
        ! allocate(edge_profiles%putNonTimed())
        allocate(edge_profiles%ggd(1)%grid%space(1))
        allocate(edge_profiles%ggd(1)%grid%space(1)%objects_per_dimension(3))
        edge_profiles%ids_properties%homogeneous_time = 1

        allocate(edge_profiles%time(1))
        edge_profiles%time(1)=time



        call ids_put(idx,"edge_profiles",edge_profiles)
        call ids_deallocate(edge_profiles)
        call imas_close(idx)

    end subroutine write_ids_edge_profiles

    !!!!!!!!!!!!!!!Opening IDS to check if writing works correctly !!!!!!!!!!!!
    subroutine read_ids(treename, shot, run, idx, username, &
                                        & machine, version)
        ! Example for readin IDS database in Fortran90
        integer                 ::  shot, run, idx
        character(len=24)       ::  treename, username, machine, version
        type(ids_edge_profiles) ::  edge_profiles
        character(len=255) :: imas_connect_url

        ! Open input datafile from local database
        write (0,*) "Started reading input IDS", idx, shot, run

        call imas_open_env('treename', shot, run, idx, username, machine, version)
        call ids_get(idx, 'edge_profiles', edge_profiles)

        write(0,*) 'homogeneous_time = ', edge_profiles%ids_properties%homogeneous_time
        call imas_close(idx)
        write (0,*) "Finished reading input IDS"

    end subroutine read_ids

end program b2_ual_write

