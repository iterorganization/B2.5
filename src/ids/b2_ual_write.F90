!-----------------------------------------------------------------------
!     .specification

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

  ! UAL Access
#ifdef IMAS
!  use ids_schemas  ! IGNORE
!  use ids_routines ! IGNORE
#else
#ifdef ITM
  use euITM_schemas  ! IGNORE
  use euITM_routines ! IGNORE
#endif
#endif

  ! B2/IDS-CPO Mapping
  use b2mod_grid_mapping
  use b2mod_ual_io_grid
  use b2mod_ual_io_data
  use ggd
  use b2mod_ual
  use b2mod_ual_io


  implicit none

  !-----------------------------------------------------------------------
  !     .documentation
  !
  !     1. purpose
  !
  !     b2_ual_write (main program) is a postprocessor for b2. It reads the
  !     plasma state and writes it as a CPO the ITM UAL.
  !
  !
  !     2. specification
  !
  !     Main program.
  !
  !
  !     3. description (see also routine b2cdca)
  !
  !     The complete program performs postprocessing of the
  !     result of a b2 calculation.
  !     This program unit opens and closes the input/output units, and
  !     may perform some other system-dependent operations. 
  !
  !     The input units are:
  !     ninp(0): formatted; provides output control parameters.
  !     ninp(1): un*formatted; provides the geometry.
  !     ninp(2): un*formatted; provides the run parameters.
  !     ninp(3): un*formatted; provides the plasma state.
  !     ninp(4): unformatted; provides the detailed plasma state.
  !     ninp(5): formatted; provides the run switches.
  !     ninp(6): un*formatted; provides the atomic data.
  !
  !     The output units are:
  !     nout(0): formatted; provides printed output.
  !
  !     (See routine b2cdca for the meaning of 'un*formatted'.)
  !
  !
  !     5. parameters (see also routine b2cdcv)
  !
  !     None.
  !
  !
  !     6. error indicators
  !
  !     In case an error condition is detected, a call is made to the
  !     routine xerrab. This causes an error message to be printed,
  !     after which the program halts.
  !
  !-----------------------------------------------------------------------
  !     .declarations

  !     ..common blocks

  !     ..local variables
  integer ninp(0:6), nout(0:2), nx, ny, ns, idum(0:9), io
  integer :: idx
  integer :: ix, ixx
#ifdef IMAS
  type (ids_edge_profiles)  :: edge_profiles
  type (ids_edge_sources)   :: edge_sources
  type (ids_edge_transport) :: edge_transport
  type (ids_core_profiles)  :: cp
#else
#ifdef ITM
  type (type_edge),pointer :: cpoedge(:) => null()
#endif
#endif
  real(B2R8) :: time
  character*120 lblgm
  logical includeGhostCells

  integer :: shot, run
  character(len=255) :: imas_connect_url

  !     ..procedures
  external prgini, prgend, xerset, xertst, cfopen, cfruin
  !..initialise input/output units
  data ninp/50,51,52,53,54,55,56/
  data nout/60,61,62/

  !-----------------------------------------------------------------------
  !     .documentation-internal
  !
  !     The following common blocks have their outermost declaration in
  !     this routine; they need not be preserved between calls.
  !
  !     Description of some local variables:
  !
  !     ninp - (0:6) integer array.
  !     ninp specifies the input unit numbers.
  !
  !     nout - (0:2) integer array.
  !     nout specifies the output unit numbers.
  !
  !     nx, ny - integer.
  !     nx and ny specify the number of interior cells along the first
  !     and the second coordinate, respectively. The total number of
  !     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
  !     It will hold that 0.le.nx and 0.le.ny.
  !
  !     ns - integer.
  !     ns specifies the number of atomic species in the calculation.
  !     The species are indexed by (0:ns-1).
  !     It will hold that 1.le.ns.
  !
  !-----------------------------------------------------------------------
  !     .computation

!  call ids_test_put
  !     ..program start-up calls
  call prgini ('b2_ual_write')
  !     ..open files
  !call cfopen (ninp(0),'b2_ual_write.dat','old','formatted')
  call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
!!  call cfopen (ninp(2),'b2fparam','old','un*formatted')
  call cfopen (ninp(3),'b2fstate','old','un*formatted')
!!  call cfopen (ninp(4),'b2fplasma','old','unformatted')
!!  call cfopen (ninp(5),'b2mn.dat','old','formatted')
!!  call cfopen (ninp(6),'b2frates','old','formatted')
  call cfopen (nout(0),'b2_ual_write.prt','new','formatted')
  !     ..mark output unit for error messages
  call xerset (0)
  !     ..obtain version numbers						 !xpb
  call cfverr (ninp(1), b2fgmtry_version) !xpb
!!  call cfverr (ninp(2), b2fparam_version) !xpb
  call cfverr (ninp(3), b2fstate_version) !xpb
!!  call cfverr (ninp(6), b2frates_version) !xpb
  !     ..obtain nx, ny, ns
  call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
  nx = idum(0)
  ny = idum(1)
  ns = idum(2)
  call xertst (0.le.nx.and.0.le.ny.and.1.le.ns, &
       & 'faulty input nx, ny, ns from plasma state file')
  call cfruin (ninp(1), 2, idum(0), 'nx,ny')
  call xertst (idum(0).eq.nx.and.idum(1).eq.ny, &
       & 'faulty input nx, ny from geometry file')

!!! moved from read_b2fplasma
  call alloc_b2mod_geo(nx,ny)
  call alloc_b2mod_plasma(nx,ny,ns)
  call alloc_b2mod_indirect(nx,ny,nncutmax) !xpb
  call alloc_b2mod_sources(nx,ny,ns)
!!! end of moved calls

!!  call read_b2fplasma(ninp, nx, ny, ns)
  call read_additional(ninp, nout, nx, ny, ns)

  ! Fill IDS or CPO
#ifdef IMAS
  call write_ids(edge_profiles,edge_sources,edge_transport)
#else
# ifdef ITM
   allocate(edge_profiles(1))
   call write_cpo(edge_profiles(1))
# endif
#endif


  idx = -666
  shot = 8148
  run = 1

#ifdef IMAS
  call getenv ("imas_connect_url", imas_connect_url)
  if (0.lt.len_trim(imas_connect_url)) then
    write(*,*) 'imas_connect_url: ', trim(imas_connect_url)
    call imas_connect(imas_connect_url) ! Use to connect to a remote IMAS database / account
  endif
#endif
  ! Write IDS or CPO to UAL
  call open_ual(idx, shot, run, time=time, doCreate=.true., useHdf5=.false., nmlFile='b2_ual_write.dat')
  write(*,*) 'After open_ual: idx = ', idx

#ifdef IMAS
  allocate(edge_profiles%ids_properties%comment(1))
  edge_profiles%ids_properties%homogeneous_time = 1
  edge_profiles%ids_properties%comment(1) = 'Test put (edge)'
  allocate(edge_profiles%time(1))
  edge_profiles%time(1)=tim
  write(0,*) 'nx, ny : ', nx, ny
  allocate(edge_profiles%profiles_1d(1))
  allocate(edge_profiles%profiles_1d(1)%electrons%temperature(-1:ny))
  edge_profiles%profiles_1d(1)%electrons%temperature(-1:ny)=te(nx/2,-1:ny)/ev
  allocate(edge_profiles%ggd(1))
  call ids_put(idx,"edge_profiles",edge_profiles)
  call ids_deallocate(edge_profiles)
#else
# ifdef ITM
   call euitm_put(idx,"edge",cpoedge)
# endif
#endif
  call close_ual(idx)

#ifdef IMAS
! trying to re-open and read
  if (0.lt.len_trim(imas_connect_url)) then
    call imas_connect(imas_connect_url) ! Use to connect to a remote IMAS database / account
  endif
#endif
  call open_ual(idx,shot,run, time=time, doCreate=.true., useHdf5=.false., nmlFile='b2_ual_write.dat')
#ifdef IMAS
  call ids_get(idx,"edge_profiles",edge_profiles)
  write(*,*) "edge_profiles"
  write(*,*) "IDS_Properties homogeneous : ", edge_profiles%IDS_Properties%homogeneous_time
  write(*,'("IDS_Properties comment(1) : ",A)') edge_profiles%IDS_Properties%comment(1)
  write(*,'("Time : ",g14.7)') edge_profiles%time(1)
  call ids_deallocate(edge_profiles)
#else
# ifdef ITM
  call euitm_get(idx,"edge",cpoedge)
# endif
#endif
  call close_ual(idx)


  !     ..close files			 !
  do io=1,size(ninp)-1
     close(ninp(io))
  enddo
  do io=0,size(nout)-1
     close(nout(io))
  enddo

  call dealloc_b2mod_geo
  call dealloc_b2mod_plasma
  call dealloc_b2mod_indirect
!!$  call dealloc_b2mod_regions
  call dealloc_b2mod_sources

  !     ..end of computation
  call prgend()
  stop 'b2_ual_write'

  !-----------------------------------------------------------------------
  !     .end b2_ual_write

contains

  subroutine read_b2fplasma(ninp, nx, ny, ns)
    integer ninp(0:6), nx, ny, ns

    ! allocate module data structures and read b2fplasma file (taken from b2md.F)

    call alloc_b2mod_rates(nx,ny,ns)
    call alloc_b2mod_residuals(nx,ny,ns)
    call alloc_b2mod_transport(nx,ny,ns)
    call alloc_b2mod_anomalous_transport(nx,ny,ns)
    call alloc_b2mod_work(nx,ny,ns)
    ! ..read geometry
    call cfruch (ninp(1), 120, lblgm, 'label')
    call cfruin (ninp(1), 1, idum, 'isymm')
    isymm = idum(0)
    call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz, &
     &  bb, vol, hx, hy, qz, qc, qcb, gs, pbs, &
     &  wbbl, wbbr, wbbv, wbbc, cell_width, cell_height, gmap)
    !     ..read plasma state
    call cfverr(ninp(4), b2fplasma_version)
    call read_b2mod_geo(nx,ny,ninp(4))
    call read_b2mod_plasma(nx,ny,ns,ninp(4))
    call read_b2mod_residuals(ninp(4))
    call read_b2mod_sources(ninp(4))
!!$    call read_b2mod_transport(nx,ny,ns,ninp(4))
    call read_b2mod_anomalous_transport(ninp(4))
    call read_b2mod_neutr_src_scaling(ninp(4),ns,nstrai)
    call alloc_b2mod_ppout(nx,ny,ns)


  end subroutine read_b2fplasma
  
  subroutine read_additional(ninp, nout, nx, ny, ns)

    ! read and compute additional data (taken from b2mddr.F)

    use b2mod_version
    use b2mod_geo
    use b2mod_plasma
    use b2mod_rates
    use b2mod_residuals
    use b2mod_sources
    use b2mod_transport
    use b2mod_ppout
    use b2mod_indirect
    use b2mod_constants
    use b2mod_geo_corner
    use b2mod_work
    use b2mod_time
    use b2mod_b2cmfs
    use b2mod_b2cmrc
    use b2mod_b2cmgs
    use b2mod_b2cmpa
    use b2mod_b2cmpb
    use b2mod_b2cmpt
    use b2mod_b2cmwg
    implicit none

    !   ..input arguments (unchanged on exit)
    integer ninp(0:6), nout(0:2), nx, ny, ns
    !   ..output arguments (unspecified on entry)
    !     (none)
    !   ..common blocks

    !   ..local variables
    integer k, is, idum(0:9), nscx, iscx(0:nscxmax-1), ismain
    real (kind=B2R8) :: BoRiS, t0, t1
    real (kind=B2R8) ::&
         & rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),&
         & rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1), &
         & rbr(-1:nx,-1:ny,0:ns-1), rcx(-1:nx,-1:ny,0:ns-1)
    character lblgm*120, lblcp*120, lblmn*120, lblrc*120
    character cnamip*80, cvalip*80
    !   ..initialisation
    save BoRiS, ismain
    data BoRiS/0.0_B2R8/, ismain/1/
    !   ..procedures
    intrinsic abs
!!    logical ltst
!!    ltst(t0,t1) = abs(t0-t1).le.1.0e-6_B2R8*(abs(t0)+abs(t1))
    external subini, subend, xertst, xerrab, cfruch, cfruin
    external b2rups, b2rugm, b2rucp, ipsetc, ipgetr
    external b2rurc, b2spel, b2sqel, b2sqcx, b2spcx, b2xvps
    external b2xxid, b2rflb, b2rfcp
    !   ..namelist
    character :: exp*128, comment*128
    integer :: shot, overwrite_shotnumber
    real (kind=B2R8) :: time
    logical :: timedep, snapshot, tallies, movies
    namelist /b2md_namelist/ exp, shot, time, comment, timedep, snapshot, tallies, movies, overwrite_shotnumber

    !   ..subprogram start-up calls
    call subini ('b2mddr')
    !   ..test ninp, nout
    call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.&
         & 1.le.ninp(3).and.1.le.ninp(5).and.1.le.ninp(6).and.&
         & 1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2),&
         & 'faulty argument ninp, nout')
    !   ..test dimensions
    call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
    call xertst (1.le.ns, 'faulty argument ns')
    call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

    ! ..computation
    !   ..obtain lblmn
    call b2xxid ('b2mn', lblmn(1:60))
!!    call b2rflb (ninp(5), lblmn(61:120))
    !   ..read revised parameters in physics common
    !     (ninp(0) need supply only parameters that are to be changed from
    !     the values read in as default physics common.)
!!    call b2rfcp (ninp(5), ns)
    !   ..read and echo code internal parameters
    write (nout(0),'(/2x,a)') 'non-default internal parameters:'
1   continue
!!    read (ninp(5),*,end=2,err=91) cnamip, cvalip
!!    if (cnamip(1:1).ne.'*') then
!!      write (nout(0),'(4x,a,2x,a)') cnamip, cvalip
!!      call ipsetc (cnamip, cvalip)
!!    endif
!!    goto 1
2   continue
!!    write (nout(0),'(2x,a)') '(end of list of internal parameters)'
!!    call ipgetr ('b2news_BoRiS', BoRiS)
    !   ..read geometry
    call cfruch (ninp(1), 120, lblgm, 'label')
    call cfruin (ninp(1), 1, idum, 'isymm')
    isymm = idum(0)
    call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz, &
         & bb, vol, hx, hy, qz, qc, qcb, gs, pbs, &
         & wbbl, wbbr, wbbv, wbbc, cell_width, cell_height, gmap)
    !   ..read physics common
!!    call cfruch (ninp(2), 120, lblcp, 'label')
!!    call b2rucp (ninp(2), b2fparam_version, nx, ns)
    !   ..read plasma state file label
    call cfruch (ninp(3), 120, lblmn, 'label')
    !   ..read species data
    call b2ruzd (ninp(3), b2fstate_version, ns, zamin, zamax, zn, am, .true.)
    !   ..read plasma state
    call b2rups (ninp(3), nx, ny, ns, ne, na, ua, uadia, te, ti, po, &
         & fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy, b2fstate_version)
    !   ..test plasma state and fluxes
    call b2xvps (nx, ny, ns, &
         & ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
    !   ..read atomic rate data
!!    call cfruch (ninp(6), 120, lblrc, 'label')
!!    adpak_used = index(lblrc,'ADPAK').gt.0 
!!    call read_b2mod_ppout('rates', ninp(6), lblrc)
!!    call b2rurc (ninp(6), b2frates_version)
    !   ..test rtns, rtzmin, rtzmax, rtzn
!!    call xertst (ns.eq.rtns, 'faulty input rtns--wrong rate table?')
!!     do is = 0, ns-1
!!        call xertst (ltst(zamin(is),rtzmin(is)), 'faulty input rtzmin--wrong rate table?')
!!        call xertst (ltst(zamax(is),rtzmax(is)), 'faulty input rtzmax--wrong rate table?')
!!        call xertst (ltst(zn(is),rtzn(is)), 'faulty input rtzn--wrong rate table?')
!!     enddo
!!     !   ..compute log-log electron rate coefficients
!!     call b2spel (nx, ny, ns, ev, te, ne, rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
!!     !   ..compute electron rate coefficients
!!     call ipgeti ('b2mndr_ismain', ismain)
!!     call xertst (0.le.ismain.and.ismain.lt.ns, 'invalid main plasma species index ismain')
!!     call xertst (.not.is_neutral(ismain), 'invalid main plasma species ismain; must not be neutral')
!!     call b2sqel (nx, ny, ns, ismain, ev, te, &
!!          & rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,&
!!          & rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
!!     !   ..find nscx, iscx
!!     !     (indices for neutral hydrogen species)
!!     nscx = 0
!!     do is = 0, ns-1
!!        if (is_neutral(is).and.zn(is).eq.1.0_B2R8) then
!!           call xertst (nscx.lt.nscxmax, 'too many neutral hydrogen species')
!!           iscx(nscx) = is
!!           nscx = nscx+1
!!        endif
!!     enddo
!!     do k = nscx, nscxmax-1
!!        iscx(k) = -1
!!     enddo
!!     !    ..compute log-log charge exchange rate coefficients
!!     do k = 0, nscx-1
!!        call b2spcx (nx, ny, ns, ev, am(iscx(k)), ti, ne, rlcx(-1,-1,0,0,k))
!!     enddo
!!     !   ..compute charge exchange rate coefficients
!!     do k = 0, nscx-1
!!        call b2sqcx (nx, ny, ns, ev, am(iscx(k)), ti, rlcx(-1,-1,0,0,k), rcx, wrk0)
!!     enddo

    call subend ()
    return

    !-----------------------------------------------------------------------
    !.error conditions

91  call xerrab ('error trying to read internal parameters')

  end subroutine read_additional

end program b2_ual_write

!!!Local Variables:
!!! mode: f90
!!! End:
