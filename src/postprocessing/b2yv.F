*-----------------------------------------------------------------------
*.specification

      program b2yv
      use b2mod_types
      use b2mod_constants
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_ma28
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_elements
#ifdef WINDOWS
! AVDef is the AViz module file
      use AVDef    ! IGNORE
      use AVVIEWER ! IGNORE
      use DFLIB    ! IGNORE
#endif
      implicit none
#ifdef USE_MPI
      include 'mpif.h'
#endif

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YV reads in an unformatted b2fplasma file and visualizes
*     b2plasma's arrays.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer :: icase, ix, iy, isymm, is
      integer :: lbc, ubc, lbr, ubr   
      real(4), allocatable :: M(:,:)
      character(1) :: key
      character lblgm*120
#ifdef WINDOWS
      integer :: hv, lbnd(2) = 0
      integer status
      character(AV_MAX_LABEL_LEN) :: axis1 = 'X', axis2 = 'Y', axis3
#else
      character*1 :: axis1 = 'X', axis2 = 'Y'
      character*40 :: axis3
#endif

      integer ninp(0:5), nout(0:7), nx, ny, ns, ns0, idum(0:9), io
#ifdef USE_MPI
      integer rank_mpi,size_mpi,ierr_mpi
#endif
      character lblmn*120, title*64
      character hlp_frm*22
      logical file_open
*   ..procedures
      external prgini, prgend, xerset, xertst, xerrab, cfopen, cfruin
*   ..initialize input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/

*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
#ifdef USE_MPI
      call mpi_init(ierr_mpi)
      call mpi_comm_rank(MPI_COMM_WORLD,rank_mpi,ierr_mpi)
      call mpi_comm_size(MPI_COMM_WORLD,size_mpi,ierr_mpi)
      if(rank_mpi .ne. 0) then
        call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
        call mpi_finalize(ierr_mpi)
        stop 'b2yv_mpi_child'
      endif
#endif
      call prgini ('b2yv')
*   ..open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
!      call cfopen (nout(0),'b2yv.prt','new','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version number						 !srv 18.09.03
      call cfverr (ninp(1), b2fgmtry_version)
      call cfverr (ninp(2), b2fparam_version)
      call cfverr (ninp(3), b2fstati_version)			 !srv 18.09.03
*   ..obtain nx, ny from geometry file
      call cfruin (ninp(1), 2, idum, 'nx,ny')
      nx = idum(0)
      ny = idum(1)
      call xertst (0.le.nx.and.0.le.ny,
     &  'faulty input nx, ny from geometry file')
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain nx, ny, ns0 from plasma state file
      call cfruin (ninp(3), 3, idum, 'nx,ny,ns')
      ns0 = idum(2)
      call xertst (idum(0).eq.nx.and.idum(1).eq.ny,
     &  'faulty input nx, ny from plasma state file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstati_version,
     & ns, zamin, zamax, zn, am, .true.)
      lbc = -1
      lbr = -1
      ubc = ny
      ubr = nx
      allocate (M(lbc:ubc, lbr:ubr))
*   ..call principal driver routine
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_residuals(nx,ny,ns)
      call alloc_b2mod_sources(nx,ny,ns)
      call alloc_b2mod_transport(nx,ny,ns)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)                       !srv 25.11.02
*   ..read geometry and metric
      call b2rugm (ninp(1), nx, ny, crx, cry, fpsi, ffbz,
     &  bb, vol, hx, hy, qz, qc, qcb, gs, pbs, wbbl, wbbr, wbbv, wbbc,
     &  cell_width, cell_height, gmap)
*   ..read plasma state
      call cfverr(ninp(4), b2fplasma_version)
      call read_b2mod_geo(nx, ny, ninp(4))
      call read_b2mod_plasma(nx, ny, ns, ninp(4))
      call read_b2mod_residuals(ninp(4))
      call read_b2mod_sources(ninp(4))
      call read_b2mod_transport(ninp(4))
      call alloc_b2mod_elements(ns,zn,am)
   10 continue
      write(hlp_frm,'(a,i3,a)') '(a,',ns,'(a2,i1,1x))'
      write(*,hlp_frm) 'To look at density of',
     &      (adjustr(is_codes(is)),nint(zamin(is)),is=0,ns-1)
      write(hlp_frm,'(a,i3,a)') '(4x,a,',ns,'(1x,I2,1x)/)'
      write(*,hlp_frm)
     &        ' - input a number', (is,is=0,ns-1)
      write(hlp_frm,'(a,i3,a)') '(a,',ns,'(a2,i1,1x))'
      write(*,hlp_frm) 'To look at parallel velocity of',
     &      (adjustr(is_codes(is)),nint(zamin(is)),is=0,ns-1)
      write(hlp_frm,'(a,i3,a)') '(14x,a,',ns,'(1x,I2,1x)/)'
      write(*,hlp_frm)
     &        ' - input a number', (is+ns,is=0,ns-1) 

      write(*,"(a,a,(1x,I2,1x))") 'To look at electron tempetature',
     &        '       - input a number', ns+ns

      write(*,"(a,a,(1x,I2,1x))") 'To look at ion tempetature     ',
     &        '       - input a number', ns+ns+1

      write(*,"(a,a,(1x,I2,1x))") 'To look at electric potential',
     &        '       - input a number', ns+ns+2

      write(*,"(a,a,(1x,I2,1x))") 'To look at radial electric field',
     &        '    - input a number', ns+ns+3

      read (*,*) icase 

      if(0.le.icase .and. icase.le.ns-1 )then
	write(title,"('Density of ',a2,i1)") 
     &     adjustr(is_codes(icase)),nint(zamin(icase))
	write(axis3,"('Density of ',a2,i1,', [1/m^3]')") 
     &     adjustr(is_codes(icase)),nint(zamin(icase))
        do iy = -1,ny
          do ix = -1,nx
	    M(iy,ix) = na(ix,iy,icase)
	  enddo
	enddo
      elseif(ns.le.icase .and. icase.le.ns-1+ns )then
	write(title,"('Parallel velocity of ',a2,i1)") 
     &     adjustr(is_codes(icase-ns)),nint(zamin(icase-ns))
	write(axis3,"('Parallel velocity of ',a2,i1,', [m/c]')") 
     &     adjustr(is_codes(icase-ns)),nint(zamin(icase-ns))
        do iy = -1,ny
          do ix = -1,nx
	    M(iy,ix) = ua(ix,iy,icase-ns)
	  enddo
	enddo
      elseif(icase .eq. ns+ns )then
	write(title,"('Temperature of electrons')") 
	write(axis3,"('Temperature of electrons',', [eV]')") 
        do iy = -1,ny
          do ix = -1,nx
	    M(iy,ix) = te(ix,iy)/ev
	  enddo
	enddo
      elseif(icase .eq. ns+ns+1 )then
	write(title,"('Temperature of ions')") 
	write(axis3,"('Temperature of ions',', [eV]')") 
        do iy = -1,ny
          do ix = -1,nx
	    M(iy,ix) = ti(ix,iy)/ev
	  enddo
	enddo
      elseif(icase .eq. ns+ns+2 )then
	write(title,"('Electric potential')") 
        do iy = -1,ny
          do ix = -1,nx
	    M(iy,ix) = po(ix,iy)
	  enddo
	enddo
      elseif(icase .eq. ns+ns+3 )then
	write(title,"('Radial electric field')") 
        M=0.0d0
	do iy=0,ny
	  do ix=-1,nx
	    M(iy,ix)=-(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))
     *                /(hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
     *                *2.0e0_R8
	  enddo
	enddo
      else
        call xerrab ('Wrong icase value')
      endif

#ifdef WINDOWS
      call faglStartWatch(M, status)
      lbnd(1:size(shape(M))) = lbound(M)
      call faglLBound(M, lbnd, status)

      call favStartViewer(hv,status)
      call favSetArray(hv,M,status)
      call favSetArrayName(hv,trim(title), status)
      call favSetUseAxisLabel (hv, x_axis, 1 , status)
      call favSetUseAxisLabel (hv, y_axis, 1 , status)
      call favSetUseAxisLabel (hv, z_axis, 1 , status)
      call favSetAxisLabel (hv, x_axis, axis1, status)
      call favSetAxisLabel (hv, y_axis, axis2, status)
      call favSetAxisLabel (hv, z_axis, trim(axis3), status)

!	call faglShow(M, status)
      call favShowWindow(hv,AV_TRUE, status)


!      call faglName(M, trim(title), status)

      write (*,*) ' Press 0 - 9 to continue',
     &            ' or press any other key to finish'
      key = GETCHARQQ()
      call favEndViewer(hv,status)
!       call faglClose(M, status)
      call faglEndWatch(M, status)  
#endif

      select case (key)
	case ('0':'9')
	  goto 10
	case default
      end select

*   ..close files
      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo

      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_sources
      call dealloc_b2mod_elements
      call dealloc_b2mod_indirect
      call dealloc_b2mod_residuals
      call dealloc_b2mod_transport
      deallocate (M)
* ..end of computation
      call prgend ()
#ifdef USE_MPI
      call mpi_barrier(MPI_COMM_WORLD,ierr_mpi)
      call mpi_finalize(ierr_mpi)
#endif
      stop 'b2yv'

*-----------------------------------------------------------------------
*.end b2yv

      end
