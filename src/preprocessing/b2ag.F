*-----------------------------------------------------------------------
*.specification

      program b2ag
      use b2mod_types
      use b2mod_geo
      use b2mod_geo2
      use b2mod_indirect
      use b2ag_parameters
	  use b2us_geo
	  use b2us_map
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AG (main program) prepares a geometry and magnetic field for
*     B2. This program calls the driver routine b2agdr for the actual
*     computation.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2ag prepares the geometry and magnetic
*     field for a B2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2agdr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides geometric parameters.
*     ninp(1): formatted; provides mhd equilibrium data.
*
*     The output units are:
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; full specification of the geometry.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
*
*
*  4. references
*
*     (Refs. /1/ and /2/ are out of date; !! to be replaced.)
*     /1/ B. J. Braams, Computational Studies in Tokamak Equilibrium
*     and Transport, thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks, NET Report Nr 68, January 1987.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     None.
*
*
*  6. error indicators
*
*     In case an error condition is detected, a call is made to the
*     routine xerrab. This causes an error message to be printed,
*     after which the program halts.
*
*
*  7. accuracy
*
*     The computation performed in this program depends on the
*     application, and it is difficult to provide a useful general
*     statement concerning the accuracy of the calculation.
*     Numerical differentiation is employed to compute the metric
*     parameters of the grid and to compute the magnetic field. So as
*     to retain sufficient accuracy when performing this numerical
*     differentiation, cells should not be too thin. Be careful about
*     the boundary cells in particular.
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      character :: id*8, id2*20
      integer ninp(0:1), nout(0:1), nx, ny, nx1, ny1, io
      logical file_open, streql
*   ..procedures
      external prgini, prgend, cfopen, streql
      external b2agx0, b2agx0_st, b2agdr
*   ..initialise input/output units
      data ninp/50,51/
      data nout/60,61/

	  type(geometry), save :: g
	  type(mapping), save :: m 

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:1) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:1) integer array.
*     nout specifies the output unit numbers.
*
*     nx, ny - integer.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     nx1, ny1 - integer.
*     nx1 and ny1 specify the mesh dimensions for an auxiliary mesh
*     that will be used to evaluate accurately the metric coefficients.
*     This auxiliary mesh has no role outside the present program.
*     nx1 will be a multiple of nx and ny1 will be a multiple of ny.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2ag')
*   ..open files
      call cfopen (ninp(0),'b2ag.dat','old','formatted')
*$$$  call cfopen (ninp(1),'b2mhdeq','old','formatted')
*     (b2mhdeq is reserved for numerical mhd equilibrium data)
      call cfopen (nout(0),'b2ag.prt','new','formatted')
      call cfopen (nout(1),'b2fgmtry','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..read format (self addition)
      read (ninp(0),'(a8)') id
      read (ninp(0),*) id2
      if ( streql(id,'*format') .AND. streql(id2,'unstructured')) then
*       run the unstructured version
	  write (*,*) 'Unstructured b2ag used'
*    ..obtain nCv, nVx, nFc
	  call b2ag_read_parameters(ninp,nout,m)
*    ..call principal driver routine
	  call alloc_mapping(m)
	  call alloc_geometry(g,m%nCv,m%nFc,m%nVx,m%nFs,m%nFt,m%nVmxCv)
	  call b2agdr (ninp, nout,g,m)
*    ..close graphics files
*      (none)
*    ..close other files
	  do io=0,size(ninp)-1
		inquire(unit=ninp(io),opened=file_open)
		if (file_open) close(ninp(io))
	  enddo
	  do io=0,size(nout)-1
		inquire(unit=nout(io),opened=file_open)
		if (file_open) close(nout(io))
	  enddo
*    ..end of computation
      call dealloc_geometry(g)
	  !call dealloc_b2mod_geo(g)
	  call dealloc_mapping(m)
	  !call dealloc_b2mod_indirect(m)

      else
	  write (*,*) 'Structured b2ag used'
*   ..rewind b2ag.dat file
	  rewind ninp(0)
*   ..obtain nx, ny, nx1, ny1
      	call b2ag_read_parameters_st(ninp,nout, nx, ny, nx1, ny1)
*   ..	call principal driver routine
      	call alloc_b2mod_indirect_st(nx,ny,nncutmax)
      	call alloc_b2mod_geo_st(nx,ny)
      	call alloc_b2mod_geo2_st(nx1,ny1,nncutmax)
      	call b2agdr_st (ninp, nout, nx, ny, nx1, ny1)
*   ..close graphics files
*     	(none)
*   ..close other files
      	do io=0,size(ninp)-1
         	inquire(unit=ninp(io),opened=file_open)
         	if (file_open) close(ninp(io))
      	enddo
      	do io=0,size(nout)-1
         	inquire(unit=nout(io),opened=file_open)
         	if (file_open) close(nout(io))
      	enddo

* ..end of computation
      	call dealloc_b2mod_geo_st
      	call dealloc_b2mod_geo2_st
      	call dealloc_b2mod_indirect_st

      endif

      call prgend ()
      stop 'b2ag'

*-----------------------------------------------------------------------
*.end b2ag

      end program b2ag
