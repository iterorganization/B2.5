*-----------------------------------------------------------------------
*.specification

      program b2rd
      use b2mod_types
      use b2mod_version
      use b2mod_plasma
      use b2mod_geo
      use b2mod_work
      use b2mod_math
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmgs
      use b2mod_xerset
      use b2mod_mdsplus
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2RD (main program) is a postprocessor for B2. It imports an
*     archived plasma state from the MDSPlus tree and outputs the
*     equivalent b2fstate and B2 geometry files.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description (see also routine b2cdca)
*
*     The complete program b2rd performs a fetch of the archived
*     result of a B2 calculation.
*     This program unit opens and closes the input/output units, and
*     may perform some other system-dependent operations. It calls
*     routine b2rddr to control the actual computation.
*
*     The input units are:
*     ninp(0): formatted; provides output control parameters.
*
*     The output units are:
*     nout(0): formatted; provides printed output.
*     nout(1): un*formatted; provides the plasma state.
*     nout(2): un*formatted; provides the geometry.
*
*     (See routine b2cdca for the meaning of 'un*formatted'.)
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
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..local variables
      integer ninp(0:0), nout(0:2)
      real (kind=R8), allocatable :: zero(:,:)
#ifdef MDSPLUS
#     include "mdslib.inc"
      integer idum(0:9), io, ishot, n2, nx, ny, ns
      integer status, retlen
      integer ix, iy, is
      character*120 lblmn, lblgm
      integer, allocatable :: idummy(:,:)
      real(kind=R8), allocatable :: rdum(:), rdummy(:,:),
     *     rfx(:,:), rfy(:,:),
     *     rfn(:,:,:), rfnx(:,:,:), rfny(:,:,:),
     *     rfn2(:,:,:,:)
      real (kind=R8), allocatable :: za(:)
      logical file_open
*   ..procedures
      external b2wups, b2wugm, b2wuzd, b2xgbs, xertst, xerrab
      external b2xvgm, b2xvps, b2xvsg
      external cfwuin, cfwuch
      external intfaceh1, intfacer1, intfacev1
#endif
      external prgini, prgend, cfopen
*   ..initialise input/output units
      data ninp/50/
      data nout/60,61,62/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     Description of some local variables:
*
*     ninp - (0:0) integer array.
*     ninp specifies the input unit numbers.
*
*     nout - (0:2) integer array.
*     nout specifies the output unit numbers.
*
*     nx, ny - integer.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2rd')
*   ..open files
      call cfopen (nout(0),'b2rd.prt','new','formatted')
      call cfopen (nout(1),'b2fstate','new','un*formatted')
      call cfopen (nout(2),'b2fgmtrt','new','un*formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain run number
#ifdef MDSPLUS
      open(unit=ninp(0),file='shotnumber.history',status='old')
 11   continue
      read(ninp(0),*,end=12) ishot
      goto 11
 12   continue
      call xertst(1.le.ishot, 'faulty run number')
      write(*,*) 'Attempting recovery from run number ',ishot
*   ..obtain archived run
      status=mdsconnect('solps-mdsplus.aug.ipp.mpg.de:8001'//char(0))
      write(*,*) 'mdsconnect'
      if(status.lt.0) then
        write(*,*) ' return status: ',status
      endif
      status=mdsopen('solps'//char(0),ishot)
      write(*,*) 'mdsopen'
      call mds_err_msg(status)

#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i0('_s=\TOP.IDENT:NX'//char(0),
     1 descr2_i0(idtype_long,0),nx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT:NX'//char(0),
     1 descr2(idtype_long,0),nx,0,retlen)
#endif
      call mds_err_msg(status)
      nx=nx-2
      write(*,*) 'Found NX = ',nx
      call xertst (0.le.nx, 'faulty parameter nx')
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i0('_s=\TOP.IDENT:NY'//char(0),
     1 descr2_i0(idtype_long,0),ny,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT:NY'//char(0),
     1 descr2(idtype_long,0),ny,0,retlen)
#endif
      call mds_err_msg(status)
      ny=ny-2
      write(*,*) 'Found NY = ',ny
      call xertst (0.le.ny, 'faulty parameter ny')
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i0('_s=\TOP.IDENT:NS'//char(0),
     1 descr2_i0(idtype_long,0),ns,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT:NS'//char(0),
     1 descr2(idtype_long,0),ns,0,retlen)
#endif
      call mds_err_msg(status)
      write(*,*) 'Found NS = ',ns
      call xertst (1.le.ns, 'faulty parameter ns')
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_str('_s=\TOP.IDENT.VERSION:b2fstate'//char(0),
     1 descr2_str0(idtype_cstring,0,10),
     2 plasmastate_version,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT.VERSION:b2fstate'//char(0),
     1 descr2(idtype_cstring,0,10),plasmastate_version,0,retlen)
#endif
      call mds_err_msg(status)
      write(*,*) 'Reading plasma state version number = ',
     1 plasmastate_version
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i0('_s=\TOP.IDENT:VERSION'//char(0),
     1 descr2_i0(idtype_long,0),dbversion,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT:VERSION'//char(0),
     1 descr2(idtype_long,0),dbversion,0,retlen)
#endif
      call mds_err_msg(status)
      write(*,*) 'MDSplus database version = ',dbversion

*   ..allocate arrays
      call alloc_b2mod_geo(nx,ny)
      call alloc_b2mod_work(nx,ny,ns)
      call alloc_b2mod_rates(nx,ny,ns)
      call alloc_b2mod_plasma(nx,ny,ns)
      call alloc_b2mod_indirect(nx,ny,nncutmax)
      allocate(idummy(-1:nx,-1:ny))
      allocate(rdummy(-1:nx,-1:ny))
      allocate(rdum(0:ns-1))
      allocate(rfx( 0:nx,-1:ny))
      allocate(rfy(-1:nx, 0:ny))
      allocate(rfn (-1:nx,-1:ny,0:ns-1))
      allocate(rfnx( 0:nx,-1:ny,0:ns-1))
      allocate(rfny(-1:nx, 0:ny,0:ns-1))
      allocate(rfn2(-1:nx,-1:ny,0:1,0:ns-1))
      allocate(zero(-1:nx,-1:ny))
      zero = 0.0_R8

*   ..read species data
      if (dbversion.ge.9) then
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r1('_s=\TOP.SNAPSHOT.GRID:ZAMIN'//char(0),
     1   descr2_i1(idtype_double,NS,0),rdum,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:ZAMIN'//char(0),
     1   descr2(idtype_double,NS,0),rdum,0,retlen)
#endif
        zamin=rdum
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r1('_s=\TOP.SNAPSHOT.GRID:ZAMAX'//char(0),
     1   descr2_i1(idtype_double,NS,0),rdum,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:ZAMAX'//char(0),
     1   descr2(idtype_double,NS,0),rdum,0,retlen)
#endif
        zamax=rdum
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2('_s=\TOP.SNAPSHOT:RZA'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:RZA'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#endif
        rza=rfn
        call mds_err_msg(status)
      else
        allocate(za(0:ns-1))
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r1('_s=\TOP.SNAPSHOT.GRID:ZA'//char(0),
     1   descr2_i1(idtype_double,NS,0),za,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:ZA'//char(0),
     1   descr2(idtype_double,NS,0),za,0,retlen)
#endif
        call mds_err_msg(status)
        zamin(0:ns-1)=za(0:ns-1)
        zamax(0:ns-1)=za(0:ns-1)
        do is = 0, ns-1
          rza(-1:nx,-1:ny,is) = za(is)
        end do
        deallocate(za)
      endif
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r1('_s=\TOP.SNAPSHOT.GRID:ZN'//char(0),
     1 descr2_i1(idtype_double,NS,0),rdum,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:ZN'//char(0),
     1 descr2(idtype_double,NS,0),rdum,0,retlen)
#endif
      zn=rdum
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r1('_s=\TOP.SNAPSHOT.GRID:AM'//char(0),
     1 descr2_i1(idtype_double,NS,0),rdum,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:AM'//char(0),
     1 descr2(idtype_double,NS,0),rdum,0,retlen)
#endif
      am=rdum
      call mds_err_msg(status)

*   ..read plasma data
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:NE'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:NE'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      ne=rdummy
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:NA'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:NA'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#endif
      na=rfn
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:UA'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:UA'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#endif
      ua=rfn
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:TE'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:TE'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      te=rdummy*qe
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:TI'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:TI'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      ti=rdummy*qe
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:PO'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:PO'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      po=rdummy
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX'//char(0),
     1 descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
      fna(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
      fna(-1,:,0,0,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
      fna(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
      fna(:,-1,1,1,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_32'//char(0),
     1 descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_32'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
      fna_32(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
      fna_32(-1,:,0,0,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_32'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_32'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
      fna_32(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
      fna_32(:,-1,1,1,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_52'//char(0),
     1 descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_52'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
      fna_52(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
      fna_52(-1,:,0,0,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_52'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_52'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
      fna_52(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
      fna_52(:,-1,1,1,:)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEX'//char(0),
     1 descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEX'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
      fhe(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
      fhe(-1,:,0,0)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEY'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEY'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
      fhe(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
      fhe(:,-1,1,1)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIX'//char(0),
     1 descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIX'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
      fhi(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
      fhi(-1,:,0,0)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIY'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIY'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
      fhi(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
      fhi(:,-1,1,1)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHX'//char(0),
     1 descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHX'//char(0),
     1 descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
      fch(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
      fch(-1,:,0,0)=0.0_R8
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHY'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHY'//char(0),
     1 descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
      fch(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
      fch(:,-1,1,1)=0.0_R8
      call mds_err_msg(status)
      if (dbversion.le.9) then
        uadia=0.0_R8
        fch_32(:,:,0,0)=fch(:,:,0,0)
        fch_32(:,:,1,1)=0.0_R8
        fch_52(:,:,0,0)=0.0_R8
        fch_52(:,:,1,1)=fch(:,:,1,1)
        fch_p(:,:,0,0)=fch(:,:,0,0)                               !srv 10.06.09
        fch_p(:,:,1,1)=0.0_R8                                   !srv 10.06.09
        do is = 0, ns-1
          kinrgy(:,:,is)=0.5_R8*am(is)*mp*ua(:,:,is)**2
        enddo
      else
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r4('_s=\TOP.SNAPSHOT:UADIA'//char(0),
     1   descr2_i4(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:UADIA'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#endif
        uadia(-1:nx,-1:ny,0,0,0:ns-1)=rfn2(-1:nx,-1:ny,0,0:ns-1)
        uadia(-1:ns,-1:ny,1,1,0:ns-1)=rfn2(-1:nx,-1:ny,1,0:ns-1)
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHX_32'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHX_32'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fch_32(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fch_32(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHY_32'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHY_32'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fch_32(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fch_32(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHX_52'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHX_52'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fch_52(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fch_52(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHY_52'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHY_52'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fch_52(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fch_52(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCH_P'//char(0),  !srv 10.06.09 {
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCH_P'//char(0),  !srv 10.06.09 {
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fch_p(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fch_p(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:KINRGY'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:KINRGY'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,NS,0),rfn,0,retlen)
#endif
        kinrgy=rfn
        call mds_err_msg(status)
      endif
      if (dbversion.ge.22) then
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_MDF'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_MDF'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fna_mdf(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fna_mdf(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_MDF'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_MDF'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fna_mdf(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fna_mdf(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEX_MDF'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEX_MDF'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhe_mdf(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhe_mdf(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEY_MDF'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEY_MDF'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhe_mdf(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhe_mdf(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIX_MDF'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIX_MDF'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhi_mdf(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhi_mdf(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIY_MDF'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIY_MDF'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhi_mdf(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhi_mdf(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_FCOR'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_FCOR'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fna_fcor(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fna_fcor(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_FCOR'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_FCOR'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fna_fcor(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fna_fcor(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_NODRIFT'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_NODRIFT'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fna_nodrift(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fna_nodrift(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_NODRIFT'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_NODRIFT'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fna_nodrift(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fna_nodrift(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_HE'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_HE'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fna_he(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fna_he(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_HE'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_HE'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fna_he(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fna_he(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAPSCHX'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAPSCHX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fnaPSch(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fnaPSch(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAPSCHY'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAPSCHY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fnaPSch(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fnaPSch(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEPSCHX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEPSCHX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhePSch(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhePSch(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEPSCHY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEPSCHY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhePSch(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhePSch(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIPSCHX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIPSCHX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhiPSch(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhiPSch(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIPSCHY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIPSCHY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhiPSch(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhiPSch(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAX_EIR'//char(0),
     1   descr2_i3(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAX_EIR'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,NS,0),rfnx,0,retlen)
#endif
        fna_eir(0:nx,-1:ny,0,0,0:ns-1)=rfnx(0:nx,-1:ny,0:ns-1)
        fna_eir(-1,:,0,0,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:FNAY_EIR'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNAY_EIR'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,NS,0),rfny,0,retlen)
#endif
        fna_eir(-1:nx,0:ny,1,1,0:ns-1)=rfny(-1:nx,0:ny,0:ns-1)
        fna_eir(:,-1,1,1,:)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEX_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEX_EIR'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fne_eir(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fne_eir(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEY_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEY_EIR'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fne_eir(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fne_eir(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEX_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEX_EIR'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhe_eir(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhe_eir(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHEY_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHEY_EIR'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhe_eir(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhe_eir(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIX_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIX_EIR'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fhi_eir(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fhi_eir(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FHIY_EIR'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FHIY_EIR'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fhi_eir(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fhi_eir(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNIX_32'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNIX_32'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fni_32(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fni_32(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNIY_32'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNIY_32'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fni_32(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fni_32(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNIX_52'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNIX_52'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fni_52(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fni_52(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNIY_52'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNIY_52'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fni_52(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fni_52(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEX_32'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEX_32'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fne_32(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fne_32(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEY_32'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEY_32'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fne_32(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fne_32(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEX_52'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEX_52'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fne_52(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fne_52(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FNEY_52'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FNEY_52'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fne_52(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fne_52(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHDIAX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHDIAX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchdia(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchdia(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHDIAY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHDIAY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchdia(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchdia(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHINX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHINX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchin(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchin(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHINY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHINY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchin(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchin(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISPARX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISPARX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchvispar(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchvispar(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISPARY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISPARY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchvispar(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchvispar(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISPERX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISPERX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchvisper(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchvisper(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISPERY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISPERY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchvisper(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchvisper(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISQX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISQX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchvisq(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchvisq(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHVISQY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHVISQY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchvisq(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchvisq(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHINERTX'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHINERTX'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        fchinert(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        fchinert(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FCHINERTY'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FCHINERTY'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        fchinert(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        fchinert(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r4('_s=\TOP.SNAPSHOT:VADIA'//char(0),
     1   descr2_i4(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:VADIA'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#endif
        vadia(-1:nx,-1:ny,0,0,0:ns-1)=rfn2(-1:nx,-1:ny,0,0:ns-1)
        vadia(-1:nx,-1:ny,1,1,0:ns-1)=rfn2(-1:nx,-1:ny,1,0:ns-1)
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r4('_s=\TOP.SNAPSHOT:WADIA'//char(0),
     1   descr2_i4(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:WADIA'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#endif
        wadia(-1:nx,-1:ny,0,0,0:ns-1)=rfn2(-1:nx,-1:ny,0,0:ns-1)
        wadia(-1:nx,-1:ny,1,1,0:ns-1)=rfn2(-1:nx,-1:ny,1,0:ns-1)
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r4('_s=\TOP.SNAPSHOT:VAECRB'//char(0),
     1   descr2_i4(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:VAECRB'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,NS,0),rfn2,0,retlen)
#endif
        vaecrb(-1:nx,-1:ny,0,0,0:ns-1)=rfn2(-1:nx,-1:ny,0,0:ns-1)
        vaecrb(-1:nx,-1:ny,1,1,0:ns-1)=rfn2(-1:nx,-1:ny,1,0:ns-1)
        call mds_err_msg(status)
        deallocate(rfn)
        allocate(rfn(-1:nx,-1:ny,0:1))
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:VEDIA'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:VEDIA'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#endif
        vedia(-1:nx,-1:ny,0,0)=rfn(-1:nx,-1:ny,0)
        vedia(-1:nx,-1:ny,1,1)=rfn(-1:nx,-1:ny,1)
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:VEECRB'//char(0),
     1   descr2_i3(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:VEECRB'//char(0),
     1   descr2(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#endif
        veecrb(-1:nx,-1:ny,0,0)=rfn(-1:nx,-1:ny,0)
        veecrb(-1:nx,-1:ny,1,1)=rfn(-1:nx,-1:ny,1)
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FLOEX_NOC'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FLOEX_NOC'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        floe_noc(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        floe_noc(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FLOEY_NOC'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FLOEY_NOC'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        floe_noc(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        floe_noc(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FLOIX_NOC'//char(0),
     1   descr2_i2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FLOIX_NOC'//char(0),
     1   descr2(idtype_double,NX+1,NY+2,0),rfx,0,retlen)
#endif
        floi_noc(0:nx,-1:ny,0,0)=rfx(0:nx,-1:ny)
        floi_noc(-1,:,0,0)=0.0_R8
        call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
        status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:FLOIY_NOC'//char(0),
     1   descr2_i2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#else
        status=mdsvalue2('_s=\TOP.SNAPSHOT:FLOIY_NOC'//char(0),
     1   descr2(idtype_double,NX+2,NY+1,0),rfy,0,retlen)
#endif
        floi_noc(-1:nx,0:ny,1,1)=rfy(-1:nx,0:ny)
        floi_noc(:,-1,1,1)=0.0_R8
        call mds_err_msg(status)
      else
        fna_mdf = fna
        fhe_mdf = fhe
        fhi_mdf = fhi
        fna_fcor = fna
        fna_nodrift = fna
        fna_he = fna
        fnaPSch = 0.0_R8
        fhePSch = 0.0_R8
        fhiPSch = 0.0_R8
        fna_eir = fna
        fhe_eir = fhe
        fhi_eir = fhi
        call b2xpfi (nx, ny, ns, fna_32, fni_32)
        call b2xpfi (nx, ny, ns, fna_52, fni_52)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_32, zero, fch_32, fne_32)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_52, zero, fch_52, fne_52)
        call b2xpfe (nx, ny, ns, 0, qe,
     1     rza, zero, fna_eir, zero, fch, fne_eir)
        fchdia = 0.0_R8
        fchin = 0.0_R8
        fchvispar = 0.0_R8
        fchvisper = 0.0_R8
        fchvisq = 0.0_R8
        fchinert = 0.0_R8
        vaecrb = 0.0_R8
        vadia = 0.0_R8
        wadia = 0.0_R8
        veecrb = 0.0_R8
        vedia = 0.0_R8
        floe_noc = 0.0_R8
        floi_noc = 0.0_R8
      endif

*   ..read geometry data
      deallocate(rfn)
      allocate(rfn(-1:nx,-1:ny,0:3))
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT.GRID:R'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:R'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#endif
      crx=rfn
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT.GRID:Z'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:Z'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#endif
      cry=rfn
      call mds_err_msg(status)
      fpsi = 0.0_R8
      ffbz = 0.0_R8
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:B'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:B'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,4,0),rfn,0,retlen)
#endif
      bb=rfn
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:VOL'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:VOL'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      vol=rdummy
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:HX'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:HX'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      hx=rdummy
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r2('_s=\TOP.SNAPSHOT:HY'//char(0),
     1 descr2_i2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:HY'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,0),rdummy,0,retlen)
#endif
      hy=rdummy
      call mds_err_msg(status)
      deallocate(rfn)
      allocate(rfn(-1:nx,-1:ny,0:1))
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_r3('_s=\TOP.SNAPSHOT:QZ'//char(0),
     1 descr2_i3(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT:QZ'//char(0),
     1 descr2(idtype_double,NX+2,NY+2,2,0),rfn,0,retlen)
#endif
      qz=rfn
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:REGVOL'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:REGVOL'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      region(:,:,0)=idummy(:,:)
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:REGFLX'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:REGFLX'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      region(:,:,1)=idummy(:,:)
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:REGFLY'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:REGFLY'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      region(:,:,2)=idummy(:,:)
      call mds_err_msg(status)
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:LEFTIX'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:LEFTIX'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      leftix = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:LEFTIY'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:LEFTIY'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      leftiy = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:RIGHTIX'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:RIGHTIX'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      rightix = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:RIGHTIY'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:RIGHTIY'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      rightiy = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:BOTTOMIX'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:BOTTOMIX'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      bottomix = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:BOTTOMIY'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:BOTTOMIY'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      bottomiy = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:TOPIX'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:TOPIX'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      topix = idummy - 1
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i2('_s=\TOP.SNAPSHOT.GRID:TOPIY'//char(0),
     1 descr2_i2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.SNAPSHOT.GRID:TOPIY'//char(0),
     1 descr2(idtype_long,NX+2,NY+2,0),idummy,0,retlen)
#endif
      call mds_err_msg(status)
      topiy = idummy - 1
      nncut = maxval(region(:,:,0)/4)
      do is = 0, 2
        nnreg(is) = maxval(region(:,:,is))
      enddo
      is = 0
      iy = -1
      do ix = -1, nx
        if (rightix(ix,iy).ne.ix+1 .and. rightix(ix,iy).ne.nx+2) then
          is = is + 1
          write(*,*)
     .     'Found cut number ',is,' at ix location ',ix, rightix(ix,iy)
          leftcut(is) = ix+1
          rightcut(is) = rightix(ix,iy)
          bottomcut(is) = -1
        endif
      enddo
      call xertst (is.eq.2*nncut, 'Did not find all the cuts!')
      do is = 1, nncut
        iy = bottomcut(is)
        ix = leftcut(is)-1
        do while(iy.le.ny .and. rightix(ix,iy).ne.ix+1)
          topcut(is) = iy+1
          iy = iy+1
        enddo
      enddo
      do iy=-1,ny
        do ix=-1,nx
          if(region(ix,iy,0).eq.0 .or.
     &     leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
            resignore(ix,iy,1)=0
            resignore(ix,iy,2)=0
          else
            resignore(ix,iy,1)=1
            resignore(ix,iy,2)=1
          endif
        enddo
      enddo
      periodic_bc = 0
      do iy = -1, ny
        do ix = -1, nx
          qc(ix,iy,0) = cosang (
     &     (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &     (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &      crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
        enddo
      enddo
      do iy = -1, ny
        do ix = -1, nx
          gs(ix,iy,2) = 0.5_R8*
     &     ((crx(ix,iy,3)-crx(ix,iy,0))*(cry(ix,iy,2)-cry(ix,iy,1))-
     &      (crx(ix,iy,2)-crx(ix,iy,1))*(cry(ix,iy,3)-cry(ix,iy,0)))
          if(leftix(ix,iy).eq.-2) then
            gs(ix,iy,0) = 0.0_R8
          else
            gs(ix,iy,0) = pi*dist(ix,iy,0,2)*(crx(ix,iy,0)+crx(ix,iy,2))
          endif
          if(bottomiy(ix,iy).eq.-2) then
            gs(ix,iy,1) = 0.0_R8
          else
            gs(ix,iy,1) = pi*dist(ix,iy,0,1)*(crx(ix,iy,0)+crx(ix,iy,1))
          endif
        enddo
      enddo
      nlreg = 0
      call intfaceh1(nx,ny,2,leftix,leftiy,vol,bb,wbbl)
      call intfacer1(nx,ny,2,rightix,rightiy,vol,bb,wbbr)
      call intfacev1(nx,ny,2,bottomix,bottomiy,vol,bb,wbbv)
      do iy = -1, ny
        do ix = -1, nx
          if (leftix(ix,iy).eq.-2) wbbl(ix,iy,0:2) =
     &     bb(ix,iy,0:2) - (wbbr(ix,iy,0:2)-bb(ix,iy,0:2))
          if (rightix(ix,iy).eq.nx+1) wbbr(ix,iy,0:2) =
     &     bb(ix,iy,0:2) - (wbbl(ix,iy,0:2)-bb(ix,iy,0:2))
          if (bottomiy(ix,iy).eq.-2) wbbv(ix,iy,0:2) = bb(ix,iy,0:2) -
     &     (wbbv(topix(ix,iy),topiy(ix,iy),0:2)-bb(ix,iy,0:2))
        enddo
      enddo
      wbbl(:,:,3)=sqrt(wbbl(:,:,0)**2+wbbl(:,:,1)**2+wbbl(:,:,2)**2)
      wbbr(:,:,3)=sqrt(wbbr(:,:,0)**2+wbbr(:,:,1)**2+wbbr(:,:,2)**2)
      wbbv(:,:,3)=sqrt(wbbv(:,:,0)**2+wbbv(:,:,1)**2+wbbv(:,:,2)**2)
      call b2xgbs (nx, ny, wbbl, wbbv, gs, qc, crx, cry, pbs)

*   .. copied from b2rugm
      cell_width = hx
      cell_height(:,:) = hy(:,:)*qz(:,:,1)
      do iy = -1, ny
        do ix = -1, nx
          qc(ix,iy,1) = cosang (
     &       crx(ix,iy,1)-crx(ix,iy,0),  cry(ix,iy,1)-cry(ix,iy,0),
     &      (crx(ix,iy,2)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,1)),
     &      (cry(ix,iy,2)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,1)))
          if(qc(ix,iy,1).le.0.0_R8) then
            write(*,'(a,1f12.6,2i5)')
     1       'Problem in b2rugm: qc(ix,iy,1) <= 0',qc(ix,iy,1),ix,iy
            write(*,'(a,4f12.6)') 'crx = ', crx(ix,iy,0:3)
            write(*,'(a,4f12.6)') 'cry = ', cry(ix,iy,0:3)
          endif
        enddo
      enddo
      do iy = -1, ny
        do ix = -1, nx
          qs(ix,iy,0) = sinang (
     &     (crx(ix,iy,1)+crx(ix,iy,3))-(crx(ix,iy,0)+crx(ix,iy,2)),
     &     (cry(ix,iy,1)+cry(ix,iy,3))-(cry(ix,iy,0)+cry(ix,iy,2)),
     &      crx(ix,iy,2)-crx(ix,iy,0), cry(ix,iy,2)-cry(ix,iy,0))
          if (bottomiy(ix,iy).eq.-2) then
            qs(ix,iy,1)=0.0_R8
          else
            qs(ix,iy,1) =
     .      (hy(ix,iy)*qz(ix,iy,0)+hy(bottomix(ix,iy),bottomiy(ix,iy))*
     .       qz(bottomix(ix,iy),bottomiy(ix,iy),0))/
     .      (hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))
          endif
        enddo
      enddo
cwdk  ensure qc(:,:,:)**2 + qs(:,:,:)**2 = 1, in possible event
cwdk  of improved computation of qc
      qs(:,:,0) = sign(sqrt(1.0_R8 - qc(:,:,0)**2),qs(:,:,0))
      qs(:,:,1) = sign(sqrt(1.0_R8 - qc(:,:,1)**2),qs(:,:,1))

*   ..test geometry and plasma state
      n2 = (nx+2)*(ny+2)
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qs, gs, pbs,
     &  cell_width, cell_height)
      call b2xvps (nx, ny, ns,
     .  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
      call b2xvsg (n2*ns, kinrgy, 1, 'kinrgy', '.ge.')
*   ..save final state
      newversion = plasmastate_version
      write(lblmn,*) 'Recovered from run number ',ishot
      idum(0) = nx
      idum(1) = ny
      idum(2) = ns
      call cfverw (nout(1), plasmastate_version)                         !xpb
      call cfwuin (nout(1), 3, idum, 'nx,ny,ns')
      call cfwuch (nout(1), 120, lblmn, 'label')
      call b2wuzd (nout(1), plasmastate_version, ns,
     &  zamin, zamax, zn, am)
      call b2wups (nout(1), nx, ny, ns, ne, na, ua, uadia, te,           !xpb
     &  ti, po, fna, fhe, fhi, fch, fch_32, fch_52, fch_p, kinrgy,
     &  fna_mdf, fhe_mdf, fhi_mdf, fna_fcor, fna_nodrift, fna_he,
     &  fnaPSch, fhePSch, fhiPSch, fna_eir, fne_eir, fhe_eir, fhi_eir,
     &  fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fchdia, fchin, fchvispar, fchvisper, fchvisq, fchinert,
     &  vaecrb, vadia, wadia, veecrb, vedia, floe_noc, floi_noc)
      write(nout(0),*) 'Converted run number ',ishot,' to b2fstate'
*   ..save geometry
      write(lblgm,*) 'Recovered geometry from run number ',ishot
      idum(0) = nx
      idum(1) = ny
      call cfverw (nout(2), newversion)
      call cfwuin (nout(2), 2, idum, 'nx,ny')
      call cfwuch (nout(2), 120, lblgm, 'label')
      idum(0) = 1
      call cfwuin (nout(2), 1, idum, 'isymm')
      call b2wugm (nout(2), nx, ny)

      open(91,file='core_boundary.txt')
      iy=-1
      do ix=-1,nx
         if(mod(region(ix,iy,0),4).eq.1) then
            write(91,*) crx(ix,iy,2),cry(ix,iy,2)
         endif
      enddo
      close(91)
      open(91,file='separatrix_boundary.txt')
      iy=minval(topcut(1:nncut))-1
      write(*,*) 'sep = ',iy
      do ix=-1,nx
         if(mod(region(ix,iy,0),4).eq.1) then
            write(91,*) crx(ix,iy,2),cry(ix,iy,2)
         endif
      enddo
      close(91)

*   ..wrapping up
      write(nout(0),*) 'Recovered geometry from run number ',ishot
      write(*,*) 'Recovered plasma state from run number ',ishot
#if ( defined(GFORTRAN) && COMPILER_MAJOR_VERSION > 9 ) || defined(G95)
      status=mdsvalue2_i0('_s=\TOP.IDENT:VERSION'//char(0),
     1 descr2_i0(idtype_long,0),dbversion,0,retlen)
#else
      status=mdsvalue2('_s=\TOP.IDENT:VERSION'//char(0),
     1 descr2(idtype_long,0),dbversion,0,retlen)
#endif
      call mds_err_msg(status)
      write(*,*) 'MDSplus database version = ',dbversion

*   ..close files
      status=mdsclose('solps'//char(0),ishot)
      write(*,*) 'mdsclose'
      call mds_err_msg(status)
      status=mdsdisconnect()
      write(*,*) 'mdsdisconnect'

      do io=0,size(ninp)-1
         inquire(unit=ninp(io),opened=file_open)
         if (file_open) close(ninp(io))
      enddo
      do io=0,size(nout)-1
         inquire(unit=nout(io),opened=file_open)
         if (file_open) close(nout(io))
      enddo

* ..end of computation
      deallocate(idummy)
      deallocate(rdum,rdummy,rfx,rfy,rfn,rfnx,rfny,rfn2)
      deallocate(zero)
      call dealloc_b2mod_geo
      call dealloc_b2mod_work
      call dealloc_b2mod_rates
      call dealloc_b2mod_plasma
      call dealloc_b2mod_indirect
#else
      call xerrab('b2rd: ' //
     1 'MDSPLUS routine called when not compiled with MDSPLUS')
#endif
      call prgend ()
      stop 'b2rd'

*-----------------------------------------------------------------------
*.end b2rd

      end program b2rd
