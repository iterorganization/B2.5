      program b2plasmastate_inspect
      use b2mod_types
      use b2mod_version
      use b2mod_time
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2us_io
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      use b2mod_xerset
      implicit none

*-----------------------------------------------------------------------
*     b2plasmastate_inspect reads in a formatted plasmastate file, and
*     displays the value of the plasma time stamp "tim"
*
*     Run as: b2plasmastate_inspect.exe plasmastate.xxxx
*-----------------------------------------------------------------------
      
#if ( defined(USE_PXFGETENV) || defined(NAGFOR) )
      integer lenval, ierror
#endif
#if ( !defined(NAGF95) && !defined(NAGFOR) )
      integer iargc
#endif
      integer :: ninp(0:5), nCv, nFc, ns, ns0, idum(0:9)
      character :: arg*256
*   ..main data structures
      type(mapping), save    :: mpg
      type(geometry), save   :: geo
      type(B2State), save    :: state
      type(switches), save   :: switch
*   ..procedures
      external xertst, cfopen, cfruin
      external b2rucp
*   ..initialise input/output unit numbers
      data ninp/50,51,52,53,54,55/

*-----------------------------------------------------------------------

      call prgini('b2plasmastate_inspect')
! Handle command args
#ifndef NAGFOR
      call xertst (iargc().ge.1,
     &  'program should be supplied with a state number')
#ifdef USE_PXFGETENV
      call pxfgetarg(1,arg,lenval,ierror)
#else
      call getarg(1,arg)
#endif
#else
      call xertst (command_argument_count().ge.1,
     &  'program should be supplied with a state number')
      call get_command_argument(1,arg,lenval,ierror)
#endif

! Open files
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fpardf','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4), trim(arg), 'old','unformatted')

! Follow b2co logic      
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version numbers                                           !xpb
      call cfverr (ninp(1), b2fgmtry_version)                          !xpb
      call cfverr (ninp(2), b2fparam_version)                          !xpb
      call cfverr (ninp(3), b2fstati_version)                          !xpb
*   ..obtain mapping/geometry dimensions from b2fgmtry file
      call cfruin (ninp(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      mpg%nCi = idum(0) 
      mpg%nCg = idum(1) 
      mpg%nCv = idum(2)
      mpg%nFc = idum(3)
      mpg%nVx = idum(4)
      mpg%nFs = idum(5)
      mpg%nFt = idum(6)
      call cfruin (ninp(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      mpg%nCmxVx = idum(0)
      mpg%nCmxFc = idum(1)
      mpg%nVmxCv = idum(2)
      mpg%nVmxFc = idum(3)
      mpg%nCmxNv = idum(4)
      call cfruin (ninp(1),1,idum,'isClassicalGrid')
      mpg%isClassicalGrid = idum(0)
      call cfruin (ninp(1),3,idum,'nx,ny,nncut')
      mpg%nx    = idum(0)
      mpg%ny    = idum(1)
      mpg%nncut = idum(2)
      call xertst (0.le.mpg%nCv.and.0.le.mpg%nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain ns from parameters file
      call cfruin (ninp(2), 1, idum, 'ns')
      ns = idum(0)
      call xertst (1.le.ns, 'faulty input ns from parameters file')
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns0 = idum(2)
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc, ns from b2fstate file')
      call xertst (1.le.ns0.and.ns0.le.ns,
     &  'faulty input ns from plasma state file')
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs, 
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
!*   ..obtain lblmn
!      call b2xxid ('b2mn', lblmn(1:60))
!      call b2rflb (ninp(0), lblmn(61:120))
!      label = lblmn
!*   ..read physics common
!      call cfruch (ninp(2), 120, lblcp, 'label')
!      call b2rucp (ninp(2), b2fparam_version, ns)
!*   ..read revised parameters in physics common
!*     (ninp(0) need supply only parameters that are to be changed from
!*     the values read in as default physics common.)
!      call b2rfcp (ninp(0), ns)
!*   ..test parameters
!      call b2xvcp (ns)
!*   ..read and echo code internal parameters
!    1 continue
!       read (ninp(0),'(a)',end=2,err=91) line
!       if (.not.is_comment(line)) then
!         read (line,*) cnamip, cvalip
!         call ipsetc (cnamip, cvalip)
!       end if
!       goto 1
!    2 continue
!*   ..read header for geometry
!      call cfruch (ninp(1), 120, lblgm, 'label')
!      call cfruin (ninp(1), 1, idum, 'isymm')
!      isymm = idum(0)
!*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
!*   ..load default values of switches
!      call read_first_switches(switch)
!      call set_defaults_switches(switch,ns)
!*   ..overwrite with values provided in b2mn.dat
!      call read_switches(switch)
!*   ..check the values of the switches
!      call check_switches(switch,ns)
!*   ..read geometry and mapping
!      call read_b2fgmtry (ninp(1), mpg, geo)
!      call init_geometry (mpg, geo, switch)
*   ..read initial state
      call cfverr (ninp(4), plasmastate_version)                        !xpb
      call read_b2fstate (ninp(4), nCv, nFc, ns, state)

! Display tim      
      write(0,'(3a,1p,1e15.8,a)') 'The plasmastate file ',trim(arg),
     & ' corresponds to plasma time ',tim,' seconds.'

* ..end of computation
      call destroyB2State (state)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      call prgend ()

      stop 'b2plasmastate_inspect'
      end program b2plasmastate_inspect

