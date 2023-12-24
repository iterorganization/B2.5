*-----------------------------------------------------------------------
*.specification

      program b2yv
      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_ma28
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_elements
      use b2mod_switches
      use b2mod_xerset
      use b2mod_ipmain
      use b2us_io
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef WINDOWS
! AVDef is the AViz module file
      use AVDef    ! IGNORE
      use AVVIEWER ! IGNORE
      use DFLIB    ! IGNORE
#endif
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YV reads in an unformatted b2fplasma file and visualizes
*     b2plasma arrays.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks

*   ..main data structures
      type(mapping)  :: mpg
      type(geometry) :: geo
      type(B2State)  :: state
      type(switches) :: switch
*   ..local variables
      integer :: icase, iCv, is
      real(kind=R4), allocatable :: M(:)
      real(kind=R8), allocatable :: pov(:), gpor(:)
      character(1) :: key
      character cnamip*80, cvalip*80, line*500
      character lblgm*120
#ifdef WINDOWS
      integer :: hv, lbnd(1) = 0
      integer status
      character(AV_MAX_LABEL_LEN) :: axis1 = 'X', axis2 = 'Y', axis3
#else
      character*40 :: axis3
#endif

      integer ninp(0:5), nout(0:7), nCv, nFc, ns, idum(0:9), io
      character lblmn*120, title*64
      character hlp_frm*22
      logical file_open
*   ..procedures
      logical is_comment
      external is_comment
      external prgini, prgend, xertst, xerrab, cfopen, cfruin
*   ..initialize input/output unit numbers
      data ninp/50,51,52,53,54,55/
      data nout/60,61,62,63,64,65,66,67/

*-----------------------------------------------------------------------
*.computation

*   ..program start-up calls
      call prgini ('b2yv')
*   ..open files
      call cfopen (ninp(0),'b2mn.dat','old','formatted')
      call cfopen (ninp(1),'b2fgmtry','old','un*formatted')
      call cfopen (ninp(2),'b2fparam','old','un*formatted')
      call cfopen (ninp(3),'b2fstati','old','un*formatted')
      call cfopen (ninp(4),'b2fplasma','old','unformatted')
!      call cfopen (nout(0),'b2yv.prt','new','formatted')
*   ..mark output unit for error messages
      call xerset (0)
*   ..obtain version number                                             !srv 18.09.03
      call cfverr (ninp(1), b2fgmtry_version)
      call cfverr (ninp(2), b2fparam_version)
      call cfverr (ninp(3), b2fstati_version)                           !srv 18.09.03
      call cfverr (ninp(4), b2fplasma_version)
*   ..obtain nCv, nFc, ns from b2fstate file
      call cfruin (ninp(3), 3, idum, 'nCv,nFc,ns')
      nCv = idum(0)
      nFc = idum(1)
      ns  = idum(2)
      call xertst (0.le.nCv.and.0.le.nFc.and.1.le.ns,
     &  'faulty input nCv, nFc, ns from plasma state file')
*   ..obtain ns from b2fparam file
      call cfruin (ninp(2), 1, idum, 'ns')
      call xertst (idum(0).eq.ns,
     &  'faulty input ns from physics parameters file')
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
      call xertst (mpg%nCv.eq.nCv.and.mpg%nFc.eq.nFc,
     &  'faulty input nCv, nFc from geometry file')
*   ..obtain parameters from b2fplasmf file
      call cfruin (ninp(4),3,idum,'nCv,nFc,ns')
      call xertst (idum(0).eq.nCv.and.idum(1).eq.nFc.and.idum(2).eq.ns,
     &  'faulty input nCv, nFc, ns from b2fplasmf file')
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*   ..read plasma state file
      call cfruch (ninp(3), 120, lblmn, 'label')
*   ..read species data
      call b2ruzd (ninp(3), b2fstati_version,
     &  ns, zamin, zamax, zn, am, .true.)
      allocate (M(nCv))
*   ..call principal driver routine
      call alloc_mapping (mpg)
      call alloc_geometry (geo, mpg%nCv, mpg%nFc, mpg%nVx, mpg%nFs,
     .  mpg%nFt, mpg%nVmxCv)
      call createB2State (nCv, nFc, ns, state)
*   ..read and echo code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       end if
       goto 1
    2 continue
      call read_switch_afn(switch)
      call set_defaults_switches(switch,ns)
      call read_switches(switch)
      call check_switches(switch,ns)
*   ..read geometry and metric
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)
*   ..read plasma state
      call read_b2fplasma(ninp(4), nCv, nFc, ns, state)
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
        M(:) = real(state%pl%na(:,icase),R4)
      elseif(ns.le.icase .and. icase.le.ns-1+ns )then
        write(title,"('Parallel velocity of ',a2,i1)")
     &     adjustr(is_codes(icase-ns)),nint(zamin(icase-ns))
        write(axis3,"('Parallel velocity of ',a2,i1,', [m/c]')")
     &     adjustr(is_codes(icase-ns)),nint(zamin(icase-ns))
        M(:) = real(state%pl%ua(:,icase-ns),R4)
      elseif(icase .eq. ns+ns )then
        write(title,"('Temperature of electrons')")
        write(axis3,"('Temperature of electrons',', [eV]')")
        do iCv = 1,nCv
          M(iCv) = real(state%pl%te(iCv)/ev,R4)
        enddo
      elseif(icase .eq. ns+ns+1 )then
        write(title,"('Temperature of ions')")
        write(axis3,"('Temperature of ions',', [eV]')")
        do iCv = 1,nCv
          M(iCv) = real(state%pl%ti(iCv)/ev,R4)
        enddo
      elseif(icase .eq. ns+ns+2 )then
        write(title,"('Electric potential')")
        do iCv = 1,nCv
          M(iCv) = real(state%pl%po(iCv),R4)
        enddo
      elseif(icase .eq. ns+ns+3 )then
        write(title,"('Radial electric field')")
        allocate(pov(mpg%nVx))
        allocate(gpor(mpg%nCv))
        call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     *               state%pl%po, pov, gpor)
        M(1:nCv)=real(-gpor(1:nCv),R4)
        deallocate(pov)
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

!       call faglShow(M, status)
      call favShowWindow(hv,AV_TRUE, status)

!      call faglName(M, trim(title), status)

      write (*,*) ' Press 0 - 9 to continue',
     &            ' or press any other key to finish'
      key = GETCHARQQ()
      call favEndViewer(hv,status)
!       call faglClose(M, status)
      call faglEndWatch(M, status)
#else
      key = ' '
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

      call destroyB2State (state)
      call dealloc_geometry (geo)
      call dealloc_mapping (mpg)
      deallocate (M)
* ..end of computation
      call prgend ()
      stop 'b2yv'

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2yv

      end program b2yv
