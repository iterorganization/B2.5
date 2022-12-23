      subroutine out_for_plot (ns,mpg, gm, st)

      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2us_plasma
      implicit none

      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2state), intent(in) :: st
      integer :: iout, igeo, ielm
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iout,file='idl.conf',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
#else
      iout = newunit()
      open (unit=iout,file='idl.conf',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
#endif
      call out_config(iout, mpg)
      close (iout)
#ifdef F2003
      open (newunit=igeo,file='geometry.npco_char',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
      open (newunit=ielm,file='geometry.elemente',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
#else
      igeo = newunit()
      open (unit=igeo,file='geometry.npco_char',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
      ielm = newunit()
      open (unit=ielm,file='geometry.npco_char',form='FORMATTED',
     .            ACCESS='SEQUENTIAL')
#endif
      call out_geometry(igeo, ielm, mpg, gm)
      close (igeo)
      close (ielm)

      call out_quantities(ns, mpg, gm, st)

      return
      end subroutine out_for_plot


      subroutine out_config(iout, mpg)
      use b2mod_types
      use b2us_map
      use b2mod_diag

      implicit none

      integer, intent(in) :: iout
      type(mapping), intent(in) :: mpg
      integer :: i1, i2, i3
      character(10) :: CDATE, CTIME
      character(8) :: CNAME

      write (iout,'(A)') ' PRODUCING PROGRAM: SOLPS-ITER'

      write (iout,'(//A,I6)') ' LEVGEO ',15

!  mesh consisting of n-cornered cells
      CNAME = 'geometry'
      write (iout,'(A,A)')  ' CASENAME ',CNAME
      write (iout,'(A,I10)') ' NVx     ',mpg%nVx
      write (iout,'(A,I10)') ' NCI     ',mpg%nCi
      write (iout,'(A,I10)') ' NCORNER ',maxval(mpg%cvVxp(1:mpg%nCi,2))
      
      write (iout,'(//A)') label

      CALL DATE_AND_TIME(CDATE,CTIME)
      READ(CDATE(1:4),*) I1
      READ(CDATE(5:6),*) I2
      READ(CDATE(7:8),*) I3
      WRITE (iout,'(//1X,A6,1X,2(I2,1X),I4)') 'DATE: ',I3,I2,I1
      READ(CTIME(1:2),*) I1
      READ(CTIME(3:4),*) I2
      READ(CTIME(5:6),*) I3
      WRITE (iout,'(1X,A6,1X,3(I2,1X))') 'TIME: ',I1,I2,I3
      
      return
      end subroutine out_config


      subroutine out_geometry(igeo, ielm, mpg, gm)
      use b2mod_types
      use b2us_map
      use b2us_geo
      use b2mod_diag

      implicit none

      integer, intent(in) :: igeo, ielm
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer :: icv, ivx, mxfc, nfcl, ifc, jfc, ihelp, ilast, iv, ichk
      integer, allocatable :: indvx(:)
      real (kind=R8) :: at
      character(6) :: ccell
      character(6) :: cform

!  write vertices onto file

      write (igeo,'(i9)') mpg%nVx
      do ivx = 1, mpg%nVx
        WRITE(igeo,'(I9,1P,2E15.7)') ivx, gm%vxX(ivx), gm%vxY(ivx)
      end do

!  for each plasma cell write the numbers of the vertices describing
!  the cell in counter-clockwise sequence
      
      mxfc = maxval(mpg%cvFcP(1:mpg%nCi,2))
      allocate (indvx(2*mxfc))

      write (ielm,'(2i9)') mpg%nCi, mxfc

      do icv = 1, mpg%nCi
        nfcl = mpg%cvFcP(icv,2)
        do jfc = 1, nfcl
          ifc = mpg%cvFc(mpg%cvFcP(icv,1)+jfc-1)
          indvx((jfc-1)*2+1) = mpg%fcVx(ifc,1)
          indvx(jfc*2) = mpg%fcVx(ifc,2)
        end do
!  check first face (first 2 points) for counter-clockwise sequence
        at = areat(gm%vxX(indvx(1)),gm%vxY(indvx(1)),
     .             gm%vxX(indvx(2)),gm%vxY(indvx(2)),
     .             gm%cvX(icv),gm%cvY(icv))
        if (at < 0._R8) then
!  wrong sequence
!  switch first two points
          ihelp = indvx(1)
          indvx(1) = indvx(2)
          indvx(2) = ihelp
        end if 

        ilast = 2
        do while(ilast < nfcl*2)
          do iv = ilast+1, nfcl*2
            if (indvx(iv) == indvx(ilast)) then
              if (mod(iv,2) == 1) then
                if (iv /= ilast+1) then
                  ihelp = indvx(ilast+1)
                  indvx(ilast+1) = indvx(iv)
                  indvx(iv) = ihelp
                  ihelp = indvx(ilast+2)
                  indvx(ilast+2) = indvx(iv+1)
                  indvx(iv+1) = ihelp
                end if
              else
                ihelp = indvx(ilast+1)
                indvx(ilast+1) = indvx(iv)
                indvx(iv) = ihelp
                if (iv-ilast-1 > 1) then
                  ihelp = indvx(ilast+2)
                  indvx(ilast+2) = indvx(iv-1)
                  indvx(iv-1) = ihelp
                end if 
              end if
              ilast = ilast + 2
              exit
            end if
          end do  ! do iv
        end do ! do while
!  now the vertices should be in correct order, check

        do iv = 2, nfcl*2, 2
          ichk = iv + 1
          if (ichk > nfcl*2) ichk = 1
          if (indvx(iv) /= indvx(ichk)) then
            write (ccell,'(i6)') icv
            call xerrab 
     .       ('out_geometry -- wrong vertex sequence,ncell='//ccell)
          end if
        end do

!  write vertex numbers to file
        write (cform,'(a,i2,a)') '(',nfcl+2,'i9)'
        write (ielm,cform) icv, nfcl, indvx(1:nfcl*2:2)
      end do

      return

      contains

      function areat(x1, y1, x2, y2, x3, y3)
      real (kind=R8) :: areat
      real (kind=R8), intent(in) :: x1, y1, x2, y2, x3, y3
        areat = 0.5_R8*(x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2))
      end function areat

      end subroutine out_geometry

      subroutine out_quantities(ns, mpg, gm, st)
      use b2mod_types
      use b2mod_constants
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_diag

      implicit none

      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2state), intent(in) :: st

      character(50) :: txtspc(0:ns-1)
      logical lfirst
      real(R8), allocatable :: dummy(:,:)
      real(R8) :: evi
      integer :: is

      do is = 0, ns-1
        txtspc = repeat(' ',50)
      end do
      evi = 1.0_R8 / ev

      allocate (dummy(mpg%nCi,0:ns-1))

      lfirst = .true.
      dummy(1:mpg%nCi,0) = st%pl%te(1:mpg%nCi) * eVi
      txtspc(0) = 'D'
      call write_quantity(1, mpg%nCi, 1, 0, dummy, lfirst,
     .                   'Te',txtspc,'eV')

      lfirst = .true.
      dummy(1:mpg%nCi,0) = st%pl%ti(1:mpg%nCi) * eVi
      txtspc(0) = 'D'
      call write_quantity(1, mpg%nCi, 2, 0, dummy, lfirst,
     .                   'Ti',txtspc,'eV')

      lfirst = .true.
      dummy(1:mpg%nCi,0:ns-1) = st%pl%na(1:mpg%nCi,0:ns-1)
      txtspc(0:ns-1) = 'D'
      call write_quantity(ns, mpg%nCi, 3, 0, dummy, lfirst,
     .                   'Ni',txtspc,'#/m**3')

      lfirst = .true.
      dummy(1:mpg%nCi,0:ns-1) = st%pl%ua(1:mpg%nCi,0:ns-1)
      txtspc(0:ns-1) = 'D'
      call write_quantity(ns, mpg%nCi, 4, 0, dummy, lfirst,
     .                   'Ua',txtspc,'m/s')

      deallocate (dummy)
      
      return
      end subroutine out_quantities


      subroutine write_quantity(ns, nCi, ital, istra, qu, lfirst, 
     .                          txttal, txtspc, txtunt)
      use b2mod_types

      implicit none
      integer, intent(in) :: ns, nCi, ital, istra
      real(R8), intent(in) :: qu(nCi,0:ns-1)
      logical, intent(inout) :: lfirst
      character(*), intent(in) :: txttal, txtspc(0:ns-1), txtunt
      integer :: iout, i, k
      CHARACTER(50) :: FNAME, FORMA, FORME, FORME2
#ifndef F2003
      integer newunit
      external newunit
#endif

      FORMA=REPEAT(' ',50)
      FORMA='(6X,   A25)'
      WRITE (FORMA(5:7),'(I3)') ns

      FORME=REPEAT(' ',50)
      FORME='(I10,   ES25.7)'
      WRITE (FORME(6:8),'(I3)') ns

      FORME2=REPEAT(' ',50)
      FORME2='(6X,   ES25.7)'
      WRITE (FORME2(5:7),'(I3)') ns
      
      FNAME = 'outtal_   '
      WRITE (FNAME(8:9),'(i0)') ITAL

      IF (LFIRST) THEN
#ifdef F2003
        OPEN (NEWUNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .        ACCESS='SEQUENTIAL')
#else
        IOUT = newunit()
        OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .        ACCESS='SEQUENTIAL')
#endif
      ELSE
#ifdef F2003
        OPEN (NEWUNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .        ACCESS='SEQUENTIAL',POSITION='APPEND')
#else
        IOUT = NEWUNIT()
        OPEN (UNIT=IOUT,FILE=FNAME,FORM='FORMATTED',
     .        ACCESS='SEQUENTIAL',POSITION='APPEND')
#endif
      END IF

      WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'
      WRITE (IOUT,'(A,I6)') 'ISTRA = ',ISTRA
      WRITE (IOUT,'(A)')
     .      '+++++++++++++++++++++++++++++++++++++++++++++++++'

      WRITE (IOUT,'(A)') TXTTAL
      WRITE (IOUT,'(A,I10)') 'NCELLS:   ',NCI
      WRITE (IOUT,'(A,I10)') 'NSPECIES: ',NS
      WRITE (IOUT,'(A)') 'SPECIES'
      WRITE (IOUT,FORMA) (TRIM(TXTSPC(K)), K=0,NS-1)
      WRITE (IOUT,'(A)') 'UNITS'
      WRITE (IOUT,FORMA) (TRIM(TXTUNT), K=0,NS-1)

      WRITE (IOUT,'(A)')
     .        'TOTAL ("UNITS*CM**3), AND MEAN VALUE ("UNITS") '
      WRITE (IOUT,'(A)') 'TOTAL'
      WRITE (IOUT,FORME2) (sum(qu(:,k)), K=0,ns-1)
      WRITE (IOUT,'(A)') 'MEAN'
      WRITE (IOUT,FORME2) (sum(qu(:,k))/nCi, K=0, NS-1)

      WRITE (IOUT,'(A)')
     .      '================================================='
      DO I=1, NCI
        IF (ANY(ABS(qu(I,0:ns-1)) > 1.e-30_R8))
     .      WRITE (IOUT,FORME) I,(qu(I,K), K=0, NS-1)
      END DO
      WRITE (IOUT,'(A)')
     .      '================================================='

      CLOSE (UNIT=IOUT)
C
      LFIRST = .FALSE.
      
      return
      end subroutine write_quantity

