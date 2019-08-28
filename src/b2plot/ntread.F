c
C RHW 18.08.92:  UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
C
*//NTREAD//
*=======================================================================
*                 S U B R O U T I N E  N T R E A D
*=======================================================================
*
*
*  1. PURPOSE
*
*     NTREAD READS IN THE ATOM AND MOLECULE DATA OF EIRENE FOR PLOT.
*
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE NTREAD
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2cmfs
      use b2mod_constants
      use b2mod_indirect
     ,  , only : triangleArea
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
      use eirmod_ctrig
      use eirmod_cestim
      use eirmod_comusr
     , , only: eirene_alloc_comusr
      use eirmod_coutau
     , , only: eirene_alloc_coutau
      use eirmod_braeir
     .   , only : pux, puy, pvx, pvy
#endif
      implicit none

*  -- LOCAL SCALARS --
#ifdef B25_EIRENE
      INTEGER is, i, j, n, io, idummy, jvft46
      real (kind=R8) :: width
      character*256 filename
      logical ex
#endif

*  -- PROCEDURES --
      external subini, subend
#ifdef B25_EIRENE
      external ipgetr, xertst
!      real (kind=R8) :: triangleArea
!      external triangleArea
*   ..initialisation
      data width /1.0_R8/
      save width
#endif

*=======================================================================
*//COMPUTATION//
C
      call subini ('ntread')

#ifdef B25_EIRENE
C
* ----------------------------------------------------------------------
*  READ EIRENE SPECIES DENSITIES AND TEMPERATURES
* ----------------------------------------------------------------------
*
*   ..read input data
      call ipgetr ('b2agmt_1d_width', width)
      call xertst (0.0_R8.lt.width, 'faulty input width')
*   ..test isymm
      call xertst (0.le.isymm.and.isymm.le.4,
     &  'faulty argument isymm')

#ifndef HPUX
      filename='fort.46'
#else
      filename='ftn46'
#endif
      call find_file(filename,ex)
      if (.not.ex) then
        call subend ()
        return
      end if

      OPEN (UNIT=33,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1983)
      OPEN (UNIT=34,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1984)
      OPEN (UNIT=35,ACCESS='SEQUENTIAL',FORM='FORMATTED',ERR=1985)
      REWIND 33
      REWIND 34
      REWIND 35
C
      READ(33,*,IOSTAT=IO,ERR=1983) NRKNOT  !VK IO
      IF(IO.NE.0) NRKNOT=0 !VK
      WRITE(*,'(a14,i6)') 'NRKNOT      = ',NRKNOT
C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
      READ(34,*,IOSTAT=IO,ERR=1984) NTRII   !VK IO
      IF(IO.NE.0) NTRII=0 !VK
      WRITE(*,'(a14,i6)') 'NTRII       = ',NTRII

      READ (35,*,IOSTAT=IO,ERR=1985) IDUMMY
      IF(IO.NE.0) IDUMMY=0 !VK
      IF (IDUMMY /= NTRII) THEN
        WRITE (*,*) ' NUMBER OF TRIANGLES DO NOT MATCH '
        WRITE (*,*) ' IN ELEMENTS AND NEIGHBOR FILES'
        WRITE (*,*) ' PLEASE CHECK THE GEOMETRY '
#ifndef HPUX
        call xerrab('Mismatch between fort.34 and fort.35 files!')
#else
        call xerrab('Mismatch between ftn34 and ftn35 files!')
#endif
      END IF

      CALL EIRENE_ALLOC_CTRIG
      CALL EIRENE_ALLOC_CESTIM(1)
      CALL EIRENE_ALLOC_COMUSR(1)
      CALL EIRENE_ALLOC_COUTAU
      CALL EIRENE_SETPRM
      CALL EIRENE_ALLOC_CESTIM(2)
      CALL EIRENE_ASSOCIATE_CESTIM
      call alloc_b2mod_b2plot_eirene(natmi,nmoli,nioni,ntrii,wklng)

      READ(33,*,IOSTAT=IO) (XTRIAN(I),I=1,NRKNOT)
      READ(33,*,IOSTAT=IO) (YTRIAN(I),I=1,NRKNOT)
      IF(IO.GT.0) THEN
C TO READ THE FILE PRODUCED BY DOS
#ifndef HPUX
        WRITE(*,*) "Cannot read fort.33 in regular way"
#else
        WRITE(*,*) "Cannot read ftn33 in regular way"
#endif
        WRITE(*,*) "Trying to read line by line..."
        REWIND 33
        READ(33,*) NRKNOT
        N=INT((NRKNOT+3)/4)
        DO I=1,N
          J=(I-1)*4+1
          READ(33,*,IOSTAT=IO) XTRIAN(J:J+4)
        END DO
        DO I=1,N
          J=(I-1)*4+1
          READ(33,*,IOSTAT=IO) YTRIAN(J:J+4)
        END DO
        WRITE(*,*) "... ready!"
      END IF

      IF(NTRII.GT.0) THEN !VK
        DO I=1,NTRII
          READ(34,*) J,NECKE(1,I),NECKE(2,I),NECKE(3,I)
        ENDDO

        DO I=1,NTRII
          READ(35,*) J,NCHBAR(1,I),NSEITE(1,I),IDUMMY,
     >                 NCHBAR(2,I),NSEITE(2,I),IDUMMY,
     >                 NCHBAR(3,I),NSEITE(3,I),IDUMMY,
     >                 IXTRI(I),IYTRI(I)
        ENDDO
      END IF !VK

      CLOSE (UNIT=33)
      CLOSE (UNIT=34)
      CLOSE (UNIT=35)
C
      do i=1,nrknot
        xtrian(i)=xtrian(i)/100.0_R8
        ytrian(i)=ytrian(i)/100.0_R8
      end do
      do i=1,ntrii
        triangle_area(i)=triangleArea(
     .   xtrian(necke(1,i)), ytrian(necke(1,i)),
     .   xtrian(necke(2,i)), ytrian(necke(2,i)),
     .   xtrian(necke(3,i)), ytrian(necke(3,i)))
      end do
      if (isymm.eq.0) then ! slab geometry
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*width
        end do
      else if (isymm.eq.1 .or. isymm.eq.2) then
*    ..case of rotational symmetry about crx=0 axis
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .     (xtrian(necke(1,i))+xtrian(necke(2,i))+xtrian(necke(3,i)))/
     .      3.0_R8
        end do
      else if (isymm.eq.3 .or. isymm.eq.4) then
*    ..case of rotational symmetry about cry=0 axis
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .     (ytrian(necke(1,i))+ytrian(necke(2,i))+ytrian(necke(3,i)))/
     .      3.0_R8
        end do
      end if
  
      rewind (46)
C
      READ(46,*,END=99,ERR=99) idummy,jvft46
      IF (IDUMMY /= NTRII) THEN
        WRITE (*,*) ' NUMBER OF TRIANGLES DO NOT MATCH '
        WRITE (*,*) ' IN GRID AND DATA FILES'
        WRITE (*,*) ' PLEASE CHECK THE GEOMETRY '
#ifndef HPUX
        call xerrab('Mismatch between fort.34 and fort.46 files!')
#else
        call xerrab('Mismatch between ftn34 and ftn46 files!')
#endif
      END IF
      READ(46,'(i4,2x,i4,2x,i4)',END=99,ERR=99) natmi,nmoli,nioni
#ifndef HPUX
      write(*,'(a,4i6)')
     . ' Dimensions in fort.46 = ',ntrii,natmi,nmoli,nioni
#else
      write(*,'(a,4i6)')
     . ' Dimensions in ftn46 = ',ntrii,natmi,nmoli,nioni
#endif
      call xertst(natmi .le. natm, 'NGREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NGREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NGREAD: NIONI > NION')
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ(46,'(A8)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ(46,'(A8)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
C
* test ions
C
      do 13 IS = 1, NIONI
        READ(46,'(A8)',END=99,ERR=99) TEXTIN(IS-1)
   13 CONTINUE
C
* data
* densities
C
      call read_title(46)
      read(46,'(6e15.7)') ((PDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((PDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((PDENI(i,j),j=1,ntrii),i=1,nioni)

* average energies
C
      call read_title(46)
      read(46,'(6e15.7)') ((EDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((EDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((EDENI(i,j),j=1,ntrii),i=1,nioni)

* X-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENI(i,j),j=1,ntrii),i=1,nioni)

* Y-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENI(i,j),j=1,ntrii),i=1,nioni)

* Z-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENI(i,j),j=1,ntrii),i=1,nioni)

* Volumes
C
      if (jvft46.ge.20170930) then
        call read_title(46)
        read(46,'(6e15.7)') (triangle_vol(j),j=1,ntrii)
        triangle_vol = triangle_vol/1.0e6_R8   ! convert to m^3
        if (.not.allocated(pux)) then
          allocate(pux(ntrii))
          allocate(puy(ntrii))
          allocate(pvx(ntrii))
          allocate(pvy(ntrii))
        end if
        call read_title(46)
        read(46,'(6e15.7)') (pux(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (puy(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (pvx(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (pvy(j),j=1,ntrii)
      end if
#endif

*  ---------------------------------------------------------------------
   99 CONTINUE
C
      call subend ()
      RETURN

*=======================================================================
*//END NTREAD//

#ifndef HPUX
1983  WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN FORT.33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.33 file!')
1984  WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN FORT.34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.34 file!')
1985  WRITE(*,*) "ERROR IN NTREAD: CANNOT OPEN FORT.35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.35 file!')
#else
1983  WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN FTN33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn33 file!')
1984  WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN FTN34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn34 file!')
1985  WRITE(*,*) "ERROR IN NTREAD: CANNOT OPEN FTN35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn35 file!')
#endif
      END



