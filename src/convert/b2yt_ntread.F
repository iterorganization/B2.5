c
C RHW 18.08.92:  UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
C
*//B2YT_NGREAD//
*=======================================================================
*                 S U B R O U T I N E  B 2 Y T _ N G R E A D
*=======================================================================
*
*
*  1. PURPOSE
*
*     B2YT_NTREAD READS IN AND CONVERTS THE ATOM AND MOLECULE DATA OF 
*     EIRENE FOR B2YT.
*
*
*=======================================================================
*//DECLARATIONS//

      subroutine b2yt_ntread(natm1, iatm1, nmol1, imol1)
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_ctrig
      use eirmod_cestim
#endif
      implicit none

*  -- INPUT ARGUMENTS --
      integer natm1, nmol1
      integer iatm1(natm), imol1(nmol)

*  -- LOCAL SCALARS --
#ifdef B25_EIRENE
      integer jvft46
      integer is, i, j, n, jvftnew, nstra1, io, idummy
      real (kind=R8) :: width
      real (kind=R8), allocatable :: field1(:,:)
      character*16, allocatable :: textan1(:), textmn1(:)
      data jvftnew / 20170930 /

*   ..initialisation
      data width /1.0_R8/
      save width

*  -- PROCEDURES --
      character*32 :: get_B25_hash
      logical streql
      external ipgeti, ipgetr, xertst, streql, get_B25_hash
#endif
      external subini, subend

*=======================================================================
*//COMPUTATION//
C
      call subini ('b2yt_ntread')

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
      CALL EIRENE_ALLOC_CESTIM(2)
      CALL EIRENE_ASSOCIATE_CESTIM

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
      rewind (46)
#ifndef HPUX
      open(unit=47,file='fort.46.2',status='unknown')
      write(*,*) 'Writing fort.46.2'
#else
      open(unit=47,file='ftn46.2',status='unknown')
      write(*,*) 'Writing ftn46.2'
#endif
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
      WRITE(47,'(I6,2X,i8,2x,a32)') ntrii,jvftnew,get_B25_hash()
      READ(46,'(i4,2x,i4,2x,i4)',END=99,ERR=99) natmi,nmoli,nioni
      WRITE(47,'(I4,2X,I4,2X,i4)')              natm1,nmol1,nioni
#ifndef HPUX
      write(*,'(a,4i6)')
     . ' Dimensions in fort.46 =   ',ntrii,natmi,nmoli,nioni
      write(*,'(a,4i6)')
     . ' Dimensions in fort.46.2 = ',ntrii,natm1,nmol1,nioni
#else
      write(*,'(a,4i6)')
     . ' Dimensions in ftn46 =   ',ntrii,natmi,nmoli,nioni
      write(*,'(a,4i6)')
     . ' Dimensions in ftn46.2 = ',ntrii,natm1,nmol1,nioni
#endif
      call xertst(natmi .le. natm, 'NGREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NGREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NGREAD: NIONI > NION')
      allocate(field1(ntrii,max(natm1,nmol1)))

c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ (46,'(A8)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
      allocate(textan1(0:natm1-1))
      do is = 1, natmi
        if (iatm1(is).ne.0) TEXTAN1(iatm1(is)-1) = TEXTAN(is-1)
      end do
      do is = 1, natm1
        WRITE(47,'(A8)') TEXTAN1(IS-1)
      end do
      deallocate(textan1)
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ (46,'(A8)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
      allocate(textmn1(0:nmol1-1))
      do is = 1, nmoli
        if (imol1(is).ne.0) TEXTMN1(imol1(is)-1) = TEXTMN(is-1)
      end do
      do is = 1, nmol1
        WRITE(47,'(A8)') TEXTMN1(IS-1)
      end do
      deallocate(textmn1)
C
* test ions
C
      do 13 IS = 1, NIONI
        READ (46,'(A8)',END=99,ERR=99) TEXTIN(IS-1)
        WRITE(47,'(A8)')               TEXTIN(IS-1)
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
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field1(1:ntrii,iatm1(is)) = pdena(1:ntrii,is)
        end if
      end do
      call write_title(47,'pdena',ntrii*natm1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,natm1)
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field1(1:ntrii,imol1(is)) = pdenm(1:ntrii,is)
        end if
      end do
      call write_title(47,'pdenm',ntrii*nmol1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,nmol1)
      call write_title(47,'pdeni',ntrii*nioni)
      write(47,'(6e15.7)') ((PDENI(i,j),j=1,ntrii),i=1,nioni)

* average energies
C
      call read_title(46)
      read(46,'(6e15.7)') ((EDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((EDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((EDENI(i,j),j=1,ntrii),i=1,nioni)
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field1(1:ntrii,iatm1(is)) = edena(1:ntrii,is)
        end if
      end do
      call write_title(47,'edena',ntrii*natm1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,natm1)
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field1(1:ntrii,imol1(is)) = edenm(1:ntrii,is)
        end if
      end do
      call write_title(47,'edenm',ntrii*nmol1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,nmol1)
      call write_title(47,'edeni',ntrii*nioni)
      write(47,'(6e15.7)') ((EDENI(i,j),j=1,ntrii),i=1,nioni)

* X-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENI(i,j),j=1,ntrii),i=1,nioni)
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field1(1:ntrii,iatm1(is)) = vxdena(1:ntrii,is)
        end if
      end do
      call write_title(47,'vxdena',ntrii*natm1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,natm1)
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field1(1:ntrii,imol1(is)) = vxdenm(1:ntrii,is)
        end if
      end do
      call write_title(47,'vxdenm',ntrii*nmol1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,nmol1)
      call write_title(47,'vxdeni',ntrii*nioni)
      write(47,'(6e15.7)') ((VXDENI(i,j),j=1,ntrii),i=1,nioni)

* Y-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENI(i,j),j=1,ntrii),i=1,nioni)
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field1(1:ntrii,iatm1(is)) = vydena(1:ntrii,is)
        end if
      end do
      call write_title(47,'vydena',ntrii*natm1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,natm1)
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field1(1:ntrii,imol1(is)) = vydenm(1:ntrii,is)
        end if
      end do
      call write_title(47,'vydenm',ntrii*nmol1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,nmol1)
      call write_title(47,'vydeni',ntrii*nioni)
      write(47,'(6e15.7)') ((VYDENI(i,j),j=1,ntrii),i=1,nioni)

* Z-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENA(i,j),j=1,ntrii),i=1,natmi)
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENM(i,j),j=1,ntrii),i=1,nmoli)
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENI(i,j),j=1,ntrii),i=1,nioni)
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field1(1:ntrii,iatm1(is)) = vzdena(1:ntrii,is)
        end if
      end do
      call write_title(47,'vzdena',ntrii*natm1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,natm1)
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field1(1:ntrii,imol1(is)) = vzdenm(1:ntrii,is)
        end if
      end do
      call write_title(47,'vzdenm',ntrii*nmol1)
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ntrii),i=1,nmol1)
      call write_title(47,'vzdeni',ntrii*nioni)
      write(47,'(6e15.7)') ((VZDENI(i,j),j=1,ntrii),i=1,nioni)

* Triangle volumes
C
      if(jvft46.lt.20170930) then
        do i=1,ntrii
          triangle_area(i)=triangleArea(
     .     xtrian(necke(1,i)), ytrian(necke(1,i)),
     .     xtrian(necke(2,i)), ytrian(necke(2,i)),
     .     xtrian(necke(3,i)), ytrian(necke(3,i)))
        end do
        if (isymm.eq.0) then ! slab geometry
          do i=1,ntrii
            triangle_vol(i)=triangle_area(i)*width
          end do
        else if (isymm.eq.1 .or. isymm.eq.2) then
*    ..case of rotational symmetry about crx=0 axis
          do i=1,ntrii
            triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .       (xtrian(necke(1,i))+xtrian(necke(2,i))+xtrian(necke(3,i)))/
     .        3.0_R8
          end do
        else if (isymm.eq.3 .or. isymm.eq.4) then
*    ..case of rotational symmetry about cry=0 axis
          do i=1,ntrii
            triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .       (ytrian(necke(1,i))+ytrian(necke(2,i))+ytrian(necke(3,i)))/
     .        3.0_R8
          end do
        end if
      else
        call read_title(46)
        read(46,'(6e15.7)') (TRIANGLE_VOL(j),j=1,ntrii)
      end if
      call write_title(47,'volumes',ntrii)
      write(47,'(6e15.7)') (TRIANGLE_VOL(j),j=1,ntrii)

* Direction normals
C
      if(jvft46.ge.20170930) then
        call read_title(46)
        read(46,'(6e15.7)') (PUX(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (PUY(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (PVX(j),j=1,ntrii)
        call read_title(46)
        read(46,'(6e15.7)') (PVY(j),j=1,ntrii)
      else
        pux = 0.0_R8
        puy = 0.0_R8
        pvx = 0.0_R8
        pvy = 0.0_R8
      end if
      call write_title(47,'pux',ntrii)
      write(47,'(6e15.7)') (PUX(j),j=1,ntrii)
      call write_title(47,'puy',ntrii)
      write(47,'(6e15.7)') (PUY(j),j=1,ntrii)
      call write_title(47,'pvx',ntrii)
      write(47,'(6e15.7)') (PVX(j),j=1,ntrii)
      call write_title(47,'pvy',ntrii)
      write(47,'(6e15.7)') (PVY(j),j=1,ntrii)

      deallocate(field1)
#endif

*  ---------------------------------------------------------------------
   99 CONTINUE
      close(46)
      close(47)
C
      call subend ()
      RETURN

*=======================================================================
*//END B2YT_NTREAD//

#ifndef HPUX
1983  WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN FORT.33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.33 file!')
1984  WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN FORT.34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.34 file!')
1985  WRITE(*,*) "ERROR IN B2YT_NTREAD: CANNOT OPEN FORT.35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open fort.35 file!')
#else
1983  WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN FTN33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn33 file!')
1984  WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN FTN34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn34 file!')
1985  WRITE(*,*) "ERROR IN B2YT_NTREAD: CANNOT OPEN FTN35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open ftn35 file!')
#endif
      END



