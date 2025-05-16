c
C RHW 18.08.92: UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
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
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_wneutrals
      use eirmod_ctrig
      use eirmod_cestim
      use eirmod_cinit
      use eirmod_cgeom
     , , only : ncltal
      use eirmod_comusr
     , , only : eirene_alloc_comusr
      use eirmod_coutau
     , , only : eirene_alloc_coutau
      use eirmod_eirbra
     , , only : ncvtot
      use eirmod_braeir
     , , only : pux, puy, pvx, pvy
#endif
      use b2mod_subsys
      use b2us_data
      implicit none

*  -- LOCAL SCALARS --
#ifdef B25_EIRENE
      INTEGER is, i, j, iCv, icoadd, n, io, idummy
      character*256 filename
      logical ex
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
      filename=fort_lc//'46'
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
      WRITE(*,'(a14,i8)') 'NRKNOT      = ',NRKNOT
C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
      READ(34,*,IOSTAT=IO,ERR=1984) NTRII   !VK IO
      IF(IO.NE.0) NTRII=0 !VK
      WRITE(*,'(a14,i8)') 'NTRII       = ',NTRII

      READ (35,*,IOSTAT=IO,ERR=1985) IDUMMY
      IF(IO.NE.0) IDUMMY=0 !VK
      IF (IDUMMY /= NTRII) THEN
        WRITE (*,*) ' NUMBER OF TRIANGLES DO NOT MATCH '
        WRITE (*,*) ' IN ELEMENTS AND NEIGHBOR FILES'
        WRITE (*,*) ' PLEASE CHECK THE GEOMETRY '
        CALL XERRAB
     .   ('Mismatch between '//fort_lc//'34 and '//fort_lc//'35 files!')
      END IF

      CALL EIRENE_ALLOC_CTRIG
      CALL EIRENE_ALLOC_CESTIM(1)
      CALL EIRENE_ALLOC_COMUSR(1)
      CALL EIRENE_ALLOC_COUTAU
      CALL EIRENE_SETPRM
      CALL EIRENE_ALLOC_CESTIM(2)
      CALL EIRENE_ASSOCIATE_CESTIM
      call alloc_b2mod_b2plot_eirene(natmi,nmoli,nioni,ntrii,wklng)

      IF (NTRFRM == 0) THEN
        READ(33,*,IOSTAT=IO) (XTRIAN(I),I=1,NRKNOT)
        READ(33,*,IOSTAT=IO) (YTRIAN(I),I=1,NRKNOT)
        IF(IO.GT.0) THEN
C TO READ THE FILE PRODUCED BY DOS
          WRITE(*,*) "Cannot read "//fort_lc//"33 in regular way"
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
      ELSE
        DO I=1,NRKNOT
          READ(33,*) J,XTRIAN(I),YTRIAN(I)
        END DO
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
      call compute_triangle_area
      call compute_triangle_vol
C
      icoadd = maxval(ixtri)
      if (.not.allocated(ncltal)) allocate(ncltal(ntrii))
      do i = 1, ntrii
        iCv = ixtri(i)
        if (iCv.gt.0) then
          ncltal(i) = iCv
        else
          icoadd = icoadd + 1
          ncltal(i) = icoadd
        end if
      end do
C
      rewind (46)
C
      READ(46,*,END=99,ERR=99) ncvtot,jvft46
      READ(46,'(i4,2x,i4,2x,i4)',END=99,ERR=99) natmi,nmoli,nioni
      write(*,'(3a,4i6)')
     . ' Dimensions in ',fort_lc,'46 = ',ncvtot,natmi,nmoli,nioni
      IF (ICOADD /= NCVTOT) THEN
        WRITE (*,*) ' NUMBER OF TRIANGLES DO NOT MATCH'
        WRITE (*,*) ' BETWEEN GEOMETRY AND DATA FILES'
        WRITE (*,*) ' PLEASE CHECK THE GEOMETRY'
        CALL XERRAB
     .   ('Mismatch between '//fort_lc//'35 and '//fort_lc//'46 files!')
      END IF
      call xertst(natmi .le. natm, 'NTREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NTREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NTREAD: NIONI > NION')
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ(46,'(A)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ(46,'(A)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
C
* test ions
C
      do 13 IS = 1, NIONI
        READ(46,'(A)',END=99,ERR=99) TEXTIN(IS-1)
   13 CONTINUE
C
* data
* densities
C
      call read_title(46)
      read(46,'(6e15.7)') ((PDENA(i,j),j=1,ncvtot),i=1,max(natmi,1))
      call eirene_extend(pdena,ntrii,ncltal,max(natmi,1))
      call read_title(46)
      read(46,'(6e15.7)') ((PDENM(i,j),j=1,ncvtot),i=1,max(nmoli,1))
      call eirene_extend(pdenm,ntrii,ncltal,max(nmoli,1))
      call read_title(46)
      read(46,'(6e15.7)') ((PDENI(i,j),j=1,ncvtot),i=1,max(nioni,1))
      call eirene_extend(pdeni,ntrii,ncltal,max(nioni,1))

* average energies
C
      call read_title(46)
      read(46,'(6e15.7)') ((EDENA(i,j),j=1,ncvtot),i=1,max(natmi,1))
      call eirene_extend(edena,ntrii,ncltal,max(natmi,1))
      call read_title(46)
      read(46,'(6e15.7)') ((EDENM(i,j),j=1,ncvtot),i=1,max(nmoli,1))
      call eirene_extend(edenm,ntrii,ncltal,max(nmoli,1))
      call read_title(46)
      read(46,'(6e15.7)') ((EDENI(i,j),j=1,ncvtot),i=1,max(nioni,1))
      call eirene_extend(edeni,ntrii,ncltal,max(nioni,1))

* X-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENA(i,j),j=1,ncvtot),i=1,max(natmi,1))
      call eirene_extend(vxdena,ntrii,ncltal,max(natmi,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENM(i,j),j=1,ncvtot),i=1,max(nmoli,1))
      call eirene_extend(vxdenm,ntrii,ncltal,max(nmoli,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENI(i,j),j=1,ncvtot),i=1,max(nioni,1))
      call eirene_extend(vxdeni,ntrii,ncltal,max(nioni,1))

* Y-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENA(i,j),j=1,ncvtot),i=1,max(natmi,1))
      call eirene_extend(vydena,ntrii,ncltal,max(natmi,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENM(i,j),j=1,ncvtot),i=1,max(nmoli,1))
      call eirene_extend(vydenm,ntrii,ncltal,max(nmoli,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENI(i,j),j=1,ncvtot),i=1,max(nioni,1))
      call eirene_extend(vydeni,ntrii,ncltal,max(nioni,1))

* Z-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENA(i,j),j=1,ncvtot),i=1,max(natmi,1))
      call eirene_extend(vzdena,ntrii,ncltal,max(natmi,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENM(i,j),j=1,ncvtot),i=1,max(nmoli,1))
      call eirene_extend(vzdenm,ntrii,ncltal,max(nmoli,1))
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENI(i,j),j=1,ncvtot),i=1,max(nioni,1))
      call eirene_extend(vzdeni,ntrii,ncltal,max(nioni,1))

* Volumes
C
      if (jvft46.ge.20170930) then
        call read_title(46)
        read(46,'(6e15.7)') (triangle_vol(j),j=1,ntrii)
        triangle_vol = triangle_vol/1.0e6_R8   ! convert to m^3
        if (.not.allocated(pux)) then
          allocate(pux(mpg%nCv))
          allocate(puy(mpg%nCv))
          allocate(pvx(mpg%nCv))
          allocate(pvy(mpg%nCv))
        end if
        call read_title(46)
        read(46,'(6e15.7)') (pux(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (puy(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (pvx(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (pvy(j),j=1,mpg%nCv)
      end if
   99 CONTINUE
#endif
      close (46)
*  ---------------------------------------------------------------------
C
      call subend ()
      RETURN

*=======================================================================
*//END NTREAD//

#ifdef B25_EIRENE
 1983 WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN "//FORT//"33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'33 file!')
 1984 WRITE(*,*)
     w     "ERROR IN NTREAD: CANNOT OPEN "//FORT//"34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'34 file!')
 1985 WRITE(*,*) "ERROR IN NTREAD: CANNOT OPEN "//FORT//"35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'35 file!')
#endif
      END SUBROUTINE NTREAD

#ifdef B25_EIRENE
      subroutine compute_triangle_area
      use b2mod_b2plot
      use b2mod_math
      use eirmod_ctrig
      implicit none
      integer i

      do i=1,ntrii
        triangle_area(i)=triangleArea(
     .   xtrian(necke(1,i)), ytrian(necke(1,i)),
     .   xtrian(necke(2,i)), ytrian(necke(2,i)),
     .   xtrian(necke(3,i)), ytrian(necke(3,i)))
      end do
      return
      end subroutine compute_triangle_area

      subroutine compute_triangle_vol
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_b2plot
      use b2mod_constants
      use eirmod_ctrig
      implicit none
      integer i
      real (kind=R8) :: width
*   ..initialisation
      data width /1.0_R8/
      save width
*   ..procedures
      external ipgetr, xertst

*   ..read input data
      call ipgetr ('b2agmt_1d_width', width)
      call xertst (0.0_R8.lt.width, 'faulty input width')
*   ..test isymm
      call xertst (0.le.isymm.and.isymm.le.4,
     &  'faulty argument isymm')

      if (isymm.eq.0) then ! slab geometry
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*width
        end do
      else if (isymm.eq.1 .or. isymm.eq.2) then
*    ..case of rotational symmetry about r=0 axis
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .     (xtrian(necke(1,i))+xtrian(necke(2,i))+xtrian(necke(3,i)))/
     .      3.0_R8
        end do
      else if (isymm.eq.3 .or. isymm.eq.4) then
*    ..case of rotational symmetry about z=0 axis
        do i=1,ntrii
          triangle_vol(i)=triangle_area(i)*2.0_R8*pi*
     .     (ytrian(necke(1,i))+ytrian(necke(2,i))+ytrian(necke(3,i)))/
     .      3.0_R8
        end do
      end if
      return
      end subroutine compute_triangle_vol
#endif
