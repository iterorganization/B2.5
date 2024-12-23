c
C RHW 18.08.92: UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
C
*//B2YT_NTREAD//
*=======================================================================
*                 S U B R O U T I N E  B 2 Y T _ N T R E A D
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

      subroutine b2yt_ntread
#ifdef B25_EIRENE
     & (mpg, natm1, iatm1, nmol1, imol1, nion1, iion1)
#endif
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_math
      use b2mod_b2plot
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
#ifdef B25_EIRENE
      use eirmod_parmmod
     ,   , only : nion, natm, nmol
      use eirmod_wneutrals
      use eirmod_ctrig
      use eirmod_cinit
      use eirmod_cestim
      use eirmod_eirbra
     ,   , only : ncvtot
      use eirmod_braeir
     ,   , only : pux, puy, pvx, pvy
#endif
      use b2mod_subsys
      use b2us_map
      implicit none

#ifdef B25_EIRENE
*  -- INPUT ARGUMENTS --
      type(mapping), intent(in) :: mpg
      integer natm1, nmol1, nion1
      integer iatm1(natm), imol1(nmol), iion1(nion)

*  -- LOCAL SCALARS --
      integer is, i, j, n, jvftnew, io, idummy
      real (kind=R8) :: width
      real (kind=R8), allocatable :: field1(:,:)
      character*16, allocatable :: textan1(:), textmn1(:), textin1(:)
      character*256 filename
      data jvftnew / 20170930 /

*   ..initialisation
      data width /1.0_R8/
      save width

*  -- PROCEDURES --
      character*32 :: get_Eir_hash
      external ipgeti, ipgetr, xertst, get_Eir_hash
#endif

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
      WRITE(*,'(a14,i8)') 'NRKNOT      = ',NRKNOT
C
C     READ IN THE NUMBER OF TRIANGLES AND ATTRIBUTES OF THE TRIANGLES
      READ(34,*,IOSTAT=IO,ERR=1984) NTRII   !VK IO
      IF(IO.NE.0) NTRII=0 !VK
      WRITE(*,'(a14,i8)') 'NTRII       = ',NTRII

      READ (35,*,IOSTAT=IO,ERR=1985) IDUMMY
      IF(IO.NE.0) IDUMMY=0 !VK
      IF (IDUMMY /= NTRII) THEN
        WRITE (*,*) ' NUMBER OF TRIANGLES DO NOT MATCH'
        WRITE (*,*) ' IN ELEMENTS AND NEIGHBOR FILES'
        WRITE (*,*) ' PLEASE CHECK THE GEOMETRY'
        call xerrab
     .   ('Mismatch between '//fort_lc//'34 and '//fort_lc//'35 files!')
      END IF

      CALL EIRENE_ALLOC_CTRIG
      CALL EIRENE_ALLOC_CESTIM(1)
      CALL EIRENE_ALLOC_CESTIM(2)
      CALL EIRENE_ASSOCIATE_CESTIM

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
      rewind (46)
      filename=fort_lc//'46.2'
      open(unit=47,file=trim(filename),status='unknown')
      write(*,*) 'Writing '//trim(filename)
C
      READ(46,*,END=99,ERR=99) ncvtot,jvft46
      WRITE(47,'(I6,2X,i8,2x,a32)') ncvtot,jvftnew,get_Eir_hash()
      READ (46,'(i4,2x,i4,2x,i4)',END=99,ERR=99) natmi,nmoli,nioni
      WRITE(47,'(I4,2X,I4,2X,i4)')               natm1,nmol1,nion1
      write(*,'(3a,4i6)')
     . ' Dimensions in ',fort_lc,'46 =   ',ntrii,natmi,nmoli,nioni
      write(*,'(3a,4i6)')
     . ' Dimensions in ',fort_lc,'46.2 = ',ntrii,natm1,nmol1,nion1
      call xertst(natmi .le. natm, 'NTREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NTREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NTREAD: NIONI > NION')
      allocate(field1(max(1,natm1,nmol1,nion1),max(ntrii,ncvtot)))

c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ (46,'(A)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
      if (natm1.gt.0) then
        allocate(textan1(0:natm1-1))
        do is = 1, natmi
          if (iatm1(is).ne.0) TEXTAN1(iatm1(is)-1) = TEXTAN(is-1)
        end do
        do is = 1, natm1
          WRITE(47,'(A)') TEXTAN1(IS-1)
        end do
        deallocate(textan1)
      end if
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ (46,'(A)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
      if (nmol1.gt.0) then
        allocate(textmn1(0:nmol1-1))
        do is = 1, nmoli
          if (imol1(is).ne.0) TEXTMN1(imol1(is)-1) = TEXTMN(is-1)
        end do
        do is = 1, nmol1
          WRITE(47,'(A)') TEXTMN1(IS-1)
        end do
        deallocate(textmn1)
      end if
C
* test ions
C
      do 13 IS = 1, NIONI
        READ (46,'(A)',END=99,ERR=99) TEXTIN(IS-1)
   13 CONTINUE
      if (nion1.gt.0) then
        allocate(textin1(0:nion1-1))
        do is = 1, nioni
          if (iion1(is).ne.0) TEXTIN1(iion1(is)-1) = TEXTIN(is-1)
        end do
        do is = 1, nion1
          WRITE(47,'(A)') TEXTIN1(IS-1)
        end do
        deallocate(textin1)
      end if
C
* data
* densities
C
      call read_title(46)
      read(46,'(6e15.7)') ((PDENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      call read_title(46)
      read(46,'(6e15.7)') ((PDENM(i,j),j=1,ncvtot),i=1,max(1,nmoli))
      call read_title(46)
      read(46,'(6e15.7)') ((PDENI(i,j),j=1,ncvtot),i=1,max(1,nioni))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = pdena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'pdena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'pdena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      if (nmol1.gt.0) then
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field1(imol1(is),1:ncvtot) = pdenm(is,1:ncvtot)
          end if
        end do
        call write_title(47,'pdenm',ncvtot*nmol1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'pdenm (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nmol1))
      if (nion1.gt.0) then
        do is = 1, nioni
          if (iion1(is).ne.0) then
            field1(iion1(is),1:ncvtot) = pdeni(is,1:ncvtot)
          end if
        end do
        call write_title(47,'pdeni',ncvtot*nion1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'pdeni (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nion1))

* average energies
C
      call read_title(46)
      read(46,'(6e15.7)') ((EDENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      call read_title(46)
      read(46,'(6e15.7)') ((EDENM(i,j),j=1,ncvtot),i=1,max(1,nmoli))
      call read_title(46)
      read(46,'(6e15.7)') ((EDENI(i,j),j=1,ncvtot),i=1,max(1,nioni))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = edena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'edena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'edena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      if (nmol1.gt.0) then
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field1(imol1(is),1:ncvtot) = edenm(is,1:ncvtot)
          end if
        end do
        call write_title(47,'edenm',ncvtot*nmol1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'edenm (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nmol1))
      if (nion1.gt.0) then
        do is = 1, nioni
          if (iion1(is).ne.0) then
            field1(iion1(is),1:ncvtot) = edeni(is,1:ncvtot)
          end if
        end do
        call write_title(47,'edeni',ncvtot*nion1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'edeni (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nion1))

* X-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENM(i,j),j=1,ncvtot),i=1,max(1,nmoli))
      call read_title(46)
      read(46,'(6e15.7)') ((VXDENI(i,j),j=1,ncvtot),i=1,max(1,nioni))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = vxdena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vxdena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vxdena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      if (nmol1.gt.0) then
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field1(imol1(is),1:ncvtot) = vxdenm(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vxdenm',ncvtot*nmol1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vxdenm (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nmol1))
      if (nion1.gt.0) then
        do is = 1, nioni
          if (iion1(is).ne.0) then
            field1(iion1(is),1:ncvtot) = vxdeni(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vxdeni',ncvtot*nion1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vxdeni (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nion1))

* Y-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENM(i,j),j=1,ncvtot),i=1,max(1,nmoli))
      call read_title(46)
      read(46,'(6e15.7)') ((VYDENI(i,j),j=1,ncvtot),i=1,max(1,nioni))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = vydena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vydena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vydena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      if (nmol1.gt.0) then
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field1(imol1(is),1:ncvtot) = vydenm(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vydenm',ncvtot*nmol1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vydenm (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nmol1))
      if (nion1.gt.0) then
        do is = 1, nioni
          if (iion1(is).ne.0) then
            field1(iion1(is),1:ncvtot) = vydeni(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vydeni',ncvtot*nion1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vydeni (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nion1))

* Z-directed velocities
C
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENM(i,j),j=1,ncvtot),i=1,max(1,nmoli))
      call read_title(46)
      read(46,'(6e15.7)') ((VZDENI(i,j),j=1,ncvtot),i=1,max(1,nioni))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = vzdena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vzdena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vzdena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      if (nmol1.gt.0) then
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field1(imol1(is),1:ncvtot) = vzdenm(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vzdenm',ncvtot*nmol1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vzdenm (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nmol1))
      if (nion1.gt.0) then
        do is = 1, nioni
          if (iion1(is).ne.0) then
            field1(iion1(is),1:ncvtot) = vzdeni(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vzdeni',ncvtot*nion1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vzdeni (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,nion1))

* Triangle volumes
C
      if(jvft46.lt.20231224) then
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
      if (.not.allocated(pux)) then
        allocate(pux(mpg%nCv))
        allocate(puy(mpg%nCv))
        allocate(pvx(mpg%nCv))
        allocate(pvy(mpg%nCv))
      end if
      if(jvft46.ge.20170930) then
        call read_title(46)
        read(46,'(6e15.7)') (PUX(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (PUY(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (PVX(j),j=1,mpg%nCv)
        call read_title(46)
        read(46,'(6e15.7)') (PVY(j),j=1,mpg%nCv)
      else
        pux = 0.0_R8
        puy = 0.0_R8
        pvx = 0.0_R8
        pvy = 0.0_R8
      end if
      call write_title(47,'pux',mpg%nCv)
      write(47,'(6e15.7)') (PUX(j),j=1,mpg%nCv)
      call write_title(47,'puy',mpg%nCv)
      write(47,'(6e15.7)') (PUY(j),j=1,mpg%nCv)
      call write_title(47,'pvx',mpg%nCv)
      write(47,'(6e15.7)') (PVX(j),j=1,mpg%nCv)
      call write_title(47,'pvy',mpg%nCv)
      write(47,'(6e15.7)') (PVY(j),j=1,mpg%nCv)

      deallocate(field1)
      deallocate(pux,puy,pvx,pvy)

* Atom condensation
      call read_title(46)
      read(46,'(6e15.7)') ((PGENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = pgena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'pgena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'pgena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      read(46,'(6e15.7)') ((VGENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = vgena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'vgena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'vgena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))
      read(46,'(6e15.7)') ((EGENA(i,j),j=1,ncvtot),i=1,max(1,natmi))
      if (natm1.gt.0) then
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            field1(iatm1(is),1:ncvtot) = egena(is,1:ncvtot)
          end if
        end do
        call write_title(47,'egena',ncvtot*natm1)
      else
        field1(1,1:ncvtot) = 0.0_R8
        call write_title(47,'egena (not computed)',ncvtot)
      end if
      write(47,'(6e15.7)') ((FIELD1(i,j),j=1,ncvtot),i=1,max(1,natm1))

   99 CONTINUE
#endif

*  ---------------------------------------------------------------------
      close(46)
      close(47)
C
      call subend ()
      return

*=======================================================================
*//END B2YT_NTREAD//
#ifdef B25_EIRENE
 1983 WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN "//FORT//"33 ",
     .     "(NODES OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'33 file!')
 1984 WRITE(*,*)
     w     "ERROR IN B2YT_NTREAD: CANNOT OPEN "//FORT//"34 ",
     .     "(TABLE OF CELLS OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'34 file!')
 1985 WRITE(*,*) "ERROR IN B2YT_NTREAD: CANNOT OPEN "//FORT//"35 ",
     .           "(TABLE OF NEIGHBORS OF TRIANGULAR GRID)"
      call xerrab('Cannot open '//fort_lc//'35 file!')
#endif
      end subroutine b2yt_ntread
