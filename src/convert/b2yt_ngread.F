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
*     B2YT_NGREAD READS IN AND CONVERTS THE ATOM AND MOLECULE DATA OF 
*     EIRENE FOR PLOT.
*
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE b2yt_NGREAD(solps_version,eirene_format)
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2cmpa
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_clogau
#endif
      implicit none

*  -- INPUT ARGUMENTS --
      character*4 solps_version
      character*8 eirene_format

*  -- LOCAL SCALARS --
      INTEGER NNX, NNY, IS, IX, IY, jvft44, istra, nplsi
      integer i, j, k, ncl, jvftnew

*  -- PROCEDURES --
      logical streql
      EXTERNAL GFSUB2
      external subini, subend, ipgeti, streql

*=======================================================================
*//COMPUTATION//
C
      call subini ('b2yt_ngread')

      nplsi = 0
      do is = 0, ns-1
        if (.not.is_neutral(is)) nplsi = nplsi + 1
      end do

#ifdef B25_EIRENE
C
* ----------------------------------------------------------------------
*  READ EIRENE SPECIES DENSITIES AND TEMPERATURES
* ----------------------------------------------------------------------
*

      rewind (44)
      open(unit=45,file='fort.44.2',status='unknown')
C
      READ (44,'(I4,2X,I4,2x,i8)',END=99,ERR=99) nnx,nny,jvft44
      if((streql(solps_version,'5.0').or.streql(solps_version,'5.2'))
     &   .and.
     &   (streql(eirene_format,'old').or.streql(eirene_format,'new')))
     & then
        jvftnew = min(jvft44,961228)
      else if (streql(solps_version,'4.3')) then
        jvftnew = min(jvft44,20080706)
      else
        jvftnew = jvft44
      endif
      WRITE(45,'(I4,2X,I4,2x,i8)')               nnx,nny,jvftnew
      READ (44,'(I4,2X,I4,2X,I4)',END=99,ERR=99) NATMI,NMOLI,NIONI
      WRITE(45,'(I4,2X,I4,2X,i4)')               NATMI,NMOLI,NIONI
      write(*,'(a,5i4)') 
     . ' Dimensions in fort.44 = ',nnx,nny,natmi,nmoli,nioni
      write(*,'(a,2i4)') ' Array inner dimensions = ', ndx , ndy
      if(natmi .gt. natm) stop 'NGREAD: NATMI > NATM'
      if(nmoli .gt. nmol) stop 'NGREAD: NMOLI > NMOL'
      if(nioni .gt. nion) stop 'NGREAD: NIONI > NION'
      if(nnx   .gt. ndx ) stop 'NGREAD:   NNX > NDX '
      if(nny   .gt. ndy ) stop 'NGREAD:   NNY > NDY '
      call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstra,natmi)
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ (44,'(A8)',END=99,ERR=99) TEXTAN(IS-1)
        WRITE(45,'(A8)')               TEXTAN(IS-1)
   11 CONTINUE
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ (44,'(A8)',END=99,ERR=99) TEXTMN(IS-1)
        WRITE(45,'(A8)')               TEXTMN(IS-1)
   12 CONTINUE
C
* test ions
C
      do 13 IS = 1, NIONI
        READ (44,'(A8)',END=99,ERR=99) TEXTIN(IS-1)
        WRITE(45,'(A8)')               TEXTIN(IS-1)
   13 CONTINUE
C
* data
* neutral atoms
C
      do 14 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,dab2(:,:,is,1))
   14 CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,dab2,ndx,ndy,natm,1,1)

      do 24 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,tab2(:,:,is,1))
   24 CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,tab2,ndx,ndy,natm,1,1)

* neutral molecules
C
      do 15 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,dmb2(:,:,is,1))
   15 CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,dmb2,ndx,ndy,nmol,1,1)

      do 25 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,tmb2(:,:,is,1))
   25 CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,tmb2,ndx,ndy,nmol,1,1)

* test ions
C
      do 16 IS = 1, NIONI
        call gfsub2(44,ndx,ndy,nnx,nny,dib2(:,:,is,1))
   16 CONTINUE
      call eirene_neutr(45,nnx,nny,nioni,dib2,ndx,ndy,nion,1,1)

      do 26 IS = 1, NIONI
        call gfsub2(44,ndx,ndy,nnx,nny,tib2(:,:,is,1))
   26 CONTINUE
      call eirene_neutr(45,nnx,nny,nioni,tib2,ndx,ndy,nion,1,1)

* neutral atom radial particle flux
C
      do 114 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,rfluxa(:,:,is,1))
 114  CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,rfluxa,ndx,ndy,natm,1,1)

* neutral molecules radial particle flux
C
      do 115 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,rfluxm(:,:,is,1))
 115  CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,rfluxm,ndx,ndy,nmol,1,1)

* neutral atom poloidal particle flux
C
      do 214 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,pfluxa(:,:,is,1))
 214  CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,pfluxa,ndx,ndy,natm,1,1)

* neutral molecules poloidal particle flux
C
      do 215 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,pfluxm(:,:,is,1))
 215  CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,pfluxm,ndx,ndy,nmol,1,1)

* neutral atom radial energy flux
C
      do 314 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,refluxa(:,:,is,1))
 314  CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,refluxa,ndx,ndy,natm,1,1)

* neutral molecules radial energy flux
C
      do 415 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,refluxm(:,:,is,1))
 415  CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,refluxm,ndx,ndy,nmol,1,1)

* neutral atom poloidal energy flux
C
      do 514 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,pefluxa(:,:,is,1))
 514  CONTINUE
      call eirene_neutr(45,nnx,nny,natmi,pefluxa,ndx,ndy,natm,1,1)

* neutral molecules poloidal energy flux
C
      do 515 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,pefluxm(:,:,is,1))
 515  CONTINUE
      call eirene_neutr(45,nnx,nny,nmoli,pefluxm,ndx,ndy,nmol,1,1)

      call gfsub2(44,ndx,ndy,nnx,nny,emiss(:,:,1,1))
      call gfsub2(44,ndx,ndy,nnx,nny,emissmol(:,:,1,1))
      call eirene_neutr(45,nnx,nny,1,emiss,ndx,ndy,1,1,1)
      call eirene_neutr(45,nnx,nny,1,emissmol,ndx,ndy,1,1,1)

      if(jvft44 .ge. 960511) then
        do is=1,nmoli
          call gfsub2(44,ndx,ndy,nnx,nny,srcml(:,:,is,1))
        end do
        do is=1,nmoli
          call gfsub2(44,ndx,ndy,nnx,nny,edissml(:,:,is,0))
        end do
      else
        srcml = 0.0_R8
        edissml = 0.0_R8
      end if
      call eirene_neutr(45,nnx,nny,nmoli,srcml,ndx,ndy,nmol,1,1)
      call eirene_neutr(45,nnx,nny,nmoli,edissml,ndx,ndy,nmol,1,1)

      if ((streql(solps_version,'5.0').or.streql(solps_version,'5.2'))
     &    .and. 
     &    (streql(eirene_format,'old').or.streql(eirene_format,'new')))
     & then
        if (jvft44 .ge. 20130210) then
          call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
          write(6,*) 'In this fort.44 version,',
     .               ' the eneutrad data is summed over species!'
          do is = 2, natmi
            eneutrad(:,:,is,0) = 0.0_R8
          end do
        end if
      end if

      if(jvft44 .ge. 960513) then
        read (44,'(3i6)') nlimi,nstsi,nstrai
        write(45,'(3i6)') nlimi,nstsi,nstrai
        print *,'nlimi,nstsi,nstrai:',nlimi,nstsi,nstrai
        call rdneutrs(44,wldnek,1)
        call   neutrs(45,wldnek,1)
        call rdneutrs(44,wldnep,1)
        call   neutrs(45,wldnep,1)
        call rdneutrs(44,wldna,natmi)
        call   neutrs(45,wldna,natmi)
        call rdneutrs(44,ewlda,natmi)
        call   neutrs(45,ewlda,natmi)
        call rdneutrs(44,wldnm,nmoli)
        call   neutrs(45,wldnm,nmoli)
        call rdneutrs(44,ewldm,nmoli)
        call   neutrs(45,ewldm,nmoli)
c*** read the wall geometry
        read (44,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlimi)
        write(45,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlimi)
      end if
cank
cank-960727
      if(jvft44 .ge. 960623) then
        print *,'natmi,nmoli,nioni,nplsi=',
     .           natmi,nmoli,nioni,nplsi
        call rdneutrs(44,wldra,natmi)
        call   neutrs(45,wldra,natmi)
        call rdneutrs(44,wldrm,nmoli)
        call   neutrs(45,wldrm,nmoli)
      end if
      if(jvft44 .ge. 960727) then
        if(nstrai .gt. 1) then
          do istra=1,nstrai
            call rdneutrs(44,wldnek(1,istra),1)
            call rdneutrs(44,wldnep(1,istra),1)
            call rdneutrs(44,wldna(:,:,istra),natmi)
            call rdneutrs(44,ewlda(:,:,istra),natmi)
            call rdneutrs(44,wldnm(:,:,istra),nmoli)
            call rdneutrs(44,ewldm(:,:,istra),nmoli)
            call rdneutrs(44,wldra(:,:,istra),natmi)
            call rdneutrs(44,wldrm(:,:,istra),nmoli)
          end do
          do istra=1,nstrai
            call neutrs(45,wldnek(1,istra),1)
            call neutrs(45,wldnep(1,istra),1)
            call neutrs(45,wldna(1,1,istra),natmi)
            call neutrs(45,ewlda(1,1,istra),natmi)
            call neutrs(45,wldnm(1,1,istra),nmoli)
            call neutrs(45,ewldm(1,1,istra),nmoli)
            call neutrs(45,wldra(1,1,istra),natmi)
            call neutrs(45,wldrm(1,1,istra),nmoli)
          end do
        end if
      end if
      if(jvft44 .ge. 961228) then
        call rdneutrs(44,wldpp,nplsi)
        call rdneutrs(44,wldpa,natmi)
        call rdneutrs(44,wldpm,nmoli)
        call rdneutrs(44,wldpeb,1)
        call rdneutrs(44,wldspt,1)
        call   neutrs(45,wldpp,nplsi)
        call   neutrs(45,wldpa,natmi)
        call   neutrs(45,wldpm,nmoli)
        call   neutrs(45,wldpeb,1)
        call   neutrs(45,wldspt,1)
        if(nstrai .gt. 1) then
          do istra=1,nstrai
            call rdneutrs(44,wldpp(:,:,istra),nplsi)
            call rdneutrs(44,wldpa(:,:,istra),natmi)
            call rdneutrs(44,wldpm(:,:,istra),nmoli)
            call rdneutrs(44,wldpeb(1,istra),1)
            call rdneutrs(44,wldspt(1,istra),1)
          end do
          do istra=1,nstrai
            call neutrs(45,wldpp(1,1,istra),nplsi)
            call neutrs(45,wldpa(1,1,istra),natmi)
            call neutrs(45,wldpm(1,1,istra),nmoli)
            call neutrs(45,wldpeb(1,istra),1)
            call neutrs(45,wldspt(1,istra),1)
          end do
        end if
      end if
cank
      if(jvftnew.le.961228) goto 99
      if(jvft44 .ge. 20000727) then
        read (44,'(18i4)') (isrftype(i),i=1,nlimi)
        read (44,'(18i4)') (isrftype(nlim+i),i=1,nstsi)
        write(45,'(18i4)') (isrftype(i),i=1,nlimi)
        write(45,'(18i4)') (isrftype(nlim+i),i=1,nstsi)
      end if
      if (NLWRMSH) then
        if(jvft44 .ge. 20051115) then
          read (44,'(6e13.5)') (wlarea(i),i=1,nlimi)
          read (44,'(6e13.5)') (wlarea(nlim+i),i=1,nstsi)
          write(45,'(6e13.5)') (wlarea(i),i=1,nlimi)
          write(45,'(6e13.5)') (wlarea(nlim+i),i=1,nstsi)
        end if
        if(jvft44 .ge. 20060206) then
          k=0
          read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,natmi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          k=k+natmi
          read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nmoli),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          k=k+nmoli
          read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nioni),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nplsi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nplsi),
     ,                          i=1,nstsi)
          k=0
          write(45,'(6a13)') (eirtxt(j+k),j=1,natmi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,i),j=1,natmi),i=1,nlimi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          k=k+natmi
          write(45,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          write(45,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nmoli),i=1,nlimi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          k=k+nmoli
          write(45,'(6a13)') (eirtxt(j+k),j=1,nioni)
          write(45,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nioni),i=1,nlimi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          write(45,'(6a13)') (eirtxt(j+k),j=1,nplsi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nplsi),i=1,nlimi)
          write(45,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nplsi),
     ,                          i=1,nstsi)
        end if

        if (jvft44 .ge. 20071209) then
          if (streql(solps_version,'4.3').or.
     &        streql(solps_version,'5.1').or.
     &        streql(eirene_format,'facelift')) then
            call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
            write(6,*) 'In this fort.44 version,',
     .                 ' the eneutrad data is summed over species!'
            do is = 2, natmi
              eneutrad(:,:,is,0) = 0.0_R8
            end do
          else
            do is = 1, natmi
              call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,is,0))
            end do
          end if
          call eirene_neutr(45,nnx,nny,natmi,eneutrad(:,:,:,0),
     .                         ndx,ndy,natm,1,1)
        end if
      
        if(jvft44 .ge. 20080706) then
          READ (44,'(12i6)') (eirdiag_nds_ind(i),i=1,NSTS+1)        
          READ (44,'(12i6)') (eirdiag_nds_typ(i),i=1,NSTS)    
          READ (44,'(12i6)') (eirdiag_nds_srf(i),i=1,NSTS)    
          READ (44,'(12i6)') (eirdiag_nds_start(i),i=1,NSTS)    
          READ (44,'(12i6)') (eirdiag_nds_end(i),i=1,NSTS) 
          WRITE(45,'(12i6)') (eirdiag_nds_ind(i),i=1,NSTS+1)        
          WRITE(45,'(12i6)') (eirdiag_nds_typ(i),i=1,NSTS)    
          WRITE(45,'(12i6)') (eirdiag_nds_srf(i),i=1,NSTS)    
          WRITE(45,'(12i6)') (eirdiag_nds_start(i),i=1,NSTS)    
          WRITE(45,'(12i6)') (eirdiag_nds_end(i),i=1,NSTS) 
      
          NCL=eirdiag_nds_ind(NSTS+1)
          IF(.NOT.ALLOCATED(wldna_res)) 
     >     ALLOCATE(wldna_res(NATM,NCL), wldnm_res(NMOL,NCL),
     .              ewlda_res(NATM,NCL), ewldm_res(NMOL,NCL), 
     .              ewldea_res(NATM,NCL), ewldem_res(NMOL,NCL), 
     .              ewldrp_res(NCL), ewldmr_res(NMOL,NCL),
     .              sarea_res(NCL) ) 
 
          READ (44,'(6e15.7)') (sarea_res(i),i=1,NCL)
          READ (44,'(6e15.7)') ((wldna_res(j,i),j=1,natm),i=1,NCL)
          READ (44,'(6e15.7)') ((wldnm_res(j,i),j=1,nmol),i=1,NCL) 
          READ (44,'(6e15.7)') ((ewlda_res(j,i),j=1,natm),i=1,NCL)
          READ (44,'(6e15.7)') ((ewldm_res(j,i),j=1,nmol),i=1,NCL)
          READ (44,'(6e15.7)') ((ewldea_res(j,i),j=1,natm),i=1,NCL)
          READ (44,'(6e15.7)') ((ewldem_res(j,i),j=1,nmol),i=1,NCL)
          READ (44,'(6e15.7)') (ewldrp_res(i),i=1,NCL) 
          READ (44,'(6e15.7)') ((ewldmr_res(j,i),j=1,nmol),i=1,NCL) 
          WRITE(45,'(6e15.7)') (sarea_res(i),i=1,NCL)
          WRITE(45,'(6e15.7)') ((wldna_res(j,i),j=1,natm),i=1,NCL)
          WRITE(45,'(6e15.7)') ((wldnm_res(j,i),j=1,nmol),i=1,NCL) 
          WRITE(45,'(6e15.7)') ((ewlda_res(j,i),j=1,natm),i=1,NCL)
          WRITE(45,'(6e15.7)') ((ewldm_res(j,i),j=1,nmol),i=1,NCL)
          WRITE(45,'(6e15.7)') ((ewldea_res(j,i),j=1,natm),i=1,NCL)
          WRITE(45,'(6e15.7)') ((ewldem_res(j,i),j=1,nmol),i=1,NCL)
          WRITE(45,'(6e15.7)') (ewldrp_res(i),i=1,NCL) 
          WRITE(45,'(6e15.7)') ((ewldmr_res(j,i),j=1,nmol),i=1,NCL) 
        end if

        if(streql(solps_version,'4.3')) goto 99
        if(jvft44 .ge. 20081111) then
          READ (44,'(6e15.7)') ((PDENA_INT(I,J),I=1,NATM),J=0,NSTRAI)
          READ (44,'(6e15.7)') ((PDENM_INT(I,J),I=1,NMOL),J=0,NSTRA)
          READ (44,'(6e15.7)') ((PDENI_INT(I,J),I=1,NION),J=0,NSTRA)
          READ (44,'(6e15.7)') ((PDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRA)
          READ (44,'(6e15.7)') ((PDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRA)
          READ (44,'(6e15.7)') ((PDENI_INT_B2(I,J),I=1,NION),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENA_INT(I,J),I=1,NATM),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENM_INT(I,J),I=1,NMOL),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENI_INT(I,J),I=1,NION),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRA)
          READ (44,'(6e15.7)') ((EDENI_INT_B2(I,J),I=1,NION),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((PDENA_INT(I,J),I=1,NATM),J=0,NSTRAI)
          WRITE(45,'(6e15.7)') ((PDENM_INT(I,J),I=1,NMOL),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((PDENI_INT(I,J),I=1,NION),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((PDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((PDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((PDENI_INT_B2(I,J),I=1,NION),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENA_INT(I,J),I=1,NATM),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENM_INT(I,J),I=1,NMOL),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENI_INT(I,J),I=1,NION),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENA_INT_B2(I,J),I=1,NATM),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENM_INT_B2(I,J),I=1,NMOL),J=0,NSTRA)
          WRITE(45,'(6e15.7)') ((EDENI_INT_B2(I,J),I=1,NION),J=0,NSTRA)
        end if
      else  ! NLWRMSH
!! FIXME: If coming from 5.1, date is 20130210 and only sum over species is provided
        if(streql(solps_version,'5.1')) then
          if (jvft44.ge.20130210) then
            call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
            write(6,*) 'In this fort.44 version,',
     .                 ' the eneutrad data is summed over species!'
            do is = 2, natmi
              eneutrad(:,:,is,0) = 0.0_R8
            end do
            call eirene_neutr(45,nnx,nny,natmi,eneutrad(:,:,:,0),
     .                           ndx,ndy,natm,1,1)
          end if
        elseif(jvft44 .ge. 20071209) then
          do is = 1, natmi
            call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,is,0))
          end do
          call eirene_neutr(45,nnx,nny,natmi,eneutrad(:,:,:,0),
     .                         ndx,ndy,natm,1,1)
        end if
      end if

#endif

*  ---------------------------------------------------------------------
   99 CONTINUE
C
      call subend ()
      RETURN

*=======================================================================
*//END B2YT_NGREAD//

      END



