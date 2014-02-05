c
C RHW 18.08.92:  UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
C
*//NGREAD//
*=======================================================================
*                 S U B R O U T I N E  N G R E A D
*=======================================================================
*
*
*  1. PURPOSE
*
*     NGREAD READS IN THE ATOM AND MOLECULE DATA OF EIRENE FOR PLOT.
*
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE NGREAD (ismain,NFEND)
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_neutrals_namelist
      use b2mod_b2cmpa
csw 02feb2011
      use eirmod_extrab25
csw
      implicit none

*  -- INPUT ARGUMENTS --
      INTEGER ismain

*  -- OUTPUT ARGUMENTS --
      INTEGER NFEND

!#     include <COUPLE/EIRDIAG.F>                                !srv 14.10.13

*  -- LOCAL SCALARS --
      INTEGER NNX, NNY, IS, IX, IY, jvft44, istra, nplsi

*  -- PROCEDURES --
      EXTERNAL GFSUB2
      external subini, subend

*=======================================================================
*//COMPUTATION//
C
      call subini ('ngread')
      NFEND = 0

      nplsi = 0
      do is = 0, ns-1
        if (.not.is_neutral(is)) nplsi = nplsi + 1
      end do
C
* ----------------------------------------------------------------------
*  READ EIRENE SPECIES DENSITIES AND TEMPERATURES
* ----------------------------------------------------------------------
*

      rewind (44)
C
      READ(44,'(I4,2X,I4,2x,i8)',END=99,ERR=99) nnx,nny,jvft44
      READ(44,'(I4,2X,I4,2X,I4)',END=99,ERR=99) NATMI,NMOLI,NIONI
      write(*,*) ' Dimensions in fort.44 = ',nnx,nny,natmi,nmoli,nioni
      write(*,*) ' Outer dimensions = ', ndxp, ndyp
      if(natmi .gt. natm) stop 'NGREAD: NATMI > NATM'
      if(nmoli .gt. nmol) stop 'NGREAD: NMOLI > NMOL'
      if(nioni .gt. nion) stop 'NGREAD: NIONI > NION'
      if(nnx   .gt. nxdd) stop 'NGREAD:   NNX > NXDD'
      if(nny   .gt. nydd) stop 'NGREAD:   NNY > NYDD'
      call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstra,natmi)
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ(44,'(A8)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ(44,'(A8)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
C
* test ions
C
      do 13 IS = 1, NIONI
        READ(44,'(A8)',END=99,ERR=99) TEXTIN(IS-1)
   13 CONTINUE
C
* data
* neutral atoms
C
      do 14 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,dab2(:,:,is,1))
   14 CONTINUE
      do is = 1,natmi                                                      !srv 04.06.09 16.02.12{
       call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,is,1),'dab2')
      enddo

      do 24 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,tab2(:,:,is,1))
   24 CONTINUE

* neutral molecules
C
      do 15 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,dmb2(:,:,is,1))
   15 CONTINUE
      do is = 1,nmoli                                                      !srv 16.02.12
       call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,is,1),'dmb2')          !srv 29.04.10
      enddo                                                                 !srv 04.06.09 }

      do 25 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,tmb2(:,:,is,1))
   25 CONTINUE

* test ions
C
      do 16 IS = 1, NIONI
        call gfsub2(44,nxdd,nydd,nnx,nny,dib2(:,:,is,1))
   16 CONTINUE

      do 26 IS = 1, NIONI
        call gfsub2(44,nxdd,nydd,nnx,nny,tib2(:,:,is,1))
   26 CONTINUE

* neutral atom radial particle flux
C
      do 114 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,rfluxa(:,:,is,1))
 114  CONTINUE

* neutral molecules radial particle flux
C
      do 115 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,rfluxm(:,:,is,1))
 115  CONTINUE

* neutral atom poloidal particle flux
C
      do 214 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,pfluxa(:,:,is,1))
 214  CONTINUE

* neutral molecules poloidal particle flux
C
      do 215 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,pfluxm(:,:,is,1))
 215  CONTINUE

* neutral atom radial energy flux
C
      do 314 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,refluxa(:,:,is,1))
 314  CONTINUE

* neutral molecules radial energy flux
C
      do 415 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,refluxm(:,:,is,1))
 415  CONTINUE

* neutral atom poloidal energy flux
C
      do 514 IS = 1, NATMI
        call gfsub2(44,nxdd,nydd,nnx,nny,pefluxa(:,:,is,1))
 514  CONTINUE

* neutral molecules poloidal energy flux
C
      do 515 IS = 1, NMOLI
        call gfsub2(44,nxdd,nydd,nnx,nny,pefluxm(:,:,is,1))
 515  CONTINUE

      call gfsub2(44,nxdd,nydd,nnx,nny,emiss)
      call gfsub2(44,nxdd,nydd,nnx,nny,emissmol)

cank 960511
      do is=1,nmoli
        do ix=0,nx+1
          do iy=0,ny+1
            srcml(ix,iy,is,1)=0.0_R8
            edissml(ix,iy,is,1)=0.0_R8
          end do
        end do
      end do
      do istra=0,nstra
        do ix=1,nlmpgs
          wldnek(ix,istra)=0.0_R8
          wldnep(ix,istra)=0.0_R8
          do is=1,natmi
            wldna(ix,is,istra)=0.0_R8
            wldra(ix,is,istra)=0.0_R8
            ewlda(ix,is,istra)=0.0_R8
          end do
          do is=1,nmoli
            wldnm(ix,is,istra)=0.0_R8
            wldrm(ix,is,istra)=0.0_R8
            ewldm(ix,is,istra)=0.0_R8
          end do
        end do
      end do
      if(jvft44 .ge. 960511) then
        do is=1,nmoli
          call gfsub2(44,nxdd,nydd,nnx,nny,srcml)
          call gfsub2(44,nxdd,nydd,nnx,nny,edissml)
        end do
        read(44,'(2i6)') nlimi,nstsi
        print *,'nlimi,nstsi:',nlimi,nstsi
        call rdneutrs(44,wldnek,1)
        call rdneutrs(44,wldnep,1)
        call rdneutrs(44,wldna,natmi)
        call rdneutrs(44,ewlda,natmi)
        call rdneutrs(44,wldnm,nmoli)
        call rdneutrs(44,ewldm,nmoli)
c*** read the wall geometry
        read (44,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlimi)
      end if
cank
cank-960727
      if(jvft44 .ge. 960623) then
        print *,'nstrai,natmi,nmoli,nplsi=',nstrai,natmi,nmoli,nplsi
        call rdneutrs(44,wldra,natmi)
        call rdneutrs(44,wldrm,nmoli)
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
        else
          do ix=1,nlmpgs
            wldnek(ix,1)=wldnek(ix,0)
            wldnep(ix,1)=wldnep(ix,0)
            do is=1,natmi
              wldna(ix,is,1)=wldna(ix,is,0)
              wldra(ix,is,1)=wldra(ix,is,0)
              ewlda(ix,is,1)=ewlda(ix,is,0)
            end do
            do is=1,nmoli
              wldnm(ix,is,1)=wldnm(ix,is,0)
              wldrm(ix,is,1)=wldrm(ix,is,0)
              ewldm(ix,is,1)=ewldm(ix,is,0)
            end do
          end do
        end if
        if(nspecies .eq. natmi) then
          do ix=1,nlmpgs
            wldnek(ix,0)=0.0_R8
            wldnep(ix,0)=0.0_R8
            do is=1,natmi
              wldna(ix,is,0)=0.0_R8
              wldra(ix,is,0)=0.0_R8
              ewlda(ix,is,0)=0.0_R8
            end do
            do is=1,nmoli
              wldnm(ix,is,0)=0.0_R8
              wldrm(ix,is,0)=0.0_R8
              ewldm(ix,is,0)=0.0_R8
            end do
            do istra=1,nstrai
              wldnek(ix,0)=wldnek(ix,0)+
     +         wldnek(ix,istra)*recyc(ismain,istra)
              wldnep(ix,0)=wldnep(ix,0)+
     +         wldnep(ix,istra)*recyc(ismain,istra)
              iy=1
              do is=1,natmi
                wldna(ix,is,0)=wldna(ix,is,0)+
     +                             wldna(ix,is,istra)*recyc(iy,istra)
                wldra(ix,is,0)=wldra(ix,is,0)+
     +                             wldra(ix,is,istra)*recyc(iy,istra)
                ewlda(ix,is,0)=ewlda(ix,is,0)+
     +                             ewlda(ix,is,istra)*recyc(iy,istra)
                iy=iy+nfluids(is)+1
              end do
              do is=1,nmoli
                if (lmolscl(is).gt.0 .and. lmolscl(is).le.natmi) then
                  iy=eb2spcr(lmolscl(is))+1
                  wldnm(ix,is,0)=wldnm(ix,is,0)+
     +             wldnm(ix,is,istra)*recyc(iy,istra)
                  wldrm(ix,is,0)=wldrm(ix,is,0)+
     +             wldrm(ix,is,istra)*recyc(iy,istra)
                  ewldm(ix,is,0)=ewldm(ix,is,0)+
     +             ewldm(ix,is,istra)*recyc(iy,istra)
                else
                  write(*,*) 'Faulty of missing LMOLSCL assignment:'
                  write(*,*) 'LMOLSCL(',is,')=',lmolscl(is)
                  write(*,*) 'Should be valued between 1 and NATM'
                  write(*,*) 'Modify in b2.neutrals.parameters'
                  write(*,*) 'See manual or b2cdcn.f file for usage'
                  write(*,*) 'Skipping scaling of molecular arrays'
                  wldnm(ix,is,0)=wldnm(ix,is,0)+wldnm(ix,is,istra)
                  wldrm(ix,is,0)=wldrm(ix,is,0)+wldrm(ix,is,istra)
                  ewldm(ix,is,0)=ewldm(ix,is,0)+ewldm(ix,is,istra)
                endif
              end do
            end do
          end do
        else
          print *,'natmi, nspecies= ',natmi,nspecies
          print *,'        --> no source rescaling!'
        end if
      end if
      if(jvft44 .ge. 961228) then
        call rdneutrs(44,wldpp,nplsi)
        call rdneutrs(44,wldpa,natmi)
        call rdneutrs(44,wldpm,nmoli)
        call rdneutrs(44,wldpeb,1)
        call rdneutrs(44,wldspt,1)
        if(nstrai .gt. 1) then
          do istra=1,nstrai
            call rdneutrs(44,wldpp(:,:,istra),nplsi)
            call rdneutrs(44,wldpa(:,:,istra),natmi)
            call rdneutrs(44,wldpm(:,:,istra),nmoli)
            call rdneutrs(44,wldpeb(1,istra),1)
            call rdneutrs(44,wldspt(1,istra),1)
          end do
        end if
      end if
cank

* -- NO ERROR --

      call extend(dab2,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(tab2,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(dmb2,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(tmb2,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(dib2,nxdd,nydd,nion,nnx,nny,nioni)
      call extend(tib2,nxdd,nydd,nion,nnx,nny,nioni)
      call extend(rfluxa,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(rfluxm,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(pfluxa,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(pfluxm,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(refluxa,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(refluxm,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(pefluxa,nxdd,nydd,natm,nnx,nny,natmi)
      call extend(pefluxm,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(emiss,nxdd,nydd,1,nnx,nny,1)
      call extend(emissmol,nxdd,nydd,1,nnx,nny,1)
cank 960511
      call extend(srcml,nxdd,nydd,nmol,nnx,nny,nmoli)
      call extend(edissml,nxdd,nydd,nmol,nnx,nny,nmoli)
cank

      NFEND = 1

*  ---------------------------------------------------------------------
   99 CONTINUE
C
      call subend ()
      RETURN

*=======================================================================
*//END NGREAD//

      END



