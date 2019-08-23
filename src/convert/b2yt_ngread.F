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

      subroutine b2yt_ngread(solps_version, eirene_format,
     . nx1, ny1, natm1, iatm1, nmol1, imol1, npls1, ipls1, nfla,
     . no_sputtering, savestratum, non_commensurate, e_bundle_contr,
     . jsep1, region1, leftcut1, rightcut1, bottomcut1, topcut1,
     . bottomiy1,
     . hx, hy, qc, hx1, hy1, qc1, gs1)
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_clogau
#endif
      implicit none

*  -- INPUT ARGUMENTS --
      integer nx1,ny1,natm1,nmol1,npls1,nfla,non_commensurate
      integer iatm1(natm), imol1(nmol), ipls1(npls)
      integer jsep1
      integer region1(-1:nx1,-1:ny1,0:2)
      integer leftcut1(nncutmax), rightcut1(nncutmax)
      integer bottomcut1(nncutmax), topcut1(nncutmax)
      integer bottomiy1(-1:nx1,-1:ny1)
      real (kind=R8) :: e_bundle_contr(npls1,nfla)
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qc(-1:nx,-1:ny),
     *  hx1(-1:nx1,-1:ny1), hy1(-1:nx1,-1:ny1), qc1(-1:nx1,-1:ny1),
     *  gs1(-1:nx1,-1:ny1,0:2)
      character*4 solps_version
      character*8 eirene_format
      logical savestratum(0:nstraid), no_sputtering

*  -- LOCAL SCALARS --
      integer is, nplsi
#ifdef B25_EIRENE
      INTEGER NNX, NNY, IX, IY, jvft44, istra, istra1
      integer i, j, k, ncl, ncl1, nclm, jvftnew, nstra1, n2
      real (kind=R8) :: wrk0(-1:nx,-1:ny1), field(0:nx+1,0:ny+1)
      real (kind=R8), allocatable :: dab1(:,:,:), tab1(:,:,:)
      real (kind=R8), allocatable :: dmb1(:,:,:), tmb1(:,:,:)
      real (kind=R8), allocatable :: dib1(:,:,:), tib1(:,:,:)
      real (kind=R8), allocatable :: rfa1(:,:,:), pfa1(:,:,:)
      real (kind=R8), allocatable :: rfm1(:,:,:), pfm1(:,:,:)
      real (kind=R8), allocatable :: refa1(:,:,:), pefa1(:,:,:)
      real (kind=R8), allocatable :: refm1(:,:,:), pefm1(:,:,:)
      real (kind=R8), allocatable :: emiss1(:,:), emissmol1(:,:)
      real (kind=R8), allocatable :: srcml1(:,:,:), edissml1(:,:,:,:)
      real (kind=R8), allocatable :: wldna1(:,:,:), ewlda1(:,:,:)
      real (kind=R8), allocatable :: wldnm1(:,:,:), ewldm1(:,:,:)
      real (kind=R8), allocatable :: wldra1(:,:,:), wldrm1(:,:,:)
      real (kind=R8), allocatable :: wldpa1(:,:,:), wldpm1(:,:,:)
      real (kind=R8), allocatable :: wldpp1(:,:,:)
      real (kind=R8), allocatable :: wlabsrp1(:,:)
      real (kind=R8), allocatable :: wlpump1(:,:)
      real (kind=R8), allocatable :: eneutrad1(:,:,:)
      real (kind=R8), allocatable :: wldna1res(:,:), wldnm1res(:,:)
      real (kind=R8), allocatable :: ewlda1res(:,:), ewldm1res(:,:)
      real (kind=R8), allocatable :: ewldea1res(:,:), ewldem1res(:,:)
      real (kind=R8), allocatable :: wldspta1res(:,:), wldsptm1res(:,:)
      real (kind=R8), allocatable :: wldspta1(:,:,:), wldsptm1(:,:,:)
      real (kind=R8), allocatable :: ewldmr1res(:,:), sarea1res(:)
      real (kind=R8), allocatable :: wlpump1res(:,:), ewldrp1res(:)
      real (kind=R8), allocatable :: wldspt1res(:)
      real (kind=R8), allocatable :: data1res(:)
      real (kind=R8), allocatable :: area_res_contr(:,:)
      real (kind=R8), allocatable ::
     *  pdena1int(:,:), pdenm1int(:,:), pdeni1int(:,:),
     *  pdena1int_b2(:,:), pdenm1int_b2(:,:), pdeni1int_b2(:,:),
     *  edena1int(:,:), edenm1int(:,:), edeni1int(:,:),
     *  edena1int_b2(:,:), edenm1int_b2(:,:), edeni1int_b2(:,:)
      real (kind=R8), allocatable :: wd_old(:), wd_new(:)
      integer, allocatable :: eirdiag_nds_ind1(:), eirdiag_nds_srf1(:),
     *                        eirdiag_nds_start1(:), eirdiag_nds_end1(:)
      character*16, allocatable :: textan1(:), textmn1(:)
      character*8, allocatable :: eirtxt1(:)
      character*3 stnum
      character*256 filename
      logical ex

*  -- PROCEDURES --
      integer convert_eirdiag, conformal_convert
      character*32 :: get_B25_hash
      logical streql
      EXTERNAL GFSUB2
      external ipgeti, streql, get_B25_hash
      external convert_eirdiag, conformal_convert
#endif
      external subini, subend

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
#ifndef HPUX
      open(unit=45,file='fort.44.2',status='unknown')
      write(*,*) 'Writing fort.44.2'
#else
      open(unit=45,file='ftn44.2',status='unknown')
      write(*,*) 'Writing ftn44.2'
#endif
C
      READ (44,'(I4,2X,I4,2x,i8)',END=99,ERR=99) nnx,nny,jvft44
      if((streql(solps_version,'5.0').or.streql(solps_version,'5.2'))
     &   .and.
     &   (streql(eirene_format,'old').or.streql(eirene_format,'new')))
     & then
        jvftnew = min(jvft44,961228)
      else if (streql(solps_version,'4.3')) then
        jvftnew = min(jvft44,20080706)
      else if (streql(solps_version,'5.1')) then
        jvftnew = min(jvft44,20130210)
      else
        jvftnew = max(jvft44,20160829)
      endif
      WRITE(45,'(I4,2X,I4,2x,i8,2x,a32)') nx1,ny1,jvftnew,get_B25_hash()
      READ (44,'(I4,2X,I4,2X,I4)',END=99,ERR=99) NATMI,NMOLI,NIONI
      WRITE(45,'(I4,2X,I4,2X,i4)')               natm1,NMOL1,NIONI
#ifndef HPUX
      write(*,'(a,5i4)')
     . ' Dimensions in fort.44 =   ',nnx,nny,natmi,nmoli,nioni
      write(*,'(a,5i4)')
     . ' Dimensions in fort.44.2 = ',nx1,ny1,natm1,nmol1,nioni
#else
      write(*,'(a,5i4)')
     . ' Dimensions in ftn44 =   ',nnx,nny,natmi,nmoli,nioni
      write(*,'(a,5i4)')
     . ' Dimensions in ftn44.2 = ',nx1,ny1,natm1,nmol1,nioni
#endif
      write(*,'(a,2i4)') ' Array inner dimensions = ', ndx , ndy
      call xertst(natmi .le. natm, 'NGREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NGREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NGREAD: NIONI > NION')
      call xertst(nnx   .le. ndx , 'NGREAD:   NNX > NDX ')
      call xertst(nny   .le. ndy , 'NGREAD:   NNY > NDY ')
      call alloc_b2mod_diag(nx,ny,ns,nnreg(0),nstra,natmi)
      n2 = nx1*ny1
      filename='input.dat'
      call find_file(filename,ex)
#ifndef HPUX
      filename='fort.33'
#else
      filename='ftn33'
#endif
      if (.not.ex) call find_file(filename,nlwrmsh)
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ (44,'(A8)',END=99,ERR=99) TEXTAN(IS-1)
   11 CONTINUE
      allocate(textan1(0:natm1-1))
      do is = 1, natmi
        if (iatm1(is).ne.0) TEXTAN1(iatm1(is)-1) = TEXTAN(is-1)
      end do
      do is = 1, natm1
        WRITE(45,'(A8)') TEXTAN1(IS-1)
      end do
      deallocate(textan1)
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ (44,'(A8)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
      allocate(textmn1(0:nmol1-1))
      do is = 1, nmoli
        if (imol1(is).ne.0) TEXTMN1(imol1(is)-1) = TEXTMN(is-1)
      end do
      do is = 1, nmol1
        WRITE(45,'(A8)') TEXTMN1(IS-1)
      end do
      deallocate(textmn1)
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
      allocate (dab1(0:nx1+1,0:ny1+1,natm1))
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = dab2(0:nx+1,0:ny+1,is,1)
          call b2xtvc(nx, ny, nx1, ny1,
     &        field, dab1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
        end if
      end do
      call write_title(45,'dab2',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,dab1,nx1,ny1,natm1,1,1)
      deallocate (dab1)

      do 24 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,tab2(:,:,is,1))
   24 CONTINUE
      allocate (tab1(0:nx1+1,0:ny1+1,natm1))
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = tab2(0:nx+1,0:ny+1,is,1)
          if (streql(solps_version,'4.3')) then
            do iy = 0, ny+1
              do ix = 0, nx+1
                if (dab2(ix,iy,is,1).eq.0.0_R8) then
                  field(ix,iy) = field(ix,iy)*ev
                end if
              end do
            end do
          end if
          call b2xtvc(nx, ny, nx1, ny1,
     &        field, tab1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
        end if
      end do
      call write_title(45,'tab2',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,tab1,nx1,ny1,natm1,1,1)
      deallocate (tab1)

* neutral molecules
C
      do 15 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,dmb2(:,:,is,1))
   15 CONTINUE
      allocate (dmb1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = dmb2(0:nx+1,0:ny+1,is,1)
          call b2xtvc(nx, ny, nx1, ny1,
     &        field, dmb1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
        end if
      end do
      call write_title(45,'dmb2',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,dmb1,nx1,ny1,nmol1,1,1)
      deallocate(dmb1)

      do 25 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,tmb2(:,:,is,1))
   25 CONTINUE
      allocate (tmb1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = tmb2(0:nx+1,0:ny+1,is,1)
          if (streql(solps_version,'4.3')) then
            do iy = 0, ny+1
              do ix = 0, nx+1
                if (dmb2(ix,iy,is,1).eq.0.0_R8) then
                  field(ix,iy) = field(ix,iy)*ev
                end if
              end do
            end do
          end if
          call b2xtvc(nx, ny, nx1, ny1,
     &        field, tmb1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
        end if
      end do
      call write_title(45,'tmb2',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,tmb1,nx1,ny1,nmol1,1,1)
      deallocate(tmb1)

* test ions
C
      allocate (dib1(0:nx1+1,0:ny1+1,nioni))
      do 16 IS = 1, NIONI
        call gfsub2(44,ndx,ndy,nnx,nny,dib2(:,:,is,1))
        field(0:nx+1,0:ny+1) = dib2(0:nx+1,0:ny+1,is,1)
        call b2xtvc(nx, ny, nx1, ny1,
     &        field, dib1(0,0,is),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
   16 CONTINUE
      call write_title(45,'dib2',n2*nioni)
      call eirene_neutr(45,nx1,ny1,nioni,dib1,nx1,ny1,nioni,1,1)
      deallocate(dib1)

      allocate (tib1(0:nx1+1,0:ny1+1,nioni))
      do 26 IS = 1, NIONI
        call gfsub2(44,ndx,ndy,nnx,nny,tib2(:,:,is,1))
        field(0:nx+1,0:ny+1) = tib2(0:nx+1,0:ny+1,is,1)
        if (streql(solps_version,'4.3')) then
          do iy = 0, ny+1
            do ix = 0, nx+1
              if (dib2(ix,iy,is,1).eq.0.0_R8) then
                field(ix,iy) = field(ix,iy)*ev
              end if
            end do
          end do
        end if
        call b2xtvc(nx, ny, nx1, ny1,
     &        field, tib1(0,0,is),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
   26 CONTINUE
      call write_title(45,'tib2',n2*nioni)
      call eirene_neutr(45,nx1,ny1,nioni,tib1,nx1,ny1,nioni,1,1)
      deallocate(tib1)

* neutral atom radial particle flux
C
      do 114 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,rfluxa(:,:,is,1))
  114 CONTINUE
      allocate (rfa1(0:nx1+1,0:ny1+1,natm1))
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = rfluxa(0:nx+1,0:ny+1,is,1)
          call b2xtfy (nx, ny, nx1, ny1,
     &        field, rfa1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hx, hx1, bottomiy1)
        end if
      end do
      call write_title(45,'rfluxa',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,rfa1,nx1,ny1,natm1,1,1)
      deallocate (rfa1)

* neutral molecules radial particle flux
C
      do 115 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,rfluxm(:,:,is,1))
  115 CONTINUE
      allocate (rfm1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = rfluxm(0:nx+1,0:ny+1,is,1)
          call b2xtfy (nx, ny, nx1, ny1,
     &        field, rfm1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hx, hx1, bottomiy1)
        end if
      end do
      call write_title(45,'rfluxm',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,rfm1,nx1,ny1,nmol1,1,1)
      deallocate (rfm1)

* neutral atom poloidal particle flux
C
      allocate (pfa1(0:nx1+1,0:ny1+1,natm1))
      do 214 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,pfluxa(:,:,is,1))
  214 CONTINUE
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = pfluxa(0:nx+1,0:ny+1,is,1)
          call b2xtfx (nx, ny, nx1, ny1,
     &        field, pfa1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
        end if
      end do
      call write_title(45,'pfluxa',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,pfa1,nx1,ny1,natm1,1,1)
      deallocate (pfa1)

* neutral molecules poloidal particle flux
C
      do 215 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,pfluxm(:,:,is,1))
  215 CONTINUE
      allocate (pfm1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = pfluxm(0:nx+1,0:ny+1,is,1)
          call b2xtfx (nx, ny, nx1, ny1,
     &        field, pfm1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
        end if
      end do
      call write_title(45,'pfluxm',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,pfm1,nx1,ny1,nmol1,1,1)
      deallocate (pfm1)

* neutral atom radial energy flux
C
      do 314 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,refluxa(:,:,is,1))
  314 CONTINUE
      allocate (refa1(0:nx1+1,0:ny1+1,natm1))
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = refluxa(0:nx+1,0:ny+1,is,1)
          call b2xtfy (nx, ny, nx1, ny1,
     &        field, refa1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hx, hx1, bottomiy1)
        end if
      end do
      call write_title(45,'refluxa',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,refa1,nx1,ny1,natm1,1,1)
      deallocate (refa1)

* neutral molecules radial energy flux
C
      do 415 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,refluxm(:,:,is,1))
  415 CONTINUE
      allocate (refm1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = refluxm(0:nx+1,0:ny+1,is,1)
          call b2xtfy (nx, ny, nx1, ny1,
     &        field, refm1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hx, hx1, bottomiy1)
        end if
      end do
      call write_title(45,'refluxm',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,refm1,nx1,ny1,nmol1,1,1)
      deallocate (refm1)

* neutral atom poloidal energy flux
C
      allocate (pefa1(0:nx1+1,0:ny1+1,natm1))
      do 514 IS = 1, NATMI
        call gfsub2(44,ndx,ndy,nnx,nny,pefluxa(:,:,is,1))
  514 CONTINUE
      do is = 1, natmi
        if (iatm1(is).ne.0) then
          field(0:nx+1,0:ny+1) = pefluxa(0:nx+1,0:ny+1,is,1)
          call b2xtfx (nx, ny, nx1, ny1,
     &        field, pefa1(0,0,iatm1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
        end if
      end do
      call write_title(45,'pefluxa',n2*natm1)
      call eirene_neutr(45,nx1,ny1,natm1,pefa1,nx1,ny1,natm1,1,1)
      deallocate (pefa1)

* neutral molecules poloidal energy flux
C
      do 515 IS = 1, NMOLI
        call gfsub2(44,ndx,ndy,nnx,nny,pefluxm(:,:,is,1))
  515 CONTINUE
      allocate (pefm1(0:nx1+1,0:ny1+1,nmol1))
      do is = 1, nmoli
        if (imol1(is).ne.0) then
          field(0:nx+1,0:ny+1) = pefluxm(0:nx+1,0:ny+1,is,1)
          call b2xtfx (nx, ny, nx1, ny1,
     &        field, pefm1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        hy, qc, hy1, qc1)
        end if
      end do
      call write_title(45,'pefluxm',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,pefm1,nx1,ny1,nmol1,1,1)
      deallocate (pefm1)

      call gfsub2(44,ndx,ndy,nnx,nny,emiss(:,:,1,1))
      call gfsub2(44,ndx,ndy,nnx,nny,emissmol(:,:,1,1))
      allocate (emiss1(0:nx1+1,0:ny1+1))
      allocate (emissmol1(0:nx1+1,0:ny1+1))
      field(0:nx+1,0:ny+1) = emiss(0:nx+1,0:ny+1,1,1)
      call b2xtvc(nx, ny, nx1, ny1,
     &        field, emiss1(0,0),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
      field(0:nx+1,0:ny+1) = emissmol(0:nx+1,0:ny+1,1,1)
      call b2xtvc(nx, ny, nx1, ny1,
     &        field, emissmol1(0,0),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
      call write_title(45,'emiss',n2)
      call eirene_neutr(45,nx1,ny1,1,emiss1,nx1,ny1,1,1,1)
      call write_title(45,'emissmol',n2)
      call eirene_neutr(45,nx1,ny1,1,emissmol1,nx1,ny1,1,1,1)
      deallocate (emiss1)
      deallocate (emissmol1)

      allocate (srcml1(0:nx1+1,0:ny1+1,nmol1))
      allocate (edissml1(0:nx1+1,0:ny1+1,nmol1,0:nstra+1))
      if(jvft44 .ge. 960511) then
        do is=1,nmoli
          call gfsub2(44,ndx,ndy,nnx,nny,srcml(:,:,is,1))
        end do
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field(0:nx+1,0:ny+1) = srcml(0:nx+1,0:ny+1,is,1)
            call b2xtvc(nx, ny, nx1, ny1,
     &        field, srcml1(0,0,imol1(is)),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
          end if
        end do
        do is=1,nmoli
          call gfsub2(44,ndx,ndy,nnx,nny,edissml(:,:,is,0))
        end do
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            field(0:nx+1,0:ny+1) = edissml(0:nx+1,0:ny+1,is,1)
            call b2xtvc(nx, ny, nx1, ny1,
     &        field, edissml1(0,0,imol1(is),0),
     &        nncutmax, non_commensurate,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
          end if
        end do
      else
        srcml1 = 0.0_R8
        edissml1 = 0.0_R8
      end if
      call write_title(45,'srcml',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,srcml1,nx1,ny1,nmol1,1,1)
      call write_title(45,'edissml',n2*nmol1)
      call eirene_neutr(45,nx1,ny1,nmol1,edissml1,nx1,ny1,nmol1,1,1)
      deallocate(srcml1,edissml1)

      if ((streql(solps_version,'5.0').or.streql(solps_version,'5.2'))
     &    .and. 
     &    (streql(eirene_format,'old').or.streql(eirene_format,'new')))
     & then
        if (jvft44 .ge. 20130210) then
          call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
#ifndef HPUX
          write(6,'(2a)') 'In this fort.44 version,',
     .                    ' the eneutrad data is summed over species!'
#else
          write(6,'(2a)') 'In this ftn44 version,',
     .                    ' the eneutrad data is summed over species!'
#endif
          do is = 2, natmi
            eneutrad(:,:,is,0) = 0.0_R8
          end do
        end if
      end if

      if(jvft44 .ge. 960513) then
        nlim1 = nlimi
        read (44,'(3i6)') nlimi,nstsi,nstrai
        if (nlim1.ne.nlimi) then
          if (nlimi.eq.nlim1+1) then
#ifndef HPUX
            write(0,*) 'NLIM in fort.44 file is 1 greater than '//
     ,                 'NLIM in Eirene input file'
#else
            write(0,*) 'NLIM in ftn44 file is 1 greater than '//
     ,                 'NLIM in Eirene input file'
#endif
            write(0,*) 'Assuming dummy wall to be removed!'
          else
#ifndef HPUX
            call xerrab
     .       ('Mismatched NLIM between fort.44 and input file')
#else
            call xerrab
     .       ('Mismatched NLIM between ftn44 and input file')
#endif
          end if
        end if
        nstra1 = 0
        do istra = 1, nstrai
          if (savestratum(istra)) nstra1 = nstra1 + 1
        end do
        allocate(wldna1(nlimps,natm1,0:nstrai+1))
        allocate(ewlda1(nlimps,natm1,0:nstrai+1))
        allocate(wldnm1(nlimps,nmol1,0:nstrai+1))
        allocate(ewldm1(nlimps,nmol1,0:nstrai+1))
        write(45,'(3i6)') nlim1,nstsi,nstra1
        nnlimi=nlim1
        nnstsi=nstsi
        print '(1x,a,4i6)','nlimi,nstsi,nstrai,nstra1:',
     *                      nlim1,nstsi,nstrai,nstra1
        call rdneutrs(44,wldnek,1)
        call rdneutrs(44,wldnep,1)
        call rdneutrs(44,wldna,natmi)
        call rdneutrs(44,ewlda,natmi)
        call rdneutrs(44,wldnm,nmoli)
        call rdneutrs(44,ewldm,nmoli)
        call write_title(45,'wldnek(0)',nlimps)
        call wrneutrs(45,wldnek,1)
        call write_title(45,'wldnep(0)',nlimps)
        call wrneutrs(45,wldnep,1)
        do is = 1, natmi
          if (iatm1(is).ne.0) then
            wldna1(:,iatm1(is),0) = wldna(:,is,0)
            ewlda1(:,iatm1(is),0) = ewlda(:,is,0)
          end if
        end do
        do is = 1, nmoli
          if (imol1(is).ne.0) then
            wldnm1(:,imol1(is),0) = wldnm(:,is,0)
            ewldm1(:,imol1(is),0) = ewldm(:,is,0)
          end if
        end do
        call write_title(45,'wldna(0)',nlimps*natm1)
        call wrneutrs(45,wldna1,natm1)
        call write_title(45,'ewlda(0)',nlimps*natm1)
        call wrneutrs(45,ewlda1,natm1)
        call write_title(45,'wldnm(0)',nlimps*nmol1)
        call wrneutrs(45,wldnm1,nmol1)
        call write_title(45,'ewldm(0)',nlimps*nmol1)
        call wrneutrs(45,ewldm1,nmol1)
c*** read the wall geometry
        call read_title(44)
        read (44,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlimi)
        call write_title(45,'wall_geometry',4*nnlimi)
        write(45,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlim1)
      end if
cank
cank-960727
      if(jvft44 .ge. 960623) then
        print *,'natmi,nmoli,nioni,nplsi=',
     .           natmi,nmoli,nioni,nplsi
        call rdneutrs(44,wldra,natmi)
        call rdneutrs(44,wldrm,nmoli)
        allocate(wldra1(nlimps,natm1,0:nstrai+1))
        allocate(wldrm1(nlimps,nmol1,0:nstrai+1))
        do is = 1, natmi
          if (iatm1(is).ne.0) wldra1(:,iatm1(is),0) = wldra(:,is,0)
        end do
        do is = 1, nmoli
          if (imol1(is).ne.0) wldrm1(:,imol1(is),0) = wldrm(:,is,0)
        end do
        call write_title(45,'wldra(0)',nlimps*natm1)
        call wrneutrs(45,wldra1,natm1)
        call write_title(45,'wldrm(0)',nlimps*nmol1)
        call wrneutrs(45,wldrm1,nmol1)
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
          istra1=0
          do istra=1,nstrai
            if (.not.savestratum(istra)) cycle
            istra1=istra1+1
            do is = 1, natmi
              if (iatm1(is).ne.0) then
                wldna1(:,iatm1(is),istra) = wldna(:,is,istra)
                ewlda1(:,iatm1(is),istra) = ewlda(:,is,istra)
                wldra1(:,iatm1(is),istra) = wldra(:,is,istra)
              end if
            end do
            do is = 1, nmoli
              if (imol1(is).ne.0) then
                wldnm1(:,imol1(is),istra) = wldnm(:,is,istra)
                ewldm1(:,imol1(is),istra) = ewldm(:,is,istra)
                wldrm1(:,imol1(is),istra) = wldrm(:,is,istra)
              end if
            end do
            write (stnum,'(i3)') istra1
            call write_title(45,'wldnek('//stnum//')',nlimps)
            call wrneutrs(45,wldnek(1,istra),1)
            call write_title(45,'wldnep('//stnum//')',nlimps)
            call wrneutrs(45,wldnep(1,istra),1)
            call write_title(45,'wldna('//stnum//')',nlimps*natm1)
            call wrneutrs(45,wldna1(1,1,istra),natm1)
            call write_title(45,'ewlda('//stnum//')',nlimps*natm1)
            call wrneutrs(45,ewlda1(1,1,istra),natm1)
            call write_title(45,'wldnm('//stnum//')',nlimps*nmol1)
            call wrneutrs(45,wldnm1(1,1,istra),nmol1)
            call write_title(45,'ewldm('//stnum//')',nlimps*nmol1)
            call wrneutrs(45,ewldm1(1,1,istra),nmol1)
            call write_title(45,'wldra('//stnum//')',nlimps*natm1)
            call wrneutrs(45,wldra1(1,1,istra),natm1)
            call write_title(45,'wldrm('//stnum//')',nlimps*nmol1)
            call wrneutrs(45,wldrm1(1,1,istra),nmol1)
          end do
        end if
        deallocate(wldna1,wldnm1)
        deallocate(ewlda1,ewldm1)
        deallocate(wldra1,wldrm1)
      end if
      if(jvft44 .ge. 961228) then
        if (streql(solps_version,'5.1')) then
          call rdneutrs(44,wldpp,nfla)
        else
          call rdneutrs(44,wldpp,nplsi)
        endif
        call rdneutrs(44,wldpa,natmi)
        call rdneutrs(44,wldpm,nmoli)
        call rdneutrs(44,wldpeb,1)
        call rdneutrs(44,wldspt,1)
        if (jvft44 .ge. 20170328) then
          call rdneutrs(44,wldspta,natmi)
          call rdneutrs(44,wldsptm,nmoli)
        end if
        allocate(wldpp1(nlimps,npls1,0:nstrai+1))
        allocate(wldpa1(nlimps,natm1,0:nstrai+1))
        allocate(wldpm1(nlimps,nmol1,0:nstrai+1))
        allocate(wd_old(nplsi))
        allocate(wd_new(npls1))
        do i = 1, nlimps
          wd_old(1:nplsi)=wldpp(i,1:nplsi,0)
          call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .     wd_old,wd_new,1)
          wldpp1(i,1:npls1,0)=wd_new(1:npls1)
        end do
        do is = 1, natmi
          if (iatm1(is).ne.0) wldpa1(:,iatm1(is),0) = wldpa(:,is,0)
        end do
        do is = 1, nmoli
          if (imol1(is).ne.0) wldpm1(:,imol1(is),0) = wldpm(:,is,0)
        end do
        call write_title(45,'wldpp(0)',nlimps*npls1)
        call wrneutrs(45,wldpp1,npls1)
        call write_title(45,'wldpa(0)',nlimps*natm1)
        call wrneutrs(45,wldpa1,natm1)
        call write_title(45,'wldpm(0)',nlimps*nmol1)
        call wrneutrs(45,wldpm1,nmol1)
        call write_title(45,'wldpeb(0)',nlimps)
        call wrneutrs(45,wldpeb,1)
        call write_title(45,'wldspt(0)',nlimps)
        if (no_sputtering) wldspt(:,0) = 0.0_R8
        call wrneutrs(45,wldspt,1)
        if (jvftnew .ge. 20170328) then
          allocate(wldspta1(nlimps,natm1,0:nstrai+1))
          allocate(wldsptm1(nlimps,nmol1,0:nstrai+1))
          if (no_sputtering) then
            wldspta1 = 0.0_R8
            wldsptm1 = 0.0_R8
          else
            do is = 1, natmi
              if (iatm1(is).ne.0)
     >          wldspta1(:,iatm1(is),0) = wldspta(:,is,0)
            end do
            do is = 1, nmoli
              if (imol1(is).ne.0)
     >          wldsptm1(:,imol1(is),0) = wldsptm(:,is,0)
            end do
          end if
          call write_title(45,'wldspta(0)',nlimps*natm1)
          call wrneutrs(45,wldspta1,natm1)
          call write_title(45,'wldsptm(0)',nlimps*nmol1)
          call wrneutrs(45,wldsptm1,nmol1)
        end if
        if(nstrai .gt. 1) then
          do istra=1,nstrai
            if (streql(solps_version,'5.1')) then
              call rdneutrs(44,wldpp(:,:,istra),nfla)
            else
              call rdneutrs(44,wldpp(:,:,istra),nplsi)
            endif
            call rdneutrs(44,wldpa(:,:,istra),natmi)
            call rdneutrs(44,wldpm(:,:,istra),nmoli)
            call rdneutrs(44,wldpeb(1,istra),1)
            call rdneutrs(44,wldspt(1,istra),1)
            if (no_sputtering) wldspt(:,istra) = 0.0_R8
            if (jvft44 .ge. 20170328) then
              call rdneutrs(44,wldspta(:,:,istra),natmi)
              call rdneutrs(44,wldsptm(:,:,istra),nmoli)
            end if
          end do
          istra1=0
          do istra=1,nstrai
            if (.not.savestratum(istra)) cycle
            istra1=istra1+1
            do i = 1, nlimps
              wd_old(1:nplsi)=wldpp(i,1:nplsi,istra)
              call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .         wd_old,wd_new,1)
              wldpp1(i,1:npls1,istra)=wd_new(1:npls1)
            end do
            do is = 1, natmi
              if (iatm1(is).ne.0) then
                wldpa1(:,iatm1(is),istra) = wldpa(:,is,istra)
                if (jvftnew .ge. 20170328 .and. .not.no_sputtering) then
                  wldspta1(:,iatm1(is),istra) = wldspta(:,is,istra)
                end if
              end if
            end do
            do is = 1, nmoli
              if (imol1(is).ne.0) then
                wldpm1(:,imol1(is),istra) = wldpm(:,is,istra)
                if (jvftnew .ge. 20170328 .and. .not.no_sputtering) then
                  wldsptm1(:,imol1(is),istra) = wldsptm(:,is,istra)
                end if
              end if
            end do
            write (stnum,'(i3)') istra1
            call write_title(45,'wldpp('//stnum//')',nlimps*npls1)
            call wrneutrs(45,wldpp1(1,1,istra),npls1)
            call write_title(45,'wldpa('//stnum//')',nlimps*natm1)
            call wrneutrs(45,wldpa1(1,1,istra),natm1)
            call write_title(45,'wldpm('//stnum//')',nlimps*nmol1)
            call wrneutrs(45,wldpm1(1,1,istra),nmol1)
            call write_title(45,'wldpeb('//stnum//')',nlimps)
            call wrneutrs(45,wldpeb(1,istra),1)
            call write_title(45,'wldspt('//stnum//')',nlimps)
            call wrneutrs(45,wldspt(1,istra),1)
            if (jvftnew .ge. 20170328) then
              call write_title(45,'wldspta('//stnum//')',nlimps*natm1)
              call wrneutrs(45,wldspta1(1,1,istra),natm1)
              call write_title(45,'wldsptm('//stnum//')',nlimps*nmol1)
              call wrneutrs(45,wldsptm1(1,1,istra),nmol1)
            end if
          end do
        end if
        deallocate(wldpp1)
        deallocate(wldpa1)
        deallocate(wldpm1)
        if (jvftnew .ge. 20170328) deallocate(wldspta1,wldsptm1)
      end if
cank
      if(jvftnew.le.961228) goto 99
      if(jvft44 .ge. 20000727) then
        call read_title(44)
        read (44,'(18i4)') (isrftype(i),i=1,nlimi)
        read (44,'(18i4)') (isrftype(nlim+i),i=1,nstsi)
        call write_title(45,'isrftype',nnlimi+nnstsi)
        write(45,'(18i4)') (isrftype(i),i=1,nlim1)
        write(45,'(18i4)') (isrftype(nlim+i),i=1,nstsi)
      end if
      if (NLWRMSH .and. (streql(solps_version,'4.3') .or.
     &                   streql(solps_version,'5.1') .or.
     &                   streql(solps_version,'iter'))) then
        if(jvft44 .ge. 20051115) then
          call read_title(44)
          read (44,'(6e13.5)') (wlarea(i),i=1,nlimi)
          read (44,'(6e13.5)') (wlarea(nlim+i),i=1,nstsi)
          call write_title(45,'wlarea',nnlimi+nnstsi)
          write(45,'(1p,6e13.5)') (wlarea(i),i=1,nlim1)
          write(45,'(1p,6e13.5)') (wlarea(nlim+i),i=1,nstsi)
        end if
        if(jvft44 .ge. 20060206) then
          k=0
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,natmi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          k=k+natmi
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nmoli),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          k=k+nmoli
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nioni),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          call read_title(44)
          if (streql(solps_version,'5.1')) then
            read (44,'(6a13)') (eirtxt(j+k),j=1,nfla)
            read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nfla),i=1,nlimi)
            read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nfla),
     ,                          i=1,nstsi)
          else
            read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
            read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nplsi),i=1,nlimi)
            read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nplsi),
     ,                          i=1,nstsi)
          endif

          allocate (eirtxt1(natm1+nmol1+nioni+npls1))
          allocate (wlabsrp1(natm1+nmol1+nioni+npls1,nlim+nstsi))
          do k = 1, natmi
            if (iatm1(k).ne.0) then
              eirtxt1(iatm1(k)) = eirtxt(k)
              wlabsrp1(iatm1(k),:) = wlabsrp(k,:)
            end if
          end do
          do k = 1, nmoli
            if (imol1(k).ne.0) then
              eirtxt1(natm1+imol1(k)) = eirtxt(natmi+k)
              wlabsrp1(natm1+imol1(k),:) = wlabsrp(natmi+k,:)
            end if
          end do
          do k = 1, nioni
            eirtxt1(natm1+nmol1+k) = eirtxt(natmi+nmoli+k)
            wlabsrp1(natm1+nmol1+k,:) = wlabsrp(natmi+nmoli+k,:)
          end do
          do k = 1, nplsi
            if (ipls1(k).ne.0) then
              eirtxt1(natm1+nmol1+nioni+ipls1(k)) =
     .         eirtxt(natmi+nmol1+nioni+k)
            end if
          end do
          do i = 1, nlimi
            wd_old(1:nplsi)=wlabsrp(natmi+nmoli+nioni+1:
     .                              natmi+nmoli+nioni+nplsi,i)
            call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .           wd_old,wd_new,2)
           wlabsrp1(natm1+nmol1+nioni+1:natm1+nmol1+nioni+npls1,i)=
     .      wd_new(1:npls1) 
          end do
          do i = 1, nstsi
            wd_old(1:nplsi)=wlabsrp(natmi+nmoli+nioni+1:
     .                              natmi+nmoli+nioni+nplsi,nlim+i)
            call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .           wd_old,wd_new,2)
           wlabsrp1(natm1+nmol1+nioni+1:
     .              natm1+nmol1+nioni+npls1,nlim+i)=wd_new(1:npls1)
          end do
          k=0
          call write_title(45,'wlabsrp(A)',(nnlimi+nnstsi)*natm1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,natm1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,i),j=1,natm1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,nlim+i),j=1,natm1),
     ,                                                       i=1,nstsi)
          k=k+natm1
          call write_title(45,'wlabsrp(M)',(nnlimi+nnstsi)*nmol1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,nmol1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,i),j=1,nmol1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,nlim+i),j=1,nmol1),
     ,                                                       i=1,nstsi)
          k=k+nmol1
          call write_title(45,'wlabsrp(I)',(nnlimi+nnstsi)*nioni)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,nioni)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,i),j=1,nioni),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          call write_title(45,'wlabsrp(P)',(nnlimi+nnstsi)*npls1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,npls1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,i),j=1,npls1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlabsrp1(j+k,nlim+i),j=1,npls1),
     ,                          i=1,nstsi)
          deallocate (wlabsrp1)
        end if

        if(jvft44 .ge. 20170328) then
          k=0
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,natmi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          k=k+natmi
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nmoli),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          k=k+nmoli
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nioni),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nplsi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nplsi),
     ,                          i=1,nstsi)

          allocate (wlpump1(natm1+nmol1+nioni+npls1,nlim+nstsi))
          do k = 1, natmi
            if (iatm1(k).ne.0) wlpump1(iatm1(k),:) = wlpump(k,:)
          end do
          do k = 1, nmoli
            if (imol1(k).ne.0)
     >       wlpump1(natm1+imol1(k),:) = wlpump(natmi+k,:)
          end do
          do k = 1, nioni
            wlpump1(natm1+nmol1+k,:) = wlpump(natmi+nmoli+k,:)
          end do
          do i = 1, nlimi
            wd_old(1:nplsi)=wlpump(natmi+nmoli+nioni+1:
     .                             natmi+nmoli+nioni+nplsi,i)
            call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .           wd_old,wd_new,2)
           wlpump1(natm1+nmol1+nioni+1:natm1+nmol1+nioni+npls1,i)=
     .      wd_new(1:npls1)
          end do
          do i = 1, nstsi
            wd_old(1:nplsi)=wlpump(natmi+nmoli+nioni+1:
     .                             natmi+nmoli+nioni+nplsi,nlim+i)
            call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .           wd_old,wd_new,2)
           wlpump1(natm1+nmol1+nioni+1:
     .             natm1+nmol1+nioni+npls1,nlim+i)=wd_new(1:npls1)
          end do
          k=0
          call write_title(45,'wlpump(A)',(nnlimi+nnstsi)*natm1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,natm1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,i),j=1,natm1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,nlim+i),j=1,natm1),
     ,                                                       i=1,nstsi)
          k=k+natm1
          call write_title(45,'wlpump(M)',(nnlimi+nnstsi)*nmol1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,nmol1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,i),j=1,nmol1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,nlim+i),j=1,nmol1),
     ,                                                       i=1,nstsi)
          k=k+nmol1
          call write_title(45,'wlpump(I)',(nnlimi+nnstsi)*nioni)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,nioni)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,i),j=1,nioni),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          k=k+nioni
          call write_title(45,'wlpump(P)',(nnlimi+nnstsi)*npls1)
          write(45,'(6a13)') (eirtxt1(j+k),j=1,npls1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,i),j=1,npls1),
     ,                                                       i=1,nlim1)
          write(45,'(1p,6e13.5)') ((wlpump1(j+k,nlim+i),j=1,npls1),
     ,                          i=1,nstsi)
          deallocate (wlpump1)
        end if

        if (jvft44 .ge. 20071209) then
          if (streql(solps_version,'4.3').or.
     &        streql(solps_version,'5.1').or.
     &        streql(eirene_format,'facelift')) then
            call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
#ifndef HPUX
            write(6,'(2a)') 'In this fort.44 version,',
     .                      ' the eneutrad data is summed over species!'
#else
            write(6,'(2a)') 'In this ftn44 version,',
     .                      ' the eneutrad data is summed over species!'
#endif
            do is = 2, natmi
              eneutrad(:,:,is,0) = 0.0_R8
            end do
          else
            do is = 1, natmi
              call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,is,0))
            end do
          end if
          allocate(eneutrad1(0:nx1+1,0:ny1+1,natm1))
          do is = 1, natmi
            if (iatm1(is).ne.0) then
              field(0:nx+1,0:ny+1) = eneutrad(0:nx+1,0:ny+1,is,0)
              call b2xtvc(nx, ny, nx1, ny1,
     &         field, eneutrad1(0,0,iatm1(is)),
     &         nncutmax, non_commensurate,
     &         nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &         leftcut, leftcut1, rightcut, rightcut1,
     &         bottomcut, bottomcut1, topcut, topcut1,
     &         wrk0)
            end if
          end do
          call write_title(45,'eneutrad',n2*natm1)
          call eirene_neutr(45,nx1,ny1,natm1,eneutrad1,
     .                         nx1,ny1,natm1,1,1)
          deallocate(eneutrad1)
        end if

        if(jvft44 .ge. 20080706) then
          call read_title(44)
          READ (44,'(12i6)') (eirdiag_nds_ind(i),i=1,nstsi+1)        
          READ (44,'(12i6)') (eirdiag_nds_typ(i),i=1,nstsi)    
          READ (44,'(12i6)') (eirdiag_nds_srf(i),i=1,nstsi)    
          READ (44,'(12i6)') (eirdiag_nds_start(i),i=1,nstsi)    
          READ (44,'(12i6)') (eirdiag_nds_end(i),i=1,nstsi)
          allocate(eirdiag_nds_ind1(nstsi+1))
          allocate(eirdiag_nds_srf1(nstsi))
          allocate(eirdiag_nds_start1(nstsi))
          allocate(eirdiag_nds_end1(nstsi))
          eirdiag_nds_ind1(1) = 1
          do i = 1, nstsi
            if (eirdiag_nds_typ(i).eq.1) then
              eirdiag_nds_srf1(i) =
     .         conformal_convert(non_commensurate,2,
     .         (eirdiag_nds_start(i)+eirdiag_nds_end(i))/2,
     .          eirdiag_nds_srf(i),nx,ny,nx1,ny1,4)
              eirdiag_nds_start1(i) =
     .         convert_eirdiag(non_commensurate,1,
     .          eirdiag_nds_start(i),eirdiag_nds_srf(i),nx,ny,nx1,ny1,1)
              eirdiag_nds_end1(i) =
     .         convert_eirdiag(non_commensurate,1,
     .          eirdiag_nds_end(i),eirdiag_nds_srf(i),nx,ny,nx1,ny1,3)
              eirdiag_nds_ind1(i+1) = eirdiag_nds_ind1(i) +
     .         (eirdiag_nds_end1(i) - eirdiag_nds_start1(i) + 1) + 1
            else if (eirdiag_nds_typ(i).eq.2) then
              eirdiag_nds_srf1(i) =
     .         conformal_convert(non_commensurate,1,
     .          eirdiag_nds_srf(i),
     .         (eirdiag_nds_start(i)+eirdiag_nds_end(i))/2,
     .          nx,ny,nx1,ny1,4)
              eirdiag_nds_start1(i) =
     .         convert_eirdiag(non_commensurate,2,
     .          eirdiag_nds_srf(i),eirdiag_nds_start(i),nx,ny,nx1,ny1,0)
              eirdiag_nds_end1(i) =
     .         convert_eirdiag(non_commensurate,2,
     .          eirdiag_nds_srf(i),eirdiag_nds_end(i),nx,ny,nx1,ny1,2)
              eirdiag_nds_ind1(i+1) = eirdiag_nds_ind1(i) +
     .         (eirdiag_nds_end1(i) - eirdiag_nds_start1(i) + 1) + 1
            else if (eirdiag_nds_typ(i).eq.0) then
              eirdiag_nds_ind1(i+1) = eirdiag_nds_ind1(i)
              eirdiag_nds_ind1(i) = -1
              eirdiag_nds_srf1(i) = 0
              eirdiag_nds_start1(i) = 0
              eirdiag_nds_end1(i) = 0
            end if
          end do
          call write_title(45,'eirdiag', 5*nstsi+1)
          WRITE(45,'(12i6)') (eirdiag_nds_ind1(i),i=1,nstsi+1)
          WRITE(45,'(12i6)') (eirdiag_nds_typ(i),i=1,nstsi)
          WRITE(45,'(12i6)') (eirdiag_nds_srf1(i),i=1,nstsi)
          WRITE(45,'(12i6)') (eirdiag_nds_start1(i),i=1,nstsi)
          WRITE(45,'(12i6)') (eirdiag_nds_end1(i),i=1,nstsi)
      
          NCL=eirdiag_nds_ind(nstsi+1)
          NCL1=eirdiag_nds_ind1(nstsi+1)
          NCLM=max(NCL,NCL1)
          IF(NCL.GT.0) THEN
            IF(.NOT.ALLOCATED(wldna_res))
     >       ALLOCATE(wldna_res(natmi,NCL), wldnm_res(nmoli,NCL),
     .                ewlda_res(natmi,NCL), ewldm_res(nmoli,NCL), 
     .                ewldea_res(natmi,NCL), ewldem_res(nmoli,NCL), 
     .                ewldrp_res(NCL), ewldmr_res(nmoli,NCL),
     .                wldspt_res(NCL), sarea_res(NCL),
     .                wldspta_res(natmi,NCL), wldsptm_res(nmoli,NCL),
     .                wlpump_res(nspz,NCL) )
 
            call read_title(44)
            READ (44,'(6e15.7)') (sarea_res(i),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((wldna_res(j,i),j=1,natmi),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((wldnm_res(j,i),j=1,nmoli),i=1,NCL) 
            call read_title(44)
            READ (44,'(6e15.7)') ((ewlda_res(j,i),j=1,natmi),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((ewldm_res(j,i),j=1,nmoli),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((ewldea_res(j,i),j=1,natmi),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((ewldem_res(j,i),j=1,nmoli),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') (ewldrp_res(i),i=1,NCL)
            call read_title(44)
            READ (44,'(6e15.7)') ((ewldmr_res(j,i),j=1,nmoli),i=1,NCL)
            if (jvft44.ge.20170228) then
              call read_title(44)
              READ (44,'(6e15.7)') (wldspt_res(i),i=1,NCL)
            end if
            if (jvft44.ge.20170328) then
              call read_title(44)
              READ (44,'(6e15.7)')
     .              ((wldspta_res(j,i),j=1,natmi),i=1,NCL)
              call read_title(44)
              READ (44,'(6e15.7)')
     .              ((wldsptm_res(j,i),j=1,nmoli),i=1,NCL)
              k=0
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
              read (44,'(6e13.5)') ((wlpump_res(j+k,i),j=1,natmi),
     ,                                i=1,NCL)
              k=k+natmi
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
              read (44,'(6e13.5)') ((wlpump_res(j+k,i),j=1,nmoli),
     ,                                i=1,NCL)
              k=k+nmoli
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
              read (44,'(6e13.5)') ((wlpump_res(j+k,i),j=1,nioni),
     ,                                i=1,NCL)
              k=k+nioni
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
              read (44,'(6e13.5)') ((wlpump_res(j+k,i),j=1,nplsi),
     ,                                i=1,NCL)
            end if
            allocate(sarea1res(NCL1))
            allocate(wldna1res(natm1,NCL1), ewlda1res(natm1,NCL1),
     .               ewldea1res(natm1,NCL1))
            allocate(wldnm1res(nmol1,NCL1), ewldm1res(nmol1,NCL1),
     .               ewldem1res(nmol1,NCL1), ewldmr1res(nmol1,NCL1))
            allocate(ewldrp1res(NCL1))
            allocate(wldspt1res(NCL1))
            allocate(data1res(NCLM))
            if (jvftnew .ge. 20170328)
     >       allocate(wldspta1res(natm1,NCL1),wldsptm1res(nmol1,NCL1))
            allocate(area_res_contr(NCL,NCL1))
            call b2yt_area_res_map ( NCL, NCL1, nstsi, nx, ny, nx1, ny1,
     .       non_commensurate, sarea_res, sarea1res, area_res_contr,
     .       eirdiag_nds_ind, eirdiag_nds_typ, eirdiag_nds_srf, 
     .       eirdiag_nds_start, eirdiag_nds_end,
     .       eirdiag_nds_ind1, eirdiag_nds_srf1, 
     .       eirdiag_nds_start1, eirdiag_nds_end1,
     .       gs1 )
            data1res(1:NCL) = ewldrp_res(1:NCL)
            call b2yt_res_map ( NCL, NCL1, NCLM,
     .       area_res_contr, data1res)
            ewldrp1res(1:NCL1) = data1res(1:NCL1)
            if (no_sputtering) then
              wldspt1res(:) = 0.0_R8
            else
              data1res(1:NCL) = wldspt_res(1:NCL)
              call b2yt_res_map ( NCL, NCL1, NCLM,
     .         area_res_contr, data1res) 
              wldspt1res(1:NCL1) = data1res(1:NCL1)
            end if
            do is = 1, natmi
              if (iatm1(is).ne.0) then
                data1res(1:NCL) = wldna_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                wldna1res(iatm1(is),1:NCL1) = data1res(1:NCL1)
                data1res(1:NCL) = ewlda_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                ewlda1res(iatm1(is),1:NCL1) = data1res(1:NCL1)
                data1res(1:NCL) = ewldea_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                ewldea1res(iatm1(is),1:NCL1) = data1res(1:NCL1)
                if (no_sputtering) then
                  wldspta1res(iatm1(is),:) = 0.0_R8
                else
                  data1res(1:NCL) = wldspta_res(is,1:NCL)
                  call b2yt_res_map ( NCL, NCL1, NCLM,
     .             area_res_contr, data1res) 
                  wldspta1res(iatm1(is),1:NCL1) = data1res(1:NCL1)
                end if
              end if
            end do
            do is = 1, nmoli
              if (imol1(is).ne.0) then
                data1res(1:NCL) = wldnm_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                wldnm1res(imol1(is),1:NCL1) = data1res(1:NCL1)
                data1res(1:NCL) = ewldm_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                ewldm1res(imol1(is),1:NCL1) = data1res(1:NCL1)
                data1res(1:NCL) = ewldem_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                ewldem1res(imol1(is),1:NCL1) = data1res(1:NCL1)
                data1res(1:NCL) = ewldmr_res(is,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                ewldmr1res(imol1(is),1:NCL1) = data1res(1:NCL1)
                if (no_sputtering) then
                  wldsptm1res(imol1(is),:) = 0.0_R8
                else
                  data1res(1:NCL) = wldsptm_res(is,1:NCL)
                  call b2yt_res_map ( NCL, NCL1, NCLM,
     .             area_res_contr, data1res) 
                  wldsptm1res(imol1(is),1:NCL1) = data1res(1:NCL1)
                end if
              end if
            end do
            call write_title(45,'sarea_res',ncl1)
            WRITE(45,'(6e15.7)') (sarea1res(i),i=1,NCL1)
            call write_title(45,'wldna_res',ncl1*natm1)
            WRITE(45,'(6e15.7)') ((wldna1res(j,i),j=1,natm1),i=1,NCL1)
            call write_title(45,'wldnm_res',ncl1*nmol1)
            WRITE(45,'(6e15.7)') ((wldnm1res(j,i),j=1,nmol1),i=1,NCL1)
            call write_title(45,'ewlda_res',ncl1*natm1)
            WRITE(45,'(6e15.7)') ((ewlda1res(j,i),j=1,natm1),i=1,NCL1)
            call write_title(45,'ewldm_res',ncl1*nmol1)
            WRITE(45,'(6e15.7)') ((ewldm1res(j,i),j=1,nmol1),i=1,NCL1)
            call write_title(45,'ewldea_res',ncl1*natm1)
            WRITE(45,'(6e15.7)') ((ewldea1res(j,i),j=1,natm1),i=1,NCL1)
            call write_title(45,'ewldem_res',ncl1*nmol1)
            WRITE(45,'(6e15.7)') ((ewldem1res(j,i),j=1,nmol1),i=1,NCL1)
            call write_title(45,'ewldrp_res',ncl1)
            WRITE(45,'(6e15.7)') (ewldrp1res(i),i=1,NCL1)
            call write_title(45,'ewldmr_res',ncl1*nmol1)
            WRITE(45,'(6e15.7)') ((ewldmr1res(j,i),j=1,nmol1),i=1,NCL1)
            if (jvft44.ge.20170228) then
              call write_title(45,'wldspt_res',ncl1)
              WRITE(45,'(6e15.7)') (wldspt1res(i),i=1,NCL1)
            end if
            if (jvft44.ge.20170328) then
              call write_title(45,'wldspta_res',ncl1*natm1)
              WRITE(45,'(6e15.7)')
     .              ((wldspta1res(j,i),j=1,natm1),i=1,NCL1)
              call write_title(45,'wldsptm_res',ncl1*nmol1)
              WRITE(45,'(6e15.7)')
     .              ((wldsptm1res(j,i),j=1,nmol1),i=1,NCL1)
              allocate (wlpump1res(natm1+nmol1+nioni+npls1,NCLM))
              do k = 1, natmi
                if (iatm1(k).ne.0) then
                  data1res(1:NCL) = wlpump_res(k,1:NCL)
                  call b2yt_res_map ( NCL, NCL1, NCLM,
     .             area_res_contr, data1res) 
                  wlpump1res(iatm1(k),1:NCL1) = data1res(1:NCL1)
                end if
              end do
              do k = 1, nmoli
                if (imol1(k).ne.0) then
                  data1res(1:NCL) = wlpump_res(natmi+k,1:NCL)
                  call b2yt_res_map ( NCL, NCL1, NCLM,
     .             area_res_contr, data1res) 
                  wlpump1res(natm1+imol1(k),1:NCL1) = data1res(1:NCL1)
                end if
              end do
              do k = 1, nioni
                data1res(1:NCL) = wlpump_res(natmi+nmoli+k,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                wlpump1res(natm1+nmol1+k,1:NCL1) = data1res(1:NCL1)
              end do
              do i = 1, NCL
                wd_old(1:nplsi)=wlpump_res(natmi+nmoli+nioni+1:
     .                                     natmi+nmoli+nioni+nplsi,i)
                call e_bundle_transmit_r(nplsi,npls1,e_bundle_contr,
     .           wd_old,wd_new,2)
                wlpump1res(natm1+nmol1+nioni+1:
     .                     natm1+nmol1+nioni+npls1,i)=
     .           wd_new(1:npls1)
              end do
              do k = 1, npls1
                data1res(1:NCL) = wlpump1res(natm1+nmol1+nioni+k,1:NCL)
                call b2yt_res_map ( NCL, NCL1, NCLM,
     .           area_res_contr, data1res) 
                wlpump1res(natm1+nmol1+nioni+k,1:NCL1) =
     .           data1res(1:NCL1)
              end do
              k=0
              call write_title(45,'wlpump_res(A)',ncl1*natm1)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,natm1)
              write(45,'(1p,6e13.5)') ((wlpump1res(j+k,i),j=1,natm1),
     ,                                                       i=1,ncl1)
              k=k+natm1
              call write_title(45,'wlpump_res(M)',ncl1*nmol1)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,nmol1)
              write(45,'(1p,6e13.5)') ((wlpump1res(j+k,i),j=1,nmol1),
     ,                                                       i=1,ncl1)
              k=k+nmol1
              call write_title(45,'wlpump_res(I)',ncl1*nioni)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,nioni)
              write(45,'(1p,6e13.5)') ((wlpump1res(j+k,i),j=1,nioni),
     ,                                                       i=1,ncl1)
              k=k+nioni
              call write_title(45,'wlpump_res(P)',ncl1*npls1)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,npls1)
              write(45,'(1p,6e13.5)') ((wlpump1res(j+k,i),j=1,npls1),
     ,                                                       i=1,ncl1)
              deallocate (wlpump1res)
            end if
            deallocate(sarea1res)
            deallocate(wldna1res, ewlda1res, ewldea1res)
            deallocate(wldnm1res, ewldm1res, ewldem1res)
            deallocate(ewldmr1res)
            deallocate(ewldrp1res)
            deallocate(wldspt1res)
            deallocate(data1res)
            deallocate(area_res_contr)
            if (jvftnew .ge. 20170328)
     >       deallocate(wldspta1res,wldsptm1res)
          ELSE
            IS = 9
            if (jvft44.ge.20170228) IS = IS + 1
            if (jvft44.ge.20170328) IS = IS + 2
            DO I = 1, IS
              call read_title(44)
              READ (44,1028)
              WRITE(45,1028)
            END DO
            if (jvft44.ge.20170328) THEN
              k=0
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,natm1)
              READ (44,1028)
              WRITE(45,1028)
              k=k+natmi
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,nmol1)
              READ (44,1028)
              WRITE(45,1028)
              k=k+nmoli
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,nioni)
              READ (44,1028)
              WRITE(45,1028)
              k=k+nioni
              call read_title(44)
              read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
              write(45,'(6a13)') (eirtxt1(j+k),j=1,npls1)
              READ (44,1028)
              WRITE(45,1028)
            END IF
          END IF
          deallocate (eirtxt1)
          deallocate (wd_old, wd_new)
          deallocate (eirdiag_nds_ind1,eirdiag_nds_srf1)
          deallocate (eirdiag_nds_start1,eirdiag_nds_end1)
        end if

        if(streql(solps_version,'4.3')) goto 99
        if(jvft44 .ge. 20081111) then
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENA_INT(I,J),I=1,natmi),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENM_INT(I,J),I=1,nmoli),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENI_INT(I,J),I=1,nioni),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENA_INT_B2(I,J),I=1,natmi),
     .                                                     J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENM_INT_B2(I,J),I=1,nmoli),
     .                                                     J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENI_INT_B2(I,J),I=1,nioni),
     .                                                     J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENA_INT(I,J),I=1,natmi),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENM_INT(I,J),I=1,nmoli),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENI_INT(I,J),I=1,nioni),J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENA_INT_B2(I,J),I=1,natmi),
     .                                                     J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENM_INT_B2(I,J),I=1,nmoli),
     .                                                     J=0,NSTRAI)
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENI_INT_B2(I,J),I=1,nioni),
     .                                                     J=0,NSTRAI)
          allocate (
     *     pdena1int(natm1,0:nstrai), pdenm1int(nmol1,0:nstrai),
     *     pdeni1int(nioni,0:nstrai),
     *     pdena1int_b2(natm1,0:nstrai), pdenm1int_b2(nmol1,0:nstrai),
     *     pdeni1int_b2(nioni,0:nstrai),
     *     edena1int(natm1,0:nstrai), edenm1int(nmol1,0:nstrai),
     *     edeni1int(nioni,0:nstrai),
     *     edena1int_b2(natm1,0:nstrai), edenm1int_b2(nmol1,0:nstrai),
     *     edeni1int_b2(nioni,0:nstrai))
          istra = -1
          do j = 0, nstrai
            if (savestratum(j)) then
              istra = istra+1
              do is = 1, natmi
                if (iatm1(is).ne.0) then
                  pdena1int(iatm1(is),istra) = pdena_int(is,j)
                  pdena1int_b2(iatm1(is),istra) = pdena_int_b2(is,j)
                  edena1int(iatm1(is),istra) = edena_int(is,j)
                  edena1int_b2(iatm1(is),istra) = edena_int_b2(is,j)
                end if
              end do
              do is = 1, nmoli
                if (imol1(is).ne.0) then
                  pdenm1int(imol1(is),istra) = pdenm_int(is,j)
                  pdenm1int_b2(imol1(is),istra) = pdenm_int_b2(is,j)
                  edenm1int(imol1(is),istra) = edenm_int(is,j)
                  edenm1int_b2(imol1(is),istra) = edenm_int_b2(is,j)
                end if
              end do
              pdeni1int(1:nioni,istra) = pdeni_int(1:nioni,j)
              pdeni1int_b2(1:nioni,istra) = pdeni_int_b2(1:nioni,j)
              edeni1int(1:nioni,istra) = edeni_int(1:nioni,j)
              edeni1int_b2(1:nioni,istra) = edeni_int_b2(1:nioni,j)
            end if
          end do
          call write_title(45,'pdena_int',(nstra1+1)*natm1)
          WRITE(45,'(6e15.7)') ((PDENA1INT(I,J),I=1,natm1),J=0,NSTRA1)
          call write_title(45,'pdenm_int',(nstra1+1)*nmol1)
          WRITE(45,'(6e15.7)') ((PDENM1INT(I,J),I=1,nmol1),J=0,NSTRA1)
          call write_title(45,'pdeni_int',(nstra1+1)*nioni)
          WRITE(45,'(6e15.7)') ((PDENI1INT(I,J),I=1,nioni),J=0,NSTRA1)
          call write_title(45,'pdena_int_b2',(nstra1+1)*natm1)
          WRITE(45,'(6e15.7)') ((PDENA1INT_B2(I,J),I=1,natm1),
     .                                                     J=0,NSTRA1)
          call write_title(45,'pdenm_int_b2',(nstra1+1)*nmol1)
          WRITE(45,'(6e15.7)') ((PDENM1INT_B2(I,J),I=1,nmol1),
     .                                                     J=0,NSTRA1)
          call write_title(45,'pdeni_int_b2',(nstra1+1)*nioni)
          WRITE(45,'(6e15.7)') ((PDENI1INT_B2(I,J),I=1,nioni),
     .                                                     J=0,NSTRA1)
          call write_title(45,'edena_int',(nstra1+1)*natm1)
          WRITE(45,'(6e15.7)') ((EDENA1INT(I,J),I=1,natm1),J=0,NSTRA1)
          call write_title(45,'edenm_int',(nstra1+1)*nmol1)
          WRITE(45,'(6e15.7)') ((EDENM1INT(I,J),I=1,nmol1),J=0,NSTRA1)
          call write_title(45,'edeni_int',(nstra1+1)*nioni)
          WRITE(45,'(6e15.7)') ((EDENI1INT(I,J),I=1,nioni),J=0,NSTRA1)
          call write_title(45,'edena_int_b2',(nstra1+1)*natm1)
          WRITE(45,'(6e15.7)') ((EDENA1INT_B2(I,J),I=1,natm1),
     .                                                     J=0,NSTRA1)
          call write_title(45,'edenm_int_b2',(nstra1+1)*nmol1)
          WRITE(45,'(6e15.7)') ((EDENM1INT_B2(I,J),I=1,nmol1),
     .                                                     J=0,NSTRA1)
          call write_title(45,'edeni_int_b2',(nstra1+1)*nioni)
          WRITE(45,'(6e15.7)') ((EDENI1INT_B2(I,J),I=1,nioni),
     .                                                     J=0,NSTRA1)
          deallocate (
     *     pdena1int, pdenm1int, pdeni1int,
     *     pdena1int_b2, pdenm1int_b2, pdeni1int_b2,
     *     edena1int, edenm1int, edeni1int,
     *     edena1int_b2, edenm1int_b2, edeni1int_b2)
        end if
      else  ! NLWRMSH
        if(streql(solps_version,'5.1')) then
          if (jvft44.ge.20130210) then
            call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,1,0))
#ifndef HPUX
            write(6,'(2a)') 'In this fort.44 version,',
     .                      ' the eneutrad data is summed over species!'
#else
            write(6,'(2a)') 'In this ftn44 version,',
     .                      ' the eneutrad data is summed over species!'
#endif
            do is = 2, natmi
              eneutrad(:,:,is,0) = 0.0_R8
            end do
            allocate(eneutrad1(0:nx1+1,0:ny1+1,natm1))
            do is = 1, natmi
              if (iatm1(is).ne.0) then
                field(0:nx+1,0:ny+1) = eneutrad(0:nx+1,0:ny+1,is,0)
                call b2xtvc(nx, ny, nx1, ny1,
     &           field, eneutrad1(0,0,iatm1(is)),
     &           nncutmax, non_commensurate,
     &           nnreg, nncut, periodic_bc,
     &           jsep, jsep1, region, region1,
     &           leftcut, leftcut1, rightcut, rightcut1,
     &           bottomcut, bottomcut1, topcut, topcut1,
     &           wrk0)
              end if
            end do
            call write_title(45,'eneutrad',n2*natm1)
            call eirene_neutr(45,nx1,ny1,natm1,eneutrad1,
     .                           nx1,ny1,natm1,1,1)
            deallocate(eneutrad1)
          end if
        elseif(streql(solps_version,'4.3')) then
          if (jvft44.ge.20071209) then
            do is = 1, natmi
              call gfsub2(44,ndx,ndy,nnx,nny,eneutrad(:,:,is,0))
            end do
            allocate(eneutrad1(0:nx1+1,0:ny1+1,natm1))
            do is = 1, natmi
              if (iatm1(is).ne.0) then
                field(0:nx+1,0:ny+1) = eneutrad(0:nx+1,0:ny+1,is,0)
                call b2xtvc(nx, ny, nx1, ny1,
     &           field, eneutrad1(0,0,iatm1(is)),
     &           nncutmax, non_commensurate,
     &           nnreg, nncut, periodic_bc, 
     &           jsep, jsep1, region, region1,
     &           leftcut, leftcut1, rightcut, rightcut1,
     &           bottomcut, bottomcut1, topcut, topcut1,
     &           wrk0)
              end if
            end do
            call write_title(45,'eneutrad',n2*natm1)
            call eirene_neutr(45,nx1,ny1,natm1,eneutrad1,
     .                           nx1,ny1,natm1,1,1)
            deallocate(eneutrad1)
          end if
        end if
      end if
      if(jvft44.ge.960513) nlimi = nlim1

#endif

*  ---------------------------------------------------------------------
   99 CONTINUE
      close(44)
      close(45)
 1028 format()
C
      call subend ()
      RETURN

*=======================================================================
*//END B2YT_NGREAD//

      END

      subroutine b2yt_area_res_map ( NCL, NCL1, nstsi, nx, ny, nx1, ny1,
     .  non_commensurate, sarea_res, sarea1res, area_res_contr,
     .  eirdiag_nds_ind, eirdiag_nds_typ, eirdiag_nds_srf, 
     .  eirdiag_nds_start, eirdiag_nds_end,
     .  eirdiag_nds_ind1, eirdiag_nds_srf1, 
     .  eirdiag_nds_start1, eirdiag_nds_end1,
     .  gs1)
      use b2mod_types
      implicit none
      integer nstsi, nx, ny, nx1, ny1, non_commensurate
      integer NCL, NCL1, NCLM
      real (kind=R8) :: sarea_res( NCL )
      real (kind=R8) :: sarea1res( NCL1 )
      real (kind=R8) :: area_res_contr( NCL, NCL1 )
      integer :: eirdiag_nds_typ(nstsi)
      integer :: eirdiag_nds_ind(nstsi+1), eirdiag_nds_ind1(nstsi+1)
      integer :: eirdiag_nds_srf(nstsi), eirdiag_nds_srf1(nstsi)
      integer :: eirdiag_nds_start(nstsi), eirdiag_nds_start1(nstsi)
      integer :: eirdiag_nds_end(nstsi), eirdiag_nds_end1(nstsi)
      real (kind=R8) :: gs1(-1:nx1,-1:ny1,0:2)
      real (kind=R8) :: frac(NCL), frac1(NCL1)
      real (kind=R8) :: lb, ub, lb1, ub1
      integer :: i, j, ind, ind1, ix, iy, stride

      real (kind=R8) :: rratio
      external rratio
      
      sarea1res(:) = 0.0_R8
      area_res_contr(:,:) = 0.0_R8
      do i = 1, nstsi
        ind = eirdiag_nds_ind(i)
        if (ind.le.0) cycle
        ind1 = eirdiag_nds_ind1(i)
        area_res_contr(ind,ind1) = 1.0_R8
        frac(ind) = 1.0_R8
        frac1(ind1) = 1.0_R8
        do j = eirdiag_nds_start1(i), eirdiag_nds_end1(i)
          ind1 = eirdiag_nds_ind1(i) + j - eirdiag_nds_start1(i) + 1
          if (eirdiag_nds_typ(i).eq.1) then
            ix = j - 1
            iy = eirdiag_nds_srf1(i)
            sarea1res(ind1) = gs1(ix,iy,1)
          else if (eirdiag_nds_typ(i).eq.2) then
            ix = eirdiag_nds_srf1(i)
            iy = j - 1
            sarea1res(ind1) = gs1(ix,iy,0)
          end if
          sarea1res(eirdiag_nds_ind1(i)) =
     .     sarea1res(eirdiag_nds_ind1(i)) + sarea1res(ind1)
        end do
cxpb frac contains the accumulated fractional area along the original surface
        do j = eirdiag_nds_start(i), eirdiag_nds_end(i)
          ind = eirdiag_nds_ind(i) + j - eirdiag_nds_start(i) + 1
          if (j.eq.eirdiag_nds_start(i)) then
            frac(ind) = sarea_res(ind)/sarea_res(eirdiag_nds_ind(i))
          else
            frac(ind) = frac(ind-1) +
     .       sarea_res(ind)/sarea_res(eirdiag_nds_ind(i))
          end if
        end do
cxpb frac1 contains the accumulated fractional area along the mapped surface
        do j = eirdiag_nds_start1(i), eirdiag_nds_end1(i)
          ind1 = eirdiag_nds_ind1(i) + j - eirdiag_nds_start1(i) + 1
          if (j.eq.eirdiag_nds_start1(i)) then
            frac1(ind1) = sarea1res(ind1)/sarea1res(eirdiag_nds_ind1(i))
          else
            frac1(ind1) = frac1(ind1-1) +
     .       sarea1res(ind1)/sarea1res(eirdiag_nds_ind1(i))
          end if
        end do
cxpb The area_res_contr(ind,ind1) array contains the fraction of the surface area from
cxpb original element ind that is to be found as part of mapped element ind1
cxpb lb and ub are lower and upper interval bounds
        if (non_commensurate.eq.0) then
          if (eirdiag_nds_typ(i).eq.1) then
            if (nx.eq.nx1) then
              do ind = eirdiag_nds_ind(i) + 1,
     .                 eirdiag_nds_ind(i) + eirdiag_nds_end(i) -
     .                                      eirdiag_nds_start(i) + 1
                ind1 = eirdiag_nds_ind1(i) + ind - eirdiag_nds_ind(i)
                area_res_contr(ind,ind1) = 1.0_R8
              end do
            else if (nx.gt.nx1) then
              stride = nx/nx1
              ind1 = eirdiag_nds_ind1(i) + 1
              do ind = eirdiag_nds_ind(i) + 1,
     .                 eirdiag_nds_ind(i) + eirdiag_nds_end(i) -
     .                 eirdiag_nds_start(i) + 1
                area_res_contr(ind,ind1) = 1.0_R8
                if (mod(ind-eirdiag_nds_ind(i),stride).eq.0)
     .           ind1 = ind1 + 1
              end do  
            else if (nx1.gt.nx) then
              stride = nx1/nx
              ind = eirdiag_nds_ind(i) + 1
              do ind1 = eirdiag_nds_ind1(i) + 1,
     .                  eirdiag_nds_ind1(i) + eirdiag_nds_end1(i) -
     .                  eirdiag_nds_start1(i) + 1
                area_res_contr(ind,ind1) = rratio(nx,nx1)
                if (mod(ind1-eirdiag_nds_ind1(i),stride).eq.0)
     .           ind = ind + 1
              end do
            end if
          else if (eirdiag_nds_typ(i).eq.2) then
            if (ny.eq.ny1) then
              do ind = eirdiag_nds_ind(i) + 1,
     .                 eirdiag_nds_ind(i) + eirdiag_nds_end(i) -
     .                 eirdiag_nds_start(i) + 1
                ind1 = eirdiag_nds_ind1(i) + ind - eirdiag_nds_ind(i) 
                area_res_contr(ind,ind1) = 1.0_R8
              end do
            else if (ny.gt.ny1) then
              stride = ny/ny1
              ind1 = eirdiag_nds_ind1(i) + 1
              do ind = eirdiag_nds_ind(i) + 1,
     .                 eirdiag_nds_ind(i) + eirdiag_nds_end(i) -
     .                 eirdiag_nds_start(i) + 1
                area_res_contr(ind,ind1) = 1.0_R8
                if (mod(ind-eirdiag_nds_ind(i),stride).eq.0)
     .           ind1 = ind1 + 1
              end do  
            else if (ny1.gt.ny) then
              stride = ny1/ny
              ind = eirdiag_nds_ind(i) + 1
              do ind1 = eirdiag_nds_ind1(i) + 1,
     .                  eirdiag_nds_ind1(i) + eirdiag_nds_end1(i) -
     .                  eirdiag_nds_start1(i) + 1
                area_res_contr(ind,ind1) = rratio(nx,nx1)
                if (mod(ind1-eirdiag_nds_ind1(i),stride).eq.0)
     .           ind = ind + 1
              end do
            end if
          end if
        else
          ind = eirdiag_nds_ind(i)+1
          lb = 0.0_R8
          ind1 = eirdiag_nds_ind1(i)+1
          lb1 = 0.0_R8
          do while (ind.le.eirdiag_nds_ind(i)+1+
     .                     eirdiag_nds_end(i)-eirdiag_nds_start(i) .and.
     .              ind1.le.eirdiag_nds_ind1(i)+1+
     .                      eirdiag_nds_end1(i)-eirdiag_nds_start1(i))
            ub = frac(ind)
            ub1 = frac1(ind1)
            if (lb.le.lb1 .and. ub.lt.ub1) then
              area_res_contr(ind,ind1) = (ub-lb1)/(ub-lb)
              ind = ind+1
              lb = ub
            else if (lb.le.lb1 .and. ub.gt.ub1) then
              area_res_contr(ind,ind1) = (ub1-lb1)/(ub-lb)
              ind1 = ind1 + 1
              lb1 = ub1
            else if (lb.le.lb1 .and. ub.eq.ub1) then
              area_res_contr(ind,ind1) = (ub-lb1)/(ub-lb)
              ind = ind + 1
              lb = ub
              ind1 = ind1 + 1
              lb1 = ub1
            else if (lb.gt.lb1 .and. ub.lt.ub1) then
              area_res_contr(ind,ind1) = 1.0_R8
              ind = ind + 1
              lb = ub
            else if (lb.gt.lb1 .and. ub.gt.ub1) then
              area_res_contr(ind,ind1) = (ub1-lb)/(ub-lb)
              ind1 = ind1 + 1
              lb1 = ub1
            else if (lb.gt.lb1 .and. ub.eq.ub1) then
              area_res_contr(ind,ind1) = 1.0_R8
              ind = ind + 1
              lb = ub
              ind1 = ind1 + 1
              lb1 = ub1
            end if
          end do
        end if       
      end do

      return
      end subroutine b2yt_area_res_map

      subroutine b2yt_res_map ( NCL, NCL1, NCLM,
     . area_res_contr, data1res )
      use b2mod_types
      implicit none
      integer :: NCL, NCL1, NCLM
      real (kind=R8) :: area_res_contr(NCL,NCL1)
      real (kind=R8) :: data1res(NCLM)
      real (kind=R8) :: source_res(NCL)
      real (kind=R8) :: mapped_res(NCL1)
      integer :: ind, ind1

      source_res(1:NCL) = data1res(1:NCL)
      mapped_res = 0.0_R8
      do ind1 = 1, NCL1
        do ind = 1, NCL
          mapped_res(ind1) = mapped_res(ind1) +
     .     source_res(ind)*area_res_contr(ind,ind1)
        end do
      end do
      data1res(1:NCL1) = mapped_res(1:NCL1)

      return
      end subroutine b2yt_res_map

      integer function convert_eirdiag
     . (non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside)
      use b2mod_geo2
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer non_commensurate,xory,ix,iy,nx,ny,nx1,ny1,iside
c xory  = 1 : x-direction
c       = 2 : y-direction
c iside = 0 : bottom
c       = 1 : left
c       = 2 : top
c       = 3 : right
c       = 4 : center
*   ..local variables
      integer ixpos, iypos, iplg, nsize, nsize1
      integer resleft(6), resright(6)
      integer resleft1(6), resright1(6)
      integer nxtl, nxtr
      integer nxtl1, nxtr1
      logical atend
*   ..procedures
      integer conformal_convert, convert_position
      external conformal_convert, convert_position, get_nxt, b2yt_nxt,
     .         xerrab

      ixpos = ix
      resleft = nx+2
      resright= nx+2
      if (nncut.eq.0) then
        resleft(1) = 1
        resright(1)= nx+1
        resleft1(1) = 1
        resright1(1)= nx1+1
      else if (nncut.eq.1) then
        resleft(1) = 1
        resright(1)= leftcut(1)
        resleft(2) = leftcut(1)+1
        resright(2)= rightcut(1)
        resleft(3) = rightcut(1)+1
        resright(3)= nx+1
        resleft1(1) = 1
        resright1(1)= leftcut1(1)
        resleft1(2) = leftcut1(1)+1
        resright1(2)= rightcut1(1)
        resleft1(3) = rightcut1(1)+1
        resright1(3)= nx1+1
      else if (nncut.eq.2) then
        call get_nxt(nx, nxtl, nxtr)
        call b2yt_nxt(nx1, nxtl1, nxtr1)
        resleft(1) = 1
        resright(1)= leftcut(1)
        resleft(2) = leftcut(1)+1
        resright(2)= leftcut(2)
        resleft(3) = leftcut(2)+1
        resright(3)= nxtl
        resleft(4) = nxtr+1
        resright(4)= rightcut(2)
        resleft(5) = rightcut(2)+1
        resright(5)= rightcut(1)
        resleft(6) = rightcut(1)+1
        resright(6)= nx+1
        resleft1(1) = 1
        resright1(1)= leftcut1(1)
        resleft1(2) = leftcut1(1)+1
        resright1(2)= leftcut1(2)
        resleft1(3) = leftcut1(2)+1
        resright1(3)= nxtl1
        resleft1(4) = nxtr1+1
        resright1(4)= rightcut1(2)
        resleft1(5) = rightcut1(2)+1
        resright1(5)= rightcut1(1)
        resleft1(6) = rightcut1(1)+1
        resright1(6)= nx1+1
      endif
      iplg = 1
      do while (ix.gt.resright(iplg).and.iplg.le.5)
        iplg = iplg+1
      enddo
      if (xory.eq.2) then
        if (iy.eq.0) then
          convert_eirdiag = 0
        else if (iy.gt.ny+1) then
          call xerrab ('Too large IY value in convert_eirdiag')
          convert_eirdiag = 0
        else
          iypos = iy - 1
          convert_eirdiag = 1 + conformal_convert
     .     (non_commensurate,2,ixpos,iypos,nx,ny,nx1,ny1,iside)
        end if
      else if (xory.eq.1) then
        if (ix.eq.0) then
          convert_eirdiag = 0
        elseif (ix.eq.resleft(iplg)) then
          convert_eirdiag = resleft1(iplg)
        elseif (ix.eq.resright(iplg)) then
          convert_eirdiag = resright1(iplg)
        elseif (ix.gt.nx+1) then
          call xerrab ('Too large IX value in convert_eirdiag')
          convert_eirdiag = 0
        else
          ixpos = ix - resleft(iplg)
          nsize = resright(iplg) - resleft(iplg)
          nsize1= resright1(iplg)- resleft1(iplg)
          atend = iside.eq.3
          convert_eirdiag = resleft1(iplg) +
     .     convert_position(ixpos,nsize,nsize1,0,atend)
        endif
      else
        call xerrab ('In convert_eirdiag = 0 : xory test failed')
        convert_eirdiag = 0
      endif

      return
      end function convert_eirdiag




