c
C RHW 18.08.92: UEBERARBEITETE VERSION VON RFS UEBERNOMMEN.
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

      SUBROUTINE NGREAD (mpg,geo,ns,switch,ismain,NFEND)
c
c  version : 06.11.97 22:11
c
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_geometry
      use b2mod_b2cmpa
      use b2mod_switches
      use b2mod_subsys
#ifdef B25_EIRENE
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use eirmod_parmmod
     .   , only : natm, nmol, nion, nstra, nlmpgs, nspz
      use eirmod_wneutrals
      use eirmod_eirbra
      use eirmod_clogau
      use eirmod_cinit
#endif
      use b2us_map
      use b2us_geo
      implicit none

*  -- INPUT ARGUMENTS --
      integer, intent(in) :: ns, ismain
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(switches), intent(in) :: switch

*  -- OUTPUT ARGUMENTS --
      INTEGER NFEND

*  -- LOCAL SCALARS --
      integer is, nplsi
#ifdef B25_EIRENE
      integer nCvp, iCv, iFc, ix, iy, istra, jvft44_in
      integer nstsi_in, nstrai_in
      integer, allocatable :: idum(:)
      integer i, j, k, ke, kw, ncl
      integer geometryType
      real (kind=R8) :: recyctmp
      character charns*3
      character*18 hlp_frm
      character*72 zeile
      character*256 filename
      logical ex, target_east, target_west, extra_time
#endif
      integer, save :: use_recyceir, iout
      data use_recyceir /1/, iout/0/

*  -- PROCEDURES --
      EXTERNAL GFSUB2
      external ipgeti, xertst
#ifdef B25_EIRENE
      integer lines_for_format
      logical streql
      external streql, lines_for_format

      interface
        subroutine rdneutrs(kard,dummy,ldmf,extra_time)
        use b2mod_types
        use eirmod_parmmod
        integer, intent(in) :: kard, ldmf
        real(kind=R8) :: dummy(nlimps,*)
        logical, intent(in), optional :: extra_time
        end subroutine rdneutrs
      end interface
#endif

*=======================================================================
*//COMPUTATION//
C
      call subini ('ngread')
      call ipgeti ('eirene_use_recyceir', use_recyceir)
      call ipgeti ('b2stbr_iout', iout)
      NFEND = 0
      b2_fnmti = 0

      nplsi = 0
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          nplsi = nplsi + 1
        elseif (switch%spatial_hybrid.ne.0) then
          if (nint(zn(is)).eq.1) nplsi = nplsi + 1
        endif
      end do
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'faulty argument ismain')

#ifdef B25_EIRENE
C
* ----------------------------------------------------------------------
*  READ EIRENE SPECIES DENSITIES AND TEMPERATURES
* ----------------------------------------------------------------------
*
      rewind (44)
C
      READ(44,'(A72)',END=99) ZEILE
      IF (streql(ZEILE(7:8),'  ')) THEN
        READ(ZEILE,'(I6,2X,I8)',ERR=99) nCvp,jvft44_in
      ELSE IF (streql(ZEILE(5:6),'  ')) THEN
        READ(ZEILE,'(I4,2X,I8)',ERR=99) nCvp,jvft44_in
      ELSE
        GOTO 99
      END IF
      READ(44,'(I4,2X,I4,2X,I4)',END=99,ERR=99) NATMI,NMOLI,NIONI
      write(*,'(3a,t27,i6,3i4)')
     . ' Dimensions in ',fort_lc,'44 = ',nCvp,natmi,nmoli,nioni
      call xertst(nCvp .le. mpg%nCv, 'NGREAD: nCvp > nCv')
      call xertst(natmi .le. natm, 'NGREAD: NATMI > NATM')
      call xertst(nmoli .le. nmol, 'NGREAD: NMOLI > NMOL')
      call xertst(nioni .le. nion, 'NGREAD: NIONI > NION')
      call alloc_b2mod_diag(mpg%nCv,mpg%nFc,ns,mpg%nnreg(0),
     .   nstrat,natmi,switch)
      filename='input.dat'
      call find_file(filename,ex)
      filename=fort_lc//'33'
      if (.not.ex) call find_file(filename,nlwrmsh)
c
* labels
* neutral atoms
C
      do 11 IS = 1, NATMI
        READ(44,'(A)',END=99,ERR=99) TEXTAN(IS-1)
        charns=textan(is-1)(1:3)
        call strip_spaces(charns)
        atoms_list(is)=charns(1:2)
   11 CONTINUE
      write(hlp_frm,'(a7,i3,a8)') '(1x,3a,',natmi,'(a2,1x))'
      write(*,hlp_frm) 'atoms_list (',fort_lc,'44) : ',
     . (trim(atoms_list(is)),is=1,natmi)
C
* neutral molecules
C
      do 12 IS = 1, NMOLI
        READ(44,'(A)',END=99,ERR=99) TEXTMN(IS-1)
   12 CONTINUE
C
* test ions
C
      do 13 IS = 1, NIONI
        READ(44,'(A)',END=99,ERR=99) TEXTIN(IS-1)
   13 CONTINUE
C
* data
* neutral atoms
C
      do 14 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,dab2(:,is,1))
        write(charns,'(i3.3)') is
        if (iout.ne.0)
     >   call my_out_us(70,mpg%nCv,0,dab2(1,is,1),'dab2_'//charns)
   14 CONTINUE

      do 24 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,tab2(:,is,1))
   24 CONTINUE

* neutral molecules
C
      do 15 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,dmb2(:,is,1))
        write(charns,'(i3.3)') is
        if (iout.ne.0)
     >   call my_out_us(70,mpg%nCv,0,dmb2(1,is,1),'dmb2_'//charns)
   15 CONTINUE

      do 25 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,tmb2(:,is,1))
   25 CONTINUE

* test ions
C
      do 16 IS = 1, NIONI
        call gfsub2(44,mpg%nCv,nCvp,dib2(:,is,1))
   16 CONTINUE

      do 26 IS = 1, NIONI
        call gfsub2(44,mpg%nCv,nCvp,tib2(:,is,1))
   26 CONTINUE

* neutral atom radial particle flux
C
      do 114 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,rfluxa(:,is,1))
  114 CONTINUE

* neutral molecules radial particle flux
C
      do 115 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,rfluxm(:,is,1))
  115 CONTINUE

* neutral atom poloidal particle flux
C
      do 214 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,pfluxa(:,is,1))
  214 CONTINUE

* neutral molecules poloidal particle flux
C
      do 215 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,pfluxm(:,is,1))
  215 CONTINUE

* neutral atom radial energy flux
C
      do 314 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,refluxa(:,is,1))
  314 CONTINUE

* neutral molecules radial energy flux
C
      do 415 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,refluxm(:,is,1))
  415 CONTINUE

* neutral atom poloidal energy flux
C
      do 514 IS = 1, NATMI
        call gfsub2(44,mpg%nCv,nCvp,pefluxa(:,is,1))
  514 CONTINUE

* neutral molecules poloidal energy flux
C
      do 515 IS = 1, NMOLI
        call gfsub2(44,mpg%nCv,nCvp,pefluxm(:,is,1))
  515 CONTINUE

      call gfsub2(44,mpg%nCv,nCvp,emiss(:,1,1))
      call gfsub2(44,mpg%nCv,nCvp,emissmol(:,1,1))

cank 960511
      do is=1,nmoli
        do iCv=1,mpg%nCv
          srcml(iCv,is,1)=0.0_R8
          edissml(iCv,is,0)=0.0_R8
        end do
      end do
      do istra=0,nstra
        do ix=1,nlimps
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
      if(jvft44_in .ge. 960511) then
        do is=1,nmoli
          call gfsub2(44,mpg%nCv,nCvp,srcml(:,is,1))
        end do
        do is=1,nmoli
          call gfsub2(44,mpg%nCv,nCvp,edissml(:,is,0))
        end do
      end if
      if(jvft44_in .ge. 960513) then
        read(44,'(3i6)') nlimi,nstsi_in,nstrai_in
        allocate(idum(0:nstsi_in+1))
        if (nsts.eq.0) then
          nstsi = nstsi_in
        else
          nstsi = nsts
        end if
        if (nstrai.eq.0) nstrai = nstrai_in
        print '(1x,a,5i6)','nlimi,nstsi_in,nstrai_in,nlimps,nlmpgs:',
     .                      nlimi,nstsi_in,nstrai_in,nlimps,nlmpgs
        call xertst(nlimi.le.nlim, 'NGREAD:  NLIMI > NLIM')
        if (nstsi.ne.nstsi_in) then ! fort.44 contains a time horizon surface
                                    ! and a time census stratum that
                                    ! are no longer wanted
          call xertst(nstsi_in .le.nsts+1, 'NGREAD: NSTSI > NSTS+1')
          call xertst(nstrai_in.le.nstra+1,'NGREAD: NSTRAI > NSTRA+1')
        else
          call xertst(nstsi .le.nsts ,'NGREAD: NSTSI > NSTS')
          call xertst(nstrai.le.nstra,'NGREAD: NSTRAI > NSTRA')
        end if
        extra_time = nstsi_in.eq.nstsi+1
        call rdneutrs(44,wldnek,1,extra_time)
        call rdneutrs(44,wldnep,1,extra_time)
        call rdneutrs(44,wldna,natmi,extra_time)
        call rdneutrs(44,ewlda,natmi,extra_time)
        call rdneutrs(44,wldnm,nmoli,extra_time)
        call rdneutrs(44,ewldm,nmoli,extra_time)
c*** read the wall geometry
        call read_title(44)
        read (44,'(8f10.4)') (p1(1,ix), p1(2,ix),
     ,                        p2(1,ix), p2(2,ix), ix=1,nlimi)
      end if
cank
cank-960727
      if(jvft44_in .ge. 960623) then
        print *,'natmi,nmoli,nioni,nplsi = ',
     .           natmi,nmoli,nioni,nplsi
        call rdneutrs(44,wldra,natmi,extra_time)
        call rdneutrs(44,wldrm,nmoli,extra_time)
      end if
      if(jvft44_in .ge. 960727) then
        if(nstrai_in .gt. 1) then
          do istra=1,nstrai_in
            call rdneutrs(44,wldnek(1,istra),1,extra_time)
            call rdneutrs(44,wldnep(1,istra),1,extra_time)
            call rdneutrs(44,wldna(:,:,istra),natmi,extra_time)
            call rdneutrs(44,ewlda(:,:,istra),natmi,extra_time)
            call rdneutrs(44,wldnm(:,:,istra),nmoli,extra_time)
            call rdneutrs(44,ewldm(:,:,istra),nmoli,extra_time)
            call rdneutrs(44,wldra(:,:,istra),natmi,extra_time)
            call rdneutrs(44,wldrm(:,:,istra),nmoli,extra_time)
          end do
        else
          do ix=1,nlimps
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
        do ix=1,nlimps
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
          do istra=1,nstrai_in
            if(use_recyceir.gt.0) then
              recyctmp=recyceir(istra)
            else
              recyctmp=recyc(ismain,istra)
            endif
            wldnek(ix,0)=wldnek(ix,0)+
     +       wldnek(ix,istra)*recyctmp
            wldnep(ix,0)=wldnep(ix,0)+
     +       wldnep(ix,istra)*recyctmp
            do is=1,natmi
              iy=eb2atcr(is)
              if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
              wldna(ix,is,0)=wldna(ix,is,0)+
     +                             wldna(ix,is,istra)*recyctmp
              wldra(ix,is,0)=wldra(ix,is,0)+
     +                             wldra(ix,is,istra)*recyctmp
              ewlda(ix,is,0)=ewlda(ix,is,0)+
     +                             ewlda(ix,is,istra)*recyctmp
            end do
            do is=1,nmoli
              if (lmolscl(is).gt.0 .and. lmolscl(is).le.natmi) then
                iy=eb2atcr(lmolscl(is))
                if(use_recyceir.eq.0) recyctmp=recyc(iy,istra)
                wldnm(ix,is,0)=wldnm(ix,is,0)+
     +           wldnm(ix,is,istra)*recyctmp
                wldrm(ix,is,0)=wldrm(ix,is,0)+
     +           wldrm(ix,is,istra)*recyctmp
                ewldm(ix,is,0)=ewldm(ix,is,0)+
     +           ewldm(ix,is,istra)*recyctmp
              else if (use_recyceir.gt.0) then !here we are rescaling by stratum, not by species
                wldnm(ix,is,0)=wldnm(ix,is,0)+
     +           wldnm(ix,is,istra)*recyctmp
                wldrm(ix,is,0)=wldrm(ix,is,0)+
     +           wldrm(ix,is,istra)*recyctmp
                ewldm(ix,is,0)=ewldm(ix,is,0)+
     +           ewldm(ix,is,istra)*recyctmp
              else
                if (ix.eq.1 .and. istra.eq.1) then
                  write(*,*) 'Faulty or missing LMOLSCL assignment:'
                  write(*,*) 'LMOLSCL(',is,')=',lmolscl(is)
                  write(*,*) 'Should be valued between 1 and NATM'
                  write(*,*) 'Modify in b2.neutrals.parameters'
                  write(*,*) 'See manual or b2cdcn.f file for usage'
                  write(*,*) 'Skipping scaling of molecular arrays'
                endif
                wldnm(ix,is,0)=wldnm(ix,is,0)+wldnm(ix,is,istra)
                wldrm(ix,is,0)=wldrm(ix,is,0)+wldrm(ix,is,istra)
                ewldm(ix,is,0)=ewldm(ix,is,0)+ewldm(ix,is,istra)
              endif
            end do
          end do
        end do
        if(nspecies .ne. natmi) then
          print *,'Warning: natmi, nspecies= ',natmi,nspecies
          print *,'Source rescaling for additional atomic species ',
     .            'may be incomplete'
        end if
      end if
      if(jvft44_in .ge. 961228) then
        call rdneutrs(44,wldpp,nplsi,extra_time)
        call rdneutrs(44,wldpa,natmi,extra_time)
        call rdneutrs(44,wldpm,nmoli,extra_time)
        call rdneutrs(44,wldpeb,1,extra_time)
        call rdneutrs(44,wldspt,1,extra_time)
        if (jvft44_in .ge. 20170328) then
          call rdneutrs(44,wldspta,natmi,extra_time)
          call rdneutrs(44,wldsptm,nmoli,extra_time)
          do j = 1, natmi
            eirspta(j) = sum(wldspta(:,j,0))
          enddo
          do j = 1, nmoli
            eirsptm(j) = sum(wldsptm(:,j,0))
          enddo
        end if
        if(nstrai_in .gt. 1) then
          do istra=1,nstrai_in
            call rdneutrs(44,wldpp(:,:,istra),nplsi,extra_time)
            call rdneutrs(44,wldpa(:,:,istra),natmi,extra_time)
            call rdneutrs(44,wldpm(:,:,istra),nmoli,extra_time)
            call rdneutrs(44,wldpeb(1,istra),1,extra_time)
            call rdneutrs(44,wldspt(1,istra),1,extra_time)
            if (jvft44_in .ge. 20170328) then
              call rdneutrs(44,wldspta(:,:,istra),natmi,extra_time)
              call rdneutrs(44,wldsptm(:,:,istra),nmoli,extra_time)
            end if
          end do
        end if
      end if
cank
      if(jvft44_in .ge. 20000727) then
        call read_title(44)
        read (44,'(18i4)') (isrftype(i),i=1,nlimi)
        read (44,'(18i4)') (isrftype(nlim+i),i=1,nstsi)
        if (nstsi_in.eq.nstsi+1 .and. mod(nstsi,18).eq.0) read(44,1028)
      end if
      if (NLWRMSH) then
        if(jvft44_in .ge. 20051115) then
          call read_title(44)
          read (44,'(6e13.5)',end=999) (wlarea(i),i=1,nlimi)
          read (44,'(6e13.5)') (wlarea(nlim+i),i=1,nstsi)
          if (nstsi_in.eq.nstsi+1 .and. mod(nstsi,6).eq.0) read(44,1028)
        end if
        if(jvft44_in .ge. 20060206) then
          k=0
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,natmi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*natmi).gt.
     .        lines_for_format(6,nstsi*natmi)) then
            do i = lines_for_format(6,nstsi*natmi)+1,
     .             lines_for_format(6,nstsi_in*natmi)
              read(44,1028)
            end do
          end if
          k=k+natmi
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nmoli),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nmoli).gt.
     .        lines_for_format(6,nstsi*nmoli)) then
            do i = lines_for_format(6,nstsi*nmoli)+1,
     .             lines_for_format(6,nstsi_in*nmoli)
              read(44,1028)
            end do
          end if
          k=k+nmoli
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nioni),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nioni).gt.
     .        lines_for_format(6,nstsi*nioni)) then
            do i = lines_for_format(6,nstsi*nioni)+1,
     .             lines_for_format(6,nstsi_in*nioni)
              read(44,1028)
            end do
          end if
          k=k+nioni
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,i),j=1,nplsi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlabsrp(j+k,nlim+i),j=1,nplsi),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nplsi).gt.
     .        lines_for_format(6,nstsi*nplsi)) then
            do i = lines_for_format(6,nstsi*nplsi)+1,
     .             lines_for_format(6,nstsi_in*nplsi)
              read(44,1028)
            end do
          end if
        end if

        if(jvft44_in .ge. 20170328) then
          k=0
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,natmi)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,natmi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,natmi),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*natmi).gt.
     .        lines_for_format(6,nstsi*natmi)) then
            do i = lines_for_format(6,nstsi*natmi)+1,
     .             lines_for_format(6,nstsi_in*natmi)
              read(44,1028)
            end do
          end if
          k=k+natmi
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nmoli)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nmoli),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nmoli),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nmoli).gt.
     .        lines_for_format(6,nstsi*nmoli)) then
            do i = lines_for_format(6,nstsi*nmoli)+1,
     .             lines_for_format(6,nstsi_in*nmoli)
              read(44,1028)
            end do
          end if
          k=k+nmoli
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nioni)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nioni),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nioni),
     ,                                                       i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nioni).gt.
     .        lines_for_format(6,nstsi*nioni)) then
            do i = lines_for_format(6,nstsi*nioni)+1,
     .             lines_for_format(6,nstsi_in*nioni)
              read(44,1028)
            end do
          end if
          k=k+nioni
          call read_title(44)
          read (44,'(6a13)') (eirtxt(j+k),j=1,nplsi)
          read (44,'(6e13.5)') ((wlpump(j+k,i),j=1,nplsi),i=1,nlimi)
          read (44,'(6e13.5)') ((wlpump(j+k,nlim+i),j=1,nplsi),
     ,                                                      i=1,nstsi)
          if (lines_for_format(6,nstsi_in*nplsi).gt.
     .        lines_for_format(6,nstsi*nplsi)) then
            do i = lines_for_format(6,nstsi*nplsi)+1,
     .             lines_for_format(6,nstsi_in*nplsi)
              read(44,1028)
            end do
          end if

          do j=1,k+nplsi
            eirpump(j) = sum(wlpump(j,:))
          end do
        end if

        if(jvft44_in .ge. 20071209) then
          do is = 1, natmi
            call gfsub2(44,mpg%nCv,nCvp,eneutrad(:,is,0))
          end do
        end if

        if(jvft44_in .ge. 20180323) then
          do is = 1, nmoli
            call gfsub2(44,mpg%nCv,nCvp,emolrad(:,is,0))
          end do
          do is = 1, nioni
            call gfsub2(44,mpg%nCv,nCvp,eionrad(:,is,0))
          end do
        end if

        if(jvft44_in .ge. 20080706) then
          call read_title(44)
          READ (44,'(12i6)') (idum(i),i=1,nstsi_in+1)
          eirdiag_nds_ind(1:nstsi) = idum(1:nstsi)
          eirdiag_nds_ind(nstsi+1) = idum(nstsi_in+1)
          NCL=eirdiag_nds_ind(nstsi+1)-1
          IF (NCL.GT.0) THEN
            IF(.NOT.ALLOCATED(fnmti_b2)) ALLOCATE(fnmti_b2(NCL))
            call read_title(44)
            READ (44,'(12i6)') (fnmti_b2(i),i=1,NCL)
            IF(.NOT.ALLOCATED(wldna_res))
     >         ALLOCATE(wldna_res(natmi,NCL), wldnm_res(nmoli,NCL),
     .                  ewlda_res(natmi,NCL), ewldm_res(nmoli,NCL),
     .                  ewldt_res(NCL),
     .                  ewldea_res(natmi,NCL), ewldem_res(nmoli,NCL),
     .                  ewldrp_res(NCL), ewldmr_res(nmoli,NCL),
     .                  wldspt_res(NCL), sarea_res(NCL),
     .                  wldspta_res(natmi,NCL), wldsptm_res(nmoli,NCL),
     .                  wlpump_res(nspz,NCL) )

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
            if (jvft44_in.ge.20170228) then
              call read_title(44)
              READ (44,'(6e15.7)') (wldspt_res(i),i=1,NCL)
            end if
            if (jvft44_in.ge.20170328) then
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
            if (jvft44_in.ge.20201006) then
              call read_title(44)
              read(44,'(6e15.7)') (ewldt_res(i),i=1,NCL)
            end if
          ELSE
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            call read_title(44)
            READ (44,1028)
            if (jvft44_in.ge.20170228) then
              call read_title(44)
              READ (44,1028)
            end if
            if (jvft44_in.ge.20170328) then
              call read_title(44)
              READ (44,1028)
              call read_title(44)
              READ (44,1028)
              call read_title(44)
              READ (44,1028)
              call read_title(44)
              READ (44,1028)
              call read_title(44)
              READ (44,1028)
              call read_title(44)
              READ (44,1028)
            end if
            if (jvft44_in.ge.20201006) then
              call read_title(44)
              READ (44,1028)
            end if
          END IF

          geometryType = geometryId( mpg, geo )
          if (geometryType.eq.GEOMETRY_LINEAR .or.
     &        geometryType.eq.GEOMETRY_CYLINDER) then
            target_east = .false.
            target_west = .false.
            do i = 1, nbc
              if(boundaryAssignments(
     &  mpg%fcReg(mpg%bcFc(mpg%bcFcP(i,1))),geometryType).eq.EAST) then
                target_east = bcene(i).eq.3 .or. bcene(i).eq.12 .or.
     &                        bcene(i).eq.15.or. target_east
              end if
              if(boundaryAssignments(
     &  mpg%fcReg(mpg%bcFc(mpg%bcFcP(i,1))),geometryType).eq.WEST) then
                target_west = bcene(i).eq.3 .or. bcene(i).eq.12 .or.
     &                        bcene(i).eq.15.or. target_west
              end if
            end do
          else if (geometryType.eq.GEOMETRY_ANNULUS) then
            target_east = .false.
            target_west = .false.
          else
            target_east = .true.
            target_west = .true.
          end if
          kw = 0
          ke = 0
          do i = 1, ntns
            if (ltns(i).eq.0) cycle
            iFc = fnmti_b2(eirdiag_nds_ind(ltns(i))+1)
            if
     &       (boundaryAssignments(mpg%fcReg(iFc),geometryType).eq.WEST)
     &       then
              kw = kw+1
              ktns(i) = 2*kw-1
            else if
     &       (boundaryAssignments(mpg%fcReg(iFc),geometryType).eq.EAST)
     &       then
              ke = ke+1
              if (mpg%nncut.eq.0) then
                if (target_west) then
                  ktns(i) = 2
                else
                  ktns(i) = 1
                end if
              else
                ktns(i) = (3-ke)*mpg%nncut
              end if
            end if
          end do
#ifdef DBG
          if(ntns.gt.0) then !{
            write(hlp_frm,'(a,i3,a)') '(1x,a,3x,',max(ntns,4),'i5)'
            write(*,hlp_frm) 'ltns = ',(ltns(i),i=1,ntns)
            write(*,hlp_frm) 'Calculated ktns :',(ktns(i),i=1,ntns)
          end if !}
#endif
        end if

        if(jvft44_in .ge. 20081111) then
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENA_INT(I,J),I=1,natmi),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*natmi).gt.
     .        lines_for_format(6,(nstrai+1)*natmi)) then
            do i = lines_for_format(6,(nstrai+1)*natmi)+1,
     .             lines_for_format(6,(nstrai_in+1)*natmi)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENM_INT(I,J),I=1,nmoli),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nmoli).gt.
     .        lines_for_format(6,(nstrai+1)*nmoli)) then
            do i = lines_for_format(6,(nstrai+1)*nmoli)+1,
     .             lines_for_format(6,(nstrai_in+1)*nmoli)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENI_INT(I,J),I=1,nioni),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nioni).gt.
     .        lines_for_format(6,(nstrai+1)*nioni)) then
            do i = lines_for_format(6,(nstrai+1)*nioni)+1,
     .             lines_for_format(6,(nstrai_in+1)*nioni)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENA_INT_B2(I,J),I=1,natmi),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*natmi).gt.
     .        lines_for_format(6,(nstrai+1)*natmi)) then
            do i = lines_for_format(6,(nstrai+1)*natmi)+1,
     .             lines_for_format(6,(nstrai_in+1)*natmi)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENM_INT_B2(I,J),I=1,nmoli),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nmoli).gt.
     .        lines_for_format(6,(nstrai+1)*nmoli)) then
            do i = lines_for_format(6,(nstrai+1)*nmoli)+1,
     .             lines_for_format(6,(nstrai_in+1)*nmoli)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((PDENI_INT_B2(I,J),I=1,nioni),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nioni).gt.
     .        lines_for_format(6,(nstrai+1)*nioni)) then
            do i = lines_for_format(6,(nstrai+1)*nioni)+1,
     .             lines_for_format(6,(nstrai_in+1)*nioni)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENA_INT(I,J),I=1,natmi),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*natmi).gt.
     .        lines_for_format(6,(nstrai+1)*natmi)) then
            do i = lines_for_format(6,(nstrai+1)*natmi)+1,
     .             lines_for_format(6,(nstrai_in+1)*natmi)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENM_INT(I,J),I=1,nmoli),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nmoli).gt.
     .        lines_for_format(6,(nstrai+1)*nmoli)) then
            do i = lines_for_format(6,(nstrai+1)*nmoli)+1,
     .             lines_for_format(6,(nstrai_in+1)*nmoli)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENI_INT(I,J),I=1,nioni),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nioni).gt.
     .        lines_for_format(6,(nstrai+1)*nioni)) then
            do i = lines_for_format(6,(nstrai+1)*nioni)+1,
     .             lines_for_format(6,(nstrai_in+1)*nioni)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENA_INT_B2(I,J),I=1,natmi),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*natmi).gt.
     .        lines_for_format(6,(nstrai+1)*natmi)) then
            do i = lines_for_format(6,(nstrai+1)*natmi)+1,
     .             lines_for_format(6,(nstrai_in+1)*natmi)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENM_INT_B2(I,J),I=1,nmoli),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nmoli).gt.
     .        lines_for_format(6,(nstrai+1)*nmoli)) then
            do i = lines_for_format(6,(nstrai+1)*nmoli)+1,
     .             lines_for_format(6,(nstrai_in+1)*nmoli)
              read(44,1028)
            end do
          end if
          call read_title(44)
          READ (44,'(6e15.7)') ((EDENI_INT_B2(I,J),I=1,nioni),
     .                                                    J=0,NSTRAI)
          if (lines_for_format(6,(nstrai_in+1)*nioni).gt.
     .        lines_for_format(6,(nstrai+1)*nioni)) then
            do i = lines_for_format(6,(nstrai+1)*nioni)+1,
     .             lines_for_format(6,(nstrai_in+1)*nioni)
              read(44,1028)
            end do
          end if
        end if
  999   continue
      else  ! NLWRMSH
        if(jvft44_in .ge. 20071209) then
          do is = 1, natmi
            call gfsub2(44,mpg%nCv,nCvp,eneutrad(:,is,0))
          end do
        end if
        if(jvft44_in .ge. 20180323) then
          do is = 1, nmoli
            call gfsub2(44,mpg%nCv,nCvp,emolrad(:,is,0))
          end do
          do is = 1, nioni
            call gfsub2(44,mpg%nCv,nCvp,eionrad(:,is,0))
          end do
        end if
      end if

cc<<<
* -- NO ERROR --
      call extend(dab2,mpg,natm)
      call extend(tab2,mpg,natm)
      call extend(dmb2,mpg,nmol)
      call extend(tmb2,mpg,nmol)
      call extend(dib2,mpg,nion)
      call extend(tib2,mpg,nion)
      call extend(rfluxa,mpg,natm)
      call extend(rfluxm,mpg,nmol)
      call extend(pfluxa,mpg,natm)
      call extend(pfluxm,mpg,nmol)
      call extend(refluxa,mpg,natm)
      call extend(refluxm,mpg,nmol)
      call extend(pefluxa,mpg,natm)
      call extend(pefluxm,mpg,nmol)
      call extend(emiss,mpg,1)
      call extend(emissmol,mpg,1)
cxpb We do not extend the volumetric losses quantities
cank 960511
c      call extend(srcml,mpg,nmol)
c      call extend(edissml,mpg,nmol)
cank
cxpb 20130210
c      call extend(eneutrad,mpg,natm)
c      call extend(emolrad,mpg,nmol)
c      call extend(eionrad,mpg,nion)
cxpb
cxpb 20180323
cxpb
      do i = 1, NCL
        if (fnmti_b2(i).le.0) cycle
        b2_fnmti(fnmti_b2(i)) = i
      end do
 1028 format()
      NFEND = 1
   99 CONTINUE
#endif
      close (44)
*  ---------------------------------------------------------------------
C
      call subend ()
      RETURN

*=======================================================================
*//END NGREAD//

      END SUBROUTINE NGREAD

#ifdef B25_EIRENE
      integer function lines_for_format( nfields_per_line, ndata )
      integer, intent(in) :: nfields_per_line, ndata

      if (ndata.eq.0) then
        lines_for_format = 0
      else if (nfields_per_line.eq.0) then
        call xerrab ('Faulty nfields_per_line argument!')
      else
        lines_for_format = ndata/nfields_per_line
        if (mod(ndata,nfields_per_line).ne.0)
     .   lines_for_format = lines_for_format + 1
      end if

      return
      end function lines_for_format
#endif
