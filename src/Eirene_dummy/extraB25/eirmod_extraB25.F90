      module eirmod_extrab25
      use eirmod_precision
      use eirmod_parmmod
      implicit none
      private

      public :: eirene_extrab25_cleanup,eirene_extrab25_wneutrals
      public :: eirene_extrab25_wneuinit, eirene_extrab25_wneufill
      public :: eirene_extrab25_wneusave, eirene_extraB25_wneuclean
      public :: eirene_extrab25_eirpbls !,eirene_extrab25_eirpbls_init !srv 03.04.13
      public :: eirene_extrab25_srfprvsl,eirene_extrab25_srfprvsl_init
      public :: eirene_extrab25_alloc_mods
      public :: eirene_extrab25_iniusr_init
      public :: eirene_extrab25_eirdiag_alloc                          !srv 14.10.13
      public :: eirene_extrab25_eirpbls_alloc                          !srv 03.04.13

      real*8, save, allocatable, dimension(:,:,:,:),public :: &
       dab2,dmb2,dib2,tab2,tmb2,tib2,rfluxa,rfluxm,refluxa,refluxm,&
       pfluxa,pfluxm,pefluxa,pefluxm,emiss,emissmol,srcml,edissml
      real*8, save, allocatable, dimension(:,:),public :: &
       wldnek,wldnep
      real*8, save, allocatable, dimension (:,:,:), public :: &
       wldna,ewlda,wldnm,ewldm,wldra,wldrm,wldpp,wldpa,wldpm
      real*8, save, allocatable, dimension (:,:),public :: &
       wldpeb,wldspt
      real*8, save, allocatable, dimension (:,:,:,:), public :: &
       eneutrad
      integer, save, public :: nnlimi,nnstsi,nnatmi,nnmoli,nnioni
      integer, save, public :: nnplsi,nns, nnstrai
      integer, save, allocatable, public :: isrftype(:)
      logical, save, public :: lhalpha=.false.,lvib=.false.


      ! eirpbls globals
      integer, save, allocatable,public :: lkindp(:),lkindm(:),lkindi(:)
      real*8, save, allocatable,public :: srccrfc(:,:), flxspci(:,:)
      real*8, save, allocatable, public :: srcstrn(:)

      ! srfprvls
      integer, save, public :: nsrfcls
      integer, public, parameter :: msrfclsx=12
      integer, save, public, allocatable :: msrfcls(:),lsrfcls(:,:)

      ! b2.5 neutrals parameters modicifcations
      integer, save, public :: bn_spcsrf
      integer, save, public, allocatable :: bl_spcsrf(:),bi_spcsrf(:)
      integer, save, public, allocatable :: bj_spcsrf(:),bsps_sgrp(:)
      real*8, save, public, allocatable :: bsps_absr(:),bsps_trno(:)
      real*8, save, public, allocatable :: bsps_mtri(:),bsps_tmpr(:)
      real*8, save, public, allocatable :: bsps_trni(:), bsps_spph(:)
      real*8, save, public, allocatable :: bsps_spch(:)
      character*8, save, public, allocatable :: bsps_mtrl(:), bsps_id(:)


      ! remaining bits and pieces from braeir common
      real*8, save, public :: chemical_sputter_yield,fchar_chemical
      integer, save, public :: igass_chemical,itsput_chemical, &
                               issput_chemical

      !pb volumes of b2.5 cells
      real*8, save, public, allocatable :: volcel(:,:)

      !flux_save
      real*8, save, public, allocatable :: flux_save(:)

      contains

      subroutine eirene_extrab25_alloc_mods(nnx,nny)
      implicit none
      integer, intent(in) :: nnx,nny
      end subroutine

      subroutine eirene_extrab25_eirdiag_alloc(ndxp,ndyp, &
                             natm,nmol,nion,nlmpgs,nstra,npls)                  !srv 14.10.13
      implicit none
      integer, intent(in) :: ndxp,ndyp,natm,nmol,nion,nlmpgs,nstra,npls
      if(.not.allocated(dab2)) then                 !srv 14.10.13
      allocate(dab2(0:ndxp,0:ndyp,natm,1))
      allocate(dmb2(0:ndxp,0:ndyp,nmol,1))
      allocate(dib2(0:ndxp,0:ndyp,nion,1))
      allocate(tab2(0:ndxp,0:ndyp,natm,1))
      allocate(tmb2(0:ndxp,0:ndyp,nmol,1))
      allocate(tib2(0:ndxp,0:ndyp,nion,1))
      allocate(rfluxa(0:ndxp,0:ndyp,natm,1))
      allocate(rfluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(refluxa(0:ndxp,0:ndyp,natm,1))
      allocate(refluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(pfluxa(0:ndxp,0:ndyp,natm,1))
      allocate(pfluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(pefluxa(0:ndxp,0:ndyp,natm,1))
      allocate(pefluxm(0:ndxp,0:ndyp,nmol,1))
      allocate(emiss(0:ndxp,0:ndyp,1,1))
      allocate(emissmol(0:ndxp,0:ndyp,1,1))
      allocate(srcml(0:ndxp,0:ndyp,nmol,1))
      allocate(edissml(0:ndxp,0:ndyp,nmol,1))
      allocate(wldnek(nlmpgs,0:nstra+1)) 
      allocate(wldnep(nlmpgs,0:nstra+1))
      allocate(wldna(nlmpgs,natm,0:nstra+1))
      allocate(ewlda(nlmpgs,natm,0:nstra+1))
      allocate(wldnm(nlmpgs,nmol,0:nstra+1))
      allocate(ewldm(nlmpgs,nmol,0:nstra+1))
      allocate(wldra(nlmpgs,natm,0:nstra+1))
      allocate(wldrm(nlmpgs,nmol,0:nstra+1))
      allocate(wldpp(nlmpgs,npls,0:nstra+1))
      allocate(wldpa(nlmpgs,natm,0:nstra+1))
      allocate(wldpm(nlmpgs,nmol,0:nstra+1))
      allocate(wldpeb(nlmpgs,0:nstra+1))
      allocate(wldspt(nlmpgs,0:nstra+1))
      allocate(eneutrad(0:ndxp,0:ndyp,natm,1))
      write(*,*) ' Allocated arrays in eirene_extrab25_eirdiag_alloc ',  &
      ' with ndxp,ndyp,natm,nmol,nion,nlmpgs,nstra,npls',             &
             ndxp,ndyp,natm,nmol,nion,nlmpgs,nstra,npls
      dab2=0
      dmb2=0
      dib2=0
      tab2=0
      tmb2=0
      tib2=0
      rfluxa=0
      rfluxm=0
      refluxa=0
      refluxm=0
      pfluxa=0
      pfluxm=0
      pefluxa=0
      pefluxm=0
      emiss=0
      emissmol=0
      srcml=0
      edissml=0
      wldnek=0
      wldnep=0
      wldna=0
      ewlda=0
      wldnm=0
      ewldm=0
      wldra=0
      wldrm=0
      wldpp=0
      wldpa=0
      wldpm=0
      wldpeb=0
      wldspt=0
      eneutrad=0
	endif                     !srv 14.10.13
      end subroutine eirene_extrab25_eirdiag_alloc

      subroutine eirene_extrab25_wneutrals
      implicit none
      integer :: istra_in
      entry eirene_extraB25_wneuinit
      entry eirene_extraB25_wneufill(istra_in)
      entry eirene_extraB25_wneusave
      return
      end subroutine

      subroutine eirene_extrab25_wneuclean
      implicit none
      return
      end subroutine


      subroutine eirene_extrab25_eirpbls(istr)
      implicit none
      integer :: istr,im,na,nm,ni,np,nst,nnfl
      entry eirene_extrab25_eirpbls_alloc(im,na,nm,ni,np,nst,nnfl)      !srv 03.04.13 14.10.13 {
      select case (im)
      case (1)
!        allocate(lkindm(nm))
!        allocate(lkindp(np))
!        allocate(lkindi(ni))
!        allocate(srccrfc(na,nst))
       if(.not.allocated(srcstrn)) then                                  !srv 14.10.13
        allocate(srcstrn(nst))
!        lkindm=0
!        lkindp=0
!        lkindi=0
!        srccrfc=0.
        srcstrn=0.
       endif                                                             !srv 14.10.13
      case (2)
!        allocate(flxspci(nnfl,nst))
!        flxspci=0.
      end select                                                         !srv 03.04.13 } 
      return
      end subroutine               

      subroutine eirene_extrab25_srfprvsl_init
      implicit none
      end subroutine

      subroutine eirene_extrab25_srfprvsl
      implicit none
      end subroutine

      subroutine eirene_extrab25_iniusr_init(n_spcsrf,l_spcsrf,&
                i_spcsrf,&
                j_spcsrf,sps_sgrp,sps_absr,sps_trno,sps_trni,&
                sps_mtri,sps_tmpr,sps_spph,sps_spch,&
                sps_mtrl,sps_id)
      implicit none
      integer, intent(in) :: n_spcsrf,l_spcsrf(:),i_spcsrf(:),&
            j_spcsrf(:),sps_sgrp(:)
      real*8, intent(in) :: sps_absr(:),sps_trno(:),sps_trni(:),&
            sps_mtri(:), sps_tmpr(:), sps_spph(:), sps_spch(:)
      character*8, intent(in) :: sps_mtrl(:),sps_id(:)
      end subroutine

      subroutine eirene_extrab25_cleanup
      implicit none
      if(allocated(dab2)) then
        deallocate(dab2)
        deallocate(dmb2)
        deallocate(dib2)
        deallocate(tab2)
        deallocate(tmb2)
        deallocate(tib2)
        deallocate(rfluxa)
        deallocate(rfluxm)
        deallocate(refluxa)
        deallocate(refluxm)
        deallocate(pfluxa)
        deallocate(pfluxm)
        deallocate(pefluxa)
        deallocate(pefluxm)
        deallocate(emiss)
        deallocate(emissmol)
        deallocate(srcml)
        deallocate(edissml)
        deallocate(wldnek) 
        deallocate(wldnep)
        deallocate(wldna)
        deallocate(ewlda)
        deallocate(wldnm)
        deallocate(ewldm)
        deallocate(wldra)
        deallocate(wldrm)
        deallocate(wldpp)
        deallocate(wldpa)
        deallocate(wldpm)
        deallocate(wldpeb)
        deallocate(wldspt)
        deallocate(isrftype)
        deallocate(eneutrad)
      endif
      if(allocated(srcstrn)) then                  !srv 03.04.13
        deallocate(srcstrn)
      endif                                        !srv 03.04.13

      end subroutine

      end module
