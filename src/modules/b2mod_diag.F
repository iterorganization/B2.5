      module b2mod_diag
      use b2mod_types
      use b2mod_elements
      use b2mod_residuals
      use b2mod_eirene_globals
      use b2mod_neutr_src_scaling
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
     &                    , only : nbcd
      implicit none
c
c  version : 27.04.2005
c
      real (kind=R8) ::
     1 totaddee, totaddei, totaddni, rad_neutr
#ifndef SOLPS4_3
      real (kind=R8) :: b2brem, b2rad
#endif
      real (kind=R8), allocatable :: rad_imp(:)
      REAL (kind=R8), allocatable :: RESHLP(:,:)
      real (kind=R8), allocatable :: rdneureg(:),
     2  b2radirreg(:), b2radlrreg(:),
     3  recycereg(:), recycireg(:),
     4  sseereg(:), sseireg(:), ssnireg(:,:),
     5  ssnastra(:,:),   ! IYS 29.09.2017
     6  sna_eir(:),      ! IYS 12.01.2018
     7  sna_stbc_neut(:,:), sna_stbr_neut(:,:),       ! IYS 04.06.2018
     8  sna_stbr_sput(:,:)
      REAL (kind=R8), allocatable ::
     .      FEEXGT(:,:), FEEYGT(:,:),
     .      FEIXGT(:,:), FEIYGT(:,:)
      real (kind=R8), allocatable :: RADARR(:,:)
      real (kind=R8), allocatable ::
     ,     src_ee(:,:)   , src_ei(:,:),
     ,     src_es(:,:)   , src_ms(:,:),
     ,     src_ni(:,:,:) , src_mo(:,:,:),
     ,     src_ne(:,:)
      real (kind=R8), allocatable ::
     ,     cor_ee(:,:)   , cor_ei(:,:),
     ,     cor_ni(:,:,:) , cor_mo(:,:,:),
     ,     cor_ne(:,:)
      real (kind=R8), allocatable ::
     ,     rescoreg(:,:) , resmoreg(:,:),
     ,     reshereg(:), reshireg(:)

c----------------------------------------------------------------------
c*** sclstra: scaling factors for neutral fluxes per species and strata
c*** sclalph: scaling factors for particle sources
c*** sclsrcn: primary neutral sources, I, from Eirene
c*** sclbase: ratio of plasma to neutral sources, S/I, from Eirene
c*** sclrtio: (S/I from B2) / (S/I from Eirene)

c----------------------------------------------------------------------
      integer lstime, ltime, lmovie, ldiag, ldiag2, ldiag3, ldiag4,
     1 nldetail, nstages
      integer, allocatable :: IZ(:,:), nfluids(:)
c----------------------------------------------------------------------
c*** Specification of the output files and units
c***
c*** ntrc    : unit numbers for tracing files
c*** lv_trc  : their versions
c*** filedata: if true, produce the file
c*** fn_ntrc :  file names for the tracing data
c***
c***  ntrc_srt  :  tracing/intshrt.trc    :  integral data (short version)
c***  ntrc_tst  :  tracing/test.trc       :  test data (can be absent)
c***  ntrc_res  :  tracing/residuals.trc  :  residuals by regions
c***  ntrc_src  :  tracing/sources.trc    :  integral sources
c***  ntrc_bln  :  tracing/blnn.trc       :  particle balance - old style by A.Kukushkin
c***  ntrc_ble  :  tracing/blne.trc       :  energy balance
c***  ntrc_int  :  tracing/integral.trc   :  additional integrals
c***  ntrc_usr  :  tracing/user.trc       :  user-specific data
c***  ntrc_blm  :  tracing/blnm.trc       :  parallel momentum balance
c***  ntrc_sep  :  tracing/sepdata.trc    :  separatrix data
c***  ntrc_bls  :  tracing/blnn_SPb.trc   :  particle balance - new SPb style compatible with speed-up schemes (Kaveeva et al., Nucl. Fusion 58 126018)
c***  ntrc_uss  :  tracing/user_SPb.trc   :  user-specific data introduced by SPb group - contains alternative to user.trc set of data and is more straightforward in programming
c***
c*** print_first_step : if false, skip the first timestep in the history
c***                    tracing

      integer len_ntrc, ntrc_used
      parameter (len_ntrc=12)
      integer ntrc_srt, ntrc_tst, ntrc_res, ntrc_src, ntrc_bln,
     ,        ntrc_ble, ntrc_int, ntrc_usr, ntrc_blm, ntrc_sep,
     ,        ntrc_bls, ntrc_uss
      parameter (ntrc_srt=1, ntrc_tst=2,  ntrc_res=3, ntrc_src=4,
     ,           ntrc_bln=5, ntrc_ble=6,  ntrc_int=7, ntrc_usr=8,
     ,           ntrc_blm=9, ntrc_sep=10, ntrc_bls=11,ntrc_uss=12 )
      integer ntrc(len_ntrc), lv_trc(len_ntrc)
      logical filedata(len_ntrc), print_first_step
      integer lenfnm, mhdr
      parameter (lenfnm=120, mhdr=24)
      character fn_ntrc(len_ntrc)*(lenfnm)
      character*2, allocatable :: atoms_list(:)
      character*2, allocatable :: species_list(:)
      integer, allocatable :: ions_list(:), fluids_list(:)
      character*120 label
      integer :: ank_tracing = 0
c
      save

      contains

      subroutine alloc_b2mod_diag(nxd,nyd,nsd,nregionv,nstra,natm)
      use b2mod_b2cmpa
      implicit none
      integer nxd,nyd,nsd,nregionv,nstra,natm
      integer i, is, iss
      character*18 hlp_frm
      logical streql
      external streql, ipgeti, xertst, xerrab

      call xertst(1.le.nxd, 'faulty input nxd')
      call xertst(1.le.nyd, 'faulty input nyd')
      call xertst(1.le.nsd, 'faulty input nsd')
      call xertst(1.le.nregionv, 'faulty input nnreg(0)')
      call xertst(0.le.nstra, 'faulty input nstra')
      call xertst(1.le.natm, 'faulty input natm')

      call ipgeti ('ank_tracing', ank_tracing)

      if (allocated(IZ)) return
      allocate(IZ(-1:NXD,-1:NYD))
      allocate(RESHLP(-1:NXD,-1:NYD))
      allocate(rad_imp(1:nspecies))
      allocate(rdneureg(0:nregionv),
     2  b2radirreg(nregionv),b2radlrreg(nregionv),
     3  recycereg(nregionv), recycireg(nregionv),
     4  sseereg(nregionv), sseireg(nregionv),
     5  ssnireg(nspecies,nregionv))
      allocate(ssnastra(natm,nstra))      ! IYS 29.09.2017
      allocate(sna_eir(1:nspecies))       ! IYS 12.01.2018
      allocate(sna_stbr_neut(1:nspecies,1:nbcd),
     &         sna_stbc_neut(1:nspecies,1:nbcd),       ! IYS 04.06.2018
     &         sna_stbr_sput(1:nspecies,1:nbcd))       ! IYS 18.09.2018
      allocate(
     .      FEEXGT(-1:NXD,-1:NYD), FEEYGT(-1:NXD,-1:NYD),
     .      FEIXGT(-1:NXD,-1:NYD), FEIYGT(-1:NXD,-1:NYD))
      allocate(RADARR(0:NXD+1,0:NYD+1))
      allocate(species_list(1:nspecies))
      allocate(atoms_list(1:natm))
      allocate(ions_list(1:nsd))
      allocate(fluids_list(0:nsd-1))
      allocate(nfluids(1:nspecies))
      allocate(
     ,     src_ee(0:nregionv,3), src_ei(0:nregionv,3),
     ,     src_es(0:nregionv,3), src_ms(0:nregionv,3),
     ,     src_ni(0:nregionv,3,nspecies),
     ,     src_mo(0:nregionv,3,nspecies),
     ,     src_ne(0:nregionv,3))
      allocate(
     ,     cor_ee(-1:nxd,-1:nyd), cor_ei(-1:nxd,-1:nyd),
     ,     cor_ni(-1:nxd,-1:nyd,nspecies),
     ,     cor_mo(-1:nxd,-1:nyd,nspecies),
     ,     cor_ne(-1:nxd,-1:nyd))
      allocate(
     ,     rescoreg(0:nregionv,0:nsd-1), resmoreg(0:nregionv,0:nsd-1),
     ,     reshereg(0:nregionv)        , reshireg(0:nregionv))

      call xertst(1.le.nspecies.and.nspecies.le.natm,
     . 'faulty parameter nspecies')

      species_list(1:nspecies)='  '
      atoms_list(1:natm)='  '
      ions_list(1:nsd)=-1
      fluids_list(0:nsd-1)=0
      nfluids=0
      i=0
      iss=0
      nstages=0
      do is = 0, nsd-1
        if(species_list(b2espcr(is)).eq.'  ') then
          species_list(b2espcr(is))=is_codes(is)
          iss=0
        endif
        if (.not.is_neutral(is)) then
          nfluids(b2espcr(is))=nfluids(b2espcr(is))+1
          i=i+1
          iss=iss+1
          ions_list(i)=is
          fluids_list(is)=i
        endif
        nstages=max(nstages,iss)
      enddo
      do is = 1, natm
        if (latmscl(is).ne.0.and.streql(atoms_list(is),' ')) then
          atoms_list(is)=species_list(latmscl(is))
        end if
      end do
      write(*,'(1x,a,3(i3,1x))') 'natm, nsd, nstages ',
     .                            natm, nsd, nstages
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i3)'
      write(*,hlp_frm) 'b2espcr : ', b2espcr(0:nsd-1)
      write(*,hlp_frm) 'b2eatcr : ', b2eatcr(0:nsd-1)
      write(hlp_frm,'(a6,i3,a8)') '(1x,a,',nsd,'(a2,1x))'
      write(*,hlp_frm) 'is_codes : ', (trim(is_codes(is)),is=0,nsd-1)
      write(hlp_frm,'(a6,i3,a8)') '(1x,a,',nspecies,'(a2,1x))'
      write(*,hlp_frm) 'species_list : ',
     .                        (trim(species_list(is)),is=1,nspecies)
      write(hlp_frm,'(a6,i3,a8)') '(1x,a,',natm,'(a2,1x))'
      write(*,hlp_frm) 'atoms_list (B2) : ',
     .                        (trim(atoms_list(is)),is=1,natm)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',i,'i4)'
      write(*,hlp_frm) 'ions_list : ', ions_list(1:i)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nspecies,'i4)'
      write(*,hlp_frm) 'nfluids : ', nfluids(1:nspecies)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i4)'
      write(*,hlp_frm) 'fluids_list : ', fluids_list(0:nsd-1)

c
c*** Tracing unit assignment (see b2file, b2blnc, etc.)
c
      NTRC(ntrc_tst) = 39
      fn_ntrc(ntrc_tst) = 'tracing/test.trc'
      NTRC(ntrc_srt) = 59
      fn_ntrc(ntrc_srt) = 'tracing/intshrt.trc'
      NTRC(ntrc_res) = 71
      fn_ntrc(ntrc_res) = 'tracing/residuals.trc'
      NTRC(ntrc_src) = 72
      fn_ntrc(ntrc_src) = 'tracing/sources.trc'
      NTRC(ntrc_bln) = 73
      fn_ntrc(ntrc_bln) = 'tracing/blnn.trc'
      NTRC(ntrc_bls) = 74                        ! IYS 07.11.2018
      fn_ntrc(ntrc_bls) = 'tracing/blnn_SPb.trc'  ! IYS 07.11.2018
      NTRC(ntrc_ble) = 75
      fn_ntrc(ntrc_ble) = 'tracing/blne.trc'
      NTRC(ntrc_int) = 76
      fn_ntrc(ntrc_int) = 'tracing/integral.trc'
      NTRC(ntrc_usr) = 77
      fn_ntrc(ntrc_usr) = 'tracing/user.trc'
      NTRC(ntrc_uss) = 78                        ! IYS 22.03.2022
      fn_ntrc(ntrc_uss) = 'tracing/user_SPb.trc'  ! IYS 22.03.2022
      NTRC(ntrc_blm) = 79
      fn_ntrc(ntrc_blm) = 'tracing/blnm.trc'
      NTRC(ntrc_sep) = 80
      fn_ntrc(ntrc_sep) = 'tracing/sepdata.trc'

      ntrc_used = 12

      if(ntrc_used.gt.len_ntrc) then
        call xerrab('default: insufficient length of ntrc')
      end if

      do i=1,len_ntrc
        lv_trc(i)=0
        filedata(i)=.true.
      end do
c<<<
      filedata(ntrc_tst)=.true.
c>>>

      rad_imp=0.0_R8
      rdneureg=0.0_R8
      sseereg=0.0_R8
      sseireg=0.0_R8
      ssnireg=0.0_R8
      ssnastra=0.0_R8
      sna_eir=0.0_R8        ! IYS 12.01.2018
      sna_stbr_sput=0.0_R8  ! IYS 28.09.2018
      sna_stbr_neut=0.0_R8  ! IYS 04.06.2018
      sna_stbc_neut=0.0_R8  ! IYS 04.06.2018
      rescoreg=0.0_R8
      resmoreg=0.0_R8
      reshereg=0.0_R8
      reshireg=0.0_R8

      end subroutine alloc_b2mod_diag

      subroutine dealloc_b2mod_diag
      implicit none

      if (.not.allocated(IZ)) return
      deallocate(IZ)
      deallocate(RESHLP)
      deallocate(rad_imp)
      deallocate(rdneureg,
     2  b2radirreg, b2radlrreg,
     3  recycereg, recycireg,
     4  sseereg, sseireg, ssnireg,
     5  ssnastra)
      deallocate(sna_eir)
      deallocate(sna_stbr_neut, sna_stbr_sput)
      deallocate(sna_stbc_neut)  ! IYS 12.01.2018, 04.06.2018
      deallocate(FEEXGT, FEEYGT, FEIXGT, FEIYGT)
      deallocate(RADARR)
      deallocate(species_list)
      deallocate(atoms_list)
      deallocate(ions_list)
      deallocate(fluids_list)
      deallocate(nfluids)
      deallocate(src_ee, src_ei, src_es, src_ms, src_ni, src_mo, src_ne)
      deallocate(cor_ee, cor_ei, cor_ni, cor_mo, cor_ne)
      deallocate(rescoreg, resmoreg, reshereg, reshireg)

      return
      end subroutine dealloc_b2mod_diag

      end module b2mod_diag
