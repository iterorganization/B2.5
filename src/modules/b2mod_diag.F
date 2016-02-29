      module b2mod_diag
      use b2mod_types
      use b2mod_elements
      use b2mod_residuals
      use b2mod_eirene_globals
      use b2mod_neutr_src_scaling
      use b2mod_neutrals_namelist
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
     4  sseereg(:), sseireg(:), ssnireg(:,:)
      REAL (kind=R8), allocatable :: 
     .      FEEXGT(:,:), FEEYGT(:,:),
     .      FEIXGT(:,:), FEIYGT(:,:)
      real (kind=R8), allocatable ::  RADARR(:,:)
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
     1 nldetail, nltrsol, nstages
      integer, allocatable :: IZ(:,:), nfluids(:)
c----------------------------------------------------------------------
c*** Specification of the output files and units
c***
c*** ntrc    : unit numbers for tracing files
c*** lv_trc  : their versions
c*** filedata: if true, produce the file, else dump into b2.log
c*** fn_ntrc :  file names for the tracing data
c***
c***  ntrc_mnp  :  B2MNPRT                :  traditional data
c***  ntrc_tst  :  tracing/test.trc       :  test data (can be absent)
c***  ntrc_res  :  tracing/residuals.trc  :  residuals by regions
c***  ntrc_src  :  tracing/sources.trc    :  integral sources
c***  ntrc_bln  :  tracing/blnn.trc       :  particle balance
c***  ntrc_ble  :  tracing/blne.trc       :  energy balance
c***  ntrc_int  :  tracing/integral.trc   :  additional integrals
c***  ntrc_usr  :  tracing/user.trc       :  user-specific data
c***  ntrc_blm  :  tracing/blnm.trc       :  parallel momentum balance
c***  ntrc_sep  :  tracing/sepdata.trc    :  separatrix data
c***
c*** print_first_step : if false, skip the first timestep in the history
c***                    tracing

      integer len_ntrc, ntrc_used
      parameter (len_ntrc=10)
      integer ntrc_mnp, ntrc_tst, ntrc_res, ntrc_src, ntrc_bln,
     ,        ntrc_ble, ntrc_int, ntrc_usr, ntrc_blm, ntrc_sep
      parameter (ntrc_mnp=1, ntrc_tst=2, ntrc_res=3, ntrc_src=4,
     ,           ntrc_bln=5, ntrc_ble=6, ntrc_int=7, ntrc_usr=8,
     ,           ntrc_blm=9, ntrc_sep=10)
      integer ntrc(len_ntrc), lv_trc(len_ntrc)
      logical filedata(len_ntrc), print_first_step
      character fn_ntrc(len_ntrc)*32
      character*2, allocatable :: atoms_list(:)
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
      allocate(rad_imp(1:natm))
      allocate(rdneureg(nregionv),
     2  b2radirreg(nregionv),b2radlrreg(nregionv),
     3  recycereg(nregionv), recycireg(nregionv),
     4  sseereg(nregionv), sseireg(nregionv), ssnireg(natm,nregionv))
      allocate( 
     .      FEEXGT(-1:NXD,-1:NYD), FEEYGT(-1:NXD,-1:NYD),
     .      FEIXGT(-1:NXD,-1:NYD), FEIYGT(-1:NXD,-1:NYD))
      allocate(RADARR(0:NXD+1,0:NYD+1))
      allocate(atoms_list(1:natm))
      allocate(ions_list(1:nsd))
      allocate(fluids_list(0:nsd-1))
      allocate(nfluids(1:natm))
      allocate(
     ,     src_ee(0:nregionv,3)   , src_ei(0:nregionv,3),
     ,     src_es(0:nregionv,3)   , src_ms(0:nregionv,3),
     ,     src_ni(0:nregionv,3,natm), src_mo(0:nregionv,3,natm),
     ,     src_ne(0:nregionv,3))
      allocate(
     ,     cor_ee(-1:nxd,-1:nyd)     , cor_ei(-1:nxd,-1:nyd),
     ,     cor_ni(-1:nxd,-1:nyd,natm), cor_mo(-1:nxd,-1:nyd,natm),
     ,     cor_ne(-1:nxd,-1:nyd))
      allocate(
     ,     rescoreg(0:nregionv,0:nsd-1), resmoreg(0:nregionv,0:nsd-1),
     ,     reshereg(0:nregionv)        , reshireg(0:nregionv))

      call xertst(1.le.nspecies.and.nspecies.le.natm,
     . 'faulty parameter nspecies')

      atoms_list(1:natm)='  '
      ions_list(1:nsd)=-1
      fluids_list(0:nsd-1)=0
      nfluids=0
      i=0
      iss=0
      nstages=0
      do is = 0, nsd-1
        if(atoms_list(b2espcr(is)).eq.'  ') then
          atoms_list(b2espcr(is))=is_codes(is)
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
      write(*,'(1x,a,3i3)') 'natm, nsd, nstages ', natm, nsd, nstages
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i3)'
      write(*,hlp_frm) 'b2espcr : ', b2espcr(0:nsd-1)
      write(hlp_frm,'(a6,i3,a8)') '(1x,a,',nsd,'(a2,1x))'
      write(*,hlp_frm) 'is_codes : ', is_codes
      write(hlp_frm,'(a6,i3,a8)') '(1x,a,',nspecies,'(a2,1x))'
      write(*,hlp_frm) 'atoms_list : ', atoms_list(1:nspecies)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',i,'i4)'
      write(*,hlp_frm) 'ions_list : ', ions_list(1:i)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nspecies,'i4)'
      write(*,hlp_frm) 'nfluids : ', nfluids(1:nspecies)
      write(hlp_frm,'(a6,i3,a3)') '(1x,a,',nsd,'i4)'
      write(*,hlp_frm) 'fluids_list : ', fluids_list(0:nsd-1)

c
c*** Tracing unit assignment (see b2file, b2blnc, etc.)
c
c<<<
      NTRC(ntrc_tst) = 70
      fn_ntrc(ntrc_tst) = 'tracing/test.trc'
c>>>
      NTRC(ntrc_mnp) = 71
      fn_ntrc(ntrc_mnp)='B2MNPRT'
      NTRC(ntrc_res) = 72
      fn_ntrc(ntrc_res) = 'tracing/residuals.trc'
      NTRC(ntrc_src) = 73
      fn_ntrc(ntrc_src) = 'tracing/sources.trc'
      NTRC(ntrc_bln) = 74
      fn_ntrc(ntrc_bln) = 'tracing/blnn.trc'
      NTRC(ntrc_ble) = 75
      fn_ntrc(ntrc_ble) = 'tracing/blne.trc'
      NTRC(ntrc_int) = 76
      fn_ntrc(ntrc_int) = 'tracing/integral.trc'
      NTRC(ntrc_usr) = 77
      fn_ntrc(ntrc_usr) = 'tracing/user.trc'
      NTRC(ntrc_blm) = 78
      fn_ntrc(ntrc_blm) = 'tracing/blnm.trc'
      NTRC(ntrc_sep) = 79
      fn_ntrc(ntrc_sep) = 'tracing/sepdata.trc'

      ntrc_used = 10

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
     4  sseereg, sseireg, ssnireg)
      deallocate( 
     .      FEEXGT, FEEYGT,
     .      FEIXGT, FEIYGT)
      deallocate(RADARR)
      deallocate(atoms_list)
      deallocate(ions_list)
      deallocate(fluids_list)
      deallocate(nfluids)
      deallocate(
     ,     src_ee, src_ei,
     ,     src_es, src_ms,
     ,     src_ni, src_mo,
     ,     src_ne)
      deallocate(
     ,     cor_ee, cor_ei,
     ,     cor_ni, cor_mo,
     ,     cor_ne)
      deallocate(
     ,     rescoreg, resmoreg,
     ,     reshereg, reshireg)

      return
      end subroutine dealloc_b2mod_diag

      end module b2mod_diag
