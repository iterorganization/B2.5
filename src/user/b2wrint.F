      subroutine b2wrint(nx,ny,ns,nregionv)
      use b2mod_types
      use b2mod_geo
      use b2mod_time
      use b2mod_diag
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_external
      use b2mod_constants
      use b2mod_numerics_namelist
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_eirbra
#else
      use b2mod_eirbra
#endif
c
c  version : 15.03.2002 16:33
c
c========================================================================
c*** Puts some integral quantities on the output file
c*** "tracing/integral.trc"
c-----------------------------------------------------------------------
      implicit none
c
#     include <DIMENSIONS.F>
c-----------------------------------------------------------------------
c*** lv    :  version number to appear in the file
c
      integer nx, ny, ns, nregionv
c
      integer lv
c      parameter (lv=970522)
c#ifndef SCL_IMP_RCI
      parameter (lv=980813)
c#else
c      parameter (lv=20010716)
c#endif
c
c-----------------------------------------------------------------------
c*** Additional integrals:
c***
c***      header     var       Legend
c***...........................................................(global)
c***     SSEE      hlp_ssee   neutral-related source of electron energy
c***     SSEI      hlp_ssei   neutral-related source of ion energy
c***     SNI_rcmb  hlp_sumn   recombination source of ions ([A]!!!)
c***     SEE_rcmb  hlp_sume   recombination source of electron energy
c***     SEI_rcmb  hlp_sumi   recombination source of ion energy
c***     Zeff_cor  hlp_zeff   volume-average Zeff at the core boundary
c***     Zeff_crs  hlp_zefs   surface-average Zeff at the core boundary
c***     Zeff_crn  hlp_zefn   calculated Zeff at the core boundary
c***     Brmsreg   hlp_brms   Bremsstrahlung
c***     imp_rad   hlp_rdim   radiation from impurity ions
c***     neu_rad   hlp_rdnu   radiation from neutrals
c***     tot_rad   hlp_rdtt   total radiation
c***
c***..........................................................(regions)
c***     Brmsreg_  brms_tmp   Bremsstrahlung
c***     imp_rad_  radim_tmp  radiation from impurity ions
c***     neu_rad_  radnu_tmp  radiation from neutrals
c***     tot_rad_  radtt_tmp  total radiation
c***     SSEEreg_  ssee_tmp   neutral-related source of electron energy
c***     SSEIreg_  ssei_tmp   neutral-related source of ion energy
c***     th_enre_  eesum      internal energy of electrons \int 3/2 n_e T_e dV
c***     th_enri_  eisum      \int 3/2 n_i T_i dV
c***     N_e_reg_  nesum      total number of electrons \int n_e dV
c***     vol_reg__ volsum     volume of whole B2.5 domain
c***
c***..........................................................(species)
c***     SSNItot_  ssni_tmp   neutral-related source of ions
c***
c***..............................................(species and regions)
c***     SSNIreg_  ssna_tmp   neutral-related source of ions
c***
c***...........................................................(strata)
c***     Stratum_  fluxt_tmp  primary sources of neutrals
c***
c***...............................................(species and strata)
c***     SSNIstr_             neutral-related source of ions
c***     sclfact_  sclf_tmp   scaling factors for neutral fluxes
c***     sclrtio_  rtio_tmp   figure of merit for scaling (should be 0)
c***
c***..........................................................(species)
c***     res_rci_  rrci_tmp   residuals in impurity balance
c-----------------------------------------------------------------------
c
c*** Global values
      integer nglb_tmp, ng_tmp
      parameter (nglb_tmp=12) ! Number of global parameters
      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     hlp_ssee, hlp_ssei, hlp_sumn, hlp_sume, hlp_sumi, hlp_zeff,
     ,     hlp_zefs, hlp_zefn, hlp_brms, hlp_rdim, hlp_rdnu, hlp_rdtt
      equivalence
     ,     (equg_tmp(1),  hlp_ssee), (equg_tmp(2),  hlp_ssei),
     ,     (equg_tmp(3),  hlp_sumn), (equg_tmp(4),  hlp_sume),
     ,     (equg_tmp(5),  hlp_sumi), (equg_tmp(6),  hlp_zeff),
     ,     (equg_tmp(7),  hlp_zefs), (equg_tmp(8),  hlp_zefn),
     ,     (equg_tmp(9),  hlp_brms), (equg_tmp(10), hlp_rdim),
     ,     (equg_tmp(11), hlp_rdnu), (equg_tmp(12), hlp_rdtt)
      character hlp_hg(nglb_tmp)*8
c
c*** Values by regions
      integer nrgn_tmp, nrg_tmp, nr_tmp
      parameter (nrgn_tmp=11) ! Number of quantities computed in each region
      real (kind=R8) ::
     * radim_tmp(NNREGMAX), radnu_tmp(NNREGMAX),
     1 radtt_tmp(NNREGMAX), ssee_tmp(NNREGMAX),
     2 ssei_tmp(NNREGMAX),  brms_tmp(NNREGMAX),
     3 eesum(NNREGMAX),     eisum(NNREGMAX),
     4 kisum(NNREGMAX),
     5 nesum(NNREGMAX),     volsum(NNREGMAX)

      real (kind=R8) ::
     *     equr_tmp(NNREGMAX,nrgn_tmp)
      equivalence
     , (equr_tmp(1,1), brms_tmp),  (equr_tmp(1,2), radim_tmp),
     , (equr_tmp(1,3), radnu_tmp), (equr_tmp(1,4), radtt_tmp),
     , (equr_tmp(1,5), ssee_tmp),  (equr_tmp(1,6), ssei_tmp),
     , (equr_tmp(1,7), eesum),     (equr_tmp(1,8), eisum),
     , (equr_tmp(1,9), kisum),
     , (equr_tmp(1,10), nesum),    (equr_tmp(1,11), volsum)
      character*2 hlp_hr(nrgn_tmp)*8
c
c*** Values by species
      integer nspc_tmp, nsp_tmp, ns_tmp
!      parameter (nspc_tmp=1)
      parameter (nspc_tmp=2) ! IYS 29.09.2017 to add total number of neutrals in whole domain
      real (kind=R8) ::
     *     equs_tmp(DEF_NATM,nspc_tmp)
      real (kind=R8) ::
     *     ssni_tmp(DEF_NATM), Natm_tmp(DEF_NATM)
      equivalence (equs_tmp(1,1), ssni_tmp), (equs_tmp(1,2), Natm_tmp)
      character*2 hlp_hs(nspc_tmp)*8

c*** added 20010716
c#ifdef SCL_IMP_RCI
c      integer nspc1_tmp, ns1_tmp
c      parameter (nspc1_tmp=1)
c      real (kind=R8) ::
c    *     equs1_tmp(nsd,nspc1_tmp)
c      real(kind=R8) ::
c    *      rrci_tmp(nsd)
c      equivalence (equs1_tmp(1,1), rrci_tmp)
c      character*2 hlp_hs1(nspc_tmp)*8
c#endif
c
c*** Values by species and regions
      integer nsar_tmp, na_tmp
!      parameter (nsar_tmp=2)
      parameter (nsar_tmp=3) ! IYS to add one more record to count the total
                             ! number of ions (charged) and total number of
                             ! neutrals
      real (kind=R8) ::
     *     equa_tmp(DEF_NATM,NNREGMAX,nsar_tmp)
      real (kind=R8) ::
     *     nasum(DEF_NATM,NNREGMAX),  ! IYS 29.09.2017 added total number of neutrals
     *     ssna_tmp(DEF_NATM,NNREGMAX),nisum(DEF_NATM,NNREGMAX)
      equivalence 
     -  (equa_tmp(1,1,1), ssna_tmp), (equa_tmp(1,1,2), nisum),
     &  (equa_tmp(1,1,3), nasum) ! IYS 29.09.2017
      character*2 hlp_ha(nsar_tmp)*8
c
c*** Values by strata
      integer nstr_tmp, nf_tmp, nfl_tmp
      parameter (nstr_tmp=1)
      real (kind=R8) ::
     *     equf_tmp(DEF_NSTRA,nstr_tmp)
      real (kind=R8) :: fluxt_tmp(DEF_NSTRA)
      character hlp_hf(nstr_tmp)*8
c
c*** Values by species and strata
      integer nssb_tmp, nb_tmp
c      parameter (nssb_tmp=2)
      parameter (nssb_tmp=4)  ! IYS 29.09.2017
      real (kind=R8) ::
     *     equb_tmp(0:DEF_NATM,DEF_NSTRA,nssb_tmp)
      real (kind=R8) ::
     *     sclf_tmp(0:DEF_NATM,DEF_NSTRA), 
     *     rtio_tmp(0:DEF_NATM,DEF_NSTRA),
     *     Na_stra_tmp(0:DEF_NATM,DEF_NSTRA),
     *     ssna_stra_tmp(0:DEF_NATM,DEF_NSTRA)
      equivalence (equb_tmp(0,1,1), sclf_tmp)
      equivalence (equb_tmp(0,1,2), rtio_tmp)
      equivalence (equb_tmp(0,1,3), Na_stra_tmp)    ! IYS 29.09.2017
      equivalence (equb_tmp(0,1,4), ssna_stra_tmp)  ! IYS 29.09.2017
      character hlp_hb(nssb_tmp)*8
c-----------------------------------------------------------------------
c*** Local variables of general use
      integer lun, ltrc
      parameter (ltrc=ntrc_int)
      real (kind=R8) ::
     *     hlp_s,hlp_v,hlp_zs,hlp_zv,hlp_zn1,hlp_zn2,zeff
      integer i, j, k, ir,is,ix,iy, nv
      logical ex, warn_tmp
      character tab, hlp_c, f*(lenfnm), hlp_frm*256, hlp_frf*400
      integer iatm,imol ! IYS 29.09.2017
c#ifdef SCL_IMP_RCI
c     integer iy
c#endif
c-----------------------------------------------------------------------
c*** Column headers:
      data hlp_hg / ! These are headers for 12 global parameters
     1          'SSEE_tot','SSEI_tot','SNI_rcmb','SEE_rcmb','SEI_rcmb',
     2          'Zeff_cor','Zeff_crs','Zeff_crn', 'Brmsreg ','imp_rad ',
     2          'neu_rad ','tot_rad '/
      data hlp_hr / ! These are headers for 10 quantities computed in each region
     1          'Brmsreg_','imp_rad_','neu_rad_','tot_rad_','SSEEreg_',
     2          'SSEIreg_','th_enre_','th_enri_',
     3          'K_en_ri_','N_e_reg_','vol_reg_'/
      data hlp_hs / 'SSNItot_', 'N_a_tot_' /
c#ifdef SCL_IMP_RCI
c      data hlp_hs1 / 'res_rci_' /
c#endif
      data hlp_ha / 'SSNIreg_', 'N_i_reg_', 'N_a_reg_' /
      data hlp_hf / 'Stratum_' /
      data hlp_hb / 'sclfact_', 'sclrtio_', 'Na_stra_', 'sna_str_' /
c-----------------------------------------------------------------------
c
      save
     +   ng_tmp,
     +  equg_tmp,
     +  hlp_ssee, hlp_ssei, hlp_sumn, hlp_sume, hlp_sumi, hlp_zeff,
     ,     hlp_zefs, hlp_zefn, hlp_brms, hlp_rdim, hlp_rdnu, hlp_rdtt,
     +  hlp_hg,
     +  nrg_tmp, nr_tmp,
     +  radim_tmp, radnu_tmp, radtt_tmp,
     ,     ssee_tmp,  ssei_tmp,  brms_tmp,
     +  equr_tmp,
     +  hlp_hr,
     +  nsp_tmp, ns_tmp, 
c#ifdef SCL_IMP_RCI
c     + ns1_tmp, equs1_tmp, rrci_tmp, hlp_hs1,
c#endif
     +  equs_tmp,
     +  ssni_tmp,
     +  Natm_tmp,
     +  Na_stra_tmp,
     +  ssna_stra_tmp,
     +  hlp_hs,
     +  na_tmp,
     +  equa_tmp,
     +  ssna_tmp,
     +  hlp_ha,
     +  nf_tmp, nfl_tmp,
     +  equf_tmp,
     +  fluxt_tmp,
     +  hlp_hf,
     +  nb_tmp,equb_tmp,hlp_hb,
     +  lun,
     +  hlp_s,hlp_v,hlp_zs,hlp_zv,hlp_zn1,hlp_zn2,
     +  i, j, k, is, ix, nv,
c#ifdef SCL_IMP_RCI
c    +  iy,
c#endif
     +  ex, warn_tmp,
     +  tab, hlp_c, f, hlp_frm, hlp_frf
c=======================================================================

      if(first_time_step) return

c*** Initialise the arrays
c
      equg_tmp=0.0_R8
      equr_tmp=0.0_R8
      equs_tmp=0.0_R8
c#ifdef SCL_IMP_RCI
c     equs1_tmp=0.0_R8
c#endif
      equa_tmp=0.0_R8
      equf_tmp=0.0_R8
      equb_tmp=0.0_R8
c
c*** Compute the data
c
csw 18nov2013 check alloc
#ifdef B25_EIRENE
      if(allocated(volsumn)) then
#endif
        hlp_sumn=sum(volsumn(1:nstrai))
        hlp_sume=sum(volsumee(1:nstrai))
        hlp_sumi=sum(volsumei(1:nstrai))
#ifdef B25_EIRENE
      endif
#endif

      hlp_brms = 0.0_R8
      hlp_rdim = 0.0_R8
      hlp_rdnu = 0.0_R8
      hlp_rdtt = 0.0_R8
      do k=1,nregionv
        brms_tmp(k)  = b2bremreg(k)
        radim_tmp(k) = b2radreg(k)
        radnu_tmp(k) = rdneureg(k)
        radtt_tmp(k) = radim_tmp(k) + radnu_tmp(k) + brms_tmp(k)
        hlp_brms     = hlp_brms     + brms_tmp(k)
        hlp_rdim     = hlp_rdim     + radim_tmp(k)
        hlp_rdnu     = hlp_rdnu     + radnu_tmp(k)
        hlp_rdtt     = hlp_rdtt     + radtt_tmp(k)
        ssee_tmp(k)  = sseereg(k)
        ssei_tmp(k)  = sseireg(k)
        hlp_ssee     = hlp_ssee     + sseereg(k)
        hlp_ssei     = hlp_ssei     + sseireg(k)
      end do
      do is = 1,nspecies
        do k = 1, nnreg(0)
          ssni_tmp(is)   = ssni_tmp(is)   + ssnireg(is,k)
          ssna_tmp(is,k) = ssna_tmp(is,k) + ssnireg(is,k)
        end do
      end do
      do iatm = 1, nnatmi
#ifdef B25_EIRENE
        is = latmscl(iatm)
        if (is.eq.0) cycle
        Natm_tmp(is) = Natm_tmp(is) + pdena_int(iatm,0)
        do imol = 1, nnmoli
          if (mlcmp(iatm,imol).gt.0)
     &      Natm_tmp(is)   = Natm_tmp(is) +
     &        mlcmp(is,imol) * pdenm_int(imol,0)
        enddo
#endif
      end do

csw 20sep2013 added
      do iy=0,ny-1
        do ix=0,nx-1
          ir=region(ix,iy,0)
          eesum(ir)  = eesum(ir)  + 1.5*te(ix,iy)*ne(ix,iy)*vol(ix,iy)
          eisum(ir)  = eisum(ir)  + 1.5*ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy)
          nesum(ir)  = nesum(ir)  + ne(ix,iy)*vol(ix,iy)
          volsum(ir) = volsum(ir) + vol(ix,iy)
          do is=0,ns-1
            nisum(b2espcr(is),ir) = nisum(b2espcr(is),ir) +
     .                          na(ix,iy,is)*vol(ix,iy)
            kisum(ir) = kisum(ir) +
     &       0.5*mp*am(is)*na(ix,iy,is)*(ua(ix,iy,is)**2)*vol(ix,iy)
            if (is_neutral(is)) then
              nasum(b2espcr(is),ir) = nasum(b2espcr(is),ir) +
     .                                na(ix,iy,is)*vol(ix,iy)
! For fluid neutrals Natm_tmp is just the number of atoms integrated over
! the whole B2.5 domain
              Natm_tmp(b2espcr(is)) = Natm_tmp(b2espcr(is)) +
     .                                na(ix,iy,is)*vol(ix,iy)
            endif
          end do
#ifdef B25_EIRENE
          do iatm = 1, nnatmi
            is = latmscl(iatm)
            if (is.eq.0) cycle
            nasum(is,ir) = nasum(is,ir) +
     .                     dab2(ix+1,iy+1,iatm,1)*vol(ix,iy)
            do imol = 1, nnmoli
              if (mlcmp(iatm,imol).gt.0)
     >          nasum(is,ir) = nasum(is,ir) +
     .           mlcmp(iatm,imol)*dmb2(ix+1,iy+1,imol,1)*vol(ix,iy)
            end do
          end do
#endif
        enddo
      enddo

      do j=1,nstrai
csw 18nov2013 check alloc
#ifdef B25_EIRENE
        if(allocated(srcstrn)) then
#endif
          fluxt_tmp(j) = srcstrn(j)
#ifdef B25_EIRENE
        endif
#endif
      end do
      hlp_s   = 0.0_R8
      hlp_v   = 0.0_R8
      hlp_zs  = 0.0_R8
      hlp_zv  = 0.0_R8
      hlp_zn1 = 0.0_R8
      hlp_zn2 = 0.0_R8
      do ix=0, nx-1
        if(on_closed_surface(ix,0)) then
          hlp_s=hlp_s+gs(ix,0,1)
          hlp_v=hlp_v+vol(ix,0)
          zeff=0.0_R8
          do is=0,ns-1
            zeff = zeff + rz2(ix,0,is)*na(ix,0,is)
          end do
          do is=0,ns_ext-1
            zeff = zeff + za2_ext(ix,0,is)*na_ext(ix,0,is)
          end do
          zeff=zeff/ne(ix,0)
          hlp_zs = hlp_zs + gs(ix,0,1)*zeff
          hlp_zv = hlp_zv + vol(ix,0)*zeff
          do is=0,ns-1
            hlp_zn1 = hlp_zn1 + vol(ix,0)*na(ix,0,is)*rz2(ix,0,is)
            hlp_zn2 = hlp_zn2 + vol(ix,0)*na(ix,0,is)*rza(ix,0,is)
          end do
          do is=0,ns_ext-1
            hlp_zn1 = hlp_zn1+vol(ix,0)*na_ext(ix,0,is)*za2_ext(ix,0,is)
            hlp_zn2 = hlp_zn2+vol(ix,0)*na_ext(ix,0,is)*za_ext(ix,0,is)
          end do
        end if
      end do
      if (hlp_v .gt. 0.0_R8) hlp_zeff = hlp_zv / hlp_v
      if (hlp_s .gt. 0.0_R8) hlp_zefs = hlp_zs / hlp_s
      if (hlp_zn2 .gt. 0.0_R8) hlp_zefn = hlp_zn1 / hlp_zn2

      do i=0,nnatmi
        do j=1,nstrai
          sclf_tmp(i,j) = sclstra(i,j)
          rtio_tmp(i,j) = sclrtio(i,j)-1.0_R8
#ifdef B25_EIRENE
          if (i.gt.0) then
            Na_stra_tmp(i,j) = pdena_int(i,j)
            ssna_stra_tmp(i,j) = ssnastra(i,j)
            do imol = 1, nnmoli
              if (mlcmp(i,imol).gt.0)
     >         Na_stra_tmp(i,j) = Na_stra_tmp(i,j) +
     .         mlcmp(i,imol)*pdenm_int(imol,j)
            end do
          end if
#endif
        end do
      end do
c#ifdef SCL_IMP_RCI
c      do i=1,nspecies !{
c        do ix=0,nx-1 !{
c          do iy=0,ny-1 !{
c            rrci_tmp(i)=rrci_tmp(i)+sclirci(ix,iy,i)*snicrcis(ix,iy,i)
c            do is=nsph(i),nsph(i)-1+nfluids(i) !{
c              rrci_tmp(i)=rrci_tmp(i)+snivrci(ix,iy,is)*na(ix,iy,is)
c            end do !}
c          end do !}
c        end do !}
c      end do !}
c#endif
c
      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then
c
c*** Write the data in a separate file for plotting...
c
        if(nfl_tmp.gt.0) then
          write(lun,hlp_frm,err=94) tim,
     ,          (equg_tmp(j),j=1,ng_tmp),
     ,         ((equr_tmp(i,j),i=1,nrg_tmp),j=1,nr_tmp),
     ,         ((equs_tmp(i,j),i=1,nsp_tmp),j=1,ns_tmp),
     ,        (((equa_tmp(i,j,k),j=1,nrg_tmp),i=1,nsp_tmp),k=1,na_tmp),
     ,         ((equf_tmp(i,j),i=1,nfl_tmp),j=1,nf_tmp),
     ,        (((equb_tmp(i,j,k),j=1,nfl_tmp),i=1,nsp_tmp),k=1,nb_tmp),
     ,         ((equb_tmp(0,j,k),j=1,nfl_tmp),k=1,nb_tmp)
        end if
c#ifdef SCL_IMP_RCI
c     ,        ,((equs1_tmp(i,j),i=1,nsp_tmp),j=1,ns1_tmp)
c#endif
      end if
      return
c=======================================================================
c
      entry b2wrinti(nx,ny,ns,nregionv)
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c-----------------------------------------------------------------------
c
      tab=char(9)
c
c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
c
      inquire(file=f,exist=ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) then
        read(lun,700,err=92) nv,ng_tmp,nr_tmp,ns_tmp,na_tmp,nf_tmp,
     ,                          nrg_tmp,nsp_tmp,nfl_tmp,nb_tmp
c#ifdef SCL_IMP_RCI
c     ,                         ,ns1_tmp
c#endif
        lv_trc(ltrc)=nv
        warn_tmp=.false.

c Check the consistency of existing file and current run

c Check version
        if(nv .ne. lv) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current version ',lv,
     ,                         '  differs from that in the file ',f
c          if(nv .lt. 20010716) then
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc) = .false.
            return
c          end if
        end if

c Check the number of global values
        if(ng_tmp.ne.nglb_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of global values ',
     ,                nglb_tmp,'  differs from that in the file ',f
          if(ng_tmp .gt. nglb_tmp) filedata(ltrc)=.false.
        end if

c Check the number of quantities computed for each region
        if(nr_tmp .ne. nrgn_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of values by region ',
     ,                nrgn_tmp,'  differs from that in the file ',f
          if(nr_tmp .gt. nrgn_tmp) filedata(ltrc) = .false.
        end if

c Check the number of regions
        if(nrg_tmp .ne. nregionv) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of regions ',
     ,                nregionv,'  differs from that in the file ',f
          filedata(ltrc) = .false.
        end if

c Check the number of quantities computed for each species
        if(ns_tmp .ne. nspc_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of values by species',
     ,                nspc_tmp,'  differs from that in the file ',f
          if(ns_tmp .gt. nspc_tmp) filedata(ltrc)=.false.
        end if
c#ifdef SCL_IMP_RCI
c       if(ns1_tmp.ne.nspc1_tmp) then
c         warn_tmp=.true.
c         write(*,*) '!!! b2wrint: current number of values by species',
c     ,          ' (ext1)',nspc1_tmp,'  differs from that in the file ',f
c         if(ns_tmp.gt.nspc_tmp) filedata(ltrc)=.false.
c       end if
c#endif

c Check the number of species
        if(nsp_tmp .ne. nspecies) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of species ',
     ,                nspecies,'  differs from that in the file ',f
          filedata(ltrc) = .false.
        end if

c Check the number quantities computed for each species in each region
        if(na_tmp .ne. nsar_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of values ',
     ,               'by species and region ',
     ,               nsar_tmp,'  differs from that in the file ',f
          if(na_tmp .gt. nsar_tmp) filedata(ltrc) = .false.
        end if

c Check the number of quantities computed for each stratum
        if(nf_tmp .ne. nstr_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of values by strata',
     ,               nstr_tmp,'  differs from that in the file ',f
          if(nf_tmp .gt. nstr_tmp) filedata(ltrc) = .false.
        end if

c Check the number of strata
        if(nfl_tmp .ne. nstrai) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of strata ',
     ,               nfl_tmp,'  differs from that in the file ',f
          filedata(ltrc) = .false.
        end if

c Check the number quantities computed for each species and for each stratum
        if(nb_tmp .ne. nssb_tmp) then
          warn_tmp = .true.
          write(*,*) '!!! b2wrint: current number of values ',
     ,               'by species and strata ',
     ,               nssb_tmp,'  differs from that in the file ',f
          if(nb_tmp .gt. nssb_tmp) filedata(ltrc) = .false.
        end if
c

c If inconsistency is found, generate an error message
        if(.not. filedata(ltrc)) then
          write(*,*) 'No tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if
        if(warn_tmp) write(*,*) '  You should remove the file!'
c*** Read to the end of the file
 50     read(lun,'(a1)',end=55) hlp_c
        go to 50
 55     backspace lun
      else
c*** New file
        ng_tmp  = nglb_tmp
        nr_tmp  = nrgn_tmp
        nrg_tmp = nregionv
        ns_tmp  = nspc_tmp
c#ifdef SCL_IMP_RCI
c       ns1_tmp=nspc1_tmp
c#endif
        nsp_tmp = nspecies
        na_tmp  = nsar_tmp
        nf_tmp  = nstr_tmp
        nfl_tmp = nstrai
        nb_tmp  = nssb_tmp
        write(lun,700) lv,ng_tmp,nr_tmp,ns_tmp,na_tmp,nf_tmp,
     ,                      nrg_tmp,nsp_tmp,nfl_tmp,nb_tmp,
c#ifdef SCL_IMP_RCI
c     ,                     ns1_tmp,
c#endif
     ,                      label

c Create a format string for the column headers
        hlp_frf='(18h(1h#,7x,4htime,4x,,i3,9h(1h'//tab//',a8),,'//
     /          'i3,14h(1h'//tab//',a8,i2.2),,'            //
     /          'i3,12h(1h'//tab//',a8,a2),,'              //
     /          'i3,17h(1h'//tab//',a8,i2.2,a2),,'         //
     /          'i3,18h(1h'//tab//',a8,1h_,i2.2),,'        //
     /          'i3,25h(1h'//tab//',a8,1h_,i2.2,1h_,a2),,' //
     /          'i3,22h(1h'//tab//',a8,1h_,i2.2,1h_)'      //
c#ifdef SCL_IMP_RCI
c     /         ',,i3,12h(1h'//tab//',a8,a2)'              //
c#endif
     /          '))'

        if (nfl_tmp.gt.0) then
         write(hlp_frm,hlp_frf) ng_tmp, nr_tmp*nrg_tmp, ns_tmp*nsp_tmp,
     ,                          na_tmp*nrg_tmp*nsp_tmp, nf_tmp*nfl_tmp,
     ,                          nb_tmp*nsp_tmp*nfl_tmp, nb_tmp*nfl_tmp
        else
          write(*,*)
     &     'Unexpected inconsistency: number of strata NSTRAI = ',
     &      nfl_tmp
          write(*,*) 'Impossible to write into ', f
        endif
c#ifdef SCL_IMP_RCI
c     ,                        ,ns1_tmp*nsp_tmp
c#endif

c<<<
c        print *,'%%%b2wrint: nb_tmp,nsp_tmp,nfl_tmp= ',
c     ,                                           nb_tmp,nsp_tmp,nfl_tmp
c        print *,'hlp_frm:'
c        print *,hlp_frm
c>>>

        if (nfl_tmp.gt.0) then
          write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp),
     ,          ((hlp_hr(j),i,i=1,nrg_tmp),j=1,nr_tmp),
     ,          ((hlp_hs(j),species_list(i),i=1,nsp_tmp),j=1,ns_tmp),
     ,         (((hlp_ha(k),j,species_list(i),j=1,nrg_tmp),
     ,                                        i=1,nsp_tmp),k=1,na_tmp),
     ,          ((hlp_hf(j),i,i=1,nfl_tmp),j=1,nf_tmp),
     ,         (((hlp_hb(k),j,species_list(i),j=1,nfl_tmp),
     ,                                        i=1,nsp_tmp),k=1,nb_tmp),
     ,          ((hlp_hb(k),j,j=1,nfl_tmp),k=1,nb_tmp)
        end if
c#ifdef SCL_IMP_RCI
c     ,        ,((hlp_hs1(j),species_list(i),i=1,nsp_tmp),j=1,ns1_tmp)
c#endif

      end if

      if (nfl_tmp.gt.0) then
        write(hlp_frm,'(10h(1p,e16.7,,i3,12h(1h'//tab//',e12.4)))')
     -      ng_tmp + nr_tmp*nrg_tmp + ns_tmp*nsp_tmp +
     +        na_tmp*nrg_tmp*nsp_tmp  + nf_tmp*nfl_tmp +
     +          (nsp_tmp+1)*nfl_tmp*nb_tmp
      end if
c#ifdef SCL_IMP_RCI
c     +          + ns1_tmp*nsp_tmp
c#endif
c      print *,'hlp_frm:' !###
c      print *,hlp_frm    !###

      return
c=======================================================================
 90   write(*,*) 'b2wrint: cannot open a file ',f
      filedata(ltrc)=.false.
      return
c
 92   write(*,*) 'b2wrint: wrong format of the file ',f
      filedata(ltrc)=.false.
      return
c
 94   write(*,*) 'b2wrint: cannot write in the file ',f
      filedata(ltrc)=.false.
      return
c=======================================================================
      entry b2wrints(nx,ny,ns,nregionv)
c-----------------------------------------------------------------------
c*** Save the values (not used for this set)
c-----------------------------------------------------------------------
c
      return
c=======================================================================
      entry b2wrintf(nx,ny,ns,nregionv)
c-----------------------------------------------------------------------
c*** Finalization step
c-----------------------------------------------------------------------
c
      return
c=======================================================================
c#ifndef SCL_IMP_RCI
  700 format('#   ',10i8/'##  ',a)
c#else
c  700 format('#   ',11i8/'##  ',a)
c#endif
c=======================================================================
      end subroutine b2wrint
