      module b2mod_wrsep
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_indirect
      implicit none

c======================================================================
c*** Puts some separatrix data on the output file "tracing/sepdata.trc"
c***  - primarily, to provide boundary conditions for core models
c----------------------------------------------------------------------

c----------------------------------------------------------------------
c*** lv    :  version number to appear in the file

      integer lv
      parameter (lv=20000828)
      integer ltrc
      parameter (ltrc=ntrc_sep)

c----------------------------------------------------------------------
c*** Separatrix data:
c***
c***      header     var       Legend
c***...........................................................(global)
c***     ne_sep     hlp_nesp    average ne
c***     Te_sep     hlp_tesp    density-weighted average Te
c***     Ti_sep     hlp_tisp    density-weighted average Ti
c***     Zef_sep    hlp_zesp    density-weighted average Zeff
c***     Qt_sep     hlp_qtsp    total heat flow from the core
c***     Qe_sep     hlp_qesp    electron heat flow from the core
c***     Qi_sep     hlp_qisp    ion heat flow from the core
c***     Qrad_SOL   hlp_qrsl    radiation from SOL and divertors
c***     Spmp_nom   hlp_spmn    nominal pumping speed
c***
c***..........................................................(species)
c***     ni_sep_    nisp_tmp    average ni for each species
c***     Zi_sep_    zisp_tmp    average Z for each ion species
c***     Ze_sep_    zesp_tmp    average Zeff for each ion species
c***     Gi_sep_    gisp_tmp    outward ion flux for each species
c***     na_sep_    nasp_tmp    average natm for each species
c***     Ta_sep_    tasp_tmp    average Tatm for each species
c***     Ga_sep_    gasp_tmp    inward neutral flux for each species
c***     Gp_sep_    gpuf_tmp    gas puffing rate for each species
c***     pr_sep_    prsr_tmp    relative total pressure for each species
c***     nr_sep_    nrsr_tmp    relative ion density for each species
c***
c***.............................................................(ions)
c***     nii_sep_   nisi_tmp    average ni for each ion
c***     Gii_sep_   gisi_tmp    outward ion flux for each ion
c----------------------------------------------------------------------

      integer, save :: ng_tmp, nsp_tmp, ns_tmp, ni_tmp, nii_tmp, nis_tmp
      integer, save :: jxi, jxa, jsep, lun
      integer, save, allocatable :: ii(:)
      character*2, save, allocatable :: hlp_a(:)
      character, save :: f*(lenfnm), hlp_frm*256

      integer nglb_tmp
      parameter (nglb_tmp=9)
      integer nspc_tmp
      parameter (nspc_tmp=10)
      integer nion_tmp
      parameter (nion_tmp=2)
      character, save :: hlp_hg(nglb_tmp)*8
      character, save :: hlp_hs(nspc_tmp)*7
      character, save :: hlp_hi(nion_tmp)*8
c*** Values by species
      real (kind=R8), allocatable ::
     *      equs_tmp(:,:)
      real (kind=R8), allocatable ::
     *      nisp_tmp(:), zisp_tmp(:), zesp_tmp(:),
     ,      gisp_tmp(:), nasp_tmp(:), tasp_tmp(:),
     ,      gasp_tmp(:), gpuf_tmp(:),
     ,      prsr_tmp(:), nrsr_tmp(:)

c*** Values by ion charge states
      real (kind=R8), allocatable ::
     *     equi_tmp(:,:)
      real (kind=R8), allocatable ::
     *     nisi_tmp(:), gisi_tmp(:)
c----------------------------------------------------------------------
c*** Column headers:
      data hlp_hg /
     1 '  ne_sep','  Te_sep','  Ti_sep',' Zef_sep',
     2 '  Qt_sep','  Qe_sep','  Qi_sep','Qrad_SOL',
     3 'Spmp_nom'/
      data hlp_hs /
     1  'ni_sep_','Zi_sep_','Ze_sep_','Gi_sep_','na_sep_',
     2  'Ta_sep_','Ga_sep_','Gp_sep_','pr_sep_','nr_sep_' /
      data hlp_hi / 'nii_sep_','Gii_sep_' /
c----------------------------------------------------------------------

      contains

      subroutine b2wrsep(nx,ny,ns,BoRiS)
      use b2mod_geo
      use b2mod_rates
      use b2mod_plasma
      use b2mod_sources
     . , only : sna
      use b2mod_tallies
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_external
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_eirbra
#else
      use b2mod_eirbra
#endif
      implicit none

      integer nx, ny, ns
      real (kind=R8) :: BoRiS
c
c*** Global values
      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     hlp_nesp, hlp_tesp, hlp_tisp, hlp_zesp, hlp_qtsp, hlp_qesp,
     ,     hlp_qisp, hlp_qrsl, hlp_spmn
      equivalence
     ,     (equg_tmp(1),  hlp_nesp), (equg_tmp(2),  hlp_tesp),
     ,     (equg_tmp(3),  hlp_tisp), (equg_tmp(4),  hlp_zesp),
     ,     (equg_tmp(5),  hlp_qtsp), (equg_tmp(6),  hlp_qesp),
     ,     (equg_tmp(7),  hlp_qisp), (equg_tmp(8),  hlp_qrsl),
     ,     (equg_tmp(9),  hlp_spmn)

c*** Local variables of general use
      real (kind=R8) ::
     *     v,hlp_i,hlp_n,hlp_p,hlp_v,zeff(-1:nx,-1:ny)
      integer i,j,k,l,ia,is,ix,iy

      save
c======================================================================

      if(first_time_step) return

      if(ank_tracing.eq.0) return

c*** Initialise the arrays
c
      equg_tmp=0.0_R8
      do i=1,nspecies !{
        nisp_tmp(i)=0.0_R8
        zisp_tmp(i)=0.0_R8
        zesp_tmp(i)=0.0_R8
        gisp_tmp(i)=0.0_R8
        nasp_tmp(i)=0.0_R8
        tasp_tmp(i)=0.0_R8
        gasp_tmp(i)=0.0_R8
        gpuf_tmp(i)=0.0_R8
        prsr_tmp(i)=0.0_R8
        nrsr_tmp(i)=0.0_R8
      end do !}
      do i=1,ns !{
        nisi_tmp(i)=0.0_R8
        gisi_tmp(i)=0.0_R8
      end do !}
      hlp_n=0.0_R8
      hlp_p=0.0_R8
      hlp_v=0.0_R8
c
c*** Compute the data
c
      call b2xzef_st(nx,ny,ns,rz2,na,ne,zeff)
      iy=jsep+1
      do ix=0,nx-1 !{
        if(on_closed_surface(ix,jsep)) then !{
c          v=vol(ix,iy)
          v=gs(ix,iy,1)
          hlp_i=0.0_R8
          do is=0,ns-1 !{
            if (.not.is_neutral(is)) hlp_i=hlp_i+na(ix,iy,is)
          end do !}
          hlp_v=hlp_v+v
          hlp_n=hlp_n+hlp_i*v
          hlp_nesp=hlp_nesp+v*ne(ix,iy)
          hlp_tesp=hlp_tesp+v*ne(ix,iy)*te(ix,iy)
          hlp_p=hlp_p+v*(ne(ix,iy)*te(ix,iy)+hlp_i*ti(ix,iy))
          do is = 0, ns_ext-1
            hlp_p=hlp_p+v*na_ext(ix,iy,is)*ta_ext(ix,iy,is)
            hlp_tisp=hlp_tisp+v*na_ext(ix,iy,is)*ta_ext(ix,iy,is)
          end do
          hlp_tisp=hlp_tisp+v*hlp_i*ti(ix,iy)
          hlp_zesp=hlp_zesp+v*ne(ix,iy)*zeff(ix,iy)
          hlp_qesp=hlp_qesp+fhe(ix,iy,1,1)
          hlp_qisp=hlp_qisp+fhi(ix,iy,1,1)
          do is=0,ns-1
            hlp_qisp=hlp_qisp+(1.0_R8-BoRiS)*fhm(ix,iy,1,1,is)
          enddo
          is=0
          do j=1,nspecies !{
            do k=1,nfluids(j) !{
              is=is+1
              nisi_tmp(is)=nisi_tmp(is)+v*na(ix,iy,ions_list(is))
              nisp_tmp(j)=nisp_tmp(j)+v*na(ix,iy,ions_list(is))
              prsr_tmp(j)=prsr_tmp(j)+v*na(ix,iy,ions_list(is))*
     *            (ti(ix,iy)+rza(ix,iy,ions_list(is))*te(ix,iy))
              zisp_tmp(j)=zisp_tmp(j)+
     *             v*na(ix,iy,ions_list(is))*rza(ix,iy,ions_list(is))
              zesp_tmp(j)=zesp_tmp(j)+
     *             v*na(ix,iy,ions_list(is))*rz2(ix,iy,ions_list(is))
              gisi_tmp(is)=gisi_tmp(is)+fna_mdf(ix,iy,1,1,ions_list(is))
              gisp_tmp(j)=gisp_tmp(j)+fna_mdf(ix,iy,1,1,ions_list(is))
              do i=0,jsep !{
                gasp_tmp(j)=gasp_tmp(j)+
     1               sna(ix,i,0,ions_list(is))+
     2               sna(ix,i,1,ions_list(is))*na(ix,i,ions_list(is))
              end do !}
            end do !}
#ifdef WG_TODO
            do k=1,nnatmi !{
              if (latmscl(k).eq.j) then !{
                nasp_tmp(j)=nasp_tmp(j)+v*dab2(ix+1,iy+1,k,1)
                tasp_tmp(j)=tasp_tmp(j)+v*dab2(ix+1,iy+1,k,1)
     .                                   *tab2(ix+1,iy+1,k,1)
              end if !}
            end do !}
#endif
          end do !}
        end if !}
      end do !}

      do i=2,nnreg(0) !{
        if(nnreg(0).eq.8.and.(i.eq.5.or.i.eq.6)) cycle
        hlp_qrsl=hlp_qrsl-(b2radreg(i)+rdneureg(i)+b2bremreg(i))
      end do !}

      is=0
      do j=1,nspecies !{
        if(hlp_nesp.gt.0.0_R8) then !{
          zesp_tmp(j)=zesp_tmp(j)/hlp_nesp
        else !}{
          zesp_tmp(j)=0.0_R8
        end if !}
        if(nisp_tmp(j).gt.0.0_R8) then !{
          zisp_tmp(j)=zisp_tmp(j)/nisp_tmp(j)
        else !}{
          zisp_tmp(j)=0.0_R8
        end if !}
        if(nasp_tmp(j).gt.0.0_R8) then !{
          tasp_tmp(j)=tasp_tmp(j)/nasp_tmp(j)/qe
        else !}{
          tasp_tmp(j)=0.0_R8
        end if !}
        if(hlp_v.gt.0.0_R8) then !{
          nisp_tmp(j)=nisp_tmp(j)/hlp_v
          nasp_tmp(j)=nasp_tmp(j)/hlp_v
        else !}{
          nisp_tmp(j)=0.0_R8
          nasp_tmp(j)=0.0_R8
        end if !}
        if(hlp_p.gt.0.0_R8) then !{
          prsr_tmp(j)=prsr_tmp(j)/hlp_p
        else !}{
          prsr_tmp(j)=0.0_R8
        end if !}
        if(hlp_v.gt.0.0_R8) then !{
          do i=1,nfluids(j) !{
            is=is+1
            nisi_tmp(is)=nisi_tmp(is)/hlp_v
          end do !}
        else !}{
          do i=1,nfluids(j) !{
            is=is+1
            nisi_tmp(is)=0.0_R8
          end do !}
        end if !}
      end do !}

      if(hlp_nesp.gt.0.0_R8) then !{
        hlp_tesp=hlp_tesp/hlp_nesp/qe
        hlp_zesp=hlp_zesp/hlp_nesp
      else !}{
        hlp_tesp=0.0_R8
        hlp_zesp=0.0_R8
      end if !}
      if(hlp_n.gt.0.0_R8) then !{
        hlp_tisp=hlp_tisp/hlp_n/qe
      else !}{
        hlp_tisp=0.0_R8
      end if !}
      if(hlp_v.gt.0.0_R8) then !{
        hlp_nesp=hlp_nesp/hlp_v
      else !}{
        hlp_nesp=0.0_R8
      end if !}
      hlp_qtsp=hlp_qisp+hlp_qesp
      hlp_spmn=Spmp_nom

      do is=1,nstrat !{
        if(crcstra(is).eq.'C') then !{
          v=0.0_R8
          do i=1,nnatmi !{
            v=v+gpfc(i,is)
          end do !}
          do i=1,nspecies !{
            do j=1,nnatmi !{
              if(latmscl(j).eq.i.and.v.gt.0.0_R8)
     -          gpuf_tmp(i)=gpuf_tmp(i)+userfluxparm(is,1)*gpfc(j,is)/v
            end do !}
          end do !}
        end if !}
      end do !}

c      print *,'@@gasp= ',(gasp_tmp(i),i=1,nspecies)              !###
      do k=1,ncns !{
        l=nnlimi+lcns(k)
        do is=1,nspecies !{
          do ia=1,nnatmi !{
            if (latmscl(ia).eq.is) then !{
c             print *,'@@gasp: k,l,is,ia,wldna=',                 !###
c     ,                        k,l,is,ia,wldna(l,ia,0)            !###
              gasp_tmp(is)=gasp_tmp(is)+wldna(l,ia,0)
              do i=1,nnmoli !{
c               if(mlcmp(ia,i).gt.0) print *,'@@i,wldnm,mlcmp=',  !###
c     ,                                i,wldnm(l,i,0),mlcmp(ia,i) !###
                gasp_tmp(is)=gasp_tmp(is)+wldnm(l,i,0)*mlcmp(ia,i)
              end do !}
            end if !}
          end do !}
        end do !}
      end do !}

      if(hlp_nesp.gt.0.0_R8) then
        do k=1,nspecies !{
          nrsr_tmp(k)=nisp_tmp(k)/hlp_nesp
        end do !}
      else !}{
        do k=1,nspecies !{
          nrsr_tmp(k)=0.0_R8
        end do !}
      end if !}

      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then !{

c*** Write the data in a separate file for plotting...
        do k=1, nspecies !{
          equs_tmp(k,1) = nisp_tmp(k)
          equs_tmp(k,2) = zisp_tmp(k)
          equs_tmp(k,3) = zesp_tmp(k)
          equs_tmp(k,4) = gisp_tmp(k)
          equs_tmp(k,5) = nasp_tmp(k)
          equs_tmp(k,6) = tasp_tmp(k)
          equs_tmp(k,7) = gasp_tmp(k)
          equs_tmp(k,8) = gpuf_tmp(k)
          equs_tmp(k,9) = prsr_tmp(k)
          equs_tmp(k,10)= nrsr_tmp(k)
        end do !}
        do k = 1, ns !{
          equi_tmp(k,1) = nisi_tmp(k)
          equi_tmp(k,2) = gisi_tmp(k)
        end do !}

        write(lun,hlp_frm,err=94) tim,
     ,          (equg_tmp(j),j=1,ng_tmp),
     ,         ((equs_tmp(i,j),i=1,nsp_tmp),j=1,ns_tmp),
     ,         ((equi_tmp(i,j),i=1,nii_tmp),j=1,ni_tmp)
      end if !}
      return
c
   94 write(*,*) 'b2wrsep: cannot write in the file ',trim(f)
      filedata(ltrc)=.false.
      return
      end subroutine b2wrsep
c======================================================================

      subroutine b2wrsepi(nx,ny,ns,switch)
      use b2mod_switches
      implicit none
      integer, intent(in) :: nx, ny, ns
      type (switches), intent (in) :: switch
      integer i, j, k, nv
      character hlp_c
      logical ex, warn_tmp
c----------------------------------------------------------------------
c*** Initialise the file for tracing data
c----------------------------------------------------------------------

      call get_jsep(nx,ny,jxi,jxa,jsep)

      if(nspecies.le.0) then !{
        write(*,*) '!!! b2wrsep: nspecies.lt.0 - no tracing data'
        filedata(ltrc)=.false.
      else !}{
        nis_tmp=0
        do i=1,nspecies !{
          nis_tmp=nis_tmp+nfluids(i)
          if(nfluids(i).le.0) then !{
            write(*,*) '!!! b2wrsep: nfluids(',i,' ) = ',nfluids(i),
     ,                                            '  - no tracing data'
            filedata(ltrc)=.false.
          end if !}
        end do !}
      end if !}
      if(nncut.ne.1 .and. nncut.ne.2) then !{
        write(*,*) '!!! b2wrsep: nncut = ',nncut,' - no tracing data'
        filedata(ltrc)=.false.
      end if !}

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)

      allocate(ii(ns))
      allocate(hlp_a(ns))
      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)
      if(ex.and.switch%b2mndr_stim.lt.0.0_R8) then
        read(lun,*,err=92) hlp_c,
     ,   nv,ng_tmp,ns_tmp,ni_tmp,nsp_tmp,nii_tmp
        lv_trc(ltrc)=nv
        warn_tmp=.false.
        if(nv.ne.lv) then
          warn_tmp=.true.
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2wrsep: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
c          if(nv .lt. 970212) then
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            return
c          end if
        end if
        if(ng_tmp.ne.nglb_tmp) then
          warn_tmp=.true.
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2wrsep: current number of global values ',nglb_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ng_tmp
          if(ng_tmp.gt.nglb_tmp) filedata(ltrc)=.false.
        end if
        if(ns_tmp.ne.nspc_tmp) then
          warn_tmp=.true.
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2wrsep: current number of values by species ',nspc_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ns_tmp
          if(ns_tmp.gt.nspc_tmp) filedata(ltrc)=.false.
        end if
        if(nsp_tmp.ne.nspecies) then
          warn_tmp=.true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2wrsep: current number of species ',nspecies,
     ,     ' differs from that in the file ',trim(f),' : ',nsp_tmp
          filedata(ltrc)=.false.
        end if
        if(ni_tmp.ne.nion_tmp) then
          warn_tmp=.true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2wrsep: current number of values',
     ,     ' by ion charge state ',nion_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ni_tmp
          if(ni_tmp.gt.nion_tmp) filedata(ltrc)=.false.
        end if
        if(nii_tmp.ne.nis_tmp) then
          warn_tmp=.true.
          write(*,'(a,i3,3a,i3)')
     ,     '!!! b2wrsep: current number of ion charge states ',
     ,     ns-nspecies,
     ,     ' differs from that in the file ',trim(f),' : ',nii_tmp
          filedata(ltrc)=.false.
        end if

        if(.not.filedata(ltrc)) then
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if
        if(warn_tmp) write(*,*) ' You should remove the file!'
c*** Read to the end of the file
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun
      else
c*** New file
        ng_tmp=nglb_tmp
        ns_tmp=nspc_tmp
        nsp_tmp=nspecies
        ni_tmp=nion_tmp
        nii_tmp=nis_tmp
        write(lun,700) lv,ng_tmp,ns_tmp,ni_tmp,nsp_tmp,nii_tmp
        write(lun,'(2a)') '#  ',trim(label)

        write(hlp_frm,'(a,i3,a,i3,a,i3,a)')
     ,   '(a1,6x,a4,5x,',ng_tmp,'(5x,a8,1x),',
     ,                   ns_tmp*nsp_tmp,'(5x,a7,a2),',
     ,                   ni_tmp*nii_tmp,'(1x,a8,a2,i2.2,''+''))'
c<<<
c        print *,'%%%b2wrsep: nb_tmp,nsp_tmp,nfl_tmp= ',
c     ,                       nb_tmp,nsp_tmp,nfl_tmp
c        print *,'hlp_frm:'
c        print *,trim(hlp_frm)
c>>>
        k=0
        do i=1,nspecies !{
          do j=1,nfluids(i) !{
            k=k+1
            ii(k)=j
            hlp_a(k)=species_list(i)
            if(hlp_a(k)(2:2).eq.' ') hlp_a(k)(2:2)='_'
          end do !}
        end do !}
        write(lun,hlp_frm) '#','time',
     ,    (hlp_hg(i),i=1,ng_tmp),
     ,   ((hlp_hs(j),species_list(i),i=1,nsp_tmp),j=1,ns_tmp),
     ,   ((hlp_hi(k),hlp_a(i),ii(i),i=1,nii_tmp),k=1,ni_tmp)

      end if
      write(hlp_frm,'(a,i3,a)') '(1p,e16.7,',
     -  ng_tmp + ns_tmp*nsp_tmp + ni_tmp*nii_tmp,'(2x,e12.4))'

      allocate(equs_tmp(nspecies,nspc_tmp))
      allocate(nisp_tmp(nspecies))
      allocate(zisp_tmp(nspecies))
      allocate(zesp_tmp(nspecies))
      allocate(gisp_tmp(nspecies))
      allocate(nasp_tmp(nspecies))
      allocate(tasp_tmp(nspecies))
      allocate(gasp_tmp(nspecies))
      allocate(gpuf_tmp(nspecies))
      allocate(prsr_tmp(nspecies))
      allocate(nrsr_tmp(nspecies))
      allocate(equi_tmp(ns,nion_tmp))
      allocate(nisi_tmp(ns))
      allocate(gisi_tmp(ns))
c<<<
c      print *,'hlp_frm:'
c      print *,trim(hlp_frm)
c>>>
      return
c======================================================================
   90 write(*,*) 'b2wrsep: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2wrsept: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return
c======================================================================
  700 format('#   ',6i8)
      end subroutine b2wrsepi
c======================================================================
c
      subroutine b2wrsepc(iline)
      implicit none
      integer, intent(out) :: iline
      logical ex
      real (kind=R8) :: tdummy
      character hlp_c
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_sep)) return
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)
      open(lun,file=f,err=190)
      rewind(lun)
      iline = 0
      read(lun,'(a1)',err=155,end=155) hlp_c
      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',err=155,end=155) hlp_c
      end do
      backspace lun
      read(lun,*,err=155,end=155) tdummy
      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,err=155,end=155) tdummy
      end do
  155 continue
!      close(lun)
      return
c======================================================================
  190 write(*,*) 'b2wrsep: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
      end subroutine b2wrsepc
c======================================================================
      subroutine b2wrseps
c----------------------------------------------------------------------
c
      return
      end subroutine b2wrseps
c======================================================================
      subroutine b2wrsepf
      implicit none
c----------------------------------------------------------------------
c     Finalization step
      logical file_open

      lun=ntrc(ltrc)
      inquire(unit=lun,opened=file_open)
      if (file_open) close(lun)

      if (.not.allocated(equs_tmp)) return
      deallocate(equs_tmp)
      deallocate(nisp_tmp)
      deallocate(zisp_tmp)
      deallocate(zesp_tmp)
      deallocate(gisp_tmp)
      deallocate(nasp_tmp)
      deallocate(tasp_tmp)
      deallocate(gasp_tmp)
      deallocate(gpuf_tmp)
      deallocate(prsr_tmp)
      deallocate(nrsr_tmp)
      deallocate(equi_tmp)
      deallocate(nisi_tmp)
      deallocate(gisi_tmp)
      deallocate(ii)
      deallocate(hlp_a)
      return
      end subroutine b2wrsepf
c======================================================================
      end module b2mod_wrsep
