      module b2mod_wrsrt
      use b2mod_types
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      implicit none

c========================================================================
c*** Puts some integral quantities on the output file
c*** "tracing/intshrt.trc"
c-----------------------------------------------------------------------
c
c*** lv    :  version number to appear in the file
c
      integer lv
      parameter (lv=20171020)
      integer ltrc
      parameter (ltrc=ntrc_srt)
      integer, save :: lun

c-----------------------------------------------------------------------
c*** Additional integrals:
c***
c***..........................................................(regions)
c***     totenre_  tesum      total plasma energy \int 3/2 n_e T_e dV plus
c***                          \int 3/2 n_i T_i dV plus \int 1/2 m_a n_a u_{a||}^2 dV
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

      integer, save :: nrg_tmp, nr_tmp, nsp_tmp, ns_tmp, na_tmp, nfl_tmp
      integer nrgn_tmp
      parameter (nrgn_tmp=1) ! Number of quantities computed in each region
      integer nspc_tmp
      parameter (nspc_tmp=1)
      integer nsar_tmp
      parameter (nsar_tmp=2)
      character, save :: hlp_hr(nrgn_tmp)*8
      character, save :: hlp_hs(nspc_tmp)*8
      character, save :: hlp_ha(nsar_tmp)*8

      character, save :: f*(lenfnm), hlp_frm*256
c-----------------------------------------------------------------------
c*** Column headers:
      data hlp_hr / 'totenre_' /
      data hlp_hs / 'N_a_tot_' /
      data hlp_ha / 'N_i_reg_', 'N_a_reg_' /
c-----------------------------------------------------------------------

      contains

      subroutine b2wrsrt(ns,geo,mpg,st)
      use b2mod_constants
      use b2mod_numerics_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_eirbra
#else
      use b2mod_eirbra
#endif
      use b2mod_dimensions
      use b2us_map
      use b2us_geo
      use b2us_plasma
      implicit none
c
      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(B2state), intent(in) :: st
c
c*** Values by regions
      real (kind=R8) :: tesum(2*DEF_NCUT)

      real (kind=R8) :: equr_tmp(2*DEF_NCUT,nrgn_tmp)
      equivalence (equr_tmp(1,1), tesum)
c
c*** Values by species
      real (kind=R8) :: equs_tmp(DEF_NATM,nspc_tmp)
      real (kind=R8) :: Natm_tmp(DEF_NATM)
      equivalence (equs_tmp(1,1), Natm_tmp)

c*** Values by species and regions
      real (kind=R8) :: equa_tmp(DEF_NATM,2*DEF_NCUT,nsar_tmp)
      real (kind=R8) ::
     *     nasum(DEF_NATM,2*DEF_NCUT), nisum(DEF_NATM,2*DEF_NCUT)     ! IYS 29.09.2017 added total number of neutrals
      equivalence (equa_tmp(1,1,1), nisum), (equa_tmp(1,1,2), nasum)  ! IYS 29.09.2017
c-----------------------------------------------------------------------
c*** Local variables of general use
      integer i, j, k, ir, is, iCv
      integer iatm ! IYS 29.09.2017
#ifdef B25_EIRENE
      integer imol
#endif
c
      save
     +  equr_tmp,
     +  equs_tmp,
     +  Natm_tmp,
     +  equa_tmp,
     +  i, j, k, is
c=======================================================================

c*** Initialise the arrays
c
      equr_tmp=0.0_R8
      equs_tmp=0.0_R8
      equa_tmp=0.0_R8
c
c*** Compute the data
c
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
      do iCv=1,mpg%nCv
        ir = mpg%cvReg(iCv)
        tesum(ir)  = tesum(ir)
     &          + 1.5_R8*st%pl%te(iCv)*st%dv%ne(iCv)*geo%cvVol(iCv)
     &          + 1.5_R8*st%pl%ti(iCv)*st%dv%ni(iCv,0)*geo%cvVol(iCv)
        do is=0,ns-1
          nisum(b2espcr(is),ir) = nisum(b2espcr(is),ir) +
     .                        st%pl%na(iCv,is)*geo%cvVol(iCv)
          tesum(ir) = tesum(ir) + 0.5_R8*mp*am(is)*
     &       st%pl%na(iCv,is)*(st%pl%ua(iCv,is)**2)*geo%cvVol(iCv)
          if (is_neutral(is)) then
            nasum(b2espcr(is),ir) = nasum(b2espcr(is),ir) +
     .                              st%pl%na(iCv,is)*geo%cvVol(iCv)
! For fluid neutrals Natm_tmp is just the number of atoms integrated over
! the whole B2.5 domain
            Natm_tmp(b2espcr(is)) = Natm_tmp(b2espcr(is)) +
     .                              st%pl%na(iCv,is)*geo%cvVol(iCv)
          endif
        end do
#ifdef B25_EIRENE
        do iatm = 1, nnatmi
          is = latmscl(iatm)
          if (is.eq.0) cycle
          nasum(is,ir) = nasum(is,ir) +
     .                   dab2(iCv,iatm,1)*geo%cvVol(iCv)
          do imol = 1, nnmoli
            if (mlcmp(iatm,imol).gt.0)
     >       nasum(is,ir) = nasum(is,ir) +
     .        mlcmp(iatm,imol)*dmb2(iCv,imol,1)*geo%cvVol(iCv)
          end do
        end do
#endif
      enddo

      if(filedata(ltrc)) then
c
c*** Write the data in a separate file for plotting...
c
        if(nfl_tmp.gt.0) then
          write(lun,hlp_frm,err=94) tim,
     ,         ((equr_tmp(i,j),i=1,nrg_tmp),j=1,nr_tmp),
     ,         ((equs_tmp(i,j),i=1,nsp_tmp),j=1,ns_tmp),
     ,        (((equa_tmp(i,j,k),j=1,nrg_tmp),i=1,nsp_tmp),k=1,na_tmp)
        end if
      end if
      return
c
   94 write(*,*) 'b2wrsrt: cannot write in the file ',trim(f)
      filedata(ltrc)=.false.
      return
      end subroutine b2wrsrt
c=======================================================================
c
      subroutine b2wrsrti(nregionv,switch)
      use b2mod_switches
      implicit none
      integer, intent(in) :: nregionv
      type(switches), intent(in) :: switch
      integer i, j, k, nv
      character hlp_c
      logical ex, warn_tmp
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c-----------------------------------------------------------------------
c
c
c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)

!      write(*,*) 'Inside b2wrsrt: lun = ',lun
!      write(*,*) 'Inside b2wrsrt: f = ', trim(f)

c
      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)

      if(ex.and.switch%b2mndr_stim.lt.0.0_R8) then ! the old file exists and is found
                                     ! and is to be appended
        read(lun,*,err=92) hlp_c, nv,
     &     nr_tmp, ns_tmp, na_tmp,
     ,     nrg_tmp, nsp_tmp, nfl_tmp
        lv_trc(ltrc)=nv
        warn_tmp=.false.

c Check the consistency of existing file and current run

c Check version
        if(nv .ne. lv) then
          warn_tmp = .true.
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2wrsrt: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
c*** It was a major change of the file format - impossible to continue
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          filedata(ltrc) = .false.
          return
        end if

c Check the number of quantities computed for each region
        if(nr_tmp .ne. nrgn_tmp) then
          warn_tmp = .true.
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2wrsrt: current number of values by region ',nrgn_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',nr_tmp
          if(nr_tmp .gt. nrgn_tmp) filedata(ltrc) = .false.
        end if

c Check the number of regions
        if(nrg_tmp .ne. nregionv) then
          warn_tmp = .true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2wrsrt: current number of regions ',nregionv,
     ,     ' differs from that in the file ',trim(f),' : ',nrg_tmp
          filedata(ltrc) = .false.
        end if

c Check the number of quantities computed for each species
        if(ns_tmp .ne. nspc_tmp) then
          warn_tmp = .true.
          write(*,'(a,i4,3a,i4)')
     ,     '!!! b2wrsrt: current number of values by species ',nspc_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',ns_tmp
          if(ns_tmp .gt. nspc_tmp) filedata(ltrc)=.false.
        end if

c Check the number of species
        if(nsp_tmp .ne. nspecies) then
          warn_tmp = .true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2wrsrt: current number of species ',nspecies,
     ,     ' differs from that in the file ',trim(f),' : ',nsp_tmp
          filedata(ltrc) = .false.
        end if

c Check the number quantities computed for each species in each region
        if(na_tmp .ne. nsar_tmp) then
          warn_tmp = .true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2wrsrt: current number of values',
     ,     ' by species and region ',nsar_tmp,
     ,     ' differs from that in the file ',trim(f),' : ',na_tmp
          if(na_tmp .gt. nsar_tmp) filedata(ltrc) = .false.
        end if

c Check the number of strata
        if(nfl_tmp .ne. nstrat) then
          warn_tmp = .true.
          write(*,'(a,i3,3a,i3)')
     ,     '!!! b2wrsrt: current number of strata ',nstrat,
     ,     ' differs from that in the file ',trim(f),' : ',nfl_tmp
          filedata(ltrc) = .false.
        end if

c If inconsistency is found, generate an error message
        if(.not. filedata(ltrc)) then
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if

c However, there exist two level of inconsistency strength!..
        if(warn_tmp) write(*,*) ' You should remove the file!'

c*** Read to the end of the file
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun

      else
c*** New file
        nr_tmp  = nrgn_tmp
        nrg_tmp = nregionv
        ns_tmp  = nspc_tmp
        nsp_tmp = nspecies
        na_tmp  = nsar_tmp
        nfl_tmp = nstrat
        write(lun,700) lv,nr_tmp,ns_tmp,na_tmp,
     ,                      nrg_tmp,nsp_tmp,nfl_tmp,
     ,                      label

c Create a format string for the column headers
        if (nfl_tmp.gt.0) then
          write(hlp_frm,'(a,i3,a,i3,a,i3,a)') '(a1,6x,a4,4x,',
     ,             nr_tmp*nrg_tmp,'(4x,a8,i2.2),',
     ,             ns_tmp*nsp_tmp,'(4x,a8,a2),1x,',
     ,     na_tmp*nrg_tmp*nsp_tmp,'(1x,a8,i2.2,''_'',a2))'
#ifdef DBG
c        print *,'%%%b2wrsrt: nsp_tmp,nfl_tmp= ',nsp_tmp,nfl_tmp
         print *,'hlp_frm:'
         print *,trim(hlp_frm)
#endif
        endif
        if (nfl_tmp.gt.0) then
          write(lun,hlp_frm) '#','time',
     ,          ((hlp_hr(j),i,i=1,nrg_tmp),j=1,nr_tmp),
     ,          ((hlp_hs(j),species_list(i),i=1,nsp_tmp),j=1,ns_tmp),
     ,         (((hlp_ha(k),j,species_list(i),j=1,nrg_tmp),
     ,                                        i=1,nsp_tmp),k=1,na_tmp)
        end if

      end if

      if (nfl_tmp.gt.0) then
        write(hlp_frm,'(a,i3,a)') '(1p,e16.7,',
     +   nr_tmp*nrg_tmp + ns_tmp*nsp_tmp + na_tmp*nrg_tmp*nsp_tmp,
     +  '(2x,e12.4))'
      end if
#ifdef DBG
      print *,'hlp_frm:'       !###
      print *,trim(hlp_frm)    !###
#endif
      return
c=======================================================================
   90 write(*,*) 'b2wrsrt: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      return
c
   92 write(*,*) 'b2wrsrt: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      return
c=======================================================================
  700 format('#   ',7i8/'#  ',a)
      end subroutine b2wrsrti
c======================================================================
c
      subroutine b2wrsrtc(iline)
      implicit none
      integer, intent(out) :: iline
      real (kind=R8) :: tdummy
      character hlp_c
      logical ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_srt)) return
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
      close(lun)
      return
c======================================================================
  190 write(*,*) 'b2wrsrt: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
      end subroutine b2wrsrtc
c=======================================================================
      subroutine b2wrsrts
c-----------------------------------------------------------------------
c*** Save the values (not used for this set)
c-----------------------------------------------------------------------
c
      return
      end subroutine b2wrsrts
c=======================================================================
      subroutine b2wrsrtf
c-----------------------------------------------------------------------
c*** Finalization step
c-----------------------------------------------------------------------
c
      logical file_open

      lun=ntrc(ltrc)
      inquire(unit=lun,opened=file_open)
      if (file_open) close(lun)
      return
      end subroutine b2wrsrtf
c=======================================================================
      end module b2mod_wrsrt
