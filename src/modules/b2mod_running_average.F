      module b2mod_running_average

      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version
      use b2mod_constants
      use b2us_plasma
      use b2mod_subsys
      implicit none

      integer, save :: naver, ntotdt, nfluids, nstrata

!! These variables are still needed for conversion from structured to unstructured.
!! The variables used by the unstructured code are saved in B2Averege
      real (kind=R8), allocatable, save ::
     .                             na_mean(:,:,:), te_mean(:,:),
     .                             ti_mean(:,:), ua_mean(:,:,:),
     .                             po_mean(:,:)
      real (kind=R8), allocatable, save ::
     .                             e_na(:,:,:), e_te(:,:),
     .                             e_ti(:,:), e_ua(:,:,:),
     .                             e_po(:,:)
      real (kind=R8), allocatable, save ::
     .                             sna_mean(:,:,:),
     .                             smo_mean(:,:,:),
     .                             she_mean(:,:),
     .                             shi_mean(:,:)
      real (kind=R8), allocatable, save ::
     .                             e_sna(:,:,:),
     .                             e_smo(:,:,:),
     .                             e_she(:,:),
     .                             e_shi(:,:)

      character lble*120, lbli*120, check*10

! For smooth restart with eirene_underrelax
      real (kind=R8), allocatable, save ::
     .          sni_underrelax(:,:,:), smo_underrelax(:,:,:),
     .          sei_underrelax(:,:,:), see_underrelax(:,:)

      contains

      subroutine running_average (nCv, ns, pl, srw, avg)
      implicit none
c input
      integer, intent(in) :: nCv, ns
      type (B2Plasma), intent (in) :: pl
      type (B2SourceWork), intent (in) :: srw
      type (B2Average), intent(inout) :: avg
c local
      integer :: itp
      external batch_average, batch_average_sq

      call subini ('running_average')

      naver = naver + 1
      ntotdt = ntotdt + 1

      itp = 0

*  calculate averaged plasma profiles and sources
      call batch_average (nCv*ns, pl%na, avg%na_mean, itp, naver)
      call batch_average (nCv,    pl%te, avg%te_mean, itp, naver)
      call batch_average (nCv,    pl%ti, avg%ti_mean, itp, naver)
      call batch_average (nCv,    pl%po, avg%po_mean, itp, naver)
      call batch_average (nCv*ns, pl%ua, avg%ua_mean, itp, naver)
      call batch_average (nCv,    pl%kt, avg%kt_mean, itp, naver)
      call batch_average (nCv,    pl%zt, avg%zt_mean, itp, naver)

      call batch_average (nCv*ns, srw%sna0_eir_tot, avg%sna_mean,
     .                                            itp, naver)
      call batch_average (nCv,    srw%she0_eir_tot, avg%she_mean,
     .                                            itp, naver)
      call batch_average (nCv,    srw%shi0_eir_tot, avg%shi_mean,
     .                                            itp, naver)
      call batch_average (nCv,    srw%shn0_eir_tot, avg%shn_mean,
     .                                            itp, naver)
      call batch_average (nCv*ns, srw%smo0_eir_tot, avg%smo_mean,
     .                                            itp, naver)

*  calculate squared averaged plasma profiles
      call batch_average_sq (nCv*ns, pl%na, avg%e_na, itp, naver)
      call batch_average_sq (nCv,    pl%te, avg%e_te, itp, naver)
      call batch_average_sq (nCv,    pl%ti, avg%e_ti, itp, naver)
      call batch_average_sq (nCv,    pl%po, avg%e_po, itp, naver)
      call batch_average_sq (nCv*ns, pl%ua, avg%e_ua, itp, naver)
      call batch_average_sq (nCv,    pl%kt, avg%e_kt, itp, naver)
      call batch_average_sq (nCv,    pl%zt, avg%e_zt, itp, naver)

      call batch_average_sq (nCv*ns, srw%sna0_eir_tot, avg%e_sna,
     .                                            itp, naver)
      call batch_average_sq (nCv,    srw%she0_eir_tot, avg%e_she,
     .                                            itp, naver)
      call batch_average_sq (nCv,    srw%shi0_eir_tot, avg%e_shi,
     .                                            itp, naver)
      call batch_average_sq (nCv,    srw%shn0_eir_tot, avg%e_shn,
     .                                            itp, naver)
      call batch_average_sq (nCv*ns, srw%smo0_eir_tot, avg%e_smo,
     .                                            itp, naver)

      call subend ()
      return
      end subroutine running_average

      subroutine run_av_init (iav_cont, nCv, ns, lpr, ler, avg)
      implicit none
c input
      integer, intent(in) :: nCv, ns, iav_cont
      logical, intent(in) :: lpr, ler
      type(B2Average), intent(inout) :: avg
c local variables
      integer :: io, iun, idum(0:1), iex
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      call subini ('run_av_init')

      if (.not.allocated(avg%na_mean)) then
         allocate (avg%na_mean(nCv,0:ns-1))
         allocate (avg%ua_mean(nCv,0:ns-1))
         allocate (avg%te_mean(nCv))
         allocate (avg%ti_mean(nCv))
         allocate (avg%po_mean(nCv))
         allocate (avg%kt_mean(nCv))
         allocate (avg%zt_mean(nCv))
         allocate (avg%e_na(nCv,0:ns-1))
         allocate (avg%e_ua(nCv,0:ns-1))
         allocate (avg%e_te(nCv))
         allocate (avg%e_ti(nCv))
         allocate (avg%e_po(nCv))
         allocate (avg%e_kt(nCv))
         allocate (avg%e_zt(nCv))

         allocate (avg%sna_mean(nCv,0:ns-1))
         allocate (avg%smo_mean(nCv,0:ns-1))
         allocate (avg%she_mean(nCv))
         allocate (avg%shi_mean(nCv))
         allocate (avg%shn_mean(nCv))
         allocate (avg%e_sna(nCv,0:ns-1))
         allocate (avg%e_smo(nCv,0:ns-1))
         allocate (avg%e_she(nCv))
         allocate (avg%e_shi(nCv))
         allocate (avg%e_shn(nCv))
      end if

      naver = 0
      ntotdt = 0

      avg % na_mean = 0.0_R8
      avg % te_mean = 0.0_R8
      avg % ti_mean = 0.0_R8
      avg % ua_mean = 0.0_R8
      avg % po_mean = 0.0_R8
      avg % kt_mean = 0.0_R8
      avg % zt_mean = 0.0_R8
      avg % e_na = 0.0_R8
      avg % e_te = 0.0_R8
      avg % e_ti = 0.0_R8
      avg % e_ua = 0.0_R8
      avg % e_po = 0.0_R8
      avg % e_kt = 0.0_R8
      avg % e_zt = 0.0_R8

      avg % sna_mean = 0.0_R8
      avg % she_mean = 0.0_R8
      avg % shi_mean = 0.0_R8
      avg % shn_mean = 0.0_R8
      avg % smo_mean = 0.0_R8
      avg % e_sna = 0.0_R8
      avg % e_she = 0.0_R8
      avg % e_shi = 0.0_R8
      avg % e_shn = 0.0_R8
      avg % e_smo = 0.0_R8

      if (iav_cont > 0) then

         naver = 0
         ntotdt = 0
#ifndef LEGACYCOMP
         open (newunit=iun,file='b2faveri',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)
#else
         iun = newunit()
         open (unit=iun,file='b2faveri',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)
#endif
         io = iex
         read (iun,'(A10)',iostat=io) check
*  read 2d averaged plasma profiles from b2fstati_run_av
         if (io == 0) then
            backspace iun
            if (lpr) then
              call cfverr (iun, b2faveri_version)
              call cfruin (iun, 2, idum, 'nCv,ns')
              call cfruch (iun, 120, lbli, 'label')
              call b2ruzd (iun, b2faveri_version, ns,
     &          zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)

              call cfruin (iun, 1, idum, 'naver')
              naver = idum(0)
              call cfruin (iun, 1, idum, 'ntotdt')
              ntotdt = idum(0)
              call cfrure (iun, nCv*ns, avg%na_mean, 'na_mean')
              call cfrure (iun, nCv,    avg%te_mean, 'te_mean')
              call cfrure (iun, nCv,    avg%ti_mean, 'ti_mean')
              call cfrure (iun, nCv,    avg%po_mean, 'po_mean')
              call cfrure (iun, nCv*ns, avg%ua_mean, 'ua_mean')
              call cfrure (iun, nCv,    avg%kt_mean, 'kt_mean')
              call cfrure (iun, nCv,    avg%zt_mean, 'zt_mean')

              call cfrure (iun, nCv*ns, avg%sna_mean, 'sna_mean')
              call cfrure (iun, nCv,    avg%she_mean, 'she_mean')
              call cfrure (iun, nCv,    avg%shi_mean, 'shi_mean')
              call cfrure (iun, nCv,    avg%shn_mean, 'shn_mean')
              call cfrure (iun, nCv*ns, avg%smo_mean, 'smo_mean')
            end if

*  read 2d averaged squared plasma profiles from b2fstati_run_var

            if (ler) then
              call cfrure (iun, nCv*ns, avg%e_na, 'e_na')
              call cfrure (iun, nCv,    avg%e_te, 'e_te')
              call cfrure (iun, nCv,    avg%e_ti, 'e_ti')
              call cfrure (iun, nCv,    avg%e_po, 'e_po')
              call cfrure (iun, nCv*ns, avg%e_ua, 'e_ua')
              call cfrure (iun, nCv,    avg%e_kt, 'e_kt')
              call cfrure (iun, nCv,    avg%e_zt, 'e_zt')

              call cfrure (iun, nCv*ns, avg%e_sna, 'e_sna')
              call cfrure (iun, nCv,    avg%e_she, 'e_she')
              call cfrure (iun, nCv,    avg%e_shi, 'e_shi')
              call cfrure (iun, nCv,    avg%e_shn, 'e_shn')
              call cfrure (iun, nCv*ns, avg%e_smo, 'e_smo')
            end if

         else

            write (*,*) ' NO FILE b2faveri FOUND'
            write (*,*) ' AVERAGING STARTS FROM ZERO'

         end if

         close (iun)

      end if

      call subend ()
      return
      end subroutine run_av_init

      subroutine run_av_fin( avg )
      implicit none
      type(B2Average), intent(inout) :: avg

      if (.not.allocated(avg%na_mean)) then

        deallocate (avg%na_mean)
        deallocate (avg%ua_mean)
        deallocate (avg%te_mean)
        deallocate (avg%ti_mean)
        deallocate (avg%po_mean)
        deallocate (avg%kt_mean)
        deallocate (avg%zt_mean)
        deallocate (avg%e_na)
        deallocate (avg%e_ua)
        deallocate (avg%e_te)
        deallocate (avg%e_ti)
        deallocate (avg%e_po)
        deallocate (avg%e_kt)
        deallocate (avg%e_zt)
        deallocate (avg%sna_mean)
        deallocate (avg%smo_mean)
        deallocate (avg%she_mean)
        deallocate (avg%shi_mean)
        deallocate (avg%shn_mean)
        deallocate (avg%e_sna)
        deallocate (avg%e_smo)
        deallocate (avg%e_she)
        deallocate (avg%e_shi)
        deallocate (avg%e_shn)

      endif

      if (.not.allocated(sni_underrelax)) then

        deallocate (sni_underrelax)
        deallocate (smo_underrelax)
        deallocate (sei_underrelax)
        deallocate (see_underrelax)

      endif

      return
      end subroutine run_av_fin

      subroutine run_av_save (savefile, nCv, ns, lpr, ler, nur, avg)
      implicit none
      integer, intent(in) :: nCv, ns, nur
      logical, intent(in) :: lpr, ler
      type(B2Average), intent(in) :: avg
      character(*) :: savefile
c local variables
      integer :: idum(0:2), iun
      external cfwuch, cfwuin, cfwure
      external b2xxid, b2wuzd
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      call subini ('run_av_save')

*  write 2d averaged plasma profiles

#ifndef LEGACYCOMP
      open (newunit=iun,file=savefile,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')
#else
      iun = newunit()
      open (unit=iun,file=savefile,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')
#endif

!     write (*,*) ' iun ',iun

      if (lpr) then

*   ..get label
         call b2xxid ('b2mn', lble(1:60))
         lble(61:120) = repeat(' ',60)

*   ..write header
         idum(0) = nCv
         idum(1) = ns
         call cfverw (iun, newversion)
         call cfwuin (iun, 2, idum, 'nCv,ns')
         call cfwuch (iun, 120, lble, 'label')
         call b2wuzd (iun, newversion, ns, zamin, zamax, zn, am)

*   ..write average state and sources
         idum(0) = naver
         call cfwuin (iun, 1, idum, 'naver')
         idum(0) = ntotdt
         call cfwuin (iun, 1, idum, 'ntotdt')

         call cfwure (iun, nCv*ns, avg%na_mean, 'na_mean')
         call cfwure (iun, nCv,    avg%te_mean, 'te_mean')
         call cfwure (iun, nCv,    avg%ti_mean, 'ti_mean')
         call cfwure (iun, nCv,    avg%po_mean, 'po_mean')
         call cfwure (iun, nCv*ns, avg%ua_mean, 'ua_mean')
         call cfwure (iun, nCv,    avg%kt_mean, 'kt_mean')
         call cfwure (iun, nCv,    avg%zt_mean, 'zt_mean')

         call cfwure (iun, nCv*ns, avg%sna_mean, 'sna_mean')
         call cfwure (iun, nCv,    avg%she_mean, 'she_mean')
         call cfwure (iun, nCv,    avg%shi_mean, 'shi_mean')
         call cfwure (iun, nCv,    avg%shn_mean, 'shn_mean')
         call cfwure (iun, nCv*ns, avg%smo_mean, 'smo_mean')

      end if

*  write 2d squared averaged plasma profiles

      if (ler) then

         call cfwure (iun, nCv*ns, avg%e_na, 'e_na')
         call cfwure (iun, nCv,    avg%e_te, 'e_te')
         call cfwure (iun, nCv,    avg%e_ti, 'e_ti')
         call cfwure (iun, nCv,    avg%e_po, 'e_po')
         call cfwure (iun, nCv*ns, avg%e_ua, 'e_ua')
         call cfwure (iun, nCv,    avg%e_kt, 'e_kt')
         call cfwure (iun, nCv,    avg%e_zt, 'e_zt')

         call cfwure (iun, nCv*ns, avg%e_sna, 'e_sna')
         call cfwure (iun, nCv,    avg%e_she, 'e_she')
         call cfwure (iun, nCv,    avg%e_shi, 'e_shi')
         call cfwure (iun, nCv,    avg%e_shn, 'e_shn')
         call cfwure (iun, nCv*ns, avg%e_smo, 'e_smo')

      end if

*  write underrelaxed Eirene sources for smooth restart with eirene_underrelaxation

      if (nur.gt.1) then
         idum(0) = nCv
         idum(1) = nfluids
         idum(2) = nstrata
         call cfwuin (iun, 3, idum, 'nCv,nfluids,nstrata')
         call cfwure (iun, nCv*nfluids*nstrata, sni_underrelax,
     .                                          'sni_underrelax')
         call cfwure (iun, nCv*nfluids*nstrata, smo_underrelax,
     .                                          'smo_underrelax')
         call cfwure (iun, nCv*nfluids*nstrata, sei_underrelax,
     .                                          'sei_underrelax')
         call cfwure (iun, nCv*nstrata, see_underrelax,
     .                                          'see_underrelax')
      end if

      close (iun)
      write (*,*) ' file ',savefile,' written '

      call subend ()
      return
      end subroutine run_av_save

      subroutine run_av_get_plasma (na_arg, te_arg, ti_arg, ua_arg,
     &                          po_arg, kt_arg, zt_arg, nCv, ns, avg)
      implicit none
      integer, intent(in) :: nCv, ns
      type(B2Average), intent(in) :: avg
      real (kind=R8), intent(inout) ::
     *  na_arg(nCv,0:ns-1), ua_arg(nCv,0:ns-1),
     &  te_arg(nCv), ti_arg(nCv), po_arg(nCv), kt_arg(nCv), zt_arg(nCv)

      call subini ('run_av_get_plasma')

      na_arg = avg%na_mean
      ua_arg = avg%ua_mean
      te_arg = avg%te_mean
      ti_arg = avg%ti_mean
      po_arg = avg%po_mean
      kt_arg = avg%kt_mean
      zt_arg = avg%zt_mean

      call subend ()
      return
      end subroutine run_av_get_plasma

      subroutine run_av_get_sources (sna_out, smo_out, she_out, shi_out,
     &                               shn_out, nCv, ns, avg)
      implicit none
      integer, intent(in) :: nCv, ns
      type(B2Average), intent(in) :: avg
      real (kind=R8), intent(out) ::
     *  sna_out(nCv,0:ns-1), smo_out(nCv,0:ns-1),
     &  she_out(nCv), shi_out(nCv), shn_out(nCv)

      call subini ('run_av_get_sources')

c only the total sources are averaged, not the linearised components.

      sna_out = avg%sna_mean
      smo_out = avg%smo_mean
      she_out = avg%she_mean
      shi_out = avg%shi_mean
      shn_out = avg%shn_mean

      call subend ()
      return
      end subroutine run_av_get_sources


      subroutine save_underrelaxed_sources (nCv,nfl,nstra,sni_s,smo_s,
     .                                      sei_s,see_s)
      implicit none
c input
      integer, intent(in) :: nCv, nfl, nstra
      real(kind=R8), intent(in) :: sni_s(nCv,nfl,nstra),
     .  smo_s(nCv,nfl,nstra), sei_s(nCv,nfl,nstra), see_s(nCv,nstra)
c     type(B2Average), intent(in) :: avg


      call subini ('save_underrelaxed_sources')

      nfluids = nfl
      nstrata = nstra

      if (.not.allocated(sni_underrelax)) then
         allocate (sni_underrelax(nCv,nfl,nstra))
         allocate (smo_underrelax(nCv,nfl,nstra))
         allocate (sei_underrelax(nCv,nfl,nstra))
         allocate (see_underrelax(nCv,nstra))
      end if

      sni_underrelax = sni_s
      smo_underrelax = smo_s
      sei_underrelax = sei_s
      see_underrelax = see_s

      call subend ()
      return
      end subroutine save_underrelaxed_sources

      subroutine read_b2faveri_st ( nx, ny, ns )
      use b2mod_rates
      implicit none
      integer, intent(in) :: nx, ny, ns
      integer iun, io, iex, is, n2, idum(0:2)
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4), dummy(-1:nx,-1:ny), zero(-1:nx,-1:ny),
     *   swrk2(-1:nx,-1:ny,0:1), swrk2s(-1:nx,-1:ny,0:1,0:ns-1),
     *   swrk4(-1:nx,-1:ny,0:3), swrk4s(-1:nx,-1:ny,0:3,0:ns-1),
     *   ne_mean(-1:nx,-1:ny), ni_mean(-1:nx,-1:ny,0:1)
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      n2 = (nx+2)*(ny+2)

      allocate (na_mean(-1:nx,-1:ny,0:ns-1))
      allocate (ua_mean(-1:nx,-1:ny,0:ns-1))
      allocate (te_mean(-1:nx,-1:ny))
      allocate (ti_mean(-1:nx,-1:ny))
      allocate (po_mean(-1:nx,-1:ny))
      allocate (e_na(-1:nx,-1:ny,0:ns-1))
      allocate (e_ua(-1:nx,-1:ny,0:ns-1))
      allocate (e_te(-1:nx,-1:ny))
      allocate (e_ti(-1:nx,-1:ny))
      allocate (e_po(-1:nx,-1:ny))

      allocate (sna_mean(-1:nx,-1:ny,0:ns-1))
      allocate (smo_mean(-1:nx,-1:ny,0:ns-1))
      allocate (she_mean(-1:nx,-1:ny))
      allocate (shi_mean(-1:nx,-1:ny))
      allocate (e_sna(-1:nx,-1:ny,0:ns-1))
      allocate (e_smo(-1:nx,-1:ny,0:ns-1))
      allocate (e_she(-1:nx,-1:ny))
      allocate (e_shi(-1:nx,-1:ny))

#ifndef LEGACYCOMP
      open (newunit=iun,file='b2faveri',access='SEQUENTIAL',
     .      form='FORMATTED',status='OLD',iostat=iex)
#else
      iun = newunit()
      open (unit=iun,file='b2faveri',access='SEQUENTIAL',
     .     form='FORMATTED',status='OLD',iostat=iex)
#endif
      io = iex
      read (iun,'(A10)',iostat=io) check
*  read 2d averaged plasma profiles from b2faveri
      if (io == 0) then
        backspace iun
        call cfverr (iun, b2faveri_version)
        call cfruin (iun, 3, idum, 'nx,ny,ns')
        call cfruch (iun, 120, lbli, 'label')
        call b2ruzd (iun, b2faveri_version, ns,
     &          zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)

        call cfruin (iun, 1, idum, 'naver')
        naver = idum(0)
        call cfruin (iun, 1, idum, 'ntotdt')
        ntotdt = idum(0)
        call cfrure (iun, n2*ns, na_mean, 'na_mean')
        call cfrure (iun, n2,    te_mean, 'te_mean')
        call cfrure (iun, n2,    ti_mean, 'ti_mean')
        call cfrure (iun, n2,    po_mean, 'po_mean')
        call cfrure (iun, n2*ns, ua_mean, 'ua_mean')

        if (b2faveri_version.ge.'03.000.008') then
          call cfrure (iun, n2*ns, sna_mean, 'sna_mean')
          call cfrure (iun, n2,    dummy,    'sne_mean')
          call cfrure (iun, n2,    she_mean, 'she_mean')
          call cfrure (iun, n2,    shi_mean, 'shi_mean')
          call cfrure (iun, n2,    dummy,    'sch_mean')
          call cfrure (iun, n2*ns, smo_mean, 'smo_mean')
        else
          zero = 0.0_R8
          call b2xpne_st ( nx, ny, ns, rza, na_mean, zero, ne_mean)
          call b2xpni ( n2, ns, na_mean, ni_mean)
          call cfrure (iun, n2*2*ns, swrk2s, 'sna_mean')
          sna_mean(:,:,:) = swrk2s(:,:,0,:) +
     .                      swrk2s(:,:,1,:) * na_mean(:,:,:)
          call cfrure (iun, n2*2,    swrk2 , 'sne_mean')
          call cfrure (iun, n2*4,    swrk4 , 'she_mean')
          she_mean(:,:) = swrk4(:,:,0) +
     .                    swrk4(:,:,1) * ne_mean(:,:) +
     .                    swrk4(:,:,2) * te_mean(:,:) +
     .                    swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
          call cfrure (iun, n2*4,    swrk4 , 'shi_mean')
          shi_mean(:,:) = swrk4(:,:,0) +
     .                    swrk4(:,:,1) * ni_mean(:,:,0) +
     .                    swrk4(:,:,2) * ti_mean(:,:) +
     .                    swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
          call cfrure (iun, n2*4,    swrk4 , 'sch_mean')
          call cfrure (iun, n2*4*ns, swrk4s, 'smo_mean')
          do is = 0, ns-1
            smo_mean(:,:,is) = swrk4s(:,:,0,is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) +
     .             swrk4s(:,:,2,is) * na_mean(:,:,is) * mp * am(is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) *
     .                                na_mean(:,:,is) * mp * am(is)
          end do
        end if

        read (iun,'(A10)',iostat=io) check
*  read 2d averaged squared plasma profiles from b2faveri if present
        if (io == 0) then
          backspace iun
          call cfrure (iun, n2*ns, e_na, 'e_na')
          call cfrure (iun, n2,    e_te, 'e_te')
          call cfrure (iun, n2,    e_ti, 'e_ti')
          call cfrure (iun, n2,    e_po, 'e_po')
          call cfrure (iun, n2*ns, e_ua, 'e_ua')

          if (b2faveri_version.ge.'03.000.008') then
            call cfrure (iun, n2*ns, e_sna, 'e_sna')
            call cfrure (iun, n2,    dummy, 'e_sne')
            call cfrure (iun, n2,    e_she, 'e_she')
            call cfrure (iun, n2,    e_shi, 'e_shi')
            call cfrure (iun, n2,    dummy, 'e_sch')
            call cfrure (iun, n2*ns, e_smo, 'e_smo')
          else
            call cfrure (iun, n2*2*ns, swrk2s, 'e_sna')
            e_sna(:,:,:) = swrk2s(:,:,0,:) +
     .                     swrk2s(:,:,1,:) * na_mean(:,:,:)
            call cfrure (iun, n2*2,    swrk2 , 'e_sne')
            call cfrure (iun, n2*4,    swrk4 , 'e_she')
            e_she(:,:) = swrk4(:,:,0) +
     .                   swrk4(:,:,1) * ne_mean(:,:) +
     .                   swrk4(:,:,2) * te_mean(:,:) +
     .                   swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
            call cfrure (iun, n2*4,    swrk4 , 'e_shi')
            e_shi(:,:) = swrk4(:,:,0) +
     .                   swrk4(:,:,1) * ni_mean(:,:,0) +
     .                   swrk4(:,:,2) * ti_mean(:,:) +
     .                   swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
            call cfrure (iun, n2*4,    swrk4 , 'e_sch')
            call cfrure (iun, n2*4*ns, swrk4s, 'e_smo')
            do is = 0, ns-1
              e_smo(:,:,is) = swrk4s(:,:,0,is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) +
     .             swrk4s(:,:,2,is) * na_mean(:,:,is) * mp * am(is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) *
     .                                na_mean(:,:,is) * mp * am(is)
            end do
          end if
        end if

      end if
      close (iun)

      return
      end subroutine read_b2faveri_st

      subroutine dealloc_avg_st
      implicit none

      deallocate (na_mean)
      deallocate (ua_mean)
      deallocate (te_mean)
      deallocate (ti_mean)
      deallocate (po_mean)
      deallocate (e_na)
      deallocate (e_ua)
      deallocate (e_te)
      deallocate (e_ti)
      deallocate (e_po)
      deallocate (sna_mean)
      deallocate (smo_mean)
      deallocate (she_mean)
      deallocate (shi_mean)
      deallocate (e_sna)
      deallocate (e_smo)
      deallocate (e_she)
      deallocate (e_shi)

      return
      end subroutine dealloc_avg_st

      subroutine read_underrelaxed_sources(nCv, nfl, nstra)
      implicit none
c input
      integer, intent(in) :: nCv, nfl, nstra
c nfl is the number of species indices for the EIRENE sources sni, smo, and sei, =/= ns in B2.5
c local variables
      integer :: io, iun, idum(0:2), iex
      character :: check*10, line*51
      character lbli*120
      character*10 :: version
      logical :: found
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

      call subini ('read_underrelaxed_sources')

      if(.not.allocated(sni_underrelax)) then
          allocate (sni_underrelax(nCv,nfl,nstra))
          allocate (smo_underrelax(nCv,nfl,nstra))
          allocate (sei_underrelax(nCv,nfl,nstra))
          allocate (see_underrelax(nCv,nstra))
      end if
      sni_underrelax = 0.0_R8
      smo_underrelax = 0.0_R8
      sei_underrelax = 0.0_R8
      see_underrelax = 0.0_R8

#ifndef LEGACYCOMP
      open (newunit=iun,file='b2faveri',
     .    access='SEQUENTIAL',form='FORMATTED',status='OLD',iostat=iex)
#else
         iun = newunit()
      open (unit=iun,file='b2faveri',
     .    access='SEQUENTIAL',form='FORMATTED',status='OLD',iostat=iex)
#endif

*  read saved underrelaxed sources for smooth restart with nrelax.gt.1
      io = iex
      read (iun,'(A10)',iostat=io) check
      if (io == 0) then
        ! Scan until we find the start of underrelaxed sources data
        do
          read(iun, '(A)', iostat=io) line
          if (io /= 0) then
            call xerrab('Underrelaxed sources not found before EOF')
          end if
          if (trim(line) ==
     .    '*cf:    int                3    nCv,nfluids,nstrata') then
            found = .true.
            exit
          end if
        end do
        backspace iun
        call cfruin (iun, 3, idum, 'nCv,nfluids,nstrata')
        if (idum(0).ne.nCv) call xerrab('nCv not consistent')
        if (idum(1).ne.nfl) call xerrab('nfl not consistent')
        if (idum(2).ne.nstra) call xerrab('nstra not consistent')
        call cfrure (iun, nCv*nfl*nstra, sni_underrelax,
     .                                  'sni_underrelax')
        call cfrure (iun, nCv*nfl*nstra, smo_underrelax,
     .                                  'smo_underrelax')
        call cfrure (iun, nCv*nfl*nstra, sei_underrelax,
     .                                  'sei_underrelax')
        call cfrure (iun, nCv*nstra,    see_underrelax,
     .                                 'see_underrelax')
        write (*,*) ' eirene_underrelaxed_sources read succesfully '
      else
        call xerrab('NO FILE b2faveri FOUND')
      end if
      close (iun)

      call subend ()
      return
      end subroutine read_underrelaxed_sources

      end module b2mod_running_average
