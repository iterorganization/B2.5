      module b2mod_running_average

      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version
      use b2mod_constants
      use b2mod_subsys
      implicit none

      integer, save :: naver, ntotdt
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
     .                             shi_mean(:,:),
     .                             sne_mean(:,:),
     .                             sch_mean(:,:)
      real (kind=R8), allocatable, save ::
     .                             e_sna(:,:,:),
     .                             e_smo(:,:,:),
     .                             e_she(:,:),
     .                             e_shi(:,:),
     .                             e_sne(:,:),
     .                             e_sch(:,:)

      character lble*120

      contains

      subroutine running_average (nx,ny,ns)
      use b2mod_version
      use b2mod_constants
      use b2mod_plasma
      use b2mod_sources
      implicit none
c input
      integer, intent(in) :: nx, ny, ns
c local
      integer :: n2, itp
      external batch_average, batch_average_sq

      call subini ('running_average')

      n2 = (nx+2)*(ny+2)

      naver = naver + 1
      ntotdt = ntotdt + 1

      itp = 0

*  calculate averaged plasma profiles and sources
      call batch_average (n2*ns, na, na_mean, itp, naver)
      call batch_average (n2,    te, te_mean, itp, naver)
      call batch_average (n2,    ti, ti_mean, itp, naver)
      call batch_average (n2,    po, po_mean, itp, naver)
      call batch_average (n2*ns, ua, ua_mean, itp, naver)

      call batch_average (n2*ns, sna0_eir_tot, sna_mean, itp, naver)
      call batch_average (n2,    sne0_eir_tot, sne_mean, itp, naver)
      call batch_average (n2,    she0_eir_tot, she_mean, itp, naver)
      call batch_average (n2,    shi0_eir_tot, shi_mean, itp, naver)
      call batch_average (n2,    sch0_eir_tot, sch_mean, itp, naver)
      call batch_average (n2*ns, smo0_eir_tot, smo_mean, itp, naver)

*  calculate squared averaged plasma profiles
      call batch_average_sq (n2*ns, na, e_na, itp, naver)
      call batch_average_sq (n2,    te, e_te, itp, naver)
      call batch_average_sq (n2,    ti, e_ti, itp, naver)
      call batch_average_sq (n2,    po, e_po, itp, naver)
      call batch_average_sq (n2*ns, ua, e_ua, itp, naver)

      call batch_average_sq (n2*ns, sna0_eir_tot, e_sna, itp, naver)
      call batch_average_sq (n2,    sne0_eir_tot, e_sne, itp, naver)
      call batch_average_sq (n2,    she0_eir_tot, e_she, itp, naver)
      call batch_average_sq (n2,    shi0_eir_tot, e_shi, itp, naver)
      call batch_average_sq (n2,    sch0_eir_tot, e_sch, itp, naver)
      call batch_average_sq (n2*ns, smo0_eir_tot, e_smo, itp, naver)

      call subend ()

      return
      end subroutine running_average

      subroutine run_av_init (iav_cont, nx, ny, ns, lpr, ler)
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_infcop    ! IGNORE
#else
      use b2_4to5          ! IGNORE
#endif
#endif
      use b2mod_rates
      implicit none
c input
      integer, intent(in) :: nx, ny, ns, iav_cont
      logical, intent(in) :: lpr, ler
c local variables
      integer :: io, is, iex, iun, n2, idum(0:2)
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4), zero(-1:nx,-1:ny),
     *   swrk2(-1:nx,-1:ny,0:1), swrk2s(-1:nx,-1:ny,0:1,0:ns-1),
     *   swrk4(-1:nx,-1:ny,0:3), swrk4s(-1:nx,-1:ny,0:3,0:ns-1),
     *   ne_mean(-1:nx,-1:ny), ni_mean(-1:nx,-1:ny,0:1)
      character check*10
      character lbli*120
      external b2ruzd, b2xpne, b2xpni
      external cfruch, cfruin, cfrure

      call subini ('run_av_init')

      if (.not.allocated(na_mean)) then
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
         allocate (sne_mean(-1:nx,-1:ny))
         allocate (smo_mean(-1:nx,-1:ny,0:ns-1))
         allocate (she_mean(-1:nx,-1:ny))
         allocate (shi_mean(-1:nx,-1:ny))
         allocate (sch_mean(-1:nx,-1:ny))
         allocate (e_sna(-1:nx,-1:ny,0:ns-1))
         allocate (e_sne(-1:nx,-1:ny))
         allocate (e_smo(-1:nx,-1:ny,0:ns-1))
         allocate (e_she(-1:nx,-1:ny))
         allocate (e_shi(-1:nx,-1:ny))
         allocate (e_sch(-1:nx,-1:ny))

      end if

      n2 = (nx+2)*(ny+2)

      naver = 0
      ntotdt = 0

      na_mean = 0.0_R8
      te_mean = 0.0_R8
      ti_mean = 0.0_R8
      ua_mean = 0.0_R8
      po_mean = 0.0_R8
      e_na = 0.0_R8
      e_te = 0.0_R8
      e_ti = 0.0_R8
      e_ua = 0.0_R8
      e_po = 0.0_R8

      sna_mean = 0.0_R8
      sne_mean = 0.0_R8
      she_mean = 0.0_R8
      shi_mean = 0.0_R8
      sch_mean = 0.0_R8
      smo_mean = 0.0_R8
      e_sna = 0.0_R8
      e_sne = 0.0_R8
      e_she = 0.0_R8
      e_shi = 0.0_R8
      e_sch = 0.0_R8
      e_smo = 0.0_R8

      if (iav_cont > 0) then

         naver = 0
         ntotdt = 0

#ifdef B25_EIRENE
         iun = newunit()
#else
         iun = 99
#endif
         open (unit=iun,file='b2faveri',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=iex)

         io = iex
         read (iun,'(A10)',iostat=io) check
*  read 2d averaged plasma profiles from b2fstati_run_av
         if (io == 0) then
            backspace iun
            if (lpr) then
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
                call cfrure (iun, n2,    sne_mean, 'sne_mean')
                call cfrure (iun, n2,    she_mean, 'she_mean')
                call cfrure (iun, n2,    shi_mean, 'shi_mean')
                call cfrure (iun, n2,    sch_mean, 'sch_mean')
                call cfrure (iun, n2*ns, smo_mean, 'smo_mean')
              else
                zero = 0.0_R8
                call b2xpne ( nx, ny, ns, rza, na_mean, zero, ne_mean)
                call b2xpni ( nx, ny, ns, na_mean, ni_mean)
                call cfrure (iun, n2*2*ns, swrk2s, 'sna_mean')
                sna_mean(:,:,:) = swrk2s(:,:,0,:) +
     .           swrk2s(:,:,1,:) * na_mean(:,:,:)
                call cfrure (iun, n2*2,    swrk2 , 'sne_mean')
                sne_mean(:,:) = swrk2(:,:,0) +
     .           swrk2(:,:,1) * ne_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'she_mean')
                she_mean(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * ne_mean(:,:) +
     .           swrk4(:,:,2) * te_mean(:,:) +
     .           swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'shi_mean')
                shi_mean(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * ni_mean(:,:,0) +
     .           swrk4(:,:,2) * ti_mean(:,:) +
     .           swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'sch_mean')
                sch_mean(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * po_mean(:,:) +
     .           swrk4(:,:,2) * ne_mean(:,:) +
     .           swrk4(:,:,3) * ne_mean(:,:) * po_mean(:,:)
                call cfrure (iun, n2*4*ns, swrk4s, 'smo_mean')
                do is = 0, ns-1
                  smo_mean(:,:,is) = swrk4s(:,:,0,is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) +
     .             swrk4s(:,:,2,is) * na_mean(:,:,is) * mp * am(is) +
     .             swrk4s(:,:,1,is) * ua_mean(:,:,is) *
     .                                na_mean(:,:,is) * mp * am(is)
                end do
              end if
            end if

*  read 2d averaged squared plasma profiles from b2fstati_run_var

            if (ler) then
              call cfrure (iun, n2*ns, e_na, 'e_na')
              call cfrure (iun, n2,    e_te, 'e_te')
              call cfrure (iun, n2,    e_ti, 'e_ti')
              call cfrure (iun, n2,    e_po, 'e_po')
              call cfrure (iun, n2*ns, e_ua, 'e_ua')

              if (b2faveri_version.ge.'03.000.008') then
                call cfrure (iun, n2*ns, e_sna, 'e_sna')
                call cfrure (iun, n2,    e_sne, 'e_sne')
                call cfrure (iun, n2,    e_she, 'e_she')
                call cfrure (iun, n2,    e_shi, 'e_shi')
                call cfrure (iun, n2,    e_sch, 'e_sch')
                call cfrure (iun, n2*ns, e_smo, 'e_smo')
              else
               call cfrure (iun, n2*2*ns, swrk2s, 'e_sna')
                e_sna(:,:,:) = swrk2s(:,:,0,:) +
     .           swrk2s(:,:,1,:) * na_mean(:,:,:)
                call cfrure (iun, n2*2,    swrk2 , 'e_sne')
                e_sne(:,:) = swrk2(:,:,0) +
     .           swrk2(:,:,1) * ne_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'e_she')
                e_she(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * ne_mean(:,:) +
     .           swrk4(:,:,2) * te_mean(:,:) +
     .           swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'e_shi')
                e_shi(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * ni_mean(:,:,0) +
     .           swrk4(:,:,2) * ti_mean(:,:) +
     .           swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
                call cfrure (iun, n2*4,    swrk4 , 'e_sch')
                e_sch(:,:) = swrk4(:,:,0) +
     .           swrk4(:,:,1) * po_mean(:,:) +
     .           swrk4(:,:,2) * ne_mean(:,:) +
     .           swrk4(:,:,3) * ne_mean(:,:) * po_mean(:,:)
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

         else

            write (*,*) ' NO FILE b2faveri FOUND '
            write (*,*) ' AVERAGING STARTS FROM ZERO '

         end if

         close (iun)

      end if

      call subend ()

      return
      end subroutine run_av_init

      subroutine run_av_fin
      implicit none

      if (.not.allocated(na_mean)) return

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
      deallocate (sne_mean)
      deallocate (smo_mean)
      deallocate (she_mean)
      deallocate (shi_mean)
      deallocate (sch_mean)
      deallocate (e_sna)
      deallocate (e_sne)
      deallocate (e_smo)
      deallocate (e_she)
      deallocate (e_shi)
      deallocate (e_sch)

      return
      end subroutine run_av_fin

      subroutine run_av_save (savefile, nx, ny, ns, lpr, ler)
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_infcop    ! IGNORE
#else
      use b2_4to5          ! IGNORE
#endif
#endif
      implicit none
      integer, intent(in) :: nx, ny, ns
      logical, intent(in) :: lpr, ler
      character(*) :: savefile
c local variables
      integer :: n2, idum(0:2), iun
      external cfwuch, cfwuin, cfwure
      external b2xxid, b2wuzd

      call subini ('run_av_save')

      n2 = (nx+2)*(ny+2)

*  write 2d averaged plasma profiles

#ifdef B25_EIRENE
      iun = newunit()
#else
      iun = 99
#endif
      open (unit=iun,file=savefile,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')

!     write (*,*) ' iun ',iun

      if (lpr) then

*   ..get label
         call b2xxid ('b2mn', lble(1:60))
         lble(61:120) =
     *  '                                                            '

*   ..write header
         idum(0) = nx
         idum(1) = ny
         idum(2) = ns
         call cfverw (iun, newversion)
         call cfwuin (iun, 3, idum, 'nx,ny,ns')
         call cfwuch (iun, 120, lble, 'label')
         call b2wuzd (iun, newversion, ns, zamin, zamax, zn, am)

*   ..write average state and sources
         idum(0) = naver
         call cfwuin (iun, 1, idum, 'naver')
         idum(0) = ntotdt
         call cfwuin (iun, 1, idum, 'ntotdt')

         call cfwure (iun, n2*ns, na_mean, 'na_mean')
         call cfwure (iun, n2,    te_mean, 'te_mean')
         call cfwure (iun, n2,    ti_mean, 'ti_mean')
         call cfwure (iun, n2,    po_mean, 'po_mean')
         call cfwure (iun, n2*ns, ua_mean, 'ua_mean')

         call cfwure (iun, n2*ns, sna_mean, 'sna_mean')
         call cfwure (iun, n2,    sne_mean, 'sne_mean')
         call cfwure (iun, n2,    she_mean, 'she_mean')
         call cfwure (iun, n2,    shi_mean, 'shi_mean')
         call cfwure (iun, n2,    sch_mean, 'sch_mean')
         call cfwure (iun, n2*ns, smo_mean, 'smo_mean')

      end if

*  write 2d squared averaged plasma profiles

      if (ler) then

         call cfwure (iun, n2*ns, e_na, 'e_na')
         call cfwure (iun, n2,    e_te, 'e_te')
         call cfwure (iun, n2,    e_ti, 'e_ti')
         call cfwure (iun, n2,    e_po, 'e_po')
         call cfwure (iun, n2*ns, e_ua, 'e_ua')

         call cfwure (iun, n2*ns, e_sna, 'e_sna')
         call cfwure (iun, n2,    e_sne, 'e_sne')
         call cfwure (iun, n2,    e_she, 'e_she')
         call cfwure (iun, n2,    e_shi, 'e_shi')
         call cfwure (iun, n2,    e_sch, 'e_sch')
         call cfwure (iun, n2*ns, e_smo, 'e_smo')

      end if
!
      close (iun)
      write (*,*) ' file ',savefile,' written '

      call subend ()

      return
      end subroutine run_av_save

      subroutine run_av_get_plasma (na_arg, te_arg, ti_arg, ua_arg,
     &                          po_arg, nx, ny, ns)
      implicit none
      integer, intent(in) :: nx, ny, ns
      real (kind=R8), intent(inout) ::
     *  na_arg(-1:nx,-1:ny,0:ns-1), ua_arg(-1:nx,-1:ny,0:ns-1),
     &  te_arg(-1:nx,-1:ny), ti_arg(-1:nx,-1:ny), po_arg(-1:nx,-1:ny)

      call subini ('run_av_get_plasma')

      na_arg = na_mean
      ua_arg = ua_mean
      te_arg = te_mean
      ti_arg = ti_mean
      po_arg = po_mean

      call subend ()

      return
      end subroutine run_av_get_plasma

      subroutine run_av_get_sources (sna0, smo0, she0, shi0,
     &                               sch0, sne0, nx, ny, ns)
      implicit none
      integer, intent(in) :: nx, ny, ns
      real (kind=R8), intent(inout) ::
     *  sna0(-1:nx,-1:ny,0:ns-1), smo0(-1:nx,-1:ny,0:ns-1),
     &  she0(-1:nx,-1:ny), shi0(-1:nx,-1:ny),
     &  sch0(-1:nx,-1:ny), sne0(-1:nx,-1:ny)

      call subini ('run_av_get_sources')

      sna0 = sna_mean
      sne0 = sne_mean
      smo0 = smo_mean
      she0 = she_mean
      shi0 = shi_mean
      sch0 = sch_mean

      call subend ()

      return
      end subroutine run_av_get_sources

      end module b2mod_running_average
