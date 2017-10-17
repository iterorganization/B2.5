
      subroutine running_average (nx,ny,ns)

      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_sources
#ifdef B25_EIRENE
      use eirmod_infcop
#endif

      implicit none

      integer, intent(in) :: nx, ny, ns, iav_cont
      logical, intent(in) :: lpr, ler

      real (kind=R8), intent(inout) ::
     *  na_arg(-1:nx,-1:ny,0:ns-1), ua_arg(-1:nx,-1:ny,0:ns-1),
     &  te_arg(-1:nx,-1:ny), ti_arg(-1:nx,-1:ny), po_arg(-1:nx,-1:ny)
      real (kind=R8), intent(inout) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)

      character(*) :: savefile

      real (kind=R8), allocatable, save ::
     .                             na_mean(:,:,:), te_mean(:,:),
     .                             ti_mean(:,:), ua_mean(:,:,:),
     .                             po_mean(:,:)
      real (kind=R8), allocatable, save ::
     .                             e_na(:,:,:), e_te(:,:),
     .                             e_ti(:,:), e_ua(:,:,:),
     .                             e_po(:,:)
      real (kind=R8), allocatable, save ::
     .                             sna_mean(:,:,:,:),
     .                             sne_mean(:,:,:),
     .                             smo_mean(:,:,:,:),
     .                             she_mean(:,:,:),
     .                             shi_mean(:,:,:),
     .                             sch_mean(:,:,:)
      real (kind=R8), allocatable, save ::
     .                             e_sna(:,:,:,:),
     .                             e_sne(:,:,:),
     .                             e_smo(:,:,:,:),
     .                             e_she(:,:,:),
     .                             e_shi(:,:,:),
     .                             e_sch(:,:,:)

      integer :: io, iun, n2, itp, idum(0:2), iex
      integer, save :: naver, ntotdt
      character lble*120, lbli*120, check*10
      character*10, save ::
     *   b2faveri_version = '00.000.000'
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)

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

      call batch_average (n2*2*ns, sna, sna_mean, itp, naver)
      call batch_average (n2*2,    sne, sne_mean, itp, naver)
      call batch_average (n2*4,    she, she_mean, itp, naver)
      call batch_average (n2*4,    shi, shi_mean, itp, naver)
      call batch_average (n2*4,    sch, sch_mean, itp, naver)
      call batch_average (n2*4*ns, smo, smo_mean, itp, naver)

*  calculate squared averaged plasma profiles
      call batch_average_sq (n2*ns, na, e_na, itp, naver)
      call batch_average_sq (n2,    te, e_te, itp, naver)
      call batch_average_sq (n2,    ti, e_ti, itp, naver)
      call batch_average_sq (n2,    po, e_po, itp, naver)
      call batch_average_sq (n2*ns, ua, e_ua, itp, naver)
      
      call batch_average_sq (n2*2*ns, sna, e_sna, itp, naver)
      call batch_average_sq (n2*2,    sne, e_sne, itp, naver)
      call batch_average_sq (n2*4,    she, e_she, itp, naver)
      call batch_average_sq (n2*4,    shi, e_shi, itp, naver)
      call batch_average_sq (n2*4,    sch, e_sch, itp, naver)
      call batch_average_sq (n2*4*ns, smo, e_smo, itp, naver)

      call subend ()

      return


      entry run_av_init (iav_cont, nx, ny, ns, lpr, ler)

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

ctb: the sources are now stored in linearized form I propose to store the total part only
ctb..as the mean of the product does not equal the product of the means etc.
         allocate (sna_mean(-1:nx,-1:ny,0:1,0:ns-1))
         allocate (sne_mean(-1:nx,-1:ny,0:1))
         allocate (smo_mean(-1:nx,-1:ny,0:3,0:ns-1))
         allocate (she_mean(-1:nx,-1:ny,0:3))
         allocate (shi_mean(-1:nx,-1:ny,0:3))
         allocate (sch_mean(-1:nx,-1:ny,0:3))
         allocate (e_sna(-1:nx,-1:ny,0:1,0:ns-1))
         allocate (e_sne(-1:nx,-1:ny,0:1))
         allocate (e_smo(-1:nx,-1:ny,0:3,0:ns-1))
         allocate (e_she(-1:nx,-1:ny,0:3))
         allocate (e_shi(-1:nx,-1:ny,0:3))
         allocate (e_sch(-1:nx,-1:ny,0:3))

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

              call cfrure (iun, n2*2*ns, sna_mean, 'sna_mean')
              call cfrure (iun, n2*2,    sne_mean, 'sne_mean')
              call cfrure (iun, n2*4,    she_mean, 'she_mean')
              call cfrure (iun, n2*4,    shi_mean, 'shi_mean')
              call cfrure (iun, n2*4,    sch_mean, 'sch_mean')
              call cfrure (iun, n2*4*ns, smo_mean, 'smo_mean')
            end if

*  read 2d averaged squared plasma profiles from b2fstati_run_var

            if (ler) then
              call cfrure (iun, n2*ns, e_na, 'e_na')
              call cfrure (iun, n2,    e_te, 'e_te')
              call cfrure (iun, n2,    e_ti, 'e_ti')
              call cfrure (iun, n2,    e_po, 'e_po')
              call cfrure (iun, n2*ns, e_ua, 'e_ua')

              call cfrure (iun, n2*2*ns, e_sna, 'e_sna')
              call cfrure (iun, n2*2,    e_sne, 'e_sne')
              call cfrure (iun, n2*4,    e_she, 'e_she')
              call cfrure (iun, n2*4,    e_shi, 'e_shi')
              call cfrure (iun, n2*4,    e_sch, 'e_sch')
              call cfrure (iun, n2*4*ns, e_smo, 'e_smo')
            end if

         else

            write (*,*) ' NO FILE b2faveri FOUND '
            write (*,*) ' AVERAGING STARTS FROM ZERO '

         end if

         close (iun)

      end if

      call subend ()

      return

      entry run_av_fin

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

      entry run_av_save (savefile, nx, ny, ns, lpr, ler)

      call subini ('run_av_save')

      n2 = (nx+2)*(ny+2)

!      write (*,*) ' in run_av_save '

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
         call cfwure (iun, n2, te_mean, 'te_mean')
         call cfwure (iun, n2, ti_mean, 'ti_mean')
         call cfwure (iun, n2, po_mean, 'po_mean')
         call cfwure (iun, n2*ns, ua_mean, 'ua_mean')

         call cfwure (iun, n2*2*ns, sna_mean, 'sna_mean')
         call cfwure (iun, n2*2, sne_mean, 'sne_mean')
         call cfwure (iun, n2*4, she_mean, 'she_mean')
         call cfwure (iun, n2*4, shi_mean, 'shi_mean')
         call cfwure (iun, n2*4, sch_mean, 'sch_mean')
         call cfwure (iun, n2*4*ns, smo_mean, 'smo_mean')

      end if

*  write 2d squared averaged plasma profiles

      if (ler) then

         call cfwure (iun, n2*ns, e_na, 'e_na')
         call cfwure (iun, n2,    e_te, 'e_te')
         call cfwure (iun, n2,    e_ti, 'e_ti')
         call cfwure (iun, n2,    e_po, 'e_po')
         call cfwure (iun, n2*ns, e_ua, 'e_ua')

         call cfwure (iun, n2*2*ns, e_sna, 'e_sna')
         call cfwure (iun, n2*2,    e_sne, 'e_sne')
         call cfwure (iun, n2*4,    e_she, 'e_she')
         call cfwure (iun, n2*4,    e_shi, 'e_shi')
         call cfwure (iun, n2*4,    e_sch, 'e_sch')
         call cfwure (iun, n2*4*ns, e_smo, 'e_smo')

      end if
!
      close (iun)
      write (*,*) ' file ',savefile,' written '

      call subend ()

      return


      entry run_av_get_plasma (na_arg, te_arg, ti_arg, ua_arg, po_arg,
     &                          nx, ny, ns)

      call subini ('run_av_get_plasma')

      na_arg = na_mean
      ua_arg = ua_mean
      te_arg = te_mean
      ti_arg = ti_mean
      po_arg = po_mean

      call subend ()

      return


      entry run_av_get_sources (sna0, smo0, she0, shi0,
     &                          sch0, sne0, nx, ny, ns)

      call subini ('run_av_get_sources')

      sna0 = sna_mean
      sne0 = sne_mean
      smo0 = smo_mean
      she0 = she_mean
      shi0 = shi_mean
      sch0 = sch_mean

      call subend ()

      end subroutine running_average
