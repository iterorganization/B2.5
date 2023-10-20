      module b2mod_batch
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version
      implicit none

      integer, save :: naver
!     integer, save :: ntotdt

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
     .                             sne_mean(:,:),
     .                             smo_mean(:,:,:),
     .                             she_mean(:,:),
     .                             shi_mean(:,:),
     .                             sch_mean(:,:)
      real (kind=R8), allocatable, save ::
     .                             e_sna(:,:,:),
     .                             e_sne(:,:),
     .                             e_smo(:,:,:),
     .                             e_she(:,:),
     .                             e_shi(:,:),
     .                             e_sch(:,:)

      contains

      subroutine batch_av_all (itim,ntim_batch,nx,ny,ns)
      use b2mod_constants
      use b2mod_plasma
      use b2mod_sources
      implicit none

      integer, intent(in) :: itim, ntim_batch, nx, ny, ns

      integer :: n2, itp

      n2 = (nx+2)*(ny+2)

      itp = itim + 1

*  calculate averaged plasma profiles and sources
      call batch_average (n2*ns, na, na_mean, itp, ntim_batch)
      call batch_average (n2,    te, te_mean, itp, ntim_batch)
      call batch_average (n2,    ti, ti_mean, itp, ntim_batch)
      call batch_average (n2,    po, po_mean, itp, ntim_batch)
      call batch_average (n2*ns, ua, ua_mean, itp, ntim_batch)

      call batch_average (n2*ns, sna, sna_mean, itp, ntim_batch)
      call batch_average (n2,    sne, sne_mean, itp, ntim_batch)
      call batch_average (n2,    she, she_mean, itp, ntim_batch)
      call batch_average (n2,    shi, shi_mean, itp, ntim_batch)
      call batch_average (n2,    sch, sch_mean, itp, ntim_batch)
      call batch_average (n2*ns, smo, smo_mean, itp, ntim_batch)

*  calculate squared averaged plasma profiles
      call batch_average_sq (n2*ns, na, e_na, itp, ntim_batch)
      call batch_average_sq (n2,    te, e_te, itp, ntim_batch)
      call batch_average_sq (n2,    ti, e_ti, itp, ntim_batch)
      call batch_average_sq (n2,    po, e_po, itp, ntim_batch)
      call batch_average_sq (n2*ns, ua, e_ua, itp, ntim_batch)

      call batch_average_sq (n2*ns, sna, e_sna, itp, ntim_batch)
      call batch_average_sq (n2,    sne, e_sne, itp, ntim_batch)
      call batch_average_sq (n2,    she, e_she, itp, ntim_batch)
      call batch_average_sq (n2,    shi, e_shi, itp, ntim_batch)
      call batch_average_sq (n2,    sch, e_sch, itp, ntim_batch)
      call batch_average_sq (n2*ns, smo, e_smo, itp, ntim_batch)

      naver = naver + 1

      return
      end subroutine batch_av_all

      subroutine batch_av_all_init(nx,ny,ns)
      implicit none
      integer, intent(in) :: nx, ny, ns

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

      naver = 0

      return
      end subroutine batch_av_all_init

      subroutine batch_av_all_fin
      implicit none

      if (.not.allocated(na_mean)) return

      deallocate (na_mean)
      deallocate (ua_mean)
      deallocate (te_mean)
      deallocate (ti_mean)
      deallocate (po_mean)
      deallocate (sna_mean)
      deallocate (sne_mean)
      deallocate (smo_mean)
      deallocate (she_mean)
      deallocate (shi_mean)
      deallocate (sch_mean)
      deallocate (e_na)
      deallocate (e_ua)
      deallocate (e_te)
      deallocate (e_ti)
      deallocate (e_po)
      deallocate (e_sna)
      deallocate (e_sne)
      deallocate (e_smo)
      deallocate (e_she)
      deallocate (e_shi)
      deallocate (e_sch)

      return
      end subroutine batch_av_all_fin

      subroutine batch_av_all_save (batch_name,nx,ny,ns)
      implicit none
      integer, intent(in) :: nx, ny, ns
      character(*), intent(in) :: batch_name
      integer iun, n2, idum(0:2)
      character lble*120
#ifndef F2003
      integer newunit
      external newunit
#endif

      n2 = (nx+2)*(ny+2)
#ifdef F2003
      open (newunit=iun,file=batch_name,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')
#else
      iun = newunit()
      open (unit=iun,file=batch_name,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')
#endif

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

      idum(0) = naver
      call cfwuin (iun, 1, idum, 'naver')
!     idum(0) = ntotdt
!     call cfwuin (iun, 1, idum, 'ntotdt')

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

      naver = 0

      close(iun)
      return
      end subroutine batch_av_all_save

      end module b2mod_batch

