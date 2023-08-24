
      module b2mod_batch_average
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_sources
      implicit none
      real (kind=R8), allocatable, save ::
     .                             na_mean(:,:), te_mean(:),
     .                             ti_mean(:), ua_mean(:,:),
     .                             po_mean(:), kt_mean(:)
      real (kind=R8), allocatable, save :: 
     .                             e_na(:,:), e_te(:),
     .                             e_ti(:), e_ua(:,:),
     .                             e_po(:), e_kt(:)
      real (kind=R8), allocatable, save ::
     .                             sna_mean(:,:),
     .                             smo_mean(:,:),
     .                             she_mean(:),
     .                             shi_mean(:)
      real (kind=R8), allocatable, save ::
     .                             e_sna(:,:),
     .                             e_smo(:,:),
     .                             e_she(:),
     .                             e_shi(:)
      integer, save :: naver, ntotdt
      contains

      subroutine batch_av_all (itim,ntim_batch,nCv,ns)
      implicit none
      integer, intent(in) :: itim, ntim_batch, nCv, ns
      integer :: itp
      real (kind=R8) ::
     *   sna_tot(nCv,0:ns-1), smo_tot(nCv,0:ns-1), 
     *   she_tot(nCv), shi_tot(nCv)

      itp = itim + 1

*  calculate averaged plasma profiles and sources
      call batch_average (nCv*ns, na, na_mean, itp, ntim_batch)
      call batch_average (nCv,    te, te_mean, itp, ntim_batch)
      call batch_average (nCv,    ti, ti_mean, itp, ntim_batch)
      call batch_average (nCv,    po, po_mean, itp, ntim_batch)
      call batch_average (nCv*ns, ua, ua_mean, itp, ntim_batch)
      call batch_average (nCv,    kt, kt_mean, itp, ntim_batch)

      call batch_average (nCv*ns, sna_tot, sna_mean, itp, ntim_batch)
      call batch_average (nCv,    she_tot, she_mean, itp, ntim_batch)
      call batch_average (nCv,    shi_tot, shi_mean, itp, ntim_batch)
      call batch_average (nCv*ns, smo_tot, smo_mean, itp, ntim_batch)

*  calculate squared averaged plasma profiles
      call batch_average_sq (nCv*ns, na, e_na, itp, ntim_batch)
      call batch_average_sq (nCv,    te, e_te, itp, ntim_batch)
      call batch_average_sq (nCv,    ti, e_ti, itp, ntim_batch)
      call batch_average_sq (nCv,    po, e_po, itp, ntim_batch)
      call batch_average_sq (nCv*ns, ua, e_ua, itp, ntim_batch)
      call batch_average_sq (nCv,    kt, e_kt, itp, ntim_batch)
      
      call batch_average_sq (nCv*ns, sna_tot, e_sna, itp, ntim_batch)
      call batch_average_sq (nCv,    she_tot, e_she, itp, ntim_batch)
      call batch_average_sq (nCv,    shi_tot, e_shi, itp, ntim_batch)
      call batch_average_sq (nCv*ns, smo_tot, e_smo, itp, ntim_batch)

      naver = naver + 1

      return

      end subroutine batch_av_all

      subroutine batch_av_all_init(nCv,ns)
      implicit none
      integer, intent(in) :: nCv, ns

      if (.not.allocated(na_mean)) then

         allocate (na_mean(1:nCv,0:ns-1))
         allocate (ua_mean(1:nCv,0:ns-1))      
         allocate (te_mean(1:nCv))        
         allocate (ti_mean(1:nCv))
         allocate (po_mean(1:nCv))
         allocate (kt_mean(1:nCv))
         allocate (e_na(1:nCv,0:ns-1))
         allocate (e_ua(1:nCv,0:ns-1))
         allocate (e_te(1:nCv))
         allocate (e_ti(1:nCv))
         allocate (e_po(1:nCv))
         allocate (e_kt(1:nCv))

ctb: the sources are now stored in linearized form I propose to store the total part only 
ctb..as the mean of the product does not equal the product of the means etc.
         allocate (sna_mean(nCv,0:ns-1))
         allocate (smo_mean(nCv,0:ns-1))
         allocate (she_mean(nCv))
         allocate (shi_mean(nCv))
         allocate (e_sna(nCv,0:ns-1))
         allocate (e_smo(nCv,0:ns-1))
         allocate (e_she(nCv))
         allocate (e_shi(nCv))

      end if

      na_mean = 0.0_R8
      te_mean = 0.0_R8
      ti_mean = 0.0_R8
      ua_mean = 0.0_R8
      po_mean = 0.0_R8
      kt_mean = 0.0_R8
      e_na = 0.0_R8
      e_te = 0.0_R8
      e_ti = 0.0_R8
      e_ua = 0.0_R8
      e_po = 0.0_R8
      e_kt = 0.0_R8
      
      sna_mean = 0.0_R8
      she_mean = 0.0_R8
      shi_mean = 0.0_R8
      smo_mean = 0.0_R8
      e_sna = 0.0_R8
      e_she = 0.0_R8
      e_shi = 0.0_R8
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
      deallocate (kt_mean)
      deallocate (e_na)
      deallocate (e_ua)
      deallocate (e_te)
      deallocate (e_ti)
      deallocate (e_po)
      deallocate (e_kt)
      deallocate (sna_mean)
      deallocate (smo_mean)
      deallocate (she_mean)
      deallocate (shi_mean)
      deallocate (e_sna)
      deallocate (e_smo)
      deallocate (e_she)
      deallocate (e_shi)

      return
      end subroutine batch_av_all_fin

      subroutine batch_av_all_save (batch_name,nCv,ns)
      implicit none
      integer, intent(in) :: nCv, ns
      character(*), intent(in) :: batch_name
      integer :: idum(0:1), iun
      character lble*120
#ifndef F2003
      integer newunit
      external newunit
#endif

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
      idum(0) = nCv
      idum(1) = ns
      call cfverw (iun, newversion)
      call cfwuin (iun, 2, idum, 'nCv,ns')
      call cfwuch (iun, 120, lble, 'label')
      call b2wuzd (iun, newversion, ns, zamin, zamax, zn, am)

      idum(0) = naver
      call cfwuin (iun, 1, idum, 'naver')
      idum(0) = ntotdt
!     call cfwuin (iun, 1, idum, 'ntotdt')
      
      call cfwure (iun, nCv*ns, na_mean, 'na_mean')
      call cfwure (iun, nCv, te_mean, 'te_mean')
      call cfwure (iun, nCv, ti_mean, 'ti_mean')
      call cfwure (iun, nCv, po_mean, 'po_mean')
      call cfwure (iun, nCv*ns, ua_mean, 'ua_mean')
      if (newversion.ge.'03.002.000')
     &  call cfwure (iun, nCv, kt_mean, 'kt_mean')
      
      call cfwure (iun, nCv*ns, sna_mean, 'sna_mean')
      call cfwure (iun, nCv, she_mean, 'she_mean')
      call cfwure (iun, nCv, shi_mean, 'shi_mean')
      call cfwure (iun, nCv*ns, smo_mean, 'smo_mean')
      
      call cfwure (iun, nCv*ns, e_na, 'e_na')
      call cfwure (iun, nCv,    e_te, 'e_te')
      call cfwure (iun, nCv,    e_ti, 'e_ti')
      call cfwure (iun, nCv,    e_po, 'e_po')
      call cfwure (iun, nCv*ns, e_ua, 'e_ua')
      if (newversion.ge.'03.002.000')
     &  call cfwure (iun, nCv,    e_kt, 'e_kt')
      
      call cfwure (iun, nCv*ns, e_sna, 'e_sna')
      call cfwure (iun, nCv,    e_she, 'e_she')
      call cfwure (iun, nCv,    e_shi, 'e_shi')
      call cfwure (iun, nCv*ns, e_smo, 'e_smo')

      naver = 0
      
      close(iun)
      return
      
      end subroutine batch_av_all_save
      
      end module b2mod_batch_average
