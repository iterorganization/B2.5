      module b2mod_running_average

      use b2mod_types
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_plasma
#ifdef B25_EIRENE
      use eirmod_infcop
#endif

      real (kind=R8), allocatable, save ::
     .                             na_mean(:,:), te_mean(:),
     .                             ti_mean(:), ua_mean(:,:),
     .                             po_mean(:)
      real (kind=R8), allocatable, save ::
     .                             e_na(:,:), e_te(:),
     .                             e_ti(:), e_ua(:,:),
     .                             e_po(:)
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

      character lble*120, lbli*120, check*10
      character*10, save ::
     *   b2faveri_version = '00.000.000'

      contains
     

      subroutine running_average (nCv,ns,pl,dv,
     &           sna0_eir, smo0_eir, she0_eir, shi0_eir)

      use b2mod_types
      use b2mod_constants
      use b2mod_plasma
      use b2us_plasma
#ifdef B25_EIRENE
      use eirmod_infcop
#endif

      implicit none
c input
      integer, intent(in) :: nCv, ns
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      real (kind=R8), intent(in) ::
     *   sna0_eir(nCv,0:1,0:ns-1), smo0_eir(nCv,0:3,0:ns-1), 
     *   she0_eir(nCv,0:3), shi0_eir(nCv,0:3)
c local
      integer :: itp, is
      real (kind=R8) ::
     *   sna_tot(nCv,0:ns-1), smo_tot(nCv,0:ns-1), 
     *   she_tot(nCv), shi_tot(nCv)

      call subini ('running_average')

      naver = naver + 1
      ntotdt = ntotdt + 1

      itp = 0

cnh   Calculation of the total sources
      do is = 0, ns-1
        sna_tot(:,is) = sna0_eir(:,0,is) + sna0_eir(:,1,is)*pl%na(:,is)
        smo_tot(:,is) = smo0_eir(:,0,is) + smo0_eir(:,1,is)*pl%ua(:,is)
     1    + smo0_eir(:,2,is)*pl%na(:,is)*mp*am(is) + 
     2    smo0_eir(:,3,is)*pl%ua(:,is)*pl%na(:,is)*mp*am(is) 
      end do
      she_tot = she0_eir(:,0) + she0_eir(:,1)*pl%te + 
     1  she0_eir(:,2)*dv%ne + she0_eir(:,3)*pl%te*dv%ne
      shi_tot = shi0_eir(:,0) + shi0_eir(:,1)*pl%ti +
     1  shi0_eir(:,2)*dv%ni(:,0) + shi0_eir(:,3)*pl%ti*dv%ni(:,0)


*  calculate averaged plasma profiles and sources
      call batch_average (nCv*ns, pl%na, na_mean, itp, naver)
      call batch_average (nCv,    pl%te, te_mean, itp, naver)
      call batch_average (nCv,    pl%ti, ti_mean, itp, naver)
      call batch_average (nCv,    pl%po, po_mean, itp, naver)
      call batch_average (nCv*ns, pl%ua, ua_mean, itp, naver)

      call batch_average (nCv*ns, sna_tot, sna_mean, itp, naver)
      call batch_average (nCv,    she_tot, she_mean, itp, naver)
      call batch_average (nCv,    shi_tot, shi_mean, itp, naver)
      call batch_average (nCv*ns, smo_tot, smo_mean, itp, naver)

*  calculate squared averaged plasma profiles
      call batch_average_sq (nCv*ns, pl%na, e_na, itp, naver)
      call batch_average_sq (nCv,    pl%te, e_te, itp, naver)
      call batch_average_sq (nCv,    pl%ti, e_ti, itp, naver)
      call batch_average_sq (nCv,    pl%po, e_po, itp, naver)
      call batch_average_sq (nCv*ns, pl%ua, e_ua, itp, naver)
      
      call batch_average_sq (nCv*ns, sna_tot, e_sna, itp, naver)
      call batch_average_sq (nCv,    she_tot, e_she, itp, naver)
      call batch_average_sq (nCv,    shi_tot, e_shi, itp, naver)
      call batch_average_sq (nCv*ns, smo_tot, e_smo, itp, naver)

      call subend ()
      return
      end subroutine running_average



      subroutine run_av_init(iav_cont, nCv, ns, lpr, ler)

      implicit none
c input
      integer, intent(in) :: nCv, ns, iav_cont
      logical, intent(in) :: lpr, ler
c local variables
      integer :: io, iun, idum(0:1), iex
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)


      call subini ('run_av_init')

      if (.not.allocated(na_mean)) then
         allocate (na_mean(nCv,0:ns-1))
         allocate (ua_mean(nCv,0:ns-1))
         allocate (te_mean(nCv))
         allocate (ti_mean(nCv))
         allocate (po_mean(nCv))
         allocate (e_na(nCv,0:ns-1))
         allocate (e_ua(nCv,0:ns-1))
         allocate (e_te(nCv))
         allocate (e_ti(nCv))
         allocate (e_po(nCv))

         allocate (sna_mean(nCv,0:ns-1))
         allocate (smo_mean(nCv,0:ns-1))
         allocate (she_mean(nCv))
         allocate (shi_mean(nCv))
         allocate (e_sna(nCv,0:ns-1))
         allocate (e_smo(nCv,0:ns-1))
         allocate (e_she(nCv))
         allocate (e_shi(nCv))

      end if

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
      she_mean = 0.0_R8
      shi_mean = 0.0_R8
      smo_mean = 0.0_R8
      e_sna = 0.0_R8
      e_she = 0.0_R8
      e_shi = 0.0_R8
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
              call cfruin (iun, 2, idum, 'nCv,ns')
              call cfruch (iun, 120, lbli, 'label')
              call b2ruzd (iun, b2faveri_version, ns,
     &          zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)

              call cfruin (iun, 1, idum, 'naver')
              naver = idum(0)
              call cfruin (iun, 1, idum, 'ntotdt')
              ntotdt = idum(0)
              call cfrure (iun, nCv*ns, na_mean, 'na_mean')
              call cfrure (iun, nCv,    te_mean, 'te_mean')
              call cfrure (iun, nCv,    ti_mean, 'ti_mean')
              call cfrure (iun, nCv,    po_mean, 'po_mean')
              call cfrure (iun, nCv*ns, ua_mean, 'ua_mean')

              call cfrure (iun, nCv*ns, sna_mean, 'sna_mean')
              call cfrure (iun, nCv,    she_mean, 'she_mean')
              call cfrure (iun, nCv,    shi_mean, 'shi_mean')
              call cfrure (iun, nCv*ns, smo_mean, 'smo_mean')
            end if

*  read 2d averaged squared plasma profiles from b2fstati_run_var

            if (ler) then
              call cfrure (iun, nCv*ns, e_na, 'e_na')
              call cfrure (iun, nCv,    e_te, 'e_te')
              call cfrure (iun, nCv,    e_ti, 'e_ti')
              call cfrure (iun, nCv,    e_po, 'e_po')
              call cfrure (iun, nCv*ns, e_ua, 'e_ua')

              call cfrure (iun, nCv*ns, e_sna, 'e_sna')
              call cfrure (iun, nCv,    e_she, 'e_she')
              call cfrure (iun, nCv,    e_shi, 'e_shi')
              call cfrure (iun, nCv*ns, e_smo, 'e_smo')
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
      deallocate (smo_mean)
      deallocate (she_mean)
      deallocate (shi_mean)
      deallocate (e_sna)
      deallocate (e_smo)
      deallocate (e_she)
      deallocate (e_shi)

      return
      end subroutine run_av_fin

      subroutine run_av_save (savefile, nCv, ns, lpr, ler)
      implicit none
      integer, intent(in) :: nCv, ns
      logical, intent(in) :: lpr, ler
      character(*) :: savefile
c local variables
      integer :: idum(0:1), iun

      call subini ('run_av_save')


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

         call cfwure (iun, nCv*ns, na_mean, 'na_mean')
         call cfwure (iun, nCv, te_mean, 'te_mean')
         call cfwure (iun, nCv, ti_mean, 'ti_mean')
         call cfwure (iun, nCv, po_mean, 'po_mean')
         call cfwure (iun, nCv*ns, ua_mean, 'ua_mean')

         call cfwure (iun, nCv*ns, sna_mean, 'sna_mean')
         call cfwure (iun, nCv, she_mean, 'she_mean')
         call cfwure (iun, nCv, shi_mean, 'shi_mean')
         call cfwure (iun, nCv*ns, smo_mean, 'smo_mean')

      end if

*  write 2d squared averaged plasma profiles

      if (ler) then

         call cfwure (iun, nCv*ns, e_na, 'e_na')
         call cfwure (iun, nCv,    e_te, 'e_te')
         call cfwure (iun, nCv,    e_ti, 'e_ti')
         call cfwure (iun, nCv,    e_po, 'e_po')
         call cfwure (iun, nCv*ns, e_ua, 'e_ua')

         call cfwure (iun, nCv*ns, e_sna, 'e_sna')
         call cfwure (iun, nCv,    e_she, 'e_she')
         call cfwure (iun, nCv,    e_shi, 'e_shi')
         call cfwure (iun, nCv*ns, e_smo, 'e_smo')

      end if
!
      close (iun)
      write (*,*) ' file ',savefile,' written '

      call subend ()
      return
      end subroutine run_av_save


      subroutine run_av_get_plasma (na_arg, te_arg, ti_arg, ua_arg, 
     &                          po_arg, nCv, ns)
      implicit none
      integer, intent(in) :: nCv, ns
      real (kind=R8), intent(inout) ::
     *  na_arg(nCv,0:ns-1), ua_arg(nCv,0:ns-1),
     &  te_arg(nCv), ti_arg(nCv), po_arg(nCv)

      call subini ('run_av_get_plasma')

      na_arg = na_mean
      ua_arg = ua_mean
      te_arg = te_mean
      ti_arg = ti_mean
      po_arg = po_mean

      call subend ()
      return
      end subroutine run_av_get_plasma


      subroutine run_av_get_sources (sna_out, smo_out, she_out, shi_out,
     &                          nCv, ns)
      implicit none
      integer, intent(in) :: nCv, ns
      real (kind=R8), intent(out) ::
     *  sna_out(nCv,0:ns-1), smo_out(nCv,0:ns-1),
     &  she_out(nCv), shi_out(nCv)

      call subini ('run_av_get_sources')

      sna_out = sna_mean
      smo_out = smo_mean
      she_out = she_mean
      shi_out = shi_mean

      call subend ()

      end subroutine run_av_get_sources

      end module b2mod_running_average


