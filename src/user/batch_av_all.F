
      subroutine batch_av_all (itim,ntim_batch,nx,ny,ns)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_sources

      implicit none

      integer, intent(in) :: itim,ntim_batch,nx, ny, ns
      character(*), intent(in) :: batch_name
 
      real (kind=R8), allocatable, save ::
     .                             na_mean(:,:,:), te_mean(:,:),
     .                             ti_mean(:,:), ua_mean(:,:,:),
     .                             po_mean(:,:)

      real (kind=R8), allocatable, save :: 
     .                             sna_mean(:,:,:,:),
     .                             sne_mean(:,:,:),
     .                             smo_mean(:,:,:,:),
     .                             she_mean(:,:,:),
     .                             shi_mean(:,:,:),
     .                             sch_mean(:,:,:)

      integer ::  iun, n2, itp
      integer, save :: naver, ntotdt

      n2 = (nx+2)*(ny+2)

      itp = itim + 1

      call batch_average (n2*ns, na, na_mean, itp, ntim_batch)
      call batch_average (n2,    te, te_mean, itp, ntim_batch)
      call batch_average (n2,    ti, ti_mean, itp, ntim_batch)
      call batch_average (n2,    po, po_mean, itp, ntim_batch)
      call batch_average (n2*ns, ua, ua_mean, itp, ntim_batch)
      
      call batch_average (n2*2*ns, sna, sna_mean, itp, ntim_batch)
      call batch_average (n2*2,    sne, sne_mean, itp, ntim_batch)
      call batch_average (n2*4,    she, she_mean, itp, ntim_batch)
      call batch_average (n2*4,    shi, shi_mean, itp, ntim_batch)
      call batch_average (n2*4,    sch, sch_mean, itp, ntim_batch)
      call batch_average (n2*4*ns, smo, smo_mean, itp, ntim_batch)
       
      return

      entry batch_av_all_init(nx,ny,ns)

      if (.not.allocated(na_mean)) then

         allocate (na_mean(-1:nx,-1:ny,0:ns-1))
         allocate (ua_mean(-1:nx,-1:ny,0:ns-1))      
         allocate (te_mean(-1:nx,-1:ny))        
         allocate (ti_mean(-1:nx,-1:ny))
         allocate (po_mean(-1:nx,-1:ny))

ctb: the sources are now stored in linearized form I propose to store the total part only 
ctb..as the mean of the product does not equal the product of the means etc.
         allocate (sna_mean(-1:nx,-1:ny,0:1,0:ns-1))
         allocate (sne_mean(-1:nx,-1:ny,0:1))
         allocate (smo_mean(-1:nx,-1:ny,0:3,0:ns-1))
         allocate (she_mean(-1:nx,-1:ny,0:3))
         allocate (shi_mean(-1:nx,-1:ny,0:3))
         allocate (sch_mean(-1:nx,-1:ny,0:3))

      end if

      na_mean = 0.d0
      te_mean = 0.d0
      ti_mean = 0.d0
      ua_mean = 0.d0
      po_mean = 0.d0
      
      sna_mean = 0.d0
      sne_mean = 0.d0
      she_mean = 0.d0
      shi_mean = 0.d0
      sch_mean = 0.d0
      smo_mean = 0.d0
      
      return

      entry batch_av_all_fin

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

      return

      entry batch_av_all_save (batch_name,nx,ny,ns)

      n2 = (nx+2)*(ny+2)
      open (newunit=iun,file=batch_name,access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')

!      call cfwuin (iun, 1, naver, 'naver')
!      call cfwuin (iun, 1, ntotdt, 'ntotdt')
      
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
      
      close(iun)
      return
      
      end subroutine batch_av_all
      
