

      subroutine read_averages(na_arg, te_arg, ti_arg, ua_arg, po_arg,
     &                         nx, ny, ns)
      
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_geo
      use b2mod_plasma
      
      implicit none
      
      integer, intent(in) :: nx, ny, ns

      real (kind=R8), intent(inout) ::
     *  na_arg(-1:nx,-1:ny,0:ns-1), ua_arg(-1:nx,-1:ny,0:ns-1),
     &  te_arg(-1:nx,-1:ny), ti_arg(-1:nx,-1:ny), po_arg(-1:nx,-1:ny)

      real (kind=R8), intent(inout) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)

      real*8, allocatable, save :: na_mean(:,:,:), te_mean(:,:),
     .                             ti_mean(:,:), ua_mean(:,:,:),
     .                             po_mean(:,:)
      real*8, allocatable, save :: sna_mean(:,:,:,:),
     .                             sne_mean(:,:,:),
     .                             smo_mean(:,:,:,:),
     .                             she_mean(:,:,:),
     .                             shi_mean(:,:,:),
     .                             sch_mean(:,:,:)
      real*8, allocatable, save :: e_na(:,:,:), e_te(:,:),
     .                             e_ti(:,:), e_ua(:,:,:),
     .                             e_po(:,:)
      real*8, allocatable, save :: e_sna(:,:,:,:),
     .                             e_sne(:,:,:),
     .                             e_smo(:,:,:,:),
     .                             e_she(:,:,:),
     .                             e_shi(:,:,:),
     .                             e_sch(:,:,:)
      real*8, allocatable, save :: sig_na(:,:,:), sig_te(:,:),
     .                             sig_ti(:,:), sig_ua(:,:,:),
     .                             sig_po(:,:)
      real*8, save :: sg_na,sg_te,sg_ti,sg_ua, sg_po
      real*8 :: sfee, sfei, sfni

      integer :: io, ix, iy, is, istra, iun, k
      integer, save :: jxi, jxa, jsep
      integer, save :: naver, ntotdt

      call get_jsep(nx,ny,jxi,jxa,jsep)

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
      allocate (sig_na(-1:nx,-1:ny,0:ns-1))
      allocate (sig_ua(-1:nx,-1:ny,0:ns-1))
      allocate (sig_te(-1:nx,-1:ny))
      allocate (sig_ti(-1:nx,-1:ny))
      allocate (sig_po(-1:nx,-1:ny))

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
      sig_na = 0.0_R8
      sig_te = 0.0_R8
      sig_ti = 0.0_R8
      sig_ua = 0.0_R8
      sig_po = 0.0_R8

      sna_mean = 0.0_R8
      sne_mean = 0.0_R8
      smo_mean = 0.0_R8
      she_mean = 0.0_R8
      shi_mean = 0.0_R8
      sch_mean = 0.0_R8
      e_sna = 0.0_R8
      e_sne = 0.0_R8
      e_smo = 0.0_R8
      e_she = 0.0_R8
      e_shi = 0.0_R8
      e_sch = 0.0_R8
 
      io = 0
      iun = 3
      
      open (newunit=iun,file='aver.dat',access='SEQUENTIAL',
     .      form='FORMATTED',status='OLD',iostat=io)

      if (io == 0) then
         read (iun,*) naver, ntotdt
         read (iun,*) 
     .        (((na_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
         read (iun,*) ((te_mean(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) ((ti_mean(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) ((po_mean(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) 
     .        (((ua_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
         read (iun,*) (((e_na(ix,iy,is),
     .                   ix=-1,nx),iy=-1,ny),is=0,ns-1)
         read (iun,*) ((e_te(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) ((e_ti(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) ((e_po(ix,iy),ix=-1,nx),iy=-1,ny)
         read (iun,*) (((e_ua(ix,iy,is),
     .                 ix=-1,nx),iy=-1,ny),is=0,ns-1)
         read (iun,*) ((((sna_mean(ix,iy,k,is),
     .                 ix=-1,nx),iy=-1,ny),k=0,1),is=0,ns-1)
         read (iun,*) (((sne_mean(ix,iy,k),
     .                 ix=-1,nx),iy=-1,ny),k=0,1)
         read (iun,*) (((she_mean(ix,iy,k),
     .                 ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) (((shi_mean(ix,iy,k),
     .                 ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) (((sch_mean(ix,iy,k),
     .                 ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) ((((smo_mean(ix,iy,k,is),
     .                 ix=-1,nx),iy=-1,ny),k=0,3),is=0,ns-1)
         read (iun,*) ((((e_sna(ix,iy,k,is),
     .                ix=-1,nx),iy=-1,ny),k=0,1),is=0,ns-1)
         read (iun,*) (((e_sne(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,1)
         read (iun,*) (((e_she(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) (((e_shi(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) (((e_sch(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
         read (iun,*) ((((e_smo(ix,iy,k,is),
     .                ix=-1,nx),iy=-1,ny),k=0,3),is=0,ns-1)
         close (iun)
      end if
      
      na_arg(-1:nx,-1:ny,0:ns-1) = na_mean(-1:nx,-1:ny,0:ns-1)
      ua_arg(-1:nx,-1:ny,0:ns-1) = ua_mean(-1:nx,-1:ny,0:ns-1)
      te_arg(-1:nx,-1:ny) = te_mean(-1:nx,-1:ny)
      ti_arg(-1:nx,-1:ny) = ti_mean(-1:nx,-1:ny)
      po_arg(-1:nx,-1:ny) = po_mean(-1:nx,-1:ny)

      te_mean = te_mean/qe
      ti_mean = ti_mean/qe
      e_te = e_te /qe/qe
      e_ti = e_ti /qe/qe

      sig_na(-1:nx,-1:ny,0:ns-1) = e_na(-1:nx,-1:ny,0:ns-1)
     .                            - na_arg(-1:nx,-1:ny,0:ns-1)**2
      sig_ua(-1:nx,-1:ny,0:ns-1) = e_ua(-1:nx,-1:ny,0:ns-1)
     .                            - ua_arg(-1:nx,-1:ny,0:ns-1)**2
      sig_te(-1:nx,-1:ny) = e_te(-1:nx,-1:ny)
     .                        - te_mean(-1:nx,-1:ny)**2
      sig_ti(-1:nx,-1:ny) = e_ti(-1:nx,-1:ny)
     .                        - ti_mean(-1:nx,-1:ny)**2

      return

      
      entry source_averages (sna0, smo0, she0, shi0,
     &                       sch0, sne0, nx, ny, ns)


      sna0(-1:nx,-1:ny,0:1,0:ns-1) = sna_mean(-1:nx,-1:ny,0:1,0:ns-1) 
      sne0(-1:nx,-1:ny,0:1)        = sne_mean(-1:nx,-1:ny,0:1) 
      smo0(-1:nx,-1:ny,0:3,0:ns-1) = smo_mean(-1:nx,-1:ny,0:3,0:ns-1) 
      she0(-1:nx,-1:ny,0:3)        = she_mean(-1:nx,-1:ny,0:3) 
      shi0(-1:nx,-1:ny,0:3)        = shi_mean(-1:nx,-1:ny,0:3) 
      sch0(-1:nx,-1:ny,0:3)        = sch_mean(-1:nx,-1:ny,0:3) 

      return


      entry write_averages (nx, ny, ns)

      iun = 4
      open (newunit=iun,file='prof.out',form='FORMATTED',
     .      access='SEQUENTIAL')

      write (iun,*) ' averaged over ',naver,' timesteps'
      write (iun,*) ' total number of timesteps done ',ntotdt
      write (iun,*)

      if (nnreg(0) > 1) then
        iy = jsep + 1
        write (iun,*) ' profiles along jsep+1 = ',iy
      else
        iy = ny/2
        write (iun,*) ' profiles along NY/2 = ',iy
      endif 

      write (iun,'(/a6,12a24)')
     .           'ix','na','sqrt(sig_na)','te','sqrt(sig_te)',
     .           'ti','sqrt(sig_ti)','ua','sqrt(sig_ua)',
     .           'po','sqrt(sig_po)',
     .           'vol', 'dx'
      do ix = -1, nx          
        write (iun,'(i6,12es24.16)') ix,
     .            na(ix,iy,1), sqrt(sig_na(ix,iy,1)),
     .            te_mean(ix,iy),sqrt(sig_te(ix,iy)),
     .            ti_mean(ix,iy),sqrt(sig_ti(ix,iy)),
     .            ua(ix,iy,1), sqrt(sig_ua(ix,iy,1)),
     .            po(ix,iy), sqrt(sig_po(ix,iy)),
     .            vol(ix,iy),hx(ix,iy)
      end do

      ix = -1
      write (iun,'(//a,i4)') ' profiles along IX = ',ix

      write (iun,'(/a6,18a24)')
     .           'iy','na','sqrt(sig_na)','te','sqrt(sig_te)',
     .           'ti','sqrt(sig_ti)','ua','sqrt(sig_ua)',
     .           'po','sqrt(sig_po)',
     .           'fna/gs','fhe/gs','fhi/gs',
     .           'vol','hy','gs'
!    .          ,'qcos(0)','qcos(1)'
      sfee = 0.0_R8
      sfei = 0.0_R8
      sfni = 0.0_R8
      do iy = -1, ny       
        write (iun,'(i6,18es24.16)') iy,
     .            na(ix,iy,1), sqrt(sig_na(ix,iy,1)),
     .            te_mean(ix,iy),sqrt(sig_te(ix,iy)),
     .            ti_mean(ix,iy),sqrt(sig_ti(ix,iy)),
     .            ua(ix,iy,1), sqrt(sig_ua(ix,iy,1)),
     .            po(ix,iy), sqrt(sig_po(ix,iy)),
     .            fna(ix+1,iy,0,1)/gs(ix+1,iy,0),
     .            fhe(ix+1,iy,0)/gs(ix+1,iy,0),
     .            fhi(ix+1,iy,0)/gs(ix+1,iy,0),
     .            vol(ix,iy), hy(ix,iy),
     .            gs(ix+1,iy,0)
!    .           ,qcos(ix,iy,0),qcos(ix,iy,1)
        SFEE=SFEE+FHE(IX+1,IY,0)
        SFEI=SFEI+FHI(IX+1,IY,0)
        SFNI=SFNI+FNA(IX+1,IY,0,1)
      end do

      write (iun,'(/1X,3(A,ES12.4))') 
     .    ' SFNI = ',SFNI,' SFEE = ',SFEE,' SFEI = ',SFEI

      ix = nx
      write (iun,'(//a,i4)') ' profiles along IX = ',ix

      write (iun,'(/a6,18a24)')
     .           'iy','na','sqrt(sig_na)','te','sqrt(sig_te)',
     .           'ti','sqrt(sig_ti)','ua','sqrt(sig_ua)',
     .           'po','sqrt(sig_po)',
     .           'fna/gs','fhe/gs','fhi/gs',
     .           'vol','hy','gs'
!    .          ,'qcos(0)','qcos(1)'
      sfee = 0.0_R8
      sfei = 0.0_R8
      sfni = 0.0_R8
      do iy = -1, ny         
        write (iun,'(i6,18es24.16)') iy,
     .            na(ix,iy,1), sqrt(sig_na(ix,iy,1)),
     .            te_mean(ix,iy),sqrt(sig_te(ix,iy)),
     .            ti_mean(ix,iy),sqrt(sig_ti(ix,iy)),
     .            ua(ix,iy,1), sqrt(sig_ua(ix,iy,1)),
     .            po(ix,iy), sqrt(sig_po(ix,iy)),
     .            fna(ix,iy,0,1)/gs(ix,iy,0),
     .            fhe(ix,iy,0)/gs(ix,iy,0),
     .            fhi(ix,iy,0)/gs(ix,iy,0),
     .            vol(ix,iy),hy(ix,iy),
     -            gs(ix,iy,0)
!    .           ,qcos(ix-1,iy,0),qcos(ix-1,iy,1)
        SFEE=SFEE-FHE(IX,IY,0)
        SFEI=SFEI-FHI(IX,IY,0)
        SFNI=SFNI-FNA(IX,IY,0,1)
      end do

      write (iun,'(/1X,3(A,ES12.4))') 
     .    ' SFNI = ',SFNI,' SFEE = ',SFEE,' SFEI = ',SFEI

      write (iun,*) 'min_sig_na ',minval(sig_na(-1:nx,-1:ny,1))
      write (iun,*) 'min_sig_te ',minval(sig_te(-1:nx,-1:ny))
      write (iun,*) 'min_sig_ti ',minval(sig_ti(-1:nx,-1:ny))
      write (iun,*) 'min_sig_ua ',minval(sig_ua(-1:nx,-1:ny,1))
      write (iun,*) 'min_sig_po ',minval(sig_po(-1:nx,-1:ny))

      sg_na = sum(sqrt(abs(sig_na(-1:nx,-1:ny,1)))/
     .            na(nx/2,-1,1))/(nx+2)/(ny+2)
      sg_te = sum(sqrt(abs(sig_te(-1:nx,-1:ny)))/
     .            te_mean(nx/2,-1))/(nx+2)/(ny+2)
      sg_ti = sum(sqrt(abs(sig_ti(-1:nx,-1:ny)))/ 
     .            ti_mean(nx/2,-1))/(nx+2)/(ny+2)
      sg_ua = sum(sqrt(abs(sig_ua(-1:nx,-1:ny,1)))/
     .            abs(ua(0,0,1)))/(nx+2)/(ny+2)
      sg_po = sum(sqrt(abs(sig_po(-1:nx,-1:ny)))/
     .            po(nx/2,-1))/(nx+2)/(ny+2)
      write (iun,'(//1x,A,5es12.4)') 'sig_na,ee,ei,ua,po ',
     .          sg_na,sg_te,sg_ti,sg_ua,sg_po
      
      close (unit=iun)

      deallocate(na_mean)
      deallocate(ua_mean)
      deallocate(te_mean)
      deallocate(ti_mean)
      deallocate(po_mean)
      deallocate(e_na)
      deallocate(e_ua)
      deallocate(e_te)
      deallocate(e_ti)
      deallocate(e_po)
      deallocate(sig_na)
      deallocate(sig_ua)
      deallocate(sig_te)
      deallocate(sig_ti)
      deallocate(sig_po)
      deallocate(sna_mean)
      deallocate(sne_mean)
      deallocate(smo_mean)
      deallocate(she_mean)
      deallocate(shi_mean)
      deallocate(sch_mean)
      deallocate(e_sna)
      deallocate(e_sne)
      deallocate(e_smo)
      deallocate(e_she)
      deallocate(e_shi)
      deallocate(e_sch)

      write (6,*) ' read_averages finished '

      return
      end subroutine read_averages
