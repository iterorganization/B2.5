
      subroutine profile_average (laststep,nx,ny,ns)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_sources

      implicit none

      integer, intent(in) :: nx, ny, ns
      logical, intent(in) :: laststep
 
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

      integer :: io, ix, iy, is, istra, iun, k
      real (kind=R8) :: gold, gnew
      integer, save :: ifirst=0, iaver_restart=0
      integer, save :: naver, ntotdt

      if (ifirst == 0) then
         call ipgeti ('average_restart', iaver_restart)
         naver = 0
         ntotdt = 0
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
         
         iun = 3
         open (newunit=iun,file='aver.dat',access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD',iostat=io)

         if (io == 0) then
            if (iaver_restart > 0) then
               read (iun,*) naver, ntotdt
               naver = 0
               close (iun)
            else
               read (iun,*) naver, ntotdt
               read (iun,*) 
     .              (((na_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
               read (iun,*) ((te_mean(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) ((ti_mean(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) ((po_mean(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) 
     .              (((ua_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
               read (iun,*) (((e_na(ix,iy,is),
     .                            ix=-1,nx),iy=-1,ny),is=0,ns-1)
               read (iun,*) ((e_te(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) ((e_ti(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) ((e_po(ix,iy),ix=-1,nx),iy=-1,ny)
               read (iun,*) (((e_ua(ix,iy,is),
     .                            ix=-1,nx),iy=-1,ny),is=0,ns-1)
               read (iun,*) ((((sna_mean(ix,iy,k,is),
     .                ix=-1,nx),iy=-1,ny),k=0,1),is=0,ns-1)
               read (iun,*) (((sne_mean(ix,iy,k),
     .                                 ix=-1,nx),iy=-1,ny),k=0,1)
               read (iun,*) (((she_mean(ix,iy,k),
     .                                 ix=-1,nx),iy=-1,ny),k=0,3)
               read (iun,*) (((shi_mean(ix,iy,k),
     .                                 ix=-1,nx),iy=-1,ny),k=0,3)
               read (iun,*) (((sch_mean(ix,iy,k),
     .                                 ix=-1,nx),iy=-1,ny),k=0,3)
               read (iun,*) ((((smo_mean(ix,iy,k,is),
     .                  ix=-1,nx),iy=-1,ny),k=0,3),is=0,ns-1)
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
         end if
         ifirst = 1
      end if

      if (.not.laststep) then
        naver = naver + 1
        ntotdt = ntotdt + 1
         
        gnew = 1.0_R8 / naver
        gold = (naver - 1) * gnew

        na_mean(-1:nx,-1:ny,0:ns-1) = gold*na_mean(-1:nx,-1:ny,0:ns-1) + 
     .                                gnew*na(-1:nx,-1:ny,0:ns-1)
        te_mean(-1:nx,-1:ny)        = gold*te_mean(-1:nx,-1:ny) + 
     .                                gnew*te(-1:nx,-1:ny)
        ti_mean(-1:nx,-1:ny)        = gold*ti_mean(-1:nx,-1:ny) + 
     .                                gnew*ti(-1:nx,-1:ny)
        po_mean(-1:nx,-1:ny)        = gold*po_mean(-1:nx,-1:ny) + 
     .                                gnew*po(-1:nx,-1:ny)
        ua_mean(-1:nx,-1:ny,0:ns-1) = gold*ua_mean(-1:nx,-1:ny,0:ns-1) + 
     .                                gnew*ua(-1:nx,-1:ny,0:ns-1)

        e_na(-1:nx,-1:ny,0:ns-1) = gold*e_na(-1:nx,-1:ny,0:ns-1) + 
     .                             gnew*na(-1:nx,-1:ny,0:ns-1)**2
        e_te(-1:nx,-1:ny)        = gold*e_te(-1:nx,-1:ny) + 
     .                             gnew*te(-1:nx,-1:ny)**2
        e_ti(-1:nx,-1:ny)        = gold*e_ti(-1:nx,-1:ny) + 
     .                             gnew*ti(-1:nx,-1:ny)**2
        e_po(-1:nx,-1:ny)        = gold*e_po(-1:nx,-1:ny) + 
     .                             gnew*po(-1:nx,-1:ny)**2
        e_ua(-1:nx,-1:ny,0:ns-1) = gold*e_ua(-1:nx,-1:ny,0:ns-1) + 
     .                             gnew*ua(-1:nx,-1:ny,0:ns-1)**2
        
        sna_mean(-1:nx,-1:ny,0:1,0:ns-1)  = 
     .        gold*sna_mean(-1:nx,-1:ny,0:1,0:ns-1) + 
     .        gnew*sna(-1:nx,-1:ny,0:1,0:ns-1) 
        sne_mean(-1:nx,-1:ny,0:1) = 
     .        gold*sne_mean(-1:nx,-1:ny,0:1) + 
     .        gnew*sne(-1:nx,-1:ny,0:1)  
        she_mean(-1:nx,-1:ny,0:3) = 
     .        gold*she_mean(-1:nx,-1:ny,0:3) + 
     .        gnew*she(-1:nx,-1:ny,0:3)  
        shi_mean(-1:nx,-1:ny,0:3) = 
     .        gold*shi_mean(-1:nx,-1:ny,0:3) + 
     .        gnew*shi(-1:nx,-1:ny,0:3) 
        sch_mean(-1:nx,-1:ny,0:3) = 
     .        gold*sch_mean(-1:nx,-1:ny,0:3) + 
     .        gnew*sch(-1:nx,-1:ny,0:3)  
        smo_mean(-1:nx,-1:ny,0:3,0:ns-1)  = 
     .        gold*smo_mean(-1:nx,-1:ny,0:3,0:ns-1) + 
     .        gnew*smo(-1:nx,-1:ny,0:3,0:ns-1) 
 
        e_sna(-1:nx,-1:ny,0:1,0:ns-1)  = 
     .        gold*e_sna(-1:nx,-1:ny,0:1,0:ns-1) + 
     .        gnew*sna(-1:nx,-1:ny,0:1,0:ns-1) 
        e_sne(-1:nx,-1:ny,0:1) = 
     .        gold*e_sne(-1:nx,-1:ny,0:1) + 
     .        gnew*sne(-1:nx,-1:ny,0:1)  
        e_she(-1:nx,-1:ny,0:3) = 
     .        gold*e_she(-1:nx,-1:ny,0:3) + 
     .        gnew*she(-1:nx,-1:ny,0:3)  
        e_shi(-1:nx,-1:ny,0:3) = 
     .        gold*e_shi(-1:nx,-1:ny,0:3) + 
     .        gnew*shi(-1:nx,-1:ny,0:3)  
        e_sch(-1:nx,-1:ny,0:3) = 
     .        gold*e_sch(-1:nx,-1:ny,0:3) + 
     .        gnew*sch(-1:nx,-1:ny,0:3) 
        e_smo(-1:nx,-1:ny,0:3,0:ns-1)  = 
     .        gold*e_smo(-1:nx,-1:ny,0:3,0:ns-1) + 
     .        gnew*smo(-1:nx,-1:ny,0:3,0:ns-1) 


!pb      if (laststep) then
      else
        iun = 3
        open (newunit=iun,file='aver.dat',access='SEQUENTIAL',
     .          form='FORMATTED',status='REPLACE')
        write (iun,'(2i9)') naver, ntotdt
        write (iun,'(6es24.16)') 
     .        (((na_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .        ((te_mean(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        ((ti_mean(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        ((po_mean(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        (((ua_mean(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .        (((e_na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .        ((e_te(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        ((e_ti(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        ((e_po(ix,iy),ix=-1,nx),iy=-1,ny)
        write (iun,'(6es24.16)') 
     .        (((e_ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .    ((((sna_mean(ix,iy,k,is),ix=-1,nx),iy=-1,ny),k=0,1),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .        (((sne_mean(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,1)
        write (iun,'(6es24.16)') 
     .        (((she_mean(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .        (((shi_mean(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .        (((sch_mean(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .    ((((smo_mean(ix,iy,k,is),ix=-1,nx),iy=-1,ny),k=0,3),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .    ((((e_sna(ix,iy,k,is),ix=-1,nx),iy=-1,ny),k=0,1),is=0,ns-1)
        write (iun,'(6es24.16)') 
     .        (((e_sne(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,1)
        write (iun,'(6es24.16)') 
     .        (((e_she(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .        (((e_shi(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .        (((e_sch(ix,iy,k),ix=-1,nx),iy=-1,ny),k=0,3)
        write (iun,'(6es24.16)') 
     .    ((((e_smo(ix,iy,k,is),ix=-1,nx),iy=-1,ny),k=0,3),is=0,ns-1)

        close(iun)

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

      end if

      return
      end subroutine profile_average
