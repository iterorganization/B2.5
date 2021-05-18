      module b2mod_impgyro

      use b2mod_external

      implicit none
      private

#ifdef USE_MPI
      include 'mpif.h'

      integer :: mpi_err,mpi_stat,mpi_dum
      integer :: mpi_rank,mpi_size,mpi_local
#else
#ifdef USE_IMPGYRO
      include 'mpif.h'

      integer :: mpi_err,mpi_stat,mpi_dum
      integer :: mpi_rank,mpi_size,mpi_local
#endif
#endif

      logical :: last
      real(kind=R8) , allocatable ::
     &     ep_ir (:,:),     ! electron particle source
     &                      !    by ionization/recombination
     &     im_cc (:,:,:),   ! ion momentum source by Coulomb collision
     &     ie_cc (:,:),     ! ion heat source by Coulomb collision
     &     ee_rd (:,:),     ! electron heat source by radiation
     &     b2znw (:,:,:),   ! impurity density
     &     b2zfw (:,:,:),   ! impurity mean velocity
     &     b2ztw (:,:,:)    ! impurity mean random energy

      public :: impgyro_init,impgyro_main,impgyro_last
      logical, public :: use_impgyro

      contains

      subroutine impgyro_init (nxd, nyd, nsd, ismain)

      use b2mod_types
      use b2mod_geo
      use b2mod_constants
      use b2mod_indirect
      use b2mod_elements
      implicit none
      integer nxd,nyd,nsd,ismain
      integer is,iss,ix,iy
      integer, save :: ikmn,ikmx
      real(kind=R8), allocatable :: pdum (:), pot(:)
      character*256 filename
      logical open_file
      external xertst, xerrab, open_file
      logical brestarted
      integer fd
 
*   ..test arguments
      call xertst (0.le.nxd.and.0.le.nyd, 'faulty argument nxd, nyd')
      call xertst (1.le.nsd, 'faulty argument nsd')
      call xertst (0.le.ismain.and.ismain.lt.nsd,
     * 'faulty argument ismain')

      filename='use_impgyro'
      call find_file(filename,use_impgyro)
      if ( .not. use_impgyro ) return

      write (*,*) 'b2ig: Entering impgyro_init'

      last = .false.

#ifdef USE_MPI
      write (*,*) 'b2ig: Calling MPI COMM FLAG=MPI'
      call mpi_comm_size (mpi_comm_world, mpi_size, mpi_err)
      call mpi_comm_rank (mpi_comm_world, mpi_rank, mpi_err)
      call mpi_comm_split (mpi_comm_world, mpi_undefined,
     &  1, mpi_local, mpi_err)
      call mpi_bcast (ismain, 1, mpi_integer4, 0, mpi_comm_world,
     &  mpi_err)

      if(mpi_rank.eq.0) then
      call mpi_recv (ikmn, 1, mpi_integer4, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ikmx, 1, mpi_integer4, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      end if

      write (*,*) 'b2ig: mpi_size =', mpi_size
      write (*,*) 'b2ig: ikmn, ikmx obtained from IMPGYRO', ikmn, ikmx
#else
#ifdef USE_IMPGYRO
      write (*,*) 'b2ig: Calling MPI COMM FLAG=IG'
      call mpi_comm_size (mpi_comm_world, mpi_size, mpi_err)
      call mpi_comm_rank (mpi_comm_world, mpi_rank, mpi_err)
      call mpi_comm_split (mpi_comm_world, mpi_undefined,
     &  1, mpi_local, mpi_err)
      call mpi_bcast (ismain, 1, mpi_integer4, 0, mpi_comm_world,
     &  mpi_err)

      if(mpi_rank.eq.0) then
      call mpi_recv (ikmn, 1, mpi_integer4, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ikmx, 1, mpi_integer4, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      end if

      write (*,*) 'b2ig: ikmn, ikmx obtained from IMPGYRO', ikmn, ikmx
#else
      ikmn = 0
      ikmx = 74
      write (*,*) 'b2ig: WARNING YOU ARE RUNNING NON-MPI IMPGYRO case'
      write (*,*) 'b2ig: You may have forgotten set_ig option!'
      write (*,*) 'b2ig: Setting default ikmn, ikmx', ikmn, ikmx
#endif
#endif

      allocate(ep_ir(0:nxd-1,0:nyd-1))
      allocate(im_cc(0:nxd-1,0:nyd-1,0:nsd-1))
      allocate(ie_cc(0:nxd-1,0:nyd-1))
      allocate(ee_rd(0:nxd-1,0:nyd-1))
      allocate(b2znw(ikmn:ikmx,0:nxd-1,0:nyd-1))
      allocate(b2zfw(ikmn:ikmx,0:nxd-1,0:nyd-1))
      allocate(b2ztw(ikmn:ikmx,0:nxd-1,0:nyd-1))

      ns_ext = ikmx - ikmn + 1

      if (allocated(zn_ext)) deallocate(zn_ext)
      allocate(zn_ext(0:ns_ext-1))
      if (allocated(za_ext)) deallocate(za_ext)
      allocate(za_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(za2_ext)) deallocate(za2_ext)
      allocate(za2_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(pt_ext)) deallocate(pt_ext)
      allocate(pt_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(na_ext)) deallocate(na_ext)
      allocate(na_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(ua_ext)) deallocate(ua_ext)
      allocate(ua_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(ta_ext)) deallocate(ta_ext)
      allocate(ta_ext(-1:nxd,-1:nyd,0:ns_ext-1))
      if (allocated(fa_ext)) deallocate(fa_ext)
      allocate(fa_ext(-1:nxd,-1:nyd,0:1,0:ns_ext-1))
      if (allocated(am_ext)) deallocate(am_ext)
      allocate(am_ext(0:ns_ext-1))
      if (allocated(text_ext)) deallocate(text_ext)
      allocate(text_ext(0:ns_ext-1))
      if (allocated(is_ext_neutral)) deallocate(is_ext_neutral)
      allocate(is_ext_neutral(0:ns_ext-1))
      if (allocated(pot)) deallocate (pot)
      allocate(pot(0:ns_ext-1))
      pot(:)=0.0_R8

cxpb
cxpb  Here we read the external data
cxpb   Species data
cxpb    za_ext, pt_ext, am_ext, text_ext

      do is = 0, ns_ext-1
        zn_ext(is) = 74.0_R8
        za_ext(:,:,is) = real(ikmn + is)
        am_ext(is) = 3.07832E-25_R8/mp
        if(ikmn+is.eq.0) then
          text_ext(is)=' W^{0}'
        else if(ikmn+is.lt.10) then
          write(text_ext(is),'(a,i1,a)') ' W^{+',ikmn+is,'}'
        else
          write(text_ext(is),'(a,i2,a)') ' W^{+',ikmn+is,'}'
        endif
        if (maxval(za_ext(-1:nxd,-1:nyd,is)).eq.0.0_R8) then
          is_ext_neutral(is)=.true.
        else
          is_ext_neutral(is)=.false.
        endif
        za2_ext(-1:nxd,-1:nyd,is)=za_ext(-1:nxd,-1:nyd,is)**2
      end do
      if(open_file(99,'ionization_potentials',
     1  'data.local|modules/B2.5/Database')) then
        do iss=1,73
          read(99,*)
        end do
        allocate(pdum(1:74))
        read(99,*) (pdum(iss),iss=1,74)
        pot(:)=0.0_R8
        do is = 0, ns_ext-1
          do iss = 1, ikmn+is
            pot(is) = pot(is) + pdum(iss)
          enddo
          pt_ext(:,:,is) = pot(is)
        enddo
        deallocate(pdum,pot)
      else
        call xerrab('Please provide ionization potentails!')
      endif

c     read the file if exists
      inquire (file = 'restart_ext', exist = brestarted)
      if (brestarted) then

        call my_open_file (fd, 'restart_ext' )
        read (fd) na_ext
        read (fd) ua_ext
        read (fd) ta_ext
        read (fd) sna_ext
        read (fd) sne_ext
        read (fd) smo_ext
        read (fd) she_ext
        read (fd) shi_ext
        read (fd) sch_ext
        close (fd)

      else

       ! every external quantities are set to 0 otherwise
        na_ext = 0.0_R8
        ua_ext = 0.0_R8
        ta_ext = 0.0_R8
        sna_ext = 0.0_R8
        sne_ext = 0.0_R8
        smo_ext = 0.0_R8
        she_ext = 0.0_R8
        shi_ext = 0.0_R8
        sch_ext = 0.0_R8

      end if

         call extend(na_ext,nxd,nyd,ns_ext,nxd,nyd,ns_ext)
         call extend(ua_ext,nxd,nyd,ns_ext,nxd,nyd,ns_ext)
         call extend(ta_ext,nxd,nyd,ns_ext,nxd,nyd,ns_ext)

cxpb
cxpb  Here are computed the derived quantities
cxpb
        ne_ext = 0.0_R8
        ni_ext = 0.0_R8
        ue_ext = 0.0_R8
        ne2_ext = 0.0_R8
        do is = 0, ns_ext-1
          ne_ext(-1:nxd,-1:nyd) = ne_ext(-1:nxd,-1:nyd) +
     .     za_ext(-1:nxd,-1:nyd,is) * na_ext(-1:nxd,-1:nyd,is)
          ne2_ext(-1:nxd,-1:nyd) = ne2_ext(-1:nxd,-1:nyd) +
     .     za2_ext(-1:nxd,-1:nyd,is) * na_ext(-1:nxd,-1:nyd,is)
          ni_ext(-1:nxd,-1:nyd,0) = ni_ext(-1:nxd,-1:nyd,0) +
     .     na_ext(-1:nxd,-1:nyd,is)
          if (.not.is_ext_neutral(is)) ni_ext(-1:nxd,-1:nyd,1) = 
     .     ni_ext(-1:nxd,-1:nyd,1) + na_ext(-1:nxd,-1:nyd,is)
cxpb Computing ua from fa is what we really want to do
        if (.false.) then
          do iy=-1,nyd
            do ix=-1,nxd
              if (leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nxd+1) then
                ua_ext(ix,iy,is) = 0.5_R8*
     .           (fa_ext(ix,iy,0,is)/pbs(ix,iy,0)+
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0))/
     .           na_ext(ix,iy,is)
              else if (leftix(ix,iy).eq.-2) then
                ua_ext(ix,iy,is) = 
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0)/
     .            na_ext(ix,iy,is)
              else if (rightix(ix,iy).eq.nxd+1) then
                ua_ext(ix,iy,is) = 
     .           fa_ext(ix,iy,0,is)/pbs(ix,iy,0)/na_ext(ix,iy,is)
              end if
            end do
          end do
        else
          do iy=-1,nyd
            do ix=-1,nxd
              if(leftix(ix,iy).ne.-2) then
                fa_ext(ix,iy,0,is)=((ua_ext(ix,iy,is)+
     .           ua_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           ((na_ext(ix,iy,is)+
     .           na_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           pbs(ix,iy,0)
                fhi_ext(ix,iy,0) = fhi_ext(ix,iy,0) +
     .           1.5_R8*fa_ext(ix,iy,0,is)*(ta_ext(ix,iy,is)+
     .           ta_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
              else
                fa_ext(ix,iy,0,is)=0.0_R8
              endif
            enddo
          enddo
          fa_ext(:,:,1,:) = 0.0_R8
        endif
cxpb
        ue_ext(-1:nxd,-1:nyd) = ue_ext(-1:nxd,-1:nyd) +
     .     za_ext(-1:nxd,-1:nyd,is) * ua_ext(-1:nxd,-1:nyd,is) *
     .     na_ext(-1:nxd,-1:nyd,is)
        end do
        do ix = -1, nxd
          do iy = -1, nyd
            if (ne_ext(ix,iy).gt.0.0_R8) then
              ue_ext(ix,iy) = ue_ext(ix,iy)/ne_ext(ix,iy)
            else
              ue_ext(ix,iy) = 0.0_R8
            endif
          enddo
        enddo

      write (*,*) 'b2ig: Leaving impgyro_init'

      end subroutine impgyro_init

      subroutine impgyro_main (nx,ny,ns,ns_ext,
     &  bb,na,ne,ua,te,ti,po,
     &  fna,fhe,fhi,fch,
     &  sna_ext,sne_ext,smo_ext,she_ext,shi_ext,sch_ext,
     &  na_ext,ua_ext,ta_ext)

      integer,intent(in) :: nx,ny,ns,ns_ext

      real(kind=R8),intent(in) ::
     &  bb (-1:nx,-1:ny,0:3),
     &  na (-1:nx,-1:ny,0:ns-1),
     &  ne (-1:nx,-1:ny),
     &  ua (-1:nx,-1:ny,0:ns-1),
     &  te (-1:nx,-1:ny),
     &  ti (-1:nx,-1:ny),
     &  po (-1:nx,-1:ny),
     &  fna (-1:nx,-1:ny,0:1,0:ns-1),
     &  fhe (-1:nx,-1:ny,0:1),
     &  fhi (-1:nx,-1:ny,0:1),
     &  fch (-1:nx,-1:ny,0:1)

      real(kind=R8),intent(out) ::
     &  sna_ext (-1:nx,-1:ny,0:ns-1),
     &  sne_ext (-1:nx,-1:ny),
     &  smo_ext (-1:nx,-1:ny,0:ns-1),
     &  she_ext (-1:nx,-1:ny),
     &  shi_ext (-1:nx,-1:ny),
     &  sch_ext (-1:nx,-1:ny),
     &  na_ext (-1:nx,-1:ny,0:ns_ext-1),
     &  ua_ext (-1:nx,-1:ny,0:ns_ext-1),
     &  ta_ext (-1:nx,-1:ny,0:ns_ext-1)

      integer :: ix ,iy ,is

      write (*,*) 'b2ig: Entering impgyro_main'

#ifdef USE_MPI
      write (*,*) 'b2ig: impgyro_main...flag=MPI'
      call mpi_bcast (last, 1, mpi_logical , 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (na, size(na), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ne, size(ne), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ua, size(ua), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (te, size (te), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ti, size (ti), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (po, size (po), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (fna, size (fna), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
#else
#ifdef USE_IMPGYRO
      write (*,*) 'b2ig: impgyro_main...flag=IG'
      call mpi_bcast (last, 1, mpi_logical , 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (na, size(na), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ne, size(ne), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ua, size(ua), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (te, size (te), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (ti, size (ti), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (po, size (po), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
      call mpi_bcast (fna, size (fna), mpi_real8, 0,
     &  mpi_comm_world, mpi_err)
#endif
#endif

      write (*,*) 'b2ig: running impgyro ...'

#ifdef USE_MPI
      if (mpi_rank.eq.0) then
      call mpi_recv (ep_ir, size(ep_ir), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (im_cc, size(im_cc), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ie_cc, size(ie_cc), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ee_rd, size(ee_rd), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2znw, size(b2znw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2zfw, size(b2zfw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2ztw, size(b2ztw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      end if
#else
#ifdef USE_IMPGYRO
      if (mpi_rank.eq.0) then
      call mpi_recv (ep_ir, size(ep_ir), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (im_cc, size(im_cc), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ie_cc, size(ie_cc), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (ee_rd, size(ee_rd), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2znw, size(b2znw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2zfw, size(b2zfw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      call mpi_recv (b2ztw, size(b2ztw), mpi_real8, 1, 0,
     &  mpi_comm_world, mpi_stat, mpi_err)
      end if
#else
      ep_ir = 0
      im_cc = 0
      ie_cc = 0
      ee_rd = 0
      b2znw = 0
      b2zfw = 0
      b2ztw = 0
      write (*,*) 'b2ig: IMPGYRO run failed!' 
      write (*,*) 'b2ig: Every external specie information set to zero.'
      write (*,*) 'b2ig: You may have forgotten set_ig option!'
#endif
#endif

      sna_ext = 0

      do is =0,ns-1
        do iy =0,ny-1
          do ix =0,nx-1
            smo_ext (ix,iy,is) = im_cc (ix,iy,is)*
     .       sign(1.0_R8,bb(ix,iy,0))
          end do
        end do
      end do

      do iy = 0, ny-1
        do ix = 0, nx-1
          sne_ext (ix,iy) = ep_ir (ix,iy)
          she_ext (ix,iy) = ee_rd (ix,iy)
          shi_ext (ix,iy) = ie_cc (ix,iy)
        end do
      end do

      sch_ext = 0

      do is = 0 , ns_ext - 1
        do iy = 0 , ny - 1
          do ix = 0 , nx - 1
            na_ext (ix,iy,is) = b2znw (is,ix,iy)
            ua_ext (ix,iy,is) = b2zfw (is,ix,iy)*
     .       sign(1.0_R8,bb(ix,iy,0))
            ta_ext (ix,iy,is) = b2ztw (is,ix,iy)
          end do
        end do
      end do

      !yamoto check whether the B2 receives na,ua,ta_ext properly or not. 
      !yamoto these summation have no physical meanings.

      write(*,*) 'b2ig: total impurity density in the system'
      write(*,*) 'b2ig: sum na_ext = ', sum(na_ext(:,:,:))
      write(*,*) 'b2ig: sum ua_ext = ', sum(ua_ext(:,:,:))
      write(*,*) 'b2ig: sum ta_ext = ', sum(ta_ext(:,:,:))

      write (*,*) 'b2ig: Leaving impgyro_main'

      end subroutine impgyro_main

      subroutine impgyro_last

      use b2mod_external

      integer :: fd

      if (.not. use_impgyro) return

      write (*,*) 'b2ig: Entering impgyro_last'

      last = .true.

#ifdef USE_MPI
      call mpi_bcast ( last , 1 , mpi_logical , 0 ,
     &     mpi_comm_world , mpi_err )
#else
#ifdef USE_IMPGYRO
      call mpi_bcast ( last , 1 , mpi_logical , 0 ,
     &     mpi_comm_world , mpi_err )
#endif
#endif

      call my_open_file ( fd , 'restart_ext' )
      write (fd) na_ext
      write (fd) ua_ext
      write (fd) ta_ext
      write (fd) sna_ext
      write (fd) sne_ext
      write (fd) smo_ext
      write (fd) she_ext
      write (fd) shi_ext
      write (fd) sch_ext
      close (fd)

      deallocate (ep_ir)
      deallocate (im_cc)
      deallocate (ie_cc)
      deallocate (ee_rd)
      deallocate (b2znw)
      deallocate (b2zfw)
      deallocate (b2ztw)

      write (*,*) 'b2ig: IMPGYRO restart files written'
      write (*,*) 'b2ig: Leaving impgyro_last'

      end subroutine impgyro_last

      subroutine my_open_file (fd ,fname)

      integer , intent(out) :: fd

      character*(*) ,intent(in) :: fname

      logical :: bopened

      do fd = 10 , 10000
        inquire (fd, opened = bopened )
        if (.not.bopened ) exit
      end do
      open (fd, file=fname, form='unformatted')

      end subroutine my_open_file

      end module b2mod_impgyro
