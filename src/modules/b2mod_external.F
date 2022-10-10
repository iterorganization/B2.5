      module b2mod_external
      use b2mod_types

      integer, save :: ns_ext
      real (kind=R8), allocatable, save ::
     *  zn_ext(:),
     *  za_ext(:,:,:),
     *  za2_ext(:,:,:),
     *  pt_ext(:,:,:),
     *  na_ext(:,:,:),
     *  ni_ext(:,:,:),
     *  am_ext(:),
     *  ua_ext(:,:,:),
     *  ue_ext(:,:),
     *  ta_ext(:,:,:),
     *  ne_ext(:,:),
     *  fa_ext(:,:,:,:),
     *  ne2_ext(:,:),
     *  fhi_ext(:,:,:)
      real (kind=R8), allocatable, save ::
     *  sne_ext(:,:),
     *  sna_ext(:,:,:),
     *  smo_ext(:,:,:),
     *  she_ext(:,:),
     *  shi_ext(:,:),
     *  sch_ext(:,:)
      character*13, allocatable, save ::
     *  text_ext(:)
      logical, allocatable, save ::
     *  is_ext_neutral(:)
      logical, save :: initialised_externals = .false.
      logical, save :: allocated_externals = .false.

      contains

      subroutine alloc_b2mod_external(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd
      external xertst

*   ..test arguments
      call xertst (0.le.nxd.and.0.le.nyd, 'faulty argument nxd, nyd')
      call xertst (1.le.nsd, 'faulty argument nsd')

      if (allocated_externals) return
      ns_ext = 0
      allocate(ne_ext(-1:nxd,-1:nyd))
      allocate(ue_ext(-1:nxd,-1:nyd))
      allocate(ne2_ext(-1:nxd,-1:nyd))
      allocate(ni_ext(-1:nxd,-1:nyd,0:1))
      allocate(fhi_ext(-1:nxd,-1:nyd,0:1))
      ne_ext = 0.0_R8
      ni_ext = 0.0_R8
      ue_ext = 0.0_R8
      ne2_ext = 0.0_R8
      fhi_ext = 0.0_R8
      allocate(sne_ext(-1:nxd,-1:nyd))
      allocate(sna_ext(-1:nxd,-1:nyd,0:nsd-1))
      allocate(smo_ext(-1:nxd,-1:nyd,0:nsd-1))
      allocate(she_ext(-1:nxd,-1:nyd))
      allocate(shi_ext(-1:nxd,-1:nyd))
      allocate(sch_ext(-1:nxd,-1:nyd))
      sne_ext = 0.0_R8
      sna_ext = 0.0_R8
      smo_ext = 0.0_R8
      she_ext = 0.0_R8
      shi_ext = 0.0_R8
      sch_ext = 0.0_R8
      allocated_externals = .true.

      return
      end subroutine alloc_b2mod_external

      subroutine alloc_b2mod_external_dummy(nxd,nyd)
      implicit none
      integer nxd,nyd

      if (initialised_externals) return
cxpb  Allocate dummy arrays needed by b2xpfe
      allocate(za_ext(-1:nxd,-1:nyd,0))
      allocate(fa_ext(-1:nxd,-1:nyd,0:1,0))
      za_ext = 0.0_R8
      fa_ext = 0.0_R8
      initialised_externals = .true.

      return
      end subroutine alloc_b2mod_external_dummy

      subroutine dealloc_b2mod_external

      if (initialised_externals) then
        deallocate(za_ext)
        deallocate(fa_ext)
        if (ns_ext.gt.0) then
          deallocate(zn_ext)
          deallocate(za2_ext)
          deallocate(pt_ext)
          deallocate(na_ext)
          deallocate(ua_ext)
          deallocate(ta_ext)
          deallocate(am_ext)
          deallocate(text_ext)
          deallocate(is_ext_neutral)
        endif
        initialised_externals = .false.
      endif
      if (allocated_externals) then
        deallocate(ne_ext)
        deallocate(ue_ext)
        deallocate(ne2_ext)
        deallocate(ni_ext)
        deallocate(fhi_ext)
        deallocate(sne_ext)
        deallocate(sna_ext)
        deallocate(smo_ext)
        deallocate(she_ext)
        deallocate(shi_ext)
        deallocate(sch_ext)
        allocated_externals = .false.
      endif
      return
      end subroutine dealloc_b2mod_external

      subroutine dealloc_b2mod_external_dummy

      if (.not.initialised_externals) return
      deallocate(za_ext)
      deallocate(fa_ext)
      initialised_externals = .false.

      return
      end subroutine dealloc_b2mod_external_dummy

      subroutine read_b2mod_external(nxd,nyd,nsd,ismain)
      use b2mod_geo
      use b2mod_constants
      use b2mod_indirect
      use b2mod_elements
      implicit none
      integer nxd,nyd,nsd,ismain
      integer is,iss,is1,ix,iy,len
      integer nx,ny
      integer, save :: ikmn,ikmx
      real (kind=R8), allocatable :: pdum(:), pot(:), dext(:,:)
      character*256 filename
      logical toread
      logical open_file
      integer lnblnk
      external xertst, xerrab, open_file, find_file, extend, lnblnk

*   ..test arguments
      call xertst (0.le.nxd.and.0.le.nyd, 'faulty argument nxd, nyd')
      call xertst (1.le.nsd, 'faulty argument nsd')
      call xertst (0.le.ismain.and.ismain.lt.nsd,
     * 'faulty argument ismain')
cxpb
cxpb  First read ns_ext
cpxb
      filename='W.nut'
      call find_file(filename,toread)
      if (.not.toread) goto 974
      open(unit=74,file=filename)
      if (.not.initialised_externals) then
        read(74,*) nx, ny, ikmn, ikmx
        call xertst(nx.eq.nxd.and.ny.eq.nyd,
     .   'Dimension mismatch in W.nut!')
        ns_ext = ikmx - ikmn + 1
      else
        read(74,*)
      endif

      if (ns_ext.gt.0) then
        allocate(dext(0:nxd-1,0:nyd-1))
        if (.not.initialised_externals) then
          allocate(zn_ext(0:ns_ext-1))
          allocate(za_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(za2_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(pt_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(na_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(ua_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(ta_ext(-1:nxd,-1:nyd,0:ns_ext-1))
          allocate(fa_ext(-1:nxd,-1:nyd,0:1,0:ns_ext-1))
          allocate(am_ext(0:ns_ext-1))
          allocate(text_ext(0:ns_ext-1))
          allocate(is_ext_neutral(0:ns_ext-1))
          allocate(pot(0:ns_ext-1))
          pot(:)=0.0_R8
cxpb
cxpb  Here we read the external data
cxpb   Species data
cxpb    za_ext, pt_ext, am_ext, text_ext
          do is = 0, ns_ext-1
            toread = is.eq.0
            if (.not.toread)
     >       toread = nint(zn_ext(is)).ne.nint(zn_ext(is-1))
            if (toread) is1 = ikmn
            zn_ext(is) = 74.0_R8
            za_ext(:,:,is) = real(is1)
            am_ext(is) = 3.07832E-25_R8/mp
            text_ext(is) = ' '//elements(nint(zn_ext(is)))
            len=lnblnk(text_ext(is))
            if(is1.eq.0) then
              text_ext(is)(len+1:len+4)='^{0}'
            else if(is1.lt.10) then
              write(text_ext(is)(len+1:len+5),'(a,i1,a)') '^{+',is1,'}'
            else
              write(text_ext(is)(len+1:len+6),'(a,i2,a)') '^{+',is1,'}'
            endif
            if (maxval(za_ext(-1:nxd,-1:nyd,is)).eq.0.0_R8) then
              is_ext_neutral(is)=.true.
            else
              is_ext_neutral(is)=.false.
            endif
            za2_ext(-1:nxd,-1:nyd,is)=za_ext(-1:nxd,-1:nyd,is)**2
            if (toread) then
              if(allocated(pdum)) deallocate(pdum)
              allocate(pdum(1:nint(zn_ext(is))))
              if(open_file(99,'ionization_potentials',
     1         'data.local|modules/B2.5/Database')) then
                do iss=1,nint(zn_ext(is))-1
                  read(99,*)
                end do
                read(99,*) (pdum(iss),iss=1,nint(zn_ext(is)))
                do iss = 1, is1
                  pot(is) = pot(is) + pdum(iss)
                enddo
                pt_ext(:,:,is) = pot(is)
              else
                call xerrab('Please provide ionisation potentials!')
              endif
            else
              pot(is) = pot(is-1) + pdum(is1)
              pt_ext(:,:,is) = pot(is)
            endif
            is1 = is1 + 1
          end do
          deallocate(pot,pdum)
        endif

cxpb   Primary quantities
cxpb    na_ext, fa_ext, ta_ext, fhi_ext
cxpb
        na_ext = 0.0_R8
        ua_ext = 0.0_R8
        ta_ext = 0.0_R8
        fa_ext = 0.0_R8
        fhi_ext = 0.0_R8
        do is = 0, ns_ext-1
          read(74,*) dext
          na_ext(0:nxd-1,0:nyd-1,is) = dext(:,:) ! tungsten ik-th charge state density
          read(74,*) dext
          ua_ext(0:nxd-1,0:nyd-1,is) = sign(1.0_R8,bb(:,:,0))*dext(:,:) ! tungsten ik-th charge state velocity
          read(74,*) dext
          ta_ext(0:nxd-1,0:nyd-1,is) = dext(:,:) ! tungsten ik-th charge state temperature
        enddo
        read(74,*) dext(:,:)     ! electron particle source
     &                           !    by ionisation/recombination == ep_ir
        sne_ext(0:nxd-1,0:nyd-1) = dext(:,:)
        read(74,*) dext(:,:)     ! ion momentum source by Coulomb collision == im_cc
        smo_ext(0:nxd-1,0:nyd-1,ismain) = dext(:,:)
        read(74,*) dext(:,:)     ! ion heat source by Coulomb collision = ie_cc
        shi_ext(0:nxd-1,0:nyd-1) = dext(:,:)
        read(74,*) dext(:,:)     ! electron heat source by radiation = ee_rd
        she_ext(0:nxd-1,0:nyd-1) = dext(:,:)
        deallocate(dext)
        close(74)
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
        if (maxval(fa_ext).ne.0.0_R8.or.minval(fa_ext).ne.0.0_R8) then
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
     .            fa_ext(ix,iy,0,is)/pbs(ix,iy,0)/na_ext(ix,iy,is)
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
              endif
            enddo
          enddo
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
      else if (.not.initialised_externals) then
        call alloc_b2mod_external_dummy(nxd,nyd)
      endif

      initialised_externals=.true.
      return
  974 continue
      if (.not.initialised_externals)
     > call alloc_b2mod_external_dummy(nxd,nyd)
      initialised_externals=.true.
      return

      end subroutine read_b2mod_external

      end module b2mod_external
