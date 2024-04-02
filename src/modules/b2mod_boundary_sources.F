      module b2mod_boundary_sources                          !srv 26.02.18 {
      use b2mod_types
      implicit none
      real (kind=R8), allocatable, save ::
     *  bv_po(:),
     *  bv_te(:),
     *  bv_ti(:),
     *  bv_ua(:,:),
     *  bv_na(:,:)

      contains

      subroutine alloc_b2mod_boundary_sources(nCv,nsd)
      implicit none
      integer nCv, nsd

      allocate(bv_po(nCv))
      allocate(bv_te(nCv))
      allocate(bv_ti(nCv))
      allocate(bv_ua(nCv,0:nsd-1))
      allocate(bv_na(nCv,0:nsd-1))

      bv_po = 0.0_R8
      bv_te = 0.0_R8
      bv_ti = 0.0_R8
      bv_ua = 0.0_R8
      bv_na = 0.0_R8

      write(*,*) 'Allocated space for b2mod_boundary_sources'
      return

      end subroutine alloc_b2mod_boundary_sources

      subroutine dealloc_b2mod_boundary_sources
      implicit none

      deallocate(bv_po)
      deallocate(bv_te)
      deallocate(bv_ti)
      deallocate(bv_ua)
      deallocate(bv_na)

      return
      end subroutine dealloc_b2mod_boundary_sources

      subroutine read_boundary_sources(nCv, ns, switch, itim)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_switches
      implicit none
      type (switches), intent (inout) :: switch
      integer nCv, ns
      integer is, iCv, i, itim
      logical b2news_solving(4)
      character*36 filename
      logical lexist
      character*4 chitim

      b2news_solving(1:4)=switch%no_solve.eq.0
      do i=1,4
         b2news_solving(i)=b2news_solving(i).or.
     &    (switch%no_solve.lt.0.and.
     &     mod(abs(switch%no_solve),2**i).ge.2**(i-1))
      enddo

      if (b2news_solving(2)) then
        if (switch%use_mms.eq.1) then
          filename='bv_na_us.dat'
          call find_file(filename,lexist)
          call xertst(lexist,'File '//trim(filename)//' is not found')
          open(70,file=trim(filename))
          read (70,"(1x,e23.16)") ((bv_na(iCv,is),iCv=1,nCv),is=0,ns-1)
          close(70)
        elseif (switch%use_mms.eq.2) then
          write (chitim,'(i4.4)') itim
          filename='bv_na_us'//chitim//'.dat'
          call find_file(filename,lexist)
          call xertst(lexist,'File '//trim(filename)//' is not found')
          open(70,file=trim(filename))
          read (70,"(1x,e23.16)") ((bv_na(iCv,is),iCv=1,nCv),is=0,ns-1)
          close(70)
        endif
      endif

      if (b2news_solving(1)) then
        filename='bv_ua_us.dat'
        call find_file(filename,lexist)
        call xertst(lexist,'File '//trim(filename)//' is not found')
        open(70,file=trim(filename))
        read (70,"(1x,e23.16)") ((bv_ua(iCv,is),iCv=1,nCv),is=0,ns-1)
        close(70)
      endif

      if (b2news_solving(4)) then
        filename='bv_te_us.dat'
        call find_file(filename,lexist)
        call xertst(lexist,'File '//trim(filename)//' is not found')
        open(70,file=trim(filename))
        read (70,"(1x,e23.16)") (bv_te(iCv),iCv=1,nCv)
        close(70)
        filename='bv_ti_us.dat'
        call find_file(filename,lexist)
        call xertst(lexist,'File '//trim(filename)//' is not found')
        open(70,file=trim(filename))
        read (70,"(1x,e23.16)") (bv_ti(iCv),iCv=1,nCv)
        close(70)
      endif

      if (b2news_solving(3)) then
        filename='bv_po_us.dat'
        call find_file(filename,lexist)
        call xertst(lexist,'File '//trim(filename)//' is not found')
        open(70,file=trim(filename))
        read (70,"(1x,e23.16)") (bv_po(iCv),iCv=1,nCv)
        close(70)
      endif

      return
      end subroutine read_boundary_sources

      end module b2mod_boundary_sources                      !srv 26.02.18 }
