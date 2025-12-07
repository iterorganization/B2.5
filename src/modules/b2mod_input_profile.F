      module b2mod_input_profile
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_time
      use b2mod_sources
     &     , only: ext_sna, ext_smo, ext_she, ext_shi, ext_sch, ext_sne
      use b2mod_balance !djm Jan2017
     &     , only: ext_sna0to1, ext_smo0to3, ext_she0to3, ext_shi0to3,
     &             balance_netcdf
      use b2mod_geo_corner
     &     , only: nmdpl
      use b2mod_subsys
#ifndef ALLOCATE_AND_NAMELIST
      use b2mod_dimensions
#endif
      implicit none
      integer nrr, nxx, nss
      integer, parameter :: nkind_source=6,nkind_data=2,nscale=10
      integer, parameter :: nkind_coeff=9
      integer, save :: csig_an_style=1
      integer :: ixref_profile, iyref_profile
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable, save :: nsdata(:,:,:)
      integer, allocatable, save :: nxdata(:,:,:)
      real(kind=R8), allocatable, save :: sdata(:,:,:,:)
      real(kind=R8), allocatable, save :: xdata(:,:,:,:)
      integer, allocatable, save :: ndata(:,:,:)
      real(kind=R8), allocatable, save :: tdata(:,:,:,:)
      integer, allocatable, save :: addspec(:,:,:)
#else
      parameter(nrr = DEF_NYD+2, nxx = DEF_NXD+2, nss = DEF_NSD)
      integer, save :: nsdata(nkind_data,nkind_source,0:nss)
      integer, save :: nxdata(nkind_data,nkind_source,0:nss)
      real(kind=R8), save :: sdata(2,nrr,nkind_source,0:nss)
      real(kind=R8), save :: xdata(2,nxx,nkind_source,0:nss)
      integer, save :: ndata(nkind_data,nkind_coeff,0:nss)
      real(kind=R8), save :: tdata(3,nrr,nkind_coeff,0:nss)
      integer, save :: addspec(nss,nkind_coeff,0:nss)
#endif

      contains

      subroutine alloc_input_profile(nx,ny,ns)
      implicit none
      integer nx, ny, ns

#ifdef ALLOCATE_AND_NAMELIST
      nrr = ny+2
      nxx = nx+2
      nss = ns

      if (.not.allocated(nsdata)) then
        allocate(nsdata(nkind_data,nkind_source,0:nss))
        allocate(nxdata(nkind_data,nkind_source,0:nss))
        allocate(sdata(2,nrr,nkind_source,0:nss))
        allocate(xdata(2,nxx,nkind_source,0:nss))
        allocate(ndata(nkind_data,nkind_coeff,0:nss))
        allocate(tdata(3,nrr,nkind_coeff,0:nss))
        allocate(addspec(nss,nkind_coeff,0:nss))
      endif
#else
      call xertst (nx.le.DEF_NXD, 'faulty input nx')
      call xertst (ny.le.DEF_NYD, 'faulty input ny')
      call xertst (ns.le.DEF_NSD, 'faulty input ns')
#endif

      return
      end subroutine alloc_input_profile

      subroutine dealloc_input_profile
      implicit none

#ifdef ALLOCATE_AND_NAMELIST
      if (allocated(nsdata)) then
        deallocate(nsdata)
        deallocate(nxdata)
        deallocate(sdata)
        deallocate(xdata)
        deallocate(ndata)
        deallocate(tdata)
        deallocate(addspec)
      endif
#endif

      return
      end subroutine dealloc_input_profile

      subroutine source_input(nx,ny,ns,
     &  sna0,smo0,she0,shi0,sch0,sne0,vol,nc)
!**************************************************************************
!*     subroutine for reading in source-profiles of density, parallel     *
!*     momentum and/or energy into B2.5 (written by Heimo Buerbaumer 1999)*
!*     Inputfile: b2.sources.profile                                      *
!**************************************************************************
      integer :: nx,ny,ns,nc
      integer :: ix,iy,i,ireg
      integer :: nanf,ndim,ifail,is,ndat,kind_data
      integer :: nxanf,nxend,nxdim,nxdat
      integer :: kind_source,spec,ifa
      real(kind=R8), save :: divheat
      real(kind=R8), save :: sources_time_mod,sources_time_switch
      character*256, save :: sources_filename
      real(kind=R8), dimension(-1:nx,-1:ny,0:3) :: she0,shi0,sch0
      real(kind=R8) :: sna0(-1:nx,-1:ny,0:1,0:ns-1)
      real(kind=R8) :: sne0(-1:nx,-1:ny,0:1)
      real(kind=R8) :: smo0(-1:nx,-1:ny,0:3,0:ns-1)
      real(kind=R8) :: vol(-1:nx,-1:ny),totpower,totna(0:ns-1)
      real(kind=R8) :: f(nrr),r(nrr),pr(nrr),pf(nrr),d(nrr)
      real(kind=R8) :: dhtot
      real(kind=R8) :: fx(nxx),x(nxx),xr(-1:nxx),xf(-1:nxx),dx(nxx)
      integer, save :: elm_count
      logical, save :: catch_up, new_files
#ifdef DBG
      logical, parameter :: debug = .true.
#else
      logical, parameter :: debug = .false.
#endif
      intrinsic int, mod

! initialise sources and data to zero
      if (nc.eq.0) then
        nsdata(:,:,:)=0
        sdata(:,:,:,:)=0.0_R8
        nxdata(:,:,:)=0
        xdata(:,:,:,:)=0.0_R8
        sna0(:,:,:,:)=0.0_R8
        sne0(:,:,:)=0.0_R8
        smo0(:,:,:,:)=0.0_R8
        she0(:,:,:)=0.0_R8
        shi0(:,:,:)=0.0_R8
        sch0(:,:,:)=0.0_R8
        divheat=0.0_R8
        sources_filename='b2.sources.profile'
        sources_time_mod=0.0_R8
        sources_time_switch=0.0_R8
        ixref_profile=-2
        iyref_profile=-2
! read profile
        new_files = .true.
        call read_source_profile
      endif

      if(nc.eq.0) then
        if(sources_time_mod.gt.0.0_R8) then
          elm_count=int(tim/sources_time_mod)
        else
          elm_count=0
        endif
        if(sources_time_mod.gt.0.0_R8) then
          catch_up=mod(tim,sources_time_mod).ge.sources_time_switch.and.
     &             sources_time_switch.gt.0.0_R8
        else if (sources_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.sources_time_switch
        else
          catch_up=.false.
        end if
        do while(catch_up)
          write(*,*) 'sources: catching up...'
          call read_source_profile
          if(sources_time_mod.gt.0.0_R8) then
            catch_up=mod(tim,sources_time_mod).ge.sources_time_switch
     &        .and.sources_time_switch.gt.0.0_R8
          else if (sources_time_switch.gt.0.0_R8) then
            catch_up=tim.ge.sources_time_switch
          else
            catch_up=.false.
          endif
        enddo
      else
        if(sources_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3e12.4,2i6)') 'SOURCES: ',
     .     tim,mod(tim,sources_time_mod),
     1     sources_time_switch,int(tim/sources_time_mod),elm_count
          catch_up=
     1    (mod(tim,sources_time_mod).ge.sources_time_switch.and.
     1     sources_time_switch.gt.0.0_R8).or.
     2     int(tim/sources_time_mod).gt.elm_count
        else if (sources_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.sources_time_switch
          if (catch_up) write(*,'(a,a)')
     .     'Switching to next sources profile datafile ',
     .      trim(sources_filename)
        else
          catch_up=.false.
        end if
        if(catch_up) then
          if(sources_time_mod.gt.0.0_R8)
     >     elm_count=int(tim/sources_time_mod)
          call read_source_profile
          new_files = .true.
        endif
      endif

      if (new_files) then
        if (maxval(nxdata).gt.0) open (96,file='dx2')
        if (maxval(nxdata).gt.0) open (95,file='dx1')
        if (maxval(nsdata).gt.0) open (94,file='da1')
        if (maxval(nsdata).gt.0) open (93,file='da2')
      end if

      do kind_source=1,nkind_source
        do spec=0,ns-1
          kind_data=0
          ndat=0
          do i=1,nkind_data
            if (nsdata(i,kind_source,spec).ne.0.and.kind_data.eq.0) then
              kind_data=i
              ndat=nsdata(kind_data,kind_source,spec)
            else if (nsdata(i,kind_source,spec).ne.0) then
              write(*,*) 'Conflicting data types:'
              write(*,*) 'Both KIND_DATA = 1 and 2 are declared for '
              write(*,*) 'source poloidal profile. '
              write(*,*) 'KIND_DATA = 2 ignored.'
              write(*,'(a,i2,a,i3)')
     .         'KIND_SOURCE = ', kind_source, ' Species = ' ,spec
            endif
          enddo
          if (ndat.eq.0) goto 1
          call xertst(ndat.le.nrr,
     .     'Increase NYD in b2mod_dimensions to match size of profile')
!       allocate (f(1:ndat))
!       allocate (r(1:ndat))
          do i=1,ndat
            r(i)=sdata(1,i,kind_source,spec)
            f(i)=sdata(2,i,kind_source,spec)
            if (new_files) write(93,*) r(i),f(i)
          enddo
!          write(*,*) 'vor sources_profile'
          call sources_profile
     .     (kind_data,ndat,r,ndim,pr,nx,ny,ifa,nanf)
!          write (*,*) 'ifa=',ifa
          if (ifa.eq.1) goto 99
!          allocate (pf(1:ndim))
          call E01BEF(ndat,r,f,d,ifail)
!          write(*,*) 'zwischen,ifail=',ifail
          call E01BFF(ndat,r,f,d,ndim,pr,pf,ifail)
!       obtain axial profile
          nxdat=0
          kind_data=0
          do i=1,nkind_data
            if (nxdata(i,kind_source,spec).ne.0.and.kind_data.eq.0) then
              kind_data=i
              nxdat=nxdata(kind_data,kind_source,spec)
            else if (nxdata(i,kind_source,spec).ne.0) then
              write(*,*) 'Conflicting data types:'
              write(*,*) 'Both KIND_DATA = 1 and 2 are declared for '
              write(*,*) 'source axial profile. '
              write(*,*) 'KIND_DATA = 1 ignored.'
              write(*,'(a,i2,a,i3)')
     .         'KIND_SOURCE = ', kind_source, ' Species = ' ,spec
            endif
          enddo
          xf(-1)=0.0_R8
          xf(0:nx-1)=1.0_R8
          xf(nx)=0.0_R8
          if (nxdat.gt.0) then
            call xertst(nxdat.le.nxx+2,
     .     'Increase NXD in b2mod_dimensions to match size of profile')
            do i=1,nxdat
              x(i)=xdata(1,i,kind_source,spec)
              fx(i)=xdata(2,i,kind_source,spec)
              if (new_files) write(95,*) x(i),fx(i)
            enddo
!            write(*,*) 'vor axial_profile'
            call axial_profile
     .       (kind_data,nxdat,x,nxdim,xr,nx,ny,ifa,nxanf,nxend)
            if (ifa.eq.1) goto 98
            call E01BEF(nxdat,x,fx,dx,ifail)
            call E01BFF(nxdat,x,fx,dx,nxdim,
     .       xr(nxanf:nxend),xf(nxanf:nxend),ifail)
          endif

          do iy=1,ndim
            do ix=-1,nx
              ireg=region(ix,nanf-3+iy,0)
              if (bottomiy(ix,nanf-3+iy).ne.-2 .and.
     &            topiy(ix,nanf-3+iy).ne.ny+1 .and.
     &            leftix(ix,nanf-3+iy).ne.-2 .and.
     &            rightix(ix,nanf-3+iy).ne.nx+1 .and.
     &            ireg.ne.0 .and.
     &            (mod(ireg,4).eq.1 .or. mod(ireg,4).eq.2)) then
                if (kind_source.eq.1) then
                  sna0(ix,nanf-3+iy,0,spec)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                elseif (kind_source.eq.2) then
                  smo0(ix,nanf-3+iy,0,spec)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                elseif (kind_source.eq.3) then
                  she0(ix,nanf-3+iy,0)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                elseif (kind_source.eq.4) then
                  shi0(ix,nanf-3+iy,0)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                elseif (kind_source.eq.5) then
                  sch0(ix,nanf-3+iy,0)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                elseif (kind_source.eq.6) then
                  sne0(ix,nanf-3+iy,0)=
     .             pf(iy)*xf(ix)*vol(ix,nanf-3+iy)
                endif
              endif
              if ((nc.eq.0.and.iy.eq.1).or.debug) then
                if (nxdat.gt.0) then
                  if (kind_source.eq.1) write(*,'(a,i2,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: sna ',spec,xr(ix),xf(ix)
                  if (kind_source.eq.2) write(*,'(a,i2,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: smo ',spec,xr(ix),xf(ix)
                  if (kind_source.eq.3) write(*,'(a,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: she ',xr(ix),xf(ix)
                  if (kind_source.eq.4) write(*,'(a,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: shi ',xr(ix),xf(ix)
                  if (kind_source.eq.5) write(*,'(a,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: sch ',xr(ix),xf(ix)
                  if (kind_source.eq.6) write(*,'(a,1p,2(1x,e14.7))')
     &             'AXIAL-SOURCES: sne ',xr(ix),xf(ix)
                  if (new_files) write (96,*) xr(ix),xf(ix)
                endif
              endif
            enddo
            if (nc.eq.0.or.debug) then
              if (ndat.gt.0) then
                if (kind_source.eq.1) write(*,'(a,i2,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: sna ',spec,pr(iy),pf(iy)
                if (kind_source.eq.2) write(*,'(a,i2,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: smo ',spec,pr(iy),pf(iy)
                if (kind_source.eq.3) write(*,'(a,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: she ',pr(iy),pf(iy)
                if (kind_source.eq.4) write(*,'(a,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: shi ',pr(iy),pf(iy)
                if (kind_source.eq.5) write(*,'(a,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: sch ',pr(iy),pf(iy)
                if (kind_source.eq.6) write(*,'(a,1p,2(1x,e14.7))')
     &           'PROFILE-SOURCES: sne ',pr(iy),pf(iy)
                if (new_files) write (94,*) pr(iy),pf(iy)
              endif
            endif
          enddo
    1     continue
        enddo
      enddo
! Calculation of the total power input
      if (nc.eq.0 .or. debug) then
        totpower=0.0_R8
        totna(:)=0.0_R8
        do iy=-1,ny
          do ix=-1,nx
            totpower=totpower+she0(ix,iy,0)+shi0(ix,iy,0)
            do is=0,ns-1
              totna(is)=totna(is)+sna0(ix,iy,0,is)
            enddo
          enddo
        enddo
        write (*,*) 'Total profile-source Power=', totpower
        do is=0,ns-1
          write(*,*) 'Total profile-source na:', is, totna(is)
        enddo
      endif

      if (divheat.GT.0.0_R8) then
        dhtot=0.0_R8
        do iy=-1,ny
          do ix=-1,nx
            if((mod(region(ix,iy,0),4).eq.3 .or.
     &         (mod(region(ix,iy,0),4).eq.0 .and. region(ix,iy,0).ne.0))
     &        .and. leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1
     &        .and. bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     >       then
              shi0(ix,iy,0)=shi0(ix,iy,0)+divheat*vol(ix,iy)
              dhtot=dhtot+divheat*vol(ix,iy)
            endif
          enddo
        enddo
        write (*,*) 'DIVHEAT=', dhtot,divheat
      endif

      goto 100
   98 write (*,*) 'Error in computing axial profile'
      goto 100
   99 write (*,*) 'Error in reading the dsp file'
  100 continue

*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smo0(-1,-1,0,is))
      enddo
      call b2xzdd (nx, ny, 1, sne0)
      call b2xzdd (nx, ny, 3, she0)
      call b2xzdd (nx, ny, 3, shi0)
      call b2xzdd (nx, ny, 3, sch0)

      ext_sna(:,:,:)=sna0(:,:,0,:)
      ext_smo(:,:,:)=smo0(:,:,0,:)
      ext_sne(:,:)=sne0(:,:,0)
      ext_she(:,:)=she0(:,:,0)
      ext_shi(:,:)=shi0(:,:,0)
      ext_sch(:,:)=sch0(:,:,0)

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        ext_sna0to1 = sna0
        ext_smo0to3 = smo0
        ext_she0to3 = she0
        ext_shi0to3 = shi0
      endif

      if (new_files) then
        if (maxval(nsdata).gt.0) close(93)
        if (maxval(nsdata).gt.0) close(94)
        if (maxval(nxdata).gt.0) close(95)
        if (maxval(nxdata).gt.0) close(96)
        new_files = .false.
      end if

      return

      contains

      subroutine read_source_profile

      character*260 filename
      logical file_ok

      namelist/profile/nsdata,nxdata,sdata,xdata,divheat,
     . ixref_profile,iyref_profile,
     . sources_time_mod,sources_time_switch,sources_filename

!  read the profile data
      filename=sources_filename
      sources_time_mod = 0.0_R8
      sources_time_switch = 0.0_R8
      call find_file(filename,file_ok)
      if(file_ok) then
        open(92,file=filename)
        write(*,*)
     .   'Reading sources profile namelist from ', trim(filename)
#ifdef USE_MPI
        read(92,profile,err=999)
#else
        read(92,profile)
#endif
        write(*,*)
     .   'Read sources profile namelist from ',trim(filename)
        close(92)
      else
        call xerrab ('No '//trim(filename))
      endif
!     write(*,profile)
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading PROFILE namelist!')
#endif
      end subroutine read_source_profile

*-----------------------------------------------------------------------

      subroutine sources_profile(kind_data,ndat,r,ndim,pr,nx,ny,ifa,
     & nfitanf)

      use b2mod_indirect
      use b2mod_geo
      use b2mod_mwti
     & , only: output_ds
      implicit none
      integer:: i,nx,ny,ifa
      integer:: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: pra(1:ny+2),prdsa(1:ny+2),r(nrr),pr(nrr),prtmp
      real (kind=R8) :: dsp(-1:ny)
      integer jxi, jxa, jsep, iyastrt, iyaend, ncall, nya
      integer ixref, iyref
      character*256 filename
      logical exist
      external ipgetr
      external get_jsep
      save jxi, jxa, jsep, ixref, iyref, iyastrt, iyaend, ncall
      data ncall/0/

      call subini('sources_profile')
      if (ncall.eq.0) then
        call ipgeti ('b2trno_csig_an_style', csig_an_style)
        call get_jsep(nx,ny,jxi,jxa,jsep)
        call xertst (-1.le.jsep.and.jsep.le.ny,
     2   'faulty parameter jsep')
        call ipgeti ('b2mwti_jxa',jxa)
        call xertst (-1.le.jxa.and.jxa.le.nx,
     2   'faulty internal parameter jxa')
      end if
      if (new_files) then
        if (ixref_profile.ne.-2) then
          ixref = ixref_profile
        else
          ixref = jxa
          call ipgeti ('set_transport_ixref', ixref)
        end if
        if (iyref_profile.ne.-2) then
          iyref = iyref_profile
        else
          iyref = jsep
          call ipgeti ('set_transport_iyref', iyref)
        end if
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument for sources profile ixref')
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument for sources profile iyref')
        call output_ds(ny,ixref,0,iyref,iyastrt,iyaend,'dsp',dsp)
      endif

      filename='dsp'
      call find_file(filename,exist)
      if (exist) then
        open(unit=91,file=filename)
      else
        ifa=1
        call subend()
        return
      end if
      i = 0
      nya = 0
   91 continue
      i = i + 1
      read (91,*,END=92,ERR=97) prtmp
      prdsa(i) = prtmp
      nya = i
      if (i.lt.ny+2) goto 91
   92 continue
      nfitanf=1
      IFA=0
      do i=1,nya
! input in metres from separatrix
        if (kind_data.eq.1) pra(i)=prdsa(i)
! input in poloidal rho
        if (kind_data.eq.2) then
          write(*,*) 'Input profile in poloidal rho not yet coded!'
          goto 97
        endif

        if (pra(i).LT.r(1)) nfitanf=nfitanf+1
        nfitend=i-1
        if (pra(i).GT.r(ndat)) exit
        nfitend=i
      enddo

      if (nfitend.LE.nfitanf) goto 97

      ndim=nfitend+1-nfitanf
!     allocate (pr(1:ndim))

      do i=1,ndim
        pr(i)=pra(nfitanf-1+i)
      enddo
      goto 108
   97 ifa=1

  108 continue
      close (91)
      ncall=ncall+1
      call subend ()
      return
      end subroutine sources_profile

*-----------------------------------------------------------------------

      subroutine axial_profile(kind_data,ndat,x,ndim,pr,nx,ny,ifa,
     & nfitanf,nfitend)

      use b2mod_indirect
      use b2mod_geo
      implicit none
      integer :: i,nx,ny,ifa
      integer :: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: x(1:ndat),pr(-1:nxx-2)
      real(kind=R8) :: l_conn(-1:nx), l_total, l_ref
      integer jxi, jxa, jsep, ncall
      integer ixref, iyref
      integer ireg
      external ipgetr
      external get_jsep
      save jxi, jxa, jsep, ixref, iyref, ncall
      data ncall/0/

      call subini('axial_profile')
      if (ncall.eq.0) then
        call get_jsep(nx,ny,jxi,jxa,jsep)
        call ipgeti ('b2mwti_jxa',jxa)
        call xertst (-1.le.jxa.and.jxa.le.nx,
     2   'faulty internal parameter jxa')
        ixref = jxa
        iyref = jsep
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
      endif

      if (kind_data.eq.1) then
! Establish relative profile in relative connection length
! For closed field lines, zero point sits at ix = ixref
! Reference connection length computed along iy = iyref
        ireg = mod(region(ixref,iyref,0),4)
        l_total = 0.0_R8
        pr(-1) = 0.0_R8
        l_conn(-1) = 0.0_R8
        do i = 0, nx-1
          l_conn(i) = hx(i,iyref)/
     1         (sqrt(bb(i,iyref,0)**2+bb(i,iyref,1)**2)/bb(i,iyref,3))
          pr(i) = pr(i-1) + 0.5_R8*(l_conn(i)+l_conn(i-1))
          l_total = l_total + l_conn(i)
        end do
        pr(nx) = pr(nx-1) + 0.5_R8*l_conn(nx-1)
        l_ref = pr(ixref)
        if (l_total.le.0.0_R8) goto 97
        if (ireg.eq.1 .and. nnreg(0).ge.4) then
          do i = -1, nx
            if (mod(region(i,iyref,0),4).eq.ireg) then
              pr(i) = (pr(i) - l_ref) / l_total
              if (pr(i).lt.0.0_R8) pr(i) = pr(i) + 1.0_R8
            else
              pr(i) = 0.0_R8
            end if
          end do
        else
          do i = -1, nx
            pr(i) = pr(i) / l_total
          end do
        endif
      else if (kind_data.eq.2) then
! Establish profile in computational space
        pr(-1) = 0.0_R8
        do i = 0, nx-1
          pr(i) = (i + 0.5_R8) / (real(nx))
        end do
        pr(nx) = 1.0_R8
      endif
      ifa = 0

      nfitanf = -1
      do i = 0, nx
        if (pr(i).lt.x(1)) nfitanf = nfitanf+1
        nfitend=i-1
        if (pr(i).gt.x(ndat)) exit
        nfitend=i
      end do
!     write (*,*) nfitanf,nfitend
      if (nfitend.LE.nfitanf) goto 97

      ndim=nfitend+1-nfitanf
!     write (*,*) 'ndim=',ndim
!     allocate (pr(1:ndim))

      goto 108
   97 ifa=1

  108 continue
      ncall=ncall+1
      call subend ()
      return
      end subroutine axial_profile

      end subroutine source_input

*-----------------------------------------------------------------------

      subroutine transport_input(nx, ny, ns, rza, na, ne, te, ti, hcib,
     &  dna0, dpa0, vla0, vsa0, hci0, hce0, sig0, alf0, ncall)

!**************************************************************************
!*     subroutine for reading in experimental transport coefficients from *
!*     the midplane and calculating transport coefficients                *
!*     for B2.5 (written by Heimo Buerbaumer 1999)                        *
!*     Inputfile: b2.transport.inputfile                                  *
!**************************************************************************
      use b2mod_time
      use b2mod_constants
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_geo
      IMPLICIT NONE
      integer :: nx,ny,ns
      integer :: iy,ix,i,ifail,spec,nspec,iiy,iir,
     &      nfitanf,ndim,is,nsp(1:ns+1),iloop,ncall,ifa,ndat,
     &      kind_coeff,kind_data,isc,idy
      integer :: ixwidth(nscale), ixmean(nscale)

      real(kind=R8), dimension(-1:nx,-1:ny,0:ns-1) ::
     &    rza,dna0,dpa0,vsa0,na
      real(kind=R8), dimension(-1:nx,-1:ny) ::
     &    hci0,hce0,sig0,alf0,ne,te,ti

      real(kind=R8), save :: scaling_width(nscale)
      real(kind=R8), save :: scaling_strength(nscale)
      real(kind=R8) :: f(nrr),r(nrr),pr(nrr),pf(nrr),d(nrr),
     &    pelm(nrr),felm(nrr),delm(nrr),prof(-1:nx,nrr),
     &    vla0(-1:nx,-1:ny,0:1,0:ns-1),hcib(-1:nx,-1:ny,0:ns-1),
     &    pelm_rescaled(-1:nx,nrr)
      character*256, save :: transport_filename
      real (kind=R8), save :: transport_time_mod
      real (kind=R8), save :: transport_time_switch
      real (kind=R8), save :: elm_time_begin(nkind_coeff), 
     1 elm_time_plateau(nkind_coeff), elm_time_recovery(nkind_coeff),
     2 elm_time_end(nkind_coeff), elm_crash_factor(nkind_coeff),
     3 elm_recovery_factor(nkind_coeff), elm_ballooning_exp(nkind_coeff)
      real (kind=R8), save :: elm_time_period
      real (kind=R8) :: frac_elm
      real (kind=R8) :: poloidal_scale_factor
      real (kind=R8) :: elm_poloidal_scaling_factor(-1:nx,nkind_coeff)
      real (kind=R8) :: thet, ilp(-1:ny), ilpb, hxface, hxpface, hxnface
      real (kind=R8) :: theta(-1:nx,-1:ny)
      integer, save :: scaling_ix_begin(nscale), scaling_ix_end(nscale)
      integer, save :: elm_ix_begin, elm_ix_end, elm_ixref
      logical, save :: region_flags(NNREGMAX,nkind_coeff), no_pflux,
     1 no_div, poloidal_scaling(nscale), elm_crash_time(nkind_coeff),
     2 elm_plateau_time(nkind_coeff), elm_recovery_time(nkind_coeff),
     3 elm_crash_region, elm_plateau_region, elm_recovery_region,
     4 elm_data, elm_dynamics_multiple
      character*260 filename
      character*26 hlp_frm
      integer, save :: elm_count
      logical catch_up, lpflux, ldiv
      logical, save :: new_files
      intrinsic sqrt
      namelist /transport/ ndata, tdata, addspec,
     1   transport_filename, transport_time_mod, transport_time_switch,
     2   scaling_ix_begin, scaling_ix_end, scaling_strength,
     3   region_flags, poloidal_scaling, scaling_width,
     4   elm_time_begin, elm_time_plateau, elm_time_recovery,
     5   elm_time_end, elm_time_period, elm_ix_begin, elm_ix_end, 
     6   elm_ixref, elm_crash_factor, elm_recovery_factor,
     7   elm_ballooning_exp, elm_dynamics_multiple, no_pflux, no_div

      data transport_filename/'b2.transport.inputfile'/,
     1 transport_time_mod/0.0_R8/, transport_time_switch/0.0_R8/

      if(ncall.eq.0) then
        call xertst(NNREGMAX.ge.nnreg(0),
     1   'Increase NNREGMAX in b2mod_numerics_namelist')
        filename=transport_filename
        new_files = .true.
        call read_transport_input
!  read the profile data
        if(transport_time_mod.gt.0.0_R8) then
          elm_count=int(tim/transport_time_mod)
        else
          elm_count=0
        endif
        if(transport_time_mod.gt.0.0_R8) then
          catch_up=mod(tim,transport_time_mod).ge.transport_time_switch
     &      .and.transport_time_switch.gt.0.0_R8
        else if (transport_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.transport_time_switch
        else
          catch_up=.false.
        end if
        do while(catch_up)
          write(*,*) 'transport_input: catching up ...'
          filename=transport_filename
          call read_transport_input
          if(transport_time_mod.gt.0.0_R8) then
            catch_up=transport_time_switch.gt.0.0_R8.and.
     .       mod(tim,transport_time_mod).ge.transport_time_switch
          else if (transport_time_switch.gt.0.0_R8) then
            catch_up=tim.ge.transport_time_switch
          else
            catch_up=.false.
          endif
        enddo
      else
        if(transport_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3(1x,e14.7),2i4)')
     .     'ELM: ',tim,mod(tim,transport_time_mod),
     1     transport_time_switch,int(tim/transport_time_mod),elm_count
          catch_up=
     1    (mod(tim,transport_time_mod).ge.transport_time_switch.and.
     1     transport_time_switch.gt.0.0_R8).or.
     2     int(tim/transport_time_mod).gt.elm_count
        else if (transport_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.transport_time_switch
          if (catch_up) write(*,'(a,a)')
     .     'Switching to next transport namelist file ',
     .      trim(transport_filename)
        else
          catch_up=.false.
        endif
        if(catch_up) then
          if(transport_time_mod.gt.0.0_R8)
     >     elm_count=int(tim/transport_time_mod)
          filename=transport_filename
          new_files = .true.
          call read_transport_input
        endif
      endif

      if (new_files) then
        open (94,file='da1')
        open (93,file='da2')
      end if

      do isc=1,nscale
        ixwidth(isc) = scaling_ix_end(isc) - scaling_ix_begin(isc) + 1
        ixmean(isc) = (scaling_ix_end(isc) + scaling_ix_begin(isc)) / 2
        if (scaling_width(isc).eq.0.0_R8 .and. ixwidth(isc).gt.1)
     >   scaling_width(isc) = 0.5_R8*
     .      sqrt((ixwidth(isc)+1.0_R8)*(ixwidth(isc)-1.0_R8)/3.0_R8)
      enddo
      if (elm_time_period.gt.0.) then
        if (elm_dynamics_multiple) then
          do kind_coeff=1,nkind_coeff
            if (elm_time_plateau(kind_coeff).gt.0. .and. 
     1              elm_time_recovery(kind_coeff).gt.0.) then
              if (elm_time_begin(kind_coeff).eq.
     1                elm_time_plateau(kind_coeff)) then
                elm_crash_time(kind_coeff)=.false.
              else
                elm_crash_time(kind_coeff)=mod(tim,elm_time_period).ge.
     &                  elm_time_begin(kind_coeff).and.
     &                  mod(tim,elm_time_period).le.
     &                  elm_time_plateau(kind_coeff)
              endif
              elm_plateau_time(kind_coeff)=mod(tim,elm_time_period).ge.
     &                  elm_time_plateau(kind_coeff).and.
     &                  mod(tim,elm_time_period).le.
     &                  elm_time_recovery(kind_coeff)
              if (elm_time_recovery(kind_coeff).eq.
     1                elm_time_end(kind_coeff)) then
                elm_recovery_time(kind_coeff)=.false.
              else
                elm_recovery_time(kind_coeff)=mod(tim,elm_time_period)
     &                 .ge.elm_time_recovery(kind_coeff).and.
     &                 mod(tim,elm_time_period).le.
     &                 elm_time_end(kind_coeff)
              endif
            else
              elm_crash_time(kind_coeff)=.false.
              elm_plateau_time(kind_coeff)=mod(tim,elm_time_period).ge.
     &                  elm_time_begin(kind_coeff)
     &                 .and.mod(tim,elm_time_period)
     &                 .le.elm_time_end(kind_coeff)
              elm_recovery_time(kind_coeff)=.false.
            endif
          enddo
        else
          if (elm_time_plateau(1).gt.0. .and. 
     1            elm_time_recovery(1).gt.0.) then
            if (elm_time_begin(1).eq.elm_time_plateau(1)) then
              elm_crash_time(1)=.false.
            else
              elm_crash_time(1)=mod(tim,elm_time_period)
     &               .ge.elm_time_begin(1).and.
     &               mod(tim,elm_time_period).le.elm_time_plateau(1)
            endif
            elm_plateau_time(1)=mod(tim,elm_time_period)
     &               .ge.elm_time_plateau(1).and.
     &               mod(tim,elm_time_period).le.elm_time_recovery(1)
            if (elm_time_recovery(1).eq.elm_time_end(1)) then
              elm_recovery_time(1)=.false.
            else
              elm_recovery_time(1)=mod(tim,elm_time_period).ge.
     &               elm_time_recovery(1).and.
     &               mod(tim,elm_time_period).le.elm_time_end(1)
            endif
          else
            elm_crash_time(1)=.false.
            elm_plateau_time(1)=mod(tim,elm_time_period).ge.
     &               elm_time_begin(1).and.mod(tim,elm_time_period)
     &               .le.elm_time_end(1)
            elm_recovery_time(1)=.false.
          endif
        endif
      else
        if (elm_dynamics_multiple) then
          do kind_coeff=1,nkind_coeff
            elm_crash_time(kind_coeff)=.false.
            elm_plateau_time(kind_coeff)=.false.
            elm_recovery_time(kind_coeff)=.false.
          enddo
        else
          elm_crash_time(1)=.false.
          elm_plateau_time(1)=.false.
          elm_recovery_time(1)=.false.
        endif
      endif
      elm_data = maxval(tdata(3,:,:,:)).gt.0.0_R8

c     Calculate poloidal angle theta (same code used in b2tqna.F)
      do iy=-1,ny
        do ix=-1,nx
          theta(ix,iy)=-1.0_R8
        enddo
      enddo
      do iy=0,ny
        thet=0._R8
        ilp(iy)=0._R8
        do ix=elm_ixref,-1,-1
           if (on_closed_surface(ix,iy)) then
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
        enddo
        do ix=nx,elm_ixref+1,-1
           if (on_closed_surface(ix,iy)) then
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
        enddo
        theta(elm_ixref,iy)=thet
        ilpb=thet
        do ix=-1,nx
           if (on_closed_surface(ix,iy)) then
            theta(ix,iy)=2.0_R8*pi*theta(ix,iy)/ilpb
           endif
        enddo
      enddo

!        write (*,*) 'vor loop'
      iloop=0
      do kind_coeff=1,nkind_coeff
        do spec=0,ns-1
          do kind_data=1,nkind_data
!           write (*,*) 'loop-begin', kind_coeff,spec,kind_data
            ndat=ndata(kind_data,kind_coeff,spec)
            if (ndat.eq.0) goto 1
            call xertst(ndat.le.nrr,
     .      'Increase NYD in b2mod_dimensions to match size of profile')
            nspec=1
            nsp(1)=spec
            do is=1,ns
              if (addspec(is,kind_coeff,spec).GE.0) then
                nsp(nspec+1)=addspec(is,kind_coeff,spec)
                nspec=nspec+1
              endif
            enddo
            do i=1,ndat
              r(i)=tdata(1,i,kind_coeff,spec)
              f(i)=tdata(2,i,kind_coeff,spec)
              if (elm_data) felm(i)=tdata(3,i,kind_coeff,spec)
!             if (new_files) write(93,*) r(i),f(i),felm(i)
            enddo

            call transport_profile(1,ndat,r,ndim,pr,nx,ny,ifa,nfitanf)
!           if (iloop.eq.0) print *,'R=',r(1:ndat)
!           if (iloop.eq.0) print *,'F=',f(1:ndat)
!           if (iloop.eq.0) print *,'PR=',pr(1:ndim)
!           write (*,*) 'ifa=',ifa
            if (ifa.eq.1) goto 99

            call E01BEF(ndat,r,f,d,ifail)
            if (elm_data) call E01BEF(ndat,r,felm,delm,ifail)
!           write(*,*) 'zwischen,ifail=',ifail
            call E01BFF(ndat,r,f,d,ndim,pr,pf,ifail)
!           if (iloop.eq.0) print *,'PF=',pf(1:ndim)
            if (elm_data)
     >      call E01BFF(ndat,r,felm,delm,ndim,pr,pelm,ifail)

            if (new_files) then

              if (elm_data) then
                do iy=1,ndat
                  write (93,*) r(iy),f(iy),felm(iy)
                enddo

                do iy=1,ndim
                  write (94,*) pr(iy),pf(iy),pelm(iy)
                enddo
              else
                do iy=1,ndat
                  write (93,*) r(iy),f(iy)
                enddo

                do iy=1,ndim
                  write (94,*) pr(iy),pf(iy)
                enddo
              endif
            endif

            do is=1,nspec
              do ix=-1,nx
                do iy=1,ndim

                  iiy=nfitanf-3+iy
                  iir=region(ix,iiy,0)
                  if(iir.lt.0.or.iir.gt.NNREGMAX) then
                    write(*,'(2a,3i4)')
     &               'transport_input: out of bounds ',
     &               'ix,iy,region(ix,iy) ',ix,iiy,iir
                    iir=NNREGMAX
                  endif
                  lpflux=iir.eq.3.and.iiy.le.topcut(1).and.nnreg(0).ge.4
                  lpflux=lpflux.or.
     &             (iir.eq.4.and.iiy.le.topcut(1).and.
     &             (nnreg(0).eq.4.or.nnreg(0).eq.5.or.nnreg(0).eq.7))
                  lpflux=lpflux.or.
     &             (iir.eq.4.and.iiy.le.topcut(2).and.nnreg(0).eq.8)
                  lpflux=lpflux.or.
     &             (iir.eq.5.and.iiy.le.topcut(2).and.nnreg(0).eq.7)
                  lpflux=lpflux.or.(iir.eq.6.and.nnreg(0).eq.7)
                  lpflux=lpflux.or.(iir.eq.7.and.nnreg(0).eq.7.and.
     &             (iiy.lt.min(topcut(1),topcut(2)).or.
     &              iiy.ge.max(topcut(1),topcut(2))))
                  lpflux=lpflux.or.
     &             (iir.eq.7.and.iiy.le.topcut(2).and.nnreg(0).eq.8)
                  lpflux=lpflux.or.(iir.eq.8.and.iiy.le.topcut(1))
                  ldiv=iir.ne.0.and.mod(iir-1,4)+1.ge.3
                  ldiv=ldiv.or.
     &             (nnreg(0).eq.7.and.(iir.eq.5.or.iir.eq.6))
                  poloidal_scale_factor = 1.0_R8
                  if (iir.gt.0) then
                   do isc=1,nscale
                    if (poloidal_scaling(isc) .and.
     $                ix.ge.scaling_ix_begin(isc) .and.
     $                    ix.le.scaling_ix_end(isc)) then
                      poloidal_scale_factor = poloidal_scale_factor*
     $                    (1.0_R8 + scaling_strength(isc)/
     $                    (sqrt(2.0_R8*pi)*scaling_width(isc))*
     $                    exp(-(ix-ixmean(isc))*(ix-ixmean(isc))/
     $                         (2.0_R8*scaling_width(isc)**2)))
                    endif
                   enddo

                   if (elm_dynamics_multiple) then
                     elm_crash_region=
     &                      elm_crash_time(kind_coeff).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                     elm_plateau_region=
     &                      elm_plateau_time(kind_coeff).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                     elm_recovery_region=
     &                      elm_recovery_time(kind_coeff).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                   else
                     elm_crash_region=
     &                      elm_crash_time(1).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                     elm_plateau_region=
     &                      elm_plateau_time(1).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                     elm_recovery_region=
     &                      elm_recovery_time(1).and.elm_data
     &                     .and.ix.ge.elm_ix_begin.and.ix.le.elm_ix_end
                   endif

                   if (theta(ix,0).eq.-1.0_R8.and.
     $              elm_ballooning_exp(kind_coeff).ne.0.0_R8) then
                     elm_poloidal_scaling_factor(ix,kind_coeff) = 0.0_R8
                   else
                     elm_poloidal_scaling_factor(ix,kind_coeff) = 
     &                 abs(cos((theta(ix,0)-2.0_R8*acos(-1.0_R8))
     &                 /2.0_R8))**elm_ballooning_exp(kind_coeff)
                   endif

                   if (elm_poloidal_scaling_factor(ix,kind_coeff)
     $                 .lt.1.0e-10_R8) then
                       elm_poloidal_scaling_factor(ix,kind_coeff)=0.0_R8
                   endif

                   do idy=1,ndim
                      pelm_rescaled(ix,idy)=
     &                  pf(idy)*(1.0_R8-
     &                  elm_poloidal_scaling_factor(ix,kind_coeff))+
     &                  pelm(idy)*
     &                  elm_poloidal_scaling_factor(ix,kind_coeff)
                   enddo

                   if(.not.
     1              (.not.region_flags(iir,kind_coeff).or.
     1              (lpflux.and.no_pflux).or.
     1              (ldiv.and.no_div))) then

                    if (elm_crash_region) then
                      if (elm_dynamics_multiple) then
                        frac_elm =
     &                     log(1.0_R8+elm_crash_factor(kind_coeff)*
     &                     ((mod(tim,elm_time_period)-
     &                     elm_time_begin(kind_coeff))/
     &                     (elm_time_plateau(kind_coeff)-
     &                     elm_time_begin(kind_coeff))))/
     &                     log(1.0_R8+elm_crash_factor(kind_coeff))
                      else
                        frac_elm =
     &                     log(1.0_R8+elm_crash_factor(1)*
     &                     ((mod(tim,elm_time_period)-
     &                     elm_time_begin(1))/
     &                     (elm_time_plateau(1)-
     &                     elm_time_begin(1))))/
     &                     log(1.0_R8+elm_crash_factor(1))
                      endif
                      do idy=1,ndim
                        prof(ix,idy)=
     &                    pf(idy)+(frac_elm)*
     &                    (pelm_rescaled(ix,idy)-pf(idy))
                      enddo
                    elseif (elm_plateau_region) then
                      prof(ix,1:ndim)=pelm_rescaled(ix,1:ndim)
                    elseif (elm_recovery_region) then
                      if (elm_dynamics_multiple) then
                        frac_elm =
     &                     (exp(-elm_recovery_factor(kind_coeff)*
     &                     ((mod(tim,elm_time_period)-
     &                     elm_time_recovery(kind_coeff))
     &                     /(elm_time_end(kind_coeff)-
     &                     elm_time_recovery(kind_coeff))))-
     &                     exp(-elm_recovery_factor(kind_coeff)))/
     &                     (1.0_R8-
     &                     exp(-elm_recovery_factor(kind_coeff)))
                      else
                        frac_elm =
     &                     (exp(-elm_recovery_factor(1)*
     &                     ((mod(tim,elm_time_period)-
     &                     elm_time_recovery(1))
     &                     /(elm_time_end(1)-
     &                     elm_time_recovery(1))))-
     &                     exp(-elm_recovery_factor(1)))/
     &                     (1.0_R8-
     &                     exp(-elm_recovery_factor(1)))
                      endif
                      do idy=1,ndim
                        prof(ix,idy)=
     &                    pf(idy)+(frac_elm)*
     &                    (pelm_rescaled(ix,idy)-pf(idy))
                      enddo
                    else
                      prof(ix,1:ndim)=pf(1:ndim)
                    endif

                    if (kind_coeff.eq.1) then
                      dna0(ix,iiy,nsp(is))=
     1                 prof(ix,iy)*poloidal_scale_factor

                    elseif (kind_coeff.eq.2) then
                      dpa0(ix,iiy,nsp(is))=
     1                 prof(ix,iy)*poloidal_scale_factor/
     &                 (rza(ix,iiy,nsp(is))*te(ix,iiy)+ti(ix,iiy))

                    elseif (kind_coeff.eq.3) then
                      hcib(ix,iiy,nsp(is))=poloidal_scale_factor*
     &                 prof(ix,iy)*na(ix,iiy,nsp(is))

                    elseif (kind_coeff.eq.4) then
                      hce0(ix,iiy)=prof(ix,iy)*ne(ix,iiy)*
     1                 poloidal_scale_factor

                    elseif (kind_coeff.eq.5) then
                      vla0(ix,iiy,0,nsp(is))=
     1                 prof(ix,iy)*poloidal_scale_factor

                    elseif (kind_coeff.eq.6) then
                      vla0(ix,iiy,1,nsp(is))=
     1                 prof(ix,iy)*poloidal_scale_factor

                    elseif (kind_coeff.eq.7) then
                      vsa0(ix,iiy,nsp(is))=poloidal_scale_factor*
     1                 prof(ix,iy)*mp*am(nsp(is))*na(ix,iiy,nsp(is))

                    elseif (kind_coeff.eq.8) then
                      if (csig_an_style.eq.0) then
                        sig0(ix,iiy)=prof(ix,iy)*qe*ne(ix,iiy)*
     1                   poloidal_scale_factor
                      elseif(csig_an_style.eq.1) then
                        sig0(ix,iiy)=prof(ix,iy)*qe*ne(nmdpl,-1)*
     1                   poloidal_scale_factor
                      endif

                    elseif (kind_coeff.eq.9) then
                      alf0(ix,iiy)=poloidal_scale_factor*
     1                 prof(ix,iy)*ne(ix,iiy)*
     2                 sqrt(qe/te(ix,iiy))

                    else
                      write(*,*)
     .                 'Option',kind_coeff,'not yet implemented'
                    endif

                   endif
                  endif
                enddo
              enddo
            enddo
            if (ncall.eq.0) then
              if (elm_data) then
                write(hlp_frm,'(a,i3,a)')
     .           '(a,1x,a,',nspec,'(1x,i3),3f14.7)'
                do iy=1,ndim
                  if (kind_coeff.eq.1)
     >              write (*,hlp_frm) 'transport=Dna','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.2)
     >              write (*,hlp_frm) 'transport=Dpa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.3)
     >              write (*,hlp_frm) 'transport=hci','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.4)
     >              write (*,'(a,3f14.7)') 'transport=hce',
     &              pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.5)
     >              write (*,hlp_frm) 'transport=vlax','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.6)
     >              write (*,hlp_frm) 'transport=vlay','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.7)
     >              write (*,hlp_frm) 'transport=vsa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.8)
     >              write (*,hlp_frm) 'transport=sig','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.9)
     >              write (*,hlp_frm) 'transport=alf','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                enddo
              else
                write(hlp_frm,'(a,i3,a)')
     .           '(a,1x,a,',nspec,'(1x,i3),2f14.7)'
                do iy=1,ndim
                  if (kind_coeff.eq.1)
     >              write (*,hlp_frm) 'transport=Dna','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.2)
     >              write (*,hlp_frm) 'transport=Dpa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.3)
     >              write (*,hlp_frm) 'transport=hci','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.4)
     >              write (*,'(a,2f14.7)') 'transport=hce',
     &              pr(iy),pf(iy)
                  if (kind_coeff.eq.5)
     >              write (*,hlp_frm) 'transport=vlax','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.6)
     >              write (*,hlp_frm) 'transport=vlay','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.7)
     >              write (*,hlp_frm) 'transport=vsa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.8)
     >              write (*,hlp_frm) 'transport=sig','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.9)
     >              write (*,hlp_frm) 'transport=alf','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                enddo
              endif
            endif

!           deallocate(f,r,pr,pf,felm,pelm)
            iloop=iloop+1
    1       continue
          enddo
        enddo
      enddo

      goto 100
   99 write (*,*) 'Error in reading the dsp file'
  100 continue

      if (new_files) then
        close(93)
        close(94)
        new_files = .false.
      end if

      return

*-----------------------------------------------------------------------

      contains

      subroutine read_transport_input
      implicit none
      integer jxa, jxi, jsep
      logical file_ok

c*** Clear data to forget previous profiles
      ndata(:,:,:)=0
      tdata(:,:,:,:)=0.0_R8
      addspec(:,:,:)=-5
      region_flags=.true.
      no_pflux=.false.
      no_div=.false.
      poloidal_scaling=.false.
      scaling_strength=0.0_R8
      scaling_width=0.0_R8
      scaling_ix_begin=-2
      scaling_ix_end=-2
      elm_dynamics_multiple=.false.
      elm_time_begin(:)=0.0_R8
      elm_time_plateau(:)=0.0_R8
      elm_time_recovery(:)=0.0_R8
      elm_time_end(:)=0.0_R8
      elm_time_period=0.0_R8
      elm_ix_begin=-2
      elm_ix_end=-2
      call get_jsep(nx,ny,jxi,jxa,jsep)
      elm_ixref = jxa
      elm_crash_factor(:)=1.0e-10_R8
      elm_recovery_factor(:)=1.0e-10_R8
      elm_ballooning_exp(:)=0.0_R8
      transport_time_mod=0.0_R8
      transport_time_switch=0.0_R8
      call find_file(filename,file_ok)
      if (file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist transport from ',trim(filename)
#ifdef USE_MPI
        read(99,transport,err=999)
#else
        read(99,transport)
#endif
        write(*,*) 'Read namelist transport from ',trim(filename)
        write(*,transport)
        close(99)
      else
        call xerrab (' No '//trim(filename))
      endif
      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading TRANSPORT inputfile!')
#endif
      end subroutine read_transport_input

*-----------------------------------------------------------------------

      subroutine transport_profile(kind_data,ndat,r,ndim,pr,nx,ny,ifa,
     & nfitanf)

      use b2mod_indirect
      use b2mod_geo
      use b2mod_mwti
     & , only: output_ds
      implicit none
      integer:: i,nx,ny,ifa
      integer:: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: pra(1:ny+2),prdsa(1:ny+2),r(nrr),pr(nrr),prtmp
      real (kind=R8) :: dsp(-1:ny)
      integer jxi, jxa, jsep, iyastrt, iyaend, ncall, nya
      integer ixref, iyref
      character*256 filename
      logical exist
      external ipgetr
      external get_jsep
      save jxi, jxa, jsep, ixref, iyref, iyastrt, iyaend, ncall
      data ncall/0/

      call subini('transport_profile')
      if (ncall.eq.0) then
        call get_jsep(nx,ny,jxi,jxa,jsep)
        call xertst (-1.le.jsep.and.jsep.le.ny,
     2   'faulty parameter jsep')
        call ipgeti ('b2mwti_jxa',jxa)
        call xertst (-1.le.jxa.and.jxa.le.nx,
     2   'faulty internal parameter jxa')
        ixref = jxa
        iyref = jsep
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
        call output_ds(ny,ixref,0,iyref,iyastrt,iyaend,'dsp',dsp)
        call ipgeti ('b2trno_csig_an_style', csig_an_style)
      endif

      filename='dsp'
      call find_file(filename,exist)
      if (exist) then
        open(unit=91,file=filename)
      else
        ifa=1
        call subend()
        return
      end if
      i = 0
      nya = 0
   91 continue
      i = i + 1
      read (91,*,END=92,ERR=97) prtmp
      prdsa(i) = prtmp
      nya = i
      if (i.lt.ny+2) goto 91
   92 continue
      nfitanf=1
      IFA=0
      do i=1,nya
! input in metres from separatrix
        if (kind_data.eq.1) pra(i)=prdsa(i)
! input in poloidal rho
        if (kind_data.eq.2) then
          write(*,*) 'Input profile in poloidal rho not yet coded!'
          goto 97
        endif

        if (pra(i).LT.r(1)) nfitanf=nfitanf+1
        nfitend=i-1
        if (pra(i).GT.r(ndat)) exit
        nfitend=i
      enddo

      if (nfitend.LE.nfitanf) goto 97

      ndim=nfitend+1-nfitanf
!     allocate (pr(1:ndim))

      do i=1,ndim
        pr(i)=pra(nfitanf-1+i)
      enddo
      goto 108
   97 ifa=1

  108 continue
      close (91)
      ncall=ncall+1
      call subend ()
      return
      end subroutine transport_profile

      end subroutine transport_input

      end module b2mod_input_profile

