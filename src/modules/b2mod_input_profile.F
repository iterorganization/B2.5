      module b2mod_input_profile
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_time
      use b2mod_sources 
     &     , only: ext_sna, ext_smo, ext_she, ext_shi, ext_sch, ext_sne
      implicit none
      integer nrr, nxx, nss
      integer, parameter :: nkind_source=4,nkind_data=2
      integer, parameter :: nkind_coeff=9
#ifdef ALLOCATE_AND_NAMELIST
      integer, allocatable, save :: nsdata(:,:,:)
      integer, allocatable, save :: nxdata(:,:,:)
      real(kind=R8), allocatable, save :: sdata(:,:,:,:)
      real(kind=R8), allocatable, save :: xdata(:,:,:,:)
      integer, allocatable, save  :: ndata(:,:,:)
      real(kind=R8), allocatable, save :: tdata(:,:,:,:)
      integer, allocatable, save :: addspec(:,:,:)
#else
#     include "DIMENSIONS.F"
      parameter(nrr = DEF_NYD+2, nxx = DEF_NXD, nss = DEF_NSD)
      integer, save :: nsdata(nkind_data,nkind_source,0:nss)
      integer, save :: nxdata(nkind_data,nkind_source,0:nss)
      real(kind=R8), save :: sdata(2,nrr,nkind_source,0:nss)
      real(kind=R8), save :: xdata(2,nxx,nkind_source,0:nss)
      integer, save :: ndata(nkind_data,nkind_coeff,0:nss)
      real(kind=R8), save :: tdata(3,nrr,nkind_coeff,0:nss)
      integer, save :: addspec(nss,nkind_coeff,0:nss)
#endif

      contains

      subroutine alloc_input_profile(ny,ns)
      implicit none
      integer ny, ns

#ifdef ALLOCATE_AND_NAMELIST
      nrr = ny+2
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
#endif

      return
      end subroutine alloc_input_profile

      subroutine source_input(nx,ny,ns,
     &  sna0,smo0,she0,shi0,sch0,sne0,vol,nc)
!**************************************************************************
!*     subroutine for reading in source-profiles of density, parallel     *
!*     momentum and/or energy into B2.5 (written by Heimo Buerbaumer 1999)*
!*     Inputfile:  b2.sources.profile                                     *
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
      logical catch_up
      logical, parameter :: debug = .false.
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
! read profile
        call read_source_profile
        open (96,file='../dx2')
        open (95,file='../dx1')
        open (94,file='../da1')
        open (93,file='../da2')
      endif    
   
      if(nc.eq.0) then
        if(sources_time_mod.gt.0.0_R8) then
          elm_count=tim/sources_time_mod
        else
          elm_count=0
        endif
        catch_up=sources_time_mod.gt.0.0_R8
        do while(catch_up)
          if(mod(tim,sources_time_mod).gt.sources_time_switch) then
            write(*,*) 'sources: catching up...'
            call read_source_profile
            catch_up=sources_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo
      else
        if(sources_time_mod.gt.0.0_R8) then
          write(*,'(a,1p,3e12.4,2i6)') 'SOURCES: ',
     .     tim,mod(tim,sources_time_mod),
     1     sources_time_switch,int(tim/sources_time_mod),elm_count
          if(mod(tim,sources_time_mod).gt.sources_time_switch.or.
     1     int(tim/sources_time_mod).gt.elm_count) then
            elm_count=tim/sources_time_mod 
            call read_source_profile
          endif
        endif
      endif
 
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
     .     'Increase NYD in DIMENSIONS.F to match size of profile')
!       allocate (f(1:ndat))
!       allocate (r(1:ndat))
          do i=1,ndat
            r(i)=sdata(1,i,kind_source,spec)
            f(i)=sdata(2,i,kind_source,spec)
            if (nc.eq.0) write(93,*) r(i),f(i)
          enddo
!          write(*,*) 'vor region_profile'
          call region_profile
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
     .       'Increase NXD in DIMENSIONS.F to match size of profile')
            do i=1,nxdat
              x(i)=xdata(1,i,kind_source,spec)
              fx(i)=xdata(2,i,kind_source,spec)
              if (nc.eq.0) write(95,*) x(i),fx(i)
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
              if (nc.eq.0 .and. iy.eq.1. .or. debug) then
                if (kind_source.eq.1) write(*,*) 'AXIAL-SOURCES: sna ',
     &           spec,xr(ix),xf(ix)
                if (kind_source.eq.2) write(*,*) 'AXIAL-SOURCES: smo ',
     &           spec,xr(ix),xf(ix)
                if (kind_source.eq.3) write(*,*) 'AXIAL-SOURCES: she ',
     &           xr(ix),xf(ix)
                if (kind_source.eq.4) write(*,*) 'AXIAL-SOURCES: shi ',
     &           xr(ix),xf(ix)
                if (kind_source.eq.5) write(*,*) 'AXIAL-SOURCES: sch ',
     &           xr(ix),xf(ix)
                if (kind_source.eq.6) write(*,*) 'AXIAL-SOURCES: sne ',
     &           xr(ix),xf(ix)
              endif
              if (nc.eq.0) write (96,*) xr(ix),xf(ix)
            enddo
            if (nc.eq.0 .or. debug) then
              if (kind_source.eq.1) write(*,*) 'PROFILE-SOURCES: sna ',
     &         spec,pr(iy),pf(iy)
              if (kind_source.eq.2) write(*,*) 'PROFILE-SOURCES: smo ',
     &         spec,pr(iy),pf(iy)
              if (kind_source.eq.3) write(*,*) 'PROFILE-SOURCES: she ',
     &         pr(iy),pf(iy)
              if (kind_source.eq.4) write(*,*) 'PROFILE-SOURCES: shi ',
     &         pr(iy),pf(iy)
              if (kind_source.eq.5) write(*,*) 'PROFILE-SOURCES: sch ',
     &         pr(iy),pf(iy)
              if (kind_source.eq.6) write(*,*) 'PROFILE-SOURCES: sne ',
     &         pr(iy),pf(iy)
            endif
            if (nc.eq.0) write (94,*) pr(iy),pf(iy)
          enddo
 1        continue
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
 98   write (*,*) 'Error in computing axial profile'
      goto 100
 99   write (*,*) 'Error in reading the dsa file'
 100  continue

      ext_sna(:,:,:)=sna0(:,:,0,:)
      ext_smo(:,:,:)=smo0(:,:,0,:)
      ext_sne(:,:)=sne0(:,:,0)
      ext_she(:,:)=she0(:,:,0)
      ext_shi(:,:)=shi0(:,:,0)
      ext_sch(:,:)=sch0(:,:,0)
       
      close(93)
      close(94)
      close(95)
      close(96)
      return

      contains

      subroutine read_source_profile

      character*260 filename
      logical file_ok

      namelist/profile/nsdata,nxdata,sdata,xdata,divheat,
     . sources_time_mod,sources_time_switch,sources_filename

!  read the profile data
      filename='../' // trim(sources_filename)
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(92,file=filename)
        write(*,*) 
     .   'Reading sources profile namelist from ', trim(filename)
        read(92,profile)
        write(*,*) 
     .   'Read sources profile namelist from ',trim(filename)
        close(92)
      else
        call xerrab ('No '//trim(filename))
      endif
!     write(*,profile)
      return

      end subroutine read_source_profile

      end subroutine source_input

*-----------------------------------------------------------------------

      subroutine transport_input(nx, ny, ns, rza, na, ne, te, ti, hcib,
     &  dna0, dpa0, vla0, vsa0, hci0, hce0, sig0, alf0, ncall)

!**************************************************************************
!*     subroutine for reading in experimental transport coefficients from *
!*     the midplane and calculating transport coefficients                *
!*     for B2.5 (written by Heimo Buerbaumer 1999)                        *
!*     Inputfile:  b2.transport.inputfile                                 *
!**************************************************************************
      use b2mod_time
      use b2mod_constants
      use b2mod_indirect
      use b2mod_b2cmpa
      IMPLICIT NONE
      integer :: nx,ny,ns
      integer :: iy,ix,i,ifail,spec,nspec,iiy,iir,
     &      nfitanf,ndim,is,nsp(1:ns+1),iloop,ncall,ifa,ndat,
     &      kind_coeff,kind_data
      integer ixwidth, ixmean

      real(kind=R8), dimension(-1:nx,-1:ny,0:ns-1) :: 
     &    rza,dna0,dpa0,vsa0,na
      real(kind=R8), dimension(-1:nx,-1:ny) :: 
     &    hci0,hce0,sig0,alf0,ne,te,ti
       
      real(kind=R8) :: scaling_width, scaling_strength
      real(kind=R8) :: f(nrr),r(nrr),pr(nrr),pf(nrr),d(nrr),
     &    pelm(nrr),felm(nrr),delm(nrr),prof(nrr),
     &    vla0(-1:nx,-1:ny,0:1,0:ns-1),hcib(-1:nx,-1:ny,0:ns-1)
      character*256, save :: transport_filename
      real (kind=R8), save :: transport_time_mod
      real (kind=R8), save :: transport_time_switch
      real (kind=R8), save :: elm_time_begin, elm_time_end
      real (kind=R8), save :: elm_time_period
      real (kind=R8) :: poloidal_scale_factor
      integer, save :: scaling_ix_begin, scaling_ix_end
      integer, save :: elm_ix_begin, elm_ix_end
      logical, save :: region_flags(NNREGMAX,nkind_coeff),
     1 no_pflux, poloidal_scaling, elm_time, elm_region, elm_data
      character*260 filename
      character*26 hlp_frm
      integer, save :: elm_count
      logical catch_up, lpflux
      intrinsic sqrt
      namelist /transport/ ndata, tdata, addspec,
     1   transport_filename, transport_time_mod, transport_time_switch,
     2   scaling_ix_begin, scaling_ix_end, scaling_strength,
     3   region_flags, no_pflux, poloidal_scaling, scaling_width,
     4   elm_time_begin, elm_time_end, elm_time_period,
     5   elm_ix_begin, elm_ix_end

      data transport_filename/'b2.transport.inputfile'/,
     1 transport_time_mod/0.0_R8/,transport_time_switch/0.0_R8/

      if(ncall.eq.0) then
        call xertst(NNREGMAX.ge.nnreg(0),
     1   'Increase NNREGMAX in b2mod_numerics_namelist')
        filename='../' // transport_filename
        open (94,file='../da1')
        open (93,file='../da2')
!       open (91,file='../dsa')
        call read_transport_input         
!  read the profile data
        if(transport_time_mod.gt.0.0_R8) then
          elm_count=tim/transport_time_mod
        else
          elm_count=0
        endif
        catch_up=transport_time_mod.gt.0.0_R8
        do while(catch_up) 
          if(mod(tim,transport_time_mod).gt.transport_time_switch)
     1     then           
            write(*,*) 'transport_input: catching up ...'
            filename='../' // transport_filename
            call read_transport_input
            catch_up=transport_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo
      else
        if(transport_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,transport_time_mod),
     1     transport_time_switch,int(tim/transport_time_mod),elm_count
          if(mod(tim,transport_time_mod).gt.transport_time_switch.or.
     1     int(tim/transport_time_mod).gt.elm_count) then
            elm_count=tim/transport_time_mod
            filename='../' // transport_filename
            call read_transport_input
          endif
        endif
      endif
      ixwidth = scaling_ix_end - scaling_ix_begin + 1
      ixmean = (scaling_ix_end + scaling_ix_begin) / 2
      if (scaling_width.eq.0.0_R8 .and. ixwidth.gt.1)
     > scaling_width = 
     .  0.5_R8*sqrt((ixwidth+1.0_R8)*(ixwidth-1.0_R8)/3.0_R8) 
      if (elm_time_period.gt.0.)then
        elm_time=mod(tim,elm_time_period).ge.elm_time_begin.and.
     &           mod(tim,elm_time_period).le.elm_time_end
      else  
        elm_time=.false.
      endif 
      elm_data = maxval(tdata(3,:,:,:)).gt.0.0_R8  

!        write (*,*) 'vor loop'
      iloop=0
      do kind_coeff=1,nkind_coeff
        do spec=0,ns-1
          do kind_data=1,nkind_data
!           write (*,*) 'loop-begin', kind_coeff,spec,kind_data
            ndat=ndata(kind_data,kind_coeff,spec)
            if (ndat.eq.0) goto 1
            call xertst(ndat.le.nrr,
     .       'Increase NYD in DIMENSIONS.F to match size of profile')
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
!             if (iloop.eq.2) write(93,*) r(i),f(i),felm(i)
            enddo
              
            call region_profile(1,ndat,r,ndim,pr,nx,ny,ifa,nfitanf)
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

            if (iloop.eq.1) then

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
                  if(iir.le.0.or.iir.gt.NNREGMAX) then
                    write(*,*) 'transport_input: out of bounds ',
     1               'ix,iy,region(ix,iy) ',ix,iiy,iir
                    iir=NNREGMAX
                  endif
                  lpflux=mod(iir-1,4)+1.ge.3.and.
     2             (((ix.lt.leftcut(1).or.ix.ge.rightcut(1)).and.
     3               iiy.lt.topcut(1)).or.(nnreg(0).eq.8.and.
     4               ix.ge.leftcut(2).and.ix.lt.rightcut(2).and.
     5               iiy.lt.topcut(2)))
                  if (poloidal_scaling .and.
     &                ix.ge.scaling_ix_begin .and. ix.le.scaling_ix_end) 
     &              then
                    poloidal_scale_factor = 1.0_R8 + scaling_strength/
     &                    (sqrt(2.0_R8*pi)*scaling_width)*
     &                    exp(-(ix-ixmean)*(ix-ixmean)/
     &                         (2.0_R8*(scaling_width**2)))
                  else
                    poloidal_scale_factor = 1.0_R8
                  endif
                  elm_region=elm_time.and.elm_data.and.
     &                     ix.ge.elm_ix_begin.and.ix.le.elm_ix_end

                  if(.not.
     1             (.not.region_flags(iir,kind_coeff).or.
     1             (lpflux.and.no_pflux))) then

                    if (elm_region) then
                      prof(1:ndim)=pelm(1:ndim)
                    else
                      prof(1:ndim)=pf(1:ndim)
                    endif

                    if (kind_coeff.eq.1) then
                      dna0(ix,iiy,nsp(is))=
     1                 prof(iy)*poloidal_scale_factor

                    elseif (kind_coeff.eq.2) then
                      dpa0(ix,iiy,nsp(is))=
     1                 prof(iy)*poloidal_scale_factor/
     &                 (rza(ix,iiy,nsp(is))*te(ix,iiy)+ti(ix,iiy))

                    elseif (kind_coeff.eq.3) then         
                      hcib(ix,iiy,nsp(is))=poloidal_scale_factor*
     &                 prof(iy)*na(ix,iiy,nsp(is))
             
                    elseif (kind_coeff.eq.4) then
                      hce0(ix,iiy)=prof(iy)*ne(ix,iiy)*
     1                 poloidal_scale_factor
                    
                    elseif (kind_coeff.eq.5) then
                      vla0(ix,iiy,0,nsp(is))=
     1                 prof(iy)*poloidal_scale_factor
                
                    elseif (kind_coeff.eq.6) then
                     vla0(ix,iiy,1,nsp(is))=
     1                prof(iy)*poloidal_scale_factor

                    elseif (kind_coeff.eq.7) then
                      vsa0(ix,iiy,nsp(is))=poloidal_scale_factor*
     1                 prof(iy)*mp*am(nsp(is))*na(ix,iiy,nsp(is))
                 
                    elseif (kind_coeff.eq.8) then
                      sig0(ix,iiy)=prof(iy)*qe*ne(ix,iiy)*
     1                 poloidal_scale_factor
                
                    elseif (kind_coeff.eq.9) then
                      alf0(ix,iiy)=poloidal_scale_factor*
     1                 prof(iy)*ne(ix,iiy)*
     2                 sqrt(qe/te(ix,iiy))
              
                    else
                      write(*,*) 
     .                 'Option',kind_coeff,'not yet implemented'
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
  1         continue
          enddo
        enddo
      enddo
       
       
      goto 100
 99   write (*,*) 'Error in reading the dsa file'
 100  continue
       
        
      close(99)
      close(93)
      close(94)
      return

*-----------------------------------------------------------------------

      contains

      subroutine read_transport_input
      implicit none
      logical file_ok

c*** Clear data to forget previous profiles
      ndata(:,:,:)=0
      tdata(:,:,:,:)=0.0_R8
      addspec(:,:,:)=-5
      region_flags=.true.
      no_pflux=.false.
      poloidal_scaling=.false.
      scaling_strength=0.0_R8
      scaling_width=0.0_R8
      scaling_ix_begin=-2
      scaling_ix_end=-2
      elm_time_begin=0.0_R8
      elm_time_end=0.0_R8
      elm_time_period=0.0_R8
      elm_ix_begin=-2
      elm_ix_end=-2
      inquire(file=filename,exist=file_ok)
      if (file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist transport from ',trim(filename)
        read(99,transport)
        write(*,*) 'Read namelist transport from ',trim(filename)
        write(*,transport)
      else
        call xerrab (' No '//trim(filename))
      endif
      return
      end subroutine read_transport_input

      end subroutine transport_input

*-----------------------------------------------------------------------

      subroutine region_profile(kind_data,ndat,r,ndim,pr,nx,ny,ifa,
     & nfitanf)

      use b2mod_indirect
      use b2mod_geo
      implicit none
      integer:: i,nx,ny,ifa
      integer:: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: pra(1:ny+2),prdsa(1:ny+2),r(nrr),pr(nrr),prtmp
      integer jxi, jxa, jsep, iyastrt, iyaend, ncall, nya
      integer ixref, iyref
      logical ex
      external ipgeti, ipgetr, xertst, subini, subend
      external output_ds, get_jsep
      save jxi, jxa, jsep, ixref, iyref, iyastrt, iyaend, ncall
      data ncall/0/

      call subini('region_profile')
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
        call output_ds(crx,cry,nx,ny,ixref,0,
     1   iyref,iyastrt,iyaend,'../dsp')
      endif

      open (91,file='../dsp')
      i = 0
      nya = 0
 91   continue
      i = i + 1
      read (91,*,END=92,ERR=97) prtmp
      prdsa(i) = prtmp
      nya = i
      if (i.lt.ny+2) goto 91
 92   continue
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
      close (91)
      goto 108
 97   ifa=1
         
 108  continue
      ncall=ncall+1
      call subend ()
      return
      end subroutine region_profile

*-----------------------------------------------------------------------

      subroutine axial_profile(kind_data,ndat,x,ndim,pr,nx,ny,ifa,
     & nfitanf,nfitend)

      use b2mod_indirect
      use b2mod_geo
      implicit none
      integer :: i,nx,ny,ifa,nxm
      integer :: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: x(1:ndat),pr(-1:nxx-2)
      real(kind=R8) :: l_conn(-1:nx), l_total, l_ref
      integer jxi, jxa, jsep, ncall
      integer ixref, iyref
      integer ireg
      external ipgeti, ipgetr, xertst, subini, subend
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
 97   ifa=1
         
 108  continue
      ncall=ncall+1
      call subend ()
      return
      end subroutine axial_profile

      end module b2mod_input_profile
      
