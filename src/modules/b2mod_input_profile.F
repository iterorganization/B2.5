      module b2mod_input_profile
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_sources 
     &     , only: ext_sna, ext_smo, ext_she, ext_shi, ext_sch
      implicit none
#     include "DIMENSIONS.F"

      contains

      subroutine source_input(nx,ny,ns,sna0,smo0,she0,shi0,sch0,vol,nc)
!**************************************************************************
!*     subroutine for reading in source-profiles of density, momentum     *
!*     and/or energy into B2.5 (written by Heimo Buerbaumer 1999)         *
!*     Inputfile:  b2.sources.profile                                     *
!**************************************************************************
      integer :: nx,ny,ns,nc
      integer,parameter :: nkind_source=4,nkind_data=2
      integer,parameter :: nrm=DEF_NYD+2,nsm=DEF_NSD
      integer :: ix,iy,i,ireg
      integer :: nanf,ndim,ifail,is,ndat,kind_data,iloop
      integer :: kind_source,spec,ifa
      integer :: nsdata(nkind_data,nkind_source,0:nsm)
      real(kind=R8), dimension(-1:nx,-1:ny,0:3) :: she0,shi0,sch0
      real(kind=R8) :: sna0(-1:nx,-1:ny,0:1,0:ns-1)
      real(kind=R8) :: smo0(-1:nx,-1:ny,0:3,0:ns-1),d(nrm)
      real(kind=R8) :: sdata(2,nrm,nkind_source,0:nsm)
      real(kind=R8) :: vol(-1:nx,-1:ny),totpower,totna(0:ns-1)
      real(kind=R8) :: f(nrm),r(nrm),pr(nrm),pf(nrm),divheat,dhtot
      intrinsic mod
     
      namelist/profile/nsdata,sdata,divheat

! initialise sources and data to zero
      nsdata(:,:,:)=0
      sdata(:,:,:,:)=0.0_R8
      sna0(:,:,:,:)=0.0_R8
      smo0(:,:,:,:)=0.0_R8
      she0(:,:,:)=0.0_R8
      shi0(:,:,:)=0.0_R8
      sch0(:,:,:)=0.0_R8
      divheat=0.0_R8
      
      open (94,file='../da1')
      open (93,file='../da2')
      open (92,file='../b2.sources.profile')
   
      iloop=0
!  read the profile data
      read(92,profile)
!        write(*,profile)
 
      do kind_source=1,nkind_source
        do spec=0,ns-1
          do kind_data=1,nkind_data   
            ndat=nsdata(kind_data,kind_source,spec)
            if (ndat.eq.0) goto 1
            call xertst(ndat.le.nrm,
     .       'Increase NYD in DIMENSIONS.F to match size of profile')
!       allocate (f(1:ndat))
!       allocate (r(1:ndat))
            do i=1,ndat
              r(i)=sdata(1,i,kind_source,spec)
              f(i)=sdata(2,i,kind_source,spec)
              if (iloop.eq.2) write(93,*) r(i),f(i)
            enddo
!        write(*,*) 'vor region_profile'
            call region_profile
     .       (kind_data,ndat,r,ndim,pr,nx,ny,ifa,nrm,nanf)
!        write (*,*) 'ifa=',ifa
            if (ifa.eq.1) goto 99
!        allocate (pf(1:ndim))
            call E01BEF(ndat,r,f,d,ifail)
!       write(*,*) 'zwischen,ifail=',ifail
            call E01BFF(ndat,r,f,d,ndim,pr,pf,ifail)

            do iy=1,ndim
              do ix=-1,nx
                ireg=region(ix,nanf-3+iy,0)
                if (bottomiy(ix,nanf-3+iy).ne.-2 .and.
     &              topiy(ix,nanf-3+iy).ne.ny+1 .and.
     &              leftix(ix,nanf-3+iy).ne.-2 .and.
     &              rightix(ix,nanf-3+iy).ne.nx+1 .and.
     &              ireg.ne.0) then
                  if (kind_source.eq.1) then
                    if (pr(iy).le.0.0_R8 .and. mod(ireg,4).eq.1) then
                      sna0(ix,nanf-3+iy,0,spec)=pf(iy)*vol(ix,nanf-3+iy)
                    else if 
     .                 (pr(iy).ge.0.0_R8 .and. mod(ireg,4).eq.2) then
                      sna0(ix,nanf-3+iy,0,spec)=pf(iy)*vol(ix,nanf-3+iy)
                    endif
                  elseif (kind_source.eq.2) then
                    if (pr(iy).le.0.0_R8 .and. mod(ireg,4).eq.1) then
                      smo0(ix,nanf-3+iy,0,spec)=pf(iy)*vol(ix,nanf-3+iy)
                    else if 
     .                 (pr(iy).ge.0.0_R8 .and. mod(ireg,4).eq.2) then
                      smo0(ix,nanf-3+iy,0,spec)=pf(iy)*vol(ix,nanf-3+iy)  
                    endif      
                  elseif (kind_source.eq.3) then 
                    if (pr(iy).le.0.0_R8 .and. mod(ireg,4).eq.1) then
                      she0(ix,nanf-3+iy,0)=pf(iy)*vol(ix,nanf-3+iy)
                    else if 
     .                 (pr(iy).ge.0.0_R8 .and. mod(ireg,4).eq.2) then
                      she0(ix,nanf-3+iy,0)=pf(iy)*vol(ix,nanf-3+iy)
                    endif
                  elseif (kind_source.eq.4 .and. mod(ireg,4).eq.1) then 
                    if (pr(iy).le.0.0_R8 .and. mod(ireg,4).eq.1) then
                      shi0(ix,nanf-3+iy,0)=pf(iy)*vol(ix,nanf-3+iy)
                    else if 
     .                 (pr(iy).ge.0.0_R8 .and. mod(ireg,4).eq.2) then
                      shi0(ix,nanf-3+iy,0)=pf(iy)*vol(ix,nanf-3+iy)
                    endif   
                  endif
                endif
              enddo
!         if (nc.eq.0) then
!         if (kind_source.eq.1) write(*,*) 'PROFILE-SOURCES: sna ',
!     &           spec,pr(iy),pf(iy)
!         if (kind_source.eq.2) write(*,*) 'PROFILE-SOURCES: smo ',
!     &           spec,pr(iy),pf(iy)
!         if (kind_source.eq.3) write(*,*) 'PROFILE-SOURCES: she ',
!     &           pr(iy),pf(iy)
!         if (kind_source.eq.4) write(*,*) 'PROFILE-SOURCES: shi ',
!     &           pr(iy),pf(iy)
!         
!         endif
!        if (iloop.eq.2) write (94,*) pr(iy),pf(iy)
            enddo
 1          iloop=iloop+1
          enddo
        enddo
      enddo
! Calculation of the total power input
      if (nc.eq.0) then
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
 99   write (*,*) 'Error in reading the dsa-file'
 100  continue

      ext_sna(:,:,:)=sna0(:,:,0,:)
      ext_smo(:,:,:)=smo0(:,:,0,:)
      ext_she(:,:)=she0(:,:,0)
      ext_shi(:,:)=shi0(:,:,0)
      ext_sch(:,:)=sch0(:,:,0)
       
      close (92)
      close(93)
      close(94)
      return
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
      integer,parameter :: nkind_coeff=9,nkind_data=2
      integer,parameter :: nrm=DEF_NYD+2,nsm=DEF_NSD
      integer :: nx,ny,ns,ndata(nkind_data,nkind_coeff,0:nsm)
      integer :: iy,ix,i,ifail,spec,nspec,iiy,iir,
     &      nfitanf,ndim,is,nsp(1:ns+1),iloop,ncall,ifa,ndat,
     &      kind_coeff,kind_data,addspec(nsm,nkind_coeff,0:nsm)
      integer ixwidth, ixmean

      real(kind=R8), dimension(-1:nx,-1:ny,0:ns-1) :: 
     &    rza,dna0,dpa0,vsa0,na
      real(kind=R8), dimension(-1:nx,-1:ny) :: 
     &    hci0,hce0,sig0,alf0,ne,te,ti
       
      real(kind=R8) :: tdata(3,nrm,nkind_coeff,0:nsm)
      real(kind=R8) :: scaling_width,scaling_strength
      real(kind=R8) :: f(nrm),r(nrm),pr(nrm),pf(nrm),d(nrm),
     &    pelm(nrm),felm(nrm),delm(nrm),prof(nrm),
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
      integer elm_count
      logical catch_up, lpflux
      intrinsic sqrt
      save ndata, tdata, addspec
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
            call xertst(ndat.le.nrm,
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
              
            call region_profile(1,ndat,r,ndim,pr,nx,ny,ifa,nrm,nfitanf)
!           write (*,*) 'ifa=',ifa
            if (ifa.eq.1) goto 99
        
            call E01BEF(ndat,r,f,d,ifail)
            if (elm_data) call E01BEF(ndat,r,felm,delm,ifail)
!           write(*,*) 'zwischen,ifail=',ifail
            call E01BFF(ndat,r,f,d,ndim,pr,pf,ifail)
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
                do iy=1,ndim
                  if (kind_coeff.eq.1) 
     >              write (*,*) 'transport=Dna','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.2) 
     >              write (*,*) 'transport=Dpa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.3) 
     >              write (*,*) 'transport=hci','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy) 
                  if (kind_coeff.eq.4) 
     >              write (*,*) 'transport=hce',
     &              pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.5) 
     >              write (*,*) 'transport=vlax','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.6) 
     >              write (*,*) 'transport=vlay','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.7) 
     >              write (*,*) 'transport=vsa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.8) 
     >              write (*,*) 'transport=sig','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                  if (kind_coeff.eq.9) 
     >              write (*,*) 'transport=alf','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy),pelm(iy)
                enddo
              else
                do iy=1,ndim
                  if (kind_coeff.eq.1) 
     >              write (*,*) 'transport=Dna','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.2) 
     >              write (*,*) 'transport=Dpa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.3) 
     >              write (*,*) 'transport=hci','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy) 
                  if (kind_coeff.eq.4) 
     >              write (*,*) 'transport=hce',
     &              pr(iy),pf(iy)
                  if (kind_coeff.eq.5) 
     >              write (*,*) 'transport=vlax','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.6) 
     >              write (*,*) 'transport=vlay','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.7) 
     >              write (*,*) 'transport=vsa','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.8) 
     >              write (*,*) 'transport=sig','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                  if (kind_coeff.eq.9) 
     >              write (*,*) 'transport=alf','spec=',
     &              (nsp(i),i=1,nspec),pr(iy),pf(iy)
                enddo
              endif
            endif

!           deallocate(f,r,pr,pf,felm,pelm)
            iloop=iloop+1
  1       enddo
        enddo
      enddo
       
      goto 100
 99   write (*,*) 'Error in reading the dsa-file'
 100  continue
       
      close(99)
      close(93)
      close(94)
      ncall=ncall+1
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
        write(*,*) ' No ',trim(filename)
      endif
      return
      end subroutine read_transport_input

      end subroutine transport_input

*-----------------------------------------------------------------------

      subroutine region_profile(kind_data,ndat,r,ndim,pr,nx,ny,ifa,nrm,
     & nfitanf)

      use b2mod_indirect
      use b2mod_geo
      implicit none
      integer:: i,nx,ny,ifa,nrm
      integer:: nfitanf,nfitend,ndim,ndat,kind_data
      real(kind=R8) :: pra(1:ny+2),prdsa(1:ny+2),r(nrm),pr(nrm),prtmp
      integer jxi, jxa, jsep, iyastrt, iyaend, ncall, nya
      logical ex
      external ipgeti, ipgetr, xertst, subini, subend
      external output_ds, get_jsep
      save jxi, jxa, jsep, iyastrt, iyaend, ncall
      data ncall/0/

      call subini('region_profile')
      if (ncall.eq.0) then
        inquire(file='../dsa',exist=ex)
        if (.not.ex) then
          call get_jsep(nx,ny,jxi,jxa,jsep)
          call xertst (-1.le.jsep.and.jsep.le.ny,
     2     'faulty parameter jsep')
          call ipgeti ('b2mwti_jxa',jxa)
          call xertst (-1.le.jxa.and.jxa.le.nx,
     2     'faulty internal parameter jxa')
          call output_ds(crx,cry,nx,ny,jxa,0,
     1     jsep,iyastrt,iyaend,'../dsa')
        endif
      endif

      open (91,file='../dsa')
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
! input in cm from separatrix
        if (kind_data.eq.1) pra(i)=prdsa(i)
! input in poloidal rho
!          if (kind_data.eq.2) pra(i)=

        if (pra(i).LT.r(1)) nfitanf=nfitanf+1
        nfitend=i-1
        if (pra(i).GT.r(ndat)) exit
        nfitend=i
      enddo


!        write (*,*) nfitanf,nfitend
      if (nfitend.LE.nfitanf) goto 97
           
      ndim=nfitend+1-nfitanf
!     write (*,*) 'ndim=',ndim
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

      end module b2mod_input_profile
      
