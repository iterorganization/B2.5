      module b2mod_sputter
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
      use b2mod_layer

      integer, save :: sput_dst, sput_dst2, sput_dst3
      integer, save :: target_intcode(3+6), target_index(92)
      integer, save :: zwall(6)
      character*2, save :: target_code(3+6)
      integer, allocatable, save ::
     *  sputter_index(:),
     *  sputter_length(:,:),
     *  reflect_length(:,:)
      integer, allocatable, save ::
     *  sputter_length2(:,:), reflect_length2(:,:)
      real (kind=R8), allocatable, save ::
     *  sputter_coefficient(:,:,:), reflect_coefficient(:,:,:),
     *  sputter_energy(:,:,:), reflect_energy(:,:,:)
      real (kind=R8), allocatable, save ::
     *  sputter_coefficient2(:,:,:), reflect_coefficient2(:,:,:),
     *  sputter_energy2(:,:,:), reflect_energy2(:,:,:)
      real (kind=R8), allocatable, save :: tgdns(:), es(:), theta(:),
     *  monolayer_l(:)
      logical, save :: sputter_namelist_used
      logical, save :: sputter_initialised
      data sputter_namelist_used /.false./
      data sputter_initialised /.false./
      data sput_dst/-1/, sput_dst2/-1/, sput_dst3/-1/

#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *  sputter_yield(:,:,:,:), sputter_yield2(:,:,:,:)
#else
#     include "DIMENSIONS.F"
      integer nxmax, nymax, nsdmax
      parameter (NXMAX=DEF_NXD,NYMAX=DEF_NYD,NSDMAX=DEF_NSD)
      real (kind=R8), save ::
     *  sputter_yield(-1:NXMAX,-1:NYMAX,0:NSDMAX-1,1:2), 
     *  sputter_yield2(-1:NXMAX,-1:NYMAX,0:NSDMAX-1,1:2)      
#endif
      integer, save :: jagged_rates=0

      integer, save :: nx_sp, ny_sp, ns_sp
      namelist /sputter_save/ nx_sp, ny_sp, ns_sp,
     *   sputter_yield, sputter_yield2

      contains

      subroutine alloc_b2mod_sputter(nxd,nyd,nsd)
      use b2mod_b2cmpa
      use b2mod_elements
      implicit none
      integer nxd, nyd, nsd
      logical errdim

      if(sputter_initialised) return
      call ipgeti('b2mod_sputter_jagged_rates', jagged_rates)
#ifdef ALLOCATE_AND_NAMELIST
      allocate(sputter_yield(-1:nxd,-1:nyd,0:nsd-1,2))
      allocate(sputter_yield2(-1:nxd,-1:nyd,0:nsd-1,2))
#else
      errdim=.false.
      if(nxd.gt.nxmax) then
        write(*,*) 'alloc_b2mod_sputter: nxd > nxmax',nxd,nxmax
        errdim=.true.
      endif
      if(nyd.gt.nymax) then
        write(*,*) 'alloc_b2mod_sputter: nyd > nymax',nyd,nymax
        errdim=.true.
      endif
      if(nsd.gt.nsdmax) then
        write(*,*) 'alloc_b2mod_sputter: nsd > nsdmax',nsd,nsdmax
        errdim=.true.
      endif
      call xertst (.not.errdim,
     1 'faulty NXMAX, NYMAX and/or NSDMAX')
#endif
      ntargsp=0
      target_index=0
      target_intcode=-1
      if(sput_dst.ge.0.and.sput_dst.lt.nsd) then
        ntargsp=ntargsp+1
        target_intcode(ntargsp)=sput_dst
        target_code(ntargsp)=is_codes(sput_dst)
        target_index(nint(zn(sput_dst)))=ntargsp
      endif
      if(sput_dst2.ge.0.and.sput_dst2.lt.nsd) then
        ntargsp=ntargsp+1
        target_intcode(ntargsp)=sput_dst2
        target_code(ntargsp)=is_codes(sput_dst2)
        target_index(nint(zn(sput_dst2)))=ntargsp
      endif
      if(sput_dst3.ge.0.and.sput_dst3.lt.nsd) then
        ntargsp=ntargsp+1
        target_intcode(ntargsp)=sput_dst3
        target_code(ntargsp)=is_codes(sput_dst3)
        target_index(nint(zn(sput_dst3)))=ntargsp
      endif
      call xertst(ntargsp.le.ntrack, 
     &     'alloc_b2mod_sputter: ntargsp > ntrack')
      allocate(monolayer_l(ntargsp))
      return
      end subroutine alloc_b2mod_sputter

      subroutine dealloc_b2mod_sputter
      implicit none
      if(.not.allocated(monolayer_l)) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(sputter_yield)
      deallocate(sputter_yield2)
#endif
      deallocate(monolayer_l)
      return
      end subroutine dealloc_b2mod_sputter

      subroutine read_b2mod_sputter_namelist(nread,edition,found)
      implicit none
      integer nread
      character*4 edition
      logical found
      character*256 filename
      logical file_ok

      if (edition.eq.'    ') then
        filename = '../b2.sputter_save.parameters'
      else
        filename = 'b2.sputter_save.parameters.'//edition
      endif
      inquire(file=trim(filename),exist=file_ok)
      if(file_ok) then
        open(nread,file=trim(filename))
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,sputter_save,err=999)
#else
        read(nread,sputter_save)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
      else
        write(*,*) 'No ', trim(filename)
      endif
      found = maxval(sputter_yield).gt.0.0_R8 .or.
     .        maxval(sputter_yield2).gt.0.0_R8
      return
 999  call xerrab ('Error while reading SPUTTER_SAVE namelist!')
      end subroutine read_b2mod_sputter_namelist

      subroutine write_b2mod_sputter_namelist(nwrite,nx,ny,ns,edition)
      implicit none
      integer nwrite, nx, ny, ns
      character*4 edition
      integer ix, iy, is, ib
      character*256 filename

      if(maxval(sputter_yield).gt.0.0_R8.or.
     .   maxval(sputter_yield2).gt.0.0_R8) then
        if (edition.eq.'    ') then
          filename = '../b2.sputter_save.parameters'
        else
          filename = 'b2.sputter_save.parameters.'//edition
        endif
        open(nwrite,file=trim(filename))
cxpb        write(nwrite,sputter_save)
cxpb  Replaced with explicit writing of nonzero elements only
        write(nwrite,'(a13)') '&SPUTTER_SAVE'
        write(nwrite,'(a8,i4,a1)') 'NX_SP = ', nx, ','
        write(nwrite,'(a8,i4,a1)') 'NY_SP = ', ny, ','
        write(nwrite,'(a8,i4,a1)') 'NS_SP = ', ns, ','
        do ix = -1, nx
          do iy = -1, ny
            do is = 0, ns-1
              if ( sputter_yield(ix,iy,is,1).ne.0.0_R8 .or.
     &             sputter_yield(ix,iy,is,2).ne.0.0_R8 .or.
     &            sputter_yield2(ix,iy,is,1).ne.0.0_R8 .or.
     &            sputter_yield2(ix,iy,is,2).ne.0.0_R8) then
                do ib = 1, 2
                  write(nwrite,'(a14,i4,a1,i4,a1,i3,a1,i1,a2,
     &                           1pe21.13e3,
     &                           a17,i4,a1,i4,a1,i3,a1,i1,a2,
     &                           1pe21.13e3,a1)')
     .             'SPUTTER_YIELD(',ix,',',iy,',',is,',',ib,')=',
     .              sputter_yield(ix,iy,is,ib),
     .           ', SPUTTER_YIELD2(',ix,',',iy,',',is,',',ib,')=',
     .              sputter_yield2(ix,iy,is,ib),','
                end do
              end if
            end do
          end do
        end do
        write(nwrite,'(a1)') '/'
        close(nwrite)
        write(*,*) 'Saved ',trim(filename)
      endif
      return
      end subroutine write_b2mod_sputter_namelist

      subroutine read_sputter_data(nsd,column,sput_model,refl_model)
      use b2mod_b2cmpa
      use b2mod_elements
      implicit none
      integer nsd, column, sput_model, refl_model
      integer nfile
#ifndef WINDOWS
      parameter (nfile=97)
#else
      parameter (nfile=87)
#endif
      character*256 filename, filename2, filename3, filename4, line
      character*2 adum
      integer idum1
      real (kind=R8) :: dummy
      integer ibnd, iss, isp, is, il, ic, i, itargsp, nmax
      logical open_file, found
      external open_file, ipgeti
      real (kind=R8), save  :: angles(9)
      data angles /0.0_R8, 15.0_R8, 30.0_R8, 45.0_R8, 55.0_R8,
     &   65.0_R8, 75.0_R8, 80.0_R8, 85.0_R8/
#ifndef PSR
      character*2 sputter_codes(0:nsd-1+6)
#else
      character*2 sputter_codes(0:100)
#endif
      external subini, subend

      if(sputter_initialised) return
      call subini ('read_sputter_data')
      allocate(sputter_index(0:nsd-1))
      sputter_index=-1
      sputter_codes='  '
      isp=-1
      do is=0,nsd-1
        do il=0,isp
          if(is_codes(is).eq.sputter_codes(il)) then
            sputter_index(is)=il
            goto 5
          endif
        enddo
        isp=isp+1
#ifdef PSR
        call xertst(isp.le.100, 
     &   'alloc_b2mod_sputter (1): increase size of sputter_codes')
#endif
        sputter_index(is)=isp
        sputter_codes(isp)=is_codes(is)
    5   continue
        write(*,'(1x,i3,1x,i3,1x,a)') is,sputter_index(is),
     1   sputter_codes(sputter_index(is))
      enddo
      nmax = 0
      zwall = 0
      do ibnd = 1, nbnd
        do is = 1, matsurf_nconstituents(ibnd)
          found = .false.
          iss = 0
          do while (.not.found .and. iss.lt.nsd)
            if (is_neutral(iss)) then
              found = nint(zn(iss)).eq.matsurf_nzconstituents(ibnd,is)
              iss = iss + 1
            else
              iss = iss + nint(zn(iss))
            endif
          enddo
          iss = 1
          do while (.not.found .and. iss.le.nmax)
            found = zwall(iss).eq.matsurf_nzconstituents(ibnd,is)
            iss = iss + 1
          enddo
          if (.not.found) then
            nmax = nmax + 1
            call xertst (nmax.le.6, 
     .       'increase size of zwall in b2mod_sputter')
            zwall(nmax) = matsurf_nzconstituents(ibnd,is)
#ifdef PSR
            call xertst(ntargsp+nmax.le.100,
     &       'alloc_b2mod_sputter (2): increase size of sputter_codes')
#endif
            target_code(ntargsp+nmax) = elements(zwall(nmax))
            target_index(zwall(nmax)) = ntargsp+nmax
          endif
        enddo
      enddo
      allocate(sputter_coefficient(30,0:isp,ntargsp+nmax))
      allocate(sputter_energy(30,0:isp,ntargsp+nmax))
      allocate(sputter_length(0:isp,ntargsp+nmax))
      allocate(sputter_coefficient2(30,0:isp,ntargsp+nmax))
      allocate(sputter_energy2(30,0:isp,ntargsp+nmax))
      allocate(sputter_length2(0:isp,ntargsp+nmax))
      allocate(reflect_coefficient(30,0:isp,ntargsp+nmax))
      allocate(reflect_energy(30,0:isp,ntargsp+nmax))
      allocate(reflect_length(0:isp,ntargsp+nmax))
      allocate(reflect_coefficient2(30,0:isp,ntargsp+nmax))
      allocate(reflect_energy2(30,0:isp,ntargsp+nmax))
      allocate(reflect_length2(0:isp,ntargsp+nmax))
      allocate(tgdns(ntargsp+nmax))
      allocate(Es(ntargsp+nmax))
      allocate(theta(ntargsp+nmax))
      do itargsp=1,ntargsp+nmax
        do is=0,isp
          filename = 
     1     sputter_codes(is)(1:len_trim(sputter_codes(is)))//
     2     target_code(itargsp)(1:len_trim(target_code(itargsp)))//'.y'
          call lower_case(filename)
          filename2 = 
     1     sputter_codes(is)(1:len_trim(sputter_codes(is)))//
     2     target_code(itargsp)(1:len_trim(target_code(itargsp)))//'.ye'
          call lower_case(filename2)
          if(sput_model.eq.1.and.open_file(99,filename,
     1     'data.local/trim.data/sputter.data|' // 
     2     'modules/B2.5/Database/trim.data/sputter.data')) then
            read(99,*)
            read(99,'(a80)') line
            read(line(index(line,'Es=')+3:index(line,'Es=')+7),*)
     1       Es(itargsp)
            write(*,*) 'Es read in from TRIM data file =', Es(itargsp)
            read(line(index(line,'rho=')+4:index(line,'rho=')+8),*)
     1       tgdns(itargsp)
            write(*,*) 'tgdns read in from TRIM data file =',
     1       tgdns(itargsp)
            read(99,*)
            read(99,*)
            read(99,*)(angles(i),i=1,9)
            theta(itargsp) = angles(column)
            write(*,*) 'You have chosen an incident angle of =',
     1       theta(itargsp)
            read(99,*)
            ic=0
   10       continue
            read(99,*,end=20) sputter_energy(ic+1,is,itargsp),
     1       (dummy,i=1,column-1),sputter_coefficient(ic+1,is,itargsp)
            write(*,*) sputter_energy(ic+1,is,itargsp),
     1       sputter_coefficient(ic+1,is,itargsp)
            sputter_energy(ic+1,is,itargsp)=
     1       sputter_energy(ic+1,is,itargsp)*ev
            ic=ic+1
            call xertst (ic.lt.30, 'increase 30 in alloc_b2mod_sputter')
            goto 10
   20       continue
            close(99)
            sputter_length(is,itargsp)=ic
            if(open_file(nfile,filename2,
     1       'data.local/trim.data/sputter.data|' // 
     2       'modules/B2.5/Database/trim.data/sputter.data')) then
              do i=1,6
                read(nfile,*)
              enddo
              ic=0
   30         continue
              read(nfile,*,end=40) sputter_energy2(ic+1,is,itargsp),
     1         (dummy,i=1,column-1),
     2         sputter_coefficient2(ic+1,is,itargsp)
              write(*,*) sputter_energy2(ic+1,is,itargsp),
     1         sputter_coefficient2(ic+1,is,itargsp)
              sputter_energy2(ic+1,is,itargsp)=
     1         sputter_energy2(ic+1,is,itargsp)*ev
              ic=ic+1
              call xertst (ic.lt.30,
     .         'increase 30 in alloc_b2mod_sputter')
              goto 30
   40         continue
              close(nfile)
              sputter_length2(is,itargsp)=ic
            else
              sputter_length2(is,itargsp)=0
            endif
          else
            if (sput_model.eq.1) write(*,*) 'No sputtering data for ',
     1       sputter_codes(is),' onto ',target_code(itargsp)
            if(open_file(98,'incatm',
     .                      'data.local|modules/B2.5/Database')) then
              do ic=1,94
                read(98,*) adum, idum1, atomic_mass(idum1), 
     1           surface_binding_energy(idum1), target_density(idum1)
              enddo
              close(98)
              if (itargsp.le.ntargsp) then
                tgdns(itargsp) =
     1           target_density(nint(zn(target_intcode(itargsp))))
                es(itargsp) = surface_binding_energy(
     1                          nint(zn(target_intcode(itargsp))))
              else
                tgdns(itargsp) = target_density(zwall(itargsp-ntargsp))
                es(itargsp) = 
     1           surface_binding_energy(zwall(itargsp-ntargsp))
              endif
              write(*,*) 'tgdns read in from "incatm" data file =',
     1         tgdns(itargsp)
              write(*,*) 'Es read in from "incatm" data file =', 
     1         Es(itargsp)
              write(*,*) 'Using empirical formula for ',
     1         sputter_codes(is),' onto ', target_code(itargsp)
              theta(itargsp) = angles(column)
              sputter_length(is,itargsp)=0
              sputter_length2(is,itargsp)=0
            else
              if (sput_model.ne.1) write(*,*) 'No sputtering data for ',
     1         sputter_codes(is),' onto ',target_code(itargsp)
              sputter_length(is,itargsp)=-1
              sputter_length2(is,itargsp)=0
            endif
          endif
c         write(*,*) filename
c         write(*,*)
c     1    (sputter_energy(i,is,itargsp)/ev,
c     1    sputter_coefficient(i,is,itargsp),
c     2    i=1,sputter_length(is,itargsp))
        enddo ! is
      enddo ! itargsp
      do itargsp=1,ntargsp+nmax
        do is=0,isp
          filename3 = 
     1     sputter_codes(is)(1:len_trim(sputter_codes(is)))//
     2     target_code(itargsp)(1:len_trim(target_code(itargsp)))//'.rn'
          call lower_case(filename3)
          filename4 = 
     1     sputter_codes(is)(1:len_trim(sputter_codes(is)))//
     2     target_code(itargsp)(1:len_trim(target_code(itargsp)))//'.re'
          call lower_case(filename4)
          if(refl_model.eq.1.and.open_file(96,filename3,
     1     'data.local/trim.data/refl.data|' // 
     2     'modules/B2.5/Database/trim.data/refl.data')) then
            read(96,*)
            if (itargsp.gt.ntargsp) then
              read(96,'(a80)') line
              read(line(index(line,'Es=')+3:index(line,'Es=')+7),*)
     1         Es(itargsp)
              write(*,*) 'Es read in from TRIM data file =', Es(itargsp)
              read(line(index(line,'rho=')+4:index(line,'rho=')+8),*)
     1         tgdns(itargsp)
              read(96,*)
              read(96,*)
              read(96,*)(angles(i),i=1,9)
              theta(itargsp) = angles(column)
              write(*,*) 'You have chosen an incident angle of =',
     1         theta(itargsp)
              read(96,*)
            else
              read(96,*)
              read(96,*)
              read(96,*)
              read(96,*)
              read(96,*)
            endif
            ic=0
   50       continue
            read(96,*,end=60) reflect_energy(ic+1,is,itargsp),
     1       (dummy,i=1,column-1),reflect_coefficient(ic+1,is,itargsp)
            write(*,*) reflect_energy(ic+1,is,itargsp),
     1       reflect_coefficient(ic+1,is,itargsp)
            reflect_energy(ic+1,is,itargsp)=
     1       reflect_energy(ic+1,is,itargsp)*ev
            ic=ic+1
            call xertst (ic.lt.30, 'increase 30 in alloc_b2mod_sputter')
            goto 50
   60       continue
            close(96)
            reflect_length(is,itargsp)=ic
            if(open_file(95,filename4,
     1       'data.local/trim.data/refl.data|' // 
     2       'modules/B2.5/Database/trim.data/refl.data')) then
              do i=1,6
                read(95,*)
              enddo
              ic=0
   70         continue
              read(95,*,end=80) reflect_energy2(ic+1,is,itargsp),
     1         (dummy,i=1,column-1),
     2         reflect_coefficient2(ic+1,is,itargsp)
              write(*,*) reflect_energy2(ic+1,is,itargsp),
     1         reflect_coefficient2(ic+1,is,itargsp)
              reflect_energy2(ic+1,is,itargsp)=
     1         reflect_energy2(ic+1,is,itargsp)*ev
              ic=ic+1
              call xertst (ic.lt.30,
     .         'increase 30 in alloc_b2mod_sputter')
              goto 70
   80         continue
              close(95)
              reflect_length2(is,itargsp)=ic
            else
              reflect_length2(is,itargsp)=0
            endif
          else
            if (refl_model.eq.1) write(*,*) 'No reflection data for ',
     1       sputter_codes(is),' onto ',target_code(itargsp)
            write(*,*) 'Using empirical formula for ',
     1       sputter_codes(is),' onto ',target_code(itargsp)
            theta(itargsp) = angles(column)
            reflect_length(is,itargsp)=0
            reflect_length2(is,itargsp)=0
          endif
c         write(*,*) filename3
c         write(*,*)
c     1    (reflect_energy(i,is,itargsp)/ev,
c     1    reflect_coefficient(i,is,itargsp),
c     2    i=1,reflect_length(is,itargsp))
        enddo       ! iss
        if (itargsp.le.ntargsp) then
          write(*,*) am(target_intcode(itargsp)),tgdns(itargsp)
          monolayer_l(itargsp)=
     &       (am(target_intcode(itargsp))/tgdns(itargsp)/avogr/1.0e6_R8)
     &       **(1.0_R8/3.0_R8)
          write(*,*) 'Monolayer thickness for ',target_code(itargsp),
     &       ' is ',monolayer_l(itargsp)
        else
          write(*,*) atomic_mass(zwall(itargsp-ntargsp)),tgdns(itargsp)
        endif
      enddo      ! itargsp
      if (ntargsp.eq.0) write(*,*) 'Sputtering models turned off!'
      sputter_initialised=.true.
      call subend ()
      return

      end subroutine read_sputter_data

      subroutine dealloc_sputter_data
      if(.not.sputter_initialised) return
      deallocate(sputter_index)
      deallocate(sputter_coefficient, sputter_energy, sputter_length)
      deallocate(sputter_coefficient2, sputter_energy2, sputter_length2)
      deallocate(reflect_coefficient, reflect_energy, reflect_length)
      deallocate(reflect_coefficient2, reflect_energy2, reflect_length2)
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(sputter_yield, sputter_yield2)
#endif
      sputter_initialised=.false.
      end subroutine dealloc_sputter_data

      function sputter(is,zb,issput,te)
      use b2mod_b2cmpa
      implicit none
      integer is, issput
      real (kind=R8) ::
     *   sputter, zb, te, yldphy_lowlim, yldphy_hilim, yldphy_te, f
      integer iss, ib, it, ih, lng, itargsp
c      write(*,*) 'sputter called',is,te
      iss=sputter_index(is)
      if(issput.eq.target_intcode(1)) then
        itargsp=1
      else if(ntargsp.ge.2.and.issput.eq.target_intcode(2)) then
        itargsp=2
      else if(ntargsp.ge.3.and.issput.eq.target_intcode(3)) then
        itargsp=3
      else
        write(0,*) 'ntargsp,issput,target_intcode', 
     *              ntargsp,issput,target_intcode(ntargsp)
        call xerrab('sputter: wrong issput')
      endif
      lng=sputter_length(iss,itargsp)
      if(lng.eq.0) then
        call physput(te/ev,zb,zn(is),zn(issput),
     *    am(is),am(issput),itargsp,sputter)
        return
      endif
      if(lng.eq.-1.or.te.le.0.0_R8) then
        sputter = 0.0_R8
        return
      endif

      if(te.le.sputter_energy(1,iss,itargsp)) then
c       write(*,*) 'te <= lower limit'
        call physput(sputter_energy(1,iss,itargsp)/ev,zb,zn(is),
     &   zn(issput),am(is),am(issput),itargsp,yldphy_lowlim)
        call physput(te/ev,zb,zn(is),
     &   zn(issput),am(is),am(issput),itargsp,yldphy_te)
        if (yldphy_lowlim.gt.0.0_R8.and.Yldphy_te.gt.0.0_R8) then
          sputter=sputter_coefficient(1,iss,itargsp)/
     1     yldphy_lowlim*yldphy_te
        else
          sputter = 0.0_R8
        endif
        return
      elseif(te.ge.sputter_energy(lng,iss,itargsp)) then
c       write(*,*) 'te >= upper limit'
        call physput(sputter_energy(lng,iss,itargsp)/ev,zb,zn(is),
     &   zn(issput),am(is),am(issput),itargsp,yldphy_hilim)
        call physput(te/ev,zb,zn(is),
     &   zn(issput),am(is),am(issput),itargsp,yldphy_te)
        sputter=sputter_coefficient(lng,iss,itargsp)/
     1   yldphy_hilim*yldphy_te
        return
      else
        ib=1
        it=lng
        do while (it-ib.gt.1)    ! new
c   10  ih=(ib+it)/2
          ih=(ib+it)/2
          if(te.lt.sputter_energy(ih,iss,itargsp)) then
            it=ih
          elseif(te.gt.sputter_energy(ih,iss,itargsp)) then
            ib=ih
          else
            sputter=sputter_coefficient(ih,iss,itargsp)
            return
          endif
c       write(*,*) it-ib
c       if(it-ib.gt.1) goto 10
        enddo    ! new
c       write(*,*) sputter_energy(ib,iss,itargsp),te,
c                  sputter_energy(it,iss,itargsp)
        f=(te-sputter_energy(ib,iss,itargsp))/
     1       (sputter_energy(it,iss,itargsp)-
     2        sputter_energy(ib,iss,itargsp))
        if(jagged_rates.eq.0) then
          sputter=f*
     1     (sputter_coefficient(it,iss,itargsp)-
     2     sputter_coefficient(ib,iss,itargsp))+
     1     sputter_coefficient(ib,iss,itargsp)
        else
          sputter=f*
     1     (sputter_coefficient(it,iss,itargsp)-
     2     sputter_coefficient(ib,iss,itargsp))+
     1     sputter_coefficient(it,iss,itargsp)
        endif
        return
      endif
      end function sputter
c
c  Function to evaluate the sputtered energy from TRIM data.
      function sputter2(is,issput,te)
      use b2mod_b2cmpa
      implicit none
      real (kind=R8) ::
     *   sputter2, te, f
      integer is, issput, iss, ib, it, ih, lng, itargsp
c      write(*,*) 'sputter2 called',is,te
      iss=sputter_index(is)
      if(issput.eq.target_intcode(1)) then
        itargsp=1
      else if(ntargsp.ge.2.and.issput.eq.target_intcode(2)) then
        itargsp=2
      else if(ntargsp.ge.3.and.issput.eq.target_intcode(3)) then
        itargsp=3
      else
        call xerrab('sputter2: wrong issput')
      endif
      lng=sputter_length2(iss,itargsp)
      if(lng.eq.0) then
        sputter2 = 0.0_R8
        return
      endif

      if(te.le.sputter_energy2(1,iss,itargsp)) then
c       write(*,*) 'te <= lower limit'
        sputter2 = 0.0_R8
        return
      elseif(te.ge.sputter_energy2(lng,iss,itargsp)) then
c       write(*,*) 'te >= upper limit'
        sputter2=sputter_coefficient2(lng,iss,itargsp)
        return
      else
        ib=1
        it=lng
        do while (it-ib.gt.1)    ! new
c   10  ih=(ib+it)/2
          ih=(ib+it)/2
          if(te.lt.sputter_energy2(ih,iss,itargsp)) then
            it=ih
          elseif(te.gt.sputter_energy2(ih,iss,itargsp)) then
            ib=ih
          else
            sputter2=sputter_coefficient2(ih,iss,itargsp)
            return
          endif
c       write(*,*) it-ib
c       if(it-ib.gt.1) goto 10
        enddo    ! new
c       write(*,*) sputter_energy2(ib,iss,itargsp),te,
c     @     sputter_energy2(it,iss,itargsp)
        f=(te-sputter_energy2(ib,iss,itargsp))/
     1       (sputter_energy2(it,iss,itargsp)-
     2        sputter_energy2(ib,iss,itargsp))
        if(jagged_rates.eq.0) then
          sputter2=f*(sputter_coefficient2(it,iss,itargsp)-
     1     sputter_coefficient2(ib,iss,itargsp))+
     1     sputter_coefficient2(ib,iss,itargsp)
        else
          sputter2=f*(sputter_coefficient2(it,iss,itargsp)-
     1     sputter_coefficient2(ib,iss,itargsp))+
     1     sputter_coefficient2(it,iss,itargsp)
        endif
        return
      endif
      end function sputter2

      function reflect(is,te,ibnd)
      use b2mod_b2cmpa
      implicit none
      integer is, ibnd
      real (kind=R8) :: te, reflect
      real (kind=R8) ::
     *   reffrac, reflrn_hilim, reflrn_lolim, reflrn_te, f
      integer iss, ib, it, ih, lng
      integer isp, itargsp
c     write(*,*) 'reflect called',is,te/ev
      if(.not.allocated(sputter_index)) then
        reflect=0.0_R8
        return
      endif
      iss=sputter_index(is)
      reflect=0.0_R8
      do isp = 1, layer_nconstituents(ibnd)
        itargsp = target_index(layer_nzconstituents(ibnd,isp))
        lng=reflect_length(iss,itargsp)
        if(lng.eq.-1.or.te.le.0.0_R8) then  ! no data
          reffrac = 0.0_R8
        elseif(lng.eq.0) then  ! analytic formula
          reffrac = rnion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1     theta(itargsp))
        elseif(te.le.reflect_energy(1,iss,itargsp)) then  ! TRIM low Te extrapolation
c         write(*,*) 'te <= lower limit'
          reflrn_lolim=rnion(reflect_energy(1,iss,itargsp)/ev,
     &     nint(zn(is)),nint(am(is)),ibnd,theta(itargsp))
          reflrn_te=rnion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1     theta(itargsp))
          reffrac=reflect_coefficient(1,iss,itargsp)/
     1     reflrn_lolim*reflrn_te
        elseif(te.ge.reflect_energy(lng,iss,itargsp)) then  ! TRIM high Te extrapolation
c         write(*,*) 'te >= upper limit'
          reflrn_hilim=rnion(reflect_energy(lng,iss,itargsp)/ev,
     &     nint(zn(is)),nint(am(is)),ibnd,theta(itargsp))
          reflrn_te=rnion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1     theta(itargsp))
          reffrac=reflect_coefficient(lng,iss,itargsp)/
     1     reflrn_hilim*reflrn_te
        else  ! TRIM data interpolation
          ib=1
          it=lng
          do while (it-ib.gt.1)
            ih=(ib+it)/2
            if(te.lt.reflect_energy(ih,iss,itargsp)) then
              it=ih
            elseif(te.gt.reflect_energy(ih,iss,itargsp)) then
              ib=ih
            else
              reffrac=reflect_coefficient(ih,iss,itargsp)
              return
            endif
          enddo
c         write(*,*) 
c    1     reflect_energy(ib,iss,itargsp)/ev,te/ev,
c    2     reflect_energy(it,iss,itargsp)/ev
          f=(te-reflect_energy(ib,iss,itargsp))/
     1         (reflect_energy(it,iss,itargsp)-
     2          reflect_energy(ib,iss,itargsp))
          if(jagged_rates.eq.0) then
            reffrac=f*
     1       (reflect_coefficient(it,iss,itargsp)-
     2        reflect_coefficient(ib,iss,itargsp))+
     1        reflect_coefficient(ib,iss,itargsp)
          else
            reffrac=f*
     1       (reflect_coefficient(it,iss,itargsp)-
     2        reflect_coefficient(ib,iss,itargsp))+
     1        reflect_coefficient(it,iss,itargsp)
          endif
        endif
        reflect=reflect+reffrac*layer_nrelconstituents(ibnd,isp)
      enddo
      return
      end function reflect
c
c  Function to evaluate the reflected energy from TRIM data.
      function reflect2(is,te,ibnd)
      use b2mod_b2cmpa
      implicit none
      integer is, ibnd
      real (kind=R8) :: te, reflect2
      real (kind=R8) :: 
     *   reffrac, f, reflre_hilim, reflre_lolim, reflre_te
      integer iss, ib, it, ih, lng, isp, itargsp
c     write(*,*) 'reflect2 called',is,te/ev
      if(.not.allocated(sputter_index)) then
        reflect2=0.0_R8
        return
      endif
      iss=sputter_index(is)
      reflect2=0.0_R8
      do isp = 1, layer_nconstituents(ibnd)
        itargsp = target_index(layer_nzconstituents(ibnd,isp))
        lng=reflect_length2(iss,itargsp)
        if(lng.eq.-1.or.te.le.0.0_R8) then  ! no data
          reffrac = 0.0_R8
        elseif(lng.eq.0) then  ! analytic formula
          reffrac = reion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1        theta(itargsp))
        elseif(te.le.reflect_energy2(1,iss,itargsp)) then  ! TRIM low Te extrapolation
c       write(*,*) 'te <= lower limit'
          reflre_lolim=reion(reflect_energy2(1,iss,itargsp)/ev,
     &     nint(zn(is)),nint(am(is)),ibnd,theta(itargsp))
          reflre_te=reion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1     theta(itargsp))
          reffrac=reflect_coefficient2(1,iss,itargsp)/
     2     reflre_lolim*reflre_te
        elseif(te.ge.reflect_energy2(lng,iss,itargsp)) then  ! TRIM high Te extrapolation
c       write(*,*) 'te >= upper limit'
          reflre_hilim=reion(reflect_energy2(lng,iss,itargsp)/ev,
     &     nint(zn(is)),nint(am(is)),ibnd,theta(itargsp))
          reflre_te=reion(te/ev,nint(zn(is)),nint(am(is)),ibnd,
     1     theta(itargsp))
          reffrac=reflect_coefficient2(lng,iss,itargsp)/
     2     reflre_hilim*reflre_te
        else  ! TRIM data interpolation
          ib=1
          it=lng
          do while (it-ib.gt.1)
            ih=(ib+it)/2
            if(te.lt.reflect_energy2(ih,iss,itargsp)) then
              it=ih
            elseif(te.gt.reflect_energy2(ih,iss,itargsp)) then
              ib=ih
            else
              reffrac=reflect_coefficient2(ih,iss,itargsp)
              return
            endif
          enddo
c       write(*,*) reflect_energy2(ib,iss,itargsp)/ev,te/ev,
c     @     reflect_energy2(it,iss,itargsp)/ev
          f=(te-reflect_energy2(ib,iss,itargsp))/
     1       (reflect_energy2(it,iss,itargsp)-
     2        reflect_energy2(ib,iss,itargsp))
          if(jagged_rates.eq.0) then
            reffrac=f*(reflect_coefficient2(it,iss,itargsp)-
     1       reflect_coefficient2(ib,iss,itargsp))+
     1       reflect_coefficient2(ib,iss,itargsp)
          else
            reffrac=f*(reflect_coefficient2(it,iss,itargsp)-
     1       reflect_coefficient2(ib,iss,itargsp))+
     1       reflect_coefficient2(it,iss,itargsp)
          endif
        endif
        reflect2=reflect2+reffrac*layer_nrelconstituents(ibnd,isp)
      enddo
      return
      end function reflect2

c  Subroutine to evaluate sputtering yields of
c  any target material due to physical sputtering.

c  NOTE : These formulae are valid only for normal incidence
c         and have to be modified for angular incidence.

c  REF-1: "Erosion processes in plasma-wall interactions"
c         C. Garcia Rosales
c         Jnl. Nucl. Mat. 211 (1994) 202-214

      subroutine physput(eo,zb,z1,z2,am1,am2,itargsp,yldphy)
      implicit none

      real (kind=R8) :: eo, zb, z1, z2, am1, am2, eth, etf, 
     @  stoppwr, yieldpar, ethbyeo, yldphy,
     @  pwr1by6, pwr5by6, pwr1by3, pwr2by3,
     @  z123xz223, z123z223sum, am2byam1,
     @  ypar1, ypar2, ypar3, eobyetf,
     @  stoppwr1, stoppwr2
      real (kind=R8) :: Ro, a, f, fs, q, eps, anu, ang_opt,
     .  eta, psi, sigma, t, angcntrb, gam
      integer iionflg, itargsp
 
c  Inputs to the subroutine
c  eo   -> energy of the incident particle.
c  z1   -> nuclear charge of incident atom.
c  z2   -> nuclear charge of target atom.
c  am1  -> atomic mass of incident atom (amu).
c  am2  -> atomic mass of target atom (amu).

c  es, theta and tgdns are passed as a module global variable
c  es   -> surface binding energy (heat of sublimation) of target (eV).
c  theta -> angle (in degrees) with normal to target of incident particle.
c  tgdns -> target density (gms/cc)

c  Output from the subroutine
c  yldphy -> The physical sputtering yield.

      if (eo.le.0.0_R8) then
        yldphy = 0.0_R8
        return
      endif

c  Initialising

      pwr1by6 = 1.0_R8 / 6.0_R8
      pwr5by6 = 5.0_R8 / 6.0_R8
      pwr1by3 = 1.0_R8 / 3.0_R8
      pwr2by3 = 2.0_R8 / 3.0_R8

      z123xz223 = z1**pwr2by3 * z2**pwr2by3
      z123z223sum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

      ypar1 = 1.633_R8 * z123xz223 * z123z223sum**pwr1by3 /
     @  es(itargsp)**pwr2by3
      ypar2 = am1**pwr5by6 * am2**pwr1by6 / (am1 + am2)
      ypar3 = (0.15_R8 + 0.05_R8*am2byam1) /
     @  (1.0_R8 + 0.05_R8*am2byam1**1.6_R8)
      yieldpar = ypar1 * ypar2 * ypar3

      eth = (7.0_R8/am2byam1**0.54_R8+
     @  0.15_R8*am2byam1**1.12_R8)*es(itargsp)
      ethbyeo = eth / eo

c  The Thomas Fermi energy etf used below is valid only for singly
c  charged ions. If it is multiply charged just multiply etf by the
c  charge of the ion to get the etf for multiply charged ions.

      etf = 30.74_R8*(am1+am2)/am2 *z1*z2*z123z223sum**0.5_R8
     @  *max(1.0_R8,zb)
c       print*,z1,' ',am1,' ',z2,' ',am2,' ',etf

      eobyetf = eo / etf
      stoppwr1=(0.5_R8*log(1.0_R8+1.2288_R8*eobyetf))
      stoppwr2=(eobyetf + 0.1728_R8*eobyetf**0.5_R8 +
     @  0.008_R8*eobyetf**0.1504_R8)
      stoppwr = stoppwr1 / stoppwr2

c  The angular dependence of the sputtering yield is included for light ion
c  incidence. The yield is valid for heavy ion scattering too, but for that
c  you have to dig out the "f" and "angl_opt" from
c  Ref-2: Institute of Plasma Physics, Nagoya University report - IPPJ-AM-26,
c         Japan.
c  ** Note that I have used etf for E_L (eqn 10 of above reference). e**2
c     is not defined in the report .. at least I did not find it. I initially
c     thought e was the charge of incident ion and used 1, i.e. singly ionised
c     ion, but that gave wrong results. Using etf solved the problem.

c  Checking whether it is light ion sputtering or heavy ion sputtering.
c  iionflg -> flag for light/heavy ion.
c  iionflg = 0 => light ion sputtering.
c  iionflg = 1 => heavy ion sputtering.
      if (am1.le.4.0) then
        iionflg = 0
      else
        iionflg = 1
      endif

      Ro = (am2/tgdns(itargsp)/avogr)**(1.0_R8/3.0_R8)*1.0e8_R8
      a = 0.4685_R8*sqrt(1.0_R8/z123z223sum)
      if (iionflg.eq.0) then
        f = (0.94_R8 - 1.33e-3_R8 * am2byam1) * sqrt(es(itargsp))
        gam = 4.0_R8*am1*am2/(am1+am2)**2
        q = sqrt(es(itargsp)/gam/eo)
        eps = eo/etf
        anu = a/Ro*sqrt(1.0_R8/2.0_R8/eps/q)
        ang_opt = 90.0_R8 - 180.0_R8*anu/pi
      else
c         For a really good quantitative study you have to input the
c         value of fs from Ref-2 below
c         I am circumventing the need of making a big database
c         by using fig-3 from Ref-2.
        if (am2byam1.le.2.0_R8) fs = 1.75_R8
        if (2.0_R8.le.am2byam1.and.am2byam1.lt.4.0_R8) fs=1.5_R8
        if (4.0_R8.le.am2byam1.and.am2byam1.lt.6.0_R8) fs=1.25_R8
        if (6.0_R8.le.am2byam1.and.am2byam1.lt.8.0_R8) fs=1.0_R8
        if (am2byam1.ge.8.0_R8) fs = 0.8_R8
        eta = 1.0_R8 - sqrt(eth/eo)
        if (eta.gt.0.0_R8) then
          f = min(10.0_R8,fs*(1.0_R8+2.5_R8*(1.0_R8-eta)/eta))
        else
          f = fs
        endif
        psi = (a/Ro)**(3.0_R8/2.0_R8)*
     @    sqrt(z1*z2/sqrt(z123z223sum)/eo)
        ang_opt = 90.0_R8-286.0_R8*psi**0.45_R8
      endif
      Sigma = cos(ang_opt*pi/180.0_R8)
      t = 1.0_R8/cos(theta(itargsp)*pi/180.0_R8)
      angcntrb = exp(f*(Sigma*(1.0_R8-t)+log(t)))
c  Calculation of the physical sputtering yield
      yldphy = yieldpar * stoppwr *
     @  (1.0_R8-ethbyeo**pwr2by3)
     @  * ( 1.0_R8 - ethbyeo ) ** 2 *
     @  angcntrb
      if (yldphy.lt.0.0e0_R8) then
c       write(*,*) 'Physical sputtering yield < 0 , Ei = ',eo
        yldphy = 0.0e0_R8
      endif

c     print*,yldphy,' ',angcntrb,' ',etf,' ',stoppwr,' ',eo
      return
      end subroutine physput


c******************************************************************

c     SUBROUTINE  chmsput

c**********************************************************************

c  Inputs to the subroutine :
c  --------------------------
c  eo -> energy of incident projectile (eV)
c  z1   -> nuclear charge of incident atom.
c  z2   -> nuclear charge of target atom.
c  am1  -> atomic mass of incident atom (amu).
c  am2  -> atomic mass of target atom (amu).
c  tp   -> temperature of the target plate (K).
c  gamma-> incident ion flux (part/m**2/sec).

c  Es is passed as part of the module global variables
c  es   -> surface binding energy (heat of sublimation) of target (eV).

c  Outputs from the subroutine :
c  -----------------------------
c  yldchm -> Chemical sputtering yield.


c  Subroutine to evaluate sputtering yields. 

      subroutine chmsput(ibnd,eo,zb,z1,z2,am1,am2,tp,gamma,yldchm)
      implicit none

      integer ibnd
      real (kind=R8) :: eo, zb, z1, z2, am1, am2, tp, gamma, 
     @  etf, stoppwr, yieldpar, yldchm
      real (kind=R8) :: pwr2by3, z123z223sum, am2byam1, tpev, eobyetf,
     .  stoppwr1, stoppwr2, flxm2, Cmy, etherm, erel, edam, edes,
     .  Csp3num, Csp3denom, Csp3, ytherm, D, Chd,
     .  ydam, ydes, ysurf, yldthmprt
      integer, save :: formula
      data formula/0/
 
      if (eo.le.0.0_R8 .or. gamma.le.0.0_R8) then
        yldchm = 0.0_R8
        return
      endif

c  Initialising
      pwr2by3 = 2.0e0_R8 / 3.0e0_R8

      z123z223sum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

c  The Thomas Fermi energy etf used below is valid only for singly
c  charged ions. If it is multiply charged just multiply etf by the
c  charge of the ion to get the etf for multiply charged ions.
      etf = 30.74e0_R8 * (am1+am2)/am2  * z1*z2*z123z223sum**0.5e0_R8
     @  *max(1.0_R8,zb)
      tpev = tp / (ev/KBolt)

      eobyetf = eo / etf
      stoppwr1 = (0.5e0_R8*log(1.0e0_R8+1.2288e0_R8*eobyetf))
      stoppwr2 = (eobyetf + 0.1728e0_R8*eobyetf**0.5e0_R8 +
     @  0.008e0_R8*eobyetf**0.1504e0_R8)
      stoppwr = stoppwr1 / stoppwr2

      flxm2 = gamma
cxpb Modified to follow new formula from Roth et al., PSI16
c      flxlmt = 1.0e30_R8*exp(-1.4/tpev)
c      if (flxm2.gt.flxlmt) then  ! Assuming high flux
c        Cmy = 1.0e0_R8/(1.0e0_R8+3.0e-23_R8*flxm2)
c      else
c        Cmy = 1.0e0_R8/(1+3.0e7_R8*exp(-1.4/tpev))
c      endif

cxpb - 16.06.05 Corrected unphysical flux dependence behaviour pointed
cxpb - out by Detlev Reiter.
      if (formula.eq.1) then
C  Haasz-Davis formula, 1998, with flx. dep from Roth, Nucl.Fus 2004
        Chd=1.0e0_R8/(1.0e0_R8+(1.67e-22_R8*flxm2)**0.54_R8)
        yldchm=Chd * yhaasz97m(eo,tp)
      else

        Cmy = 1.0e0_R8/((1.0e0_R8+(flxm2/6.0e21_R8)**0.54_R8)
     .                 *(1.0e0_R8+3.0e7_R8*exp(-1.6_R8/tpev)))
        etherm = 1.7e0_R8
        erel = matsurf_chemsput_erel(ibnd)    ! erel =1.8 (for pure Graphite);
                                              !      =1.5 (for Si,Ti,W doped graphite)
                                              !      =1.2 (for B doped graphite)
        edam = 15.0e0_R8
        edes = 2.0e0_R8
        Csp3num = Cmy*(2.0e-32_R8*flxm2+exp(-etherm/tpev))
        Csp3denom = 2.0e-32_R8*flxm2+
     @    (1.0e0_R8+2.0e29_R8/flxm2*exp(-erel/tpev))*
     @    exp(-etherm/tpev)
        Csp3 = Csp3num/Csp3denom
        ytherm = Csp3*0.033*exp(-etherm/tpev) /
     @    (2.0e-32_R8*flxm2+exp(-etherm/tpev))
        if (nint(am1).eq.1) then
          yieldpar = 0.035
        elseif (nint(am1).eq.2) then
          yieldpar = 0.1
        elseif (nint(am1).eq.3) then
          yieldpar = 0.12
        endif
        D = 250.0e0_R8/am1
        ydam = yieldpar*stoppwr*(1.0e0_R8-(edam/eo)**pwr2by3)*
     @    (1.0e0_R8-edam/eo)**2.0e0_R8
        ydes = yieldpar*stoppwr*(1.0e0_R8-(edes/eo)**pwr2by3)*
     @    (1.0e0_R8-edes/eo)**2.0e0_R8
        ysurf = csp3*ydes/(1.0e0_R8+exp((eo-65.0e0_R8)/40.0e0_R8))
        yldthmprt = ytherm*(1.0e0_R8+D*ydam)
        yldchm = yldthmprt+ysurf
      endif ! Haasz-Davis or Roth formula
      if (yldchm.le.0.0_R8) then
c       write(*,*) 'Chemical sputtering yield < 0, Ei = ',eo
        yldchm = 0.0_R8
      endif

      return
      end subroutine chmsput

c  Subroutine to evaluate sputtering yields of
c  Graphite due to Radiation Enhanced Sublimation.

c  NOTE : These formula are valid only for normal incidence
c         and have to be modified for angular incidence.

      subroutine ressput(eo,zb,z1,z2,am1,am2,tp,gamma_in,issput,
     @  yldres)

      implicit none
      real (kind=R8) :: eo, zb, z1, z2, am1, am2, tp, gamma_in, 
     @  yldres

      real (kind=R8) :: pwr1by6, pwr5by6, pwr1by3, pwr2by3, 
     @  z123xz223, z123z223sum, 
     @  tpev, am2byam1, stoppwr, stoppwr1, stoppwr2,
     @  yldparres, eobyetf, ethbyeo, eth, etf, gamma
        integer issput, itargsp
c  Inputs to the code

c  eo   -> incident ion energy (in eV).
c  z1   -> nuclear charge of incident atom.
c  z2   -> nuclear charge of target atom.
c  am1  -> atomic mass of incident atom (amu).
c  am2  -> atomic mass of target atom (amu).
c  tp   -> temperature of the target plate (eV).
c  gamma -> incident ion flux (part/m**2/sec).

c  es is passed as a module global variable
c  es   -> surface binding energy (heat of sublimation) of target (eV).

      if(issput.eq.target_intcode(1)) then
        itargsp=1
      else if(ntargsp.ge.2.and.issput.eq.target_intcode(2)) then
        itargsp=2
      else if(ntargsp.ge.3.and.issput.eq.target_intcode(3)) then
        itargsp=3
      else
        call xerrab('ressput: wrong issput')
      endif

      gamma = gamma_in
      if (eo.le.0.0_R8 .or. gamma.le.0.0_R8 .or. nint(z2).ne.6) then
        yldres = 0.0_R8
        return
      endif

      pwr1by6 = 1.0e0_R8 / 6.0e0_R8
      pwr5by6 = 5.0e0_R8 / 6.0e0_R8
      pwr1by3 = 1.0e0_R8 / 3.0e0_R8
      pwr2by3 = 2.0e0_R8 / 3.0e0_R8
      gamma = gamma / 1.0e4_R8

      z123xz223 = z1**pwr2by3 * z2**pwr2by3
      z123z223sum = z1**pwr2by3 + z2**pwr2by3
      am2byam1 = am2 / am1

      eth = ( 7.0e0_R8/am2byam1**0.54e0_R8 +
     @  0.15e0_R8*am2byam1**1.12e0_R8) * es(itargsp)

c  The Thomas Fermi energy etf used below is valid only for singly
c  charged ions. If it is multiply charged just multiply etf by the
c  charge of the ion to get the etf for multiply charged ions.

      etf = 30.74e0_R8 * (am1+am2)/am2  * z1*z2*
     @  z123z223sum**0.5e0_R8 *max(1.0_R8,zb)

      tpev = tp / (ev/KBolt)

      yldparres = 54.0e0_R8 * am1**1.18e0_R8 *
     @  exp((-0.78e0_R8)/tpev)*(gamma/1.0e16_R8)**(-0.1)

      ethbyeo = eth / eo
      eobyetf = eo / etf
      stoppwr1 = (0.5e0_R8*
     @  log(1.0e0_R8+1.2288e0_R8*eobyetf))
      stoppwr2 = (eobyetf+0.1728e0_R8*eobyetf**0.5e0_R8+
     @  0.008e0_R8*eobyetf**0.1504e0_R8)
      stoppwr = stoppwr1 / stoppwr2

c  Calculation of the physical sputtering yield

c  Calculation of the RES sputtering yield
c  assuming that the flux factor is unity.

      yldres = yldparres * stoppwr *
     @  (1.0e0_R8-ethbyeo**pwr2by3)
     @  * ( 1.0e0_R8 - ethbyeo ) ** 2

      gamma = gamma * 1.0e4_R8
      if (yldres.le.0.0_R8) then
c       write(*,*) 'Radiation Enhanced Sublimation < 0, Ei = ', eo
        yldres = 0.0_R8
      endif

      return
      end subroutine ressput


C   Subroutine to evaluate the thermal evaporation of Graphite.

c   Ref: 1) C. Garc'ia Rosales, Erosion processes in plasma-wall
c           interactions, J. Nucl. Mater. 211 (1994) Pg.202-214.
c        2) R.A. Langley, Evaporation, Data Compendium for Plasma
c           Surface Interactions, Special Issue, NUCLEAR FUSION,
c           IAEA, Vienna, (1984) 55-61.

c   Written by Manoj Warrier on 17-05-96.

*********************************************************************************

c  This routine takes the temperature of the graphite plate ( tp ), and the
c  atomic mass ( am2 ) of the target as inputs and outputs the flux of carbon
c  atoms evaporated from the graphite plate ( flxthev ) per meter**2 per sec.
c
c  Table of A and B for some first wall materials in Fusion research:
c  (Ref (2) above).
c  (note that A and B are provided such that the vappres is in Pascals.
c  element    |  A       |  B      |
c  Aluminum   | -16165   | 10.91   |
c  Beryllium  | -16720   | 11.61   |
c  Graphite   | -40181   | 14.80   |
c  Copper     | -17079   | 11.26   |
c  Molybdenum | -33069   | 12.02   |
c  Titanium   | -23340   | 9.79    |
c  Tungsten   | -44485   | 12.74   |
c  Titanium - |          |         |
c     Carbide | -39117.9 | 15.23   |
c  Stainless -|          |         |
c     Steel   | -14321.9 | 9.44    |

      subroutine thermev(ibnd,is,tp,flxthev,qthev)
      use b2mod_b2cmpa
      implicit none
      integer ibnd, is, i, itrack, nz
      real (kind=R8) :: tp, A, B, am2, vappres, alpha
      real (kind=R8) :: flxthev, qthev
      logical bulk

c Inputs:
c  tp       ->  temperature of the plate in Kelvin.
c  am2      ->  atomic mass of the target ( in amu )
c  A        ->  Parameter in typical Vapor pressure curve.
c  B        ->  Parameter in typical Vapor pressure curve.
c  alpha    ->  Sticking Coefficient ( = 0.6-0.9 for metals and
c               0.05 for C_3 molecules and 0.4 for C )

c Output:
c  flxthev  ->  Flux of target atoms evaporated ( /m**2/sec )
c  qthev    ->  Total heat of sublimation (Joules / m^2 / s)

c  vappres  ->  Vapour Pressure of the target in Pascals.

      itrack = track_index(is)
      flxthev = 0.0_R8
      qthev = 0.0_R8
      nz = matsurf_nconstituents(ibnd)

cxpb Contribution from bulk if species is part of bulk
      bulk = .false.
      do i = 1, nz
        if (zn(is).eq.matsurf_nzconstituents(ibnd,i)) then
          A = matsurf_vapprs_A(ibnd)
          B = matsurf_vapprs_B(ibnd)
          am2 = atomic_mass(matsurf_nzconstituents(ibnd,i))
          alpha = matsurf_stick_coef(ibnd)
          bulk = .true.
        end if
      end do
      
      if (bulk) then
        vappres = 10.0e0_R8**(A/tp + B)
        flxthev = alpha*vappres/sqrt(2.0_R8*pi*am2*mp*tp*KBolt)*
     .     (layer_nrelconstituents(ibnd,1)/
     .    matsurf_nrelconstituents(ibnd,1))
        if (flxthev.le.0.0_R8) then
c       write(*,*) 'Thermal evaporation < 0, Surface temp = ',tp
          flxthev = 0.0_R8
        endif
        qthev = flxthev * am(is) * mp * matsurf_heat_sublim(ibnd)
      endif

      if (.not.(itrack.gt.0.and.itrack.le.ntrack)) return   ! this is not a tracked species
      if (nz.eq.layer_nconstituents(ibnd)) return           ! we are not doing mixed materials

cxpb Contribution from re-deposited species
cxpb if using the mixed matarials model
cxpb Also need to consider alloy contributions
      A = species_vapprs_A(itrack)
      B = species_vapprs_B(itrack)
      alpha = species_stick_coef(itrack)
      am2 = atomic_mass(layer_nzconstituents(ibnd,nz+itrack))
      vappres = 10.0e0_R8**(A/tp + B)
      flxthev = flxthev+alpha*vappres/sqrt(2.0_R8*pi*am2*mp*tp*KBolt)*
     .  layer_nrelconstituents(ibnd,nz+itrack)
      qthev = qthev + flxthev * am(is) * mp * matsurf_heat_sublim(ibnd)

      return
      end subroutine thermev


c  Subroutine to find the surface temperature of a plasma facing
c  surface assuming a 1-D. steady state, heat diffusion with a
c  thermal conductivity that varies with temperature as
c   K(T) = 1/(KA * T + KB)

      subroutine heatdiff0D(ibnd,qin,qel,Tp)
      implicit none
      integer ibnd
      real (kind=R8) :: qin, qel, Qtot, TL, L, K, KB, KA, Tp, Tc, Td
      real (kind=R8), save :: cutll, cutlo
      logical, save  :: lfirst, coating
      real (kind=R8) :: t0
      real (kind=R8) :: machsfr, expu2, layer_Kappa
      external machsfr, layer_Kappa
      intrinsic dim, exp, min, max
      expu2(t0) = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
*     (expu2 is the exponential function with underflow set to 0 and
*     large values limited at a somewhat arbitrary level.)
*   ..initialisation
      data lfirst /.true./

c  Inputs:
c  qin -> incident ion heat flux (Joules/m^2/sec)
c  qel -> incident electron heat flux (Joules/m^2/sec)
c  TL -> Temperature at x=L (x=0 is plasma facing surface,
c                            x=L is the cold end).
c  L  -> Thickness of the target plate (m).
c  KB -> parameter in the expression for conductivity (see above)
c  KA -> parameter in the expression for conductivity (see above)

c  Outputs:
c  Tp -> Surface temperature (at x=0) (in Kelvins).

c  Used module variables:
c  layer_thickness : deposited - eroded
c  layer_depth : sum of deposition, i.e. minimum depth over which material mixing occurs

c  Intermediate variables:
c  Tc -> Temperature at interface between coating and bulk (in Kelvins).
c  Td -> Temperature at interface between mixing layer and original wall material.

c  Remarks:
c  Consider adding the option of layer_thermal_contact

*   ..set cutll, cutlo
      if (lfirst) then
        cutll = log(machsfr())
        cutlo = max(machsfr(),exp(cutll))
        lfirst = .false.
      end if

      TL = backplate_temp(ibnd)
      L = plate_thickness(ibnd)
      if (L.eq.0.0_R8) then
        Tp = TL
        return
      endif
      Qtot = qel + qin
      coating = coating_thickness(ibnd).gt.0.0_R8 .and.
     . (coating_thickness(ibnd)+layer_thickness(ibnd)).gt.0.0_R8
      if (.not. coating) then
        KA = matbulk_Kappa_KA(ibnd)
        KB = matbulk_Kappa_KB(ibnd)
        if (KA.ne.0.0_R8) then
          Td = max(TL, min( matsurf_vapor_temp(ibnd), (TL + KB/KA) * 
     .     expu2(KA*Qtot*(min(L,L+layer_thickness(ibnd)))) - KB/KA )) 
        else
          Td = max(TL, min( matsurf_vapor_temp(ibnd), 
     .     TL + min(L,L+layer_thickness(ibnd))*Qtot*KB))
        endif
      else if (coating) then
        KA = matbulk_Kappa_KA(ibnd)
        KB = matbulk_Kappa_KB(ibnd)
        if (KA.ne.0.0_R8) then
          Tc = max(TL, min( matsurf_vapor_temp(ibnd), (TL + KB/KA) * 
     .     expu2(KA*Qtot*(L-coating_thickness(ibnd))) - KB/KA ))
        else
          Tc = max(TL, min( matsurf_vapor_temp(ibnd), 
     .     TL + (L-coating_thickness(ibnd))*Qtot*KB))
        endif
        KA = coating_Kappa_KA(ibnd)
        KB = coating_Kappa_KB(ibnd)
        if (KA.ne.0.0_R8) then
          Td = max(Tc, min( matsurf_vapor_temp(ibnd), (Tc + KB/KA) * 
     .     expu2(KA*Qtot*min(coating_thickness(ibnd),
     .      coating_thickness(ibnd)+layer_thickness(ibnd)))
     .      - KB/KA )) 
        else
          Td = max(Tc, min( matsurf_vapor_temp(ibnd), Tc + 
     .     min(coating_thickness(ibnd),
     .         coating_thickness(ibnd)+layer_thickness(ibnd))*Qtot*KB))
        endif
      endif
      if (layer_depth(ibnd).gt.0.0_R8) then
        K = layer_Kappa(ibnd,Td)
        Tp = min(matsurf_vapor_temp(ibnd),Td+Qtot*layer_depth(ibnd)/K)
      else
        Tp = Td
      endif

      return
      end subroutine heatdiff0D

c  A FORTRAN FUNCTION to evaluate empirical formulae for number and
c  energy BACKSCATTERING COEFFICIENTS of light ions incident on
c  Elemental and Compound targets.
c
c***********************************************************************
c
c  program copied from
c
c  NAYOGA UNIVERSITY REPORT #  IPPJ-AM-41 ( Data On The Backscattering
c  Coefficients Of Light Ions From Solids ) - R.Ito, T.Tabata, N.Itoh,
c  K.Morita, T.Kato, H.Tawara.
c
c  Modifications :
c
c       i)  Converted into double precision by Manoj.
c
c       ii) Include heavy particle backscattering as described in the
c           paper " Present Status of Data Compilation on Ion
c           Backscattering. - Tatsuo Tabata and Rinsuke Ito
c           Osaka Prefectural Radiation Research Institute,
c           Shinke-cho, Sakai, Osaka 593 " IPPJ-AM-64, Nagoya
c           University Report.
c
c***********************************************************************
c
c  call RNION for the number backscattering coefficient and
c  REION for the energy backscattering coefficient.
c
c  valid up to any mass number for incident ions, provided m1/m2 < 4.5,
c  where m1 is the mass of the incident ion and m2 is the mass of the
c  target atom.
c
      function rnion (energy,nz1,m1,ibnd,theta)
c
c  number backscattering coefficient of light ions.
c
c meaning of input variables
c
c       theta   = angle of incidence in degrees.
c       energy  = incident kinetic energy in eV.
c       nz1     = atomic number of the projectile.
c       m1      = mass number of the projectile.
c       ne      = number of constituent elements in the target.
c       nz2     = array for atomic numbers of the constituents.
c       nw      = array for relative numbers of the constituents.
c
c       example : for Tritium ions of 1 KeV energy incident on the TiO2 
c               (Titanium Dioxide) target, Energy=1000, nz1=1,
c               m1=3, ne=2, nz2(1)=22, nw(1)=1, nz2(2)=8 and 
c               nw(2)=2.
c
c  tables prepared
c       a1t     = table of masses of projectiles
c       a2t     = table of atomic weights of elemental targets
c               IUPAC (1983), C12=12.
c       d       = table of correction factors for the low energy
c               electronic stopping cross-section.
c
      implicit none
      real (kind=R8) :: energy, rnion, reion, theta
      real (kind=R8), save :: a1t(95), a2t(92), d(92)
      integer ibnd, nz1, m1

      integer i, n
      integer, save :: icount
      real (kind=R8) :: th, z1, z2, e, a1, a2, r, r0, as1, as2, sl,
     *  a2eff, sum, delta, z1p, z2p, z2eff, saeff, sbeff, wi, fmu, aa,
     *  zz, eps, sa, sb, w, fi0, tau, beta2, epsb, ch, redengy, psinum,
     *  psiden1, psiden2, psi, sklfact, r0heavy, feps
      real (kind=R8) :: f, fneps, rne
      intrinsic sqrt, tan
c
c  data table a1t : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data a1t /
c  z1 = 1 (for Hydrogen, Deuterium and Tritium)
     . 1.0072766_R8,2.0135536_R8,3.1055011_R8,
c  z1 = 2 (for Helium-3 and Helium-4)
     . 3.0149325_R8,4.0015059_R8, 
c  z1 = 3-10
     . 6.941_R8,9.01218_R8,10.811_R8,
     . 12.011_R8,14.0067_R8,15.9994_R8,18.998403_R8,20.179_R8,
c  z1 = 11 - 20
     . 22.98977_R8,24.305_R8,26.98154_R8,28.0855_R8,30.97376_R8,
     . 32.066_R8,35.453_R8,39.948_R8,39.0893_R8,40.078_R8,
c  z1 = 21 - 30
     . 44.95591_R8,47.88_R8,50.9415_R8,51.9961_R8,54.9380_R8,
     . 55.847_R8,58.9332_R8,58.69_R8,63.546_R8,65.39_R8,
c  z1 = 31 - 40
     . 69.723_R8,72.59_R8,74.9216_R8,78.96_R8,79.904_R8,
     . 83.80_R8,85.4678_R8,87.62_R8,88.9059_R8,91.224_R8,
c  z1 = 41 - 50
     . 92.9064_R8,95.94_R8,99.0_R8,101.07_R8,102.9055_R8,
     . 106.42_R8,107.8682_R8,112.41_R8,114.82_R8,118.710_R8,
c  z1 = 51 - 60
     . 121.75_R8,127.60_R8,126.9045_R8,131.29_R8,132.9054_R8,
     . 137.33_R8,138.9055_R8,140.12_R8,140.9077_R8,144.24_R8,
c  z1 = 61 - 70
     . 145.0_R8,150.36_R8,151.96_R8,157.25_R8,158.9254_R8,
     . 162.50_R8,164.9304_R8,167.26_R8,168.9342_R8,173.04_R8,
c  z1 = 71 - 80
     . 174.967_R8,178.49_R8,180.9479_R8,183.85_R8,186.207_R8,
     . 190.2_R8,192.22_R8,195.08_R8,196.9655_R8,200.59_R8,
c  z1 = 81 - 90
     . 204.383_R8,207.2_R8,208.9804_R8,210.0_R8,210.0_R8,
     . 222.0_R8,223.0_R8,226.0_R8,227.0_R8,232.0_R8,
c  z12 = 91 - 92
     . 231.0_R8,238.0289_R8 /
c
c  data table a2t : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data a2t /
c  z2 = 1 - 10
     . 1.00794_R8,4.002602_R8,6.941_R8,9.01218_R8,10.811_R8,
     . 12.011_R8,14.0067_R8,15.9994_R8,18.998403_R8,20.179_R8,
c  z2 = 11 - 20
     . 22.98977_R8,24.305_R8,26.98154_R8,28.0855_R8,30.97376_R8,
     . 32.066_R8,35.453_R8,39.948_R8,39.0893_R8,40.078_R8,
c  z2 = 21 - 30
     . 44.95591_R8,47.88_R8,50.9415_R8,51.9961_R8,54.9380_R8,
     . 55.847_R8,58.9332_R8,58.69_R8,63.546_R8,65.39_R8,
c  z2 = 31 - 40
     . 69.723_R8,72.59_R8,74.9216_R8,78.96_R8,79.904_R8,
     . 83.80_R8,85.4678_R8,87.62_R8,88.9059_R8,91.224_R8,
c  z2 = 41 - 50
     . 92.9064_R8,95.94_R8,99.0_R8,101.07_R8,102.9055_R8,
     . 106.42_R8,107.8682_R8,112.41_R8,114.82_R8,118.710_R8,
c  z2 = 51 - 60
     . 121.75_R8,127.60_R8,126.9045_R8,131.29_R8,132.9054_R8,
     . 137.33_R8,138.9055_R8,140.12_R8,140.9077_R8,144.24_R8,
c  z2 = 61 - 70
     . 145.0_R8,150.36_R8,151.96_R8,157.25_R8,158.9254_R8,
     . 162.50_R8,164.9304_R8,167.26_R8,168.9342_R8,173.04_R8,
c  z2 = 71 - 80
     . 174.967_R8,178.49_R8,180.9479_R8,183.85_R8,186.207_R8,
     . 190.2_R8,192.22_R8,195.08_R8,196.9655_R8,200.59_R8,
c  z2 = 81 - 90
     . 204.383_R8,207.2_R8,208.9804_R8,210.0_R8,210.0_R8,
     . 222.0_R8,223.0_R8,226.0_R8,227.0_R8,232.0_R8,
c  z2 = 91 - 92
     . 231.0_R8,238.0289_R8 /
c
c  data table d : * * * * * * * * * * * * * * * * * * * * * * * * * * 
c
      data d / 
c  z2 = 1 - 10
     . 0.9341_R8,0.6693_R8,0.6654_R8,0.9712_R8,1.007_R8,
     . 1.024_R8,1.111_R8,0.9699_R8,0.7357_R8,0.6842_R8,
c  z2 = 11 - 20
     . 0.8769_R8,1.290_R8,1.395_R8,1.378_R8,1.063_R8,
     . 1.123_R8,1.632_R8,1.839_R8,1.642_R8,1.749_R8,
c  z2 = 21 - 30
     . 1.638_R8,1.523_R8,1.396_R8,1.236_R8,1.072_R8,
     . 1.083_R8,0.9624_R8,1.085_R8,1.125_R8,1.277_R8,
c  z2 = 31 - 40
     . 1.525_R8,1.675_R8,1.601_R8,1.762_R8,1.679_R8,
     . 1.914_R8,1.696_R8,1.884_R8,1.9_R8,1.993_R8,
c  z2 = 41 - 50
     . 2.039_R8,1.894_R8,2.001_R8,1.795_R8,1.738_R8,
     . 1.534_R8,1.644_R8,1.698_R8,1.816_R8,1.866_R8,
c  z2 = 51 - 60
     . 2.181_R8,2.027_R8,2.240_R8,2.384_R8,2.108_R8,
     . 2.283_R8,2.321_R8,2.159_R8,2.1_R8,2.042_R8,
c  z2 = 61 - 70
     . 1.986_R8,1.932_R8,1.879_R8,1.931_R8,1.779_R8,
     . 1.578_R8,1.492_R8,1.448_R8,1.406_R8,1.365_R8,
c  z2 = 71 - 80
     . 1.394_R8,1.431_R8,1.348_R8,1.3_R8,1.477_R8,
     . 1.439_R8,1.403_R8,1.269_R8,1.376_R8,1.22_R8,
c  z2 = 81 - 90
     . 1.336_R8,1.504_R8,1.683_R8,1.739_R8,1.751_R8,
     . 1.744_R8,1.959_R8,2.115_R8,2.167_R8,2.170_R8,
c  z2 = 91 - 92
     . 2.084_R8,2.050_R8 /

      data icount/0/

      f(z1,z2) = sqrt ( z1 ** (2.0_R8/3.0_R8) + z2 ** (2.0_R8/3.0_R8) )
      fneps(e,z1,a1,z2,a2) = 0.032534_R8 * e / 
     .  (z1 * z2 * ( 1.0_R8 + a1/a2 ) * f(z1,z2))
      rne(th) = r0 + ( 1.0_R8 - r0 ) / ( 1.0_R8 + as1 / tan(th) **
     .  ( 2.0_R8 * as2 ) )
 600  format('  undefined projectile is called  ')
      n = 1
      go to 1

      entry reion(energy,nz1,m1,ibnd,theta)

c  energy backscattering coefficient of light ions

      n = 0
 1    if (nz1.eq.1) go to 3
      if (nz1.eq.2) go to 4
      if (nz1.gt.2) go to 5
 2    write(6,600)
      reion = 0.0_R8
      return
 3    if (m1.lt.1.or.m1.gt.3) go to 2
      a1 = a1t(m1)
      go to 6
 4    if (m1.lt.3.or.m1.gt.4) go to 2
      a1 = a1t( m1+1 )
      go to 6
 5    a1 = a1t( nz1+3 )
c     write(*,*)'atomic # of projectile > 2'
c  finding effective a2 :
cxpb   in mixed materials scenarios, the layer_ arrays may vary over time
      a2eff = 0.0_R8
      sum = 0.0_R8
      do i = 1, layer_nconstituents(ibnd)
        a2 = a2t(layer_nzconstituents(ibnd,i))
        a2eff = a2eff + layer_nrelconstituents(ibnd,i) * a2
        sum = sum + layer_nrelconstituents(ibnd,i)
      enddo
      a2 = a2eff / sum
      delta = a1 / a2
c     write(*,*)'delta = ',delta
      if (delta.gt.4.5_R8 .and. icount.le.100) then 
        icount = icount + 1
        write(*,'(a,f7.2,a)') 
     .   'delta projectile/target > 4.5', delta, '. Continuing...'
      endif
 6    z1 = nz1
      z1p = z1 ** (2.0_R8/3.0_R8)
      z2eff = 0.0_R8
      a2eff = 0.0_R8
      saeff = 0.0_R8
      sbeff = 0.0_R8
      sum = 0.0_R8

      do 70 i = 1, layer_nconstituents(ibnd)
        z2 = layer_nzconstituents(ibnd,i)
        a2 = a2t(layer_nzconstituents(ibnd,i))
        wi = layer_nrelconstituents(ibnd,i)
        z2eff = z2eff + wi * z2
        a2eff = a2eff + wi * a2
        sum = sum + wi
        fmu = a1/a2
        aa = 1.0_R8 + fmu
        z2p = z2 ** (2.0_R8/3.0_R8)
        zz = sqrt ( 1.0_R8 + z1p/z2p )
        eps = fneps(energy,z1,a1,z2,a2)
        sa = 0.0793 * sqrt (eps) / fmu
        w = z1 * zz * aa / z2p
        wi = wi * z1 / w / a2
        saeff = saeff + wi * sa
        sl = d(layer_nzconstituents(ibnd,i)) * z1p * 
     .   ( aa/zz ) ** 1.5_R8 * sa / sqrt(a1)
        if ( z2 - 12.9_R8 ) 10 , 10 , 20
 10     fi0 = 12.0_R8 + 7.0_R8 / z2
        go to 30
 20     fi0 = 9.76_R8 + 58.5_R8 / z2 ** 1.19_R8
 30     tau = eps * z2 ** 2.0_R8 * w / a1 / 3.03055e7_R8
        beta2 = tau * ( tau + 2.0_R8 ) / ( tau + 1.0_R8 ) ** 2
        epsb = 1.022e6_R8 * beta2 / z2 / fi0
        if ( z1 - 2.9_R8 ) 40,40,50
 40     ch = 100.0_R8 * z1 / z2
        go to 60
 50     ch = 5.0_R8
 60     sb = 61.474_R8 / fmu * w * 
     .   ( log ( epsb/(1.0_R8-beta2) + 1.0_R8 + ch/epsb ) - beta2 ) 
     .   / fi0 / epsb
        sb = 2.0_R8 * eps / rangen(eps) + sl * sb / (0.25_R8 * sl 
     .   + sb )
        sbeff = sbeff + wi * sb
 70   continue
      z2 = z2eff / sum
      a2 = a2eff / sum
      sa = saeff
      sb = sbeff
      eps = fneps(energy,z1,a1,z2,a2)
      if (delta.gt.0.5_R8) then
        redengy = 1.0_R8 /
     @   ( 1+(eps/0.104_R8)**0.577_R8 + (eps/0.73_R8)**1.5_R8 )
        fmu = a1/a2
        psinum = 1.0_R8+24.1_R8*fmu**(-3.995_R8)
        psiden1 = fmu / ( (1+fmu) * (eps/1.84_R8+1.0_R8) )
        psiden2 = (1.0_R8+fmu)**3 * eps / ( fmu**3 * (1.0_R8 -
     @   3.0_R8/(2.0_R8*fmu) + 0.9_R8/fmu**2 +
     @   1.0_R8/(2.0_R8*fmu**3) ) * (eps+13.3_R8) )
        psi = psinum / ( psiden1 + psiden2 )

c  sbeff/saeff corresponds to the rho_a/rho-t in the reference paper (of course
c  sbeff is not equal to rho_a nor is saeff equal to rho_t).

        sklfact = z1**(2.0_R8/3.0_R8) * sqrt(a1) * sb/sa *
     @   fmu**2 / (1.0_R8+fmu)**2 / psi
        r0heavy = redengy / sklfact
        r = 1.0_R8 / ( 1.0_R8 + (eps/0.133_R8 ) ** 0.285_R8 ) + 0.530_R8 
     .   / ( 1.0_R8 + (85.0_R8/eps) **1.46_R8 )
        th = theta * pi / 180.0_R8
c       write(*,*) 'did heavy projectile calculation for delta > 0.5 '
        if ( n.eq.1 ) go to 80
        if (th.eq.0.0_R8) then
c         write(*,*) 'did theta=0 calculation of reion'
          reion = r0heavy * r
          go to 79
        else
c         write(*,*) 'did theta > 0 calculation for reion'
          r0 = r0heavy * r
          as1 = 17.9_R8 * eps ** 0.453_R8
          as2 = 0.771_R8 / eps ** 0.014_R8
          reion = rne(th)
          go to 79
        endif
      else
        continue
      endif

c     write(*,*) 'did light projectile calculation delta < 0.5 '

      feps = (2.0_R8 * eps - 3.0_R8 ) / ( eps + 1.0_R8 )
      r0 = 0.705_R8 * sa / sb / ( 1.0_R8 + a1/a2 ) ** feps
      r0 = r0 / (1.0_R8 + (eps/0.047_R8) ** 0.597_R8 + 
     . ( eps / 0.619_R8 ) ** 1.5_R8 )
      th = theta * pi / 180.0_R8
      if (n.eq.1) go to 80
      if (th.eq.0.0_R8) go to 78
      as1 = 17.9_R8 * eps ** 0.453_R8
      as2 = 0.771_R8 / eps ** 0.014_R8
      reion = rne(th)
      go to 79
 78   reion = r0
 79   return
 80   r = 1.0_R8 / ( 1.0_R8 + (eps/0.133_R8 ) ** 0.285 ) + 0.530 / 
     . ( 1.0_R8 + (85.0_R8/eps) ** 1.46_R8 )
      if (delta.gt.0.5_R8) then
c       write(*,*) 'did heavy projectile calculation for rnion'
        if (th.ne.0.0_R8) then
c         write(*,*) 'did theta <> 0 calculation for rnion'
          r0 = r0heavy
          as1 = 7.38_R8 * eps ** 0.359_R8
          as2 = 0.836_R8 /eps ** 0.087_R8
          rnion = rne(th)
          return
        else
c         write(*,*) 'did theta = 0 calculation for rnion'
          rnion = r0heavy
c         write(*,*) r0heavy, rnion, r, reion
          return
        endif
      else
c       write(*,*) 'did light projectile calculation for rnion'
        r0 = r0 / r
        if (th.eq.0.0_R8) go to 88
        as1 = 7.38_R8 * eps ** 0.359_R8
        as2 = 0.836_R8/ eps ** 0.087_R8
        rnion = rne(th)
        return
  88    rnion = r0
        return
      endif

      end function rnion
c
      function rangen(epsiln)
c  reduced range of ions for nuclear stopping only
c  W.D.Wilson , L.G.Haggmark  &  J.P.Biersack : PHYS. REV. B 15 , 2458
c  (1977).
      implicit none
      real (kind=R8) :: epsiln, rangen
      real (kind=R8) :: x, x1, x2
      real (kind=R8), save :: aa, bb, cc
      data aa,bb,cc / 0.56258_R8 , 1.1776_R8 , 0.62680_R8 /
      x = log ( bb * epsiln )
      x1 = ( cc-1.0_R8 ) * x
      x2 = -2.0_R8 * x
      rangen = ( expint(x1) - expint(x2) ) / aa / bb
      return
      end function rangen
c
      function expint(x)
      implicit none
c  incomplete gamma function : gamma ( 0 , x )
c  gamma( 0 , x ) = - ei ( - x )        for x > 0
c  gamma( 0 , x ) = - eibar( - x )      for x < 0
      real (kind=R8) :: expint, x
      real (kind=R8) :: f, p, q, t, z
      real (kind=R8), save :: a(4), b(4), conver
      data a / 0.2677737343_R8 , 8.6347608925_R8 , 18.059016973_R8 ,
     .  8.5733287401_R8 /
      data b / 3.9584969228_R8 , 21.0996530827_R8 , 25.6329561486_R8 ,
     .  9.5733223454_R8 /
      data conver / 1.0e-7_R8 /
      z = dble(x)
      if (x-1.0_R8) 10,40,40
 10   p = 0.57721566490153_R8
      if (x) 12,30,15
 12   if ( x + 83.5_R8 ) 13,15,15
 13   expint = - 1.0e38_R8
      return
 15   p = p + log(abs(z))
      f = 1.0_R8
      t = - z
 20   q = t + p
      if ( abs(p-q) - conver ) 30,30,22
 22   p = q
      t = - z * f * t
      f = f + 1.0_R8
      t = t / f ** 2
      go to 20
 30   expint = - p
      return
 40   if ( x - 84.0_R8 ) 43,43,41
 41   expint = 1.0e-38_R8
      return
 43   p = (((( z + a(4) ) * z + a(3)) * z + a(2)) * z + a(1) ) /
     .  (((( z + b(4) ) * z + b(3)) * z + b(2)) * z + b(1) )
      expint = exp ( -z ) * p / z
      return
      end function expint

      FUNCTION YHAASZ97(E0,TEMP)
c  This is the Haasz/Davis chemical sputtering yield (1997).
c  program provided via Toronto group (P.C. Stangeby/D. Elder)      
c  D. Reiter: single --> double prec.
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  E0  (eV)       -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *                                                                   *
C  *********************************************************************
C            
      use b2mod_types
      IMPLICIT NONE

      REAL (kind=R8), INTENT(IN) :: E0
      REAL (kind=R8), INTENT(IN) :: TEMP
      REAL (kind=R8) :: 
     >           FITC300(4),FITC350(4),FITC400(4),FITC450(4),FITC500(4),
     >           FITC550(4),FITC600(4),FITC650(4),FITC700(4),FITC750(4),
     >           FITC800(4),FITC850(4),FITC900(4),FITC950(4),FITC1000(4)
      REAL(kind=R8) :: POLY_C(4),POLY_H(4),YFIT,FITE0,FRAC
      REAL(kind=R8) :: YHAASZ97,YDAVIS98
      INTEGER I
C
C     Poly. fit c. /       a0,         a1,         a2,         a3
C
      DATA FITC300 / -0.03882_R8, 0.07432_R8,-0.03470_R8, 0.00486_R8/
      DATA FITC350 / -0.05185_R8, 0.10126_R8,-0.05065_R8, 0.00797_R8/
      DATA FITC400 / -0.06089_R8, 0.12186_R8,-0.06240_R8, 0.01017_R8/
      DATA FITC450 / -0.08065_R8, 0.16884_R8,-0.09224_R8, 0.01625_R8/
      DATA FITC500 / -0.08872_R8, 0.19424_R8,-0.10858_R8, 0.01988_R8/
      DATA FITC550 / -0.08728_R8, 0.20002_R8,-0.11420_R8, 0.02230_R8/
      DATA FITC600 / -0.05106_R8, 0.13146_R8,-0.07514_R8, 0.01706_R8/
      DATA FITC650 /  0.07373_R8,-0.13263_R8, 0.09571_R8,-0.01672_R8/
      DATA FITC700 /  0.02722_R8,-0.03599_R8, 0.02064_R8, 0.00282_R8/
      DATA FITC750 /  0.09052_R8,-0.18253_R8, 0.12362_R8,-0.02109_R8/
      DATA FITC800 /  0.02604_R8,-0.05480_R8, 0.04025_R8,-0.00484_R8/
      DATA FITC850 /  0.03478_R8,-0.08537_R8, 0.06883_R8,-0.01404_R8/
      DATA FITC900 /  0.02173_R8,-0.06399_R8, 0.05862_R8,-0.01380_R8/
      DATA FITC950 / -0.00086_R8,-0.01858_R8, 0.02897_R8,-0.00829_R8/
      DATA FITC1000/ -0.01551_R8, 0.01359_R8, 0.00600_R8,-0.00353_R8/
C
      real (kind=R8) :: m1,m2,m3,reducf

      DATA m1/602.39_R8/, m2/202.24_R8/, m3/43.561_R8/
c
C Find right polynomial fit coefficients for a given temperature
C
c
      IF      (TEMP.LE.350.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC300(I)
           POLY_H(I) = FITC350(I)
         ENDDO
         IF (TEMP.LE.300.0_R8) THEN
           FRAC = 1.0_R8
         ELSE
           FRAC = (350.0_R8 - TEMP) / 50.0_R8
         ENDIF
      ELSE IF (TEMP.GT.350.0_R8 .AND. TEMP.LE.400.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC350(I)
           POLY_H(I) = FITC400(I)
         ENDDO
         FRAC = (400.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.400.0_R8 .AND. TEMP.LE.450.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC400(I)
           POLY_H(I) = FITC450(I)
         ENDDO
         FRAC = (450.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.450.0_R8 .AND. TEMP.LE.500.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC450(I)
           POLY_H(I) = FITC500(I)
         ENDDO
         FRAC = (500.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.500.0_R8 .AND. TEMP.LE.550.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC500(I)
           POLY_H(I) = FITC550(I)
         ENDDO
         FRAC = (550.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.550.0_R8 .AND. TEMP.LE.600.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC550(I)
           POLY_H(I) = FITC600(I)
         ENDDO
         FRAC = (600.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.600.0_R8 .AND. TEMP.LE.650.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC600(I)
           POLY_H(I) = FITC650(I)
         ENDDO
         FRAC = (650.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.650.0_R8 .AND. TEMP.LE.700.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC650(I)
           POLY_H(I) = FITC700(I)
         ENDDO
         FRAC = (700.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.700.0_R8 .AND. TEMP.LE.750.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC700(I)
           POLY_H(I) = FITC750(I)
         ENDDO
         FRAC = (750.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.750.0_R8 .AND. TEMP.LE.800.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC750(I)
           POLY_H(I) = FITC800(I)
         ENDDO
         FRAC = (800.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.800.0_R8 .AND. TEMP.LE.850.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC800(I)
           POLY_H(I) = FITC850(I)
         ENDDO
         FRAC = (850.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.850.0_R8 .AND. TEMP.LE.900.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC850(I)
           POLY_H(I) = FITC900(I)
         ENDDO
         FRAC = (900.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.900.0_R8 .AND. TEMP.LE.950.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC900(I)
           POLY_H(I) = FITC950(I)
         ENDDO
         FRAC = (950.0_R8 - TEMP) / 50.0_R8
      ELSE IF (TEMP.GT.950.0_R8) THEN
         DO I = 1,4
           POLY_C(I) = FITC950(I)
           POLY_H(I) = FITC1000(I)
         ENDDO
         IF (TEMP.GT.1000.0_R8) THEN
           FRAC = 0.0_R8
         ELSE
           FRAC = (1000.0_R8 - TEMP) / 50.0_R8
         ENDIF
      ENDIF
C
C Calculate chemical yield according to the 3th poly. fit
C
      IF      (E0.LT.10.0_R8)  THEN
         FITE0 = 10.0_R8
      ELSE IF (E0.GT.200.0_R8) THEN
         FITE0 = 200.0_R8
      ELSE
         FITE0 = E0
      ENDIF 
C
      YFIT = 0.0_R8
      DO I = 1,4
        YFIT = YFIT +         FRAC *POLY_C(I)*LOG10(FITE0)**(I-1)
        YFIT = YFIT + (1.0_R8-FRAC)*POLY_H(I)*LOG10(FITE0)**(I-1)
      ENDDO

      IF (TEMP.LT.300.0_R8) THEN
        YDAVIS98 = 1.0_R8/(m2*(( 300.0_R8/m1)**2 - 1.0_R8)**2 + m3)
        reducf = YFIT/YDAVIS98
        YHAASZ97 = reducf/(m2*((     TEMP/m1)**2 - 1.0_R8)**2 + m3)
      ELSE IF (TEMP.GT.1000.0_R8) THEN
        YDAVIS98 = 1.0_R8/(m2*((1000.0_R8/m1)**2 - 1.0_R8)**2 + m3)
        reducf = YFIT/YDAVIS98
        YHAASZ97 = reducf/(m2*((     TEMP/m1)**2 - 1.0_R8)**2 + m3)
      ELSE
        YHAASZ97 = YFIT
      ENDIF

CW    WRITE(6,*) 'YHAASZ97 = ',YHAASZ97

      RETURN
      END FUNCTION YHAASZ97
c
      FUNCTION YHAASZ97M(E0,TEMP)
c  This is the Haasz/Davis chemical sputtering yield (1998).
c  program provided via Toronto group (P.C. Stangeby/D. Elder)      
c  D. Reiter: single --> double prec.
C
C  *********************************************************************
C  *                                                                   *
C  *  CHEMICAL SPUTTERING FROM Haasz's NEW DATA (February 1997)        *
C  *  - poly. fit: Y = a0 + a1*log(E) + a2*log(E)^2 + a3*log(E)^3      *
C  *  with the addition of a new fit below 10 eV as suggested by       *
C  *  J. Davis and parametrised by G. Porter; now interpolates between *
C  *  5 and 10 eV to lower value (YDAVIS98), and is fixed below 5 eV   *
C  *  E0  (eV)       -  Ion or neutral incident energy                 *
C  *  TEMP (K)       -  Temperature at target or wall                  *
C  *                                                                   *
C  *********************************************************************
C            
      use b2mod_types
      IMPLICIT NONE

      REAL (kind=R8), INTENT(IN) :: E0, TEMP
      real (kind=R8) :: YHAASZ97M, YDAVIS98
      real (kind=R8) :: m1,m2,m3,reducf,FRAC

      DATA m1/602.39_R8/, m2/202.24_R8/, m3/43.561_R8/, reducf/0.2_R8/
C
      IF (E0 .GE. 10.0_R8) THEN
         YHAASZ97M = YHAASZ97(E0,TEMP)
      ELSEIF (E0 .LT. 10.0_R8 .AND. E0 .GE. 5.0_R8) THEN
         FRAC = (E0-5.0_R8)/5.0_R8
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1.0_R8)**2 + m3)
         YHAASZ97M = FRAC*YHAASZ97(E0,TEMP) + (1.0_R8-FRAC)*YDAVIS98
      ELSEIF (E0 .LT. 5.0_R8) THEN
         YDAVIS98 = reducf/(m2*((TEMP/m1)**2 - 1.0_R8)**2 + m3)
         YHAASZ97M = YDAVIS98   
      ENDIF

      RETURN
      END FUNCTION YHAASZ97M

      end module b2mod_sputter
