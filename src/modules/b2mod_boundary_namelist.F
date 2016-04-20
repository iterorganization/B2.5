      module b2mod_boundary_namelist
      use b2mod_types
      use b2mod_indirect
      use b2mod_version
cank-20020724 {
c* Parameters controlling the source scaling procedure are added:
c*  lcbs          : boundary segments representing the core boundary in B2
c*  ncbs          : number of these boundary segments
cank }
      integer NBCD, NNISOD
      parameter (NBCD=20,NNISOD=10)
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   CONPAR(:,:,:), 
     *   MOMPAR(:,:,:), 
     *   ENEPAR(:,:), 
     *   ENIPAR(:,:), 
     *   POTPAR(:,:), 
     *   NIISO(:)
      integer, allocatable, save ::
     *   BCPOS(:), 
     *   BCSTART(:), 
     *   BCEND(:), 
     *   BCCON(:,:), 
     *   BCMOM(:,:), 
     *   BCENE(:), 
     *   BCENI(:), 
     *   BCPOT(:), 
     *   NXISO1(:), 
     *   NXISO2(:), 
     *   NYISO1(:), 
     *   NYISO2(:),
     *   lcbs(:),
     *   BC_LIST_SIZE(:),
     *   BC_LIST_X(:,:),
     *   BC_LIST_Y(:,:)
      character*1, allocatable, save ::
     *   BCCHAR(:)
#else
#     include "DIMENSIONS.F"
      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)
      real (kind=R8), save :: CONPAR(0:nsdmax-1,NBCD,3)
      real (kind=R8), save :: MOMPAR(0:nsdmax-1,NBCD,2)
      real (kind=R8), save :: ENEPAR(NBCD,2)
      real (kind=R8), save :: ENIPAR(NBCD,2)
      real (kind=R8), save :: POTPAR(NBCD,2)
      integer, save :: BCPOS(NBCD)
      integer, save :: BCSTART(NBCD)
      integer, save :: BCEND(NBCD)
      integer, save :: BC_LIST_SIZE(NBCD)
      integer, save :: BC_LIST_X(2*(DEF_NXD+DEF_NYD),NBCD)
      integer, save :: BC_LIST_Y(2*(DEF_NXD+DEF_NYD),NBCD)
      integer, save :: BCCON(0:nsdmax-1,NBCD)
      integer, save :: BCMOM(0:nsdmax-1,NBCD)
      integer, save :: BCENE(NBCD)
      integer, save :: BCENI(NBCD)
      integer, save :: BCPOT(NBCD)
      character*1, save :: BCCHAR(NBCD)
      integer, save :: NXISO1(NNISOD)
      integer, save :: NXISO2(NNISOD)
      integer, save :: NYISO1(NNISOD)
      integer, save :: NYISO2(NNISOD)
      real (kind=R8), save :: NIISO(0:nsdmax-1)
      integer, save :: lcbs(NBCD)
#endif
csw 18feb2014 added for correct implementation of boundary conditions
      integer, save :: l_fbound=0, l_enebound=0, l_timbound=0
      integer, allocatable, save ::
     .  iffboundx(:,:), iffboundy(:,:),
     .  ifenebound(:,:), iftimbound(:,:)

      integer, save :: NBC, NNISO, ncbs
      real (kind=R8) , save :: gammai, gammae, TEISO, TIISO, PHIISO
      logical, save :: lbndusr,lfeedback
      logical, private, save :: write_nml = .true.
      character*256, save :: boundary_filename
      real (kind=R8), save :: boundary_time_mod
      real (kind=R8), save :: boundary_time_switch
      integer, save :: bc_elm_count

      private :: boundary

      namelist /boundary/ NBC,CONPAR,MOMPAR,ENEPAR,ENIPAR,POTPAR,
     1 BCPOS,BCSTART,BCEND,BCCON,BCMOM,BCENE,BCENI,BCPOT,
     1 BC_LIST_SIZE,BC_LIST_X,BC_LIST_Y,
     2 gammai,gammae,lbndusr,lfeedback,
     3 NNISO,NIISO,TEISO,TIISO,PHIISO,NXISO1,NXISO2,NYISO1,NYISO2,
     4 bcchar,
     5 boundary_filename,boundary_time_mod,boundary_time_switch,
     6 lcbs, write_nml
      
      contains

      subroutine initialize_b2mod_boundary(nx,ny,ns,tim)
      implicit none
      integer nx ,ny, ns
      real (kind=R8) :: tim
      integer ix, iy, is, ib
      logical done, catch_up

      call alloc_b2mod_boundary(nx,ny,ns,nbcd,nnisod)
      call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
      call xertst (nbc.le.nbcd, 'increase size in b2stbc for nbcd')
      call xertst (nniso.le.nnisod, 
     1     'increase size in b2stbc for nnisod')
      call write_b2mod_boundary_namelist
      if(boundary_time_mod.gt.0.0_R8) then
         bc_elm_count=tim/boundary_time_mod
      else
         bc_elm_count=0
      endif
      catch_up=boundary_time_mod.gt.0.0_R8
      do while(catch_up)
         if(mod(tim,boundary_time_mod).gt.boundary_time_switch)
     1        then
            write(*,*) 'b2stbc: catching up ...'
            call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
            call xertst (nbc.le.nbcd, 
     1           'increase size in b2stbc for nbcd')
            call xertst (nniso.le.nnisod, 
     1           'increase size in b2stbc for nnisod')
            call write_b2mod_boundary_namelist
            catch_up=boundary_time_mod.gt.0.0_R8
         else
            catch_up=.false.
         endif
      enddo
      if(nniso.gt.0) then
         done=.true.
         do ib = 1, nniso
            do ix = nxiso1(ib), nxiso2(ib)
               do iy = nyiso1(ib), nyiso2(ib)
                  if(region(ix,iy,0).ne.0) then
                     done=.false.
                     region(ix,iy,0) = 0
                     resignore(ix,iy,1)=0
                     resignore(ix,iy,2)=0
                     write(*,*)
     .                    'REGION set to zero for dead cell ',ix,iy
                  endif
               enddo
            enddo
         enddo
         if(.not.done) then
            write(*,*) 'Isolated region was unknown from geometry'
            do ib = 1, nniso
               do ix = nxiso1(ib), nxiso2(ib)
                  bottomiy(ix,nyiso1(ib))=nyiso2(ib)
                  write(*,*) 'ix,iy,iyb',ix,nyiso1(ib),nyiso2(ib)
                  topiy(ix,nyiso2(ib))=nyiso1(ib)
                  write(*,*) 'ix,iy,iyt',ix,nyiso2(ib),nyiso1(ib)
               enddo
               do iy = nyiso1(ib),nyiso2(ib)
                  leftix(nxiso1(ib),iy)=nxiso2(ib)
                  write(*,*) 'ix,iy,ixl',nxiso1(ib),iy,nxiso2(ib)
                  rightix(nxiso2(ib),iy)=nxiso1(ib)
                  write(*,*) 'ix,iy,ixr',nxiso2(ib),iy,nxiso1(ib)
               enddo
            enddo
         endif
      endif
      do ib = 1, nbc
         do is = 0, ns-1
            if (bccon(is,ib).ne.0)
     >           resignore(bc_list_x(1:bc_list_size(ib),ib),
     .           bc_list_y(1:bc_list_size(ib),ib),1) = 0
            if (bcmom(is,ib).ne.0)
     >           resignore(bc_list_x(1:bc_list_size(ib),ib),
     .           bc_list_y(1:bc_list_size(ib),ib),2) = 0
         end do
         if (bcene(ib).ne.0 .or. bceni(ib).ne.0)
     >        resignore(bc_list_x(1:bc_list_size(ib),ib),
     .        bc_list_y(1:bc_list_size(ib),ib),1) = 0
      enddo

csw 21feb2014, setup some additional work variables for boundary conditions
cxpb overwrite the defaults with the namelist data
      iffboundx=0
      iffboundy=0
      if(l_fbound .ne. 0) then
       do ib=1,nbc
         if(bcchar(ib).eq.'S' .or. bcchar(ib).eq.'N') then
           do ix=1, bc_list_size(ib)
             iffboundx(bc_list_x(ix,ib),bc_list_y(ix,ib))=1
           end do
         else if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
           do ix = 1, bc_list_size(ib)
             iffboundy(bc_list_x(ix,ib),bc_list_y(ix,ib))=1
           enddo
         endif
       enddo
      endif

      ifenebound=0
      if(l_enebound .ne. 0) then
       do ib=1,nbc
         if(bcchar(ib).eq.'S' .or. bcchar(ib).eq.'N' .or.
     &      bcchar(ib).eq.'W' .or. bcchar(ib).eq.'E') then
           do ix=1, bc_list_size(ib)
             ifenebound(bc_list_x(ix,ib),bc_list_y(ix,ib))=1
           end do
         end if
       enddo
      endif

      iftimbound=0
      if(l_timbound .ne. 0) iftimbound(:,:)=ifenebound(:,:)

      return
      end subroutine initialize_b2mod_boundary

      subroutine catchup_b2mod_boundary(tim)
      Implicit None
      real (kind=R8) :: tim
      Integer :: nx, ny, ns

      if(boundary_time_mod.gt.0.0_R8) then
         write(*,'(a,1p,3(1x,e14.7),2i4)') 
     .       'ELM: ', tim, mod(tim,boundary_time_mod),
     1        boundary_time_switch,
     2        int(tim/boundary_time_mod),
     3        bc_elm_count
         if(mod(tim,boundary_time_mod).gt.boundary_time_switch.or.
     1        int(tim/boundary_time_mod).gt.bc_elm_count) then
            bc_elm_count=tim/boundary_time_mod
            call read_b2mod_boundary_namelist(nx, ny, ns, .true.)
            call xertst (nbc.le.nbcd, 
     1           'increase size in b2stbc for nbcd')
            call xertst (nniso.le.nnisod, 
     1           'increase size in b2stbc for nnisod')
            call write_b2mod_boundary_namelist
         endif
      endif

      end subroutine catchup_b2mod_boundary

      subroutine alloc_b2mod_boundary(nxd,nyd,nsd,nbcd_in,nnisod_in)
      Implicit None
      integer nxd, nyd, nsd, nbcd_in, nnisod_in, ix, iy
#ifndef ALLOCATE_AND_NAMELIST
      external xertst
#endif

#ifdef ALLOCATE_AND_NAMELIST
      NBCD=nbcd_in
      NNISOD=nnisod_in
      allocate(CONPAR(0:nsd-1,NBCD,3))
      allocate(MOMPAR(0:nsd-1,NBCD,2))
      allocate(ENEPAR(NBCD,2))
      allocate(ENIPAR(NBCD,2))
      allocate(POTPAR(NBCD,2))
      allocate(BCPOS(NBCD))
      allocate(BCSTART(NBCD))
      allocate(BCEND(NBCD))
      allocate(BC_LIST_SIZE(NBCD))
      allocate(BC_LIST_X(2*(NXD+NYD),NBCD))
      allocate(BC_LIST_Y(2*(NXD+NYD),NBCD))
      allocate(BCCON(0:nsd-1,NBCD))
      allocate(BCMOM(0:nsd-1,NBCD))
      allocate(BCENE(NBCD))
      allocate(BCENI(NBCD))
      allocate(BCPOT(NBCD))
      allocate(BCCHAR(NBCD))
      allocate(NXISO1(NNISOD))
      allocate(NXISO2(NNISOD))
      allocate(NYISO1(NNISOD))
      allocate(NYISO2(NNISOD))
      allocate(NIISO(0:nsd-1))
      allocate(lcbs(NBCD))
#else
      call xertst (nxd.le.DEF_NXD, 'faulty parameter nxd')
      call xertst (nyd.le.DEF_NYD, 'faulty parameter nyd')
      call xertst (nsd.le.nsdmax, 'faulty parameter nsd')
      call xertst (nbcd_in.le.nbcd, 'faulty parameter nbcd_in')
      call xertst (nnisod_in.le.nnisod, 'faulty parameter nnisod_in')
#endif
      allocate(iffboundx(-1:nxd,-1:nyd))
      allocate(iffboundy(-1:nxd,-1:nyd))
      allocate(ifenebound(-1:nxd,-1:nyd))
      allocate(iftimbound(-1:nxd,-1:nyd))
      nbc=0
      ncbs=0
      lcbs=0
      bcchar=' '
      bcpos=-2
      bcstart=-2
      bcend=-2
      bc_list_size=0
      bc_list_x=-2
      bc_list_y=-2
      bccon=0
      bcmom=0
      bcene=0
      bceni=0
      bcpot=0
      conpar=0.0_R8
      mompar=0.0_R8
      enepar=0.0_R8
      enipar=0.0_R8
      potpar=0.0_R8
      gammai=5.0_R8/3.0_R8
      gammae=0.5_R8
      lbndusr=.false.
!DPC-2003-10-23 set as a ipget parameter      lfeedback=.false.
      nniso=0
      nxiso1=-2
      nxiso2=-2
      nyiso1=-2
      nyiso2=-2
      niiso=1.0e10_R8
      teiso=1.0_R8
      tiiso=1.0_R8
      phiiso=0.0_R8
      boundary_filename='b2.boundary.parameters'
      boundary_time_mod=0.0_R8
      boundary_time_switch=0.0_R8

csw 21feb2014, setup some additional work variables for boundary conditions
cxpb Have to do this here already in case the namelist is not used
      call ipgeti ('b2news_l_fbound', l_fbound)
      call ipgeti ('b2news_l_enebound', l_enebound)
      call ipgeti ('b2news_l_timbound', l_timbound)
      iffboundx=0
      iffboundy=0
      if(l_fbound .ne. 0) then
       do ix=-1,nxd
        do iy=-1,nyd
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nxd+1) then
            iffboundy(ix,iy)=1
          endif
          if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.nyd+1) then
            iffboundx(ix,iy)=1
          endif
        enddo
       enddo
      endif

      ifenebound=0
      if(l_enebound .ne. 0) then
       do ix=-1,nxd
        do iy=-1,nyd
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nxd+1) then
            ifenebound(ix,iy)=1
          endif
        enddo
       enddo
       do iy=-1,nyd,nyd+1
         ifenebound(:,iy)=1
       enddo
      endif

      iftimbound=0
      if(l_timbound .ne. 0) iftimbound(:,:)=ifenebound(:,:)

      return
      end subroutine alloc_b2mod_boundary

      subroutine dealloc_b2mod_boundary

      if (.not.allocated(iffboundx)) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(CONPAR)
      deallocate(MOMPAR)
      deallocate(ENEPAR)
      deallocate(ENIPAR)
      deallocate(POTPAR)
      deallocate(BCPOS)
      deallocate(BCSTART)
      deallocate(BCEND)
      deallocate(BC_LIST_SIZE)
      deallocate(BC_LIST_X)
      deallocate(BC_LIST_Y)
      deallocate(BCCON)
      deallocate(BCMOM)
      deallocate(BCENE)
      deallocate(BCENI)
      deallocate(BCPOT)
      deallocate(BCCHAR)
      deallocate(NXISO1)
      deallocate(NXISO2)
      deallocate(NYISO1)
      deallocate(NYISO2)
      deallocate(NIISO)
      deallocate(LCBS)
#endif
      deallocate(iffboundx)
      deallocate(iffboundy)
      deallocate(ifenebound)
      deallocate(iftimbound)

      return
      end subroutine dealloc_b2mod_boundary

      subroutine read_b2mod_boundary_namelist(nx, ny, ns, dotest)
      Implicit None
      integer nx, ny, ns
      logical dotest
      integer is, ib, ic, ix, iy, ixx, iyy
      logical covered, corner, guard_cell, found, side
      external xertst, xerrab

      logical :: file_ok, todo, assigned
      logical, allocatable :: old_style(:), done(:)
      character*260 filename
      filename='../' // boundary_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,boundary_version)
        write(*,*) 'Reading namelist boundary from ',trim(filename)
        read(99,boundary)
        write(*,*) 'Read namelist boundary from ',trim(filename)
      else
        call xerrab ('No '//trim(filename))
      endif
      if (nbc.gt.0) then
        allocate(old_style(nbc))
        do ib = 1, nbc
          if (bc_list_size(ib).eq.0) then
            old_style(ib) = .true.
            bc_list_size(ib) = bcend(ib) - bcstart(ib) + 1
            if (bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S' .or. 
     &          bcchar(ib).eq.'X') then
              do ix = bcstart(ib), bcend(ib)
                bc_list_x(ix-bcstart(ib)+1,ib) =  ix
              end do
              bc_list_y(1:bc_list_size(ib),ib) = bcpos(ib)
            else if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W' .or. 
     &               bcchar(ib).eq.'Y') then
              bc_list_x(1:bc_list_size(ib),ib) = bcpos(ib)
              do iy = bcstart(ib), bcend(ib)
                bc_list_y(iy-bcstart(ib)+1,ib) =  iy
              end do
            endif
          else
            old_style(ib) = .false.
          end if
        end do
c.. Make a list of core boundaries
        todo=maxval(lcbs).gt.0
        if (.not.todo) then
          ncbs=0
        else
          ncbs=1
          do while (lcbs(ncbs+1).gt.0 .and. ncbs.lt.nbc)
            ncbs=ncbs+1
          end do
        end if
        allocate(done(nbc))
        done=.false.
        do ib=1,nbc
          if(bcchar(ib).eq.'S') then
            if(old_style(ib)) then
              if(mod(region(bcstart(ib),bcpos(ib),0),4).eq.1 .and.
     .           mod(region(bcend(ib),bcpos(ib),0),4).eq.1) then
                if (.not.todo) then
                  ncbs=ncbs+1
                  lcbs(ncbs)=ib
                  done(ncbs)=.true.
                else
                  ic = 1
                  assigned = .false.
                  do while (ic.le.ncbs)
                    if (lcbs(ic).eq.ib) then
                      call xertst(.not.done(ic),
     .                 'Core boundary already done!')
                      assigned = .true.
                      done(ic) = .true.
                    end if
                    ic = ic + 1
                  end do
                  if (.not.assigned) then
                    ncbs = ncbs + 1
                    done(ncbs) = .true.
                    lcbs(ncbs) = ib
                    write(*,*) 'Added boundary ',ib,
     .               ' to list of core boundaries'
                  end if
                end if
              endif
            else
              found = .true.
              do ix = 1, bc_list_size(ib)
                found = found .and.
     .           mod(region(bc_list_x(ix,ib),bc_list_y(ix,ib),0),4).eq.1
              end do
              if (found) then
                if (.not.todo) then
                  ncbs=ncbs+1
                  lcbs(ncbs)=ib
                  done(ncbs)=.true.
                else
                  ic = 1
                  assigned = .false.
                  do while (ic.le.ncbs)
                    if (lcbs(ic).eq.ib) then
                      call xertst(.not.done(ic),
     .                 'Core boundary already done!')
                      assigned = .true.
                      done(ic) = .true.
                    end if
                    ic = ic + 1
                  end do
                  if (.not.assigned) then
                    ncbs = ncbs + 1
                    done(ncbs) = .true.
                    lcbs(ncbs) = ib
                    write(*,*) 'Added boundary ',ib,
     .               ' to list of core boundaries'
                  end if
                endif
              endif
            endif
          endif
        enddo
        call xertst (1.le.ncbs.or.ny.eq.1, 'faulty parameter ncbs')
        do ib = 1, ncbs
          call xertst (done(ib), 
     .     'Did not find all declared core boundaries!')
          call xertst (1.le.lcbs(ib).and.lcbs(ib).le.nbc,
     .     'faulty parameter lcbs')
        enddo
        deallocate(done)
      endif
      if ((.not.file_ok).or.(.not.dotest)) then
        if (.not.file_ok) write(*,*) 'No ',trim(filename)
        write(*,*) 'Testing skipped!'
        return
      endif

*   ..test input
      call xertst (0.le.nbc, 'faulty input nbc')
      call xertst (0.le.nniso, 'faulty input nniso')
      call xertst (0.0_R8.le.gammae, 'faulty input gammae')
      call xertst (0.0_R8.le.gammai, 'faulty input gammai')
      call xertst (0.0_R8.lt.teiso, 'faulty input teiso')
      call xertst (0.0_R8.lt.tiiso, 'faulty input tiiso')
      call xertst (0.0_R8.le.boundary_time_mod, 
     . 'faulty input boundary_time_mod')
      call xertst (0.0_R8.le.boundary_time_switch, 
     . 'faulty input boundary_time_switch')
      do is = 0, ns-1
        call xertst (0.0_R8.lt.niiso(is), 'faulty input niiso')
      enddo
      do ib = 1, nbc
        if (bcchar(ib).eq.'N') then
          if(old_style(ib)) then
            call xertst (1.le.bcpos(ib).and.bcpos(ib).le.ny+1,
     .       'faulty input bcpos (North)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .       'faulty input bcstart (North)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .       'faulty input bcend (North)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (1.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny+1,
     .                     'faulty input bc_list_y (North)') 
              call xertst (-1.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx+1,
     .                     'faulty input bc_list_x (North)') 
            end do
          endif
        else if (bcchar(ib).eq.'S') then
          if(old_style(ib)) then
            call xertst (-1.le.bcpos(ib).and.bcpos(ib).le.ny-1,
     .       'faulty input bcpos (South)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .       'faulty input bcstart (South)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .       'faulty input bcend (South)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (-1.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny-1,
     .                     'faulty input bc_list_y (South)') 
              call xertst (-1.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx+1,
     .                     'faulty input bc_list_x (South)') 
            end do
          endif
        else if (bcchar(ib).eq.'E') then
          if(old_style(ib)) then
            call xertst (1.le.bcpos(ib).and.bcpos(ib).le.nx+1,
     .       'faulty input bcpos (East)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .       'faulty input bcstart (East)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .       'faulty input bcend (East)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (-1.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny+1,
     .                     'faulty input bc_list_y (East)') 
              call xertst (1.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx+1,
     .                     'faulty input bc_list_x (East)') 
            end do
          endif
        else if (bcchar(ib).eq.'W') then
          if(old_style(ib)) then
            call xertst (-1.le.bcpos(ib).and.bcpos(ib).le.nx-1,
     .       'faulty input bcpos (West)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .       'faulty input bcstart (West)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .       'faulty input bcend (West)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (-1.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny+1,
     .                     'faulty input bc_list_y (West)') 
              call xertst (-1.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx-1,
     .                     'faulty input bc_list_x (West)') 
            end do
          endif
        else if (bcchar(ib).eq.'X') then
          if(old_style(ib)) then
            call xertst (0.le.bcpos(ib).and.bcpos(ib).le.ny,
     .       'faulty input bcpos (X)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .       'faulty input bcstart (X)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .       'faulty input bcend (X)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (0.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny,
     .                     'faulty input bc_list_y (X)') 
              call xertst (-1.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx+1,
     .                     'faulty input bc_list_x (X)') 
            end do
          endif
        else if (bcchar(ib).eq.'Y') then
          if(old_style(ib)) then
            call xertst (0.le.bcpos(ib).and.bcpos(ib).le.nx,
     .       'faulty input bcpos (Y)')
            call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .       'faulty input bcstart (Y)')
            call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .       'faulty input bcend (Y)')
          else
            do ic = 1, bc_list_size(ib)
              call xertst (-1.le.bc_list_y(ic,ib).and.
     .                     bc_list_y(ic,ib).le.ny+1,
     .                     'faulty input bc_list_y (Y)') 
              call xertst (0.le.bc_list_x(ic,ib).and.
     .                     bc_list_x(ic,ib).le.nx,
     .                     'faulty input bc_list_x (Y)') 
            end do
          endif
        endif
        if (old_style(ib)) call xertst (bcstart(ib).le.bcend(ib),
     .   'faulty inputs bcstart > bcend')
        do is = 0, ns-1
          call xertst (0.le.bccon(is,ib), 'faulty input bccon')
          call xertst (0.le.bcmom(is,ib), 'faulty input bcmom')
        enddo
        call xertst (0.le.bcene(ib), 'faulty input bcene')
        call xertst (0.le.bceni(ib), 'faulty input bcene')
        call xertst (0.le.bcpot(ib), 'faulty input bcpot')
      enddo
      do ib = 1, nniso
        call xertst (-1.le.nxiso1(ib).and.nxiso1(ib).le.nx,
     .   'faulty input nxiso1')
        call xertst (-1.le.nxiso2(ib).and.nxiso2(ib).le.nx,
     .   'faulty input nxiso2')
        call xertst (-1.le.nyiso1(ib).and.nyiso1(ib).le.ny,
     .   'faulty input nyiso1')
        call xertst (-1.le.nyiso2(ib).and.nyiso2(ib).le.ny,
     .   'faulty input nyiso2')
        do ix = nxiso1(ib), nxiso2(ib)
          do iy = nyiso1(ib), nyiso2(ib)
            region(ix,iy,0) = 0
            resignore(ix,iy,1)=0
            resignore(ix,iy,2)=0
            do ixx = -1, nx
              do iyy = -1, ny
                if(region(ixx, iyy, 0) .ne. 0) then
                  if(leftix(ixx, iyy) .eq. ix .and.
     &               leftiy(ixx, iyy) .eq. iy) then
                    write(*,*) 'updated left neighbour of ', ixx, iyy
                    leftix(ixx, iyy) = -2
                    leftiy(ixx, iyy) = -2
                  endif
                  if(rightix(ixx, iyy) .eq. ix .and.
     &               rightiy(ixx, iyy) .eq. iy) then
                    write(*,*) 'updated right neighbour of ', ixx, iyy
                    rightix(ixx, iyy) = nx+1
                    rightiy(ixx, iyy) = ny+1
                  endif
                  if(bottomix(ixx, iyy) .eq. ix .and.
     &               bottomiy(ixx, iyy) .eq. iy) then
                    write(*,*) 'updated bottom neighbour of ', ixx, iyy
                    bottomix(ixx, iyy) = -2
                    bottomiy(ixx, iyy) = -2
                  endif
                  if(topix(ixx, iyy) .eq. ix .and.
     &               topiy(ixx, iyy) .eq. iy) then
                    write(*,*) 'updated top neighbour of ', ixx, iyy
                    topix(ixx, iyy) = nx+1
                    topiy(ixx, iyy) = ny+1
                  endif
                endif
              enddo
            enddo
          enddo
        enddo
      enddo
c..  test that all boundary cells to the real domain are covered at least once
c..  we exclude corners from this test
      do ix = -1, nx
        do iy = -1, ny
          guard_cell = region(ix,iy,0).ne.0 .and.
     &     (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &      bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)
          corner = 
     &     (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) .and.
     &     (bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)
          covered = .false.
          do ib = 1, nbc
            if (old_style(ib)) then
              if (bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S') 
     >         covered = covered .or. (iy.eq.bcpos(ib) .and.
     &         ix.ge.bcstart(ib) .and. ix.le.bcend(ib))
              if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') 
     >         covered = covered .or. (ix.eq.bcpos(ib) .and.
     &         iy.ge.bcstart(ib) .and. iy.le.bcend(ib))
            else
              do ic = 1, bc_list_size(ib)
                covered = covered .or. (iy.eq.bc_list_y(ic,ib) .and.
     &           ix.eq.bc_list_x(ic,ib))
              end do
            endif
          enddo
          if (guard_cell .and. .not.covered .and. .not.corner) then
            write(*,'(a,i4,a,i4,a)') 'Guard cell (',ix,',',iy,
     &       ') not included in b2.boundary.parameters!'
            call xerrab('Incomplete boundary definition')
          endif
        enddo
      enddo
c..  test that boundary cells are guard cells
      do ib = 1, nbc
        if (bcchar(ib).eq.'N') then
          if (old_style(ib)) then
            do ix = bcstart(ib), bcend(ib)
              call xertst(topiy(ix,bcpos(ib)).eq.ny+1, 
     &         'inconsistent North boundary!')
            enddo
          else
            do ic = 1, bc_list_size(ib)
              call xertst(
     .         topiy(bc_list_x(ic,ib),bc_list_y(ic,ib)).eq.ny+1, 
     &         'inconsistent North boundary!')
            end do
          endif
        endif
        if (bcchar(ib).eq.'S') then
          if (old_style(ib)) then
            do ix = bcstart(ib), bcend(ib)
              call xertst(bottomiy(ix,bcpos(ib)).eq.-2, 
     &         'inconsistent South boundary!')
            enddo
          else
            do ic = 1, bc_list_size(ib)
              call xertst(
     .         bottomiy(bc_list_x(ic,ib),bc_list_y(ic,ib)).eq.-2, 
     &         'inconsistent South boundary!')
            end do
          endif
        endif
        if (bcchar(ib).eq.'E') then
          if (old_style(ib)) then
            do iy = bcstart(ib), bcend(ib)
              call xertst(rightix(bcpos(ib),iy).eq.nx+1, 
     &         'inconsistent East boundary!')
            enddo
          else
            do ic = 1, bc_list_size(ib)
              call xertst(
     .         rightix(bc_list_x(ic,ib),bc_list_y(ic,ib)).eq.nx+1, 
     &         'inconsistent East boundary!')
            end do
          endif
        endif
        if (bcchar(ib).eq.'W') then
          if (old_style(ib)) then
            do iy = bcstart(ib), bcend(ib)
              call xertst(leftix(bcpos(ib),iy).eq.-2, 
     &         'inconsistent West boundary!')
            enddo
          else
            do ic = 1, bc_list_size(ib)
              call xertst(
     .         leftix(bc_list_x(ic,ib),bc_list_y(ic,ib)).eq.-2, 
     &         'inconsistent West boundary!')
            end do
          endif
        endif
      enddo
c.. Test consistency with region boundaries
      do ic = 1, nncut
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          if (old_style(ib)) then
            found = bcend(ib).eq.leftcut(ic)-1 .and.
     &              bcchar(ib).eq.'S'
          else
            ix = 1
            do while (.not.found .and. ix.le.bc_list_size(ib))
              found = bc_list_x(ix,ib).eq.leftcut(ic)-1 .and.
     .                bc_list_y(ix,ib).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found, 
     .   'Mismatch (1) with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          if (old_style(ib)) then
            found = bcstart(ib).eq.leftcut(ic) .and.
     &              bcchar(ib).eq.'S'
          else
            ix = 1
            do while (.not.found .and. ix.le.bc_list_size(ib))
              found = bc_list_x(ix,ib).eq.leftcut(ic) .and.
     .                bc_list_y(ix,ib).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found, 
     .   'Mismatch (2) with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          if (old_style(ib)) then
            found = bcend(ib).eq.rightcut(ic)-1 .and.
     &              bcchar(ib).eq.'S'
          else
            ix = 1
            do while (.not.found .and. ix.le.bc_list_size(ib))
              found = bc_list_x(ix,ib).eq.rightcut(ic)-1 .and.
     .                bc_list_y(ix,ib).le.topcut(ic)
              ix = ix + 1
            end do
          end if   
        end do
        call xertst (found, 
     .   'Mismatch (3) with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          if (old_style(ib)) then
            found = bcstart(ib).eq.rightcut(ic) .and.
     &              bcchar(ib).eq.'S'
          else
            ix = 1
            do while (.not.found .and. ix.le.bc_list_size(ib))
              found = bc_list_x(ix,ib).eq.rightcut(ic) .and.
     .                bc_list_y(ix,ib).le.topcut(ic)
              ix = ix + 1
            end do
          end if
        end do
        call xertst (found, 
     .   'Mismatch (4) with region boundaries (South) !')
      end do
      return

      end subroutine read_b2mod_boundary_namelist

      subroutine write_b2mod_boundary_namelist
      Implicit None
      If (write_nml) write(*,boundary)
      end subroutine write_b2mod_boundary_namelist

      function get_phi_applied(boundary_namelist,ix,iy)
      use b2mod_types
      use b2mod_plasma
      implicit none
      real (kind=R8) :: get_phi_applied
      integer boundary_namelist, ix, iy
      integer i, ib
      logical found

      get_phi_applied = 0.0_R8
      if (boundary_namelist.eq.1) then
        ib = 1
        found = .false.
        do while (ib.le.nbc .and. .not.found)
          if (bcchar(ib).eq.'S' .or. bcchar(ib).eq.'N' .or.
     .        bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
            i = 1
            do while (i.le.bc_list_size(ib) .and. .not.found)
              found = bc_list_x(i,ib).eq.ix .and.
     &                bc_list_y(i,ib).eq.iy
              i = i + 1
            end do
          endif
          if (.not.found) ib = ib + 1
        enddo
        if (found) then
          select case(bcpot(ib))
          case(1,4)
            get_phi_applied = POTPAR(ib,1)
          case(3)
            get_phi_applied = POTPAR(ib,2)
          end select
        endif
      endif
      return

      end function get_phi_applied

      end module b2mod_boundary_namelist
