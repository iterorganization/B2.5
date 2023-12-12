      module b2mod_boundary_namelist
      use b2mod_types
      use b2mod_indirect
      use b2mod_version
      use b2us_map
      use b2us_geo
      use carre_constants
      use b2mod_ad
     & , only : nsdmax
      implicit none
cank-20020724 {
c* Parameters controlling the source scaling procedure are added:
c*  lcbs          : boundary segments representing the core boundary in B2
c*  ncbs          : number of these boundary segments
cank }
#     include "DIMENSIONS.F"
#ifndef DEF_NBC
#define DEF_NBC 20
#endif
#ifndef DEF_NISO
#define DEF_NISO 10
#endif
      integer NBCD, NNISOD
      parameter (NBCD=DEF_NBC,NNISOD=DEF_NISO)
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   CONPAR(:,:,:),
     *   MOMPAR(:,:,:),
     *   ENEPAR(:,:),
     *   ENIPAR(:,:),
     *   POTPAR(:,:),
     *   ENKPAR(:,:),
     *   ENZPAR(:,:),
     *   NIISO(:),
     *   BC_FACE_ORI(:,:)
      integer, allocatable, save ::
     *   BCPOS(:),
     *   BCSTART(:),
     *   BCEND(:),
     *   ABCSTART(:),
     *   ABCEND(:),
     *   BCCON(:,:),
     *   BCMOM(:,:),
     *   BCENE(:),
     *   BCENI(:),
     *   BCPOT(:),
     *   BCENK(:),
     *   BCENZ(:),
     *   NXISO1(:),
     *   NXISO2(:),
     *   NYISO1(:),
     *   NYISO2(:),
     *   lcbs(:),
     *   BC_LIST_SIZE(:),
     *   BC_LIST_X(:,:),
     *   BC_LIST_Y(:,:),
     *   NBCFACES(:),
     *   BC_FACES(:,:)
      character*1, allocatable, save ::
     *   BCCHAR(:)
     *   ACCHAR(:)
     *   BC_LIST_CHAR(:,:)
#else
      real (kind=R8), save :: CONPAR(0:nsdmax-1,NBCD,3)
      real (kind=R8), save :: MOMPAR(0:nsdmax-1,NBCD,2)
      real (kind=R8), save :: ENEPAR(NBCD,2)
      real (kind=R8), save :: ENIPAR(NBCD,2)
      real (kind=R8), save :: POTPAR(NBCD,2)
      real (kind=R8), save :: ENKPAR(NBCD,2)
      real (kind=R8), save :: ENZPAR(NBCD,2)
      integer, save :: BCPOS(NBCD)
      integer, save :: BCSTART(NBCD)
      integer, save :: BCEND(NBCD)
      integer, save :: ABCSTART(NBCD)
      integer, save :: ABCEND(NBCD)
      integer, save :: BC_LIST_SIZE(NBCD)
      integer, save :: BC_LIST_X(2*(DEF_NXD+DEF_NYD),NBCD)
      integer, save :: BC_LIST_Y(2*(DEF_NXD+DEF_NYD),NBCD)
      character, save :: BC_LIST_CHAR(2*(DEF_NXD+DEF_NYD),NBCD)
      integer, save :: BCCON(0:nsdmax-1,NBCD)
      integer, save :: BCMOM(0:nsdmax-1,NBCD)
      integer, save :: BCENE(NBCD)
      integer, save :: BCENI(NBCD)
      integer, save :: BCPOT(NBCD)
      integer, save :: BCENK(NBCD)
      integer, save :: BCENZ(NBCD)
      character*1, save :: BCCHAR(NBCD)
      character*1, save :: ACCHAR(NBCD)
      integer, save :: NXISO1(NNISOD)
      integer, save :: NXISO2(NNISOD)
      integer, save :: NYISO1(NNISOD)
      integer, save :: NYISO2(NNISOD)
      real (kind=R8), save :: NIISO(0:nsdmax-1)
      integer, save :: lcbs(NBCD)
      integer, save :: NBCFACES(NBCD)
      integer, save :: BC_FACES(2*(DEF_NXD+DEF_NYD),NBCD)
      real (kind=R8), save :: BC_FACE_ORI(2*(DEF_NXD+DEF_NYD),NBCD)
#endif
#ifdef WG_OUT
csw 18feb2014 added for correct implementation of boundary conditions
      integer, save :: l_fbound=0, l_enebound=0, l_timbound=0
      integer, allocatable, save ::
     .  iffboundx(:,:), iffboundy(:,:),
     .  ifenebound(:,:), iftimbound(:,:)
#endif

      integer, save :: NBC, NNISO, ncbs
      real (kind=R8) , save :: gammai, gammae, TEISO, TIISO, PHIISO
      logical, save :: lbndusr, lfeedback
      logical, save :: write_nml_bnd = .true.
      character*256, save :: boundary_filename
      real (kind=R8), save :: boundary_time_mod
      real (kind=R8), save :: boundary_time_switch
      integer, save :: bc_elm_count
      logical, allocatable, save :: isOnFluxBoundary(:,:)

      interface read_b2mod_boundary_namelist
        module procedure read_b2mod_boundary_namelist_st
        module procedure read_b2mod_boundary_namelist_us
      end interface

      private :: boundary

      namelist /boundary/ NBC,CONPAR,MOMPAR,ENEPAR,ENIPAR,POTPAR,ENKPAR,
     1 ENZPAR,
     1 BCPOS,BCSTART,BCEND,BCCON,BCMOM,BCENE,BCENI,BCPOT,BCENK,BCENZ,
     1 BC_LIST_SIZE,BC_LIST_X,BC_LIST_Y,BC_LIST_CHAR,
     2 gammai,gammae,lbndusr,lfeedback,
     3 NNISO,NIISO,TEISO,TIISO,PHIISO,NXISO1,NXISO2,NYISO1,NYISO2,
     4 bcchar,
     5 boundary_filename,boundary_time_mod,boundary_time_switch,
     6 lcbs, write_nml_bnd

      contains

      subroutine initialize_b2mod_boundary(ns,mpg,tim)
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent(inout) :: mpg
      real (kind=R8), intent(in) :: tim
#ifdef WG_TODO
      integer ix, iy, is, ib
      logical done
#endif
      logical catch_up
      external xertst

!WG_TODO need to clean up here;
!WG_TODO for now: call alloc_b2mod_boundary with dummy nx=1,ny=1;
      call alloc_b2mod_boundary(1,1,ns,nbcd,nnisod)
      call read_b2mod_boundary_namelist(ns, mpg, .true.)
      call xertst (nbc.le.nbcd, 'increase DEF_NBC in DIMENSIONS.F')
      call xertst (nniso.le.nnisod,
     1     'increase DEF_NISO in DIMENSIONS.F')
      call write_b2mod_boundary_namelist
      if(boundary_time_mod.gt.0.0_R8) then
         bc_elm_count=int(tim/boundary_time_mod)
      else
         bc_elm_count=0
      endif
      if(boundary_time_mod.gt.0.0_R8) then
        catch_up=mod(tim,boundary_time_mod).ge.boundary_time_switch.and.
     1           boundary_time_switch.gt.0.0_R8
      else if (boundary_time_switch.gt.0.0_R8) then
        catch_up=tim.ge.boundary_time_switch
      else
        catch_up=.false.
      end if
      do while(catch_up)
        write(*,*) 'b2stbc: catching up ...'
        call read_b2mod_boundary_namelist(ns, mpg, .true.)
        call xertst (nbc.le.nbcd,
     1           'increase DEF_NBC in DIMENSIONS.F')
        call xertst (nniso.le.nnisod,
     1           'increase DEF_NISO in DIMENSIONS.F')
        call write_b2mod_boundary_namelist
        if(boundary_time_mod.gt.0.0_R8) then
          catch_up=boundary_time_switch.gt.0.0_R8.and.
     1     mod(tim,boundary_time_mod).ge.boundary_time_switch
        else if (boundary_time_switch.gt.0.0_R8) then
          catch_up=tim.ge.boundary_time_switch
        else
          catch_up=.false.
        endif
      enddo
#ifdef WG_TODO
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
#endif
#ifdef WG_TODO
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
#endif

      return
      end subroutine initialize_b2mod_boundary

      subroutine catchup_b2mod_boundary(ns, mpg, tim)
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent(inout) :: mpg
      real (kind=R8), intent(in) :: tim
      logical catch_up
      external xertst

      if(boundary_time_mod.gt.0.0_R8) then
        write(*,'(a,1p,3(1x,e14.7),2i4)')
     .       'ELM: ', tim, mod(tim,boundary_time_mod),
     1        boundary_time_switch,
     2        int(tim/boundary_time_mod),
     3        bc_elm_count
        catch_up=
     1   (mod(tim,boundary_time_mod).ge.boundary_time_switch.and.
     1    boundary_time_switch.gt.0.0_R8).or.
     2    int(tim/boundary_time_mod).gt.bc_elm_count
      else if (boundary_time_switch.gt.0.0_R8) then
        catch_up=tim.ge.boundary_time_switch
        if(catch_up) write(*,'(a,a)')
     .   'Switching to next boundary namelist file ',
     .    trim(boundary_filename)
      else
        catch_up=.false.
      end if
      if(catch_up) then
        if(boundary_time_mod.gt.0.0_R8)
     >   bc_elm_count=int(tim/boundary_time_mod)
        call read_b2mod_boundary_namelist(ns, mpg, .true.)
        call xertst (nbc.le.nbcd, 'increase DEF_NBC in DIMENSIONS.F')
        call xertst (nniso.le.nnisod,
     .                           'increase DEF_NISO in DIMENSIONS.F')
        call write_b2mod_boundary_namelist
      endif
      return

      end subroutine catchup_b2mod_boundary

      subroutine alloc_b2mod_boundary(nxd,nyd,nsd,nbcd_in,nnisod_in)
      implicit none
      integer, intent(in) :: nxd, nyd, nsd, nbcd_in, nnisod_in
      external ipgeti, xertst

#ifdef ALLOCATE_AND_NAMELIST
      allocate(CONPAR(0:nsd-1,nbcd_in,3))
      allocate(MOMPAR(0:nsd-1,nbcd_in,2))
      allocate(ENEPAR(nbcd_in,2))
      allocate(ENIPAR(nbcd_in,2))
      allocate(POTPAR(nbcd_in,2))
      allocate(ENKPAR(nbcd_in,2))
      allocate(ENZPAR(nbcd_in,2))
      allocate(BCPOS(nbcd_in))
      allocate(BCSTART(nbcd_in))
      allocate(BCEND(nbcd_in))
      allocate(ABCSTART(nbcd_in))
      allocate(ABCEND(nbcd_in))
      allocate(BC_LIST_SIZE(nbcd_in))
      allocate(BC_LIST_X(2*(NXD+NYD),nbcd_in))
      allocate(BC_LIST_Y(2*(NXD+NYD),nbcd_in))
      allocate(BCCON(0:nsd-1,nbcd_in))
      allocate(BCMOM(0:nsd-1,nbcd_in))
      allocate(BCENE(nbcd_in))
      allocate(BCENI(nbcd_in))
      allocate(BCPOT(nbcd_in))
      allocate(BCENK(nbcd_in))
      allocate(BCENZ(nbcd_in))
      allocate(BCCHAR(nbcd_in))
      allocate(ACCHAR(nbcd_in))
      allocate(NXISO1(nnisod_in))
      allocate(NXISO2(nnisod_in))
      allocate(NYISO1(nnisod_in))
      allocate(NYISO2(nnisod_in))
      allocate(NIISO(0:nsd-1))
      allocate(lcbs(nbcd_in))
#else
#ifdef WG_OUT
      call xertst (nxd.le.DEF_NXD, 'faulty parameter nxd')
      call xertst (nyd.le.DEF_NYD, 'faulty parameter nyd')
#endif
      call xertst (nsd.le.nsdmax, 'faulty parameter nsd')
#endif
      call xertst (nbcd_in.le.nbcd, 'faulty parameter nbcd_in')
      call xertst (nnisod_in.le.nnisod, 'faulty parameter nnisod_in')
      allocate(isOnFluxBoundary(-1:nxd,-1:nyd))

      nbc=0
      ncbs=0
      lcbs=0
      bcchar=' '
      acchar=' '
      bcpos=-2
      bcstart=-2
      bcend=-2
      abcstart=-2
      abcend=-2
      bc_list_size=0
      bc_list_x=-2
      bc_list_y=-2
      bc_list_char=' '
      bccon=0
      bcmom=0
      bcene=0
      bceni=0
      bcpot=0
      bcenk=0
      bcenz=0
      conpar=0.0_R8
      mompar=0.0_R8
      enepar=0.0_R8
      enipar=0.0_R8
      potpar=0.0_R8
      enkpar=0.0_R8
      enzpar=0.0_R8
      gammai=1.0_R8
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

      return
      end subroutine alloc_b2mod_boundary

      subroutine dealloc_b2mod_boundary
      implicit none

      if (.not.allocated(isOnFluxBoundary)) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(CONPAR)
      deallocate(MOMPAR)
      deallocate(ENEPAR)
      deallocate(ENIPAR)
      deallocate(POTPAR)
      deallocate(ENKPAR)
      deallocate(ENZPAR)
      deallocate(BCPOS)
      deallocate(BCSTART)
      deallocate(BCEND)
      deallocate(ACSTART)
      deallocate(ACEND)
      deallocate(BC_LIST_SIZE)
      deallocate(BC_LIST_X)
      deallocate(BC_LIST_Y)
      deallocate(BCCON)
      deallocate(BCMOM)
      deallocate(BCENE)
      deallocate(BCENI)
      deallocate(BCPOT)
      deallocate(BCENK)
      deallocate(BCENZ)
      deallocate(BCCHAR)
      deallocate(ACCHAR)
      deallocate(NXISO1)
      deallocate(NXISO2)
      deallocate(NYISO1)
      deallocate(NYISO2)
      deallocate(NIISO)
      deallocate(LCBS)
#endif
      deallocate(isOnFluxBoundary)

      return
      end subroutine dealloc_b2mod_boundary

      subroutine read_b2mod_boundary_namelist_st(nx, ny, ns, dotest)
      implicit none
      integer, intent(in) :: nx, ny, ns
      logical, intent(in) :: dotest
      integer is, ib, ic, ix, iy, ixx, iyy, iyn
      logical covered, corner, guard_cell, found
      logical is_flux_surface

      logical :: file_ok, todo, assigned
      logical, allocatable :: old_style(:), done(:)
      character*260 filename
      external find_file, xerrab, xertst

      filename=boundary_filename
      boundary_time_mod = 0.0_R8
      boundary_time_switch = 0.0_R8
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,boundary_version)
        write(*,*) 'Reading namelist boundary from ',trim(filename)
#ifdef USE_MPI
        read(99,boundary,err=999)
#else
        read(99,boundary)
#endif
        write(*,*) 'Read namelist boundary from ',trim(filename)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif

      acchar = bcchar
      abcstart = bcstart
      abcend = bcend

      if (nbc.gt.0) then
        allocate(old_style(nbc))
        do ib = 1, nbc
          if (bc_list_size(ib).eq.0) then
            if (BCCHAR(ib)=="A") then
               ! Automatic creation of boundary condition list from structures
               old_style(ib) = .false.
               call build_boundary_list( bcstart(ib), bcend(ib),
     $              bc_list_size(ib),
     $              bc_list_x(:,ib), bc_list_y(:,ib),
     $              bc_list_char(:,ib))
            else
              ! Automatic creation of boundary condition list from old style cell block definition
              old_style(ib) = .true.
              bc_list_size(ib) = bcend(ib) - bcstart(ib) + 1
              bc_list_char(1:bc_list_size(ib),ib) = bcchar(ib)
              if (bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S' .or.
     &            bcchar(ib).eq.'X') then
                do ix = bcstart(ib), bcend(ib)
                  bc_list_x(ix-bcstart(ib)+1,ib) =  ix
                end do
                bc_list_y(1:bc_list_size(ib),ib) = bcpos(ib)
              else if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W' .or.
     &                 bcchar(ib).eq.'Y') then
                bc_list_x(1:bc_list_size(ib),ib) = bcpos(ib)
                do iy = bcstart(ib), bcend(ib)
                  bc_list_y(iy-bcstart(ib)+1,ib) =  iy
                end do
              endif
            endif
          else
            old_style(ib) = .false.
          end if
        end do
      endif
      isOnFluxBoundary = .false.
      do ib = 1, nbc
        IS_FLUX_SURFACE = BCCHAR(IB).EQ.'N' .or. BCCHAR(IB).EQ.'S'
        IF (BCCHAR(IB).eq.'A') THEN
          IX = 1
          IS_FLUX_SURFACE = BC_LIST_CHAR(IX,IB).EQ.'N' .or.
     .                      BC_LIST_CHAR(IX,IB).EQ.'S'
          DO WHILE (IS_FLUX_SURFACE .AND. IX.LT.BC_LIST_SIZE(IB))
            IF (BC_LIST_CHAR(IX,IB).EQ.'N') THEN
              IY = BOTTOMIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ELSEIF (BC_LIST_CHAR(IX,IB).EQ.'S') THEN
              IY = TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB))
            ENDIF
            IF (BC_LIST_CHAR(IX+1,IB).EQ.'N') THEN
              IYN = BOTTOMIY(BC_LIST_X(IX+1,IB),BC_LIST_Y(IX+1,IB))
            ELSEIF (BC_LIST_CHAR(IX+1,IB).EQ.'S') THEN
              IYN = TOPIY(BC_LIST_X(IX+1,IB),BC_LIST_Y(IX+1,IB))
            ENDIF
            IS_FLUX_SURFACE = IS_FLUX_SURFACE .AND.
     .       IY.EQ.IYN .AND.
     .       BC_LIST_CHAR(IX,IB).EQ.BC_LIST_CHAR(IX+1,IB)
            IX = IX+1
          END DO
        END IF
        IF (IS_FLUX_SURFACE) THEN
          DO IX = 1, BC_LIST_SIZE(IB)
            isOnFluxBoundary(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)) = .true.
          END DO
        END IF
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
          if(bcchar(ib).eq.'S' .or. bcchar(ib).eq.'A') then
            if(old_style(ib)) then
              if(on_closed_surface(bcstart(ib),bcpos(ib)) .and.
     .           on_closed_surface(bcend(ib),bcpos(ib))) then
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
     .           on_closed_surface(bc_list_x(ix,ib),bc_list_y(ix,ib))
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
        if (ny.gt.1) then
          if (nnreg(0).gt.1 .or. periodic_bc.eq.1) then
            call xertst(1.le.ncbs, 'faulty parameters ncbs')
          end if
        end if
        do ib = 1, ncbs
          call xertst (1.le.lcbs(ib).and.lcbs(ib).le.nbc,
     .     'faulty parameter lcbs')
          call xertst (done(ib),
     .     'Did not find all declared core boundaries!')
        enddo
        deallocate(done)
!WG      endif
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
        call xertst (0.le.bcenk(ib), 'faulty input bcenk')
        call xertst (0.le.bcenz(ib), 'faulty input bcenz')
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
          enddo
        enddo
      enddo
      do ib = 1, nniso
        do ix = nxiso1(ib), nxiso2(ib)
          do iy = nyiso1(ib), nyiso2(ib)
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
#ifdef USE_MPI
  999 call xerrab ('Error while reading BOUNDARY namelist!')
#endif
      end subroutine read_b2mod_boundary_namelist_st

      subroutine read_b2mod_boundary_namelist_us(ns, m, dotest)
      implicit none
      integer, intent(in) :: ns
      type(mapping), intent(inout) :: m
      logical, intent(in) :: dotest
      integer is, ib, ic, noss, iss, indss, idb, nbcfc, inbc
      integer :: ian, ien, nfaces
      character*3 ss
#ifdef WG_OUT
      logical covered, corner, guard_cell, found
#else
#ifdef WG_TODO
      integer ix, iy, ixx, iyy
      logical found
      logical, allocatable :: old_style(:)
#endif
#endif

      logical file_ok, todo, assigned
      logical, allocatable :: done(:)
      character*260 filename
#ifndef F2003
      integer newunit
      external newunit
#endif
      external find_file, xerrab, check_boundary_labels

      filename=boundary_filename
      boundary_time_mod = 0.0_R8
      boundary_time_switch = 0.0_R8
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        call cfverr(99,boundary_version)
        write(*,*) 'Reading namelist boundary from ',trim(filename)
#ifdef USE_MPI
        read(99,boundary,err=999)
#else
        read(99,boundary)
#endif
        write(*,*) 'Read namelist boundary from ',trim(filename)
        close(99)
      else
        call xerrab ('No '//trim(filename))
      endif

      acchar = bcchar
      abcstart = bcstart
      abcend = bcend

      if (nbc > 0) then

        m%nBc = nbc
        m%mxnBc = count(m%fcLbl /= 0)
        call alloc_mapping_bc (m)
        m%bcFcP = 0

#ifdef F2003
        open (newunit=idb,file='debug_bound.out')
#else
        idb=newunit()
        open (unit=idb,file='debug_bound.out')
#endif
        inbc = 0
        do ib = 1, nbc
          if (bcchar(ib) == 'A') then
            nbcfc = 0
            m%bcFcP(ib,1) = inbc +1
            noss = bcend(ib) - bcstart(ib) + 1
            do iss = 1, noss
              indss = bcstart(ib)+iss-1
              write (*,*) ' surface index ',indss
! find all faces belonging to surface structure INDSS
              nfaces = 0
              call find_faces (indss, inbc, nbcfc,
     >             m%mxnbc,m%bcFc, m%bcFcOr, m, idb, nfaces)
              write(ss,"(I0)") indss
              call xertst(nfaces.gt.0, 
     &         'No faces found for boundary face label = '//ss)
            end do              ! iss
            m%bcFcP(ib,2) = nbcfc
          else
            call xerrab ('BCCHAR /= A ')
          endif
        end do
        close(idb)
      end if
!      call bc_to_struct (m)

      ! now check that all boundary faces are assigned to a BC
      call check_boundary_labels(m, nBc, bcstart(1:nbc), bcend(1:nbc),
     & .true.)

      if (any(bc_list_size > 0))
     .  call xerrab
     .  ('read_b2mod_boundary_namelist_us -- bc_list_size > 0')

      call bc_to_struct(m)

      if (nbc.gt.0) then
#ifdef WG_TODO
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
#endif
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
          if(bcchar(ib).eq.'S' .or. bcchar(ib).eq.'A') then
            ian = m%bcCvP(ib,1)
            ien = m%bcCvP(ib,1) + m%bcCvP(ib,2) - 1
            if (m%cvonClosedSurface(m%bcCv(ian,2)) .and.
     .          m%cvonClosedSurface(m%bcCv(ien,2))) then
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
     .                   'Core boundary already done!')
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
     .                 ' to list of core boundaries'
                end if
              end if
            end if
#ifdef WG_TODO
            if(old_style(ib)) then
            else
              found = .true.
              do ix = 1, bc_list_size(ib)
                found = found .and.
     .           on_closed_surface(bc_list_x(ix,ib),bc_list_y(ix,ib))
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
#endif
          endif
        enddo
        if (maxval(m%fsVxP(:,2)).gt.2) then
          if (m%nnreg(0).gt.1 .or. m%periodic_bc.eq.1) then
            call xertst(1.le.ncbs, 'faulty parameters ncbs')
          end if
        end if
        do ib = 1, ncbs
          call xertst (1.le.lcbs(ib).and.lcbs(ib).le.nbc,
     .     'faulty parameter lcbs')
          call xertst (done(ib),
     .     'Did not find all declared core boundaries!')
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
#ifdef WG_TODO
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
#endif
        call xertst (bcstart(ib).le.bcend(ib),
     .   'faulty inputs bcstart > bcend')
        do is = 0, ns-1
          call xertst (0.le.bccon(is,ib), 'faulty input bccon')
          call xertst (0.le.bcmom(is,ib), 'faulty input bcmom')
        enddo
        call xertst (0.le.bcene(ib), 'faulty input bcene')
        call xertst (0.le.bceni(ib), 'faulty input bcene')
        call xertst (0.le.bcpot(ib), 'faulty input bcpot')
        call xertst (0.le.bcenk(ib), 'faulty input bcenk')
        call xertst (0.le.bcenz(ib), 'faulty input bcenz')
      enddo
#ifdef WG_OUT
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
          enddo
        enddo
      enddo
      do ib = 1, nniso
        do ix = nxiso1(ib), nxiso2(ib)
          do iy = nyiso1(ib), nyiso2(ib)
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
#endif
      return

#ifdef USE_MPI
  999 call xerrab ('Error while reading BOUNDARY namelist!')
#endif
      end subroutine read_b2mod_boundary_namelist_us

      subroutine write_b2mod_boundary_namelist
      implicit none
      if (write_nml_bnd) write(*,boundary)
      end subroutine write_b2mod_boundary_namelist

      function get_phi_applied(boundary_namelist,ix,iy)
      use b2mod_types
      use b2mod_plasma
      implicit none
      real (kind=R8) :: get_phi_applied
      integer, intent(in) :: boundary_namelist, ix, iy
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

      function getPhiApplied(boundary_namelist,mpg,iCv)
      use b2mod_types
      use b2us_map
      implicit none
      real (kind=R8) :: getPhiApplied
      integer boundary_namelist, iCv
      type (mapping), intent (in) :: mpg
      integer i, ib
      logical found

      getPhiApplied = 0.0_R8
      if (boundary_namelist.eq.1) then
        ib = 1
        found = .false.
        do while (ib.le.nbc .and. .not.found)
          if (bcchar(ib).eq.'S' .or. bcchar(ib).eq.'N' .or.
     .        bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') then
            i = 1
            do while (i.le.mpg%bcCvP(ib,2) .and. .not.found)
              found = mpg%bcCv(mpg%bcCvP(ib,1)+i-1,1) .eq. iCv
              i = i + 1
            end do
          endif
          if (.not.found) ib = ib + 1
        enddo
        if (found) then
          select case(bcpot(ib))
          case(1,4)
            getPhiApplied = POTPAR(ib,1)
          case(3)
            getPhiApplied = POTPAR(ib,2)
          end select
        endif
      endif
      return

      end function getPhiApplied

***********************************************************************

      subroutine conv_bcparms (ns, m)
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent (inout) :: m
      integer :: ib, ibc, is, inbc, noss, iss, indss,
     .           i, iFc, iCv, idb, nfaces
      integer :: bcs(nbc), bce(nbc), bcl(nbc)
      integer :: minlbl(nbc), maxlbl(nbc)
      integer, allocatable :: ifchit(:)
      character(1) :: lchar(nbc)
#ifndef F2003
      integer newunit
      external newunit
#endif
      intrinsic huge

      m%nBc   = nbc     ! number of boundary conditions
      m%mxnBc = 0       ! number of cells with imposed BC
      do ib = 1, nbc
        bcs(ib) = bcstart(ib)
        bce(ib) = bcend(ib)
        if (bcchar(ib).eq.'A') then
! find the number of faces for this boundary condition
cwdk need more complex test here; present scheme double counts;
cwdk also need to correct for degenerate faces?
          bcl(ib) = count((cflags(:,:,2:5) >= bcstart(ib)) .and.
     >                    (cflags(:,:,2:5) <= bcend(ib)))
        else
        ! exclude corner cells if present
        if (bcchar(ib).eq.'N') then
          if (m%imapCv(bottomix(bcstart(ib),bcpos(ib)),
     >                 bottomiy(bcstart(ib),bcpos(ib))) .gt. m%nCi)
     >      bcs(ib) = bcs(ib) + 1
          if (m%imapCv(bottomix(bcend(ib),bcpos(ib)),
     >                 bottomiy(bcend(ib),bcpos(ib))) .gt. m%nCi)
     >      bce(ib) = bce(ib) - 1
        elseif (bcchar(ib).eq.'S') then
          if (m%imapCv(topix(bcstart(ib),bcpos(ib)),
     >                 topiy(bcstart(ib),bcpos(ib))) .gt. m%nCi)
     >      bcs(ib) = bcs(ib) + 1
          if (m%imapCv(topix(bcend(ib),bcpos(ib)),
     >                 topiy(bcend(ib),bcpos(ib))) .gt. m%nCi)
     >      bce(ib) = bce(ib) - 1
        elseif (bcchar(ib).eq.'E') then
          if (m%imapCv(leftix(bcpos(ib),bcstart(ib)),
     >                 leftiy(bcpos(ib),bcstart(ib))) .gt. m%nCi)
     >      bcs(ib) = bcs(ib) + 1
          if (m%imapCv(leftix(bcpos(ib),bcend(ib)),
     >                 leftiy(bcpos(ib),bcend(ib))) .gt. m%nCi)
     >      bce(ib) = bce(ib) - 1
        elseif (bcchar(ib).eq.'W') then
          if (m%imapCv(rightix(bcpos(ib),bcstart(ib)),
     >                 rightiy(bcpos(ib),bcstart(ib))) .gt. m%nCi)
     >      bcs(ib) = bcs(ib) + 1
          if (m%imapCv(rightix(bcpos(ib),bcend(ib)),
     >                 rightiy(bcpos(ib),bcend(ib))) .gt. m%nCi)
     >      bce(ib) = bce(ib) - 1
        else
          call xerrab ('bcchar not yet available for WG')
        endif
        bcl(ib) = bce(ib) - bcs(ib) + 1
        endif
        m%mxnBc = m%mxnBc + bcl(ib)
      enddo
      m%mxnBc = m%nCg
      if (m%mxnBc .ne. m%nCg )
     >  call xerrab ('conv_bcparms -- wrong mxnBc')

      call alloc_mapping_bc (m)
#ifdef F2003
      open (newunit=idb,file='bound.out')
#else
      idb = newunit()
      open (unit=idb,file='bound.out')
#endif

      minlbl = huge(1)
      maxlbl = -huge(1)
      lchar = ' '

      inbc = 0 ! index in bcCv
      do ib = 1, nbc
        m%bcFcP(ib,1) = inbc + 1
        m%bcFcP(ib,2) = bcl(ib)
        if (bcchar(ib).eq.'S') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
    ! number of the corresponding domain face
            m%bcFc(inbc) =
     >         m%cvFc(m%cvFcP(m%imapCv(bcs(ib)+ibc-1,bcpos(ib)),1))
            m%bcFcOr(inbc) = -1.0_R8
            minlbl(ib) = min(minlbl(ib),m%fcLbl(m%bcFc(inbc)))
            maxlbl(ib) = max(maxlbl(ib),m%fcLbl(m%bcFc(inbc)))
            lchar(ib) = 'A'
          enddo
        elseif (bcchar(ib).eq.'N') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) =
     >         m%cvFc(m%cvFcP(m%imapCv(bcs(ib)+ibc-1,bcpos(ib)),1))
            m%bcFcOr(inbc) = 1.0_R8
            minlbl(ib) = min(minlbl(ib),m%fcLbl(m%bcFc(inbc)))
            maxlbl(ib) = max(maxlbl(ib),m%fcLbl(m%bcFc(inbc)))
            lchar(ib) = 'A'
          enddo
        elseif (bcchar(ib).eq.'W') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) =
     >         m%cvFc(m%cvFcP(m%imapCv(bcpos(ib),bcs(ib)+ibc-1),1))
            m%bcFcOr(inbc) = -1.0_R8
            minlbl(ib) = min(minlbl(ib),m%fcLbl(m%bcFc(inbc)))
            maxlbl(ib) = max(maxlbl(ib),m%fcLbl(m%bcFc(inbc)))
            lchar(ib) = 'A'
          enddo
        elseif (bcchar(ib).eq.'E') then
          do ibc = 1, bcl(ib)
            inbc = inbc + 1
            m%bcFc(inbc) =
     >         m%cvFc(m%cvFcP(m%imapCv(bcpos(ib),bcs(ib)+ibc-1),1))
            m%bcFcOr(inbc) = 1.0_R8
            minlbl(ib) = min(minlbl(ib),m%fcLbl(m%bcFc(inbc)))
            maxlbl(ib) = max(maxlbl(ib),m%fcLbl(m%bcFc(inbc)))
            lchar(ib) = 'A'
          enddo
        elseif (bcchar(ib).eq.'A') then
          minlbl(ib) = bcstart(ib)
          maxlbl(ib) = bcend(ib)
          lchar(ib) = 'A'
          bcl(ib) = 0
          noss = bcend(ib) - bcstart(ib) + 1
          do iss = 1, noss
            indss = bcstart(ib)+iss-1
            write (*,*) ' surface index ',indss
!  find all faces belonging to surface structure INDSS
            call find_faces (indss, inbc, bcl(ib),
     >                       m%mxnbc,m%bcFc, m%bcFcOr, m, idb, nfaces)
          end do  ! iss
          m%bcFcP(ib,2) = bcl(ib)
        else
          call xerrab ('conv_bcparms -- wrong BCCHAR')
        endif
      enddo
      close (idb)

      acchar(1:nbc) = lchar(1:nbc)
      abcstart(1:nbc) = minlbl(1:nbc)
      abcend(1:nbc) = maxlbl(1:nbc)

      ! now check that all boundary faces are assigned to a BC
      call check_boundary_labels(m, nBc, abcstart(1:nbc), abcend(1:nbc),
     & .true.)

      if (inbc .ne. m%nCg ) then
        write (*,*) 'inbc, nCg ',inbc,m%nCg
! try to find missing or doubled faces
        allocate (ifchit(m%nFC))
        ifchit = 0
! for all guard cells mark face to plasma cell
        do i = 1, m%nCg
          iCv = m%nCi + i
          iFc = m%cvFc(m%cvFcP(iCv,1))
          ifchit(iFc) = ifchit(iFc) + 1
        end do
! remove mark from boundary faces
        do ib = 1, m%nbc
          do i = 1,  m%bcFcP(ib,2)
            iFc = m%bcFc(i)
            ifchit(iFc) = ifchit(iFc) - 1
          end do
        end do
        if (count(ifchit /= 0) > 0) then
          write (*,*) 'suspicious faces '
          do i = 1, m%nFc
            if (ifchit(i) /= 0) then
              write (*,*) 'iFc ',i,ifchit(i)
              write (*,*) 'fcCv ',m%fcCv(i,1),m%fcCv(i,2)
            end if
          end do
        end if

        call xerrab ('conv_bcparms -- problem setting up bc arrays')
      end if

! now convert feedback BC: BCCON 11/12-->8, BCENE/I/POT 10-->8
      do ib = 1, m%nBc
        if (bcene(ib) .eq. 10) bcene(ib) = 8
        if (bceni(ib) .eq. 10) bceni(ib) = 8
        if (bcpot(ib) .eq. 10) bcpot(ib) = 8
        do is = 0, ns-1
          if (bccon(is,ib) .eq. 11 .or. bccon(is,ib) .eq. 12)
     &      bccon(is,ib) = 8
        enddo
      enddo

      return
      end subroutine conv_bcparms

***********************************************************************

      subroutine write_bcparms (ns, m)
      use b2mod_version
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent (inout) :: m
      integer :: ib, iout, is
      character(1) :: ap =''''
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open(newunit=iout,file='b2us.boundary.parameters')
#else
      iout = newunit()
      open(unit=iout,file='b2us.boundary.parameters')
#endif
      call cfverw(iout,newversion)

      write (iout,'(1x,a)') '&boundary'
      write (iout,'(1x,a,i0)') 'NBC     = ',m%nBc
      write (iout,'(1x,a,(1x,T14,10(2x,4a1)))')
     .       'BCCHAR  =  ',(ap,acchar(ib),ap,',',ib=1,m%nBc)

      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'BCSTART =  ',(abcstart(ib),',',ib=1,m%nBc)
      write(iout,'(1x,a,(1x,T14,10(i5,a1)))')
     .       'BCEND  =  ',(abcend(ib),',',ib=1,m%nBc)

      write (iout,'(1x,a,(1x,T14,10(1x,i9,a1)))')
     .       'BCENE   =   ',(bcene(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENEPAR(1,1)=',(enepar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENEPAR(1,2)=',(enepar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(1x,T14,10(1x,i9,a1)))')
     .       'BCENI   =   ',(bceni(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENIPAR(1,1)=',(enipar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENIPAR(1,2)=',(enipar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(1x,T14,10(1x,i9,a1)))')
     .      'BCENK   =   ',(bcenk(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENKPAR(1,1)=',(enkpar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENKPAR(1,2)=',(enkpar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(1x,T14,10(1x,i9,a1)))')
     .      'BCENZ   =   ',(bcenz(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENZPAR(1,1)=',(enzpar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'ENZPAR(1,2)=',(enzpar(ib,2),',',ib=1,m%nBc)

      write (iout,'(1x,a,(1x,T14,10(1x,i9,a1)))')
     .      'BCPOT   =   ',(bcpot(ib),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'POTPAR(1,1)=',(potpar(ib,1),',',ib=1,m%nBc)
      write (iout,'(1x,a,(1x,T14,10(es10.2,a1)))')
     .       'POTPAR(1,2)=',(potpar(ib,2),',',ib=1,m%nBc)

      do ib = 1, m%nBc
        write (iout,'(1x,a,i2,a,(1x,T17,10(1x,i9,a1)))')
     .       'BCCON(0,',ib,') =  ',(bccon(is,ib),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(1x,T17,10(es10.2,a1)))')
     .       'CONPAR(0,',ib,',1)=',(conpar(is,ib,1),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(1x,T17,10(es10.2,a1)))')
     .       'CONPAR(0,',ib,',2)=',(conpar(is,ib,2),',',is=0,ns-1)
      end do

      do ib = 1, m%nBc
        write (iout,'(1x,a,i2,a,(1x,T17,10(1x,i9,a1)))')
     .       'BCMOM(0,',ib,') =  ',(bcmom(is,ib),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(1x,T17,10(es10.2,a1)))')
     .       'MOMPAR(0,',ib,',1)=',(mompar(is,ib,1),',',is=0,ns-1)
        write (iout,'(1x,a,i2,a,(1x,T17,10(es10.2,a1)))')
     .       'MOMPAR(0,',ib,',2)=',(mompar(is,ib,2),',',is=0,ns-1)
      end do

      write (iout,'(1x,a,f10.7,a1)') 'GAMMAI=',gammai,','
      write (iout,'(1x,a,f10.7,a1)') 'GAMMAE=',gammae,','
      write (iout,'(1x,a,l1,a1)') 'LBNDUSR=',lbndusr,','
      write (iout,'(1x,a,l1,a1)') 'LFEEDBACK=',lfeedback,','

      write (iout,'(1x,a,i0,a1)') 'NNISO = ',nniso,','

      if (nniso > 0)
     .   call xerrab ('write_bcparms -- nniso > 0')

      write (iout,'(1x,a,es10.2,a1)') 'TEISO = ',teiso,','
      write (iout,'(1x,a,es10.2,a1)') 'TIISO = ',tiiso,','
      write (iout,'(1x,a,es10.2,a1)') 'PHIISO = ',phiiso,','
      write (iout,'(1x,a,(10(es10.2,a1)/))')
     .       'NIISO(0)=',(niiso(is),',',is=0,ns)

      if (ncbs > 0)
     .  write (iout,'(1x,a,(10(1x,i5,a1)/))')
     .      'lcbs=',(lcbs(ib),',',ib=1,ncbs)

      write (iout,'(1x,5a)') 'BOUNDARY_FILENAME=',ap,
     .                       trim(boundary_filename),ap,','
      write (iout,'(1x,a,es10.2,a1)')
     .       'BOUNDARY_TIME_MOD = ',boundary_time_mod,','
      write (iout,'(1x,a,es10.2,a1)')
     .       'BOUNDARY_TIME_SWITCH = ',boundary_time_switch,','
      write (iout,'(1x,a,l1,a1)') 'WRITE_NML_BND=',write_nml_bnd,','

      write (iout,'(1x,a)') '/'
      close (iout)

      end subroutine write_bcparms


      subroutine replace_outdated_bcparms (ns, m)
      implicit none
      integer, intent(in) :: ns
      type (mapping), intent (inout) :: m
      integer :: ib, is
      logical changed

cc Here outdated BCs that will not be revived in the WG code can be automatically replaced by new versions.
cc This should only be done if the change is backward compatible! For example, the non-drift sheath BC BCMOM = 3 can be automaticaly replaced by
cc BCMOM = 13, because it is drift-compatible but gives the same result as the old BCMOM = 3 for cases without drifts.

      changed = .false.
      do ib=1,m%nBc
        do is=0,ns-1
cc BCMOM
          if (bcmom(is,ib).eq.3) then
            bcmom(is,ib) = 13
            write(*,*) 'Automatically changing BCMOM from 3 to 13'
            changed = .true.
          endif 

        enddo
cc BCENE
        if (bcene(ib).eq.3) then 
          bcene(ib) = 15
          write(*,*) 'Automatically changing BCENE from 3 to 15'
          changed = .true.
        endif 
        if (bcene(ib).eq.12) then 
          bcene(ib) = 15
          write(*,*) 'Automatically changing BCENE from 12 to 15'
          changed = .true.
        endif 
cc BCENI
        if (bceni(ib).eq.3) then 
          bceni(ib) = 15
          write(*,*) 'Automatically changing BCENI from 3 to 15'
          changed = .true.
        endif 
        if (bceni(ib).eq.12) then 
          bceni(ib) = 15
          write(*,*) 'Automatically changing BCENI from 12 to 15'
          changed = .true.
        endif 

cc BCPOT

        if (bcpot(ib).eq.3) then 
          bcpot(ib) = 11
          write(*,*) 'Automatically changing BCPOT from 3 to 11'
          changed = .true.
        endif 

      enddo

      if (changed) then
          write(*,*) 'Some outdated BCs' //
     .    ' in WG have been replaced automatically.'
          write(*,*) 'See lines above.'
      endif

      end subroutine replace_outdated_bcparms

***********************************************************************

      subroutine bc_to_struct (m)
      type(mapping), intent(inout) :: m
      integer :: ib, inbc, ibc, iFc, ic1, ic2
      integer :: iout, ian, ien, i, j, k
#ifndef F2003
      integer newunit
      external newunit
#endif

      m%nBc = nbc
!pb   m%mxnBc = sum(nbcfaces(1:nbc))       ! number of cells with imposed BC
!pb   m%mxnBc = sum(m%bcFcP(1:m%nbc,2))    ! number of cells with imposed BC
      m%mxnBc = count(m%fcLbl /= 0)

      call alloc_mapping_bc(m)

      inbc = 0 ! index in bcCv
      do ib = 1, m%nBc
        m%bcCvP(ib,1) = inbc + 1
!pb     m%bcCvP(ib,2) = nbcfaces(ib)
        m%bcCvP(ib,2) = m%bcFcP(ib,2)
!pb     do ibc = 1, nbcfaces(ib)
        do ibc = 1, m%bcFcP(ib,2)
          inbc = inbc + 1
!pb       iFc = bc_faces(ibc,ib)
          iFc = m%bcFc(m%bcFcP(ib,1)+ibc-1)
          ic1 = min(m%fcCv(iFc,1), m%fcCv(iFc,2))
          ic2 = max(m%fcCv(iFc,1), m%fcCv(iFc,2))
          if ((ic1 <= m%nCi) .and. (ic2 >m%nCi)) then
            m%bcCv(inbc,1) = ic2
            m%bcCv(inbc,2) = ic1
!pb         m%bcFcOr(inbc) = bc_face_ori(ibc,ib)
          else
            call xerrab('bc_to_struct -- problem with boundary face')
          end if
        end do
      end do

      if (allocated(m%strDiv)) then
        do i = 1, maxval(m%strDiv)
          do j = m%divFcP(i,1), m%divFcP(i,1) + m%divFcP(i,2) - 1
            iFc = m%divFc(j)
            do k = 1, size(m%bcFc)
              if (m%bcFc(k).eq.iFc) m%divFcOr(j) = m%bcFcOr(k)
            end do
          end do
        end do
      end if
#ifdef WG_DEBUG
#ifdef F2003
      open (newunit=iout,file='bc.testout')
#else
      iout = newunit()
      open (unit=iout,file='bc.testout')
#endif
      write (iout,*) 'nbc =',m%nBc
      write (iout,'(a,/(10i6))') 'bcCvP(1) ',(m%bcCvP(ib,1),ib=1,m%nBc)
      write (iout,'(a,/(10i6))') 'bcCvP(2) ',(m%bcCvP(ib,2),ib=1,m%nBc)
      do ib = 1, m%nBc
        ian = m%bcCvP(ib,1)
        ien = m%bcCvP(ib,1) + m%bcCvP(ib,2) - 1
        write (iout,*) 'ib = ',ib
        write (iout,'(a,/(10i6))') 'bcCv(1) ',(m%bcCv(i,1),i=ian,ien)
        write (iout,'(a,/(10i6))') 'bcCv(2) ',(m%bcCv(i,2),i=ian,ien)
        write (iout,'(a,/(10f6.2))') 'bcCv(2) ',(m%bcFcOr(i),i=ian,ien)
      end do
      close (iout)
#endif
      end subroutine bc_to_struct

      end module b2mod_boundary_namelist
