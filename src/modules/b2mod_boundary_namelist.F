      module b2mod_boundary_namelist
      use b2mod_types
      use b2mod_indirect
cank-20020724 {
c* Parameters controlling the source scaling procedure are added:
c*  lcbs          : boundary segments representing the core boundary in B2
c*  ncbs          : number of these boundary segments
cank }
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
     *   lcbs(:)
      character*1, allocatable, save ::
     *   BCCHAR(:)
#else
#     include "DIMENSIONS.F"
      integer NBCD, NNISOD, NSDMAX
      parameter (NBCD=20,NNISOD=10,NSDMAX=DEF_NSD)
      real (kind=R8), save :: CONPAR(0:nsdmax-1,NBCD,2)
      real (kind=R8), save :: MOMPAR(0:nsdmax-1,NBCD,2)
      real (kind=R8), save :: ENEPAR(NBCD,2)
      real (kind=R8), save :: ENIPAR(NBCD,2)
      real (kind=R8), save :: POTPAR(NBCD,2)
      integer, save :: BCPOS(NBCD)
      integer, save :: BCSTART(NBCD)
      integer, save :: BCEND(NBCD)
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

      integer, save :: NBC, NNISO, ncbs
      real (kind=R8) , save :: gammai, gammae, TEISO, TIISO, PHIISO
      logical, save :: lbndusr,lfeedback
      character*256, save :: boundary_filename
      real (kind=R8), save :: boundary_time_mod
      real (kind=R8), save :: boundary_time_switch

      contains

      subroutine alloc_b2mod_boundary(nsd,nbcd_in,nnisod_in)
      integer nsd, nbcd_in, nnisod_in
#ifndef ALLOCATE_AND_NAMELIST
      external xertst
#endif

#ifdef ALLOCATE_AND_NAMELIST
      NBCD=nbcd_in
      NNISOD=nnisod_in
      allocate(CONPAR(0:nsd-1,NBCD,2))
      allocate(MOMPAR(0:nsd-1,NBCD,2))
      allocate(ENEPAR(NBCD,2))
      allocate(ENIPAR(NBCD,2))
      allocate(POTPAR(NBCD,2))
      allocate(BCPOS(NBCD))
      allocate(BCSTART(NBCD))
      allocate(BCEND(NBCD))
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
      call xertst (nsd.le.nsdmax, 'faulty parameter nsd')
      call xertst (nbcd_in.le.nbcd, 'faulty parameter nbcd_in')
      call xertst (nnisod_in.le.nnisod, 'faulty parameter nnisod_in')
#endif
      nbc=0
      ncbs=0
      lcbs=0
      bcchar=' '
      bcpos=-2
      bcstart=-2
      bcend=-2
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

      end subroutine alloc_b2mod_boundary

      subroutine read_b2mod_boundary_namelist(nx, ny, ns, dotest)
      integer nx, ny, ns
      logical dotest
      integer is, ib, ic, ix, iy
      logical covered, guard_cell, found
      external xertst, xerrab

      namelist /boundary/ NBC,CONPAR,MOMPAR,ENEPAR,ENIPAR,POTPAR,
     1 BCPOS,BCSTART,BCEND,BCCON,BCMOM,BCENE,BCENI,BCPOT,
     2 gammai,gammae,lbndusr,lfeedback,
     3 NNISO,NIISO,TEISO,TIISO,PHIISO,NXISO1,NXISO2,NYISO1,NYISO2,
     4 bcchar,
     5 boundary_filename,boundary_time_mod,boundary_time_switch,
     6 lcbs
      logical :: file_ok
      character*260 filename

      filename='../' // boundary_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist boundary from ',trim(filename)
        read(99,boundary)
        write(*,*) 'Read namelist boundary from ',trim(filename)
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
          call xertst (1.le.bcpos(ib).and.bcpos(ib).le.ny+1,
     .     'faulty input bcpos (North)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .     'faulty input bcstart (North)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .     'faulty input bcend (North)')
        else if (bcchar(ib).eq.'S') then
          call xertst (-1.le.bcpos(ib).and.bcpos(ib).le.ny-1,
     .     'faulty input bcpos (South)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .     'faulty input bcstart (South)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .     'faulty input bcend (South)')
        else if (bcchar(ib).eq.'E') then
          call xertst (1.le.bcpos(ib).and.bcpos(ib).le.nx+1,
     .     'faulty input bcpos (East)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .     'faulty input bcstart (East)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .     'faulty input bcend (East)')
        else if (bcchar(ib).eq.'W') then
          call xertst (-1.le.bcpos(ib).and.bcpos(ib).le.nx-1,
     .     'faulty input bcpos (West)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .     'faulty input bcstart (West)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .     'faulty input bcend (West)')
        else if (bcchar(ib).eq.'X') then
          call xertst (0.le.bcpos(ib).and.bcpos(ib).le.ny,
     .     'faulty input bcpos (X)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.nx+1,
     .     'faulty input bcstart (X)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.nx+1,
     .     'faulty input bcend (X)')
        else if (bcchar(ib).eq.'Y') then
          call xertst (0.le.bcpos(ib).and.bcpos(ib).le.nx,
     .     'faulty input bcpos (Y)')
          call xertst (-1.le.bcstart(ib).and.bcstart(ib).le.ny+1,
     .     'faulty input bcstart (Y)')
          call xertst (-1.le.bcend(ib).and.bcend(ib).le.ny+1,
     .     'faulty input bcend (Y)')
        endif
        call xertst (bcstart(ib).le.bcend(ib),
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
        call xertst (-1.le.nxiso1(ib).and.nxiso1(ib).le.nx+1,
     .   'faulty input nxiso1')
        call xertst (-1.le.nxiso2(ib).and.nxiso2(ib).le.nx+1,
     .   'faulty input nxiso2')
        call xertst (-1.le.nyiso1(ib).and.nyiso1(ib).le.ny+1,
     .   'faulty input nyiso1')
        call xertst (-1.le.nyiso2(ib).and.nyiso2(ib).le.ny+1,
     .   'faulty input nyiso2')
        do ix = nxiso1(ib), nxiso2(ib)
          do iy = nyiso1(ib), nyiso2(ib)
            region(ix,iy,0) = 0
            resignore(ix,iy,1)=0
            resignore(ix,iy,2)=0
          enddo
        enddo
      enddo
c..  test that all boundary cells to the real domain are covered at least once
      do ix = -1, nx
        do iy = -1, ny
          guard_cell = region(ix,iy,0).ne.0 .and.
     &     (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &      bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)
          covered = .false.
          do ib = 1, nbc
            if (bcchar(ib).eq.'N' .or. bcchar(ib).eq.'S') 
     >       covered = covered .or. (iy.eq.bcpos(ib) .and.
     &        ix.ge.bcstart(ib) .and. ix.le.bcend(ib))
            if (bcchar(ib).eq.'E' .or. bcchar(ib).eq.'W') 
     >       covered = covered .or. (ix.eq.bcpos(ib) .and.
     &        iy.ge.bcstart(ib) .and. iy.le.bcend(ib))
          enddo
          if (guard_cell .and. .not.covered) then
            write(*,*) 'Guard cell',ix,iy,
     &       'not included in b2.boundary.parameters!'
            call xerrab('Incomplete boundary definition')
          endif
        enddo
      enddo
c..  test that boundary cells are guard cells
      do ib = 1, nbc
        if (bcchar(ib).eq.'N') then
          do ix = bcstart(ib), bcend(ib)
            call xertst(topiy(ix,bcpos(ib)).eq.ny+1, 
     &       'inconsistent North boundary!')
          enddo
        endif
        if (bcchar(ib).eq.'S') then
          do ix = bcstart(ib), bcend(ib)
            call xertst(bottomiy(ix,bcpos(ib)).eq.-2, 
     &       'inconsistent South boundary!')
          enddo
        endif
        if (bcchar(ib).eq.'E') then
          do iy = bcstart(ib), bcend(ib)
            call xertst(rightix(bcpos(ib),iy).eq.nx+1, 
     &       'inconsistent East boundary!')
          enddo
        endif
        if (bcchar(ib).eq.'W') then
          do iy = bcstart(ib), bcend(ib)
            call xertst(leftix(bcpos(ib),iy).eq.-2, 
     &       'inconsistent West boundary!')
          enddo
        endif
      enddo
c.. Test consistency with region boundaries
      do ic=1, nncut
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          found = bcend(ib).eq.leftcut(ic)-1 .and.
     &            bcchar(ib).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          found = bcstart(ib).eq.leftcut(ic) .and.
     &            bcchar(ib).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          found = bcend(ib).eq.rightcut(ic)-1 .and.
     &            bcchar(ib).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
        found = .false.
        ib = 0
        do while (ib.lt.nbc .and. .not.found)
          ib = ib + 1
          found = bcstart(ib).eq.rightcut(ic) .and.
     &            bcchar(ib).eq.'S'
        end do
        call xertst (found, 
     .   'Mismatch with region boundaries (South) !')
      end do
      if (nncut.ge.2) then
        do ic = 1, nniso
          found = .false.
          ib = 0
          do while (ib.lt.nbc .and. .not.found)
            ib = ib + 1
            found = bcchar(ib).eq.'S' .and. 
     &       (bcend(ib).eq.nxiso1(ic)-1.or.bcend(ib).eq.nxiso1(ic))
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (South) !')
          found = .false.
          ib = 0
          do while (ib.lt.nbc .and. .not.found)
            ib = ib + 1
            found = bcchar(ib).eq.'N' .and. 
     &       (bcend(ib).eq.nxiso1(ic)-1.or.bcend(ib).eq.nxiso1(ic))
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (North)  !')
          found = .false.
          ib = 0
          do while (ib.lt.nbc .and. .not.found)
            ib = ib + 1
            found = bcchar(ib).eq.'S' .and. 
     &       (bcstart(ib).eq.nxiso2(ic).or.bcstart(ib).eq.nxiso2(ic)+1)
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (South) !')
          found = .false.
          ib = 0
          do while (ib.lt.nbc .and. .not.found)
            ib = ib + 1
            found = bcchar(ib).eq.'N' .and. 
     &       (bcstart(ib).eq.nxiso2(ic).or.bcstart(ib).eq.nxiso2(ic)+1)
          end do
          call xertst (found, 
     .     'Mismatch with isolating boundaries (North)  !')
        end do
      end if
c.. Make a list of core boundaries
      ncbs=0
      lcbs=0
      do ib=1,nbc
        if(bcchar(ib).eq.'S') then                             !srv 28.11.08 {
         if(mod(region(bcstart(ib),bcpos(ib),0),4).eq.1 .and.
     .      mod(region(bcend(ib),bcpos(ib),0),4).eq.1) then
             ncbs=ncbs+1
             lcbs(ncbs)=ib
         endif
	  endif                                                  !srv 28.11.08 }
      enddo 
      call xertst (1.le.ncbs.or.ny.eq.1, 'faulty parameter ncbs')
      do ib = 1, ncbs
        call xertst (1.le.lcbs(ib).and.lcbs(ib).le.nbc, 
     .   'faulty parameter lcbs')
      enddo
      return

      end subroutine read_b2mod_boundary_namelist

      subroutine write_b2mod_boundary_namelist

      namelist /boundary/ NBC,CONPAR,MOMPAR,ENEPAR,ENIPAR,POTPAR,
     1 BCPOS,BCSTART,BCEND,BCCON,BCMOM,BCENE,BCENI,BCPOT,
     2 gammai,gammae,lbndusr,lfeedback,
     3 NNISO,NIISO,TEISO,TIISO,PHIISO,NXISO1,NXISO2,NYISO1,NYISO2,
     4 bcchar,
     5 boundary_filename,boundary_time_mod,boundary_time_switch,
     6 lcbs

      write(*,boundary)
      end subroutine write_b2mod_boundary_namelist

      end module b2mod_boundary_namelist
