      subroutine conv_surface_models (nout, mpg, gm,
     .           bndmdl, trifc, elem_cv, nelem, node_vert, nnode)

      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_indirect
      use b2us_map
      use b2us_geo

      implicit none

      integer, intent(in) :: nout, nelem, nnode
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer, intent(in) :: trifc(3,nelem), elem_cv(nelem),
     .           node_vert(nnode)
      integer, intent(inout) :: bndmdl(3,nelem)
      character(500) :: line
      integer :: ineir, iouteir, nstsi, ists, i, ier, icut, nd, nsf,
     .           isf, j, ifc, ib, ntrg, ihelp(12), isum, nadmod, nasmod,
     .           itrg, ix, iy, ic, nt, it, is
      integer, allocatable :: idim(:), ipos(:), istart(:), iend(:)
      integer, allocatable :: ipos_b2(:), indelm(:), tris(:)
      integer, allocatable :: idim_eir(:), ipos_eir(:),
     .                        istart_eir(:), iend_eir(:), idtrg(:),
     .                        istart_dum(:), iend_dum(:)
      logical :: found,lhelp(60)
      integer :: artificial_slab=0, spatial_hybrid = 0
      integer :: find_loc
#ifndef F2003
      integer newunit
      external newunit
#endif
      intrinsic pack

      call ipgeti ('artificial_slab', artificial_slab)
      call ipgeti ('b2mn_spatial_hybrid', spatial_hybrid)

      nd = mpg%nRc + mpg%nBc
      if (artificial_slab > 0) nd = nd + 4*nncut
      allocate (idim(nd))
      allocate (ipos(nd))
      allocate (ipos_b2(nd))
      allocate (istart(nd))
      allocate (iend(nd))
      allocate (idtrg(nd))
      idtrg = 0

      call find_eirene_numbering(mpg, idim, ipos, istart, iend,
     .                           nd, nsf, ipos_b2)
#ifdef F2003
      open (newunit=ineir,file='input.dat',iostat=ier)
#else
      ineir = newunit()
      open (unit=ineir,file='input.dat',iostat=ier)
#endif
      call xertst (ier.eq.0, 'input.dat not found')
#ifdef F2003
      open (newunit=iouteir,file='input.dat.us')
#else
      iouteir = newunit()
      open (unit=iouteir,file='input.dat.us')
#endif

      do
        read (ineir,'(a500)') line
        write (iouteir,'(a)') trim(line)
        if (index(line,'*** 2') > 0) exit
      end do
!  input block 2 found
!  indgrd
      read (ineir,'(a500)') line
      write (iouteir,'(a)') trim(line)
!  nlrad
      read (ineir,'(a500)') line
      write (iouteir,'(a)') trim(line)

      read (ineir,'(12(5l1,1x))') lhelp(1:8)
      if (lhelp(5)) then
        lhelp(5) = .false.
        lhelp(6) = .true.
      end if
      write (iouteir,'(12(5l1,1x))') lhelp(1:8)

      do
        read (ineir,'(a500)') line
        write (iouteir,'(a)') trim(line)
        if (index(line,'*** 3') > 0) exit
      end do

!  input block 3 found
!      write (iouteir,'(a)') trim(line)
      read (ineir,'(i6)') nstsi
      write (iouteir,'(i6)') nstsi

      allocate (idim_eir(nstsi))
      allocate (ipos_eir(nstsi))
      allocate (istart_eir(nstsi))
      allocate (iend_eir(nstsi))
      allocate (istart_dum(nstsi))
      allocate (iend_dum(nstsi))
      ists = 0
      do
        do
          read (ineir,'(a500)') line
          write (iouteir,'(a)') trim(line)
          if (line(1:1) /= '*') exit
        end do
        ists = ists + 1
        read (line,'(12i6)') i,idim_eir(ists),ipos_eir(ists),
     .                         istart_eir(ists), iend_eir(ists),
     .                         istart_dum(ists), iend_dum(ists)
        !! Accommodate different input formats
        if (idim_eir(ists).eq.1 .and.
     .      istart_eir(ists).eq.iend_eir(ists) .and.
     .      istart_eir(ists).eq.ipos_eir(ists) .and.
     .      istart_dum(ists).ne.0 .and. iend_dum(ists).ne.0) then
           ! range specified in 'IRPTA2' and 'IRPTE2'
           istart_eir(ists) = istart_dum(ists)
           iend_eir(ists) = iend_dum(ists)
        endif
        do
          read (ineir,'(a500)') line
          write (iouteir,'(a)') trim(line)
          if (line(1:1) == '*') exit
        end do
        if (line(1:3) == '***') exit
      end do

!  finshed reading input block 3

!  find triangles on nondefault standard surfaces

      allocate (indelm(nelem))
      indelm = (/ (i, i=1,nelem) /)

      ntrg = 0
      do ists = 1, nstsi
        found = .false.

!  check neutral source regions and boundary surfaces for matching definitions
        do isf = 1, nsf
          if (idim(isf) < 0) cycle
          if (idim_eir(ists) == idim(isf)) then
            if ((ipos_eir(ists) == ipos(isf)) .and.
     .          (istart_eir(ists) == istart(isf)) .and.
     .          (iend_eir(ists) == iend(isf))) then
!  matching surfaces identified
              if ( isf <= mpg%nRc) then
!  neutral source
                found = .true.
                do i = 1, mpg%rcFcP(isf,2)
                  j = mpg%rcFcP(isf,1)+i-1
                  ifc = mpg%rcFc(j) ! face on neutral source
                  where (abs(trifc) == ifc)
                    bndmdl = -ists
                  end where
                end do
                ntrg = ntrg + 1
                idtrg(isf) = -ists
!                exit  ! surface marked
              elseif ( isf <= mpg%nRc + mpg%nBc ) then
!  boundary
                found = .true.
                ib = isf - mpg%nRc
                do i = 1, mpg%bcFcP(ib,2)
                  j = mpg%bcFcP(ib,1)+i-1
                  ifc = mpg%bcFc(j) ! face on boundary
                  where (abs(trifc) == ifc)
                    bndmdl = -ists
                  end where
                end do
                exit            ! surface marked
              else if (artificial_slab .eq. 1) then
                found = .true.
                ix = ipos_b2(isf)
                icut = (isf - mpg%nRc - mpg%nBc) / 4 + 1
                if (mod(isf - mpg%nRc - mpg%nBc,4) == 0) icut = icut - 1
                do iy = bottomcut(icut)+1, topcut(icut)-1
                  if (ix <= leftcut(icut)) then
                    ifc = mpg%imapFcx(leftcut(icut), iy)
                  else
                    ifc = mpg%imapFcx(rightcut(icut), iy)
                  end if
                  ic = mpg%imapCv(ix,iy)
                  nt = count(elem_cv == ic)
                  if (nt > 0) then
                    allocate(tris(nt))
                    tris = pack(indelm, elem_cv == ic)
                    do it = 1, nt
                      do is = 1, 3
                        if (abs(trifc(is,tris(it))) == ifc) then
                          bndmdl(is,tris(it)) = -ists
                          exit
                        end if
                      end do    ! is
                    end do      ! it
                    deallocate (tris)
                  end if
                end do          ! iy
              end if
            end if
          end if
        end do                  ! isf

        if (.not. found) then
          write (*,*) ' Eirene nondefault standard surface ',ists,
     .                ' could not be identified '
        end if

      end do

!  copy input file until block 14
      do
        read (ineir,'(a500)') line
        write (iouteir,'(a)') trim(line)
        if (index(line,'*** 14') > 0) exit
      end do

!  copy line of logicals
      read (ineir,'(a500)') line
      write (iouteir,'(a)') trim(line)

!  nfla ntrfrm ...
!  ntrfrm=1
      read (ineir,'(12i6)') ihelp
      write (iouteir,'(12i6)') ihelp(1),1,ihelp(5:12)

!  iflb ...
      do
        read (ineir,'(a500)') line
        if (scan(line,'eEdD-') > 0) then
          write (iouteir,'(a)') trim(line)
        else
!  ndxa, ndya
          write (iouteir,'(3i6)') mpg%nCv, mpg%nCi, mpg%nFc
         exit
        end if
      end do

!  number of targets
      read (ineir,'(12i6)') ihelp(1)
      itrg = 0
      do i = 1, ntrg
        if (index('EWSN',crcstra(i)) > 0) itrg = itrg + 1
      end do

      if (spatial_hybrid.eq.0) then !! default, no spatially hybrid
        if (ihelp(1) /= itrg)
     .      call xerrab ('number of targets mismatch')
      else !! spatially hybrid neutrals are used
        !! SH_TODO
      endif
      write (iouteir,'(12i6)') itrg
!  number of target parts
      isum = 0
      do i = 1, ntrg, 12
        ihelp = 0
        read (ineir,'(12i6)') ihelp
        write (iouteir,'(12i6)') ihelp
        isum = isum + sum(ihelp)
      end do
!  write target information
      do i = 1, isum
        read (ineir,'(12i6)') ihelp
        if (i <= ntrg) then
          write (iouteir,'(12i6)') ihelp(1), idtrg(i), ihelp(8:12)
        else
!  nothing to be done
        end if
      end do

!  copy input file until block 15 or end of file
      do
        read (ineir,'(a500)',end=99) line
        write (iouteir,'(a)') trim(line)
        if (index(line,'*** 15') > 0) exit
      end do

!  adapt block 15
      read (ineir,'(12i6)') nadmod, nasmod
      write (iouteir,'(12i6)') nadmod, nasmod

      do i = 1, nadmod
        read (ineir,'(a500)') line
        write (iouteir,'(a)') trim(line)
      end do

      if (nasmod /= 0) then
        if (nasmod.eq.2 .and. count(mpg%icornVx > 0) == 4) then
          !! limiter or linear case: only NW and NE corners
          j = find_loc(node_vert,nnode,mpg%icornVx(2))
          read (ineir,'(2i6)') ihelp(1), ihelp(2)
          write (iouteir,'(3i6)') ihelp(1), ihelp(2), j
          j = find_loc(node_vert,nnode,mpg%icornVx(4))
          read (ineir,'(2i6)') ihelp(1), ihelp(2)
          write (iouteir,'(3i6)') ihelp(1), ihelp(2), j
        elseif (count(mpg%icornVx > 0) == nasmod) then
          !! SN of DN case: 4 or 8 corners
          do i = 1, nasmod
            j = find_loc(node_vert,nnode,mpg%icornVx(i))
            read (ineir,'(2i6)') ihelp(1), ihelp(2)
            write (iouteir,'(3i6)') ihelp(1), ihelp(2), j
          end do
        else
         call xerrab('conv_surface_models -- nasmod')
        endif
      else
         write (*,*) "Warning: found NASMOD = 0 in input block 15"
      endif

!  copy input file until end of file
      do
        read (ineir,'(a500)',end=99) line
        write (iouteir,'(a)') trim(line)
        if (index(line,'*** 15') > 0) exit
      end do

 99   continue
      close (iouteir)
      close (ineir)

      deallocate (idim)
      deallocate (ipos)
      deallocate (ipos_b2)
      deallocate (istart)
      deallocate (iend)
      deallocate (idim_eir)
      deallocate (ipos_eir)
      deallocate (istart_eir)
      deallocate (iend_eir)
      deallocate (istart_dum)
      deallocate (iend_dum)
      deallocate (indelm)

      end subroutine conv_surface_models

      subroutine find_nondef_surfaces (nout, mpg, gm,
     .           elem, bndmdl, elem_cv,nelem, node_vert, nnode)

      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_indirect
      use b2us_map
      use b2us_geo

      implicit none

      integer, intent(in) :: nout, nelem, nnode
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      integer, intent(in) :: elem(3,nelem), elem_cv(nelem),
     .                       node_vert(nnode)
      integer, intent(inout) :: bndmdl(3,nelem)
      integer, allocatable :: indelm(:), tris(:), jsts(:)
      integer :: istra, imodel, i, j, ifc, ifcmin, ifcmax, ic, nt, it,
     .           jt, is, is1, itmin, itmax, ibold, ib,
     .           ntarg, nblk3, nblk14, nd, nsf, istb
      integer, allocatable :: idim(:), ipos(:), istart(:), iend(:)
      integer, allocatable :: ipos_b2(:)
      integer :: artificial_slab=0
      logical :: exists
#ifndef F2003
      integer newunit
      external newunit
#endif

      call ipgeti ('artificial_slab', artificial_slab)

!  set recycling boundary model indicators
      bndmdl = 0

      allocate (indelm(nelem))
      indelm = (/ (i, i=1,nelem) /)

      imodel = 0
      do istra = 1, mpg%nRc
!  skip neutral sources without associated faces
!  CRCSTRA = 'V','C','T'
        if (mpg%rcFcP(istra,2) == 0) cycle
!nh     check if the surface is not already written for a different species
        exists = .false.
        do i = 1,istra-1
          if (rcstart(istra).eq.rcstart(i)) exists = .true.
        end do
        if (exists) cycle
        imodel = imodel + 1
        write (nout,'(a15,i6,a,i6)') 'recyc. surface ',istra,
     .   ' is associated with nondefault standard surface ',imodel
        do i = 1, mpg%rcFcP(istra,2)
          j = mpg%rcFcP(istra,1)+i-1
          ifc = mpg%rcFc(j) ! face on neutral source
          ifcmin = min(mpg%fcVx(ifc,1),mpg%fcVx(ifc,2))
          ifcmax = max(mpg%fcVx(ifc,1),mpg%fcVx(ifc,2))
          ic = mpg%rcCv(j,2) ! cell inside plasma
!  identify triangle side that represents the face ifc
!  therefore check all triangles of former cell ic
          nt = count(elem_cv == ic)
          if (nt > 0) then
            allocate(tris(nt))
            tris = pack(indelm, elem_cv == ic)
            itloop: do it = 1, nt
              jt = tris(it)
              do is = 1, 3
                is1 = is + 1
                if (is1 > 3) is1 = 1
                itmin = min(node_vert(elem(is,jt)),
     .                      node_vert(elem(is1,jt)))
                itmax = max(node_vert(elem(is,jt)),
     .                      node_vert(elem(is1,jt)))
                if ((ifcmin == itmin) .and. (ifcmax == itmax)) then
                  bndmdl(is,jt) = -imodel
                  exit itloop
                end if
              end do
            end do itloop
            deallocate (tris)
          end if
        end do
        
      end do

      ntarg = imodel

!  now check for triangle sides on non recycling boundaries

      ibold = -1
      do ib = 1, mpg%nBc
        do i = 1, mpg%bcFcP(ib,2)
          j = mpg%bcFcP(ib,1)+i-1
          ifc = mpg%bcFc(j) ! face on boundary
          ifcmin = min(mpg%fcVx(ifc,1),mpg%fcVx(ifc,2))
          ifcmax = max(mpg%fcVx(ifc,1),mpg%fcVx(ifc,2))
          ic = mpg%bcCv(j,2) ! cell inside plasma
!  identify triangle side that represents the face ifc
!  therefore check all triangles of former cell ic
          nt = count(elem_cv == ic)
          if (nt > 0) then
            allocate(tris(nt))
            tris = pack(indelm, elem_cv == ic)
            ibloop: do it = 1, nt
              jt = tris(it)
              do is = 1, 3
                is1 = is + 1
                if (is1 > 3) is1 = 1
                itmin = min(node_vert(elem(is,jt)),
     .                      node_vert(elem(is1,jt)))
                itmax = max(node_vert(elem(is,jt)),
     .                      node_vert(elem(is1,jt)))
                if ((ifcmin == itmin) .and. (ifcmax == itmax)) then
                  if (bndmdl(is,jt) == 0) then
                    if (ib /= ibold) then
                      imodel = imodel + 1
                      ibold = ib
                      write (nout,'(a15,i6,a,i6)') 'bound. surface ',ib,
     .   ' is associated with nondefault standard surface ',imodel
                   end if
                    bndmdl(is,jt) = -imodel
                  end if
                  exit ibloop
                end if
              end do
            end do ibloop
            deallocate (tris)
          end if
        end do
      end do

      deallocate (indelm)

!  write nondefault standard surface block for Eirene input

#ifdef F2003
      open (newunit=nblk3,file='eirene_block3.dat')
#else
      nblk3 = newunit()
      open (unit=nblk3,file='eirene_block3.dat')
#endif
      write (nblk3,'(a)')
     .  '*** 3A DATA FOR NONDEFAULT STANDARD SURFACES, NSTSI='
      write (nblk3,'(i6)') imodel

      nd = mpg%nRc + mpg%nBc
      if (artificial_slab > 0) nd = nd + 4*nncut
      allocate (idim(nd))
      allocate (ipos(nd))
      allocate (ipos_b2(nd))
      allocate (istart(nd))
      allocate (iend(nd))

      call find_eirene_numbering(mpg, idim, ipos, istart, iend,
     .                           nd, nsf, ipos_b2)

!  ... from neutral sources
      imodel = 0
      allocate (jsts(ntarg))
      jsts = 0

      do istra = 1, mpg%nRc
        if (idim(istra) < 0) cycle
!nh     check if the surface is not already written for a different species
        exists = .false.
        do i = 1,istra-1
          if (rcstart(istra).eq.rcstart(i)) exists=.true.
        end do
        if (exists) cycle
        imodel = imodel + 1
        jsts(imodel) = -imodel

        write (nblk3,'(a1,i3,a2,a)')
     .    '*',imodel,crcstra(istra),' surface from B2.5'
        write (nblk3,'(5i6)')
     .    imodel, idim(istra), ipos(istra), istart(istra), iend(istra)
        write (nblk3,'(12i6)') 1, 2, 0, 0, 0, imodel, 0, 0, 0
        write (nblk3,'(12i6)') 1, 0
        write (nblk3,'(6es12.4)') 1206., 2., 0., 0., 0., 2.8
        write (nblk3,'(6es12.4)') 1., 1., 1., 1., 0.5, 1.
      end do

!  ... from boundaries
      istb = mpg%nRc
      bound: do ib = 1, mpg%nBc
        do istra = 1, mpg%nRc
          if ((bcchar(ib) == crcstra(istra)) .and.
     .        (bcpos(ib) == rcpos(istra)) .and.
     .        (bcstart(ib) <= rcstart(istra)) .and.
     .        (bcend(ib) >= rcend(istra))) cycle bound
        end do

        imodel = imodel + 1
        istb = istb + 1
        write (nblk3,'(a1,i3,a2,a)')
     .    '*',imodel,bcchar(ib),' boundary from B2.5'
        write (nblk3,'(5i6)')
     .    imodel, idim(istb), ipos(istb), istart(istb), iend(istb)
        write (nblk3,'(12i6)') 2, 2, 0, 0, 0, imodel, 0, 0, 0
        write (nblk3,'(12i6)') 1, 0
        write (nblk3,'(6es12.4)') 1206., 2., 0., 0., 0., 2.8
        write (nblk3,'(6es12.4)') 1., 1., 1., 1., 0.5, 1.
      end do bound

      close (nblk3)

!  write target information for Eirene input block 14

#ifdef F2003
      open (newunit=nblk14,file='eirene_block14.dat')
#else
      nblk14 = newunit()
      open (unit=nblk14,file='eirene_block14.dat')
#endif
      write (nblk14,'(3i6)') mpg%nCv, mpg%nCi, mpg%nFc
!nh   determine the number of surface strata
      ntarg = 0
      imodel = 0
      do istra = 1, mpg%nRc
!  skip neutral sources without associated faces
!  CRCSTRA = 'V','C','T'
        if (mpg%rcFcP(istra,2) == 0) cycle
        ntarg = ntarg+1
      end do
      write (nblk14,'(i6)') ntarg
      write (nblk14,'(12i6)') (1, i=1,ntarg)
      do i = 1, ntarg
        imodel = imodel + 1
        do j = 1, i-1
          if (rcstart(i).eq.rcstart(j)) then
            imodel = imodel - 1 ! surface is already defined
            exit
          end if
        end do
        write (nblk14,'(12i6)') i, -imodel, 0, 0, 1, 1
      end do
      close (nblk14)

      deallocate(jsts)
      deallocate (idim)
      deallocate (ipos)
      deallocate (istart)
      deallocate (iend)

      return
      end subroutine find_nondef_surfaces


      subroutine find_eirene_numbering (mpg,
     .           idim, ipos, istart, iend, nd, nsf, ipos_b2)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_indirect
      use b2us_map
      use b2us_geo

      implicit none

      type(mapping), intent(in) :: mpg
      integer, intent(in) :: nd
      integer, intent(out) :: nsf
      integer, intent(out) :: idim(nd), ipos(nd), istart(nd),
     .                        iend(nd), ipos_b2(nd)
      integer :: icut, istra, ib, ix
      integer :: artificial_slab=0

      nsf = 0
      idim = -1

! find grid lines numbers for neutral source surfaces
      do istra = 1, mpg%nRc
        nsf = nsf + 1
        if (index('WESNA',crcstra(istra)) == 0) cycle
        select case (crcstra(istra))
        case ('S')
          idim(nsf) = 1
          ipos(nsf) = rcpos(istra) + 2
          ipos_b2(nsf) = rcpos(istra)
          istart(nsf) = rcstart(istra) + 1
          iend(nsf) = rcend(istra) + 2
          do icut = 1, nncut
            if (rcstart(istra) >= leftcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcstart(istra) >= rightcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcend(istra) >= leftcut(icut))
     .        iend(nsf) = iend(nsf) + 1
            if (rcend(istra) >= rightcut(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (rcstart(istra) >= mpg%ictcntrl(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcend(istra) >= mpg%ictcntrl(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
        case ('N')
          idim(nsf) = 1
          ipos(nsf) = rcpos(istra) + 1
          ipos_b2(nsf) = rcpos(istra)
          istart(nsf) = rcstart(istra) + 1
          iend(nsf) = rcend(istra) + 2
          do icut = 1, nncut
            if (rcstart(istra) >= leftcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcstart(istra) >= rightcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcend(istra) >= leftcut(icut))
     .        iend(nsf) = iend(nsf) + 1
            if (rcend(istra) >= rightcut(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (rcstart(istra) >= mpg%ictcntrl(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (rcend(istra) >= mpg%ictcntrl(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
        case ('W')
          idim(nsf) = 2
          ipos(nsf) = rcpos(istra) + 2
          ipos_b2(nsf) = rcpos(istra)
          istart(nsf) = rcstart(istra) + 1
          iend(nsf) = rcend(istra) + 2
          do icut = 1, nncut
            if (rcpos(istra) > leftcut(icut))
     .        ipos(nsf) = ipos(nsf) + 1
            if (rcpos(istra) >= rightcut(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (rcpos(istra) >= mpg%ictcntrl(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
        case ('E')
          idim(nsf) = 2
          ipos(nsf) = rcpos(istra) + 1
          ipos_b2(nsf) = rcpos(istra)
          istart(nsf) = rcstart(istra) + 1
          iend(nsf) = rcend(istra) + 2
          do icut = 1, nncut
            if (rcpos(istra) > leftcut(icut))
     .        ipos(nsf) = ipos(nsf) + 1
            if (rcpos(istra) >= rightcut(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (rcpos(istra) >= mpg%ictcntrl(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
        case ('A')
          idim(nsf) = 1
          ipos(nsf) = 1
          ipos_b2(nsf) = 1
          istart(nsf) = 1
          iend(nsf) = 1
        case ('C','V','T')
        case default
          call xerrab ('Wrong CRCSTRA or not yet available for WG')
        end select
      end do

!  ... from boundaries
      bound: do ib = 1, mpg%nBc
        nsf = nsf + 1
        select case (bcchar(ib))
        case ('S')
          idim(nsf) = 1
          ipos(nsf) = bcpos(ib) + 2
          ipos_b2(nsf) = bcpos(ib)
          istart(nsf) = bcstart(ib) + 1
          iend(nsf) = bcend(ib) + 2
          do icut = 1, nncut
            if (bcstart(ib) >= leftcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcstart(ib) >= rightcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcend(ib) >= leftcut(icut)) iend(nsf) = iend(nsf) + 1
            if (bcend(ib) >= rightcut(icut)) iend(nsf) = iend(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (bcstart(ib) >= mpg%ictcntrl(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcend(ib) >= mpg%ictcntrl(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
        case ('N')
          idim(nsf) = 1
          ipos(nsf) = bcpos(ib) + 1
          ipos_b2(nsf) = bcpos(ib)
          istart(nsf) = bcstart(ib) + 1
          iend(nsf) = bcend(ib) + 2
          do icut = 1, nncut
            if (bcstart(ib) >= leftcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcstart(ib) >= rightcut(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcend(ib) >= leftcut(icut)) iend(nsf) = iend(nsf) + 1
            if (bcend(ib) >= rightcut(icut)) iend(nsf) = iend(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (bcstart(ib) >= mpg%ictcntrl(icut))
     .        istart(nsf) = istart(nsf) + 1
            if (bcend(ib) >= mpg%ictcntrl(icut))
     .        iend(nsf) = iend(nsf) + 1
          end do
        case ('W')
          idim(nsf) = 2
          ipos(nsf) = bcpos(ib) + 2
          ipos_b2(nsf) = bcpos(ib)
          istart(nsf) = bcstart(ib) + 1
          iend(nsf) = bcend(ib) + 2
          do icut = 1, nncut
            if (bcpos(ib) > leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (bcpos(ib) >= rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (bcpos(ib) >= mpg%ictcntrl(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
        case ('E')
          idim(nsf) = 2
          ipos(nsf) = bcpos(ib) + 1
          ipos_b2(nsf) = bcpos(ib)
          istart(nsf) = bcstart(ib) + 1
          iend(nsf) = bcend(ib) + 2
          do icut = 1, nncut
            if (bcpos(ib) > leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (bcpos(ib) >= rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
          do icut = 1, mpg%nctcnt
            if (bcpos(ib) >= mpg%ictcntrl(icut))
     .        ipos(nsf) = ipos(nsf) + 1
          end do
        case ('A')
          idim(nsf) = 1
          ipos(nsf) = 1
          ipos_b2(nsf) = 1
          istart(nsf) = 1
          iend(nsf) = 1
        case default
          call xerrab ('Wrong BCCHAR or not yet available for WG')
        end select
      end do bound

!     ... from cuts in case of artificial_slab
      call ipgeti ('artificial_slab', artificial_slab)

      if (artificial_slab > 0) then
        do ib = 1, nncut
!  right face of cell left of left cut
          nsf = nsf + 1
          idim(nsf) = 2
          ix = leftcut(ib) - 1
          ipos(nsf) = ix + 2
          ipos_b2(nsf) = ix
          istart(nsf) = bottomcut(ib) + 2
          iend(nsf) = topcut(ib) + 1
          do icut = 1, nncut
            if (ix >= leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (ix >= rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
!  left face of cell right of left  cut
          nsf = nsf + 1
          idim(nsf) = 2
          ix = leftcut(ib)
          ipos(nsf) = ix + 2
          ipos_b2(nsf) = ix
          istart(nsf) = bottomcut(ib) + 2
          iend(nsf) = topcut(ib) + 1
          do icut = 1, nncut
            if (ix > leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (ix >= rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
!  right face of cell left of right cut
          nsf = nsf + 1
          idim(nsf) = 2
          ix = rightcut(ib) - 1
          ipos(nsf) = ix + 2
          ipos_b2(nsf) = ix
          istart(nsf) = bottomcut(ib) + 2
          iend(nsf) = topcut(ib) + 1
          do icut = 1, nncut
            if (ix >= leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (ix >= rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
!  left face of cell right of right  cut
          nsf = nsf + 1
          idim(nsf) = 2
          ix = rightcut(ib)
          ipos(nsf) = ix + 2
          ipos_b2(nsf) = ix
          istart(nsf) = bottomcut(ib) + 2
          iend(nsf) = topcut(ib) + 1
          do icut = 1, nncut
            if (ix > leftcut(icut)) ipos(nsf) = ipos(nsf) + 1
            if (ix > rightcut(icut)) ipos(nsf) = ipos(nsf) + 1
          end do
        end do
      end if

      return
      end subroutine find_eirene_numbering


