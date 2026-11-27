      subroutine b2pedges(lun)
c
c  version : 31.03.2005 17:35
c=======================================================================
c***  Produces a table of ion fluxes and plasma temperatures at the
c***  outer edge of the grid for sputtering calculations
c***  for G. Strohmayer
c***
c*** Output files:
c***
c***  gridedge.dat with the data for G. Strohmayer
c=======================================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_elements
      use b2mod_constants
      use b2us_data
      implicit none
      integer lun

c*** Local variables
c***
c***  nf_hlp  : number of edge faces
c***  hlp_flx : fluxes of different species across the cell boundary
c***  mz_hlp  : ion charge
c***  hlp_al  : species name
c***  hlp_area: area of the corresponding cell face

      integer nf_hlp
      integer, allocatable :: mz_hlp(:)
      integer, allocatable :: old_face_list(:), new_face_list(:),
     *                        tmp_face_list(:)
      integer :: i,j,m,is,iCv,iCn,iFc,inew,iVa,iVb,ifmin,flmin
      real (kind=R8), allocatable :: hlp_flx(:), hlp_fly(:)
      real (kind=R8) :: ra,rb,re,za,zb,ze,zz,z,d,hlp_area
      character hlp_frmr*100
      character*2, allocatable :: hlp_al(:)
      logical :: left_match_found, right_match_found

c=======================================================================

      allocate(mz_hlp(0:state%ns-1),hlp_al(0:state%ns-1))
      allocate(hlp_flx(0:state%ns-1),hlp_fly(0:state%ns-1))
      do is = 0, state%ns-1
        mz_hlp(is) = nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is) = is_codes(is)
      end do
      allocate(old_face_list(mpg%nCg))
      allocate(new_face_list(mpg%nCg))
      allocate(tmp_face_list(mpg%nCg))
      nf_hlp=0
      old_face_list=0
      new_face_list=0
      tmp_face_list=0
      do iCv = mpg%nCi+1, mpg%nCi+mpg%nCg
        iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
        iCn = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
        if (mpg%cvOnClosedSurface(iCn)) cycle
        nf_hlp = nf_hlp+1
        old_face_list(nf_hlp) = iFc
      end do
cxpb Now order the faces to get a continuous contour
      new_face_list(1) = old_face_list(1)
      old_face_list(1) = 0
      iVa = mpg%fcVx(new_face_list(1),1)
      iVb = mpg%fcVx(new_face_list(1),2)
      inew = 1
      do while (inew.lt.nf_hlp)
        left_match_found = .false.
        do j = 1, nf_hlp
          if (left_match_found) cycle
          if (old_face_list(j).eq.0) cycle
          if (iVa.eq.mpg%fcVx(old_face_list(j),1)) then
            left_match_found = .true.
            iVa = mpg%fcVx(old_face_list(j),2)
          else if (iVa.eq.mpg%fcVx(old_face_list(j),2)) then
            left_match_found = .true.
            iVa = mpg%fcVx(old_face_list(j),1)
          end if
          if (left_match_found) then
            tmp_face_list(1:inew) = new_face_list(1:inew)
            new_face_list(1) = old_face_list(j)
            new_face_list(2:inew+1) = tmp_face_list(1:inew)
            old_face_list(j) = 0
            inew = inew + 1
          end if
        end do
        right_match_found = .false.
        do j = 1, nf_hlp
          if (right_match_found) cycle
          if (old_face_list(j).eq.0) cycle
          if (iVb.eq.mpg%fcVx(old_face_list(j),1)) then
            right_match_found = .true.
            iVb = mpg%fcVx(old_face_list(j),2)
          else if (iVb.eq.mpg%fcVx(old_face_list(j),2)) then
            right_match_found = .true.
            iVb = mpg%fcVx(old_face_list(j),1)
          end if
          if (right_match_found) then
            new_face_list(inew+1) = old_face_list(j)
            old_face_list(j) = 0
            inew = inew + 1
          end if
        end do
        if (.not.(left_match_found.or.right_match_found).and.
     &      inew.lt.nf_hlp) then
          do j = 1, inew
            iFc = new_face_list(j)
            iVa = mpg%fcVx(iFc,1)
            iVb = mpg%fcVx(iFc,2)
            write(6,'(a,4i6)')
     &       'Chain of faces and vertices in new list ',j,iFc,iVa,iVb
          end do
          do j = 1, nf_hlp
            if (old_face_list(j).eq.0) cycle
            iFc = old_face_list(j)
            iVa = mpg%fcVx(iFc,1)
            iVb = mpg%fcVx(iFc,2)
            write(6,'(a,4i6)')
     &       'Leftover faces and vertices in old list',j,iFc,iVa,iVb
          end do
          call xerrab ('Broken face chain along contour')
        end if
      end do
cxpb Now choose the start of the list to be the first face from the
cxpb set with the smallest absolute face label
      flmin = 99
      ifmin = 0
      do i = 1, nf_hlp
        if (abs(mpg%fcLbl(new_face_list(i))).lt.flmin) then
          flmin = abs(mpg%fcLbl(new_face_list(i)))
          ifmin = i
        end if
      end do
      tmp_face_list(1:nf_hlp-ifmin+1) = new_face_list(ifmin:nf_hlp)
      tmp_face_list(nf_hlp-ifmin+2:nf_hlp) = new_face_list(1:ifmin-1)
      new_face_list(1:nf_hlp) = tmp_face_list(1:nf_hlp)
      if (abs(mpg%fcLbl(new_face_list(nf_hlp))).eq.flmin .and.
     &    maxval(abs(mpg%fcLbl(new_face_list(1:nf_hlp)))).gt.flmin) then
        ifmin = nf_hlp
        do while (abs(mpg%fcLbl(new_face_list(ifmin-1))).eq.flmin.and.
     &            ifmin.gt.1)
          ifmin = ifmin-1
        end do
        tmp_face_list(1:nf_hlp-ifmin+1) = new_face_list(ifmin:nf_hlp)
        tmp_face_list(nf_hlp-ifmin+2:nf_hlp) = new_face_list(1:ifmin-1)
        new_face_list(1:nf_hlp) = tmp_face_list(1:nf_hlp)
      end if

c*** Write the header

      open(lun,file='gridedge.dat')

      write (lun,'(2a)') '# Ion fluxes across the outer grid edge.',
     ,                                         '   Version : 27.11.2026'
      write (lun,'(2a)') '# ',trim(label)
      write (lun,*)
      write (lun,'(a,3i6)') '# nCv,nFc,ns = ',mpg%nCv,mpg%nFc,state%ns
      write (hlp_frmr,'(a,i2,a)') '(a,i3,',nspecies,'i3)'
      if(nspecies.gt.0) write (lun,hlp_frmr) '# nspecies,nfluids = ',
     ,                                nspecies,(nfluids(i),i=1,nspecies)

      write (lun,*)
      write (lun,'(2a)') '# il     : cell index in the list;'
      write (lun,'(2a)')
     ,  '# iFc    : face indices in the grid;'
      write (lun,'(2a)')
     ,  '# R , Z  : R and Z coordinates of the ',
     ,                                  'centres of the cell faces [m];'
      write (lun,'(2a)')
     ,  '# x      : reference distance along the grid edge [m];'
      write (lun,'(2a)')
     ,  '# area   : area of the corresponding cell face [m**2];'
      write (lun,'(2a)')
     ,  '# Te, Ti : electron and ion temperatures [eV];'
      write (lun,'(2a)')
     ,  '# flx... : normal fluxes of the ion species [1/m**2/s];'
      write (lun,'(2a)')
     ,  '# fly... : tangential fluxes of the ion species [1/m**2/s];'
      write (lun,*)

      write(hlp_frmr,'(a,i3,a)')
     , '(a17,7x,a1,13x,a1,13x,a1,11x,a4,11x,a2,10x,a2,6x,',
     ,   2*state%ns,'(a5,i2.2,a2,3x))'
c      print *,hlp_frmr               !###
      write (lun,hlp_frmr) '#  il   iFc label',
     ,                     'R','Z','x','area','Ti','Te',
     ,                    (' flx_',mz_hlp(i),hlp_al(i),i=0,state%ns-1),
     ,                    (' fly_',mz_hlp(i),hlp_al(i),i=0,state%ns-1)

      write(hlp_frmr,720) 2*state%ns+2
c      print *,hlp_frmr               !###
  720 format('(1p,i5,i6,1x,i3,2x,4e14.6,',i3,'e12.4)')
      zz=0.0_R8
      do m=1,nf_hlp !{
        hlp_area=geo%fcS(new_face_list(m))
        iCn = min(mpg%fcCv(new_face_list(m),1),
     &            mpg%fcCv(new_face_list(m),2))
        do j=0,state%ns-1 !{
          hlp_flx(j)=state%dv%fna(new_face_list(m),0,j)/hlp_area*
     &     mpg%bcFcOr(mpg%fcBc(new_face_list(m),1))
          hlp_fly(j)=state%dv%fna(new_face_list(m),1,j)/hlp_area*
     &     mpg%bcFcOr(mpg%fcBc(new_face_list(m),1))
        end do !}
        ra=geo%cvX(mpg%fcVx(new_face_list(m),1))
        za=geo%cvY(mpg%fcVx(new_face_list(m),1))
        rb=geo%cvX(mpg%fcVx(new_face_list(m),2))
        zb=geo%cvY(mpg%fcVx(new_face_list(m),2))
        d=geo%fcHt(new_face_list(m))
        re=0.5_R8*(ra+rb)
        ze=0.5_R8*(za+zb)
        z=zz+0.5_R8*d
        zz=zz+d

        write(lun,hlp_frmr)
     ,    m,new_face_list(m),mpg%fcLbl(new_face_list(m)),
     ,    re,ze,z,hlp_area,
     ,    state%pl%ti(iCn)/ev,state%pl%te(iCn)/ev,
     ,   (hlp_flx(j),j=0,state%ns-1),
     ,   (hlp_fly(j),j=0,state%ns-1)
      end do !}

      close (lun)
      deallocate(mz_hlp,hlp_al)
      deallocate(hlp_flx,hlp_fly)
      deallocate(old_face_list)
      deallocate(new_face_list)
      deallocate(tmp_face_list)

      return
c=======================================================================
      end subroutine b2pedges
