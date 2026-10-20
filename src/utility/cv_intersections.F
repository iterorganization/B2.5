      subroutine cv_intersections (geo, m, segm, list, nncv,
     &  nmax, mode, ifail)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nmax
      type (geometry), intent (in) :: geo
      type(mapping), intent(in) :: m
      real(kind=r8) :: segm(2,2)
      integer :: list(nmax), nncv, mode
      integer :: ic, ifc, ii, kk
      real(kind=r8) :: p1(2), q1(2), p2(2), q2(2)
      real(kind=r8) :: sgn, key
      real(kind=r8), allocatable :: pos(:)
      logical found, intersects, ifail, slab_test
      external intersects, xertst

***
*     This routine finds all the cv intersecting with the segment specified in segm
*     It is assumed that segm(1,:) contains x-coordinate and segm(2,:) the y-coordinates
*     list will contain all the CVs intersecting (without repetitions) and nncv the number of
*     CVs. The CVs will be sorted according to mode, mode=1 sorted along x-dir, mode=2 sorted along y-dir
***

      call xertst(mode.ge.1.and.mode.le.2,'mode value not encoded')
      nncv = 1 !number of CVs in the list
      list = 0
      ifail = .false.
      p1=segm(:,1)
      q1=segm(:,2)
      do ifc = 1, m%nFc
        p2(1) = geo%vxX(m%fcVx(ifc,1))
        p2(2) = geo%vxY(m%fcVx(ifc,1))
        q2(1) = geo%vxX(m%fcVx(ifc,2))
        q2(2) = geo%vxY(m%fcVx(ifc,2))
!! for artificial slab cases it sometimes happens that vxY and vxX are wrong at the cuts, leading to incorrect intersections.
!! fcHt will be correct though. Therefore we check here if the distance between the vertex points corresponds to the face length. If not, no intersection is counted.
        slab_test = .true.                                      
        if (sqrt((p2(1)-q2(1))**2+(p2(2)-q2(2))**2).gt.         
     &                              geo%fcHt(ifc)*1.1_R8) then 
          slab_test = .false.
        endif
        if (intersects(p1,q1,p2,q2).and.slab_test) then
          do ii=1,2
            ic=m%fcCv(ifc,ii)
            found = .false.
            kk = 1
            do while (.not.found .and. kk.le.nncv)
              if (ic.ne.0.and.list(kk).eq.ic) found = .true.
              kk = kk + 1
            enddo
            if (.not.found) then
              if (nncv.gt.nmax) then
                ifail = .true.
                return
              endif
              list(nncv) = ic
              nncv = nncv + 1
            endif
          enddo
        endif
      enddo
      nncv = nncv-1
*     Now sort the obtained list along the segment: the sort key of a cell is
*     the position of its chord midpoint on the segment (cv_segment_positions)..
      if (nncv.gt.1) then
        allocate(pos(nncv))
        call cv_segment_positions (geo, m, segm, list, nncv, pos)
        sgn = 1.0_R8
        if (segm(mode,2).lt.segm(mode,1)) sgn = -1.0_R8
        do ii=2,nncv
          ic = list(ii)
          key = sgn*pos(ii)
          do kk=ii-1,1,-1
            if (sgn*pos(kk).le.key) exit
            list(kk+1) = list(kk)
            pos(kk+1) = pos(kk)
          enddo
          list(kk+1) = ic
          pos(kk+1) = key*sgn
        enddo
        deallocate(pos)
      endif

      return
      end subroutine cv_intersections


      subroutine cv_segment_positions (geo, m, segm, list, nncv, pos)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      type (geometry), intent (in) :: geo
      type(mapping), intent(in) :: m
      real(kind=r8), intent(in) :: segm(2,2)
      integer, intent(in) :: nncv
      integer, intent(in) :: list(nncv)
      real(kind=r8), intent(out) :: pos(nncv)
      integer :: ii, ic, i, ifc, nhit
      real(kind=r8) :: d(2), e(2), w(2), dlen, den, t, u, tmin, tmax
      real(kind=r8), parameter :: tol = 1.0e-9_R8

      d = segm(:,2) - segm(:,1)
      dlen = sqrt(d(1)**2 + d(2)**2)
      do ii = 1, nncv
        ic = list(ii)
        nhit = 0
        tmin = 0.0_R8
        tmax = 0.0_R8
        do i = m%cvFcP(ic,1), m%cvFcP(ic,1) + m%cvFcP(ic,2) - 1
          ifc = m%cvFc(i)
          e(1) = geo%vxX(m%fcVx(ifc,2)) - geo%vxX(m%fcVx(ifc,1))
          e(2) = geo%vxY(m%fcVx(ifc,2)) - geo%vxY(m%fcVx(ifc,1))
          den = d(1)*e(2) - d(2)*e(1)
          if (abs(den).le.tol*dlen*sqrt(e(1)**2+e(2)**2)) cycle
          w(1) = geo%vxX(m%fcVx(ifc,1)) - segm(1,1)
          w(2) = geo%vxY(m%fcVx(ifc,1)) - segm(2,1)
          t = (w(1)*e(2) - w(2)*e(1))/den
          u = (w(1)*d(2) - w(2)*d(1))/den
          if (t.lt.-tol .or. t.gt.1.0_R8+tol) cycle
          if (u.lt.-tol .or. u.gt.1.0_R8+tol) cycle
          if (nhit.eq.0) then
            tmin = t
            tmax = t
          else
            tmin = min(tmin, t)
            tmax = max(tmax, t)
          endif
          nhit = nhit + 1
        enddo
        if (nhit.gt.0) then
          pos(ii) = 0.5_R8*(tmin + tmax)*dlen
        else
          pos(ii) = ((geo%cvX(ic) - segm(1,1))*d(1) +
     &               (geo%cvY(ic) - segm(2,1))*d(2))/dlen
        endif
      enddo

      return
      end subroutine cv_segment_positions
