      subroutine cv_intersections (nCv, geo, m, segm, list, nncv, 
     &  nmax, mode, ifail)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nmax
      type (geometry), intent (in) :: geo
      type(mapping), intent(in) :: m
      real(kind=r8) :: segm(2,2)
      integer :: list(nmax), nncv, mode
      integer :: ic, ifc, ii, kk
      real(kind=r8) :: p1(2), q1(2), p2(2), q2(2), stemp
      logical found, intersects, ifail
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
        if (intersects(p1,q1,p2,q2)) then
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
*     Now sort the obtained list
      do ii=2,nncv
        ic = list(ii)
        do kk=ii-1,1,-1
          if (mode.eq.1) then
            if (geo%cVx(list(kk)).le.geo%cVx(ic)) exit
          else
            if (geo%cVy(list(kk)).le.geo%cVy(ic)) exit
          endif
            list(kk+1)=list(kk)
        enddo
        list(kk+1)=ic
      enddo

      return
      end subroutine cv_intersections
