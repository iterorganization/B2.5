      subroutine cv_intersections (geo, m, segm, list, list_ci,
     &   list_fc, nncv, nnci, nnfc, nmax, mode, ifail)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nmax
      type (geometry), intent (in) :: geo
      type(mapping), intent(in) :: m
      real(kind=r8) :: segm(2,2)
      integer :: list(nmax), list_ci(m%nCi), nncv, nnci, mode
      integer :: list_fc(m%nFc), nnfc
      integer :: ic, ifc, ii, kk
      real(kind=r8) :: p1(2), q1(2), p2(2), q2(2)
      logical found, intersects, ifail
      external intersects, xertst

***
*     This routine finds all the cv intersecting with the segment specified in segm
*     It is assumed that segm(1,:) contains x-coordinate and segm(2,:) the y-coordinates
*     list will contain all the CVs intersecting (without repetitions) and nncv the number of 
*     CVs. list_ci will only contain the internal CVs intersecting, and nnci the number of CVs.
*     The CVs will be sorted according to mode, mode=1 sorted along x-dir, mode=2 sorted along y-dir
***

      call xertst(mode.ge.1.and.mode.le.2,'mode value not encoded')
      nncv = 1 !number of CVs in the list
      nnci = 0
      nnfc = 0
      list = 0
      call ifill (m%nCi, 0, list_ci, 1)
      call ifill (m%nFc, 0, list_fc, 1)
      ifail = .false.
      p1=segm(:,1)
      q1=segm(:,2)
      do ifc = 1, m%nFc
        p2(1) = geo%vxX(m%fcVx(ifc,1))
        p2(2) = geo%vxY(m%fcVx(ifc,1))
        q2(1) = geo%vxX(m%fcVx(ifc,2))
        q2(2) = geo%vxY(m%fcVx(ifc,2))
        if (intersects(p1,q1,p2,q2)) then
          nnfc = nnfc + 1
          list_fc(nnfc) = ifc
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
            if (geo%cvX(list(kk)).le.geo%cvX(ic)) exit
          else
            if (geo%cvY(list(kk)).le.geo%cvY(ic)) exit
          endif
            list(kk+1)=list(kk)
        enddo
        list(kk+1)=ic
      enddo
*     Now create the nCi list
      kk = 1
      do ic = 1,nncv
        if(list(ic).le.m%nCi) then
          list_ci(kk) = list(ic)
          kk=kk+1
        endif
      enddo
      nnci = kk-1

*     Now sort the faces list
c      nnfc = nnfc-1
      do ifc=2,nnfc
        ic = list_fc(ifc)
        do kk=ifc-1,1,-1
c          if (mode.eq.1) then
c            if (geo%cvX(list_fc(kk)).le.geo%cvX(ic)) exit
c          else
c            if (geo%cvY(list_fc(kk)).le.geo%cvY(ic)) exit
c          endif
            list_fc(kk+1)=list_fc(kk)
        enddo
        list_fc(kk+1)=ic
      enddo

      return
      end subroutine cv_intersections
