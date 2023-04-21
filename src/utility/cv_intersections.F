      recursive subroutine cv_intersections (nCv, geo, m, segm, 
     &  list, nncv, nmax, mode, ifail)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nmax
      type (geometry), intent (in) :: geo
      type(mapping), intent(in) :: m
      real(kind=r8) :: segm(2,2)
      integer :: list(nmax), nncv, mode, list_sorted(nmax)
      integer :: ic, ifc, ivx, ii, kk, i, j, k, l, icn
      integer :: icsep, ic1, ic2, ifc1, ifc2, ifc11, ifc22, ifcc
      real(kind=r8) :: p1(2), q1(2), p2(2), q2(2), stemp
      logical found, intersects, ifail, in_list, in_list_sorted
      external intersects, xertst


***
*     This routine finds all the cv intersecting with the segment specified in segm
*     It is assumed that segm(1,:) contains x-coordinate and segm(2,:) the y-coordinates
*     list will contain all the CVs intersecting (without repetitions) and nncv the number of 
*     CVs. The CVs will be sorted according to connectivity.
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
*     Now sort the obtained list based on connectivity
*     First find a boundary cell to start
       list_sorted = 0     
       do i = 1,nncv
         ic = list(i)
         if (ic.gt.m%nCi) then
            if (list_sorted(1).eq.0) then
              list_sorted(1) = ic
            elseif (mode.eq.1 .and.
     &      geo%cvX(ic).lt.geo%cvX(list_sorted(1)) ) then
                list_sorted(1) = ic   
            elseif (mode.eq.2 .and. 
     &      geo%cvY(ic).lt.geo%cvY(list_sorted(1)) ) then
                list_sorted(1) = ic 
              endif              
           endif
       enddo

*    Start sorting itself
       do i = 1,nncv-1
         ic = list_sorted(i)
         ! First through poloidal faces         
         do j = 1,m%cvFcP(ic,2)
           ifc = m%cvFc(m%cvFcP(ic,1)+j-1)
           if (geo%fcQalf(ifc,0).ne.0) then
             do l = 1,2
               in_list = .false.
               in_list_sorted = .false.
               icn = m%fcCv(ifc,l)
               if (icn.ne.ic .and. list_sorted(i+1).eq.0) then
                 ! Check whether icn is in list but not in list_sorted
                 do k = 1,nncv
                   if (list(k).eq.icn) in_list = .true.
                 enddo
                 do k = 1,i
                   if (list_sorted(k).eq.icn) in_list_sorted = .true.
                 enddo
                 if (in_list .and. .not.in_list_sorted ) then
                   list_sorted(i+1) = icn
                 endif
               endif
             enddo
           endif
         enddo

         ! Second through radial faces
         if (list_sorted(i+1).eq.0) then
         do j = 1,m%cvFcP(ic,2)
           ifc = m%cvFc(m%cvFcP(ic,1)+j-1)
           if (geo%fcQalf(ifc,0).eq.0) then
             do l = 1,2
               in_list = .false.
               in_list_sorted = .false.
               icn = m%fcCv(ifc,l)
               if (icn.ne.ic .and. list_sorted(i+1).eq.0) then
                 ! Check whether icn is in list but not in list_sorted
                 do k = 1,nncv
                   if (list(k).eq.icn) in_list = .true.
                 enddo
                 do k = 1,i
                   if (list_sorted(k).eq.icn) in_list_sorted = .true.
                 enddo
                 if (in_list .and. .not.in_list_sorted ) then
                   list_sorted(i+1) = icn
                 endif
               endif
             enddo
           endif
         enddo
         endif

         ! Third try to sort via vertex connectivity
         if (list_sorted(i+1).eq.0) then
           do j = 1,m%cvVxP(ic,2)
             ivx = m%cvVx(m%cvVxP(ic,1)+j-1)
             do l = 1,m%vxCvP(ivx,2)
               in_list = .false.
               in_list_sorted = .false.
               icn = m%vxCv(m%vxCvP(ivx,1)+l-1)
               if ( icn.ne.ic .and. list_sorted(i+1).eq.0 ) then
                 ! Check whether icn is in list but not in list_sorted
                 do k = 1,nncv
                   if (list(k).eq.icn) in_list = .true.
                 enddo
                 do k = 1,i
                   if (list_sorted(k).eq.icn) in_list_sorted = .true.
                 enddo
                 if ( in_list .and. .not.in_list_sorted ) then
                   list_sorted(i+1) = icn
                 endif                
               endif
             enddo
           enddo
         endif
       enddo
       ! Save
       list = list_sorted

       ! Check whether separatrix can be identified
       icsep = 0
       ic = 1
       found = .false.
       do while (ic.le.nncv .and..not.found)
         if (m%cvReg(list(ic)).ne.m%cvReg(list(ic+1))) then
           found = .true.
           icsep = ic+1
         endif
       ic = ic + 1 
       end do
       if (.not.found) write(*,*) ('Did not find separatrix location. '
     &  // 'Check specification OMP and IMP in b2.user.paramters.')
       ! Try to get the common face
       ic1 = list(icsep-1) !cv of flux tube just inside sep
       ic2 = list(icsep) !cv of flux tube just outside sep
       ifc1 = m%cvFcP(ic1,1)
       ifc2 = ifc1 + m%cvFcP(ic1,2) - 1
       ifc11 = m%cvFcP(ic2,1)
       ifc22 = ifc11 + m%cvFcP(ic2,2) - 1
       ! find face belonging to both ic1 and ic2
       found = .false.
       do i = ifc1, ifc2
         ifc = m%cvFc(i)
         do j = ifc11, ifc22
           ifcc = m%cvFc(j)
           if (ifc.eq.ifcc) then
             found = .true.
           endif
         enddo
       end do  
       
       if (.not.found) then
       ! Call cv_intersects again with slightly shifted segm
         segm = segm + 1.0E-6
         call cv_intersections (nCv, geo, m, segm, 
     &  list, nncv, nmax, mode, ifail)
       endif     


      return
      end subroutine cv_intersections
