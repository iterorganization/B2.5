

!> Computes additional connectivity data based
!! on the connectivity data of the traduitoutb2us file
      subroutine calc_add_connectivity(g,m)
      use b2mod_types
      !use b2mod_indirect   ! causes issue
      !use b2mod_geo        ! causes issue
      use b2us_map
      use b2us_geo
      implicit none

      !!.. arguments
      type(geometry) :: g
      type(mapping) :: m

      !!.. local variables
      integer :: i, ic, ivx, ifc, n, vxFc_end, vxCv_end
      integer, allocatable :: indFc(:), indCv(:), fh(:), ch(:), fch(:), &
                        & indcvfc(:)
      intrinsic pack, count

      !!Read out fcVx, cvFc, cvFcP, cvVx and cvVxP
      !!Construct vxCv, vxCvP, vxFcP, fcCv !(because of ghostcell, every face has two cells)
  

      !Construct vxFcP
      allocate (indFc(2*m%nFc))
      indFc(1:m%nFc)          = (/ (i, i=1,m%nFc) /)
      indFc(m%nFc+1:2*m%nFc)    = (/ (i, i=1,m%nFc) /)

      allocate (fh(2*m%nFc))
      fh(1:m%nFc)             = m%fcVx(1:m%nFc,1)
      fh(m%nFc+1:2*m%nFc)       = m%fcVx(1:m%nFc,2)

      allocate (indCv(m%nCmxVx))
      do ic = 1,m%nCv
         indCv(m%cvVxP(ic,1):m%cvVxP(ic,1)+m%cvVxP(ic,2)-1) = ic
      enddo

      allocate (ch(m%nCmxVx))
      ch(1:m%nCmxVx) = m%cvVx(1:m%nCmxVx)

      m%vxFcP(1,1) = 1
      m%vxCvP(1,1) = 1
      do ivx = 1, m%nVx-1
         m%vxFcP(ivx,2) = count (fh == ivx)
         m%vxFcP(ivx+1,1) = m%vxFcP(ivx,1) + m%vxFcP(ivx,2)
         if (m%vxFcP(ivx,2) > 0) then
             m%vxFc(m%vxFcP(ivx,1):m%vxFcP(ivx+1,1)-1) = &
           &  pack(indFc,fh == ivx)
         endif

         m%vxCvP(ivx,2) = count(ch == ivx)
         m%vxCvP(ivx+1,1) = m%vxCvP(ivx,1) + m%vxCvP(ivx,2)
         if (m%vxCvP(ivx,2) > 0 ) then
            m%vxCv(m%vxCvP(ivx,1):m%vxCvP(ivx+1,1)-1) = &
           &   pack( indCv, ch == ivx ) 
         endif
      enddo

      ivx = m%nVx
      m%vxFcP(ivx,2) = count(fh == ivx)
      vxFc_end = m%vxFcP(ivx,1) + m%vxFcP(ivx,2) -1
      if (m%vxFcP(ivx,2) > 0) then
          m%vxFc(m%vxFcP(ivx,1):vxFc_end) = &
         &  pack( indFc, fh == ivx )
      end if
      m%vxCvP(ivx,2) = count(ch == ivx)
      vxCv_end = m%vxCvP(ivx,1) + m%vxCvP(ivx,2) - 1
      if (m%vxCvP(ivx,2) > 0) then
         m%vxCv(m%vxCvP(ivx,1):vxCv_end) = &
        &  pack( indCv, ch == ivx )
      endif

      !Try analog for fcCv(nFc,2), similar to vxCv - check this in debugger

      allocate (indcvfc(m%nCmxFc))
      do ic = 1,m%nCv
         indcvfc(m%cvFcP(ic,1):m%cvFcP(ic,1)+m%cvFcP(ic,2)-1) = ic
      enddo

      allocate (fch(m%nCmxFc))
      fch(1:m%nCmxFc) = m%cvFc(1:m%nCmxFc)

      !hier loop over de faces

      do ifc = 1,m%nFc
           n = count(fch == ifc)
           if (n .eq. 2) then
             m%fcCv(ifc,1:2) = pack(indcvfc, fch == ifc)
           else
             call xerrab &
             & ('calc_add_connectivity: Wrong # of cell per face')
           endif
      enddo


      deallocate (indFc)
      deallocate (indCv)
      deallocate (indcvfc)
      deallocate (fh)
      deallocate (ch)
      deallocate (fch)


      end subroutine calc_add_connectivity

! For the ghost cell, get the cell index
! of the boundary cell to which the ghost cell is connected
      integer function ghostGetBndCell(iCv,g,m)
      use b2mod_connectivity
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer :: iCv, face_g, cells_g(0:1)
      type(geometry_ag) :: g
      type(mapping_ag) :: m

       if (g%cvCflags(iCv,1).eq.GRID_GUARD) then
          !steps: 
          !- ghostcell has one face => get that face via cvFc
          !- get the two cell of that face
          ! cell which is not equal to the ghostcell itself is the needed boundary cell

          face_g = m%cvFc(m%cvFcP(iCv,1))
          cells_g = m%fcCv(face_g,:)

          if ((cells_g(0).eq.iCv) .and. &
     &         (g%cvCflags(cells_g(1),1).eq.GRID_BOUNDARY)) then
              ghostGetBndCell = cells_g(1)
          elseif  ((cells_g(1).eq.iCv) .and. &
     &         (g%cvCflags(cells_g(0),1).eq.GRID_BOUNDARY)) then
              ghostGetBndCell = cells_g(0)
          else
              call xerrab  &
     &         ('ghostGetBndCell: ghostGetBndCell not allocated')
          endif

       else
          call xerrab('ghostGetBndCell: cell is no ghost cell')
       endif



      end function ghostGetBndCell
