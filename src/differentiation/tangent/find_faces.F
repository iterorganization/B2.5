      subroutine find_faces(indss, inbc, bcl, ndfc, Fc, FcOr,
     .                      m, idb)
      
      use b2mod_types
      use b2us_map_diffv
      use b2us_geo_diffv
      use carre_constants
      implicit none
      integer, intent(in) :: indss, ndfc, idb
      integer, intent(inout) :: inbc, bcl
      integer, intent(inout) :: Fc(ndfc)
      real(R8), intent(inout) :: FcOr(ndFc)
      
      type (mapping), intent (inout) :: m

      integer :: ncoss, icoss,
     .           iFc1, i, imn, imx,
     .           ncout, istart
      integer, allocatable :: fcss(:), indfc(:)
      logical, allocatable :: lout(:)
      logical, allocatable :: lfchit(:)
      
      type fcpart
        integer :: fcno
        integer :: vx1, vx2
        type(fcpart), pointer :: prev, next
      end type fcpart

      type(fcpart), pointer :: head, tail
      type(fcpart), pointer :: face
      intrinsic pack
      
!  find all faces belonging to surface structure INDSS
      ncoss = count(m%fcLbl == indss)
      if (ncoss.eq.0) return
      allocate (fcss(ncoss))
      allocate (indfc(m%nFc))
      indfc = [(i,i=1,m%nFc)]
      fcss(1:ncoss) = pack(indfc, m%fcLbl == indss)
      deallocate (indfc)
       
      bcl = bcl + ncoss

! sort face for a continuous surface
      istart = 1
      allocate (lout(ncoss))
      lout = .false.
      ncout = ncoss

      do while (count(lout) < ncoss)
        allocate(face)
        nullify(face%prev)
        nullify(face%next)
        face%fcno = fcss(istart)
        lout(istart) =.true.
        face%vx1 = m%fcVx(face%fcno,1)
        face%vx2 = m%fcVx(face%fcno,2)
        head => face
        tail => face

        do while (ncout > 0)
          do icoss = 2, ncoss
            if (lout(icoss)) cycle
            imn = min(m%fcVX(fcss(icoss),1),m%fcVX(fcss(icoss),2))
            imx = max(m%fcVX(fcss(icoss),1),m%fcVX(fcss(icoss),2))
            if ((head%vx1 == imn) .or. (head%vx1 == imx)) then
              allocate(face) 
              nullify(face%prev)
              nullify(face%next)
              face%fcno = fcss(icoss)
              face%vx1 = imn+imx-head%vx1
              face%vx2 = head%vx1
              face%next => head
              head%prev => face
              head => face
              lout(icoss) = .true.
              ncout = ncout - 1
              exit
            end if
            if ((tail%vx2 == imn) .or. (tail%vx2 == imx)) then
              allocate(face) 
              nullify(face%prev)
              nullify(face%next)
              face%fcno = fcss(icoss)
              face%vx1 = tail%vx2
              face%vx2 = imn+imx-tail%vx2
              face%prev => tail
              tail%next => face
              tail => face
              lout(icoss) = .true.
              ncout = ncout - 1
              exit
            end if
          end do                ! icoss
          if (icoss > ncoss) then
            write (*,*) 'face list, indss = ',indss
            write (idb,*) 'face list, indss = ',indss
            do while (associated(head))
              write (*,*) head%fcno, head%vx1, head%vx2
!              write (idb,'(i6,2es16.7)')
!     .             head%vx1,g%vxX(head%vx1),g%vxY(head%vx1)
!              write (idb,'(i6,2es16.7)')
!     .             head%vx2,g%vxX(head%vx2),g%vxY(head%vx2)
              inbc = inbc + 1
              iFc1 = head%fcno
              Fc(inbc) = iFc1
              if (m%fcCv(iFc1,1) .lt.  m%fcCv(iFc1,2))
     >             FcOr(inbc) = 1.0_R8
              if (m%fcCv(iFc1,1) .gt.  m%fcCv(iFc1,2))
     >             FcOr(inbc) = -1.0_R8
              tail => head
              head => head%next
              deallocate (tail)
            end do

            do i = 1, ncoss
              if (.not.lout(i)) then
                istart = i
                exit
              end if
            end do
            exit
          end if
        end do
      end do              
      deallocate(lout)
!     
!     store faces in m%bcFc(inbc)...
      do while (associated(head))
        inbc = inbc + 1
        iFc1 = head%fcno
        Fc(inbc) = iFc1
        if (m%fcCv(iFc1,1) .lt.  m%fcCv(iFc1,2))
     >       FcOr(inbc) = 1.0_R8
        if (m%fcCv(iFc1,1) .gt.  m%fcCv(iFc1,2))
     >       FcOr(inbc) = -1.0_R8
        tail => head
        head => head%next
        deallocate (tail)
      end do

      deallocate (fcss)
      if (allocated(lfchit)) deallocate (lfchit)
      
      end subroutine find_faces
