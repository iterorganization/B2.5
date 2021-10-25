      subroutine sort_faces(nFcList,fcList,nFcLbl,fcLbl,iLbl,m)
      !! Routine that looks for all faces in fcLbl that have a particular
      !! face label/index iLbl, and returns them as a sorted list fcList
      !!
      !! Implicit assumption: faces connect into a single line
      !! nFcList         : number of faces with label iLbl
      !! fcList(nFcList) : sorted list of faces (on output)
      !! nFcLbl          : number of faces in fcLbl
      !! fcLbl(nFcLbl)   : face labels assigned to all faces
      !! iLbl            : label to look for within fcLbl
      !! m               : mapping
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none
      integer, intent (in)        :: nFcList, nFcLbl, iLbl
      integer, intent (inout)     :: fcList(nFcList), fcLbl(nFcLbl)
      type (mapping), intent (in) :: m

      integer :: ncoss, icoss,
     .           iFc, iFcList, i, imn, imx,
     .           ncout, istart
      integer, allocatable :: fcss(:), indfc(:)
      logical, allocatable :: lout(:)
      
      type fcpart
        integer :: fcno
	integer :: vx1, vx2
	type(fcpart), pointer :: prev, next
      end type fcpart

      type(fcpart), pointer :: head, tail
      type(fcpart), pointer :: face
      intrinsic pack
      
!     find all faces in fcLbl with label iLbl
      ncoss = count(fcLbl == iLbl)
      if (ncoss.ne.nFcList) call xerrab ('Inconsistent input')
      if (ncoss.eq.0) return

      allocate (fcss(ncoss))
      allocate (indfc(nFcLbl))
      indfc = [(i,i=1,nFcLbl)]
      fcss(1:ncoss) = pack(indfc, fcLbl == iLbl)
      deallocate (indfc)
       
!     sort faces for a continuous surface
      iFcList = 0
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
            imn = min(m%fcVx(fcss(icoss),1),m%fcVx(fcss(icoss),2))
            imx = max(m%fcVx(fcss(icoss),1),m%fcVx(fcss(icoss),2))
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
!            write (*,*) 'face list, iLbl = ',iLbl
            do while (associated(head))
!              write (*,*) head%fcno, head%vx1, head%vx2
              iFcList = iFcList + 1
              iFc = head%fcno
              fcList(iFcList) = iFc
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

      do while (associated(head))
        iFcList = iFcList + 1
        iFc = head%fcno
        fcList(iFcList) = iFc
        tail => head
        head => head%next
        deallocate (tail)
      end do

      deallocate (fcss)
      
      end subroutine sort_faces
