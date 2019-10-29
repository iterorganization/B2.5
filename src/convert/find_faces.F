#ifdef CFLABEL
      subroutine find_faces(indss, inbc, nx, ny, ns,
     .                      bcl, ndfc, Fc, FcOr,
     .                      m, idb)
#else
      subroutine find_faces(indss, inbc, bcl, ndfc, Fc, FcOr,
     .                      m, idb)
#endif
      
      use b2mod_types
      use b2mod_indirect
      use b2us_map
      use b2us_geo
      use carre_constants
      implicit none
#ifdef CFLABEL
      integer, intent(in) :: indss, nx, ny, ns, ndfc, idb
#else
      integer, intent(in) :: indss, ndfc, idb
#endif
      integer, intent(inout) :: inbc, bcl
      integer, intent(inout) :: Fc(ndfc)
      real(R8), intent(inout) :: FcOr(ndFc)
      
      type (mapping), intent (inout) :: m

      integer :: ncoss, icoss,
     .           ix, iy, iFc1, i, imn, imx, iFc, 
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

      
!  find all faces belonging to surface structure INDSS
      
#ifdef CFLABEL
      ncoss = count(cflags(:,:,2:5) == indss)
      if (ncoss.eq.0) return
      allocate (fcss(ncoss))
      if (.not.allocated(lfchit)) allocate (lfchit(m%nFc))
      lfchit = .false.
      icoss = 0
      do iy = -1, ny
        do ix = 0, nx
          if (count(cflags(ix,iy,2:5) == indss) > 0) then
! left face 
            if (cflags(ix,iy,CELLFLAG_LEFTFACE) == indss) then
              iFc = m%imapFcx(ix,iy)
              if (iFc == 0) then
                write (*,*)
                write (*,*) 'iFc = 0 in '
                write (*,*) 'ix, iy ',ix,iy
                write (*,*) 'cflags ',cflags(ix,iy,:)
                cycle
              end if
              if (.not.lfchit(iFc)) then
                icoss = icoss + 1
                fcss(icoss) = iFc
                lfchit(ifc) =.true.
              end if
            end if             
! right face
            if (cflags(ix,iy,CELLFLAG_RIGHTFACE) == indss) then
              if (rightix(ix,iy).ne.nx+1) then
                iFc = m%imapFcx(rightix(ix,iy),rightiy(ix,iy))
                if (iFc == 0) then
                  write (*,*)
                  write (*,*) 'iFc = 0 in '
                  write (*,*) 'ix, iy ',ix,iy
                  write (*,*) 'right ix,iy ',
     .                 rightix(ix,iy),rightiy(ix,iy)
                  write (*,*) 'cflags ',cflags(ix,iy,:)
                  cycle
                end if
                if (.not.lfchit(iFc)) then
                  icoss = icoss + 1
                  fcss(icoss) = iFc
                  lfchit(ifc) =.true.
                end if
              endif
            end if
! bottom face                   
            if (cflags(ix,iy,CELLFLAG_BOTTOMFACE) == indss) then
              iFc = m%imapFcy(ix,iy)
              if (iFc == 0) then
                write (*,*)
                write (*,*) 'iFc = 0 in '
                write (*,*) 'ix, iy ',ix,iy
                write (*,*) 'cflags ',cflags(ix,iy,:)
                cycle
              end if
              if (.not.lfchit(iFc)) then
                icoss = icoss + 1
                fcss(icoss) = iFc
                lfchit(ifc) =.true.
              end if
            end if
! top face                  
            if (cflags(ix,iy,CELLFLAG_TOPFACE) == indss) then
              if (topiy(ix,iy).ne.ny+1) then
                iFc = m%imapFcy(topix(ix,iy),topiy(ix,iy))
                if (iFc == 0) then
                  write (*,*)
                  write (*,*) 'iFc = 0 in '
                  write (*,*) 'ix, iy ',ix,iy
                  write (*,*) 'top ix,iy ',
     .                 topix(ix,iy),topiy(ix,iy)
                  write (*,*) 'cflags ',cflags(ix,iy,:)
                  cycle
                end if
                if (.not.lfchit(iFc)) then
                  icoss = icoss + 1
                  fcss(icoss) = iFc
                  lfchit(ifc) =.true.
                end if
              endif
            end if             
          end if
        end do                  ! ix
      end do                    ! iy
     
      do i = 1, icoss
        write (*,'(a,4i6,4es16.7)')' i,fcss(i),vx1,vx2 ',
     .             i,fcss(i),m%fcVx(fcss(i),1), m%fcVx(fcss(i),2)
!     .           , g%vxX(m%fcVx(fcss(i),1)),g%vxY(m%fcVx(fcss(i),1)),
!     .             g%vxX(m%fcVx(fcss(i),2)),g%vxY(m%fcVx(fcss(i),2))
      end do
      ncoss = icoss
#else
      ncoss = count(m%fcLbl == indss)
      if (ncoss.eq.0) return
      allocate (fcss(ncoss))
      allocate (indfc(m%nFc))
      indfc = [(i,i=1,m%nFc)]
      fcss(1:ncoss) = pack(indfc, m%fcLbl == indss)
      deallocate (indfc)
#endif
       
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
