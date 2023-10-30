      subroutine write_triangles (nnodes, nelem, xco, yco, elem,
     .    neighbor, bndmdl, elem_cv, trifc)

      use b2mod_types
      use b2us_map
      use b2us_geo

      implicit none

      integer, intent(in) :: nnodes, nelem
      integer, intent(in) :: elem(3,nelem), neighbor(3,2,nelem),
     .    bndmdl(3,nelem), elem_cv(nelem), trifc(3,nelem)
      real(R8), intent(in) :: xco(nnodes), yco(nnodes)

      ! write triangle.coordinates
      call write_npco (nnodes, xco, yco, 1.0_R8)

      ! write triangle.elements
      call write_elemente (nelem, elem)

      ! write triangle.neighbors
      call write_neighbors (nelem, neighbor, bndmdl, elem_cv, trifc)

      return
      end subroutine write_triangles


      subroutine write_npco (nco, xco, yco, sclfac)

      use b2mod_types

      implicit none
      integer, intent(in) :: nco
      real(R8), intent(in) :: xco(nco), yco(nco), sclfac
      integer :: iun, i
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.coordinates',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.coordinates',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nco
      do i = 1, nco
        write (iun,'(i9,2es24.16)') i, xco(i)*sclfac,
     .                                 yco(i)*sclfac
      end do
      close (iun)

      return
      end subroutine write_npco


      subroutine write_elemente (nel, elem)

      use b2mod_types

      implicit none
      integer, intent(in) :: nel, elem(3,nel)
      integer :: iun, i, is
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.elements',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.elements',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nel
      do i = 1, nel
        write (iun,'(4i9)') i, (elem(is,i), is = 1,3)
      end do
      close (iun)

      return
      end subroutine write_elemente


      subroutine write_neighbors (nel, neighbor, bndmdl, elem_cv, triFc)

      use b2mod_types

      implicit none
      integer, intent(in) :: nel, neighbor(3,2,nel),bndmdl(3,nel),
     .                       elem_cv(nel), triFc(3,nel)
      integer :: iun, i, is, js
#ifndef F2003
      integer newunit
      external newunit
#endif

#ifdef F2003
      open (newunit=iun,file='triangle.neighbors',
     .           access='SEQUENTIAL',form='FORMATTED')
#else
      iun = newunit()
      open (unit=iun,file='triangle.neighbors',
     .           access='SEQUENTIAL',form='FORMATTED')
#endif
      write (iun,'(i9)') nel
      do i = 1, nel
        write (iun,'(i9,3x,3(3i9,3x),4i9)')
     .     i, (neighbor(is,1,i), neighbor(is,2,i),
     .         bndmdl(is,i), is = 1,3), elem_cv(i),
     .        (triFc(js,i), js = 1,3)
      end do
      close (iun)

      return
      end subroutine write_neighbors
