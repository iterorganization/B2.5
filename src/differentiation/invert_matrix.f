
      subroutine invert_matrix (A, Ainv, nnucl, nz2n) 
      use b2mod_types
      implicit none
      integer :: nnucl
      real(kind=R8) :: A(0:2*nnucl-1,0:2*nnucl-1), nz2n(0:nnucl-1),
     &                 Ainv(0:2*nnucl-1,0:2*nnucl-1)
      integer :: dnnucl, i, j

      dnnucl = 2*nnucl
      print *, nz2n
      do i = 0, dnnucl-1
        do j = 0, dnnucl-1
          Ainv(i,j) = 1.0_R8/A(i,j)
        end do
      end do

      return
      end subroutine invert_matrix


