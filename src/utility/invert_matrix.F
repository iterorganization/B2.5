
      subroutine invert_matrix (A, Ainv, nnucl, nz2n) 
      use b2mod_types
      implicit none
      integer :: nnucl
      real(kind=R8) :: A(0:2*nnucl-1,0:2*nnucl-1), nz2n(0:nnucl-1),
     &                 Ainv(0:2*nnucl-1,0:2*nnucl-1)
      real(kind=R8), dimension(size(A,1)) :: work  ! work array for LAPACK
      integer, dimension(size(A,1)) :: ipiv   ! pivot indices
      integer :: info, dnnucl, i, j
      external DGETRF
      external DGETRI

*    ..the matrix inversion
      dnnucl=2*nnucl
      Ainv=A
      ! DGETRF computes an LU factorization of a general M-by-N matrix A
      ! using partial pivoting with row interchanges.
      call DGETRF(dnnucl, dnnucl, Ainv, dnnucl, ipiv, info)

      if (info /= 0) then
         print *, nz2n
         print *
         do i = 0, dnnucl-1
           print *, (A(i,j), j = 0, dnnucl-1 )
         enddo
         print *
#ifndef NAGFOR
         call sleep(3)
#endif
         flush(6)
         call xerrab ('Matrix is numerically singular!')
      end if

      ! DGETRI computes the inverse of a matrix using the LU factorization
      ! computed by DGETRF.
      call DGETRI(dnnucl, Ainv, dnnucl, ipiv, work, dnnucl, info)

      call xertst(info.eq.0, 'Matrix inversion failed!')

      return
      end subroutine invert_matrix


