      subroutine invert_matrix_dv (A, Ad, Ainv, Ainvd, nnucl, nz2n,
     & nbdirs)
      use b2mod_types
      use b2mod_diffsizes
      implicit none
      integer :: nnucl, nbdirs
      real(kind=R8) :: A(0:2*nnucl-1,0:2*nnucl-1), nz2n(0:nnucl-1),
     &                 Ad(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1)
      real (kind=R8) :: Ainvd(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
      real (kind=R8) :: dum(nbdirsmax,0:2*nnucl-1,0:2*nnucl-1)
      integer :: i, j, nd

*    ..recalculate inverse of A
      call invert_matrix(A, Ainv, nnucl, nz2n)

*    ..discrete version of d(inv(A)):
*      A*Ainv = I
*      d(A*Ainv) = Ad*Ainv + A*Ainvd = 0
*      Ainvd = -Ainv*Ad*Ainv
      Ainvd = 0.0_R8
      dum = 0.0_R8
      do nd = 1, nbdirs
*     ..first do dum = -Ainv*Ad
        do i = 0, 2*nnucl-1
          do j = 0, 2*nnucl-1
            dum(nd,i,j) = -sum(Ainv(i,:)*Ad(nd,:,j))
          end do
        end do
*     ..then do Ainvd = dum*Ainv
        do i = 0, 2*nnucl-1
          do j = 0, 2*nnucl-1
            Ainvd(nd,i,j) = sum(dum(nd,i,:)*Ainv(:,j))
          end do
        end do
      end do

      return
      end subroutine invert_matrix_dv

