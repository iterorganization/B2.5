      subroutine invert_matrix_b (A, Ab, Ainv, Ainvb, nnucl, nz2n) 
      use b2mod_types
      implicit none
      integer :: nnucl
      real(kind=R8) :: A(0:2*nnucl-1,0:2*nnucl-1), nz2n(0:nnucl-1),
     &                 Ainv(0:2*nnucl-1,0:2*nnucl-1)
      real(kind=R8) :: Ab(0:2*nnucl-1,0:2*nnucl-1),
     &                 Ainvb(0:2*nnucl-1,0:2*nnucl-1),
     &                 dum(0:2*nnucl-1,0:2*nnucl-1)
      integer :: i, j

*    ..recalculate inverse of A
      call invert_matrix(A, Ainv, nnucl, nz2n)

*    ..discrete adjoint to get Ab from Ainv and Ainvb
*      Ab = -A**-T*Cb*A**-T = -Ainv**T*Ainvb*Ainv**T
*    ..first do dum = -Ainv**T*Ainvb
      dum = 0.0_R8
      Ab = 0.0_R8
      do i = 0, 2*nnucl-1
        do j = 0, 2*nnucl-1
          dum(i,j) = -sum(Ainv(:,i)*Ainvb(:,j))
        end do
      end do
*    ..then do Ab = dum*Ainv**T
      do i = 0, 2*nnucl-1
        do j = 0, 2*nnucl-1
          Ab(i,j) = sum(dum(i,:)*Ainv(j,:))
        end do
      end do
      Ainvb = 0.0_R8

      return
      end subroutine invert_matrix_b

