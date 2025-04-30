
      subroutine solve_covariance(n, z, C, res, LdetC, U, icf)
      use b2mod_types
      use b2mod_subsys
      implicit none
      integer, intent (in) :: n !leading dimensions of vectors/matrices
      real(kind=R8), intent (in) :: z(n) !RHS
      real(kind=R8), intent (in) :: C(n*(n+1)/2) !Covariance matrix in packed format
      real(kind=R8), intent (out) :: res(n) !result of the system
      real(kind=R8), intent (out) :: LdetC !log-determinant of C
      real(kind=R8), intent (out) :: U(n*(n+1)/2) !Cholesky factorization of C (not actually used later on if not for forward/adjoint AD)
      integer, intent (in) :: icf !cost function number (only used for error purposes)
      real(kind=r8) :: dummy
      integer :: ii
!csc This is a dummy routine to be fed to Tapenade isntead of the real solve_covariance, 
!    such that only the necessary dependencies are created
      do ii=1,n
        res(ii) = z(ii)*C(ii)
      end do
      dummy=1.0_R8
      do ii=1,n*(n+1)/2
        dummy = dummy*C(ii)
      enddo
      LdetC = log(dummy)
      U = C
      return
      end subroutine solve_covariance


