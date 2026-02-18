      subroutine solve_covariance(n, z, C, res, LdetC, U, icf)
      use b2mod_types
      use b2mod_subsys
      implicit none
!     This routines solves the system res=C^-1*z
!     that arises in the Gaussian Likelihood function.
!     C must be a positive definite matrix stored in packed
!     format according to https://netlib.org/lapack/explore-html/de/d49/group__double_o_t_h_e_rsolve_ga59a299a18cb0af017ee2c972ae7defb5.html
      integer, intent (in) :: n !leading dimension of vectors/matrices
      real(kind=R8), intent (in) :: z(n) !RHS
      real(kind=R8), intent (in) :: C(n*(n+1)/2) !Covariance matrix in packed format
      real(kind=R8), intent (out) :: res(n) !result of the system
      real(kind=R8), intent (out) :: LdetC !log-determinant of C
      real(kind=R8), intent (out) :: U(n*(n+1)/2) !Cholesky factorization of C (not actually used later on if not for adjoint AD)
      integer, intent (in) :: icf !cost function number (only used for error purposes)
!     local variables to avoid aliasing in dppsv
      real(kind=R8) :: dd(n)
      integer :: info, ii

      call subini('solve_covariance')

      ! perform cholesky factorization C=U'*U and solve C*res = zz
      res = z
      U = C
      call dppsv('U',n,1,U,res,n,info)
      write(*,*) 'DPPSV outcome = ',info
      if (.not. info.eq.0) then
        write(*,*) "Covariance matrix for CF",icf
        write(*,*) "Format C(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;"
        do ii = 1, (n*(n+1)/2)
          write(*,*) C(ii)
        end do
        call xerrab('solve_covariance: unable to solve '//
     &   'system of equations with covariance matrix')
      endif

      ! calculate determinant of C from Cholesky matrix U
      ! set dd = diagonal of U
      !  if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
      do ii = 1, n
        dd(ii) = U(ii + (ii-1)*ii/2)
      end do
      ! log(detC) = log(PROD(U_ii**2))
      LdetC = sum(2.0_R8*log(dd))

      call subend ()
      return
      end subroutine solve_covariance
