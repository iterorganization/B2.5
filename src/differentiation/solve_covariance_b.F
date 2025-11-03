      subroutine solve_covariance_b(n, z, zb, C, Cb, res, resb,
     &  LdetC, LdetCb, U, icf)
      use b2mod_types
      use b2mod_subsys
      implicit none
!     This routines solves the system res=C^-1*z in its adjoint form 
      integer, intent (in) :: n !leading dimensions of vectors/matrices
      real(kind=R8), intent (in)    :: z(n) !RHS
      real(kind=R8), intent (inout) :: zb(n) !RHS diff'ed
      real(kind=R8), intent (in)    :: C(n*(n+1)/2) !Covariance matrix in packed format
      real(kind=R8), intent (inout) :: Cb(n*(n+1)/2) !Covariance diff matrix in packed format
      real(kind=R8), intent (in)    :: res(n) !result of the original system, unchanged from original call
      real(kind=R8), intent (inout) :: resb(n) !result of the diff'ed system
      real(kind=R8), intent (in)    :: LdetC !log-determinant of C, unchanged from original call
      real(kind=R8), intent (inout) :: LdetCb !log-determinant of C
      real(kind=R8), intent (in)    :: U(n*(n+1)/2) !Cholesky factorization of C, unchanged from original call
      integer, intent (in)          :: icf !cost function number (only used for error purposes)
!     local variables to avoid aliasing in Lapack routines
      real(kind=R8) :: invC(n*(n+1)/2), dd(n), bb(n)
      real(kind=r8), parameter :: eps_c = 1.0e-30
      integer :: info, ii, jj

      call subini('solve_covariance_b')

      ! in principle the adjoint requires first adjoint of q = log(d) --> db = qb/d
      ! which means detCb = LdetCb/detC. However down below we actually already do
      ! detCb*detC = LdetCb. --> no need for recalculating the determinant. This
      ! avoids issues with detC close to zero or too large

      ! initialize Cb
      Cb = 0.0_R8

      ! adjoint piece coming from determinant d=det(C) --> Cb = db*d*C^-T
      ! first calculate inverse of C using U
      invC = U
      call dpptri('U',n,invC,info)
      write(*,*) 'dpptri outcome = ',info
      call xertst(info.eq.0,'solve_covariance_b: unable to invert '//
     & 'covariance matrix')
      ! Now perform Cb = db*d*C^-T (C transpose == C since it is symmetric!)
      do ii = 1, (n*(n+1)/2)
        Cb(ii) = Cb(ii) + LdetCb*invC(ii)
      end do

      ! now adjoint of matrix solve A*x=b --> bb = (A^-1)^T*xb and Ab(i,j) = -x(j)*bb(i)
      ! with C*res=bb  --> bb = invC^T*resb = invC*resb..
      do ii = 1, n
        ! store in dd the line of coefficients in invC, starting from end
        do jj = n, ii, -1
          dd(jj) = invC(ii + (jj-1)*jj/2)
        enddo
        ! these following coefficients are swapped, we store only upper triangle of invC
        do jj = 1, ii-1 
          dd(jj) = invC(jj + (ii-1)*ii/2)
        enddo
        bb(ii) = sum(dd*resb)
      end do
      zb = zb + bb
      ! ..and now Cb(i,j) = Cb(i,j) + res(j)*bb(i) with  A(i,j)= AP(i + (j-1)*j/2)  for 1<=i<=j
      do jj = 1, n
        do ii = 1, jj
          Cb(ii + (jj-1)*jj/2) = Cb(ii + (jj-1)*jj/2) - res(jj)*bb(ii)
          ! IMPORTANT: we need to increase by factor 2 the off-diagonal elements of Cb because only the upper triangular part
          ! is stored and used to build CC, otherwise the contribution of the lower triangle to the sigma/corr_length will not the accounted for!
          if (.not. (ii.eq.jj) ) 
     &      Cb(ii + (jj-1)*jj/2) = 2.0_R8*Cb(ii + (jj-1)*jj/2)
        end do
      end do

      ! finalize resb = 0
      resb = 0.0_R8

      call subend ()
      return
      end subroutine solve_covariance_b
