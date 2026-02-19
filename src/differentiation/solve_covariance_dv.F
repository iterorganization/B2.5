      subroutine solve_covariance_dv(n, z, zd, C, Cd, res, resd,
     &  LdetC, LdetCd, U, icf, nbdirs)
      use b2mod_types
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
!     This routines solves the system res=C^-1*z also in its forward
!     differentiated form 
      integer, intent (in) :: n !leading dimension of vectors/matrices
      real(kind=R8), intent (in) :: z(n) !RHS
      real(kind=R8), intent (in) :: zd(nbdirsmax,n) !RHS diff'ed
      real(kind=R8), intent (in) :: C(n*(n+1)/2) !Covariance matrix in packed format
      real(kind=R8), intent (in) :: Cd(nbdirsmax, n*(n+1)/2) !Covariance diff matrix in packed format
      real(kind=R8), intent (out) :: res(n) !result of the system
      real(kind=R8), intent (out) :: resd(nbdirsmax, n) !result of the diff'ed system
      real(kind=R8), intent (out) :: LdetC !log-determinant of C
      real(kind=R8), intent (out) :: LdetCd(nbdirsmax) !log-determinant of C
      real(kind=R8), intent (out) :: U(n*(n+1)/2) !Cholesky factorization of C (not actually used later on if not for adjoint AD)
      integer, intent (in) :: icf !cost function number (only used for error purposes)
      integer, intent (in) :: nbdirs !number of tangent directions
!     local variables to avoid aliasing in Lapack routines
      real(kind=R8) :: dd(n), yy(n), invC(n*(n+1)/2)
      integer :: info, ii, jj, nd

      call subini('solve_covariance_dv')

      ! perform cholesky factorization CC=U'*U, we need it later!
      res = 0.0_R8
      LdetC = 0.0_R8
      LdetCd = 0.0_R8
      U = C
      call dpptrf('U',n,U,info)
      write(*,*) 'dpptrf outcome = ',info
      if (.not. info.eq.0) then
        write(*,*) "Covariance matrix for CF",icf
        write(*,*) "Format C(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;"
        do ii = 1, (n*(n+1)/2)
          write(*,*) C(ii)
        end do
        call xerrab('solve_covariance_dv: unable to calculate'//
     &   ' Cholesky factorization of covariance matrix')
      endif
      
      ! solve original system C*res = zz
      res = z
      call dpptrs('U',n,1,U,res,n,info)
      write(*,*) 'dpptrs outcome = ',info
      call xertst(info.eq.0,'solve_covariance_dv: unable to solve '//
     & 'system of equations with covariance matrix')

      ! get inverse of covariance matrix for derivative of log-determinant
      invC = U
      call dpptri('U',n,invC,info)
      write(*,*) 'dpptri outcome = ',info
      call xertst(info.eq.0,'solve_covariance_dv: unable to invert '//
     & 'covariance matrix')

      ! calculate determinant of C
      ! set dd = diagonal of U
      !  if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
      do ii = 1, n
        dd(ii) = U(ii + (ii-1)*ii/2)
      end do
      ! log(detC) = log(PROD(U_ii**22))
      LdetC = sum(2.0_R8*log(dd))
      write(*,*) 'MAP LdetC', LdetC

      resd = 0.0_R8
      do nd = 1, nbdirs
        ! Now solve the differentiated form Cd*res + C*resd = zd for resd
        ! which becomes a system in the form C*resd = zd - Cd*res = yy
        yy = zd(nd,:)
        ! The covariance matrix can be affected only by perturbation
        ! in sigma or corr_length (and shift?). For other parameters
        ! the perturbation is zero
        do ii = 1, n
          ! store in dd the line of coefficients in Cd, starting from end
          do jj = n, ii, -1
            dd(jj) = Cd(nd, ii + (jj-1)*jj/2)
          enddo
          ! these following coefficients are swapped, we store only upper triangle of Cd
          do jj = 1, ii-1
            dd(jj) = Cd(nd, jj + (ii-1)*ii/2)
          enddo
          yy(ii) = yy(ii) - sum(dd*res)
        end do
        call dpptrs('U',n,1,U,yy,n,info)
        write(*,*) 'dpptrs outcome = ',info, 'for direction ',nd
        call xertst(info.eq.0,'solve_covariance_dv: unable to solve '//
     &   'diffed system of equations with covariance matrix')
        resd(nd,:) = yy
        ! Now we still need derivative of log(det(C))
        ! from M.Giles Report no. 08/01 "An extended collection of matrix ..."
        ! p = det(C) --> pd = p*Tr(C**-1*Cd)
        ! q = log(p) --> qd = pd/p = Tr(C**-1*Cd)
        ! only sigma and corr_length have effect of covariance!
        if (any(abs(Cd(nd,:)) .gt. 0.0_R8)) then
          ! since matrices are symmetric we can take directly the trace as 1strow*1strow, etc
          do ii = 1, n
            do jj = n, ii, -1
              LdetCd(nd) = LdetCd(nd) + invC(ii + (jj-1)*jj/2)*
     &                                  Cd(nd, ii + (jj-1)*jj/2)
            enddo
            ! these following coefficients are swapped
            do jj = 1, ii-1
              LdetCd(nd) = LdetCd(nd) + invC(jj + (ii-1)*ii/2)*
     &                                  Cd(nd, jj + (ii-1)*ii/2)
            enddo
          enddo
        endif
        write(*,*) 'MAP LdetCd', LdetCd(nd)
      end do

      call subend ()
      return
      end subroutine solve_covariance_dv
