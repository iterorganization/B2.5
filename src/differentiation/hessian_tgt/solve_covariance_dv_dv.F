      subroutine solve_covariance_dv_dv(n, z, zd0, zd, zdd, C, Cd0, Cd,
     &  Cdd, res, resd0, resd, resdd, LdetC, LdetCd, LdetCdd,
     &  U, Ud0, icf, nbdirs, nbdirs0)
      use b2mod_types
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
!     This routines solves the system res=C^-1*z also in its forward
!     differentiated form 
      integer, intent (in) :: n !leading dimensions of vectors/matrices
      real(kind=R8), intent (in) :: z(n) !RHS
      real(kind=R8), intent (in) :: zd0(nbdirsmax0,n) !RHS diff'ed
      real(kind=R8), intent (in) :: zd(nbdirsmax,n) !RHS diff'ed
      real(kind=R8), intent (in) :: zdd(nbdirsmax0,nbdirsmax,n) !RHS diff'ed
      real(kind=R8), intent (in) :: C(n*(n+1)/2) !Covariance matrix in packed format
      real(kind=R8), intent (in) :: Cd0(nbdirsmax0,n*(n+1)/2) !Covariance diff matrix in packed format
      real(kind=R8), intent (in) :: Cd(nbdirsmax,n*(n+1)/2) !Covariance diff matrix in packed format
      real(kind=R8), intent (in) :: Cdd(nbdirsmax0,nbdirsmax,n*(n+1)/2) !Covariance diff matrix in packed format
      real(kind=R8), intent (out) :: res(n) !result of the system
      real(kind=R8), intent (out) :: resd0(nbdirsmax0,n) !result of the diff'ed system
      real(kind=R8), intent (out) :: resd(nbdirsmax,n) !result of the diff'ed system
      real(kind=R8), intent (out) :: resdd(nbdirsmax0,nbdirsmax,n) !result of the diff'ed system
      real(kind=R8), intent (out) :: LdetC !log-determinant of C
      real(kind=R8), intent (out) :: LdetCd(nbdirsmax) !log-determinant of C
      real(kind=R8), intent (out) :: LdetCdd(nbdirsmax0,nbdirsmax) !log-determinant of C
      real(kind=R8), intent (out) :: U(n*(n+1)/2) !Cholesky factorization of C (not actually used later on if not for adjoint AD)
      real(kind=R8), intent (out) :: Ud0(nbdirsmax0,n*(n+1)/2)
      integer, intent (in) :: icf !cost function number (only used for error purposes)
      integer, intent (in) :: nbdirs, nbdirs0 !number of tangent directions
!     local variables to avoid aliasing in Lapack routines
      real(kind=R8) :: rr(n), cc(n), yy(nbdirs,n), invC(n*(n+1)/2),
     & yyd0(nbdirsmax0,nbdirs,n)
      real(kind=R8) :: dummy1, dummy2(nbdirsmax0), dummy3(n*(n+1)/2),
     & dum(nbdirsmax, n, n), dum1(nbdirsmax0, nbdirsmax, n, n)
      integer :: info, ii, jj, nd, nd0, j1, j2, i1, i2

      call subini('solve_covariance_dv_dv')
      ! solve_covariance_dv solves the diff'ed form d(C*res=z)
      ! first it solves C*res=z
      ! then C*resd = zd - Cd*res = yy
      ! differentiated again it becomes
      ! d(C*res)=d(z) --> solve_covariance_dv
      ! and then d(C*resd) = d(yy) --> solve_covariance_dv
      ! with d(yy) = yyd0 = zdd - Cdd*res - Cd*resd0

      ! first part: solve the derivative of d(C*res)=d(z)
      call solve_covariance_dv(n, z, zd0, C, Cd0, res, resd0,
     &  LdetC, LdetCd, U, icf, nbdirs0)
      resd(1:nbdirs,:) = resd0(1:nbdirs0,:)

      ! second part, differential of yy = zd - Cd*res
      ! which is yyd0 = zdd - Cdd*res - Cd*resd0
      yy = 0.0_R8
      yyd0 = 0.0_R8
      do nd = 1, nbdirs
        do nd0 = nd, nbdirs0
          yyd0(nd0,nd,:) = zdd(nd0,nd,:)
          ! - Cdd*res
          do ii = 1, n
            ! store in dd the line of coefficients in Cd, starting from end
            do jj = n, ii, -1
              rr(jj) = Cdd(nd0, nd, ii + (jj-1)*jj/2)
            enddo
            ! these following coefficients are swapped, we store only upper triangle of Cd
            do jj = 1, ii-1 
              rr(jj) = Cdd(nd0, nd, jj + (ii-1)*ii/2)
            enddo
            yyd0(nd0,nd,ii) = yyd0(nd0,nd,ii) - sum(rr*res)
          end do
          ! - Cd*resd0
          do ii = 1, n
            ! store in dd the line of coefficients in Cd, starting from end
            do jj = n, ii, -1
              rr(jj) = Cd(nd, ii + (jj-1)*jj/2)
            enddo
            ! these following coefficients are swapped, we store only upper triangle of Cd
            do jj = 1, ii-1 
              rr(jj) = Cd(nd, jj + (ii-1)*ii/2)
            enddo
            yyd0(nd0,nd,ii) = yyd0(nd0,nd,ii) - sum(rr*resd0(nd0,:))
          end do
        end do ! nd0
        yy(nd,:) = zd(nd,:)
        do ii = 1, n
          ! store in dd the line of coefficients in Cd, starting from end
          do jj = n, ii, -1
            rr(jj) = Cd(nd, ii + (jj-1)*jj/2)
          enddo
          ! these following coefficients are swapped, we store only upper triangle of Cd
          do jj = 1, ii-1 
            rr(jj) = Cd(nd, jj + (ii-1)*ii/2)
          enddo
          yy(nd,ii) = yy(nd,ii) - sum(rr*res)
        end do
      end do

      ! third part, solve d(C*resd)=d(yy)
      ! C -->  C
      ! Cd --> Cd0
      ! res --> resd
      ! resd --> resdd
      ! z --> yy
      ! zd --> yyd0
      do nd = 1, nbdirs
        call solve_covariance_dv(n, yy(nd,:), yyd0(:,nd,:), C, Cd0,
     &    resd(nd,:), resdd(:,nd,:), dummy1, dummy2, dummy3, icf,
     &    nbdirs0)
      end do

      LdetCdd = 0.0_R8
!     The 1st derivative of q=log(det(C)) -> p=det(C); q = log(p)
!     pd = d(det(C)) = p*Tr(C^-1*Cd)
!     qd = d(log(p)) = 1/p*pd = Tr(C^-1*Cd)
!     setting A = C^-1*Cd --> qd = Tr(A)
!     The second derivative qdd = d(qd) = d(Tr(A)) = Tr(Ad)
!     From report of M. Giles Ad = C^-1*(Cdd - Cd0*A) and thus
!     qdd = Tr(C^-1*Cdd) - Tr(C^-1*Cd0*C^-1*Cd)

      ! get inverse of covariance matrix for derivative of log-determinant
      invC = U
      call dpptri('U',n,invC,info)
      write(*,*) 'dpptri outcome = ',info
      call xertst(info.eq.0,'solve_covariance_dv: unable to invert '//
     & 'covariance matrix')

      ! Calculate qdd in three steps
      ! 1 - dum = C^-1*Cd
      ! 2 - dum1 = Cd0*dum
      ! 3 - qdd =Tr(C^-1*Cdd) - Tr(C^-1*dum1)

      ! step 1, dum = C^-1*Cd
      ! the result may not be symmetric anymore
      dum = 0.0_R8
      do nd = 1, nbdirs
        do i1 = 1, n
          rr = 0.0_R8 ! row of left matrix
          do j1 = 1, i1-1
            rr(j1) = invC(j1 + (i1-1)*i1/2)
          end do
          do j1 = i1, n
            rr(j1) = invC(i1 + (j1-1)*j1/2)
          end do
          do i2 = 1, n
            cc = 0.0_R8 ! column of right matrix
            do j2 = 1, i2-1
              cc(j2) = Cd(nd, j2 + (i2-1)*i2/2)
            end do
            do j2 = i2, n
              cc(j2) = Cd(nd, i2 + (j2-1)*j2/2)
            end do
            dum(nd, i1, i2) = sum(rr*cc)
          end do ! i2
        end do ! i1
      end do ! nbdirs

      ! step 2, dum1 = Cd0*dum
      ! dum1 may also be non-symmetric
      dum1 = 0.0_R8
      do nd = 1, nbdirs
        do nd0 = 1, nbdirs0
          do i1 = 1, n
            rr = 0.0_R8 ! row of left matrix
            do j1 = 1, i1-1
              rr(j1) = Cd0(nd0, j1 + (i1-1)*i1/2)
            end do
            do j1 = i1, n
              rr(j1) = Cd0(nd0, i1 + (j1-1)*j1/2)
            end do
            ! dum in principle not symmetric anymore
            do j2 = 1, n
              dum1(nd0, nd, i1, j2) = sum(rr*dum(nd, 1:n, j2))
            end do ! i2
          end do ! i1
        end do ! nbdirs0
      end do ! nbdirs

      ! step 3, qdd =Tr(C^-1*Cdd) - Tr(C^-1*dum1)
      LdetCdd = 0.0_R8
      do nd = 1, nbdirs
        do nd0 = 1, nbdirs0
          do ii = 1, n
            rr = 0.0_R8 ! row of invC
            do jj = 1, ii-1
              rr(jj) = invC(jj + (ii-1)*ii/2)
              LdetCdd(nd0,nd) = LdetCdd(nd0,nd) +
     &          invC(jj + (ii-1)*ii/2)*Cdd(nd0, nd, jj + (ii-1)*ii/2)
            enddo
            do jj = ii, n
              rr(jj) = invC(ii + (jj-1)*jj/2)
              LdetCdd(nd0,nd) = LdetCdd(nd0,nd) +
     &          invC(ii + (jj-1)*jj/2)*Cdd(nd0, nd, ii + (jj-1)*jj/2)
            enddo
            ! dum1 not symmetric anymore
            LdetCdd(nd0,nd) = LdetCdd(nd0,nd) - 
     &        sum(rr*dum1(nd0, nd, 1:n, ii))
          end do
          write(*,*) 'MAP LdetCdd', LdetCdd(nd0,nd)
        end do
      end do

      ! set Ud0 to zero (not used in tangent)
      Ud0 = 0.0_R8
      
      call subend ()
      return
      end subroutine solve_covariance_dv_dv
