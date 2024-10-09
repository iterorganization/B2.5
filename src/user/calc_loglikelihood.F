      subroutine calc_loglikelihood(nn, ff, xx, yy, ss, lll, isigma,
     &  imean, icf)
      use b2mod_types
      use b2mod_par_opt
     & , only : sigma, scale_sigma, mean, nmean, corr_model,
     &          corr_length, corr_cutoff
      use b2mod_constants
      use b2mod_subsys
      implicit none
      integer, intent (in) :: nn, isigma, imean, icf
      real(kind=R8), intent (in) :: ff(nn) !SOLPS results
      real(kind=R8), intent (in) :: xx(nn) !Spatial points
      real(kind=R8), intent (in) :: yy(nn) !Data
      real(kind=R8), intent (in) :: ss(nn) !Standard deviation
      real(kind=R8), intent (inout) :: lll !loglikelihood
      integer :: ii, jj, ind
      real(kind=R8) :: zz(nn), RR(nn), invC_z(nn), CC(nn*(nn+1)/2),
     &                 val, LdetC, UU(nn*(nn+1)/2)
      intrinsic sqrt, log

! sc  Routine based on R. De Wolf et al 2021 Nucl. Fusion 61 046048
      call subini('calc_loglikelihood')

      if (scale_sigma(isigma)) then
        RR = (yy*sigma(isigma) + ss)**2 ! diagonal elements
      else
        RR = (sigma(isigma) + ss)**2
      endif
      call xertst(any(RR .gt. 0.0_R8), ! FIXME add some lower threshold? 
     & 'calc_loglikelihood: zero or negative variances not allowed!')
      zz = ff - yy
      if (nmean.gt.0) zz = zz - mean(imean)
!     solution of quadratic problem z'*C^-1*z
      select case (corr_model(icf))

      case (0) ! simple diagonal covariance matrix
        do ii = 1, nn
          invC_z(ii) = zz(ii)**2/RR(ii)
        end do
        LdetC = sum(log(RR))

      case (1) ! we have a covariance matrix C
        ! create covariance C packed columnwise in a linear array. 
        ! The j-th column of C is stored in the array as follows
        ! if UPLO = 'U', AP(i + (j-1)*j/2) = A(i,j) for 1<=i<=j;
        ! Two-dimensional storage of the symmetric matrix A when N = 4, UPLO = 'U':
        !     a11 a12 a13 a14
        !         a22 a23 a24
        !             a33 a34     (aij = conjg(aji))
        !                 a44
        !  Packed storage of the upper triangle of A:
        !  AP = [ a11, a12, a22, a13, a23, a33, a14, a24, a34, a44 ]
        CC = 0.0_R8
        ! start from the end and go backwards is easier
        ind = nn*(nn+1)/2
        do ii = nn, 1 , -1
          !main diagonal element
          CC(ind) = RR(ii)
          ind = ind - 1
          ! off diagonal elements
          do jj = ii-1, 1, -1
#ifdef TEST_MAP
            val = exp(-abs(xx(ii)-xx(jj))/corr_length(icf))
#else
            val = exp(-abs(xx(ii)-xx(jj))/corr_length(icf)*1.0e3_R8) !corr length is in [mm]!
#endif
            if (val .lt. corr_cutoff) val = 0.0_R8
            CC(ind) = val*sqrt(RR(ii)*RR(jj))
            ind = ind - 1
          end do
        enddo
        ! solve C*invC_z = zz and get log-determinant
        call solve_covariance(nn, zz, CC, invC_z, LdetC, UU, icf)
        ! multiply invC_z = z'*invC_z
        invC_z = invC_z*zz

      case default
       write(*,*) 'icf, corr_model', icf, corr_model(icf)
       call xerrab ('corr_model out of bounds')

      end select

!     Calculate log-likelihood
      lll = -0.5_R8*(nn*log(2.0_R8*pi) + LdetC + sum(invC_z))

      call subend ()
      return
      end subroutine calc_loglikelihood
