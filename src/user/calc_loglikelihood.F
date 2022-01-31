      subroutine calc_loglikelihood(nn, ff, yy, lll, isigma)
      use b2mod_types
      use b2mod_user_namelist
     & , only : sigma, scale_sigma
      use b2mod_constants
      implicit none
      integer, intent (in) :: nn, isigma
      real(kind=R8), intent (in) :: ff(nn) !SOLPS results
      real(kind=R8), intent (in) :: yy(nn) !Data
      real(kind=R8), intent (inout) :: lll !loglikelihood
      integer :: ii
      real(kind=R8) :: zz(nn), RR(nn), invC_z(nn)
      intrinsic sqrt, exp, log
      external subini, subend
! sc  Routine based on R. De Wolf et al 2021 Nucl. Fusion 61 046048
      call subini('calc_loglikelihood')
! sc  TO BE DONE: proper covariance matrix should be built here
!     for the moment no correlation => only diagonal elements
      if (scale_sigma) then
        RR = (yy*sigma(isigma))**2 ! diagonal elements
      else
        RR = sigma(isigma)**2
      endif
      zz = ff - yy ! FIXME in principle also -mu
!     solution of quadratic problem for DIAGONAL matrix
! FIXME in principle invC_z = Cov\zz = SIGMA^-1*zz
      do ii = 1, nn
        invC_z(ii) = zz(ii)**2/RR(ii)
      end do 
!     rest = Z'*invCov_z
!     Calculate log-likelihood
      RR = sqrt(RR) ! FIXME, should go away when real cholezky factorization is there
      lll = -0.5_R8*(nn*log(2.0_R8*pi) + sum(2.0_R8*log(RR)) 
     &               + sum(invC_z))
      
     
      call subend ()
      return
      end subroutine calc_loglikelihood
