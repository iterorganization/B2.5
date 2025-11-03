      subroutine calc_log_prior(prior,inrange)
      use b2mod_types
      use b2mod_par_opt
      use b2mod_constants
      use b2mod_subsys
      implicit none
      real(kind=R8), intent (inout) :: prior
      logical, intent (inout) :: inrange
      integer :: ii, isigma, imean, ind
      real(kind=R8) :: aa, bb
#ifndef LEGACYCOMP
      intrinsic log_gamma
#endif
! sc  Routine based on R. De Wolf et al 2021 Nucl. Fusion 61 046048
      call subini('calc_log_prior')
! sc  Initialize prior
      prior = 0.0_R8
      inrange = .true.
! sc  Loop on number of plasma parameters (no sigma, no mean, no shift)
      do ii = 1, npar_opt - nsigma_opt - nmean_opt - nshift_opt -
     &  ncorr_opt
        select case (prior_type(ii))
        case (0) !Uniform distribution
          prior = prior !(+=log(1.0))
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ii,2)) - 0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((par_opt_phys(ii)-prior_par(ii,1))/
     &            prior_par(ii,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          ! first define shape alpha (aa) and rate beta (bb) using mean=aa/bb and variance=aa/bb**2
          ! --> aa = (mu/sigma)**2; bb = mu/sigma**2
          aa = (prior_par(ii,1)/prior_par(ii,2))**2
          bb = prior_par(ii,1)/prior_par(ii,2)**2
          prior = prior + (aa-1)*log(par_opt_phys(ii)) -
     &            bb*par_opt_phys(ii) + aa*log(bb) - log_gamma(aa)
        case default
          write(*,*) prior_type(ii)
          call xerrab ('prior_type out of bounds')
        end select
        if (par_opt_phys(ii).lt.prior_range(ii,1) .or.
     &      par_opt_phys(ii).gt.prior_range(ii,2)) then
          prior = 1.0E+60_R8
          inrange = .false.
        endif
      end do

! Now loop only on sigma vector
      isigma = npar_opt - nsigma_opt-nmean_opt-nshift_opt-ncorr_opt
      ind = 1 !index for prior vectors (may differ from ii due to not optimizing some sigmas)
      do ii = 1, nsigma
        if (.not.sigma_opt(ii)) cycle
        select case (prior_type(ind+isigma))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ind+isigma,2)) -
     &            0.5_R8*log(2.0_R8*pi) - 
     &            0.5_R8*((sigma(ii)-prior_par(ind+isigma,1))/
     &            prior_par(ind+isigma,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          aa = (prior_par(ind+isigma,1)/prior_par(ind+isigma,2))**2
          bb = prior_par(ind+isigma,1)/prior_par(ind+isigma,2)**2
          prior = prior + (aa-1)*log(sigma(ii)) -
     &            bb*sigma(ii) + aa*log(bb) - log_gamma(aa)
        case (3) !Jeffrey's 1/sigma
          prior = prior - log(sigma(ii))
        case default
          write(*,*) prior_type(ind+isigma)
          call xerrab ('prior_type out of bounds for sigma')
        end select
        if (sigma(ii).lt.prior_range(ind+isigma,1) .or.
     &      sigma(ii).gt.prior_range(ind+isigma,2)) then
          prior = 1.0E+60_R8
          inrange = .false.
        endif
        ind = ind + 1
      end do

! Now loop only on mean vector
      imean = npar_opt - nmean_opt - nshift_opt - ncorr_opt
      ind = 1 !index for prior vectors (may differ from ii due to not optimizing some means)
      do ii = 1, nmean
        if (.not.mean_opt(ii)) cycle
        select case (prior_type(ind+imean))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ind+imean,2)) -
     &            0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((mean(ii)-prior_par(ind+imean,1))/
     &            prior_par(ind+imean,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          aa = (prior_par(ind+imean,1)/prior_par(ind+imean,2))**2
          bb = prior_par(ind+imean,1)/prior_par(ind+imean,2)**2
          prior = prior + (aa-1)*log(mean(ii)) -
     &            bb*mean(ii) + aa*log(bb) - log_gamma(aa)
        case default
          write(*,*) prior_type(ind+imean)
          call xerrab ('prior_type out of bounds for mean')
        end select
        if (mean(ii).lt.prior_range(ind+imean,1) .or.
     &      mean(ii).gt.prior_range(ind+imean,2)) then
          prior = 1.0E+60_R8
          inrange = .false.
        endif
        ind = ind + 1
      end do

! Now do shift priors
      ind = 1 !index for prior vectors
      do ii = 1, nshift
        if (.not.shiftopt(ii)) cycle
        select case (shift_prior_type(ind))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(shift_prior_par(ind,2)) -
     &            0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((shift(ii)-shift_prior_par(ind,1))/
     &            shift_prior_par(ind,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          aa = (shift_prior_par(ind,1)/shift_prior_par(ind,2))**2
          bb = shift_prior_par(ind,1)/shift_prior_par(ind,2)**2
          prior = prior + (aa-1)*log(shift(ii)) -
     &            bb*shift(ii) + aa*log(bb) - log_gamma(aa)
        case default
          write(*,*) shift_prior_type(ind)
          call xerrab ('prior_type out of bounds for shift')
        end select
        if (shift(ii).lt.shift_prior_range(ind,1) .or.
     &      shift(ii).gt.shift_prior_range(ind,2)) then
          prior = 1.0E+60_R8
          inrange = .false.
        endif
        ind = ind + 1
      end do

! Now do corr priors
      ind = 1 !index for prior vectors
      do ii = 1, ncf
        if (.not.corr_opt(ii)) cycle
        select case (corr_prior_type(ind))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(corr_prior_par(ind,2)) -
     &            0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((corr_length(ii)-corr_prior_par(ind,1))/
     &            corr_prior_par(ind,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          aa = (corr_prior_par(ind,1)/corr_prior_par(ind,2))**2
          bb = corr_prior_par(ind,1)/corr_prior_par(ind,2)**2
          prior = prior + (aa-1)*log(corr_length(ii)) -
     &            bb*corr_length(ii) + aa*log(bb) - log_gamma(aa)
        case default
          write(*,*) corr_prior_type(ind)
          call xerrab ('prior_type out of bounds for corr')
        end select
        if (corr_length(ii).lt.corr_prior_range(ind,1) .or.
     &      corr_length(ii).gt.corr_prior_range(ind,2)) then
          prior = 1.0E+60_R8
          inrange = .false.
        endif
        ind = ind + 1
      end do


      call subend ()
      return

#ifdef LEGACYCOMP
      contains

      function log_gamma(xx)
      implicit none
      real(kind=R8) :: log_gamma
      real(kind=R8), intent(in) :: xx

C+
C Project     :	Course FORTRAN Codes
C
C Name        :	GAMMLN
C
C Purpose     :	Numerical Recipes natural log of the gamma function.
C
C Explanation :	This function calculates the natural log of the gamma
C               function for XX > 0, where the gamma function, GAMMA(a),
C               is given by:
C
C                   GAMMA(XX) = INT_0^INFTY t**(a-1) EXP(-t) dt.
C
C Use         :	Result = GAMMLN(XX)
C
C Inputs      : XX:     Scalar containing the value passed to the gamma function.
C
C Outputs     :	Result: ln(GAMMA(XX)).
C
C Calls       :	None.
C
C Common      :	None.
C
C Restrictions:	FORTRAN 77 coding.
C
C Side effects:	None.
C
C Category    :	Data fitting.
C
C Prev. Hist. :	Based on the GAMMLN function of Numerical Recipes for FORTRAN 77.
C
C Written     :	Donald G. Luttermoser, ETSU/Physics, 2 October 2013.
C
C Modified    :	Version 1, Donald G. Luttermoser, ETSU/Physics, 2 Oct 2013
C			Initial program.
C
C Version     :	Version 1,  2 October 2013.
C
C-
      INTEGER J
      REAL(KIND=R8) :: SER, STP, TMP, X, Y, COF(6)
      SAVE COF, STP
C
      DATA COF, STP/76.18009172947146_R8, -86.50532032941677_R8,
     +    24.01409824083091_R8, -1.231739572450155_R8,
     +    0.1208650973866179E-2_R8, -0.5395239384953E-5_R8,
     +    2.5066282746310005_R8/
C
      X = XX
      Y = X
      TMP = X + 5.5_R8
      TMP = (X + 0.5_R8) * LOG(TMP) - TMP
      SER = 1.000000000190015_R8
C
      DO 11 J = 1, 6
          Y = Y + 1.0_R8
          SER = SER + COF(J)/Y
 11   CONTINUE
C
      log_gamma = TMP + LOG(STP * SER / X)
C
      RETURN
      end function log_gamma
#endif
      end subroutine calc_log_prior
