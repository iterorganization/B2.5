      subroutine calc_log_prior(prior,inrange)
      use b2mod_types
      use b2mod_par_opt
      use b2mod_constants
      use b2mod_subsys
      implicit none
      real(kind=R8), intent (inout) :: prior
      logical, intent (inout) :: inrange
      integer :: ii, isigma, imean
      intrinsic sqrt, exp, log_gamma
! sc  Routine based on R. De Wolf et al 2021 Nucl. Fusion 61 046048
      call subini('calc_log_prior')
! sc  Initialize prior
      prior = 0.0_R8
      inrange = .true.
! sc  Loop on number of plasma parameters (no sigma, no mean)
      do ii = 1, npar_opt - nsigma_opt - nmean_opt
        select case (prior_type(ii))
        case (0) !Uniform distribution
          prior = prior !(+=log(1.0))
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ii,2)) - 0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((par_opt_phys(ii)-prior_par(ii,1))/
     &            prior_par(ii,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          prior = prior + (prior_par(ii,1)-1)*log(par_opt_phys(ii)) -
     &            prior_par(ii,2)*par_opt_phys(ii) + 
     &            prior_par(ii,1)*log(prior_par(ii,2)) -
     &            log_gamma(prior_par(ii,1))
        case default
          write(*,*) prior_type(ii)
          call xerrab ('prior_type out of bounds')
        end select
        if (par_opt_phys(ii).lt.prior_range(ii,1) .or.
     &      par_opt_phys(ii).gt.prior_range(ii,2)) then
          prior = log(0.0_R8)
          inrange = .false.
        endif
      end do

! Now loop only on sigma vector
      isigma = npar_opt - nsigma_opt - nmean_opt
      do ii = 1, nsigma_opt
        select case (prior_type(ii+isigma))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ii+isigma,2)) -
     &            0.5_R8*log(2.0_R8*pi) - 
     &            0.5_R8*((sigma(ii)-prior_par(ii+isigma,1))/
     &            prior_par(ii+isigma,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          prior = prior + (prior_par(ii+isigma,1)-1)*log(sigma(ii)) -
     &            prior_par(ii+isigma,2)*sigma(ii) + 
     &            prior_par(ii+isigma,1)*log(prior_par(ii+isigma,2)) -
     &            log_gamma(prior_par(ii,1))
        case (3) !Jeffrey's 1/sigma
          prior = prior - log(sigma(ii))
        case default
          write(*,*) prior_type(ii+isigma)
          call xerrab ('prior_type out of bounds for sigma')
        end select
        if (sigma(ii).lt.prior_range(ii+isigma,1) .or.
     &      sigma(ii).gt.prior_range(ii+isigma,2)) then
          prior = log(0.0_R8)
          inrange = .false.
        endif
      end do

! Now loop only on mean vector
      imean = npar_opt - nmean_opt
      do ii = 1, nmean_opt
        select case (prior_type(ii+imean))
        case (0) !Uniform distribution
          prior = prior
        case (1) !Uninformative, proper, Gaussian prior.
          prior = prior - log(prior_par(ii+imean,2)) -
     &            0.5_R8*log(2.0_R8*pi) -
     &            0.5_R8*((mean(ii)-prior_par(ii+imean,1))/
     &            prior_par(ii+imean,2))**2
        case (2) ! gamma distribution x**(alpha-1)*exp(-beta*xx)*beta**alpha/Gamma(alpha)
          prior = prior + (prior_par(ii+imean,1)-1)*log(mean(ii)) -
     &            prior_par(ii+imean,2)*mean(ii) + 
     &            prior_par(ii+imean,1)*log(prior_par(ii+imean,2)) -
     &            log_gamma(prior_par(ii+imean,1))
        case default
          write(*,*) prior_type(ii+imean)
          call xerrab ('prior_type out of bounds for mean')
        end select
        if (mean(ii).lt.prior_range(ii+imean,1) .or.
     &      mean(ii).gt.prior_range(ii+imean,2)) then
          prior = log(0.0_R8)
          inrange = .false.
        endif
      end do
      
     
      call subend ()
      return
      end subroutine calc_log_prior
