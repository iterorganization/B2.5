      subroutine calc_prior(prior,inrange)
      use b2mod_types
      use b2mod_par_opt
      use b2mod_constants
      use b2mod_subsys
      implicit none
      real(kind=R8), intent (inout) :: prior
      logical, intent (inout) :: inrange
      integer :: ii, isigma, nn
      intrinsic sqrt, exp
! sc  Routine based on R. De Wolf et al 2021 Nucl. Fusion 61 046048
      call subini('calc_prior')
! sc  Initialize prior
      prior = 1.0_R8
      inrange = .true.
! sc  Loop on number of plasma parameters (no sigma)
      nn = nsigma
      if (nsigma_opt.gt.0) nn = nsigma_opt
      do ii = 1, npar_opt - nn
        select case (prior_type(ii))
        case (0) !Uniform distribution
          prior = prior*1.0_R8
        case (1) !Uninformative, proper, Gaussian prior. Factor two because only positive parameters are allowed.
          prior = prior*2.0_R8/(prior_par(ii,2)*sqrt(2.0_R8*pi))*
     &            exp(-0.5_R8*((par_opt_phys(ii)-prior_par(ii,1))
     &            /prior_par(ii,2))**2)
        case default
          write(*,*) prior_type(ii)
          call xerrab ('prior_type out of bounds')
        end select
        if (par_opt_phys(ii).lt.prior_range(ii,1) .or.
     &      par_opt_phys(ii).gt.prior_range(ii,2)) then
          prior = 0.0_R8
          inrange = .false.
        endif
      end do

! Now loop only on sigma vector
      isigma = max(0,npar_opt-nn)
      do ii = 1, nsigma_opt
        select case (prior_type(ii+isigma))
        case (0) !Uniform distribution
          prior = prior*1.0_R8
        case (1) !Uninformative, proper, Gaussian prior. Factor two because only positive parameters are allowed.
          prior = prior*2.0_R8/(prior_par(ii+isigma,2)*sqrt(2.0_R8*pi))*
     &            exp(-0.5_R8*((sigma(ii)-prior_par(ii+isigma,1))
     &            /prior_par(ii+isigma,2))**2)
        case (2) !Jeffrey's 1/sigma
          prior = prior/sigma(ii)
        case default
          write(*,*) prior_type(ii+isigma)
          call xerrab ('prior_type out of bounds')
        end select
        if (sigma(ii).lt.prior_range(ii+isigma,1) .or.
     &      sigma(ii).gt.prior_range(ii+isigma,2)) then
          prior = 0.0_R8
          inrange = .false.
        endif
      end do
      
     
      call subend ()
      return
      end subroutine calc_prior
