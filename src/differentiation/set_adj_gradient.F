      subroutine set_adj_gradient(npar,grad,switchdiff)
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_transport_namelist_diff
      use b2mod_input_profile_diff
      use b2mod_boundary_namelist_diff
      use b2mod_switches_diff
      implicit none
      integer :: npar
      real(kind=R8) :: grad(npar)
      type(switches), intent(inout) :: switchdiff
      integer :: idir, ii, ip, isigma

      grad = 0.0_R8
      idir = 1 !this indicates the different directions of adjoint mode
      do ip = 1, nnvar - nsigma_opt
        if (spatial_dep(ip)) then
!         spatially dependent coefficient
          do ii = 1, spatial_points(ip)
            grad(idir) = tdatab(2,ii,partype(ip),paris(ip))
            grad(idir) = grad(idir) + par_opt_physb(idir)
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          case (1) ! dna
            grad(idir) = parm_dnab(paris(ip))
          case (2) ! dpa
            grad(idir) = parm_dpab(paris(ip))
          case (3) ! hci
            grad(idir) = parm_hcib(paris(ip))
          case (4) ! hce
            grad(idir) = parm_hceb
          case (5) ! vla
            grad(idir) = parm_vlab(paris(ip))
          case (7) ! vsa
            grad(idir) = parm_vsab(paris(ip))
          case (8) ! sig
            grad(idir) = parm_sigb
          case (9) ! alf
            grad(idir) = parm_alfb
          case (10) ! enkpar(1,1)
            grad(idir) = enkparb(parib(ip),1)
          case (11) ! b2sikt_fac_sheath
            grad(idir) = switchdiff%b2sikt_fac_sheath
          case (12) ! b2sikt_fac_sheath_core
            grad(idir) = switchdiff%b2sikt_fac_sheath_core
          case (13) ! keps_cd
            grad(idir) = switchdiff%keps_cd
          case (14) ! keps_heat
            grad(idir) = switchdiff%keps_heat
          case (15) ! keps_heat_i
            grad(idir) = switchdiff%keps_heat_i
          case (16) ! enepar
            grad(idir) = eneparb(parib(ip),1)
          case (17) ! enipar
            grad(idir) = eniparb(parib(ip),1) 
          case (18) ! conpar
            grad(idir) = conparb(paris(ip),parib(ip),1)
          case (19) ! mompar
            grad(idir) = momparb(paris(ip),parib(ip),1)
          case (20) ! potpar
            grad(idir) = potparb(parib(ip),1) 
          case default
            write(*,*) partype(ip)
            call xerrab ('partype out of bounds')
          end select
          grad(idir) = grad(idir) + par_opt_physb(idir)
          idir = idir + 1
        endif
      end do
      if (nsigma_opt.gt.0) then
        do ip = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ip)) then
            grad(idir) = sigmab(ip)
            idir = idir + 1
          endif
        end do
      endif

      return
      end subroutine set_adj_gradient

