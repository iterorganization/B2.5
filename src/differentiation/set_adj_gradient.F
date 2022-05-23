      subroutine set_adj_gradient(npar,grad,switchdiff)
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_transport_namelist_diff
      use b2mod_input_profile_diff
      use b2mod_boundary_namelist_diff
      use b2mod_neutrals_namelist_diff
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
          ! transport coefficients
          case (1)
            grad(idir) = parm_dnab(paris(ip))
          case (2)
            grad(idir) = parm_dpab(paris(ip))
          case (3)
            grad(idir) = parm_hcib(paris(ip))
          case (4)
            grad(idir) = parm_hceb
          case (5)
            grad(idir) = parm_vlab(paris(ip))
          case (7)
            grad(idir) = parm_vsab(paris(ip))
          case (8)
            grad(idir) = parm_sigb
          case (9)
            grad(idir) = parm_alfb
          ! boundary conditions
          case (10)
            grad(idir) = eneparb(parib(ip),1)
          case (11)
            grad(idir) = eniparb(parib(ip),1) 
          case (12)
            grad(idir) = conparb(paris(ip),parib(ip),1)
          case (13)
            grad(idir) = momparb(paris(ip),parib(ip),1)
          case (14)
            grad(idir) = potparb(parib(ip),1) 
          case (15)
            grad(idir) = enkparb(parib(ip),1)
          ! recycling
          case (16)
            grad(idir) = b2recycb(paris(ip),parib(ip))
          ! ballooning
          case (17)
            grad(idir) = switchdiff%b2tqna_ballooning
          case (18)
            grad(idir) = switchdiff%b2tqna_ballooning_rescale
          ! k-model switches
          case (19)
            grad(idir) = switchdiff%keps_cd
          case (20)
            grad(idir) = switchdiff%keps_heat
          case (21)
            grad(idir) = switchdiff%keps_heat_i
          case (22)
            grad(idir) = switchdiff%keps_sig
          case (23)
            grad(idir) = switchdiff%keps_alf
          case (24)
            grad(idir) = switchdiff%keps_visc
          case (25)
            grad(idir) = switchdiff%keps_dkt
          case (26)
            grad(idir) = switchdiff%keps_dzt
          case (27)
            grad(idir) = switchdiff%b2sikt_fac_diss
          case (28)
            grad(idir) = switchdiff%b2sikt_fac_diss_core 
          case (29)
            grad(idir) = switchdiff%b2sikt_fac_sheath
          case (30)
            grad(idir) = switchdiff%b2sikt_fac_sheath_core
          case (31)
            grad(idir) = switchdiff%keps_shear
          case (32)
            grad(idir) = switchdiff%b2tfhi_fsigkt
          case (33)
            grad(idir) = switchdiff%b2sikt_fac_vis_RS
          case (34)
            grad(idir) = switchdiff%b2tfhi_fflokt
          case (35)
            grad(idir) = switchdiff%b2tfhi_fconkt
          case (36)
            grad(idir) = switchdiff%b2tfhi_fflozt
          case (37)
            grad(idir) = switchdiff%b2tfhi_fconzt
          case (38)
            grad(idir) = switchdiff%b2tfhi_fkt_hie
          case (39)
            grad(idir) = switchdiff%b2tfhe_vis_kt
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

