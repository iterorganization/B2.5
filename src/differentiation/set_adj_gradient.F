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
*-----------------------------------------------------------------------
*     set_adj_gradient is used in b2optim_tao and b2optim_ipopt
*     to correctly set and pass the gradient from adjoint AD to the
*     optimization routines
*-----------------------------------------------------------------------
      integer :: idir, ii, ip, isigma, imean
      character*3 ss1
      character*3 ss2
      character*70 key
      grad = 0.0_R8
      idir = 1 !this indicates the different directions of adjoint mode
      do ip = 1, nnvar - nsigma_opt - nmean_opt
        if (spatial_dep(ip)) then
!         spatially dependent coefficient
          select case (partype(ip))
            case(1)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_dna, species ", ss1
            case(2)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_dpa, species ", ss1
            case(3)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_hci, species ", ss1
            case(4)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_hce, species ", ss1
            case(5)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_vla (poloidal), species ", ss1
            case(6)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_vla (radial), species ", ss1
            case(7)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_vsa, species ", ss1
            case(8)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_sig, species ", ss1
            case(9)
              write(ss1,"(I3)") paris(ip)
              write (key, *) "parm_alf, species ", ss1
            case default
              write(*,*) 'partype = ',partype(ip)
              call xerrab ('partype incompatible with spatial points')
          end select
          do ii = 1, spatial_points(ip)
            grad(idir) = tdatab(2,ii,partype(ip),paris(ip))
            grad(idir) = grad(idir) + par_opt_physb(idir)
            write (*,*) '  Cost function gradient 1: ',
     &        grad(idir),'wrt',trim(key)
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          ! transport coefficients
          case (1)
            grad(idir) = parm_dnab(paris(ip))
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dna, species ", ss1
          case (2)
            grad(idir) = parm_dpab(paris(ip))
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dpa, species ", ss1
          case (3)
            grad(idir) = parm_hcib(paris(ip))
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hci, species ", ss1
          case (4)
            grad(idir) = parm_hceb
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hce, species ", ss1
          case (6)
            grad(idir) = parm_vlab(paris(ip))
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vla, species ", ss1
          case (7)
            grad(idir) = parm_vsab(paris(ip))
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vsa, species ", ss1
          case (8)
            grad(idir) = parm_sigb
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_sig, species ", ss1
          case (9)
            grad(idir) = parm_alfb
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_alf, species ", ss1
          ! boundary conditions
          case (10)
            grad(idir) = eneparb(parib(ip),1)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enepar, boundary ",ss2
          case (11)
            grad(idir) = eniparb(parib(ip),1)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enipar, boundary ",ss2 
          case (12)
            grad(idir) = conparb(paris(ip),parib(ip),1)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "conpar, species ",ss1,", boundary ",ss2
          case (13)
            grad(idir) = momparb(paris(ip),parib(ip),1)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "mompar, species ",ss1,", boundary ",ss2
          case (14)
            grad(idir) = potparb(parib(ip),1)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "potpar, boundary ",ss2
          case (15)
            grad(idir) = enkparb(parib(ip),1)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enkpar, boundary ",ss2
          ! recycling
          case (16)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "b2recyc, species ",ss1,", stratum ",ss2
          ! ballooning
          case (17)
            grad(idir) = switchdiff%b2tqna_ballooning
            write (key, *) "b2tqna_ballooning"
          case (18)
            grad(idir) = switchdiff%b2tqna_ballooning_rescale
            write (key, *) "b2tqna_ballooning_rescale"
          ! k-model switches
          case (19)
            grad(idir) = switchdiff%keps_cd
            write (key, *) "keps_cd"
          case (20)
            grad(idir) = switchdiff%keps_heat
            write (key, *) "keps_heat"
          case (21)
            grad(idir) = switchdiff%keps_heat_i
            write (key, *) "keps_heat_i"
          case (22)
            grad(idir) = switchdiff%keps_sig
            write (key, *) "keps_sig"
          case (23)
            grad(idir) = switchdiff%keps_alf
            write (key, *) "keps_alf"
          case (24)
            grad(idir) = switchdiff%keps_visc
            write (key, *) "keps_visc"
          case (25)
            grad(idir) = switchdiff%keps_dkt
            write (key, *) "keps_dkt"
          case (26)
            grad(idir) = switchdiff%keps_dzt
            write (key, *) "keps_dzt"
          case (27)
            grad(idir) = switchdiff%b2sikt_fac_diss
            write (key, *) "b2sikt_fac_diss"
          case (28)
            grad(idir) = switchdiff%b2sikt_fac_diss_core 
            write (key, *) "b2sikt_fac_diss_core"
          case (29)
            grad(idir) = switchdiff%b2sikt_fac_sheath
            write (key, *) "b2sikt_fac_sheath"
          case (30)
            grad(idir) = switchdiff%b2sikt_fac_sheath_core
            write (key, *) "b2sikt_fac_sheath_core"
          case (31)
            grad(idir) = switchdiff%keps_shear
            write (key, *) "keps_shear"
          case (32)
            grad(idir) = switchdiff%b2tfhi_fsigkt
            write (key, *) "b2tfhi_fsigkt"
          case (33)
            grad(idir) = switchdiff%b2sikt_fac_vis_RS
            write (key, *) "b2sikt_fac_vis_RS"
          case (34)
            grad(idir) = switchdiff%b2tfhi_fflokt
            write (key, *) "b2tfhi_fflokt"
          case (35)
            grad(idir) = switchdiff%b2tfhi_fconkt
            write (key, *) "b2tfhi_fconkt"
          case (36)
            grad(idir) = switchdiff%b2tfhi_fflozt
            write (key, *) "b2tfhi_fflozt"
          case (37)
            grad(idir) = switchdiff%b2tfhi_fconzt
            write (key, *) "b2tfhi_fconzt"
          case (38)
            grad(idir) = switchdiff%b2tfhi_fkt_hie
            write (key, *) "b2tfhi_fkt_hie"
          case (39)
            grad(idir) = switchdiff%b2tfhe_vis_kt
            write (key, *) "b2tfhe_vis_kt"
          case default
            write(*,*) 'partype = ',partype(ip)
            call xerrab ('partype out of bounds')
          end select
          grad(idir) = grad(idir) + par_opt_physb(idir)
          write (*,*) '  Cost function gradient 1: ',
     &      grad(idir),'wrt',trim(key)
          idir = idir + 1
        endif
      end do

      if (nsigma_opt.gt.0) then
        do ip = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ip)) then
            grad(idir) = sigmab(ip)
            write(ss2,"(I3)") ip
            write (key, *) "sigma",ss2
            write (*,*) '  Cost function gradient 1: ',
     &        grad(idir),'wrt',trim(key)
            idir = idir + 1
          endif
        end do
      endif

      if (nmean_opt.gt.0) then
        do ip = 1, nmean
!         only if mean is being optimized!
          if (mean_opt(ip)) then
            grad(idir) = meanb(ip)
            write(ss2,"(I3)") ip
            write (key, *) "mean",ss2
            write (*,*) '  Cost function gradient 1: ',
     &        grad(idir),'wrt',trim(key)
            idir = idir + 1
          endif
        end do
      endif

      return
      end subroutine set_adj_gradient

