      subroutine set_tgt_perturbation(switchdiff)
      use b2mod_types
      use b2mod_par_opt_diffv
      use b2mod_transport_namelist_diffv
      use b2mod_input_profile_diffv
      use b2mod_boundary_namelist_diffv
      use b2mod_neutrals_namelist_diffv
      use b2mod_switches_diffv
      implicit none
      type(switches_diffv), intent(inout) :: switchdiff
*-----------------------------------------------------------------------
*     set_tgt_perturbation is used in b2optim_tao and b2optim_ipopt
*     to correctly set and pass the input perturbation to tangent AD
*     coming from the optimization routines
*-----------------------------------------------------------------------
      integer :: idir, ii, ip, isigma
!     FIXME initialization to 0.0_R8 should be done elsewhere?
!     FIXME if not differentiated wrt certain variables, then their derivative here will not exist and such the lines should be manually commented-out
      tdatad = 0.0_R8
      parm_dnad = 0.0_R8
      parm_dpad = 0.0_R8
      parm_hcid = 0.0_R8
      parm_hced = 0.0_R8
      parm_vlad = 0.0_R8
      parm_vsad = 0.0_R8
      parm_sigd = 0.0_R8
      parm_alfd = 0.0_R8
      conpard = 0.0_R8
      mompard = 0.0_R8
      potpard = 0.0_R8
      enepard = 0.0_R8
      enipard = 0.0_R8
      enkpard = 0.0_R8
      switchdiff%b2tqna_ballooning = 0.0_R8
      switchdiff%b2tqna_ballooning_rescale = 0.0_R8
      switchdiff%keps_cd = 0.0_R8
      switchdiff%keps_heat = 0.0_R8
      switchdiff%keps_heat_i = 0.0_R8
      switchdiff%keps_sig = 0.0_R8
      switchdiff%keps_alf = 0.0_R8
      switchdiff%keps_visc = 0.0_R8
      switchdiff%keps_dkt = 0.0_R8
      switchdiff%keps_dzt = 0.0_R8
      switchdiff%b2sikt_fac_diss = 0.0_R8
      switchdiff%b2sikt_fac_diss_core = 0.0_R8 
      switchdiff%b2sikt_fac_sheath = 0.0_R8
      switchdiff%b2sikt_fac_sheath_core = 0.0_R8
      switchdiff%keps_shear = 0.0_R8
      switchdiff%b2tfhi_fsigkt = 0.0_R8
      switchdiff%b2sikt_fac_vis_RS = 0.0_R8
      switchdiff%b2tfhi_fflokt = 0.0_R8
      switchdiff%b2tfhi_fconkt = 0.0_R8
      switchdiff%b2tfhi_fflozt = 0.0_R8
      switchdiff%b2tfhi_fconzt = 0.0_R8
      switchdiff%b2tfhi_fkt_hie = 0.0_R8
      switchdiff%b2tfhe_vis_kt = 0.0_R8
      idir = 1 !this indicates the different directions of multidirectional derivative mode
      do ip = 1, nnvar - nsigma_opt
        if (spatial_dep(ip)) then
!         spatially dependent coefficient
          do ii = 1, spatial_points(ip)
            tdatad(idir,2,ii,partype(ip),paris(ip)) = 1.0_R8
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          ! transport coefficients
          case (1)
            parm_dnad(idir,paris(ip)) = 1.0_R8 
          case (2)
            parm_dpad(idir,paris(ip)) = 1.0_R8 
          case (3)
            parm_hcid(idir,paris(ip)) = 1.0_R8
          case (4)
            parm_hced(idir) = 1.0_R8
          case (5)
            parm_vlad(idir,paris(ip)) = 1.0_R8
          case (7)
            parm_vsad(idir,paris(ip)) = 1.0_R8
          case (8)
            parm_sigd(idir) = 1.0_R8
          case (9)
            parm_alfd(idir) = 1.0_R8
          ! boundary conditions
          case (10)
            enepard(idir,parib(ip),1) = 1.0_R8     
          case (11)
            enipard(idir,parib(ip),1) = 1.0_R8      
          case (12)
            conpard(idir,paris(ip),parib(ip),1) = 1.0_R8      
          case (13)
            mompard(idir,paris(ip),parib(ip),1) = 1.0_R8      
          case (14)
            potpard(idir,parib(ip),1) = 1.0_R8      
          case (15) 
            enkpard(idir,parib(ip),1) = 1.0_R8     
          ! recycling
          case (16)
            b2recycd(idir,paris(ip),parib(ip)) = 1.0_R8     
          ! ballooning switches
          case (17)
            switchdiff%b2tqna_ballooning(idir) = 1.0_R8     
          case (18)
            switchdiff%b2tqna_ballooning_rescale(idir) = 1.0_R8     
          ! k-model switches
          case (19)
            switchdiff%keps_cd(idir) = 1.0_R8     
          case (20)
            switchdiff%keps_heat(idir) = 1.0_R8     
          case (21)
            switchdiff%keps_heat_i(idir) = 1.0_R8     
          case (22)
            switchdiff%keps_sig(idir) = 1.0_R8     
          case (23)
            switchdiff%keps_alf(idir) = 1.0_R8     
          case (24)
            switchdiff%keps_visc(idir) = 1.0_R8     
          case (25)
            switchdiff%keps_dkt(idir) = 1.0_R8     
          case (26)
            switchdiff%keps_dzt(idir) = 1.0_R8     
          case (27)
            switchdiff%b2sikt_fac_diss(idir) = 1.0_R8     
          case (28)
            switchdiff%b2sikt_fac_diss_core(idir)  = 1.0_R8     
          case (29)
            switchdiff%b2sikt_fac_sheath(idir) = 1.0_R8     
          case (30)
            switchdiff%b2sikt_fac_sheath_core(idir) = 1.0_R8     
          case (31)
            switchdiff%keps_shear(idir) = 1.0_R8     
          case (32)
            switchdiff%b2tfhi_fsigkt(idir) = 1.0_R8     
          case (33)
            switchdiff%b2sikt_fac_vis_RS(idir) = 1.0_R8     
          case (34)
            switchdiff%b2tfhi_fflokt(idir) = 1.0_R8     
          case (35)
            switchdiff%b2tfhi_fconkt(idir) = 1.0_R8     
          case (36)
            switchdiff%b2tfhi_fflozt(idir) = 1.0_R8     
          case (37)
            switchdiff%b2tfhi_fconzt(idir) = 1.0_R8     
          case (38)
            switchdiff%b2tfhi_fkt_hie(idir) = 1.0_R8     
          case (39)
            switchdiff%b2tfhe_vis_kt(idir) = 1.0_R8     
          case default
            write(*,*) partype(ip)
            call xerrab ('partype out of bounds')
          end select
          idir = idir + 1
        endif
      end do

      if (nsigma_opt.gt.0) then
        sigmad = 0.0_R8
        isigma = npar_opt - nsigma_opt + 1
        do ip = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ip)) then
            sigmad(isigma,ip) = 1.0_R8
            isigma = isigma + 1
          endif
        end do
      endif
      par_opt_physd = 0.0_R8
      do ip = 1, npar_opt
        par_opt_physd(ip,ip) = 1.0_R8
      enddo

      return
      end subroutine set_tgt_perturbation

