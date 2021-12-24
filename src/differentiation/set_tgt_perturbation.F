      subroutine set_tgt_perturbation(swtichdiff)
      use b2mod_types
      use b2mod_par_opt_diffv
      use b2mod_transport_namelist_diffv
      use b2mod_input_profile_diffv
      use b2mod_boundary_namelist_diffv
      use b2mod_switches_diffv
      implicit none
      type(switches), intent(inout) :: switchdiff
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
      switchdiff%b2sikt_fac_sheath = 0.0_R8
      switchdiff%b2sikt_fac_sheath_core = 0.0_R8
      switchdiff%keps_cd = 0.0_R8
      switchdiff%keps_heat = 0.0_R8
      switchdiff%keps_heat_i = 0.0_R8
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
          case (1) ! dna
            parm_dnad(idir,paris(ip)) = 1.0_R8 
          case (2) ! dpa
            parm_dpad(idir,paris(ip)) = 1.0_R8 
          case (3) ! hci
            parm_hcid(idir,paris(ip)) = 1.0_R8
          case (4) ! hce
            parm_hced(idir) = 1.0_R8
          case (5) ! vla
            parm_vlad(idir,paris(ip)) = 1.0_R8
          case (7) ! vsa
            parm_vsad(idir,paris(ip)) = 1.0_R8
          case (8) ! sig
            parm_sigd(idir) = 1.0_R8
          case (9) ! alf
            parm_alfd(idir) = 1.0_R8
          case (10) ! enkpar(1,1)
            enkpard(idir,parib(ip),1) = 1.0_R8
          case (11) ! b2sikt_fac_sheath
            switchdiff%b2sikt_fac_sheath(idir) = 1.0_R8
          case (12) ! b2sikt_fac_sheath_core
            switchdiff%b2sikt_fac_sheath_core(idir) = 1.0_R8
          case (13) ! keps_cd
            switchdiff%keps_cd(idir) = 1.0_R8
          case (14) ! keps_heat
            switchdiff%keps_heat(idir) = 1.0_R8
          case (15) ! keps_heat_i
            switchdiff%keps_heat_i(idir) = 1.0_R8
          case (16) ! enepar
            enepard(idir,parib(ip),1) = 1.0_R8
          case (17) ! enipar
            enipard(idir,parib(ip),1) = 1.0_R8 
          case (18) ! conpar
            conpard(idir,paris(ip),parib(ip),1) = 1.0_R8 
          case (19) ! mompar
            mompard(idir,paris(ip),parib(ip),1) = 1.0_R8 
          case (20) ! potpar
            potpard(idir,parib(ip),1) = 1.0_R8 
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

