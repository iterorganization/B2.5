      subroutine set_parameters(switch)
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_transport_namelist_diff
      use b2mod_input_profile_diff
      use b2mod_boundary_namelist_diff
      use b2mod_switches_diff
      implicit none
      type(switches), intent(inout) :: switch
      integer :: idir, ii, ip

      if (.not.flag_optim) return

      idir = 1 !this indicates the different parameters active
      do ip = 1, nnvar - nsigma_opt
        if (spatial_dep(ip)) then
!         spatially dependent coefficient
          do ii = 1, spatial_points(ip)
            tdata(2,ii,partype(ip),paris(ip)) = par_opt_phys(idir)
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          case (1) ! dna
            parm_dna(paris(ip)) = par_opt_phys(idir) 
          case (2) ! dpa
            parm_dpa(paris(ip)) = par_opt_phys(idir) 
          case (3) ! hci
            parm_hci(paris(ip)) = par_opt_phys(idir)
          case (4) ! hce
            parm_hce = par_opt_phys(idir)
          case (5) ! vla
            parm_vla(paris(ip)) = par_opt_phys(idir)
          case (7) ! vsa
            parm_vsa(paris(ip)) = par_opt_phys(idir)
          case (8) ! sig
            parm_sig = par_opt_phys(idir)
          case (9) ! alf
            parm_alf = par_opt_phys(idir)
          case (10) ! enkpar
            enkpar(parib(ip),1) = par_opt_phys(idir)
          case (11) ! b2sikt_fac_sheath
            switch%b2sikt_fac_sheath = par_opt_phys(idir)
          case (12) ! b2sikt_fac_sheath_core
            switch%b2sikt_fac_sheath_core = par_opt_phys(idir)
          case (13) ! keps_cd
            switch%keps_cd = par_opt_phys(idir)
          case (14) ! keps_heat
            switch%keps_heat = par_opt_phys(idir)
          case (15) ! keps_heat_i
            switch%keps_heat_i = par_opt_phys(idir)
          case (16) ! enepar
            enepar(parib(ip),1) = par_opt_phys(idir)
          case (17) ! enip
            enipar(parib(ip),1) = par_opt_phys(idir) 
          case (18) ! conpar
            conpar(paris(ip),parib(ip),1) = par_opt_phys(idir) 
          case (19) ! mompar
            mompar(paris(ip),parib(ip),1) = par_opt_phys(idir) 
          case (20) ! potpar
            potpar(parib(ip),1) = par_opt_phys(idir) 
          case default
            write(*,*) partype(ip)
            call xerrab ('partype out of bounds')
          end select
          idir = idir + 1
        endif
      end do

      return
      end subroutine set_parameters

