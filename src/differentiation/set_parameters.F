      subroutine set_parameters(switch)
      use b2mod_types
      use b2mod_par_opt_diff
      use b2mod_transport_namelist_diff
      use b2mod_input_profile_diff
      use b2mod_boundary_namelist_diff
      use b2mod_switches_diff
      use b2mod_neutrals_namelist_diff
      implicit none
      type(switches), intent(inout) :: switch
*-----------------------------------------------------------------------
*     set_parameters is added to the diff'ed version of b2mn_step
*     to assign the optimization parameters obtained from optimization
*     to the input variables actually used by B2.5.
*-----------------------------------------------------------------------
      integer :: idir, ii, ip

      if (.not.flag_optim) return

      idir = 1 !this indicates the different parameters active
      do ip = 1, nnvar - nsigma_opt
        if (spatial_dep(ip)) then
!         spatially dependent coefficient
          if (partype(ip).lt.1 .or. partype(ip).gt.9) then
            write(*,*) partype(ip)
            call xerrab ('partype incompatible with spatial points')
          endif
          do ii = 1, spatial_points(ip)
            tdata(2,ii,partype(ip),paris(ip)) = par_opt_phys(idir)
            idir = idir + 1
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          ! transport coefficients
          case (1) 
            parm_dna(paris(ip)) = par_opt_phys(idir) 
          case (2)
            parm_dpa(paris(ip)) = par_opt_phys(idir) 
          case (3)
            parm_hci(paris(ip)) = par_opt_phys(idir)
          case (4)
            parm_hce = par_opt_phys(idir)
          case (5)
            parm_vla(paris(ip)) = par_opt_phys(idir)
          case (7)
            parm_vsa(paris(ip)) = par_opt_phys(idir)
          case (8)
            parm_sig = par_opt_phys(idir)
          case (9)
            parm_alf = par_opt_phys(idir)
          ! boundary conditions
          case (10)
            enepar(parib(ip),1) = par_opt_phys(idir)
          case (11)
            enipar(parib(ip),1) = par_opt_phys(idir) 
          case (12)
            conpar(paris(ip),parib(ip),1) = par_opt_phys(idir) 
          case (13)
            mompar(paris(ip),parib(ip),1) = par_opt_phys(idir) 
          case (14)
            potpar(parib(ip),1) = par_opt_phys(idir) 
          case (15) 
            enkpar(parib(ip),1) = par_opt_phys(idir)
          ! recycling
          case (16)
            b2recyc(paris(ip),parib(ip)) = par_opt_phys(idir)
          ! ballooning switches
          case (17)
            switch%b2tqna_ballooning = par_opt_phys(idir)
          case (18)
            switch%b2tqna_ballooning_rescale = par_opt_phys(idir)
          ! k-model switches
          case (19)
            switch%keps_cd = par_opt_phys(idir)
          case (20)
            switch%keps_heat = par_opt_phys(idir)
          case (21)
            switch%keps_heat_i = par_opt_phys(idir)
          case (22)
            switch%keps_sig = par_opt_phys(idir)
          case (23)
            switch%keps_alf = par_opt_phys(idir)
          case (24)
            switch%keps_visc = par_opt_phys(idir)
          case (25)
            switch%keps_dkt = par_opt_phys(idir)
          case (26)
            switch%keps_dzt = par_opt_phys(idir)
          case (27)
            switch%b2sikt_fac_diss = par_opt_phys(idir)
          case (28)
            switch%b2sikt_fac_diss_core  = par_opt_phys(idir)
          case (29)
            switch%b2sikt_fac_sheath = par_opt_phys(idir)
          case (30)
            switch%b2sikt_fac_sheath_core = par_opt_phys(idir)
          case (31)
            switch%keps_shear = par_opt_phys(idir)
          case (32)
            switch%b2tfhi_fsigkt = par_opt_phys(idir)
          case (33)
            switch%b2sikt_fac_vis_RS = par_opt_phys(idir)
          case (34)
            switch%b2tfhi_fflokt = par_opt_phys(idir)
          case (35)
            switch%b2tfhi_fconkt = par_opt_phys(idir)
          case (36)
            switch%b2tfhi_fflozt = par_opt_phys(idir)
          case (37)
            switch%b2tfhi_fconzt = par_opt_phys(idir)
          case (38)
            switch%b2tfhi_fkt_hie = par_opt_phys(idir)
          case (39)
            switch%b2tfhe_vis_kt = par_opt_phys(idir)
          case default
            write(*,*) partype(ip)
            call xerrab ('partype out of bounds')
          end select
          idir = idir + 1
        endif
      end do

      return
      end subroutine set_parameters

