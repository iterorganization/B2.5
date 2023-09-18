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
      integer :: idir, ii, ip, ipar
      character*3 ss1
      character*3 ss2
      character*70 key
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
      do ip = 1, nnvar - nsigma_opt - nmean_opt ! not subtracting nshift_opt because it is not part of nnvar!
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
              call xerrab ('partype not compatible with spatial points')
          end select
          do ii = 1, spatial_points(ip)
            tdatad(idir,2,ii,partype(ip),paris(ip)) = 1.0_R8
            idir = idir + 1
            write (*,*) '  Calculating gradient of',trim(key)
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          ! transport coefficients
          case (1)
            parm_dnad(idir,paris(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dna, species ", ss1
          case (2)
            parm_dpad(idir,paris(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dpa, species ", ss1
          case (3)
            parm_hcid(idir,paris(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hci, species ", ss1
          case (4)
            parm_hced(idir) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hce, species ", ss1
          case (6)
            parm_vlad(idir,paris(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vla, species ", ss1
          case (7)
            parm_vsad(idir,paris(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vsa, species ", ss1
          case (8)
            parm_sigd(idir) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_sig, species ", ss1
          case (9)
            parm_alfd(idir) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_alf, species ", ss1
          ! boundary conditions
          case (10)
            enepard(idir,parib(ip),1) = 1.0_R8
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enepar, boundary ",ss2
          case (11)
            enipard(idir,parib(ip),1) = 1.0_R8
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enipar, boundary ",ss2
          case (12)
            conpard(idir,paris(ip),parib(ip),1) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "conpar, species ",ss1,", boundary ",ss2
          case (13)
            mompard(idir,paris(ip),parib(ip),1) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "mompar, species ",ss1,", boundary ",ss2
          case (14)
            potpard(idir,parib(ip),1) = 1.0_R8
            write(ss2,"(I3)") parib(ip)
            write (key, *) "potpar, boundary ",ss2
          case (15) 
            enkpard(idir,parib(ip),1) = 1.0_R8
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enkpar, boundary ",ss2
          ! recycling
          case (16)
            b2recycd(idir,paris(ip),parib(ip)) = 1.0_R8
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "b2recyc, species ",ss1,", stratum ",ss2
          ! ballooning switches
          case (17)
            switchdiff%b2tqna_ballooning(idir) = 1.0_R8
            write (key, *) "b2tqna_ballooning"
          case (18)
            switchdiff%b2tqna_ballooning_rescale(idir) = 1.0_R8
            write (key, *) "b2tqna_ballooning_rescale"
          ! k-model switches
          case (19)
            switchdiff%keps_cd(idir) = 1.0_R8
            write (key, *) "keps_cd"
          case (20)
            switchdiff%keps_heat(idir) = 1.0_R8
            write (key, *) "keps_heat"
          case (21)
            switchdiff%keps_heat_i(idir) = 1.0_R8
            write (key, *) "keps_heat_i"
          case (22)
            switchdiff%keps_sig(idir) = 1.0_R8
            write (key, *) "keps_sig"
          case (23)
            switchdiff%keps_alf(idir) = 1.0_R8
            write (key, *) "keps_alf"
          case (24)
            switchdiff%keps_visc(idir) = 1.0_R8
            write (key, *) "keps_visc"
          case (25)
            switchdiff%keps_dkt(idir) = 1.0_R8
            write (key, *) "keps_dkt"
          case (26)
            switchdiff%keps_dzt(idir) = 1.0_R8
            write (key, *) "keps_dzt"
          case (27)
            switchdiff%b2sikt_fac_diss(idir) = 1.0_R8
            write (key, *) "b2sikt_fac_diss"
          case (28)
            switchdiff%b2sikt_fac_diss_core(idir)  = 1.0_R8
            write (key, *) "b2sikt_fac_diss_core"
          case (29)
            switchdiff%b2sikt_fac_sheath(idir) = 1.0_R8
            write (key, *) "b2sikt_fac_sheath"
          case (30)
            switchdiff%b2sikt_fac_sheath_core(idir) = 1.0_R8
            write (key, *) "b2sikt_fac_sheath_core"
          case (31)
            switchdiff%keps_shear(idir) = 1.0_R8
            write (key, *) "keps_shear"
          case (32)
            switchdiff%b2tfhi_fsigkt(idir) = 1.0_R8
            write (key, *) "b2tfhi_fsigkt"
          case (33)
            switchdiff%b2sikt_fac_vis_RS(idir) = 1.0_R8
            write (key, *) "b2sikt_fac_vis_RS"
          case (34)
            switchdiff%b2tfhi_fflokt(idir) = 1.0_R8
            write (key, *) "b2tfhi_fflokt"
          case (35)
            switchdiff%b2tfhi_fconkt(idir) = 1.0_R8
            write (key, *) "b2tfhi_fconkt"
          case (36)
            switchdiff%b2tfhi_fflozt(idir) = 1.0_R8
            write (key, *) "b2tfhi_fflozt"
          case (37)
            switchdiff%b2tfhi_fconzt(idir) = 1.0_R8
            write (key, *) "b2tfhi_fconzt"
          case (38)
            switchdiff%b2tfhi_fkt_hie(idir) = 1.0_R8
            write (key, *) "b2tfhi_fkt_hie"
          case (39)
            switchdiff%b2tfhe_vis_kt(idir) = 1.0_R8
            write (key, *) "b2tfhe_vis_kt"
          case default
            write(*,*) 'partype = ',partype(ip)
            call xerrab ('partype out of bounds')
          end select
          write (*,*) '  Calculating gradient of',trim(key)
          idir = idir + 1
        endif
      end do

      ipar = nnvar - nmean_opt - nsigma_opt ! not subtracting nshift_opt because it is not part of nnvar!

      if (nsigma_opt.gt.0) then
        sigmad = 0.0_R8
        do ip = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ip)) then
            ipar = ipar + 1
            if (.not. (partype(ipar).eq.0)) then
              write(*,*) ' ipar, partype(ipar) = ', ipar, partype(ipar)
              write(*,*) ' this was supposed to be a sigma, partype = 0'
              write(ss1,"(I3)") ip
              write(ss2,"(I3)") ipar
              call xerrab('set_tgt_perturbation: inconsistency '//
     &         'between sigma '//trim(ss1)//'and partype '//trim(ss2))
            endif
            sigmad(idir,ip) = 1.0_R8
            idir = idir + 1
            write(ss2,"(I3)") ip
            write (key, *) "sigma",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      if (nmean_opt.gt.0) then
        meand = 0.0_R8
        do ip = 1, nmean
!         only if mean is being optimized!
          if (mean_opt(ip)) then
            ipar = ipar + 1
            if (.not. (partype(ipar).eq.-1)) then
              write(*,*) ' ipar, partype(ipar) = ', ipar, partype(ipar)
              write(*,*) ' this was supposed to be a mean, partype = -1'
              write(ss1,"(I3)") ip
              write(ss2,"(I3)") ipar
              call xerrab('set_tgt_perturbation: inconsistency '//
     &         'between mean '//trim(ss1)//'and partype '//trim(ss2))
            endif
            meand(idir,ip) = 1.0_R8
            idir = idir + 1
            write(ss2,"(I3)") ip
            write (key, *) "mean",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      if (nshift_opt.gt.0) then
        shiftd = 0.0_R8
        do ip = 1, nshift
!         only if shift is being optimized!
          if (shiftopt(ip)) then
            shiftd(idir,ip) = 1.0_R8
            idir = idir + 1
            write(ss2,"(I3)") ip
            write (key, *) "shift",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      par_opt_physd = 0.0_R8
      do ip = 1, npar_opt
        par_opt_physd(ip,ip) = 1.0_R8
      enddo

      return
      end subroutine set_tgt_perturbation

