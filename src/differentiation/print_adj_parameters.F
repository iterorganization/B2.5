      subroutine print_adj_parameters()
      use b2mod_par_opt_diff
      implicit none
*-----------------------------------------------------------------------
*     print_adj_parameters is used in b2mn_b and b2optim_* to print which
*     parameters are used for sensitivity calculation and optimization
*-----------------------------------------------------------------------
      integer :: icf, ii, ip, ipar
      character*3 ss1
      character*3 ss2
      character*70 key
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
              call xerrab ('partype incompatible with spatial points')
          end select
          do ii = 1, spatial_points(ip)
            write (*,*) '  Calculating gradient of',trim(key)
          end do
        else
!         non spatially dependent --> specified in b2tqna
          select case (partype(ip))
          ! transport coefficients
          case (1)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dna, species ", ss1
          case (2)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_dpa, species ", ss1
          case (3)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hci, species ", ss1
          case (4)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_hce, species ", ss1
          case (6)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vla, species ", ss1
          case (7)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_vsa, species ", ss1
          case (8)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_sig, species ", ss1
          case (9)
            write(ss1,"(I3)") paris(ip)
            write (key, *) "parm_alf, species ", ss1
          ! boundary conditions
          case (10)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enepar, boundary ",ss2
          case (11)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enipar, boundary ",ss2
          case (12)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "conpar, species ",ss1,", boundary ",ss2
          case (13)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "mompar, species ",ss1,", boundary ",ss2
          case (14)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "potpar, boundary ",ss2
          case (15) 
            write(ss2,"(I3)") parib(ip)
            write (key, *) "enkpar, boundary ",ss2
          ! recycling
          case (16)
            write(ss1,"(I3)") paris(ip)
            write(ss2,"(I3)") parib(ip)
            write (key, *) "b2recyc, species ",ss1,", stratum ",ss2
          ! ballooning switches
          case (17)
            write (key, *) "b2tqna_ballooning"
          case (18)
            write (key, *) "b2tqna_ballooning_rescale"
          ! k-model switches
          case (19)
            write (key, *) "keps_cd"
          case (20)
            write (key, *) "keps_heat"
          case (21)
            write (key, *) "keps_heat_i"
          case (22)
            write (key, *) "keps_sig"
          case (23)
            write (key, *) "keps_alf"
          case (24)
            write (key, *) "keps_visc"
          case (25)
            write (key, *) "keps_dkt"
          case (26)
            write (key, *) "keps_dzt"
          case (27)
            write (key, *) "b2sikt_fac_diss"
          case (28)
            write (key, *) "b2sikt_fac_diss_core"
          case (29)
            write (key, *) "b2sikt_fac_sheath"
          case (30)
            write (key, *) "b2sikt_fac_sheath_core"
          case (31)
            write (key, *) "keps_shear"
          case (32)
            write (key, *) "b2tfhi_fsigkt"
          case (33)
            write (key, *) "b2sikt_fac_vis_RS"
          case (34)
            write (key, *) "b2tfhi_fflokt"
          case (35)
            write (key, *) "b2tfhi_fconkt"
          case (36)
            write (key, *) "b2tfhi_fflozt"
          case (37)
            write (key, *) "b2tfhi_fconzt"
          case (38)
            write (key, *) "b2tfhi_fkt_hie"
          case (39)
            write (key, *) "b2tfhe_vis_kt"
          case default
            write(*,*) 'partype = ',partype(ip)
            call xerrab ('partype out of bounds')
          end select
          write (*,*) '  Calculating gradient of',trim(key)
        endif
      end do !nbdirs-nsigma_opt-nmean_opt

      ipar = nnvar - nmean_opt - nsigma_opt ! not subtracting nshift_opt because it is not part of nnvar!

      if (nsigma_opt.gt.0) then
        do ip = 1, nsigma
!         only if sigma is being optimized!
          if (sigma_opt(ip)) then
            ipar = ipar + 1
            if (.not. (partype(ipar).eq.0)) then
              write(*,*) ' ipar, partype(ipar) = ', ipar, partype(ipar)
              write(*,*) ' this was supposed to be a sigma, partype = 0'
              write(ss1,"(I3)") ip
              write(ss2,"(I3)") ipar
              call xerrab('print_adj_parameters: inconsistency '//
     &         'between sigma '//trim(ss1)//'and partype '//trim(ss2))
            endif
            write(ss2,"(I3)") ip
            write (key, *) "sigma",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      if (nmean_opt.gt.0) then
        do ip = 1, nmean
!         only if mean is being optimized!
          if (mean_opt(ip)) then
            ipar = ipar + 1
            if (.not. (partype(ipar).eq.-1)) then
              write(*,*) ' ipar, partype(ipar) = ', ipar, partype(ipar)
              write(*,*) ' this was supposed to be a mean, partype = -1'
              write(ss1,"(I3)") ip
              write(ss2,"(I3)") ipar
              call xerrab('print_adj_parameters: inconsistency '//
     &         'between mean '//trim(ss1)//'and partype '//trim(ss2))
            endif
            write(ss2,"(I3)") ip
            write (key, *) "mean",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      if (nshift_opt.gt.0) then
        do ip = 1, nshift
!         only if shift is being optimized!
          if (shiftopt(ip)) then
            ipar = ipar + 1
            write(ss2,"(I3)") ip
            write (key, *) "shift",ss2
            write (*,*) '  Calculating gradient of',trim(key)
          endif
        end do
      endif

      return
      end subroutine print_adj_parameters

