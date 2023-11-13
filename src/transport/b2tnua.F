*-----------------------------------------------------------------------
*.specification

      subroutine b2tnua (nCv, ns, switch, pl, dv, rt, nuab, nua)
      use b2mod_types
      use b2mod_math
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_rates
      use b2mod_transport_namelist
      use b2mod_switches
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none

*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      integer is, is2, iCv
      real (kind=R8) :: lnab(nCv), vtha(nCv,0:ns)
      real (kind=R8) :: nuab(nCv,0:ns,0:ns), nua(nCv,0:ns)
c     parameters
c     conversion of constants to CGS units for the collision frequencies
      real (kind=R8) :: erg, cm, gauss, gramm
      parameter (erg=1.0E+7_R8, cm=1.0E+2_R8,
     &          gauss=1.0E+4_R8, gramm=1.0E+3_R8)
      real (kind=R8) :: eV_cgs, e_cgs, c_cgs
      parameter (eV_cgs=ev*erg, e_cgs=qe*c*10.0_R8, c_cgs=c*cm)
      character :: chns*3

*   ..subprogram start-up calls
      call subini ('b2tnua')

c Initialise
      call sfill (nCv*ns, 1.0E+5_R8, nua, 1)
      call sfill (nCv*ns*ns, 1.0E+5_R8, nuab, 1)
      
*   ..calculate thermal velocities
      call thermal_velocity(nCv, ns, pl, switch, vtha)

*   ..calculate collision frequency nuab of species a and b
*   ..ion-ion collisions
      do is=0,ns-1
        if (.not. is_neutral(is)) then
          do is2=0,ns-1
            if (.not. is_neutral(is2)) then
              call b2tlnl(nCv, switch, 1, pl%te, pl%ti, dv%ne, lnab)
              do iCv = 1,nCv
                nuab(iCv,is,is2)=16.0_R8*sqrt(pi)*pl%na(iCv,is2)/
     &             1.0E+6_R8*rt%rza(iCv,is)*rt%rza(iCv,is2)*e_cgs**4/
     &             (3.0_R8*am(is)*mp*gramm)*lnab(iCv)*
     &             (1.0_R8/(am(is)*mp*gramm)+1.0_R8/
     &             (am(is2)*mp*gramm))/(vtha(iCv,is)**2
     &             +vtha(iCv,is2)**2)**1.5_R8
              enddo !! end loop over volumes
            endif
          enddo !! end inner loop over species
        endif
      enddo !! end outer over species

*   ..ion-electron collisions and electron-ion collisions
      do is=0,ns-1
        if (.not. is_neutral(is)) then
          call b2tlnl(nCv, switch, 3, pl%te, pl%ti, dv%ne, lnab)
          do iCv = 1,nCv
            nuab(iCv,is,ns)=16.0_R8*sqrt(pi)*dv%ne(iCv)/
     &        1.0E+6_R8*rt%rza(iCv,is)*
     &        e_cgs**4/(3.0_R8*am(is)*mp*gramm)*lnab(iCv)*
     &        (1.0_R8/(am(is)*mp*gramm)+1.0_R8/(me*gramm))/
     &        (vtha(iCv,is)**2+vtha(iCv,ns-1)**2)**1.5_R8

            nuab(iCv,ns,is)=16.0_R8*sqrt(pi)*pl%na(iCv,is)/
     &        1.0E+6_R8*rt%rza(iCv,is)
     &        *e_cgs**4/(3.0_R8*me*gramm)*lnab(iCv)*
     &        (1.0_R8/(am(is)*mp*gramm)+1.0_R8/(me*gramm))/
     &        (vtha(iCv,is)**2+vtha(iCv,ns-1)**2)**1.5_R8
          enddo !! end loop over volumes
        endif
      enddo !! end loop over species

*   ..electron-electron collisions
      call b2tlnl(nCv, switch, 2, pl%te, pl%ti, dv%ne, lnab)
      do iCv = 1,nCv
        nuab(iCv,ns,ns)=16.0_R8*sqrt(pi)*dv%ne(iCv)/1.0E+6_R8
     &    *e_cgs**4/(3.0_R8*me*gramm)*lnab(iCv)*(2.0_R8/(me*gramm))/
     &    (2.0_R8*vtha(iCv,ns)**2)**1.5_R8
      enddo

*   ..calculate total collision frequency nua of species a
      do is=0,ns
        do is2=0,ns
          do iCv = 1,nCv
            nua(iCv,is)=nua(iCv,is)+nuab(iCv,is,is2)
          enddo !! end loop over volumes
        enddo !! end inner loop over species
      enddo !! end outer loop over species

      if (switch%bohm_iout.eq.1) then
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,nua(:,is),'b2tqna_bohm_nua'//chns)
        enddo
!        write (chns,'(i3.3)') ns electrons
        call my_out_us(70,nCv,0,nua(:,ns),'b2tqna_bohm_nua_electrons')
      endif

*   ..return
      call subend ()
      return
      end subroutine b2tnua


*-----------------------------------------------------------------------
*.specification

      subroutine thermal_velocity (nCv, ns, pl, switch, vtha)
      use b2mod_types
      use b2mod_math
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_rates
      use b2mod_transport_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none

*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (B2Plasma), intent (in) :: pl
      type (switches), intent (in) :: switch
*   ..input/output arguments
      integer is
      real (kind=R8), intent (out) :: vtha(nCv,0:ns)
      real (kind=R8) :: erg, gramm
      parameter (erg=1.0E+7_R8, gramm=1.0E+3_R8)
      character :: chns*3

*   ..subprogram start-up calls
      call subini ('thermal_velocities')

      call sfill (nCv*ns, 0.0_R8, vtha, 1)

      do is=0,ns-1
        vtha(:,is) = sqrt(2.0_R8*erg*pl%ti/(am(is)*mp*gramm)) ! ion + neutral species
      enddo
      vtha(:,ns) = sqrt(2.0_R8*erg*pl%te/(me*gramm)) ! electrons

      if (switch%bohm_iout.eq.1.or.switch%iout_b2wdat.eq.4) then
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,vtha(:,is),'b2tqna_bohm_thvel'//chns)
        enddo
        write (chns,'(i3.3)') ns
        call my_out_us(70,nCv,0,vtha(:,ns),'b2tqna_bohm_thvel'//chns)
      endif

*   ..return
      call subend ()
      return
      end subroutine thermal_velocity
