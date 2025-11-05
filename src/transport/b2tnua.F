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
      real (kind=R8) :: mred, kk
      character :: chns*3

*   ..subprogram start-up calls
      call subini ('b2tnua')

c Initialise
      call sfill (nCv*(ns+1), 0.0_R8, nua, 1)
      call sfill (nCv*(ns+1)*(ns+1), 0.0_R8, nuab, 1)
      
*   ..calculate thermal velocities
      call thermal_velocity(nCv, ns, pl, switch, vtha)

*   ..calculate collision frequency nuab of species a and b
*   ..ion-electron collisions and electron-ion collisions
      do is=0,ns-1
        if (.not. is_neutral(is)) then
          mred = (am(is)*mp*me)/(am(is)*mp+me)
          call b2tlnl(nCv, switch, 3, pl%te, pl%ti, dv%ne, lnab)
          do iCv = 1,nCv
c            nuab(iCv,ns,is) = (1.0_R8/(4.0_R8*pi)*qe**4*pl%na(iCv,is)/
c     &                        (eps0**2*me*mred)*lnab(iCv))*
c     &                        1.0_R8/(vtha(iCv,ns)**3+1.3_R8*
c     &                        vtha(iCv,is)**3) ! electron-ion Freidberg p.216 formula 9.48

c            nuab(iCv,is,ns) = (me/(am(is)*mp))*nuab(iCv,ns,is) ! ion-electron Freidberg p.217 formula 9.53

            nuab(iCv,is,ns) = pl%na(iCv,is)*(4.8032E-10_R8**4.0_R8)*
     &                        8.0_R8*sqrt(2.0_R8*pi)*lnab(iCv)/(3.0_R8*
     &                        me*1000*am(is)*mp*1000*(pl%te(iCv)/ev*
     &                        1.602E-12_R8/me/1000 + pl%ti(iCv)/ev*
     &                        1.602E-12_R8/am(is)/mp/1000)**1.5_R8) ! TOMATOR1D

           nuab(iCv,ns,ns) = nuab(iCv,ns,ns) + nuab(iCv,is,ns) ! TOMATOR1D
          enddo !! end loop over volumes
        endif
      enddo !! end loop over species

*   ..ion-ion collisions
      do is=0,ns-1
        if (.not. is_neutral(is)) then
          do is2=0,ns-1
            if (.not. is_neutral(is2)) then
              mred = (am(is)*mp*am(is2)*mp)/(am(is)*mp+am(is2)*mp)
              call b2tlnl(nCv, switch, 1, pl%te, pl%ti, dv%ne, lnab)
              do iCv = 1,nCv
                nuab(iCv,is,is2) = pl%na(iCv,is)*4.8032E-10_R8**4*
     &                             8.0_R8*sqrt(2.0_R8*pi)*lnab(iCv)/(
     &                             3.0_R8*(am(is)*mp*1000)**2*(2*
     &                             pl%ti(iCv)/ev*1.602E-12_R8/am(is)/
     &                             mp/1000)**1.5_R8)

c                nuab(iCv,is,is2) = (me/(am(is)*mp))**(1.0_R8/2.0_R8)*
c     &                         (1.0_R8/(4.0_R8*pi)*qe**4*pl%na(iCv,is)/
c     &                         (eps0**2*me*mred)*lnab(iCv))*
c     &                         1.0_R8/(vtha(iCv,ns)**3+1.3_R8*
c     &                         vtha(iCv,is)**3) ! Freidberg
              enddo !! end loop over volumes
            endif
          enddo !! end inner loop over species
        endif
      enddo !! end outer loop over species

*   ..electron-electron collisions
      mred = (me*me)/(me+me)
      call b2tlnl(nCv, switch, 2, pl%te, pl%ti, dv%ne, lnab)
      do iCv = 1,nCv
         nuab(iCv,ns,ns) = 2.91E-6_R8*dv%ne(iCv)*lnab(iCv)
     &                    *pl%te(iCv)**(-3.0_R8/2.0_R8) ! Tom's code

c        nuab(iCv,ns,ns) = (1.0_R8/(4.0_R8*pi)*qe**4*pl%na(iCv,is-1)/
c     &                    (eps0**2*me*mred)*lnab(iCv))*
c     &                    1.0_R8/(vtha(iCv,ns)**3+1.3_R8*
c     &                    vtha(iCv,is)**3) ! Freidberg
      enddo

*   ..neutral collisions ! TODO test
      do is=0,ns-1
        if (is_neutral(is)) then
          do iCv = 1,nCv
! ALL the kk-s are taken from TOMATOR1D
c neutral H with electrons (ionization)
            kk = 2.91E-8_R8*(13.6_R8/pl%te(iCv)**0.39_R8*
     &           exp(-13.6_R8/pl%te(iCv))/(0.232_R8+13.6_R8/pl%te(iCv)))
            nuab(iCv,is,ns) = kk*dv%ne(iCv)
c neutral H (charge exchange)
            kk = 7.829E-9_R8*(pl%ti(iCv) + pl%tn(iCv)/2.0_R8)**0.41_R8
            nuab(iCv,is,is+1) = kk*pl%na(iCv,is+1)
          enddo
        endif
      enddo

*   ..calculate total collision frequency nua of species a
      do is=0,ns-1
        do is2=0,ns-1
          do iCv = 1,nCv
            nua(iCv,is)=nua(iCv,is)+nuab(iCv,is,is2)
          enddo !! end loop over volumes
        enddo !! end inner loop over species
      enddo !! end outer loop over species

      if (switch%bohm_iout.eq.1.or.switch%iout_b2wdat.eq.4) then
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,nua(:,is),'b2tqna_bohm_nua_'//chns)
        enddo
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
      call subini ('thermal_velocity')

      call sfill (nCv*ns, 0.0_R8, vtha, 1)

      do is=0,ns-1
        vtha(:,is) = sqrt(2.0_R8*pl%ti/(am(is)*mp)) ! ion + neutral species
      enddo
      vtha(:,ns) = sqrt(2.0_R8*pl%te/me) ! electrons

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
