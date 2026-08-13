*-----------------------------------------------------------------------
*.specification

      subroutine integrate_elastic_collisions (mt, mb, rtt, nt, 
     .    rttr, ntr, rt_mcoll, n_mcoll, coeff, 
     .            N_int, rate_mom, rate_E)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nt, ntr, n_mcoll, N_int
      real (kind=R8), intent(in) ::
     *  rtt(0:nt), rttr(0:ntr), rt_mcoll(0:n_mcoll),
     *  coeff(0:8,0:8), mt, mb
      real (kind=R8), intent(out) ::
     *  rate_mom(0:nt,0:ntr,0:n_mcoll), rate_E(0:nt,0:ntr)
*-----------------------------------------------------------------------
*.documentation

*
*  1. purpose
*
*     INTEGRATE_ELASTIC_COLLISIONS creates rate coefficients for momentum
*     and energy transfer due to elastic collisions between molecules and ions, 
*     consitently with how EIRENE treats them. 
*
*
*  3. description
*
*     rtt represents temperature or temperature per amu; units are eV.
*
*  4. references
*
*  5. parameters
*
*     rtt - (0:nt) real*8 array, input.
*     rtt(0:nt) specifies the sequence of values of the temperature
*     variable, expressed in eV, for which the rate coefficients are
*     evaluated. rtt(0:nt) must be positive.
*
*     nt - integer, input.
*     nt specifies the upper dimension bound for the array rtt.
*     It will hold that 0.le.nt.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind = R8) ::P(1:9), Tt, Tb, dV_coll,
     .     VXN, VYN, VZN, ZRG, VXI, VYI, VZI, VXDR, VYDR, VZDR,
     .     VXNp, VYNp, VZNp, VXIp, VYIp, VZIp, !! post-collisional velocities
     .     delta_mom_x, delta_mom_x_sum, delta_E, delta_E_sum,
     .     AMU, CVEL2A, FG1, FG2, FG3, FCRG, cou, poly(1:9),
     .     rc1min, rc1max, rc2min, rc2max, vth, 
     .     ADDT, ADDTL, TII, CVELI2, TMASS,
     .     PVELQ, SIGVEL, EEFEL, RCMIN, RCMAX, ELB, EMINL, FP(1:6),
     .     b2_sngl_poly, CVELAA, mr,  rate_mom0, rate_E0, denom_E,
     .     delta_E_out(0:nt,0:ntr,-1:n_mcoll), dE, fric_tot, ratio, 
     .     delta_mom_x_out(0:nt,0:ntr,0:n_mcoll), dE_tot, dE_T,
     .     ratio_out(0:nt,0:ntr,0:n_mcoll), wrk, rel_err, fric,
     .     rel_err_out(0:nt,0:ntr,0:n_mcoll), dE_num, dE_reconstr
      integer :: i_n, i, j, k, neg_rate_cnt_m, neg_rate_cnt_E, N_mult,
     .           N_final

*   ..initialisation
      data P /2.70000000E+00_R8, 3.00000000E+00_R8,1.00000000E+00_R8, !! Reaction H0 0.3T p + H2 ! p + H2 potential, binary elastic TODO read from AMJUEL
     .        2.83550000E+00_R8, 2.18038000E+00_R8,3.49068700E+00_R8, !! Morse potential 
     .        0.00000000E+00_R8,-2.02500000E+00_R8,0.00000000E+00_R8/

      external xerrab

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('integrate_elastic_collisions')
      call xertst (0.le.nt, 'faulty argument nt')


#if ( defined(_OPENMP) )
!$omp parallel default(none)
cc shared variables
!$omp& shared(amu, cvel2a, cvelaa, cveli2, tmass, eminl, mr)
!$omp& shared(rc1min, rc1max, rc2min, rc2max, addt, addtl)
!$omp& shared(nt, ntr, n_mcoll, mt, mb)
!$omp& shared(rtt, rttr, rt_mcoll, coeff, P, N_int)
!$omp& shared(rcmin, rcmax, fp)
!$omp& shared(rate_mom, rate_E, neg_rate_cnt_m, neg_rate_cnt_E)
cc private variables
!$omp& private(j, k, Tb, Tt, vth, dv_coll, tii)
!$omp& private(denom_E, N_mult, N_final)
!$omp& private(delta_mom_x, delta_mom_x_sum, delta_E, delta_E_sum)
!$omp& private(FCRG, ZRG, FG1, FG2, FG3, VXDR, VYDR, VZDR)
!$omp& private(VXN, VYN, VZN, VXI, VYI, VZI)
!$omp& private(VXNp, VYNp, VZNp, VXIp, VYIp, VZIp)
!$omp& private(EEFEL, ELB, SIGVEL, PVELQ, cou, poly, rate_mom0, rate_E0)
!$omp& private(delta_E_out, delta_mom_x_out)

#endif

      rate_mom = 0.0_R8

cc set constants
      AMU=1.6606D-27 !!    ATOMIC MASS UNIT       [kg]
      CVEL2A=SQRT(10000.0_R8*ev/AMU)           !! [cm/s] (a constant to translate from CGS to SI)
      CVELAA = sqrt(2.0_R8)*CVEL2A             !! [cm/s] (a constant to translate from CGS to SI)
      CVELI2 = 1.0_R8/CVELAA**2                !! [s^2/cm^2]
      TMASS = 2.0146000000000002_R8            !! [-]
      EMINL = -2.3_R8                          !! hard-coded cut-off for EBEAM parameter in H.3 fits
      mr = (mb*mt)/(mb+mt)                     !! reduced mass [kg]
      rc1min = -20.0_R8                        !! database cutoffs TODO generalize? -> also these values in EIRENE
      rc1max =  20.0_R8  
      rc2min = -20.0_R8  
      rc2max =  100.0_R8
      ADDT=  0.5_R8                            !! PMASS/RMASSP(IPL), TODO generalize
      ADDTL=LOG(ADDT)
      FP = 0.0_R8                              !! ?? TODO
      RCMIN = -HUGE(1.0_R8)                    !! database cutoffs TODO generalize?
      RCMAX = HUGE(1.0_R8)

      neg_rate_cnt_m = 0                       !! initialize negative rate counter
      neg_rate_cnt_E = 0                       !! initialize negative rate counter


cc FIRST DO ONLY MOMENTUM: 3D table
#if (defined(_OPENMP) )
!$omp do 
#endif
      do i = 0, nt
        Tb = rtt(i)                          !! Background temperature (ion) [ev]
        !! pre-evaluate the double polynomial for the ion temperature
        TII = log(Tb)+ADDTL                  !! [log(ev)], ADDTL = due to use of D while tables are for H
        call b2_dbl_poly(coeff, TII, 0.0_R8, cou, poly, 
     .             rc1min, rc1max, rc2min, rc2max)
        do j = 0, ntr
          Tt = Tb*rttr(j)                        !! 'Target' temperature (neutral molecule) [ev]
          vth = sqrt(2.0_R8*Tb*ev/(mp*2.0_R8))   !! Thermal velocity TODO generalize species, mass   [m/s]
          do k = -1, n_mcoll                     !! k = -1 -> energy term, k = 0:n_mcoll -> momentum
            if (k.eq.-1) then  
              delta_E_sum = 0.0_R8    
              denom_E = 1.5_R8*ev*(Tt-Tb)      
              dV_coll = 0.0_R8           
              N_mult = 10                        !! take 10x more particles for energy
            else
              delta_E_sum = 0.0_R8    
              delta_mom_x_sum = 0.0_R8       
              dV_coll = rt_mcoll(k)*vth          !! velocity difference in the direction of the momentum transfer [m/s]
              N_mult = 1
            endif 
cc calculate contributions with Monte Carlo sampling
            N_final = N_int*N_mult
            do i_n = 1,N_final
              !! sample a molecule particle velocity
              call b2_fgauss(FG1,FG2,FG3)
              !! drift velocities (DR)
              VXDR = dV_coll*100.0_R8             !! [m/s -> cm/s]
              VYDR = 0.0_R8                       
              VZDR = 0.0_R8
              !! sample random velocity from distribution
              FCRG=CVEL2A/SQRT(TMASS*2.0_R8) !! TODO generalize masses
              ZRG = FCRG*SQRT(Tt)
              !! neutral molecule particle velocity
              VXN=FG1*ZRG+VXDR
              VYN=FG2*ZRG+VYDR
              VZN=FG3*ZRG+VZDR

              !! sample an ion particle velocity (w/o bulk velocity, this is captured in the neutral)
              call b2_fgauss(FG1,FG2,FG3)     !! units_TODO (no dimension?)
              FCRG=CVEL2A/SQRT(TMASS)         !! [cm/s/sqrt(kg)]
              ZRG = FCRG*SQRT(Tb)             !! units_TODO ?? -> Tb = joule = kg m^2/s^2
              VXI=FG1*ZRG
              VYI=FG2*ZRG
              VZI=FG3*ZRG

              !! evaluate for the particle velocity, calculate SIGMA (carefull, density not included here, as opposed to EIRENE)
              PVELQ = VXN**2+VYN**2+VZN**2 !! Square of the velocity difference between the kinetic D2 particle and the ion bulk velocity. 
                                             !! Because the bulk velocity difference is fully assigned to the neutral here above, 
                                             !! this is the square of the neutral particle velocity
                                             !! [cm^2/s^2]
              EEFEL = LOG(CVELI2*TMASS)    !! [log(s^2/cm^2)]
              ELB = MAX(EMINL,LOG(PVELQ)+EEFEL) !! units_TODO no dimension?
              SIGVEL = b2_sngl_poly(poly,ELB,RCMIN,RCMAX,FP,0,0) !! Caution with name conventions and units: in EIRENE, this quantity
                                                                 !! contains the ion density as well. Here, it is a reaction rate coefficient [cm^3/s]
              SIGVEL = SIGVEL/1.0e6 !! convert [cm^3/s -> m^3/s]

              !! execute the elastic collision event
              call execute_elastic_collision (P, mt/mp, mb/mp, cveli2, 
     .             VXN,  VYN,  VZN,  VXI,  VYI,  VZI,
     .              VXNp,  VYNp,  VZNp,  VXIp,  VYIp,  VZIp) !! both input and output is in cm/s

              if (k.eq.-1) then
                !! calculate energy transferred
                delta_E = 0.5_R8*mt*( (VXN**2+VYN**2+VZN**2)-
     .                (VXNp**2+VYNp**2+VZNp**2))*SIGVEL/1.0e4_R8
                delta_E_sum = delta_E_sum + delta_E
              else
                !! calculate momentum transferred
                delta_mom_x = mt*(VXN-VXNp)*SIGVEL/100.0_R8 !! x = direction of interest + convert back to m/s
                delta_mom_x_sum = delta_mom_x_sum + delta_mom_x
                !! calculate energy transferred
                delta_E = 0.5_R8*mt*( (VXN**2+VYN**2+VZN**2)-
     .                (VXNp**2+VYNp**2+VZNp**2))*SIGVEL/1.0e4_R8
                delta_E_sum = delta_E_sum + delta_E
              endif
            enddo !! end Monte Carlo loop
cc compute final data entries
            if (k.eq. -1) then
              !! energy
              delta_E = delta_E_sum / N_final          !! average over Monte Carlo samples
              delta_E_out(i,j,k) = delta_E
              rate_E0 =delta_E/denom_E
              if (rate_E0.gt.0) then
                rate_E(i,j) = rate_E0                  !! save result
              else
                !! prevent negative rate
                write(*,*) 'negative E rate obtained !!'
                rate_E(i,j) = 1.0e-20_R8  
                neg_rate_cnt_E = neg_rate_cnt_E + 1
              endif
            else
              !! momentum
              delta_mom_x = delta_mom_x_sum / N_final  !! average over Monte Carlo samples
              delta_mom_x_out(i,j,k) = delta_mom_x
              rate_mom0 = delta_mom_x/mr/dV_coll
              if (rate_mom0.gt.0) then
                rate_mom(i,j,k) = rate_mom0             !! save result
              else
                !! prevent negative rate
                write(*,*) 'negative mom rate obtained !!'
                rate_mom(i,j,k) = 1.0e-20_R8  
                neg_rate_cnt_m = neg_rate_cnt_m + 1
              endif
              delta_E = delta_E_sum / N_final          !! average over Monte Carlo samples
              delta_E_out(i,j,k) = delta_E
            endif
cc end loops over table dimensions   
          enddo 
        enddo
      enddo
#if ( defined(_OPENMP) )
!$omp end do 
!$omp end parallel
#endif

      write(*,*) 'negative mom rate obtained for'
      write(*,*) neg_rate_cnt_m, 'out of'
      write(*,*) (nt+1)*(ntr+1)*(n_mcoll+1), 'entries'

      write(*,*) 'negative E rate obtained for'
      write(*,*) neg_rate_cnt_E, 'out of'
      write(*,*) (nt+1)*(ntr+1), 'entries'

!! put both in single loop, document everything
      do i = 0, nt
        Tb = rtt(i)                          !! Background temperature (ion) [ev]
        vth = sqrt(2.0_R8*Tb*ev/(mp*2.0_R8))   !! Thermal velocity TODO generalize species, mass   [m/s]
        do j = 0, ntr
          Tt = Tb*rttr(j)                        !! 'Target' temperature (neutral molecule) [ev]
          do k = 0, n_mcoll                     !! k = -1 -> energy term, k = 0:n_mcoll -> momentum
            write(*,*) 'i, j, k', i, j, k
            dV_coll = rt_mcoll(k)*vth          !! velocity difference in the direction of the momentum transfer [m/s]
            dE_tot = -delta_E_out(i,j,k)        !! total energy gain molecules
            write(*,*) 'dE_tot', dE_tot
            dE_T = -delta_E_out(i,j,-1)          !! energy gain molecules without bulk velocity
            write(*,*) 'dE_T', dE_T
            wrk = -dV_coll*delta_mom_x_out(i,j,k) !! drag work on molecules
            write(*,*) 'wrk', wrk 
            dE = dE_tot-dE_T-wrk
            write(*,*) 'diff dE', dE
            fric_tot = abs(wrk)
            write(*,*) 'fric_tot', fric_tot 
            ratio = dE/fric_tot 
            write(*,*) 'ratio', ratio
            ratio_out(i,j,k) = ratio
            dE_reconstr = dE_T + wrk + fric_tot/3.0_R8
            write(*,*) 'dE reconstr', dE_reconstr
            fric = wrk + fric_tot/3.0_R8
            rel_err = abs(dE_tot-dE_reconstr)/
     .       max(max(abs(dE_tot),abs(fric)),
     .           max(abs(dE_reconstr),abs(dE_T))) 
            rel_err_out(i,j,k) = rel_err 
            write(*,*) 'rel error', rel_err
          enddo
        enddo
      enddo


      call my_out_us(70,(nt+1)*(ntr+1)*(n_mcoll+1),0,delta_mom_x_out,
     .                                            'delta_mom_x_out')
      call my_out_us(70,(nt+1)*(ntr+1)*(n_mcoll+2),0,delta_E_out,
     .                                            'delta_E_out')
      call my_out_us(70,(nt+1)*(ntr+1)*(n_mcoll+1),0,ratio_out,
     .                                            'ratio_out')
      call my_out_us(70,(nt+1)*(ntr+1)*(n_mcoll+1),0,rel_err_out,
     .                                            'rel_err')


* ..return
      call subend ()
      return

*----------------------------------------------------------------------
*.end integrate_elastic_collisions

      end subroutine integrate_elastic_collisions


cc this routine is only for debugging: calculate the elastic rate coefficients for one specific cell
      subroutine integrate_elastic_collisions_cell (mt, mb, 
     *           tm, ti, ui_p, ui_r, ui_d, un_p, un_r, un_d)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      use b2mod_eirreac
      implicit none
*   ..input arguments (unchanged on exit)
      real (kind=R8), intent(in) ::
     *  mt, mb, tm, ti, ui_p, ui_r, ui_d,
     *  un_p, un_r, un_d
      integer, save :: N_int = 100000000
      integer :: i_n, I, J

*   ..local variables
      real (kind = R8) ::P(1:9), Tt, Tb, dV_coll,
     .                   VXN, VYN, VZN, ZRG, VXI, VYI, VZI,
     .                   VXDR, VYDR, VZDR,
     .                   VXNp, VYNp, VZNp, VXIp, VYIp, VZIp, !! post-collisional velocities
     .                   delta_mom_x, delta_mom_x_sum, 
     .                   delta_mom_y, delta_mom_y_sum, 
     .                   delta_mom_z, delta_mom_z_sum, 
     .                   delta_E, delta_E_sum, AMU, CVEL2A,
     .                   FG1, FG2, FG3, FCRG, cou, poly(1:9),
     .                   rc1min, rc1max, rc2min, rc2max, vth, 
     .                   ADDT, ADDTL, TII, CVELI2, TMASS,
     .                   PVELQ, SIGVEL, EEFEL, RCMIN, RCMAX, ELB,
     .                   EMINL, FP(1:6), b2_sngl_poly, CVELAA, mr,
     .                   rate_mom_x, rate_mom_y, rate_mom_z, 
     .                   rate_E0, denom_E, denom_E0, Qfr, Qfri,
     .                   denom_E1, rate_E0nb, rate_E0s, coeff(0:8,0:8),
     .                   delta_E_noV, wrki, wrkn, Vcom_x, Vcom_y, 
     .                   Vcom_z, Qt

*   ..initialisation
      data P /2.70000000E+00_R8, 3.00000000E+00_R8,1.00000000E+00_R8, !! Reaction H0 0.3T p + H2 ! p + H2 potential, binary elastic TODO read from AMJUEL
     .        2.83550000E+00_R8, 2.18038000E+00_R8,3.49068700E+00_R8, !! Morse potential 
     .        0.00000000E+00_R8,-2.02500000E+00_R8,0.00000000E+00_R8/

      external xerrab


* ..preliminaries
      call subini ('integrate_elastic_collisions_cell')

cc get data
      CALL SLREAC_b2plot(1,'AMJUEL  ','H.3 ','0.3T      ','EL ')
      do J=1,9
       do I=1,9
         coeff(J-1,I-1)=CREAC(J,I,1)
       enddo
      enddo
     
cc set constants
      AMU=1.6606D-27 !!    ATOMIC MASS UNIT       [kg]
      CVEL2A=SQRT(10000.0_R8*ev/AMU)           !! [cm/s] (a constant to translate from CGS to SI)
      CVELAA = sqrt(2.0_R8)*CVEL2A             !! [cm/s] (a constant to translate from CGS to SI)
      CVELI2 = 1.0_R8/CVELAA**2                !! [s^2/cm^2]
      TMASS = 2.0146000000000002_R8            !! [-]
      EMINL = -2.3_R8                          !! hard-coded cut-off for EBEAM parameter in H.3 fits
      mr = (mb*mt)/(mb+mt)                     !! reduced mass [kg]
      rc1min = -20.0_R8                        !! database cutoffs TODO generalize? -> also these values in EIRENE
      rc1max =  20.0_R8  
      rc2min = -20.0_R8  
      rc2max =  100.0_R8
      ADDT=  0.5_R8                            !! PMASS/RMASSP(IPL), TODO generalize
      ADDTL=LOG(ADDT)
      FP = 0.0_R8                              !! ?? TODO
      RCMIN = -HUGE(1.0_R8)                    !! database cutoffs TODO generalize?
      RCMAX = HUGE(1.0_R8)

      Tb = ti                           !! Background temperature (ion) [ev]
      !! pre-evaluate the double polynomial for the ion temperature
      TII = log(Tb)+ADDTL                  !! [log(ev)], ADDTL = due to use of D while tables are for H
      call b2_dbl_poly(coeff, TII, 0.0_R8, cou, poly, 
     .             rc1min, rc1max, rc2min, rc2max)
      Tt = tm                        !! 'Target' temperature (neutral molecule) [ev]


cc 1) WITHOUT VELOCITY DIFFERENCE
cc calculate contributions with Monte Carlo sampling          
      delta_E_sum = 0.0_R8    
      denom_E0 = 1.5_R8*ev*(Tt-Tb)
      denom_E1 = 0.0_R8 
      denom_E = denom_E0 + denom_E1
cc take more Monte Carlo samples if denominator is small
      do i_n = 1,N_int
        !! sample a molecule particle velocity
         call b2_fgauss(FG1,FG2,FG3)
         !! sample random velocity from distribution
         FCRG=CVEL2A/SQRT(TMASS*2.0_R8) !! TODO generalize masses
         ZRG = FCRG*SQRT(Tt)
        !! neutral molecule particle velocity
         VXN=FG1*ZRG
         VYN=FG2*ZRG
         VZN=FG3*ZRG
         !! sample an ion particle velocity 
         call b2_fgauss(FG1,FG2,FG3)     !! units_TODO (no dimension?)
         FCRG=CVEL2A/SQRT(TMASS)         !! [cm/s/sqrt(kg)]
         ZRG = FCRG*SQRT(Tb)             !! units_TODO ?? -> Tb = joule = kg m^2/s^2
         VXI=FG1*ZRG
         VYI=FG2*ZRG
         VZI=FG3*ZRG
         !! evaluate for the particle velocity, calculate SIGMA (carefull, density not included here, as opposed to EIRENE)
         PVELQ = VXN**2+VYN**2+VZN**2 !! Square of the velocity difference between the kinetic D2 particle and the ion bulk velocity. 
                                             !! Because the bulk velocity difference is fully assigned to the neutral here above, 
                                             !! this is the square of the neutral particle velocity
                                             !! [cm^2/s^2]
         EEFEL = LOG(CVELI2*TMASS)    !! [log(s^2/cm^2)]
         ELB = MAX(EMINL,LOG(PVELQ)+EEFEL) !! units_TODO no dimension?
         SIGVEL = b2_sngl_poly(poly,ELB,RCMIN,RCMAX,FP,0,0) !! Caution with name conventions and units: in EIRENE, this quantity
                                                            !! contains the ion density as well. Here, it is a reaction rate coefficient [cm^3/s]
         SIGVEL = SIGVEL/1.0e6 !! convert [cm^3/s -> m^3/s]
         !! execute the elastic collision event
         call execute_elastic_collision (P, mt/mp, mb/mp, cveli2, 
     .             VXN,  VYN,  VZN,  VXI,  VYI,  VZI,
     .              VXNp,  VYNp,  VZNp,  VXIp,  VYIp,  VZIp) !! both input and output is in cm/s
         !! calculate momentum and energy transferred
         delta_E = 0.5_R8*mt*( (VXN**2+VYN**2+VZN**2)-
     .                (VXNp**2+VYNp**2+VZNp**2))*SIGVEL/1.0e4_R8
         delta_E_sum = delta_E_sum + delta_E
      enddo !! end Monte Carlo loop
       !! average over Monte Carlo samples
      delta_E = delta_E_sum / N_int
      delta_E_noV = delta_E
cc compute final data entries
      write(*,*) 'Energy transfer rate no dV', delta_E  !! first check if this one is correct to begin with ...

cc 2) NOW WITH VELOCITY DIFFERENCE
cc calculate contributions with Monte Carlo sampling
      delta_mom_x_sum = 0.0_R8   
      delta_mom_y_sum = 0.0_R8       
      delta_mom_z_sum = 0.0_R8           
      delta_E_sum = 0.0_R8    
      denom_E0 = 1.5_R8*ev*(Tt-Tb)
      denom_E1 = 0.0_R8 
      denom_E = denom_E0 + denom_E1
cc take more Monte Carlo samples if denominator is small
      do i_n = 1,N_int
        !! sample a molecule particle velocity
         call b2_fgauss(FG1,FG2,FG3)
         !! drift velocities (DR)
         VXDR = un_p*100.0_R8             !! [m/s -> cm/s]
         VYDR = un_r*100.0_R8                       
         VZDR = un_d*100.0_R8
         !! sample random velocity from distribution
         FCRG=CVEL2A/SQRT(TMASS*2.0_R8) !! TODO generalize masses
         ZRG = FCRG*SQRT(Tt)
        !! neutral molecule particle velocity
         VXN=FG1*ZRG+VXDR
         VYN=FG2*ZRG+VYDR
         VZN=FG3*ZRG+VZDR
         !! sample an ion particle velocity 
         call b2_fgauss(FG1,FG2,FG3)     !! units_TODO (no dimension?)
         VXDR = ui_p*100.0_R8             !! [m/s -> cm/s]
         VYDR = ui_r*100.0_R8                       
         VZDR = ui_d*100.0_R8
         FCRG=CVEL2A/SQRT(TMASS)         !! [cm/s/sqrt(kg)]
         ZRG = FCRG*SQRT(Tb)             !! units_TODO ?? -> Tb = joule = kg m^2/s^2
         VXI=FG1*ZRG+VXDR
         VYI=FG2*ZRG+VYDR
         VZI=FG3*ZRG+VZDR
         !! evaluate for the particle velocity, calculate SIGMA (carefull, density not included here, as opposed to EIRENE)
         PVELQ = VXN**2+VYN**2+VZN**2 !! Square of the velocity difference between the kinetic D2 particle and the ion bulk velocity. 
                                             !! Because the bulk velocity difference is fully assigned to the neutral here above, 
                                             !! this is the square of the neutral particle velocity
                                             !! [cm^2/s^2]
         EEFEL = LOG(CVELI2*TMASS)    !! [log(s^2/cm^2)]
         ELB = MAX(EMINL,LOG(PVELQ)+EEFEL) !! units_TODO no dimension?
         SIGVEL = b2_sngl_poly(poly,ELB,RCMIN,RCMAX,FP,0,0) !! Caution with name conventions and units: in EIRENE, this quantity
                                                            !! contains the ion density as well. Here, it is a reaction rate coefficient [cm^3/s]
         SIGVEL = SIGVEL/1.0e6 !! convert [cm^3/s -> m^3/s]
         !! execute the elastic collision event
         call execute_elastic_collision (P, mt/mp, mb/mp, cveli2, 
     .             VXN,  VYN,  VZN,  VXI,  VYI,  VZI,
     .              VXNp,  VYNp,  VZNp,  VXIp,  VYIp,  VZIp) !! both input and output is in cm/s
         !! calculate momentum and energy transferred
         delta_mom_x = mt*(VXN-VXNp)*SIGVEL/100.0_R8 !! x = direction of interest + convert back to m/s
         delta_mom_y = mt*(VYN-VYNp)*SIGVEL/100.0_R8 
         delta_mom_z = mt*(VZN-VZNp)*SIGVEL/100.0_R8 
         delta_E = 0.5_R8*mt*( (VXN**2+VYN**2+VZN**2)-
     .                (VXNp**2+VYNp**2+VZNp**2))*SIGVEL/1.0e4_R8
         delta_mom_x_sum = delta_mom_x_sum + delta_mom_x
         delta_mom_y_sum = delta_mom_y_sum + delta_mom_y
         delta_mom_z_sum = delta_mom_z_sum + delta_mom_z
         delta_E_sum = delta_E_sum + delta_E
       enddo !! end Monte Carlo loop
       !! average over Monte Carlo samples
       delta_mom_x = delta_mom_x_sum / N_int
       delta_mom_y = delta_mom_y_sum / N_int
       delta_mom_z = delta_mom_z_sum / N_int
       delta_E = delta_E_sum / N_int
cc compute final data entries
       !! momentum
       rate_mom_x = delta_mom_x/mr/(un_p-ui_p)
       rate_mom_y = delta_mom_y/mr/(un_r-ui_r)
       rate_mom_z = delta_mom_z/mr/(un_d-ui_d)
       Vcom_x = (mt*un_p+mb*ui_p)/(mt+mb)
       Vcom_y = (mt*un_r+mb*ui_r)/(mt+mb)
       Vcom_z = (mt*un_d+mb*ui_d)/(mt+mb)
       write(*,*) 'momentum rate x', delta_mom_x !! = force on ions
       write(*,*) 'momentum rate y', delta_mom_y
       write(*,*) 'momentum rate z', delta_mom_z
       write(*,*) 'momentum rate coeff x', rate_mom_x
       write(*,*) 'momentum rate coeff y', rate_mom_y
       write(*,*) 'momentum rate coeff z', rate_mom_z
       write(*,*) 'Energy transfer rate', delta_E  !! first check if this one is correct to begin with ...
cc write input quantities
       write(*,*) 'Ti', Tb
       write(*,*) 'Tm', Tt
       write(*,*) 'ui_p', ui_p
       write(*,*) 'un_p', un_p
       write(*,*) 'un_r', un_r
       write(*,*) 'un_d', un_d
       wrki = delta_mom_x*ui_p
       write(*,*) 'mech. work ions', wrki
       wrkn = -delta_mom_x*un_p-delta_mom_y*un_r-delta_mom_z*un_d
       write(*,*) 'mech. work mols', wrkn
       write(*,*) 'mech. work tot', wrkn + wrki
       Qfr = delta_mom_x*(un_p-ui_p) + 
     .       delta_mom_y*(un_r-ui_r) + 
     .       delta_mom_z*(un_d-ui_d) 
       Qfri = Qfr*mt/(mb+mt)
       Qt = Vcom_x*delta_mom_x + Vcom_y*delta_mom_y + Vcom_z*delta_mom_z
       write(*,*) 'total friction heating', Qfr
       write(*,*) 'ion friction heating', Qfri
       write(*,*) 'Energy transfer T + Kin', Qfri + wrki + delta_E_noV
       write(*,*) 'Energy transfer T + Qt', delta_E_noV + Qt



* ..return
      call subend ()
      return

*----------------------------------------------------------------------
*.end integrate_elastic_collisions_cell

      end subroutine integrate_elastic_collisions_cell




      subroutine execute_elastic_collision (P, mt, mb, cveli2,
     .         vxt_in,  vyt_in,  vzt_in,  vxb_in,  vyb_in,  vzb_in,
     .         vxt_out, vyt_out, vzt_out, vxb_out, vyb_out, vzb_out)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      !! t = target particle, b = background particle, in = before collision, out = after collision
      real (kind=R8), intent(in) ::
     *  P(1:9), mt, mb, vxt_in, vyt_in, vzt_in, vxb_in, vyb_in, vzb_in,
     *  cveli2
      real (kind=R8), intent(out) ::
     *  vxt_out, vyt_out, vzt_out, vxb_out, vyb_out, vzb_out


*-----------------------------------------------------------------------
*.documentation
*
*     Execute a single elastic collision event. This function mimmicks 
*     what happens in EIRMOD_VELOEL
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind = R8) :: RMN, RMI, RMSI, RLMS,
     . VX, VY, VZ, VXI, VYI, VZI, VRELX, VRELY, VRELZ,
     . VRQYZ, VRELQ, VREL, VRYZ, VSX, VSY, VSZ,
     . PMASS, DEFEL, ELAB, EPS60, CEL, PIAI, B, BMAX, rb,
     . RS, ER, RES, CHI, CCHI, PH, CPH, SPH, EPS, CEPS, SEPS, 
     . VRSX, VRSY, VRSZ, VELX, VELY, VELZ, VELQ, VEL, randx,
     . B2_RSTERN, b2_cross
      integer :: IFLAG

*   ..initialisation
      data EPS60/1E-60_R8/
      external xerrab

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('execute_elastic_collision')

cc set masses
      RMN  = mt !! mass D2
      RMI = mb !! mass D+
      RMSI = 1.0_R8/(RMN+RMI)
      RLMS = RMN*RMI*RMSI !! reduced mass = mr in Kotov
cc set velocities (switch to EIRENE naming convention)
      VX = vxt_in 
      VY = vyt_in
      VZ = vzt_in
      VXI = vxb_in
      VYI = vyb_in
      VZI = vzb_in
cc relative velocities between ion and neutral are calculated
      VRELX=VX-VXI
      VRELY=VY-VYI
      VRELZ=VZ-VZI
cc additional quantities related to the relative velocities are calculated
      VRQYZ = VRELY**2+VRELZ**2
      VRELQ=VRELX**2+VRQYZ !! kinetic energy associated with relative velocity
      VREL=SQRT(VRELQ) !! = vr in Kotov
      VRYZ=SQRT(VRQYZ+EPS60)
cc TODO what is the purpose of VRQYZ? How do Y and Z differ from X?
cc I assume this has to do with some arbitrary reference frame?

cc calculate center-of-mass velocities (= v_c in Kotov (3.23))
      VSX=(RMI*VXI+RMN*VX)*RMSI
      VSY=(RMI*VYI+RMN*VY)*RMSI
      VSZ=(RMI*VZI+RMN*VZ)*RMSI

cc compute ELAB (Here, by abuse of notation, ELAB is already log(ELAB))
      PMASS = 1.0073000000000001               !! [-]
      DEFEL=LOG(CVELI2*PMASS)
      ELAB=LOG(VRELQ)+DEFEL

cc evaluate cross-section for ELAB
      !! Coefficients are read from “1.2.7 Reaction 0.3T p + H2 -> p + H2 total cross section”
      !! Polynomial is evaluated for log ELAB (see Eq. 3.19 in Kotov)
      CEL = b2_cross(elab) !! total cross-section !! TODO




cc determine impact parameter b
      PIAI = 1/pi
      rb = 0.52917E-8_R8  !! = Bohr radius (in cm)
      BMAX=SQRT(CEL*PIAI)/rb !! from Kotov p132: “Here bmax is found from the condition: σ^t=Πb_{max}^2 “
      call random_number(randx)
      B= SQRT(randx)*BMAX

cc TODO iflag
      IFLAG = 2

cc find distance off shortest approach
      ER = RLMS*VRELQ*CVELI2 !! = E_r in Kotov 
      RS=B2_RSTERN(ER,B,IFLAG,P) !! This function solves Eq. B.4 for the distance of closest approach, r*

cc Kotov p132: “To find the deflection function _ the integral (3.5) is calculated using a Gauss-Mehler type quadrature rule with 16 nodes (subroutine GAUMEH)”
      call B2_GAUMEH (RS,ER,B,IFLAG,P,10,1,RES)
      CHI=pi-2.0_R8*B/RS*RES !! determine Chi from Eq. 3.5 in Kotov

cc determine the scattering angle Θ (PH) from Chi (Eq 3.5)
      CCHI=COS(CHI)
      PH=ACOS(CCHI)
      CPH=COS(PH)  !!= cos⁡(Θ)
      SPH=SQRT(1.0-CPH*CPH) !! = sin⁡(Θ)

cc a random number between  0 and 2 Pi is chosen for the polar angle
      call random_number(randx)
      EPS=2.0_R8*pi*randx

cc convert post-collisional velocities back into lab frame
      CEPS=COS(EPS)
      SEPS=SIN(EPS)
      !! VRS = Post-collisional relative velocities in center-off-mass system:
      !! VS = center-off-mass velocity (does not change after collision!)
      VRSX=VRELX*CPH+SPH*SEPS*VRYZ
      VRSY=VRELY*CPH+SPH*(VREL*VRELZ*CEPS-VRELX*VRELY*SEPS)/VRYZ
      VRSZ=VRELZ*CPH-SPH*(VREL*VRELY*CEPS+VRELX*VRELZ*SEPS)/VRYZ
      !! Post-collisional velocity in lab-frame (see Eq 3.24 in Kotov):
      VELX=VSX+RLMS/RMN*VRSX
      VELY=VSY+RLMS/RMN*VRSY
      VELZ=VSZ+RLMS/RMN*VRSZ
      VELQ=VELX*VELX+VELY*VELY+VELZ*VELZ
      VEL=SQRT(VELQ)

cc output
      vxt_out = VELX  
      vyt_out = VELY
      vzt_out = VELZ
      vxb_out = 0.0_R8
      vyb_out = 0.0_R8
      vzb_out = 0.0_R8

* ..return
      call subend ()
      return


      end subroutine execute_elastic_collision



C #####################################################################################################
C The functions below are copies of EIRENE functions. Code duplication could not be avoided in this case,
C as these functions are needed on the standalone B2.5 side. 
C #####################################################################################################
C
      function b2_rstern (ER,B,IFLAG,P)
C
C     IFLAG=1: H+ + H, PURELY REPULSIVE POTENTIAL
C     IFLAG=2: MORSE POTENTIAL, He+ + He, H+ + Noble Gases, H+ + H2
C
C               PARAMETERS: P(1):
C                           P(2):
C                           P(3):
C
C  RSTERN(ER,B) IS THE LARGEST ROOT OF THE EQUATION:
C
C     FI(R):=1.-V(R)/ER-(B/R)**2=0.
C
C   HERE: ER COLLISION ENERGY (EV)
C         B  IMPACT PARAMETER
C         V  INTERACTION POTENTIAL (EV)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      real(kind = R8), intent(in) :: ER, B, P(*)
      integer, intent(in) :: IFLAG

      real(kind = R8) :: FIW, RV, B2_RTSAF, BQ, TOL, B2_RSTERN,
     .            B2_FI, DFI,
     .        FITEST, rup, rlw, vw, rw, r0
      external xerrab

      IF (IFLAG.EQ.1) THEN
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAF
C  FOR IFLAG=1, RLW=B IS A LOWER BOUND AND RUP=INF IS AN UPPER BOUND.
C  FIRST TRY TO FIND BETTER BOUNDS:
        RUP=B
C
   10   RLW=RUP
        RUP=RLW*2.
        FITEST=B2_FI(RUP,ER,B,IFLAG,P,DFI)
        IF (FITEST.LT.0.0_R8) GOTO 10
C
        TOL=1.0e-7_R8*(RUP-RLW)
        B2_RSTERN=B2_RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
C
      ELSEIF (IFLAG.EQ.2) THEN
C  GENERALISED MORSE POTENTIAL
C  FIND UPPER AND LOWER BOUND FOR NUMERICAL ROOT FINDER RTSAFE
C  FOR IFLAG=2 AND IFLAG=3, THE FOLLOWING INTERVAL CAN BE USED:
C  R0 IS THE ROOT OF THE INTERACTION POTENTIAL V(R): V(R0)=0.
C  RM IS THE RADIUS OF THE MINIMUM OF V(R)           V(RM)= VM = -EPS
C  RW IS THE RADIUS OF THE POINT OF INFLECTION       V(RM)= VW
C        RM=P(4)
        R0=P(5)
        RW=P(6)
        VW=P(8)
C
        BQ=B*B
        RV=R0

        IF (RV.GT.B) THEN !!  CASE 1: B < R0,
          RLW=B
          RUP=RV
        ELSEIF (RV.LT.B) THEN !! CASE 2: R0 < B,
          RLW=RV
          RUP=B
C  TRY TO IMPROVE LOWER BOUND RLW
          FIW=1.0_R8-VW/ER-BQ/(RW**2)
          IF (FIW.LT.0.0_R8) RLW=RW
        ELSEIF (RV.eq.0.0_R8.and.B.eq.0.0_R8) THEN
          call xerrab('B can not be 0')
        ELSE
          call xerrab('RV exactly equal to B??')
        ENDIF
C
        TOL=1.0e-7_R8*(RUP-RLW)
        B2_RSTERN=B2_RTSAF(RLW,RUP,TOL,ER,B,IFLAG,P)
      ELSE
        call xerrab('ERROR IN FUNCTION B2_RSTERN ')
      ENDIF
C
      RETURN

      END FUNCTION b2_rstern

      FUNCTION b2_fi(R,ER,B,IFLAG,P,DFI)

C  EVALUATE EFFECTIVE POTENTIAL FUNCTION FI AT R
C  EVALUATE DFI(R)/DR AT R
C  RETURN FI=FI(R), DFI=DFI(R)/DR
C     --------------
C  IFLAG=1: H+ + H
C  IFLAG=2: H+ + NOBLE GASES, H+ + H2, HE+ + HE
C

      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      REAL(kind=R8), INTENT(IN) :: R, ER, B, P(*)
      REAL(kind=R8), INTENT(OUT) :: DFI
      INTEGER, INTENT(IN) :: IFLAG

      REAL(kind=R8) :: DSS, U, DU, SS, RI, RIQ, G1, REFF, RR,
     .            EMRR, G2, EX2, B2, V, DV, EX, RLOW, R2, R3,
     .            B2_FI, G, RMI, EPS

      external xerrab

      B2=B*B

      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.0e-7_R8
        R2=R*R
        R3=R2*R
C  FIND V=V(R) and DV=DV(R)/DR
        IF (R.GT.160_R8) THEN
          V=0.0_R8
          DV=0.0_R8
        ELSEIF (R.LT.RLOW) THEN
          R2=RLOW*RLOW
          R3=R2*RLOW
          EX=EXP(-RLOW)
          EX2=EX*EX
          SS=(1.0_R8+RLOW+R2/3.0_R8)*EX-1.0_R8
          RI=1.0_R8/RLOW
          V=27.211_R8*((RI-(1.0_R8+RI)*EX2-(1.0_R8+RLOW)*EX)/SS+RI)
          DV=0.0_R8
        ELSE
          EX=EXP(-R)
          EX2=EX*EX
          SS=(1.0_R8+R+R2/3.0_R8)*EX-1.0_R8
          RI=1.0_R8/R
          V=27.211_R8*((RI-(1.0_R8+RI)*EX2-(1.0_R8+R)*EX)/SS+RI)
          RIQ=RI*RI
          DSS=-R/3.0_R8*(1.0_R8+R)*EX
          U=RI-(1.0_R8+RI)*EX2-(1.0_R8+R)*EX
          DU=-RIQ+(RIQ+2.0_R8*RI+2.0_R8)*EX2+R*EX
          DV=27.211_R8*((DU*SS-U*DSS)/(SS*SS)-RIQ)
        ENDIF
C  FIND FI=FI(R) AND DFI=DFI(R)/DR
        B2_FI=1.0_R8-V/ER-B2/R2
        DFI=-DV/ER+2.0_R8*B2/R3
C
      ELSEIF (IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1.0_R8/P(4)
        R2=R*R
        R3=R2*R
C
        RR=R*RMI
        EMRR=1.0_R8-RR
C
C       G2=1.00+(1.0-G2)*MAX(0.D0,-EMRR)/EMRR
C       REFF=-G1*G2*EMRR
        IF (RR.LT.1.0_R8) THEN
          G=G1
        ELSE
          G=G1*G2
        ENDIF
        REFF=-G*EMRR
C
        IF (REFF.GT.160.D0) THEN
          V=0.0_R8
          DV=0.0_R8
        ELSE
          EX=EXP(-REFF)
          EX2=EX*EX
          V=EPS*(EX2-(EX+EX))
          DV=-2.0_R8*EPS*RMI*G*(EX2-EX)
        ENDIF
        B2_FI=1.0_R8-V/ER-B2/R2
        DFI=-DV/ER+2.0_R8*B2/R3
C
      ELSE
        call xerrab('ERROR IN FUNCTION B2_FI. IFLAG INVALID.')


      ENDIF
C
      RETURN
      END FUNCTION b2_fi

      function b2_rtsaf(X1,X2,XACC,ER,B,IFLAG,P)

C  TAKEN FROM: NUMERICAL RECIPES, W.H.PRESS ET AL.,
C              CAMBRIDGE UNIV. PRESS, 1989, P258
C  MODIFIED, TO FIND LARGEST ROOT, IN CASE MORE THAN ONE ROOTS
C  AND TO SPEED UP

      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      REAL(kind=R8), INTENT(IN) :: X1, X2, XACC, ER, B, P(*)
      INTEGER, INTENT(IN) :: IFLAG
      INTEGER, PARAMETER :: MAXIT=100
      REAL(kind=R8) :: DF, XDIST, F, B2_FI, TEMP, DX, DXOLD, XH,
     .            B2_RTSAF, XL, POI
      INTEGER :: J

      XL=X1
      XH=X2
C
      F=B2_FI(XH,ER,B,IFLAG,P,DF)
      XDIST=(XH-XL)
      DX=MIN(XDIST,F/DF)
      DXOLD=DX
      B2_RTSAF=XH-DX
      POI = B2_RTSAF
C
      F=B2_FI(POI,ER,B,IFLAG,P,DF)

      IF(F.LT.0.0_R8) THEN
        XL=B2_RTSAF
      ELSE
        XH=B2_RTSAF
      ENDIF
C
      DO J=1,MAXIT
        IF(((B2_RTSAF-XH)*DF-F)*((B2_RTSAF-XL)*DF-F).GT.0.
     *      .OR. ABS(F+F).GT.ABS(DXOLD*DF) ) THEN
          DXOLD=DX
          DX=0.5_R8*(XH-XL)
          B2_RTSAF=XL+DX
          POI=B2_RTSAF
        ELSE
          DXOLD=DX
          DX=F/DF
          TEMP=B2_RTSAF
          B2_RTSAF=B2_RTSAF-DX
          POI=B2_RTSAF
        ENDIF
C
        IF(ABS(DX).LT.XACC) RETURN
C
        F=B2_FI(POI,ER,B,IFLAG,P,DF)
        IF(F.LT.0.0_R8) THEN
          XL=B2_RTSAF
        ELSE
          XH=B2_RTSAF
        ENDIF
      ENDDO
      WRITE(*,*) 'RTSAF EXCEEDING MAXIMUM ITERATIONS'
      RETURN
      END FUNCTION B2_RTSAF


      SUBROUTINE B2_GAUMEH(RS,ER,B,IFLAG,P,N,IGAUS,RES)
C
C  IGAUS=1
C  GAUSS MEHLER QUADRATURE, N=5,10,20, W(X)=1./SQRT(1+X)/SQRT(1-X)
C                           INTEGRATION FROM A=-1 TO B=1
C                           TRANSFORMED TO A=0,B=1, AND
C                           W(X)=1./SQRT(X)/SQRT(1-X)
C
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      REAL(kind=R8), INTENT(IN) :: RS, ER, B, P(*)
      REAL(kind=R8), INTENT(OUT) :: RES
      INTEGER, INTENT(IN) :: IFLAG, N, IGAUS

      REAL(kind=R8) :: AR(128), AFI(128)
      INTEGER :: NFI
      REAL(kind=R8) :: X5A(5),X5B(5),X10A(10),X10B(10),X20A(20),X20B(20)
      REAL(kind=R8) :: X5(5),X10(10),X20(20)
C
      REAL(kind=R8) :: XG10(10),WG10(10),XG10A(10),XG10B(10)
      REAL(kind=R8) :: SUM, F, X, W5, W10, W20
      INTEGER :: I, IFIRST, IFI
C
      SAVE
      DATA X5 /0.95105654E+00_R8, 0.58778542E+00_R8,0.31391647E-06_R8,
     .       - 0.58778363E+00_R8,-0.95105600E+00_R8/,
     .     W5 /0.62831837E+00_R8/
      DATA X10/9.876884E-01_R8, 8.910065E-01_R8, 7.071068E-01_R8,
     .         4.539905E-01_R8, 1.564344E-01_R8,-1.564345E-01_R8,
     .        -4.539905E-01_R8,-7.071068E-01_R8,-8.910065E-01_R8,
     .        -9.876884E-01_R8/,
     .     W10/3.141593E-01_R8/
      DATA X20/9.969173E-01_R8, 9.723699E-01_R8, 9.238795E-01_R8,
     .         8.526402E-01_R8, 7.604060E-01_R8, 6.494480E-01_R8,
     .         5.224985E-01_R8, 3.826834E-01_R8, 2.334453E-01_R8,
     .         7.845905E-02_R8,-7.845914E-02_R8,-2.334454E-01_R8,
     .        -3.826835E-01_R8,-5.224986E-01_R8,-6.494481E-01_R8,
     .        -7.604060E-01_R8,-8.526402E-01_R8,-9.238796E-01_R8,
     .        -9.723700E-01_R8,-9.969174E-01_R8/,
     .     W20/1.570796E-01_R8/
C
C  IGAUS=2
C  GAUSS MEHLER QUADRATURE, N=10, W(X)=1./SQRT(1-X)
C                           INTEGRATION FROM A=0 TO B=1
      DATA XG10/
     .    .013695585480651072296_R8,
     .    .070758123420104485088_R8,
     .    .16782834791297652828_R8,
     .    .29588270759990962383_R8,
     .    .44298868539955669668_R8,
     .    .59543571523425259033_R8,
     .    .73901490631777356542_R8,
     .    .86034375925702257262_R8,
     .    .94811360601967719799_R8,
     .    .99414369156320387287_R8/
      DATA WG10/
     .    .035228014278304236584_R8,
     .    .081202859600773882686_R8,
     .    .12534409666821812627_R8,
     .    .16655348315340949082_R8,
     .    .20386023963448093937_R8,
     .    .23638906392303662412_R8,
     .    .26337727689835336360_R8,
     .    .28419221863676482994_R8,
     .    .29834597294520589473_R8,
     .    .30550677426145261144_R8/
      DATA IFIRST /0/
      IF (IFIRST.EQ.0) THEN
C SET ROOTS AND WEIGHTS FOR GAUSS QUADRATURE RULES
        IFIRST=1
        DO I=1,5
          X=0.5_R8*(1.0_R8+X5(I))
          X5A(I)=SQRT(X*(1.0_R8-X))
          X5B(I)=1.0_R8/X
        END DO
        DO I=1,10
          X=0.5_R8*(1.0_R8+X10(I))
          X10A(I)=SQRT(X*(1.0_R8-X))
          X10B(I)=1.0_R8/X
        END DO
        DO I=1,10
          X=XG10(I)
          XG10A(I)=SQRT(1.0_R8-X)
          XG10B(I)=1.0_R8/X
        END DO
        DO I=1,20
          X=0.5_R8*(1.0_R8+X20(I))
          X20A(I)=SQRT(X*(1.0_R8-X))
          X20B(I)=1.0_R8/X
        END DO
      ENDIF

      NFI=N
      IF (N.EQ.5) THEN
        DO IFI=1,5
          AR(IFI)=RS*X5B(IFI)
        END DO

        CALL B2_FIVEC(ER,B,IFLAG,P,NFI,AR,AFI)
        SUM=0.0_R8
        DO IFI=1,5
          F=AFI(IFI)
          IF (F.GT.0.0_R8)
     .    SUM=SUM+X5A(IFI)/SQRT(F)
        END DO
        RES=W5*SUM

      ELSEIF (N.EQ.10) THEN
        IF (IGAUS.EQ.1) THEN
          DO IFI=1,10
            AR(IFI)=RS*X10B(IFI)
          END DO

          CALL B2_FIVEC(ER,B,IFLAG,P,NFI,AR,AFI)
          SUM=0.0_R8
          DO IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.0_R8)
     .      SUM=SUM+X10A(IFI)/SQRT(F)
          END DO
          RES=W10*SUM
        ELSEIF (IGAUS.EQ.2) THEN
          DO IFI=1,10
            AR(IFI)=RS*XG10B(IFI)
          END DO

          CALL B2_FIVEC(ER,B,IFLAG,P,NFI,AR,AFI)
          SUM=0.0_R8
          DO IFI=1,10
            F=AFI(IFI)
            IF (F.GT.0.D0)
     .      SUM=SUM+XG10A(IFI)/SQRT(F)*WG10(IFI)
          END DO
          RES=SUM
        ENDIF

      ELSEIF (N.EQ.20) THEN
        DO IFI=1,20
          AR(IFI)=RS*X20B(IFI)
        END DO

        CALL B2_FIVEC(ER,B,IFLAG,P,NFI,AR,AFI)
        SUM=0.0_R8
        DO IFI=1,20
          F=AFI(IFI)
          IF (F.GT.0.0_R8)
     .    SUM=SUM+X20A(IFI)/SQRT(F)
        END DO
        RES=W20*SUM

      ELSE
        call xerrab('ERROR IN GAUMEH, WRONG PARAMETER N')
      ENDIF
      RETURN
      END SUBROUTINE B2_GAUMEH


      SUBROUTINE B2_FIVEC(ER,B,IFLAG,P,NFI,AR,AFI)

C  VECTORIZED VERSION OF FUNCTION FI FOR GAUSS MEHLER QUADRATURE
C  EVALUATE EFFECTIVE POTENTIAL FUNCTION AT AR(I),I=1,NFI
C  NOTE: NFI.LE.128 IS NOT CHECKED, BUT USED
C  RETURN FI(AR(I)) IN THE ARRAY AFI(I),I=1,NFI
C     --------------
C  IFLAG=1: H+ + H
C  IFLAG=2: H+ + NOBLE GASES, H+ + H2, HE+ + HE
C

      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      REAL(kind=R8), INTENT(IN) :: ER, B, P(*), AR(128)
      INTEGER, INTENT(IN) :: IFLAG, NFI
      REAL(kind=R8), INTENT(OUT) :: AFI(128)
      REAL(kind=R8) :: RI, G1, SS, EX, EX2, RR, REFF, EMRR, G2, V,
     .          R, R2, B2, RLOW, G, RMI, EPS
      INTEGER :: IFI


      B2=B*B

      IF(IFLAG.EQ.1) THEN
C  INTERACTION POTENTIAL V(R): H+ + H
C  R IN A0, V IN EV
        RLOW=1.0e-7_R8
        DO IFI=1,NFI
          R=AR(IFI)
          R2=R*R
C  FIND V=V(R)
          IF (R.GT.160.D0) THEN
            V=0.0_R8
          ELSEIF (R.LT.RLOW) THEN
            R2=RLOW*RLOW
            EX=EXP(-RLOW)
            EX2=EX*EX
            SS=(1.0_R8+RLOW+R2/3.0_R8)*EX-1.0_R8
            RI=1.0_R8/RLOW
            V=27.211_R8*((RI-(1.0_R8+RI)*EX2-(1.0_R8+RLOW)*EX)/SS+RI)
          ELSE
            EX=EXP(-R)
            EX2=EX*EX
            SS=(1.0_R8+R+R2/3.0_R8)*EX-1.0_R8
            RI=1.0_R8/R
            V=27.211_R8*((RI-(1.0_R8+RI)*EX2-(1.0_R8+R)*EX)/SS+RI)
          ENDIF
          AFI(IFI)=1.0_R8-V/ER-B2/R2
        END DO
C
      ELSEIF(IFLAG.EQ.2) THEN
C  INTERACTION POTENTIAL V(R): H+ + NOBLE GASES, (MORSE LIKE POTENTIAL)
C     R IN A0, V IN EV
        EPS=P(1)
        G1=P(2)
        G2=P(3)
        RMI=1.0/P(4)
        DO IFI=1,NFI
          R=AR(IFI)
          R2=R*R
          RR=R*RMI
          EMRR=1.0_R8-RR
          IF (RR.LT.1) THEN
            G=G1
          ELSE
            G=G1*G2
          ENDIF
          REFF=-G*EMRR
C
          IF (REFF.GT.160.0_R8) THEN
            V=0.0_R8
          ELSE
            EX=EXP(-REFF)
            EX2=EX*EX
            V=EPS*(EX2-(EX+EX))
          ENDIF
          AFI(IFI)=1.0_R8-V/ER-B2/R2
        END DO
C
      ELSE
        call xerrab('ERROR IN FUNCTION EIRENE_FIVEC. IFLAG INVALID.')
      ENDIF
C
      RETURN
      end subroutine b2_fivec


      SUBROUTINE B2_FGAUSS(FG1,FG2,FG3)
C
C  SAMPLE IRNDVC TRIPLES OF RANDOM NUMBER FROM A GAUSSIAN
C  AND RETURN THEM IN THE ARRAYS FG1,FG2,FG3
C  IT FIRST GENERATES IRNDVC PAIRS, AND THEN IRNDVC/2 PAIRS
C  THIS ROUTINE IS CALLED IF INIV2=0. IT RETURNS INIV2=IRNDVC

      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
cc output variables
      REAL(kind=R8) :: FG1, FG2, FG3
cc local variables
      REAL(kind=R8) :: ZCOS, ZSIN, ZT, AR, ZZ, xrand

C    THIS FUNCTION SAMPLES FROM A GAUSSIAN OF THE
C    FORM EXP(-(X-ZA)**2/(2.*ZS**2))/(ZS*SQRT(2*PI))
C    ZA=0.
C    ZS=1.
C
C  IT IS THE BOX-MUELLER METHOD
C
C   1. COMPUTE THE SINE AND COSINE OF 2*PI*RAN(1)
C
      call random_number(xrand)
      ZZ=2.0_R8*pi*xrand
      ZSIN=SIN(ZZ)
      ZCOS=COS(ZZ)
      call random_number(xrand)
      AR=LOG(xrand)
      ZT=SQRT(-2.0_R8*AR)
      FG1=ZT*ZSIN
      FG2=ZT*ZCOS

      call random_number(xrand)
      ZZ=2.0_R8*pi*xrand
      ZSIN=SIN(ZZ)
      ZCOS=COS(ZZ)
      call random_number(xrand)
      AR=LOG(xrand)
      ZT=SQRT(-2.0_R8*AR)
      call random_number(xrand)
      if (xrand.gt.0.5_R8) then
        FG3=ZT*ZSIN
      else
        FG3=ZT*ZCOS
      endif

      RETURN
      END SUBROUTINE B2_FGAUSS



      subroutine b2_dbl_poly (cf, al1, al2, cou, dum,
     .           rc1min, rc1max, rc2min, rc2max)

c  used for evaluating double polynomial fit.
c  return        COU=fit2(AL1,AL2)
C  INTERNAL INTERPRETATION: AL1=LOG(PARM1),
C                           AL2=LOG(PARM2), WITH PHYSICAL PARAMETERS PARM1, PARM2
c  and return also:
c  reduced fit coefficients dum(1,9) for AL2 dependency at fixed first parameter AL1.
c  Then COU=fit1(AL2)=sum_1^9 [dum(j) log(parm2)^(j-1)]

c  input:
c  nd1, nd2: orders of double fit polynomial: nd1-1,nd2-1 (here still: nd1=nd2=9, merging tbd.)
c  cf    : fit coefficients for 2D polyn. fit. We require data: cf(1:9,1:9)
c          f(parm1,parm2)=sum_1^9 [sum_1^9 (cf(i,j) log(parm1)^(i-1)] log(parm2)^(j-1)
c  al1   : 1st argument of fit, log(parm1)  e.g. Te, Ti, ...
c  al2   : 2nd argument of fit, log(parm2)  e.g. ne, E0, ...

c  rc1min : left boundary of valid range for first parameter al1
c  rc1max : right boundary of valid range for first parameter al1
c  fp1    : parameters for extrapolation from valid range

c  rc2min : bottom boundary of valid range for second parameter al2
c  rc2max : top    boundary of valid range for second parameter al2
c  fp2    : parameters for extrapolation from valid range

c  output:
c  cou   : log of value of fit(al1,al2)
c  dum   : dum(1,9), summation over 1st parameter is done,
c          dum provides 1D fit coefficients for 2nd parameter dependency (e.g. ne, E0)
c          for reduced 1D fit evaluation "on the fly", at fixed parameter AL1.
c          Also: dum(1) is the fit value at al2=0, e.g. at density=1e8 for amjuel.


      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      real(kind=R8), intent(in) :: cf(9,9)
      real(kind=R8), intent(in) :: al1, al2, rc1min, rc1max, 
     .                                       rc2min, rc2max
      real(kind=R8), intent(out) :: cou, dum(9)

      real(kind=R8) :: p1, p2
      integer :: kk, jj

      dum = 0.0_R8
cc if outside of valid range, evaluate at the boundary (no extrapolation)
      p1 = max(min(rc1max,al1),rc1min)
      p2 = max(min(rc2max,al2),rc2min)

cc  split the fit evaluation into two steps.
cc  first:  p1=al1 dependence, to provide collapsed fit dum(..)
cc  second: p2=al2 dependence using dum(..)

cc  each of the 9 coefficients dum(jj) for the al2 dependence
cc  is obtained by summing over all 9 terms
      do jj = 9, 1, -1
        dum(jj) = cf(9,jj)
        do kk = 8, 1, -1
          dum(jj) = dum(jj) * p1 + cf(kk,jj)
        end do
      end do

cdr  this second evaluation may not be needed, if only a collapsed fit is wanted,
cdr  or if al2=0.0 (as in H.4, H.10, AMJUEL fits at 1e8 density), for automatic corona limit.

cdr  H.4, H.10, H.12 fits from AMJUEL: corona at p2 <= log(ne/10**8) = rc2min = 0.0
cdr  if p2.le.0.0, just return cou=dum(1) = fit2(AL1,AL2)

      cou = dum(9)
      do jj = 8, 1, -1
        cou = cou * p2 + dum(jj)
      end do

      return
      end subroutine b2_dbl_poly




      FUNCTION B2_CROSS(AL)

      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none


      REAL(kind=R8), INTENT(IN) :: AL
      REAL(kind=R8), save :: POLY(9), FP(6), FPL(3), FPR(3),
     .                       rcmin, rcmax
      REAL(kind=R8) :: B2_CROSS,
     .            RES, B2_EXTRAP, B2_SNGL_POLY
      INTEGER :: JFEX1MN, JFEX1MX

cc 1.2.7 Reaction 0.3T p + H2 ! p + H2 total cross section. Hardcoded in testing phase TODO generalize, read from AMJUEL
      data poly 
     . /-3.452141819446e+01_R8, 1.092015526305e+01_R8,
     .  -2.732690257819e+01_R8, 3.466297654768e+01_R8,
     .  -2.524607958646e+01_R8, 1.092376446349e+01_R8, 
     .  -2.770065796605e+00_R8, 3.796353200921e-01_R8,
     .  -2.168988142310e-02_R8/
      data FPL /-3.275286840950e+01_R8,
     . -2.351764912137e-01_R8, -1.045602118569e-02_R8/
      data FPR /-3.537275807146e+01_R8,
     .  2.144573517210e-01_R8, -4.643079956637e-02_R8/

      FP(1:3) = FPL
      FP(4:6) = FPR

      JFEX1MN = 5 !! type of extrap. scheme
      JFEX1MX = 5

cc ELABMIN and ELABMAX from AMJUEL 1.2.7 Reaction 0.3T p + H2 ! p + H2 total cross section
      rcmin = log( 1.55980_R8)
      rcmax = log( 6.18164E+01_R8)
      
      B2_CROSS = B2_SNGL_POLY(poly,AL,
     .                            rcmin,rcmax,FP,
     .                            JFEX1MN,JFEX1MX)
      RETURN
      END FUNCTION B2_CROSS


      function b2_sngl_poly (cf, al, rcmin, rcmax, fpp, ifexmn, ifexmx)
     .         result(cou)
c  input:
c  cf      : fit coefficients for fit POLY=f(parm=)=sum_1^9 (cf(i) log(parm)^(i-1))
c  al      : argument of fit, al=log(parm)
c  rcmin   : left boundary of valid range of PARM
c  rcmax   : right boundary of valid range of PARM
c  fpp(1:6): parameters for extrapolation from valid range
c            1:3  left (low PARM values) extrapolation.
c            4:6  right (high PARM values) extrapolation
c output: return exp(POLY)


      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      real(kind=R8), intent(in) :: cf(9), fpp(6)
      real(kind=R8), intent(in) :: al, rcmin, rcmax
      integer, intent(in) :: ifexmn, ifexmx
      real(kind=R8) :: p1, parm, cou, fp(6), 
     .            almin, almax, coumin, coumax,
     .            b2_extrap
      integer :: ii, ifex

      p1=al

      if (p1 < rcmin.and.ifexmn.gt.0) then

        coumin = cf(9)
        do ii = 8, 1, -1
          coumin = coumin * RCMIN + cf(ii)
        end do
        ifex = ifexmn
cc  determine parameter and fit value at left boundary. 
        ALMIN= RCMIN
        COUMIN=EXP(COUMIN)
        FP(1:3) = FPP(1:3)
        COU=B2_EXTRAP(P1,IFEX,FP(1),FP(2),FP(3))
          cou = log(cou)

      elseif (p1 > rcmax.and.ifexmx.gt.0) then

        COUMAX = cf(9)
        do ii = 8, 1, -1
          coumax = coumax * RCMAX + cf(ii)
        end do
        ifex = ifexmx
C  determine parameter and fit value at right boundary. may be needed by fct. extrap
        ALMAX= RCMAX
        COUMAX=EXP(COUMAX)
        FP(4:6) = FPP(4:6)

        COU=B2_EXTRAP(P1,IFEX,FP(4),FP(5),FP(6))
        cou = log(cou)

      else

        p1 = min(max(p1,rcmin),rcmax)

        cou = cf(9)
        do ii = 8, 1, -1
          cou = cou * p1 + cf(ii)
        end do

      endif 

      cou = exp(max(-100.0_R8, cou))

      return
      end function b2_sngl_poly


      FUNCTION B2_EXTRAP(PARM,IFLAG,FP1,FP2,FP3)
C
C  FUNCTION FOR EXTRAPOLATING SINGLE PARAMETER FITS "FF(PP)" BEYOND THEIR
C  RANGE OF VALIDITY, WHEN LN(PP)=PARM IS OUT OF RANGE. RETURN EXTRAP=FF(PARM)
C  For now only iflag=5 implemented here
C  IFLAG=5 : RETURN EXTRAP=EXP(FP1+FP2*PARM+FP3*PARM**2), 2ND ORDER ON LOG SCALE
C
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none

      REAL(kind=R8), INTENT(IN) :: PARM, FP1, FP2, FP3
      INTEGER, INTENT(IN)  :: IFLAG
      REAL(kind=R8)             :: B2_EXTRAP


      IF (IFLAG.eq.0) THEN 
      ELSEIF (IFLAG.EQ.5) THEN
C  QUADRATIC EXTRAPOLATION OF LN(FF) IN PARM=LN(PP)
CNR FLOOR AT EPS30 TO AVOID LOG(0) IN POLY ROUTINES
        B2_EXTRAP=EXP(MAX(-100.0_R8,FP1+FP2*PARM+FP3*PARM**2))
      ELSE
        call xerrab( 'ERROR IN EXTRAP. IFLAG OUT OF RANGE.')
      ENDIF
      RETURN

      END FUNCTION B2_EXTRAP

