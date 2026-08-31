*-----------------------------------------------------------------------
*.specification

      subroutine b2smcu_3D (nCv, nFc, nVx, mb, switch, geo, mpg,   
     *           nb, u_pol, u_rad, u_tor, sm_pol, sm_rad, sm_tor)                   
      use b2mod_types
      use b2mod_constants
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx                               
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   mb, nb(nCv), u_pol(nCv), u_rad(nCv), u_tor(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   sm_pol(nCv), sm_rad(nCv), sm_tor(nCv)                       
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SMCU_3D computes terms from the divergence of the stress tensor
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv
      real (kind=R8) :: wrks(nCv), zeroes_Vx(nVx)
      real (kind=R8) ::
     *   wrk0(nCv), wrk1(nCv), wrk2(nCv), wrk3(nCv), wrk4(nCv),
     *   wrk5(nCv), t0, t1, t2, t3, t4, t5, t6
*   ..procedures
      external xertst
      external b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2smcu_3D')
*   ..set internal parameters on first call
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test mb
      call xertst (0.le.mb, 'faulty argument mb')

      sm_pol = 0.0_R8 
      sm_rad = 0.0_R8 
      sm_tor = 0.0_R8
      zeroes_Vx = 0.0_R8

* ..compute centrifugal force
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,geo%cvHz,zeroes_Vx,wrk0)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,geo%cvHx,zeroes_Vx,wrk1)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,geo%cvHy,zeroes_Vx,wrk2)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,geo%cvHz,zeroes_Vx,wrk3)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,geo%cvHx,zeroes_Vx,wrk4)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,geo%cvHy,zeroes_Vx,wrk5)

*   .. compute curvature terms
        do iCv = 1, mpg%nCi
          t0 = geo%cvVol(iCv)*mb*u_rad(iCv)*nb(iCv)*
     .         (wrk4(iCv)/geo%cvHx(iCv)*u_pol(iCv) -
     .          wrk2(iCv)/geo%cvHy(iCv)*u_rad(iCv))
          t1 = geo%cvVol(iCv)*mb*u_tor(iCv)*nb(iCv)*
     .         (-wrk0(iCv)/geo%cvHz(iCv)*u_tor(iCv))    
          t2 = geo%cvVol(iCv)*mb*u_pol(iCv)*nb(iCv)*
     .         (wrk2(iCv)/geo%cvHy(iCv)*u_rad(iCv) -
     .          wrk4(iCv)/geo%cvHx(iCv)*u_pol(iCv))  
          t3 = geo%cvVol(iCv)*mb*u_tor(iCv)*nb(iCv)*
     .         (-wrk3(iCv)/geo%cvHz(iCv)*u_tor(iCv))  
          sm_pol(iCv) = t0 + t1 
          sm_rad(iCv) = t2 + t3 
        enddo





c      if (switch%b2smcu_3D_iout.ne.0) then                         
c#ifdef _OPENMP
c!$OMP MASTER
c        write(*,*) 'B2SMCU_3D OpenMP warning: no file output in ',
c     &             'parallel mode'
c!$OMP END MASTER
c#else
c        call my_out_us(70,nCv,0,sm_rad,'b2smcu_rad')
c#endif
c      endif                                                          

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sicf

      end subroutine b2smcu_3D
