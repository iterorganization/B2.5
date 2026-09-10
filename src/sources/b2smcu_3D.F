*-----------------------------------------------------------------------
*.specification

      subroutine b2smcu_3D (nCv, nFc, nVx, mb, switch, geo, mpg,   
     *          nb, u_pol, u_rad, u_tor, vs_mol, sm_pol, sm_rad, sm_tor)                   
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
     *   mb, nb(nCv), u_pol(nCv), u_rad(nCv), u_tor(nCv), vs_mol(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   sm_pol(nCv,0:1), sm_rad(nCv,0:1), sm_tor(nCv,0:1)                       
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SMCU_3D computes terms from the divergence of the stress tensor
*     This routines use the equations from the Master thesis of 
*     N. Horsten & P. Bauweraerts, Appdendix B2, for the computation of
*     curvature terms from the convective term and the divergence of the 
*     stress tensor
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, iCv1, iCv2, iCv_guard, iCv_int
      real (kind=R8) :: wrks(nCv), zeroes_Vx(nVx)
      real (kind=R8) ::
     *   dhzdx(nCv), dhxdx(nCv), dhydx(nCv), dhzdy(nCv), dhxdy(nCv),
     *   dhydy(nCv), t0, t1, t2, t3, t4, t5, t6,
     *   dupdp(nCv), durdp(nCv), dutdp(nCv), 
     *   dupdr(nCv), durdr(nCv), dutdr(nCv),
     *   sigma_p_p(nCv), sigma_p_r(nCv), sigma_p_t(nCV),
     *   sigma_r_r(nCv), sigma_r_t(nCv), sigma_t_t(nCv),
     *   vol(nCv), hz(nCv), hx(nCv), hy(nCv),
     *   termGx(nCv), termGy(nCv), dtermGx(nCv), dtermGy(nCv),
     *   term_p_p(nCv), term_p_r(nCv), term_r_p(nCv), term_r_r(nCv),
     *   term_t_p(nCV), term_t_r(nCV),
     *   dterm_p_p(nCv), dterm_p_r(nCv), dterm_r_p(nCv), dterm_r_r(nCv),
     *   dterm_t_p(nCV), dterm_t_r(nCV),
     *   div_Pi_p(nCv), div_Pi_r(nCv), div_Pi_t(nCv),
     *   curv_conv_pol(nCv), curv_conv_rad(nCV), curv_conv_tor(nCv)
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

*   ..local copies for variables that are used a lot
      vol = geo%cvVol 
      hz = geo%cvHz 
      hy = geo%cvHy 
      hx = geo%cvHx
      !! apparently hy is not defined in guard cells -> this gives issues later on when computing gradients->manually recompute here 
cc      do iCv = mpg%nCi+1, mpg%nCv 
cc        hy(iCv) = vol(iCv)/hx(iCV)/hz(iCv)
cc      enddo 
      !! set the guard cell values to the internal ones, to avoid strong gradients at the boundaries 
      do iFc = 1, nFc 
        iCv1 = mpg%fcCv(iFc,1)
        iCv2 = mpg%fcCv(iFc,2)
        iCv_guard = max(iCv1,iCv2)
        iCv_int = min(iCv1,iCv2)
        if (iCv_guard.gt.mpg%nCi) then 
          hx(iCv_guard) = hx(iCv_int)
          hy(iCv_guard) = hy(iCv_int)
          hz(iCv_guard) = hz(iCv_int)
          vol(iCv_guard) = vol(iCv_int)
        endif
      enddo
      !! FM_TODO CHECK !!! -> rename wrk arrays to avoid confusion
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,hz,zeroes_Vx,dhzdx)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,hx,zeroes_Vx,dhxdx)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,hy,zeroes_Vx,dhydx)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,hz,zeroes_Vx,dhzdy)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,hx,zeroes_Vx,dhxdy)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,hy,zeroes_Vx,dhydy)

      call gradc_p(nCv,nFc,nVx,0,geo,mpg,u_pol,zeroes_Vx,dupdp)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,u_rad,zeroes_Vx,durdp)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,u_tor,zeroes_Vx,dutdp)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,u_pol,zeroes_Vx,dupdr)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,u_rad,zeroes_Vx,durdr)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,u_tor,zeroes_Vx,dutdr)


*   .. compute curvature pieces of convective terms
      do iCv = 1, mpg%nCi
        t0 = vol(iCv)*mb*u_rad(iCv)*nb(iCv)*
     .         (dhxdy(iCv)/hx(iCv)*u_pol(iCv) -
     .          dhydx(iCv)/hy(iCv)*u_rad(iCv))
        t1 = vol(iCv)*mb*u_tor(iCv)*nb(iCv)*
     .         (-dhzdx(iCv)/hz(iCv)*u_tor(iCv))    
        t2 = vol(iCv)*mb*u_pol(iCv)*nb(iCv)*
     .         (dhydx(iCv)/hy(iCv)*u_rad(iCv) -
     .          dhxdy(iCv)/hx(iCv)*u_pol(iCv))  
        t3 = vol(iCv)*mb*u_tor(iCv)*nb(iCv)*
     .         (-dhzdy(iCv)/hz(iCv)*u_tor(iCv))  
        t4 = vol(iCv)*mb*u_pol(iCv)*nb(iCv)*
     .         (dhzdx(iCv)/hz(iCv)*u_tor(iCv))  
        t5 = vol(iCv)*mb*u_rad(iCv)*nb(iCv)*
     .         (dhzdy(iCv)/hz(iCv)*u_tor(iCv)) 
        sm_pol(iCv,0) = -t0 - t1 
        sm_rad(iCv,0) = -t2 - t3 
        sm_tor(iCv,0) = -t4 - t5
        curv_conv_pol(iCv) = -t0 - t1 
        curv_conv_rad(iCv) = -t2 - t3 
        curv_conv_tor(iCv) = -t4 - t5
      enddo

*   .. compute gradients of sqrt(g)/hx or hy
      termGx = vol/hx
      termGy = vol/hy 
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,termGx,zeroes_Vx,dtermGx)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,termGy,zeroes_Vx,dtermGy)
      termGx = dtermGx/vol*hx
      termGy = dtermGy/vol*hy




*   .. compute components of stress tensor
!! FM_TODO: multiply with density / cell volume / mass ??
      sigma_p_p = vs_mol*(
     *   4.0_R8/3.0_R8*dupdp - 2.0_R8/3.0*durdr 
     *  -2.0_R8/3.0_R8*u_pol*termGx  
     * + 2.0_R8*u_rad*(dhxdy/hx - 1.0_R8/3.0_R8*termGy))

      sigma_p_r = vs_mol*(
     * dupdr + durdp - u_pol/hx*dhxdy - u_rad/hy*dhydx)

      sigma_p_t = vs_mol*(
     * dutdp - u_tor/hx*dhzdx)

      sigma_r_r = vs_mol*(
     *   4.0_R8/3.0_R8*durdr - 2.0_R8/3.0_R8*dupdp 
     * + 2.0_R8*u_pol*(dhydx/hy-1.0_R8/3.0_R8*termGx)
     *  -2.0_R8/3.0_R8*u_rad*termGy )

      sigma_r_t = vs_mol*(
     *  dutdr - u_tor/hz*dhzdy)

      sigma_t_t = vs_mol*(
     *  -2.0_R8/3.0_R8*dupdp-2.0_R8/3.0_R8*durdr
     * + 2.0_R8*u_pol*(dhzdx/hz-1.0_R8/3.0_R8*termGx)
     * + 2.0_R8*u_rad*(dhzdy/hz-1.0_R8/3.0_R8*termGy))

*   .. compute terms of the form 1/sqrt(g) d/dx (sqrt(g)/hx sigma)
cc first compute the term inside the derivative, then the derivative, then add the prefactor
      term_p_p = vol/hx*sigma_p_p 
      term_p_r = vol/hy*sigma_p_r 
      term_r_p = vol/hx*sigma_p_r 
      term_r_r = vol/hy*sigma_r_r 
      term_t_p = vol/hx*sigma_p_t 
      term_t_r = vol/hy*sigma_r_t 
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,term_p_p,zeroes_Vx,dterm_p_p)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,term_p_r,zeroes_Vx,dterm_p_r)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,term_r_p,zeroes_Vx,dterm_r_p)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,term_r_r,zeroes_Vx,dterm_r_r)
      call gradc_p(nCv,nFc,nVx,0,geo,mpg,term_t_p,zeroes_Vx,dterm_t_p)
      call gradc_r(nCv,nFc,nVx,0,geo,mpg,term_t_r,zeroes_Vx,dterm_t_r)
      term_p_p = dterm_p_p*hx/vol
      term_p_r = dterm_p_r*hy/vol
      term_r_p = dterm_r_p*hx/vol
      term_r_r = dterm_r_r*hy/vol
      term_t_p = dterm_t_p*hx/vol
      term_t_r = dterm_t_r*hy/vol


*   .. compute terms related to divergence of the stress tensor (+ multiply by cell volume to get source per cell)
      div_Pi_p = vol*(term_p_p + term_p_r + sigma_p_r/hx*dhxdy 
     *           - sigma_r_r/hy*dhydx - sigma_t_t/hz*dhzdx) 
      div_Pi_r = vol*(term_r_p + term_r_r + sigma_p_r/hy*dhydx 
     *           - sigma_p_p/hx*dhxdy - sigma_t_t/hz*dhzdy) 
      div_Pi_t = vol*(term_t_p + term_t_r + sigma_p_t/hz*dhzdx 
     *           + sigma_r_t/hz*dhzdy)


      do iCv = 1, mpg%nCi
        sm_pol(iCv,0) = sm_pol(iCv,0) - div_Pi_p(iCv)
        sm_rad(iCv,0) = sm_rad(iCv,0) - div_Pi_r(iCv)
        sm_tor(iCv,0) = sm_tor(iCv,0) - div_Pi_t(iCv)
      enddo

cc linearize
      do iCv = 1, mpg%nCi
        sm_rad(iCv,0) = sm_rad(iCv,0)+
     &     max(0.0_R8,sm_rad(iCv,1))*u_rad(iCv)
        sm_rad(iCv,1) = min(sm_rad(iCv,1),0.0_R8)

        sm_pol(iCv,0) = sm_pol(iCv,0)+
     &     max(0.0_R8,sm_pol(iCv,1))*u_pol(iCv)
        sm_pol(iCv,1) = min(sm_pol(iCv,1),0.0_R8)

        sm_tor(iCv,0) = sm_tor(iCv,0)+
     &     max(0.0_R8,sm_tor(iCv,1))*u_tor(iCv)
        sm_tor(iCv,1) = min(sm_tor(iCv,1),0.0_R8)
      enddo


      if (switch%b2smcu_3D_iout.ne.0) then                         
#ifdef _OPENMP
!$OMP MASTER
        write(*,*) 'B2SMCU_3D OpenMP warning: no file output in ',
     &             'parallel mode'
!$OMP END MASTER
#else
      call my_out_us(70,nCv,0,div_Pi_p,'b2smcu_3D_div_Pi_p')
      call my_out_us(70,nCv,0,div_Pi_r,'b2smcu_3D_div_Pi_r')
      call my_out_us(70,nCv,0,div_Pi_t,'b2smcu_3D_div_Pi_t')
      call my_out_us(70,nCv,0,curv_conv_pol,'b2smcu_3D_curv_conv_pol')
      call my_out_us(70,nCv,0,curv_conv_rad,'b2smcu_3D_curv_conv_rad')
      call my_out_us(70,nCv,0,curv_conv_tor,'b2smcu_3D_curv_conv_tor')
      call my_out_us(70,nCv,0,sigma_p_p,'b2smcu_3D_sigma_p_p')
      call my_out_us(70,nCv,0,sigma_p_r,'b2smcu_3D_sigma_p_r')
      call my_out_us(70,nCv,0,sigma_p_t,'b2smcu_3D_sigma_p_t')
      call my_out_us(70,nCv,0,sigma_r_r,'b2smcu_3D_sigma_r_r')
      call my_out_us(70,nCv,0,sigma_r_t,'b2smcu_3D_sigma_r_t')
      call my_out_us(70,nCv,0,sigma_t_t,'b2smcu_3D_sigma_t_t')
      call my_out_us(70,nCv,0,term_p_p,'b2smcu_3D_term_p_p')
      call my_out_us(70,nCv,0,term_p_r,'b2smcu_3D_term_p_r')
      call my_out_us(70,nCv,0,term_r_p,'b2smcu_3D_term_r_p')
      call my_out_us(70,nCv,0,term_r_r,'b2smcu_3D_term_r_r')
      call my_out_us(70,nCv,0,term_t_p,'b2smcu_3D_term_t_p')
      call my_out_us(70,nCv,0,term_t_r,'b2smcu_3D_term_t_r')
      call my_out_us(70,nCv,0,hz,'b2smcu_3D_hz')
      call my_out_us(70,nCv,0,hx,'b2smcu_3D_hx')
      call my_out_us(70,nCv,0,hy,'b2smcu_3D_hy')
      call my_out_us(70,nCv,0,dhzdx,'b2smcu_3D_dhzdx')
      call my_out_us(70,nCv,0,dhxdx,'b2smcu_3D_dhxdx')
      call my_out_us(70,nCv,0,dhydx,'b2smcu_3D_dhydx')
      call my_out_us(70,nCv,0,dhzdy,'b2smcu_3D_dhzdy')
      call my_out_us(70,nCv,0,dhxdy,'b2smcu_3D_dhxdy')
      call my_out_us(70,nCv,0,dhydy,'b2smcu_3D_dhydy')
      call my_out_us(70,nCv,0,vs_mol,'b2smcu_3D_vs_mol')
      call my_out_us(70,nCv,0,termGx,'b2smcu_3D_termGx')
      call my_out_us(70,nCv,0,termGy,'b2smcu_3D_termGy')
      call my_out_us(70,nCv,0,dupdp,'b2smcu_3D_dupdp')
      call my_out_us(70,nCv,0,dupdr,'b2smcu_3D_dupdr')
      call my_out_us(70,nCv,0,durdp,'b2smcu_3D_durdp')
      call my_out_us(70,nCv,0,durdr,'b2smcu_3D_durdr')
      call my_out_us(70,nCv,0,dutdp,'b2smcu_3D_dutdp')
      call my_out_us(70,nCv,0,dutdr,'b2smcu_3D_dutdr')
      call my_out_us(70,nCv,0,u_pol,'b2smcu_3D_upol')
      call my_out_us(70,nCv,0,u_rad,'b2smcu_3D_urad')
      call my_out_us(70,nCv,0,u_tor,'b2smcu_3D_utor')


#endif
      endif                                                          

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sicf

      end subroutine b2smcu_3D
