*-----------------------------------------------------------------------
*.specification

      subroutine get_mult_dt(style, dtim, dv, geo, mpg, nCv, nFc, ns, 
     &                     te, ti, na, ua, te_min, te_max, mult_max, 
     &                     Cmax, mult_dt, dt_var4_ns,dt_var4_min,
     &                     nfilter)
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_plasma
      use b2us_map
      use b2mod_b2cmpa
      implicit none
*   ..input arguments
      type (B2Derivatives), intent (in) :: dv
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      integer, intent(in) :: style, nCv, nFc, ns, nfilter
      real (kind=R8), intent(in) ::
     *   dtim, te_min, te_max, mult_max, te(nCv), ti(nCv), 
     *   na(nCv,0:ns-1), ua(nCv,0:ns-1), Cmax
*   ..output arguments
      real (kind=R8), intent(out) ::
     *   mult_dt(nCv), dt_var4_ns(nCv,0:ns-1), dt_var4_min(nCv)
*-----------------------------------------------------------------------
*.documentation
*
*     get_mult_dt calculates a local (cell-based) multiplier to the timestep. 
*     The style switch determines the model.
*
*     style = 0: no multiplier, dt_mult = 1.0
*
*     style = 1: based on Te. In low Te regions (typically more interaction with neutrals) a small timestep is preferred. 
*                In high Te regions there is typically less neutral interaction, or this interaction
*                is a least less critical from a numerical point of view. 
*
*     style = 2: based on the magnitude of the sources compared to the particle/momentum/energy content of the cell. 
*
*     style = 3: CFL-like criterium for drift velocities
*
*     style = 4: Combination of option 2 and 3. The minimum of both options is taken.
*
*     style = 5: CFL-like criterium for drift velocities, based on final fluxes
*
*     Other styles can be added in the future.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: iCv, is, iFc, i1, i2, ifilter
      real (kind=R8) :: te_min_J, te_max_J, exp1, hx, hy, Vx, Vy, 
     &            Vx_a, Vy_a, Vx_e, Vy_e, dt3,
     &            dtx, dty, wrk0(nCv,0:ns-1), wrk1(nCv,0:ns-1),
     &            uadia_c0(nCv,0:ns-1), uadia_c1(nCv,0:ns-1),
     &            vedia_c0(nCv), vedia_c1(nCv), 
     &            veecrb_c0(nCv), veecrb_c1(nCv),
     &            wrk00(0:ns-1), wrk11(0:ns-1), dt4,
     &            wrk(nCv), mult_dt_tmp(nCv), min_mult
      real (kind=R8), parameter :: eps = 1.0e-20_R8
*   ..procedures
      external xertst
*-----------------------------------------------------------------------
*.computation  


      if (style.eq.0) then
        mult_dt = 1.0_R8
      elseif (style.eq.1) then
        te_min_J = te_min*ev !! convert from eV to joule
        te_max_J = te_max*ev
        do iCv = 1, nCv
          if (te(iCv).ge.te_max_J) then
            mult_dt(iCv) = mult_max
          elseif (te(iCv).le.te_min_J) then
            mult_dt(iCv) = 1.0_R8
          else
            exp1 = log10(mult_max)*(log10(te(iCv))-log10(te_min_J)) / 
     &                              (log10(te_max_J)-log10(te_min_J))
            mult_dt(iCv) = 10**exp1
          endif
        enddo
      elseif (style.eq.2) then !! This option is worked out in b2srdt itself.
        mult_dt = 1.0_R8 
      elseif (style.eq.3) then 
        !! interpolate drift velocities to cell centers
        call intcell1(nFc, nCv, ns-1, mpg, mpg%intcellP,
     &               dv%uadia(:,0,:),uadia_c0)
        call intcell1(nFc, nCv, ns-1, mpg, mpg%intcellR,
     &               dv%uadia(:,1,:),uadia_c1)
        call intcell(nFc, nCv, mpg, mpg%intcellP,
     &               dv%vedia(:,0),vedia_c0)
        call intcell(nFc, nCv, mpg, mpg%intcellR,
     &               dv%vedia(:,1),vedia_c1)
        call intcell(nFc, nCv, mpg, mpg%intcellP,
     &               dv%veecrb(:,0),veecrb_c0)
        call intcell(nFc, nCv, mpg, mpg%intcellR,
     &               dv%veecrb(:,1),veecrb_c1)
        do iCv = 1, nCv
          hy = geo%cvVol(iCv)/geo%cvHx(iCv)/geo%cvHz(iCv) !! hx and hy are approximations for non-quadrilateral cells
          hx = geo%cvHz(iCv)
          wrk00 = uadia_c0(iCv,:)+
     &            geo%cvBb(iCv,0)/geo%cvBb(iCv,3)*ua(iCv,:)
          wrk11 = uadia_c1(iCv,:)
          Vx_a = maxval(abs(wrk00))
          Vx_e = abs(vedia_c0(iCv) + veecrb_c0(iCv) + 
     &           geo%cvBb(iCv,0)/geo%cvBb(iCv,3)*dv%ue(iCv))
          Vy_a = maxval(abs(wrk11))
          Vy_e = abs(vedia_c1(iCv) + veecrb_c1(iCv))
          dtx = Cmax*hx/max(Vx_a,Vx_e)
          dty = Cmax*hy/max(Vy_a,Vy_e)
          dt3 = min(dtx,dty)
          dt3 = min(max(dtim,dt3),dtim*mult_max) !! limit with minimum (= baseline) and maximum (= baseline*mult_max) timestep
          mult_dt(iCv) = dt3/dtim
        enddo
        do ifilter = 1, nfilter
          mult_dt_tmp = mult_dt
          do iFc = 1, nFc
            i1 = mpg%fcCv(iFc,1)
            i2 = mpg%fcCv(iFc,2)
            min_mult = min(mult_dt_tmp(i1),mult_dt_tmp(i2))
            mult_dt(i1) = min_mult
            mult_dt(i2) = min_mult
          enddo
        enddo
      elseif (style.eq.4) then 
        mult_dt = 1.0_R8 !! For this option mult_dt is not used but other pre-compuations are performed here
        !! interpolate drift velocities to cell centers
        call intcell1(nFc, nCv, ns-1, mpg, mpg%intcellP,
     &               dv%uadia(:,0,:),wrk0)
        call intcell1(nFc, nCv, ns-1, mpg, mpg%intcellR,
     &               dv%uadia(:,1,:),wrk1)
        do iCv = 1, nCv
          hy = geo%cvVol(iCv)/geo%cvHx(iCv)/geo%cvHz(iCv) !! hx and hy are approximations for non-quadrilateral cells
          hx = geo%cvHz(iCv)
          wrk00 = wrk0(iCv,:)
          wrk11 = wrk1(iCv,:)
          do is=0,ns-1
            Vx = abs(wrk00(is))
            Vy = abs(wrk11(is))
            dtx = Cmax*hx/(Vx+0.001_R8)
            dty = Cmax*hy/(Vy+0.001_R8)
            dt4 = min(dtx,dty)
            dt_var4_ns(iCv,is) = dt4
          enddo
          dt_var4_min(iCv) = minval(dt_var4_ns(iCv,:))
        enddo

      else
        call xerrab('This b2mndt_var_style is not available')
      endif

*   ..return
      return

*-----------------------------------------------------------------------


      end subroutine get_mult_dt
