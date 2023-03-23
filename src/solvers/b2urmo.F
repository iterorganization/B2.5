*-----------------------------------------------------------------------
*.specification

      subroutine b2urmo (nCv, nFc, nVx, switch, geo, mpg,  isb, ub, rob, 
     .  smb, flub, flubv, cvsb, resmb, fmb, fkb, cvsb_eff)                      !srv 02.06.11
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
!WG_TODO      use b2mod_balance                                                 !djm Jan2017
!WG_TODO     & , only : fmo_flua, fmo_cvsa, balance_netcdf
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2urmo
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   ub(nCv), rob(nCv), smb(nCv,0:3), flub(nFc,0:1), flubv(nFc,0:1),
     *   cvsb(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resmb(nCv), fmb(nFc,0:1), fkb(nFc,0:1), cvsb_eff(nFc,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URMO computes the residual of the parallel momentum balance
*     equation for one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of the parallel momentum
*     balance equation for one species, integrated over the cells of
*     the B2 mesh. The equation has the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resmb = rhs-div(flux) .
*     The flux will be computed in the array fmb; it has the form:
*       fmbx = flubx*ub-cvsbx*(d/dx)ub
*       fmby = fluby*up-cvsby*(d/dy)ub
*     The right hand side has a linearised representation with
*     coefficients given by smb. This linearisation may include a
*     contribution from an implicit discretisation in time.
*
*     The input arguments smb, flub and cvsb contain appropriate
*     geometric coefficients so that in the present routine the
*     expression for the integrated residual on any cell does not
*     involve any geometric factors. However, the correct computation
*     of (radial) differences of ub across the faces requires knowledge
*     of geometric factors.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
c      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0             !srv 06.07.06
      character chns*3                                                  !srv 02.06.11
      real (kind=R8) ::
     *   flf(nFc,0:1), flv(nFc,0:1),
     *   wrk(nCv), wrkf(nFc,0:1), wrkv(nVx)
      real (kind=R8), parameter :: eps = 1.0e-60_R8
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, calcflow, calccoef
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2urmo')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2urmo.lt.3) then
*    ..test sign of rob, smb1, smb3
       call b2xvsg (nCv, rob, 1, 'rob', '.gt.')
       call b2xvsg (nCv, smb(1,1), 1, 'smb1', '.le.')
       call b2xvsg (nCv, smb(1,3), 1, 'smb3', '.le.')
*    ..test sign of cvsb
       wrkf(:,0) = cvsb(:,0)*geo%fcQalf(:,0)
       wrkf(:,1) = cvsb(:,1)*geo%fcQalf(:,1)
       call b2xvsg (nFc, wrkf(:,0), 1, 'cvsb0', '.ge.')
       call b2xvsg (nFc, wrkf(:,1), 1, 'cvsb1', '.ge.')
      endif

* ..compute flux
cwdk note: compared to original code, flv is now based on "effective" 
cwdk       viscosity, and will be zero if cvsb is zero
      call calcflow(nCv,nFc,nVx, switch%b2npmo_discr_meth,
     &                       geo,mpg,ub,flub,cvsb,fmb,flf,flv)
cwdk store effective viscosity, for use in viscous heating
      call diff(nCv, nFc, nVx, 0, geo, mpg, ub, wrkv, wrkf)
      cvsb_eff  = -flv/(wrkf + eps)

!WG_TODO      if (balance_netcdf.ne.0) then !djm Jan2017
!WG_TODO        fmo_cvsa(:,:,isb)=flv
!WG_TODO        fmo_flua(:,:,isb)=flf
!WG_TODO      endif

      if(switch%b2npmo_iout.eq.1 .or. switch%iout_b2wdat.eq.4) then
       write (chns,'(i3.3)') isb
       call my_out_us(70,nFc,1,flv(1,0),'b2urmo_etaPat_gradua_th'//chns)
       call my_out_us(70,nFc,1,flv(1,1),'b2urmo_etaPat_gradua_r'//chns)
      endif  

* ..compute kinetic energy flux
      ! Note: flf and flv used as dummy arrays in this call !!
      ! To be checked:
      ! - consistency in case of drifts (in particular: fna_cor?)
      ! - should contribution from flubv be in here?
      wrk  = 0.5_R8*ub**2
      wrkf = flub + flubv
      call calcflow(nCv,nFc,nVx,switch%b2npmo_discr_ meth,
     &                     geo,mpg,wrk,wrkf,cvsb,fkb,flf,flv)
      fkb(:,0) = fkb(:,0)/geo%fcHz
      fkb(:,1) = fkb(:,1)/geo%fcHz

* ..compute residual
*   ..contribution from right hand side
      resmb = smb(:,0) + smb(:,1)*ub + smb(:,2)*rob + smb(:,3)*rob*ub

*   ..contribution from fmb
      do iFc = 1, nFc
        resmb(mpg%fcCv(iFc,1)) = 
     &    resmb(mpg%fcCv(iFc,1)) - fmb(iFc,0) - fmb(iFc,1)
        resmb(mpg%fcCv(iFc,2)) = 
     &    resmb(mpg%fcCv(iFc,2)) + fmb(iFc,0) + fmb(iFc,1)
      enddo

* ..return
      ncall_b2urmo = ncall_b2urmo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2urmo

      end subroutine b2urmo
