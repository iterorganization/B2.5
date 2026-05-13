*-----------------------------------------------------------------------
*.specification

      subroutine b2sihs_mol (nCv, nFc, nVx, switch, geo, mpg,
     &  cvs_mol, cvshz_mol_eff, sh0_mol)
      use b2mod_types
      use b2mod_tallies
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa                 
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      use b2mod_geometry
      use b2mod_fluid_molecules
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx                          
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg

      real (kind=R8) :: cvshz_mol_eff(nFc,0:1), cvs_mol(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   sh0_mol(nCv,0:3)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIHS_MOL computes the source coefficients for a set of terms in the
*     fluid molecule heat equations.
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants

*   ..local variables
      integer iCv, iFc
      real (kind=R8), parameter :: eps = 1.0e-60_R8
      real (kind=R8) ::
     &  wrkc(nCv), wrkf(nFc,0:1), wrkf1(nFc,0:1),
     &  wrkv(nVx), visc(nCv,0:1)
      real (kind=R8) ::
     *  sh_du_mol(nCv,0:3),                          
     *  sh_va_mol(nCv,0:3), sh_vc_mol(nCv,0:3)
      real (kind=R8) :: t0, t1
      real (kind=R8) :: nm(nCv), um(nCv), tm(nCv)
*   ..procedures
      intrinsic min, max, sqrt
      external xertst, ipgeti, ipgetr, sfill
      external b2xvfv, b2xvsg

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sihs_mol')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')

* ..compute heat source terms
*   ..initialize sources to zero
      call sfill (nCv*4, 0.0_R8, sh0_mol, 1)
      call sfill (nCv*4, 0.0_R8, sh_du_mol, 1)
      call sfill (nCv*4, 0.0_R8, sh_va_mol, 1)
      call sfill (nCv*4, 0.0_R8, sh_vc_mol, 1)    
      um = u_mol(:,0)
      nm = n_mol(:,0)
      tm = t_mol(:,0)                      
                    

*   ..contribution from divergence um
      if(switch%b2sihs_style_int_vel.eq.0) then
        call intface(nCv,nFc,mpg%fcCv,geo%fcHc,um,wrkf(:,0))
      else
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,um,wrkf(:,0))
      endif
      wrkf(:,0) = wrkf(:,0)*geo%fcPbs
      call div0(nCv, nFc, mpg, wrkf(:,0), wrkc)
      do iCv = 1, mpg%nCi
        t0 = -wrkc(iCv)     
        sh_du_mol(iCv,0) = sh_du_mol(iCv,0)+
     &        (max(0.0_R8,t0)+switch%b2sihs__rf1*abs(t0))*nm(iCv)*
     &         tm(iCv)
        sh_du_mol(iCv,3) = sh_du_mol(iCv,3)+min(0.0_R8,t0)
     &                     -switch%b2sihs__rf1*abs(t0)
      enddo


*   ..contribution from viscous heating

      wrkf1(:,0) = cvshz_mol_eff(:,0)/((cvs_mol(:,0)+eps)*geo%fcHz)
      wrkf1(:,1) = cvshz_mol_eff(:,1)/((cvs_mol(:,1)+eps)*geo%fcHz)

      call diff(nCv, nFc, nVx, 0, geo, mpg, um, wrkv, wrkf)
      wrkf = cvs_mol(:,:)*wrkf1*wrkf**2/(geo%fcQalf+eps)
      call intcell(nFc, nCv, mpg, mpg%intcellP,wrkf(:,0),visc(:,0))
      call intcell(nFc, nCv, mpg, mpg%intcellR,wrkf(:,1),visc(:,1))

      do iCv = 1, mpg%nCi
        t0 = visc(iCv,0)
        t1 = visc(iCv,1)
        sh_vc_mol(iCv,0) = sh_vc_mol(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t0
        sh_vc_mol(iCv,3) = sh_vc_mol(iCv,3)-switch%b2sihs__rf2*t0/
     &                     (nm(iCv)*tm(iCv))
        sh_va_mol(iCv,0) = sh_va_mol(iCv,0)+
     &                     (1.0_R8+switch%b2sihs__rf2)*t1
        sh_va_mol(iCv,3) = sh_va_mol(iCv,3)-switch%b2sihs__rf2*t1/
     &                     (nm(iCv)*tm(iCv))                             
      enddo

      if (switch%b2sihs_mol.eq.1) then
        call b2saxpy (nCv*4, 1.0_R8, sh_du_mol , 1, sh0_mol, 1)
        call b2saxpy (nCv*4, 1.0_R8, sh_vc_mol , 1, sh0_mol, 1)
        call b2saxpy (nCv*4, 1.0_R8, sh_va_mol , 1, sh0_mol, 1)
      endif



      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sihs_mol

      end subroutine b2sihs_mol
