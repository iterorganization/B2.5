*-----------------------------------------------------------------------
*.specification

      subroutine b2tinnt (nCv, nFc, nVx, ns, switch, geo, mpg, rt,
     &       facdrift, fac_ExB, pl, fna,
     &       csigin, fchin)                                              !srv 08.03.99 16.06.09
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef B25_EIRENE
      use eirmod_wneutrals
     & , only : dab2, rfluxa, pfluxa, tfluxa
#else
      use b2mod_eirdiag
     & , only : dab2, rfluxa, pfluxa, tfluxa
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer :: nCv, nFc, nVx, ns
      real (kind=R8) ::
     & facdrift(nFc), fac_ExB(nFc),
     & fna(nFc,0:1,0:ns-1),
     & csigin(nFc,0:1,0:ns-1,0:ns-1)
      type (B2Plasma), intent (in) :: pl
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Rates), intent (in) :: rt
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   fchin(nFc,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINNT computes the ion-neutral current.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, is0, iCv
#ifdef WG_TODO
      integer iamain
#ifdef B25_EIRENE
      integer ia, im
      character chns*3
#endif
#endif
      integer, save :: ncall = 0
      real (kind=R8), parameter :: eps = 1.0e-60_R8
      real (kind=R8) ::
     &   fnac(nCv,0:1,0:ns-1), vntrc(nCv,0:1),                            !srv 08.07.09
     &   vntrc_fluid(nCv,0:1), vntrc_kin(nCv,0:1),
     &   na_fluid(nCv,0:ns-1), na_kin(nCv,0:ns-1), na_tot(nCv,0:ns-1),
     &   vntrf(nFc,0:1), tif(nFc), naf(nFc,0:ns-1), dpo(nFc,0:1),
     &   fnan(nFc,0:1,0:ns-1),!WG_TODO work(-1:nx,-1:ny,0:1,0:ns-1),
     &   wrk(nCv), wrk0(nFc,0:1,0:ns-1), wrk1(nFc,0:1), wrk2(nFc),
     &   wrkVx(nVx), pit(nCv), pit2(nCv), vntrc_par(nCv)
      real (kind=R8) ::
#ifdef WG_TODO
     &   sigia(nCv), sec(nCv), csigia(nFc,0:1), sigi(nCv),               !srv 15.02.11
#endif
     &   facdriftm, fac_ExBm                                             !srv 20.10.17
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c      integer, save :: mode = 1
c      integer, save :: iout = 0
c      integer, save :: fchin_in_core = 0                                !srv 16.11.17
*   ..procedures
*   ..initialisation
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tinnt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgeti ('b2tinnt_iout', iout)
c        call ipgeti ('b2tral_mode', mode)
c        call ipgeti ('b2tinnt_fchin_in_core', fchin_in_core)            !srv 16.11.17
      endif
*   ..test input
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      fchin     = 0.0_R8                                                 !srv 09.07.99
      facdriftm = maxval(facdrift)                                       !srv 20.10.17
      fac_ExBm  = maxval(fac_ExB)

      if(switch%pot_eq.eq.1 .and.(facdriftm.gt.0.0_R8 .or.
     . fac_ExBm.gt.0.0_R8).and.switch%no_current.eq.0
     .            .and. switch%fch_ion_neutral .ne. 0.0_R8) then         !srv 26.07.06 01.07.99

*   ..extensive tests on first few calls
          if (ncall.lt.3) then
*    ..test sign of te
            call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
*    ..test sign of csigin
            do is = 0, ns-1
              do is0 = 0, ns-1
                wrk2 = csigin(:,0,is,is0)*geo%fcQalf(:,0)
                call b2xvsg (nFc, wrk2, 1, 'csigin0', '.ge.')
                wrk2 = csigin(:,1,is,is0)*geo%fcQalf(:,1)
                call b2xvsg (nFc, wrk2, 1, 'csigin1', '.ge.')
              enddo
            enddo
          endif

*    ..assemble densities (for fluid, kinetic, or hybrid fluid-kinetic neutrals)
          do is = 0, ns-1
            na_fluid(:,is) = pl%na(:,is)
            if (is_neutral(is).and.switch%use_eirene.ne.0) then
              na_kin(:,is) = (dab2(:,b2eatcr(is),1)+
     &                          switch%b2mndr_na_min)  !! -> to prevent dividing by zero for the velocities below
            else
              na_kin(:,is) = 0.0_R8
            endif
            na_tot(:,is) = na_kin(:,is) + na_fluid(:,is)
*   ..interpolate na to cell faces
            call intface(nCv, nFc, mpg%fcCv, geo%fcVol, na_tot(:,is),
     &          naf(:,is))
          enddo
*    ..computation of gradients and interpolations
          if(switch%b2tinnt_style.eq.0) then
          !! b2tinnt_style = 0
*   ..interpolate ti to cell faces
            call intface(nCv, nFc, mpg%fcCv, geo%fcVol, pl%ti, tif)
*   ..pre-computation of differences
            call diff(nCv, nFc, nVx, 0, geo, mpg, pl%po, wrkVx, dpo)
            do is = 0, ns-1
              wrk = log(na_tot(:,is))                                    !srv 12.05.18
              call diff(nCv, nFc, nVx, 0, geo, mpg, wrk, wrkVx,
     &          wrk0(:,:,is))
            enddo
            call diff(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkVx, wrk1)
          else
          !! b2tinnt_style = 1 (default)
            do is = 0, ns-1
*   ..pre-computation of differences
              wrk = pl%na(:,is)*pl%ti
              call diff(nCv, nFc, nVx, 0, geo, mpg, wrk,
     &          wrkVx, wrk0(:,:,is))
            enddo
            call diff(nCv, nFc, nVx, 0, geo, mpg, pl%po, wrkVx, dpo)
          endif
*   ..calculate neutral velocities
!! Remark: when the perpendicular momentum sources are scored in EIRENE and available in B2.5, this routine should be rewritten. This would be more self-consistent, whereas now
!! it is a mix of post-processed EIRENE atom velocities and B2.5 CX rates.
          do is0 = 0, ns-1
            if (is_neutral(is0).and.nint(zn(is0)).eq.1) then             !srv 20.10.17
*      ..get fluxes divided by area
              fnan(:,0,is0) = fna(:,0,is0)/(geo%fcS*geo%fcQalf(:,0)+eps)
              fnan(:,1,is0) = fna(:,1,is0)/(geo%fcS*geo%fcQalf(:,1)+eps)
*       ..get fluxes in the cell centre
              call intcell(nFc,nCv,mpg,mpg%intcellP,fnan(:,0,is0),
     &          fnac(:,0,is0))
              call intcell(nFc,nCv,mpg,mpg%intcellR,fnan(:,1,is0),
     &          fnac(:,1,is0))
*      ..get velocities in the cell centre
*        .. fluid velocities
              vntrc_fluid(:,0) = fnac(:,0,is0)/na_fluid(:,is0)  !! poloidal
              vntrc_fluid(:,1) = fnac(:,1,is0)/na_fluid(:,is0)  !! radial
              !! correct from poloidal to diamagnetic
              pit = abs(geo%cvBb(:,0))/geo%cvBb(:,3)  !! pitch, b_theta / B
              pit2 = abs(geo%cvBb(:,2))/geo%cvBb(:,3) !! b_phi / B
              vntrc_fluid(:,0) =(vntrc_fluid(:,0)-pl%ua(:,is0)*pit)/pit2
*        .. kinetic velocities
              if (switch%use_eirene.ne.0) then
                vntrc_kin(:,0) = pfluxa(:,b2eatcr(is0),1)/na_kin(:,is0) !! poloidal
                vntrc_kin(:,1) = rfluxa(:,b2eatcr(is0),1)/na_kin(:,is0) !! radial
              !! correct from poloidal to diamagnetic
                pit = abs(geo%cvBb(:,0))/geo%cvBb(:,3)  !! pitch, b_theta / B
                pit2 = abs(geo%cvBb(:,2))/geo%cvBb(:,3) !! b_phi / B
                vntrc_par = tfluxa(:,b2eatcr(is0),1)/na_kin(:,is0) !! approximate parallel velocity as toroidal velocity/b_phi
     &                       /pit2
                vntrc_kin(:,0) = (vntrc_kin(:,0)-vntrc_par*pit)/pit2
              else
                vntrc_kin(:,0) = 0.0_R8
                vntrc_kin(:,1) = 0.0_R8
                vntrc_par = 0.0_R8
              endif
!! density-weighted average of the velocities, needed for spatially hybrid fluid-kinetic neutrals
              do iCv=1,nCv
                vntrc(iCv,0) = (vntrc_fluid(iCv,0)*na_fluid(iCv,is0) +
     &                         vntrc_kin(iCv,0)*na_kin(iCv,is0))/
     &                         na_tot(iCv,is0)
                vntrc(iCv,1) = (vntrc_fluid(iCv,1)*na_fluid(iCv,is0) +
     &                         vntrc_kin(iCv,1)*na_kin(iCv,is0))/
     &                         na_tot(iCv,is0)
              enddo
*      ..interpolate velocities in the cell face
              call intface(nCv, nFc, mpg%fcCv, geo%fcVol, vntrc(:,0),
     &          vntrf(:,0))
              call intface(nCv, nFc, mpg%fcCv, geo%fcVol, vntrc(:,1),
     &          vntrf(:,1))
*      ..components of ion-neutral current
              do is = 0, ns-1
                if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then      !srv 20.10.17
                  wrk  = qe*rt%rza(:,is)
                  call intface(nCv,nFc,mpg%fcCv,geo%fcVol,wrk,wrk2)
                  if(switch%b2tinnt_style.eq.0) then
                    fchin(:,0) = fchin(:,0) + csigin(:,0,is,is0)*(
     &                (-dpo(:,0) - (wrk1(:,0) + tif*wrk0(:,0,is))/wrk2)*
     &                (geo%fcBb(:,2)/geo%fcBb(:,3))**2 +
     &                geo%fcBb(:,2)*vntrf(:,1)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                    fchin(:,1) = fchin(:,1) + csigin(:,1,is,is0)*(
     &                (-dpo(:,1) - (wrk1(:,1) + tif*wrk0(:,1,is))/wrk2)-
     &                geo%fcBb(:,3)*vntrf(:,0)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                  else
                    fchin(:,0) = fchin(:,0) + csigin(:,0,is,is0)*(
     &                (-dpo(:,0) - wrk0(:,0,is)/(naf(:,is)*wrk2))*
     &                (geo%fcBb(:,2)/geo%fcBb(:,3))**2 +
     &                geo%fcBb(:,2)*vntrf(:,1)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                    fchin(:,1) = fchin(:,1) + csigin(:,1,is,is0)*(
     &                (-dpo(:,1) - wrk0(:,1,is)/(naf(:,is)*wrk2))-
     &                geo%fcBb(:,3)*vntrf(:,0)*geo%fcQgam(:,0)*
     &                (geo%fcHc(:,1) + geo%fcHc(:,2)))
                  endif !! end b2tinnt_style
                endif !! end condition for ion species
              enddo !! end loop over ion species
            endif !! end condition for neutral species
          enddo  !! end loop over neutral species                                           !srv 08.07.09 }
      endif  !! end conditions to calculated i-n current

#ifdef B25_EIRENE
#ifdef WG_TODO
      if(iout.ne.0.and.switch%use_eirene.ne.0) then
       do ia = 1,natm
        write (chns,'(i3.3)') ia
        call my_out(70,nx,ny,rfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_rfluxa_eir'//chns)
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_pfluxa_eir'//chns)
        call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,ia,1),
     &                'b2tinnt_dab2_eir'//chns)
       enddo
       do im = 1,nmol
        write (chns,'(i3.3)') im
        call my_out(70,nx,ny,rfluxm(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_rfluxm_eir'//chns)
        call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_pfluxm_eir'//chns)
        call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,im,1),
     &              'b2tinnt_dmb2_eir'//chns)
       enddo
!       call my_out(70,nx,ny,dab2_,'b2tinnt_dab2_01')
       call my_out(70,nx,ny,sigia,'b2tinnt_sigia')
       call my_out(70,nx,ny,csigia(-1,-1,0,0),'b2tinnt_csigiax')
       call my_out(70,nx,ny,csigia(-1,-1,1,1),'b2tinnt_csigiay')
       call my_out(70,nx,ny,csigin(-1,-1,0,0,ismain,ismain0),
     &                                      'b2tinnt_csiginx')
       call my_out(70,nx,ny,csigin(-1,-1,1,1,ismain,ismain0),
     &                                      'b2tinnt_csiginy')
       call my_out(70,nx,ny,fnan(-1,-1,0,0,ismain0),'b2tinnt_fnanx')
       call my_out(70,nx,ny,fnan(-1,-1,1,1,ismain0),'b2tinnt_fnany')
       call my_out(70,nx,ny,fchin(-1,-1,0,0),'b2tinnt_fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1,1),'b2tinnt_fchiny')
!       call my_out(70,nx,ny,vntrc(-1,-1,0),'b2tinnt_vntrcx')
!       call my_out(70,nx,ny,vntrc(-1,-1,1),'b2tinnt_vntrcy')
      endif
#endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tinnt

      end subroutine b2tinnt
