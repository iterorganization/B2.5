*-----------------------------------------------------------------------
*.specification

      subroutine b2npht (nx, ny, ns, nregionv, itcnt, ismain, solving,   !srv 05.07.17
     &  solvereg, solvireg, solvtreg, solvpreg, solvmreg,                !srv 22.05.18
     &  rxf, dtim, BoRiS, fac_ExB, hx, hy, qz, qc, qs, vol, gs, bb,      !srv 07.07.21
     &  pbs, lnlam, floe, cone, floi, coni, ceqp, pccm, she, shi,
     &  fne, fch, fchstoch, cvsa, cvsa_cl, cvsa_drho,                    !som 16.08.21 !srv 24.08.21
     &  fchdia, uadia, wadia, vadia, vaecrb,
     &  na, ua, ue, te, ti, po, ne, ni, fhe, fhi,                        !srv 09.07.08
     &  shef, shif, reshe, reshi, resht, corte, corti, ierr)
      use b2mod_types
      use b2mod_indirect
      use b2mod_tallies
      use b2mod_b2cmpa
      use b2mod_boundary_namelist
      use b2mod_balance !djm Jan2017
     & , only : b2npht_shei_bal, update_balance_sources_heat,
     &          balance_netcdf
      use b2mod_subsys
      use b2mod_sources
     & , only : shedt, shidt                                             !srv 31.01.21
      use b2mod_numerics_namelist
     & , only : dtee, dtei, time_factor
      use b2mod_plasma
     , , only : tep, tip, nep, nip, nap, kinrgyp
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv, itcnt, ismain                        !srv 05.07.17
      real (kind=R8) ::
     *  rxf, dtim, BoRiS, fac_ExB(-1:nx,-1:ny), vol(-1:nx,-1:ny),        !srv 07.07.21
     *  pbs(-1:nx,-1:ny,0:1), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     *  gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny,0:1),
     *  qz(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     *  lnlam(-1:nx,-1:ny)                                               !srv 01.07.09
      real (kind=R8) ::
     *  floe(-1:nx,-1:ny,0:1,0:1), floi(-1:nx,-1:ny,0:1,0:1),
     &  cone(-1:nx,-1:ny,0:1,0:1,0:2), coni(-1:nx,-1:ny,0:1,0:1,0:2),    !srv 02.01.07
     &  ceqp(-1:nx,-1:ny), pccm(-1:nx,-1:ny,0:1,0:1),
     &  she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *  fne(-1:nx,-1:ny,0:1,0:1), fchstoch(-1:nx,-1:ny,0:1),             !srv 09.07.08 24.08.21
     &  fch(-1:nx,-1:ny,0:1,0:1), cvsa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_cl(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     &  cvsa_drho(-1:nx,-1:ny,0:1,0:1,0:ns-1),                           !som 16.08.21
     1  fchdia(-1:nx,-1:ny,0:1,0:1), uadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     2  bb(-1:nx,-1:ny,0:3), wadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     3  vadia(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     4  vaecrb(-1:nx,-1:ny,0:1,0:1,0:ns-1)                               !srv 07.07.09
      logical solving(4),
     *  solvereg(0:nregionv), solvireg(0:nregionv),
     *  solvtreg(0:nregionv),                                            !srv 22.05.18
     *  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), ue(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     &  reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny), resht(-1:nx,-1:ny),
     &  corte(-1:nx,-1:ny), corti(-1:nx,-1:ny),
     &  shef(-1:nx,-1:ny,0:3), shif(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPHT performs one iteration of a relaxation process for the
*     solution of the equations for internal electron and atom energy.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equations that
*     govern the temperatures te and ti.
*     An outline of the procedure follows.
*     ..compute auxiliaries
*     ..compute residuals
*     ..compute corrections
*     ..apply corrections
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is
      integer, save :: iout = 0, iout_b2wdat = 0                       !srv 17.07.05
      integer, save :: style = 1                                       !srv 09.07.08
      integer, save :: get_residuals = 0                               !srv 22.05.18
      integer, save :: style_b2mndt = 1                                !srv 07.07.21
      character chns*3                                                 !srv 17.07.05
      logical warning_e, warning_i
      real (kind=R8) ::
     *  rxg, pcm0, pcm1, t0, ttim                                      !srv 07.07.21
      real (kind=R8) ::
     *  aa5(-1:nx,-1:ny,-2:2), wrkf(-1:nx,-1:ny,0:1,0:1),
     1  wrkc(-1:nx,-1:ny,0:1,0:1,0:2), cortt(-1:nx,-1:ny),                 !srv 25.11.08
     2  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     3  wrk0(-1:nx,-1:ny), shei(-1:nx,-1:ny)                           !srv 18.09.09
      save ncall, rxg, pcm0, pcm1
*   ..procedures
      intrinsic sqrt
      external xertst, ipgetr, b2saxpy
      external b2xvsg, b2xvff, b2xvfx, b2sihs, b2ursd, b2usht, b2upht
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     rxg - real.
*     rxg specifies a special under-relaxation parameter. A term of the
*     form "abs(residual)/rxg" is added to the diagonal of the matrix
*     of the correction equation, with the effect of limiting the
*     computed correction. rxg is dimensionless and of order unity;
*     smaller values of rxg imply stronger damping.
*     Not used when running in time-dependent mode.
*
*     pcm0 - real.
*     ...
*
*     pcm1 - real.
*     ...
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npht')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2mndt_style', style_b2mndt)                     !srv 07.07.21 {
       if(style_b2mndt.eq.2) then
        rxg = huge(1.0_R8)
        pcm0 = 0.0_R8
        pcm1 = 0.0_R8
        warning_e = .false.
        warning_i = .false.
        do iy = -1, ny
         do ix = -1, nx
           warning_e = warning_e .or.
     .      rxf*dtee(region(ix,iy,0))*time_factor(ix,iy).ne.1.0_R8
           warning_i = warning_i .or.
     .      rxf*dtei(region(ix,iy,0))*time_factor(ix,iy).ne.1.0_R8
         enddo
        enddo
        if (warning_e) write(6,*) 'rxf*dtee*time_factor should be 1 !'
        if (warning_i) write(6,*) 'rxf*dtei*time_factor should be 1 !'
        if (warning_e .or. warning_i)
     >   write(6,*) 'Time derivative terms may be inacurrate !'
       else
        rxg = 1.0_R8
        pcm0 = 1.0_R8
        pcm1 = 1.0_R8
        call ipgetr ('b2npht_rxg', rxg)
        call ipgetr ('b2npht_pcm0', pcm0)
        call ipgetr ('b2npht_pcm1', pcm1)
       endif                                                          !srv 07.07.21 }
       call xertst (0.0_R8.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.0_R8.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.0_R8.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2npht_iout', iout)                              !srv 02.01.07
       call ipgeti ('b2wdat_iout', iout_b2wdat)                       !iys 17.04.17
       call ipgeti ('b2npht_style', style)                            !srv 09.07.08
       call ipgeti ('get_residuals', get_residuals)                   !srv 22.05.18
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     & 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of transport coefficients
       call b2xvsg (n2, cone(-1,-1,0,0,0), 1, 'cone', '.ge.')
       call b2xvsg (n2, cone(-1,-1,1,1,0), 1, 'cone', '.ge.')
       call b2xvsg (n2, coni(-1,-1,0,0,0), 1, 'coni', '.ge.')
       call b2xvsg (n2, coni(-1,-1,1,1,0), 1, 'coni', '.ge.')
       call b2xvsg (n2, ceqp, 1, 'ceqp', '.ge.')
*    ..test edge values of transport coefficients
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0,0), 'cone0')               !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,0,2), 'coni2')               !srv 02.01.07 }
*    ..test sign of she, shi
       call b2xvsg (n2, she(-1,-1,0), 1, 'she0', '.ge.')
       call b2xvsg (n2, she(-1,-1,1), 1, 'she1', '.le.')
       call b2xvsg (n2, she(-1,-1,2), 1, 'she2', '.ge.')
       call b2xvsg (n2, she(-1,-1,3), 1, 'she3', '.le.')
       call b2xvsg (n2, shi(-1,-1,0), 1, 'shi0', '.ge.')
       call b2xvsg (n2, shi(-1,-1,1), 1, 'shi1', '.le.')
       call b2xvsg (n2, shi(-1,-1,2), 1, 'shi2', '.ge.')
       call b2xvsg (n2, shi(-1,-1,3), 1, 'shi3', '.le.')
*    ..test sign of te, ti
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test edge values of the fluxes
       call b2xvff (nx, ny, fhe, 'fhe')
       call b2xvff (nx, ny, fhi, 'fhi')
*    ..test fac_ExB
       call xertst (0.0_R8.le.minval(fac_ExB) .and.
     &       1.0_R8.ge.maxval(fac_ExB), 'faulty argument fac_ExB')
      endif

* ..diagnostic of sources (incoming)
      b2she=0.0_R8
      b2shi=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
            b2she(region(ix,iy,0))=b2she(region(ix,iy,0))+
     1       she(ix,iy,0)+she(ix,iy,1)*te(ix,iy)+
     1       she(ix,iy,2)*ne(ix,iy)+she(ix,iy,3)*ne(ix,iy)*te(ix,iy)
            b2shi(region(ix,iy,0))=b2shi(region(ix,iy,0))+
     1       shi(ix,iy,0)+shi(ix,iy,1)*ti(ix,iy)+
     1       shi(ix,iy,2)*ni(ix,iy,0)+
     1       shi(ix,iy,3)*ni(ix,iy,0)*ti(ix,iy)
          endif
        enddo
      enddo
* ..main computation
      tep = te                                                           !srv 07.07.21 {
      tip = ti
      if(style_b2mndt.eq.2) then
*  .. calculate time-dependent source
        call b2shdt (nx, ny, ns, dtim, BoRiS,
     &    vol, nap, tep, tip, ne, ni, kinrgyp, shedt, shidt, .false.)
*  .. add contribution from time-dependent source
        she = she+shedt
        she(:,:,0)=she(:,:,0)-1.5_R8*(ne-nep)*tep/dtim*vol
        shi = shi+shidt
        shi(:,:,0)=shi(:,:,0)-1.5_R8*(ni(:,:,0)-nip(:,:,0))*tip/dtim*vol
      endif                                                              !srv 207.07.21 }
*   ..compute final heat sources
      if(style.eq.0) then
       call b2sihs (nx, ny, ns, BoRiS, ismain, fac_ExB, hx, hy,          !xpb !srv 05.07.17
     &  vol, bb, pbs, lnlam, na, ua, ue, te, ti, po, ne, ni,             !srv 01.07.09
     &  fch, fchstoch, cvsa, fchdia, uadia, wadia, vaecrb, she0, shi0)   !xpb !srv 24.08.21
      else
       call b2sihs_ (nx, ny, ns, BoRiS, ismain, fac_ExB,                 !srv 07.06.00 02.06.99 28.06.03 05.07.17
     *  hx, hy, qz, qc, qs, vol, gs, bb, pbs, lnlam, na, ua, ue,
     *  te, ti, po, ne, ni, fne, fch, fchstoch, cvsa, cvsa_cl,           !som 16.08.21 !srv 24.08.21
     *  cvsa_drho, she0, shi0)                                           !som 16.08.21 !srv 07.06.00 02.06.99 07.07.09
      endif
csw 21feb2014 remove sources in ghostcells (SOLPS4)
      do ix=-1,nx
       do iy=-1,ny
        if(ifenebound(ix,iy) .ne. 0) then
         she0(ix,iy,:) = 0.0_R8
         shi0(ix,iy,:) = 0.0_R8
        endif
       enddo
      enddo

* ..diagnostic of sources (b2sihs)
      b2she0=0.0_R8
      b2shi0=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &     bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &     region(ix,iy,0).ne.0) then
            b2she0(region(ix,iy,0))=b2she0(region(ix,iy,0))+
     1       she0(ix,iy,0)+she0(ix,iy,1)*te(ix,iy)+
     1       she0(ix,iy,2)*ne(ix,iy)+she0(ix,iy,3)*ne(ix,iy)*te(ix,iy)
            b2shi0(region(ix,iy,0))=b2shi0(region(ix,iy,0))+
     1       shi0(ix,iy,0)+shi0(ix,iy,1)*ti(ix,iy)+
     1       shi0(ix,iy,2)*ni(ix,iy,0)+
     1       shi0(ix,iy,3)*ni(ix,iy,0)*ti(ix,iy)
          endif
        enddo
      enddo
      call b2saxpy (n2*4, 1.0_R8, she, 1, she0, 1)
      call b2saxpy (n2*4, 1.0_R8, shi, 1, shi0, 1)

cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        call update_balance_sources_heat(nx,ny,ne,ni,te,ti,ns)
      endif

      if (iout.eq.1 .or. iout_b2wdat.eq.4) then                         !srv 07.07.21
        call my_out(70,nx,ny,te,'b2npht_te_inp')
        call my_out(70,nx,ny,ti,'b2npht_ti_inp')
      endif                                                             !srv 07.07.21

*   ..compute provisional reshe, reshi
*     (first ignore equipartition)
*   ..save final electron heat source ignoring equipartition            !srv 10.11.02 {
      shef = she0
*   ..save final ion heat source ignoring equipartition
      shif = shi0                                                       !srv 10.11.02 }
      if ((solving(4).and.any(solvereg(0:nregionv))).or.
     &                    get_residuals.eq.1) then                      !srv 22.05.18
        call b2ursd (nx, ny, ne, te, she0, fhe, reshe)
      endif                                                             !srv 22.05.18
      if ((solving(4).and.any(solvireg(0:nregionv))).or.
     &                   get_residuals.eq.1) then                       !srv 22.05.18
        call b2ursd (nx, ny, ni, ti, shi0, fhi, reshi)
      endif                                                             !srv 22.05.18
*   ..compute heat source due to electron-ions collisions               !srv 18.09.09 30.01.10 {
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1 .or.
     &     bottomiy(ix,iy).eq.-2.or.topiy(ix,iy).eq.ny+1) then
         shei(ix,iy) = 0.0_R8
        else
         shei(ix,iy) = ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/
     &                (te(ix,iy)*sqrt(te(ix,iy)))*(te(ix,iy)-ti(ix,iy))
        endif
       enddo
      enddo

csw 21feb2014 remove source in ghostcells (SOLPS4)
      do ix=-1,nx
       do iy=-1,ny
        if(ifenebound(ix,iy) .ne. 0) then
         shei(ix,iy) = 0.0_R8
        endif
       enddo
      enddo

      if (balance_netcdf.ne.0) then
        b2npht_shei_bal(-1:nx,-1:ny) = shei(-1:nx,-1:ny) !djm Jan2017
      endif

*   ..compute final resht, reshe, reshi
      b2qie=0.0_R8
      do iy = -1, ny
       do ix = -1, nx
        b2qie(region(ix,iy,0))=b2qie(region(ix,iy,0))+shei(ix,iy)
       enddo
      enddo

      if ((solving(4).and.any(solvereg(0:nregionv)))                  !srv 22.05.18 {
     &    .or.get_residuals.eq.1) reshe = reshe - shei
      if ((solving(4).and.any(solvireg(0:nregionv)))
     &    .or.get_residuals.eq.1) reshi = reshi + shei
      if (solving(4).and.any(solvereg(0:nregionv)).and.               ! if total equation is not solved then resht = 0
     &   any(solvireg(0:nregionv)).and.any(solvtreg(0:nregionv)))
     &   resht = reshe + reshi                                        !srv 22.05.18 }

      if(iout .eq. 1 .or. iout_b2wdat.eq.4) then                      !srv 17.07.05 {
       call my_out(70,nx,ny,shei,'b2npht_shei')
       call my_out(70,nx,ny,te,'b2npht_te')
       call my_out(70,nx,ny,ti,'b2npht_ti')
       call my_out(70,nx,ny,reshe,'b2npht_reshe')
       call my_out(70,nx,ny,fhe(-1,-1,0,0),'b2npht_fhex')               !srv 10.03.10
       call my_out(70,nx,ny,fhe(-1,-1,1,1),'b2npht_fhey')               !srv 10.03.10
       call my_out(70,nx,ny,she0(-1,-1,0),'b2npht_she0')
       call my_out(70,nx,ny,she0(-1,-1,1),'b2npht_she1')
       call my_out(70,nx,ny,she0(-1,-1,2),'b2npht_she2')
       call my_out(70,nx,ny,she0(-1,-1,3),'b2npht_she3')
       wrk0=she0(:,:,0)+she0(:,:,1)*te+she0(:,:,2)*ne+
     &                  she0(:,:,3)*te*ne
       call my_out(70,nx,ny,wrk0,'b2npht_she')
       call my_out(70,nx,ny,reshi,'b2npht_reshi')
       call my_out(70,nx,ny,fhi(-1,-1,0,0),'b2npht_fhix')               !srv 10.03.10
       call my_out(70,nx,ny,fhi(-1,-1,1,1),'b2npht_fhiy')               !srv 10.03.10
       call my_out(70,nx,ny,shi0(-1,-1,0),'b2npht_shi0')
       call my_out(70,nx,ny,shi0(-1,-1,1),'b2npht_shi1')
       call my_out(70,nx,ny,shi0(-1,-1,2),'b2npht_shi2')
       call my_out(70,nx,ny,shi0(-1,-1,3),'b2npht_shi3')
       wrk0=shi0(:,:,0)+shi0(:,:,1)*ti+shi0(:,:,2)*ni(:,:,0)+
     &                  shi0(:,:,3)*ti*ni(:,:,0)
       call my_out(70,nx,ny,wrk0,'b2npht_shi')
      endif                                                           !srv 17.07.05 }
      if (solving(4)) then
*   ..compute correction
        call b2usht (nx, ny, ns, nregionv, itcnt,
     &   rxg, pcm1, pbs, ne, ni, te, ti,
     &   floe, cone, floi, coni, pccm(-1,-1,1,1), she0, shi0,
     &   ceqp, resht, reshe, reshi, cortt, corte, corti,
     &   wrkf, wrkc, aa5, 'b2npht_')                                  !srv 18.03.10
*   ..apply correction
        call b2upht (nx, ny, ns, nregionv, solving(3), solving(1),
     &   solvereg, solvireg, solvtreg, solvpreg, solvmreg,            !srv 22.05.18
     &   rxf, pcm0, cortt, corte, corti,
     &   pccm(-1,-1,0,0), na, ua, te, ti, po, ne, ni, wrk0, vol)
      endif
      if(iout .eq. 1) then                                            !srv 17.07.05 {
       do is = 0,ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,na(-1,-1,is),'b2npht_na'//chns)
        call my_out(70,nx,ny,ua(-1,-1,is),'b2npht_ua'//chns)
       enddo
       call my_out(70,nx,ny,te,'b2npht_te')
       call my_out(70,nx,ny,ti,'b2npht_ti')
       call my_out(70,nx,ny,po,'b2npht_po')
      elseif (iout_b2wdat.eq.4) then
       call my_out(70,nx,ny,te,'b2npht_te')
       call my_out(70,nx,ny,ti,'b2npht_ti')
       if(style_b2mndt.eq.2) then                                     !srv 07.07.21
         wrk0=1.5_R8*vol*ne*(te-tep)/dtim
         call my_out(70,nx,ny,wrk0,'b2npht_nedtedt')
         wrk0=1.5_R8*vol*tep*(ne-nep)/dtim
         call my_out(70,nx,ny,wrk0,'b2npht_tednedt')
         wrk0=1.5_R8*vol*ni(:,:,0)*(ti-tip)/dtim
         call my_out(70,nx,ny,wrk0,'b2npht_nidtidt')
         wrk0=1.5_R8*vol*tip*(ni(:,:,0)-nip(:,:,0))/dtim
         call my_out(70,nx,ny,wrk0,'b2npht_tidnidt')
       else
         do iy = -1,ny
          do ix = -1,nx
           ttim=dtim*dtee(region(ix,iy,0))*time_factor(ix,iy)
           t0 = 1.5_R8*vol(ix,iy)/ttim
           wrk0(ix,iy)=(ne(ix,iy)*te(ix,iy)-nep(ix,iy)*tep(ix,iy))*t0
          enddo
         enddo
         call my_out(70,nx,ny,wrk0,'b2npht_dnetedt')
         do iy = -1,ny
          do ix = -1,nx
           ttim=dtim*dtei(region(ix,iy,0))*time_factor(ix,iy)
           t0 = 1.5_R8*vol(ix,iy)/ttim
           wrk0(ix,iy)=(ni(ix,iy,0)*ti(ix,iy)-
     &                  nip(ix,iy,0)*tip(ix,iy))*t0
          enddo
         enddo
         call my_out(70,nx,ny,wrk0,'b2npht_dnitidt')
       endif                                                          !srv 07.07.21 }
      end if
*   ..set error parameter and return
*     (at present all errors cause an abort through xerrab)
      ierr = 0

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npht

      end subroutine b2npht
