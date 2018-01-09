*-----------------------------------------------------------------------
*.specification

      subroutine b2nph9 (nx, ny, ns, nregionv, itcnt, ismain, solving,  !srv 05.07.17
     &  solvereg, solvireg, solvpreg, solvmreg,
     &  rxf, BoRiS, fac_ExB, hx, hy, qz, qs, vol, bb, pbs, lnlam,
     &  floe, cone, floi, coni, ceqp, pccm, she, shi, fne, fch, cvsa,   !srv 28.06.03
     1  fchdia, uadia, wadia, vadia, vaecrb, na, ua, ue,                !srv 07.07.09
     &  te, ti, po, ne, ni, fhe, fhi, shef, shif, reshe, reshi, resht,  !srv 10.11.02
     &  corte, corti, ierr)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     & , only : b2npht_shei_bal, update_balance_sources_heat,
     &          balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv, itcnt, ismain                       !srv 05.07.17
      real (kind=R8) ::
     *  rxf, BoRiS, fac_ExB(-1:nx,-1:ny), vol(-1:nx,-1:ny),             !srv 07.06.00 02.06.99
     *  pbs(-1:nx,-1:ny,0:1),hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), 
     &  qz(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1), lnlam(-1:nx,-1:ny)    !srv 07.05.08 01.07.09
      real (kind=R8) ::
     *  floe(-1:nx,-1:ny,0:1,0:1), floi(-1:nx,-1:ny,0:1,0:1),
     &  cone(-1:nx,-1:ny,0:1,0:1,0:2), coni(-1:nx,-1:ny,0:1,0:1,0:2),   !srv 02.01.07
     &  ceqp(-1:nx,-1:ny), pccm(-1:nx,-1:ny,0:1,0:1),
     &  she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3),
     &  fhe(-1:nx,-1:ny,0:1,0:1), fhi(-1:nx,-1:ny,0:1,0:1),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), fne(-1:nx,-1:ny,0:1,0:1),
     &  fch(-1:nx,-1:ny,0:1,0:1), cvsa(-1:nx,-1:ny,0:1,0:1,0:ns-1),
     1  fchdia(-1:nx,-1:ny,0:1,0:1), uadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  bb(-1:nx,-1:ny,0:3), wadia(-1:nx,-1:ny,0:1,0:ns-1), 
     3  vadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1)   !srv 07.07.09
      logical solving(4),
     *  solvereg(0:nregionv), solvireg(0:nregionv),
     *  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  ua(-1:nx,-1:ny,0:ns-1), ue(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *   reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny), resht(-1:nx,-1:ny),
     &  corte(-1:nx,-1:ny), corti(-1:nx,-1:ny),
     &  shef(-1:nx,-1:ny,0:3), shif(-1:nx,-1:ny,0:3)                    !srv 10.11.02
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPH9 performs one iteration of a relaxation process for the
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
      real (kind=R8) ::
     *   rxg, pcm0, pcm1
      real (kind=R8) ::
     * aa9(-1:nx,-1:ny,-6:6), wrkf(-1:nx,-1:ny,0:1,0:1),                !srv 02.01.07
     1 wrkc(-1:nx,-1:ny,0:1,0:1,0:2), cortt(-1:nx,-1:ny),               !srv 02.01.07
     2 she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     3 wrk0(-1:nx,-1:ny), shei(-1:nx,-1:ny)                             !srv 18.09.09
      integer, save :: iout = 0, iout_b2wdat = 0                        !srv 17.07.05
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0             !srv 06.07.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0             !srv 06.07.06
      integer, save :: style = 1                                        !srv 09.07.08
      character chns*3                                                  !srv 17.07.05
      save ncall, rxg, pcm0, pcm1
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, ipgetr, lwget, b2saxpy
      external b2xvsg, b2xvff, b2xvfx, b2sihs, b2ursd, b2ush9, b2upht   !srv 02.01.07
*   ..initialization
      data ncall/0/, rxg/1.0e0_R8/, pcm0/1.0e0_R8/, pcm1/1.0e0_R8/

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
      call subini ('b2nph9')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2nph9_rxg', rxg)                                  !srv 02.01.07 {
       call ipgetr ('b2nph9_pcm0', pcm0)
       call ipgetr ('b2nph9_pcm1', pcm1)                                !srv 02.01.07 }
       call xertst (0.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                 !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)                 !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                 !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)                 !srv 06.07.06
       call ipgeti ('b2nph9_iout', iout)                                !srv 02.01.07
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       call ipgeti ('b2nph9_style', style)                              !srv 09.07.08
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of pbs
       call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbs')
       call b2xvfy (nx, ny, pbs(-1,-1,1), 'pbs')
*    ..test sign of transport coefficients
       if(no_Ptncr_x_he.eq.0 .and. no_Ptncr_y_he.eq.0) then             !srv 06.07.06 {
        call b2xvsg (n2, cone(-1,-1,0,0,0), 1, 'cone', '.ge.')          !srv 02.01.07
        call b2xvsg (n2, cone(-1,-1,1,1,0), 1, 'cone', '.ge.')
       endif                                                            !srv 06.07.06 }
       if(no_Ptncr_x_hi.eq.0 .and. no_Ptncr_y_hi.eq.0) then             !srv 06.07.06 {
       call b2xvsg (n2, coni(-1,-1,0,0,0), 1, 'coni', '.ge.')
       call b2xvsg (n2, coni(-1,-1,1,1,0), 1, 'coni', '.ge.')           !srv 02.01.07
       endif                                                            !srv 06.07.06 }
       call b2xvsg (n2, ceqp, 1, 'ceqp', '.ge.')
*    ..test edge values of transport coefficients
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0,0), 'cone0')                 !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,0,2), 'coni2')                 !srv 02.01.07 }
!       call b2xvfx (nx, ny, pccm(-1,-1,0,0), 'pccm0')                  !srv 15.04.05 {
!       call b2xvfx (nx, ny, pccm(-1,-1,1,1), 'pccm1')                  !srv 15.04.05 }
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
      endif

* ..main computation
*   ..compute final heat sources
      if(style.eq.0) then
       call b2sihs (nx, ny, ns, BoRiS, ismain, fac_ExB,
     &  hx, hy, qz, qs, vol, bb, pbs,                                   !xpb !srv 05.07.17
     &  lnlam, na, ua, ue, te, ti, po, ne, ni,                          !srv 01.07.09
     &  fch, cvsa, fchdia, uadia, wadia, vaecrb,                        !xpb
     &  she0, shi0)                                                     !xpb
      else
       call b2sihs_ (nx, ny, ns, BoRiS, ismain, fac_ExB,                !srv 07.06.00 02.06.99 28.06.03 05.07.17
     *  hx, hy, qz, qs, vol, bb, pbs, lnlam,
     *  na, ua, ue, te, ti, po, ne, ni, fne, fch, cvsa, she0, shi0)     !srv 07.06.00 02.06.99 07.07.09
      endif
      call b2saxpy (n2*4, 1.0e0_R8, she, 1, she0, 1)
      call b2saxpy (n2*4, 1.0e0_R8, shi, 1, shi0, 1)

cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        call update_balance_sources_heat(nx,ny,ne,ni,te,ti)
      endif

*   ..compute provisional reshe, reshi
*     (first ignore equipartition)
*   ..save final electron heat source ignoring equipartition            !srv 10.11.02 {
      shef = she0
*   ..save final ion heat source ignoring equipartition
      shif = shi0                                                       !srv 10.11.02 }
      call b2ursd (nx, ny, ne, te, she0, fhe, reshe)
      call b2ursd (nx, ny, ni, ti, shi0, fhi, reshi)
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

      if (balance_netcdf.ne.0) then
        b2npht_shei_bal(-1:nx,-1:ny) = shei(-1:nx,-1:ny)                !djm Jan2017
      endif

*   ..compute final resht, reshe, reshi
      do iy = -1, ny
       do ix = -1, nx
        resht(ix,iy) = reshe(ix,iy) + reshi(ix,iy)
        reshe(ix,iy) = reshe(ix,iy) - shei(ix,iy)
        reshi(ix,iy) = reshi(ix,iy) + shei(ix,iy)
       enddo
      enddo                                                             !srv 18.09.09 30.01.09 }

      if(iout .eq. 1 .or. iout_b2wdat .eq. 4) then                      !srv 17.07.05 01.11.10 {
        call my_out(70,nx,ny,shei,'b2nph9_shei')
!        do is = 0,ns-1
!         write (chns,'(i3.3)') is
!         call my_out(70,nx,ny,na(-1,-1,is),'b2nph9_na'//chns//'_input')
!         call my_out(70,nx,ny,ua(-1,-1,is),'b2nph9_ua'//chns//'_input')
!        enddo
!        call my_out(70,nx,ny,te,'b2nph9_te_input')
!        call my_out(70,nx,ny,ti,'b2nph9_ti_input')
!        call my_out(70,nx,ny,po,'b2nph9_po_input')
	call my_out(70,nx,ny,reshe,'b2nph9_reshe')
	call my_out(70,nx,ny,fhe(-1,-1,0,0),'b2nph9_fhex')
	call my_out(70,nx,ny,fhe(-1,-1,1,1),'b2nph9_fhey')
!	 call my_out(70,nx,ny,she0(-1,-1,0),'b2nph9_she0')
!	 call my_out(70,nx,ny,she0(-1,-1,1),'b2nph9_she1')
!	 call my_out(70,nx,ny,she0(-1,-1,2),'b2nph9_she2')
!	 call my_out(70,nx,ny,she0(-1,-1,3),'b2nph9_she3')
        wrk0=she0(:,:,0)+she0(:,:,1)*te+she0(:,:,2)*ne+
     &                   she0(:,:,3)*te*ne
	call my_out(70,nx,ny,wrk0,'b2nph9_she')
	call my_out(70,nx,ny,reshi,'b2nph9_reshi')
	call my_out(70,nx,ny,fhi(-1,-1,0,0),'b2nph9_fhix')
	call my_out(70,nx,ny,fhi(-1,-1,1,1),'b2nph9_fhiy')
!	 call my_out(70,nx,ny,shi0(-1,-1,0),'b2nph9_shi0')
!	 call my_out(70,nx,ny,shi0(-1,-1,1),'b2nph9_shi1')
!	 call my_out(70,nx,ny,shi0(-1,-1,2),'b2nph9_shi2')
!	 call my_out(70,nx,ny,shi0(-1,-1,3),'b2nph9_shi3')
        wrk0=shi0(:,:,0)+shi0(:,:,1)*ti+
     &       shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrk0,'b2nph9_shi')
      endif                                                             !srv 17.07.05 01.11.10 }
*   ..compute correction
      if (solving(4)) then                                              !srv 11.01.13
        call b2ush9 (nx, ny, ns, nregionv, itcnt,
     &   rxg, pcm1, pbs, ne, ni, te, ti,
     &   floe, cone, floi, coni, pccm(-1,-1,1,1), she0, shi0, ceqp,
     &   resht, reshe, reshi, cortt, corte, corti, wrkf, wrkc, aa9)
*   ..apply correction
        call b2upht (nx, ny, ns, nregionv, solving(3), solving(1),
     &   solvereg, solvireg, solvpreg, solvmreg,
     &   rxf, pcm0, cortt, corte, corti,
     &   pccm(-1,-1,0,0), na, ua, te, ti, po, ne, ni, wrk0, vol)
      endif
*   ..set error parameter and return
*     (at present all errors cause an abort through xerrab)
      ierr = 0

      if (iout.eq.1) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out(70,nx,ny,na(-1,-1,is),'b2nph9_na'//chns)
          call my_out(70,nx,ny,ua(-1,-1,is),'b2nph9_ua'//chns)
        end do
        call my_out(70,nx,ny,te,'b2nph9_te')
        call my_out(70,nx,ny,ti,'b2nph9_ti')
        call my_out(70,nx,ny,po,'b2nph9_po')
      elseif (iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,te,'b2nph9_te')
        call my_out(70,nx,ny,ti,'b2nph9_ti')
      end if

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nph9                                                            !srv 02.01.07 }

      end
