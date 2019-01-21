*-----------------------------------------------------------------------
*.specification

      subroutine b2npht (nCv, nFc, nVx, ns, geo, mpg, nregionv, itcnt,  !srv 05.07.17
     &  ismain, solving, solvereg, solvireg, solvpreg, solvmreg,
     &  rxf, BoRiS, pl, dv, co, sr, ierr)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2npht_shei_bal, update_balance_sources_heat,
!WG_TODO     &          balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nregionv, itcnt, ismain                !srv 05.07.17
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxf, BoRiS
      logical solving(4),
     *  solvereg(0:nregionv), solvireg(0:nregionv),
     *  solvpreg(0:nregionv), solvmreg(0:ns-1,0:nregionv)
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2Source), intent (inout) :: sr
*   ..output arguments (unspecified on entry)
      integer ierr
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
      integer ncall, iCv, is
      real (kind=R8) ::
     *   rxg, pcm0, pcm1
      real (kind=R8) ::
     * aa(mpg%nCmxNv), wrkf(nFc,0:1),                !srv 02.01.07
     1 wrkc(nFc,0:1,0:2), cortt(nCv),               !srv 02.01.07
     2 she0(nCv,0:3), shi0(nCv,0:3),
     3 wrk0(nCv), shei(nCv),                            !srv 18.09.09
     4 shef(nCv,0:3), shif(nCv,0:3)                     !srv 10.11.02
      integer, save :: iout = 0, iout_b2wdat = 0                        !srv 17.07.05
      integer, save :: mdf_fhe = 0, mdf_fhi = 0
      integer, save :: style = 1                                        !srv 09.07.08
      character chns*3                                                  !srv 17.07.05
      save ncall, rxg, pcm0, pcm1
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, ipgetr, lwget, b2saxpy
      external b2xvsg, b2xvff, b2xvfx, b2ursd, b2ush9, b2upht   !srv 02.01.07
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
      call subini ('b2npht')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2npht_rxg', rxg)                                  !srv 02.01.07 {
       call ipgetr ('b2npht_pcm0', pcm0)
       call ipgetr ('b2npht_pcm1', pcm1)                                !srv 02.01.07 }
       call xertst (0.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2npht_iout', iout)                                !srv 02.01.07
       call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
       call ipgeti ('b2npht_style', style)                              !srv 09.07.08
       call ipgeti ('b2tfhe_mdf', mdf_fhe)
       call ipgeti ('b2tfhi_mdf', mdf_fhi) 
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of transport coefficients
       call b2xvsg (nCv, co%ceqp, 1, 'ceqp', '.ge.')
*    ..test sign of she, shi
       call b2xvsg (nCv, sr%she(1,0), 1, 'she0', '.ge.')
       call b2xvsg (nCv, sr%she(1,1), 1, 'she1', '.le.')
       call b2xvsg (nCv, sr%she(1,2), 1, 'she2', '.ge.')
       call b2xvsg (nCv, sr%she(1,3), 1, 'she3', '.le.')
       call b2xvsg (nCv, sr%shi(1,0), 1, 'shi0', '.ge.')
       call b2xvsg (nCv, sr%shi(1,1), 1, 'shi1', '.le.')
       call b2xvsg (nCv, sr%shi(1,2), 1, 'shi2', '.ge.')
       call b2xvsg (nCv, sr%shi(1,3), 1, 'shi3', '.le.')
*    ..test sign of te, ti
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
      endif

* ..main computation
*   ..compute final heat sources
      if(style.eq.0) then
       call xerrab('b2npht -- style 0 not available for WG')
#ifdef WG_TODO
       call b2sihs (nCv, nFc, nVx, ns, geo, mpg, BoRiS, ismain, 
     &  dv%fac_ExB, dv%lnlam, pl%na, pl%ua, dv%ue, pl%te, pl%ti, pl%po, 
     &  dv%ne, dv%ni, dv%fch, co%cvsa, dv%fchdia, dv%uadia, dv%wadia, 
     &  dv%vaecrb, she0, shi0)
#endif
      else
       call b2sihs_ (nCv, nFc, nVx, ns, geo, mpg, BoRiS, ismain, 
     &  dv%fac_ExB, dv%lnlam, pl%na, pl%ua, dv%ue, pl%te, pl%ti, pl%po, 
     &  dv%ne, dv%ni, dv%fne, dv%fch, co%cvsa, she0, shi0)
      endif
      call b2saxpy (nCv*4, 1.0e0_R8, sr%she, 1, she0, 1)
      call b2saxpy (nCv*4, 1.0e0_R8, sr%shi, 1, shi0, 1)
#ifdef WG_TODO
cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        call update_balance_sources_heat(nx,ny,ne,ni,te,ti)
      endif
#endif
*   ..compute provisional reshe, reshi
*     (first ignore equipartition)
*   ..save final electron heat source ignoring equipartition            !srv 10.11.02 {
      shef = she0
*   ..save final ion heat source ignoring equipartition
      shif = shi0                                                       !srv 10.11.02 }
      if (mdf_fhe.eq.0) then
        call b2ursd (nCv, nFc, geo, mpg, dv%ne, pl%te, she0, dv%fhe,
     &     dv%reshe)
      else
        call b2ursd (nCv, nFc, geo, mpg, dv%ne, pl%te, she0, dv%fhe_mdf,
     &     dv%reshe)
      endif
      if (mdf_fhi.eq.0) then
        call b2ursd (nCv, nFc, geo, mpg, dv%ni, pl%ti, shi0, dv%fhi, 
     &     dv%reshi)
      else
        call b2ursd (nCv, nFc, geo, mpg, dv%ni, pl%ti, shi0, dv%fhi_mdf, 
     &     dv%reshi)
      endif
      do iCv = 1, mpg%nCi
        shei(iCv) = co%ceqp(iCv)*dv%ne(iCv)*dv%ni(iCv,0)/
     &     (pl%te(iCv)*sqrt(pl%te(iCv)))*(pl%te(iCv)-pl%ti(iCv))
      enddo
      do iCv = mpg%nCi+1, nCv
        shei(iCv) = 0.0_R8
      enddo
#ifdef WG_TODO
      if (balance_netcdf.ne.0) then
        b2npht_shei_bal(-1:nx,-1:ny) = shei(-1:nx,-1:ny)                !djm Jan2017
      endif
#endif
*   ..compute final resht, reshe, reshi
      dv%resht = dv%reshe + dv%reshi
      dv%reshe = dv%reshe - shei
      dv%reshi = dv%reshi + shei

      if(iout .eq. 1 .or. iout_b2wdat .eq. 4) then                      !srv 17.07.05 01.11.10 {
        call my_out_us(70,nCv,0,shei,'b2npht_shei')
        call my_out_us(70,nCv,0,co%ceqp,'b2npht_ceqp')
!        do is = 0,ns-1
!         write (chns,'(i3.3)') is
!         call my_out_us(70,nCv,0,na(1,is),'b2npht_na'//chns//'_input')
!         call my_out_us(70,nCv,0,ua(1,is),'b2npht_ua'//chns//'_input')
!        enddo
!        call my_out_us(70,nCv,0,te,'b2npht_te_input')
!        call my_out_us(70,nCv,0,ti,'b2npht_ti_input')
!        call my_out_us(70,nCv,0,po,'b2npht_po_input')
	call my_out_us(70,nCv,0,dv%reshe,'b2npht_reshe')
	call my_out_us(70,nFc,1,dv%fhe(1,0),'b2npht_fhex')
	call my_out_us(70,nFc,1,dv%fhe(1,1),'b2npht_fhey')
!	 call my_out_us(70,nCv,0,she0(1,0),'b2npht_she0')
!	 call my_out_us(70,nCv,0,she0(1,1),'b2npht_she1')
!	 call my_out_us(70,nCv,0,she0(1,2),'b2npht_she2')
!	 call my_out_us(70,nCv,0,she0(1,3),'b2npht_she3')
        wrk0=she0(:,0)+she0(:,1)*pl%te+she0(:,2)*dv%ne+
     &       she0(:,3)*pl%te*dv%ne
	call my_out_us(70,nCv,0,wrk0,'b2npht_she')
	call my_out_us(70,nCv,0,dv%reshi,'b2npht_reshi')
	call my_out_us(70,nFc,1,dv%fhi(1,0),'b2npht_fhix')
	call my_out_us(70,nFc,1,dv%fhi(1,1),'b2npht_fhiy')
!	 call my_out_us(70,nCv,0,shi0(1,0),'b2npht_shi0')
!	 call my_out_us(70,nCv,0,shi0(1,1),'b2npht_shi1')
!	 call my_out_us(70,nCv,0,shi0(1,2),'b2npht_shi2')
!	 call my_out_us(70,nCv,0,shi0(1,3),'b2npht_shi3')
        wrk0=shi0(:,0)+shi0(:,1)*pl%ti+
     &       shi0(:,2)*dv%ni(:,0)+shi0(:,3)*pl%ti*dv%ni(:,0)
        call my_out_us(70,nCv,0,wrk0,'b2npht_shi')
      endif                                                             !srv 17.07.05 01.11.10 }
*   ..compute correction
      if (solving(4)) then                                              !srv 11.01.13
        call b2usht (nCv, nFc, nVx, ns, geo, mpg, nregionv, itcnt,
     &  rxg, pcm1, dv%ne, dv%ni, pl%te, pl%ti,
     &  dv%floe, dv%cone, dv%floi, dv%coni, dv%pccm(1,1), she0, shi0, 
     &  co%ceqp, dv%resht, dv%reshe, dv%reshi, dv%cortt, dv%corte, 
     &  dv%corti, wrkf, wrkc, aa)
*   ..apply correction
        call b2upht (nCv, nFc, nVx, ns, geo, mpg, nregionv, 
     &  solving(3), solving(1), solvereg, solvireg, solvpreg, solvmreg,
     &  rxf, pcm0, dv%cortt, dv%corte, dv%corti,
     &  dv%pccm(1,0), pl%na, pl%ua, pl%te, pl%ti, pl%po, dv%ne, dv%ni)
      endif
*   ..set error parameter and return
*     (at present all errors cause an abort through xerrab)
      ierr = 0

      if (iout.eq.1) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,pl%na(1,is),'b2npht_na'//chns)
          call my_out_us(70,nCv,0,pl%ua(1,is),'b2npht_ua'//chns)
        end do
        call my_out_us(70,nCv,0,pl%te,'b2npht_te')
        call my_out_us(70,nCv,0,pl%ti,'b2npht_ti')
        call my_out_us(70,nCv,0,pl%po,'b2npht_po')
        call my_out_us(70,nCv,0,dv%corte,'b2npht_corte')
        call my_out_us(70,nCv,0,dv%corti,'b2npht_corti')
        call my_out_us(70,nCv,0,dv%cortt,'b2npht_cortt')
      elseif (iout_b2wdat.eq.4) then
        call my_out_us(70,nCv,0,pl%te,'b2npht_te')
        call my_out_us(70,nCv,0,pl%ti,'b2npht_ti')
      end if

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npht                                                            !srv 02.01.07 }

      end
