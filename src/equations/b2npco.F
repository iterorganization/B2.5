*-----------------------------------------------------------------------
*.specification

      subroutine b2npco (nCv, nFc, nVx, nregionv, solving,               !srv 22.05.18
     &  solvcreg, ua_solve, solvmreg, itcnt, isb, rxf, dtim, switch,
     &  geo, mpg, pl, dv, sr, psnl)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_numerics_namelist                                        !srv 07.07.21 {
     & , only : dtco, time_factor                                        !srv 07.07.21 }
#ifndef FIXED_POINT
      use b2mod_balance
     & ,only : update_balance_sources_particles,
     &         update_balance_eirene_sources_particles,
     &         balance_netcdf,                                           !djm Jan2017
     &         tot_sna0to1, fna_tot                                      !djm Mar2019
#endif
      use b2mod_ad
     & , only : ncall_b2npco, rxg_npco
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2upco, ncall_b2usco, my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                               !srv 27.06.21
#endif
     &          ncall_b2ursc
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt, isb                        !srv 05.07.08
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Source), intent (inout) :: sr
      type (B2PlasmaSnapshot), intent (inout) :: psnl
      real (kind=R8) ::
     *  rxf, dtim
      logical solving, solvcreg(0:nregionv), ua_solve,                   !srv 22.05.18
     *  solvmreg(0:nregionv)
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPCO performs one iteration of a relaxation process for the
*     solution of the continuity equation for one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equation that
*     governs the particle density, nb, for a single atomic species.
*     An outline of the procedure follows.
*     ..compute residual
*     ..compute correction
*     ..apply correction
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      character chns*3                                             !srv 04.07.08
      real (kind=R8) ::
     *   aa(mpg%nCmxNv), wrk0(nCv), ts_factor
      logical warning
      integer :: iCv
*   ..procedures
      external xertst, ipgetr
      external b2xvsg, b2ursc, b2usco, b2upco
*   ..initialisation

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
*     pcm0, pcm1 - real.
*     pcm0 specifies an additional under-relaxation factor that is
*     applied to the velocity correction.
*     It is required that 0.le.pcm0.and.0.le.pcm1.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npco')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')
*   ..set internal parameters on first call
      if (ncall_b2npco.eq.0) then
        if(switch%b2mndt_style.eq.2) then
          rxg_npco = huge(1.0_R8)
          if (switch%b2npco_pcm0.gt.0.0_R8) then
            call xerrab('b2npco_pcm0 must be 0 for time-dependent mode')
          endif
          if (switch%b2npco_pcm1.gt.0.0_R8) then
            call xerrab('b2npco_pcm1 must be 0 for time-dependent mode')
          endif
          warning = .false.
          do iCv=1, nCv
            warning = warning .or.
     &        rxf*dtco(isb,mpg%cvReg(iCv))*time_factor(iCv).ne.1.0_R8
          enddo
          if (warning) then
            write(6,'(a,i3,a)')
     &        'rxf*dtco*time_factor for species ',isb,' should be 1 !'
            write(6,*) 'Time derivative terms may be inacurrate !'
          endif
        else
          rxg_npco = switch%b2npco_rxg
        endif
      endif
*   ..extensive tests on first few calls
      if (ncall_b2npco.lt.3) then
*    ..test sign of conb
       call b2xvsg (nFc, dv%conb(1,0,0), 1, 'conb00', '.ge.')
       call b2xvsg (nFc, dv%conb(1,0,1), 1, 'conb01', '.ge.')
       call b2xvsg (nFc, dv%conb(1,0,2), 1, 'conb02', '.ge.')
*    ..test sign of snbc, snbv
       call b2xvsg (nCv, sr%sna(1,0,isb), 1, 'snbc', '.ge.')
       call b2xvsg (nCv, sr%sna(1,1,isb), 1, 'snbv', '.le.')
*    ..test sign of pb, nb
       call b2xvsg (nCv, dv%pa(1,isb), 1, 'pb', '.gt.')
       call b2xvsg (nCv, pl%na(1,isb), 1, 'nb', '.gt.')
      endif
#ifndef FIXED_POINT
cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        tot_sna0to1(:,:,isb)=sr%sna(:,:,isb)    !djm Mar2019
        fna_tot(:,:,isb)=dv%fna(:,:,isb)        !djm Mar2019
        call update_balance_sources_particles(nCv,isb,pl%na(:,isb))
      endif
      if (balance_netcdf.ne.0 .or.
     &    switch%write_eirene_sources.ne.0) then
        call update_balance_eirene_sources_particles
     &   (nCv,isb,pl%na(:,isb))
      endif
#endif

      !Compute residual without time-stepping sources (needed at next time-step for ckn scheme)
      if (switch%b2mndt_ckn.eq.2) then
        if((solving.and.any(solvcreg(0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then      
          call b2ursc (nCv, nFc, mpg, pl%na(:,isb), sr%sna(:,:,isb),
     &      dv%fna(:,:,isb), dv%resco0(:,isb))
        endif
      endif

      if(switch%b2mndt_style.eq.2) then                                 !srv 07.07.21 {
*  .. calculate time-dependent source
        call b2scdt (nCv, mpg%nCi, isb, dtim,
     &   geo%cvVol, psnl%na(:,isb), switch, mpg, sr%snadt(1,0,isb))
*  .. add contribution from time-dependent source
        sr%sna(:,:,isb) = sr%sna(:,:,isb) + sr%snadt(:,:,isb)
      endif                                                             !srv 07.07.21 }

* ..particle balance
*   ..compute residual
      if((solving.and.any(solvcreg(0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then                   !srv 22.05.18
       if (switch%mdf_fnb.eq.0) then
        call b2ursc (nCv, nFc, mpg, pl%na(:,isb), sr%sna(:,:,isb),
     &    dv%fna(:,:,isb), dv%resco(:,isb))
       else
        call b2ursc (nCv, nFc, mpg, pl%na(:,isb), sr%sna(:,:,isb),
     &    dv%fna_mdf(:,:,isb), dv%resco(:,isb))
       endif
      endif                                                             !srv 22.05.18

*   ..modify residual according to time integration scheme
      if (switch%b2mndt_ckn.eq.2) then
        if((solving.and.any(solvcreg(0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then
          dv%resco(1:mpg%nCi,isb) = dv%resco(1:mpg%nCi,isb) +
     &      psnl%resco0(1:mpg%nCi,isb)
        endif
      elseif (switch%b2mndt_bdf.eq.2) then
        if((solving.and.any(solvcreg(0:mpg%nnreg(0))))
     &            .or.switch%get_residuals.eq.1) then
          dv%resco(1:mpg%nCi,isb) = dv%resco(1:mpg%nCi,isb) +
     &      0.5_R8*psnl%dnadt(1:mpg%nCi,isb) ! +1/2*V*(dn_i/dt) (where i denotes the previous time-step)
        endif
      endif

*   ..compute correction
      write (chns,'(i3.3)') isb
      if(solving.and.any(solvcreg(0:mpg%nnreg(0)))) then                !srv 22.05.18
        call b2usco (nCv, nFc, nVx, isb, switch, geo, mpg, nregionv,    !srv 13.10.06 26.09.12
     &   solvcreg, itcnt, rxg_npco, pl%na(:,isb), dv%pa(:,isb),
     &   sr%sna(:,1,isb), dv%flob, dv%conb,
     &   dv%resco(:,isb), dv%corpa(:,isb), aa, 'b2usco_'//chns)
*   ..apply correction
        call b2upco (nCv, nFc, nVx, geo, mpg, switch, nregionv,
     &   ua_solve, solvmreg, rxf, switch%b2npco_pcm0, dv%corpa(:,isb),
     &   dv%pcca(:,0,isb),
     &   dv%pa(:,isb), pl%na(:,isb), pl%ua(:,isb), isb, psnl)
      endif                                                             !srv 22.05.18

*    ..compute time derivative
      ! to obtain the true physical derivative V(dn/dt), we need to divide by ts_factor, which was previously multiplied to snadt inside b2scdt
      if (switch%b2mndt_ckn.eq.2) then
        ts_factor = 2.0_R8
      elseif (switch%b2mndt_bdf.eq.2) then
        ts_factor = 3.0_R8/2.0_R8
      else !default implicit Euler scheme
        ts_factor = 1.0_R8
      endif
      dv%dnadt(:,isb) = (-sr%snadt(:,1,isb)*pl%na(:,isb)
     &                  -sr%snadt(:,0,isb))/ts_factor               ! dnadt = V(dna/dt)

*    ..output quantities for continuity equation
      if (switch%b2npco_iout.ne.0) then                                 !srv 12.07.02 {
        write (chns,'(i3.3)') isb
        call my_out_us(70,nCv,0,dv%pa(1,isb),'b2npco_pa'//chns)
        call my_out_us(70,nCv,0,pl%na(1,isb),'b2npco_na'//chns)
        call my_out_us(70,nCv,0,pl%ua(1,isb),'b2npco_ua'//chns)
        call my_out_us(70,nCv,0,dv%resco(1,isb),'b2npco_resco'//chns)
        call my_out_us(70,nCv,0,sr%sna(1,0,isb),'b2npco_sna0'//chns)
        call my_out_us(70,nCv,0,sr%sna(1,1,isb),'b2npco_sna1'//chns)
        wrk0=sr%sna(:,0,isb)+pl%na(:,isb)*sr%sna(:,1,isb)
        call my_out_us(70,nCv,0,wrk0,'b2npco_sna'//chns)
        !wrk0=-sr%snadt(:,1,isb)*pl%na(:,isb)-sr%snadt(:,0,isb)
        call my_out_us(70,nCv,0,dv%dnadt(1,isb),'b2npco_dnadt'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2npco_fna_th'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2npco_fna_r'//chns)
        call my_out_us(70,nFc,1,dv%fna_mdf(1,0,isb),
     .    'b2npco_fna_mdf_th'//chns)
        call my_out_us(70,nFc,1,dv%fna_mdf(1,1,isb),
     .    'b2npco_fna_mdf_r'//chns)
      elseif (switch%iout_b2wdat.eq.4) then                             !iys 17.04.17
        write (chns,'(i3.3)') isb
        call my_out_us(70,nCv,0,pl%na(1,isb),'b2npco_na'//chns)
        call my_out_us(70,nCv,0,dv%resco(1,isb),'b2npco_resco'//chns)
        wrk0=sr%sna(:,0,isb)+pl%na(:,isb)*sr%sna(:,1,isb)
        call my_out_us(70,nCv,0,wrk0,'b2npco_sna'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2npco_fna_th'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2npco_fna_r'//chns)
        call my_out_us(70,nFc,1,dv%fna_mdf(1,0,isb),
     .    'b2npco_fna_mdf_th'//chns)
        call my_out_us(70,nFc,1,dv%fna_mdf(1,1,isb),
     .    'b2npco_fna_mdf_r'//chns)
      endif

* ..return
      ncall_b2npco = ncall_b2npco+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npco

      end subroutine b2npco
