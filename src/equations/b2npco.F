*-----------------------------------------------------------------------
*.specification

      subroutine b2npco (nCv, nFc, nVx, nregionv, solvcreg, ua_solve, 
     &  solvmreg, itcnt, isb, rxf, geo, mpg, pl, dv, sr, ierr)
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance
!WG_TODO     & ,only : update_balance_sources_particles, balance_netcdf      !djm Jan2017
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt, isb                       !srv 05.07.08
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Source), intent (in) :: sr
      real (kind=R8) ::
     *  rxf
      logical solvcreg(0:nregionv), ua_solve, solvmreg(0:nregionv)
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
*   ..output arguments (unspecified on entry)
      integer ierr
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
      integer ncall
      integer, save :: iout = 0, iout_b2wdat = 0                     !srv 04.07.08
      character chns*3                                               !srv 04.07.08
      real (kind=R8) ::
     *   rxg, pcm0, pcm1
      real (kind=R8) ::
     *   aa(mpg%nCmxNv), wrk0(nCv)
      save ncall, rxg, pcm0, pcm1
*   ..procedures
      external subini, subend, xertst, ipgetr, lwget
      external b2xvsg, b2xvff, b2ursc, b2usc9, b2upco
*   ..initialization
      data ncall/0/, rxg/1.0_R8/, pcm0/1.0_R8/, pcm1/1.0_R8/

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
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2npco_rxg', rxg)
       call ipgetr ('b2npco_pcm0', pcm0)
       call ipgetr ('b2npco_pcm1', pcm1)
       call xertst (0.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2npco_iout', iout)                              !srv 13.07.05
       call ipgeti ('b2wdat_iout', iout_b2wdat)                       !iys 17.04.17
      endif
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of conb
       call b2xvsg (nFc, dv%conb(1,0,0), 1, 'conb00', '.ge.')
       call b2xvsg (nFc, dv%conb(1,1,0), 1, 'conb10', '.ge.')
       call b2xvsg (nFc, dv%conb(1,0,1), 1, 'conb01', '.ge.')
       call b2xvsg (nFc, dv%conb(1,1,1), 1, 'conb11', '.ge.')
       call b2xvsg (nFc, dv%conb(1,0,2), 1, 'conb02', '.ge.')
       call b2xvsg (nFc, dv%conb(1,1,2), 1, 'conb12', '.ge.')
*    ..test sign of snbc, snbv
       call b2xvsg (nCv, sr%sna(1,0,isb), 1, 'snbc', '.ge.')
       call b2xvsg (nCv, sr%sna(1,1,isb), 1, 'snbv', '.le.')
*    ..test sign of pb, nb
       call b2xvsg (nCv, dv%pa(1,isb), 1, 'pb', '.gt.')
       call b2xvsg (nCv, pl%na(1,isb), 1, 'nb', '.gt.')
      endif
#ifdef WG_TODO
cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        call update_balance_sources_particles(nx,ny,isb,nb)
      endif
#endif
* ..particle balance
*   ..compute residual
      call b2ursc (nCv, nFc, geo, mpg, pl%na(1,isb), sr%sna(1,0,isb), 
     &  dv%fna(1,0,isb), dv%resco(1,isb))
*   ..compute correction
      write (chns,'(i3.3)') isb
      call b2usco (nCv, nFc, nVx, isb, geo, mpg, nregionv,            !srv 13.10.06 26.09.12
     &  solvcreg, itcnt, rxg, pcm1, pl%na(1,isb), dv%pa(1,isb), 
     &  sr%sna(1,1,isb), dv%flob, dv%conb, dv%pcca(1,1,isb), 
     &  dv%resco(1,isb), dv%corpa, aa, 'b2usco_'//chns)
*   ..apply correction
      call b2upco (nCv, nFc, nVx, geo, mpg, nregionv, ua_solve,
     &  solvmreg, rxf, pcm0, dv%corpa(1,isb), dv%pcca(1,0,isb), 
     &  dv%pa(1,isb), pl%na(1,isb), pl%ua(1,isb))
      if (iout.ne.0) then                                             !srv 12.07.02 {
        write (chns,'(i3.3)') isb
        call my_out_us(70,nCv,0,dv%pa(1,isb),'b2npco_pa'//chns)
        call my_out_us(70,nCv,0,pl%na(1,isb),'b2npco_na'//chns)
        call my_out_us(70,nCv,0,pl%ua(1,isb),'b2npco_ua'//chns)
        call my_out_us(70,nCv,0,dv%resco(1,isb),'b2npco_resco'//chns)
        call my_out_us(70,nCv,0,sr%sna(1,0,isb),'b2npco_sna0'//chns)
        call my_out_us(70,nCv,0,sr%sna(1,1,isb),'b2npco_sna1'//chns)
        wrk0=sr%sna(:,0,isb)+pl%na(:,isb)*sr%sna(:,1,isb)
        call my_out_us(70,nCv,0,wrk0,'b2npco_sna'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2npco_fnax'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2npco_fnay'//chns)
      elseif (iout_b2wdat.eq.4) then                                  !iys 17.04.17
        write (chns,'(i3.3)') isb
        call my_out_us(70,nCv,0,pl%na(1,isb),'b2npco_na'//chns)
        call my_out_us(70,nCv,0,dv%resco(1,isb),'b2npco_resco'//chns)
        wrk0=sr%sna(:,0,isb)+pl%na(:,isb)*sr%sna(:,1,isb)
        call my_out_us(70,nCv,0,wrk0,'b2npco_sna'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,0,isb),'b2npco_fnax'//chns)
        call my_out_us(70,nFc,1,dv%fna(1,1,isb),'b2npco_fnay'//chns)
      endif                                                          

*   ..set error parameter
*     (at present all errors cause an abort through xerrab)
      ierr = 0

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npco

      end
