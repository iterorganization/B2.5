*-----------------------------------------------------------------------
*.specification

      subroutine b2npc7 (nx, ny, fac_vis, guard_flows, zb,              !srv 05.07.08 02.04.18 {
     &  nregionv, solving, solvcreg, ua_solve, solvmreg, itcnt, isb,    !srv 22.05.18
     &  rxf, dtim, vol, pbs, flob, conb, snb,                           !srv 07.07.21
     &  pccb, ti, pb, nb, ub, fnb, rescb, corpb, ierr)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_balance
     & ,only : update_balance_sources_particles,                        !djm Jan2017
     &         update_balance_eirene_sources_particles,
     &         balance_netcdf,                                          !djm Jan2017
     &         tot_sna0to1, fna_tot                                     !djm Mar2019
      use b2mod_subsys
      use b2mod_sources
     , , only : snadt                                                   !srv 07.07.21 {
      use b2mod_numerics_namelist
     & , only : dtco, time_factor
       use b2mod_plasma
     & , only : nap                                                     !srv 07.07.21 }
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, guard_flows, nregionv, itcnt, isb                 !srv 05.07.08
      real (kind=R8) ::
     *  fac_vis(-1:nx,-1:ny), rxf, dtim, pbs(-1:nx,-1:ny,0:1),          !srv 07.07.21
     &  flob(-1:nx,-1:ny,0:1,0:1), conb(-1:nx,-1:ny,0:1,0:1,0:4),       !srv 13.10.06
     &  snb(-1:nx,-1:ny,0:1), pccb(-1:nx,-1:ny,0:1,0:1),
     &  vol(-1:nx,-1:ny), ti(-1:nx,-1:ny), zb(-1:nx,-1:ny)
      logical solving, solvcreg(0:nregionv), ua_solve,                  !srv 22.05.18
     &  solvmreg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  pb(-1:nx,-1:ny), nb(-1:nx,-1:ny), ub(-1:nx,-1:ny),
     &  fnb(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *  rescb(-1:nx,-1:ny), corpb(-1:nx,-1:ny)
*   ..common blocks


*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPC7 performs one iteration of a relaxation process for the
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
      integer ncall, n2, ix, iy, use_eirene
      integer, save :: iout = 0, iout_b2wdat = 0                   !srv 04.07.08
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0        !srv 06.07.06
      integer, save :: get_residuals = 0                           !srv 22.05.18
      integer, save :: style_b2mndt = 1                            !srv 07.07.21
      integer, save :: write_eirene_sources = 0
      character chns*3                                             !srv 04.07.08
      logical warning
#ifdef _OPENMP
!$OMP THREADPRIVATE(iout, iout_b2wdat)
!$OMP THREADPRIVATE(no_Ptncr_x_co, no_Ptncr_y_co)
#endif
      real (kind=R8) ::
     *  rxg, pcm0, pcm1, ttim, t0                                  !srv 07.07.21
      real (kind=R8) ::
     *  aa7(-1:nx,-1:ny,-3:3), wrk0(-1:nx,-1:ny)
      save ncall, rxg, pcm0, pcm1, use_eirene
#ifdef _OPENMP
!$OMP THREADPRIVATE(ncall, rxg, pcm0, pcm1, use_eirene)
#endif
*   ..procedures
      external xertst, ipgetr, ipgeti
      external b2xvsg, b2xvff, b2ursc, b2usco, b2upco
*   ..initialisation
      data ncall/0/, use_eirene/0/

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
*     Not used if running in time-dependent mode.
*
*     pcm0, pcm1 - real.
*     pcm0 specifies an additional under-relaxation factor that is
*     applied to the velocity correction.
*     It is required that 0.le.pcm0.and.0.le.pcm1.
*     Not used if running in time-dependent mode.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npc7')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2mndt_style', style_b2mndt)                       !srv 07.07.21 {
       if(style_b2mndt.eq.2) then
        rxg = huge(1.0_R8)
        pcm0 = 0.0_R8
        pcm1 = 0.0_R8
        warning = .false.
        do iy = -1, ny
         do ix = -1, nx
          warning = warning .or.
     .     rxf*dtco(isb,region(ix,iy,0))*time_factor(ix,iy).ne.1.0_R8
         enddo
        enddo
        if (warning) then
          write(6,'(a,i3,a)')
     .     'rxf*dtco*time_factor for species ',isb,' should be 1 !'
          write(6,*) 'Time derivative terms may be inacurrate !'
        end if
       else
        rxg = 1.0_R8
        pcm0 = 1.0_R8
        pcm1 = 1.0_R8
        call ipgetr ('b2npc7_rxg', rxg)
        call ipgetr ('b2npc7_pcm0', pcm0)
        call ipgetr ('b2npc7_pcm1', pcm1)
       endif                                                          !srv 07.07.21 }
       call xertst (0.0_R8.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.0_R8.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.0_R8.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)               !srv 06.07.06
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)               !srv 06.07.06
       call ipgeti ('b2npc7_iout', iout)                              !srv 04.07.08
       call ipgeti ('b2wdat_iout', iout_b2wdat)                       !iys 17.04.17
       call ipgeti ('b2mndr_eirene', use_eirene)
       call ipgeti ('get_residuals', get_residuals)                   !srv 22.05.18
       call ipgeti ('cdfmovie_eirene_sources',write_eirene_sources)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then

*    ..test sign of conb
       if(no_Ptncr_x_co.eq.0 .and. no_Ptncr_y_co.eq.0) then           !srv 06.07.06 {
        call b2xvsg (n2, conb(-1,-1,0,0,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,2), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,2), 1, 'conb', '.ge.')
       endif                                                          !srv 06.07.06 }
*    ..test edge values of flob, conb
       call b2xvff (nx, ny, flob, 'flob')
       call b2xvff (nx, ny, conb(-1,-1,0,0,0), 'conb0')               !srv 13.10.06 {
       call b2xvff (nx, ny, conb(-1,-1,0,0,1), 'conb1')
       call b2xvff (nx, ny, conb(-1,-1,0,0,2), 'conb2')               !srv 13.10.06 }
*    ..test sign of snbc, snbv
       call b2xvsg (n2, snb(-1,-1,0), 1, 'snbc', '.ge.')
       call b2xvsg (n2, snb(-1,-1,1), 1, 'snbv', '.le.')

*    ..test sign of pb, nb
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')

*    ..test edge values of fnb
       call b2xvff (nx, ny, fnb, 'fnb')
      endif

      if(style_b2mndt.eq.2) then                                        !srv 07.07.21 {
*  .. calculate time-dependent source
        call b2scdt (nx, ny, isb, guard_flows, dtim,
     &   vol, nap(-1,-1,isb), snadt(-1,-1,0,isb), .false.)
*  .. add contribution from time-dependent source
        snb(:,:,:) = snb(:,:,:) + snadt(:,:,:,isb)
      endif                                                             !srv 24.05.21 }

cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        tot_sna0to1(:,:,:,isb)=snb !djm Mar2019
        fna_tot(:,:,:,:,isb)=fnb   !djm Mar2019
        call update_balance_sources_particles(nx,ny,isb,nb)
      endif
      if (balance_netcdf.ne.0 .or. write_eirene_sources.ne.0) then
        call update_balance_eirene_sources_particles(nx,ny,isb,nb)
      endif

* ..particle balance
*   ..compute residual
      if((solving.and.any(solvcreg(0:nnreg(0))))
     &            .or.get_residuals.eq.1) then                          !srv 22.05.18
       call b2ursc (nx, ny, nb, snb, fnb, rescb)
       if(guard_flows.eq.2 .or.                                          !srv 06.07.08 {
     &    (is_neutral(isb) .and. use_eirene.ne.0)) then
        do iy = -1, ny, ny+1                                             !srv 27.06.08 {
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           rescb(ix,iy) = 0.0_R8
          endif
         enddo
        enddo                                                            !srv 27.06.08 }
       endif
      endif                                                             !srv 22.05.18
csw 11dec2012, check for solv flags
      write (chns,'(i3.3)') isb
      if(solving.and.any(solvcreg(0:nnreg(0)))) then                    !srv 22.05.18
*   ..compute correction
       call b2usc7 (nx, ny, fac_vis, isb, zb, nregionv, solvcreg,
     &  itcnt, rxg, pcm1, vol, pbs, ti, nb, pb, snb(-1,-1,1), flob,
     &  conb, pccb(-1,-1,1,1), rescb, corpb, aa7, 'b2usc7_'//chns)      !srv 02.04.18
*   ..apply correction
       call b2upco (nx, ny, nregionv, ua_solve, solvmreg, rxf, pcm0,
     &  vol, corpb, pccb(-1,-1,0,0), pb, nb, ub, isb)                   ! IYS 19.12.2017
       if(guard_flows.eq.2 .or.                                         !srv 06.07.08 {
     &    (is_neutral(isb) .and. use_eirene.ne.0)) then
*   ..replace the incorrect corner values by the interpolated values
        iy=-1                                                           !srv 27.06.08 {
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2)
     &     nb(ix,iy) = exp(
     .      log(nb(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(nb(topix(ix,iy),topiy(ix,iy))) -
     &      log(nb(topix(rightix(ix,iy),rightiy(ix,iy)),
     &             topiy(rightix(ix,iy),rightiy(ix,iy)))))
          if(rightix(ix,iy).eq.nx+1)
     &     nb(ix,iy) = exp(
     .      log(nb(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(nb(topix(ix,iy),topiy(ix,iy))) -
     &      log(nb(topix(leftix(ix,iy),leftiy(ix,iy)),
     &             topiy(leftix(ix,iy),leftiy(ix,iy)))))
        enddo
        iy=ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2)
     &     nb(ix,iy) = exp(
     .      log(nb(rightix(ix,iy),rightiy(ix,iy))) +
     .      log(nb(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(nb(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &             bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
          if(rightix(ix,iy).eq.nx+1)
     &     nb(ix,iy) = exp(
     .      log(nb(leftix(ix,iy),leftiy(ix,iy))) +
     .      log(nb(bottomix(ix,iy),bottomiy(ix,iy))) -
     &      log(nb(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &             bottomiy(leftix(ix,iy),leftiy(ix,iy)))))
        enddo                                                           !srv 27.06.08 }
       endif                                                            !srv 06.07.08 }
      endif

      if (iout.ne.0) then                                            !srv 12.07.02 {
       if (in_parallel()) then
#ifdef _OPENMP
!$OMP MASTER
#endif
         write(*,*) 'B2NPC7 will not write to file in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
       else
        call my_out(70,nx,ny,pb,'b2npc7_pa'//chns)
        call my_out(70,nx,ny,nb,'b2npc7_na'//chns)
        call my_out(70,nx,ny,ub,'b2npc7_ua'//chns)
        call my_out(70,nx,ny,rescb,'b2npc7_resco'//chns)
        call my_out(70,nx,ny,snb(-1,-1,0),'b2npc7_sna0'//chns)
        call my_out(70,nx,ny,snb(-1,-1,1),'b2npc7_sna1'//chns)
        wrk0=snb(:,:,0)+nb*snb(:,:,1)
        call my_out(70,nx,ny,wrk0,'b2npc7_sna'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,0),'b2npc7_fnax'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,1),'b2npc7_fnay'//chns)
        wrk0=-snadt(:,:,1,isb)*nb(:,:)-snadt(:,:,0,isb)              !srv 09.03.21
        call my_out(70,nx,ny,wrk0,'b2npc7_dnadt'//chns)              !srv 09.03.21
       endif
      endif                                                          !srv 12.07.02 }

      if (iout_b2wdat.eq.4) then
       if (in_parallel()) then
#ifdef _OPENMP
!$OMP MASTER
#endif
         write(*,*) 'B2NPC7 will not write to file in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
       else
        call my_out(70,nx,ny,nb,'b2npc7_na'//chns)
        call my_out(70,nx,ny,rescb,'b2npc7_resco'//chns)
        wrk0=snb(:,:,0)+nb*snb(:,:,1)
        call my_out(70,nx,ny,wrk0,'b2npc7_sna'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,0),'b2npc7_fnax'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,1),'b2npc7_fnay'//chns)
        if(style_b2mndt.eq.2) then                                     !srv 07.07.21 {
         wrk0 = vol*(nb-nap(:,:,isb))/dtim
        else
         wrk0 = 0.0_R8
         do iy = -1, ny
          do ix = -1, nx
           if((guard_flows.eq.2 .or.
     &        (is_neutral(isb) .and. use_eirene.ne.0)) .and.
     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
           ttim=dtim*dtco(isb,region(ix,iy,0))*time_factor(ix,iy)
           t0=vol(ix,iy)/ttim
           wrk0(ix,iy) = t0*(nb(ix,iy)-nap(ix,iy,isb))
          enddo
         enddo
        endif
        call my_out(70,nx,ny,wrk0,'b2npc7_dnadt'//chns)                 !srv 24.05.21 }
       endif
      endif

*   ..set error parameter
*     (at present all errors cause an abort through xerrab)
      ierr = 0
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npc7

      end subroutine b2npc7                                             !srv 02.04.18 }
