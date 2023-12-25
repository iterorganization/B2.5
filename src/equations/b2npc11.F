*-----------------------------------------------------------------------
*.specification

      subroutine b2npc11 (nx, ny, fac_vis, guard_flows, zb,              !srv 05.07.08 02.04.18 {
     &  nregionv, solving, solvcreg, ua_solve, solvmreg, itcnt, isb,
     &  rxf, dtim, vol, pbs, flob, conb, snb,                            !srv 22.01.02 13.10.06 07.07.21
     &  pccb, po, ti, pb, nb, ub, fnb, rescb, corpb, ierr)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_balance
     & ,only : update_balance_sources_particles, balance_netcdf,     !djm Jan2017
     &         tot_sna0to1, fna_tot                                  !djm Mar2019
      use b2mod_subsys
      use b2mod_lwmain
      use b2mod_sources
     , , only : snadt                                                !srv 31.01.21
      use b2mod_numerics_namelist                                    !srv 07.07.21 {
     & , only : dtco, time_factor
       use b2mod_plasma                                              !srv 07.07.21 }
     & , only : nap
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, guard_flows, nregionv, itcnt, isb              !srv 05.07.08
      real (kind=R8) ::
     *  rxf, dtim, vol(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),           !srv 07.07.21
     &  fac_vis(-1:nx,-1:ny), flob(-1:nx,-1:ny,0:1,0:1),
     &  conb(-1:nx,-1:ny,0:1,0:1,0:4), snb(-1:nx,-1:ny,0:1),         !srv 13.10.06
     &  pccb(-1:nx,-1:ny,0:1,0:1), fnb(-1:nx,-1:ny,0:1,0:1),
     &  po(-1:nx,-1:ny), ti(-1:nx,-1:ny), zb(-1:nx,-1:ny)
      logical solving, solvcreg(0:nregionv), ua_solve,               !srv 22.05.18
     &  solvmreg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *   pb(-1:nx,-1:ny), nb(-1:nx,-1:ny), ub(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *   rescb(-1:nx,-1:ny), corpb(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPC11 performs one iteration of a relaxation process for the
*     solution of the continuity equation for one atomic species using
*     an 11-point stencil.
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
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0          !srv 06.07.06
      integer, save :: iout = 0, iout_b2wdat = 0                     !srv 04.07.08
      integer, save :: get_residuals = 0                             !srv 22.05.18
      integer, save :: style_b2mndt = 1                              !srv 07.07.21
      integer, save :: use_9pt_stencil = 1
      character chns*3                                               !srv 04.07.08
      logical warning
      !$OMP THREADPRIVATE(iout, iout_b2wdat)
      !$OMP THREADPRIVATE(no_Ptncr_x_co, no_Ptncr_y_co)
      real (kind=R8) ::
     *   rxg, pcm0, pcm1, ttim, t0                                   !srv 07.07.21
      real (kind=R8) ::
     *   aa11(-1:nx,-1:ny,-7:7), wrk0(-1:nx,-1:ny)
      save ncall, rxg, pcm0, pcm1, use_eirene
      !$OMP THREADPRIVATE(ncall, rxg, pcm0, pcm1, use_eirene)
*   ..procedures
      external xertst, ipgetr
      external b2xvsg, b2xvff, b2ursc, b2usc9, b2upco
*   ..initialization
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
      call subini ('b2npc11')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2mndt_style', style_b2mndt)                       !srv 07.07.21 {
       call ipgeti ('b2mndr_use_9pt_stencil', use_9pt_stencil)
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
        call ipgetr ('b2npc11_rxg', rxg)
        call ipgetr ('b2npc11_pcm0', pcm0)
        call ipgetr ('b2npc11_pcm1', pcm1)
       endif                                                           !srv 07.07.21 }
       call xertst (0.lt.rxg, 'faulty internal parameter rxg')
       call xertst (0.le.pcm0, 'faulty internal parameter pcm0')
       call xertst (0.le.pcm1, 'faulty internal parameter pcm1')
       call ipgeti ('b2npc11_iout', iout)                             !srv 13.07.05
       call ipgeti ('b2wdat_iout', iout_b2wdat)                       !iys 17.04.17
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)               !srv 06.07.06
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)               !srv 06.07.06
       call ipgeti ('get_residuals', get_residuals)                   !srv 22.05.18
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf
      call xertst (0.le.rxf.and.rxf.le.1, 'faulty argument rxf')
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

      if(style_b2mndt.eq.2) then                                         !srv 07.07.21 {
*  .. calculate time-dependent source
        call b2scdt (nx, ny, isb, guard_flows, dtim,
     &   vol, nap(-1,-1,isb), snadt(-1,-1,0,isb), .false.)
*  .. add contribution from time-dependent source
        snb(:,:,:) = snb(:,:,:) + snadt(:,:,:,isb)
      endif                                                              !srv 24.05.21 }

cdjm Jan2017 update particle source components for balance
      if (balance_netcdf.ne.0) then
        tot_sna0to1(:,:,:,isb)=snb !djm Mar2019
        fna_tot(:,:,:,:,isb)=fnb   !djm Mar2019
        call update_balance_sources_particles(nx,ny,isb,nb)
      endif

      if (iout.ne.0) then                                             !srv 12.07.02 {
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,pb,'b2npc11_pa0_'//chns)
        call my_out(70,nx,ny,nb,'b2npc11_na0_'//chns)
        wrk0=snb(:,:,0)+nb*snb(:,:,1)
        call my_out(70,nx,ny,wrk0,'b2npc11_sna_input'//chns)
      endif

* ..particle balance
*   ..compute residual
      if((solving.and.any(solvcreg(0:nnreg(0))))
     &           .or.get_residuals.eq.1) then                          !srv 22.05.18
       call b2ursc (nx, ny, nb, snb, fnb, rescb)
       if(guard_flows.eq.2 .or.                                        !srv 06.07.08 {
     &    (is_neutral(isb) .and. use_eirene.ne.0)) then                !srv 15.02.10 { this part should be re-written
        do iy=-1, ny, ny+1                                             !srv 27.06.08 {
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 ) then
            rescb(ix,iy) = 0.0_R8
          endif
         enddo
        enddo                                                          !srv 27.06.08 }
       endif                                                           !srv 15.02.10 }
      endif                                                           !srv 22.05.18
*   ..compute correction
      write (chns,'(i3.3)') isb
      if(solving.and.any(solvcreg(0:nnreg(0)))) then                  !srv 22.05.18
        call b2usc11 (nx, ny, fac_vis, isb, zb, nregionv, solvcreg,   !srv 13.10.06 26.09.12
     &    itcnt, rxg, pcm1,
     &    vol, pbs, po, ti, nb, pb, snb(-1,-1,1), flob,
     &    conb, pccb(-1,-1,1,1), rescb, corpb, aa11,
     &    'b2usc11_'//chns)                                           !srv 18.05.02 02.04.18
*   ..apply correction
        call b2upco (nx, ny, nregionv, ua_solve, solvmreg, rxf, pcm0,
     &   vol, corpb, pccb(-1,-1,0,0), pb, nb, ub, isb)                !IYS 19.12.2017
        if(guard_flows.eq.2 .or.                                      !srv 06.07.08 {
     &    (is_neutral(isb) .and. use_eirene.ne.0)) then               !srv 06.07.08 {
*   ..replace the incorrect corner values by the interpolated values  !srv 30.04.08 {
          iy=-1                                                       !srv 27.06.08 {
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2) then
            if(use_9pt_stencil.ne.0) then
             nb(ix,iy) = nb(rightix(ix,iy),rightiy(ix,iy))
            else
             nb(ix,iy) = exp(
     .         log(nb(rightix(ix,iy),rightiy(ix,iy))) +
     .         log(nb(topix(ix,iy),topiy(ix,iy))) -
     &         log(nb(topix(rightix(ix,iy),rightiy(ix,iy)),
     &                topiy(rightix(ix,iy),rightiy(ix,iy)))))
            end if
           end if
           if(rightix(ix,iy).eq.nx+1) then
            if(use_9pt_stencil.ne.0) then
             nb(ix,iy) = nb(leftix(ix,iy),leftiy(ix,iy))
            else
             nb(ix,iy) = exp(
     .          log(nb(leftix(ix,iy),leftiy(ix,iy))) +
     .          log(nb(topix(ix,iy),topiy(ix,iy))) -
     &          log(nb(topix(leftix(ix,iy),leftiy(ix,iy)),
     &                 topiy(leftix(ix,iy),leftiy(ix,iy)))))
            end if
           end if
          enddo
          iy=ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2) then
             if(use_9pt_stencil.ne.0) then
              nb(ix,iy) = nb(rightix(ix,iy),rightiy(ix,iy))
             else
              nb(ix,iy) = exp(
     .          log(nb(rightix(ix,iy),rightiy(ix,iy))) +
     .          log(nb(bottomix(ix,iy),bottomiy(ix,iy))) -
     &          log(nb(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &                 bottomiy(rightix(ix,iy),rightiy(ix,iy)))))
             end if
            end if
            if(rightix(ix,iy).eq.nx+1) then
             if(use_9pt_stencil.ne.0) then
              nb(ix,iy) = nb(leftix(ix,iy),leftiy(ix,iy))
             else
              nb(ix,iy) = exp(
     .          log(nb(leftix(ix,iy),leftiy(ix,iy))) +
     .          log(nb(bottomix(ix,iy),bottomiy(ix,iy))) -
     &          log(nb(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &                 bottomiy(leftix(ix,iy),leftiy(ix,iy)))))
             end if
            end if
          enddo                                                       !srv 27.06.08 }
        endif                                                         !srv 06.07.08 }
      endif                                                           !srv 22.05.18
      if (iout.ne.0) then                                             !srv 12.07.02 {
       if (in_parallel()) then
         !$OMP MASTER
         write(*,*) 'B2NPC11 will not write to file in parallel mode'
         !$OMP END MASTER
       else
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,pb,'b2npc11_pa'//chns)
        call my_out(70,nx,ny,nb,'b2npc11_na'//chns)
        call my_out(70,nx,ny,ub,'b2npc11_ua'//chns)
        call my_out(70,nx,ny,rescb,'b2npc11_resco'//chns)
        call my_out(70,nx,ny,snb(-1,-1,0),'b2npc11_sna0'//chns)
        call my_out(70,nx,ny,snb(-1,-1,1),'b2npc11_sna1'//chns)
        wrk0=snb(:,:,0)+nb*snb(:,:,1)
        call my_out(70,nx,ny,wrk0,'b2npc11_sna'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,0),'b2npc11_fnax'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,1),'b2npc11_fnay'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,1),'b2npc11_fnaxy'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,0),'b2npc11_fnayx'//chns)
       endif
      elseif (iout_b2wdat.eq.4) then                                  !iys 17.04.17
       if (in_parallel()) then
         !$OMP MASTER
         write(*,*) 'B2NPC11 will not write to file in parallel mode'
         !$OMP END MASTER
       else
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,nb,'b2npc11_na'//chns)
        call my_out(70,nx,ny,rescb,'b2npc11_resco'//chns)
        wrk0=snb(:,:,0)+nb*snb(:,:,1)
        call my_out(70,nx,ny,wrk0,'b2npc11_sna'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,0),'b2npc11_fnax'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,1),'b2npc11_fnay'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,0,1),'b2npc11_fnaxy'//chns)
        call my_out(70,nx,ny,fnb(-1,-1,1,0),'b2npc11_fnayx'//chns)
        if(style_b2mndt.eq.2) then                                    !srv 24.05.21 {
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
        call my_out(70,nx,ny,wrk0,'b2npc11_dnadt'//chns)              !srv 24.05.21 }
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
*.end b2npc11

      end subroutine b2npc11                                          !srv 02.04.18 }
