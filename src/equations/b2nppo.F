*-----------------------------------------------------------------------
*.specification

      subroutine b2nppo (nx, ny, nregionv, solving, solvereg, itcnt,     !srv 22.05.18
     &  rxf, vol, ne, te, conc, sch, po, fch, respo, corpo, ierr)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv, itcnt
      real (kind=R8) ::
     *  rxf, ne(-1:nx,-1:ny), te(-1:nx,-1:ny),
     &  conc(-1:nx,-1:ny,0:1,0:1), sch(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny)
      logical solving, solvereg(0:nregionv)                              !srv 22.05.18
*   ..input/output arguments
      real (kind=R8) ::
     *  po(-1:nx,-1:ny), fch(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *  respo(-1:nx,-1:ny), corpo(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPPO performs one iteration of a relaxation process for the
*     solution of the charge conservation equation, which governs the
*     electric potential.
*
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equation that
*     governs the electric potential po.
*     An outline of the procedure follows.
*     ..compute residual
*     ..compute correction
*     ..apply correction
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2
      integer, save :: iout = 0, iout_b2wdat = 0
      integer, save :: get_residuals = 0                              !srv 22.05.18
      real (kind=R8) ::
     *  aa9(-1:nx,-1:ny,-6:6), wrk0(-1:nx,-1:ny)
      save ncall
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, b2ursd, b2uspo, b2uppo
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nppo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2nppo_iout', iout)                            !srv 08.07.08
       call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
       call ipgeti ('get_residuals', get_residuals)                 !srv 22.05.18
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
*    ..test sign of ne, te
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
*    ..test sign of conc
       call b2xvsg (n2, conc(-1,-1,0,0), 1, 'conc00', '.ge.')
       call b2xvsg (n2, conc(-1,-1,1,1), 1, 'conc11', '.ge.')
*    ..test edge values of conc
       call b2xvff (nx, ny, conc, 'conc')
*    ..test sign of sch
       call b2xvsg (n2, sch(-1,-1,1), 1, 'sch1', '.le.')
       call b2xvsg (n2, sch(-1,-1,3), 1, 'sch3', '.le.')
*    ..test edge values of fch
       call b2xvff (nx, ny, fch, 'fch')
      endif

* ..main computation
*   ..compute residual
      if((solving.and.any(solvereg(0:nregionv))).or.
     &                get_residuals.eq.1) then                      !srv 22.05.18
       call b2ursd (nx, ny, ne, po, sch, fch, respo)
      endif                                                         !srv 22.05.18
c-rfs      write(6,*) 'ix=20, iy=-1, po ',po(20,-1), po(21,-1)
c-rfs      write(6,*) 'ix=20, iy=-1, sch ',
c-rfs     .  sch(20,-1,0), sch(20,-1,1),
c-rfs     .  sch(20,-1,2), sch(20,-1,3)
c-rfs      write(6,*) 'ix=20, iy=-1, fch, respo',
c-rfs     .  fch(20,-1,0), fch(20,-1,1),
c-rfs     .  fch(20,0,1), respo(20,-1)
c-rfs      write(6,*) 'ix=21, iy=-1, fch, respo',
c-rfs     .  fch(21,-1,0), fch(21,-1,1),
c-rfs     .  fch(21,0,1), respo(21,-1)
*   ..compute correction
      if(solving.and.any(solvereg(0:nregionv))) then                !srv 22.05.18
        call b2uspo (nx, ny, nregionv, solvereg, itcnt, vol,
     &    ne, te, conc, sch, respo, corpo, aa9, 'b2nppo')           !srv 18.03.10
*   ..apply correction
        call b2uppo (nx, ny, rxf, corpo, po)
      endif                                                         !srv 22.05.18

      if(iout .eq. 1 .or. iout_b2wdat.eq.4) then                    !srv 11.09.09 {
!       call my_out(70,nx,ny,te,'b2nppo_te')
!       call my_out(70,nx,ny,ne,'b2nppo_ne')
        call my_out(70,nx,ny,po,'b2nppo_po')
        call my_out(70,nx,ny,respo,'b2nppo_respo')
        call my_out(70,nx,ny,fch(-1,-1,0,0),'b2nppo_fchx')
        call my_out(70,nx,ny,fch(-1,-1,1,1),'b2nppo_fchy')
!       call my_out(70,nx,ny,sch(-1,-1,0),'b2nppo_sch0')
!       call my_out(70,nx,ny,sch(-1,-1,1),'b2nppo_sch1')
!       call my_out(70,nx,ny,sch(-1,-1,2),'b2nppo_sch2')
!       call my_out(70,nx,ny,sch(-1,-1,3),'b2nppo_sch3')
         wrk0=sch(:,:,0)+sch(:,:,1)*po+sch(:,:,2)*ne+
     &        sch(:,:,3)*ne*po
        call my_out(70,nx,ny,wrk0,'b2nppo_sch')
      endif                                                       !srv 11.09.09 }

*   ..set error parameter and return
*     (at present all errors cause an abort through xertst)
      ierr = 0
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nppo

      end subroutine b2nppo
