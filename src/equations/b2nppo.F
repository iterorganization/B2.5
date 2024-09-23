*-----------------------------------------------------------------------
*.specification

      subroutine b2nppo (nCv, nFc, nVx, nregionv, solving, solvereg,     !srv 22.05.18
     &  itcnt, rxf, switch, geo, mpg, pl, dv, sr, ierr)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_ad
     & , only : ncall_b2nppo
#ifndef NAGFOR
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2uspo,
     &          my_out_folder,
#ifdef WG_DEBUG
     &          my_out_sf,                                               !srv 27.06.21
#endif
     &          ncall_b2ursd
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, itcnt
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  rxf
      logical solving, solvereg(0:nregionv)                              !srv 22.05.18
*   ..input/output arguments
      type (B2Plasma), intent (inout) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Source), intent (inout) :: sr
*   ..output arguments (unspecified on entry)
      integer ierr
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
      real (kind=R8) ::
     *  aa(mpg%nCmxNv), wrk0(nCv), wrk
      real (kind=R8), save :: b2nppo_restr_po = 0.0_R8
      integer :: iCv
*   ..procedures
      external xertst
      external b2xvsg, b2ursd, b2uspo, b2uppo
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nppo')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     . 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall_b2nppo.lt.3) then
*    ..test sign of ne, te
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
*    ..test sign of sch
       call b2xvsg (nCv, sr%sch(1,1), 1, 'sch1', '.le.')
       call b2xvsg (nCv, sr%sch(1,3), 1, 'sch3', '.le.')
      endif
      if (ncall_b2nppo.eq.0) then
            call ipgetr ('b2nppo_restr_po', b2nppo_restr_po)
      endif

* ..main computation
*   ..compute residual
      if((solving.and.any(solvereg(0:nregionv))).or.
     &    switch%get_residuals.eq.1) then                           !srv 22.05.18
        call b2ursd (nCv, nFc, mpg, dv%ne, pl%po, sr%sch, dv%fch,
     &    dv%respo)
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
        call b2uspo (nCv, nFc, nVx, nregionv, solvereg, itcnt, switch,
     &   geo, mpg, dv%ne, pl%te, dv%conc, sr%sch, dv%respo, dv%corpo,aa,
     &  'b2nppo')
*   ..apply correction
        call b2uppo (nCv, rxf, dv%corpo, pl%po)
        if (abs(b2nppo_restr_po).gt.1.0_R8) then
         do iCv=1,mpg%nCv
            if (mpg%cvOnClosedSurface(iCv)) cycle
            if (b2nppo_restr_po.gt.0.0_R8) then
              if ((abs(qe*pl%po(iCv)/(pl%te(iCv)))
     &                   .gt. b2nppo_restr_po )) then
                    wrk = pl%po(iCv)
                    pl%po(iCv) = sign(b2nppo_restr_po
     &                   *pl%te(iCv)/qe,pl%po(iCv))
                   print *, "Change po",iCv, wrk,pl%po(iCv)
              endif
            else
             if ((abs(qe*pl%po(iCv)/(pl%te(iCv)))
     &                   .gt. abs(b2nppo_restr_po)).and.
     &             (pl%po(iCv).lt.0.0_R8)) then
                    wrk = pl%po(iCv)
                    pl%po(iCv) = sign(b2nppo_restr_po
     &                   *pl%te(iCv)/qe,pl%po(iCv))
                   print *, "Change po",iCv, wrk,pl%po(iCv)
             endif
            endif
         enddo
       endif
      endif                                                         !srv 22.05.18

      if(switch%b2nppo_iout .eq. 1 .or. switch%iout_b2wdat.eq.4) then    !srv 11.09.09 {
!       call my_out_us(70,nCv,0,pl%te,'b2nppo_te')
!       call my_out_us(70,nCv,0,dv%ne,'b2nppo_ne')
        call my_out_us(70,nCv,0,pl%po,'b2nppo_po')
        call my_out_us(70,nCv,0,dv%respo,'b2nppo_respo')
        call my_out_us(70,nFc,1,dv%fch(1,0),'b2nppo_fch_th')
        call my_out_us(70,nFc,1,dv%fch(1,1),'b2nppo_fch_r')
!       call my_out_us(70,nCv,0,sr%sch(1,0),'b2nppo_sch0')
!       call my_out_us(70,nCv,0,sr%sch(1,1),'b2nppo_sch1')
!       call my_out_us(70,nCv,0,sr%sch(1,2),'b2nppo_sch2')
!       call my_out_us(70,nCv,0,sr%sch(1,3),'b2nppo_sch3')
         wrk0=sr%sch(:,0)+sr%sch(:,1)*pl%po+sr%sch(:,2)*dv%ne+
     &        sr%sch(:,3)*dv%ne*pl%po
        call my_out_us(70,nCv,0,wrk0,'b2nppo_sch')
      endif                                                       !srv 11.09.09 }

*   ..set error parameter and return
*     (at present all errors cause an abort through xertst)
      ierr = 0
* ..return
      ncall_b2nppo = ncall_b2nppo+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nppo

      end subroutine b2nppo
