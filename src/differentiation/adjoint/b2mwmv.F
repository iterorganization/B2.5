*-----------------------------------------------------------------------
*.specification

      subroutine b2mwmv (nput, nCv, state, itim, tim, dtim)
      use b2mod_types
      use b2mod_subsys
      use b2us_plasma_diff
      use b2us_work
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nput, nCv, itim
      type(B2State), intent(in) :: state
      real (kind=R8), intent(in) :: tim, dtim
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWMV writes to an un*formatted file a selection of data
*     describing the state of the plasma.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) data that are suitable for
*     producing movie output. Unit nput must be connected to an
*     un*formatted file. Presently, only the changes in the plasma
*     state are written out.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nCv - integer, input.
*     nCv specifies the number of control volumes in the plasma grid.
*
*     state - plasma state data structure, input.
*     See the calling routine for description.
*
*     itim - integer, input.
*     itim specifies the timestep counter. It is initialised to 0 and
*     incremented by 1 following each successful call to b2mndt.
*
*     tim - real*8 input.
*     tim specifies the simulation time for the present call to b2mndt.
*     It will hold that 0.le.tim.
*
*     dtim - real, input.
*     dtim specifies the timestep for the present call to b2mndt.
*     It will hold that 0.lt.dtim.
*
*     na0, ua0, te0, ti0, po0 - real*8 array, input.
*     Previous state of the plasma.
*     See the calling routine for description.
*
*     na, ua, te, ti, po - real*8 array, input.
*     Present state of the plasma.
*     See the calling routine for description.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, idum(0:9)
      save ncall
      real (kind=R8) ::
     *   rdum(0:9)
*   ..procedures
      external xertst, b2scopy_nodiff, b2saxpy_nodiff, cfwuin, cfwure
      external b2xvsg
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwmv')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nCv
      call xertst (1.le.nCv, 'faulty argument nCv')
*   ..test dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvsg (nCv*state%ns, state%pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv, state%pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, state%pl%ti, 1, 'ti', '.gt.')
*    ..test old state
       call b2xvsg (nCv*state%ns, state%psnl%na, 1, 'na0', '.gt.')
       call b2xvsg (nCv, state%psnl%te, 1, 'te0', '.gt.')
       call b2xvsg (nCv, state%psnl%ti, 1, 'ti0', '.gt.')
      endif
*   ..compute change in plasma state
      call b2scopy_nodiff (nCv*state%ns, state%pl%na, 1, delna, 1)
      call b2saxpy_nodiff (nCv*state%ns, -1.0_R8, state%psnl%na, 1,
     &  delna, 1)
      call b2scopy_nodiff (nCv*state%ns, state%pl%ua, 1, delua, 1)
      call b2saxpy_nodiff (nCv*state%ns, -1.0_R8, state%psnl%ua, 1,
     &  delua, 1)
      call b2scopy_nodiff (nCv, state%pl%te, 1, delte, 1)
      call b2saxpy_nodiff (nCv, -1.0_R8, state%psnl%te, 1, delte, 1)
      call b2scopy_nodiff (nCv, state%pl%ti, 1, delti, 1)
      call b2saxpy_nodiff (nCv, -1.0_R8, state%psnl%ti, 1, delti, 1)
      call b2scopy_nodiff (nCv, state%pl%po, 1, delpo, 1)
      call b2saxpy_nodiff (nCv, -1.0_R8, state%psnl%po, 1, delpo, 1)

* ..save data for movie
*   ..save itim, tim, dtim
      idum(0) = itim
      call cfwuin (nput, 1, idum, 'itim')
      rdum(0) = tim
      rdum(1) = dtim
      call cfwure (nput, 2, rdum, 'tim,dtim')
*   ..save change in plasma state
      call cfwure (nput, nCv*state%ns, delna, 'delna')
      call cfwure (nput, nCv*state%ns, delua, 'delua')
      call cfwure (nput, nCv, delte, 'delte')
      call cfwure (nput, nCv, delti, 'delti')
      call cfwure (nput, nCv, delpo, 'delpo')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwmv

      end subroutine b2mwmv
