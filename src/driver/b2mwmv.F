*-----------------------------------------------------------------------
*.specification

      subroutine b2mwmv (nput, nx, ny, ns, itim, tim, dtim)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_indirect
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny, ns, itim
      real (kind=R8) ::
     *   tim, dtim
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks

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
*     nx, ny, ns - integer, input.
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
*     (parameters in common:)
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, idum(0:9)
      save ncall
      real (kind=R8) ::
     *   rdum(0:9)
*   ..procedures
      external xertst, b2scopy, b2saxpy, cfwuin, cfwure
      external b2xvsg, b2xvfx
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwmv')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test dtim
      call xertst (0.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test old state
       call b2xvsg (n2*ns, na0, 1, 'na0', '.gt.')
       call b2xvsg (n2, te0, 1, 'te0', '.gt.')
       call b2xvsg (n2, ti0, 1, 'ti0', '.gt.')
      endif
*   ..compute change in plasma state
      call b2scopy (n2*ns, na, 1, delna, 1)
      call b2saxpy (n2*ns, -1.0_R8, na0, 1, delna, 1)
      call b2scopy (n2*ns, ua, 1, delua, 1)
      call b2saxpy (n2*ns, -1.0_R8, ua0, 1, delua, 1)
      call b2scopy (n2, te, 1, delte, 1)
      call b2saxpy (n2, -1.0_R8, te0, 1, delte, 1)
      call b2scopy (n2, ti, 1, delti, 1)
      call b2saxpy (n2, -1.0_R8, ti0, 1, delti, 1)
      call b2scopy (n2, po, 1, delpo, 1)
      call b2saxpy (n2, -1.0_R8, po0, 1, delpo, 1)

* ..save data for movie
*   ..save itim, tim, dtim
      idum(0) = itim
      call cfwuin (nput, 1, idum, 'itim')
      rdum(0) = tim
      rdum(1) = dtim
      call cfwure (nput, 2, rdum, 'tim,dtim')
*   ..save change in plasma state
      call cfwure (nput, n2*ns, delna, 'delna')
      call cfwure (nput, n2*ns, delua, 'delua')
      call cfwure (nput, n2, delte, 'delte')
      call cfwure (nput, n2, delti, 'delti')
      call cfwure (nput, n2, delpo, 'delpo')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwmv

      end subroutine b2mwmv
