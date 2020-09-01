*-----------------------------------------------------------------------
*.specification

      subroutine b2mwq0 (nput, ns, switch)
      use b2mod_types
      use b2mod_switches
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, ns
      type (switches), intent (in) :: switch


*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWQ0 writes a header for the output from b2mwqt.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the header information for
*     the 'quick trace' evolution data written by subroutine b2mwqt.
*     Unit (nput) must be connected to an un*formatted file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer idum(0:9), nfld, k, is, iix, iiy, ll(0:1)
      character chtit*60, chlab*12, id*3, idr*7
*   ..procedures
      external subini, subend, xertst, ipgeti, cfwuin, cfwuch
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     style - integer.
*     style specifies the amount of output produced on each call.
*     See the companion routine b2mwqt for a precise description.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwq0')
*   ..test nput
      call xertst (1.le.nput, 'faulty argument nput')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set nfld
      nfld = 0
      if (1.le.switch%b2mwqt_style) then
       nfld = nfld+2+2*(2*ns+4)
      endif
      if (2.le.switch%b2mwqt_style) then
       nfld = nfld+36*(ns+3)
      endif
      if (3.le.switch%b2mwqt_style) then
       nfld = nfld+40*(ns+3)
      endif
      if (4.le.switch%b2mwqt_style) then
       nfld = nfld+72*(ns+3)
      endif
      nfld = nfld + 1
*   ..set ll
*     (silly time; ll is needed as an argument to cfwuin)
      ll(0) = 0
      ll(1) = 1
*   ..write nfld
      idum(0) = nfld
      call cfwuin (nput, 1, idum(0), 'nfld')

* ..write titles and labels
      k = 0
      if (1.le.switch%b2mwqt_style) then
*    ..iteration counters
       chtit = 'time step counter'
       chlab = 'itim'
       call cfwuin (nput, 1, ll(0), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'iteration counter'
       chlab = 'itcnt'
       call cfwuin (nput, 1, ll(0), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
*     ..norms of residuals
       do is = 0, ns-1
        write (id,'(i3)') is
        chtit =
     &    'norm of residual of continuity equation, species '//id
        chlab = 'aresco('//id//')'
        call cfwuin (nput, 1, ll(1), 'linlog')
        call cfwuch (nput, 60, chtit, 'title')
        call cfwuch (nput, 12, chlab, 'label')
        k = k+1
       enddo
       do is = 0, ns-1
        write (id,'(i3)') is
        chtit =
     &    'norm of residual of parallel momentum equation, species '
     &    //id
        chlab = 'aresmo('//id//')'
        call cfwuin (nput, 1, ll(1), 'linlog')
        call cfwuch (nput, 60, chtit, 'title')
        call cfwuch (nput, 12, chlab, 'label')
        k = k+1
       enddo
       chtit = 'norm of residual of total parallel momentum equation'
       chlab = 'aresmt'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of residual of electron heat equation'
       chlab = 'areshe'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of residual of all atom heat equation'
       chlab = 'areshi'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of residual of the fluid neutral heat equation'
       chlab = 'areshn'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of residual of electric potential equation'
       chlab = 'arespo'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
*    ..norms of corrections
       do is = 0, ns-1
        write (id,'(i3)') is
        chtit =
     &    'norm of correction of pressure, species '//id
        chlab = 'acorpa('//id//')'
        call cfwuin (nput, 1, ll(1), 'linlog')
        call cfwuch (nput, 60, chtit, 'title')
        call cfwuch (nput, 12, chlab, 'label')
        k = k+1
       enddo
       do is = 0, ns-1
        write (id,'(i3)') is
        chtit =
     &    'norm of correction of velocity, species '//id
        chlab = 'acorua('//id//')'
        call cfwuin (nput, 1, ll(1), 'linlog')
        call cfwuch (nput, 60, chtit, 'title')
        call cfwuch (nput, 12, chlab, 'label')
        k = k+1
       enddo
       chtit = 'norm of correction of combined velocity'
       chlab = 'acorut'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of correction of electron temperature'
       chlab = 'acorte'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of correction of all atom temperature'
       chlab = 'acorti'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
       chtit = 'norm of correction of electric potential'
       chlab = 'acorpo'
       call cfwuin (nput, 1, ll(1), 'linlog')
       call cfwuch (nput, 60, chtit, 'title')
       call cfwuch (nput, 12, chlab, 'label')
       k = k+1
      endif
      if (2.le.switch%b2mwqt_style) then
*    ..total content on regions
       do is = 0, ns-1
        write (id,'(i3)') is
        do iiy = -1, 4
         do iix = -1, 4
          write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
          chtit = 'total particle content, species '//id//
     &      ', region '//idr
          chlab = 'int-na('//id//')'
          call cfwuin (nput, 1, ll(0), 'linlog')
          call cfwuch (nput, 60, chtit, 'title')
          call cfwuch (nput, 12, chlab, 'label')
          k = k+1
         enddo
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'total electron heat content, region '//idr
         chlab = 'int-ee'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'total all atom heat content, region '//idr
         chlab = 'int-ei'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'total electric potential (qe*ne*po), region '//idr
         chlab = 'int-po'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
      endif
      if (3.le.switch%b2mwqt_style) then
*    ..integrated fluxes
       do is = 0, ns-1
        write (id,'(i3)') is
        do iiy = 0, 3
         do iix = 0, 4
          write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
          chtit = 'particle flux along x, species '//id//
     &      ', region '//idr
          chlab = 'rfnax'
          call cfwuin (nput, 1, ll(0), 'linlog')
          call cfwuch (nput, 60, chtit, 'title')
          call cfwuch (nput, 12, chlab, 'label')
          k = k+1
         enddo
        enddo
        do iiy = 0, 4
         do iix = 0, 3
          write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
          chtit = 'particle flux along y, species '//id//
     &      ', region '//idr
          chlab = 'rfnay'
          call cfwuin (nput, 1, ll(0), 'linlog')
          call cfwuch (nput, 60, chtit, 'title')
          call cfwuch (nput, 12, chlab, 'label')
          k = k+1
         enddo
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'electron heat flux along x, region '//idr
         chlab = 'rfhex'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'electron heat flux along y, region '//idr
         chlab = 'rfhey'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'all atom heat flux along x, region '//idr
         chlab = 'rfhix'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'all atom heat flux along y, region '//idr
         chlab = 'rfhiy'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = 0, 3
        do iix = 0, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'electric charge flux along x, region '//idr
         chlab = 'rfchx'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = 0, 4
        do iix = 0, 3
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'electric charge flux along y, region '//idr
         chlab = 'rfchy'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
      endif
      if (4.le.switch%b2mwqt_style) then
*    ..min and max on regions
       do is = 0, ns-1
        write (id,'(i3)') is
        do iiy = -1, 4
         do iix = -1, 4
          write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
          chtit = 'minimum density, species '//id//', region '//idr
          chlab = 'na0('//id//')'
          call cfwuin (nput, 1, ll(0), 'linlog')
          call cfwuch (nput, 60, chtit, 'title')
          call cfwuch (nput, 12, chlab, 'label')
          k = k+1
          chtit = 'maximum density, species '//id//', region '//idr
          chlab = 'na1('//id//')'
          call cfwuin (nput, 1, ll(0), 'linlog')
          call cfwuch (nput, 60, chtit, 'title')
          call cfwuch (nput, 12, chlab, 'label')
          k = k+1
         enddo
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'minimum electron temperature, region '//idr
         chlab = 'te0'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
         chtit = 'maximum electron temperature, region '//idr
         chlab = 'te1'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'minimum all atom temperature, region '//idr
         chlab = 'ti0'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
         chtit = 'maximum all atom temperature, region '//idr
         chlab = 'ti1'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
       do iiy = -1, 4
        do iix = -1, 4
         write (idr,'(a,i2,a,i2,a)') '(', iix, ',', iiy, ')'
         chtit = 'minimum electric potential, region '//idr
         chlab = 'po0'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
         chtit = 'maximum electric potential, region '//idr
         chlab = 'po1'
         call cfwuin (nput, 1, ll(0), 'linlog')
         call cfwuch (nput, 60, chtit, 'title')
         call cfwuch (nput, 12, chlab, 'label')
         k = k+1
        enddo
       enddo
      endif
*   ..test k
      call xertst (k.eq.nfld, 'faulty parameter nfld')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwq0

      end
