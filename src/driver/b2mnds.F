*-----------------------------------------------------------------------
*.specification

      subroutine b2mnds (ninp, nout, nCv, nFc, ns, ns0, switch)
      use b2mod_types
      use b2mod_math
      use b2mod_version
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
      use b2mod_ipmain
      use b2mod_ppout
      use b2mod_diag
     . , only : label
      use b2mod_switches
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:6), nout(0:10), nCv, nFc, ns, ns0
      type (switches), intent (inout) :: switch
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MNDS performs initialisation functions for a B2 calculation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine performs initialisation for a B2 calculation. It
*     reads and checks header records from the input files, writes
*     header records to the output files, and initialises variables in
*     common.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:5) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides physics and numerics parameters.
*     ninp(1): un*formatted; provides the geometry.
*     ninp(2): un*formatted; provides default physics parameters.
*     ninp(3): un*formatted; provides the initial state.
*     ninp(4): un*formatted; provides atomic rate data.
*     ninp(5): un*formatted; provides source term data.
*     ninp(6): un*formatted; provides the initial differentiated state.
*
*     nout - (0:5) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary print output.
*     nout(1): un*formatted; provides actual physics parameters.
*     nout(2): un*formatted; provides the final state.
*     nout(3): un*formatted; may provide data for movie.
*     nout(4): un*formatted; provides data on the numerical evolution.
*     nout(5): un*formatted; may provide detailed convergence data.
*     nout(9): un*formatted; provides the final differentiated state.
*     nout(10): un*formatted; provides data on the numerical evolution of differentiated quantities.
*
*     nCv, nFc - integer, input.
*     nCv and nFc specify the number of cells and faces in the
*     unstructured data format.
*     It will hold that 0.le.nCv and 0.le.nFc.
*
*     ns - integer, input.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     ns0 - integer, input.
*     ns0 specifies the number of atomic species in the data file
*     for the initial state.  The initial data for species (0:ns0-1)
*     are taken from that file, while the initial densities for species
*     (ns0:ns-1) are set to some very small value.
*     It will hold that 1.le.ns0.le.ns.
*
*     (parameters in common:)
*
*     /b2cmgs/ specifies the singular regions in the geometry.
*     /b2cmpa/ specifies basic physics parameters.
*     /b2cmpb/ specifies physics parameters for boundary conditions.
*     /b2cmpt/ specifies physics parameters for transport.
*     /b2cmrc/ contains a table of atomic rate coefficients.
*     /b2cmsr/ specifies source terms from neutral gas model.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, idum(0:9)
      character lblgm*120, lblcp*120, lblps*120, lblrc*120, lblsr*120,
     &  lblmn*120, cnamip*80, cvalip*80, line*80, ligne*80!, lblpsd*120
      real (kind=R8) ::
     *   zpdum(0:nsdecl-1,4)
*   ..procedures
      logical is_comment
      external xertst, xerrab, cfruch, cfruin, cfwuch, cfwuin
      external b2rucp, b2rurc, b2rusr, b2rflb, b2rfcp, b2xxid, b2mwq0,
     &  b2wfcp, b2wucp, b2ruzd, b2wuzd, b2xvcp
      external is_comment

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lblgm - character*120.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     lblcp - character*120.
*     lblcp identifies the instance of the program employed to prepare
*     the table of physics parameters.
*
*     lblps - character*120.
*     lblps identifies the instance of the program employed to prepare
*     the initial state.
*
*     lblrc - character*120.
*     lblrc identifies the instance of the program employed to prepare
*     the table of atomic rate coefficients.
*
*     lblsr - character*120.
*     lblsr identifies the instance of the program employed to prepare
*     the neutral gas source data.
*
*     lblmn - character*120.
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mnds')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.ninp(2).and.
     &             1.le.ninp(3).and.1.le.ninp(4).and.1.le.ninp(5).and.
     &             1.le.ninp(6).and.
     &             1.le.nout(0).and.1.le.nout(1).and.1.le.nout(2).and.
     &             1.le.nout(3).and.1.le.nout(4).and.1.le.nout(5).and.
     &             1.le.nout(9).and.1.le.nout(10),
     &    'faulty argument ninp, nout')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (ns.le.nsdecl, 'faulty parameter nsdecl')

* ..initialisation
*   ..obtain lblmn
      call b2xxid ('b2mn', lblmn(1:60))
      call b2rflb (ninp(0), lblmn(61:120))
*   ..read header for geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
#if ( defined(JET) || defined(IMAS) )
      labgeo=lblgm
      labgeo=lblgm
#endif
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
*     (isymm is now stored in the b2cmfs common)               ! xpb 10.08.00
*   ..read default physics common
      call cfruch (ninp(2), 120, lblcp, 'label')
      call b2rucp (ninp(2), b2fparam_version, ns)
*   ..read revised parameters in physics common
*     (ninp(0) need supply only parameters that are to be changed from
*     the values read in as default physics common.)
      call b2rfcp (ninp(0), ns)
*   ..test parameters
      call b2xvcp (ns)
*   ..read header for initial state
*     (zamin, zamax, zn, am from ninp(3) will be checked against the data
*     already read in)
      call cfruch (ninp(3), 120, lblps, 'label')
      call b2ruzd (ninp(3), b2fstati_version, ns0, 
     &  zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)
      do is = 0, ns0-1
       call xertst (ltst(zpdum(is,1),zamin(is)), 
     &  'faulty data for zamin')
       call xertst (ltst(zpdum(is,2),zamax(is)), 
     &  'faulty data for zamax')
       call xertst (ltst(zpdum(is,3),zn(is)), 'faulty data for zn')
       call xertst (ltst(zpdum(is,4),am(is)), 'faulty data for am')
      enddo
*   ..read atomic rate data
      call cfverr (ninp(4), b2frates_version)
      call cfruch (ninp(4), 120, lblrc, 'label')
      adpak_used = index(lblrc,'ADPAK').gt.0 
      call read_b2mod_ppout('rates', ninp(4), lblrc)
      call b2rurc (ninp(4), b2frates_version)
*   ..test rtns, rtzmin, rtzmax, rtzn
      call xertst (ns.eq.rtns, 'faulty input rtns--wrong rate table?')
      do is = 0, ns-1
       call xertst (ltst(zamin(is),rtzmin(is)),
     &   'faulty input rtzmin--wrong rate table?')
       call xertst (ltst(zamax(is),rtzmax(is)),
     &   'faulty input rtzmax--wrong rate table?')
       call xertst (ltst(zn(is),rtzn(is)),
     &   'faulty input rtzn--wrong rate table?')
      enddo
*   ..read source data
*!!   (at present we have no external source data)
*$$$  call cfruch (ninp(5), 120, lblsr, 'label')
      lblsr(1:60) = 'external source terms ...'
      lblsr(61:120) = '... not defined'
      call b2rusr (ninp(5))
*   ..test source data
*     (no data, no test)
*   ..write all labels
      write (nout(0),'(/2x,a)') 'labels for this calculation:'
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblgm:',
     &  lblgm(1:60), lblgm(61:120)
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblcp:',
     &  lblcp(1:60), lblcp(61:120)
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblps:',
     &  lblps(1:60), lblps(61:120)
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblrc:',
     &  lblrc(1:60), lblrc(61:120)
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblsr:',
     &  lblsr(1:60), lblsr(61:120)
      write (nout(0),'(2x,a/4x,a/4x,a)') 'lblmn:',
     &  lblmn(1:60), lblmn(61:120)
      label = lblmn
*   ..read and echo code internal parameters
      write (nout(0),'(/2x,a)') 'non-default internal parameters:'
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         ligne = line
         call strip_spaces(ligne)
         if (ligne(1:1).eq.'''') then
           read (line,*) cnamip, cvalip
           write (nout(0),'(4x,a,2x,a)') cnamip, trim(cvalip)
           call ipsetc (cnamip, cvalip)
         endif
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
      write (nout(0),'()')
*   ..check for prohibited switches
      call ipprhp ()
*   ..test /b2cmpa/, /b2cmpb/, /b2cmpt/
      call b2xvcp (ns)
*   ..write physics parameters
      call b2wfcp (nout(0), 0, ns)
*   ..write header to nout(1)
      idum(0) = ns
      call cfverw (nout(1), newversion)                                 !xpb
      call cfwuin (nout(1), 1, idum, 'ns')
      call cfwuch (nout(1), 120, lblmn, 'label')
*   ..save physics parameters to nout(1)
      call b2wucp (nout(1), newversion, ns)
*   ..write header to nout(2), nout(3) and nout(5)
      idum(0) = nCv
      idum(1) = nFc
      idum(2) = ns
      call cfverw (nout(2), newversion)                                 !xpb
      call cfwuin (nout(2), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(2), 120, lblmn, 'label')
      call b2wuzd (nout(2), newversion, ns, zamin, zamax, zn, am)
      call cfverw (nout(3), newversion)                                 !xpb
      call cfwuin (nout(3), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(3), 120, lblmn, 'label')
      call b2wuzd (nout(3), newversion, ns, zamin, zamax, zn, am)
      call cfverw (nout(5), newversion)                                 !xpb
      call cfwuin (nout(5), 3, idum, 'nCv,nFc,ns')
      call cfwuch (nout(5), 120, lblmn, 'label')
      call b2wuzd (nout(5), newversion, ns, zamin, zamax, zn, am)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c The switches are now also loaded in a module b2mod_switches.
c Hence, they can be accessed in the classical way using ipgeti/ipgetr
c or in a new way using this module.
c The b2mod_switches is mainly intended for switches that are used in
c different functions/subroutines to reduce code duplication
c and to centralize consistency checks on the switches.

c read main model switches to correctly set the corresponding defaults
      call read_first_switches(switch)
c load default values of switches
      call set_defaults_switches(switch,ns)
c overwrite with values provided in b2mn.dat
      call read_switches(switch)
c check the values of the switches
      call check_switches(switch,ns)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

*   ..write header to nout(4)
      idum(0) = ns
      call cfwuin (nout(4), 1, idum, 'ns')
      call cfwuch (nout(4), 120, lblmn, 'label')
      call b2mwq0 (nout(4), ns, switch)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2mnds

      end subroutine b2mnds
