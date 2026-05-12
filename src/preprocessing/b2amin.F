*-----------------------------------------------------------------------
*.specification

      subroutine b2amin (nrid, nms_max, tlo, thi, nlo, nhi,
     &  nms, tailep)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nms_max
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   tlo, thi, nlo, nhi, tailep
      integer nms

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AMIN reads from a formatted file the list of (hydrogenic) 
*     molecular species
*     for which the rate coefficients are to be
*     computed.
*
*
*  3. description 
*
*     b2amin reads from unit (nrid) the variables: nms and tailep.
*     nms is the number of molecular species. tailep is a smoothing parameter.
*     Unit nrid must be connected to a formatted file. The input file
*     must be positioned properly. See the code to discern the precise
*     format required of the input file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nms_max - integer, input.
*     nms_max specifies the maximum allowed number of moleuclar species,
*     i.e. the maximum output value of nms.
*
*     tlo, thi - real, output.
*     tlo and thi specify the range of temperature over which the rate
*     coefficients are to be computed. tlo and thi are expressed in ev.
*
*     nlo, nhi - real, output.
*     nlo and nhi specify the range of density over which the rate
*     coefficients are to be computed. nlo and nhi are expressed in
*     1/m**3.
*
*     nms - integer, output.
*     nms specifies the number of molecular species.
*     The relation (1.le.nms.le.nms_max) will hold.
*
*     tailep - real.
*     tailep specifies a smoothing parameter for the temperature
*     dependence of some rate coefficients. It must be non-negative and
*     small compared to 1; tailep.eq.0.0 implies no smoothing.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      character id*8
      integer ims !! index of the molecular species
*   ..procedures
      logical streql
      external xertst, xerrab, streql

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2amin')
*   ..read heading for nms
      read (nrid,'(a8)',err=93) id
      call xertst (streql(id,'*nms'),
     &  'faulty input id instead of *nms, id='//id)
*   ..read nms
      read (nrid,*,err=93) nms
      call xertst (0.lt.nms.and.nms.le.nms_max,
     &  'faulty input nms')
*   ..read heading for tlo, thi
      read (nrid,'(a8)',err=91) id
      call xertst (streql(id,'*tlohi'),
     &  'faulty input id instead of *tlohi, id='//id)
*   ..read tlo, thi
      read (nrid,*,err=91) tlo, thi
      call xertst (1.0e-3_R8.le.tlo.and.tlo.lt.thi.and.thi.le.1.0e6_R8,
     &  'unreasonable input tlo, thi')
*   ..read heading for nlo, nhi
      read (nrid,'(a8)',err=92) id
      call xertst (streql(id,'*nlohi'),
     &  'faulty input id instead of *nlohi, id='//id)
*   ..read nlo, nhi
      read (nrid,*,err=92) nlo, nhi
      call xertst (1.0e12_R8.le.nlo.and.nlo.lt.nhi.and.nhi.le.1.0e30_R8,
     &  'unreasonable input nlo, nhi')
*   ..read heading for tailep
      read (nrid,'(a8)',err=96) id
      call xertst (streql(id,'*tailep'),
     &  'faulty input id instead of *tailep, id='//id)
*   ..read tailep
      read (nrid,*,err=96) tailep
      call xertst (0.le.tailep, 'faulty input tailep')
      call xertst (tailep.lt.1, 'unreasonable input tailep')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read tlo, thi')
   92 call xerrab ('error trying to read nlo, nhi')
   93 call xerrab ('error trying to read nnuc')
   96 call xerrab ('error trying to read tailep')

*-----------------------------------------------------------------------
*.end b2amin

      end subroutine b2amin
