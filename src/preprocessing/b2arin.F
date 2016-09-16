*-----------------------------------------------------------------------
*.specification

      subroutine b2arin (nrid, nnucm, nzmax, tlo, thi, nlo, nhi,
     &  nnuc, nz, izlo, izhi, flag, tailep)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nnucm, nzmax
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   tlo, thi, nlo, nhi, tailep
      integer nnuc, nz(0:nnucm-1), izlo(0:nnucm-1), izhi(0:nnucm-1)
      character flag*8

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2ARIN reads from a formatted file the list of nuclear species
*     and atomic states for which the rate coefficients are to be
*     computed.
*
*
*  3. description (see also routine b2cdca)
*
*     b2arin reads from unit (nrid) the variables: nnuc, nz(0:nnuc-1),
*     izlo(0:nnuc-1), izhi(0:nnuc-1), flag and tailep.
*     nnuc is the number of nuclear species, nz() are the nuclear
*     charges, izlo() and izhi() specify the range of atomic states,
*     and flag specifies which procedure is to be used to calculate the
*     rate coefficients. tailep is a smoothing parameter.
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
*     nnucm - integer, input.
*     nnucm specifies the maximum allowed number of nuclear species,
*     i.e. the maximum output value of nnuc.
*
*     nzmax - integer, input.
*     nzmax specifies the maximum allowed value of the nuclear charge
*     of any species.
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
*     nnuc - integer, output.
*     nnuc specifies the number of nuclear species.
*     The relation (1.le.nnuc.le.nnucm) will hold.
*
*     nz - (0:nnucm-1) integer array, output.
*     For inuc in (0:nnuc-1), nz(inuc) specifies the nuclear charge of
*     the nuclear species inuc. It will hold that 1.le.nz(inuc).
*     For inuc in (nnuc:nnucm-1), nz(inuc) will hold 0.
*
*     izlo, izhi - (0:nnucm-1) integer array, output.
*     For inuc in (0:nnuc-1), (izlo(inuc):izhi(inuc)) specifies the
*     range of atomic charge state represented in b2 for the nuclear
*     species inuc. It will hold that
*       0.le.izlo(inuc).le.izhi(inuc)+1 and izhi(inuc).le.nz(inuc).
*     For inuc in (nnuc:nnucm-1), (izlo(inuc),izhi(inuc)) will hold
*     (0,-1).
*
*     flag - character*8, output.
*     flag specifies which code is to be used to calculate the rate
*     coefficients. Its value may be 'adpak' or 'strahl'.
*
*     tailep - real.
*     tailep specifies a smoothing parameter for the temperature
*     dependence of some rate coefficients. It must be non-negative and
*     small compared to 1; tailep.eq.0.0 implies no smoothing.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      character id*8, flag0*8
      integer inuc
*   ..procedures
      logical streql
      external subini, subend, xertst, xerrab, streql

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2arin')
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
*   ..read heading for nnuc
      read (nrid,'(a8)',err=93) id
      call xertst (streql(id,'*numnuc'),
     &  'faulty input id instead of *numnuc, id='//id)
*   ..read nnuc
      read (nrid,*,err=93) nnuc
      call xertst (0.lt.nnuc.and.nnuc.le.nnucm,
     &  'faulty input nnuc')
*   ..read heading for nz
      read (nrid,'(a8)',err=94) id
      call xertst (streql(id,'*nucspec'),
     &  'faulty input id instead of *nucspec, id='//id)
*   ..read nz, izlo, izhi
      do inuc = 0, nnuc-1
       read (nrid,*,err=94) nz(inuc), izlo(inuc), izhi(inuc)
       call xertst (0.lt.nz(inuc).and.nz(inuc).le.nzmax,
     &   'faulty input nz')
       call xertst (0.le.izlo(inuc).and.
     &   izlo(inuc).le.izhi(inuc)+1.and.izhi(inuc).le.nz(inuc),
     &   'faulty input izlo, izhi')
      enddo
*   ..set remaining nz, izlo, izhi
      do inuc = nnuc, nnucm-1
       nz(inuc) = 0
       izlo(inuc) = 0
       izhi(inuc) = -1
      enddo
*   ..read heading for flag
      read (nrid,'(a8)',err=95) id
      call xertst (streql(id,'*flag'),
     &  'faulty input id instead of *flag, id='//id)
*   ..read flag
      read (nrid,*,err=95) flag0
      if (streql(flag0,'adpak')) then
       flag = 'adpak'
      else if (streql(flag0,'strahl')) then
       flag = 'strahl'
      else if (streql(flag0,'adas')) then
       flag = 'adas'
#ifdef ITM
      else if (streql(flag0,'amns')) then
       flag = 'amns'
#endif
      else
       call xerrab ('faulty input flag, flag='//flag0)
      endif
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
   94 call xerrab ('error trying to read nz, izlo, izhi')
   95 call xerrab ('error trying to read flag')
   96 call xerrab ('error trying to read tailep')

*-----------------------------------------------------------------------
*.end b2arin

      end
