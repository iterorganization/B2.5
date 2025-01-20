*-----------------------------------------------------------------------
*.specification

      subroutine b2mxac (nCv, ns, dv, diag)
      use b2mod_types
      use b2us_plasma
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (B2Derivatives), intent (in) :: dv
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic), intent (inout) :: diag

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MXAC computes the (maximum-)norms of the corrections.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes maximum-norms of the corrections computed
*     towards the solution of the equations for continuity, parallel
*     momentum balance, electron heat balance, atom heat balance, and
*     charge balance.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is
*   ..procedures
      real (kind=R8) ::
     *   damax
      external xertst, damax
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxac')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')

* ..compute norms of the corrections
*   ..compute acorpa, acorua
      do is = 0, ns-1
       diag%acorpa(is) = damax(nCv,dv%corpa(1,is),1)
       diag%acorua(is) = damax(nCv,dv%corua(1,is),1)
      enddo
*   ..compute acorut, acorte, acorti, acorpo
      diag%acorut = damax(nCv,dv%corut,1)
      diag%acorte = damax(nCv,dv%corte,1)
      diag%acorti = damax(nCv,dv%corti,1)
      diag%acortn = damax(nCv,dv%cortn,1)
      diag%acorpo = damax(nCv,dv%corpo,1)
      diag%acorkt = damax(nCv,dv%corkt,1)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxac

      end subroutine b2mxac
