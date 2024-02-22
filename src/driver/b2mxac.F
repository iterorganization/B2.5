*-----------------------------------------------------------------------
*.specification

      subroutine b2mxac (nx, ny, ns, corpa, corua, corut, corte, corti,
     &  corpo, acorpa, acorua, acorut, acorte, acorti, acorpo)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  corpa(-1:nx,-1:ny,0:ns-1), corua(-1:nx,-1:ny,0:ns-1),
     &  corut(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny),
     &  corpo(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  acorpa(0:ns-1), acorua(0:ns-1), acorut, acorte, acorti,
     &  acorpo

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
      integer n2, is
*   ..procedures
      real (kind=R8) ::
     *   damax
      external xertst, damax
      external b2xvfx
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxac')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')


* ..compute norms of the corrections
*   ..compute acorpa, acorua
      do is = 0, ns-1
       acorpa(is) = damax(n2,corpa(-1,-1,is),1)
       acorua(is) = damax(n2,corua(-1,-1,is),1)
      enddo
*   ..compute acorut, acorte, acorti, acorpo
      acorut = damax(n2,corut,1)
      acorte = damax(n2,corte,1)
      acorti = damax(n2,corti,1)
      acorpo = damax(n2,corpo,1)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxac

      end subroutine b2mxac
