*-----------------------------------------------------------------------
*.specification

      subroutine b2mxac_diff (nCv, ns, dv, diag)
      use b2mod_types
      use b2us_plasma_diff
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (B2Derivatives), intent (in) :: dv
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic), intent (inout) :: diag
*-----------------------------------------------------------------------
*.documentation
* This is a modified version of b2mxac used in the B2.5 code differentiated
* with TAPENADE. The main difference wrt b2mxac is the B2Derivatives_diff
* type used for input argument dv, which differs from the non-diff one.
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, is
      save ncall
*   ..procedures
      real (kind=R8) ::
     *   samax
      external xertst, samax
      external b2xvfx
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxac')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')


* ..compute norms of the corrections
*   ..compute acorpa, acorua
      do is = 0, ns-1
       diag%acorpa(is) = samax(nCv,dv%corpa(1,is),1)
       diag%acorua(is) = samax(nCv,dv%corua(1,is),1)
      enddo
*   ..compute acorut, acorte, acorti, acorpo
      diag%acorut = samax(nCv,dv%corut,1)
      diag%acorte = samax(nCv,dv%corte,1)
      diag%acorti = samax(nCv,dv%corti,1)
      diag%acortn = samax(nCv,dv%cortn,1)
      diag%acorpo = samax(nCv,dv%corpo,1)
      diag%acorkt = samax(nCv,dv%corkt,1)

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxac

      end
