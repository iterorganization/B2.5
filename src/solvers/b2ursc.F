*-----------------------------------------------------------------------
*.specification

      subroutine b2ursc (nCv, nFc, geo, mpg, fun, sfn, ffn, resfn)
      use b2mod_types
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2ursc
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   fun(nCv), sfn(nCv,0:1), ffn(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resfn(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URSC computes the residual of a scalar conservation equation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of a scalar conservation
*     equation, integrated over cells of the b2 mesh. The equation has
*     the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resfn = rhs-div(flux) .
*     The flux is specified by ffn, and the right hand side has the
*     linearised representation
*       rhs() = sfn(,0)+sfn(,1)*fun() .
*     This linearisation may include a contribution from an implicit
*     discretisation in time.
*
*     The input argument ffn specifies fluxes integrated over cell
*     faces, and sfn specifies source coefficients integrated over
*     cell volumes. Therefore, the expression for the integrated
*     residual on any cell does not involve any geometric factors.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
*   ..procedures
      external subini, subend, xertst
      external b2xvsg, b2xvff
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ursc')
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2ursc.lt.3) then
*    ..test sign of sfn
       call b2xvsg (nCv, sfn(1,1), 1, 'sfn', '.le.')
      endif

* ..compute residual

*   ..contribution from right hand side
      resfn = sfn(:,0)+sfn(:,1)*fun

*   ..contribution from ffn
      do iFc = 1, nFc
        resfn(mpg%fcCv(iFc,1)) = 
     &    resfn(mpg%fcCv(iFc,1)) - ffn(iFc,0) - ffn(iFc,1)
        resfn(mpg%fcCv(iFc,2)) = 
     &    resfn(mpg%fcCv(iFc,2)) + ffn(iFc,0) + ffn(iFc,1)
      enddo

* ..return
      ncall_b2ursc = ncall_b2ursc+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ursc

      end subroutine b2ursc
