*-----------------------------------------------------------------------
*.specification

      subroutine b2ursd (nCv, nFc, mpg, den, fun, sdf, fdf, resdf)
      use b2mod_types
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2ursd
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   den(nCv), fun(nCv), sdf(nCv,0:3), fdf(nFc,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resdf(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URSD computes the residual of a scalar conservation equation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of a scalar conservation
*     equation, integrated over cells of the B2 mesh. The equation has
*     the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resdf = rhs-div(flux) .
*     The flux is specified by fdf, and the right hand side has the
*     linearised representation
*       rhs() = sdf(,0)+
*         sdf(,1)*fun()+
*         sdf(,2)*den()+
*         sdf(,3)*den()*fun() .
*     This linearisation may include a contribution from an implicit
*     discretisation in time.
*
*     The input argument fdf specifies fluxes integrated over cell
*     faces, and sdf specifies source coefficients integrated over
*     cell volumes. Therefore, the expression for the integrated
*     residual on any cell does not involve any geometric factors.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
*   ..procedures
      external xertst
      external b2xvsg, b2xvff
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ursd')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2ursd.lt.3) then
*    ..test sign of den, sdf
       call b2xvsg (nCv, den, 1, 'den', '.ge.')
       call b2xvsg (nCv, sdf(1,1), 1, 'sdf1', '.le.')
       call b2xvsg (nCv, sdf(1,3), 1, 'sdf3', '.le.')
      endif

* ..compute residual
*   ..contribution from right hand side
      resdf = sdf(:,0) + sdf(:,1)*fun + sdf(:,2)*den + sdf(:,3)*den*fun

*   ..contribution from fdf
      do iFc = 1, nFc
        resdf(mpg%fcCv(iFc,1)) = 
     &    resdf(mpg%fcCv(iFc,1)) - fdf(iFc,0) - fdf(iFc,1)
        resdf(mpg%fcCv(iFc,2)) = 
     &    resdf(mpg%fcCv(iFc,2)) + fdf(iFc,0) + fdf(iFc,1)
      enddo

* ..return
      ncall_b2ursd = ncall_b2ursd+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ursd

      end subroutine b2ursd
