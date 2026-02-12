*-----------------------------------------------------------------------
*.specification

      subroutine b2uxus_dv(nCv, mpg, aa, aad, itcnt, resfun,  
     .  resfund, corfun, corfund, name, style, nbdirs)
      use b2mod_types
      use b2us_map_diffv
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, itcnt, nbdirs, style
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      real (kind=R8) ::
     *   aad(nbdirsmax, 13*nCv), resfund(nbdirsmax, nCv)
      character*(*) name 
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(nCv)
      real (kind=R8) ::
     *   corfund(nbdirsmax, nCv)
*   ..common blocks

*-----------------------------------------------------------------------
*     Rearranged routine which solves a linear system of equation A*x=b
*     also in the differentiated form Ad*x + A*xd = bd, using b2uxus
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nCv - integer, input.
*     See the calling routine for description.
*
*     aa - () real*8 array, input.
*     aa contains the stencil.
*
*     itcnt - integer, input.
*     itcnt specifies the iteration counter in the calling routine.
*
*     resfun - (nCv) real*8 array, input.
*     resfun specifies the right hand side.
*
*     corfun - (nCv) real*8 array, output.
*     corfun contains the computed solution.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) :: yy(nCv)
      integer :: iCv, nd
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxus_d')

*     First solve the original system of equation
      call b2uxus(nCv, mpg, aa, itcnt, resfun, corfun, name, style)

      do nd=1,nbdirs
*     Now solve the differentiated form Ad*x + A*xd = bd for xd
*     which becomes a system in the form A*xd = bd - Ad*x = yy
        yy = 0.0_R8
        do iCv = 1, nCv
          yy(iCv) = resfund(nd,iCv) - 
     .     sum(aad(nd,mpg%cvNvP(iCv,1):mpg%cvNvP(iCv,1)
     .     +mpg%cvNvP(iCv,2)-1)*corfun(mpg%cvNv(mpg%cvNvP(iCv,1)
     .     :mpg%cvNvP(iCv,1)+mpg%cvNvP(iCv,2)-1)))
        end do

        call b2uxus(nCv, mpg, aa, itcnt, yy, corfund(nd,:), name, style)
      end do

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uxus_dd

      end subroutine b2uxus_dv
