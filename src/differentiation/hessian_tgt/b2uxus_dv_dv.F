*-----------------------------------------------------------------------
*.specification

      subroutine b2uxus_dv_dv(nCv, mpg, aa, aad0, aad, aadd, itcnt,
     .  resfun, resfund0, resfund, resfundd, corfun, corfund0, corfund,
     .  corfundd, name, style, nbdirs, nbdirs0)
      use b2mod_types
      use b2us_map_diffv_diffv
      use b2mod_subsys
      use b2mod_diffsizes
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, itcnt, nbdirs, nbdirs0, style
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      real (kind=R8) ::
     *   aad(nbdirsmax, 13*nCv), resfund(nbdirsmax, nCv),
     *   aad0(nbdirsmax0, 13*nCv), resfund0(nbdirsmax0, nCv),
     *   aadd(nbdirsmax0, nbdirsmax, 13*nCv),
     *   resfundd(nbdirsmax0, nbdirsmax, nCv)
      character*(*) name 
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corfun(nCv)
      real (kind=R8) ::
     *   corfund(nbdirsmax, nCv), corfund0(nbdirsmax0, nCv),
     *   corfundd(nbdirsmax0, nbdirsmax, nCv)
*   ..common blocks

*-----------------------------------------------------------------------
*     Rearranged routine which solves a linear system of equation A*x=b
*     differentiated two times in tangent AD mode, using b2uxus
*     First tangent derivative routine solves d(A*x)=d(b) giving:
*     Ad*x + A*xd = bd
*     A*xd = bd - Ad*x = yy
*     The above is done via SOLVE_D(A,Ad,b,bd,x,xd)
*     So the command to differentiate again (with d=d0) are
*     A*x=b
*     yy = bd - Ad*x
*     A*xd = yy
*     This results in:
*     d0(A*x)=d0(b) --> SOLVE_D(A,Ad0,b,bd0,x,xd0)
*     yyd = bdd - Add*x - Ad*xd0
*     And finally:
*     d0(A*xd) = d0(yy) --> SOLVE_D(A,Ad0,yy,yyd,xd,xdd)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      real (kind=R8) :: yy(nCv), yyd(nbdirsmax0, nCv)
      integer :: iCv, nd, nd0
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxus_dv_dv')
     
*     First solve derivative of A*x=b
      call b2uxus_dv (nCv, mpg, aa, aad0, itcnt, resfun, resfund0,
     &  corfun, corfund0, name, style, nbdirs0)

      do nd = 1, nbdirs
*     Now evaluate yy = bd - Ad*x and yyd = bdd - Add*x - Ad*xd0
        yyd = 0.0_R8
        yy = 0.0_R8
        do iCv = 1, nCv
          yy(iCv) = resfund(nd,iCv) - 
     .     sum(aad(nd,mpg%cvNvP(iCv,1):mpg%cvNvP(iCv,1)
     .     +mpg%cvNvP(iCv,2)-1)*corfun(mpg%cvNv(mpg%cvNvP(iCv,1)
     .     :mpg%cvNvP(iCv,1)+mpg%cvNvP(iCv,2)-1)))
          do nd0 = nd, nbdirs0
            yyd(nd0,iCv) = resfundd(nd0,nd,iCv) - 
     .       sum(aadd(nd0,nd,mpg%cvNvP(iCv,1):mpg%cvNvP(iCv,1)
     .       +mpg%cvNvP(iCv,2)-1)*corfun(mpg%cvNv(mpg%cvNvP(iCv,1)
     .       :mpg%cvNvP(iCv,1)+mpg%cvNvP(iCv,2)-1))) -
     .       sum(aad(nd,mpg%cvNvP(iCv,1):mpg%cvNvP(iCv,1)
     .       +mpg%cvNvP(iCv,2)-1)*corfund0(nd0,mpg%cvNv(mpg%cvNvP(iCv,1)
     .       :mpg%cvNvP(iCv,1)+mpg%cvNvP(iCv,2)-1)))
          end do
        end do
*       And finally solve for derivative of A*xd=yy
        call b2uxus_dv (nCv, mpg, aa, aad0, itcnt, yy, yyd,
     &    corfund(nd, :), corfundd(:, nd, :), name, style, nbdirs0)
      end do

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uxus_dd_dd

      end subroutine b2uxus_dv_dv
