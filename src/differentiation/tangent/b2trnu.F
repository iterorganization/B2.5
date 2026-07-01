*.specification

      subroutine b2trnu (ns, nnucl, znmax)                             !som 05.02.2021
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ns
*   ..output arguments (unspecified on entry)
      integer :: nnucl, znmax
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRNU computes the number of isonuclear sequences (ions with
*     different nuclei) and the maximal atom number Z in the plasma
*     composition.
*
*  3. description
*
*     This routine computes:
*     - nnucl: the number of isonuclear sequences
*     - znmax: the maximal atom number Z in the plasma composition
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer is
      real (kind=R8) :: zn_temp, am_temp
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trnu')
*   ..main
      zn_temp=zn(0)
      am_temp=am(0)
      nnucl=1
      znmax=0
      do is = 0, ns-1
        if (nint(zn(is)).ne.nint(zn_temp) .or.
     &      nint(am(is)).ne.nint(am_temp)) then
          zn_temp=zn(is)
          am_temp=am(is)
          nnucl=nnucl+1
        endif
        znmax=max(znmax,nint(zn_temp))
      enddo
      call subend ()
      return

      end subroutine b2trnu
