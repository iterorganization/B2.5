*.specification

      subroutine b2trnu (ns,nnucl,zmaxmax)                 !som 05.02.2021  
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ns
*   ..output arguments (unspecified on entry)
      integer :: nnucl, zmaxmax

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRNU computes the number of different types of species (ions with the 
*     different nucleus) and the maximal atom number Z in the plasma composition
*
*  3. description 
*
*     This routine computes:
*     - nnucl: the number of different types of species 
*     - zmaxmax: the maximal atom number Z in the plasma composition
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer is
      real (kind=R8) :: zamax_temp, am_temp
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trnu')
*   ..main
        zamax_temp=zamax(0)
        am_temp=am(0)
        nnucl=1 
        do is = 0, ns-1
!       print *, nnucl
          if (zamax(is).ne.zamax_temp .and. am(is).ne.am_temp) then
           zamax_temp=zamax(is)
           am_temp=am(is)
           nnucl=nnucl+1 
          endif 
        zmaxmax=max(zmaxmax,nint(zamax_temp))
        enddo        
      call subend ()
      return


      end subroutine b2trnu
