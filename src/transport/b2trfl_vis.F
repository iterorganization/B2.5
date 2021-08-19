*.specification

      subroutine b2trfl_vis (ix,iy,ns,eta,pa,gamma,
     & cF_limv) 
      use b2mod_types
      use b2mod_b2cmpt
      use b2mod_indirect
      use b2mod_zhfrtf , only: gavm_uc
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ix,iy,ns
      real (kind=R8) :: gamma                             !som 20.03.21
      real (kind=R8) :: eta(0:ns-1), pa(0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: 
     & cF_limv(0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2trfl_vis computes the flux limiting coefficient for the stress-visosity  
*
*  3. description 
*
*     This routine computes:
*     - cF_limv: the flux limiting coefficient for the thermal force coefficient
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer is
      real (kind=R8) :: t0,t1
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfl_vis')
*   ..main
      if (cflim(2).ne.0.0_R8.and..not.on_closed_surface(ix,iy)) then
        do is=0, ns-1
          t0=abs(gavm_uc(ix,iy))
          t1=cflim(2)*pa(is)
        
          cF_limv(is)=
     & 1.0_R8/(1.0_R8+(eta(is)*t0/t1)**gamma)**(1.0_R8/gamma)
        enddo
      else
        cF_limv(0:ns-1)=1.0_R8
      endif 

      call subend ()
      return

      end subroutine b2trfl_vis
