*-----------------------------------------------------------------------
*.specification

      subroutine b2trcs(nspecies, mum, mdm, nz2n, c_pi_uB, c_pi_hB)
      use b2mod_types
      use b2mod_subsys
      use b2mod_Zzh
      implicit none
*   ..input arguments (unchanged on exit)
      integer nspecies
      real (kind=R8) ::
     &  nz2n(0:nspecies-1),
     &  mum(0:nspecies-1,0:nspecies-1), mdm(0:nspecies-1,0:nspecies-1)
*   ..output arguments
      real (kind=R8) ::
     &  c_pi_uB(0:nspecies-1), c_pi_hB(0:nspecies-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCS calculates charge state corrections for the viscous
*     stress tensor.
*
*  3. description
*
*     This routine computes:
*     - c_pi_uB: correction for the velocity-dependent part
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i
      real (kind=R8) :: sqrt2
      real (kind=R8) ::
     &  Z11pi(0:nspecies-1), Z12pi(0:nspecies-1), Z22pi(0:nspecies-1),
     &  rz2n(0:nspecies-1), Dpipart_al(0:nspecies-1),
     &  mumi(0:nspecies-1), mdmi(0:nspecies-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trcs')
*   ..main
*    ..Z variables
      sqrt2=sqrt(2.0_R8)
      do i = 0, nspecies-1
        mdmi(0:nspecies-1)=mdm(i,0:nspecies-1)
        mumi(0:nspecies-1)=mum(i,0:nspecies-1)
        rz2n(0:nspecies-1)=nz2n(0:nspecies-1)/nz2n(i)
        Z11pi(i)=b2mod_Zzh_Z11pi(nspecies,i,mumi,mdmi,rz2n)
        Z12pi(i)=b2mod_Zzh_Z12pi(nspecies,i,mumi,mdmi,rz2n)
        Z22pi(i)=b2mod_Zzh_Z22pi(nspecies,i,mumi,mdmi,rz2n)
      enddo
*    ..denominator of the corrections
      do i = 0, nspecies-1
        Dpipart_al(i)=1+
     &   159.0_R8*sqrt2/167.0_R8*Z11pi(i)+
     &   816.0_R8*sqrt2/835.0_R8*Z22pi(i)-
     &   108.0_R8*sqrt2/167.0_R8*Z12pi(i)+
     &  1224.0_R8/835.0_R8*Z11pi(i)*Z22pi(i)-
     &   648.0_R8/835.0_R8*Z12pi(i)**2
      enddo

*    ..the correction for the velocity-dependent part
      do i = 0, nspecies-1
        c_pi_uB(i)=265.0_R8/334.0_R8*
     &  (1+204.0_R8*sqrt2/265.0_R8*Z22pi(i))/
     &   Dpipart_al(i)
      enddo

*    ..the correction for the heat flux dependent part
      do i = 0, nspecies-1
        c_pi_hB(i)=475.0_R8/334.0_R8*
     &  (1+204.0_R8*sqrt2/475.0_R8*Z22pi(i)+
     &     252.0_R8*sqrt2/475.0_R8*Z12pi(i))/
     &   Dpipart_al(i)
      enddo
      call subend ()
      return

      end subroutine b2trcs
