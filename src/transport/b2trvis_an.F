*-----------------------------------------------------------------------
*.specification

      subroutine b2trvis_an (nspecies, mum, mdm, nz2n,
     &  c_pi_uA_an, c_pi_hA_an)
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
     &  c_pi_uA_an(0:nspecies-1),
     &  c_pi_hA_an(0:nspecies-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRVIS_AN computes the viscous stress tensor transport coefficient
*     analytically.
*
*  3. description
*     - c_pi_uA_an: the velocity-dependent viscosity transport coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i
      real (kind=R8) :: sqrt2
      real (kind=R8) ::
     &  Z11pi(0:nspecies-1), Z12pi(0:nspecies-1), Z22pi(0:nspecies-1),
     &  rz2n(0:nspecies-1),
     &  Deltapi_al(0:nspecies-1)

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trvis_an')
*   ..main
      sqrt2=sqrt(2.0_R8)
*    ..Z-variables
      do i = 0, nspecies-1
        rz2n(0:nspecies-1)=nz2n(0:nspecies-1)/nz2n(i)
        Z11pi(i)=b2mod_Zzh_Z11pi(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
        Z12pi(i)=b2mod_Zzh_Z12pi(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
        Z22pi(i)=b2mod_Zzh_Z22pi(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
      enddo

*    ..denominator of the analytical transport coefficients
      do i = 0, nspecies-1
        Deltapi_al(i)=1+
     &   205.0_R8*sqrt2/178.0_R8*Z11pi(i)+
     &   102.0_R8*sqrt2/89.0_R8*Z22pi(i)-
     &    54.0_R8*sqrt2/89.0_R8*Z12pi(i)+
     &   204.0_R8/89.0_R8*Z11pi(i)*Z22pi(i)-
     &   108.0_R8/89.0_R8*Z12pi(i)**2
      enddo

*    ..the velocity-dependent viscosity transport coefficient
      do i = 0, nspecies-1
        c_pi_uA_an(i)=1025.0_R8/1068.0_R8*
     &   (1.0_R8+204.0_R8*sqrt2/205.0_R8*Z22pi(i))/Deltapi_al(i)
      enddo

*    ..the velocity-dependent viscosity transport coefficient
      do i = 0, nspecies-1
        c_pi_hA_an(i)=1655.0_R8/1068.0_R8*
     &   (1.0_R8+204.0_R8*sqrt2/331.0_R8*Z22pi(i)+
     &           252.0_R8*sqrt2/331.0_R8*Z12pi(i))/Deltapi_al(i)
      enddo

      call subend ()
      return

      end subroutine b2trvis_an
