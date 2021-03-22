
      module b2mod_Zzh
      use b2mod_types
      implicit none
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*     b2mod_Zzh contains functions to calculate Z-variables
*
*-----------------------------------------------------------------------
      contains

*.........................................................................
      function b2mod_Zzh_Z2s(nspecies,i,mumal,rz2n)
      implicit none
      integer :: i,j
      integer nspecies
      real (kind=R8) :: b2mod_Zzh_Z2s,
     & mumal(0:nspecies-1),rz2n(0:nspecies-1)

      b2mod_Zzh_Z2s=0.0_R8
      do j = 0, nspecies-1
        if (i.ne.j) then
          b2mod_Zzh_Z2s=b2mod_Zzh_Z2s+mumal(j)*sqrt(mumal(j))*rz2n(j)
        endif
      enddo
      end function b2mod_Zzh_Z2s

*.........................................................................
      function b2mod_Zzh_Z5s(nspecies,i,mumal,mdmal,rz2n)
      implicit none
      integer :: i,j
      integer nspecies
      real (kind=R8) :: b2mod_Zzh_Z5s,
     & mumal(0:nspecies-1),mdmal(0:nspecies-1),rz2n(0:nspecies-1)

      b2mod_Zzh_Z5s=0.0_R8
      do j = 0, nspecies-1
        if (i.ne.j) then
          b2mod_Zzh_Z5s=b2mod_Zzh_Z5s+
     & (1+16.0_R8/13.0_R8*mdmal(j)+30.0_R8/13.0_R8*mdmal(j)**2)*
     & mumal(j)**2*sqrt(mumal(j))*rz2n(j)
        endif
      enddo
      end function b2mod_Zzh_Z5s

*.........................................................................
      function b2mod_Zzh_Z8s(nspecies,i,mumal,rz2n)
      implicit none
      integer :: i,j
      integer nspecies
      real (kind=R8) :: b2mod_Zzh_Z8s,
     & mumal(0:nspecies-1),rz2n(0:nspecies-1)

      b2mod_Zzh_Z8s=0.0_R8
      do j = 0, nspecies-1
        if (i.ne.j) then
          b2mod_Zzh_Z8s=b2mod_Zzh_Z8s+
     & mumal(j)**2*sqrt(mumal(j))*rz2n(j)
        endif
      enddo
      end function b2mod_Zzh_Z8s

*.........................................................................
      function b2mod_Zzh_Z9s(nspecies,i,mumal,mdmal,rz2n)
      implicit none
      integer :: i,j
      integer nspecies
      real (kind=R8) :: b2mod_Zzh_Z9s,
     & mumal(0:nspecies-1),mdmal(0:nspecies-1),rz2n(0:nspecies-1)

      b2mod_Zzh_Z9s=0.0_R8
      do j = 0, nspecies-1
        if (i.ne.j) then
          b2mod_Zzh_Z9s=b2mod_Zzh_Z9s+
     & (1+32.0_R8/23.0_R8*mdmal(j)+84.0_R8/23.0_R8*mdmal(j)**2)*
     & mumal(j)**3*sqrt(mumal(j))*rz2n(j)
        endif
      enddo
      end function b2mod_Zzh_Z9s

*.........................................................................
      function b2mod_Zzh_Z11s(nspecies,i,mumal,mdmal,rz2n)
      implicit none
      integer :: i,j
      integer nspecies
      real (kind=R8) :: b2mod_Zzh_Z11s,
     & mumal(0:nspecies-1),mdmal(0:nspecies-1),rz2n(0:nspecies-1)

      b2mod_Zzh_Z11s=0.0_R8
      do j = 0, nspecies-1
        if (i.ne.j) then
          b2mod_Zzh_Z11s=b2mod_Zzh_Z11s+
     & (1+1088.0_R8/433.0_R8*mdmal(j)+3672.0_R8/433.0_R8*mdmal(j)**2+
     & 1792.0_R8/433.0_R8*mdmal(j)**3+1400.0_R8/433.0_R8*mdmal(j)**4)*
     & mumal(j)**4*sqrt(mumal(j))*rz2n(j)
        endif
      enddo
      end function b2mod_Zzh_Z11s
      end module b2mod_Zzh
