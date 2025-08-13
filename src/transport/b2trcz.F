*-----------------------------------------------------------------------
*.specification

      subroutine b2trcz(nspecies, mum, mdm, nz2n,
     &  c_R_TB, c_R_wB, c_hwB, c_rwB, c_hTB, c_rTB)
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
      real (kind=R8) :: c_R_TB(0:nspecies-1),
     &  c_R_wB(0:nspecies-1,0:nspecies-1), c_hwB(0:nspecies-1),
     &  c_rwB(0:nspecies-1), c_hTB(0:nspecies-1),
     &  c_rTB(0:nspecies-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCZ calculates charge state corrections to the thermal and
*     friction force terms.
*
*  3. description 
*
*     This routine computes:
*     - c_R_TB: correction for the thermal force
*     - c_R_wB: correction for the friction force
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i, j
      real (kind=R8) :: sqrt2
      real (kind=R8) :: Z2s(0:nspecies-1), Z5s(0:nspecies-1),
     &  Z8s(0:nspecies-1), Z9s(0:nspecies-1), Z11s(0:nspecies-1),
     &  rz2n(0:nspecies-1), Dpart_al(0:nspecies-1),
     &  mumi(0:nspecies-1), mdmi(0:nspecies-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trcz')
*   ..main
*    ..Z variables
      sqrt2=sqrt(2.0_R8) 
      do i = 0, nspecies-1
        mumi(0:nspecies-1)=mum(i,0:nspecies-1)
        mdmi(0:nspecies-1)=mdm(i,0:nspecies-1)
        rz2n(0:nspecies-1)=nz2n(0:nspecies-1)/nz2n(i)
        Z2s (i)=b2mod_Zzh_Z2s (nspecies,i,mumi,rz2n)
        Z5s (i)=b2mod_Zzh_Z5s (nspecies,i,mumi,mdmi,rz2n)
        Z8s (i)=b2mod_Zzh_Z8s (nspecies,i,mumi,rz2n)
        Z9s (i)=b2mod_Zzh_Z9s (nspecies,i,mumi,mdmi,rz2n)
        Z11s(i)=b2mod_Zzh_Z11s(nspecies,i,mumi,mdmi,rz2n)
      enddo
*    ..denominator of the corrections
      do i = 0, nspecies-1
        Dpart_al(i)=1+
     &   204376.0_R8*sqrt2/160413.0_R8*Z11s(i)+
     &    72670.0_R8*sqrt2/53471.0_R8*Z5s(i)-
     &    76728.0_R8*sqrt2/53471.0_R8*Z9s(i)+
     &   360256.0_R8/160413.0_R8*Z5s(i)*Z11s(i)-
     &   101568.0_R8/53471.0_R8*Z9s(i)**2
      enddo

      do i = 0, nspecies-1
*    ..the correction for the thermal force
        c_R_TB(i)=31500.0_R8/53471.0_R8*
     &   (559.0_R8/420.0_R8*(1+6928.0_R8*sqrt2/8385.0_R8*Z11s(i))*
     &   (1+2.0_R8*sqrt2*Z2s(i))-139.0_R8/420.0_R8*
     &   (1+184.0_R8*sqrt2/139.0_R8*Z9s(i))*(1+4.0_R8*sqrt2*Z8s(i)))/
     &    Dpart_al(i)

*    ..the correction for the heat flux (temperature)
        c_hTB(i)=139750.0_R8/53471.0_R8*
     &   (1+6928.0_R8*sqrt2/8385.0_R8*Z11s(i))/
     &    Dpart_al(i)

*    ..the correction for the R-moment (temperature)
        c_rTB(i)=194600.0_R8/53471.0_R8*
     &   (1+184.0_R8*sqrt2/139.0_R8*Z9s(i))/
     &    Dpart_al(i)

*    ..the correction for the heat flux (velocity)
        c_hwB(i)=31500.0_R8/53471.0_R8*
     &  (1-139.0_R8*sqrt2/105.0_R8*Z8s(i)-
     &      46.0_R8*sqrt2/105.0_R8*Z9s(i)+
     &    1732.0_R8*sqrt2/1575.0_R8*Z11s(i)+
     &     559.0_R8*sqrt2/210.0_R8*Z2s(i)-
     &     368.0_R8/105.0_R8*Z8s(i)*Z9s(i)+
     &    6928.0_R8/1575.0_R8*Z2s(i)*Z11s(i))/
     &   Dpart_al(i)

*    ..the correction for the R-moment (velocity)
        c_rwB(i)=17080.0_R8/53471.0_R8*
     &   (1+276.0_R8*sqrt2/61.0_R8*Z9s(i)+
     &      417.0_R8*sqrt2/61.0_R8*Z2s(i)-
     &      590.0_R8*sqrt2/61.0_R8*Z8s(i)-
     &      130.0_R8*sqrt2/61.0_R8*Z5s(i)+
     &     1104.0_R8/61.0_R8*Z2s(i)*Z9s(i)-
     &     1040.0_R8/61.0_R8*Z5s(i)*Z8s(i))/
     &   Dpart_al(i)

*    ..the correction for the thermal force
        do j = 0, nspecies-1
          c_R_wB(j,i)=1-0.6_R8*mum(i,j)*c_hwB(i)+
     &     3.0_R8/14.0_R8*mum(i,j)**2*c_rwB(i)
        enddo
      enddo

      call subend ()
      return

      end subroutine b2trcz
