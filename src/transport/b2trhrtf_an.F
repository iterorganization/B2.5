*-----------------------------------------------------------------------
*.specification

      subroutine b2trhrtf_an(nspecies, mum, mdm, nz2n,
     &  c_hwA_an, c_rwA_an, c_hTA_an, c_rTA_an, c_R_TA_an)
      use b2mod_types
      use b2mod_subsys
      use b2mod_b2zhco
     , , only: zh_analyt_mdf
      use b2mod_Zzh
      implicit none
*   ..input arguments (unchanged on exit)
      integer nspecies
      real (kind=R8) ::
     &  nz2n(0:nspecies-1),
     &  mum(0:nspecies-1,0:nspecies-1), mdm(0:nspecies-1,0:nspecies-1)
*   ..output arguments
      real (kind=R8) ::
     &  c_R_TA_an(0:nspecies-1,0:nspecies-1),
     &  c_hwA_an(0:nspecies-1,0:nspecies-1),
     &  c_rwA_an(0:nspecies-1,0:nspecies-1),
     &  c_hTA_an(0:nspecies-1), c_rTA_an(0:nspecies-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRHRTF_AN computes the heat flux, the R-moment and the thermal
*     force transport coefficient analytically.
*
*  3. description
*     - c_R_TA_an: the thermal force transport coefficient (not the same as c_R_TA)
*     - c_hwA_an: the heat flux velocity-dependent part of the transport coefficient
*     - c_rwA_an: the r-moment velocity-dependent part of the transport coefficient
*     - c_hTA_an: the heat flux temperature-dependent of the part transport coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i, j
      real (kind=R8) :: sqrt2
      real (kind=R8) :: Z5s(0:nspecies-1),
     &  Z9s(0:nspecies-1), Z11s(0:nspecies-1),
     &  rz2n(0:nspecies-1), Delta_al(0:nspecies-1)

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trhrtf_an')
*   ..main
      sqrt2=sqrt(2.0_R8)
*    ..Z-variables
      do i = 0, nspecies-1
        rz2n(0:nspecies-1)=nz2n(0:nspecies-1)/nz2n(i)
        Z5s(i)=b2mod_Zzh_Z5s(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
        Z9s(i)=b2mod_Zzh_Z9s(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
        Z11s(i)=b2mod_Zzh_Z11s(nspecies,i,mum(i,0:nspecies-1),
     &   mdm(i,0:nspecies-1),rz2n)
      enddo

*    ..denominator of the analytical transport coefficients
      do i = 0, nspecies-1
        Delta_al(i)=1+
     &   433.0_R8*sqrt2/288.0_R8*Z11s(i)+
     &    65.0_R8*sqrt2/32.0_R8*Z5s(i)-
     &    23.0_R8*sqrt2/16.0_R8*Z9s(i)+
     &  5629.0_R8/1152.0_R8*Z5s(i)*Z11s(i)-
     &   529.0_R8/128.0_R8*Z9s(i)**2
      enddo

*    ..the thermal force transport coefficient
      do i = 0, nspecies-1
        c_R_TA_an(i,0:nspecies-1)=25.0_R8*sqrt2/16.0_R8*
     &    mum(i,0:nspecies-1)*sqrt(mum(i,0:nspecies-1))*
     &   (1.5_R8+433.0_R8*sqrt2/240.0_R8*Z11s(i)-
     &    0.5_R8*mum(i,0:nspecies-1)*
     &   (1.0_R8+23.0_R8*sqrt2/8.0_R8*Z9s(i)))/Delta_al(i)
      enddo

*    ..the heat flux temperature-dependent part of the transport coefficient
      do i = 0, nspecies-1
        c_hTA_an(i)=125.0_R8/32.0_R8*
     &   (1.0_R8+433.0_R8*sqrt2/360.0_R8*Z11s(i))
     &   /Delta_al(i)
      enddo

*    ..the R-moment temperature-dependent part of the transport coefficient
      do i = 0, nspecies-1
        c_rTA_an(i)=175.0_R8/48.0_R8*
     &   (1.0_R8+23.0_R8*sqrt2/8.0_R8*Z9s(i))
     &   /Delta_al(i)
      enddo

*    ..the heat flux velocity-dependent part of the transport coefficient
      do i = 0, nspecies-1
        c_hwA_an(0:nspecies-1,i)=c_R_TA_an(i,0:nspecies-1)*
     &   (nz2n(0:nspecies-1)/nz2n(i))
      enddo

*    ..the R-moment velocity-dependent part of the transport coefficient
      do i = 0, nspecies-1
        c_rwA_an(0:nspecies-1,i)=-35.0_R8*sqrt2/24.0_R8*
     &    mum(i,0:nspecies-1)*sqrt(mum(i,0:nspecies-1))*
     &   (-1.5_R8-69.0_R8*sqrt2/16.0_R8*Z9s(i)+
     &     2.5_R8*mum(i,0:nspecies-1)*
     &   (1.0_R8+13.0_R8*sqrt2/8.0_R8*Z5s(i)))/Delta_al(i)*
     &   (nz2n(0:nspecies-1)/nz2n(i))
      enddo

*    ..the mdf correction
      if (zh_analyt_mdf.eq.1) then
        do j = 0, nspecies-1
          do i = 0, nspecies-1
            if (mdm(i,j).gt.1.0_R8) then
              c_hwA_an(j,i)=0.0_R8
              c_rwA_an(j,i)=0.0_R8
            endif
          enddo
        enddo
      endif

      call subend ()
      return

      end subroutine b2trhrtf_an
