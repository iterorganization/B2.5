*-----------------------------------------------------------------------
*.specification

      subroutine b2trcn (nnucl, amnucl, N2, N3, N4, N5, N6, N8, N9,  
     & N10, N11, N12, N13, N14, N15, N16, mu, rmu, rmd2, mum, mdm)      !som 02.08.2021
      use b2mod_types
!     use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
      integer nnucl
      real (kind=R8) :: amnucl(0:nnucl-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::N2(0:nnucl-1,0:nnucl-1), N3(0:nnucl-1,0:nnucl-1),
     & N4(0:nnucl-1,0:nnucl-1), N5(0:nnucl-1,0:nnucl-1),
     & N6(0:nnucl-1,0:nnucl-1), N8(0:nnucl-1,0:nnucl-1),
     & N9(0:nnucl-1,0:nnucl-1),
     & N10(0:nnucl-1,0:nnucl-1), N11(0:nnucl-1,0:nnucl-1),
     & N12(0:nnucl-1,0:nnucl-1),
     & N13(0:nnucl-1,0:nnucl-1), N14(0:nnucl-1,0:nnucl-1),
     & N15(0:nnucl-1,0:nnucl-1), N16(0:nnucl-1,0:nnucl-1),
     & mu(0:nnucl-1,0:nnucl-1), mum(0:nnucl-1,0:nnucl-1),
     & mdm(0:nnucl-1,0:nnucl-1), rmu(0:nnucl-1,0:nnucl-1),
     & rmd2(0:nnucl-1)



*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCN computes the N-coefficients and all mass related variables
*
*  3. description 
*
*     This routine computes:
*     - mu: the reduced mass (in atom mass numbers)
*     - rmu: the square root of the reduced mass (in atom mass numbers)
*     - rmd2: the square root of the mass devided by 2 (in atom mass numbers)
*     - kappa: the energy transfer factor
*     - mdm: the mass to mass ratio
*     - mum: the reduced mass to mass ratio
*     - N*: the N-coefficients, which are used in the G-objects
*
*-----------------------------------------------------------------------
*.declarations
      integer i,j
      real (kind=R8) :: kappa(0:nnucl-1,0:nnucl-1)
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trcn')
*   ..main
      do i = 0, nnucl-1
        do j = 0, nnucl-1
          mu(i,j)=amnucl(i)*amnucl(j)/(amnucl(i)+amnucl(j))
          rmu(i,j)=sqrt(mu(i,j))
          kappa(i,j)=amnucl(i)*amnucl(j)/(amnucl(i)+amnucl(j))**2
          mdm(i,j)=amnucl(i)/amnucl(j)
        enddo
      enddo
      
      do i = 0, nnucl-1
        rmd2(i)=rmu(i,i)
        do j = 0, nnucl-1
          mum(i,j)=mu(i,j)/amnucl(i)
          N2(i,j)=3.0_R8/5.0_R8*sqrt(mu(i,j))
          N3(i,j)=-2.0_R8*(1.0_R8
     &    +3.0_R8/5.0_R8*amnucl(j)/amnucl(i))
     &    /(amnucl(i)+amnucl(j))*sqrt(mu(i,j))
          N4(i,j)=4.0_R8/5.0_R8
     &    /(amnucl(i)+amnucl(j))*sqrt(mu(i,j))
          N8(i,j)=-3.0_R8/14.0_R8*sqrt(mu(i,j))
          N5(i,j)=-(13.0_R8/10.0_R8*amnucl(j)/amnucl(i)+8.0_R8/5.0_R8
     &    +3.0_R8*amnucl(i)/amnucl(j))*kappa(i,j)*sqrt(mu(i,j))
          N6(i,j)=27.0_R8/10.0_R8
     &    *kappa(i,j)*sqrt(mu(i,j))
          N9(i,j)=3.0_R8/5.0_R8*(23.0_R8/28.0_R8*amnucl(j)/amnucl(i)
     &    +8.0_R8/7.0_R8+3.0_R8*amnucl(i)/amnucl(j))
     &    *kappa(i,j)*sqrt(mu(i,j))*mum(i,j)
          N10(i,j)=-45.0_R8/28.0_R8
     &    *kappa(i,j)*sqrt(mu(i,j))*mum(i,j)
          N11(i,j)=-(433.0_R8/280.0_R8*(amnucl(j)/amnucl(i))**2
     &    +136.0_R8/35.0_R8*amnucl(j)/amnucl(i)+
     &    +459.0_R8/35.0_R8+32.0_R8/5.0_R8*amnucl(i)/amnucl(j)+
     &    5.0_R8*(amnucl(i)/amnucl(j))**2)*kappa(i,j)**2*sqrt(mu(i,j))
          N12(i,j)=75.0_R8/8.0_R8
     &    *kappa(i,j)**2*sqrt(mu(i,j))
          N13(i,j)=(18.0_R8/35.0_R8*amnucl(j)/amnucl(i)
     &    +6.0_R8/5.0_R8)
     &    /amnucl(i)*kappa(i,j)*sqrt(mu(i,j))
          N14(i,j)=-24.0_R8/35.0_R8
     &    /amnucl(i)*kappa(i,j)*sqrt(mu(i,j))
          N15(i,j)=-(51.0_R8/35.0_R8*(amnucl(j)/amnucl(i))**2
     &    +37.0_R8/7.0_R8*amnucl(j)/amnucl(i)+
     &    +22.0_R8/5.0_R8+4.0_R8*amnucl(i)/amnucl(j))
     &    /(amnucl(i)+amnucl(j))*kappa(i,j)*sqrt(mu(i,j))
          N16(i,j)=24.0_R8/7.0_R8
     &    /mdm(i,j)/(amnucl(i)+amnucl(j))*kappa(i,j)*sqrt(mu(i,j))
        enddo
      enddo        
      call subend ()
      return

      end subroutine b2trcn
