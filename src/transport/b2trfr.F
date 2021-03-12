*-----------------------------------------------------------------------
*.specification

      subroutine b2trfr(nnucl, G2, G8, mum, rmu,
     & z2n, Ainv, Ainv_mdf, c_R_w)
      use b2mod_types
      use b2mod_subsys
!     use b2mod_b2zhco, only: c_R_wA
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) :: amti(0:nnucl-1),
     & z2n(0:nnucl-1),
     & G2(0:nnucl-1,0:nnucl-1), G8(0:nnucl-1,0:nnucl-1),
     & mum(0:nnucl-1,0:nnucl-1), rmu(0:nnucl-1,0:nnucl-1)
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1),
     & Ainv_mdf(0:nnucl-1,0:2*nnucl-1) 
*   ..output arguments 
      real (kind=R8) :: c_R_w(0:nnucl-1,0:nnucl-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRFR computes the friction force transport coefficient 
*
*  3. description 
*
*     This routine computes:
*     -c_R_wA: the friction force transport coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,k
      real (kind=R8) :: temp, rel
      real (kind=R8) :: c_hwAm(0:nnucl-1,0:nnucl-1,0:nnucl-1),
     & c_rwAm(0:nnucl-1,0:nnucl-1,0:nnucl-1),
     & c_hwA(0:nnucl-1,0:nnucl-1), c_rwA(0:nnucl-1,0:nnucl-1),
     & c_R_wA1, c_R_wA2, c_R_wA3,
     & c_R_wA(0:nnucl-1,0:nnucl-1)

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfr')
*   ..main
*   ..components of the velocity dependent parts  
      do j = 0, nnucl-1
        do k = 0, nnucl-1
          do i = 0, nnucl-1
*   ..of the heat flux
            c_hwAm(i,k,j)=-2.5_R8*mum(k,j)*G2(k,j)*Ainv(i,k)-
     & 17.5_R8*mum(k,j)**2*G8(k,j)*Ainv(i,k+nnucl)           
*   ..of the r-moment
            c_rwAm(i,k,j)=-2.5_R8*mum(k,j)*G2(k,j)*Ainv_mdf(i,k)-
     & 17.5_R8*mum(k,j)**2*G8(k,j)*Ainv_mdf(i,k+nnucl)           
          enddo
        enddo
      enddo
 
!     do k = 0, nnucl-1
!     do i=0,nnucl-1
!     print *, (c_hwAm(i,j,k), j=0,nnucl-1 )
!     enddo
!     print *
!     enddo
!     do k = 0, nnucl-1
!     do i=0,nnucl-1
!     print *, (c_rwAm(i,j,k), j=0,nnucl-1 )
!     enddo
!     print *
!     enddo

*   ..velocity dependent parts 
      c_hwA=0
      c_rwA=0
      do i = 0, nnucl-1
        do j = 0, nnucl-1
          do k = 0, nnucl-1
*   ..of the heat flux
            c_hwA(j,i)=c_hwA(j,i)+c_hwAm(i,k,j)-c_hwAm(i,j,k)
*   ..of the r-moment 
            c_rwA(j,i)=c_rwA(j,i)+c_rwAm(i,k,j)-c_rwAm(i,j,k)
          enddo
        enddo
      enddo

!     do i=0,nnucl-1
!     print *, (c_hwA(i,j), j=0,nnucl-1 )
!     enddo
!     print *
!     do i=0,nnucl-1
!     print *, (c_rwA(i,j), j=0,nnucl-1 )
!     enddo
!     print *

*   ..the friction force transport coefficient
      rel=3.0_R8/14.0_R8
      do i = 0, nnucl-1
        do j = 0, nnucl-1
          c_R_wA1=0
          c_R_wA2=0
          c_R_wA3=0
          do k = 0, nnucl-1
            if (k.ne.i) then
              temp=(rmu(i,k)/rmu(i,j))*(z2n(k)/z2n(j))
*   ..first part (main for the light species)
              c_R_wA1=c_R_wA1-mum(i,k)*temp*
     & (0.6_R8*c_hwA(j,i)-rel*mum(i,k)*c_rwA(j,i))
*   ..second part (main for the heavy species)
              c_R_wA2=c_R_wA2+mum(k,i)*temp*
     & (0.6_R8*c_hwA(j,k)-rel*mum(k,i)*c_rwA(j,k))
            endif
*   ..third part (main for the heavy species)
            c_R_wA3=c_R_wA3-mum(j,i)*
     & (0.6_R8*c_hwA(k,j)-rel*mum(j,i)*c_rwA(k,j))
          enddo
*   ..finally
          c_R_wA(j,i)=1+c_R_wA1+c_R_wA2+c_R_wA3
        enddo
      enddo

!     print *
!     do, i=0,nnucl-1
!     print *, (c_R_wA(ix,iy,i,j), j=0,nnucl-1 )
!     enddo
!     print *
        
      c_R_w=c_R_wA
!     do i = 0, nnucl-1
!       temp=0
!       do k = 0, nnucl-1
!         temp=temp+(rmu(i,k)/rmu(i,i))*(z2n(k)/z2n(i))*
!    & (1-c_R_wA(k,i))
!       enddo
!       c_R_w(i,i)=c_R_w(i,i)+temp
!     enddo

      do i = 0, nnucl-1
        c_R_w(i,i)=1
      enddo
      call subend ()
      return


      end subroutine b2trfr
