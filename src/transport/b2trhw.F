*-----------------------------------------------------------------------
*.specification

      subroutine b2trhw(nnucl, G2, G8, mum, rmu,
     & z2n, Ainv, Ainv_mdf, c_hwA, c_rwA)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) ::
     & z2n(0:nnucl-1),
     & G2(0:nnucl-1,0:nnucl-1), G8(0:nnucl-1,0:nnucl-1),
     & mum(0:nnucl-1,0:nnucl-1), rmu(0:nnucl-1,0:nnucl-1)
      real (kind=R8) :: Ainv(0:2*nnucl-1,0:2*nnucl-1),
     & Ainv_mdf(0:nnucl-1,0:2*nnucl-1) 
*   ..output arguments 
      real (kind=R8) :: 
     & c_hwA(0:nnucl-1,0:nnucl-1), c_rwA(0:nnucl-1,0:nnucl-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2trhw computes the heat flux and r-moment velocity dependent part transport coefficient 
*
*  3. description 
*
*     This routine computes:
*     -c_hwA: the heat flux velocity dependent transport coefficient
*     -c_rwA: the r-moment velocity dependent transport coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,k
      real (kind=R8) :: temp, rel
      real (kind=R8) :: c_hwAm(0:nnucl-1,0:nnucl-1,0:nnucl-1),
     & c_rwAm(0:nnucl-1,0:nnucl-1,0:nnucl-1)

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trhw')
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

      call subend ()
      return


      end subroutine b2trhw


