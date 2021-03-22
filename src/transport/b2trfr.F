*-----------------------------------------------------------------------
*.specification

      subroutine b2trfr(nnucl, mum, rmu,
     & z2n,c_hwA,c_rwA, c_R_wB,c_R_wA,c_R_w)
      use b2mod_types
      use b2mod_subsys
      use b2mod_b2zhco, only: zhcscorr
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) :: 
     & z2n(0:nnucl-1),
     & mum(0:nnucl-1,0:nnucl-1), rmu(0:nnucl-1,0:nnucl-1)
      real (kind=R8) :: 
     & c_hwA(0:nnucl-1,0:nnucl-1), c_rwA(0:nnucl-1,0:nnucl-1),
     & c_R_wB(0:nnucl-1,0:nnucl-1)
*   ..output arguments 
      real (kind=R8) :: c_R_w(0:nnucl-1,0:nnucl-1),
     & c_R_wA(0:nnucl-1,0:nnucl-1)
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
*     - c_R_wA: the friction force transport coefficient for isonuclear ions
*     - c_R_w: the total friction force transport coefficient 
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,k
      real (kind=R8) :: temp, rel
      real (kind=R8) :: 
     & c_R_wA1, c_R_wA2, c_R_wA3

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfr')
*   ..main
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

*   ..apply charge state correction 
      c_R_w=c_R_wA
      do i = 0, nnucl-1
        if (zhcscorr.eq.1) then
          temp=0
          do k = 0, nnucl-1
            temp=temp+(rmu(i,k)/rmu(i,i))*(z2n(k)/z2n(i))*
     &   (c_R_wB(k,i)-c_R_wA(k,i))
          enddo
          c_R_w(i,i)=c_R_w(i,i)+temp
        else
          c_R_w(i,i)=1
        endif
      enddo
      call subend ()
      return

      end subroutine b2trfr
