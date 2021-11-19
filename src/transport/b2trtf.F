*-----------------------------------------------------------------------
*.specification

      subroutine b2trtf (nnucl, G2, G8, mum, mdm,
     &  Ainv, Ainv_mdf, c_R_TA)                                         !som 05.02.2021
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nnucl
      real (kind=R8) ::
     &  amti(0:nnucl-1), nz2n(0:nnucl-1),
     &  G2(0:nnucl-1,0:nnucl-1), G8(0:nnucl-1,0:nnucl-1),
     &  mum(0:nnucl-1,0:nnucl-1), mdm(0:nnucl-1,0:nnucl-1)
      real (kind=R8) ::
     &  Ainv(0:2*nnucl-1,0:2*nnucl-1),
     &  Ainv_mdf(0:nnucl-1,0:2*nnucl-1)
*   ..output arguments
      real (kind=R8) ::
     &  c_R_TA(0:nnucl-1,0:nnucl-1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRTF computes the thermal force transport coefficient.
*
*  3. description
*
*     This routine computes:
*     -c_R_TA: the thermal force transport coefficient
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer :: i,j,k
      real (kind=R8) :: temp
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trtf')
*   ..main
      do j = 0, nnucl-1
        do i = 0, nnucl-1
          temp = 0.0_R8
          do k = 0, nnucl-1
            temp=temp-2.5_R8*(mum(i,k)*G2(i,k)*(Ainv(i,j)-
     &       mdm(i,k)*Ainv(k,j))+mum(i,k)**2*G8(i,k)*(Ainv_mdf(i,j)-
     &       mdm(i,k)**2*Ainv_mdf(k,j)))
          enddo
          c_R_TA(i,j)=temp
        enddo
      enddo

      call subend ()
      return

      end subroutine b2trtf
