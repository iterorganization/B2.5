      subroutine strcopy (l0, l1, str0, str1)
      implicit none
      integer l0, l1
      character str0*(l0), str1*(l1)
*     ------------------------------------------------------------------
*     STRCOPY copies the first l0 characters of string str1 into str0.
*     ------------------------------------------------------------------
      integer i
*     ------------------------------------------------------------------
      do i = 1, min(l0,l1)
        str0(i:i) = str1(i:i)
      end do
      return
*     ------------------------------------------------------------------
      end subroutine strcopy
