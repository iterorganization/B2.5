      subroutine b2xpne (nCv, ns, rza, na, ne_ext, ne)
      use b2mod_types
      implicit none
      integer nCv, ns
      real (kind=R8) ::
     *  rza(1:nCv,0:ns-1), na(1:nCv,0:ns-1), ne_ext(1:nCv), ne(1:nCv)
*     ------------------------------------------------------------------
*     B2XPNE computes the electron density, ne:
*       ne() = (sum is :: rza(,is)*na(,is)) .
*     A contribution from external charged species not followed by B2
*     can be included via use of the ne_ext array.
*     ------------------------------------------------------------------
      integer is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpne')
      call sfill (nCv, 0.0_R8, ne, 1)
      do is = 0, ns-1
       ne = ne + rza(:,is)*na(:,is)
      enddo
      ne = ne + ne_ext
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpne
