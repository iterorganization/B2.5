      subroutine b2xzmf (nCv, ns, ismain, rz2, na, zmff)             !srv 13.01.17 {
      use b2mod_types
      implicit none
      integer nCv, ns, ismain
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1), zmff(nCv)
*     ------------------------------------------------------------------
*     B2XZMF computes zmff:
*       zmff() = (sum over is /= ismain :: rz2(,is)*na(,is))/
*                  na(,ismain)*rz2(,ismain)
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzmf')
      zmff = 0.0_R8
      do is = 0, ns-1
       if(is.ne.ismain) then
         zmff = zmff + rz2(:,is)*na(:,is)
        endif
      enddo
      zmff=zmff/(rz2(:,ismain)*na(:,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzmf                                             !srv 13.01.17 }
