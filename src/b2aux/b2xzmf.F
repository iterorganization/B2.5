      subroutine b2xzmf (nx, ny, ns, ismain, rz2, na, zmff)             !srv 13.01.17 {
      use b2mod_types
      implicit none
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *  rz2(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1),
     *  zmff(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2XZMF computes zmff:
*       zmff(,) = (sum over is /= ismain :: rz2(,,is)*na(,,is))/
*                  na(,,ismain)*rz2(,,ismain)
*     ------------------------------------------------------------------
      integer n2, is
      external sfill
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzmf')
      n2 = (nx+2)*(ny+2)
      call sfill (n2, 0.0_R8, zmff, 1)
      do is = 0, ns-1
       if(is.ne.ismain) then
         zmff(-1:nx,-1:ny) = zmff(-1:nx,-1:ny) +
     .       rz2(-1:nx,-1:ny,is)*na(-1:nx,-1:ny,is)
        endif
      enddo
      zmff=zmff/(rz2(-1:nx,-1:ny,ismain)*na(-1:nx,-1:ny,ismain))
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzmf                                             !srv 13.01.17 }
