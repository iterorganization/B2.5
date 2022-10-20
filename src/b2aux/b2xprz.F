      subroutine b2xprz (nCv, ns, mp, am, na, rz, st_ext)
      use b2mod_types
      use b2mod_b2cmpa
     * , only : is_neutral
      use b2us_plasma
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *   mp, am(0:ns-1), na(nCv,0:ns-1), rz(nCv)
*     ------------------------------------------------------------------
*     B2XPRZ computes the mass density of the charged particles, rz:
*       rz() = (sum is :: am(is)*mp*na(,is)) .
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xprz')
      rz = 0.0_R8
      do is = 0, ns-1
       if (.not.is_neutral(is)) then
        rz = rz + mp*am(is)*na(:,is)
       endif
      enddo
      do is = 0, st_ext%ns-1
       if (.not.st_ext%is_neutral(is)) then
        rz = rz + mp*st_ext%am(is)*st_ext%na(:,is)
       endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xprz
