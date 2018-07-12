      subroutine b2xpro (nCv, ns, mp, am, na, ro, st_ext)
      use b2mod_types
      use b2us_plasma
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *   mp, am(0:ns-1), na(nCv,0:ns-1), ro(nCv)
*     ------------------------------------------------------------------
*     B2XPRO computes the mass density, ro:
*       ro() = (sum is :: am(is)*mp*na(,is)) .
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpro')
      ro = 0.0_R8
      do is = 0, ns-1
        ro = ro + mp*am(is)*na(:,is)
      enddo
      do is = 0, st_ext%ns-1
        ro = ro + mp*st_ext%am(is)*st_ext%na(:,is)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpro
