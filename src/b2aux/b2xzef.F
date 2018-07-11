      subroutine b2xzef (nCv, ns, rz2, na, ne, zeff, st_ext)                 !srv 13.01.17 {
      use b2mod_types
      use b2us_plasma
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  rz2(nCv,0:ns-1), na(nCv,0:ns-1), ne(nCv), zeff(nCv)
*     ------------------------------------------------------------------
*     B2XZEF computes zeff:
*       zeff() = (sum is :: rz2(,is)*na(,is))/ne .
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xzef')
      zeff = 0.0_R8
      do is = 0, ns-1
        zeff = zeff + rz2(:,is)*na(:,is)
      enddo
      do is = 0, st_ext%ns-1
        zeff = zeff + st_ext%za2(:,is)*st_ext%na(:,is)
      enddo
      zeff=zeff/ne
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xzef                                             !srv 13.01.17 }
