      subroutine b2xppz (nCv, ns, ne, na, te, ti, pz, st_ext)         !sv 18.01.02 {
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_plasma
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  ne(nCv), na(nCv,0:ns-1), te(nCv), ti(nCv), pz(nCv)
*     ------------------------------------------------------------------
*     B2XPPZ computes pressure of charged particles, pz:
*       pz() = ne()*te() + (sum is :: na(,is)*ti() for za(is).ne.0)
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppz')
      pz = ne*te
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          pz = pz + na(:,is)*ti
        endif
      enddo
      do is = 0, st_ext%ns-1
        if (.not.st_ext%is_neutral(is)) then
          pz = pz + st_ext%na(:,is)*st_ext%ta(:,is)
        endif
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                        !sv 18.01.02 }
