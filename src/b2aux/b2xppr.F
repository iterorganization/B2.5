      subroutine b2xppr (nCv, ns, ne, na, te, ti, tn, pr, st_ext)
      use b2mod_types
      use b2us_plasma
      use b2mod_b2cmpa
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  ne(nCv), na(nCv,0:ns-1), te(nCv), ti(nCv), tn(nCv), pr(nCv)
*     ------------------------------------------------------------------
*     B2XPPR computes the total pressure, pr:
*       pr() = ne()*te() + (sum is :: na(,is)*ti()) .
*     ------------------------------------------------------------------
      integer is
*     ------------------------------------------------------------------
*$$$  call subini ('b2xppr')
      pr = ne*te
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          pr = pr + na(:,is)*ti
        else
          pr = pr + na(:,is)*tn
        endif
      enddo
      do is = 0, st_ext%ns-1
        pr = pr + st_ext%na(:,is)*st_ext%ta(:,is)
      enddo
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xppr
