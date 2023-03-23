      logical function esseq (x, y)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   x, y
*     ------------------------------------------------------------------
*     ESSEQ determines whether two real*8 numbers differ by no more than
*     a small multiple of their machine precision.
*     ------------------------------------------------------------------
      real (kind=R8) ::
     *   t0, t1, t2
      intrinsic abs
      real (kind=R8) ::
     *   eval
      external eval
*     ------------------------------------------------------------------
      t0 = eval(abs(x)+abs(y))
      t1 = eval(t0+abs(x-y)/16)
      t2 = eval(t0-abs(x-y)/16)
      esseq = t0.eq.t1.or.t0.eq.t2
*     (The purpose of comparing t0 both to t1 and to t2 is to deal with
*     unidirectional rounding.)
      return
*     ------------------------------------------------------------------
      end function esseq
