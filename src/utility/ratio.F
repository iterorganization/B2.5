      function rratio (m, n)
      use b2mod_types
      implicit none
      real (kind=R8) ::
     *   rratio
      integer m, n
*     ------------------------------------------------------------------
*     RATIO evaluates very carefully the real*8 ratio of two integers.
*     rratio(m,n) evaluates the mathematical expression real(m)/real(n).
*     It is ensured that the expressions rratio(m,n) and rratio(k*m,k*n)
*     give identical numerical results for any nonzero integer k.
*     ------------------------------------------------------------------
      integer m0, n0, k
      intrinsic max, abs, mod, real
      external xertst
*     ------------------------------------------------------------------
      call xertst (n.ne.0, 'ratio--error, n.eq.0')
*     (compute  k = gcd(abs(m),abs(n)))
      m0 = abs(m)
      n0 = abs(n)
    1 if (m0.ne.0.and.m0.le.n0) then
       n0 = mod(n0,m0)
       goto 1
      else if (n0.ne.0.and.n0.le.m0) then
       m0 = mod(m0,n0)
       goto 1
      endif
      k = max(m0,n0)
*     (compute ratio)
      rratio = real(m/k)/real(n/k)
      return
*     ------------------------------------------------------------------
      end
