      subroutine calc_res_fp(ns, state, state0, res)
      use b2mod_types
      use b2us_plasma
      implicit none
      integer, intent (in) :: ns
      type(b2state), intent(in) :: state
      type(b2state), intent(in) :: state0
      real(kind=r8), intent (inout) :: res
*-----------------------------------------------------------------------
*     calc_res_fp calculates a maximum residual for the fixed point
*     form of the B2.5 time iteration scheme (while-loop instad of
*     for-loop). This is evaluated as the maximum among value of
*     t0 = sqrt(sum((q-q0)**2))/sqrt(sum(q**2)), where q is a plasma
*     state variable at time t and q0 the same variable at t-dt.
*-----------------------------------------------------------------------
      real(kind=r8) :: t0
      integer :: is

      res = 0.0_R8

      do is = 0, ns-1
        t0 = sqrt(sum((state%pl%na(:,is)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &   sqrt(sum((state%pl%na(:,is)-state0%pl%na(:,is))**2))/t0)
        t0 = sqrt(sum((state%pl%ua(:,is)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((state%pl%ua(:,is)-state0%pl%ua(:,is))**2))/t0)
      enddo
      t0 = sqrt(sum((state%pl%te**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     & sqrt(sum((state%pl%te-state0%pl%te)**2))/t0)
      t0 = sqrt(sum((state%pl%ti**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((state%pl%ti-state0%pl%ti)**2))/t0)
      t0 = sqrt(sum((state%pl%tn**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((state%pl%tn-state0%pl%tn)**2))/t0)
      t0 = sqrt(sum((state%pl%po**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((state%pl%po-state0%pl%po)**2))/t0)
      t0 = sqrt(sum((state%pl%kt**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((state%pl%kt-state0%pl%kt)**2))/t0)

      return
      end subroutine calc_res_fp



