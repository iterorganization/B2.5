      subroutine b2mxnu (n, funp, fun0, funq)
      use b2mod_types
      use b2mod_math
      use b2mod_subsys
      implicit none
      integer n
      real (kind=R8) ::
     *   funp(0:n-1), fun0(0:n-1), funq(0:n-1)
*     ------------------------------------------------------------------
*     B2MXNU extrapolates an indefinite quantity to the next time.
*     It performs the operation
*       (funp,fun0,funq) = (fun0,funq,extrap(funp,fun0,funq)) ,
*     in which extrap refers to a certain form of extrapolation.
*     This is used to update the history over three timelevels of the
*     state of an indefinite (positive or negative) quantity fun.
*     The interpretation is that funp is the previous state, fun0 is
*     the state at the start of this timestep, and funq is the state at
*     the end of the timestep.
*     The new value of funq is a conservative prediction for fun at the
*     end of the next timestep.
*     ------------------------------------------------------------------
      integer i
      real (kind=R8) ::
     *   eps, t0, t1
      parameter (eps=0.0_R8)
      intrinsic abs
      external b2scopy
*     ------------------------------------------------------------------
*     Must set eps to a small non-negative value. eps.eq.0 is safest,
*     positive values of eps are suitable only if the inner iterations
*     are carried through to sufficient accuracy. If this is not done,
*     then even eps.eq.1.0e-2 can be disastrous.
*     ------------------------------------------------------------------
      call subini ('b2mxnu')
      if (eps.ne.0.0_R8) then
       do i = 0, n-1
        t0 = med(0.0_R8,funq(i)-fun0(i),fun0(i)-funp(i))
        t1 = max(abs(funp(i)),abs(fun0(i)),abs(funq(i)))
        funp(i) = fun0(i)
        fun0(i) = funq(i)
        funq(i) = funq(i)+med(-eps*t1,t0,eps*t1)
       enddo
      else
*      (faster code for common special case)
!      call b2scopy (n, fun0, 1, funp, 1)
       call b2scopy (n, funq, 1, funp, 1)
       call b2scopy (n, funq, 1, fun0, 1)
*      (leave funq unchanged)
      endif
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2mxnu
