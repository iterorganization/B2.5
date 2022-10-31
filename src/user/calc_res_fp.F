      subroutine calc_res_fp(nCv, ns, state, state0, res)
      use b2mod_types
      use b2us_plasma
      implicit none
      integer, intent (in) :: nCv, ns
      type(b2state), intent(in) :: state
      type(b2state), intent(in) :: state0
      real(kind=r8), intent (inout) :: res
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



