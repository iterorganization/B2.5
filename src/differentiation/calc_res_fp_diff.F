      subroutine calc_res_fp_diff(nCv, ns, stated, stated0, res)
      use b2mod_types
      use b2us_plasma_diff
      implicit none
      integer, intent (in) :: nCv, ns
      type(b2state_diff), intent(in) :: stated
      type(b2state_diff), intent(in) :: stated0
      real(kind=r8), intent (inout) :: res
      real(kind=r8) :: t0
      integer :: icv, is

      res = 0.0_R8
      
      do is = 0, ns-1
        t0 = sqrt(sum((stated%pl%na(:,is)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &   sqrt(sum((stated%pl%na(:,is)-stated0%pl%na(:,is))**2))/t0)
        t0 = sqrt(sum((stated%pl%ua(:,is)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%ua(:,is)-stated0%pl%ua(:,is))**2))/t0)
      enddo
      t0 = sqrt(sum((stated%pl%te**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     & sqrt(sum((stated%pl%te-stated0%pl%te)**2))/t0)
      t0 = sqrt(sum((stated%pl%ti**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((stated%pl%ti-stated0%pl%ti)**2))/t0)
      t0 = sqrt(sum((stated%pl%tn**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((stated%pl%tn-stated0%pl%tn)**2))/t0)
      t0 = sqrt(sum((stated%pl%po**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &  sqrt(sum((stated%pl%po-stated0%pl%po)**2))/t0)
      t0 = sqrt(sum((stated%pl%kt**2)))
      if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%kt-stated0%pl%kt)**2))/t0)

      return
      end subroutine calc_res_fp_diff



