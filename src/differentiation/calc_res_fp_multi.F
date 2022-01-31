      subroutine calc_res_fp_multi(nbdirs, nCv, ns, stated, stated0,
     & res)
      use b2mod_types
      use b2us_plasma_diffv
      implicit none
      integer, intent (in) :: nCv, ns, nbdirs
      type(b2state_diffv), intent(in) :: stated
      type(b2state_diffv), intent(in) :: stated0
      real(kind=r8), intent (inout) :: res
      real(kind=r8) :: t0
      integer :: icv, is, nd

      res = 0.0_R8
      
      do nd = 1, nbdirs
        do is = 0, ns-1
          t0 = sqrt(sum((stated%pl%na(nd,:,is)**2)))
          if (t0 .gt. 0.0_R8) res = max(res,
     &     sqrt(sum((stated%pl%na(nd,:,is)-stated0%pl%na(nd,:,is))**2))
     &     /t0)
          t0 = sqrt(sum((stated%pl%ua(nd,:,is)**2)))
          if (t0 .gt. 0.0_R8) res = max(res,
     &     sqrt(sum((stated%pl%ua(nd,:,is)-stated0%pl%ua(nd,:,is))**2))
     &     /t0)
        enddo
        t0 = sqrt(sum((stated%pl%te(nd,:)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &   sqrt(sum((stated%pl%te(nd,:)-stated0%pl%te(nd,:))**2))/t0)
        t0 = sqrt(sum((stated%pl%ti(nd,:)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%ti(nd,:)-stated0%pl%ti(nd,:))**2))/t0)
        t0 = sqrt(sum((stated%pl%tn(nd,:)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%tn(nd,:)-stated0%pl%tn(nd,:))**2))/t0)
        t0 = sqrt(sum((stated%pl%po(nd,:)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%po(nd,:)-stated0%pl%po(nd,:))**2))/t0)
        t0 = sqrt(sum((stated%pl%kt(nd,:)**2)))
        if (t0 .gt. 0.0_R8) res = max(res,
     &    sqrt(sum((stated%pl%kt(nd,:)-stated0%pl%kt(nd,:))**2))/t0)
      enddo

      return
      end subroutine calc_res_fp_multi



