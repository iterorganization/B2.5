      subroutine calc_res_fp_multi(nbdirs, nCv, ns, tn_style, 
     & solve_keps, diagd, res_max)
      use b2mod_types
      use b2us_plasma_diffv
      use b2mod_numerics_namelist_diffv
      implicit none
      integer, intent (in) :: nCv, ns, nbdirs, tn_style, solve_keps 
      type(B2Diagnostic_diffv), intent(in) :: diagd
      real(kind=r8), intent (inout) :: res_max
*-----------------------------------------------------------------------
*     calc_res_fp_multi calculates a maximum residual for the vector   
*     tangent fixed point form of the B2.5 time iteration scheme.
*-----------------------------------------------------------------------
      integer :: is, nd

      res_max = 0.0_R8
      
      do nd = 1, nbdirs
        do is = 0, ns-1
          if (any(solveco(is,:)))
     &       res_max = max(res_max,diagd%aresco(nd, is))
          if (any(solvemo(is,:)))
     &      res_max = max(res_max,diagd%aresmo(nd, is))
        end do
        if (any(solvepo))
     &    res_max = max(res_max,diagd%arespo(nd))
        if (any(solveee))
     &    res_max = max(res_max,diagd%areshe(nd))
        if (any(solveei))
     &    res_max = max(res_max,diagd%areshi(nd))
        if (any(solveen).and.tn_style.eq.2)
     &    res_max = max(res_max,diagd%areshn(nd))
        if (any(solvekt).and.solve_keps.ge.1)
     &    res_max = max(res_max,diagd%areskt(nd))
        if (any(solvezt).and.solve_keps.eq.2)
     &    res_max = max(res_max,diagd%areszt(nd))
      enddo

      return
      end subroutine calc_res_fp_multi

