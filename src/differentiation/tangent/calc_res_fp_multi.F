      subroutine calc_res_fp_multi(nbdirs, nCv, ns, tn_style, 
     & solve_keps, diagd, res_max)
      use b2mod_types
      use b2us_plasma_diffv
      use b2mod_numerics_namelist_diffv
      use B2MOD_PAR_OPT_DIFFV, only: par_rescale
      implicit none
      integer, intent (in) :: nCv, ns, nbdirs, tn_style, solve_keps 
      type(B2Diagnostic_diffv), intent(in) :: diagd
      real(kind=r8), intent (inout) :: res_max
*-----------------------------------------------------------------------
*     calc_res_fp_multi calculates a maximum residual for the vector   
*     tangent fixed point form of the B2.5 time iteration scheme,
*     normalized by the par_rescale factors from b2.optimization.parameters.
*-----------------------------------------------------------------------
      integer :: is, nd

      res_max = 0.0_R8
      
      do nd = 1, nbdirs
        do is = 0, ns-1
          if (any(solveco(is,:)))
     &       res_max = max(res_max,diagd%aresco(nd,is)/par_rescale(nd))
          if (any(solvemo(is,:)))
     &      res_max = max(res_max,diagd%aresmo(nd,is)/par_rescale(nd))
        end do
        if (any(solvepo))
     &    res_max = max(res_max,diagd%arespo(nd)/par_rescale(nd))
        if (any(solveee))
     &    res_max = max(res_max,diagd%areshe(nd)/par_rescale(nd))
        if (any(solveei))
     &    res_max = max(res_max,diagd%areshi(nd)/par_rescale(nd))
        if (any(solveen).and.tn_style.eq.2)
     &    res_max = max(res_max,diagd%areshn(nd)/par_rescale(nd))
        if (any(solvekt).and.solve_keps.ge.1)
     &    res_max = max(res_max,diagd%areskt(nd)/par_rescale(nd))
        if (any(solvezt).and.solve_keps.eq.2)
     &    res_max = max(res_max,diagd%areszt(nd)/par_rescale(nd))
      enddo

      return
      end subroutine calc_res_fp_multi
