      subroutine b2txcc (nCv, nFc, mode, geo, mpg, fcVol, fcS, kap, 
     .  con)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nFc, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   fcVol(nFc,2), fcS(nFc), kap(nCv), con(nFc)
*     ------------------------------------------------------------------
*     B2TXCC computes normal conductances between neighbouring cells, 
*     given the conductivities.
*     ------------------------------------------------------------------
      integer iFc
      real (kind=R8) ::
     *   t0, t1
      intrinsic max, sqrt
      external subini, subend, xertst
      external b2xxgs
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txcc')
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.3,
     &  'b2txcc--faulty argument mode')

*   ..compute con
      if (mode.eq.-1) then
*    ..use harmonic averaging
       do iFc = 1, nFc
        t0 = fcVol(iFc,1)*kap(mpg%fcCv(iFc,1))/geo%fcHc(iFc,1)**2
        t1 = fcVol(iFc,2)*kap(mpg%fcCv(iFc,2))/geo%fcHc(iFc,2)**2
        if (t0*t1.gt.0.0_R8) then
         con(iFc) = 
     &    t0*t1/(t0+t1)/geo%fcQgam(iFc,0)**2
        else
         con(iFc) = 0.0_R8
        endif
       enddo
      else if (mode.eq.0) then
*    ..use geometric averaging
       do iFc = 1, nFc
        con(iFc) = 
     &    1.0_R8/geo%fcQgam(iFc,0)**2*sqrt(max(0.0_R8,
     &   (fcVol(iFc,1)*kap(mpg%fcCv(iFc,1))/geo%fcHc(iFc,1))*
     &   (fcVol(iFc,2)*kap(mpg%fcCv(iFc,2))/geo%fcHc(iFc,2))))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))
       enddo
      else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
       do iFc = 1, nFc
        con(iFc) = 1.0_R8/geo%fcQgam(iFc,0)**2*
     &   (fcVol(iFc,1)*kap(mpg%fcCv(iFc,1)) + 
     &    fcVol(iFc,2)*kap(mpg%fcCv(iFc,2)))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))**2
       enddo
      else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
       do iFc = 1, nFc
        con(iFc) = 
     &    fcS(iFc)/geo%fcQgam(iFc,0)*
     &   (geo%fcHc(iFc,1)*kap(mpg%fcCv(iFc,1)) + 
     &    geo%fcHc(iFc,2)*kap(mpg%fcCv(iFc,2)))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))**2
       enddo
      else if (mode.eq.3) then
*    ..use arithmetic averaging (drift formulation)
       do iFc = 1, nFc
        con(iFc) = 
     &    fcS(iFc)/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))*
     &   (geo%fcVol(iFc,1)*kap(mpg%fcCv(iFc,1)) + 
     &    geo%fcVol(iFc,2)*kap(mpg%fcCv(iFc,2)))/
     &   (geo%fcVol(iFc,1) + geo%fcVol(iFc,2))
       enddo
      endif
*   ..reset null regions
!      call b2xxgs (nx, ny, 0.0_R8, con(-1,-1,0), 1)
!      call b2xxgs (nx, ny, 0.0_R8, con(-1,-1,1), 2)
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end
