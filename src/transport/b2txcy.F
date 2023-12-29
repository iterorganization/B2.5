! $AD NOCHECKPOINT
      subroutine b2txcy (nCv, nFc, mode, geo, mpg, fcVol, fcS, kapy, 
     .  cony)
      use b2mod_types
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      integer nCv, nFc, mode
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &   fcVol(nFc,2), fcS(nFc), kapy(nCv), cony(nFc)
*     ------------------------------------------------------------------
*     B2TXCY computes radial conductances between neighbouring cells,
*     given the conductivities.
*     ------------------------------------------------------------------
      integer iFc
      real (kind=R8) ::
     *   t0, t1
      intrinsic max, sqrt
      external xertst
*     ------------------------------------------------------------------
*   ..subprogram start-up calls
      call subini ('b2txcy')
*   ..test mode
      call xertst (-1.le.mode.and.mode.le.3,
     &  'b2txcy--faulty argument mode')

cwdk factor 0.5 (vs. 2.0 in original code)!
cwdk note that fcHc is typically hx/2 or hy/2!

*   ..compute cony
      if (mode.eq.-1) then
*    ..use harmonic averaging
       do iFc = 1, nFc
        t0 = fcVol(iFc,1)*kapy(mpg%fcCv(iFc,1))/geo%fcHc(iFc,1)**2
        t1 = fcVol(iFc,2)*kapy(mpg%fcCv(iFc,2))/geo%fcHc(iFc,2)**2
        if (t0*t1.gt.0.0_R8) then
         cony(iFc) = 
     &    t0*t1/(t0+t1)*geo%fcQalf(iFc,1)/geo%fcQgam(iFc,0)**2
        else
         cony(iFc) = 0.0_R8
        endif
       enddo
      else if (mode.eq.0) then
*    ..use geometric averaging
       do iFc = 1, nFc
        cony(iFc) = 
     &    geo%fcQalf(iFc,1)/geo%fcQgam(iFc,0)**2*sqrt(max(0.0_R8,
     &   (fcVol(iFc,1)*kapy(mpg%fcCv(iFc,1))/geo%fcHc(iFc,1))*
     &   (fcVol(iFc,2)*kapy(mpg%fcCv(iFc,2))/geo%fcHc(iFc,2))))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))
       enddo
      else if (mode.eq.1) then
*    ..use arithmetic averaging (SOLPS5.0 formulation)
       do iFc = 1, nFc
        cony(iFc) = geo%fcQalf(iFc,1)/geo%fcQgam(iFc,0)**2*
     &   (fcVol(iFc,1)*kapy(mpg%fcCv(iFc,1)) + 
     &    fcVol(iFc,2)*kapy(mpg%fcCv(iFc,2)))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))**2
       enddo
      else if (mode.eq.2) then
*    ..use arithmetic averaging (SOLPS4.0 formulation)
       do iFc = 1, nFc
        cony(iFc) = 
     &    fcS(iFc)*geo%fcQalf(iFc,1)/geo%fcQgam(iFc,0)*
     &   (geo%fcHc(iFc,1)*kapy(mpg%fcCv(iFc,1)) + 
     &    geo%fcHc(iFc,2)*kapy(mpg%fcCv(iFc,2)))/
     &   (geo%fcHc(iFc,1) + geo%fcHc(iFc,2))**2
       enddo
      else if (mode.eq.3) then
*    ..use arithmetic averaging (drift formulation)
       do iFc = 1, nFc
        cony(iFc) = 
     &    fcS(iFc)*geo%fcQalf(iFc,1)/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))*
     &   (geo%fcVol(iFc,1)*kapy(mpg%fcCv(iFc,1)) + 
     &    geo%fcVol(iFc,2)*kapy(mpg%fcCv(iFc,2)))/
     &   (geo%fcVol(iFc,1) + geo%fcVol(iFc,2))
       enddo
      endif
*   ..return
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2txcy
