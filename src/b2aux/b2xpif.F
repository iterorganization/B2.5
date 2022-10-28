      subroutine b2xpif (nCv, nFc, switch, geo, mpg, fch_pi_c, fch_pi_f)      !srv 13.01.17
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nFc
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  fch_pi_c(nCv,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fch_pi_f(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPIF computes the ion parallel current on cell faces.
*     ------------------------------------------------------------------
      integer, save :: ncall = 0
*     ------------------------------------------------------------------

*  .. interpolate fch_pi_c(,0) to faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,fch_pi_c,fch_pi_f)

*  .. multiply with poloidal face area
      fch_pi_f(:,0) = fch_pi_f(:,0)*geo%fcS(:)*geo%fcQalf(:,0)
     &                                         *switch%prl_cur

*  .. set radial component zero
      fch_pi_f(:,1) = 0.0_R8

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpif
