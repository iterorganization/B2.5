      subroutine b2xpve (nCv, nFc, ns, switch, geo, mpg, qe, rza,
     &                rz2, alfx_c,                                     !srv 17.11.16 {
     &                sigx_c, fch_pi_c, po, ne, te, na, ua, st_ext, ve)
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none
      integer, intent (in) :: nCv, nFc, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     &   qe, rza(nCv,0:ns-1), rz2(nCv,0:ns-1), alfx_c(nCv), sigx_c(nCv),
     &   te(nCv), po(nCv), ne(nCv), na(nCv,0:ns-1), ua(nCv,0:ns-1),
     &   fch_pi_c(nCv,0:1)
      real (kind=R8) ::
     *   ve(nCv)
*     ------------------------------------------------------------------
*     B2XPVE computes the parallel electron velocity on cell centers.
*     ------------------------------------------------------------------
      integer is
      integer, save :: ncall = 0                                        !srv 03.04.15
      real (kind=R8) ::
     &   wrksig(nFc), wrkalf(nFc), wrksig_c(nCv), wrkalf_c(nCv),
     &   gpop(nFc), gpep(nFc), gtep(nFc), dumm(mpg%nVx), nef(nFc),
     &   zeff(nCv), wght(nFc,2), pe(nCv)
*     ------------------------------------------------------------------
      if(switch%b2sigp_style.eq.2) then                                 !srv 13.01.17 {
        call b2xzef (nCv, ns, rz2, na, ne, zeff, st_ext)
        call b2xpic (nCv, ns, qe, rza, geo%cvBb, na, ua, zeff,
     &    st_ext, fch_pi_c)
      else
        fch_pi_c = 0.0_R8
      endif                                                             !srv 13.01.17 }
*   ..compute ve, electron parallel velocity                            !iyv 19.11.13 {
      wght = 1.0_R8
      pe = ne*te
      call grad_p(nCv, nFc, mpg%nVx, 0, geo, mpg, po, dumm, gpop)
      call grad_p(nCv, nFc, mpg%nVx, 0, geo, mpg, pe, dumm, gpep)
      call grad_p(nCv, nFc, mpg%nVx, 0, geo, mpg, te, dumm, gtep)
      call intface(nCv,nFc,mpg%fcCv,wght,ne,nef)
      wrksig = -gpop + gpep/(qe*nef)
      wrkalf =  gtep

      call intcell(nFc,nCv,mpg,mpg%intcellP,wrksig,wrksig_c)
      call intcell(nFc,nCv,mpg,mpg%intcellP,wrkalf,wrkalf_c)
      ve = -sigx_c*wrksig_c + switch%b2sifr_phm2*alfx_c*wrkalf_c

      ve = (ve+fch_pi_c(:,0))/(geo%cvBb(:,0)/geo%cvBb(:,3))/qe*
     &                                               switch%prl_cur      !srv 21.12.15
      do is=0,ns-1
        ve = ve + rza(:,is)*ua(:,is)*na(:,is)
      enddo
      ve = ve + st_ext%ue*st_ext%ne
      ve = ve/ne

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpve                                             !srv 17.11.16 }
