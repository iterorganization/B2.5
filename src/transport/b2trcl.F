*-----------------------------------------------------------------------
*.specification

      subroutine b2trcl (nCv, nFc, nVx, ns, nscx, iscx, ismain,          !srv 23.11.10
     1  geo, mpg, pl, dv, rt, co, st_ext, cvsahz,                        !srv 09.01.01 03.06.03
     3  cvsa, chce, chve, chci, chvi, csig, calf, csigin,                !xpb
     4  cthe, cthi, vsaf_cl, cvsa_cl, cvsahz_cl, fllimvisc,                       !srv 09.01.01
     5  csig_cl, calf_cl)                                                !xpb
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, ismain
      integer iscx(0:nscx-1)
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2Coeff), intent (in) :: co
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
cwdk to consider: create an output struct co0 of type (B2Coeff) for this?
      real (kind=R8) ::
     &  cvsa(nFc,0:1,0:ns-1),
     &  chce(nFc,0:1), chve(nFc,0:1),
     &  chci(nFc,0:1), chvi(nFc,0:1),
     &  csig(nFc,0:1), calf(nFc,0:1),
     &  csigin(nFc,0:1,0:ns-1,0:ns-1),				 !xpb
     &  cthe(nFc,0:ns-1), cthi(nFc,0:ns-1),
     &  vsaf_cl(nFc,0:1,0:ns-1),
     &  cvsa_cl(nFc,0:1,0:ns-1),
     &  cvsahz_cl(nFc,0:1,0:ns-1),                               !srv 09.01.01
     &  fllimvisc(nFc,0:ns-1), 
     &  cvsahz(nFc,0:1,0:ns-1),   !srv 09.01.01
     &  csig_cl(nFc,0:1), calf_cl(nFc,0:1),
     &  fllim_ke(nFc), fllim_ki(nFc), fllim_al(nFc), fllim_al_c(nCv),
     &  f_luc_ke(nFc), f_luc_ki(nFc), f_luc_et(nFc), f_luc_sg(nFc), 
     &  f_luc_al(nFc), alfx_c(nCv), sigx_c(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRCL computes the classical transport coefficients, including
*     geometric factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete classical transport
*     coefficients.
*     - cvsa: viscosity.
*     - chce: electron thermal conductivity.
*     - chve: electron thermal strange velocity.
*     - chci: atom thermal conductivity.
*     - chvi: atom thermal strange velocity.
*     - csig: electrical conductivity.
*     - csigin: ion-neutral conductivity.
*     - calf: thermo-electric coefficient.
*
*     The local differential transport coefficients, defined on cells,
*     are obtained by calls to b2tqca and b2tqce. To compute from these
*     the discrete transport coefficients involves multiplication by
*     geometric factors (vol, vol/hx**2, vol/hy**2), and, usually,
*     interpolation to cell faces.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, k, is, is0, lthf21, lvis21, lluciani, iCv, iFc
      integer, save :: conductive_limit = 1                            !srv 16.06.08
      integer, save :: style = 2                                       !srv 13.01.17
      integer, save :: iout = 0, iout_b2wdat = 0,                      !srv 17.06.08
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0,
     &                 integral_modeling = 1, tn_style = 0
      character(3) chns, chns0                                        !srv 17.06.08
      integer mintopcut                                               !lk 19.09.07  
      integer mode
      real (kind=R8) ::
     *  cflmv, cthiv, cthev, vis_par, jemax                           !srv 13.01.17
      real (kind=R8) ::                                               !srv 03.06.03
     *   t0, dte(nFc), dti(nFc), flomx
      real (kind=R8) :: tauia(nCv,0:ns-1), nu2(nCv),  !lk 19.09.07
     &    k2(nCv), tt0(mpg%nFs), tt1(mpg%nFs), eps1(mpg%nFs)          !lk 19.09.07
      real (kind=R8) :: tau(nCv), nu1(nCv),                           !lk 12.05.11
     &    k1(nCv), t1, t2                                             !lk 12.05.11
      real (kind=R8), save :: cvsa_mltpl = 1.0_R8                     !iyv 06.12.13
      real (kind=R8), save :: cutlo = 1.0e-10_R8
      real (kind=R8) ::
     *  hce(nCv,0:1), sig(nCv,0:1),
     1  alf(nCv,0:1), vsa(nCv,0:1,0:ns-1),
     2  hci(nCv,0:1), sigin(nCv,0:1,0:ns-1,0:nscx-1),
     3  wrk(nFc), wrk1(nFc,2), wrkv(nVx)                              !srv 01.07.99 09.01.01
csw 13mar2012
      integer, save :: core_cond_limit=0
*   ..procedures
      real (kind=R8) ::
     *  smin, smax                                                    !srv 22.04.16
      external subini, subend, xertst, xerrab, sfill, smin, smax
      external b2tqca, b2tqce, b2tqin, b2txcx, b2txcy, b2tlmv,        !xpb
     &  ipgeti, ipgetr, b2ttia                                        !lk 19.09.07
*   ..initialisation
      save ncall, lthf21, lvis21, cthiv, cthev, lluciani, mode
      data ncall/0/, lthf21/0/, lvis21/0/, mode/1/,
     1  cthiv/2.65_R8/, cthev/0.0_R8/, vis_par/0.0_R8/,
     1  lluciani/3/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trcl')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tral_mode',mode)
        call ipgeti ('b2trcl_lthf21',lthf21)
        call ipgeti ('b2trcl_lvis21',lvis21)
#ifdef NO_21
        call xertst (lthf21.eq.0.and.lvis21.eq.0,
     &   'Recompile without NO_21 option for 21-moment approach!')
#endif
        call ipgetr ('b2tfhe_vis_par',vis_par)
        call xertst (.not.(lvis21.ne.0.and.vis_par.ne.0.0_R8),
     &   '21-moment viscosity already includes parallel viscosity!')
        call ipgeti ('b2trcl_lluciani',lluciani)
        call ipgetr ('b2trcl_cthi', cthiv)
        call xertst (0.0_R8.le.cthiv,
     &	 'faulty internal parameter cthiv')
        call ipgetr ('b2trcl_cthe', cthev)
        call xertst (0.0_R8.le.cthev,
     &   'faulty internal parameter cthev')
        call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)                !iyv 06.12.13
        call ipgeti ('b2trcl_conductive_limit', conductive_limit)    !srv 16.06.08
        call ipgeti ('b2trcl_iout', iout)                            !srv 17.06.08
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 11.09.09 }
csw 13mar2012
        call ipgeti ('b2trcl_core_cond_limit', core_cond_limit)
        call ipgeti ('b2trcl_integral_modeling', integral_modeling)  !srv 01.11.11
        call ipgeti ('b2sigp_style', style)                          !srv 13.01.17
        call ipgeti ('b2mn_tn_style', tn_style)                         
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
        call b2xvsg (nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*    ..test sign of vol
        call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range fcQgam(,0)')
*    ..test sign of na, te, ti
        call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
        call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
        call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
      endif

* ..compute transport coefficients
*   ..compute differential classical transport coefficients
      call b2tqca (nCv, ns, geo, pl, dv, rt, st_ext, vsa, hci)        !srv 01.07.09 20.09.11
      call b2tqce (nCv, geo, mpg, pl, dv, hce, sig, alf)              !srv 01.07.09 20.09.11
      call b2tqin (nCv, ns, nscx, iscx, geo, pl, rt, sigin)           !xpb
      call sfill (nFc*ns, cthev, cthe,1)
      call sfill (nFc*ns, cthiv, cthi,1)
      call sfill (nFc, 1.0_R8, fllim_ki, 1)                            !srv 13.01.17 {
      call sfill (nFc, 1.0_R8, fllim_ke, 1)
      call sfill (nFc, 1.0_R8, fllim_al, 1)
      call sfill (nCv, 1.0_R8, fllim_al_c, 1)
      alfx_c(:) = alf(:,0)  !/abs((bb(:,:,0))/bb(:,:,3))          !iys 19.11.13 !srv 08.10.14 { !srv 13.01.17
      sigx_c(:) = sig(:,0)  !/abs((bb(:,:,0))/bb(:,:,3))          !             }               !srv 13.01.17
*   ..interpolate vsa to faces
      !wdk need to write general routine for volume interpolation like this...
      !wdk intface_v? intface with mode switch?
      do is = 0, ns-1
        do iFc = 1, nFc
          vsaf_cl(iFc,0,is) =
     &      (geo%fcVol(iFc,1)*vsa(mpg%fcCv(iFc,1),0,is)+
     &       geo%fcVol(iFc,2)*vsa(mpg%fcCv(iFc,2),0,is))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
          vsaf_cl(iFc,1,is) =
     &      (geo%fcVol(iFc,1)*vsa(mpg%fcCv(iFc,1),1,is)+
     &       geo%fcVol(iFc,2)*vsa(mpg%fcCv(iFc,2),1,is))/
     &      (geo%fcVol(iFc,1)+geo%fcVol(iFc,2))
        enddo
      enddo
*   ..incorporate geometry with hz                                    !srv 09.01.01
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,geo%cvHz,wrk)
      wrk = geo%fcS*wrk 
      do iFc = 1, nFc
        wrk1(iFc,1) = geo%fcVol(iFc,1)*geo%cvHz(mpg%fcCv(iFc,1))
        wrk1(iFc,2) = geo%fcVol(iFc,2)*geo%cvHz(mpg%fcCv(iFc,2))
      enddo
      do is = 0, ns-1
       call b2txcx (nCv, nFc, mode, geo, mpg, wrk1, wrk,              !srv 01.07.05 {
     &   vsa(1,0,is), cvsahz(1,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, wrk1, wrk,
     &   vsa(1,1,is), cvsahz(1,1,is))                                 !srv 01.07.05 }
      enddo                                                           !srv 09.01.01
*   ..incorporate geometry
      call sfill (nFc*2*ns*ns, 0.0_R8, csigin(1,0,0,0), 1)
      do is = 0, ns-1
       call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   vsa(1,0,is), cvsa(1,0,is))
       call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &   vsa(1,1,is), cvsa(1,1,is))
       do k = 0, nscx-1
         is0 = iscx(k)
         if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
           call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,   !sv
     &       sigin(1,0,is,k), csigin(1,0,is,is0))                       !sv 01.07.99
           call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,   !sv
     &       sigin(1,1,is,k), csigin(1,1,is,is0))                       !sv
        endif
       enddo
      enddo

      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  hce(1,0), chce(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  hce(1,1), chce(1,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  hci(1,0), chci(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  hci(1,1), chci(1,1))
      call b2txcx (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  sig(1,0), csig(1,0))
      call b2txcy (nCv, nFc, mode, geo, mpg, geo%fcVol, geo%fcS,
     &  sig(1,1), csig(1,1))
      call b2txcx (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS,
     &  alf(1,0), calf(1,0))
      call b2txcy (nCv, nFc,   -1, geo, mpg, geo%fcVol, geo%fcS,
     &  alf(1,1), calf(1,1))

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        if (.not. is_neutral(is)) then
          call my_out_us(70,nFc,1,cvsahz(1,0,is),'b2trcl_cvsahzx'//
     &      chns)
        endif
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                   .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),'b2trcl_chcex')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_chcey')
        call my_out_us(70,nFc,1,calf(1,0),'b2trcl_calfx')
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_calfy')
        call my_out_us(70,nFc,1,alfx_c,'b2trcl_alfx_c')                    !srv 05.07.17
        call my_out_us(70,nFc,1,chci(1,0),'b2trcl_chcix')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_chciy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                   .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,csig(1,0),'b2trcl_csigx')
!        call my_out_us(70,nFc,1,csig(1,1),'b2trcl_csigy')
        call my_out_us(70,nFc,1,sigx_c,'b2trcl_sigx_c')                    !srv 05.07.17
      endif

*   ..call 21 moment routines
#ifdef WG_TODO
#ifndef NO_21
      if (lthf21.ne.0) call b2t21t (nx, ny, ns, lnlam, bb, vol, hx,     !srv 01.07.09
     1 pbs, ne, na, ne2, te, ti, ua, fch_p, chce, chci, cthe, cthi)
      if (lvis21.ne.0) call b2t21v (nx, ny, ns, lnlam,                  !srv 01.07.09
     1 bb, vol, hx, pbs, ne, na, te, ti, ua, cvsa)
#endif
#endif
*   ..luciani
      f_luc_ke = 1.0_R8                                                 !srv 13.01.17 {
      f_luc_ki = 1.0_R8
      f_luc_et = 1.0_R8
      f_luc_sg = 1.0_R8
      f_luc_al = 1.0_R8                                                 !srv 13.01.17 }
#ifdef WG_TODO
      if(lluciani.eq.1) then   ! connection length formulation
        do iy=-1,ny
          if(conn(iy).gt.0.0_R8) then
            do ix=0,nx
              if(on_closed_surface(ix,iy)) then
                f_luc_ke(ix,iy)=1.0_R8/
     1           (1.0_R8+7.5e16_R8/conn(iy)/2.0_R8*
     2           ((te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3           (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
                f_luc_ki(ix,iy)=1.0_R8/
     1           (1.0_R8+7.5e16_R8/conn(iy)/2.0_R8*
     2           ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3           (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
                f_luc_sg(ix,iy)=f_luc_ke(ix,iy)
                f_luc_al(ix,iy)=f_luc_ke(ix,iy)
                f_luc_et(ix,iy)=f_luc_ki(ix,iy)
                do is=0,ns-1
                  if(.not.is_neutral(is)) then
                    cvsa(ix,iy,0,0,is)=
     .               cvsa(ix,iy,0,0,is)*f_luc_et(ix,iy)
                    cvsahz(ix,iy,0,0,is)=
     .               cvsahz(ix,iy,0,0,is)*f_luc_et(ix,iy)                 !srv 02.07.08
                    !WG_TODO flux limit vsaf_cl here!
                  endif
                enddo
                chce(ix,iy,0,0)=chce(ix,iy,0,0)*f_luc_ke(ix,iy)
                csig(ix,iy,0,0)=csig(ix,iy,0,0)*f_luc_sg(ix,iy)
                calf(ix,iy,0,0)=calf(ix,iy,0,0)*f_luc_al(ix,iy)
                chci(ix,iy,0,0)=chci(ix,iy,0,0)*f_luc_ki(ix,iy)
                sigx_c(ix,iy)=sigx_c(ix,iy)*f_luc_sg(ix,iy)
                alfx_c(ix,iy)=alfx_c(ix,iy)*f_luc_al(ix,iy)
              endif
            enddo
          endif
        enddo
      else if (lluciani.eq.2) then  ! SOLPS4.0 formulation
        do iy=-1,ny
          do ix=0,nx
            if (on_closed_surface(ix,iy)) then
              f_luc_ke(ix,iy)=1.0_R8/
     1             (1.0_R8+0.3e16_R8/2.0_R8*
     2           ((te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3           (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
              f_luc_ki(ix,iy)=1.0_R8/
     1             (1.0_R8+0.3e16_R8/2.0_R8*
     2           ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3           (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
              f_luc_sg(ix,iy)=f_luc_ke(ix,iy)
              f_luc_al(ix,iy)=f_luc_ke(ix,iy)
              f_luc_et(ix,iy)=f_luc_ki(ix,iy)
              do is=0,ns-1
                if(.not.is_neutral(is)) then
                  cvsa(ix,iy,0,0,is)=cvsa(ix,iy,0,0,is)*f_luc_et(ix,iy)
                  cvsahz(ix,iy,0,0,is)=
     .             cvsahz(ix,iy,0,0,is)*f_luc_et(ix,iy)
                  !WG_TODO flux limit vsaf_cl here!
                endif
              enddo
              chce(ix,iy,0,0)=chce(ix,iy,0,0)*f_luc_ke(ix,iy)
              csig(ix,iy,0,0)=csig(ix,iy,0,0)*f_luc_sg(ix,iy)
              calf(ix,iy,0,0)=calf(ix,iy,0,0)*f_luc_al(ix,iy)
              chci(ix,iy,0,0)=chci(ix,iy,0,0)*f_luc_ki(ix,iy)
            endif
          enddo
        enddo
      else if (lluciani.eq.3) then  ! new formulation
*   .. plateau and banana correction for ion heat conduction              !lk 19.09.07 {
        call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)                   !srv 01.07.09
*   .. calculate inverse aspect ratio and collisionality parameter nu_star   
        tt0 = 0.0e0_R8
        tt1 = 0.0e0_R8
        eps1 = 1.0e0_R8

!wdk init mintopcut to -2 to skip loops in case nncut = 0
!pbwdk        mintopcut = ny
        mintopcut = -2
        if (nncut > 0) mintopcut = minval(topcut(1:nncut))

        do iy = -1,mintopcut-1    ! for inner surfaces                    !srv 22.04.16
         do ix = -1,nx
          if(on_closed_surface(ix,iy)) then
           tt0(iy) = tt0(iy) + hx(ix,iy)
           tt1(iy) = tt1(iy) + vol(ix,iy)/hy1(ix,iy)                      !srv 22.04.16
          endif
         enddo
        enddo

        do iy = -1,mintopcut-1                                            !srv 22.04.16
          eps1(iy) = tt0(iy)**2/tt1(iy)
        enddo

        nu1 = 0.0e0_R8                                                    !lk 12.05.11{
        k1 = 0.58e0_R8
        nu2 = 0.0e0_R8
        k2 = 1.6e0_R8
        t1 = (4.0_R8*pi*eps0/qe)/qe
        t2 = 3.0_R8/(4.0_R8*sqrt(2.0_R8*pi))*t1*t1*ev*sqrt(ev)*sqrt(me)
     
        do iy = -1,mintopcut-1                                            !srv 22.04.16
         do ix = -1,nx
          if(on_closed_surface(ix,iy)) then
           tau(ix,iy) = t2/lnlam(ix,iy)*te(ix,iy)/ev
     &      *sqrt(te(ix,iy)/ev)/ne2(ix,iy)
           nu1(ix,iy) = conn(iy)/(tau(ix,iy)
     &      *sqrt(te(ix,iy)/me)
     &      *2.0e0_R8*pi*eps1(iy)**1.5e0_R8)

           k1(ix,iy) = 1.0e0_R8/(1.0e0_R8+2.01e0_R8*sqrt(nu1(ix,iy))
     &      +1.53e0_R8*nu1(ix,iy))
     &      +0.52e0_R8*(eps1(iy)**3)*nu1(ix,iy)
     &      /(1.0e0_R8+0.89e0_R8*sqrt(eps1(iy)**3)*nu1(ix,iy))
          
           nu2(ix,iy) = conn(iy)/(tauia(ix,iy,ismain)
     &      *sqrt(ti(ix,iy)/(mp*am(ismain)))
     &      *2.0e0_R8*pi*eps1(iy)**1.5e0_R8)
  
           k2(ix,iy) = (0.66e0_R8+1.88e0_R8*sqrt(eps1(iy))
     &      -1.54e0_R8*eps1(iy))
     &      /(1.0e0_R8+1.03e0_R8*sqrt(nu2(ix,iy))+0.31e0_R8*nu2(ix,iy))
     &      +1.17e0_R8*(eps1(iy)**3)*nu2(ix,iy)
     &      /(1.0e0_R8+0.74e0_R8*sqrt(eps1(iy)**3)*nu2(ix,iy))
          endif
         enddo
        enddo                                                             !lk 19.09.07 }
        do iy = -1,mintopcut-1                                            !srv 22.04.16
         do ix = 0,nx
          if(on_closed_surface(ix,iy)) then                               !srv 26.11.02
           f_luc_ke(ix,iy)=1.0_R8/k1(ix,iy)*0.58e0_R8*sqrt(eps1(iy)**3)
           f_luc_ki(ix,iy)=1.0_R8/k2(ix,iy)*1.60e0_R8*sqrt(eps1(iy)**3)
           f_luc_et(ix,iy)=cvsa_mltpl/
     1        (1.0e0_R8+15.12e16_R8/conn(iy)/2.0e0_R8*
     2        ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3        (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
     4        /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     5        /conn(iy)/2.0e0_R8*
     6        ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     7        (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
           f_luc_sg(ix,iy)=f_luc_ke(ix,iy)
           f_luc_al(ix,iy)=f_luc_ke(ix,iy)
           do is=0,ns-1
            if(.not.is_neutral(is)) then                                  !srv 10.07.02
             cvsa(ix,iy,0,0,is)=cvsa(ix,iy,0,0,is)*f_luc_et(ix,iy)
             cvsahz(ix,iy,0,0,is)=cvsahz(ix,iy,0,0,is)*f_luc_et(ix,iy)
             !WG_TODO flux limit vsaf_cl here!
            endif                                                         !srv 10.07.02
           enddo
! May need to use flux limiting factors calculated on cell faces
           chce(ix,iy,0,0)=chce(ix,iy,0,0)*f_luc_ke(ix,iy)
           csig(ix,iy,0,0)=csig(ix,iy,0,0)*f_luc_sg(ix,iy)
           calf(ix,iy,0,0)=calf(ix,iy,0,0)*f_luc_al(ix,iy)
           chci(ix,iy,0,0)=chci(ix,iy,0,0)*f_luc_ki(ix,iy)                    !lk 20.11.07
           sigx_c(ix,iy)=sigx_c(ix,iy)*f_luc_sg(ix,iy)
           alfx_c(ix,iy)=alfx_c(ix,iy)*f_luc_al(ix,iy)
          endif                                                           !lk 12.05.11}
         enddo
        enddo
      endif
#ifdef ASTRA
      bbaverage2 = 0.0e0_R8                                               !srv 28.07.11 {
      chci_koeff = 0.0e0_R8
      tt0 = 0.0e0_R8
      tt1 = 0.0e0_R8
      if(integral_modeling.ne.0) then
       do iy = 0,mintopcut-1    ! for inner surfaces
        inv_hi_astra(iy) =  -(Ti_astra(iy)-Ti_astra(iy-1))*1.0e3_R8/
     &   (qi_neo_astra(iy)*1.0e6_R8)/(Rmajc_astra(iy)-Rmajc_astra(iy-1))  
        do ix=0,nx                                                  
         if(on_closed_surface(ix,iy)) then
          tt0(iy) = tt0(iy) +bb(ix,iy,3)*bb(ix,iy,3)*vol(ix,iy)
          tt1(iy) = tt1(iy) +vol(ix,iy)
         endif
        enddo
        bbaverage2(iy) = tt0(iy)/tt1(iy)
        tt0(iy) = 0.0e0_R8
        do ix=0,nx                                                  
         if(on_closed_surface(ix,iy)) then
          tt0(iy) = tt0(iy) +vol(ix,iy)/bb(ix,iy,3)/bb(ix,iy,3)*
     &     (1.0_R8-bb(ix,iy,3)*bb(ix,iy,3)/bbaverage2(iy))**2 
         endif
        enddo
        do ix=0,nx                                                  
         if(on_closed_surface(ix,iy)) then
          chci_koeff(ix, iy) = 25.0_R8/4.0_R8*hy1(nmdpl,iy)/qe*tt0(iy)*   !srv 22.04.16
     &     0.5_R8*((bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)/hx(ix,iy))**2      !srv 22.04.16
     &     *vol(ix,iy)+(bb(leftix(ix,iy),iy,2)/bb(leftix(ix,iy),iy,3)
     &     /hy1(leftix(ix,iy),iy)/hx(leftix(ix,iy),iy))**2                !srv 22.04.16
     &     *vol(leftix(ix,iy),iy))*inv_hi_astra(iy) !     chi_N_astra(iy) 
         endif
        enddo
       enddo
!       do iy = 0, iy_Ratb    ! for the region of overlap with ASTRA
       do iy =0,mintopcut-1    
        do ix=0,nx                                                  
         if(on_closed_surface(ix,iy)) then
          chci(ix,iy,0,0) = chci_koeff(ix, iy)* 
     &     0.5_R8*((ti(ix,iy)*na(ix,iy,ismain))**2
     &     +(ti(leftix(ix,iy),iy)*na(leftix(ix,iy),iy,ismain))**2)
         endif
        enddo
       enddo                                   !srv 28.07.11 }
      endif
#endif
#endif
*   ..set strange velocities chve, chvi
      call sfill (nFc*2, 0.0_R8, chve, 1)
      call sfill (nFc*2, 0.0_R8, chvi, 1)

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        if(.not. is_neutral(is)) 
     &  call my_out_us(70,nFc,1,cvsahz(1,0,is),
     &              'b2trcl_luciani_cvsahzx'//chns)
!        call my_out_us(70,nFc,1,cvsahz(1,1,is),'b2trcl_luciani_cvsahzy'
!     &              //chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),'b2trcl_luciani_chcex')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_luciani_chcey')
        call my_out_us(70,nFc,1,calf(1,0),'b2trcl_luciani_calfx')
        call my_out_us(70,nFc,1,alfx_c,'b2trcl_luciani_alfx_c')              !srv 05.07.17
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_luciani_calfy')
        call my_out_us(70,nFc,1,chci(1,0),'b2trcl_luciani_chcix')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_luciani_chciy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,csig(1,0),'b2trcl_luciani_csigx')
!        call my_out_us(70,nFc,1,csig(1,1),'b2trcl_luciani_csigy')
        call my_out_us(70,nFc,1,sigx_c,'b2trcl_luciani_sigx_c')              !srv 05.07.17
      endif

*   .. save chce(:,0) after applying Luciani
      wrk = chce(:,0)                                                   !srv 13.01.1

      if(conductive_limit.ne.0) then                                      !srv 16.06.08 {
!!        AFN_TODO to be discussed what should happen here
!!        if (tn_style.ne.0) then
!!          call xerrab('flux limiters not yet developed for tn_style>0')
!!        endif
*   .. pre-compute poloidal differences of te, ti
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, dte)
       call diff_p(nCv, nFc, nVx, 0, geo, mpg, pl%ti, wrkv, dti)
*   .. apply the conductive limit of parallel electron heat flux          !srv 03.06.03 {
       if(cflim(0).ne.0.0e0_R8) then
        do iFc = 1, nFc
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c            if(.not.on_closed_surface(ix,iy)) then                         !lk 06.05.07
         if(core_cond_limit.ne.0.or.
     &      .not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
          flomx = co%chvemx(iFc)*(dv%ne(mpg%fcCv(iFc,1))+
     &                            dv%ne(mpg%fcCv(iFc,2)))/2.0_R8
          t0 = abs(chce(iFc,0)*dte(iFc))/(flomx*max(
     &      pl%te(mpg%fcCv(iFc,1)),pl%te(mpg%fcCv(iFc,2))) + cutlo)
          fllim_ke(iFc)=1.0_R8/(1.0e0_R8+t0)
          chce(iFc,0)=chce(iFc,0)*fllim_ke(iFc)
         endif
        enddo
       endif
*   .. apply the conductive limit of parallel ion heat flux
       if(cflim(1).ne.0.0e0_R8) then
        do iFc = 1,nFc
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c          if(.not.on_closed_surface(ix,iy)) then                           !lk 06.05.07
         if(core_cond_limit.ne.0.or.
     &      .not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
          flomx = co%chvimx(iFc)*(dv%ni(mpg%fcCv(iFc,1),0)
     &                           +dv%ni(mpg%fcCv(iFc,2),0))/2.0_R8
          t0 = abs(chci(iFc,0)*dti(iFc))/(flomx*max(
     &      pl%ti(mpg%fcCv(iFc,1)),pl%ti(mpg%fcCv(iFc,2))) + cutlo)
          fllim_ki(iFc) = 1.0_R8/(1.0e0_R8+t0)
          chci(iFc,0)=chci(iFc,0)*fllim_ki(iFc)
         endif
        enddo
       endif
      endif                                                               !srv 16.06.08 }
*   ..apply flux limit to calf                                            !srv 13.01.17 {
      if(style.ne.2) then
       if(cflim(0).ne.0.0e0_R8) then
        do iFc = 1, nFc
         if(.not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
            flomx = co%chvemx(iFc)*(dv%ne(mpg%fcCv(iFc,1))+
     & 	                            dv%ne(mpg%fcCv(iFc,2)))/2.0_R8
            t0 = abs(wrk(iFc)*dte(iFc))/(flomx*max(
     & 	      pl%te(mpg%fcCv(iFc,1)),pl%te(mpg%fcCv(iFc,2))) + cutlo)
            fllim_al(iFc)=1.0_R8/(1.0e0_R8+t0)
            calf(iFc,0)=calf(iFc,0)*fllim_al(iFc)            !srv 13.01.17
!WG_TODO            alfx_c(iCv)=alfx_c(iCv)*fllim_al(iCv)                !iyv 19.11.13 !srv 08.10.14
         endif
        enddo
       endif
      else
       if(cflim(3).ne.0.0e0_R8) then
*   ..apply flux limit to calf at cell faces
        do iFc = 1, nFc
         if(.not.(mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &            mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then
            jemax = sqrt(
     &       (pl%te(mpg%fcCv(iFc,1))+pl%te(mpg%fcCv(iFc,2)))/2.0_R8/me)*
     &        qe*(dv%ne(mpg%fcCv(iFc,1))+dv%ne(mpg%fcCv(iFc,2)))/2.0_R8*
     &        abs(geo%fcPbs(iFc))
            t0 = abs(calf(iFc,0)*dte(iFc))/(cflim(3)*jemax + cutlo)
            fllim_al(iFc)=1.0_R8/(1.0e0_R8+t0)
            calf(iFc,0)=calf(iFc,0)*fllim_al(iFc)                !srv 13.01.17
         endif
        enddo
*   ..apply flux limit to alfx_c at cell centers
#ifdef WG_TODO
        call gradc_p(nCv, nFc, nVx, 0, geo, mpg, pl%te, wrkv, gtec)
        do iCv = 1, nCv
!WG_TODO           if(.not.on_closed_surface(ix,iy)) then
!             if(leftix(ix,iy).eq.-2) then
!               dte =
!     &          (te(rightix(ix,iy),rightiy(ix,iy))-te(ix,iy))/
!     &          (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))/2.0_R8
!             elseif(rightix(ix,iy).eq.nx+1) then
!               dte =
!     &          (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
!     &          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
!             else
!               dte =
!     &          ((te(rightix(ix,iy),rightiy(ix,iy))-te(ix,iy))/
!     &           (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))/2.0_R8+
!     &           (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
!     &           (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8)/
!     &          2.0_R8
!             endif
             jemax = sqrt(pl%te(iCv)/me)*
     &               qe*dv%ne(iCv)*abs(geo%cvBb(iCv,0)/geo%cvBb(iCv,3))
             t0 = abs(alfx_c(iCv))*abs(gtec(iCv))/(cflim(3)*jemax + cutlo)
             fllim_al_c(iCv)=1.0_R8/(1.0e0_R8+t0)
             alfx_c(iCv)=alfx_c(iCv)*fllim_al_c(iCv)                !iyv 19.11.13 !srv 08.10.14
!WG_TODO           endif
        enddo
#endif
       endif
      endif
*   ..apply flux limit to cvsa
      cflmv = cflim(2)
      do is = 0, ns-1
       call b2tlmv (nCv, nFc, nVx, is, cflmv, geo, mpg, pl%na(1,is),
     &   pl%ti, pl%ua(1,is), co%vsaf_cl(1,0,is), cvsa(1,0,is),
     &   cvsahz(1,0,is), fllimvisc(1,is))
      enddo
*   .. save the classical transport coefficients
      csig_cl = csig
      calf_cl = calf
      cvsa_cl=cvsa
      cvsahz_cl=cvsahz                                                    !srv 09.01.01

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4) then
       do is = 0, ns-1
       write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,cvsahz(1,0,is),
     &              'b2trcl_luciani_fllim_cvsahzx'//chns)
!        call my_out_us(70,nFc,1,cvsahz(1,1,is),
!     &              'b2trcl_luciani_fllim_cvsahzy'//chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &       .or. iout_b2wdat.eq.4) then
        call my_out_us(70,nFc,1,chce(1,0),
     &             'b2trcl_luciani_fllim_chcex')
!        call my_out_us(70,nFc,1,chce(1,1),'b2trcl_luciani_fllim_chcey')
        call my_out_us(70,nFc,1,calf(1,0),
     &             'b2trcl_luciani_fllim_calfx')
!        call my_out_us(70,nFc,1,calf(1,1),'b2trcl_luciani_fllim_calfy')
        call my_out_us(70,nFc,1,alfx_c,'b2trcl_luciani_fllim_alfx_c')         !srv 05.07.17
        call my_out_us(70,nFc,1,chci(1,0),
     &             'b2trcl_luciani_fllim_chcix')
!        call my_out_us(70,nFc,1,chci(1,1),'b2trcl_luciani_fllim_chciy')
      endif
                                         
      if(iout.ne.0) then                             !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is 
        call my_out_us(70,nFc,1,cvsa(1,0,is),'cl_cvsax'//chns)
        call my_out_us(70,nFc,1,cvsa(1,1,is),'cl_cvsay'//chns)
        call my_out_us(70,nFc,1,cvsahz(1,0,is),'cl_cvsahzx'//chns)
        call my_out_us(70,nFc,1,cvsahz(1,1,is),'cl_cvsahzy'//chns)
        call my_out_us(70,nFc,1,cvsa_cl(1,0,is),'cl_cvsa_clx'//chns)
        call my_out_us(70,nFc,1,cvsa_cl(1,1,is),'cl_cvsa_cly'//chns)
        call my_out_us(70,nFc,1,cvsahz_cl(1,0,is),
     &              'cl_cvsahz_clx'//chns)
        call my_out_us(70,nFc,1,cvsahz_cl(1,1,is),
     &              'cl_cvsahz_cly'//chns)
        call my_out_us(70,nFc,1,cthe(1,is),'cl_cthe'//chns)
        call my_out_us(70,nFc,1,cthi(1,is),'cl_cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i3.3)') is0 
         call my_out_us(70,nFc,1,csigin(1,0,is,is0),'cl_csiginx'//
     &        chns//'_'//chns0)
         call my_out_us(70,nFc,1,csigin(1,1,is,is0),'cl_csiginy'//
     &        chns//'_'//chns0)
        enddo
       enddo
       call my_out_us(70,nFc,1,chce(1,0),'cl_chcex')
       call my_out_us(70,nFc,1,chce(1,1),'cl_chcey')
       call my_out_us(70,nFc,1,chci(1,0),'cl_chcix')
       call my_out_us(70,nFc,1,chci(1,1),'cl_chciy')
       call my_out_us(70,nFc,1,chve(1,0),'cl_chvex')
       call my_out_us(70,nFc,1,chve(1,1),'cl_chvey')
       call my_out_us(70,nFc,1,chvi(1,0),'cl_chvix')
       call my_out_us(70,nFc,1,chvi(1,1),'cl_chviy')
       call my_out_us(70,nFc,1,csig(1,0),'cl_csigx')
       call my_out_us(70,nFc,1,csig(1,1),'cl_csigy')
       call my_out_us(70,nFc,1,calf(1,0),'cl_calfx')
       call my_out_us(70,nFc,1,calf(1,1),'cl_calfy')
       call my_out_us(70,nFc,1,sigx_c,'cl_sigx_c')                      !srv 05.07.17
       call my_out_us(70,nFc,1,alfx_c,'cl_alfx_c')                      !srv 05.07.17
      endif                                                          !srv 17.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Still must consider how to deal with a non-orthogonal grid.

*-----------------------------------------------------------------------
*.end b2trcl

      end
