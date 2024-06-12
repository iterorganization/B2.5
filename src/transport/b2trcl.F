*-----------------------------------------------------------------------
*.specification

      subroutine b2trcl (nx, ny, ns, nscx, iscx, ismain,
     &  bb, conn, vol, gs, hx, hy, hy1, hz, qz, qc, pbs,                 !srv 23.11.10
     1  fch_p, na, ua, te, ti, po, ne, ni, ne2,                          !srv 09.01.01 01.07.09
     2  lnlam, chvemx, chvimx, cvsahz,                                   !srv 09.01.01 03.06.03
     3  cvsa, chce, chve, chci, chvi, csig, calf, csigin,                !xpb
     4  cthe, cthi, cvsa_cl, cvsahz_cl, fllimvisc,                       !srv 09.01.01
     5  csig_cl, calf_cl)                                                !xpb
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_rates                                                    !som 22.06.2021
      use b2mod_transport
     &                    , only :
     &                     f_luc_ke, f_luc_ki, f_luc_et, f_luc_sg,
     &                     f_luc_al, fllim_ki, fllim_ke, fllim_al,       !srv 13.01.17
     &                     fllim_al_c, alfx_c, sigx_c,                   !srv 08.10.14
     &                     alpha1, kt_neo                                !som 05.10.21
      use b2mod_transport_namelist
#ifdef ASTRA
      use b2mod_geo_corner
     &                    , only : nmdpl
      use b2mod_b2_to_astra
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, ismain
      integer iscx(0:nscx-1)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  hy1(-1:nx,-1:ny), hz(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne2(-1:nx,-1:ny),
     &  fch_p(-1:nx,-1:ny,0:1), conn(-1:ny), qc(-1:nx,-1:ny),            !srv 10.06.08
     &  chvemx(-1:nx,-1:ny), chvimx(-1:nx,-1:ny), lnlam(-1:nx,-1:ny)     !srv 03.06.03 01.07.09
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     &  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  csig(-1:nx,-1:ny,0:1), calf(-1:nx,-1:ny,0:1),
     &  csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),                           !xpb
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     &  cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz_cl(-1:nx,-1:ny,0:1,0:ns-1),                               !srv 09.01.01
     &  fllimvisc(-1:nx,-1:ny,0:ns-1),
     &  cvsahz(-1:nx,-1:ny,0:1,0:ns-1),                                  !srv 09.01.01
     &  csig_cl(-1:nx,-1:ny,0:1), calf_cl(-1:nx,-1:ny,0:1)
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
      integer ncall, n2, k, is, is0, lthf21, lvis21, lluciani, ix, iy
      integer, save :: conductive_limit = 1                              !srv 16.06.08
      integer, save :: style = 2                                         !srv 13.01.17
      integer, save :: Zhdanov_closure = 0                               !som 22.06.21
      integer, save :: Zhdanov_test = 0                                  !som 22.06.21
      integer, save :: Zhdanov_cond = 1                                  !som 22.06.21
      integer, save :: luc_et_form = 1                                   !som 01.09.22
      integer, save :: iout = 0, iout_b2wdat = 0,                        !srv 17.06.08
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0,
     &                 iout_b2npmo = 0, integral_modeling = 1
      character(3) chns, chns0                                           !srv 17.06.08
      integer mintopcut                                                  !lk 19.09.07
      integer mode
      real (kind=R8) ::
     *  cthiv, cthev, vis_par, jemax                                     !srv 13.01.17
      real (kind=R8) ::                                                  !srv 03.06.03
     *  t0, dte, dti, flomx
      real (kind=R8) :: tauia(-1:nx,-1:ny,0:ns-1), nu2(-1:nx,-1:ny),     !lk 19.09.07
     &  k2(-1:nx,-1:ny), tt0(-1:ny), tt1(-1:ny), eps1(-1:ny)             !lk 19.09.07
      real (kind=R8) :: tau(-1:nx,-1:ny), nu1(-1:nx,-1:ny),              !lk 12.05.11
     &  k1(-1:nx,-1:ny), t1, t2, t3, tauVeff                             !lk 12.05.11 !som 22.06.21
      real (kind=R8), save :: cvsa_mltpl = 1.0_R8, csig_mltpl=1.0_R8     !iyv 06.12.13
      real (kind=R8) ::
     *  hce(-1:nx,-1:ny,0:1), sig(-1:nx,-1:ny,0:1),
     1  alf(-1:nx,-1:ny,0:1), vsa(-1:nx,-1:ny,0:1,0:ns-1),
     2  hci(-1:nx,-1:ny,0:1), sigin(-1:nx,-1:ny,0:1,0:ns-1,0:nscx-1),
     3  cfl(-1:nx,-1:ny), wrk(-1:nx,-1:ny)                               !srv 01.07.99 09.01.01
csw 13mar2012
      integer, save :: core_cond_limit=0
*   ..procedures
      real (kind=R8) ::
     *  smin, smax                                                       !srv 22.04.16
      external xertst, xerrab, sfill, smin, smax
      external b2tqca, b2tqce, b2tqin, b2txcx, b2txcy, b2tlmv,           !xpb
     &  ipgeti, ipgetr, b2ttia                                           !lk 19.09.07
*   ..initialisation
      save ncall, lthf21, lvis21, cthiv, cthev, lluciani, mode
      data ncall/0/, lthf21/0/, lvis21/0/, mode/1/,
     1  cthiv/2.65_R8/, cthev/0.0_R8/, vis_par/1.0_R8/,
     1  lluciani/3/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trcl')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
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
     &   'faulty internal parameter cthiv')
        call ipgetr ('b2trcl_cthe', cthev)
        call xertst (0.0_R8.le.cthev,
     &   'faulty internal parameter cthev')
        call ipgetr ('b2trcl_cvsa_mltpl', cvsa_mltpl)                    !iyv 06.12.13
        call ipgetr ('b2trcl_csig_mltpl', csig_mltpl)                    !iyv 06.12.13
        call ipgeti ('b2trcl_conductive_limit', conductive_limit)        !srv 16.06.08
        call ipgeti ('b2trcl_iout', iout)                                !srv 17.06.08
        call ipgeti ('b2wdat_iout', iout_b2wdat)                         !iys 17.04.17
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                         !srv 11.09.09 }
csw 13mar2012
        call ipgeti ('b2trcl_core_cond_limit', core_cond_limit)
        call ipgeti ('b2trcl_integral_modeling', integral_modeling)      !srv 01.11.11
        call ipgeti ('b2sigp_style', style)                              !srv 13.01.17
        call ipgeti ('b2tral_Zhdanov_closure', Zhdanov_closure)          !som 22.06.21
        call ipgeti ('b2sifr_Zhdanov_test', Zhdanov_test)                !som 22.06.21
        call ipgeti ('b2tral_Zhdanov_cond', Zhdanov_cond)                !som 22.06.21
        call ipgeti ('b2trcl_luc_et_form', luc_et_form)                  !som 01.09.22
        call xertst (luc_et_form.eq.1.or.nint(zn(ismain)).eq.1,
     .   'b2trcl_luc_et_form value not appropriate for '//
     .   'non-hydrogenic plasmas!')
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
        call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test sign of vol, hx, hy, hy1
        call b2xvsg (n2, vol, 1, 'vol', '.gt.')
        call b2xvsg (n2, hx, 1, 'hx', '.gt.')
        call b2xvsg (n2, hy, 1, 'hy', '.gt.')
        call b2xvsg (n2, hy1, 1, 'hy1', '.gt.')
*    ..test range of qz
        call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &               smax(n2,qz(-1,-1,0),1).lt.1.0_R8,
     &               'faulty argument range qz(,,0)')
        call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &               smax(n2,qz(-1,-1,1),1).le.1.0_R8,
     &               'faulty argument range qz(,,1)')
*    ..test sign of na, te, ti
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2, te, 1, 'te', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

* ..compute transport coefficients
*   ..compute differential classical transport coefficients
      call b2tqca (nx, ny, ns, lnlam, na, te, ti, ne, ne2, vsa, hci)     !srv 01.07.09 20.09.11
      call b2tqce (nx, ny, lnlam, te, ti, ne, ne2, hce, sig, alf)        !srv 01.07.09 20.09.11
      call b2tqin (nx, ny, ns, nscx, iscx, bb, ti, na, sigin)            !xpb
      call sfill (n2*ns, cthev, cthe,1)
      call sfill (n2*ns, cthiv, cthi,1)
      call sfill (n2, 1.0_R8, fllim_ki, 1)                               !srv 13.01.17 {
      call sfill (n2, 1.0_R8, fllim_ke, 1)
      call sfill (n2, 1.0_R8, fllim_al, 1)
      call sfill (n2, 1.0_R8, fllim_al_c, 1)
      alfx_c(:,:) = alf(:,:,0)  !/abs((bb(:,:,0))/bb(:,:,3))             !iys 19.11.13 !srv 08.10.14 { !srv 13.01.17
      sigx_c(:,:) = sig(:,:,0)  !/abs((bb(:,:,0))/bb(:,:,3))         ! } !srv 13.01.17
*   ..incorporate geometry with hz                                       !srv 09.01.01
      wrk = hz*vol
      do is = 0, ns-1
       call b2txcx (nx, ny, mode, wrk, gs, hx, qz,                       !srv 01.07.05 {
     &   vsa(-1,-1,0,is), cvsahz(-1,-1,0,is))
       call b2txcy (nx, ny, mode, wrk, gs, hy1,
     &   vsa(-1,-1,1,is), cvsahz(-1,-1,1,is))                            !srv 01.07.05 }
      enddo                                                              !srv 09.01.01
*   ..incorporate geometry
      do is = 0, ns-1
        do is0 = 0, ns-1
          if(.not.is_neutral(is).and.is_neutral(is0)) then
            call sfill (n2*4, 0.0_R8, csigin(-1,-1,0,is,is0), 1)
          endif
        enddo
      enddo
      do is = 0, ns-1
       call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     &   vsa(-1,-1,0,is), cvsa(-1,-1,0,is))
       call b2txcy (nx, ny, mode, vol, gs, hy1,
     &   vsa(-1,-1,1,is), cvsa(-1,-1,1,is))
       do k = 0, nscx-1
         is0 = iscx(k)
         if (.not.is_neutral(is).and.nint(zn(is)).eq.1) then
           call b2txcx (nx, ny, mode, vol, gs, hx, qz,                   !srv
     &       sigin(-1,-1,0,is,k), csigin(-1,-1,0,is,is0))                !srv 01.07.99
           call b2txcy (nx, ny, mode, vol, gs, hy1,                      !srv
     &       sigin(-1,-1,1,is,k), csigin(-1,-1,1,is,is0))                !srv
         endif
       enddo
      enddo

      call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     & hce(-1,-1,0), chce(-1,-1,0))
      call b2txcy (nx, ny, mode, vol, gs, hy1,
     & hce(-1,-1,1), chce(-1,-1,1))
      call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     & hci(-1,-1,0), chci(-1,-1,0))
      call b2txcy (nx, ny, mode, vol, gs, hy1,
     & hci(-1,-1,1), chci(-1,-1,1))
      call b2txcx (nx, ny, mode, vol, gs, hx, qz,
     & sig(-1,-1,0), csig(-1,-1,0))
      call b2txcy (nx, ny, mode, vol, gs, hy1,
     & sig(-1,-1,1), csig(-1,-1,1))
      call b2txcx (nx, ny, -1, vol, gs, hx, qz,
     & alf(-1,-1,0), calf(-1,-1,0))
      call b2txcy (nx, ny, -1, vol, gs, hy1,
     & alf(-1,-1,1), calf(-1,-1,1))

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4 ) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if (.not. is_neutral(is)) then
          call my_out(70,nx,ny,cvsahz(-1,-1,0,is),
     &               'b2trcl_cvsahzx'//chns)
        endif
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                   .or. iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,chce(-1,-1,0),'b2trcl_chcex')
!        call my_out(70,nx,ny,chce(-1,-1,1),'b2trcl_chcey')
        call my_out(70,nx,ny,calf(-1,-1,0),'b2trcl_calfx')
!        call my_out(70,nx,ny,calf(-1,-1,1),'b2trcl_calfy')
        call my_out(70,nx,ny,alfx_c,'b2trcl_alfx_c')                    !srv 05.07.17
        call my_out(70,nx,ny,chci(-1,-1,0),'b2trcl_chcix')
!        call my_out(70,nx,ny,chci(-1,-1,1),'b2trcl_chciy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                   .or. iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,csig(-1,-1,0),'b2trcl_csigx')
!        call my_out(70,nx,ny,csig(-1,-1,1),'b2trcl_csigy')
        call my_out(70,nx,ny,sigx_c,'b2trcl_sigx_c')                    !srv 05.07.17
      endif

*   ..call 21 moment routines
#ifndef NO_21
      if (lthf21.ne.0) call b2t21t (nx, ny, ns, lnlam, bb, vol, hx,     !srv 01.07.09
     1 pbs, ne, na, ne2, te, ti, ua, fch_p, chce, chci, cthe, cthi)
      if (lvis21.ne.0) call b2t21v (nx, ny, ns, lnlam,                  !srv 01.07.09
     1 bb, vol, hx, pbs, ne, na, te, ti, ua, cvsa)
#endif
*   ..luciani
      f_luc_ke = 1.0_R8                                                 !srv 13.01.17 {
      f_luc_ki = 1.0_R8
      f_luc_et = 1.0_R8
      f_luc_sg = 1.0_R8
      f_luc_al = 1.0_R8                                                 !srv 13.01.17 }
      alpha1 = 1.0_R8                                                   !som 05.10.21
      kt_neo = 2.7_R8                                                   !som 05.10.21
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
                    cvsa(ix,iy,0,is)=
     .               cvsa(ix,iy,0,is)*f_luc_et(ix,iy)
                    cvsahz(ix,iy,0,is)=
     .               cvsahz(ix,iy,0,is)*f_luc_et(ix,iy)                 !srv 02.07.08
                  endif
                enddo
                chce(ix,iy,0)=chce(ix,iy,0)*f_luc_ke(ix,iy)
                csig(ix,iy,0)=csig(ix,iy,0)*f_luc_sg(ix,iy)
                calf(ix,iy,0)=calf(ix,iy,0)*f_luc_al(ix,iy)
                chci(ix,iy,0)=chci(ix,iy,0)*f_luc_ki(ix,iy)
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
                  cvsa(ix,iy,0,is)=cvsa(ix,iy,0,is)*f_luc_et(ix,iy)
                  cvsahz(ix,iy,0,is)=
     .             cvsahz(ix,iy,0,is)*f_luc_et(ix,iy)
                endif
              enddo
              chce(ix,iy,0)=chce(ix,iy,0)*f_luc_ke(ix,iy)
              csig(ix,iy,0)=csig(ix,iy,0)*f_luc_sg(ix,iy)
              calf(ix,iy,0)=calf(ix,iy,0)*f_luc_al(ix,iy)
              chci(ix,iy,0)=chci(ix,iy,0)*f_luc_ki(ix,iy)
            endif
          enddo
        enddo
      else if (lluciani.eq.3) then  ! new formulation
*   .. plateau and banana correction for ion heat conduction             !lk 19.09.07 {
        call b2ttia (nx, ny, ns, ti, ne2, lnlam, tauia)                  !srv 01.07.09
*   .. calculate inverse aspect ratio and collisionality parameter nu_star
        tt0 = 0.0e0_R8
        tt1 = 0.0e0_R8
        eps1 = 1.0e0_R8

!wdk init mintopcut to -2 to skip loops in case nncut = 0
!pbwdk        mintopcut = ny
        mintopcut = -2
        if (nncut > 0) then
          if (nnreg(0).eq.7) then
            mintopcut=topcut(1)
          else
            mintopcut=minval(topcut(1:nncut))
          endif
        endif

        do iy = -1,mintopcut-1    ! for inner surfaces                   !srv 22.04.16
         do ix = -1,nx
          if(on_closed_surface(ix,iy)) then
           tt0(iy) = tt0(iy) + hx(ix,iy)
           tt1(iy) = tt1(iy) + vol(ix,iy)/hy1(ix,iy)                     !srv 22.04.16
          endif
         enddo
        enddo

        do iy = -1,mintopcut-1                                           !srv 22.04.16
          eps1(iy) = tt0(iy)**2/tt1(iy)
        enddo

        nu1 = 0.0e0_R8                                                   !lk 12.05.11{
        k1 = 0.58e0_R8
        nu2 = 0.0e0_R8
        k2 = 1.6e0_R8
        t1 = (4.0_R8*pi*eps0/qe)/qe
        t2 = 3.0_R8/(4.0_R8*sqrt(2.0_R8*pi))*t1*t1*ev*sqrt(ev)*sqrt(me)

        do iy = -1,mintopcut-1                                           !srv 22.04.16
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

           if (Zhdanov_closure.eq.1 .and. Zhdanov_test.eq.0 .and.
     &         Zhdanov_cond.eq.1) then                                   !som 22.06.21

             tauVeff = 0.0e0_R8
             do is = 0, ns-1
               if(.not.is_neutral(is)) then
                 t3=sum(sqrt(2.0e0_R8*am(0:ns-1)/(am(is)+am(0:ns-1)))*
     &            rz2(ix,iy,0:ns-1)*na(ix,iy,0:ns-1))
                 tauVeff=tauVeff+(na(ix,iy,is)/ni(ix,iy,1))*
     &            tauia(ix,iy,is)*(ne2(ix,iy)/t3)*sqrt(ti(ix,iy)/
     &           (mp*am(is)))
               endif
             enddo

             nu2(ix,iy) = conn(iy)/(tauVeff
     &        *2.0e0_R8*pi*eps1(iy)**1.5e0_R8)
           else

             nu2(ix,iy) = conn(iy)/(tauia(ix,iy,ismain)
     &        *sqrt(ti(ix,iy)/(mp*am(ismain)))
     &        *2.0e0_R8*pi*eps1(iy)**1.5e0_R8)
           endif

           k2(ix,iy) = (0.66e0_R8+1.88e0_R8*sqrt(eps1(iy))
     &      -1.54e0_R8*eps1(iy))
     &      /(1.0e0_R8+1.03e0_R8*sqrt(nu2(ix,iy))+0.31e0_R8*nu2(ix,iy))
     &      +1.17e0_R8*(eps1(iy)**3)*nu2(ix,iy)
     &      /(1.0e0_R8+0.74e0_R8*sqrt(eps1(iy)**3)*nu2(ix,iy))
           alpha1(ix,iy) = ((-0.17e0_R8+1.05e0_R8*sqrt(nu2(ix,iy))+
     &           2.7e0_R8*nu2(ix,iy)**2*eps1(iy)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu2(ix,iy))+
     &           nu2(ix,iy)**2*eps1(iy)**3)-1.0e0_R8)/1.7e0_R8           !som 05.10.21
           kt_neo(ix,iy) =
     &           ((-0.17e0_R8+1.05e0_R8*sqrt(nu2(ix,iy))
     &           +2.7e0_R8*(nu2(ix,iy))**2*eps1(iy)**3)/(1.0e0_R8+
     &           0.7e0_R8*sqrt(nu2(ix,iy))+
     &           (nu2(ix,iy))**2*eps1(iy)**3))                           !som 05.10.21
          endif
         enddo
        enddo                                                            !lk 19.09.07 }
        do iy = -1,mintopcut-1                                           !srv 22.04.16
         do ix = 0,nx
          if(on_closed_surface(ix,iy)) then                              !srv 26.11.02
           f_luc_ke(ix,iy)=1.0_R8/k1(ix,iy)*0.58e0_R8*sqrt(eps1(iy)**3)
           f_luc_ki(ix,iy)=1.0_R8/k2(ix,iy)*1.60e0_R8*sqrt(eps1(iy)**3)
           if(luc_et_form.eq.0) then ! for hydrogen isotopes + imp plasmas only
             f_luc_et(ix,iy)=cvsa_mltpl/
     1        (1.0e0_R8+15.12e16_R8/conn(iy)/2.0e0_R8*
     2        ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     3        (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
     4        /(1.0e0_R8+15.12e16_R8*eps1(iy)**1.5e0_R8
     5        /conn(iy)/2.0e0_R8*
     6        ((ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/ev)**2/
     7        (ne2(ix,iy)+ne2(leftix(ix,iy),leftiy(ix,iy))))
           else if (luc_et_form.eq.1) then ! generalized formulation     !som 01.09.22
             f_luc_et(ix,iy)=cvsa_mltpl/
     1        (1.0e0_R8+2.0e0_R8/eps1(iy)**1.5e0_R8/
     2        (nu2(ix,iy)+nu2(leftix(ix,iy),leftiy(ix,iy))))
     3        /(1.0e0_R8+2.0e0_R8/
     4        (nu2(ix,iy)+nu2(leftix(ix,iy),leftiy(ix,iy))))
           endif
           f_luc_sg(ix,iy)=f_luc_ke(ix,iy)
           f_luc_al(ix,iy)=f_luc_ke(ix,iy)
           do is=0,ns-1
            if(.not.is_neutral(is)) then                                 !srv 10.07.02
             cvsa(ix,iy,0,is)=cvsa(ix,iy,0,is)*f_luc_et(ix,iy)
             cvsahz(ix,iy,0,is)=cvsahz(ix,iy,0,is)*f_luc_et(ix,iy)
            endif                                                        !srv 10.07.02
           enddo
! May need to use flux-limiting factors calculated on cell faces
           chce(ix,iy,0)=chce(ix,iy,0)*f_luc_ke(ix,iy)
           csig(ix,iy,0)=csig(ix,iy,0)*f_luc_sg(ix,iy)
           calf(ix,iy,0)=calf(ix,iy,0)*f_luc_al(ix,iy)
           chci(ix,iy,0)=chci(ix,iy,0)*f_luc_ki(ix,iy)                   !lk 20.11.07
           sigx_c(ix,iy)=sigx_c(ix,iy)*f_luc_sg(ix,iy)
           alfx_c(ix,iy)=alfx_c(ix,iy)*f_luc_al(ix,iy)
           if( bottomiy(ix,iy).eq.-2 ) then
            csig(ix,iy,0) = csig_mltpl*csig(ix,iy,0)                     !iyv
            sigx_c(ix,iy) = csig_mltpl*sigx_c(ix,iy)
           endif
          endif                                                          !lk 12.05.11}
         enddo
        enddo
      endif
#ifdef ASTRA
      bbaverage2 = 0.0e0_R8                                              !srv 28.07.11 {
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
          chci_koeff(ix,iy) = 25.0_R8/4.0_R8*hy1(nmdpl,iy)/qe*tt0(iy)*   !srv 22.04.16
     &     0.5_R8*((bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)/hx(ix,iy))**2     !srv 22.04.16
     &     *vol(ix,iy)+(bb(leftix(ix,iy),iy,2)/bb(leftix(ix,iy),iy,3)
     &     /hy1(leftix(ix,iy),iy)/hx(leftix(ix,iy),iy))**2               !srv 22.04.16
     &     *vol(leftix(ix,iy),iy))*inv_hi_astra(iy) !     chi_N_astra(iy)
         endif
        enddo
       enddo
!       do iy = 0, iy_Ratb    ! for the region of overlap with ASTRA
       do iy=0,mintopcut-1
        do ix=0,nx
         if(on_closed_surface(ix,iy)) then
          chci(ix,iy,0) = chci_koeff(ix, iy)*
     &     0.5_R8*((ti(ix,iy)*na(ix,iy,ismain))**2
     &     +(ti(leftix(ix,iy),iy)*na(leftix(ix,iy),iy,ismain))**2)
         endif
        enddo
       enddo                                                             !srv 28.07.11 }
      endif
#endif
*   ..set strange velocities chve, chvi
      call sfill (n2*2, 0.0_R8, chve, 1)
      call sfill (n2*2, 0.0_R8, chvi, 1)

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        if(.not. is_neutral(is))
     &  call my_out(70,nx,ny,cvsahz(-1,-1,0,is),
     &              'b2trcl_luciani_cvsahzx'//chns)
!        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),
!     &              'b2trcl_luciani_cvsahzy'//chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,chce(-1,-1,0),'b2trcl_luciani_chcex')
!        call my_out(70,nx,ny,chce(-1,-1,1),'b2trcl_luciani_chcey')
        call my_out(70,nx,ny,calf(-1,-1,0),'b2trcl_luciani_calfx')
        call my_out(70,nx,ny,alfx_c,'b2trcl_luciani_alfx_c')             !srv 05.07.17
!        call my_out(70,nx,ny,calf(-1,-1,1),'b2trcl_luciani_calfy')
        call my_out(70,nx,ny,chci(-1,-1,0),'b2trcl_luciani_chcix')
!        call my_out(70,nx,ny,chci(-1,-1,1),'b2trcl_luciani_chciy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0
     &                  .or. iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,csig(-1,-1,0),'b2trcl_luciani_csigx')
!        call my_out(70,nx,ny,csig(-1,-1,1),'b2trcl_luciani_csigy')
        call my_out(70,nx,ny,sigx_c,'b2trcl_luciani_sigx_c')             !srv 05.07.17
      endif

*   .. save chce(:,:,0) after applying Luciani
      wrk = chce(:,:,0)                                                  !srv 13.01.17

      if(conductive_limit.ne.0) then                                     !srv 16.06.08 {
*   .. apply the conductive limit of parallel electron heat flux         !srv 03.06.03 {
       if(cflim(0).ne.0.0e0_R8.or.cflme_core.ne.cflme_sol) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c            if(.not.on_closed_surface(ix,iy)) then                      !lk 06.05.07
           if(core_cond_limit.ne.0.or.
     &        .not.on_closed_surface(ix,iy)) then
            dte = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
            flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &       ne(ix,iy))/2.0_R8
            t0 = chce(ix,iy,0)*abs(dte)/(flomx*max(
     &       te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
            fllim_ke(ix,iy)=1.0_R8/(1.0e0_R8+t0)
            chce(ix,iy,0)=chce(ix,iy,0)*fllim_ke(ix,iy)
           endif
          endif
         enddo
        enddo
       endif
*   .. apply the conductive limit of parallel ion heat flux
       if(cflim(1).ne.0.0e0_R8.or.cflmi_core.ne.cflmi_sol) then
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).ne.-2) then
csw 19oct2011, also in core for SOLPS4.3 comparison!!!
c          if(.not.on_closed_surface(ix,iy)) then                        !lk 06.05.07
           if(core_cond_limit.ne.0.or.
     &        .not.on_closed_surface(ix,iy)) then
            dti = ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy))
            flomx = chvimx(ix,iy)*(ni(leftix(ix,iy),leftiy(ix,iy),0)
     &                          +ni(ix,iy,0))/2.0_R8
            t0 = chci(ix,iy,0)*abs(dti)/
     &        (flomx*max(ti(leftix(ix,iy),leftiy(ix,iy)),ti(ix,iy)))
            fllim_ki(ix,iy) = 1.0_R8/(1.0e0_R8+t0)
            chci(ix,iy,0)=chci(ix,iy,0)*fllim_ki(ix,iy)
           endif
          endif
         enddo
        enddo
       endif
      endif                                                              !srv 16.06.08 }
*   ..apply flux limit to calf                                           !srv 13.01.17 {
      if(style.eq.1) then
       if(cflim(0).ne.0.0e0_R8.or.cflme_core.ne.cflme_sol) then
        do iy = -1, ny
         do ix = -1, nx
           if(leftix(ix,iy).ne.-2) then
             if(.not.on_closed_surface(ix,iy)) then
               dte = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
               flomx = chvemx(ix,iy)*(ne(leftix(ix,iy),leftiy(ix,iy))+
     &                                ne(ix,iy))/2.0_R8
               t0 = wrk(ix,iy)*abs(dte)/(flomx*max(
     &              te(leftix(ix,iy),leftiy(ix,iy)),te(ix,iy)))
               fllim_al(ix,iy)=1.0_R8/(1.0e0_R8+t0)
               calf(ix,iy,0)=calf(ix,iy,0)*fllim_al(ix,iy)               !srv 13.01.17
               alfx_c(ix,iy)=alfx_c(ix,iy)*fllim_al(ix,iy)               !iyv 19.11.13 !srv 08.10.14
             endif
           endif
         enddo
        enddo
       endif
      else if (style.eq.2) then
       if(cflim(3).ne.0.0e0_R8) then
*   ..apply flux limit to calf at cell faces
        do iy = -1, ny
         do ix = -1, nx
           if(.not.on_closed_surface(ix,iy)) then
             if(leftix(ix,iy).ne.-2) then
               dte = te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy))
               jemax = sqrt(
     &          (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8/me)*
     &           qe*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))/2.0_R8*
     &           abs(pbs(ix,iy,0))
               t0 = abs(calf(ix,iy,0)*dte)/(cflim(3)*jemax)
               fllim_al(ix,iy)=1.0_R8/(1.0e0_R8+t0)
               calf(ix,iy,0)=calf(ix,iy,0)*fllim_al(ix,iy)               !srv 13.01.17
             endif
           endif
         enddo
        enddo
*   ..apply flux limit to alfx_c at cell centers
        do iy = -1, ny
         do ix = -1, nx
           if(.not.on_closed_surface(ix,iy)) then
             if(leftix(ix,iy).eq.-2) then
               dte =
     &          (te(rightix(ix,iy),rightiy(ix,iy))-te(ix,iy))/
     &          (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))/2.0_R8
             elseif(rightix(ix,iy).eq.nx+1) then
               dte =
     &          (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
     &          (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8
             else
               dte =
     &          ((te(rightix(ix,iy),rightiy(ix,iy))-te(ix,iy))/
     &           (hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy))/2.0_R8+
     &           (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))/
     &           (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8)/
     &          2.0_R8
             endif
             jemax = sqrt(te(ix,iy)/me)*
     &               qe*ne(ix,iy)*abs(bb(ix,iy,0)/bb(ix,iy,3))
             t0 = abs(alfx_c(ix,iy))*abs(dte)/(cflim(3)*jemax)
             fllim_al_c(ix,iy)=1.0_R8/(1.0e0_R8+t0)
             alfx_c(ix,iy)=alfx_c(ix,iy)*fllim_al_c(ix,iy)               !iyv 19.11.13 !srv 08.10.14
           endif
         enddo
        enddo
       endif
      endif
*   ..apply flux limit to cvsa
      call cfl_tanh_profile(nx, ny, cflim(2),
     .  cflmv_core, cflmv_sol, cflmv_tanh_a, cflmv_tanh_b, cfl)
      do is = 0, ns-1
       call b2tlmv (nx, ny, is, cfl, pbs(-1,-1,0), qc, na(-1,-1,is),     !srv 12.04.08
     &   ti, ua(-1,-1,is), cvsa(-1,-1,0,is), cvsahz(-1,-1,0,is),
     &   fllimvisc(-1,-1,is))
      enddo
*   ..save the classical transport coefficients
      csig_cl = csig
      calf_cl = calf
      cvsa_cl = cvsa
      cvsahz_cl = cvsahz                                                 !srv 09.01.01

      if(iout_b2npmo.ne.0 .or. iout_b2wdat.eq.4) then
       do is = 0, ns-1
       write (chns,'(i3.3)') is
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),
     &              'b2trcl_luciani_fllim_cvsahzx'//chns)
!        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),
!     &              'b2trcl_luciani_fllim_cvsahzy'//chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &       .or. iout_b2wdat.eq.4) then
        call my_out(70,nx,ny,chce(-1,-1,0),
     &             'b2trcl_luciani_fllim_chcex')
!        call my_out(70,nx,ny,chce(-1,-1,1),'b2trcl_luciani_fllim_chcey')
        call my_out(70,nx,ny,calf(-1,-1,0),
     &             'b2trcl_luciani_fllim_calfx')
!        call my_out(70,nx,ny,calf(-1,-1,1),'b2trcl_luciani_fllim_calfy')
        call my_out(70,nx,ny,alfx_c,'b2trcl_luciani_fllim_alfx_c')       !srv 05.07.17
        call my_out(70,nx,ny,chci(-1,-1,0),
     &             'b2trcl_luciani_fllim_chcix')
!        call my_out(70,nx,ny,chci(-1,-1,1),'b2trcl_luciani_fllim_chciy')
      endif

      if(iout.ne.0) then                                                 !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'cl_cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'cl_cvsay'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'cl_cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'cl_cvsahzy'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'cl_cvsa_clx'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,1,is),'cl_cvsa_cly'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,0,is),
     &              'cl_cvsahz_clx'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,1,is),
     &              'cl_cvsahz_cly'//chns)
        call my_out(70,nx,ny,cthe(-1,-1,is),'cl_cthe'//chns)
        call my_out(70,nx,ny,cthi(-1,-1,is),'cl_cthi'//chns)
        do is0 = 0, ns-1
         write (chns0,'(i3.3)') is0
         call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'cl_csiginx'//
     &        chns//'_'//chns0)
         call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'cl_csiginy'//
     &        chns//'_'//chns0)
        enddo
       enddo
       call my_out(70,nx,ny,chce(-1,-1,0),'cl_chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'cl_chcey')
       call my_out(70,nx,ny,chci(-1,-1,0),'cl_chcix')
       call my_out(70,nx,ny,chci(-1,-1,1),'cl_chciy')
       call my_out(70,nx,ny,chve(-1,-1,0),'cl_chvex')
       call my_out(70,nx,ny,chve(-1,-1,1),'cl_chvey')
       call my_out(70,nx,ny,chvi(-1,-1,0),'cl_chvix')
       call my_out(70,nx,ny,chvi(-1,-1,1),'cl_chviy')
       call my_out(70,nx,ny,csig(-1,-1,0),'cl_csigx')
       call my_out(70,nx,ny,csig(-1,-1,1),'cl_csigy')
       call my_out(70,nx,ny,calf(-1,-1,0),'cl_calfx')
       call my_out(70,nx,ny,calf(-1,-1,1),'cl_calfy')
       call my_out(70,nx,ny,sigx_c,'cl_sigx_c')                          !srv 05.07.17
       call my_out(70,nx,ny,alfx_c,'cl_alfx_c')                          !srv 05.07.17
      endif                                                              !srv 17.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Still must consider how to deal with a non-orthogonal grid.

*-----------------------------------------------------------------------
*.end b2trcl

      end subroutine b2trcl
