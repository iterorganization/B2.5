*-----------------------------------------------------------------------
*.specification

      subroutine b2stbc (nx, ny, ns, ismain, ismain0,
     &  BoRiS, facdrift, fac_ExB,                                       !xpb
     &  na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     &  te, ti, po, ne, ni, floe_noc, floi_noc,                         !srv 28.10.02
     1  fna, fna_mdf, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,
     &  fhe, fhi, fhe_mdf, fhi_mdf, fhm, fch, fne, fni, fna_fcor,       !srv 16.02.12
     &  fna_nodrift, fna_he, fna_eir, fnaPSch, fhePSch, fhiPSch,        !srv 25.02.00 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,  !srv 21.10.02
     &  fchanml,
     2  cdna, chce, chve, chci, chvi, calf,
     &  sna0, smo0, she0, shi0, sch0, sne0, wrong_flow, main_call)      !xpb
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_indirect
      use b2mod_boundary_namelist
      use b2mod_geo_corner
      use b2mod_constants
      use b2mod_sources
     1 , only : b2stbc_sna,b2stbc_smo,b2stbc_she,b2stbc_shi,b2stbc_sch,
     2          b2stbc_sne
      use b2mod_balance !djm Jan 2017
     1 , only : b2stbc_sna0to1,b2stbc_smo0to3,b2stbc_she0to3,
     2          b2stbc_shi0to3,balance_netcdf
      use b2mod_feedback
      use b2mod_b2cmpa
      use b2mod_external
      use b2mod_time
      use b2mod_diag
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0
      real (kind=R8) ::
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1),
     2  vedia(-1:nx,-1:ny,0:1), vadia(-1:nx,-1:ny,0:1,0:ns-1),
     3  veecrb(-1:nx,-1:ny,0:1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fna_32(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_52(-1:nx,-1:ny,0:1,0:ns-1), fna_mdf(-1:nx,-1:ny,0:1,0:ns-1),
     &  fni_32(-1:nx,-1:ny,0:1), fni_52(-1:nx,-1:ny,0:1),
     &  fne_32(-1:nx,-1:ny,0:1), fne_52(-1:nx,-1:ny,0:1),
     &  fhe(-1:nx,-1:ny,0:1), fhe_mdf(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fhi_mdf(-1:nx,-1:ny,0:1),
     &  fch(-1:nx,-1:ny,0:1), fhm(-1:nx,-1:ny,0:1,0:ns-1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),              !srv 11.10.02 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &  fchinert(-1:nx,-1:ny,0:1),                                       !srv 21.10.02 }
     &  fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),                                !srv 25.02.00
     &  fna_nodrift(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_he(-1:nx,-1:ny,0:1,0:ns-1),
     &  fna_eir(-1:nx,-1:ny,0:1,0:ns-1),
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1), fchanml(-1:nx,-1:ny,0:1),       !srv 13.10.06
     &  fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),              !srv 13.10.06
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1)             !srv 28.10.02
      real (kind=R8) ::
     1  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     2  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     3  calf(-1:nx,-1:ny,0:1), cdna(-1:nx,-1:ny,0:1,0:ns-1)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, regcnt, ixx
      real (kind=R8) ::
     *   sna0ep, she0ep, shi0ep, t0, testqe, testqi,
     *   b2stbc_neoclassical, b2stbc_fchy_dia,
     &   praverage(2), prsum(2), gssum(2), tesum(2), tisum(2),
     1   nasum(0:ns-1,2),
     1   nef, nif, tef, tif, naf(0:ns-1), zaf(0:ns-1),
     x   nxf, taf, zxf,
     1   volsum(2), qeaverage(2), qiaverage(2), qesum(2), qisum(2),
     1   sum_sch0, sum_vol, integral_current
      integer diagno, b2stbc_fchy_dia_coreonly,
     *   area_fix, sheath_drift_fix
      integer element_done(-1:nx)
      integer b2stbc_feedback, secmodel
      save ncall, sna0ep, she0ep, shi0ep, diagno,
     *   b2stbc_neoclassical, b2stbc_fchy_dia, integral_current,
     1   b2stbc_feedback, b2stbc_fchy_dia_coreonly, secmodel,
     *   area_fix, sheath_drift_fix
*   ..procedures
      intrinsic nint
      real (kind=R8) ::
     *   b2ssum
      external xertst, sfill, b2ssum, b2saxpy, ipgetr, ipgeti
      external b2xvsg, b2xvff, b2xvps, b2xxgs,
     1   b2stbc_phys, b2stbc_bas, b2stbc_fb_pre, b2stbc_fb_post
*   ..initialisation
      data ncall/0/, sna0ep/1.0e-36_R8/, she0ep/1.0e-36_R8/,
     1   shi0ep/1.0e-36_R8/, diagno/0/,
     2   b2stbc_neoclassical/0.0_R8/, b2stbc_fchy_dia/0.0_R8/,
     3   integral_current/0.0_R8/, b2stbc_feedback/0/, area_fix/3/,
     4   b2stbc_fchy_dia_coreonly/1/, secmodel/0/,
     5   sheath_drift_fix/1/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_sna0ep', sna0ep)
       call xertst (0.0_R8.le.sna0ep,
     .  'faulty internal parameter sna0ep')
       call ipgetr ('b2stbc_she0ep', she0ep)
       call xertst (0.0_R8.le.she0ep,
     .  'faulty internal parameter she0ep')
       call ipgetr ('b2stbc_shi0ep', shi0ep)
       call xertst (0.0_R8.le.shi0ep,
     .  'faulty internal parameter shi0ep')
       call get_boundary_namelist
       call ipgeti ('b2stbc_diagno',diagno)                              !xpb
       call ipgetr ('b2stbc_neoclassical', b2stbc_neoclassical)          !xpb
       call xertst (0.0_R8.le.b2stbc_neoclassical .and.                  !xpb
     &              b2stbc_neoclassical.le.1.0_R8,                       !xpb
     1              'faulty internal parameter b2stbc_neoclassical')     !xpb
       call ipgetr ('b2stbc_fchy_dia', b2stbc_fchy_dia)                  !xpb
       call xertst (0.0_R8.le.b2stbc_fchy_dia .and.                      !xpb
     &              b2stbc_fchy_dia.le.1.0_R8,                           !xpb
     1              'faulty internal parameter b2stbc_fchy_dia')         !xpb
       call ipgetr ('b2stbc_integral_current', integral_current)         !dpc
       call xertst (.not.(integral_current.ne.0.0_R8.and.
     .  ((maxval(facdrift).gt.0.0_R8.and.b2stbc_fchy_dia.gt.0.0_R8).or.
     .   b2stbc_neoclassical.gt.0.0_R8)),
     .  'inconsistent boundary settings for current equation')
       call ipgeti ('b2stbc_feedback', b2stbc_feedback)
       lfeedback=b2stbc_feedback.ne.0                          ! dpc-2003-10-23
       call ipgeti ('b2stbc_fchy_dia_coreonly',
     .               b2stbc_fchy_dia_coreonly)
       call ipgeti ('b2stbc_secmodel', secmodel)
       call ipgeti ('b2news_area_fix', area_fix)
       call ipgeti ('b2stbc_sheath_drift_fix', sheath_drift_fix)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main neutral species ismain; must not be neutral')
      call xertst (0.le.ismain0.and.ismain0.lt.ns,
     . 'invalid main neutral species index ismain0')
      call xertst (is_neutral(ismain0).or.ismain.eq.ismain0,
     . 'invalid main neutral species ismain0; must be neutral')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvps (nx, ny, ns,
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*    ..test edge values of fne, fni, chce, chve, chci, chvi, calf
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
       call b2xvff (nx, ny, calf, 'calf')
*    ..test sign of chce, chci
       call b2xvsg (n2*2, chce, 1, 'chce', '.ge.')
       call b2xvsg (n2*2, chci, 1, 'chci', '.ge.')
      endif
*   ..test facdrift           !xpb
      call xertst (0.0_R8.le.minval(facdrift) .and.
     & 1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..check kinetic secondary electron emission model
      if (ncall.eq.0 .and. secmodel.eq.1) then
        call xertst (nint(zn(ismain)).eq.1,
     &   'Kinetic secondary emission model only for hydrogen plasmas!')
        write (*,*) 'Kinetic secondary emission model turned on'
        write (*,*) 'Graphite target plates are assumed'
        if (nint(zamax(ismain)).ne.1) write(*,*)
     .   'WARNING! Kinetic secondary emission expects HYDROGEN plasmas!'
      endif
*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0_R8, smo0, 1)
      call sfill (n2*2, 0.0_R8, sne0, 1)
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)

* ..compute standard form volume sources
*   (This code placed here for want of a better location)
*   ..low-level particle source
      do is = 0, ns-1
       call b2saxpy (n2, sna0ep, vol, 1, sna0(-1,-1,0,is), 1)
      enddo
      call b2saxpy (n2, she0ep, vol, 1, she0(-1,-1,0), 1)
      call b2saxpy (n2, shi0ep, vol, 1, shi0(-1,-1,0), 1)

*   ..find the average pressure for inner flux surface           !sv 24.03.99
      iy = -1                                                    !sv 24.03.99
      tesum(1)=0.0_R8
      tisum(1)=0.0_R8
      nasum(:,1)=0.0_R8
      prsum(1)=0.0_R8                                                    !xpb
      qesum(1)=0.0_R8
      qisum(1)=0.0_R8
      gssum(1)=0.0_R8                                                    !xpb
      volsum(1)=0.0_R8
      do ix = -1, nx                                                     !xpb
        if(on_closed_surface(ix,iy)) then                                !xpb
          tef = (te(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           te(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ti(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ne(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          nif = (ni(ix,iy,1)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ni(topix(ix,iy),topiy(ix,iy),1)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           na(topix(ix,iy),topiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          zaf(:) = (rza(ix,iy,:)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           rza(topix(ix,iy),topiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tesum(1)=tesum(1) + tef*gs(topix(ix,iy),topiy(ix,iy),1)
          tisum(1)=tisum(1) + tif*gs(topix(ix,iy),topiy(ix,iy),1)
          nasum(:,1)=nasum(:,1) + naf(:)*gs(topix(ix,iy),topiy(ix,iy),1)
          prsum(1)=prsum(1) +
     1     (nef*tef+nif*tif)*vol(ix,iy)
     2     *wbbv(topix(ix,iy),topiy(ix,iy),2)*vol(ix,iy)/gs(ix,iy,2)
          qesum(1)=qesum(1)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     nef*tef**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/qe*
     2     gs(topix(ix,iy),topiy(ix,iy),1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(1)=qisum(1)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         naf(is)*tif**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/
     5         (qe*zaf(is))*gs(topix(ix,iy),topiy(ix,iy),1)
            endif
          enddo
          do is = 0, ns_ext-1
            nxf = (na_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &             na_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &            (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            taf = (ta_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &             ta_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &            (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
            prsum(1)=prsum(1) + nxf*taf*vol(ix,iy)
     2        *wbbv(topix(ix,iy),topiy(ix,iy),2)*vol(ix,iy)/gs(ix,iy,2)
            if (.not.is_ext_neutral(is)) then
              zxf = (za_ext(ix,iy,is)*hy1(topix(ix,iy),topiy(ix,iy))+
     &               za_ext(topix(ix,iy),topiy(ix,iy),is)*hy1(ix,iy))/
     &              (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
              qisum(1)=qisum(1)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         nxf*taf**2*wbbv(topix(ix,iy),topiy(ix,iy),2)/
     5         (qe*zxf)*gs(topix(ix,iy),topiy(ix,iy),1)
            endif
          enddo
          gssum(1)=gssum(1) + gs(topix(ix,iy),topiy(ix,iy),1)            !xpb
          volsum(1)=volsum(1) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                      !sv 24.03.99
      if (volsum(1).gt.0.0_R8) praverage(1)=prsum(1)/volsum(1)           !xpb
      if (gssum(1).gt.0.0_R8) then
        qeaverage(1)=qesum(1)/gssum(1)
        qiaverage(1)=qisum(1)/gssum(1)
      end if
*   ..find the average pressure for island center flux surface           !xpb
      iy = ny                                                            !xpb
      tesum(2)=0.0_R8
      tisum(2)=0.0_R8
      nasum(:,2)=0.0_R8
      prsum(2)=0.0_R8                                                    !xpb
      qesum(2)=0.0_R8
      qisum(2)=0.0_R8
      gssum(2)=0.0_R8                                                    !xpb
      volsum(2)=0.0_R8
      praverage(2)=0.0_R8                                                !xpb
      qeaverage(2)=0.0_R8                                                !xpb
      qiaverage(2)=0.0_R8                                                !xpb
      do ix = -1, nx                                                     !xpb
        if(on_closed_surface(ix,iy)) then                                !xpb
          tef = (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ne(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          nif = (ni(ix,iy,1)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ni(bottomix(ix,iy),bottomiy(ix,iy),1)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           na(bottomix(ix,iy),bottomiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          zaf(:) = (rza(ix,iy,:)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           rza(bottomix(ix,iy),bottomiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tesum(2)=tesum(2) + tef*gs(ix,iy,1)                            !xpb
          tisum(2)=tisum(2) + tif*gs(ix,iy,1)                            !xpb
          nasum(:,2)=nasum(:,2) + naf(:)*gs(ix,iy,1)                     !xpb
          prsum(2)=prsum(2) +
     1     (nef*tef+nif*tif)*vol(ix,iy)*
     2     wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)*
     3     vol(ix,iy)/gs(ix,iy,2)
          qesum(2)=qesum(2)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     nef*tef**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/qe*
     2     gs(ix,iy,1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(2)=qisum(2)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         naf(is)*tif**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/
     5         (qe*zaf(is))*gs(ix,iy,1)
            endif
          enddo
          do is = 0, ns_ext-1
            nxf = (na_ext(ix,iy,is)*
     .             hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &             na_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .             hy1(ix,iy))/
     &            (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            taf = (ta_ext(ix,iy,is)*
     .             hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &             ta_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .             hy1(ix,iy))/
     &            (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            prsum(2)=prsum(2) + nxf*taf*vol(ix,iy)
     2        *wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)
     3        *vol(ix,iy)/gs(ix,iy,2)
            if (.not.is_ext_neutral(is)) then
              zxf = (za_ext(ix,iy,is)*
     .               hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &               za_ext(bottomix(ix,iy),bottomiy(ix,iy),is)*
     .               hy1(ix,iy))/
     &              (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
              qisum(2)=qisum(2)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         nxf*taf**2*wbbv(bottomix(ix,iy),bottomiy(ix,iy),2)/
     5         (qe*zxf)*gs(ix,iy,1)
            endif
          enddo
          gssum(2)=gssum(2) + gs(ix,iy,1)                                !xpb
          volsum(2)=volsum(2) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                              !xpb
      if (volsum(2).gt.0.0_R8) praverage(2)=prsum(2)/volsum(2)           !xpb
      if (gssum(2).gt.0.0_R8) then
        qeaverage(2)=qesum(2)/gssum(2)                                   !xpb
        qiaverage(2)=qisum(2)/gssum(2)                                   !xpb
      endif
      if (ncall.eq.0) lbndusr=.false.                                    !xpb

      if(boundary_namelist.ge.1) then

         if(lfeedback.or.b2stbc_feedback.ne.0)
     1     call b2stbc_fb_pre(nx, ny, ns, ismain, BoRiS,
     1      vol, gs, na, ne, ni, te, ti,
     &      fna, fna_mdf, fna_fcor, fna_nodrift,
     &      fna_he, fnaPSch, fna_eir, fna_32, fna_52,
     &      fhe, fhi, fhm, fch,
     &      sna0, smo0, she0, shi0, sch0, sne0,
     &      main_call, diagno)

        if(ncall.eq.0) write(*,*) ' b2stbc_phys called'
        call b2stbc_phys(nx, ny, ns, ismain, ismain0,
     1   BoRiS, facdrift, fac_ExB,                                       !srv 24.08.08
     1   na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     2   te, ti, po, ne, ni, floe_noc, floi_noc, fna, fna_mdf,           !srv 28.10.02
     1   fnaPSch, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,        !srv 13.10.06
     &   fhe, fhi, fhm, fhe_mdf, fhi_mdf, fna_fcor, fhePSch, fhiPSch,    !srv 16.02.12
     &   fch, fne, fni,                                                  !srv 13.10.06
     &   fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,  !srv 21.10.02
     &   fchanml,
     &   cdna, chce, chve, chci, chvi, calf,                             !srv 25.09.09
     &   sna0, smo0, she0, shi0, sch0, sne0,
     &   wrong_flow, main_call, diagno,
     &   b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &   b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     &   sheath_drift_fix)
        if(wrong_flow) write(*,*)
     &        'b2stbc: wrong_flow returned from b2stbc_phys'
        if(lbndusr) then
          if(ncall.eq.0) write(*,*) 'BNDUSR called'
          call b2stbc_bas(nx, ny, ns, ismain,
     1     facdrift, na, ua, uadia, vaecrb, te, ti, po, ne, ni,
     1     fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &     sna0, smo0, she0, shi0, sch0, sne0,
     &     wrong_flow, main_call, diagno,
     &     b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &     b2stbc_fchy_dia_coreonly, secmodel, lbndusr,
     &     area_fix, sheath_drift_fix)
          if(wrong_flow) write(*,*)
     &           'b2stbc: wrong_flow returned from b2stbc_bas'
        endif      !lbndusr
      else      !boundary_namelist.ge.1

        call b2stbc_bas(nx, ny, ns, ismain,
     1   facdrift, na, ua, uadia, vaecrb, te, ti, po, ne, ni,
     1   fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &   sna0, smo0, she0, shi0, sch0, sne0,
     &   wrong_flow, main_call, diagno,
     &   b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &   b2stbc_fchy_dia_coreonly, secmodel, lbndusr,
     &   area_fix, sheath_drift_fix)
        if(wrong_flow) write(*,*)
     &       'b2stbc: wrong_flow returned from b2stbc_bas'
      endif
      if(lfeedback.or.b2stbc_feedback.ne.0)
     1 call b2stbc_fb_post(nx, ny, ns, ismain,
     1  na, ne, ni, te, ti, fna, fhe, fhi, fch,
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  main_call, diagno)

* ..compute conditions on the cuts
*     (no flow)
*   ..compute suitable large coefficient
      t0 = 1.0e20_R8*b2ssum(n2,vol,1)
*     (the dimension of t0 is m**3/s)
*   ..loop over species
      do is = 0, ns-1
*    ..set parallel momentum source
       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,0,is), 1)
       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,1,is), 1)
       call b2xxgs (nx, ny, 0.0_R8, smo0(-1,-1,2,is), 1)
       call b2xxgs (nx, ny, -t0, smo0(-1,-1,3,is), 1)
      enddo

c added 2001.11.26 dpc/xpb
c to fix problems with the core energy boundary condition arising
c from diamagnetic drifts
      testqe=0.0_R8
      testqi=0.0_R8
      if(maxval(facdrift).ne.0.0_R8) then
        do ix=-1,nx
          iy=0                                        ! south side
          if(.not.on_closed_surface(ix,iy) .and.
     &       b2stbc_fchy_dia_coreonly.eq.0) then  ! code for private flux
            do is=0,ns-1
              if(.not.is_neutral(is)) then
                t0 =
     &           -2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(rza(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            *na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            *te(bottomix(ix,iy),bottomiy(ix,iy))
     &            /ti(bottomix(ix,iy),bottomiy(ix,iy))
     &           +rza(ix,iy,is)*na(ix,iy,is)*te(ix,iy)/ti(ix,iy))/
     1           2.0_R8*vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                she0(ix,-1,0) = she0(ix,-1,0) + t0
                testqe = testqe + t0
                t0 =
     &           +2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                shi0(ix,-1,0) = shi0(ix,-1,0) + t0
                testqi = testqi + t0
              endif
            enddo
          elseif(on_closed_surface(ix,iy)) then
            if (diagno.ge.2) write(*,*) 'Before she0,shi0 =',
     &            ix,she0(ix,-1,0),shi0(ix,-1,0)
            she0(ix,-1,0)=she0(ix,-1,0)+qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            testqe=testqe+qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            shi0(ix,-1,0)=shi0(ix,-1,0)+qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            testqi=testqi+qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            if (diagno.ge.2) write(*,*) 'Electron term =',
     .       +qeaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            if (diagno.ge.2) write(*,*) 'Ion term =',
     .       +qiaverage(1)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
          endif
          iy = ny
          if(on_closed_surface(ix,iy)) then  ! code for island center
            she0(ix,iy,0)=she0(ix,iy,0)-qeaverage(2)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
            shi0(ix,iy,0)=shi0(ix,iy,0)-qiaverage(2)*facdrift(ix,iy)*
     7       (OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     8        OnedBsqc(ix,iy))*b2stbc_fchy_dia
          elseif(region(ix,iy,0).ne.5 .and.
     &           b2stbc_fchy_dia_coreonly.eq.0) then  ! North boundary
            do is = 0, ns-1
              if (.not.is_neutral(is)) then
                t0 =
     &           +2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(rza(bottomix(ix,iy),bottomiy(ix,iy),is)
     &           *na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &           *te(bottomix(ix,iy),bottomiy(ix,iy))
     &           /ti(bottomix(ix,iy),bottomiy(ix,iy))
     &           +rza(ix,iy,is)*na(ix,iy,is)*te(ix,iy)/ti(ix,iy))/
     1           2.0_R8*vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                she0(ix,iy,0)=she0(ix,iy,0)+t0
                testqe=testqe+t0
                t0 =
     &           -2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
                shi0(ix,iy,0)=shi0(ix,iy,0)+t0
                testqi=testqi+t0
              endif
            enddo
          endif
        enddo
      endif
      if (diagno.ge.2) write(*,*) 'testqe, testqi ', testqe, testqi
c end 2001.11.26 dpc/xpb addition

c corrections for N & S current boundary conditions
c if we want the integral current to be 0
      if(integral_current.gt.0.0_R8) then
        regcnt=0
        do iy=-1,ny,ny+1                         ! do S + N boundaries
          element_done=0
          ix=-1
          do while(ix.le.nx)                     ! do each ix
            if(element_done(ix).eq.0) then       ! have we done this chain yet?
              regcnt=regcnt+1
              sum_sch0=sch0(ix,iy,0)+
     +         sch0(ix,iy,1)*po(ix,iy)+
     +         sch0(ix,iy,2)*ne(ix,iy)+
     +         sch0(ix,iy,3)*ne(ix,iy)*po(ix,iy)
              sum_vol=vol(ix,iy)
              element_done(ix)=regcnt
              ixx=rightix(ix,iy)
              do while(ixx.ge.-1.and.ixx.le.nx)  ! do chain of elements connected to this element
                if(element_done(ixx).ne.0) then
                  ixx=nx+1                       ! get here for periodic case
                else
                  sum_sch0=sum_sch0+sch0(ixx,iy,0)+
     +             sch0(ixx,iy,1)*po(ixx,iy)+
     +             sch0(ixx,iy,2)*ne(ixx,iy)+
     +             sch0(ixx,iy,3)*ne(ixx,iy)*po(ixx,iy)
                  sum_vol=sum_vol+vol(ixx,iy)
                  element_done(ixx)=regcnt
                  ixx=rightix(ixx,iy)
                endif
              enddo
              write(*,*) 'Subtracting ',sum_sch0,' for region ',regcnt
              do ixx=-1,nx
                if(element_done(ixx).eq.regcnt)
     1           sch0(ixx,iy,0)=sch0(ixx,iy,0)-
     2           sum_sch0*vol(ixx,iy)/sum_vol*integral_current
              enddo
            endif
            ix=ix+1
          enddo
          if (ncall.eq.0 .or. diagno.ge.2)
     >     write(*,'(i5/(10(1x,10i1)))') iy,element_done
        enddo
      endif

* ..return
      if((ncall.eq.0 .and. diagno.ge.2) .or. diagno.ge.3) then
        do ix=-1,nx
          write(*,'(a,i4,1p,6(1x,e14.7))') 'S:ix,ti,ni,shi0 ',
     1     ix,ti(ix,-1)/ev,ni(ix,-1,0),(shi0(ix,-1,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'S:ix,te,ne,she0 ',
     1     ix,te(ix,-1)/ev,ne(ix,-1),(she0(ix,-1,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'N:ix,ti,ni,shi0 ',
     1     ix,ti(ix,ny)/ev,ni(ix,ny,0),(shi0(ix,ny,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'N:ix,te,ne,she0 ',
     1     ix,te(ix,ny)/ev,ne(ix,ny),(she0(ix,ny,is),is=0,3)
        enddo
        do iy=-1,ny
          write(*,'(a,i4,1p,6(1x,e14.7))') 'W:iy,ti,ni,shi0 ',
     1     iy,ti(-1,iy)/ev,ni(-1,iy,0),(shi0(-1,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'W:iy,te,ne,she0 ',
     1     iy,te(-1,iy)/ev,ne(-1,iy),(she0(-1,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'E:iy,ti,ni,shi0 ',
     1     iy,ti(nx,iy)/ev,ni(nx,iy,0),(shi0(nx,iy,is),is=0,3)
          write(*,'(a,i4,1p,6(1x,e14.7))') 'E:iy,te,ne,she0 ',
     1     iy,te(nx,iy)/ev,ne(nx,iy),(she0(nx,iy,is),is=0,3)
        enddo
      endif
c      write(*,*) 'Copying sna0 to b2stbc_sna'
      b2stbc_sna(-1:nx,-1:ny,0:ns-1)=
     1 sna0(-1:nx,-1:ny,0,0:ns-1)+
     2 sna0(-1:nx,-1:ny,1,0:ns-1)*na(-1:nx,-1:ny,0:ns-1)
c      write(*,*) 'Copying smo0 to b2stbc_smo'
      do is=0,ns-1
        b2stbc_smo(-1:nx,-1:ny,is)=
     1   smo0(-1:nx,-1:ny,0,is)+
     2   smo0(-1:nx,-1:ny,1,is)*ua(-1:nx,-1:ny,is)+
     3   smo0(-1:nx,-1:ny,2,is)*am(is)*mp*na(-1:nx,-1:ny,is)+
     4   smo0(-1:nx,-1:ny,3,is)*am(is)*mp*na(-1:nx,-1:ny,is)*
     5   ua(-1:nx,-1:ny,is)
      enddo
c      write(*,*) 'Copying she0 to b2stbc_she'
      b2stbc_she(-1:nx,-1:ny)=
     1 she0(-1:nx,-1:ny,0)+
     2 she0(-1:nx,-1:ny,1)*te(-1:nx,-1:ny)+
     3 she0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 she0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*te(-1:nx,-1:ny)
c      write(*,*) 'Copying shi0 to b2stbc_shi'
      b2stbc_shi(-1:nx,-1:ny)=
     1 shi0(-1:nx,-1:ny,0)+
     2 shi0(-1:nx,-1:ny,1)*ti(-1:nx,-1:ny)+
     3 shi0(-1:nx,-1:ny,2)*ni(-1:nx,-1:ny,0)+
     4 shi0(-1:nx,-1:ny,3)*ni(-1:nx,-1:ny,0)*ti(-1:nx,-1:ny)
c      write(*,*) 'Copying sch0 to b2stbc_sch'
      b2stbc_sch(-1:nx,-1:ny)=
     1 sch0(-1:nx,-1:ny,0)+
     2 sch0(-1:nx,-1:ny,1)*po(-1:nx,-1:ny)+
     3 sch0(-1:nx,-1:ny,2)*ne(-1:nx,-1:ny)+
     4 sch0(-1:nx,-1:ny,3)*ne(-1:nx,-1:ny)*po(-1:nx,-1:ny)
c      write(*,*) 'Copying sne0 to b2stbc_sne'
      b2stbc_sne(-1:nx,-1:ny)=
     1 sne0(-1:nx,-1:ny,0)+
     2 sne0(-1:nx,-1:ny,1)*ne(-1:nx,-1:ny)
c      write(*,*) 'Copied s..0 to b2stbc_s..'

! { IYS 04.06.2018
      sna_stbc_neut = 0.0_R8
      if (ank_tracing .gt. 0) then
        do is = 1,nspecies  ! is stands here for the ls - the index of isonuclear sequences.  Do not mix with the fluid index!!!
          do iy = 1, nbc  ! iy stands here for the boundary index
            do ix = 1, BC_LIST_SIZE(iy) ! ix stands here for the index which conuts cells along the boundary
!           avoid corner cells
              if ( .not. (
     & ( leftix(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.-2 .and.
     &       bottomiy(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.-2) .or.
     & ( leftix(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.-2 .and.
     &       topiy(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.ny+1) .or.
     & ( rightix(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.nx+1 .and.
     &       bottomiy(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.-2) .or.
     & ( rightix(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.nx+1 .and.
     &       topiy(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy)).eq.ny+1) ) )
     &        sna_stbc_neut(is,iy) = sna_stbc_neut(is,iy) +
     &           sna0(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy),0,EB2SPCR(is)) +
     &           sna0(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy),1,EB2SPCR(is)) *
     &           na(BC_LIST_X(ix,iy),BC_LIST_Y(ix,iy),EB2SPCR(is))
            enddo
          enddo
        enddo
      endif     ! IYS 04.06.2018 }

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stbc_sna0to1 = sna0
        b2stbc_smo0to3 = smo0
        b2stbc_she0to3 = she0
        b2stbc_shi0to3 = shi0
      endif

      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc

      end subroutine b2stbc

*-----------------------------------------------------------------------
*.specification
      function bcint(fun,nx,ny,il,ib)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_boundary_namelist
      implicit none
*   ..output result (undefined upon entry)
      real (kind=R8) :: bcint
*   ..input arguments (unchanged on exit)
      integer nx, ny, il, ib
      real (kind=R8) :: fun(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     BCINT provides interpolated face-centred values for cell-centred
*     quantities (array fun), for use in boundary conditions.
*     The interpolation is done according to distance away from the
*     cell face. The boundary index is ib, and the cell index along
*     the boundary is il.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..procedures
      external xertst, xerrab

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*  ..test input parameters
      call xertst (0.le.nx.and.0.le.ny, 'faulty parameter nx, ny')
      call xertst (0.lt.ib.and.nbc.ge.ib, 'faulty parameter ib')
      call xertst (1.le.il.and.bc_list_size(ib).ge.il,
     .  'faulty parameter il')

* ..main computation
      if (bcchar(ib).eq.'S') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hy1(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hy1(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hy1(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hy1(topix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        topiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'N') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hy1(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hy1(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hy1(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hy1(bottomix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        bottomiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'E') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hx(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       leftiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        leftiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hx(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hx(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hx(leftix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       leftiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      elseif (bcchar(ib).eq.'W') then
        bcint =
     .   (fun(bc_list_x(il,ib),bc_list_y(il,ib))*
     .    hx(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       rightiy(bc_list_x(il,ib),bc_list_y(il,ib)))+
     .    fun(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .        rightiy(bc_list_x(il,ib),bc_list_y(il,ib)))*
     .    hx(bc_list_x(il,ib),bc_list_y(il,ib)))
     .  /(hx(bc_list_x(il,ib),bc_list_y(il,ib))+
     .    hx(rightix(bc_list_x(il,ib),bc_list_y(il,ib)),
     .       rightiy(bc_list_x(il,ib),bc_list_y(il,ib))))
      else
        call xerrab ('bcint : Unknown orientation!')
      endif

* ..return
      return

*-----------------------------------------------------------------------
*.end bcint

      end function bcint
