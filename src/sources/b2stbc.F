*-----------------------------------------------------------------------
*.specification

      subroutine b2stbc (nx, ny, ns, ismain, BoRiS, facdrift, fac_ExB,   !xpb
     &  vol, crx, hx, hy, hz, qz, qc, gs, pbs, pbshz, bb,                !xpb !srv 09.01.01
     &  na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     &  te, ti, po, ne, ni, kinrgy, floe_noc, floi_noc,                  !srv 28.10.02
     1  fna, fna_mdf, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52, 
     &  fhe, fhi, fhe_mdf, fhi_mdf, fch, fne, fni, fna_fcor,             !srv 16.02.12
     &  fnaPSch, fhePSch, fhiPSch,                                       !srv 25.02.00 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,   !srv 21.10.02
     2  cdna, chce, chve, chci, chvi, calf,
     &  sna0, smo0, she0, shi0, sch0, wrong_flow, main_call)             !xpb
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_boundary_namelist
      use b2mod_geo_corner
      use b2mod_constants
      use b2mod_sources
     1 , only : b2stbc_sna,b2stbc_smo,b2stbc_she,b2stbc_shi,b2stbc_sch
      use b2mod_feedback
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *  BoRiS, facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3), hx(-1:nx,-1:ny),
     &  hy(-1:nx,-1:ny), hz(-1:nx,-1:ny), pbshz(-1:nx,-1:ny,0:1),       !srv 09.01.01
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     1  bb(-1:nx,-1:ny,0:3), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1), kinrgy(-1:nx,-1:ny,0:ns-1),
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
     &  fch(-1:nx,-1:ny,0:1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &  fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),              !srv 11.10.02 21.10.02 {
     &  fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &  fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &  fchinert(-1:nx,-1:ny,0:1),                                       !srv 21.10.02 }
     &  fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),                                !srv 25.02.00
     &  fnaPSch(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 13.10.06
     &  fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),              !srv 13.10.06
     &  floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1)             !srv 28.10.02
      real (kind=R8) ::
     1  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     2  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     3  calf(-1:nx,-1:ny,0:1), cdna(-1:nx,-1:ny,0:1,0:ns-1)
      logical main_call                                                  !xpb
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3)
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC provides linearised expressions for the sources of
*     density, momentum, electron heat, atom heat, and charge, due to
*     the boundary conditions.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, regcnt, ixx
      real (kind=R8) ::
     *   sna0ep, she0ep, shi0ep, t0, testqe, testqi,
     *   b2stbc_neoclassical, b2stbc_fchy_dia,
     &   praverage(2), prsum(2), gssum(2), tesum(2), tisum(2),
     1   nasum(0:ns-1,2), teave(2), tiave(2),
     1   naave(0:ns-1,2), nef, tef, tif, naf(0:ns-1),
     1   volsum(2), qeaverage(2), qiaverage(2), qesum(2), qisum(2),
     1   sum_sch0, sum_vol, integral_current
      integer boundary_namelist, diagno, b2stbc_fchy_dia_coreonly,
     *   sound_model, area_fix
      integer element_done(-1:nx)
      integer b2stbc_feedback, secmodel
      save ncall, sna0ep, she0ep, shi0ep, boundary_namelist, diagno,
     *   b2stbc_neoclassical, b2stbc_fchy_dia, integral_current,
     1   b2stbc_feedback, b2stbc_fchy_dia_coreonly, secmodel,
     *   sound_model, area_fix
*   ..procedures
      intrinsic mod, nint
      real (kind=R8) ::
     *   ssum, hy1
      external subini, subend, xertst, sfill, ssum, saxpy,
     *   ipgetr, ipgeti
      external b2xvsg, b2xvff, b2xvps, b2xxgs,
     1   b2stbc_phys, b2stbc_bas, b2stbc_fb
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, sna0ep/1.0e-36_R8/, she0ep/1.0e-36_R8/,
     1   shi0ep/1.0e-36_R8/, boundary_namelist/0/,diagno/0/,
     2   b2stbc_neoclassical/0.0_R8/, b2stbc_fchy_dia/1.0_R8/,
     3   integral_current/0.0_R8/, b2stbc_feedback/0/, area_fix/3/,
     4   b2stbc_fchy_dia_coreonly/1/, secmodel/0/, sound_model/0/
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
       call ipgeti ('b2stbc_boundary_namelist', boundary_namelist)
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
       call ipgeti ('b2stbc_sound_model', sound_model)
       call ipgeti ('b2news_area_fix', area_fix)
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
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
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
!      if (ncall.eq.1 .and. main_call) then                               !xpb !srv 18.05.09
!        if (maxval(facdrift).eq.0.0_R8) then                             !xpb
!          write (*,*) 'b2stbc_neoclassical =', b2stbc_neoclassical       !xpb
!        else                                                             !xpb
!          write (*,*) 'b2stbc_fchy_dia =', b2stbc_fchy_dia               !xpb
!        endif                                                            !xpb
!      endif
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
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)

* ..compute standard form volume sources
*   (This code placed here for want of a better location)
*   ..low-level particle source
      do is = 0, ns-1
       call saxpy (n2, sna0ep, vol, 1, sna0(-1,-1,0,is), 1)
      enddo
      call saxpy (n2, she0ep, vol, 1, she0(-1,-1,0), 1)
      call saxpy (n2, shi0ep, vol, 1, shi0(-1,-1,0), 1)

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
        if(mod(region(ix,iy,0),4).eq.1) then                             !xpb
          tef = (te(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           te(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ti(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           ne(topix(ix,iy),topiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(topix(ix,iy),topiy(ix,iy))+
     &           na(topix(ix,iy),topiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(topix(ix,iy),topiy(ix,iy))+hy1(ix,iy))
          tesum(1)=tesum(1) + tef*gs(topix(ix,iy),topiy(ix,iy),1)
          tisum(1)=tisum(1) + tif*gs(topix(ix,iy),topiy(ix,iy),1)
          nasum(:,1)=nasum(:,1) + naf(:)*gs(topix(ix,iy),topiy(ix,iy),1)
          prsum(1)=prsum(1) +
     1     (ne(ix,iy)*te(ix,iy)+ni(ix,iy,1)*ti(ix,iy))*vol(ix,iy)
     2     *bb(ix,iy,2)*vol(ix,iy)/gs(ix,iy,2)
          qesum(1)=qesum(1)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     ne(ix,iy)*te(ix,iy)**2*bb(ix,iy,2)/qe*
     2     gs(topix(ix,iy),topiy(ix,iy),1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(1)=qisum(1)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         na(ix,iy,is)*ti(ix,iy)**2*bb(ix,iy,2)/
     5         (qe*rza(ix,iy,is))*gs(topix(ix,iy),topiy(ix,iy),1)
            endif
          enddo
          gssum(1)=gssum(1) + gs(topix(ix,iy),topiy(ix,iy),1)            !xpb
          volsum(1)=volsum(1) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                      !sv 24.03.99
      praverage(1)=prsum(1)/volsum(1)                                    !xpb
      qeaverage(1)=qesum(1)/gssum(1)
      qiaverage(1)=qisum(1)/gssum(1)
      teave(1)=tesum(1)/gssum(1)
      tiave(1)=tisum(1)/gssum(1)
      naave(:,1)=nasum(:,1)/gssum(1)
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
        if(region(ix,iy,0).eq.5) then                                    !xpb
          tef = (te(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           te(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tif = (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          nef = (ne(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           ne(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          naf(:) = (na(ix,iy,:)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &           na(bottomix(ix,iy),bottomiy(ix,iy),:)*hy1(ix,iy))/
     &          (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
          tesum(2)=tesum(2) + tef*gs(ix,iy,1)                            !xpb
          tisum(2)=tisum(2) + tif*gs(ix,iy,1)                            !xpb
          nasum(:,2)=nasum(:,2) + naf(:)*gs(ix,iy,1)                     !xpb
          prsum(2)=prsum(2) +
     1     (ne(ix,iy)*te(ix,iy)+ni(ix,iy,1)*ti(ix,iy))*vol(ix,iy)
     2     *bb(ix,iy,2)*vol(ix,iy)/gs(ix,iy,2)
          qesum(2)=qesum(2)+2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     1     ne(ix,iy)*te(ix,iy)**2*bb(ix,iy,2)/qe*gs(ix,iy,1)
          do is = 0, ns-1
            if (.not.is_neutral(is)) then
              qisum(2)=qisum(2)-2.5_R8*vol(ix,iy)/gs(ix,iy,2)*
     3         na(ix,iy,is)*ti(ix,iy)**2*bb(ix,iy,2)/
     5         (qe*rza(ix,iy,is))*gs(ix,iy,1)
            endif
          enddo
          gssum(2)=gssum(2) + gs(ix,iy,1)                                !xpb
          volsum(2)=volsum(2) + vol(ix,iy)
        endif                                                            !xpb
      enddo                                                              !xpb
      if (gssum(2).gt.0.0_R8) then
        praverage(2)=prsum(2)/volsum(2)                                  !xpb
        qeaverage(2)=qesum(2)/gssum(2)                                   !xpb
        qiaverage(2)=qisum(2)/gssum(2)                                   !xpb
        teave(2)=tesum(2)/gssum(2)
        tiave(2)=tisum(2)/gssum(2)
        naave(:,2)=nasum(:,2)/gssum(2)
      endif
      if (ncall.eq.0) lbndusr=.false.                                    !xpb

      if(boundary_namelist.ge.1) then
        if(ncall.eq.0) write(*,*) ' b2stbc_phys called'
        call b2stbc_phys(nx, ny, ns, ismain,
     1   BoRiS, facdrift, fac_ExB, vol, crx, hx, hy, hz, qz, qc, gs,   !srv 24.08.08
     1   pbs, pbshz, bb, na, ua, uadia, vedia, vadia, veecrb, vaecrb,
     2   te, ti, po, ne, ni, kinrgy,  floe_noc, floi_noc, fna, fna_mdf, !srv 28.10.02
     1   fnaPSch, fna_32, fna_52, fni_32, fni_52, fne_32, fne_52,      !srv 13.10.06 
     &   fhe, fhi, fhe_mdf, fhi_mdf, fna_fcor, fhePSch, fhiPSch, fch,  !srv 16.02.12
     &   fne, fni,                                                     !srv 13.10.06
     &   fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,!srv 21.10.02
     &   cdna, chce, chve, chci, chvi, calf,                           !srv 25.09.09
     &   sna0, smo0, she0, shi0, sch0, wrong_flow, main_call, diagno,
     &   b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &   b2stbc_fchy_dia_coreonly, secmodel, sound_model, area_fix)
        if(lbndusr) then
          if(ncall.eq.0) write(*,*) 'BNDUSR called'
          call b2stbc_bas(nx, ny, ns, ismain,
     1     facdrift, vol, hx, hy, qz, qc, gs,
     1     pbs, bb, na, ua, uadia, te, ti, po, ne, ni,
     1     fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &     sna0, smo0, she0, shi0, sch0, wrong_flow, main_call, diagno,
     &     b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &     b2stbc_fchy_dia_coreonly, secmodel, lbndusr, sound_model, 
     &     area_fix)
        endif      !lbndusr
      else      !boundary_namelist.ge.1
        call b2stbc_bas(nx, ny, ns, ismain,
     1   facdrift, vol, hx, hy, qz, qc, gs,
     1   pbs, bb, na, ua, uadia, te, ti, po, ne, ni,
     1   fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &   sna0, smo0, she0, shi0, sch0, wrong_flow, main_call, diagno,
     &   b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &   b2stbc_fchy_dia_coreonly, secmodel, lbndusr, sound_model, 
     &   area_fix)
      endif
      if(lfeedback.or.b2stbc_feedback.ne.0)
     1 call b2stbc_fb(nx, ny, ns, ismain, BoRiS,
     1  vol, gs, na, ne, ni, te, ti, fna, fhe, fhi, fch, kinrgy,
     &  main_call, diagno)

* ..compute conditions on the cuts
*     (no flow)
*   ..compute suitable large coefficient
      t0 = 1.0e20_R8*ssum(n2,vol,1)
*     (the dimension of t0 is m**3/s)
*   ..loop over species
      do is = 0, ns-1
*    ..set momentum source
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
          if(mod(region(ix,-1,0),4).ne.1 .and.
     &       b2stbc_fchy_dia_coreonly.eq.0) then  ! code for private flux 
            do is=0,ns-1
              if(.not.is_neutral(is)) then
                if(ncall.eq.0)
     >           vadia(ix,iy,1,is) = ti(ix,iy)*bb(ix,iy,2)/            !xpb
     &            (qe*rza(ix,iy,is)*hx(ix,iy))*(2.0_R8/bb(ix,iy,3)**3) !xpb
     &            *(wbbr(ix,iy,3)-wbbl(ix,iy,3))*facdrift(ix,iy)       !xpb
                she0(ix,-1,0)=she0(ix,-1,0)
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
                shi0(ix,-1,0)=shi0(ix,-1,0)
     &           +2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
              endif
            enddo
          elseif(mod(region(ix,-1,0),4).eq.1) then
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
          if(region(ix,iy,0).eq.5 .and. periodic_bc.eq.1 .and.
     &       b2stbc_fchy_dia_coreonly.eq.0) then  ! code for island center
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
                if(ncall.eq.0)
     >           vadia(ix,iy,1,is) = 
     .            ti(bottomix(ix,iy),bottomiy(ix,iy))*
     &            bb(bottomix(ix,iy),bottomiy(ix,iy),2)                  !xpb
     &            /(qe*rza(bottomix(ix,iy),bottomiy(ix,iy),is)*
     &              hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &             (2.0_R8/bb(bottomix(ix,iy),bottomiy(ix,iy),3)**3)     !xpb
     &            *(wbbr(bottomix(ix,iy),bottomiy(ix,iy),3)-
     &              wbbl(bottomix(ix,iy),bottomiy(ix,iy),3))*
     &              facdrift(ix,iy)                                      !xpb
                she0(ix,iy,0)=she0(ix,iy,0)
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
                shi0(ix,iy,0)=shi0(ix,iy,0)
     &           -2.5_R8*
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
     &           /(hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
     &           *(na(bottomix(ix,iy),bottomiy(ix,iy),is)
     &            +na(ix,iy,is))/2.0_R8
     &           *vadia(ix,iy,1,is)*b2stbc_fchy_dia*
     1           (ti(ix,iy)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     &            ti(bottomix(ix,iy),bottomiy(ix,iy))*hy1(ix,iy))/
     &           (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
              endif
            enddo
          endif
        enddo
      endif
      if (diagno.ge.2) write(*,*) 'testqe, testqi ',testqe, testqi
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
          write(*,*) 'S:ix,ti,ni,shi0 ',
     1     ix,ti(ix,-1)/ev,ni(ix,-1,0),shi0(ix,-1,:)
          write(*,*) 'S:ix,te,ne,she0 ',
     1     ix,te(ix,-1)/ev,ne(ix,-1),she0(ix,-1,:)
          write(*,*) 'N:ix,ti,ni,shi0 ',
     1     ix,ti(ix,ny)/ev,ni(ix,ny,0),shi0(ix,ny,:)
          write(*,*) 'N:ix,te,ne,she0 ',
     1     ix,te(ix,ny)/ev,ne(ix,ny),she0(ix,ny,:)
        enddo
        do iy=-1,ny
          write(*,*) 'W:iy,ti,ni,shi0 ',
     1     iy,ti(-1,iy)/ev,ni(-1,iy,0),shi0(-1,iy,:)
          write(*,*) 'W:iy,te,ne,she0 ',
     1     iy,te(-1,iy)/ev,ne(-1,iy),she0(-1,iy,:)
          write(*,*) 'E:iy,ti,ni,shi0 ',
     1     iy,ti(nx,iy)/ev,ni(nx,iy,0),shi0(nx,iy,:)
          write(*,*) 'E:iy,te,ne,she0 ',
     1     iy,te(nx,iy)/ev,ne(nx,iy),she0(nx,iy,:)
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
c      write(*,*) 'Copied s..0 to b2stbc_s..'
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc

      end subroutine b2stbc
