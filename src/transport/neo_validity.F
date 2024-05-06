*-----------------------------------------------------------------------
*.specification

      subroutine neo_validity(mpg, geo, state, ns, safetyfactor,  !cjk
     &  nuab, nua, nustar, rhol, lnn, ltt, deltal, wb, deltaw,
     &  deltawn, deltawt, omega, deltab, deltapar, vtha)
      use b2mod_types
      use b2mod_constants
      use b2mod_ranges
      use b2mod_neoclassical
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2plot_debug
      use b2mod_subsys
      use b2us_map
      use b2us_geo
      use b2us_plasma
      implicit none
*   ..input arguments (unchanged on exit)
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(B2State), intent(in) :: state
      integer ns
      real (kind=R8) ::
     *   safetyfactor(mpg%nFt)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   nuab(mpg%nCv,0:ns,0:ns), nua(mpg%nCv,0:ns),
     &   nustar(mpg%nCv,0:ns), rhol(mpg%nCv,0:ns),
     &   lnn(mpg%nCv,0:ns), ltt(mpg%nCv,0:ns),
     &   deltal(mpg%nCv,0:ns), wb(mpg%nCv,0:ns),
     &   deltaw(mpg%nCv,0:ns), deltawn(mpg%nCv,0:ns),
     &   deltawt(mpg%nCv,0:ns), omega(mpg%nCv,0:ns),
     &   deltab(mpg%nCv,0:ns), deltapar(mpg%nCv,0:ns),
     &   vtha(mpg%nCv,0:ns)

      real (kind=R8) :: min_value
c     minimum value for real arrays to avoid division by zero
      parameter (min_value=1.e-8_R8)
c     conversion of constants to CGS units
      real (kind=R8) :: erg, cm, gauss, gramm
      parameter (erg=1.e7_R8, cm=1.e2_R8, gauss=1.e4_R8, gramm=1.e3_R8)
      real (kind=R8) :: eV_cgs, e_cgs, c_cgs
      parameter (eV_cgs=ev*erg,
     &           e_cgs=qe*c*10._R8, c_cgs=c*cm)
      real (kind=R8) ::
     *   bav, b2av, bi2a, rbt, bgradp, dpsidr, rnq, fc, gclass,
     &   fm(0:nmaxgr-1), r2i, cosbl(0:nmaxgr-1),
     &   cosb2(0:nmaxgr-1), sinbl(0:nmaxgr-1), sinb2(0:nmaxgr-1)
      common /b2ngeom/ bav,b2av,bi2a,rbt,bgradp,dpsidr,rnq,fc,
     *                 gclass,fm,r2i,cosbl,cosb2,sinbl,sinb2

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     NEO_VALIDITY computes diagnostics parameters for the neoclassical
*     transport
*
*
*
*  3. description
*
*     This routine computes the following diagnostics parameters:
*
*     - vtha   : the thermal velocity of species a
*     - nuab   : the collision frequency of species a with species b
*     - nua    : the sum of nuab over all species b
*     - nustar : the normalized collisionality
*     - rhol   : the Larmor radius
*     - lnn    : the density gradient length (in y-direction)
*     - ltt    : the temperature gradient length (in y-direction)
*     - deltal : rhol/min{lnn,ltt}
*     - wb     : the radial width of the Banana orbit
*     - deltaw : wb/min{lnn,ltt}
*     - deltawn: wb/lnn
*     - deltawt: wb/ltt
*     - omega  : the cyclotron frequency
*     - deltab : nua/omega
*     - deltapar: ft^4*nustar*wb*R*nua/(min{lnn,ltt}^2*omega)
*                 (a measure for the density perturbations along a flux
*                  surface, should be << 1 for neoclassical theory)
*     all of the above parameters are calculated for the charged particles
*     only and set to zero for all neutrals
*     wb is valid only in the Banana regime (nustar<1)
*     deltapar is important only in the PS regime (nustar>eps^-3/2)
*     all delta ratios are only sensible in the outboard midplane due
*     to the one-dimensional nature of the implemented neoclassical model
*-----------------------------------------------------------------------
*   ..local variables
      integer iCv, iFt, i, is, is2
      real (kind=R8) ::
     &   gradient(mpg%nCv), funv(mpg%nVx), lnab, minlnlt, crmin, crmax
      real (kind=R8) ::
     &   rmaj_(mpg%nFt)

*   ..procedures
      external xertst, b2xvsg

*   ..subprogram start-up calls
      call subini ('neo_validity')
*   ..test ns
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
*   ..test bb
      call b2xvsg (mpg%nCv, geo%cvBb(1,3), 1, 'bb3', '.gt.')
*   ..test state
      call b2xvsg (mpg%nCv*ns, state%pl%na, 1, 'na', '.gt.')
      call b2xvsg (mpg%nCv, state%pl%te, 1, 'te', '.gt.')
      call b2xvsg (mpg%nCv, state%pl%ti, 1, 'ti', '.gt.')
      call b2xvsg (mpg%nCv*2, state%dv%ni, 1, 'ni', '.gt.')
      call b2xvsg (mpg%nCv, state%dv%ne, 1, 'ne', '.gt.')

*   ..main

*   ..intialize arrays (set to zero for neutrals)
      vtha(:,:)=0._R8
      nuab(:,:,:)=0._R8
      nua(:,:)=0._R8
      nustar(:,:)=0._R8
      rhol(:,:)=0._R8
      lnn(:,:)=0._R8
      ltt(:,:)=0._R8
      deltal(:,:)=0._R8
      wb(:,:)=0._R8
      deltaw(:,:)=0._R8
      deltawn(:,:)=0._R8
      deltawt(:,:)=0._R8
      omega(:,:)=0._R8
      deltab(:,:)=0._R8
      deltapar(:,:)=0._R8
      rmaj_(:)=0._R8

*   ..calculate major radius of closed flux surfaces
      do iFt=1,mpg%nFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if(.not.mpg%cvOnClosedSurface(iCv)) cycle
        crmin=geo%cvX(iCv)
        crmax=crmin
        do i=2,mpg%ftCvP(iFt,2)
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          crmin=min(crmin,geo%cvX(iCv))
          crmax=max(crmax,geo%cvX(iCv))
        enddo
        rmaj_(iFt)=0.5_R8*(crmin+crmax)
      enddo

*   ..calculate gradient lengths lnn and ltt
      do is=0,ns-1
       call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     &  state%pl%na(:,is), funv, gradient)
       lnn(:,is)=max(state%pl%na(:,is),min_value)/
     &           max(abs(gradient(:)),min_value)
       if (is_neutral(is).and.nint(zn(is)).eq.1) then
         call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     &    state%pl%tn(:), funv, gradient)
         ltt(:,is)=max(state%pl%tn(:),min_value*ev)/
     &             max(abs(gradient(:)),min_value*ev)
       else
         call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     &    state%pl%ti(:), funv, gradient)
         ltt(:,is)=max(state%pl%ti(:),min_value*ev)/
     &             max(abs(gradient(:)),min_value*ev)
       end if
      end do
      call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     & state%dv%ne(:), funv, gradient)
      lnn(:,ns)=max(state%dv%ne(:),min_value)/
     &          max(abs(gradient(:)),min_value)
      call gradc_r(mpg%nCv, mpg%nFc, mpg%nVx, 0, geo, mpg,
     & state%pl%te(:), funv, gradient)
      ltt(:,ns)=max(state%pl%te(:),min_value*ev)/
     &          max(abs(gradient(:)),min_value*ev)

*   ..calculate themal velocities
      do is=0,ns-1
        if(is_neutral(is).and.nint(zn(is)).eq.1) then
          vtha(:,is)=sqrt(erg*state%pl%tn(:)/(am(is)*mp*gramm))
        else
          vtha(:,is)=sqrt(erg*state%pl%ti(:)/(am(is)*mp*gramm))
        endif
      enddo
      vtha(:,ns)=sqrt(erg*state%pl%te(:)/(me*gramm))

*   ..calculate collision frequency nuab of species a and b
*      ion-ion collisions
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do is2=0,ns-1
         if (.not. is_neutral(is2)) then
          do iCv=1,mpg%nCv
*      Coulomb logarithm
            lnab=17.3-0.5*log(state%dv%ne(iCv)/1.e20)
     &               +1.5*log(state%pl%ti(iCv)*erg/(eV_cgs*1.e3))
            nuab(iCv,is,is2)=16.*sqrt(pi)*state%pl%na(iCv,is2)/1.e6
     &        *state%rt%rz2(iCv,is)*state%rt%rz2(iCv,is2)*e_cgs**4./
     &        (3.*am(is)*mp*gramm)*lnab*
     &        (1./(am(is)*mp*gramm)+1./(am(is2)*mp*gramm))/
     &        (vtha(iCv,is)**2+vtha(iCv,is2)**2)**1.5
          enddo
         endif
        enddo
       endif
      enddo
*      ion-electron collisions and electron-ion collisions
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do iCv=1,mpg%nCv
*      Coulomb logarithm
          lnab=15.2-0.5*log(state%dv%ne(iCv)/1.e20)
     &             +log(state%pl%te(iCv)*erg/(eV_cgs*1.e3))
          nuab(iCv,is,ns)=16.*sqrt(pi)*state%dv%ne(iCv)/1.e6
     &      *state%rt%rz2(iCv,is)
     &      *e_cgs**4./(3.*am(is)*mp*gramm)*lnab*
     &      (1./(am(is)*mp*gramm)+1./(me*gramm))/
     &      (vtha(iCv,is)**2+vtha(iCv,ns)**2)**1.5
          nuab(iCv,ns,is)=16.*sqrt(pi)*state%pl%na(iCv,is)/1.e6
     &      *state%rt%rz2(iCv,is)
     &      *e_cgs**4./(3.*me*gramm)*lnab*
     &      (1./(am(is)*mp*gramm)+1./(me*gramm))/
     &      (vtha(iCv,is)**2+vtha(iCv,ns)**2)**1.5
        enddo
       endif
      enddo
*      electron-electron collisions
      do iCv=1,mpg%nCv
*      Coulomb logarithm
        lnab=14.9-0.5*log(state%dv%ne(iCv)/1.e20)
     &           +log(state%pl%te(iCv)*erg/(eV_cgs*1.e3))
        nuab(iCv,ns,ns)=16.*sqrt(pi)*state%dv%ne(iCv)/1.e6
     &    *e_cgs**4./(3.*me*gramm)*lnab*(2./(me*gramm))/
     &    (2.*vtha(iCv,ns)**2)**1.5
      enddo

*   ..calculate total collision frequency nua of species a
      do is=0,ns
       do is2=0,ns
        nua(:,is)=nua(:,is)+nuab(:,is,is2)
       enddo
      enddo

*   ..calculate collisionality nustar of species a
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do iCv=1,mpg%nCv
         iFt=mpg%cvFt(iCv)
         if (iFt.gt.0.and.mpg%cvOnClosedSurface(iCv)) then
           nustar(iCv,is)=nua(iCv,is)*safetyfactor(iFt)*
     &      rmaj_(iFt)*cm/(vtha(iFt,is)*(1.-fc_(iFt))**3)
         else
           nustar(iCv,is)=0.0_R8
         endif
        enddo
       endif
      enddo
      do iCv=1,mpg%nCv
        iFt=mpg%cvFt(iCv)
        if (iFt.gt.0.and.mpg%cvOnClosedSurface(iCv)) then
         nustar(iCv,ns)=nua(iCv,ns)*safetyfactor(iFt)*
     &    rmaj_(iFt)*cm/(vtha(iCv,ns)*(1.-fc_(iFt))**3)
        else
         nustar(iCv,ns)=0.0_R8
        endif
      enddo

*   ..calculate validity tests
      do iCv=1,mpg%nCv
        iFt=mpg%cvFt(iCv)
        do is=0,ns-1
         if (.not. is_neutral(is)) then
*      minimum gradient length
          if (lnn(iCv,is).lt.ltt(iCv,is)) then
           minlnlt=lnn(iCv,is)
          else
           minlnlt=ltt(iCv,is)
          endif
*      Larmor radius
          rhol(iCv,is)=am(is)*mp*gramm*c_cgs*vtha(iCv,is)/
     &      (state%rt%rza(iCv,is)*e_cgs*geo%cvBb(iCv,3)*gauss*cm)
*      ratio of Larmor radius to minimum gradient length
          deltal(iCv,is)=rhol(iCv,is)/minlnlt
*      width of Banana orbit
          if (iFt.gt.0)
     >     wb(iCv,is)=am(is)*mp*gramm*c_cgs*(1.-fc_(iFt))*
     &      vtha(iCv,is)/(state%rt%rza(iCv,is)*e_cgs*
     &                    geo%cvBb(iCv,0)*gauss*cm)
*      ratio of Banana orbit to minimum gradient length
          deltaw(iCv,is)=wb(iCv,is)/minlnlt
*      ratio of Banana orbit to density gradient length
          deltawn(iCv,is)=wb(iCv,is)/lnn(iCv,is)
*      ratio of Banana orbit to temperature gradient length
          deltawt(iCv,is)=wb(iCv,is)/ltt(iCv,is)
*      ion gyro-frequency
          omega(iCv,is)=state%rt%rza(iCv,is)*e_cgs*geo%cvBb(iCv,3)*
     &     gauss/(am(is)*mp*gramm*c_cgs)
*      ratio of the total collision frequency to the gyro-frequency
          deltab(iCv,is)=nua(iCv,is)/omega(iCv,is)
*      measure for the density perturbation in a collision dominated
*      regime
          if (iFt.gt.0)
     >     deltapar(iCv,is)=(1.0_R8-fc_(iFt))**4*nustar(iCv,is)*
     &      wb(iCv,is)*rmaj_(iFt)*nua(iCv,is)/(minlnlt**2*
     &      omega(iCv,is))
         endif
        enddo
*      minimum gradient length
        if (lnn(iCv,ns).lt.ltt(iCv,ns)) then
         minlnlt=lnn(iCv,ns)
        else
         minlnlt=ltt(iCv,ns)
        endif
*      Larmor radius
        rhol(iCv,ns)=me*gramm*c_cgs*vtha(iCv,ns)/
     &    (e_cgs*geo%cvBb(iCv,3)*gauss*cm)
*      ratio of Larmor radius to minimum gradient length
        deltal(iCv,ns)=rhol(iCv,ns)/minlnlt
*      width of Banana orbit
        if (iFt.gt.0)
     >   wb(iCv,ns)=me*gramm*c_cgs*(1.0_R8-fc_(iFt))*
     &    vtha(iCv,ns)/(e_cgs*geo%cvBb(iCv,0)*gauss*cm)
*      ratio of Banana orbit to minimum gradient length
        deltaw(iCv,ns)=wb(iCv,ns)/minlnlt
*      ratio of Banana orbit to density gradient length
        deltawn(iCv,ns)=wb(iCv,ns)/lnn(iCv,ns)
*      ratio of Banana orbit to temperature gradient length
        deltawt(iCv,ns)=wb(iCv,ns)/ltt(iCv,ns)
*      electron gyro-frequency
        omega(iCv,ns)=e_cgs*geo%cvBb(iCv,3)*gauss/(me*gramm*c_cgs)
*      ratio of the total collision frequency to the gyro-frequency
        deltab(iCv,ns)=nua(iCv,ns)/omega(iCv,ns)
*      measure for the density perturbation in a collision dominated
*      regime
        if (iFt.gt.0)
     >   deltapar(iCv,ns)=(1.-fc_(iFt))**4*nustar(iCv,ns)*
     &    wb(iCv,ns)*rmaj_(iFt)*nua(iCv,ns)/(minlnlt**2*
     &    omega(iCv,ns))
      enddo

*   ..return
      call subend ()
      return
      end subroutine neo_validity

      subroutine chang_hinton(mpg,geo,state,state_ext,switch,ns,
     1 safetyfactor,transport)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2plot_debug
      use b2mod_subsys
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_switches
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      type(B2State), intent(in) :: state
      type(B2StateExt), intent(in) :: state_ext
      type(switches), intent(in) :: switch
      integer, intent(in) :: ns
      real (kind=R8) ::
     3 safetyfactor(mpg%nFt),
     2 transport(mpg%nCv,0:ns-1)

      integer iCv,iFt,i,is,icmax,is_fuel
      real (kind=R8) :: crmin,crmax,czmin,czmax,
     2 bave,epsilon,vthi,omegaci,
     1 rhoi,alpha,deltap,f1,f2,taui,nui,nustari,mustari,chi,
     2 bracket1,bracket2,ammin
      real (kind=R8) :: lnlam(mpg%nCv)
      external xertst

      call subini ('chang_hinton')
      transport=0.0_R8
      call b2tlnl(mpg%nCv,switch,1,state%pl%te,state%pl%ti,state%dv%ne,
     .  lnlam)
      do iFt=1,mpg%nFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if(.not.mpg%cvOnClosedSurface(iCv)) cycle
        if (b2plot_debug.ge.8) write(*,*) 'iFt: ',iFt
        icmax=1
        crmin=geo%cvX(iCv)
        crmax=crmin
        czmin=geo%cvY(iCv)
        czmax=czmin
        do i=2,mpg%ftCvP(iFt,2)
          iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
          if(geo%cvX(iCv).gt.crmax) icmax=iCv
          crmin=min(crmin,geo%cvX(iCv))
          crmax=max(crmax,geo%cvX(iCv))
          czmin=min(czmin,geo%cvY(iCv))
          czmax=max(czmax,geo%cvY(iCv))
        enddo
        if (b2plot_debug.ge.8)
     >   write(*,*) 'crmin,crmax,czmin,czmax',
     1   crmin,crmax,czmin,czmax
        iCv=icmax
        bave=geo%cvX(iCv)*2.0_R8/(crmax+crmin)*geo%cvBb(iCv,3)
        if (b2plot_debug.ge.8) write(*,*) 'ne,te,ti,b,q: ',
     1   state%dv%ne(iCv),state%pl%te(iCv)/qe,state%pl%ti(iCv)/qe,
     2   bave,safetyfactor(iFt)
        epsilon=(crmax-crmin)/(crmax+crmin)
        if (b2plot_debug.ge.8) write(*,*) 'epsilon: ',epsilon
        if (b2plot_debug.ge.8) write(*,*) 'lnlam: ',lnlam(iCv)
        deltap=-2._R8
        if (b2plot_debug.ge.8) write(*,*) 'deltap: ',deltap
        f1=(1.0_R8+1.5_R8*(epsilon**2+epsilon*deltap)+
     1   0.375_R8*epsilon**3*deltap)/(1_R8+0.5_R8*epsilon*deltap)
        if (b2plot_debug.ge.8) write(*,*) 'f1: ',f1
        f2=sqrt(1_R8-epsilon**2)*(1.0_R8+0.5_R8*epsilon*deltap)/
     1   (1_R8+deltap/epsilon*(sqrt(1_R8-epsilon**2)-1.0_R8))
        if (b2plot_debug.ge.8) write(*,*) 'f2: ',f2
*   ..calculate alpha
        alpha=0._R8
        ammin=54._R8
        is_fuel=1
        do is=0,ns-1
          if ((.not. is_neutral(is)).and.(am(is).lt.ammin)) then
            is_fuel=is
            ammin=am(is)
          endif
        enddo
        do is=0,ns-1
          if ((.not. is_neutral(is)).and.(is.ne.is_fuel)) then
            alpha=alpha+state%pl%na(iCv,is)*state%rt%rz2(iCv,is)
          endif
        enddo
        do is=0,state_ext%ns-1
          if (.not. state_ext%is_neutral(is)) then
            alpha=alpha+state_ext%na(iCv,is)*state_ext%za2(iCv,is)
          endif
        enddo
        alpha=alpha/(state%pl%na(iCv,is_fuel)*state%rt%rz2(iCv,is_fuel))
        if (b2plot_debug.ge.8) write(*,*) 'alpha: ', alpha
        do is=0,ns-1
          if(.not. is_neutral(is)) then
            if (b2plot_debug.ge.8) write(*,*) 'rza,am,na: ',
     >       state%rt%rza(iCv,is),am(is),state%pl%na(iCv,is)
            vthi=sqrt(state%pl%ti(iCv)/mp/am(is))
            if (b2plot_debug.ge.8) write(*,*) 'vthi: ',vthi
            omegaci=state%rt%rza(iCv,is)*qe*bave/mp/am(is)
            if (b2plot_debug.ge.8) write(*,*) 'omegaci: ',omegaci
            rhoi=sqrt(2.0_R8)*vthi/omegaci
            if (b2plot_debug.ge.8) write(*,*) 'rhoi: ',rhoi
            taui=12.0_R8*pi**1.5_R8*eps0**2*
     1       (am(is)*mp)**0.5_R8*state%pl%ti(iCv)**1.5_R8/
     2       (state%pl%na(iCv,is)*state%rt%rz2(iCv,is)**2*qe**4*
     3        lnlam(iCv))
            if (b2plot_debug.ge.8) write(*,*) 'taui: ',taui
            nui=1.0_R8/taui
            nustari=nui*(crmax+crmin)/2.0_R8*safetyfactor(iFt)/
     1       epsilon**(1.5_R8)/vthi
            if (b2plot_debug.ge.8) write(*,*) 'nustari: ',nustari
            mustari=nustari*(1.0_R8+1.54_R8*alpha)
            if (b2plot_debug.ge.8) write(*,*) 'mustari: ',mustari
            bracket1=(0.66_R8*(1.0_R8+1.54_R8*alpha)+
     3               (1.88_R8*epsilon**0.5_R8-1.54_R8*epsilon)*
     4               (1.0_R8+3.75_R8*alpha))/
     5               (1.0_R8+1.03_R8*mustari**0.5_R8+0.31*mustari)
            bracket2=0.59_R8*mustari*epsilon/
     6         (1_R8+0.74_R8*mustari*epsilon**1.5_R8)*
     7         (1_R8+
     8           1.33_R8*alpha*(1.0_R8+0.60_R8*alpha)/
     9           (1.0_R8+1.79_R8*alpha))
            if (b2plot_debug.ge.8) write(*,*) 'bracket1 =', bracket1
            if (b2plot_debug.ge.8) write(*,*) 'bracket2 =', bracket2
            chi=epsilon**(-1.5_R8)*
     1       safetyfactor(iFt)**2*rhoi**2/taui*
     2       ((0.66_R8*(1.0_R8+1.54_R8*alpha)+
     3               (1.88_R8*epsilon**0.5_R8-1.54_R8*epsilon)*
     4               (1_R8+3.75_R8*alpha))/
     5               (1.0_R8+1.03_R8*mustari**0.5_R8+0.31*mustari)*f1
     5        +0.59_R8*mustari*epsilon/
     6         (1_R8+0.74_R8*mustari*epsilon**1.5_R8)*
     7         (1_R8+
     8           1.33_R8*alpha*(1.0_R8+0.60_R8*alpha)/
     9           (1.0_R8+1.79_R8*alpha))*(f1-f2))
            if (b2plot_debug.ge.8) write(*,*) 'chi: ',chi
            do i=1,mpg%ftCvP(iFt,2)
              iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
              transport(iCv,is)=chi
            enddo
          endif
        enddo
      enddo
      call subend ()
      return
      end subroutine chang_hinton















