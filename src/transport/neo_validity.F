*-----------------------------------------------------------------------
*.specification

      subroutine neo_validity(nx, ny, ns, jsep, jxa, safetyfactor,  !cjk
     &  nuab, nua, nustar, rhol, lnn, ltt, deltal, wb, deltaw,
     &  deltawn, deltawt, omega, deltab, deltapar, vtha)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_ranges
      use b2mod_rates
      use b2mod_neoclassical
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2plot_debug
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, jsep, jxa
      real (kind=R8) ::
     *   safetyfactor(-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   nuab(-1:nx,-1:ny,0:ns,0:ns), nua(-1:nx,-1:ny,0:ns),
     &   nustar(-1:nx,-1:ny,0:ns), rhol(-1:nx,-1:ny,0:ns),
     &   lnn(-1:nx,-1:ny,0:ns), ltt(-1:nx,-1:ny,0:ns),
     &   deltal(-1:nx,-1:ny,0:ns), wb(-1:nx,-1:ny,0:ns),
     &   deltaw(-1:nx,-1:ny,0:ns), deltawn(-1:nx,-1:ny,0:ns),
     &   deltawt(-1:nx,-1:ny,0:ns), omega(-1:nx,-1:ny,0:ns),
     &   deltab(-1:nx,-1:ny,0:ns), deltapar(-1:nx,-1:ny,0:ns),
     &   vtha(-1:nx,-1:ny,0:ns)

      real min_value
c     minimum value for real arrays to avoid division by zero
      parameter (min_value=1.e-8)
c     conversion of constants to CGS units
      real erg, cm, gauss, gramm
      parameter (erg=1.e7, cm=1.e2, gauss=1.e4, gramm=1.e3)
      real eV_cgs, e_cgs, c_cgs
      parameter (eV_cgs=ev*erg,
     &           e_cgs=qe*c*10., c_cgs=c*cm)
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
      integer n2, ix, iy, is, is2, ixref, iyref
      save ixref, iyref
      real (kind=R8) ::
     &   gradient, lnab, minlnlt, crmin, crmax
      real (kind=R8) ::
     &   rmaj_(-1:ny)

*   ..procedures
      external ipgeti, ipgetr, xertst, b2xvsg

*   ..subprogram start-up calls
      call subini ('neo_validity')
      if(nnreg(0).eq.4) then
        if (LSN) then
          ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
        else
          ixref=leftcut(1)+(rightcut(1)-leftcut(1))/4
        endif
        iyref=2*topcut(1)/3
      elseif(nnreg(0).eq.5) then
        ixref=leftcut(1)+(rightcut(1)-leftcut(1))/2
        iyref=2*topcut(1)/3
      elseif(nnreg(0).eq.8) then
        ixref=(rightcut(1)+rightcut(2))/2
        iyref=2*min(topcut(1),topcut(2))/3
      elseif(nnreg(0).eq.2) then
        ixref=nx/2
        iyref=2*topcut(1)/3
      else
        ixref=nx/2
        iyref=ny/2
      endif
      call ipgeti ('set_transport_ixref', ixref)
      call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
      call ipgeti ('set_transport_iyref', iyref)
      call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
*   ..test bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*   ..test state
      call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
      call b2xvsg (n2, te, 1, 'te', '.gt.')
      call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
      call b2xvsg (n2, ne, 1, 'ne', '.gt.')

*   ..main

*   ..intialize arrays (set to zero for neutrals)
      vtha(:,:,:)=0._R8
      nuab(:,:,:,:)=0._R8
      nua(:,:,:)=0._R8
      nustar(:,:,:)=0._R8
      rhol(:,:,:)=0._R8
      lnn(:,:,:)=0._R8
      ltt(:,:,:)=0._R8
      deltal(:,:,:)=0._R8
      wb(:,:,:)=0._R8
      deltaw(:,:,:)=0._R8
      deltawn(:,:,:)=0._R8
      deltawt(:,:,:)=0._R8
      omega(:,:,:)=0._R8
      deltab(:,:,:)=0._R8
      deltapar(:,:,:)=0._R8
      rmaj_(:)=0._R8

*   ..calculate major radius of closed flux surfaces
      do iy=-1,jsep
        ix=jxa
        call xertst (on_closed_surface(ix,iy), 'Not in core!')
        crmin=cr(ix,iy)
        crmax=crmin
        ix=rightix(ix,iy)
        do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
          crmin=min(crmin,cr(ix,iy))
          crmax=max(crmax,cr(ix,iy))
          ix=rightix(ix,iy)
        enddo
        call xertst (ix.eq.jxa, 'Core not continuous!')
        rmaj_(iy)=0.5_R8*(crmin+crmax)
      enddo

*   ..calculate gradient lengths lnn and ltt
      do ix=-1,nx
       do iy=0,ny-1
        do is=0,ns-1
         if (abs(na(topix(ix,iy),topiy(ix,iy),is)-na(bottomix(ix,iy),
     &     bottomiy(ix,iy),is)).lt.min_value) then
          if (b2plot_debug.ge.6)
     &     write(*,*) 'ion density gradient too small'
          gradient=min_value/(hy1(ix,iy)+0.5*(hy1(topix(ix,iy),
     &     topiy(ix,iy))+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
         else
          gradient=abs(na(topix(ix,iy),topiy(ix,iy),is)-
     &     na(bottomix(ix,iy),bottomiy(ix,iy),is))/(hy1(ix,iy)+0.5*
     &     (hy1(topix(ix,iy),topiy(ix,iy))+hy1(bottomix(ix,iy),
     &     bottomiy(ix,iy))))
         endif
         if (na(ix,iy,is).lt.min_value) then
          if (b2plot_debug.ge.6) write(*,*) 'ion density too small'
          lnn(ix,iy,is)=min_value/gradient
         else
          lnn(ix,iy,is)=na(ix,iy,is)/gradient
         endif
         if (abs(ti(topix(ix,iy),topiy(ix,iy))-ti(bottomix(ix,iy),
     &     bottomiy(ix,iy))).lt.(min_value*ev)) then
          if (b2plot_debug.ge.6)
     &     write(*,*) 'ion temperature gradient too small'
          gradient=min_value*ev/(hy1(ix,iy)+0.5*(hy1(topix(ix,iy),
     &     topiy(ix,iy))+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
         else
          gradient=abs(ti(topix(ix,iy),topiy(ix,iy))-
     &     ti(bottomix(ix,iy),bottomiy(ix,iy)))/(hy1(ix,iy)+0.5*
     &     (hy1(topix(ix,iy),topiy(ix,iy))+hy1(bottomix(ix,iy),
     &     bottomiy(ix,iy))))
         endif
         if (ti(ix,iy).lt.(min_value*ev)) then
          if (b2plot_debug.ge.6) write(*,*) 'ion temperature too low'
          ltt(ix,iy,is)=min_value*ev/gradient
         else
          ltt(ix,iy,is)=ti(ix,iy)/gradient
         endif
        enddo
        if (abs(ne(topix(ix,iy),topiy(ix,iy))-ne(bottomix(ix,iy),
     &    bottomiy(ix,iy))).lt.min_value) then
         if (b2plot_debug.ge.6)
     &    write(*,*) 'electron density gradient too small'
         gradient=min_value/(hy1(ix,iy)+0.5*(hy1(topix(ix,iy),
     &    topiy(ix,iy))+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
        else
         gradient=abs(ne(topix(ix,iy),topiy(ix,iy))-
     &    ne(bottomix(ix,iy),bottomiy(ix,iy)))/(hy1(ix,iy)+0.5*
     &    (hy1(topix(ix,iy),topiy(ix,iy))+hy1(bottomix(ix,iy),
     &    bottomiy(ix,iy))))
        endif
        if (ne(ix,iy).lt.min_value) then
         if (b2plot_debug.ge.6) write(*,*) 'electron density too small'
         lnn(ix,iy,ns)=min_value/gradient
        else
         lnn(ix,iy,ns)=ne(ix,iy)/gradient
        endif
        if (abs(te(topix(ix,iy),topiy(ix,iy))-te(bottomix(ix,iy),
     &    bottomiy(ix,iy))).lt.(min_value*ev)) then
         if (b2plot_debug.ge.6)
     &    write(*,*) 'electron temperature gradient too small'
         gradient=min_value*ev/(hy1(ix,iy)+0.5*(hy1(topix(ix,iy),
     &    topiy(ix,iy))+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
        else
         gradient=abs(te(topix(ix,iy),topiy(ix,iy))-
     &    te(bottomix(ix,iy),bottomiy(ix,iy)))/(hy1(ix,iy)+0.5*
     &    (hy1(topix(ix,iy),topiy(ix,iy))+hy1(bottomix(ix,iy),
     &    bottomiy(ix,iy))))
        endif
        if (te(ix,iy).lt.(min_value*ev)) then
         if (b2plot_debug.ge.6)
     &    write(*,*) 'electron temperature too low'
         ltt(ix,iy,ns)=min_value*ev/gradient
        else
         ltt(ix,iy,ns)=te(ix,iy)/gradient
        endif
       enddo
       do is=0,ns
        lnn(ix,-1,is)=lnn(ix,0,is)
        ltt(ix,-1,is)=ltt(ix,0,is)
        lnn(ix,ny,is)=lnn(ix,ny-1,is)
        ltt(ix,ny,is)=ltt(ix,ny-1,is)
       enddo
      enddo

*   ..calculate themal velocities
      do is=0,ns-1
       do ix=-1,nx
        do iy=-1,ny
         vtha(ix,iy,is)=sqrt(erg*ti(ix,iy)/(am(is)*mp*gramm))
        enddo
       enddo
      enddo
      do ix=-1,nx
       do iy=-1,ny
        vtha(ix,iy,ns)=sqrt(erg*te(ix,iy)/(me*gramm))
       enddo
      enddo

*   ..calculate collision frequency nuab of species a and b
*      ion-ion collisions
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do is2=0,ns-1
         if (.not. is_neutral(is2)) then
          do ix=-1,nx
           do iy=-1,ny
*      Coulomb logarithm
            lnab=17.3-0.5*log(ne(ix,iy)/1.e20)+1.5*log(ti(ix,iy)*
     &           erg/(eV_cgs*1.e3))
            nuab(ix,iy,is,is2)=16.*sqrt(pi)*na(ix,iy,is2)/1.e6
     &        *rz2(ix,iy,is)*rz2(ix,iy,is2)*e_cgs**4./
     &        (3.*am(is)*mp*gramm)*lnab*
     &        (1./(am(is)*mp*gramm)+1./(am(is2)*mp*gramm))/
     &        (vtha(ix,iy,is)**2.+vtha(ix,iy,is2)**2.)**1.5
           enddo
          enddo
         endif
        enddo
       endif
      enddo
*      ion-electron collisions and electron-ion collisions
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do ix=-1,nx
         do iy=-1,ny
*      Coulomb logarithm
          lnab=15.2-0.5*log(ne(ix,iy)/1.e20)+log(te(ix,iy)*
     &         erg/(eV_cgs*1.e3))
          nuab(ix,iy,is,ns)=16.*sqrt(pi)*ne(ix,iy)/1.e6*rz2(ix,iy,is)
     &      *e_cgs**4./(3.*am(is)*mp*gramm)*lnab*
     &      (1./(am(is)*mp*gramm)+1./(me*gramm))/
     &      (vtha(ix,iy,is)**2.+vtha(ix,iy,ns)**2.)**1.5
          nuab(ix,iy,ns,is)=16.*sqrt(pi)*na(ix,iy,is)/1.e6
     &      *rz2(ix,iy,is)
     &      *e_cgs**4./(3.*me*gramm)*lnab*
     &      (1./(am(is)*mp*gramm)+1./(me*gramm))/
     &      (vtha(ix,iy,is)**2.+vtha(ix,iy,ns)**2.)**1.5
         enddo
        enddo
       endif
      enddo
*      electron-electron collisions
      do ix=-1,nx
       do iy=-1,ny
*      Coulomb logarithm
        lnab=14.9-0.5*log(ne(ix,iy)/1.e20)+log(te(ix,iy)*
     &         erg/(eV_cgs*1.e3))
        nuab(ix,iy,ns,ns)=16.*sqrt(pi)*ne(ix,iy)/1.e6
     &    *e_cgs**4./(3.*me*gramm)*lnab*(2./(me*gramm))/
     &    (2.*vtha(ix,iy,ns)**2)**1.5
       enddo
      enddo

*   ..calculate total collision frequency nua of species a
      do is=0,ns
       do is2=0,ns
        do ix=-1,nx
         do iy=-1,ny
          nua(ix,iy,is)=nua(ix,iy,is)+nuab(ix,iy,is,is2)
         enddo
        enddo
       enddo
      enddo

*   ..calculate collisionality nustar of species a
      do is=0,ns-1
       if (.not. is_neutral(is)) then
        do ix=-1,nx
         do iy=-1,ny
          if (on_closed_surface(ix,iy)) then
           nustar(ix,iy,is)=nua(ix,iy,is)*safetyfactor(iy)*
     &      rmaj_(iy)*cm/(vtha(ix,iy,is)*(1.-fc_(iy))**3)
          else
           nustar(ix,iy,is)=0.0_R8
          endif
         enddo
        enddo
       endif
      enddo
      do ix=-1,nx
       do iy=-1,ny
        if (on_closed_surface(ix,iy)) then
         nustar(ix,iy,ns)=nua(ix,iy,ns)*safetyfactor(iy)*
     &    rmaj_(iy)*cm/(vtha(ix,iy,ns)*(1.-fc_(iy))**3)
        else
         nustar(ix,iy,ns)=0.0_R8
        endif
       enddo
      enddo

*   ..calculate validity tests
      do ix=-1,nx
       do iy=-1,ny
        do is=0,ns-1
         if (.not. is_neutral(is)) then
*      minimum gradient length
          if (lnn(ix,iy,is).lt.ltt(ix,iy,is)) then
           minlnlt=lnn(ix,iy,is)
          else
           minlnlt=ltt(ix,iy,is)
          endif
*      Larmor radius
          rhol(ix,iy,is)=am(is)*mp*gramm*c_cgs*vtha(ix,iy,is)/
     &      (rza(ix,iy,is)*e_cgs*bb(ix,iy,3)*gauss*cm)
*      ratio of Larmor radius to minimum gradient length
          deltal(ix,iy,is)=rhol(ix,iy,is)/minlnlt
*      width of Banana orbit
          wb(ix,iy,is)=am(is)*mp*gramm*c_cgs*(1.-fc_(iy))*
     &      vtha(ix,iy,is)/(rza(ix,iy,is)*e_cgs*bb(ix,iy,0)*gauss*cm)
*      ratio of Banana orbit to minimum gradient length
          deltaw(ix,iy,is)=wb(ix,iy,is)/minlnlt
*      ratio of Banana orbit to density gradient length
          deltawn(ix,iy,is)=wb(ix,iy,is)/lnn(ix,iy,is)
*      ratio of Banana orbit to temperature gradient length
          deltawt(ix,iy,is)=wb(ix,iy,is)/ltt(ix,iy,is)
*      ion gyro-frequency
          omega(ix,iy,is)=rza(ix,iy,is)*e_cgs*bb(ix,iy,3)*gauss/
     &      (am(is)*mp*gramm*c_cgs)
*      ratio of the total collision frequency to the gyro-frequency
          deltab(ix,iy,is)=nua(ix,iy,is)/omega(ix,iy,is)
*      measure for the density perturbation in a collision dominated
*      regime
          deltapar(ix,iy,is)=(1.-fc_(iy))**4.*nustar(ix,iy,is)*
     &      wb(ix,iy,is)*rmaj_(iy)*nua(ix,iy,is)/(minlnlt**2.*
     &      omega(ix,iy,is))
         endif
        enddo
*      minimum gradient length
        if (lnn(ix,iy,ns).lt.ltt(ix,iy,ns)) then
         minlnlt=lnn(ix,iy,ns)
        else
         minlnlt=ltt(ix,iy,ns)
        endif
*      Larmor radius
        rhol(ix,iy,ns)=me*gramm*c_cgs*vtha(ix,iy,ns)/
     &    (e_cgs*bb(ix,iy,3)*gauss*cm)
*      ratio of Larmor radius to minimum gradient length
        deltal(ix,iy,ns)=rhol(ix,iy,ns)/minlnlt
*      width of Banana orbit
        wb(ix,iy,ns)=me*gramm*c_cgs*(1.-fc_(iy))*
     &    vtha(ix,iy,ns)/(e_cgs*bb(ix,iy,0)*gauss*cm)
*      ratio of Banana orbit to minimum gradient length
        deltaw(ix,iy,ns)=wb(ix,iy,ns)/minlnlt
*      ratio of Banana orbit to density gradient length
        deltawn(ix,iy,ns)=wb(ix,iy,ns)/lnn(ix,iy,ns)
*      ratio of Banana orbit to temperature gradient length
        deltawt(ix,iy,ns)=wb(ix,iy,ns)/ltt(ix,iy,ns)
*      electron gyro-frequency
        omega(ix,iy,ns)=e_cgs*bb(ix,iy,3)*gauss/(me*gramm*c_cgs)
*      ratio of the total collision frequency to the gyro-frequency
        deltab(ix,iy,ns)=nua(ix,iy,ns)/omega(ix,iy,ns)
*      measure for the density perturbation in a collision dominated
*      regime
        deltapar(ix,iy,ns)=(1.-fc_(iy))**4.*nustar(ix,iy,ns)*
     &    wb(ix,iy,ns)*rmaj_(iy)*nua(ix,iy,ns)/(minlnlt**2.*
     &    omega(ix,iy,ns))
       enddo
      enddo


*   ..return
      call subend ()
      return
      end subroutine neo_validity

      subroutine chang_hinton(nx,ny,ns,jsep,jxa,lambda,safetyfactor,
     1 cr,cz,crx,bb,ne,te,ti,na,transport)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_b2plot_debug
      use b2mod_external
      implicit none
      integer nx,ny,ns,jsep,jxa
      real (kind=R8) :: lambda,
     1 cr(-1:nx,-1:ny), cz(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     2 bb(-1:nx,-1:ny,0:3),
     2 ne(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     3 na(-1:nx,-1:ny,0:ns-1), safetyfactor(-1:ny),
     2 transport(-1:nx,-1:ny,0:ns-1)

      integer ix,iy,is,irmin,irmax, is_fuel
      real (kind=R8) :: crmin,crmax,czmin,czmax,
     1 cr2max,cr2min,cr1max,cr1min,
     2 bave,epsilon,lnlam,vthi,omegaci,
     1 rhoi,alpha,deltap,f1,f2,taui,nui,nustari,mustari,chi,
     2 bracket1, bracket2, ammin
      external subini, subend, xertst

      call subini ('chang_hinton')
      transport=0.0_R8
      do iy=-1,jsep
        if (b2plot_debug.ge.8) write(*,*) 'iy: ',iy
        ix=jxa
        call xertst (on_closed_surface(ix,iy), 'Not in core')
        irmin=jxa
        irmax=jxa
        crmin=cr(ix,iy)
        crmax=crmin
        czmin=cz(ix,iy)
        czmax=czmin
        ix=rightix(ix,iy)
        do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
          if(cr(ix,iy).lt.crmin) irmin=ix
          if(cr(ix,iy).gt.crmax) irmax=ix
          crmin=min(crmin,cr(ix,iy))
          crmax=max(crmax,cr(ix,iy))
          czmin=min(czmin,cz(ix,iy))
          czmax=max(czmax,cz(ix,iy))
          ix=rightix(ix,iy)
        enddo
        call xertst (ix.eq.jxa, 'Core not continuous')
        if (b2plot_debug.ge.8)
     >   write(*,*) 'crmin,crmax,czmin,czmax',
     1   crmin,crmax,czmin,czmax
        bave=cr(ix,iy)*2.0_R8/(crmax+crmin)*bb(ix,iy,3)
        if (b2plot_debug.ge.8) write(*,*) 'ne,te,ti,b,q: ',
     1   ne(ix,iy),te(ix,iy)/qe,ti(ix,iy)/qe,bave,safetyfactor(iy)
        epsilon=(crmax-crmin)/(crmax+crmin)
        if (b2plot_debug.ge.8) write(*,*) 'epsilon: ',epsilon
        if (lambda.lt.0.0_R8) then
          lnlam=max(-lambda,17.3_R8-
     1     log(ne(ix,iy)/1e20_R8)/2.0_R8+
     2     1.5_R8*log(ti(ix,iy)/ev/1000.0_R8))
        else
          lnlam = lambda
        endif
        if (b2plot_debug.ge.8) write(*,*) 'lnlam: ',lnlam
        cr2max=0.5_R8*(crx(irmax,iy,1)+crx(irmax,iy,3))
        cr2min=0.5_R8*(crx(irmin,iy,1)+crx(irmin,iy,3))
        cr1max=0.5_R8*(crx(irmax,iy,0)+crx(irmax,iy,2))
        cr1min=0.5_R8*(crx(irmin,iy,0)+crx(irmin,iy,2))
        if (b2plot_debug.ge.8)
     >   write(*,*) 'cr2max,cr2min,cr1max,cr1min: ',
     1   cr2max,cr2min,cr1max,cr1min
        deltap=((cr2max+cr2min)-(cr1max+cr1min))/
     1   ((cr2max-cr2min)-(cr1max-cr1min))
        deltap=-2._R8
        if (b2plot_debug.ge.8) write(*,*) 'deltap: ',deltap
        f1=(1_R8+1.5_R8*(epsilon**2+epsilon*deltap)+
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
            alpha=alpha+na(ix,iy,is)*rz2(ix,iy,is)
          endif
        enddo
        do is=0,ns_ext-1
          if (.not. is_ext_neutral(is)) then
            alpha=alpha+na_ext(ix,iy,is)*za2_ext(ix,iy,is)
          endif
        enddo
        alpha=alpha/(na(ix,iy,is_fuel)*rz2(ix,iy,is_fuel))
        if (b2plot_debug.ge.8) write(*,*) 'alpha: ', alpha
        do is=0,ns-1
          if(.not. is_neutral(is)) then
            if (b2plot_debug.ge.8) write(*,*) 'rza,am,na: ',
     >       rza(ix,iy,is),am(is),na(ix,iy,is)
            vthi=sqrt(ti(ix,iy)/mp/am(is))
            if (b2plot_debug.ge.8) write(*,*) 'vthi: ',vthi
            omegaci=rza(ix,iy,is)*qe*bave/mp/am(is)
            if (b2plot_debug.ge.8) write(*,*) 'omegaci: ',omegaci
            rhoi=sqrt(2.0_R8)*vthi/omegaci
            if (b2plot_debug.ge.8) write(*,*) 'rhoi: ',rhoi
            taui=12.0_R8*pi**1.5_R8*eps0**2*
     1       (am(is)*mp)**0.5_R8*ti(ix,iy)**1.5_R8/
     2       (na(ix,iy,is)*rz2(ix,iy,is)**2*qe**4*lnlam)
            if (b2plot_debug.ge.8) write(*,*) 'taui: ',taui
            nui=1.0_R8/taui
            nustari=nui*(crmax+crmin)/2.0_R8*safetyfactor(iy)/
     1       epsilon**(1.5_R8)/vthi
            if (b2plot_debug.ge.8) write(*,*) 'nustari: ',nustari
            mustari=nustari*(1.0_R8+1.54_R8*alpha)
            if (b2plot_debug.ge.8) write(*,*) 'mustari: ',mustari
            bracket1=(0.66_R8*(1.0_R8+1.54_R8*alpha)+
     3               (1.88_R8*epsilon**0.5_R8-1.54_R8*epsilon)*
     4               (1_R8+3.75_R8*alpha))/
     5               (1.0_R8+1.03_R8*mustari**0.5_R8+0.31*mustari)
            bracket2=0.59_R8*mustari*epsilon/
     6         (1_R8+0.74_R8*mustari*epsilon**1.5_R8)*
     7         (1_R8+
     8           1.33_R8*alpha*(1.0_R8+0.60_R8*alpha)/
     9           (1.0_R8+1.79_R8*alpha))
            if (b2plot_debug.ge.8) write(*,*) 'bracket1 =', bracket1
            if (b2plot_debug.ge.8) write(*,*) 'bracket2 =', bracket2
            chi=epsilon**(-1.5_R8)*
     1       safetyfactor(iy)**2*rhoi**2/taui*
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
            transport(ix,iy,is)=chi
            ix=rightix(ix,iy)
            do while(ix.ne.jxa.and.ix.ge.-1.and.ix.le.nx)
              transport(ix,iy,is)=chi
              ix=rightix(ix,iy)
            enddo
          endif
        enddo
      enddo
      call subend ()
      return
      end subroutine chang_hinton















