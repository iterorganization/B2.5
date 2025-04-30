*-----------------------------------------------------------------------
*.specification

      subroutine neo_diagnostics(mpg, geo, state, ns, ic,
     &  dn_neo, chi_neo, vn_neo, vw_neo)
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
*   ..set the contribution ic
*     0 --- classical particle flux
*     1 --- banana plateau contribution
*     2 --- Pfirsch Schlueter contribution
*     3 --- both banana and PS
*     4 --- all contributions
*     mind that B2 already calculates the classical transport !
*     avoid double transport, 0+4 for cross-checks only !
      integer ic
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   dn_neo(mpg%nCv,0:ns), chi_neo(mpg%nCv,0:ns),
     &   vn_neo(mpg%nCv,0:ns), vw_neo(mpg%nCv,0:ns)

*     nints=number of intervals for integration
      integer nints
      parameter (nints=10)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     NEO_DIAGNOSTICS computes the neoclassical transport coefficients
*
*
*  3. description
*
*     This routine computes the following neoclassical transport
*     coefficients
*
*     - dn_neo:  particle transport coefficients for grad(na).
*     - chi_neo: heat conductivity
*     - vn_neo:  particle pinch velocity
*     - vw_neo:  heat pinch velocity
*-----------------------------------------------------------------------
*   ..local variables
      integer iCv, iCm, iCn, iFc, iFt, iFm, iFp, iVx, iVx1, iVx2
      integer i, j, is, is2
*   ..local variables for neoclassical transport
      integer :: istate, ispecies
      integer :: isorth(0:ns,0:nzdecl-1), nclh(0:ns), nslh
      real (kind=R8) ::
     *   zah(0:ns-1), amh(0:ns-1), t(0:ns),
     &   den(0:ns,0:nzdecl-1), ds(0:ns,0:nzdecl-1,2),
     &   denl(0:ns,0:nzdecl-1), tl(0:ns),
     &   denr(0:ns,0:nzdecl-1), tr(0:ns),
     &   vtn(0:ns,0:nzdecl-1,mpg%nFt),
     &   kappa(0:ns,0:nzdecl-1,mpg%nFt),
     &   vw(0:ns,0:nzdecl-1,mpg%nFt),
     &   dfhcb, dfhce, hthet(mpg%nCv), quant(mpg%nCv),
     &   flux_width(mpg%nFt)
      real (kind=R8) ::
     *   zacomp, amcomp, epsn
*   ..set the accuracy
      parameter (epsn=1.e-3)
      real (kind=R8) ::
     *   thet, ilpb, ilp(mpg%nFt), bmax, nlambda, nfactor,
     &   eparr, rho, rhomax, hy1
      save eparr
      real (kind=R8) ::
     *   fluxav, sfluxav, tfluxav
*   ..for consistency check
      real (kind=R8) ::
     *   densav(0:ns,0:nzdecl-1,mpg%nFt),
     &   gradnsav(0:ns,0:nzdecl-1,mpg%nFt),
     &   tforce(0:ns,0:nzdecl-1,1:2,mpg%nFt)

*   ..procedures
      external xertst, b2xvsg, fluxav, sfluxav,
     &         tfluxav, neodv, neochi

*   ..subprogram start-up calls
      call subini ('neo_diagnostics')
*   ..neoclassical transport only valid for single charge states
*     bundled charge states yet to be implemented
      write(6,*) 'Allocated neo_diagnostics'
      write(*,*) 'NEOART option ', ic, ' has been used'
      if(ic.eq.0)
     1 write(*,*) '  NEOART: classical'
      if(ic.eq.1)
     1 write(*,*) '  NEOART: banana plateau'
      if(ic.eq.2)
     1 write(*,*) '  NEOART: Pfirsch-Schlueter'
      if(ic.eq.3)
     1 write(*,*) '  NEOART: banana plateau + Pfirsch-Schlueter'
      if(ic.eq.4)
     1 write(*,*) '  NEOART: all'
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

*   ..initialize densities and temperatures
      do ispecies=0,ns
        t(ispecies)=0.0_R8
        do istate=0,nzdecl-1
          den(ispecies,istate)=0.0_R8
          zsp(ispecies,istate)=0.0_R8
        enddo
      enddo

*   ..species reordering and calculation of geometric parameters

*    ..loop voltage eparr (to be included in a separate routine!)
*     default value:
      eparr=0.0_R8
*     read scalar from b2mn.dat
      call ipgetr ('b2tqna_eparr', eparr)
      if (eparr.ne.0.0_R8 .or. b2plot_debug.ge.3)
     > write(*,*) 'Loop voltage: eparr =', eparr

*   ..initialize geometric parameters
      do iFt=1,mpg%nFt
        bav_(iFt)=0.0_R8
        b2av_(iFt)=0.0_R8
        bi2a_(iFt)=0.0_R8
        rbt_(iFt)=0.0_R8
        bgradp_(iFt)=0.0_R8
        dpsidr_(iFt)=0.0_R8
        rnq_(iFt)=0.0_R8
        fc_(iFt)=0.0_R8
        gclass_(iFt)=0.0_R8
        do is=1,mmx
          fm_(is,iFt)=0.0_R8
        enddo
        r2i_(iFt)=0.0_R8
        rho_(iFt)=0.0_R8
      enddo

*   ..split species into species and charge states
      do ispecies=0,ns
        ncl(ispecies)=0
        do istate=0,nzdecl-1
          isort(ispecies,istate)=0
        enddo
      enddo
      do is=0,ns-1
        zah(is)=(zamin(is)+zamax(is))/2.0_R8
        if (zamin(is).ne.zamax(is)) then
          write(6,*) 'Warning: '//
     .     'neoclassical diagnostics may not be accurate '//
     .     'for bundled charge states!'
        end if
        amh(is)=am(is)
      enddo
      do is=0,ns-2
        zacomp=zah(is)
        amcomp=amh(is)
        ispecies=is
        do is2=is+1,ns-1
          if((abs(amh(is2)-amcomp).lt.0.5_R8).and.
     &       (zah(is2).lt.zacomp)) then
            zacomp=zah(is2)
            ispecies=is2
          elseif(amh(is2).lt.amcomp) then
            amcomp=amh(is2)
            zacomp=zah(is2)
            ispecies=is2
          endif
        enddo
        zah(ispecies)=zah(is)
        amh(ispecies)=amh(is)
        zah(is)=zacomp
        amh(is)=amcomp
      enddo
      ispecies=1
      istate=1
      amcomp=amh(0)
      do is=1,ns-1
        if((amh(is)-amcomp).gt.0.5_R8) then
          ncl(ispecies-1)=istate
          ispecies=ispecies+1
          amcomp=amh(is)
          istate=1
        elseif(abs(amh(is)-amcomp).lt.0.5_R8) then
          istate=istate+1
        endif
      enddo

      ncl(ispecies-1)=istate

      nsl=ispecies

      is=0
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          do is2=0,ns-1
            if((abs(am(is2)-amh(is)).lt.0.5_R8).and.
     &         (abs((zamin(is2)+zamax(is2))/2.0_R8-zah(is)).lt.
     .          0.5_R8)) then
              isort(ispecies,istate)=is2
            endif
          enddo
          is=is+1
        enddo
      enddo

*   ..remove neutrals from list
      do ispecies=0,nsl-1
        nclh(ispecies)=ncl(ispecies)
      enddo
      nslh=nsl
      isorth=isort
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          if (is_neutral(isort(ispecies,istate))) then
            nclh(ispecies)=nclh(ispecies)-1
            if (nclh(ispecies).eq.0) nslh=nslh-1
          endif
        enddo
      enddo
      is=0
      do ispecies=0,nsl-1
        if (nclh(ispecies).ne.0) then
           nclh(is)=nclh(ispecies)
          is=is+1
        endif
      enddo
      is=0
      is2=0
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          if (.not. is_neutral(isorth(ispecies,istate))) then
            isort(is,is2)=isorth(ispecies,istate)
            is2=is2+1
            if (is2.gt.(nclh(is)-1)) then
              is=is+1
              is2=0
            endif
          endif
        enddo
      enddo

      do ispecies=0,nsl-1
        if (ispecies.lt.nslh) then
          ncl(ispecies)=nclh(ispecies)
        else
          ncl(ispecies)=0
        endif
      enddo
      nsl=nslh

*   ..mass
      do ispecies=0,nsl-1
        m(ispecies)=am(isort(ispecies,0))*mp
      enddo

*   ..add electrons as a separate species
      nsl=nsl+1
      ncl(nsl-1)=1
      m(nsl-1)=me
      zsp(nsl-1,0)=-1.0_R8

*   ..calculate generalized poloidal angle theta
*      default value 0 for grid points outside the core and on
*      the inner boundary
      theta=0.0_R8

*    ..determine poloidal angle
      do iFt=1,mpg%nFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if(mpg%cvOnClosedSurface(iCv)) then
          iCm=iCv
          thet=0.5_R8*geo%cvBb(iCv,3)/geo%cvBb(iCv,0)*geo%cvHx(iCv)
          ilp(iFt)=geo%cvBb(iCv,0)*geo%cvHx(iCv)
          do i=2,mpg%ftCvP(iFt,2)
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            iCn=mpg%ftCv(mpg%ftCvP(iFt,1)+i-2)
            if(geo%cvX(iCv).gt.geo%cvX(iCm)) iCm=iCv
            thet=thet
     .       +0.5_R8*geo%cvBb(iCn,3)/geo%cvBb(iCn,0)*geo%cvHx(iCn)
     .       +0.5_R8*geo%cvBb(iCv,3)/geo%cvBb(iCv,0)*geo%cvHx(iCv)
            theta(iCv)=thet
            ilp(iFt)=ilp(iFt)+geo%cvBb(iCv,0)*geo%cvHx(iCv)
          enddo
          thet=thet
     .       +0.5_R8*geo%cvBb(iCv,3)/geo%cvBb(iCv,0)*geo%cvHx(iCv)
          ilpb=thet
          do i=1,mpg%ftCvP(iFt,2)
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            theta(iCv)=theta(iCv)-theta(iCm)
            theta(iCv)=2.0_R8*pi*mod(theta(iCv)/ilpb,1.0_R8)
            if (theta(iCv).lt.0.0_R8) theta(iCv)=theta(iCv)+2.0_R8*pi
          enddo
          is=0
          flux_width(iFt)=0.0_R8
          do i=1,mpg%cvFcP(iCm,2)
            iFc=mpg%cvFc(mpg%cvFcP(iCm,1)+i-1)
            iVx1=mpg%fcVx(iFc,1)
            iVx2=mpg%fcVx(iFc,2)
            if (mpg%vxFs(iVx1).ne.0.and.mpg%vxFs(iVx2).ne.0.and.
     &          mpg%vxFs(iVx1).ne.mpg%vxFs(iVx2)) then
              is=is+1
              flux_width(iFt)=flux_width(iFt)+
     &          geo%fcHt(iFc)*geo%fcQgam(iFc,0)
            end if
          end do
          if(is.gt.0) flux_width(iFt)=flux_width(iFt)/real(is,R8)
        endif
      enddo

*   ..calculate geometric parameters

      rhomax=0.0_R8
*   ..loop over flux tubes
      do iFt=1,mpg%nFt

*   ..calculate neoclassical transport only for closed flux tubes

        iCv = mpg%ftCv(mpg%ftCvP(iFt,1))
        if (mpg%cvOnClosedSurface(iCv)) then
          iCm=iCv
          bmax=geo%cvBb(iCv,3)
          do i=2,mpg%ftCvP(iFt,2)
            iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
            if(geo%cvX(iCv).gt.geo%cvX(iCm)) iCm=iCv
            bmax=max(bmax,geo%cvBb(iCv,3))
          end do

*   ..flux averages

          quant(:)=geo%cvBb(:,3)
          bav_(iFt)=fluxav(quant,mpg,geo,iFt)
          quant(:)=geo%cvBb(:,3)**2
          b2av_(iFt)=fluxav(quant,mpg,geo,iFt)
          quant(:)=1.0_R8/geo%cvBb(:,3)**2
          bi2a_(iFt)=fluxav(quant,mpg,geo,iFt)
          quant(:)=geo%cvBb(:,2)*geo%cvX(:)
          rbt_(iFt)=fluxav(quant,mpg,geo,iFt)
          bgradp_(iFt)=2.0_R8*pi/(bav_(iFt)*ilp(iFt))
          dpsidr_(iFt)=geo%cvBb(iCm,0)*geo%cvX(iCm)
          rnq_(iFt)=geo%ftConn(iFt)/(2.0_R8*pi)

          nlambda=1.0_R8/nints

          fc_(iFt)=0.0_R8
          do is2=0,nints
            if ((is2.eq.0).or.(is2.eq.nints)) then
              nfactor=0.5_R8
            else
              nfactor=1.0_R8
            endif
            quant(:)=geo%cvBb(:,3)
            fc_(iFt)=fc_(iFt)+is2*nlambda/
     &       sfluxav(quant,mpg,geo,iFt,is2*nlambda,bmax)*nlambda*nfactor
          enddo
          fc_(iFt)=0.75_R8*b2av_(iFt)/bmax**2*fc_(iFt)
          quant(:)=(geo%cvX(:)*geo%cvBb(:,0)/geo%cvBb(:,3))**2
          gclass_(iFt)=fluxav(quant,mpg,geo,iFt)

          do is2=1,mmx
            hthet(:)=is2*theta(:)
            quant(:)=cos(hthet(:))*geo%cvBb(:,3)*log(geo%cvBb(:,3))
            cosbl(is2)=fluxav(quant,mpg,geo,iFt)
            quant(:)=cos(hthet(:))*geo%cvBb(:,3)**2
            cosb2(is2)=fluxav(quant,mpg,geo,iFt)
            quant(:)=sin(hthet(:))*geo%cvBb(:,3)*log(geo%cvBb(:,3))
            sinbl(is2)=fluxav(quant,mpg,geo,iFt)
            quant(:)=sin(hthet(:))*geo%cvBb(:,3)**2
            sinb2(is2)=fluxav(quant,mpg,geo,iFt)
          enddo
          do is2=1,mmx
            fm_(is2,iFt)=2.0_R8*(is2*bgradp_(iFt))**2/
     &                (b2av_(iFt)*bav_(iFt))*
     &                (cosb2(is2)*cosbl(is2)+sinb2(is2)*sinbl(is2))
          enddo
          quant(:)=(1.0_R8/geo%cvX(:))**2
          r2i_(iFt)=fluxav(quant,mpg,geo,iFt)

*    ..flux tube label
          rho_(iFt)=geo%cvX(iCm)

          do i=1,mpg%cvVxP(iCm,2)
            rhomax=max(rhomax,geo%vxX(mpg%cvVx(mpg%cvVxP(iCm,1)+i-1)))
          end do

        endif

      enddo

*   ..now calculate the neoclassical transport

*   ..loop over flux tubes
       do iFt=1,mpg%nFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if (iCv.gt.mpg%nCi) cycle

*   ..calculate neoclassical transport only for closed flux tubes
        if (mpg%cvOnClosedSurface(iCv)) then
          iFm=0
          iFp=0
          do i=1,mpg%cvVxP(iCv,2)
            iVx=mpg%cvVx(mpg%cvVxP(iCv,1)+i-1)
            do j=1,mpg%vxCvP(iVx,2)
              iCn=mpg%vxCv(mpg%vxCvP(iVx,1)+j-1)
              if (mpg%cvFt(iCn).eq.iFt) cycle
              if (mpg%cvFt(iCn).eq.iFm) cycle
              if (mpg%cvFt(iCn).eq.iFp) cycle
              if (iFm.eq.0) then
                if (rho_(mpg%cvFt(iCn)).lt.rho_(iFt)) iFm=mpg%cvFt(iCn)
              end if
              if (iFp.eq.0) then
                if (rho_(mpg%cvFt(iCn)).ne.0.0_R8) then
                  if (rho_(mpg%cvFt(iCn)).gt.rho_(iFt))
     >             iFp=mpg%cvFt(iCn)
                else if (mod(mpg%cvReg(iCn),4).eq.2) then
                  iFp=mpg%cvFt(iCn)
                end if
              end if
            end do
          end do

*   ..flux averages
*     temperatures in keV, densities in 10^19 m^-3
          do ispecies=0,nsl-2
            t(ispecies)=fluxav(state%pl%ti,mpg,geo,iFt)/
     &                  (1.0e3_R8*ev)
            tl(ispecies)=fluxav(state%pl%ti,mpg,geo,iFm)/
     &                  (1.0e3_R8*ev)
            tr(ispecies)=tfluxav(state%pl%ti,mpg,geo,iFp)/
     &                  (1.0e3_R8*ev)
            do istate=0,ncl(ispecies)-1
              quant(:)=state%pl%na(:,isort(ispecies,istate))
              den(ispecies,istate)=
     &          fluxav(quant,mpg,geo,iFt)/1.0e19_R8
              denl(ispecies,istate)=
     &          fluxav(quant,mpg,geo,iFm)/
     &          1.0e19_R8
              denr(ispecies,istate)=
     &          tfluxav(quant,mpg,geo,iFp)/1.0e19_R8
              quant(:)=state%rt%rza(:,isort(ispecies,istate))
              zsp(ispecies,istate)=fluxav(quant,mpg,geo,iFt)
            enddo
          enddo
          t(nsl-1)=fluxav(state%pl%te,mpg,geo,iFt)/(1.0e3_R8*ev)
          den(nsl-1,0)=fluxav(state%dv%ne,mpg,geo,iFt)/1.0e19_R8
          tl(nsl-1)=fluxav(state%pl%te,mpg,geo,iFm)/(1.e3_R8*ev)
          denl(nsl-1,0)=fluxav(state%dv%ne,mpg,geo,iFm)/1.e19_R8
          tr(nsl-1)=tfluxav(state%pl%te,mpg,geo,iFp)/(1.0e3_R8*ev)
          denr(nsl-1,0)=tfluxav(state%dv%ne,mpg,geo,iFp)/1.0e19_R8

*   ..thermodynamic forces
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              if (mpg%cvOnClosedSurface(mpg%ftCv(mpg%ftCvP(iFp,1))))
     >         then
                ds(ispecies,istate,1)=(log(denr(ispecies,istate)*
     &           tr(ispecies))-log(denl(ispecies,istate)*
     &           tl(ispecies)))/(rho_(iFp)-rho_(iFm))
                ds(ispecies,istate,2)=(log(tr(ispecies))-
     &           log(tl(ispecies)))/(rho_(iFp)-rho_(iFm))
              else
                ds(ispecies,istate,1)=(log(denr(ispecies,istate)*
     &           tr(ispecies))-log(denl(ispecies,istate)*
     &           tl(ispecies)))/(rhomax-rho_(iFm))
                ds(ispecies,istate,2)=(log(tr(ispecies))-
     &           log(tl(ispecies)))/(rhomax-rho_(iFm))
              endif

*   ..for consistency check
              densav(ispecies,istate,iFt)=den(ispecies,istate)
              gradnsav(ispecies,istate,iFt)=den(ispecies,istate)*
     &            (ds(ispecies,istate,1)-ds(ispecies,istate,2))
              tforce(ispecies,istate,1,iFt)=ds(ispecies,istate,1)
              tforce(ispecies,istate,2,iFt)=ds(ispecies,istate,2)
            enddo
          enddo

*   ..recover geometric parameters for the flux tube iFt
          bav=bav_(iFt)
          b2av=b2av_(iFt)
          bi2a=bi2a_(iFt)
          rbt=rbt_(iFt)
          bgradp=bgradp_(iFt)
          dpsidr=dpsidr_(iFt)
          rnq=rnq_(iFt)
          fc=fc_(iFt)
          gclass=gclass_(iFt)
          do is2=1, mmx
            fm(is2)=fm_(is2,iFt)
          enddo
          r2i=r2i_(iFt)
          rho=rho_(iFt)

*    ..calculate the diffusion coefficients
          call neodv(ns,nzdecl,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &               epsn,dd(:,:,iFt),vvn(:,:,iFt),vtn(:,:,iFt))

*    ..calculate the heat conduction coefficients
          call neochi(ns,nzdecl,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &                epsn,kappa(:,:,iFt),chi(:,:,iFt),vw(:,:,iFt))

*    ..calculate thermal conductivity and heat pinch from
*        q(i,j) = - kappa(i,j)*grad.T(i) + vw(i,j)*T(i)
*               = - kappa(i,j)*grad.T(i) + den(i,j)*vwt(i,j)*T(i)
          do ispecies=0, nsl-1
            do istate=0, ncl(ispecies)-1
              vwt(ispecies,istate,iFt)=vw(ispecies,istate,iFt)/
     &                                (den(ispecies,istate)*1e19)
           enddo
          enddo

        else
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              dd(ispecies,istate,iFt)=0.0_R8
              vvn(ispecies,istate,iFt)=0.0_R8
              chi(ispecies,istate,iFt)=0.0_R8
              kappa(ispecies,istate,iFt)=0.0_R8
              vwt(ispecies,istate,iFt)=0.0_R8
              densav(ispecies,istate,iFt)=0.0_R8
              gradnsav(ispecies,istate,iFt)=0.0_R8
              tforce(ispecies,istate,1,iFt)=0.0_R8
              tforce(ispecies,istate,2,iFt)=0.0_R8
            enddo
          enddo
        endif

       enddo

*    ..set transport coefficients on boundary by symmetric condition
       do iFt=1,mpg%nFt
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
        if (iCv.le.mpg%nCi) cycle
        iFc=mpg%cvFc(mpg%cvFcP(iCv,1))
        if (mpg%fcCv(iFc,1).eq.iCv) then
          iCn=mpg%fcCv(iFc,2)
        else
          iCn=mpg%fcCv(iFc,1)
        endif
        if(.not.mpg%cvOnClosedSurface(iCn)) cycle
        iFp=mpg%cvFt(iCn)

        do ispecies=0,nsl-1
         do istate=0,ncl(ispecies)-1
           densav(ispecies,istate,iFt)=densav(ispecies,istate,iFp)
           gradnsav(ispecies,istate,iFt)=gradnsav(ispecies,istate,iFp)
           tforce(ispecies,istate,1,iFt)=tforce(ispecies,istate,1,iFp)
           tforce(ispecies,istate,2,iFt)=tforce(ispecies,istate,2,iFp)
           dd(ispecies,istate,iFt)=dd(ispecies,istate,iFp)
           vvn(ispecies,istate,iFt)=vvn(ispecies,istate,iFp)
           chi(ispecies,istate,iFt)=chi(ispecies,istate,iFp)
           kappa(ispecies,istate,iFt)=kappa(ispecies,istate,iFp)
           vwt(ispecies,istate,iFt)=vwt(ispecies,istate,iFp)
         enddo
        enddo
       enddo

*    ..determine the neoclassical contribution to the transport coefficients
*       assuming the following structure for the particle and heat flux
*         gamma(i,j) = - dd(i,j)*grad.n(i,j) + vvn(i,j)*den(i,j)
*           already including the ware pinch, and
*         q(i,j) = - kappa(i,j)*grad.T(i) + den(i,j)*vwt(i,j)*T(i)
*       also including the ware pinch
*       coefficients: dn_neo, chi_neo, vn_neo, and vw_neo

*    ..initialize transport coefficients
      dn_neo(:,:)=0.0_R8
      chi_neo(:,:)=0.0_R8
      vn_neo(:,:)=0.0_R8
      vw_neo(:,:)=0.0_R8

      do iCv=1,mpg%nCv
        if (mpg%cvOnClosedSurface(iCv)) then
         iFt=mpg%cvFt(iCv)
         is=0
         hy1=0.0_R8
         do i=1,mpg%cvFcP(iCv,2)
           iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
           iVx1=mpg%fcVx(iFc,1)
           iVx2=mpg%fcVx(iFc,2)
           if (mpg%vxFs(iVx1).ne.0.and.mpg%vxFs(iVx2).ne.0.and.
     &         mpg%vxFs(iVx1).ne.mpg%vxFs(iVx2)) then
             is=is+1
             hy1=hy1+geo%fcHt(iFc)*geo%fcQgam(iFc,0)
          end if
         end do
         if(is.gt.0) hy1=hy1/real(is,R8)
         do ispecies=0,nsl-2
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dn_neo(iCv,is)=dd(ispecies,istate,iFt)*
     &                      (hy1/flux_width(iFt))**2
           vn_neo(iCv,is)=vvn(ispecies,istate,iFt)*
     &                      (hy1/flux_width(iFt))
           dfhcb=chi(ispecies,istate,iFt)*
     &                      (hy1/flux_width(iFt))**2
           chi_neo(iCv,is)=dfhcb
           vw_neo(iCv,is)=vwt(ispecies,istate,iFt)*
     &                      (hy1/flux_width(iFt))
          enddo
         enddo
         dfhce=chi(nsl-1,0,iFt)*(hy1/flux_width(iFt))**2
         if (dfhce.lt.0.0_R8) then
           write(*,*) 'dfhce(',iCv,')=',dfhce
         endif
         dn_neo(iCv,ns)=dd(nsl-1,0,iFt)*
     &                    (hy1/flux_width(iFt))**2
         vn_neo(iCv,ns)=vvn(nsl-1,0,iFt)*
     &                    (hy1/flux_width(iFt))
         chi_neo(iCv,ns)=dfhce
         vw_neo(iCv,ns)=vwt(nsl-1,0,iFt)*(hy1/flux_width(iFt))

*   ..set transport to zero in the SOL and divertor regions
        elseif (.not.mpg%cvOnClosedSurface(iCv)) then
         do ispecies=0,nsl-1
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dn_neo(iCv,is)=0.0_R8
           vn_neo(iCv,is)=0.0_R8
           chi_neo(iCv,is)=0.0_R8
           vw_neo(iCv,is)=0.0_R8
          enddo
         enddo
        endif

      enddo

*   ..return
      call subend ()
      return
      end subroutine neo_diagnostics
