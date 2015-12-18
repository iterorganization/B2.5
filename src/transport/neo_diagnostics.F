*-----------------------------------------------------------------------
*.specification

      subroutine neo_diagnostics(nx, ny, ns, nstages,            ! cjk
     &  dn_neo, chi_neo, vn_neo, vw_neo)
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
      integer nx, ny, ns, nstages
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   dn_neo(-1:nx,-1:ny,0:ns), chi_neo(-1:nx,-1:ny,0:ns),
     &   vn_neo(-1:nx,-1:ny,0:ns), vw_neo(-1:nx,-1:ny,0:ns)

      integer ::
     *  nar,nzm
*     nints=number of intervals for integration
      integer nints
      parameter (nints=10)
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
      integer n2, ix, iy, is, is2, ixref, iyref
      save ixref, iyref
*   ..local variables for neoclassical transport
      integer ic
      integer :: istate,ispecies,ishift
      integer :: isorth(0:ns,0:nstages-1), nclh(0:ns), nslh
      real (kind=R8) ::
     *   zah(0:ns-1), amh(0:ns-1), t(0:ns),
     &   den(0:ns,0:nstages-1), ds(0:ns,0:nstages-1,2),
     &   denl(0:ns,0:nstages-1), tl(0:ns),
     &   denr(0:ns,0:nstages-1), tr(0:ns),
     &   vtn(0:ns,0:nstages-1,-1:ny), 
     &   kappa(0:ns,0:nstages-1,-1:ny),
     &   chi(0:ns,0:nstages-1,-1:ny), vw(0:ns,0:nstages-1,-1:ny),
     &   xi(0:ns,0:nstages-1), ntau(0:ns,0:ns), 
     &   dfhcb, dfhce, hthet(-1:nx,-1:ny),
     &   cellface(0:ns,0:nstages-1,-1:ny,1:5),
     &   hy1(-1:nx,-1:ny), vwt(0:ns,0:nstages-1,-1:ny)
      real (kind=R8) ::
     *   zacomp, amcomp, epsn
*   ..set the accuracy
      parameter (epsn=1.e-3)
*   ..set the contribution ic
*     0 --- classical particle flux
*     1 --- banana plateau contribution
*     2 --- Pfirsch Schlueter contribution
*     3 --- both banana and PS
*     4 --- all contributions
*     mind that B2 already calculates the classical transport ! 
*     avoid double transport, 0+4 for cross checks only !
      parameter (ic=3)
      real (kind=R8) ::
     *   thet, ilpb, hxface, ilp(-1:ny), bmax, nlambda, nfactor,
     &   hxpface, hxnface, eparr, rho, prefactor, postfactor
      save eparr
      real (kind=R8) ::
     *   fluxav, sfluxav, tfluxav
      logical closedf
*   ..for consistency check
      real (kind=R8) ::
     *   densav(0:ns,0:nstages-1,-1:ny),
     &   gradnsav(0:ns,0:nstages-1,-1:ny),
     &   tforce(0:ns,0:nstages-1,1:2,-1:ny)

*   ..procedures
      external ipgeti, ipgetr, xertst, b2xvsg, fluxav, sfluxav, 
     &         tfluxav, neodv, neochi

*   ..subprogram start-up calls
      call subini ('neo_diagnostics')
*   ..neoclassical transport only valid for single charge states
*     bundled charge states yet to be implemented
      do is=0,ns-1
        call xertst (zamax(is).eq.zamin(is),
     &   'ERROR: neoclassical transport cannot be calculated')
c     1     // ' with bundled charge states'
      enddo
      call alloc_b2mod_neoclassical(nx,ny,ns,nstages,mmx)
      if(nnreg(0).gt.1 .and. nnreg(0).le.5) then
        ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
        iyref=2*topcut(1)/3
      elseif(nnreg(0).eq.8) then
        ixref=(rightcut(1)+rightcut(2))/2
        iyref=2*min(topcut(1),topcut(2))/3
      else
        ixref=3*nx/4
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

*   ..initialize densities and temperatures
      do is=0,ns
        t(is)=0.0_R8
        do istate=0,nstages-1
          den(is,istate)=0.0_R8
        enddo
      enddo

*   ..calculate hy-projection
      do ix=-1,nx
        do iy=-1,ny
          hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
        enddo
      enddo

*   ..species reordering and calculation of geometric parameters

*    ..loop voltage eparr (to be included in a separate routine!)
*       default value:
        eparr=0.0_R8                                  
*       read scalar from b2mn.dat
        call ipgetr ('b2tqna_eparr', eparr)
        if (eparr.ne.0.0_R8 .or. b2plot_debug.ge.3) 
     >   write(*,*) 'eparr =', eparr

*   ..initialize geometric parameters
        do iy=-1,ny
          bav_(iy)=0.0_R8
          b2av_(iy)=0.0_R8
          bi2a_(iy)=0.0_R8
          rbt_(iy)=0.0_R8
          bgradp_(iy)=0.0_R8
          dpsidr_(iy)=0.0_R8
          rnq_(iy)=0.0_R8
          fc_(iy)=0.0_R8
          gclass_(iy)=0.0_R8
          do is=1,mmx
            fm_(is,iy)=0.0_R8
          enddo
          r2i_(iy)=0.0_R8
          rho_(iy)=0.0_R8
        enddo 

*   ..split species into species and charge states
        nar=ns
        nzm=nstages-1
        do is=0,ns
          ncl(is)=0
          do istate=0,nstages-1
            isort(is,istate)=0
          enddo
        enddo
        do is=0,ns-1
          zah(is)=rza(-1,-1,is)  ! assume constant charge state
          amh(is)=am(is)
        enddo
        do is=0,ns-2
          zacomp=zah(is)
          amcomp=amh(is)
          ispecies=is
          do is2=is+1,ns-1
            if((abs(amh(is2)-amcomp).lt.(0.5)).and.
     &         (zah(is2).lt.zacomp)) then
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
         if((amh(is)-amcomp).gt.(0.5)) then
           ncl(ispecies-1)=istate
           ispecies=ispecies+1
           amcomp=amh(is)
           istate=1
         elseif(abs(amh(is)-amcomp).lt.(0.5)) then
           istate=istate+1
         endif
        enddo

        ncl(ispecies-1)=istate

        nsl=ispecies

        is=0
        do ispecies=0,nsl-1
          do istate=0,ncl(ispecies)-1
            do is2=0,ns-1 
              if((abs(am(is2)-amh(is)).lt.(0.5)).and.
     &           (abs(rza(-1,-1,is2)-zah(is)).lt.0.5)) then
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
            
*   ..mass and charges
        do ispecies=0,nsl-1
          m(ispecies)=am(isort(ispecies,0))*mp
          do istate=0,ncl(ispecies)-1
            zsp(ispecies,istate)=rza(-1,-1,isort(ispecies,istate))
          enddo
        enddo

*   ..add electrons as a separate species
        nsl=nsl+1
        ncl(nsl-1)=1
        m(nsl-1)=me
        zsp(nsl-1,0)=-1.0_R8

*   ..calculate generalized poloidal angle theta
*      default value 0 for grid points outside the core and on
*      the inner boundary
        do iy=-1,ny
          do ix=-1,nx
            theta(ix,iy)=0._R8
          enddo
        enddo
        do iy=0,ny
          thet=0._R8
          ilp(iy)=0._R8
          do ix=ixref,-1,-1
           if (mod(region(ix,iy,0),4).eq.1) then 
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
          enddo
          do ix=nx,ixref+1,-1
           if (mod(region(ix,iy,0),4).eq.1) then 
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
          enddo
          theta(ixref,iy)=thet
          ilpb=thet
          do ix=-1,nx
           if (mod(region(ix,iy,0),4).eq.1) then 
            theta(ix,iy)=2._R8*pi*theta(ix,iy)/ilpb
           endif
          enddo
        enddo

*   ..calculate geometric parameters

*   ..loop over flux surfaces
        do iy=0,ny

*   ..calculate neoclassical transport only for closed flux surfaces
         closedf=.false.
         do ix=-1,nx
          if (mod(region(ix,iy,0),4).eq.1) then
            closedf=.true.
          endif
         enddo

*   ..calculate the connection length
         conn(iy)=0
         do ix=-1,nx
            if(mod(region(ix,iy,0),4).eq.1) conn(iy)=conn(iy)+
     1           hx(ix,iy)*bb(ix,iy,3)/abs(bb(ix,iy,0))
         enddo

         if (closedf) then
           closedf=.false.

*   ..flux surface averages

           bav_(iy)=fluxav(bb(:,:,3),bb,hx,vol,nx,ny,iy)
           b2av_(iy)=fluxav(bb(:,:,3)**2._R8,bb,hx,vol,nx,ny,iy)
           bi2a_(iy)=fluxav(bb(:,:,3)**(-2._R8),bb,hx,vol,nx,ny,iy)
           rbt_(iy)=fluxav(bb(:,:,2)*(crx(:,:,0)+crx(:,:,1)+
     &           crx(:,:,2)+crx(:,:,3))/4._R8,bb,hx,vol,nx,ny,iy)
           bgradp_(iy)=2._R8*pi/(bav_(iy)*ilp(iy))
           dpsidr_(iy)=bb(ixref,iy,0)*(crx(ixref,iy,0)+
     &            crx(ixref,iy,2))/2._R8
           rnq_(iy)=conn(iy)/(2._R8*pi)

           bmax=0._R8
           do ix=-1,nx
            if (mod(region(ix,iy,0),4).eq.1) then          
             if (bb(ix,iy,3).gt.bmax) then
              bmax=bb(ix,iy,3)
             endif
            endif
           enddo
           nlambda=1._R8/nints

           fc_(iy)=0
           do is2=0,nints
            if ((is2.eq.0).or.(is2.eq.nints)) then
               nfactor=0.5_R8
             else
               nfactor=1._R8
            endif
            fc_(iy)=fc_(iy)+is2*nlambda/sfluxav(bb(:,:,3),bb,hx,vol,
     &            nx,ny,iy,is2*nlambda,bmax)*nlambda*nfactor
           enddo
           fc_(iy)=0.75_R8*b2av_(iy)/bmax**2._R8*fc_(iy)

           gclass_(iy)=fluxav(((crx(:,:,0)+crx(:,:,1)+crx(:,:,2)+
     &            crx(:,:,3))/4._R8*bb(:,:,0)/bb(:,:,3))**2._R8,
     &            bb,hx,vol,nx,ny,iy)

           do is2=1,mmx
            do is=-1,ny
             do ix=-1,nx
               hthet(ix,is)=is2*theta(ix,iy)
             enddo
            enddo
            cosbl(is2)=fluxav(cos(hthet)*bb(:,:,3)*log(bb(:,:,3)),
     &                        bb,hx,vol,nx,ny,iy)
            cosb2(is2)=fluxav(cos(hthet)*bb(:,:,3)**2._R8,
     &                        bb,hx,vol,nx,ny,iy)
            sinbl(is2)=fluxav(sin(hthet)*bb(:,:,3)*log(bb(:,:,3)),
     &                        bb,hx,vol,nx,ny,iy)
            sinb2(is2)=fluxav(sin(hthet)*bb(:,:,3)**2._R8,
     &                        bb,hx,vol,nx,ny,iy)
           enddo
           do is2=1,mmx
             fm_(is2,iy)=2._R8*(is2*bgradp_(iy))**2._R8/
     &               (b2av_(iy)*bav_(iy))*
     &               (cosb2(is2)*cosbl(is2)+sinb2(is2)*sinbl(is2))
           enddo
           r2i_(iy)=fluxav(((crx(:,:,0)+crx(:,:,1)+crx(:,:,2)+
     &         crx(:,:,3))/4._R8)**(-2._R8),bb,hx,vol,nx,ny,iy)

*    ..flux surface label
           rho_(iy)=(crx(ixref,iy,0)+crx(ixref,iy,2))/2._R8

          else
         endif

        enddo


*   ..now calculate the neoclassical transport

*   ..loop over flux surfaces
       do iy=0,ny

*   ..calculate neoclassical transport only for closed flux surfaces
        closedf=.false.
        do ix=-1,nx
         if (mod(region(ix,iy,0),4).eq.1) then
           closedf=.true.
         endif
        enddo
        if (closedf) then
          closedf=.false.

*   ..flux surfaces averages
*     temperatures in keV, densities in 10^19 m^-3
          do ispecies=0,nsl-2
            t(ispecies)=fluxav(ti,bb,hx,vol,nx,ny,iy)/
     &                  (1.e3_R8*ev)
            tl(ispecies)=fluxav(ti,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &                  (1.e3_R8*ev)
            tr(ispecies)=tfluxav(ti,bb,hx,vol,nx,ny,iy)/
     &                  (1.e3_R8*ev)
            do istate=0,ncl(ispecies)-1
              den(ispecies,istate)=
     &          fluxav(na(:,:,isort(ispecies,istate)),bb,hx,vol,
     &          nx,ny,iy)/1.e19_R8
              denl(ispecies,istate)=
     &          fluxav(na(:,:,isort(ispecies,istate)),bb,hx,vol,
     &          nx,ny,bottomiy(ixref,iy))/1.e19_R8
              denr(ispecies,istate)=
     &          tfluxav(na(:,:,isort(ispecies,istate)),bb,hx,vol,
     &          nx,ny,iy)/1.e19_R8
            enddo
          enddo
          t(nsl-1)=fluxav(te,bb,hx,vol,nx,ny,iy)/(1.e3_R8*ev)
          den(nsl-1,0)=fluxav(ne,bb,hx,vol,nx,ny,iy)/1.e19_R8
          tl(nsl-1)=fluxav(te,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &               (1.e3_R8*ev)
          denl(nsl-1,0)=fluxav(ne,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &                1.e19_R8
          tr(nsl-1)=tfluxav(te,bb,hx,vol,nx,ny,iy)/(1.e3_R8*ev)
          denr(nsl-1,0)=tfluxav(ne,bb,hx,vol,nx,ny,iy)/1.e19_R8


*   ..thermodynamic forces
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              ds(ispecies,istate,1)=(log(denr(ispecies,istate)*
     &           tr(ispecies))-log(denl(ispecies,istate)*
     &           tl(ispecies)))/(hy1(bottomix(ixref,iy),
     &           bottomiy(ixref,iy))+hy1(ixref,iy))
              ds(ispecies,istate,2)=(log(tr(ispecies))-
     &           log(tl(ispecies)))/(hy1(bottomix(ixref,iy),
     &           bottomiy(ixref,iy))+hy1(ixref,iy))


*   ..for consistency check
              densav(ispecies,istate,iy)=den(ispecies,istate)
              gradnsav(ispecies,istate,iy)=den(ispecies,istate)*
     &            (ds(ispecies,istate,1)-ds(ispecies,istate,2))
              tforce(ispecies,istate,1,iy)=ds(ispecies,istate,1)
              tforce(ispecies,istate,2,iy)=ds(ispecies,istate,2)
            enddo
          enddo

*   ..recover geometric parameters for the flux surface iy
          bav=bav_(iy)
          b2av=b2av_(iy)
          bi2a=bi2a_(iy)
          rbt=rbt_(iy)
          bgradp=bgradp_(iy)
          dpsidr=dpsidr_(iy)
          rnq=rnq_(iy)
          fc=fc_(iy)
          gclass=gclass_(iy)
          do is2=1, mmx
            fm(is2)=fm_(is2,iy)
          enddo
          r2i=r2i_(iy)
          rho=rho_(iy)

*    ..calculate the diffusion coefficients
          call neodv(ns,nstages,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &               epsn,dd(:,:,iy),vvn(:,:,iy),vtn(:,:,iy))

*    ..calculate the heat conduction coefficients
          call neochi(ns,nstages,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &                epsn,kappa(:,:,iy),chi(:,:,iy),vw(:,:,iy))

*    ..calculate thermal conductivity and heat pinch from
*        q(i,j) = - kappa(i,j)*grad.T(i) + vw(i,j)*T(i)
*               = - kappa(i,j)*grad.T(i) + den(i,j)*vwt(i,j)*T(i)
          do ispecies=0, nsl-1
            do istate=0, ncl(ispecies)-1
              vwt(ispecies,istate,iy)=vw(ispecies,istate,iy)/
     &                                (den(ispecies,istate)*1e19)
           enddo
          enddo

         else
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              dd(ispecies,istate,iy)=0._R8
              vvn(ispecies,istate,iy)=0._R8
              chi(ispecies,istate,iy)=0._R8
              kappa(ispecies,istate,iy)=0._R8
              vwt(ispecies,istate,iy)=0._R8
              densav(ispecies,istate,iy)=0._R8
              gradnsav(ispecies,istate,iy)=0._R8
              tforce(ispecies,istate,1,iy)=0._R8
              tforce(ispecies,istate,2,iy)=0._R8
            enddo
          enddo
        endif

       enddo

*    ..set transport coefficients on boundary by symmetric condition
       do ispecies=0,nsl-1
         do istate=0,ncl(ispecies)-1
           densav(ispecies,istate,-1)=densav(ispecies,istate,0)
           gradnsav(ispecies,istate,-1)=gradnsav(ispecies,istate,0)
           tforce(ispecies,istate,1,-1)=tforce(ispecies,istate,1,0)
           tforce(ispecies,istate,2,-1)=tforce(ispecies,istate,2,0)
           dd(ispecies,istate,-1)=dd(ispecies,istate,0)
           vvn(ispecies,istate,-1)=vvn(ispecies,istate,0)
           chi(ispecies,istate,-1)=chi(ispecies,istate,0)
           kappa(ispecies,istate,-1)=kappa(ispecies,istate,0)
           vwt(ispecies,istate,-1)=vwt(ispecies,istate,0)
         enddo
       enddo

*    ..calculate the transport coefficients on cell centers
       do iy=-1,ny
         do ispecies=0,nsl-1
           do istate=0,ncl(ispecies)-1
             cellface(ispecies,istate,iy,1)=dd(ispecies,istate,iy)
             cellface(ispecies,istate,iy,2)=vvn(ispecies,istate,iy)
             cellface(ispecies,istate,iy,3)=chi(ispecies,istate,iy)
             cellface(ispecies,istate,iy,4)=kappa(ispecies,istate,iy)
             cellface(ispecies,istate,iy,5)=vwt(ispecies,istate,iy)
           enddo
         enddo
       enddo

       do iy=0,ny
         if (mod(region(ixref,iy,0),4).eq.1) then
           if (iy.eq.ny) then
             prefactor=-hy(ixref,iy)/(2._R8*hy(ixref,
     &                  bottomiy(ixref,iy)))
             postfactor=1._R8+hy(ixref,iy)/(2._R8*hy(ixref,
     &                  bottomiy(ixref,iy)))
             ishift=bottomiy(ixref,iy)
            else
             if (mod(region(ixref,topiy(ixref,iy),0),4).eq.1) then
               prefactor=0.5_R8
               postfactor=0.5_R8
               ishift=iy
              else
               prefactor=-hy(ixref,iy)/(2._R8*hy(ixref,
     &                    bottomiy(ixref,iy)))
               postfactor=1._R8+hy(ixref,iy)/(2._R8*hy(ixref,
     &                    bottomiy(ixref,iy)))
               ishift=bottomiy(ixref,iy)
             endif
           endif
           do ispecies=0, nsl-1
             do istate=0, ncl(ispecies)-1
               dd(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,1)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),1)*postfactor
               vvn(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,2)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),2)*postfactor
               chi(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,3)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),3)*postfactor
               kappa(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,4)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),4)*postfactor
               vwt(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,5)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),5)*postfactor
             enddo
           enddo
         endif
       enddo

*     ..calculate transport coefficients at inner boundary by
*       a mixture of extrapolation and symmetry
       prefactor=0.5_R8*(hy(ixref,-1)+hy(ixref,topiy(ixref,-1)))/
     &     (hy(ixref,topiy(ixref,-1))+hy(ixref,topiy(ixref,
     &      topiy(ixref,-1))))
       do ispecies=0, nsl-1
         do istate=0, ncl(ispecies)-1
           if (mod(region(ixref,-1,0),4).eq.1) then
             dd(ispecies,istate,-1)=dd(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(dd(ispecies,istate,
     &          topiy(ixref,-1))-dd(ispecies,istate,
     &          topiy(ixref,topiy(ixref,-1))))
             vvn(ispecies,istate,-1)=vvn(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(vvn(ispecies,istate,
     &          topiy(ixref,-1))-vvn(ispecies,istate,
     &          topiy(ixref,topiy(ixref,-1))))
             chi(ispecies,istate,-1)=chi(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(chi(ispecies,istate,
     &          topiy(ixref,-1))-chi(ispecies,istate,
     &          topiy(ixref,topiy(ixref,-1))))
             kappa(ispecies,istate,-1)=kappa(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(kappa(ispecies,istate,
     &          topiy(ixref,-1))-kappa(ispecies,istate,
     &          topiy(ixref,topiy(ixref,-1))))
             vwt(ispecies,istate,-1)=vwt(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(vwt(ispecies,istate,
     &          topiy(ixref,-1))-vwt(ispecies,istate,
     &          topiy(ixref,topiy(ixref,-1))))
            else
             dd(ispecies,istate,-1)=0.
             vvn(ispecies,istate,-1)=0.
             chi(ispecies,istate,-1)=0.
             kappa(ispecies,istate,-1)=0.
             vwt(ispecies,istate,-1)=0.
           endif
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
      dn_neo(:,:,:)=0._R8
      chi_neo(:,:,:)=0._R8
      vn_neo(:,:,:)=0._R8
      vw_neo(:,:,:)=0._R8

      do ix=-1,nx
       do iy=-1,ny
        if (mod(region(ix,iy,0),4).eq.1) then
         do ispecies=0,nsl-2
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dn_neo(ix,iy,is)=dd(ispecies,istate,iy)*
     &                      (hy(ix,iy)/hy(ixref,iy))**2.
           vn_neo(ix,iy,is)=vvn(ispecies,istate,iy)*
     &                      (hy(ix,iy)/hy(ixref,iy))
           dfhcb=chi(ispecies,istate,iy)*
     &           (hy(ix,iy)/hy(ixref,iy))**2.
           chi_neo(ix,iy,is)=dfhcb
           vw_neo(ix,iy,is)=vwt(ispecies,istate,iy)*
     &                      (hy(ix,iy)/hy(ixref,iy))
          enddo
         enddo
         dfhce=chi(nsl-1,0,iy)*(hy(ix,iy)/hy(ixref,iy))**2.
         if (dfhce.lt.0) then
           write(*,*) 'dfhce(',ix,',',iy,')=',dfhce
         endif
         dn_neo(ix,iy,ns)=dd(nsl-1,0,iy)*
     &                    (hy(ix,iy)/hy(ixref,iy))**2.
         vn_neo(ix,iy,ns)=vvn(nsl-1,0,iy)*
     &                    (hy(ix,iy)/hy(ixref,iy))
         chi_neo(ix,iy,ns)=dfhce
         vw_neo(ix,iy,ns)=vwt(nsl-1,0,iy)*(hy(ix,iy)/hy(ixref,iy))

*   ..set transport to zero in the SOL and divertor regions
        elseif (mod(region(ix,iy,0),4).gt.1) then
         do ispecies=0,nsl-1
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dn_neo(ix,iy,is)=0._R8
           vn_neo(ix,iy,is)=0._R8
           chi_neo(ix,iy,is)=0._R8
           vw_neo(ix,iy,is)=0._R8
          enddo
         enddo
        endif

       enddo
      enddo
*   ..return
      call subend ()
      return
      end
