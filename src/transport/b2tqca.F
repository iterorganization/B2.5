*-----------------------------------------------------------------------
*.specification

      subroutine b2tqca (nx, ny, ns, lnlam, na, te, ti, ne, ne2,       !srv 01.07.09
     1 vsa, hci)
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_external
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  lnlam(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),                        !srv 01.07.09
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &  ne2(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vsa(-1:nx,-1:ny,0:1,0:ns-1), hci(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TQCA computes the classical contributions to the atom
*     transport coefficients.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the classical contributions to the following
*     atom transport coefficients.
*     - vsa: viscosity.
*     - hci: atom thermal conductivity.
*     These coefficients occur in the equations as indicated below.
*     Parallel momentum transport.
*       momx(is) = -vsax(is)*(d/dx)up(is)
*       momy(is) = -vsay(is)*(d/dy)up(is)
*     Atom heat flux.
*       fhix = -hcix*(d/dx)ti
*       fhiy = -hciy*(d/dy)ti
*     The particle and parallel momentum transport is a property of
*     each atomic species, the other transport laws involve summation
*     over all atomic species.
*
*     The equations indicated above do not describe all the relevant
*     transport properties of a plasma. The most significant omission
*     is the absence of diamagnetic transport terms. Also, one might
*     define temperatures and heat fluxes for each individual atom
*     species. These tasks are left for future work.
*
*
*  4. references
*
*     /1/ R. Balescu: Transport Processes in Plasmas, Vol. 1.
*     North Holland, Amsterdam, 1988, pp. 230-239.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, model, style, n2, ix, iy, is, js
      integer, save :: icase = 0                                      !srv 20.09.11
      integer, save :: iout = 0, iout_b2wdat=0,                       !srv 16.06.08
     &                 iout_b2nph9 = 0, iout_b2npht = 0, iout_b2npmo = 0
      character chns*3                                                !srv 16.06.08
      real (kind=R8) :: 
     *  fki, fxi, eta, ctaup, ctaui, ci, ft,
     *  t, t0, t1, tau(-1:nx,-1:ny), tnp, phm0, w                     !srv 11.09.09
      parameter (fki=2.253_R8, eta=0.96_R8, ci = 3.9e0_R8)
      save ncall, model, phm0, style
*   ..procedures
      intrinsic sqrt
      external xertst, sfill, ipgetr, ipgeti
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, model/1/, style/1/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lnlam is the log of the plasma parameter.
*     fki is the dimensionless parallel atom thermal conductivity.
*     eta is the dimensionless parallel atom viscosity.
*     (see table 4.1 of Balescu for fki and eta.)
*
*     Description of some local procedures:
*
*     pitx(ix,iy) is the local pitch of the magnetic field.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqca')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read input parameter
      if (ncall.eq.0) then
        call ipgeti ('b2tqca_model', model)
        call ipgetr ('b2tqca_phm0', phm0)
        call xertst (0.0_R8.le.phm0, 'faulty input phm0')
        call ipgeti ('b2tral_cvsa_style', style)
        call ipgeti ('b2tqca_iout', iout)                            !srv 16.06.08                       
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npmo_iout', iout_b2npmo)                     !srv 30.01.10
        call ipgeti ('b2tlnl_ii', icase)
        call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
        call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
        call b2xvsg (n2, ti, 1, 'ti', '.gt.')
        call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaup
*     Following Balescu; Braginskii has ctaup different by a factor
*     sqrt(2). The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaup = t0*sqrt(mp)/lnlam
*     The numerical value in SI units is: ctaup = 1.4743e13/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
!xpb The lnlam term has been removed and is put back in later so it may
!xpb vary with space and species or be set by user.
      if (style.eq.1) then
        ft = 1.0_R8
      else
        ft = 4.0_R8/3.0_R8
      endif
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      ctaup = t0*sqrt(mp)                                          !xpb

*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)               !srv 01.07.09  20.09.11
* ..compute vsa, hci
*   ..compute vsax, hcix
*     (initialise hcix to 0)
      call sfill (n2, 0.0_R8, hci(-1,-1,0), 1)
      do is = 0, ns-1
       if (is_neutral(is)) then
*       (neutral vsax and hcix are strictly anomalous)
        call sfill (n2, 0.0_R8, vsa(-1,-1,0,is), 1)
       else
        do iy = -1, ny
         do ix = -1, nx
          t0 = ti(ix,iy)/ev
          tau(ix,iy) = (ctaup/lnlam(ix,iy))*sqrt(am(is))*t0*sqrt(t0)/             !xpb !srv 01.07.09
     &      (rz2(ix,iy,is)*ne2(ix,iy))
          tnp = tau(ix,iy)*na(ix,iy,is)*pitx(ix,iy)**2
          if (model.eq.1) then  ! SOLPS5.0 model : Balescu 29-moment result
            hci(ix,iy,0) = hci(ix,iy,0)+
     &       2.5_R8*(ti(ix,iy)/(mp*am(is)))*fki*tnp
c DPC/XPB changed to 4/3*  at the strong recommendation of Rozhansky
c WDK/SRV 21.07.2016 factor 4/3 removed from vsa, but added explicitly to cvsahz (not cvsa),
c         according to old model.eq.3 treatment
            vsa(ix,iy,0,is) = 
     &       phm0*ft*eta*sqrt(2.0_R8)*ti(ix,iy)*tnp
          else if (model.eq.2) then ! old SOLPS4.0 model : Braginskii
            w = 0.0_R8
            t = 0.0_R8
            do js = 0, ns-1
              w = w + rz2(ix,iy,js) * na(ix,iy,js) * 
     &            sqrt(2.0_R8*mp*am(is)*am(js)/(am(is)+am(js))) 
              t = t + rz2(ix,iy,js) * na(ix,iy,js) / 
     &            sqrt((am(is)+am(js))/2.0_R8) 
            enddo
            do js = 0, ns_ext-1
              w = w + za2_ext(ix,iy,js) * na_ext(ix,iy,js) *
     &            sqrt(2.0_R8*mp*am(is)*am_ext(js)/(am(is)+am_ext(js)))
              t = t + za2_ext(ix,iy,js) * na_ext(ix,iy,js) /
     &            sqrt((am(is)+am_ext(js))/2.0_R8)
            enddo
cxpb            ctaui = ctaup*sqrt(2.0_R8)/lnlam/rz2(ix,iy,is)  ! lower precision as per 4.0
            ctaui = 2.1e13_R8 / (lnlam(ix,iy)*rz2(ix,iy,is))             !srv 01.07.09
            fxi = ci * ctaui / (ev*sqrt(ev*mp))
            hci(ix,iy,0) = hci(ix,iy,0) + fxi*(na(ix,iy,is)/w)*
     &       pitx(ix,iy)**2*(ti(ix,iy)**2*sqrt(ti(ix,iy)))
cwdk      no factor 4/3 here, for consistent treatment with model.eq.1
            vsa(ix,iy,0,is) = 
     &       phm0*eta*ctaui/(ev*sqrt(ev))*(na(ix,iy,is)/t)*
     &       pitx(ix,iy)**2*(ti(ix,iy)**2*sqrt(ti(ix,iy)))
          else
            call xerrab ('b2tqca--invalid value b2tqca_model > 2. ')
          endif
         enddo
        enddo
        if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0
     &                      .or. iout_b2wdat.eq.4) then
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,tau,'b2tqca_taua'//chns)
        endif
       endif
      enddo
*   ..compute vsay, hciy
*!!   (classical perpendicular transport is set to 0)
      do is = 0, ns-1
       call sfill (n2, 0.0_R8, vsa(-1,-1,1,is), 1)
      enddo
      call sfill (n2, 0.0_R8, hci(-1,-1,1), 1)

      if(iout_b2npmo.ne.0) then                                        !srv 16.06.08 30.01.10 {
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,vsa(-1,-1,0,is),'b2tqca_vsax'//chns)
        call my_out(70,nx,ny,vsa(-1,-1,1,is),'b2tqca_vsay'//chns)
!        call my_out(70,nx,ny,na(-1,-1,is),'b2tqca_na'//chns)
       enddo
!        call my_out(70,nx,ny,lnlam,'b2tqca_lnlam')
!        call my_out(70,nx,ny,ti,'b2tqca_ti')
!        call my_out(70,nx,ny,ne,'b2tqca_ne')
!        call my_out(70,nx,ny,ne2,'b2tqca_ne2')
      endif                                                          !srv 16.06.08 }

      if(iout_b2wdat.eq.4) then                                      !srv 16.06.08 30.01.10 {
        call my_out(70,nx,ny,lnlam,'b2tqca_lnlam')
      endif

      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
        call my_out(70,nx,ny,hci(-1,-1,0),'b2tqca_hcix')
        call my_out(70,nx,ny,hci(-1,-1,1),'b2tqca_hciy')
      endif                                                          !srv 16.06.08 30.01.10 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqca

      end subroutine b2tqca
