*-----------------------------------------------------------------------
*.specification

      subroutine b2tqce (nx, ny, lnlam, bb, te, ti, ne, ne2, hce, sig,
     &                   alf)  !srv 01.07.09
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *  lnlam(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), te(-1:nx,-1:ny),        !srv 01.07.09
     &  ti(-1:nx,-1:ny), ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  hce(-1:nx,-1:ny,0:1), sig(-1:nx,-1:ny,0:1),
     &  alf(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TQCE computes the classical contributions to the electron
*     transport coefficients.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the classical contributions to the following
*     electron transport coefficients.
*     - hce: electron thermal conductivity.
*     - sig: electrical conductivity.
*     - alf: thermo-electric coefficient.
*     These coefficients occur in the equations as indicated below.
*     Electron heat flux.
*       fhex = alfx*te*ehx-hcex*(d/dx)te
*       fhey = alfy*te*ehy-hcey*(d/dy)te
*     Electric current.
*       fchx = sigx*ehx-alfx*(d/dx)te
*       fchy = sigy*ehy-alfy*(d/dy)te
*     Here, ehx and ehy are the components of the modified electric
*     field: ehx = -(d/dx)po+(1/(qe*ne))*(d/dx)(ne*te),
*     ehy = -(d/dy)po+(1/(ne*te))*(d/dy)(ne*te). See Balescu /1/.
*
*     The equations indicated above do not describe all the relevant
*     transport properties of a plasma. The most significant omission
*     is the absence of diamagnetic transport terms.
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
      integer ncall, n2, ix, iy, model
	integer, save :: icase = 0                                      !srv 20.09.11
	integer, save :: iout = 0,                                      !srv 16.06.08
     &                 iout_b2nph9 = 0, iout_b2npht = 0, fke_Zhdanov = 1
      integer, save :: style_guard_cells = 0                          !srv 19.04.13
      character(23) hlp_frm                                           !srv 11.03.13
	real (kind=R8), save :: fch_pTe = 0.0_R8                        !srv 30.12.09
      real (kind=R8) ::
     *  ctaue, ce, t0, t1, z, tau(-1:nx,-1:ny), tnp                   !srv 11.09.09
!xpb      parameter (lnlam=12.0_R8)
      parameter (ce = 3.16e0_R8)
      save ncall, model
*   ..procedures
      intrinsic abs, sqrt, max, log
      real (kind=R8) ::
     *  pitx, fke, fse, fal, fke_Zh
      external subini, subend, xertst, sfill
      external b2xvsg, ipgeti
      pitx(ix,iy) = abs(bb(ix,iy,0))/bb(ix,iy,3)
      fke(z) =  2.16_R8*z/(1.0_R8+0.27_R8*z)
      fse(z) =  4.95_R8*z/(1.0_R8+1.54_R8*z)
      fal(z) = -1.25_R8*z/(1.0_R8+0.41_R8*z)
      fke_Zh(z) = (3.9_R8+2.3_R8/z)/(0.31_R8+1.2_R8/z+0.41_R8/z**2)   !srv 02.12.10
      real (kind=R8) c11e, c13e, c15e
      parameter (c11e = 1.0000_R8, c13e = 0.9487_R8, c15e =-0.8964_R8)
      real (kind=R8) c31e
      parameter (c31e = c13e)
      real (kind=R8) c33e, c35e, c55e
      c33e(z) = 1.3000_R8 + 0.5657_R8/z
      c35e(z) =-1.3040_R8 - 0.3207_R8/z
      c55e(z) = 1.5464_R8 + 0.9091_R8/z
      real (kind=R8) c53e
      c53e(z) = c35e(z)
      real (kind=R8) d13e, d15e, d35e
      d13e(z) = c11e   *c33e(z) - c13e   *c13e
      d15e(z) = c11e   *c55e(z) - c15e   *c15e
      d35e(z) = c33e(z)*c55e(z) - c35e(z)*c35e(z)
      real (kind=R8) h135e, h153e, h315e
      h135e(z) = c15e   *c33e(z) - c13e   *c35e(z)
      h153e(z) = c13e   *c55e(z) - c15e   *c53e(z)
      h315e(z) = c35e(z)*c11e    - c31e   *c15e
      real (kind=R8) f135e
      f135e(z) = c55e(z)*d13e(z) - c35e(z)*h315e(z) - c15e   *h135e(z)
*   ..initialisation
      data ncall/0/, model/1/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lnlam is the log of the plasma parameter.
*
*     Description of some local procedures:
*
*     fke(z) is a fit to the dimensionless parallel electron thermal
*     conductivity as function of the atomic charge.
*     fse(z) is a fit to the dimensionless parallel electrical
*     conductivity as function of the atomic charge.
*     fal(z) is a fit to the dimensionless parallel thermo-electric
*     coefficient as function of the atomic charge.
*     (see fig. 4.1 of Balescu for fke, fse, fal).
*     pitx(ix,iy) is the local pitch of the magnetic field.
*     The coefficients and functions c11e, c13e, c15e, c33e, c35e,
*     c53e, c55e, d13e, d15e, d35e, h135e, h153e, h315e, and f135e
*     are the corresponding expressions from Balescu from computation
*     of his 13-moment and 21-moment expressions for the electron
*     transport coefficient. See /1/.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqce')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..read input data
      if (ncall.eq.0) then
        call ipgeti ('b2tqce_model', model)
        call ipgeti ('b2tqce_iout', iout)                            !srv 16.06.08                       
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                           !srv 30.12.09
        call ipgeti ('b2tqce_fke_Zhdanov', fke_Zhdanov)                   !srv 02.12.10
        call ipgeti ('b2tlnl_ee', icase)
        call xertst (0.eq.icase .or. 2.eq.icase, 
     &               'faulty input b2tlnl_ee')
        call ipgeti ('b2tqce_style_guard_cells', style_guard_cells)       !srv 19.04.13
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif
*   ..compute ctaue
*     Following Balescu. The mathematical expression is:
*       t0 = (3/(4*sqrt(2*pi)))*(4*pi*eps0/qe**2)**2*ev*sqrt(ev)
*       ctaue = t0*sqrt(me)/lnlam
*     The numerical value in SI units is: ctaue = 3.4407e11/lnlam.
*     The following calculation of t0 is designed to avoid overflow
*     even for IEEE single precision arithmetic.
!xpb The lnlam term has been removed and is put back in later so it may
!xpb vary with space and species or be set by user.
      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      if (model.eq.1.or.model.eq.3) then
        ctaue = t0*sqrt(me)						 !xpb
      else if (model.eq.2) then
        ctaue = 3.5e11_R8             ! lower precision as per SOLPS4.0
      endif

*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)               !srv 01.07.09  20.09.11
* ..compute classical hce, sig, alf
*   ..compute hcex, sigx, alfx.
      if(style_guard_cells.eq.0) then
       do iy = -1, ny
        do ix = -1, nx
         t0 = te(ix,iy)/ev
         tau(ix,iy) = (ctaue/lnlam(ix,iy))*t0*sqrt(t0)/ne2(ix,iy)	   !xpb !srv 01.07.09
         tnp = tau(ix,iy)*ne(ix,iy)*pitx(ix,iy)**2
         z = ne2(ix,iy)/ne(ix,iy)
         if (model.eq.1) then  ! SOLPS5.0 model : Balescu fit
	    if(fch_pTe .eq. 1.0_R8.and.fke_Zhdanov.eq.1) then            !srv 02.12.10
           hce(ix,iy,0) =     (te(ix,iy)/me)*fke_zh(z)*tnp
	    else
           hce(ix,iy,0) = 2.5_R8*(te(ix,iy)/me)*fke(z)*tnp
	    endif
           sig(ix,iy,0) = (qe**2/me)*fse(z)*tnp
           alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*fal(z)*tnp
         else if (model.eq.2) then ! old SOLPS4.0 Braginskii model for hce
           hce(ix,iy,0) = ce*(te(ix,iy)/me)*tnp
           sig(ix,iy,0) = (qe**2/me)*fse(z)*tnp
           alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*fal(z)*tnp
         else if (model.eq.3) then  ! Balescu exact 21-moment expressions
          hce(ix,iy,0) = 2.5_R8*(te(ix,iy)/me)*(d15e(z)/f135e(z))*tnp
          sig(ix,iy,0) = (qe**2/me)*(d35e(z)/f135e(z))*tnp
          alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*(-h153e(z)/f135e(z))*tnp
         endif
        enddo
       enddo
      elseif(style_guard_cells.eq.1) then
       call sfill (n2, 0.0_R8, hce(-1,-1,0), 1)  ! IYS and lk 21.11.2012
       call sfill (n2, 0.0_R8, sig(-1,-1,0), 1)  ! IYS and lk 21.11.2012
       call sfill (n2, 0.0_R8, alf(-1,-1,0), 1)  ! IYS and lk 21.11.2012

       do iy = -1, ny
        do ix = -1, nx
         t0 = te(ix,iy)/ev
         tau(ix,iy) = (ctaue/lnlam(ix,iy))*t0*sqrt(t0)/ne2(ix,iy)	   !xpb !srv 01.07.09
         tnp = tau(ix,iy)*ne(ix,iy)*pitx(ix,iy)**2
         z = ne2(ix,iy)/ne(ix,iy)
         if (model.eq.1) then  ! SOLPS5.0 model : Balescu fit
	    if(fch_pTe .eq. 1.0_R8.and.fke_Zhdanov.eq.1) then            !srv 02.12.10
           hce(ix,iy,0) =     (te(ix,iy)/me)*fke_zh(z)*tnp
	    else
           hce(ix,iy,0) = 2.5_R8*(te(ix,iy)/me)*fke(z)*tnp
	    endif
          if ((iy .ne. -1 .and. iy .ne. ny) 
     &                .or. mod(region(ix,iy,0),4).eq.1) then   ! IYS and lk 21.11.2012
            sig(ix,iy,0) = (qe**2/me)*fse(z)*tnp
            alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*fal(z)*tnp
          endif
         else if (model.eq.2) then ! old SOLPS4.0 Braginskii model for hce
          hce(ix,iy,0) = ce*(te(ix,iy)/me)*tnp
          if ((iy .ne. -1 .and. iy .ne. ny) 
     &                .or. mod(region(ix,iy,0),4).eq.1) then   ! IYS and lk 21.11.2012
            sig(ix,iy,0) = (qe**2/me)*fse(z)*tnp
            alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*fal(z)*tnp
          endif
         else if (model.eq.3) then  ! Balescu exact 21-moment expressions
          hce(ix,iy,0) = 2.5_R8*(te(ix,iy)/me)*(d15e(z)/f135e(z))*tnp
          if ((iy .ne. -1 .and. iy .ne. ny) 
     &                .or. mod(region(ix,iy,0),4).eq.1) then   ! IYS and lk 21.11.2012
            sig(ix,iy,0) = (qe**2/me)*(d35e(z)/f135e(z))*tnp
            alf(ix,iy,0) = sqrt(2.5_R8)*(qe/me)*(-h153e(z)/f135e(z))*tnp
          endif
         endif
        enddo
       enddo
	else
        CALL XERRAB ('B2TQCE - incorrect value of style_guard_cells')
	endif
*   ..compute hcey, sigy, alfy
*!!   (classical perpendicular transport is set to 0)
      call sfill (n2, 0.0_R8, hce(-1,-1,1), 1)
      call sfill (n2, 0.0_R8, sig(-1,-1,1), 1)
      call sfill (n2, 0.0_R8, alf(-1,-1,1), 1)

      if(iout.ne.0) then                                             !srv 16.06.08 {
        call my_out(70,nx,ny,tau,'taue')
        call my_out(70,nx,ny,hce(-1,-1,0),'hcex')
        call my_out(70,nx,ny,hce(-1,-1,1),'hcey')
        call my_out(70,nx,ny,sig(-1,-1,0),'sigx')
        call my_out(70,nx,ny,sig(-1,-1,1),'sigy')
        call my_out(70,nx,ny,alf(-1,-1,0),'alfx')
        call my_out(70,nx,ny,alf(-1,-1,1),'alfy')
        open (70,file='b2tqce_fke.dat')                         
        open (71,file='b2tqce_fke_Zh.dat')                         
        write(hlp_frm,"(a10,i4,a9)") '(1x,3x,1x,',nx+2,'(es14.6))'  !srv 11.03.13
        do iy = -1, ny
          z = ne2(ix,iy)/ne(ix,iy)
	    write(70,hlp_frm) (fke(z), ix = -1, nx)                    !srv 11.03.13
	    write(71,hlp_frm) (fke_Zh(z), ix = -1, nx)                 !srv 11.03.13
        enddo
	close(71)
	close(70) 
      endif                                                          !srv 16.06.08 }
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
        call my_out(70,nx,ny,tau,'b2tqce_taue')
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqce

      end
