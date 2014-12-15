*-----------------------------------------------------------------------
*.specification

      subroutine b2npmo (nx, ny, ns, itcnt, solving, rxf, hx, hy, hz, 
     &  vol, bb, pbs, pbshz, lnlam, na, te, ti, po, fna_fcor,        !srv 09.01.01 01.07.09
     &  cvsahz, cvsa, cvsa_cl, cvsahz_cl, smo, smq, cthe, cthi, chci, 
     &  ua, smof, wadia, vaecrb, resmo, resmt, corua, corut,         !srv 10.11.02
     &  pcca, pccm, ierr)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_sources
     1 , only : b2npmo_smaf, b2npmo_smag, b2npmo_smav
      use b2mod_plasma
     & , only : fmo
      use b2mod_numerics_namelist
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, itcnt
      real (kind=R8) ::
     *  rxf, 
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),           !srv 09.01.01
     &  vol(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), pbs(-1:nx,-1:ny,0:1),
     &  pbshz(-1:nx,-1:ny,0:1), lnlam(-1:nx,-1:ny),                  !srv 01.07.09
     &  na(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  po(-1:nx,-1:ny), fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvsahz(-1:nx,-1:ny,0:1,0:ns-1), cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     1  cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1), 
     &  cvsahz_cl(-1:nx,-1:ny,0:1,0:ns-1),                            !srv 09.01.01
     2  smo(-1:nx,-1:ny,0:3,0:ns-1), smq(-1:nx,-1:ny,0:3,0:ns-1),
     &  smof(-1:nx,-1:ny,0:3,0:ns-1), chci(-1:nx,-1:ny,0:1),          !srv 10.11.02
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1)
      logical solving
*   ..input/output arguments
      real (kind=R8) ::
     *  ua(-1:nx,-1:ny,0:ns-1)
*   ..output arguments (unspecified on entry)
      integer ierr
      real (kind=R8) ::
     *  resmo(-1:nx,-1:ny,0:ns-1), resmt(-1:nx,-1:ny),
     &  corua(-1:nx,-1:ny,0:ns-1), corut(-1:nx,-1:ny),
     &  pcca(-1:nx,-1:ny,0:1,0:ns-1), pccm(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NPMO performs one iteration of a relaxation process for the
*     solution of the b2 equations for parallel velocity.
*
*  3. description (see also routine b2cdca)
*
*     For a detailed description of the equations and the method of
*     solution, see the references cited in the main program unit.
*
*     This routine performs one inner iteration for the equations that
*     govern the parallel velocity ua.
*     An outline of the procedure follows.
*     ..compute auxiliaries
*     ..do is=0,ns-1
*      ..compute residual for species (is)
*      ..compute correction for species (is)
*      ..compute contribution to total residual
*     ..enddo
*     ..compute simultaneous correction for all species
*     ..apply correction to all species
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, diagno, ireg, k,ii
      integer, save :: modvis = 0                                    !srv 19.01.05
      integer, save :: iout = 0, b2sifr_ver=0                        !srv 04.07.08 26.10.08
      character chns*2, chk*1                                        !srv 04.07.08
      real (kind=R8) ::
     *  rxg, t0, t1, t2, vis_par
      real (kind=R8) ::
     *  smbch(-1:nx,-1:ny,0:3), smbtf(-1:nx,-1:ny,0:3),              !srv 11.09.09
     *  smbgp(-1:nx,-1:ny,0:3),
     *  smbvh(-1:nx,-1:ny,0:3), smbvv(-1:nx,-1:ny,0:3),              !srv 23.07.99 04.10.99
     *  smbcf(-1:nx,-1:ny,0:3),                                      !srv 09.07.01
     1  pb(-1:nx,-1:ny), cvct(-1:nx,-1:ny,0:1),
     2  rob(-1:nx,-1:ny), flub(-1:nx,-1:ny,0:1), flubv(-1:nx,-1:ny,0:1), !srv 11.09.09
     3  wrkf(-1:nx,-1:ny,0:1), aa5(-1:nx,-1:ny,-2:2),
     4  smb(-1:nx,-1:ny,0:3), flcb(-1:nx,-1:ny,0:1),
     5  cvcb(-1:nx,-1:ny,0:1), ro(-1:nx,-1:ny),
     6  pr(-1:nx,-1:ny), ne2(-1:nx,-1:ny),
     7  ni(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     8  ehxp(-1:nx,-1:ny), smt(-1:nx,-1:ny,0:3),
     9  flct(-1:nx,-1:ny,0:1), wrk1(-1:nx,-1:ny),
     1  smqdu(-1:nx,-1:ny), wrk0(-1:nx,-1:ny), wrk2(-1:nx,-1:ny),
     2  ctcfb(-1:nx,-1:ny), ctcft(-1:nx,-1:ny)                         !srv 06.04.07
      logical solvedum(0:NNREGMAX)
      save ncall, rxg, vis_par, diagno, solvedum
*   ..procedures
      intrinsic min, max
      real (kind=R8) ::
     *  med
      external subini, subend, xertst, ipgetr, sfill, scopy,
     &  sscal, saxpy
      external b2xvsg, b2xvff, b2xpne, b2xpro, b2xppr,
     &  b2xppb, b2nxdu, b2nxdv, b2nxfx, b2nxfm, b2nxfc, b2nxst,
     &  b2sifr, b2sigp, b2siav, b2sicf, b2nxcm, b2urmo, b2usmo
      med(t0,t1,t2) = min(max(t0,t1),max(t1,t2),max(t2,t0))
*   ..initialisation
      data ncall/0/, rxg/1.0e6_R8/, vis_par/0.0_R8/, diagno/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2npmo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2npmo_rxg', rxg)
        call xertst (0.0_R8.lt.rxg, 'faulty internal parameter rxg')
        call ipgeti ('b2npmo_diagno', diagno)
        call ipgetr ('b2tfhe_vis_par', vis_par)
        call xertst (0.0_R8.le.vis_par, 'faulty input vis_par')
        call ipgeti ('b2npmo_modvis', modvis)                        !srv 19.01.05
        call ipgeti ('b2npmo_iout', iout)                            !srv 04.07.08
        call ipgeti ('b2npmo_b2sifr_', b2sifr_ver)                   !srv 04.07.08
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..test rxf
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8, 
     . 'faulty argument rxf')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of na, te, ti
       do is = 0, ns-1
        call b2xvsg (n2, na(-1,-1,is), 1, 'na', '.gt.')
       enddo
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test edge values of fna_fcor
       do is = 0, ns-1
        call b2xvff (nx, ny, fna_fcor(-1,-1,0,is), 'fna_fcor')
       enddo
*    ..test sign of cvsa
       call b2xvsg (n2*2*ns, cvsa, 1, 'cvsa', '.ge.')
*    ..test edge values of cvsa
       do is = 0, ns-1
        call b2xvff (nx, ny, cvsa(-1,-1,0,is), 'cvsa')
       enddo
*    ..test sign of smo, smq
       do is = 0, ns-1
        call b2xvsg (n2, smo(-1,-1,1,is), 1, 'smo 1', '.le.')
        call b2xvsg (n2, smo(-1,-1,3,is), 1, 'smo 3', '.le.')
        call b2xvsg (n2, smq(-1,-1,1,is), 1, 'smq 1', '.le.')
        call b2xvsg (n2, smq(-1,-1,3,is), 1, 'smq 3', '.le.')
       enddo
      endif
*   ..compute ne, ne2, ro, pr
      call b2xpni (nx, ny, ns, na, ni)
      call b2xpne (nx, ny, ns, rza, na, ne)
      call b2xpne (nx, ny, ns, rz2, na, ne2)
      call b2xpro (nx, ny, ns, mp, am, na, ro)
      call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
*   ..compute ehxp
      call b2nxfx (nx, ny, qe, vol, hx, bb, ne, te, po, 
     &  ehxp, wrk0, wrk1)
*   ..compute smqdu
      call b2nxdu (nx, ny, ns, na, ua, smq, smqdu, wrk0, wrk1)
      if (iout.ne.0) then
       call my_out(70,nx,ny,smqdu,'b2npmo_smqdu')
      endif
*   ..initialise
      call sfill (n2*4, 0.0_R8, smt, 1)
      call sfill (n2*2, 0.0_R8, flct, 1)
      call sfill (n2*2, 0.0_R8, cvct, 1)
      call sfill (n2, 0.0_R8, resmt, 1)
      call sfill (n2, 0.0e0_R8, ctcft, 1)                              !srv 06.04.07

* ..main computation
*   ..compute individual species corrections
      do is = 0, ns-1
*    ..compute rob, pb
       call scopy (n2, na(-1,-1,is), 1, rob, 1)
       call sscal (n2, am(is)*mp, rob, 1)
       call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)
*    ..compute extended friction term
       if(b2sifr_ver.eq.0) then                                       !srv 26.10.08 {
        call b2sifr (nx, ny, ns, is, vol, hx, hz, bb, na, ua, te, ti, !srv 09.01.01
     &   ni, ne, ne2, ehxp, lnlam, cthe, cthi, smbch, smbtf,       !srv 11.09.09
     &   wrk0, wrk1, wrk2)                                             !srv 02.07.08 01.07.09
       else
        call b2sifr_ (nx, ny, ns, is, hx, hz, vol, bb, na, ua, te, ti, !srv 09.01.01 18.10.02
     &   ni, ne, ne2, ehxp, lnlam, cthe, cthi, smbch, smbtf,       !srv 11.09.09
     &   wrk0, wrk1, wrk2)                                             !srv 01.07.09
       endif                                                          !srv 26.10.08 }
*    ..compute pressure gradient term
       call b2sigp (nx, ny, is, rza(-1,-1,is), vol, hx, hz, bb,  
     &   na(-1,-1,is), ti, ne, te, po, smbgp, wrk0, wrk1)              !srv 02.07.08
*    ..compute flub
       call scopy (n2*2, fna_fcor(-1,-1,0,is), 1, flub, 1)
       call sscal (n2*2, am(is)*mp, flub, 1)
!       call b2nxfm (nx, ny, am(is)*mp, fna_fcor(-1,-1,0,is), flub)
       if (modvis.eq.1) then                                           !srv 19.01.05 11.09.09 {
*    ..add convective term due to new form of ion viscosity term only
        call b2nxfv (nx, ny, is, hz, cvsahz(-1,-1,0,is), flubv)        !srv 01.07.05 11.09.09
        call saxpy (n2*2, 1.0_R8, flubv, 1, flub, 1)                   !srv 11.09.09
       endif                                                           !srv 19.01.05 }
*    ..compute flcb, cvcb
       call b2nxfc (nx, ny, is, vol, pbshz, rob, ua(-1,-1,is), flub,
     &   cvsahz(-1,-1,0,is), flcb, cvcb, wrk0)                         !srv 09.01.01
*    ..compute additional viscosity term
        call b2siav (nx, ny, ns, is, modvis, am(is)*mp, hx, hy, hz, vol, 
     &   bb, ua(-1,-1,is), wadia(-1,-1,0,is), vaecrb(-1,-1,0,is),  
     &   te, ti, na(-1,-1,is), po, ni(-1,-1,1), ne2, lnlam, chci,      !srv 01.07.09
     &   cvsa(-1,-1,0,is), cvsa_cl(-1,-1,0,is), cvsahz(-1,-1,0,is), 
     &   cvsahz_cl(-1,-1,0,is), smbvh, smbvv, wrk0, wrk1, wrk2)         !srv 23.07.99 04.10.99 17.01.01 09.01.01 31.07.02 01.07.05 
*    ..compute centrifugial force term                                 !srv 09.07.01
       call b2sicf (nx, ny, is, am(is)*mp, bb, vol, hx, hz, 
     *     na(-1,-1,is), ua(-1,-1,is), smbcf, ctcfb, wrk0)              !srv 09.07.01 31.07.02 22.07.05} 06.04.07
*    ..compute total momentum source
       call scopy (n2*4, smo(-1,-1,0,is), 1, smb, 1)
       call saxpy (n2*4, 1.0_R8, smq(-1,-1,0,is), 1, smb, 1)
       b2npmo_smaf(:,:,:,is)=smbch(:,:,:)+smbtf(:,:,:)            !srv 11.09.09 {
       call saxpy (n2*4, 1.0_R8, smbch, 1, smb, 1)
       call saxpy (n2*4, 1.0_R8, smbtf, 1, smb, 1)                  !srv 11.09.09 }
       b2npmo_smag(:,:,:,is)=smbgp(:,:,:)
       call saxpy (n2*4, 1.0_R8, smbgp, 1, smb, 1)
       b2npmo_smav(:,:,:,is)=smbvh(:,:,:)+smbvv(:,:,:)                 !srv 03.07.08 {
       call saxpy (n2*4, 1.0_R8, smbvh,  1, smb, 1)
       call saxpy (n2*4, 1.0_R8, smbvv, 1, smb, 1)
       call saxpy (n2*4, 1.0_R8, smbcf,  1, smb, 1)                    !srv 03.07.08 }
       call b2nxdv (nx, ny, rob, smq(-1,-1,0,is), smqdu, smb)         !srv ???
*    ..save final source for species (is)                             !srv 10.11.02
       smof(:,:,:,is) = smb                                           !srv 10.11.02
*    ..compute residual
       call b2urmo (nx, ny, is, ua(-1,-1,is), rob, smb,
     &   flub, cvsahz(-1,-1,0,is), resmo(-1,-1,is), wrkf)           !srv 04.07.08
       fmo(:,:,:,is)=wrkf
       write (chns,'(i2.2)') is                                     !srv 18.03.10
       if (iout.ne.0) then                                          !srv 12.07.02 {
        call my_out(70,nx,ny,ua(-1,-1,is),'b2npmo_ua_input'//chns)
        call my_out(70,nx,ny,resmo(-1,-1,is),'b2npmo_resmo'//chns)
        call my_out(70,nx,ny,fmo(-1,-1,0,is),'b2npmo_fmox'//chns)
        call my_out(70,nx,ny,fmo(-1,-1,1,is),'b2npmo_fmoy'//chns)
        if(.not. is_neutral(is))
     &   call my_out(70,nx,ny,flubv(-1,-1,0),'b2npmo_flubvx'//chns)
!	   call my_out(70,nx,ny,flubv(-1,-1,1),'b2npmo_flubvy'//chns)
!	   do k=0,3
!	    write (chk,'(i1)') k  
!	    call my_out(70,nx,ny,smb(-1,-1,k),'b2npmo_smo'//chk//chns)
!	   enddo
        call my_out(70,nx,ny,smb(:,:,0)+smb(:,:,1)*ua(:,:,is)+
     &        smb(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smb(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smo'//chns)
        call my_out(70,nx,ny,smo(:,:,0,is)+smo(:,:,1,is)*ua(:,:,is)+
     &        smo(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smo(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smb'//chns)
        call my_out(70,nx,ny,smbch(:,:,0)+smbch(:,:,1)*ua(:,:,is)+
     &        smbch(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbch(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smoch'//chns)
        call my_out(70,nx,ny,smbtf(:,:,0)+smbtf(:,:,1)*ua(:,:,is)+
     &        smbtf(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbtf(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smotf'//chns)
        call my_out(70,nx,ny,smbgp(:,:,0)+smbgp(:,:,1)*ua(:,:,is)+
     &        smbgp(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbgp(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smogp'//chns)
        call my_out(70,nx,ny,smbvv(:,:,0)+smbvv(:,:,1)*ua(:,:,is)+
     &        smbvv(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbvv(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smovv'//chns)
        call my_out(70,nx,ny,smbvh(:,:,0)+smbvh(:,:,1)*ua(:,:,is)+
     &        smbvh(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbvh(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smovh'//chns)
        call my_out(70,nx,ny,smbcf(:,:,0)+smbcf(:,:,1)*ua(:,:,is)+
     &        smbcf(:,:,2)*na(:,:,is)*mp*am(is)+
     &        smbcf(:,:,3)*ua(:,:,is)*na(:,:,is)*mp*am(is),
     &        'b2npmo_smocf'//chns)
       endif                                                           !srv 12.07.02 }
       if (solving.and.all(solvemo(is,0:nnreg(0)))) then                !srv 06.07.10
         solvedum(0:nnreg(0))=solvemo(is,0:nnreg(0))
*    ..compute correction
         call b2usmo (nx, ny, is, nnreg(0), solvedum, itcnt, rxg,        !srv 04.07.08
     &   vol, hx, hz, bb, rob, pb, ua(-1,-1,is),                        !srv 09.01.01
     &   smb, flcb, cvcb, resmo(-1,-1,is), ctcfb, corua(-1,-1,is),      !srv 06.04.07
     &   pcca(-1,-1,0,is), aa5, wrk1, 'b2usmo_'//chns)
csw 21feb2014 added this else branch
       else
        pcca(:,:,:,is)  = 0.0_R8
        resmo(:,:,is) = 0.0_R8
        corua(:,:,is) = 0.0_R8
       endif
*    ..compute contribution to total residual
*      (now the friction cancels out)

csw 21feb2014 and do this only for the solved species, i.e. added next branch
       if (solving.and.all(solvemo(is,0:nnreg(0)))) then 
        call scopy (n2*4, smo(-1,-1,0,is), 1, smb, 1)
        call saxpy (n2*4, 1.0_R8, smbgp, 1, smb, 1)
        call saxpy (n2*4, 1.0_R8, smbvh, 1, smb, 1)
        call saxpy (n2*4, 1.0_R8, smbvv, 1, smb, 1)                     !srv 03.07.08
        call saxpy (n2*4, 1.0_R8, smbcf, 1, smb, 1)                     !srv 09.07.01
        call b2urmo (nx, ny, is, ua(-1,-1,is), rob, smb,                !srv 09.01.01
     &   flub, cvsahz(-1,-1,0,is), wrk0, wrkf)
*    ..accumulate residual
        call saxpy (n2, 1.0_R8, wrk0, 1, resmt, 1)
*    ..accumulate smb
        call b2nxst (nx, ny, rob, smb, ro, smt)
*    ..accumulate transport coefficients
        call saxpy (n2*2, 1.0_R8, flcb, 1, flct, 1)
        call saxpy (n2*2, 1.0_R8, cvcb, 1, cvct, 1)
        call saxpy (n2, 1.0_R8, ctcfb, 1, ctcft, 1)                  !srv 06.04.07
       endif

      enddo
csw 21feb2014, execute simultaneous corrections always
c      if (solving.and.all(solvemo(0:ns-1,0:nnreg(0)))) then           !srv 06.07.10
      if (solving) then    
*   ..compute simultaneous correction
       call sfill (n2, 0.0_R8, wrk0, 1)
       solvedum(0:nnreg(0))=.true.
csw 14feb2014 warning: variable 'is' unspecified here, (but no use of 'is' in b2usmo)
       call b2usmo (nx, ny, is, nnreg(0), solvedum, itcnt, rxg,       !srv 04.07.08
     &  vol, hx, hz, bb, ro, pr, wrk0, smt, flct, cvct,               !srv 09.01.01
     &  resmt, ctcft, corut, pccm, aa5, wrk1, 'b2usmo_total')         !srv 06.04.07
      else
       pccm  = 0.0_R8
       resmt = 0.0_R8
       corut = 0.0_R8
      endif                                                           !srv 06.07.10
*   ..apply correction
       do is = 0, ns-1
        do iy = -1, ny
         do ix = -1, nx
          ireg = region(ix,iy,0)
          if (solvemo(is,ireg)) then
            if((leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) .or.
     &       (bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) ) then !srv 27.01.11
               ua(ix,iy,is) = ua(ix,iy,is) + corua(ix,iy,is)
            else
               ua(ix,iy,is) = ua(ix,iy,is)+
     &         rxf*med(corua(ix,iy,is),corut(ix,iy),
     &         corua(ix,iy,is)+corut(ix,iy))
            endif
          endif
         enddo
        enddo
       enddo
*   ..correct pcca
       call b2nxcm (nx, ny, ns, pccm, pcca)
*   ..set error parameter
*     (at present all errors cause an abort through xerrab)
      ierr = 0
      if (iout.ne.0) then                                          !srv 12.07.02 {
       do is = 0, ns-1
        write (chns,'(i2.2)') is                                    !srv 18.05.02
         call my_out(70,nx,ny,ua(-1,-1,is),'b2npmo_ua'//chns)
!	   call my_out(70,nx,ny,pcca(-1,-1,0,is),'b2npmo_pccax'//chns)
!	   call my_out(70,nx,ny,pcca(-1,-1,1,is),'b2npmo_pccay'//chns)
!	   call my_out(70,nx,ny,fna_fcor(-1,-1,0,is),'b2npmo_fna_fcorx'//chns)
!	   call my_out(70,nx,ny,fna_fcor(-1,-1,1,is),'b2npmo_fna_fcory'//chns)
       enddo
       call my_out(70,nx,ny,pccm(-1,-1,0),'b2npmo_pccmx')          !srv 11.01.13
       call my_out(70,nx,ny,pccm(-1,-1,1),'b2npmo_pccmy')          !srv 11.01.13
      endif                                                   !srv 12.07.02 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2npmo

      end subroutine b2npmo
