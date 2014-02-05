*-----------------------------------------------------------------------
*.specification

      subroutine b2tral (nx, ny, ns, nscx, nscxmax, iscx, ismain,    !srv 23.11.10
     1 bb, conn, vol, gs, hx, hy, hz, qz, qc, pbs, bzb, lnlam,       !srv 09.01.01 01.07.09
     1 fch_p, na, ua, te, ti, po, ne, ni, ne2, chvemx, chvimx,       !srv 03.06.03
     2 cdna, cdpa, cddi, cvla, cvsa, chce, chve, chci,               !srv 13.10.06
     3 chvi, csig, csigin, calf, cthe, cthi, 
     4 cdnahz, cdpahz, cvlahz, cvmahz, cvsahz, cvsa_cl, cvsahz_cl,   !srv 09.01.01 15.12.05
     4 vsal, fllim0fhi, fllimvisc)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_anomalous_transport
      use b2mod_rates, only : rza                                    !srv 16.06.08
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, nscxmax, iscx(0:nscxmax-1), ismain    !srv 23.11.10
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), 
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), hz(-1:nx,-1:ny),            !srv 09.01.01
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1), bzb(-1:nx,-1:ny),  !srv 13.10.06
     &  lnlam(-1:nx,-1:ny),                                           !srv 01.07.09
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne2(-1:nx,-1:ny),
     &  fch_p(-1:nx,-1:ny,0:1), conn(-1:ny), qc(-1:nx,-1:ny),         !srv 10.06.08
     &  chvemx(-1:nx,-1:ny), chvimx(-1:nx,-1:ny)                      !srv 03.06.03
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cdna(-1:nx,-1:ny,0:1,0:ns-1), cdpa(-1:nx,-1:ny,0:1,0:ns-1),
     *  cddi(-1:nx,-1:ny,0:1,0:ns-1),                                   !srv 13.10.06
     &  cvla(-1:nx,-1:ny,0:1,0:ns-1), cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  cdnahz(-1:nx,-1:ny,0:1,0:ns-1), cdpahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvlahz(-1:nx,-1:ny,0:1,0:ns-1), cvsahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvsahz_cl(-1:nx,-1:ny,0:1,0:ns-1),                              !srv 09.01.01
     &  cvmahz(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 15.12.05
     &  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     &  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  csig(-1:nx,-1:ny,0:1), csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),    !xpb
     &  calf(-1:nx,-1:ny,0:1),
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1),
     1  cvsa_cl(-1:nx,-1:ny,0:1,0:ns-1), vsal(-1:nx,-1:ny,0:1,0:ns-1), 
     1  fllim0fhi(-1:nx,-1:ny,0:1,0:ns-1), 
     2  fllimvisc(-1:nx,-1:ny,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRAL computes the transport coefficients, including geometric
*     factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete transport
*     coefficients.
*     - cdna: particle transport proportional to grad na.
*     - cdpa: particle transport proportional to grad pa.
*     - cvla: anomalous velocity.
*     - cvsa: viscosity.
*     - chce: electron thermal conductivity.
*     - chve: electron thermal strange velocity.
*     - chci: atom thermal conductivity.
*     - chvi: atom thermal strange velocity.
*     - csig: electrical conductivity.
*     - csigin: ion-neutral conductivity.
*     - calf: thermo-electric coefficient.
*
*     The routine operates as a simple interface to routines b2trcl
*     and b2trno, for the classical and the anomalous contributions,
*     respectively.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ncall, ix, iy, is, is0
      integer, save :: iout = 0, model = 1,                             !srv 17.06.08
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0
      character(2) chns, chns0                                          !srv 17.06.08
      real (kind=R8) ::
     1  chvi0(-1:nx,-1:ny,0:1), csig0(-1:nx,-1:ny,0:1),
     2  calf0(-1:nx,-1:ny,0:1), chci0(-1:nx,-1:ny,0:1),
     3  cvsa0(-1:nx,-1:ny,0:1,0:ns-1), cvsahz0(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     4  chce0(-1:nx,-1:ny,0:1), chve0(-1:nx,-1:ny,0:1),
     1  visc(-1:nx,-1:ny,0:ns-1), pol_anom_scale,
     *    wrk(-1:nx,-1:ny), wrkmin
      character*72 message
      real (kind=R8), save :: addvis1=0.0_R8                            !srv 31.07.02
      save ncall, pol_anom_scale
*   ..procedures
      external subini, subend, xertst, ipgetr, saxpy
      external b2trcl, b2trno
*   ..initialisation
      data ncall/0/, pol_anom_scale/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tral')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall.eq.0) then
        call ipgetr ('b2trno_pol_anom_scale', pol_anom_scale)
        call ipgetr ('b2siav_addvis1', addvis1)                      !srv 31.07.02
        call ipgeti ('b2tral_iout', iout)                            !srv 17.06.08                    
        call ipgeti ('b2tqca_model', model)
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 11.09.09 }
      end if

* ..compute transport coefficients
*   ..compute anomalous terms
      call b2trno (nx, ny, ns, nscx, nscxmax, iscx, 
     &  bb, vol, gs, hx, hy, hz, qz, pbs, na, ua,                     !srv 09.01.01
     &  te, ti, po, ne, ni, ne2, cdna, cdpa, cvla, cvsa,
     &  cdnahz, cdpahz, cvlahz, cvmahz, cvsahz,                       !srv 09.01.01 15.12.05
     &  chce, chve, chci, chvi, csig, calf, fllim0fhi)
*   ..compute classical terms
      call b2trcl (nx, ny, ns, ismain, bb, conn, vol, gs, hx, hy, hz, !srv 23.11.10
     1  qz, qc, pbs, fch_p, na, ua, te, ti, po, ne, ni, ne2, lnlam,   !srv 01.07.09
     2  chvemx, chvimx, cvsahz0,                                      !srv 09.01.01 03.06.03
     3  cvsa0, chce0, chve0, chci0, chvi0, csig0, calf0, csigin,
     4  cthe, cthi, cvsa_cl, cvsahz_cl, vsal, fllimvisc)              !srv 09.01.01
*   ..add to previous contributions
      call saxpy (n2*2*ns, 1.0_R8, cvsahz0, 1, cvsahz, 1)              !sv 09.01.01
      call saxpy (n2*2*ns, 1.0_R8, cvsa0, 1, cvsa, 1)
      visc=vsa0*pol_anom_scale
      call saxpy (n2*ns, 1.0_R8, visc, 1, vsal(-1,-1,0,0), 1)
      call saxpy (n2*ns, 1.0_R8, vsa0(-1,-1,0), 1, vsal(-1,-1,1,0), 1)
      call saxpy (n2*2, 1.0_R8, chce0, 1, chce, 1)
      call saxpy (n2*2, 1.0_R8, chve0, 1, chve, 1)
      call saxpy (n2*2, 1.0_R8, chci0, 1, chci, 1)
      call saxpy (n2*2, 1.0_R8, chvi0, 1, chvi, 1)
      call saxpy (n2*2, 1.0_R8, csig0, 1, csig, 1)
      call saxpy (n2*2, 1.0_R8, calf0, 1, calf, 1)
*   ..multiply x-component cvsahz by 4/3 only base ion                 !srv 31.07.02 {
!      if (addvis1.ne.0.0_R8.and.model.eq.3) then 
!      if (addvis1.ne.0.0_R8) then 
      if (model.eq.3) then                                             !srv 20.10.09
       call sscal (n2, 4.0e0_R8/3.0e0_R8, cvsahz(-1,-1,0,ismain), 1)   !srv 23.11.10
      endif                                                            !srv 31.07.02 }

*   ..compute terms from diamagnetic parts                             !srv 13.10.06 {
      do is = 0, ns-1
       if(is_neutral(is)) then
        cddi(:,:,0,is) = 0.0_R8
        cddi(:,:,1,is) = 0.0_R8
       else
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
           cddi(ix,iy,0,is) = 0.0_R8
          else
           cddi(ix,iy,0,is) =
     &      +2.0e0_R8*(vol(leftix(ix,iy),leftiy(ix,iy))*
     &         bzb(leftix(ix,iy),leftiy(ix,iy))
     &         +vol(ix,iy)*bzb(ix,iy))/
     &       ((hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*
     &        (hy(leftix(ix,iy),leftiy(ix,iy))+hy(ix,iy)))/
     &        (rza(ix,iy,is)*qe)
          endif
         enddo
        enddo

        do iy = -1, ny
         do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) then
           cddi(ix,iy,1,is) = 0.0e0_R8
          else
           cddi(ix,iy,1,is) =
     &       2.0e0_R8*(vol(bottomix(ix,iy),bottomiy(ix,iy))*
     &         bzb(bottomix(ix,iy),bottomiy(ix,iy))
     &         +vol(ix,iy)*bzb(ix,iy))/
     &       ((hx(bottomix(ix,iy),bottomiy(ix,iy))+hx(ix,iy))*
     &        (hy(bottomix(ix,iy),bottomiy(ix,iy))+hy(ix,iy)))/
     &        (rza(ix,iy,is)*qe)
          endif
         enddo
        enddo

       endif 
      enddo                                                           !srv 13.10.06 }

*  ..check that we have some particle diffusion
      if (ncall.le.3) then                                            !srv 12.03.12 {
        do is = 0, ns-1
          wrk(:,:) = cdna(:,:,0,is)+cdpa(:,:,0,is)
          wrkmin = huge(1.0_R8)
          do iy = -1, ny
            do ix = -1, nx
              if (leftix(ix,iy).ne.-2) then
                wrkmin = min(wrkmin, wrk(ix,iy))
              endif 
            enddo
          enddo
          write(message,'(a, i2, a)') 
     .     'No particle diffusion for species ',is,
     .     ' in poloidal direction'
          call xertst (wrkmin.gt.0.0_R8 .or. pol_anom_scale.eq.0.0_R8,
     .     message)
          wrk(:,:) = cdna(:,:,1,is)+cdpa(:,:,1,is)
          wrkmin = huge(1.0_R8)
          do iy = -1, ny
            do ix = -1, nx
              if (bottomiy(ix,iy).ne.-2) then
                wrkmin = min(wrkmin, wrk(ix,iy))
              endif 
            enddo
          enddo
          write(message,'(a, i2, a)') 
     .     'No particle diffusion for species ',is,
     .     ' in radial direction'
          call xertst (wrkmin.gt.0.0_R8, message)
        enddo
      endif                                                           !srv 12.03.12 }

      if(iout.ne.0) then                                              !srv 17.06.08 {
       do is = 0, ns-1
        write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'cdpay'//chns)
        call my_out(70,nx,ny,cddi(-1,-1,0,is),'cddix'//chns)
        call my_out(70,nx,ny,cddi(-1,-1,1,is),'cddiy'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,0,is),'cvlax'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,1,is),'cvlay'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'cvsay'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,0,is),'cdnahzx'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,1,is),'cdnahzy'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,0,is),'cdpahzx'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,1,is),'cdpahzy'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,0,is),'cvlahzx'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,1,is),'cvlahzy'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'cvsahzy'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,0,is),'cvmahzx'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,1,is),'cvmahzy'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'cvsa_clx'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,1,is),'cvsa_cly'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,0,is),'cvsahz_clx'//chns)
        call my_out(70,nx,ny,cvsahz_cl(-1,-1,1,is),'cvsahz_cly'//chns)
        call my_out(70,nx,ny,cthe(-1,-1,is),'cthe'//chns)
        call my_out(70,nx,ny,cthi(-1,-1,is),'cthi'//chns)
        do is0=0,ns-1
         write (chns0,'(i2.2)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then 
          call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'csiginx'//chns//
     &        '_'//chns0)
          call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'csiginy'//chns//
     &        '_'//chns0)
         endif
        enddo
       enddo
       call my_out(70,nx,ny,chce(-1,-1,0),'chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'chcey')
       call my_out(70,nx,ny,chci(-1,-1,0),'chcix')
       call my_out(70,nx,ny,chci(-1,-1,1),'chciy')
       call my_out(70,nx,ny,chve(-1,-1,0),'chvex')
       call my_out(70,nx,ny,chve(-1,-1,1),'chvey')
       call my_out(70,nx,ny,chvi(-1,-1,0),'chvix')
       call my_out(70,nx,ny,chvi(-1,-1,1),'chviy')
       call my_out(70,nx,ny,csig(-1,-1,0),'csigx')
       call my_out(70,nx,ny,csig(-1,-1,1),'csigy')
       call my_out(70,nx,ny,calf(-1,-1,0),'calfx')
       call my_out(70,nx,ny,calf(-1,-1,1),'calfy')
      endif                                                          !srv 17.06.08 }

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0) then                !srv 11.09.09 {
       do is = 0, ns-1
        write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'b2tral_cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'b2tral_cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'b2tral_cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'b2tral_cdpay'//chns)
        if(.not.is_neutral(is)) then
         call my_out(70,nx,ny,cddi(-1,-1,0,is),'b2tral_cddix'//chns)
         call my_out(70,nx,ny,cddi(-1,-1,1,is),'b2tral_cddiy'//chns)
        endif
!        call my_out(70,nx,ny,cvla(-1,-1,0,is),'b2tral_cvlax'//chns)
!        call my_out(70,nx,ny,cvla(-1,-1,1,is),'b2tral_cvlay'//chns)
       enddo
      endif
        
      if(iout_b2npmo.ne.0) then
       do is = 0, ns-1
        write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'b2tral_cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'b2tral_cvsahzy'//chns)
!        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'b2tral_cvsa_clx'//
!     &              chns)
!        call my_out(70,nx,ny,cvsa_cl(-1,-1,1,is),'b2tral_cvsa_cly'//
!     &              chns)
!        call my_out(70,nx,ny,cvsahz_cl(-1,-1,0,is),'b2tral_cvsahz_clx'//
!     &             chns)
!        call my_out(70,nx,ny,cvsahz_cl(-1,-1,1,is),'b2tral_cvsahz_cly'//
!     &             chns)
        call my_out(70,nx,ny,cthe(-1,-1,is),'b2tral_cthe'//chns)
        call my_out(70,nx,ny,cthi(-1,-1,is),'b2tral_cthi'//chns)
       enddo
      endif
      if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
       call my_out(70,nx,ny,chce(-1,-1,0),'b2tral_chcex')
       call my_out(70,nx,ny,chce(-1,-1,1),'b2tral_chcey')
       call my_out(70,nx,ny,chci(-1,-1,0),'b2tral_chcix')
       call my_out(70,nx,ny,chci(-1,-1,1),'b2tral_chciy')
!       call my_out(70,nx,ny,chve(-1,-1,0),'b2tral_chvex')
!       call my_out(70,nx,ny,chve(-1,-1,1),'b2tral_chvey')
!       call my_out(70,nx,ny,chvi(-1,-1,0),'b2tral_chvix')
!       call my_out(70,nx,ny,chvi(-1,-1,1),'b2tral_chviy')
      endif
      if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
       do is = 0, ns-1
        write (chns,'(i2.2)') is 
        do is0=0,ns-1
         write (chns0,'(i2.2)') is0
         if(is0.ne.is .and. .not.is_neutral(is)) then 
          call my_out(70,nx,ny,csigin(-1,-1,0,is,is0),'b2tral_csiginx'
     &        //chns//'_'//chns0)
          call my_out(70,nx,ny,csigin(-1,-1,1,is,is0),'b2tral_csiginy'
     &        //chns//'_'//chns0)
         endif
        enddo
       enddo
       call my_out(70,nx,ny,csig(-1,-1,0),'b2tral_csigx')
       call my_out(70,nx,ny,csig(-1,-1,1),'b2tral_csigy')
      endif                                                          !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tral

      end subroutine b2tral
