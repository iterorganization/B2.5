*-----------------------------------------------------------------------
*.specification

      subroutine b2trno (nx, ny, ns, nscx, nscxmax, iscx, 
     &  bb, vol, gs, hx, hy, hz, qz, pbs,                            !srv 09.01.01
     &  na, ua, te, ti, po, ne, ni, ne2, cdna, cdpa, cvla, cvsa,
     &  cdnahz, cdpahz, cvlahz, cvmahz, cvsahz,                      !srv 09.01.01 15.12.05
     &  chce, chve, chci, chvi, csig, calf,
     4	fllim0fhi)
      use b2mod_types
      use b2mod_anomalous_transport
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, nscxmax, iscx(0:nscxmax-1)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), 
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), 
     &  pbs(-1:nx,-1:ny,0:1), hz(-1:nx,-1:ny),                       !srv 09.01.01
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), ne2(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cdna(-1:nx,-1:ny,0:1,0:ns-1), cdpa(-1:nx,-1:ny,0:1,0:ns-1),
     &  cvla(-1:nx,-1:ny,0:1,0:ns-1), cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  cdnahz(-1:nx,-1:ny,0:1,0:ns-1), cdpahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvlahz(-1:nx,-1:ny,0:1,0:ns-1), cvsahz(-1:nx,-1:ny,0:1,0:ns-1), !srv 09.01.01
     &  cvmahz(-1:nx,-1:ny,0:1,0:ns-1),                                 !srv 15.12.05
     &  chce(-1:nx,-1:ny,0:1), chve(-1:nx,-1:ny,0:1),
     &  chci(-1:nx,-1:ny,0:1), chvi(-1:nx,-1:ny,0:1),
     &  csig(-1:nx,-1:ny,0:1), calf(-1:nx,-1:ny,0:1),
     1	fllim0fhi(-1:nx,-1:ny,0:1,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TRNO computes the anomalous transport coefficients, including
*     geometric factors such as vol/hx**2 and vol/hy**2.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the following discrete anomalous transport
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
*     - calf: thermo-electric coefficient.
*
*     The local differential transport coefficients, defined on cells,
*     are obtained by a calls to b2tqna. To compute from these the
*     discrete transport coefficients involves multiplication by
*     geometric factors (vol, vol/hx**2, vol/hy**2), and, usually,
*     interpolation to cell faces.
*
*     We treat the anomalous electrical conductivity separately, and
*     do not add it to the classical conductivity. See treatment in
*     b2tfch.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is
      integer hcimx_flag, mode
	integer, save :: iout = 0,                                      !srv 17.06.08
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npp7 = 0, iout_b2nppo = 0, iout_b2npmo =0
	character(2) chns                                            !srv 17.06.08
	integer, save :: flux_limit_to_dpa = 0, csig_an_style = 0      !srv 15.06.08
      real (kind=R8) ::
     *  wrk(-1:nx,-1:ny), hcimx(-1:nx,-1:ny,0:1,0:ns-1),
     1	pol_anom_scale, 
     2  work(-1:nx,-1:ny), dpa0mx(-1:nx,-1:ny,0:1,0:ns-1)            !srv 09.07.03
      save ncall, pol_anom_scale, hcimx_flag, mode
*   ..procedures
      real (kind=R8) ::
     *  smin, smax
      external subini, subend, xertst, sfill, smin, smax, scopy
      external ipgetr, ipgeti
      external b2tqna, b2txcx, b2txcy, b2txcv, b2txsx, b2txsy
*   ..initialisation
      data ncall/0/, pol_anom_scale/1.0_R8/, hcimx_flag/1/, mode/1/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2trno')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..read data on first call
      if (ncall.eq.0) then
        call ipgetr ('b2trno_pol_anom_scale', pol_anom_scale)
        call ipgeti ('b2trno_hcimx_flag', hcimx_flag)
        call ipgeti ('b2tral_mode', mode)
        call ipgeti ('b2trno_flux_limit_to_dpa', flux_limit_to_dpa)
	  call ipgeti ('b2trno_csig_an_style',csig_an_style)           !srv 15.06.08
        call ipgeti ('b2trno_iout', iout)                            !srv 17.06.08                    
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 19.05.10 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)                     !srv 19.05.10 }
      end if
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test sign of pbs
!       call b2xvsg (n2*2, pbs, 1, 'pbs', '.ge.')                    !srv 11.01.13
*    ..test edge values of pbs
       call b2xvff (nx, ny, pbs, 'pbs')
*    ..test range of qz
       call xertst (-1.0_R8.lt.smin(n2,qz(-1,-1,0),1).and.
     &              smax(n2,qz(-1,-1,0),1).lt.1.0_R8, 
     .              'faulty argument range qz(,,0)')
       call xertst (0.0_R8.lt.smin(n2,qz(-1,-1,1),1).and.
     &              smax(n2,qz(-1,-1,1),1).le.1.0_R8, 
     .              'faulty argument range qz(,,1)')
*    ..test sign of na, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

* ..compute transport coefficients
*   ..compute differential anomalous transport coefficients
      call b2tqna (nx, ny, ns, nscx, iscx, bb, na, te, ti, 
     & ni, ne, ne2, hx, hy, qz, dna0, dpa0, vla0, vma0,             !srv 15.12.05
     & vsa0, hci0, hcib, hce0, sig0, alf0, wrk)
*   ..apply flux limit to neutrals contribution to hci
      call b2tlh0 (nx, ny, ns, hx, hy, na, ti, hcib, hcimx, fllim0fhi)      
      if(hcimx_flag.eq.1) then
         hcib(:,:,:)=hcimx(:,:,1,:)
         do ix=-1,nx
            do iy=-1,ny
               hci0(ix,iy)=0.0_R8
               do is=0,ns-1
                  hci0(ix,iy)=hci0(ix,iy)+hcib(ix,iy,is)
               enddo
            enddo
         enddo
      endif
	if(flux_limit_to_dpa .ne. 0) then                              !srv 15.06.08 {
*   .. apply flux limit to neutrals contribution to dpa0
       call b2tlc0 (nx, ny, ns, hx, hy, na, te, ti, dpa0, dpa0mx)    !srv 09.07.03 
	endif                                                          !srv 15.06.08 }
*   ..incorporate geometry
      do is = 0, ns-1
       call b2txcx (nx, ny, mode, vol, gs, hx, qz, dna0(-1,-1,is),
     &   cdna(-1,-1,0,is))
       call b2txcy (nx, ny, mode, vol, gs, hy, qz, dna0(-1,-1,is),
     &   cdna(-1,-1,1,is))
	 if(flux_limit_to_dpa .ne. 0) then                             !srv 15.06.08
        call b2txcx (nx, ny, mode, vol, gs, hx, qz, dpa0mx(-1,-1,0,is),
     &   cdpa(-1,-1,0,is))
        call b2txcy (nx, ny, mode, vol, gs, hy, qz, dpa0mx(-1,-1,1,is),
     &   cdpa(-1,-1,1,is))
	 else                                                          !srv 15.06.08 {
        call b2txcx (nx, ny, mode, vol, gs, hx, qz, dpa0(-1,-1,is),
     &   cdpa(-1,-1,0,is))
        call b2txcy (nx, ny, mode, vol, gs, hy, qz, dpa0(-1,-1,is),
     &   cdpa(-1,-1,1,is))
	 endif                                                         !srv 15.06.08 }
       call b2txsx (nx, ny, vol, hx, qz, vla0(-1,-1,0,is),
     &   cvla(-1,-1,0,is))
       call b2txsy (nx, ny, vol, hy, qz, vla0(-1,-1,1,is),
     &   cvla(-1,-1,1,is))
       call b2txcx (nx, ny, mode, vol, gs, hx, qz, vsa0(-1,-1,is),
     &   cvsa(-1,-1,0,is))
       call b2txcy (nx, ny, mode, vol, gs, hy, qz, vsa0(-1,-1,is),
     &   cvsa(-1,-1,1,is))
      enddo
*   ..incorporate geometry and hz                                     !srv 09.01.01
      work = hz*vol
      do is = 0, ns-1
       call b2txcx (nx, ny, mode, work, gs, hx, qz, dna0(-1,-1,is),
     &   cdnahz(-1,-1,0,is))
       call b2txcy (nx, ny, mode, work, gs, hy, qz, dna0(-1,-1,is),
     &   cdnahz(-1,-1,1,is))
	 if(flux_limit_to_dpa .ne. 0) then                             !srv 15.06.08
        call b2txcx (nx, ny, mode, work, gs, hx, qz, dpa0mx(-1,-1,0,is),
     &   cdpahz(-1,-1,0,is))                                         !srv 09.07.03 {
        call b2txcy (nx, ny, mode, work, gs, hy, qz, dpa0mx(-1,-1,1,is),
     &   cdpahz(-1,-1,1,is))                                         !srv 09.07.03 }
	 else                                                          !srv 15.06.08 {
        call b2txcx (nx, ny, mode, work, gs, hx, qz, dpa0(-1,-1,is),
     &   cdpahz(-1,-1,0,is))
        call b2txcy (nx, ny, mode, work, gs, hy, qz, dpa0(-1,-1,is),
     &   cdpahz(-1,-1,1,is))
	 endif                                                         !srv 15.06.08 }
       call b2txsx (nx, ny, work, hx, qz, vla0(-1,-1,0,is),
     &   cvlahz(-1,-1,0,is))
       call b2txsy (nx, ny, work, hy, qz, vla0(-1,-1,1,is),
     &   cvlahz(-1,-1,1,is))
       call b2txsx (nx, ny, work, hx, qz, vma0(-1,-1,0,is),         !srv 15.12.05
     &   cvmahz(-1,-1,0,is))
       call b2txsy (nx, ny, work, hy, qz, vma0(-1,-1,1,is),         !srv 15.12.05
     &   cvmahz(-1,-1,1,is))
       call b2txcx (nx, ny, mode, work, gs, hx, qz, vsa0(-1,-1,is), !srv 01.07.05 {
     &   cvsahz(-1,-1,0,is))
       call b2txcy (nx, ny, mode, work, gs, hy, qz, vsa0(-1,-1,is),
     &   cvsahz(-1,-1,1,is))                                        !srv 01.07.05 }
      enddo                                                         !srv 09.01.01
*   ..save the anomalous transport coefficients 					  !srv 04.10.99
      cvsa_an = cvsa                                                !srv 04.10.99
      call b2txcx (nx, ny, mode, vol, gs, hx, qz, hce0, chce(-1,-1,0))
      call b2txcy (nx, ny, mode, vol, gs, hy, qz, hce0, chce(-1,-1,1))
      call b2txfx (nx, ny, ns, mode, vol, gs, hx, qz, hcimx, 
     & chci(-1,-1,0))
      call b2txfy (nx, ny, ns, mode, vol, gs, hy, qz, hcimx, 
     & chci(-1,-1,1))
      call b2txcx (nx, ny, -1, vol, gs, hx, qz, alf0, calf(-1,-1,0))
      call b2txcy (nx, ny, -1, vol, gs, hy, qz, alf0, calf(-1,-1,1))
	if(csig_an_style. eq.0) then                                !srv 15.06.08
       call b2txcx (nx, ny, mode, vol, gs, hx, qz, sig0, 
     &      csig_an(-1,-1,0))
       call b2txcy (nx, ny, mode, vol, gs, hy, qz, sig0, 
     &      csig_an(-1,-1,1))
*   .. set unused contribution to the electrical conductivity to zero
       call sfill (n2*2, 0.0_R8, csig, 1)
	else 					                                     !srv
       call b2txcx (nx, ny, mode, vol, gs, hx, qz, sig0, csig(-1,-1,0))
       call b2txcy (nx, ny, mode, vol, gs, hy, qz, sig0, csig(-1,-1,1))
*   .. save the anomalous transport coefficients
       csig_an = csig 
*   .. set csigx = 0 and calfx = 0
	 call sfill (n2, 0.0e0_R8, csig(-1,-1,0), 1)
       call sfill (n2, 0.0e0_R8, calf(-1,-1,0), 1)
	endif					                                     !srv
*   ..set strange velocities chve, chvi
      call sfill (n2*2, 0.0_R8, chve, 1)
      call sfill (n2*2, 0.0_R8, chvi, 1)
*   ..scale parallel coefficients if desired
      if (pol_anom_scale.ne.1.0_R8) then
        cdna(:,:,0,:) = cdna(:,:,0,:) * pol_anom_scale
        cdpa(:,:,0,:) = cdpa(:,:,0,:) * pol_anom_scale
        cvla(:,:,0,:) = cvla(:,:,0,:) * pol_anom_scale
        cvsa(:,:,0,:) = cvsa(:,:,0,:) * pol_anom_scale
        chce(:,:,0) = chce(:,:,0) * pol_anom_scale
        chci(:,:,0) = chci(:,:,0) * pol_anom_scale
        calf(:,:,0) = calf(:,:,0) * pol_anom_scale
        csig_an(:,:,0) = csig_an(:,:,0) * pol_anom_scale
      end if

      if(iout_b2npc9.ne.0 .or. iout_b2npco.ne.0) then                !srv 19.05.10 {
       do is = 0, ns-1
	  write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'b2trno_cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'b2trno_cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'b2trno_cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'b2trno_cdpay'//chns)
!        call my_out(70,nx,ny,cvla(-1,-1,0,is),'b2trno_cvlax'//chns)
!        call my_out(70,nx,ny,cvla(-1,-1,1,is),'b2trno_cvlay'//chns)
       enddo
	endif
	if(iout_b2npmo.ne.0) then
       do is = 0, ns-1
	  write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'b2trno_cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'b2trno_cvsay'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'b2trno_cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'b2trno_cvsahzy'//chns)
       enddo
	endif
	if(iout_b2nph9.ne.0 .or. iout_b2npht.ne.0) then
        call my_out(70,nx,ny,chce(-1,-1,0),'b2trno_chcex')
        call my_out(70,nx,ny,chce(-1,-1,1),'b2trno_chcey')
        call my_out(70,nx,ny,chci(-1,-1,0),'b2trno_chcix')
        call my_out(70,nx,ny,chci(-1,-1,1),'b2trno_chciy')
!        call my_out(70,nx,ny,chve(-1,-1,0),'b2trno_chvex')
!        call my_out(70,nx,ny,chve(-1,-1,1),'b2trno_chvey')
!        call my_out(70,nx,ny,chvi(-1,-1,0),'b2trno_chvix')
!        call my_out(70,nx,ny,chvi(-1,-1,1),'b2trno_chviy')
	endif
	if(iout_b2npp7.ne.0 .or. iout_b2nppo.ne.0) then
        call my_out(70,nx,ny,csig(-1,-1,0),'b2trno_csigx')
        call my_out(70,nx,ny,csig(-1,-1,1),'b2trno_csigy')
	endif                                                          !srv 19.05.10 }
      if(iout.ne.0) then                                        !srv 17.06.08 {
       do is = 0, ns-1
	  write (chns,'(i2.2)') is 
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'b2trno_cdnax'//chns)
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'b2trno_cdnay'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'b2trno_cdpax'//chns)
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'b2trno_cdpay'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,0,is),'b2trno_cvlax'//chns)
        call my_out(70,nx,ny,cvla(-1,-1,1,is),'b2trno_cvlay'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'b2trno_cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'b2trno_cvsay'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,0,is),'b2trno_dnahzx'//chns)
        call my_out(70,nx,ny,cdnahz(-1,-1,1,is),'b2trno_dnahzy'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,0,is),'b2trno_cdpahzx'//chns)
        call my_out(70,nx,ny,cdpahz(-1,-1,1,is),'b2trno_cdpahzy'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,0,is),'b2trno_cvlahzx'//chns)
        call my_out(70,nx,ny,cvlahz(-1,-1,1,is),'b2trno_cvlahzy'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,0,is),'b2trno_cvsahzx'//chns)
        call my_out(70,nx,ny,cvsahz(-1,-1,1,is),'b2trno_cvsahzy'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,0,is),'b2trno_cvmahzx'//chns)
        call my_out(70,nx,ny,cvmahz(-1,-1,1,is),'b2trno_cvmahzy'//chns)
      enddo
        call my_out(70,nx,ny,chce(-1,-1,0),'b2trno_chcex')
        call my_out(70,nx,ny,chce(-1,-1,1),'b2trno_chcey')
        call my_out(70,nx,ny,chci(-1,-1,0),'b2trno_chcix')
        call my_out(70,nx,ny,chci(-1,-1,1),'b2trno_chciy')
        call my_out(70,nx,ny,chve(-1,-1,0),'b2trno_chvex')
        call my_out(70,nx,ny,chve(-1,-1,1),'b2trno_chvey')
        call my_out(70,nx,ny,chvi(-1,-1,0),'b2trno_chvix')
        call my_out(70,nx,ny,chvi(-1,-1,1),'b2trno_chviy')
        call my_out(70,nx,ny,csig(-1,-1,0),'b2trno_csigx')
        call my_out(70,nx,ny,csig(-1,-1,1),'b2trno_csigy')
        call my_out(70,nx,ny,calf(-1,-1,0),'b2trno_calfx')
        call my_out(70,nx,ny,calf(-1,-1,1),'b2trno_calfy')
        call my_out(70,nx,ny,csig_an(-1,-1,0),'b2trno_csig_anx')
        call my_out(70,nx,ny,csig_an(-1,-1,1),'b2trno_csig_any')
	endif                                                          !srv 17.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.scribbles

*!!   Consider to apply a flux limit to the anomalous particle transport.
*!!   Still must consider how to deal with a non-orthogonal grid.

*-----------------------------------------------------------------------
*.end b2trno

      end

