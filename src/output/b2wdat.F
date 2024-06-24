      subroutine b2wdat (nx, ny, ns, ismain, ismain0)                    !srv 05.10.01
      use b2mod_types
      use b2mod_constants
      use b2mod_geo
      use b2mod_geo_corner
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_anomalous_transport
      use b2mod_transport
      use b2mod_math
      use b2mod_work
      use b2mod_indirect
      use b2mod_external
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
      use b2mod_elements
!      use b2mod_sign_of_magnetic_field
      use b2mod_b2cmgs
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmpt
      use b2mod_b2cmrc
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WDAT is an output routine that writes out B2 arrays into ASCII
*     data files.
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer ncall, iout, is, ix, iy, ib, ic, iw, iout_iout_ASTRA !srv 05.10.01 03.04.13 19.04.13
      integer jsep, jxi, jxa, ndat, k, mintopcut, ixx,             !srv 12.07.10
     *        nxtl, nxtr, ix1, ix2                                 !iyv 04.02.11
      integer use_astra
      integer, save :: ishift=0                                    !srv 05.11.09
      integer, save :: iappend = 0
      real (kind=R8) ::                                            !srv
     *  Total_current(-1:ny), Total_flux(-1:ny),                   !srv 31.08.00 19.02.02
     *  Total_abscurrent(-1:ny), Total_absflux(-1:ny),
     *  Total_energy_flux_i(-1:ny), Total_energy_flux_e(-1:ny),    !srv 23.11.00
     *  wrk2(-1:nx,-1:ny), wrkt(-1:ny,-1:nx),                      !srv 23.11.00
     *  Total_energy_flux(-1:ny), power,                           !srv 23.11.00
     *  Total_energy_flux_mdf(-1:ny),
     *  Total_energy_flux_i_mdf(-1:ny),
     *  Total_energy_flux_e_mdf(-1:ny),
!     *  kintmp,
     *  sfuadia(-1:ny,0:ns-1), sfvadia(-1:ny,0:ns-1),
     *  sfvaecrb(-1:ny,0:ns-1), sfna(-1:ny,0:ns-1),
     &  sfna_mdf(-1:ny,0:ns-1), sfn(-1:ny), sfn_mdf(-1:ny),
     &  boundary_flux(1:nbc,0:ns-1), sum_sdiv, sum_ssna,
     &  sdiv(0:ns-1), ssna(0:ns-1),                                !srv 12.04.12
     &  phi_app,
     &  pz(-1:nx,-1:ny), rz(-1:nx,-1:ny), cst_a(-1:nx,-1:ny,0:ns-1),
     &  pr(-1:nx,-1:ny), ro(-1:nx,-1:ny), cst(-1:nx,-1:ny)         !srv 30.01.15
      real (kind=R8) :: bsp(-1:ny), csp(-1:ny), dsp(-1:ny), dy,
     &     crd(4*(ny+1)+1), dseval, wdim(8*ns)                     !srv 03.04.13
      real (kind=R8) ::
     *     hlp, hlp_s, s
      real (kind=R8)  ::                                           !iyv 15.10.09
     *   qe_out_sum, qi_out_sum, q_out_sum,                        !iyv 04.02.11
     *   qe_inn_sum, qi_inn_sum, q_inn_sum,
     *   qe_oup_sum, qi_oup_sum, q_oup_sum,
     *   qe_iup_sum, qi_iup_sum, q_iup_sum,
     *   uvbtot, voltot                                            !iyv 29.11.13
      character chns*3                                             !srv 27.04.01
      character*32 hlp_frm                                         !srv 03.04.13
      character*26 boundary_location(nbcd)
      character(7), save :: folder                                 !srv 26.09.17
#ifdef B25_EIRENE
      integer js, ind, ias
#endif
*   ..procedures
      integer len_of_digits
      external len_of_digits
      external ipgeti
*   ..procedures                                                   !srv
      intrinsic abs                                                !srv
*   ..initialisation
      save ncall, iout, use_astra, jxi, jxa, jsep, nxtl, nxtr
      data ncall/0/, iout/0/, use_astra/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wdat')
*   ..initialise files and common

*   ..set internal parameters on first (and only) call
      if (ncall.eq.0) then
       call ipgeti('b2wdat_iout', iout)                           !srv 05.10.01
       call ipgeti('b2wdat_ishift', ishift)                       !srv 05.10.01
       call ipgeti('b2wdat_append', iappend)
       call ipgeti('b2mndr_astra', use_astra)
       call get_boundary_namelist
       call get_jsep(nx,ny,jxi,jxa,jsep)
       call get_nxt(nx, nxtl, nxtr) ! position of additional targets in DN configuration
*    ..try to create a file in optional output folder
       open(70,file='output/__1__.dat',err=10)
*    ..folder exists: future output will be done in folder
       folder='output/'
       write(70,*) 'Test file to check on output folder'
       close(70)
       go to 20
*    ..folder does not exist: future output will be done in current folder
   10  folder=''
   20  continue
      endif
* ..initialisation
      if (use_astra.ne.0) then
        iout_iout_ASTRA = 1
      else
        iout_iout_ASTRA = 0
      endif

      mintopcut = ny                                              !srv 25.02.16
      if(nncut > 0) then
        if (nnreg(0).eq.7) then
          mintopcut = topcut(1)
        else
          mintopcut = minval(topcut(1:nncut))                     !iyv 29.11.13 25.02.16
        endif
      endif
      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
      call b2xprz (nx, ny, ns, mp, am, na, rz)                    !srv 30.01.15 {
      cst = sqrt(pz/rz)
      do is=0,ns-1
        cst_a(:,:,is) = sqrt((rza(:,:,is)*te+ti)/(mp*am(is)))     !srv 24.03.10 {
      end do

*   ..write unchanging geometrical quantities only once
      if (ncall.eq.0) then

       if (iout.eq.1) then

        call my_out(70,-1,ny,conn,'conn')                         !srv

        wrk0=vol/hx**2                                             !srv
        call my_out(70,nx,ny,wrk0,'voldhx2')                       !srv
        wrk0=vol/(hy*qz(:,:,1))**2                                 !srv
        call my_out(70,nx,ny,wrk0,'voldhy2')                       !srv
        call my_out(70,nx,ny,OnedBsq,'OnedBsq')                    !srv
        call my_out(70,nx+1,ny+1,OnedBsqc,'OnedBsqc')              !srv

        wrk1=bb(:,:,0)/bb(:,:,3)
        call my_out(70,nx,ny,wrk1,'bx')                            !srv 18.12.14 }

        if (isymm.eq.1 .or. isymm.eq.2) then
         do iy=-1,ny
           do ix=-1,nx
             wrk0(ix,iy)=0.5_R8*pi*
     &        (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
           enddo
         enddo
         call my_out(70,nx,ny,wrk0,'2piR')                         !srv 15.09.14
         wrk1=wrk0-hz
         call my_out(70,nx,ny,wrk1,'2piRmhz')                      !srv 15.09.14
        else if (isymm.eq.3 .or. isymm.eq.4) then
         do iy=-1,ny
           do ix=-1,nx
             wrk0(ix,iy)=0.5_R8*pi*
     &        (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))
           enddo
         enddo
         call my_out(70,nx,ny,wrk0,'2piZ')                         !srv 15.09.14
         wrk1=wrk0-hz
         call my_out(70,nx,ny,wrk1,'2piZmhz')                      !srv 15.09.14
        end if

!       call my_out(70,nx,ny,pbs(-1,-1,1),'pbsy')

        call my_outi(70,nx,ny,region(-1,-1,0),'region0')
        call my_outi(70,nx,ny,region(-1,-1,1),'region1')
        call my_outi(70,nx,ny,region(-1,-1,2),'region2')

       end if

       if(iout.eq.2) then

        call b2wdat_open(70,'separatrix_location_on_plate')        !srv 26.09.17
        ix = -1
        if (LSN) then
         write(70,"(1x,a,12es14.6)")
     &    ' Separatrix location on Inner Lower plate =',
     &      cry(ix+ishift,mintopcut,0)
        else
         write(70,"(1x,a,12es14.6)")
     &    ' Separatrix location on Outer Upper plate =',
     &      cry(ix+ishift,mintopcut,0)
        endif
        ix = nx
        if (LSN) then
         write(70,"(1x,a,12es14.6)")
     &    ' Separatrix location on Outer Lower plate =',
     &      crx(ix-ishift,mintopcut,0)
        else
         write(70,"(1x,a,12es14.6)")
     &    ' Separatrix location on Inner Upper plate =',
     &      crx(ix-ishift,mintopcut,0)
        endif
        do ix = -1,nx
         if(leftix(ix,mintopcut).eq.-2 .and. ix.gt.-1) then
           write(70,"(1x,a,12es14.6)")
     &      ' Separatrix location on Outer Upper plate =',
     &        crx(ix+ishift,mintopcut,0)
         endif
         if(rightix(ix,mintopcut).eq.nx+1 .and. ix.lt.nx) then
           write(70,"(1x,a,12es14.6)")
     &      ' Separatrix location on Inner Upper plate =',
     &        cry(ix-ishift,mintopcut,0)
         endif
        enddo
        ix = nmdpl
        write(70,"(1x,a,12es14.6)")
     &   ' Separatrix location on Outer Midplane =',
     &     crx(ix-ishift,mintopcut,0)
        close(70)

       end if

       if(iout.eq.4 .or. iout.eq.3 .or. iout.eq.1) then           !iys 17.04.17
#ifdef B25_EIRENE
        if(iout.ne.4) then
         call b2wdat_open(70,'eirdiag')                           !srv 26.09.17
         write(70,"(a/6i5)") 'nlimps natm nmol npls nstra nstsid',
     &                        nlimps, natm, nmol, npls, nstra, nstsid
         write(70,"(a/6i5)") 'natmid nmolid ', natmid, nmolid
         write(70,"(a/6i5)") 'nnatmi nnmoli ', nnatmi, nnmoli
         close(70)
        end if
#endif
        call my_outi(70,nx,ny,bottomix,'bottomix')
        call my_outi(70,nx,ny,bottomiy,'bottomiy')
        call my_outi(70,nx,ny,topix,'topix')
        call my_outi(70,nx,ny,topiy,'topiy')
        call my_outi(70,nx,ny,leftix,'leftix')
        call my_outi(70,nx,ny,leftiy,'leftiy')
        call my_outi(70,nx,ny,rightix,'rightix')
        call my_outi(70,nx,ny,rightiy,'rightiy')

        call my_out(70,nx,ny,lengthx,'lengthx')                   !srv 13.02.17
        wrkt=transpose(lengthx)
        call my_out(70,ny,nx,wrkt,'lengthxTrans')                 !srv 13.02.17
        call my_out(70,nx,ny,lengthy,'lengthy')                   !srv
        wrk0=0.0d0                                                !srv 18.03.09 26.04.09 {
        do ix=-1,nx
         wrk0(ix,-1)=hy1(ix,-1)
         do iy=0,ny
           wrk0(ix,iy)=wrk0(ix,iy-1)+hy1(ix,iy)
         enddo
        enddo
        call my_out(70,nx,ny,wrk0,'lengthy1')                     !srv

        call my_out(70,nx,ny,len_surf,'len_surf')

        call my_out(70,nx,ny,crdy_surf,'crdy_surf')
        call my_out(70,nx,ny,crdy_c_cell,'crdy_c_cell')           !srv 18.03.09 26.04.09}

        call b2wdat_open(70,'info')                               !srv 26.09.17
        write(70,"(' Number of outer midplane cell =',i4)") nmdpl !srv 24.03.06 }
        write (70,*) ' nncut    =', nncut, ' nnreg =',nnreg       !srv
        write (70,*) ' leftcut  =',(leftcut(ix),ix=1,nncut)       !srv
        write (70,*) ' rightcut =',(rightcut(ix),ix=1,nncut)      !srv
        write (70,*) ' topcut   =',(topcut(ix),ix=1,nncut)        !srv
     &   , ' is the radial index of the ring directly above the X-point'
        write (70,*) ' bottomcut=',(bottomcut(ix),ix=1,nncut)     !srv
        write (70,*) ' Number of dead (or isolated) regions NNISO',NNISO
        if (NNISO.gt.0) then
         write (70,*) ' Column number of bottom left corner ',
     &    'of the dead region NXISO1', (NXISO1(IS),IS=1,NNISO)
         write (70,*) ' Column number of top right corner ',
     &    'of the dead region NXISO2', (NXISO2(IS),IS=1,NNISO)
         write (70,*) ' Row number of bottom left corner ',
     &    'of the dead region NYISO1', (NYISO1(IS),IS=1,NNISO)
         write (70,*) ' Row number of top right corner ',
     &    'of the dead region NYISO2', (NYISO2(IS),IS=1,NNISO)
        end if
        close (70)                                                !srv

        call my_out(70,nx,ny,crx(-1,-1,0),'crx0')
        call my_out(70,nx,ny,crx(-1,-1,1),'crx1')
        call my_out(70,nx,ny,crx(-1,-1,2),'crx2')
        call my_out(70,nx,ny,crx(-1,-1,3),'crx3')
        call my_out(70,nx,ny,cry(-1,-1,0),'cry0')
        call my_out(70,nx,ny,cry(-1,-1,1),'cry1')
        call my_out(70,nx,ny,cry(-1,-1,2),'cry2')
        call my_out(70,nx,ny,cry(-1,-1,3),'cry3')
        call my_out(70,nx,ny,bb(-1,-1,0),'bbx')
        call my_out(70,nx,ny,bb(-1,-1,2),'bbz')
        call my_out(70,nx,ny,bb(-1,-1,3),'bb')
        call my_out(70,nx,ny,hx,'hx')
        call my_out(70,nx,ny,hy,'hy')
        call my_out(70,nx,ny,hz,'hz')
        call my_out(70,nx,ny,qc,'qc')
        call my_out(70,nx,ny,hy1,'hy1')                          !srv 07.11.12
        call my_out(70,nx,ny,vol,'vol')
        call my_out(70,nx,ny,pbs(-1,-1,0),'pbsx')

        call my_out(70,nx,ny,wbbl(-1,-1,0),'bbxl')
        call my_out(70,nx,ny,wbbl(-1,-1,2),'bbzl')
        call my_out(70,nx,ny,wbbl(-1,-1,3),'bbl')
        call my_out(70,nx,ny,wbbr(-1,-1,0),'bbxr')
        call my_out(70,nx,ny,wbbr(-1,-1,2),'bbzr')
        call my_out(70,nx,ny,wbbr(-1,-1,3),'bbr')
        call my_out(70,nx,ny,wbbv(-1,-1,0),'bbxv')
        call my_out(70,nx,ny,wbbv(-1,-1,2),'bbzv')
        call my_out(70,nx,ny,wbbv(-1,-1,3),'bbv')
        wrk0(-1:nx,-1:ny)=wbbc(-1:nx,-1:ny,0)
        call my_out(70,nx,ny,wrk0,'bbxc')
        wrk0(-1:nx,-1:ny)=wbbc(-1:nx,-1:ny,2)
        call my_out(70,nx,ny,wrk0,'bbzc')
        wrk0(-1:nx,-1:ny)=wbbc(-1:nx,-1:ny,3)
        call my_out(70,nx,ny,wrk0,'bbc')
        call my_out(70,nx,ny,hxb,'hxb')
        call my_out(70,nx,ny,hyl,'hyl')
        call my_out(70,nx,ny,hzl,'hzl')
        call my_out(70,nx,ny,hzb,'hzb')
        call my_out(70,nx,ny,gs(-1,-1,0),'gsx')
        call my_out(70,nx,ny,gs(-1,-1,1),'gsy')
        call my_out(70,nx,ny,gs(-1,-1,2),'gsz')
        call my_out(70,nx,ny,qz(-1,-1,0),'qz0sint')
        call my_out(70,nx,ny,qz(-1,-1,1),'qz1cost')

c     calculate coordinates of cell centre
        wrk0=0.25_R8*(crx(:,:,0)+crx(:,:,1)+crx(:,:,2)+crx(:,:,3))  !srv 18.12.14 {
        call my_out(70,nx,ny,wrk0,'crc_')
        call b2wdat_open(70,'crc')   ! R-coordinate of the (ix,iy) cell centre !srv 26.09.17
        write(hlp_frm,"(a1,i4,a9)") '(',nx+2,'es24.16))'            !srv 07.08.15
        do iy = -1, ny
         write(70,hlp_frm) (wrk0(ix,iy),ix = -1, nx)                !srv 07.08.15
        enddo
        close(70)
        wrk1=0.25_R8*(cry(:,:,0)+cry(:,:,1)+cry(:,:,2)+cry(:,:,3))
        call my_out(70,nx,ny,wrk1,'czc_')
        call b2wdat_open(70,'czc')   ! Z-coordinate of the (ix,iy) cell centre !srv 26.09.17
        write(hlp_frm,"(a1,i4,a9)") '(',nx+2,'es24.16))'            !srv 07.08.15
        do iy = -1, ny
         write(70,hlp_frm) (wrk1(ix,iy),ix = -1, nx)                !srv 07.08.15
        enddo
        close(70)

       end if

      end if

      if (iout.eq.1) then                                           !srv begin if for iout

       Total_current=0.0_R8
       Total_abscurrent=0.0_R8
       Total_flux=0.0_R8
       Total_absflux=0.0_R8
       wrk2=0.0_R8
       call b2xpro (nx, ny, ns, mp, am, na, ro)                     !srv 30.01.15 {
       call b2xppr (nx, ny, ns, ne, na, te, ti, pr)
       call my_out(70,nx,ny,cst,'usound_cst')
       call my_out(70,nx,ny,ro,'ro')
       call my_out(70,nx,ny,rz,'rz')
       call my_out(70,nx,ny,pr,'pr')                                !srv 30.01.15 }
       call my_out(70,nx,ny,pz,'pz')                                !srv 30.01.15 }
       call my_out(70,nx,ny,ni(-1,-1,0),'ni0')                      !srv
       call my_out(70,nx,ny,ni(-1,-1,1),'ni1')                      !srv

       do is = 0,ns-1
        write (chns,'(i3.3)') is
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, wrk0)
        call my_out(70,nx,ny,wrk0,'pb'//chns)
        call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)               !srv
        call my_out(70,nx,ny,cst_a,'usound'//chns)
        wrk1 = ua(:,:,is)/wrk0
        call my_out(70,nx,ny,wrk1,'ua_M_'//chns)                    !srv 24.03.10 }
        call my_out(70,nx,ny,uadia(-1,-1,0,is),'uadiax'//chns)      !srv
        call my_out(70,nx,ny,uadia(-1,-1,1,is),'uadiay'//chns)      !srv
        call my_out(70,nx,ny,vaecrb(-1,-1,0,is),'vaecrbx'//chns)    !srv
        call my_out(70,nx,ny,vaecrb(-1,-1,1,is),'vaecrby'//chns)    !srv
        call my_out(70,nx,ny,vadia(-1,-1,0,is),'vadiax'//chns)      !srv
        call my_out(70,nx,ny,vadia(-1,-1,1,is),'vadiay'//chns)      !srv
        call my_out(70,nx,ny,wadia(-1,-1,0,is),'wadiax'//chns)      !srv
        call my_out(70,nx,ny,wadia(-1,-1,1,is),'wadiay'//chns)      !srv
        call my_out(70,nx,ny,fna(-1,-1,0,is),'fnax'//chns)          !srv
        call my_out(70,nx,ny,fna(-1,-1,1,is),'fnay'//chns)          !srv
        call my_out(70,nx,ny,fna_mdf(-1,-1,0,is),'fnax_mdf'//chns)  !srv 13.10.06
        call my_out(70,nx,ny,fna_mdf(-1,-1,1,is),'fnay_mdf'//chns)  !srv 13.10.06
        call my_out(70,nx,ny,fna_eir(-1,-1,0,is),'fna_eirx'//chns)  !srv 14.07.10
        call my_out(70,nx,ny,fna_eir(-1,-1,1,is),'fna_eiry'//chns)  !srv 13.07.10
        call my_out(70,nx,ny,fnaPSch(-1,-1,0,is),'fnaPSchx'//chns)  !srv 13.10.06
        call my_out(70,nx,ny,fnaPSch(-1,-1,1,is),'fnaPSchy'//chns)  !srv 13.10.06
        call my_out(70,nx,ny,fna_he(-1,-1,0,is),'fna_hex'//chns)
        call my_out(70,nx,ny,fna_he(-1,-1,1,is),'fna_hey'//chns)
        call my_out(70,nx,ny,fna_fcor(-1,-1,0,is),'fna_fcorx'//chns)
        call my_out(70,nx,ny,fna_fcor(-1,-1,1,is),'fna_fcory'//chns)
        call my_out(70,nx,ny,fna_nodrift(-1,-1,0,is),
     *                      'fna_nodriftx'//chns)
        call my_out(70,nx,ny,fna_nodrift(-1,-1,1,is),
     *                      'fna_nodrifty'//chns)
        call my_out(70,nx,ny,sna(-1,-1,0,is),'sna0'//chns)          !srv 25.02.16
        call my_out(70,nx,ny,sna(-1,-1,1,is),'sna1'//chns)          !srv 25.02.16

        call my_out(70,nx,ny,cdna(-1,-1,0,is),'cdnax'//chns)        !srv
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'cdnay'//chns)        !srv
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'cdpax'//chns)        !srv
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'cdpay'//chns)        !srv
        call my_out(70,nx,ny,cddi(-1,-1,0,is),'cddix'//chns)        !srv 13.10.06
        call my_out(70,nx,ny,cddi(-1,-1,1,is),'cddiy'//chns)        !srv 13.10.06
        call my_out(70,nx,ny,cvsa_an(-1,-1,0,is),'cvsax_an'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'cvsax_cl'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'cvsay'//chns)
        call my_out(70,nx,ny,kinrgy(-1,-1,is),'kinrgy'//chns)       !srv 08.10.04
        call my_out(70,nx,ny,vsa0(-1,-1,is),'vsaAN'//chns)          !srv 25.02.16

        wrk0=0.0d0
        do iy = -1,ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) wrk0(ix,iy)=
     1        wbbl(ix,iy,2)/wbbl(ix,iy,3)
     4       *vaecrb(ix,iy,0,is)*
     5       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/         !srv 01.03.00
     6       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))            !srv 01.03.00
     7       *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
          enddo
        enddo
        call my_out(70,nx,ny,wrk0,'fvaecrbx'//chns)                 !srv

        wrk1=0.0d0
        do iy = -1,ny
          do ix = -1, nx
            if (leftix(ix,iy).ne.-2) wrk1(ix,iy)=
     1        wbbl(ix,iy,2)/wbbl(ix,iy,3)
     4       *vadia(ix,iy,0,is)*
     5       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/         !srv 01.03.00
     6       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))            !srv 01.03.00
     7       *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
          enddo
        enddo
        call my_out(70,nx,ny,wrk1,'fvadiax'//chns)                  !srv

        wrk0=fna(:,:,0,is)-wrk0-wrk1
        call my_out(70,nx,ny,wrk0,'fdifx'//chns)                    !srv

        wrk0=0.0d0
        do iy = -1, ny
          do ix = -1, nx
            if (bottomiy(ix,iy).ne.-2)
     >        wrk0(ix,iy)=vaecrb(ix,iy,1,is)*(na(ix,iy,is)
     &          +na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8*
     5       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/     !srv 01.03.00
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))      !srv 01.03.00
          enddo
        enddo
        call my_out(70,nx,ny,wrk0,'fvaecrby'//chns)                 !srv
        wrk1=0.0d0
        do iy = -1, ny
          do ix = -1, nx
            if (bottomiy(ix,iy).ne.-2)
     >        wrk1(ix,iy)=vadia(ix,iy,1,is)*(na(ix,iy,is)
     &          +na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8*
     5       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/     !srv 01.03.00
     6       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))      !srv 01.03.00
          enddo
        enddo
        call my_out(70,nx,ny,wrk1,'fvadiay'//chns)                  !srv
        wrk0=fna(:,:,1,is)-wrk0-wrk1
        call my_out(70,nx,ny,wrk0,'fdify'//chns)                    !srv

        call my_out(70,nx,ny,fmo(-1,-1,0,is),'fmox'//chns)          !srv 24.02.10
        call my_out(70,nx,ny,fmo(-1,-1,1,is),'fmoy'//chns)

        call b2wdat_open(70,'Total_flux_pm'//chns)                  !srv 19.02.02 { 26.09.17
        write (70,*) '  N      Total_flux    Total_|flux| '
        do iy = mintopcut,0,-1
          Total_flux(iy) = 0.0d0
          Total_absflux(iy) = 0.0d0
          do ix=-1,nx
            if(bottomiy(ix,iy).eq.-2) cycle
            if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then ! for CORE
              Total_flux(iy) = Total_flux(iy) + fmo(ix,iy,1,is)
              Total_absflux(iy)=Total_absflux(iy)+abs(fmo(ix,iy,1,is))
            endif
          enddo
          write (70,'(1x,i4,2(2x,es14.7))') iy, Total_flux(iy),
     &           Total_absflux(iy)
        enddo
        close (70)                                                !srv 19.02.02 }

        call my_out(70,nx,ny,rqahe(-1,-1,is),'rqahe'//chns)       !srv
        call my_out(70,nx,ny,rqrad(-1,-1,is),'rqrad'//chns)       !xpb
        call my_out(70,nx,ny,rqbrm(-1,-1,is),'rqbrm'//chns)       !xpb
        call my_out(70,nx,ny,pcca(-1,-1,0,is),'pccax'//chns)      !srv 11.01.13
        call my_out(70,nx,ny,pcca(-1,-1,1,is),'pccay'//chns)      !srv 11.01.13

        call my_out(70,nx,ny,fchanml_a(-1,-1,0,is),               !srv 08.09.21 {
     &                      'fchanml_ax'//chns)
        call my_out(70,nx,ny,fchanml_a(-1,-1,1,is),
     &                      'fchanml_ay'//chns)
        call my_out(70,nx,ny,fchinert_a(-1,-1,0,is),
     &                      'fchinert_ax'//chns)
        call my_out(70,nx,ny,fchinert_a(-1,-1,1,is),
     &                      'fchinert_ay'//chns)
        call my_out(70,nx,ny,fchvispar_a(-1,-1,0,is),
     &                      'fchvispar_ax'//chns)
        call my_out(70,nx,ny,fchvispar_a(-1,-1,1,is),
     &                      'fchvispar_ay'//chns)
        call my_out(70,nx,ny,fchvisper_a(-1,-1,0,is),
     &                      'fchvisper_ax'//chns)
        call my_out(70,nx,ny,fchvisper_a(-1,-1,1,is),
     &                      'fchvisper_ay'//chns)
        call my_out(70,nx,ny,fchvisq_a(-1,-1,0,is),
     &                      'fchvisq_ax'//chns)
        call my_out(70,nx,ny,fchvisq_a(-1,-1,1,is),
     &                      'fchvisq_ay'//chns)                   !srv 08.09.21 }
       enddo

       call my_out(70,nx,ny,pccm(-1,-1,0),'pccmx')                !srv 11.01.13
       call my_out(70,nx,ny,pccm(-1,-1,1),'pccmy')                !srv 11.01.13
       call my_out(70,nx,ny,fhe(-1,-1,0),'fhex')                  !srv 23.11.00
       call my_out(70,nx,ny,fhe(-1,-1,1),'fhey')                  !srv 23.11.00
       call my_out(70,nx,ny,fhi(-1,-1,0),'fhix')                  !srv 23.11.00
       call my_out(70,nx,ny,fhi(-1,-1,1),'fhiy')                  !srv 23.11.00

       call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'fhe_mdfx')          !srv 13.10.06 {
       call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'fhe_mdfy')
       call my_out(70,nx,ny,fhePSch(-1,-1,0),'fhePSchx')
       call my_out(70,nx,ny,fhePSch(-1,-1,1),'fhePSchy')
       call my_out(70,nx,ny,fhi_mdf(-1,-1,0),'fhi_mdfx')
       call my_out(70,nx,ny,fhi_mdf(-1,-1,1),'fhi_mdfy')
       call my_out(70,nx,ny,fhiPSch(-1,-1,0),'fhiPSchx')
       call my_out(70,nx,ny,fhiPSch(-1,-1,1),'fhiPSchy')          !srv 13.10.06 }

       call my_out(70,nx,ny,te,'te')                              !srv
       call my_out(70,nx,ny,ti,'ti')                              !srv
       call my_out(70,nx,ny,fch(-1,-1,0),'fchx')                  !srv
       call my_out(70,nx,ny,fch(-1,-1,1),'fchy')                  !srv
       call my_out(70,nx,ny,fch_p(-1,-1,0),'fchx_p')              !srv 10.06.08
       call my_out(70,nx,ny,fch_p(-1,-1,1),'fchy_p')              !srv 10.06.08
       call my_out(70,nx,ny,fchvispar(-1,-1,0),'fchvisparx')
       call my_out(70,nx,ny,fchvispar(-1,-1,1),'fchvispary')
       call my_out(70,nx,ny,fchvisper(-1,-1,0),'fchvisperx')
       call my_out(70,nx,ny,fchvisper(-1,-1,1),'fchvispery')
       call my_out(70,nx,ny,fchvisq(-1,-1,0),'fchvisqx')
       call my_out(70,nx,ny,fchvisq(-1,-1,1),'fchvisqy')
       call my_out(70,nx,ny,fchinert(-1,-1,0),'fchinertx')        !srv 07.07.99
       call my_out(70,nx,ny,fchinert(-1,-1,1),'fchinerty')        !srv 07.07.99
       call my_out(70,nx,ny,fchdia(-1,-1,0),'fchdiax')
       call my_out(70,nx,ny,fchdia(-1,-1,1),'fchdiay')
       call my_out(70,nx,ny,fchin(-1,-1,0),'fchinx')
       call my_out(70,nx,ny,fchin(-1,-1,1),'fchiny')
       call my_out(70,nx,ny,fchanml(-1,-1,0),'fchanmlx')
       call my_out(70,nx,ny,fchanml(-1,-1,1),'fchanmly')
       call b2tlnl (nx, ny, 0, te, ti, ne, wrk0)                  !srv 02.07.13 {
       call my_out(70,nx,ny,wrk0,'lnlam_Br')
       call b2tlnl (nx, ny, 1, te, ti, ne, wrk0)
       call my_out(70,nx,ny,wrk0,'lnlam_We_ii')
       call b2tlnl (nx, ny, 2, te, ti, ne, wrk0)
       call my_out(70,nx,ny,wrk0,'lnlam_We_ee')
       call b2tlnl (nx, ny, 3, te, ti, ne, wrk0)                  !srv 02.07.13 }
       call my_out(70,nx,ny,wrk0,'lnlam_We_ei')
       call my_out(70,nx,ny,csig(-1,-1,0),'csig_x')               !srv
       call my_out(70,nx,ny,csig(-1,-1,1),'csig_y')               !srv
       call my_out(70,nx,ny,csig_cl(-1,-1,0),'csig_clx')          !srv
       call my_out(70,nx,ny,csig_cl(-1,-1,1),'csig_cly')          !srv
       call my_out(70,nx,ny,csig_an(-1,-1,0),'csig_anx')          !srv
       call my_out(70,nx,ny,csig_an(-1,-1,1),'csig_any')          !srv
       if (ismain0.ne.ismain) then
        call my_out(70,nx,ny,csigin(-1,-1,0,ismain,ismain0),
     &                      'csig_inx')                           !srv
        call my_out(70,nx,ny,csigin(-1,-1,1,ismain,ismain0),
     &                      'csig_iny')                           !srv
       endif
       call my_out(70,nx,ny,cvsa(-1,-1,0,ismain),'cvsa')          !srv
       call my_out(70,nx,ny,cvsa_cl(-1,-1,0,ismain),'cvsa_cl')    !srv
       call my_out(70,nx,ny,calf(-1,-1,0),'calfx')                !srv
       call my_out(70,nx,ny,calf(-1,-1,1),'calfy')                !srv
       call my_out(70,nx,ny,calf_cl(-1,-1,0),'calf_clx')          !srv
       call my_out(70,nx,ny,calf_cl(-1,-1,1),'calf_cly')          !srv
       call my_out(70,nx,ny,calf_an(-1,-1,0),'calf_anx')          !srv
       call my_out(70,nx,ny,calf_an(-1,-1,1),'calf_any')          !srv
       call my_out(70,nx,ny,chce(-1,-1,0),'chcex')                !srv
       call my_out(70,nx,ny,chce(-1,-1,1),'chcey')                !srv
       call my_out(70,nx,ny,chci(-1,-1,0),'chcix')                !srv
       call my_out(70,nx,ny,chci(-1,-1,1),'chciy')                !srv

       do is = 0,ns-1                                             !srv 16.11.14 {
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,resco(-1,-1,is),'resco'//chns)
        call my_out(70,nx,ny,resmo(-1,-1,is),'resmo'//chns)
       enddo                                                      !srv 16.11.14 }
       call my_out(70,nx,ny,respo,'respo')                        !srv
       call my_out(70,nx,ny,reshe,'reshe')                        !srv
       call my_out(70,nx,ny,reshi,'reshi')                        !srv
       call my_out(70,nx,ny,she(-1,-1,0),'she0')                  !srv
       call my_out(70,nx,ny,she(-1,-1,1),'she1')                  !srv
       call my_out(70,nx,ny,she(-1,-1,2),'she2')                  !srv
       call my_out(70,nx,ny,she(-1,-1,3),'she3')                  !srv
       call my_out(70,nx,ny,shi(-1,-1,0),'shi0')                  !srv
       call my_out(70,nx,ny,shi(-1,-1,1),'shi1')                  !srv
       call my_out(70,nx,ny,shi(-1,-1,2),'shi2')                  !srv
       call my_out(70,nx,ny,shi(-1,-1,3),'shi3')                  !srv
       call my_out(70,nx,ny,sch(-1,-1,0),'sch0')                  !srv
       call my_out(70,nx,ny,sch(-1,-1,1),'sch1')                  !srv
       call my_out(70,nx,ny,sch(-1,-1,2),'sch2')                  !srv
       call my_out(70,nx,ny,sch(-1,-1,3),'sch3')                  !srv
       wrk0=sna(:,:,0,ismain)+sna(:,:,1,ismain)*na(:,:,ismain)
       call my_out(70,nx,ny,wrk0,'sour_ni')                       !srv

       call b2wdat_open(70,'Total_current')                       !srv 31.08.00 23.07.10 { 26.09.17
       write (70,*) '  N      Total_current  Total_|current| '
       do iy = mintopcut,0,-1
        Total_current(iy) = 0.0_R8
        Total_abscurrent(iy) = 0.0_R8
        do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) cycle
          if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then ! for CORE
            Total_current(iy) = Total_current(iy) + fch(ix,iy,1)
            Total_abscurrent(iy)=Total_current(iy)+abs(fch(ix,iy,1))
          endif
        enddo
        write (70,'(1x,i4,2(2x,es14.7))') iy, Total_current(iy),
     &                                        Total_abscurrent(iy)
       enddo
       ix = jxa ; iy = 1
       write (70,'(1x,a,i4,a,i4,a,g14.7)')
     &         'potential in (',ix,',',iy,')=',po(ix,iy)
       close (70)

*   ..write physics parameters after changing                     !srv 19.12.07 {
       call b2wdat_open(70,'b2mn_chk')                            !srv 26.09.17
       call b2wfcp (70, 0, ns)
       close(70)                                                  !srv 19.12.07 }
*   ..write "experimental" data
       ndat = 1
       crd(ndat) = crdy_c_cell(nmdpl,-1)
       do iy = -1, ny-1
         dy = (crdy_c_cell(nmdpl,iy+1)-crdy_c_cell(nmdpl,iy))/4.0_R8
         do k=1,4
           crd(ndat+1)=crd(ndat)+dy
           ndat = ndat+1
         enddo
       enddo

       call dspline(ny+2,crdy_c_cell(nmdpl,:),ne(nmdpl,:),bsp,csp,dsp)
       call b2wdat_open(70,'splne')                               !srv 26.09.17
       do iy=1,ndat
         write(70,"(1x,2(es14.6))") crd(iy),
     &    dseval(ny+2,crd(iy),crdy_c_cell(nmdpl,:),ne(nmdpl,:),
     &           bsp,csp,dsp)
       enddo
       close(70)

       call dspline(ny+2,crdy_c_cell(nmdpl,:),te(nmdpl,:),bsp,csp,dsp)
       call b2wdat_open(70,'splte')                               !srv 26.09.17
       do iy=1,ndat
         write(70,"(1x,2(es14.6))") crd(iy),
     &    dseval(ny+2,crd(iy),crdy_c_cell(nmdpl,:),te(nmdpl,:),
     &           bsp,csp,dsp)
       enddo
       close(70)

       call dspline(ny+2,crdy_c_cell(nmdpl,:),ti(nmdpl,:),bsp,csp,dsp)
       call b2wdat_open(70,'splti')                               !srv 26.09.17
       do iy=1,ndat
         write(70,"(1x,2(es14.6))") crd(iy),
     &    dseval(ny+2,crd(iy),crdy_c_cell(nmdpl,:),ti(nmdpl,:),
     &           bsp,csp,dsp)
       enddo
       close(70)

       call b2wdat_open(70,'expne')                               !srv 26.09.17
       do iy=-1,ny
        write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), ne(nmdpl,iy)
       enddo
       close(70)
       call b2wdat_open(70,'expte')                               !srv 26.09.17
       do iy=-1,ny
        write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), te(nmdpl,iy)
       enddo
       close(70)
       call b2wdat_open(70,'expti')                               !srv 26.09.17
       do iy=-1,ny
        write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), ti(nmdpl,iy)
       enddo
       close(70)

*   ..calculate total core energy flux
       power = 0.0_R8
       if (boundary_namelist.eq.1) then
         do ib = 1, nbc
           if(core_boundary(ib)) then
             do ix = 1, bc_list_size(ib)
               power = power +
     .         (fhe(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     .          fhi(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .              topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))
             enddo
           endif
         enddo
       else
         do ix = -1, nx
           if (on_closed_surface(ix,0)) then
             power = power + (fhe(ix,0,1)+fhi(ix,0,1))
           endif
         enddo
       endif
       call b2wdat_open(70,'exppw')                                !srv 26.09.17
       write(70,"(es14.6)") power
       close(70)

*   ..calculate integrated particle fluxes for checking
       call b2wdat_open(70,'Total_particle_fluxes_in_Core')        !srv 26.09.17
       do is = 0,ns-1
         write (70,"(a,1x,a2,i1)") 'Total_particle_flux ',
     &     adjustr(is_codes(is)),nint(zamin(is))
         write(70,'(6a)')
     &     ' iy ','   sfuadia  ','   sfvadia  ','   sfvaecrb ',
     &     '    sfna   ','    sfna_mdf '
         do iy = ny,-1, -1
           sfuadia(iy,is)  =0.0_R8
           sfvadia(iy,is)  =0.0_R8
           sfvaecrb(iy,is) =0.0_R8
           sfna(iy,is)     =0.0_R8
           sfna_mdf(iy,is) =0.0_R8
           do ix = -1, nx
             if(on_closed_surface(ix,iy).and.
     &          bottomiy(ix,iy).ne.-2) then  ! for CORE            !srv 01.02.12
               sfuadia(iy,is) = sfuadia(iy,is) + uadia(ix,iy,1,is)*
     5          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &          (na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &           na(ix,iy,is))/2.0_R8
               sfvadia(iy,is) = sfvadia(iy,is) + vadia(ix,iy,1,is)*
     5          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &          (na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &           na(ix,iy,is))/2.0_R8
               sfvaecrb(iy,is) = sfvaecrb(iy,is) + vaecrb(ix,iy,1,is)*
     5          (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     6          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &          (na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &           na(ix,iy,is))/2.0_R8
               sfna(iy,is) = sfna(iy,is) + fna(ix,iy,1,is)
               sfna_mdf(iy,is) = sfna_mdf(iy,is) + fna_mdf(ix,iy,1,is)
             endif
           enddo
           write(70,"(i4,5es15.7)") iy, sfuadia(iy,is), sfvadia(iy,is),
     &      sfvaecrb(iy,is), sfna(iy,is), sfna_mdf(iy,is)
         enddo
       enddo

       write (70,"(a)") 'Total particle flux over all particles '  !srv 12.04.12 {
       write(70,"(' iy ','   sfn    ','   sfn_mdf ')")
       do iy = ny,-1, -1
         sfn    (iy) = 0.0_R8
         sfn_mdf(iy) = 0.0_R8
         do is = 0,ns-1
           sfn(iy)     =  sfn(iy)     + sfna(iy,is)
           sfn_mdf(iy) =  sfn_mdf(iy) + sfna_mdf(iy,is)
         enddo
         write(70,"(i4,5es15.7)") iy, sfn(iy), sfn_mdf(iy)
       enddo
       close(70)

      endif              !srv   endif for iout=1

      if(iout.eq.1 .or. iout_iout_ASTRA.eq.1
     &                .or. iout.eq.3 .or. iout.eq.4) then                       !srv 19.04.13
!______________________________________________________________________________ !iyv 15.10.09 {
! target energy fluxes
        wrk0=0.0e0_R8
        wrk1=0.0e0_R8

        if(nncut.eq.2) then            !iyv 04.02.11
          ix1=nxtl    !inner upper
          ix2=nxtr    !outer upper
        endif

        do iy=-1,ny
          do is = 0, ns-1
            if(.not.is_neutral(is)) then          !iyv 07.02.11
! pr_West
              wrk1(-1,iy)=wrk1(-1,iy)+na(-1,iy,is)*
     &         (rza(-1,iy,is)*te(-1,iy)+ti(-1,iy))
! pr_East
              wrk1(0,iy)=wrk1(0,iy)+na(nx,iy,is)*
     &         (rza(nx,iy,is)*te(nx,iy)+ti(nx,iy))
! ro_West
              wrk1(1,iy) = wrk1(1,iy)+am(is)*mp*na(-1,iy,is)
! ro_East
              wrk1(2,iy) = wrk1(2,iy)+am(is)*mp*na(nx,iy,is)

              if(nncut.eq.2) then !{    !iyv 04.02.11
! pr_inner_upper
                wrk1(3,iy)=wrk1(3,iy)+na(ix1,iy,is)*
     &           (rza(ix1,iy,is)*te(ix1,iy)+ti(ix1,iy))
! pr_outer_upper
                wrk1(4,iy)=wrk1(4,iy)+na(ix2,iy,is)*
     &           (rza(ix2,iy,is)*te(ix2,iy)+ti(ix2,iy))
! ro_inner_upper
                wrk1(5,iy) = wrk1(5,iy)+am(is)*mp*na(ix1,iy,is)
! ro_outer_upper
                wrk1(6,iy) = wrk1(6,iy)+am(is)*mp*na(ix2,iy,is)
              endif
            endif
          enddo

! qe West
          phi_app = get_phi_applied(boundary_namelist,-1,iy)
          wrk0(-1,iy)=hy1(-1,iy)/hy(-1,iy)*abs(bb(-1,iy,0)/bb(-1,iy,3))*    !srv 19.04. 13
     &     ( te(-1,iy)*(2.0e0_R8)+qe*po(-1,iy)*(1.0_R8-gammae) )*
     &     ( ne(-1,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(-1,iy)/me)*
     &     expu2(max(-50.0_R8, -qe*(po(-1,iy)-phi_app)/te(-1,iy))) )
! qe East
          phi_app = get_phi_applied(boundary_namelist,nx,iy)
          wrk0(0,iy)=hy1(nx,iy)/hy(nx,iy)*abs(bb(nx,iy,0)/bb(nx,iy,3))*    !srv 19.04. 13
     &     ( te(nx,iy)*(2.0e0_R8)+qe*po(nx,iy)*(1.0_R8-gammae) )*
     &     ( ne(nx,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(nx,iy)/me)*
     &     expu2(max(-50.0_R8, -qe*(po(nx,iy)-phi_app)/te(nx,iy))) )

          do is = 0, ns-1
            if(.not.is_neutral(is)) then
              wrk0(1,iy)=wrk0(1,iy)+na(-1,iy,is)
              wrk0(2,iy)=wrk0(2,iy)+na(nx,iy,is)
           endif
         enddo
! qi West
         wrk0(3,iy)=2.5e0_R8*
     &               hy1(-1,iy)/hy(-1,iy)*abs(bb(-1,iy,0)/bb(-1,iy,3))*  ! IYS 22.03.2017 changed bb(:,:,0) to abs(bb(:,:,0))
     &               ti(-1,iy)*wrk0(1,iy)*sqrt(wrk1(-1,iy)/wrk1(1,iy))
!  IYS added the heat flux carried by recycling neutrals
     &        - b2stbr_shi(-1,iy)/hy(-1,iy)/hz(-1,iy)                ! IYS 21.11.2012
!  IYS believes that even if the coupling to EIRENE is active, the EIRENE
!  does not disturb the guard cells and all the fluid sources related to neutrals
!  are multiplied by a very small factor (like 1.0e-10) in the coupling mode

! qi East
         wrk0(4,iy)=2.5e0_R8*
     &               hy1(nx,iy)/hy(nx,iy)*abs(bb(nx,iy,0)/bb(nx,iy,3))*  ! IYS 22.03.2017 changed bb(:,:,0) to abs(bb(:,:,0))
     &               ti(nx,iy)*wrk0(2,iy)*sqrt(wrk1(0,iy)/wrk1(2,iy))
!  IYS added the heat flux carried by recycling neutrals
     &        - b2stbr_shi(nx,iy)/hy(nx,iy)/hz(nx,iy)  ! IYS 21.11.2012
!  IYS believes that even if the coupling to EIRENE is active, the EIRENE
!  does not disturb the guard cells and all the fluid sources related to neutrals
!  are multiplied by a very small factor (like 1.0e-10) in the coupling mode

         if(nncut.eq.2) then !{      !iyv 04.02.11

! qe inner upper
           phi_app = get_phi_applied(boundary_namelist,ix1,iy)
           wrk0(7,iy)=hy1(ix1,iy)/hy(ix1,iy)*                             !srv 19.04. 13
     &        abs(bb(ix1,iy,0)/bb(ix1,iy,3))*
     &      ( te(ix1,iy)*(2.0e0_R8)+qe*po(ix1,iy)*(1.0_R8-gammae) )*
     &      ( ne(ix1,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(ix1,iy)/me)*
     &       expu2(max(-50.0_R8, -qe*(po(ix1,iy)-phi_app)/te(ix1,iy))) )
! qe outer upper
           phi_app = get_phi_applied(boundary_namelist,ix2,iy)
           wrk0(8,iy)=hy1(ix2,iy)/hy(ix2,iy)*                             !srv 19.04. 13
     &        abs(bb(ix2,iy,0)/bb(ix2,iy,3))*
     &      ( te(ix2,iy)*(2.0e0_R8)+qe*po(ix2,iy)*(1.0_R8-gammae) )*
     &      ( ne(ix2,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(ix2,iy)/me)*
     &       expu2(max(-50.0_R8, -qe*(po(ix2,iy)-phi_app)/te(ix2,iy))) )

           do is = 0, ns-1
             if(.not.is_neutral(is)) then
               wrk0(9,iy) =wrk0(9,iy) +na(ix1,iy,is)
               wrk0(10,iy)=wrk0(10,iy)+na(ix2,iy,is)
             endif
           enddo
! qi inner upper
           wrk0(11,iy)=2.5e0_R8*
     &           hy1(ix1,iy)/hy(ix1,iy)*abs(bb(ix1,iy,0)/bb(ix1,iy,3))*  ! IYS 22.03.2017 changed bb(:,:,0) to abs(bb(:,:,0))
     &               ti(ix1,iy)*wrk0(9,iy)*sqrt(wrk1(3,iy)/wrk1(5,iy))
!  IYS added the heat flux carried by recycling neutrals
     &        - b2stbr_shi(ix1,iy)/hy(ix1,iy)/hz(ix1,iy)  ! IYS 21.11.2012
!  IYS believes that even if the coupling to EIRENE is active, the EIRENE
!  does not disturb the guard cells and all the fluid sources related to neutrals
!  are multiplied by a very small factor (like 1.0e-10) in the coupling mode

! qi outer upper
           wrk0(12,iy)=2.5e0_R8*
     &           hy1(ix2,iy)/hy(ix2,iy)*abs(bb(ix2,iy,0)/bb(ix2,iy,3))*  ! IYS 22.03.2017 changed bb(:,:,0) to abs(bb(:,:,0))
     &               ti(ix2,iy)*wrk0(10,iy)*sqrt(wrk1(4,iy)/wrk1(6,iy))
!  IYS added the heat flux carried by recycling neutrals
     &        - b2stbr_shi(ix2,iy)/hy(ix2,iy)/hz(ix2,iy)  ! IYS 21.11.2012
!  IYS believes that even if the coupling to EIRENE is active, the EIRENE
!  does not disturb the guard cells and all the fluid sources related to neutrals
!  are multiplied by a very small factor (like 1.0e-10) in the coupling mode

         endif
       enddo

! next part is copied from b2usrtrc
!________________________________________________________________

      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        s=gs(0,iy,0)
        hlp=fht(0,iy,0)
        hlp=abs(hlp)/s
        wrk0(5,iy) = hlp
        hlp_s=hlp_s+s
      end do !}
#ifdef B25_EIRENE
      do iy=0,ny-1 !{
        do js=1,nsts !{
          if (eirdiag_nds_typ(js).ne.2) cycle
          if (eirdiag_nds_srf(js).ne.0) cycle
          if (eirdiag_nds_start(js).gt.(iy+1)) cycle
          if (eirdiag_nds_end(js).lt.(iy+1)) cycle
          ind=eirdiag_nds_ind(js)
          ias=ind+1-(eirdiag_nds_start(js)-1)
          wrk0(5,iy)=wrk0(5,iy)+ewldt_res(iy+ias)
        enddo !}
      enddo !}
#endif
      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        s=gs(nx,iy,0)
        hlp=fht(nx,iy,0)
        hlp=abs(hlp)/s
        wrk0(6,iy) = hlp
        hlp_s=hlp_s+s
      end do !}
#ifdef B25_EIRENE
      do iy=0,ny-1 !{
        do js=1,nsts !{
          if (eirdiag_nds_typ(js).ne.2) cycle
          if (eirdiag_nds_srf(js).ne.nx) cycle
          if (eirdiag_nds_start(js).gt.(iy+1)) cycle
          if (eirdiag_nds_end(js).lt.(iy+1)) cycle
          ind=eirdiag_nds_ind(js)
          ias=ind+1-(eirdiag_nds_start(js)-1)
          wrk0(6,iy)=wrk0(6,iy)+ewldt_res(iy+ias)
        enddo !}
      enddo !}
#endif

! for double-null geometry     !iyv 04.02.11
! outer upper target

      if(nncut.eq.2) then !{

        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          s=gs(ix2+1,iy,0)
          hlp=fht(ix2+1,iy,0)
          hlp=abs(hlp)/s
          wrk0(14,iy) = hlp
          hlp_s=hlp_s+s
        end do !}
#ifdef B25_EIRENE
        do iy=0,ny-1 !{
          do js=1,nsts !{
            if (eirdiag_nds_typ(js).ne.2) cycle
            if (eirdiag_nds_srf(js).ne.ix2+1) cycle
            if (eirdiag_nds_start(js).gt.(iy+1)) cycle
            if (eirdiag_nds_end(js).lt.(iy+1)) cycle
            ind=eirdiag_nds_ind(js)
            ias=ind+1-(eirdiag_nds_start(js)-1)
            wrk0(14,iy)=wrk0(14,iy)+ewldt_res(iy+ias)
          enddo !}
        enddo !}
#endif
! inner upper target
        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          s=gs(ix1,iy,0)
          hlp=fht(ix1,iy,0)
          hlp=abs(hlp)/s
          wrk0(13,iy) = hlp
          hlp_s=hlp_s+s
        end do !}
#ifdef B25_EIRENE
        do iy=0,ny-1 !{
          do js=1,ntns !{
            if (eirdiag_nds_typ(js).ne.2) cycle
            if (eirdiag_nds_srf(js).ne.ix1) cycle
            if (eirdiag_nds_start(js).gt.(iy+1)) cycle
            if (eirdiag_nds_end(js).lt.(iy+1)) cycle
            ind=eirdiag_nds_ind(js)
            ias=ind+1-(eirdiag_nds_start(js)-1)
            wrk0(13,iy)=wrk0(13,iy)+ewldt_res(iy+ias)
          enddo !}
        end do !}
#endif
      end if !}

!________________________________________________________________

      qe_out_sum = 0.0d0
      qe_inn_sum = 0.0d0
      qi_out_sum = 0.0d0
      qi_inn_sum = 0.0d0
      q_inn_sum  = 0.0d0
      q_out_sum  = 0.0d0

      do iy = -1, ny
        qe_out_sum = qe_out_sum + wrk0(0,iy) *hy(nx,iy)*hz(nx,iy)
        qe_inn_sum = qe_inn_sum + wrk0(-1,iy)*hy(-1,iy)*hz(-1,iy)
        qi_out_sum = qi_out_sum + wrk0(4,iy) *hy(nx,iy)*hz(nx,iy)
        qi_inn_sum = qi_inn_sum + wrk0(3,iy) *hy(-1,iy)*hz(-1,iy)
        q_out_sum  = q_out_sum  + wrk0(6,iy) *hy(nx,iy)*hz(nx,iy)
        q_inn_sum  = q_inn_sum  + wrk0(5,iy) *hy(-1,iy)*hz(-1,iy)
      enddo

      call b2wdat_open(70,'q_plates.dat')                            !srv 26.09.17

      if(nncut.eq.2) then !{
        qe_oup_sum = 0.0d0
        qe_iup_sum = 0.0d0
        qi_oup_sum = 0.0d0
        qi_iup_sum = 0.0d0
        q_oup_sum  = 0.0d0
        q_iup_sum  = 0.0d0

        do iy = -1, ny
          qe_oup_sum = qe_oup_sum + wrk0(8,iy)  *hy(ix2,iy)*hz(ix2,iy)
          qe_iup_sum = qe_iup_sum + wrk0(7,iy)  *hy(ix1,iy)*hz(ix1,iy)
          qi_oup_sum = qi_oup_sum + wrk0(12,iy) *hy(ix2,iy)*hz(ix2,iy)
          qi_iup_sum = qi_iup_sum + wrk0(11,iy) *hy(ix1,iy)*hz(ix1,iy)
          q_oup_sum  = q_oup_sum  + wrk0(14,iy) *hy(ix2,iy)*hz(ix2,iy)
          q_iup_sum  = q_iup_sum  + wrk0(13,iy) *hy(ix1,iy)*hz(ix1,iy)
        enddo

        write (70,"(5x, 'nxtr, ix2:', 2i4, 5x, 'nxtl, ix1:', 2i4 )")
     &                   nxtr, ix2, nxtl, ix1

        write (70,"(a7, 16a11)") 'iy',
     &              'qsum_olo','qsum_oup', 'qsum_ilo','qsum_iup',
     &              'qe_olo', 'qe_oup', 'qi_olo',  'qi_oup',
     &              'qei_olo','qei_oup',
     &              'qe_ilo', 'qe_iup', 'qi_ilo',  'qi_iup',
     &              'qei_ilo','qei_iup'
        do iy = ny, -1, -1
          write(70,"(5x,i2,16(1x,es10.3))") iy,
     &             wrk0( 6,iy), wrk0(14,iy), wrk0(5,iy), wrk0(13,iy),
     &             wrk0( 0,iy), wrk0( 8,iy), wrk0(4,iy), wrk0(12,iy),
     &             wrk0( 0,iy)+wrk0(4,iy),   wrk0(8,iy)+wrk0(12,iy),
     &             wrk0(-1,iy), wrk0( 7,iy), wrk0(3,iy), wrk0(11,iy),
     &             wrk0(-1,iy)+wrk0(3,iy),   wrk0(7,iy)+wrk0(11,iy)
        enddo
        write(70,41) qe_oup_sum, qi_oup_sum, qe_out_sum, qi_out_sum,
     &               qe_iup_sum, qi_iup_sum, qe_inn_sum, qi_inn_sum,
     &               q_oup_sum, q_out_sum, q_iup_sum, q_inn_sum

   41   format(/'Energy e-flux to outer upper plate     =',es11.3/
     &          'Energy i-flux to outer upper plate     =',es11.3/
     &          'Energy e-flux to outer lower plate     =',es11.3/
     &          'Energy i-flux to outer lower plate     =',es11.3//

     &          'Energy e-flux to inner lower plate     =',es11.3/
     &          'Energy i-flux to inner lower plate     =',es11.3/
     &          'Energy e-flux to inner upper plate     =',es11.3/
     &          'Energy i-flux to inner upper plate     =',es11.3//

     &          'Total Energy flux to outer upper plate =',es11.3/
     &          'Total Energy flux to outer lower plate =',es11.3/
     &          'Total Energy flux to inner upper plate =',es11.3/
     &          'Total Energy flux to inner lower plate =',es11.3/ )

      else

        write (70,
     &   "(5x,a2,5x,a6,5x,a6,4x,a7,4x,a7,5x,a6,5x,a6,4x,a7,4x,a7)")
     &   'iy','qe_out','qi_out','qei_out','q_outtg','qe_inn',
     &   'qi_inn','qei_inn','q_inntg'
        do iy = ny, -1, -1
          write(70,"(5x,i2,8(1x,es10.3))") iy,
     &             wrk0( 0,iy), wrk0( 4,iy), wrk0(0,iy)+wrk0(4,iy),
     &             wrk0( 6,iy), wrk0(-1,iy), wrk0(3,iy),
     &             wrk0(-1,iy)+wrk0(3,iy),   wrk0(5,iy)
        enddo
        write(70,51) qe_out_sum, qi_out_sum, q_out_sum,
     &               qe_inn_sum, qi_inn_sum, q_inn_sum

   51   format(/'Energy e-flux to Eastern plate      =',es11.3/
     &          'Energy i-flux to Eastern plate      =',es11.3/
     &          'Total energy flux to Eastern plate  =',es11.3//
     &          'Energy e-flux to Western plate      =',es11.3/
     &          'Energy i-flux to Western plate      =',es11.3/
     &          'Total energy flux to Western plate  =',es11.3/)
      endif
      close(70)

      ishift = 0

      if (LSN) then
        call b2wdat_open(70,'table_lower_inner')                   !srv 26.09.17
      else
        call b2wdat_open(70,'table_lower_outer')                   !srv 26.09.17
      endif
      write(70,"(1x,9a14)")
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
      do iy = 0,ny
        ixx=-1 !max(-1,ix-1)
        write(70,"(1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx+ishift,iy),
     &     0.5_R8*( crx(ixx+ishift,iy,0) + crx(ixx+ishift,iy,2) ),
     &     te(ixx+ishift,iy)/ev, ti(ixx+ishift,iy)/ev,
     &     ne(ixx+ishift,iy),
     &     qe * hy1(ixx+ishift,iy)/hy(ixx+ishift,iy) *
     &     abs(bb(ixx+ishift,iy,0))/bb(ixx+ishift,iy,3) *
     &     na(ixx+ishift,iy,ismain) *
     &     sqrt((te(ixx+ishift,iy)+ti(ixx+ishift,iy))/(am(ismain)*mp)),
     &     wrk0(-1,iy)/1.0e6_R8,wrk0(3,iy)/1.0e6_R8,wrk0(5,iy)/1.0e6_R8
      enddo
      close(70)

      if(LSN) then
        call b2wdat_open(70,'table_lower_outer')                   !srv 26.09.17
      else
        call b2wdat_open(70,'table_lower_inner')                   !srv 26.09.17
      endif
      write(70,"(1x,9a14)")
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
      do iy = 0,ny
        ixx=nx !min(nx,ix+1)
        write(70,"(1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx-ishift,iy),
     &     0.5_R8*( crx(ixx-ishift,iy,0) + crx(ixx-ishift,iy,2) ),
     &     te(ixx-ishift,iy)/ev, ti(ixx-ishift,iy)/ev,
     &     ne(ixx-ishift,iy),
     &     qe * hy1(ixx-ishift,iy)/hy(ixx-ishift,iy) *
     &     abs(bb(ixx-ishift,iy,0))/bb(ixx-ishift,iy,3) *
     &     na(ixx-ishift,iy,ismain) *
     &     sqrt((te(ixx-ishift,iy)+ti(ixx-ishift,iy))/(am(ismain)*mp)),
     &     wrk0(0,iy)/1.0e6_R8, wrk0(4,iy)/1.0e6_R8, wrk0(6,iy)/1.0e6_R8
      enddo
      close(70)                                       !srv 19.09.11 }

      if(nncut.eq.2) then

        call b2wdat_open(70,'table_upper_inner')      !srv 26.09.17
        write(70,"(1x,9a14)")
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
        do iy = 0,ny
          ixx=ix1
          write(70,"( 1x,9es14.6)")                              !srv 03.04.13
     &     crdy_c_cell(ixx-ishift,iy),
     &     0.5_R8*( crx(ixx-ishift,iy,0) + crx(ixx-ishift,iy,2) ),
     &     te(ixx-ishift,iy)/ev, ti(ixx-ishift,iy)/ev,
     &     ne(ixx-ishift,iy),
     &     qe * hy1(ixx-ishift,iy)/hy(ixx-ishift,iy) *
     &     abs(bb(ixx-ishift,iy,0))/bb(ixx-ishift,iy,3) *
     &     na(ixx-ishift,iy,ismain) *
     &     sqrt((te(ixx-ishift,iy)+ti(ixx-ishift,iy))/(am(ismain)*mp)),
     &     wrk0(7,iy)/1.0e6_R8, wrk0(11,iy)/1.0e6_R8,
     &     wrk0(13,iy)/1.0e6_R8
        enddo
        close(70)

        call b2wdat_open(70,'table_upper_outer')                   !srv 26.09.17
        write(70,"(1x,9a14)")
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
        do iy = 0,ny
          ixx=ix2
          write(70,"( 1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx+ishift,iy),
     &     0.5_R8*( crx(ixx+ishift,iy,0) + crx(ixx+ishift,iy,2) ),
     &     te(ixx+ishift,iy)/ev, ti(ixx+ishift,iy)/ev,
     &     ne(ixx+ishift,iy),
     &     qe * hy1(ixx+ishift,iy)/hy(ixx+ishift,iy) *
     &     abs(bb(ixx+ishift,iy,0))/bb(ixx+ishift,iy,3) *
     &     na(ixx+ishift,iy,ismain) *
     &     sqrt((te(ixx+ishift,iy)+ti(ixx+ishift,iy))/(am(ismain)*mp)),
     &     wrk0(8,iy)/1.0e6_R8, wrk0(12,iy)/1.0e6_R8,
     &     wrk0(14,iy)/1.0e6_R8
        enddo
        close(70)                                         !srv 19.09.11 }
       endif
      endif      !  iout.eq.1 or. iout_iout_ASTRA.eq.1    !srv 19.04.13
!______________________________________________________________________________ !iyv 15.10.09 }

      if(iout.eq.2) then                                  !srv 06.05.09 {
*   ..write crdy_c_cell, crdy_c_cell_abs [mm], ne, te, ti, js at Western edge
        if(LSN) then
          call b2wdat_open(70,'inner')                    !srv 26.09.17
        else
          call b2wdat_open(70,'outer')                    !srv 26.09.17
        endif
        do iy = -1,ny                                     !srv 03.04.13 {
          iw = 0
          do ix = -1,nx
            if(leftix(ix,iy).eq.-2) then
              ixx=max(-1,ix-1)
              wdim(iw+1)= crdy_c_cell(ixx+ishift,iy)
              wdim(iw+2)= ne(ixx+ishift,iy)
              wdim(iw+3)= te(ixx+ishift,iy)/ev
              wdim(iw+4)= ti(ixx+ishift,iy)/ev
              iw = iw + 4
            endif
          enddo
          write(70,"(2(1x,4es14.6))") (wdim(ix),ix=1,iw)
        enddo                                             !srv 03.04.13 }
        close(70)

        call b2wdat_open(70,'midpln')                     !srv 26.09.17
        ix=nmdpl
        do iy=-1,ny
          write(70,"(1x,5es14.6)") crdy_c_cell(ix,iy),
     &     ne(ix,iy), te(ix,iy)/ev, ti(ix,iy)/ev
        enddo
        close(70)

        if (LSN) then
          call b2wdat_open(70,'outer')                    !srv 26.09.17
        else
          call b2wdat_open(70,'inner')                    !srv 26.09.17
        endif
        do iy = -1,ny                                     !srv 03.04.13 {
         iw = 0
         do ix = -1,nx
           if(rightix(ix,iy).eq.nx+1) then
             ixx=min(nx,ix+1)
             wdim(iw+1)= crdy_c_cell(ixx-ishift,iy)
             wdim(iw+2)= ne(ixx-ishift,iy)
             wdim(iw+3)= te(ixx-ishift,iy)/ev
             wdim(iw+4)= ti(ixx-ishift,iy)/ev
             iw = iw + 4
           endif
         enddo
         write(70,"(2(1x,4es14.6))") (wdim(ix),ix=1,iw)
       enddo                                              !srv 03.04.13 }
       close(70)

       call b2wdat_open(70,'midplnUi')                    !srv 26.09.17
       ix=nmdpl
       do iy=-1,ny
         write(70,"(1x,2es14.6)") crdy_c_cell(ix,iy), ua(ix,iy,1)
       enddo
       close(70)

       call b2wdat_open(70,'midplnEy')                    !srv 26.09.17
       ix=nmdpl
       do iy=0,ny
         write(70,"(1x,2es14.6)")
     *    crdy_surf(ix,iy),
     *     -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*2.0e0_R8
       enddo
       close(70)

       call b2wdat_open(70,'OMP_parallel_electron_flux')
       do iy=-1,ny
         write(70,"(1x,2es14.6)") crdy_c_cell(jxa,iy),
     *     fne(jxa,iy,0)/pitx(jxa,iy)
       enddo
       close(70)

       call b2wdat_open(70,'IMP_parallel_electron_flux')
       do iy=-1,ny
         write(70,"(1x,2es14.6)") crdy_c_cell(jxi,iy),
     *     fne(jxi,iy,0)/pitx(jxi,iy)
       enddo
       close(70)

*   ..calculate total core energy flux
       power = 0.0_R8
       if (boundary_namelist.eq.1) then
         do ib = 1, nbc
           if(core_boundary(ib)) then
             do ix = 1, bc_list_size(ib)
               power = power +
     .          (fhe(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1)+
     .           fhi(topix(bc_list_x(ix,ib),bc_list_y(ix,ib)),
     .               topiy(bc_list_x(ix,ib),bc_list_y(ix,ib)),1))
             enddo
           endif
         enddo
       else
         do ix = -1, nx
           if (on_closed_surface(ix,0)) then
             power = power + (fhe(ix,0,1)+fhi(ix,0,1))
           endif
         enddo
       endif
       call b2wdat_open(70,'power')                                !srv 26.09.17
       write(70,"('Power =',es14.6)") power
       close(70)

*   ..write crdy_c_cell, na, ua, cst at Western edges
       if (LSN) then
         call b2wdat_open(70,'inner_na_ua_cst')                   !srv 03.04.13 { 26.09.17
       else
         call b2wdat_open(70,'outer_na_ua_cst')                   !srv 03.04.13 { 26.09.17
       endif
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+2,'es14.6))'
       do iy = -1,ny
         iw = 0
         do ix = -1,nx
          if(leftix(ix,iy).eq.-2) then
            ixx=max(-1,ix-1)
            iw = iw +1
            wdim(iw)= crdy_c_cell(ixx+ishift,iy)
            do is=0,ns-1
              iw = iw + 1
              wdim(iw)= na(ixx+ishift,iy,is)
            enddo
            do is=0,ns-1
              iw = iw + 1
              wdim(iw)= ua(ixx+ishift,iy,is)
            enddo
            iw = iw + 1
            wdim(iw) = cst(ixx+ishift,iy)
          endif
        enddo
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, na, ua, cst at Eastern edges
       if (LSN) then
        call b2wdat_open(70,'outer_na_ua_cst')                  !srv 26.09.17
       else
        call b2wdat_open(70,'inner_na_ua_cst')                  !srv 26.09.17
       endif
       do iy = -1,ny
        iw = 0
        do ix = -1,nx
          if(rightix(ix,iy).eq.nx+1) then
            ixx=min(nx,ix+1)
            iw = iw +1
            wdim(iw)= crdy_c_cell(ixx-ishift,iy)
            do is=0,ns-1
              iw = iw + 1
              wdim(iw)= na(ixx-ishift,iy,is)
            enddo
            do is=0,ns-1
              iw = iw + 1
              wdim(iw)= ua(ixx-ishift,iy,is)
            enddo
            iw = iw + 1
            wdim(iw) = cst(ixx-ishift,iy)
          endif
        enddo
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)                                               !srv 03.04.13 }

*   ..write crdy_c_cell, na, ua, cst at Inner Midplane
       call b2wdat_open(70,'IMP_na_ua_cst')
       do iy = -1,ny
        iw = 1
        wdim(iw)= crdy_c_cell(jxi,iy)
        do is=0,ns-1
          iw = iw + 1
          wdim(iw)= na(jxi,iy,is)
        enddo
        do is=0,ns-1
          iw = iw + 1
          wdim(iw)= ua(jxi,iy,is)
        enddo
        iw = iw + 1
        wdim(iw) = cst(jxi,iy)
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, parallel power flux densities at Inner Midplane
       call b2wdat_open(70,'IMP_parallel_power')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+5,'es14.6))'
       do iy = -1,ny
        iw = 1
        wdim(iw)= crdy_c_cell(jxi,iy)
        iw = 2
        wdim(iw) = fhe(jxi,iy,0)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        iw = 3
        wdim(iw) = fhi(jxi,iy,0)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        iw = 4
        wdim(iw) = fhj(jxi,iy,0)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        do is=0,ns-1
         iw = iw + 1
         wdim(iw)= fhm(jxi,iy,0,is)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        enddo
        do is=0,ns-1
         iw = iw + 1
         wdim(iw)= fhp(jxi,iy,0,is)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        enddo
        iw = iw + 1
        wdim(iw) = fht(jxi,iy,0)/gs(jxi,iy,0)/qc(jxi,iy)
     .           /pitx(jxi,iy)
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, na, ua, cst at Outer Midplane
       call b2wdat_open(70,'OMP_na_ua_cst')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+2,'es14.6))'
       do iy = -1,ny
        iw = 1
        wdim(iw)= crdy_c_cell(jxa,iy)
        do is=0,ns-1
          iw = iw + 1
          wdim(iw)= na(jxa,iy,is)
        enddo
        do is=0,ns-1
          iw = iw + 1
          wdim(iw)= ua(jxa,iy,is)
        enddo
        iw = iw + 1
        wdim(iw) = cst(jxa,iy)
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, parallel power flux densities at Outer Midplane
       call b2wdat_open(70,'OMP_parallel_power')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+5,'es14.6))'
       do iy = -1,ny
        iw = 1
        wdim(iw)= crdy_c_cell(jxa,iy)
        iw = 2
        wdim(iw) = fhe(jxa,iy,0)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        iw = 3
        wdim(iw) = fhi(jxa,iy,0)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        iw = 4
        wdim(iw) = fhj(jxa,iy,0)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        do is=0,ns-1
         iw = iw + 1
         wdim(iw)= fhm(jxa,iy,0,is)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        enddo
        do is=0,ns-1
         iw = iw + 1
         wdim(iw)= fhp(jxa,iy,0,is)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        enddo
        iw = iw + 1
        wdim(iw) = fht(jxa,iy,0)/gs(jxa,iy,0)/qc(jxa,iy)
     .           /pitx(jxa,iy)
        write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, na, ua, cst at Eastern divertor entrances
       call b2wdat_open(70,'Eastern_div_entry_na_ua_cst')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+2,'es14.6))'
       do ic = 1, nncut
         if (ic.eq.1) then
           ixx = rightcut(1)
         else if (ic.eq.2) then
           ixx = leftcut(2)
         end if
         do iy = topcut(ic), ny
           iw = 1
           wdim(iw)= crdy_c_cell(ixx,iy)
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= na(ixx,iy,is)
           enddo
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= ua(ixx,iy,is)
           enddo
           iw = iw + 1
           wdim(iw) = cst(ixx,iy)
           write(70,hlp_frm) (wdim(ix),ix=1,iw)
         enddo
       enddo
       close(70)

*   ..write crdy_c_cell, parallel power flux densities at Eastern divertor entrances
       call b2wdat_open(70,'Eastern_div_entry_parallel_power')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+5,'es14.6))'
       do ic = 1, nncut
         if (ic.eq.1) then
           ixx = rightcut(1)
         else if (ic.eq.2) then
           ixx = leftcut(2)
         end if
         do iy = topcut(ic), ny
           iw = 1
           wdim(iw)= crdy_c_cell(ixx,iy)
           iw = 2
           wdim(iw) = fhe(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           iw = 3
           wdim(iw) = fhi(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           iw = 4
           wdim(iw) = fhj(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           do is=0,ns-1
             iw = iw + 1
             wdim(iw)= fhm(ixx,iy,0,is)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           enddo
           do is=0,ns-1
             iw = iw + 1
             wdim(iw)= fhp(ixx,iy,0,is)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           enddo
           iw = iw + 1
           wdim(iw) = fht(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           write(70,hlp_frm) (wdim(ix),ix=1,iw)
         enddo
       enddo
       close(70)

*   ..write crdy_c_cell, na, ua, cst at Western divertor entrances
       call b2wdat_open(70,'Western_div_entry_na_ua_cst')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+2,'es14.6))'
       do ic = 1, nncut
         if (ic.eq.1) then
           ixx = leftcut(1)-1
         else if (ic.eq.2) then
           ixx = rightcut(2)-1
         end if
         do iy = topcut(ic), ny
           iw = 1
           wdim(iw)= crdy_c_cell(ixx,iy)
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= na(ixx,iy,is)
           enddo
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= ua(ixx,iy,is)
           enddo
           iw = iw + 1
           wdim(iw) = cst(ixx,iy)
           write(70,hlp_frm) (wdim(ix),ix=1,iw)
         enddo
       enddo
       close(70)

*   ..write crdy_c_cell, parallel power flux densities at Western divertor entrances
       call b2wdat_open(70,'Western_div_entry_parallel_power')
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+5,'es14.6))'
       do ic = 1, nncut
         if (ic.eq.1) then
           ixx = leftcut(1)-1
         else if (ic.eq.2) then
           ixx = rightcut(2)-1
         end if
         do iy = topcut(ic), ny
           iw = 1
           wdim(iw)= crdy_c_cell(ixx,iy)
           iw = 2
           wdim(iw) = fhe(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           iw = 3
           wdim(iw) = fhi(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           iw = 4
           wdim(iw) = fhj(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= fhm(ixx,iy,0,is)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           enddo
           do is = 0, ns-1
             iw = iw + 1
             wdim(iw)= fhp(ixx,iy,0,is)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           enddo
           iw = iw + 1
           wdim(iw) = fht(ixx,iy,0)/gs(ixx,iy,0)/qc(ixx,iy)/
     .                                             pitx(ixx,iy)
           write(70,hlp_frm) (wdim(ix),ix=1,iw)
         enddo
       enddo
       close(70)

*   ..write crdy_c_cell, j_sat at Western edges
       if (LSN) then
         call b2wdat_open(70,'inner_jsat')                   !srv 03.04.13 { 26.09.17
       else
         call b2wdat_open(70,'outer_jsat')                   !srv 03.04.13 { 26.09.17
       endif
       write(hlp_frm,'(a)') '(2(1x,2es14.6))'
       do iy = -1,ny
         wdim = 0.0_R8
         iw = 0
         do ix = -1,nx
           if(leftix(ix,iy).eq.-2) then
             ixx=max(-1,ix-1)
             iw = iw +1
             wdim(iw)= crdy_c_cell(ixx+ishift,iy)
             iw = iw +1
             do is = 0, ns-1
               wdim(iw)= wdim(iw)+qe*rza(ixx+ishift,iy,is)*
     .                   na(ixx+ishift,iy,is)*cst_a(ixx+ishift,iy,is)
             enddo
           endif
         enddo
         write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write crdy_c_cell, j_sat at Eastern edges
       if (LSN) then
         call b2wdat_open(70,'outer_jsat')                   !srv 03.04.13 { 26.09.17
       else
         call b2wdat_open(70,'inner_jsat')                   !srv 03.04.13 { 26.09.17
       endif
       do iy = -1,ny
         wdim = 0.0_R8
         iw = 0
         do ix = -1,nx
           if(rightix(ix,iy).eq.nx+1) then
             ixx=min(nx,ix+1)
             iw = iw +1
             wdim(iw)= crdy_c_cell(ixx-ishift,iy)
             iw = iw +1
             do is = 0, ns-1
               wdim(iw)= wdim(iw)+qe*rza(ixx-ishift,iy,is)*
     .                   na(ixx-ishift,iy,is)*cst_a(ixx-ishift,iy,is)
             enddo
           endif
         enddo
         write(70,hlp_frm) (wdim(ix),ix=1,iw)
       enddo
       close(70)

*   ..write poloidal coordinates, Up upper separatrix
       iy = mintopcut+4                                        !srv 13.02.17 {
       write (chns,'(i3.3)') iy
       call b2wdat_open(70,'poloidal_ua_at_'//chns)            !srv 26.09.17
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+1,'es14.6))'
       do ix = -1,nx
        iw = 1
        wdim(iw)= lengthx(ix,iy)
        do is = 0, ns-1
          iw = iw + 1
          wdim(iw)= ua(ix,iy,is)
        enddo
        write(70,hlp_frm) (wdim(is),is=1,iw)
       enddo
       close(70)

*   ..write poloidal coordinates, na upper separatrix
       iy = mintopcut+4
       write (chns,'(i3.3)') iy
       call b2wdat_open(70,'poloidal_na_at_'//chns)            !srv 26.09.17
       write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+1,'es14.6))'
       do ix = -1,nx
        iw = 1
        wdim(iw)= lengthx(ix,iy)
        do is=0,ns-1
          iw = iw + 1
          wdim(iw)= na(ix,iy,is)
        enddo
        write(70,hlp_frm) (wdim(is),is=1,iw)
       enddo
       close(70)                                               !srv 13.02.17 }

      endif                     !srv   endif for iout=2        !srv 06.05.09 }

      if(iout.eq.4 .or. iout.eq.3 .or. iout.eq.1) then         !iys 17.04.17

       Total_energy_flux=0.0_R8
       Total_energy_flux_mdf=0.0_R8                            !srv 30.04.15 {
       Total_energy_flux_e=0.0_R8
       Total_energy_flux_i=0.0_R8
       Total_energy_flux_e_mdf=0.0_R8                          !srv 30.04.15 }
       Total_energy_flux_i_mdf=0.0_R8

#ifdef B25_EIRENE
       if (iout.ne.4) then
        do is = 1,nnatmi                                            !srv 04.06.09 16.02.12{
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,is,1),'dab2'//chns)
         call my_out(70,nx,ny,tab2(0:nx+1,0:ny+1,is,1),'tab2'//chns)
         call my_out(70,nx,ny,rfluxa(0:nx+1,0:ny+1,is,1),'rfluxa'//chns)
         call my_out(70,nx,ny,refluxa(0:nx+1,0:ny+1,is,1),
     .                       'refluxa'//chns)
         call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,is,1),'pfluxa'//chns)
         call my_out(70,nx,ny,pefluxa(0:nx+1,0:ny+1,is,1),
     .                       'pefluxa'//chns)
        enddo
        do is = 1,nnmoli                                            !srv 16.02.12
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,is,1),'dmb2'//chns) !srv 29.04.10
         call my_out(70,nx,ny,tmb2(0:nx+1,0:ny+1,is,1),'tmb2'//chns) !srv 29.04.10
         call my_out(70,nx,ny,rfluxm(0:nx+1,0:ny+1,is,1),'rfluxm'//chns)
         call my_out(70,nx,ny,refluxm(0:nx+1,0:ny+1,is,1),
     .                       'refluxm'//chns)
         call my_out(70,nx,ny,pfluxm(0:nx+1,0:ny+1,is,1),'pfluxm'//chns)
         call my_out(70,nx,ny,pefluxm(0:nx+1,0:ny+1,is,1),
     .                       'pefluxm'//chns)
        enddo                                                       !srv 04.06.09 }
       endif
#endif

       call my_out(70,nx,ny,po,'po')                                !srv

       wrk0=0.0e0_R8
       do iy=0,ny
         do ix=-1,nx
           wrk0(ix,iy)=
     *         -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *          (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))    !srv 08.05.07
     *          *2.0e0_R8
         enddo
       enddo
       call my_out(70,nx,ny,wrk0,'Eyf')                              !srv 15.06.04 }
       wrk0=0.0e0_R8                                                 !iyv 15.10.09 {
       do iy=0,ny
         do ix=-1,nx
           if(topiy(ix,iy).ne.ny+1) then
             wrk0(ix,iy)=
     *        -(po(topix(ix,iy),topiy(ix,iy))-
     *          po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *         (hy1(ix,iy)+0.5e0_R8*
     *         (hy1(bottomix(ix,iy),bottomiy(ix,iy))) +
     *          hy1(topix(ix,iy),topiy(ix,iy)))
           else
             wrk0(ix,iy)=
     *        -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy))) !
     *        *2.0e0_R8
           endif
         enddo
       enddo
       call my_out(70,nx,ny,wrk0,'Eyc')

*     calculate shear
       wrk1=0.0e0_R8
       do iy=-1,ny
         do ix=-1,nx
           wrk1(ix,iy)= wrk0(ix,iy)/(hz(ix,iy)*bb(ix,iy,0))
         enddo
       enddo

       wrk0=0.0e0_R8
       do iy=-1,ny
         do ix=-1,nx
           if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     *       wrk0(ix,iy)=abs(hz(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*
     *              (wrk1(topix(ix,iy),topiy(ix,iy))-
     *               wrk1(bottomix(ix,iy),bottomiy(ix,iy)))/
     *                (hy1(ix,iy)+0.5e0_R8*
     *                (hy1(bottomix(ix,iy),bottomiy(ix,iy))) +
     *                 hy1(topix(ix,iy),topiy(ix,iy))))
         enddo
       enddo
       call my_out(70,nx,ny,wrk0,'Wsh')                              !iyv 15.10.09 }

       uvbtot=0.0e0_R8
       voltot=0.0e0_R8
       wrk0=0.0e0_R8
       wrk1=0.0e0_R8
       wrk2=0.0e0_R8
       do iy = 0,mintopcut-1
         do ix = 0,nx
           if(on_closed_surface(ix,iy)) then
             uvbtot=uvbtot+
     &              ua(ix,iy,ismain)*vol(ix,iy)*bb(ix,iy,3)
             voltot=voltot+vol(ix,iy)
             wrk0(ix,iy)=
     &           (ne(ix,iy)-ne(bottomix(ix,iy),bottomiy(ix,iy)))/
     &           (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
             wrk1(ix,iy)=
     &           (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))/
     &           (hy(ix,iy)+hy(bottomix(ix,iy),bottomiy(ix,iy)))
           endif
         enddo
       enddo
       do iy = 0,mintopcut-1
         do ix = 0,nx
           if(on_closed_surface(ix,iy)) then
             wrk2(ix,iy) = -uvbtot/voltot*bb(ix,iy,0)/bb(ix,iy,3)+
     &           (ti(ix,iy)/ne(ix,iy)*
     &           (wrk0(ix,iy)+wrk0(topix(ix,iy),topiy(ix,iy)))+
     &            kt_neo(ix,iy)*
     &           (wrk1(ix,iy)+wrk1(topix(ix,iy),topiy(ix,iy))))/qe
           endif
         enddo
       enddo
       call my_out(70,nx,ny,kt_neo,'kt_neo')
       call my_out(70,nx,ny,wrk2,'Ey_neo')                         !iyv 15.10.09 }

       call b2xzef(nx,ny,ns,rz2,na,ne,wrk0)                        !srv 13.01.17
       call my_out(70,nx,ny,wrk0,'Zeff')                           !iyv 12.12.13 }

       if (iout.ne.4) then
         wrk2=0.0_R8
         do is = 0,ns-1
           write (chns,'(i3.3)') is
           call my_out(70,nx,ny,na(-1,-1,is),'na'//chns)           !srv
           call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)           !srv
           wrk2 = ua(:,:,is)/cst                                   !srv 30.01.15
           call my_out(70,nx,ny,wrk2,'ua_Mcst_'//chns)             !srv 30.01.15
         end do
       end if

       call my_out(70,nx,ny,ne,'ne')                               !srv
       wrk0=te/ev
       call my_out(70,nx,ny,wrk0,'te_eV')
       wrk0=ti/ev
       call my_out(70,nx,ny,wrk0,'ti_eV')
       call b2wdat_open(70,'Total_energy_flux_i')                  !srv 26.09.17
       write (70,*) ' N  Total_energy_flux_i Total_energy_flux_i_mdf '

       do iy = mintopcut,0,-1
        Total_energy_flux_i(iy) = 0.0_R8
        Total_energy_flux_i_mdf(iy) = 0.0_R8
        do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) cycle
          if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then ! for CORE
            Total_energy_flux_i(iy) =
     &       Total_energy_flux_i(iy) + fhi(ix,iy,1)
            Total_energy_flux_i_mdf(iy) =
     &       Total_energy_flux_i_mdf(iy) + fhi_mdf(ix,iy,1)
          endif
        enddo
        write (70,'(1x,i4,2es18.7)')
     &          iy, Total_energy_flux_i(iy), Total_energy_flux_i_mdf(iy)
       enddo
       close (70)

       call b2wdat_open(70,'Total_energy_flux_e')               !srv 26.09.17
       write (70,*) '  N  Total_energy_flux_e Total_energy_flux_e_mdf'

       do iy = mintopcut,0,-1
        Total_energy_flux_e(iy) = 0.0_R8
        Total_energy_flux_e_mdf(iy) = 0.0_R8
        do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) cycle
          if(on_closed_surface(bottomix(ix,iy),bottomiy(ix,iy))) then ! for CORE
            Total_energy_flux_e(iy) =
     &       Total_energy_flux_e(iy) + fhe(ix,iy,1)
            Total_energy_flux_e_mdf(iy) =
     &       Total_energy_flux_e_mdf(iy) + fhe_mdf(ix,iy,1)
          endif
        enddo
        write (70,'(1x,i4,2es18.7)')
     &          iy, Total_energy_flux_e(iy), Total_energy_flux_e_mdf(iy)
       enddo
       close (70)

       call b2wdat_open(70,'Total_energy_flux')                    !srv 26.09.17
       write (70,*) '  N  Total_energy_flux  Total_energy_flux_mdf'

       do iy = mintopcut,0,-1
        Total_energy_flux(iy) = Total_energy_flux_e(iy)
     *                         +Total_energy_flux_i(iy)
        Total_energy_flux_mdf(iy) = Total_energy_flux_e_mdf(iy)
     *                             +Total_energy_flux_i_mdf(iy)
        write (70,'(1x,i4,2es18.7)')
     &          iy, Total_energy_flux(iy), Total_energy_flux_mdf(iy)
       enddo
       close (70)                                                 !srv 31.08.00 23.07.10 }

       call b2wdat_open(70,'Total_boundary_particle_fluxes')      !srv 26.09.17

       sum_sdiv = 0.0_R8
       sum_ssna= 0.0_R8

       do is=0, ns-1
        write (hlp_frm,'(a,i1,a)')
     .   '(a,1x,a2,i',len_of_digits(nint(zamin(is))),')'
        write (70,hlp_frm) 'Total_boundary_particle_flux ',
     &               adjustr(is_codes(is)),nint(zamin(is))

        sdiv(is) = 0.0_R8
        ssna(is) = 0.0_R8

        do ib=1, nbc
         boundary_flux(ib,is)=0.0e0_R8
         IF(BCCHAR(IB).EQ.'S') THEN
           DO IX = 1, BC_LIST_SIZE(IB)
             boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &            TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
           end do
           sdiv(is) = sdiv(is) - boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).EQ.'N') THEN
           DO IX = 1, BC_LIST_SIZE(IB)
            boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
           end do
           sdiv(is) = sdiv(is) + boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).EQ.'W') THEN
           DO IY = 1, BC_LIST_SIZE(IB)
            boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &            rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
           end do
           sdiv(is) = sdiv(is) - boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY = 1, BC_LIST_SIZE(IB)
             boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
           end do
           sdiv(is) = sdiv(is) + boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).NE.'X' .AND. BCCHAR(IB).NE.'Y') THEN
           CALL XERRAB
     1         ('B2WDAT: INVALID OPTION')
         ENDIF
         if (bcpos(ib).eq.-2) then
           write(boundary_location(ib),'(a23,i2)')
     .        ' on boundary cell list ',IB
         else
           write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .        ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
         endif
         WRITE(70,"(a,a,a,es15.7)") BCCHAR(IB),
     2         boundary_location(ib),
     &         ' flux =',boundary_flux(ib,is)
        enddo

        do iy = -1,ny
          do ix = -1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &        bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
              ssna(is)=ssna(is)+sna(ix,iy,0,is)+
     &                        sna(ix,iy,1,is)*na(ix,iy,is)
            endif
          enddo
        enddo

        sum_ssna = sum_ssna + ssna(is)

        sum_sdiv = sum_sdiv + sdiv(is)
        if (nint(zamin(is)).lt.10) then
          write(hlp_frm,'(a)') '(a,9x,a2,i1,4x,3(a,es15.7,4x))'
        else
          write(hlp_frm,'(a)') '(a,9x,a2,i2,3x,3(a,es15.7,4x))'
        endif
        write (70,hlp_frm) 'Balance for equation ',
     &     adjustr(is_codes(is)),nint(zamin(is)),
     &     ' Sdiv =', sdiv(is),' Ssna =', ssna(is),
     &     '     Ssna-Sdiv =', ssna(is) - sdiv(is)
       enddo

       write(70,"(a/33x,4(a,es15.7))") 'Balance over all equations:',
     &      ' sum_Sdiv =', sum_sdiv,' sum_Ssna =', sum_ssna,
     &      ' sum_Ssna-sum_Sdiv =', sum_ssna - sum_sdiv
       close(70)

       call b2wdat_open(70,'Total_boundary_particle_fluxes_modified')  !srv 26.09.17
       sum_sdiv = 0.0_R8
       sum_ssna = 0.0_R8

       do is=0, ns-1
         write (hlp_frm,'(a,i1,a)')
     .    '(a,1x,a2,i',len_of_digits(nint(zamin(is))),')'
         write (70,hlp_frm) 'Total_boundary_particle_flux ',
     &               adjustr(is_codes(is)),nint(zamin(is))

         sdiv(is) = 0.0_R8
         ssna(is) = 0.0_R8

         do ib=1, nbc
          boundary_flux(ib,is)=0.0e0_R8
          IF(BCCHAR(IB).EQ.'S') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              boundary_flux(ib,is) = boundary_flux(ib,is) +
     &         fna_mdf(TOPIX(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),
     &                 TOPIY(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB)),1,is)
            end do
            sdiv(is) = sdiv(is) - boundary_flux(ib,is)
          ELSEIF(BCCHAR(IB).EQ.'N') THEN
            DO IX = 1, BC_LIST_SIZE(IB)
              boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna_mdf(BC_LIST_X(IX,IB),BC_LIST_Y(IX,IB),1,is)
            end do
            sdiv(is) = sdiv(is) + boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).EQ.'W') THEN
           DO IY = 1, BC_LIST_SIZE(IB)
            boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna_mdf(rightix(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),
     &                rightiy(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB)),0,is)
           end do
           sdiv(is) = sdiv(is) - boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).EQ.'E') THEN
           DO IY = 1, BC_LIST_SIZE(IB)
             boundary_flux(ib,is) = boundary_flux(ib,is) +
     &        fna_mdf(BC_LIST_X(IY,IB),BC_LIST_Y(IY,IB),0,is)
           end do
           sdiv(is) = sdiv(is) + boundary_flux(ib,is)
         ELSEIF(BCCHAR(IB).NE.'X' .AND. BCCHAR(IB).NE.'Y') THEN
           CALL XERRAB('B2WDAT: INVALID BOUNDARY OPTION')
         ENDIF
         WRITE(70,"(a,a,a,es15.7)") BCCHAR(IB),
     2         boundary_location(ib),
     &         ' flux =',boundary_flux(ib,is)
        enddo

        do iy = -1,ny
          do ix = -1,nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &        bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
              ssna(is)=ssna(is)+sna(ix,iy,0,is)+
     &                        sna(ix,iy,1,is)*na(ix,iy,is)
            endif
          enddo
        enddo

        sum_ssna = sum_ssna + ssna(is)

        sum_sdiv = sum_sdiv + sdiv(is)
        if (nint(zamin(is)).lt.10) then
          write(hlp_frm,'(a)') '(a,9x,a2,i1,4x,3(a,es15.7,4x))'
        else
          write(hlp_frm,'(a)') '(a,9x,a2,i2,3x,3(a,es15.7,4x))'
        endif
        write (70,hlp_frm) 'Balance for equation ',
     &     adjustr(is_codes(is)),nint(zamin(is)),
     &     ' Sdiv =', sdiv(is),' Ssna =', ssna(is),
     &     '     Ssna-Sdiv =', ssna(is) - sdiv(is)
       enddo

       write(70,"(a/33x,4(a,es15.7))") 'Balance over all equations:',
     &      ' sum_Sdiv =', sum_sdiv,' sum_Ssna =', sum_ssna,
     &      ' sum_Ssna-sum_Sdiv =', sum_ssna - sum_sdiv
       close(70)                                                          !srv 12.04.12 }

      endif              !srv   endif for iout=1 .or. iout=3

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wdat

      contains

      subroutine b2wdat_open(nwrite,first_name)
      integer, intent(in) :: nwrite
      character*(*), intent(in) :: first_name
      character*256 :: filename

      filename=trim(folder)//trim(first_name)//'.dat'
      if(iappend.eq.0) then
        open (nwrite,file=trim(filename))
      else
        open (nwrite,file=trim(filename),position='append')
      end if
      return
      end subroutine b2wdat_open

      end subroutine b2wdat
