      subroutine b2wdat (nx, ny, ns, ismain, ismain0, BoRiS)     !srv 05.10.01
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
      use eirmod_extrab25 
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, ismain0
      real (kind=R8) :: BoRiS
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
      integer ncall, iout, is, ix, iy, ic, ib, iw, iout_iout_ASTRA   !srv 05.10.01 03.04.13 19.04.13
      integer n2, jsep, jxi, jxa, ndat, k, mintopcut, ixx,         !srv 12.07.10
     *        js, nxtl, nxtr, ix1, ix2                             !iyv 04.02.11
      integer use_astra
      integer, save :: ishift=0                                    !srv 05.11.09
      real (kind=R8) ::                                            !srv
     *  Total_current(-1:ny), B0out, Total_flux(-1:ny),            !srv 31.08.00 19.02.02
     *  Total_abscurrent(-1:ny), Total_absflux(-1:ny),
     *  Total_energy_flux_i(-1:ny), Total_energy_flux_e(-1:ny), !srv 23.11.00
     *  wrk2(-1:nx,-1:ny), wrk3(-1:nx,-1:ny),                   !srv 23.11.00
     *  Total_energy_flux(-1:ny), power,                        !srv 23.11.00
     *  kintmp,
     *  sfuadia(-1:ny,0:ns-1), sfvadia(-1:ny,0:ns-1), 
     *  sfvaecrb(-1:ny,0:ns-1), sfna(-1:ny,0:ns-1),
     &  sfna_mdf(-1:ny,0:ns-1), sfn(-1:ny), sfn_mdf(-1:ny),
     &  boundary_flux(1:nbc,0:ns-1), sum_sdiv, sum_ssna,
     &  sdiv(0:ns-1), ssna(0:ns-1),                                !srv 12.04.12
     &  cutll, cutlo, t0
      real (kind=R8) ::  bsp(-1:ny), csp(-1:ny), dsp(-1:ny), dy,
     &     cr(4*(ny+1)+1), dseval, wdim(8*ns)                      !srv 03.04.13
      real (kind=R8) ::
     *     hlp_ppwi, hlp_ppwo, hlp, hlp_s, hlp_neu_peak, s
      real (kind=R8)  ::                                           !iyv 15.10.09 
     *   qe_out_sum, qi_out_sum, q_out_sum,                        !iyv 04.02.11
     *   qe_inn_sum, qi_inn_sum, q_inn_sum, 
     *   qe_oup_sum, qi_oup_sum, q_oup_sum, 
     *   qe_iup_sum, qi_iup_sum, q_iup_sum 
      character chns*3                                             !srv 27.04.01
      character(23) hlp_frm                                        !srv 03.04.13
      character*26 boundary_location(nbcd)
*   ..procedures
      external subini, subend, ipgeti, machsfr
      real (kind=R8)  ::                                           !srv
     *   con, flo, dif                                             !srv
*   ..procedures                                                   !srv
      intrinsic max, abs                                           !srv
      real (kind=R8) ::                                            !srv
     *   hy1, expu2, machsfr                                       !srv
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      expu2(t0) = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
*   ..initialisation
      save ncall, iout, use_astra, hlp_neu_peak
      data ncall/0/, iout/0/, use_astra/0/
      data hlp_neu_peak / 5.0_R8 /                                 !iyv 15.10.09

*-----------------------------------------------------------------------
*.documentation-internal
*
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wdat')
*   ..initialise files and common

*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
*   ..set internal parameters on first (and only) call
      if (ncall.eq.0) then
       call ipgeti('b2wdat_iout', iout)                           !srv 05.10.01
       call ipgeti('b2wdat_ishift', ishift)                       !srv 05.10.01
       call ipgeti('b2mndr_astra', use_astra)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)

* ..initialisation
      if (use_astra.ne.0) then
        iout_iout_ASTRA = 1
      else
        iout_iout_ASTRA = 0
      endif

      if (iout.eq.1) then                                          !srv begin if for iout

#ifdef B25_EIRENE
       open(70,file='eirdiag.dat')
	  write(70,"(a/6i5)") 'nlmpgs natm nmol npls nstra nstsid',
     &                     nlmpgs, natm, nmol, npls, nstra, nstsid
	  write(70,"(a/6i5)") 'natmid nmolid ', natmid, nmolid
	  write(70,"(a/6i5)") 'nnatmi nnmoli ', nnatmi, nnmoli
	 close(70)
      do is = 1,nnatmi                                                      !srv 04.06.09 16.02.12{
	 write (chns,'(i2.2)') is
       call my_out(70,nx,ny,dab2(0:nx+1,0:ny+1,is,1),'dab2'//chns)
       call my_out(70,nx,ny,tab2(0:nx+1,0:ny+1,is,1),'tab2'//chns)
       call my_out(70,nx,ny,rfluxa(0:nx+1,0:ny+1,is,1),'rfluxa'//chns)
       call my_out(70,nx,ny,refluxa(0:nx+1,0:ny+1,is,1),'refluxa'//chns)
       call my_out(70,nx,ny,pfluxa(0:nx+1,0:ny+1,is,1),'pfluxa'//chns)
       call my_out(70,nx,ny,pefluxa(0:nx+1,0:ny+1,is,1),'pefluxa'//chns)
      enddo
      do is = 1,nnmoli                                                      !srv 16.02.12
	 write (chns,'(i2.2)') is
       call my_out(70,nx,ny,dmb2(0:nx+1,0:ny+1,is,1),'dmb2'//chns)          !srv 29.04.10
       call my_out(70,nx,ny,tmb2(0:nx+1,0:ny+1,is,1),'tmb2'//chns)          !srv 29.04.10
       call my_out(70,nx,ny,rfluxm(0:nx+1,0:ny+1,is,1),'rfluxm'//chns)
       call my_out(70,nx,ny,refluxm(0:nx+1,0:ny+1,is,1),'refluxm'//chns)
       call my_out(70,nx,ny,pfluxm(0:nx+1,0:ny+1,is,1),'pfluxm'//chns)
       call my_out(70,nx,ny,pefluxm(0:nx+1,0:ny+1,is,1),'pefluxm'//chns)
      enddo                                                                !srv 04.06.09 }
#endif
      call get_jsep(nx,ny,jxi,jxa,jsep)
      Total_current=0.0_R8
      Total_abscurrent=0.0_R8
      B0out=0.0_R8
      Total_flux=0.0_R8
      Total_absflux=0.0_R8
      Total_energy_flux_i=0.0_R8
      Total_energy_flux_e=0.0_R8
      wrk2=0.0_R8
      wrk3=0.0_R8
      Total_energy_flux=0.0_R8
      call my_outi(70,nx,ny,bottomix,'bottomix')
      call my_outi(70,nx,ny,bottomiy,'bottomiy')
      call my_outi(70,nx,ny,topix,'topix')
      call my_outi(70,nx,ny,topiy,'topiy')
      call my_outi(70,nx,ny,leftix,'leftix')
      call my_outi(70,nx,ny,leftiy,'leftiy')
      call my_outi(70,nx,ny,rightix,'rightix')
      call my_outi(70,nx,ny,rightiy,'rightiy')
      call my_out(70,nx,ny,po,'po')                                !srv
      call my_out(70,nx,ny,ni(-1,-1,0),'ni0')                      !srv
      call my_out(70,nx,ny,ni(-1,-1,1),'ni1')                      !srv
      wrk0=0.0e0_R8
	do iy=0,ny
	 do ix=-1,nx
	   wrk0(ix,iy)=-(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *                (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy))) !srv 08.05.07
     *                *2.0e0_R8
	 enddo
	enddo
      call my_out(70,nx,ny,wrk0,'Eyf')                              !srv 15.06.04 }

	wrk0=0.0e0_R8                                                !iyv 15.10.09 {
	do iy=0,ny
	 do ix=-1,nx
         if(topiy(ix,iy).ne.ny+1) then
         wrk0(ix,iy)=-(po(topix(ix,iy),topiy(ix,iy))-
     *                  po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *                (hy1(ix,iy)+0.5e0_R8*
     *                (hy1(bottomix(ix,iy),bottomiy(ix,iy))) +
     *                 hy1(topix(ix,iy),topiy(ix,iy))) 
         else
	   wrk0(ix,iy)=-(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *                (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy))) !
     *                *2.0e0_R8
         endif	                     
	 enddo
	enddo
      call my_out(70,nx,ny,wrk0,'Eyc')                              

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
     *   wrk0(ix,iy)=abs(hz(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)* 
     *              (wrk1(topix(ix,iy),topiy(ix,iy))-
     *               wrk1(bottomix(ix,iy),bottomiy(ix,iy)))/
     *                (hy1(ix,iy)+0.5e0_R8*
     *                (hy1(bottomix(ix,iy),bottomiy(ix,iy))) +
     *                 hy1(topix(ix,iy),topiy(ix,iy)))) 
	  enddo
      enddo
      call my_out(70,nx,ny,wrk0,'Wsh')                              !iyv 15.10.09 }

      do is = 0,ns-1
	 write (chns,'(i3.3)') is
       call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, wrk0)
       call my_out(70,nx,ny,wrk0,'pb'//chns)
       call my_out(70,nx,ny,na(-1,-1,is),'na'//chns)             !srv
       call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)             !srv
       call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)             !srv
	 wrk0 = sqrt((rza(:,:,is)*te+ti)/(mp*am(is)))              !srv 24.03.10 {
       call my_out(70,nx,ny,wrk0,'usound'//chns)
	 wrk1 = ua(:,:,is)/wrk0
       call my_out(70,nx,ny,wrk1,'ua_M_'//chns)                  !srv 24.03.10 }
       call my_out(70,nx,ny,uadia(-1,-1,0,is),'uadiax'//chns)    !srv
       call my_out(70,nx,ny,uadia(-1,-1,1,is),'uadiay'//chns)    !srv
       call my_out(70,nx,ny,vaecrb(-1,-1,0,is),'vaecrbx'//chns)            !srv
       call my_out(70,nx,ny,vaecrb(-1,-1,1,is),'vaecrby'//chns)            !srv
       call my_out(70,nx,ny,vadia(-1,-1,0,is),'vadiax'//chns)    !srv
       call my_out(70,nx,ny,vadia(-1,-1,1,is),'vadiay'//chns)    !srv
       call my_out(70,nx,ny,wadia(-1,-1,0,is),'wadiax'//chns)    !srv
       call my_out(70,nx,ny,wadia(-1,-1,1,is),'wadiay'//chns)    !srv
       call my_out(70,nx,ny,fna(-1,-1,0,is),'fnax'//chns)        !srv
       call my_out(70,nx,ny,fna(-1,-1,1,is),'fnay'//chns)        !srv
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
     *                     'fna_nodriftx'//chns)
       call my_out(70,nx,ny,fna_nodrift(-1,-1,1,is),
     *                     'fna_nodrifty'//chns)
               
        call my_out(70,nx,ny,cdna(-1,-1,0,is),'cdnax'//chns)      !srv
        call my_out(70,nx,ny,cdna(-1,-1,1,is),'cdnay'//chns)      !srv
        call my_out(70,nx,ny,cdpa(-1,-1,0,is),'cdpax'//chns)      !srv
        call my_out(70,nx,ny,cdpa(-1,-1,1,is),'cdpay'//chns)      !srv
        call my_out(70,nx,ny,cddi(-1,-1,0,is),'cddix'//chns)      !srv 13.10.06
        call my_out(70,nx,ny,cddi(-1,-1,1,is),'cddiy'//chns)      !srv 13.10.06
        call my_out(70,nx,ny,cvsa_an(-1,-1,0,is),'cvsax_an'//chns)
        call my_out(70,nx,ny,cvsa_cl(-1,-1,0,is),'cvsax_cl'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,0,is),'cvsax'//chns)
        call my_out(70,nx,ny,cvsa(-1,-1,1,is),'cvsay'//chns)
        call my_out(70,nx,ny,kinrgy(-1,-1,is),'kinrgy'//chns)     !srv 08.10.04

       wrk0=0.0d0
       do iy = -1,ny
         do ix = -1, nx
           if (leftix(ix,iy).ne.-2) wrk0(ix,iy)=
     1       wbbl(ix,iy,2)/wbbl(ix,iy,3)
     4       *vaecrb(ix,iy,0,is)*
     5       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/        !srv 01.03.00
     6       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))           !srv 01.03.00
     7       *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
         enddo
       enddo
       call my_out(70,nx,ny,wrk0(-1,-1),'fvaecrbx'//chns)        !srv

       wrk1=0.0d0
       do iy = -1,ny
         do ix = -1, nx
           if (leftix(ix,iy).ne.-2) wrk1(ix,iy)=
     1       wbbl(ix,iy,2)/wbbl(ix,iy,3)
     4       *vadia(ix,iy,0,is)*
     5       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/        !srv 01.03.00
     6       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))           !srv 01.03.00
     7       *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
         enddo
       enddo
       call my_out(70,nx,ny,wrk1(-1,-1),'fvadiax'//chns)         !srv

       wrk0=fna(:,:,0,is)-wrk0-wrk1
       call my_out(70,nx,ny,wrk0(-1,-1),'fdifx'//chns)           !srv

       wrk0=0.0d0
       do iy = -1, ny
         do ix = -1, nx
           if (bottomiy(ix,iy).ne.-2)
     >      wrk0(ix,iy)=vaecrb(ix,iy,1,is)*(na(ix,iy,is)
     &          +na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8*
     5       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/    !srv 01.03.00
     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 01.03.00
         enddo
       enddo
       call my_out(70,nx,ny,wrk0(-1,-1),'fvaecrby'//chns)        !srv
       wrk1=0.0d0
       do iy = -1, ny
         do ix = -1, nx
           if (bottomiy(ix,iy).ne.-2)
     >       wrk1(ix,iy)=vadia(ix,iy,1,is)*(na(ix,iy,is)
     &          +na(bottomix(ix,iy),bottomiy(ix,iy),is))/2.0_R8*
     5       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/    !srv 01.03.00
     6     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 01.03.00
         enddo
       enddo
       call my_out(70,nx,ny,wrk1(-1,-1),'fvadiay'//chns)         !srv
       wrk0=fna(:,:,1,is)-wrk0-wrk1
       call my_out(70,nx,ny,wrk0(-1,-1),'fdify'//chns)           !srv

       call my_out(70,nx,ny,fmo(-1,-1,0,is),'fmox'//chns)           !srv 24.02.10
       call my_out(70,nx,ny,fmo(-1,-1,1,is),'fmoy'//chns)

       open(70,file='Total_flux_pm'//chns//'.dat')               !srv 19.02.02 {
       write (70,*) '  N      Total_flux    Total_|flux| '
       mintopcut = minval(topcut(1:nncut))
         do iy = mintopcut,0,-1 
           Total_flux(iy) = 0.0d0
           Total_absflux(iy) = 0.0d0
           do ix=-1,nx
	      if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
             Total_flux(iy) = Total_flux(iy) + fmo(ix,iy,1,is)
             Total_absflux(iy)=Total_absflux(iy)+abs(fmo(ix,iy,1,is))
	      endif
           enddo                    
           write (70,'(1x,i4,2(2x,es14.7))') iy, Total_flux(iy),
     &           Total_absflux(iy)
         enddo
       close (70)                                                  !srv 19.02.02 }

       call my_out(70,nx,ny,rqahe(-1,-1,is),'rqahe'//chns)       !srv
       call my_out(70,nx,ny,rqrad(-1,-1,is),'rqrad'//chns)       !xpb
       call my_out(70,nx,ny,rqbrm(-1,-1,is),'rqbrm'//chns)       !xpb
       call my_out(70,nx,ny,pcca(-1,-1,0,is),'pccax'//chns)           !srv 11.01.13
       call my_out(70,nx,ny,pcca(-1,-1,1,is),'pccay'//chns)           !srv 11.01.13

      enddo

      call my_out(70,nx,ny,pccm(-1,-1,0),'pccmx')           !srv 11.01.13
      call my_out(70,nx,ny,pccm(-1,-1,1),'pccmy')           !srv 11.01.13
      call my_out(70,nx,ny,fhe(-1,-1,0),'fhex')                    !srv 23.11.00
      call my_out(70,nx,ny,fhe(-1,-1,1),'fhey')                    !srv 23.11.00
      call my_out(70,nx,ny,fhi(-1,-1,0),'fhix')                    !srv 23.11.00
      call my_out(70,nx,ny,fhi(-1,-1,1),'fhiy')                    !srv 23.11.00

      call my_out(70,nx,ny,fhe_mdf(-1,-1,0),'fhe_mdfx')            !srv 13.10.06 {
      call my_out(70,nx,ny,fhe_mdf(-1,-1,1),'fhe_mdfy')
      call my_out(70,nx,ny,fhePSch(-1,-1,0),'fhePSchx')
      call my_out(70,nx,ny,fhePSch(-1,-1,1),'fhePSchy')
      call my_out(70,nx,ny,fhi_mdf(-1,-1,0),'fhi_mdfx')
      call my_out(70,nx,ny,fhi_mdf(-1,-1,1),'fhi_mdfy')
      call my_out(70,nx,ny,fhiPSch(-1,-1,0),'fhiPSchx')
      call my_out(70,nx,ny,fhiPSch(-1,-1,1),'fhiPSchy')            !srv 13.10.06 }

      call my_out(70,nx,ny,ne,'ne')                                !srv
      call my_out(70,nx,ny,te,'te')                                !srv
      call my_out(70,nx,ny,te/ev,'te_eV')
      call my_out(70,nx,ny,ti,'ti')                                !srv
      call my_out(70,nx,ny,ti/ev,'ti_eV')
      call my_out(70,nx,ny,fch(-1,-1,0),'fchx')                    !srv
      call my_out(70,nx,ny,fch(-1,-1,1),'fchy')                    !srv
      call my_out(70,nx,ny,fch_p(-1,-1,0),'fchx_p')                !srv 10.06.08
      call my_out(70,nx,ny,fch_p(-1,-1,1),'fchy_p')                !srv 10.06.08
      call my_out(70,nx,ny,fchvispar(-1,-1,0),'fchvisparx')
      call my_out(70,nx,ny,fchvispar(-1,-1,1),'fchvispary')
      call my_out(70,nx,ny,fchvisper(-1,-1,0),'fchvisperx')
      call my_out(70,nx,ny,fchvisper(-1,-1,1),'fchvispery')
      call my_out(70,nx,ny,fchvisq(-1,-1,0),'fchvisqx')
      call my_out(70,nx,ny,fchvisq(-1,-1,1),'fchvisqy')
      call my_out(70,nx,ny,fchinert(-1,-1,0),'fchinertx')          !srv 07.07.99
      call my_out(70,nx,ny,fchinert(-1,-1,1),'fchinerty')          !srv 07.07.99
      call my_out(70,nx,ny,fchdia(-1,-1,0),'fchdiax')
      call my_out(70,nx,ny,fchdia(-1,-1,1),'fchdiay')
      call my_out(70,nx,ny,fchin(-1,-1,0),'fchinx')
      call my_out(70,nx,ny,fchin(-1,-1,1),'fchiny')
      call my_out(70,nx,ny,fchanml(-1,-1,0),'fchanmlx')
      call my_out(70,nx,ny,fchanml(-1,-1,1),'fchanmly')
      call b2tlnl (nx, ny, 0, te, ti, ne, wrk0)                   !srv 02.07.13 {
      call my_out(70,nx,ny,wrk0,'lnlam_Br')
      call b2tlnl (nx, ny, 1, te, ti, ne, wrk0)
      call my_out(70,nx,ny,wrk0,'lnlam_We_ii')
      call b2tlnl (nx, ny, 2, te, ti, ne, wrk0)
      call my_out(70,nx,ny,wrk0,'lnlam_We_ee')
      call b2tlnl (nx, ny, 3, te, ti, ne, wrk0)                   !srv 02.07.13 }
      call my_out(70,nx,ny,wrk0,'lnlam_We_ei')
      call my_out(70,-1,ny,conn,'conn')                            !srv
      call my_out(70,nx,ny,csig(-1,-1,0),'csig_x')                 !srv
      call my_out(70,nx,ny,csig(-1,-1,1),'csig_y')                 !srv
      call my_out(70,nx,ny,csig_an(-1,-1,0),'csig_anx')            !srv
      call my_out(70,nx,ny,csig_an(-1,-1,1),'csig_any')            !srv
      if (ismain0.ne.ismain) then
        call my_out(70,nx,ny,csigin(-1,-1,0,ismain,ismain0),'csig_inx') !srv
        call my_out(70,nx,ny,csigin(-1,-1,1,ismain,ismain0),'csig_iny') !srv
      endif
      call my_out(70,nx,ny,vsa0(-1,-1,0),'vsa0')                   !srv
      call my_out(70,nx,ny,vsa0(-1,-1,1),'vsa1')                   !srv
      call my_out(70,nx,ny,cvsa(-1,-1,0,ismain),'cvsa')            !srv
      call my_out(70,nx,ny,cvsa_cl(-1,-1,0,ismain),'cvsa_cl')      !srv
      call my_out(70,nx,ny,calf(-1,-1,0),'calfx')                 !srv
      call my_out(70,nx,ny,calf(-1,-1,1),'calfy')                 !srv
      call my_out(70,nx,ny,chce(-1,-1,0),'chcex')                 !srv
      call my_out(70,nx,ny,chce(-1,-1,1),'chcey')                 !srv
      call my_out(70,nx,ny,chci(-1,-1,0),'chcix')                 !srv
      call my_out(70,nx,ny,chci(-1,-1,1),'chciy')                 !srv

      call my_out(70,nx,ny,lengthy,'lengthy')                      !srv
	wrk0=0.0d0                                                   !srv 18.03.09 26.04.09 {
	do ix=-1,nx
	 wrk0(ix,-1)=hy1(ix,-1)
	 do iy=0,ny
	  wrk0(ix,iy)=wrk0(ix,iy-1)+hy1(ix,iy)
	 enddo
	enddo
      call my_out(70,nx,ny,wrk0,'lengthy1')                              !srv

      call my_out(70,nx,ny,len_surf,'len_surf')

      call my_out(70,nx,ny,crdy_surf,'crdy_surf')
      call my_out(70,nx,ny,crdy_c_cell,'crdy_c_cell')             !srv 18.03.09 26.04.09}

       call my_out(70,nx,ny,gs(-1,-1,0),'gsx')
       call my_out(70,nx,ny,gs(-1,-1,1),'gsy')
       call my_out(70,nx,ny,gs(-1,-1,2),'gsz')

       call my_out(70,nx,ny,wbbl(-1,-1,0),'bbxl')
       call my_out(70,nx,ny,wbbl(-1,-1,2),'bbzl')
       call my_out(70,nx,ny,wbbl(-1,-1,3),'bbl')
       call my_out(70,nx,ny,wbbr(-1,-1,0),'bbxr')
       call my_out(70,nx,ny,wbbr(-1,-1,2),'bbzr')
       call my_out(70,nx,ny,wbbr(-1,-1,3),'bbr')
       call my_out(70,nx,ny,wbbv(-1,-1,0),'bbxv')
       call my_out(70,nx,ny,wbbv(-1,-1,2),'bbzv')
       call my_out(70,nx,ny,wbbv(-1,-1,3),'bbv')

      wrk0=vol/hx**2                                               !srv
      call my_out(70,nx,ny,wrk0,'voldhx2')                         !srv
      wrk0=vol/(hy*qz(:,:,1))**2                                   !srv
      call my_out(70,nx,ny,wrk0,'voldhy2')                         !srv
      call my_out(70,nx,ny,OnedBsq,'OnedBsq')                      !srv
      call my_out(70,nx+1,ny+1,OnedBsqc,'OnedBsqc')                !srv

      call my_out(70,nx,ny,resco(-1,-1,0),'resco0')                !srv
      call my_out(70,nx,ny,resco(-1,-1,1),'resco1')                !srv
      call my_out(70,nx,ny,resmo(-1,-1,0),'resmo0')                !srv
      call my_out(70,nx,ny,resmo(-1,-1,1),'resmo1')                !srv
      call my_out(70,nx,ny,respo,'respo')                          !srv
      call my_out(70,nx,ny,reshe,'reshe')                          !srv
      call my_out(70,nx,ny,reshi,'reshi')                          !srv
      call my_out(70,nx,ny,she(-1,-1,0),'she0')                    !srv
      call my_out(70,nx,ny,she(-1,-1,1),'she1')                    !srv
      call my_out(70,nx,ny,she(-1,-1,2),'she2')                    !srv
      call my_out(70,nx,ny,she(-1,-1,3),'she3')                    !srv
      call my_out(70,nx,ny,shi(-1,-1,0),'shi0')                    !srv
      call my_out(70,nx,ny,shi(-1,-1,1),'shi1')                    !srv
      call my_out(70,nx,ny,shi(-1,-1,2),'shi2')                    !srv
      call my_out(70,nx,ny,shi(-1,-1,3),'shi3')                    !srv
      call my_out(70,nx,ny,sch(-1,-1,0),'sch0')                    !srv
      call my_out(70,nx,ny,sch(-1,-1,1),'sch1')                    !srv
      call my_out(70,nx,ny,sch(-1,-1,2),'sch2')                    !srv
      call my_out(70,nx,ny,sch(-1,-1,3),'sch3')                    !srv
      call my_out(70,nx,ny,sna(-1,-1,0,1),'sna10')                 !srv
      call my_out(70,nx,ny,sna(-1,-1,1,1),'sna11')                 !srv
      wrk0=sna(:,:,0,ismain)+sna(:,:,1,ismain)*na(:,:,ismain)
      call my_out(70,nx,ny,wrk0(-1,-1),'sour_ni')                  !srv

       open(70,file='Total_current'//'.dat')                       !srv 31.08.00 23.07.10 {
        write (70,*) '  N      Total_current  Total_|current| '
        mintopcut = minval(topcut(1:nncut))
        do iy = mintopcut,0,-1
          Total_current(iy) = 0.0_R8
          Total_abscurrent(iy) = 0.0_R8
          do ix = -1, nx 
	     if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
            Total_current(iy) = Total_current(iy) + fch(ix,iy,1)
            Total_abscurrent(iy)=Total_current(iy)+abs(fch(ix,iy,1))
	     endif
          enddo
          write (70,'(1x,i4,2(2x,es14.7))') iy, Total_current(iy),
     &                                      Total_abscurrent(iy)
        enddo
        ix = jxa ; iy = 1
        write (70,'(1x,a,i4,a,i4,a,g14.7)')
     &         'potential in (',ix,',',iy,')=',po(ix,iy)
       close (70)

       open(70,file='Total_energy_flux_i'//'.dat')
        write (70,*) '  N      Total_energy_flux_i '

        do iy = mintopcut,0,-1
          Total_energy_flux_i(iy) = 0.0_R8
          do ix = -1, nx
	     if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
            Total_energy_flux_i(iy) = 
     &       Total_energy_flux_i(iy) + fhi(ix,iy,1)
!            do is = 0, ns-1                                       ! ask Rozhansky
!              kintmp = 
!     .         (kinrgy(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
!     .          kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
!     .         (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
!              Total_energy_flux_i(iy) = 
!     &         Total_energy_flux_i(iy) + kintmp*(1.0_R8-BoRiS)
!            enddo
	     endif
          enddo
          write (70,'(1x,i4,2x,es14.7)') iy, Total_energy_flux_i(iy)
        enddo
       close (70)

       open(70,file='Total_energy_flux_e'//'.dat')
        write (70,*) '  N      Total_energy_flux_e '

        do iy = mintopcut,0,-1
          Total_energy_flux_e(iy) = 0.0_R8
          do ix = -1, nx
	     if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
            Total_energy_flux_e(iy) = 
     &       Total_energy_flux_e(iy) + fhe(ix,iy,1)
	     endif
          enddo
          write (70,'(1x,i4,2x,es14.7)') iy, Total_energy_flux_e(iy)
        enddo
       close (70)

       open(70,file='Total_energy_flux.dat')
        write (70,*) '  N      Total_energy_flux '

        do iy = mintopcut,0,-1
         Total_energy_flux(iy) = Total_energy_flux_e(iy)
     *                          +Total_energy_flux_i(iy)
          write (70,'(1x,i4,2x,es14.7)') iy, Total_energy_flux(iy)
        enddo
       close (70)                                                  !srv 31.08.00  23.07.10 }

       open(70,file='info.dat')                                     !srv
        write(70,"(' Number of outer midplane cell =',i4)") nmdpl   !srv 24.03.06 }
	  write (70,*) ' nncut    =', nncut, ' nnreg =',nnreg         !srv
        write (70,*) ' leftcut  =',(leftcut(ix),ix=1,nncut)         !srv
        write (70,*) ' rightcut =',(rightcut(ix),ix=1,nncut)        !srv
        write (70,*) ' topcut   =',(topcut(ix),ix=1,nncut)          !srv
     &           , ' is the number of first outer cells in SOL'
        write (70,*) ' bottomcut=',(bottomcut(ix),ix=1,nncut)       !srv
        write (70,*) ' Number of dead (or isolated) regions NNISO',NNISO
        write (70,*) ' Column number of bottom left corner ', 
     &   'of the dead region NXISO1', (NXISO1(IS),IS=1,NNISO)
        write (70,*) ' Column number of top right corner ', 
     &   'of the dead region NXISO2', (NXISO2(IS),IS=1,NNISO)
        write (70,*) ' Row number of bottom left corner ', 
     &   'of the dead region NYISO1', (NYISO1(IS),IS=1,NNISO)
        write (70,*) ' Row number of top right corner ', 
     &   'of the dead region NYISO2', (NYISO2(IS),IS=1,NNISO)
	 close (70)                                                   !srv
*   ..write physics parameters after changing                      !srv 19.12.07 {
      open(70,file='b2mn_chk.prt')
       call b2wfcp (70, 0, ns)
      close(70)                                                    !srv 19.12.07 }
*   ..write "experimental" data                 
	ndat = 1
	cr(ndat) = crdy_c_cell(nmdpl,-1)
	do iy = -1, ny-1
	 dy = (crdy_c_cell(nmdpl,iy+1)-crdy_c_cell(nmdpl,iy))/4.0_R8
	 do k=1,4
	  cr(ndat+1)=cr(ndat)+dy
        ndat = ndat+1
	 enddo
	enddo

      call dspline(ny+2,crdy_c_cell(nmdpl,:),ne(nmdpl,:),bsp,csp,dsp)
      open(70,file='splne.dat')
	 do iy=1,ndat
	  write(70,"(1x,2(es14.6))") cr(iy),
     &  dseval(ny+2,cr(iy),crdy_c_cell(nmdpl,:),ne(nmdpl,:),bsp,csp,dsp)
	 enddo
	close(70)

      call dspline(ny+2,crdy_c_cell(nmdpl,:),te(nmdpl,:),bsp,csp,dsp)
      open(70,file='splte.dat')
	 do iy=1,ndat
	  write(70,"(1x,2(es14.6))") cr(iy),
     &  dseval(ny+2,cr(iy),crdy_c_cell(nmdpl,:),te(nmdpl,:),bsp,csp,dsp)
	 enddo
	close(70)

      call dspline(ny+2,crdy_c_cell(nmdpl,:),ti(nmdpl,:),bsp,csp,dsp)
      open(70,file='splti.dat')
	 do iy=1,ndat
	  write(70,"(1x,2(es14.6))") cr(iy),
     &  dseval(ny+2,cr(iy),crdy_c_cell(nmdpl,:),ti(nmdpl,:),bsp,csp,dsp)
	 enddo
	close(70)

      open(70,file='expne.dat')
	 do iy=-1,ny
	  write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), ne(nmdpl,iy)
	 enddo
	close(70)
      open(70,file='expte.dat')
	 do iy=-1,ny
	  write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), te(nmdpl,iy)
	 enddo
	close(70)
      open(70,file='expti.dat')
	 do iy=-1,ny
	  write(70,"(1x,2(es14.6))") crdy_c_cell(nmdpl,iy), ti(nmdpl,iy)
	 enddo
	close(70)
*   ..calculate total energy flux
	 iy = 0
       power=0.0_R8
       do ix = -1, nx
	  if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
	    power = power + (fhe(ix,iy,1)+fhi(ix,iy,1))
	  endif
	 enddo
	 open(70,file='exppw.dat')
	  write(70,"(es14.6)") power
	 close(70)
*   ..calculate integrated particle fluxes for checking
	 open(70,file='Total_particle_fluxes_in_Core.dat')
       do is = 0,ns-1
	 write (70,"(a,1x,a2,i1)") 'Total_particle_flux ',
     &     adjustr(is_codes(is)),nint(zamin(is))
	  write(70,"(' iy ','   sfuadia  ','   sfvadia  ','   sfvaecrb ',
     &  '    sfna   ','    sfna_mdf ')")
       do iy = ny,-1, -1
        sfuadia(iy,is)  =0.0_R8
        sfvadia(iy,is)  =0.0_R8
        sfvaecrb(iy,is) =0.0_R8
        sfna(iy,is)     =0.0_R8
        sfna_mdf(iy,is) =0.0_R8
        do ix = -1, nx
	   if(mod(region(ix,iy,0),4).eq.1.and.bottomiy(ix,iy).ne.-2) then  ! for CORE  !srv 01.02.12
	    sfuadia(iy,is) = sfuadia(iy,is) + uadia(ix,iy,1,is)*
     5     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/     
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))/2.0_R8
	    sfvadia(iy,is) = sfvadia(iy,is) + vadia(ix,iy,1,is)*
     5     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/     
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))/2.0_R8
	    sfvaecrb(iy,is) = sfvaecrb(iy,is) + vaecrb(ix,iy,1,is)*
     5     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/     
     6      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*
     &     (na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))/2.0_R8
	    sfna(iy,is) = sfna(iy,is) + fna(ix,iy,1,is)
	    sfna_mdf(iy,is) = sfna_mdf(iy,is) + fna_mdf(ix,iy,1,is)
	   endif
	  enddo
	  write(70,"(i4,5es15.7)") iy, sfuadia(iy,is), sfvadia(iy,is),
     &     sfvaecrb(iy,is), sfna(iy,is), sfna_mdf(iy,is)
	 enddo
      enddo

	 write (70,"(a)") 'Total particle flux over all particles '      !srv 12.04.12 {
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

       open(70,file='Total_boundary_particle_fluxes.dat')

	 sum_sdiv = 0.0_R8
	 sum_ssna= 0.0_R8

       do is=0, ns-1
	write (70,"(a,1x,a2,i1)") 'Total_boundary_particle_flux ',
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
     .      ' on boundary cell list ',IB
         else
           write(boundary_location(ib),'(a4,i4,a6,i4,a4,i4)')
     .      ' at ',BCPOS(IB),' from ',BCSTART(IB),' to ',BCEND(IB)
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

	sum_sdiv = sum_sdiv +sdiv(is)
	write (70,"(a,9x,a2,i1,3(4x,a,es15.7))")'Balance for equation ',
     &     adjustr(is_codes(is)),nint(zamin(is)),
     &  ' Sdiv =',sdiv(is), ' Ssna =', ssna(is),
     &  '     Ssna-Sdiv =', ssna(is) - sdiv(is)
       enddo

       write(70,"(a/33x,4(a,es15.7))") 'Balance over all equations:',
     &      ' sum_Sdiv =', sum_sdiv, ' sum_Ssna =', sum_ssna,
     &      ' sum_Ssna-sum_Sdiv =', sum_ssna - sum_sdiv 
	close(70)

        open(70,file='Total_boundary_particle_fluxes_modified.dat')
	sum_sdiv = 0.0_R8
	sum_ssna= 0.0_R8

        do is=0, ns-1
	 write (70,"(a,1x,a2,i1)") 'Total_boundary_particle_flux ',
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
           CALL XERRAB
     1         ('B2WDAT: INVALID OPTION')
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

	sum_sdiv = sum_sdiv +sdiv(is)
	write (70,"(a,9x,a2,i1,3(4x,a,es15.7))")'Balance for equation ',
     &     adjustr(is_codes(is)),nint(zamin(is)),
     &  ' Sdiv =',sdiv(is), ' Ssna =', ssna(is),
     &  '     Ssna-Sdiv =', ssna(is) - sdiv(is)
       enddo

       write(70,"(a/33x,4(a,es15.7))") 'Balance over all equations:',
     &      ' sum_Sdiv =', sum_sdiv, ' sum_Ssna =', sum_ssna,
     &      ' sum_Ssna-sum_Sdiv =', sum_ssna - sum_sdiv 
	close(70)                                                          !srv 12.04.12 }

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
	 do iy=-1,ny
	  do ix=-1,nx
         wrk0(ix,iy)=hy1(ix,iy)
	  enddo
	 enddo
       call my_out(70,nx,ny,wrk0,'hy1')                             !srv 07.11.12                            
       call my_out(70,nx,ny,vol,'vol')
       call my_out(70,nx,ny,pbs(-1,-1,0),'pbsx')
!       call my_out(70,nx,ny,pbs(-1,-1,1),'pbsy')

       call my_out(70,nx,ny,gs(-1,-1,0),'gsx')
       call my_out(70,nx,ny,gs(-1,-1,1),'gsy')
       call my_out(70,nx,ny,gs(-1,-1,2),'gsz')

!	 call my_outi(70,nx,ny,region(-1,-1,0),'region0')
!	 call my_outi(70,nx,ny,region(-1,-1,1),'region1')
!	 call my_outi(70,nx,ny,region(-1,-1,2),'region2')

      endif              !srv   endif for iout=1

      if(iout.eq.1 .or. iout_iout_ASTRA.eq.1) then                                !srv 19.04.13
!______________________________________________________________________________ !iyv 15.10.09 {
! target energy fluxes                                  
      wrk0=0.0e0_R8
	wrk1=0.0e0_R8

      if(nncut.eq.2)  then           !iyv 04.02.11 
	  call get_nxt(nx, nxtl, nxtr) ! position of additional targets in DN configuration
        ix1=nxtl    !inner upper
        ix2=nxtr+1  !outer upper 
      endif

	do iy=-1,ny
        do is = 0, ns-1          
          if(.not.is_neutral(is)) then          !iyv 07.02.11
! pr_inner
          wrk1(-1,iy)=wrk1(-1,iy)+na(-1,iy,is)*(rza(-1,iy,is)*te(-1,iy)+
     &                                                 ti(-1,iy))
! pr_outer
          wrk1(0,iy)=wrk1(0,iy)+na(nx,iy,is)*(rza(nx,iy,is)*te(nx,iy)+
     &                                                 ti(nx,iy))
! ro_inner
          wrk1(1,iy) = wrk1(1,iy)+am(is)*mp*na(-1,iy,is)
! ro_outer
          wrk1(2,iy) = wrk1(2,iy)+am(is)*mp*na(nx,iy,is)

          if(nncut.eq.2) then !{    !iyv 04.02.11
! pr_inner_upper
            wrk1(3,iy)=wrk1(3,iy)+na(ix1,iy,is)*(rza(ix1,iy,is)*
     &                                        te(ix1,iy)+ti(ix1,iy))
! pr_outer_upper
            wrk1(4,iy)=wrk1(4,iy)+na(ix2,iy,is)*(rza(ix2,iy,is)*
     &                                        te(ix2,iy)+ti(ix2,iy))
! ro_inner_upper
            wrk1(5,iy) = wrk1(5,iy)+am(is)*mp*na(ix1,iy,is)
! ro_outer_upper
            wrk1(6,iy) = wrk1(6,iy)+am(is)*mp*na(ix2,iy,is)
          endif
         endif
        enddo

! qe inner
        wrk0(-1,iy)=hy1(-1,iy)/hy(-1,iy)*abs(bb(-1,iy,0)/bb(-1,iy,3))*    !srv 19.04. 13
     &   ( te(-1,iy)*(2.0e0_R8)+qe*po(-1,iy)*(1.0_R8-gammae) )*
     &   ( ne(-1,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(-1,iy)/me)* 
     &    expu2(max(-50.0_R8, -qe*(po(-1,iy)-0.0_R8)/te(-1,iy))) )   !instead of 0.0_R8 should be potpar(ib,2)
! qe outer
         wrk0(0,iy)=hy1(nx,iy)/hy(nx,iy)*abs(bb(nx,iy,0)/bb(nx,iy,3))*    !srv 19.04. 13
     &     ( te(nx,iy)*(2.0e0_R8)+qe*po(nx,iy)*(1.0_R8-gammae) )*
     &     ( ne(nx,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(nx,iy)/me)* 
     &    expu2(max(-50.0_R8, -qe*(po(nx,iy)-0.0_R8)/te(nx,iy))) )   !instead of 0.0_R8 should be potpar(ib,2)

         do is = 0, ns-1
           if(.not.is_neutral(is)) then
              wrk0(1,iy)=wrk0(1,iy)+na(-1,iy,is)
              wrk0(2,iy)=wrk0(2,iy)+na(nx,iy,is)
           endif
	   enddo 	  		  	  
! qi inner
         wrk0(3,iy)=2.5e0_R8*
     &               hy1(-1,iy)/hy(-1,iy)*bb(-1,iy,0)/bb(-1,iy,3)*
     &               ti(-1,iy)*wrk0(1,iy)*sqrt(wrk1(-1,iy)/wrk1(1,iy))
! qi outer
	   wrk0(4,iy)=2.5e0_R8*
     &               hy1(nx,iy)/hy(nx,iy)*bb(nx,iy,0)/bb(nx,iy,3)*
     &               ti(nx,iy)*wrk0(2,iy)*sqrt(wrk1(0,iy)/wrk1(2,iy))
 
       if(nncut.eq.2) then !{      !iyv 04.02.11

! qe inner upper
        wrk0(7,iy)=hy1(ix1,iy)/hy(ix1,iy)*                             !srv 19.04. 13
     &     abs(bb(ix1,iy,0)/bb(ix1,iy,3))*
     &   ( te(ix1,iy)*(2.0e0_R8)+qe*po(ix1,iy)*(1.0_R8-gammae) )*
     &   ( ne(ix1,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(ix1,iy)/me)* 
     &    expu2(max(-50.0_R8, -qe*(po(ix1,iy)-0.0_R8)/te(ix1,iy))) )   !instead of 0.0_R8 should be potpar(ib,2)
! qe outer upper
        wrk0(8,iy)=hy1(ix2,iy)/hy(ix2,iy)*                             !srv 19.04. 13
     &     abs(bb(ix2,iy,0)/bb(ix2,iy,3))*
     &   ( te(ix2,iy)*(2.0e0_R8)+qe*po(ix2,iy)*(1.0_R8-gammae) )*
     &   ( ne(ix2,iy)*sqrt(1.0_R8/(2.0_R8*pi))*sqrt(te(ix2,iy)/me)* 
     &    expu2(max(-50.0_R8, -qe*(po(ix2,iy)-0.0_R8)/te(ix2,iy))) )   !instead of 0.0_R8 should be potpar(ib,2)

         do is = 0, ns-1
           if(.not.is_neutral(is)) then
              wrk0(9,iy) =wrk0(9,iy) +na(ix1,iy,is)
              wrk0(10,iy)=wrk0(10,iy)+na(ix2,iy,is)
           endif
	   enddo 	  		  	  
! qi inner upper
         wrk0(11,iy)=2.5e0_R8*
     &               hy1(ix1,iy)/hy(ix1,iy)*bb(ix1,iy,0)/bb(ix1,iy,3)*
     &               ti(ix1,iy)*wrk0(9,iy)*sqrt(wrk1(3,iy)/wrk1(5,iy))
! qi outer upper
	   wrk0(12,iy)=2.5e0_R8*
     &               hy1(ix2,iy)/hy(ix2,iy)*bb(ix2,iy,0)/bb(ix2,iy,3)*
     &               ti(ix2,iy)*wrk0(10,iy)*sqrt(wrk1(4,iy)/wrk1(6,iy))
       
	  endif
	enddo

! next part is copied from b2usrtrc                                                      
!________________________________________________________________	        

      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        s=gs(0,iy,0)
        hlp=0.0_R8
        do is=0,ns-1 !{
          hlp=hlp+abs(fna(0,iy,0,is)*(
     .     (rpt(0,iy,is)*hx(-1,iy)+rpt(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(0,iy,is)*hx(-1,iy)+kinrgy(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy))))
        end do !}
        hlp=hlp/s
        wrk0(5,iy) = hlp + wrk0(-1,iy)+wrk0(3,iy) 
        hlp_s=hlp_s+s
      end do !}
#ifdef B25_EIRENE
      js=ltns(1)+nlim
      do iy=0,ny-1 !{
        wrk0(5,iy)=wrk0(5,iy)+(wldnek(js,0)+wldnep(js,0))*
     &                                hlp_neu_peak/hlp_s
      enddo
#endif
      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        hlp=0.0_R8
        s=gs(nx,iy,0)
        do is=0,ns-1 !{
          hlp=hlp+abs(fna(nx,iy,0,is)*(
     .     (rpt(nx,iy,is)*hx(nx-1,iy)+rpt(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(nx,iy,is)*hx(nx-1,iy)+kinrgy(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy))))
        end do !}
        hlp=hlp/s
        wrk0(6,iy) = hlp + wrk0(0,iy)+wrk0(4,iy) 
        hlp_s=hlp_s+s
      end do !}
#ifdef B25_EIRENE
      js=ltns(nncut)+nlim
      do iy=0,ny-1 !{
      wrk0(6,iy)=wrk0(6,iy)+(wldnek(js,0)+wldnep(js,0))*
     &                              hlp_neu_peak/hlp_s
      enddo
#endif

! for double_null geometry     !iyv 04.02.11
! outer upper target

      if(nncut.eq.2) then !{

        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          s=gs(ix2,iy,0)
          hlp=0.0_R8
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(ix2,iy,0,is)*(
     .       (rpt(ix2,iy,is)*hx(ix2-1,iy)+
     .        rpt(ix2-1,iy,is)*hx(ix2,iy))/
     .       (hx(ix2-1,iy)+hx(ix2,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(ix2,iy,is)*hx(ix2-1,iy)+
     .        kinrgy(ix2-1,iy,is)*hx(ix2,iy))/
     .       (hx(ix2-1,iy)+hx(ix2,iy))))
          end do !}
          hlp=hlp/s
          wrk0(14,iy) = hlp + wrk0(8,iy)+wrk0(12,iy) 
          hlp_s=hlp_s+s
        end do !}
#ifdef B25_EIRENE
        js=ltns(3)+nlim
        do iy=0,ny-1 !{
          wrk0(14,iy) = wrk0(14,iy) +(wldnek(js,0)+wldnep(js,0))*
     &                                       hlp_neu_peak/hlp_s
        enddo
#endif
! inner upper target
        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          hlp=0.0_R8
          s=gs(ix1,iy,0)
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(ix1,iy,0,is)*(
     .       (rpt(ix1,iy,is)*hx(ix1-1,iy)+
     .        rpt(ix1-1,iy,is)*hx(ix1,iy))/
     .       (hx(ix1-1,iy)+hx(ix1,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(ix1,iy,is)*hx(ix1-1,iy)+
     .        kinrgy(ix1-1,iy,is)*hx(ix1,iy))/
     .       (hx(ix1-1,iy)+hx(ix1,iy))))
          end do !}
          hlp=hlp/s
          wrk0(13,iy) = hlp + wrk0(7,iy)+wrk0(11,iy) 
          hlp_s=hlp_s+s
        end do !}
#ifdef B25_EIRENE
        js=ltns(2)+nlim
        do iy=0,ny-1 !{
          wrk0(13,iy)=wrk0(13,iy)+(wldnek(js,0)+wldnep(js,0))*
     &                                 hlp_neu_peak/hlp_s
        enddo
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

      open (70,file='q_plates.dat') 


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
     &                                                             
	  do iy = ny, -1, -1
           write(70,"(5x,i2,16(1x,es10.3))") iy, 
     &             wrk0(6,iy), wrk0(14,iy), wrk0(5,iy), wrk0(13,iy),
     &             wrk0(0,iy), wrk0(8,iy), wrk0(4,iy), wrk0(12,iy),
     &             wrk0(0,iy)+wrk0(4,iy),  wrk0(8,iy)+wrk0(12,iy),
     &             wrk0(-1,iy), wrk0(7,iy), wrk0(3,iy), wrk0(11,iy),
     &             wrk0(-1,iy)+wrk0(3,iy),  wrk0(7,iy)+wrk0(11,iy)
                 
	  enddo
	  write(70,41) qe_oup_sum, qi_oup_sum, qe_out_sum, qi_out_sum, 
     &               qe_iup_sum, qi_iup_sum, qe_inn_sum, qi_inn_sum, 
     &               q_oup_sum, q_out_sum, q_iup_sum, q_inn_sum
                  
   41    format(/'Energy e-flux to outer upper plate     =',es11.3/
     &           'Energy i-flux to outer upper plate     =',es11.3/
     &           'Energy e-flux to outer lower plate     =',es11.3/
     &           'Energy i-flux to outer lower plate     =',es11.3//
   
     &           'Energy e-flux to inner low plate       =',es11.3/
     &           'Energy i-flux to inner low plate       =',es11.3/
     &           'Energy e-flux to inner upper plate     =',es11.3/
     &       	   'Energy i-flux to inner upper plate     =',es11.3// 

     &           'Total Energy flux to outer upper plate =',es11.3/
     &           'Total Energy flux to outer low plate   =',es11.3/
     &           'Total Energy flux to inner upper plate =',es11.3/ 
     &           'Total Energy flux to inner low plate   =',es11.3/ )

	else

	  write (70,"(5x,'iy',5x,'qe_out',5x,'qi_out',4x,'qei_out',
     &                    4x,'q_outtg',
     &                    5x,'qe_inn',5x,'qi_inn',4x,'qei_inn',
     &                    4x,'q_inntg')") 
	  do iy = ny, -1, -1
           write(70,"(5x,i2,8(1x,es10.3))") iy, 
     &             wrk0(0,iy), wrk0(4,iy), wrk0(0,iy) + wrk0(4,iy), 
     &             wrk0(6,iy), wrk0(-1,iy),wrk0(3,iy), 
     &             wrk0(-1,iy)+ wrk0(3,iy),wrk0(5,iy)
	  enddo
	  write(70,51) qe_out_sum, qi_out_sum, q_out_sum,
     &               qe_inn_sum, qi_inn_sum, q_inn_sum

   51    format(/'Energy e-flux to outer plate      =',es11.3/
     &           'Energy i-flux to outer plate      =',es11.3/
     &           'Total Energy flux to outer plate  =',es11.3/
     &           'Energy e-flux to inner plate      =',es11.3/
     &           'Energy i-flux to inner plate      =',es11.3//
     &           'Total Energy flux to inner plate  =',es11.3/)
      endif
	close(70)                                                    

	 ishift = 0

       open(70,file='table_low_inner.dat')
	  write(70,"(1x,9a14)") 
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
	  do iy = 0,ny
	     ixx=-1 !max(-1,ix-1) 
	     write(70,"(1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx+ishift,iy),
     &     0.5_R8*( crx(ixx+ishift,iy,0) + crx(ixx+ishift,iy,2) ),
     &     te(ixx+ishift,iy)/ev, ti(ixx+ishift,iy)/ev,ne(ixx+ishift,iy),

     &     qe * hy1(ixx+ishift,iy)/hy(ixx+ishift,iy) *
     &     bb(ixx+ishift,iy,0)/bb(ixx+ishift,iy,3) *
     &     na(ixx+ishift,iy,ismain) *
     &     sqrt((te(ixx+ishift,iy)+ti(ixx+ishift,iy))/(am(ismain)*mp)),
     &     wrk0(-1,iy)/1.0e6_R8,wrk0(3,iy)/1.0e6_R8,wrk0(5,iy)/1.0e6_R8        
	  enddo
	 close(70)

       open(70,file='table_low_outer.dat')
	  write(70,"(1x,9a14)") 
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
	  do iy = 0,ny
	     ixx=nx !min(nx,ix+1) 
	     write(70,"(1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx-ishift,iy), 
     &     0.5_R8*( crx(ixx-ishift,iy,0) + crx(ixx-ishift,iy,2) ),
     &     te(ixx-ishift,iy)/ev, ti(ixx-ishift,iy)/ev,ne(ixx-ishift,iy),
     &     qe * hy1(ixx-ishift,iy)/hy(ixx-ishift,iy) *
     &     bb(ixx+ishift,iy,0)/bb(ixx-ishift,iy,3) *
     &     na(ixx-ishift,iy,ismain) *
     &     sqrt((te(ixx-ishift,iy)+ti(ixx-ishift,iy))/(am(ismain)*mp)),
     &     wrk0(0,iy)/1.0e6_R8, wrk0(4,iy)/1.0e6_R8, wrk0(6,iy)/1.0e6_R8
	  enddo
	 close(70)                                       !srv 19.09.11 }

       if(nncut.eq.2) then 

        open(70,file='table_upper_inner.dat')
	   write(70,"(1x,9a14)") 
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
	   do iy = 0,ny
	     ixx=ix1
	     write(70,"( 1x,9es14.6)")                              !srv 03.04.13
     &     crdy_c_cell(ixx+ishift,iy),
     &     0.5_R8*( crx(ixx+ishift,iy,0) + crx(ixx+ishift,iy,2) ),
     &     te(ixx+ishift,iy)/ev, ti(ixx+ishift,iy)/ev,ne(ixx+ishift,iy),

     &     qe * hy1(ixx+ishift,iy)/hy(ixx+ishift,iy) *
     &     bb(ixx+ishift,iy,0)/bb(ixx+ishift,iy,3) *
     &     na(ixx+ishift,iy,ismain) *
     &     sqrt((te(ixx+ishift,iy)+ti(ixx+ishift,iy))/(am(ismain)*mp)),
     &   wrk0(7,iy)/1.0e6_R8, wrk0(11,iy)/1.0e6_R8, wrk0(13,iy)/1.0e6_R8        
	   enddo
	  close(70)

        open(70,file='table_upper_outer.dat')
	   write(70,"(1x,9a14)") 
     &       'l(m)', 'Rc(m)', 'Te(eV)', 'Ti(eV)',
     &       'ne(1/m^3)', 'jsat(A/m^2)', 'qe(MW/m^2)','qi(MW/m^2)',
     &       'qtot(MW/m^2)'
	   do iy = 0,ny
	     ixx=ix2 
	     write(70,"( 1x,9es14.6)")                               !srv 03.04.13
     &     crdy_c_cell(ixx-ishift,iy), 
     &     0.5_R8*( crx(ixx-ishift,iy,0) + crx(ixx-ishift,iy,2) ),
     &     te(ixx-ishift,iy)/ev, ti(ixx-ishift,iy)/ev,ne(ixx-ishift,iy),
     &     qe * hy1(ixx-ishift,iy)/hy(ixx-ishift,iy) *
     &     bb(ixx+ishift,iy,0)/bb(ixx-ishift,iy,3) *
     &     na(ixx-ishift,iy,ismain) *
     &     sqrt((te(ixx-ishift,iy)+ti(ixx-ishift,iy))/(am(ismain)*mp)),
     &   wrk0(8,iy)/1.0e6_R8, wrk0(12,iy)/1.0e6_R8, wrk0(14,iy)/1.0e6_R8
	   enddo
	  close(70)                                       !srv 19.09.11 }
       endif
      endif      !  iout.eq.1 or. iout_iout_ASTRA.eq.1    !srv 19.04.13
!______________________________________________________________________________ !iyv 15.10.09 }

	if(iout.eq.2) then                                !srv 06.05.09 {
       mintopcut = minval(topcut(1:nncut))
*   ..write crdy_c_cell, crdy_c_cell_abs [mm], ne, te, ti, js at inner
       open(70,file='inner.dat')
	  do iy = -1,ny                                        !srv 03.04.13 {
         iw = 0      
	   do ix = -1,nx
	    if(leftix(ix,iy).eq.-2) then
	     ixx=max(-1,ix-1) 
           wdim(iw+1 )= crdy_c_cell(ixx+ishift,iy)
           wdim(iw+2)= ne(ixx+ishift,iy)
           wdim(iw+3)= te(ixx+ishift,iy)/ev
           wdim(iw+4)= ti(ixx+ishift,iy)/ev
           iw = iw + 4
!	     write(70,"(\1x,12es14.6)")
!     &     crdy_c_cell(ixx+ishift,iy),
!!     &     (crdy_c_cell(ixx+ishift,iy)+cry(ixx+ishift,mintopcut,0))*1e3,
!     &     ne(ixx+ishift,iy),te(ixx+ishift,iy)/ev, ti(ixx+ishift,iy)/ev
!!     &    ,ne(ixx+ishift,iy),te(ixx+ishift,iy), ti(ixx+ishift,iy),
!!     &    ,qe*ne(ixx+ishift,iy)*bb(ixx+ishift,iy,0)/bb(ixx+ishift,iy,3)*
!!     &     sqrt((te(ixx+ishift,iy)+ti(ixx+ishift,iy))/(mp)) 
	    endif
	   enddo
	   write(70,"(2(1x,4es14.6))") (wdim(ix),ix=1,iw)
	  enddo                                              !srv 03.04.13 }
	 close(70)

       open(70,file='midpln.dat')
	 ix=nmdpl
	  do iy=-1,ny 
	   write(70,"(1x,5es14.6)") 
     *    crdy_c_cell(ix,iy), 
!     &    crx(ix,iy,0)*1e3, 
     &    ne(ix,iy), te(ix,iy)/ev, ti(ix,iy)/ev 
	  enddo
	 close(70)

       open(70,file='outer.dat')
	  do iy = -1,ny                                        !srv 03.04.13 {
         iw = 0
	   do ix = -1,nx
	    if(rightix(ix,iy).eq.nx+1) then
	     ixx=min(nx,ix+1) 
           wdim(iw+1 )= crdy_c_cell(ixx-ishift,iy)
           wdim(iw+2)= ne(ixx-ishift,iy)
           wdim(iw+3)= te(ixx-ishift,iy)/ev
           wdim(iw+4)= ti(ixx-ishift,iy)/ev
           iw = iw + 4
!	     write(70,"(\1x,12es14.6)") 
!     &     crdy_c_cell(ixx-ishift,iy), 
!!     &     (crdy_c_cell(ixx-ishift,iy)+crx(ixx-ishift,mintopcut,0))*1e3,
!     &     ne(ixx-ishift,iy),te(ixx-ishift,iy)/ev, ti(ixx-ishift,iy)/ev 
!!     &    ,ne(ixx-ishift,iy),te(ixx-ishift,iy), ti(ixx-ishift,iy) 
!!     &    ,qe*ne(ixx-ishift,iy)*bb(ixx-ishift,iy,0)/bb(ixx-ishift,iy,3)*
!!     &     sqrt((te(ixx-ishift,iy)+ti(ixx-ishift,iy))/(mp)) 
	    endif
	   enddo
	   write(70,"(2(1x,4es14.6))") (wdim(ix),ix=1,iw)
	  enddo                                              !srv 03.04.13 }
	 close(70)

       open(70,file='midplnUi.dat')
	 ix=nmdpl
	 do iy=-1,ny
	  write(70,"(1x,2es14.6)") 
     *   crdy_c_cell(ix,iy), ua(ix,iy,1)
     	 enddo
	 close(70)

       open(70,file='midplnEy.dat')
	 ix=nmdpl
	 do iy=0,ny
	  write(70,"(1x,2es14.6)") 
     *   crdy_surf(ix,iy), 
     *     -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))/
     *     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*2.0e0_R8
     	 enddo
	 close(70)
*   ..calculate total energy flux
	 iy = 0
       power=0.0_R8
       do ix = -1, nx
	  if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
	    power = power + (fhe(ix,iy,1)+fhi(ix,iy,1))
	  endif
	 enddo
	 open(70,file='power.dat')
	  write(70,"('Power =',es14.6)") power
	 close(70)

       open(70,file='separatrix_location_on_plate.dat')
	  ix = -1
	  write(70,"(1x,' Separatrix Locationon on Inner Lower =',
     &             12es14.6)") cry(ix+ishift,mintopcut,0)
	  ix = nx
	  write(70,"(1x,' Separatrix Locationon on Outer Lower =',
     &             12es14.6)") crx(ix-ishift,mintopcut,0)
	  do ix = -1,nx
	   if(leftix(ix,mintopcut).eq.-2 .and. ix.gt.-1) then
	    write(70,"(1x,' Separatrix Locationon on Outer Upper =',
     &             12es14.6)") crx(ix+ishift,mintopcut,0)
	    endif
	   if(rightix(ix,mintopcut).eq.nx+1 .and. ix.lt.nx) then
	    write(70,"(1x,' Separatrix Locationon on Inner Upper =',
     &             12es14.6)") cry(ix-ishift,mintopcut,0)
	    endif
	  enddo
	  ix = nmdpl
	  write(70,"(1x,' Separatrix Locationon on Outer Mdpln =',
     &             12es14.6)") crx(ix-ishift,mintopcut,0)
	 close(70)

*   ..write crdy_c_cell, crdy_c_cell_abs [m], na, ua at inner
       open(70,file='inner_na_ua.dat')                         !srv 03.04.13 {
	  write(hlp_frm,"(a6,i4,a8)") '(2(1x,',2*ns+1,'es14.6))'
	  do iy = -1,ny                               
         iw = 0      
	   do ix = -1,nx
	    if(leftix(ix,iy).eq.-2) then
	     ixx=max(-1,ix-1) 
           iw = iw +1
           wdim(iw )= crdy_c_cell(ixx+ishift,iy)
           do is=0,ns-1
            iw = iw + 1
            wdim(iw)= na(ixx+ishift,iy,is)
	     enddo
           do is=0,ns-1
            iw = iw + 1
            wdim(iw)= ua(ixx+ishift,iy,is)
           enddo
!	     write(70,"(\1x,12es14.6)")
!     &     crdy_c_cell(ixx+ishift,iy),
!     &     (na(ixx+ishift,iy,is),is=0,ns-1), 
!     &     (ua(ixx+ishift,iy,is),is=0,ns-1)
	    endif
	   enddo
	   write(70,hlp_frm) (wdim(ix),ix=1,iw)
!	   write(70,"(<2*ns+1>es14.6)")
	  enddo  
	 close(70)
*   ..write crdy_c_cell, crdy_c_cell_abs [m], na, ua at outer
       open(70,file='outer_na_ua.dat')
	  do iy = -1,ny
         iw = 0      
	   do ix = -1,nx
	    if(rightix(ix,iy).eq.nx+1) then
	     ixx=min(nx,ix+1) 
           iw = iw +1
           wdim(iw )= crdy_c_cell(ixx-ishift,iy)
           do is=0,ns-1
            iw = iw + 1
            wdim(iw)= na(ixx-ishift,iy,is)
	     enddo
           do is=0,ns-1
            iw = iw + 1
            wdim(iw)= ua(ixx-ishift,iy,is)
           enddo
!	     write(70,"(\1x,12es14.6)") 
!     &     crdy_c_cell(ixx-ishift,iy), 
!     &     (na(ixx-ishift,iy,is),is=0,ns-1),
!     &     (ua(ixx-ishift,iy,is),is=0,ns-1) 
	    endif
	   enddo
	   write(70,hlp_frm) (wdim(ix),ix=1,iw)
!	   write(70,"(<2*ns+1>es14.6)")
	  enddo
	 close(70)                                                !srv 03.04.13 }

      endif              !srv   endif for iout=2        !srv 06.05.09 }

	if(iout.eq.3) then


	endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wdat

      end subroutine b2wdat
