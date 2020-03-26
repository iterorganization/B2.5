*-----------------------------------------------------------------------
*.specification
      subroutine b2srsm (nCv, ns, dtim, geo, mpg, na, ua, te, ti, ne,
     &  ni, sr, main_call)
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     1 , only : b2srsm_sna0to1, b2srsm_smo0to3, b2srsm_she0to3,
!WG_TODO     2          b2srsm_shi0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  dtim
      real (kind=R8) ::
     *  na(nCv,0:ns-1), 
     *  ua(nCv,0:ns-1), ni(nCv,0:1),
     &  te(nCv), ti(nCv), ne(nCv)
      logical main_call 
*   ..input/output arguments
      type (B2Source), intent (inout) :: sr
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRSM smoothes expressions the sources strengths of
*     density, parallel momentum, electron heat, atom heat, and charge.
*     The idea is to limit the source strength by (value)/(dt*t._f._r.)
*     - and to do this smoothly!
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, iCv, ncall, diagno, enable
      integer k
      integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0, tn_style = 0                                
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *  timescale, invtimescale, timefactor, hlp
      real (kind=R8) :: !djm Jan2017
     &  snatmp(nCv,0:1,0:ns-1), smotmp(nCv,0:3,0:ns-1), 
     &  shetmp(nCv,0:3), shitmp(nCv,0:3)
*   ..procedures
      real (kind=R8) :: roxa
      roxa(iCv,is) = am(is)*mp*na(iCv,is)
      external subini, subend, xertst, ipgeti, b2xvsg
      intrinsic abs, min, mod, int
*   ..initialisation
      save ncall, diagno, enable
      data ncall/0/, diagno/0/, enable/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srsm')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test dtim
      call xertst (0.lt.dtim, 'faulty argument dtim')
*   ..test sign of na, ni, ne, te, ti
      do is = 0, ns-1
        call b2xvsg (nCv, na(1,is), 1, 'na', '.gt.')
      enddo
      do is = 0, 1
        call b2xvsg (nCv, ni(1,is), 1, 'ni', '.gt.')
      enddo
      call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      call b2xvsg (nCv, te, 1, 'te', '.gt.')
      call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
*   ..obtain input
      if (ncall.eq.0) then
        call ipgeti ('b2srsm_iout', iout)                         !srv 11.09.09 {
        call ipgeti ('b2npc9_iout', iout_b2npc9) 
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)                  !srv 11.09.09 }
        call ipgeti ('b2srsm_diagno', diagno)
        call ipgeti ('b2srsm_enable', enable)
        call ipgeti ('b2mn_tn_style', tn_style)
        if(enable.eq.0) then
          write(*,*) 'b2srsm disabled'
        endif
        if(eirene_mod.le.0) eirene_mod=1
      endif

cdjm Jan2017 Keep sources at the start to calculate total change for balance routines
      snatmp = sr%sna
      smotmp = sr%smo
      shetmp = sr%she
      shitmp = sr%shi

      if(enable.ne.0) then
        call xerrab ('b2srsm not adapted for WG')
#ifdef WG_TODO
c*** calculate the integral sources for all regions:
c*** before correction (1), after correction (2)
        timefactor = 1.0_R8

*   ..smooth particle and parallel momentum sources
        do is = 0, ns-1

          do ix = -1, nx
            do iy = -1, ny
              if (resignore(ix,iy,1).ne.0) then
c check particles
                invtimescale =
     1           abs(sr%sna(ix,iy,0,is)+sr%sna(ix,iy,1,is)*na(ix,iy,is))
     2           /(na(ix,iy,is)*vol(ix,iy))  
                hlp = invtimescale*time_factor_required*
     1           dtim*dtco(is,region(ix,iy,0))*time_factor(ix,iy)
                if (hlp.gt.1.0_R8) then
                  timescale  = 1.0_R8/invtimescale
                  timefactor = min(timefactor,1.0_R8/hlp)
                  if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                    write(*,*) 'WARNING: particle sources very strong @'
                    write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,               timescale,
     2               dtim*dtco(is,region(ix,iy,0))*time_factor(ix,iy)
                    write(*,'(a,2i4,i3)') 'ix,iy,is = ', ix,iy,is
                    write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1               te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                  endif
                  sr%sna(ix,iy,0,is) = sr%sna(ix,iy,0,is)/hlp
                  sr%sna(ix,iy,1,is) = sr%sna(ix,iy,1,is)/hlp
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_ni(region(ix,iy,0),2,b2espcr(is)) = 
     &                src_ni(region(ix,iy,0),2,b2espcr(is)) +
     &                cor_ni(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                    src_ni(0,2,b2espcr(is)) = 
     &                src_ni(0,2,b2espcr(is)) +
     &                cor_ni(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  endif
                endif
              endif
            enddo    
          enddo
          do ix = -1, nx
            do iy = -1, ny
              if (resignore(ix,iy,2).ne.0) then
c check parallel momentum
                invtimescale = abs(sr%smo(ix,iy,0,is) +
     1           sr%smo(ix,iy,1,is)*ua(ix,iy,is) +
     &           sr%smo(ix,iy,2,is)*roxa(ix,iy,is) + 
     &           sr%smo(ix,iy,3,is)*roxa(ix,iy,is)*ua(ix,iy,is))
     2           /(sqrt(ti(ix,iy)/
     3           (am(is)*mp))*roxa(ix,iy,is)*vol(ix,iy))  
                hlp = invtimescale*time_factor_required*
     1           dtim*dtmo(is,region(ix,iy,0))*time_factor(ix,iy)
                if (hlp.gt.1.0_R8) then
                  timescale  = 1.0_R8/invtimescale
                  timefactor = min(timefactor,1.0_R8/hlp)
                  if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                    write(*,*) 
     1               'WARNING: parallel momentum sources very strong @'
                    write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,               timescale,
     2               dtim*dtmo(is,region(ix,iy,0))*time_factor(ix,iy)
                    write(*,'(a,2i4,i3)') 'ix,iy,is = ', ix,iy,is
                    write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1               te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                  endif
                  sr%smo(ix,iy,0,is) = sr%smo(ix,iy,0,is)/hlp
                  sr%smo(ix,iy,1,is) = sr%smo(ix,iy,1,is)/hlp
                  sr%smo(ix,iy,2,is) = sr%smo(ix,iy,2,is)/hlp
                  sr%smo(ix,iy,3,is) = sr%smo(ix,iy,3,is)/hlp
                  if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                    src_mo(region(ix,iy,0),2,b2espcr(is)) = 
     &                src_mo(region(ix,iy,0),2,b2espcr(is)) +
     &                cor_mo(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                    src_mo(0,2,b2espcr(is)) = 
     &                src_mo(0,2,b2espcr(is)) +
     &                cor_mo(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                    src_ms(region(ix,iy,0),2) = 
     &                src_ms(region(ix,iy,0),2) +
     &                cor_mo(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                    src_ms(0,2) = src_ms(0,2) +
     &                cor_mo(ix,iy,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  endif
                endif
              endif
            enddo    
          enddo
        enddo
*   ..smooth electron heat, ion heat and electron particle sources
        do ix = -1, nx
          do iy = -1, ny
            if (resignore(ix,iy,1).ne.0) then
c check ion energy
              invtimescale =
     1         abs(sr%shi(ix,iy,0) + sr%shi(ix,iy,1)*ti(ix,iy) + 
     &         sr%shi(ix,iy,2)*ni(ix,iy,0) +
     &         sr%shi(ix,iy,3)*ti(ix,iy)*ni(ix,iy,0))
     2         /(ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy))
              hlp = invtimescale*time_factor_required*
     .         dtim*dtei(region(ix,iy,0))*time_factor(ix,iy)
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                  write(*,*) 'WARNING: ion heat sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,
     2             dtim*dtei(region(ix,iy,0))*time_factor(ix,iy)
                  write(*,'(a,2i4)') 'ix,iy = ', ix,iy
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                endif
                sr%shi(ix,iy,0) = sr%shi(ix,iy,0)/hlp
                sr%shi(ix,iy,1) = sr%shi(ix,iy,1)/hlp
                sr%shi(ix,iy,2) = sr%shi(ix,iy,2)/hlp
                sr%shi(ix,iy,3) = sr%shi(ix,iy,3)/hlp
                if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                  src_ei(region(ix,iy,0),2) = 
     &              src_ei(region(ix,iy,0),2) +
     &              cor_ei(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_ei(0,2) = src_ei(0,2) +
     &              cor_ei(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_es(region(ix,iy,0),2) = 
     &              src_es(region(ix,iy,0),2) +
     &              cor_ei(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_es(0,2) = src_es(0,2) +
     &              cor_ei(ix,iy)*(1.0_R8/hlp-1.0_R8)
                endif
              endif
c check neutral energy
          if (tn_style.ne.0) then
            call xerrab('b2srsm not developed for tn_style.ne.0')
          endif
c check electron energy
              invtimescale =
     1         abs(sr%she(ix,iy,0) + sr%she(ix,iy,1)*te(ix,iy) + 
     &         sr%she(ix,iy,2)*ne(ix,iy) +
     &         sr%she(ix,iy,3)*te(ix,iy)*ne(ix,iy))
     2         /(te(ix,iy)*ne(ix,iy)*vol(ix,iy))
              hlp = invtimescale*time_factor_required*
     .         dtim*dtee(region(ix,iy,0))*time_factor(ix,iy)
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                  write(*,*) 
     ,             'WARNING: electron heat sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,
     2             dtim*dtee(region(ix,iy,0))*time_factor(ix,iy)
                  write(*,'(a,2i4)') 'ix,iy = ', ix,iy
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                endif
                sr%she(ix,iy,0) = sr%she(ix,iy,0)/hlp
                sr%she(ix,iy,1) = sr%she(ix,iy,1)/hlp
                sr%she(ix,iy,2) = sr%she(ix,iy,2)/hlp
                sr%she(ix,iy,3) = sr%she(ix,iy,3)/hlp
                if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                  src_ee(region(ix,iy,0),2) = 
     &              src_ee(region(ix,iy,0),2) +
     &              cor_ee(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_ee(0,2) = src_ee(0,2) +
     &              cor_ee(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_es(region(ix,iy,0),2) = 
     &              src_es(region(ix,iy,0),2) +
     &              cor_ee(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_es(0,2) = src_es(0,2) +
     &              cor_ee(ix,iy)*(1.0_R8/hlp-1.0_R8)
                endif
              endif
c check electron particle number
              invtimescale =
     1         abs(sr%sne(ix,iy,0) + sr%sne(ix,iy,1)*ne(ix,iy))
     2         /(ne(ix,iy)*vol(ix,iy))
              hlp = invtimescale*time_factor_required*
     1         dtim*time_factor(ix,iy)
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                  write(*,*) 
     ,             'WARNING: electron particle sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,dtim*time_factor(ix,iy)
                  write(*,'(a,2i4)') 'ix,iy = ', ix,iy
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                endif
                sr%sne(ix,iy,0) = sr%sne(ix,iy,0)/hlp
                sr%sne(ix,iy,1) = sr%sne(ix,iy,1)/hlp
                if (main_call.and.mod(ncall,eirene_mod).eq.0) then
                  src_ne(region(ix,iy,0),2) = 
     &              src_ne(region(ix,iy,0),2) +
     &              cor_ne(ix,iy)*(1.0_R8/hlp-1.0_R8)
                  src_ne(0,2) = src_ne(0,2) +
     &              cor_ne(ix,iy)*(1.0_R8/hlp-1.0_R8)
                endif
              endif
            endif
          enddo    
        enddo
#endif
      endif

cdjm Store the total change in linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2srsm_sna0to1 = sr%sna-snatmp
!WG_TODO        b2srsm_smo0to3 = sr%smo-smotmp
!WG_TODO        b2srsm_she0to3 = sr%she-shetmp
!WG_TODO        b2srsm_shi0to3 = sr%shi-shitmp
!WG_TODO      endif

      if (iout_b2npmo.eq.1) then                                     !srv 11.09.09 {
	 do is = 0, ns-1          
	  write (chns,'(i3.3)') is
!	  call my_out(70,nCv,sr%smo(:,:,0,is)+                        !srv 07.10.09
!     &        sr%smo(:,:,1,is)*ua(:,:,is)+
!     &        sr%smo(:,:,2,is)*na(:,:,is)*mp*am(is)+
!     &        sr%smo(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
!     &        'b2srsm_smosm'//chns)
	 enddo
      endif
      if (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1) then
        do is = 0, ns-1
	  write (chns,'(i3.3)') is
!        call my_out(70,nCv,sr%sna(:,:,0,is)+                       !srv 07.10.09
!     &             sr%sna(:,:,1,is)*na(:,:,is),'b2srsm_snasm'//chns)
        enddo
      endif
      if (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1) then
!	 call my_out(70,nCv,sr%she(:,:,0)+sr%she(:,:,1)*te+            !srv 07.10.09
!     &   sr%she(:,:,2)*ne + sr%she(:,:,3)*te*ne,'b2srsm_shesm')
!	 call my_out(70,nCv,sr%shi(:,:,0)+sr%shi(:,:,1)*ti+
!     &   sr%shi(:,:,2)*ni(:,:,0) + 
!     &   sr%shi(:,:,3)*ti*ni(:,:,0),'b2srsm_shism')
      endif            

      if (iout.eq.1) then
        do is = 0, ns-1
	  write (chns,'(i3.3)') is
	  do k=0,3
	   write (chk,'(i1)') k
	   call my_out_us(70,nCv,0,sr%smo(1,k,is),
     &                       'b2srsm_smosm'//chk//chns)
	  enddo
	  do k=0,1
	   write (chk,'(i1)') k
	   call my_out_us(70,nCv,0,sr%sna(1,k,is),
     &                       'b2srsm_snasm'//chk//chns)
	  enddo
	enddo
        do k=0,3
	  call my_out_us(70,nCv,0,sr%she(1,k),'b2srsm_shesm'//chk)
	  call my_out_us(70,nCv,0,sr%shi(1,k),'b2srsm_shism'//chk)
	enddo
	do k=0,1
	  call my_out_us(70,nCv,0,sr%sne(1,k),'b2srsm_snesm'//chk)
	enddo
      endif                                                          !srv 11.09.09 }

* ..return
      if (ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srsm

      end
