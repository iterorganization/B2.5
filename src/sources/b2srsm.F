*-----------------------------------------------------------------------
*.specification
      subroutine b2srsm (nx, ny, ns, dtim, vol, na, ua, te, ti, ne, ni,
     &  snasm, smosm, shesm, shism, main_call)
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  dtim
      real (kind=R8) ::
     *  vol(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1), 
     *  ua(-1:nx,-1:ny,0:ns-1), ni(-1:nx,-1:ny,0:1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny)
      logical main_call 
*   ..input/output arguments
      real (kind=R8) ::
     *  snasm(-1:nx,-1:ny,0:1,0:ns-1), smosm(-1:nx,-1:ny,0:3,0:ns-1), 
     &  shesm(-1:nx,-1:ny,0:3), shism(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRSM smoothes expressions the sources strengths of
*     density, momentum, electron heat, atom heat, and charge.
*     The idea is to limit the source strength by (value)/(dt*t._f._r.)
*     - and to do this smoothly!
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, ix, iy, n2, ncall, diagno, enable, numerics_namelist
      integer elm_count, k
	integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0                                
	character chns*2, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *  timescale, invtimescale, timefactor, hlp
      logical catch_up
      save elm_count
*   ..procedures
      real (kind=R8) :: roxa
      roxa(ix,iy,is) = am(is)*mp*na(ix,iy,is)
      external subini, subend, xertst, ipgeti, b2xvsg
      intrinsic abs, min, mod, int
*   ..initialisation
      save ncall, diagno, enable, numerics_namelist
      data ncall/0/, diagno/0/, enable/0/, numerics_namelist/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srsm')
*   ..compute n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test dtim
      call xertst (0.lt.dtim, 'faulty argument dtim')
*   ..test vol
      call b2xvsg (n2, vol, 1, 'vol', '.gt.')
*   ..test sign of na, ni, ne, te, ti
      do is = 0, ns-1
        call b2xvsg (n2, na(-1,-1,is), 1, 'na', '.gt.')
      enddo
      do is = 0, 1
        call b2xvsg (n2, ni(-1,-1,is), 1, 'ni', '.gt.')
      enddo
      call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      call b2xvsg (n2, te, 1, 'te', '.gt.')
      call b2xvsg (n2, ti, 1, 'ti', '.gt.')
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
        if(enable.eq.0) then
          write(*,*) 'b2srsm disabled'
        endif
        call ipgeti ('b2srdt_numerics_namelist', numerics_namelist)
        call alloc_b2mod_numerics(ns)
        if(numerics_namelist.ge.1) then
          call read_b2mod_numerics_namelist
          call write_b2mod_numerics_namelist
        endif
        if(eirene_mod.le.0) eirene_mod=1
        if(numerics_time_mod.gt.0.0_R8) then
          elm_count=tim/numerics_time_mod
        else
          elm_count=0
        endif
        catch_up=numerics_time_mod.gt.0.0_R8
        do while(catch_up)
          if(mod(tim,numerics_time_mod).gt.numerics_time_switch)
     1     then
            write(*,*) 'b2srsm: catching up ...' 
            call read_b2mod_numerics_namelist
            call write_b2mod_numerics_namelist
            catch_up=numerics_time_mod.gt.0.0_R8
          else
            catch_up=.false.
          endif
        enddo
      else
        if(numerics_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,numerics_time_mod),
     1     numerics_time_switch,int(tim/numerics_time_mod),elm_count
          if(mod(tim,numerics_time_mod).gt.numerics_time_switch.or.
     1     int(tim/numerics_time_mod).gt.elm_count) then
            elm_count=tim/numerics_time_mod
            call read_b2mod_numerics_namelist
            call write_b2mod_numerics_namelist
          endif
        endif
      endif

      if(enable.ne.0) then
c*** calculate the integral sources for nregion regions:
c*** before correction (1), after correction (2)
        timefactor = 1.0_R8

*   ..smooth particle and momentum sources
        do is = 0, ns-1
          do ix = -1, nx
            do iy = -1, ny
              if (resignore(ix,iy,1).ne.0) then
c check particles
                invtimescale =
     1           abs(snasm(ix,iy,0,is)+snasm(ix,iy,1,is)*na(ix,iy,is))
     2           /(na(ix,iy,is)*vol(ix,iy))  
                hlp = invtimescale*
     1           time_factor_required*dtim*dtco(region(ix,iy,0),is)
                if (hlp.gt.1.0_R8) then
                  timescale  = 1.0_R8/invtimescale
                  timefactor = min(timefactor,1.0_R8/hlp)
                  if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                    write(*,*) 'WARNING: particle sources very strong @'
                    write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,               timescale,dtim*dtco(region(ix,iy,0),is)
                    write(*,'(a,3i4)') 'ix,iy,is = ', ix,iy,is
                    write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1               te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                  endif
                  snasm(ix,iy,0,is) = snasm(ix,iy,0,is)/hlp
                  snasm(ix,iy,1,is) = snasm(ix,iy,1,is)/hlp
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
c check momentum
                invtimescale = abs(smosm(ix,iy,0,is) +
     1           smosm(ix,iy,1,is)*ua(ix,iy,is) +
     &           smosm(ix,iy,2,is)*roxa(ix,iy,is) + 
     &           smosm(ix,iy,3,is)*roxa(ix,iy,is)*ua(ix,iy,is))
     2           /(sqrt(ti(ix,iy)/
     3           (am(is)*mp))*roxa(ix,iy,is)*vol(ix,iy))  
                hlp = invtimescale*
     1           time_factor_required*dtim*dtmo(region(ix,iy,0),is)
                if (hlp.gt.1.0_R8) then
                  timescale  = 1.0_R8/invtimescale
                  timefactor = min(timefactor,1.0_R8/hlp)
                  if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                    write(*,*) 'WARNING: momentum sources very strong @'
                    write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,               timescale,dtim*dtmo(region(ix,iy,0),is)
                    write(*,'(a,6i4)') 'ix,iy,is = ', ix,iy,is
                    write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1               te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                  endif
                  smosm(ix,iy,0,is) = smosm(ix,iy,0,is)/hlp
                  smosm(ix,iy,1,is) = smosm(ix,iy,1,is)/hlp
                  smosm(ix,iy,2,is) = smosm(ix,iy,2,is)/hlp
                  smosm(ix,iy,3,is) = smosm(ix,iy,3,is)/hlp
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
*   ..smooth electron heat, ion heat and charge sources
        do ix = -1, nx
          do iy = -1, ny
            if (resignore(ix,iy,1).ne.0) then
c check ion energy
              invtimescale =
     1         abs(shism(ix,iy,0) + shism(ix,iy,1)*ti(ix,iy) + 
     &         shism(ix,iy,2)*ni(ix,iy,0) +
     &         shism(ix,iy,3)*ti(ix,iy)*ni(ix,iy,0))
     2         /(ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy))
              hlp = invtimescale*time_factor_required*
     .         dtim*dtei(region(ix,iy,0))
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                  write(*,*) 'WARNING: ion heat sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,dtim*dtei(region(ix,iy,0))
                  write(*,'(a,6i4)') 'ix,iy,is = ', ix,iy,is
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                endif
                shism(ix,iy,0) = shism(ix,iy,0)/hlp
                shism(ix,iy,1) = shism(ix,iy,1)/hlp
                shism(ix,iy,2) = shism(ix,iy,2)/hlp
                shism(ix,iy,3) = shism(ix,iy,3)/hlp
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
c check electron energy
              invtimescale =
     1         abs(shesm(ix,iy,0) + shesm(ix,iy,1)*te(ix,iy) + 
     &         shesm(ix,iy,2)*ne(ix,iy) +
     &         shesm(ix,iy,3)*te(ix,iy)*ne(ix,iy))
     2         /(te(ix,iy)*ne(ix,iy)*vol(ix,iy))
              hlp = invtimescale*time_factor_required*
     .         dtim*dtee(region(ix,iy,0))
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                  write(*,*) 
     ,             'WARNING: electron heat sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,dtim*dtee(region(ix,iy,0))
                  write(*,'(a,6i4)') 'ix,iy,is = ', ix,iy,is
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy)
                endif
                shesm(ix,iy,0) = shesm(ix,iy,0)/hlp
                shesm(ix,iy,1) = shesm(ix,iy,1)/hlp
                shesm(ix,iy,2) = shesm(ix,iy,2)/hlp
                shesm(ix,iy,3) = shesm(ix,iy,3)/hlp
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
            endif
          enddo    
        enddo
      endif

      if (iout_b2npmo.eq.1) then                                     !srv 11.09.09 {
	 do is = 0, ns-1          
	  write (chns,'(i2.2)') is
!	  call my_out(70,nx,ny,smosm(:,:,0,is)+                        !srv 07.10.09
!     &        smosm(:,:,1,is)*ua(:,:,is)+
!     &        smosm(:,:,2,is)*na(:,:,is)*mp*am(is)+
!     &        smosm(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is),
!     &        'b2srsm_smosm'//chns)
	 enddo
      endif
      if (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1) then
	 do is = 0, ns-1
	  write (chns,'(i2.2)') is
!        call my_out(70,nx,ny,snasm(:,:,0,is)+                       !srv 07.10.09
!     &             snasm(:,:,1,is)*na(:,:,is),'b2srsm_snasm'//chns)
	 enddo
	endif
	if (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1) then
!	 call my_out(70,nx,ny,shesm(:,:,0)+shesm(:,:,1)*te+            !srv 07.10.09
!     &   shesm(:,:,2)*ne + shesm(:,:,3)*te*ne,'b2srsm_shesm')
!	 call my_out(70,nx,ny,shism(:,:,0)+shism(:,:,1)*ti+
!     &   shism(:,:,2)*ni(:,:,0) + 
!     &   shism(:,:,3)*ti*ni(:,:,0),'b2srsm_shism')
	endif            

	if (iout.eq.1) then
       do is = 0, ns-1
	  write (chns,'(i2.2)') is
	  do k=0,3
	   write (chk,'(i1)') k
	   call my_out(70,nx,ny,smosm(-1,-1,k,is),
     &                       'b2srsm_smosm'//chk//chns)
	  enddo
	  do k=0,1
	   write (chk,'(i1)') k
	   call my_out(70,nx,ny,snasm(-1,-1,k,is),
     &                       'b2srsm_snasm'//chk//chns)
	  enddo
	 enddo
	 do k=0,3
	  call my_out(70,nx,ny,shesm(-1,-1,k),'b2srsm_shesm'//chk)
	  call my_out(70,nx,ny,shism(-1,-1,k),'b2srsm_shism'//chk)
	 enddo
	endif                                                          !srv 11.09.09 }

* ..return
      if (ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srsm

      end
