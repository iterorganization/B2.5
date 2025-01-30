*-----------------------------------------------------------------------
*.specification
      subroutine b2srsm (nCv, ns, dtim, switch, geo, mpg, na, ua,
     &   te, ti, ne, ni, nn, sr, main_call)
      use b2mod_types
      use b2mod_numerics_namelist
      use b2mod_neutrals_namelist
      use b2mod_constants
      use b2mod_diag
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     1 , only : b2srsm_sna0to1, b2srsm_smo0to3, b2srsm_she0to3,
     2          b2srsm_shi0to3, balance_netcdf
#endif
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2srsm
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  dtim
      real (kind=R8) ::
     *  na(nCv,0:ns-1),
     *  ua(nCv,0:ns-1), ni(nCv,0:1),
     &  te(nCv), ti(nCv), ne(nCv), nn(nCv)
      logical main_call
*   ..input/output arguments
      type (B2Source), intent (inout) :: sr
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
      integer is, iCv
      integer k
      character chns*3, chk*1                                        !srv 11.09.09 }
      integer diagno
      real (kind=R8) ::
     *  timescale, invtimescale, timefactor, hlp
      real (kind=R8) :: !djm Jan2017
     &  snatmp(nCv,0:1,0:ns-1), smotmp(nCv,0:3,0:ns-1),
     &  shetmp(nCv,0:3), shitmp(nCv,0:3)
*   ..procedures
      external xertst, ipgeti, b2xvsg
*   ..initialisation
      save diagno
      data diagno/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srsm')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
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
      call b2xvsg (nCv, nn, 1, 'nn', '.gt.')
      call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      call b2xvsg (nCv, te, 1, 'te', '.gt.')
      call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
      call xertst (maxval(dabs(ua)).lt.c, 'Supra-luminal velocities !')
*   ..obtain input
      if (ncall_b2srsm.eq.0) then
c The following switches are only used in 'WG-TODO' blocks, i.e. not yet converted to wide grid functionality
c        call ipgeti ('b2srsm_diagno', diagno)
        if(switch%b2srsm_enable.eq.0) then
          write(*,*) 'b2srsm disabled'
        endif
        if(eirene_mod.le.0) eirene_mod=1
      endif

cdjm Jan2017 Keep sources at the start to calculate total change for balance routines
      snatmp = sr%sna
      smotmp = sr%smo
      shetmp = sr%she
      shitmp = sr%shi

      if(switch%b2srsm_enable.ne.0) then
!WG_TODO
        call xerrab ('b2srsm not adapted for WG')
c*** calculate the integral sources for all regions:
c*** before correction (1), after correction (2)
        timefactor = 1.0_R8

*   ..smooth particle and parallel momentum sources
        do is = 0, ns-1
          do iCv = 1, mpg%nCv
!WG_TODO            if (resignore(iCv,1).ne.0) then
c check particles
              invtimescale =
     1         abs(sr%sna(iCv,0,is)+sr%sna(iCv,1,is)*na(iCv,is))
     2         /(na(iCv,is)*geo%cvVol(iCv))
              hlp = invtimescale*time_factor_required*
     1         dtim*dtco(is,mpg%cvReg(iCv))*time_factor(iCv)
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                  write(*,*) 'WARNING: particle sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,
     2             dtim*dtco(is,mpg%cvReg(iCv))*time_factor(iCv)
                  write(*,'(a,2i4,i3)') 'iCv,is = ', iCv,is
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(iCv)/ev,ti(iCv)/ev,ne(iCv)
                endif
                sr%sna(iCv,0,is) = sr%sna(iCv,0,is)/hlp
                sr%sna(iCv,1,is) = sr%sna(iCv,1,is)/hlp
                if (main_call.and.
     &           mod(ncall_b2srsm,eirene_mod).eq.0) then
                  src_ni(mpg%cvReg(iCv),2,b2espcr(is)) =
     &             src_ni(mpg%cvReg(iCv),2,b2espcr(is)) +
     &             cor_ni(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  src_ni(0,2,b2espcr(is)) =
     &             src_ni(0,2,b2espcr(is)) +
     &             cor_ni(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                endif
              endif
!WG_TODO            endif
          enddo
          do iCv = 1, mpg%nCv
!WG_TODO            if (resignore(iCv,2).ne.0) then
c check parallel momentum
              invtimescale = abs(sr%smo(iCv,0,is) +
     1         sr%smo(iCv,1,is)*ua(iCv,is) +
     &         sr%smo(iCv,2,is)*roxa(iCv,is) +
     &         sr%smo(iCv,3,is)*roxa(iCv,is)*ua(iCv,is))
     2         /(sqrt(ti(iCv)/
     3         (am(is)*mp))*roxa(iCv,is)*geo%cvVol(iCv))
              hlp = invtimescale*time_factor_required*
     1         dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)
              if (hlp.gt.1.0_R8) then
                timescale  = 1.0_R8/invtimescale
                timefactor = min(timefactor,1.0_R8/hlp)
                if ((main_call.and.diagno.ge.1).or.diagno.ge.2) then
                  write(*,*)
     1             'WARNING: parallel momentum sources very strong @'
                  write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,             timescale,
     2             dtim*dtmo(is,mpg%cvReg(iCv))*time_factor(iCv)
                  write(*,'(a,2i4,i3)') 'iCv,is = ', iCv,is
                  write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1             te(iCv)/ev,ti(iCv)/ev,ne(iCv)
                endif
                sr%smo(iCv,0,is) = sr%smo(iCv,0,is)/hlp
                sr%smo(iCv,1,is) = sr%smo(iCv,1,is)/hlp
                sr%smo(iCv,2,is) = sr%smo(iCv,2,is)/hlp
                sr%smo(iCv,3,is) = sr%smo(iCv,3,is)/hlp
                if (main_call.and.
     &            mod(ncall_b2srsm,eirene_mod).eq.0) then
                  src_mo(mpg%cvReg(iCv),2,b2espcr(is)) =
     &             src_mo(mpg%cvReg(iCv),2,b2espcr(is)) +
     &             cor_mo(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  src_mo(0,2,b2espcr(is)) =
     &             src_mo(0,2,b2espcr(is)) +
     &             cor_mo(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  src_ms(mpg%cvReg(iCv),2) =
     &             src_ms(mpg%cvReg(iCv),2) +
     &             cor_mo(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                  src_ms(0,2) = src_ms(0,2) +
     &             cor_mo(iCv,b2espcr(is))*(1.0_R8/hlp-1.0_R8)
                endif
              endif
!WG_TODO            endif
          enddo
        enddo
*   ..smooth electron heat, ion heat and electron particle sources
        do iCv = 1, mpg%nCv
!WG_TODO          if (resignore(iCv,1).ne.0) then
c check ion energy
            invtimescale =
     1       abs(sr%shi(iCv,0) + sr%shi(iCv,1)*ti(iCv) +
     &           sr%shi(iCv,2)*ni(iCv,0) +
     &           sr%shi(iCv,3)*ti(iCv)*ni(iCv,0))
     2         /(ti(iCv)*ni(iCv,0)*geo%cvVol(iCv))
            hlp = invtimescale*time_factor_required*
     .       dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
            if (hlp.gt.1.0_R8) then
              timescale  = 1.0_R8/invtimescale
              timefactor = min(timefactor,1.0_R8/hlp)
              if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                write(*,*) 'WARNING: ion heat sources very strong @'
                write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,           timescale,
     2           dtim*dtei(mpg%cvReg(iCv))*time_factor(iCv)
                write(*,'(a,2i4)') 'iCv = ', iCv
                write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1           te(iCv)/ev,ti(iCv)/ev,ne(iCv)
              endif
              sr%shi(iCv,0) = sr%shi(iCv,0)/hlp
              sr%shi(iCv,1) = sr%shi(iCv,1)/hlp
              sr%shi(iCv,2) = sr%shi(iCv,2)/hlp
              sr%shi(iCv,3) = sr%shi(iCv,3)/hlp
              if (main_call.and.
     &         mod(ncall_b2srsm,eirene_mod).eq.0) then
                src_ei(mpg%cvReg(iCv),2) =
     &           src_ei(mpg%cvReg(iCv),2) +
     &           cor_ei(iCv)*(1.0_R8/hlp-1.0_R8)
                src_ei(0,2) = src_ei(0,2) +
     &           cor_ei(iCv)*(1.0_R8/hlp-1.0_R8)
                src_es(mpg%cvReg(iCv),2) =
     &           src_es(mpg%cvReg(iCv),2) +
     &           cor_ei(iCv)*(1.0_R8/hlp-1.0_R8)
                src_es(0,2) = src_es(0,2) +
     &           cor_ei(iCv)*(1.0_R8/hlp-1.0_R8)
              endif
            endif
c check neutral energy
            if (switch%tn_style.ne.0) then
              call xerrab ('b2srsm not developed for tn_style.ne.0')
            endif
c check electron energy
            invtimescale =
     1       abs(sr%she(iCv,0) + sr%she(iCv,1)*te(iCv) +
     &           sr%she(iCv,2)*ne(iCv) +
     &           sr%she(iCv,3)*te(iCv)*ne(iCv))
     2         /(te(iCv)*ne(iCv)*geo%cvVol(iCv))
            hlp = invtimescale*time_factor_required*
     .       dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
            if (hlp.gt.1.0_R8) then
              timescale  = 1.0_R8/invtimescale
              timefactor = min(timefactor,1.0_R8/hlp)
              if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                write(*,*)
     ,           'WARNING: electron heat sources very strong @'
                write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,           timescale,
     2           dtim*dtee(mpg%cvReg(iCv))*time_factor(iCv)
                write(*,'(a,2i4)') 'iCv = ', iCv
                write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1           te(iCv)/ev,ti(iCv)/ev,ne(iCv)
              endif
              sr%she(iCv,0) = sr%she(iCv,0)/hlp
              sr%she(iCv,1) = sr%she(iCv,1)/hlp
              sr%she(iCv,2) = sr%she(iCv,2)/hlp
              sr%she(iCv,3) = sr%she(iCv,3)/hlp
              if (main_call.and.
     &         mod(ncall_b2srsm,eirene_mod).eq.0) then
                src_ee(mpg%cvReg(iCv),2) =
     &           src_ee(mpg%cvReg(iCv),2) +
     &           cor_ee(iCv)*(1.0_R8/hlp-1.0_R8)
                src_ee(0,2) = src_ee(0,2) +
     &           cor_ee(iCv)*(1.0_R8/hlp-1.0_R8)
                src_es(mpg%cvReg(iCv),2) =
     &           src_es(mpg%cvReg(iCv),2) +
     &           cor_ee(iCv)*(1.0_R8/hlp-1.0_R8)
               src_es(0,2) = src_es(0,2) +
     &           cor_ee(iCv)*(1.0_R8/hlp-1.0_R8)
              endif
            endif
c check electron particle number
            invtimescale =
     1       abs(sr%sne(iCv,0) + sr%sne(iCv,1)*ne(iCv))
     2         /(ne(iCv)*geo%cvVol(iCv))
            hlp = invtimescale*time_factor_required*
     1       dtim*time_factor(iCv)
            if (hlp.gt.1.0_R8) then
              timescale  = 1.0_R8/invtimescale
              timefactor = min(timefactor,1.0_R8/hlp)
              if ((main_call .and. diagno.ge.1) .or. diagno.ge.2) then
                write(*,*)
     ,           'WARNING: electron particle sources very strong @'
                write(*,'(a,2e10.3)') 'timescale,dt = ',
     ,           timescale,dtim*time_factor(iCv)
                write(*,'(a,2i4)') 'iCv = ', iCv
                write(*,'(a,3e10.3)') 'te,ti,ne = ',
     1           te(iCv)/ev,ti(iCv)/ev,ne(iCv)
              endif
              sr%sne(iCv,0) = sr%sne(iCv,0)/hlp
              sr%sne(iCv,1) = sr%sne(iCv,1)/hlp
              if (main_call.and.
     &          mod(ncall_b2srsm,eirene_mod).eq.0) then
                src_ne(mpg%cvReg(iCv),2) =
     &           src_ne(mpg%cvReg(iCv),2) +
     &           cor_ne(iCv)*(1.0_R8/hlp-1.0_R8)
                src_ne(0,2) = src_ne(0,2) +
     &           cor_ne(iCv)*(1.0_R8/hlp-1.0_R8)
              endif
            endif
!WG_TODO          endif
        enddo
      endif

#ifndef FIXED_POINT
cdjm Store the total change in linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2srsm_sna0to1 = sr%sna-snatmp
        b2srsm_smo0to3 = sr%smo-smotmp
        b2srsm_she0to3 = sr%she-shetmp
        b2srsm_shi0to3 = sr%shi-shitmp
      endif
#endif

      if (switch%b2npmo_iout.eq.1) then                                  !srv 11.09.09 {
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,sr%smo(:,0,is)+                        !srv 07.10.09
     &     sr%smo(:,1,is)*ua(:,is)+
     &     sr%smo(:,2,is)*na(:,is)*mp*am(is)+
     &     sr%smo(:,3,is)*ua(:,is)*na(:,is)*mp*am(is),
     &    'b2srsm_smosm'//chns)
        enddo
      endif
      if (switch%b2npco_iout.eq.1) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,sr%sna(:,0,is)+                        !srv 07.10.09
     &             sr%sna(:,1,is)*na(:,is),'b2srsm_snasm'//chns)
        enddo
      endif
      if (switch%b2npht_iout.eq.1) then
        call my_out_us(70,nCv,0,sr%she(:,0)+sr%she(:,1)*te+              !srv 07.10.09
     &   sr%she(:,2)*ne + sr%she(:,3)*te*ne,'b2srsm_shesm')
        call my_out_us(70,nCv,0,sr%shi(:,0)+sr%shi(:,1)*ti+
     &   sr%shi(:,2)*ni(:,0) +
     &   sr%shi(:,3)*ti*ni(:,0),'b2srsm_shism')
      endif

      if (switch%b2srsm_iout.eq.1) then
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
      if (ncall_b2srsm.eq.0 .or. main_call)
     &  ncall_b2srsm = ncall_b2srsm+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srsm

      contains

      function roxa(iCv, is)
      implicit none
      integer, intent(in) :: iCv, is
      real (kind=R8) :: roxa

      roxa = am(is)*mp*na(iCv,is)
      return
      end function roxa

      end subroutine b2srsm
