*-----------------------------------------------------------------------
*.specification

      subroutine b2sigp (nCv, nFc, nVx, isb, switch, geo, mpg, rzb,
     &  nb, ti, tn, ne, te, po, smbgp, wrk0, wrk1)
      use b2mod_types
      use b2mod_constants                                         !srv 23.03.10
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2mod_openmp
      use b2mod_b2cmpa
#ifndef FIXED_POINT
      use b2mod_balance                                           !djm Jan2017
     & , only : b2sigp_smogp0to3, b2sigp_pstat0to3, balance_netcdf,
     &          b2sigp_pstati0to3, b2sigp_pstate0to3, !djm Jul2018
     &          b2sigp_smogpi0to3, b2sigp_smogpe0to3, b2sigp_smogpgr0to3
#endif
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2sigp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) ::
     *  rzb(nCv),
     &  nb(nCv), ti(nCv), tn(nCv),
     &  ne(nCv), te(nCv), po(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8), intent(out) ::
     *   smbgp(nCv,0:3)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     * wrk0(nCv), wrk1(nCv), wrk3(nCv),
#ifdef WG_TODO
     * wrk4(nCv),                                  !djm Jan2017
#endif
     * wrk5(nCv), wrk6(nCv),                       !srv 28.03.17
     * tb(nCv), nete(nCv),
     * smbgp_b4gr(nCv,0:3)                         !djm Jul2018
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIGP computes the pressure gradient term on the right hand
*     of the parallel momentum balance equation for one species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc!, k
      character chns*3 !, chk*1                                       !srv 17.06.02
      real (kind=R8) ::
     *   t0, wrkf(nFc), wrkv(nVx), nbf(nFc), tbf(nFc), rzbf(nFc),
     *   tef(nFc)
*   ..procedures
      external xertst
      external b2xvsg
*   ..initialisation


*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sigp')
*   ..test nCv, nFv
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..extensive tests on first few calls
      if (ncall_b2sigp.lt.3) then
*    ..test sign of nb, ti, tn, ne, te, rzb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, rzb, 1, 'rzb', '.ge.')
      endif
*   .. check output writing
      if ((switch%b2sigp_iout.ne.0.or.switch%iout_b2wdat.eq.4)) then                       !srv 17.06.02 {
        if (in_parallel() .and. ncall_b2sigp == 0) then
#ifdef _OPENMP
!$OMP MASTER
#endif
           if (ncall_b2sigp==0) then
             write(*,*) 'b2sigp OpenMP warning: no file output ',
     &                 'in parallel mode'
           end if
#ifdef _OPENMP
!$OMP END MASTER
#endif
        end if 
      end if

* ..choose the correct temperature
      if (switch%tn_style.eq.2.and.is_neutral(isb).and.
     &      nint(zn(isb)).eq.1)then
        tb = tn !! tn = ti if no separate neutral energy equation is solved
      else
        tb = ti
      endif

* ..interpolate some variables to cell faces
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,nb,nbf)
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,tb,tbf)
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,rzb,rzbf)
      call intface(nCv,nFc,mpg%fcCv,geo%fcVol,te,tef)

      if (switch%b2sigp_phm0.ne.0.0_R8) then                         !srv 26.06.18
* ..compute pressure gradient term
*   ..compute smbgp
        if(switch%b2sigp_style.eq.0) then                              !srv 02.07.08
          call xerrab('b2sigp_style=0 not available for WG')
        elseif(switch%b2sigp_style.eq.1) then                           !srv 02.07.08 {
*   ..  compute gradients in cell centers
          wrk3 = nb*tb
          nete = ne*te
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk3, wrkv, wrk0)
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, nete, wrkv, wrk1)
* ..compute pressure gradient term
          smbgp = 0.0_R8
          wrk5  = 0.0_R8
          wrk6  = 0.0_R8
*   ..  compute smbgp
          do iCv = 1, mpg%nCi
           t0=switch%b2sigp_phm0*geo%cvHz(iCv)*geo%cvVol(iCv)*
     &        geo%cvBb(iCv,0)/geo%cvBb(iCv,3) !srv 18.10.02 10.06.05
           wrk5(iCv) = -t0*wrk0(iCv)
           wrk6(iCv) = -t0*rzb(iCv)*(nb(iCv)/ne(iCv))*wrk1(iCv)
           smbgp(iCv,0) = wrk5(iCv) + wrk6(iCv)
          enddo
#ifndef FIXED_POINT
cdjm compute b2sigp_pstat0to3 on cell face
          if (balance_netcdf.ne.0) then
            b2sigp_pstat0to3 = 0.0_R8
            b2sigp_pstati0to3 = 0.0_R8
!$omp critical
            b2sigp_pstate0to3 = 0.0_R8
            do iFc = 1, nFc
              t0 = switch%b2sigp_phm0*geo%fcHz(iFc)*geo%fcS(iFc)*
     &         geo%fcBb(iFc,0)/geo%fcBb(iFc,3)*geo%fcQalf(iFc,0)
              b2sigp_pstat0to3(iFc,0,isb) = t0*(nbf(iFc)*tbf(iFc)+
     &         rzbf(iFc)*nbf(iFc)*tef(iFc))
              b2sigp_pstati0to3(iFc,0,isb) = t0*nbf(iFc)*tbf(iFc)
              b2sigp_pstate0to3(iFc,0) = t0*rzbf(iFc)*nbf(iFc)*tef(iFc)
            enddo
!$omp end critical
          endif
#endif
        else                                                            !srv 23.03.10 {
          call xertst(switch%pot_eq.eq.1,
     &   'b2sigp_style.eq.2 requires pot_eq.eq.1!')
* ..compute pressure gradient term
          smbgp = 0.0_R8
          wrk3  = nb*tb
          wrk5  = 0.0_R8
          wrk6  = 0.0_R8
*   ..  compute gradients in cell centers
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, wrk3, wrkv, wrk0)
          call gradc_p(nCv, nFc, nVx, 0, geo, mpg, po, wrkv, wrk1)
*   ..  compute smbgp
          do iCv = 1, mpg%nCi
            t0=switch%b2sigp_phm0*geo%cvHz(iCv)*geo%cvVol(iCv)*
     &         geo%cvBb(iCv,0)/geo%cvBb(iCv,3)
            wrk5(iCv) = -t0*wrk0(iCv)
            wrk6(iCv) = -t0*qe*rzb(iCv)*nb(iCv)*
     &                   wrk1(iCv)
            smbgp(iCv,0) = wrk5(iCv) + wrk6(iCv)
          enddo                                                         !srv 23.03.10 }
#ifndef FIXED_POINT
cdjm compute b2sigp_pstat0to3 on cell face
          if (balance_netcdf.ne.0) then
            b2sigp_pstat0to3 = 0.0_R8
!WG_TODO          do iFc = 1,mpg%nFc
!WG_TODO             t0=switch%b2sigp_phm0*geo%fcHz(iFc)*geo%fcS(iFc)*geo%fcBb(iFc,0)/
!WG_TODO      &        geo%fcBb(iFc,3)*geo%fcQalf(iFc,0)
!WG_TODO             b2sigp_pstat0to3(iFc,0,isb) =
!WG_TODO      &         -t0*wrk0(ix,iy)
!WG_TODO      &         -t0*qe*wrk1(ix,iy)*
!WG_TODO      &         (rzb(ix,iy)*nb(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy))+
!WG_TODO      &          rzb(leftix(ix,iy),leftiy(ix,iy))*
!WG_TODO      &          nb(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy))/
!WG_TODO      &         (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
!WG_TODO               b2sigp_pstat0to3(ix,iy,1,isb) = 0.0_R8
!WG_TODO               b2sigp_pstat0to3(ix,iy,2,isb) = 0.0_R8
!WG_TODO               b2sigp_pstat0to3(ix,iy,3,isb) = 0.0_R8
!WG_TODO          enddo
          endif
#endif
        endif                                                           !srv 02.07.08 }

        smbgp_b4gr=smbgp !djm Jul2018
cxpb apply restriction to pressure gradient for minority ions
        if (switch%b2sigp_pressure_restriction.ne.0) then
          wrk0 = nb*tb
          call intface(nCv,nFc,mpg%fcCv,geo%fcHc,wrk0,wrkf)
          call intcell(nFc,nCv,mpg,mpg%intcellP,wrkf,wrk1)
          do iCv = 1, mpg%nCi
            if (nb(iCv)/ne(iCv) .lt. 0.1_R8) then
              smbgp(iCv,0) = smbgp(iCv,0)*
     &        min(1.0_R8,wrk0(iCv)/wrk1(iCv))
            endif
          enddo                                                          !lk 23.06.17}
        endif

      else          !!switch%b2sigp_phm0.eq.0.0_R8                                                  !srv 26.06.18 {
        smbgp = 0.0_R8
#ifdef WG_TODO
        b2sigp_pstat0to3(:, :, isb) = 0.0_R8
#endif
        wrk5 = 0.0_R8
        wrk6 = 0.0_R8
      endif                                                           !srv 26.06.18 }

      if (switch%b2sigp_iout.ne.0.or.switch%iout_b2wdat.eq.4) then    !srv 17.06.02 {
        if (.not. in_parallel()) then 
          write (chns,'(i3.3)') isb
    !      do k=0,3
    !        write (chk,'(i1)') k
    !        call my_out_us(70,nCv,0,smbgp(-1,-1,k),'b2sigp_smogp'//chk//chns)
    !      enddo
          call my_out_us(70,nCv,0,smbgp(1,0),'b2sigp_smogp'//chns)
        if(switch%b2sigp_style.eq.0 .or. switch%b2sigp_style.eq.1 ) then
              call my_out_us(70,nCv,0,wrk5,'b2sigp_smogpi'//chns)
              call my_out_us(70,nCv,0,wrk6,'b2sigp_smogpe'//chns)
          elseif(switch%b2sigp_style.eq.2) then
              call my_out_us(70,nCv,0,wrk5,'b2sigp_smogpi'//chns)
              call my_out_us(70,nCv,0,wrk6,'b2sigp_smogpo'//chns)
          endif
        endif
      endif                                                           !srv 17.06.02 }
#ifndef FIXED_POINT
cdjm Jan2017 Store linearised source for balance
      if (balance_netcdf.ne.0) then
        b2sigp_smogp0to3(1:nCv,0:3,isb)=smbgp
        b2sigp_smogpi0to3(1:nCv,0,isb)=wrk5 !djm Jul2018
        b2sigp_smogpe0to3(1:nCv,0,isb)=wrk6 !djm Jul2018
        b2sigp_smogpgr0to3(1:nCv,0:3,isb)=smbgp-smbgp_b4gr
      endif
#endif

* ..return
      ncall_b2sigp = ncall_b2sigp+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sigp

      end subroutine b2sigp
