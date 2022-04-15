*-----------------------------------------------------------------------
*.specification

      subroutine b2srdt (nx, ny, ns, guard_flows, dtim, BoRiS, hz,   !srv 09.01.01 01.07.08
     &  vol, na0, ua0, te0, ti0, ne0, ni0, kinrgy0,
     &  na, ua, te, ti, ne, ni, kinrgy,                              !srv 16.07.10
     &  snadt, smodt, shedt, shidt, schdt, snedt, lout)              !srv 11.09.09
      use b2mod_types
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_boundary_namelist
     , , only : iftimbound
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_sna0to1,b2srdt_smo0to3,b2srdt_she0to3,
     2          b2srdt_shi0to3, balance_netcdf
      use b2mod_subsys
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      logical(1) :: lout                                             !srv 16.07.10
      integer nx, ny, ns, guard_flows                                !srv 01.07.08
      real (kind=R8) ::
     *  dtim, BoRiS, vol(-1:nx,-1:ny), na0(-1:nx,-1:ny,0:ns-1),
     &  ua0(-1:nx,-1:ny,0:ns-1), te0(-1:nx,-1:ny), ti0(-1:nx,-1:ny),
     &  ne0(-1:nx,-1:ny), ni0(-1:nx,-1:ny,0:1), hz(-1:nx,-1:ny),     !srv 09.01.01
     &  kinrgy0(-1:nx,-1:ny,0:ns-1),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),              !srv 16.07.10
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &  ni(-1:nx,-1:ny,0:1), kinrgy(-1:nx,-1:ny,0:ns-1)
*   ..input/output arguments
      real (kind=R8) ::
     *  snadt(-1:nx,-1:ny,0:1,0:ns-1), smodt(-1:nx,-1:ny,0:3,0:ns-1),
     &  shedt(-1:nx,-1:ny,0:3), shidt(-1:nx,-1:ny,0:3),
     &  schdt(-1:nx,-1:ny,0:3), snedt(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SRDT modifies the source coefficients to account for the
*     implicit time discretisation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, k, use_eirene
      integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,            !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo = 0, iout_b2npp7 = 0,
     &                 iout_b2nppo = 0, iout_b2wdat = 0
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *   phm0, phm1, phm3, phm4, phm5, t0, ttim
      real (kind=R8) :: wrk0(-1:nx,-1:ny)
      save ncall, phm0, phm1, phm3, phm4, phm5, use_eirene           !xpb
      real (kind=R8) :: !djm Jan2017
     &  snatmp(-1:nx,-1:ny,0:1,0:ns-1), smotmp(-1:nx,-1:ny,0:3,0:ns-1), 
     &  shetmp(-1:nx,-1:ny,0:3), shitmp(-1:nx,-1:ny,0:3)
*   ..procedures
      external xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, phm1/1.0_R8/, phm3/1.0_R8/,
     1  phm4/0.0_R8/, phm5/1.0_R8/, use_eirene/0/                    !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2srdt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2srdt_phm0', phm0)
        call ipgetr ('b2srdt_phm1', phm1)
        call ipgetr ('b2srdt_phm3', phm3)
        call ipgetr ('b2srdt_phm4', phm4)                            !xpb
        call ipgetr ('b2srdt_phm5', phm5)                            !xpb
        call xertst (0.0_R8.le.phm0.and.0.0_R8.le.phm1               !xpb
     &   .and.0.0_R8.le.phm3.and.0.0_R8.le.phm4.and.0.0_R8.le.phm5,  !xpb
     &   'faulty internal parameter phm0, phm1, phm3, phm4, phm5')   !xpb
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2srdt_iout', iout)                            !srv 30.06.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npc7_iout', iout_b2npc7)                     !srv 02.04.18
        call ipgeti ('b2npc11_iout', iout_b2npc11)                   !srv 02.04.18
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2npp7_iout', iout_b2npp7)
        call ipgeti ('b2nppo_iout', iout_b2nppo)
        call ipgeti ('b2nph9_iout', iout_b2nph9)
        call ipgeti ('b2npht_iout', iout_b2npht)                     !srv 11.09.09 }
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (n2*ns, na0, 1, 'na0', '.gt.')
       call b2xvsg (n2, te0, 1, 'te0', '.gt.')
       call b2xvsg (n2, ti0, 1, 'ti0', '.gt.')
       call b2xvsg (n2, ne0, 1, 'ne0', '.gt.')
       call b2xvsg (2*n2, ni0, 1, 'ni0', '.gt.')
       call b2xvsg (n2*ns, kinrgy0, 1, 'kinrgy0', '.ge.')
      endif

#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
#ifndef NO_OPENMP_B2SRDT
!$OMP PARALLEL
!$OMP& default(none)
!$OMP& private(is,ix,iy,ttim,t0)
!$OMP& shared(nx,ny,ns,iftimbound,is_neutral,time_factor)
!$OMP& shared(leftix,rightix,bottomiy,topiy,region)
!$OMP& shared(guard_flows,use_eirene,dtim,dtco,dtmo,dtee,dtei)
!$OMP& shared(phm0,phm1,phm3,phm4,phm5,BoRiS,kinrgy0)
!$OMP& shared(vol,na0,ne0,ni0,ua0,ti0,te0,hz,am)
! these are the result arrays:
!$OMP& shared(snadt,smodt,shedt,shidt,schdt,snedt)
!$OMP& shared(snatmp,smotmp,shetmp,shitmp)
#endif
cdjm Jan2017 Keep sources at the start to calculate total change for balance routines
      snatmp = snadt
      smotmp = smodt
      shetmp = shedt
      shitmp = shidt
* ..include the contributions from dtim
*   ..modify snadt

#ifndef NO_OPENMP_B2SRDT
!$OMP DO
#endif
      do is = 0, ns-1
        do iy = -1, ny
          do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
           if(iftimbound(ix,iy) .ne. 0) cycle
!           if((guard_flows.eq.2 .or.
!     &        (is_neutral(is) .and. use_eirene.ne.0)) .and.
!     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
           ttim=dtim*dtco(is,region(ix,iy,0))*time_factor(ix,iy)
           t0=(phm0/ttim)*vol(ix,iy)                                     !srv 06.01.21
           snadt(ix,iy,0,is) = snadt(ix,iy,0,is)+t0*na0(ix,iy,is)
           snadt(ix,iy,1,is) = snadt(ix,iy,1,is)-t0
          enddo
        enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO NOWAIT
#endif

*   ..modify smodt
#ifndef NO_OPENMP_B2SRDT
!$OMP DO
#endif
      do is = 0, ns-1
        do iy = -1, ny
*       smodt(-1,iy,0:3,is) = (unchanged at 0.0)
          do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
            if(iftimbound(ix,iy) .ne. 0) cycle

            ttim=dtim*dtmo(is,region(ix,iy,0))*time_factor(ix,iy)
            t0 = (phm1/ttim)*vol(ix,iy)*hz(ix,iy)                        !srv 09.01.01
            smodt(ix,iy,0,is) = smodt(ix,iy,0,is)+
     1       t0*am(is)*mp*ua0(ix,iy,is)*
     &       na0(ix,iy,is)
*        smodt(ix,iy,1,is) = (unchanged)
*        smodt(ix,iy,2,is) = (unchanged)
            smodt(ix,iy,3,is) = smodt(ix,iy,3,is)-t0
        enddo
       enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO
#endif

*   ..modify shedt, shidt, schdt, snedt
#ifndef NO_OPENMP_B2SRDT
!$OMP DO
#endif
      do iy = -1, ny
        do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
          if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtee(region(ix,iy,0))*time_factor(ix,iy)
          t0 = (phm3/ttim)*vol(ix,iy)
          shedt(ix,iy,0) = shedt(ix,iy,0)+
     &     1.5_R8*t0*ne0(ix,iy)*te0(ix,iy)
*         shedt(ix,iy,1) = (unchanged)
*         shedt(ix,iy,2) = (unchanged)
          shedt(ix,iy,3) = shedt(ix,iy,3)-1.5_R8*t0
       enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO nowait

!$OMP DO
#endif
      do iy = -1, ny
        do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
          if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtei(region(ix,iy,0))*time_factor(ix,iy)
          t0 = (phm3/ttim)*vol(ix,iy)
          shidt(ix,iy,0) = shidt(ix,iy,0)+
     &     1.5_R8*t0*ni0(ix,iy,0)*ti0(ix,iy)
          if (BoRiS.eq.1.0_R8) then
            do is = 0, ns-1
              shidt(ix,iy,0) = shidt(ix,iy,0)+t0*
     &         kinrgy0(ix,iy,is)*na0(ix,iy,is)
              shidt(ix,iy,2) = shidt(ix,iy,2)-t0*
     &         kinrgy0(ix,iy,is)*na0(ix,iy,is)/
     &         ni0(ix,iy,0)
            enddo
          endif
*         shidt(ix,iy,1) = (unchanged)
*         shidt(ix,iy,2) = (otherwise unchanged)
          shidt(ix,iy,3) = shidt(ix,iy,3)-1.5_R8*t0
       enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO nowait

!$OMP DO
#endif
      do iy = -1, ny
        do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
          if(iftimbound(ix,iy) .ne. 0) cycle


cxpb  SCH should normally remain unchanged (phm4.eq.0.0)
cxpb  This treatment added as a fix to get the ExB runs going
          t0 = (phm4/dtim)*vol(ix,iy)*qe*time_factor(ix,iy)              !xpb
          schdt(ix,iy,0) = schdt(ix,iy,0)+t0*(ni0(ix,iy,0)-ne0(ix,iy))   !xpb
*         schdt(ix,iy,1) = (unchanged)
*         schdt(ix,iy,2) = (unchanged)
          schdt(ix,iy,3) = schdt(ix,iy,3)-t0                             !xpb
       enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO nowait

!$OMP DO
#endif
      do iy = -1, ny
        do ix = -1, nx
          if(iftimbound(ix,iy) .ne. 0) cycle
          t0 = (phm5/dtim)*vol(ix,iy)*time_factor(ix,iy)                 !xpb
          snedt(ix,iy,0) = snedt(ix,iy,0)+t0*ne0(ix,iy)
          snedt(ix,iy,1) = snedt(ix,iy,1)-t0
       enddo
      enddo
#ifndef NO_OPENMP_B2SRDT
!$OMP END DO

!$OMP END PARALLEL
#endif
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif

cdjm Store the total change in linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2srdt_sna0to1 = snadt-snatmp
        b2srdt_smo0to3 = smodt-smotmp
        b2srdt_she0to3 = shedt-shetmp
        b2srdt_shi0to3 = shidt-shitmp
      endif

      if ((iout_b2npmo.eq.1 .or. iout_b2wdat .eq. 4) .and. lout) then                                     !srv 11.09.09 {
       do is = 0, ns-1 
         write (chns,'(i3.3)') is
         wrk0=smodt(:,:,0,is)+
     &        smodt(:,:,1,is)*ua(:,:,is)+
     &        smodt(:,:,2,is)*na(:,:,is)*mp*am(is)+
     &        smodt(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
         call my_out(70,nx,ny,wrk0,'b2srdt_smodt'//chns)
       enddo
      endif
      if ((iout_b2npc9.eq.1 .or. iout_b2npco.eq.1 .or. iout_b2wdat.eq.4)
     &             .and. lout) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          wrk0=snadt(:,:,0,is)+snadt(:,:,1,is)*na(:,:,is)
          call my_out(70,nx,ny,wrk0,'b2srdt_snadt'//chns)
        enddo
      endif
      if ((iout_b2nph9.eq.1 .or. iout_b2npht.eq.1 .or. iout_b2wdat.eq.4)
     &             .and. lout) then
        wrk0=shedt(:,:,0)+shedt(:,:,1)*te+
     &       shedt(:,:,2)*ne + shedt(:,:,3)*te*ne
        call my_out(70,nx,ny,wrk0,'b2srdt_shedt')
        wrk0=shidt(:,:,0)+shidt(:,:,1)*ti+
     &       shidt(:,:,2)*ni(:,:,0)+
     &       shidt(:,:,3)*ti*ni(:,:,0)
        call my_out(70,nx,ny,wrk0,'b2srdt_shidt')
      endif

      if (iout.eq.1 .and. lout) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
         write (chk,'(i1)') k
         call my_out(70,nx,ny,smodt(-1,-1,k,is),
     &                       'b2srdt_smodt'//chk//chns)
        enddo
        do k=0,1
         write (chk,'(i1)') k
         call my_out(70,nx,ny,snadt(-1,-1,k,is),
     &                       'b2srdt_snadt'//chk//chns)
        enddo
       enddo
       do k=0,3
        write (chk,'(i1)') k
        call my_out(70,nx,ny,shedt(-1,-1,k),'b2srdt_shedt'//chk)
        call my_out(70,nx,ny,shidt(-1,-1,k),'b2srdt_shidt'//chk)
        call my_out(70,nx,ny,schdt(-1,-1,k),'b2srdt_schdt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2srdt

      end subroutine b2srdt
