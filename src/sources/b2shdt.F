*-----------------------------------------------------------------------
*.specification

      subroutine b2shdt (nx, ny, ns, dtim, BoRiS,                       !srv 07.07.21
     &  vol, na0, tep, ti0, ne0, ni0, kinrgy0, shedt, shidt, lout)
      use b2mod_types 
      use b2mod_indirect 
      use b2mod_numerics_namelist 
      use b2mod_constants 
      use b2mod_b2cmpa 
      use b2mod_boundary_namelist
     , , only : iftimbound
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_she0to3, b2srdt_shi0to3, balance_netcdf
      use b2mod_subsys
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      logical(1) :: lout                                             !srv 16.07.10
      integer nx, ny, ns
      real (kind=R8) ::
     *  dtim, BoRiS, vol(-1:nx,-1:ny), na0(-1:nx,-1:ny,0:ns-1),
     &  tep(-1:nx,-1:ny), ti0(-1:nx,-1:ny),
     &  ne0(-1:nx,-1:ny), ni0(-1:nx,-1:ny,0:1),
     &  kinrgy0(-1:nx,-1:ny,0:ns-1)
*   ..input/output arguments
      real (kind=R8) ::
     &  shedt(-1:nx,-1:ny,0:3), shidt(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SHDT modifies the source coefficients to account for the
*     implicit time discretisation for Te and Ti equations.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, k, use_eirene
      integer, save :: iout = 0, iout_b2nph9 = 0, iout_b2npht = 0,   !srv 11.09.09
     &                 iout_b2wdat = 0
      character chk*1                                                !srv 11.09.09 }
      real (kind=R8) ::
     *   phm3, t0, ttim
      save ncall, phm3, use_eirene                                   !xpb
*   ..procedures
      external xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm3/1.0_R8/, use_eirene/0/                    !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2shdt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2srdt_phm3', phm3)
        call xertst (0.0_R8.le.phm3,'faulty internal parameter phm3')   !xpb
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2shdt_iout', iout)                            !srv 30.06.08
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
       call b2xvsg (n2, tep, 1, 'tep', '.gt.')
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
!$OMP& shared(use_eirene,dtim,dtee,dtei)
!$OMP& shared(phm3,BoRiS,kinrgy0)
!$OMP& shared(vol,na0,ne0,ni0,ti0,tep,am)
! these are the result arrays:
!$OMP& shared(shedt,shidt)
#endif

*   ..modify shedt, shidt
#ifndef NO_OPENMP_B2SRDT
!$OMP DO
#endif
      do iy = -1, ny
        do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
          if(iftimbound(ix,iy) .ne. 0) cycle

          ttim=dtim*dtee(region(ix,iy,0))*time_factor(ix,iy)
          t0 = (phm3/ttim)*vol(ix,iy)
          shedt(ix,iy,0) = 1.5_R8*t0*ne0(ix,iy)*tep(ix,iy)               !srv 06.01.21 {
          shedt(ix,iy,1) = 0.0_R8
          shedt(ix,iy,2) = 0.0_R8
          shedt(ix,iy,3) = -1.5_R8*t0                                    !srv 06.01.21 }
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
          shidt(ix,iy,0) = 1.5_R8*t0*ni0(ix,iy,0)*ti0(ix,iy)             !srv 06.01.21 {
          shidt(ix,iy,1) = 0.0_R8
          shidt(ix,iy,2) = 0.0_R8
          shidt(ix,iy,3) = -1.5_R8*t0                                    !srv 06.01.21 }
          if (BoRiS.eq.1.0_R8) then
            do is = 0, ns-1
              shidt(ix,iy,0) = shidt(ix,iy,0)+t0*
     &         (kinrgy0(ix,iy,is)+ti0(ix,iy))*na0(ix,iy,is)
              shidt(ix,iy,2) = shidt(ix,iy,2)-t0*
     &         (kinrgy0(ix,iy,is)+ti0(ix,iy))*na0(ix,iy,is)/
     &         ni0(ix,iy,0)
            enddo
          endif
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
        b2srdt_she0to3 = shedt
        b2srdt_shi0to3 = shidt
      endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated. 
      if (iout.eq.1 .and. lout) then
       do k=0,3
        write (chk,'(i1)') k
        call my_out(70,nx,ny,shedt(-1,-1,k),'b2shdt_shedt'//chk)
        call my_out(70,nx,ny,shidt(-1,-1,k),'b2shdt_shidt'//chk)
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2shdt

      end subroutine b2shdt
