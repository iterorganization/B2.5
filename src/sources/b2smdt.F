*-----------------------------------------------------------------------
*.specification

      subroutine b2smdt (nx, ny, ns, dtim, hz,   !srv 09.01.01 01.07.08
     &  vol, na, uap, smodt, lout)               !srv 11.09.09 10.10.22
      use b2mod_types 
      use b2mod_indirect 
      use b2mod_numerics_namelist 
      use b2mod_constants 
      use b2mod_b2cmpa 
      use b2mod_boundary_namelist
     , , only : iftimbound
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_smo0to3, balance_netcdf
      use b2mod_subsys
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      logical :: lout                                               !srv 16.07.10
      integer nx, ny, ns
      real (kind=R8) ::
     *  dtim, vol(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),             !srv 10.10.22
     &  uap(-1:nx,-1:ny,0:ns-1), hz(-1:nx,-1:ny)                    !srv 09.01.01 10.10.22
*   ..input/output arguments
      real (kind=R8) ::
     *  smodt(-1:nx,-1:ny,0:3,0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SMDT modifies the source coefficients to account for the
*     implicit time discretisation for velocity equations.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, k, use_eirene
      integer, save :: iout = 0, iout_b2npmo = 0, iout_b2wdat = 0
      character chns*3, chk*1
      real (kind=R8) ::
     *   phm1, t0, ttim
      save ncall, phm1, use_eirene           !xpb
*   ..procedures
      external xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm1/1.0_R8/, use_eirene/0/                    !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2smdt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2srdt_phm1', phm1)
        call xertst (0.0_R8.le.phm1,'faulty internal parameter phm1')   !xpb
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2smdt_iout', iout)                            !srv 30.06.08
        call ipgeti ('b2npmo_iout', iout_b2npmo)
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
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')                          !srv 10.10.22
      endif
    
#ifdef WITH_FTIMINGS
       call timer%start("parallel")
#endif
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$OMP PARALLEL
!$OMP& default(none)
!$OMP& private(is,ix,iy,ttim,t0)
!$OMP& shared(nx,ny,ns,iftimbound,is_neutral,time_factor)
!$OMP& shared(leftix,rightix,bottomiy,topiy,region)
!$OMP& shared(use_eirene,dtim,dtmo)
!$OMP& shared(phm1)
!$OMP& shared(vol,na,uap,hz,am)
! these are the result arrays:
!$OMP& shared(smodt)
#endif
* ..include the contributions from dtim

*   ..modify smodt
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
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
            smodt(ix,iy,0,is) = t0*am(is)*mp*uap(ix,iy,is)*na(ix,iy,is)  !srv 06.01.21 { 10.10.22
            smodt(ix,iy,1,is) = 0.0_R8
            smodt(ix,iy,2,is) = 0.0_R8
            smodt(ix,iy,3,is) = -t0                                      !srv 06.01.21 }
        enddo
       enddo
      enddo
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$OMP END DO
#endif

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$OMP END PARALLEL
#endif
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
cdjm Store the total change in linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2srdt_smo0to3 = smodt
      endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated. 
      if (iout.eq.1 .and. lout) then
       do is = 0, ns-1
        write (chns,'(i3.3)') is
        do k=0,3
         write (chk,'(i1)') k
         call my_out(70,nx,ny,smodt(-1,-1,k,is),
     &                       'b2smdt_smodt'//chk//chns)
        enddo
       enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2smdt

      end subroutine b2smdt
