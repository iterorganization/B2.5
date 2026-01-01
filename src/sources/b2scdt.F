*-----------------------------------------------------------------------
*.specification

      subroutine b2scdt (nx, ny, isb, guard_flows, dtim,      !srv 07.07.21
     &  vol, nap, snadt, lout)                                !srv 10.10.22
      use b2mod_types 
      use b2mod_indirect 
      use b2mod_numerics_namelist 
      use b2mod_constants 
      use b2mod_b2cmpa 
      use b2mod_boundary_namelist
     , , only : iftimbound
      use b2mod_balance !djm Jan2017
     1 , only : b2srdt_sna0to1, balance_netcdf
      use b2mod_subsys
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      logical :: lout                                                 !srv 16.07.10
      integer nx, ny, isb, guard_flows                                !srv 01.07.08
      real (kind=R8) ::
     *  dtim, vol(-1:nx,-1:ny), nap(-1:nx,-1:ny)                      !srv 10.10.22
*   ..input/output arguments
      real (kind=R8) ::
     *  snadt(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SCDT modifies the source coefficients to account for the
*     implicit time discretisation for the isb-th density equation.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, k, use_eirene
      integer, save :: iout = 0, iout_b2npc9 = 0, iout_b2npco = 0,   !srv 11.09.09 {
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,            !srv 02.04.18
     &                 iout_b2wdat = 0
      character chns*3, chk*1                                        !srv 11.09.09 }
      real (kind=R8) ::
     *   phm0, t0, ttim
      save ncall, phm0, use_eirene           !xpb
*   ..procedures
      external xertst, ipgetr
      external b2xvsg
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, use_eirene/0/                    !xpb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2scdt')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgetr ('b2srdt_phm0', phm0)
        call xertst (0.0_R8.le.phm0,'faulty internal parameter phm0')   !xpb
        call ipgeti ('b2mndr_eirene', use_eirene)
        call ipgeti ('b2scdt_iout', iout)                            !srv 30.06.08
        call ipgeti ('b2npc9_iout', iout_b2npc9)                     !srv 11.09.09 {
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npc7_iout', iout_b2npc7)                     !srv 02.04.18
        call ipgeti ('b2npc11_iout', iout_b2npc11)                   !srv 02.04.18
        call ipgeti ('b2wdat_iout', iout_b2wdat)                     !iys 17.04.17
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.le.isb, 'faulty argument isb')
*   ..test sign of dtim
      call xertst (0.0_R8.lt.dtim, 'faulty argument dtim')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test old state
       call b2xvsg (n2, nap, 1, 'nap', '.gt.')                       !srv 10.10.22
      endif
    
#ifdef WITH_FTIMINGS
       call timer%start("parallel")
#endif
#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$OMP PARALLEL
!$OMP& default(none)
!$OMP& private(ix,iy,ttim,t0)
!$OMP& shared(nx,ny,isb,iftimbound,is_neutral,time_factor)
!$OMP& shared(leftix,rightix,bottomiy,topiy,region)
!$OMP& shared(guard_flows,use_eirene,dtim,dtco)
!$OMP& shared(phm0)
!$OMP& shared(vol,nap)
! these are the result arrays:
!$OMP& shared(snadt)
#endif
* ..include the contributions from dtim
*   ..modify snadt

#if ( defined(_OPENMP) && !defined(NO_OPENMP_B2SRDT) )
!$OMP DO
#endif
        do iy = -1, ny
          do ix = -1, nx
csw 24feb2014 check for guard cell (SOLPS4)
           if(iftimbound(ix,iy) .ne. 0) cycle
           if((guard_flows.eq.2 .or.
     &        (is_neutral(isb) .and. use_eirene.ne.0)) .and. 
     &       (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1)) cycle
           ttim=dtim*dtco(isb,region(ix,iy,0))*time_factor(ix,iy)
           t0=(phm0/ttim)*vol(ix,iy)                                     !srv 06.01.21
           snadt(ix,iy,0) = t0*nap(ix,iy)                                !srv 06.01.21 10.10.22
           snadt(ix,iy,1) =-t0                                           !srv 06.01.21
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
        b2srdt_sna0to1(:,:,:,isb) = snadt
      endif

*    The diagnostic output of the time derivative will be done
*    in the subroutines where the solution is calculated. 
      if (iout.eq.1 .and. lout) then
        write (chns,'(i3.3)') isb
        do k=0,1
         write (chk,'(i1)') k
         call my_out(70,nx,ny,snadt(-1,-1,k),
     &                       'b2scdt_snadt'//chk//chns)
        enddo
      endif                                                           !srv 11.09.09 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2scdt

      end subroutine b2scdt
