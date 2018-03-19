*-----------------------------------------------------------------------
*.specification

      subroutine b2sian (nx, ny, isb, bb, vol, hx, hy, hz, qz, fchanml, !srv 25.09.17 {
     &           smban)
      use b2mod_types
      use b2mod_constants
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_balance !djm Jan2017
     & , only : b2sian_smo0to3, balance_netcdf
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     &   bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),
     &   hy(-1:nx,-1:ny), hz(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1),
     &   fchanml(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smban(-1:nx,-1:ny,0:3)
*   ..workspace
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SIAN computes contribution to the parallel momentum source due
*     to the anomalous current
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, k
      integer, save :: iout = 0
      character chns*3, chark*1
      real (kind=R8) ::
     *    phm0
      save ncall, phm0
*   ..procedures
      real (kind=R8) ::
     *   hy1
      external subini, subend, xertst, ipgeti, ipgetr
      external b2xvff, b2xvsg
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*   ..initialization
      data ncall/0/, phm0/1.0e0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sian')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2sian_phm0', phm0)
       call ipgeti ('b2sian_iout', iout)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      smban = 0.0_R8
      if(nint(zn(isb)).eq.1.and.nint(zamax(isb)).eq.1) then
*   ..  compute smban
        do iy = -1, ny
          do ix = -1, nx
           if(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1 .or.
     &       bottomiy(ix,iy).eq.-2.or.topiy(ix,iy).eq.ny+1) then
            smban(ix,iy,0) = 0.0e0_R8
            smban(ix,iy,1) = 0.0e0_R8
            smban(ix,iy,2) = 0.0e0_R8
            smban(ix,iy,3) = 0.0e0_R8
           else
            smban(ix,iy,0) = -phm0*hz(ix,iy)*hy1(ix,iy)*
     &          bb(ix,iy,0)/bb(ix,iy,2)*bb(ix,iy,3)*
     &          (fchanml(ix,iy,1,1)+
     &           fchanml(topix(ix,iy),topiy(ix,iy),1,1))*0.5_R8
            smban(ix,iy,1) = 0.0e0_R8
            smban(ix,iy,2) = 0.0e0_R8
            smban(ix,iy,3) = 0.0e0_R8
           endif !if(left...
          enddo
        enddo      !iyv }
      endif
cdjm Store linearised source for balance
      if (balance_netcdf.ne.0) then
        b2sian_smo0to3(-1:nx,-1:ny,0:3,isb)=smban
      endif

      if (iout.ne.0) then
        write (chns,'(i3.3)') isb
        do k=0,3
          write (chark,'(i1)') k
          call my_out(70,nx,ny,smban(-1,-1,k),
     .                         'b2sian_smoan'//chark//chns)
        enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sian

      end subroutine b2sian
