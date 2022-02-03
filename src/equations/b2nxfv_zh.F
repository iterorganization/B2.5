*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfv_zh (nx, ny, rho_rel, hz, cvsbhz_cl, flubv)
      use b2mod_types
      use b2mod_indirect
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) :: rho_rel(-1:nx,-1:ny),
     *   hz(-1:nx,-1:ny), cvsbhz_cl(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   flubv(-1:nx,-1:ny,0:1)

*----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFV_ZH calculates the convective term due to new form of ion
*     viscosity term in the momentum correction equation if Zhdanov
*     closure is applied.
*
*----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy
      integer, save :: zh_vis_type = 0                              !som 13.01.22
      real (kind=R8) ::
     *   ft
      save ncall
      !$OMP THREADPRIVATE(ncall)
      external xertst, ipgeti, ipgetr, b2scopy
      external b2xvff, b2xvfx, b2xvfv, b2xvsg
*   ..initialization
      data ncall/0/

*----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfv_zh')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2tral_zh_vis_type', zh_vis_type)              !som 13.01.22
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif

* ..compute convection term
      ft = 4.0_R8/3.0_R8
*   ..compute flubvx
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).ne.-2) then
          if (zh_vis_type.eq.0) then                          !som 13.01.22
           flubv(ix,iy,0) = ft*cvsbhz_cl(ix,iy,0)*
     &     (log(hz(ix,iy))-log(hz(leftix(ix,iy),leftiy(ix,iy))))
          else if (zh_vis_type.eq.1) then
           flubv(ix,iy,0) = ft*cvsbhz_cl(ix,iy,0)*
     &     ((log(hz(ix,iy))-log(hz(leftix(ix,iy),leftiy(ix,iy))))-
     &      (log(rho_rel(ix,iy))-
     &       log(rho_rel(leftix(ix,iy),leftiy(ix,iy)))))
          endif
         else
          flubv(ix,iy,0) = 0.0e0_R8
         endif
        enddo
       enddo

       flubv(:,:,1) = 0.0_R8

* ..return
      ncall = ncall+1
      call subend ()
      return

*----------------------------------------------------------------------
*.end b2nxfv_zh

      end subroutine b2nxfv_zh
