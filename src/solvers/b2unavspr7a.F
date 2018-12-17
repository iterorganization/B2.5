      subroutine b2unavspr7a(nx, ny, nregionv, solvereg, isb, zb,       !srv 02.04.18 {
     &              fac_vis, ti, nb, pb, aa7)
      use b2mod_types
      use b2mod_indirect
      use b2mod_geo
      use b2mod_constants
      use b2mod_openmp
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv, isb
      logical solvereg(0:nregionv)
      real (kind=R8) ::
     &   zb(-1:nx,-1:ny), fac_vis(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), nb(-1:nx,-1:ny), pb(-1:nx,-1:ny)
*   ..output arguments (unchanged on exit)
      real (kind=R8) ::
     &   aa7(-1:nx,-1:ny,-3:3)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UNAVSPR7A evaluates the coefficients of matrix of the contribution
*     from divergence of the viscosity perpendicular current by analytical
*     manipulation with respect to main ion density using a 7-point stencil.
*     The coefficients will be stored in the matrix.
*
*   ..local variables
      integer, save :: ncall = 0, iout = 0
      integer ix, iy, ixt, ixb, iyt, iyb, ireg
      real (kind=R8), save :: vis_per = 1.0_R8
      real (kind=R8) ::
     &   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),
     &   cvisper(-1:nx,-1:ny,0:1), dvisper(-1:nx,-1:ny,0:1),
     &   kmpr(-1:nx,-1:ny)
      real (kind=R8) ::
     &   aa7p2, aa7m2, coef
!$OMP THREADPRIVATE(ncall, iout, vis_per)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2unavspr7a')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2unavspr7a_iout', iout)
       call ipgetr ('b2tfnb_vis_per', vis_per)
       if (in_parallel() .and. iout.ne.0) then
!$OMP MASTER
         write(*,*)
     .    'b2unavspr7a OpenMP warning: no file IO in parallel mode'
!$OMP END MASTER
         iout = 0
       endif
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

*   ..calculate the coefficients in difference approximation
*     of viscosity perpendicular current
      kmpr = nb/pb
      call b2txvsprco(nx, ny, isb, fac_vis, vol, hx, hy, qz,            !srv 02.04.18
     &    bb, zb, ti, kmpr, nb, avisper, bvisper, cvisper, dvisper)

      aa7 = 0.0_R8

      if(vis_per.ne.0.0_R8) then
*   ..compute the seven-point matrix                               !srv 26.01.00
      do iy = -1, ny
       do ix = -1, nx                                              !srv 26.11.02
         ireg = region(ix,iy,0)
         ixt = topix(ix,iy)
         iyt = topiy(ix,iy)
         ixb = bottomix(ix,iy)
         iyb = bottomiy(ix,iy)
         coef = 1.0_R8/qe  ! we divide by qe because the contribution is fchvisper/qe
         if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &      iyt.eq.ny+1 .or. iyb.eq.-2 .or. .not.solvereg(ireg)) then   !srv 26.11.02 12.12.17 {
*   .. do nothing
         else
           if(bottomiy(ixb,iyb).eq.-2) then
              aa7(ix,iy, 3) = coef*
     &         (avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)*
     &          dvisper(topix(ixt,iyt),topiy(ixt,iyt),1))
             aa7m2 =   coef*(
     &         avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*dvisper(ix ,iy ,1))
             aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7m2
             aa7p2 =     coef*(
     &         avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)*dvisper(ixt,iyt,1)
     &        +avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*dvisper(ixt,iyt,1))
             aa7(ix,iy,2) = aa7(ix,iy,2) - aa7(ix,iy,3) - aa7p2
             aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2
           elseif(topiy(ixt,iyt).eq.ny+1) then
              aa7(ix,iy,-3) = coef*
     &        (avisper(ix,iy,1)*bvisper(ixb,iyb,1)*dvisper(ixb,iyb,1))
             aa7m2 =   coef*(
     &        +avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*dvisper(ix ,iy ,1)
     &        +avisper(ix ,iy ,1)*bvisper(ixb,iyb,1)*dvisper(ix ,iy ,1))
             aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7(ix,iy,-3)  - aa7m2
             aa7p2 =     coef*(
     &        avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*dvisper(ixt,iyt,1))
             aa7(ix,iy,2) = aa7(ix,iy,2) - aa7p2
             aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2
           else
             aa7(ix,iy,-3) = coef*
     &        (avisper(ix,iy,1)*bvisper(ixb,iyb,1)*dvisper(ixb,iyb,1))
             aa7(ix,iy, 3) = coef*
     &         (avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)*
     &          dvisper(topix(ixt,iyt),topiy(ixt,iyt),1))
             aa7m2 =   coef*(
     &         avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*dvisper(ix ,iy ,1)
     &        +avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*dvisper(ix ,iy ,1)
     &        +avisper(ix ,iy ,1)*bvisper(ixb,iyb,1)*dvisper(ix ,iy ,1))
             aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7(ix,iy,-3)  - aa7m2
             aa7p2 =     coef*(
     &         avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)*dvisper(ixt,iyt,1)
     &        +avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*dvisper(ixt,iyt,1)
     &        +avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*dvisper(ixt,iyt,1))
             aa7(ix,iy,2) = aa7(ix,iy,2) - aa7(ix,iy,3) - aa7p2
             aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2
           endif
         endif                                                     !srv 26.11.02 12.12.17 }
       enddo
      enddo

      endif   !   vis_per.ne.0.0_R8

      if(iout.ne.0) then
       call my_out(70,nx,ny,aa7(-1,-1,  3),'b2unavspr7a_aa7p3')
       call my_out(70,nx,ny,aa7(-1,-1,  2),'b2unavspr7a_aa7p2')
       call my_out(70,nx,ny,aa7(-1,-1,  0),'b2unavspr7a_aa7m0')
       call my_out(70,nx,ny,aa7(-1,-1, -2),'b2unavspr7a_aa7m2')
       call my_out(70,nx,ny,aa7(-1,-1, -3),'b2unavspr7a_aa7m3')
      endif

      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2unavspr7a

      end subroutine b2unavspr7a                           !srv 02.04.18 }
