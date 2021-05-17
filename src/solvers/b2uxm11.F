*-----------------------------------------------------------------------
*.specification

      subroutine b2uxm11 (nx, ny, isb, zb, fac_vis, nregionv, itcnt,
     &                    solvereg, po, ti, nb, pb, conb, aa11, name)    !srv 02.04.18
      use b2mod_types
      use b2mod_indirect
      use b2mod_openmp
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, nregionv, itcnt
      real (kind=R8) ::
     *   fac_vis(-1:nx,-1:ny), zb(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &   ti(-1:nx,-1:ny), nb(-1:nx,-1:ny), pb(-1:nx,-1:ny),
     &   conb(-1:nx,-1:ny,0:1,0:1)
      logical solvereg(0:nregionv)
      character*(*) name
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   aa11(-1:nx,-1:ny,-7:7)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UXM11 computes a matrix which corresponds of approximation mixed
*     derivatives in equations on 11-point stencil.
*
*-----------------------------------------------------------------------
!                                                           7
!                                                           |
!                                                           |
!       2_____3_____4          5___________6           _____3_____
!       |     |     |          |     |     |          |     |     |
!       |     |     |          |     |     |          |     |     |
!      -1_____0_____1          |_____|_____|          |_____0_____|
!       |     |     |          |     |     |          |     |     |
!       |     |     |          |     |     |          |     |     |
!      -4____-3____-2         -6_____|____-5          |____-3_____|
!                                                           |
!                                                           |
!                                                          -7
!
!                                 (ix,iy+2)
!                                    |
!                                    |
!                                    |
!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
!                                    |
!                                    |
!                                    |
!                                (ix,iy-2)
!
*.declarations

*   ..local variables
      integer ncall, ix, iy, ixl, iyl, ixr, iyr, ixb, iyb, ixt, iyt
      integer ireg
      integer, save :: iout = 0                                          !srv 13.10.06
      real (kind=R8) ::
     *   aaf(-1:nx,-1:ny,-2:2), aa7(-1:nx,-1:ny,-3:3)
      save ncall
      !$OMP THREADPRIVATE(ncall, iout)
*   ..procedures
      external subini, subend, xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux5p
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxm11')
*   ..set internal parameters on first call                              !srv 25.08.99
      if (ncall.eq.0) then                                               !srv 06.07.06 {
       call ipgeti ('b2uxm11_iout', iout)                                 !srv 13.10.06
       if (in_parallel() .and. iout.ne.0) then
         !$OMP MASTER
         write(*,*)
     .    'b2uxm11 OpenMP warning: no file IO in parallel mode'
         !$OMP END MASTER
         iout = 0
       endif
      endif                                                              !srv 06.07.06 }
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of conb
       call b2xvff (nx, ny, conb(-1,-1,0,0), 'conb0')                    !srv 13.10.06 {
       call b2xvff (nx, ny, conb(-1,-1,0,1), 'conb1')
      endif

*   ..compute contribution in x-part of five-point stencil
       aaf = 0.0e0_R8
*     from P.Sch x-part of nine-point stencil
       do iy = -1, ny
        do ix = -1, nx
         ireg=region(ix,iy,0)
*       Case of left neighbor
         if(leftix(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then

         else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            aaf(ix,iy,-1) = conb(ixl,iyl,0,0) -
     &                      conb(topix(ixl,iyl),topiy(ixl,iyl),0,1)
         endif
*       Case of right neighbor
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then

         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            aaf(ix,iy,1) = - conb(ixr,iyr,0,0) +
     &                       conb(topix(ixr,iyr),topiy(ixr,iyr),0,1)
         endif
        enddo
       enddo

       if(iout.ne.0) then
          call my_out(70,nx,ny,conb(-1,-1,0,0),'conbx3'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,0,1),'conbx4'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,1,0),'conby3'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,1,1),'conby4'//trim(name))

          call my_out(70,nx,ny,aaf(-1,-1,-1),'aafm1'//trim(name))
          call my_out(70,nx,ny,aaf(-1,-1, 1),'aafp1'//trim(name))
          endif

*     from P.Sch y-part of nine-point stencil
       do iy = -1, ny
        do ix = -1, nx
         ireg=region(ix,iy,0)
*       Case of bottom neighbor
         if(bottomiy(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &      .not.solvereg(ireg)) then

         else
          ixb = bottomix(ix,iy)
          iyb = bottomiy(ix,iy)
          aaf(ix,iy,-2) = -conb(ixb,iyb,1,0) +
     &                     conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,1)
         endif
*       Case of top neighbor
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then

         else
          ixt = topix(ix,iy)
          iyt = topiy(ix,iy)
          aaf(ix,iy,2) =  conb(ixt,iyt,1,0) -
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1)
         endif
        enddo
       enddo

       if(iout.ne.0) then
         call my_out(70,nx,ny,aaf(-1,-1,-2),'aafm2'//trim(name))
         call my_out(70,nx,ny,aaf(-1,-1, 2),'aafp2'//trim(name))
       endif

       aa11(:,:,-3) = aa11(:,:,-3) + aaf(:,:,-2)
       aa11(:,:,-1) = aa11(:,:,-1) + aaf(:,:,-1)
       aa11(:,:, 1) = aa11(:,:, 1) + aaf(:,:, 1)
       aa11(:,:, 3) = aa11(:,:, 3) + aaf(:,:, 2)

*   ..compute additional P.Sch part of nine-point matrix.
       do iy = -1, ny
        do ix = -1, nx
         ireg=region(ix,iy,0)
*       Case (ix-1,iy-1)
         if(leftix(ix,iy).eq.-2 .or. bottomiy(ix,iy).eq.-2 .or.
     &      .not.solvereg(ireg)) then
            aa11(ix,iy,-4) = 0.0e0_R8
            aa11(ix,iy,-6) = 0.0e0_R8
         else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if(bottomix(ixl,iyl) .eq. leftix(ixb,iyb) .and.
     &         bottomiy(ixl,iyl) .eq. leftiy(ixb,iyb)) then
              aa11(ix,iy,-4) = aa11(ix,iy,-4)                           !srv 28.05.08
     &               - conb(bottomix(ixl,iyl),bottomiy(ixl,iyl),0,1)
     &               + conb(leftix(ixb,iyb),leftiy(ixb,iyb),1,0)
               aa11(ix,iy,-6) = 0.0e0_R8
            else
               aa11(ix,iy,-4) = aa11(ix,iy,-4)                          !srv 28.05.08
     &              - conb(bottomix(ixl,iyl),bottomiy(ixl,iyl),0,1)
               aa11(ix,iy,-6) = aa11(ix,iy,-6) +                        !srv 28.05.08
     &                conb(leftix(ixb,iyb),leftiy(ixb,iyb),1,0)
            endif
         endif
*       Case (ix-1,iy+1)
         if(leftix(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then
          aa11(ix,iy,2) = 0.0e0_R8
          aa11(ix,iy,5) = 0.0e0_R8
         else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if(leftix(ixt,iyt) .eq. topix(ixl,iyl) .and.
     &         leftiy(ixt,iyt) .eq. topiy(ixl,iyl)) then
              aa11(ix,iy,2) = aa11(ix,iy,2) +                           !srv 28.05.08
     &                    conb(topix(ixl,iyl),topiy(ixl,iyl),0,0) -
     &                    conb(leftix(ixt,iyt),leftiy(ixt,iyt),1,1)
             aa11(ix,iy,5) = 0.0e0_R8
            else
             aa11(ix,iy,2) = aa11(ix,iy,2) +                            !srv 28.05.08
     &                      conb(topix(ixl,iyl),topiy(ixl,iyl),0,0)
             aa11(ix,iy,5) = aa11(ix,iy,5) -                            !srv 28.05.08
     &                      conb(leftix(ixt,iyt),leftiy(ixt,iyt),1,1)
            endif
         endif
*       Case (ix+1,iy-1)
         if(rightix(ix,iy).eq.nx+1 .or. bottomiy(ix,iy).eq.-2 .or.
     &      .not.solvereg(ireg)) then
          aa11(ix,iy,-2) = 0.0e0_R8
          aa11(ix,iy,-5) = 0.0e0_R8
         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if(bottomix(ixr,iyr) .eq. rightix(ixb,iyb) .and.
     &       bottomiy(ixr,iyr) .eq. rightiy(ixb,iyb)) then
              aa11(ix,iy,-2) = aa11(ix,iy,-2) +                         !srv 28.05.08
     &                   conb(bottomix(ixr,iyr),bottomiy(ixr,iyr),0,1)
     &                 - conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,0)
              aa11(ix,iy,-5) = 0.0e0_R8
            else
              aa11(ix,iy,-2) = aa11(ix,iy,-2) +                         !srv 28.05.08
     &                   conb(bottomix(ixr,iyr),bottomiy(ixr,iyr),0,1)
              aa11(ix,iy,-5) = aa11(ix,iy,-5) -                         !srv 28.05.08
     &                   conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,0)
            endif
         endif
*       Case (ix+1,iy+1)
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then
           aa11(ix,iy,4) = 0.0e0_R8
           aa11(ix,iy,6) = 0.0e0_R8
         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if(topix(ixr,iyr) .eq. rightix(ixt,iyt) .and.
     &       topiy(ixr,iyr) .eq. rightiy(ixt,iyt)) then
              aa11(ix,iy,4) = aa11(ix,iy,4) -                           !srv 28.05.08
     &                    conb(topix(ixr,iyr),topiy(ixr,iyr),0,0)+
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1)
              aa11(ix,iy,6) = 0.0e0_R8
            else
              aa11(ix,iy,4) = aa11(ix,iy,4) -                           !srv 28.05.08
     &                    conb(topix(ixr,iyr),topiy(ixr,iyr),0,0)
              aa11(ix,iy,6) = aa11(ix,iy,6) +                           !srv 28.05.08
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1)
            endif
         endif

        enddo
      enddo
*  .. calculate the coefficients of 7-point stencil from
*     viscosity perpendicular current by numerically
!      call b2unavspr7n(nx, ny, isb, zb, fac_vis, po, ti, nb, pb, aa7)  !srv 02.04.18

      call b2unavspr7a(nx, ny, nregionv, solvereg, isb, zb,              !srv 02.04.18 {
     &              fac_vis, ti, nb, pb, aa7)
*  .. calculate the coefficients of 7-point stencil from
*     viscosity perpendicular current by analytically

*  .. add coefficients
      aa11(:,:,-7) = aa11(:,:,-7) + aa7(:,:,-3)
      aa11(:,:,-3) = aa11(:,:,-3) + aa7(:,:,-2)
      aa11(:,:, 0) = aa11(:,:, 0) + aa7(:,:, 0)
      aa11(:,:, 3) = aa11(:,:, 3) + aa7(:,:, 2)
      aa11(:,:, 7) = aa11(:,:, 7) + aa7(:,:, 3)

      if(iout.ne.0) then
        call my_out(70,nx,ny,aa11(-1,-1,-7),'aa11xym2'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-6),'aa11ym1xm1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-5),'aa11ym1xp1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-4),'aa11xm1ym1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-3),'aa11xym1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-2),'aa11xp1ym1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1,-1),'aa11xm1y'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 0),'aa11xy'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 1),'aa11xp1y'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 2),'aa11xm1yp1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 3),'aa11xyp1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 4),'aa11xp1yp1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 5),'aa11yp1xm1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 6),'aa11yp1xp1'//trim(name))
        call my_out(70,nx,ny,aa11(-1,-1, 7),'aa11xyp2'//trim(name))
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uxm11

      end subroutine b2uxm11                            !srv 02.04.18 }
