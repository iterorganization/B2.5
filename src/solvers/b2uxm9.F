*-----------------------------------------------------------------------
*.specification

      subroutine b2uxm9 (nx, ny, nregionv, itcnt, solvereg, conb,
     . aa9, name)                                                       !sv 13.10.06
      use b2mod_types
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv, itcnt
      real (kind=R8) ::
     *   conb(-1:nx,-1:ny,0:1,0:1,0:1)
      logical solvereg(0:nregionv)
      character*(*) name
*   ..input/output arguments
      real (kind=R8) ::
     *   aa9(-1:nx,-1:ny,-6:6)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UXM9 computes a matrix which corresponds of approximation mixed
*     derivatives in equations on nine-point stencil.
*
*-----------------------------------------------------------------------
!                 2_____3_____4                     5___________6
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -1_____0_____1                     |_____|_____|
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -4____-3____-2                    -6_____|____-5

!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
*.declarations

*   ..local variables
      integer ncall, ix, iy, ixl, iyl, ixr, iyr, ixb, iyb, ixt, iyt
      integer ireg
      integer, save :: iout = 0                                          !sv 13.10.06
      real (kind=R8) ::
     *   aaf(-1:nx,-1:ny,-2:2)
      save ncall
      !$OMP THREADPRIVATE(ncall, iout)
*   ..procedures
      external xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux5p
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uxm9')
*   ..set internal parameters on first call                              !sv 25.08.99
      if (ncall.eq.0) then                                               !sv 06.07.06 {
       call ipgeti ('b2uxm9_iout', iout)                                 !sv 13.10.06
       if (in_parallel() .and. iout.ne.0) then
         !$OMP MASTER
         write(*,*) 'b2uxm9 OpenMP warning: no file IO in parallel mode'
         !$OMP END MASTER
         iout = 0
       endif
      endif                                                              !sv 06.07.06 }
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of conb
       call b2xvff (nx, ny, conb(-1,-1,0,0,0), 'conb0')               !sv 13.10.06 {
       call b2xvff (nx, ny, conb(-1,-1,0,0,1), 'conb1')
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
            aaf(ix,iy,-1) = conb(ixl,iyl,0,0,0) -
     &                      conb(topix(ixl,iyl),topiy(ixl,iyl),0,0,1)
         endif
*       Case of right neighbor
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then

         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            aaf(ix,iy,1) = - conb(ixr,iyr,0,0,0) +
     &                       conb(topix(ixr,iyr),topiy(ixr,iyr),0,0,1)
         endif
        enddo
       enddo

       if(iout.ne.0) then
          call my_out(70,nx,ny,conb(-1,-1,0,0,0),'conbx3'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,0,0,1),'conbx4'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,1,1,0),'conby3'//trim(name))
          call my_out(70,nx,ny,conb(-1,-1,1,1,1),'conby4'//trim(name))

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
          if (rightix(ixb,iyb).ne.nx+1)
     >     aaf(ix,iy,-2) = -conb(ixb,iyb,1,1,0) +
     &      conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,1,1)
         endif
*       Case of top neighbor
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then

         else
          ixt = topix(ix,iy)
          iyt = topiy(ix,iy)
          aaf(ix,iy,2) =  conb(ixt,iyt,1,1,0) -
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1,1)
         endif
        enddo
       enddo

       if(iout.ne.0) then
         call my_out(70,nx,ny,aaf(-1,-1,-2),'aafm2'//trim(name))
         call my_out(70,nx,ny,aaf(-1,-1, 2),'aafp2'//trim(name))
       endif

       aa9(:,:,-3) = aa9(:,:,-3) + aaf(:,:,-2)
       aa9(:,:,-1) = aa9(:,:,-1) + aaf(:,:,-1)
       aa9(:,:, 1) = aa9(:,:, 1) + aaf(:,:, 1)
       aa9(:,:, 3) = aa9(:,:, 3) + aaf(:,:, 2)

*   ..compute additional P.Sch part of nine-point matrix.
       do iy = -1, ny
        do ix = -1, nx
         ireg=region(ix,iy,0)
*       Case (ix-1,iy-1)
         if(leftix(ix,iy).eq.-2 .or. bottomiy(ix,iy).eq.-2 .or.
     &      .not.solvereg(ireg)) then
            aa9(ix,iy,-4) = 0.0e0_R8
            aa9(ix,iy,-6) = 0.0e0_R8
         else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if(bottomix(ixl,iyl) .eq. leftix(ixb,iyb) .and.
     &         bottomiy(ixl,iyl) .eq. leftiy(ixb,iyb)) then
              aa9(ix,iy,-4) = aa9(ix,iy,-4)                           !sv 28.05.08
     &               - conb(bottomix(ixl,iyl),bottomiy(ixl,iyl),0,0,1)
     &               + conb(leftix(ixb,iyb),leftiy(ixb,iyb),1,1,0)
               aa9(ix,iy,-6) = 0.0e0_R8
            else if (leftix(ixb,iyb).ne.-2) then
               aa9(ix,iy,-4) = aa9(ix,iy,-4)                          !sv 28.05.08
     &              - conb(bottomix(ixl,iyl),bottomiy(ixl,iyl),0,0,1)
               aa9(ix,iy,-6) = aa9(ix,iy,-6) +                        !sv 28.05.08
     &                conb(leftix(ixb,iyb),leftiy(ixb,iyb),1,1,0)
            endif
         endif
*       Case (ix-1,iy+1)
         if(leftix(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then
            aa9(ix,iy,2) = 0.0e0_R8
            aa9(ix,iy,5) = 0.0e0_R8
         else
            ixl = leftix(ix,iy)
            iyl = leftiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if(leftix(ixt,iyt) .eq. topix(ixl,iyl) .and.
     &         leftiy(ixt,iyt) .eq. topiy(ixl,iyl)) then
              aa9(ix,iy,2) = aa9(ix,iy,2) +                           !sv 28.05.08
     &                    conb(topix(ixl,iyl),topiy(ixl,iyl),0,0,0) -
     &                    conb(leftix(ixt,iyt),leftiy(ixt,iyt),1,1,1)
             aa9(ix,iy,5) = 0.0e0_R8
            else
             aa9(ix,iy,2) = aa9(ix,iy,2) +                            !sv 28.05.08
     &                      conb(topix(ixl,iyl),topiy(ixl,iyl),0,0,0)
             aa9(ix,iy,5) = aa9(ix,iy,5) -                            !sv 28.05.08
     &                      conb(leftix(ixt,iyt),leftiy(ixt,iyt),1,1,1)
            endif
         endif
*       Case (ix+1,iy-1)
         if(rightix(ix,iy).eq.nx+1 .or. bottomiy(ix,iy).eq.-2 .or.
     &      .not.solvereg(ireg)) then
            aa9(ix,iy,-2) = 0.0e0_R8
            aa9(ix,iy,-5) = 0.0e0_R8
         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixb = bottomix(ix,iy)
            iyb = bottomiy(ix,iy)
            if(bottomix(ixr,iyr) .eq. rightix(ixb,iyb) .and.
     &         bottomiy(ixr,iyr) .eq. rightiy(ixb,iyb)) then
              aa9(ix,iy,-2) = aa9(ix,iy,-2) +                         !sv 28.05.08
     &                   conb(bottomix(ixr,iyr),bottomiy(ixr,iyr),0,0,1)
     &                 - conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,1,0)
              aa9(ix,iy,-5) = 0.0e0_R8
            else if (rightix(ixb,iyb).ne.nx+1) then
              aa9(ix,iy,-2) = aa9(ix,iy,-2) +                         !sv 28.05.08
     &                   conb(bottomix(ixr,iyr),bottomiy(ixr,iyr),0,0,1)
              aa9(ix,iy,-5) = aa9(ix,iy,-5) -                         !sv 28.05.08
     &                   conb(rightix(ixb,iyb),rightiy(ixb,iyb),1,1,0)
            endif
         endif
*       Case (ix+1,iy+1)
         if(rightix(ix,iy).eq.nx+1 .or. topiy(ix,iy).eq.ny+1 .or.
     &      .not.solvereg(ireg)) then
            aa9(ix,iy,4) = 0.0e0_R8
            aa9(ix,iy,6) = 0.0e0_R8
         else
            ixr = rightix(ix,iy)
            iyr = rightiy(ix,iy)
            ixt = topix(ix,iy)
            iyt = topiy(ix,iy)
            if(topix(ixr,iyr) .eq. rightix(ixt,iyt) .and.
     &         topiy(ixr,iyr) .eq. rightiy(ixt,iyt)) then
              aa9(ix,iy,4) = aa9(ix,iy,4) -                           !sv 28.05.08
     &                    conb(topix(ixr,iyr),topiy(ixr,iyr),0,0,0)+
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1,1)
              aa9(ix,iy,6) = 0.0e0_R8
            else
              aa9(ix,iy,4) = aa9(ix,iy,4) -                           !sv 28.05.08
     &                    conb(topix(ixr,iyr),topiy(ixr,iyr),0,0,0)
              aa9(ix,iy,6) = aa9(ix,iy,6) +                           !sv 28.05.08
     &                    conb(rightix(ixt,iyt),rightiy(ixt,iyt),1,1,1)
            endif
         endif

        enddo
       enddo

       if(iout.ne.0) then
          call my_out(70,nx,ny,aa9(-1,-1,-4),'aa9xm1ym1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1,-6),'aa9ym1xm1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1, 4),'aa9xp1yp1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1, 6),'aa9yp1xp1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1,-2),'aa9xp1ym1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1,-5),'aa9ym1xp1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1, 2),'aa9xm1yp1'//trim(name))
          call my_out(70,nx,ny,aa9(-1,-1, 5),'aa9yp1xm1'//trim(name))
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uxm9

      end subroutine b2uxm9
