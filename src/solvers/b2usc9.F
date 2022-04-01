*-----------------------------------------------------------------------
*.specification

      subroutine b2usc9 (nx, ny, isb, nregionv, solvereg, itcnt, rxg,   !srv 13.10.06 26.09.12
     &   pcm1, vol, pbs, nb, pb, snbv, flob, conb, pccb1, rescb, corpb,
     &   aa9, name)                                                     !srv 18.05.02
      use b2mod_types
      use b2mod_indirect
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt, isb, nregionv
      logical solvereg(0:nregionv)
      real (kind=R8) ::
     *  rxg, pcm1, pbs(-1:nx,-1:ny,0:1), nb(-1:nx,-1:ny),
     &  pb(-1:nx,-1:ny), snbv(-1:nx,-1:ny), flob(-1:nx,-1:ny,0:1,0:1),
     &  conb(-1:nx,-1:ny,0:1,0:1,0:4), pccb1(-1:nx,-1:ny),
     &  rescb(-1:nx,-1:ny), vol(-1:nx,-1:ny)
      character*(*) name
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpb(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa9(-1:nx,-1:ny,-6:6), wrk0(-1:nx,-1:ny)                        !srv 13.10.06

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USC9 computes a correction towards the solution of the
*     continuity equation for one atomic species using a
*     nine-point stencil.
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
      integer ncall, n2, ix, iy, ireg, ileft, ibot
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0             !srv 06.07.06
      integer, save :: iadd_9 = 1, iout = 0                             !srv 13.10.06
      !$OMP THREADPRIVATE(no_Ptncr_x_co, no_Ptncr_y_co, iadd_9, iout)
      real (kind=R8) ::
     *   t0, t1
      character chns*3
      save ncall
      !$OMP THREADPRIVATE(ncall)
*   ..procedures
      intrinsic abs
      external xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux9p, b2uxm9                   !srv 02.01.07
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usc9')
*   ..set internal parameters on first call                             !srv 25.08.99
      if (ncall.eq.0) then                                              !srv 06.07.06 {
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)
       call ipgeti ('b2usc9_add_9', iadd_9)                             !srv 13.10.06
       call ipgeti ('b2usc9_iout', iout)                                !srv 13.10.06
       if (in_parallel() .and. iout.ne.0) then
         !$OMP MASTER
         write(*,*) 'b2usc9 OpenMP warning: no file IO in parallel mode'
         !$OMP END MASTER
         iout = 0
       endif
      endif                                                             !srv 06.07.06 }
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt, rxg, pcm1
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.lt.rxg, 'faulty argument rxg')
      call xertst (0.le.pcm1, 'faulty argument pcm1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of pbs
       call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbsx')
       call b2xvfy (nx, ny, pbs(-1,-1,1), 'pbsy')
*    ..test sign of nb, pb
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
*    ..test sign of snbv, conb
       call b2xvsg (n2, snbv, 1, 'snbv', '.le.')
       if(no_Ptncr_x_co.eq.0 .and. no_Ptncr_y_co.eq.0) then             !srv 06.07.06 {
        call b2xvsg (n2, conb(-1,-1,0,0,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,2), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,2), 1, 'conb', '.ge.')
       endif                                                            !srv 06.07.06 }
*    ..test edge values of flob, conb
       call b2xvff (nx, ny, flob, 'flob')
       call b2xvff (nx, ny, conb(-1,-1,0,0,0), 'conb0')                 !srv 13.10.06 {
       call b2xvff (nx, ny, conb(-1,-1,0,0,1), 'conb1')
       call b2xvff (nx, ny, conb(-1,-1,0,0,2), 'conb2')                 !srv 13.10.06 }
!*    ..test edge values of pccb1                                       !srv 15.04.05
!       call b2xvfx (nx, ny, pccb1, 'pccb1')                            !srv 15.04.05
      endif

      aa9 = 0.0e0_R8

*   ..compute wrk0
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
            wrk0(ix,iy) = 0.0_R8
          else
*!!     Maybe use volume weighting or inverse-volume weighting here.
            wrk0(ix,iy) = pbs(ix,iy,0)*
     &       (pccb1(leftix(ix,iy),leftiy(ix,iy))*
     &        nb(leftix(ix,iy),leftiy(ix,iy))+
     &        pccb1(ix,iy)*nb(ix,iy))/2.0_R8
          endif
        enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)

* ..compute the nine-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
          aa9(ix,iy,0) = (-snbv(ix,iy)*nb(ix,iy)/pb(ix,iy)+             !srv 13.10.06 {
     &     abs(rescb(ix,iy))/(rxg*pb(ix,iy)))*rscal(pb(ix,iy))          !srv 13.10.06 }
        else
          aa9(ix,iy,0) = rscal(pb(ix,iy))
        endif
       enddo
      enddo

      do iy = -1, ny
       do ix = -1, nx                                                   !srv 26.11.02 {
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
          aa9(ix,iy,-1) = 0.0e0_R8
        else                                                            !srv 26.11.02 }
          t0 =  flob(ix,iy,0,0)/2.0_R8
     &        - (conb(ix,iy,0,0,0)+pcm1*wrk0(ix,iy))
     &        - conb(ix,iy,0,0,1)                                       !srv 13.10.06 
          t1 = -flob(ix,iy,0,0)/2.0_R8
     &        - (conb(ix,iy,0,0,0)+pcm1*wrk0(ix,iy))
     &        - conb(ix,iy,0,0,2)                                       !srv 13.10.06
          ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
          if (solvereg(ileft)) then
           aa9(leftix(ix,iy),leftiy(ix,iy),1) =
     &       aa9(leftix(ix,iy),leftiy(ix,iy),1) + t0*rscal(pb(ix,iy))
           aa9(leftix(ix,iy),leftiy(ix,iy),0) =
     &       aa9(leftix(ix,iy),leftiy(ix,iy),0) -
     &       t1*rscal(pb(leftix(ix,iy),leftiy(ix,iy)))
          endif
          if (solvereg(ireg)) then
           aa9(ix,iy,0)  = aa9(ix,iy,0) - t0*rscal(pb(ix,iy))
           aa9(ix,iy,-1) = aa9(ix,iy,-1) +
     &       t1*rscal(pb(leftix(ix,iy),leftiy(ix,iy)))
          endif
cwdk add 9-point contributions
cwdk todo: rscal(pb)-logic
         if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
          t0  = (0.5_R8*flob(ix,iy,0,1)-conb(ix,iy,0,1,0)
     &                                 -conb(ix,iy,0,1,1))/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          t1  = (-0.5_R8*flob(ix,iy,0,1)-conb(ix,iy,0,1,0)
     &                                  -conb(ix,iy,0,1,2))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          if (solvereg(ireg)) then
           aa9(ix,iy,-4) = aa9(ix,iy,-4)           + t1/
     &       vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &           bottomiy(leftix(ix,iy),leftiy(ix,iy)))
           aa9(ix,iy,-3) = aa9(ix,iy,-3)           + t1/
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
           aa9(ix,iy,-1) = aa9(ix,iy,-1) - (t0 - t1)/
     &       vol(leftix(ix,iy),leftiy(ix,iy))
           aa9(ix,iy, 0) = aa9(ix,iy, 0) - (t0 - t1)/
     &       vol(ix,iy)
           aa9(ix,iy, 2) = aa9(ix,iy, 2) -  t0/
     &       vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &           leftiy(topix(ix,iy),topiy(ix,iy)))
           aa9(ix,iy, 3) = aa9(ix,iy, 3) -  t0/
     &       vol(topix(ix,iy),topiy(ix,iy))
          endif
          if (solvereg(ileft)) then
           aa9(leftix(ix,iy),leftiy(ix,iy),-3) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy),-3)       - t1/
     &       vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &           bottomiy(leftix(ix,iy),leftiy(ix,iy)))
           aa9(leftix(ix,iy),leftiy(ix,iy),-2) =
     &       aa9(leftix(ix,iy),leftiy(ix,iy),-2)       - t1/
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
           aa9(leftix(ix,iy),leftiy(ix,iy), 0) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 0) + (t0 - t1)/
     &       vol(leftix(ix,iy),leftiy(ix,iy))
           aa9(leftix(ix,iy),leftiy(ix,iy), 1) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 1) + (t0 - t1)/
     &       vol(ix,iy)
           aa9(leftix(ix,iy),leftiy(ix,iy), 3) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 3) + t0/
     &       vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &           leftiy(topix(ix,iy),topiy(ix,iy)))
           aa9(leftix(ix,iy),leftiy(ix,iy), 4) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 4) + t0 /
     &       vol(topix(ix,iy),topiy(ix,iy))
          endif
         endif
        endif                                                           !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg))
     >   aa9(ix,iy,1) = 0.0e0_R8                                        !srv 26.11.02
       enddo
      enddo

      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = -1, ny
       do ix = -1, nx
         ireg = region(ix,iy,0)
         if(bottomiy(ix,iy).eq.-2) then
          aa9(ix,iy,-3) = 0.0_R8
         else
          t0 = 0.5_R8*flob(ix,iy,1,1)
     &       - conb(ix,iy,1,1,0) - conb(ix,iy,1,1,1)                    !srv 13.10.06
          t1 =-0.5_R8*flob(ix,iy,1,1)
     &       - conb(ix,iy,1,1,0) - conb(ix,iy,1,1,2)                    !srv 13.10.06
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot)) then
           aa9(bottomix(ix,iy),bottomiy(ix,iy),3) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),3) + t0*rscal(pb(ix,iy))
           aa9(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),0)
     &      -t1*rscal(pb(bottomix(ix,iy),bottomiy(ix,iy)))
          endif
          if (solvereg(ireg)) then
           aa9(ix,iy,0)  = aa9(ix,iy,0) - t0*rscal(pb(ix,iy))
           aa9(ix,iy,-3) = aa9(ix,iy,-3) +
     &      t1*rscal(pb(bottomix(ix,iy),bottomiy(ix,iy)))
          endif
cwdk add 9-point contributions
cwdk todo: add rscal-logic
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           t0 = (0.5_R8*flob(ix,iy,1,0)-conb(ix,iy,1,0,0)
     &                                 -conb(ix,iy,1,0,1))/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
           t1 = (-0.5_R8*flob(ix,iy,1,0)-conb(ix,iy,1,0,0)
     &                                  -conb(ix,iy,1,0,2))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           if (solvereg(ireg)) then
            aa9(ix,iy,-4) = aa9(ix,iy,-4)       + t1/
     &        vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))
            aa9(ix,iy,-3) = aa9(ix,iy,-3) - (t0 - t1)/
     &        vol(bottomix(ix,iy),bottomiy(ix,iy))
            aa9(ix,iy,-2) = aa9(ix,iy,-2) -  t0/
     &        vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))
            aa9(ix,iy,-1) = aa9(ix,iy,-1)       + t1/
     &        vol(leftix(ix,iy),leftiy(ix,iy))
            aa9(ix,iy, 0) = aa9(ix,iy, 0) - (t0 - t1)/
     &        vol(ix,iy)
            aa9(ix,iy, 1) = aa9(ix,iy, 1) -  t0/
     &        vol(rightix(ix,iy),rightiy(ix,iy))
           endif
           if (solvereg(ibot)) then
            aa9(bottomix(ix,iy),bottomiy(ix,iy),-1) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy),-1)       - t1/
     &        vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 0) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 0) + (t0 - t1)/
     &        vol(bottomix(ix,iy),bottomiy(ix,iy)) 
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 1) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 1) +  t0/
     &        vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 2) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 2)      - t1/
     &        vol(leftix(ix,iy),leftiy(ix,iy))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 3) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 3) + (t0 - t1)/
     &        vol(ix,iy)
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 4) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 4) +  t0/
     &        vol(rightix(ix,iy),rightiy(ix,iy))
           endif
          endif
         endif
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,ny,3), 1)

*   ..add matrix from contribution nine-point stencil
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nnreg(0), itcnt, solvereg,
     &   conb(-1,-1,0,0,3), aa9, name)                                  !srv 02.01.07
      endif

      do iy=-1,ny
       do ix=-1,nx
        ireg = region(ix,iy,0)
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
          aa9(ix,iy,0) = rscal(pb(ix,iy))
          rescb(ix,iy) = 0.0_R8
        end if
        if (.not.solvereg(ireg)) rescb(ix,iy) = 0.0_R8
       enddo
      enddo
*   ..temporary fix to stabilize 9-point stencil in corners:
*   ..set guard cell value equal to neighboring cell
      aa9(-1,-1,:) = 0.0_R8
      aa9(-1,-1,0) = 1.00_R8*nb(-1,-1)/pb(-1,-1)
!      aa9(-1,-1,3) =-1.00_R8
!      rescb(-1,-1) = pb(-1,0) - pb(-1,-1)
      aa9(-1,-1,1) =-1.00_R8*nb(0,-1)/pb(0,-1)
      rescb(-1,-1) = nb(0,-1) - nb(-1,-1)

      aa9(-1,ny,:)  = 0.0_R8
      aa9(-1,ny,0)  = 1.00_R8*nb(-1,ny)/pb(-1,ny)
!      aa9(-1,ny, 1) =-0.50_R8
!      aa9(-1,ny,-3) =-0.50_R8
!      rescb(-1,ny)  = 0.5_R8*(pb(0,ny) + pb(-1,ny-1)) - pb(-1,ny)
!      aa9(-1,ny,-3) =-1.00_R8
!      rescb(-1,ny)  = pb(-1,ny-1) - pb(-1,ny)
      aa9(-1,ny,1) =-1.00_R8*nb(0,ny)/pb(0,ny)
      rescb(-1,ny) = nb(0,ny) - nb(-1,ny)

      aa9(nx,-1,:)  = 0.0_R8
      aa9(nx,-1,0)  = 1.00_R8*nb(nx,-1)/pb(nx,-1)
!      aa9(nx,-1,3)  =-1.00_R8 
!      rescb(nx,-1)  = pb(nx,0) - pb(nx,-1)
      aa9(nx,-1,-1) =-1.00_R8*nb(nx-1,-1)/pb(nx-1,-1) 
      rescb(nx,-1)  = nb(nx-1,-1) - nb(nx,-1)

      aa9(nx,ny,:)  = 0.0_R8
      aa9(nx,ny,0)  = 1.00_R8*nb(nx,ny)/pb(nx,ny)
!      aa9(nx,ny,-3) =-1.00_R8 
!      rescb(nx,ny)  = pb(nx,ny-1) - pb(nx,ny)
      aa9(nx,ny,-1) =-1.00_R8*nb(nx-1,ny)/pb(nx-1,ny) 
      rescb(nx,ny)  = nb(nx-1,ny) - nb(nx,ny)

*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, rescb, corpb, name)              !srv 18.05.02
*   ..redo the scaling
      do iy = -1, ny
       do ix = -1, nx
        corpb(ix,iy) = corpb(ix,iy)*rscal(pb(ix,iy))
       enddo
      enddo
      if (iout.ne.0) then                                               !srv 17.06.02 {
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,aa9(-1,-1,-6),name//'_aa9m6'//chns)
        call my_out(70,nx,ny,aa9(-1,-1,-5),name//'_aa9m5'//chns)
        call my_out(70,nx,ny,aa9(-1,-1,-4),name//'_aa9m4'//chns)
        call my_out(70,nx,ny,aa9(-1,-1,-3),name//'_aa9m3'//chns)
        call my_out(70,nx,ny,aa9(-1,-1,-2),name//'_aa9m2'//chns)
        call my_out(70,nx,ny,aa9(-1,-1,-1),name//'_aa9m1'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 0),name//'_aa9cc'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 1),name//'_aa9p1'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 2),name//'_aa9p2'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 3),name//'_aa9p3'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 4),name//'_aa9p4'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 5),name//'_aa9p5'//chns)
        call my_out(70,nx,ny,aa9(-1,-1, 6),name//'_aa9p6'//chns)
        call my_out(70,nx,ny,corpb,name//'_corpb'//chns)
        call my_out(70,nx,ny,rescb,name//'_rescb'//chns)
        call my_out(70,nx,ny,nb,name//'_na'//chns)
      endif                                                             !srv 17.06.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usc9

      contains

      function rscal(t0)
      implicit none
      real(kind=R8) :: rscal, t0
!      rscal = t0
      rscal = 1.0e0_R8                                                  !srv 13.10.06
      return
      end function rscal

      end subroutine b2usc9
