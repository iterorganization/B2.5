*-----------------------------------------------------------------------
*.specification

      subroutine b2ush9 (nx, ny, ns, nregionv, itcnt,
     &  rxg, pcm1, vol, pbs, ne, ni, te, ti,
     &  floe, cone, floi, coni, pccm1, she, shi, ceqp,
     &  resht, reshe, reshi, cortt, corte, corti, flo0, con0, aa9)
      use b2mod_types
      use b2mod_indirect
      use b2mod_ilutern_9
      use b2mod_numerics_namelist
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv, itcnt
      real (kind=R8) ::
     *  rxg, pcm1, pbs(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  ni(-1:nx,-1:ny,0:1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  floe(-1:nx,-1:ny,0:1,0:1), floi(-1:nx,-1:ny,0:1,0:1),
     &  cone(-1:nx,-1:ny,0:1,0:1,0:2), coni(-1:nx,-1:ny,0:1,0:1,0:2),   !srv 02.01.07
     &  pccm1(-1:nx,-1:ny), vol(-1:nx,-1:ny)
      real (kind=R8) ::
     *  she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3),
     &  ceqp(-1:nx,-1:ny), resht(-1:nx,-1:ny), reshe(-1:nx,-1:ny),
     &  reshi(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cortt(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(-1:nx,-1:ny,0:1,0:1), con0(-1:nx,-1:ny,0:1,0:1,0:2),       !srv 02.01.07
     &  aa9(-1:nx,-1:ny,-6:6)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USH9 computes a correction to the fields te and ti.
*
*-----------------------------------------------------------------------
!                 2_____3_____4                     5___________6
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -1_____0_____1                     |_____|_____|
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -4____-3____-2                    -6_____|____-5
!
!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0           !srv 06.07.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0           !srv 06.07.06
      integer, save :: iadd_9 = 1                                     !srv 02.01.07
      integer, save :: style_b2mndt = 1                               !srv 07.07.21
      integer, save :: iout = 0
      integer, save :: use_9pt_stencil = 1
      real (kind=R8), save :: cor_diag = 1.0_R8                       !srv 07.07.21
      !$OMP THREADPRIVATE(no_Ptncr_x_he, no_Ptncr_y_he, iadd_9, iout)
      !$OMP THREADPRIVATE(no_Ptncr_x_hi, no_Ptncr_y_hi, style_b2mndt)
      !$OMP THREADPRIVATE(use_9pt_stencil, cor_diag)
      integer ireg, ileft, ibot
      save ncall
      !$OMP THREADPRIVATE(ncall)
      real (kind=R8) ::
     *   t0, t1
      logical new_matrix
*   ..procedures
      intrinsic abs, sqrt, dim
      external xertst, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux9p, b2uxm9                 !srv 02.01.07
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ush9')
      if (ncall.eq.0) then
       call ipgeti ('b2mndr_use_9pt_stencil', use_9pt_stencil)
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)               !srv 06.07.06 {
       call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)
       call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)               !srv 06.07.06 }
       call ipgeti ('b2ush9_add_9', iadd_9)                           !srv 02.01.07
       call ipgetr ('b2ush9_cor_diag', cor_diag)                      !srv 07.07.21
       call ipgeti ('b2mndt_style', style_b2mndt)                     !srv 07.07.21
       call ipgeti ('b2ush9_iout', iout)
       if (in_parallel() .and. iout.ne.0) then
         !$OMP MASTER
         write(*,*) 'b2ush9 OpenMP warning: no file IO in parallel mode'
         !$OMP END MASTER
         iout = 0
       endif
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt, rxg, pcm1
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.lt.rxg, 'faulty argument rxg')
      call xertst (0.le.pcm1, 'faulty argument pcm1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test edge values of pbs
       call b2xvfx (nx, ny, pbs(-1,-1,0), 'pbsx')
       call b2xvfy (nx, ny, pbs(-1,-1,1), 'pbsy')
*    ..test sign of ne, ni, te, ti
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test sign of cone, coni
       if(no_Ptncr_x_he.eq.0 .and. no_Ptncr_y_he.eq.0) then           !srv 06.07.06 {
        call b2xvsg (n2, cone(-1,-1,0,0,0), 1, 'cone', '.ge.')        !srv 02.01.07
        call b2xvsg (n2, cone(-1,-1,1,1,0), 1, 'cone', '.ge.')
       endif                                                          !srv 06.07.06 }
       if(no_Ptncr_x_hi.eq.0 .and. no_Ptncr_y_hi.eq.0) then           !srv 06.07.06 {
        call b2xvsg (n2, coni(-1,-1,0,0,0), 1, 'coni', '.ge.')
        call b2xvsg (n2, coni(-1,-1,1,1,0), 1, 'coni', '.ge.')
       endif
*    ..test edge values of floe, floi, cone, coni
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0,0), 'cone0')               !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,0,2), 'coni2')               !srv 02.01.07 }
*    ..test sign of she, shi, ceqp
       if(style_b2mndt.ne.2) then                                     !srv 07.07.21
         call b2xvsg (n2, she(-1,-1,0), 1, 'she0', '.ge.')
       endif                                                          !srv 07.07.21
       call b2xvsg (n2, she(-1,-1,1), 1, 'she1', '.le.')
       call b2xvsg (n2, she(-1,-1,2), 1, 'she2', '.ge.')
       call b2xvsg (n2, she(-1,-1,3), 1, 'she3', '.le.')
       if(style_b2mndt.ne.2) then                                     !srv 07.07.21
         call b2xvsg (n2, shi(-1,-1,0), 1, 'shi0', '.ge.')
       endif                                                          !srv 07.07.21
       call b2xvsg (n2, shi(-1,-1,1), 1, 'shi1', '.le.')
       call b2xvsg (n2, shi(-1,-1,2), 1, 'shi2', '.ge.')
       call b2xvsg (n2, shi(-1,-1,3), 1, 'shi3', '.le.')
       call b2xvsg (n2, ceqp, 1, 'ceqp', '.ge.')
      endif

      if (any(solveee(0:nnreg(0))).and.any(solveei(0:nnreg(0))).and.
     &    any(solveet(0:nnreg(0)))) then                                 !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveet(ireg).neqv.last_solve_9(ireg)
      end do
      if (new_matrix.or.use_9pt_stencil.ne.0) call restart_ma28_for_9
* ..total heat correction equation
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0,0:1)   = 0.0e0_R8
         con0(ix,iy,0,0:1,:) = 0.0e0_R8                               !srv 02.01.07
        else
         if (solveee(ireg) .and. solveei(ireg)) then
          t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*                              !srv 02.01.07 {
     &     (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &     (ne(leftix(ix,iy),leftiy(ix,iy))*
     &      te(leftix(ix,iy),leftiy(ix,iy))+
     &      ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &      ti(leftix(ix,iy),leftiy(ix,iy)))+
     &      pccm1(ix,iy)*
     &     (ne(ix,iy)*te(ix,iy)+
     &       ni(ix,iy,0)*ti(ix,iy)))/2.0_R8                           !srv 02.01.07 }
          flo0(ix,iy,0,0) = floe(ix,iy,0,0)+floi(ix,iy,0,0)
     &     -t0*(ne(ix,iy)+ni(ix,iy,0)-ne(leftix(ix,iy),leftiy(ix,iy)) !srv 02.01.07 {
     &             -ni(leftix(ix,iy),leftiy(ix,iy),0))                !srv 02.01.07 }
          con0(ix,iy,0,0,0) = cone(ix,iy,0,0,0)+coni(ix,iy,0,0,0)     !srv 02.01.07
     &     +(t0/2.0_R8)*(ne(leftix(ix,iy),leftiy(ix,iy))              !srv 02.01.07 {
     &      +ne(ix,iy)+ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0)) !srv 02.01.07 }
cwdk add linearization for 9-point terms
          flo0(ix,iy,0,1)   = floe(ix,iy,0,1)   + floi(ix,iy,0,1)
          con0(ix,iy,0,1,0) = cone(ix,iy,0,1,0) + coni(ix,iy,0,1,0)
         else if (.not.solveee(ireg) .and. solveei(ireg)) then
          t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*                              !srv 02.01.07 {
     &     (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &      ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &      ti(leftix(ix,iy),leftiy(ix,iy))+
     &      pccm1(ix,iy)*ni(ix,iy,0)*ti(ix,iy))/2.0_R8                !srv 02.01.07 }
          flo0(ix,iy,0,0) = floi(ix,iy,0,0)
     &     -t0*(ni(ix,iy,0)-ni(leftix(ix,iy),leftiy(ix,iy),0))        !srv 02.01.07 }
          con0(ix,iy,0,0,0) = coni(ix,iy,0,0,0)                       !srv 02.01.07
     &     +(t0/2.0_R8)*(ni(leftix(ix,iy),leftiy(ix,iy),0)+
     &                   ni(ix,iy,0))                                 !srv 02.01.07 }
cwdk add linearization for 9-point terms
          flo0(ix,iy,0,1)   = floi(ix,iy,0,1)
          con0(ix,iy,0,1,0) = coni(ix,iy,0,1,0)
         else if (solveee(ireg) .and. .not.solveei(ireg)) then
          t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*                              !srv 02.01.07 {
     &     (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &      ne(leftix(ix,iy),leftiy(ix,iy))*
     &      te(leftix(ix,iy),leftiy(ix,iy))+
     &      pccm1(ix,iy)*ne(ix,iy)*te(ix,iy))/2.0_R8                  !srv 02.01.07 }
          flo0(ix,iy,0,0) = floe(ix,iy,0,0)
     &     -t0*(ne(ix,iy)-ne(leftix(ix,iy),leftiy(ix,iy)))            !srv 02.01.07 }
          con0(ix,iy,0,0,0) = cone(ix,iy,0,0,0)                       !srv 02.01.07
     &     +(t0/2.0_R8)*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))   !srv 02.01.07 }
cwdk add linearization for 9-point terms
          flo0(ix,iy,0,1)   = floe(ix,iy,0,1)
          con0(ix,iy,0,1,0) = cone(ix,iy,0,1,0)
         endif
        endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1,0:1) = 0.0e0_R8
       con0(ix,-1,1,0:1,:) = 0.0e0_R8                                 !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveee(ireg) .and. solveei(ireg)) then
          flo0(ix,iy,1,1)   = floe(ix,iy,1,1)+floi(ix,iy,1,1)
          con0(ix,iy,1,1,0) = cone(ix,iy,1,1,0)+coni(ix,iy,1,1,0)     !srv 02.01.07
cwdk add linearization for 9-point terms
          flo0(ix,iy,1,0)   = floe(ix,iy,1,0)   + floi(ix,iy,1,0)
          con0(ix,iy,1,0,0) = cone(ix,iy,1,0,0) + coni(ix,iy,1,0,0)
        else if (.not.solveee(ireg) .and. solveei(ireg)) then
          flo0(ix,iy,1,1)   = floi(ix,iy,1,1)
          con0(ix,iy,1,1,0) = coni(ix,iy,1,1,0)                       !srv 02.01.07
cwdk add linearization for 9-point terms
          flo0(ix,iy,1,0)   = floi(ix,iy,1,0)
          con0(ix,iy,1,0,0) = coni(ix,iy,1,0,0)
        else if (solveee(ireg) .and. .not.solveei(ireg)) then
          flo0(ix,iy,1,1)   = floe(ix,iy,1,1)
          con0(ix,iy,1,1,0) = cone(ix,iy,1,1,0)                       !srv 02.01.07
cwdk add linearization for 9-point terms
          flo0(ix,iy,1,0)   = floe(ix,iy,1,0)
          con0(ix,iy,1,0,0) = cone(ix,iy,1,0,0)
        endif
       enddo
      enddo
      aa9 = 0.0e0_R8
*   ..compute the nine-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveee(ireg) .and. solveei(ireg)) then
          aa9(ix,iy,0) =
     &    -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)-
     &     shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)
     &    +abs(resht(ix,iy))/(rxg*min(te(ix,iy),ti(ix,iy)))           !srv 02.01.07
        else if (solveee(ireg) .and. .not.solveei(ireg)) then
          aa9(ix,iy,0) =
     &    -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)
     &    +abs(resht(ix,iy))/(rxg*te(ix,iy))                          !srv 02.01.07
        else if (.not.solveee(ireg) .and. solveei(ireg)) then
          aa9(ix,iy,0) =
     &    -shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)
     &    +abs(resht(ix,iy))/(rxg*ti(ix,iy))                          !srv 02.01.07
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
         aa9(ix,iy,-4) = 0.0e0_R8
         aa9(ix,iy, 2) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         t1 = -0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveee(ireg) .or. solveei(ireg)) then
          aa9(ix,iy,0)  = aa9(ix,iy,0)  - t0
          aa9(ix,iy,-1) = aa9(ix,iy,-1) + t1
         endif
         if (solveee(ileft) .or. solveei(ileft)) then
          aa9(leftix(ix,iy),leftiy(ix,iy),1) =
     &      aa9(leftix(ix,iy),leftiy(ix,iy),1) + t0
          aa9(leftix(ix,iy),leftiy(ix,iy),0) =
     &      aa9(leftix(ix,iy),leftiy(ix,iy),0) - t1
         endif
         if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
          if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
           t0  = (0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           t1  = (-0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           if (solveee(ireg) .or. solveei(ireg)) then
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
           if (solveee(ileft) .or. solveei(ileft)) then
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
         endif
        endif                                                         !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        t0 =  0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                !srv 02.01.07
        t1 = -0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                !srv 02.01.07
        ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
        if (solveee(ireg) .or. solveei(ireg)) then
         aa9(ix,iy,0)  = aa9(ix,iy,0)  - t0
         aa9(ix,iy,-3) = aa9(ix,iy,-3) + t1
        endif
        if (solveee(ibot) .or. solveei(ibot)) then
         aa9(bottomix(ix,iy),bottomiy(ix,iy),3) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),3) + t0
         aa9(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),0) - t1
        endif
        if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           t0 = (0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
           t1 = (-0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          if (solveee(ireg) .or. solveei(ireg)) then
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
          if (solveee(ibot) .or. solveei(ibot)) then
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
       call b2uxm9 (nx, ny, nregionv, itcnt, solveee, cone(-1,-1,0,0,1),
     &  aa9, 'b2ush9_totale')                                         !srv 02.01.07
      endif
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nregionv, itcnt, solveei, coni(-1,-1,0,0,1),
     &  aa9, 'b2ush9_totali')                                         !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
         aa9(ix,iy,0) = 1.0_R8
         resht(ix,iy) = 0.0_R8
        end if
       enddo
      enddo
*   ..temporary fix to stabilize 9-point stencil in corners:
*   ..set guard cell value equal to neighboring cell
      if (use_9pt_stencil.ne.0) then
      do iy=-1,ny
       do ix=-1,nx
        if (leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:) = 0.0_R8
         aa9(ix,iy,0) = 1.00_R8
!        aa9(ix,iy,3) =-1.00_R8
!        rescb(ix,iy) = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         resht(ix,iy) = 0.5_R8*
     &    (te(rightix(ix,iy),rightiy(ix,iy)) - te(ix,iy) +
     &     ti(rightix(ix,iy),rightiy(ix,iy)) - ti(ix,iy))
        endif
        if (leftix(ix,iy).eq.-2.and.topiy(ix,iy).eq.ny+2) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy, 1) =-0.50_R8
!        aa9(ix,iy,-3) =-0.50_R8
!        rescb(ix,iy)  = 0.5_R8*(pb(rightix(ix,iy),rightiy(ix,iy)) +
!                                pb(bottomix(ix,iy),bottomiy(ix,iy))) -
!                        pb(ix,iy)
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(-1,ny)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         resht(ix,iy) = 0.5_R8*
     &    (te(rightix(ix,iy),rightiy(ix,iy)) - te(ix,iy) +
     &     ti(rightix(ix,iy),rightiy(ix,iy)) - ti(ix,iy))
        endif
        if (rightix(ix,iy).eq.nx+1.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,3)  =-1.00_R8
!        rescb(ix,iy)  = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         resht(ix,iy)  = 0.5_R8*
     &    (te(leftix(ix,iy),leftiy(ix,iy)) - te(ix,iy) +
     &     ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
        endif
        if (rightix(ix,iy).eq.nx+1.and.topiy(ix,iy).eq.ny+1) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(ix,iy)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         resht(ix,iy)  = 0.5_R8*
     &    (te(leftix(ix,iy),leftiy(ix,iy)) - te(ix,iy) +
     &     ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy))
        endif
       enddo
      enddo
      endif ! 9-point stencil

*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, resht, cortt, 'b2ush9_total')  !srv 18.05.02

      if (iout.ne.0) then
       call my_out(70,nx,ny,aa9(-1,-1,-6),'b2ush9_total_aa9m6')
       call my_out(70,nx,ny,aa9(-1,-1,-5),'b2ush9_total_aa9m5')
       call my_out(70,nx,ny,aa9(-1,-1,-4),'b2ush9_total_aa9m4')
       call my_out(70,nx,ny,aa9(-1,-1,-3),'b2ush9_total_aa9m3')
       call my_out(70,nx,ny,aa9(-1,-1,-2),'b2ush9_total_aa9m2')
       call my_out(70,nx,ny,aa9(-1,-1,-1),'b2ush9_total_aa9m1')
       call my_out(70,nx,ny,aa9(-1,-1, 0),'b2ush9_total_aa9m0')
       call my_out(70,nx,ny,aa9(-1,-1, 1),'b2ush9_total_aa9p1')
       call my_out(70,nx,ny,aa9(-1,-1, 2),'b2ush9_total_aa9p2')
       call my_out(70,nx,ny,aa9(-1,-1, 3),'b2ush9_total_aa9p3')
       call my_out(70,nx,ny,aa9(-1,-1, 4),'b2ush9_total_aa9p4')
       call my_out(70,nx,ny,aa9(-1,-1, 5),'b2ush9_total_aa9p5')
       call my_out(70,nx,ny,aa9(-1,-1, 6),'b2ush9_total_aa9p6')
       call my_out(70,nx,ny,resht,'b2ush9_resht')
       call my_out(70,nx,ny,cortt,'b2ush9_cortt')
      endif

      endif                                                           !srv 22.05.18

* ..electron temperature correction equation
      if (any(solveee(0:nnreg(0)))) then                              !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveee(ireg).neqv.solveet(ireg)
      end do
      if (new_matrix.or.use_9pt_stencil.ne.0) call restart_ma28_for_9
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0,0:1) = 0.0e0_R8
         con0(ix,iy,0,0:1,:) = 0.0e0_R8                               !srv 02.01.07
        else
         t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*                               !srv 02.01.07 {
     &    (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &     ne(leftix(ix,iy),leftiy(ix,iy))*
     &     te(leftix(ix,iy),leftiy(ix,iy))+
     &     pccm1(ix,iy)*ne(ix,iy)*te(ix,iy))/2.0_R8                   !srv 02.01.07 }
         flo0(ix,iy,0,0) = floe(ix,iy,0,0)
     &     -t0*(ne(ix,iy)-ne(leftix(ix,iy),leftiy(ix,iy)))            !srv 02.01.07
         con0(ix,iy,0,0,0) = cone(ix,iy,0,0,0)                        !srv 02.01.07
     &     +(t0/2.0_R8)*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))   !srv 02.01.07
cwdk add linearization for 9-point terms
          flo0(ix,iy,0,1)   = floe(ix,iy,0,1)
          con0(ix,iy,0,1,0) = cone(ix,iy,0,1,0)
        endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1,0:1) = 0.0_R8
       con0(ix,-1,1,0:1,:) = 0.0_R8                                   !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        flo0(ix,iy,1,1)   = floe(ix,iy,1,1)
        con0(ix,iy,1,1,0) = cone(ix,iy,1,1,0)                         !srv 02.01.07
cwdk add linearization for 9-point terms
        flo0(ix,iy,1,0)   = floe(ix,iy,1,0)
        con0(ix,iy,1,0,0) = cone(ix,iy,1,0,0)
       enddo
      enddo
*   ..compute the nine-point matrix
      aa9 = 0.0e0_R8                                                  !srv 28.05.08
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveee(ireg))
     >   aa9(ix,iy,0) = -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)+
     &    cor_diag*dim(1.0_R8,1.5_R8*(te(ix,iy)-ti(ix,iy))/te(ix,iy))*  !srv 02.01.07 07.07.21
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))
     &    +abs(reshe(ix,iy))/(rxg*te(ix,iy))                          !srv 02.01.07
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
         aa9(ix,iy,-4) = 0.0e0_R8
         aa9(ix,iy, 2) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         t1 = -0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveee(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),1) = t0
         if (solveee(ireg))
     >    aa9(ix,iy,0) = aa9(ix,iy,0)-t0
         if (solveee(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),0) =
     &       aa9(leftix(ix,iy),leftiy(ix,iy),0)-t1
         if (solveee(ireg))
     >    aa9(ix,iy,-1) = t1
         if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
          if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
           t0  = (0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           t1  = (-0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           if (solveee(ireg)) then
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
           if (solveee(ileft)) then
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
         endif
        endif                                                         !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        t0 =  0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                !srv 02.01.07
        t1 = -0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                !srv 02.01.07
        ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
        if (solveee(ireg)) then
         aa9(ix,iy,0)  = aa9(ix,iy,0)  - t0
         aa9(ix,iy,-3) = aa9(ix,iy,-3) + t1
        endif
        if (solveee(ibot)) then
         aa9(bottomix(ix,iy),bottomiy(ix,iy),3) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),3) + t0
         aa9(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),0) - t1
        endif
        if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           t0 = (0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
           t1 = (-0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          if (solveee(ireg)) then
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
          if (solveee(ibot)) then
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
       call b2uxm9 (nx, ny, nregionv, itcnt, solveee, cone(-1,-1,0,0,1),
     &  aa9, 'b2ush9_Te')                                             !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
         aa9(ix,iy,0) = 1.0_R8
         reshe(ix,iy) = 0.0_R8
        end if
       enddo
      enddo
*   ..temporary fix to stabilize 9-point stencil in corners:
*   ..set guard cell value equal to neighboring cell
      if (use_9pt_stencil.ne.0) then
      do iy=-1,ny
       do ix=-1,nx
        if (leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:) = 0.0_R8
         aa9(ix,iy,0) = 1.00_R8
!        aa9(ix,iy,3) =-1.00_R8
!        rescb(ix,iy) = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         reshe(ix,iy) = te(rightix(ix,iy),rightiy(ix,iy)) - te(ix,iy)
        endif
        if (leftix(ix,iy).eq.-2.and.topiy(ix,iy).eq.ny+1) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy, 1) =-0.50_R8
!        aa9(ix,iy,-3) =-0.50_R8
!        rescb(ix,iy)  = 0.5_R8*(pb(rightix(ix,iy),rightiy(ix,iy)) +
!    .                           pb(bottomix(ix,iy),bottomiy(ix,iy))) -
!    .                   pb(ix,iy)
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(ix,iy)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         reshe(ix,iy) = te(rightix(ix,iy),rightiy(ix,iy)) - te(ix,iy)
        endif
        if (rightix(ix,iy).eq.ny+1.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,3)  =-1.00_R8
!        rescb(ix,iy)  = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         reshe(ix,iy)  = te(leftix(ix,iy),leftiy(ix,iy)) - te(ix,iy)
        endif
        if (rightix(ix,iy).eq.ny+1.and.topiy(ix,iy).eq.ny+1) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(ix,iy)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         reshe(ix,iy)  = te(leftix(ix,iy),leftiy(ix,iy)) - te(ix,iy)
        endif
       enddo
      enddo
      endif ! 9-point stencil

*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, reshe, corte, 'b2ush9_Te')     !srv 18.05.02

      if (iout.ne.0) then
       call my_out(70,nx,ny,aa9(-1,-1,-6),'b2ush9_he_aa9m6')
       call my_out(70,nx,ny,aa9(-1,-1,-5),'b2ush9_he_aa9m5')
       call my_out(70,nx,ny,aa9(-1,-1,-4),'b2ush9_he_aa9m4')
       call my_out(70,nx,ny,aa9(-1,-1,-3),'b2ush9_he_aa9m3')
       call my_out(70,nx,ny,aa9(-1,-1,-2),'b2ush9_he_aa9m2')
       call my_out(70,nx,ny,aa9(-1,-1,-1),'b2ush9_he_aa9m1')
       call my_out(70,nx,ny,aa9(-1,-1, 0),'b2ush9_he_aa9m0')
       call my_out(70,nx,ny,aa9(-1,-1, 1),'b2ush9_he_aa9p1')
       call my_out(70,nx,ny,aa9(-1,-1, 2),'b2ush9_he_aa9p2')
       call my_out(70,nx,ny,aa9(-1,-1, 3),'b2ush9_he_aa9p3')
       call my_out(70,nx,ny,aa9(-1,-1, 4),'b2ush9_he_aa9p4')
       call my_out(70,nx,ny,aa9(-1,-1, 5),'b2ush9_he_aa9p5')
       call my_out(70,nx,ny,aa9(-1,-1, 6),'b2ush9_he_aa9p6')
       call my_out(70,nx,ny,reshe,'b2ush9_reshe')
       call my_out(70,nx,ny,corte,'b2ush9_corte')
      endif

      endif                                                           !srv 22.05.18

* ..atom temperature correction equation
      if (any(solveei(0:nnreg(0)))) then                              !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveei(ireg).neqv.solveee(ireg)
      end do
      if (new_matrix.or.use_9pt_stencil.ne.0) call restart_ma28_for_9
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0,0:1) = 0.0e0_R8
         con0(ix,iy,0,0:1,:) = 0.0e0_R8                               !srv 02.01.07
        else
         t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*                               !srv 02.01.07 {
     &    (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &     ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &     ti(leftix(ix,iy),leftiy(ix,iy))+
     &     pccm1(ix,iy)*ni(ix,iy,0)*ti(ix,iy))/2.0_R8                 !srv 02.01.07 }
         flo0(ix,iy,0,0) = floi(ix,iy,0,0)
     &     -t0*(ni(ix,iy,0)-ni(leftix(ix,iy),leftiy(ix,iy),0))        !srv 02.01.07
         con0(ix,iy,0,0,0) = coni(ix,iy,0,0,0)                        !srv 02.01.07
     &      +(t0/2.0_R8)*(ni(leftix(ix,iy),leftiy(ix,iy),0)+
     &                    ni(ix,iy,0))                                !srv 02.01.07
cwdk add linearization for 9-point terms
         flo0(ix,iy,0,1)   = floi(ix,iy,0,1)
         con0(ix,iy,0,1,0) = coni(ix,iy,0,1,0)
        endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1,0:1) = 0.0_R8
       con0(ix,-1,1,0:1,:) = 0.0_R8                                   !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        flo0(ix,iy,1,1)   = floi(ix,iy,1,1)
        con0(ix,iy,1,1,0) = coni(ix,iy,1,1,0)                             !srv 02.01.07
cwdk add linearization for 9-point terms
        flo0(ix,iy,1,0)   = floi(ix,iy,1,0)
        con0(ix,iy,1,0,0) = coni(ix,iy,1,0,0)
       enddo
      enddo
*   ..compute the nine-point matrix
      aa9 = 0.0e0_R8                                                  !srv 28.05.08
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveei(ireg))
     >   aa9(ix,iy,0) = (-shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)+
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))
     &   +abs(reshi(ix,iy))/(rxg*ti(ix,iy)))*scalr(ti(ix,iy))         !srv 02.01.07
!     &                                      )*scalr(ti(ix,iy))       !srv 02.01.07
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
         aa9(ix,iy,-4) = 0.0e0_R8
         aa9(ix,iy, 2) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         t1 = -0.5_R8*flo0(ix,iy,0,0)-con0(ix,iy,0,0,0)               !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveei(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),1) = t0*scalr(ti(ix,iy))
         if (solveei(ireg))
     >    aa9(ix,iy,0) = aa9(ix,iy,0)-t0*scalr(ti(ix,iy))
         if (solveei(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),0) =
     &      aa9(leftix(ix,iy),leftiy(ix,iy),0)
     &     -t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
         if (solveei(ireg))
     >    aa9(ix,iy,-1) = t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
         if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
cwdk todo: scalr()-functionality; is 1.0 for now
          if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
           t0  = (0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           t1  = (-0.5_R8*flo0(ix,iy,0,1)-con0(ix,iy,0,1,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           if (solveei(ireg)) then
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
           if (solveei(ileft)) then
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
         endif
        endif                                                         !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        t0 =  0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                     !srv 02.01.07
        t1 = -0.5_R8*flo0(ix,iy,1,1)-con0(ix,iy,1,1,0)                     !srv 02.01.07
        ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
        if (solveei(ireg)) then
         aa9(ix,iy,0)  = aa9(ix,iy,0)  - t0*scalr(ti(ix,iy))
         aa9(ix,iy,-3) = aa9(ix,iy,-3) + t1*
     &     scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
        endif
        if (solveei(ibot)) then
         aa9(bottomix(ix,iy),bottomiy(ix,iy),3) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),3) + t0*scalr(ti(ix,iy))
         aa9(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),0) - t1*
     &      scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
        endif
        if (use_9pt_stencil.ne.0) then
cwdk add 9-point contributions
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
          t0 = (0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
          t1 = (-0.5_R8*flo0(ix,iy,1,0)-con0(ix,iy,1,0,0))/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          if (solveei(ireg)) then
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
          if (solveei(ibot)) then
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
       call b2uxm9 (nx, ny, nregionv, itcnt, solveei, coni(-1,-1,0,0,1),
     &  aa9, 'b2ush9_Ti')                                             !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
         aa9(ix,iy,0) = 1.0_R8
         reshi(ix,iy) = 0.0_R8
        end if
       enddo
      enddo
*   ..temporary fix to stabilize 9-point stencil in corners:
*   ..set guard cell value equal to neighboring cell
      if (use_9pt_stencil.ne.0) then
      do iy=-1,ny
       do ix=-1,nx
        if (leftix(ix,iy).eq.-2.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:) = 0.0_R8
         aa9(ix,iy,0) = 1.00_R8
!        aa9(ix,iy,3) =-1.00_R8
!        rescb(ix,iy) = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         reshi(ix,iy) = ti(rightix(ix,iy),rightiy(ix,iy)) - ti(ix,iy)
        endif
        if (leftix(ix,iy).eq.-2.and.topiy(ix,iy).eq.ny+1) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy, 1) =-0.50_R8
!        aa9(ix,iy,-3) =-0.50_R8
!        rescb(ix,iy)  = 0.5_R8*(pb(rightix(ix,iy),rightiy(ix,iy)) +
!    .                           pb(bottomix(ix,iy),bottomiy(ix,iy))) -
!    .                   pb(ix,iy)
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(ix,iy)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,1) =-1.00_R8
         reshi(ix,iy) = ti(rightix(ix,iy),rightiy(ix,iy)) - ti(ix,iy)
        endif
        if (rightix(ix,iy).eq.nx+1.and.bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,3)  =-1.00_R8
!        rescb(ix,iy)  = pb(topix(ix,iy),topiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         reshi(ix,iy)  = ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy)
        endif
        if (rightix(ix,iy).eq.nx+1.and.topiy(ix,iy).eq.ny+1) then
         aa9(ix,iy,:)  = 0.0_R8
         aa9(ix,iy,0)  = 1.00_R8
!        aa9(ix,iy,-3) =-1.00_R8
!        rescb(ix,iy)  = pb(bottomix(ix,iy),bottomiy(ix,iy)) - pb(ix,iy)
         aa9(ix,iy,-1) =-1.00_R8
         reshi(ix,iy)  = ti(leftix(ix,iy),leftiy(ix,iy)) - ti(ix,iy)
        endif
       enddo
      enddo
      endif ! 9-point stencil

*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, reshi, corti, 'b2ush9_Ti')
      last_solve_9(0:nnreg(0)) = solveei(0:nnreg(0))
      do iy = -1, ny
       do ix = -1, nx
        corti(ix,iy) = corti(ix,iy)*scalr(ti(ix,iy))
       enddo
      enddo

      if (iout.ne.0) then
       call my_out(70,nx,ny,aa9(-1,-1,-6),'b2ush9_hi_aa9m6')
       call my_out(70,nx,ny,aa9(-1,-1,-5),'b2ush9_hi_aa9m5')
       call my_out(70,nx,ny,aa9(-1,-1,-4),'b2ush9_hi_aa9m4')
       call my_out(70,nx,ny,aa9(-1,-1,-3),'b2ush9_hi_aa9m3')
       call my_out(70,nx,ny,aa9(-1,-1,-2),'b2ush9_hi_aa9m2')
       call my_out(70,nx,ny,aa9(-1,-1,-1),'b2ush9_hi_aa9m1')
       call my_out(70,nx,ny,aa9(-1,-1, 0),'b2ush9_hi_aa9m0')
       call my_out(70,nx,ny,aa9(-1,-1, 1),'b2ush9_hi_aa9p1')
       call my_out(70,nx,ny,aa9(-1,-1, 2),'b2ush9_hi_aa9p2')
       call my_out(70,nx,ny,aa9(-1,-1, 3),'b2ush9_hi_aa9p3')
       call my_out(70,nx,ny,aa9(-1,-1, 4),'b2ush9_hi_aa9p4')
       call my_out(70,nx,ny,aa9(-1,-1, 5),'b2ush9_hi_aa9p5')
       call my_out(70,nx,ny,aa9(-1,-1, 6),'b2ush9_hi_aa9p6')
       call my_out(70,nx,ny,reshi,'b2ush9_reshi')
       call my_out(70,nx,ny,corti,'b2ush9_corti')
      endif

      endif                                                           !srv 22.05.18

*   ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ush9

      contains

      function scalr(t0)
      implicit none
      real(kind=R8) :: scalr, t0
      scalr = 1.0e0_R8  ! t0                                          !srv 02.01.07
      return
      end function scalr

      end subroutine b2ush9
