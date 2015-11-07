*-----------------------------------------------------------------------
*.specification

      subroutine b2ush9 (nx, ny, itcnt, rxg, pcm1, pbs, ne, ni, te, ti,
     &  floe, cone, floi, coni, pccm1, she, shi, ceqp,
     &  resht, reshe, reshi, cortt, corte, corti, flo0, con0, aa9)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      real (kind=R8) ::
     *   rxg, pcm1, pbs(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  ni(-1:nx,-1:ny,0:1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  floe(-1:nx,-1:ny,0:1), floi(-1:nx,-1:ny,0:1),
     &  cone(-1:nx,-1:ny,0:1,0:2), coni(-1:nx,-1:ny,0:1,0:2),         !srv 02.01.07
     &  pccm1(-1:nx,-1:ny)
      real (kind=R8) ::
     *   she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3),
     &  ceqp(-1:nx,-1:ny), resht(-1:nx,-1:ny), reshe(-1:nx,-1:ny),
     &  reshi(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   cortt(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(-1:nx,-1:ny,0:1), con0(-1:nx,-1:ny,0:1,0:2),             !srv 02.01.07
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
	integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0            !srv 06.07.06
	integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0            !srv 06.07.06
	integer, save :: iadd_9 = 1                                      !srv 02.01.07
      save ncall
      real (kind=R8) ::
     *   t0, t1
*   ..procedures
      intrinsic abs, sqrt, dim
      real (kind=R8) ::
     *   scalr
      external subini, subend, xertst, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux9p, b2uxm9                 !srv 02.01.07
      scalr(t0) = 1.0e0_R8  ! t0                                      !srv 02.01.07
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ush9')
      if (ncall.eq.0) then
       call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)               !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)               !srv 06.07.06
       call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)               !srv 06.07.06
       call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)               !srv 06.07.06
       call ipgeti ('b2ush9_add_9', iadd_9)                           !srv 02.01.07
      endif
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
       call b2xvff (nx, ny, pbs, 'pbs')
*    ..test sign of ne, ni, te, ti
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test sign of cone, coni
       if(no_Ptncr_x_he.eq.0 .and. no_Ptncr_y_he.eq.0) then         !srv 06.07.06 {
        call b2xvsg (n2*2, cone, 1, 'cone', '.ge.')               !srv 02.01.07
	 endif                                                        !srv 06.07.06 }
       if(no_Ptncr_x_hi.eq.0 .and. no_Ptncr_y_hi.eq.0) then         !srv 06.07.06 {
        call b2xvsg (n2*2, coni, 1, 'coni', '.ge.')               !srv 02.01.07
	 endif                                                        !srv 06.07.06 }
*    ..test edge values of floe, floi, cone, coni
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0), 'cone0')               !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,2), 'coni2')               !srv 02.01.07 }
*    ..test sign of she, shi, ceqp
       call b2xvsg (n2, she(-1,-1,0), 1, 'she0', '.ge.')
       call b2xvsg (n2, she(-1,-1,1), 1, 'she1', '.le.')
       call b2xvsg (n2, she(-1,-1,2), 1, 'she2', '.ge.')
       call b2xvsg (n2, she(-1,-1,3), 1, 'she3', '.le.')
       call b2xvsg (n2, shi(-1,-1,0), 1, 'shi0', '.ge.')
       call b2xvsg (n2, shi(-1,-1,1), 1, 'shi1', '.le.')
       call b2xvsg (n2, shi(-1,-1,2), 1, 'shi2', '.ge.')
       call b2xvsg (n2, shi(-1,-1,3), 1, 'shi3', '.le.')
       call b2xvsg (n2, ceqp, 1, 'ceqp', '.ge.')
      endif

* ..total heat correction equation
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0e0_R8
         con0(ix,iy,0,:) = 0.0e0_R8                                    !srv 02.01.07
        else
         t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*pccm1(ix,iy)*                  !srv 02.01.07 {
     &    (ne(leftix(ix,iy),leftiy(ix,iy))*
     &     te(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy)*te(ix,iy)+
     &     ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &     ti(leftix(ix,iy),leftiy(ix,iy))+ni(ix,iy,0)*ti(ix,iy))/2   !srv 02.01.07 }
         flo0(ix,iy,0) = floe(ix,iy,0)+floi(ix,iy,0)
     &    -t0*(ne(ix,iy)+ni(ix,iy,0)-ne(leftix(ix,iy),leftiy(ix,iy))  !srv 02.01.07 {
     &             -ni(leftix(ix,iy),leftiy(ix,iy),0))                !srv 02.01.07 }
         con0(ix,iy,0,0) = cone(ix,iy,0,0)+coni(ix,iy,0,0)               !srv 02.01.07
     &    +(t0/2)*(ne(leftix(ix,iy),leftiy(ix,iy))                    !srv 02.01.07 {
     &     +ne(ix,iy)+ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))  !srv 02.01.07 }
	  endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1) = 0
       con0(ix,-1,1,:) = 0                                            !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        flo0(ix,iy,1) = floe(ix,iy,1)+floi(ix,iy,1)
        con0(ix,iy,1,0) = cone(ix,iy,1,0)+coni(ix,iy,1,0)     !srv 02.01.07
       enddo
      enddo
	aa9 = 0.0e0_R8
*   ..compute the nine-point matrix
      do iy = -1, ny
       do ix = -1, nx
        aa9(ix,iy,0) =
     &    -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)-
     &     shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)
     &    +abs(resht(ix,iy))/(rxg*min(te(ix,iy),ti(ix,iy)))          !srv 02.01.07
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         aa9(leftix(ix,iy),leftiy(ix,iy),1) = t0
         aa9(ix,iy,0) = aa9(ix,iy,0)-t0
         aa9(leftix(ix,iy),leftiy(ix,iy),0) = 
     &      aa9(leftix(ix,iy),leftiy(ix,iy),0)-t1
         aa9(ix,iy,-1) = t1
        endif                                                         !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        t0 =  flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        t1 = -flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        aa9(bottomix(ix,iy),bottomiy(ix,iy),3) = t0
        aa9(ix,iy,0) = aa9(ix,iy,0)-t0
        aa9(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &      aa9(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
        aa9(ix,iy,-3) = t1
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,ny,3), 1)
*   ..add matrix from contribution nine point stencil
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, itcnt, cone(-1,-1,0,1), aa9,
     &                                               'b2ush9_totale')   !srv 02.01.07
      endif
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, itcnt, coni(-1,-1,0,1), aa9,
     &                                               'b2ush9_totali')   !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8) aa9(ix,iy,0) = 1.0_R8
       enddo
      enddo
*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, resht, cortt, 'b2ush9_total')   !srv 18.05.02

* ..electron temperature correction equation
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0e0_R8
         con0(ix,iy,0,:) = 0.0e0_R8     !srv 02.01.07
        else
         t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*pccm1(ix,iy)*                !srv 02.01.07 {
     &    (ne(leftix(ix,iy),leftiy(ix,iy))*
     &     te(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy)*te(ix,iy))/2   !srv 02.01.07 }
         flo0(ix,iy,0) = floe(ix,iy,0)
     &     -t0*(ne(ix,iy)-ne(leftix(ix,iy),leftiy(ix,iy)))          !srv 02.01.07
         con0(ix,iy,0,0) = cone(ix,iy,0,0)     !srv 02.01.07
     &     +(t0/2)*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))      !srv 02.01.07
	  endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1) = 0
       con0(ix,-1,1,:) = 0     !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        flo0(ix,iy,1) = floe(ix,iy,1)
        con0(ix,iy,1,0) = cone(ix,iy,1,0)     !srv 02.01.07
       enddo
      enddo
*   ..compute the nine-point matrix
	aa9 = 0.0e0_R8                                                 !srv 28.05.08
      do iy = -1, ny
       do ix = -1, nx
        aa9(ix,iy,0) = -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)+
     &    dim(1.0_R8,1.5_R8*(te(ix,iy)-ti(ix,iy))/te(ix,iy))*        !srv 02.01.07
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))
     &    +abs(reshe(ix,iy))/(rxg*te(ix,iy))                         !srv 02.01.07
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         aa9(leftix(ix,iy),leftiy(ix,iy),1) = t0
         aa9(ix,iy,0) = aa9(ix,iy,0)-t0
         aa9(leftix(ix,iy),leftiy(ix,iy),0) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy),0)-t1
         aa9(ix,iy,-1) = t1
	  endif                                                         !srv 26.11.02 {
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        t0 =  flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        t1 = -flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        aa9(bottomix(ix,iy),bottomiy(ix,iy),3) = t0
        aa9(ix,iy,0) = aa9(ix,iy,0)-t0
        aa9(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &       aa9(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
        aa9(ix,iy,-3) = t1
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,ny,3), 1)
*   ..add matrix from contribution nine point stencil
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, itcnt, cone(-1,-1,0,1), aa9, 'b2ush9_Te')     !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8) aa9(ix,iy,0) = 1.0_R8
       enddo
      enddo
*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, reshe, corte, 'b2ush9_Te')     !srv 18.05.02

* ..atom temperature correction equation
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0e0_R8
         con0(ix,iy,0,:) = 0.0e0_R8     !srv 02.01.07
        else
         t0 = 2.5_R8*pbs(ix,iy,0)*pcm1*pccm1(ix,iy)*                 !srv 02.01.07 {
     &    (ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &     ti(leftix(ix,iy),leftiy(ix,iy))+ni(ix,iy,0)*ti(ix,iy))/2  !srv 02.01.07 }
         flo0(ix,iy,0) = floi(ix,iy,0)
     &     -t0*(ni(ix,iy,0)-ni(leftix(ix,iy),leftiy(ix,iy),0))       !srv 02.01.07
         con0(ix,iy,0,0) = coni(ix,iy,0,0)     !srv 02.01.07
     &      +(t0/2)*(ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))  !srv 02.01.07
	  endif                                                         !srv 26.11.02 {
       enddo
      enddo
*   ..compute flo0y, con0y
      do ix = -1, nx
       flo0(ix,-1,1) = 0
       con0(ix,-1,1,:) = 0     !srv 02.01.07
      enddo
      do iy = 0, ny
       do ix = -1, nx
        flo0(ix,iy,1) = floi(ix,iy,1)
        con0(ix,iy,1,0) = coni(ix,iy,1,0)     !srv 02.01.07
       enddo
      enddo
*   ..compute the nine-point matrix
	aa9 = 0.0e0_R8                                                 !srv 28.05.08
      do iy = -1, ny
       do ix = -1, nx
        aa9(ix,iy,0) = (-shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)+
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))
     &   +abs(reshi(ix,iy))/(rxg*ti(ix,iy)))*scalr(ti(ix,iy))        !srv 02.01.07
!     &                                      )*scalr(ti(ix,iy))        !srv 02.01.07
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !srv 26.11.02 {
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0e0_R8
        else                                                          !srv 26.11.02 }
         t0 =  flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2-con0(ix,iy,0,0)     !srv 02.01.07
         aa9(leftix(ix,iy),leftiy(ix,iy),1) = t0*scalr(ti(ix,iy))
         aa9(ix,iy,0) = aa9(ix,iy,0)-t0*scalr(ti(ix,iy))
         aa9(leftix(ix,iy),leftiy(ix,iy),0) = 
     &      aa9(leftix(ix,iy),leftiy(ix,iy),0)
     &     -t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
         aa9(ix,iy,-1) = t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
	  endif                                                         !srv 26.11.02 {
        if(rightix(ix,iy).eq.nx+1) aa9(ix,iy,1) = 0.0e0_R8            !srv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,-1,-3), 1)
      do iy = 0, ny
       do ix = -1, nx
        t0 =  flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        t1 = -flo0(ix,iy,1)/2-con0(ix,iy,1,0)     !srv 02.01.07
        aa9(bottomix(ix,iy),bottomiy(ix,iy),3) = t0*scalr(ti(ix,iy))
        aa9(ix,iy,0) = aa9(ix,iy,0)-t0*scalr(ti(ix,iy))
        aa9(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &    aa9(bottomix(ix,iy),bottomiy(ix,iy),0)
     &    -t1*scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
        aa9(ix,iy,-3) = t1*scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
       enddo
      enddo
      call sfill (nx+2, 0.0_R8, aa9(-1,ny,3), 1)
*   ..add matrix from contribution nine point stencil
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, itcnt, coni(-1,-1,0,1), aa9, 'b2ush9_Ti')   !srv 02.01.07
      endif
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8) aa9(ix,iy,0) = 1.0_R8
       enddo
      enddo
*   ..solve the correction equation
      call b2ux9p (nx, ny, aa9, itcnt, reshi, corti, 'b2ush9_Ti')
      do iy = -1, ny
       do ix = -1, nx
        corti(ix,iy) = corti(ix,iy)*scalr(ti(ix,iy))
       enddo
      enddo

*   ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ush9

      end
