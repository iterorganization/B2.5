*-----------------------------------------------------------------------
*.specification

      subroutine b2usp7_ (nx, ny, itcnt, ne, te, conc,                !sv 27.01.00
     &    avisper, bvisper, cvisper,                                  !sv 27.01.00
     &    sch, respo, corpo, aa7)                                     !sv 27.01.00
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, itcnt
      real (kind=R8) ::
     *   ne(-1:nx,-1:ny), te(-1:nx,-1:ny), conc(-1:nx,-1:ny,0:1),
     3   avisper(-1:nx,-1:ny,0:1), bvisper(-1:nx,-1:ny,0:1),        !sv 27.01.00
     4   cvisper(-1:nx,-1:ny,0:1),                                  !sv 27.01.00
     &   sch(-1:nx,-1:ny,0:3), respo(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corpo(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   aa7(-1:nx,-1:ny,-3:3)                               !sv 11.01.00
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USP7_ computes a correction towards the solution of the
*     potential equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ixb, iyb, ixt, iyt, k                !sv 26.01.00
      integer, save :: iout = 0                                       !sv 08.07.08
      character charns*3, chark*1                                     !sv 08.07.08
      real (kind=R8) :: vis_per, aa7p2, aa7m2                         !sv 27.01.00
      save ncall, vis_per                                             !sv 27.01.00
*   ..procedures
      intrinsic min, abs
      external subini, subend, xertst, sfill
      external b2xvsg, b2xvff, b2ux7p                                 !sv 11.01.00
*   ..initialization
      data ncall/0/, vis_per/0.0e0_R8/                                !sv 27.01.00

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usp7_')                                          !sv 06.07.00
*   ..set internal parameters on first call                           !sv
      if (ncall.eq.0) then                                            !sv
       call ipgetr ('b2tfhe_vis_per', vis_per)                        !sv 27.01.00 08.07.08
       call ipgeti ('b2usp7_iout', iout)                              !sv 08.07.08
      endif                                                           !sv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, conc, sch
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2*2, conc, 1, 'conc', '.ge.')
       call b2xvsg (n2, sch(-1,-1,1), 1, 'sch1', '.le.')
       call b2xvsg (n2, sch(-1,-1,3), 1, 'sch3', '.le.')
*    ..test edge values of conc
       call b2xvff (nx, ny, conc, 'conc')
      endif

* ..compute the correction
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        aa7(ix,iy,0) = -sch(ix,iy,1)-sch(ix,iy,3)*ne(ix,iy)
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                                 !sv 26.11.02
        if(leftix(ix,iy).eq.-2) then
         aa7(ix,iy,-1) = 0.0e0_R8
        else
         aa7(leftix(ix,iy),leftiy(ix,iy),1) = -conc(ix,iy,0)
         aa7(ix,iy,0) = aa7(ix,iy,0)+conc(ix,iy,0)
         aa7(leftix(ix,iy),leftiy(ix,iy),0) = 
     &       aa7(leftix(ix,iy),leftiy(ix,iy),0)+conc(ix,iy,0)
         aa7(ix,iy,-1) = -conc(ix,iy,0)
        endif                                                         !sv 26.11.02
        if(rightix(ix,iy).eq.nx+1) aa7(ix,iy,1) = 0.0e0_R8            !sv 26.11.02
       enddo
      enddo
      call sfill (nx+2, 0.0e0_R8, aa7(-1,-1,-2), 1)
      do iy = 0, ny
       do ix = -1, nx
        aa7(bottomix(ix,iy),bottomiy(ix,iy),2) = -conc(ix,iy,1)
        aa7(ix,iy,0) = aa7(ix,iy,0)+conc(ix,iy,1)
        aa7(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &       aa7(bottomix(ix,iy),bottomiy(ix,iy),0)+conc(ix,iy,1)
        aa7(ix,iy,-2) = -conc(ix,iy,1)
       enddo
      enddo
      call sfill (nx+2, 0.0e0_R8, aa7(-1,ny,2), 1)
	aa7(:,:,-3) = 0.0e0_R8
	aa7(:,:, 3) = 0.0e0_R8
*   ..compute the seven-point matrix                                !sv 26.01.00
      call sfill ((nx+2)*2, 0.0e0_R8, aa7(-1,-1,-3), 1)
      aa7(:,-1,3) = 0.0e0_R8

	iy = 0
	ix = -1
      if(leftix(ix,iy).eq.-2) aa7(ix,iy,3) = 0.0e0_R8                 !sv 26.11.02
	do ix = 0, nx-1
	  ixt = topix(ix,iy)
	  iyt = topiy(ix,iy)
        aa7(ix,iy,3) = vis_per*(                                       !sv 27.01.00
     &    avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)
     &   *cvisper(topix(ixt,iyt),topiy(ixt,iyt),1) )                   !sv 27.01.00
        aa7p2 =     vis_per*(                                          !sv 27.01.00
     &    avisper(ixt,iyt,1)*bvisper(ixt,iyt,1)*cvisper(ixt,iyt,1)     !sv 27.01.00
     &   +avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*cvisper(ixt,iyt,1) )   !sv 27.01.00
        aa7(ix,iy,2) = aa7(ix,iy,2) - aa7(ix,iy,3) - aa7p2             !sv 27.01.00
!        aa7m2 =   vis_per*(                                            !sv 27.01.00 06.05.02
!     &    avisper(ixt,iyt,1)*bvisper(ix ,iy ,1)*cvisper(ix ,iy ,1) )   !sv 27.01.00 06.05.02
        aa7m2 = 0.0e0_R8                                               !sv 06.05.02
        aa7(ix,iy,-2) = aa7(ix,iy,-2)  - aa7m2                         !sv 27.01.00
        aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2                     !sv 27.01.00
	enddo
	ix=nx
      if(rightix(ix,iy).eq.nx+1) aa7(nx,iy,3) = 0.0e0_R8              !sv 26.11.02

      do iy = 1, ny-2
       do ix = -1, nx                                                 !sv 26.11.02
	  if(leftix(ix,iy).eq.-2.or.rightix(ix,iy).eq.nx+1) then        !sv 26.11.02
         aa7(ix,iy, 3) = 0.0e0_R8                                      !sv 27.01.00
         aa7(ix,iy,-3) = 0.0e0_R8                                      !sv 27.01.00
	  else
	   ixt = topix(ix,iy)
	   iyt = topiy(ix,iy)
	   ixb = bottomix(ix,iy)
	   iyb = bottomiy(ix,iy)
         if (iy.lt.ny-2) then                                           !sv 06.05.02
          aa7(ix,iy,3) = vis_per*(                                      !sv 27.01.00
     &     avisper(ixt,iyt,1)                                           !sv 27.01.00
     &    *bvisper(ixt,iyt,1)*cvisper(topix(ixt,iyt),topiy(ixt,iyt),1)  !sv 27.01.00
     &                                                     )            !sv 27.01.00
         else                                                           !sv 06.05.02
          aa7(ix,iy,3) = 0.0e0_R8                                       !sv 06.05.02
	   endif                                                          !sv 06.05.02
	   if (iy.ge.2) then                                             !sv 06.05.02
	    aa7(ix,iy,-3) = vis_per*(                                    !sv 27.01.00
     &     avisper(ix,iy,1)*bvisper(ixb,iyb,1)                         !sv 27.01.00
     &    *cvisper(ixb,iyb,1)                                          !sv 27.01.00
     &                                           )                     !sv 27.01.00
         else                                                          !sv 06.05.02
	    aa7(ix,iy,-3) = 0.0e0_R8  !sv because of boundary condition  !sv 06.05.02
	   endif                                                         !sv 06.05.02
	   aa7p2 =     vis_per*(                                          !sv 27.01.00
     &     avisper(ixt,iyt,1)                                           !sv 27.01.00
     &    *bvisper(ixt,iyt,1)*cvisper(ixt,iyt,1)                        !sv 27.01.00
     &    +avisper(ixt,iyt,1)                                           !sv 27.01.00
     &    *bvisper(ix ,iy ,1)*cvisper(ixt,iyt,1)                        !sv 27.01.00
     &    +avisper(ix ,iy ,1)                                           !sv 27.01.00
     &    *bvisper(ix ,iy ,1)*cvisper(ixt,iyt,1)                        !sv 27.01.00
     &                                                     )            !sv 27.01.00
         aa7(ix,iy,2) = aa7(ix,iy,2) - aa7(ix,iy,3) - aa7p2             !sv 27.01.00

         aa7m2 =   vis_per*(                                            !sv 27.01.00
     &     avisper(ixt,iyt,1)                                           !sv 27.01.00
     &    *bvisper(ix ,iy ,1)*cvisper(ix ,iy ,1)                        !sv 27.01.00
     &    +avisper(ix ,iy ,1)                                           !sv 27.01.00
     &    *bvisper(ix ,iy ,1)*cvisper(ix ,iy ,1)                        !sv 27.01.00
     &    +avisper(ix ,iy ,1)                                           !sv 27.01.00
     &    *bvisper(ixb,iyb,1)*cvisper(ix ,iy ,1)                        !sv 27.01.00
     &                                                     )            !sv 27.01.00
         aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7(ix,iy,-3) - aa7m2          !sv 27.01.00

         aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2                     !sv 27.01.00
        endif                                                         !sv 26.11.02
       enddo
      enddo

      aa7(:,ny,-3) = 0.0e0_R8

	iy = ny-1
	ix=-1                                                           !sv 26.11.02
	if(leftix(ix,iy).eq.-2) aa7(ix,iy,-3) = 0.0e0_R8                !sv 26.11.02
       
	do ix = 0, nx-1
	  ixt = topix(ix,iy)
	  iyt = topiy(ix,iy)
	  ixb = bottomix(ix,iy)
	  iyb = bottomiy(ix,iy)
        aa7(ix,iy,-3) = vis_per*(                                     !sv 27.01.00
     &    avisper(ix,iy,1)*bvisper(ixb,iyb,1)*cvisper(ixb,iyb,1)      !sv 27.01.00
     &                                          )                     !sv 27.01.00
        aa7m2 =   vis_per*(                                            !sv 27.01.00
     &    avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*cvisper(ix ,iy ,1)     !sv 27.01.00
     &   +avisper(ix ,iy ,1)*bvisper(ixb,iyb,1)*cvisper(ix ,iy ,1)     !sv 27.01.00
     &                                                    )            !sv 27.01.00
        aa7(ix,iy,-2) = aa7(ix,iy,-2)  - aa7(ix,iy,-3) - aa7m2   !!!   !sv 27.01.00  06.05.02 
!        aa7p2 =     vis_per*(                                          !sv 27.01.00 06.05.02
!     &    avisper(ix ,iy ,1)*bvisper(ix ,iy ,1)*cvisper(ixt,iyt,1)     !sv 27.01.00 06.05.02
!     &                                                    )            !sv 27.01.00 06.05.02
        aa7p2 = 0.0e0_R8                                               !sv 06.05.02
        aa7(ix,iy,2) = aa7(ix,iy,2)  - aa7p2                           !sv 27.01.00

        aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2                     !sv 27.01.00
	enddo

	ix=nx                                                          !sv 26.11.02
	if(rightix(ix,iy).eq.nx+1) aa7(ix,iy,3) = 0.0e0_R8             !sv 26.11.02
       
      call sfill ((nx+2)*2, 0.0e0_R8, aa7(-1,ny-1,3), 1)               !sv 26.01.00
	if (iout.ne.0) then                                           !sv 08.07.08 {
	 do k=0,3
	  write (chark,'(i1)') k  
	  call my_out(70,nx,ny,sch(-1,-1,k),'b2usp7__sch'//chark)
	 enddo
	  call my_out(70,nx,ny,aa7(-1,-1,-3),'b2usp7__aa7m3')
	  call my_out(70,nx,ny,aa7(-1,-1,-2),'b2usp7__aa7m2')
	  call my_out(70,nx,ny,aa7(-1,-1,-1),'b2usp7__aa7m1')
	  call my_out(70,nx,ny,aa7(-1,-1, 0),'b2usp7__aa7m0')
	  call my_out(70,nx,ny,aa7(-1,-1, 1),'b2usp7__aa7p1')
	  call my_out(70,nx,ny,aa7(-1,-1, 2),'b2usp7__aa7p2')
	  call my_out(70,nx,ny,aa7(-1,-1, 3),'b2usp7__aa7p3')
	  call my_out(70,nx,ny,respo,'b2usp7__respo')
	  call my_out(70,nx,ny,corpo,'b2usp7__corpo')
	endif                                                         !sv 08.07.08 }

*   ..solve the correction equation
      call b2ux7p (nx, ny, aa7, itcnt, respo, corpo, 'b2usp7_')       !sv 11.01.00

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usp7_

      end
