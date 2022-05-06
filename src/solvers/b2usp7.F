*-----------------------------------------------------------------------
*.specification

      subroutine b2usp7 (nx, ny, ns, ismain, nregionv, solvereg, itcnt,  !srv 27.01.00
     &    ne, te, conc, avisper_a, bvisper_a, cvisper_a,                 !srv 27.01.00
     &    sch, respo, corpo, aa7)                                        !srv 27.01.00
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, nregionv, itcnt
      real (kind=R8) ::
     *   ne(-1:nx,-1:ny), te(-1:nx,-1:ny), conc(-1:nx,-1:ny,0:1),
     &   avisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   bvisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   cvisper_a(-1:nx,-1:ny,0:1,0:ns-1),
     &   sch(-1:nx,-1:ny,0:3), respo(-1:nx,-1:ny)
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corpo(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   aa7(-1:nx,-1:ny,-3:3)                                     !srv 11.01.00
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USP7 computes a correction towards the solution of the
*     potential equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ixb, iyb, ixt, iyt, k         !srv 26.01.00
      integer, save :: iout = 0                                    !srv 08.07.08
      integer, save :: visper_all_ions = 1
      character chark*1                                            !srv 08.07.08
      integer ireg, ileft, ibot
      real (kind=R8) :: aa7p2, aa7m2, aa7p3, aa7m3                 !srv 27.01.00
      save ncall                                                   !srv 27.01.00
*   ..procedures
      external xertst, sfill
      external b2xvsg, b2xvff, b2ux7p                              !srv 11.01.00
*   ..initialisation
      data ncall/0/                                                !srv 27.01.00

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usp7')                                       !srv 06.07.00
*   ..set internal parameters on first call                        !srv
      if (ncall.eq.0) then                                         !srv
       call ipgeti ('b2usp7_iout', iout)                           !srv 08.07.08
       call ipgeti ('b2tvspr_visper_all_ions', visper_all_ions)
      endif                                                        !srv
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
      aa7 = 0.0_R8
*   ..compute contribution from the five-point stencil
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         aa7(ix,iy,0) = -sch(ix,iy,1)-sch(ix,iy,3)*ne(ix,iy)
        else
         aa7(ix,iy,0) = 1.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx                                              !srv 26.11.02
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa7(ix,iy,-1) = 0.0_R8
        else
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solvereg(ileft))
     >    aa7(leftix(ix,iy),leftiy(ix,iy),1) = -conc(ix,iy,0)
         if (solvereg(ireg))
     >    aa7(ix,iy,0) = aa7(ix,iy,0)+conc(ix,iy,0)
         if (solvereg(ileft))
     >    aa7(leftix(ix,iy),leftiy(ix,iy),0) =
     &       aa7(leftix(ix,iy),leftiy(ix,iy),0)+conc(ix,iy,0)
         if (solvereg(ireg))
     >    aa7(ix,iy,-1) = -conc(ix,iy,0)
        endif                                                      !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg))
     >   aa7(ix,iy,1) = 0.0_R8                                     !srv 26.11.02
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
          aa7(ix,iy,-2) = 0.0_R8
        else
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot)) then
           aa7(bottomix(ix,iy),bottomiy(ix,iy),2) =
     &       aa7(bottomix(ix,iy),bottomiy(ix,iy),2)-conc(ix,iy,1)
           aa7(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &       aa7(bottomix(ix,iy),bottomiy(ix,iy),0)+conc(ix,iy,1)
          endif
          if (solvereg(ireg)) then
           aa7(ix,iy, 0) = aa7(ix,iy, 0)+conc(ix,iy,1)
           aa7(ix,iy,-2) = aa7(ix,iy,-2)-conc(ix,iy,1)
         endif
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solvereg(ireg)) then
          aa7(ix,iy,2) = 0.0_R8
        endif
       enddo
      enddo

*   ..compute contribution from the seven-point y-stencil
      
      do is = 0,ns-1
       if( (visper_all_ions.ne.0 .and. .not.is_neutral(is)) .or.
     &     (visper_all_ions.eq.0 .and. is.eq.ismain)) then
        do iy = -1, ny
         do ix = -1, nx                                              !srv 26.11.02
           ireg = region(ix,iy,0)
           ixt = topix(ix,iy)
           iyt = topiy(ix,iy)
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &      iyt.eq.ny+1 .or. iyb.eq.-2 .or. .not.solvereg(ireg)) then   !srv 26.11.02 12.12.17 {
*   .. do nothing
           else
            if(bottomiy(ixb,iyb).eq.-2) then
              aa7p3 =
     &           avisper_a(ixt,iyt,1,is)*bvisper_a(ixt,iyt,1,is)*
     &           cvisper_a(topix(ixt,iyt),topiy(ixt,iyt),1,is)
              aa7(ix,iy, 3) = aa7(ix,iy, 3) + aa7p3
              aa7m2 =avisper_a(ixt,iyt,1,is)*
     &               bvisper_a(ix,iy,1,is)*cvisper_a(ix,iy,1,is)
              aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7m2
              aa7p2 =avisper_a(ixt,iyt,1,is)*
     &          bvisper_a(ixt,iyt,1,is)*cvisper_a(ixt,iyt,1,is)
     &         +avisper_a(ixt,iyt,1,is)*
     &          bvisper_a(ix ,iy ,1,is)*cvisper_a(ixt,iyt,1,is)
              aa7(ix,iy,2) = aa7(ix,iy,2) - aa7p3 - aa7p2
              aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2
            elseif(topiy(ixt,iyt).eq.ny+1) then
              aa7m3 = avisper_a(ix,iy,1,is)*
     &          bvisper_a(ixb,iyb,1,is)*cvisper_a(ixb,iyb,1,is)
              aa7(ix,iy,-3) = aa7(ix,iy,-3) + aa7m3
              aa7m2 = avisper_a(ix,iy,1,is)*
     &          bvisper_a(ix,iy,1,is)*cvisper_a(ix,iy,1,is)
     &         +avisper_a(ix,iy,1,is)*
     &          bvisper_a(ixb,iyb,1,is)*cvisper_a(ix,iy,1,is)
              aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7m3 - aa7m2
              aa7p2 = avisper_a(ix,iy,1,is)*
     &         bvisper_a(ix,iy,1,is)*cvisper_a(ixt,iyt,1,is)
              aa7(ix,iy,2) = aa7(ix,iy,2) - aa7p2
              aa7(ix,iy,0) = aa7(ix,iy,0)+ aa7p2 + aa7m2
            else
              aa7m3 = avisper_a(ix,iy,1,is)*
     &          bvisper_a(ixb,iyb,1,is)*cvisper_a(ixb,iyb,1,is)
              aa7(ix,iy,-3) = aa7(ix,iy,-3) + aa7m3
              aa7p3 = avisper_a(ixt,iyt,1,is)*bvisper_a(ixt,iyt,1,is)*
     &          cvisper_a(topix(ixt,iyt),topiy(ixt,iyt),1,is)
              aa7(ix,iy, 3) = aa7(ix,iy, 3) + aa7p3
              aa7m2 = avisper_a(ixt,iyt,1,is)*
     &          bvisper_a(ix,iy,1,is)*cvisper_a(ix,iy,1,is)
     &         +avisper_a(ix,iy,1,is)*
     &          bvisper_a(ix,iy,1,is)*cvisper_a(ix,iy,1,is)
     &         +avisper_a(ix,iy,1,is)*
     &          bvisper_a(ixb,iyb,1,is)*cvisper_a(ix,iy,1,is)
              aa7(ix,iy,-2) = aa7(ix,iy,-2) - aa7m3  - aa7m2
              aa7p2 = avisper_a(ixt,iyt,1,is)*
     &          bvisper_a(ixt,iyt,1,is)*cvisper_a(ixt,iyt,1,is)
     &         +avisper_a(ixt,iyt,1,is)*
     &          bvisper_a(ix,iy,1,is)*cvisper_a(ixt,iyt,1,is)
     &         +avisper_a(ix,iy,1,is)*
     &          bvisper_a(ix,iy,1,is)*cvisper_a(ixt,iyt,1,is)
              aa7(ix,iy,2) = aa7(ix,iy,2) - aa7p3 - aa7p2
              aa7(ix,iy,0) = aa7(ix,iy,0) + aa7p2 + aa7m2
            endif
           endif                                                     !srv 26.11.02 12.12.17 }
         enddo
        enddo
       endif
      enddo

*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (aa7(ix,iy,0).eq.0.0_R8 .and.
     &      aa7(ix,iy,1).eq.0.0_R8 .and. aa7(ix,iy,-1).eq.0.0_R8 .and.
     &      aa7(ix,iy,2).eq.0.0_R8 .and. aa7(ix,iy,-2).eq.0.0_R8 .and.
     &      aa7(ix,iy,3).eq.0.0_R8 .and. aa7(ix,iy,-3).eq.0.0_R8) then
          aa7(ix,iy,0) = 1.0_R8
          respo(ix,iy) = 0.0_R8
        end if
        if (.not.solvereg(ireg)) respo(ix,iy) = 0.0_R8
       enddo
      enddo

*   ..solve the correction equation
      call b2ux7p (nx, ny, aa7, itcnt, respo, corpo, 'b2usp7')     !srv 11.01.00

      if (iout.ne.0) then                                          !srv 08.07.08 {
        do k=0,3
          write (chark,'(i1)') k
          call my_out(70,nx,ny,sch(-1,-1,k),'b2usp7_sch'//chark)
        enddo
        call my_out(70,nx,ny,aa7(-1,-1,-3),'b2usp7_aa7m3')
        call my_out(70,nx,ny,aa7(-1,-1,-2),'b2usp7_aa7m2')
        call my_out(70,nx,ny,aa7(-1,-1,-1),'b2usp7_aa7m1')
        call my_out(70,nx,ny,aa7(-1,-1, 0),'b2usp7_aa7m0')
        call my_out(70,nx,ny,aa7(-1,-1, 1),'b2usp7_aa7p1')
        call my_out(70,nx,ny,aa7(-1,-1, 2),'b2usp7_aa7p2')
        call my_out(70,nx,ny,aa7(-1,-1, 3),'b2usp7_aa7p3')
        call my_out(70,nx,ny,respo,'b2usp7_respo')
        call my_out(70,nx,ny,corpo,'b2usp7_corpo')
      endif                                                        !srv 08.07.08 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usp7

      end subroutine b2usp7
