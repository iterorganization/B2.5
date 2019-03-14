*-----------------------------------------------------------------------
*.specification

      subroutine b2uspo (nx, ny, nregionv, solvereg, itcnt, vol,
     &  ne, te, conc, sch, respo, corpo, aa9, name)                  !srv 18.03.10
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv, itcnt
      character*(*) name                                             !srv 18.03.10
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), te(-1:nx,-1:ny), 
     &  conc(-1:nx,-1:ny,0:1,0:1),
     &  sch(-1:nx,-1:ny,0:3), respo(-1:nx,-1:ny),
     &  vol(-1:nx,-1:ny)
      logical solvereg(0:nregionv)      
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpo(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa9(-1:nx,-1:ny,-6:6)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USPO computes a correction towards the solution of the
*     potential equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop, k
      integer, save :: iout = 0
      character chark*1
      real (kind=R8) :: t0, t1
      save ncall
*   ..procedures
      intrinsic min
      external subini, subend, xertst, sfill
      external b2xvsg, b2xvff, b2ux5p
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2uspo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2uspo_iout', iout)
      endif  
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
       call b2xvsg (n2, conc(-1,-1,0,0), 1, 'conc00', '.ge.')
       call b2xvsg (n2, conc(-1,-1,1,1), 1, 'conc11', '.ge.')
       call b2xvsg (n2, sch(-1,-1,1), 1, 'sch1', '.le.')
       call b2xvsg (n2, sch(-1,-1,3), 1, 'sch3', '.le.')
*    ..test edge values of conc
       call b2xvff (nx, ny, conc, 'conc')
      endif

* ..compute the correction
      aa9 = 0.0_R8
*   ..compute the nine-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         aa9(ix,iy,0) = -sch(ix,iy,1)-sch(ix,iy,3)*ne(ix,iy)
        else
         aa9(ix,iy,0) = 1.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-4) = 0.0_R8
         aa9(ix,iy,-1) = 0.0_R8
         aa9(ix,iy, 2) = 0.0_R8
        else
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solvereg(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),1) = -conc(ix,iy,0,0)
         if (solvereg(ireg))
     >    aa9(ix,iy,0) = aa9(ix,iy,0)+conc(ix,iy,0,0)
         if (solvereg(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),0) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy),0)+conc(ix,iy,0,0)
         if (solvereg(ireg))
     >    aa9(ix,iy,-1) = -conc(ix,iy,0,0)
cwdk add 9-point contributions
         if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
          t0  = -conc(ix,iy,0,1)/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          t1  = -conc(ix,iy,0,1)/
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
        endif                                                        !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg)) 
     >   aa9(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
          aa9(ix,iy,-3) = 0.0_R8
        else
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot))
     >     aa9(bottomix(ix,iy), bottomiy(ix,iy),3) = -conc(ix,iy,1,1)
          if (solvereg(ireg))
     >     aa9(ix,iy,0) = aa9(ix,iy,0)+conc(ix,iy,1,1)
          if (solvereg(ibot))
     >     aa9(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &       aa9(bottomix(ix,iy),bottomiy(ix,iy),0)+conc(ix,iy,1,1)
          if (solvereg(ireg))
     >     aa9(ix,iy,-3) = -conc(ix,iy,1,1)
cwdk add 9-point contributions
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           t0 = -conc(ix,iy,1,0)/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
           t1 = -conc(ix,iy,1,0)/
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
cxpb
          if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &       topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &	     iy.ne.ny .and. aa9(ix,iy,3).eq.0.0_R8) then
           if(ncall.lt.3) write(*,*) 'Executing xpb:b2uspo for ',ix,iy
           itop = region(topix(ix,iy),topiy(ix,iy),0)
           if (solvereg(ireg))
     >      aa9(ix,iy,3) = -conc(topix(ix,iy),topiy(ix,iy),1,1)
           if (solvereg(itop))
     >      aa9(topix(ix,iy),topiy(ix,iy),0) = 
     &      aa9(topix(ix,iy),topiy(ix,iy),0)+
     &      conc(topix(ix,iy),topiy(ix,iy),1,1)
           if (solvereg(ireg))
     >      aa9(ix,iy,0) = 
     &      aa9(ix,iy,0)+conc(topix(ix,iy),topiy(ix,iy),1,1)
           if (solvereg(itop))
     >      aa9(topix(ix,iy),topiy(ix,iy),3) =
     &      -conc(topix(ix,iy),topiy(ix,iy),1,1) 
          endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solvereg(ireg)) 
     >   aa9(ix,iy,3) = 0.0_R8
       enddo
      enddo
*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
          aa9(ix,iy,0) = 1.0_R8
          respo(ix,iy) = 0.0_R8
        end if
        if (.not.solvereg(ireg)) respo(ix,iy) = 0.0_R8
       enddo
      enddo

*   ..solve the correction equation
      write(*,*) 'Calling b2ux9p from '//trim(name)
      call b2ux9p (nx, ny, aa9, itcnt, respo, corpo, name)           !srv 18.03.10

      if (iout.ne.0) then                                          !srv 08.07.08 {
        do k=0,3
          write (chark,'(i1)') k
          call my_out(70,nx,ny,sch(-1,-1,k),'b2uspo_sch'//chark)
        enddo
        call my_out(70,nx,ny,aa9(-1,-1,-6),'b2uspo_aa9m6')
        call my_out(70,nx,ny,aa9(-1,-1,-5),'b2uspo_aa9m5')
        call my_out(70,nx,ny,aa9(-1,-1,-4),'b2uspo_aa9m4')
        call my_out(70,nx,ny,aa9(-1,-1,-3),'b2uspo_aa9m3')
        call my_out(70,nx,ny,aa9(-1,-1,-2),'b2uspo_aa9m2')
        call my_out(70,nx,ny,aa9(-1,-1,-1),'b2uspo_aa9m1')
        call my_out(70,nx,ny,aa9(-1,-1, 0),'b2uspo_aa9m0')
        call my_out(70,nx,ny,aa9(-1,-1, 1),'b2uspo_aa9p1')
        call my_out(70,nx,ny,aa9(-1,-1, 2),'b2uspo_aa9p2')
        call my_out(70,nx,ny,aa9(-1,-1, 3),'b2uspo_aa9p3')
        call my_out(70,nx,ny,aa9(-1,-1, 4),'b2uspo_aa9p4')
        call my_out(70,nx,ny,aa9(-1,-1, 5),'b2uspo_aa9p5')
        call my_out(70,nx,ny,aa9(-1,-1, 6),'b2uspo_aa9p6')
        call my_out(70,nx,ny,respo,'b2uspo_respo')
        call my_out(70,nx,ny,corpo,'b2uspo_corpo')
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uspo

      end
