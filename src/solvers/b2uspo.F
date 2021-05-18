*-----------------------------------------------------------------------
*.specification

      subroutine b2uspo (nx, ny, nregionv, solvereg, itcnt,
     &  ne, te, conc, sch, respo, corpo, aa5, name)                  !srv 18.03.10
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
     &  conc(-1:nx,-1:ny,0:1),
     &  sch(-1:nx,-1:ny,0:3), respo(-1:nx,-1:ny)
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpo(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa5(-1:nx,-1:ny,-2:2)
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
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop
      !$OMP THREADPRIVATE(ncall)
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
      call sfill (n2*5, 0.0_R8, aa5, 1)
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         aa5(ix,iy,0) = -sch(ix,iy,1)-sch(ix,iy,3)*ne(ix,iy)
        else
         aa5(ix,iy,0) = 1.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa5(ix,iy,-1) = 0.0_R8
        else
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solvereg(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),1) = -conc(ix,iy,0)
         if (solvereg(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)+conc(ix,iy,0)
         if (solvereg(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),0) =
     &       aa5(leftix(ix,iy),leftiy(ix,iy),0)+conc(ix,iy,0)
         if (solvereg(ireg))
     >    aa5(ix,iy,-1) = -conc(ix,iy,0)
        endif                                                        !srv 26.11.02
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg))
     >   aa5(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
          aa5(ix,iy,-2) = 0.0_R8
        else
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy), bottomiy(ix,iy),2) = -conc(ix,iy,1)
          if (solvereg(ireg))
     >     aa5(ix,iy,0) = aa5(ix,iy,0)+conc(ix,iy,1)
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &       aa5(bottomix(ix,iy),bottomiy(ix,iy),0)+conc(ix,iy,1)
          if (solvereg(ireg))
     >     aa5(ix,iy,-2) = -conc(ix,iy,1)
cxpb
          if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &       topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &       iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
           if(ncall.lt.3) write(*,*) 'Executing xpb:b2uspo for ',ix,iy
           itop = region(topix(ix,iy),topiy(ix,iy),0)
           if (solvereg(ireg))
     >      aa5(ix,iy,2) = -conc(topix(ix,iy),topiy(ix,iy),1)
           if (solvereg(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),0) =
     &      aa5(topix(ix,iy),topiy(ix,iy),0)+
     &      conc(topix(ix,iy),topiy(ix,iy),1)
           if (solvereg(ireg))
     >      aa5(ix,iy,0) =
     &      aa5(ix,iy,0)+conc(topix(ix,iy),topiy(ix,iy),1)
           if (solvereg(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),2) =
     &      -conc(topix(ix,iy),topiy(ix,iy),1)
          endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solvereg(ireg))
     >   aa5(ix,iy,2) = 0.0_R8
       enddo
      enddo
*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (aa5(ix,iy,0).eq.0.0_R8 .and.
     &      aa5(ix,iy,1).eq.0.0_R8 .and. aa5(ix,iy,-1).eq.0.0_R8 .and.
     &      aa5(ix,iy,2).eq.0.0_R8 .and. aa5(ix,iy,-2).eq.0.0_R8) then
          aa5(ix,iy,0) = 1.0_R8
          respo(ix,iy) = 0.0_R8
        end if
        if (.not.solvereg(ireg)) respo(ix,iy) = 0.0_R8
       enddo
      enddo

*   ..solve the correction equation
      write(*,'(1x,a,a)') 'Calling b2ux5p from '//trim(name)
      call b2ux5p (nx, ny, aa5, itcnt, respo, corpo, name)           !srv 18.03.10

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2uspo

      end subroutine b2uspo
