*-----------------------------------------------------------------------
*.specification

      subroutine b2usht (nx, ny, ns, nregionv, itcnt,
     &  rxg, pcm1, pbs, ne, ni, te, ti,
     &  floe, cone, floi, coni, pccm1, she, shi, ceqp,
     &  resht, reshe, reshi, cortt, corte, corti, flo0, con0, aa5, name)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_ilutern
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_boundary_namelist
     , , only : iffboundx, iffboundy
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv, itcnt
      character*(*) name                                                !srv 18.03.10
      real (kind=R8) ::
     *  rxg, pcm1, pbs(-1:nx,-1:ny,0:1), ne(-1:nx,-1:ny),
     &  ni(-1:nx,-1:ny,0:1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  floe(-1:nx,-1:ny,0:1), floi(-1:nx,-1:ny,0:1),
     &  cone(-1:nx,-1:ny,0:1,0:2), coni(-1:nx,-1:ny,0:1,0:2),           !srv 02.01.07
     &  pccm1(-1:nx,-1:ny)
      real (kind=R8) ::
     *  she(-1:nx,-1:ny,0:3), shi(-1:nx,-1:ny,0:3),
     &  ceqp(-1:nx,-1:ny), resht(-1:nx,-1:ny), reshe(-1:nx,-1:ny),
     &  reshi(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cortt(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  flo0(-1:nx,-1:ny,0:1), con0(-1:nx,-1:ny,0:1,0:2),               !srv 02.01.07
     &  aa5(-1:nx,-1:ny,-2:2)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USHT computes a correction to the fields te and ti.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop
      integer, save :: style_b2mndt = 1                                 !srv 07.07.21
      logical new_matrix
      save ncall
      real (kind=R8) ::
     *   t0, t1
      real (kind=R8), save :: cor_diag = 1.0_R8                         !srv 07.07.21
      !$OMP THREADPRIVATE(ncall)
*   ..procedures
      intrinsic abs, sqrt, dim
      external xertst, sfill
      external b2xvsg, b2xvff, b2ux5p
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usht')
      if (ncall.eq.0) then                                              !srv 07.07.21
       call ipgetr ('b2usht_cor_diag', cor_diag)
       call ipgeti ('b2mndt_style', style_b2mndt)
      endif                                                             !srv 07.07.21
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test itcnt, rxg, pcm1
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.0_R8.lt.rxg, 'faulty argument rxg')
      call xertst (0.0_R8.le.pcm1, 'faulty argument pcm1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, ni, te, ti
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
*    ..test sign of cone, coni
       call b2xvsg (n2*2, cone, 1, 'cone', '.ge.')
       call b2xvsg (n2*2, coni, 1, 'coni', '.ge.')
*    ..test edge values of floe, floi, cone, coni
       call b2xvff (nx, ny, floe, 'floe')
       call b2xvff (nx, ny, floi, 'floi')
       call b2xvff (nx, ny, cone(-1,-1,0,0), 'cone0')                   !srv 02.01.07 {
       call b2xvff (nx, ny, cone(-1,-1,0,1), 'cone1')
       call b2xvff (nx, ny, cone(-1,-1,0,2), 'cone2')
       call b2xvff (nx, ny, coni(-1,-1,0,0), 'coni0')
       call b2xvff (nx, ny, coni(-1,-1,0,1), 'coni1')
       call b2xvff (nx, ny, coni(-1,-1,0,2), 'coni2')                   !srv 02.01.07 }
*    ..test sign of she, shi, ceqp
       if(style_b2mndt.ne.2) then                             !srv 07.07.21
         call b2xvsg (n2, she(-1,-1,0), 1, 'she0', '.ge.')
       endif                                                  !srv 07.07.21
       call b2xvsg (n2, she(-1,-1,1), 1, 'she1', '.le.')
       call b2xvsg (n2, she(-1,-1,2), 1, 'she2', '.ge.')
       call b2xvsg (n2, she(-1,-1,3), 1, 'she3', '.le.')
       if(style_b2mndt.ne.2) then                             !srv 07.07.21
         call b2xvsg (n2, shi(-1,-1,0), 1, 'shi0', '.ge.')
       endif                                                  !srv 07.07.21
       call b2xvsg (n2, shi(-1,-1,1), 1, 'shi1', '.le.')
       call b2xvsg (n2, shi(-1,-1,2), 1, 'shi2', '.ge.')
       call b2xvsg (n2, shi(-1,-1,3), 1, 'shi3', '.le.')
       call b2xvsg (n2, ceqp, 1, 'ceqp', '.ge.')
      endif

      if (any(solveee(0:nnreg(0))).and.any(solveei(0:nnreg(0))).and.
     &    any(solveet(0:nnreg(0)))) then                                !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveet(ireg).neqv.last_solve_5(ireg)
      end do
      if (new_matrix) call restart_ma28
* ..total heat correction equation
      call sfill (n2*5, 0.0_R8, aa5, 1)
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0_R8
         con0(ix,iy,0,:) = 0.0_R8                                       !srv 02.01.07
        else
         t0 = 2.5_R8*pcm1*pbs(ix,iy,0)*
     &      (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &       (ne(leftix(ix,iy),leftiy(ix,iy))*
     &        te(leftix(ix,iy),leftiy(ix,iy))+
     &        ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &        ti(leftix(ix,iy),leftiy(ix,iy)))+
     &       pccm1(ix,iy)*(ne(ix,iy)*te(ix,iy)+ni(ix,iy,0)*ti(ix,iy)))/
     &      2.0_R8
         flo0(ix,iy,0) = floe(ix,iy,0)+floi(ix,iy,0)-
     &    t0*(ne(ix,iy)+ni(ix,iy,0)-ne(leftix(ix,iy),leftiy(ix,iy))
     &             -ni(leftix(ix,iy),leftiy(ix,iy),0))
         con0(ix,iy,0,0) = cone(ix,iy,0,0)+coni(ix,iy,0,0)+             !srv 02.01.07
     &    (t0/2.0_R8)*(ne(leftix(ix,iy),leftiy(ix,iy))
     &     +ne(ix,iy)+ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))
        endif
       enddo
      enddo
*   ..compute flo0y, con0y
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         flo0(ix,iy,1) = 0.0_R8
         con0(ix,iy,1,:) = 0.0_R8                                       !srv 02.01.07
        else
         flo0(ix,iy,1) = floe(ix,iy,1)+floi(ix,iy,1)
         con0(ix,iy,1,0) = cone(ix,iy,1,0)+coni(ix,iy,1,0)              !srv 02.01.07
        endif
       enddo
      enddo
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveet(ireg)) then
         aa5(ix,iy,0) =
     &     -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)-
     &      shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)+
     &    abs(resht(ix,iy))/(rxg*min(te(ix,iy),ti(ix,iy)))
        else
         aa5(ix,iy,0) = 1.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundx(ix,iy).ne.0) then
          aa5(ix,iy,-1) = 0.0_R8
          aa5(ix,iy,1) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa5(ix,iy,-1) = 0.0_R8
        else
         t0 =  flo0(ix,iy,0)/2.0_R8-con0(ix,iy,0,0)                     !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2.0_R8-con0(ix,iy,0,0)                     !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveet(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),1) = t0
         if (solveet(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0
         if (solveet(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),0) = 
     &      aa5(leftix(ix,iy),leftiy(ix,iy),0)-t1
         if (solveet(ireg))
     >    aa5(ix,iy,-1) = t1
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solveet(ireg))
     >   aa5(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundy(ix,iy).ne.0) then
          aa5(ix,iy,-2) = 0.0_R8
          aa5(ix,iy,2) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
         aa5(ix,iy,-2) = 0.0_R8
        else
         t0 =  flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         t1 = -flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
         if (solveet(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),2) = t0
         if (solveet(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0
         if (solveet(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &     aa5(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
         if (solveet(ireg))
     >    aa5(ix,iy,-2) = t1
cxpb
         if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &      topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &      iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
           if(ncall.lt.3) 
     >      write(*,*) 'Executing xpb:b2usht(1) for ',ix,iy
           itop = region(topix(ix,iy),topiy(ix,iy),0)
           t0 =  flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           t1 = -flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           if (solveet(ireg))
     >      aa5(ix,iy,2) = t0
           if (solveet(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),0) = 
     &      aa5(topix(ix,iy),topiy(ix,iy),0)-t0
           if (solveet(ireg))
     >      aa5(ix,iy,0) = aa5(ix,iy,0)-t1
           if (solveet(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),2) = t1
         endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solveet(ireg))
     >   aa5(ix,iy,2) = 0.0_R8
       enddo
      enddo
*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        if (aa5(ix,iy,0).eq.0.0_R8 .and.
     &      aa5(ix,iy,1).eq.0.0_R8 .and. aa5(ix,iy,-1).eq.0.0_R8 .and.
     &      aa5(ix,iy,2).eq.0.0_R8 .and. aa5(ix,iy,-2).eq.0.0_R8) then
          aa5(ix,iy,0) = 1.0_R8
          resht(ix,iy) = 0.0_R8
        end if
       enddo
      enddo

*   ..solve the correction equation
        write(*,'(1x,a,a)') 'Calling b2ux5p from ', name//'tt'
        call b2ux5p (nx, ny, aa5, itcnt, resht, cortt, name//'tt')      !srv 18.03.10
      endif                                                             !srv 22.05.18

* ..electron temperature correction equation
      if (any(solveee(0:nnreg(0)))) then                                !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveee(ireg).neqv.solveet(ireg)
      end do
      if (new_matrix) call restart_ma28
      call sfill (n2*5, 0.0_R8, aa5, 1)
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0_R8
         con0(ix,iy,0,:) = 0.0_R8                                       !srv 02.01.07
        else
         t0 = 2.5_R8*pcm1*pbs(ix,iy,0)*
     &      (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &       ne(leftix(ix,iy),leftiy(ix,iy))*
     &       te(leftix(ix,iy),leftiy(ix,iy))+
     &       pccm1(ix,iy)*ne(ix,iy)*te(ix,iy))/2.0_R8
         flo0(ix,iy,0) = floe(ix,iy,0)
     &     -t0*(ne(ix,iy)-ne(leftix(ix,iy),leftiy(ix,iy)))
         con0(ix,iy,0,0) = cone(ix,iy,0,0)+(t0/2.0_R8)*                 !srv 02.01.07
     &     (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))
        endif
       enddo
      enddo
*   ..compute flo0y, con0y
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         flo0(ix,iy,1) = 0.0_R8
         con0(ix,iy,1,:) = 0.0_R8                                       !srv 02.01.07
        else
         flo0(ix,iy,1) = floe(ix,iy,1)
         con0(ix,iy,1,0) = cone(ix,iy,1,0)                              !srv 02.01.07
        endif
       enddo
      enddo
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveee(ireg)) then
         aa5(ix,iy,0) = -she(ix,iy,1)-she(ix,iy,3)*ne(ix,iy)+
     &    cor_diag*dim(1.0_R8,1.5_R8*(te(ix,iy)-ti(ix,iy))/te(ix,iy))*   !srv 07.07.21
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))+
     &    abs(reshe(ix,iy))/(rxg*te(ix,iy))
        else
         aa5(ix,iy,0) = 1.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundx(ix,iy).ne.0) then
          aa5(ix,iy,-1) = 0.0_R8
          aa5(ix,iy,1) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa5(ix,iy,-1) = 0.0_R8
        else
         t0 =  flo0(ix,iy,0)/2.0_R8-con0(ix,iy,0,0)                     !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2.0_R8-con0(ix,iy,0,0)                     !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveee(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),1) = t0
         if (solveee(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0
         if (solveee(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),0) = 
     &       aa5(leftix(ix,iy),leftiy(ix,iy),0)-t1
         if (solveee(ireg))
     >    aa5(ix,iy,-1) = t1
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solveee(ireg))
     >   aa5(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundy(ix,iy).ne.0) then
          aa5(ix,iy,-2) = 0.0_R8
          aa5(ix,iy,2) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
         aa5(ix,iy,-2) = 0.0_R8
        else
         t0 =  flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         t1 = -flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
         if (solveee(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),2) = t0
         if (solveee(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0
         if (solveee(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &       aa5(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
         if (solveee(ireg))
     >    aa5(ix,iy,-2) = t1
cxpb
         if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &      topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &      iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
           if(ncall.lt.3) 
     >      write(*,*) 'Executing xpb:b2usht(2) for ',ix,iy
           itop = region(topix(ix,iy),topiy(ix,iy),0)
           t0 =  flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           t1 = -flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           if (solveee(ireg))
     >      aa5(ix,iy,2) = t0
           if (solveee(ireg))
     >      aa5(topix(ix,iy),topiy(ix,iy),0) = 
     &      aa5(topix(ix,iy),topiy(ix,iy),0)-t0
           if (solveee(ireg))
     >      aa5(ix,iy,0) = aa5(ix,iy,0)-t1
           if (solveee(ireg))
     >      aa5(topix(ix,iy),topiy(ix,iy),2) = t1
         endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solveee(ireg))
     >   aa5(ix,iy,2) = 0.0_R8
       enddo
      enddo
*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        if (aa5(ix,iy,0).eq.0.0_R8 .and.
     &      aa5(ix,iy,1).eq.0.0_R8 .and. aa5(ix,iy,-1).eq.0.0_R8 .and.
     &      aa5(ix,iy,2).eq.0.0_R8 .and. aa5(ix,iy,-2).eq.0.0_R8) then
          aa5(ix,iy,0) = 1.0_R8
          reshe(ix,iy) = 0.0_R8
        end if
       enddo
      enddo

*   ..solve the correction equation
        write(*,'(1x,a,a)') 'Calling b2ux5p from ', name//'te'
        call b2ux5p (nx, ny, aa5, itcnt, reshe, corte, name//'te')      !srv 18.03.10
      endif                                                             !srv 22.05.18

* ..atom temperature correction equation
      if (any(solveei(0:nnreg(0)))) then                                !srv 22.05.18
      new_matrix = .false.
      do ireg = 0, nregionv
        new_matrix = new_matrix.or.
     &   solveei(ireg).neqv.solveee(ireg)
      end do
      if (new_matrix) call restart_ma28
      call sfill (n2*5, 0.0_R8, aa5, 1)
*   ..compute flo0x, con0x
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         flo0(ix,iy,0) = 0.0_R8
         con0(ix,iy,0,:) = 0.0_R8                                       !srv 02.01.07
        else
         t0 = 2.5_R8*pcm1*pbs(ix,iy,0)*
     &      (pccm1(leftix(ix,iy),leftiy(ix,iy))*
     &       ni(leftix(ix,iy),leftiy(ix,iy),0)*
     &       ti(leftix(ix,iy),leftiy(ix,iy))+
     &       pccm1(ix,iy)*ni(ix,iy,0)*ti(ix,iy))/2.0_R8
         flo0(ix,iy,0) = floi(ix,iy,0)-t0*(ni(ix,iy,0)
     &     -ni(leftix(ix,iy),leftiy(ix,iy),0))
         con0(ix,iy,0,0) = coni(ix,iy,0,0)+(t0/2.0_R8)*                 !srv 02.01.07
     &      (ni(leftix(ix,iy),leftiy(ix,iy),0)+ni(ix,iy,0))
        endif
       enddo
      enddo
*   ..compute flo0y, con0y
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         flo0(ix,iy,1) = 0.0_R8
         con0(ix,iy,1,:) = 0.0_R8                                       !srv 02.01.07
        else
         flo0(ix,iy,1) = floi(ix,iy,1)
         con0(ix,iy,1,0) = coni(ix,iy,1,0)                              !srv 02.01.07
        endif
       enddo
      enddo
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solveei(ireg)) then
         aa5(ix,iy,0) = (-shi(ix,iy,1)-shi(ix,iy,3)*ni(ix,iy,0)+
     &    ceqp(ix,iy)*ne(ix,iy)*ni(ix,iy,0)/(te(ix,iy)*sqrt(te(ix,iy)))+
     &    abs(reshi(ix,iy))/(rxg*ti(ix,iy)))*scalr(ti(ix,iy))
        else
         aa5(ix,iy,0) = scalr(ti(ix,iy))
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundx(ix,iy).ne.0) then
          aa5(ix,iy,-1) = 0.0_R8
          aa5(ix,iy,1) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa5(ix,iy,-1) = 0.0_R8
        else
         t0 = flo0(ix,iy,0)/2-con0(ix,iy,0,0)                           !srv 02.01.07
         t1 = -flo0(ix,iy,0)/2-con0(ix,iy,0,0)                          !srv 02.01.07
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solveei(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),1) = t0*scalr(ti(ix,iy))
         if (solveei(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0*scalr(ti(ix,iy))
         if (solveei(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),0) = 
     &    aa5(leftix(ix,iy),leftiy(ix,iy),0)
     &    -t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
         if (solveei(ireg))
     >    aa5(ix,iy,-1) = t1*scalr(ti(leftix(ix,iy),leftiy(ix,iy)))
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solveei(ireg))
     >   aa5(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundy(ix,iy).ne.0) then
          aa5(ix,iy,-2) = 0.0_R8
          aa5(ix,iy,2) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
         aa5(ix,iy,-2) = 0.0_R8
        else
         t0 =  flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         t1 = -flo0(ix,iy,1)/2.0_R8-con0(ix,iy,1,0)                     !srv 02.01.07
         ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
         if (solveei(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),2) = t0*scalr(ti(ix,iy))
         if (solveei(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0*scalr(ti(ix,iy))
         if (solveei(ibot))
     >    aa5(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &    aa5(bottomix(ix,iy),bottomiy(ix,iy),0)
     &    -t1*scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
         if (solveei(ireg))
     >    aa5(ix,iy,-2) = t1*scalr(ti(bottomix(ix,iy),bottomiy(ix,iy)))
cxpb
         if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &      topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &      iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
           if(ncall.lt.3) 
     >      write(*,*) 'Executing xpb:b2usht(3) for ',ix,iy
           itop = region(topix(ix,iy),topiy(ix,iy),0)
           t0 =  flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           t1 = -flo0(topix(ix,iy),topiy(ix,iy),1)/2.0_R8
     &          -con0(topix(ix,iy),topiy(ix,iy),1,0)                    !srv 02.01.07
           if (solveei(ireg))
     >      aa5(ix,iy,2) = t0*scalr(ti(topix(ix,iy),topiy(ix,iy)))
           if (solveei(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),0) = 
     &      aa5(topix(ix,iy),topiy(ix,iy),0)-
     &      t0*scalr(ti(topix(ix,iy),topiy(ix,iy)))
           if (solveei(ireg))
     >      aa5(ix,iy,0) = aa5(ix,iy,0)-t1*scalr(ti(ix,iy))
           if (solveei(itop))
     >      aa5(topix(ix,iy),topiy(ix,iy),2) = t1*scalr(ti(ix,iy))
         endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solveei(ireg))
     >   aa5(ix,iy,2) = 0.0_R8
       enddo
      enddo
*   ..correct RHS if not solving
      do iy = -1, ny
       do ix = -1, nx
        if (aa5(ix,iy,0).eq.0.0_R8 .and.
     &      aa5(ix,iy,1).eq.0.0_R8 .and. aa5(ix,iy,-1).eq.0.0_R8 .and.
     &      aa5(ix,iy,2).eq.0.0_R8 .and. aa5(ix,iy,-2).eq.0.0_R8) then
          aa5(ix,iy,0) = 1.0_R8
          reshi(ix,iy) = 0.0_R8
        end if
       enddo
      enddo

*   ..solve the correction equation
      write(*,'(1x,a,a)') 'Calling b2ux5p from ', name//'ti'
      call b2ux5p (nx, ny, aa5, itcnt, reshi, corti, name//'ti')        !srv 18.03.10
      last_solve_5(0:nnreg(0)) = solveei(0:nnreg(0))
      do iy = -1, ny
        do ix = -1, nx
          corti(ix,iy) = corti(ix,iy)*scalr(ti(ix,iy))
        enddo
      enddo
      endif                                                             !srv 22.05.18

*   ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usht

      contains

      function scalr(t0)
      implicit none
      real(kind=R8) :: scalr, t0
      scalr = t0
      return
      end function scalr

      end subroutine b2usht
