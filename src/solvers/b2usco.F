*-----------------------------------------------------------------------
*.specification

      subroutine b2usco (nx, ny, isb, nregionv, solvereg, itcnt, rxg, !srv 06.07.08
     &  pcm1, vol, pbs, nb, pb,
     &  snbv, flob, conb, pccb1, rescb, corpb, aa5, name)             !srv 18.05.02
      use b2mod_types
      use b2mod_indirect
      use b2mod_boundary_namelist
     , , only : iffboundx, iffboundy
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, nregionv, itcnt
      real (kind=R8) ::
     *  rxg, pcm1, pbs(-1:nx,-1:ny,0:1), nb(-1:nx,-1:ny),
     &  pb(-1:nx,-1:ny), snbv(-1:nx,-1:ny), flob(-1:nx,-1:ny,0:1,0:1),
     &  conb(-1:nx,-1:ny,0:1,0:1,0:4), pccb1(-1:nx,-1:ny),            !srv 13.10.06
     &  rescb(-1:nx,-1:ny), vol(-1:nx,-1:ny)
      logical solvereg(0:nregionv)      
      character*(*) name                                              !srv 18.03.10
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpb(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa5(-1:nx,-1:ny,-2:2), wrk0(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USCO computes a correction towards the solution of the
*     continuity equation for one atomic species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !srv 06.07.06
      integer, save :: style_b2mndt = 1                               !srv 07.07.21
      integer, save :: iout = 0                                       !srv 17.06.02
      character chns*3                                                !srv 17.06.02
      real (kind=R8) ::
     *  t0, t1
      save ncall
#ifdef _OPENMP
!$OMP THREADPRIVATE(ncall, no_Ptncr_x_co, no_Ptncr_y_co, iout)
#endif
*   ..procedures
      intrinsic abs
      external xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2xxgs, b2ux5p
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usco')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                            !srv 06.07.08 {
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)
       call ipgeti ('b2mndt_style', style_b2mndt)                     !srv 07.07.21
       call ipgeti ('b2usco_iout', iout)                              !srv 17.06.02
       if (in_parallel() .and. iout.ne.0) then
#ifdef _OPENMP
!$OMP MASTER
#endif
         write(*,*) 'b2usco OpenMP warning: no file IO in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
         iout = 0
       endif
      endif                                                           !srv 06.07.08 }
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt, rxg, pcm1
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.0_R8.lt.rxg, 'faulty argument rxg')
      call xertst (0.0_R8.le.pcm1, 'faulty argument pcm1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, nb, pb
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
*    ..test sign of snbv, conb
       call b2xvsg (n2, snbv, 1, 'snbv', '.le.')
       if(no_Ptncr_x_co.eq.0 .and. no_Ptncr_y_co.eq.0) then           !srv 06.07.06 {
        call b2xvsg (n2, conb(-1,-1,0,0,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,0), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,1), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,0,0,2), 1, 'conb', '.ge.')
        call b2xvsg (n2, conb(-1,-1,1,1,2), 1, 'conb', '.ge.')
       endif                                                          !srv 06.07.06 }
*    ..test edge values of flob, conb
       call b2xvff (nx, ny, flob, 'flob')
       call b2xvff (nx, ny, conb(-1,-1,0,0,0), 'conb0')               !srv 13.10.06 {
       call b2xvff (nx, ny, conb(-1,-1,0,0,1), 'conb1')
       call b2xvff (nx, ny, conb(-1,-1,0,0,2), 'conb2')               !srv 13.10.06 }
      endif

* ..compute the correction
      call sfill (n2*5, 0.0_R8, aa5, 1)
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
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         aa5(ix,iy,0) = -snbv(ix,iy)*nb(ix,iy)/pb(ix,iy)
         if (style_b2mndt.lt.2) aa5(ix,iy,0) = aa5(ix,iy,0)+
     &    abs(rescb(ix,iy))/(rxg*pb(ix,iy))
        else
         aa5(ix,iy,0) = 1.0_R8
        endif
        aa5(ix,iy,0) = aa5(ix,iy,0)*rscal(pb(ix,iy))
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
          t0 =  flob(ix,iy,0,0)/2.0_R8 -
     &         (conb(ix,iy,0,0,0)+pcm1*wrk0(ix,iy)) - conb(ix,iy,0,0,1) !srv 13.10.06
          t1 = -flob(ix,iy,0,0)/2.0_R8 -
     &         (conb(ix,iy,0,0,0)+pcm1*wrk0(ix,iy)) - conb(ix,iy,0,0,2) !srv 13.10.06
          ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
          if (solvereg(ileft))
     >     aa5(leftix(ix,iy),leftiy(ix,iy),1) = t0*rscal(pb(ix,iy))
          if (solvereg(ireg))
     >     aa5(ix,iy,0) = aa5(ix,iy,0)-t0*rscal(pb(ix,iy))
          if (solvereg(ileft))
     >     aa5(leftix(ix,iy),leftiy(ix,iy),0) = 
     &       aa5(leftix(ix,iy),leftiy(ix,iy),0)
     &       -t1*rscal(pb(leftix(ix,iy),leftiy(ix,iy)))
          if (solvereg(ireg))
     >     aa5(ix,iy,-1) = t1*rscal(pb(leftix(ix,iy),leftiy(ix,iy)))
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg)) 
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
          t0 =  flob(ix,iy,1,1)/2.0_R8 - conb(ix,iy,1,1,0)
     &                                 - conb(ix,iy,1,1,1)             !srv 13.10.06
          t1 = -flob(ix,iy,1,1)/2.0_R8 - conb(ix,iy,1,1,0)
     &                                 - conb(ix,iy,1,1,2)             !srv 13.10.06
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy),bottomiy(ix,iy),2) = t0*rscal(pb(ix,iy))
          if (solvereg(ireg))
     >     aa5(ix,iy,0) = aa5(ix,iy,0)-t0*rscal(pb(ix,iy))
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &      aa5(bottomix(ix,iy),bottomiy(ix,iy),0)
     &      -t1*rscal(pb(bottomix(ix,iy),bottomiy(ix,iy)))
          if (solvereg(ireg))
     >     aa5(ix,iy,-2) = t1*rscal(pb(bottomix(ix,iy),bottomiy(ix,iy)))
cxpb
          if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &       topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &       iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
            if(ncall.lt.3) write(*,*) 'Executing xpb:b2usco for ',ix,iy
            t0 = flob(topix(ix,iy),topiy(ix,iy),1,1)/2-
     &           conb(topix(ix,iy),topiy(ix,iy),1,1,0)-                !srv 13.10.06
     &           conb(topix(ix,iy),topiy(ix,iy),1,1,1)                 !srv 13.10.06
            t1 =-flob(topix(ix,iy),topiy(ix,iy),1,1)/2-
     &           conb(topix(ix,iy),topiy(ix,iy),1,1,0)-                !srv 13.10.06
     &           conb(topix(ix,iy),topiy(ix,iy),1,1,2)                 !srv 13.10.06
            itop = region(topix(ix,iy),topiy(ix,iy),0)
            if (solvereg(ireg))
     >       aa5(ix,iy,2) = t0*rscal(pb(topix(ix,iy),topiy(ix,iy)))
            if (solvereg(itop))
     >       aa5(topix(ix,iy),topiy(ix,iy),0) = 
     &       aa5(topix(ix,iy),topiy(ix,iy),0)-
     &       t0*rscal(pb(topix(ix,iy),topiy(ix,iy)))
            if (solvereg(ireg))
     >       aa5(ix,iy,0) = aa5(ix,iy,0)-t1*rscal(pb(ix,iy))
            if (solvereg(itop))
     >       aa5(topix(ix,iy),topiy(ix,iy),2) = t1*rscal(pb(ix,iy))
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
          aa5(ix,iy,0) = rscal(pb(ix,iy))
          rescb(ix,iy) = 0.0_R8
        end if
        if (.not.solvereg(ireg)) rescb(ix,iy) = 0.0_R8
       enddo
      enddo
      
*   ..solve the correction equation
      write(*,'(1x,a,a)') 'Calling b2ux5p from ',name
      call b2ux5p (nx, ny, aa5, itcnt, rescb, corpb, name)             !srv 18.03.10
*   ..redo the scaling
      do iy = -1, ny
       do ix = -1, nx
        corpb(ix,iy) = corpb(ix,iy)*rscal(pb(ix,iy))
       enddo
      enddo

      if (iout.ne.0) then                                              !srv 17.06.02 {
       if (in_parallel()) then
        if (ncall==0) then
#ifdef _OPENMP
!$OMP MASTER
#endif
         write(*,*)'OpenMP warning from b2usco: not writing files ',
     &        'in parallel mode'
#ifdef _OPENMP
!$OMP END MASTER
#endif
        end if
       else
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,aa5(-1,-1,-2),name//'_aa5m2'//chns)
        call my_out(70,nx,ny,aa5(-1,-1,-1),name//'_aa5m1'//chns)
        call my_out(70,nx,ny,aa5(-1,-1, 0),name//'_aa5m0'//chns)
        call my_out(70,nx,ny,aa5(-1,-1, 1),name//'_aa5p1'//chns)
        call my_out(70,nx,ny,aa5(-1,-1, 2),name//'_aa5p2'//chns)
        call my_out(70,nx,ny,corpb,name//'_corpb'//chns)
        call my_out(70,nx,ny,rescb,name//'_rescb'//chns)
        call my_out(70,nx,ny,nb,name//'_nb'//chns)
       endif
      endif                                                            !srv 17.06.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usco

      contains

      function rscal(t0)
      implicit none
      real(kind=R8) :: rscal, t0
! if uncomment rscal(t0) = t0 then uncomment nb(ix,iy) = nb(ix,iy)+rxf*corpb(ix,iy)/pb(ix,iy) in b2upco
      rscal = t0                                                  !srv 13.10.06
!      rscal = 1.0e0_R8                                           !srv 13.10.06
      return
      end function rscal

      end subroutine b2usco
