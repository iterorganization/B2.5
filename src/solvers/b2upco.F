*-----------------------------------------------------------------------
*.specification

      subroutine b2upco (nx, ny, nregionv, ua_solve, solvereg, 
     &  rxf, pcm0, vol, corpb, pccb0, pb, nb, ub)
      use b2mod_indirect
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nregionv
      real (kind=R8) ::
     *  rxf, pcm0, vol(-1:nx,-1:ny), corpb(-1:nx,-1:ny),
     *  pccb0(-1:nx,-1:ny)
      logical ua_solve, solvereg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  pb(-1:nx,-1:ny), nb(-1:nx,-1:ny), ub(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPCO corrects the state of the plasma towards the solution of
*     continuity equation for one atomic species.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy, ireg
!      integer, save :: mdf_fnb = 1                                       !srv 13.10.06
      save ncall
*     ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
      real (kind=R8) ::
     *  t0
*   ..procedures
      external subini, subend, xertst, b2xxgs
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upco')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8, 
     .  'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')

* ..implement the correction
*   ..compute wrk0
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         wrk0(ix,iy) = 0.0_R8
        else
         wrk0(ix,iy) = corpb(ix,iy)-corpb(leftix(ix,iy),leftiy(ix,iy))
        endif
       enddo
      enddo
      call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
*   ..compute new pb and nb
      do iy = -1, ny
       do ix = -1, nx
        if (0.0_R8.le.corpb(ix,iy)) then
         t0 = 1.0_R8+rxf*corpb(ix,iy)/pb(ix,iy)
        else
         t0 = 1.0_R8/(1.0_R8-rxf*corpb(ix,iy)/pb(ix,iy))
        endif
        nb(ix,iy) = t0*nb(ix,iy)
        pb(ix,iy) = t0*pb(ix,iy)
       enddo
      enddo
*   ..compute new ub
      if (ua_solve) then
       do iy = -1, ny
        do ix = -1, nx
*       (note the further damping on the velocity correction)
         ireg = region(ix,iy,0)
         if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &      solvereg(ireg))
     &    ub(ix,iy) = ub(ix,iy)-pcm0*rxf*pccb0(ix,iy)*vol(ix,iy)*
     &    (wrk0(ix,iy)/(vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))+ 
     &     wrk0(rightix(ix,iy),rightiy(ix,iy))/(vol(ix,iy)+
     &      vol(rightix(ix,iy),rightiy(ix,iy))))
        enddo
       enddo
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upco

      end
