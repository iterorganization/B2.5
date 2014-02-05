*-----------------------------------------------------------------------
*.specification

      subroutine b2sqcx (nx, ny, ns, ev, amh0, th0, rlcx0, rcx, wrk0)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   ev, amh0, th0(-1:nx,-1:ny), rlcx0(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   rcx(-1:nx,-1:ny,0:ns-1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SQCX computes on all cells of the mesh the rate coefficients
*     for charge exchange with one neutral hydrogen species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine employs the log-log linear representation provided
*     in the array rlcx0 to evaluate on all cells of the mesh the rate
*     coefficients for charge exchange with one neutral hydrogen
*     species, characterised by atomic mass amh0.
*
*
*  8. further comments
*
*     In case a rate coefficient vanishes, the representation of the
*     log of the rate coefficient will provide a large negative value,
*     for which the exponential function should underflow. The present
*     code tests for this condition and explicitly sets those rate
*     coefficients to 0.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, styl0, n2, ix, iy, is
      real (kind=R8) ::
     *   phm0, cutlo, cutll, t0
      save ncall, styl0, phm0
*   ..procedures
      intrinsic log, exp, max, dim
      real (kind=R8) ::
     *   machsfr, expu
      external subini, subend, xertst, ipgeti, ipgetr, machsfr
      external b2xvsg
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
*     (expu is the exponential function with underflow set to 0)
*   ..initialisation
      data ncall/0/, styl0/0/, phm0/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sqcx')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sqcx_styl0', styl0)
       call xertst (0.le.styl0, 'faulty internal parameter styl0')
       call ipgetr ('b2sqcx_phm0', phm0)
       call xertst (0.0_R8.le.phm0, 'faulty internal parameter phm0')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev, amh0
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
      call xertst (0.0_R8.lt.amh0.and.amh0.le.3.0_R8,
     &  'faulty argument range amh0')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test th0
       call b2xvsg (n2, th0, 1, 'th0', '.gt.')
      endif
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
*   ..compute wrk0(,) = log(th0(,)/(amh0*ev))
      do iy = -1, ny
       do ix = -1, nx
        wrk0(ix,iy) = log(th0(ix,iy)/(amh0*ev))
       enddo
      enddo

* ..main computation
*   ..compute rate coefficients
      do is = 0, ns-1
      if (styl0.eq.0.or.is.lt.2) then
       do iy = -1, ny
        do ix = -1, nx
*      ..compute charge exchange rate coefficient
         rcx(ix,iy,is) = phm0*
     &     expu(rlcx0(ix,iy,0,is)+rlcx0(ix,iy,1,is)*wrk0(ix,iy))
        enddo
       enddo
      else
       do iy = -1, ny
        do ix = -1, nx
*      ..compute charge exchange rate coefficient
         rcx(ix,iy,is) = 0.0_R8
        enddo
       enddo
      endif
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sqcx

      end subroutine b2sqcx
