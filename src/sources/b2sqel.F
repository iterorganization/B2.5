*-----------------------------------------------------------------------
*.specification

      subroutine b2sqel (nx, ny, ns, ismain, ev, te, 
     &  rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi,
     &  rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi)
      use b2mod_types
      use b2mod_math
      use b2mod_b2cmpa
#ifdef WITH_FTIMINGS
      use timings
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      real (kind=R8) ::
     *  ev, te(-1:nx,-1:ny), rlsa(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlra(-1:nx,-1:ny,0:1,0:ns-1), rlqa(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlrd(-1:nx,-1:ny,0:1,0:ns-1), rlbr(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlza(-1:nx,-1:ny,0:1,0:ns-1), rlz2(-1:nx,-1:ny,0:1,0:ns-1),
     &  rlpt(-1:nx,-1:ny,0:1,0:ns-1), rlpi(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  rsa(-1:nx,-1:ny,0:ns-1), rra(-1:nx,-1:ny,0:ns-1),
     &  rqa(-1:nx,-1:ny,0:ns-1), rrd(-1:nx,-1:ny,0:ns-1),
     &  rbr(-1:nx,-1:ny,0:ns-1),
     &  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1),
     &  rpt(-1:nx,-1:ny,0:ns-1), rpi(-1:nx,-1:ny,0:ns-1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SQEL computes on all cells of the mesh the rate coefficients
*     for ionisation, recombination, electron heat loss and radiation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine employs the log-log linear representation provided
*     in the arrays rlsa, rlra, rlqa, rlrd and rlbr to evaluate on all cells 
*     of the mesh the rate coefficients for ionisation, recombination, 
*     electron heat loss, and radiation.
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
      integer ncall, n2, ix, iy, is
      real (kind=R8) ::
     *   phm0, phm1, phm2, art_rad, tnorm
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
      save ncall, phm0, phm1, phm2, art_rad
      integer, save :: use_eirene=0
*   ..procedures
      intrinsic log, min
      external subini, subend, xertst, ipgetr, machsfr
*   ..initialisation
      data ncall/0/, phm0/1.0_R8/, phm1/1.0_R8/, phm2/1.0_R8/,
     1 art_rad/0.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sqel')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2sqel_phm0', phm0)
       call ipgetr ('b2sqel_phm1', phm1)
       call ipgetr ('b2sqel_phm2', phm2)
       call xertst (0.0_R8.le.min(phm0,phm1,phm2),
     &  'faulty internal parameter phm0, phm1, phm2')
       call ipgetr ('b2sqel_artificial_radiation', art_rad)
       call ipgeti ('b2mndr_eirene', use_eirene)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test ev
      call xertst (0.0_R8.lt.ev, 'faulty argument ev')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif
*   ..set cutll, cutlo
      if (.not.b2mod_math_initialised) call b2mod_math_init
      if ((art_rad.ne.0.0_R8) .and. (ncall.eq.0)) then
          art_rad=abs(art_rad)
          write(*,*)
     1     'Artificial radiation model active, art_rad = ',art_rad
      endif
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
#ifndef NO_OPENMP_B2SQEL
!$OMP PARALLEL DEFAULT(none)
!$OMP& private(ix, iy, is),
!$OMP& shared(rsa, rra, rqa, rrd, rbr, rza, rz2, rpt, rpi, wrk0)
!$OMP& shared(rlsa, rlra, rlqa, rlrd, rlbr, rlza, rlz2, rlpt, rlpi)
!$OMP& shared(zamin, zamax, phm0, phm1, phm2, ev, ns, nx, ny)
!$OMP& shared(art_rad, ismain, ncall, TNORM, TE)
!$OMP DO
#endif
*   ..compute wrk0(,)=log(te(,)/ev)
      do iy = -1, ny
       do ix = -1, nx
        wrk0(ix,iy) = log(te(ix,iy)/ev)
       enddo
      enddo
#ifndef NO_OPENMP_B2SQEL
!$OMP END DO
#endif

* ..main computation
*   ..compute rate coefficients

#ifndef NO_OPENMP_B2SQEL
!$OMP DO
#endif
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
*      ..compute ionisation rate coefficient
         rsa(ix,iy,is) = phm0*
     &     expu(rlsa(ix,iy,0,is)+rlsa(ix,iy,1,is)*wrk0(ix,iy))
*      ..compute recombination rate coefficient
         rra(ix,iy,is) = phm1*
     &     expu(rlra(ix,iy,0,is)+rlra(ix,iy,1,is)*wrk0(ix,iy))
*      ..compute heat loss rate coefficient
         rqa(ix,iy,is) = phm2*
     &     ev*expu(rlqa(ix,iy,0,is)+rlqa(ix,iy,1,is)*wrk0(ix,iy))
*      ..compute line radiation rate coefficient
         rrd(ix,iy,is) = phm2*
     &     ev*expu(rlrd(ix,iy,0,is)+rlrd(ix,iy,1,is)*wrk0(ix,iy))
*      ..compute bremsstrahlung radiation rate coefficient
         rbr(ix,iy,is) = phm2*
     &     ev*expu(rlbr(ix,iy,0,is)+rlbr(ix,iy,1,is)*wrk0(ix,iy))
*      ..compute effective charge state
         rza(ix,iy,is) = max(zamin(is),min(zamax(is), 
     &     rlza(ix,iy,0,is)+rlza(ix,iy,1,is)*wrk0(ix,iy)))
*      ..compute effective square charge
         rz2(ix,iy,is) = max(zamin(is)**2,min(zamax(is)**2,
     &     rlz2(ix,iy,0,is)+rlz2(ix,iy,1,is)*wrk0(ix,iy)))
*      ..compute cumulative ionisation potential
         rpt(ix,iy,is) = 
     &     rlpt(ix,iy,0,is)+rlpt(ix,iy,1,is)*wrk0(ix,iy)
*      ..compute effective ionisation potential
         rpi(ix,iy,is) = 
     &     rlpi(ix,iy,0,is)+rlpi(ix,iy,1,is)*wrk0(ix,iy)
        enddo
       enddo
      enddo
#ifndef NO_OPENMP_B2SQEL
!$OMP END DO
#endif

* ..artificial radiation                DPC 2000.08.25
      if(art_rad.ne.0.0_R8) then
        if(ncall.eq.0) then
          art_rad=abs(art_rad)
          write(*,*)
     1     'Artificial radiation model active, art_rad = ',art_rad
        endif
        is=ismain    ! XPB 2007.11.02 was ns-1
#ifndef NO_OPENMP_B2SQEL
!$OMP DO
#endif
        do iy = -1, ny
          do ix = -1, nx
*      ..compute artificial radiation rate coefficient
            TNORM = MAX(0.000001_R8,TE(IX,IY)/(15.0_R8*EV))
            rqa(ix,iy,is) = rqa(ix,iy,is) +
     1       art_rad*2.0E-31_R8/(TNORM**1.5_R8+1.0_R8/TNORM**3)
            rrd(ix,iy,is) = rrd(ix,iy,is) +
     1       art_rad*2.0E-31_R8/(TNORM**1.5_R8+1.0_R8/TNORM**3)
          enddo
        enddo
#ifndef NO_OPENMP_B2SQEL
!$OMP END DO
#endif
      endif
#ifndef NO_OPENMP_B2SQEL
!$OMP END PARALLEL
#endif
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sqel

      end subroutine b2sqel
