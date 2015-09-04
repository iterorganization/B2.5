*-----------------------------------------------------------------------
*.specification

      subroutine b2wgcr (nprt, idout, nx, ny, ns, resco, resmo, resmt,
     &  reshe, reshi, respo, corpa, corua, corut, corte, corti, corpo)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
      real (kind=R8) ::
     *   resco(-1:nx,-1:ny,0:ns-1), resmo(-1:nx,-1:ny,0:ns-1),
     &  resmt(-1:nx,-1:ny), reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny),
     &  respo(-1:nx,-1:ny)
      real (kind=R8) ::
     *   corpa(-1:nx,-1:ny,0:ns-1), corua(-1:nx,-1:ny,0:ns-1),
     &  corut(-1:nx,-1:ny), corte(-1:nx,-1:ny), corti(-1:nx,-1:ny),
     &  corpo(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGCR provides graphical output of the residuals and the
*     corrections.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the residuals and the corrections
*     computed in the course of a b2 calculation.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     idout - character*(*), input.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     resco, resmo, resmt, reshe, reshi, respo - real*8 array, input.
*     These arrays specify the most recent residuals of the equations
*     for (respectively) continuity (each species), parallel momentum
*     (each species), parallel momentum (all species), electron heat,
*     atom heat, and electric potential.
*     See the calling routine for description.
*
*     corpa, corua, corut, corte, corti, corpo - real*8 array, input.
*     These arrays specify the most recent computed corrections for
*     (respectively) pressure (each species), parallel velocity (each
*     species), parallel velocity (all species), electron temperature,
*     atom temperature, and electric potential.
*     See the calling routine for description.
*
*     wrk0 - (-1:nx,-1:ny) workspace
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, ix, iy, is
      real (kind=R8) ::
     *   tthi, t0
      character fmth*(*), fmth2*(*), fmtp*(*), header*60
      parameter (
     &  fmth='(a,t44,''inc=(5*)'',1pe9.2)',
     &  fmth2='(a,t35,''is ='',i3,t44,''inc=(5*)'',1pe9.2)',
     &  fmtp='(2x,a6,4x,a)')
      real (kind=R8) ::
     *   xf(-1:nx),xc(-1:nx),yf(-1:ny),yc(-1:ny),zero(-1:nx,-1:ny)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)
*   ..procedures
      intrinsic max
      real (kind=R8) ::
     *   samax, scal10
      logical strmas
      external subini, subend, xertst, samax, b2scopy, sfill,
     &  scal10, strmas
      external b2xpne, b2xpni, b2xxgs, b2wgxp

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgcr')
*   ..set n2 and zero
      n2 = (nx+2)*(ny+2)
      zero = 0.0_R8
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..compute tthi
      tthi = ev*scal10(0.8_R8*
     &  max(samax(n2,corte,1)/ev,samax(n2,corti,1)/ev))
*   ..set coordinates
      do ix = -1, nx
       xf(ix) = ix
       xc(ix) = ix+0.5_R8
      enddo
      do iy = -1, ny
       yf(iy) = iy
       yc(iy) = iy+0.5_R8
      enddo

* ..make plots of residuals
*   ..make contour plots of (ne) residual
      if (strmas(idout,'cpne')) then
*      (ignore boundaries and singular regions)
       call b2xpne (nx, ny, ns, rza, resco, zero, wrk0)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or. 
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'ELECTRON CONTINUITY RESID (1/s).', 
     &  t0/50.0_R8
       write (nprt,fmtp) '+cpne', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plots of (ni) residual
      if (strmas(idout,'cpni')) then
*      (ignore boundaries and singular regions)
       call b2xpni (nx, ny, ns, resco, wrk0)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or. 
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'ALL ATOM CONTINUITY RESID (1/s).', 
     &  t0/50.0_R8
       write (nprt,fmtp) '+cpni', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plot of (te) residual
      if (strmas(idout,'cpte')) then
*      (ignore boundaries and singular regions)
       call b2scopy (n2, reshe, 1, wrk0, 1)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'ELECTRON HEAT RESIDUAL (W).', t0/50.0_R8
       write (nprt,fmtp) '+cpte', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plot of (ti) residual
      if (strmas(idout,'cpti')) then
*      (ignore boundaries and singular regions)
       call b2scopy (n2, reshi, 1, wrk0, 1)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'ATOM HEAT RESIDUAL (W).', t0/50.0_R8
       write (nprt,fmtp) '+cpti', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plot of (po) residual
      if (strmas(idout,'cppo')) then
*      (ignore boundaries and singular regions)
       call b2scopy (n2, respo, 1, wrk0, 1)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'EL. POTENTIAL RESIDUAL (A).', t0/50.0_R8
       write (nprt,fmtp) '+cppo', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plots of (na) residual
      if (strmas(idout,'cpna')) then
       do is = 0, ns-1
*       (ignore boundaries and singular regions)
        call b2scopy (n2, resco(-1,-1,is), 1, wrk0, 1)
        call b2xxgs (nx, ny, 0.0_R8, wrk0, 0)
        call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
        call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
        call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
        call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,1).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
        t0 = scal10(0.8_R8*samax(n2,wrk0,1))
        write (header,fmth2) 'CONTINUITY RESIDUAL (1/s).', is, 
     &   t0/50.0_R8
        write (nprt,fmtp) '+cpna', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &    wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	  24, -t0/10.0_R8, t0/10.0_R8, 10)
       enddo
      endif
*   ..make contour plots of (ua) residual
      if (strmas(idout,'cpua')) then
       do is = 0, ns-1
*       (ignore boundaries and singular regions)
        call b2scopy (n2, resmo(-1,-1,is), 1, wrk0, 1)
        call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
        call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
        call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
        call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
        call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
        do ix = 0, nx-1
          do iy = 0, ny-1
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &         bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &         region(ix,iy,0).eq.0 .or. resignore(ix,iy,2).eq.0)
     >       wrk0(ix,iy)=0.0_R8
          enddo
        enddo
        t0 = scal10(0.8_R8*samax(n2,wrk0,1))
        write (header,fmth2) 'PARALLEL MOMENTUM RESIDUAL (J/m).',
     &    is, t0/50.0_R8
        write (nprt,fmtp) '+cpua', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &    wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	  24, -t0/10.0_R8, t0/10.0_R8, 10)
       enddo
      endif
*   ..make contour plot of (ua) residual, all species
      if (strmas(idout,'cpua')) then
*      (ignore boundaries and singular regions)
       call b2scopy (n2, resmt, 1, wrk0, 1)
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,-1), 1)
       call sfill (nx+2, 0.0_R8, wrk0(-1,ny), 1)
       call sfill (ny, 0.0_R8, wrk0(-1,0), nx+2)
       call sfill (ny, 0.0_R8, wrk0(nx,0), nx+2)
       do ix = 0, nx-1
         do iy = 0, ny-1
           if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &        bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 .or.
     &        region(ix,iy,0).eq.0 .or. resignore(ix,iy,2).eq.0)
     >       wrk0(ix,iy)=0.0_R8
         enddo
       enddo
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth)
     &  'PARALLEL MOMENTUM RESIDUAL (J/m), all species', t0/50.0_R8
       write (nprt,fmtp) '+cpua', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif

* ..make plots of corrections
*   ..make contour plots of (pr) correction
      if (strmas(idout,'cppr')) then
       call b2xpni (nx, ny, ns, corpa, wrk0)
       t0 = scal10(0.8_R8*samax(n2,wrk0,1))
       write (header,fmth) 'TOTAL PRESSURE CORRECTION (J/m**3).', 
     &  t0/50.0_R8
       write (nprt,fmtp) '+cppr', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   wrk0, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plot of (te) correction
      if (strmas(idout,'cpte')) then
       write (header,fmth) 'EL. TEMP. CORRECTION (eV).', 
     &  tthi/(50.0_R8*ev)
       write (nprt,fmtp) '+cpte', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   corte, 0, ev, -1.2_R8*tthi, 1.2_R8*tthi, 
     &   24, -tthi/10.0_R8, tthi/10.0_R8, 10)
      endif
*   ..make contour plot of (ti) correction
      if (strmas(idout,'cpti')) then
       write (header,fmth) 'ATOM TEMP. CORRECTION (eV).', 
     &  tthi/(50.0_R8*ev)
       write (nprt,fmtp) '+cpti', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   corti, 0, ev, -1.2_R8*tthi, 1.2_R8*tthi,  
     &   24, -tthi/10.0_R8, tthi/10.0_R8, 10)
      endif
*   ..make contour plot of (po) correction
      if (strmas(idout,'cppo')) then
       write (header,fmth) 'EL. POTENTIAL CORRECTION (V).', 
     &  tthi/(50.0_R8*ev)
       write (nprt,fmtp) '+cppo', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   corpo, 0, 1.0_R8, -1.2_R8*tthi/ev, 1.2_R8*tthi/ev,  
     &   24, -tthi/(10.0_R8*ev), tthi/(10.0_R8*ev), 10)
      endif
*   ..make contour plots of (pa) correction
      if (strmas(idout,'cpna')) then
       do is = 0, ns-1
        t0 = scal10(0.8_R8*samax(n2,corpa(-1,-1,is),1))
        write (header,fmth2) 'PRESSURE CORRECTION (J/m**3).', is, 
     &   t0/50.0_R8
        write (nprt,fmtp) '+cpna', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &    corpa(-1,-1,is), 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 24,
     &    -t0/10.0_R8, t0/10.0_R8, 10)
       enddo
      endif
*   ..make contour plots of (ua) correction
      if (strmas(idout,'cpua')) then
       do is = 0, ns-1
        t0 = scal10(0.8_R8*samax(n2,corua(-1,-1,is),1))
        write (header,fmth2) 'VELOCITY CORRECTION (m/s).', is, 
     &   t0/50.0_R8
        write (nprt,fmtp) '+cpua', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &    corua(-1,-1,is), 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 24,
     &    -t0/10.0_R8, t0/10.0_R8, 10)
       enddo
      endif
*   ..make contour plots of (ua) correction, all species
      if (strmas(idout,'cpua')) then
       t0 = scal10(0.8_R8*samax(n2,corut,1))
       write (header,fmth) 'VELOCITY CORRECTION (m/s), all species',
     &   t0/50.0_R8
       write (nprt,fmtp) '+cpua', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   corut, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1	 24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgcr

      end subroutine b2wgcr
