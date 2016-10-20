*-----------------------------------------------------------------------
*.specification

      subroutine b2wgpj (nprt, idout, nx, ny, ns, na, ua, te, ti, po,
     &  fna, fhe, fhi, fch, ne, ni, pr, ro, fne, fni, wrk0)
      use b2mod_types
      use b2mod_rates
      use b2mod_constants
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), pr(-1:nx,-1:ny),
     &  ro(-1:nx,-1:ny), fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGPJ provides graphical output of the plasma state on a uniform
*     cartesian mesh.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the plasma state that results from
*     a b2 calculation. Arrays na, ua, te, ti, and po provide the state
*     of the plasma. Arrays fna, fhe, fhi, and fch provide the particle
*     flux, electron heat flux, atom heat flux, and the electric
*     current. This routine does not require any metric information.
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
*     na, ua, te, ti, po, fna, fhe, fhi, fch - real*8 array, input.
*     ne, ni, pr, ro, fne, fni - real*8 array, input.
*     See the calling routine for description.
*
*     wrk0 - (-1:nx,-1:ny) real*8 array, workspace.
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
      integer nmax, n2, ix, iy, is
      parameter (nmax=300)
      real (kind=R8) ::
     *   sc, t0
      real (kind=R8) ::
     *   xc(-1:nmax), xf(-1:nmax), yc(-1:nmax), yf(-1:nmax),
     &  wxf1(-1:nmax), wxf2(-1:nmax), wyf1(-1:nmax), wyf2(-1:nmax)
      character fmth*(*), fmth2*(*), fmtl*(*), fmtl2*(*), fmtp*(*),
     &  header*60
      parameter (
     &  fmth='(a,t44,''inc=(5*)'',1pe9.2)',
     &  fmth2='(a,t35,''is ='',i3,t44,''inc=(5*)'',1pe9.2)',
     &  fmtl='(a,t44,''logarithmic scale'')',
     &  fmtl2='(a,t35,''is ='',i3,t44,''logarithmic scale'')',
     &  fmtp='(2x,a6,4x,a)')
*   ..procedures
      intrinsic max
      real (kind=R8) ::
     *   scal10, samax, b2ssum
      logical strmas
      external subini, subend, xertst, scal10, samax, b2ssum,
     &  strmas, pagini, pagend, gsln, gpl_d
      external b2xppe, b2xpmc, b2wgxp, b2wgxq

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgpj')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test if sufficient workspace
      call xertst (nx.le.nmax.and.ny.le.nmax,
     &  'faulty parameter nmax')

*   ..set coordinates
      do ix = -1, nx
       xf(ix) = ix
       xc(ix) = ix+0.5_R8
      enddo
      do iy = -1, ny
       yf(iy) = iy
       yc(iy) = iy+0.5_R8
      enddo

* ..make two-dimensional plots
*   ..make log-contour plots of electron density
      if (strmas(idout,'cpne')) then
       write (header,fmtl) 'ELECTRON DENSITY (/m**3).'
       write (nprt,fmtp) '+cpne', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, ne, 0, 1.0_R8)
      endif
*   ..make log-contour plots of total atom density
      if (strmas(idout,'cpni')) then
       write (header,fmtl) 'TOTAL ATOM DENSITY (/m**3).'
       write (nprt,fmtp) '+cpni', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, ni, 0, 1.0_R8)
      endif
*   ..make log-contour plots of electron temperature
      if (strmas(idout,'cpte')) then
       write (header,fmtl) 'ELECTRON TEMPERATURE (eV).'
       write (nprt,fmtp) '+cpte', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, te, 0, ev)
      endif
*   ..make log-contour plots of all atom temperature
      if (strmas(idout,'cpti')) then
       write (header,fmtl) 'ALL ATOM TEMPERATURE (eV).'
       write (nprt,fmtp) '+cpti', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, ti, 0, ev)
      endif
*   ..make log-contour plots of total pressure
      if (strmas(idout,'cppr')) then
       write (header,fmtl) 'TOTAL PRESSURE (N/m**2).'
       write (nprt,fmtp) '+cppr', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, pr, 0, 1.0_R8)
      endif
*   ..make contour plots of electric potential
      if (strmas(idout,'cppo')) then
       t0 = scal10(0.8_R8*samax(n2,po,1))
       write (header,fmth) 'ELECTRIC POTENTIAL (V).', t0/100.
       write (nprt,fmtp) '+cppo', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf, po, 0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10)
      endif
*   ..make log-contour plots of fractional atom density
      if (strmas(idout,'cpfn')) then
       do is = 0, ns-1
        do iy = -1, ny
         do ix = -1, nx
          wrk0(ix,iy) = na(ix,iy,is)/(ni(ix,iy,0)+ni_ext(ix,iy,0))
         enddo
        enddo
        write (header,fmtl2) 'FRACTIONAL ATOM DENSITY (unit).', is
        write (nprt,fmtp) '+cpfn', header
        call b2wgxq (nx, ny, header, xc, xf, yc, yf, wrk0, 0, 1.0_R8)
       enddo
      endif
*   ..make log-contour plots of electron pressure
      if (strmas(idout,'cppe')) then
       call b2xppe (nx, ny, ne, te, wrk0)
       write (header,fmtl) 'ELECTRON PRESSURE (N/m**2).'
       write (nprt,fmtp) '+cppe', header
       call b2wgxq (nx, ny, header, xc, xf, yc, yf, wrk0, 0, 1.0_R8)
      endif
*   ..make contour plots of parallel velocity
      if (strmas(idout,'cpua')) then
       do is = 0, ns-1
        t0 = scal10(0.8_R8*samax(n2,ua(-1,-1,is),1))
        write (header,fmth2) 'PARALLEL VELOCITY (m/s).', is, t0/100.
        write (nprt,fmtp) '+cpua', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf, ua(-1,-1,is),
     &    0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     &    -t0/20.0_R8, t0/20.0_R8, 10)
       enddo
      endif
*   ..make contour plots of mach number
      if (strmas(idout,'cpmc')) then
       do is = 0, ns-1
        call b2xpmc (nx, ny, ua(-1,-1,is), pr, ro, wrk0)
        t0 = 1
        write (header,fmth2) 'MACH NUMBER.', is, t0/100.
        write (nprt,fmtp) '+cpmc', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf, wrk0,
     &    0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     &    -t0/20.0_R8, t0/20.0_R8, 10)
       enddo
      endif
*   ..make contour plots of electron heat flux
      if (strmas(idout,'cphe')) then
       t0 = scal10(0.8_R8*samax(n2,fhe(-1,-1,0),1))
       write (header,fmth)
     1	'x-directed electron heat flux (unit).', t0/100.
       write (nprt,fmtp) '+cphe', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf, fhe(-1,-1,0),
     1	 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     2   -t0/20.0_R8, t0/20.0_R8, 10)
       t0 = scal10(0.8_R8*samax(n2,fhe(-1,-1,1),1))
       write (header,fmth)
     1	'y-directed electron heat flux (unit).', t0/100.
       write (nprt,fmtp) '+cphe', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf, fhe(-1,-1,1),
     1	 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     2   -t0/20.0_R8, t0/20.0_R8, 10)
      endif
*   ..make contour plots of all atom heat flux
      if (strmas(idout,'cphi')) then
       t0 = scal10(0.8_R8*samax(n2,fhi(-1,-1,0),1))
       write (header,fmth)
     1	'x-directed all atom heat flux (unit).', t0/100.
       write (nprt,fmtp) '+cphi', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf, fhi(-1,-1,0),
     1	 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     2   -t0/20.0_R8, t0/20.0_R8, 10)
       t0 = scal10(0.8_R8*samax(n2,fhi(-1,-1,1),1))
       write (header,fmth)
     1	'y-directed all atom heat flux (unit).', t0/100.
       write (nprt,fmtp) '+cphi', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf, fhi(-1,-1,1),
     1	 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 48,
     2   -t0/20.0_R8, t0/20.0_R8, 10)
      endif

* ..make one-dimensional plots
*   ..compute electron and atom flux as function of x
      do ix = 0, nx
       wxf1(ix) = b2ssum(ny+2,fne(ix,-1,0),nx+2)
       wxf2(ix) = b2ssum(ny+2,fni(ix,-1,0),nx+2)
      enddo
*   ..compute common scale for particle fluxes
      t0 = max(samax(nx+1,wxf1(0),1),samax(nx+1,wxf2(0),1))
      if (t0.eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(t0)
      endif
*   ..make plots of particle fluxes as function of x
      if (strmas(idout,'fnex')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRON (ATOM) FLUX', 'AS FUNCTION OF X'
       write (nprt,fmtp) '+fnex', header
       call pagini (0, xc(-1), xc(nx), -sc, sc, 'x', 'fnex',
     &   header)
       call gsln (1)
       call gpl_d (nx+1, xf(0), wxf1(0))
       call gsln (2)
       call gpl_d (nx+1, xf(0), wxf2(0))
       call pagend ()
      endif
*   ..compute heat fluxes as function of x
      do ix = 0, nx
       wxf1(ix) = b2ssum(ny+2,fhe(ix,-1,0),nx+2)
       wxf2(ix) = b2ssum(ny+2,fhi(ix,-1,0),nx+2)
      enddo
*   ..compute common scale for two heat fluxes
      t0 = max(samax(nx+1,wxf1(0),1),samax(nx+1,wxf2(0),1))
      if (t0.eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(t0)
      endif
*   ..make plots of heat fluxes as function of x
      if (strmas(idout,'fhex')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRON (ATOM) HEAT FLUX', 'AS FUNCTION OF X'
       write (nprt,fmtp) '+fhex', header
       call pagini (0, xc(-1), xc(nx), -sc, sc, 'x', 'fhex',
     &   header)
       call gsln (1)
       call gpl_d (nx+1, xf(0), wxf1(0))
       call gsln (2)
       call gpl_d (nx+1, xf(0), wxf2(0))
       call pagend ()
      endif
*   ..compute y-directed heat fluxes as function of x on southern boundary
      do ix = 0, nx
       wxf1(ix) = fhe(ix,0,1)
       wxf2(ix) = fhi(ix,0,1)
      enddo
*   ..compute common scale for two heat fluxes
      t0 = max(samax(nx+1,wxf1(0),1),samax(nx+1,wxf2(0),1))
      if (t0.eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(t0)
      endif
*   ..make plots of heat fluxes as function of x
      if (strmas(idout,'fhex')) then
       write (header,'(a,t44,a)')
     &   'Y-DIRECTED ELECTRON (ATOM) HEAT FLUX (SO)',
     1	 'AS FUNCTION OF X'
       write (nprt,fmtp) '+fhex', header
       call pagini (0, xc(-1), xc(nx), -sc, sc, 'x', 'fhex',
     &   header)
       call gsln (1)
       call gpl_d (nx+1, xf(0), wxf1(0))
       call gsln (2)
       call gpl_d (nx+1, xf(0), wxf2(0))
       call pagend ()
      endif
*   ..compute electric current as function of x
      do ix = 0, nx
       wxf1(ix) = b2ssum(ny+2,fch(ix,-1,0),nx+2)
      enddo
*   ..scale for electric current
      if (samax(nx+1,wxf1(0),1).eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(samax(nx+1,wxf1(0),1))
      endif
*   ..make plots of electric current as function of x
      if (strmas(idout,'fchx')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRIC CURRENT', 'AS FUNCTION OF X'
       write (nprt,fmtp) '+fchx', header
       call pagini (0, xc(-1), xc(nx), -sc, sc, 'x', 'fchx',
     &   header)
       call gpl_d (nx+1, xf(0), wxf1(0))
       call pagend ()
      endif
*   ..make plots of atom flux as function of x
      do is = 0, ns-1
       if (strmas(idout,'fnax')) then
        do ix = 0, nx
         wxf1(ix) = b2ssum(ny+2,fna(ix,-1,0,is),nx+2)
        enddo
        if (samax(nx+1,wxf1(0),1).eq.0.0_R8) then
         sc = 1.0_R8
        else
         sc = scal10(samax(nx+1,wxf1(0),1))
        endif
        write (header,'(a,'' ('',i3,'')'',t44,a)')
     &    'INTEGRATED PARTICLE FLUX', is, 'AS FUNCTION OF X'
        write (nprt,fmtp) '+fnax', header
        call pagini (0, xc(-1), xc(nx), -sc, sc, 'x', 'fnax',
     &    header)
        call gpl_d (nx+1, xf(0), wxf1(0))
        call pagend ()
       endif
      enddo
*   ..compute electron and atom flux as function of y
      do iy = 0, ny
       wyf1(iy) = b2ssum(nx+2,fne(-1,iy,1),1)
       wyf2(iy) = b2ssum(nx+2,fni(-1,iy,1),1)
      enddo
*   ..cmpute commoon scale for particle fluxes
      t0 = max(samax(ny+1,wyf1(0),1),samax(ny+1,wyf2(0),1))
      if (t0.eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(t0)
      endif
*   ..make plots of particle fluxes as function of y
      if (strmas(idout,'fney')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRON (ATOM) FLUX', 'AS FUNCTION OF Y'
       write (nprt,fmtp) '+fney', header
       call pagini (0, yc(-1), yc(ny), -sc, sc, 'y', 'fney',
     &   header)
       call gsln (1)
       call gpl_d (ny+1, yf(0), wyf1(0))
       call gsln (2)
       call gpl_d (ny+1, yf(0), wyf2(0))
       call pagend ()
      endif
*   ..compute heat fluxes as function of y
      do iy = 0, ny
       wyf1(iy) = b2ssum(nx+2,fhe(-1,iy,1),1)
       wyf2(iy) = b2ssum(nx+2,fhi(-1,iy,1),1)
      enddo
*   ..compute common scale for two heat fluxes
      t0 = max(samax(ny+1,wyf1(0),1),samax(ny+1,wyf2(0),1))
      if (t0.eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(t0)
      endif
*   ..make plots of heat fluxes as function of y
      if (strmas(idout,'fhey')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRON (ATOM) HEAT FLUX', 'AS FUNCTION OF Y'
       write (nprt,fmtp) '+fhey', header
       call pagini (0, yc(-1), yc(ny), -sc, sc, 'y', 'fhey',
     &   header)
       call gsln (1)
       call gpl_d (ny+1, yf(0), wyf1(0))
       call gsln (2)
       call gpl_d (ny+1, yf(0), wyf2(0))
       call pagend ()
      endif
*   ..compute electric current as function of y
      do iy = 0, ny
       wyf1(iy) = b2ssum(nx+2,fch(-1,iy,1),1)
      enddo
*   ..scale for electric current
      if (samax(ny+1,wyf1(0),1).eq.0.0_R8) then
       sc = 1.0_R8
      else
       sc = scal10(samax(ny+1,wyf1(0),1))
      endif
*   ..make plots of electric current as function of y
      if (strmas(idout,'fchy')) then
       write (header,'(a,t44,a)')
     &   'INTEGRATED ELECTRIC CURRENT', 'AS FUNCTION OF Y'
       write (nprt,fmtp) '+fchy', header
       call pagini (0, yc(-1), yc(ny), -sc, sc, 'y', 'fchy',
     &   header)
       call gpl_d (ny+1, yf(0), wyf1(0))
       call pagend ()
      endif
*   ..make plots of atom flux as function of y
      do is = 0, ns-1
       if (strmas(idout,'fnay')) then
        do iy = 0, ny
         wyf1(iy) = b2ssum(nx+2,fna(-1,iy,1,is),1)
        enddo
        if (samax(ny+1,wyf1(0),1).eq.0.0_R8) then
         sc = 1.0_R8
        else
         sc = scal10(samax(ny+1,wyf1(0),1))
        endif
        write (header,'(a,'' ('',i3,'')'',t44,a)')
     &    'INTEGRATED PARTICLE FLUX', is, 'AS FUNCTION OF Y'
        write (nprt,fmtp) '+fnay', header
        call pagini (0, yc(-1), yc(ny), -sc, sc, 'y', 'fnay',
     &    header)
        call gpl_d (ny+1, yf(0), wyf1(0))
        call pagend ()
       endif
      enddo

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgpj

      end subroutine b2wgpj
