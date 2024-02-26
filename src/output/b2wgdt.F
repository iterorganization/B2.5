*-----------------------------------------------------------------------
*.specification

      subroutine b2wgdt (nprt, idout, nx, ny, ns, delne, delni, delna, 
     &  delua, delte, delti, delpo)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
      real (kind=R8) ::
     *  delne(-1:nx,-1:ny), delni(-1:nx,-1:ny,0:1), 
     &  delna(-1:nx,-1:ny,0:ns-1), delua(-1:nx,-1:ny,0:ns-1),
     &  delte(-1:nx,-1:ny), delti(-1:nx,-1:ny), delpo(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGDT provides graphical output of the change in the state.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the change in the plasma state
*     computed in the course of a B2 calculation.
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
*     delne, delni, delna, delua, delte, delti, delpo - real*8 array, input.
*     These arrays specify the changes in (respectively) density,
*     parallel velocity, electron temperature, atom temperature, and
*     electric potential.
*     See the calling routine for description.
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
     *   xf(-1:nx),xc(-1:nx),yf(-1:ny),yc(-1:ny)
*   ..procedures
      intrinsic max
      real (kind=R8) ::
     *   scal10, damax
      logical strmas
      external xertst, scal10, damax, strmas
      external b2wgxp

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgdt')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..compute tthi
      tthi = ev*scal10(0.8_R8*
     &  max(damax(n2,delte,1)/ev,damax(n2,delti,1)/ev))
*   ..set coordinates
      do ix = -1, nx
       xf(ix) = ix
       xc(ix) = ix+0.5_R8
      enddo
      do iy = -1, ny
       yf(iy) = iy
       yc(iy) = iy+0.5_R8
      enddo

* ..make plots
*   ..make contour plots of (ne) change
      if (strmas(idout,'cpne')) then
       t0 = scal10(0.8_R8*damax(n2,delne,1))
       write (header,fmth) 'ELECTRON DENSITY CHANGE (/m**3).', 
     &  t0/50.0_R8
       write (nprt,fmtp) '+cpne', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   delne, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1   24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plots of (ni) change
      if (strmas(idout,'cpni')) then
        t0 = scal10(0.8_R8*damax(n2,delni(-1,-1,0),1))
       write (header,fmth) 'ALL ATOM DENSITY CHANGE (/m**3).', 
     &  t0/50.0_R8
       write (nprt,fmtp) '+cpni', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   delni(-1,-1,0), 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0,
     1   24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plots of (te) change
      if (strmas(idout,'cpte')) then
       write (header,fmth) 'ELECTRON TEMP. CHANGE (eV).', 
     &  tthi/(50.0_R8*ev)
       write (nprt,fmtp) '+cpte', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   delte, 0, ev, -1.2_R8*tthi, 1.2_R8*tthi, 
     &   24, -tthi/10.0_R8, tthi/10.0_R8, 10)
      endif
*   ..make contour plots of (ti) change
      if (strmas(idout,'cpti')) then
       write (header,fmth) 'ATOM TEMPERATURE CHANGE (eV).', 
     &  tthi/(50.0_R8*ev)
       write (nprt,fmtp) '+cpti', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   delti, 0, ev, -1.2_R8*tthi, 1.2_R8*tthi, 
     &   24, -tthi/10.0_R8, tthi/10.0_R8, 10)
      endif
*   ..make contour plots of (po) change
      if (strmas(idout,'cppo')) then
       t0 = scal10(0.8_R8*damax(n2,delpo,1))
       write (header,fmth) 'ELECTRIC POTENTIAL CHANGE (V).', t0/50.0_R8
       write (nprt,fmtp) '+cppo', header
       call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &   delpo, 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 
     &   24, -t0/10.0_R8, t0/10.0_R8, 10)
      endif
*   ..make contour plots of (na) change
      if (strmas(idout,'cpna')) then
       do is = 0, ns-1
        t0 = scal10(0.8_R8*damax(n2,delna(-1,-1,is),1))
        write (header,fmth2) 'DENSITY CHANGE (/m**3).', is, t0/50.0_R8
        write (nprt,fmtp) '+cpna', header
        call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &    delna(-1,-1,is), 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 24,
     &    -t0/10.0_R8, t0/10.0_R8, 10)
       enddo
      endif
*   ..make contour plots of (ua) change
       if (strmas(idout,'cpua')) then
        do is = 0, ns-1
         t0 = scal10(0.8_R8*damax(n2,delua(-1,-1,is),1))
         write (header,fmth2) 'VELOCITY CHANGE (m/s).', is, t0/50.0_R8
         write (nprt,fmtp) '+cpua', header
         call b2wgxp (nx, ny, header, xc, xf, yc, yf,
     &     delua(-1,-1,is), 0, 1.0_R8, -1.2_R8*t0, 1.2_R8*t0, 24,
     &     -t0/10.0_R8, t0/10.0_R8, 10)
        enddo
       endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgdt

      end subroutine b2wgdt
