*-----------------------------------------------------------------------
*.specification

      subroutine b2wggm (nprt, idout, nx, ny, title, crx, cry,
     &  fpsi, ffbz, bb)
      use b2mod_types
      use b2mod_indirect
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny
      character idout*(*), title*(*)
      real (kind=R8) ::
     *   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     &  fpsi(-1:nx,-1:ny,0:3), ffbz(-1:nx,-1:ny,0:3),
     &  bb(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGGM produces graphical displays of the B2 geometry and
*     magnetic field.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the B2 geometry and of the
*     magnetic field. Arrays crx and cry provide the coordinates
*     of the four corners of all the cells in the mesh. Arrays fpsi
*     and ffbz provide the flux functions at the corner locations.
*     Array bb provides the magnetic field at cell centers.
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
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     title - character*(*), input.
*     title specifies a label to identify the geometry.
*
*     crx, cry, fpsi, ffbz, bb - real*8 array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2
      real (kind=R8) ::
     *   t0, t1
      real (kind=R8) ::
     *   wrkx(-1:nx,-1:ny), wrky(-1:nx,-1:ny)
      character fmth*(*), fmtp*(*), header*60, headr2*120
      parameter (
     &  fmth='(a,t44,''inc=(5*)'',1pe9.2)',
     &  fmtp='(2x,a6,4x,a)')
*   ..procedures
      real (kind=R8) ::
     *   scal10, smin, smax, damax
      logical strmas
      external xertst, scal10, smax, damax, strmas
      external b2wgx0, b2wgx1, b2wgx2, b2wgx3

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wggm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..make plots
*   ..plot the geometry and mesh
      headr2 = title
      if (strmas(idout,'mpcr')) then
       write (nprt,fmtp) '+mpcr', headr2
*      (cartesian representation)
       call b2wgx0 (nx, ny, headr2, crx, cry)
*      (physical representation)
       call b2wgx1 (nx, ny, headr2, crx, cry)
      endif
*   ..contour plots of fpsi
*     (plots use only one of the four sets of vertex data)
      t0 = smin(n2*4,fpsi,1)
      t1 = smax(n2*4,fpsi,1)
      write (header,fmth) 'FLUX FUNCTION FPSI (Wb).', (t1-t0)/20.
      if (strmas(idout,'mpff')) then
       write (nprt,fmtp) '+mpff', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, fpsi(-1,-1,0),
     1   3, 1.0_R8,
     &   t0, t1, 4, t0, t1, 20)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, fpsi(-1,-1,0),
     1   3, 1.0_R8,
     &   t0, t1, 4, t0, t1, 20, wrkx, wrky)
      endif
*   ..contour plots of ffbz
*     (plots use only one of the four sets of vertex data)
      t0 = smin(n2*4,ffbz,1)
      t1 = smax(n2*4,ffbz,1)
      write (header,fmth) 'FLUX FUNCTION FFBZ (Wb/m).', (t1-t0)/20.
      if (strmas(idout,'mpff')) then
       write (nprt,fmtp) '+mpff', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, ffbz(-1,-1,0),
     1   3, 1.0_R8,
     &   t0, t1, 4, t0, t1, 20)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, ffbz(-1,-1,0),
     1   3, 1.0_R8,
     &   t0, t1, 4, t0, t1, 20, wrkx, wrky)
      endif
*   ..contour plots of bx
      t0 = scal10(0.8_R8*damax(n2,bb(-1,-1,0),1))
      write (header,fmth) 'MAGNETIC FIELD BX (T).', t0/100.
      if (strmas(idout,'mpbb')) then
       write (nprt,fmtp) '+mpbb', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, bb(-1,-1,0),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, bb(-1,-1,0),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10,
     &   wrkx, wrky)
      endif
*   ..contour plots of by
      t0 = scal10(0.8_R8*damax(n2,bb(-1,-1,1),1))
      write (header,fmth) 'MAGNETIC FIELD BY (T).', t0/100.
      if (strmas(idout,'mpbb')) then
       write (nprt,fmtp) '+mpbb', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, bb(-1,-1,1),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, bb(-1,-1,1),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10,
     &   wrkx, wrky)
      endif
*   ..contour plots of bz
      t0 = scal10(0.8_R8*damax(n2,bb(-1,-1,2),1))
      write (header,fmth) 'MAGNETIC FIELD BZ (T).', t0/100.
      if (strmas(idout,'mpbb')) then
       write (nprt,fmtp) '+mpbb', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, bb(-1,-1,2),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, bb(-1,-1,2),
     1   0, 1.0_R8,
     &   -1.2_R8*t0, 1.2_R8*t0, 48, -t0/20.0_R8, t0/20.0_R8, 10,
     &   wrkx, wrky)
      endif
*   ..contour plots of abs(bb)
      t0 = scal10(0.8_R8*smax(n2,bb(-1,-1,3),1))
      write (header,fmth) 'MAGNETIC FIELD BTOT (T).', t0/100.
      if (strmas(idout,'mpbb')) then
       write (nprt,fmtp) '+mpbb', header
*      (cartesian representation)
       call b2wgx2 (nx, ny, header, bb(-1,-1,3),
     1   0, 1.0_R8,
     &   0.0_R8, 1.2_R8*t0, 24, 0.0_R8, t0/20.0_R8, 5)
*      (physical representation)
       call b2wgx3 (nx, ny, header, crx, cry, bb(-1,-1,3),
     1   0, 1.0_R8,
     &   0.0_R8, 1.2_R8*t0, 24, 0.0_R8, t0/20.0_R8, 5, wrkx, wrky)
      endif

* ..return
      flush (nprt)
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wggm

      end subroutine b2wggm
