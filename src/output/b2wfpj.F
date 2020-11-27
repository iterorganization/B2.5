*-----------------------------------------------------------------------
*.specification

      subroutine b2wfpj (nprt, idout, nx, ny, ns, pr, ro, rz, wrk0)
      use b2mod_types
      use b2mod_rates
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
      real (kind=R8) ::
     &  pr(-1:nx,-1:ny), ro(-1:nx,-1:ny), rz(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFPJ writes a summary of the plasma state to a formatted file.
*     This routine does not make use of any geometric information.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit nprt summary information describing
*     the plasma state. Unit nprt must be connected to a formatted file.
*     This routine does not make use of any geometric information.
*     The extent of the output is specified by the parameter idout.
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
*     pr, ro, rz - real*8 array, input.
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
      integer n2, ix, iy, is, ix0, iy0, ix1, iy1
*   ..procedures
      logical strmas
      external subini, subend, xertst, strmas
      external b2xvsg, b2xvff, b2xvps, b2xxmm, b2xppb

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfpj')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nprt, nx, ny
      call xertst (1.le.nprt, 'faulty argument nprt')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test plasma state
      call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*   ..test sign of pr, ro, rz
      call b2xvsg (n2, pr, 1, 'pr', '.gt.')
      call b2xvsg (n2, ro, 1, 'ro', '.gt.')
      call b2xvsg (n2, rz, 1, 'rz', '.gt.')

* ..write the state
*   ..write pzmm segment?
      if (strmas(idout,'pzmm')) then
*    ..write heading
       write (nprt,'(/2x,a6)') '+pzmm'
       write (nprt,'(/2x,a/)') 'Minimum and maximum values and'//
     &   ' their locations'
*    ..write min and max of na
       do is = 0, ns-1
        call b2xxmm (nx, ny, na(-1,-1,is), 0, ix0, iy0, ix1, iy1)
        write (nprt,103) 'na', is, na(ix0,iy0,is), ix0, iy0,
     &    na(ix1,iy1,is), ix1, iy1
       enddo
*    ..write min and max of ne
       call b2xxmm (nx, ny, ne, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'ne', ne(ix0,iy0), ix0, iy0,
     &   ne(ix1,iy1), ix1, iy1
*    ..write min and max of ni
       call b2xxmm (nx, ny, ni, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'ni', ni(ix0,iy0,0), ix0, iy0,
     &   ni(ix1,iy1,0), ix1, iy1
*    ..write min and max of te
       call b2xxmm (nx, ny, te, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'te', te(ix0,iy0)/ev, ix0, iy0,
     &   te(ix1,iy1)/ev, ix1, iy1
*    ..write min and max of ti
       call b2xxmm (nx, ny, ti, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'ti', ti(ix0,iy0)/ev, ix0, iy0,
     &   ti(ix1,iy1)/ev, ix1, iy1
*     ..write min and max of po
       call b2xxmm (nx, ny, po, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'po', po(ix0,iy0), ix0, iy0,
     &   po(ix1,iy1), ix1, iy1
*    ..write min and max of pa
       do is = 0, ns-1
        call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, wrk0)
        call b2xxmm (nx, ny, wrk0, 0, ix0, iy0, ix1, iy1)
        write (nprt,103) 'pa', is, wrk0(ix0,iy0), ix0, iy0,
     &    wrk0(ix1,iy1), ix1, iy1
       enddo
*    ..write min and max of pr
       call b2xxmm (nx, ny, pr, 0, ix0, iy0, ix1, iy1)
       write (nprt,104) 'pr', pr(ix0,iy0), ix0, iy0,
     &   pr(ix1,iy1), ix1, iy1
      endif
*    ..write min and max of ua
       do is = 0, ns-1
        call b2xxmm (nx, ny, ua(-1,-1,is), 0, ix0, iy0, ix1, iy1)
        write (nprt,103) 'ua', is, ua(ix0,iy0,is), ix0, iy0,
     &    ua(ix1,iy1,is), ix1, iy1
       enddo
*   ..write pfnx segment?
      if (strmas(idout,'pfnx')) then
       write (nprt,'(/2x,a6)') '+pfnx'
*    ..write reduced fnax
       write (nprt,105) (is,is=0,ns-1)
       do ix = 0, nx
        write (nprt,'(2x,i4,1p,(t10,10e11.3))')
     &    ix, (rfnax(nx,ny,ix,is),is=0,ns-1)
       enddo
      endif
*   ..write pfzx segment?
      if (strmas(idout,'pfzx')) then
       write (nprt,'(/2x,a6)') '+pfzx'
*    ..write reduced fnex, fnix, fhex, fhix, fchx
       write (nprt,106)
       do ix = 0, nx
        write (nprt,'(2x,i4,1p,(t10,5e11.3))')
     &    ix, rfnex(nx,ny,ix), rfnix(nx,ny,ix), rfhex(nx,ny,ix),
     &        rfhix(nx,ny,ix), rfchx(nx,ny,ix)
       enddo
      endif
*   ..write pfny segment?
      if (strmas(idout,'pfny')) then
       write (nprt,'(/2x,a6)') '+pfny'
*    ..write reduced fnay
       write (nprt,107) (is,is=0,ns-1)
       do iy = 0, ny
        write (nprt,'(2x,i4,1p,(t10,10e11.3))')
     &    iy, (rfnay(nx,iy,is),is=0,ns-1)
       enddo
      endif
*   ..write pfzy segment?
      if (strmas(idout,'pfzy')) then
       write (nprt,'(/2x,a6)') '+pfzy'
*   ..write reduced fney, fniy, fhey, fhiy, fchy
       write (nprt,108)
       do iy = 0, ny
        write (nprt,'(2x,i4,1p,(t10,5e11.3))')
     &    iy, rfney(nx,iy), rfniy(nx,iy),
     &        rfhey(nx,iy), rfhiy(nx,iy), rfchy(nx,iy)
       enddo
      endif
*   ..write pzso segment?
      if (strmas(idout,'pzso')) then
       write (nprt,'(/2x,a6)') '+pzso'
*    ..write state on south boundary
       write (nprt,111)
       do ix = -1, nx
        write (nprt,'(2x,i4)') ix
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(ix,-1)/ev, ti(ix,-1)/ev, ne(ix,-1), po(ix,-1),
     &    (na(ix,-1,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(ix,0)/ev, ti(ix,0)/ev, ne(ix,0), po(ix,0),
     &    (na(ix,0,is),is=0,ns-1)
       enddo
      endif
*   ..write pzno segment?
      if (strmas(idout,'pzno')) then
       write (nprt,'(/2x,a6)') '+pzno'
*    ..write state on north boundary
       write (nprt,112)
       do ix = -1, nx
        write (nprt,'(2x,i4)') ix
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(ix,ny)/ev, ti(ix,ny)/ev, ne(ix,ny), po(ix,ny),
     &    (na(ix,ny,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(ix,ny-1)/ev, ti(ix,ny-1)/ev, ne(ix,ny-1), po(ix,ny-1),
     &    (na(ix,ny-1,is),is=0,ns-1)
       enddo
      endif
*   ..write pzwe segment?
      if (strmas(idout,'pzwe')) then
       write (nprt,'(/2x,a6)') '+pzwe'
*    ..write state on west boundary
       write (nprt,113)
       do iy = -1, ny
        write (nprt,'(2x,i4)') iy
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(-1,iy)/ev, ti(-1,iy)/ev, ne(-1,iy), po(-1,iy),
     &    (na(-1,iy,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(0,iy)/ev, ti(0,iy)/ev, ne(0,iy), po(0,iy),
     &    (na(0,iy,is),is=0,ns-1)
       enddo
      endif
*   ..write pzea segment?
      if (strmas(idout,'pzea')) then
       write (nprt,'(/2x,a6)') '+pzea'
*    ..write state on east boundary
       write (nprt,114)
       do iy = -1, ny
        write (nprt,'(2x,i4)') iy
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(nx,iy)/ev, ti(nx,iy)/ev, ne(nx,iy), po(nx,iy),
     &    (na(nx,iy,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    te(nx-1,iy)/ev, ti(nx-1,iy)/ev, ne(nx-1,iy), po(nx-1,iy),
     &    (na(nx-1,iy,is),is=0,ns-1)
       enddo
      endif
*   ..write blank line
      write (nprt,'()')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  103 format (2x,a2,'(',i3,'):',2x,1pe11.2,' at (',i4,',',i4,');',
     &  2x,1pe11.2,' at (',i4,',',i4,').')
  104 format (2x,a2,'     :',2x,1pe11.2,' at (',i4,',',i4,');',
     &  2x,1pe11.2,' at (',i4,',',i4,').')
  105 format (/'  Flux of particles fnax(ix,,is),',
     &  ' integrated over y as function of x:'/
     &  t4, 'ix', (t9,10(5x,'is=',i3:)))
  106 format (/'  Flux of particles, heat, current',
     &  ' integrated over y as function of x:'/
     &  t4, 'ix', t14, 'fnex', t25, 'fnix', t36, 'fhex', t47, 'fhix',
     &  t58, 'fchx')
  107 format (/'  Flux of particles fnay(,iy,is),',
     &  ' integrated over x as function of y:'/
     &  t4, 'iy', (t9,10(5x,'is=',i3:)))
  108 format (/'  Flux of particles, heat, current',
     &  ' integrated over x as function of y:'/
     &  t4, 'iy', t14, 'fney', t25, 'fniy', t36, 'fhey', t47, 'fhiy',
     &  t58, 'fchy')
  111 format (/'  State of the plasma on the south boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)')
  112 format (/'  State of the plasma on the north boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)')
  113 format (/'  State of the plasma on the west boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)')
  114 format (/'  State of the plasma on the east boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)'/
     &  t6, 'te', t17, 'ti', t28, 'ne', t39, 'po', t50, 'na( 0)',
     &  t61, 'na(..)')

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: Also write total energy fluxes. Total energy flux has
*     factor (5/2) in the convective term and includes parallel kinetic
*     energy and recombination energy.

*-----------------------------------------------------------------------
*.end b2wfpj

      end
