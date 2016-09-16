*-----------------------------------------------------------------------
*.specification

      subroutine b2wfpj (nprt, idout, nx, ny, ns, na, ua, te, ti, po,
     &  fna, fhe, fhi, fch, ne, ni, pr, ro, fne, fni, wrk0)
      use b2mod_types
      use b2mod_rates
      use b2mod_constants
      use b2mod_b2cmpa
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
     *   ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), pr(-1:nx,-1:ny),
     &  ro(-1:nx,-1:ny), fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)
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
      integer n2, ix, iy, is, ix0, iy0, ix1, iy1
*   ..procedures
      real (kind=R8) ::
     *   ssum, rfnax, rfnay, rfnix, rfniy, rfnex, rfney,
     &  rfhex, rfhey, rfhix, rfhiy, rfchx, rfchy
      logical strmas
      external subini, subend, xertst, ssum, strmas
      external b2xvsg, b2xvff, b2xvps, b2xxmm, b2xppb
      rfnax(ix,is) = ssum(ny+2,fna(ix,-1,0,is),nx+2)
      rfnay(iy,is) = ssum(nx+2,fna(-1,iy,1,is),1)
      rfhex(ix) = ssum(ny+2,fhe(ix,-1,0),nx+2)
      rfhey(iy) = ssum(nx+2,fhe(-1,iy,1),1)
      rfhix(ix) = ssum(ny+2,fhi(ix,-1,0),nx+2)
      rfhiy(iy) = ssum(nx+2,fhi(-1,iy,1),1)
      rfchx(ix) = ssum(ny+2,fch(ix,-1,0),nx+2)
      rfchy(iy) = ssum(nx+2,fch(-1,iy,1),1)
      rfnex(ix) = ssum(ny+2,fne(ix,-1,0),nx+2)
      rfney(iy) = ssum(nx+2,fne(-1,iy,1),1)
      rfnix(ix) = ssum(ny+2,fni(ix,-1,0),nx+2)
      rfniy(iy) = ssum(nx+2,fni(-1,iy,1),1)

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
*   ..test sign of ni, pr, ro
      call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
      call b2xvsg (n2, pr, 1, 'pr', '.gt.')
      call b2xvsg (n2, ro, 1, 'ro', '.gt.')
*   ..test edge values of fne, fni
      call b2xvff (nx, ny, fne, 'fne')
      call b2xvff (nx, ny, fni, 'fni')

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
     &    ix, (rfnax(ix,is),is=0,ns-1)
       enddo
      endif
*   ..write pfzx segment?
      if (strmas(idout,'pfzx')) then
       write (nprt,'(/2x,a6)') '+pfzx'
*    ..write reduced fnex, fnix, fhex, fhix, fchx
       write (nprt,106)
       do ix = 0, nx
        write (nprt,'(2x,i4,1p,(t10,5e11.3))')
     &    ix, rfnex(ix), rfnix(ix), rfhex(ix), rfhix(ix), rfchx(ix)
       enddo
      endif
*   ..write pfny segment?
      if (strmas(idout,'pfny')) then
       write (nprt,'(/2x,a6)') '+pfny'
*    ..write reduced fnay
       write (nprt,107) (is,is=0,ns-1)
       do iy = 0, ny
        write (nprt,'(2x,i4,1p,(t10,10e11.3))')
     &    iy, (rfnay(iy,is),is=0,ns-1)
       enddo
      endif
*   ..write pfzy segment?
      if (strmas(idout,'pfzy')) then
       write (nprt,'(/2x,a6)') '+pfzy'
*   ..write reduced fney, fniy, fhey, fhiy, fchy
       write (nprt,108)
       do iy = 0, ny
        write (nprt,'(2x,i4,1p,(t10,5e11.3))')
     &    iy, rfney(iy), rfniy(iy), rfhey(iy), rfhiy(iy), rfchy(iy)
       enddo
      endif
*   ..write pzso segment?
      if (strmas(idout,'pzso')) then
       write (nprt,'(/2x,a6)') '+pzso'
*    ..write state on south boundary
       write (nprt,111)
       do ix = -1, nx
        write (nprt,'(2x,i4)') ix
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
     &    te(ix,-1)/ev, ti(ix,-1)/ev, ne(ix,-1), po(ix,-1),
     &    (na(ix,-1,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
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
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
     &    te(ix,ny)/ev, ti(ix,ny)/ev, ne(ix,ny), po(ix,ny),
     &    (na(ix,ny,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
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
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
     &    te(-1,iy)/ev, ti(-1,iy)/ev, ne(-1,iy), po(-1,iy),
     &    (na(-1,iy,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
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
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
     &    te(nx,iy)/ev, ti(nx,iy)/ev, ne(nx,iy), po(nx,iy),
     &    (na(nx,iy,is),is=0,ns-1)
        write (nprt,'(1p,1x,4e11.3,2x,(t48,6e11.3))')
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
  104 format (4x,a3,' :',2x,1pe11.2,' at (',i4,',',i4,');',
     &  2x,1pe11.2,' at (',i4,',',i4,').')
  105 format (/'  Flux of particles fnax(ix,,is),',
     &  ' integrated over y as function of x:'/
     &  t4, 'ix', (t8,10(6x,'is=',i3:)))
  106 format (/'  Flux of particles, heat, current, and energy',
     &  ' integrated over y as function of x:'/
     &  t4, 'ix', t15, 'fnex', t26, 'fnix', t37, 'fhex', t48, 'fhix',
     &  t59, 'fchx', t70, 'feex', t81, 'feix')
  107 format (/'  Flux of particles fnay(,iy,is),',
     &  ' integrated over x as function of y:'/
     &  t4, 'iy', (t8,10(6x,'is=',i3:)))
  108 format (/'  Flux of particles, heat, current, and energy',
     &  ' integrated over x as function of y:'/
     &  t4, 'iy', t15, 'fney', t26, 'fniy', t37, 'fhey', t48, 'fhiy',
     &  t59, 'fchy', t70, 'feey', t81, 'feiy')
  111 format (/'  State of the plasma on the south boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)')
  112 format (/'  State of the plasma on the north boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)')
  113 format (/'  State of the plasma on the west boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)')
  114 format (/'  State of the plasma on the east boundary:'/
     &  '  (second line for te, ti, etc. gives state just',
     &  ' inside the boundary.)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)'/
     &  t8, 'te', t19, 'ti', t30, 'ne', t41, 'po', t51, 'na( 0)',
     &  t63, 'na(..)')

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: Also write total energy fluxes. Total energy flux has
*     factor (5/2) in the convective term and includes parallel kinetic
*     energy and recombination energy.

*-----------------------------------------------------------------------
*.end b2wfpj

      end
