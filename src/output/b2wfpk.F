*-----------------------------------------------------------------------
*.specification

      subroutine b2wfpk (nprt, idout, nx, ny, ns)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      character idout*(*)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WFPK writes supplementary information about the plasma state to
*     a formatted file.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit nprt some supplementary information
*     describing the plasma state. Unit nprt must be connected to a
*     formatted file. The extent of the output is specified by the
*     parameter idout.
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
*     gs - real*8 array, input.
*     See the calling routine for description.
*
*     qc - real*8 array, input.
*     See the calling routine for description.
*
*     fna, fhe, fhi, fch, fne, fni - real*8 array, input.
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
      integer ix, iy, is
*   ..procedures
      logical strmas
      external xertst, strmas
      external b2xvff

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wfpk')
*   ..test nprt, nx, ny
      call xertst (1.le.nprt, 'faulty argument nprt')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

* ..write the flux densities
*   ..write pzso segment?
      if (strmas(idout,'pzso')) then
       write (nprt,'(/2x,a6)') '+pzso'
*    ..write flux densities on south boundary
       write (nprt,111)
       do ix = -1, nx
        write (nprt,'(2x,i4)') ix
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    fhe(ix,0,1)/gs(ix,0,1), fhi(ix,0,1)/gs(ix,0,1),
     &    fne(ix,0,1)/gs(ix,0,1), fch(ix,0,1)/gs(ix,0,1),
     &   (fna(ix,0,1,is)/gs(ix,0,1),is=0,ns-1)
       enddo
      endif
*   ..write pzno segment?
      if (strmas(idout,'pzno')) then
       write (nprt,'(/2x,a6)') '+pzno'
*    ..write flux densities on north boundary
       write (nprt,112)
       do ix = -1, nx
        write (nprt,'(2x,i4)') ix
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    fhe(ix,ny,1)/gs(ix,ny,1), fhi(ix,ny,1)/gs(ix,ny,1),
     &    fne(ix,ny,1)/gs(ix,ny,1), fch(ix,ny,1)/gs(ix,ny,1),
     &   (fna(ix,ny,1,is)/gs(ix,ny,1),is=0,ns-1)
       enddo
      endif
*   ..write pzwe segment?
      if (strmas(idout,'pzwe')) then
       write (nprt,'(/2x,a6)') '+pzwe'
*    ..write flux densities on west boundary
       write (nprt,113)
       do iy = -1, ny
        write (nprt,'(2x,i4)') iy
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    fhe(0,iy,0)/gs(0,iy,0)/qc(0,iy), 
     &    fhi(0,iy,0)/gs(0,iy,0)/qc(0,iy),
     &    fne(0,iy,0)/gs(0,iy,0)/qc(0,iy), 
     &    fch(0,iy,0)/gs(0,iy,0)/qc(0,iy),
     &   (fna(0,iy,0,is)/gs(0,iy,0)/qc(0,iy),is=0,ns-1)
       enddo
      endif
*   ..write pzea segment?
      if (strmas(idout,'pzea')) then
       write (nprt,'(/2x,a6)') '+pzea'
*    ..write flux densities on east boundary
       write (nprt,114)
       do iy = -1, ny
        write (nprt,'(2x,i4)') iy
        write (nprt,'(1p,1x,4e11.3,2x,(t46,6e11.3))')
     &    fhe(nx,iy,0)/gs(nx,iy,0)/qc(nx,iy), 
     &    fhi(nx,iy,0)/gs(nx,iy,0)/qc(nx,iy),
     &    fne(nx,iy,0)/gs(nx,iy,0)/qc(nx,iy), 
     &    fch(nx,iy,0)/gs(nx,iy,0)/qc(nx,iy),
     &   (fna(nx,iy,0,is)/gs(nx,iy,0)/qc(nx,iy),is=0,ns-1)
       enddo
      endif
*   ..write blank line
      write (nprt,'()')

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.formats

  111 format (/'  Flux densities on the south boundary:'/
     &  t5, 'fhey/sy', t16, 'fhiy/sy', t27, 'fney/sy',
     &  t38,'fchy/sy', t49, 'fnay/sy', t60, 'fnay/sy')
  112 format (/'  Flux densities on the north boundary:'/
     &  t5, 'fhey/sy', t16, 'fhiy/sy', t27, 'fney/sy',
     &  t38,'fchy/sy', t49, 'fnay/sy', t60, 'fnay/sy')
  113 format (/'  Flux densities on the west boundary:'/
     &  t5, 'fhex/sx', t16, 'fhix/sx', t27, 'fnex/sx',
     &  t38,'fchx/sx', t49, 'fnax/sx', t60, 'fnax/sx')
  114 format (/'  Flux densities on the east boundary:'/
     &  t5, 'fhex/sx', t16, 'fhix/sx', t27, 'fnex/sx',
     &  t38,'fchx/sx', t49, 'fnax/sx', t60, 'fnax/sx')

*-----------------------------------------------------------------------
*.scribbles

*!!   Note: Also write total energy flux densities. Total energy
*     flux has factor (5/2) in the convective term and includes kinetic
*     energy and recombination energy.

*-----------------------------------------------------------------------
*.end b2wfpk

      end subroutine b2wfpk
