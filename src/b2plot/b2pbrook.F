      subroutine b2pbrooks(lun)
c
c----------------------------------------------------------------------
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, drift velocity)
c***  in these points - for J. Brooks
c***
c*** Output files:
c***
c***  brooks.dat with the data for J. Brooks
c----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_plasma
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_constants
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,is,j
      integer, allocatable :: mz_hlp(:)
      real (kind=R8) :: vr,vz,vt,u
      character hlp_frmr*600
      character*2, allocatable :: hlp_al(:)
c----------------------------------------------------------------------
c*** Write the header
c
      allocate(mz_hlp(0:ns-1),hlp_al(0:ns-1))
      do is=0,ns-1
        mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is)=is_codes(is)
        if (hlp_al(is)(2:2).eq.' ') hlp_al(is)(2:2) = '_'
      end do
c
      open(lun,file='brooks.dat')
      write (lun,*) 'ix, iy   : x and y indices in the grid;'
      write (lun,*) 'R , Z    : R and Z coordinates of the cell ',
     ,                                                       'centres;'
      write (lun,*) 'ne       : electron density;'
      write (lun,*) 'Te, Ti   : electron and ion temperatures;'
      write (lun,*) 'vr,vz,vt : R, Z, and toroidal cosines of B;'
      write (lun,*) 'n...     : densities of the ion species;'
      write (lun,*) 'vp...    : parallel velocities of the ion species;'
      write (lun,*)
      write(hlp_frmr,710) ns
  710 format('(9h   ix  iy,5x,5hR (m),7x,5hZ (m),4x,10hne (m**-3),4x,',
     ,       '7hTe (eV),5x,7hTi (eV),7x,2hvr,10x,2hvz,10x,2hvt,5x,',
     ,                      i2,'(4x,1hn,a2,i2.2,6x,2hvp,a2,i2.2,3x)/)')
      write (lun,hlp_frmr) ((hlp_al(is),mz_hlp(is),j=1,2),is=0,ns-1)
c
      do ix=0,nx-1
        do iy=0,ny-1
          u=atan2(cz(ix+1,iy)-cz(ix-1,iy),cr(ix+1,iy)-cr(ix-1,iy))
          vr=cos(u)*abs(bb(ix,iy,0))/bb(ix,iy,3)
          vz=sin(u)*abs(bb(ix,iy,0))/bb(ix,iy,3)
          vt=sqrt(max(0.0_R8,1.0_R8-(bb(ix,iy,0)/bb(ix,iy,3))**2))
          write (lun,700) ix,iy,cr(ix,iy),cz(ix,iy),ne(ix,iy),
     ,                    te(ix,iy)/ev,ti(ix,iy)/ev,vr,vz,vt,
     ,                    (na(ix,iy,is),ua(ix,iy,is),is=0,ns-1)
  700     format(1x,1p,2i4,48e12.3)
        end do
      end do
      close (lun)
      deallocate(mz_hlp,hlp_al)
      return
c
      end subroutine b2pbrooks
