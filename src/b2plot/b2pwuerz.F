      subroutine b2pwuerz(lun)

c version : 22.09.2000 00:03

c=======================================================================
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, fluid velocity)
c***  in the outer divertor - for H. Wuerz
c***
c*** Output files:
c***
c***  wuerz-d.dat with the data
c***  wuerz-g.dat with the grid corners
c***
c*** cell count starts from the target PFR corner
c-----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_elements
      use b2mod_indirect
      use b2mod_constants
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,is,i,j,k
      integer, allocatable :: mz_hlp(:)
      real (kind=R8) :: vr,vz,vt,u
      real (kind=R8) :: pit
      character hlp_frmr*600
      character*2, allocatable :: hlp_al(:)
c-----------------------------------------------------------------------
      pit(ix,iy)=abs(bb(ix,iy,0))/bb(ix,iy,3)
c-----------------------------------------------------------------------
c*** The plasma data at the cell centers
c
      allocate(mz_hlp(0:ns-1),hlp_al(0:ns-1))
      do is=0,ns-1
        mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is)=is_codes(is)
      end do
c
      open(lun,file='wuerz-d.dat')
      write (lun,'(a)')
     ,  '# ix, iy   : poloidal and radial indices in the grid;'
      write (lun,'(a)')
     ,  '# R , Z    : R and Z coordinates of the cell centres [m];'
      write (lun,'(a)')
     ,  '# pit      : pitch angle;'
      write (lun,'(a)')
     ,  '# ne       : electron density [m^-3];'
      write (lun,'(a)')
     ,  '# Te, Ti   : electron and ion temperatures [eV];'
      write (lun,'(a)')
     ,  '# n...     : densities of the ion species [m^-3];'
      write (lun,'(a)')
     ,  '# vp...    : parallel velocities of the ion species [m/s]' 
      write (lun,*)
      write(hlp_frmr,710) ns
  710 format('(9h#  ix  iy,7x,1hR,11x,1hZ,10x,3hpit,9x,',
     ,  '2hne,11x,2hTe,10x,2hTi,4x,',
     ,  i3,'(4x,1hn,a2,i2.2,6x,2hvp,a2,i2.2,3x)/)')
      write (lun,hlp_frmr) ((hlp_al(i),mz_hlp(i),j=1,2),i=0,ns-1)

      do iy=0,ny-1
        i=0
        write (lun,*)
        do ix=nx-1,rightcut(1),-1
          write (lun,700) i,iy,cr(ix,iy),cz(ix,iy),pit(ix,iy),
     ,                    ne(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,
     ,                    (na(ix,iy,k),ua(ix,iy,k),k=0,ns-1)
  700     format(1x,1p,2i4,48e12.4)
          i=i+1
        end do
      end do
      close (lun)
c-----------------------------------------------------------------------
c*** The grid geometry

      open(lun,file='wuerz-g.dat')
      write (lun,'(2a)')
     ,  '# Coordinates of the grid cell corners'
      write (lun,'(2a)')
     ,  '# ix,iy : poloidal and radial indices in the grid;'
      write (lun,'(2a)')
     ,  '# R , Z : R and Z coordinates of the cell ',
     ,  'corners along corresponding poloidal surfaces'
      write (lun,'(a,6x,a1,11x,a1)') '#  ix  iy','R','Z'

      do iy=0,ny !{
        i=0
        write (lun,*)
        do j=nx-1,rightcut(1),-1 !{
          ix=j
          if(j.eq.rightcut(1).and.iy.le.topcut(1)) ix=leftcut(1)-1
          write (lun,'(1x,1p,2i4,2e12.4)')
     .     i,iy,crx(ix,iy,1),cry(ix,iy,1)
          i=i+1
        end do !}
      end do !}
      close (lun)
      deallocate(mz_hlp,hlp_al)
      return
c=======================================================================
      end
