      subroutine b2krieger(lun)

c version : 16.07.2004 16:20

c=======================================================================
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, fluid velocity)
c***  - for K. Krieger
c***
c*** Output files:
c***
c***  krieger.dat with the data
c***
c*** cell count starts from the target PFR corner
c-----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_b2cmpa
      use b2mod_plasma
      use b2mod_elements
      use b2mod_constants
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,is,i,j,k,l
      integer, allocatable :: mz_hlp(:)
      real (kind=R8), allocatable :: n_hlp(:),v_hlp(:)
      character hlp_frmr*600, hlp_frm*200
      character*2, allocatable :: hlp_al(:)
c-----------------------------------------------------------------------
c*** The plasma data at the cell centers
c
      allocate(mz_hlp(nspecies+ns),hlp_al(nspecies+ns))
      allocate(n_hlp(nspecies),v_hlp(nspecies))
      do is=0,ns-1 !{
        mz_hlp(nspecies+1+is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(nspecies+1+is)=is_codes(is)
        hlp_al(b2espcr(is))=hlp_al(nspecies+1+is)
        mz_hlp(b2espcr(is))=99
      end do !}
      do is=1,nspecies+ns !{
        if (hlp_al(is)(2:2).eq.' ') hlp_al(is)(2:2) = '_'
      end do !}

      open(lun,file='krieger.dat')
      write (lun,'(2a)')
     ,  '# ix, iy   : poloidal and radial indices in the grid;'
      write (lun,'(2a)')
     ,  '# R , Z    : R and Z coordinates of the cell centres [m];'
      write (lun,'(2a)')
     ,  '# pit      : pitch angle;'
      write (lun,'(2a)')
     ,  '# ne       : electron density [m^-3];'
      write (lun,'(2a)')
     ,  '# Te, Ti   : electron and ion temperatures [eV];'
      write (lun,'(2a)')
     ,  '# n...     : densities of the ion species [m^-3];'
      write (lun,'(2a)')
     ,  '# vp...    : parallel velocities of the ion species [m/s]' 
      write (lun,'(2a)')
     ,  '# n...99 and vp...99 stand for the ',
     ,  'total density and average velocity of the ion species'
      write (lun,*)
      write (hlp_frmr,'(a,i3,a)')
     ,  '(a9,7x,a1,11x,a1,10x,a3,9x,a2,11x,a2,10x,a2,4x,',
     ,    nspecies+ns,'(4x,a3,i2.2,6x,a4,i2.2,3x)/)'
      write (lun,hlp_frmr)
     ,  '#  ix  iy','R','Z','pit','ne','Te','Ti',
     ,  ('n'//hlp_al(i),mz_hlp(i),
     ,  'vp'//hlp_al(i),mz_hlp(i),i=1,nspecies+ns)
      write (hlp_frm,'(a,i3,a)')
     ,  '(1x,1p,2i4,',6+2*(nspecies+ns),'e12.4)'

      do iy=0,ny-1 !{
        write (lun,*)
        do ix=0,nx-1 !{
          k=0
          do j=1,nspecies !{
            n_hlp(j)=0.0_R8
            v_hlp(j)=0.0_R8
            do l=1,nfluids(j)+1 !{
              n_hlp(j)=n_hlp(j)+na(ix,iy,k)
              v_hlp(j)=v_hlp(j)+na(ix,iy,k)*ua(ix,iy,k)
              k=k+1
            end do !}
            v_hlp(j)=v_hlp(j)/n_hlp(j)
          end do !}
          write (lun,hlp_frm) ix,iy,cr(ix,iy),cz(ix,iy),pit(ix,iy),
     ,                        ne(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,
     ,                       (n_hlp(j),v_hlp(j),j=1,nspecies),
     ,                       (na(ix,iy,k),ua(ix,iy,k),k=0,ns-1)
        end do !}
      end do !}
      close (lun)
      deallocate(mz_hlp,hlp_al)
      deallocate(n_hlp,v_hlp)
      return
c=======================================================================
      end subroutine b2krieger
