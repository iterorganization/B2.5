      subroutine b2pblmtr(lun)

c version : 18.01.2002 00:19

c=======================================================================
c***  Produces a file with coordinates of the grid cell centers and
c***  the radiation intensity (total and components)
c***  in these points - for bolometry optimisation
c***
c*** Output files:
c***
c***  radsources.dat with the data
c***
c*** cell count starts from the target PFR corner
c-----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_b2cmpa
      use b2mod_sources
      use b2mod_elements
      use b2mod_neutrals_namelist
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,is,i,j
      real (kind=R8), allocatable :: hlp_rad(:)
      real (kind=R8) :: hlp, hlp_brm
      character hlp_frmr*600
      character*2, allocatable :: hlp_al(:)
c-----------------------------------------------------------------------
c*** The plasma data at the cell centers
c
      allocate(hlp_al(nspecies),hlp_rad(nspecies))
      do j=1,nspecies !{
        hlp_al(j)=elements(nfluids(j))
        if(nfluids(j).eq.1) then !{
          if(nint(am(eb2spcr(j))).eq.2) then !{
            hlp_al(j)='D'
          else if(nint(am(eb2spcr(j))).eq.3) then !}{
            hlp_al(j)='T'
          end if !}
        end if !}
      end do !}

      open(lun,file='radsources.dat')
      write (lun,'(2a)') '# ix, iy  : poloidal and radial indices in ',
     ,                                                      'the grid;'
      write (lun,'(2a)') '# R , Z   : R and Z coordinates of the ',
     ,                                              'cell centres [m];'
      write (lun,'(a)') '# rad_tot : total radiation [MW/m^-3];'
      write (lun,'(a)') '# rad_*   : radiation by species [MW/m^-3];'
      write (lun,'(a)') '# rad_brm : Bremsstrahlung [MW/m^-3];'
      write (lun,'(a)')

      write(hlp_frmr,710) nspecies
  710 format('(9h#  ix  iy,7x,1hR,11x,1hZ,7x,7hrad_tot,2x,',
     ,  i2,'(4x,4hrad_,a2,2x),3x,7hrad_brm,2x)')
      write (lun,hlp_frmr) (hlp_al(i),i=1,nspecies)

      do iy=0,ny-1 !{
        write (lun,*)
        do ix=0,nx-1 !{
          hlp_brm=0.0_R8
          hlp_rad=0.0_R8
          do is=0,ns-1
            i=b2espcr(is)
            hlp_brm=hlp_brm+rqbrm(ix,iy,is)
            hlp_rad(i)=hlp_rad(i)+rqrad(ix,iy,is)
          end do !}
          hlp=hlp_brm
          do i=1,nspecies
            hlp=hlp+hlp_rad(i)
          end do
          write (lun,700) ix,iy,cr(ix,iy),cz(ix,iy),
     ,      -1.0e-6_R8*hlp/vol(ix,iy),
     ,     (-1.0e-6_R8*hlp_rad(i)/vol(ix,iy),i=1,nspecies),
     ,      -1.0e-6_r8*hlp_brm/vol(ix,iy)
  700     format(1x,1p,2i4,48e12.4)
        end do !}
      end do !}
      close (lun)
      deallocate(hlp_al,hlp_rad)
      return
c=======================================================================
      end
