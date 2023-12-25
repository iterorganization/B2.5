      subroutine b2poleg(lun,ig)
c
c----------------------------------------------------------------------
c***
c***  Produces a file with plasma parameters (density, temperature,
c***  parallel flow velocity, parallel heat flux) against the
c***  distance along the magnetic field along field line index
c***  ig - for Oleg
c***
c*** Output files:
c***
c***  oleg.dat with the data for Oleg
c----------------------------------------------------------------------
      use b2mod_types
      use b2mod_geo
      use b2mod_rates
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_sources
      use b2mod_external
      use b2mod_constants
      implicit none
      integer ig
      integer lun
c
c*** Local variables
c
      integer ix, k
      real (kind=R8) :: x, y, z, hlp_q, hlp_qp, hlp_si, hlp_je
c----------------------------------------------------------------------
c*** Check the ig value
c
      if(ig.lt.0 .or. ig.gt.ny) then
        print *,'*** b2poleg: wrong value of ig = ',ig
        return
      end if
c
c*** Write the header
c
      open(lun,file='oleg.dat')
      write (lun,710) 'x','Xpol','ne','Te','qpol_e','qpar_e',
     .                'S_ion','Zeff','jpar_e'
  710 format(6x,a1,9x,a4,9x,a2,10x,a2,8x,a6,6x,a6,
     ,       6x,a5,8x,a4,7x,a6/)
c
      x=0.0_R8
      z=0.0_R8
      do ix=0,nx-1
        x=x+0.5_R8*hx(ix,ig)
        y=0.5_R8*hx(ix,ig)/pit(ix,ig)
        z=z+y
        hlp_q=fhe(ix,ig,0,0)/area(ix,ig)
        hlp_qp=hlp_q/pit(ix,ig)
        hlp_si=0.0_R8
        hlp_je=0.0_R8
        do k=0,ns-1
          hlp_je=hlp_je+fna(ix,ig,0,0,k)*rza(ix,ig,k)
          hlp_si=hlp_si+b2stbr_sna(ix,ig,k)
        end do
        do k=0,ns_ext-1
          hlp_je=hlp_je+fa_ext(ix,ig,0,0,k)*za_ext(ix,ig,k)
        end do
        hlp_je=hlp_je/(area(ix,ig)*pit(ix,ig))
        hlp_si=hlp_si/vol(ix,ig)
c
        write(lun,700) z,x,ne(ix,ig),te(ix,ig)/ev,hlp_q,hlp_qp,hlp_si,
     ,                                             zeff(ix,ig),hlp_je
        x=x+0.5_R8*hx(ix,ig)
        z=z+y
  700   format(1p,48e12.4)
      end do
      close (lun)
      return
c
      contains

      function area(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real(kind=R8) :: area
      area = gs(ix,iy,0)*qc(ix,iy,0)
      return
      end function area

      end subroutine b2poleg
