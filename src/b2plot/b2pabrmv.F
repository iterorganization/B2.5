      subroutine b2pabrmv(lun)
c
c----------------------------------------------------------------------
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density)
c***  in these points - for V. Abramov
c***
c*** Output files:
c***
c***  abramov.dat with the data for V. Abramov
c----------------------------------------------------------------------
      use b2mod_types
      use b2us_data
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_elements
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      implicit none
      integer lun
c
c*** Local variables
c
      integer iCv,is
      integer, allocatable :: mz_hlp(:)
      real (kind=R8) :: vr,vz,vt,u
      character*2, allocatable :: hlp_al(:)
c----------------------------------------------------------------------
c*** Write the header
c
      allocate(mz_hlp(0:ns-1),hlp_al(0:ns-1))
      do is=0,ns-1
        mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is)=is_codes(is)
      end do
c
      open(lun,file='abramov.dat')
  710 format(a6,8x,a5,7x,a5,4x,a10,4x,
     ,       a7,5x,a7,7x,a2,10x,a2,10x,a2,
     ,       10x,a2,10x,a3/)
      write (lun,710) '   iCv','R (m)','Z (m)','ne (m**-3)',
     ,                'Te (eV)','Ti (eV)','vr','vz','vt','nH','nH2'
c
      do iCv=1,mpg%nCv
#ifdef WG_TODO
        u=atan2(cz(ix+1,iy)-cz(ix-1,iy),cr(ix+1,iy)-cr(ix-1,iy))
#endif
        vr=cos(u)*abs(geo%cvBb(iCv,0))/geo%cvBb(iCv,3)
        vz=sin(u)*abs(geo%cvBb(iCv,0))/geo%cvBb(iCv,3)
        vt=sqrt(max(0.0_R8,1.0_R8-(geo%cvBb(iCv,0)/geo%cvBb(iCv,3))**2))
        write (lun,700) iCv,geo%cvX(iCv),geo%cvY(iCv),state%dv%ne(iCv),
     ,                  state%pl%te(iCv)/ev,state%pl%ti(iCv)/ev,vr,vz,vt
     ,                 ,dab2(iCv,j_h_at,1),dmb2(iCv,1,1)
  700   format(1x,1p,i6,10e12.3)
      end do
      close (lun)
      deallocate(mz_hlp,hlp_al)
c
      return
      end subroutine b2pabrmv
