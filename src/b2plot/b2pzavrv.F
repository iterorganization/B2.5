      subroutine b2pzavrv(lun)
c
c----------------------------------------------------------------------
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density)
c***  in these points - for V. Zaveryaev
c***
c*** Output files:
c***
c***  zaveraev.dat with the data for V. Zaveryaev
c----------------------------------------------------------------------
      use b2mod_diag
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,i,j,in,k1,k0
      real (kind=R8), allocatable :: hlp_f(:)
      character*16 stmp
      character*256 hlp_frmr
      character*1 isotope(3), hlp_at(3)
      logical streql
      external streql
      data isotope / 'H', 'D', 'T'/
c----------------------------------------------------------------------
c*** Find helium
c
      k1=0
      k0=0
      do j=1,nspecies
        i=nfluids(j)
        if(i.eq.2) then
          k1=k1+1
          go to 10
        end if
        k1=k1+i+1
      end do
   10 continue
      do j=1,natmi
        print *,j,' :',textan(j-1)
        stmp=textan(j-1)
        call lower_case(stmp)
        if(streql(stmp,' he')) k0=j
      end do
      print *,'Helium (ion, atom): ',k1,k0
      if(k1.eq.0 .or. k0.eq.0) then
        write(*,*) 'b2pzavrv: no helium found: ',k1,k0
        return
      end if
      in=0
      do j=1,3
        if (j_h_at(j).ne.0) then
          in=in+1
          hlp_at(in)=isotope(j)
        end if
      end do
c
      allocate(hlp_f(in))
      open(lun,file='zaveraev.dat')
      write(hlp_frmr,'(a,i2,a)') '(a9,5x,a5,7x,a5,5x,',
     . in,'(a10,2x),3x,a2,8x,a7,5x,a7,6x,a5,7x,a4,8x,a3)'
      write (lun,hlp_frmr)
     .  '   ix  iy','R (m)','Z (m)',
     . ('n'//hlp_at(i)//' (m**-3)',i=1,in),
     .  'ne','Te (eV)','Ti (eV)','nHe++','nHe+','nHe'
c
      do ix=0,nx-1
        do iy=0,ny-1
          i=0
          do j=1,3
            if (j_h_at(j).ne.0) then
              i=i+1
              hlp_f(i)=dab2(ix+1,iy+1,j_h_at(j),1)
            end if
          end do
          write (lun,700) ix,iy,cr(ix,iy),cz(ix,iy),
     ,        (hlp_f(i),i=1,in),
     ,         ne(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,
     ,         na(ix,iy,k1+1),na(ix,iy,k1),dab2(ix+1,iy+1,k0,1)
  700     format(1x,1p,2i4,48e12.3)
        end do
      end do
      close (lun)
      deallocate(hlp_f)
c
      return
      end subroutine b2pzavrv
