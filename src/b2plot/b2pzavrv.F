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
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,i,j,k1,k0
      character*16 stmp
      logical streql
      external streql
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
 10   continue
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
c
      open(lun,file='zaveraev.dat')
      write (lun,710)
  710 format(9h   ix  iy,5x,5hR (m),7x,5hZ (m),5x,10hnH (m**-3),
     ,       5x,2hne,8x,7hTe (eV),5x,7hTi (eV),
     ,       6x,5hnHe++,7x,4hnHe+,8x,3hnHe/)
c
      do ix=0,nx-1
        do iy=0,ny-1
          write (lun,700) ix,iy,cr(ix,iy),cz(ix,iy),
     ,         dab2(ix+1,iy+1,j_h_at,1),
     ,         ne(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,
     ,         na(ix,iy,k1+1),na(ix,iy,k1),dab2(ix+1,iy+1,k0,1)
  700     format(1x,1p,2i4,48e12.3)
        end do
      end do
      close (lun)
c
      return
      end
