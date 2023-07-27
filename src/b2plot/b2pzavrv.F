      subroutine b2pzavrv(lun,mpg,gm,st)
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
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      use b2us_map
      use b2us_geo
      use b2us_plasma
      implicit none
      integer, intent(in) :: lun
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2State), intent(in) :: st
c
c*** Local variables
c
      integer iCv,i,j,in,k1,k0
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
      write(hlp_frmr,'(a,i2,a)') '(a6,3x,a5,7x,a5,5x,',
     . in,'(a10,2x),3x,a2,8x,a7,5x,a7,6x,a5,7x,a4,8x,a3)'
      write (lun,hlp_frmr)
     .  '   iCv','R (m)','Z (m)',
     . ('n'//hlp_at(i)//' (m**-3)',i=1,in),
     .  'ne','Te (eV)','Ti (eV)','nHe++','nHe+','nHe'
c
      do iCv=1,mpg%nCi
        i=0
        do j=1,3
          if (j_h_at(j).ne.0) then
            i=i+1
            hlp_f(i)=dab2(iCv,j_h_at(j),1)
          end if
        end do
        write (lun,700) iCv,gm%cvX(iCv),gm%cvY(iCv),
     ,        (hlp_f(i),i=1,in),
     ,         st%dv%ne(iCv),st%pl%te(iCv)/ev,st%pl%ti(iCv)/ev,
     ,         st%pl%na(iCv,k1+1),st%pl%na(iCv,k1),
     ,         dab2(iCv,k0,1)
      end do
      close (lun)
      deallocate(hlp_f)
c
  700 format(1x,1p,i6,48e12.3)
      return
      end subroutine b2pzavrv
