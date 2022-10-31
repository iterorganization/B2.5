      subroutine b2phortn(lun)
c
c----------------------------------------------------------------------
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density)
c***  in these points - for L. Horton
c***
c*** Output files:
c***
c***  horton.dat with the data for L. Horton
c----------------------------------------------------------------------
      use b2mod_types
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
      integer i,ix,is,in,im,nf
#ifdef WG_TODO
      integer iy
#endif
      integer, allocatable :: mz_hlp(:),imol(:)
      real (kind=R8), allocatable :: hlp_f(:)
      character hlp_frmr*600
      character*2, allocatable :: hlp_al(:),hlp_ml(:)
      character*1 hlp_at(3),isotope(3)
      logical hydrogen
      data isotope / 'H', 'D', 'T'/
c----------------------------------------------------------------------
c*** Write the header
c
      allocate(mz_hlp(0:ns-1),hlp_al(0:ns-1),hlp_ml(nmoli),imol(nmoli))
      do is=0,ns-1
        mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is)=is_codes(is)
      end do
      in=0
      do is=1,3
        if (j_h_at(is).ne.0) then
          in=in+1
          hlp_at(in)=isotope(is)
        end if
      end do
      im=0
      imol=0
      do is=1,nmoli
        if (l_h_mol(is,1)+l_h_mol(is,2)+l_h_mol(is,3).eq.2) then
          hydrogen=.true.
          do ix=1,natmi
            if (ix.eq.j_h_at(1)) cycle
            if (ix.eq.j_h_at(2)) cycle
            if (ix.eq.j_h_at(3)) cycle
            hydrogen = mlcmp(ix,is).eq.0
          end do
          if (hydrogen) then
            im=im+1
            hlp_ml(im)=textmn(is-1)(2:3)
            imol(im)=is
          end if
        end if
      end do
c
      open(lun,file='horton.dat')
      write (lun,*) 'ix, iy : x and y indices in the grid;'
      write (lun,*) 'R , Z  : R and Z coordinates of the cell ',
     ,                                                       'centres;'
      write (lun,*) 'ne     : electron density;'
      write (lun,*) 'Te, Ti : electron and ion temperatures;'
      if (j_h_at(1).ne.0)
     > write (lun,*) 'nH, TH : density and temperature of H atoms;'
      if (j_h_at(2).ne.0)
     > write (lun,*) 'nD, TD : density and temperature of D atoms;'
      if (j_h_at(3).ne.0)
     > write (lun,*) 'nT, TT : density and temperature of T atoms;'
      do is=1,im
       write (lun,*) 'n'//hlp_ml(is)//'    : density of '
     .                  //hlp_ml(is)//' molecules.'
      end do
      write (lun,*) '     All the densities in 1/m**3'
      write (lun,*)
      write (hlp_frmr,'(a,i2,a,i2,a)')
     .    '(a9,5x,a5,7x,a5,4x,a10,4x,a7,5x,a7,7x,',
     . in,'(a2,10x,a2,10x),',im,'(a3,9x))'
      write (lun,hlp_frmr)
     .  '   ix  iy','R (m)','Z (m)','ne (m**-3)','Te (eV)','Ti (eV)',
     . ('n'//hlp_at(i),'T'//hlp_at(i),i=1,in),
     . ('n'//hlp_ml(i),i=1,im)
c
      nf=2*in+im
      allocate(hlp_f(nf))
      write(hlp_frmr,'(a,i2,a)') '(1x,1p,2i4,',5+nf,'e12.3)'
#ifdef WG_TODO
      do ix=0,nx-1
        do iy=0,ny-1
          is=0
          do i=1,3
            if (j_h_at(i).ne.0) then
              hlp_f(is+1)=dab2(ix+1,iy+1,j_h_at(i),1)
              hlp_f(is+2)=tab2(ix+1,iy+1,j_h_at(i),1)/ev
              is=is+2
            end if
          end do
          do i=1,im
            hlp_f(is+1)=dmb2(ix+1,iy+1,imol(im),1)
            is=is+1
          end do
          write (lun,hlp_frmr) ix,iy,cr(ix,iy),cz(ix,iy),ne(ix,iy),
     ,                         te(ix,iy)/ev,ti(ix,iy)/ev,
     ,                        (hlp_f(i),i=1,nf)
        end do
      end do
#endif
      close (lun)
      deallocate(mz_hlp,hlp_al,hlp_ml,hlp_f,imol)
c
      return
      end subroutine b2phortn
