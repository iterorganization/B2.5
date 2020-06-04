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
      use b2mod_plasma
      use b2mod_elements
      use b2mod_constants
      use b2mod_neutrals_namelist
      implicit none
      integer lun
c
c*** Local variables
c
      integer ix,iy,is,i,j
      integer use_eirene
      real (kind=R8), allocatable :: hlp_rad(:), hlp_ni(:)
      real (kind=R8) :: hlp, hlp_brm
      character hlp_frmr*600, species_name*13
      character*6, allocatable :: hlp_al(:)
      character*10, allocatable :: hlp_is(:)
      character*9, allocatable :: hlp_mn(:), hlp_in(:)
      data use_eirene /0/
c-----------------------------------------------------------------------
c*** The plasma data at the cell centers
c
      call ipgeti ('b2mndr_eirene', use_eirene)
      allocate(hlp_al(nspecies),hlp_rad(nspecies))
      if (nmoli.gt.0) allocate(hlp_mn(0:nmoli-1))
      if (nioni.gt.0) allocate(hlp_in(0:nioni-1))
      allocate(hlp_ni(0:ns-1))
      allocate(hlp_is(0:ns-1))
      do j=1,nspecies !{
        hlp_al(j)='rad_'//elements(nfluids(j))
        if(nfluids(j).eq.1) then !{
          if(nint(am(eb2spcr(j))).eq.2) then !{
            hlp_al(j)='rad_D '
          else if(nint(am(eb2spcr(j))).eq.3) then !}{
            hlp_al(j)='rad_T '
          end if !}
        end if !}
      end do !}
      do is=0,ns-1 !{
        call species(is,species_name,.false.)
        hlp_is(is)='n_'//species_name(2:9)
      enddo !}
      if (nmoli.gt.0) then !{
        do is=0,nmoli-1 !{
          hlp_mn(is)='n_'//textmn(is)(2:8)
        enddo !}
      endif !}
      if (nioni.gt.0) then !{
        do is=0,nioni-1 !{
          hlp_in(is)='n_'//textin(is)(2:8)
        enddo !}
      end if !}

      open(lun,file='radsources.dat')
      write (lun,'(2a)') '# ix, iy  : poloidal and radial indices in ',
     ,                                                      'the grid;'
      write (lun,'(2a)') '# R , Z   : R and Z coordinates of the ',
     ,                                              'cell centers [m];'
      write (lun,'(a)') '# Te : electron temperature [eV];'
      write (lun,'(a)') '# Ti : ion temperature [eV];'
      write (lun,'(a)') '# ne : electron density [m^⁻3];'
      write (lun,'(a)') '# n_* : species density [m^⁻3];'
      write (lun,'(a)') '# rad_tot : total radiation [MW/m^3];'
      write (lun,'(a)')
     ,       '# rad_*   : line radiation by species [MW/m^3];'
      write (lun,'(a)') '# rad_brm : bremsstrahlung [MW/m^3];'
      write (lun,'(a)')

      if (nmoli.gt.0 .and. nioni.gt.0) then
        write(hlp_frmr,'(a,i3,a,i2,a,i2,a,i3,a)')
     ,   '(a9,7x,a1,11x,a1,10x,a2,10x,a2,10x,a2,7x,',
     .    ns,'(1x,a10,1x),',nmoli,'(1x,a9,2x),',
     ,    nioni,'(1x,a9,2x),a7,2x,',nspecies,'(4x,a6,2x),3x,a7)'
        write(lun,hlp_frmr) '#  ix  iy','R','Z','Te','Ti','ne',
     ,   (hlp_is(i),i=0,ns-1),(hlp_mn(i),i=0,nmoli-1),
     ,   (hlp_in(i),i=0,nioni-1),
     ,   'rad_tot',(hlp_al(i),i=1,nspecies),'rad_brm'
      else if (nmoli.gt.0 .and. nioni.eq.0) then
        write(hlp_frmr,'(a,i3,a,i2,a,i3,a)')
     ,   '(a9,7x,a1,11x,a1,10x,a2,10x,a2,10x,a2,7x,',
     .    ns,'(1x,a10,1x),',
     ,    nmoli,'(1x,a9,2x),a7,2x,',nspecies,'(4x,a6,2x),3x,a7)'
        write(lun,hlp_frmr) '#  ix  iy','R','Z','Te','Ti','ne',
     ,   (hlp_is(i),i=0,ns-1),(hlp_mn(i),i=0,nmoli-1),
     ,   'rad_tot',(hlp_al(i),i=1,nspecies),'rad_brm'
      else if (nmoli.eq.0 .and. nioni.gt.0) then
        write(hlp_frmr,'(a,i3,a,i2,a,i3,a)')
     ,   '(a9,7x,a1,11x,a1,10x,a2,10x,a2,10x,a2,7x,',
     .    ns,'(1x,a10,1x),',
     ,    nioni,'(1x,a9,2x),a7,2x,',nspecies,'(4x,a6,2x),3x,a7)'
        write(lun,hlp_frmr) '#  ix  iy','R','Z','Te','Ti','ne',
     ,   (hlp_is(i),i=0,ns-1),
     ,   (hlp_in(i),i=0,nioni-1),
     ,   'rad_tot',(hlp_al(i),i=1,nspecies),'rad_brm'
      else if (nmoli.eq.0 .and. nioni.eq.0) then
        write(hlp_frmr,'(a,i3,a,i3,a)')
     ,   '(a9,7x,a1,11x,a1,10x,a2,10x,a2,10x,a2,7x,',
     .    ns,'(1x,a10,1x),a7,2x,',nspecies,'(4x,a6,2x),3x,a7)'
        write(lun,hlp_frmr) '#  ix  iy','R','Z','Te','Ti','ne',
     ,   (hlp_is(i),i=0,ns-1),
     ,   'rad_tot',(hlp_al(i),i=1,nspecies),'rad_brm'
      end if
      write(hlp_frmr,'(a,i3,a)')
     , '(1x,1p,2i4,',7+ns+nmoli+nioni+nspecies,'e12.4)'

      do iy=0,ny-1 !{
        write (lun,*)
        do ix=0,nx-1 !{
          hlp_brm=0.0_R8
          hlp_rad=0.0_R8
          hlp_ni=0.0_R8
          do is=0,ns-1
            i=b2espcr(is)
            hlp_brm=hlp_brm+rqbrm(ix,iy,is)
            hlp_rad(i)=hlp_rad(i)+rqrad(ix,iy,is)
            if (.not.is_neutral(is).or.use_eirene.eq.0)
     >       hlp_ni(is)=na(ix,iy,is)
          end do !}
          do is=1,natmi !{
            i=latmscl(is)
            hlp_rad(i)=hlp_rad(i)-eneutrad(ix+1,iy+1,is,0)
            i=eb2atcr(is)
            hlp_ni(i)=hlp_ni(i)+dab2(ix+1,iy+1,is,1)
          end do !}
          do is=1,nmoli !{
            i=lmolscl(is)
            if (i.eq.0) cycle
            hlp_rad(i)=hlp_rad(i)-emolrad(ix+1,iy+1,is,0)
          end do !}
          do is=1,nioni !{
            i=lionscl(is)
            if (i.eq.0) cycle
            hlp_rad(i)=hlp_rad(i)-eionrad(ix+1,iy+1,is,0)
          end do !}
          hlp=hlp_brm
          do i=1,nspecies !{
            hlp=hlp+hlp_rad(i)
          end do !}
          if (nmoli.gt.0 .and. nioni.gt.0) then !{
            write (lun,hlp_frmr) ix,iy,
     ,       cr(ix,iy),cz(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy),
     ,      (hlp_ni(is),is=0,ns-1),(dmb2(ix+1,iy+1,is,1),is=1,nmoli),
     ,      (dib2(ix+1,iy+1,is,1),is=1,nioni),
     ,       1.0e-6_R8*hlp/vol(ix,iy),
     ,      (1.0e-6_R8*hlp_rad(i)/vol(ix,iy),i=1,nspecies),
     ,       1.0e-6_R8*hlp_brm/vol(ix,iy)
          else if (nmoli.gt.0 .and. nioni.eq.0) then !}{
            write (lun,hlp_frmr) ix,iy,
     ,       cr(ix,iy),cz(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy),
     ,      (hlp_ni(is),is=0,ns-1),(dmb2(ix+1,iy+1,is,1),is=1,nmoli),
     ,       1.0e-6_R8*hlp/vol(ix,iy),
     ,      (1.0e-6_R8*hlp_rad(i)/vol(ix,iy),i=1,nspecies),
     ,       1.0e-6_R8*hlp_brm/vol(ix,iy)
          else if (nmoli.eq.0 .and. nioni.gt.0) then !}{
            write (lun,hlp_frmr) ix,iy,
     ,       cr(ix,iy),cz(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy),
     ,      (hlp_ni(is),is=0,ns-1),
     ,      (dib2(ix+1,iy+1,is,1),is=1,nioni),
     ,       1.0e-6_R8*hlp/vol(ix,iy),
     ,      (1.0e-6_R8*hlp_rad(i)/vol(ix,iy),i=1,nspecies),
     ,       1.0e-6_R8*hlp_brm/vol(ix,iy)
          else if (nmoli.eq.0 .and. nioni.eq.0) then !}{
            write (lun,hlp_frmr) ix,iy,
     ,       cr(ix,iy),cz(ix,iy),te(ix,iy)/ev,ti(ix,iy)/ev,ne(ix,iy),
     ,      (hlp_ni(is),is=0,ns-1),
     ,       1.0e-6_R8*hlp/vol(ix,iy),
     ,      (1.0e-6_R8*hlp_rad(i)/vol(ix,iy),i=1,nspecies),
     ,       1.0e-6_R8*hlp_brm/vol(ix,iy)
          end if !}
        end do !}
      end do !}
      close (lun)
      deallocate(hlp_al,hlp_rad,hlp_ni,hlp_is)
      if (nmoli.gt.0) deallocate(hlp_mn)
      if (nioni.gt.0) deallocate(hlp_in)
      return
c=======================================================================
      end
