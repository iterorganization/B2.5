      subroutine b2pedges(lun)
c
c  version : 31.03.2005 17:35
c=======================================================================
c***  Produces a table of ion fluxes and plasma temperatures at the 
c***  outer edge of the grid for sputtering calculations
c***  for G. Strohmayer
c***
c*** Output files:
c***
c***  gridedge.dat with the data for G. Strohmayer
c=======================================================================
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_elements
      use b2mod_indirect
      use b2mod_constants
      implicit none
      integer lun

c*** Local variables
c***
c***  nx_hlp  : max. number of edge sections
c***  nl_hlp  : max. number of cells in an edge section
c***  ns_hlp  : number of edge sections (5 for SN, 10 for DN)
c***  ld_hlp  : direction of the fluxes in the order S W N E
c***  nc_hlp  : number of cells in a section
c***  ix, iy  : cell indices for T
c***  jx, jy  : cell indices for fluxes
c***  hlp_flx : fluxes of different species across the cell boundary
c***  mz_hlp  : ion charge
c***  hlp_al  : species name
c***  hlp_area: area of the corresponding cell face

      integer nx_hlp,nl_hlp
      parameter (nx_hlp=10)
      integer, allocatable :: mz_hlp(:)
      integer, allocatable :: ix(:,:), iy(:,:), jx(:,:), jy(:,:)
      integer i,j,k,l,m,n,i1,i2,il,is,ns_hlp
      integer ld_hlp(nx_hlp),nc_hlp(nx_hlp)
      real (kind=R8) ::  d,re,ze,ra,za,rb,zb,z,zz,hlp_area
      real (kind=R8), allocatable :: hlp_flx(:)
      character hlp_frmr*100, tab
      character*2, allocatable :: hlp_al(:)
c=======================================================================

      tab=char(9)
      nl_hlp=max(nx,ny)
      allocate(mz_hlp(0:ns-1),hlp_al(0:ns-1))
      allocate(hlp_flx(0:ns-1))
      allocate(ix(nl_hlp,nx_hlp),iy(nl_hlp,nx_hlp))
      allocate(jx(nl_hlp,nx_hlp),jy(nl_hlp,nx_hlp))
      do is=0,ns-1
        mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
        hlp_al(is)=is_codes(is)
      end do

      if(nncut.eq.2) then !{                        === single null ===

        ns_hlp=5

        m=1
        nc_hlp(m)=leftcut(1)
        ld_hlp(m)=1
        k=leftcut(1)-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

        m=2
        nc_hlp(m)=ny
        ld_hlp(m)=2
        k=0
        do i=1,nc_hlp(m) !{
          ix(i,m)=0
          iy(i,m)=k
          jx(i,m)=0
          jy(i,m)=k
          k=k+1
        end do !}

        m=3
        nc_hlp(m)=nx
        ld_hlp(m)=3
        k=0
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=ny-1
          jx(i,m)=k
          jy(i,m)=ny
          k=k+1
        end do !}

        m=4
        nc_hlp(m)=ny
        ld_hlp(m)=4
        k=ny-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=nx-1
          iy(i,m)=k
          jx(i,m)=nx
          jy(i,m)=k
          k=k-1
        end do !}

        m=5
        nc_hlp(m)=nx-rightcut(1)
        ld_hlp(m)=1
        k=nx-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

      else if(nncut.eq.4) then !}{ === double null ===
        
        ns_hlp=10

        m=1
        nc_hlp(m)=leftcut(1)
        ld_hlp(m)=1
        k=leftcut(1)-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

        m=2
        nc_hlp(m)=ny
        ld_hlp(m)=2
        k=0
        do i=1,nc_hlp(m) !{
          ix(i,m)=1
          iy(i,m)=k
          jx(i,m)=0
          jy(i,m)=k
          k=k+1
        end do !}

        m=3
        nc_hlp(m)=nxtl
        ld_hlp(m)=3
        k=1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=ny-1
          jx(i,m)=k
          jy(i,m)=ny
          k=k+1
        end do !}

        m=4
        nc_hlp(m)=ny
        ld_hlp(m)=4
        k=ny-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=nxtl-1
          iy(i,m)=k
          jx(i,m)=nxtl
          jy(i,m)=k
          k=k-1
        end do !}

        m=5
        nc_hlp(m)=nxtl-leftcut(2)-1
        ld_hlp(m)=1
        k=nxtl-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

        m=6
        nc_hlp(m)=rightcut(2)-nxtr-1
        ld_hlp(m)=1
        k=rightcut(2)-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

        m=7
        nc_hlp(m)=ny
        ld_hlp(m)=2
        k=0
        do i=1,nc_hlp(m) !{
          ix(i,m)=nxtr+1
          iy(i,m)=k
          jx(i,m)=nxtr+1
          jy(i,m)=k
          k=k+1
        end do !}

        m=8
        nc_hlp(m)=nx-nxtr-1
        ld_hlp(m)=3
        k=nxtr+1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=ny-1
          jx(i,m)=k
          jy(i,m)=ny
          k=k+1
        end do !}

        m=9
        nc_hlp(m)=ny
        ld_hlp(m)=4
        k=ny-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=nx-1
          iy(i,m)=k
          jx(i,m)=nx
          jy(i,m)=k
          k=k-1
        end do !}

        m=10
        nc_hlp(m)=nx-rightcut(1)
        ld_hlp(m)=1
        k=nx-1
        do i=1,nc_hlp(m) !{
          ix(i,m)=k
          iy(i,m)=0
          jx(i,m)=k
          jy(i,m)=0
          k=k-1
        end do !}

      else !}{
        write(0,*) 'b2pedges: wrong value of nncut (must be 2 or 4):',
     ,                nncut
        return
      end if !}
c!###{
c        print *
c        print *,'ns_hlp=',ns_hlp
c        if(ns_hlp.gt.0) then !{
c          print '(1x,a8,20i4/(9x,20i4))','ld_hlp: ',
c     ,                                          (ld_hlp(i),i=1,ns_hlp)
c          print '(1x,a8,20i4/(9x,20i4))','nc_hlp: ',
c     ,                                          (nc_hlp(i),i=1,ns_hlp)
c          k=0
c          print '(1x,5a4)','k','ix','iy','jx','jy'
c          do m=1,ns_hlp !{
c            print *,'                           m=',m,',  n=',nc_hlp(m)
c            if(nc_hlp(m).gt.0) then !{
c              do j=1,nc_hlp(m) !{
c                k=k+1
c                print '(1x,5i4)',k,ix(j,m),iy(j,m),jx(j,m),jy(j,m)
c              end do !}
c            end if !}
c          end do !}
c        end if !}
c        print *
c!###}

c*** Write the header

      open(lun,file='gridedge.dat')

      write (lun,*) '# Ion fluxes across the outer grid edge.',
     ,                                         '   Version : 05.05.2004'
      write (lun,*) '# ',label
      write (lun,*)
      write (lun,*) '# nx,ny,ns = ',nx,ny,ns
      if(nspecies.gt.0) write (lun,*) '# nspecies,nfluids = ',
     ,                                nspecies,(nfluids(i),i=1,nspecies)
      write (lun,*) '# nncut,leftcut,rightcut = ',
     ,                 nncut,(leftcut(i),i=1,nncut),
     ,                       (rightcut(i),i=1,nncut)

      write (lun,*)
      write (lun,*) '# il     : cell index in the list;'
      write (lun,*) '# ldir   : direction of the fluxes in the order ',
     ,                                                         'S W N E'
      write (lun,*) '# ix, iy : x and y indices in the grid;'
      write (lun,*) '# R , Z  : R and Z coordinates of the ',
     ,                                  'centres of the cell faces [m];'
      write (lun,*) '# x      : reference distance along the grid ',
     ,                                                       'edge [m];'
      write (lun,*) '# area   : area of the corresponding cell face ',
     ,                                                         '[m**2];'
      write (lun,*) '# Te, Ti : electron and ion temperatures [eV];'
      write (lun,*) '# flx... : fluxes of the ion species [1/m**2/s];'
      write (lun,*)

      write(hlp_frmr,710) ns
      print *,hlp_frmr               !###
  710 format('(15h#  il   ix   iy,7x,1hR,13x,1hZ,13x,1hx,11x,4harea,',
     ,        '11x,2hTi,10x,2hTe,5x,',i3,'(6h  flx_,i2.2,a2,2x))')
      write (lun,hlp_frmr) (mz_hlp(i),hlp_al(i),i=0,ns-1)

      write(hlp_frmr,720) ns+2
      print *,hlp_frmr               !###
  720 format('(1p,3i5,4e14.6,',i3,'e12.4)')
      k=0
      zz=0.0_R8
      do m=1,ns_hlp !{
        do l=1,nc_hlp(m) !{
          k=k+1
          do j=0,ns-1 !{
            if(ld_hlp(m).eq.1) then !{                            === S:
              hlp_area=gs(jx(l,m),jy(l,m),1)
              hlp_flx(j)=-fna(jx(l,m),jy(l,m),1,j)/hlp_area
            else if(ld_hlp(m).eq.2) then !}{                      === W:
              hlp_area=gs(jx(l,m),jy(l,m),0)
              hlp_flx(j)=-fna(jx(l,m),jy(l,m),0,j)/hlp_area
            else if(ld_hlp(m).eq.3) then !}{                      === N:
              hlp_area=gs(jx(l,m),jy(l,m),1)
              hlp_flx(j)=fna(jx(l,m),jy(l,m),1,j)/hlp_area
            else if(ld_hlp(m).eq.4) then !}{                      === E:
              hlp_area=gs(jx(l,m),jy(l,m),0)
              hlp_flx(j)=fna(jx(l,m),jy(l,m),0,j)/hlp_area
            else !}{
              write(0,*) 'b2pedges: internal error. m,l,k,ld_hlp(m)= ',
     ,                                                   m,l,k,ld_hlp(m)
              close (lun)
              return
            end if !}
            hlp_flx(j)=max(hlp_flx(j),1.0e6_R8)
          end do !}
          i1=((ld_hlp(m)-1)/2)*3
          i2=(mod(ld_hlp(m),4)/2)+1
          ra=crx(ix(l,m),iy(l,m),i1)
          za=cry(ix(l,m),iy(l,m),i1)
          rb=crx(ix(l,m),iy(l,m),i2)
          zb=cry(ix(l,m),iy(l,m),i2)
          d=sqrt((ra-rb)**2+(za-zb)**2)
          re=0.5_R8*(ra+rb)
          ze=0.5_R8*(za+zb)
          z=zz+0.5_R8*d
          zz=zz+d

          write(lun,hlp_frmr) k,ix(l,m),iy(l,m),re,ze,z,hlp_area,
     ,               ti(ix(l,m),iy(l,m))/ev,te(ix(l,m),iy(l,m))/ev,
     ,                                               (hlp_flx(j),j=1,ns)
        end do !}
      end do !}

      close (lun)
      deallocate(mz_hlp,hlp_al)
      deallocate(hlp_flx)
      deallocate(ix,iy,jx,jy)
      return
c=======================================================================
      end
