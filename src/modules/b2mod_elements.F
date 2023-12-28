      module b2mod_elements
      use b2mod_types
      implicit none

      character*2, save :: elements(92)
      data elements /
     & 'H' ,'He','Li','Be','B' ,'C' ,'N' ,'O' ,'F' ,'Ne',
     1 'Na','Mg','Al','Si','P' ,'S' ,'Cl','Ar','K' ,'Ca',
     2 'Sc','Ti','V' ,'Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3 'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y' ,'Zr',
     4 'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5 'Sb','Te','I' ,'Xe','Cs','Ba','La','Ce','Pr','Nd',
     6 'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7 'Lu','Hf','Ta','W' ,'Re','Os','Ir','Pt','Au','Hg',
     8 'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9 'Pa','U '/
      integer, save :: nuclear_mass(92)
      data nuclear_mass /
     &    1,   4,   7,   9,  11,  12,  14,  16,  18,  20,
     1   23,  24,  27,  28,  31,  32,  35,  40,  39,  40,
     2   45,  48,  51,  52,  55,  56,  59,  58,  63,  64,
     3   70,  73,  75,  80,  80,  84,  85,  88,  89,  90,
     4   93,  96,  98, 100, 103, 106, 108, 112, 115, 118,
     5  122, 128, 127, 132, 133, 138, 139, 140, 141, 144,
     6  145, 152, 152, 156, 159, 162, 165, 167, 169, 174,
     7  175, 180, 181, 184, 187, 192, 193, 195, 197, 200,
     8  205, 208, 209, 209, 210, 222, 223, 226, 227, 232,
     9  231, 238/

#ifndef PSR
      character*2, allocatable, save ::
     *  is_codes(:)
#else
      character*2, save ::
     *  is_codes(0:100)
#endif

      contains

      subroutine alloc_b2mod_elements(nsd,zn,am)
      implicit none
      integer nsd
      real (kind=R8) ::
     *   zn(0:nsd-1)
      real (kind=R8) ::
     *   am(0:nsd-1)
      integer is, iz, im

#ifndef PSR
      allocate(is_codes(0:nsd-1))
#else
      call xertst(nsd.le.100, 
     *  'alloc_b2mod_elements: increase size of is_codes')
#endif
      
      do is=0,nsd-1
        iz=nint(zn(is))
        if(iz.le.1) then
          im=nint(am(is))
          if(im.eq.1) is_codes(is)='H'
          if(im.eq.2) is_codes(is)='D'
          if(im.eq.3) is_codes(is)='T'
        else
          is_codes(is)=elements(iz)
        end if
        write(*,*) is,zn(is),am(is),trim(is_codes(is))
      end do
     
      return
      end subroutine alloc_b2mod_elements

      subroutine dealloc_b2mod_elements
      implicit none

#ifndef PSR
      deallocate(is_codes)
#endif
      
      return
      end subroutine dealloc_b2mod_elements

      integer function get_atomic_number(name)
      implicit none
      character*(*) name
      character*8 species_name
      integer i
      logical skip, openp
      integer lnblnk
      logical streql, isadigit, isplusorminus
      external streql, lnblnk, isadigit, isplusorminus, strip_spaces

      get_atomic_number = 0
      species_name = repeat(' ',8)
      openp = .false.
      do i = 1, lnblnk(name)
        if (.not.openp) then
          skip = .false.
          skip = skip .or. isadigit(name(i:i))
          skip = skip .or. isplusorminus(name(i:i))
          skip = skip .or. name(i:i).eq.'*'
          openp = name(i:i).eq.'('
        else
          openp = name(i:i).ne.')'
          skip = .true.
        endif
        skip = skip .or. openp
        if (.not.skip) species_name(i:i) = name(i:i)
      end do
      call strip_spaces(species_name)
      do i = 1, 92
        if(streql(species_name(1:2),elements(i))) then
          get_atomic_number = i
          return
        end if
      end do
      if (streql(species_name,'D').or.
     &    streql(species_name,'T')) get_atomic_number=1

      return
      end function get_atomic_number

      end module b2mod_elements
