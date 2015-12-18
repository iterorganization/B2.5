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
        write(*,*) is,zn(is),am(is),is_codes(is)
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

      end module b2mod_elements
