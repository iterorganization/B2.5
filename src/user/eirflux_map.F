c
cvk interfaces to surface flux data calculated by Eirene

c  version: 11.07.2008 18:53

c
c calculating the index shift in the Eirene flux arrays.
c to be used in interface subroutines to Eirene flux data
c
c typ    : is the type of the surface which has to match
c          (for targets 2: radial surface)
c ind    : index of the radial surface which has to match
c istart,: starting and ending indices on this surface which have to match
c   iend :
c shift  : resulting shift which has to be used in accessing the arrays
c
      subroutine b2neut_ind(typ,ind,istart,iend,shift)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
      use b2mod_dimensions
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      integer,intent(in) :: typ,ind,istart,iend
      integer,intent(out) :: shift
      integer :: is

      do is=1,nnstsi
        if(typ.eq.eirdiag_nds_typ(is).and.
     f     ind.eq.eirdiag_nds_srf(is)) then
         if(istart.ge.eirdiag_nds_start(is).and.
     f      iend.le.eirdiag_nds_end(is)) then
c the surface has been found
           shift=eirdiag_nds_ind(is)+istart-eirdiag_nds_start(is)
           goto 20
         end if
        end if
      end do
      shift=-1
   20 return
      end subroutine b2neut_ind


c=======================================================================
c check if surface areas in B2 and Eirene are the same
c
c  s     : array of surface areas
c  t     : index of target
c  n     : size of array s
c  shift : shift which has to be used to access the Eirene array
c
      subroutine b2neut_surfcheck(s,t,n,shift)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_wneutrals
#else
      use b2mod_eirdiag
      use b2mod_dimensions
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      integer,intent(in) :: t,n,shift
      real (kind=R8),intent(in) :: s(n)
      integer :: i,ii
      real (kind=R8) :: ds,s2,sum1,sum2
      real (kind=R8),parameter :: eps=1e-3_R8
      logical :: iferr

      if(shift.lt.0) return

      iferr=.false.
      sum1=0.0_R8
      sum2=0.0_R8
      do i=1,n
        ii=i+shift
        s2=abs(s(i)+sarea_res(ii))
        ds=abs(s(i)-sarea_res(ii))
        sum1=sum1+s(i)
        sum2=sum2+sarea_res(ii)
        if(ds.gt.eps*s2) then
          iferr=.true.
        end if
      end do
      if(iferr) then
        write(*,*) "inconsistency of Eirene and B2 geometries ",
     w             "for target ", t
        write(*,*) " surface area in B2 ", sum1
        write(*,*) " surface area in Eirene ", sum2
      end if

      return
      end subroutine b2neut_surfcheck


c=======================================================================
cvk power flux delivered to the surface with neutrals
c
c  wneut : (kinetic energy of incident neutrals) -
c          (kinetic energy of reflected neutrals) +
c          (energy release from recombination of atoms into molecules)
c  wnpr  : kinetic energy of reflected neutrals originated from ions
c          (including sputtered neutrals)
c  n     : length of arrays wneut and wnpr
c  shift : data offset in the Eirene arrays



      subroutine b2neut_eng(wneut,wnpr,n,shift)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
      use b2mod_dimensions
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      integer,intent(in) :: n,shift
      real (kind=R8),intent(out) :: wneut(n),wnpr(n)
      integer :: j,i,ii

      if(shift.lt.0) return

      do j=1,n
       ii=j+shift
       wneut(j)=0.0_R8
       do i=1,nnatmi
        wneut(j)=wneut(j)+ewlda_res(i,ii)-ewldea_res(i,ii)
       end do
       do i=1,nnmoli
        wneut(j)=wneut(j)+ewldm_res(i,ii)-ewldem_res(i,ii)+
     +                    ewldmr_res(i,ii)
       end do
       wnpr(j)=ewldrp_res(ii)
      end do
      return
      end subroutine b2neut_eng

c=======================================================================
cvk incident neutral flux density
c
c  sfa   : incident flux density of atom
c  n     : length of array sfa wneut and wnpr
c  ia    : index of atom to be taken
c          if ia>natm then ia-natm is the index of molecule
c  shift : shift which has to be used to access the Eirene array
c
      subroutine b2neut_flux(sfa,n,ia,shift)
      use b2mod_types
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
      use b2mod_dimensions
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      integer,intent(in) :: n,ia,shift
      real (kind=R8),intent(out) :: sfa(n)
      integer :: j,ia0

      if(shift.lt.0) return
      if(ia.le.nnatmi) then
c atoms
       do j=1,n
        sfa(j)=wldna_res(ia,j+shift)
       end do
      else
c molecules
       ia0=ia-nnatmi
       do j=1,n
        sfa(j)=wldnm_res(ia0,j+shift)
       end do
      end if
      return
      end subroutine b2neut_flux



c=======================================================================
cak average energy of impinging neutrals
c
c  eav   : average energy of impinging neutrals [eV]
c  n     : length of array eav
c  ia    : index of atom to be taken
c          if ia>natm then ia-natm is the index of molecule
c  shift : shift which has to be used to access the Eirene array
c
      subroutine b2neut_enr(eav,n,ia,shift)
      use b2mod_types
      use b2mod_constants
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_extrab25
#else
      use b2mod_eirdiag
      use b2mod_dimensions
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      integer,intent(in) :: n,ia,shift
      real (kind=R8),intent(out) :: eav(n)
      integer :: j,ia0

      if(shift.lt.0) return
      if(ia.le.nnatmi) then
c atoms
       do j=1,n
        eav(j)=ewlda_res(ia,j+shift)/
     .     max(wldna_res(ia,j+shift),1.0_R8)/ev
       end do
      else
c molecules
       ia0=ia-nnatmi
       do j=1,n
        eav(j)=ewldm_res(ia0,j+shift)/
     .     max(wldnm_res(ia0,j+shift),1.0_R8)/ev
       end do
      end if
      return
      end subroutine b2neut_enr
