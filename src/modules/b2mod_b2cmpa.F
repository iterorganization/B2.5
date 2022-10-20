      module b2mod_b2cmpa
      use b2mod_types
#include "DIMENSIONS.F"
*     (/b2cmpa/ contains basic physics parameters)
* DPC: constant declaration now in b2mod_constants
      integer nsdecl
      parameter (nsdecl=DEF_NSD)
      integer nspecies
      integer partition_number(0:nsdecl-1)
      real (kind=R8), save ::
     *  zn(0:nsdecl-1), am(0:nsdecl-1),
     &  zamax(0:nsdecl-1), zamin(0:nsdecl-1),
     &  potmin(0:nsdecl-1), potmax(0:nsdecl-1),
     &  potimin(0:nsdecl-1), potimax(0:nsdecl-1)
      logical is_neutral(0:nsdecl-1)
      real (kind=R8) :: amtol
      parameter (amtol=1.0e-2_R8)

c This module contains definitions for bundled stages
c     zn      : Nuclear charge of bundled stage (is)
c     am      : Atomic mass of bundled stage (is) (in amu)
c     zamax   : Highest ionisation stage of bundle (is)
c     zamin   : Lowest ionisation stage of bundle (is)
c     potmin  : Cumulative ionisation cost to reach zamin(is) (in eV)
c     potmax  : Cumulative ionisation cost to reach zamax(is) (in eV)
c     potimin : Ionisation cost to reach zamin(is) from zamax(is-1) (in eV)
c     potimax : Ionisation cost to reach zamax(is) from zamin(is-1) (in eV)
c     is_neutral(is) is true if a stage is a neutral species
c     partition_number indicates the stage number within its isonuclear sequence
c     
c     It is not allowed to bundle neutral and charged states together.

      contains

      subroutine initialise_potentials(ns)
      use b2mod_constants
      implicit none
      logical, save :: initialised=.false.
      real (kind=R8), allocatable :: pdum(:)
      logical open_file
      integer ns,is,iss,js
      external open_file
      intrinsic nint

      if(.not.initialised) then
        if(open_file(99,'ionization_potentials',
     1    'data.local|modules/B2.5/Database')) then
          is=0
          do while(is .lt. ns)
            if(.not.is_neutral(is)) then
              rewind(99)
              allocate(pdum(1:nint(zn(is))))
              do iss=1,nint(zn(is))-1
                read(99,*)
              end do
              read(99,*) (pdum(iss),iss=1,nint(zn(is)))
              if(nint(zn(is)) .eq. 1 .and. .not.is_neutral(is)) then
                pdum(1)=pdum(1)*(am(is)*(me+mp))/(me+am(is)*mp)
              endif
              potimin(is)=0.0_R8
              potmin(is)=0.0_R8
              potimax(is)=0.0_R8
              potmax(is)=0.0_R8
              do iss = 1, nint(zamin(is))
                potimin(is) = potimin(is)+pdum(iss)
                potmin(is)=potmin(is)+pdum(iss)
                potimax(is) = potimax(is)+pdum(iss)
                potmax(is)=potmax(is)+pdum(iss)
              enddo
              do iss = nint(zamin(is))+1, nint(zamax(is))
                potimax(is) = potimax(is)+pdum(iss)
                potmax(is)=potmax(is)+pdum(iss)
              enddo
              iss = is + 1
              do while (lnext(is,min(iss,ns-1)).and.iss.lt.ns)
                potimin(iss)=0.0_R8
                potimax(iss)=0.0_R8
                do js = nint(zamin(iss-1)), nint(zamin(iss))
                  potmin(iss)=potmin(iss-1)+pdum(js)
                  potmax(iss)=potmax(iss-1)+pdum(js)
                enddo
                do js = nint(zamin(iss)+1), nint(zamax(iss))
                  potmax(iss)=potmax(iss)+pdum(js)
                enddo
                do js = nint(zamax(iss-1)+1), nint(zamin(iss))
                  potimin(iss)=potimin(iss)+pdum(js)
                enddo
                do js = nint(zamin(iss-1)+1), nint(zamax(iss))
                  potimax(iss)=potimax(iss)+pdum(js)
                enddo
                iss = iss + 1
              enddo
              is=is+nint(zn(is))
              deallocate(pdum)
            else
              potimin(is)=0.0_R8
              potmin(is)=potimin(is)
              potimax(is)=0.0_R8
              potmax(is)=potimax(is)
              is=is+1
            end if
          end do
          close(99)
        else
          do is=0,ns-1
            potimin(is)=0.0_R8
            potmin(is)=potimin(is)
            potimax(is)=0.0_R8
            potmax(is)=potimax(is)
            if(nint(zn(is)) .eq. 1 .and. .not.is_neutral(is)) then
              potimin(is)=13.598434_R8*(am(is)*(me+mp))/(me+am(is)*mp)
              potmin(is)=potimin(is)
              potimax(is)=potimin(is)
              potmax(is)=potimax(is)
            endif
          end do
        end if
        do is=0,ns-1
          write(*,'(i3,1x,f5.1,1x,f5.1,1x,f5.1,4f11.3)')
     1     is,zamin(is),zamax(is),zn(is),
     2     potimin(is),potimax(is),potmin(is),potmax(is)
        end do
        initialised=.true.
      endif

      return
      end subroutine initialise_potentials

      logical function lnext(is,js)
      implicit none
      integer, intent(in) :: is, js
      lnext = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.amtol.and.is.ne.js
      return
      end function lnext

      end module b2mod_b2cmpa
