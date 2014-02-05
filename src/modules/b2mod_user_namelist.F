      module b2mod_user_namelist
      use b2mod_types
#     include "DIMENSIONS.F"
      integer nhetrgtx, npfrgrpd
      parameter (nhetrgtx=DEF_NLIM, npfrgrpd=6)
      real (kind=R8) ::
     ,  fusion_power,Spmp_He_to_D,Spmp_nom
      integer lhetrgts(nhetrgtx),lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,
     ,  j_he_at,j_h_at,l_h_mol(DEF_NMOL),lpfrs_pmp,
     ,  npfrgrp,lpfrgrp(DEF_NLIM),ipfrgrp(npfrgrpd),jpfrgrp(npfrgrpd)
      character gpfrgrp(npfrgrpd)*8
      character*80 user_filename
      logical user_namelist_used

c----------------------------------------------------------------------
c*** npfrgrpd     : max. number of surface groups for PFR flows
c*** fusion_power : fusion power (including neutrons) [MW]
c*** Spmp_He_to_D : ratio of He to DT pumping speeds (typically, 0.8)
c*** lhetrgts  :  list of surface indices (EIRENE notation) which
c***              are used for calculation of helium enrichment
c*** lpfrb_[io]:  B2 x-cell indices corresponding to the bottom of the
c***              bypass between the inner and outer divertor (counted
c***              from the bottom, if non-positive)
c*** lpfrt_[io]:  B2 x-cell indices corresponding to the top of the
c***              bypass between the inner and outer divertor (counted
c***              from the x-point, if non-positive)
c*** j_he_at   :  Helium index in EIRENE atoms list
c*** j_h_at    :  Hydrogen index in EIRENE atoms list
c*** l_h_mol   :  Hydrogen content in the molecules
c*** lpfrs_pmp :  Location of the pump.
c***              0 no pump at all, 1 - lower PFR, 2 - lower outer
c*** npfrgrp   :  actual number of surface groups for PFR flows
c*** lpfrgrp   :  list of surface segments for PFR flows
c*** [ij]pfrgrp:  first and last positions in this list for each group
c----------------------------------------------------------------------

      contains

      subroutine read_b2mod_user_namelist(nx)
      use b2mod_b2cmpa
      implicit none
      integer nx

      namelist /user/ 
     ,  lhetrgts,lpfrb_i,lpfrb_o,lpfrt_i,lpfrt_o,j_he_at,j_h_at,
     ,  l_h_mol,fusion_power,Spmp_He_to_D,lpfrs_pmp,
     ,  npfrgrp,lpfrgrp,ipfrgrp,jpfrgrp,gpfrgrp,
     ,  Spmp_nom,
     ,  user_filename
      logical :: file_ok
      character*260 filename
      integer i
      external subini, subend, xertst

      call subini ('readuser')

*  -- USER --
c*** PFR specification - for He enrichment and neutral pressure
c***  Standard setup:
      lpfrs_pmp=0       ! no pumping
      fusion_power=0    ! no helium
      Spmp_He_to_D=1    ! same pumping speed for He and D
      Spmp_nom=0
      j_h_at = 0
      l_h_mol=0
      j_he_at=0
      lpfrb_i=-2
      lpfrb_o=-2
      lpfrt_i=-2
      lpfrt_o=-2
      npfrgrp=0
      ipfrgrp=0
      jpfrgrp=0
      do i=1,nhetrgtx
        lhetrgts(i)=0
      end do
      lhetrgts(1)=-2
      lhetrgts(2)=-3
      gpfrgrp=' '
      user_filename='b2.user.parameters'

      filename='../' // user_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist user from ',trim(filename)
        read(99,user)
        write(*,*) 'Read namelist user from ',trim(filename)
        user_namelist_used = .true.
      else
        write(*,*) 'No ',trim(filename)
        user_namelist_used = .false.
      endif

*   ..test input
      call xertst(0.le.fusion_power,
     . 'faulty input fusion_power')
      call xertst(0.le.Spmp_nom,
     . 'faulty input Spmp_nom')
      call xertst(0.le.Spmp_He_to_D,
     . 'faulty input Spmp_He_to_D')
      if (lpfrb_i.ne.-2) call xertst (0.le.lpfrb_i.and.lpfrb_i.lt.nx,
     . 'faulty input lpfrb_i')
      if (lpfrb_o.ne.-2) call xertst (0.le.lpfrb_o.and.lpfrb_o.lt.nx,
     . 'faulty input lpfrb_o')
      if (lpfrt_i.ne.-2) call xertst (0.le.lpfrt_i.and.lpfrt_i.lt.nx,
     . 'faulty input lpfrt_i')
      if (lpfrt_o.ne.-2) call xertst (0.le.lpfrt_o.and.lpfrt_o.lt.nx,
     . 'faulty input lpfrt_o')
      call xertst (0.le.npfrgrp,
     . 'faulty input npfrgrp')
      do i = 1, npfrgrp
        call xertst(1.le.ipfrgrp(i).and.ipfrgrp(i).le.jpfrgrp(i),
     .   'faulty inputs ipfrgrp, jpfrgrp')
      end do
      do i = 2, npfrgrp
        call xertst(jpfrgrp(i-1).lt.ipfrgrp(i),
     .   'faulty inputs ipfrgrp, jpfrgrp')
      end do

      call subend ()

      entry write_b2mod_user_namelist

      if (user_namelist_used) write(*,user)
      return

      end subroutine read_b2mod_user_namelist

      end module b2mod_user_namelist
