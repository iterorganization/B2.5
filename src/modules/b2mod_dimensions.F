      module b2mod_dimensions
      implicit none
!  Common dimensions
!
!  version : 01.12.98 21:42
!
#ifndef DIMENSIONS_DEFINED
#define DIMENSIONS_DEFINED

!
! parameters that are unique to B2
!

      integer, parameter :: DEF_NSD = 42
! DEF_NCUT is common to Eirene as well but needs to be defined here first
      integer, parameter :: DEF_NCUT = 4
      integer, parameter :: DEF_NBC = 6*(DEF_NCUT/2+1)
      integer, parameter :: DEF_NISO = 10

!
! parameters that are common to Eirene and B2
!

! NOTE: DEF_NXD should not include the additional cells to handle the cuts
      integer, parameter :: DEF_NXD = 128
      integer, parameter :: DEF_ISOEXTRA = 0
      integer, parameter :: DEF_NYD = 36
      integer, parameter :: DEF_NPHID = 1
      integer, parameter :: DEF_NATM = 6
      integer, parameter :: DEF_NMOL = 3
      integer, parameter :: DEF_NION = 3
      integer, parameter :: DEF_NSTRA = (5*DEF_NCUT/2+2)*DEF_NATM+1
      integer, parameter :: DEF_NFL = DEF_NSD
      integer, parameter :: DEF_NPLS = (DEF_NSD-DEF_NATM)+
     &                                 (DEF_NATM+DEF_NMOL)*
     &                                 (DEF_NATM+DEF_NMOL)
      integer, parameter :: DEF_NLIM = 300
      integer, parameter :: DEF_NSTS = 50
      integer, parameter :: DEF_NTOR = 60
!*** Max. number of groups of Eirene surfaces for which the data can
!*** be transferred from B2 (DG specification "Surface special")
      integer, parameter :: DEF_NSPCSRFG = 40

! new! [2002.04.22]
      integer, parameter :: DEF_NADD = 4
      integer, parameter :: DEF_NPRNL = 50000
! new! [2002.06.14]
      integer, parameter :: DEF_NSRFS = 4

!
! parameters that are unique to Eirene
!

      integer, parameter :: DEF_NGSTAL = 1

!
! parameters needed by uinp
!

      integer, parameter :: DEF_NREI = DEF_NATM+DEF_NMOL+DEF_NION+
     &                                 DEF_NFL
      integer, parameter :: DEF_NREC = DEF_NATM
      integer, parameter :: DEF_NRCX = DEF_NATM*(DEF_NION+DEF_NFL)
      integer, parameter :: DEF_NREL = (DEF_NATM+DEF_NMOL)*
     &                                 (DEF_NATM+DEF_NMOL)
      integer, parameter :: DEF_NREAC = DEF_NREI+DEF_NREC+DEF_NRCX+
     &                                  DEF_NREL

#endif
      end module b2mod_dimensions
