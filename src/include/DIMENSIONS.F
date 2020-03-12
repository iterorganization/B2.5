!  Common dimensions
!
!  version : 01.12.98 21:42
!
#ifndef DIMENSIONS_DEFINED
#define DIMENSIONS_DEFINED

!
! parameters that are common to Eirene and B2
!

! NOTE: DEF_NXD should not include the additional cells to handle the cuts
#define DEF_NXD 128
#define DEF_NCUT 4
#define DEF_ISOEXTRA 0
#define DEF_NYD 36
#define DEF_NPHID 1
#define DEF_NSTRA (5*DEF_NCUT/2+2)*DEF_NATM+1
#define DEF_NFL DEF_NSD
#define DEF_NPLS (DEF_NSD-DEF_NATM)+(DEF_NATM+DEF_NMOL)*(DEF_NATM+DEF_NMOL)
#define DEF_NATM 6
#define DEF_NMOL 3
#define DEF_NION 3
#define DEF_NLIM 300
#define DEF_NSTS 50
#define DEF_NTOR 60
!*** Max. number of groups of Eirene surfaces for which the data can
!*** be transferred from B2 (DG specification "Surface special")
#define DEF_NSPCSRFG 40

! new! [2002.04.22]
#define DEF_NADD 4
#define DEF_NPRNL 50000
! new! [2002.06.14]
#define DEF_NSRFS 4

!
! parameters that are unique to B2
!

#define DEF_NSD 42
#define DEF_NBC 6*(DEF_NCUT/2+1)
#define DEF_NISO 10

!
! parameters that are unique to Eirene
!

#define DEF_NGSTAL 1

!
! parameters needed by uinp
!

#define DEF_NREAC 60
#define DEF_NREI 70
#define DEF_NREC DEF_NATM
#define DEF_NREL (DEF_NATM+DEF_NMOL)*(DEF_NATM+DEF_NMOL)

#endif
