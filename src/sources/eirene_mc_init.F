      subroutine eirene_mc_init
      use b2mod_types
csw 17feb2011 extrab25
      use eirmod_extrab25
      implicit none
#     include <COUPLE/KOPPLDIM.F>
csw#     include <COUPLE/EIRDIAG.F>

#ifndef USE_EIRENE
      call eirene_extrab25_eirdiag_alloc(nxdd,nydd,
     .                     natm,nmol,nion,nlmpgs,nstra,npls)          !srv 14.10.13
      call eirene_extrab25_eirpbls_alloc(1,natm,nmol,nion,npls,       !srv 14.10.13
     &                                       nstra,nfl)               !srv 03.04.13
#endif
      return
      end
