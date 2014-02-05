      subroutine eirene_mc_init
      use b2mod_types
csw 17feb2011 extrab25
      use eirmod_extrab25
      implicit none
#     include <COUPLE/KOPPLDIM.F>
csw#     include <COUPLE/EIRDIAG.F>

c      call eirene_extrab25_dummyalloc(nxdd,nydd,
c     .                     natm,nmol,nion,nlmpgs,nstra)
c      call eirene_extrab25_dummyeirpbls_init(1,natm,nmol,nion,npls,
c     &                                       nstra,nfl)               !srv 03.04.13
      return
      end
