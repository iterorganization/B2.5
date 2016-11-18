      module b2mod_plasma
      use b2mod_types
      use b2mod_version
      implicit none

      real (kind=R8), allocatable, save ::
     *  fch(:,:,:),  fch0(:,:,:),  fchp(:,:,:),  fch_32(:,:,:),
     *  fch_52(:,:,:),
     *  fna_32(:,:,:,:),  fna_52(:,:,:,:),
     *  fni_32(:,:,:),  fni_52(:,:,:),
     *  fne_32(:,:,:),  fne_52(:,:,:),
     *  fhe(:,:,:),  fhe0(:,:,:),  fhep(:,:,:),  fhet(:,:,:),
     *  fhi(:,:,:),  fhi0(:,:,:),  fhip(:,:,:),  fhit(:,:,:),
     *  fna(:,:,:,:),  fna0(:,:,:,:),  fnap(:,:,:,:),
     *  fne(:,:,:),  fni(:,:,:),
     *  na(:,:,:),  na0(:,:,:),  nap(:,:,:),
     *  ne(:,:),  ne0(:,:),  ne2(:,:),  nep(:,:),
     *  ni(:,:,:),  ni0(:,:,:),
     *  pb(:,:),
     *  po(:,:),  po0(:,:),  pop(:,:),
     *  te(:,:),  te0(:,:),  tep(:,:),
     *  ti(:,:),  ti0(:,:),  tip(:,:),
     *  ua(:,:,:),  ua0(:,:,:),  uap(:,:,:),
     *  uadia(:,:,:,:),  fchdia(:,:,:),
     *  fch_p(:,:,:), uahel(:,:,:),                                      !srv 23.12.03 10.06.08 {
     *  fhm(:,:,:,:), fhp(:,:,:,:), fhj(:,:,:), fht(:,:,:),
     *  fhe_mdf(:,:,:), fhePSch(:,:,:), fhi_mdf(:,:,:), fhiPSch(:,:,:),  !srv 13.10.06
     *  fna_mdf(:,:,:,:), fnaPSch(:,:,:,:), 
     *  fna_fcor(:,:,:,:), fna_he(:,:,:,:),                              !srv 01.10.99 21.01.00
     *  fna_eir(:,:,:,:), fne_eir(:,:,:), fhe_eir(:,:,:),                !srv 14.07.10
     *  fhi_eir(:,:,:), 
     *  fna_32nd(:,:,:,:), fna_52nd(:,:,:,:)                             !srv 18.06.08

      real (kind=R8), allocatable, save ::
     *  fchvispar(:,:,:), fchvisper(:,:,:), fchvisq(:,:,:),              !srv 11.07.99
     *  fchinert(:,:,:), fchin(:,:,:), fchanml(:,:,:),                   !srv 07.07.99 24.03.07 10.06.08 }
     *  vadia(:,:,:,:), vedia(:,:,:), wadia(:,:,:,:),                    !xpb
     *  vaecrb(:,:,:,:), veecrb(:,:,:),
     *  wedia(:,:,:), ve(:,:),                                           !srv 14.07.10
     *  fna_nodrift(:,:,:,:), 
     *  kinrgy(:,:,:), kinrgy0(:,:,:), kinrgyp(:,:,:),                   !xpb
     *  fmo(:,:,:,:)
      real (kind=R8), allocatable, save ::
     *  facdrift(:,:), fac_ExB(:,:), fac_vis(:,:)

      real (kind=R8), allocatable, save ::                               !djm 17/10/16
     *  fna_pinch(:,:,:,:), fna_pll(:,:,:,:), fna_drift(:,:,:,:),
     *  fna_ch(:,:,:,:), fna_nanom(:,:,:,:), fna_panom(:,:,:,:)

      contains

      subroutine alloc_b2mod_plasma(nxd,nyd,nsd)
      implicit none
      integer nxd,nyd,nsd

      allocate(fch(-1:nxd,-1:nyd,0:1))
      allocate(fch0(-1:nxd,-1:nyd,0:1))
      allocate(fchp(-1:nxd,-1:nyd,0:1))
      allocate(fch_32(-1:nxd,-1:nyd,0:1))
      allocate(fch_52(-1:nxd,-1:nyd,0:1))
      allocate(fna_32(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_52(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fni_32(-1:nxd,-1:nyd,0:1))
      allocate(fni_52(-1:nxd,-1:nyd,0:1))
      allocate(fne_32(-1:nxd,-1:nyd,0:1))
      allocate(fne_52(-1:nxd,-1:nyd,0:1))
      allocate(fna_32nd(-1:nxd,-1:nyd,0:1,0:nsd-1))                !srv 18.06.08
      allocate(fna_52nd(-1:nxd,-1:nyd,0:1,0:nsd-1))                !srv 18.06.08
      allocate(fna_pinch(-1:nxd,-1:nyd,0:1,0:nsd-1))               !djm 17/10/16
      allocate(fna_pll(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_drift(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_ch(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_nanom(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna_panom(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fhe(-1:nxd,-1:nyd,0:1))
      allocate(fhe0(-1:nxd,-1:nyd,0:1))
      allocate(fhep(-1:nxd,-1:nyd,0:1))
      allocate(fhet(-1:nxd,-1:nyd,0:1))
      allocate(fhe_mdf(-1:nxd,-1:nyd,0:1))                         !srv 13.10.06
      allocate(fhePSch(-1:nxd,-1:nyd,0:1))                         !srv 13.10.06
      allocate(fhi(-1:nxd,-1:nyd,0:1))
      allocate(fhi0(-1:nxd,-1:nyd,0:1))
      allocate(fhip(-1:nxd,-1:nyd,0:1))
      allocate(fhit(-1:nxd,-1:nyd,0:1))
      allocate(fhi_mdf(-1:nxd,-1:nyd,0:1))                         !srv 13.10.06
      allocate(fhiPSch(-1:nxd,-1:nyd,0:1))                         !srv 13.10.06
      allocate(fna_he(-1:nxd,-1:nyd,0:1,0:nsd-1))                  !srv 21.01.00
      allocate(fna_fcor(-1:nxd,-1:nyd,0:1,0:nsd-1))                !srv 01.10.99
      allocate(fnaPSch(-1:nxd,-1:nyd,0:1,0:nsd-1))                 !srv 13.10.06
      allocate(fna_mdf(-1:nxd,-1:nyd,0:1,0:nsd-1))                 !srv 13.10.06
      allocate(fna_eir(-1:nxd,-1:nyd,0:1,0:nsd-1))                 !srv 14.07.10 {
      allocate(fne_eir(-1:nxd,-1:nyd,0:1))
      allocate(fhe_eir(-1:nxd,-1:nyd,0:1))
      allocate(fhi_eir(-1:nxd,-1:nyd,0:1))
      fna_he = 0.0_R8
      fnaPSch = 0.0_R8
      fhePSch = 0.0_R8
      fhiPSch = 0.0_R8
      fna_mdf = 0.0_R8
      fhe_mdf = 0.0_R8
      fhi_mdf = 0.0_R8
      fna_eir = 0.0_R8
      fne_eir = 0.0_R8
      fhe_eir = 0.0_R8
      fhi_eir = 0.0_R8  
      fna_fcor = 0.0_R8                                            !srv 14.07.10 }
      allocate(fna(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fna0(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fnap(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fne(-1:nxd,-1:nyd,0:1))
      allocate(fni(-1:nxd,-1:nyd,0:1))
      allocate(fhm(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fhp(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(fhj(-1:nxd,-1:nyd,0:1))
      allocate(fht(-1:nxd,-1:nyd,0:1))
      allocate(na(-1:nxd,-1:nyd,0:nsd-1))
      allocate(na0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(nap(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ne(-1:nxd,-1:nyd))
      allocate(ne0(-1:nxd,-1:nyd))
      allocate(ne2(-1:nxd,-1:nyd))
      allocate(nep(-1:nxd,-1:nyd))
      allocate(ni(-1:nxd,-1:nyd,0:1))
      allocate(ni0(-1:nxd,-1:nyd,0:1))
      allocate(pb(-1:nxd,-1:nyd))
      allocate(po(-1:nxd,-1:nyd))
      allocate(po0(-1:nxd,-1:nyd))
      allocate(pop(-1:nxd,-1:nyd))
      allocate(te(-1:nxd,-1:nyd))
      allocate(te0(-1:nxd,-1:nyd))
      allocate(tep(-1:nxd,-1:nyd))
      allocate(ti(-1:nxd,-1:nyd))
      allocate(ti0(-1:nxd,-1:nyd))
      allocate(tip(-1:nxd,-1:nyd))
      allocate(ua(-1:nxd,-1:nyd,0:nsd-1))
      allocate(ua0(-1:nxd,-1:nyd,0:nsd-1))
      allocate(uap(-1:nxd,-1:nyd,0:nsd-1))
      allocate(uadia(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(uahel(-1:nxd,-1:nyd,0:nsd-1))                         !srv 23.12.03
      allocate(fchdia(-1:nxd,-1:nyd,0:1))
      fchdia=0.0_R8
      allocate(fch_p(-1:nxd,-1:nyd,0:1))                             !xpb srv 10.06.08
      allocate(fchvispar(-1:nxd,-1:nyd,0:1))                         !xpb
      fchvispar=0.0_R8
      allocate(fchvisper(-1:nxd,-1:nyd,0:1))                         !xpb
      fchvisper=0.0_R8
      allocate(fchvisq(-1:nxd,-1:nyd,0:1))                           !srv 11.07.99
      fchvisq=0.0_R8
      allocate(fchinert(-1:nxd,-1:nyd,0:1))                          !srv 07.07.99
      fchinert=0.0_R8
      allocate(fchin(-1:nxd,-1:nyd,0:1))                             !srv 30.09.02
      fchin=0.0_R8
      allocate(fchanml(-1:nxd,-1:nyd,0:1))                           !srv 24.03.07
      allocate(vadia(-1:nxd,-1:nyd,0:1,0:nsd-1))                     !xpb
      vadia=0.0_R8
      allocate(vedia(-1:nxd,-1:nyd,0:1))                             !xpb
      allocate(wadia(-1:nxd,-1:nyd,0:1,0:nsd-1))                     !xpb
      allocate(wedia(-1:nxd,-1:nyd,0:1))                             !srv 14.07.10
      allocate(ve(-1:nxd,-1:nyd))                                    !srv 14.07.10
      allocate(vaecrb(-1:nxd,-1:nyd,0:1,0:nsd-1))                    !xpb
      vaecrb=0.0_R8
      allocate(veecrb(-1:nxd,-1:nyd,0:1))                            !xpb
      allocate(fna_nodrift(-1:nxd,-1:nyd,0:1,0:nsd-1))               !xpb
      fna_nodrift=0.0_R8
      allocate(kinrgy(-1:nxd,-1:nyd,0:nsd-1))                        !xpb
      allocate(kinrgy0(-1:nxd,-1:nyd,0:nsd-1))                       !xpb
      allocate(kinrgyp(-1:nxd,-1:nyd,0:nsd-1))                       !xpb
      allocate(fmo(-1:nxd,-1:nyd,0:1,0:nsd-1))
      allocate(facdrift(-1:nxd,-1:nyd))                              !dpc
      facdrift=0.0_R8
      allocate(fac_ExB(-1:nxd,-1:nyd))                               !dpc
      fac_ExB=0.0_R8
      allocate(fac_vis(-1:nxd,-1:nyd))                               !dpc
      fac_vis=0.0_R8

      return
      end subroutine alloc_b2mod_plasma

      subroutine dealloc_b2mod_plasma
      implicit none

      deallocate(fch)
      deallocate(fch0)
      deallocate(fchp)
      deallocate(fch_32)
      deallocate(fch_52)
      deallocate(fna_32)
      deallocate(fna_52)
      deallocate(fna_32nd)                                         !srv 18.06.08
      deallocate(fna_52nd)                                         !srv 18.06.08
      deallocate(fni_32)
      deallocate(fni_52)
      deallocate(fne_32)
      deallocate(fne_52)
      deallocate(fna_pinch) !djm 17/10/16
      deallocate(fna_pll)
      deallocate(fna_drift)
      deallocate(fna_ch)
      deallocate(fna_nanom)
      deallocate(fna_panom)
      deallocate(fhe)
      deallocate(fhe0)
      deallocate(fhep)
      deallocate(fhet)
      deallocate(fhi)
      deallocate(fhi0)
      deallocate(fhip)
      deallocate(fhit)
      deallocate(fna)
      deallocate(fna0)
      deallocate(fnap)
      deallocate(fne)
      deallocate(fni)
      deallocate(fhm)
      deallocate(fhp)
      deallocate(fhj)
      deallocate(fht)
      deallocate(na)
      deallocate(na0)
      deallocate(nap)
      deallocate(ne)
      deallocate(ne0)
      deallocate(ne2)
      deallocate(nep)
      deallocate(ni)
      deallocate(ni0)
      deallocate(pb)
      deallocate(po)
      deallocate(po0)
      deallocate(pop)
      deallocate(te)
      deallocate(te0)
      deallocate(tep)
      deallocate(ti)
      deallocate(ti0)
      deallocate(tip)
      deallocate(ua)
      deallocate(ua0)
      deallocate(uap)
      deallocate(uadia)
      deallocate(fchdia)
      deallocate(fch_p)
      deallocate(fchvispar)
      deallocate(fchvisper)
      deallocate(fchvisq)                           !srv 11.07.99
      deallocate(fchinert)                          !srv 07.07.99
      deallocate(fchin)                             !srv 30.09.02
      deallocate(fchanml)                           !srv 24.03.07
      deallocate(vadia)
      deallocate(vedia)
      deallocate(wadia)
      deallocate(vaecrb)
      deallocate(veecrb)
      deallocate(fna_nodrift)
      deallocate(kinrgy)
      deallocate(kinrgy0)
      deallocate(kinrgyp)
      deallocate(fmo)
      deallocate(facdrift)
      deallocate(fac_ExB)
      deallocate(fac_vis)
      deallocate(fhe_mdf)                           !srv 13.10.06
      deallocate(fhePSch)                           !srv 13.10.06
      deallocate(fhi_mdf)                           !srv 13.10.06
      deallocate(fhiPSch)                           !srv 13.10.06
      deallocate(fna_he)                            !srv 21.01.00
      deallocate(fna_fcor)                          !srv 01.10.99
      deallocate(fnaPSch)                           !srv 13.10.06
      deallocate(fna_mdf)                           !srv 13.10.06
      deallocate(uahel  )                           !srv 23.12.03
      deallocate(fna_eir)                           !srv 14.07.10 {
      deallocate(fne_eir)
      deallocate(fhe_eir)
      deallocate(fhi_eir)
      deallocate(wedia)
      deallocate(ve)                                !srv 14.07.10 }

      return
      end subroutine dealloc_b2mod_plasma

      subroutine read_b2mod_plasma(nxd,nyd,nsd,nread)
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer nxd,nyd,nsd,nread
      integer is

      call cfrure (nread, size(fch), fch, 'fch')
      call cfrure (nread, size(fch0), fch0, 'fch0')
      call cfrure (nread, size(fchp), fchp, 'fchp')
      call cfrure (nread, size(fhe), fhe, 'fhe')
      call cfrure (nread, size(fhe0), fhe0, 'fhe0')
      call cfrure (nread, size(fhep), fhep, 'fhep')
      if(b2fplasma_version.gt.'01.001.003') then
        call cfrure (nread, size(fhet), fhet, 'fhet')
      else
        fhet = fhe
      endif
      call cfrure (nread, size(fhi), fhi, 'fhi')
      call cfrure (nread, size(fhi0), fhi0, 'fhi0')
      call cfrure (nread, size(fhip), fhip, 'fhip')
      if(b2fplasma_version.gt.'01.001.003') then
        call cfrure (nread, size(fhit), fhit, 'fhit')
      else
        fhit = fhi
      endif
      call cfrure (nread, size(fna), fna, 'fna')
      call cfrure (nread, size(fna0), fna0, 'fna0')
      call cfrure (nread, size(fnap), fnap, 'fnap')
      call cfrure (nread, size(fne), fne, 'fne')
      call cfrure (nread, size(fni), fni, 'fni')
      call cfrure (nread, size(na), na, 'na')
      call cfrure (nread, size(na0), na0, 'na0')
      call cfrure (nread, size(nap), nap, 'nap')
      call cfrure (nread, size(ne), ne, 'ne')
      call cfrure (nread, size(ne0), ne0, 'ne0')
      call cfrure (nread, size(ne2), ne2, 'ne2')
      call cfrure (nread, size(nep), nep, 'nep')
      call cfrure (nread, size(ni), ni, 'ni')
      call cfrure (nread, size(ni0), ni0, 'ni0')
      call cfrure (nread, size(pb), pb, 'pb')
      call cfrure (nread, size(po), po, 'po')
      call cfrure (nread, size(po0), po0, 'po0')
      call cfrure (nread, size(pop), pop, 'pop')
      call cfrure (nread, size(te), te, 'te')
      call cfrure (nread, size(te0), te0, 'te0')
      call cfrure (nread, size(tep), tep, 'tep')
      call cfrure (nread, size(ti), ti, 'ti')
      call cfrure (nread, size(ti0), ti0, 'ti0')
      call cfrure (nread, size(tip), tip, 'tip')
      call cfrure (nread, size(ua), ua, 'ua')
      call cfrure (nread, size(ua0), ua0, 'ua0')
      call cfrure (nread, size(uap), uap, 'uap')
      call cfrure (nread, size(uadia), uadia, 'uadia')
      call cfrure (nread, size(fchdia), fchdia, 'fchdia')
      if(b2fplasma_version.ge.'01.001.008') then
        call cfrure (nread, size(fmo), fmo, 'fmo')
        call cfrure (nread, size(fna_32), fna_32, 'fna_32')
        call cfrure (nread, size(fna_52), fna_52, 'fna_52')
        call cfrure (nread, size(fni_32), fni_32, 'fni_32')
        call cfrure (nread, size(fni_52), fni_52, 'fni_52')
        call cfrure (nread, size(fne_32), fne_32, 'fne_32')
        call cfrure (nread, size(fne_52), fne_52, 'fne_52')
      else
        do is = 0, nsd-1
          fmo(-1:nxd,-1:nyd,0,is) = am(is) * mp * ua(-1:nxd,-1:nyd,is) *
     .     fna(-1:nxd,-1:nyd,0,is)
          fmo(-1:nxd,-1:nyd,1,is) = am(is) * mp * ua(-1:nxd,-1:nyd,is) *
     .     fna(-1:nxd,-1:nyd,1,is)
        enddo
        fna_32 = fna
        fna_52 = 0.0_R8
        fni_32 = fni
        fni_52 = 0.0_R8
        fne_32 = fne
        fne_52 = 0.0_R8
      endif
      if(b2fplasma_version.ge.'03.000.006') then !djm 17/10/16
        call cfrure (nread, size(fna_pinch), fna_pinch, 'fna_pinch')
        call cfrure (nread, size(fna_pll), fna_pll, 'fna_pll')
        call cfrure (nread, size(fna_drift),fna_drift , 'fna_drift')
        call cfrure (nread, size(fna_ch), fna_ch, 'fna_ch')
        call cfrure (nread, size(fna_nanom), fna_nanom, 'fna_nanom')
        call cfrure (nread, size(fna_panom), fna_panom, 'fna_panom')
      endif
      if(b2fplasma_version.ge.'01.001.010') then
        call cfrure (nread, size(wadia), wadia, 'wadia')
        call cfrure (nread, size(vaecrb), vaecrb, 'vaecrb')
      else
        wadia = uadia
        vaecrb = 0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.012') then
        call cfrure (nread, size(facdrift), facdrift, 'facdrift')
        call cfrure (nread, size(fac_ExB), fac_ExB, 'fac_ExB')
      else
        facdrift = 0.0_R8
        fac_ExB = 0.0_R8
      endif
      if(b2fplasma_version.ge.'01.001.014') then
        call cfrure (nread, size(fchvispar), fchvispar, 'fchvispar')
        call cfrure (nread, size(fchvisper), fchvisper, 'fchvisper')
        call cfrure (nread, size(fchin), fchin, 'fchin')
        call cfrure (nread, size(fna_nodrift), fna_nodrift,
     1   'fna_nodrift')
      else
        fchvispar = 0.0_R8
        fchvisper = 0.0_R8
        fchin = 0.0_R8
        fna_nodrift = fna
      endif
      if(b2fplasma_version.ge.'01.001.015') then
        call cfrure (nread, size(fac_vis), fac_vis, 'fac_vis')
      else
        fac_vis = 0.0_R8
      endif
      if(b2fplasma_version.ge.'03.000.003') then
        call cfrure (nread, size(fna_mdf), fna_mdf, 'fna_mdf')
        call cfrure (nread, size(fhe_mdf), fhe_mdf, 'fhe_mdf')
        call cfrure (nread, size(fhi_mdf), fhi_mdf, 'fhi_mdf')
        call cfrure (nread, size(fnaPSch), fnaPSch, 'fnaPSch')
        call cfrure (nread, size(fhePSch), fhePSch, 'fhePSch')
        call cfrure (nread, size(fhiPSch), fhiPSch, 'fhiPSch')
        call cfrure (nread, size(fna_fcor), fna_fcor, 'fna_fcor')
        call cfrure (nread, size(fna_he), fna_he, 'fna_he')
        call cfrure (nread, size(fchvisq), fchvisq, 'fchvisq')
        call cfrure (nread, size(fchinert), fchinert, 'fchinert')
      endif

      return
      end subroutine read_b2mod_plasma

      subroutine write_b2mod_plasma(nwrite)
      implicit none
      integer nwrite

      call cfwure (nwrite, size(fch), fch, 'fch')
      call cfwure (nwrite, size(fch0), fch0, 'fch0')
      call cfwure (nwrite, size(fchp), fchp, 'fchp')
      call cfwure (nwrite, size(fhe), fhe, 'fhe')
      call cfwure (nwrite, size(fhe0), fhe0, 'fhe0')
      call cfwure (nwrite, size(fhep), fhep, 'fhep')
      if(newversion.gt.'01.001.003') then
        write(*,*) newversion,'fhet --- fhit'
        call cfwure (nwrite, size(fhet), fhet, 'fhet')
      endif
      call cfwure (nwrite, size(fhi), fhi, 'fhi')
      call cfwure (nwrite, size(fhi0), fhi0, 'fhi0')
      call cfwure (nwrite, size(fhip), fhip, 'fhip')
      if(newversion.gt.'01.001.003') then
        call cfwure (nwrite, size(fhit), fhit, 'fhit')
      endif
      call cfwure (nwrite, size(fna), fna, 'fna')
      call cfwure (nwrite, size(fna0), fna0, 'fna0')
      call cfwure (nwrite, size(fnap), fnap, 'fnap')
      call cfwure (nwrite, size(fne), fne, 'fne')
      call cfwure (nwrite, size(fni), fni, 'fni')
      call cfwure (nwrite, size(na), na, 'na')
      call cfwure (nwrite, size(na0), na0, 'na0')
      call cfwure (nwrite, size(nap), nap, 'nap')
      call cfwure (nwrite, size(ne), ne, 'ne')
      call cfwure (nwrite, size(ne0), ne0, 'ne0')
      call cfwure (nwrite, size(ne2), ne2, 'ne2')
      call cfwure (nwrite, size(nep), nep, 'nep')
      call cfwure (nwrite, size(ni), ni, 'ni')
      call cfwure (nwrite, size(ni0), ni0, 'ni0')
      call cfwure (nwrite, size(pb), pb, 'pb')
      call cfwure (nwrite, size(po), po, 'po')
      call cfwure (nwrite, size(po0), po0, 'po0')
      call cfwure (nwrite, size(pop), pop, 'pop')
      call cfwure (nwrite, size(te), te, 'te')
      call cfwure (nwrite, size(te0), te0, 'te0')
      call cfwure (nwrite, size(tep), tep, 'tep')
      call cfwure (nwrite, size(ti), ti, 'ti')
      call cfwure (nwrite, size(ti0), ti0, 'ti0')
      call cfwure (nwrite, size(tip), tip, 'tip')
      call cfwure (nwrite, size(ua), ua, 'ua')
      call cfwure (nwrite, size(ua0), ua0, 'ua0')
      call cfwure (nwrite, size(uap), uap, 'uap')
      call cfwure (nwrite, size(uadia), uadia, 'uadia')
      call cfwure (nwrite, size(fchdia), fchdia, 'fchdia')
      if(newversion.ge.'01.001.008') then
        call cfwure (nwrite, size(fmo), fmo, 'fmo')
        call cfwure (nwrite, size(fna_32), fna_32, 'fna_32')
        call cfwure (nwrite, size(fna_52), fna_52, 'fna_52')
        call cfwure (nwrite, size(fni_32), fni_32, 'fni_32')
        call cfwure (nwrite, size(fni_52), fni_52, 'fni_52')
        call cfwure (nwrite, size(fne_32), fne_32, 'fne_32')
        call cfwure (nwrite, size(fne_52), fne_52, 'fne_52')
      endif
      if(newversion.ge.'03.000.006') then !djm 17/10/16
        call cfwure (nwrite, size(fna_pinch), fna_pinch, 'fna_pinch')
        call cfwure (nwrite, size(fna_pll), fna_pll, 'fna_pll')
        call cfwure (nwrite, size(fna_drift),fna_drift , 'fna_drift')
        call cfwure (nwrite, size(fna_ch), fna_ch, 'fna_ch')
        call cfwure (nwrite, size(fna_nanom), fna_nanom, 'fna_nanom')
        call cfwure (nwrite, size(fna_panom), fna_panom, 'fna_panom')
      endif
      if(newversion.ge.'01.001.010') then
        call cfwure (nwrite, size(wadia), wadia, 'wadia')
        call cfwure (nwrite, size(vaecrb), vaecrb, 'vaecrb')
      endif
      if(newversion.ge.'01.001.012') then
        call cfwure (nwrite, size(facdrift), facdrift, 'facdrift')
        call cfwure (nwrite, size(fac_ExB), fac_ExB, 'fac_ExB')
      endif
      if(newversion.ge.'01.001.014') then
        call cfwure (nwrite, size(fchvispar), fchvispar, 'fchvispar')
        call cfwure (nwrite, size(fchvisper), fchvisper, 'fchvisper')
        call cfwure (nwrite, size(fchin), fchin, 'fchin')
        call cfwure (nwrite, size(fna_nodrift), fna_nodrift,
     1   'fna_nodrift')
      endif
      if(newversion.ge.'01.001.015') then
        call cfwure (nwrite, size(fac_vis), fac_vis, 'fac_vis')
      endif
      if(newversion.ge.'03.000.003') then
        call cfwure (nwrite, size(fna_mdf), fna_mdf, 'fna_mdf')
        call cfwure (nwrite, size(fhe_mdf), fhe_mdf, 'fhe_mdf')
        call cfwure (nwrite, size(fhi_mdf), fhi_mdf, 'fhi_mdf')
        call cfwure (nwrite, size(fnaPSch), fnaPSch, 'fnaPSch')
        call cfwure (nwrite, size(fhePSch), fhePSch, 'fhePSch')
        call cfwure (nwrite, size(fhiPSch), fhiPSch, 'fhiPSch')
        call cfwure (nwrite, size(fna_fcor), fna_fcor, 'fna_fcor')
        call cfwure (nwrite, size(fna_he), fna_he, 'fna_he')
        call cfwure (nwrite, size(fchvisq), fchvisq, 'fchvisq')
        call cfwure (nwrite, size(fchinert), fchinert, 'fchinert')
      endif

      return
      end subroutine write_b2mod_plasma

      end module b2mod_plasma
