      subroutine tallies(nx,ny,ns,BoRiS)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_work
      use b2mod_tallies
      use b2mod_indirect
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_time
      use b2mod_diag
      use b2mod_b2cmpa
      implicit none
      integer nx,ny,ns
      real (kind=R8) :: BoRiS
      real (kind=R8), save :: potave(0:nsdecl-1)
#ifndef NO_CDF
#     include <netcdf.inc>
      integer status, ncid
      real (kind=R8) :: tim_arr(1)
      real (kind=R8) :: naregtmp(0:nnreg(0),0:ns-1)
#endif
*   ..common blocks
      integer ncall
      integer ix,iy,is,ir,tallies_netcdf
      character*256 tally_format, species_format
      character*10, allocatable, save :: species_name(:)
      save ncall,tallies_netcdf
      data ncall/0/,tallies_netcdf/0/
#ifdef NO_LONG_LINES
      data tally_format, species_format
     *  /'(a,1p,8e12.4,/,(''+ '',10e12.4))',
     *   '(a,1x,i3,18(1x,a),/,(''+ '',20(1x,a)))'/
#else
      data tally_format, species_format
     *  /'(a,1p,(10000e12.4))','(a,1x,i3,999(1x,a))'/
#endif
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('tallies')
*    ..initialise species_name
      if (ncall.eq.0) then
        call initialise_potentials(ns)
        allocate(species_name(0:ns-1))
        do is = 0, ns-1
          potave(is) = (potmin(is)+potmax(is))/2.0_R8
          call species(is,species_name(is))
        enddo
        call ipgeti('tallies_netcdf',tallies_netcdf)
      endif
*    ..produce tallies
      fnaxreg=0.0_R8
      fnayreg=0.0_R8
      fhpxreg=0.0_R8
      fhpyreg=0.0_R8
      fhmxreg=0.0_R8
      fhmyreg=0.0_R8
      fhixreg=0.0_R8
      fhiyreg=0.0_R8
      fhexreg=0.0_R8
      fheyreg=0.0_R8
      fchxreg=0.0_R8
      fchyreg=0.0_R8
      fhtxreg=0.0_R8                                               !cjk
      fhtyreg=0.0_R8
      nareg=0.0_R8
      tereg=0.0_R8
      nereg=0.0_R8
      ne2reg=0.0_R8
      tireg=0.0_R8
      nireg=0.0_R8
      poreg=0.0_R8
      volreg=0.0_R8
      qconvixreg=0.0_R8                                            !cjk
      qconviyreg=0.0_R8
      qconvexreg=0.0_R8
      qconveyreg=0.0_R8
      b2stbr_sna_reg=0.0_R8
      b2stbr_smo_reg=0.0_R8
      b2stbr_she_reg=0.0_R8
      b2stbr_shi_reg=0.0_R8
      b2stbr_sch_reg=0.0_R8
      b2stbc_sna_reg=0.0_R8
      b2stbc_smo_reg=0.0_R8
      b2stbc_she_reg=0.0_R8
      b2stbc_shi_reg=0.0_R8
      b2stbc_sch_reg=0.0_R8
      resco_reg=0.0_R8
      reshe_reg=0.0_R8
      reshi_reg=0.0_R8
      resmo_reg=0.0_R8
      resmt_reg=0.0_R8
      respo_reg=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          do is=0,ns-1
            fnaxreg(region(ix,iy,1),is)=fnaxreg(region(ix,iy,1),is)+
     1       fna(ix,iy,0,is)
            fnayreg(region(ix,iy,2),is)=fnayreg(region(ix,iy,2),is)+
     1       fna(ix,iy,1,is)
            nareg(region(ix,iy,0),is)=nareg(region(ix,iy,0),is)+
     1       na(ix,iy,is)*vol(ix,iy)
            b2stbr_sna_reg(region(ix,iy,0),is)=
     1       b2stbr_sna_reg(region(ix,iy,0),is)+
     1       b2stbr_sna(ix,iy,is)
            b2stbr_smo_reg(region(ix,iy,0),is)=
     1       b2stbr_smo_reg(region(ix,iy,0),is)+
     1       b2stbr_smo(ix,iy,is)
            b2stbc_sna_reg(region(ix,iy,0),is)=
     1       b2stbc_sna_reg(region(ix,iy,0),is)+
     1       b2stbc_sna(ix,iy,is)
            b2stbc_smo_reg(region(ix,iy,0),is)=
     1       b2stbc_smo_reg(region(ix,iy,0),is)+
     1       b2stbc_smo(ix,iy,is)
            resco_reg(region(ix,iy,0),is)=resco_reg(region(ix,iy,0),is)+
     1       resco(ix,iy,is)
            resmo_reg(region(ix,iy,0),is)=resmo_reg(region(ix,iy,0),is)+
     1       resmo(ix,iy,is)
            if (leftix(ix,iy).ne.-2) then
              fhmxreg(region(ix,iy,1),is)=fhmxreg(region(ix,iy,1),is)+
     1         0.5_R8*(kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,0,is)
              fhpxreg(region(ix,iy,1),is)=fhpxreg(region(ix,iy,1),is)+
     1         0.5_R8*(rpt(leftix(ix,iy),leftiy(ix,iy),is)+
     2                 rpt(ix,iy,is))*fna(ix,iy,0,is)*ev
              fhixreg(region(ix,iy,1))=fhixreg(region(ix,iy,1))+
     .         (1.0_R8-BoRiS)*
     1         0.5_R8*(kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,0,is)
              fhtxreg(region(ix,iy,1))=fhtxreg(region(ix,iy,1))+         !cjk
     1         0.5_R8*(kinrgy(leftix(ix,iy),leftiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,0,is)+
     1         0.5_R8*(rpt(leftix(ix,iy),leftiy(ix,iy),is)+
     2                 rpt(ix,iy,is))*fna(ix,iy,0,is)*ev
              qconvixreg(region(ix,iy,1),is)=                            !cjk
     1         qconvixreg(region(ix,iy,1),is)+
     1         (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))/2.0_R8*
     1         (1.5_R8*fna_32(ix,iy,0,is)+2.5_R8*fna_52(ix,iy,0,is))         
              qconvexreg(region(ix,iy,1))=                               !cjk
     1         qconvexreg(region(ix,iy,1))+
     1         (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))/2.0_R8*
     1         (1.5_R8*fna_32(ix,iy,0,is)+
     1          2.5_R8*fna_52(ix,iy,0,is))*
     2         (rza(leftix(ix,iy),leftiy(ix,iy),is)+
     2          rza(ix,iy,is))/2.0_R8         
            endif
            if (bottomiy(ix,iy).ne.-2) then
              fhmyreg(region(ix,iy,2),is)=fhmyreg(region(ix,iy,2),is)+
     1         0.5_R8*(kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,1,is)
              fhpyreg(region(ix,iy,2),is)=fhpyreg(region(ix,iy,2),is)+
     1         0.5_R8*(rpt(bottomix(ix,iy),bottomiy(ix,iy),is)+
     2                 rpt(ix,iy,is))*fna(ix,iy,1,is)*ev
              fhiyreg(region(ix,iy,2))=fhiyreg(region(ix,iy,2))+
     .         (1.0_R8-BoRiS)*
     1         0.5_R8*(kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,1,is)
              fhtyreg(region(ix,iy,2))=fhtyreg(region(ix,iy,2))+
     1         0.5_R8*(kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)+
     2                 kinrgy(ix,iy,is))*fna(ix,iy,1,is)+
     1         0.5_R8*(rpt(bottomix(ix,iy),bottomiy(ix,iy),is)+
     2                 rpt(ix,iy,is))*fna(ix,iy,1,is)*ev
              qconviyreg(region(ix,iy,2),is)=
     1         qconviyreg(region(ix,iy,2),is)+
     1         (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/2.0_R8*
     1         (1.5_R8*fna_32(ix,iy,1,is)+2.5_R8*fna_52(ix,iy,1,is))         
              qconveyreg(region(ix,iy,2))=
     1         qconveyreg(region(ix,iy,2))+
     1         (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/2.0_R8*
     1         (1.5_R8*fna_32(ix,iy,1,is)+
     1          2.5_R8*fna_52(ix,iy,1,is))*
     2         (rza(bottomix(ix,iy),bottomiy(ix,iy),is)+rza(ix,iy,is))
     2         /2.0_R8                  
            endif
          enddo
          fhixreg(region(ix,iy,1))=fhixreg(region(ix,iy,1))+
     1     fhi(ix,iy,0)
          fhiyreg(region(ix,iy,2))=fhiyreg(region(ix,iy,2))+
     1     fhi(ix,iy,1)
          fhexreg(region(ix,iy,1))=fhexreg(region(ix,iy,1))+
     1     fhe(ix,iy,0)
          fheyreg(region(ix,iy,2))=fheyreg(region(ix,iy,2))+
     1     fhe(ix,iy,1)
          fchxreg(region(ix,iy,1))=fchxreg(region(ix,iy,1))+
     1     fch(ix,iy,0)
          fchyreg(region(ix,iy,2))=fchyreg(region(ix,iy,2))+
     1     fch(ix,iy,1)
          fhtxreg(region(ix,iy,1))=fhtxreg(region(ix,iy,1))+             !cjk
     1     fhi(ix,iy,0)+fhe(ix,iy,0)
          fhtyreg(region(ix,iy,2))=fhtyreg(region(ix,iy,2))+
     1     fhi(ix,iy,1)+fhe(ix,iy,1)
          tereg(region(ix,iy,0))=tereg(region(ix,iy,0))+
     1     te(ix,iy)*ne(ix,iy)*vol(ix,iy)
          nereg(region(ix,iy,0))=nereg(region(ix,iy,0))+
     1     ne(ix,iy)*vol(ix,iy)
          ne2reg(region(ix,iy,0))=ne2reg(region(ix,iy,0))+
     1     ne2(ix,iy)*vol(ix,iy)
          tireg(region(ix,iy,0))=tireg(region(ix,iy,0))+
     1     ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy)
          nireg(region(ix,iy,0))=nireg(region(ix,iy,0))+
     1     ni(ix,iy,0)*vol(ix,iy)
          poreg(region(ix,iy,0))=poreg(region(ix,iy,0))+
     1     po(ix,iy)*vol(ix,iy)
          volreg(region(ix,iy,0))=volreg(region(ix,iy,0))+
     1     vol(ix,iy)
          b2stbr_she_reg(region(ix,iy,0))=
     1     b2stbr_she_reg(region(ix,iy,0))+
     1     b2stbr_she(ix,iy)
          b2stbr_shi_reg(region(ix,iy,0))=
     1     b2stbr_shi_reg(region(ix,iy,0))+
     1     b2stbr_shi(ix,iy)
          b2stbr_sch_reg(region(ix,iy,0))=
     1     b2stbr_sch_reg(region(ix,iy,0))+
     1     b2stbr_sch(ix,iy)
          b2stbc_she_reg(region(ix,iy,0))=
     1     b2stbc_she_reg(region(ix,iy,0))+
     1     b2stbc_she(ix,iy)
          b2stbc_shi_reg(region(ix,iy,0))=
     1     b2stbc_shi_reg(region(ix,iy,0))+
     1     b2stbc_shi(ix,iy)
          b2stbc_sch_reg(region(ix,iy,0))=
     1     b2stbc_sch_reg(region(ix,iy,0))+
     1     b2stbc_sch(ix,iy)
          reshe_reg(region(ix,iy,0))=reshe_reg(region(ix,iy,0))+
     1     reshe(ix,iy)
          reshi_reg(region(ix,iy,0))=reshi_reg(region(ix,iy,0))+
     1     reshi(ix,iy)
          resmt_reg(region(ix,iy,0))=resmt_reg(region(ix,iy,0))+
     1     resmt(ix,iy)
          respo_reg(region(ix,iy,0))=respo_reg(region(ix,iy,0))+
     1     respo(ix,iy)
        enddo
      enddo
      do is=0,ns-1
        rsanareg(0,is)=0.0_R8
        rsahireg(0,is)=0.0_R8
        rsamoreg(0,is)=0.0_R8
        rranareg(0,is)=0.0_R8
        rrahireg(0,is)=0.0_R8
        rramoreg(0,is)=0.0_R8
        rqahereg(0,is)=0.0_R8
        rqradreg(0,is)=0.0_R8
        rqbrmreg(0,is)=0.0_R8
        rcxnareg(0,is)=0.0_R8
        rcxhireg(0,is)=0.0_R8
        rcxmoreg(0,is)=0.0_R8
        b2stbr_sna_reg(0,is)=0.0_R8
        b2stbr_smo_reg(0,is)=0.0_R8
        b2stbc_sna_reg(0,is)=0.0_R8
        b2stbc_smo_reg(0,is)=0.0_R8
        nareg(0,is)=0.0_R8
        resco_reg(0,is)=0.0_R8
        resmo_reg(0,is)=0.0_R8
        do ir=1,nnreg(0)
          rsanareg(0,is)=rsanareg(0,is)+rsanareg(ir,is)
          rsahireg(0,is)=rsahireg(0,is)+rsahireg(ir,is)
          rsamoreg(0,is)=rsamoreg(0,is)+rsamoreg(ir,is)
          rranareg(0,is)=rranareg(0,is)+rranareg(ir,is)
          rrahireg(0,is)=rrahireg(0,is)+rrahireg(ir,is)
          rramoreg(0,is)=rramoreg(0,is)+rramoreg(ir,is)
          rqahereg(0,is)=rqahereg(0,is)+rqahereg(ir,is)
          rqradreg(0,is)=rqradreg(0,is)+rqradreg(ir,is)
          rqbrmreg(0,is)=rqbrmreg(0,is)+rqbrmreg(ir,is)
          rcxnareg(0,is)=rcxnareg(0,is)+rcxnareg(ir,is)
          rcxhireg(0,is)=rcxhireg(0,is)+rcxhireg(ir,is)
          rcxmoreg(0,is)=rcxmoreg(0,is)+rcxmoreg(ir,is)
          nareg(0,is)=nareg(0,is)+nareg(ir,is)
          b2stbr_sna_reg(0,is)=b2stbr_sna_reg(0,is)+
     1     b2stbr_sna_reg(ir,is)
          b2stbr_smo_reg(0,is)=b2stbr_smo_reg(0,is)+
     1     b2stbr_smo_reg(ir,is)
          b2stbc_sna_reg(0,is)=b2stbc_sna_reg(0,is)+
     1     b2stbc_sna_reg(ir,is)
          b2stbc_smo_reg(0,is)=b2stbc_smo_reg(0,is)+
     1     b2stbc_smo_reg(ir,is)
          resco_reg(0,is)=resco_reg(0,is)+resco_reg(ir,is)
          resmo_reg(0,is)=resmo_reg(0,is)+resmo_reg(ir,is)
        enddo
      enddo
      do is=0,ns-1
        fnaxreg(0,is)=0.0_R8
        fnayreg(0,is)=0.0_R8
        fhpxreg(0,is)=0.0_R8
        fhpyreg(0,is)=0.0_R8
        fhmxreg(0,is)=0.0_R8
        fhmyreg(0,is)=0.0_R8
        qconvixreg(0,is)=0.0_R8                                             !cjk
        qconviyreg(0,is)=0.0_R8
      enddo
      fhixreg(0)=0.0_R8
      fhiyreg(0)=0.0_R8
      fhexreg(0)=0.0_R8
      fheyreg(0)=0.0_R8
      fhtxreg(0)=0.0_R8
      fhtyreg(0)=0.0_R8
      fchxreg(0)=0.0_R8
      fchyreg(0)=0.0_R8
      qconvexreg(0)=0.0_R8                                                  !cjk
      qconveyreg(0)=0.0_R8
      tereg(0)=0.0_R8
      nereg(0)=0.0_R8
      ne2reg(0)=0.0_R8
      tireg(0)=0.0_R8
      nireg(0)=0.0_R8
      poreg(0)=0.0_R8
      volreg(0)=0.0_R8
      b2stbr_she_reg(0)=0.0_R8
      b2stbr_shi_reg(0)=0.0_R8
      b2stbr_sch_reg(0)=0.0_R8
      b2stbc_she_reg(0)=0.0_R8
      b2stbc_shi_reg(0)=0.0_R8
      b2stbc_sch_reg(0)=0.0_R8
      b2divue(0)=0.0_R8
      b2divua(0)=0.0_R8
      b2exbe(0)=0.0_R8
      b2exba(0)=0.0_R8
      b2visa(0)=0.0_R8
      b2joule(0)=0.0_R8
      b2fraa(0)=0.0_R8
      b2str(0)=0.0_R8
      b2wrong1(0)=0.0_R8
      b2wrong2(0)=0.0_R8
      b2wrong3(0)=0.0_R8
      reshe_reg(0)=0.0_R8
      reshi_reg(0)=0.0_R8
      resmt_reg(0)=0.0_R8
      respo_reg(0)=0.0_R8
      do ir=1,nnreg(0)
        tereg(0)=tereg(0)+tereg(ir)
        nereg(0)=nereg(0)+nereg(ir)
        ne2reg(0)=ne2reg(0)+ne2reg(ir)
        tireg(0)=tireg(0)+tireg(ir)
        nireg(0)=nireg(0)+nireg(ir)
        poreg(0)=poreg(0)+poreg(ir)
        volreg(0)=volreg(0)+volreg(ir)
        b2stbr_she_reg(0)=b2stbr_she_reg(0)+b2stbr_she_reg(ir)
        b2stbr_shi_reg(0)=b2stbr_shi_reg(0)+b2stbr_shi_reg(ir)
        b2stbr_sch_reg(0)=b2stbr_sch_reg(0)+b2stbr_sch_reg(ir)
        b2stbc_she_reg(0)=b2stbc_she_reg(0)+b2stbc_she_reg(ir)
        b2stbc_shi_reg(0)=b2stbc_shi_reg(0)+b2stbc_shi_reg(ir)
        b2stbc_sch_reg(0)=b2stbc_sch_reg(0)+b2stbc_sch_reg(ir)
        b2qie(0)=b2qie(0)+b2qie(ir)
        b2divue(0)=b2divue(0)+b2divue(ir)
        b2divua(0)=b2divua(0)+b2divua(ir)
        b2exbe(0)=b2exbe(0)+b2exbe(ir)
        b2exba(0)=b2exba(0)+b2exba(ir)
        b2visa(0)=b2visa(0)+b2visa(ir)
        b2joule(0)=b2joule(0)+b2joule(ir)
        b2fraa(0)=b2fraa(0)+b2fraa(ir)
        b2str(0)=b2str(0)+b2str(ir)
        b2wrong1(0)=b2wrong1(0)+b2wrong1(ir)
        b2wrong2(0)=b2wrong2(0)+b2wrong2(ir)
        b2wrong3(0)=b2wrong3(0)+b2wrong3(ir)
        b2she(0)=b2she(0)+b2she(ir)
        b2shi(0)=b2shi(0)+b2shi(ir)
        b2she0(0)=b2she0(0)+b2she0(ir)
        b2shi0(0)=b2shi0(0)+b2shi0(ir)
        reshe_reg(0)=reshe_reg(0)+reshe_reg(ir)
        reshi_reg(0)=reshi_reg(0)+reshi_reg(ir)
        resmt_reg(0)=resmt_reg(0)+resmt_reg(ir)
        respo_reg(0)=respo_reg(0)+respo_reg(ir)
      enddo
*    ..write tallies
      write(*,tally_format) 'Start tallies by region'
      write(*,'(a,3(1x,i2))') 'nnreg',nnreg
      write(*,species_format) 'species',ns,
     * (trim(species_name(is)),is=0,ns-1)
      if(ncall.eq.0) then
        write(*,tally_format) 'pot',potave(0:ns-1)
      endif
      write(*,tally_format) 'rsanareg',rsanareg
      write(*,tally_format) 'rsahireg',rsahireg
      write(*,tally_format) 'rsamoreg',rsamoreg
      write(*,tally_format) 'rranareg',rranareg
      write(*,tally_format) 'rrahireg',rrahireg
      write(*,tally_format) 'rramoreg',rramoreg
      write(*,tally_format) 'rqahereg',rqahereg
      write(*,tally_format) 'rqradreg',rqradreg
      write(*,tally_format) 'rqbrmreg',rqbrmreg
      write(*,tally_format) 'rcxnareg',rcxnareg
      write(*,tally_format) 'rcxhireg',rcxhireg
      write(*,tally_format) 'rcxmoreg',rcxmoreg
      write(*,tally_format) 'fnaxreg',fnaxreg
      write(*,tally_format) 'fnayreg',fnayreg
      write(*,tally_format) 'fhixreg',fhixreg
      write(*,tally_format) 'fhiyreg',fhiyreg
      write(*,tally_format) 'fhexreg',fhexreg
      write(*,tally_format) 'fheyreg',fheyreg
      write(*,tally_format) 'fhpxreg',fhpxreg
      write(*,tally_format) 'fhpyreg',fhpyreg
      write(*,tally_format) 'fhmxreg',fhmxreg
      write(*,tally_format) 'fhmyreg',fhmyreg
      write(*,tally_format) 'fchxreg',fchxreg
      write(*,tally_format) 'fchyreg',fchyreg
      write(*,tally_format) 'fhtxreg',fhtxreg                          !cjk
      write(*,tally_format) 'fhtyreg',fhtyreg
      write(*,tally_format) 'qconvixreg',qconvixreg
      write(*,tally_format) 'qconviyreg',qconviyreg
      write(*,tally_format) 'qconvexreg',qconvexreg
      write(*,tally_format) 'qconveyreg',qconveyreg
      write(*,tally_format) 'b2stbr_sna_reg',b2stbr_sna_reg
      write(*,tally_format) 'b2stbr_smo_reg',b2stbr_smo_reg
      write(*,tally_format) 'b2stbr_she_reg',b2stbr_she_reg
      write(*,tally_format) 'b2stbr_shi_reg',b2stbr_shi_reg
      write(*,tally_format) 'b2stbr_sch_reg',b2stbr_sch_reg
      write(*,tally_format) 'b2stbc_sna_reg',b2stbc_sna_reg
      write(*,tally_format) 'b2stbc_smo_reg',b2stbc_smo_reg
      write(*,tally_format) 'b2stbc_she_reg',b2stbc_she_reg
      write(*,tally_format) 'b2stbc_shi_reg',b2stbc_shi_reg
      write(*,tally_format) 'b2stbc_sch_reg',b2stbc_sch_reg
      write(*,tally_format) 'nareg',nareg
      write(*,tally_format) 'tereg',tereg
      write(*,tally_format) 'nereg',nereg
      write(*,tally_format) 'ne2reg',ne2reg
      write(*,tally_format) 'tireg',tireg
      write(*,tally_format) 'nireg',nireg
      write(*,tally_format) 'poreg',poreg
      write(*,tally_format) 'volreg',volreg
      write(*,tally_format) 'b2brem',b2brem
      write(*,tally_format) 'b2rad',b2rad
      write(*,tally_format) 'b2qie',b2qie
      write(*,tally_format) 'b2vdp',b2vdp
      write(*,tally_format) 'b2divue',b2divue
      write(*,tally_format) 'b2divua',b2divua
      write(*,tally_format) 'b2exbe',b2exbe
      write(*,tally_format) 'b2exba',b2exba
      write(*,tally_format) 'b2visa',b2visa
      write(*,tally_format) 'b2joule',b2joule
      write(*,tally_format) 'b2fraa',b2fraa
      write(*,tally_format) 'b2str',b2str
      write(*,tally_format) 'b2wrong1',b2wrong1
      write(*,tally_format) 'b2wrong2',b2wrong2
      write(*,tally_format) 'b2wrong3',b2wrong3
      write(*,tally_format) 'b2she',b2she
      write(*,tally_format) 'b2shi',b2shi
      write(*,tally_format) 'b2she0',b2she0
      write(*,tally_format) 'b2shi0',b2shi0
      write(*,tally_format) 'resco_reg',resco_reg
      write(*,tally_format) 'resmo_reg',resmo_reg
      write(*,tally_format) 'reshe_reg',reshe_reg
      write(*,tally_format) 'reshi_reg',reshi_reg
      write(*,tally_format) 'resmt_reg',resmt_reg
      write(*,tally_format) 'respo_reg',respo_reg
      write(*,tally_format) 'End tallies by region'

#ifndef NO_CDF
      if(tallies_netcdf.ne.0) then
        if(ncall.eq.0) then
          call define_cdf_b2mod_tallies(nnreg,ns)
        endif
        write(*,tally_format) 'Start (CDF) tallies by region'
        status=nf_open('b2tallies.nc',NF_WRITE,ncid)
        call check_cdf_status(status)

        tim_arr=tim
        call write_cdf_variable(ncid,ncall+1,
     *   'times',tim_arr)
        call write_cdf_variable(ncid,ncall+1,
     *   'rsanareg',rsanareg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rsahireg',rsahireg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rsamoreg',rsamoreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rranareg',rranareg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rrahireg',rrahireg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rramoreg',rramoreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rqahereg',rqahereg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rqradreg',rqradreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rqbrmreg',rqbrmreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rcxnareg',rcxnareg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rcxhireg',rcxhireg)
        call write_cdf_variable(ncid,ncall+1,
     *   'rcxmoreg',rcxmoreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fnaxreg',fnaxreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fnayreg',fnayreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhixreg',fhixreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhiyreg',fhiyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhexreg',fhexreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fheyreg',fheyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhpxreg',fhpxreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhpyreg',fhpyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhmxreg',fhmxreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhmyreg',fhmyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fchxreg',fchxreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fchyreg',fchyreg)
        call write_cdf_variable(ncid,ncall+1,                          !cjk
     *   'fhtxreg',fhtxreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'fhtyreg',fhtyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'qconvixreg',qconvixreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'qconviyreg',qconviyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'qconvexreg',qconvexreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'qconveyreg',qconveyreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbr_sna_reg',b2stbr_sna_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbr_smo_reg',b2stbr_smo_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbr_she_reg',b2stbr_she_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbr_shi_reg',b2stbr_shi_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbr_sch_reg',b2stbr_sch_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbc_sna_reg',b2stbc_sna_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbc_smo_reg',b2stbc_smo_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbc_she_reg',b2stbc_she_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbc_shi_reg',b2stbc_shi_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2stbc_sch_reg',b2stbc_sch_reg)
	do is=0,ns-1
	  naregtmp(:,is)=nareg(:,is)/volreg(:)
	enddo
        call write_cdf_variable(ncid,ncall+1,
     *   'nareg',naregtmp)
        call write_cdf_variable(ncid,ncall+1,
     *   'tereg',tereg/nereg/qe)
        call write_cdf_variable(ncid,ncall+1,
     *   'nereg',nereg/volreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'ne2reg',ne2reg/volreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'tireg',tireg/nireg/qe)
        call write_cdf_variable(ncid,ncall+1,
     *   'nireg',nireg/volreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'poreg',poreg/volreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'volreg',volreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2bremreg',b2bremreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2radreg',b2radreg)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2qie',b2qie)
!       call write_cdf_variable(ncid,ncall+1,
!     *  'b2vdp',b2vdp)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2divue',b2divue)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2divua',b2divua)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2exbe',b2exbe)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2exba',b2exba)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2visa',b2visa)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2joule',b2joule)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2fraa',b2fraa)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2str',b2str)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2wrong1',b2wrong1)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2wrong2',b2wrong2)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2wrong3',b2wrong3)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2she',b2she)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2shi',b2shi)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2she0',b2she0)
        call write_cdf_variable(ncid,ncall+1,
     *   'b2shi0',b2shi0)
        call write_cdf_variable(ncid,ncall+1,
     *   'resco_reg',resco_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'resmo_reg',resmo_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'reshe_reg',reshe_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'reshi_reg',reshi_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'resmt_reg',resmt_reg)
        call write_cdf_variable(ncid,ncall+1,
     *   'respo_reg',respo_reg)
        write(*,tally_format) 'End (CDF) tallies by region'

        status=nf_close(ncid)
        call check_cdf_status(status)
      endif
#endif
* ..return
      ncall = ncall+1
      call subend ()
      return
      end subroutine tallies
