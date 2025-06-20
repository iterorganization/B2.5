      subroutine cdfmovie
#ifndef NO_CDF
     . (ncid,nCv,nCi,ns,geo,state,switch)
      use b2mod_types
      use b2mod_time
      use b2us_geo
      use b2us_plasma
      use b2mod_dimensions
      use b2mod_movies
      use b2mod_constants
#ifdef B25_EIRENE
      use eirmod_parmmod
     ,   , ncv_eir => ncv
      use eirmod_wneutrals
      use eirmod_extrab25
      use eirmod_cestim
      use b2mod_balance
      use eirmod_eirbra
#else
      use b2mod_eirdiag
#endif
      use b2mod_switches
      use b2mod_subsys
      implicit none
      integer ncid, nCv, nCi, ns
      type (geometry), intent(in) :: geo
      type (B2State), intent(in) :: state
      type (switches), intent(in) :: switch
#     include <netcdf.inc>
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      real (kind=R8) :: time(1), dummy(1:nCv), dummy3(1:nCv, 1:ns)
#ifdef B25_EIRENE
      integer nstrad
      integer iatmi, imoli, istra
      real (kind=R8) :: dummytriatm(1:ncvtot-1, 1:nnatmi),
     & dummytrimol(1:ncvtot-1, 1:nnmoli),
     & dummy3atm(1:nCv, 1:nnatmi), dummy3mol(1:nCv, 1:nnmoli)
      character*256, save :: filename_tri
#endif
      integer is
      integer ncall
      integer status, ndims, nvars, natts, unlimid
      integer, save :: write_b25_state_variables = 1
      integer, save :: write_b25_sources = 0
      integer, save :: write_eirene_state_variables = 0
      integer, save :: write_eirene_sources = 0
      integer, save :: write_triangles = 0
      character*256, save :: filename
      logical ex
      save ncall
      data ncall/0/

      call ipgeti ('cdfmovie_b25_state_variables',
     1  write_b25_state_variables)
      call ipgeti ('cdfmovie_b25_sources', write_b25_sources)
      call ipgeti ('cdfmovie_eirene_state_variables',
     1  write_eirene_state_variables)
      call ipgeti ('cdfmovie_eirene_sources', write_eirene_sources)
      call ipgeti ('cdfmovie_triangles', write_triangles)

      call subini ('cdfmovie')
      if (ncall.eq.0) then
        filename='b2movies.nc'
        call find_file(filename,ex)
        if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
          write(6,'(a)') trim(filename)//' will be appended'
          status=nf_open(filename,or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(status)
          status=nf_inq(ncid,ndims,nvars,natts,unlimid)
          call check_cdf_status(status)
        else
          if (ex) then
            write(6,'(a)') trim(filename)//' will be replaced'
          else
            write(6,'(a)') trim(filename)//' will be created'
          end if
#ifdef B25_EIRENE
          call create_cdf_start(ncid, filename,
     &     nCv, ns, nnatmi, nnmoli, ncvtot)
#else
          call create_cdf_start(ncid, filename,
     &     nCv, ns, nnatmi, nnmoli, nCi+1)
#endif
          if (write_b25_state_variables.ge.1) then
            call create_cdf_field3(ncid, 'na', 1,
     1       'Ion density', 'm^-3')
            call create_cdf_field2(ncid, 'ne',
     1       'Electron density', 'm^-3')
            call create_cdf_field2(ncid, 'te',
     1       'Electron temperature', 'eV')
            call create_cdf_field2(ncid, 'ti',
     1       'Ion temperature', 'eV')
            call create_cdf_field2(ncid, 'po',
     1       'Electric potential', 'V')
          endif
          if (write_b25_sources.ge.1) then
            call create_cdf_field3(ncid, 'rsana', 1,
     1       'Volumetric ionisation rate', 'm^-3 s^-1')
            call create_cdf_field3(ncid, 'rsahi', 1,
     1       'Volumetric ionisation energy source/sink', 'W m^-3')
            call create_cdf_field2(ncid, 'rsahisum',
     1       'Total volumetric ionisation energy source/sink',
     2       'W m^-3')
            call create_cdf_field3(ncid, 'rrana', 1,
     1       'Volumetric recombination rate', 'm^-3 s^-1')
            call create_cdf_field3(ncid, 'rrahi', 1,
     1       'Volumetric recombination energy source/sink',
     2       'W m^-3')
            call create_cdf_field2(ncid, 'rrahisum',
     1       'Total volumetric recombination energy source/sink',
     2       'W m^-3')
            call create_cdf_field3(ncid, 'rcxna', 1,
     1       'Volumetric charge exchange rate', 'm^-3 s^-1')
            call create_cdf_field3(ncid, 'rcxhi', 1,
     1       'Volumetric charge exchange energy source/sink',
     2       'W m^-3')
            call create_cdf_field2(ncid, 'rcxhisum',
     1       'Total volumetric charge exchange energy source/sink',
     2       'W m^-3')
            call create_cdf_field3(ncid, 'rqahe', 1,
     1       'Volumetric electron cooling rate', 'W m^-3')
            call create_cdf_field2(ncid, 'rqahesum',
     1       'Total volumetric electron cooling rate', 'W m^-3')
            call create_cdf_field3(ncid, 'rqrad', 1,
     1       'Volumetric line radiation rate', 'W m^-3')
            call create_cdf_field2(ncid, 'rqradsum',
     1       'Total volumetric line radiation rate', 'W m^-3')
          endif
#ifdef B25_EIRENE
          if (write_eirene_state_variables.ge.1) then
            if(nnatmi.gt.0) then
              call create_cdf_field3(ncid, 'dab2', 2,
     1         'Atom density', 'm^-3')
              call create_cdf_field3(ncid, 'tab2', 2,
     1         'Atom temperature', 'eV')
            endif
            if(nnmoli.gt.0) then
              call create_cdf_field3(ncid, 'dmb2', 3,
     1         'Molecular density', 'm^-3')
              call create_cdf_field3(ncid, 'tmb2', 3,
     1         'Molecule temperature', 'eV')
            endif
          endif
          if (write_eirene_sources.ge.1) then
            if(nnatmi.gt.0) then
              call create_cdf_field3(ncid, 'eirene_papl_sna', 1,
     1         'Particle source from atom-plasma interactions',
     2         'm^-3 s^-1')
              call create_cdf_field3(ncid, 'eirene_pppl_sna', 1,
     1         'Particle source from recombinations',
     2         'm^-3 s^-1')
              call create_cdf_field2(ncid, 'eirene_eael_she',
     1         'Electron energy source from atom-plasma interactions',
     2         'W m^-3')
              call create_cdf_field2(ncid, 'eirene_epel_she',
     1         'Electron energy source from recombinations',
     2         'W m^-3')
              call create_cdf_field2(ncid, 'eirene_eapl_shi',
     1         'Ion energy source from atom-plasma interactions',
     2         'W m^-3')
              call create_cdf_field2(ncid, 'eirene_eppl_shi',
     1         'Ion energy source from recombinations',
     2         'W m^-3')
              call create_cdf_field3(ncid, 'eneutrad', 2,
     1         'Volumetric radiation rate due to atoms', 'W m^-3')
              call create_cdf_field2(ncid, 'eneutradsum',
     1         'Total volumetric radiation rate due to atoms',
     2         'W m^-3')
            endif
            if(nnmoli.gt.0) then
              call create_cdf_field3(ncid, 'eirene_pmpl_sna', 1,
     1         'Particle source from molecule-plasma interactions',
     2         'm^-3 s^-1')
              call create_cdf_field2(ncid, 'eirene_emel_she',
     1         'Electron energy source from molecule-plasma ' //
     2         'interactions', 'W m^-3')
              call create_cdf_field2(ncid, 'eirene_empl_shi',
     1         'Ion energy source from molecule-plasma ' //
     2         'interactions', 'W m^-3')
              call create_cdf_field3(ncid, 'emolrad', 3,
     1         'Volumetric radiation rate due to molecules',
     2         'W m^-3')
             call create_cdf_field2(ncid, 'emolradsum',
     1         'Total volumetric radiation rate due to molecules',
     2         'W m^-3')
            endif
          endif
#endif
          call create_cdf_end(ncid)
        endif
      endif

#ifdef B25_EIRENE
      if (ncall.eq.0 .and. write_triangles.ge.1) then
        filename_tri='eirenemovies.nc'
        call find_file(filename_tri,ex)
        if (ex.and.switch%b2mndr_stim.lt.0.0_R8) then
          write(6,'(a)') trim(filename_tri)//' will be appended'
          status=nf_open(filename_tri,or(NF_WRITE,NF_SHARE),ncid)
          call check_cdf_status(status)
          status=nf_inq(ncid,ndims,nvars,natts,unlimid)
          call check_cdf_status(status)
        else
          if (ex) then
            write(6,'(a)') trim(filename_tri)//' will be replaced'
          else
            write(6,'(a)') trim(filename_tri)//' will be created'
          end if
          call create_cdf_start(ncid, filename_tri,
     &     nCv, ns, nnatmi, nnmoli, ncvtot-1)
          if (nnatmi.gt.0) then
            call create_cdf_field_tri(ncid, 'pdena', 1,
     1       'Atom particle density', 'm^-3')
            call create_cdf_field_tri(ncid, 'edena', 1,
     1       'Atom energy density', 'eV m^-3')
            call create_cdf_field_tri(ncid, 'vxdena', 1,
     1       'Atom x-directed momentum density', 'kg s^-1 m^-2')
            call create_cdf_field_tri(ncid, 'vydena', 1,
     1       'Atom y-directed momentum density', 'kg s^-1 m^-2')
            call create_cdf_field_tri(ncid, 'vzdena', 1,
     1       'Atom z-directed momentum density', 'kg s^-1 m^-2')
          endif
          if(nnmoli.gt.0) then
            call create_cdf_field_tri(ncid, 'pdenm', 2,
     1       'Molecular particle density', 'm^-3')
            call create_cdf_field_tri(ncid, 'edenm', 2,
     1       'Molecular energy density', 'eV m^-3')
            call create_cdf_field_tri(ncid, 'vxdenm', 2,
     1       'Molecular x-directed momentum density', 'kg s^-1 m^-2')
            call create_cdf_field_tri(ncid, 'vydenm', 2,
     1       'Molecular y-directed momentum density', 'kg s^-1 m^-2')
            call create_cdf_field_tri(ncid, 'vzdenm', 2,
     1       'Molecular z-directed momentum density', 'kg s^-1 m^-2')
          endif
          call create_cdf_end(ncid)
        endif
      endif
#endif

      time(1) = tim

      write(*,*) 'put_cdf_start'
      call put_cdf_start(ncid,filename,time)
      if (write_b25_state_variables.ge.1) then
        write(*,*) 'put_cdf_field3: na'
        call put_cdf_field3(ncid,'na',state%pl%na,nCv,ns)
        write(*,*) 'put_cdf_field2: ne'
        call put_cdf_field2(ncid,'ne',state%dv%ne,nCv)
        write(*,*) 'put_cdf_field2: te'
        dummy(:) = state%pl%te(:)/qe
        call put_cdf_field2(ncid,'te',dummy,nCv)
        write(*,*) 'put_cdf_field2: ti'
        dummy(:) = state%pl%ti(:)/qe
        call put_cdf_field2(ncid,'ti',dummy,nCv)
        write(*,*) 'put_cdf_field2: po'
        call put_cdf_field2(ncid,'po',state%pl%po,nCv)
      endif
      if (write_b25_sources.ge.1) then
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rsana(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rsana'
        call put_cdf_field3(ncid,'rsana',dummy3,nCv,ns)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rsahi(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rsahi'
        call put_cdf_field3(ncid,'rsahi',dummy3,nCv,ns)
        dummy = 0.0_R8
        do is=0,ns-1
          dummy(:) = dummy(:)+state%srw%rsahi(:,is)
        enddo
        dummy(:)=dummy(:)/geo%cvVol(:)
        write(*,*) 'put_cdf_field2: rsahisum'
        call put_cdf_field2(ncid,'rsahisum',dummy,nCv)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rcxna(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rcxna'
        call put_cdf_field3(ncid,'rcxna',dummy3,nCv,ns)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rrahi(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rrahi'
        call put_cdf_field3(ncid,'rrahi',dummy3,nCv,ns)
        dummy = 0.0_R8
        do is=0,ns-1
          dummy(:) = dummy(:)+state%srw%rrahi(:,is)
        enddo
        dummy(:)=dummy(:)/geo%cvVol(:)
        write(*,*) 'put_cdf_field2: rrahisum'
        call put_cdf_field2(ncid,'rrahisum',dummy,nCv)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rsana(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rsana'
        call put_cdf_field3(ncid,'rsana',dummy3,nCv,ns)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rcxhi(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rcxhi'
        call put_cdf_field3(ncid,'rcxhi',dummy3,nCv,ns)
        dummy = 0.0_R8
        do is=0,ns-1
          dummy(:)=dummy(:)+state%srw%rcxhi(:,is)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field2: rcxhisum'
        call put_cdf_field2(ncid,'rcxhisum',dummy,nCv)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rqahe(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rqahe'
        call put_cdf_field3(ncid,'rqahe',dummy3,nCv,ns)
        dummy = 0.0_R8
        do is=0,ns-1
          dummy(:) = dummy(:) + state%srw%rqahe(:,is)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field2: rqahesum'
        call put_cdf_field2(ncid,'rqahesum',dummy,nCv)
        dummy3 = 0.0_R8
        do is=1,ns
          dummy3(:,is) = state%srw%rqrad(:,is-1)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field3: rqrad'
        call put_cdf_field3(ncid,'rqrad',dummy3,nCv,ns)
        dummy = 0.0_R8
        do is=0,ns-1
          dummy(:) = dummy(:) + state%srw%rqrad(:,is)/geo%cvVol(:)
        enddo
        write(*,*) 'put_cdf_field2: rqradsum'
        call put_cdf_field2(ncid,'rqradsum',dummy,nCv)
      endif
#ifdef B25_EIRENE
      if (write_eirene_state_variables.ge.1) then
        if(nnatmi.gt.0) then
          dummy3atm = 0.0_R8
          dummy3atm(1:nCv,1:nnatmi) = dab2(1:nCv,1:nnatmi,1)
          write(*,*) 'put_cdf_field3: dab2'
          call put_cdf_field3(ncid,'dab2',dummy3atm,nCv,nnatmi)
          dummy3atm = 0.0_R8
          dummy3atm(1:nCv,1:nnatmi) = tab2(1:nCv,1:nnatmi,1)/qe
          write(*,*) 'put_cdf_field3: tab2'
          call put_cdf_field3(ncid,'tab2',dummy3atm,nCv,nnatmi)
        endif
        if(nnmoli.gt.0) then
          dummy3mol = 0.0_R8
          dummy3mol(1:nCv,1:nnmoli) = dmb2(1:nCv,1:nnmoli,1)
          write(*,*) 'put_cdf_field3: dmb2'
          call put_cdf_field3(ncid,'dmb2',dummy3mol,nCv,nnmoli)
          dummy3mol = 0.0_R8
          dummy3mol(1:nCv,1:nnmoli) = tmb2(1:nCv,1:nnmoli,1)/qe
          write(*,*) 'put_cdf_field3: tmb2'
          call put_cdf_field3(ncid,'tmb2',dummy3mol,nCv,nnmoli)
        endif
      endif
      if (write_eirene_sources.ge.1) then
        if(nnatmi.gt.0) then
          nstrad=size(eirene_mc_papl_sna_bal,3)
          dummy3 = 0.0_R8
          do istra = 1, nstrad
            dummy3(1:nCv,1:ns) = dummy3(1:nCv,1:ns) +
     .       eirene_mc_papl_sna_bal(1:nCv,0:ns-1,istra)
          enddo
          do is=1,ns
            dummy3(1:nCv,is) = dummy3(1:nCv,is)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field3: eirene_papl_sna'
          call put_cdf_field3(ncid,'eirene_papl_sna',dummy3,nCv,ns)
          dummy3 = 0.0_R8
          do istra = 1, nstrad
            dummy3(1:nCv,1:ns) = dummy3(1:nCv,1:ns) +
     .       eirene_mc_pppl_sna_bal(1:nCv,0:ns-1,istra)
          enddo
          do is=1,ns
            dummy3(1:nCv,is) = dummy3(1:nCv,is)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field3: eirene_pppl_sna'
          call put_cdf_field3(ncid,'eirene_pppl_sna',dummy3,nCv,ns)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_eael_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_eael_she'
          call put_cdf_field2(ncid,'eirene_eael_she',dummy,nCv)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_epel_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_epel_she'
          call put_cdf_field2(ncid,'eirene_epel_she',dummy,nCv)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_eapl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_eapl_shi'
          call put_cdf_field2(ncid,'eirene_eapl_shi',dummy,nCv)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_eppl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_eppl_shi'
          call put_cdf_field2(ncid,'eirene_eppl_shi',dummy,nCv)
          dummy3atm = 0.0_R8
          do iatmi = 1, nnatmi
            dummy3atm(1:nCv,iatmi) =
     .       eneutrad(1:nCv,iatmi,0)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field3: eneutrad'
          call put_cdf_field3(ncid,'eneutrad',dummy3atm,nCv,
     1     nnatmi)
          dummy = 0.0_R8
          do iatmi = 1, nnatmi
            dummy(:) = dummy(:) + eneutrad(:,iatmi,0)/geo%cvVol(:)
          enddo
          write(*,*) 'put_cdf_field2: eneutradsum'
          call put_cdf_field2(ncid,'eneutradsum',dummy,nCv)
        endif
        if(nnmoli.gt.0) then
          nstrad=size(eirene_mc_pmpl_sna_bal,3)
          dummy3 = 0.0_R8
          do istra = 1, nstrad
            dummy3(1:nCv,1:ns) = dummy3(1:nCv,1:ns) +
     .       eirene_mc_pmpl_sna_bal(1:nCv,0:ns-1,istra)
          enddo
          do is=1,ns
            dummy3(:,is) = dummy3(:,is)/geo%cvVol(:)
          enddo
          write(*,*) 'put_cdf_field3: eirene_pmpl_sna'
          call put_cdf_field3(ncid,'eirene_pmpl_sna',dummy3,nCv,ns)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_emel_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_emel_she'
          call put_cdf_field2(ncid,'eirene_emel_she',dummy3,nCv)
          dummy = 0.0_R8
          do istra = 1, nstrad
            dummy(1:nCv) = dummy(1:nCv) +
     .       eirene_mc_empl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
          enddo
          write(*,*) 'put_cdf_field2: eirene_empl_shi'
          call put_cdf_field2(ncid,'eirene_empl_shi',dummy3,nCv)
          dummy3mol = 0.0_R8
          do imoli = 1, nnmoli
            dummy3mol(:,imoli) = emolrad(:,imoli,0)/geo%cvVol(:)
          enddo
          write(*,*) 'put_cdf_field3: emolrad'
          call put_cdf_field3(ncid,'emolrad',dummy3mol,nCv,nnmoli)
          dummy = 0.0_R8
          do imoli = 1, nnmoli
            dummy(:) = dummy(:) + emolrad(:,imoli,0)/geo%cvVol(:)
          enddo
          write(*,*) 'put_cdf_field2: emolradsum'
          call put_cdf_field2(ncid,'emolradsum',dummy,nCv)
        endif
      endif
#endif
      write(*,*) 'put_cdf_end'
      call put_cdf_end(ncid)

#ifdef B25_EIRENE
      if (write_triangles.ge.1) then
        write(*,*) 'put_cdf_start'
        call put_cdf_start(ncid,filename_tri,time)
        if (nnatmi.gt.0) then
          dummytriatm = 0.0_R8
          do iatmi = 1, nnatmi
            dummytriatm(1:ncvtot-1,iatmi) =
     .       pdena(iatmi,1:ncvtot-1)*1e6_R8
          enddo
          write(*,*) 'put_cdf_field_tri: pdena'
          call put_cdf_field_tri
     .     (ncid,'pdena',dummytriatm,ncvtot-1,nnatmi)
          dummytriatm = 0.0_R8
          do iatmi = 1, nnatmi
            dummytriatm(1:ncvtot-1,iatmi) =
     .       edena(iatmi,1:ncvtot-1)*1e6_R8
          enddo
          write(*,*) 'put_cdf_field_tri: edena'
          call put_cdf_field_tri
     .     (ncid,'edena',dummytriatm,ncvtot-1,nnatmi)
          dummytriatm = 0.0_R8
          do iatmi = 1, nnatmi
            dummytriatm(1:ncvtot-1,iatmi) = vxdena(iatmi,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vxdena'
          call put_cdf_field_tri
     .     (ncid,'vxdena',dummytriatm,ncvtot-1,nnatmi)
          dummytriatm = 0.0_R8
          do iatmi = 1, nnatmi
            dummytriatm(1:ncvtot-1,iatmi) = vydena(iatmi,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vydena'
          call put_cdf_field_tri
     .     (ncid,'vydena',dummytriatm,ncvtot-1,nnatmi)
          dummytriatm = 0.0_R8
          do iatmi = 1, nnatmi
            dummytriatm(1:ncvtot-1,iatmi) = vzdena(iatmi,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vzdena'
          call put_cdf_field_tri
     .     (ncid,'vzdena',dummytriatm,ncvtot-1,nnatmi)
        endif
        if(nnmoli.gt.0) then
          dummytrimol = 0.0_R8
          do imoli = 1, nnmoli
            dummytrimol(1:ncvtot-1,imoli) =
     .       pdenm(imoli,1:ncvtot-1)*1e6_R8
          enddo
          write(*,*) 'put_cdf_field_tri: pdenm'
          call put_cdf_field_tri
     .     (ncid,'pdenm',dummytrimol,ncvtot-1,nnmoli)
          dummytrimol = 0.0_R8
          do imoli = 1, nnmoli
            dummytrimol(1:ncvtot-1,imoli) =
     .       edenm(imoli,1:ncvtot-1)*1e6_R8
          enddo
          write(*,*) 'put_cdf_field_tri: edenm'
          call put_cdf_field_tri
     .     (ncid,'edenm',dummytrimol,ncvtot-1,nnmoli)
          dummytrimol = 0.0_R8
          do imoli = 1, nnmoli
            dummytrimol(1:ncvtot-1,imoli) = vxdenm(imoli,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vxdenm'
          call put_cdf_field_tri
     .     (ncid,'vxdenm',dummytrimol,ncvtot-1,nnmoli)
          dummytrimol = 0.0_R8
          do imoli = 1, nnmoli
            dummytrimol(1:ncvtot-1,imoli) = vydenm(imoli,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vydenm'
          call put_cdf_field_tri
     .     (ncid,'vydenm',dummytrimol,ncvtot-1,nnmoli)
          dummytrimol = 0.0_R8
          do imoli = 1, nnmoli
            dummytrimol(1:ncvtot-1,imoli) = vzdenm(imoli,1:ncvtot-1)
          enddo
          write(*,*) 'put_cdf_field_tri: vzdenm'
          call put_cdf_field_tri
     .     (ncid,'vzdenm',dummytrimol,ncvtot-1,nnmoli)
        endif
        write(*,*) 'put_cdf_end'
        call put_cdf_end(ncid)
      endif
#endif

* ..return
      ncall = ncall+1
      call subend ()
#endif
      return
      end subroutine cdfmovie
