      subroutine cdfchord
#ifndef NO_CDF
     & (ncid,nCv,nFc,ns,geo,mpg,state,switch)
      use b2mod_types
      use b2mod_time
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_dimensions
      use b2mod_chords
      use b2mod_chord_interp
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_user_namelist
     & , only : rzchords, nchords, nchordmax
      use b2mod_switches
      use b2mod_subsys
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_extrab25
      use eirmod_cestim
      use b2mod_balance
      use eirmod_eirbra
      use eirmod_ctrig
     & , only : NTRII
      use b2mod_chord_interp_tri
     & , only : build_tri_field, sample_tri_at_points, mask_inside_tri
#endif
      implicit none
      integer ncid, nCv, nFc, ns
      type (geometry), intent(in) :: geo
      type (mapping), intent(in) :: mpg
      type (B2State), intent(in) :: state
      type (switches), intent(in) :: switch
#     include <netcdf.inc>
      integer is, ncall
      integer ni_lb, m, ioff, grp
c   ..b2chords.nc geometry: sample points kept inside the B2.5 grid, one
c   ..netCDF-4 group per surviving chord
      integer, save :: npt = 0, nchord_b = 0
      logical, save :: have_b25 = .false.
      real (kind=R8), allocatable, save :: crx(:), cry(:), ds(:)
      integer,        allocatable, save :: ichord(:), chordid(:)
      real (kind=R8), allocatable, save :: cr1(:), cz1(:)
      real (kind=R8), allocatable, save :: cr2(:), cz2(:)
      integer,        allocatable, save :: off(:), cnt(:), grpids(:)
c   ..full (unfiltered) sample set for all chords, built once
      integer :: npt0
      real (kind=R8), allocatable :: crx0(:), cry0(:), ds0(:)
      integer,        allocatable :: ichord0(:)
      logical,        allocatable :: inside(:)
      real (kind=R8), allocatable :: valg(:), valg3(:,:), valg3b(:,:)
      real (kind=R8), allocatable :: wrk(:), wrk3(:,:), cpol(:), crad(:)
      real (kind=R8) :: time(1), dl, res_old
      real (kind=R8) :: mrs(nchordmax), mzs(nchordmax), mre(nchordmax),
     & mze(nchordmax)
      integer :: cid_o(nchordmax), npt_o(nchordmax), nchord_old
      integer :: natm_b, nmol_b
      logical :: ex, okmeta, samechords
      character*256, save :: filename
#ifdef B25_EIRENE
c   ..eirenechords.nc geometry: points kept inside the triangular grid
      integer, save :: ncid_tri
      integer, save :: npt_t = 0, nchord_t = 0
      logical, save :: have_tri = .false.
      real (kind=R8), allocatable, save :: crx_t(:), cry_t(:), ds_t(:)
      integer,        allocatable, save :: ichord_t(:), chordid_t(:)
      real (kind=R8), allocatable, save :: cr1_t(:), cz1_t(:)
      real (kind=R8), allocatable, save :: cr2_t(:), cz2_t(:)
      integer, allocatable, save :: off_t(:), cnt_t(:), grpids_t(:)
      character*256, save :: filename_tri
      integer :: iatmi, imoli, istra, nstrad
      logical :: ex_tri, samechords_tri
      real (kind=R8), allocatable :: ftri(:), valgtri(:,:)
      real (kind=R8), allocatable :: valgatm(:,:), valgmol(:,:)
      real (kind=R8), allocatable :: wrkatm(:,:), wrkmol(:,:)
#endif
      save ncall
      data ncall/0/

      call subini ('cdfchord')
      dl = switch % cdfchord_resolution
      natm_b = 0
      nmol_b = 0
#ifdef B25_EIRENE
      natm_b = nnatmi
      nmol_b = nnmoli
#endif

      if (nchords .lt. 1) then
        if (ncall.eq.0) write(6,'(a)')
     &    'cdfchord: no chords defined in rzchords; b2chords.nc skipped'
        ncall = ncall + 1
        call subend ()
        return
      endif

      if (ncall.eq.0) then
c   ..full concatenated sample set for all chords (grid-independent)
        call chords_total_points(rzchords, nchords, dl, npt0)
        allocate(crx0(npt0), cry0(npt0), ds0(npt0), ichord0(npt0))
        allocate(inside(npt0))
        call chords_geometry(rzchords, nchords, dl, npt0,
     &    crx0, cry0, ds0, ichord0)
c   ..b2chords.nc: keep only points inside the B2.5 grid, dropping any
c   ..chord that lies entirely outside it
        filename = 'b2chords.nc'
        call mask_inside_cv(geo, mpg, crx0, cry0, ichord0, npt0, inside)
        call compact_chords(crx0, cry0, ds0, ichord0, npt0, nchords,
     &    inside, crx, cry, ds, ichord, npt, chordid, nchord_b,
     &    cr1, cz1, cr2, cz2)
        have_b25 = (nchord_b .gt. 0)
        if (have_b25) then
c   ..per-chord offsets/counts in the concatenated arrays
          allocate(off(nchord_b), cnt(nchord_b), grpids(nchord_b))
          ioff = 0
          do m = 1, nchord_b
            off(m) = ioff
            cnt(m) = count(ichord(1:npt).eq.m)
            ioff = ioff + cnt(m)
          enddo
c   ..decide between append (same chords) and replace/create
          call find_file(filename, ex)
          samechords = .false.
          if (ex .and. switch%b2mndr_stim.lt.0.0_R8) then
            call get_cdf_chord_meta(filename, nchordmax, okmeta,
     &        nchord_old, cid_o, npt_o, mrs, mzs, mre, mze, res_old)
            if (okmeta .and. nchord_old.eq.nchord_b .and.
     &          abs(res_old-dl).le.1.0e-12_R8) then
              samechords = .true.
              do m = 1, nchord_b
                if (cid_o(m).ne.chordid(m) .or. npt_o(m).ne.cnt(m) .or.
     &            abs(mrs(m)-cr1(m)).gt.1.0e-9_R8 .or.
     &            abs(mzs(m)-cz1(m)).gt.1.0e-9_R8 .or.
     &            abs(mre(m)-cr2(m)).gt.1.0e-9_R8 .or.
     &            abs(mze(m)-cz2(m)).gt.1.0e-9_R8) samechords = .false.
              enddo
            endif
          endif
          if (samechords) then
            write(6,'(a)') trim(filename)//' will be appended'
          else
            if (ex) then
              write(6,'(a)') trim(filename)//' will be replaced'
            else
              write(6,'(a)') trim(filename)//' will be created'
            endif
            call create_cdf_chord_start(ncid, filename, ns, natm_b,
     &        nmol_b)
            do m = 1, nchord_b
              call create_cdf_chord_group(ncid, chordid(m), cnt(m), grp)
              call define_b25_fields(grp)
            enddo
            call create_cdf_chord_end(ncid, nchord_b, dl)
            call put_cdf_chord_coords(filename, nchord_b, chordid, off,
     &        cnt, npt, crx, cry, ds, cr1, cz1, cr2, cz2)
          endif
        else
          write(6,'(a)') 'cdfchord: no chord intersects the B2.5 grid;'
     &      //' b2chords.nc skipped'
        endif
#ifdef B25_EIRENE
c   ..eirenechords.nc: EIRENE triangle quantities along the chords,
c   ..keeping only the points inside the triangular grid
        if (switch%write_chord_triangles.ge.1) then
          filename_tri = 'eirenechords.nc'
          call mask_inside_tri(crx0, cry0, ichord0, npt0, inside)
          call compact_chords(crx0, cry0, ds0, ichord0, npt0, nchords,
     &      inside, crx_t, cry_t, ds_t, ichord_t, npt_t, chordid_t,
     &      nchord_t, cr1_t, cz1_t, cr2_t, cz2_t)
          have_tri = (nchord_t .gt. 0)
          if (have_tri) then
            allocate(off_t(nchord_t), cnt_t(nchord_t))
            allocate(grpids_t(nchord_t))
            ioff = 0
            do m = 1, nchord_t
              off_t(m) = ioff
              cnt_t(m) = count(ichord_t(1:npt_t).eq.m)
              ioff = ioff + cnt_t(m)
            enddo
            call find_file(filename_tri, ex_tri)
            samechords_tri = .false.
            if (ex_tri .and. switch%b2mndr_stim.lt.0.0_R8) then
              call get_cdf_chord_meta(filename_tri, nchordmax, okmeta,
     &          nchord_old, cid_o, npt_o, mrs, mzs, mre, mze, res_old)
              if (okmeta .and. nchord_old.eq.nchord_t .and.
     &            abs(res_old-dl).le.1.0e-12_R8) then
                samechords_tri = .true.
                do m = 1, nchord_t
                  if (cid_o(m).ne.chordid_t(m) .or.
     &                npt_o(m).ne.cnt_t(m) .or.
     &                abs(mrs(m)-cr1_t(m)).gt.1.0e-9_R8 .or.
     &                abs(mzs(m)-cz1_t(m)).gt.1.0e-9_R8 .or.
     &                abs(mre(m)-cr2_t(m)).gt.1.0e-9_R8 .or.
     &                abs(mze(m)-cz2_t(m)).gt.1.0e-9_R8)
     &                samechords_tri = .false.
                enddo
              endif
            endif
            if (samechords_tri) then
              write(6,'(a)') trim(filename_tri)//' will be appended'
            else
              if (ex_tri) then
                write(6,'(a)') trim(filename_tri)//' will be replaced'
              else
                write(6,'(a)') trim(filename_tri)//' will be created'
              endif
              call create_cdf_chord_start(ncid_tri, filename_tri, ns,
     &          nnatmi, nnmoli)
              do m = 1, nchord_t
                call create_cdf_chord_group(ncid_tri, chordid_t(m),
     &            cnt_t(m), grp)
                call define_tri_fields(grp)
              enddo
              call create_cdf_chord_end(ncid_tri, nchord_t, dl)
              call put_cdf_chord_coords(filename_tri, nchord_t,
     &          chordid_t, off_t, cnt_t, npt_t, crx_t, cry_t, ds_t,
     &          cr1_t, cz1_t, cr2_t, cz2_t)
            endif
          else
            write(6,'(a)')
     &        'cdfchord: no chord intersects the triangular grid;'
     &        //' eirenechords.nc skipped'
          endif
        endif
#endif
        deallocate(crx0, cry0, ds0, ichord0, inside)
      endif

      time(1) = tim

c   ..b2chords.nc frame
      if (have_b25) then
        allocate(valg(npt), valg3(npt,ns), valg3b(npt,ns))
        allocate(wrk(nCv), wrk3(nCv,ns), cpol(nCv), crad(nCv))
#ifdef B25_EIRENE
        allocate(wrkatm(nCv,max(natm_b,1)), wrkmol(nCv,max(nmol_b,1)))
        allocate(valgatm(npt,max(natm_b,1)), valgmol(npt,max(nmol_b,1)))
#endif
        call put_cdf_chord_start(ncid, filename, time)
        do m = 1, nchord_b
          call get_cdf_chord_group(ncid, chordid(m), grpids(m))
        enddo

        if (switch%write_chord_b25_state_variables.ge.1) then
          call put_chord_field3('na', state%pl%na)
          call put_chord_field2('ne', state%dv%ne)
          call put_chord_field3('ua', state%pl%ua)
          wrk(:) = state%pl%te(:)/qe
          call put_chord_field2('te', wrk)
          wrk(:) = state%pl%ti(:)/qe
          call put_chord_field2('ti', wrk)
          call put_chord_field2('po', state%pl%po)
        endif

        if (switch%write_chord_b25_fluxes.ge.1) then
c   ..particle flux density, per species (poloidal fnax, radial fnay)
          do is = 1, ns
            call facefield_to_cellcenters(geo, mpg,
     &        state%dv%fna(:,:,is-1), cpol, crad)
            call sample_cv_at_points(geo,mpg,cpol,crx,cry,ichord,npt,
     &        valg3(:,is))
            call sample_cv_at_points(geo,mpg,crad,crx,cry,ichord,npt,
     &        valg3b(:,is))
          enddo
          do m = 1, nchord_b
            call put_cdf_chord_field3(grpids(m),'fnax',
     &        valg3(off(m)+1:off(m)+cnt(m),:), cnt(m), ns)
            call put_cdf_chord_field3(grpids(m),'fnay',
     &        valg3b(off(m)+1:off(m)+cnt(m),:), cnt(m), ns)
          enddo
c   ..electron heat flux density
          call facefield_to_cellcenters(geo,mpg,state%dv%fhe,cpol,crad)
          call put_chord_field2('fhex', cpol)
          call put_chord_field2('fhey', crad)
c   ..ion heat flux density
          call facefield_to_cellcenters(geo,mpg,state%dv%fhi,cpol,crad)
          call put_chord_field2('fhix', cpol)
          call put_chord_field2('fhiy', crad)
c   ..electric current density
          call facefield_to_cellcenters(geo,mpg,state%dv%fch,cpol,crad)
          call put_chord_field2('fchx', cpol)
          call put_chord_field2('fchy', crad)
        endif

        if (switch%write_chord_b25_sources.ge.1) then
          do is=1,ns
            wrk3(:,is) = state%srw%rsana(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rsana', wrk3)
          do is=1,ns
            wrk3(:,is) = state%srw%rsahi(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rsahi', wrk3)
          wrk = 0.0_R8
          do is=0,ns-1
            wrk(:) = wrk(:)+state%srw%rsahi(:,is)
          enddo
          wrk(:) = wrk(:)/geo%cvVol(:)
          call put_chord_field2('rsahisum', wrk)
          do is=1,ns
            wrk3(:,is) = state%srw%rcxna(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rcxna', wrk3)
          do is=1,ns
            wrk3(:,is) = state%srw%rrahi(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rrahi', wrk3)
          wrk = 0.0_R8
          do is=0,ns-1
            wrk(:) = wrk(:)+state%srw%rrahi(:,is)
          enddo
          wrk(:) = wrk(:)/geo%cvVol(:)
          call put_chord_field2('rrahisum', wrk)
          do is=1,ns
            wrk3(:,is) = state%srw%rsana(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rsana', wrk3)
          do is=1,ns
            wrk3(:,is) = state%srw%rcxhi(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rcxhi', wrk3)
          wrk = 0.0_R8
          do is=0,ns-1
            wrk(:) = wrk(:)+state%srw%rcxhi(:,is)/geo%cvVol(:)
          enddo
          call put_chord_field2('rcxhisum', wrk)
          do is=1,ns
            wrk3(:,is) = state%srw%rqahe(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rqahe', wrk3)
          wrk = 0.0_R8
          do is=0,ns-1
            wrk(:) = wrk(:) + state%srw%rqahe(:,is)/geo%cvVol(:)
          enddo
          call put_chord_field2('rqahesum', wrk)
          do is=1,ns
            wrk3(:,is) = state%srw%rqrad(:,is-1)/geo%cvVol(:)
          enddo
          call put_chord_field3('rqrad', wrk3)
          wrk = 0.0_R8
          do is=0,ns-1
            wrk(:) = wrk(:) + state%srw%rqrad(:,is)/geo%cvVol(:)
          enddo
          call put_chord_field2('rqradsum', wrk)
        endif

        if (switch%write_chord_b25_transport.ge.1) then
          call put_chord_field3('dna', state%co%dna0)
          do is=1,ns
            if (switch%tn_style.eq.0) then
              wrk3(:,is) = state%co%dpa0(:,is-1)*(state%rt%rza(:,is-1)*
     &              state%pl%te(:)+state%pl%ti(:))
            else
              wrk3(:,is) = state%co%dpa0(:,is-1)*state%pl%tn(:)
            endif
          enddo
          call put_chord_field3('dpa', wrk3)
          ni_lb = lbound(state%dv%ni,2)
          if (switch%tn_style.eq.0) then
            wrk(:) = state%co%hci0(:)/ state%dv%ni(:,ni_lb)
          else if (switch%tn_style.eq.1) then
            wrk(:) = state%co%hci0(:)/state%dv%ni(:,ni_lb+1)
          else if (switch%tn_style.eq.2) then
            wrk(:) = state%co%hci0(:) /
     &                  (state%dv%ni(:,ni_lb) - state%dv%nn(:))
          end if
          call put_chord_field2('hci', wrk)
          wrk(:) = state%co%hce0(:)/state%dv%ne(:)
          call put_chord_field2('hce', wrk)
          do is=1,ns
            wrk3(:,is) = state%co%vla0(:,0,is-1)
          enddo
          call put_chord_field3('vlax', wrk3)
          do is=1,ns
            wrk3(:,is) = state%co%vla0(:,1,is-1)
          enddo
          call put_chord_field3('vlay', wrk3)
          do is=1,ns
            wrk3(:,is) = state%co%vsa0(:,is-1)/(mp*am(is-1)*
     &          state%pl%na(:,is-1))
          enddo
          call put_chord_field3('vsa', wrk3)
        endif

#ifdef B25_EIRENE
        if (switch%write_chord_eirene_state_variables.ge.1) then
          if (nnatmi.gt.0) then
            wrkatm(1:nCv,1:nnatmi) = dab2(1:nCv,1:nnatmi,1)
            call put_chord_fieldn('dab2', wrkatm, nnatmi, valgatm)
            wrkatm(1:nCv,1:nnatmi) = tab2(1:nCv,1:nnatmi,1)/qe
            call put_chord_fieldn('tab2', wrkatm, nnatmi, valgatm)
          endif
          if (nnmoli.gt.0) then
            wrkmol(1:nCv,1:nnmoli) = dmb2(1:nCv,1:nnmoli,1)
            call put_chord_fieldn('dmb2', wrkmol, nnmoli, valgmol)
            wrkmol(1:nCv,1:nnmoli) = tmb2(1:nCv,1:nnmoli,1)/qe
            call put_chord_fieldn('tmb2', wrkmol, nnmoli, valgmol)
          endif
        endif
        if (switch%write_chord_eirene_fluxes.ge.1) then
          if (nnatmi.gt.0) then
            wrkatm(1:nCv,1:nnatmi) = rfluxa(1:nCv,1:nnatmi,1)
            call put_chord_fieldn('rfluxa', wrkatm, nnatmi, valgatm)
            wrkatm(1:nCv,1:nnatmi) = pfluxa(1:nCv,1:nnatmi,1)
            call put_chord_fieldn('pfluxa', wrkatm, nnatmi, valgatm)
            wrkatm(1:nCv,1:nnatmi) = refluxa(1:nCv,1:nnatmi,1)
            call put_chord_fieldn('refluxa', wrkatm, nnatmi, valgatm)
            wrkatm(1:nCv,1:nnatmi) = pefluxa(1:nCv,1:nnatmi,1)
            call put_chord_fieldn('pefluxa', wrkatm, nnatmi, valgatm)
          endif
          if (nnmoli.gt.0) then
            wrkmol(1:nCv,1:nnmoli) = rfluxm(1:nCv,1:nnmoli,1)
            call put_chord_fieldn('rfluxm', wrkmol, nnmoli, valgmol)
            wrkmol(1:nCv,1:nnmoli) = pfluxm(1:nCv,1:nnmoli,1)
            call put_chord_fieldn('pfluxm', wrkmol, nnmoli, valgmol)
            wrkmol(1:nCv,1:nnmoli) = refluxm(1:nCv,1:nnmoli,1)
            call put_chord_fieldn('refluxm', wrkmol, nnmoli, valgmol)
            wrkmol(1:nCv,1:nnmoli) = pefluxm(1:nCv,1:nnmoli,1)
            call put_chord_fieldn('pefluxm', wrkmol, nnmoli, valgmol)
          endif
        endif
        if (switch%write_chord_eirene_sources.ge.1) then
          if (nnatmi.gt.0) then
            nstrad = size(eirene_mc_papl_sna_bal,3)
            wrk3 = 0.0_R8
            do istra = 1, nstrad
              wrk3(1:nCv,1:ns) = wrk3(1:nCv,1:ns) +
     &         eirene_mc_papl_sna_bal(1:nCv,0:ns-1,istra)
            enddo
            do is = 1, ns
              wrk3(1:nCv,is) = wrk3(1:nCv,is)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field3('eirene_papl_sna', wrk3)
            wrk3 = 0.0_R8
            do istra = 1, nstrad
              wrk3(1:nCv,1:ns) = wrk3(1:nCv,1:ns) +
     &         eirene_mc_pppl_sna_bal(1:nCv,0:ns-1,istra)
            enddo
            do is = 1, ns
              wrk3(1:nCv,is) = wrk3(1:nCv,is)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field3('eirene_pppl_sna', wrk3)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_eael_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_eael_she', wrk)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_epel_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_epel_she', wrk)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_eapl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_eapl_shi', wrk)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_eppl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_eppl_shi', wrk)
            do iatmi = 1, nnatmi
              wrkatm(1:nCv,iatmi) =
     &         eneutrad(1:nCv,iatmi,0)/geo%cvVol(1:nCv)
            enddo
            call put_chord_fieldn('eneutrad', wrkatm, nnatmi, valgatm)
            wrk = 0.0_R8
            do iatmi = 1, nnatmi
              wrk(:) = wrk(:) + eneutrad(:,iatmi,0)/geo%cvVol(:)
            enddo
            call put_chord_field2('eneutradsum', wrk)
          endif
          if (nnmoli.gt.0) then
            nstrad = size(eirene_mc_pmpl_sna_bal,3)
            wrk3 = 0.0_R8
            do istra = 1, nstrad
              wrk3(1:nCv,1:ns) = wrk3(1:nCv,1:ns) +
     &         eirene_mc_pmpl_sna_bal(1:nCv,0:ns-1,istra)
            enddo
            do is = 1, ns
              wrk3(1:nCv,is) = wrk3(1:nCv,is)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field3('eirene_pmpl_sna', wrk3)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_emel_she_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_emel_she', wrk)
            wrk = 0.0_R8
            do istra = 1, nstrad
              wrk(1:nCv) = wrk(1:nCv) +
     &         eirene_mc_empl_shi_bal(1:nCv,istra)/geo%cvVol(1:nCv)
            enddo
            call put_chord_field2('eirene_empl_shi', wrk)
            do imoli = 1, nnmoli
              wrkmol(1:nCv,imoli) =
     &         emolrad(1:nCv,imoli,0)/geo%cvVol(1:nCv)
            enddo
            call put_chord_fieldn('emolrad', wrkmol, nnmoli, valgmol)
            wrk = 0.0_R8
            do imoli = 1, nnmoli
              wrk(:) = wrk(:) + emolrad(:,imoli,0)/geo%cvVol(:)
            enddo
            call put_chord_field2('emolradsum', wrk)
          endif
        endif
#endif

        call put_cdf_chord_end(ncid)
        deallocate(valg, valg3, valg3b, wrk, wrk3, cpol, crad)
#ifdef B25_EIRENE
        deallocate(wrkatm, wrkmol, valgatm, valgmol)
#endif
      endif

#ifdef B25_EIRENE
c   ..eirenechords.nc frame (EIRENE triangle quantities along the chords)
      if (have_tri) then
        allocate(ftri(NTRII), valgtri(npt_t, max(nnatmi,nnmoli,1)))
        call put_cdf_chord_start(ncid_tri, filename_tri, time)
        do m = 1, nchord_t
          call get_cdf_chord_group(ncid_tri, chordid_t(m), grpids_t(m))
        enddo
        if (nnatmi.gt.0) then
          call put_chord_tri('pdena', pdena, nnatmi, 1.0e6_R8)
          call put_chord_tri('edena', edena, nnatmi, 1.0e6_R8)
          call put_chord_tri('vxdena', vxdena, nnatmi, 1.0_R8)
          call put_chord_tri('vydena', vydena, nnatmi, 1.0_R8)
          call put_chord_tri('vzdena', vzdena, nnatmi, 1.0_R8)
        endif
        if (nnmoli.gt.0) then
          call put_chord_tri('pdenm', pdenm, nnmoli, 1.0e6_R8)
          call put_chord_tri('edenm', edenm, nnmoli, 1.0e6_R8)
          call put_chord_tri('vxdenm', vxdenm, nnmoli, 1.0_R8)
          call put_chord_tri('vydenm', vydenm, nnmoli, 1.0_R8)
          call put_chord_tri('vzdenm', vzdenm, nnmoli, 1.0_R8)
        endif
        call put_cdf_chord_end(ncid_tri)
        deallocate(ftri, valgtri)
      endif
#endif

      ncall = ncall+1
      call subend ()
#endif
      return

#ifndef NO_CDF
      contains

      subroutine define_b25_fields(g)
c     Define all b2chords fields (selected by the switches) in group g.
      implicit none
      integer g
      if (switch%write_chord_b25_state_variables.ge.1) then
        call create_cdf_chord_field3(g, 'na', 1,
     1   'Ion density', 'm^-3')
        call create_cdf_chord_field2(g, 'ne',
     1   'Electron density', 'm^-3')
        call create_cdf_chord_field3(g, 'ua', 1,
     1   'Ion parallel velocity', 'm s^-1')
        call create_cdf_chord_field2(g, 'te',
     1   'Electron temperature', 'eV')
        call create_cdf_chord_field2(g, 'ti',
     1   'Ion temperature', 'eV')
        call create_cdf_chord_field2(g, 'po',
     1   'Electric potential', 'V')
      endif
      if (switch%write_chord_b25_fluxes.ge.1) then
        call create_cdf_chord_field3(g, 'fnax', 1,
     1   'Poloidal particle flux density', 'm^-2 s^-1')
        call create_cdf_chord_field3(g, 'fnay', 1,
     1   'Radial particle flux density', 'm^-2 s^-1')
        call create_cdf_chord_field2(g, 'fhex',
     1   'Poloidal electron heat flux density', 'W m^-2')
        call create_cdf_chord_field2(g, 'fhey',
     1   'Radial electron heat flux density', 'W m^-2')
        call create_cdf_chord_field2(g, 'fhix',
     1   'Poloidal ion heat flux density', 'W m^-2')
        call create_cdf_chord_field2(g, 'fhiy',
     1   'Radial ion heat flux density', 'W m^-2')
        call create_cdf_chord_field2(g, 'fchx',
     1   'Poloidal electric current density', 'A m^-2')
        call create_cdf_chord_field2(g, 'fchy',
     1   'Radial electric current density', 'A m^-2')
      endif
      if (switch%write_chord_b25_sources.ge.1) then
        call create_cdf_chord_field3(g, 'rsana', 1,
     1   'Volumetric ionisation rate', 'm^-3 s^-1')
        call create_cdf_chord_field3(g, 'rsahi', 1,
     1   'Volumetric ionisation energy source/sink', 'W m^-3')
        call create_cdf_chord_field2(g, 'rsahisum',
     1   'Total volumetric ionisation energy source/sink', 'W m^-3')
        call create_cdf_chord_field3(g, 'rrana', 1,
     1   'Volumetric recombination rate', 'm^-3 s^-1')
        call create_cdf_chord_field3(g, 'rrahi', 1,
     1   'Volumetric recombination energy source/sink', 'W m^-3')
        call create_cdf_chord_field2(g, 'rrahisum',
     1   'Total volumetric recombination energy source/sink', 'W m^-3')
        call create_cdf_chord_field3(g, 'rcxna', 1,
     1   'Volumetric charge exchange rate', 'm^-3 s^-1')
        call create_cdf_chord_field3(g, 'rcxhi', 1,
     1   'Volumetric charge exchange energy source/sink', 'W m^-3')
        call create_cdf_chord_field2(g, 'rcxhisum',
     1   'Total volumetric charge exchange energy source/sink',
     2   'W m^-3')
        call create_cdf_chord_field3(g, 'rqahe', 1,
     1   'Volumetric electron cooling rate', 'W m^-3')
        call create_cdf_chord_field2(g, 'rqahesum',
     1   'Total volumetric electron cooling rate', 'W m^-3')
        call create_cdf_chord_field3(g, 'rqrad', 1,
     1   'Volumetric line radiation rate', 'W m^-3')
        call create_cdf_chord_field2(g, 'rqradsum',
     1   'Total volumetric line radiation rate', 'W m^-3')
      endif
      if (switch%write_chord_b25_transport.ge.1) then
        call create_cdf_chord_field3(g, 'dna', 1,
     1   'Particle density-driven diffusivity', 'm^-2 s^-1')
        call create_cdf_chord_field3(g, 'dpa', 1,
     1   'Particle pressure-driven diffusivity', 'm^-2 s^-1')
        call create_cdf_chord_field2(g, 'hci',
     1   'Ion thermal anomalous diffusivity', 'm^-2 s^-1')
        call create_cdf_chord_field2(g, 'hce',
     1   'Electron thermal anomalous diffusivity', 'm^-2 s^-1')
        call create_cdf_chord_field3(g, 'vlax', 1,
     1   'Poloidal anomalous pinch velocity', 'm s^-1')
        call create_cdf_chord_field3(g, 'vlay', 1,
     1   'Radial anomalous pinch velocity', 'm s^-1')
        call create_cdf_chord_field3(g, 'vsa', 1,
     1   'Anomalous viscosity', 'm kg^-1 s^-1')
      endif
#ifdef B25_EIRENE
      if (switch%write_chord_eirene_state_variables.ge.1) then
        if (nnatmi.gt.0) then
          call create_cdf_chord_field3(g, 'dab2', 2,
     1     'Atom density', 'm^-3')
          call create_cdf_chord_field3(g, 'tab2', 2,
     1     'Atom temperature', 'eV')
        endif
        if (nnmoli.gt.0) then
          call create_cdf_chord_field3(g, 'dmb2', 3,
     1     'Molecular density', 'm^-3')
          call create_cdf_chord_field3(g, 'tmb2', 3,
     1     'Molecule temperature', 'eV')
        endif
      endif
      if (switch%write_chord_eirene_fluxes.ge.1) then
        if (nnatmi.gt.0) then
          call create_cdf_chord_field3(g, 'rfluxa', 2,
     1     'Atomic radial flux density', 'm^-2 s^-1')
          call create_cdf_chord_field3(g, 'pfluxa', 2,
     1     'Atomic poloidal flux density', 'm^-2 s^-1')
          call create_cdf_chord_field3(g, 'refluxa', 2,
     1     'Atomic radial energy flux density', 'W m^-2')
          call create_cdf_chord_field3(g, 'pefluxa', 2,
     1     'Atomic poloidal energy flux density', 'W m^-2')
        endif
        if (nnmoli.gt.0) then
          call create_cdf_chord_field3(g, 'rfluxm', 3,
     1     'Molecular radial flux density', 'm^-2 s^-1')
          call create_cdf_chord_field3(g, 'pfluxm', 3,
     1     'Molecular poloidal flux density', 'm^-2 s^-1')
          call create_cdf_chord_field3(g, 'refluxm', 3,
     1     'Molecular radial energy flux density', 'W m^-2')
          call create_cdf_chord_field3(g, 'pefluxm', 3,
     1     'Molecular poloidal energy flux density', 'W m^-2')
        endif
      endif
      if (switch%write_chord_eirene_sources.ge.1) then
        if (nnatmi.gt.0) then
          call create_cdf_chord_field3(g, 'eirene_papl_sna', 1,
     1     'Particle source from atom-plasma interactions',
     2     'm^-3 s^-1')
          call create_cdf_chord_field3(g, 'eirene_pppl_sna', 1,
     1     'Particle source from recombinations', 'm^-3 s^-1')
          call create_cdf_chord_field2(g, 'eirene_eael_she',
     1     'Electron energy source from atom-plasma interactions',
     2     'W m^-3')
          call create_cdf_chord_field2(g, 'eirene_epel_she',
     1     'Electron energy source from recombinations', 'W m^-3')
          call create_cdf_chord_field2(g, 'eirene_eapl_shi',
     1     'Ion energy source from atom-plasma interactions',
     2     'W m^-3')
          call create_cdf_chord_field2(g, 'eirene_eppl_shi',
     1     'Ion energy source from recombinations', 'W m^-3')
          call create_cdf_chord_field3(g, 'eneutrad', 2,
     1     'Volumetric radiation rate due to atoms', 'W m^-3')
          call create_cdf_chord_field2(g, 'eneutradsum',
     1     'Total volumetric radiation rate due to atoms', 'W m^-3')
        endif
        if (nnmoli.gt.0) then
          call create_cdf_chord_field3(g, 'eirene_pmpl_sna', 1,
     1     'Particle source from molecule-plasma interactions',
     2     'm^-3 s^-1')
          call create_cdf_chord_field2(g, 'eirene_emel_she',
     1     'Electron energy source from molecule-plasma '//
     2     'interactions', 'W m^-3')
          call create_cdf_chord_field2(g, 'eirene_empl_shi',
     1     'Ion energy source from molecule-plasma '//
     2     'interactions', 'W m^-3')
          call create_cdf_chord_field3(g, 'emolrad', 3,
     1     'Volumetric radiation rate due to molecules', 'W m^-3')
          call create_cdf_chord_field2(g, 'emolradsum',
     1     'Total volumetric radiation rate due to molecules',
     2     'W m^-3')
        endif
      endif
#endif
      return
      end subroutine define_b25_fields

      subroutine put_chord_field2(name, cellf)
c     Interpolate a cell-centered field onto the chords and write each
c     chord's segment into its own group.
      implicit none
      character*(*) name
      real (kind=R8) :: cellf(nCv)
      integer mm
      call sample_cv_at_points(geo, mpg, cellf, crx, cry, ichord, npt,
     &  valg)
      do mm = 1, nchord_b
        call put_cdf_chord_field2(grpids(mm), name,
     &    valg(off(mm)+1:off(mm)+cnt(mm)), cnt(mm))
      enddo
      return
      end subroutine put_chord_field2

      subroutine put_chord_field3(name, cellf3)
c     Interpolate a per-species cell-centered field onto the chords and
c     write each chord's segment into its own group. cellf3 columns map
c     to species 1..ns (for a 0:ns-1 array).
      implicit none
      character*(*) name
      real (kind=R8) :: cellf3(nCv, ns)
      integer js, mm
      do js = 1, ns
        call sample_cv_at_points(geo, mpg, cellf3(:,js), crx, cry,
     &    ichord, npt, valg3(:,js))
      enddo
      do mm = 1, nchord_b
        call put_cdf_chord_field3(grpids(mm), name,
     &    valg3(off(mm)+1:off(mm)+cnt(mm),:), cnt(mm), ns)
      enddo
      return
      end subroutine put_chord_field3

#ifdef B25_EIRENE
      subroutine put_chord_fieldn(name, cellf, ncol, valgn)
c     Interpolate an (nCv,ncol) cell-centered field onto the chords and
c     write each chord's segment into its own group. Used for the EIRENE
c     atom/molecule-indexed quantities (ncol = nnatmi or nnmoli).
      implicit none
      character*(*) name
      integer :: ncol
      real (kind=R8) :: cellf(nCv, ncol)
      real (kind=R8) :: valgn(npt, ncol)
      integer :: jc, mm
      do jc = 1, ncol
        call sample_cv_at_points(geo, mpg, cellf(:,jc), crx, cry,
     &    ichord, npt, valgn(:,jc))
      enddo
      do mm = 1, nchord_b
        call put_cdf_chord_field3(grpids(mm), name,
     &    valgn(off(mm)+1:off(mm)+cnt(mm),:), cnt(mm), ncol)
      enddo
      return
      end subroutine put_chord_fieldn

      subroutine put_chord_tri(name, src, nsp, scale)
c     Map a per-EIRENE-cell field src(species,cell) to per-triangle
c     values, scale it, interpolate onto the chords for each of the nsp
c     species, and write each chord's segment into its eirenechords group.
      implicit none
      character*(*) name
      real (kind=R8) :: src(:,:)
      integer :: nsp
      real (kind=R8) :: scale
      integer :: isp, mm
      do isp = 1, nsp
        call build_tri_field(src(isp,1:ncvtot), ncvtot, ftri)
        if (scale.ne.1.0_R8) ftri = ftri*scale
        call sample_tri_at_points(ftri, crx_t, cry_t, ichord_t, npt_t,
     &    valgtri(:,isp))
      enddo
      do mm = 1, nchord_t
        call put_cdf_chord_field3(grpids_t(mm), name,
     &    valgtri(off_t(mm)+1:off_t(mm)+cnt_t(mm),:), cnt_t(mm), nsp)
      enddo
      return
      end subroutine put_chord_tri

      subroutine define_tri_fields(g)
c     Define all eirenechords (triangle-grid) fields in group g.
      implicit none
      integer g
      if (nnatmi.gt.0) then
        call create_cdf_chord_field3(g, 'pdena', 2,
     1   'Atom particle density', 'm^-3')
        call create_cdf_chord_field3(g, 'edena', 2,
     1   'Atom energy density', 'eV m^-3')
        call create_cdf_chord_field3(g, 'vxdena', 2,
     1   'Atom x-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_chord_field3(g, 'vydena', 2,
     1   'Atom y-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_chord_field3(g, 'vzdena', 2,
     1   'Atom z-directed momentum density', 'kg s^-1 m^-2')
      endif
      if (nnmoli.gt.0) then
        call create_cdf_chord_field3(g, 'pdenm', 3,
     1   'Molecular particle density', 'm^-3')
        call create_cdf_chord_field3(g, 'edenm', 3,
     1   'Molecular energy density', 'eV m^-3')
        call create_cdf_chord_field3(g, 'vxdenm', 3,
     1   'Molecular x-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_chord_field3(g, 'vydenm', 3,
     1   'Molecular y-directed momentum density', 'kg s^-1 m^-2')
        call create_cdf_chord_field3(g, 'vzdenm', 3,
     1   'Molecular z-directed momentum density', 'kg s^-1 m^-2')
      endif
      return
      end subroutine define_tri_fields
#endif
#endif

      end subroutine cdfchord
