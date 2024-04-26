      subroutine copy_background(nCv, nFc, ns, switch, mpg, gm,
     &                           pl, dv, st_ext)
      use b2mod_types
      use b2mod_switches
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_diag
     &     , only: fluids_list
      use b2mod_subsys
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_parmmod
     .   , ncv_eir => ncv
      use eirmod_braeir
      use eirmod_extrab25
      use b2mod_neutrals_namelist
     & , only: rcfe, rcfi, mol
#else
      use parmmod   ! IGNORE
     .   , ncv_eir => ncv
      use braeir    ! IGNORE
#endif
#else
      use b2mod_braeir
      use b2mod_neutrals_namelist
     & , only : rcfe, rcfi, mol
#endif
      use b2mod_recycle
     & , only : fna_mol
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer nCv, nFc, ns
      type (switches), intent(in) :: switch
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2Plasma), intent(in) :: pl
      type(B2Derivatives), intent(in) :: dv
      type(B2StateExt), intent (in) :: st_ext

      character chns*3                                                 !srv 16.09.02
      character(2) :: ccall
      integer, save :: icall = 0                                       !srv 21.06.08
      integer is,ncall
      external ipgetr, ipgeti, xerrab
      real (kind=R8) ::
     * te_min_val, te_max_val, ti_min_val, ti_max_val,
     * na_min_val, na_max_val, ua_min_val, ua_max_val,
     * ne_min_val, ne_max_val
      data ncall/0/
      real (kind=R8) :: kintmp(nFc,0:ns-1)
      real (kind=R8) ::
     &  vaecrbc(nCv,0:1,0:ns-1), wadiac(nCv,0:1,0:ns-1),
     &  vadiac(nCv,0:1,0:ns-1)
      real (kind=R8) :: naf(nFc,0:ns-1), fna_tmp(nFc,0:ns-1),
     &  uuc(nCv,0:ns-1),vvc(nCv,0:ns-1),nafsqalf(nFc,0:1,0:ns-1)
      integer :: ic, iFc, ian, ien, ntrgdat, ifl, i, istra, ic2,
     &           iv1, iv2, iCp
      real (kind=R8) :: evi, ebx, eby, ebn, uu, vv, ww, vtot, vmax,
     1 pvxs, pvys, puxs, puys
      real (kind=R8) :: wght(nFc,2)
      real (kind=R8) :: fnn_inc_b(nFc,1:ns), fnn_mol_b(nFc,1:ns)
      real (kind=R8) :: cs(nCv), rz(nCv), pz(nCv)

      call subini ('copy_background')

c Initialization
      fnn_inc_b = 0.0_R8
      fnn_mol_b = 0.0_R8

      evi = 1.0_R8 / ev
      if(switch%eirene_print_minmax.ne.0) then
         te_min_val=minval(pl%te(:))*evi
         te_max_val=maxval(pl%te(:))*evi
         ti_min_val=minval(pl%ti(:))*evi
         ti_max_val=maxval(pl%ti(:))*evi
         na_min_val=minval(pl%na(:,0:ns-1))
         na_max_val=maxval(pl%na(:,0:ns-1))
         ne_min_val=minval(dv%ne(:))
         ne_max_val=maxval(dv%ne(:))
         ua_min_val=minval(pl%ua(:,0:ns-1))
         ua_max_val=maxval(pl%ua(:,0:ns-1))
         write(*,'(a,1p,10(1x,g12.3))') 'minmax: te, ti, na, ua, ne',
     &     te_min_val, te_max_val, ti_min_val, ti_max_val,
     &     na_min_val, na_max_val, ua_min_val, ua_max_val,
     &     ne_min_val, ne_max_val
      endif

* get drift velocities in cell centers
      do is=0,ns-1
        call intcell(nFc,nCv,mpg,mpg%intcellP,
     &               dv%vaecrb(:,0,is),vaecrbc(:,0,is))
        call intcell(nFc,nCv,mpg,mpg%intcellP,
     &               dv%wadia(:,0,is),wadiac(:,0,is))
        call intcell(nFc,nCv,mpg,mpg%intcellP,
     &               dv%vadia(:,0,is),vadiac(:,0,is))
        call intcell(nFc,nCv,mpg,mpg%intcellR,
     &               dv%vaecrb(:,1,is),vaecrbc(:,1,is))
        call intcell(nFc,nCv,mpg,mpg%intcellR,
     &               dv%wadia(:,1,is),wadiac(:,1,is))
        call intcell(nFc,nCv,mpg,mpg%intcellR,
     &               dv%vadia(:,1,is),vadiac(:,1,is))
      enddo

      if(switch%copy_background_iout.ne.0) then                                        !srv 18.06.08 {
        do is = 0,ns-1
         write (chns,'(i3.3)') is
         write (ccall,'(a1,i1)') '_', icall
         call my_out_us(70,nFc,1,dv%fna(1,0,is),'copy_fna_th'//chns)
         call my_out_us(70,nFc,1,dv%fna(1,1,is),'copy_fna_r'//chns)
         call my_out_us(70,nFc,1,dv%fna_eir(1,0,is),
     &     'copy_fna_eir_th'//chns)
         call my_out_us(70,nFc,1,dv%fna_eir(1,1,is),
     &     'copy_fna_eir_r'//chns)
         call my_out_us(70,nFc,1,dv%vaecrb(1,0,is),
     &     'copy_vecrb_th'//chns)
         call my_out_us(70,nFc,1,dv%vaecrb(1,1,is),'copy_vecrb_r'//chns)
         call my_out_us(70,nFc,1,dv%wadia(1,0,is),'copy_wadia_th'//chns)
         call my_out_us(70,nFc,1,dv%wadia(1,1,is),'copy_wadia_r'//chns)
         call my_out_us(70,nFc,1,dv%vadia(1,0,is),'copy_vadia_th'//chns)
         call my_out_us(70,nFc,1,dv%vadia(1,1,is),'copy_vadia_r'//chns)
         call my_out_us(70,nCv,0,vaecrbc(1,0,is),'copy_vecrbc_th'//chns)
         call my_out_us(70,nCv,0,vaecrbc(1,1,is),'copy_vecrbc_r'//chns)
         call my_out_us(70,nCv,0,wadiac(1,0,is),'copy_wadiac_th'//chns)
         call my_out_us(70,nCv,0,wadiac(1,1,is),'copy_wadiac_r'//chns)
         call my_out_us(70,nCv,0,vadiac(1,0,is),'copy_vadiac_th'//chns)
         call my_out_us(70,nCv,0,vadiac(1,1,is),'copy_vadiac_r'//chns)
         call my_out_us(70,nCv,0,pl%ua(1,is),'copy_ua'//chns//ccall)
         call my_out_us(70,nCv,0,pl%na(1,is),'copy_na'//chns)
        enddo
        call my_out_us(70,nFc,1,dv%fne_eir(1,0),'copy_fne_eir_th')
        call my_out_us(70,nFc,1,dv%fne_eir(1,1),'copy_fne_eir_r')
        call my_out_us(70,nFc,1,dv%fhe_eir(1,0),'copy_fhe_eir_th')
        call my_out_us(70,nFc,1,dv%fhe_eir(1,1),'copy_fhe_eir_r')
        call my_out_us(70,nFc,1,dv%fhi_eir(1,0),'copy_fhi_eir_th')
        call my_out_us(70,nFc,1,dv%fhi_eir(1,1),'copy_fhi_eir_r')
        call my_out_us(70,nCv,0,dv%ne,'copy_ne')
        call my_out_us(70,nCv,0,pl%te,'copy_te')
        call my_out_us(70,nCv,0,pl%ti,'copy_ti')
        call my_out_us(70,nCv,0,pl%po,'copy_po')
        call my_out_us(70,nCv,0,gm%cvBb(1,0),'copy_bb_0')
        call my_out_us(70,nCv,0,gm%cvBb(1,1),'copy_bb_1')
        call my_out_us(70,nCv,0,gm%cvBb(1,2),'copy_bb_2')
        call my_out_us(70,nCv,0,gm%cvBb(1,3),'copy_bb_3')
      endif

      ! precompute sound speed for velocity limiting
      call b2xprz (nCv, ns, mp, am, pl%na, rz, st_ext)
      call b2xppz (nCv, ns, dv%ne, pl%na, pl%te, pl%ti, pz, st_ext)
      cs = sqrt(pz/rz)

#ifndef SOLPS4_3
* cell centred
      do ic=1,nCv
        neb(ic)=max(min(dv%ne(ic),switch%eir_ne_max),switch%eir_ne_min)
        teb(ic)=max(min(pl%te(ic)*evi,switch%eir_te_max),
     &                                 switch%eir_te_min)
        tib(ic)=max(min(pl%ti(ic)*evi,switch%eir_ti_max),
     &                                 switch%eir_ti_min)
        tnb(ic)=max(min(pl%tn(ic)*evi,switch%eir_ti_max),
     &                                 switch%eir_ti_min)
        pob(ic)=pl%po(ic)
        prb(ic)=teb(ic)*max(min(dv%ne(ic),switch%eir_na_max),
     &                                 switch%eir_na_min)
        rrb(ic)=gm%cvBb(ic,0)/gm%cvBb(ic,3)
        volb(ic)=gm%cvVol(ic)
        bbxb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,0)
        bbyb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,1)
        bbzb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,2)
        bfeldb(ic)=gm%cvBb(ic,3)
      enddo

* faces
!      call intface(nCv,nFc,mpg%fcCv,gm%fcVol,dv%kinrgy,kintmp)
      call intface1(nCv,nFc,ns-1,mpg%fcCv,gm%fcHc,dv%kinrgy,kintmp)
      do ifc = 1, nFc
        feib(ifc) = dv%fhi(ifc,0) + dv%fhi(ifc,1)
        feeb(ifc) = dv%fhe(ifc,0) + dv%fhe(ifc,1)
        do is = 0, ns-1
          feib(ifc) = feib(ifc) + kintmp(ifc,is)*(1.0_R8-switch%BoRiS)
     .                *(dv%fna(ifc,0,is) + dv%fna(ifc,1,is))
        end do
      end do
#endif
* by ns

*   ..interpolate density to cell faces
!pb      call intface1(nCv,nFc,ns-1,mpg%fcCv,gm%fcVol,pl%na,naf)
      wght = 1.0_R8
      call intface1(nCv,nFc,ns-1,mpg%fcCv,wght,pl%na,naf)

      if (switch%eirene_uub_style.eq.1) then
        do is = 0,ns-1
          nafsqalf(1:nFc,0,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,0)*
     .                           naf(1:nFc,is)+1.e-30_R8
        end do
        if (switch%eirene_background.eq.0) then
          fna_tmp(1:nFc,0:ns-1) = dv%fna(1:nFc,0,0:ns-1)/
     .                            nafsqalf(1:nFc,0,0:ns-1)
        else
          fna_tmp(1:nFc,0:ns-1) = dv%fna_eir(1:nFc,0,0:ns-1)/
     .                            nafsqalf(1:nFc,0,0:ns-1)
        endif
        call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellP,fna_tmp,uuc)
      end if

      do is = 0,ns-1
        nafsqalf(1:nFc,1,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,1)*
     .                         naf(1:nFc,is)+1.e-30_R8
        if (switch%eirene_background.eq.0) then
          fna_tmp(1:nFc,is) = dv%fna(1:nFc,1,is)/
     .                        nafsqalf(1:nFc,1,is)
          do ic = mpg%nCi+1, nCv   ! loop over all guard cells
            iFc = mpg%cvFc(mpg%cvFcP(ic,1))
            iCp = mpg%fcCv(iFc,1) + mpg%fcCv(iFc,2) - ic  ! index of plasma cell
            fna_tmp(iFc,is) = dv%fna(iFc,1,is)/
     .                       (gm%fcS(iFc)*gm%fcQalf(iFc,1)*
     .                        pl%na(ic,is)+1.e-30_R8)
          end do
        else
          fna_tmp(1:nFc,is) = dv%fna_eir(1:nFc,1,is)/
     .                        nafsqalf(1:nFc,1,is)
        endif
      end do
      call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellR,fna_tmp,vvc)

c for neutrals (needed on Eirene side in case of spatially hybrid neutrals),
c always use fna instead of fna_eir, regardless of eirene_background switch
      if (switch%spatial_hybrid.ne.0) then
        do is=0,ns-1
          if (is_neutral(is)) then
            fna_tmp(1:nFc,is) = dv%fna(1:nFc,1,is)/
     .                          nafsqalf(1:nFc,1,is)
            call intcell(nFc,nCv,mpg,mpg%intcellR,
     .              fna_tmp(1:nFc,is),vvc(1:nFc,is))
          endif
        enddo
      endif

#ifndef SOLPS4_3
      nflai=0
      do is=0,ns-1
!!        if(.not.is_neutral(is)) then
!!          nflai=nflai+1
          if (fluids_list(is) == 0) cycle
          nflai = fluids_list(is)
          call xertst (nflai.le.NFL, 'eirene_mc: increase size of nfl')
          amb(nflai)=am(is)
          znb(nflai)=zn(is)
          if (zamin(is) == zamax(is)) then
            zab(nflai)=zamin(is)
          else
            call xerrab ('eirene_mc: do not use bundled ' //
     .                   'species for Eirene')
          end if
* cell centred
          do ic = 1, nCv
            dnib(ic,nflai)=max(min(pl%na(ic,is),switch%eir_na_max),
     1                                           switch%eir_na_min)
            upb(ic,nflai)=pl%ua(ic,is)
            uudiab(ic,nflai)=wadiac(ic,0,is)+vaecrbc(ic,0,is)
            vvdiab(ic,nflai)=wadiac(ic,1,is)+vaecrbc(ic,1,is)
            if(.not.is_neutral(is)) then
              prb(ic)=prb(ic)+tib(ic)*dnib(ic,nflai)
            else
              prb(ic)=prb(ic)+tnb(ic)*dnib(ic,nflai)
            endif
c     Velocities uub, vvb, wwb
c     * cell-centered
c     * internal cell values used for plasma background properties in infcop
c     * guard cell values used for source sampling in infcop
            if(switch%eirene_uub_style.eq.0) then
              if (switch%eirene_background.eq.0.or.ic.gt.mpg%nCi) then
c     internal cell, old treatment -> grad-B drift only
c     for boundary cells, -> always grad-B drift only, regardless of eirene_background switch
                uub(ic,nflai)=
     &               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     &               switch%xwdia*vadiac(ic,0,is)+
     &               switch%xvecrb*vaecrbc(ic,0,is)
              else
c     internal cell -> full diamagnetic drift
                uub(ic,nflai)=
     &               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     &               switch%xwdia*wadiac(ic,0,is)+
     &               switch%xvecrb*vaecrbc(ic,0,is)
                 endif
            elseif(switch%eirene_uub_style.eq.1) then
cwdk  uub computed from face fluxes, and interpolated to cell centers
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at N/S boundaries questionable at present, because (certain/all) terms
c     in fna may have been set to 0. Should consider using the eirene_uub_style.eq.0 treatment
c     in these cells!
              uub(ic,nflai)=uuc(ic,is)
            end if

cwdk  vvb computed from face fluxes, and interpolated to cell centers
c     Compatible with SOLPS4.3.
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at E/W boundaries questionable, because (certain/all) terms
c     in fna may have been set to 0.
            vvb(ic,nflai)=vvc(ic,is)
cwdk  wwb directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but may not be fully consistent with respect to drifts.
c     In particular: convective anomalous contribution (if present) is missing
c                    (diffusive anomalous contribution zero due to toroidal symm.)
            if (switch%eirene_background.eq.0.or.ic.gt.mpg%nCi) then
c     internal cell, old treatment -> grad-B drift only (fna)
c     for boundary cells, -> always grad-B drift only, regardless of eirene_background switch
                wwb(ic,nflai)=
     &            upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     &            (switch%xwdia*vadiac(ic,0,is)+
     &             switch%xvecrb*vaecrbc(ic,0,is))*
     &            gm%cvBb(ic,0)/gm%cvBb(ic,2)
            else
c     internal cell -> full diamagnetic drift
                 wwb(ic,nflai)=
     &             upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     &             (switch%xwdia*wadiac(ic,0,is)+
     &              switch%xvecrb*vaecrbc(ic,0,is))*
     &             gm%cvBb(ic,0)/gm%cvBb(ic,2)
            end if

            ebx = gm%cvEb(ic,0)
            eby = gm%cvEb(ic,1)
            ebn = sqrt(ebx*ebx + eby*eby)
            if (ebn > 0.0_R8) then
              ebx = ebx / ebn
              eby = eby / ebn
            endif
            pux(ic) = ebx
            puy(ic) = eby
            if (gm%cvEb(ic,2).ge.0.0_R8) then
              pvx(ic) = -eby
              pvy(ic) = ebx
            else
              pvx(ic) = eby
              pvy(ic) = -ebx
            endif

            ! limit the total velocity, and then scale
            ! all individual components
            vtot = sqrt(uub(ic,nflai)**2 + vvb(ic,nflai)**2 +
     &        wwb(ic,nflai)**2)
            vmax = min(switch%eir_ua_max,abs(switch%eir_ua_min),
     &        switch%eir_M_max*cs(ic))
            if (vtot.gt.vmax) then
              uub(ic,nflai) = uub(ic,nflai)*vmax/vtot
              vvb(ic,nflai) = vvb(ic,nflai)*vmax/vtot
              wwb(ic,nflai) = wwb(ic,nflai)*vmax/vtot
              upb(ic,nflai) = upb(ic,nflai)*vmax/vtot
              uudiab(ic,nflai) = uudiab(ic,nflai)*vmax/vtot
              vvdiab(ic,nflai) = vvdiab(ic,nflai)*vmax/vtot
            endif

            vxb(ic,nflai) = uub(ic,nflai)*pux(ic) +
     1                      vvb(ic,nflai)*pvx(ic)
            vyb(ic,nflai) = uub(ic,nflai)*puy(ic) +
     1                      vvb(ic,nflai)*pvy(ic)
            vzb(ic,nflai) = wwb(ic,nflai)*sign(1.0_R8,gm%cvEb(ic,2))
          enddo
* cell faces
          do iFc = 1, nFc
            fnib(iFc,nflai) = dv%fna(iFc,0,is) + dv%fna(iFc,1,is)
            if (is_neutral(is)) then
              fnn_inc_b(iFc,nflai) = dv%fnn_inc(iFc,is)
              fnn_mol_b(iFc,nflai) = fna_mol(iFc,is)
            endif
          end do

!!        endif
      enddo

!pb  05.03.18
      nflai = count(fluids_list > 0,1)

      if(switch%eirene_extrap.eq.0) then
        ! do nothing
      elseif(switch%eirene_extrap.eq.1) then
c we derive the guard cell value from its neighbouring interior cell
        do ic = mpg%nCi+1, mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(ic,1))
          ic2 = min(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))   ! index of plasma cell
          teb(ic)=teb(ic2)
          tib(ic)=tib(ic2)
          tnb(ic)=tnb(ic2)
          pob(ic)=pob(ic2)
          dnib(ic,:)=dnib(ic2,:)
        enddo
      else
        write(*,*) switch%eirene_extrap,
     1   ' not available as an eirene_extrap option'
      endif

      if(switch%eirene_print_minmax.ne.0) then
         te_min_val=minval(teb(:))
         te_max_val=maxval(teb(:))
         ti_min_val=minval(tib(:))
         ti_max_val=maxval(tib(:))
         na_min_val=minval(dnib(:,1:nflai))
         na_max_val=maxval(dnib(:,1:nflai))
         ua_min_val=minval(uub(:,1:nflai))
         ua_max_val=maxval(uub(:,1:nflai))
         write(*,'(a,1p,8(1x,g12.3))') 'minmax:teb,tib,nab,uub ',
     &     te_min_val, te_max_val, ti_min_val, ti_max_val,
     &     na_min_val, na_max_val, ua_min_val, ua_max_val
      endif

      if(switch%copy_background_iout.ne.0) then
       nflai=0
       do is=0,ns-1
        if(.not.is_neutral(is)) then
         nflai=nflai+1
         write (chns,'(i3.3)') is
         write (ccall,'(a1,i1)') '_', icall
         call my_out_us(70,nCv,0,uub(1,nflai),'copy_uub'//chns//ccall)
         call my_out_us(70,nCv,0,vvb(1,nflai),'copy_vvb'//chns//ccall)
         call my_out_us(70,nCv,0,wwb(1,nflai),'copy_wwb'//chns//ccall)
         call my_out_us(70,nCv,0,upb(1,nflai),'copy_upb'//chns//ccall)
        endif
       enddo
       call my_out_us(70,nCv,0,teb,'copy_teb')
       call my_out_us(70,nCv,0,tib,'copy_tib')
       call my_out_us(70,nCv,0,tnb,'copy_tnb')
       call my_out_us(70,nCv,0,pob,'copy_pob')
       icall = icall + 1
      endif

!  .. define target data

      if (.not.allocated(trgt)) allocate(trgt(mpg%nRc))
      do istra = 1, mpg%nRc
        ian = mpg%rcFcP(istra,1)
        ien = ian + mpg%rcFcP(istra,2) - 1
        ntrgdat = ien-ian+1
        if (ien-ian+1 <= 0) then
            trgt(istra)%ntrgdat = 0
            cycle
        endif

#ifdef B25_EIRENE
        call eirene_alloc_target_data(istra,ntrgdat,nflai)
#else
        call alloc_target_data(istra,ntrgdat,nflai)
#endif

        do i = 1, ntrgdat
          iFc = mpg%rcFc(ian+i-1)
          trgt(istra)%faces(i) = iFc
          trgt(istra)%fcori(i) = mpg%rcFcOr(ian+i-1)
          iv1 = mpg%fcVx(ifc,1)
          iv2 = mpg%fcVx(ifc,2)
          ic = mpg%rcCv(ian+i-1,1)   ! index of guard cell
          trgt(istra)%flength(i) = sqrt((gm%vxX(iv1)-gm%vxX(iv2))**2 +
     >                            (gm%vxY(iv1)-gm%vxY(iv2))**2)*100._R8
          trgt(istra)%tet(i) = teb(ic)
          trgt(istra)%tit(i) = tib(ic)
          trgt(istra)%tnt(i) = tnb(ic)
          trgt(istra)%fe(i) = rcfe(istra)
          trgt(istra)%fsh(i) = delta_sheathb(iFc)
          do ifl = 1, nflai
            trgt(istra)%dnit(i,ifl) = dnib(ic,ifl)
            if (zab(ifl).ne.0.and.mol(istra).eq.0) then !! ion species
              trgt(istra)%flux(i,ifl) = fnib(iFc,ifl)*
     &                mpg%rcFcOr(ian+i-1)*dv%kin_frac_hyb(iFc)
            elseif (zab(ifl).eq.0.and.mol(istra).eq.0) then !! neutral atom species
              trgt(istra)%flux(i,ifl) = fnn_inc_b(iFc,ifl)*
     &                mpg%rcFcOr(ian+i-1)*dv%kin_frac_hyb(iFc)
              fnib(iFc,ifl) = fnn_inc_b(iFc,ifl)
            elseif (zab(ifl).eq.0.and.mol(istra).eq.1) then     !! neutral molecule species
              trgt(istra)%flux(i,ifl) = 2.0_R8*fnn_mol_b(iFc,ifl)*
     &                                  dv%fluid_frac_hyb(iFc)
              fnib(iFc,ifl) = 0.0_R8
            else
              if (fnn_mol_b(iFc,ifl).gt.0) then
                call xerrab ('no molecular contribution expected here')
              endif
            endif
            if (trgt(istra)%flux(i,ifl).lt.0.0_R8) then
              write(*,*) 'negative flux !! istra, i, ifl = ',istra,i,ifl
              trgt(istra)%flux(i,ifl) = 0.0_R8
              fnib(iFc,ifl) = 0.0_R8
            endif
            trgt(istra)%fi(i,ifl) = rcfi(istra)

            ebx = gm%fcEb(iFc,0)
            eby = gm%fcEb(iFc,1)
            ebn = sqrt(ebx*ebx + eby*eby)
            if (ebn > 0) then
              ebx = ebx / ebn
              eby = eby / ebn
            endif
            puxs = ebx
            puys = eby
            if (gm%fcEb(iFc,2).ge.0.0_R8) then
              pvxs = -eby
              pvys = ebx
            else
              pvxs = eby
              pvys = -ebx
            endif

            uu = uub(ic,ifl)
            vv = vvb(ic,ifl)
            ww = wwb(ic,ifl)
            trgt(istra)%vxt(i,ifl) = puxs*uu + pvxs*vv
            trgt(istra)%vyt(i,ifl) = puys*uu + pvys*vv
            trgt(istra)%vzt(i,ifl) = ww*sign(1.0_R8,gm%fcEb(iFc,2))
            trgt(istra)%fel(i,ifl) = 0.0_R8
            trgt(istra)%vpart(i,ifl) = 0.0_R8
            trgt(istra)%mach(i,ifl) = 1.0_R8
            trgt(istra)%usr(i,ifl) = 0.0_R8
            trgt(istra)%zit(i,ifl) = zab(ifl)
          end do
        end do
      end do
#endif

      ncall=ncall+1
      call subend ()
      return
      end subroutine copy_background

      subroutine write_f30(geo,mpg)
      use b2us_map
      use b2us_geo
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: geo
      integer :: i, nB, iLbl, minLbl, maxLbl, nBFaces, nBVerts
      integer, allocatable :: bVerts(:)

      write(*,*) 'Writing fort.30 file'

      ! open the file
      open(unit=30,file='fort.30')

      ! identify range of face labels
      minLbl = minval(mpg%fcLbl)
      maxLbl = maxval(mpg%fcLbl)

      ! first pass: identify number of boundary segments to be written
      nB = 0
      do iLbl = minLbl,maxLbl

        ! only treat faces with non-zero fclbl
        if (iLbl.eq.0) cycle

        ! determine number of faces with this label
        nBFaces = count(mpg%fcLbl==iLbl)

        if (nBFaces.gt.0) nB = nB + 1

      enddo

      ! write file
      write(30,'(i3)') nB

      do iLbl = minLbl,maxLbl

        ! only treat faces with non-zero fclbl
        if (iLbl.eq.0) cycle

        ! determine number of faces with this label
        nBFaces = count(mpg%fcLbl==iLbl)

        ! corresponding number of vertices
        nBVerts = nBFaces + 1

        if (nBFaces.gt.0) then

          allocate (bVerts(nBVerts))

          ! sort vertices of the boundary faces into a polygon
          call sort_vertices(nBVerts,bVerts,mpg%nFc,mpg%fcLbl,iLbl,mpg)

          ! write into fort.30
          write(30,'(i3)') iLbl
          write(30,'(i12)') nBVerts
          do i = 1,nBVerts
            write(30,'(2(2X,E21.14))') 
     .        geo%vxX(bVerts(i)), geo%vxY(bVerts(i))
          end do

          deallocate (bVerts)

        end if
      end do

      ! close file
      close(30)

      return
      end subroutine write_f30
c
      subroutine write_f31(mpg)
      use b2mod_types
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_parmmod
      use eirmod_braeir
      use eirmod_extrab25
      use eirmod_wneutrals
#else
      use parmmod   ! IGNORE
      use braeir    ! IGNORE
#endif
#else
      use b2mod_braeir
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      type(mapping), intent(in) :: mpg
      integer :: i, ic, ifl, ntrgdat, istra

!  .. use approach to generalized user interface

#ifndef SOLPS4_3
!  .. plasma data
      rewind(31)
      write (31,'(a)') '* EIRENE PLASMA FILE '
      write (31,'(a)') '*'
      write (31,'(a)') '*** MISC PLASMA DATA'
      write (31,'(a)') '*'
      write (31,'(a1,a6,15a15)') '*','Index','Te','ne',
     .      'Bx','By','Bz','Ex','Ey','Ez','pot','psi','f',
     .      'Pux','Puy','Pvx','Pvy'
      write (31,'(a1,a6,15a15)') '*',' ','(eV)','(1/m^3)',
     .      'Tesla','Tesla','Tesla','(V/m)','(V/m)','(V/m)',
     .      'V','??','??',' ',' ',' ',' '

      write (31,'(3i6)') mpg%nCv, mpg%nCi, mpg%nFc

      do ic = 1, mpg%nCi
        write (31,'(i7,15es15.7)') ic, teb(ic), neb(ic),
     .      bbxb(ic), bbyb(ic), bbzb(ic), 0.0_R8, 0.0_R8, 0.0_R8,
     .      pob(ic), 0.0_R8, 0.0_R8,
     .      pux(ic), puy(ic), pvx(ic), pvy(ic)
      end do

      do ifl = 1, nflai

        write (31,'(a)') '*'
        write (31,'(a,i0,a)') '*** ION #',ifl,' PLASMA DATA'
        write (31,'(a)') '*'
        write (31,'(a1,a6,5a15)') '*','Index','Ti','ni',
     .      'Vx','Vy','Vz'
        write (31,'(a1,a6,5a15)') '*',' ','(eV)','(1/m^3)',
     .      'm/s','m/s','m/s'

        write (31,'(3es15.7)') amb(ifl), znb(ifl), zab(ifl)

        write (31,'(i6)') mpg%nCi

        do ic = 1, mpg%nCi
          write (31,'(i7,5es15.7)') ic, tib(ic), dnib(ic,ifl),
     .      vxb(ic,ifl), vyb(ic,ifl), vzb(ic,ifl)
        end do

      end do

!  .. target data
      write (31,'(a)') '*'
      write (31,'(a)') '*** TARGET DATA'
      write (31,'(a)') '*'
!      write (31,'(i6)') mpg%nRc
      write (31,'(i6)') count(trgt(1:mpg%nRc)%ntrgdat > 0)

      do istra = 1, mpg%nRc
        ntrgdat = trgt(istra)%ntrgdat
        if (ntrgdat == 0) cycle
        write (31,'(a)') '*'
        write (31,'(a,i0)') '*** MISC TARGET DATA #',istra
        write (31,'(a)') '*'
        write (31,'(a1,2a6,6a15)') '*','Index','face','ori','Te',
     .      'fe','fsh','length'
        write (31,'(a1,2a6,6a15)') '*',' ',' ',' ','(eV)',' ',' ',' '
        write (31,'(i6)') ntrgdat

        do i = 1, ntrgdat
          write (31,'(i7,i6,6es15.7)')
     .      i,trgt(istra)%faces(i),trgt(istra)%fcori(i),
     .        trgt(istra)%tet(i),
     .        trgt(istra)%fe(i),trgt(istra)%fsh(i),
     .        trgt(istra)%flength(i)
        end do

        do ifl = 1, nflai

          write (31,'(a)') '*'
          write (31,'(a,i0,a,i0)') '*** ION #',ifl,
     .                             ' TARGET DATA #',istra
          write (31,'(a)') '*'

          write (31,'(a1,a6,12a15)') '*','Index','Flux','Ti','ni',
     .      'Vx','Vy','Vz','fi','fel','vpart','mach','usr','Zi'
          write (31,'(a1,a6,12a15)') '*',' ',' ','(eV)','(1/m^3)',
     .      'm/s','m/s','m/s',' ',' ',' ',' ',' '

          write (31,'(3es15.7)') amb(ifl), znb(ifl), zab(ifl)
          write (31,'(i6)') ntrgdat

          do i = 1, ntrgdat
            write (31,'(i7,12es15.7)')
     .        i,trgt(istra)%flux(i,ifl),trgt(istra)%tit(i),
     .          trgt(istra)%dnit(i,ifl),trgt(istra)%vxt(i,ifl),
     .          trgt(istra)%vyt(i,ifl),trgt(istra)%vzt(i,ifl),
     .          trgt(istra)%fi(i,ifl),trgt(istra)%fel(i,ifl),
     .          trgt(istra)%vpart(i,ifl),trgt(istra)%mach(i,ifl),
     .          trgt(istra)%usr(i,ifl),trgt(istra)%zit(i,ifl)
          end do
        end do ! ifl
      end do   ! istra
      flush(31)
#endif
      return
      end subroutine write_f31
      
