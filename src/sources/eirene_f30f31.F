      subroutine copy_background(nCv, nFc, ns, switch, mpg, gm, 
     &                           pl, dv, BoRiS, kin_frac_hyb,
     &                           fluid_frac_hyb,fnn_inc)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_switches
      use b2mod_geo
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_diag
     &     , only: fluids_list
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_parmmod
     .   , ncv_eir => ncv
      use eirmod_braeir
      use eirmod_extrab25
      use b2mod_neutrals_namelist
     & , only: crcstra, rcfe, rcfi, mol
#else
      use parmmod   ! IGNORE
      use braeir    ! IGNORE
#endif
#else
      use b2mod_braeir
      use b2mod_neutrals_namelist
     & , only : crcstra, rcfe, rcfi, mol
#endif
      use b2mod_recycle
     & , only : fna_mol
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer nCv, nFc, ns
      real (kind=R8) :: BoRiS
      type (switches), intent (in) :: switch
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2Plasma), intent(in) :: pl
      type(B2Derivatives), intent(in) :: dv

      character chns*3                                                  !srv 16.09.02
      character(2) :: ccall
      integer, save :: icall = 0             !srv 21.06.08
      integer ix,iy,is,ncall,in
      external ipgetr, ipgeti, subini, subend, xerrab
      real (kind=R8) ::
     * te_min, te_max, ti_min, ti_max, na_min, na_max, ua_min, ua_max,
     * te_min_val, te_max_val, ti_min_val, ti_max_val,
     * na_min_val, na_max_val, ua_min_val, ua_max_val,
     * ne_min_val, ne_max_val, ne_min, ne_max
      data ncall/0/
      real (kind=R8) :: kintmp(nFc,0:ns-1)
      real (kind=R8) :: 
     &  vaecrbc(nCv,0:1,0:ns-1), wadiac(nCv,0:1,0:ns-1),
     &  vadiac(nCv,0:1,0:ns-1), wrk(nCv)
      real (kind=R8) :: naf(nFc,0:ns-1), fna_tmp(nFc,0:1,0:ns-1),
     &  uuc(nCv,0:ns-1),vvc(nCv,0:ns-1),nafsqalf(nFc,0:1,0:ns-1)
      integer :: ic, ifc, ian, ien, ntrgdat, ifl, i, istra, ic2,
     &           iv1, iv2, icp
      integer, save :: check_target=1
      real (kind=R8) :: intface_s
      real (kind=R8) :: evi, ebx, eby, ebn, uu, vv, ww, pnorm,
     1 pvxs, pvys, pupv, puxs, puys, puxic, puyic, pvxic, pvyic   
      real (kind=R8) :: wght(nFc,2)
      real (kind=R8), intent(in) ::kin_frac_hyb(nFc),
     &                   fluid_frac_hyb(nFc),fnn_inc(nFc,0:ns-1)
      real (kind=R8) :: fnn_inc_b(nFc,1:ns), fnn_mol_b(nFc,1:ns)

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
         call my_out_us(70,nFc,1,dv%fna(1,0,is),'copy_fnax'//chns)
         call my_out_us(70,nFc,1,dv%fna(1,1,is),'copy_fnay'//chns)
         call my_out_us(70,nFc,1,dv%fna_eir(1,0,is),
     &     'copy_fna_eirx'//chns)
         call my_out_us(70,nFc,1,dv%fna_eir(1,1,is),
     &     'copy_fna_eiry'//chns)
         call my_out_us(70,nFc,1,dv%vaecrb(1,0,is),'copy_vecrbx'//chns)
         call my_out_us(70,nFc,1,dv%vaecrb(1,1,is),'copy_vecrby'//chns)
         call my_out_us(70,nFc,1,dv%wadia(1,0,is),'copy_wadiax'//chns)
         call my_out_us(70,nFc,1,dv%wadia(1,1,is),'copy_wadiay'//chns)
         call my_out_us(70,nFc,1,dv%vadia(1,0,is),'copy_vadiax'//chns)
         call my_out_us(70,nFc,1,dv%vadia(1,1,is),'copy_vadiay'//chns)
         call my_out_us(70,nCv,0,vaecrbc(1,0,is),'copy_vecrbcx'//chns)
         call my_out_us(70,nCv,0,vaecrbc(1,1,is),'copy_vecrbcy'//chns)
         call my_out_us(70,nCv,0,wadiac(1,0,is),'copy_wadiacx'//chns)
         call my_out_us(70,nCv,0,wadiac(1,1,is),'copy_wadiacy'//chns)
         call my_out_us(70,nCv,0,vadiac(1,0,is),'copy_vadiacx'//chns)
         call my_out_us(70,nCv,0,vadiac(1,1,is),'copy_vadiacy'//chns)
         call my_out_us(70,nCv,0,pl%ua(1,is),'copy_ua'//chns//ccall)
         call my_out_us(70,nCv,0,pl%na(1,is),'copy_na'//chns)
        enddo
        call my_out_us(70,nFc,1,dv%fne_eir(1,0),'copy_fne_eirx')
        call my_out_us(70,nFc,1,dv%fne_eir(1,1),'copy_fne_eiry')
        call my_out_us(70,nFc,1,dv%fhe_eir(1,0),'copy_fhe_eirx')
        call my_out_us(70,nFc,1,dv%fhe_eir(1,1),'copy_fhe_eiry')
        call my_out_us(70,nFc,1,dv%fhi_eir(1,0),'copy_fhi_eirx')
        call my_out_us(70,nFc,1,dv%fhi_eir(1,1),'copy_fhi_eiry')
        call my_out_us(70,nCv,0,dv%ne,'copy_ne')
        call my_out_us(70,nCv,0,pl%te,'copy_te')
        call my_out_us(70,nCv,0,pl%ti,'copy_ti')
        call my_out_us(70,nCv,0,pl%po,'copy_po')
        call my_out_us(70,nCv,0,gm%cvBb(1,0),'copy_bb_0')
        call my_out_us(70,nCv,0,gm%cvBb(1,1),'copy_bb_1')
        call my_out_us(70,nCv,0,gm%cvBb(1,2),'copy_bb_2')
        call my_out_us(70,nCv,0,gm%cvBb(1,3),'copy_bb_3')
      endif 
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
          feib(ifc) = feib(ifc) + kintmp(ifc,is)*(1.0_R8-BoRiS)
     .                *(dv%fna(ifc,0,is) + dv%fna(ifc,1,is))
        end do
      end do

* by ns

*   ..interpolate density to cell faces
!pb      call intface1(nCv,nFc,ns-1,mpg%fcCv,gm%fcVol,pl%na,naf)
      wght = 1.0_R8
      call intface1(nCv,nFc,ns-1,mpg%fcCv,wght,pl%na,naf)

      if (switch%eirene_uub_style.eq.1) then
        do is = 0,ns-1
          nafsqalf(1:nFc,0,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,0)*
     .                           naf(1:nFc,is)+1.E-30_R8
        end do
        if (switch%eirene_background.eq.0) then
          fna_tmp(1:nFc,0,0:ns-1) = dv%fna(1:nFc,0,0:ns-1)/
     .                              nafsqalf(1:nFc,0,0:ns-1)
        else
          fna_tmp(1:nFc,0,0:ns-1) = dv%fna_eir(1:nFc,0,0:ns-1)/
     .                              nafsqalf(1:nFc,0,0:ns-1)
        endif
        call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellP,
     .                fna_tmp(:,0,:),uuc)
      end if

      do is = 0,ns-1
        nafsqalf(1:nFc,1,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,1)*
     .                         naf(1:nFc,is)+1.E-30_R8
        if (switch%eirene_background.eq.0) then
          fna_tmp(1:nFc,1,is) = dv%fna(1:nFc,1,is)/
     .                          nafsqalf(1:nFc,1,is)
          do ic = mpg%nCi+1, nCv   ! loop over all guard cells
            iFc = mpg%cvFc(mpg%cvFcP(ic,1))
            icp = mpg%fcCv(ifc,1) + mpg%fcCv(iFc,2) - ic  ! index of plasma cell
            fna_tmp(iFc,1,is) = dv%fna(iFc,1,is)/
     .                         (gm%fcS(iFc)*gm%fcQalf(iFc,1)*
     .                          pl%na(ic,is)+1.E-30_R8)
          end do
        else
          fna_tmp(1:nFc,1,is) = dv%fna_eir(1:nFc,1,is)/
     .                            nafsqalf(1:nFc,1,is)
        endif
      end do
      call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellR,
     .              fna_tmp(:,1,:),vvc)

c for neutrals (needed on Eirene side in case of spatially hybrid neutrals),
c always use fna instead of fna_eir, regardless of eirene_background switch
      if (switch%spatial_hybrid.ne.0) then 
        do is=0,ns-1
          if (is_neutral(is)) then
            fna_tmp(1:nFc,1,is) = dv%fna(1:nFc,1,is)/
     .                            nafsqalf(1:nFc,1,is)
            call intcell(nFc,nCv,mpg,mpg%intcellR,
     .              fna_tmp(1:nFc,1,is),vvc(1:nCv,is))
          endif
        enddo
      endif

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
            upb(ic,nflai)=max(min(pl%ua(ic,is),switch%eir_ua_max),
     1                                           switch%eir_ua_min)
            uudiab(ic,nflai)=max(min(
     1         wadiac(ic,0,is)+vaecrbc(ic,0,is),switch%eir_ua_max),
     2                                            switch%eir_ua_min)
            vvdiab(ic,nflai)=max(min(
     1         wadiac(ic,1,is)+vaecrbc(ic,1,is),switch%eir_ua_max),
     2                                            switch%eir_ua_min)
            if(.not.is_neutral(is)) then
              prb(ic)=prb(ic)+tib(ic)*dnib(ic,nflai)
            else
              prb(ic)=prb(ic)+tnb(ic)*dnib(ic,nflai)
            endif
c     Velocities uub, vvb, wwb
c     * cell centered
c     * internal cell values used for plasma background properties in infcop
c     * guard cell values used for source sampling in infcop
            if(switch%eirene_uub_style.eq.0) then
              if (switch%eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only
                uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               vadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               switch%eir_ua_max),switch%eir_ua_min)
              else
c     internal cell -> full diamagnetic drift
                uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               wadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               switch%eir_ua_max),switch%eir_ua_min)
                 endif
            elseif(switch%eirene_uub_style.eq.1) then
cwdk  uub computed from face fluxes, and interpolated to cell centers
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at N/S boundaries questionable at present, because (certain/all) terms
c     in fna may have been set 0. Should consider using the eirene_uub_style.eq.0 treatment 
c     in these cells!
              uub(ic,nflai)=max(min(uuc(ic,is),switch%eir_ua_max),
     1                                           switch%eir_ua_min)
            end if

cwdk  vvb computed from face fluxes, and interpolated to cell centers
c     Compatible with SOLPS4.3.
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at E/W boundaries questionable, because (certain/all) terms
c     in fna may have been set 0. 
            vvb(ic,nflai)=max(min(vvc(ic,is),switch%eir_ua_max),
     1                                         switch%eir_ua_min)
cwdk  wwb directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but may not be fully consistent with respect to drifts.
c     In particular: convective anomalous contribution  (if present) is missing
c                    (diffusive anomalous contribution zero due to toroidal symm.)
              if (switch%eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                wwb(ic,nflai)=max(min(
     .            upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     2            (vadiac(ic,0,is)+vaecrbc(ic,0,is))*
     3            gm%cvBb(ic,0)/gm%cvBb(ic,2),
     4            switch%eir_ua_max),switch%eir_ua_min)
              else
c     internal cell -> full diamagnetic drift
                 wwb(ic,nflai)=max(min(
     1             upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     2             (wadiac(ic,0,is)+vaecrbc(ic,0,is))*
     3             gm%cvBb(ic,0)/gm%cvBb(ic,2),
     4             switch%eir_ua_max),switch%eir_ua_min)
              end if
              
              ebx = gm%cvEb(ic,0)
              eby = gm%cvEb(ic,1)
              ebn = sqrt(ebx*ebx + eby*eby)
              ebx = ebx / ebn
              eby = eby / ebn
              pux(ic) = ebx
              puy(ic) = eby
              pvx(ic) = -eby
              pvy(ic) = ebx
              vxb(ic,nflai) = uub(ic,nflai)*pux(ic) +
     1                        vvb(ic,nflai)*pvx(ic)
              vyb(ic,nflai) = uub(ic,nflai)*puy(ic) +
     1                        vvb(ic,nflai)*pvy(ic)
              vzb(ic,nflai) = wwb(ic,nflai)*sign(1.0_R8,gm%cvEb(ic,2))

          enddo
* cell faces
          do ifc = 1, nFc
            fnib(ifc,nflai) = dv%fna(ifc,0,is) + dv%fna(ifc,1,is)
            if (is_neutral(is)) then
              fnn_inc_b(ifc,nflai) = fnn_inc(ifc,is)
              fnn_mol_b(ifc,nflai) = fna_mol(ifc,is)
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
     1   ' not avaible as an eirene_extrap option'
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
        if (ien-ian+1 <= 0) cycle
        
#ifdef B25_EIRENE
        call eirene_alloc_target_data(istra,ntrgdat,nflai,ifoff)
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
     &                mpg%rcFcOr(ian+i-1)*kin_frac_hyb(iFc)
            elseif (zab(ifl).eq.0.and.mol(istra).eq.0) then !! neutral atom species
              trgt(istra)%flux(i,ifl) = fnn_inc_b(iFc,ifl)*
     &                mpg%rcFcOr(ian+i-1)*kin_frac_hyb(iFc)
            elseif (zab(ifl).eq.0.and.mol(istra).eq.1) then     !! neutral molecule species
              trgt(istra)%flux(i,ifl) = 2.0_R8*fnn_mol_b(iFc,ifl)*
     &                                  fluid_frac_hyb(iFc)
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
            ebx = ebx / ebn
            eby = eby / ebn
            puxs = ebx
            puys = eby
            pvxs = -eby
            pvys = ebx
            uu = uub(ic,ifl)
            vv = vvb(ic,ifl)
            ww = wwb(ic,ifl)
            trgt(istra)%vxt(i,ifl) = puxs*uu + pvxs*vv
            trgt(istra)%vyt(i,ifl) = puys*uu + pvys*vv
            trgt(istra)%vzt(i,ifl) = ww
            trgt(istra)%fel(i,ifl) = 0.0_R8
            trgt(istra)%vpart(i,ifl) = 0.0_R8
            trgt(istra)%mach(i,ifl) = 1.0_R8
            trgt(istra)%usr(i,ifl) = 0.0_R8
          end do
        end do
      end do

      ncall=ncall+1
      call subend ()
      return
      end subroutine copy_background

      subroutine write_f30
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
      real*8, allocatable :: dummy(:,:,:)
      integer :: i, ic, iFc, ifl, ntrgdat, istra

!  .. use approach to generalized user interface

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

      write (31,'(2i6)') mpg%nCv, mpg%nCi, mpg%nFc
      
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
     .        trgt(istra)%fe(i),trgt(istra)%flength(i)
        end do

        do ifl = 1, nflai

          write (31,'(a)') '*'
          write (31,'(a,i0,a,i0)') '*** ION #',ifl,
     .                             ' TARGET DATA #',istra
          write (31,'(a)') '*'

          write (31,'(a1,a6,11a15)') '*','Index','Flux','Ti','ni',
     .      'Vx','Vy','Vz','fi','fel','vpart','mach','usr'
          write (31,'(a1,a6,11a15)') '*',' ',' ','(eV)','(1/m^3)',
     .      'm/s','m/s','m/s',' ',' ',' ',' '

          write (31,'(3es15.7)') amb(ifl), znb(ifl), zab(ifl)
          write (31,'(i6)') ntrgdat
 
          do i = 1, ntrgdat
            write (31,'(i7,11es15.7)')
     .        i,trgt(istra)%flux(i,ifl),trgt(istra)%tit(i),
     .          trgt(istra)%dnit(i,ifl),trgt(istra)%vxt(i,ifl),
     .          trgt(istra)%vyt(i,ifl),trgt(istra)%vzt(i,ifl),
     .          trgt(istra)%fi(i,ifl),trgt(istra)%fel(i,ifl),
     .          trgt(istra)%vpart(i,ifl),trgt(istra)%mach(i,ifl),
     .          trgt(istra)%usr(i,ifl)
          end do
        end do ! ifl
      end do   ! istra
      return
      end subroutine write_f31
      
