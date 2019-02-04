*//GFSUB3//
*=======================================================================
*                 S U B R O U T I N E  G F S U B 3
*=======================================================================
*
*
*  1. PURPOSE
*
*     GFSUB3 WRITES AN ARRAY ONTO AN OUTPUT FILE
*
*
*=======================================================================
*//DECLARATIONS//

      SUBROUTINE GFSUB3(KARD,NX,NY,NDIMX,NDIMY,NDIMS,DUMMY)
      use b2mod_types
      IMPLICIT NONE
      INTEGER KARD,NX,NY,NDIMX,NDIMY, NDIMS, ND1, LIM, IS, IY, IX, III
      REAL (kind=R8) ::
     *   DUMMY(0:NDIMX+1,0:NDIMY+1,1:NDIMS)
      ND1 = NX + 2
      LIM = (ND1/5)*5 - 4
      DO  120  IS = 1,NDIMS
        DO  110  IY = 0,NY+1
          DO  100  IX = 1,LIM,5
            WRITE(KARD,910) (DUMMY(IX-2+III,IY,IS),III = 1,5)
  100     CONTINUE
          IF( (LIM+4).EQ.ND1 )   GOTO 110
          WRITE(KARD,910) (DUMMY(IX-1,IY,IS),IX = LIM+5,ND1)
  110   CONTINUE
  120 CONTINUE
      RETURN
  910 FORMAT(1P,5(E16.8))
*=======================================================================
*//END GFSUB3//

      END SUBROUTINE GFSUB3
c
      subroutine copy_background(nCv, nFc, ns, mpg, gm, pl, dv, BoRiS)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_constants
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
      integer nCv, nFc, ns
      real (kind=R8) :: BoRiS
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2Plasma), intent(in) :: pl
      type(B2Derivatives), intent(in) :: dv

      character chns*3                                                 !srv 16.09.02
      integer, save :: drift_style = 1, iout = 0                        !srv 21.06.08
      integer ix,iy,is,ncall,eirene_extrap,in,
     1 eirene_uub_style, print_minmax, eirene_background
      external ipgetr, ipgeti, subini, subend, xerrab
      real (kind=R8) ::
     * te_min, te_max, ti_min, ti_max, na_min, na_max, ua_min, ua_max,
     * te_min_val, te_max_val, ti_min_val, ti_max_val,
     * na_min_val, na_max_val, ua_min_val, ua_max_val,
     * ne_min_val, ne_max_val, ne_min, ne_max
      save te_min, te_max, ti_min, ti_max,
     1 na_min, na_max, ua_min, ua_max, ncall, ne_min, ne_max,
     2 eirene_extrap, eirene_uub_style, print_minmax
      data te_min/0.0_R8/, te_max/1.0e30_R8/,
     1 ti_min/0.0_R8/, ti_max/1.0e30_R8/,
     1 na_min/0.0_R8/, na_max/1.0e30_R8/,
     1 ne_min/0.0_R8/, ne_max/1.0e30_R8/,
     2 ncall/0/, eirene_extrap/1/,
     & eirene_uub_style/0/, print_minmax/0/, eirene_background/1/
      real (kind=R8) :: kintmp(nFc)
      real (kind=R8) :: 
     &  vaecrbc(nCv,0:1,0:ns-1), wadiac(nCv,0:1,0:ns-1),
     &  vadiac(nCv,0:1,0:ns-1), wrk(nCv)
      real (kind=R8) :: naf(nFc,0:ns-1), fna_tmp(nFc,0:1,0:ns-1),
     &  uuc(nCv,0:ns-1),vvc(nCv,0:ns-1),nafsqalf(nFc,0:1,0:ns-1)
      integer :: ic, ifc, ian, ien, ntrgdat, ifl, i, istra
      real (kind=R8) :: intface_s
      real (kind=R8) :: evi

      call subini ('copy_background')
      if(ncall.eq.0) then
        ua_min = -c
        ua_max = +c
        call ipgetr('eirene_te_min',te_min)
        call ipgetr('eirene_te_max',te_max)
        call ipgetr('eirene_ti_min',ti_min)
        call ipgetr('eirene_ti_max',ti_max)
        call ipgetr('eirene_na_min',na_min)
        call ipgetr('eirene_na_max',na_max)
        call ipgetr('eirene_ne_min',ne_min)
        call ipgetr('eirene_ne_max',ne_max)
        call ipgetr('eirene_ua_min',ua_min)
        call ipgetr('eirene_ua_max',ua_max)
        call ipgeti('eirene_extrap',eirene_extrap)
        call ipgeti('eirene_uub_style',eirene_uub_style)
        call ipgeti('eirene_background',eirene_background)
        call ipgeti('eirene_print_minmax',print_minmax)
        call ipgeti('b2tfnb_drift_style', drift_style)                 !srv 21.06.08
        call ipgeti('copy_background_iout', iout)                      !srv 14.07.10
      endif 

      evi = 1.0_R8 / ev
      if(print_minmax.ne.0) then
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
      if (drift_style.eq.1) then
        do is=0,ns-1
          call intcell(nFc,nCv,mpg,mpg%intcellP,
     &                 dv%vaecrb(:,0,is),vaecrbc(:,0,is))
          call intcell(nFc,nCv,mpg,mpg%intcellP,
     &                 dv%wadia(:,0,is),wadiac(:,0,is))
          call intcell(nFc,nCv,mpg,mpg%intcellP,
     &                 dv%vadia(:,0,is),vadiac(:,0,is))
          call intcell(nFc,nCv,mpg,mpg%intcellR,
     &                 dv%vaecrb(:,1,is),vaecrbc(:,1,is))
          call intcell(nFc,nCv,mpg,mpg%intcellR,
     &                 dv%wadia(:,1,is),wadiac(:,1,is))
          call intcell(nFc,nCv,mpg,mpg%intcellR,
     &                 dv%vadia(:,1,is),vadiac(:,1,is))
        enddo
      else
        call xerrab ('copy_background -- drift_style .ne. 1')
      endif

      if(iout.ne.0) then                                               !srv 18.06.08 {
        do is = 0,ns-1
         write (chns,'(i3.3)') is
         call my_out_us(70,nFc,dv%fna(1,0,is),'copy_fnax'//chns)
         call my_out_us(70,nFc,dv%fna(1,1,is),'copy_fnay'//chns)
         call my_out_us(70,nFc,dv%fna_eir(1,0,is),'copy_fna_eirx'//chns)
         call my_out_us(70,nFc,dv%fna_eir(1,1,is),'copy_fna_eiry'//chns)
         call my_out_us(70,nFc,dv%vaecrb(1,0,is),'copy_vecrbx'//chns)
         call my_out_us(70,nFc,dv%vaecrb(1,1,is),'copy_vecrby'//chns)
         call my_out_us(70,nFc,dv%wadia(1,0,is),'copy_wadiax'//chns)
         call my_out_us(70,nFc,dv%wadia(1,1,is),'copy_wadiay'//chns)
         call my_out_us(70,nFc,dv%vadia(1,0,is),'copy_vadiax'//chns)
         call my_out_us(70,nFc,dv%vadia(1,1,is),'copy_vadiay'//chns)
         call my_out_us(70,nCv,vaecrbc(1,0,is),'copy_vecrbcx'//chns)
         call my_out_us(70,nCv,vaecrbc(1,1,is),'copy_vecrbcy'//chns)
         call my_out_us(70,nCv,wadiac(1,0,is),'copy_wadiacx'//chns)
         call my_out_us(70,nCv,wadiac(1,1,is),'copy_wadiacy'//chns)
         call my_out_us(70,nCv,vadiac(1,0,is),'copy_vadiacx'//chns)
         call my_out_us(70,nCv,vadiac(1,1,is),'copy_vadiacy'//chns)
         call my_out_us(70,nCv,pl%ua(1,is),'copy_ua'//chns)
         call my_out_us(70,nCv,pl%na(1,is),'copy_na'//chns)
        enddo
        call my_out_us(70,nFc,dv%fne_eir(1,0),'copy_fne_eirx')
        call my_out_us(70,nFc,dv%fne_eir(1,1),'copy_fne_eiry')
        call my_out_us(70,nFc,dv%fhe_eir(1,0),'copy_fhe_eirx')
        call my_out_us(70,nFc,dv%fhe_eir(1,1),'copy_fhe_eiry')
        call my_out_us(70,nFc,dv%fhi_eir(1,0),'copy_fhi_eirx')
        call my_out_us(70,nFc,dv%fhi_eir(1,1),'copy_fhi_eiry')
        call my_out_us(70,nCv,dv%ne,'copy_ne')
        call my_out_us(70,nCv,pl%te,'copy_te')
        call my_out_us(70,nCv,pl%ti,'copy_ti')
        call my_out_us(70,nCv,pl%po,'copy_po')
      endif 
* cell centred
      do ic=1,nCv
        neb(ic)=max(min(dv%ne(ic),ne_max),ne_min)
        teb(ic)=max(min(pl%te(ic)*evi,te_max),te_min)
        tib(ic)=max(min(pl%ti(ic)*evi,ti_max),ti_min)
        pob(ic)=pl%po(ic)
        prb(ic)=teb(ic)*max(min(dv%ne(ic),na_max),na_min)
        rrb(ic)=abs(gm%cvBb(ic,0)/gm%cvBb(ic,3))
        volb(ic)=gm%cvVol(ic)
!        bpolb(ic)=gm%cvBb(ic,0)
!        bradb(ic)=gm%cvBb(ic,1)
!        btorb(ic)=gm%cvBb(ic,2)
!        bfeldb(ic)=gm%cvBb(ic,3)
        bbxb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,0)
        bbyb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,1)
        bbzb(ic)=gm%cvBb(ic,3) * gm%cvEb(ic,2)
        bfeldb(ic)=gm%cvBb(ic,3)
      enddo

* faces
      call intface(nCv,nFc,mpg%fcCv,gm%fcVol,dv%kinrgy,kintmp)
      do ifc = 1, nFc
        feib(ifc) = dv%fhi(ifc,0) + dv%fhi(ifc,1)
        feeb(ifc) = dv%fhe(ifc,0) + dv%fhe(ifc,1)
        do is = 0, ns-1
          feib(ifc) = feib(ifc) + kintmp(ifc)*(1.0_R8-BoRiS)
     .                *(dv%fna(ifc,0,is) + dv%fna(ifc,1,is))
        end do
      end do
#ifdef WG_TODO
* cell face x
      do iy=-1,ny
        do ix=0,nx
          if(rightix(ix-1,iy).ne.nx+1) then
            feixb(ix,iy+1)=fhi(rightix(ix-1,iy),rightiy(ix-1,iy),0,0)+
     .                     fhi(rightix(ix-1,iy),rightiy(ix-1,iy),0,1)
            do is=0,ns-1
              kintmp =
     .         (kinrgy(rightix(ix-1,iy),rightiy(ix-1,iy),is)*
     .          hx(ix-1,iy)+
     .          kinrgy(ix-1,iy,is)*
     .          hx(rightix(ix-1,iy),rightiy(ix-1,iy)))/
     .         (hx(ix-1,iy)+hx(rightix(ix-1,iy),rightiy(ix-1,iy)))
              feixb(ix,iy+1)=feixb(ix,iy+1)+kintmp*(1.0_R8-BoRiS)
     .         *(fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,0,is) + 
     .           fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,1,is))
            enddo
            feexb(ix,iy+1)=fhe(rightix(ix-1,iy),rightiy(ix-1,iy),0,0)+
     .                     fhe(rightix(ix-1,iy),rightiy(ix-1,iy),0,1)
          else
            feixb(ix,iy+1)=0.0_R8
            feexb(ix,iy+1)=0.0_R8
          endif
        enddo
      enddo
      feixb(nx+1,:)=0.0_R8
      feexb(nx+1,:)=0.0_R8
* cell face y
      do iy=0,ny
        do ix=-1,nx
          feiyb(ix+1,iy)=fhi(ix,iy,1,0)+fhi(ix,iy,1,1)
          do is=0,ns-1
            if (bottomiy(ix,iy).ne.-2) then
              kintmp =
     .         (kinrgy(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     .          kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
     .         (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            else
              kintmp = kinrgy(ix,iy,is)
            endif
            feiyb(ix+1,iy)=feiyb(ix+1,iy)+kintmp*(1.0_R8-BoRiS)
     .       *(fna(ix,iy,1,0,is)+fna(ix,iy,1,1,is))
          enddo
          feeyb(ix+1,iy)=fhe(ix,iy,1,0)+fhe(ix,iy,1,1)
        enddo
      enddo
      feiyb(:,ny+1)=0.0_R8
      feeyb(:,ny+1)=0.0_R8
#endif

* by ns
      if (eirene_uub_style.eq.1) then
*   ..interpolate density to cell faces
        call intface1(nCv,nFc,ns-1,mpg%fcCv,gm%fcVol,pl%na,naf)
        do is = 0,ns-1
          nafsqalf(1:nFc,0,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,0)*
     .                           naf(1:nFc,is)+1.E-30_R8
          nafsqalf(1:nFc,1,is) = gm%fcS(1:nFc)*gm%fcQalf(1:nFc,1)*
     .                           naf(1:nFc,is)+1.E-30_R8
        end do
        if (eirene_background.eq.0) then
!         fna_tmp(1:nFc,0,0:ns-1) = dv%fna(1:nFc,0,0:ns-1)/
!    .      (gm%fcS(1:nFc)*gm%fcQalf(1:nFc,0)*naf(1:nFc,0:ns-1)+
!    .       1.E-30_R8)
          fna_tmp(1:nFc,0,0:ns-1) = dv%fna(1:nFc,0,0:ns-1)/
     .                              nafsqalf(1:nFc,0,0:ns-1)
!         fna_tmp(1:nFc,1,0:ns-1) = dv%fna(1:nFc,1,0:ns-1)/
!    .      (gm%fcS(1:nFc)*gm%fcQalf(1:nFc,1)*naf(1:nFc,0:ns-1)+
!    .      1.E-30_R8)
          fna_tmp(1:nFc,1,0:ns-1) = dv%fna(1:nFc,1,0:ns-1)/
     .                              nafsqalf(1:nFc,1,0:ns-1)
        else
!         fna_tmp(1:nFc,0,0:ns-1) = dv%fna_eir(1:nFc,0,0:ns-1)/
!    .      (gm%fcS(1:nFc)*gm%fcQalf(1:nFc,0)*naf(1:nFc,0:ns-1)+
!    .      1.E-30_R8)
          fna_tmp(1:nFc,0,0:ns-1) = dv%fna_eir(1:nFc,0,0:ns-1)/
     .                              nafsqalf(1:nFc,0,0:ns-1)
!         fna_tmp(1:nFc,1,0:ns-1) = dv%fna_eir(1:nFc,1,0:ns-1)/
!    .      (gm%fcS(1:nFc)*gm%fcQalf(1:nFc,1)*naf(1:nFc,0:ns-1)+
!    .      1.E-30_R8)
          fna_tmp(1:nFc,1,0:ns-1) = dv%fna_eir(1:nFc,1,0:ns-1)/
     .                              nafsqalf(1:nFc,1,0:ns-1)
        endif
        call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellP,
     .                fna_tmp(:,0,:),uuc)
        call intcell1(nFc,nCv,ns-1,mpg,mpg%intcellR,
     .                fna_tmp(:,1,:),vvc)
      end if

      nflai=0
      do is=0,ns-1
        if(.not.is_neutral(is)) then
          nflai=nflai+1
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
            dnib(ic,nflai)=max(min(pl%na(ic,is),na_max),na_min)
            upb(ic,nflai)=max(min(pl%ua(ic,is),ua_max),ua_min)
            uudiab(ic,nflai)=max(min(
     1         wadiac(ic,0,is)+vaecrbc(ic,0,is),ua_max),ua_min)
            vvdiab(ic,nflai)=max(min(
     1         wadiac(ic,1,is)+vaecrbc(ic,1,is),ua_max),ua_min)
            prb(ic)=prb(ic)+tib(ic)*dnib(ic,nflai)
c     Velocities uub, vvb, wwb
c     * cell centered
c     * internal cell values used for plasma background properties in infcop
c     * guard cell values used for source sampling in infcop
            if(eirene_uub_style.eq.0) then
              if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only
                uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               vadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               ua_max),ua_min)
              else
c     internal cell -> full diamagnetic drift
                uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               wadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               ua_max),ua_min)
                 endif
            elseif(eirene_uub_style.eq.1) then
cwdk  uub computed from face fluxes, and interpolated to cell centers
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at N/S boundaries questionable at present, because (certain/all) terms
c     in fna may have been set 0. Should consider using the eirene_uub_style.eq.0 treatment 
c     in these cells!
              uub(ic,nflai)=max(min(uuc(ic,is),ua_max),ua_min)
            end if
#ifdef WG_TODO
              if(eirene_uub_style.eq.0) then
cwdk  uub directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but not fully consistent with respect to anomalous transport.
c     In particular: possible anomalous contributions to parallel and 
c     diamagnetic velocities are missing (case b2trno_pol_anom_scale.ne.0.0)
               if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .           bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only
                 uub(ic,nflai)=max(min(
     1             upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2             vadiac(ic,0,is)+vaecrbc(ic,0,is),
     3             ua_max),ua_min)
               else
                 if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only
                   uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               vadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               ua_max),ua_min)
                 else
c     internal cell -> full diamagnetic drift
                   uub(ic,nflai)=max(min(
     1               upb(ic,nflai)*gm%cvBb(ic,0)/gm%cvBb(ic,3)+
     2               wadiac(ic,0,is)+vaecrbc(ic,0,is),
     3               ua_max),ua_min)
                 endif
               endif
              elseif(eirene_uub_style.eq.1) then
cwdk  uub computed from face fluxes, and interpolated to cell centers
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at N/S boundaries questionable at present, because (certain/all) terms
c     in fna may have been set 0. Should consider using the eirene_uub_style.eq.0 treatment 
c     in these cells!
!              if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
!    .           bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
               if (ic > nCi) then
c     guard cell -> use grad-B drift only (fna)
                 if (leftix(ix,iy).eq.-2) then
                   uub(ix+1,iy+1,nflai)=max(min(
     1               fna(rightix(ix,iy),rightiy(ix,iy),0,0,is)/(
     3               gs(rightix(ix,iy),rightiy(ix,iy),0)*
     4               na(ix,iy,is)*qc(rightix(ix,iy),rightiy(ix,iy),0)),
     5               ua_max),ua_min)
                 elseif (rightix(ix,iy).eq.nx+1) then
                   uub(ix+1,iy+1,nflai)=max(min(
     1               fna(ix,iy,0,0,is)/(gs(ix,iy,0)*
     4               na(ix,iy,is)*qc(ix,iy,0)),
     5               ua_max),ua_min)
                 else
                  uub(ix+1,iy+1,nflai)=max(min(
     1             fna(ix,iy,0,0,is)/(gs(ix,iy,0)*qc(ix,iy,0)*
     2             (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3             fna(rightix(ix,iy),rightiy(ix,iy),0,0,is)/(
     4             gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5             qc(rightix(ix,iy),rightiy(ix,iy),0)*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7             ua_max),ua_min)
                 endif
               else
                 if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                  uub(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,0,0,is)/(gs(ix,iy,0)*qc(ix,iy,0)*
     2              (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3              fna(rightix(ix,iy),rightiy(ix,iy),0,0,is)/(
     4              gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5              qc(rightix(ix,iy),rightiy(ix,iy),0)*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7              ua_max),ua_min)
                 else
c     internal cell -> full diamagnetic drift (fna_eir)
                  uub(ix+1,iy+1,nflai)=max(min(
     1              fna_eir(ix,iy,0,0,is)/(gs(ix,iy,0)*qc(ix,iy,0)*
     2              (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3              fna_eir(rightix(ix,iy),rightiy(ix,iy),0,0,is)/(
     4              gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5              qc(rightix(ix,iy),rightiy(ix,iy),0)*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7              ua_max),ua_min)
                 endif
               endif
              endif
#endif
cwdk  vvb computed from face fluxes, and interpolated to cell centers
c     Compatible with SOLPS4.3.
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at E/W boundaries questionable, because (certain/all) terms
c     in fna may have been set 0. 
            vvb(ic,nflai)=max(min(vvc(ic,is),ua_max),ua_min)
#ifdef WG_TODO
              if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .          bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only (fna)
                if (bottomiy(ix,iy).eq.-2) then
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(topix(ix,iy),topiy(ix,iy),1,1,is)/(
     2              gs(topix(ix,iy),topiy(ix,iy),1)*na(ix,iy,is)),
     3              ua_max),ua_min)
                elseif (topiy(ix,iy).eq.ny+1) then
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,1,is)/(gs(ix,iy,1)*na(ix,iy,is)),
     5              ua_max),ua_min)
                else
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3               na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna(topix(ix,iy),topiy(ix,iy),1,1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                endif
              else
                if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3              na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna(topix(ix,iy),topiy(ix,iy),1,1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                else
c     internal cell -> full diamagnetic drift (fna_eir)
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna_eir(ix,iy,1,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3              na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna_eir(topix(ix,iy),topiy(ix,iy),1,1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                endif
              endif
#endif
cwdk  wwb directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but may not be fully consistent with respect to drifts.
c     In particular: convective anomalous contribution  (if present) is missing
c                    (diffusive anomalous contribution zero due to toroidal symm.)
              if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                wwb(ic,nflai)=max(min(
     .            upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     2            (vadiac(ic,0,is)+vaecrbc(ic,0,is))*
     3            gm%cvBb(ic,0)/gm%cvBb(ic,2),
     4            ua_max),ua_min)
              else
c     internal cell -> full diamagnetic drift
                 wwb(ic,nflai)=max(min(
     1             upb(ic,nflai)*gm%cvBb(ic,2)/gm%cvBb(ic,3)-
     2             (wadiac(ic,0,is)+vaecrbc(ic,0,is))*
     3             gm%cvBb(ic,0)/gm%cvBb(ic,2),
     4             ua_max),ua_min)
              end if
              
              vxb(ic,nflai) = gm%cvEb(ic,0)*(uub(ic,nflai) +
     1                        vvb(ic,nflai))
              vyb(ic,nflai) = gm%cvEb(ic,1)*(uub(ic,nflai) +
     1                        vvb(ic,nflai))
              vzb(ic,nflai) = wwb(ic,nflai)

#ifdef WG_TODO
              if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .          bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only
               wwb(ix+1,iy+1,nflai)=max(min(
     1           upb(ix+1,iy+1,nflai)*bb(ix,iy,2)/bb(ix,iy,3)-
     2           (vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is))*
     3           bb(ix,iy,0)/bb(ix,iy,2),
     4           ua_max),ua_min)
              else
                if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                 wwb(ix+1,iy+1,nflai)=max(min(
     1             upb(ix+1,iy+1,nflai)*bb(ix,iy,2)/bb(ix,iy,3)-
     2             (vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is))*
     3             bb(ix,iy,0)/bb(ix,iy,2),
     4             ua_max),ua_min)
                else
c     internal cell -> full diamagnetic drift
                 wwb(ix+1,iy+1,nflai)=max(min(
     1             upb(ix+1,iy+1,nflai)*bb(ix,iy,2)/bb(ix,iy,3)-
     2             (wadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is))*
     3             bb(ix,iy,0)/bb(ix,iy,2),
     4             ua_max),ua_min)
                endif
              endif
#endif
          enddo
* cell faces
          do ifc = 1, nFc
            fnib(ifc,nflai) = dv%fna(ifc,0,is) + dv%fna(ifc,1,is)
          end do
#ifdef WG_TODO
* cell face x
          do iy=-1,ny
            do ix=0,nx
              if(rightix(ix-1,iy).ne.nx+1) then
               if(leftix(ix,iy).ne.-2) then
                  fnixb(ix,iy+1,nflai)=
     1             fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,0,is)+
     2             fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,1,is)
               else
                 fnixb(ix,iy+1,nflai)=0.0_R8
               endif
              else
               fnixb(ix,iy+1,nflai)=0.0_R8
              endif
            enddo
          enddo
          fnixb(nx+1,:,nflai)=0.0_R8
* cell face y
          do iy=0,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).ne.-2) then
                fniyb(ix+1,iy,nflai)=fna(ix,iy,0,1,is)+fna(ix,iy,1,1,is)
              else
                fniyb(ix+1,iy,nflai)=0.0_R8
              endif
            enddo
          enddo
          fniyb(:,ny+1,nflai)=0.0_R8
        endif
#endif
        endif
      enddo

#ifdef WG_TODO
      if(eirene_extrap.eq.0) then
        ! do nothing
      elseif(eirene_extrap.eq.1) then
c we derive the guard cell value from its neighbouring interior cell
        do in=1,2    ! do twice so we get the corner cells as well
          do iy=-1,ny
            do ix=-1,nx
c the ix+1,iy+1 form arises because the eirene quantity arrays are declared
c 0:nx+1,0:ny+1 for compatibility with the old B2, and the B2.5 arrays are
c declared -1:nx,-1:ny
              if(leftix(ix,iy).eq.-2) then
                teb(ix+1,iy+1)=teb(rightix(ix,iy)+1,rightiy(ix,iy)+1)
                tib(ix+1,iy+1)=tib(rightix(ix,iy)+1,rightiy(ix,iy)+1)
                pob(ix+1,iy+1)=pob(rightix(ix,iy)+1,rightiy(ix,iy)+1)
                dnib(ix+1,iy+1,:)=
     1           dnib(rightix(ix,iy)+1,rightiy(ix,iy)+1,:)
              endif
              if(rightix(ix,iy).eq.nx+1) then
                teb(ix+1,iy+1)=teb(leftix(ix,iy)+1,leftiy(ix,iy)+1)
                tib(ix+1,iy+1)=tib(leftix(ix,iy)+1,leftiy(ix,iy)+1)
                pob(ix+1,iy+1)=pob(leftix(ix,iy)+1,leftiy(ix,iy)+1)
                dnib(ix+1,iy+1,:)=
     1           dnib(leftix(ix,iy)+1,leftiy(ix,iy)+1,:)
              endif
              if(bottomiy(ix,iy).eq.-2) then
                teb(ix+1,iy+1)=teb(topix(ix,iy)+1,topiy(ix,iy)+1)
                tib(ix+1,iy+1)=tib(topix(ix,iy)+1,topiy(ix,iy)+1)
                pob(ix+1,iy+1)=pob(topix(ix,iy)+1,topiy(ix,iy)+1)
                dnib(ix+1,iy+1,:)=
     1           dnib(topix(ix,iy)+1,topiy(ix,iy)+1,:)
              endif
              if(topiy(ix,iy).eq.ny+1) then
                teb(ix+1,iy+1)=teb(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
                tib(ix+1,iy+1)=tib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
                pob(ix+1,iy+1)=pob(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
                dnib(ix+1,iy+1,:)=
     1           dnib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1,:)
              endif
            enddo
          enddo
        enddo
      else
        write(*,*) eirene_extrap,
     1   ' not avaible as an eirene_extrap option'
      endif
#endif

      if(print_minmax.ne.0) then
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
      
      if(iout.ne.0) then
       nflai=0
       do is=0,ns-1
        if(.not.is_neutral(is)) then
         nflai=nflai+1
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,uub(1,nflai),'copy_uub'//chns)
         call my_out_us(70,nCv,vvb(1,nflai),'copy_vvb'//chns)
         call my_out_us(70,nCv,wwb(1,nflai),'copy_wwb'//chns)
         call my_out_us(70,nCv,upb(1,nflai),'copy_upb'//chns)
        endif
       enddo
       call my_out_us(70,nCv,teb,'copy_teb')
       call my_out_us(70,nCv,tib,'copy_tib')
       call my_out_us(70,nCv,pob,'copy_pob')
      endif 

!  .. define target data

      if (.not.allocated(trgt)) allocate(trgt(mpg%nRc))
      do istra = 1, mpg%nRc
        ian = mpg%rcFcP(istra,1)
        ien = ian + mpg%rcFcP(istra,2) - 1
        ntrgdat = ien-ian+1
        if (ien-ian+1 <= 0) cycle
        
        call alloc_target_data(istra,ntrgdat,nflai)

        do i = 1, ntrgdat
          iFc = mpg%rcFc(ian+i-1)
          trgt(istra)%faces(i) = iFc
          trgt(istra)%tet(i) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,pl%te)*evi
          trgt(istra)%tit(i) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,pl%ti)*evi
          trgt(istra)%fe(i) = 2.0_R8     ! to be specified correctly
          trgt(istra)%fsh(i) = delta_sheathb(iFc)
          do ifl = 1, nflai
            trgt(istra)%flux(i,ifl) = fnib(iFc,ifl)*mpg%rcFcOr(ian+i-1)
            trgt(istra)%dnit(i,ifl) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,dnib(:,ifl))
            trgt(istra)%vxt(i,ifl) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,vxb(:,ifl))
            trgt(istra)%vyt(i,ifl) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,vyb(:,ifl))
            trgt(istra)%vzt(i,ifl) = intface_s(iFc,nCv,nFc,mpg%fcCv,
     >                                   gm%fcVol,vzb(:,ifl))
            trgt(istra)%fi(i,ifl) = 2.5_R8
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
      write (31,'(a1,a6,11a15)') '*','Index','Te','ne',
     .      'Bx','By','Bz','Ex','Ey','Ez','pot','psi','f'
      write (31,'(a1,a6,11a15)') '*',' ','(eV)','(1/m^3)',
     .      'Tesla','Tesla','Tesla','(V/m)','(V/m)','(V/m)',
     .      'V','??','??'

      write (31,'(2i6)') mpg%nCi, mpg%nFc
      
      do ic = 1, mpg%nCi
        write (31,'(i7,11es15.7)') ic, teb(ic), neb(ic), 
     .      bbxb(ic), bbyb(ic), bbzb(ic), 0.0_R8, 0.0_R8, 0.0_R8,
     .      pob(ic), 0.0_R8, 0.0_R8
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
        write (31,'(a1,2a6,4a15)') '*','Index','face','Te',
     .      'fe','fsh'
        write (31,'(a1,2a6,5a15)') '*',' ',' ','(eV)',' ',
     .      ' ',' '
        write (31,'(i6)') ntrgdat
 
        do i = 1, ntrgdat
          write (31,'(i7,i6,3es15.7)')
     .      i,trgt(istra)%faces(i),trgt(istra)%tet(i),
     .        trgt(istra)%fe(i),trgt(istra)%fsh(i)
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
        

      
      
#ifdef WG_TODO
      write(*,*) 'NX, NY, NFLAI = ',nx,ny,nflai
      rewind(31)
c ion density
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,DNIB(0,0,1))
c poloidal velocity
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,UUB(0,0,1))
c radial velocity
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,VVB(0,0,1))
c toroidal velocity
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,WWB(0,0,1))
c electron temperature
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,TEB(0,0))
c ion temperature
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,TIB(0,0))
c pressure
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,PRB(0,0))
c parallel velocity (ua(,,,))
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,UPB(0,0,1))
c pitch angle (bx/btot)
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,RRB(0,0))
c poloidal ion flow on the left face
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,FNIXB(0,0,1))
c radial ion flow on the bottom face
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,FNIYB(0,0,1))
c poloidal ion heat flux on the left face
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,FEIXB(0,0))
c radial ion heat flux on the bottom face
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,FEIYB(0,0))
c poloidal electron heat flux on the left face
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,FEEXB(0,0))
c radial electron heat flux on the bottom face
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,FEEYB(0,0))
c total ion drift velocity in diamagnetic direction
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,UUDIAB(0,0,1))
c total ion drift velocity in radial direction
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,VVDIAB(0,0,1))
!pb preparation for future use
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,POB(0,0))
c cell volumes
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,VOLB(0,0))
c magnitude and components of the magnetic field
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,BFELDB(0,0))
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,BPOLB(0,0))
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,BRADB(0,0))
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,BTORB(0,0))
csw 13apr2011 dummies for vparx,...,deltae_...
      allocate(dummy(0:nxx+1,0:nyy+1,nflai))
      dummy=0.d0
      call gfsub3(31,nx,ny,nxx,nyy,nflai,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,nflai,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,nflai,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,nflai,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,delta_sheathxb)
      call gfsub3(31,nx,ny,nxx,nyy,1,delta_sheathyb)



      deallocate(dummy)

      write(*,*) 'Background plasma written'
      close(31)
#endif
      return
      end subroutine write_f31
