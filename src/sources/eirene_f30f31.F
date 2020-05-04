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
      subroutine copy_background(nx, ny, ns, BoRiS)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_parmmod
      use eirmod_braeir
      use eirmod_extrab25
#else
      use parmmod   ! IGNORE
      use braeir    ! IGNORE
      use eirmod_cinit, only: fort_lc
#endif
#else
      use b2mod_braeir
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/KOPPLDIM.F>
#endif
      integer nx,ny,ns
      real (kind=R8) :: BoRiS
      character chns*3                                                 !srv 16.09.02
      integer, save :: drift_style = 1, iout = 0                        !srv 21.06.08
      integer ix,iy,is,ncall,eirene_extrapolate,in,
     1 eirene_uub_style, print_minmax, eirene_background
      external ipgetr, ipgeti, subini, subend, xerrab
      real (kind=R8) ::
     * te_min, te_max, ti_min, ti_max, na_min, na_max, ua_min, ua_max,
     * te_min_val, te_max_val, ti_min_val, ti_max_val,
     * na_min_val, na_max_val, ua_min_val, ua_max_val
      save te_min, te_max, ti_min, ti_max,
     1 na_min, na_max, ua_min, ua_max, ncall,
     2 eirene_extrapolate, eirene_uub_style, print_minmax
      data te_min/0.0_R8/, te_max/1.0e30_R8/,
     1 ti_min/0.0_R8/, ti_max/1.0e30_R8/,
     1 na_min/0.0_R8/, na_max/1.0e30_R8/,
     2 ncall/0/, eirene_extrapolate/1/,
     & eirene_uub_style/0/, print_minmax/0/, eirene_background/1/
      real (kind=R8) :: vdrift,
     &  vaecrbc(-1:nx,-1:ny,0:1,0:ns-1), wadiac(-1:nx,-1:ny,0:1,0:ns-1),
     &  vadiac(-1:nx,-1:ny,0:1,0:ns-1), wrk(-1:nx,-1:ny)

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
        call ipgetr('eirene_ua_min',ua_min)
        call ipgetr('eirene_ua_max',ua_max)
        call ipgeti('eirene_extrap',eirene_extrapolate)
        call ipgeti('eirene_uub_style',eirene_uub_style)
        call ipgeti('eirene_background',eirene_background)
        call ipgeti('eirene_print_minmax',print_minmax)
        call ipgeti('b2tfnb_drift_style', drift_style)                 !srv 21.06.08
        call ipgeti('copy_background_iout', iout)                      !srv 14.07.10
      endif 

      call xertst (nx.lt.NDXP, 'eirene_mc: increase size of ndxp')
      call xertst (ny.lt.NDYP, 'eirene_mc: increase size of ndyp')
      if(print_minmax.ne.0) then
         te_min_val=minval(te(-1:nx,-1:ny))/ev
         te_max_val=maxval(te(-1:nx,-1:ny))/ev
         ti_min_val=minval(ti(-1:nx,-1:ny))/ev
         ti_max_val=maxval(ti(-1:nx,-1:ny))/ev
         na_min_val=minval(na(-1:nx,-1:ny,0:ns-1))
         na_max_val=maxval(na(-1:nx,-1:ny,0:ns-1))
         ua_min_val=minval(ua(-1:nx,-1:ny,0:ns-1))
         ua_max_val=maxval(ua(-1:nx,-1:ny,0:ns-1))
         write(*,'(a,1p,8(1x,g12.3))') 'minmax: te, ti, na, ua ',
     &     te_min_val, te_max_val, ti_min_val, ti_max_val,
     &     na_min_val, na_max_val, ua_min_val, ua_max_val
      endif

* get drift velocities in cell centers
      vaecrbc = vaecrb
      wadiac  = wadia
      vadiac  = vadia
      if (drift_style.eq.1) then
        do is=0,ns-1
          call mapx(nx,ny,1,vaecrbc(:,:,0,is),wrk)
          call mapy(nx,ny,1,vaecrbc(:,:,1,is),wrk)
          call mapx(nx,ny,1,wadiac(:,:,0,is),wrk)
          call mapy(nx,ny,1,wadiac(:,:,1,is),wrk)
          call mapx(nx,ny,1,vadiac(:,:,0,is),wrk)
          call mapy(nx,ny,1,vadiac(:,:,1,is),wrk)
        enddo
      endif

      if(iout.ne.0) then                                               !srv 18.06.08 {
        do is = 0,ns-1
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,fna(-1,-1,0,is),'copy_fnax'//chns)
         call my_out(70,nx,ny,fna(-1,-1,1,is),'copy_fnay'//chns)
         call my_out(70,nx,ny,fna_eir(-1,-1,0,is),'copy_fna_eirx'//chns)
         call my_out(70,nx,ny,fna_eir(-1,-1,1,is),'copy_fna_eiry'//chns)
         call my_out(70,nx,ny,vaecrb(-1,-1,0,is),'copy_vecrbx'//chns)
         call my_out(70,nx,ny,vaecrb(-1,-1,1,is),'copy_vecrby'//chns)
         call my_out(70,nx,ny,wadia(-1,-1,0,is),'copy_wadiax'//chns)
         call my_out(70,nx,ny,wadia(-1,-1,1,is),'copy_wadiay'//chns)
         call my_out(70,nx,ny,vadia(-1,-1,0,is),'copy_vadiax'//chns)
         call my_out(70,nx,ny,vadia(-1,-1,1,is),'copy_vadiay'//chns)
         call my_out(70,nx,ny,vaecrbc(-1,-1,0,is),'copy_vecrbcx'//chns)
         call my_out(70,nx,ny,vaecrbc(-1,-1,1,is),'copy_vecrbcy'//chns)
         call my_out(70,nx,ny,wadiac(-1,-1,0,is),'copy_wadiacx'//chns)
         call my_out(70,nx,ny,wadiac(-1,-1,1,is),'copy_wadiacy'//chns)
         call my_out(70,nx,ny,vadiac(-1,-1,0,is),'copy_vadiacx'//chns)
         call my_out(70,nx,ny,vadiac(-1,-1,1,is),'copy_vadiacy'//chns)
         call my_out(70,nx,ny,ua(-1,-1,is),'copy_ua'//chns)
         call my_out(70,nx,ny,na(-1,-1,is),'copy_na'//chns)
        enddo
        call my_out(70,nx,ny,fne_eir(-1,-1,0),'copy_fne_eirx')
        call my_out(70,nx,ny,fne_eir(-1,-1,1),'copy_fne_eiry')
        call my_out(70,nx,ny,fhe_eir(-1,-1,0),'copy_fhe_eirx')
        call my_out(70,nx,ny,fhe_eir(-1,-1,1),'copy_fhe_eiry')
        call my_out(70,nx,ny,fhi_eir(-1,-1,0),'copy_fhi_eirx')
        call my_out(70,nx,ny,fhi_eir(-1,-1,1),'copy_fhi_eiry')
        call my_out(70,nx,ny,ne,'copy_ne')
        call my_out(70,nx,ny,te,'copy_te')
        call my_out(70,nx,ny,ti,'copy_ti')
        call my_out(70,nx,ny,po,'copy_po')
      endif 

* cell centred
      do iy=-1,ny
        do ix=-1,nx
          teb(ix+1,iy+1)=max(min(te(ix,iy),te_max*ev),te_min*ev)
          tib(ix+1,iy+1)=max(min(ti(ix,iy),ti_max*ev),ti_min*ev)
          pob(ix+1,iy+1)=po(ix,iy)
          prb(ix+1,iy+1)=teb(ix+1,iy+1)*
     1     max(min(ne(ix,iy),na_max),na_min)
          rrb(ix+1,iy+1)=abs(bb(ix,iy,0)/bb(ix,iy,3))
          volb(ix+1,iy+1)=vol(ix,iy)
          bpolb(ix+1,iy+1)=bb(ix,iy,0)
          bradb(ix+1,iy+1)=bb(ix,iy,1)
          btorb(ix+1,iy+1)=bb(ix,iy,2)
          bfeldb(ix+1,iy+1)=bb(ix,iy,3)
        enddo
      enddo
* cell face x
      do iy=-1,ny
        do ix=0,nx
          if(rightix(ix-1,iy).ne.nx+1) then
            feixb(ix,iy+1)=fhi(rightix(ix-1,iy),rightiy(ix-1,iy),0)
            do is=0,ns-1
              feixb(ix,iy+1)=feixb(ix,iy+1)+(1.0_R8-BoRiS)
     .         *fhm(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)
            enddo
            feexb(ix,iy+1)=fhe(rightix(ix-1,iy),rightiy(ix-1,iy),0)
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
          feiyb(ix+1,iy)=fhi(ix,iy,1)
          do is=0,ns-1
            feiyb(ix+1,iy)=feiyb(ix+1,iy)+(1.0_R8-BoRiS)*fhm(ix,iy,1,is)
          enddo
          feeyb(ix+1,iy)=fhe(ix,iy,1)
        enddo
      enddo
      feiyb(:,ny+1)=0.0_R8
      feeyb(:,ny+1)=0.0_R8
* by ns
      nflai=0
      do is=0,ns-1
        if(.not.is_neutral(is)) then
          nflai=nflai+1
          call xertst (nflai.le.NFL, 'eirene_mc: increase size of nfl')
* cell centred
          do iy=-1,ny
            do ix=-1,nx
              dnib(ix+1,iy+1,nflai)=max(min(na(ix,iy,is),na_max),na_min)
              upb(ix+1,iy+1,nflai)=max(min(ua(ix,iy,is),ua_max),ua_min)
              uudiab(ix+1,iy+1,nflai)=max(min(
     1         wadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is),ua_max),ua_min)
              vvdiab(ix+1,iy+1,nflai)=max(min(
     1         wadiac(ix,iy,1,is)+vaecrbc(ix,iy,1,is),ua_max),ua_min)
              prb(ix+1,iy+1)=prb(ix+1,iy+1)+tib(ix+1,iy+1)*
     1         dnib(ix+1,iy+1,nflai)
c     Velocities uub, vvb, wwb
c     * cell centered
c     * internal cell values used for plasma background properties in infcop
c     * guard cell values used for source sampling in infcop
              if(eirene_uub_style.eq.0) then
cwdk  uub directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but not fully consistent with respect to anomalous transport.
c     In particular: possible anomalous contributions to parallel and 
c     diamagnetic velocities are missing (case b2trno_pol_anom_scale.ne.0.0)
               if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .           bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only
                 uub(ix+1,iy+1,nflai)=max(min(
     1             upb(ix+1,iy+1,nflai)*bb(ix,iy,0)/bb(ix,iy,3)+
     2             vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is),
     3             ua_max),ua_min)
               else
                 if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only
                   uub(ix+1,iy+1,nflai)=max(min(
     1               upb(ix+1,iy+1,nflai)*bb(ix,iy,0)/bb(ix,iy,3)+
     2               vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is),
     3               ua_max),ua_min)
                 else
c     internal cell -> full diamagnetic drift
                   uub(ix+1,iy+1,nflai)=max(min(
     1               upb(ix+1,iy+1,nflai)*bb(ix,iy,0)/bb(ix,iy,3)+
     2               wadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is),
     3               ua_max),ua_min)
                 endif
               endif
              elseif(eirene_uub_style.eq.1) then
cwdk  uub computed from face fluxes, and interpolated to cell centers
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at N/S boundaries questionable at present, because (certain/all) terms
c     in fna may have been set 0. Should consider using the eirene_uub_style.eq.0 treatment 
c     in these cells!
               if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .           bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only (fna)
                 if (leftix(ix,iy).eq.-2) then
                   uub(ix+1,iy+1,nflai)=max(min(
     1               fna(rightix(ix,iy),rightiy(ix,iy),0,is)/(
     3               gs(rightix(ix,iy),rightiy(ix,iy),0)*
     4               na(ix,iy,is)*qc(rightix(ix,iy),rightiy(ix,iy))),
     5               ua_max),ua_min)
                 elseif (rightix(ix,iy).eq.nx+1) then
                   uub(ix+1,iy+1,nflai)=max(min(
     1               fna(ix,iy,0,is)/(gs(ix,iy,0)*
     4               na(ix,iy,is)*qc(ix,iy)),
     5               ua_max),ua_min)
                 else
                  uub(ix+1,iy+1,nflai)=max(min(
     1             fna(ix,iy,0,is)/(gs(ix,iy,0)*qc(ix,iy)*
     2             (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3             fna(rightix(ix,iy),rightiy(ix,iy),0,is)/(
     4             gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5             qc(rightix(ix,iy),rightiy(ix,iy))*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7             ua_max),ua_min)
                 endif
               else
                 if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                  uub(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,0,is)/(gs(ix,iy,0)*qc(ix,iy)*
     2              (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3              fna(rightix(ix,iy),rightiy(ix,iy),0,is)/(
     4              gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5              qc(rightix(ix,iy),rightiy(ix,iy))*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7              ua_max),ua_min)
                 else
c     internal cell -> full diamagnetic drift (fna_eir)
                  uub(ix+1,iy+1,nflai)=max(min(
     1              fna_eir(ix,iy,0,is)/(gs(ix,iy,0)*qc(ix,iy)*
     2              (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)))+
     3              fna_eir(rightix(ix,iy),rightiy(ix,iy),0,is)/(
     4              gs(rightix(ix,iy),rightiy(ix,iy),0)*
     5              qc(rightix(ix,iy),rightiy(ix,iy))*
     6             (na(ix,iy,is)+na(rightix(ix,iy),rightiy(ix,iy),is))),
     7              ua_max),ua_min)
                 endif
               endif
              endif
cwdk  vvb computed from face fluxes, and interpolated to cell centers
c     Compatible with SOLPS4.3.
c     Should be fully compatible with anom. transport, including all contributions to the flow.
c     Note, computation at E/W boundaries questionable, because (certain/all) terms
c     in fna may have been set 0. 
              if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .          bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only (fna)
                if (bottomiy(ix,iy).eq.-2) then
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(topix(ix,iy),topiy(ix,iy),1,is)/(
     2              gs(topix(ix,iy),topiy(ix,iy),1)*na(ix,iy,is)),
     3              ua_max),ua_min)
                elseif (topiy(ix,iy).eq.ny+1) then
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,is)/(gs(ix,iy,1)*na(ix,iy,is)),
     5              ua_max),ua_min)
                else
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3               na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna(topix(ix,iy),topiy(ix,iy),1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                endif
              else
                if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna(ix,iy,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3              na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna(topix(ix,iy),topiy(ix,iy),1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                else
c     internal cell -> full diamagnetic drift (fna_eir)
                  vvb(ix+1,iy+1,nflai)=max(min(
     1              fna_eir(ix,iy,1,is)/(gs(ix,iy,1)*
     2              (na(ix,iy,is)+
     3              na(bottomix(ix,iy),bottomiy(ix,iy),is)))+
     4              fna_eir(topix(ix,iy),topiy(ix,iy),1,is)/
     5              (gs(topix(ix,iy),topiy(ix,iy),1)*
     6              (na(ix,iy,is)+na(topix(ix,iy),topiy(ix,iy),is))),
     7              ua_max),ua_min)
                endif
              endif
cwdk  wwb directly computed from cell-centered velocities.
c     Compatible with SOLPS4.3, and consistent with sign convention ua,
c     but may not be fully consistent with respect to drifts.
c     In particular: convective anomalous contribution  (if present) is missing
c                    (diffusive anomalous contribution zero due to toroidal symm.)
              if (leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     .          bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy)  .eq.ny+1) then
c     guard cell -> use grad-B drift only
                vdrift=vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is)
              else if (eirene_background.eq.0) then
c     internal cell, old treatment -> grad-B drift only (fna)
                vdrift=vadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is)
              else
c     internal cell -> full diamagnetic drift
                vdrift=wadiac(ix,iy,0,is)+vaecrbc(ix,iy,0,is)
              endif
              if (bb(ix,iy,2).ne.0.0_R8)
     >         vdrift=vdrift*bb(ix,iy,0)/bb(ix,iy,2)
              wwb(ix+1,iy+1,nflai)=max(min(
     1          upb(ix+1,iy+1,nflai)*bb(ix,iy,2)/bb(ix,iy,3)-vdrift,
     2          ua_max),ua_min)
            enddo
          enddo
* cell face x
          do iy=-1,ny
            do ix=0,nx
              if(rightix(ix-1,iy).ne.nx+1) then
               if(leftix(ix,iy).ne.-2) then
                  fnixb(ix,iy+1,nflai)=
     1             fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)
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
                fniyb(ix+1,iy,nflai)=fna(ix,iy,1,is)
              else
                fniyb(ix+1,iy,nflai)=0.0_R8
              endif
            enddo
          enddo
          fniyb(:,ny+1,nflai)=0.0_R8
        endif
      enddo

      if(eirene_extrapolate.eq.0) then
        ! do nothing
      elseif(eirene_extrapolate.eq.1) then
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
        write(*,*) eirene_extrapolate,
     1   ' not avaible as an eirene_extrap option'
      endif

      if(print_minmax.ne.0) then
         te_min_val=minval(teb(0:nx+1,0:ny+1))/ev
         te_max_val=maxval(teb(0:nx+1,0:ny+1))/ev
         ti_min_val=minval(tib(0:nx+1,0:ny+1))/ev
         ti_max_val=maxval(tib(0:nx+1,0:ny+1))/ev
         na_min_val=minval(dnib(0:nx+1,0:ny+1,1:nflai))
         na_max_val=maxval(dnib(0:nx+1,0:ny+1,1:nflai))
         ua_min_val=minval(uub(0:nx+1,0:ny+1,1:nflai))
         ua_max_val=maxval(uub(0:nx+1,0:ny+1,1:nflai))
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
         call my_out(70,nx,ny,uub(0:nx+1,0:ny+1,nflai),'copy_uub'//chns)
         call my_out(70,nx,ny,vvb(0:nx+1,0:ny+1,nflai),'copy_vvb'//chns)
         call my_out(70,nx,ny,wwb(0:nx+1,0:ny+1,nflai),'copy_wwb'//chns)
         call my_out(70,nx,ny,upb(0:nx+1,0:ny+1,nflai),'copy_upb'//chns)
        endif
       enddo
       call my_out(70,nx,ny,teb(0:nx+1,0:ny+1),'copy_teb')
       call my_out(70,nx,ny,tib(0:nx+1,0:ny+1),'copy_tib')
       call my_out(70,nx,ny,pob(0:nx+1,0:ny+1),'copy_pob')
      endif 


      ncall=ncall+1
      call subend ()
      return
      end subroutine copy_background
c
      subroutine write_f30(nx, ny)
      use b2mod_types
      use b2mod_indirect
      use b2mod_geo
#ifdef B25_EIRENE
#ifndef SOLPS4_3
      use eirmod_ctrig
#else
      use ctrig    ! IGNORE
#endif
      use eirmod_cinit, only: fort_lc
#else
      use b2mod_eirbra, only: fort_lc
#endif
      implicit none
      integer nx,ny
      integer ix,iy,i,nniso

      integer nxcut1(nncutmax),nxcut2(nncutmax),
     &        nycut1(nncutmax),nycut2(nncutmax),
     &        nxiso1(nncutmax),nxiso2(nncutmax),
     &        nyiso1(nncutmax),nyiso2(nncutmax)

      character*7 filename
      character*54 format_string(4)
      integer, save :: new_format
      logical ex

#ifdef B25_EIRENE
      integer idum
      integer, save :: exp_location, dot_location, spacing
      character*1 :: ldum(0:3)
      character*80 line
      data new_format/1/
      logical open_file, streql
      external open_file, streql
#else
      data new_format/2/
#endif
      data format_string/ '(1P,4E15.7)', '(1P,4E16.8)',
     .                       '(4E15.7)',    '(4E16.8)'/

      filename=fort_lc//'30'
      inquire(file=filename,exist=ex)
      if (ex) then
        write(*,*) 'Skipping writing '//trim(filename)//
     .   ' because file already present'
        return
      endif
      nniso=nncut-1
      if(nncut.eq.1) then
        if (nnreg(0).eq.1) then
          nxcut1=0
          nxcut2=0
          nycut1=0
          nycut2=0
        else
          nxcut1(1)=leftcut(1)
          nxcut1(2)=rightcut(1)
          nxcut2(1)=rightcut(1)+1
          nxcut2(2)=leftcut(1)+1
          nycut1(1)=bottomcut(1)+1
          nycut1(2)=bottomcut(1)+1
          nycut2(1)=topcut(1)
          nycut2(2)=topcut(1)
        endif
      elseif(nncut.eq.2) then
        nxcut1(1)=leftcut(1)
        nxcut1(2)=leftcut(2)
        nxcut1(3)=rightcut(2)
        nxcut1(4)=rightcut(1)
        nxcut2(1)=rightcut(1)+1
        nxcut2(2)=rightcut(2)+1
        nxcut2(3)=leftcut(2)+1
        nxcut2(4)=leftcut(1)+1
        nycut1(1)=bottomcut(1)+1
        nycut1(2)=bottomcut(2)+1
        nycut1(3)=bottomcut(2)+1
        nycut1(4)=bottomcut(1)+1
        nycut2(1)=topcut(1)
        nycut2(2)=topcut(2)
        nycut2(3)=topcut(2)
        nycut2(4)=topcut(1)
        ix=0
        do while (rightix(ix,ny/2).ne.nx+1 .and. ix.lt.nx)
          ix=ix+1
        enddo
        if(ix.ne.nx) nxiso1(1)=ix+1
        do while (leftix(ix,ny/2).ne.-2 .and. ix.lt.nx)
          ix=ix+1
        enddo
        if(ix.ne.nx) nxiso2(1)=ix
        nyiso1(1)=0
        nyiso2(1)=ny
      endif

#ifdef B25_EIRENE
cxpb Deduce from fort.33 triangle vertices files which format precision to use
      filename=fort_lc//'33'
      if (open_file(33,trim(filename),'')) then
        rewind(33)
        read(33,*,end=933,err=933) idum
        if (nrknot.eq.0) nrknot = idum
        if (nrknot.gt.0) then
          read(33,'(A80)') line
          exp_location = 81
          dot_location = 0
          if (index(line,'E').ne.0)
     .     exp_location = min(exp_location,index(line,'E'))
          if (index(line,'e').ne.0)
     .     exp_location = min(exp_location,index(line,'e'))
          if (index(line,'d').ne.0)
     .     exp_location = min(exp_location,index(line,'d'))
          if (index(line,'D').ne.0)
     .     exp_location = min(exp_location,index(line,'D'))
          if (exp_location.ne.81)
     .     dot_location = index(line(1:exp_location),'.')
          if (dot_location.ne.0) then
            if (exp_location-dot_location.eq.8) then
              new_format = 1
            else if (exp_location-dot_location.eq.9) then
              new_format = 2
            else
              call xerrab
     .         ('Unrecognized format (1) in '//fort_lc//'33 file')
            end if
          else
            call xerrab
     .       ('Unrecognized format (2) in '//fort_lc//'33 file')
          end if
cxpb Now deal with older files without the 1P in the format string
cxpb If the 1P was missing in the format, the numbers will start with 0
          spacing=index(line(dot_location+1:len(line)),'.')
          do i = 0, 3
            ldum(i) = line(dot_location-1+i*spacing:
     .                     dot_location-1+i*spacing)
            if (.not.streql(ldum(i),'0')) exit
            if (i.eq.3) new_format = new_format+2
          end do
        else
          new_format = 2
        end if
        goto 939
  933   new_format = 2
  939   continue
        close(33)
      else
        new_format = 2
      end if
#endif
      write(*,'(4a)')
     . 'Writing ',fort_lc,'30 with ', format_string(new_format)
      rewind(30)
      write (30,*) 'Geometry from B2'
      write (30,*) ' '
      if (nnreg(0).eq.1) then
        write (30,*) nx,ny,' 0 '
        write (30,*) ' 0 0 0 0 '
      else
        write (30,*) nx,ny,nncut*2
        write (30,*)
     1   (nxcut1(i),nxcut2(i),nycut1(i),nycut2(i),i=1,nncut*2)
      endif
      if(nniso.gt.0) then
        write (30,*) nniso
        write (30,*)
     1   (nxiso1(i),nxiso2(i),nyiso1(i),nyiso2(i),i=1,nniso)
      endif
      write (30,*) ' '
      do ix=0,nx-1
        do iy=0,ny-1
          write (30,format_string(new_format))
     1     crx(ix,iy,1),crx(ix,iy,3),crx(ix,iy,2),crx(ix,iy,0)
          write (30,format_string(new_format))
     1     cry(ix,iy,1),cry(ix,iy,3),cry(ix,iy,2),cry(ix,iy,0)
        end do
      end do
      close (30)

      return
 1028 format()
      end subroutine write_f30
c
      subroutine write_f31(nx, ny)
      use b2mod_types
      use b2mod_eirene_globals
      use b2mod_b2cmpa
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
      integer nx,ny
      integer nxx,nyy
      real*8, allocatable :: dummy(:,:,:)

csw 2013 corrected
c      nxx=nxdd
c      nyy=nydd
#ifdef B25_EIRENE
      nxx=ndx
      nyy=ndy
#else
      nxx=nx
      nyy=ny
#endif

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

      return
      end subroutine write_f31
