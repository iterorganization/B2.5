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
  910 FORMAT(5(E16.8))
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
      use eirmod_extrab25
      use eirmod_braeir
      implicit none
#     include <COUPLE/KOPPLDIM.F>
      integer nx,ny,ns
      real (kind=R8) :: BoRiS
      integer ix,iy,is,ncall,eirene_extrap,in,eirene_fixmap,
     1 eirene_fixuub
      external ipgetr, ipgeti
      real (kind=R8) ::
     * te_min, te_max, ti_min, ti_max, na_min, na_max, ua_min, ua_max
      save te_min, te_max, ti_min, ti_max,
     1 na_min, na_max, ua_min, ua_max, ncall,
     2 eirene_extrap, eirene_fixmap, eirene_fixuub
      data te_min/0.0_R8/, te_max/1.0e30_R8/,
     1 ti_min/0.0_R8/, ti_max/1.0e30_R8/,
     1 na_min/0.0_R8/, na_max/1.0e30_R8/,
     1 ua_min/-2.99792458e8_R8/, ua_max/+2.99792458e8_R8/,
     2 ncall/0/, eirene_extrap/1/, eirene_fixmap/1/, eirene_fixuub/0/   !srv 28.07.10
      real (kind=R8) :: kintmp
      real (kind=R8) :: hy1
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)

      if(ncall.eq.0) then
	call ipgetr('eirene_te_min',te_min)
	call ipgetr('eirene_te_max',te_max)
	call ipgetr('eirene_ti_min',ti_min)
	call ipgetr('eirene_ti_max',ti_max)
	call ipgetr('eirene_na_min',na_min)
	call ipgetr('eirene_na_max',na_max)
	call ipgetr('eirene_ua_min',ua_min)
	call ipgetr('eirene_ua_max',ua_max)
	call ipgeti('eirene_extrap',eirene_extrap)
	call ipgeti('eirene_fixmap',eirene_fixmap)
	call ipgeti('eirene_fixuub',eirene_fixuub)
      endif

      if(nx.ge.NDXP) call xerrab('eirene_mc: increase size of ndxp')
      if(ny.ge.NDYP) call xerrab('eirene_mc: increase size of ndyp')
* cell centred
      do iy=-1,ny
        do ix=-1,nx
          teb(ix+1,iy+1)=max(min(te(ix,iy),te_max*ev),te_min*ev)
          tib(ix+1,iy+1)=max(min(ti(ix,iy),ti_max*ev),ti_min*ev)
          pob(ix+1,iy+1)=po(ix,iy)
	    prb(ix+1,iy+1)=teb(ix+1,iy+1)*
     1	max(min(ne(ix,iy),na_max),na_min)
          rrb(ix+1,iy+1)=abs(bb(ix,iy,0)/bb(ix,iy,3))
          volb(ix+1,iy+1)=vol(ix,iy)
          bfeldb(ix+1,iy+1)=bb(ix,iy,3)
          qqc(ix+1,iy+1)=qc(ix,iy)
        enddo
      enddo
* cell face x
      do iy=-1,ny
	  do ix=0,nx
	    if(eirene_fixmap.eq.0) then
	    feixb(ix,iy+1)=fhi(ix,iy,0)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .          kinrgy(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .         (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
              feixb(ix,iy+1)=feixb(ix,iy+1)+kintmp*(1.0_R8-BoRiS)
     .         *fna(ix,iy,0,is)
            enddo
	    feexb(ix,iy+1)=fhe(ix,iy,0)
	  else
	    if(rightix(ix-1,iy).ne.nx+1) then
	      feixb(ix,iy+1)=fhi(rightix(ix-1,iy),rightiy(ix-1,iy),0)
              do is=0,ns-1
                kintmp = 
     .           (kinrgy(rightix(ix-1,iy),rightiy(ix-1,iy),is)*
     .            hx(ix-1,iy)+
     .            kinrgy(ix-1,iy,is)*
     .            hx(rightix(ix-1,iy),rightiy(ix-1,iy)))/
     .           (hx(ix-1,iy)+hx(rightix(ix-1,iy),rightiy(ix-1,iy)))
                feixb(ix,iy+1)=feixb(ix,iy+1)+kintmp*(1.0_R8-BoRiS)
     .           *fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)
              enddo
	      feexb(ix,iy+1)=fhe(rightix(ix-1,iy),rightiy(ix-1,iy),0)
	    else
	      feixb(ix,iy+1)=0.0_R8
	      feexb(ix,iy+1)=0.0_R8
	    endif
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
            kintmp = 
     .       (kinrgy(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     .        kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
     .       (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            feiyb(ix+1,iy)=feiyb(ix+1,iy)+kintmp*(1.0_R8-BoRiS)
     .       *fna(ix,iy,1,is)
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
          if(nflai.gt.NFL) stop 'eirene_mc: increase size of nfl'
* cell centred
          do iy=-1,ny
            do ix=-1,nx
	       dnib(ix+1,iy+1,nflai)=max(min(na(ix,iy,is),na_max),na_min)
              upb(ix+1,iy+1,nflai)=max(min(ua(ix,iy,is),ua_max),ua_min)
	       uudiab(ix+1,iy+1,nflai)=max(min(
     1         wadia(ix,iy,0,is)+vaecrb(ix,iy,0,is),ua_max),ua_min)
	       vvdiab(ix+1,iy+1,nflai)=max(min(
     1	       wadia(ix,iy,1,is)+vaecrb(ix,iy,1,is),ua_max),ua_min)
	       prb(ix+1,iy+1)=prb(ix+1,iy+1)+tib(ix+1,iy+1)*
     1	       dnib(ix+1,iy+1,nflai)
	       if(eirene_fixuub.eq.1) then
		    if(gs(ix,iy,0).ne.0_R8.and.rightix(ix,iy).ne.nx+1) then
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   (fna(ix,iy,0,is)/gs(ix,iy,0)+
     2		   fna(rightix(ix,iy),rightiy(ix,iy),0,is)/
     3		   gs(rightix(ix,iy),rightiy(ix,iy),0))/2.0_R8/
     4		   na(ix,iy,is)/qz(ix,iy,1),
     5		   ua_max),ua_min)
		    elseif(
     1		 gs(ix,iy,0).eq.0_R8.and.rightix(ix,iy).ne.nx+1) then
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   fna(rightix(ix,iy),rightiy(ix,iy),0,is)/
     3		   gs(rightix(ix,iy),rightiy(ix,iy),0)/
     4		   na(ix,iy,is)/qc(rightix(ix,iy),rightiy(ix,iy)),
     5		   ua_max),ua_min)
		    else
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   fna(ix,iy,0,is)/gs(ix,iy,0)/
     4		   na(ix,iy,is)/qc(ix,iy),
     5		   ua_max),ua_min)
		    endif
             elseif(eirene_fixuub.eq.2) then
                uub(ix+1,iy+1,nflai)=max(min(
     1           upb(ix+1,iy+1,nflai)*rrb(ix+1,iy+1),
     2           ua_max),ua_min)
	       endif
            enddo
          enddo
* cell face x
          do iy=-1,ny
            do ix=0,nx
	       if(eirene_fixmap.eq.0) then
		    if(leftix(ix,iy).ne.-2) then
		     fnixb(ix,iy+1,nflai)=fna(ix,iy,0,is)
		     if(eirene_fixuub.eq.0) then
		      uub(ix,iy+1,nflai)=max(min(
     1		    fna(ix,iy,0,is)/gs(ix,iy,0)/qc(ix,iy)*2.0_R8/
     1		    (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)),
     3		    ua_max),ua_min)
		     endif
		    else
		     fnixb(ix,iy+1,nflai)=0.0_R8
		     if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		    endif
	       else
		    if(rightix(ix-1,iy).ne.nx+1) then
		     if(leftix(ix,iy).ne.-2) then
		      fnixb(ix,iy+1,nflai)=
     1		     fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)
		      if(eirene_fixuub.eq.0) then
		       uub(ix,iy+1,nflai)=max(min(
     1		       fna(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)/
     2		       gs(rightix(ix-1,iy),rightiy(ix-1,iy),0)*2.0_R8/
     1		       (na(rightix(ix-1,iy),rightiy(ix-1,iy),is)+
     2		       na(ix-1,iy,is))/
     3		       qc(rightix(ix-1,iy),rightiy(ix-1,iy)),
     3		       ua_max),ua_min)
		      endif
		     else
		      write(*,*) 'WARNING: ',
     1		     'You should not have executed this block!'
		      fnixb(ix,iy+1,nflai)=0.0_R8
		      if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		     endif
		    else
		     fnixb(ix,iy+1,nflai)=0.0_R8
		     if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		    endif
	       endif  !eirene_fixmap
	      enddo
	    enddo
          fnixb(nx+1,:,nflai)=0.0_R8
          if(eirene_fixuub.eq.0) uub(nx+1,:,nflai)=0.0_R8
* cell face y
          do iy=0,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).ne.-2) then
                fniyb(ix+1,iy,nflai)=fna(ix,iy,1,is)
                vvb(ix+1,iy,nflai)=max(min(
     1           fna(ix,iy,1,is)/gs(ix,iy,1)*2.0_R8/
     1           (na(ix,iy,is)+na(bottomix(ix,iy),bottomiy(ix,iy),is)),
     3           ua_max),ua_min)
              else
                fniyb(ix+1,iy,nflai)=0.0_R8
                vvb(ix+1,iy,nflai)=0.0_R8
              endif
            enddo
          enddo
          fniyb(:,ny+1,nflai)=0.0_R8
          vvb(:,ny+1,nflai)=0.0_R8
        endif
      enddo

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
     1		 dnib(rightix(ix,iy)+1,rightiy(ix,iy)+1,:)
	      endif
	      if(rightix(ix,iy).eq.nx+1) then
		teb(ix+1,iy+1)=teb(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(leftix(ix,iy)+1,leftiy(ix,iy)+1,:)
	      endif
	      if(bottomiy(ix,iy).eq.-2) then
		teb(ix+1,iy+1)=teb(topix(ix,iy)+1,topiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(topix(ix,iy)+1,topiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(topix(ix,iy)+1,topiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(topix(ix,iy)+1,topiy(ix,iy)+1,:)
	      endif
	      if(topiy(ix,iy).eq.ny+1) then
		teb(ix+1,iy+1)=teb(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1,:)
	      endif
	    enddo
	  enddo
	enddo
      else
	write(*,*) eirene_extrap,
     1	 ' not avaible as an eirene_extrap option'
      endif
      ncall=ncall+1
      return
      end subroutine copy_background
c
      subroutine copy_background_(nx, ny, ns, BoRiS)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_geo
      use b2mod_constants
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use eirmod_extrab25
      use eirmod_braeir
      implicit none
#     include <COUPLE/KOPPLDIM.F>
      integer nx,ny,ns
      real (kind=R8) :: BoRiS
      character chns*2                                                 !srv 16.09.02
      integer, save :: drift_style = 0, iout = 0                        !srv 21.06.08
      integer ix,iy,is,ncall,eirene_extrap,in,eirene_fixmap,
     1 eirene_fixuub
      external ipgetr, ipgeti
      real (kind=R8) ::
     * te_min, te_max, ti_min, ti_max, na_min, na_max, ua_min, ua_max
      save te_min, te_max, ti_min, ti_max,
     1 na_min, na_max, ua_min, ua_max, ncall,
     2 eirene_extrap, eirene_fixmap, eirene_fixuub
      data te_min/0.0_R8/, te_max/1.0e30_R8/,
     1 ti_min/0.0_R8/, ti_max/1.0e30_R8/,
     1 na_min/0.0_R8/, na_max/1.0e30_R8/,
     1 ua_min/-2.99792458e8_R8/, ua_max/+2.99792458e8_R8/,
     2 ncall/0/, eirene_extrap/1/, eirene_fixmap/1/, eirene_fixuub/0/   !srv 28.07.10
      real (kind=R8) :: kintmp
      real (kind=R8) :: hy1
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)

      if(ncall.eq.0) then
	call ipgetr('eirene_te_min',te_min)
	call ipgetr('eirene_te_max',te_max)
	call ipgetr('eirene_ti_min',ti_min)
	call ipgetr('eirene_ti_max',ti_max)
	call ipgetr('eirene_na_min',na_min)
	call ipgetr('eirene_na_max',na_max)
	call ipgetr('eirene_ua_min',ua_min)
	call ipgetr('eirene_ua_max',ua_max)
	call ipgeti('eirene_extrap',eirene_extrap)
	call ipgeti('eirene_fixmap',eirene_fixmap)
	call ipgeti('eirene_fixuub',eirene_fixuub)
      call ipgeti('b2tfnb_drift_style', drift_style)                 !srv 21.06.08
      call ipgeti('copy_background_iout', iout)                      !srv 14.07.10
      endif
	if(iout.ne.0) then                                             !srv 18.06.08 {
	 do is = 0,ns-1
	 write (chns,'(i2.2)') is
       call my_out(70,nx,ny,fna_eir(-1,-1,0,is),'copy_fna_eirx'//chns)
       call my_out(70,nx,ny,fna_eir(-1,-1,1,is),'copy_fna_eiry'//chns)
	 call my_out(70,nx,ny,vaecrb(-1,-1,0,is),'copy_vecrbx'//chns)
	 call my_out(70,nx,ny,vaecrb(-1,-1,1,is),'copy_vecrby'//chns)
	 call my_out(70,nx,ny,wadia(-1,-1,0,is),'copy_wadiax'//chns)
	 call my_out(70,nx,ny,wadia(-1,-1,1,is),'copy_wadiay'//chns)
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
	endif                                                          !srv 18.06.08 }

      if(nx.ge.NDXP) call xerrab('eirene_mc: increase size of ndxp')
      if(ny.ge.NDYP) call xerrab('eirene_mc: increase size of ndyp')

* cell centred
      do iy=-1,ny
        do ix=-1,nx
          teb(ix+1,iy+1)=max(min(te(ix,iy),te_max*ev),te_min*ev)
          tib(ix+1,iy+1)=max(min(ti(ix,iy),ti_max*ev),ti_min*ev)
          pob(ix+1,iy+1)=po(ix,iy)
	  prb(ix+1,iy+1)=teb(ix+1,iy+1)*
     1	   max(min(ne(ix,iy),na_max),na_min)
          rrb(ix+1,iy+1)=abs(bb(ix,iy,0)/bb(ix,iy,3))
          volb(ix+1,iy+1)=vol(ix,iy)
          bfeldb(ix+1,iy+1)=bb(ix,iy,3)
        enddo
      enddo

* cell face x
      do iy=-1,ny
	do ix=0,nx
	  if(eirene_fixmap.eq.0) then
	    feixb(ix,iy+1)=fhi_eir(ix,iy,0)                               !srv 14.07.10
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,iy,is)*hx(leftix(ix,iy),leftiy(ix,iy))+
     .          kinrgy(leftix(ix,iy),leftiy(ix,iy),is)*hx(ix,iy))/
     .         (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
              feixb(ix,iy+1)=feixb(ix,iy+1)+kintmp*(1.0_R8-BoRiS)
     .         *fna_eir(ix,iy,0,is)                                     !srv 14.07.10
            enddo
	    feexb(ix,iy+1)=fhe_eir(ix,iy,0)                               !srv 14.07.10
	  else
	    if(rightix(ix-1,iy).ne.nx+1) then
	      feixb(ix,iy+1)=fhi_eir(rightix(ix-1,iy),rightiy(ix-1,iy),0) !srv 14.07.10
              do is=0,ns-1
                kintmp = 
     .           (kinrgy(rightix(ix-1,iy),rightiy(ix-1,iy),is)*
     .            hx(ix-1,iy)+
     .            kinrgy(ix-1,iy,is)*
     .            hx(rightix(ix-1,iy),rightiy(ix-1,iy)))/
     .           (hx(ix-1,iy)+hx(rightix(ix-1,iy),rightiy(ix-1,iy)))
                feixb(ix,iy+1)=feixb(ix,iy+1)+kintmp*(1.0_R8-BoRiS)
     .           *fna_eir(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)       !srv 14.07.10
              enddo
	      feexb(ix,iy+1)=fhe_eir(rightix(ix-1,iy),rightiy(ix-1,iy),0) !srv 14.07.10
	    else
	      feixb(ix,iy+1)=0.0_R8
	      feexb(ix,iy+1)=0.0_R8
	    endif
	  endif
	enddo
      enddo
      feixb(nx+1,:)=0.0_R8
      feexb(nx+1,:)=0.0_R8

* cell face y
      do iy=0,ny
        do ix=-1,nx
          feiyb(ix+1,iy)=fhi_eir(ix,iy,1)                               !srv 14.07.10
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ix,iy,is)*hy1(bottomix(ix,iy),bottomiy(ix,iy))+
     .        kinrgy(bottomix(ix,iy),bottomiy(ix,iy),is)*hy1(ix,iy))/
     .       (hy1(bottomix(ix,iy),bottomiy(ix,iy))+hy1(ix,iy))
            feiyb(ix+1,iy)=feiyb(ix+1,iy)+kintmp*(1.0_R8-BoRiS)
     .       *fna_eir(ix,iy,1,is)                                       !srv 14.07.10
          enddo
          feeyb(ix+1,iy)=fhe_eir(ix,iy,1)                               !srv 14.07.10
        enddo
      enddo
      feiyb(:,ny+1)=0.0_R8
      feeyb(:,ny+1)=0.0_R8
* by ns
      nflai=0
      do is=0,ns-1
        if(.not.is_neutral(is)) then
          nflai=nflai+1
          if(nflai.gt.NFL) stop 'eirene_mc: increase size of nfl'
* cell centred
          do iy=-1,ny
            do ix=-1,nx
	       dnib(ix+1,iy+1,nflai)=max(min(na(ix,iy,is),na_max),na_min)
              upb(ix+1,iy+1,nflai)=max(min(ua(ix,iy,is),ua_max),ua_min)
	       if(drift_style.eq.0) then                                  !srv 21.06.08 14.07.10 {
	        uudiab(ix+1,iy+1,nflai)=max(min(
     1         wadia(ix,iy,0,is)+vaecrb(ix,iy,0,is),ua_max),ua_min)
	        vvdiab(ix+1,iy+1,nflai)=max(min(
     1	   wadia(ix,iy,1,is)+vaecrb(ix,iy,1,is),ua_max),ua_min)
	       else
              if(leftix(ix,iy).eq.-2) then
	         uudiab(ix+1,iy+1,nflai)=max(min(
     1          wadia(rightix(ix,iy),rightiy(ix,iy),0,is)+
     &          vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is),
     &          ua_max),ua_min)
              elseif(rightix(ix,iy).eq.nx+1) then
	         uudiab(ix+1,iy+1,nflai)=max(min(
     1          wadia(ix,iy,0,is)+vaecrb(ix,iy,0,is),ua_max),ua_min)
              else
	         uudiab(ix+1,iy+1,nflai)=max(min(0.5_R8*
     &          (wadia(rightix(ix,iy),rightiy(ix,iy),0,is)+
     &           vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is)+
     1           wadia(ix,iy,0,is)+vaecrb(ix,iy,0,is)),ua_max),ua_min)
              endif      
              if(bottomiy(ix,iy).eq.-2) then
	         vvdiab(ix+1,iy+1,nflai)=max(min(
     1	    wadia(topix(ix,iy),topiy(ix,iy),1,is)+
     &          vaecrb(topix(ix,iy),topiy(ix,iy),1,is),ua_max),ua_min)
              elseif(topiy(ix,iy).eq.ny+1) then
	         vvdiab(ix+1,iy+1,nflai)=max(min(
     1	    wadia(ix,iy,1,is)+vaecrb(ix,iy,1,is),ua_max),ua_min)
              else
	         vvdiab(ix+1,iy+1,nflai)=max(min(0.5_R8*
     &          (wadia(topix(ix,iy),topiy(ix,iy),1,is)+
     &           vaecrb(topix(ix,iy),topiy(ix,iy),1,is)+
     1	     wadia(ix,iy,1,is)+vaecrb(ix,iy,1,is)),ua_max),ua_min)
              endif      
	       endif                                                      !srv 21.06.08 14.07.10 }
	       prb(ix+1,iy+1)=prb(ix+1,iy+1)+tib(ix+1,iy+1)*
     1	       dnib(ix+1,iy+1,nflai)
	       if(eirene_fixuub.eq.1) then
		    if(gs(ix,iy,0).ne.0_R8.and.rightix(ix,iy).ne.nx+1) then
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   (fna_eir(ix,iy,0,is)/gs(ix,iy,0)+
     2		   fna_eir(rightix(ix,iy),rightiy(ix,iy),0,is)/       !srv 14.07.10
     3		   gs(rightix(ix,iy),rightiy(ix,iy),0))/2.0_R8/
     4		   na(ix,iy,is)/qz(ix,iy,1),
     5		   ua_max),ua_min)
		    elseif(
     1		 gs(ix,iy,0).eq.0_R8.and.rightix(ix,iy).ne.nx+1) then
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   fna_eir(rightix(ix,iy),rightiy(ix,iy),0,is)/       !srv 14.07.10
     3		   gs(rightix(ix,iy),rightiy(ix,iy),0)/
     4		   na(ix,iy,is)/qc(rightix(ix,iy),rightiy(ix,iy)),
     5		   ua_max),ua_min)
		    else
		     uub(ix+1,iy+1,nflai)=max(min(
     1		   fna_eir(ix,iy,0,is)/gs(ix,iy,0)/                   !srv 14.07.10
     4		   na(ix,iy,is)/qc(ix,iy),
     5		   ua_max),ua_min)
		    endif
             elseif(eirene_fixuub.eq.2) then
                uub(ix+1,iy+1,nflai)=max(min(
     1           upb(ix+1,iy+1,nflai)*rrb(ix+1,iy+1),
     2           ua_max),ua_min)
	       endif
            enddo
          enddo

* cell face x
          do iy=-1,ny
            do ix=0,nx
	       if(eirene_fixmap.eq.0) then
		    if(leftix(ix,iy).ne.-2) then
		     fnixb(ix,iy+1,nflai)=fna_eir(ix,iy,0,is)               !srv 14.07.10
		     if(eirene_fixuub.eq.0) then
		      uub(ix,iy+1,nflai)=max(min(
     1		    fna_eir(ix,iy,0,is)/gs(ix,iy,0)/qc(ix,iy)*2.0_R8/ !srv 14.07.10
     1		    (na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is)),
     3		    ua_max),ua_min)
		     endif
		    else
		     fnixb(ix,iy+1,nflai)=0.0_R8
		     if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		    endif
	       else
		    if(rightix(ix-1,iy).ne.nx+1) then
		     if(leftix(ix,iy).ne.-2) then
		      fnixb(ix,iy+1,nflai)=
     1		     fna_eir(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)  !srv 14.07.10
		      if(eirene_fixuub.eq.0) then
		       uub(ix,iy+1,nflai)=max(min(
     1		    fna_eir(rightix(ix-1,iy),rightiy(ix-1,iy),0,is)/  !srv 14.07.10
     2		    gs(rightix(ix-1,iy),rightiy(ix-1,iy),0)*2.0_R8/
     1		    (na(rightix(ix-1,iy),rightiy(ix-1,iy),is)+
     2		    na(ix-1,iy,is))/
     3		    qc(rightix(ix-1,iy),rightiy(ix-1,iy)),
     3		    ua_max),ua_min)
		      endif
		     else
		      write(*,*) 'WARNING: ',
     1		     'You should not have executed this block!'
		      fnixb(ix,iy+1,nflai)=0.0_R8
		      if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		     endif
		    else
		     fnixb(ix,iy+1,nflai)=0.0_R8
		     if(eirene_fixuub.eq.0) uub(ix,iy+1,nflai)=0.0_R8
		    endif
	       endif  !eirene_fixmap
	      enddo
	    enddo
          fnixb(nx+1,:,nflai)=0.0_R8
          if(eirene_fixuub.eq.0) uub(nx+1,:,nflai)=0.0_R8

* cell face y
          do iy=0,ny
            do ix=-1,nx
              if(bottomiy(ix,iy).ne.-2) then
                fniyb(ix+1,iy,nflai)=fna_eir(ix,iy,1,is)               !srv 14.07.10
                vvb(ix+1,iy,nflai)=max(min(
     1           fna_eir(ix,iy,1,is)/gs(ix,iy,1)*2.0_R8/               !srv 14.07.10
     1           (na(ix,iy,is)+na(bottomix(ix,iy),bottomiy(ix,iy),is)),
     3           ua_max),ua_min)
              else
                fniyb(ix+1,iy,nflai)=0.0_R8
                vvb(ix+1,iy,nflai)=0.0_R8
              endif
            enddo
          enddo
          fniyb(:,ny+1,nflai)=0.0_R8
          vvb(:,ny+1,nflai)=0.0_R8
        endif
      enddo

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
     1		 dnib(rightix(ix,iy)+1,rightiy(ix,iy)+1,:)
	      endif
	      if(rightix(ix,iy).eq.nx+1) then
		teb(ix+1,iy+1)=teb(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(leftix(ix,iy)+1,leftiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(leftix(ix,iy)+1,leftiy(ix,iy)+1,:)
	      endif
	      if(bottomiy(ix,iy).eq.-2) then
		teb(ix+1,iy+1)=teb(topix(ix,iy)+1,topiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(topix(ix,iy)+1,topiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(topix(ix,iy)+1,topiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(topix(ix,iy)+1,topiy(ix,iy)+1,:)
	      endif
	      if(topiy(ix,iy).eq.ny+1) then
		teb(ix+1,iy+1)=teb(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		tib(ix+1,iy+1)=tib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		pob(ix+1,iy+1)=pob(bottomix(ix,iy)+1,bottomiy(ix,iy)+1)
		dnib(ix+1,iy+1,:)=
     1		 dnib(bottomix(ix,iy)+1,bottomiy(ix,iy)+1,:)
	      endif
	    enddo
	  enddo
	enddo
      else
	write(*,*) eirene_extrap,
     1	 ' not avaible as an eirene_extrap option'
      endif
      ncall=ncall+1
      return
      end subroutine copy_background_
c
      subroutine write_f30(nx, ny)
      use b2mod_types
      use b2mod_indirect
      use b2mod_geo
      implicit none
      integer nx,ny
      integer ix,iy,i,nniso

      integer nxcut1(nncutmax),nxcut2(nncutmax),
     &        nycut1(nncutmax),nycut2(nncutmax),
     &        nxiso1(nncutmax),nxiso2(nncutmax),
     &        nyiso1(nncutmax),nyiso2(nncutmax)

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
          write (30,'(4e15.7)')
     1     crx(ix,iy,1),crx(ix,iy,3),crx(ix,iy,2),crx(ix,iy,0)
          write (30,'(4e15.7)')
     1     cry(ix,iy,1),cry(ix,iy,3),cry(ix,iy,2),cry(ix,iy,0)
        end do
      end do
      close (30)

      return
      end subroutine write_f30
c
      subroutine write_f31(nx, ny)
      use b2mod_types
      use b2mod_eirene_globals
      use b2mod_b2cmpa
      use eirmod_extrab25
      use eirmod_braeir
      implicit none
#     include <COUPLE/KOPPLDIM.F>
      integer nx,ny
      integer nxx,nyy
      real*8, allocatable :: dummy(:,:,:)

csw 2013 corrected
c      nxx=nxdd
c      nyy=nydd
      nxx=nx+2
      nyy=ny

      write(*,*) 'NX, NY, NFLAI = ',nx,ny,nflai
      rewind(31)
c ion density
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,DNIB(0,0,1))
c poloidal velocity (derived from fna(,,0,))
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,UUB(0,0,1))
c radial velocity (derived from fna(,,1))
      CALL GFSUB3(31,NX,NY,NXx,NYy,NFLAI,VVB(0,0,1))
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
c cell volumes
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,VOLB(0,0))
c magnitude of the magnetic field
      CALL GFSUB3(31,NX,NY,NXx,NYy,1,BFELDB(0,0))
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
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      call gfsub3(31,nx,ny,nxx,nyy,1,dummy)
      deallocate(dummy)

c electric potential
c      CALL GFSUB3(31,NX,NY,NXx,NYy,1,POB(0,0))
      write(*,*) 'Background plasma written'
      close(31)

      return
      end subroutine write_f31
