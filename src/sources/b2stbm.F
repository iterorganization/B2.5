*-----------------------------------------------------------------------
*.specification

      subroutine b2stbm (nCv, nFc, nVx, ns, ismain, dtim, switch,
     &  geo, mpg, pl, dv, st_ext, sr, srw, main_call)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
      use b2mod_impgyro
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
#ifdef WG_RM
      use b2mod_sources
     1 , only : b2stbm_sna,b2stbm_smo,b2stbm_she,b2stbm_shi,b2stbm_sch,
     2          b2stbm_sne
#endif
      use b2mod_balance !djm Jan2017
     1 , only : balance_netcdf
     2        , b2stbm_sna0to1, b2stbm_smo0to3, b2stbm_she0to3,
     3          b2stbm_shi0to3
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      logical main_call
      real (kind=R8) :: dtim
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2StateExt), intent (inout) :: st_ext
*   ..input/output arguments 
      type (B2Derivatives), intent (inout) :: dv
      type (B2Source), intent (inout) :: sr
      type (B2SourceWork), intent (inout) :: srw
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBM provides linearised expressions for the sources of
*     density, parallel momentum, and heat, due to miscellaneous
*     application specific processes.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ia, is, is0, iCv, iFc, ireg
      real (kind=R8) :: t0, vthe, ua_tmp, dua
      real (kind=R8) :: n0_tmp, nr_tmp, nt_tmp, p0_tmp, pr_tmp, r0_tmp,
     *                  pe_tmp, re_tmp
      real (kind=R8) ::
     *  na_tot(nCv), na0(nCv,0:ns-1),
     *  ua0(nCv,0:ns-1), fna0(nFc,0:1,0:ns-1),
     *  fhi0(nFc,0:1), ti0(nCv),
     *  fhe0(nFc,0:1), fch0(nFc,0:1)
      integer, save :: ncall = 0
*   ..procedures
      real (kind=R8) ::
     * smin, smax
      external xertst, sfill,
     * smin, smax, b2xvsg, b2xvff, ipgeti, ipgetr
      intrinsic abs, mod, sqrt

*-----------------------------------------------------------------------
*.documentation-internal
*
*     This version of b2stbm imports the sources computed from an
*     external code such as IMPGYRO.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbm')
*   ..test nCv, nFc, ns, ismain
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     & 'faulty argument ismain')
*   ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nFc*2, geo%fcVol, 1, 'vol', '.gt.')
       call b2xvsg (nCv, geo%cvHx, 1, 'hx', '.gt.')
       call b2xvsg (nFc, geo%fcHt, 1, 'hy', '.gt.')
*    ..test range of cvQgam, fcQgam
       call xertst (-1.0_R8.lt.smin(nCv,geo%cvQgam(1,1),1).and.
     &   smax(nCv,geo%cvQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range cvQgam(,1)')
       call xertst (0.0_R8.lt.smin(nCv,geo%cvQgam(1,0),1).and.
     &   smax(nCv,geo%cvQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range cvQgam(,0)')
       call xertst (-1.0_R8.lt.smin(nFc,geo%fcQgam(1,1),1).and.
     &   smax(nFc,geo%fcQgam(1,1),1).lt.1.0_R8, 
     &   'faulty argument range fcQgam(,1)')
       call xertst (0.0_R8.lt.smin(nFc,geo%fcQgam(1,0),1).and.
     &   smax(nFc,geo%fcQgam(1,0),1).le.1.0_R8, 
     &   'faulty argument range fcQgam(,0)')
*    ..test state
       call b2xvps (nCv, nFc, ns, pl, dv)

*   ..test edge values of fluxes
#ifdef WG_TODO
       do is = 0, ns-1
        call b2xvff (nx, ny, fna(-1,-1,0,0,is), 'fna')
       enddo
       call b2xvff (nx, ny, fhe, 'fhe')
       call b2xvff (nx, ny, fhi, 'fhi')
       call b2xvff (nx, ny, fch, 'fch')
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
#endif
      endif

*   ..initialise sources to 0
      srw%sna0 = 0.0_R8
      srw%smo0 = 0.0_R8
      srw%sne0 = 0.0_R8
      srw%she0 = 0.0_R8
      srw%shi0 = 0.0_R8
      srw%sch0 = 0.0_R8
      srw%b2stbm_sna = 0.0_R8
      srw%b2stbm_smo = 0.0_R8
      srw%b2stbm_sne = 0.0_R8
      srw%b2stbm_she = 0.0_R8
      srw%b2stbm_shi = 0.0_R8
      srw%b2stbm_sch = 0.0_R8

* ..obtain sources from IMPGYRO
      if (switch%b2stbm_impgyro_mod.eq.0) goto 999

cxpb First provide plasma background (include neutral quantities)
cxpb Must convert fluxes to flux densities
      ti0 = 0.0_R8
      fhi0 = 0.0_R8
      na_tot = 0.0_R8
      do is = 0, ns-1
        if (is_neutral(is)) then
          is0 = b2espcr(is)
        else
          is0 = 0
        endif
        if (is0.ne.0) then
          do iCv = 1, nCv
            n0_tmp = 0.0_R8
            nr_tmp = 0.0_R8
            nt_tmp = 0.0_R8
            p0_tmp = 0.0_R8
            pe_tmp = 0.0_R8
            pr_tmp = 0.0_R8
            r0_tmp = 0.0_R8
            re_tmp = 0.0_R8
            do ia = 1, nnatmi
              if (latmscl(ia).eq.is0) then
                n0_tmp = n0_tmp + dab2(iCv,ia,1)
!WG_TODO        nr_tmp = nr_tmp + dab2(ixr+1,iyr+1,ia,1)
                nt_tmp = nt_tmp + dab2(iCv,ia,1)*
     .                            tab2(iCv,ia,1)
                p0_tmp = p0_tmp + pfluxa(iCv,ia,1)
                pe_tmp = pe_tmp + pefluxa(iCv,ia,1)
!WG_TODO        pr_tmp = pr_tmp + pfluxa(ixr+1,iyr+1,ia,1)
                r0_tmp = r0_tmp + rfluxa(iCv,ia,1)
                re_tmp = re_tmp + refluxa(iCv,ia,1)
              end if
            end do
            if (n0_tmp.gt.0.0_R8 .and.
     .          nr_tmp.gt.0.0_R8) then
              ua_tmp = (p0_tmp/n0_tmp+pr_tmp/nr_tmp)/2.0_R8
            elseif (n0_tmp.gt.0.0_R8) then
              ua_tmp = p0_tmp/n0_tmp
            elseif (nr_tmp.gt.0.0_R8) then
              ua_tmp = pr_tmp/nr_tmp
            else
              ua_tmp = 0.0_R8
            endif  
            na_tot(iCv) = na_tot(iCv) + n0_tmp
            na0(iCv,is) = pl%na(iCv,is) + n0_tmp
            ti0(iCv) = ti0(iCv) + nt_tmp
            ua0(iCv,is) = (pl%ua(iCv,is)*pl%na(icV,is)+ua_tmp*n0_tmp)
     /                    /na0(iCv,is)
          end do
          do iFc = 1, nFc
            fna0(iFc,0,is) = 
     .         dv%fna(iFc,0,is)/geo%fcS(iFc)/geo%fcQgam(iFc,0) + p0_tmp
            fna0(iFc,1,is) = dv%fna(iFc,1,is)/geo%fcS(iFc) +r0_tmp
            fhi0(iFc,0) = fhi0(iFc,0) + pe_tmp
            fhi0(iFc,1) = fhi0(iFc,1) + re_tmp
          end do
        else
          na0(:,is) = pl%na(:,is)
          ua0(:,is) = pl%ua(:,is)
          fna0(:,0,is) = dv%fna(:,0,is)/geo%fcS(:)/geo%fcQgam(iFc,0)
          fna0(:,1,is) = dv%fna(:,1,is)/geo%fcS(:)
        end if
      end do
      do is = 1, nnmoli
        do iCv = 1, nCv
          na_tot(iCv) = na_tot(iCv) + dmb2(iCv,is,1)
          ti0(iCv) = ti0(iCv) + 
     .       dmb2(iCv,is,1)*tmb2(iCv,is,1)
        end do
      end do
      do is = 1, nnioni
        do iCv = 1, nCv
          na_tot(iCv) = na_tot(iCv) + dib2(iCv,is,1)
          ti0(iCv) = ti0(iCv) + 
     .       dib2(iCv,is,1)*tib2(iCv,is,1)
        end do
      end do
      do iCv = 1, nCv
        ti0(iCv) = (ti0(iCv) + 
     .                  pl%ti(iCv)*dv%ni(iCv,0)) /
     .                 (na_tot(iCv) + dv%ni(iCv,0))
      end do
      do iFc = 1, nFc
        fhi0(iFc,0) = fhi0(iFc,0) + 
     .     dv%fhi(iFc,0)/geo%fcS(iFc)/geo%fcQgam(iFc,0)
        fhi0(iFc,1) = fhi0(iFc,1) + 
     .     dv%fhi(iFc,1)/geo%fcS(iFc)
        fhe0(iFc,0) = dv%fhe(iFc,0)/geo%fcS(iFc)/geo%fcQgam(iFc,0)
        fhe0(iFc,1) = dv%fhe(iFc,1)/geo%fcS(iFc)
        fch0(iFc,0) = dv%fch(iFc,0)/geo%fcS(iFc)/geo%fcQgam(iFc,0)
        fch0(iFc,1) = dv%fch(iFc,1)/geo%fcS(iFc)
      end do

      if (ncall.eq.0 .or.
     & (main_call .and. mod(ncall,switch%b2stbm_impgyro_mod).eq.0)) then
cxpb The call below obtains the following quantities from IMPGYRO
cxpb Main W data: na_ext, ua_ext, ta_ext
cxpb Sources due to W: sne_ext, smo_ext, she_ext, shi_ext
cxpb Derived quantities: ne_ext, ne2_ext, ni_ext, fa_ext, fhi_ext, ue_ext
c        call read_b2mod_external(nx,ny,ns,ismain)
cxpb In the current arrangement, IMPGYRO provides: na_ext, ua_ext, ta_ext
cxpb It would be best to have instead as primary quantities:
cxpb na_ext, fa_ext, ta_ext, fhi_ext
cxpb And then derive: ua_ext

       st_ext%sna = 0.0_R8
       st_ext%sne = 0.0_R8
       st_ext%smo = 0.0_R8
       st_ext%she = 0.0_R8
       st_ext%shi = 0.0_R8
       st_ext%sch = 0.0_R8
       st_ext%na = 0.0_R8
       st_ext%ua = 0.0_R8
       st_ext%ta = 0.0_R8

#ifdef WG_TODO
       call impgyro_main(
     *  nx, ny, ns, ns_ext,
     *  bb, na0, ne, ua0, te, ti0, po,
     *  fna0, fhe0, fhi0, fch0,
     *  sna_ext, sne_ext, smo_ext,
     *  she_ext, shi_ext, sch_ext,
     *  na_ext, ua_ext, ta_ext)

       call extend(na_ext,nx,ny,ns_ext,nx,ny,ns_ext)
       call extend(ua_ext,nx,ny,ns_ext,nx,ny,ns_ext)
       call extend(ta_ext,nx,ny,ns_ext,nx,ny,ns_ext)
        
       ne_ext = 0.0_R8
       ni_ext = 0.0_R8
       ue_ext = 0.0_R8
       ne2_ext = 0.0_R8
       do is = 0, ns_ext-1
         ne_ext(-1:nx,-1:ny) = ne_ext(-1:nx,-1:ny) +
     .     za_ext(-1:nx,-1:ny,is) * na_ext(-1:nx,-1:ny,is)
         ne2_ext(-1:nx,-1:ny) = ne2_ext(-1:nx,-1:ny) +
     .     za2_ext(-1:nx,-1:ny,is) * na_ext(-1:nx,-1:ny,is)
         ni_ext(-1:nx,-1:ny,0) = ni_ext(-1:nx,-1:ny,0) +
     .     na_ext(-1:nx,-1:ny,is)
         if (.not.is_ext_neutral(is)) ni_ext(-1:nx,-1:ny,1) =
     .     ni_ext(-1:nx,-1:ny,1) + na_ext(-1:nx,-1:ny,is)
cxpb Computing ua from fa is what we really want to do
         if (maxval(fa_ext).ne.0.0_R8.or.minval(fa_ext).ne.0.0_R8) then
          do iy=-1,ny
            do ix=-1,nx
              if (leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                ua_ext(ix,iy,is) = 0.5_R8*
     .           (fa_ext(ix,iy,0,0,is)/pbs(ix,iy,0)+
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0))/
     .           na_ext(ix,iy,is)
              else if (leftix(ix,iy).eq.-2) then
                ua_ext(ix,iy,is) =
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0)/
     .            na_ext(ix,iy,is)
              else if (rightix(ix,iy).eq.nx+1) then
                ua_ext(ix,iy,is) =
     .           fa_ext(ix,iy,0,0,is)/pbs(ix,iy,0)/na_ext(ix,iy,is)
              end if
            end do
          end do
         else
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2) then
                fa_ext(ix,iy,0,0,is)=((ua_ext(ix,iy,is)+
     .           ua_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           ((na_ext(ix,iy,is)+
     .           na_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           pbs(ix,iy,0)
                 fhi_ext(ix,iy,0,0) = fhi_ext(ix,iy,0,0) +
     .            1.5_R8*fa_ext(ix,iy,0,0,is)*(ta_ext(ix,iy,is)+
     .            ta_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
               else
                 fa_ext(ix,iy,0,0,is)=0.0_R8
               endif
            enddo
          enddo
          fa_ext(:,:,1,1,:) = 0.0_R8
         endif
cxpb
         ue_ext(-1:nx,-1:ny) = ue_ext(-1:nx,-1:ny) +
     .     za_ext(-1:nx,-1:ny,is) * ua_ext(-1:nx,-1:ny,is) *
     .     na_ext(-1:nx,-1:ny,is) *
     .     sign(1.0_R8,bb(-1:nx,-1:ny,0))
       end do
       do ix = -1, nx
         do iy = -1, ny
           if (ne_ext(ix,iy).gt.0.0_R8) then
             ue_ext(ix,iy) = ue_ext(ix,iy)/ne_ext(ix,iy)
           else
             ue_ext(ix,iy) = 0.0_R8
           endif
         enddo
       enddo
#endif
      end if

* ..linearise external sources
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (switch%b2stbm_linearisation.gt.0.0_R8) then
*  ...electron particle source
          t0 = st_ext%sne(iCv)
          if (t0.ge.0.0_R8) then
            srw%sne0(iCv,0) = t0
          else
            srw%sne0(iCv,1) = t0/dv%ne(iCv)
          end if
          srw%sne0(iCv,0) = srw%sne0(iCv,0)+
     .       switch%b2stbm_linearisation*abs(t0)
          srw%sne0(iCv,1) = srw%sne0(iCv,1)-
     .       switch%b2stbm_linearisation*abs(t0)/dv%ne(iCv)
          do is = 0, ns-1
*  ...heavy particles sources
            t0 = st_ext%sna(iCv,is)
            if (t0.ge.0.0_R8) then
              srw%sna0(iCv,0,is) = t0
            else
              srw%sna0(iCv,1,is) = t0/pl%na(iCv,is)
            end if
            srw%sna0(iCv,0,is) = srw%sna0(iCv,0,is)+
     .         switch%b2stbm_linearisation*abs(t0)
            srw%sna0(iCv,1,is) = srw%sna0(iCv,1,is)-
     .         switch%b2stbm_linearisation*abs(t0)/pl%na(iCv,is)
*  ...parallel momentum sources
            vthe=sqrt(pl%te(iCv)/me)
            vthe=vthe/sqrt(am(is)*mp/me)
            t0 = st_ext%smo(iCv,is)
            srw%smo0(iCv,0,is) = srw%smo0(iCv,0,is)+
     .         t0+switch%b2stbm_linearisation*abs(t0)*pl%ua(iCv,is)/vthe
            srw%smo0(iCv,1,is) = srw%smo0(iCv,1,is)-
     .         switch%b2stbm_linearisation*abs(t0)/vthe
          end do
*  ...electron energy source
          t0 = st_ext%she(iCv)
          if (t0.ge.0.0_R8) then
            srw%she0(iCv,0) = t0
          else
            srw%she0(iCv,1) = t0/pl%te(iCv)
          end if
          srw%she0(iCv,0) = srw%she0(iCv,0)+
     .       switch%b2stbm_linearisation*abs(t0)
          srw%she0(iCv,1) = srw%she0(iCv,1)-
     .       switch%b2stbm_linearisation*abs(t0)/pl%te(iCv)
*  ...ion energy source
          t0 = st_ext%shi(iCv)
          if (switch%b2stbm_internal_energy_sources.eq.0) then
            do is = 0, ns-1
              if (st_ext%smo(iCv,is).eq.0.0_R8 .and.
     .            st_ext%sna(iCv,is).eq.0.0_R8) cycle
              ua_tmp = (ua0(iCv,is)*na0(iCv,is) +
     .           st_ext%smo(iCv,is)*dtim*dtmo(is,ireg)/(am(is)*mp))/
     .          (na0(iCv,is)+st_ext%sna(iCv,is)*dtim*dtco(is,ireg))
              dua = ua_tmp - ua0(iCv,is)
              t0 = t0 - ua0(iCv,is)*st_ext%smo(iCv,is) + 
     .           0.5_R8*(am(is)*mp)*(ua0(iCv,is)**2)*
     .           st_ext%sna(iCv,is) -
     .          (na0(iCv,is)+st_ext%sna(iCv,is)*dtim*dtco(is,ireg))/
     .          (dtim*dtei(ireg))*0.5_R8*(am(is)*mp)*dua**2 
            end do
          end if
          if (t0.ge.0.0_R8) then
            srw%shi0(iCv,0) = t0
          else
            srw%shi0(iCv,1) = t0/pl%ti(iCv)
          end if
          srw%shi0(iCv,0) = srw%shi0(iCv,0)+
     .       switch%b2stbm_linearisation*abs(t0)
          srw%shi0(iCv,1) = srw%shi0(iCv,1)-
     .       switch%b2stbm_linearisation*abs(t0)/pl%ti(iCv)
*  ...electric charge source
          t0 = st_ext%sch(iCv)
          if (t0.ge.0.0_R8) then
            srw%sch0(iCv,0) = t0
          else
            srw%sch0(iCv,1) = t0/pl%po(iCv)
          end if
          srw%sch0(iCv,0) = srw%sch0(iCv,0)+
     .      switch%b2stbm_linearisation*abs(t0)
          srw%sch0(iCv,1) = srw%sch0(iCv,1)-
     .       switch%b2stbm_linearisation*abs(t0)/pl%po(iCv)
        else
*  ...raw sources: no linearisation
          srw%sne0(iCv,0) = st_ext%sne(iCv)
          srw%sna0(iCv,0,0:ns-1) = st_ext%sna(iCv,0:ns-1)
          srw%smo0(iCv,0,0:ns-1) = st_ext%smo(iCv,0:ns-1)
          srw%she0(iCv,0) = st_ext%she(iCv)
          srw%sch0(iCv,0) = st_ext%sch(iCv)
          srw%shi0(iCv,0) = st_ext%shi(iCv)
          if (switch%b2stbm_internal_energy_sources.eq.0) then
            do is = 0, ns-1
              if (st_ext%smo(iCv,is).eq.0.0_R8 .and.
     .            st_ext%sna(iCv,is).eq.0.0_R8) cycle
              ua_tmp = (ua0(iCv,is)*na0(iCv,is) +
     .           st_ext%smo(iCv,is)*dtim*dtmo(is,ireg)/(am(is)*mp))/
     .          (na0(iCv,is)+st_ext%sna(iCv,is)*dtim*dtco(is,ireg))
              dua = ua_tmp - ua0(iCv,is)
              srw%shi0(iCv,0) = srw%shi0(iCv,0) - 
     .           ua0(iCv,is)*st_ext%smo(iCv,is) + 
     .           0.5_R8*(am(is)*mp)*(ua0(iCv,is)**2)*
     .           st_ext%sna(iCv,is) -
     .          (na0(iCv,is)+st_ext%sna(iCv,is)*dtim*dtco(is,ireg))/
     .          (dtim*dtei(ireg))*0.5_R8*(am(is)*mp)*dua**2
            end do
          end if
        end if
      end do

      srw%b2stbm_sna(:,:)=srw%sna0(:,0,:)+srw%sna0(:,1,:)*pl%na(:,:)
      srw%b2stbm_smo(:,:)=srw%smo0(:,0,:)+srw%smo0(:,1,:)*pl%ua(:,:)
      srw%b2stbm_sne(:)=srw%sne0(:,0)+srw%sne0(:,1)*dv%ne(:)
      srw%b2stbm_she(:)=srw%she0(:,0)+srw%she0(:,1)*pl%te(:)
      srw%b2stbm_shi(:)=srw%shi0(:,0)+srw%shi0(:,1)*pl%ti(:)
      srw%b2stbm_sch(:)=srw%sch0(:,0)+srw%sch0(:,1)*pl%po(:)

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stbm_sna0to1 = srw%sna0
        b2stbm_smo0to3 = srw%smo0
        b2stbm_she0to3 = srw%she0
        b2stbm_shi0to3 = srw%shi0
      endif

* ..return
  999 if (main_call.or.ncall.eq.0) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbm

      end subroutine b2stbm
