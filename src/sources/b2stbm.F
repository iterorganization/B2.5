*-----------------------------------------------------------------------
*.specification

      subroutine b2stbm (nx, ny, ns, ismain, dtim,
     &  vol, gs, hx, hy, qz, qc, pbs, bb,
     &  na, ua, te, ti, po, ne, ni, fna, fhe, fhi, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0, sne0, main_call)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
      use b2mod_impgyro
      use b2mod_neutrals_namelist
      use b2mod_numerics_namelist
      use b2mod_sources
     1 , only : b2stbm_sna,b2stbm_smo,b2stbm_she,b2stbm_shi,b2stbm_sch,
     2          b2stbm_sne
      use b2mod_balance !djm Jan2017
     1 , only : b2stbm_sna0to1,b2stbm_smo0to3,b2stbm_she0to3,
     2          b2stbm_shi0to3, balance_netcdf
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain
      logical main_call
      real (kind=R8) :: dtim,
     *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1),
     &  gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &  fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
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
      integer n2, ia, is, is0, ix, iy, ixr, iyr, ireg
      real (kind=R8) :: t0, vthe, ua_tmp, dua
      real (kind=R8) :: n0_tmp, nr_tmp, nt_tmp, p0_tmp, pr_tmp, r0_tmp,
     *                  pe_tmp, re_tmp
      real (kind=R8) ::
     *  na_tot(-1:nx,-1:ny), na0(-1:nx,-1:ny,0:ns-1),
     *  ua0(-1:nx,-1:ny,0:ns-1), fna0(-1:nx,-1:ny,0:1,0:ns-1),
     *  fhi0(-1:nx,-1:ny,0:1), ti0(-1:nx,-1:ny),
     *  fhe0(-1:nx,-1:ny,0:1), fch0(-1:nx,-1:ny,0:1)
      integer, save :: ncall, impgyro_mod
      real (kind=R8), save :: linearisation
      logical, save :: internal_energy_sources
*   ..procedures
      real (kind=R8) ::
     * smin, smax
      external subini, subend, xertst, sfill,
     * smin, smax, b2xvsg, b2xvff, ipgeti, ipgetr
      intrinsic abs, mod, sqrt
*   ..initialisation
      data ncall/0/, linearisation/1.0_R8/, impgyro_mod/0/, 
     * internal_energy_sources/.false./

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
*   ..set internal parameters on first call
      if(ncall.eq.0) then !{
        call ipgetr ('b2stbm_linearisation', linearisation)
        if (linearisation.eq.1.0_R8)   ! for backward compatibility
     >   call ipgetr ('b2stbm_linearization', linearisation)
        call xertst 
     .   (linearisation.le.1.0_R8 .and. linearisation.ge.0.0_R8,
     .    'faulty input parameter linearisation')
        call ipgeti ('b2stbm_impgyro_mod', impgyro_mod)
      endif !}
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns, ismain
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     & 'faulty argument ismain')
*   ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test range of qz
       call xertst (-1.lt.smin(n2,qz(-1,-1,0),1).and.
     &   smax(n2,qz(-1,-1,0),1).lt.1, 'faulty argument range qz(,,0)')
       call xertst (0.lt.smin(n2,qz(-1,-1,1),1).and.
     &   smax(n2,qz(-1,-1,1),1).le.1, 'faulty argument range qz(,,1)')
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')

*   ..test edge values of fluxes
       do is = 0, ns-1
        call b2xvff (nx, ny, fna(-1,-1,0,is), 'fna')
       enddo
       call b2xvff (nx, ny, fhe, 'fhe')
       call b2xvff (nx, ny, fhi, 'fhi')
       call b2xvff (nx, ny, fch, 'fch')
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
      endif

*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0_R8, smo0, 1)
      call sfill (n2*2, 0.0_R8, sne0, 1)
      call sfill (n2*4, 0.0_R8, she0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
      call sfill (n2*4, 0.0_R8, sch0, 1)
      call sfill (n2*ns, 0.0_R8, b2stbm_sna, 1)
      call sfill (n2*ns, 0.0_R8, b2stbm_smo, 1)
      call sfill (n2, 0.0_R8, b2stbm_sne, 1)
      call sfill (n2, 0.0_R8, b2stbm_she, 1)
      call sfill (n2, 0.0_R8, b2stbm_shi, 1)
      call sfill (n2, 0.0_R8, b2stbm_sch, 1)

* ..obtain sources from IMPGYRO
      if (impgyro_mod.eq.0) goto 999

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
          do ix = -1, nx
            do iy = -1, ny
              if (bottomiy(ix,iy).eq.-2) cycle
              if (topiy(ix,iy).eq.ny+1) cycle
              ixr = rightix(ix,iy)
              iyr = rightiy(ix,iy)
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
                  n0_tmp = n0_tmp + dab2(ix+1,iy+1,ia,1)
                  nr_tmp = nr_tmp + dab2(ixr+1,iyr+1,ia,1)
                  nt_tmp = nt_tmp + dab2(ix+1,iy+1,ia,1)*
     .                              tab2(ix+1,iy+1,ia,1)
                  p0_tmp = p0_tmp + pfluxa(ix+1,iy+1,ia,1)
                  pe_tmp = pe_tmp + pefluxa(ix+1,iy+1,ia,1)
                  pr_tmp = pr_tmp + pfluxa(ixr+1,iyr+1,ia,1)
                  r0_tmp = r0_tmp + rfluxa(ix+1,iy+1,ia,1)
                  re_tmp = re_tmp + refluxa(ix+1,iy+1,ia,1)
                end if
              end do
              if (leftix(ix,iy).ne.-2 .and. n0_tmp.gt.0.0_R8 .and.
     .            rightix(ix,iy).ne.nx+1 .and. nr_tmp.gt.0.0_R8) then
                ua_tmp = (p0_tmp/n0_tmp+pr_tmp/nr_tmp)/2.0_R8
              elseif (leftix(ix,iy).ne.-2 .and. n0_tmp.gt.0.0_R8) then
                ua_tmp = p0_tmp/n0_tmp
              elseif (rightix(ix,iy).ne.nx+1 .and.
     .                                          nr_tmp.gt.0.0_R8) then
                ua_tmp = pr_tmp/nr_tmp
              else
                ua_tmp = 0.0_R8
              endif  
              na_tot(ix,iy) = na_tot(ix,iy) + n0_tmp
              na0(ix,iy,is) = na(ix,iy,is) + n0_tmp
              ti0(ix,iy) = ti0(ix,iy) + nt_tmp
              ua0(ix,iy,is) = (ua(ix,iy,is)*na(ix,iy,is)*
     .         sign(1.0_R8,bb(ix,iy,0)) + ua_tmp*n0_tmp)/na0(ix,iy,is)
              fna0(ix,iy,0,is) = 
     .         fna(ix,iy,0,is)/gs(ix,iy,0)/qc(ix,iy) + p0_tmp
              fna0(ix,iy,1,is) = fna(ix,iy,1,is)/gs(ix,iy,1) + r0_tmp
              fhi0(ix,iy,0) = fhi0(ix,iy,0) + pe_tmp
              fhi0(ix,iy,1) = fhi0(ix,iy,1) + re_tmp
            end do
          end do
        else
          na0(:,:,is) = na(:,:,is)
          ua0(:,:,is) = ua(:,:,is)*sign(1.0_R8,bb(:,:,0))
          fna0(:,:,0,is) = fna(:,:,0,is)/gs(:,:,0)/qc(:,:)
          fna0(:,:,1,is) = fna(:,:,1,is)/gs(:,:,1)
        end if
      end do
      do is = 1, nnmoli
        do ix = -1, nx
          do iy = -1, ny
            na_tot(ix,iy) = na_tot(ix,iy) + dmb2(ix+1,iy+1,is,1)
            ti0(ix,iy) = ti0(ix,iy) + 
     .       dmb2(ix+1,iy+1,is,1)*tmb2(ix+1,iy+1,is,1)
          end do
        end do
      end do
      do is = 1, nnioni
        do ix = -1, nx
          do iy = -1, ny
            na_tot(ix,iy) = na_tot(ix,iy) + dib2(ix+1,iy+1,is,1)
            ti0(ix,iy) = ti0(ix,iy) + 
     .       dib2(ix+1,iy+1,is,1)*tib2(ix+1,iy+1,is,1)
          end do
        end do
      end do
      do ix = -1, nx
        do iy = -1, ny
          ti0(ix,iy) = (ti0(ix,iy) + 
     .                  ti(ix,iy)*ni(ix,iy,0)) /
     .                 (na_tot(ix,iy) + ni(ix,iy,0))
          fhi0(ix,iy,0) = fhi0(ix,iy,0) + 
     .     fhi(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
          fhi0(ix,iy,1) = fhi0(ix,iy,1) + 
     .     fhi(ix,iy,1)/gs(ix,iy,1)
          fhe0(ix,iy,0) = fhe(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
          fhe0(ix,iy,1) = fhe(ix,iy,1)/gs(ix,iy,1)
          fch0(ix,iy,0) = fch(ix,iy,0)/gs(ix,iy,0)/qc(ix,iy)
          fch0(ix,iy,1) = fch(ix,iy,1)/gs(ix,iy,1)
        end do
      end do

      if (ncall.eq.0 .or.
     & (main_call .and. mod(ncall,impgyro_mod).eq.0)) then
cxpb The call below obtains the following quantities from IMPGYRO
cxpb Main W data: na_ext, ua_ext, ta_ext
cxpb Sources due to W: sne_ext, smo_ext, she_ext, shi_ext
cxpb Derived quantities: ne_ext, ne2_ext, ni_ext, fa_ext, fhi_ext, ue_ext
c        call read_b2mod_external(nx,ny,ns,ismain)
cxpb In the current arrangement, IMPGYRO provides: na_ext, ua_ext, ta_ext
cxpb It would be best to have instead as primary quantities:
cxpb na_ext, fa_ext, ta_ext, fhi_ext
cxpb And then derive: ua_ext

       sna_ext = 0.0_R8
       sne_ext = 0.0_R8
       smo_ext = 0.0_R8
       she_ext = 0.0_R8
       shi_ext = 0.0_R8
       sch_ext = 0.0_R8
       na_ext = 0.0_R8
       ua_ext = 0.0_R8
       ta_ext = 0.0_R8

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
         if (.false.) then
          do iy=-1,ny
            do ix=-1,nx
              if (leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1) then
                ua_ext(ix,iy,is) = 0.5_R8*
     .           (fa_ext(ix,iy,0,is)/pbs(ix,iy,0)+
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0))/
     .           na_ext(ix,iy,is)
              else if (leftix(ix,iy).eq.-2) then
                ua_ext(ix,iy,is) =
     .            fa_ext(rightix(ix,iy),rightiy(ix,iy),0,is)/
     .            pbs(rightix(ix,iy),rightiy(ix,iy),0)/
     .            na_ext(ix,iy,is)
              else if (rightix(ix,iy).eq.nx+1) then
                ua_ext(ix,iy,is) =
     .           fa_ext(ix,iy,0,is)/pbs(ix,iy,0)/na_ext(ix,iy,is)
              end if
            end do
          end do
         else
          do iy=-1,ny
            do ix=-1,nx
              if(leftix(ix,iy).ne.-2) then
                fa_ext(ix,iy,0,is)=((ua_ext(ix,iy,is)+
     .           ua_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           ((na_ext(ix,iy,is)+
     .           na_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8)*
     .           pbs(ix,iy,0)
                 fhi_ext(ix,iy,0) = fhi_ext(ix,iy,0) +
     .            1.5_R8*fa_ext(ix,iy,0,is)*(ta_ext(ix,iy,is)+
     .            ta_ext(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
               else
                 fa_ext(ix,iy,0,is)=0.0_R8
               endif
            enddo
          enddo
          fa_ext(:,:,1,:) = 0.0_R8
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
      end if

* ..linearise external sources
      do ix = 0, nx-1
        do iy = 0, ny-1
          ireg = region(ix,iy,0)
          if (linearisation.gt.0.0_R8) then
*  ...electron particle source
            t0 = sne_ext(ix,iy)
            if (t0.ge.0.0_R8) then
              sne0(ix,iy,0) = t0
            else
              sne0(ix,iy,1) = t0/ne(ix,iy)
            end if
            sne0(ix,iy,0) = sne0(ix,iy,0)+linearisation*abs(t0)
            sne0(ix,iy,1) = sne0(ix,iy,1)-
     .       linearisation*abs(t0)/ne(ix,iy)
            do is = 0, ns-1
*  ...heavy particles sources
              t0 = sna_ext(ix,iy,is)
              if (t0.ge.0.0_R8) then
                sna0(ix,iy,0,is) = t0
              else
                sna0(ix,iy,1,is) = t0/na(ix,iy,is)
              end if
              sna0(ix,iy,0,is) = sna0(ix,iy,0,is)+linearisation*abs(t0)
              sna0(ix,iy,1,is) = sna0(ix,iy,1,is)-
     .         linearisation*abs(t0)/na(ix,iy,is)
*  ...parallel momentum sources
              vthe=sqrt(te(ix,iy)/me)
              vthe=vthe/sqrt(am(is)*mp/me)
              t0 = smo_ext(ix,iy,is)*
     .         sign(1.0_R8,bb(ix,iy,0))
              smo0(ix,iy,0,is) = smo0(ix,iy,0,is)+
     .         t0+linearisation*abs(t0)*ua(ix,iy,is)/vthe
              smo0(ix,iy,1,is) = smo0(ix,iy,1,is)-
     .         linearisation*abs(t0)/vthe
            end do
*  ...electron energy source
            t0 = she_ext(ix,iy)
            if (t0.ge.0.0_R8) then
              she0(ix,iy,0) = t0
            else
              she0(ix,iy,1) = t0/te(ix,iy)
            end if
            she0(ix,iy,0) = she0(ix,iy,0)+linearisation*abs(t0)
            she0(ix,iy,1) = she0(ix,iy,1)-
     .       linearisation*abs(t0)/te(ix,iy)
*  ...ion energy source
            t0 = shi_ext(ix,iy)
            if (.not.internal_energy_sources) then
              do is = 0, ns-1
                if (smo_ext(ix,iy,is).eq.0.0_R8 .and.
     .              sna_ext(ix,iy,is).eq.0.0_R8) cycle
                ua_tmp = (ua0(ix,iy,is)*na0(ix,iy,is) +
     .           smo_ext(ix,iy,is)*dtim*dtmo(is,ireg)/(am(is)*mp))/
     .          (na0(ix,iy,is)+sna_ext(ix,iy,is)*dtim*dtco(is,ireg))
                dua = ua_tmp - ua0(ix,iy,is)
                t0 = t0 - ua0(ix,iy,is)*smo_ext(ix,iy,is) + 
     .           0.5_R8*(am(is)*mp)*(ua0(ix,iy,is)**2)*
     .           sna_ext(ix,iy,is) -
     .          (na0(ix,iy,is)+sna_ext(ix,iy,is)*dtim*dtco(is,ireg))/
     .          (dtim*dtei(ireg))*0.5_R8*(am(is)*mp)*dua**2 
              end do
            end if
            if (t0.ge.0.0_R8) then
              shi0(ix,iy,0) = t0
            else
              shi0(ix,iy,1) = t0/ti(ix,iy)
            end if
            shi0(ix,iy,0) = shi0(ix,iy,0)+linearisation*abs(t0)
            shi0(ix,iy,1) = shi0(ix,iy,1)-
     .       linearisation*abs(t0)/ti(ix,iy)
*  ...electric charge source
            t0 = sch_ext(ix,iy)
            if (t0.ge.0.0_R8) then
              sch0(ix,iy,0) = t0
            else
              sch0(ix,iy,1) = t0/po(ix,iy)
            end if
            sch0(ix,iy,0) = sch0(ix,iy,0)+linearisation*abs(t0)
            sch0(ix,iy,1) = sch0(ix,iy,1)-
     .       linearisation*abs(t0)/po(ix,iy)
          else
*  ...raw sources: no linearisation
            sne0(ix,iy,0) = sne_ext(ix,iy)
            sna0(ix,iy,0,0:ns-1) = sna_ext(ix,iy,0:ns-1)
            smo0(ix,iy,0,0:ns-1) = smo_ext(ix,iy,0:ns-1)*
     .         sign(1.0_R8,bb(ix,iy,0))
            she0(ix,iy,0) = she_ext(ix,iy)
            sch0(ix,iy,0) = sch_ext(ix,iy)
            shi0(ix,iy,0) = shi_ext(ix,iy)
            if (.not.internal_energy_sources) then
              do is = 0, ns-1
                if (smo_ext(ix,iy,is).eq.0.0_R8 .and.
     .              sna_ext(ix,iy,is).eq.0.0_R8) cycle
                ua_tmp = (ua0(ix,iy,is)*na0(ix,iy,is) +
     .           smo_ext(ix,iy,is)*dtim*dtmo(is,ireg)/(am(is)*mp))/
     .          (na0(ix,iy,is)+sna_ext(ix,iy,is)*dtim*dtco(is,ireg))
                dua = ua_tmp - ua0(ix,iy,is)
                shi0(ix,iy,0) = shi0(ix,iy,0) - 
     .           ua0(ix,iy,is)*smo_ext(ix,iy,is) + 
     .           0.5_R8*(am(is)*mp)*(ua0(ix,iy,is)**2)*
     .           sna_ext(ix,iy,is) -
     .          (na0(ix,iy,is)+sna_ext(ix,iy,is)*dtim*dtco(is,ireg))/
     .          (dtim*dtei(ireg))*0.5_R8*(am(is)*mp)*dua**2
              end do
            end if
          end if
        end do
      end do

      b2stbm_sna(:,:,:)=sna0(:,:,0,:)+sna0(:,:,1,:)*na(:,:,:)
      b2stbm_smo(:,:,:)=smo0(:,:,0,:)+smo0(:,:,1,:)*ua(:,:,:)
      b2stbm_sne(:,:)=sne0(:,:,0)+sne0(:,:,1)*ne(:,:)
      b2stbm_she(:,:)=she0(:,:,0)+she0(:,:,1)*te(:,:)
      b2stbm_shi(:,:)=shi0(:,:,0)+shi0(:,:,1)*ti(:,:)
      b2stbm_sch(:,:)=sch0(:,:,0)+sch0(:,:,1)*po(:,:)

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stbm_sna0to1 = sna0
        b2stbm_smo0to3 = smo0
        b2stbm_she0to3 = she0
        b2stbm_shi0to3 = shi0
      endif

* ..return
  999 if (main_call.or.ncall.eq.0) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbm

      end subroutine b2stbm
