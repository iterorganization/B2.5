      subroutine b2usrtrc(nx,ny,ns,BoRiS)
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_time
      use b2mod_rates
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_boundary_namelist
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_eirbra
#endif

c  version : 03.04.2001 15:55

c=======================================================================
c*** Puts user-specific quantities on the output file
c*** "tracing/user.trc"
c-----------------------------------------------------------------------
      implicit none

#     include <DIMENSIONS.F>
#ifndef B25_EIRENE
#     include <COUPLE/EIRBRA.F>
#endif

      integer nx, ny, ns
      real (kind=R8) :: BoRiS
c-----------------------------------------------------------------------
c*** lv    :  version number to appear in the file
c*** lvmj  :  version number referring to the last major change

      integer lv,lvmj
      parameter (lv=20000727, lvmj=970212)

c-----------------------------------------------------------------------
c*** User-specific data:
c***
c***      header     var       Legend
c***                                                           (global)
c***     Pneu_PFR  hlp_ppfr   neutral pressure in PFR - mean
c***     Pneu_PFi  hlp_ppfi   neutral pressure in PFR - inner
c***     Pneu_PFo  hlp_ppfo   neutral pressure in PFR - outer
c***     He_fract  hlp_hefr   helium fraction in the pump
c***     He_enrch  hlp_henr   helium "enrichment"
c***     He_p_prs  hlp_hepr   helium partial pressure in PFR
c***     pk_pwr_i  hlp_ppwi   peak power load onto the inner target
c***     pk_pwr_o  hlp_ppwo   peak power load onto the outer target
c***     pk_pw_iu  hlp_pwiu   peak power load on the inner upper target
c***     pk_pw_ou  hlp_pwou   peak power load on the outer upper target
c***                          (the peaks do not include radiation load)
c***     Sp_DT     hlp_spdt   effective pumping speed for DT
c***     Sp_He     hlp_sphe   effective pumping speed for He
c***     Sp_He_sc  hlp_sshe   re-scaled pumping speed for He
c***     He_src    hlp_hesr   helium source
c***     He_sc_cr  hlp_hesc   rescaled average nHe/ne at the CEI
c***     He_sc_sp  hlp_hess   rescaled average nHe/ne at the separatrix
c***     He_og_cr  hlp_heoc   original average nHe/ne at the CEI
c***     He_og_sp  hlp_heos   original average nHe/ne at the separatrix
c***     He_enrcs  hlp_hens   helium enrichment to ns from throughputs
c***     He_enrcc  hlp_henc   helium enrichment to nc from throughputs
c***     DT_thrpt  hlp_dttp   DT throughput
c***     He_thrpt  hlp_hetp   He throughput
c***                                        (PFR - by group and species)
c***     PFRflx    hlp_pfrf   net nuclei flux through each surface group
c-----------------------------------------------------------------------

c*** Global values
      integer nglb_tmp, ng_tmp
      parameter (nglb_tmp=22)
      real (kind=R8) ::
     *     equg_tmp(nglb_tmp)
      real (kind=R8) ::
     *     hlp_ppfr, hlp_ppfi, hlp_ppfo, hlp_hefr, hlp_henr, hlp_hepr,
     ,     hlp_ppwi, hlp_ppwo, hlp_pwiu, hlp_pwou, hlp_spdt, hlp_sphe,
     ,     hlp_sshe, hlp_hesr, hlp_hesc, hlp_hess, hlp_heoc, hlp_heos,
     ,     hlp_hens, hlp_henc, hlp_dttp, hlp_hetp
      equivalence
     ,     (equg_tmp(1),  hlp_ppfr), (equg_tmp(2),  hlp_ppfi),
     ,     (equg_tmp(3),  hlp_ppfo), (equg_tmp(4),  hlp_hefr),
     ,     (equg_tmp(5),  hlp_henr), (equg_tmp(6),  hlp_hepr),
     ,     (equg_tmp(7),  hlp_ppwi), (equg_tmp(8),  hlp_ppwo),
     ,     (equg_tmp(9),  hlp_pwiu), (equg_tmp(10), hlp_pwou),
     ,     (equg_tmp(11), hlp_spdt), (equg_tmp(12), hlp_sphe),
     ,     (equg_tmp(13), hlp_sshe), (equg_tmp(14), hlp_hesr),
     ,     (equg_tmp(15), hlp_hesc), (equg_tmp(16), hlp_hess),
     ,     (equg_tmp(17), hlp_heoc), (equg_tmp(18), hlp_heos),
     ,     (equg_tmp(19), hlp_hens), (equg_tmp(20), hlp_henc),
     ,     (equg_tmp(21), hlp_dttp), (equg_tmp(22), hlp_hetp)
      character hlp_hg(nglb_tmp)*8
c*** PFR flow values
      integer npff_tmp, np_tmp, ns_tmp, na_tmp, lnh_tmp, i0
      parameter (npff_tmp=1, lnh_tmp=20)
      real (kind=R8) ::
     *     equp_tmp(0:npfrgrpd,DEF_NATM,npff_tmp)
      real (kind=R8) ::
     *     hlp_pfrf(0:npfrgrpd,DEF_NATM,npff_tmp)
      equivalence (equp_tmp(0,1,1), hlp_pfrf)
      character hlp_hp(npff_tmp)*8,hlp_hb*8,
     ,  hlp_ha(0:npfrgrpd,ns,npff_tmp)*(lnh_tmp)
c-----------------------------------------------------------------------

      integer lun, ltrc
      parameter (ltrc=ntrc_usr)
      real (kind=R8) ::
     *     hlp, hlp_p, hlp_s, hlp_pi, hlp_si, hlp_po, hlp_so, hlp_tot,
     ,     hlp_he, hlp_hep, hlp_h, hlp_neu_peak, s, hlpnw,
     ,     hlp_ne_he_c, hlp_ne_he_s
      integer i, j, k, l, m, n, im, is, ix, iy, iyp, js, ls, nv, j_he_i
      integer jxi, jxa, jsep
      logical ex, warn_tmp, first, double_null
      character tab, hlp_c, f*40, hlp_frm*120, hlp_frf*240, hlp_tpl*2
      data hlp_neu_peak / 5.0_R8 /
      data first / .true. /, double_null / .false. /
c-----------------------------------------------------------------------
c*** Column headers:
      data hlp_hg /
     1      'Pneu_PFR', 'Pneu_PFi', 'Pneu_PFo', 'He_fract', 'He_enrch',
     2      'He_p_prs', 'pk_pwr_i', 'pk_pwr_o', 'pk_pw_iu', 'pk_pw_ou',
     3      'Sp_DT   ', 'Sp_He   ', 'Sp_He_sc', 'He_src  ', 'He_sc_cr',
     4      'He_sc_sp', 'He_og_cr', 'He_og_sp', 'He_enrcs', 'He_enrcc',
     4      'DT_thrpt', 'He_thrpt'/
      data hlp_hp /'PFR_flow'/
c-----------------------------------------------------------------------

      save
     l   first, double_null, lun,
     i   j_he_i, nv, jsep, ng_tmp,
     c   hlp_tpl, hlp_hb, hlp_frm
c=======================================================================
      if(.not.filedata(ltrc)) return
c=======================================================================

      if (first) then !{

c*** Check the topology

        hlp_tpl='UD'
        if(nncut.eq.1) then !{
c*** Single null
          hlp_tpl='SN'
        else if(nncut.eq.2) then !}{
          if(nniso.gt.0) then !{
            do i=1,nniso !{
              if(nyiso1(i).le.-1 .and. nyiso2(i).ge.ny) then !{
c*** Double null
                hlp_tpl='DN'
              end if !}
            end do !}
          end if !}
        else if(nncut.eq.0) then !}{
c*** Limiter case
          hlp_tpl='LM'
        end if !}
        if(hlp_tpl.eq.'SN') then !{
          double_null=.false.
        else if(hlp_tpl.eq.'DN') then !}{
          double_null=.true.
        else !}{
          write(0,*) 'b2usrtrc: only works for SN or DN topology.'
          filedata(ltrc)=.false.
          return
        end if !}

c*** Check the input data for consistency and prepare indices

        j_he_i=0
        do j=1,nspecies !{
          if(nfluids(j).eq.1) then !{
            if(j_h_at.le.0 .or. j_h_at.gt.nnatmi) then !{
c*** Assume the same order of species in B2 and Eirene
              j_h_at=j
            end if !}
          end if !}
          if(nfluids(j).eq.2) then !{
            j_he_i=eb2spcr(j)
            if(j_he_at.le.0 .or. j_he_at.gt.nnatmi) then !{
c*** Assume the same order of species in B2 and Eirene
              j_he_at=j
            end if !}
          end if !}
        end do !}

        print *,'%% b2usrtr: j_he_at,j_h_at,j_he_i,double_null = ', !###
     ,                       j_he_at,j_h_at,j_he_i,double_null      !###

        first=.false.
      end if !}
c=======================================================================
c*** Initialise the arrays

      equg_tmp=0.0_R8
      equp_tmp=0.0_R8

c*** Re-scale the neutral fluxes

      do i=1,nnlimi+nnstsi !{
        do j=1,nnatmi !{
          wldna(i,j,0)=0.0_R8
          wldra(i,j,0)=0.0_R8
          wldpa(i,j,0)=0.0_R8
          do k=1,nstrai !{
            wldna(i,j,0)=wldna(i,j,0)+wldna(i,j,k)*sclstra(j,k)
            wldra(i,j,0)=wldra(i,j,0)+wldra(i,j,k)*sclstra(j,k)
            wldpa(i,j,0)=wldpa(i,j,0)+wldpa(i,j,k)*sclstra(j,k)
          end do !}
        end do !}
        do j=1,nnmoli !{
          l=lmolscl(j)
          wldnm(i,j,0)=0.
          wldrm(i,j,0)=0.
          wldpm(i,j,0)=0.
          do k=1,nstrai !{
            if(l.gt.0) then !{
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)*sclstra(l,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)*sclstra(l,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)*sclstra(l,k)
            else !}{
              wldnm(i,j,0)=wldnm(i,j,0)+wldnm(i,j,k)
              wldrm(i,j,0)=wldrm(i,j,0)+wldrm(i,j,k)
              wldpm(i,j,0)=wldpm(i,j,0)+wldpm(i,j,k)
            end if !}
          end do !}
        end do !}
      end do !}
c=======================================================================
      if(j_he_at.gt.0) then !{

c*** Calculate the helium-related data. First, average densities
c*** for the CEI...

        hlp_pi=0.0_R8
        hlp_he=0.0_R8
        hlp_si=0.0_R8
        hlp_ne_he_c=0.0_R8

c*** ... and separatrix
        hlp_po=0.0_R8
        hlp_hep=0.0_R8
        hlp_so=0.0_R8
        hlp_ne_he_s=0.0_R8
        iy=jsep+1

c      print *,'%% b2usrtrc: j_he_i,iy = ',j_he_i,iy !###

        do ix=0,nx-1 !{
          if (mod(region(ix,iy-1,0),4).eq.1) then  !{
            hlp_si=hlp_si+gs(ix,0,1)
            hlp_pi=hlp_pi+gs(ix,0,1)*ne(ix,0)
            hlp=0.0_R8
            do j=1,2 !{
              hlp=hlp+na(ix,0,j_he_i+j)
              hlp_ne_he_c=hlp_ne_he_c+na(ix,0,j_he_i+j)*j*gs(ix,0,1)
            end do !}
            hlp_he=hlp_he+gs(ix,0,1)*hlp

            hlp_so=hlp_so+gs(ix,iy,1)
            hlp_po=hlp_po+gs(ix,iy,1)*ne(ix,iy)
            hlp=0.0_R8
            do j=1,2 !{
              hlp=hlp+na(ix,iy,j_he_i+j)
              hlp_ne_he_s=hlp_ne_he_s+na(ix,iy,j_he_i+j)*j*gs(ix,iy,1)
            end do !}
            hlp_hep=hlp_hep+gs(ix,iy,1)*hlp
          end if !}
        end do !}
        if(hlp_pi.gt.0.0_R8) then !{
          hlp_heoc=hlp_he/hlp_pi
          hlp_ne_he_c=hlp_ne_he_c/hlp_pi
        else !}{
          hlp_heoc=0.0_R8
          hlp_ne_he_c=0.0_R8
        end if !}
        if(hlp_po.gt.0.0_R8) then !{
          hlp_heos=hlp_hep/hlp_po
          hlp_ne_he_s=hlp_ne_he_s/hlp_po
        else !}{
          hlp_heos=0.0_R8
          hlp_ne_he_s=0.0_R8
        end if !}
      end if !}
c=======================================================================
      if(lpfrs_pmp.gt.0 .and. lpfrs_pmp.le.3) then !{

c-----------------------------------------------------------------------
c*** Calculate the pumping-related data.
c***
c*** First, estimate the neutral pressure in the duct - as the average
c*** over the plasma cells adjacent to the duct, weighted with the
c*** surface area. The pump location is specified by parameters
c***
c***  lpfrs_pmp:  1 - in PFR, 2 - on the outer wall, 3 - in upper PFR
c***  lpfrb_[io]: cell indices corresponding to the bottom (PFR) or
c***              just beginning of the region open to the duct
c***  lpfrt_[io]: cell indices corresponding to the top (PFR) or just
c***              end of the region open to the duct
c***
c*** In case of PFR, the bottom cells are counted back from the targets
c*** and the top cells down from the x-point if the corresponding value
c*** is negative.
c*** In case of outer wall, negative values of lpfr[tb]_o mean counting
c*** back from the last (in terms of B2) target; lpfr[tb]_i is not used
c***
c***  j_he_at  is the helium index in the Eirene atoms list. If 0, then
c***           the B2 species is taken assuming the same order.
c***  j_h_at   is the same for the hydrogen species, and
c***  l_h_mol  gives the number of hydrogen nuclei in each molecule
c***  lhetrgts is the list on Eirene surfaces used to estimate the
c***           He enrichment from the single-sided fluxes in Eirene.
c-----------------------------------------------------------------------

        hlp_p=0.0_R8
        hlp_s=0.0_R8
        hlp_hep=0.0_R8

c      print *,'%%## j_he_at,nnatmi,lpfrb_i,lpfrt_i,lpfrb_o,lpfrt_o',!###
c     ,              j_he_at,nnatmi,lpfrb_i,lpfrt_i,lpfrb_o,lpfrt_o  !###

        hlp_pi=0.0_R8
        hlp_si=0.0_R8
        if(lpfrs_pmp.ne.2) then !{

c*** Inner divertor

          m=lpfrb_i
          if(m.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              m=1-m
            else !}{
              m=nxiso1(1)-1+m
            end if !}
          end if !}
          n=lpfrt_i
          if(n.le.0) then !{
            if(lpfrs_pmp.eq.1) then !{
              n=leftcut(1)-1+n
            else !}{
              n=rightcut(1)-1+n
            end if !}
          end if !}
          do j=m,n !{
            hlp_si=hlp_si+gs(j,0,1)
            hlp=0.0_R8
            do k=1,nnatmi !{
              hlp=hlp+dab2(j+1,1,k,1)*tab2(j+1,1,k,1)
              if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +          dab2(j+1,1,k,1)*tab2(j+1,1,k,1)*gs(j,0,1)
            end do !}
            do k=1,nnmoli !{
              hlp=hlp+dmb2(j+1,1,k,1)*tmb2(j+1,1,k,1)
            end do !}
            hlp_pi=hlp_pi+hlp*gs(j,0,1)
          end do !}

c          print '(a,1p,2i4,20e12.3)',' %%## m,n,hlp_si,hlp_hep: ',  !###
c     ,                                      m,n,hlp_si,hlp_hep      !###

        end if !}

c*** Outer divertor

        hlp_po=0.0_R8
        hlp_so=0.0_R8
        m=lpfrt_o
        if(m.le.0) then !{
          if(lpfrs_pmp.eq.1) then !{
            m=rightcut(1)-m
          else if(lpfrs_pmp.eq.2) then !}{
            m=nx+m
          else !}{
            m=rightcut(1)-m
          end if !}
        end if !}
        n=lpfrb_o
        if(n.le.0) then !{
          if(lpfrs_pmp.eq.1) then !{
            n=nx+n
          else if(lpfrs_pmp.eq.2) then !}{
            n=nx+n
          else !}{
            n=nxiso1(1)+2+m
          end if !}
        end if !}
        iy=0
        iyp=1
        if(lpfrs_pmp.eq.2) then !{
          iy=ny
          iyp=ny
        endif !}
        do j=m,n !{
          hlp_so=hlp_so+gs(j,iy,1)
          hlp=0.0_R8
          do k=1,nnatmi !{
            hlp=hlp+dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)
            if(k.eq.j_he_at) hlp_hep=hlp_hep+
     +            dab2(j+1,iyp,k,1)*tab2(j+1,iyp,k,1)*gs(j,iy,1)
          end do !}
          do k=1,nnmoli !{
            hlp=hlp+dmb2(j+1,iyp,k,1)*tmb2(j+1,iyp,k,1)
          end do !}
          hlp_po=hlp_po+hlp*gs(j,iy,1)
        end do !}

        hlp_s=hlp_si+hlp_so
        hlp_p=hlp_pi+hlp_po
c<<<
c      print '(a,1p,2i4,20e12.3)',' %%## m,n,hlp_so,hlp_hep,hlp_p: ',
c     ,                                  m,n,hlp_so,hlp_hep,hlp_p
c>>>
        if(hlp_s.gt.0.0_R8) then !{
          hlp_ppfr=hlp_p/hlp_s
          hlp_hepr=hlp_hep/hlp_s
        end if !}
        if(hlp_si.gt.0.0_R8) hlp_ppfi=hlp_pi/hlp_si
        if(hlp_so.gt.0.0_R8) hlp_ppfo=hlp_po/hlp_so

c*** Calculate throughputs, pumping speed, and re-scale He if necessary

        ls=j_h_at
        do l=1,2 !{
          if(ls.gt.0) then !{

c*** Neutral flux across the side walls

            hlpnw=0.0_R8
            do js=1,nnlimi !{
              hlpnw=hlpnw+wldna(js,ls,0)-wldra(js,ls,0)
              do j=1,nnmoli !{
                if(mlcmp(ls,j).gt.0) hlpnw=hlpnw+
     +                        mlcmp(ls,j)*(wldnm(js,j,0)-wldrm(js,j,0))
              end do !}
            end do !}
          end if !}
          equg_tmp(20+l)=hlpnw/2.7e20_R8
          ls=j_he_at
        end do !}
        hlp_dttp=0.5_R8*hlp_dttp     ! assume molecular DT pumping
        if(hlp_ppfr.gt.0.0_R8) hlp_spdt=hlp_dttp/hlp_ppfr
        if(hlp_hepr.gt.0.0_R8) hlp_sphe=hlp_hetp/hlp_hepr
        hlp_sshe=spmp_he_to_d*hlp_spdt

c*** Calculate the helium-related data

c       print *,'%%## j_he_at,nhetrgtx,fusion_power =',  !###
c     ,                            j_he_at,nhetrgtx,fusion_power  !###
        if(fusion_power.gt.0.0_R8) then !{

c*** Re-scale the helium density upstream to meet both throughput and
c*** pumping speed

          hlp_hesr=3.55e17_R8/2.7e20_R8*fusion_power
          hlp=hlp_sshe*hlp_hepr
          if(hlp.gt.0.0_R8) then !{
            hlp=hlp_hesr/hlp
            if(hlp_ne_he_c.lt.1.0_R8) then !{
              hlp_hesc=
     .         1.0_R8-1.0_R8/(1.0_R8+hlp_heoc*hlp/(1.0_R8-hlp_ne_he_c))
            else !}{
              hlp_hesc=1
            end if !}
            if(hlp_ne_he_s.lt.1.0_R8) then !{
              hlp_hess=
     .         1.0_R8-1.0_R8/(1.0_R8+hlp_heos*hlp/(1.0_R8-hlp_ne_he_s))
            else !}{
              hlp_hess=1.0_R8
            end if !}
          end if !}
          hlp=2.0_R8*hlp_dttp*hlp_hesc
          if(hlp.gt.0.0_R8) hlp_henc=hlp_hesr/hlp
          hlp=2.0_R8*hlp_dttp*hlp_hess
          if(hlp.gt.0.0_R8) hlp_hens=hlp_hesr/hlp
        end if !}
      end if !}

      if(j_he_at.gt.0) then !{

c*** Helium fraction and related data - flux estimation

        hlp_tot=0.0_R8
        hlp_he=0.0_R8
        hlp_h=0.0_R8
        do i=1,nhetrgtx !{
          k=lhetrgts(i)
          if(k.eq.0) go to 10
          if(k.lt.0) k=nnlimi-k
          hlp_he=hlp_he+wldna(k,j_he_at,0)
          hlp_h=hlp_h+wldna(k,j_h_at,0)
          do j=1,nnatmi !{
            hlp_tot=hlp_tot+wldna(k,j,0)
          end do !}
          do j=1,nnmoli !{
            hlp_tot=hlp_tot+wldnm(k,j,0)
            hlp_h=hlp_h+wldnm(k,j,0)*l_h_mol(j)
          end do !}
        end do !}
 10     if(hlp_h.gt.0.0_R8) hlp_hefr=hlp_he/(hlp_h+hlp_he)
        if(hlp_heoc.gt.0.0_R8) hlp_henr=hlp_hefr/hlp_heoc
      end if !}

c*** Calculate the peak powers to the targets
c***  (assuming a peaking factor of hlp_neu_peak for neutrals)

      hlp_ppwi=0.0_R8
      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        s=gs(0,iy,0)
        hlp=abs(fhe(0,iy,0)+fhi(0,iy,0))
        do is=0,ns-1 !{
          hlp=hlp+abs(fna(0,iy,0,is)*(
     .     (rpt(0,iy,is)*hx(-1,iy)+rpt(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(0,iy,is)*hx(-1,iy)+kinrgy(-1,iy,is)*hx(0,iy))/
     .     (hx(-1,iy)+hx(0,iy))))
        end do !}
        hlp=hlp/s
        hlp_ppwi=max(hlp_ppwi,hlp)
        hlp_s=hlp_s+s
      end do !}
      js=ltns(1)+nnlimi
      hlp_ppwi=hlp_ppwi+(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s

      hlp_ppwo=0.0_R8
      hlp_s=0.0_R8
      do iy=0,ny-1 !{
        hlp=abs(fhe(nx,iy,0)+fhi(nx,iy,0))
        s=gs(nx,iy,0)
        do is=0,ns-1 !{
          hlp=hlp+abs(fna(nx,iy,0,is)*(
     .     (rpt(nx,iy,is)*hx(nx-1,iy)+rpt(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy))*ev+(1.0_R8-BoRiS)*
     .     (kinrgy(nx,iy,is)*hx(nx-1,iy)+kinrgy(nx-1,iy,is)*hx(nx,iy))/
     .     (hx(nx-1,iy)+hx(nx,iy))))
        end do !}
        hlp=hlp/s
        hlp_ppwo=max(hlp_ppwo,hlp)
        hlp_s=hlp_s+s
      end do !}
      js=ltns(nncut)+nnlimi
      hlp=(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s
      hlp_ppwo=hlp_ppwo+(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s

      if(double_null) then !{
        hlp_pwou=0.0_R8
        ix=nxiso2(1)+1
        hlp_s=0.0_R8
        do iy=0,ny-1 !{
          s=gs(ix,iy,0)
          hlp=abs(fhe(ix,iy,0)+fhi(ix,iy,0))
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(ix,iy,0,is)*(
     .       (rpt(ix,iy,is)*hx(ix-1,iy)+
     .        rpt(ix-1,iy,is)*hx(ix,iy))/
     .       (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(ix,iy,is)*hx(ix-1,iy)+
     .        kinrgy(ix-1,iy,is)*hx(ix,iy))/
     .       (hx(ix-1,iy)+hx(ix,iy))))
          end do !}
          hlp=hlp/s
          hlp_pwou=max(hlp_pwou,hlp)
          hlp_s=hlp_s+s
        end do !}
        js=ltns(3)+nnlimi
        hlp_pwou=hlp_pwou+(wldnek(js,0)+
     +                                 wldnep(js,0))*hlp_neu_peak/hlp_s

        hlp_pwiu=0.0_R8
        ix=nxiso1(1)
        hlp_s=0.0_R8
        do iy=1,ny !{
          hlp=abs(fhe(ix,iy,0)+fhi(ix,iy,0))
          s=gs(ix,iy,0)
          do is=0,ns-1 !{
            hlp=hlp+abs(fna(ix,iy,0,is)*(
     .       (rpt(ix,iy,is)*hx(ix-1,iy)+
     .        rpt(ix-1,iy,is)*hx(ix,iy))/
     .       (hx(ix-1,iy)+hx(ix,iy))*ev+(1.0_R8-BoRiS)*
     .       (kinrgy(ix,iy,is)*hx(ix-1,iy)+
     .        kinrgy(ix-1,iy,is)*hx(ix,iy))/
     .       (hx(ix-1,iy)+hx(ix,iy))))
          end do !}
          hlp=hlp/s
          hlp_pwiu=max(hlp_pwiu,hlp)
          hlp_s=hlp_s+s
        end do !}
        js=ltns(2)+nnlimi
        hlp=(wldnek(js,0)+wldnep(js,0))*hlp_neu_peak/hlp_s
        hlp_pwiu=hlp_pwiu+(wldnek(js,0)+
     +                                 wldnep(js,0))*hlp_neu_peak/hlp_s
      end if !}

c*** Calculate the neutral flows in PFR 

      do k=1,npfrgrp !{
        do j=1,nspecies !{
          hlp=1.0_R8/2.7e20_R8
          if(atoms_list(j).eq.'H' .or. atoms_list(j).eq.'D' .or.
     ,                                atoms_list(j).eq.'T') hlp=0.5*hlp
          do i=ipfrgrp(k),jpfrgrp(k) !{
            m=lpfrgrp(i)
            if(m.lt.0) m=nnlimi-m
            if(isrftype(m).ne.-3) then !{
              write (*,*) '#b2usrtrc: wrong type of diagnostic ',
     ,                  'surface. Surface:',m,',  type: ',isrftype(m)
              write (*,*) 'Check the DG model!!!'
              hlp_pfrf(k,j,1)=0.0_R8
              go to 30
            end if !}
            hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)+wldna(m,j,0)
            do im=1,nnmoli !{
              hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)+
     +                                      wldnm(m,im,0)*mlcmp(j,im)
            end do !}
          end do !}
          hlp_pfrf(k,j,1)=hlp_pfrf(k,j,1)*hlp
        end do !}
      end do !}
 30   do k=1,npfrgrp !{
        do j=1,nspecies !{
          hlp_pfrf(0,j,1)=hlp_pfrf(0,j,1)+hlp_pfrf(k,j,1)
        end do !}
      end do !}
c=======================================================================
      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then !{

c*** Write the data in a separate file for plotting...

        if(npfrgrp.gt.0) then !{
          i0=0
          if(npfrgrp.le.1) i0=1
          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp),
     ,      (((equp_tmp(i,j,k),i=i0,npfrgrp),j=1,nspecies),k=1,npff_tmp)
        else !}{
          write(lun,hlp_frm,err=94) tim,(equg_tmp(j),j=1,ng_tmp)
        end if !}
      end if !}
      return
c=======================================================================

      entry b2usrtri(nx,ny,ns,BoRiS)
c-----------------------------------------------------------------------
c*** Initialise the file for tracing data
c-----------------------------------------------------------------------

      tab=char(9)
      call get_jsep(nx,ny,jxi,jxa,jsep)

c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name

      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      i0=0                                                           !srv 05.11.10
      if(npfrgrp.le.1) i0=1                                          !srv 05.11.10

      inquire(file=f,exist=ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) then !{
        na_tmp=0
        np_tmp=0
        ns_tmp=0
        read(lun,700,err=92) nv,ng_tmp
        lv_trc(ltrc)=nv
        if(nv.ge.20000727) then !{
          rewind(lun)
          read(lun,700,err=92) nv,ng_tmp,np_tmp,ns_tmp,na_tmp
        end if !}
        warn_tmp=.false.
        if(nv.ne.lv) then !{
          warn_tmp=.true.
          write(*,*) '!!! b2usrtrc: current version ',lv,
     ,                         '  differs from that in the file ',f
          if(nv .lt. lvmj) then !{
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No tracing file - ',
     ,                              'remove the existing file first: ',f
            filedata(ltrc)=.false.
            return
          end if !}
          if(nv .lt. 20000727 .and. npfrgrp.gt.0) then !{
            write(*,*) 'No tracing file - ',
     ,                              'remove the existing file first: ',f
            filedata(ltrc)=.false.
            return
          end if !}
        end if !}
        if(na_tmp.ne.npfrgrp) then !{
          write(*,*) '!!! b2usrtrc: current number of PFR surface ',
     ,          'groups ',npfrgrp,'  differs from that in the file ',f
          filedata(ltrc)=.false.
        end if !}
        if(npfrgrp.gt.0) then !{
          if(ng_tmp.ne.nglb_tmp) then !{
            write(*,*) '!!! b2usrtrc: current number of global values ',
     ,                  nglb_tmp,'  differs from that in the file ',f
            filedata(ltrc)=.false.
          end if !}
          if(np_tmp.ne.npff_tmp) then !{
            write(*,*) '!!! b2usrtrc: current number of PFR values ',
     ,                  npff_tmp,'  differs from that in the file ',f
            filedata(ltrc)=.false.
          end if !}
          if(ns_tmp.ne.nspecies) then !{
            write(*,*) '!!! b2usrtrc: current number of species ',
     ,                  nspecies,'  differs from that in the file ',f
            filedata(ltrc)=.false.
          end if !}
        end if !}

        if(.not.filedata(ltrc)) then !{
          write(*,*) 'No tracing file - ',
     ,                                'remove the existing file first!'
          return
        end if !}
        if(warn_tmp) write(*,*) '  You should remove the file!'
c*** Read to the end of the file
 50     read(lun,'(a1)',end=55) hlp_c
        go to 50
 55     backspace lun
      else !}{
c*** New file
        ng_tmp=nglb_tmp
        np_tmp=npff_tmp
        ns_tmp=nspecies
        na_tmp=npfrgrp
        do i=0,npfrgrp !{
          hlp_hb='total'
          if(i.gt.0) hlp_hb=gpfrgrp(i)
          do j=1,nspecies !{
            do k=1,npff_tmp !{
              hlp_ha(i,j,k)=hlp_hp(k)//'_'//hlp_hb//'_'//atoms_list(j)
              call lefta(hlp_ha(i,j,k),lnh_tmp,l)
            end do !}
          end do !}
        end do !}
        write(lun,700) lv,ng_tmp,np_tmp,ns_tmp,na_tmp
        write(lun,'(2a)') '# ',label
c       hlp_frf='(9h(5h#time,,i3,8h(1h'//tab//',a),,'//
c     /                                      'i3,8h(1h'//tab//',a)))'
        hlp_frf='(9h(5h#time,,i3,8h(1h'//tab//',a)))'
        write(hlp_frm,hlp_frf) ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp
        if(npfrgrp.gt.0) then !{
          write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp),
     ,        (((hlp_ha(i,j,k),i=i0,npfrgrp),j=1,nspecies),k=1,npff_tmp)
        else !}{
          write(lun,hlp_frm) (hlp_hg(i),i=1,ng_tmp)
        end if !}
c        print *,'b2usrtrc frm:',hlp_frm  !###
      end if !}
      write(hlp_frm,'(10h(1p,e16.7,,i3,12h(1h'//tab//',e12.3)))')
     -                                ng_tmp+np_tmp*(na_tmp+1-i0)*ns_tmp
      return
c=======================================================================
 90   write(*,*) 'b2usrtrc: cannot open a file ',f
      filedata(ltrc)=.false.
      return

 92   write(*,*) 'b2usrtrc: wrong format of the file ',f
      filedata(ltrc)=.false.
      return

 94   write(*,*) 'b2usrtrc: cannot write in the file ',f
      filedata(ltrc)=.false.
      return
c=======================================================================
      entry b2usrtrs(nx,ny,ns,BoRiS)

      return
c=======================================================================
      entry b2usrtrf(nx,ny,ns,BoRiS)
c     Finalization step
      return
c=======================================================================
  700 format('#   ',5i8)
c=======================================================================
      end subroutine b2usrtrc
