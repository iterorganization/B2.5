      subroutine b2pwldw(lun,trgl,izn,iln,BoRiS)
      use b2us_data
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_diag
      use b2mod_external
      use b2mod_constants
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_wneutrals
      use eirmod_ccona
      use eirmod_cinit
#else
      use b2mod_eirbra
#endif
      use b2mod_dimensions
c
c  version : 12.01.2002 01:19
c
c=======================================================================
c*** Calculate and write down the power loading, particle fluxes, and
c*** so on, onto the segment of a plot zone or a target.
c*** Power in W/m**2, fluxes in 1/m**2/s, erosion in mm/y
c***
c*** Input:
c***  lun   : output unit number (must be opened in advance)
c***  trgl  : target (some extra data with respect to a plot zone)
c***  izn   : zone number or target number
c***  iln   : segment number within the zone
c***
c*** Output (written to the lun unit):
c***  Reference distance
c***  Total power
c***  Radiation load
c***  Total power with particles
c***  Total power with plasma (electrons and ions)
c***  Total power with neutrals
c***  Heat total
c***  Potential energy total
c***  Impinging ion flux     \
c***  Average ion energy      |  for each species
c***  Impinging atom flux     |
c***  Average atom energy    /
c***  Impinging molecule flux    \ for each molecular species
c***  Average molecule energy    /
c***  Net erosion             |   for each kind of sputtered material
c***  Gross erosion           |
c***
c*** Extra output for targets:
c***  ne
c***  Te
c***  Ti
c***  Midplane coordinate
c***
c*** Apart from this, some extra output is written to the common block
c*** for the targets:
c***  Peak power load - total
c***  Peak power with radiation
c***  Peak power with plasma
c***  Peak power with neutrals
c***  Peak heat
c***  Peak potential
c***  Radiation power at the total peak
c***  Plasma power at the total peak
c***  Neutral power at the total peak
c***  Heat at the total peak
c***  Potential at the total peak
c***  Power total
c***  Radiation total
c***  Plasma total
c***  Neutrals total
c***  Heat total
c***  Potential total
c***  Sputtered particles - gross | for each kind of sputtered material
c***  Sputtered particles - net   |
c***  Ion flux total        \
c***  Neutral flux total     |  for each species
c***  Net particle flux     /
c=======================================================================
      implicit none
#ifndef B25_EIRENE
      real (kind=R8) :: EionH2
      parameter (EionH2=15.4_R8)
      character*8 TEXTA(DEF_NATM)
#else
      integer :: use_eirene = 0
      real (kind=R8) :: pump, pumptot
      character*13 species_name
#endif
      integer lun,izn,iln
      real (kind=R8) :: BoRiS
      logical trgl,ft44

c*** Local variables
c***  ils     : the segment index in the whole list
c***  jtrg    : target number jtrg prepared
c***  izns    : plot zone izns started
c***  iznt    : plot zone of the last done target
c***  zz      : reference distance
c***  i_ersmtr: index of eroded material for the wall segments

      integer i,j,k,kk,l,m,i1,i2,ib,is,iss,iy,im,js,ils,itg,
     ,  jtrg,izns,iznt,
     ,  iCb,iCv,iFc,iFs,iFs1,iFs2,iVx1,iVx2,
     ,  i_ersmtr(nsgmx),l_ersmtr(0:nsgmx),
     ,  iatm,imol,ia,ims,IAS,shift,
     ,  sptlistrg(ntrgsx)
      integer :: nyd
      parameter (nyd=DEF_NYD)
      integer midplane_list(nromp),nmp

c*** Arrays to store the target profiles
c***  xtrg    : coordinate along the target
c***  xmp_hlp : the same projected to the midplane
c***  s       : ring area
c***  netrg   : electron density
c***  tetrg   : electron tempearture
c***  titrg   : ion tempearture
c***  wptnu   : neutral recombination
c***  wptpl   : plasma recombination
c***  whtnu   : neutral heat
c***  whtpl   : plasma heat
c***  wpttrg  : particle power to the target
c***  wpltrg  : plasma power to the target
c***  wnutrg  : neutral power to the target
c***  wnmtrg  : power with molecules to the target
c***  wnatrg  : power with atoms to the target
c***  wrdtrg  : radiation power to the target
c***  wtttrg  : total power to the target
c***  fiitrg  : impinging ion flux
c***  eiitrg  : average energy of ions
c***  faitrg  : impinging atom flux
c***  eaitrg  : average energy of atoms
c***  fmitrg  : impinging molecules flux
c***  emitrg  : average energy of molecules
c***  pmitrg  : pressure of molecules
c***  nmitrg  : density of molecules
c***  naitrg  : density of neutral atoms
c***  paitrg  : pressure of neutral atoms
c***  erntrg  : erosion [mm/y]
c***  spttrg  : sputtering [g/y]

c*** Arrays to store the peak values
c***  wpk_tt  : peak power - total
c***  wpk_rd  : peak power - radiation
c***  wpk_pl  : peak power - plasma
c***  wpk_nu  : peak power - neutrals
c***  wlc_tt  : location of total peak
c***  wlc_rd  : location of radiation peak
c***  wlc_pl  : location of plasma peak
c***  wlc_nu  : location of neutral peak
c***  wrd_pk  : radiation at total peak
c***  wpl_pk  : plasma at total peak
c***  wnu_pk  : neutrals at total peak

      real (kind=R8) :: h,p,t,u,wo,wwo,mion,
     ,  frac,psit1,psit2,psim1,psim2,fs_ds_mp(mpg%nFs),
     ,  wi(ntrgsx),wwi(ntrgsx),z_hlp(ntrgsx),
     ,  s(0:NYD-1,ntrgsx),xmp_hlp(0:NYD-1,ntrgsx),xtrg(0:NYD-1,ntrgsx),
     ,  netrg(0:NYD-1,ntrgsx), tetrg(0:NYD-1,ntrgsx),
     ,  titrg(0:NYD-1,ntrgsx),
     ,  wtttrg(0:NYD-1,ntrgsx), wptnu(0:NYD-1,ntrgsx),
     ,  wptpl(0:NYD-1,ntrgsx), wnutpr(0:NYD-1,ntrgsx),
     ,  whtnu(0:NYD-1,ntrgsx), whtpl(0:NYD-1,ntrgsx),
     ,  wpttrg(0:NYD-1,ntrgsx),
     ,  wpltrg(0:NYD-1,ntrgsx), wnutrg(0:NYD-1,ntrgsx),
     ,  wrdtrg(0:NYD-1,ntrgsx,0:DEF_NATM),
     ,  fiitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  eiitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  faitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  eaitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  naitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  paitrg(0:NYD-1,ntrgsx,DEF_NATM),
     ,  fmitrg(0:NYD-1,ntrgsx,DEF_NMOL),
     ,  emitrg(0:NYD-1,ntrgsx,DEF_NMOL),
     ,  pmitrg(0:NYD-1,ntrgsx,DEF_NMOL),
     ,  nmitrg(0:NYD-1,ntrgsx,DEF_NMOL),
     ,  erntrg(0:NYD-1,ntrgsx), spttrg(0:NYD-1,ntrgsx),
     ,  sjs(0:NYD-1,ntrgsx), sjs_t(ntrgsx),
     ,  sjs_w(0:ntrgsx), sptmass, sptflx(DEF_NLIM+DEF_NSTS)
c      real (kind=R8) :: zz

c*** Recap of fluxes from Eirene
c***    wldnek  :   net kinetic energy deposited by neutrals. Evaluated as
c***                   EOTAT+EOTML+EOTIO-ERFAAT-ERFMAT-ERFIAT-ERFAML-ERFMML-ERFIML
c***    wldnep  :   potential energy released by neutrals. Evaluated as
c***                   diss_pot_H2*(PRFAML+PRFPML)
c***    wldna   :   flux of atoms impinging onto the surface. Equal to POTAT
c***    ewlda   :   their average energy. Equal to EOTAT/POTAT
c***    wldnm   :   flux of molecules impinging onto the surface. Equal to POTML
c***    ewldm   :   their average energy. Equal to EOTML/POTML
c***    wldra   :   flux of reflected atoms. Equal to PRFAAT+PRFMAT+PRFIAT
c***    wldrm   :   flux of reflected molecules. Equal to PRFAML+PRFMML+PRFIML
c***    wldpp   :   flux of plasma ions impinging onto the surface. Equal to POTPL
c***    wldpa   :   flux of emitted atoms from BI. Equal to PRFPAT
c***    wldpm   :   flux of emitted molecules from BI. Equal to PRFPML
c***    wldpeb  :   power carried away by these atoms and molecules. Equal to
c***                  ERFPAT+ERFPML+ERFPIO
c***    wldspt  :   flux of sputtered wall material
c***    wldspta :   sputtered flux by emitted atom
c***    wldsptm :   sputtered flux by emitted molecule
c***    wlarea  :   areas of the surface segments from Eirene
c***    wlabsrp :   absorption at the surfaces (1-RECYCT from Eirene)
c***    wlpump  :   pumped flux for each pumped particle
c***    ewldea_res  :   emitted energy flux of atoms. Equal to
c***                      ERFAAT+ERFMAT+ERFIAT+ERFPHAT
c***    ewldem_res  :   emitted energy flux due of molecules. Equal to
c***                      ERFAML+ERFMML+ERFIML+ERFPHML
c***    ewldrp_res  :   kinetic energy of reflected atoms+molecules originated from ions.
c***                      ERFPAT+ERFPML
c***    ewldmr_res  :   energy due to recombination of atoms and atomic ions into molecules.
c***                     (PRFAML+PRFPML)*H_EDISS

c*** Output data headers
c***  m_hdpr    : max. length of a quantity name
c***  m_hdln    : max. length of the header line
c***  n_hdpr_c  : number of common quantities
c***  n_hdpr_t  : number of target quantities
c***  n_hdpr_a  : number of quantities by atoms
c***  n_hdpr_ml : number of quantities by molecules
c***  n_hdpr_s  : number of quantities by species
c***  n_hdpr_sp : number of quantities for sputtering
c***  n_hdpr_er : number of quantities for erosion
c***  n_sfpr_s  : number of active species
c***  n_sfpr_sp : number of different sputtered materials
c***  n_sfpr_er : number of different eroded materials
c***  hdpr_c    : headers for common quantities
c***  hdpr_t    : headers for target quantities
c***  hdpr_a    : headers for the quantities by atoms
c***  hdpr_s    : headers for the quantities by species
c***  hdpr_sp   : headers for the sputtered quantities
c***  hdpr_er   : headers for the eroded quantities
c***  hdpr_ml   : headers for the quantities by molecules
c***  sfpr_s    : list of active species
c***  sfpr_a    : list of active atoms
c***  sfpr_m    : list of active materials

      integer n_hdpr_c,n_hdpr_t,n_hdpr_s,n_hdpr_a,n_hdpr_sp,n_sfpr_s,
     ,  m_hdpr,m_hdln,n_hdpr_ml,n_hdpr_er,n_srfs_m,nclm_na,nclm_nm
      parameter
     , (n_hdpr_c=10, n_hdpr_t=4, n_hdpr_s=2, n_hdpr_a=4, n_hdpr_sp=2,
     ,  n_hdpr_er=2, m_hdpr=5, n_hdpr_ml=4, n_srfs_m=8,
     ,  m_hdln=(n_hdpr_c+n_hdpr_t)*(m_hdpr+1)+
     +         (n_hdpr_s+n_hdpr_a+n_hdpr_sp)*DEF_NATM*(m_hdpr+4)+
     +          n_hdpr_ml*DEF_NMOL*(m_hdpr+4)+
     *          n_hdpr_er*n_srfs_m*(m_hdpr+4))
      integer n_line_x, n_sfpr_sp, n_sfpr_er
      integer do_spt_wll(DEF_NLIM+DEF_NSTS)
      integer i_sfpr_sp(DEF_NATM+DEF_NMOL)
      real (kind=R8), allocatable :: w_line(:)
      real (kind=R8) :: nsecyr
      parameter (nsecyr=24.0_R8*3600.0_R8*365.2425_R8)
      logical :: first,do_trg,match_found,do_spt
      character*(m_hdpr) tab*1,hdpr_c(n_hdpr_c),hdpr_t(n_hdpr_t),
     ,  hdpr_s(n_hdpr_s),hdpr_a(n_hdpr_a),hdpr_sp(n_hdpr_sp),
     ,  sfpr_s(DEF_NATM)*2,sfpr_a(DEF_NATM)*2,
     ,  sfpr_m(nsptmx)*(lnwlmtrl),hdr_line*(m_hdln),frmt_hlp*24,
     ,  hdpr_ml(n_hdpr_ml), hdpr_er(n_hdpr_er),
     ,  sfpr_sp(DEF_NATM+DEF_NMOL)*3
      character*16 ssmol,ssion
      data hdpr_c /'x', 'Wtot','Wrad','Wpart','Wpls',
     2  'Wneu','Wheat','Wpot', 'Whtpl', 'Wptpl' /
      data hdpr_t /'ne', 'Te','Ti','xMP'/
      data hdpr_s /'flxi','Eavi'/
      data hdpr_a /'flxa','Eava','na','pa'/
      data hdpr_ml /'flxm','Eavm','nm','pm'/
      data hdpr_sp /'sptne','sptng'/
      data hdpr_er /'ersne','ersng'/
      data jtrg, izns, iznt, n_sfpr_s /4*0/
      data first /.true./
      data sfpr_s, sfpr_a, sfpr_m /
     /     DEF_NATM*' ', DEF_NATM*' ', nsptmx*' ' /

c*** Data for erosion calculations
c***  n_srfs_m  : number of possible wall materials
c***  srfs_m    : list of possible wall materials
c***  srfs_x    : corresponding species (say, Fe for SS)
c***  srfs_d    : specific masses of possible wall materials [g/cm**3]
c***  srfs_a    : atomic weight of possible wall materials
c***  wllersf   : ratio of erosion to sputtering for different materials
      integer map_mtr(n_srfs_m)
      real (kind=R8) :: srfs_d(n_srfs_m),srfs_a(n_srfs_m),
     .                  wllersf(nsptmx)
      character srfs_m(n_srfs_m)*4,srfs_x(n_srfs_m)*2
      character*256 filename
      logical streql
      external streql
      data srfs_m /'Be','C','CFC','Fe','SS','ss','Mo' ,'W'/
      data srfs_x /'Be','C','C'  ,'Fe','Fe','Fe','Mo' ,'W'/
      data srfs_a /9.0_R8,12.0_R8,12.0_R8,56.0_R8,56.0_R8,56.0_R8,
     .            96.0_R8,184.0_R8/
      data srfs_d /1.8_R8, 2.0_R8, 2.0_R8, 7.8_R8, 7.8_R8, 7.8_R8,
     .            10.2_R8,19.3_R8/
      save
c=======================================================================

c      print *,'Starting b2pwldw...',lun,trgl,izn,iln  !###

c*** At present, only a crude estimate for the neutrals at the targets:
c*** power flux proportional to n*T**(3/2),
c*** particle flux proportional to n*T**(1/2),
c*** sputtering proportional to the D particle flux (chem. dominant),
c*** re-deposition proportional to the C particle flux.
c*** Here n and T are the density and temperature of the hydrogenic
c*** atoms near the target.
c*** The proportionality factor is adjusted to match the integrals
c*** which are calculated in Eirene.
c-----------------------------------------------------------------------

      if(first) then !{
c        print *,'Starting b2pwldw first time',lun,trgl,izn,iln  !###
        tab=char(9)
        filename=fort_lc//'44'
        call find_file(filename,ft44)
        do_trg=.true.

        n_sfpr_s=1
        if(nspecies.gt.1) n_sfpr_s=nspecies
        do is=1,n_sfpr_s !{
          sfpr_s(is)=adjustl(species_list(is))
        end do !}
#ifndef B25_EIRENE
        nclm_na=nspecies
        nclm_nm=0
        do ia=1,nspecies
          TEXTA(ia)=species_list(ia)
        end do
#else
        call ipgeti ('b2mndr_eirene', use_eirene)
        if (use_eirene.eq.0) then
          nclm_na=nspecies
          nclm_nm=0
          if (.not.allocated(texta)) allocate(texta(nspecies))
          do ia=1,nspecies
            TEXTA(ia)=species_list(ia)
          end do
        else
          nclm_na=natmi
          nclm_nm=nmoli
        endif
#endif
        do ia=1,nclm_na
          sfpr_a(ia)=TEXTA(ia)(1:2)
        end do
c        print *,'First call. format:',frmt_hlp,':'                  !###
c        print *,'n_sfpr_s=',n_sfpr_s                                !###
c        if(n_sfpr_s.gt.0) print *,(' ',sfpr_s(i),i=1,n_sfpr_s)      !###

c*** data for sputtering and erosion

c        print *,'Input wall materials'                              !###
c        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=1,lppp)            !###
c        print *,'Input target materials'                            !###
c        print '(1h:,20(a,1h:))',(mtrtrgsg(itrg(i)),i=1,ntrgts)      !###
        k=lppp
        do j=1,ntrgts !{
          if (itrg(j).eq.0) cycle
          do i=k+1,k+mpg%divFcP(itrg(j),2) !{
c carlis insert wall material for each target element
            mtrwldsg(i)=mtrtrgsg(itrg(j))
          end do !}
          k=k+mpg%divFcP(itrg(j),2)
        end do !}
        kk=k
c        print *,'k,lppp,ntrgts=',k,lppp,ntrgts                      !###
c        print *,'Target segment materials'                          !###
c        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=lppp+1,k)          !###

        n_sfpr_er=0
        do k=1,kk !{
          l=0
          if(.not.lskip(k).and.
     &       .not.streql(mtrwldsg(k),' ').and.
     &       .not.streql(mtrwldsg(k),'dummy')) then !{
            do i=1,n_sfpr_er !{
c carlis l is index of wall material referred to list of ACTIVE materials
              if(streql(sfpr_m(i),mtrwldsg(k))) l=i
            end do !}
            if(l.eq.0) then !{
              if(n_sfpr_er.lt.nsptmx) then !{
                n_sfpr_er=n_sfpr_er+1
                l=n_sfpr_er
c carlis added if not in the list of ACTIVE materials
                sfpr_m(l)=mtrwldsg(k)
c              print *,'Added :',sfpr_m(l),':',k,l                   !###
              else !}{
              end if !}
            end if !}
          end if !}
c carlis list of indices of wall materials referred to list of ACTIVE materials
          i_ersmtr(k)=l
        end do !}
c        print *,'Input materials: ',sfpr_m                        !###
c        if(n_sfpr_m.gt.0) print *,':',(sfpr_m(i),':',i=1,n_sfpr_m)  !###
        l_ersmtr(0)=0
        map_mtr(:)=0
        m=0
        do k=1,kk !{
          l=0
          if(i_ersmtr(k).ne.0) then !{
            do i=1,n_srfs_m !{
c carlis l is index of wall material in ACTIVE list referred to list of POSSIBLE materials
              if(streql(srfs_m(i),sfpr_m(i_ersmtr(k)))) then !{
                l=i
c carlis in list of ACTIVE materials we put the 'real' material e.g. CFC becomes C
c        the list is not ordered as the list of possible materials
                if(map_mtr(i_ersmtr(k)).eq.0) then !{
                  m=m+1
                  sfpr_m(m)=srfs_x(l)
c carlis     [s/y*m2*mm]   =[s/y]*[m2/cm2*mm/cm*]*[g/mol]/[1/mol]/[g/cm3]
                  wllersf(m)=nsecyr*1.e-3*srfs_a(l)/(avogr*srfs_d(l))
                  map_mtr(i_ersmtr(k))=m
                else !}{
                  m=map_mtr(i_ersmtr(k))
                end if !}
              end if !}
            end do !}
            if(l.eq.0) write(*,*)
     .       'b2pwldw: unknown wall material ',sfpr_m(i_ersmtr(k))
c carlis list of indices that tells the new ordering of wall materials in active list wrt the ordering of i_ersmtr
            l_ersmtr(k)=m
          end if !}
        end do !}
        spttrg(:,:)=0.0_R8
        erntrg(:,:)=0.0_R8
        n_sfpr_er=m
        do i=1,kk !{
c carlis substitute indices of wall materials with the new ordering of ACTIVE list
          i_ersmtr(i)=l_ersmtr(i)
        end do !}
        write(*,*) 'Number of wall materials: ',n_sfpr_er
        if(n_sfpr_er.gt.0) then !{
          write(*,*) 'Wall materials : ',(sfpr_m(i),i=1,n_sfpr_er)
          write(*,'(a,1p,10e10.2)') ' Erosion factors: ',
     ,                                        (wllersf(i),i=1,n_sfpr_er)
        end if !}
c        print *,'i_sptmtr'                                    !###
c        print '(1x,50i2)',(i_sptmtr(i),i=1,lppp+kk)           !###

        do_spt_wll(:)=0
        i_sfpr_sp(:)=0
        sptflx(:)=0.0_R8
        i=0
#ifdef B25_EIRENE
c carlis evaluate number of sputtered atoms and molecules
        do ia=1,natmi
          match_found=.false.
          j=1
          do j=1,nlimps
            if (wldspta(j,ia,0).gt.0.0_R8) then
              sptflx(j)=wldspta(j,ia,0)
              do_spt_wll(j)=ia
              match_found=.true.
            end if
          end do
          if (match_found) then
            i=i+1
            i_sfpr_sp(i)=ia
            sfpr_sp(i)=trim(TEXTA(ia))
          end if
        end do
        do imol=1,nmoli
          match_found=.false.
          j=1
          do j=1,nlimps
            if (wldsptm(j,imol,0).gt.0.0_R8) then
              sptflx(j)=wldsptm(j,imol,0)
              do_spt_wll(j)=imol+natmi
              match_found=.true.
            end if
          end do
          if (match_found) then
            i=i+1
            i_sfpr_sp(i)=imol+natmi
            sfpr_sp(i)=trim(adjustl(textmn(imol-1)))
          end if
        end do
#endif
        n_sfpr_sp=i
        write(*,*) 'Number of sputtered species: ',n_sfpr_sp
        if(n_sfpr_sp.gt.0) then !{
          write(*,*) 'Sputtered species : ',(sfpr_sp(i),i=1,n_sfpr_sp)
        end if !}
        n_line_x=n_hdpr_c+n_hdpr_t+
     ,           n_hdpr_s*n_sfpr_s+nclm_na*n_hdpr_a+
     ,           n_hdpr_ml*nclm_nm+n_hdpr_sp*n_sfpr_sp+
     ,           n_hdpr_er*n_sfpr_er
        write(frmt_hlp,'(a7,i5,a11)') '(1p,1x,',n_line_x,'(e12.5,a1))'
c carlis calculation of total pumped fluxes
#ifdef B25_EIRENE
        print *,'Pumped fluxes'
        k=0
        pumptot=0.0_R8
        do ia=1,natmi
          pump=0.0_R8
          do j=1,nlimps
            if (wlabsrp(k+ia,j).gt.0.0) then
              pump=pump+wlpump(k+ia,j)
            end if
          end do
          pumptot=pumptot+pump
          print *,'Pumped ',trim(TEXTA(ia)),': ',pump,' (1/s)'
        end do

        k=k+natmi
        do ia=1,nmoli
          pump=0.0_R8
          do j=1,nlimps
            if (wlabsrp(k+ia,j).gt.0.0) then
              pump=pump+wlpump(k+ia,j)
            end if
          end do
          pumptot=pumptot+pump
          print *,'Pumped ',trim(adjustl(textmn(ia-1))),
     ,': ',pump,' (1/s)'
        end do

        k=k+nmoli
        do ia=1,nioni
          pump=0.0_R8
          do j=1,nlimps
            if (wlabsrp(k+ia,j).gt.0.0) then
              pump=pump+wlpump(k+ia,j)
            end if
          end do
          pumptot=pumptot+pump
          print *,'Pumped ',trim(adjustl(textin(ia-1))),
     ,': ',pump,' (1/s)'
        end do

        k=k+nioni
        ia=0
        do is = 0, state%ns-1
          if (.not.is_neutral(is)) then
            ia=ia+1
            pump=0.0_R8
            do j=1,nlimps
              if (wlabsrp(k+ia,j).gt.0.0) then
                pump=pump+wlpump(k+ia,j)
              end if
            end do
            pumptot=pumptot+pump
            call species(is,species_name,.false.)
            print *,'Pumped ',trim(adjustl(species_name)),
     ,              ': ',pump,' (1/s)'
          end if
        end do
        print *,'Total pumped flux: ',pumptot,' (1/s)'
#endif

        first=.false.
      end if !} end first call
      ils=iabs(lwldlis(iln,izn))
c      print *,'lwldlis,ils,lppp=',lwldlis(iln,izn),ils,lppp       !###

c*** Check whether a target segment is requested

      if(ils.le.lppp) then !{
c*** Wall segment. Reset the target pointer
        itg=0
      else !}{
        if(trgl) then !{
          itg=izn
        else !}{
          k=ils-lppp
          itg=0
          do j=1,ntrgts !{
            if (itg.gt.0) cycle
            if (itrg(j).eq.0) cycle
            if (k.le.mpg%divFcP(itrg(j),2)) then
              itg=itrg(j)
            else
              k=k-mpg%divFcP(itrg(j),2)
            end if
          end do !}
        end if !}
c        print *,'itg,iznt,jtrg=',itg,iznt,jtrg  !###
c*** Check whether this target has been prepared
        if(itg.ne.jtrg) then !{

          if(jtrg.eq.0 .and. do_trg) then !{

c*** This is the first request for a target segment.
c*** Prepare the necessary data for all the targets - just to make
c*** the modification from the old version easier.
c*** Here we need to calculate all totals and peaks, and also to
c*** evaluate the neutral fluxes which are not directly available
c*** from Eirene.
c*** All the extra output for targets is also prepared here.

c*** Determine topology and calculate radiation load onto the targets
c*** Indeed, b2ptrdl works only once

            call b2ptrdl(wrdtrg)

c            print *,'After b2ptrdl: ntrgts = ',ntrgts   !###

c*** Target geometry

            do j=1,ntrgts !{
              z_hlp(j)=0.0_R8
              do iy=1,lseppos(itrg(j))-1 !{
                iFc = mpg%divFc(mpg%divFcP(itrg(j),1) + iy - 1)
                z_hlp(j)=z_hlp(j) - geo%fcHt(iFc)
              end do !}
              do iy=1,mpg%divFcP(itrg(j),2) !{
                iFc = mpg%divFc(mpg%divFcP(itrg(j),1) + iy - 1)
                s(iy,itrg(j))=geo%fcS(iFc)
                z_hlp(j)=z_hlp(j)+0.5_R8*geo%fcHt(iFc)
                xtrg(iy,itrg(j))=z_hlp(j)
                z_hlp(j)=z_hlp(j)+0.5_R8*geo%fcHt(iFc)
              end do !}
            end do !}

            if(trgl) then !{

c*** Midplane coordinates
              do j=1,2*ntrgts
                if (j.le.ntrgts/2 .and. nimp.gt.0) then
                  midplane_list = imp
                  fs_ds_mp = fs_ds_imp
                  im = icsepimp
                  nmp = nimp
                  l = -1
                else
                  midplane_list = omp
                  fs_ds_mp = fs_ds_omp
                  im = icsepomp
                  nmp = nomp
                  l = 1
                end if
                u=0.0_R8
                if (itrg(j).ne.0) then !{
                  do i=1,lseppos(itrg(j))-1 !{
                    iFc = mpg%divFc(mpg%divFcP(itrg(j),1) + i - 1)
                    h=0.5_R8*geo%fcHt(iFc)
                    u=u+h
                    xmp_hlp(i,itrg(j))=u
                    u=u+h
                  end do !}
                  do i=1,lseppos(itrg(j))-1 !{
                    xmp_hlp(i,itrg(j))=xmp_hlp(i,itrg(j))-u
                  end do !}
                  u=0.0_R8
                  do i=lseppos(itrg(j)),mpg%divFcP(itrg(j),2) !{
! May be a partial mapping only!
! First, find to which flux tube the target surface belongs
! Second, determine how much of the flux tube it intersects
! Third, calculate the distance according to the midplane distances
                    iFc = mpg%divFc(mpg%divFcP(itrg(j),1) + i - 1)
                    iVx1 = mpg%fcVx(iFc,1)
                    iVx2 = mpg%fcVx(iFc,2)
                    iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                    iFs1 = mpg%vxFs(iVx1)
                    iFs2 = mpg%vxFs(iVx2)
                    if (iFs1.eq.0.or.iFs2.eq.0) then
                      Psit1 = minval(
     .                 geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1):
     .                                     mpg%cvVxP(iCv,1)+
     .                                     mpg%cvVxP(iCv,2)-1)))
                      Psit2 = maxval(
     .                 geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1):
     .                                     mpg%cvVxP(iCv,1)+
     .                                     mpg%cvVxP(iCv,2)-1)))
                      Psim1 = minval(
     .         geo%vxFpsi(mpg%cvVx(mpg%cvVxP(midplane_list(im),1):
     .                             mpg%cvVxP(midplane_list(im),1)+
     .                             mpg%cvVxP(midplane_list(im),2)-1)))
                      Psim2 = maxval(
     .         geo%vxFpsi(mpg%cvVx(mpg%cvVxP(midplane_list(im),1):
     .                             mpg%cvVxP(midplane_list(im),1)+
     .                             mpg%cvVxP(midplane_list(im),2)-1)))
                      do while (im.gt.0.and.im.lt.nmp.and.
     .                 (Psit1.lt.Psim1.or.Psit1.gt.Psim2.or.
     .                  Psit2.lt.Psim1.or.Psit2.gt.Psim2))
                        im = im+l
                        Psim1 = minval(
     .           geo%vxFpsi(mpg%cvVx(mpg%cvVxP(midplane_list(im),1):
     .                               mpg%cvVxP(midplane_list(im),1)+
     .                               mpg%cvVxP(midplane_list(im),2)-1)))
                        Psim2 = maxval(
     .           geo%vxFpsi(mpg%cvVx(mpg%cvVxP(midplane_list(im),1):
     .                               mpg%cvVxP(midplane_list(im),1)+
     .                               mpg%cvVxP(midplane_list(im),2)-1)))
                      end do
                      call xertst (.not.(im.le.0.or.im.gt.nmp),
     .                 'Stepped outside midplane range of cells!')
                      iFs1 = 0
                      iFs2 = 0
                      do k = mpg%cvVxP(midplane_list(im),1),
     .                       mpg%cvVxP(midplane_list(im),1) +
     .                       mpg%cvVxP(midplane_list(im),2) - 1
                        iFs = mpg%vxFs(k)
                        if (iFs.eq.0) cycle
                        if (iFs1.eq.0) then
                          iFs1 = iFs
                        else if (iFs2.eq.0) then
                          iFs2 = iFs
                        end if
                      end do
                      frac = (Psit2-Psit1)/(Psim2-Psim1)
                    else
                      frac = 1.0_R8
                    end if
                    h=0.5_R8*frac*abs(fs_ds_mp(iFs2)-fs_ds_mp(iFs1))
                    u=u+h
                    xmp_hlp(i,itrg(j))=u
                    u=u+h
                  end do !}
                end if !}
              end do !}
            end if !}

c*** Power loading:

c*** Plasma and neutrals

c            print *,'Starting target power. ntrgts=',ntrgts   !###
            do j=1,2*max(1,mpg%nncut) !{ Do not use ntrgts here to allow for sign flips
              if (itrg(j).gt.0) then !{
                do iy=1,mpg%divFcP(itrg(j),2) !{
                  iFc = mpg%divFc(mpg%divFcP(itrg(j),1) + iy - 1)
                  iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                  wo=0.0_R8
                  wnutrg(iy,itrg(j))=0.0_R8
                  whtpl(iy,itrg(j))=
     .             (state%dv%fhe(iFc,0)+state%dv%fhe(iFc,1)+
     .              state%dv%fhi(iFc,0)+state%dv%fhi(iFc,1)+
     .              state_ext%fhi(iFc,0)+state_ext%fhi(iFc,1))
                  do is=0,state%ns-1 !{
                    if (.not.is_neutral(is).or.ft44) then
                      wo=wo+(state%dv%fhp(iFc,0,is)+
     .                       state%dv%fhp(iFc,1,is))
                    else
                      if (nint(zn(is)).eq.1) then
                        wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     *                   (1.0_R8-BoRiS)*
     *          (state%dv%fhm(iFc,0,is)+state%dv%fhm(iFc,1,is)+
     .          (state%dv%fna(iFc,0,is)+state%dv%fna(iFc,1,is))*
     .           state%pl%tn(iCv))
                      else
                        wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     *                   (1.0_R8-BoRiS)*
     *          (state%dv%fhm(iFc,0,is)+state%dv%fhm(iFc,1,is)+
     .          (state%dv%fna(iFc,0,is)+state%dv%fna(iFc,1,is))*
     .           state%pl%ti(iCv))
                      end if
                    end if
                    if (nint(zn(is)).eq.1.and.is_neutral(is)) then
                      whtpl(iy,itrg(j))=
     .                 whtpl(iy,itrg(j))+(1.0_R8-BoRiS)*
     *           (state%dv%fhm(iFc,0,is)+state%dv%fhm(iFc,1,is)+
     .           (state%dv%fna(iFc,0,is)+state%dv%fna(iFc,1,is))*
     .            state%pl%tn(iCv))
                    else
                      whtpl(iy,itrg(j))=
     .                 whtpl(iy,itrg(j))+(1.0_R8-BoRiS)*
     *           (state%dv%fhm(iFc,0,is)+state%dv%fhm(iFc,1,is)+
     .           (state%dv%fna(iFc,0,is)+state%dv%fna(iFc,1,is))*
     .            state%pl%ti(iCv))
                    end if
c carlis 31jan2017  added ti/tn as suggested by sw in b2blne.
                  end do !}
                  do is=0,state_ext%ns-1 !{
                    if (.not.state_ext%is_neutral(is).or.ft44) then
                      wo=wo+
     .                 (state_ext%fa(iFc,0,is)+state_ext%fa(iFc,1,is))*
     .                  state_ext%pt(iCv,is)*ev
                    else
                      wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .                (state_ext%fa(iFc,0,is)+state_ext%fa(iFc,1,is))
     .               *(1.0_R8-BoRiS)*(0.5_R8*state_ext%am(is)*mp*
     .                 state_ext%ua(iCv,is)**2+state_ext%ta(iCv,is))
                    end if
                    whtpl(iy,itrg(j))=whtpl(iy,itrg(j))+
     .              (state_ext%fa(iFc,0,is)+state_ext%fa(iFc,1,is))
     .             *(1.0_R8-BoRiS)*(state_ext%ta(iCv,is)+0.5_R8*
     .               state_ext%am(is)*mp*state_ext%ua(iCv,is)**2)
c carlis 31jan2017  added ta as suggested by sw in b2blne.
                  end do !}
                  whtpl(iy,itrg(j))=whtpl(iy,itrg(j))+
     .              (state%dv%fne(iFc,0)+state%dv%fne(iFc,1))
     .             *(1.0_R8-BoRiS)*state%pl%te(iCv)
c carlis 31jan2017  added electron contribution as suggested by sw in b2blne.
                  u = mpg%divFcOr(mpg%divFcP(itrg(j),1)+iy-1)
                  wptpl(iy,itrg(j))=wo/s(iy,itrg(j))*u
                  if (.not.ft44)
     >             wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))*u/s(iy,itrg(j))
                  whtpl(iy,itrg(j))=whtpl(iy,itrg(j))/s(iy,itrg(j))*u
                  wpltrg(iy,itrg(j))=wptpl(iy,itrg(j))+whtpl(iy,itrg(j))
                  netrg(iy,itrg(j))=state%dv%ne(iCv)
                  tetrg(iy,itrg(j))=state%pl%te(iCv)/ev
                  titrg(iy,itrg(j))=state%pl%ti(iCv)/ev
                  if (ft44) then
                    ias = b2_fnmti(iFc)
                    call xertst(ias.gt.0,
     .               'Eirene boundary face not found!')
                    do is=1,natmi
c                   carlis Using values directly calculated from Eirene
                      wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))
     .                 +(ewlda_res(is,IAS)-ewldea_res(is,IAS))*
     .                  sarea_res(IAS)
                    end do
                    do is=1,nmoli
c                   carlis Using values directly calculated from Eirene
                      wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))
     .                 +(ewldm_res(is,IAS)-ewldem_res(is,IAS)
     .                                    +ewldmr_res(is,IAS))*
     .                  sarea_res(IAS)
                    enddo
                    wnutpr(iy,itrg(j))=ewldrp_res(IAS)
                    whtpl(iy,itrg(j))=
     >               whtpl(iy,itrg(j))-wnutpr(iy,itrg(j))
                    wpltrg(iy,itrg(j))=
     >               wpltrg(iy,itrg(j))-wnutpr(iy,itrg(j))
                    do is=1,nioni
                      mion=0.0_R8
                      imol=1
                      ssion=textin(is-1)
                      ssion=trim(adjustl(ssion))
                      match_found=.false.
                      do while (.not.match_found.and.imol.le.nmoli)
                        ssmol=textmn(imol-1)
                        ssmol=trim(adjustl(ssmol))//'+'
                        if (streql(ssion,ssmol)) then
                          match_found=.true.
                          do ims = 1,natmi
                            mion=mion+mlcmp(ims,imol)*am(eb2atcr(ims))
                          enddo
                        endif
                        imol=imol+1
                      enddo
                      t=tib2(iCv,is,1)/ev
                      if (mion.gt.0.0_R8)
     >                 wnutrg(iy,itrg(j))=wnutrg(iy,itrg(j))+
     .                  s(iy,itrg(j))*dib2(iCv,is,1)*
     *                  sqrt(max(t/mion,0.0_R8))*0.25_R8*1.0e4_R8*
     .                 (2.0_R8*tib2(iCv,is,1)+EionH2*ev)
                    enddo
                  end if
                end do !} loop over boundary faces
              end if !}
            end do !} loop over targets

c*** Re-scale the neutral target loading and calculate the sums

            do j=1,ntrgts !{
              wi(itrg(j))=0.0_R8
              do i=1,mpg%divFcP(itrg(j),2) !{
                wi(itrg(j))=wi(itrg(j))+wnutrg(i,itrg(j))
              end do !}
            end do !}
c            print *,'ntrgts, wi: ',ntrgts,(wi(itrg(i)),i=1,ntrgts)  !###
            do j=1,ntrgts !{
              if (ft44) then
                if(count(ktns(:)==itrg(j)).gt.0) then !{
                  wo=0.0_R8
                  do iss=1,ntns !{
                    if (ktns(iss).ne.itrg(j)) cycle
                    js=ltns(iss)+nlim
                    wo=wo+wldnek(js,0)+wldnep(js,0)
c                    print *,'j,js,wo = ',j,js,wo    !###
                  end do !}
                  if(wi(itrg(j)).gt.0.0_R8) then !{
                    wi(itrg(j))=wo/wi(itrg(j))
                  else !}{
                    wi(itrg(j))=0.0_R8
                  end if !}
                else !}{
                  write(*,*) 'b2pwldw: no Eirene surface for target ',j,
     ,              ' - no neutral load onto the target calculated'
                  write(*,*) '  ==> check specification of ltns in ',
     ,              'b2.neutrals.parameters'
                  wi(itrg(j))=0.0_R8
                end if !}
              endif
              do iy=1,mpg%divFcP(itrg(j),2) !{
                if (ft44) then
                  wnutrg(iy,itrg(j))=
     >             wnutrg(iy,itrg(j))*wi(itrg(j))/s(iy,itrg(j))
                  if(abs(wo).gt.0.0_R8) then !{
                    do iss=1,ntns !{
                      if (ktns(iss).ne.itrg(j)) cycle
                      js=ltns(iss)+nlim
                      wptnu(iy,itrg(j))=wnutrg(iy,itrg(j))*
     >                 wldnep(js,0)/(wldnek(js,0)+wldnep(js,0))
                    end do !}
                    whtnu(iy,itrg(j))=
     >               wnutrg(iy,itrg(j))-wptnu(iy,itrg(j))
                  else !}{
                    wptnu(iy,itrg(j))=0.0_R8
                    whtnu(iy,itrg(j))=0.0_R8
                  end if !}
                end if
                wpttrg(iy,itrg(j))=wnutrg(iy,itrg(j))+wpltrg(iy,itrg(j))
                wtttrg(iy,itrg(j))=
     >           wpttrg(iy,itrg(j))+wrdtrg(iy,itrg(j),0)
              end do !}
            end do !}
            if (ntrgts.gt.0) write(*,'(1x,a,1p,16e12.3)')
     ,       'Scaling factors for Wn : ', (wi(itrg(i)),i=1,ntrgts)

c*** Particle fluxes and energies
            i2=0
            do is=1,n_sfpr_s !{
              i1=eb2spcr(is)+1
              i2=eb2spcr(is)+nfluids(is)
              do j=1,2*max(1,mpg%nncut) !{ Do not use ntrgts here to allow for sign flips
                if (itrg(j).gt.0) then !{
c!###{
c                  print '(a,7i6)','is,j,i1,i2,k,ib=',
c     ,                       is,itrg(j),i1,i2,k,ib
c                  zz=0.0_R8
c                  print '(5x,7(2x,a8))','phi','fti','wo','wwo',
c     ,             'fiitrg','eiitrg','s'
c!###}
                  do iy=1,mpg%divFcP(itrg(j),2) !{
                    iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                    iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                    iCb = max(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                    ib = mpg%fcBc(iFc,2)
                    wo=0.0_R8
                    wwo=0.0_R8
                    p=(state%pl%po(iCb)
     .     -get_phi_applied(switch%b2stbc_boundary_namelist,mpg,iFc))*qe
                    h=state%pl%ti(iCb)*
     .               (enipar(ib,1)+2.0_R8-enipar(ib,2))
                    do i=i1,i2 !{
                      wo=wo+
     .                 (state%dv%fna(iFc,0,i)+state%dv%fna(iFc,1,i))
                      wwo=wwo+(h+p*state%rt%rza(iCv,i))*
     .                 (state%dv%fna(iFc,0,i)+state%dv%fna(iFc,1,i))
                    end do !}
c                    zz=zz+wwo*u                 !###
                    fiitrg(iy,itrg(j),is)=u*wo
                    if(abs(wo).gt.0.0_R8) then !{
                      eiitrg(iy,itrg(j),is)=wwo/wo
                    else !}{
                      eiitrg(iy,itrg(j),is)=0.0_R8
                    end if !}
c!###{
c                    print '(1x,i4,1p,8e10.2)',iy,p/ev,h/ev,wo,wwo,
c     ,               fiitrg(iy,itrg(j),is),eiitrg(iy,itrg(j),is)/ev,
c     ,               s(iy,itrg(j))
c!###}
                  end do !} loop over boundary faces
c                print *,'j,is,zz=',itrg(j),is,zz  !###
                end if !}
              end do !} loop over targets
c              write(*,*) 'Species',is
c              write(*,*) 'Flux scaling'
c              write(*,'(1p,12e10.2)') (wi(itrg(i)),i=1,ntrgts)
            end do !} loop over active species

c*** Atomic fluxes and energies
            do iatm=1,nclm_na !{
              do j=1,2*max(1,mpg%nncut) !{ Do not use ntrgts here to allow for sign flips
                if (itrg(j).gt.0) then !{
                  if (ft44) then
c!###{
c                    print '(a,7i6)','is,j,i1,i2,k,ib=',
c     ,                         is,itrg(j),i1,i2,k,ib
c                    print '(5x,4(2x,a8))',
c     ,                         'faitrg','eaitrg','dab2','tab2'
c!###}
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                      iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                      naitrg(iy,itrg(j),iatm)=dab2(iCv,iatm,1)
                      paitrg(iy,itrg(j),iatm)=dab2(iCv,iatm,1)*
     *                                        tab2(iCv,iatm,1)
c carlis using values directly from Eirene
                      ias = b2_fnmti(iFc)
                      call xertst(ias.gt.0,
     .                 'Eirene boundary face not found!')
                      faitrg(iy,itrg(j),iatm)=wldna_res(iatm,IAS)*
     .                       sarea_res(IAS) !units in 1/m2/s*m2
                      eaitrg(iy,itrg(j),iatm)=ewlda_res(iatm,IAS)*
     .                       sarea_res(IAS) !units in W/m2*m2
c!###{
c                    print '(1x,i4,1p,4e10.2)',iy,
c     ,               faitrg(iy,itrg(j),iatm),eaitrg(iy,itrg(j),iatm),
c     ,               dab2(iCv,iatm,1),tab2(iCv,iatm,1)/ev
c!###}
                    end do !} loop over boundary faces
                    wo=0.0_R8
                    wwo=0.0_R8
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      wo=wo+faitrg(iy,itrg(j),iatm) ! units 1/s
                      wwo=wwo+eaitrg(iy,itrg(j),iatm) ! units W
                    end do !}
                    wi(itrg(j))=0.0_R8
                    wwi(itrg(j))=0.0_R8
                    do iss=1,ntns !{
                      if (ktns(iss).ne.itrg(j)) cycle
                      js=ltns(iss)+nlim
                      wi(itrg(j))=wi(itrg(j))+wldna(js,iatm,0)
                      wwi(itrg(j))=wwi(itrg(j))+
     .                     ewlda(js,iatm,0)*wldna(js,iatm,0)
c                    print *,'j,iatm,ltns,js=',itrg(j),iatm,ltns(iss),js    !###
c                    print *,'wldna,ewlda=',wldna(js,iatm,0),
c     ,                                     ewlda(js,iatm,0)                !###
c                    print *,'wo,wwo=',wo,wwo                               !###
                    end do !}
                    if(wo.gt.0.0_R8) then !{
                      wi(itrg(j))=wi(itrg(j))/wo !units [1/s]/[1/s]
                    else !}{
                      wi(itrg(j))=0.0_R8
                    end if !}
                    if(wwo.gt.0.0_R8) then !{
                      wwi(itrg(j))=wwi(itrg(j))/wwo ! units [J]*[1/s]/[W]
                    else !}{
                      wwi(itrg(j))=0.0_R8
                    end if !}
                    do iy=1,mpg%divFcp(itrg(j),2) !{
                      faitrg(iy,itrg(j),iatm)=
     .                 faitrg(iy,itrg(j),iatm)*wi(itrg(j))
                      if(faitrg(iy,itrg(j),iatm).gt.0.0_R8) then !{
c carlis             [J]=[J/s]/[1/s]
                        eaitrg(iy,itrg(j),iatm)=
     .                   eaitrg(iy,itrg(j),iatm)*wwi(itrg(j))/
     /                   faitrg(iy,itrg(j),iatm)
                      else !}{
                        eaitrg(iy,itrg(j),iatm)=0.0_R8
                      end if !}
                    end do !} loop over boundary faces
                  else
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                      iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                      u = mpg%divFcOr(mpg%divFcP(itrg(j),1)+iy-1)
                      naitrg(iy,itrg(j),iatm)=
     .                 state%pl%na(iCv,eb2atcr(iatm))
                      paitrg(iy,itrg(j),iatm)=
     .                 state%pl%na(iCv,eb2atcr(iatm))*state%pl%ti(iCv)
                      faitrg(iy,itrg(j),iatm)=u*
     .                 (state%dv%fna(iFc,0,eb2atcr(iatm))+
     .                  state%dv%fna(iFc,1,eb2atcr(iatm)))
                      eaitrg(iy,itrg(j),iatm)=2.0_R8*state%pl%ti(iCv)/ev
                    end do !} loop over boundary faces
                  end if
                end if !}
              end do !} loop over targets
c              write(*,*) 'Atoms',iatm
c              write(*,*) 'Flux scaling'
c              write(*,'(1p,12e10.2)') (wi(itrg(i)),i=1,ntrgts)
c              write(*,*) 'Energy scaling'
c              write(*,'(1p,12e10.2)') (wwi(itrg(i)),i=1,ntrgts)
            end do !} loop over atoms

c*** Molecules fluxes and energies
            do imol=1,nmoli !{
              do j=1,2*max(1,mpg%nncut) !{ Do not use ntrgts here to allow for sign flips
                if (itrg(j).gt.0) then !{
                  if (ft44) then !{
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                      iCv = min(mpg%fcCv(iFc,1),mpg%fcCv(iFc,2))
                      ias = b2_fnmti(iFc)
                      fmitrg(iy,itrg(j),imol)=wldnm_res(imol,IAS)*
     .                 sarea_res(IAS)
                      pmitrg(iy,itrg(j),imol)=dmb2(iCv,imol,1)*
     *                 tmb2(iCv,imol,1)
                      nmitrg(iy,itrg(j),imol)=dmb2(iCv,imol,1)
                      emitrg(iy,itrg(j),imol)=ewldm_res(imol,iy)*
     .                 sarea_res(iy) !units in W/m2*m2
                    end do !} loop over boundary faces
                    wo=0.0_R8
                    wwo=0.0_R8
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      wo=wo+fmitrg(iy,itrg(j),imol) ! units 1/s
                      wwo=wwo+emitrg(iy,itrg(j),imol) ! units W
                    end do !}
                    wi(itrg(j))=0.0_R8
                    wwi(itrg(j))=0.0_R8
                    do iss=1,ntns !{
                      if (ktns(iss).ne.itrg(j)) cycle
                      js=ltns(iss)+nlim
                      wi(itrg(j))=wi(itrg(j))+wldnm(js,imol,0)
                      wwi(itrg(j))=wwi(itrg(j))+
     .                     ewldm(js,imol,0)*wldnm(js,imol,0)
                    end do !}
                    if(wo.gt.0.0_R8) then !{
                      wi(itrg(j))=wi(itrg(j))/wo !units [1/s]/[1/s]
                    else !}{
                      wi(itrg(j))=0.0_R8
                    end if !}
                    if(wwo.gt.0.0_R8) then !{
                      wwi(itrg(j))=wwi(itrg(j))/wwo ! units [J]*[1/s]/[W]
                    else !}{
                      wwi(itrg(j))=0.0_R8
                    end if !}
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      fmitrg(iy,itrg(j),imol)=
     .                 fmitrg(iy,itrg(j),imol)*wi(j)
                      if(fmitrg(iy,itrg(j),imol).gt.0.0_R8) then !{
                        emitrg(iy,itrg(j),imol)=
     .                   emitrg(iy,itrg(j),imol)*wwi(itrg(j))/
     /                   fmitrg(iy,itrg(j),imol)
                      else !}{
                        emitrg(iy,itrg(j),imol)=0.0_R8
                      end if !}
                    end do !}
                  else !}{
                    do iy=1,mpg%divFcP(itrg(j),2) !{
                      fmitrg(iy,itrg(j),imol)=0.0_R8
                      pmitrg(iy,itrg(j),imol)=0.0_R8
                      nmitrg(iy,itrg(j),imol)=0.0_R8
                      emitrg(iy,itrg(j),imol)=0.0_R8
                    end do !}
                  end if !}
                end if !}
              end do !} loop over targets
            end do !} loop over molecules

c*** Sputtering data
c carlis wldspta and wldsptm give sputtered flux for each emitted atom and molecule
            if (ntargsp.gt.0.and.ft44) then
              shift=lppp+1
              do j=1,ntrgts !{
                sjs_t(itrg(j))=0.0_R8
                sjs(:,itrg(j))=0.0_R8
                sptlistrg(itrg(j))=0
                do iy=1,mpg%divFcP(itrg(j),2) !{
                  iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                  ias = b2_fnmti(iFc)
                  do ia=1,natmi
                    if (wldspta_res(ia,IAS).gt.0.0_R8) then
                      sptlistrg(itrg(j))=ia
                      sjs(iy,itrg(j))=wldspta_res(ia,IAS)*
     .                 sarea_res(IAS)
                      sjs_t(itrg(j))=sjs_t(itrg(j))+sjs(iy,itrg(j))
                    end if
                  end do
                  do imol=1,nmoli
                    if (wldsptm_res(imol,IAS).gt.0.0_R8) then
                      sptlistrg(itrg(j))=imol+natmi
                      sjs(iy,itrg(j))=wldsptm_res(imol,IAS)*
     .                 sarea_res(IAS)
                      sjs_t(itrg(j))=sjs_t(itrg(j))+sjs(iy,itrg(j))
                    end if
                  end do
                end do

                if(sptlistrg(itrg(j)).ne.0) then
                  sjs_w(itrg(j))=0.0_R8
                  do iss=1,ntns
                    if (ktns(iss).ne.itrg(j)) cycle
                    js=ltns(iss)+nlim
                    if (sptlistrg(itrg(j)).le.natmi) then
                      sjs_w(itrg(j))=sjs_w(itrg(j))+
     .                 wldspta(js,sptlistrg(itrg(j)),0)
                    else
                      sjs_w(itrg(j))=sjs_w(itrg(j))+
     .                 wldsptm(js,sptlistrg(itrg(j))-natmi,0)
                    end if
                  end do
                  if(sjs_t(itrg(j)).gt.0.0_R8) then
                    sjs_t(itrg(j))=sjs_w(itrg(j))/sjs_t(itrg(j))
                  else
                    sjs_t(itrg(j))=0.0_R8
                  end if
                else
                  sjs_t(itrg(j))=0.0_R8
                end if

                if (sjs_t(itrg(j)).gt.0.0_R8) then
                  do_spt=.true.
                  do iy=1,mpg%divFcP(itrg(j),2) !{
                    iFc = mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
                    ias = b2_fnmti(iFc)
                    sjs(iy,itrg(j))=sjs_t(itrg(j))*sjs(iy,itrg(j))
c carlis The emitted amount in g/y is evaluated using the sputtered atom or molecule
c carlis        [g/y]=[1/s]*[s/y]/[1/mol]*[g/mol]
                    if (sptlistrg(itrg(j)).le.natmi) then
                      sptmass=am(eb2atcr(sptlistrg(itrg(j))))
                    else
                      imol = sptlistrg(itrg(j))-natmi
                      sptmass=0.0_R8
                      do ia = 1, natmi
                        if (mlcmp(ia,imol).gt.0) then
                          sptmass=sptmass+mlcmp(ia,imol)*am(eb2atcr(ia))
                        end if
                      end do
                    end if
                    spttrg(iy,itrg(j))=
     .               sjs(iy,itrg(j))*nsecyr/avogr*sptmass
c carlis The erosion of the wall is evaluated using the wall material and not the emitted atom/molecule
c carlis           [mm/y] = [1/s]*[s*mm*m2/y]*
                    erntrg(iy,itrg(j))=sjs(iy,itrg(j))*
     .                wllersf(i_ersmtr(iy+shift))/sarea_res(IAS)
                  end do
                  if (sptlistrg(itrg(j)).le.natmi) then
                    ssmol=texta(sptlistrg(itrg(j)))
                  else
                    ssmol=textmn(sptlistrg(itrg(j))-natmi-1)
                    ssmol=adjustl(ssmol)
                  end if
                  print *,'target material = ', mtrtrgsg(itrg(j))
                  print *,'sputtered atom from target = ',
     ,                ssmol
                  write(*,'(1x,a,1p,16e12.3)')
     .             'Scaling factors for '//mtrtrgsg(itrg(j))//' Sp : ',
     ,             (sjs_t(i),i=1,ntrgts)
                else ! no sputtering
                  do_spt=.false.
                  print *,
     ,             'b2pwldw, zero sputtering flux for target j=',j
                  sjs_w(itrg(j))=0.0_R8
                end if
                shift=shift+mpg%divFcP(itrg(j),2)
              end do !} loop on targets for sputtering
            end if !if sputtering happening
          end if !} jtrg.eq.0 .and. do_trg
          jtrg=itg
          iznt=izn
          iy=ils-lppp
          do j=1,itg-1
            if (itrg(j).eq.0) cycle
            iy=iy-mpg%divFcP(itrg(j),2)
          end do
        else if (iznt.ne.izn) then !}{ itg.ne.jtrg
          iy=ils-lppp
          do j=1,itg-1 !{
            if (itrg(j).eq.0) cycle
            iy=iy-mpg%divFcP(itrg(j),2)
          end do !}
          iznt=izn
        end if !} iznt.ne.izn
      end if !} ils.gt.lppp

c      print *,'izns,izn=',izns,izn      !###
      if(izns.ne.izn) then !{
        izns=izn
c*** Print the column headers
        hdr_line='#'
        k=1
        do j=1,n_hdpr_c !{
          l=len_trim(hdpr_c(j))
          hdr_line=hdr_line(1:k)//hdpr_c(j)(1:l)//tab
          k=k+l+1
        end do !}
        if(trgl) then !{
          do j=1,n_hdpr_t !{
            l=len_trim(hdpr_t(j))
            hdr_line=hdr_line(1:k)//hdpr_t(j)(1:l)//tab
            k=k+l+1
          end do !}
        end if !}
        do i=1,n_sfpr_s !{
          m=len_trim(sfpr_s(i))
          do j=1,n_hdpr_s !{
            l=len_trim(hdpr_s(j))
            hdr_line=hdr_line(1:k)//hdpr_s(j)(1:l)//'_'//
     /                              sfpr_s(i)(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
        do i=1,nclm_na !{
          m=len_trim(sfpr_a(i))
          do j=1,n_hdpr_a !{
            l=len_trim(hdpr_a(j))
            hdr_line=hdr_line(1:k)//hdpr_a(j)(1:l)//'_'//
     /                              sfpr_a(i)(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
        do i=1,nclm_nm !{
          ssmol=textmn(i-1)
          ssmol=adjustl(ssmol)
          m=len_trim(ssmol)
          do j=1,n_hdpr_ml !{
            l=len_trim(hdpr_ml(j))
            hdr_line=hdr_line(1:k)//hdpr_ml(j)(1:l)//'_'//
     /                              ssmol(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
        do j=1,n_sfpr_sp !{
          m=len_trim(sfpr_sp(j))
          do i=1,n_hdpr_sp !{
            l=len_trim(hdpr_sp(i))
            hdr_line=hdr_line(1:k)//hdpr_sp(i)(1:l)//'_'//
     /                              sfpr_sp(j)(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
        do j=1,n_sfpr_er !{
          m=len_trim(sfpr_m(j))
          do i=1,n_hdpr_er !{
            l=len_trim(hdpr_er(i))
            hdr_line=hdr_line(1:k)//hdpr_er(i)(1:l)//'_'//
     /                              sfpr_m(j)(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
c        print *,'About to write a header:',k-1,' symbols'   !###
c        print *,hdr_line(1:k-1)
                           !###
        write(lun,740) 'x',tab,'coordinate along plot zone [m].'
        write(lun,740) 'Wtot',tab,'total power flux [W/m^2]'
        write(lun,740) 'Wrad',tab,'power flux from radiation [W/m^2]'
        write(lun,740) 'Wpart',tab,'power flux from particles [W/m^2]'
        write(lun,740) 'Wpls',tab,
     .   'power flux from charged particles [W/m^2]'
        write(lun,740) 'Wneut',tab,
     .   'power flux from neutral particles [W/m^2]'
        write(lun,740) 'Wheat',tab,
     .   'power flux due to particles impinging on wall [W/m^2]'
        write(lun,740) 'Wpot',tab,
     .   'power flux due to recombination [W/m^2]'
        write(lun,740) 'Whtpl',tab,
     .   'power flux due to plasma particles impinging on wall [W/m^2]'
        write(lun,740) 'Wptpl',tab,
     .   'power flux due to plasma recombination [W/m^2]'
        write(lun,740) 'flxi_XX',tab,
     .   'flux of ions of species XX [1/s]'
        write(lun,740) 'Eavi_XX',tab,
     .   'average energy of ions of species XX [eV]'
        write(lun,740) 'flxa_XX',tab,
     .   'flux of atoms of species XX [1/s]'
        write(lun,740) 'Eava_XX',tab,
     .   'average energy of atoms of species XX [eV]'
        write(lun,740) 'na_XX',tab,
     .   'density of neutral atoms of species XX [1/m^3]'
        write(lun,740) 'pa_XX',tab,
     .   'pressure of neutral atoms of species XX [Pa]'
        if (nclm_nm.gt.0) then
          write(lun,740) 'flxm_XX',tab,
     .     'incident flux of molecules XX [1/m^2/s]'
          write(lun,740) 'Eavm_XX',tab,
     .     'average energy of molecules XX [eV]'
          write(lun,740) 'pm_XX',tab,
     .     'density of molecules XX [1/m^3]'
          write(lun,740) 'pm_XX',tab,
     .     'pressure of molecules XX [Pa]'
        end if
        write(lun,740) 'sptne_XX',tab,
     .   'net sputtering of wall material XX [g/y]'
        write(lun,740) 'sptng_XX',tab,
     .   'gross sputtering of wall material XX [g/y]'
        write(lun,740) 'ersne_XX',tab,
     .   'net erosion of wall material XX [mm/y]'
        write(lun,740) 'ersng_XX',tab,
     .   'gross erosion of wall material XX [mm/y]'
  740   format('# ',2a,' : ',a)
        write(lun,'(a)') hdr_line(1:k-1)
      end if !}
c*** Print the line
      allocate(w_line(n_line_x))
      m=lelmind(ils)
      k=n_hdpr_c
      if(itg.eq.0) then !{ WALL not target
        w_line(10)=0.0 !recombination due to plasma
        w_line(9)=0.0 !heat due to plasma
        w_line(8)=wldnep(m,0)/segm_area(ils)
        w_line(7)=wldnek(m,0)/segm_area(ils)
        w_line(6)=w_line(7)+w_line(8)
        w_line(5)=0 !load from plasma
        w_line(4)=w_line(6)
        w_line(3)=wldrdld(ils)
        w_line(2)=w_line(3)+w_line(4)
        w_line(1)=wldrfdst(iln,izn)
        do is=1,n_sfpr_s !{
          w_line(k+1)=0.0 !ions flux
          w_line(k+2)=0.0 !average ion energy
          k=k+n_hdpr_s
        end do !}
        do is=1,nclm_na !{
          w_line(k+1)=wldna(m,is,0)/segm_area(ils)
          w_line(k+2)=ewlda(m,is,0)
          w_line(k+3)=0.0 !atom density
          w_line(k+4)=0.0 !atom pressure
          k=k+n_hdpr_a
        end do !}
        do is=1,nclm_nm !{
          w_line(k+1)=wldnm(m,is,0)/segm_area(ils) !molecule flux
          w_line(k+2)=ewldm(m,is,0) !molecule energy
          w_line(k+3)=0.0 !molecule density
          w_line(k+4)=0.0 !molecule pressure
          k=k+n_hdpr_ml
        end do !}
        if (do_spt_wll(m).gt.0) then
          do is=1,n_sfpr_sp
            if(do_spt_wll(m).eq.i_sfpr_sp(is)) then !{
c carlis For output sputtering I assume that the sputtered flux is acting on the sputtered atom
c carlis        [g/y]=[1/s]*[s/y]/[1/mol]*[g/mol]
              if (do_spt_wll(m).le.natmi) then
                sptmass=am(eb2atcr(do_spt_wll(m)))
              else
                imol = do_spt_wll(m)-natmi
                sptmass=0.0_R8
                do ia = 1, natmi
                  if (mlcmp(ia,imol).gt.0) then
                    sptmass=sptmass+mlcmp(ia,imol)*am(eb2atcr(ia))
                  end if
                end do
              end if
              w_line(k+1)=0.0_R8 !net sputtering
              w_line(k+2)=sptflx(m)*nsecyr/avogr*sptmass !gross sputtering
              l=i_sfpr_sp(l_ersmtr(ils))
              if(l.gt.0) then !{
                wo=0.0_R8
                do i=1,natmi !{
                  if (latmscl(i).eq.l)
     >             wo=wo-wldna(m,i,0)+wldra(m,i,0)+wldpa(m,i,0)
                end do !}
                do i=1,nmoli !{
                  wo=wo-mlcmp(l,i)*(wldnm(m,i,0)-
     -                  wldrm(m,i,0)-wldpm(m,i,0))
                end do !}
                i1=0
                do i=1,l-1 !{
                  i1=i1+nfluids(i)
                end do !}
                do i=i1+1,i1+nfluids(l) !{
                  wo=wo-wldpp(m,i,0)
                end do !}
                w_line(k+1)=wo*nsecyr/avogr*sptmass
              end if
            else
              w_line(k+1)=0.0_R8
              w_line(k+2)=0.0_R8
            end if
            k=k+n_hdpr_sp
          end do
          do is=1,n_sfpr_er
            if(i_ersmtr(ils).eq.is) then !{
c carlis For output erosion I assume that the sputtered flux is acting on the target material
c carlis        [mm/y] = [1/s]*[s*mm*m2/y]/[m2]
              w_line(k+1)=0.0_R8 !net erosion
              w_line(k+2)=sptflx(m)*wllersf(is)/segm_area(ils) !gross erosion
              l=i_sfpr_sp(l_ersmtr(ils))
              if(l.gt.0) then !{
                wo=0.0_R8
                do i=1,natmi !{
                  if (latmscl(i).eq.l)
     >             wo=wo-wldna(m,i,0)+wldra(m,i,0)+wldpa(m,i,0)
                end do !}
                do i=1,nmoli !{
                  do ia=1,natmi !{
                    if (latmscl(ia).eq.l.and.mlcmp(ia,i).ne.0)
     >               wo=wo-mlcmp(ia,i)*(wldnm(m,i,0)-
     -                                  wldrm(m,i,0)-wldpm(m,i,0))
                  end do !}
                end do !}
                i1=0
                do i=1,l-1 !{
                  i1=i1+nfluids(i)
                end do !}
                do i=i1+1,i1+nfluids(l) !{
                  wo=wo-wldpp(m,i,0)
                end do !}
                w_line(k+1)=wo*wllersf(is)/segm_area(ils)
              end if
            else
              w_line(k+1)=0.0_R8
              w_line(k+2)=0.0_R8
            end if
            k=k+n_hdpr_er
          end do
        else
c carlis sputtering not happening on this wall but possibly on other walls of the same plot zone
          do is=1,n_sfpr_sp
            w_line(k+1)=0.0_R8
            w_line(k+2)=0.0_R8
            k=k+n_hdpr_sp
          end do
          do is=1,n_sfpr_er
            w_line(k+1)=0.0_R8
            w_line(k+2)=0.0_R8
            k=k+n_hdpr_er
          end do
        end if
      else if(do_trg) then !}{
c        print *,'jtrg,iy=',jtrg,iy
        w_line(10)=wptpl(iy,jtrg) !recombination due to plasma
        w_line(9)=whtpl(iy,jtrg) !heat due to plasma                           !###
        w_line(8)=wptnu(iy,jtrg)+wptpl(iy,jtrg)
        w_line(7)=whtnu(iy,jtrg)+whtpl(iy,jtrg)
        w_line(6)=wnutrg(iy,jtrg)
        w_line(5)=wpltrg(iy,jtrg)
        w_line(4)=wpttrg(iy,jtrg)
        w_line(3)=wrdtrg(iy,jtrg,0)
        w_line(2)=wtttrg(iy,jtrg)
        w_line(1)=wldrfdst(iln,izn)
        do is=1,n_sfpr_s !{
          w_line(k+1)=fiitrg(iy,jtrg,is)
          w_line(k+2)=eiitrg(iy,jtrg,is)/eV
          k=k+n_hdpr_s
        end do !}
        do is=1,nclm_na !{
          w_line(k+1)=faitrg(iy,jtrg,is)
          w_line(k+2)=eaitrg(iy,jtrg,is)
          w_line(k+3)=naitrg(iy,jtrg,is)
          w_line(k+4)=paitrg(iy,jtrg,is)
          k=k+n_hdpr_a
        end do !}
        do is=1,nclm_nm !{
          w_line(k+1)=fmitrg(iy,jtrg,is) !molecule flux
          w_line(k+2)=emitrg(iy,jtrg,is) !molecule energy
          w_line(k+3)=nmitrg(iy,jtrg,is) !molecule density
          w_line(k+4)=pmitrg(iy,jtrg,is) !molecule pressure
          k=k+n_hdpr_ml
        end do !}
        if (do_spt) then
          do is=1,n_sfpr_sp
            if(sptlistrg(jtrg).eq.i_sfpr_sp(is)) then !{
              w_line(k+1)=0.0_R8 !net sputtering
              w_line(k+2)=spttrg(iy,jtrg) !gross sputtering
            else
              w_line(k+1)=0.0_R8
              w_line(k+2)=0.0_R8
            end if
            k=k+n_hdpr_sp
          end do
          do is=1,n_sfpr_er
            if(i_ersmtr(ils).eq.is) then !{
              w_line(k+1)=0.0_R8 !net erosion
              w_line(k+2)=erntrg(iy,jtrg) !gross erosion
            else
              w_line(k+1)=0.0_R8
              w_line(k+2)=0.0_R8
            end if
            k=k+n_hdpr_er
          end do
        else
c carlis sputtering not happening on this target but possibly on other walls of the same plot zone
          do is=1,n_sfpr_sp
            w_line(k+1)=0.0_R8
            w_line(k+2)=0.0_R8
            k=k+n_hdpr_sp
          end do
          do is=1,n_sfpr_er
            w_line(k+1)=0.0_R8
            w_line(k+2)=0.0_R8
            k=k+n_hdpr_er
          end do
        end if
        if(trgl) then !{
          w_line(1)=xtrg(iy,jtrg)
          w_line(k+1)=netrg(iy,jtrg)
          w_line(k+2)=tetrg(iy,jtrg)
          w_line(k+3)=titrg(iy,jtrg)
          w_line(k+4)=xmp_hlp(iy,jtrg)
          k=k+n_hdpr_t
        end if !}
        iy=iy+1
      else !}{
        k=n_hdpr_c+n_sfpr_s+n_sfpr_sp
        if(trgl) then !{
          k=k+n_hdpr_t
          w_line(1)=xtrg(iy,jtrg)
        else !}{
          w_line(1)=wldrfdst(iln,izn)
        end if !}
        do i=2,k !{
          w_line(i)=0.0_R8
        end do !}
      end if !}
c      print *,'About to write a line:',k,' values'   !###
      if(k.gt.0) write(lun,frmt_hlp) (w_line(i),tab,i=1,k)
      deallocate(w_line)
c      print *,'TEXTMN*',adjustl(TEXTMN),'*'
c      print *,'TEXTIN*',adjustl(TEXTIN),'*'
      return

c=======================================================================
      end subroutine b2pwldw
