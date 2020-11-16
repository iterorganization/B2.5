      subroutine b2pwldw(lun,trgl,izn,iln,BoRiS)
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_rates
      use b2mod_diag
      use b2mod_external
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_extrab25
      use eirmod_ccona
      use eirmod_cinit
#else
      use b2mod_eirbra
#endif
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
#include <DIMENSIONS.F>
#ifndef B25_EIRENE
      real (kind=R8) :: EionH2
      parameter (EionH2=15.4_R8)
      character*8 TEXTA(DEF_NATM)
#else
      real (kind=R8) :: pump, pumptot
      character*13 species_name
#endif
      integer lun,izn,iln
      real (kind=R8) :: BoRiS
      logical trgl,ft44

c*** Local variables
c***  ils     : the segment index in the whole list
c***  itrg    : target number itrg prepared
c***  izns    : plot zone izns started
c***  iznt    : plot zone of the last done target
c***  zz      : reference distance
c***  i_ersmtr: index of eroded material for the wall segments

      integer i,j,k,l,m,n,i1,i2,ib,is,iss,iy,im,js,kk,ils,itg,
     ,  itrg,izns,iznt,
     ,  i_ersmtr(nsgmx),l_ersmtr(0:nsgmx),
     ,  iatm,imol,ia,ims,ind,IAS,shift,
     ,  sptlistrg(0:nsgmx)

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

      real (kind=R8) :: h,p,t,u,wo,wwo,eni_hlp(ntrgsx),mion,
     ,  phi_hlp(ntrgsx),wi(ntrgsx),wwi(ntrgsx),z_hlp(ntrgsx),
     ,  s(0:NY-1,ntrgsx),xmp_hlp(0:NY-1,ntrgsx),xtrg(0:NY-1,ntrgsx),
     ,  netrg(0:NY-1,ntrgsx),  tetrg(0:NY-1,ntrgsx),
     ,  titrg(0:NY-1,ntrgsx),
     ,  wtttrg(0:NY-1,ntrgsx), wptnu(0:NY-1,ntrgsx),
     ,  wptpl(0:NY-1,ntrgsx), wnutpr(0:NY-1,ntrgsx),
     ,  whtnu(0:NY-1,ntrgsx),  whtpl(0:NY-1,ntrgsx),
     ,  wpttrg(0:NY-1,ntrgsx),
     ,  wpltrg(0:NY-1,ntrgsx), wnutrg(0:NY-1,ntrgsx),
     ,  wrdtrg(0:NY-1,ntrgsx,0:DEF_NATM),
     ,  fiitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  eiitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  faitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  eaitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  naitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  paitrg(0:NY-1,ntrgsx,DEF_NATM),
     ,  fmitrg(0:NY-1,ntrgsx,DEF_NMOL),
     ,  emitrg(0:NY-1,ntrgsx,DEF_NMOL),
     ,  pmitrg(0:NY-1,ntrgsx,DEF_NMOL),
     ,  nmitrg(0:NY-1,ntrgsx,DEF_NMOL),
     ,  erntrg(0:NY-1,ntrgsx), spttrg(0:NY-1,ntrgsx),
     ,  sjs(0:NY-1,ntrgsx), sjs_t(ntrgsx),
     ,  sjs_w(0:ntrgsx),sptmass,sptflx(DEF_NLIM+DEF_NSTS)
c     real (kind=R8) :: zz

c*** Recap of fluxes from Eirene
c***    wldnek  :   heat transferred with neutrals. Evaluated as EOTAT+EOTML+EOTIO
c***                   -erfaat-erfmat-erfiat-erfaml-erfmml-erfiml
c***    wldnep  :   potential energy released by neutrals. Evaluated as diss_pot_H2*prfaml
c***    wldna   :   flux of atoms impinging onto the surface. Equal to POTAT
c***    ewlda   :   their average energy. Equal to EOTAT/POTAT
c***    wldnm   :   flux of molecules impinging onto the surface. Equal to POTML
c***    ewldm   :   their average energy. Equal to EOTML/POTML
c***    wldra   :   flux of reflected atoms. Equal to PRFAAT+PRFMAT+PRFIAT
c***    wldrm   :   flux of reflected molecules. Equal to PRFAML+PRFMML+PRFIML
c***    wldpp   :   flux of plasma ions impinging onto the surface. Equal to POTPL
c***    wldpa   :   flux of emitted atoms from BI. Equal to PRFPAT
c***    wldpm   :   flux of emitted molecules from BI. Equal to PRFPML
c***    wldpeb  :   power carried away by these atoms and molecules. Equal to ERFPAT+ERFPML+ERFPIO
c***    wldspt  :   flux of sputtered wall material
c***    wldspta :   sputtered flux by emitted atom
c***    wldsptm :   sputtered flux by emitted molecule
c***    wlarea  :   areas of the surface segments from Eirene
c***    wlabsrp :   absorption at the surfaces (1-recyct from Eirene)
c***    wlpump  :   pumped flux for each pumped particle
c***    ewldea_res  :   emitted energy flux of atoms. Equal to ERFAAT+ERFMAT+ERFIAT+ERFPHAT
c***    ewldem_res  :   emitted energy flux due of molecules. Equal to ERFAML+ERFMML+
c***                      ERFIML+ERFPHML
c***    ewldrp_res  :   kinetic energy of reflected atoms+molecules originated from ions. ERFPAT+ERFPML
c***    ewldmr_res  :   energy due to recombination of atoms into molecules. ()PRFAML+PRFPML)*H_EDISS

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
      real (kind=R8), allocatable :: w_line(:)
      real (kind=R8) :: nsecyr
      parameter (nsecyr=24.0_R8*3600.0_R8*365.25_R8*0.9925_R8)
      logical :: first,do_trg,match_found,do_spt
      character*(m_hdpr) tab*1,hdpr_c(n_hdpr_c),hdpr_t(n_hdpr_t),
     ,  hdpr_s(n_hdpr_s),hdpr_a(n_hdpr_a),hdpr_sp(n_hdpr_sp),
     ,  sfpr_s(DEF_NATM)*2,sfpr_a(DEF_NATM)*2,
     ,  sfpr_m(nsptmx)*(lnwlmtrl),hdr_line*(m_hdln),frmt_hlp*24,
     ,  hdpr_ml(n_hdpr_ml), hdpr_er(n_hdpr_er),
     ,  sfpr_sp(DEF_NATM+DEF_NMOL)*3
      character*10 ssmol,ssion
      data hdpr_c /'x', 'Wtot','Wrad','Wpart','Wpls',
     2  'Wneu','Wheat','Wpot', 'Whtpl', 'Wptpl' /
      data hdpr_t /'ne', 'Te','Ti','xMP'/
      data hdpr_s /'flxi','Eavi'/
      data hdpr_a /'flxa','Eava','na','pa'/
      data hdpr_ml /'flxm','Eavm','nm','pm'/
      data hdpr_sp /'sptne','sptng'/
      data hdpr_er /'ersne','ersng'/
      data itrg, izns, iznt, n_sfpr_s /4*0/
      data first /.true./
      data sfpr_s, sfpr_a, sfpr_m /
     /     DEF_NATM*' ', DEF_NATM*' ', nsptmx*' ' /

c*** Data for erosion calculations
c***  n_srfs_m  : number of possible wall materials
c***  srfs_m    : list of possible wall materials
c***  srfs_x    : corresponding species (say, Fe for SS)
c***  srfs_d    : specific masses of possible wall materials [g/cm**3]
c***  srfs_a    : atomic weight of possible wall materials
c***  wllersf   : ratio of erosion to sputtering for different meterials
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
        nclm_na=natmi
        nclm_nm=nmoli
#endif
        do ia=1,nclm_na
          sfpr_a(ia)=adjustl(TEXTA(ia))
        end do
c        print *,'First call. format:',frmt_hlp,':'                  !###
c        print *,'n_sfpr_s=',n_sfpr_s                                !###
c        if(n_sfpr_s.gt.0) print *,(' ',sfpr_s(i),i=1,n_sfpr_s)      !###
        k=0
        do i=1,nbc !{
          l=bcpos(i)
          if(bcchar(i).eq.'W') then !{
            k=k+1
            iy=(bcstart(i)+bcend(i))/2
            if(l.eq.-1) then !{
              phi_hlp(1)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(1)=enipar(i,1)+2.0_R8-enipar(i,2)
            else !}{
              phi_hlp(3)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(3)=enipar(i,1)+2.0_R8-enipar(i,2)
            end if !}
          else if(bcchar(i).eq.'E') then !}{
            k=k+1
            iy=(bcstart(i)+bcend(i))/2
            if(l.eq.nx) then !{
              if(ntrgts.eq.1) then !{
                phi_hlp(1)=get_phi_applied(boundary_namelist,l,iy)
                eni_hlp(1)=enipar(i,1)+2.0_R8-enipar(i,2)
              else if(ntrgts.eq.2) then !{
                phi_hlp(2)=get_phi_applied(boundary_namelist,l,iy)
                eni_hlp(2)=enipar(i,1)+2.0_R8-enipar(i,2)
              else if (ntrgts.eq.4) then !}{
                phi_hlp(4)=get_phi_applied(boundary_namelist,l,iy)
                eni_hlp(4)=enipar(i,1)+2.0_R8-enipar(i,2)
              end if !}
            else if(ntrgts.eq.4) then !}{
              phi_hlp(2)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(2)=enipar(i,1)+2.0_R8-enipar(i,2)
            end if !}
          end if !}
        end do !}
c        print '(1x,a,1p,4e12.3)','phi_hlp:',(phi_hlp(i),i=1,ntrgts) !###
c        print '(1x,a,1p,4e12.3)','eni_hlp:',(eni_hlp(i),i=1,ntrgts) !###
        if(k.ne.ntrgts) then !{
          write(0,*) 'b2pwldw: inconsistent target specification: ',
     ,                                                         ntrgts,k
          do_trg=.false.
        end if !}

c*** data for sputtering and erosion

c        print *,'Input wall materials'                              !###
c        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=1,lppp)            !###
c        print *,'Input target materials'                            !###
c        print '(1h:,20(a,1h:))',(mtrtrgsg(i),i=1,ntrgts)            !###
        k=lppp
        do j=1,ntrgts !{
          do i=k+1,k+ny !{
c carlis insert wall material for each target element
            mtrwldsg(i)=mtrtrgsg(j)
          end do !}
          k=k+ny
        end do !}
c        print *,'k,lppp,ny,ntrgts=',k,lppp,ny,ntrgts                !###
c        print *,'Target segment materials'                          !###
c        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=lppp+1,k)          !###

        n_sfpr_er=0
        do k=1,lppp+ny*ntrgts !{
          l=0
          do i=1,n_sfpr_er !{
c carlis l is index of wall material referred to list of ACTIVE materials
            if(streql(mtrwldsg(k),'dummy')) cycle
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
c carlis list of indexes of wall materials referred to list of ACTIVE materials
          i_ersmtr(k)=l
        end do !}
c        print *,'Input materials: ',sfpr_m                        !###
c        if(n_sfpr_m.gt.0) print *,':',(sfpr_m(i),':',i=1,n_sfpr_m)  !###
        l_ersmtr(0)=0
        m=0
        do k=1,n_sfpr_er !{
          l=0
          do i=1,n_srfs_m !{
c carlis l is index of wall material in ACTIVE list referred to list of POSSIBLE materials
            if(streql(srfs_m(i),sfpr_m(k))) l=i
          end do !}
          if(l.eq.0) then !{
            write(*,*) 'b2pwldw: unknown wall material ',sfpr_m(k)
            l_ersmtr(k)=0

          else !}{
            m=m+1
c carlis in list of ACTIVE materials we put the 'real' material e.g. CFC becomes C
c        the list is not erdered as the list of possible materials
            sfpr_m(m)=srfs_x(l)
c carlis     [s/y*m2*mm]   =[s/y]*[m2/cm2*mm/cm*]*[g/mol]/[1/mol]/[g/cm3]
            wllersf(m)=nsecyr*1.e-3*srfs_a(l)/(avogr*srfs_d(l))
c carlis list of indexes that tells the new ordering of wall materials in active list wrt the ordering of i_ersmtr
            l_ersmtr(k)=m
          end if !}
        end do !}
        spttrg(:,:)=0.0_R8
        erntrg(:,:)=0.0_R8
        n_sfpr_er=m
        do i=1,lppp+ny*ntrgts !{
c carlis substitute indexes of wall materials with the new ordering of ACTIVE list
          i_ersmtr(i)=l_ersmtr(i_ersmtr(i))
        end do !}
        write(*,*) 'Number of wall materials: ',n_sfpr_er
        if(n_sfpr_er.gt.0) then !{
          write(*,*) 'Wall materials : ',(sfpr_m(i),i=1,n_sfpr_er)
          write(*,'(a,1p,10e10.2)') ' Erosion factors: ',
     ,                                        (wllersf(i),i=1,n_sfpr_er)
        end if !}
c        print *,'i_sptmtr'                                    !###
c        print '(1x,50i2)',(i_sptmtr(i),i=1,lppp+ny*ntrgts)    !###

        do_spt_wll(:)=0
        sptflx(j)=0.0_R8
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
            sfpr_sp(i)=trim(TEXTA(ia))
          end if
        end do
        do imol=1,nmoli
          match_found=.false.
          j=1
          do j=1,nlimps
            if (wldsptm(j,imol,0).gt.0.0_R8) then
              sptflx(j)=wldsptm(j,imol,0)
              do_spt_wll(j)=imol+natm
              match_found=.true.
            end if
          end do
          if (match_found) then
            i=i+1
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
        do is = 0, ns-1
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
          itg=(ils-lppp-1)/ny+1
        end if !}
c        print *,'itg,iznt,itrg=',itg,iznt,itrg  !###
c*** Check whether this target has been prepared
        if(itg.ne.itrg) then !{

          if(itrg.eq.0 .and. do_trg) then !{

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
              k=ltrgpos(j)
              do iy=0,lseppos(j) !{
                z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                                   (cry(k,iy,3)-cry(k,iy,1))**2)
              end do !}
              do iy=0,ny-1 !{
                h=sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                 (cry(k,iy,3)-cry(k,iy,1))**2)
                s(iy,j)=gs(k+1,iy,0)
                z_hlp(j)=z_hlp(j)+0.5_R8*h
                xtrg(iy,j)=z_hlp(j)
                z_hlp(j)=z_hlp(j)+0.5_R8*h
              end do !}
            end do !}

            if(trgl) then !{

c*** Midplane coordinates

              l=jxi
              m=0
              do j=1,ntrgts !{
                u=0.0_R8
                k=l
                n=2*m
                do i=1,lseppos(j) !{
                  h=0.5_R8*sqrt((crx(k,i,n+1)-crx(k,i,n))**2+
     +                          (cry(k,i,n+1)-cry(k,i,n))**2)
                  u=u+h
                  xmp_hlp(i,j)=u
                  u=u+h
                end do !}
                do i=1,lseppos(j) !{
                  xmp_hlp(i,j)=xmp_hlp(i,j)-u
                end do !}
                u=0.0_R8
                do i=lseppos(j)+1,ny-1 !{
                  h=0.5_R8*sqrt((crx(k,i,n+1)-crx(k,i,n))**2+
     +                          (cry(k,i,n+1)-cry(k,i,n))**2)
                  u=u+h
                  xmp_hlp(i,j)=u
                  u=u+h
                end do !}
                if(j.eq.ntrgts/2) l=jxa
                m=1-m
              end do !}
            end if !}

c*** Power loading:

c*** Plasma and neutrals

c            print *,'Starting target power. ntrgts=',ntrgts   !###
            l=1
            u=-1
            do j=1,ntrgts !{
              if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
                u=1
                l=0
              end if
              k=ltrgpos(j)
              kk=k+l
              ISS=ltns(j)
              if (ISS.eq.0) cycle
              IND=eirdiag_nds_ind(ISS)
              IA=IND+1
              do iy=0,ny-1 !{
                wo=0.0_R8
                wnutrg(iy,j)=0.0_R8
                whtpl(iy,j)=fhe(k+1-l,iy,0)+fhi(k+1-l,iy,0)+
     .                      fhi_ext(k+1-l,iy,0)
                do is=0,ns-1 !{
                  if (.not.is_neutral(is).or.ft44) then
                    wo=wo+fhp(k+1-l,iy,0,is)
                  else
                    wnutrg(iy,j)=wnutrg(iy,j)+(1.0_R8-BoRiS)*
     *               (fhm(k+1-l,iy,0,is)+
     *                fna(k+1-l,iy,0,is)*ti(kk,iy))
                  end if
                  whtpl(iy,j)=whtpl(iy,j)+(1.0_R8-BoRiS)*
     *               (fhm(k+1-l,iy,0,is)+
     *                fna(k+1-l,iy,0,is)*ti(kk,iy))
c carlis 31jan2017  added ti as suggested by sw in b2blne.
                end do !}
                do is=0,ns_ext-1 !{
                  if (.not.is_ext_neutral(is).or.ft44) then
                    wo=wo+fa_ext(k+1-l,iy,0,is)*pt_ext(kk,iy,is)*ev
                  else
                    wnutrg(iy,j)=wnutrg(iy,j)+fa_ext(k+1-l,iy,0,is)
     .               *(1.0_R8-BoRiS)*(0.5_R8*am_ext(is)*mp*
     .                 ua_ext(kk,iy,is)**2+ta_ext(kk,iy,is))
                  end if
                  whtpl(iy,j)=whtpl(iy,j)+fa_ext(k+1-l,iy,0,is)
     .             *(1.0_R8-BoRiS)*(ta_ext(kk,iy,is)+0.5_R8*
     .             am_ext(is)*mp*ua_ext(kk,iy,is)**2)
c carlis 31jan2017  added ta as suggested by sw in b2blne.
                end do !}
                whtpl(iy,j)=whtpl(iy,j)+fne(k+1-l,iy,0)
     .           *(1.0_R8-BoRiS)*te(kk,iy)
c carlis 31jan2017  added electron contribution as suggested by sw in b2blne.
                wptpl(iy,j)=wo/s(iy,j)*u
                if (.not.ft44) wnutrg(iy,j)=wnutrg(iy,j)*u/s(iy,j)
                whtpl(iy,j)=whtpl(iy,j)/s(iy,j)*u
                wpltrg(iy,j)=wptpl(iy,j)+whtpl(iy,j)
                netrg(iy,j)=ne(kk,iy)
                tetrg(iy,j)=te(kk,iy)/ev
                titrg(iy,j)=ti(kk,iy)/ev
                if (ft44) then
                  do is=1,natmi
c                   carlis Using values directly calculated from Eirene
                    wnutrg(iy,j)=wnutrg(iy,j)
     .                           +(ewlda_res(is,iy+IA)
     .                           -ewldea_res(is,iy+IA))*
     .                            sarea_res(iy+IA)
                  end do
                  do is=1,nmoli
c                   carlis Using values directly calculated from Eirene
                    wnutrg(iy,j)=wnutrg(iy,j)
     .                           +(ewldm_res(is,iy+IA)
     .                           -ewldem_res(is,iy+IA)
     .                           +ewldmr_res(is,iy+IA))*
     .                            sarea_res(iy+IA)
                  enddo
                  wnutpr(iy,j)=ewldrp_res(iy+IA)
                  whtpl(iy,j)=whtpl(iy,j)-wnutpr(iy,j)
                  wpltrg(iy,j)=wpltrg(iy,j)-wnutpr(iy,j)
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
                    t=tib2(kk+1,iy+1,is,1)/ev
                    if (mion.gt.0.0_R8)
     >               wnutrg(iy,j)=wnutrg(iy,j)+
     .                s(iy,j)*dib2(kk+1,iy+1,is,1)*
     *                sqrt(max(t/mion,0.0_R8))*0.25_R8*1.0e4_R8*
     .                (2.0_R8*tib2(kk+1,iy+1,is,1)+EionH2*ev)
                  enddo
                end if
              end do !} loop over B2 cells
              l=1-l
              u=-u
            end do !} loop over targets

c*** Re-scale the neutral target loading and calculate the sums

            do j=1,ntrgts !{
              wi(j)=0.0_R8
              do i=0,ny-1 !{
                wi(j)=wi(j)+wnutrg(i,j)
              end do !}
            end do !}
c            print *,'ntrgts, wi: ',ntrgts,(wi(i),i=1,ntrgts)  !###
            do j=1,ntrgts !{
              if (ft44) then
                if(ntns.eq.ntrgts) then !{
                  js=ltns(j)+nlim
                  if (ltns(j).eq.0) cycle
                  wo=wldnek(js,0)+wldnep(js,0)
c                  print *,'j,js,wo = ',j,js,wo    !###
                  if(wi(j).gt.0.0_R8) then !{
                    wi(j)=wo/wi(j)
                  else !}{
                    wi(j)=0.0_R8
                  end if !}
                else !}{
                  write(*,*) 'b2pwldw: ntns.ne.ntrgts (',
     ,              ntns,' vs.',ntrgts,
     ,              ' ) - no neutral load onto the targets calculated'
                  write(*,*) '  ==> check specification of ltns in ',
     ,              'b2.neutrals.parameters'
                  wi(j)=0.0_R8
                end if !}
              endif
              do iy=0,ny-1 !{
                if (ft44) then
                  wnutrg(iy,j)=wnutrg(iy,j)*wi(j)/s(iy,j)
                  if(abs(wo).gt.0.0_R8) then !{
                    wptnu(iy,j)=wnutrg(iy,j)*wldnep(js,0)/wo
                    whtnu(iy,j)=wnutrg(iy,j)-wptnu(iy,j)
                  else !}{
                    wptnu(iy,j)=0.0_R8
                    whtnu(iy,j)=0.0_R8
                  end if !}
                end if
                wpttrg(iy,j)=wnutrg(iy,j)+wpltrg(iy,j)
                wtttrg(iy,j)=wpttrg(iy,j)+wrdtrg(iy,j,0)
              end do !}
            end do !}
            write(*,'(1x,a,1p,16e12.3)') 'Scaling factors for Wn : ',
     ,                                    (wi(i),i=1,ntrgts)

c*** Particle fluxes and energies
            i2=0
            do is=1,n_sfpr_s !{
              i1=eb2spcr(is)+1
              i2=eb2spcr(is)+nfluids(is)
              l=1
              u=-1
              do j=1,ntrgts !{
                if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
                  u=1
                  l=0
                end if
                k=ltrgpos(j)
                kk=k+l
                ib=k+1-l
                ISS=ltns(j)
                if (ISS.eq.0) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1
c!###{
c                print *,'is,j,i1,i2,k,kk,ib=',is,j,i1,i2,k,kk,ib
c                zz=0.0_R8
c                print '(3x,12(2x,a8))','phi','fti','tn','wo','wwo',
c     ,            'fiitrg','eiitrg','faitrg','eaitrg','dab2',
c     ,            'tab2','s'
c!###}
                do iy=0,ny-1 !{
                  wo=0.0_R8
                  wwo=0.0_R8
                  p=(po(ib,iy)-phi_hlp(j))*qe
                  h=ti(ib,iy)*eni_hlp(j)
                  do i=i1,i2 !{
                    wo=wo+fna(ib,iy,0,i)
                    wwo=wwo+(h+p*rza(kk,iy,i))*fna(ib,iy,0,i)
                  end do !}
c                  zz=zz+wwo*u                 !###
                  fiitrg(iy,j,is)=u*wo
                  if(abs(wo).gt.0.0_R8) then !{
                    eiitrg(iy,j,is)=wwo/wo
                  else !}{
                    eiitrg(iy,j,is)=0.0_R8
                  end if !}
c!###{
c                  print '(1x,i2,1p,8e10.2)',iy,p/ev,h/ev,t,wo,wwo,
c     ,              fiitrg(iy,j,is),eiitrg(iy,j,is),s(iy,j)
c!###}
                end do !} loop over B2 cells
c                print *,'j,is,zz=',j,is,zz  !###
                l=1-l
                u=-u
              end do !} loop over targets
c              write(*,*) 'Species',is
c              write(*,*) 'Flux scaling'
c              write(*,'(1p,12e10.2)') (wi(i),i=1,ntrgts)
c              write(*,*) 'Energy scaling'
c              write(*,'(1p,12e10.2)') (wwi(i),i=1,ntrgts)
            end do !} loop over active species

c*** Atomic fluxes and energies
            do iatm=1,nclm_na !{
              l=1
              u=-1
              do j=1,ntrgts !{
                if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
                  u=1
                  l=0
                end if
                k=ltrgpos(j)
                kk=k+l
                ib=k+1-l
                ISS=ltns(j)
                if (ISS.eq.0) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1
                do iy=0,ny-1 !{
                  if (ft44) then
                    naitrg(iy,j,iatm)=dab2(kk+1,iy+1,iatm,1)
                    paitrg(iy,j,iatm)=dab2(kk+1,iy+1,iatm,1)*
     *                                tab2(kk+1,iy+1,iatm,1)
c carlis using values directly from Eirene
                    faitrg(iy,j,iatm)=wldna_res(iatm,iy+IAS)*
     .                       sarea_res(iy+IAS) !units in 1/m2/s*m2
                    eaitrg(iy,j,iatm)=ewlda_res(iatm,iy+IAS)*
     .                       sarea_res(iy+IAS) !units in W/m2*m2
                  else
                    naitrg(iy,j,iatm)=na(kk,iy,eb2atcr(iatm))
                    paitrg(iy,j,iatm)=na(kk,iy,eb2atcr(iatm))*ti(kk,iy)
                    faitrg(iy,j,iatm)=u*fna(ib,iy,0,eb2atcr(iatm))
                    eaitrg(iy,j,iatm)=2.0_R8*ti(kk,iy)
                  end if
c!###{
c                  print '(1x,i2,1p,4e10.2)',iy,
c     ,              faitrg(iy,j,iatm),eaitrg(iy,j,iatm)/ev,
c     ,              dab2(kk+1,iy+1,iatm,1),tab2(kk+1,iy+1,iatm,1)/ev
c!###}
                end do !} loop over B2 cells
                wo=0.0_R8
                wwo=0.0_R8
                do iy=0,ny-1 !{
                  wo=wo+faitrg(iy,j,iatm) ! units 1/s
                  wwo=wwo+eaitrg(iy,j,iatm) ! units W
                end do !}
                js=ltns(j)+nlim
c                print *,'j,iatm,ltns,js=',j,iatm,ltns(j),js              !###
c                print *,'wldna,ewlda=',wldna(js,iatm,0),ewlda(js,iatm,0) !###
c                print *,'wo,wwo=',wo,wwo                                 !###
                if(wo.gt.0.0_R8) then !{
                  wi(j)=wldna(js,iatm,0)/wo !units [1/s]/[1/s]
                else !}{
                  wi(j)=0.0_R8
                end if !}
                if(wwo.gt.0.0_R8) then !{
                  wwi(j)=ewlda(js,iatm,0)*wldna(js,iatm,0)/wwo ! units [J]*[1/s]/[W]
                else !}{
                  wwi(j)=0.0_R8
                end if !}
                do iy=0,ny-1 !{
                  faitrg(iy,j,iatm)=faitrg(iy,j,iatm)*wi(j)
                  if(faitrg(iy,j,iatm).gt.0.0_R8) then !{
c carlis            [J]=[J/s]/[1/s]
                    eaitrg(iy,j,iatm)=eaitrg(iy,j,iatm)*wwi(j)/
     /                                faitrg(iy,j,iatm)
                  else !}{
                    eaitrg(iy,j,iatm)=0.0_R8
                  end if !}
                end do !} loop over B2 cells
                l=1-l
                u=-u
              end do !} loop over targets
c              write(*,*) 'Atoms',iatm
c              write(*,*) 'Flux scaling'
c              write(*,'(1p,12e10.2)') (wi(i),i=1,ntrgts)
c              write(*,*) 'Energy scaling'
c              write(*,'(1p,12e10.2)') (wwi(i),i=1,ntrgts)
            end do !} loop over atoms

c*** Molecules fluxes and energies
            do imol=1,nmoli !{
              l=1
              u=-1
              do j=1,ntrgts !{
                if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
                  u=1
                  l=0
                endif
                k=ltrgpos(j)
                kk=k+l
                ib=k+1-l
                ISS=ltns(j)
                if (ISS.eq.0) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1
                do iy=0,ny-1 !{
                  if (ft44) then
                    fmitrg(iy,j,imol)=wldnm_res(imol,iy+IAS)*
     .                     sarea_res(iy+IAS)
                    pmitrg(iy,j,imol)=dmb2(kk+1,iy+1,imol,1)*
     *                     tmb2(kk+1,iy+1,imol,1)
                    nmitrg(iy,j,imol)=dmb2(kk+1,iy+1,imol,1)
                    emitrg(iy,j,imol)=ewldm_res(imol,iy+IAS)*
     .                     sarea_res(iy+IAS) !units in W/m2*m2
                  else
                    fmitrg(iy,j,imol)=0.0
                    pmitrg(iy,j,imol)=0.0
                    nmitrg(iy,j,imol)=0.0
                    emitrg(iy,j,imol)=0.0
                  endif
                end do !} loop over B2 cells
                wo=0.0_R8
                wwo=0.0_R8
                do iy=0,ny-1 !{
                  wo=wo+fmitrg(iy,j,imol) ! units 1/s
                  wwo=wwo+emitrg(iy,j,imol) ! units W
                end do !}
                js=ltns(j)+nlim
                if(wo.gt.0.0_R8) then !{
                  wi(j)=wldnm(js,imol,0)/wo !units [1/s]/[1/s]
                else !}{
                  wi(j)=0.0_R8
                end if !}
                if(wwo.gt.0.0_R8) then !{
                  wwi(j)=ewldm(js,imol,0)*wldnm(js,imol,0)/wwo ! units [J]*[1/s]/[W]
                else !}{
                  wwi(j)=0.0_R8
                end if !}
                do iy=0,ny-1 !{
                  fmitrg(iy,j,imol)=fmitrg(iy,j,imol)*wi(j)
                  if(fmitrg(iy,j,imol).gt.0.0_R8) then !{
                    emitrg(iy,j,imol)=emitrg(iy,j,imol)*wwi(j)/
     /                                                fmitrg(iy,j,imol)
                  else !}{
                    emitrg(iy,j,imol)=0.0_R8
                  end if !}
                end do !}
                l=1-l
                u=-u
              end do !} loop over targets
            end do !} loop over molecules

c*** Sputtering data
c carlis wldspta and wldsptm give sputtered flux for each emitted atom and molecule
            if (ntargsp.gt.0.and.ft44) then
              shift=lppp+1
              do j=1,ntrgts !{
                ISS=ltns(j)
                if (ISS.eq.0) cycle
                IND=eirdiag_nds_ind(ISS)
                IAS=IND+1
                sjs_t(j)=0.0_R8
                sjs(:,j)=0.0_R8
                sptlistrg(j)=0
                do iy=0,ny-1 !{
                  do ia=1,natmi
                    if (wldspta_res(ia,iy+IAS).gt.0.0_R8) then
                      sptlistrg(j)=ia
                      sjs(iy,j)=wldspta_res(ia,iy+IAS)*
     .                    sarea_res(iy+IAS)
                      sjs_t(j)=sjs_t(j)+sjs(iy,j)
                    end if
                  end do
                  do imol=1,nmoli
                    if (wldsptm_res(imol,iy+IAS).gt.0.0_R8) then
                      sptlistrg(j)=imol+natm
                      sjs(iy,j)=wldsptm_res(imol,iy+IAS)*
     .                    sarea_res(iy+IAS)
                      sjs_t(j)=sjs_t(j)+
     .                                 sjs(iy,j)
                    end if
                  end do
                end do

                if((ntns.eq.2 .or. ntns.eq.4) .and.
     .            (sptlistrg(j).ne.0)) then
                  js=ltns(j)+nlim
                  if (sptlistrg(j).le.natm) then
                    sjs_w(j)=wldspta(js,sptlistrg(j),0)
                  else
                    sjs_w(j)=wldsptm(js,sptlistrg(j)-natm,0)
                  end if
                  if(sjs_t(j).gt.0.0_R8) then
                    sjs_t(j)=sjs_w(j)/sjs_t(j)
                  else
                    sjs_t(j)=0.0_R8
                  end if
                else
                  sjs_t(j)=0.0_R8
                end if

                if (sjs_t(j).gt.0.0_R8) then
                  do_spt=.true.
                  do iy=lytg,ny-1!{
                    sjs(iy,j)=sjs_t(j)*sjs(iy,j)
c carlis The emitted amount in g/y is evaluated using the sputtered atom or molecule
c carlis        [g/y]=[1/s]*[s/y]/[1/mol]*[g/mol]
                    if (sptlistrg(j).le.natm) then
                      sptmass=am(eb2atcr(sptlistrg(j)))
                    else
                      imol = sptlistrg(j)-natm
                      sptmass=0.0_R8
                      do ia = 1, natmi
                        if (mlcmp(ia,imol).gt.0) then
                          sptmass=sptmass+mlcmp(ia,imol)*am(eb2atcr(ia))
                        end if
                      end do
                    end if
                    spttrg(iy,j)=sjs(iy,j)*nsecyr/avogr*sptmass
c carlis The erosion of the wall is evaluated using the wall material and not the emitted atom/molecule
c carlis           [mm/y] = [1/s]*[s*mm*m2/y]*
                    erntrg(iy,j)=sjs(iy,j)*
     .                wllersf(i_ersmtr(iy+shift))/sarea_res(iy+IAS)
                  end do
                  if (sptlistrg(j).le.natm) then
                    ssmol=texta(sptlistrg(j))
                  else
                    ssmol=textmn(sptlistrg(j)-natm-1)
                    ssmol=adjustl(ssmol)
                  end if
                  print *,'target material = ', mtrtrgsg(j)
                  print *,'sputtered atom from target = ',
     ,                ssmol
                  write(*,'(1x,a,1p,16e12.3)')
     .             'Scaling factors for '//mtrtrgsg(j)//' Sp : ',
     ,             (sjs_t(i),i=1,ntrgts)
                else ! no sputtering
                  do_spt=.false.
                  print *,
     ,             'b2pwldw, zero sputtering flux for target j=',j
                  sjs_w(j)=0.0
                end if
                shift=shift+ny
              end do !} loop on targets for sputtering
            end if !if sputtering happening
          end if !} itrg.eq.0 .and. do_trg
          itrg=itg
          iznt=izn
          if((lwlddir(izn) .and. (itg.eq.1 .or. itg.eq.3)) .or.
     .     (.not.lwlddir(izn) .and. (itg.eq.2 .or. itg.eq.4))) then !{
            iy= ny*itg - (ils-lppp) + (ny-1)*(itg-1)
            im=-1
          else !}{
            iy= ny+itg-1-(ils-lppp)
            im=1
          end if !}
        else if (iznt.ne.izn) then !}{ itg.ne.itrg
          if((lwlddir(izn) .and. (itg.eq.1 .or. itg.eq.3)) .or.
     .     (.not.lwlddir(izn) .and. (itg.eq.2 .or. itg.eq.4))) then !{
            iy= ny*itg - (ils-lppp) + (ny-1)*(itg-1)
            im=-1
          else !}{
            iy= ny+itg-1-(ils-lppp)
            im=1
          end if !}
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
        write(lun,740) 'Wrad',tab,'power flux with radiation [W/m^2]'
        write(lun,740) 'Wpart',tab,'power flux with particles [W/m^2]'
        write(lun,740) 'Wpls',tab,
     .   'power flux with charged particles [W/m^2]'
        write(lun,740) 'Wneut',tab,
     .   'power flux with neutral particles [W/m^2]'
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
#ifdef B25_EIRENE
        write(lun,740) 'flxm_XX',tab,
     .   'incident flux of molecules XX [1/m^2/s]'
        write(lun,740) 'Eavm_XX',tab,
     .   'average energy of molecules XX [eV]'
        write(lun,740) 'pm_XX',tab,
     .   'density of molecules XX [1/m^3]'
        write(lun,740) 'pm_XX',tab,
     .   'pressure of molecules XX [Pa]'
#endif
        write(lun,740) 'ersne_XX',tab,
     .   'net erosion of wall material XX [mm/y]'
        write(lun,740) 'ersng_XX',tab,
     .   'gross erosion of wall material XX [mm/y]'
        write(lun,740) 'sptne_XX',tab,
     .   'net sputtering of wall material XX [g/y]'
        write(lun,740) 'sptng_XX',tab,
     .   'gross sputtering of wall material XX [g/y]'
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
            if(do_spt_wll(m).eq.is) then !{
c carlis For output sputtering I assume that the sputtered flux is acting on the sputtered atom
c carlis        [g/y]=[1/s]*[s/y]/[1/mol]*[g/mol]
              if (do_spt_wll(m).le.natm) then
                sptmass=am(eb2atcr(do_spt_wll(m)))
              else
                imol = do_spt_wll(m)-natm
                sptmass=0.0_R8
                do ia = 1, natmi
                  if (mlcmp(ia,imol).gt.0) then
                    sptmass=sptmass+mlcmp(ia,imol)*am(eb2atcr(ia))
                  end if
                end do
              end if
              w_line(k+1)=0.0_R8 !net sputtering
              w_line(k+2)=sptflx(m)*nsecyr/avogr*sptmass !gross sputtering
              l=l_ersmtr(ils)
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
              l=l_ersmtr(ils)
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
c        print *,'itrg,iy=',itrg,iy
        w_line(10)=wptpl(iy,itrg) !recombination due to plasma
        w_line(9)=whtpl(iy,itrg) !heat due to plasma                           !###
        w_line(8)=wptnu(iy,itrg)+wptpl(iy,itrg)
        w_line(7)=whtnu(iy,itrg)+whtpl(iy,itrg)
        w_line(6)=wnutrg(iy,itrg)
        w_line(5)=wpltrg(iy,itrg)
        w_line(4)=wpttrg(iy,itrg)
        w_line(3)=wrdtrg(iy,itrg,0)
        w_line(2)=wtttrg(iy,itrg)
        w_line(1)=wldrfdst(iln,izn)
        do is=1,n_sfpr_s !{
          w_line(k+1)=fiitrg(iy,itrg,is)
          w_line(k+2)=eiitrg(iy,itrg,is)/eV
          k=k+n_hdpr_s
        end do !}
        do is=1,nclm_na !{
          w_line(k+1)=faitrg(iy,itrg,is)
          w_line(k+2)=eaitrg(iy,itrg,is)/eV
          w_line(k+3)=naitrg(iy,itrg,is)
          w_line(k+4)=paitrg(iy,itrg,is)
          k=k+n_hdpr_a
        end do !}
        do is=1,nclm_nm !{
          w_line(k+1)=fmitrg(iy,itrg,is) !molecule flux
          w_line(k+2)=emitrg(iy,itrg,is)/eV !molecule energy
          w_line(k+3)=nmitrg(iy,itrg,is) !molecule density
          w_line(k+4)=pmitrg(iy,itrg,is) !molecule pressure
          k=k+n_hdpr_ml
        end do !}
        if (do_spt) then
          do is=1,n_sfpr_sp
            if(sptlistrg(ils).eq.is) then !{
              w_line(k+1)=0.0_R8 !net sputtering
              w_line(k+2)=spttrg(iy,itrg) !gross sputtering
            else
              w_line(k+1)=0.0_R8
              w_line(k+2)=0.0_R8
            end if
            k=k+n_hdpr_sp
          end do
          do is=1,n_sfpr_er
            if(i_ersmtr(ils).eq.is) then !{
              w_line(k+1)=0.0_R8 !net erosion
              w_line(k+2)=erntrg(iy,itrg) !gross erosion
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
          w_line(1)=xtrg(iy,itrg)
          w_line(k+1)=netrg(iy,itrg)
          w_line(k+2)=tetrg(iy,itrg)
          w_line(k+3)=titrg(iy,itrg)
          w_line(k+4)=xmp_hlp(iy,itrg)
          k=k+n_hdpr_t
        end if !}
        iy=iy+im
      else !}{
        k=n_hdpr_c+n_sfpr_s+n_sfpr_sp
        if(trgl) then !{
          k=k+n_hdpr_t
          w_line(1)=xtrg(iy,itrg)
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

c=======================================================================
      end subroutine b2pwldw
