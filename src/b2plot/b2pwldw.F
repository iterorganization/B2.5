      subroutine b2pwldw(lun,trgl,izn,iln,boundary_namelist)
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_rates
      use b2mod_diag
      use b2mod_constants
      use b2mod_neutrals_namelist
      use b2mod_boundary_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_extrab25
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
c***  Net erosion             |   for each kind of sputtered material
c***  Gross erosion           |
c***
c*** Extra output for targets:
c***  ne
c***  Te
c***  Ti
c***  Mid-plane co-ordinate
c***
c*** Apart from this, some extra output is written to the common-block 
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

      integer lun,izn,iln,boundary_namelist
      logical trgl

c*** Local variables
c***  ils     : the segment index in the whole list
c***  itrg    : target number itrg prepared
c***  izns    : plot zone izns started
c***  zz      : reference distance
c***  i_sptmtr: index of sputtered material for the wall segments

      integer i,j,k,l,m,n,i1,i2,i3,ib,is,iy,im,jj,js,kk,ils,itg,jcf,jci,
     ,  itrg,izns,nclm_ld,nclm_er,nclm_pk,nx_hlp,m_hlp(ntrgsx),
     ,  i_sptmtr(nsgmx),l_sptmtr(0:nsgmx)

c*** Arrays to store the target profiles
c***  xtrg    : co-ordinate along the target
c***  xmp_hlp : the same projected to the mid-plane
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
c***  wrdtrg  : radiation power to the target
c***  wtttrg  : total power to the target
c***  fiitrg  : impinging ion flux 
c***  eiitrg  : average energy of ions
c***  faitrg  : impinging atom flux
c***  eaitrg  : average energy of atoms
c***  erntrg  : erosion net
c***  ergtrg  : erosion gross

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

      real (kind=R8) :: h,p,t,u,wo,zz,wwo,eni_hlp(ntrgsx),
     ,  phi_hlp(ntrgsx),wi(ntrgsx),wwi(ntrgsx),z_hlp(ntrgsx),
     ,  s(0:ny-1,ntrgsx),xmp_hlp(0:ny-1,ntrgsx),xtrg(0:ny-1,ntrgsx),  
     ,  netrg(0:ny-1,ntrgsx) , tetrg(0:ny-1,ntrgsx),
     ,  titrg(0:ny-1,ntrgsx), 
     ,  wtttrg(0:ny-1,ntrgsx), wptnu(0:ny-1,ntrgsx),
     ,   wptpl(0:ny-1,ntrgsx),
     ,  whtnu(0:ny-1,ntrgsx) , whtpl(0:ny-1,ntrgsx),
     ,  wpttrg(0:ny-1,ntrgsx),
     ,  wpltrg(0:ny-1,ntrgsx),wnutrg(0:ny-1,ntrgsx),
     ,  wrdtrg(0:ny-1,ntrgsx),
     ,  fiitrg(0:ny-1,ntrgsx,DEF_NATM), eiitrg(0:ny-1,ntrgsx,DEF_NATM),
     ,  faitrg(0:ny-1,ntrgsx,DEF_NATM), eaitrg(0:ny-1,ntrgsx,DEF_NATM),
     ,  erntrg(0:ny-1,ntrgsx,nsptmx)  , ergtrg(0:ny-1,ntrgsx,nsptmx),
     ,  wpk_tt(0:ntrgsx),
     ,  wpk_rd(0:ntrgsx),wpk_pl(0:ntrgsx),wpk_nu(0:ntrgsx),
     ,  wrd_pk(0:ntrgsx),wpl_pk(0:ntrgsx),wnu_pk(0:ntrgsx),
     ,  wlc_tt(ntrgsx),wlc_rd(ntrgsx),wlc_pl(ntrgsx),wlc_nu(ntrgsx)

c*** Output data headers
c***  m_hdpr    : max. length of a quantity name
c***  m_hdln    : max. length of the header line
c***  n_hdpr_c  : number of common quantities
c***  n_hdpr_t  : number of target quantities
c***  n_hdpr_s  : number of quantities by species
c***  n_hdpr_m  : number of quantities by material
c***  n_sfpr_s  : number of active species
c***  n_sfpr_m  : number of different sputtered materials for erosion
c***  hdpr_c    : headers for common quantities
c***  hdpr_t    : headers for target quantities
c***  hdpr_s    : headers for the quantities by species
c***  hdpr_m    : headers for the quantities by material
c***  sfpr_s    : list of active species
c***  sfpr_m    : list of active materials

      integer n_hdpr_c,n_hdpr_t,n_hdpr_s,n_hdpr_m,n_sfpr_s,n_sfpr_m,
     ,  m_hdpr,m_hdln,n_line_x
      parameter (n_hdpr_c=8, n_hdpr_t=4, n_hdpr_s=4, n_hdpr_m=2,
     ,  m_hdpr=5, 
     ,  n_line_x=n_hdpr_c+n_hdpr_t+(n_hdpr_s+n_hdpr_m)*DEF_NATM,
     ,  m_hdln=(n_hdpr_c+n_hdpr_t)*(m_hdpr+1)+
     +         (n_hdpr_s+n_hdpr_m)*DEF_NATM*(m_hdpr+4))
      real (kind=R8) :: w_line(n_line_x)
      real (kind=R8) :: s_hlp
      logical :: first,do_trg
      character*(m_hdpr) tab*1,hdpr_c(n_hdpr_c),hdpr_t(n_hdpr_t),
     ,  hdpr_s(n_hdpr_s),hdpr_m(n_hdpr_m),sfpr_s(DEF_NATM)*2,
     ,  sfpr_m(nsptmx)*(lnwlmtrl),hdr_line*(m_hdln),frmt_hlp*24

      data hdpr_c /'x', 'Wtot','Wrad','Wpart','Wpls',
     2  'Wneu','Wheat','Wpot' /
      data hdpr_t /'ne', 'Te','Ti','xMP'/
      data hdpr_s /'flxi','Eavi','flxa','Eava'/
      data hdpr_m /'ersne','ersng'/
      data itrg, izns, n_sfpr_s /3*0/
      data first /.true./
      data sfpr_s, sfpr_m / DEF_NATM*' ', nsptmx*' ' /

c*** Data for erosion calculations
c***  n_srfs_m  : number of possible wall materials
c***  srfs_m    : list of possible wall materials
c***  srfs_x    : corresponding species (say, Fe for SS)
c***  srfs_d    : specific masses of possible wall materials [g/cm**3]
c***  srfs_a    : atomic weight of possible wall materials
c***  wllersf   : ratio of erosion to sputtering for different meterials
      integer n_srfs_m
      parameter (n_srfs_m=7)
      real (kind=R8) :: srfs_d(n_srfs_m),srfs_a(n_srfs_m),
     .                  wllersf(nsptmx)
      character srfs_m(n_srfs_m)*4,srfs_x(n_srfs_m)*2
      logical streql
      external streql
      data srfs_m /'Be','C','CFC','Fe','SS','Mo','W'/
      data srfs_x /'Be','C','C'  ,'Fe','Fe','Mo','W'/
      data srfs_a /  9.,12., 12. , 56., 56., 96.,184./
      data srfs_d /1.8 ,2.0, 2.0 , 7.8, 7.8,10.2,19.3/
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
        print *,'Starting b2pwldw first time',lun,trgl,izn,iln  !###
        tab=char(9)
        do_trg=.true.
        write(frmt_hlp,'(7h(1p,1x,,i5,11h(e12.5,a1)))') n_line_x
        
        n_sfpr_s=1
        if(nspecies.gt.1) n_sfpr_s=nspecies
        if(n_sfpr_s.ne.natmi) then !{
          write(0,*) '!!! Inconsistency found in b2pwldw: numbers',
     ,      ' of atomic and plasma species differ: ',natmi,n_sfpr_s
          write(0,*) '==> No particle fluxes in the output'
          n_sfpr_s=0
        end if !}
        do is=1,n_sfpr_s !{
          sfpr_s(is)=atoms_list(is)
        end do !}
        print *,'First call. format:',frmt_hlp,':'                  !###
        print *,'n_sfpr_s=',n_sfpr_s                                !###
        if(n_sfpr_s.gt.0) print *,(' ',sfpr_s(i),i=1,n_sfpr_s)      !###
        k=0
        do i=1,nbc !{
          l=bcpos(i)
          if(bcchar(i).eq.'W') then !{
            k=k+1
            iy=(bcstart(i)+bcend(i))/2
            if(l.eq.-1) then !{
              phi_hlp(1)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(1)=enipar(i,1)
            else !}{
              phi_hlp(3)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(3)=enipar(i,1)
            end if !}
          else if(bcchar(i).eq.'E') then !}{
            k=k+1
            iy=(bcstart(i)+bcend(i))/2
            if(l.eq.nx) then !{
              if(ntrgts.eq.2) then !{
                phi_hlp(2)=get_phi_applied(boundary_namelist,l,iy)
                eni_hlp(2)=enipar(i,1)
              else !}{
                phi_hlp(4)=get_phi_applied(boundary_namelist,l,iy)
                eni_hlp(4)=enipar(i,1)
              end if !}
            else if(ntrgts.eq.4) then !}{
              phi_hlp(2)=get_phi_applied(boundary_namelist,l,iy)
              eni_hlp(2)=enipar(i,1)
            end if !}
          end if !}
        end do !}
        print '(1x,a,1p,4e12.3)','phi_hlp:',(phi_hlp(i),i=1,ntrgts) !###
        print '(1x,a,1p,4e12.3)','eni_hlp:',(eni_hlp(i),i=1,ntrgts) !###
        if(k.ne.ntrgts) then !{
          write(0,*) 'b2pwldw: inconsistent target specification: ',
     ,                                                         ntrgts,k
          do_trg=.false.
        end if !}

c*** data for sputtering and erosion

        print *,'Input wall materials'                              !###
        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=1,lppp)            !###
        print *,'Input target materials'                            !###
        print '(1h:,20(a,1h:))',(mtrtrgsg(i),i=1,ntrgts)            !###
        k=lppp
        do j=1,ntrgts !{
          do i=k+1,k+ny !{
            mtrwldsg(i)=mtrtrgsg(j)
          end do !}
          k=k+ny
        end do !}
        print *,'k,lppp,ny,ntrgts=',k,lppp,ny,ntrgts                !###
        print *,'Target segment materials'                          !###
        print '((1h:,20(a,1h:)))',(mtrwldsg(i),i=lppp+1,k)          !###

        n_sfpr_m=0
        do k=1,lppp+ny*ntrgts !{
          l=0
          do i=1,n_sfpr_m !{
            if(streql(sfpr_m(i),mtrwldsg(k))) l=i
          end do !}
          if(l.eq.0) then !{
            if(n_sfpr_m.lt.nsptmx) then !{
              n_sfpr_m=n_sfpr_m+1
              l=n_sfpr_m
              sfpr_m(l)=mtrwldsg(k)
              print *,'Added :',sfpr_m(l),':',k,l                   !###
            else !}{
            end if !}
          end if !}
          i_sptmtr(k)=l
        end do !}
        print *,'Input materials: ',n_sfpr_m                        !###
        if(n_sfpr_m.gt.0) print *,':',(sfpr_m(i),':',i=1,n_sfpr_m)  !###
        l_sptmtr(0)=0
        m=0
        do k=1,n_sfpr_m !{
          l=0
          do i=1,n_srfs_m !{
            if(streql(srfs_m(i),sfpr_m(k))) l=i
          end do !}
          if(l.eq.0) then !{
            write(*,*) 'b2pwldw: unknown wall material ',sfpr_m(k)
            l_sptmtr(k)=0
            
          else !}{
            m=m+1
            sfpr_m(m)=srfs_x(l)
            wllersf(m)=3.15e7*1.e-3*srfs_a(l)/(avogr*srfs_d(l))
            l_sptmtr(k)=m
          end if !}
        end do !}
        n_sfpr_m=m
        do i=1,lppp+ny*ntrgts !{
          i_sptmtr(i)=l_sptmtr(i_sptmtr(i))
        end do !}
        write(*,*) 'Sputtered materials: ',n_sfpr_m
        if(n_sfpr_m.gt.0) then !{
          write(*,*) 'Eroded species : ',(sfpr_m(i),i=1,n_sfpr_m)
          write(*,'(a,1p,10e10.2)') ' Erosion factors: ',
     ,                                        (wllersf(i),i=1,n_sfpr_m)
        end if !}
        print *,'i_sptmtr'                                    !###
        print '(1x,50i2)',(i_sptmtr(i),i=1,lppp+ny*ntrgts)    !###
        first=.false.
      end if !}

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
c        print *,'itg,itrg=',itg,itrg  !###
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
                if (mod(j,2).eq.1) then ! West targets
                  z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                                     (cry(k,iy,3)-cry(k,iy,1))**2)
                else ! East targets
                  z_hlp(j)=z_hlp(j) - sqrt((crx(k,iy,2)-crx(k,iy,0))**2+
     +                                     (cry(k,iy,2)-cry(k,iy,0))**2)
                endif
              end do !}
              do iy=0,ny-1 !{
                if (mod(j,2).eq.1) then ! West targets
                  h=sqrt((crx(k,iy,3)-crx(k,iy,1))**2+
     +                   (cry(k,iy,3)-cry(k,iy,1))**2)
                  s(iy,j)=gs(k,iy,0)
                else ! East targets
                  h=sqrt((crx(k,iy,2)-crx(k,iy,0))**2+
     +                   (cry(k,iy,2)-cry(k,iy,0))**2)
                  s(iy,j)=gs(k,iy+1,0)
                endif
                z_hlp(j)=z_hlp(j)+0.5_R8*h
                xtrg(iy,j)=z_hlp(j)
                z_hlp(j)=z_hlp(j)+0.5_R8*h
              end do !}
            end do !}

            if(trgl) then !{

c*** Mid-plane co-ordinates

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
              k=ltrgpos(j)
              kk=k+l
              do iy=0,ny-1 !{
                wo=0.0_R8
                do is=0,ns-1 !{
                  wo=wo+rpt(k,iy,is)*ev*fna(k,iy,0,is)
                end do !}
                wptpl(iy,j)=wo/s(iy,j)*u
                whtpl(iy,j)=(fhe(k,iy,0)+fhi(k,iy,0))/s(iy,j)*u
                wpltrg(iy,j)=wptpl(iy,j)+whtpl(iy,j)
                netrg(iy,j)=ne(kk,iy)
                tetrg(iy,j)=te(kk,iy)/ev
                titrg(iy,j)=ti(kk,iy)/ev
                wnutrg(iy,j)=0.0_R8
                do is=1,natmi
                  t=tab2(kk+1,iy+1,is,1)/ev
                  wnutrg(iy,j)=wnutrg(iy,j)+s(iy,j)*
     *                       dab2(kk+1,iy+1,is,1)*tab2(kk+1,iy+1,is,1)*
     *                       sqrt(max(t/am(eb2spcr(is)),0.0_R8))*
     *                       0.25_R8*2.0e4_R8
                end do
              end do !}
              l=1-l
              u=-u
            end do !}

c*** Re-scale the neutral target loading and calculate the sums

            do j=1,ntrgts !{
              wi(j)=0.0_R8
              do i=0,ny-1 !{
                wi(j)=wi(j)+wnutrg(i,j)
              end do !}
            end do !}
c            print *,'ntrgts, wi: ',ntrgts,(wi(i),i=1,ntrgts)  !###
            do j=1,ntrgts !{
              if(ntns.eq.ntrgts) then !{
                js=ltns(j)+nlim
                wo=wldnek(js,0)+wldnep(js,0)
c                print *,'j,js,wo = ',j,js,wo    !###
                if(wi(j).gt.0.0_R8) then !{
                  wi(j)=wo/wi(j)
                else !}{
                  wi(j)=0.0_R8
                end if !}
              else !}{
                write(*,*) 'b2pwldw: ntns.ne.ntrgts (',
     ,            ntns,' vs.',ntrgts,
     ,            ' ) - no neutral load onto the targets calculated'
                write(*,*) '  ==> check specification of ltns in ',
     ,            'b2.parameters (section NEUTRALS)'
                wi(j)=0.0_R8
              end if !}

              do iy=0,ny-1 !{
                wnutrg(iy,j)=wnutrg(iy,j)*wi(j)/s(iy,j)
                if(abs(wo).gt.0.0_R8) then !{
                  wptnu(iy,j)=wnutrg(iy,j)*wldnep(js,0)/wo
                  whtnu(iy,j)=wnutrg(iy,j)-wptnu(iy,j)
                else !}{
                  wptnu(iy,j)=0.0_R8
                  whtnu(iy,j)=0.0_R8
                end if !}
                wpttrg(iy,j)=wnutrg(iy,j)+wpltrg(iy,j)
                wtttrg(iy,j)=wpttrg(iy,j)+wrdtrg(iy,j)
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
                k=ltrgpos(j)
                kk=k+l
                ib=k+1-l
c!###{
c                print *,'is,j,i1,i2,k,kk,ib=',is,j,i1,i2,k,kk,ib
c                zz=0
c                print '(3x,12(2x,a8))','phi','fti','tn','wo','wwo',
c     ,            'fiitrg','eiitrg','faitrg','eaitrg','dab2',
c     ,            'tab2','s'
c!###}
                do iy=0,ny-1 !{
                  wo=0.0_R8
                  wwo=0.0_R8
                  p=(po(kk,iy)-phi_hlp(j))*qe
                  h=ti(kk,iy)*eni_hlp(j)
                  do i=i1,i2 !{
                    wo=wo+fna(k,iy,0,i)
                    wwo=wwo+(h+p*rza(kk,iy,i))*fna(k,iy,0,i)
                  end do !}
c                  zz=zz+wwo*u                 !###
                  fiitrg(iy,j,is)=u*wo
                  if(abs(wo).gt.0.0_R8) then !{
                    eiitrg(iy,j,is)=wwo/wo/ev
                  else !}{
                    eiitrg(iy,j,is)=0.0_R8
                  end if !}
                  t=tab2(kk+1,iy+1,is,1)/ev
                  faitrg(iy,j,is)=s(iy,2)*dab2(kk+1,iy+1,is,1)*
     *             sqrt(max(t/am(eb2spcr(is)),0.))*0.25_R8*1.0e4_R8
                  eaitrg(iy,j,is)=s(iy,2)*dab2(kk+1,iy+1,is,1)*
     *             sqrt(max(t/am(eb2spcr(is)),0.0_R8))*
     *             tab2(kk+1,iy+1,is,1)*0.25_R8*2.0e4_R8
c!###{
c                  print '(1x,i2,1p,12e10.2)',iy,p/ev,h/ev,t,wo,wwo,
c     ,              fiitrg(iy,j,is),eiitrg(iy,j,is),
c     ,              faitrg(iy,j,is),eaitrg(iy,j,is),
c     ,              dab2(kk+1,iy+1,is,1),tab2(kk+1,iy+1,is,1)/ev,s(iy,2)
c!###}
                end do !}
c                print *,'j,is,zz=',j,is,zz  !###
                wo=0.0_R8
                wwo=0.0_R8
                do iy=0,ny-1 !{
                  wo=wo+faitrg(iy,j,is)
                  wwo=wwo+eaitrg(iy,j,is)
                end do !}
                js=ltns(j)+nlim
c                print *,'j,is,ltns,js=',j,is,ltns(j),js  !###
c                print *,'wldna,ewlda=',wldna(js,is,0),ewlda(js,is,0)!###
c                print *,'wo,wwo=',wo,wwo                   !###
                if(wo.gt.0.0_R8) then !{
                  wi(j)=wldna(js,is,0)/wo
                else !}{
                  wi(j)=0.0_R8
                end if !}
                if(wwo.gt.0.0_R8) then !{
                  wwi(j)=ewlda(js,is,0)*wldna(js,is,0)*eV/wwo
                else !}{
                  wwi(j)=0.0_R8
                end if !}
                do iy=0,ny-1 !{
                  faitrg(iy,j,is)=faitrg(iy,j,is)*wi(j)
                  if(faitrg(iy,j,is).gt.0.0_R8) then !{
                    eaitrg(iy,j,is)=eaitrg(iy,j,is)*wwi(j)/ev/
     /                                                  faitrg(iy,j,is)
                  else !}{
                  end if !}
                end do !}
                l=1-l
                u=-u
              end do !}
              write(*,*) 'Species',is
              write(*,*) 'Flux scaling'
              write(*,'(1p,12e10.2)') (wi(i),i=1,ntrgts)
              write(*,*) 'Energy scaling'
              write(*,'(1p,12e10.2)') (wwi(i),i=1,ntrgts)
            end do !}

c*** Sputtering data: zeroes for the time being

            do j=1,ntrgts !{
              do iy=0,ny-1 !{
                do is=1,n_sfpr_m !{
                  erntrg(iy,j,is)=0.0_R8
                  ergtrg(iy,j,is)=0.0_R8
                end do !}
              end do !}
            end do !}
          end if !}
          itrg=itg
          if((lwlddir(izn) .and. (itg.eq.1 .or. itg.eq.3)) .or.
     .      (.not.lwlddir(izn) .and. (itg.eq.2 .or. itg.eq.4))) then !{
            iy=ny-1
            im=-1
          else !}{
            iy=0
            im=1
          end if !}
        end if !}
      end if !}

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
        do j=1,n_hdpr_m !{
          l=len_trim(hdpr_m(j))
          do i=1,n_sfpr_m !{
            m=len_trim(sfpr_m(i))
            hdr_line=hdr_line(1:k)//hdpr_m(j)(1:l)//'_'//
     /                              sfpr_m(i)(1:m)//tab
            k=k+l+m+2
          end do !}
        end do !}
c        print *,'About to write a header:',k-1,' symbols'   !###
c        print *,hdr_line(1:k-1)                             !###
        write(lun,'(a)') hdr_line(1:k-1)
      end if !}

c*** Print the line

      m=lelmind(ils)
      k=n_hdpr_c
      if(itg.eq.0) then !{
        w_line(8)=wldnep(m,0)/segm_area(ils)
        w_line(7)=wldnek(m,0)/segm_area(ils)
        w_line(6)=w_line(7)+w_line(8)
        w_line(5)=0
        w_line(4)=w_line(6)
        w_line(3)=wldrdld(ils)
        w_line(2)=w_line(3)+w_line(4)
        w_line(1)=wldrfdst(iln,izn)
        do is=1,n_sfpr_s !{
          w_line(k+1)=0
          w_line(k+2)=0
          w_line(k+3)=wldna(m,is,0)/segm_area(ils)
          w_line(k+4)=ewlda(m,is,0)
          k=k+n_hdpr_s
        end do !}
        do is=1,n_sfpr_m !{
          if(i_sptmtr(ils).eq.is) then !{
            w_line(k+2)=wldspt(m,0)*wllersf(is)/segm_area(ils)
            w_line(k+1)=0
            l=l_sptmtr(ils)
            if(l.gt.0) then !{
              wo=-wldna(m,l,0)+wldra(m,l,0)+wldpa(m,l,0)
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
              w_line(k+1)=wo*wllersf(is)/segm_area(ils)
            end if !}
          else !}{
            w_line(k+2)=0
            w_line(k+1)=0
          end if !}
          k=k+n_hdpr_m
        end do !}
      else if(do_trg) then !}{
c        print *,'itrg,iy=',itrg,iy                                !###
        w_line(8)=wptnu(iy,itrg)+wptpl(iy,itrg)
        w_line(7)=whtnu(iy,itrg)+whtpl(iy,itrg)
        w_line(6)=wnutrg(iy,itrg)
        w_line(5)=wpltrg(iy,itrg)
        w_line(4)=wpttrg(iy,itrg)
        w_line(3)=wrdtrg(iy,itrg)
        w_line(2)=wtttrg(iy,itrg)
        w_line(1)=wldrfdst(iln,izn)
        do is=1,n_sfpr_s !{
          w_line(k+1)=fiitrg(iy,itrg,is)
          w_line(k+2)=eiitrg(iy,itrg,is)
          w_line(k+3)=faitrg(iy,itrg,is)
          w_line(k+4)=eaitrg(iy,itrg,is)
          k=k+n_hdpr_s
        end do !}
        do is=1,n_sfpr_m !{
          w_line(k+2)=ergtrg(iy,itrg,is)
          w_line(k+1)=erntrg(iy,itrg,is)
          k=k+n_hdpr_m
        end do !}
        if(trgl) then !{
          w_line(1)=xtrg(iy,itrg)
          w_line(k+1)=netrg(iy,itrg)
          w_line(k+2)=tetrg(iy,itrg)
          w_line(k+3)=titrg(iy,itrg)
          w_line(k+4)=xmp_hlp(iy,itrg)
          k=k+n_hdpr_t
        end if !}
      else !}{
        k=n_hdpr_c+n_sfpr_s+n_sfpr_m
        if(trgl) then !{
          k=k+n_hdpr_t
          w_line(1)=xtrg(iy,itrg)
        else !}{
          w_line(1)=wldrfdst(iln,izn)
        end if !}
        do i=2,k !{
          w_line(i)=0
        end do !}
      end if !}
c      print *,'About to write a line:',k,' values'   !###

      if(k.gt.0) write(lun,frmt_hlp) (w_line(i),tab,i=1,k)
      
c=======================================================================
      end
