      module b2mod_b2plot

      use b2mod_types
      use b2mod_layer
     . , only : ntrack
      implicit none
#include <DIMENSIONS.F>

      real (kind=R8), allocatable, save ::
     *  work(:,:,:,:),
     *  wrk(:),
     *  arrowx(:,:,:),
     *  arrowy(:,:,:),
     *  p1(:,:),
     *  p2(:,:),
     *  xx(:,:),
     *  yy(:,:),
     *  cr(:,:),
     *  cz(:,:),
     *  nh(:,:),
     *  uutraj(:,:,:),
     *  vvtraj(:,:,:),
     *  targetx(:,:),
     *  targety(:,:),
     *  targetr(:,:),
     *  targetz(:,:),
     *  wwrk(:),
     *  wrk_hlp(:),
     *  wwrk_hlp(:),
     *  twrk_hlp(:)

      real (kind=R8), save ::
     *  pp1(2,DEF_NLIM+2*DEF_NYD),
     *  pp2(2,DEF_NLIM+2*DEF_NYD)

      integer, allocatable, save ::
     *  ind_hlp(:), indt_hlp(:)

      real (kind=R8), allocatable, save ::
     *  wall_stack(:,:,:),
     *  wall_flux(:,:),
     *  wall_area(:),
     *  incident_energy(:,:),
     *  wall_heatflux(:),
     *  wall_temp(:),
     *  target_stack(:,:,:)

      integer wklng
      parameter (wklng=10)

#ifndef PSR
      character*16, allocatable, save ::
     *  textan(:),
     *  textmn(:),
     *  textin(:),
     *  textpl(:),
     *  textplel(:),
     *  textplel2(:,:),
     *  textsp(:),
     *  textwk(:,:),
     *  textwwk(:,:),
     *  textbl(:),
     *  textni(:),
     *  textcomponents(:),
     *  text0123(:),
     *  textwl(:),
     *  textcorner(:),
     *  texttrack(:),
     *  texttargsp(:)    
#else
      character*16, save ::
     *  textan(0:DEF_NATM-1),
     *  textmn(0:DEF_NMOL-1),
     *  textin(0:DEF_NION-1),
     *  textpl(0:DEF_NPLS-1),
     *  textplel(0:DEF_NPLS),
     *  textplel2(0:DEF_NPLS,0:DEF_NPLS),
     *  textsp(0:DEF_NATM-1),
     *  textwl(0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  textwk(0:(DEF_NSD+1)*(DEF_NSD+1),wklng),
     *  textbl(0:0),
     *  textni(0:1),
     *  textcomponents(0:3),
     *  text0123(0:3),
     *  textcorner(0:3),
     *  texttrack(0:ntrack),
     *  texttargsp(ntrack)
#endif

      integer npart
c typeofplot=1 physical domain only, =2 computational domain only, =3 both
      parameter(npart=100000)

      real (kind=R8), save ::
     * arrowl,arroww,sepwidth,lwidth,wwidth,owidth,vwidth,ewidth,zwidth
      integer, save ::
     * nidim,linlog,wlinlog,tlinlog,trglinlog
      integer, save ::
     *  ndim(3,wklng),nwdim(2,wklng),
     *  typeofplot,iwrk,iwwrk,itwrk,movie_skip,average
cank 960511
      integer, save ::
     *  nlimi,nstsi,natmi,nmoli,nioni
cank
      integer, save ::
     *  mcol,scol,vcol,ocol,ecol,acol,zcol,y_transform_l,
     *  dcol,wcol,tcol,hcol,gcol,bcol,xcol,lfsz,expsymbol

      logical, save ::
     *  writing_to_file,color_set,use_compile,
     *  set_resignore_regions_to_zero,plotting_wall,map_to_midplane,
     *  xparallel,expline

      character*256, save ::
     *  graphlabel(wklng), wall_label(wklng), target_label(wklng)
      character*256, save ::
     *  solpstop,moviename,expfile,wallfile
      character*80, save ::
     *  b2version

      integer, save :: nx, ny, ns, nsw, jxa, jxi, jsep, ixref, iyref,
     *  nxtr, nxtl

      integer, save :: neo_ic

      contains

      subroutine alloc_b2mod_b2plot(nxd,nyd,nsd,ns_ext,wklngd,
     1 natmd,nmold,niond,nlim)
      implicit none
      integer nxd,nyd,nsd,ns_ext,wklngd,natmd,nmold,niond,nlim
      integer nsw
      external ipgeti, xertst, get_jsep, get_nxt

      nx=nxd
      ny=nyd
      ns=nsd
      nsw=max(nsd+3,(nsd+1)*(nsd+1)-1,ns_ext)

      call get_jsep(nx,ny,jxi,jxa,jsep)
      ixref=jxa
      iyref=-2
      call ipgeti ('set_transport_ixref', ixref)
      call xertst (-1.le.ixref.and.ixref.le.nxd,
     1   'faulty argument set_transport_ixref')
      call ipgeti ('set_transport_iyref', iyref)
      call xertst (-2.le.iyref.and.iyref.le.nyd,
     1   'faulty argument set_transport_iyref')
      call get_nxt(nx,nxtl,nxtr)

      allocate(work(-1:nxd,-1:nyd,0:nsw,1:wklngd))
      allocate(wrk(-1:nxd+nyd+1))
      allocate(wrk_hlp(1:nxd*nyd))
      allocate(ind_hlp(1:nxd*nyd))
      allocate(arrowx(-1:nxd,-1:nyd,0:nsw))
      allocate(arrowy(-1:nxd,-1:nyd,0:nsw))
      allocate(p1(2,nlim))
      allocate(p2(2,nlim))
      allocate(xx(-1:nx,2))
      allocate(yy(-1:ny,2))     
      allocate(cr(-1:nxd,-1:nyd))
      allocate(cz(-1:nxd,-1:nyd))
      allocate(nh(-1:nxd,-1:nyd))
      allocate(uutraj(-1:nxd,-1:nyd,0:nsd-1))
      allocate(vvtraj(-1:nxd,-1:nyd,0:nsd-1))
#ifndef PSR
      allocate(textan(0:natmd-1))
      allocate(textmn(0:nmold-1))
      allocate(textin(0:niond-1))
      allocate(textpl(0:nsd-1))
      allocate(textplel(0:nsd))
      allocate(textplel2(0:nsd,0:nsd))
      allocate(textsp(0:natmd-1))
      allocate(textwl(0:nsd-1+natmd+nmold+niond))
      allocate(textwk(0:nsw,wklngd))
      allocate(textbl(0:0))
      allocate(textni(0:1))
      allocate(textcomponents(0:3))
      allocate(text0123(0:3))
      allocate(textcorner(0:3))
      allocate(texttrack(0:ntrack))
      allocate(texttargsp(ntrack))
#else
      call xertst(natmd .le.DEF_NATM,'increase size of textam array')
      call xertst(nmold .le.DEF_NMOL,'increase size of textml array')
      call xertst(niond .le.DEF_NION,'increase size of textin array')
      call xertst(nsd   .le.DEF_NPLS,'increase size of textpl array')
      call xertst(wklngd.le.wklng   ,'increase size of textwk array')
#endif
      allocate(textwwk(0:nsd-1+natmd+nmold+niond,wklngd))

      return
      end subroutine alloc_b2mod_b2plot

      subroutine alloc_b2mod_b2plot_wall(nsd,natmd,nmold,niond,
     .  nwall,ndepth,wklngd)
      implicit none
      integer nsd,natmd,nmold,niond,nwall,ndepth,wklngd

      allocate(wwrk_hlp(1:nwall))
      allocate(twrk_hlp(1:nwall*ndepth))
      allocate(indt_hlp(1:nwall*ndepth))
      allocate(wall_stack(1:nwall,0:nsd-1+natmd+nmold+niond,1:wklngd))
      allocate(wall_flux(1:nwall,0:nsd-1+natmd+nmold+niond))
      allocate(wall_area(1:nwall))
      allocate(incident_energy(1:nwall,0:nsd-1+natmd+nmold+niond))
      allocate(wall_heatflux(1:nwall))
      allocate(wall_temp(1:nwall))
      allocate(targetx(1:nwall,2))
      allocate(targety(1:ndepth,2))
      allocate(targetr(1:nwall,1:ndepth))
      allocate(targetz(1:nwall,1:ndepth))
      allocate(wwrk(1:nwall+ndepth))
      allocate(target_stack(1:nwall,1:ndepth,wklngd))

      return
      end subroutine alloc_b2mod_b2plot_wall

      subroutine read_b2mod_b2plot(nread)
      implicit none
      integer nread
      call xertst(0.lt.nread, 'faulty input nread')

      return
      end subroutine read_b2mod_b2plot

      subroutine write_b2mod_b2plot(nwrite)
      implicit none
      integer nwrite
      call xertst(0.lt.nwrite, 'faulty input nwrite')

      return
      end subroutine write_b2mod_b2plot

      end module b2mod_b2plot
