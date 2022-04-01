      module b2mod_b2plot

      use b2mod_types
      use b2mod_layer
     . , only : ntrack
      implicit none
#include <DIMENSIONS.F>

      real (kind=R8), allocatable, save ::
     *  work(:,:,:,:),
     *  wrk(:),
     *  vxa(:,:),
     *  vya(:,:),
     *  vza(:,:),
     *  vxm(:,:),
     *  vym(:,:),
     *  vzm(:,:),
     *  vxi(:,:),
     *  vyi(:,:),
     *  vzi(:,:),
     *  arrowx(:,:,:),
     *  arrowy(:,:,:),
     *  earrowx(:,:),
     *  earrowy(:,:),
     *  p1(:,:),
     *  p2(:,:),
     *  xx(:,:),
     *  yy(:,:),
     *  cr(:,:),
     *  cz(:,:),
     *  nh(:,:),
     *  psi(:,:),
     *  fbz(:,:),
     *  zeff(:,:),
     *  uutraj(:,:,:),
     *  vvtraj(:,:,:),
     *  targetx(:,:),
     *  targety(:,:),
     *  targetr(:,:),
     *  targetz(:,:),
     *  wwrk(:),
     *  wrk_hlp(:),
     *  ewrk_hlp(:),
     *  wwrk_hlp(:),
     *  twrk_hlp(:)

      real (kind=R8), save ::
     *  pp1(2,DEF_NLIM+2*DEF_NYD),
     *  pp2(2,DEF_NLIM+2*DEF_NYD)

      integer, allocatable, save ::
     *  ind_hlp(:), indt_hlp(:)

      real (kind=R8), allocatable, save ::
     *  eirene_stack(:,:,:),
     *  triangle_area(:),
     *  triangle_vol(:),
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
      character*19, allocatable, save ::
     *  textpl(:),
     *  textplel(:),
     *  textplel2(:,:),
     *  textsp(:),
     *  textwk(:,:),
     *  textewk(:,:),
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
      character*19, save ::
     *  textpl(0:DEF_NPLS-1),
     *  textplel(0:DEF_NPLS),
     *  textplel2(0:DEF_NPLS,0:DEF_NPLS),
     *  textsp(0:DEF_NATM-1),
     *  textwl(0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION),
     *  textwk(0:(DEF_NSD+1)*(DEF_NSD+1),wklng),
     *  textewk(0:max(DEF_NATM,DEF_NMOL,DEF_NION)-1,wklng),
     *  textwwk(0:DEF_NSD-1+DEF_NATM+DEF_NMOL+DEF_NION,wklng),
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
     * arrowl,earrowl,arroww,psi_sep,
     * sepwidth,lwidth,wwidth,owidth,vwidth,ewidth,zwidth
      integer, save ::
     * nidim,linlog,wlinlog,tlinlog,trglinlog
      integer, save ::
     *  ndim(3,wklng),nedim(2,wklng),nwdim(2,wklng),
     *  typeofplot,iwrk,iewrk,iwwrk,itwrk,movie_skip,average
cank 960511
      integer, save ::
     *  nlimi,nstsi,natmi,nmoli,nioni,nlim1
cank
      integer, save ::
     *  mcol,scol,vcol,ocol,ecol,acol,zcol,y_transform_l,
     *  dcol,wcol,tcol,hcol,gcol,bcol,xcol,lfsz,expsymbol

      logical, save ::
     *  writing_to_file,color_set,use_compile,
     *  plot_vs_full_psi,plot_vs_norm_psi,plot_vs_root_psi,
     *  set_resignore_regions_to_zero,plotting_wall,map_to_midplane,
     *  xparallel,plotting_sep,expline

      character*256, save ::
     *  graphlabel(wklng),
     *  eirene_label(wklng), wall_label(wklng), target_label(wklng)
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
      external get_jsep, get_nxt, ipgeti, xertst

      nx=nxd
      ny=nyd
      ns=nsd
      nsw=max(nsd+3,(nsd+1)*(nsd+1)-1,ns_ext)

      call get_jsep(nx,ny,jxi,jxa,jsep)
      ixref=jxa
      iyref=jsep
      call ipgeti ('set_transport_ixref', ixref)
      call xertst (-1.le.ixref.and.ixref.le.nxd,
     1   'faulty argument set_transport_ixref')
      call ipgeti ('set_transport_iyref', iyref)
      call xertst (-1.le.iyref.and.iyref.le.nyd,
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
      allocate(psi(-1:nxd,-1:nyd))
      allocate(fbz(-1:nxd,-1:nyd))
      allocate(zeff(-1:nxd,-1:nyd))
      allocate(uutraj(-1:nxd,-1:nyd,0:nsd-1))
      allocate(vvtraj(-1:nxd,-1:nyd,0:nsd-1))
#ifndef PSR
      allocate(textpl(0:nsd-1))
      allocate(textplel(0:nsd))
      allocate(textplel2(0:nsd,0:nsd))
      allocate(textsp(0:max(1,natmd)-1))
      allocate(textwl(0:nsd-1+natmd+nmold+niond))
      allocate(textwk(0:nsw,wklngd))
      allocate(textewk(0:max(1,natmd,nmold,niond)-1,wklngd))
      allocate(textwwk(0:nsd-1+natmd+nmold+niond,wklngd))
      allocate(textbl(0:0))
      allocate(textni(0:1))
      allocate(textcomponents(0:3))
      allocate(text0123(0:3))
      allocate(textcorner(0:3))
      allocate(texttrack(0:ntrack))
      allocate(texttargsp(ntrack))
#else
      call xertst(natmd .le.DEF_NATM,
     . 'increase size of DEF_NATM in DIMENSIONS.F')
      call xertst(nmold .le.DEF_NMOL,
     . 'increase size of DEF_NMOL in DIMENSIONS.F')
      call xertst(niond .le.DEF_NION,
     . 'increase size of DEF_NION in DIMENSIONS.F')
      call xertst(nsd   .le.DEF_NPLS,
     . 'increase size of DEF_NPLS in DIMENSIONS.F')
      call xertst(wklngd.le.wklng,
     . 'increase size of textwk array')
#endif

      return
      end subroutine alloc_b2mod_b2plot

      subroutine alloc_b2mod_b2plot_eirene(natmd,nmold,niond,
     .  ntrii,wklngd)
      implicit none
      integer natmd,nmold,niond,ntrii,wklngd

      allocate(eirene_stack(ntrii,0:max(1,natmd,nmold,niond)-1,wklngd))
      allocate(triangle_area(ntrii))
      allocate(triangle_vol(ntrii))
      allocate(vxa(max(1,natmd),ntrii))
      allocate(vya(max(1,natmd),ntrii))
      allocate(vza(max(1,natmd),ntrii))
      allocate(vxm(max(1,nmold),ntrii))
      allocate(vym(max(1,nmold),ntrii))
      allocate(vzm(max(1,nmold),ntrii))
      allocate(vxi(max(1,niond),ntrii))
      allocate(vyi(max(1,niond),ntrii))
      allocate(vzi(max(1,niond),ntrii))
      allocate(earrowx(ntrii,0:max(1,natmd,nmold,niond)-1))
      allocate(earrowy(ntrii,0:max(1,natmd,nmold,niond)-1))
      allocate(ewrk_hlp(ntrii))

      return
      end subroutine alloc_b2mod_b2plot_eirene

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

      subroutine dealloc_b2mod_b2plot
      implicit none

      if (.not.allocated(work)) return
      deallocate(work)
      deallocate(wrk)
      deallocate(wrk_hlp)
      deallocate(ind_hlp)
      deallocate(arrowx)
      deallocate(arrowy)
      deallocate(p1)
      deallocate(p2)
      deallocate(xx)
      deallocate(yy)     
      deallocate(cr)
      deallocate(cz)
      deallocate(nh)
      deallocate(psi)
      deallocate(fbz)
      deallocate(zeff)
      deallocate(uutraj)
      deallocate(vvtraj)
#ifndef PSR
      deallocate(textpl)
      deallocate(textplel)
      deallocate(textplel2)
      deallocate(textsp)
      deallocate(textwl)
      deallocate(textwk)
      deallocate(textewk)
      deallocate(textwwk)
      deallocate(textbl)
      deallocate(textni)
      deallocate(textcomponents)
      deallocate(text0123)
      deallocate(textcorner)
      deallocate(texttrack)
      deallocate(texttargsp)
#endif

      return
      end subroutine dealloc_b2mod_b2plot

      subroutine dealloc_b2mod_b2plot_eirene
      implicit none

      deallocate(eirene_stack)
      deallocate(triangle_area)
      deallocate(triangle_vol)
      deallocate(vxa)
      deallocate(vya)
      deallocate(vza)
      deallocate(vxm)
      deallocate(vym)
      deallocate(vzm)
      deallocate(vxi)
      deallocate(vyi)
      deallocate(vzi)
      deallocate(earrowx)
      deallocate(earrowy)
      deallocate(ewrk_hlp)

      return
      end subroutine dealloc_b2mod_b2plot_eirene

      subroutine dealloc_b2mod_b2plot_wall
      implicit none

      deallocate(wwrk_hlp)
      deallocate(twrk_hlp)
      deallocate(indt_hlp)
      deallocate(wall_stack)
      deallocate(wall_flux)
      deallocate(wall_area)
      deallocate(incident_energy)
      deallocate(wall_heatflux)
      deallocate(wall_temp)
      deallocate(targetx)
      deallocate(targety)
      deallocate(targetr)
      deallocate(targetz)
      deallocate(wwrk)
      deallocate(target_stack)

      return
      end subroutine dealloc_b2mod_b2plot_wall

      subroutine read_b2mod_b2plot(nread)
      implicit none
      integer nread
      external xertst

      call xertst(0.lt.nread, 'faulty input nread')

      return
      end subroutine read_b2mod_b2plot

      subroutine write_b2mod_b2plot(nwrite)
      implicit none
      integer nwrite
      external xertst

      call xertst(0.lt.nwrite, 'faulty input nwrite')

      return
      end subroutine write_b2mod_b2plot

      function ix_e2b(ixe)
      use b2mod_indirect
      implicit none
      integer ix_e2b
      integer ixe
      integer ic, i

      if (ixe.eq.-1 .or. ixe.eq.0) then
        ix_e2b = -2
      else
        i = 1
        ix_e2b = 0
        do while (i.lt.ixe)
          i = i + 1
          ix_e2b = ix_e2b + 1
          do ic = 1, nncut
            if (ix_e2b.eq.leftcut(ic)) i = i + 1
            if (ix_e2b.eq.rightcut(ic)) i = i + 1
          end do
          if (nncut.eq.2 .and. ix_e2b.eq.nxtr) i = i + 1
        end do
      end if
      return
      end function ix_e2b

      function ix_b2e(ixb)
      use b2mod_indirect
      implicit none
      integer ix_b2e
      integer ixb
      integer i

      if (ixb.eq.-2) then
        ix_b2e = -1
      else if (ixb.gt.nx) then
        ix_b2e = -1
      else if (ixb.eq.-1 .or. ixb.eq.nx) then
        ix_b2e = -1
      else if (nncut.eq.2 .and.
     .        (ixb.ge.nxtl .and. ixb.le.nxtr)) then
        ix_b2e = -1
      else
        ix_b2e = ixb + 1
        do i = 1, nncut
          if (ixb.ge.leftcut(i)) ix_b2e = ix_b2e + 1
          if (ixb.ge.rightcut(i)) ix_b2e = ix_b2e + 1
        end do
        if (nncut.eq.2 .and. ixb.gt.nxtr) ix_b2e = ix_b2e + 1
      end if
      return
      end function ix_b2e

      end module b2mod_b2plot
