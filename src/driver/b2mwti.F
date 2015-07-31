*-----------------------------------------------------------------------
*.specification

#ifdef IBMaix
@process noextchk
#endif
      subroutine b2mwti (itim, tim, nx, ny, ns, ismain, ismain0, BoRiS)
      use b2mod_types
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_residuals
      use b2mod_sources
      use b2mod_transport
      use b2mod_anomalous_transport
      use b2mod_neutrals_namelist
      use b2mod_work
      use b2mod_indirect
      use b2mod_constants
      use b2mod_tallies
      use b2mod_wall
      use b2mod_b2cmpa
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include <COUPLE/EIRDIAG.F>
#endif
*   ..input arguments (unchanged on exit)
      integer itim, nx, ny, ns, ismain, ismain0
      real (kind=R8) ::
     *  tim, BoRiS
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
#ifndef NO_CDF
#     include <netcdf.inc>
#endif
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MWTI writes to an un*formatted file a selection of data
*     describing the state of the plasma.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) data that are suitable for
*     producing movie output. Unit nput must be connected to an
*     un*formatted file. Presently, only the changes in the plasma
*     state are written out.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ntstep
      real (kind=R8) ::
     *      fnixip(nncutmax), feexip(nncutmax), feixip(nncutmax),
     *      fnixap(nncutmax), feexap(nncutmax), feixap(nncutmax),
     *      pwmxip(nncutmax), fchxip(nncutmax), fchxap(nncutmax),
     *      fetxip(nncutmax), fetxap(nncutmax)
      real (kind=R8) ::
     *      fniyip(nncutmax), feeyip(nncutmax), feiyip(nncutmax),
     *      fniyap(nncutmax), feeyap(nncutmax), feiyap(nncutmax),
     *      pwmxap(nncutmax), fchyip(nncutmax), fchyap(nncutmax),
     *      fetyip(nncutmax), fetyap(nncutmax)
      real (kind=R8) ::
     *      nemxip(nncutmax), temxip(nncutmax), timxip(nncutmax),
     *      nemxap(nncutmax), temxap(nncutmax), timxap(nncutmax),
     *      pomxip(nncutmax), pomxap(nncutmax), tpmxip(nncutmax),
     *      tpmxap(nncutmax)
      real (kind=R8) ::
     *      fnisip(nncutmax), feesip(nncutmax), feisip(nncutmax),
     *      fnisap(nncutmax), feesap(nncutmax), feisap(nncutmax),
     *      fchsip(nncutmax), fchsap(nncutmax),
     *      fetsip(nncutmax), fetsap(nncutmax)
      real (kind=R8) ::
     *      fnisipp(nncutmax), feesipp(nncutmax), feisipp(nncutmax),
     *      fnisapp(nncutmax), feesapp(nncutmax), feisapp(nncutmax),
     *      fchsipp(nncutmax), fchsapp(nncutmax),
     *      fetsipp(nncutmax), fetsapp(nncutmax)
      real (kind=R8) ::
     *      tmne(1),tmte(1),tmti(1),tmvol,kintmp,rpttmp,fettmp

      integer iy, ix, ic, is, ixtl, ixtr, jsep
      integer jxi, jxa, target_offset
      integer iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt,
     *        iyaend,  iyiend,  iylend,  iyrend,  iytlend,  iytrend, 
     *        nybl, nybr, nytl, nytr, nya, nyi, nc
*   ..procedures
      external subini, subend, xertst, ipgeti, output_ds
#ifndef NO_CDF
      integer imap(maxvdims), iret, ncid
      real (kind=R8) ::
     *      nesepi(nncutmax), tesepi(nncutmax), tisepi(nncutmax),
     *      nesepm(nncutmax), tesepm(nncutmax), tisepm(nncutmax),
     *      nesepa(nncutmax), tesepa(nncutmax), tisepa(nncutmax),
     *      posepi(nncutmax), posepm(nncutmax), posepa(nncutmax),
     1      dnsepm(nncutmax), dpsepm(nncutmax), kesepm(nncutmax),
     2      kisepm(nncutmax), vxsepm(nncutmax), vysepm(nncutmax),
     3      vssepm(nncutmax), tpsepi(nncutmax), tpsepa(nncutmax)
      real (kind=R8) ::
     *      tmhacore(1), tmhasol(1), tmhadiv(1), slice(-1:ny), tstepn(1)
      real (kind=R8) ::
     *      timesa(1)
      character*5 rw
      external b2crtimecdf, rwcdf, rwcdf_settime, ncclos
#endif
      intrinsic max, sqrt
*   ..initialisation
      save ncall, ntstep, jxi, jxa, jsep, ixtl, ixtr, target_offset,
     * iyastrt, iyistrt, iylstrt, iyrstrt, iytlstrt, iytrstrt,
     * iyaend,  iyiend,  iylend,  iyrend,  iytlend,  iytrend,
     * nybl, nybr, nytl, nytr, nya, nyi, nc
      data ncall/0/,target_offset/1/

      real (kind=R8) ::
     *   hy1
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mwti')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns, 
     . 'invalid main plasma species index ismain')
*   ..extensive tests on first few calls
      if (ncall.eq.0) then
*    ..test state
        call get_jsep(nx,ny,jxi,jxa,jsep)
        call ipgeti ('b2mwti_target_offset',target_offset)
        call xertst (0.le.target_offset.and.target_offset.le.1,
     2    'faulty internal parameter target_offset')
        write(*,*) 'target_offset ', target_offset
        call output_ds(crx,cry,nx,ny, -1,+target_offset,
     1   jsep,iylstrt,iylend,'../dsl')
        call output_ds(crx,cry,nx,ny,jxi,0,
     1   jsep,iyistrt,iyiend,'../dsi')
        call output_ds(crx,cry,nx,ny,jxa,0,
     1   jsep,iyastrt,iyaend,'../dsa')
        call output_ds(crx,cry,nx,ny, nx,-target_offset,
     1   jsep,iyrstrt,iyrend,'../dsr')
        if (nnreg(0).eq.8) then
          ixtl = 0
          do while (rightix(ixtl,max(topcut(1),topcut(2))).ne.nx+1
     &     .and.ixtl.lt.nx)
            ixtl=ixtl+1
          enddo
          ixtr = ixtl
          do while (leftix(ixtr,max(topcut(1),topcut(2))).ne.-2 .and.
     &              ixtr.lt.nx)
            ixtr=ixtr+1
          enddo
          call output_ds(crx,cry,nx,ny,ixtl,-target_offset,
     1     jsep,iytlstrt,iytlend,'../dstl')
          call output_ds(crx,cry,nx,ny,ixtr,+target_offset,
     1     jsep,iytrstrt,iytrend,'../dstr')
          nytl = iytlend - iytlstrt + 1
          nytr = iytrend - iytrstrt + 1
        else
          nytl = 0
          nytr = 0
        endif
        nybl = iylend  - iylstrt  + 1
        nybr = iyrend  - iyrstrt  + 1
        nya  = iyaend  - iyastrt  + 1
        nyi  = iyiend  - iyistrt  + 1
        nc = max(nncut,1)
c Target areas
        open(99,file='../dsL')
        do iy=-1,ny
          if(region(-1,iy,0).ne.0)
     >     write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)
        enddo
        close(99)
        open(99,file='../dsR')
        do iy=-1,ny
          if(region(nx,iy,0).ne.0)
     >     write(99,*) gs(nx,iy,0)
        enddo
        close(99)
        if (nnreg(0).eq.8) then
          open(99,file='../dsTL')
          do iy=iytlstrt,iytlend
            write(99,*) gs(ixtl,iy,0)
          enddo
          close(99)
          open(99,file='../dsTR')
          do iy=iytrstrt,iytrend
            write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          enddo
          close(99)
        endif
c Poloidal contact areas
        open(99,file='../dsLP')
        do iy=-1,ny
          if(region(-1,iy,0).ne.0)
     >     write(99,*) gs(rightix(-1,iy),rightiy(-1,iy),0)
     $          *qc(rightix(-1,iy),rightiy(-1,iy))
        enddo
        close(99)
        open(99,file='../dsRP')
        do iy=-1,ny
          if(region(nx,iy,0).ne.0)
     >     write(99,*) gs(nx,iy,0)*qc(nx,iy)
        enddo
        close(99)
        if (nnreg(0).eq.8) then
          open(99,file='../dsTLP')
          do iy=iytlstrt,iytlend
            write(99,*) gs(ixtl,iy,0)*qc(ixtl,iy)
          enddo
          close(99)
          open(99,file='../dsTRP')
          do iy=iytrstrt,iytrend
            write(99,*) gs(rightix(ixtr,iy),rightiy(ixtr,iy),0)
     $            *qc(rightix(ixtr,iy),rightiy(ixtr,iy))
          enddo
          close(99)
        endif
#ifndef NO_CDF
        call b2crtimecdf(nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc,
     1   ns, iret)
        rw='write'
        ncid=ncopn('b2time.nc',ncwrite,iret)
c       write(*,*) 'After ncopn: ', ncid, ncwrite, iret
        ntstep=0
        imap(1)=1
        tstepn(1) = ntstep
        call rwcdf (rw, ncid, 'ntstep', imap, tstepn, iret)
        call ncclos (ncid, iret)
#endif
      endif
      if(ncall.lt.3) then
        call xertst (0.le.itim, 'faulty parameter itim')
        call xertst (0.0_R8.le.tim, 'faulty parameter tim')
      endif
*   ..compute change in plasma state
c
      ntstep = ntstep + 1
c     write(*,*) 'ntstep = ',ntstep
#ifndef NO_CDF
      rw = 'write'
      ncid = ncopn ('b2time.nc', ncwrite, iret)
      call rwcdf_settime ('time', ntstep)
      timesa(1) = tim
#endif
c
c    total flows to the divertor plates
c
      fnixip = 0.0_R8
      feexip = 0.0_R8
      feixip = 0.0_R8
      fetxip = 0.0_R8
      fchxip = 0.0_R8
      fnixap = 0.0_R8
      feexap = 0.0_R8
      feixap = 0.0_R8
      fetxap = 0.0_R8
      fchxap = 0.0_R8
      nemxip = 0.0_R8
      temxip = 0.0_R8
      timxip = 0.0_R8
      tpmxip = 0.0_R8
      pomxip = 0.0_R8
      nemxap = 0.0_R8
      temxap = 0.0_R8
      timxap = 0.0_R8
      tpmxap = 0.0_R8
      pomxap = 0.0_R8
      pwmxip = 0.0_R8
      pwmxap = 0.0_R8
      fnisip = 0.0_R8
      feesip = 0.0_R8
      feisip = 0.0_R8
      fetsip = 0.0_R8
      fchsip = 0.0_R8
      fnisap = 0.0_R8
      feesap = 0.0_R8
      feisap = 0.0_R8
      fetsap = 0.0_R8
      fchsap = 0.0_R8
      fnisipp = 0.0_R8
      feesipp = 0.0_R8
      feisipp = 0.0_R8
      fetsipp = 0.0_R8
      fchsipp = 0.0_R8
      fnisapp = 0.0_R8
      feesapp = 0.0_R8
      feisapp = 0.0_R8
      fetsapp = 0.0_R8
      fchsapp = 0.0_R8
c
      do iy=iylstrt,iylend
        fnixip(1) = fnixip(1) + 
     1   fna(rightix(-1,iy),rightiy(-1,iy),0,ismain)
        feexip(1) = feexip(1) + fhe(rightix(-1,iy),rightiy(-1,iy),0)
        feixip(1) = feixip(1) + fhi(rightix(-1,iy),rightiy(-1,iy),0)
        fchxip(1) = fchxip(1) + fch(rightix(-1,iy),rightiy(-1,iy),0)
        nemxip(1) = max (nemxip(1), ne(-1+target_offset,iy))
        temxip(1) = max (temxip(1), te(-1+target_offset,iy))
        timxip(1) = max (timxip(1), ti(-1+target_offset,iy))
        pomxip(1) = max (pomxip(1), po(-1+target_offset,iy))
        fettmp = fhe(rightix(-1,iy),rightiy(-1,iy),0) +
     2           fhi(rightix(-1,iy),rightiy(-1,iy),0)
        do is=0,ns-1
          kintmp = (kinrgy(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     .              kinrgy(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     .             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
          rpttmp = (rpt(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     .              rpt(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     .             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
          fettmp = fettmp + (rpttmp*ev + kintmp*(1.0_R8-BoRiS))*
     .     fna(rightix(-1,iy),rightiy(-1,iy),0,is)
        enddo
        fetxip(1) = fettmp
        pwmxip(1) = max (pwmxip(1),
     1   -fettmp/gs(rightix(-1,iy),rightiy(-1,iy),0))
        if (bottomiy(-1,iy).ne.-2 .and. topiy(-1,iy).ne.ny+1)
     >   tpmxip(1) = max (tpmxip(1), target_temp(xymap(-1,iy),1))
      enddo
      do iy=iyrstrt,iyrend
        fnixap(1) = fnixap(1) + fna(nx,iy,0,ismain)
        feexap(1) = feexap(1) + fhe(nx,iy,0)
        feixap(1) = feixap(1) + fhi(nx,iy,0)
        fchxap(1) = fchxap(1) + fch(nx,iy,0)
        nemxap(1) = max (nemxap(1), ne(nx-target_offset,iy))
        temxap(1) = max (temxap(1), te(nx-target_offset,iy))
        timxap(1) = max (timxap(1), ti(nx-target_offset,iy))
        pomxap(1) = max (pomxap(1), po(nx-target_offset,iy))
        fettmp = fhe(nx,iy,0) + fhi(nx,iy,0)
        do is=0,ns-1
          kintmp = (kinrgy(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     .              kinrgy(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     .             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
          rpttmp = (rpt(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     .              rpt(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     .             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
          fettmp = fettmp + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))
     .     *fna(nx,iy,0,is)
        enddo
        fetxap(1) = fettmp
        pwmxap(1) = max (pwmxap(1), fettmp/gs(nx,iy,0))
        if (bottomiy(nx,iy).ne.-2 .and. topiy(nx,iy).ne.ny+1)
     >   tpmxap(1) = max (tpmxap(1), target_temp(xymap(nx,iy),1))
      enddo
      if(nncut.ge.2) then
        do iy = iytrstrt, iytrend
          fnixip(2) = fnixip(2) + 
     .     fna(rightix(ixtr,iy),rightiy(ixtr,iy),0,ismain)
          feexip(2) = feexip(2) +  
     .     fhe(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          feixip(2) = feixip(2) +  
     .     fhi(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          fchxip(2) = fchxip(2) +  
     .     fch(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          nemxip(2) = max (nemxip(2), ne(ixtr+target_offset,iy))
          temxip(2) = max (temxip(2), te(ixtr+target_offset,iy))
          timxip(2) = max (timxip(2), ti(ixtr+target_offset,iy))
          pomxip(2) = max (pomxip(2), po(ixtr+target_offset,iy))
          fettmp = fhe(rightix(ixtr,iy),rightiy(ixtr,iy),0) +
     .             fhi(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ixtr,iy,is)*hx(rightix(ixtr,iy),rightiy(ixtr,iy))+
     .        kinrgy(rightix(ixtr,iy),rightiy(ixtr,iy),is)*hx(ixtr,iy))/
     .       (hx(rightix(ixtr,iy),rightiy(ixtr,iy))+hx(ixtr,iy))
            rpttmp = 
     .       (rpt(ixtr,iy,is)*hx(rightix(ixtr,iy),rightiy(ixtr,iy))+
     .        rpt(rightix(ixtr,iy),rightiy(ixtr,iy),is)*hx(ixtr,iy))/
     .       (hx(rightix(ixtr,iy),rightiy(ixtr,iy))+hx(ixtr,iy))
            fettmp = fettmp + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .       fna(rightix(ixtr,iy),rightiy(ixtr,iy),0,is)
          enddo
          fetxip(2) = fettmp
          pwmxip(2) = max (pwmxip(2),
     1     -fettmp/gs(rightix(ixtr,iy),rightiy(ixtr,iy),0))
          if (bottomiy(ixtr,iy).ne.-2 .and. topiy(ixtr,iy).ne.ny+1)
     >     tpmxip(2) = max (tpmxip(2), target_temp(xymap(ixtr,iy),1))
        enddo
        do iy = iytlstrt, iytlend
          fnixap(2) = fnixap(2) + fna(ixtl,iy,0,ismain)
          feexap(2) = feexap(2) + fhe(ixtl,iy,0)
          feixap(2) = feixap(2) + fhi(ixtl,iy,0)
          fchxap(2) = fchxap(2) + fch(ixtl,iy,0)
          nemxap(2) = max (nemxap(2), ne(ixtl-target_offset,iy))
          temxap(2) = max (temxap(2), te(ixtl-target_offset,iy))
          timxap(2) = max (timxap(2), ti(ixtl-target_offset,iy))
          pomxap(2) = max (pomxap(2), po(ixtl-target_offset,iy))
          fettmp = fhe(ixtl,iy,0) + fhi(ixtl,iy,0)
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ixtl,iy,is)*hx(leftix(ixtl,iy),leftiy(ixtl,iy))+
     .        kinrgy(leftix(ixtl,iy),leftiy(ixtl,iy),is)*hx(ixtl,iy))/
     .       (hx(leftix(ixtl,iy),leftiy(ixtl,iy))+hx(ixtl,iy))
            rpttmp = 
     .       (rpt(ixtl,iy,is)*hx(leftix(ixtl,iy),leftiy(ixtl,iy))+
     .        rpt(leftix(ixtl,iy),leftiy(ixtl,iy),is)*hx(ixtl,iy))/
     .       (hx(leftix(ixtl,iy),leftiy(ixtl,iy))+hx(ixtl,iy))
            fettmp = fettmp + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .       fna(ixtl,iy,0,is)
          enddo
          fetxap(2) = fettmp
          pwmxap(2) = max (pwmxap(2), fettmp/gs(ixtl,iy,0))
          if (bottomiy(ixtl,iy).ne.-2 .and. topiy(ixtl,iy).ne.ny+1)
     >     tpmxap(2) = max (tpmxap(2), target_temp(xymap(ixtl,iy),1))
        enddo
      endif
      if(nnreg(0).ge.3) then
        do iy=-1,jsep
          do ic = 1, nncut
            fnisipp(ic) = fnisipp(ic) + fna(leftcut(ic),iy,0,ismain)
            feesipp(ic) = feesipp(ic) + fhe(leftcut(ic),iy,0)
            feisipp(ic) = feisipp(ic) + fhi(leftcut(ic),iy,0)
            fetsipp(ic) = fetsipp(ic) + fhe(leftcut(ic),iy,0)
     .                                + fhi(leftcut(ic),iy,0)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(leftix(leftcut(ic),iy),
     .                 leftiy(leftcut(ic),iy),is)*hx(leftcut(ic),iy)+
     .          kinrgy(leftcut(ic),iy,is)*
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))/
     .         (hx(leftcut(ic),iy)+
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))
              rpttmp = 
     .         (rpt(leftix(leftcut(ic),iy),
     .              leftiy(leftcut(ic),iy),is)*hx(leftcut(ic),iy)+
     .          rpt(leftcut(ic),iy,is)*
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))/
     .         (hx(leftcut(ic),iy)+
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))
              fetsipp(ic) = fetsipp(ic) + 
     .         (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(leftcut(ic),iy,0,is)
            enddo
            fchsipp(ic) = fchsipp(ic) + fch(leftcut(ic),iy,0)
            fnisapp(ic) = fnisapp(ic) + fna(rightcut(ic),iy,0,ismain)
            feesapp(ic) = feesapp(ic) + fhe(rightcut(ic),iy,0)
            feisapp(ic) = feisapp(ic) + fhi(rightcut(ic),iy,0)
            fetsapp(ic) = fetsapp(ic) + fhe(rightcut(ic),iy,0)
     .                                + fhi(rightcut(ic),iy,0)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(leftix(rightcut(ic),iy),
     .                 leftiy(rightcut(ic),iy),is)*hx(rightcut(ic),iy)+
     .          kinrgy(rightcut(ic),iy,is)*
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))/
     .         (hx(rightcut(ic),iy)+
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))
              rpttmp = 
     .         (rpt(leftix(rightcut(ic),iy),
     .              leftiy(rightcut(ic),iy),is)*hx(rightcut(ic),iy)+
     .          rpt(rightcut(ic),iy,is)*
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))/
     .         (hx(rightcut(ic),iy)+
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))
              fetsapp(ic) = fetsapp(ic) + 
     .         (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(rightcut(ic),iy,0,is)
            enddo
            fchsapp(ic) = fchsapp(ic) + fch(rightcut(ic),iy,0)
          enddo
        enddo
        do iy=jsep+1,ny
          do ic = 1, nncut
            fnisip(ic) = fnisip(ic) + fna(leftcut(ic),iy,0,ismain)
            feesip(ic) = feesip(ic) + fhe(leftcut(ic),iy,0)
            feisip(ic) = feisip(ic) + fhi(leftcut(ic),iy,0)
            fetsip(ic) = fetsip(ic) + fhe(leftcut(ic),iy,0)
     .                              + fhi(leftcut(ic),iy,0)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(leftix(leftcut(ic),iy),
     .                 leftiy(leftcut(ic),iy),is)*hx(leftcut(ic),iy)+
     .          kinrgy(leftcut(ic),iy,is)*
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))/
     .         (hx(leftcut(ic),iy)+
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))
              rpttmp = 
     .         (rpt(leftix(leftcut(ic),iy),
     .              leftiy(leftcut(ic),iy),is)*hx(leftcut(ic),iy)+
     .          rpt(leftcut(ic),iy,is)*
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))/
     .         (hx(leftcut(ic),iy)+
     .          hx(leftix(leftcut(ic),iy),leftiy(leftcut(ic),iy)))
              fetsip(ic) = fetsip(ic) + 
     .         (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(leftcut(ic),iy,0,is)
            enddo
            fchsip(ic) = fchsip(ic) + fch(leftcut(ic),iy,0)
            fnisap(ic) = fnisap(ic) + fna(rightcut(ic),iy,0,ismain)
            feesap(ic) = feesap(ic) + fhe(rightcut(ic),iy,0)
            feisap(ic) = feisap(ic) + fhi(rightcut(ic),iy,0)
            fetsap(ic) = fetsap(ic) + fhe(rightcut(ic),iy,0)
     .                              + fhi(rightcut(ic),iy,0)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(leftix(rightcut(ic),iy),
     .                 leftiy(rightcut(ic),iy),is)*hx(rightcut(ic),iy)+
     .          kinrgy(rightcut(ic),iy,is)*
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))/
     .         (hx(rightcut(ic),iy)+
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))
              rpttmp = 
     .         (rpt(leftix(rightcut(ic),iy),
     .              leftiy(rightcut(ic),iy),is)*hx(rightcut(ic),iy)+
     .          rpt(rightcut(ic),iy,is)*
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))/
     .         (hx(rightcut(ic),iy)+
     .          hx(leftix(rightcut(ic),iy),leftiy(rightcut(ic),iy)))
              fetsap(ic) = fetsap(ic) + 
     .         (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(rightcut(ic),iy,0,is)
            enddo
            fchsap(ic) = fchsap(ic) + fch(rightcut(ic),iy,0)
          enddo
        enddo
      endif
      
      fniyip = 0.0_R8
      feeyip = 0.0_R8
      feiyip = 0.0_R8
      fetyip = 0.0_R8
      fchyip = 0.0_R8
      fniyap = 0.0_R8
      feeyap = 0.0_R8
      feiyap = 0.0_R8
      fetyap = 0.0_R8
      fchyap = 0.0_R8
      
      if(nnreg(0).eq.4) then
        do ix=-1,nx
          if(region(ix,ny,0).eq.2) then
            fniyip(1) = fniyip(1) + fna(ix,ny,1,ismain)
            feeyip(1) = feeyip(1) + fhe(ix,ny,1)
            feiyip(1) = feiyip(1) + fhi(ix,ny,1)
            fchyip(1) = fchyip(1) + fch(ix,ny,1)
            fetyip(1) = fetyip(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyip(1)=fetyip(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,ny,0).eq.3.or.region(ix,ny,0).eq.4) then
            fniyap(1) = fniyap(1) + fna(ix,ny,1,ismain)
            feeyap(1) = feeyap(1) + fhe(ix,ny,1)
            feiyap(1) = feiyap(1) + fhi(ix,ny,1)
            fchyap(1) = fchyap(1) + fch(ix,ny,1)
            fetyap(1) = fetyap(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyap(1)=fetyap(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
            fniyap(1) = fniyap(1) - 
     &       fna(topix(ix,-1),topiy(ix,-1),1,ismain)
            feeyap(1) = feeyap(1) - fhe(topix(ix,-1),topiy(ix,-1),1)
            feiyap(1) = feiyap(1) - fhi(topix(ix,-1),topiy(ix,-1),1)
            fchyap(1) = fchyap(1) - fch(topix(ix,-1),topiy(ix,-1),1)
            fetyap(1) = fetyap(1) -
     1       fhe(topix(ix,-1),topiy(ix,-1),1) -
     2       fhi(topix(ix,-1),topiy(ix,-1),1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          kinrgy(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              rpttmp = 
     .         (rpt(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          rpt(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              fetyap(1)=fetyap(1) - (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     1         fna(topix(ix,-1),topiy(ix,-1),1,is)
            enddo
          endif
        enddo
      elseif(nnreg(0).eq.5) then
        do ix=-1,nx
          if(region(ix,ny,0).eq.5) then
            fniyip(1) = fniyip(1) + fna(ix,ny,1,ismain)
            feeyip(1) = feeyip(1) + fhe(ix,ny,1)
            feiyip(1) = feiyip(1) + fhi(ix,ny,1)
            fchyip(1) = fchyip(1) + fch(ix,ny,1)
            fetyip(1) = fetyip(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyip(1)=fetyip(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,-1,0).eq.3.or.region(ix,-1,0).eq.4) then
            fniyap(1) = fniyap(1) - 
     &       fna(topix(ix,-1),topiy(ix,-1),1,ismain)
            feeyap(1) = feeyap(1) - fhe(topix(ix,-1),topiy(ix,-1),1)
            feiyap(1) = feiyap(1) - fhi(topix(ix,-1),topiy(ix,-1),1)
            fchyap(1) = fchyap(1) - fch(topix(ix,-1),topiy(ix,-1),1)
            fetyap(1) = fetyap(1) -
     1       fhe(topix(ix,-1),topiy(ix,-1),1) -
     2       fhi(topix(ix,-1),topiy(ix,-1),1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          kinrgy(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              rpttmp = 
     .         (rpt(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          rpt(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              fetyap(1)=fetyap(1) - (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     1         fna(topix(ix,-1),topiy(ix,-1),1,is)
            enddo
          endif
        enddo
      elseif (nnreg(0).eq.8) then
        do ix=-1,nx
          if(mod(region(ix,ny,0),4).eq.2) then
            fniyip(region(ix,ny,0)/4+1) = fniyip(region(ix,ny,0)/4+1) + 
     &       fna(ix,ny,1,ismain)
            feeyip(region(ix,ny,0)/4+1) = feeyip(region(ix,ny,0)/4+1) +
     &       fhe(ix,ny,1)
            feiyip(region(ix,ny,0)/4+1) = feiyip(region(ix,ny,0)/4+1) +
     &       fhi(ix,ny,1)
            fchyip(region(ix,ny,0)/4+1) = fchyip(region(ix,ny,0)/4+1) +
     &       fch(ix,ny,1)
            fetyip(region(ix,ny,0)/4+1) = fetyip(region(ix,ny,0)/4+1) +
     1       fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyip(region(ix,ny,0)/4+1) = fetyip(region(ix,ny,0)/4+1)+
     1         (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,ny,0).eq.3 .or. region(ix,ny,0).eq.8) then
            fniyap(1) = fniyap(1) + fna(ix,ny,1,ismain)
            feeyap(1) = feeyap(1) + fhe(ix,ny,1)
            feiyap(1) = feiyap(1) + fhi(ix,ny,1)
            fchyap(1) = fchyap(1) + fch(ix,ny,1)
            fetyap(1) = fetyap(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyap(1)=fetyap(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,-1,0).eq.3 .or. region(ix,-1,0).eq.8) then
            fniyap(1) = fniyap(1) - 
     &       fna(topix(ix,-1),topiy(ix,-1),1,ismain)
            feeyap(1) = feeyap(1) - fhe(topix(ix,-1),topiy(ix,-1),1)
            feiyap(1) = feiyap(1) - fhi(topix(ix,-1),topiy(ix,-1),1)
            fchyap(1) = fchyap(1) - fch(topix(ix,-1),topiy(ix,-1),1)
            fetyap(1) = fetyap(1) -
     1       fhe(topix(ix,-1),topiy(ix,-1),1) -
     2       fhi(topix(ix,-1),topiy(ix,-1),1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          kinrgy(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              rpttmp = 
     .         (rpt(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          rpt(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              fetyap(1)=fetyap(1) - (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     1         fna(topix(ix,-1),topiy(ix,-1),1,is)
            enddo
          endif
          if(region(ix,ny,0).eq.4 .or. region(ix,ny,0).eq.7) then
            fniyap(2) = fniyap(2) + fna(ix,ny,1,ismain)
            feeyap(2) = feeyap(2) + fhe(ix,ny,1)
            feiyap(2) = feiyap(2) + fhi(ix,ny,1)
            fchyap(2) = fchyap(2) + fch(ix,ny,1)
            fetyap(2) = fetyap(2) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyap(2)=fetyap(2) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
          if(region(ix,-1,0).eq.4 .or. region(ix,-1,0).eq.7) then
            fniyap(2) = fniyap(2) - 
     &       fna(topix(ix,-1),topiy(ix,-1),1,ismain)
            feeyap(2) = feeyap(2) - fhe(topix(ix,-1),topiy(ix,-1),1)
            feiyap(2) = feiyap(2) - fhi(topix(ix,-1),topiy(ix,-1),1)
            fchyap(2) = fchyap(2) - fch(topix(ix,-1),topiy(ix,-1),1)
            fetyap(2) = fetyap(2) -
     1       fhe(topix(ix,-1),topiy(ix,-1),1) -
     2       fhi(topix(ix,-1),topiy(ix,-1),1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          kinrgy(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              rpttmp = 
     .         (rpt(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     .          rpt(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     .         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
              fetyap(2)=fetyap(2) - (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     1         fna(topix(ix,-1),topiy(ix,-1),1,is)
            enddo
          endif
        enddo
      else if (nnreg(0).eq.2) then
        do ix=-1,nx
          if(region(ix,ny,0).eq.2) then
            fniyip(1) = fniyip(1) + fna(ix,ny,1,ismain)
            feeyip(1) = feeyip(1) + fhe(ix,ny,1)
            feiyip(1) = feiyip(1) + fhi(ix,ny,1)
            fchyip(1) = fchyip(1) + fch(ix,ny,1)
            fetyip(1) = fetyip(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
            do is=0,ns-1
              kintmp = 
     .         (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              rpttmp = 
     .         (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .          rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .         (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
              fetyip(1)=fetyip(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .         fna(ix,ny,1,is)
            enddo
          endif
        enddo
      else
        do ix=-1,nx
          fniyip(1) = fniyip(1) + fna(ix,ny,1,ismain)
          feeyip(1) = feeyip(1) + fhe(ix,ny,1)
          feiyip(1) = feiyip(1) + fhi(ix,ny,1)
          fchyip(1) = fchyip(1) + fch(ix,ny,1)
          fetyip(1) = fetyip(1) + fhe(ix,ny,1) + fhi(ix,ny,1)
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .        kinrgy(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
            rpttmp = 
     .       (rpt(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     .        rpt(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     .       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
            fetyip(1)=fetyip(1) + (rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .       fna(ix,ny,1,is)
          enddo
        enddo
      endif
c
c    other quantities related to the target plates
c
#ifndef NO_CDF
      if(nnreg(0).ne.2) then
        nesepi(1) = 0.5_R8 *
     1  (ne(-1+target_offset,jsep)+
     2   ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
        tesepi(1) = 0.5_R8/ev *
     1   (te(-1+target_offset,jsep)+
     2    te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
        tisepi(1) = 0.5_R8/ev *
     1   (ti(-1+target_offset,jsep)+
     2    ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
        tpsepi(1) = 0.5_R8 *
     1   (target_temp(xymap(-1,jsep),1)+
     2    target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1))
        posepi(1) = 0.5_R8 *
     1   (po(-1+target_offset,jsep)+
     2    po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep)))
      else
        nesepi(1) = 
     1   ne(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
        tesepi(1) = 1.0_R8/ev *
     1   te(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
        tisepi(1) = 1.0_R8/ev *
     1   ti(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
        tpsepi(1) = 
     1   target_temp(xymap(topix(-1,jsep),topiy(-1,jsep)),1)
        posepi(1) = 
     1   po(topix(-1+target_offset,jsep),topiy(-1+target_offset,jsep))
      endif
      nesepm(1) = 0.5_R8 * (ne(jxa,jsep)+
     1 ne(topix(jxa,jsep),topiy(jxa,jsep)))
      tesepm(1) = 0.5_R8 * (te(jxa,jsep)+
     1 te(topix(jxa,jsep),topiy(jxa,jsep)))/ev
      tisepm(1) = 0.5_R8 * (ti(jxa,jsep)+
     1 ti(topix(jxa,jsep),topiy(jxa,jsep)))/ev
      posepm(1) = 0.5_R8 * (po(jxa,jsep)+
     1 po(topix(jxa,jsep),topiy(jxa,jsep)))
      dnsepm(1) = 0.5_R8 * (dna0(jxa,jsep,ismain)+
     1 dna0(topix(jxa,jsep),topiy(jxa,jsep),ismain))
      dpsepm(1) = 0.5_R8 * (
     1 dpa0(jxa,jsep,ismain0)*(
     2 rza(jxa,jsep,ismain0)*te(jxa,jsep)+ti(jxa,jsep))+
     3 dpa0(topix(jxa,jsep),topiy(jxa,jsep),ismain0)*(
     4 rza(topix(jxa,jsep),topiy(jxa,jsep),ismain0)*
     5 te(topix(jxa,jsep),topiy(jxa,jsep))+
     6 ti(topix(jxa,jsep),topiy(jxa,jsep))))
      kesepm(1) = 0.5_R8 * (hce0(jxa,jsep)/ne(jxa,jsep)+
     1 hce0(topix(jxa,jsep),topiy(jxa,jsep))/
     2 ne(topix(jxa,jsep),topiy(jxa,jsep)))
      kisepm(1) = 0.5_R8 * (hci0(jxa,jsep)/ni(jxa,jsep,0)+
     1 hci0(topix(jxa,jsep),topiy(jxa,jsep))/
     2 ni(topix(jxa,jsep),topiy(jxa,jsep),0))
      vxsepm(1) = 0.5_R8 * (vla0(jxa,jsep,0,ismain)+
     1 vla0(topix(jxa,jsep),topiy(jxa,jsep),0,ismain))
      vysepm(1) = 0.5_R8 * (vla0(jxa,jsep,1,ismain)+
     1 vla0(topix(jxa,jsep),topiy(jxa,jsep),1,ismain))
      vssepm(1) = 0.5_R8 * (
     1 vsa0(jxa,jsep,ismain)/(mp*am(ismain)*na(jxa,jsep,ismain))+
     2 vsa0(topix(jxa,jsep),topiy(jxa,jsep),ismain)/
     3 (mp*am(ismain)*na(topix(jxa,jsep),topiy(jxa,jsep),ismain)))
      if(nnreg(0).ne.2) then
        nesepa(1) = 0.5_R8 *
     1   (ne(nx-target_offset,jsep)+
     2    ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
        tesepa(1) = 0.5_R8/ev *
     1   (te(nx-target_offset,jsep)+
     2    te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
        tisepa(1) = 0.5_R8/ev *
     1   (ti(nx-target_offset,jsep)+
     2    ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
        tpsepa(1) = 0.5_R8 *
     1   (target_temp(xymap(nx,jsep),1)+
     2    target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1))
        posepa(1) = 0.5_R8 *
     1   (po(nx-target_offset,jsep)+
     2    po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep)))
      else
        nesepa(1) = 
     1    ne(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
        tesepa(1) = 1.0_R8/ev *
     1    te(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
        tisepa(1) = 1.0_R8/ev *
     1    ti(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
        tpsepa(1) = 
     1   target_temp(xymap(topix(nx,jsep),topiy(nx,jsep)),1)
        posepa(1) = 
     1    po(topix(nx-target_offset,jsep),topiy(nx-target_offset,jsep))
      endif
      if(nncut.eq.2) then
        nesepi(2) = 0.5_R8 * (ne(ixtr+target_offset,jsep)+
     1   ne(topix(ixtr+target_offset,jsep),
     2      topiy(ixtr+target_offset,jsep)))
        tesepi(2) = 0.5_R8 * (te(ixtr+target_offset,jsep)+
     1   te(topix(ixtr+target_offset,jsep),
     2      topiy(ixtr+target_offset,jsep)))/ev
        tisepi(2) = 0.5_R8 * (ti(ixtr+target_offset,jsep)+
     1   ti(topix(ixtr+target_offset,jsep),
     2      topiy(ixtr+target_offset,jsep)))/ev
        tpsepi(2) = 0.5_R8 *
     1   (target_temp(xymap(ixtr,jsep),1)+
     2    target_temp(xymap(topix(ixtr,jsep),topiy(ixtr,jsep)),1))
        posepi(2) = 0.5_R8 * (po(ixtr+target_offset,jsep)+
     1   po(topix(ixtr+target_offset,jsep),
     2      topiy(ixtr+target_offset,jsep)))
        nesepm(2) = 0.5_R8 * (ne(jxi,jsep)+
     1   ne(topix(jxi,jsep),topiy(jxi,jsep)))
        tesepm(2) = 0.5_R8 * (te(jxi,jsep)+
     1   te(topix(jxi,jsep),topiy(jxi,jsep)))/ev
        tisepm(2) = 0.5_R8 * (ti(jxi,jsep)+
     1   ti(topix(jxi,jsep),topiy(jxi,jsep)))/ev
        posepm(2) = 0.5_R8 * (po(jxi,jsep)+
     1   po(topix(jxi,jsep),topiy(jxi,jsep)))
        dnsepm(2) = 0.5_R8 * (dna0(jxi,jsep,ismain)+
     1   dna0(topix(jxi,jsep),topiy(jxi,jsep),ismain))
        dpsepm(2) = 0.5_R8 * (
     1   dpa0(jxi,jsep,ismain0)*(
     2   rza(jxi,jsep,ismain0)*te(jxi,jsep)+ti(jxi,jsep))+
     3   dpa0(topix(jxi,jsep),topiy(jxi,jsep),ismain0)*(
     4   rza(topix(jxi,jsep),topiy(jxi,jsep),ismain0)*
     5   te(topix(jxi,jsep),topiy(jxi,jsep))+
     6   ti(topix(jxi,jsep),topiy(jxi,jsep))))
        kesepm(2) = 0.5_R8 * (hce0(jxi,jsep)/ne(jxi,jsep)+
     1   hce0(topix(jxi,jsep),topiy(jxi,jsep))/
     2   ne(topix(jxi,jsep),topiy(jxi,jsep)))
        kisepm(2) = 0.5_R8 * (hci0(jxi,jsep)/ni(jxi,jsep,0)+
     1   hci0(topix(jxi,jsep),topiy(jxi,jsep))/
     2   ni(topix(jxi,jsep),topiy(jxi,jsep),0))
        vxsepm(2) = 0.5_R8 * (vla0(jxi,jsep,0,ismain)+
     1   vla0(topix(jxi,jsep),topiy(jxi,jsep),0,ismain))
        vysepm(2) = 0.5_R8 * (vla0(jxi,jsep,1,ismain)+
     1   vla0(topix(jxi,jsep),topiy(jxi,jsep),1,ismain))
        vssepm(2) = 0.5_R8 * (
     1   vsa0(jxi,jsep,ismain)/(mp*am(ismain)*na(jxi,jsep,ismain))+
     2   vsa0(topix(jxi,jsep),topiy(jxi,jsep),ismain)/
     3   (mp*am(ismain)*na(topix(jxi,jsep),topiy(jxi,jsep),ismain)))
        nesepa(2) = 0.5_R8 * (ne(ixtl-target_offset,jsep)+
     1   ne(topix(ixtl-target_offset,jsep),
     2      topiy(ixtl-target_offset,jsep)))
        tesepa(2) = 0.5_R8 * (te(ixtl-target_offset,jsep)+
     1   te(topix(ixtl-target_offset,jsep),
     2      topiy(ixtl-target_offset,jsep)))/ev
        tisepa(2) = 0.5_R8 * (ti(ixtl-target_offset,jsep)+
     1   ti(topix(ixtl-target_offset,jsep),
     2      topiy(ixtl-target_offset,jsep)))/ev
        tpsepa(2) = 0.5_R8 *
     1   (target_temp(xymap(ixtl,jsep),1)+
     2    target_temp(xymap(topix(ixtl,jsep),topiy(ixtl,jsep)),1))
        posepa(2) = 0.5_R8 * (po(ixtl-target_offset,jsep)+
     1   po(topix(ixtl-target_offset,jsep),
     2      topiy(ixtl-target_offset,jsep)))
      endif
#endif
c
      do ic = 1, nc
        nemxip(ic) = nemxip(ic)
        temxip(ic) = temxip(ic)/ev
        timxip(ic) = timxip(ic)/ev
        pomxip(ic) = pomxip(ic)
        nemxap(ic) = nemxap(ic)
        temxap(ic) = temxap(ic)/ev
        timxap(ic) = timxap(ic)/ev
        pomxap(ic) = pomxap(ic)
        pwmxip(ic) = pwmxip(ic)
        pwmxap(ic) = pwmxap(ic)
      enddo
      
      tmne(1)=0.0_R8
      tmte(1)=0.0_R8
      tmti(1)=0.0_R8
      tmvol=0.0_R8
      do iy=-1,ny
        do ix=-1,nx
          if(region(ix,iy,0).ne.0) then
            tmne(1)=tmne(1)+ne(ix,iy)*vol(ix,iy)
            tmte(1)=tmte(1)+te(ix,iy)*ne(ix,iy)*vol(ix,iy)
            tmti(1)=tmti(1)+ti(ix,iy)*ni(ix,iy,0)*vol(ix,iy)
            tmvol=tmvol+vol(ix,iy)
          endif
        enddo
      enddo
      tmne(1)=tmne(1)
      tmte(1)=tmte(1)/ev
      tmti(1)=tmti(1)/ev

#ifndef NO_CDF      
      tmhacore(1)=0.0_R8
      tmhasol(1)=0.0_R8
      tmhadiv(1)=0.0_R8
#ifdef B25_EIRENE
c note no emission in guard cells and offset of array by 1
      do iy=-1,ny
        do ix=-1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &       bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 ) then
            if(mod(region(ix,iy,0),4).eq.1) then
              tmhacore(1)=tmhacore(1)+
     1         (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
            elseif(mod(region(ix,iy,0),4).eq.3 .or.
     1       (mod(region(ix,iy,0),4).eq.0 .and.
     2       region(ix,iy,0).ne.0)) then
              tmhadiv(1)=tmhadiv(1)+
     1         (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
            elseif(mod(region(ix,iy,0),4).eq.2) then
              tmhasol(1)=tmhasol(1)+
     1         (emiss(ix+1,iy+1,1,1)+emissmol(ix+1,iy+1,1,1))*vol(ix,iy)
            elseif(region(ix,iy,0).ne.0) then
              write(*,*) 'b2mwti: unknown region @ ',
     1         ix,iy,region(ix,iy,0)
            endif
          endif
        enddo
      enddo
#endif
      imap(1)=1
      tstepn(1) = ntstep
      call rwcdf(rw,ncid,'ntstep',imap,tstepn,iret)
      call rwcdf(rw,ncid,'timesa',imap,timesa,iret)

      imap(1)=1
      imap(2)=nncutmax
      call rwcdf(rw,ncid,'fnixip',imap,fnixip,iret)
      call rwcdf(rw,ncid,'feexip',imap,feexip,iret)
      call rwcdf(rw,ncid,'feixip',imap,feixip,iret)
      call rwcdf(rw,ncid,'fetxip',imap,fetxip,iret)
      call rwcdf(rw,ncid,'fchxip',imap,fchxip,iret)
      call rwcdf(rw,ncid,'fnixap',imap,fnixap,iret)
      call rwcdf(rw,ncid,'feexap',imap,feexap,iret)
      call rwcdf(rw,ncid,'feixap',imap,feixap,iret)
      call rwcdf(rw,ncid,'fetxap',imap,fetxap,iret)
      call rwcdf(rw,ncid,'fchxap',imap,fchxap,iret)

      call rwcdf(rw,ncid,'nesepi',imap,nesepi,iret)
      call rwcdf(rw,ncid,'tesepi',imap,tesepi,iret)
      call rwcdf(rw,ncid,'tisepi',imap,tisepi,iret)
      call rwcdf(rw,ncid,'posepi',imap,posepi,iret)
      call rwcdf(rw,ncid,'nesepm',imap,nesepm,iret)
      call rwcdf(rw,ncid,'tesepm',imap,tesepm,iret)
      call rwcdf(rw,ncid,'tisepm',imap,tisepm,iret)
      call rwcdf(rw,ncid,'posepm',imap,posepm,iret)
      call rwcdf(rw,ncid,'dnsepm',imap,dnsepm,iret)
      call rwcdf(rw,ncid,'dpsepm',imap,dpsepm,iret)
      call rwcdf(rw,ncid,'kesepm',imap,kesepm,iret)
      call rwcdf(rw,ncid,'kisepm',imap,kisepm,iret)
      call rwcdf(rw,ncid,'vxsepm',imap,vxsepm,iret)
      call rwcdf(rw,ncid,'vysepm',imap,vysepm,iret)
      call rwcdf(rw,ncid,'vssepm',imap,vssepm,iret)
      call rwcdf(rw,ncid,'nesepa',imap,nesepa,iret)
      call rwcdf(rw,ncid,'tesepa',imap,tesepa,iret)
      call rwcdf(rw,ncid,'tisepa',imap,tisepa,iret)
      call rwcdf(rw,ncid,'posepa',imap,posepa,iret)
      call rwcdf(rw,ncid,'nemxip',imap,nemxip,iret)
      call rwcdf(rw,ncid,'temxip',imap,temxip,iret)
      call rwcdf(rw,ncid,'timxip',imap,timxip,iret)
      call rwcdf(rw,ncid,'pomxip',imap,pomxip,iret)
      call rwcdf(rw,ncid,'nemxap',imap,nemxap,iret)
      call rwcdf(rw,ncid,'temxap',imap,temxap,iret)
      call rwcdf(rw,ncid,'timxap',imap,timxap,iret)
      call rwcdf(rw,ncid,'pomxap',imap,pomxap,iret)
      call rwcdf(rw,ncid,'fniyip',imap,fniyip,iret)
      call rwcdf(rw,ncid,'feeyip',imap,feeyip,iret)
      call rwcdf(rw,ncid,'feiyip',imap,feiyip,iret)
      call rwcdf(rw,ncid,'fetyip',imap,fetyip,iret)
      call rwcdf(rw,ncid,'fchyip',imap,fchyip,iret)
      call rwcdf(rw,ncid,'fniyap',imap,fniyap,iret)
      call rwcdf(rw,ncid,'feeyap',imap,feeyap,iret)
      call rwcdf(rw,ncid,'feiyap',imap,feiyap,iret)
      call rwcdf(rw,ncid,'fetyap',imap,fetyap,iret)
      call rwcdf(rw,ncid,'fchyap',imap,fchyap,iret)
      call rwcdf(rw,ncid,'pwmxip',imap,pwmxip,iret)
      call rwcdf(rw,ncid,'pwmxap',imap,pwmxap,iret)

      imap(1)=1
      call rwcdf(rw,ncid,'tmne',imap,tmne,iret)
      call rwcdf(rw,ncid,'tmte',imap,tmte,iret)
      call rwcdf(rw,ncid,'tmti',imap,tmti,iret)
      call rwcdf(rw,ncid,'tmhacore',imap,tmhacore,iret)
      call rwcdf(rw,ncid,'tmhasol',imap,tmhasol,iret)
      call rwcdf(rw,ncid,'tmhadiv',imap,tmhadiv,iret)

      imap(1)=1
      imap(2)=nncutmax
      call rwcdf(rw,ncid,'fnisip',imap,fnisip,iret)
      call rwcdf(rw,ncid,'feesip',imap,feesip,iret)
      call rwcdf(rw,ncid,'feisip',imap,feisip,iret)
      call rwcdf(rw,ncid,'fetsip',imap,fetsip,iret)
      call rwcdf(rw,ncid,'fchsip',imap,fchsip,iret)
      call rwcdf(rw,ncid,'fnisap',imap,fnisap,iret)
      call rwcdf(rw,ncid,'feesap',imap,feesap,iret)
      call rwcdf(rw,ncid,'feisap',imap,feisap,iret)
      call rwcdf(rw,ncid,'fetsap',imap,fetsap,iret)
      call rwcdf(rw,ncid,'fchsap',imap,fchsap,iret)
      call rwcdf(rw,ncid,'fnisipp',imap,fnisipp,iret)
      call rwcdf(rw,ncid,'feesipp',imap,feesipp,iret)
      call rwcdf(rw,ncid,'feisipp',imap,feisipp,iret)
      call rwcdf(rw,ncid,'fetsipp',imap,fetsipp,iret)
      call rwcdf(rw,ncid,'fchsipp',imap,fchsipp,iret)
      call rwcdf(rw,ncid,'fnisapp',imap,fnisapp,iret)
      call rwcdf(rw,ncid,'feesapp',imap,feesapp,iret)
      call rwcdf(rw,ncid,'feisapp',imap,feisapp,iret)
      call rwcdf(rw,ncid,'fetsapp',imap,fetsapp,iret)
      call rwcdf(rw,ncid,'fchsapp',imap,fchsapp,iret)
c
      imap(1)=nx+2     ! bl
      imap(2)=nx+2
      call rwcdf
     . (rw,ncid,'ne3dl',imap,ne(-1+target_offset,iylstrt),iret)
      call rwcdf
     . (rw,ncid,'te3dl',imap,te(-1+target_offset,iylstrt),iret)
      call rwcdf
     . (rw,ncid,'ti3dl',imap,ti(-1+target_offset,iylstrt),iret)
      call rwcdf
     . (rw,ncid,'po3dl',imap,po(-1+target_offset,iylstrt),iret)
      call rwcdf(rw,ncid,'fn3dl',imap,fna(0,iylstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dl',imap,fhe(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dl',imap,fhi(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dl',imap,fch(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dl',imap,fne(0,iylstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dl',imap,fni(0,iylstrt,0),iret)
      imap(1)=nx+2     ! i
      imap(2)=nx+2
      call rwcdf(rw,ncid,'ne3di',imap,ne(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'te3di',imap,te(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'ti3di',imap,ti(jxi,iyistrt),iret)
      call rwcdf(rw,ncid,'po3di',imap,po(jxi,iyistrt),iret)
      imap(1)=nx+2     ! a
      imap(2)=nx+2
      call rwcdf(rw,ncid,'ne3da',imap,ne(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'te3da',imap,te(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'ti3da',imap,ti(jxa,iyastrt),iret)
      call rwcdf(rw,ncid,'po3da',imap,po(jxa,iyastrt),iret)
      imap(1)=nx+2     ! br
      imap(2)=nx+2
      call rwcdf
     . (rw,ncid,'ne3dr',imap,ne(nx-target_offset,iyrstrt),iret)
      call rwcdf
     . (rw,ncid,'te3dr',imap,te(nx-target_offset,iyrstrt),iret)
      call rwcdf
     . (rw,ncid,'ti3dr',imap,ti(nx-target_offset,iyrstrt),iret)
      call rwcdf
     . (rw,ncid,'po3dr',imap,po(nx-target_offset,iyrstrt),iret)
      call rwcdf(rw,ncid,'fn3dr',imap,fna(nx,iyrstrt,0,ismain),iret)
      call rwcdf(rw,ncid,'fe3dr',imap,fhe(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fi3dr',imap,fhi(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fc3dr',imap,fch(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fl3dr',imap,fne(nx,iyrstrt,0),iret)
      call rwcdf(rw,ncid,'fo3dr',imap,fni(nx,iyrstrt,0),iret)
      if (nnreg(0).ge.8) then
        imap(1)=nx+2     ! tr
        imap(2)=nx+2
        call rwcdf
     .   (rw,ncid,'ne3dtr',imap,ne(ixtr+target_offset,iytrstrt),iret)
        call rwcdf
     .   (rw,ncid,'te3dtr',imap,te(ixtr+target_offset,iytrstrt),iret)
        call rwcdf
     .   (rw,ncid,'ti3dtr',imap,ti(ixtr+target_offset,iytrstrt),iret)
        call rwcdf
     .   (rw,ncid,'po3dtr',imap,po(ixtr+target_offset,iytrstrt),iret)
        call rwcdf
     .   (rw,ncid,'fn3dtr',imap,fna(ixtr+1,iytrstrt,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtr',imap,fhe(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fi3dtr',imap,fhi(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fc3dtr',imap,fch(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fl3dtr',imap,fne(ixtr+1,iytrstrt,0),iret)
        call rwcdf(rw,ncid,'fo3dtr',imap,fni(ixtr+1,iytrstrt,0),iret)
        imap(1)=nx+2     ! tl
        imap(2)=nx+2
        call rwcdf
     .   (rw,ncid,'ne3dtl',imap,ne(ixtl-target_offset,iytlstrt),iret)
        call rwcdf
     .   (rw,ncid,'te3dtl',imap,te(ixtl-target_offset,iytlstrt),iret)
        call rwcdf
     .   (rw,ncid,'ti3dtl',imap,ti(ixtl-target_offset,iytlstrt),iret)
        call rwcdf
     .   (rw,ncid,'po3dtl',imap,po(ixtl-target_offset,iytlstrt),iret)
        call rwcdf
     .   (rw,ncid,'fn3dtl',imap,fna(ixtl,iytlstrt,0,ismain),iret)
        call rwcdf(rw,ncid,'fe3dtl',imap,fhe(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fi3dtl',imap,fhi(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fc3dtl',imap,fch(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fl3dtl',imap,fne(ixtl,iytlstrt,0),iret)
        call rwcdf(rw,ncid,'fo3dtl',imap,fni(ixtl,iytlstrt,0),iret)
      endif
c
      imap(1)=1
      imap(2)=nx+2
      slice=0.0_R8
      slice(iylstrt:iylend)=na(-1+target_offset,iylstrt:iylend,0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(1,1:ny,1,1)
      call rwcdf(rw,ncid,'an3dl',imap,slice(iylstrt),iret)
      slice=0.0_R8
      slice(iyistrt:iyiend)=na(jxi,iyistrt:iyiend,0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(jxi+1,1:ny,1,1)
      call rwcdf(rw,ncid,'an3di',imap,slice(iyistrt),iret)
      slice=0.0_R8
      slice(iyastrt:iyaend)=na(jxa,iyastrt:iyaend,0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(jxa+1,1:ny,1,1)
      call rwcdf(rw,ncid,'an3da',imap,slice(iyastrt),iret)
      slice=0.0_R8
      slice(iyrstrt:iyrend)=na(nx-target_offset,iyrstrt:iyrend,0)
      slice(0:ny-1)=slice(0:ny-1)+dab2(nx,1:ny,1,1)
      call rwcdf(rw,ncid,'an3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.8) then
        slice(iytlstrt:iytlend)=
     .   na(ixtl-target_offset,iytlstrt:iytlend,0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(ixtl,1:ny,1,1)
        call rwcdf(rw,ncid,'an3dtl',imap,slice(iytlstrt),iret)
        slice(iytrstrt:iytrend)=
     .   na(ixtr+target_offset,iytrstrt:iytrend,0)
        slice(0:ny-1)=slice(0:ny-1)+dab2(ixtr+1,1:ny,1,1)
        call rwcdf(rw,ncid,'an3dtr',imap,slice(iytrstrt),iret)
      endif
      slice=0
      slice(0:ny-1)=dmb2(1,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3dl',imap,slice,iret)
      slice(0:ny-1)=dmb2(jxi+1,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3di',imap,slice,iret)
      slice(0:ny-1)=dmb2(jxa+1,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3da',imap,slice,iret)
      slice(0:ny-1)=dmb2(nx,1:ny,1,1)
      call rwcdf(rw,ncid,'mn3dr',imap,slice,iret)
      if (nnreg(0).ge.8) then
        slice(0:ny-1)=dmb2(ixtl,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dtl',imap,slice,iret)
        slice(0:ny-1)=dmb2(ixtr+1,1:ny,1,1)
        call rwcdf(rw,ncid,'mn3dtr',imap,slice,iret)
      endif
c      
      do iy = iylstrt, iylend
        slice(iy)=
     .   fhe(rightix(-1,iy),rightiy(-1,iy),0)+
     .   fhi(rightix(-1,iy),rightiy(-1,iy),0)
        do is=0,ns-1
          kintmp = (kinrgy(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     .              kinrgy(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     .             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
          rpttmp = (rpt(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     .              rpt(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     .             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
          slice(iy)=slice(iy)+(rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .     fna(rightix(-1,iy),rightiy(-1,iy),0,is)
        enddo
      enddo
      call rwcdf(rw,ncid,'ft3dl',imap,slice(iylstrt),iret)
      slice(iyrstrt:iyrend)=fhe(nx,iyrstrt:iyrend,0)+
     .                      fhi(nx,iyrstrt:iyrend,0)
      do iy=iyrstrt,iyrend
        do is=0,ns-1
          kintmp = (kinrgy(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     .              kinrgy(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     .             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
          rpttmp = (rpt(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     .              rpt(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     .             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
          slice(iy)=slice(iy)+(rpttmp*ev+kintmp*(1.0_R8-BoRiS))
     .     *fna(nx,iy,0,is)
        enddo
      enddo
      call rwcdf(rw,ncid,'ft3dr',imap,slice(iyrstrt),iret)
      if (nnreg(0).ge.8) then
        slice(iytlstrt:iytlend)=fhe(ixtl,iytlstrt:iytlend,0)+
     .                          fhi(ixtl,iytlstrt:iytlend,0)
        do iy=iytlstrt,iytlend
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ixtl,iy,is)*hx(leftix(ixtl,iy),leftiy(ixtl,iy))+
     .        kinrgy(leftix(ixtl,iy),leftiy(ixtl,iy),is)*hx(ixtl,iy))/
     .       (hx(leftix(ixtl,iy),leftiy(ixtl,iy))+hx(ixtl,iy))
            rpttmp = 
     .       (rpt(ixtl,iy,is)*hx(leftix(ixtl,iy),leftiy(ixtl,iy))+
     .        rpt(leftix(ixtl,iy),leftiy(ixtl,iy),is)*hx(ixtl,iy))/
     .       (hx(leftix(ixtl,iy),leftiy(ixtl,iy))+hx(ixtl,iy))
            slice(iy)=slice(iy)+(rpttmp*ev+kintmp*(1.0_R8-BoRiS))
     .       *fna(ixtl,iy,0,is)
          enddo
        enddo
        call rwcdf(rw,ncid,'ft3dtl',imap,slice(iytlstrt),iret)
        do iy = iytrstrt, iytrend
          slice(iy)=
     .     fhe(rightix(ixtr,iy),rightiy(ixtr,iy),0)+
     .     fhi(rightix(ixtr,iy),rightiy(ixtr,iy),0)
          do is=0,ns-1
            kintmp = 
     .       (kinrgy(ixtr,iy,is)*hx(rightix(ixtr,iy),rightiy(ixtr,iy))+
     .        kinrgy(rightix(ixtr,iy),rightiy(ixtr,iy),is)*hx(ixtr,iy))/
     .       (hx(rightix(ixtr,iy),rightiy(ixtr,iy))+hx(ixtr,iy))
            rpttmp = 
     .       (rpt(ixtr,iy,is)*hx(rightix(ixtr,iy),rightiy(ixtr,iy))+
     .        rpt(rightix(ixtr,iy),rightiy(ixtr,iy),is)*hx(ixtr,iy))/
     .       (hx(rightix(ixtr,iy),rightiy(ixtr,iy))+hx(ixtr,iy))
            slice(iy)=slice(iy)+(rpttmp*ev+kintmp*(1.0_R8-BoRiS))*
     .       fna(rightix(ixtr,iy),rightiy(ixtr,iy),0,is)
          enddo
        enddo
        call rwcdf(rw,ncid,'ft3dtr',imap,slice(iytrstrt),iret)
      endif
      slice(-1:ny)=dna0(jxi,-1:ny,ismain)
      call rwcdf(rw,ncid,'dn3di',imap,slice,iret)
      slice(-1:ny)=dna0(jxa,-1:ny,ismain)
      call rwcdf(rw,ncid,'dn3da',imap,slice,iret)
      slice(-1:ny)=dpa0(jxi,-1:ny,ismain0)*
     1           (rza(jxi,-1:ny,ismain0)*te(jxi,-1:ny)+ti(jxi,-1:ny))
      call rwcdf(rw,ncid,'dp3di',imap,slice,iret)
      slice(-1:ny)=dpa0(jxa,-1:ny,ismain0)*
     1           (rza(jxa,-1:ny,ismain0)*te(jxa,-1:ny)+ti(jxa,-1:ny))
      call rwcdf(rw,ncid,'dp3da',imap,slice,iret)
      slice(-1:ny)=fllim0fhi(jxi,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'lh3di',imap,slice,iret)
      slice(-1:ny)=fllim0fhi(jxa,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'lh3da',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxi,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'ln3di',imap,slice,iret)
      slice(-1:ny)=fllim0fna(jxa,-1:ny,1,ismain0)
      call rwcdf(rw,ncid,'ln3da',imap,slice,iret)
      slice(-1:ny)=hce0(jxi,-1:ny)/ne(jxi,-1:ny)
      call rwcdf(rw,ncid,'ke3di',imap,slice,iret)
      slice(-1:ny)=hce0(jxa,-1:ny)/ne(jxa,-1:ny)
      call rwcdf(rw,ncid,'ke3da',imap,slice,iret)
      slice(-1:ny)=hci0(jxi,-1:ny)/ni(jxi,-1:ny,0)
      call rwcdf(rw,ncid,'ki3di',imap,slice,iret)
      slice(-1:ny)=hci0(jxa,-1:ny)/ni(jxa,-1:ny,0)
      call rwcdf(rw,ncid,'ki3da',imap,slice,iret)
      slice(-1:ny)=vla0(jxi,-1:ny,0,ismain)
      call rwcdf(rw,ncid,'vx3di',imap,slice,iret)
      slice(-1:ny)=vla0(jxa,-1:ny,0,ismain)
      call rwcdf(rw,ncid,'vx3da',imap,slice,iret)
      slice(-1:ny)=vla0(jxi,-1:ny,1,ismain)
      call rwcdf(rw,ncid,'vy3di',imap,slice,iret)
      slice(-1:ny)=vla0(jxa,-1:ny,1,ismain)
      call rwcdf(rw,ncid,'vy3da',imap,slice,iret)
      slice(-1:ny)=vsa0(jxi,-1:ny,ismain)/
     1 (mp*am(ismain)*na(jxi,-1:ny,ismain))
      call rwcdf(rw,ncid,'vs3di',imap,slice,iret)
      slice(-1:ny)=vsa0(jxa,-1:ny,ismain)/
     1 (mp*am(ismain)*na(jxa,-1:ny,ismain))
      call rwcdf(rw,ncid,'vs3da',imap,slice,iret)
c
      imap(1)=1
      imap(2)=nncutmax
      call rwcdf(rw,ncid,'tpsepi',imap,tpsepi,iret)
      call rwcdf(rw,ncid,'tpsepa',imap,tpsepa,iret)
      call rwcdf(rw,ncid,'tpmxip',imap,tpmxip,iret)
      call rwcdf(rw,ncid,'tpmxap',imap,tpmxap,iret)
      imap(1)=1
      imap(2)=nx+2
      slice=0.0_R8
      slice(-1)=target_temp(xymap(-1,0),1)
      slice(0:ny-1)=target_temp(xymap(-1,0:ny-1),1)
      slice(ny)=target_temp(xymap(-1,ny-1),1)
      call rwcdf(rw,ncid,'tp3dl',imap,slice,iret)
      slice=0.0_R8
      slice(-1)=target_temp(xymap(nx,0),1)
      slice(0:ny-1)=target_temp(xymap(nx,0:ny-1),1)
      slice(ny)=target_temp(xymap(nx,ny-1),1)
      call rwcdf(rw,ncid,'tp3dr',imap,slice,iret)
      if (nnreg(0).ge.8) then
        slice=0.0_R8
        slice(-1)=target_temp(xymap(ixtl,0),1)
        slice(0:ny-1)=target_temp(xymap(ixtl,0:ny-1),1)
        slice(ny)=target_temp(xymap(ixtl,ny-1),1)
        call rwcdf(rw,ncid,'tp3dtl',imap,slice,iret)
        slice=0.0_R8
        slice(-1)=target_temp(xymap(ixtr,0),1)
        slice(0:ny-1)=target_temp(xymap(ixtr,0:ny-1),1)
        slice(ny)=target_temp(xymap(ixtr,ny-1),1)
        call rwcdf(rw,ncid,'tp3dtr',imap,slice,iret)
      endif
c      
      call ncclos (ncid, iret)
#endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mwti

      end subroutine b2mwti

#ifndef NO_CDF
#ifdef IBMaix
@process noextchk
#endif
      subroutine b2crtimecdf(nx, ny, nybl, nytl, nytr, nybr, nya, nyi, 
     1  nc, ns, iret)
      use b2mod_types
      use b2mod_constants
#     include <netcdf.inc>
      integer nx, ny, nybl, nytl, nytr, nybr, nya, nyi, nc, ns, iret
* netcdf id
      integer  ncid
* dimension ids
      integer  nxdim, nydim, nsdim, timedim,
     1 nybldim,nytldim,nytrdim,nybrdim,nyadim,nyidim,ncdim
* variable ids
      integer  ntstepid, timesaid, fnixipid, feexipid, feixipid, 
     1 fnixapid, feexapid, feixapid, nesepiid, tesepiid, tisepiid,
     2 nesepmid, tesepmid, tisepmid, nesepaid, tesepaid, tisepaid,
     3 nemxipid, temxipid, timxipid, nemxapid, temxapid, timxapid,
     4 fniyipid, feeyipid, feiyipid, fniyapid, feeyapid, feiyapid,
     5 pwmxipid, pwmxapid, tmneid, tmteid, tmtiid, tmhacoreid,
     6 tmhasolid, tmhadivid, fnisipid, feesipid, feisipid, fnisapid,
     7 feesapid, feisapid, fnisippid, feesippid, feisippid, fnisappid,
     8 feesappid, feisappid, ne3dlid, te3dlid, ti3dlid, ne3diid,
     9 te3diid, ti3diid, ne3daid, te3daid, ti3daid, ne3drid, te3drid,
     A ti3drid, an3dlid, mn3dlid, an3diid, mn3diid, an3daid, mn3daid,
     B an3drid, mn3drid, fn3dlid, fe3dlid, fi3dlid, fn3drid, fe3drid,
     C fi3drid, ne3dtlid, te3dtlid, ti3dtlid, ne3dtrid, te3dtrid,
     D ti3dtrid, an3dtlid, mn3dtlid, 
     E an3dtrid, mn3dtrid, fn3dtlid, fe3dtlid, fi3dtlid, fn3dtrid, 
     F fe3dtrid, fi3dtrid, fetxipid, fetxapid, fetyipid, fetyapid,
     G fetsipid, fetsapid, fetsippid, fetsappid
      integer  fchxipid, fchxapid, posepiid, posepmid, posepaid,
     3 pomxipid, pomxapid, fchyipid, fchyapid,
     5 fchsipid, fchsapid, fchsippid, fchsappid, po3dlid, 
     9 po3diid, po3daid, po3drid, fc3dlid, fc3drid,
     1 fl3dlid, fl3drid, fo3dlid, fo3drid,
     1 ft3dlid, ft3drid, po3dtlid, po3dtrid, fc3dtlid, fc3dtrid,
     1 fl3dtlid, fl3dtrid, fo3dtlid, fo3dtrid,
     1 ft3dtlid, ft3dtrid,
     2 dn3diid, dn3daid, dp3diid, dp3daid, ke3diid, ke3daid,
     3 ki3diid, ki3daid, vx3diid, vx3daid, vy3diid, vy3daid,
     4 vs3diid, vs3daid, lh3diid, lh3daid, ln3diid, ln3daid,
     5 dnsepmid, dpsepmid, kesepmid, kisepmid,
     6 vxsepmid, vysepmid, vssepmid,
     7 tpmxipid, tpmxapid, tp3drid, tp3dlid, tp3dtlid, tp3dtrid,
     8 tpsepiid, tpsepaid
* variable shapes
      integer dims(2)
* attribute vectors
      real (kind=R8) ::
     *   doubleval(1)
* enter define mode
      ncid = nccre ('b2time.nc', ncclob, iret)
* define dimensions
      nxdim = ncddef(ncid, 'nx', nx+2, iret)
      nydim = ncddef(ncid, 'ny', ny+2, iret)
      nybldim = ncddef (ncid, 'nybl', nybl, iret)
      if(nytl.gt.0) nytldim = ncddef (ncid, 'nytl', nytl, iret)
      if(nytr.gt.0) nytrdim = ncddef (ncid, 'nytr', nytr, iret)
      nybrdim = ncddef (ncid, 'nybr', nybr, iret)
      nyidim = ncddef (ncid, 'nyi', nyi, iret)
      nyadim = ncddef (ncid, 'nya', nya, iret)
      ncdim = ncddef(ncid, 'nc', nc, iret)
      nsdim = ncddef(ncid, 'ns', ns, iret)
      timedim = ncddef(ncid, 'time', ncunlim, iret)
* define variables
      dims(1) = 1
      ntstepid = ncvdef (ncid, 'ntstep', NCDOUBLE, 1, dims, iret)
      dims(1) = timedim
      timesaid = ncvdef (ncid, 'timesa', NCDOUBLE, 1, dims, iret)
      dims(1) = ncdim
      dims(2) = timedim
      fnixipid = ncvdef (ncid, 'fnixip', NCDOUBLE, 2, dims, iret)
      feexipid = ncvdef (ncid, 'feexip', NCDOUBLE, 2, dims, iret)
      feixipid = ncvdef (ncid, 'feixip', NCDOUBLE, 2, dims, iret)
      fetxipid = ncvdef (ncid, 'fetxip', NCDOUBLE, 2, dims, iret)
      fchxipid = ncvdef (ncid, 'fchxip', NCDOUBLE, 2, dims, iret)
      fnixapid = ncvdef (ncid, 'fnixap', NCDOUBLE, 2, dims, iret)
      feexapid = ncvdef (ncid, 'feexap', NCDOUBLE, 2, dims, iret)
      feixapid = ncvdef (ncid, 'feixap', NCDOUBLE, 2, dims, iret)
      fetxapid = ncvdef (ncid, 'fetxap', NCDOUBLE, 2, dims, iret)
      fchxapid = ncvdef (ncid, 'fchxap', NCDOUBLE, 2, dims, iret)
      nesepiid = ncvdef (ncid, 'nesepi', NCDOUBLE, 2, dims, iret)
      tesepiid = ncvdef (ncid, 'tesepi', NCDOUBLE, 2, dims, iret)
      tisepiid = ncvdef (ncid, 'tisepi', NCDOUBLE, 2, dims, iret)
      tpsepiid = ncvdef (ncid, 'tpsepi', NCDOUBLE, 2, dims, iret)
      posepiid = ncvdef (ncid, 'posepi', NCDOUBLE, 2, dims, iret)
      nesepmid = ncvdef (ncid, 'nesepm', NCDOUBLE, 2, dims, iret)
      tesepmid = ncvdef (ncid, 'tesepm', NCDOUBLE, 2, dims, iret)
      tisepmid = ncvdef (ncid, 'tisepm', NCDOUBLE, 2, dims, iret)
      posepmid = ncvdef (ncid, 'posepm', NCDOUBLE, 2, dims, iret)
      dnsepmid = ncvdef (ncid, 'dnsepm', NCDOUBLE, 2, dims, iret)
      dpsepmid = ncvdef (ncid, 'dpsepm', NCDOUBLE, 2, dims, iret)
      kesepmid = ncvdef (ncid, 'kesepm', NCDOUBLE, 2, dims, iret)
      kisepmid = ncvdef (ncid, 'kisepm', NCDOUBLE, 2, dims, iret)
      vxsepmid = ncvdef (ncid, 'vxsepm', NCDOUBLE, 2, dims, iret)
      vysepmid = ncvdef (ncid, 'vysepm', NCDOUBLE, 2, dims, iret)
      vssepmid = ncvdef (ncid, 'vssepm', NCDOUBLE, 2, dims, iret)
      nesepaid = ncvdef (ncid, 'nesepa', NCDOUBLE, 2, dims, iret)
      tesepaid = ncvdef (ncid, 'tesepa', NCDOUBLE, 2, dims, iret)
      tisepaid = ncvdef (ncid, 'tisepa', NCDOUBLE, 2, dims, iret)
      tpsepaid = ncvdef (ncid, 'tpsepa', NCDOUBLE, 2, dims, iret)
      posepaid = ncvdef (ncid, 'posepa', NCDOUBLE, 2, dims, iret)
      nemxipid = ncvdef (ncid, 'nemxip', NCDOUBLE, 2, dims, iret)
      temxipid = ncvdef (ncid, 'temxip', NCDOUBLE, 2, dims, iret)
      timxipid = ncvdef (ncid, 'timxip', NCDOUBLE, 2, dims, iret)
      tpmxipid = ncvdef (ncid, 'tpmxip', NCDOUBLE, 2, dims, iret)
      pomxipid = ncvdef (ncid, 'pomxip', NCDOUBLE, 2, dims, iret)
      nemxapid = ncvdef (ncid, 'nemxap', NCDOUBLE, 2, dims, iret)
      temxapid = ncvdef (ncid, 'temxap', NCDOUBLE, 2, dims, iret)
      timxapid = ncvdef (ncid, 'timxap', NCDOUBLE, 2, dims, iret)
      tpmxapid = ncvdef (ncid, 'tpmxap', NCDOUBLE, 2, dims, iret)
      pomxapid = ncvdef (ncid, 'pomxap', NCDOUBLE, 2, dims, iret)
      fniyipid = ncvdef (ncid, 'fniyip', NCDOUBLE, 2, dims, iret)
      feeyipid = ncvdef (ncid, 'feeyip', NCDOUBLE, 2, dims, iret)
      feiyipid = ncvdef (ncid, 'feiyip', NCDOUBLE, 2, dims, iret)
      fetyipid = ncvdef (ncid, 'fetyip', NCDOUBLE, 2, dims, iret)
      fchyipid = ncvdef (ncid, 'fchyip', NCDOUBLE, 2, dims, iret)
      fniyapid = ncvdef (ncid, 'fniyap', NCDOUBLE, 2, dims, iret)
      feeyapid = ncvdef (ncid, 'feeyap', NCDOUBLE, 2, dims, iret)
      feiyapid = ncvdef (ncid, 'feiyap', NCDOUBLE, 2, dims, iret)
      fetyapid = ncvdef (ncid, 'fetyap', NCDOUBLE, 2, dims, iret)
      fchyapid = ncvdef (ncid, 'fchyap', NCDOUBLE, 2, dims, iret)
      pwmxipid = ncvdef (ncid, 'pwmxip', NCDOUBLE, 2, dims, iret)
      pwmxapid = ncvdef (ncid, 'pwmxap', NCDOUBLE, 2, dims, iret)
      dims(1) = timedim
      tmneid = ncvdef (ncid, 'tmne', NCDOUBLE, 1, dims, iret)
      tmteid = ncvdef (ncid, 'tmte', NCDOUBLE, 1, dims, iret)
      tmtiid = ncvdef (ncid, 'tmti', NCDOUBLE, 1, dims, iret)
      tmhacoreid = ncvdef (ncid, 'tmhacore', NCDOUBLE, 1, dims, iret)
      tmhasolid = ncvdef (ncid, 'tmhasol', NCDOUBLE, 1, dims, iret)
      tmhadivid = ncvdef (ncid, 'tmhadiv', NCDOUBLE, 1, dims, iret)
      dims(1) = ncdim
      dims(2) = timedim
      fnisipid = ncvdef (ncid, 'fnisip', NCDOUBLE, 2, dims, iret)
      feesipid = ncvdef (ncid, 'feesip', NCDOUBLE, 2, dims, iret)
      feisipid = ncvdef (ncid, 'feisip', NCDOUBLE, 2, dims, iret)
      fetsipid = ncvdef (ncid, 'fetsip', NCDOUBLE, 2, dims, iret)
      fchsipid = ncvdef (ncid, 'fchsip', NCDOUBLE, 2, dims, iret)
      fnisapid = ncvdef (ncid, 'fnisap', NCDOUBLE, 2, dims, iret)
      feesapid = ncvdef (ncid, 'feesap', NCDOUBLE, 2, dims, iret)
      feisapid = ncvdef (ncid, 'feisap', NCDOUBLE, 2, dims, iret)
      fetsapid = ncvdef (ncid, 'fetsap', NCDOUBLE, 2, dims, iret)
      fchsapid = ncvdef (ncid, 'fchsap', NCDOUBLE, 2, dims, iret)
      fnisippid = ncvdef (ncid, 'fnisipp', NCDOUBLE, 2, dims, iret)
      feesippid = ncvdef (ncid, 'feesipp', NCDOUBLE, 2, dims, iret)
      feisippid = ncvdef (ncid, 'feisipp', NCDOUBLE, 2, dims, iret)
      fetsippid = ncvdef (ncid, 'fetsipp', NCDOUBLE, 2, dims, iret)
      fchsippid = ncvdef (ncid, 'fchsipp', NCDOUBLE, 2, dims, iret)
      fnisappid = ncvdef (ncid, 'fnisapp', NCDOUBLE, 2, dims, iret)
      feesappid = ncvdef (ncid, 'feesapp', NCDOUBLE, 2, dims, iret)
      feisappid = ncvdef (ncid, 'feisapp', NCDOUBLE, 2, dims, iret)
      fetsappid = ncvdef (ncid, 'fetsapp', NCDOUBLE, 2, dims, iret)
      fchsappid = ncvdef (ncid, 'fchsapp', NCDOUBLE, 2, dims, iret)
      dims(2) = timedim
      dims(1) = nybldim
      fn3dlid = ncvdef (ncid, 'fn3dl', NCDOUBLE, 2, dims, iret)
      fe3dlid = ncvdef (ncid, 'fe3dl', NCDOUBLE, 2, dims, iret)
      fi3dlid = ncvdef (ncid, 'fi3dl', NCDOUBLE, 2, dims, iret)
      ft3dlid = ncvdef (ncid, 'ft3dl', NCDOUBLE, 2, dims, iret)
      fc3dlid = ncvdef (ncid, 'fc3dl', NCDOUBLE, 2, dims, iret)
      fl3dlid = ncvdef (ncid, 'fl3dl', NCDOUBLE, 2, dims, iret)
      fo3dlid = ncvdef (ncid, 'fo3dl', NCDOUBLE, 2, dims, iret)
      ne3dlid = ncvdef (ncid, 'ne3dl', NCDOUBLE, 2, dims, iret)
      te3dlid = ncvdef (ncid, 'te3dl', NCDOUBLE, 2, dims, iret)
      ti3dlid = ncvdef (ncid, 'ti3dl', NCDOUBLE, 2, dims, iret)
      po3dlid = ncvdef (ncid, 'po3dl', NCDOUBLE, 2, dims, iret)
      an3dlid = ncvdef (ncid, 'an3dl', NCDOUBLE, 2, dims, iret)
      mn3dlid = ncvdef (ncid, 'mn3dl', NCDOUBLE, 2, dims, iret)
      tp3dlid = ncvdef (ncid, 'tp3dl', NCDOUBLE, 2, dims, iret)
      if(nytl.gt.0) then
        dims(2) = timedim
        dims(1) = nytldim
        fn3dtlid = ncvdef (ncid, 'fn3dtl', NCDOUBLE, 2, dims, iret)
        fe3dtlid = ncvdef (ncid, 'fe3dtl', NCDOUBLE, 2, dims, iret)
        fi3dtlid = ncvdef (ncid, 'fi3dtl', NCDOUBLE, 2, dims, iret)
        ft3dtlid = ncvdef (ncid, 'ft3dtl', NCDOUBLE, 2, dims, iret)
        fc3dtlid = ncvdef (ncid, 'fc3dtl', NCDOUBLE, 2, dims, iret)
        fl3dtlid = ncvdef (ncid, 'fl3dtl', NCDOUBLE, 2, dims, iret)
        fo3dtlid = ncvdef (ncid, 'fo3dtl', NCDOUBLE, 2, dims, iret)
        ne3dtlid = ncvdef (ncid, 'ne3dtl', NCDOUBLE, 2, dims, iret)
        te3dtlid = ncvdef (ncid, 'te3dtl', NCDOUBLE, 2, dims, iret)
        ti3dtlid = ncvdef (ncid, 'ti3dtl', NCDOUBLE, 2, dims, iret)
        po3dtlid = ncvdef (ncid, 'po3dtl', NCDOUBLE, 2, dims, iret)
        an3dtlid = ncvdef (ncid, 'an3dtl', NCDOUBLE, 2, dims, iret)
        mn3dtlid = ncvdef (ncid, 'mn3dtl', NCDOUBLE, 2, dims, iret)
        tp3dtlid = ncvdef (ncid, 'tp3dtl', NCDOUBLE, 2, dims, iret)
      endif
      dims(2) = timedim
      dims(1) = nyidim
      ne3diid = ncvdef (ncid, 'ne3di', NCDOUBLE, 2, dims, iret)
      te3diid = ncvdef (ncid, 'te3di', NCDOUBLE, 2, dims, iret)
      ti3diid = ncvdef (ncid, 'ti3di', NCDOUBLE, 2, dims, iret)
      po3diid = ncvdef (ncid, 'po3di', NCDOUBLE, 2, dims, iret)
      an3diid = ncvdef (ncid, 'an3di', NCDOUBLE, 2, dims, iret)
      mn3diid = ncvdef (ncid, 'mn3di', NCDOUBLE, 2, dims, iret)
      dn3diid = ncvdef (ncid, 'dn3di', NCDOUBLE, 2, dims, iret)
      dp3diid = ncvdef (ncid, 'dp3di', NCDOUBLE, 2, dims, iret)
      lh3diid = ncvdef (ncid, 'lh3di', NCDOUBLE, 2, dims, iret)
      ln3diid = ncvdef (ncid, 'ln3di', NCDOUBLE, 2, dims, iret)
      ke3diid = ncvdef (ncid, 'ke3di', NCDOUBLE, 2, dims, iret)
      ki3diid = ncvdef (ncid, 'ki3di', NCDOUBLE, 2, dims, iret)
      vx3diid = ncvdef (ncid, 'vx3di', NCDOUBLE, 2, dims, iret)
      vy3diid = ncvdef (ncid, 'vy3di', NCDOUBLE, 2, dims, iret)
      vs3diid = ncvdef (ncid, 'vs3di', NCDOUBLE, 2, dims, iret)
      dims(2) = timedim
      dims(1) = nyadim
      ne3daid = ncvdef (ncid, 'ne3da', NCDOUBLE, 2, dims, iret)
      te3daid = ncvdef (ncid, 'te3da', NCDOUBLE, 2, dims, iret)
      ti3daid = ncvdef (ncid, 'ti3da', NCDOUBLE, 2, dims, iret)
      po3daid = ncvdef (ncid, 'po3da', NCDOUBLE, 2, dims, iret)
      an3daid = ncvdef (ncid, 'an3da', NCDOUBLE, 2, dims, iret)
      mn3daid = ncvdef (ncid, 'mn3da', NCDOUBLE, 2, dims, iret)
      dn3daid = ncvdef (ncid, 'dn3da', NCDOUBLE, 2, dims, iret)
      dp3daid = ncvdef (ncid, 'dp3da', NCDOUBLE, 2, dims, iret)
      lh3daid = ncvdef (ncid, 'lh3da', NCDOUBLE, 2, dims, iret)
      ln3daid = ncvdef (ncid, 'ln3da', NCDOUBLE, 2, dims, iret)
      ke3daid = ncvdef (ncid, 'ke3da', NCDOUBLE, 2, dims, iret)
      ki3daid = ncvdef (ncid, 'ki3da', NCDOUBLE, 2, dims, iret)
      vx3daid = ncvdef (ncid, 'vx3da', NCDOUBLE, 2, dims, iret)
      vy3daid = ncvdef (ncid, 'vy3da', NCDOUBLE, 2, dims, iret)
      vs3daid = ncvdef (ncid, 'vs3da', NCDOUBLE, 2, dims, iret)
      dims(2) = timedim
      dims(1) = nybrdim
      ne3drid = ncvdef (ncid, 'ne3dr', NCDOUBLE, 2, dims, iret)
      te3drid = ncvdef (ncid, 'te3dr', NCDOUBLE, 2, dims, iret)
      ti3drid = ncvdef (ncid, 'ti3dr', NCDOUBLE, 2, dims, iret)
      po3drid = ncvdef (ncid, 'po3dr', NCDOUBLE, 2, dims, iret)
      an3drid = ncvdef (ncid, 'an3dr', NCDOUBLE, 2, dims, iret)
      mn3drid = ncvdef (ncid, 'mn3dr', NCDOUBLE, 2, dims, iret)
      fn3drid = ncvdef (ncid, 'fn3dr', NCDOUBLE, 2, dims, iret)
      fe3drid = ncvdef (ncid, 'fe3dr', NCDOUBLE, 2, dims, iret)
      fi3drid = ncvdef (ncid, 'fi3dr', NCDOUBLE, 2, dims, iret)
      ft3drid = ncvdef (ncid, 'ft3dr', NCDOUBLE, 2, dims, iret)
      fc3drid = ncvdef (ncid, 'fc3dr', NCDOUBLE, 2, dims, iret)
      fl3drid = ncvdef (ncid, 'fl3dr', NCDOUBLE, 2, dims, iret)
      fo3drid = ncvdef (ncid, 'fo3dr', NCDOUBLE, 2, dims, iret)
      tp3drid = ncvdef (ncid, 'tp3dr', NCDOUBLE, 2, dims, iret)
      if(nytr.gt.0) then
        dims(2) = timedim
        dims(1) = nytrdim
        ne3dtrid = ncvdef (ncid, 'ne3dtr', NCDOUBLE, 2, dims, iret)
        te3dtrid = ncvdef (ncid, 'te3dtr', NCDOUBLE, 2, dims, iret)
        ti3dtrid = ncvdef (ncid, 'ti3dtr', NCDOUBLE, 2, dims, iret)
        po3dtrid = ncvdef (ncid, 'po3dtr', NCDOUBLE, 2, dims, iret)
        an3dtrid = ncvdef (ncid, 'an3dtr', NCDOUBLE, 2, dims, iret)
        mn3dtrid = ncvdef (ncid, 'mn3dtr', NCDOUBLE, 2, dims, iret)
        fn3dtrid = ncvdef (ncid, 'fn3dtr', NCDOUBLE, 2, dims, iret)
        fe3dtrid = ncvdef (ncid, 'fe3dtr', NCDOUBLE, 2, dims, iret)
        fi3dtrid = ncvdef (ncid, 'fi3dtr', NCDOUBLE, 2, dims, iret)
        ft3dtrid = ncvdef (ncid, 'ft3dtr', NCDOUBLE, 2, dims, iret)
        fc3dtrid = ncvdef (ncid, 'fc3dtr', NCDOUBLE, 2, dims, iret)
        fl3dtrid = ncvdef (ncid, 'fl3dtr', NCDOUBLE, 2, dims, iret)
        fo3dtrid = ncvdef (ncid, 'fo3dtr', NCDOUBLE, 2, dims, iret)
        tp3dtrid = ncvdef (ncid, 'tp3dtr', NCDOUBLE, 2, dims, iret)
      endif
* assign attributes
      call ncaptc(ncid, timesaid, 'long_name', NCCHAR, 4, 'time', iret)
      call ncaptc(ncid, timesaid, 'units', NCCHAR, 2, 's ', iret)
      call ncaptc(ncid, fnixipid, 'long_name', NCCHAR, 51, 
     1'integrated poloidal particle flux, inboard divertor', iret)
      call ncaptc(ncid, fnixipid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feexipid, 'long_name', NCCHAR, 58,
     1'integrated poloidal electron energy flux, inboard divertor',
     2 iret)
      call ncaptc(ncid, feexipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feixipid, 'long_name', NCCHAR, 53,
     1'integrated poloidal ion energy flux, inboard divertor', iret)
      call ncaptc(ncid, feixipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetxipid, 'long_name', NCCHAR, 53,
     1'integrated poloidal tot energy flux, inboard divertor', iret)
      call ncaptc(ncid, fetxipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchxipid, 'long_name', NCCHAR, 45,
     1'integrated poloidal current, inboard divertor', iret)
      call ncaptc(ncid, fchxipid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, fnixapid, 'long_name', NCCHAR, 52,
     1'integrated poloidal particle flux, outboard divertor', iret)
      call ncaptc(ncid, fnixapid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feexapid, 'long_name', NCCHAR, 59,
     1'integrated poloidal electron energy flux, outboard divertor',
     2 iret)
      call ncaptc(ncid, feexapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feixapid, 'long_name', NCCHAR, 54,
     1'integrated poloidal ion energy flux, outboard divertor', iret)
      call ncaptc(ncid, feixapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetxapid, 'long_name', NCCHAR, 54,
     1'integrated poloidal tot energy flux, outboard divertor', iret)
      call ncaptc(ncid, fetxapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchxapid, 'long_name', NCCHAR, 46,
     1'integrated poloidal current, outboard divertor', iret)
      call ncaptc(ncid, fchxapid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, nesepiid, 'long_name', NCCHAR, 45,
     1'separatrix electron density, inboard divertor', iret)
      call ncaptc(ncid, nesepiid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, tesepiid, 'long_name', NCCHAR, 49,
     1'separatrix electron temperature, inboard divertor', iret)
      call ncaptc(ncid, tesepiid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tisepiid, 'long_name', NCCHAR, 44,
     1'separatrix ion temperature, inboard divertor', iret)
      call ncaptc(ncid, tisepiid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tpsepiid, 'long_name', NCCHAR, 46,
     1'separatrix plate temperature, inboard divertor', iret)
      call ncaptc(ncid, tpsepiid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, posepiid, 'long_name', NCCHAR, 38,
     1'separatrix potential, inboard divertor', iret)
      call ncaptc(ncid, posepiid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, nesepmid, 'long_name', NCCHAR, 43,
     1'separatrix electron density, outer midplane', iret)
      call ncaptc(ncid, nesepmid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, tesepmid, 'long_name', NCCHAR, 47,
     1'separatrix electron temperature, outer midplane', iret)
      call ncaptc(ncid, tesepmid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tisepmid, 'long_name', NCCHAR, 42,
     1'separatrix ion temperature, outer midplane', iret)
      call ncaptc(ncid, tisepmid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, posepmid, 'long_name', NCCHAR, 36,
     1'separatrix potential, outer midplane', iret)
      call ncaptc(ncid, posepmid, 'units', NCCHAR, 2, 'V ', iret)

      call ncaptc(ncid, dnsepmid, 'long_name', NCCHAR, 37,
     1'diffusion coefficient, outer midplane', iret)
      call ncaptc(ncid, dnsepmid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, dpsepmid, 'long_name', NCCHAR, 46,
     1'pressure diffusion coefficient, outer midplane', iret)
      call ncaptc(ncid, dpsepmid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, kesepmid, 'long_name', NCCHAR, 44,
     1'electron thermal diffusivity, outer midplane', iret)
      call ncaptc(ncid, kesepmid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, kisepmid, 'long_name', NCCHAR, 39,
     1'ion thermal diffusivity, outer midplane', iret)
      call ncaptc(ncid, kisepmid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, vxsepmid, 'long_name', NCCHAR, 39,
     1'poloidal pinch velocity, outer midplane', iret)
      call ncaptc(ncid, vxsepmid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vysepmid, 'long_name', NCCHAR, 37,
     1'radial pinch velocity, outer midplane', iret)
      call ncaptc(ncid, vysepmid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vssepmid, 'long_name', NCCHAR, 37,
     1'viscosity coefficient, outer midplane', iret)
      call ncaptc(ncid, vssepmid, 'units', NCCHAR, 12, 'm.kg^-1.s^-1',
     1 iret)


      call ncaptc(ncid, nesepaid, 'long_name', NCCHAR, 46,
     1'separatrix electron density, outboard divertor', iret)
      call ncaptc(ncid, nesepaid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, tesepaid, 'long_name', NCCHAR, 50,
     1'separatrix electron temperature, outboard divertor', iret)
      call ncaptc(ncid, tesepaid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tisepaid, 'long_name', NCCHAR, 45,
     1'separatrix ion temperature, outboard divertor', iret)
      call ncaptc(ncid, tisepaid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tpsepaid, 'long_name', NCCHAR, 47,
     1'separatrix plate temperature, outboard divertor', iret)
      call ncaptc(ncid, tpsepaid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, posepaid, 'long_name', NCCHAR, 43,
     1'separatrix potential, outboard divertor', iret)
      call ncaptc(ncid, posepaid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, nemxipid, 'long_name', NCCHAR, 42,
     1'maximum electron density, inboard divertor', iret)
      call ncaptc(ncid, nemxipid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, temxipid, 'long_name', NCCHAR, 46,
     1'maximum electron temperature, inboard divertor', iret)
      call ncaptc(ncid, temxipid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, timxipid, 'long_name', NCCHAR, 41,
     1'maximum ion temperature, inboard divertor', iret)
      call ncaptc(ncid, timxipid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tpmxipid, 'long_name', NCCHAR, 43,
     1'maximum plate temperature, inboard divertor', iret)
      call ncaptc(ncid, tpmxipid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, pomxipid, 'long_name', NCCHAR, 35,
     1'maximum potential, inboard divertor', iret)
      call ncaptc(ncid, pomxipid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, nemxapid, 'long_name', NCCHAR, 43,
     1'maximum electron density, outboard divertor', iret)
      call ncaptc(ncid, nemxapid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, temxapid, 'long_name', NCCHAR, 47,
     1'maximum electron temperature, outboard divertor', iret)
      call ncaptc(ncid, temxapid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, timxapid, 'long_name', NCCHAR, 43,
     1'maximum ion temperature, outboard divertor', iret)
      call ncaptc(ncid, timxapid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tpmxapid, 'long_name', NCCHAR, 45,
     1'maximum plate temperature, outboard divertor', iret)
      call ncaptc(ncid, tpmxapid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, pomxapid, 'long_name', NCCHAR, 37,
     1'maximum potential, outboard divertor', iret)
      call ncaptc(ncid, pomxapid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, fniyipid, 'long_name', NCCHAR, 45,
     1'integrated radial particle flux, main chamber', iret)
      call ncaptc(ncid, fniyipid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feeyipid, 'long_name', NCCHAR, 52,
     1'integrated radial electron energy flux, main chamber',
     2 iret)
      call ncaptc(ncid, feeyipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feiyipid, 'long_name', NCCHAR, 47,
     1'integrated radial ion energy flux, main chamber', iret)
      call ncaptc(ncid, feiyipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetyipid, 'long_name', NCCHAR, 47,
     1'integrated radial tot energy flux, main chamber', iret)
      call ncaptc(ncid, fetyipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchyipid, 'long_name', NCCHAR, 39,
     1'integrated radial current, main chamber', iret)
      call ncaptc(ncid, fchyipid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, fniyapid, 'long_name', NCCHAR, 48,
     1'integrated radial particle flux, divertor region', iret)
      call ncaptc(ncid, fniyapid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feeyapid, 'long_name', NCCHAR, 55,
     1'integrated radial electron energy flux, divertor region'
     2,iret)
      call ncaptc(ncid, feeyapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feiyapid, 'long_name', NCCHAR, 50,
     1'integrated radial ion energy flux, divertor region',
     2 iret)
      call ncaptc(ncid, feiyapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetyapid, 'long_name', NCCHAR, 50,
     1'integrated radial tot energy flux, divertor region',
     2 iret)
      call ncaptc(ncid, fetyapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchyapid, 'long_name', NCCHAR, 42,
     1'integrated radial current, divertor region',
     2 iret)
      call ncaptc(ncid, fchyapid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, pwmxipid, 'long_name', NCCHAR, 42,
     1'maximum total power flux, inboard divertor', iret)
      call ncaptc(ncid, pwmxipid, 'units', NCCHAR, 6, 'W.m^-2', iret)
      call ncaptc(ncid, pwmxapid, 'long_name', NCCHAR, 43,
     1'maximum total power flux, outboard divertor', iret)
      call ncaptc(ncid, pwmxapid, 'units', NCCHAR, 6, 'W.m^-2', iret)
      call ncaptc(ncid, tmneid, 'long_name', NCCHAR, 25,
     1'total number of particles', iret)
      call ncaptc(ncid, tmneid, 'units', NCCHAR, 2, '  ', iret)
      call ncaptc(ncid, tmteid, 'long_name', NCCHAR, 21,
     1'total electron energy', iret)
      call ncaptc(ncid, tmteid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tmtiid, 'long_name', NCCHAR, 16,
     1'total ion energy', iret)
      call ncaptc(ncid, tmtiid, 'units', NCCHAR, 2, 'eV', iret)
      call ncaptc(ncid, tmhacoreid, 'long_name', NCCHAR, 24,
     1'H-alpha emissivity, core', iret)
      call ncaptc(ncid, tmhacoreid, 'units', NCCHAR, 19,
     1'photons.m^-2.sr^-1?', iret)
      call ncaptc(ncid, tmhasolid, 'long_name', NCCHAR, 41,
     1'H-alpha emissivity, SOL above the x-point', iret)
      call ncaptc(ncid, tmhasolid, 'units', NCCHAR, 19,
     1'photons.m^-2.sr^-1?', iret)
      call ncaptc(ncid, tmhadivid, 'long_name', NCCHAR, 35,
     1'H-alpha emissivity, divertor region', iret)
      call ncaptc(ncid, tmhadivid, 'units', NCCHAR, 19,
     1'photons.m^-2.sr^-1?', iret)
      call ncaptc(ncid, fnisipid, 'long_name', NCCHAR, 54,
     1'poloidal particle flux, into inboard separatrix throat', iret)
      call ncaptc(ncid, fnisipid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feesipid, 'long_name', NCCHAR, 61,
     1'poloidal electron energy flux, into inboard separatrix throat'
     2,iret)
      call ncaptc(ncid, feesipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feisipid, 'long_name', NCCHAR, 56,
     1'poloidal ion energy flux, into inboard separatrix throat',
     2 iret)
      call ncaptc(ncid, feisipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetsipid, 'long_name', NCCHAR, 58,
     1'poloidal total energy flux, into inboard separatrix throat',
     2 iret)
      call ncaptc(ncid, fetsipid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchsipid, 'long_name', NCCHAR, 48,
     1'poloidal current, into inboard separatrix throat',
     2 iret)
      call ncaptc(ncid, fchsipid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, fnisapid, 'long_name', NCCHAR, 55,
     1'poloidal particle flux, into outboard separatrix throat',
     2 iret)
      call ncaptc(ncid, fnisapid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feesapid, 'long_name', NCCHAR, 62,
     1'poloidal electron energy flux,'//
     2' into outboard separatrix throat', iret)
      call ncaptc(ncid, feesapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feisapid, 'long_name', NCCHAR, 57,
     1'poloidal ion energy flux, into outboard separatrix throat', iret)
      call ncaptc(ncid, feisapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetsapid, 'long_name', NCCHAR, 59,
     1'poloidal total energy flux, into outboard separatrix throat',
     2 iret)
      call ncaptc(ncid, fetsapid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchsapid, 'long_name', NCCHAR, 49,
     1'poloidal current, into outboard separatrix throat', iret)
      call ncaptc(ncid, fchsapid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, fnisippid, 'long_name', NCCHAR, 40,
     1'poloidal particle flux, core x-pt region', iret)
      call ncaptc(ncid, fnisippid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feesippid, 'long_name', NCCHAR, 52,
     1'poloidal electron energy flux, core x-pt flux region',
     2 iret)
      call ncaptc(ncid, feesippid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feisippid, 'long_name', NCCHAR, 47,
     1'poloidal ion energy flux, core x-pt flux region', iret)
      call ncaptc(ncid, feisippid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetsippid, 'long_name', NCCHAR, 49,
     1'poloidal total energy flux, core x-pt flux region',
     2 iret)
      call ncaptc(ncid, fetsippid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchsippid, 'long_name', NCCHAR, 39,
     1'poloidal current, core x-pt flux region', iret)
      call ncaptc(ncid, fchsippid, 'units', NCCHAR, 2, 'A ', iret)
      call ncaptc(ncid, fnisappid, 'long_name', NCCHAR, 48,
     1'poloidal particle flux, x-pt private flux region', iret)
      call ncaptc(ncid, fnisappid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, feesappid, 'long_name', NCCHAR, 55,
     1'poloidal electron energy flux, x-pt private flux region',
     2iret)
      call ncaptc(ncid, feesappid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, feisappid, 'long_name', NCCHAR, 50,
     1'poloidal ion energy flux, x-pt private flux region', iret)
      call ncaptc(ncid, feisappid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fetsappid, 'long_name', NCCHAR, 52,
     1'poloidal total energy flux, x-pt private flux region',
     2 iret)
      call ncaptc(ncid, fetsappid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fchsappid, 'long_name', NCCHAR, 42,
     1'poloidal current, x-pt private flux region', iret)
      call ncaptc(ncid, fchsappid, 'units', NCCHAR, 2, 'A ', iret)

c inboard divertor quantities
      call ncaptc(ncid, ne3dlid, 'long_name', NCCHAR, 34,
     1'electron density, inboard divertor', iret)
      call ncaptc(ncid, ne3dlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3dlid, 'long_name', NCCHAR, 38,
     1'electron temperature, inboard divertor', iret)
      call ncaptc(ncid, te3dlid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3dlid, 'long_name', NCCHAR, 33,
     1'ion temperature, inboard divertor', iret)
      call ncaptc(ncid, ti3dlid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, tp3dlid, 'long_name', NCCHAR, 35,
     1'plate temperature, inboard divertor', iret)
      call ncaptc(ncid, tp3dlid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, po3dlid, 'long_name', NCCHAR, 26,
     1'potential, inboard divertor', iret)
      call ncaptc(ncid, po3dlid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3dlid, 'long_name', NCCHAR, 30,
     1'atom density, inboard divertor', iret)
      call ncaptc(ncid, an3dlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3dlid, 'long_name', NCCHAR, 34,
     1'molecule density, inboard divertor', iret)
      call ncaptc(ncid, mn3dlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, fn3dlid, 'long_name', NCCHAR, 40,
     1'poloidal spec one flux, inboard divertor', iret)
      call ncaptc(ncid, fn3dlid, 'units', NCCHAR, 4, 's^-1', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fn3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fl3dlid, 'long_name', NCCHAR, 40,
     1'poloidal electron flux, inboard divertor', iret)
      call ncaptc(ncid, fl3dlid, 'units', NCCHAR, 4, 's^-1', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fl3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fo3dlid, 'long_name', NCCHAR, 35,
     1'poloidal ion flux, inboard divertor', iret)
      call ncaptc(ncid, fo3dlid, 'units', NCCHAR, 4, 's^-1', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fo3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fe3dlid, 'long_name', NCCHAR, 47,
     1'poloidal electron energy flux, inboard divertor', iret)
      call ncaptc(ncid, fe3dlid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fe3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fi3dlid, 'long_name', NCCHAR, 42,
     1'poloidal ion energy flux, inboard divertor', iret)
      call ncaptc(ncid, fi3dlid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fi3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ft3dlid, 'long_name', NCCHAR, 44,
     1'poloidal total energy flux, inboard divertor', iret)
      call ncaptc(ncid, ft3dlid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, ft3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fc3dlid, 'long_name', NCCHAR, 35,
     1'poloidal current, inboard divertor', iret)
      call ncaptc(ncid, fc3dlid, 'units', NCCHAR, 2, 'A ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fc3dlid, 'scale', NCDOUBLE, 1, doubleval, iret)
c inboard midplane quantities
      call ncaptc(ncid, ne3diid, 'long_name', NCCHAR, 34,
     1'electron density, inboard midplane', iret)
      call ncaptc(ncid, ne3diid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3diid, 'long_name', NCCHAR, 38,
     1'electron temperature, inboard midplane', iret)
      call ncaptc(ncid, te3diid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3diid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3diid, 'long_name', NCCHAR, 33,
     1'ion temperature, inboard midplane', iret)
      call ncaptc(ncid, ti3diid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3diid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, po3diid, 'long_name', NCCHAR, 27,
     1'potential, inboard midplane', iret)
      call ncaptc(ncid, po3diid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3diid, 'long_name', NCCHAR, 30,
     1'atom density, inboard midplane', iret)
      call ncaptc(ncid, an3diid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3diid, 'long_name', NCCHAR, 34,
     1'molecule density, inboard midplane', iret)
      call ncaptc(ncid, mn3diid, 'units', NCCHAR, 4, 'm^-3', iret)

      call ncaptc(ncid, dn3diid, 'long_name', NCCHAR, 39,
     1'diffusion coefficient, inboard midplane', iret)
      call ncaptc(ncid, dn3diid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, dp3diid, 'long_name', NCCHAR, 48,
     1'pressure diffusion coefficient, inboard midplane', iret)
      call ncaptc(ncid, dp3diid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, lh3diid, 'long_name', NCCHAR, 48,
     1'perp neutral heat flux limiter, inboard midplane', iret)
      call ncaptc(ncid, lh3diid, 'units', NCCHAR, 1, ' ', iret)
      call ncaptc(ncid, ln3diid, 'long_name', NCCHAR, 48,
     1'perp neutral part flux limiter, inboard midplane', iret)
      call ncaptc(ncid, ln3diid, 'units', NCCHAR, 1, ' ', iret)
      call ncaptc(ncid, ke3diid, 'long_name', NCCHAR, 46,
     1'electron thermal diffusivity, inboard midplane', iret)
      call ncaptc(ncid, ke3diid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, ki3diid, 'long_name', NCCHAR, 41,
     1'ion thermal diffusivity, inboard midplane', iret)
      call ncaptc(ncid, ki3diid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, vx3diid, 'long_name', NCCHAR, 41,
     1'poloidal pinch velocity, inboard midplane', iret)
      call ncaptc(ncid, vx3diid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vy3diid, 'long_name', NCCHAR, 39,
     1'radial pinch velocity, inboard midplane', iret)
      call ncaptc(ncid, vy3diid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vs3diid, 'long_name', NCCHAR, 39,
     1'viscosity coefficient, inboard midplane', iret)
      call ncaptc(ncid, vs3diid, 'units', NCCHAR, 12, 'm.kg^-1.s^-1',
     1 iret)

c upper inboard divertor quantities
      if(nytl.gt.0) then
      call ncaptc(ncid, ne3dtlid, 'long_name', NCCHAR, 40,
     1'electron density, upper inboard divertor', iret)
      call ncaptc(ncid, ne3dtlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3dtlid, 'long_name', NCCHAR, 44,
     1'electron temperature, upper inboard divertor', iret)
      call ncaptc(ncid, te3dtlid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3dtlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3dtlid, 'long_name', NCCHAR, 39,
     1'ion temperature, upper inboard divertor', iret)
      call ncaptc(ncid, ti3dtlid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3dtlid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, tp3dtlid, 'long_name', NCCHAR, 41,
     1'plate temperature, upper inboard divertor', iret)
      call ncaptc(ncid, tp3dtlid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, po3dtlid, 'long_name', NCCHAR, 32,
     1'potential, upper inboard divertor', iret)
      call ncaptc(ncid, po3dtlid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3dtlid, 'long_name', NCCHAR, 36,
     1'atom density, upper inboard divertor', iret)
      call ncaptc(ncid, an3dtlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3dtlid, 'long_name', NCCHAR, 40,
     1'molecule density, upper inboard divertor', iret)
      call ncaptc(ncid, mn3dtlid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, fn3dtlid, 'long_name', NCCHAR, 46,
     1'poloidal spec one flux, upper inboard divertor', iret)
      call ncaptc(ncid, fn3dtlid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, fl3dtlid, 'long_name', NCCHAR, 46,
     1'poloidal electron flux, inboard divertor', iret)
      call ncaptc(ncid, fl3dtlid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, fo3dtlid, 'long_name', NCCHAR, 41,
     1'poloidal ion flux, upper inboard divertor', iret)
      call ncaptc(ncid, fo3dtlid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, fe3dtlid, 'long_name', NCCHAR, 53,
     1'poloidal electron energy flux, upper inboard divertor', iret)
      call ncaptc(ncid, fe3dtlid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fi3dtlid, 'long_name', NCCHAR, 48,
     1'poloidal ion energy flux, upper inboard divertor', iret)
      call ncaptc(ncid, fi3dtlid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, ft3dtlid, 'long_name', NCCHAR, 50,
     1'poloidal total energy flux, upper inboard divertor', iret)
      call ncaptc(ncid, ft3dtlid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fc3dtlid, 'long_name', NCCHAR, 41,
     1'poloidal current, upper inboard divertor', iret)
      call ncaptc(ncid, fc3dtlid, 'units', NCCHAR, 2, 'A ', iret)
      endif
c outboard midplane quantities
      call ncaptc(ncid, ne3daid, 'long_name', NCCHAR, 35,
     1'electron density, outboard midplane', iret)
      call ncaptc(ncid, ne3daid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3daid, 'long_name', NCCHAR, 39,
     1'electron temperature, outboard midplane', iret)
      call ncaptc(ncid, te3daid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3daid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3daid, 'long_name', NCCHAR, 34,
     1'ion temperature, outboard midplane', iret)
      call ncaptc(ncid, ti3daid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3daid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, po3daid, 'long_name', NCCHAR, 28,
     1'potential, outboard midplane', iret)
      call ncaptc(ncid, po3daid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3daid, 'long_name', NCCHAR, 31,
     1'atom density, outboard midplane', iret)
      call ncaptc(ncid, an3daid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3daid, 'long_name', NCCHAR, 35,
     1'molecule density, outboard midplane', iret)
      call ncaptc(ncid, mn3daid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, dn3daid, 'long_name', NCCHAR, 40,
     1'diffusion coefficient, outboard midplane', iret)
      call ncaptc(ncid, dn3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, dp3daid, 'long_name', NCCHAR, 49,
     1'pressure diffusion coefficient, outboard midplane', iret)
      call ncaptc(ncid, dp3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, lh3daid, 'long_name', NCCHAR, 49,
     1'perp neutral heat flux limiter, outboard midplane', iret)
      call ncaptc(ncid, lh3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, ln3daid, 'long_name', NCCHAR, 49,
     1'perp neutral part flux limiter, outboard midplane', iret)
      call ncaptc(ncid, ln3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, ke3daid, 'long_name', NCCHAR, 47,
     1'electron thermal diffusivity, outboard midplane', iret)
      call ncaptc(ncid, ke3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, ki3daid, 'long_name', NCCHAR, 42,
     1'ion thermal diffusivity, outboard midplane', iret)
      call ncaptc(ncid, ki3daid, 'units', NCCHAR, 8, 'm^2.s^-1', iret)
      call ncaptc(ncid, vx3daid, 'long_name', NCCHAR, 42,
     1'poloidal pinch velocity, outboard midplane', iret)
      call ncaptc(ncid, vx3daid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vy3daid, 'long_name', NCCHAR, 40,
     1'radial pinch velocity, outboard midplane', iret)
      call ncaptc(ncid, vy3daid, 'units', NCCHAR, 6, 'm.s^-1', iret)
      call ncaptc(ncid, vs3daid, 'long_name', NCCHAR, 40,
     1'viscosity coefficient, outboard midplane', iret)
      call ncaptc(ncid, vs3daid, 'units', NCCHAR, 12, 'm.kg^-1.s^-1',
     1 iret)
c outboard divertor quantities
      call ncaptc(ncid, ne3drid, 'long_name', NCCHAR, 35,
     1'electron density, outboard divertor', iret)
      call ncaptc(ncid, ne3drid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3drid, 'long_name', NCCHAR, 39,
     1'electron temperature, outboard divertor', iret)
      call ncaptc(ncid, te3drid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3drid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3drid, 'long_name', NCCHAR, 34,
     1'ion temperature, outboard divertor', iret)
      call ncaptc(ncid, ti3drid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3drid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, tp3drid, 'long_name', NCCHAR, 36,
     1'plate temperature, outboard divertor', iret)
      call ncaptc(ncid, tp3drid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, po3drid, 'long_name', NCCHAR, 28,
     1'potential, outboard divertor', iret)
      call ncaptc(ncid, po3drid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3drid, 'long_name', NCCHAR, 31,
     1'atom density, outboard divertor', iret)
      call ncaptc(ncid, an3drid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3drid, 'long_name', NCCHAR, 35,
     1'molecule density, outboard divertor', iret)
      call ncaptc(ncid, mn3drid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, fn3drid, 'long_name', NCCHAR, 41,
     1'poloidal spec one flux, outboard divertor', iret)
      call ncaptc(ncid, fn3drid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fl3drid, 'long_name', NCCHAR, 41,
     1'poloidal electron flux, outboard divertor', iret)
      call ncaptc(ncid, fl3drid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fo3drid, 'long_name', NCCHAR, 36,
     1'poloidal ion flux, outboard divertor', iret)
      call ncaptc(ncid, fo3drid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fe3drid, 'long_name', NCCHAR, 48,
     1'poloidal electron energy flux, outboard divertor', iret)
      call ncaptc(ncid, fe3drid, 'units', NCCHAR, 4, 's^-1', iret)
      call ncaptc(ncid, fi3drid, 'long_name', NCCHAR, 43,
     1'poloidal ion energy flux, outboard divertor', iret)
      call ncaptc(ncid, fi3drid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, ft3drid, 'long_name', NCCHAR, 45,
     1'poloidal total energy flux, outboard divertor', iret)
      call ncaptc(ncid, ft3drid, 'units', NCCHAR, 2, 'W ', iret)
      call ncaptc(ncid, fc3drid, 'long_name', NCCHAR, 36,
     1'poloidal current, outboard divertor', iret)
      call ncaptc(ncid, fc3drid, 'units', NCCHAR, 2, 'A ', iret)
c upper outboard divertor quantities
      if(nytr.gt.0) then
      call ncaptc(ncid, ne3dtrid, 'long_name', NCCHAR, 41,
     1'electron density, upper outboard divertor', iret)
      call ncaptc(ncid, ne3dtrid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, te3dtrid, 'long_name', NCCHAR, 45,
     1'electron temperature, upper outboard divertor', iret)
      call ncaptc(ncid, te3dtrid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, te3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ti3dtrid, 'long_name', NCCHAR, 40,
     1'ion temperature, upper outboard divertor', iret)
      call ncaptc(ncid, ti3dtrid, 'units', NCCHAR, 2, 'eV', iret)
      doubleval(1) = 1.0_R8/ev
      call ncapt(ncid, ti3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, tp3dtrid, 'long_name', NCCHAR, 42,
     1'plate temperature, upper outboard divertor', iret)
      call ncaptc(ncid, tp3dtrid, 'units', NCCHAR, 2, 'K ', iret)
      call ncaptc(ncid, po3dtrid, 'long_name', NCCHAR, 34,
     1'potential, upper outboard divertor', iret)
      call ncaptc(ncid, po3dtrid, 'units', NCCHAR, 2, 'V ', iret)
      call ncaptc(ncid, an3dtrid, 'long_name', NCCHAR, 37,
     1'atom density, upper outboard divertor', iret)
      call ncaptc(ncid, an3dtrid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, mn3dtrid, 'long_name', NCCHAR, 41,
     1'molecule density, upper outboard divertor', iret)
      call ncaptc(ncid, mn3dtrid, 'units', NCCHAR, 4, 'm^-3', iret)
      call ncaptc(ncid, fn3dtrid, 'long_name', NCCHAR, 47,
     1'poloidal spec one flux, upper outboard divertor', iret)
      call ncaptc(ncid, fn3dtrid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fn3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fl3dtrid, 'long_name', NCCHAR, 47,
     1'poloidal electron flux, upper outboard divertor', iret)
      call ncaptc(ncid, fl3dtrid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fl3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fo3dtrid, 'long_name', NCCHAR, 42,
     1'poloidal ion flux, upper outboard divertor', iret)
      call ncaptc(ncid, fo3dtrid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fo3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fe3dtrid, 'long_name', NCCHAR, 54,
     1'poloidal electron energy flux, upper outboard divertor', iret)
      call ncaptc(ncid, fe3dtrid, 'units', NCCHAR, 4, 's^-1', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fe3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fi3dtrid, 'long_name', NCCHAR, 49,
     1'poloidal ion energy flux, upper outboard divertor', iret)
      call ncaptc(ncid, fi3dtrid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fi3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, ft3dtrid, 'long_name', NCCHAR, 51,
     1'poloidal total energy flux, upper outboard divertor', iret)
      call ncaptc(ncid, ft3dtrid, 'units', NCCHAR, 2, 'W ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, ft3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      call ncaptc(ncid, fc3dtrid, 'long_name', NCCHAR, 42,
     1'poloidal current, upper outboard divertor', iret)
      call ncaptc(ncid, fc3dtrid, 'units', NCCHAR, 2, 'A ', iret)
      doubleval(1) = -1.0_R8
      call ncapt(ncid, fc3dtrid, 'scale', NCDOUBLE, 1, doubleval, iret)
      endif

* leave define mode
      call ncendf (ncid, iret)
      call ncclos (ncid, iret)
      return
      end subroutine b2crtimecdf

#ifdef IBMaix
@process noextchk
#endif
      subroutine rwcdf(rw,ncid,data_name,imap,data_set,iret)
c
      use b2mod_types
#     include <netcdf.inc>
c
      character*(*) rw,data_name
      integer ncid,imap(*),iret,i
      real (kind=R8) ::
     *   data_set(*), rcode
      character*(*) datac(*)
      character*(maxncnam) varnam,dimnam
      integer vartyp,nvdims,vdims(maxvdims),nvatts,
     1 start(maxvdims),count(maxvdims),lmap(maxvdims),
     2 noofbytes,one(maxvdims)
c
      character*(*) timnam
      character*(maxncnam) timsav
      integer ntsav,ntstep
      save timsav,ntsav
      data timsav /'!!!! INVALID NAME !!!!'/
      data one /maxvdims*1/
c
      integer lnblnk
      external lnblnk
      integer data_id
c
c     write(*,*) 'Processing ',data_name,' with ncid ',ncid
      data_id=ncvid(ncid,data_name,iret)
      if(iret.ne.0) then
        write(*,*) data_name,' not declared'
        stop
      else
        call ncvinq(ncid,data_id,varnam,vartyp,nvdims,vdims,
     1   nvatts,rcode)
        if(vartyp.eq.ncbyte.or.vartyp.eq.ncchar) then
          noofbytes=1
        elseif(vartyp.eq.ncshort) then
          noofbytes=2
        elseif(vartyp.eq.nclong.or.vartyp.eq.ncfloat) then
          noofbytes=4
        elseif(vartyp.eq.ncdouble) then
          noofbytes=8
        else
          write(*,*) 'Unknown data type ',vartyp,' for ',data_name
          stop
        endif
        do i=1,nvdims
          call ncdinq(ncid,vdims(i),dimnam,count(i),iret)
c        write(*,*) 'Dimension ',i,' is ',dimnam(1:lnblnk(dimnam)),
c    1    ' with length ',count(i)
          if(dimnam.eq.timsav) then
            start(i)=ntsav
            count(i)=1
c           write(*,*) 'TIME dimension'
          else
            start(i)=1
          endif
          lmap(i)=imap(i)*noofbytes
        enddo
        if(rw.eq.'read') then
          call ncvgtg(ncid,data_id,start,count,one,lmap,data_set,iret)
        elseif(rw.eq.'write') then
c         write(*,*) 'ncid : ', ncid
c         write(*,*) 'data_id :', data_id
c         write(*,*) 'start :', start(1:nvdims)
c         write(*,*) 'count :', count(1:nvdims)
c         write(*,*) 'one :', one(1:nvdims)
c         write(*,*) 'lmap :', lmap(1:nvdims)
c         write(*,*) 'data_set :', data_set(1:nvdims)
          call ncvptg(ncid,data_id,start,count,one,lmap,data_set,iret)
        else
          write(*,*) 'Either "read" or "write" must be chosen, not',
     1     rw
        endif
      endif
      return
c
      entry rwcdfc(rw,ncid,data_name,imap,datac,iret)
c
c      write(*,*) 'Processing ',data_name,' (char)'
      data_id=ncvid(ncid,data_name,iret)
      if(iret.ne.0) then
        write(*,*) data_name,' not declared'
        stop
      else
        call ncvinq(ncid,data_id,varnam,vartyp,nvdims,vdims,
     1   nvatts,rcode)
        if(vartyp.eq.ncbyte.or.vartyp.eq.ncchar) then
          noofbytes=1
        elseif(vartyp.eq.ncshort) then
          noofbytes=2
        elseif(vartyp.eq.nclong.or.vartyp.eq.ncfloat) then
          noofbytes=4
        elseif(vartyp.eq.ncdouble) then
          noofbytes=8
        else
          write(*,*) 'Unknown data type ',vartyp,' for ',data_name
          stop
        endif
        do i=1,nvdims
          call ncdinq(ncid,vdims(i),dimnam,count(i),iret)
c         write(*,*) 'Dimension ',i,' is ',dimnam(1:lnblnk(dimnam)),
c     1    ' with length ',count(i)
          if(dimnam.eq.timsav) then
            start(i)=ntsav
            count(i)=1
c           write(*,*) 'TIME dimension'
          else
            start(i)=1
          endif
          lmap(i)=imap(i)*noofbytes
        enddo
        if(rw.eq.'read') then
          call ncvggc(ncid,data_id,start,count,one,lmap,datac,iret)
        elseif(rw.eq.'write') then
          call ncvpgc(ncid,data_id,start,count,one,lmap,datac,iret)
        else
          write(*,*) 'Either "read" or "write" must be chosen, not',
     1     rw
        endif
      endif
      return
c
      entry rwcdf_settime(timnam,ntstep)
      write(*,*) 'Saving ',timnam(1:lnblnk(timnam)),
     1 ' as the time dimension'
      write(*,*) 'ntstep = ',ntstep
      timsav=timnam
      ntsav=ntstep
      return
      end subroutine rwcdf
#endif
c
      integer function lnblnk(string)
      implicit none
      integer i
      intrinsic len
c
c returns the position of the last non-blank character in "string"
c
      character*(*) string
      do i=len(string),1,-1
        if(string(i:i).ne.' ') then
          lnblnk=i
          return
        endif
      enddo
      lnblnk=1
      return
      end function lnblnk
c
      subroutine output_ds(crx,cry,nx,ny,iref,target_offset,
     1 jsep,iystart,iyend,filename)
      use b2mod_types
      use b2mod_indirect
      implicit none
      integer nx,ny,iref,jsep,iystart,iyend,target_offset
      real (kind=R8) ::
     *      crx(-1:nx,-1:ny,0:3),cry(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *      ds(-1:ny), ds_offset
      character*(*) filename
      integer ix,iy
      external subini, subend, xertst
      intrinsic sqrt
      real (kind=R8) ::
     *   cr,cz

      cr(ix,iy)=0.25_R8*
     1 (crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)+crx(ix,iy,3))
      cz(ix,iy)=0.25_R8*
     1 (cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)+cry(ix,iy,3))

      call subini ('output_ds')
      iystart=-1
      do while (region(iref,iystart,0).eq.0 .and. iystart.lt.ny)
        iystart=iystart+1
      enddo
      call xertst(iystart.le.ny, 'faulty parameter iystart')
      iyend=ny
      do while (region(iref,iyend,0).eq.0 .and. iyend.gt.-1)
        iyend=iyend-1
      enddo
      call xertst(iyend.ge.iystart, 'faulty parameter iyend')
      ds(iystart)=
     1 sqrt((cr(iref+target_offset,iystart)-
     1          0.5_R8*(crx(iref+target_offset,iystart,0)+
     2                  crx(iref+target_offset,iystart,1)))**2+
     2      (cz(iref+target_offset,iystart)-
     3          0.5_R8*(cry(iref+target_offset,iystart,0)+
     4                  cry(iref+target_offset,iystart,1)))**2)
      do iy=iystart+1,iyend
        ds(iy)=ds(iy-1)+
     1   sqrt((cr(iref+target_offset,iy)-
     2         cr(iref+target_offset,iy-1))**2+
     2        (cz(iref+target_offset,iy)-
     3         cz(iref+target_offset,iy-1))**2)
      enddo
      if(iystart.le.jsep.and.iyend.gt.jsep) then
        ds_offset=(ds(jsep)+ds(jsep+1))/2.0_R8
        do iy=iystart,iyend
          ds(iy)=ds(iy)-ds_offset
        enddo
      endif
      open(99,file=filename)
      do iy=iystart,iyend
        write(99,*) ds(iy)
      enddo
      close(99)
      call subend ()
      return
      end subroutine output_ds
