c
c
c
      subroutine process(in_nx,in_ny,in_ns,
     1 in_crx,in_cry,in_bb,in_vol,in_hx,in_hy,in_qz,in_gs,
     2 in_za,in_zn,in_am,in_na,
     3 in_ua,in_uadia,in_te,in_ti,in_po,
     4 in_fna,in_fhe,in_fhi,in_fch,in_fchdia,
     3 in_resco,in_reshe,in_reshi,in_resmo,in_resmt,in_respo,
     4 in_nncut,in_leftcut,in_rightcut,in_topcut,in_bottomcut)

      implicit none

      integer in_nx,in_ny,in_ns
      integer in_nncut,
     1 in_leftcut(in_nncut),in_rightcut(in_nncut),
     2 in_topcut(in_nncut),in_bottomcut(in_nncut)
      real*8
     * in_crx(0:in_nx+1,0:in_ny+1,0:3),
     . in_cry(0:in_nx+1,0:in_ny+1,0:3),
     . in_bb(0:in_nx+1,0:in_ny+1,0:3),
     . in_vol(0:in_nx+1,0:in_ny+1),
     . in_hx(0:in_nx+1,0:in_ny+1),
     . in_hy(0:in_nx+1,0:in_ny+1),
     . in_qz(0:in_nx+1,0:in_ny+1,0:1),
     . in_gs(0:in_nx+1,0:in_ny+1,0:2),
     . in_za(1:in_ns),
     . in_zn(1:in_ns),
     . in_am(1:in_ns),
     . in_na(0:in_nx+1,0:in_ny+1,in_ns),
     . in_ua(0:in_nx+1,0:in_ny+1,in_ns),
     . in_uadia(0:in_nx+1,0:in_ny+1,0:1,in_ns),
     . in_te(0:in_nx+1,0:in_ny+1),
     . in_ti(0:in_nx+1,0:in_ny+1),
     . in_po(0:in_nx+1,0:in_ny+1)
      real*8
     . in_fna(0:in_nx+1,0:in_ny+1,0:1,in_ns),
     . in_fhe(0:in_nx+1,0:in_ny+1,0:1),
     . in_fhi(0:in_nx+1,0:in_ny+1,0:1),
     . in_fch(0:in_nx+1,0:in_ny+1,0:1),
     . in_fchdia(0:in_nx+1,0:in_ny+1,0:1),
     . in_resco(0:in_nx+1,0:in_ny+1,in_ns),
     . in_resmo(0:in_nx+1,0:in_ny+1,in_ns),
     . in_reshe(0:in_nx+1,0:in_ny+1),
     . in_reshi(0:in_nx+1,0:in_ny+1),
     . in_resmt(0:in_nx+1,0:in_ny+1),
     . in_respo(0:in_nx+1,0:in_ny+1)

#ifdef WANT_PROCESS
#     include <COUPLE/EIRDIAG.F>
#     include <COMDIM.F>
#     include <COMLAB.F>
#     include <COMGEO.F>
#     include <COMPLASMA.F>
#     include <COMEQ.F>
#     include <COMDIAG.F>
#     include <COMNUM.F>
#     include <COMINIT.F>
#     include <COMBOUND.F>
#     include <CCUT.F>
#     include <IMPDIM.F>
#     include <COMNEUT.F>
#     include <COMTRA.F>
#     include <COMILU.F>
#     include <COMELM.F>
#     include <COMUSER.F>
#     include <COMFEEDBACK.F>

      integer ix,iy,is,in

      real*8
     *   b2ssum,sdot,b2s2sum
      external b2ssum,sdot,b2s2sum

      real*8
     *   pitmin,pitmax

      nx=in_nx
      ny=in_ny
      ns=in_ns
      if(in_topcut(1).lt.-1) then
        nncut=0
        jxi=3*nx/8
        jxa=5*nx/8
        jsep=0
        write(*,*) 'jxi,jxa,jsep ',jxi,jxa,jsep
      else
        nncut=2*in_nncut
        nniso=in_nncut-1
        if(in_nncut.eq.1) then
          jxi=in_leftcut(1)+(in_rightcut(1)-in_leftcut(1))/4
          jxa=in_rightcut(1)-(in_rightcut(1)-in_leftcut(1))/4
          jsep=in_topcut(1)
          write(*,*) 'jxi,jxa,jsep ',jxi,jxa,jsep
          nxcut1(1)=in_leftcut(1)
          nxcut1(2)=in_rightcut(1)
          nxcut2(1)=in_rightcut(1)+1
          nxcut2(2)=in_leftcut(1)+1
          nycut1(1)=in_bottomcut(1)+1
          nycut1(2)=in_bottomcut(1)+1
          nycut2(1)=in_topcut(1)
          nycut2(2)=in_topcut(1)
        elseif(in_nncut.eq.2) then
          jxi=(in_rightcut(1)+in_rightcut(2))/4
          jxa=(in_leftcut(1)+in_leftcut(2))/4
          jsep=min(in_topcut(1),in_topcut(2))
          write(*,*) 'jxi,jxa,jsep ',jxi,jxa,jsep
          nxcut1(1)=in_leftcut(1)
          nxcut1(2)=in_leftcut(2)
          nxcut1(3)=in_rightcut(2)
          nxcut1(4)=in_rightcut(1)
          nxcut2(1)=in_rightcut(1)+1
          nxcut2(2)=in_rightcut(2)+1
          nxcut2(3)=in_leftcut(2)+1
          nxcut2(4)=in_leftcut(1)+1
          nycut1(1)=in_bottomcut(1)+1
          nycut1(2)=in_bottomcut(2)+1
          nycut1(3)=in_bottomcut(2)+1
          nycut1(4)=in_bottomcut(1)+1
          nycut2(1)=in_topcut(1)
          nycut2(2)=in_topcut(2)
          nycut2(3)=in_topcut(2)
          nycut2(4)=in_topcut(1)
          ix=0
          do while (rightix(ix,in_ny/2).ne.in_nx+1 .and. ix.lt.in_nx)
            ix=ix+1
          enddo
          if(ix.ne.in_nx) nxiso1(1)=ix+1
          do while (leftix(ix,in_ny/2).ne.-2 .and. ix.lt.in_nx)
            ix=ix+1
          enddo
          if(ix.ne.in_nx) nxiso2(1)=ix
          nyiso1(1)=0
          nyiso2(1)=in_ny
        endif
      endif
      write(*,*) 'nncut ',nncut
c cell centred quantities
      do iy=0,in_ny+1
        do ix=0,in_nx+1
          br(ix,iy,1)=in_crx(ix,iy,1)
          br(ix,iy,2)=in_crx(ix,iy,3)
          br(ix,iy,3)=in_crx(ix,iy,2)
          br(ix,iy,4)=in_crx(ix,iy,0)
          bz(ix,iy,1)=in_cry(ix,iy,1)
          bz(ix,iy,2)=in_cry(ix,iy,3)
          bz(ix,iy,3)=in_cry(ix,iy,2)
          bz(ix,iy,4)=in_cry(ix,iy,0)
          cr(ix,iy)=0.25*
     1     (br(ix,iy,1)+br(ix,iy,2)+br(ix,iy,3)+br(ix,iy,4))
          cz(ix,iy)=0.25*
     1     (bz(ix,iy,1)+bz(ix,iy,2)+bz(ix,iy,3)+bz(ix,iy,4))
          bb(ix,iy)=in_bb(ix,iy,3)
          pit(ix,iy)=in_bb(ix,iy,0)/in_bb(ix,iy,3)
          vol(ix,iy)=in_vol(ix,iy)
          dx(ix,iy)=in_hx(ix,iy)
          dy(ix,iy)=in_hy(ix,iy)
          qcos(ix,iy,0)=in_qz(ix,iy,0)
          qcos(ix,iy,1)=in_qz(ix,iy,1)
          te(ix,iy)=in_te(ix,iy)
          ti(ix,iy)=in_ti(ix,iy)
          phi(ix,iy)=in_po(ix,iy)
          resee(ix,iy)=in_reshe(ix,iy)
          resei(ix,iy)=in_reshi(ix,iy)
          respo(ix,iy)=in_respo(ix,iy)
        enddo
      enddo
      do iy=0,nyd+1
        do ix=0,nxd+1
          sx(ix,iy)=0
          feex(ix,iy)=0
          feix(ix,iy)=0
          fchx(ix,iy)=0
          sy(ix,iy)=0
          feey(ix,iy)=0
          feiy(ix,iy)=0
        enddo
      enddo
      do is=1,nsd
        do iy=0,nyd+1
          do ix=0,nxd+1
            up(ix,iy,is)=0
            fnix(ix,iy,is)=0
            fniy(ix,iy,is)=0
          enddo
        enddo
      enddo
c x-direction cell face quantities
      do iy=0,in_ny+1
        do ix=0,in_nx
          sx(ix,iy)=in_gs(ix+1,iy,0)*
     .     (in_qz(leftix(ix+1),leftiy(iy),1)*in_vol(ix+1,iy)+
     .      in_qz(ix+1,iy,1)*in_vol(leftix(ix+1),leftiy(iy)))/
     .     (in_vol(ix+1,iy)+in_vol(leftix(ix+1),leftiy(iy)))
          feex(ix,iy)=in_fhe(ix+1,iy,0)
          feix(ix,iy)=in_fhi(ix+1,iy,0)
          fchx(ix,iy)=in_fch(ix+1,iy,0)
          cd(ix,iy)=in_fchdia(ix+1,iy,0)
        enddo
      enddo
c y-direction cell face quantities
      do iy=0,in_ny
        do ix=0,in_nx+1
          sy(ix,iy)=in_gs(ix,iy+1,1)
          feey(ix,iy)=in_fhe(ix,iy+1,1)
          feiy(ix,iy)=in_fhi(ix,iy+1,1)
          fchy(ix,iy)=in_fch(ix,iy+1,1)
          dd(ix,iy)=in_fchdia(ix,iy+1,1)
        enddo
      enddo
      do is=1,in_ns
        zi(is)=in_za(is)
        mi(is)=in_am(is)
c cell centred quantities
        do iy=0,in_ny+1
          do ix=0,in_nx+1
            ni(ix,iy,is)=in_na(ix,iy,is)
            resco(ix,iy,is)=in_resco(ix,iy,is)
            resmo(ix,iy,is)=in_resmo(ix,iy,is)
          enddo
        enddo
c x-direction cell face quantities
        do iy=0,in_ny+1
          do ix=0,in_nx
            up(ix,iy,is)=in_ua(ix+1,iy,is)
            fnix(ix,iy,is)=in_fna(ix+1,iy,0,is)
            ud(ix,iy,is)=in_uadia(ix+1,iy,0,is)
          enddo
        enddo
c y-direction cell face quantities
        do iy=0,in_ny
          do ix=0,in_nx+1
            fniy(ix,iy,is)=in_fna(ix,iy+1,1,is)
            vd(ix,iy,is)=in_uadia(ix,iy+1,1,is)
          enddo
        enddo
      enddo

      do iy=0,in_ny+1
        do ix=0,in_nx+1
          ne(ix,iy)=0
          pr(ix,iy)=0
          zeff(ix,iy)=0
          do is=1,ns
            ne(ix,iy)=ne(ix,iy)+ni(ix,iy,is)*zi(is)
            pr(ix,iy)=pr(ix,iy)+ni(ix,iy,is)*ti(ix,iy)
            zeff(ix,iy)=zeff(ix,iy)+zi(is)**2*ni(ix,iy,is)
          enddo
          pr(ix,iy)=pr(ix,iy)+ne(ix,iy)*te(ix,iy)
          zeff(ix,iy)=zeff(ix,iy)/ne(ix,iy)
        enddo
      enddo

*  -- compute 0d reductions --
*     -- total volume --
         ff(1) = b2s2sum (nx+2, ny+2, vol, 1, nxd+2)
*     -- global dx --
c         ff(2) = b2ssum (nx+2, xx(0,4), 1)
*     -- global dy --
c         ff(3) = b2ssum (ny+2, yy(0,4), 1)
*     -- global sx --
c         ff(4) = ff(1) * ff(2)
*     -- global sy --
c         ff(5) = ff(1) * ff(3)
*     -- global pit --
c         ff(6) = sdot (ny+2, yy(0,6), 1, yy(0,3), 1) /
c     .         b2ssum (ny+2, yy(0,3), 1)
         do in=2,6
           ff(in) = 1.0
         enddo
*  -- x-coordinate labels --
*     cell centers in xx(0:nx+1,1), cell faces in xx(0:nx,2);
       do ix = 0,nx+1
        xx(ix,2) = (real(ix)/nx)
        xx(ix,1) = ((ix-0.5)/nx)
      enddo
*  -- y-coordinate labels --
*     cell centers in yy(0:ny+1,1), cell faces in yy(0:ny,2);
      do iy = 0, ny+1
        yy(iy,2) = (real(iy)/ny)
        yy(iy,1) = ((iy-0.5)/ny)
      enddo
      output_format=1
c      output_format=4
      version=newversion
      pitmin=pit(0,0)
      pitmax=pit(0,0)
      do iy=0,ny+1
        do ix=0,nx+1
          pitmin=min(pitmin,pit(ix,iy))
          pitmax=max(pitmax,pit(ix,iy))
        enddo
      enddo
      write(*,*) 'min(pit),max(pit)',pitmin,pitmax
#ifndef NO_CDF
      call b2rw(.false.)
#else
      write(*,*) 'NO_CDF defined, no reading possible!'
#endif

#endif
      return
      end subroutine process
