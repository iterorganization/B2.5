      subroutine get_jsep(nx,ny,jxi,jxa,jsep)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      implicit none
      integer nx,ny,jxi,jxa,jsep
      integer ix,iy
      real (kind=R8) ::
     * geom_match_dist
      intrinsic min
      external subini, xerrab, xertst, ipgeti, subend
      data geom_match_dist/1.0e-6_R8/
      real (kind=R8) ::
     *   dist
      integer ix1,iy1,ip1,ix2,iy2,ip2
      dist(ix1,iy1,ip1,ix2,iy2,ip2)=
     1 sqrt((crx(ix1,iy1,ip1)-crx(ix2,iy2,ip2))**2+
     2      (cry(ix1,iy1,ip1)-cry(ix2,iy2,ip2))**2)
      logical match
      match(ix1,iy1,ix2,iy2)=
     1 (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist

      call subini ('get_jsep')
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')
      if(nnreg(0).eq.1) then
        jxi=nx/4
        jxa=3*nx/4
        jsep=ny/2
        call ipgeti ('b2mwti_jsep',jsep)
      elseif(nnreg(0).eq.4 .or. nnreg(0).eq.5) then
        jxi=leftcut(1)+(rightcut(1)-leftcut(1))/4
        jxa=rightcut(1)-(rightcut(1)-leftcut(1))/4
        jsep=topcut(1)-1
      elseif(nnreg(0).eq.8) then
        jxi=(leftcut(1)+leftcut(2))/2
        jxa=(rightcut(1)+rightcut(2))/2
        jsep=min(topcut(1),topcut(2))-1
      elseif(nnreg(0).eq.2) then
        jxi=nx/4
        jxa=nx/2
        jsep=-2
        do iy=-1,ny
          do ix=2,0,-1
            if(match(nx-ix,iy,ix-1,iy)) then
              jsep=iy
              exit
            endif
          enddo
        enddo
      else
        call xerrab ('Invalid number of regions!')
      endif
      call xertst (-1.le.jsep.and.jsep.le.ny,
     2  'faulty parameter jsep')
      call ipgeti ('b2mwti_jxi',jxi)
      call xertst (-1.le.jxi.and.jxi.le.nx,
     2 'faulty internal parameter jxi')
      call ipgeti ('b2mwti_jxa',jxa)
      call xertst (-1.le.jxa.and.jxa.le.nx,
     2 'faulty internal parameter jxa')
      write(*,*) 'jxi, jxa, jsep : ', jxi, jxa, jsep

      call subend ()
      return
      end subroutine get_jsep

      subroutine get_nxt(nx, nxtl, nxtr)
      use b2mod_indirect
      implicit none
      integer nx, nxtl, nxtr
      external subini, subend

      call subini ('get_nxt')

cxpb  Compute position of additional targets in DN configuration
      nxtl = -2
      nxtr = nx+1
      if (nnreg(0).eq.8) then
        nxtl = 0
        do while (rightix(nxtl,max(topcut(1),topcut(2))).ne.nx+1
     &   .and.nxtl.lt.nx)
          nxtl=nxtl+1
        enddo
        nxtr = nxtl
        do while (leftix(nxtr,max(topcut(1),topcut(2))).ne.-2 .and.
     &            nxtr.lt.nx)
          nxtr=nxtr+1
        enddo
      endif

      call subend ()
      return
      end subroutine get_nxt
