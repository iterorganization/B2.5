      subroutine b2yt_jsep(nx1,ny1,jxi1,jxa1,jsep1,geom_match_dist)
      use b2mod_types
      use b2mod_geo2
      use b2mod_indirect
      implicit none
      integer nx1,ny1,jxi1,jxa1,jsep1
      integer ix,iy
      real (kind=R8) ::
     * geom_match_dist
      intrinsic min
      external subini, xerrab, xertst, ipgetr, subend
      real (kind=R8) ::
     *   dist
      integer ix1,iy1,ip1,ix2,iy2,ip2
      dist(ix1,iy1,ip1,ix2,iy2,ip2)=
     1 sqrt((crx1(ix1,iy1,ip1)-crx1(ix2,iy2,ip2))**2+
     2      (cry1(ix1,iy1,ip1)-cry1(ix2,iy2,ip2))**2)
      logical match
      match(ix1,iy1,ix2,iy2)=
     1 (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist

      call subini ('b2yt_jsep')
      if(nnreg(0).eq.1) then
        jxi1=nx1/4
        jxa1=3*nx1/4
        jsep1=ny1/2
      elseif(nnreg(0).eq.4 .or. nnreg(0).eq.5) then
        jxi1=leftcut1(1)+(rightcut1(1)-leftcut1(1))/4
        jxa1=rightcut1(1)-(rightcut1(1)-leftcut1(1))/4
        jsep1=topcut1(1)-1
      elseif(nnreg(0).eq.8) then
        jxi1=(leftcut1(1)+leftcut1(2))/2
        jxa1=(rightcut1(1)+rightcut1(2))/2
        jsep1=min(topcut1(1),topcut1(2))-1
      elseif(nnreg(0).eq.2) then
        jxi1=nx1/4
        jxa1=nx1/2
        jsep1=-2
        do iy=-1,ny1
          do ix=2,0,-1
            if(match(nx1-ix,iy,ix-1,iy)) then
              jsep1=iy
              exit
            endif
          enddo
        enddo
      else
        call xerrab ('Invalid number of regions!')
      endif
      call xertst (-1.le.jsep1.and.jsep1.le.ny1,
     2  'faulty parameter jsep')
      call xertst (-1.le.jxi1.and.jxi1.le.nx1,
     2 'faulty internal parameter jxi')
      call xertst (-1.le.jxa1.and.jxa1.le.nx1,
     2 'faulty internal parameter jxa')
      write(*,*) 'jxi1, jxa1, jsep1 : ', jxi1, jxa1, jsep1

      call subend ()
      return
      end subroutine b2yt_jsep

      subroutine b2yt_nxt(nx1, nxtl, nxtr)
      use b2mod_geo2
      use b2mod_indirect
      implicit none
      integer nx1, nxtl, nxtr
      external subini, subend

      call subini ('b2yt_nxt')

cxpb  Compute position of additional targets in DN configuration
      nxtl = -2
      nxtr = nx1+1
      if (nnreg(0).eq.8) then
        nxtl = 0
        do while (rightix1(nxtl,max(topcut1(1),topcut1(2))).ne.nx1+1
     &   .and.nxtl.lt.nx1)
          nxtl=nxtl+1
        enddo
        nxtr = nxtl
        do while (leftix1(nxtr,max(topcut1(1),topcut1(2))).ne.-2 .and.
     &            nxtr.lt.nx1)
          nxtr=nxtr+1
        enddo
      endif

      call subend ()
      return
      end subroutine b2yt_nxt

