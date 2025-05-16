      subroutine b2yt_jsep(nx1,ny1,jxi1,jxa1,jsep1)
      use b2mod_types
      use b2mod_geo
     , ,only: geom_match_dist
      use b2mod_geo2
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx1,ny1,jxi1,jxa1,jsep1
      integer ix,iy
      intrinsic min
      external xerrab, xertst, ipgetr

      call subini ('b2yt_jsep')
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')
      if(nnreg(0).eq.1) then
        jxi1=nx1/4
        jxa1=nx1/2
        jsep1=ny1/2
        call ipgeti ('b2ytdr_jsep1',jsep1)
      elseif(nnreg(0).eq.4) then
        if (cry1(leftcut1(1)+1,topcut1(1),0).lt.
     .      cry1(leftcut1(1)+1,0,0).or.
     .     (cry1(leftcut1(1)+1,topcut1(1),0).eq.
     .      cry1(leftcut1(1)+1,0,0).and.
     .      crx1(leftcut1(1)+1,topcut1(1),0).lt.
     .      crx1(leftcut1(1)+1,0,0))) then ! LSN topology
          jxi1=leftcut1(1)+(rightcut1(1)-leftcut1(1))/4
          jxa1=rightcut1(1)-(rightcut1(1)-leftcut1(1))/4
        else ! USN topology
          jxi1=rightcut1(1)-(rightcut1(1)-leftcut1(1))/4
          jxa1=leftcut1(1)+(rightcut1(1)-leftcut1(1))/4
        end if
        jsep1=topcut1(1)-1
      elseif(nnreg(0).eq.5) then
        jxi1=leftcut1(1)+(rightcut1(1)-leftcut1(1))/4
        jxa1=leftcut1(1)+(rightcut1(1)-leftcut1(1))/2
        jsep1=topcut1(1)-1
      elseif(nnreg(0).eq.7) then
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
        topcut1(1)=jsep1+1
      else
        call xerrab ('Invalid number of regions!')
      endif
      write(*,*) 'jxi1, jxa1, jsep1 : ', jxi1, jxa1, jsep1
      call xertst (-1.le.jsep1.and.jsep1.le.ny1,
     2  'faulty parameter jsep1')
      call xertst (-1.le.jxi1.and.jxi1.le.nx1,
     2 'faulty internal parameter jxi1')
      call xertst (-1.le.jxa1.and.jxa1.le.nx1,
     2 'faulty internal parameter jxa1')

      call subend ()
      return

      contains

      logical function match(ix1,iy1,ix2,iy2)
      use b2mod_geo2
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      match=
     1 (distc1(ix1,iy1,1,ix2,iy2,0)+distc1(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist
      return
      end function match

      end subroutine b2yt_jsep

      subroutine b2yt_nxt(nx1, nxtl, nxtr)
      use b2mod_geo2
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx1, nxtl, nxtr

      call subini ('b2yt_nxt')

cxpb  Compute position of additional targets in DN configuration
      nxtl = -2
      nxtr = nx1+1
      if (nnreg(0).ge.7) then
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
