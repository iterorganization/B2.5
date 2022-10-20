      subroutine get_jsep(nx,ny,jxi,jxa,jsep)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx,ny,jxi,jxa,jsep
      integer ix,iy
      intrinsic min
      external xerrab, xertst, ipgeti

      call subini ('get_jsep')
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')
      if(nnreg(0).eq.1) then
        jxi=nx/4
        jxa=nx/2
        jsep=ny/2
      elseif(nnreg(0).eq.4) then
        if (cry(leftcut(1)+1,topcut(1),0).lt.
     .      cry(leftcut(1)+1,0,0).or.
     .     (cry(leftcut(1)+1,topcut(1),0).eq.
     .      cry(leftcut(1)+1,0,0).and.
     .      crx(leftcut(1)+1,topcut(1),0).lt.
     .      crx(leftcut(1)+1,0,0))) then ! LSN topology
          jxi=leftcut(1)+(rightcut(1)-leftcut(1))/4
          jxa=rightcut(1)-(rightcut(1)-leftcut(1))/4
        else ! USN topology
          jxi=rightcut(1)-(rightcut(1)-leftcut(1))/4
          jxa=leftcut(1)+(rightcut(1)-leftcut(1))/4
        end if
        jsep=topcut(1)-1
      elseif(nnreg(0).eq.5) then
        jxi=leftcut(1)+(rightcut(1)-leftcut(1))/4
        jxa=leftcut(1)+(rightcut(1)-leftcut(1))/2
        jsep=topcut(1)-1
      elseif(nnreg(0).eq.7) then
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
        topcut(1)=jsep+1
      else
        call xerrab ('Invalid number of regions!')
      endif
      write(*,'(a,3i5)') 'Default settings for jxi, jxa, jsep : ',
     1                                         jxi, jxa, jsep
      if (nnreg(0).eq.1) call ipgeti ('b2mwti_jsep',jsep)
      call xertst (-1.le.jsep.and.jsep.le.ny,
     2  'faulty parameter jsep')
      call ipgeti ('b2mwti_jxi',jxi)
      call xertst (-1.le.jxi.and.jxi.le.nx,
     2 'faulty internal parameter jxi')
      call ipgeti ('b2mwti_jxa',jxa)
      call xertst (-1.le.jxa.and.jxa.le.nx,
     2 'faulty internal parameter jxa')
      write(*,'(a,3i5)')
     .           'User-corrected settings for jxi, jxa, jsep : ',
     1                                        jxi, jxa, jsep

      call subend ()
      return
      end subroutine get_jsep

      subroutine get_nxt(nx, nxtl, nxtr)
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx, nxtl, nxtr

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
        write(*,*) 'nxtl, nxtr : ', nxtl, nxtr
      endif

      call subend ()
      return
      end subroutine get_nxt
