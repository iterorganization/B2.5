c
      subroutine cont
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      implicit none
      character*32 tmpstring
      integer is,lnblnk
      external lnblnk
      if(b2plot_debug.ge.10) write(*,*) 'Plotting array: nx, ny, nz = ',
     1 ndim(1,iwrk),ndim(2,iwrk),ndim(3,iwrk)
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
      call map(nx,ny,vol,cr,work(-1,-1,ns+1,iwrk))
      call map(nx,ny,vol,cz,work(-1,-1,ns+2,iwrk))
      do is=0,ndim(3,iwrk)-1
        tmpstring=textwk(is,iwrk)
        if(typeofplot .ge. 2) then
          if(b2plot_debug.ge.7) write(*,*) 'computational domain'
          call init_page(
     1	   graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//tmpstring,
     1	   values(1),values(2),values(3),values(4),1,.false.,.true.,
     3     .false.)
          call plotcontour(nx,ny,nx,ny,is,
     2	   xx,yy,
     3	   work(-1,-1,is,iwrk),
     4	   values(9),values(10),.false.,linlog)
          call end_page
        endif
        if(mod(typeofplot,2) .eq. 1) then
          if(b2plot_debug.ge.7) write(*,*) 'physical domain'
          call map(nx,ny,
     1         vol,work(-1,-1,is,iwrk),work(-1,-1,ns+3,iwrk))
          call init_page(
     1	   graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//tmpstring,
     1	   values(5),values(6),values(7),values(8),1,.true.,.true.,
     2     .true.)
          call plotcontour(nx,ny,nx,nx,is,
c     2      work(-1,-1,ns+1,iwrk),work(-1,-1,ns+2,iwrk),
     2	   crx(-1,-1,3),cry(-1,-1,3),
     3	   work(-1,-1,ns+3,iwrk),
     4	   values(9),values(10),.true.,linlog)
          if(arrowl .gt. 0.0_R8) call plot_arrows(is)
          call end_page
        endif
      enddo
c      call gridon(.false.)
c      call separatrixon(.false.)
c      call vesselon(.false.)
      return
      end
