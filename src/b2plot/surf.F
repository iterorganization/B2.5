c
      subroutine surf(igraphs)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      character*256 tmpstring,titlestring
      integer is,lnblnk,igraphs
      external lnblnk

      if(b2plot_debug.ge.10) write(*,*) 'Plotting array: nx, ny, nz = ',
     1 ndim(1,iwrk),ndim(2,iwrk),ndim(3,iwrk)
      if(ndim(1,iwrk) .lt. 0 .or. ndim(2,iwrk) .lt. 0) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1     'must be specified'
        return
      endif

      do is=0,ndim(3,iwrk)-1
        tmpstring=graphlabel(iwrk)(1:lnblnk(graphlabel(iwrk)))//
     1   textwk(is,iwrk)
        call customise_label(tmpstring,titlestring)
        if(typeofplot .ge. 2) then
          if(b2plot_debug.ge.7) write(*,*) 'computational domain'
          call init_page(titlestring,
     1     values(1),values(2),values(3),values(4),1,.false.,.true.,
     2     .false.)
          call plotpatch(NX,NY,is,
     .     crx,cry,region,work(:,:,is,iwrk),
     4     values(9),values(10),.false.,linlog,igraphs)
          call end_page
        endif
        if(mod(typeofplot,2) .eq. 1) then
          if(b2plot_debug.ge.7) write(*,*) 'physical domain'
          call init_page(titlestring,
     1     values(5),values(6),values(7),values(8),1,.true.,.true.,
     2     .true.)
          call plotpatch(nx,ny,is,
     2     crx,cry,region,
     3     work(:,:,is,iwrk),
     4     values(9),values(10),.true.,linlog,igraphs)
          if(arrowl .gt. 0.0_R8) call plot_arrows(is)
          call end_page
        endif
      enddo
      return
      end
