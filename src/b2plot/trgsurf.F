c
      subroutine trgsurf(nCv,igraphs)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      integer nCv,lnblnk,igraphs
      external lnblnk

      if(b2plot_debug.ge.10) 
     > write(*,*) 'Plotting target array: nbnd, ndepth_nml = ',
     1 nbnd,ndepth_nml
      if(nbnd .lt. 1 .or. ndepth_nml .lt. 1) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1     'must be specified'
        return
      endif

      if(typeofplot .ge. 2) then
        if(b2plot_debug.ge.7) write(*,*) 'computational domain'
        call init_page(
     1     target_label(itwrk)(1:lnblnk(target_label(itwrk))),
     1     values(24),values(25),values(26),values(27),1,.false.,.true.,
     2     .false.)
        call plottarget(nCv,nbnd,ndepth_nml,
     .     targetcr(1:nbnd,1:ndepth_nml,1:2),
     .     targetcz(1:nbnd,1:ndepth_nml,1:2),
     .     target_stack(1:nbnd,1:ndepth_nml,itwrk),
     4     values(32),values(33),.false.,trglinlog,igraphs)
        call end_page
      endif
      if(mod(typeofplot,2) .eq. 1) then
        if(b2plot_debug.ge.7) write(*,*) 'physical domain'
        call init_page(
     1     target_label(itwrk)(1:lnblnk(target_label(itwrk))),
     1     values(28),values(29),values(30),values(31),1,.true.,.true.,
     2     .true.)
        call plottarget(nCv,nbnd,ndepth_nml,
     2     targetcr(1:nbnd,1:ndepth_nml,1:2),
     .     targetcz(1:nbnd,1:ndepth_nml,1:2),
     3     target_stack(1:nbnd,1:ndepth_nml,itwrk),
     4     values(32),values(33),.true.,trglinlog,igraphs)
        call end_page
      endif
      return
      end subroutine trgsurf
