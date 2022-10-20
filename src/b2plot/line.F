c
      subroutine line(idim)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer idim
      nidim=idim
      if(idim .eq. 1) then
        call init_page(graphlabel(iwrk),
     1   values(1),values(2),values(9),values(10),linlog,
     2   .false.,.false.,.false.)
      else if(idim .eq. 2) then
        call init_page(graphlabel(iwrk),
     1   values(3),values(4),values(9),values(10),linlog,
     2   .false.,.false.,.false.)
      else
        write(*,*) 'Either x or y must be chosen'
        nidim=0
      endif
      return
      end subroutine line
