
      subroutine set_same_page
      use b2mod_types
      use b2mod_graphics
      implicit none

      if(same_page.ne.0) then
        same_page=0
        call end_page
      endif
      same_page=2
      write(*,*) 'same_page set to ',same_page
      return
      end subroutine set_same_page
