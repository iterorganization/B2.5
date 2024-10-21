
      subroutine b2page(ix,iy)
      use b2mod_types
      use b2mod_graphics
      implicit none
      integer ix,iy
      nx_page=ix
      ny_page=iy
      ix_page=1
      iy_page=1
      new_page=.true.
      return
      end subroutine b2page
