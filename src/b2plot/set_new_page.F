      subroutine set_new_page
      use b2mod_types
#include <graphics_vars.h>
      same_page=0
      write(*,*) 'same_page set to ',same_page
      call end_page
      return
      end
