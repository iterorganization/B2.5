
      subroutine clear_extralabel(label)
      use b2mod_types
      use b2mod_graphics
      implicit none
      character*(*) label
      extralabel=label
      return
      end subroutine clear_extralabel

      subroutine clear_xlabel(label)
      use b2mod_types
      use b2mod_graphics
      implicit none
      character*(*) label
      xlabel=label
      return
      end subroutine clear_xlabel

      subroutine clear_ylabel(label)
      use b2mod_types
      use b2mod_graphics
      implicit none
      character*(*) label
      ylabel=label
      return
      end subroutine clear_ylabel
