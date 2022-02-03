
      subroutine clear_extralabel(label)
      use b2mod_types
#include <graphics_vars.h>
      character*(*) label
      extralabel=label
      return
      end

      subroutine clear_xlabel(label)
      use b2mod_types
#include <graphics_vars.h>
      character*(*) label
      xlabel=label
      return
      end

      subroutine clear_ylabel(label)
      use b2mod_types
#include <graphics_vars.h>
      character*(*) label
      ylabel=label
      return
      end
