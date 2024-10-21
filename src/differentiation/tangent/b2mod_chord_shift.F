      module b2mod_chord_shift
      use b2mod_types
      implicit none

      real (kind=R8) ::
     *   rshift, zshift, r0, z0, a0, b0, tt0

      contains

      subroutine init_chord_shift
      implicit none

      rshift = 0.0_R8
      zshift = 0.0_R8
      r0 = 0.0_R8
      z0 = 0.0_R8
      a0 = 0.0_R8
      b0 = 0.0_R8
      tt0 = 0.0_R8

      return
      end subroutine init_chord_shift

      end module b2mod_chord_shift
