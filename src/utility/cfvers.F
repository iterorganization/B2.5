      function get_B25_hash()
      implicit none
      character*32 get_B25_hash
#     include <git_version_B25.h>

      get_B25_hash = git_version_B25

      return
      end function get_B25_hash

      function get_ADAS_hash()
      implicit none
      character*32 get_ADAS_hash
#     include <git_version_B25.h>

      get_ADAS_hash = git_version_ADAS

      return
      end function get_ADAS_hash

      function get_SOLPS_hash()
      implicit none
      character*32 get_SOLPS_hash
#     include <git_version_B25.h>

      get_SOLPS_hash = git_version_SOLPS

      return
      end function get_SOLPS_hash
