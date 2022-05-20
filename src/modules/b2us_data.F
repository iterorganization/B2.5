      module b2us_data
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_switches
      implicit none
c**********************************************************
      type(mapping), save    :: mpg
      type(geometry), save   :: geo
      type(B2State), save    :: state
      type(B2StateExt), save :: state_ext
      type(switches), save   :: switch
c**********************************************************
      end module b2us_data
