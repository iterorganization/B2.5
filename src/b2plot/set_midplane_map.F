c
      subroutine set_midplane_mapping(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        map_to_midplane=.not.map_to_midplane
        write (*,*) 'mapping to midplane option changed'
      else
        map_to_midplane=.false.
        write (*,*) 'Not mapping to midplane'
      endif
      return
      end subroutine set_midplane_mapping
