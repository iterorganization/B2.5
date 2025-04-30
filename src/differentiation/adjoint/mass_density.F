c....................................................................
c                      function mass_density
c  function mass_density(ibnd,i) gives the mass density of the target
c  element ibnd at depth i.

      function mass_density(ibnd,i)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer ibnd, i
      real (kind=R8) :: mass_density

      if (i.gt.ndepth_nml) then
        mass_density = matbulk_target_density(ibnd)
      else if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .         target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        mass_density = matbulk_target_density(ibnd)
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        mass_density = coating_target_density(ibnd)
      else
        mass_density = layer_mass_density(ibnd)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          mass_density = mass_density+
     .      matbulk_target_density(ibnd)*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          mass_density = mass_density+
     .      coating_target_density(ibnd)*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function mass_density
