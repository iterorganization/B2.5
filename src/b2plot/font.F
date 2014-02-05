c
      subroutine font(charid)
      use b2mod_types
      integer charid
#ifndef PLPLOT
      CALL PCSETI('FN',charid)  ! =====>  set font to ...
      CALL PCSETI('CC(0)',1)    ! =====>  in black
#endif
      return
      end
