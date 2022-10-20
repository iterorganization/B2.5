      function ipos(vals,nvals,val)
      use b2mod_types
      implicit none
      integer ipos,nvals
      real (kind=R8) ::
     *   vals(1:nvals),val
      integer is,ie,im
      if(val.lt.vals(1)) then
        ipos=0
      elseif(val.gt.vals(nvals)) then
        ipos=nvals+1
      else
        is=1
        ie=nvals
   10   im=(is+ie)/2
        if(val.lt.vals(im)) then
          ie=im
        elseif(val.gt.vals(im)) then
          is=im
        else
          ipos=im
          return
        endif
        if((ie-is).gt.1) goto 10
        ipos=is
      endif
      return
      end function ipos
