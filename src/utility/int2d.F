      subroutine int2d(vals,mode,tdim,ndim,temp,ntemp,dens,ndens,
     1 rtemp,rdens,rval)
      use b2mod_types
      implicit none
      integer mode,tdim,ndim,ntemp,ndens
      real (kind=R8) ::
     *   vals(1:tdim,1:ndim),temp(1:ntemp),dens(1:ndens),
     1   rtemp,rdens,rval,rval_edge
      real (kind=R8) ::
     *   ftemp,fdens,fratio
      integer it,in
      integer ipos
      logical extrap

c Three extrapolation modes are possible:
c mode.eq.0: simple 2-D linear extrapolation in log-log space
c            (may lead to spurious results)
c mode.eq.1: use the extrapolated result, but growing no faster
c            than the parameter variation to the 3/2 power (default)
c mode.eq.2: use the nearest domain edge value

      it=ipos(temp,ntemp,rtemp)
      in=ipos(dens,ndens,rdens)
      it=max(1,min(it,ntemp-1))
      in=max(1,min(in,ndens-1))
      ftemp=(rtemp-temp(it))/(temp(it+1)-temp(it))
      fdens=(rdens-dens(in))/(dens(in+1)-dens(in))
      rval=(1.0_R8-ftemp)*(1.0_R8-fdens)*vals(it,in)+
     1             ftemp *(1.0_R8-fdens)*vals(it+1,in)+
     2             ftemp *        fdens *vals(it+1,in+1)+
     3     (1.0_R8-ftemp)*        fdens *vals(it,in+1)
      if (mode.gt.0) then
        extrap = ftemp.lt.0.0_R8 .or. ftemp.gt.1.0_R8 .or.
     &           fdens.lt.0.0_R8 .or. fdens.gt.1.0_R8
        if (extrap) then
          fratio = 0.0_R8
          if (ftemp.lt.0.0_R8) fratio = fratio - ftemp
          if (ftemp.gt.1.0_R8) fratio = fratio +(ftemp-1.0_R8)
          if (fdens.lt.0.0_R8) fratio = fratio - fdens
          if (fdens.gt.1.0_R8) fratio = fratio +(fdens-1.0_R8)
          ftemp = min(max(ftemp,0.0_R8),1.0_R8)
          fdens = min(max(fdens,0.0_R8),1.0_R8)
          rval_edge = (1.0_R8-ftemp)*(1.0_R8-fdens)*vals(it,in)+
     1                        ftemp *(1.0_R8-fdens)*vals(it+1,in)+
     2                        ftemp *        fdens *vals(it+1,in+1)+
     3                (1.0_R8-ftemp)*        fdens *vals(it,in+1)
          if (mode.eq.1) then
            rval = min(rval,rval_edge+1.5_R8*fratio)
          else if (mode.eq.2) then
            rval = rval_edge
          end if
        end if
      endif
c      write(*,*) it,in,rval

      return
      end subroutine int2d
