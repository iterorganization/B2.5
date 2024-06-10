      subroutine int2d(vals,tdim,ndim,temp,ntemp,dens,ndens,
     1 rtemp,rdens,rval)
      use b2mod_types
      implicit none
      integer tdim,ndim,ntemp,ndens
      real (kind=R8) ::
     *   vals(1:tdim,1:ndim),temp(1:ntemp),dens(1:ndens),
     1   rtemp,rdens,rval
      real (kind=R8) ::
     *   ftemp,fdens
      integer it,in
      integer ipos
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
c      write(*,*) it,in,rval

      return
      end subroutine int2d

