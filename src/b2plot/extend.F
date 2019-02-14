
      subroutine extend(field,nCv,mpg,ns)
      use b2mod_types
      use b2mod_indirect
      use b2us_map
      integer nCv, ns
      type(mapping) :: mpg
      real (kind=R8) ::
     *   field(nCv,0:ns-1)
      integer ic,ifc,icin
      do is=0,ns-1
        do ic=mpg%nCi+1,nCv 
          ifc = mpg%cvFc(ic)
          icin = min(mpg%fcCv(ifc,1), mpg%fcCv(ifc,2))
          field(ic,is)=field(icin,is)
        enddo
      enddo
      return
      end
