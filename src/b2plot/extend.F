
      subroutine extend(field,nCv,mpg,ns)
      use b2mod_types
      use b2mod_indirect
      use b2us_map
      implicit none
      integer nCv, ns
      type(mapping) :: mpg
      real (kind=R8) ::
     *   field(nCv,0:ns-1)
      integer is,ic,ifc,icin
      do is=0,ns-1
        do ic=mpg%nCi+1,nCv 
          ifc = mpg%cvFc(mpg%cvFcP(ic,1))
          icin = min(mpg%fcCv(ifc,1), mpg%fcCv(ifc,2))
          field(ic,is)=field(icin,is)
        enddo
      enddo
      return
      end subroutine extend
