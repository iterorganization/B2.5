
      subroutine extend(field,mpg,ns)
      use b2mod_types
      use b2mod_indirect
      use b2us_map
      implicit none
      integer, intent(in) :: ns
      type(mapping), intent(in) :: mpg
      real (kind=R8) ::
     *   field(mpg%nCv,0:ns-1)
      integer is,iCv,iFc,iCin

      do is=0,ns-1
        do iCv=mpg%nCi+1,mpg%nCv
          iFc = mpg%cvFc(mpg%cvFcP(iCv,1))
          iCin = min(mpg%fcCv(iFc,1), mpg%fcCv(iFc,2))
          field(iCv,is)=field(iCin,is)
        enddo
      enddo

      return
      end subroutine extend
