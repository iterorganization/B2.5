
      subroutine extend(field,mpg,ns)
      use b2mod_types
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

#ifdef B25_EIRENE
      subroutine eirene_extend(field,ntrii,ncltal,ns)
cxpb  Routine to map Eirene data defined
cxpb  on the control volumes and void triangles
cxpb  to all the Eirene triangles
      use b2mod_types
      implicit none
      integer, intent(in) :: ntrii, ns
      integer, intent(in) :: ncltal(ntrii)
      real (kind=R8) ::
     *   field(ns,ntrii), field_copy(ns,ntrii)
      integer i,iCv

      field_copy = 0.0_R8
      do i = 1, ntrii
        field_copy(1:ns,i) = field(1:ns,ncltal(i))
      end do
      field = field_copy

      return
      end subroutine eirene_extend
#endif
