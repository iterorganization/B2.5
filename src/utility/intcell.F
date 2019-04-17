      subroutine intcell1(nFc,nCv,nz,mpg,wght,face,centre)
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv,nFc,nz
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   wght(mpg%nCmxFc),face(nFc,0:nz),centre(nCv,0:nz)
      integer is
      do is=0,nz
        call intcell(nFc,nCv,mpg,wght,face(1,is),centre(1,is))
      enddo
      return
      end subroutine intcell1
c
      subroutine intcell(nFc,nCv,mpg,wght,face,centre)
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv,nFc
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   wght(mpg%nCmxFc),face(nFc),centre(nCv)
  
      integer iFc,iCv
      real (kind=R8) :: wghtsum

      centre=0.0_R8
      do iCv = 1, nCv
        wghtsum = 0.0_R8
        do iFc = 1,mpg%cvFcP(iCv,2)
          wghtsum = wghtsum + wght(mpg%cvFcP(iCv,1) + iFc - 1)
          centre(iCv) = centre(iCv) + 
     &      face(mpg%cvFc(mpg%cvFcP(iCv,1) + iFc - 1))*
     &      wght(mpg%cvFcP(iCv,1) + iFc - 1)
        enddo
        centre(iCv) = centre(iCv)/wghtsum
      enddo

      return
      end
