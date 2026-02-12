
      subroutine b2uxus (nCv, mpg, aa, itcnt, resfun, corfun, name, style) 
      use b2mod_types
      use b2us_map_diffv
      implicit none
      integer nCv, itcnt, iCv, style
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   aa(13*nCv), resfun(nCv)
      character*(*) name 
      real (kind=R8) ::
     *   corfun(nCv)
!csc This is a dummy routine to be fed to Tapenade isntead of the real b2uxus, 
!    such that only the necessary dependencies are created
      do iCv=1,nCv
        corfun(iCv) = resfun(iCv)*aa(iCv)
      end do
      return
      end subroutine b2uxus
