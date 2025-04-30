c Calculate face distances and total distance
      subroutine calc_face_dist(geo, mpg, nCv, segment, list_cv, nncv,
     &   list_ci, nnci, list_fc, nnfc, segmdist, totaldist, ifail)
*   ..other modules called 
      use b2mod_types
      use b2mod_constants
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_b2cmpa
      use b2mod_subsys
      use b2mod_switches
      use b2mod_math
      implicit none
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind = R8), intent(in) :: segment(2,2)
      integer, intent(in):: nCv
c     Output variables
      real (kind=R8), intent(out) :: segmdist(mpg%nCi), totaldist
      integer, intent(out) :: list_ci(mpg%nCi), nnci, list_fc(mpg%nFc),
     &                        nnfc, list_cv(nCv), nncv
c     Internal variables
      real (kind=R8) :: XI(mpg%nFc), YI(mpg%nFc)
      integer :: face1, face2, iFc, iCv
      logical :: ifail

      call subini ('calc_face_dist')

      call ifill (nCv, 0, list_cv, 1)
      call ifill (mpg%nCi, 0, list_ci, 1)
      call ifill (mpg%nFc, 0, list_fc, 1)
      nncv = 0
      nnci = 0
      nnfc = 0
      call sfill (mpg%nCi, 0.0_R8, segmdist, 1)
      totaldist = 0.0_R8

      call cv_intersections (geo, mpg, segment, list_cv, list_ci,
     &   list_fc, nncv, nnci, nnfc, nCv, XI, YI, 1, ifail)

      do iFc = 1, nnfc-1
        iCv = list_ci(iFc)
        face1 = list_fc(iFc)
        face2 = list_fc(iFc+1)
        segmdist(iCv) = sqrt((XI(face1)-XI(face2))**2+
     &                       (YI(face1)-YI(face2))**2)
        totaldist = totaldist + segmdist(iCv)
      enddo

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end calc_face_dist

      end subroutine calc_face_dist
