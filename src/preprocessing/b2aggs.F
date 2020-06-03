*-----------------------------------------------------------------------
*.specification

      subroutine b2aggs (ixlo, ixhi, iylo, iyhi, iloc)
      use b2mod_types
      use b2mod_b2cmgs
      implicit none
*   ..input arguments (unchanged on exit)
      integer ixlo, ixhi, iylo, iyhi, iloc
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGGS marks a singular region in the geometry.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine enters into the common block /b2cmgs/ parameters
*     describing one elementary singular region in the B2 geometry.
*     Such a singular region may be a cut associated with an x-point,
*     or it may be a larger domain on which the calculation is
*     special.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ixlo, ixhi, iylo, iyhi - integer, input.
*     ixlo, ixhi, iylo, and iyhi specify the bounds of one singular
*     region in the B2 geometry. Mesh locations (ix,iy) that satisfy
*     ixlo.le.ix.le.ixhi and iylo.le.iy.le.iyhi make up the region.
*
*     iloc - integer, input.
*     iloc specifies the class of quantities for which the given domain
*     is special, as follows:
*       iloc = 0: cell center quantities.
*       iloc = 1: cell x-face quantities.
*       iloc = 2: cell y-face quantities.
*       iloc = 3: cell corner quantities.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input/output.
*     /b2cmgs/ specifies the singular regions in the geometry.
*      nlrmx - integer constant.
*     nlrmx is an upper bound for nlreg, and is used in dimension
*     specifications.
*      nlreg - integer.
*     nlreg specifies the number of singular regions. It will hold that
*     0.le.nlreg.le.nlrmx.
*      nlxlo, nlxhi, nlylo, nlyhi - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.nlreg, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i) are
*     the bounds of the i-th singular region. Mesh locations (ix,iy)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (-1:nx,-1:ny).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (-1:nx,-1:ny).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
*
*
*  8. further comments
*
*     It is intended that /b2cmgs/ should either be read in and never
*     modified, or it should be initialised with nlreg.eq.0 and only be
*     modified through calls to b2aggs. An internal 'saved' variable
*     (nlrg0) in this routine is used to verify that the number of
*     singular regions is equal to the number of calls made to b2aggs.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nlrg0
      save nlrg0
*   ..procedures
      external subini, subend, xertst
*   ..initialisation
      data nlrg0/0/

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2aggs')
*   ..test iloc, nlreg
      call xertst (0.le.iloc.and.iloc.lt.4, 'faulty argument iloc')
      call xertst (nlreg.lt.nlrmx,
     &  'faulty call, too many null regions')
      call xertst (nlrg0.eq.nlreg, 'faulty call, nlrg0.ne.nlreg')
*   ..mark null region
      nlxlo(nlreg) = ixlo
      nlxhi(nlreg) = ixhi
      nlylo(nlreg) = iylo
      nlyhi(nlreg) = iyhi
      nlloc(nlreg) = iloc
      nlreg = nlreg+1
      nlrg0 = nlrg0+1

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2aggs

      end subroutine b2aggs
