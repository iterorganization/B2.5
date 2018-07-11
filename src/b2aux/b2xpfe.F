      subroutine b2xpfe (nCv, nFc, ns, ns_ext, geo, mpg, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
      integer nCv, nFc, ns, ns_ext
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  qe, rza(nCv,0:ns-1), fna(nFc,0:1,0:ns-1),
     *  za_ext(nCv,0:max(ns_ext,1)-1), 
     *  fa_ext(nFc,0:1,0:max(ns_ext,1)-1),
     &  fch(nFc,0:1), fne(nFc,0:1)
*     ------------------------------------------------------------------
*     B2XPFE computes the electron particle flux, fne:
*       fne(,) = (sum is :: rza(,is)*fna(,,is))-(1/qe)*fch(,) .
*     A contribution from external species can also be included.
*     ------------------------------------------------------------------
      integer iFc, is
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfe')
      call sfill (nFc*2, 0.0_R8, fne, 1)
      do is = 0, ns-1
        do iFc = 1, nFc
          fne(iFc,0:1) = fne(iFc,0:1) + fna(iFc,0:1,is)*
     .      (rza(mpg%fcCv(iFc,1),is)*geo%cvVol(mpg%fcCv(iFc,2))+
     .       rza(mpg%fcCv(iFc,2),is)*geo%cvVol(mpg%fcCv(iFc,1)))/
     .      (geo%cvVol(mpg%fcCv(iFc,1))+geo%cvVol(mpg%fcCv(iFc,2)))
        enddo
      enddo
      do is = 0, ns_ext-1
        do iFc = 1, nFc
          fne(iFc,0:1) = fne(iFc,0:1) + fa_ext(iFc,0:1,is)*
     .      (za_ext(mpg%fcCv(iFc,1),is)*geo%cvVol(mpg%fcCv(iFc,2))+
     .       za_ext(mpg%fcCv(iFc,2),is)*geo%cvVol(mpg%fcCv(iFc,1)))/
     .      (geo%cvVol(mpg%fcCv(iFc,1))+geo%cvVol(mpg%fcCv(iFc,2)))
        enddo
      enddo
      call b2saxpy (nFc*2, -1.0_R8/qe, fch, 1, fne, 1)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfe
