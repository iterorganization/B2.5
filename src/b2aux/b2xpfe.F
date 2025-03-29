      subroutine b2xpfe (nCv, nFc, ns, ns_ext, geo, mpg, qe, 
     &  rza, za_ext, fna, fa_ext, fch, fne)
      use b2mod_types
      use b2us_geo
      use b2us_map
#ifdef _OPENMP
      use b2mod_openmp
#endif
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
      integer is
      real (kind=R8) :: wrkf(nFc)
      external sfill, b2saxpy
*     ------------------------------------------------------------------
*$$$  call subini ('b2xpfe')
      call sfill (nFc*2, 0.0_R8, fne, 1)

#if (defined(_OPENMP) && !defined(NO_OPENMP_B2XPFE)) 
! fixed form fortran deserves to die imho
!$omp parallel default(none) private(is, wrkf) 
!$omp& shared(nCv, nFc, mpg, geo, rza, fna, ns, ns_ext, za_ext, 
!$omp& fa_ext)
!$omp& reduction(+:fne)
!$omp do
#endif 
      do is = 0, ns-1
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,rza(:,is),wrkf)
        fne(:,0) = fne(:,0) + fna(:,0,is)*wrkf
        fne(:,1) = fne(:,1) + fna(:,1,is)*wrkf
      enddo
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2XPFE)) 
!$omp end do
!$omp do
#endif 
      do is = 0, ns_ext-1
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,za_ext(:,is),wrkf)
        fne(:,0) = fne(:,0) + fa_ext(:,0,is)*wrkf
        fne(:,1) = fne(:,1) + fa_ext(:,1,is)*wrkf
      enddo
#if (defined(_OPENMP) && !defined(NO_OPENMP_B2XPFE))
!$omp end do
!$omp end parallel 
#endif 
      call b2saxpy (nFc*2, -1.0_R8/qe, fch, 1, fne, 1)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2xpfe
