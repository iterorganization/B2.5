*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar_diff (nCv, ns, nregionv, BoRiS, geo, mpg, pl, 
     &  dv, diag)
      use b2mod_types
      use b2mod_constants
!      use b2mod_rates
      use b2mod_b2cmpa_diff
      use b2us_geo_diff
      use b2us_map_diff
      use b2us_plasma_diff
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nregionv
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma_diff), intent (in) :: pl
      type (B2Derivatives_diff), intent (in) :: dv
      real (kind=R8) :: BoRiS
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic), intent (inout) :: diag
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
* This is a modified version of b2mxar used in the B2.5 code differentiated
* with TAPENADE. The main difference wrt b2mxar is that some differentiated
* quantities can be .le. 0.0, therefore sqrt gives nans and results cannot
* be plotted against log scale. Therefore all sums are done using absolute
* values of differentiated quantities.

*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, iCv, is
      real (kind=R8) ::
     *  navi, movi, mtvi, nevi, prvi, t0, t1, t2, t3
      save ncall
*   ..procedures
      intrinsic sqrt, abs
      real (kind=R8) ::
     *  b2sasum
      external subini, subend, xertst, b2sasum
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxar')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')

* ..compute scaled norms
*   ..compute aresco, aresmo, mtvi
      mtvi = 0.0_R8
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iCv = 1, mpg%nCi
          navi = navi+geo%cvVol(iCv)*pl%na(iCv,is)
          movi = movi+geo%cvVol(iCv)*
     &     pl%na(iCv,is)*sqrt(abs(pl%ti(iCv)*am(is)*mp))
!WG_TODO         if (resignore(ix,iy,1).ne.0) 
!WG_TODO     >    navi = navi+geo%cvVol(iCv)*pl%na(iCv,is)
!WG_TODO         if (resignore(ix,iy,2).ne.0)
!WG_TODO     >    movi = movi+geo%cvVol(iCv)*
!WG_TODO     &     pl%na(iCv,is)*sqrt(pl%ti(iCv)*am(is)*mp)
       enddo
       mtvi = mtvi+movi
       t0 = 0.0_R8
       t1 = 0.0_R8
       do iCv = 1, mpg%nCi
        t0 = t0+abs(dv%resco(iCv,is))
        t1 = t1+abs(dv%resmo(iCv,is))
       enddo
       diag%aresco(is) = abs(t0/navi)
       diag%aresmo(is) = abs(t1/movi)
      enddo
*   ..compute aresmt
      t0 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resmt(iCv))
      enddo
      diag%aresmt = abs(t0/mtvi)
*   ..compute nevi, prvi
      nevi = 0.0_R8
      prvi = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(ix,iy,1).eq.0) cycle
        nevi = nevi+geo%cvVol(iCv)*dv%ne(iCv)
        prvi = prvi+geo%cvVol(iCv)*dv%ne(iCv)*pl%te(iCv)
        do is = 0, ns-1
         prvi = prvi+geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv)
        enddo
      enddo
*   ..compute areshe, areshi, arespo
      t0 = 0.0_R8
      t1 = 0.0_R8
      t2 = 0.0_R8
      t3 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%reshe(iCv))
       t1 = t1+abs(dv%reshi(iCv))
       t2 = t2+abs(dv%respo(iCv))
       t3 = t3+abs(dv%reshn(iCv))
      enddo
      diag%areshe = abs(t0/((1.5_R8+BoRiS)*prvi))
      diag%areshi = abs(t1/((1.5_R8+BoRiS)*prvi))
      diag%arespo = abs(t2/(qe*nevi))
      diag%areshn = abs(t3/((1.5_R8+BoRiS)*prvi))

      !! AFN_TODO split prvi into prvi and prvn
      !! AFN_TODO add to plot routines

      write(*,*) 'areshn =', diag%areshn

* ..compute regional scaled sums
      diag%rescoreg = 0.0_R8
      diag%resmoreg = 0.0_R8
*   ..compute navi, movi
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        navi = navi+geo%cvVol(iCv)*pl%na(iCv,is)
        movi = movi+geo%cvVol(iCv)*
     &    pl%na(iCv,is)*sqrt(abs(pl%ti(iCv)*am(is)*mp))
       enddo
*   ..compute rescoreg, resmoreg
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%rescoreg(mpg%cvReg(iCv),is)=
     .     diag%rescoreg(mpg%cvReg(iCv),is)+abs(dv%resco(iCv,is))
        diag%resmoreg(mpg%cvReg(iCv),is)=
     .     diag%resmoreg(mpg%cvReg(iCv),is)+abs(dv%resmo(iCv,is))
       enddo
       diag%rescoreg(:,is) = abs(diag%rescoreg(:,is)/navi)
       diag%resmoreg(:,is) = abs(diag%resmoreg(:,is)/movi)
      enddo
*   ..compute prvi
      prvi = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(ix,iy,1).eq.0) cycle
        prvi = prvi+geo%cvVol(iCv)*dv%ne(iCv)*pl%te(iCv)
        do is = 0, ns-1
         if (.not.is_neutral(is)) 
     >    prvi = prvi+geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv)
        enddo
      enddo
*   ..compute reshereg, reshireg
      diag%reshereg = 0.0_R8
      diag%reshireg = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%reshereg(mpg%cvReg(iCv)) = 
     .    diag%reshereg(mpg%cvReg(iCv))+abs(dv%reshe(iCv))
        diag%reshireg(mpg%cvReg(iCv)) = 
     .    diag%reshireg(mpg%cvReg(iCv))+abs(dv%reshi(iCv))
      enddo
      diag%reshereg = abs(diag%reshereg/((1.5_R8+BoRiS)*prvi))
      diag%reshireg = abs(diag%reshireg/((1.5_R8+BoRiS)*prvi))

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar

      end subroutine b2mxar_diff
