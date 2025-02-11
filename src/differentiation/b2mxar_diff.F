*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar_diff (nCv, ns, nregionv, BoRiS, switch, geo,
     &  mpg, pl, dv, diag)
      use b2mod_types
      use b2mod_constants
!      use b2mod_rates
      use b2mod_b2cmpa_diff
      use b2mod_switches_diff
      use b2us_geo_diff
      use b2us_map_diff
      use b2us_plasma_diff
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma_diff), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
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
     *  navi, movi, mtvi, nevi, prvi, prvn, t0, t1, t2, t3, t4
      save ncall
*   ..procedures
      intrinsic sqrt, abs
      external xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxar')
*   ..set internal parameters on first call
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')

* ..compute scaled norms
*   ..compute aresco, aresmo, mtvi
      mtvi = 0.0_R8
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iCv = 1, mpg%nCi
          navi = navi+abs(geo%cvVol(iCv)*pl%na(iCv,is))
          movi = movi+abs(geo%cvVol(iCv)*
     &     pl%na(iCv,is)*sqrt(abs(pl%ti(iCv)*am(is)*mp)))
!WG_TODO         if (resignore(iCv,1).ne.0)
!WG_TODO     >    navi = navi+geo%cvVol(iCv)*pl%na(iCv,is)
!WG_TODO         if (resignore(iCv,2).ne.0)
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
       if (navi .gt. 0.0_R8) diag%aresco(is) = abs(t0/navi)
       if (movi .gt. 0.0_R8) diag%aresmo(is) = abs(t1/movi)
      enddo
*   ..compute aresmt
      t0 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resmt(iCv))
      enddo
      if (mtvi .gt. 0.0_R8) diag%aresmt = abs(t0/mtvi)
*   ..compute nevi, prvi
      nevi = 0.0_R8
      prvi = 0.0_R8
      prvn = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(iCv,1).eq.0) cycle
        nevi = nevi+abs(geo%cvVol(iCv)*dv%ne(iCv))
        prvi = prvi+abs(geo%cvVol(iCv)*dv%ne(iCv)*pl%te(iCv))
        if (switch%tn_style.eq.0) then
         do is = 0, ns-1
          prvi = prvi+abs(geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv))
         enddo
         prvn = prvi
        else
         do is = 0, ns-1
          if (.not.is_neutral(is).or.nint(zn(is)).ne.1) then
           prvi = prvi+abs(geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv))
          else
           prvn = prvn+abs(geo%cvVol(iCv)*pl%na(iCv,is)*pl%tn(iCv))
          endif
         enddo
        endif
      enddo
*   ..compute areshe, areshi, arespo
      t0 = 0.0_R8
      t1 = 0.0_R8
      t2 = 0.0_R8
      t3 = 0.0_R8
      t4 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%reshe(iCv))
       t1 = t1+abs(dv%reshi(iCv))
       t2 = t2+abs(dv%respo(iCv))
       t3 = t3+abs(dv%reshn(iCv))
       t4 = t4+abs(dv%reskt(iCv))
      enddo
      if (prvi .gt. 0.0_R8) diag%areshe = abs(t0/((1.5_R8+BoRiS)*prvi))
      if (prvi .gt. 0.0_R8) diag%areshi = abs(t1/((1.5_R8+BoRiS)*prvi))
      if (nevi .gt. 0.0_R8) diag%arespo = abs(t2/(qe*nevi))
      if (prvn .gt. 0.0_R8) diag%areshn = abs(t3/((1.5_R8+BoRiS)*prvn))
      if (prvi .gt. 0.0_R8) diag%areskt = abs(t4/(prvi))

* ..compute regional scaled sums
      diag%rescoreg = 0.0_R8
      diag%resmoreg = 0.0_R8
*   ..compute navi, movi
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        navi = navi+abs(geo%cvVol(iCv)*pl%na(iCv,is))
        movi = movi+abs(geo%cvVol(iCv)*
     &    pl%na(iCv,is)*sqrt(abs(pl%ti(iCv)*am(is)*mp)))
       enddo
*   ..compute rescoreg, resmoreg
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%rescoreg(mpg%cvReg(iCv),is)=
     .     diag%rescoreg(mpg%cvReg(iCv),is)+abs(dv%resco(iCv,is))
        diag%resmoreg(mpg%cvReg(iCv),is)=
     .     diag%resmoreg(mpg%cvReg(iCv),is)+abs(dv%resmo(iCv,is))
       enddo
       if (navi .gt. 0.0_R8)
     .  diag%rescoreg(:,is) = abs(diag%rescoreg(:,is)/navi)
       if (movi .gt. 0.0_R8)
     .  diag%resmoreg(:,is) = abs(diag%resmoreg(:,is)/movi)
      enddo
*   ..compute prvi
      prvi = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(iCv,1).eq.0) cycle
        prvi = prvi+abs(geo%cvVol(iCv)*dv%ne(iCv)*pl%te(iCv))
        do is = 0, ns-1
         if (.not.is_neutral(is))
     >    prvi = prvi+abs(geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv))
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
      if (prvi .gt. 0.0_R8)
     . diag%reshereg = abs(diag%reshereg/((1.5_R8+BoRiS)*prvi))
      if (prvi .gt. 0.0_R8)
     . diag%reshireg = abs(diag%reshireg/((1.5_R8+BoRiS)*prvi))

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar

      end subroutine b2mxar_diff
