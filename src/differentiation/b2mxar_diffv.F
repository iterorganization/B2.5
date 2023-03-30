*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar_diffv (nCv, ns, nregionv, BoRiS, switch, geo, 
     &  mpg, pl, dv, diag)
      use b2mod_types
      use b2mod_constants
!      use b2mod_rates
      use b2mod_b2cmpa_diffv
      use b2mod_switches_diffv
      use b2us_geo_diffv
      use b2us_map_diffv
      use b2us_plasma_diffv
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma_diffv), intent (in) :: pl
      type (B2Derivatives_diffv), intent (in) :: dv
      real (kind=R8) :: BoRiS
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic_diffv), intent (inout) :: diag
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
      external subini, subend, xertst
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxar')
*   ..set internal parameters on first call
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
          navi = navi+abs(geo%cvVol(iCv)*pl%na(1,iCv,is))
          movi = movi+abs(geo%cvVol(iCv)*
     &     pl%na(1,iCv,is)*sqrt(abs(pl%ti(1,iCv)*am(is)*mp)))
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
        t0 = t0+abs(dv%resco(1,iCv,is))
        t1 = t1+abs(dv%resmo(1,iCv,is))
       enddo
       diag%aresco(1,is) = abs(t0/navi)
       diag%aresmo(1,is) = abs(t1/movi)
      enddo
*   ..compute aresmt
      t0 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resmt(1,iCv))
      enddo
      diag%aresmt(1) = abs(t0/mtvi)
*   ..compute nevi, prvi
      nevi = 0.0_R8
      prvi = 0.0_R8
      prvn = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(ix,iy,1).eq.0) cycle
        nevi = nevi+abs(geo%cvVol(iCv)*dv%ne(1,iCv))
        prvi = prvi+abs(geo%cvVol(iCv)*dv%ne(1,iCv)*pl%te(1,iCv))
        if (switch%tn_style.eq.0) then
         do is = 0, ns-1
          prvi = prvi+abs(geo%cvVol(iCv)*pl%na(1,iCv,is)*pl%ti(1,iCv))
         enddo
         prvn = prvi
        else
         do is = 0, ns-1
          if (.not.is_neutral(is).or.(switch%tn_style.eq.2.and.
     &       (nint(zn(is)).ne.1.or.is.ne.ismain0))) then
           prvi = prvi+abs(geo%cvVol(iCv)*pl%na(1,iCv,is)*pl%ti(1,iCv))
          else
           prvn = prvn+abs(geo%cvVol(iCv)*pl%na(1,iCv,is)*pl%tn(1,iCv))
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
       t0 = t0+abs(dv%reshe(1,iCv))
       t1 = t1+abs(dv%reshi(1,iCv))
       t2 = t2+abs(dv%respo(1,iCv))
       t3 = t3+abs(dv%reshn(1,iCv))
       t4 = t4+abs(dv%reskt(1,iCv))
      enddo
      diag%areshe = abs(t0/((1.5_R8+BoRiS)*prvi))
      diag%areshi = abs(t1/((1.5_R8+BoRiS)*prvi))
      diag%arespo = abs(t2/(qe*nevi))
      diag%areshn = abs(t3/((1.5_R8+BoRiS)*prvn))
      diag%areskt = abs(t4/(prvi))

* ..compute regional scaled sums
      diag%rescoreg = 0.0_R8
      diag%resmoreg = 0.0_R8
*   ..compute navi, movi
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        navi = navi+abs(geo%cvVol(iCv)*pl%na(1,iCv,is))
        movi = movi+abs(geo%cvVol(iCv)*
     &    pl%na(1,iCv,is)*sqrt(abs(pl%ti(1,iCv)*am(is)*mp)))
       enddo
*   ..compute rescoreg, resmoreg
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%rescoreg(1,mpg%cvReg(iCv),is)=
     .     diag%rescoreg(1,mpg%cvReg(iCv),is)+abs(dv%resco(1,iCv,is))
        diag%resmoreg(1,mpg%cvReg(iCv),is)=
     .     diag%resmoreg(1,mpg%cvReg(iCv),is)+abs(dv%resmo(1,iCv,is))
       enddo
       diag%rescoreg(1,:,is) = abs(diag%rescoreg(1,:,is)/navi)
       diag%resmoreg(1,:,is) = abs(diag%resmoreg(1,:,is)/movi)
      enddo
*   ..compute prvi
      prvi = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(ix,iy,1).eq.0) cycle
        prvi = prvi+abs(geo%cvVol(iCv)*dv%ne(1,iCv)*pl%te(1,iCv))
        do is = 0, ns-1
         if (.not.is_neutral(is)) 
     >    prvi = prvi+abs(geo%cvVol(iCv)*pl%na(1,iCv,is)*pl%ti(1,iCv))
        enddo
      enddo
*   ..compute reshereg, reshireg
      diag%reshereg = 0.0_R8
      diag%reshireg = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%reshereg(1,mpg%cvReg(iCv)) = 
     .    diag%reshereg(1,mpg%cvReg(iCv))+abs(dv%reshe(1,iCv))
        diag%reshireg(1,mpg%cvReg(iCv)) = 
     .    diag%reshireg(1,mpg%cvReg(iCv))+abs(dv%reshi(1,iCv))
      enddo
      diag%reshereg(1,:) = abs(diag%reshereg(1,:)/((1.5_R8+BoRiS)*prvi))
      diag%reshireg(1,:) = abs(diag%reshireg(1,:)/((1.5_R8+BoRiS)*prvi))

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar

      end subroutine b2mxar_diffv
