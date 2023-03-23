*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar (nCv, ns, nregionv, BoRiS, switch, geo, mpg, 
     &  pl, dv, diag)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      real (kind=R8) :: BoRiS
*   ..output arguments (unspecified on entry)
      type (B2Diagnostic), intent (inout) :: diag
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MXAR computes scaled (L1-)norms of the residuals.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes scaled L1-norms of the residuals of the
*     equations for continuity, parallel momentum balance, electron heat
*     balance, atom heat balance, and charge balance. The scaling is
*     such that the resulting quantity has a dimension (1/time).
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, is
      real (kind=R8) ::
     *  navi, movi, mtvi, nevi, prvi, prvn, t0, t1, t2, t3, t4
*   ..procedures
      intrinsic sqrt, abs
      external xertst
*   ..initialisation

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
          navi = navi+geo%cvVol(iCv)*pl%na(iCv,is)
          movi = movi+geo%cvVol(iCv)*
     &     pl%na(iCv,is)*sqrt(pl%ti(iCv)*am(is)*mp)
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
       diag%aresco(is) = t0/navi
       diag%aresmo(is) = t1/movi
      enddo
*   ..compute aresmt
      t0 = 0.0_R8
      do iCv = 1, mpg%nCi
       t0 = t0+abs(dv%resmt(iCv))
      enddo
      diag%aresmt = t0/mtvi
*   ..compute nevi, prvi
      nevi = 0.0_R8
      prvi = 0.0_R8
      prvn = 0.0_R8
      do iCv = 1, mpg%nCi
!WG_TODO        if (resignore(ix,iy,1).eq.0) cycle
        nevi = nevi+geo%cvVol(iCv)*dv%ne(iCv)
        prvi = prvi+geo%cvVol(iCv)*dv%ne(iCv)*pl%te(iCv)
        if (switch%tn_style.eq.0) then
         do is = 0, ns-1
          prvi = prvi+geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv)
         enddo
         prvn = prvi !dummy
        elseif (switch%tn_style.eq.1) then
         do is = 0, ns-1
          if (.not.is_neutral(is)) then
           prvi = prvi+geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv)
          endif
         enddo
         prvn = prvi !dummy
        else
         do is = 0, ns-1
          if (.not.is_neutral(is).or.zn(is).ne.1) then
           prvi = prvi+geo%cvVol(iCv)*pl%na(iCv,is)*pl%ti(iCv)
          else
           prvn = prvn+geo%cvVol(iCv)*pl%na(iCv,is)*pl%tn(iCv)
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
      diag%areshe = t0/((1.5_R8+BoRiS)*prvi)
      diag%areshi = t1/((1.5_R8+BoRiS)*prvi)
      diag%arespo = t2/(qe*nevi)
      diag%areshn = t3/((1.5_R8+BoRiS)*prvn)
      diag%areskt = t4/(prvi)

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
     &    pl%na(iCv,is)*sqrt(pl%ti(iCv)*am(is)*mp)
       enddo
*   ..compute rescoreg, resmoreg
       do iCv = 1, mpg%nCi
!WG_TODO: resignore!
        diag%rescoreg(mpg%cvReg(iCv),is)=
     .     diag%rescoreg(mpg%cvReg(iCv),is)+abs(dv%resco(iCv,is))
        diag%resmoreg(mpg%cvReg(iCv),is)=
     .     diag%resmoreg(mpg%cvReg(iCv),is)+abs(dv%resmo(iCv,is))
       enddo
       diag%rescoreg(:,is) = diag%rescoreg(:,is)/navi
       diag%resmoreg(:,is) = diag%resmoreg(:,is)/movi
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
      diag%reshereg = diag%reshereg/((1.5_R8+BoRiS)*prvi)
      diag%reshireg = diag%reshireg/((1.5_R8+BoRiS)*prvi)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar

      end subroutine b2mxar
