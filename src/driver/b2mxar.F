*-----------------------------------------------------------------------
*.specification

      subroutine b2mxar (nx, ny, ns, nregionv, BoRiS, 
     &  vol, ne, na, te, ti, 
     &  resco, resmo, resmt, reshe, reshi, respo, 
     &  rescoreg, resmoreg, reshereg, reshireg,
     &  aresco, aresmo, aresmt, areshe, areshi, arespo)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nregionv
      real (kind=R8) :: BoRiS,
     *  vol(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1),
     &  ne(-1:nx,-1:ny), te(-1:nx,-1:ny), ti(-1:nx,-1:ny),
     &  resco(-1:nx,-1:ny,0:ns-1), resmo(-1:nx,-1:ny,0:ns-1),
     &  resmt(-1:nx,-1:ny), reshe(-1:nx,-1:ny), reshi(-1:nx,-1:ny),
     &  respo(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  rescoreg(0:nregionv,0:ns-1), resmoreg(0:nregionv,0:ns-1), 
     *  reshereg(0:nregionv), reshireg(0:nregionv)
      real (kind=R8) ::
     *  aresco(0:ns-1), aresmo(0:ns-1), aresmt, areshe, areshi,
     &  arespo
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
      integer ncall, ix, iy, is
      real (kind=R8) ::
     *  navi, movi, mtvi, nevi, prvi, t0, t1, t2
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
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')

* ..compute scaled norms
*   ..compute aresco, aresmo, navi, movi, mtvi
      mtvi = 0.0_R8
      do is = 0, ns-1
       t0 = 0.0_R8
       t1 = 0.0_R8
       navi = 0.0_R8
       movi = 0.0_R8
       do iy = -1, ny
        do ix =-1, nx
         if (resignore(ix,iy,1).ne.0) then
           t0 = t0 + abs(resco(ix,iy,is))
           navi = navi+vol(ix,iy)*na(ix,iy,is)
         end if
         if (resignore(ix,iy,2).ne.0) then
           t1 = t1 + abs(resmo(ix,iy,is))
           movi = movi+vol(ix,iy)*na(ix,iy,is)*sqrt(ti(ix,iy)*am(is)*mp)
         end if
        enddo
       enddo
       mtvi = mtvi+movi
       aresco(is) = t0/navi
       aresmo(is) = t1/movi
      enddo
*   ..compute aresmt
      t0 = 0.0_R8
      do iy = -1, ny
       do ix =-1, nx
        if (resignore(ix,iy,2).ne.0)t0 = t0 + abs(resmt(ix,iy))
       enddo
      enddo
      aresmt = t0/mtvi
*   ..compute areshe, areshi, arespo, nevi, prvi
      t0 = 0.0_R8
      t1 = 0.0_R8
      t2 = 0.0_R8
      nevi = 0.0_R8
      prvi = 0.0_R8
      do iy = -1, ny
       do ix = -1, nx
        if (resignore(ix,iy,1).eq.0) cycle
        t0 = t0 + abs(reshe(ix,iy))
        t1 = t1 + abs(reshi(ix,iy))
        t2 = t2 + abs(respo(ix,iy))
        nevi = nevi+vol(ix,iy)*ne(ix,iy)
        prvi = prvi+vol(ix,iy)*ne(ix,iy)*te(ix,iy)
        do is = 0, ns-1
         prvi = prvi+vol(ix,iy)*na(ix,iy,is)*ti(ix,iy)
        enddo
       enddo
      enddo
      areshe = t0/((1.5_R8+BoRiS)*prvi)
      areshi = t1/((1.5_R8+BoRiS)*prvi)
      arespo = t2/(qe*nevi)

* ..compute regional scaled sums
      rescoreg = 0.0_R8
      resmoreg = 0.0_R8
*   ..compute navi, movi
      do is = 0, ns-1
       navi = 0.0_R8
       movi = 0.0_R8
       do iy = 0, ny-1
        do ix = 0, nx-1
         if (resignore(ix,iy,1).ne.0)
     >    navi = navi+vol(ix,iy)*na(ix,iy,is)
         if (resignore(ix,iy,2).ne.0)
     >    movi = movi+vol(ix,iy)*
     &     na(ix,iy,is)*sqrt(ti(ix,iy)*am(is)*mp)
        enddo
       enddo
*   ..compute rescoreg, resmoreg
       do iy = 0, ny-1
        do ix = 0, nx-1
         if (resignore(ix,iy,1).ne.0)
     >    rescoreg(region(ix,iy,0),is)=
     .     rescoreg(region(ix,iy,0),is)+abs(resco(ix,iy,is))
         if (resignore(ix,iy,2).ne.0)
     >    resmoreg(region(ix,iy,0),is)=
     .     resmoreg(region(ix,iy,0),is)+abs(resmo(ix,iy,is))
        enddo
       enddo
       rescoreg(:,is) = rescoreg(:,is)/navi
       resmoreg(:,is) = resmoreg(:,is)/movi
      enddo
*   ..compute prvi
      prvi = 0.0_R8
      do iy = -1, ny
       do ix = -1, nx
        if (resignore(ix,iy,1).eq.0) cycle
        prvi = prvi+vol(ix,iy)*ne(ix,iy)*te(ix,iy)
        do is = 0, ns-1
         if (.not.is_neutral(is))
     >    prvi = prvi+vol(ix,iy)*na(ix,iy,is)*ti(ix,iy)
        enddo
       enddo
      enddo
*   ..compute reshereg, reshireg
      reshereg = 0.0_R8
      reshireg = 0.0_R8
      do iy = 0, ny-1
       do ix = 0, nx-1
        if (resignore(ix,iy,1).ne.0)
     >   reshereg(region(ix,iy,0)) =
     .    reshereg(region(ix,iy,0))+abs(reshe(ix,iy))
        if (resignore(ix,iy,1).ne.0)
     >   reshireg(region(ix,iy,0)) =
     .    reshireg(region(ix,iy,0))+abs(reshi(ix,iy))
       enddo
      enddo
      reshereg = reshereg/((1.5_R8+BoRiS)*prvi)
      reshireg = reshireg/((1.5_R8+BoRiS)*prvi)

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxar

      end subroutine b2mxar
