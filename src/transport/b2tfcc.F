*-----------------------------------------------------------------------
*.specification

      subroutine b2tfcc (nCv, nFc, nVx, ns, switch, geo, mpg, pl, 
     &                   dv, co, rt)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tfcc
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type(B2Plasma), intent (in) :: pl 
      type(B2Coeff), intent (in) :: co
      type(B2Rates), intent (in) :: rt
*   ..input/output arguments
      type(B2Derivatives), intent (inout) :: dv
*   ..common blocks
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, iFc
      real (kind=R8) ::
     *   nw(nCv), nwf(nFc),
     *   dnw(nFc,0:1), tef(nFc), tif(nFc), tnf(nFc),
     *   wght(nFc,2), wrkVx(nVx),
     *   flo(nFc,0:1), con(nFc,0:1), flo0(nFc,0:1), con0(nFc,0:1)
*   ..procedures
      external subini, subend, xertst
      external b2xvsg, b2xvff, diff, intface
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfcc')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2tfcc.lt.3) then
*    ..test sign of cdpa
       do is = 0, ns-1
        call b2xvsg (nFc, co%cdpa(:,0,is)*geo%fcQalf(:,0),1,'cdpa0',
     .    '.ge.')
        call b2xvsg (nFc, co%cdpa(:,1,is)*geo%fcQalf(:,1),1,'cdpa1',
     .    '.ge.')
       enddo
      endif

*   ..interpolate te and ti to cell faces
      wght = 1.0_R8
      call intface(nCv,nFc,mpg%fcCv,wght,pl%te,tef)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%ti,tif)
      call intface(nCv,nFc,mpg%fcCv,wght,pl%tn,tnf)

* ..correct floe, cone

      do is = 0, ns-1

*    ..pre-compute radial differences
       nw = pl%na(:,is)*(rt%rza(:,is) + 1.0_R8)*rt%rza(:,is)
       call diff(nCv, nFc, nVx, 0, geo, mpg, nw, wrkVx, dnw)

*    ..interpolate nwe to cell faces
       call intface(nCv,nFc,mpg%fcCv,wght,nw,nwf)

*    ..compute poloidal/radial coefficients for pressure diffusion term
       flo(:,0) = -switch%xfac*2.5_R8*co%cdpa(:,0,is)*dnw(:,0)*tef
       flo(:,1) = -switch%xfac*2.5_R8*co%cdpa(:,1,is)*dnw(:,1)*tef
       con(:,0) =  switch%xfac*2.5_R8*co%cdpa(:,0,is)*nwf*tef
       con(:,1) =  switch%xfac*2.5_R8*co%cdpa(:,1,is)*nwf*tef

*    ..compute coefficients for stencil
       call calccoef(nCv,nFc,nVx,0,geo,mpg,flo,con,flo0,con0)

*    ..add to previous contributions
       dv%floe        = dv%floe        + flo0
       dv%cone(:,:,0) = dv%cone(:,:,0) + con0

      enddo !! end loop over species

*  ..add contributions from Rhie-Chow term

*   ..pre-compute radial differences
      call diff(nCv, nFc, nVx, 0, geo, mpg, dv%ne, wrkVx, dnw)

*   ..interpolate nwe and nwi to cell faces
      call intface(nCv,nFc,mpg%fcCv,wght,dv%ne,nwf)

      do iFc = 1, nFc
        con(:,0) = switch%b2npht_pcm1*1.25_R8*geo%fcPbs(iFc)*
     &   (dv%pccm(mpg%fcCv(iFc,1),1)*
     &    dv%ne(mpg%fcCv(iFc,1))*pl%te(mpg%fcCv(iFc,1)) +
     &    dv%pccm(mpg%fcCv(iFc,2),1)*
     &    dv%ne(mpg%fcCv(iFc,2))*pl%te(mpg%fcCv(iFc,2)) )/
     &   (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
      enddo

*   ..compute poloidal/radial coefficients for pressure diffusion term
      flo(:,0) = -con(:,0)*dnw(:,0)
      con(:,0) =  con(:,0)*nwf

      flo(:,1) = 0.0_R8
      con(:,1) = 0.0_R8

*   ..compute coefficients for stencil
      call calccoef(nCv,nFc,nVx,0,geo,mpg,flo,con,flo0,con0)

*   ..add to previous contributions
      dv%floe        = dv%floe        + flo0
      dv%cone(:,:,0) = dv%cone(:,:,0) + con0

* ..correct floi, coni, flon, conn
      do is = 0, ns-1

*    ..pre-compute radial differences
       nw = pl%na(:,is)*(rt%rza(:,is) + 1.0_R8)
       call diff(nCv, nFc, nVx, 0, geo, mpg, nw, wrkVx, dnw)

*    ..interpolate nwe to cell faces
       call intface(nCv,nFc,mpg%fcCv,wght,nw,nwf)

*    ..compute poloidal/radial coefficients for pressure diffusion term
       if (.not.is_neutral(is)) then
         flo(:,0) = -switch%xfac*2.5_R8*co%cdpa(:,0,is)*dnw(:,0)*tif
         flo(:,1) = -switch%xfac*2.5_R8*co%cdpa(:,1,is)*dnw(:,1)*tif
         con(:,0) =  switch%xfac*2.5_R8*co%cdpa(:,0,is)*nwf*tif
         con(:,1) =  switch%xfac*2.5_R8*co%cdpa(:,1,is)*nwf*tif
       else
         flo(:,0) = -switch%xfac*2.5_R8*co%cdpa(:,0,is)*dnw(:,0)*tnf !! in case of default energy equation Tn = Ti
         flo(:,1) = -switch%xfac*2.5_R8*co%cdpa(:,1,is)*dnw(:,1)*tnf
         con(:,0) =  switch%xfac*2.5_R8*co%cdpa(:,0,is)*nwf*tnf
         con(:,1) =  switch%xfac*2.5_R8*co%cdpa(:,1,is)*nwf*tnf
       endif

*    ..compute coefficients for stencil
       call calccoef(nCv,nFc,nVx,0,geo,mpg,flo,con,flo0,con0)

*    ..add to previous contributions
c the logic statements could be simplified but for now this is not done to improve readability
       if (switch%tn_style.eq.0) then !! default combined ion-neutral energy equation, include both ions and neutrals in "ion" terms.
         dv%floi        = dv%floi        + flo0
         dv%coni(:,:,0) = dv%coni(:,:,0) + con0
       elseif (switch%tn_style.eq.1) then !! pure ion energy equation, only add ion contributions
         if (.not.is_neutral(is)) then
           dv%floi        = dv%floi        + flo0
           dv%coni(:,:,0) = dv%coni(:,:,0) + con0
         endif
       elseif (switch%tn_style.eq.2) then !! separate ion and neutral energy equation, split contributions
         if (.not.is_neutral(is)) then
           dv%floi        = dv%floi        + flo0
           dv%coni(:,:,0) = dv%coni(:,:,0) + con0
         else
           dv%flon        = dv%flon        + flo0
           dv%conn(:,:,0) = dv%conn(:,:,0) + con0
         endif
       endif

*  ..add contributions from Rhie-Chow term

*    ..pre-compute radial differences
       call diff(nCv, nFc, nVx, 0, geo, mpg, pl%na(1,is), wrkVx, dnw)

*    ..interpolate nwe to cell faces
       call intface(nCv,nFc,mpg%fcCv,wght,nw,nwf)

       if (.not.is_neutral(is)) then
        do iFc = 1, nFc
          con(:,0) = switch%b2npht_pcm1*1.25_R8*geo%fcPbs(iFc)*
     &     (dv%pccm(mpg%fcCv(iFc,1),1)*
     &      pl%na(mpg%fcCv(iFc,1),is)*pl%ti(mpg%fcCv(iFc,1))  +
     &      dv%pccm(mpg%fcCv(iFc,2),1)*
     &      pl%na(mpg%fcCv(iFc,2),0)*pl%ti(mpg%fcCv(iFc,2)) )/
     &     (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
        enddo
       else
        do iFc = 1, nFc
          con(:,0) = switch%b2npht_pcm1*1.25_R8*geo%fcPbs(iFc)*
     &     (dv%pccm(mpg%fcCv(iFc,1),1)*
     &      pl%na(mpg%fcCv(iFc,1),0)*pl%tn(mpg%fcCv(iFc,1))  +
     &      dv%pccm(mpg%fcCv(iFc,2),1)*
     &      pl%na(mpg%fcCv(iFc,2),0)*pl%tn(mpg%fcCv(iFc,2)) )/
     &     (geo%fcQgam(iFc,0)*(geo%fcHc(iFc,1) + geo%fcHc(iFc,2)))
        enddo
       endif

*    ..compute poloidal/radial coefficients for pressure diffusion term
       flo(:,0) = -con(:,0)*dnw(:,0)
       con(:,0) =  con(:,0)*nwf

       flo(:,1) = 0.0_R8
       con(:,1) = 0.0_R8

*    ..compute coefficients for stencil
       call calccoef(nCv,nFc,nVx,0,geo,mpg,flo,con,flo0,con0)

*    ..add to previous contributions
c the logic statements could be simplified but for now this is not done to improve readability
       if (switch%tn_style.eq.0) then !! default combined ion-neutral energy equation, include both ions and neutrals in "ion" terms.
         dv%floi        = dv%floi        + flo0
         dv%coni(:,:,0) = dv%coni(:,:,0) + con0
       elseif (switch%tn_style.eq.1) then !! pure ion energy equation, only add ion contributions
         if (.not.is_neutral(is)) then
           dv%floi        = dv%floi        + flo0
           dv%coni(:,:,0) = dv%coni(:,:,0) + con0
         endif
       elseif (switch%tn_style.eq.2) then !! separate ion and neutral energy equation, split contributions
         if (.not.is_neutral(is)) then
           dv%floi        = dv%floi        + flo0
           dv%coni(:,:,0) = dv%coni(:,:,0) + con0
         else
           dv%flon        = dv%flon        + flo0
           dv%conn(:,:,0) = dv%conn(:,:,0) + con0
         endif
       endif

      enddo

* ..return
      ncall_b2tfcc = ncall_b2tfcc+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfcc

      end
