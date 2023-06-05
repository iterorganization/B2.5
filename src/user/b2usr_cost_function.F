      subroutine b2usr_cost_function(nCv, nFc, nVx, ns, geo, mpg, st,
     & st_ext, BoRiS, J)
      use b2mod_types
      use b2mod_constants
      use b2mod_user_namelist
     & , only : omp, nomp, icsepomp
      use b2mod_b2cmpa
     & , only : am
      use b2us_map
      use b2us_plasma
      use b2us_geo
      use b2mod_par_opt
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : vold, cfnorm, nncf, b2rr, b2volOnCf, b2data, b2dataOnCf
#ifdef TEST_MAP
      use b2mod_transport_namelist
     & , only : parm_dna
#endif
      use b2mod_subsys
      implicit none
      type(mapping), intent (in)    :: mpg
      type(geometry), intent (in)   :: geo
      type(B2State), intent (in)    :: st
      type(B2StateExt), intent (in) :: st_ext
      real(kind=R8) :: BoRiS, J(nncf)
      integer :: nCv, nFc, nVx, ns
      integer :: icf, ic1, ic2, iCv, n1, n2,
     &           ifc, is, iicf
      logical, save :: first_call=.true.
      logical :: inrange
      real (kind=R8) :: qq, prior, lll, lll_cum
      real (kind=R8) :: gradr(nCv), funv(nVx), pz(nCv), rz(nCv), cs(nCv)
#ifdef TEST_MAP
      integer, parameter :: nn=50
      real (kind=R8), save :: xx(nn), dd(nn), ff(nn)
      integer :: ii, ij, jj
#endif
      integer :: isigma, imean
      external calc_dist

      call subini('b2usr_cost_function')

      if (first_call .and. ncf .gt. 0) then
        call xertst(ncf.eq.mpg%nCf,'ncf=/mpg%nCf')
        call xertst(mpg%isClassicalGrid.eq.1,' Cost function CVs and '//
     &  'FCs ordering  and interpolation in case of unstructured grid'//
     &  ' needs to be revised')
#ifdef TEST_MAP
        ij = 1
        do ii=1,5
          do jj=1,10
            xx(ij) = real(jj,kind=r8) 
            ij = ij + 1
          end do
        end do
#else
        cfnorm = 0.0_R8
        vold = 0.0_R8
        do icf = 1, ncf
          if (cfread(icf)) then
! WARNING! For now we assume that the data read from the cost function is only radially dependent
            write (*,*) 'WARNING, b2usr_cost_function assumes data '//
     &                  'read for CF is only radially dependent'
            ic1 = mpg%cfRegP(icf,1)
            ic2 = ic1 + mpg%cfRegP(icf,2) - 1
            n1 = mpg%cfRegP(icf,2)
            n2 = ncfdata(icf)
            call calc_dist(mpg,geo,mpg%cfReg(ic1:ic2),n1,icsepomp-1, !csc icsep-1 because non-internal cells have been excluded!
     &        b2rr(icf,1:n1))
            b2data(1:n1) = geo%cvVol(mpg%cfReg(ic1:ic2))
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2volOnCf(icf,1:n2))
            vold(icf) = sum(b2volOnCf(icf,1:n2))
            do iCv = 1, n2
              cfnorm(icf) = cfnorm(icf) + 
     &                      b2volOnCf(icf,iCv)*cfdata(icf,2,iCv)
            end do
            if (mapToOMP(icf)) then
            ! FIXME to be considerably improved!
            !for now assume radial data rigidly moved at the OMP (like in ST code)
            ! simply use as b2rr the OMP radial coordinates.
              call xertst(mpg%isClassicalGrid.eq.1,
     &         'Mapping to OMP need revision for extended grids')
              ! check on, and then use, only internal Cvs in OMP list
              call xertst(n1.eq.nomp-2,
     &         'number of Cv in CF does not equal Cvs in OMP')
              call calc_dist(mpg,geo,omp(2:nomp-1),n1,icsepomp-1,
     &          b2rr(icf,1:n1))
            endif
          endif
        enddo 
#endif
        first_call = .false.
      endif

      J = 0.0_R8
#ifdef TEST_MAP
      isigma = 1
      imean = 1
      iicf = 2
      do ii=1,nn
        ff(ii) = parm_dna(1)*(xx(ii)**1.5)
      end do
      call calc_log_prior(prior, inrange)
      n1 = ncfdata(iicf)
      call calc_loglikelihood(n1, ff, cfdata(iicf,2,1:n1),
     & lll_cum, isigma, imean)
      if (inrange) then
        J(1) = -(prior + lll_cum)*cfweight(1)
      else
        J(1) = inf_opt
      endif
#else
      if (ncf .gt. 0) then
        do icf = 1, ncf
          ic1 = mpg%cfRegP(icf,1)
          ic2 = ic1 + mpg%cfRegP(icf,2) - 1
          n1 = mpg%cfRegP(icf,2)
          n2 = ncfdata(icf)
          select case (cftype(icf))
          case (0) !sum of other cost function. It is skipped at this stage because it MUST be J(1) for optimization
          case (1) !electron density on desired CVs
            b2data(1:n1) = st%dv%ne(mpg%cfReg(ic1:ic2))/1.0e19_R8
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (2) !electron temperature on desired CVs
            b2data(1:n1) = st%pl%te(mpg%cfReg(ic1:ic2))/ev
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (3) !ion temperature on desired CVs
            b2data(1:n1) = st%pl%ti(mpg%cfReg(ic1:ic2))/ev
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (4) !max electron temperature on desired CVs
            do iCv = ic1, ic2
              J(icf) = max(J(icf),st%pl%te(mpg%cfReg(iCv)))
            enddo
            J(icf) = J(icf)/ev
          case (5) !max heat flux density on desired FCs
            do iFc = ic1, ic2
              qq = 0.0
              iCv = min(mpg%fcCv(mpg%cfReg(iFc),1),
     &                  mpg%fcCv(mpg%cfReg(iFc),2))
              if (iCv.gt.mpg%nCi) iCv = max(mpg%fcCv(mpg%cfReg(iFc),1),
     &          mpg%fcCv(mpg%cfReg(iFc),2))
              qq = st%dv%fhe(mpg%cfReg(iFc),0) + 
     &             st%dv%fhi(mpg%cfReg(iFc),0) + 
     &             (1-BoRiS)*st%dv%fne(mpg%cfReg(iFc),0)*st%pl%te(iCv)
              do is = 0, ns-1
                qq = qq + st%dv%fna(mpg%cfReg(iFc),0,is)*
     &               (st%rt%rpt(iCv,is)*ev + 
     &               (1-BoRiS)*(st%dv%kinrgy(iCv,is)+st%pl%ti(iCv)))
              end do
              qq = qq*mpg%cfFcOr(iFc)/geo%fcS(mpg%cfReg(iFc))
              J(icf) = max(J(icf),qq)*cfweight(icf)
            enddo
          case (6) !loglikelihood + log prior for Bayesian MAP estimate
            call calc_log_prior(prior, inrange)
            write(*,*) 'MAP log prior value:', prior
            lll_cum = 0.0_R8
            isigma = 1
            imean = 1
            do iicf = cfstart(1), cfend(1)
              n1 = mpg%cfRegP(iicf,2)
              n2 = ncfdata(iicf)
              ic1 = mpg%cfRegP(iicf,1)
              ic2 = ic1 + n1 - 1
              select case (cftype(iicf))
              case (1) !electron density
                b2data(1:n1) = st%dv%ne(mpg%cfReg(ic1:ic2))/1.0e19_R8
              case (2) !electron temperature
                b2data(1:n1) = st%pl%te(mpg%cfReg(ic1:ic2))/ev
              case (3) !ion temperature
                b2data(1:n1) = st%pl%ti(mpg%cfReg(ic1:ic2))/ev
              case (7) !electron density radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%dv%ne, 
     &            funv, gradr)
                b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/1.0e19_R8
              case (8) !electron temperature radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%te, 
     &            funv, gradr)
                b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/ev
              case (9) !ion temperature radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%ti, 
     &            funv, gradr)
                b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/ev
              case (10) !ion saturation current (kA/m**2)
                call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
                call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te,
     &            st%pl%ti, pz, st_ext)
                cs   = sqrt(pz/rz)
                do icv=1, n1
                  b2data(icv) = st%dv%ne(mpg%cfReg(ic1+icv-1))*
     &              cs(mpg%cfReg(ic1+icv-1))*eV/1.0e3_R8
                enddo
              case default
               write(*,*) cftype(iicf)
               call xerrab ('cftype for sum out of bounds for cftype=6')
              end select
              call interp1d(n1,n2,b2rr(iicf,1:n1),cfdata(iicf,1,1:n2),
     &          b2data(1:n1),b2dataOnCf(1:n2))
              call calc_loglikelihood(n2, b2dataOnCf(1:n2), 
     &          cfdata(iicf,2,1:n2), cfdata(iicf,3,1:n2), lll, isigma,
     &          imean)
              lll_cum = lll_cum + lll
              isigma = isigma + 1
              imean = imean + 1
            end do
            if (inrange) then
              J(icf) = -(prior + lll_cum)*cfweight(icf)
            else
              J(icf) = inf_opt
            endif
          case (7) !electron density radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%dv%ne, 
     &        funv, gradr)
            b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/1.0e19_R8
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (8) !electron temperature radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%te, 
     &        funv, gradr)
            b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/ev
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (9) !ion temperature radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%ti, 
     &        funv, gradr)
            b2data(1:n1) = gradr(mpg%cfReg(ic1:ic2))/ev
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case (10) !ion saturation current (kA/m**2)
            call b2xprz (nCv, ns, mp, am, st%pl%na, rz, st_ext)
            call b2xppz (nCv, ns, st%dv%ne, st%pl%na, st%pl%te,
     &        st%pl%ti, pz, st_ext)
            cs   = sqrt(pz/rz)
            ! Calculate saturation current as ne*cs*eV
            do icv = 1, n1
              b2data(icv) = st%dv%ne(mpg%cfReg(ic1+icv-1))*
     &          cs(mpg%cfReg(ic1+icv-1))*eV/1.0e3_R8
            enddo
            call interp1d(n1,n2,b2rr(icf,1:n1),cfdata(icf,1,1:n2),
     &        b2data(1:n1),b2dataOnCf(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + b2volOnCf(icf,iCv)*
     &        (b2dataOnCf(iCv)-cfdata(icf,2,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)/
     &               2.0_R8
          case default
            write(*,*) cftype(icf)
            call xerrab ('cftype out of bounds')
          end select
        enddo
        ! now sum of cost functions for cftype=0
        if (cftype(1).eq.0) then
          do icf = cfstart(1), cfend(1)
            J(1) = J(1) + J(icf)
          end do
        endif
      endif  
#endif   
      call subend ()
      return

      end subroutine b2usr_cost_function

