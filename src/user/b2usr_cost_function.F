      subroutine b2usr_cost_function(nCv, nFc, nVx, ns, mpg, st, geo, J)
      use b2mod_types
      use b2mod_constants
      use b2mod_user_namelist
      use b2us_map
      use b2us_plasma
      use b2us_geo
      use b2mod_par_opt
     & , only : npar_opt
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : vold, cfnorm, nncf, rrdata, vvdata, TMP1, TMP2
#ifdef TEST_MAP
      use b2mod_transport_namelist
     & , only : parm_dna
#endif
      implicit none
      type(mapping), intent (in)    :: mpg
      type(geometry), intent (in)   :: geo
      type(B2State), intent (in)    :: st
      real(kind=R8) :: BoRiS, J(nncf)
      integer :: nCv, nFc, nVx, ns
      integer :: icf, ic1, ic2, iCv, ncall, n1, n2, ii,
     &           ifc, ifc1, ifc2, is, iicf
      logical, save :: first_call=.true.
      logical :: inrange
      real (kind=R8) :: qq, prior, lll, lll_cum
      real (kind=R8) :: gradr(nCv), funv(nVx)
#ifdef TEST_MAP
      integer, parameter :: nn=50
      real (kind=R8), save :: xx(nn), yy(nn), ff(nn)
      integer :: ij, jj
#endif
      integer :: nele_jac, ncon, isigma
      external subini, subend

      call subini('b2usr_cost_function')

      BoRiS = 0.0_R8 !!FIXME 

      if (first_call .and. ncf .gt. 0) then
        call xertst(ncf.eq.mpg%nCf,'ncf=/mpg%nCf')
#ifdef TEST_MAP
        ij = 1
        do ii=1,5
          do jj=1,10
            xx(ij) = real(jj,kind=r8) 
            ij = ij + 1
          end do
        end do
#endif    
        if (cftype(1).eq.6) then
          call xertst(npar_opt+nsigma.gt.0,'cftype 6 requires at '//
     &    'least one sigma in b2.user.parameters')
          call xertst(cfend(1)-cfstart(1)+1.eq.nsigma, 'Each cost '//
     &    'function summed in cftype 6 must have a sigma value '//
     &    'assigned in b2.optimization.parameters')
        endif
        cfnorm = 0.0_R8
        vold = 0.0_R8
        do icf = 1, ncf
          if (cfread(icf)) then
! WARNING! For now we assume that the data read from the cost function is only radially dependent
            write (*,*) 'WARNING, b2usr_cost_function assumes data '//
     &                  'read for CF is only radially dependent'
            ic1 = mpg%cfCvP(icf,1)
            ic2 = ic1 + mpg%cfCvP(icf,2) - 1
            n1 = mpg%cfCvP(icf,2)
            n2 = ncfdata(icf)
            rrdata(icf,1:n1) = geo%cvX(mpg%cfCv(ic1:ic2))
            TMP1(1:n1) = geo%cvVol(mpg%cfCv(ic1:ic2))
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        TMP1(1:n1),vvdata(icf,1:n2))
            vold(icf) = sum(vvdata(icf,1:n2))
            do iCv = 1, n2
              cfnorm(icf) = cfnorm(icf) + 
     &                      vvdata(icf,iCv)*cfdata(icf,3,iCv)
            end do
          endif
        enddo 
        first_call = .false.
      endif

      J = 0.0_R8

      if (ncf .gt. 0) then
        do icf = 1, ncf
          ic1 = mpg%cfCvP(icf,1)
          ic2 = ic1 + mpg%cfCvP(icf,2) - 1
          ifc1 = mpg%cfFcP(icf,1)
          ifc2 = ifc1 + mpg%cfFcP(icf,2) - 1
          n1 = mpg%cfCvP(icf,2)
          n2 = ncfdata(icf)
          select case (cftype(icf))
          case (0) !sum of other cost function. It is skipped at this stage because it MUST be J(1) for optimization
          case (1) !electron density on desired CVs
            TMP1(1:n1) = st%dv%ne(mpg%cfCv(ic1:ic2))
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        TMP1(1:n1),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
          case (2) !electron temperature on desired CVs
            TMP1(1:n1) = st%pl%te(mpg%cfCv(ic1:ic2))/ev
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        TMP1(1:n1),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
          case (3) !ion temperature on desired CVs
            TMP1(1:n1) = st%pl%ti(mpg%cfCv(ic1:ic2))/ev
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        TMP1(1:n1),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
          case (4) !max electron temperature on desired CVs
            do iCv = ic1, ic2
              J(icf) = max(J(icf),st%pl%te(mpg%cfCv(iCv)))
            enddo
            J(icf) = J(icf)/ev
          case (5) !max heat flux on desired FCs
            do iFc = ifc1, ifc2
              qq = 0.0
              iCv = mpg%fcCv(mpg%cfFc(iFc),1)
              qq = st%dv%fhe(mpg%cfFc(iFc),0) + 
     &             st%dv%fhi(mpg%cfFc(iFc),0) + 
     &             (1-BoRiS)*st%dv%fne(mpg%cfFc(iFc),0)*st%pl%te(iCv)
              do is = 0, ns-1
                qq = qq + st%dv%fna(mpg%cfFc(iFc),0,is)*
     &               (st%rt%rpt(iCv,is)*ev + 
     &               (1-BoRiS)*(st%dv%kinrgy(iCv,is)+st%pl%ti(iCv)))
              end do
              qq = abs(qq)
              J(icf) = max(J(icf),qq)*cfweight(icf) !!FIXME to be better difined with cffcor...
            enddo
          case (6) !loglikelihood + log prior for Bayesian MAP estimate
#ifdef TEST_MAP
            isigma = 1
            iicf = 2
            do ii=1,nn
              ff(ii) = parm_dna(1)*(xx(ii)**1.5)
            end do
            call calc_prior(prior, inrange)
            n1 = ncfdata(iicf)
            call calc_loglikelihood(n1, ff, cfdata(iicf,3,1:n1),
     &       lll_cum, isigma)
#else
            call calc_prior(prior, inrange)
            lll_cum = 0.0_R8
            isigma = 1
            do iicf = cfstart(1), cfend(1)
              n1 = mpg%cfCvP(iicf,2)
              n2 = ncfdata(iicf)
              ic1 = mpg%cfCvP(iicf,1)
              ic2 = ic1 + n2 - 1
              select case (cftype(iicf))
              case (1) !electron density
                TMP1(1:n1) = st%dv%ne(mpg%cfCv(ic1:ic2))
              case (2) !electron temperature
                TMP1(1:n1) = st%pl%te(mpg%cfCv(ic1:ic2))/ev
              case (3) !ion temperature
                TMP1(1:n1) = st%pl%ti(mpg%cfCv(ic1:ic2))/ev
              case (7) !electron density radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%dv%ne, 
     &            funv, gradr)
                TMP1(1:n1) = gradr(mpg%cfCv(ic1:ic2))
              case (8) !electron temperature radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%te, 
     &            funv, gradr)
                TMP1(1:n1) = gradr(mpg%cfCv(ic1:ic2))/ev
              case (9) !ion temperature radial gradient
                call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%ti, 
     &            funv, gradr)
                TMP1(1:n1) = gradr(mpg%cfCv(ic1:ic2))/ev
              case default
               write(*,*) cftype(iicf)
               call xerrab ('cftype for sum out of bounds for cftype=6')
              end select
              call interp1d(n1,n2,rrdata(iicf,1:n1),cfdata(iicf,1,1:n2),
     &          TMP1(1:n1),TMP2(1:n2))
              call calc_loglikelihood(n2, TMP2(1:n2), 
     &          cfdata(iicf,3,1:n2), lll, isigma)
              lll_cum = lll_cum + lll
              isigma = isigma + 1
            end do
#endif
            if (inrange) then
              J(icf) = -(log(prior) + lll_cum)*cfweight(icf)
            else
              J(icf) = inf_opt
            endif
          case (7) !electron density radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%dv%ne, 
     &        funv, gradr)
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        gradr(mpg%cfCv(ic1:ic2)),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
          case (8) !electron temperature radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%te, 
     &        funv, gradr)
            gradr = gradr/ev
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        gradr(mpg%cfCv(ic1:ic2)),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
          case (9) !ion temperature radial gradient
            call gradc_r(nCv, nFc, nVx, 0, geo, mpg, st%pl%ti, 
     &        funv, gradr)
            gradr = gradr/ev
            call interp1d(n1,n2,rrdata(icf,1:n1),cfdata(icf,1,1:n2),
     &        gradr(mpg%cfCv(ic1:ic2)),TMP2(1:n2))
            do iCv = 1, n2
              J(icf) = J(icf) + vvdata(icf,iCv)*
     &        (TMP2(iCv)-cfdata(icf,3,iCv))**2
            end do
            J(icf) = J(icf)*cfweight(icf)*vold(icf)/(cfnorm(icf)**2)
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
      call subend ()
      return
#ifdef TEST_MAP
      contains

      function normal_random(mu,ss)
      use b2mod_constants
      implicit none

      real(kind=R8):: normal_random, u(2), mu, ss, xx

      call random_number(u)
      xx = sqrt(-2.0_R8*log(u(1)))*cos(2.0_R8*pi*u(2))
      normal_random = mu + xx*ss
      return
      end function
#endif
      end subroutine b2usr_cost_function

