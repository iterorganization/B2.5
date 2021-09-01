
      subroutine b2ye_sterr (meanbatch, nn, eps_s, tmax, T1, T2,
     .                      meansample)

      use b2mod_types

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YE_STERR computes the statistical error of the input array of
*     batch means
*
*  2. input/output variables
*
*     meanbatch(:): input variable, 1-D array with batch means
*     NN: number of batches
*     nototdt: total number of iterations for averaging
*     eps_s: output variable, statistical error
*     eps_s1: output variable, statistical error
*     eps_s2: output variable, statistical error
*     T: local variable, scalar, correlation time
*
*  3. description
*
*     For a detailed description of the equations see K. Ghoos et al.,
*     Numerical error estimation in Random Noise coupled plasma edge
*     simulations in nuclear fusion reactors, M&C2017
*
*     This routine performs the following steps.
*     ..compute T based on batch means
*     ..compute T based on autocorrelation
*     ..compute statistical error based on maximum value for T
*
*-----------------------------------------------------------------------
*.declarations

      implicit none

      integer, intent(in) :: nn
      real (kind=R8), intent(in) :: meanbatch(NN)
      real (kind=R8), intent(out) :: eps_s, meansample, tmax, T1, T2

      real (kind=R8), allocatable, save :: meangroup(:),
     .                             meangroup2(:), meanbatch2(:),
     .                             rho(:), tau(:), logrho(:),
     .                             rhonew(:), taunew(:)
      real (kind=R8) :: STAT, meantot, varbatch, eps_s1, eps_s2,
     .          meansamplebis, varbatchbis, vargroup, fit(20)
      integer :: I, M, IFAIL, ncutoff, nfit, ncutoffnew
      logical :: check
      logical, save :: first_call = .true.
      external G13ABF, G02CAF

      if (.not.allocated(meangroup)) then
        allocate (meangroup(nn))
        allocate (meangroup2(nn))
        allocate (meanbatch2(nn))
        allocate (rho(nn))
        allocate (logrho(nn))
        allocate (tau(nn))
        meangroup = 0.0_R8
        meangroup2 = 0.0_R8
        meanbatch2 = 0.0_R8
        rho = 0.0_R8
        logrho = 0.0_R8
        tau = 0.0_R8
      end if


*** Compute T based on batch means ensembles
*     compute mean and variance of samples and meanbatches
*            (meansample = mean of meanbatch !)
      meansample = sum(meanbatch)/NN
      meanbatch2 = meanbatch**2
      varbatch = (SUM(meanbatch2(1:NN))-NN*meansample**2)/(NN-1)

*     compute mean and variance of group of batches
      M = NN/30
      if (first_call) then
        write (*,*) ' Size of groups overlapping batches: ',M
        if (NN.lt.30) then
          write (*,*)' This should be at least 1, please continue',
     ,               'averaging procedure'
          call xerrab('Insufficient number of batches')
        endif
        if (M.lt.30) then
          write (*,*) ' Preferably, this should be larger than 30, ',
     ,                'consider continue the averaging'
        endif
        first_call = .false.
      endif
      DO I=1,NN-M+1
        meangroup(I)= SUM(meanbatch(I:I+M-1))/M
        meangroup2(I)=meangroup(I)**2
      ENDDO
      meantot = sum(meangroup(1:NN-M+1))/(NN-M+1)
      vargroup = (SUM(meangroup2(1:NN-M+1))-(NN-M+1)*meantot**2)/(NN-M)

*** Compute correlation time for batches

      T1 = vargroup*M/varbatch

*** Compute T based on autocorrelation
      ncutoff = 30
      ifail = 0
      check = .true.
      nfit = 0
      CALL G13ABF(meanbatch, NN, NN-ncutoff, meansamplebis, varbatchbis,
     .            rho, STAT, IFAIL )
*     --- for fit only take samples before the first one that is
*         within the 1 sigma band ---
      DO I=1,NN-ncutoff
        IF(abs(rho(I)).LT.(1.0_R8/SQRT(real(NN))).and.check) THEN
          nfit=I
          check=.false.
        ENDIF
        logrho(I)=log(rho(I))
        tau(I)=real(I)
      ENDDO

      IF (nfit.gt.2) then
        CALL G02CAF(nfit,tau,logrho,fit,ifail)
        ncutoffnew = int((log(0.001_R8)-fit(7))/fit(6))
        allocate (rhonew(ncutoffnew))
        allocate (taunew(ncutoffnew))
        taunew = (/ ( i, i=1, ncutoffnew ) /)
        rhonew(1:ncutoffnew)=exp(fit(6)*taunew(1:ncutoffnew)+fit(7))
        T2=1+2*sum(rhonew(1:ncutoffnew))
        deallocate (rhonew)
        deallocate (taunew)
      ELSE
        T2=0.0_R8
      ENDIF


*** Compute statistical error

      tmax = max(T1,T2)
      eps_s = sqrt(varbatch*tmax/NN)
      eps_s1 = sqrt(varbatch*T1/NN)
      eps_s2 = sqrt(varbatch*T2/NN)


      return
      end subroutine b2ye_sterr


* -----------------------------------------------------------------------------



