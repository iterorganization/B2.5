      subroutine b2usr_cost_function(ns, mpg, st, geo, J)
      use b2mod_types
      use b2mod_constants
      use b2mod_user_namelist
      use b2us_map
      use b2us_plasma
      use b2us_geo
      implicit none
      type(mapping), intent (in)    :: mpg
      type(geometry), intent (in)   :: geo
      type(B2State), intent (in)    :: st
      real(kind=R8) :: BoRiS, J(nncf)
      integer :: ns
      integer :: icf, ic1, ic2, iCv, ncall, ndata1, ndata2, ii,
     &           ifc, ifc1, ifc2, is
      logical, save :: first_call=.true.
      real (kind=R8), allocatable, save :: ffdata(:,:)
      real (kind=R8), allocatable ::  xx1(:), yy1(:)
      real(kind=R8) :: xcheck(2), ycheck(2), eps, machsfr, qq
      external machsfr, subini, subend

      call subini('b2usr_cost_function')

      BoRiS = 0.0_R8 !!FIXME 

      if (first_call .and. ncf .gt. 0) then
        call xertst(ncf.eq.mpg%nCf,'ncf=/mpg%nCf')
        ndata1 = maxval(mpg%cfCvP(1:ncf,2),1)
        eps = machsfr()
        allocate(ffdata(ncf,ndata1))
        ffdata = 0.0_R8
        do icf = 1, ncf
          if (cfread(icf)) then
            ic1 = mpg%cfCvP(icf,1)
!csc first interpolate cfdata onto B2.5 grid.  !FIXME safety for non monotonic!!
! In the and ffdata will contain only the data, no information on x-y is needed anymore
            ndata1= ncfdata(icf) !number of data points read
            ndata2= mpg%cfCvP(icf,2) !number of CV in cost function
            allocate(xx1(ndata2),yy1(ndata2))
            do iCv = 1, ndata2
              ! evaluate x and y data of CV in cost function
              xx1(iCv) = geo%cvX(mpg%cfCv(ic1 + iCv -1))
              yy1(iCv) = geo%cvY(mpg%cfCv(ic1 + iCv -1))
            enddo
!csc check that x and y values are both monotonic, otherwise perform simple 1d interpolation
           xcheck = 0.0_R8
           ycheck = 0.0_R8
           do ii = 2, ndata1 
            xcheck(1) = xcheck(1) + cfdata(icf,1,ii)-cfdata(icf,1,ii-1)
            ycheck(1) = ycheck(1) + cfdata(icf,2,ii)-cfdata(icf,2,ii-1)
           end do
           do ii = 2, ndata2 
            xcheck(2) = xcheck(2) + xx1(ii)-xx1(ii-1)
            ycheck(2) = ycheck(2) + yy1(ii)-yy1(ii-1)
           end do
           if (xcheck(1).lt.eps .or. xcheck(2).lt.eps) then
!csc either xx or xx1 are constant --> perform 1d interpolation and issue warning
             write(*,*) ' Warning: xx or xx1 seem to be constant,'//
     >                  ' performing 1d interpolation for cost function'
             call interp1d(ndata1,ndata2,cfdata(icf,2,1:ndata1),
     &            yy1,cfdata(icf,3,1:ndata1),ffdata(icf,1:ndata2))
           elseif (ycheck(1).lt.eps .or. ycheck(2).lt.eps) then
!csc either yy or yy1 are constant --> perform 1d interpolation and issue warning
             write(*,*) ' Warning: yy or yy1 seem to be constant,'//
     >                  ' performing 1d interpolation for cost function'
             call interp1d(ndata1,ndata2,cfdata(icf,1,1:ndata1),
     &            xx1,cfdata(icf,3,1:ndata1),ffdata(icf,1:ndata2))
           else
             call xerrab 
     >         ('2d interpolation not yet available for cost function')
           endif
           deallocate(xx1,yy1)
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
          select case (cftype(icf))
          case (1) !electron density on desired CVs
            do iCv = ic1, ic2
              J(icf) = J(icf) + (st%dv%ne(mpg%cfCv(iCv))
     &                - ffdata(icf,iCv-ic1+1))**2 ! FIXME put normalization from user 
            enddo
            J(icf) = J(icf)/1.0E19
          case (2) !electron temperature on desired CVs
            do iCv = ic1, ic2
              J(icf) = J(icf) + (st%pl%te(mpg%cfCv(iCv))
     &                - ffdata(icf,iCv-ic1+1))**2
            enddo
            J(icf) = J(icf)/ev
          case (3) !max electron temperature on desired CVs
            do iCv = ic1, ic2
              J(icf) = max(J(icf),st%pl%te(mpg%cfCv(iCv)))
            enddo
            J(icf) = J(icf)/ev
          case (4) !max heat flux on desired FCs
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
              J(icf) = max(J(icf),qq) !!FIXME to be better difined with cffcor...
            enddo
          case default
            write(*,*) cftype(icf)
            call xerrab ('cftype out of bounds')
          end select
        enddo
      endif     
      call subend ()
      return
      end subroutine b2usr_cost_function
