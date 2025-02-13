      subroutine b2poleg(lun,iFt)
c
c----------------------------------------------------------------------
c***
c***  Produces a file with plasma parameters (density, temperature,
c***  parallel flow velocity, parallel heat flux) against the
c***  distance along the magnetic field along field tube index
c***  iFt - for Oleg
c***
c*** Output files:
c***
c***  oleg.dat with the data for Oleg
c----------------------------------------------------------------------
      use b2mod_types
      use b2mod_b2plot
      use b2mod_constants
      use b2us_data
      implicit none
      integer iFt
      integer lun
c
c*** Local variables
c
      integer i, j, iCn, iCp, iCv, iFc, iFn, iFt1, iFt2, k
      real (kind=R8) :: x, y, z, hlp_je, hlp_q, hlp_qp, hlp_si
c----------------------------------------------------------------------
c*** Check the ig value
c
      if(iFt.le.0 .or. iFt.gt.mpg%nFt) then
        print *,'*** b2poleg: wrong value of iFt = ',iFt
        return
      end if
c
c*** Write the header
c
      open(lun,file='oleg.dat')
      write (lun,710) 'x','Xpol','ne','Te','qpol_e','qpar_e',
     .                'S_ion','Zeff','jpar_e'
  710 format(6x,a1,9x,a4,9x,a2,10x,a2,8x,a6,6x,a6,
     ,       6x,a5,8x,a4,7x,a6/)
c
      x=0.0_R8
      z=0.0_R8
      i=1
      iCv=mpg%ftCv(mpg%ftCvP(iFt,1))
      iCp=iCv
      do while (i.le.mpg%ftCvP(iFt,2))
        iCn=0
        do j=1,mpg%cvFcP(iCv,2)
          iFc=mpg%cvFc(mpg%cvFcP(iCv,1)+j-1)
          iFt1=mpg%cvFt(mpg%fcCv(iFc,1))
          iFt2=mpg%cvFt(mpg%fcCv(iFc,2))
          if (iFt1.eq.iFt.and.iFt2.eq.iFt) then
            if (mpg%fcCv(iFc,1).eq.iCv.and.
     &          mpg%fcCv(iFc,2).ne.iCp) then
              iCn=mpg%fcCv(iFc,2)
              iFn=iFc
              iCp=iCv
            else if (mpg%fcCv(iFc,2).eq.iCv.and.
     &          mpg%fcCv(iFc,1).ne.iCp) then
              iCn=mpg%fcCv(iFc,1)
              iFn=iFc
              iCp=iCv
            end if
          end if
        end do
        i=i+1
        iCv=iCn
        if (iCn.gt.mpg%nCi.or.iCn.eq.0) cycle
        x=x+0.5_R8*geo%cvHx(iCn)
        y=0.5_R8*geo%cvHx(iCn)/pitch(geo,iCn)
        z=z+y
        hlp_q=state%dv%fhe(iFn,0)/area(iFn)
        hlp_qp=hlp_q/pitch(geo,iCn)
        hlp_si=0.0_R8
        hlp_je=0.0_R8
        do k=0,state%ns-1
          hlp_je=hlp_je+state%dv%fna(iFn,0,k)*state%rt%rza(iCn,k)
          hlp_si=hlp_si+state%srw%b2stbr_sna(iCn,k)
        end do
        do k=0,state_ext%ns-1
          hlp_je=hlp_je+state_ext%fa(iFn,0,k)*state%rt%rza(iCn,k)
        end do
        hlp_je=hlp_je/(area(iFn)*pitch(geo,iCn))
        hlp_si=hlp_si/geo%cvVol(iCn)
c
        write(lun,700) z,x,state%dv%ne(iCn),state%pl%te(iCn)/ev,
     ,    hlp_q,hlp_qp,hlp_si,zeff(iCn),hlp_je
        x=x+0.5_R8*geo%cvHx(iCn)
        z=z+y
  700   format(1p,48e12.4)
      end do
      close (lun)
      return
c
      contains

      function area(iFc)
      implicit none
      integer, intent(in) :: iFc
      real(kind=R8) :: area
      area = geo%fcS(iFc)*geo%fcQgam(iFc,0)
      return
      end function area

      end subroutine b2poleg
