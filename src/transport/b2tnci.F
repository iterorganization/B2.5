      subroutine b2tnci (mpg, geo, switch, nCv, ns, ti, te, rz2, ne, 
     &  ne2, na, lnlam, ismain, st_ext, collisnumf, collisnumc)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_plasma
      use b2us_map
      use b2us_geo
      use b2mod_switches
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
    !   use b2mod_ad
    !  & , only : ncall_b2ttia
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns, ismain
      real (kind=R8) ::
     &   ti(nCv), ne2(nCv), lnlam(nCv), rz2(nCv,0:ns-1), na(nCv,ns-1),
     &   ne(nCv), te(nCv)
      type (B2StateExt), intent (in) :: st_ext
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      type (switches), intent (in) :: switch
      real (kind=R8) ::
     *   tauia(nCv,0:ns-1), rz(nCv), pz(nCv), wrk1(mpg%nFc,2), 
     *   tauiaf(mpg%nFc), csf(mpg%nFc), collisnum(mpg%nFt)
      integer :: iFt, k, iFc, iCv
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: collisnumf(mpg%nFc),
     *  collisnumc(mpg%nCv)

      call b2ttia (nCv, ns, ti, rz2, ne2, lnlam, tauia)
      call b2xprz (nCv, ns, mp, am, na, rz, st_ext)
      call b2xppz (nCv, ns,  ne, na, te, ti, pz, st_ext)
      wrk1 = 1.0_R8
      call intface(nCv,mpg%nFc,mpg%fcCv,wrk1,tauia(:,ismain),tauiaf)
      call intface(nCv,mpg%nFc,mpg%fcCv,wrk1,sqrt(pz/rz),csf)
      collisnum = 0.0_R8
      do iFt=1,mpg%nFt
        do k=mpg%ftFcP(iFt,1),mpg%ftFcP(iFt,1)
     &     +mpg%ftFcP(iFt,2)-1
          iFc = mpg%ftFc(k)
          collisnum(iFt) = collisnum(iFt)+1/(tauiaf(iFc)*csf(iFc))*
     &     (geo%fcHc(iFc,1)+geo%fcHc(iFc,2))*
     &      geo%fcBb(iFc,3)/geo%fcBb(iFc,0)
        enddo
      enddo

      collisnumf = 0.0_R8
      do iFc = 1,mpg%nFc
        iFt = mpg%cvFt(mpg%fcCv(iFc,1))
        if (iFt .eq. 0)  iFt = mpg%cvFt(mpg%fcCv(iFc,2))
        collisnumf(iFc) = collisnum(iFt)
      enddo

      collisnumc = 0.0_R8
      do iCv = 1,mpg%nCv
        iFt = mpg%cvFt(iCv)
        if (iFt.eq.0) cycle
        collisnumc(iCv) = collisnum(iFt)
      enddo

      if (switch%iout_b2wdat.eq.4) then
        call my_out_us(70,mpg%nFt,1,collisnum,'b2tnci_numcoli')
        call my_out_us(70,mpg%nFc,1,collisnumf,'b2tnci_numcolif')
        call my_out_us(70,mpg%nCv,1,collisnumc,'b2tnci_numcolic')
      endif
      end subroutine