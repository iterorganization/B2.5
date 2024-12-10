      subroutine b2tnce (mpg, geo, switch, nCv, te, 
     &  ne2, lnlam, collisnumf, collisnumc)
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
      integer nCv
      real (kind=R8) :: ne2(nCv), lnlam(nCv), te(nCv)
      type (mapping), intent (in) :: mpg
      type (geometry), intent (in) :: geo
      type (switches), intent (in) :: switch
      real (kind=R8) :: wrk1(mpg%nFc,2), tauf(mpg%nFc),
     *   vte(mpg%nFc), collisnum(mpg%nFt),
     *   tau(nCv), t0, t1, ctaue
      integer :: iFt, k, iFc, iCv
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: collisnumf(mpg%nFc),
     *  collisnumc(mpg%nCv)

      t1 = (4.0_R8*pi*eps0/qe)/qe
      t0 = (3.0_R8/(4.0_R8*sqrt(2.0_R8*pi)))*t1*(t1*ev*sqrt(ev))
      if (switch%b2tqce_model.eq.1.or.switch%b2tqce_model.eq.3) then
        ctaue = t0*sqrt(me)                                                !xpb
      else if (switch%b2tqce_model.eq.2) then
        ctaue = 3.5e11_R8             ! lower precision as per SOLPS4.0
      endif
      do iCv=1, nCv
        t0 = te(iCv)/ev
        tau(iCv) = (ctaue/lnlam(iCv))*t0*sqrt(t0)/ne2(iCv)  
      enddo

      wrk1 = 1.0_R8
      call intface(nCv,mpg%nFc,mpg%fcCv,wrk1,tau,tauf)
      call intface(nCv,mpg%nFc,mpg%fcCv,wrk1,sqrt(te/me),vte)
      collisnum = 0.0_R8
      do iFt=1,mpg%nFt
        do k=mpg%ftFcP(iFt,1),mpg%ftFcP(iFt,1)
     &     +mpg%ftFcP(iFt,2)-1
          iFc = mpg%ftFc(k)
          collisnum(iFt) = collisnum(iFt)+1/(tauf(iFc)*vte(iFc))*
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
        call my_out_us(70,mpg%nCv,1,tau,'b2tnce_tau')
        call my_out_us(70,mpg%nFt,1,collisnum,'b2tnce_numcole')
        call my_out_us(70,mpg%nFc,1,collisnumf,'b2tnce_numcolef')
        call my_out_us(70,mpg%nCv,1,collisnumc,'b2tnce_numcolec')
      endif
      end subroutine