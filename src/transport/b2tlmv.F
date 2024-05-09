      subroutine b2tlmv (nCv, nFc, nVx, isb, cflmv, switch, geo, mpg,
     &    nb, tb, ub, vsbfx, cvsbx, cvsbxhz, fl)                       !srv 09.01.01
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tlmv
#ifndef NAGFOR
      use b2mod_ad
     & , only : my_out_folder
#endif
      use b2mod_subsys
      implicit none
      integer nCv, nFc, nVx, isb
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  cflmv, nb(nCv), tb(nCv), ub(nCv), vsbfx(nFc), cvsbx(nFc), 
     *  cvsbxhz(nFc), fl(nFc)
*     ------------------------------------------------------------------
*     B2TLMV applies a flux limit to the viscosity.
*     ------------------------------------------------------------------
      integer iFc
      character charns*3                                             !srv 12.04.04
      real (kind=R8) ::
     *   qcl, qfl, tfl, ubv(nVx), dub(nFc), cflmv_loc(nFc), nbf(nFc)
      intrinsic abs
      external xertst
      external b2xvsg
*   ..initialisation
*     ------------------------------------------------------------------
      call subini ('b2tlmv')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test cflmv
      call xertst (0.0_R8.le.cflmv, 'faulty argument cflmv')
*   ..extensive tests on first few calls
      if (ncall_b2tlmv.lt.3) then
*    ..test sign of nb, tb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
       call b2xvsg (nCv, tb, 1, 'tb', '.gt.')
      endif

      ! Init local flux limit to constant value
      cflmv_loc = cflmv

      ! Enhanced flux limit for low-density regions
      if (switch%b2tlmv_far_sol .eq. 1) then
        ! Attempt enhanced flux limit for low-density regions
        call intface(nCv,nFc,mpg%fcCv,geo%fcVol,nb,nbf)
        do iFc = 1, nFc
          if (nbf(iFc) .le. switch%b2tlmv_ni_min) then
           ! Set at min value
           cflmv_loc(iFc) = switch%b2tlmv_ni_min
          elseif (nbf(iFc) .lt. switch%b2tlmv_ni_max ) then
           ! Interpolate for smooth transition
           cflmv_loc(iFc) = 
     .      (cflmv*(nbf(iFc) - switch%b2tlmv_ni_min) + 
     .       switch%b2tlmv_cflmv_min*(switch%b2tlmv_ni_max - nbf(iFc)))/
     .      (switch%b2tlmv_ni_max - switch%b2tlmv_ni_min)
          endif
        enddo
      endif

*   ..compute differences of ub at faces
      call diff_p(nCv, nFc, nVx, 0, geo, mpg, ub, ubv, dub)

*   ..main computation
      if(switch%b2tlmv_style.eq.0) then                              !srv 30.05.08
       do iFc = 1, nFc
          qcl = -cvsbx(iFc)*dub(iFc)
          qfl = cflmv_loc(iFc)*abs(geo%fcPbs(iFc))*                  !srv 11.01.13
     &     (nb(mpg%fcCv(iFc,1))*tb(mpg%fcCv(iFc,1))+
     &      nb(mpg%fcCv(iFc,2))*tb(mpg%fcCv(iFc,2)))/2.0_R8  
          if (qfl.gt.0.0_R8) then 
            tfl = 1.0_R8 + abs(qcl/qfl)
          else
            tfl = 1.0_R8
          endif
          cvsbx(iFc) = cvsbx(iFc)/tfl
          cvsbxhz(iFc) = cvsbxhz(iFc)/tfl                            !srv 02.07.08
          vsbfx(iFc) = vsbfx(iFc)/tfl
          fl(iFc) = tfl
       enddo
csw 18sep2014 TESTING only (SOLPS4 discr), i.e. switch off mom.flux
c             limit
c      else if(style.eq.2) then
c       fl = 1.0_R8
      else                                                           !srv 30.05.08
       do iFc = 1, nFc
          if(.not.(.not.is_neutral(isb) .and.
     &              mpg%cvOnClosedSurface(mpg%fcCv(iFc,1)) .and.
     &              mpg%cvOnClosedSurface(mpg%fcCv(iFc,2)))) then    !srv 12.04.08
           qcl = -cvsbx(iFc)*dub(iFc)
           qfl = cflmv_loc(iFc)*abs(geo%fcPbs(iFc))*                 !srv 11.01.13
     &          (nb(mpg%fcCv(iFc,1))*tb(mpg%fcCv(iFc,1))+
     &           nb(mpg%fcCv(iFc,2))*tb(mpg%fcCv(iFc,2)))/2.0_R8
           if (qfl.gt.0.0_R8) then
             tfl = 1.0_R8 + abs(qcl/qfl)
           else
             tfl = 1.0_R8
           endif
           cvsbx(iFc) = cvsbx(iFc)/tfl
           cvsbxhz(iFc) = cvsbxhz(iFc)/tfl
           vsbfx(iFc) = vsbfx(iFc)/tfl
           fl(iFc) = tfl
         else
           fl(iFc) = 1.0_R8
         endif
       enddo
      endif                                                          !srv 30.05.08

      if(switch%b2tlmv_iout.ne.0) then                               !srv 12.04.04 {
       write (charns,'(i3.3)') isb
       call my_out_us(70,nFc,1,fl,'cflmv-1'//charns)
      endif                                                          !srv 12.04.04 }

      ncall_b2tlmv = ncall_b2tlmv+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2tlmv
