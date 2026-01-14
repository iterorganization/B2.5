*-----------------------------------------------------------------------
*.specification

      subroutine b2usco (nCv, nFc, nVx, isb, switch, geo, mpg,    !srv 13.10.06 26.09.12
     &   nregionv, solvereg, itcnt, rxg, nb, pb, snbv, flob, conb,
     &   rescb, corpb, aa, name)                                  !srv 18.05.02
      use b2mod_types
      use b2mod_switches
      use b2us_geo
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2usco
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, itcnt, isb, nregionv
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      logical solvereg(0:nregionv)
      real (kind=R8) ::
     *   rxg, nb(nCv), pb(nCv), snbv(nCv), flob(nFc,0:1),
     &   conb(nFc,0:1,0:4), rescb(nCv)
      character*(*) name
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corpb(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa(mpg%nCmxNv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USCO computes a correction towards the solution of the
*     continuity equation for one atomic species on an unstructured
*     mesh.
*
*-----------------------------------------------------------------------
!                 2_____3_____4                     5___________6
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -1_____0_____1                     |_____|_____|
!                 |     |     |                     |     |     |
!                 |     |     |                     |     |     |
!                -4____-3____-2                    -6_____|____-5

!                (ix-1,iy+1)______(ix,iy+1)____(ix+1,iy+1)
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy  )______(ix,iy)______(ix+1,iy  )
!                     |              |              |
!                     |              |              |
!                     |              |              |
!                (ix-1,iy-1)______(ix,iy-1)____(ix+1,iy-1)
*.declarations

*   ..local variables
      integer iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv, ireg, k
c      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !srv 06.07.06
c      integer, save :: iadd_9 = 1                                     !srv 13.10.06
      real (kind=R8) ::
     *   t1, t2, t3, t11, t22, wrk(nCv)
      character chns*3, charnv*2
      logical lnonzero
*   ..procedures
      intrinsic abs
      external xertst, ipgetr, sfill
      external b2xvsg                                                  !srv 02.01.07

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usco')
*   ..set internal parameters on first call
      if (ncall_b2usco.eq.0) then                                     !srv 06.07.06 {
c       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)
c       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)
c       call ipgeti ('b2usco_add_9', iadd_9)                          !srv 13.10.06
      endif                                                           !srv 06.07.06 }
*   ..test nCv, nFc, nVx
      call xertst (0.lt.nCv.and.0.lt.nFc.and.0.lt.nVx,
     .  'faulty argument nCv, nFc, nVx')
*   ..test itcnt, rxg
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.lt.rxg, 'faulty argument rxg')
*   ..extensive tests on first few calls
      if (ncall_b2usco.lt.3) then
*    ..test sign of nb, pb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
       call b2xvsg (nCv, pb, 1, 'pb', '.gt.')
*    ..test sign of snbv
       call b2xvsg (nCv, snbv, 1, 'snbv', '.le.')
      endif

*   ..initialize matrix
      aa = 0.0_R8

*   ..compute the matrix
! $AD II-LOOP
      do iCv = 1, nCv
        ireg  = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
          if (switch%b2mndt_style.lt.2) then                          !srv 13.10.06 {
            aa(mpg%cvNvP(iCv,1)) = (-snbv(iCv)*nb(iCv)/pb(iCv)+
     &       abs(rescb(iCv))/(rxg*pb(iCv)))*rscal(pb(iCv))
          else
            aa(mpg%cvNvP(iCv,1)) = (-snbv(iCv)*nb(iCv)/pb(iCv))*
     &          rscal(pb(iCv))
          endif                                                       !srv 13.10.06 }
        else
          aa(mpg%cvNvP(iCv,1)) = rscal(pb(iCv))
        endif
      enddo
! $AD II-LOOP
      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flob(:,0), conb(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flob(iFc,0) - conb(iFc,0,0) - conb(iFc,0,2)
        t2 =   0.5_R8*flob(iFc,0) - conb(iFc,0,0) - conb(iFc,0,1)
        if (solvereg(mpg%cvReg(iCv1))) then
         aa(mpg%cvNvP(iCv1,1)) = aa(mpg%cvNvP(iCv1,1)) -
     &     t1*rscal(pb(iCv1))
         where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                   mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv2)
          aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &     aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) +
     &     t2*rscal(pb(iCv2))
         end where
        endif
        if (solvereg(mpg%cvReg(iCv2))) then
         aa(mpg%cvNvP(iCv2,1)) = aa(mpg%cvNvP(iCv2,1)) -
     &     t2*rscal(pb(iCv2))
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv1)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t1*rscal(pb(iCv1))
         end where
        endif
*     ..contributions from flob(:,1), conb(:,1,:) ('9-point' contributions)
        iVx1 = mpg%fcVx(iFc,1) ! number of v1
        iVx2 = mpg%fcVx(iFc,2) ! number of v2
        t1 = - 0.5_R8*flob(iFc,1) - conb(iFc,1,0) - conb(iFc,1,2)
        t2 =   0.5_R8*flob(iFc,1) - conb(iFc,1,0) - conb(iFc,1,1)
        t11 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx1,2)
          t11 = t11 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx1,1)+inCv-1)
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1.0_R8/geo%vxVol(mpg%vxCvP(iVx2,1)+inCv-1)
        enddo
        ! contribution to equation of iCv1
        if (solvereg(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)-
     &      t3*t1*rscal(pb(iCv))
          end where
         enddo
         ! contribution from iVx2
         do inCv = 1, mpg%vxCvP(iVx2,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
          t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
          where (mpg%cvNv(mpg%cvNvP(iCv1,1):
     &                    mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) == iCv)
           aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1) =
     &      aa(mpg%cvNvP(iCv1,1):mpg%cvNvP(iCv1,1)+mpg%cvNvP(iCv1,2)-1)+
     &      t3*t2*rscal(pb(iCv))
          end where
        enddo
       endif
       ! contribution to equation of iCv2
       if (solvereg(mpg%cvReg(iCv2))) then
        ! contribution from iVx1
        do inCv = 1, mpg%vxCvP(iVx1,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx1,1) + inCv - 1)*t11)
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) +
     &     t3*t1*rscal(pb(iCv))
         end where
        enddo
        ! contribution from iVx2
        do inCv = 1, mpg%vxCvP(iVx2,2)
         iCv = mpg%vxCv(mpg%vxCvP(iVx2,1) + inCv - 1)
         t3  = 1.0_R8/(geo%vxVol(mpg%vxCvP(iVx2,1) + inCv - 1)*t22)
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) =
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) -
     &     t3*t2*rscal(pb(iCv))
         end where
        enddo
       endif
      enddo
*   ..check for empty matrix lines and replace with dummy equation
! $AD II-LOOP
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1, mpg%cvNvP(iCv,2)
          lnonzero = lnonzero.or.(aa(mpg%cvNvP(iCv,1)+inCv-1).ne.0.0_R8)
        enddo
        if (.not.lnonzero) then
          aa(mpg%cvNvP(iCv,1)) = rscal(pb(iCv))
          rescb(iCv) = 0.0_R8
        endif
        ireg = mpg%cvReg(iCv)
        if (.not.solvereg(ireg)) rescb(iCv) = 0.0_R8
      enddo

*   ..solve the correction equation
      call b2uxus (nCv, mpg, aa, itcnt, rescb, corpb, name,            !srv 18.05.02
     .             switch%use_PARDISO)
*   ..redo the scaling
! $AD II-LOOP
      do iCv = 1, nCv
        corpb(iCv) = corpb(iCv)*rscal(pb(iCv))
      enddo
      if (switch%b2usco_iout.ne.0) then                                !srv 17.06.02 {
        write (chns,'(i3.3)') isb
        do k=1,maxval(mpg%cvNvP(1:nCv,2))
          write (charnv,'(i2.2)') k
          do iCv = 1, nCv
            if (k.le.mpg%cvNvP(iCv,2)) then
              wrk(iCv) = aa(mpg%cvNvP(iCv,1)+k-1)
            else
              wrk(iCv) = 0.0_R8
            end if
          end do
          call my_out_us(70,nCv,0,wrk,name//'_aam'//charnv//'_'//chns)
        end do
        call my_out_us(70,nCv,0,corpb,name//'_corpb'//chns)
        call my_out_us(70,nCv,0,rescb,name//'_rescb'//chns)
        call my_out_us(70,nCv,0,nb,name//'_na'//chns)
      endif                                                            !srv 17.06.02 }

* ..return
      ncall_b2usco = ncall_b2usco+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usco

      contains

      function rscal(t0)
      implicit none
      real(kind=R8) :: rscal, t0
! if uncomment rscal(t0) = t0 then uncomment nb(ix,iy) = nb(ix,iy)+rxf*corpb(ix,iy)/pb(ix,iy) in b2upco
      rscal = t0                                                  !srv 13.10.06
!      rscal = 1.0_R8                                             !srv 13.10.06
      return
      end function rscal

      end subroutine b2usco
