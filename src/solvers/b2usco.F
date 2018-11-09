*-----------------------------------------------------------------------
*.specification

      subroutine b2usco (nCv, nFc, nVx, isb, geo, mpg, nregionv,            !srv 13.10.06 26.09.12
     &   solvereg, itcnt, rxg, pcm1, nb, pb, snbv, flob, conb,
     &   pccb1, rescb, corpb, aa, name)                              !srv 18.05.02
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, itcnt, isb, nregionv
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      logical solvereg(0:nregionv)      
      real (kind=R8) ::
     *   rxg, pcm1, nb(nCv), pb(nCv), snbv(nCv), flob(nFc,0:1),
     &   conb(nFc,0:1,0:4), pccb1(nCv), rescb(nCv)
      character*(*) name
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   corpb(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   aa(mpg%nCmxNv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2usco computes a correction towards the solution of the
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
      integer ncall, iCv, iFc, iCv1, iCv2, iVx1, iVx2, inCv, ireg
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !srv 06.07.06
      integer, save :: iadd_9 = 1, iout = 0                           !srv 13.10.06
      real (kind=R8) ::
     *   t1, t2, t3, t11, t22, wrk0(nFc)
      character chns*3
      logical lnonzero
      save ncall
*   ..procedures
      intrinsic abs
      real (kind=R8) ::
     *   rscal
      external subini, subend, xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2xvfx, b2ux9p, b2uxm9                  !srv 02.01.07
!      rscal(t1) = t1
      rscal(t1) = 1.0e0_R8                                             !srv 13.10.06
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usco')
*   ..set internal parameters on first call                              !srv 25.08.99
      if (ncall.eq.0) then                                               !srv 06.07.06 {
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)
       call ipgeti ('b2usco_add_9', iadd_9)                              !srv 13.10.06
       call ipgeti ('b2usco_iout', iout)                                 !srv 13.10.06
      endif                                                              !srv 06.07.06 }
*   ..test nCv, nFc
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
*   ..test itcnt, rxg, pcm1
      call xertst (0.le.itcnt, 'faulty argument itcnt')
      call xertst (0.lt.rxg, 'faulty argument rxg')
      call xertst (0.le.pcm1, 'faulty argument pcm1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, pb
       call b2xvsg (nCv, nb, 1, 'nb', '.gt.')
       call b2xvsg (nCv, pb, 1, 'pb', '.gt.')
*    ..test sign of snbv
       call b2xvsg (nCv, snbv, 1, 'snbv', '.le.')
      endif

      aa = 0.0e0_R8

*   ..compute wrk0
      do iFc = 1, nFc
*!!     Maybe use volume weighting or inverse-volume weighting here.
cwdk need calccoeff here...
        wrk0(iFc) = 0.5_R8*geo%fcPbs(iFc)*
     &    (pccb1(mpg%fcCv(iFc,1))*nb(mpg%fcCv(iFc,1)) +
     &     pccb1(mpg%fcCv(iFc,2))*nb(mpg%fcCv(iFc,2)))
      enddo
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)

* ..compute the matrix
      do iCv = 1, nCv
        ireg  = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
          aa(mpg%cvNvP(iCv,1)) = (-snbv(iCv)*nb(iCv)/pb(iCv)+         !srv 13.10.06 {
     &     abs(rescb(iCv))/(rxg*pb(iCv)))*rscal(pb(iCv))        !srv 13.10.06 }
        else
          aa(mpg%cvNvP(iCv,1)) = rscal(pb(iCv))
        endif
      enddo

      do iFc = 1, nFc
        iCv1 = mpg%fcCv(iFc,1) ! cell number of c1
        iCv2 = mpg%fcCv(iFc,2) ! cell number of c2
*     ..contributions from flob(:,0), conb(:,0,:) ('5-point' contributions)
        t1 = - 0.5_R8*flob(iFc,0) - (conb(iFc,0,0)+pcm1*wrk0(iFc))
     &       - conb(iFc,0,2)                                       !srv 13.10.06
        t2 =   0.5_R8*flob(iFc,0) - (conb(iFc,0,0)+pcm1*wrk0(iFc))
     &       - conb(iFc,0,1)                                       !srv 13.10.06 
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
          t11 = t11 + 1/geo%cvVol(mpg%vxCv(mpg%vxCvP(iVx1,1)+inCv-1))
        enddo
        t22 = 0.0_R8
        do inCv = 1, mpg%vxCvP(iVx2,2)
          t22 = t22 + 1/geo%cvVol(mpg%vxCv(mpg%vxCvP(iVx2,1)+inCv-1))
        enddo
        ! contribution to equation of iCv1
        if (solvereg(mpg%cvReg(iCv1))) then
         ! contribution from iVx1
         do inCv = 1, mpg%vxCvP(iVx1,2)
          iCv = mpg%vxCv(mpg%vxCvP(iVx1,1) + inCv - 1)
          t3  = 1/(geo%cvVol(iCv)*t11)
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
         t3  = 1/(geo%cvVol(iCv)*t22)
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
         t3  = 1/(geo%cvVol(iCv)*t11)
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
         t3  = 1/(geo%cvVol(iCv)*t22)
         where (mpg%cvNv(mpg%cvNvP(iCv2,1):
     &                   mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) == iCv)
          aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) = 
     &     aa(mpg%cvNvP(iCv2,1):mpg%cvNvP(iCv2,1)+mpg%cvNvP(iCv2,2)-1) - 
     &     t3*t2*rscal(pb(iCv))
         end where
        enddo
       endif
      enddo

#ifdef WG_TODO
*   ..add matrix from contribution nine point stencil
      if(iadd_9 .eq. 1) then
       call b2uxm9 (nx, ny, nnreg(0), itcnt, solvereg,
     &   conb(-1,-1,0,0,3), aa, name)                             !srv 02.01.07
      endif
#endif
      do iCv = 1, nCv
        lnonzero=.false.
        do inCv = 1,mpg%cvNvP(iCv,2)
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
      call b2uxus (nCv, mpg, aa, itcnt, rescb, corpb, name)                 !srv 18.05.02
*   ..redo the scaling
      do iCv = 1, nCv
        corpb(iCv) = corpb(iCv)*rscal(pb(iCv))
      enddo
      if (iout.ne.0) then                                              !srv 17.06.02 {
        write (chns,'(i3.3)') isb
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-6),name//'_aam6'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-5),name//'_aam5'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-4),name//'_aam4'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-3),name//'_aam3'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-2),name//'_aam2'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1,-1),name//'_aam1'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 0),name//'_aacc'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 1),name//'_aap1'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 2),name//'_aap2'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 3),name//'_aap3'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 4),name//'_aap4'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 5),name//'_aap5'//chns)
!WG_TODO        call my_out_us(70,nx,ny,aa(-1,-1, 6),name//'_aap6'//chns)
        call my_out_us(70,nCv,0,corpb,name//'_corpb'//chns)
        call my_out_us(70,nCv,0,rescb,name//'_rescb'//chns)
        call my_out_us(70,nCv,0,nb,name//'_na'//chns)
      endif                                                            !srv 17.06.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usco

      end
