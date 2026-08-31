*-----------------------------------------------------------------------
*.specification

      subroutine b2upco_3D (nCv, nFc, nVx, geo, mpg, nregionv,
     &   um_solve, solvereg, rxf, pcm0, corpb, pccb0_pol, pccb0_rad, 
     &   pb, nb, ub_pol, ub_rad, psnl, iout)
      use b2mod_types
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, nregionv, iout
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2PlasmaSnapshot), intent (inout) :: psnl
      real (kind=R8) ::
     *  rxf, pcm0, corpb(nCv), pccb0_pol(nCv), pccb0_rad(nCv)
      logical um_solve, solvereg(0:nregionv)
*   ..input/output arguments
      real (kind=R8) ::
     *  pb(nCv), nb(nCv), ub_pol(nCv), ub_rad(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2UPCO_3D corrects the state of the fluid molecules towards the 
*     solution of continuity equation when using a 3D momentum equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer i, iCv, iFt, ireg
*     ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) :: wrk(nCv,0:1), wrkv(nVx,0:1)
      real (kind=R8) :: t0, pb_(nCv), nb_(nCv), corr_u_pol(nCv),
     .                  corr_u_rad(nCv)                                                          
*   ..procedures
      external xertst, ipgetr
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2upco_3D')
*   ..test nCv, nFc
      call xertst (0.lt.nCv.and.0.lt.nFc, 'faulty argument nCv, nFc')
*   ..test rxf, pcm0
      call xertst (0.0_R8.le.rxf.and.rxf.le.1.0_R8,
     .  'faulty argument rxf')
      call xertst (0.0_R8.le.pcm0, 'faulty argument pcm0')

* ..implement the correction
*   ..compute wrk
      call gradc(nCv, nFc, nVx, 0, geo, mpg, corpb, wrkv, wrk)
*   ..compute new pb and nb
      do iCv = 1, nCv
        if (0.0_R8.le.corpb(iCv)) then
         t0 = 1.0_R8+rxf*corpb(iCv)/pb(iCv)
        else
         t0 = 1.0_R8/(1.0_R8-rxf*corpb(iCv)/pb(iCv))
        endif
        nb_(iCv) = t0*nb(iCv)                        
        pb_(iCv) = t0*pb(iCv)                       
      enddo
      nb = nb_
      pb = pb_
*   ..compute new ub
      if (um_solve) then
       do iCv = 1, mpg%nCi ! internal cells only
*       (note the further damping on the velocity correction)
         ireg = mpg%cvReg(iCv)
         if(solvereg(ireg)) then
           t0 = pcm0*rxf*pccb0_pol(iCv)*wrk(iCv,0)
           corr_u_pol(iCv) = t0
           ub_pol(iCv) = ub_pol(iCv)-t0
           t0 = pcm0*rxf*pccb0_rad(iCv)*wrk(iCv,1)
           corr_u_rad(iCv) = t0
           ub_rad(iCv) = ub_rad(iCv)-t0
         endif
       enddo
      endif

      if (iout.ne.0) then 
        call my_out_us(70,nCv,0,corr_u_pol,'b2upco_corr_u_pol')
        call my_out_us(70,nCv,0,corr_u_rad,'b2upco_corr_u_rad')
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2upco

      end subroutine b2upco_3D
