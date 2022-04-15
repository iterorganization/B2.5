*-----------------------------------------------------------------------
*.specification

      subroutine b2nxfc (nx, ny, isb, vol, pbs, rob, ub, flub, cvsb,
     &  flcb, cvcb, wrk0)
      use b2mod_types
      use b2mod_indirect
      use b2mod_boundary_namelist
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     *  pbs(-1:nx,-1:ny,0:1), rob(-1:nx,-1:ny), ub(-1:nx,-1:ny),
     &  vol(-1:nx,-1:ny)
*   ..inout arguments
      real (kind=R8) ::
     *  cvsb(-1:nx,-1:ny,0:1,0:1), flub(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  flcb(-1:nx,-1:ny,0:1,0:1), cvcb(-1:nx,-1:ny,0:1,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)                                          !srv 02.07.08

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXFC computes the transport terms in the parallel momentum
*     correction equation.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, n2, ix, iy, ncallout, ib
      integer, save :: iout = 0                                     !srv 17.06.02
      character charns*3                                            !srv 17.06.02
      !$OMP THREADPRIVATE(ncall, style, ncallout, iout)
      save ncall, style, ncallout
*   ..procedures
      external xertst, ipgeti, ipgetr, b2scopy
      external b2xvff, b2xvsg, hbout, intfaceh
*   ..initialisation
      data ncall/0/, style/1/, ncallout/-1/
*!!   (style is a numerical switch; needs experiments)

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     style - integer.
*     style specifies the precise form of interpolation used in the
*     computation of flcb and cvcb.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxfc')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2nxfc_style', style)
       call xertst (0.le.style.and.style.lt.2,
     &   'faulty internal parameter style')
       call ipgeti ('b2news_ncallout', ncallout)
       call ipgeti ('b2nxfc_iout', iout)                             !srv 17.06.02
       if ((iout.ne.0) .and. in_parallel()) then
        !$OMP MASTER
        write(*,*) 'B2NXFC OpenMP warning: no file output in ',
     &   'parallel mode'
        !$OMP END MASTER
        iout = 0
       end if
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*   ..test edge values of flub
       call b2xvff (nx, ny, flub, 'flub')
       if (ncall.eq.ncallout) then
        write(*,*) 'ncall',ncall
        call hbout('b2nxfc_xvsg ',nx,ny,2,rob,'rob',ub,'ub',cvsb,
     &  'cvsb',cvsb,'ua',cvsb,'na',cvsb,'ua')
       endif
*    ..test sign of rob
       call b2xvsg (n2, rob, 1, 'rob', '.gt.')
      endif

* ..compute convection factors
*   ..compute wrk0
      call intfaceh(nx,ny,leftix,leftiy,vol,ub,wrk0)
*   ..compute flcbx, cvcbx
      do iy = -1, ny
       if (style.eq.0) then
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          flcb(ix,iy,0,0) = 0.0_R8
          flcb(ix,iy,0,1) = 0.0_R8
          cvcb(ix,iy,0,0) = 0.0_R8
          cvcb(ix,iy,0,1) = 0.0_R8
         else
          flcb(ix,iy,0,0) = flub(ix,iy,0,0)+pbs(ix,iy,0)*(vol(ix,iy)*
     &      rob(leftix(ix,iy),leftiy(ix,iy))*
     &      ub(leftix(ix,iy),leftiy(ix,iy))+
     &      vol(leftix(ix,iy),leftiy(ix,iy))*rob(ix,iy)*ub(ix,iy))/
     &     (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
          flcb(ix,iy,0,1) = flub(ix,iy,0,1)
          cvcb(ix,iy,0,0) = cvsb(ix,iy,0,0)
          cvcb(ix,iy,0,1) = cvsb(ix,iy,0,1)
         endif
        enddo
       elseif(style.eq.1) then
        do ix = -1, nx                                                 !srv 26.11.02
         if(leftix(ix,iy).eq.-2) then                                  !srv 26.11.02 {
          flcb(ix,iy,0,0) = 0.0_R8
          flcb(ix,iy,0,1) = 0.0_R8
          cvcb(ix,iy,0,0) = 0.0_R8
          cvcb(ix,iy,0,1) = 0.0_R8
         else
          flcb(ix,iy,0,0) = flub(ix,iy,0,0)+pbs(ix,iy,0)*wrk0(ix,iy)*
     &      (rob(ix,iy)+rob(leftix(ix,iy),leftiy(ix,iy)))*0.5e0_R8     !srv 01.07.05 }
          flcb(ix,iy,0,1) = flub(ix,iy,0,1)
          cvcb(ix,iy,0,0) = cvsb(ix,iy,0,0)
          cvcb(ix,iy,0,1) = cvsb(ix,iy,0,1)
         endif                                                         !srv 26.11.02 }
        enddo
       else
        call xerrab ('Only style 0 or 1 in b2nxfc are included')
       endif
      enddo
*   ..compute flcby, cvcby
      call b2scopy (n2, flub(-1,-1,1,0), 1, flcb(-1,-1,1,0), 1)
      call b2scopy (n2, flub(-1,-1,1,1), 1, flcb(-1,-1,1,1), 1)
      call b2scopy (n2, cvsb(-1,-1,1,0), 1, cvcb(-1,-1,1,0), 1)
      call b2scopy (n2, cvsb(-1,-1,1,1), 1, cvcb(-1,-1,1,1), 1)

csw 18feb2014 remove fluxes between guard cells
        do ix=-1,nx
         do iy=-1,ny
           if(iffboundx(ix,iy).ne.0) then
             flub(ix,iy,0,0:1)=0.0_R8
             cvsb(ix,iy,0,0:1)=0.0_R8
             flcb(ix,iy,0,0:1)=0.0_R8
             cvcb(ix,iy,0,0:1)=0.0_R8
           endif
           if(iffboundy(ix,iy).ne.0) then
             flub(ix,iy,1,0:1)=0.0_R8
             cvsb(ix,iy,1,0:1)=0.0_R8
             flcb(ix,iy,1,0:1)=0.0_R8
             cvcb(ix,iy,1,0:1)=0.0_R8
           endif
         enddo
        enddo

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            flcb(ix,iy,:,:)=0.0_R8
            cvcb(ix,iy,:,:)=0.0_R8
            flub(ix,iy,:,:)=0.0_R8
            cvsb(ix,iy,:,:)=0.0_R8
          end do
        end do
      end do

      if (iout.ne.0) then                                              !sv 17.06.02 {
        write (charns,'(i3.3)') isb
        call my_out(70,nx,ny,flub(-1,-1,0,0),'b2nxfc_flubx'//charns)
        call my_out(70,nx,ny,flub(-1,-1,0,1),'b2nxfc_flubxy'//charns)
        call my_out(70,nx,ny,flub(-1,-1,1,0),'b2nxfc_flubyx'//charns)
        call my_out(70,nx,ny,flub(-1,-1,1,1),'b2nxfc_fluby'//charns)
        call my_out(70,nx,ny,flcb(-1,-1,0,0),'b2nxfc_flcbx'//charns)
        call my_out(70,nx,ny,flcb(-1,-1,0,1),'b2nxfc_flcbxy'//charns)
        call my_out(70,nx,ny,flcb(-1,-1,1,0),'b2nxfc_flcbyx'//charns)
        call my_out(70,nx,ny,flcb(-1,-1,1,1),'b2nxfc_flcby'//charns)
        call my_out(70,nx,ny,cvcb(-1,-1,0,0),'b2nxfc_cvcbx'//charns)
        call my_out(70,nx,ny,cvcb(-1,-1,0,1),'b2nxfc_cvcbxy'//charns)
        call my_out(70,nx,ny,cvcb(-1,-1,1,0),'b2nxfc_cvcbyx'//charns)
        call my_out(70,nx,ny,cvcb(-1,-1,1,1),'b2nxfc_cvcby'//charns)
        call my_out(70,nx,ny,rob,'b2nxfc_rob'//charns)
      endif                                                            !sv 17.06.02 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxfc

      end subroutine b2nxfc
