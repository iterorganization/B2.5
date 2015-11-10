*-----------------------------------------------------------------------
*.specification

      subroutine b2usmo (nx, ny, isb, nregionv, solvereg, itcnt,     !srv 04.07.08
     &  rxg, vol, hx, hz, bb, rob, pb, ub,                           !srv 09.01.01
     &  smb, flcb, cvsb, resmb, ctcfb, corub, pccb, aa5, aad0, name) !srv 06.04.07
      use b2mod_indirect
      use b2mod_types
      use b2mod_boundary_namelist
     , , only : iffboundx, iffboundy
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, nregionv, itcnt
      character*(*) name                                             !srv 18.05.02
      real (kind=R8) ::
     *  rxg, bb(-1:nx,-1:ny,0:3), rob(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hz(-1:nx,-1:ny),          !srv 09.01.01
     &  pb(-1:nx,-1:ny), ub(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3),
     &  flcb(-1:nx,-1:ny,0:1), cvsb(-1:nx,-1:ny,0:1),
     &  resmb(-1:nx,-1:ny), ctcfb(-1:nx,-1:ny)                       !srv 06.04.07
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corub(-1:nx,-1:ny), pccb(-1:nx,-1:ny,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa5(-1:nx,-1:ny,-2:2), aad0(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2USMO computes a correction towards the solution of the
*     parallel momentum balance equation for one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the velocity correction, corub, for one
*     species.
*     ...
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop, k
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !sv 06.07.06
      integer, save :: iout = 0                                       !srv 17.06.02
      character chns*3                                                !srv 17.06.02
      real (kind=R8) ::
     *  cfc0, t0, t1, flo, dif, csb, csb1
      save ncall, cfc0
*   ..procedures
      intrinsic max, abs, sqrt
      real (kind=R8) ::
     *  hybr
      external subini, subend, xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2ux5p
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
*   ..initialisation
      data ncall/0/, cfc0/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usmo')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2usmo_cfc0', cfc0)
       call xertst (0.0_R8.le.cfc0, 'faulty internal parameter cfc0')
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)               !sv 06.07.06
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)               !sv 06.07.06
       call ipgeti ('b2usmo_iout', iout)                             !srv 17.06.02
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test itcnt
      call xertst (0.le.itcnt, 'faulty argument itcnt')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of rob, pb, smb1, smb3, cvsb
       call b2xvsg (n2, rob, 1, 'rob', '.gt.')
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
       call b2xvsg (n2, smb(-1,-1,1), 1, 'smb1', '.le.')
       call b2xvsg (n2, smb(-1,-1,3), 1, 'smb3', '.le.')
       call b2xvsg (n2*2, cvsb, 1, 'cvsb', '.ge.')
*    ..test edge values of flcb, cvsb
       call b2xvff (nx, ny, flcb, 'flcb')
       call b2xvff (nx, ny, cvsb, 'cvsb')
      endif

* ..compute the correction
      call sfill (n2*5, 0.0_R8, aa5, 1)
*   ..compute the five-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         csb = sqrt(pb(ix,iy)/rob(ix,iy))
         csb1 = sqrt(ub(ix,iy)**2+csb**2)
         aad0(ix,iy) = -smb(ix,iy,1)-
     &    smb(ix,iy,3)*rob(ix,iy)+
     &    cfc0*abs(vol(ix,iy)*hz(ix,iy)/hx(ix,iy)*bb(ix,iy,0)/       !srv 09.01.01
     &    bb(ix,iy,3)*csb)*rob(ix,iy)
         aa5(ix,iy,0) = aad0(ix,iy)+abs(resmb(ix,iy))/(rxg*csb1)
        else
         aa5(ix,iy,0) = 1.0_R8
         aad0(ix,iy) = 0.0_R8
        endif
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundx(ix,iy).ne.0) then
          aa5(ix,iy,-1) = 0.0_R8
          aa5(ix,iy,1) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa5(ix,iy,-1) = 0.0_R8
        else
         if(no_Ptncr_x_co.eq.0) then                                  !sv 06.07.06 {
          t0 =  flcb(ix,iy,0)/2-hybr(flcb(ix,iy,0),cvsb(ix,iy,0))
          t1 = -flcb(ix,iy,0)/2-hybr(flcb(ix,iy,0),cvsb(ix,iy,0))
         else
          t0 =  flcb(ix,iy,0)/2-cvsb(ix,iy,0)
          t1 = -flcb(ix,iy,0)/2-cvsb(ix,iy,0)
         endif                                                        !sv 06.07.06 }
         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solvereg(ileft)) 
     >    aa5(leftix(ix,iy),leftiy(ix,iy),1) = t0
         if (solvereg(ireg))
     >    aa5(ix,iy,0) = aa5(ix,iy,0)-t0
         if (solvereg(ileft))
     >    aa5(leftix(ix,iy),leftiy(ix,iy),0) = 
     &    aa5(leftix(ix,iy),leftiy(ix,iy),0)-t1
         if (solvereg(ireg))
     >    aa5(ix,iy,-1) = t1
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg)) 
     >   aa5(ix,iy,1) = 0.0_R8
       enddo
      enddo
      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
        if(iffboundy(ix,iy).ne.0) then
          aa5(ix,iy,-2) = 0.0_R8
          aa5(ix,iy,2) = 0.0_R8
          cycle
        endif

        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
          aa5(ix,iy,-2) = 0.0_R8
        else
         if(no_Ptncr_y_co.eq.0) then                                  !srv 06.07.06 {
          t0 =  flcb(ix,iy,1)/2 - hybr(flcb(ix,iy,1),cvsb(ix,iy,1))
          t1 = -flcb(ix,iy,1)/2 - hybr(flcb(ix,iy,1),cvsb(ix,iy,1))
         else
          t0 =  flcb(ix,iy,1)/2 - cvsb(ix,iy,1)
          t1 = -flcb(ix,iy,1)/2 - cvsb(ix,iy,1)
         endif                                                        !srv 06.07.06 }
          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy),bottomiy(ix,iy),2) = t0
          if (solvereg(ireg))
     >     aa5(ix,iy,0) = aa5(ix,iy,0)-t0
          if (solvereg(ibot))
     >     aa5(bottomix(ix,iy),bottomiy(ix,iy),0) = 
     &              aa5(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
          if (solvereg(ireg))
     >     aa5(ix,iy,-2) = t1
cxpb
          if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &       topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &	     iy.ne.ny .and. aa5(ix,iy,2).eq.0.0_R8) then
           if(ncall.lt.3) write(*,*) 'Executing xpb:b2usmo for ',ix,iy
            t0 = flcb(topix(ix,iy),topiy(ix,iy),1)/2-
     &           hybr(flcb(topix(ix,iy),topiy(ix,iy),1),
     &	              cvsb(topix(ix,iy),topiy(ix,iy),1))
            t1 = -flcb(topix(ix,iy),topiy(ix,iy),1)/2-
     &            hybr(flcb(topix(ix,iy),topiy(ix,iy),1),
     &                 cvsb(topix(ix,iy),topiy(ix,iy),1))
            itop = region(topix(ix,iy),topiy(ix,iy),0)
            if (solvereg(ireg))
     >       aa5(ix,iy,2) = t0
            if (solvereg(itop))
     >       aa5(topix(ix,iy),topiy(ix,iy),0) = 
     &	     aa5(topix(ix,iy),topiy(ix,iy),0)-t0
            if (solvereg(ireg))
     >       aa5(ix,iy,0) = aa5(ix,iy,0)-t1
            if (solvereg(itop))
     >       aa5(topix(ix,iy),topiy(ix,iy),2) = t1
          endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solvereg(ireg)) 
     >   aa5(ix,iy,2) = 0.0_R8
       enddo
      enddo
*   ..compute pccb
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         pccb(ix,iy,0) = (vol(ix,iy)*hz(ix,iy)/hx(ix,iy))*            !srv 09.01.01
     .    bb(ix,iy,0)/bb(ix,iy,3)/aa5(ix,iy,0)                        !srv 11.01.13
         pccb(ix,iy,1) = (vol(ix,iy)*hz(ix,iy)/hx(ix,iy))*            !srv 09.01.01
     .    bb(ix,iy,0)/bb(ix,iy,3)/aad0(ix,iy)                         !srv 11.01.13
        else
         pccb(ix,iy,0) = 0.0_R8
         pccb(ix,iy,1) = 0.0_R8   
         resmb(ix,iy) = 0.0_R8   ! take opportunity to correct RHS
        endif
       enddo
      enddo
*   ..contribution from linearization of centrifugal force
      aa5(:,:,0) = aa5(:,:,0) - min(0.0e0_R8,ctcfb)                   !srv 06.04.07
      do iy=-1,ny
       do ix=-1,nx
        if (aa5(ix,iy,0).eq.0.0_R8 .and.
     &      aa5(ix,iy,1).eq.0.0_R8 .and. aa5(ix,iy,-1).eq.0.0_R8 .and.
     &      aa5(ix,iy,2).eq.0.0_R8 .and. aa5(ix,iy,-2).eq.0.0_R8) then
          aa5(ix,iy,0) = 1.0_R8
          resmb(ix,iy) = 0.0_R8
        end if
       enddo
      enddo
*   ..solve the correction equation
      write(*,'(a,a)') 'Calling b2ux5p from ', name
      call b2ux5p (nx, ny, aa5, itcnt, resmb, corub, name)            !srv 18.03.10

      if (iout.ne.0) then                                             !srv 17.06.02 {
       write (chns,'(i3.3)') isb    
       call my_out(70,nx,ny,aa5(-1,-1,-2),name//'_aa5m2'//chns)
       call my_out(70,nx,ny,aa5(-1,-1,-1),name//'_aa5m1'//chns)
       call my_out(70,nx,ny,aa5(-1,-1, 0),name//'_aa5m0'//chns)
       call my_out(70,nx,ny,aa5(-1,-1, 1),name//'_aa5p1'//chns)
       call my_out(70,nx,ny,aa5(-1,-1, 2),name//'_aa5p2'//chns)
       call my_out(70,nx,ny,resmb,name//'_resmb'//chns)
       call my_out(70,nx,ny,corub,name//'_corub'//chns)
       call my_out(70,nx,ny,pccb(-1,-1,0),name//'_pccbx'//chns)
       call my_out(70,nx,ny,pccb(-1,-1,1),name//'_pccby'//chns)
      endif                                                           !srv 17.06.02 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usmo

      end
