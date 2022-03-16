*-----------------------------------------------------------------------
*.specification

      subroutine b2usm9 (nx, ny, isb, nregionv, solvereg, itcnt,     !srv 04.07.08
     &  rxg, vol, hx, hy, hz, qz, qs, bb, rob, rzb, pb, pz, ub,      !srv 09.01.01
     &  smb, flcb, cvsb, resmb, ctcfb, corub, pccb, aa9, aad0, name) !srv 06.04.07
      use b2mod_indirect
      use b2mod_types
      use b2mod_boundary_namelist
     , , only : iffboundx, iffboundy
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, nregionv, itcnt
      character*(*) name                                             !srv 18.05.02
      real (kind=R8) ::
     *  rxg, bb(-1:nx,-1:ny,0:3), 
     &  rob(-1:nx,-1:ny), rzb(-1:nx,-1:ny),
     &  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),          !srv 09.01.01
     &  hz(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), qs(-1:nx,-1:ny,0:1),
     &  pb(-1:nx,-1:ny), pz(-1:nx,-1:ny),
     &  ub(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3),
     &  flcb(-1:nx,-1:ny,0:1,0:1), cvsb(-1:nx,-1:ny,0:1,0:1),
     &  resmb(-1:nx,-1:ny), ctcfb(-1:nx,-1:ny)                       !srv 06.04.07
      logical solvereg(0:nregionv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  corub(-1:nx,-1:ny), pccb(-1:nx,-1:ny,0:1,0:1)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  aa9(-1:nx,-1:ny,-6:6), aad0(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2usm9 computes a correction towards the solution of the
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
      integer ncall, n2, ix, iy, ireg, ileft, ibot, itop
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !srv 06.07.06
      integer, save :: iout = 0                                       !srv 17.06.02
      real (kind=R8) ::
     *  cfc0, t0, t1, t2, t3, t4, flo, dif, csb, csb1
      save ncall, cfc0
*   ..procedures
      intrinsic max, abs, sqrt
      real (kind=R8) ::
     *  hybr, upwind
      external xertst, ipgetr, sfill
      external b2xvsg, b2xvff, b2ux5p, intfaceh, intfacev
      hybr(flo,dif) = max(dif,abs(flo/2.0_R8),
     .                    sqrt(dif**2+flo**2/4.0_R8))
      upwind(flo,t0,t1) = max(flo,0.0_R8)*t0+min(flo,0.0_R8)*t1
*   ..initialisation
      data ncall/0/, cfc0/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2usm9')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2usm9_cfc0', cfc0)
       call xertst (0.0_R8.le.cfc0, 'faulty internal parameter cfc0')
       call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)               !srv 06.07.06
       call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)               !srv 06.07.06
       call ipgeti ('b2usm9_iout', iout)                              !srv 17.06.02
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
       call b2xvsg (n2, rzb, 1, 'rzb', '.gt.')
       call b2xvsg (n2, pb, 1, 'pb', '.gt.')
       call b2xvsg (n2, pz, 1, 'pz', '.gt.')
       call b2xvsg (n2, smb(-1,-1,1), 1, 'smb1', '.le.')
       call b2xvsg (n2, smb(-1,-1,3), 1, 'smb3', '.le.')
       call b2xvsg (n2, cvsb(-1,-1,0,0), 1, 'cvsb', '.ge.')
       call b2xvsg (n2, cvsb(-1,-1,1,1), 1, 'cvsb', '.ge.')
*    ..test edge values of flcb, cvsb
       call b2xvff (nx, ny, flcb, 'flcb')
       call b2xvff (nx, ny, cvsb, 'cvsb')
      endif

* ..compute the correction
      aa9 = 0.0_R8
*   ..compute the nine-point matrix
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         csb = sqrt(pz(ix,iy)/rzb(ix,iy))
         csb1 = sqrt(ub(ix,iy)**2+csb**2)
         aad0(ix,iy) = -smb(ix,iy,1)-
     &    smb(ix,iy,3)*rob(ix,iy)+
     &    cfc0*abs(vol(ix,iy)*hz(ix,iy)/hx(ix,iy)*bb(ix,iy,0)/       !srv 09.01.01
     &    bb(ix,iy,3)*csb)*rob(ix,iy)
         aa9(ix,iy,0) = aad0(ix,iy)+abs(resmb(ix,iy))/(rxg*csb1)
        else
         aa9(ix,iy,0) = 1.0_R8
         aad0(ix,iy) = 0.0_R8
        endif
       enddo
      enddo

      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
!        if(iffboundx(ix,iy).ne.0) then
!          aa9(ix,iy,-1) = 0.0_R8
!         aa9(ix,iy,1) = 0.0_R8
!          cycle
!        endif

        ireg = region(ix,iy,0)
        if(leftix(ix,iy).eq.-2) then
         aa9(ix,iy,-1) = 0.0_R8
        else
         if(no_Ptncr_x_co.eq.0) then                                  !sv 06.07.06 {
          t0 =  0.5_R8*flcb(ix,iy,0,0)-
     .          hybr(flcb(ix,iy,0,0),cvsb(ix,iy,0,0))
          t1 = -0.5_R8*flcb(ix,iy,0,0)-
     .          hybr(flcb(ix,iy,0,0),cvsb(ix,iy,0,0))
         else
          t0 =  flcb(ix,iy,0,0)/2.0_R8-cvsb(ix,iy,0,0)
          t1 = -flcb(ix,iy,0,0)/2.0_R8-cvsb(ix,iy,0,0)
         endif                                                        !sv 06.07.06 }
cwdk compute contributions from 9-point stencil
         t2 =  flcb(ix,iy,0,1)
         t3 =  sign(hybr(t2,cvsb(ix,iy,0,1)) - 
     &         0.5_R8*abs(t2),cvsb(ix,iy,0,1))
         t4 = -(hy(ix,iy) + hy(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy) + hx(leftix(ix,iy),leftiy(ix,iy)))*
     &         qs(ix,iy,0)
         t0 = t0 + upwind(t2,0.0_R8,1.0_R8) - t3*t4
         t1 = t1 - upwind(t2,1.0_R8,0.0_R8) - t3*t4

         ileft = region(leftix(ix,iy),leftiy(ix,iy),0)
         if (solvereg(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),1) = 
     >      aa9(leftix(ix,iy),leftiy(ix,iy),1) + t0
         if (solvereg(ireg))
     >    aa9(ix,iy,0) = aa9(ix,iy,0)-t0
         if (solvereg(ileft))
     >    aa9(leftix(ix,iy),leftiy(ix,iy),0) =
     &    aa9(leftix(ix,iy),leftiy(ix,iy),0)-t1
         if (solvereg(ireg))
     >    aa9(ix,iy,-1) = aa9(ix,iy,-1) + t1
cwdk add 9-point contributions to matrix
         if(bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1) then
          t0  = -t3/
     &       (1/vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &              leftiy(topix(ix,iy),topiy(ix,iy))) +
     &        1/vol(topix(ix,iy),topiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          t1  = -t3/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
          if (solvereg(ireg)) then
           aa9(ix,iy,-4) = aa9(ix,iy,-4)           + t1/
     &       vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &           bottomiy(leftix(ix,iy),leftiy(ix,iy)))
           aa9(ix,iy,-3) = aa9(ix,iy,-3)           + t1/
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
           aa9(ix,iy,-1) = aa9(ix,iy,-1) - (t0 - t1)/
     &       vol(leftix(ix,iy),leftiy(ix,iy))
           aa9(ix,iy, 0) = aa9(ix,iy, 0) - (t0 - t1)/
     &       vol(ix,iy)
           aa9(ix,iy, 2) = aa9(ix,iy, 2) -  t0/
     &       vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &           leftiy(topix(ix,iy),topiy(ix,iy)))
           aa9(ix,iy, 3) = aa9(ix,iy, 3) -  t0/
     &       vol(topix(ix,iy),topiy(ix,iy))
          endif
          if (solvereg(ileft)) then
           aa9(leftix(ix,iy),leftiy(ix,iy),-3) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy),-3)       - t1/
     &       vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &           bottomiy(leftix(ix,iy),leftiy(ix,iy)))
           aa9(leftix(ix,iy),leftiy(ix,iy),-2) =
     &       aa9(leftix(ix,iy),leftiy(ix,iy),-2)       - t1/
     &       vol(bottomix(ix,iy),bottomiy(ix,iy))
           aa9(leftix(ix,iy),leftiy(ix,iy), 0) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 0) + (t0 - t1)/
     &       vol(leftix(ix,iy),leftiy(ix,iy))
           aa9(leftix(ix,iy),leftiy(ix,iy), 1) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 1) + (t0 - t1)/
     &       vol(ix,iy)
           aa9(leftix(ix,iy),leftiy(ix,iy), 3) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 3) + t0/
     &       vol(leftix(topix(ix,iy),topiy(ix,iy)),
     &           leftiy(topix(ix,iy),topiy(ix,iy)))
           aa9(leftix(ix,iy),leftiy(ix,iy), 4) = 
     &       aa9(leftix(ix,iy),leftiy(ix,iy), 4) + t0 /
     &       vol(topix(ix,iy),topiy(ix,iy))
          endif
         endif
        endif
        if(rightix(ix,iy).eq.nx+1 .or. .not.solvereg(ireg))
     >   aa9(ix,iy,1) = 0.0_R8
       enddo
      enddo

      do iy = -1, ny
       do ix = -1, nx
csw 18feb2014 corrected BCs (SOLPS4)
!        if(iffboundy(ix,iy).ne.0) then
!          aa9(ix,iy,-3) = 0.0_R8
!          aa9(ix,iy,3) = 0.0_R8
!          cycle
!        endif

        ireg = region(ix,iy,0)
        if(bottomiy(ix,iy).eq.-2) then
         aa9(ix,iy,-3) = 0.0_R8
        else
         if(no_Ptncr_y_co.eq.0) then                                 !srv 06.07.06 {
          t0 =  0.5_R8*flcb(ix,iy,1,1) - 
     &            hybr(flcb(ix,iy,1,1),cvsb(ix,iy,1,1))
          t1 = -0.5_R8*flcb(ix,iy,1,1)- 
     &            hybr(flcb(ix,iy,1,1),cvsb(ix,iy,1,1))
          t2 = flcb(ix,iy,1,1)
          t3 = hybr(flcb(ix,iy,1,1),cvsb(ix,iy,1,1))
         else
          t0 =  0.5_R8*flcb(ix,iy,1,1) - cvsb(ix,iy,1,1)
          t1 = -0.5_R8*flcb(ix,iy,1,1) - cvsb(ix,iy,1,1)
          t2 =  0.0_R8
          t3 =  cvsb(ix,iy,1,1)
         endif                                                        !srv 06.07.06 }

          ibot = region(bottomix(ix,iy),bottomiy(ix,iy),0)
          if (solvereg(ibot))
     >     aa9(bottomix(ix,iy),bottomiy(ix,iy),3) = 
     >      aa9(bottomix(ix,iy),bottomiy(ix,iy),3) + t0
          if (solvereg(ireg))
     >     aa9(ix,iy,0) = aa9(ix,iy,0)-t0
          if (solvereg(ibot))
     >     aa9(bottomix(ix,iy),bottomiy(ix,iy),0) =
     &              aa9(bottomix(ix,iy),bottomiy(ix,iy),0)-t1
          if (solvereg(ireg))
     >     aa9(ix,iy,-3) = aa9(ix,iy,-3) + t1
cwdk add contributions from 9-point stencil
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
           t4 = -(hy(ix,iy) + hy(bottomix(ix,iy),bottomiy(ix,iy)))/
     &           (hx(ix,iy) + hx(bottomix(ix,iy),bottomiy(ix,iy)))*
     &           qs(ix,iy,1)
           t0 = -(t3 - 0.5_R8*abs(t2))*t4/
     &       (1/vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(rightix(ix,iy),rightiy(ix,iy)))
           t1 = -(t3 - 0.5_R8*abs(t2))*t4/
     &       (1/vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy))) +
     &        1/vol(bottomix(ix,iy),bottomiy(ix,iy))+
     &        1/vol(ix,iy)+
     &        1/vol(leftix(ix,iy),leftiy(ix,iy)))
           if (solvereg(ireg)) then
            aa9(ix,iy,-4) = aa9(ix,iy,-4)       + t1/
     &        vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))
            aa9(ix,iy,-3) = aa9(ix,iy,-3) - (t0 - t1)/
     &        vol(bottomix(ix,iy),bottomiy(ix,iy))
            aa9(ix,iy,-2) = aa9(ix,iy,-2) -  t0/
     &        vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))
            aa9(ix,iy,-1) = aa9(ix,iy,-1)       + t1/
     &        vol(leftix(ix,iy),leftiy(ix,iy))
            aa9(ix,iy, 0) = aa9(ix,iy, 0) - (t0 - t1)/
     &        vol(ix,iy)
            aa9(ix,iy, 1) = aa9(ix,iy, 1) -  t0/
     &        vol(rightix(ix,iy),rightiy(ix,iy))
           endif
           if (solvereg(ibot)) then
            aa9(bottomix(ix,iy),bottomiy(ix,iy),-1) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy),-1)       - t1/
     &        vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &            bottomiy(leftix(ix,iy),leftiy(ix,iy)))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 0) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 0) + (t0 - t1)/
     &        vol(bottomix(ix,iy),bottomiy(ix,iy)) 
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 1) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 1) +  t0/
     &        vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &            bottomiy(rightix(ix,iy),rightiy(ix,iy)))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 2) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 2)      - t1/
     &        vol(leftix(ix,iy),leftiy(ix,iy))
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 3) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 3) + (t0 - t1)/
     &        vol(ix,iy)
            aa9(bottomix(ix,iy),bottomiy(ix,iy), 4) = 
     &        aa9(bottomix(ix,iy),bottomiy(ix,iy), 4) +  t0/
     &        vol(rightix(ix,iy),rightiy(ix,iy))
           endif
          endif
cxpb
          if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).eq.iy .and.
     &       topix(topix(ix,iy),min(topiy(ix,iy),ny)).eq.ix .and.
     &       iy.ne.ny .and. aa9(ix,iy,3).eq.0.0_R8) then
           if(ncall.lt.3) write(*,*) 'Executing xpb:b2usm9 for ',ix,iy
           t0 =  flcb(topix(ix,iy),topiy(ix,iy),1,1)/2.0_R8-
     &           hybr(flcb(topix(ix,iy),topiy(ix,iy),1,1),
     &                cvsb(topix(ix,iy),topiy(ix,iy),1,1))
            t1 = -flcb(topix(ix,iy),topiy(ix,iy),1,1)/2.0_R8-
     &            hybr(flcb(topix(ix,iy),topiy(ix,iy),1,1),
     &                 cvsb(topix(ix,iy),topiy(ix,iy),1,1))
            itop = region(topix(ix,iy),topiy(ix,iy),0)
            if (solvereg(ireg))
     >       aa9(ix,iy,3) = t0
            if (solvereg(itop))
     >       aa9(topix(ix,iy),topiy(ix,iy),0) = 
     &       aa9(topix(ix,iy),topiy(ix,iy),0)-t0
            if (solvereg(ireg))
     >       aa9(ix,iy,0) = aa9(ix,iy,0)-t1
            if (solvereg(itop))
     >       aa9(topix(ix,iy),topiy(ix,iy),3) = t1
          endif
cxpb
        endif
        if(topiy(ix,iy).eq.ny+1 .or. .not.solvereg(ireg))
     >   aa9(ix,iy,3) = 0.0_R8
       enddo
      enddo
*   ..compute pccb
      pccb = 0.0_R8
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         pccb(ix,iy,0,0) = (vol(ix,iy)*hz(ix,iy)/hx(ix,iy))*            !srv 09.01.01
     .    bb(ix,iy,0)/bb(ix,iy,3)/aa9(ix,iy,0)                        !srv 11.01.13
         pccb(ix,iy,1,1) = (vol(ix,iy)*hz(ix,iy)/hx(ix,iy))*            !srv 09.01.01
     .    bb(ix,iy,0)/bb(ix,iy,3)/aad0(ix,iy)                         !srv 11.01.13
        else
         resmb(ix,iy) = 0.0_R8   ! take opportunity to correct RHS
        endif
       enddo
      enddo
*   ..contribution from linearization of centrifugal force
      aa9(:,:,0) = aa9(:,:,0) - min(0.0e0_R8,ctcfb)                   !srv 06.04.07
      do iy=-1,ny
       do ix=-1,nx
        if (aa9(ix,iy,0).eq.0.0_R8 .and.
     &      aa9(ix,iy,1).eq.0.0_R8 .and. aa9(ix,iy,-1).eq.0.0_R8 .and.
     &      aa9(ix,iy,2).eq.0.0_R8 .and. aa9(ix,iy,-2).eq.0.0_R8 .and.
     &      aa9(ix,iy,3).eq.0.0_R8 .and. aa9(ix,iy,-3).eq.0.0_R8 .and.
     &      aa9(ix,iy,4).eq.0.0_R8 .and. aa9(ix,iy,-4).eq.0.0_R8 .and.
     &      aa9(ix,iy,5).eq.0.0_R8 .and. aa9(ix,iy,-5).eq.0.0_R8 .and.
     &      aa9(ix,iy,6).eq.0.0_R8 .and. aa9(ix,iy,-6).eq.0.0_R8) then
          aa9(ix,iy,0) = 1.0_R8
          resmb(ix,iy) = 0.0_R8
        end if
       enddo
      enddo
*   ..solve the correction equation
      write(*,'(a,a)') 'Calling b2ux9p from ', name
      call b2ux9p (nx, ny, aa9, itcnt, resmb, corub, name)            !srv 18.03.10

      if (iout.ne.0) then                                             !srv 17.06.02 {
       call my_out(70,nx,ny,aa9(-1,-1,-4),name//'_aa9m4')
       call my_out(70,nx,ny,aa9(-1,-1,-3),name//'_aa9m3')
       call my_out(70,nx,ny,aa9(-1,-1,-2),name//'_aa9m2')
       call my_out(70,nx,ny,aa9(-1,-1,-1),name//'_aa9m1')
       call my_out(70,nx,ny,aa9(-1,-1, 0),name//'_aa9m0')
       call my_out(70,nx,ny,aa9(-1,-1, 1),name//'_aa9p1')
       call my_out(70,nx,ny,aa9(-1,-1, 2),name//'_aa9p2')
       call my_out(70,nx,ny,aa9(-1,-1, 3),name//'_aa9p3')
       call my_out(70,nx,ny,aa9(-1,-1, 4),name//'_aa9p4')
       call my_out(70,nx,ny,resmb,name//'_resmb')
       call my_out(70,nx,ny,corub,name//'_corub')
       call my_out(70,nx,ny,pccb(-1,-1,0,0),name//'_pccbx')
       call my_out(70,nx,ny,pccb(-1,-1,1,1),name//'_pccby')
      endif                                                           !srv 17.06.02 }
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2usm9

      end
