*-----------------------------------------------------------------------
*.specification

      subroutine b2txnpsch (nx, ny, isb, guard_flows, facdrift,    !xpb !iyv 06.02.14 02.04.18 {
     &        kmprti, cddb, conb)                                       !srv 09.01.01 13.10.06
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, guard_flows                                  !xpb !srv 08.09.01 !lk 02.06.11
      real (kind=R8) ::
     &  kmprti(-1:nx,-1:ny,0:1), facdrift(-1:nx,-1:ny),
     &  cddb(-1:nx,-1:ny,0:1)                                           !srv 13.10.06 !iyv 06.02.14
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  conb(-1:nx,-1:ny,0:1,0:4) 
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TXNPSCH computes the coefficients of the Pfirsch-Schlueter particle
*     fluxes approximation
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the coefficients conb that linearise the
*     response of the Pfirsch-Schlueter particle flux.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, ix, iy, area_fix
      integer, save :: iout = 0                                         !srv 02.04.18
      integer, save :: PSch_style                                       !srv 12.11.08
      character chns*3                                                  !srv 16.09.02
      real (kind=R8), save :: PSch = 1.0_R8                             !srv 17.01.07
      save ncall, area_fix
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr                   !xpb
      external b2xvsg, b2xvff, intfaceh, intfacev
*   ..initialisation
      data ncall/0/, area_fix/3/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tnpsch')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tnpsch_iout', iout)                             !srv 02.04.18
        call ipgetr ('b2tfnb_PSch', PSch)                               !srv 17.01.07
        if (redef_gmtry.eq.0) then
         PSch_style = 1
        else
         PSch_style = 2
        end if
        call ipgeti ('b2tfnb_PSch_style', PSch_style)                   !srv 12.11.08
      endif

      write(chns,'(i3.3)') isb
      conb  = 0.0e0_R8                                                  !srv 13.10.06

*   ..compute P.Sch x-flux                                              !srv 13.10.06 {
      if(PSch_style.eq.0) then                                          !srv 12.11.08
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then       !srv 26.05.08
           conb(ix,iy,0,3) = 0.0e0_R8                                   !srv 26.05.08
           conb(ix,iy,0,4) = 0.0e0_R8
          else
!           conb(ix,iy,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1)*ti(ix,iy)    !srv 21.10.10 {
!           conb(ix,iy,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,3) = 0.25_R8*PSch*cddb(ix,iy,1)*
     &                       kmprti(ix,iy,1)
           conb(ix,iy,0,4) = conb(ix,iy,0,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
      else
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. 
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then       !srv 26.05.08
           conb(ix,iy,0,3) = 0.0e0_R8                                   !srv 26.05.08
           conb(ix,iy,0,4) = 0.0e0_R8
          else
!            conb(ix,iy,0,3) = 0.25e0_R8*PSch*cddb(ix,iy,1)*ti(ix,iy)   !srv 21.10.10 {
!            conb(ix,iy,0,4) = 0.25e0_R8*PSch*cddb(ix,iy,1)*
!     &                      ti(bottomix(ix,iy),bottomiy(ix,iy))
           conb(ix,iy,0,3) = 0.25_R8*PSch*cddb(ix,iy,1)*
     &                       kmprti(ix,iy,1)
           conb(ix,iy,0,4) = conb(ix,iy,0,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
      endif

*   ..compute P.Sch y-flux
      iy = -1
      if(PSch_style.eq.0) then                                          !srv 12.11.08
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           conb(ix,iy,1,3) = 0.0e0_R8                                   !srv 26.05.08
           conb(ix,iy,1,4) = 0.0e0_R8
          else
!           conb(ix,iy,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0)*ti(ix,iy)    !srv 21.10.10 {
!           conb(ix,iy,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,3) = 0.25_R8*PSch*cddb(ix,iy,0)*
     &                       kmprti(ix,iy,0)
           conb(ix,iy,1,4) = conb(ix,iy,1,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
      else
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           conb(ix,iy,1,3) = 0.0e0_R8                                   !srv 26.05.08
           conb(ix,iy,1,4) = 0.0e0_R8
          else
!           conb(ix,iy,1,3) = 0.25e0_R8*PSch*cddb(ix,iy,0)*ti(ix,iy)    !srv 21.10.10 {
!           conb(ix,iy,1,4) = 0.25e0_R8*PSch*cddb(ix,iy,0)*
!     &                        ti(leftix(ix,iy),leftiy(ix,iy))
           conb(ix,iy,1,3) = 0.25_R8*PSch*cddb(ix,iy,0)*
     &                       kmprti(ix,iy,0)
           conb(ix,iy,1,4) = conb(ix,iy,1,3)                            !srv 21.10.10 }
          endif
         enddo
        enddo
      endif

      do iy = -1, ny                                                   !srv 02.01.07 {
        do ix = -1, nx
          conb(ix,iy,0:1,3) = conb(ix,iy,0:1,3)*facdrift(ix,iy)
          conb(ix,iy,0:1,4) = conb(ix,iy,0:1,4)*facdrift(ix,iy)
        enddo
      enddo                                                            !srv 02.01.07 }
      if(guard_flows.eq.0) then                                        !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            conb(ix,iy,1,:)=0.0_R8
          endif
         enddo
        enddo
      endif                                                             !xpb

csw 18feb2014 remove fluxes between guard cells 
      do ix=-1,nx
        do iy=-1,ny
          if(iffboundx(ix,iy).ne.0) then
            conb(ix,iy,0,:)=0.0_R8
          endif
          if(iffboundy(ix,iy).ne.0) then
            conb(ix,iy,1,:)=0.0_R8
          endif
        enddo
      enddo

      if(iout.ne.0) then                                                !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out(70,nx,ny,conb(-1,-1,0,0),'b2txnpsch_conbx0'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,1),'b2txnpsch_conbx1'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,2),'b2txnpsch_conbx2'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,3),'b2txnpsch_conbx3'//chns)
       call my_out(70,nx,ny,conb(-1,-1,0,4),'b2txnpsch_conbx4'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,0),'b2txnpsch_conby0'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,1),'b2txnpsch_conby1'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,2),'b2txnpsch_conby2'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,3),'b2txnpsch_conby3'//chns)
       call my_out(70,nx,ny,conb(-1,-1,1,4),'b2txnpsch_conby4'//chns)
      endif                                                             !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2txnpsch

      end subroutine b2txnpsch                                          !srv 02.04.18 }
