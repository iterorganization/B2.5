*-----------------------------------------------------------------------
*.specification

      subroutine b2tnpsch (nx, ny, isb, guard_flows, facdrift, zb,   !xpb !iyv 06.02.14 02.04.18 {
     &   nb, ti, cddb, fnbPSch)                       !srv 09.01.01 13.10.06
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb, guard_flows                                   !xpb !srv 08.09.01 !lk 02.06.11
      real (kind=R8) ::
 !    *  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  nb(-1:nx,-1:ny),!  hzl(-1:nx,-1:ny),hzb(-1:nx,-1:ny),              !srv 21.10.10
!     &  qz(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2), qc(-1:nx,-1:ny),
     &  facdrift(-1:nx,-1:ny), zb(-1:nx,-1:ny),
     &  ti(-1:nx,-1:ny), cddb(-1:nx,-1:ny,0:1)                           !srv 13.10.06 !iyv 06.02.14
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fnbPSch(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TNPSCH computes the Pfirsch-Schlueter particle fluxes of one atomic
*     species between neighbouring cells
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle flux, fnbPSch, of one atomic
*     species
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, ib, area_fix,
     &        ixl, iyl, ixb, iyb                                         !srv 13.10.06
      integer, save :: iout = 0                                          !srv 02.04.18
      integer, save :: PSch_style                                    !srv 12.11.08
      character chns*3                                                   !srv 16.09.02
      real (kind=R8), save :: PSch = 1.0_R8                              !srv 17.01.07
      real (kind=R8) ::
     *   nbcor(-1:nx,-1:ny),  !srv 01.10.99 !lk 09.06.11
     *   nbti(-1:nx,-1:ny), nbtic(-1:nx,-1:ny)                           !srv 08.10.04
      real (kind=R8) :: weight(-1:nx,-1:ny)                              !lk
      save ncall, area_fix
*   ..procedures
      external xertst, ipgeti, ipgetr                                    !xpb
      external b2xvsg, b2xvff, intfaceh, intfacev
*   ..initialisation
      data ncall/0/, area_fix/3/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tnpsch')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2news_area_fix', area_fix)
        call ipgeti ('b2tnpsch_iout', iout)                              !srv 02.04.18
        call ipgetr ('b2tfnb_PSch', PSch)                                !srv 17.01.07
        if (redef_gmtry.eq.0) then
         PSch_style = 1
        else
         PSch_style = 2
        end if
        call ipgeti ('b2tfnb_PSch_style', PSch_style)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb, pb, cdnb, cdpb
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
      endif

      write(chns,'(i3.3)') isb
      fnbPSch = 0.0_R8                                                   !srv 13.10.06

*   ..compute P.Sch x-flux                                               !srv 13.10.06 {
      if(PSch_style.eq.0) then                                          !srv 12.11.08
        nbti=nb*ti
        call intcor(nx,ny,vol,nbti,nbtic)                                !srv 08.10.04
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or.
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then        !srv 26.05.08
           fnbPSch(ix,iy,0)=0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
           fnbPSch(ix,iy,0) = 0.25e0_R8*PSch*(
     &      cddb(topix(ix,iy),topiy(ix,iy),1)*
     &       (nb(topix(ix,iy),topiy(ix,iy))*
     &        ti(topix(ix,iy),topiy(ix,iy))- nb(ix,iy)*ti(ix,iy))+
     &      cddb(topix(ixl,iyl),topiy(ixl,iyl),1)*
     &       (nb(topix(ixl,iyl),topiy(ixl,iyl))*
     &        ti(topix(ixl,iyl),topiy(ixl,iyl))-
     &        nb(ixl,iyl)*ti(ixl,iyl))+
     &      cddb(ixl,iyl,1)*
     &       (nb(ixl,iyl)*ti(ixl,iyl)
     &                   -nb(bottomix(ixl,iyl),bottomiy(ixl,iyl))*
     &                    ti(bottomix(ixl,iyl),bottomiy(ixl,iyl)))+
     &      cddb(ix,iy,1)*
     &       (nb(ix,iy)*ti(ix,iy)-nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &                         ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &                                       )
          endif
         enddo
        enddo
      else if(PSch_style.eq.1) then
        weight = 1.0e0_R8                                               !lk
        nbti=nb*ti
        call intcor(nx,ny,weight,nbti,nbtic)                            !lk 12.12.07
        call intcor(nx,ny,weight,nb,nbcor)                              !lk 09.06.11
        do iy = 0, ny-1
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or.
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then       !srv 26.05.08
           fnbPSch(ix,iy,0)=0.0e0_R8
          else
           ixl = leftix(ix,iy)
           iyl = leftiy(ix,iy)
!           fnbPSch(ix,iy,0) = PSch*                                    !lk {
!     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixl,iyl)
!     &      +bzb(ixl,iyl)*hz(ixl,iyl)*vol(ix,iy))
!     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                  !lk }
           fnbPSch(ix,iy,0) = PSch*                                     !lk {
     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixl,iyl)
     &      +bzb(ixl,iyl)*(vol(ixl,iyl)/gs(ixl,iyl,2))*vol(ix,iy))
     &      /(vol(ixl,iyl)+vol(ix,iy))/(qe*zb(ix,iy))                   !lk }
     &      *(nb(ix,iy)+nb(ixl,iyl))/(nbcor(ix,iy)+                     !lk 09.06.11
     &      nbcor(topix(ix,iy),topiy(ix,iy)))                           !lk 09.06.11
          endif
         enddo
        enddo
      else                                                              !nn
        nbti=nb*ti
        weight = 1.0e0_R8
        call intcor(nx,ny,weight,nbti,nbtic)
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or.
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1 ) then
              fnbPSch(ix,iy,0)= 0.0e0_R8
            else
              fnbPSch(ix,iy,0) = PSch*                                  !
     &      (nbtic(topix(ix,iy),topiy(ix,iy))-nbtic(ix,iy))*bzbl(ix,iy)/
     &      (qe*zb(ix,iy))*hzl(ix,iy)
            endif
          enddo
        enddo
      endif

*   ..compute P.Sch y-flux
      iy = -1
      fnbPSch(:,iy,1)=0.0e0_R8
      if(PSch_style.eq.0) then                                          !srv 12.11.08
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1)= 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
           fnbPSch(ix,iy,1) = -0.25e0_R8*PSch*(
     &    cddb(rightix(ix,iy),rightiy(ix,iy),0)*
     &     (nb(rightix(ix,iy),rightiy(ix,iy))*
     &      ti(rightix(ix,iy),rightiy(ix,iy))-nb(ix,iy)*ti(ix,iy))+
     &    cddb(ix,iy,0)*
     &    (nb(ix,iy)*ti(ix,iy)-nb(leftix(ix,iy),leftiy(ix,iy))*
     &                         ti(leftix(ix,iy),leftiy(ix,iy)))+
     &    cddb(bottomix(ix,iy),bottomiy(ix,iy),0)*
     &     (nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &      ti(bottomix(ix,iy),bottomiy(ix,iy))
     &                         -nb(leftix(ixb,iyb),leftiy(ixb,iyb))*
     &                          ti(leftix(ixb,iyb),leftiy(ixb,iyb)))+
     &    cddb(rightix(ixb,iyb),rightiy(ixb,iyb),0)*
     &     (nb(rightix(ixb,iyb),rightiy(ixb,iyb))*
     &      ti(rightix(ixb,iyb),rightiy(ixb,iyb))
     &                  -nb(bottomix(ix,iy),bottomiy(ix,iy))*
     &                   ti(bottomix(ix,iy),bottomiy(ix,iy)))
     &                                         )

          endif
         enddo
        enddo
      else if(PSch_style.eq.1) then
        do iy = 0, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
           fnbPSch(ix,iy,1)= 0.0e0_R8
          else
           ixb = bottomix(ix,iy)
           iyb = bottomiy(ix,iy)
!           fnbPSch(ix,iy,1) = -PSch*                                   !lk {
!     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
!     &      (bzb(ix,iy)*hz(ix,iy)*vol(ixb,iyb)
!     &      +bzb(ixb,iyb)*hz(ixb,iyb)*vol(ix,iy))
!     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                  !lk }
           fnbPSch(ix,iy,1) = -PSch*                                    !lk {
     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
     &      (bzb(ix,iy)*(vol(ix,iy)/gs(ix,iy,2))*vol(ixb,iyb)
     &      +bzb(ixb,iyb)*(vol(ixb,iyb)/gs(ixb,iyb,2))*vol(ix,iy))
     &      /(vol(ixb,iyb)+vol(ix,iy))/(qe*zb(ix,iy))                   !lk }
     &      *(nb(ix,iy)+nb(ixb,iyb))/(nbcor(ix,iy)+                     !lk 09.06.11
     &      nbcor(rightix(ix,iy),rightiy(ix,iy)))                       !lk 09.06.11
          endif
         enddo
        enddo
      else                                                              !nn
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fnbPSch(ix,iy,1)= 0.0e0_R8
            else
              fnbPSch(ix,iy,1) = -PSch*                                    !
     &      (nbtic(rightix(ix,iy),rightiy(ix,iy))-nbtic(ix,iy))*
     &      bzbv(ix,iy)/(qe*zb(ix,iy))*hzb(ix,iy)
            endif
          enddo
        enddo
      endif

      do iy = -1, ny                                                    !srv 02.01.07 {
        do ix = -1, nx
          fnbPSch(ix,iy,0:1) = fnbPSch(ix,iy,0:1)*facdrift(ix,iy)
        enddo
      enddo                                                             !srv 02.01.07 }

      if(guard_flows.eq.0) then                                         !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fnbPSch(ix,iy,1)=0.0_R8
          endif
         enddo
        enddo

cxpb turn off poloidal flows in core boundary
        if (boundary_namelist.eq.1) then
          do ib = 1, nbc
            if(core_boundary(ib)) then
              do ix = 1, bc_list_size(ib)
                fnbPSch(bc_list_x(ix,ib),bc_list_y(ix,ib),0) = 0.0_R8
              enddo
            endif
          enddo
        else
          do ix = -1, nx
            if (on_closed_surface(ix,-1)) then
              fnbPSch(ix,-1,0) = 0.0_R8
            endif
          enddo
        endif

      elseif(guard_flows.eq.1) then                                     !srv 20.06.08 {
*   .. do nothing
      elseif(guard_flows.eq.2) then
*     turn off flows between boundary cells in west and east
        do iy = -1, ny
         do ix = -1, nx
          if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
            fnbPSch(ix,iy,1)        = 0.0_R8
          endif
         enddo
        enddo
*     turn off flows between boundary cells in south and in north
        do iy = -1, ny, ny+1
          fnbPSch(:,iy,0)         = 0.0_R8
        enddo                                                           !srv 20.06.08 }
      endif                                                             !xpb

csw 18feb2014 remove fluxes between guard cells
      do ix=-1,nx
        do iy=-1,ny
          if(iffboundx(ix,iy).ne.0) then
            fnbPSch(ix,iy,0)=0.0_R8
          endif
          if(iffboundy(ix,iy).ne.0) then
            fnbPSch(ix,iy,1)=0.0_R8
          endif
        enddo
      enddo

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            fnbpsch(ix,iy,0)=0.0_R8
          end do
        end do
      end do

      if(iout.ne.0) then                                                !srv 18.06.08 {
       write (chns,'(i3.3)') isb
       call my_out(70,nx,ny,fnbPSch(-1,-1,0),'b2tnpsch_fnbPSchx'//chns)
       call my_out(70,nx,ny,fnbPSch(-1,-1,1),'b2tnpsch_fnbPSchy'//chns)
      endif                                                             !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tnpsch

      end subroutine b2tnpsch                                            !srv 02.04.18 }
