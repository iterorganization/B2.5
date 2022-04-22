*-----------------------------------------------------------------------
*.specification

      subroutine b2tfed (nx, ny, facdrift, fac_ExB,
     & te, po, vedia, veecrb)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo
      use b2mod_geo_corner
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     * facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny),
     * te(-1:nx,-1:ny), po(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  vedia(-1:nx,-1:ny,0:1), veecrb(-1:nx,-1:ny,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFED computes the electron drift velocities.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      integer drift_style
      real (kind=R8) ::
     * poleldr, poldidr, radeldr, raddidr
      real (kind=R8) ::
     * wpo(-1:nx,-1:ny), weight(-1:nx,-1:ny),
     * poc(-1:nx,-1:ny)
      save ncall, drift_style, poleldr, poldidr, radeldr, raddidr
*   ..procedures
      intrinsic min
      external ipgeti, ipgetr, xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/,
     * poleldr/1.0_R8/, poldidr/1.0_R8/,
     * radeldr/1.0_R8/, raddidr/1.0_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfed')
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        if (redef_gmtry.eq.0) then
          drift_style = 1
        else
          drift_style = 2
        end if
        call ipgeti ('b2tfnb_drift_style', drift_style)
        call ipgetr ('b2tfnb_poleldr', poleldr)                         !srv 25.04.00 {
        call ipgetr ('b2tfnb_radeldr', radeldr)
        call ipgetr ('b2tfnb_poldidr', poldidr)
        call ipgetr ('b2tfnb_raddidr', raddidr)                         !srv 25.04.00 }
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of te
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif
*   ..test facdrift
      call xertst (0.0_R8.le.minval(facdrift) .and.
     &     1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')

c Drift terms
      vedia=0.0_R8                                                       !xpb
      veecrb=0.0_R8                                                      !xpb
      if(maxval(facdrift).ne.0.0_R8 .or. maxval(fac_ExB).ne.0.0_R8) then !xpb

        if(drift_style.eq.0) then ! 5.0 treatment

          call intfacev(nx,ny,bottomix,bottomiy,vol,po,wpo)
          do iy = -1, ny
            do ix = -1, nx                                               !xpb
cxpb  Must check that we are not at the island X-point or at the mesh edge
              if(topiy(topix(ix,iy),min(topiy(ix,iy),ny)).ne.iy .and.
     &           topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &           rightix(ix,iy).ne.nx+1 .and. leftiy(ix,iy).ne.-2 .and.
     &           region(ix,iy,0).ne.0) then                              !xpb
cxpb  Standard case
                veecrb(ix,iy,0) = - 1.0_R8/hy1(ix,iy)*       ! signmf    !xpb !srv 11.01.13
     1            bb(ix,iy,2)/bb(ix,iy,3)**2*                            !xpb
     2           (wpo(topix(ix,iy),topiy(ix,iy))-wpo(ix,iy))*
     &            poleldr*fac_ExB(ix,iy)                                 !xpb
                vedia(ix,iy,0) = + bb(ix,iy,2)/bb(ix,iy,3)/hy1(ix,iy)*   !xpb
     1            1.0_R8/bb(ix,iy,3)*                                    !xpb
     2            2.0_R8*te(ix,iy)/bb(ix,iy,3)/qe*                       !xpb
     3           (wbbv(topix(ix,iy),topiy(ix,iy),3)-wbbv(ix,iy,3))*
     &            poldidr*facdrift(ix,iy)                                !xpb
              else                                                       !xpb
                veecrb(ix,iy,0) = 0.0_R8                                 !xpb
                vedia(ix,iy,0) = 0.0_R8                                  !xpb
              endif                                                      !xpb
            enddo
          enddo
        
          call intfaceh(nx,ny,leftix,leftiy,vol,po,wpo)
          do iy = -1, ny                                                 !xpb
            do ix = -1, nx
              if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &           bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &           region(ix,iy,0).ne.0) then
                veecrb(ix,iy,1) = 1.0_R8/hx(ix,iy)*                      !xpb
     1            bb(ix,iy,2)/bb(ix,iy,3)**2*                            !xpb
     2           (wpo(rightix(ix,iy),rightiy(ix,iy))-wpo(ix,iy))*
     &            radeldr*fac_ExB(ix,iy)                                 !xpb
                vedia(ix,iy,1) = - te(ix,iy)*bb(ix,iy,2)/                !xpb
     &           (qe*hx(ix,iy))*(2.0_R8/bb(ix,iy,3)**3)*                 !xpb
     &           (wbbr(ix,iy,3)-wbbl(ix,iy,3))*raddidr*facdrift(ix,iy)   !xpb
              endif
            enddo
          enddo

* drift velocities at walls are extrapolated from those at neighbouring cells
* drift velocities at iy=ny are extrapolated from those at iy=ny-1       !xpb
* drift velocities at iy=-1 are extrapolated from those at iy=0          !xpb
          do ix = -1, nx                                                 !xpb
            do iy = -1, ny                                               !xpb
              if(leftix(ix,iy).eq.-2) then
                veecrb(ix,iy,:) =
     .            veecrb(rightix(ix,iy),rightiy(ix,iy),:)                !xpb
                vedia(ix,iy,:) =
     .            vedia(rightix(ix,iy),rightiy(ix,iy),:)                 !xpb
              elseif(rightix(ix,iy).eq.nx+1) then                        !xpb
                veecrb(ix,iy,:) =
     .            veecrb(leftix(ix,iy),leftiy(ix,iy),:)                  !xpb
                vedia(ix,iy,:) =
     .            vedia(leftix(ix,iy),leftiy(ix,iy),:)                   !xpb
              endif
              if(topiy(ix,iy).eq.ny+1) then                              !xpb
                veecrb(ix,iy,:) = 
     .            veecrb(bottomix(ix,iy),bottomiy(ix,iy),:)              !xpb
                vedia(ix,iy,:) =
     .            vedia(bottomix(ix,iy),bottomiy(ix,iy),:)               !xpb
              elseif(bottomiy(ix,iy).eq.-2) then
                veecrb(ix,iy,:) =
     .            veecrb(topix(ix,iy),topiy(ix,iy),:)                    !xpb
                vedia(ix,iy,:) = vedia(topix(ix,iy),topiy(ix,iy),:)      !xpb
              endif
            enddo
          enddo

* turn off side drift flows in guard cells                               !xpb
          do iy = -1, ny                                                 !xpb
            do ix = -1, nx
              if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
                veecrb(ix,iy,1) = 0.0_R8                                 !xpb
                vedia(ix,iy,1) = 0.0_R8                                  !xpb
              endif
              if(bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
                veecrb(ix,iy,0) = 0.0_R8                                 !xpb
                vedia(ix,iy,0) = 0.0_R8                                  !xpb
              endif
            enddo
          enddo                                                          !xpb

       elseif (drift_style.eq.1 .or. redef_gmtry.eq.0) then     ! 5.2 drift treatment

         weight = 1.0e0_R8                                               !lk   
         call intcor(nx,ny,weight,po,poc)                                !lk 12.12.07
         do iy = 0, ny-1
           do ix = 0, nx
             if(leftix(ix,iy).ne.-2) then                                !srv 26.11.02
               veecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)/               !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1          (hy1(ix,iy)+hy1(leftix(ix,iy),leftiy(ix,iy)))* ! signmf* !srv
     2          (bb(ix,iy,2)/bb(ix,iy,3)**2                              !srv 05.02.05
     2          +bb(leftix(ix,iy),leftiy(ix,iy),2)/
     2           bb(leftix(ix,iy),leftiy(ix,iy),3)**2)*                  !srv 05.02.05
     2          (poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))              !srv 13.05.16 }
               vedia(ix,iy,0) = -poldidr*facdrift(ix,iy)* ! signmf*      !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &           (te(leftix(ix,iy),leftiy(ix,iy))                        !srv 05.03.99
     &           *bb(leftix(ix,iy),leftiy(ix,iy),2)                      !srv 05.02.05
     &          /(qe*hy1(leftix(ix,iy),leftiy(ix,iy)))                   !srv
     &           +te(ix,iy)*bb(ix,iy,2)/(qe*hy1(ix,iy)))*0.5_R8          !srv 05.02.05
     &          *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))   !srv 05.03.99
             endif
           enddo
         enddo

         do iy = 0, ny 
           do ix = 0, nx-1
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
               veecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)/                 !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     1          (hx(ix,iy)+hx(bottomix(ix,iy),bottomiy(ix,iy)))* !signmf*!srv
     2          (bb(ix,iy,2)/bb(ix,iy,3)**2                              !srv 05.02.05
     2          +bb(bottomix(ix,iy),bottomiy(ix,iy),2)/
     2           bb(bottomix(ix,iy),bottomiy(ix,iy),3)**2)*              !srv 05.02.05
     2          (poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))          !srv 13.05.16 }
               vedia(ix,iy,1) =  raddidr*facdrift(ix,iy)* ! signmf*      !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &           (te(bottomix(ix,iy),bottomiy(ix,iy))                    !srv 05.03.99
     &           *bb(bottomix(ix,iy),bottomiy(ix,iy),2)                  !srv
     &          /(qe*hx(bottomix(ix,iy),bottomiy(ix,iy)))                !srv
     &           +te(ix,iy)*bb(ix,iy,2)/(qe*hx(ix,iy)))*0.5_R8           !srv
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 05.03.99
             endif
           enddo
         enddo

       elseif (drift_style.eq.2 .and. redef_gmtry.eq.1) then             ! IYS 11.02.2020

         weight = 1.0e0_R8                                               !lk   
         call intcor(nx,ny,weight,po,poc)                                !lk 12.12.07
         do iy = 0, ny-1
           do ix = 0, nx
             if(leftix(ix,iy).ne.-2) then                                !srv 26.11.02
               veecrb(ix,iy,0) = - poleldr*fac_ExB(ix,iy)                !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     2              *(poc(topix(ix,iy),topiy(ix,iy))-poc(ix,iy))*        !srv 13.05.16 }
     &           wbbl(ix,iy,2)/(wbbl(ix,iy,3))**2/hyl(ix,iy)/qc(ix,iy)   ! IYS 04.03.2019
               vedia(ix,iy,0) = -poldidr*facdrift(ix,iy)* ! signmf*      !srv 07.06.00 17.05.99 25.04.00 05.02.05 11.01.13
     &           0.5_R8 * (te(leftix(ix,iy),leftiy(ix,iy)) + te(ix,iy))
     &           * wbbl(ix,iy,2) / (qe * hyl(ix,iy) * qc(ix,iy))         ! IYS 04.03.2019
     &           *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))  !srv 05.03.99
             endif
           enddo
         enddo

         do iy = 0, ny 
           do ix = 0, nx-1
             if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 26.11.02
               veecrb(ix,iy,1) = radeldr*fac_ExB(ix,iy)                  !srv 17.05.99 25.08.99 07.06.00 05.02.05 13.05.16 {
     2           *(poc(rightix(ix,iy),rightiy(ix,iy))-poc(ix,iy))        !srv 13.05.16 }
     &           * wbbv(ix,iy,2) / (wbbv(ix,iy,3))**2 / hxb(ix,iy)       ! IYS 04.03.2019
               vedia(ix,iy,1) =  raddidr*facdrift(ix,iy)* ! signmf*      !srv 11.01.13 07.06.00 17.05.99 25.04.00
     &           (te(bottomix(ix,iy),bottomiy(ix,iy)) + te(ix,iy))*
     &           0.5_R8*wbbv(ix,iy,2) / (qe * hxb(ix,iy))                ! IYS 04.03.2019
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !srv 05.03.99
             endif
           enddo
         enddo

       else         ! drift_style                                    ! IYS 11.02.2020 {                            

         write(*,fmt="(2(a20,i4))") 'drift_style = ',drift_style,
     &      'redef_gmtry = ',redef_gmtry
         call xerrab (
     &   'faulty definition of drift_style and redef_gmtry, ' // 
     &   'or logical mistake occurred')

       endif

      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfed

      end subroutine b2tfed
