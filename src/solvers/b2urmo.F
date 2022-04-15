*-----------------------------------------------------------------------
*.specification

      subroutine b2urmo (nx, ny, hx, hy, vol, qs,  isb, ub, rob,
     .  smb, flub, cvsb, resmb, fmb)                                    !srv 02.06.11
      use b2mod_types
      use b2mod_math
      use b2mod_indirect
      use b2mod_balance                                                 !djm Jan2017
     & , only : fmo_flua, fmo_cvsa, fmo_hybr, balance_netcdf
      use b2mod_openmp
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     *   ub(-1:nx,-1:ny), rob(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3),
     &  flub(-1:nx,-1:ny,0:1,0:1), cvsb(-1:nx,-1:ny,0:1,0:1),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), vol(-1:nx,-1:ny), 
     &  qs(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resmb(-1:nx,-1:ny), fmb(-1:nx,-1:ny,0:1,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URMO computes the residual of the parallel momentum balance
*     equation for one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of the parallel momentum
*     balance equation for one species, integrated over the cells of
*     the B2 mesh. The equation has the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resmb = rhs-div(flux) .
*     The flux will be computed in the array fmb; it has the form:
*       fmbx = flubx*ub-cvsbx*(d/dx)ub
*       fmby = fluby*up-cvsby*(d/dy)ub
*     The right hand side has a linearised representation with
*     coefficients given by smb. This linearisation may include a
*     contribution from an implicit discretisation in time.
*
*     The input arguments smb, flub and cvsb contain appropriate
*     geometric coefficients so that in the present routine the
*     expression for the integrated residual on any cell does not
*     involve any geometric factors. However, the correct computation
*     of (radial) differences of ub across the faces requires knowledge
*     of geometric factors.
*
*     Since the introduction of the 9-point stencil, the dimensions of 
*     the flux arrays have increase to f(-1:nx,-1:ny,0:1,0:1), where
*     f(-1:nx,-1:ny,0,0) specifies the poloidal flow across an x-face
*     f(-1:nx,-1:ny,0,1) specifies the radial   flow across an x-face
*     f(-1:nx,-1:ny,1,0) specifies the poloidal flow across an y-face
*     f(-1:nx,-1:ny,1,1) specifies the radial   flow across an y-face
*
*     Similarly, the conductivity arrays have increased dimensions, so
*     cvsb(-1:nx,-1:ny,0,0) relates to (d/dx)ub on an x-face
*     cvsb(-1:nx,-1:ny,0,1) relates to (d/dy)ub on an x-face
*     cvsb(-1:nx,-1:ny,1,0) relates to (d/dx)ub on an y-face
*     cvsb(-1:nx,-1:ny,1,1) relates to (d/dy)ub on an y-face
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0             !srv 06.07.06
      integer, save :: iout_b2npmo = 0, iout_b2wdat = 0                 !srv 11.09.09
      character chns*3                                                  !srv 02.06.11
      save ncall
      !$OMP THREADPRIVATE(ncall, no_Ptncr_x_co, no_Ptncr_y_co)
      !$OMP THREADPRIVATE(iout_b2wdat, iout_b2npmo)
      real (kind=R8) ::
     *   con(-1:nx,-1:ny,0:1,0:1), flv(-1:nx,-1:ny,0:1,0:1), flo,       !srv 11.09.09
     *   flo9, con9, dubr(-1:nx,-1:ny,0:1),
     *   wght(-1:nx,-1:ny), ubc(-1:nx,-1:ny),
     *   conhyb !djm Jan2018 viscosity coefficient associated with hybrid correction
*   ..procedures
      external xertst
      external b2xvsg, b2xvff, diff_r
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2urmo')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgeti ('b2wdat_iout', iout_b2wdat)                              !iys 17.04.17
        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)
        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of rob, smb1, smb3
       call b2xvsg (n2, rob, 1, 'rob', '.gt.')
       call b2xvsg (n2, smb(-1,-1,1), 1, 'smb1', '.le.')
       call b2xvsg (n2, smb(-1,-1,3), 1, 'smb3', '.le.')
*    ..test sign of cvsb
       call b2xvsg (n2, cvsb(-1,-1,0,0), 1, 'cvsb', '.ge.')
       call b2xvsg (n2, cvsb(-1,-1,1,1), 1, 'cvsb', '.ge.')
*    ..test edge values of flub, cvsb
       call b2xvff (nx, ny, flub, 'flub')
       call b2xvff (nx, ny, cvsb, 'cvsb')
      endif

*   ..compute radial differences of ub for 9-point terms
      wght = 1.0_R8
      wght = vol
      call diff_r(nx, ny, 0, hx, hy, wght, qs, ub, ubc, dubr)

* ..compute flux
*   ..compute parallel momentum transport fmbx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         fmb(ix,iy,0,0:1) = 0.0_R8
         con(ix,iy,0,0:1) = 0.0_R8
         flv(ix,iy,0,0:1) = 0.0_R8
        else
         flo = flub(ix,iy,0,0)
         if(no_Ptncr_x_co.eq.0) then                                  !sv 06.07.06 {
          con(ix,iy,0,0) = hybr(flo,cvsb(ix,iy,0,0))
         else
          con(ix,iy,0,0) = cvsb(ix,iy,0,0)
         endif                                                        !sv 06.07.06 }
cwdk     add computation of 9-point quantities
         flo9 = flub(ix,iy,0,1)
         con9 = sign(hybr(flo9,cvsb(ix,iy,0,1)) - 0.5_R8*abs(flo9),
     &                         cvsb(ix,iy,0,1))
         flv(ix,iy,0,0) = - con(ix,iy,0,0)*
     &         (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
         fmb(ix,iy,0,0) = flo*(ub(ix,iy)
     &         +ub(leftix(ix,iy),leftiy(ix,iy)))/2.0_R8 
     &         +flv(ix,iy,0,0)
         fmb(ix,iy,0,1) = upwind(flo9,
     &          ub(leftix(ix,iy),leftiy(ix,iy)),ub(ix,iy)) -
     &          con9*dubr(ix,iy,0)
        endif
        if (balance_netcdf.ne.0) then !djm Jan2018
         if(leftix(ix,iy).eq.-2) then
          fmo_cvsa(ix,iy,0,:,isb) = 0.0_R8
          fmo_flua(ix,iy,0,:,isb) = 0.0_R8
         else
          if(no_Ptncr_x_co.eq.0) then
           conhyb = con(ix,iy,0,0)-abs(flo)/2.0_R8-cvsb(ix,iy,0,0)
           fmo_flua(ix,iy,0,0,isb) =
     &       (flo-abs(flo))/2.0_R8*ub(ix,iy)+
     &       (flo+abs(flo))/2.0_R8*ub(leftix(ix,iy),leftiy(ix,iy))
           fmo_flua(ix,iy,0,1,isb) = 0.0_R8
           fmo_cvsa(ix,iy,0,0,isb) = -cvsb(ix,iy,0,0)*
     &       (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
           fmo_cvsa(ix,iy,0,0,isb) = -cvsb(ix,iy,0,1)*
     &       (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
           fmo_hybr(ix,iy,0,0,isb) = -conhyb*
     &       (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
           fmo_hybr(ix,iy,0,1,isb) = 0.0_R8
          else
           fmo_flua(ix,iy,0,0,isb) = fmb(ix,iy,0,0)-flv(ix,iy,0,0)
           fmo_flua(ix,iy,0,1,isb) = fmb(ix,iy,0,1)-flv(ix,iy,0,1)
           fmo_cvsa(ix,iy,0,0,isb) = flv(ix,iy,0,0)
           fmo_cvsa(ix,iy,0,1,isb) = flv(ix,iy,0,1)
           fmo_hybr(ix,iy,0,:,isb) = 0.0_R8
          endif
         endif
        endif
       enddo
      enddo
*   ..compute parallel momentum transport fmby
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         fmb(ix,iy,1,0:1) = 0.0_R8
         con(ix,iy,1,0:1) = 0.0_R8
         flv(ix,iy,1,0:1) = 0.0_R8
        else
         flo = flub(ix,iy,1,1)
         if(no_Ptncr_y_co.eq.0) then                                  !sv 06.07.06 {
          con(ix,iy,1,1) = hybr(flo,cvsb(ix,iy,1,1))-0.5_R8*abs(flo)
         else
          con(ix,iy,1,1) = cvsb(ix,iy,1,1)
         endif                                                        !sv 06.07.06 }
         flv(ix,iy,1,1) = - con(ix,iy,1,1)*dubr(ix,iy,1)
         fmb(ix,iy,1,1) = upwind(flo,
     &          ub(bottomix(ix,iy),bottomiy(ix,iy)),ub(ix,iy)) +
     &          flv(ix,iy,1,1)
         fmb(ix,iy,1,0) = 0.0_R8
        endif
        if (balance_netcdf.ne.0) then !djm Jan2018
         if(bottomiy(ix,iy).eq.-2) then
          fmo_cvsa(ix,iy,1,:,isb) = 0.0_R8
          fmo_flua(ix,iy,1,:,isb) = 0.0_R8
         else
          if(no_Ptncr_y_co.eq.0) then
           conhyb = con(ix,iy,1,1)-abs(flo)/2.0_R8-cvsb(ix,iy,1,1)
           fmo_flua(ix,iy,1,1,isb) =
     &       (flo-abs(flo))/2.0_R8*ub(ix,iy)+
     &       (flo+abs(flo))/2.0_R8*ub(bottomix(ix,iy),bottomiy(ix,iy))
           fmo_cvsa(ix,iy,1,1,isb) = -cvsb(ix,iy,1,1)*
     &       (ub(ix,iy)-ub(bottomix(ix,iy),bottomiy(ix,iy)))
           fmo_hybr(ix,iy,1,1,isb) = -conhyb*
     &       (ub(ix,iy)-ub(bottomix(ix,iy),bottomiy(ix,iy)))
          else
           fmo_flua(ix,iy,1,0,isb) = fmb(ix,iy,1,0)-flv(ix,iy,1,0)
           fmo_flua(ix,iy,1,1,isb) = fmb(ix,iy,1,1)-flv(ix,iy,1,1)
           fmo_cvsa(ix,iy,1,0,isb) = flv(ix,iy,1,0)
           fmo_cvsa(ix,iy,1,1,isb) = flv(ix,iy,1,1)
           fmo_hybr(ix,iy,1,:,isb) = 0.0_R8
          endif
         endif
        endif
       enddo
      enddo
      if(iout_b2npmo.eq.1 .or. iout_b2wdat.eq.4) then
       if (in_parallel()) then
        !$OMP MASTER
        if (ncall==0) then
          write(*,*) 'b2urmo OpenMP warning: no file output ',
     &               'in parallel mode'
        end if
        !$OMP END MASTER
       else
        write (chns,'(i3.3)') isb
        call my_out(70,nx,ny,flv(-1,-1,0,0),
     &      'b2urmo_etaPat_graduax'//chns)
        call my_out(70,nx,ny,flv(-1,-1,1,1),
     &      'b2urmo_etaPat_graduay'//chns)
!       call my_out(70,nx,ny,con(-1,-1,0,0),'b2urmo_etaPatx')
!       call my_out(70,nx,ny,con(-1,-1,1,1),'b2urmo_etaPaty')
       end if
      endif      
* ..compute residual
*   ..contribution from right hand side
      do iy = -1, ny
       do ix = -1, nx
        resmb(ix,iy) = smb(ix,iy,0)+
     &    smb(ix,iy,1)*ub(ix,iy)+
     &    smb(ix,iy,2)*rob(ix,iy)+
     &    smb(ix,iy,3)*rob(ix,iy)*ub(ix,iy)
       enddo
      enddo
*   ..contribution from fmbx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         resmb(ix,iy) = resmb(ix,iy)
     &    -fmb(rightix(ix,iy),rightiy(ix,iy),0,0)
     &    -fmb(rightix(ix,iy),rightiy(ix,iy),0,1)
        elseif(rightix(ix,iy).eq.nx+1) then
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,0,0)+fmb(ix,iy,0,1)
        else
         resmb(ix,iy) = resmb(ix,iy)
     &    +fmb(ix,iy,0,0)-fmb(rightix(ix,iy),rightiy(ix,iy),0,0)
     &    +fmb(ix,iy,0,1)-fmb(rightix(ix,iy),rightiy(ix,iy),0,1)
        endif
       enddo
      enddo
*   ..contribution from fmby
      do iy = -1, ny
       do ix = -1, nx
        if (bottomiy(ix,iy).eq.-2) then
cxpb  cells with no bottom neighbour
          resmb(ix,iy) = resmb(ix,iy)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1,0)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1,1)
        elseif (topiy(ix,iy).eq.ny+1) then
cxpb  cells with no top neighbour
          resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1,0)+fmb(ix,iy,1,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then
cxpb  standard case
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1,0)+fmb(ix,iy,1,1)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1,0)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1)
     &   then
cxpb  diverted island X-point case
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1,0)+fmb(ix,iy,1,1)
        endif
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2urmo

      end subroutine b2urmo
