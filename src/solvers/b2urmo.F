*-----------------------------------------------------------------------
*.specification

      subroutine b2urmo (nx,ny,isb,ub, rob, smb, flub, cvsb, resmb, fmb)  !srv 02.06.11
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb
      real (kind=R8) ::
     *   ub(-1:nx,-1:ny), rob(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3),
     &  flub(-1:nx,-1:ny,0:1), cvsb(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resmb(-1:nx,-1:ny), fmb(-1:nx,-1:ny,0:1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URMO computes the residual of the momentum balance equation for
*     one atomic species.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of the momentum balance
*     equation for one species, integrated over the cells of
*     the b2 mesh. The equation has the form
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
*     involve any geometric factors.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0           !sv 06.07.06
      integer, save :: iout_b2npmo = 0                            !srv 11.09.09
      character chns*2                                            !srv 02.06.11
      save ncall
      real (kind=R8) ::
     *   con(-1:nx,-1:ny,0:1), flv(-1:nx,-1:ny,0:1), flo, dif     !srv 11.09.09
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   hybr
      external subini, subend, xertst
      external b2xvsg, b2xvff
      hybr(flo,dif) = max(dif,abs(flo/2),sqrt(dif**2+flo**2/4)) !!!
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
       call b2xvsg (n2*2, cvsb, 1, 'cvsb', '.ge.')
*    ..test edge values of flub, cvsb
       call b2xvff (nx, ny, flub, 'flub')
       call b2xvff (nx, ny, cvsb, 'cvsb')
      endif

* ..compute flux
*   ..compute momentum transport fmbx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         fmb(ix,iy,0) = 0
         con(ix,iy,0) = 0
         flv(ix,iy,0) = 0
        else
         flo = flub(ix,iy,0)
         if(no_Ptncr_x_co.eq.0) then                                  !sv 06.07.06 {
          con(ix,iy,0) = hybr(flo,cvsb(ix,iy,0))
         else
          con(ix,iy,0) = cvsb(ix,iy,0)
         endif                                                        !sv 06.07.06 }
         flv(ix,iy,0) = - con(ix,iy,0)*
     &         (ub(ix,iy)-ub(leftix(ix,iy),leftiy(ix,iy)))
         fmb(ix,iy,0) = flo*(ub(ix,iy)
     &         +ub(leftix(ix,iy),leftiy(ix,iy)))/2 + flv(ix,iy,0)
        endif
       enddo
      enddo
*   ..compute momentum transport fmby
      do iy = -1, ny
       do ix = -1, nx
        if(bottomiy(ix,iy).eq.-2) then
         fmb(ix,iy,1) = 0
         con(ix,iy,1) = 0
         flv(ix,iy,1) = 0
        else
         flo = flub(ix,iy,1)
         if(no_Ptncr_y_co.eq.0) then                                  !sv 06.07.06 {
          con(ix,iy,1) = hybr(flo,cvsb(ix,iy,1))
         else
          con(ix,iy,1) = cvsb(ix,iy,1)
         endif                                                        !sv 06.07.06 }
         flv(ix,iy,1) = - con(ix,iy,1)*
     &            (ub(ix,iy)-ub(bottomix(ix,iy),bottomiy(ix,iy)))
         fmb(ix,iy,1) = flo*(ub(bottomix(ix,iy),bottomiy(ix,iy))
     &         +ub(ix,iy))/2 + flv(ix,iy,1)
        endif
       enddo
      enddo
      if(iout_b2npmo.eq.1) then
       write (chns,'(i2.2)') isb
       call my_out(70,nx,ny,flv(-1,-1,0),'b2urmo_etaPat_graduax'//chns)
       call my_out(70,nx,ny,flv(-1,-1,1),'b2urmo_etaPat_graduay'//chns)
!       call my_out(70,nx,ny,con(-1,-1,0),'b2urmo_etaPatx')
 !      call my_out(70,nx,ny,con(-1,-1,1),'b2urmo_etaPaty')
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
     &    -fmb(rightix(ix,iy),rightiy(ix,iy),0)
        elseif(rightix(ix,iy).eq.nx+1) then
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,0)
        else
         resmb(ix,iy) = resmb(ix,iy)
     &    +fmb(ix,iy,0)-fmb(rightix(ix,iy),rightiy(ix,iy),0)
        endif
       enddo
      enddo
*   ..contribution from fmby
      do iy = -1, ny
       do ix = -1, nx
        if (bottomiy(ix,iy).eq.-2) then
cxpb  cells with no bottom neighbour
          resmb(ix,iy) = resmb(ix,iy)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1)
        elseif (topiy(ix,iy).eq.ny+1) then
cxpb  cells with no top neighbour
          resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).ne.iy)  then
cxpb  standard case
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1)
     &     -fmb(topix(ix,iy),topiy(ix,iy),1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1)
     &   then
cxpb  diverted island X-point case
         resmb(ix,iy) = resmb(ix,iy)+fmb(ix,iy,1)
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
