*-----------------------------------------------------------------------
*.specification

      subroutine b2ursd (nx, ny, den, fun, sdf, fdf, resdf)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   den(-1:nx,-1:ny), fun(-1:nx,-1:ny), sdf(-1:nx,-1:ny,0:3),
     &  fdf(-1:nx,-1:ny,0:1,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resdf(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URSD computes the residual of a scalar conservation equation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of a scalar conservation
*     equation, integrated over cells of the b2 mesh. The equation has
*     the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resdf = rhs-div(flux) .
*     The flux is specified by fdf, and the right hand side has the
*     linearised representation
*       rhs(,) = sdf(,,0)+
*         sdf(,,1)*fun(,)+
*         sdf(,,2)*den(,)+
*         sdf(,,3)*den(,)*fun(,) .
*     This linearisation may include a contribution from an implicit
*     discretisation in time.
*
*     The input argument fdf specifies fluxes integrated over cell
*     faces, and sdf specifies source coefficients integrated over
*     cell volumes. Therefore, the expression for the integrated
*     residual on any cell does not involve any geometric factors.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg, b2xvff
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2ursd')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of den, sdf
       call b2xvsg (n2, den, 1, 'den', '.ge.')
       call b2xvsg (n2, sdf(-1,-1,1), 1, 'sdf1', '.le.')
       call b2xvsg (n2, sdf(-1,-1,3), 1, 'sdf3', '.le.')
*    ..test edge values of fdf
       call b2xvff (nx, ny, fdf, 'fdf')
      endif

* ..compute residual
*   ..contribution from right hand side
      do iy = -1, ny
       do ix = -1, nx
        resdf(ix,iy) = sdf(ix,iy,0)+
     &    sdf(ix,iy,1)*fun(ix,iy)+
     &    sdf(ix,iy,2)*den(ix,iy)+
     &    sdf(ix,iy,3)*den(ix,iy)*fun(ix,iy)
       enddo
      enddo
*   ..contribution from fdfx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         resdf(ix,iy) = resdf(ix,iy)
     &    -fdf(rightix(ix,iy),rightiy(ix,iy),0,0)
        elseif(rightix(ix,iy).eq.nx+1) then    
         resdf(ix,iy) = resdf(ix,iy)+fdf(ix,iy,0,0)
        else
         resdf(ix,iy) = resdf(ix,iy)+fdf(ix,iy,0,0)
     &    -fdf(rightix(ix,iy),rightiy(ix,iy),0,0)
        endif
       enddo
      enddo
*   ..contribution from fdfy
      do iy = -1, ny-1
       do ix = -1, nx
        if (bottomiy(ix,iy).eq.-2) then
cxpb  cells with no bottom neighbour
         resdf(ix,iy) = resdf(ix,iy)
     &     -fdf(topix(ix,iy),topiy(ix,iy),1,1)
        elseif (topiy(ix,iy).eq.ny+1) then
cxpb  cells with no top neighbour
          resdf(ix,iy) = resdf(ix,iy)+fdf(ix,iy,1,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then
cxpb  standard case
         resdf(ix,iy) = resdf(ix,iy)+fdf(ix,iy,1,1)
     &      -fdf(topix(ix,iy),topiy(ix,iy),1,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1) 
     &   then
cxpb  diverted island X-point case
         resdf(ix,iy) = resdf(ix,iy)+fdf(ix,iy,1,1)
        endif
       enddo
      enddo
      do ix = -1, nx
       if (bottomiy(ix,ny).ne.-2)
     >  resdf(ix,ny) = resdf(ix,ny)+fdf(ix,ny,1,1)
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ursd

      end subroutine b2ursd
