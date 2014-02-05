*-----------------------------------------------------------------------
*.specification

      subroutine b2ursc (nx, ny, fun, sfn, ffn, resfn)
      use b2mod_types
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   fun(-1:nx,-1:ny), sfn(-1:nx,-1:ny,0:1), 
     &  ffn(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   resfn(-1:nx,-1:ny)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URSC computes the residual of a scalar conservation equation.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the residual of a scalar conservation
*     equation, integrated over cells of the b2 mesh. The equation has
*     the form
*       div(flux) = rhs ,
*     and the residual is defined as
*       resfn = rhs-div(flux) .
*     The flux is specified by ffn, and the right hand side has the
*     linearised representation
*       rhs(,) = sfn(,,0)+sfn(,,1)*fun(,) .
*     This linearisation may include a contribution from an implicit
*     discretisation in time.
*
*     The input argument ffn specifies fluxes integrated over cell
*     faces, and sfn specifies source coefficients integrated over
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
      call subini ('b2ursc')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of sfn
       call b2xvsg (n2, sfn(-1,-1,1), 1, 'sfn', '.le.')
*    ..test edge values of ffn
       call b2xvff (nx, ny, ffn, 'ffn')
      endif

* ..compute residual
*   ..contribution from right hand side
      do iy = -1, ny
       do ix = -1, nx
        resfn(ix,iy) = sfn(ix,iy,0)+sfn(ix,iy,1)*fun(ix,iy)
       enddo
      enddo
*   ..contribution from ffnx
      do iy = -1, ny
       do ix = -1, nx
        if(leftix(ix,iy).eq.-2) then
         resfn(ix,iy) = resfn(ix,iy)
     &                  -ffn(rightix(ix,iy),rightiy(ix,iy),0)
        elseif(rightix(ix,iy).eq.nx+1) then
         resfn(ix,iy) = resfn(ix,iy)+ffn(ix,iy,0)
        else
         resfn(ix,iy) = resfn(ix,iy)+ffn(ix,iy,0)
     &                  -ffn(rightix(ix,iy),rightiy(ix,iy),0)
        endif
       enddo
      enddo
*   ..contribution from ffny
      do iy = -1, ny
       do ix = -1, nx
        if (bottomiy(ix,iy).eq.-2) then
cxpb  cells with no bottom neighbour
         resfn(ix,iy) = resfn(ix,iy)
     &     -ffn(topix(ix,iy),topiy(ix,iy),1)
        elseif (topiy(ix,iy).eq.ny+1) then
cxpb  cells with no top neighbour
          resfn(ix,iy) = resfn(ix,iy)+ffn(ix,iy,1)
        elseif (topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then
cxpb  standard case
          resfn(ix,iy) = resfn(ix,iy)+ffn(ix,iy,1)
     &                  -ffn(topix(ix,iy),topiy(ix,iy),1)
         elseif (topiy(topix(ix,iy),topiy(ix,iy)).eq.iy .and. ny.gt.1) 
     &    then
cxpb  diverted island X-point case
          resfn(ix,iy) = resfn(ix,iy)+ffn(ix,iy,1)
        endif
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2ursc

      end subroutine b2ursc
