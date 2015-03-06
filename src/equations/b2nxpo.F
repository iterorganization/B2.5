*-----------------------------------------------------------------------
*.specification

      subroutine b2nxpo (nx, ny, ns, hx, pbs, 
     &  te, ne, fna, fch, respo, po) 
      use b2mod_types
      use b2mod_indirect
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_rates
      use b2mod_b2cmpb
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1), te(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), fna(-1:nx,-1:ny,0:1,0:ns-1), 
     &  fch(-1:nx,-1:ny,0:1), respo(-1:nx,-1:ny)
*   ..input/output arguments (changed on exit)
      real (kind=R8) ::
     *  po(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXPO solves the potential due to the boundary conditions exactly
*     when linearised solutions from B2STBC are not enough.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ib, iy, is, ireg
      real (kind=R8) ::
     1  vte, fchi, t0, t1, zarg, nef, tef, zaf
*   ..procedures
      intrinsic sqrt, log, max, abs
      real (kind=R8) ::
     *  bcint
      external subini, subend, xertst, bcint

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxpo')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')       

* -- SHEATH BOUNDARY CONDITIONS FOR POTENTIAL EQUATION --
* -- CASE WHERE THE BOUNDARY NAMELIST IS USED --      
      DO IB=1,NBC
        IF(BCCHAR(IB).EQ.'W') THEN
          t1 = 0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            t1 = max(t1,abs(respo(bc_list_x(iy,ib),bc_list_y(iy,ib))))
          ENDDO
          if(t1.gt.1.0_R8) then
          DO IY = 1, BC_LIST_SIZE(IB)
              write(*,*) ' Old potential, residual (W) : ',
     .	       bc_list_x(iy,ib),bc_list_y(iy,ib),
     .         po(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .         respo(bc_list_x(iy,ib),bc_list_y(iy,ib))
              vte = sqrt(bcint(te,nx,ny,iy,ib)/me)*
     .         abs(pbs(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .                 rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0))
              t0 = 0.0_R8
              do is = 0, ns-1
                t0 = t0-bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     .           fna(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .               rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0,is)
              enddo
              fchi = qe*t0
              zarg = sqrt(2.0_R8*pi)/(1.0_R8-gammae)*(fchi+
     &          fch(rightix(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     &              rightiy(bc_list_x(iy,ib),bc_list_y(iy,ib)),0))/
     &          (qe*bcint(ne,nx,ny,iy,ib)*vte)
              if(zarg.gt.1.0e-5_R8)
     >         po(bc_list_x(iy,ib),bc_list_y(iy,ib)) = potpar(ib,2) - 
     &          (bcint(te,nx,ny,iy,ib)/qe)*log(zarg)
              write(*,*) ' New potential : ', 
     .         po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            enddo
          endif
        ELSEIF(BCCHAR(IB).EQ.'E') THEN
          t1 = 0.0_R8
          DO IY = 1, BC_LIST_SIZE(IB)
            t1 = max(t1,abs(respo(bc_list_x(iy,ib),bc_list_y(iy,ib))))
          ENDDO
          if(t1.gt.1.0_R8) then
            DO IY = 1, BC_LIST_SIZE(IB)
              write(*,*) ' Old potential, residual (E) : ',
     .	       bc_list_x(iy,ib),bc_list_y(iy,ib),
     .         po(bc_list_x(iy,ib),bc_list_y(iy,ib)),
     .         respo(bc_list_x(iy,ib),bc_list_y(iy,ib))
              vte = sqrt(bcint(te,nx,ny,iy,ib)/me)*
     &         abs(pbs(bc_list_x(iy,ib),bc_list_y(iy,ib),0))
              t0 = 0.0_R8
              do is = 0, ns-1
                t0 = t0+bcint(rza(-1,-1,is),nx,ny,iy,ib)*
     &           fna(bc_list_x(iy,ib),bc_list_y(iy,ib),0,is)
              enddo
              fchi = qe*t0
              zarg = sqrt(2.0_R8*pi)/(1.0_R8-gammae)*
     &          (fchi-fch(bc_list_x(iy,ib),bc_list_y(iy,ib),0))/
     &          (qe*bcint(ne,nx,ny,iy,ib)*vte)
              if(zarg.gt.1.0e-5_R8)
     >         po(bc_list_x(iy,ib),bc_list_y(iy,ib)) = potpar(ib,2) -
     &          (bcint(te,nx,ny,iy,ib)/qe)*log(zarg)
              write(*,*) ' New potential : ', 
     .         po(bc_list_x(iy,ib),bc_list_y(iy,ib))
            enddo
          endif
        ENDIF
      ENDDO
      if(nbc.gt.0) goto 9000

*..default case
* ..compute conditions at west boundary
*   ..initialise ireg
      ireg = cbirwe
*   ..loop over boundary cells
      t1 = 0.0_R8
      do iy = -1, ny
        if(region(-1,iy,0).ne.0) t1 = max(t1,abs(respo(-1,iy)))
      enddo
      do iy = -1, ny
*    ..update ireg
    3   if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
          call xertst (ireg.lt.cbirwe+cbnrwe-1,
     &     'error in computing region on west boundary')
          ireg = ireg+1
          goto 3
        endif
        if(region(-1,iy,0).ne.0 .and. t1.gt.1.0_R8) then
          write(*,*) ' Old potential, residual (west) : ',
     .	   iy,po(-1,iy),respo(-1,iy)
*    ..compute cell-faced quantities
          tef = (te(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &           te(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &          (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
          nef = (ne(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &           ne(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &          (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
*    ..compute vte
          vte = sqrt(tef/me)*abs(pbs(topix(-1,iy),topiy(-1,iy),0))
*    ..compute Bohm fchi
          t0 = 0.0_R8
          do is = 0, ns-1
            zaf = (rza(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &             rza(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &            (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
            t0 = t0-zaf*
     &              fna(rightix(-1,iy),rightiy(-1,iy),0,is)
          enddo
          fchi = qe*t0
          zarg = sqrt(2.0_R8*pi)/(1.0_R8-gammae)*
     &     (fchi+fch(rightix(-1,iy),rightiy(-1,iy),0))/(qe*nef*vte)
          if(zarg.gt.1.0e-5_R8) po(-1,iy) = -(tef/qe)*log(zarg)
          write(*,*) ' New potential : ', po(-1,iy)
        endif
      enddo

* ..compute conditions at east boundary
*   ..initialise ireg
      ireg = cbirea
*   ..loop over boundary cells
      t1 = 0.0_R8
      do iy = -1, ny
        if(region(nx,iy,0).ne.0) t1 = max(t1,abs(respo(nx,iy)))
      enddo
      do iy = -1, ny
*    ..update ireg
    4   if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
          call xertst (ireg.lt.cbirea+cbnrea-1,
     &     'error in computing region on east boundary')
          ireg = ireg+1
          goto 4
        endif
        if(region(nx,iy,0).ne.0 .and. t1.gt.1.0_R8) then
          write(*,*) ' Old potential, residual (east) : ',
     .	   iy,po(nx,iy),respo(nx,iy)
*    ..compute cell-faced quantities
          tef = (te(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &           te(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &          (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
          nef = (ne(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &           ne(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &          (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
*    ..compute vte
          vte = sqrt(tef/me)*abs(pbs(nx,iy,0))
*    ..compute Bohm fchi
          t0 = 0.0_R8
          do is = 0, ns-1
            zaf = (rza(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &             rza(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &            (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
            t0 = t0+zaf*fna(nx,iy,0,is)
          enddo
          fchi = qe*t0
          zarg = sqrt(2.0_R8*pi)/(1.0_R8-gammae)*
     &     (fchi-fch(nx,iy,0))/(qe*nef*vte)
          if(zarg.gt.1.0e-5_R8) po(nx,iy) = -(tef/qe)*log(zarg)
          write(*,*) ' New potential : ', po(nx,iy)
        endif
      enddo

* ..return
 9000 continue
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxpo

      end subroutine b2nxpo

