*-----------------------------------------------------------------------
*.specification

      subroutine b2mxzr (nCv, ns, mpg, resco, resmo, resmt, reshe,
     &  reshi, reshn, respo, reskt, reszt)
      use b2mod_types
      use b2us_map
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, ns
      type (mapping), intent (in) :: mpg
*   ..input/output arguments
      real (kind=R8) ::
     *  resco(nCv,0:ns-1), resmo(nCv,0:ns-1),
     &  resmt(nCv), reshe(nCv), reshi(nCv),
     &  respo(nCv), reshn(nCv), reskt(nCv), reszt(nCv)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2MXZR sets the residuals to zero along the boundaries of the
*     mesh and on the special regions.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine sets the residuals to zero along the boundaries
*     and on the special regions. The residuals on the boundary cells
*     and on the special regions often do not have any physical
*     meaning, and it is desirable to eliminate them from the
*     computation of the norms and from the output.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, iCv
*   ..procedures
      external xertst, sfill
      external b2xvfx, b2xxgs
*   ..initialisation

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxzr')
*   ..test nCv, ns
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..zero-out parts of (na) residual
      do is = 0, ns-1
!WG_TODO       call b2xxgs (nx, ny, 0.0_R8, resco(-1,-1,is), 0)
       do iCv = mpg%nCi+1, nCv
        resco(iCv,is) = 0.0_R8
       enddo
      enddo
*   ..zero-out parts of (ua) residual
      do is = 0, ns-1
!WG_TODO       call b2xxgs (nx, ny, 0.0_R8, resmo(-1,-1,is), 0)
       do iCv = mpg%nCi+1, nCv
        resmo(iCv,is) = 0.0_R8
       enddo
      enddo
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, resmt, 0)
       do iCv = mpg%nCi+1, nCv
        resmt(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (te) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, reshe, 0)
       do iCv = mpg%nCi+1, nCv
        reshe(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (ti) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, reshi, 0)
       do iCv = mpg%nCi+1, nCv
        reshi(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (tn) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, reshn, 0)
       do iCv = mpg%nCi+1, nCv
        reshn(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (po) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, respo, 0)
       do iCv = mpg%nCi+1, nCv
        respo(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (kt) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, reskt, 0)
       do iCv = mpg%nCi+1, nCv
        reskt(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (zt) residual
!WG_TODO      call b2xxgs (nx, ny, 0.0_R8, reszt, 0)
       do iCv = mpg%nCi+1, nCv
        reszt(iCv) = 0.0_R8
       enddo

#ifdef WG_TODO
      do ix = 0, nx-1
        do iy = 0, ny-1
*   ..zero out interior boundaries and resignore regions
          if(resignore(ix,iy,1).eq.0 .or. leftix(ix,iy).eq.-2 .or.
     &       rightix(ix,iy).eq.nx+1) then
            do is = 0, ns-1
              resco(ix,iy,is) = 0.0_R8
              resmo(ix,iy,is) = 0.0_R8
            enddo
            resmt(ix,iy) = 0.0_R8
            reshe(ix,iy) = 0.0_R8
            reshi(ix,iy) = 0.0_R8
            reshn(ix,iy) = 0.0_R8
            respo(ix,iy) = 0.0_R8
            reskt(ix,iy) = 0.0_R8
            reszt(ix,iy) = 0.0_R8
          endif
        enddo
      enddo
#endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxzr

      end subroutine b2mxzr
