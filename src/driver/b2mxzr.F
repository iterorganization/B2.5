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

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2mxzr')
*   ..test nCv, ns
      call xertst (0.lt.nCv, 'faulty argument nCv')
      call xertst (1.le.ns, 'faulty argument ns')

*   ..zero-out parts of (na) residual
      do is = 0, ns-1
       do iCv = mpg%nCi+1, nCv
        resco(iCv,is) = 0.0_R8
       enddo
      enddo
*   ..zero-out parts of (ua) residual
      do is = 0, ns-1
       do iCv = mpg%nCi+1, nCv
        resmo(iCv,is) = 0.0_R8
       enddo
      enddo
       do iCv = mpg%nCi+1, nCv
        resmt(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (te) residual
       do iCv = mpg%nCi+1, nCv
        reshe(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (ti) residual
       do iCv = mpg%nCi+1, nCv
        reshi(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (tn) residual
       do iCv = mpg%nCi+1, nCv
        reshn(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (po) residual
       do iCv = mpg%nCi+1, nCv
        respo(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (kt) residual
       do iCv = mpg%nCi+1, nCv
        reskt(iCv) = 0.0_R8
       enddo
*   ..zero-out parts of (zt) residual
       do iCv = mpg%nCi+1, nCv
        reszt(iCv) = 0.0_R8
       enddo

#ifdef WG_TODO
      do iCv = 1, mpg%nCi
*   ..zero out resignore regions
        if(resignore(iCv,1).eq.0) then
          do is = 0, ns-1
            resco(iCv,is) = 0.0_R8
            resmo(iCv,is) = 0.0_R8
          enddo
          resmt(iCv) = 0.0_R8
          reshe(iCv) = 0.0_R8
          reshi(iCv) = 0.0_R8
          reshn(iCv) = 0.0_R8
          respo(iCv) = 0.0_R8
          reskt(iCv) = 0.0_R8
          reszt(iCv) = 0.0_R8
        endif
      enddo
#endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2mxzr

      end subroutine b2mxzr
