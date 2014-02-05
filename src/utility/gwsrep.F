      subroutine gwsrep (iwks, ityp)
      use b2mod_types
      implicit none
      integer iwks, ityp
*     ------------------------------------------------------------------
*     GWSREP defines the workstation representation for GKS. It should
*     normally be called just once for each workstation, immediately
*     after opening the workstation.
*     This routine calls GKS routines to define the representations
*     of polylines, polymarkers, fill area, text, and color on the
*     workstation identified by (iwks). The content of this routine
*     may depend on the environment. The variable (ityp) may be used
*     to distinguish different types of workstation.
*     ------------------------------------------------------------------
*     The present version of gwsrep ignores ityp, and only specifies
*     a color representation. The representations of polylines,
*     polymarkers, fill area, and text are left as by default.
*!!   Note 0: Parameter stylbw defines the assignment to the colors 0
*     and 1 (black-white or white-black).
*!!   Note 1: parameters ncd0 and ncd1 define the resolution of the
*     color map; they may be changed if that is desired. ncd1 might be
*     set to a small multiple of ncd0, if they are equal then the 'fine'
*     selection of colors is the same as the 'coarse' selection. I have
*     found that ncd0+ncd1 must be at most 20 for NCAR on the Cray at
*     NERSC, else we get an error message "too many colors allocated"
*     from the ctrans/ictrans viewing utility.
*     ------------------------------------------------------------------
*   ..local variables
      integer iasf(13), stylbw, ncd0, ncd1, i
      parameter (stylbw=2, ncd0=5, ncd1=5)
      save iasf
      real (kind=R8) ::
     *   t0
*   ..procedures
      intrinsic real, sqrt
      external gsasf, gscr_d
*   ..initialisation
      data iasf/13*1/
*     ------------------------------------------------------------------
*   ..set aspect source flags to 'individual'
*     (as recommended in NCARG3 manual)
      call gsasf (iasf)
*   ..set colors 0 and 1
      if (stylbw.eq.0) then
*    ..leave colors 0 and 1 as per default
*      (this is either 0=black, 1=white or 0=white, 1=black.)
      else if (stylbw.eq.1) then
*    ..set colors 0=black and 1=white
       call gscr_d (iwks, 0, 0.0_R8, 0.0_R8, 0.0_R8)
       call gscr_d (iwks, 1, 1.0_R8, 1.0_R8, 1.0_R8)
      else if (stylbw.eq.2) then
*    ..set colors 0=white and 1=black
       call gscr_d (iwks, 0, 1.0_R8, 1.0_R8, 1.0_R8)
       call gscr_d (iwks, 1, 0.0_R8, 0.0_R8, 0.0_R8)
      endif
*   ..set colors 2..2+6*ncd0-1 to a coarse selection
      do i = 0, ncd0-1
       t0 = real(i)/ncd0
       call gscr_d (iwks, 2+0*ncd0+i, 1.0_R8, sqrt(t0), 0.0_R8)
       call gscr_d (iwks, 2+1*ncd0+i, sqrt(1-t0), 1.0_R8, 0.0_R8)
       call gscr_d (iwks, 2+2*ncd0+i, 0.0_R8, 1.0_R8, sqrt(t0))
       call gscr_d (iwks, 2+3*ncd0+i, 0.0_R8, sqrt(1-t0), 1.0_R8)
       call gscr_d (iwks, 2+4*ncd0+i, sqrt(t0), 0.0_R8, 1.0_R8)
       call gscr_d (iwks, 2+5*ncd0+i, 1.0_R8, 0.0_R8, sqrt(1-t0))
      enddo
*   ..set colors 2+6*ncd0..2+6*(ncd0+ncd1)-1 to a fine selection
      do i = 0, ncd1-1
       t0 = real(i)/ncd1
       call gscr_d (iwks, 2+6*ncd0+0*ncd1+i, 1.0_R8, sqrt(t0), 0.0_R8)
       call gscr_d (iwks, 2+6*ncd0+1*ncd1+i, sqrt(1-t0), 1.0_R8, 0.0_R8)
       call gscr_d (iwks, 2+6*ncd0+2*ncd1+i, 0.0_R8, 1.0_R8, sqrt(t0))
       call gscr_d (iwks, 2+6*ncd0+3*ncd1+i, 0.0_R8, sqrt(1-t0), 1.0_R8)
       call gscr_d (iwks, 2+6*ncd0+4*ncd1+i, sqrt(t0), 0.0_R8, 1.0_R8)
       call gscr_d (iwks, 2+6*ncd0+5*ncd1+i, 1.0_R8, 0.0_R8, sqrt(1-t0))
      enddo
      return
*     ------------------------------------------------------------------
      end subroutine gwsrep
