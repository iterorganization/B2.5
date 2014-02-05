*-----------------------------------------------------------------------
*.specification

      subroutine b2wgyp (n0, n1, n0d, fun, scal, fa0, fa1, ncnta,
     &  fb0, fb1, ncntb)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer n0, n1, n0d, ncnta, ncntb
      real (kind=R8) ::
     *   fun(0:n0d-1,0:n1-1), scal, fa0, fa1, fb0, fb1
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGYP performs the core NCAR operations to produce a contour
*     plot.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine performs the core NCAR operations needed to produce
*     a contour plot.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     n0, n1 - integer, input.
*     n0 and n1 specify the dimensions of the function along the first
*     and the second coordinate direction.
*     It is required that 1.le.n0 and 1.le.n1.
*
*     n0d - integer, input.
*     n0d specifies the declared leading dimension of the array fun.
*     It is required that n0.le.n0d.
*
*     fun - (0:n0d-1,0:n1-1) real*8 array, input.
*     fun specifies the function of which a contour plot is to be made.
*     Only the section fun(0:n0-1,0:n1-1) will be used.
*
*     scal - real, input.
*     scal specifies a scaling factor for the function fun. The
*     contours will be labelled with the value of fun/scal.
*     It is required that 0.lt.scal.
*
*     fa0, fa1 - real, input.
*     fa0 and fa1 specify the extreme values of fun for the solid
*     contour region.
*
*     ncnta - integer, input.
*     ncnta specifies the number of solid contour intervals.
*     If 0.lt.ncnta a solid contour will be plotted at the values
*       fun = ((ncnta-i)*fa0+i*fa1)/ncnta, for 0.le.i.le.ncnta.
*     It is required that 0.le.ncnta.
*
*     fb0, fb1 - real, input.
*     fb0 and fb1 specify the extreme values of fun for the dashed
*     contour region.
*
*     ncntb - integer, input.
*     ncntb specifies the number of dashed contour intervals.
*     If 0.lt.ncntb a dashed contour will be plotted at the values
*       fun = ((ncntb-i)*fb0+i*fb1)/ncntb, for 0.le.i.le.ncntb.
*     It is required that 0.le.ncntb.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer niaia, nacra, niama, nrwrk, niwrk
      parameter (niaia=10)
      logical lgerr
      integer i, ic, iftp, irdum
      integer iaia(0:niaia-1), igia(0:niaia-1)
      character cllbs*9
      real (kind=R8) ::
     *   cont
      real (kind=R4), allocatable :: rwrk(:), acrx(:), acry(:)
      integer, allocatable :: iwrk(:), iama(:)
      real (kind=R4) fun4(0:n0d-1,0:n1-1)
*   ..procedures
      integer nerro
      external subini, subend, xertst, cprect, cpseti,
     &  cpsetr_d, cpsetc, cpclam, cplbam, cplbdr, cpdrpl, cpcldm, 
     &  cpcldr,
     &  arinam, arpram, arscam, gncolr, lbfilq, entsr, nerro, errof

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgyp')
*   ..test n0, n1, n0d
      call xertst (1.le.n0.and.1.le.n1.and.n0.le.n0d,
     &  'faulty argument n0, n1, n0d')
*   ..test scal
      call xertst (0.lt.scal, 'faulty argument scal')
*   ..test ncnta, ncntb
      call xertst (0.le.ncnta.and.0.le.ncntb,
     &  'faulty argument ncnta, ncntb')
*   ..set dimension parameters
*!!   Note: The following assignments may be changed.
*     We do not know very well how to estimate the amount of workspace
*     needed by the AREAS and CONPACK sub-packages of NCAR graphics.
*     It is to be expected that the required workspace increases with
*     the size of the mesh and with the complexity of the contour plots.
      nrwrk = 50000+50*n0*n1
      niwrk = 50000+50*n0*n1
      nacra = 10000+50*n0*n1
      niama = 100000+100*n0*n1
*   ..obtain real*8 work space
      allocate(rwrk(1:nrwrk))
      allocate(acrx(1:nacra))
      allocate(acry(1:nacra))
*   ..obtain integer work space
      allocate(iwrk(1:niwrk))
      allocate(iama(1:niama))

* ..make plot
      fun4=fun
      call cpsetr_d ('SFS', scal)
      call cprect (fun4, n0d, n0, n1, rwrk, nrwrk, iwrk, niwrk)
*   ..primary divisions
      call cpseti ('NCL', ncnta+1)
      ic = 0
      do i = 0, ncnta
       cont = ((ncnta-i)*fa0+i*fa1)/ncnta
       write (cllbs,'(1pe9.2)') cont/scal
       ic = ic+1
       call cpseti ('PAI', ic)
       call cpseti ('AIA', 2+i)
       call cpseti ('AIB', 1+i)
       call cpsetr_d ('CLV', cont)
       call cpseti ('CLU', 3)
       call cpseti ('CLD', 65535)
       call cpsetc ('LLT', cllbs)
      enddo
*   ..inquire area fill style
      call lbfilq (iftp)
*   ..produce the area map
*     (This computation can run out of space. It is attempted to trap
*     that error and continue with a useful computation.)
      call arinam (iama, niama)
      call entsr (irdum, 1)
      lgerr = .false.
      if (iftp.ne.0.and..not.lgerr) then
*    ..add contour lines to the area map
       call cpclam (fun4, rwrk, iwrk, iama)
       if (nerro(irdum).ne.0) then
        lgerr = .true.
        call errof ()
       endif
      endif
      if (.not.lgerr) then
*    ..add label boxes to the area map
       call cplbam (fun4, rwrk, iwrk, iama)
       if (nerro(irdum).ne.0) then
        lgerr = .true.
        call errof ()
       endif
      endif
      if (.not.lgerr) then
*    ..pre-process the area map
       call arpram (iama, 0, 0, 0)
       if (nerro(irdum).ne.0) then
        lgerr = .true.
        call errof ()
       endif
      endif
      if (iftp.ne.0.and..not.lgerr) then
*    ..scan the area map and fill the areas
       call arscam (iama, acrx, acry, nacra, iaia, igia,
     &   niaia, gncolr)
       if (nerro(irdum).ne.0) then
        lgerr = .true.
        call errof ()
       endif
      endif
      call entsr (irdum, 2)
*   ..produce output
      call cplbdr (fun4, rwrk, iwrk)
      if (.not.lgerr) then
       call cpcldm (fun4, rwrk, iwrk, iama,
     &   cpdrpl)
      else
       call cpcldr (fun4, rwrk, iwrk)
      endif
*     (Note: at one time we also produced a labelbar.  This required
*     messy code, and I now prefer simply to have the labels on the
*     contours.)
*   ..secondary divisions
      if (0.lt.ncntb) then
       call cpseti ('NCL', ncntb+1)
       ic = 0
       do i = 0, ncntb
        cont = ((ncntb-i)*fb0+i*fb1)/ncntb
        ic = ic+1
        call cpseti ('PAI', ic)
        call cpseti ('AIA', 0)
        call cpseti ('AIB', 0)
        call cpsetr_d ('CLV', cont)
        call cpseti ('CLU', 1)
*$$$    call cpseti ('CLD - dash 0101010101010101', 21845)
        call cpseti ('CLD - dash 0000111100001111', 3855)
       enddo
       if (.not.lgerr) then
        call cpcldm (fun4, rwrk, iwrk, iama, cpdrpl)
       else
        call cpcldr (fun4, rwrk, iwrk)
       endif
      endif

* ..return
      deallocate(rwrk,acrx,acry,iwrk,iama)
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgyp

      end subroutine b2wgyp
