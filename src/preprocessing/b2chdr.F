*-----------------------------------------------------------------------
*.specification

      subroutine b2chdr (ninp, nout, ns)
      use b2mod_diag
      use b2mod_ipmain
      use b2mod_subsys
      use b2mod_switches
      use b2mod_version
      use b2mod_types
      use b2us_data
      use b2us_geo
      use b2us_io
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:2), ns
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CHDR computes and writes to external files the chords which
*     mimic the path of ECH rays, taking into account their bouncing on
*     vessel walls. It also provides the location of the resonance layer
*     and the control volumes crossed by the rays, in which the ECH
*     power will be deposited, and the respective fraction of the power
*     they will receive.
*
*
*  2. specification
*
*     None.
*
*
*  3. description (see also routine b2cdca)
*
*     ...
*
*  4. references
*
*     ...
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:2) integer array, input.
*     ninp specifies the input units as follows.
*     ninp(0): formatted; provides chords parameters.
*     ninp(1): un*formatted; provides the geometry.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides printed output.
*     nout(1): formatted; provides the chords and their labels.
*     nout(2): un*formatted; provides the location of the resonance
*              layer, the list of control volumes and the relative
*              fraction of deposited power. !!TODO
*
*
*  6. error indicators
*
*     In case an error condition is detected, a call is made to the
*     routine xerrab. This causes an error message to be printed,
*     after which the program halts.
*
*
*  7. accuracy
*
*     ...
*
*
*  8. further comments
*
*     ...
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer idum(0:9)
      character lblgm*120, lblcp*120, lblmn*120, lblrc*120
      character cnamip*80, cvalip*80, line*500, temp*120
*   ..procedures
      logical is_comment, cfendf
      integer :: multirays
      real (kind=R8) :: multiwidth, multitheta
      integer :: iFc, numcoord, fileunit
      real :: x1, y1, x2, y2, f1, f2, t, bex
      real, dimension(10000) :: xcoords, ycoords
      logical :: flag
      external xertst, xerrab, is_comment
      external cfruch, cfruin, ipgeti, cfendf
      external b2xxid, b2rflb, strip_spaces

*   ..subprogram start-up calls
      call subini ('b2chdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.1.le.nout(0).and.
     &             1.le.nout(1).and.1.le.nout(2),
     &  'faulty argument ninp, nout')

* ..computation
*   ..obtain lblmn
      call b2xxid ('b2ch', lblmn(1:60))
      call b2rflb (ninp(0), lblmn(61:120))
      label = lblmn
*   ..read code internal parameters
    1 continue
       read (ninp(0),'(a)',end=2,err=91) line
       if (.not.is_comment(line)) then
         read (line,*) cnamip, cvalip
         call ipsetc (cnamip, cvalip)
       endif
       goto 1
    2 continue
      write (nout(0),'(2x,a)') '(end of list of internal parameters)'
*   ..read geometry
      call cfruch (ninp(1), 120, lblgm, 'label')
      call cfruin (ninp(1), 1, idum, 'isymm')
      isymm = idum(0)
      call read_b2fgmtry (ninp(1), mpg, geo)
      call init_geometry (mpg, geo, switch)
*   ..internal parameters
      call ipgeti('b2chdr_multiple_rays', multirays)
      call ipgetr('b2chdr_multiple_width', multiwidth)
      call ipgetr('b2chdr_multiple_angle', multitheta)

      ! open the output file
      open(unit=10, file='outtest.dat', status='replace',
     &     action='write', iostat=fileunit)
      if (fileunit /= 0) then
         print *, "error: unable to open file outtest.dat"
         stop
      endif

      numcoord = 0
      bex = 1.54

      do iFc = 1, mpg%nFc-1
        ! Loop through each edge of the cell (assuming a quadrilateral with 4 edges)
        ! Edge 1: between vertex 1 and vertex 2 (bottom edge)
        f1 = geo%fcBb(iFc,3)  ! Field value at left face (or first face)
        f2 = geo%fcBb(iFc+1,3)     ! Field value at right face (or second face)

        ! Check if contour level crosses this edge
        if (min(f1, f2) <= bex .and. max(f1, f2) >= bex) then
          numcoord = numcoord + 1
          t = (bex - f1) / (f2 - f1)  ! Linear interpolation parameter
          xcoords(numcoord) = geo%vxX(mpg%fcVx(iFc,1)) + t *
     &            (geo%vxX(mpg%fcVx(iFc+1,1)) -geo%vxX(mpg%fcVx(iFc,1)))
          ycoords(numcoord) = geo%vxY(mpg%fcVx(iFc,1)) + t *
     &            (geo%vxY(mpg%fcVx(iFc+1,1)) -geo%vxY(mpg%fcVx(iFc,1)))
          write(10, '(f10.6,1x,f10.6)') xcoords(numcoord),
     &                                  ycoords(numcoord)
        end if
      enddo

    ! Output the number of contour points found
      print *, "Number of contour points found: ", numcoord

      ! close the file
      close(unit=10)

      print *, "contour coordinates saved to outtest.dat"

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read internal parameters')

*-----------------------------------------------------------------------
*.end b2chdr

      end subroutine b2chdr
