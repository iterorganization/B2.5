*-----------------------------------------------------------------------
*.specification

      subroutine b2yqdr (ninp, nout, nscr)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1), nscr(0:0)
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YQDR is a post-processor for b2. It provides a quick display of
*     the evolution data.
*
*
*  3. description (see also routine b2cdca)
*
*     This subroutine provides a quick display of the evolution of
*     a b2 calculation.
*
*     ...
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:1) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides output control.
*     ninp(1): un*formatted; provides the evolution data.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): graphical; displays the evolution data.
*
*     nscr - (0:0) integer array, input.
*     nscr specifies the scratch unit numbers.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nfldmx, nfld, ns
      parameter (nfldmx=2000)
      integer linlog(0:nfldmx-1), idum(0:9), k, npts, ipts, indx, indy,
     &  k0, k1, style
      real (kind=R8) ::
     *   z0(0:nfldmx-1), z1(0:nfldmx-1), z(0:nfldmx-1),
     &  u0(0:nfldmx-1), u1(0:nfldmx-1), zprev(0:nfldmx-1),
     &  cutlo, x0, x1, y0, y1, t0, t1, t2, osc
      character lblmn*120, chtit(0:nfldmx-1)*60,
     &  chlab(0:nfldmx-1)*12
      real (kind=R8), allocatable :: wx(:),wy(:)
*   ..procedures
      intrinsic min, max, abs, nint, log10
      logical cfendf
      real (kind=R8) ::
     *   scal10, machsfr
      external subini, subend, xertst, xerrab, cfruin, cfruch,
     &  cfrure, cfendf, scal10, machsfr, gacwk, gdawk, pagini, pagend,
     &  gpm_d

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     ns - integer.
*     ns specifies the number of atomic species in the calculation.
*     The species are indexed by (0:ns-1).
*     It will hold that 1.le.ns.
*
*     nfld - integer.
*     nfld specifies the number of data fields specified for each
*     iteration or timestep.
*
*     lblmn - character*120.
*     lblmn identifies the instance of the program employed for the main
*     calculation.
*
*     ...
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2yqdr')
*   ..test ninp, nout, nscr
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &  1.le.nout(0).and.1.le.nout(1).and.1.le.nscr(0),
     &  'faulty argument ninp, nout, nscr')
*   ..set cutlo to the machine safe range
      cutlo = machsfr()

* ..first pass through data
*   ..read header information
      call cfruin (ninp(1), 1, idum, 'ns')
      ns = idum(0)
      call xertst (0.le.ns, 'faulty input ns')
      call cfruch (ninp(1), 120, lblmn, 'label')
      call cfruin (ninp(1), 1, idum, 'nfld')
      nfld = idum(0)
      call xertst (0.le.nfld, 'faulty input nfld')
      do k = 0, nfld-1
       call cfruin (ninp(1), 1, idum(0), 'linlog')
       linlog(k) = idum(0)
       call xertst (linlog(k).eq.0.or.linlog(k).eq.1,
     &   'faulty input linlog')
       call cfruch (ninp(1), 60, chtit(k), 'title')
       call cfruch (ninp(1), 12, chlab(k), 'label')
      enddo
*   ..read and store data
      call cfrure (ninp(1), nfld, z, 'data')
      write (nscr(0)) (z(k),k=0,nfld-1)
      do k = 0, nfld-1
       z0(k) = z(k)
       z1(k) = z(k)
       u0(k) = abs(z(k))
       u1(k) = 0
       zprev(k) = z(k)
      enddo
      npts = 1
    1 if (.not.cfendf(ninp(1))) then
       call cfrure (ninp(1), nfld, z, 'data')
       write (nscr(0)) (z(k),k=0,nfld-1)
       do k = 0, nfld-1
        z0(k) = min(z0(k),z(k))
        z1(k) = max(z1(k),z(k))
        u0(k) = u0(k)+abs(z(k))
        u1(k) = u1(k)+abs(z(k)-zprev(k))
        zprev(k) = z(k)
       enddo
       npts = npts+1
       goto 1
      endif
      endfile (nscr(0))
*   ..echo basic data
      write (nout(0),
     &  '(2x,a,i3/2(2x,a,i6/)/t7,a,t21,a,t33,a,t45,a,t53,a)')
     &  'number of atomic species:', ns,
     &  'number of data fields:', nfld,
     &  'number of data points:', npts,
     &  'label', 'min', 'max', 'osc', 'title'
      do k = 0, nfld-1
       if (u0(k).eq.0) then
        osc = 0
       else
        osc = u1(k)/u0(k)
       endif
       write (nout(0),'(1x,i3,2x,a12,1p,2e12.4,2x,0pf6.1,2x,a60)')
     &   k, chlab(k), z0(k), z1(k), 1000*osc, chtit(k)
      enddo
*   ..test signs
      do k = 0, nfld-1
       if (linlog(k).eq.1) then
        call xertst (0.le.z0(k), 'faulty data for logarithmic scale')
       endif
      enddo
*   ..allocate workspace wx and wy
      allocate(wx(0:npts))    ! dpc maybe only npts-1
      allocate(wy(0:npts))    ! dpc maybe only npts-1

* ..make plots
*   ..activate graphics on nout(1)
      call gacwk (nout(1))
*   ..start new plot
    2 continue
       read (ninp(0),*,end=3,err=3) indx, indy
       call xertst (0.le.indx.and.indx.lt.nfld, 'faulty input indx')
       call xertst (0.le.indy.and.indy.lt.nfld, 'faulty input indy')
       x0 = z0(indx)
       x1 = z1(indx)
       if (linlog(indy).eq.0) then
        t0 = scal10(z1(indy)-z0(indy))
        t1 = scal10(max(abs(z0(indy)),abs(z1(indy))))
        if (t1.eq.0) then
         t2 = 1
        else
         t2 = max(t0,1.0e-6_R8*t1)
        endif
        y0 = nint((z0(indy)+z1(indy))/(2*t2))*t2-t2
        y1 = y0+2*t2
       else
        k0 = nint(log10(z0(indy)+cutlo)-0.5_R8)
        k1 = nint(log10(z1(indy)+cutlo)+0.5_R8)
        y0 = 10.0_R8**max(k0,k1-12)
        y1 = 10.0_R8**k1
       endif
*    ..start new page
       style = linlog(indy)
       call pagini (style, x0, x1, y0, y1, chlab(indx), chlab(indy),
     &   chtit(indy))
*    ..read data points
       rewind (nscr(0))
       do ipts = 0, npts-1
        read (nscr(0)) (z(k),k=0,nfld-1)
        if (linlog(indx).eq.0) then
         wx(ipts) = z(indx)
        else
         wx(ipts) = log10(z(indx)+cutlo)
        endif
        if (linlog(indy).eq.0) then
         wy(ipts) = z(indy)
        else
         wy(ipts) = log10(z(indy)+cutlo)
        endif
       enddo
*    ..plot data points
       call gpm_d (npts, wx, wy)
*    ..close page and loop
       call pagend ()
       goto 2
    3  continue
*   ..deactivate graphics on nout(1)
      call gdawk (nout(1))

* ..return
      deallocate(wx,wy)
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2yqdr

      end
