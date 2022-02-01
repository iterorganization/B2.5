      module b2mod_aggs
      use b2mod_types
      use b2mod_b2cmgs
      use b2mod_constants
      use b2mod_subsys
      implicit none

      private

      public :: b2agsp, b2agsi

      integer, save :: istyle
      real (kind=R8), save :: xm, ym, w0y,
     *  axm, wid,
     *  sdiva0, sdiva1, sdivd0, sdivd1,
     *  dxb, dxc, dxm, dya, dyb, dym, dxa0, dxa1, dxd0, dxd1
      complex*16, save :: w0, w1, wm0, wm1
*     (Following data are meant to fit some Sonnet equilibrium for
*     Asdex-Upgrade.  The fit is only approximate.)
      real (kind=R8) :: cpsi, xmid, btmid, snl, ycore,
     *  xslo, yslo, xshi, yshi, sqz
      parameter (cpsi=0.41_R8, xmid=1.65_R8, btmid=2.00_R8,
     *  snl=0.94_R8, ycore=4.0e-1_R8, xslo=1.53_R8, yslo=-0.83_R8,
     *  xshi=1.65_R8, yshi=1.32_R8, sqz=0.0_R8)
      complex*16 :: ii, shift1, shift2
      parameter (ii=(0.0_R8,1.0_R8),
     *  shift1=(0.14_R8,0.0_R8), shift2=(0.0_R8,0.0_R8))
      logical, save :: init

*   ..initialisation
      data init/.false./

      contains
*-----------------------------------------------------------------------
*.specification

      subroutine b2agsp (ax, ay, kx, ky, crx0, cry0, psi0, fbz0)
      implicit none
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *   ax, ay
      integer kx, ky
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   crx0, cry0, psi0, fbz0
*   ..common block
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGSP evaluates the coordinates and the flux functions at one
*     single point on the mesh. This routine is application dependent,
*     and defines the specific geometry and magnetic field for a B2
*     calculation.
*
*
*  3. description (see also routine b2cdca)
*
*!!   The content of this routine is problem dependent.
*     This routine defines the geometry and the magnetic field for a
*     specific B2 calculation.
*     The B2 mesh is specified by way of a smooth function that maps
*     the unit square to the physical plane; the mesh in the physical
*     plane is the result of applying that function to a uniform mesh
*     on the unit square (boundary cells are handled differently, see
*     the calling routine). This routine implements that function.
*     The transformation may depend on certain numerical parameters
*     that have been passed to the subroutine b2agsi in the array
*     parg.
*
*     The principal input argument to b2agsp is the representation
*     (ax,ay) of a point in the unit square, or just outside the unit
*     square in the case of the boundary cells. (ax,ay) are the
*     coordinates that make the mesh appear uniform and square.
*     The transformation may have branch points, corresponding to an
*     x-point in the real geometry. The arguments (kx,ky) of b2agsp are
*     intended to resolve the ambiguity in the case that (ax,ay) lies
*     on a cut. If (ax,ay) does not lie on a cut, then (kx,ky) is
*     immaterial. If (ax,ay) does lie on a cut, then (ax,ay) is
*     interpreted as (ax+kx*eps,ay+ky*eps), where eps is a positive
*     infinitesimal. Legitimate values for kx and ky are -1, 0, and 1.
*     b2agsp returns the physical coordinates (crx0,cry0), and the
*     values of the poloidal flux function, psi0, and of the toroidal
*     field function, fbz0, at that point.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ax, ay - real, input.
*     (ax,ay) are the coordinates of a point inside or just outside the
*     unit square.
*
*     kx, ky - integer, input.
*     kx and ky resolve the ambiguity across the cut, as described in
*     Section 3. kx and ky must each have one of the values (-1,0,1).
*
*     crx0, cry0 - real, output.
*     crx0 and cry0 specify the horizontal and the vertical coordinate
*     of the point in the physical plane that corresponds to the input
*     arguments.
*
*     psi0 - real, output.
*     psi0 specifies the value of the poloidal flux function at the
*     point (crx0,cry0).
*     In the case of straight symmetry,
*       bx =  'dpsi/dy' ,
*       by = -'dpsi/dx' .
*     In the case of toroidal symmetry,
*       bx =  'dpsi/dy'/(2*pi*crx0) ,
*       by = -'dpsi/dx'/(2*pi*crx0) .
*     In the case of cylindrical symmetry,
*       bx =  'dpsi/dy'/(2*pi*cry0) ,
*       by = -'dpsi/dx'/(2*pi*cry0) .
*
*     fbz0 - real, output.
*     fbz0 specifies the value of the toroidal field function at the
*     point (crx0,cry0).
*     In the case of straight symmetry,
*       bz = fbz0 .
*     In the case of toroidal symmetry,
*       bz = fbz0/(2*pi*crx0) .
*     In the case of cylindrical symmetry,
*       bz = fbz0/(2*pi*cry0) .
*
*     (parameters for subroutine b2agsi)
*
*     nx, ny - integer, input.
*     nx and ny specify the number of interior cells along the first
*     and the second coordinate, respectively. The total number of
*     cells is (nx+2)*(ny+2); they are indexed by (-1:nx,-1:ny).
*     It will hold that 0.le.nx and 0.le.ny.
*
*     parg - (0:*) real*8 array, input
*     parg(0:*) contains a set of numerical parameters that specify
*     details of the geometry. The interpretation of these parameters
*     is problem-dependent; see the code.
*
*     isymm - integer, output.
*     isymm specifies the type of symmetry of the geometry: isymm.eq.0
*     implies a slab geometry, isymm.eq.1.or.isymm.eq.2 imply rotational 
*     symmetry about the crx=0 axis, and isymm.eq.3.or.isymm.eq.4 
*     indicate rotational symmetry about the cry=0 axis. 
*     Other values are not allowed. 
*     isymm.eq.1.or.isymm.eq.3 mean a toroidal geometry case, i.e. 
*     with a toroidal (out-of-plane) magnetic field component.
*     isymm.eq.2.or.isymm.eq.4 mean a cylindrical geometry case, i.e. 
*     no toroidal magnetic field component.
*
*     lblgm1 - character*60, output.
*     lblgm1 specifies a label to identify the geometry.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input/output.
*     /b2cmgs/ specifies the singular regions in the geometry.
*      nlrmx - integer constant.
*     nlrmx is an upper bound for nlreg, and is used in dimension
*     specifications.
*      nlreg - integer.
*     nlreg specifies the number of singular regions. It will hold that
*     0.le.nlreg.le.nlrmx.
*      nlxlo, nlxhi, nlylo, nlyhi - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.nlreg, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i) are
*     the bounds of the i-th singular region. Mesh locations (ix,iy)
*     in (nlxlo(i):nlxhi(i),nlylo(i):nlyhi(i)) make up the region, which
*     will be included in (-1:nx,-1:ny).
*     For nlreg.le.i.lt.nlrmx, nlxlo(i), nlxhi(i), nlylo(i) and nlyhi(i)
*     are unspecified.
*      nlloc - (0:nlrmx-1) integer array.
*     For 0.le.i.lt.4, nlloc(i) specifies the class of quantities for
*     which the i-th region is special, as follows:
*       nlloc(i) = 0: cell centers, indices (-1:nx,-1:ny).
*       nlloc(i) = 1: cell x-faces, indices (0:nx,-1:ny).
*       nlloc(i) = 2: cell y-faces, indices (-1:nx,0:ny).
*       nlloc(i) = 3: cell corners, indices (0:nx,0:ny).
*     For nlreg.le.i.lt.nlrmx, nlloc(i) is unspecified.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical iscut0, iscut1, iscut2
      integer type
      real (kind=R8) :: w0x
      complex*16 w2, w3, w4, w5, w6
*   ..procedures
      intrinsic real, cmplx, sqrt, exp, log, aimag
      logical esseq
      external xertst, xerrab, esseq, b2aggs

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     init - logical.
*     init is true iff the initialisation subroutine b2agsi has been called.
*
*     The values of shift2 and sqz are computed in order to get the
*     correct shape of the separatrix. The parameters shift1 and cpsi
*     determine the outward shift and the poloidal magnetic field.
*
*     Description of some statement functions:
*
*     fm2(t0,fa,fb) defines a weighted average of fa and fb.
*     fm3(t0;fa,dfa,fb) enjoys the properties:
*       fm3(0.0)=fa, dfm3/dt(0.0)=dfa, fm3(1.0)=fb
*     and fm3 is monotone if sgn(dfa).eq.sgn(fb-fa).
*     fm4(t0,ta,fa,dfa,tb,fb,dfb) enjoys the properties:
*       fm4(ta)=fa, dfm4/dt(ta)=dfa, fm4(tb)=fb, dfm4/dt(tb)=dfb.
*     Beware: it is required that fm4 be monotone as function of t0.
*     For this it is clearly necessary that dfa, dfb, and fb-fa be of
*     one sign, but that is not sufficient.
*
*     fracr(w0,s0) is a fractional linear transformation of w0; fracr
*     sends 0 to s0 and sends the unit circle to itself.
*     fraci(w0,s0) is a fractional linear transformation of w0; fraci
*     sends 0 to s0, ii to ii, and -ii to -ii.
*
*-----------------------------------------------------------------------
!.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agsp')
*   ..test init
      call xertst (init, 'b2agsi has not been called')
*   ..test ax, ay
      call xertst (-0.5_R8.le.ax.and.ax.le.1.5_R8 .and.
     &  -0.5_R8.le.ay.and.ay.le.1.5_R8, 'faulty argument ax or ay')
*   ..compute axm
      if (istyle.eq.0) then
       axm = ax
      else if (istyle.eq.1) then
       axm = ax/2.0_R8
      else if (istyle.eq.2) then
       axm = (1.0_R8+ax)/2.0_R8
      else if (istyle.eq.3) then
       axm = (0.5_R8+ax)/8.0_R8
      else if (istyle.eq.4) then
       axm = (1.5_R8+ax)/8.0_R8
      else if (istyle.eq.5) then
       axm = (2.5_R8+ax)/8.0_R8
      else if (istyle.eq.6) then
       axm = (3.5_R8+ax)/8.0_R8
      else if (istyle.eq.7) then
       axm = (4.5_R8+ax)/8.0_R8
      else if (istyle.eq.8) then
       axm = (5.5_R8+ax)/8.0_R8
      else if (istyle.eq.9) then
       axm = (6.5_R8+ax)/8.0_R8
      else if (istyle.eq.10) then
       axm = (7.5_R8+ax)/8.0_R8
      endif
*   ..set iscut0, iscut1, iscut2
      iscut0 = esseq(axm,0.25_R8).and.ay.lt.0.375_R8.and.
     1 .not.esseq(ay,0.375_R8)
      iscut1 = esseq(axm,0.5_R8)
      iscut2 = esseq(axm,0.75_R8).and.ay.lt.0.375_R8.and.
     1 .not.esseq(ay,0.375_R8)
*   ..test kx, ky
      call xertst (.not.(iscut0.or.iscut1.or.iscut2) .or.
     &  kx.eq.-1 .or. kx.eq.1, 'faulty argument kx on the cut')
      call xertst (-1.le.kx .and. kx.le.1, 'faulty argument kx')
      call xertst (-1.le.ky .and. ky.le.1, 'faulty argument ky')

* ..compute physical coordinates
*   ..note
*     (axm,ay) in the unit square must be mapped to a part of the
*     cross-section in the poloidal plane.
*     ay=0.375 is the separatrix.
*     ay=0.375 and axm=0.25: x-point from inner side.
*     ay=0.375 and axm=0.50: near top of plasma.
*     ay=0.375 and axm=0.75: x-point from outer side.
*     ay.lt.0.375 and 0.25.lt.axm.lt.0.75: main plasma region.
*     ay.lt.0.375 and axm.lt.0.25: inner divertor region below null.
*     ay.lt.0.375 and 0.75.lt.axm: outer divertor region below null.
*     0.375.lt.ay and 0.25.lt.axm.lt.0.75: main scrape-off.
*     0.375.lt.ay and axm.lt.0.25: inner divertor region.
*     0.375.lt.ay and 0.75.lt.axm: outer divertor region.
*     axm.eq.0.0: inner divertor target.
*     axm.eq.1.0: outer divertor target.
*     The parameters kx and ky serve to resolve ambiguities on the
*     cut. If (axm,ay) is not on the cut then (kx,ky) is immaterial.
*     If (axm,ay) is on the cut, then it is to be interpreted as
*     (axm+kx*eps,ay). Since the cut is along the y-coordinate, the
*     ambiguity pertains only to axm.
*   ..(axm,ay)->w0
*     The domain of w0 is a rectangle.
      if (axm.le.0.125_R8) then
       w0x = f0xa0(axm)
      else if (axm.le.0.25_R8) then
       w0x = f0xa1(axm)
      else if (axm.le.0.5_R8) then
       w0x = f0xb(axm)
      else if (axm.lt.0.75_R8) then
       w0x = f0xc(axm)
      else if (axm.lt.0.875_R8) then
       w0x = f0xd1(axm)
      else
       w0x = f0xd0(axm)
      endif
      if (ay.le.0.375_R8) then
       w0y = wid*f0ya(ay)
      else
       w0y = wid*f0yb(ay)
      endif
      w0 = cmplx(w0x,w0y,R8)
*   ..compute type
      if ((.not.iscut0.and.axm.lt.0.25_R8).or.
     1 (iscut0.and.kx.eq.-1)) then
       if (w0y.le.0.0_R8) then
        type = 2
       else
        type = 3
       endif
      else if ((.not.iscut2.and.0.75_R8.lt.axm).or.
     &   (iscut2.and.kx.eq.+1)) then
       if (w0y.le.0.0_R8) then
        type = 2
       else
        type = 1
       endif
      else if ((.not.iscut1.and.axm.lt.0.5_R8).or.
     &   (iscut1.and.kx.eq.-1)) then
       if (w0y.le.(4.0_R8*axm-1.0_R8)*ym) then
        type = 0
       else
        type = 3
       endif
      else if ((.not.iscut1.and.0.5_R8.lt.axm).or.
     &   (iscut1.and.kx.eq.+1)) then
       if (w0y.le.(3.0_R8-4.0_R8*axm)*ym) then
        type = 0
       else
        type = 1
       endif
      endif
*   ..w0->w1
*     the interval (0 to 4) is mapped to a circle.
      w1 = -ii*snl*(exp(-pi*ii*(w0+wm0)/2.0_R8)-1)-ii
*   ..w1->w2
*     a fractional linear map that preserves the unit circle.
*     the parameter shift1 controls the relative width of the inboard
*     and outboard scrape-off (in a way, the magnetic axis shift).
      w2 = fracr(w1,shift1)
*   ..w2->w3
*     the idea is simply w3=(2/pi)*arsinh(w2); this maps ii->ii and
*     -ii->-ii, and has branch cuts at those two points. the features
*     of the x-points are thereby introduced.
*     one must be careful with the cut for w2->w3.
      if (type.eq.0) then
*      (core plasma)
       w3 = (2.0_R8/pi)*log(w2+sqrt(1.0_R8+w2**2))
      else if (type.eq.1) then
*      (outer sol/div)
       w3 = (2.0_R8/pi)*(log(1+sqrt(1.0_R8+1.0_R8/w2**2))+log(w2))
      else if (type.eq.2) then
*      (lower private flux region)
       w3 =-(2.0_R8/pi)*(log(w2+sqrt(1.0_R8+w2**2))+pi*ii)
      else if (type.eq.3) then
*      (inner sol/div)
       w3 =-(2.0_R8/pi)*(log(1.0_R8+sqrt(1.0_R8+1.0_R8/w2**2))+log(-w2))
      endif
*   ..w3->w4
*     squeeze the radial coordinate a bit; the parameter sqz is found
*     empirically to give the correct elongation of the plasma.
      w4 = (w3-sqz/w3)/(1.0_R8+sqz)
*   ..w4->w5
*     another fractional linear transformation, preserving the points
*     ii and -ii. the parameter shift2 is found empirically to give
*     the proper outward shift.
      w5 = fraci(w4,shift2)
*   ..w5->w6
*     linear map that sends ii to the upper saddle point and -ii to the
*     lower saddle point.
      w6 = cmplx(xshi-xslo,yshi-yslo,R8)*w5/(2.0_R8*ii)+
     &     cmplx(xslo+xshi,yslo+yshi,R8)/2.0_R8
*   ..poloidal coordinates
      crx0 = real(w6)
      cry0 = aimag(w6)
*   ..magnetic field
      psi0 = gpsi(w0y)
      fbz0 = gfbz(w0y)

* ..return
      call subend ()
      return
      end subroutine b2agsp

*-----------------------------------------------------------------------
*.end b2agsp

      subroutine b2agsi (nx, ny, parg, isymm, lblgm1, ostyle)
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     *   parg(0:*)
*   ..output arguments (unspecified on entry)
      integer isymm, ostyle
      character lblgm1*60

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2agsi')
*   ..set isymm
      isymm = 1
      call ipgeti('b2agsi_isymm',isymm)
      call xertst (0.le.isymm .and. isymm.le.4,'faulty argument isymm')
*   ..set parameters that depend on parg(0)
      if (nint(parg(0)).eq.parg(0)) then
       lblgm1 = 'ASDEX-UPGRADE(1.65m,2.0T)SN, 1997 CONFORMAL MAP'
       istyle = nint(parg(0))
       ostyle = istyle
       dxa0 = parg(1)
       dxa1 = parg(2)
       dxb = parg(3)
       dxm = parg(4)
       dxc = parg(5)
       dxd1 = parg(6)
       dxd0 = parg(7)
       dya = parg(8)
       dym = parg(9)
       dyb = parg(10)
       wid = parg(11)
       sdiva0 = parg(12)
       sdiva1 = parg(13)
       sdivd0 = parg(14)
       sdivd1 = parg(15)
      else
       call xerrab ('faulty argument parg(0)')
      endif
*   ..invert w2=-ii to find wm0
      w1 = (-ii-shift1)/(1.0_R8+conjg(shift1)*ii)
      wm0 = (2.0_R8/pi)*ii*log(1.0_R8+(w1+ii)*ii/snl)
*   ..invert w2=ii to find wm1, xm, ym
      w1 = (ii-shift1)/(1.0_R8-conjg(shift1)*ii)
      wm1 = (2.0_R8/pi)*ii*log(1.0_R8+(w1+ii)*ii/snl)
      xm = real(wm1-wm0)
      ym = aimag(wm1-wm0)
*   ..live with a Sun f77 (1992, SC1.0) compiler error
*     (!! using 'f77 -fast' or 'f77 -O' on the Sparc it turns out that,
*     without the following no-op call to dummyz, the variable wm0 is
*     mistaken for a scalar temporary and is not stored, even though it
*     occurs in a save statement. The first call to dummyz ensures that
*     wm0 will be stored. The second call is for reasons of symmetry.)
      call dummyz (wm0)
      call dummyz (wm1)
*   ..initialise /b2cmgs/
      nlreg = 0
*   ..introduce cuts
      if (istyle.eq.0) then
       call b2aggs (nx/4, nx/4, -1, 3*ny/8-1, 1)
       call b2aggs (nx/4, nx/4, 0, 3*ny/8, 3)
       call b2aggs (3*nx/4, 3*nx/4, -1, 3*ny/8-1, 1)
       call b2aggs (3*nx/4, 3*nx/4, 0, 3*ny/8, 3)
      else if (istyle.eq.1) then
       call b2aggs (3*nx/4, 3*nx/4, -1, 3*ny/8-1, 1)
       call b2aggs (3*nx/4, 3*nx/4, 0, 3*ny/8, 3)
      else if (istyle.eq.2) then
       call b2aggs (nx/4, nx/4, -1, 3*ny/8-1, 1)
       call b2aggs (nx/4, nx/4, 0, 3*ny/8, 3)
      endif
*   ..set init
      init = .true.

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agsi

      end subroutine b2agsi

      function fm2(t0,fa,fb)
      implicit none
      real (kind=R8) :: fm2
      real (kind=R8) :: t0, fa, fb
      fm2 = (fa*(1.0_R8-t0)**4+fb*t0**4)/(t0**4+(1.0_R8-t0)**4)
      return
      end function fm2

      function fm3(t0,fa,dfa,fb)
      implicit none
      real (kind=R8) :: fm3
      real (kind=R8) :: t0, fa, dfa, fb
      fm3 = fa + (fb-fa)*dfa*t0/(dfa*t0+(fb-fa)*(1.0_R8-t0))
      return
      end function fm3

      function fm4(t0,ta,fa,dfa,tb,fb,dfb)
      implicit none
      real (kind=R8) :: fm4
      real (kind=R8) :: t0, ta, tb, fa, fb, dfa, dfb
      fm4 = fm2 ((t0-ta)/(tb-ta),
     &  fm3((t0-ta)/(tb-ta),fa,dfa*(tb-ta),fb),
     &  fm3((tb-t0)/(tb-ta),fb,-dfb*(tb-ta),fa))
      return
      end function fm4

      complex*16 function fracr(w0,s0)
      implicit none
      complex*16 w0, s0
      intrinsic conjg
      fracr = (w0+s0)/(conjg(s0)*w0+1.0_R8)
      return
      end function fracr

      complex*16 function fraci(w0,s0)
      implicit none
      complex*16 w0, s0
      fraci = (w0+s0)/(1.0_R8-s0*w0)
      return
      end function fraci

      function f0xa0(t0)
      implicit none
      real (kind=R8) :: f0xa0, t0
      f0xa0 = fm4(t0, 0.0_R8, -sdiva0, dxa0, 0.125_R8, -sdiva1, dxa1)
      return
      end function f0xa0

      function f0xa1(t0)
      implicit none
      real (kind=R8) :: f0xa1, t0
      f0xa1 = fm4(t0, 0.125_R8, -sdiva1, dxa1, 0.25_R8, 0.0_R8, dxb)
      return
      end function f0xa1

      function f0xb(t0)
      implicit none
      real (kind=R8) :: f0xb, t0
      f0xb = fm4(t0, 0.25_R8, 0.0_R8, dxb, 0.5_R8, xm, dxm)
      return
      end function f0xb

      function f0xc(t0)
      implicit none
      real (kind=R8) :: f0xc, t0
      f0xc = fm4(t0, 0.5_R8, xm, dxm, 0.75_R8, 4.0_R8, dxc)
      return
      end function f0xc

      function f0xd0(t0)
      implicit none
      real (kind=R8) :: f0xd0, t0
      f0xd0 = fm4(t0, 0.875_R8, 4.0_R8+sdivd1, dxd1,
     .                  1.0_R8, 4.0_R8+sdivd0, dxd0)
      return
      end function f0xd0

      function f0xd1(t0)
      implicit none
      real (kind=R8) :: f0xd1, t0
      f0xd1 =
     . fm4(t0, 0.75_R8, 4.0_R8, dxc, 0.875_R8, 4.0_R8+sdivd1, dxd1)
      return
      end function f0xd1

      function f0ya(t0)
      implicit none
      real (kind=R8) :: f0ya, t0
      f0ya = fm4(t0, 0.0_R8, -ycore, dya, 0.375_R8, 0.0_R8, dym)
      return
      end function f0ya

      function f0yb(t0)
      implicit none
      real (kind=R8) :: f0yb, t0
      f0yb = fm4(t0, 0.375_R8, 0.0_R8, dym, 1.0_R8, 1.0_R8, dyb)
      return
      end function f0yb

      function gfbz(w0y)
      implicit none
      real (kind=R8) :: gfbz, w0y
      gfbz = 2.0_R8*pi*btmid*xmid
      return
      end function gfbz

      function gpsi(w0y)
      implicit none
      real (kind=R8) :: gpsi, w0y
      gpsi = 2.0_R8*pi*cpsi*w0y
      return
      end function gpsi

      subroutine dummyz (z)
      use b2mod_types
      implicit none
      complex*16 z
*     ------------------------------------------------------------------
*     This subroutine does nothing.
*     ------------------------------------------------------------------
      return
      end subroutine dummyz

      end module b2mod_aggs
