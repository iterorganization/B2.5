*-----------------------------------------------------------------------
*.specification

      subroutine b2agdr (ninp, nout,g,m)
      use b2mod_types
      use b2mod_version
      use b2mod_b2cmfs
      use b2mod_b2cmgs
      use b2ag_parameters
      use b2us_geo
      use b2us_map
      use b2us_io
      implicit none
*   ..input arguments (unchanged on exit)
      integer ninp(0:1), nout(0:1) 
      type(geometry) :: g
      type(mapping) :: m
*   ..output arguments (unspecified on entry)
*     (none)

*-----------------------------------------------------------------------
!WP  b2ag_us
*.documentation 
*
*  1. purpose
*
*     B2AGDR computes and writes to an external file the 
*     specification of a b2 geometry and magnetic field.
*
*
*  3. description (see also routine b2cdca)
*
*     This subroutine computes and writes out a specification of a b2
*     geometry and magnetic field. Summary information will be written
*     to unit nout(0), which must be connected to a formatted file.
*     A complete specification of the geometry and the magnetic field
*     will be written to unit nout(1), which must be connected to an
*     un*formatted file.
*
*     B2 employs a two-dimensional, curvilinear, topologically
*     rectangular mesh. The coordinates of the mesh are labelled 'x'
*     and 'y'. 
*
*     The primary data describing the mesh are the coordinates of the
*     four vertices of all the cells of the mesh. These data are stored
*     in the arrays vxX and vxY. The metric parameters (cvVvol, fcQalf,
*     fcQbet, fcQgam; described elsewhere) are derived from vxX and vxY.
*
*     The primary data describing the magnetic field are the flux
*     functions vxFpsi and vxFfbz at the cell vertices. The magnetic 
*     field components at cell centers are derived from vxFpsi and vxFfbz.
*
*     Singular regions in the geometry are specified by the data in the
*     commmon block /b2cmgs/.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     ninp - (0:1) integer array, input.
*     ninp specifies the input unit numbers, as follows.
*     ninp(0): formatted; provides geometric parameters.
*     ninp(1): formatted; provides mhd equilibrium data.
*
*     nout - (0:1) integer array, input.
*     nout specifies the output unit numbers, as follows.
*     nout(0): formatted; provides summary printed output.
*     nout(1): un*formatted; full specification of the geometry and
*     field.
*
*     g - type, input.
*     g specifies the geometry type as described module b2us_geo
*
*     m - type, input.
*     m specifies the mapping type as described module b2us_map
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..common blocks
*   ..local variables
      integer idum(0:9), ix, iy, redef_gmtry
      integer nxiso1, nxiso2, nyiso1, nyiso2, artificial_slab
      character lblgm*120, idout*80, id*8, cnamip*80, cvalip*80
      character line*500
      integer istyle
*   ..initialisation
      data nxiso1/-2/, nxiso2/-2/, nyiso1/-2/, nyiso2/-2/, 
     . artificial_slab/0/
*   ..procedures
      logical streql, is_comment
      external xertst, xerrab, streql, cfwuin, cfwuch, cfwure
      external ipgeti, intfaceh1, intfacer1, intfacev1, is_comment
      external b2xxid, b2agsi, b2agfs, b2agcf, b2agmt, b2agmx, b2agbb,
     &  b2xgbs1, b2wfgm, b2wugm

*-----------------------------------------------------------------------
*.documentation-internal
*
*     The workspace array w is used as a stack for storage of dynamic
*     real*8 arrays. An array 'name' in w(0:*) will be stored starting at
*     location w(l'name').
*
*     The following common blocks have their outermost declaration in
*     this routine; they need not be preserved between calls.
*
*     /b2cmfs/ is reserved for numerical mhd equilibrium data.
*     (Presently empty--to be used in conjunction with grid generation
*     based on a numerical mhd equilibrium.)
*
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
*     Description of some local variables (including arrays in w(0:*):
*
*     crx, cry - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), crx(ix,iy,0:3) specify the
*     horizontal coordinates and cry(ix,iy,0:3) specify the vertical
*     coordinates (in the physical plane) of the four corners of the
*     (ix,iy) cell, mapped as follows:
*     (,,0): lower left corner,  (,,1): lower right corner,
*     (,,2): upper left corner,  (,,3): upper right corner.
*     One expects the following identities to hold away from the
*     boundaries:
*     crx(ix,iy,0)=crx(ix-1,iy,1)=crx(ix,iy-1,2)=crx(ix-1,iy-1,3),
*     cry(ix,iy,0)=cry(ix-1,iy,1)=cry(ix,iy-1,2)=cry(ix-1,iy-1,3),
*     However, there may be cuts in the geometry, in which case these
*     identities need not hold.
*
*     fpsi - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), fpsi(ix,iy,i) specifies the
*     potential for the magnetic field components (bx,by) at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bx=dpsi/dy and by=-dpsi/dx.
*     In the case of toroidal or cylindrical symmetry,
*       bx=(dpsi/dy)/(2*pi*R) and by=-(dpsi/dx)/(2*pi*R),
*     where R is the distance to the axis of symmetry.)
*
*     ffbz - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy,i) in (-1:nx,-1:ny,0:3), ffbz(ix,iy,i) specifies the
*     field function for the magnetic field component bz at the
*     (ix,iy,i) position.
*     (In the case of straight symmetry, bz=ffbz. In the case of
*     toroidal or cylindrical symmetry, bz=ffbz/(2*pi*R), where R is
*     the distance to the axis of symmetry.)
*
*     bb - (-1:nx,-1:ny,0:3) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), bb(ix,iy,0:3) specifies the
*     magnetic field at the center of the (ix,iy) cell.
*     bb(ix,iy,0:2) are the (x,y,z)-components and bb(ix,iy,3) is
*     the absolute magnetic field strength.
*     (bb(,,3)=sqrt(bb(,,0)**2+bb(,,1)**2+bb(,,2)**2.)
*     It will hold that 0.lt.bb(,,3).
*
*     vol - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), vol(ix,iy) specifies the volume of
*     the (ix,iy) cell. (In the case of toroidal symmetry this volume
*     includes a factor 2*pi*crx, while in the case of cylindrical
*     symmetry it includes a factor 2*pi*cry.)
*     It will hold that 0.lt.vol(,).
*
*     hx, hy - (-1:nx,-1:ny) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), hx(ix,iy) and hy(ix,iy) specify the
*     x-diameter and the y-diameter, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.hx(,) and 0.lt.hy(,).
*
*     cell_width, cell_height - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), cell_width(ix,iy) and cell_height(ix,iy)
*     specify the average width and height, respectively, of the (ix,iy) cell.
*     It will hold that 0.lt.cell_width(,) and 0.lt.cell_height(,).
*
*     qz - (-1:nx,-1:ny,0:1) real*8 array.
*     For (ix,iy) in (-1:nx,-1:ny), qz(ix,iy,0) specifies sin(t) and
*     qz(ix,iy,1) specifies cos(t), where t is the complementary angle
*     between the x-direction and the y-direction on the (ix,iy) cell.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. (qz(ix,iy,0) is the inner product of the x and y unit
*     vectors.)
*     It will hold that -1.lt.qz(,,0).lt.1 and 0.lt.qz(,,1).le.1.
*
*     qc - (-1:nx,-1:ny) real*8 array, output.
*     For (ix,iy) in (-1:nx,-1:ny), qc(ix,iy) specifies cos(t), where 
*     t is the complementary angle between the x-direction of the 
*     (ix,iy) cell and its left face.
*     t is in the range (-pi/2 to pi/2) and vanishes for an orthogonal
*     geometry. 
*     It will hold that 0.lt.qc(,).le.1.
*
*     gs - (-1:nx,-1:ny,0:2) real*8 array.
*     (Let gsx(,)=gs(,,0), gsy(,)=gs(,,1) and gsz(,)=gs(,,2).)
*     For (ix,iy) in (0:nx,-1:ny), gsx(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its left neighbor. For ix.eq.-1,
*     gsx(ix,-1:ny) will hold 0.
*     For (ix,iy) in (-1:nx,0:ny), gsy(ix,iy) specifies the area of the
*     face between the (ix,iy) cell and its bottom neighbor. For iy.eq.-1,
*     gsy(-1:nx,iy) will hold 0.
*     For (ix,iy) in (-1:nx,-1:ny), gsz(ix,iy) specifies the area of the
*     (ix,iy) cell normal to the third (ignorable) coordinate.
*
*     pbs - (-1:nx,-1:ny,0:1) real*8 array.
*     (Let pbsx(,)=pbs(,,0) and pbsy(,)=pbs(,,1).)
*     For (ix,iy) in (-1:nx,-1:ny), pbsx(ix,iy) specifies the product
*     (bx/bb)*sx and pbsy(ix,iy) specifies (by/bb)*sy, or (magnetic 
*     field pitch)*(surface area), namely the parallel contact area, 
*     on the (ix,iy) cell left and bottom faces respectively.  
*     See routine b2xgbs for the form of interpolation employed.
*
*     crx1, cry1, fpsi1, ffbz1, vol1, 
*     hx1, hy1, qz1, qc1, gs1, pbs1 - real*8 arrays.
*     These arrays characterise the geometry on the auxiliary mesh.
*
*     lblgm - character*120.
*     lblgm identifies the instance of the program employed to prepare
*     the geometry.
*
*     idout - character*80.
*     idout specifies the desired selection of output segments for
*     printed and graphical output. Associated with each segment is a
*     four-character word, which is identified on the printed output.
*     (The same word may be associated with several segments.)
*     A particular segment will be produced only if its associated word
*     has a match, as defined by the routine strmas, in idout.
*
*     isymm - integer.
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
*     parg - (0:99) real*8 array.
*     parg(0:*) contains a set of numerical parameters that specify
*     details of the geometry. The interpretation of these parameters
*     is problem-dependent; see the code.
*
*     periodic_bc - integer
*     periodic_bc tells whether the cuts correspond to periodic boundary
*     conditions (for example in the case of limited plasmas or
*     a diverted stellarator island in a toroidally averaged geometry). 
*     periodic_bc.eq.0 is the default.
*     periodic_bc.eq.1 indicates that periodic boundary conditions are
*     being used.
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agdr')
*   ..test ninp, nout
      call xertst (1.le.ninp(0).and.1.le.ninp(1).and.
     &             1.le.nout(0).and.1.le.nout(1),
     &  'faulty argument ninp, nout')
*   ..test nCv, nFc, nVx
      call xertst (1.le.m%nCv.and.1.le.m%nFc.and.1.le.m%nVx,
     & 'faulty argument nCv, nFc, nVx')
*   ..set output control
      idout = '....'

* ..initialise the geometry
      redef_gmtry = 1 ! Since version 3.0.6, we apply the pbsx correction by default
      call ipgeti ('b2agdr_redef_gmtry', redef_gmtry)
      call ipgeti ('b2agdr_artificial_slab', artificial_slab)
      call xertst (artificial_slab.eq.0 .or. artificial_slab.eq.1, 
     & 'faulty argument artificial_slab' )
*   ..use parg, obtain isymm and lblgm
       call b2xxid ('b2ag', lblgm(1:60))
       call b2agsi (parg, isymm, lblgm(61:120), istyle)
*   ..obtain grid and mhd equilibrium data from file
      call b2agfs(ninp(1),istyle,lblgm,g,m,artificial_slab)

*   ..echo lblgm
      write (nout(0),'(/2x,a)') 'label for geometry'
      write (nout(0),'(/4x,a/4x,a)') lblgm(1:60), lblgm(61:120)


*   ..compute metric coefficients on auxiliary mesh
       call b2agmt (g, m, isymm, artificial_slab)

*   ..for artificial_slab, correct the connectivity for the cuts
       if (artificial_slab.eq.1) then
        !set_connectivity_us
        call set_connectivity_us(g,m)
       endif

      write(*,*) 'b2agdr1: g%cvBb(m%nCv,0) ',g%cvBb(m%nCv,0)
      write(*,*) 'b2agdr2: g%cvBb(m%nCv,0) ',g%cvBb(m%nCv,0)
      write(*,'()')
      if (g%vxFfBz(m%nVx).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has POSITIVE toroidal magnetic field.'
        if (isymm.ne.0) write(*,'(a)')
     .   'This means that the Bx(grad B) direction points DOWN.'
        write(*,'(a)')
     .   'Bz is out of the page towards you.'
        write(nout(0),'(/2x,a)')
     .   'This case has POSITIVE toroidal magnetic field.'
        if (isymm.ne.0) write(nout(0),'(2x,a)')
     .   'This means that the Bx(grad B) direction points DOWN.'
        write(nout(0),'(2x,a)')
     .   'Bz is out of the page towards you.'
      else if (g%vxFfBz(m%nVx).lt.0.0_R8) then
         write(*,'(a)')
     .   'This case has NEGATIVE toroidal magnetic field.'
        if (isymm.ne.0) write(*,'(a)')
     .   'This means that the Bx(grad B) direction points UP.'
        write(*,'(a)')
     .   'Bz is into the page away from you.'
        write(nout(0),'(/2x,a)')
     .   'This case has NEGATIVE toroidal magnetic field.'
        if (isymm.ne.0) write(nout(0),'(2x,a)')
     .   'This means that the Bx(grad B) direction points UP.'
        write(nout(0),'(2x,a)')
     .   'Bz is into the page away from you.'
      else
        write(*,'(a)')
     .   'This case has ZERO toroidal magnetic field.'
        write(nout(0),'(/2x,a)')
     .   'This case has ZERO toroidal magnetic field.'
      endif
      if (g%cvBb(m%nCv,2).ne.0.0_R8) then
        write(*,'(a)')
     .   'If this is not correct, please set :'
        write(*,'(a)')
     .   '"b2agfs_Bt_reversal" to "1" in "b2ag.dat" and rerun.'
        write(*,'(a)')
     .   'If already set to "1", set to "0" instead.'
      endif
      write(*,'()')
      if (g%cvBb(m%nCv,0).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has POSITIVE poloidal magnetic field.'
        write(*,'(a)')
     .   'This means the poloidal field points in the '//
     .   'direction of INCREASING poloidal index ix.'
        write(*,'(a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'INCREASING poloidal index ix.'
        write(*,'(a)')
     .   'This is the SAME convention as for poloidal '//
     .   'quantities.'
        write(nout(0),'(2x,a)')
     .   'This case has POSITIVE poloidal magnetic field.'
        write(nout(0),'(2x,a)')
     .   'This means the poloidal field points in the '//
     .   'direction of INCREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'INCREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'This is the SAME convention as for poloidal '//
     .   'quantities.'
      else
        write(*,'(a)')
     .   'This case has NEGATIVE poloidal magnetic field.'
        write(*,'(a)')
     .   'This means the poloidal field points in the '//
     .   'direction of DECREASING poloidal index ix.'
        write(*,'(a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'DECREASING poloidal index ix.'
        write(*,'(a)')
     .   'This is the OPPOSITE convention as for '//
     .   'poloidal quantities.'
        write(nout(0),'(2x,a)')
     .   'This case has NEGATIVE poloidal magnetic field.'
        write(nout(0),'(2x,a)')
     .   'This means the poloidal field points in the '//
     .   'direction of DECREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'The parallel velocities ua and ue will be '//
     .   'positive when pointing in the direction of '//
     .   'DECREASING poloidal index ix.'
        write(nout(0),'(2x,a)')
     .   'This is the OPPOSITE convention as for '//
     .   'poloidal quantities.'
      end if
      if (g%cvBb(m%nCv,0)*g%cvBb(m%nCv,2).lt.0.0_R8) then
        write(*,'(a)')
     .   'This case has CO-directed plasma current and '//
     .   'toroidal field.'
        write(nout(0),'(2x,a)')
     .   'This case has CO-directed plasma current and '//
     .   'toroidal field.'
      else if (g%cvBb(m%nCv,0)*g%cvBb(m%nCv,2).gt.0.0_R8) then
        write(*,'(a)')
     .   'This case has COUNTER-directed plasma current and '//
     .   'toroidal field.'
        write(nout(0),'(2x,a)')
     .   'This case has COUNTER-directed plasma current and '//
     .   'toroidal field.'
      endif
      write(*,'(a)')
     . 'If this is not correct, please set :'
       write(*,'(a)')
     . '"b2agfs_pit_rescale" to "-1.0" in "b2ag.dat" and rerun.'
      write(*,'(a)')
     . 'If already set to "-1.0", set to "1.0" instead.'
      write(*,'()')

       !Take absolute values of magnetic field to comply with structured b2ag
       g%cvBb = abs(g%cvBb)
       g%fcBb = abs(g%fcBb)
       g%vxBb = abs(g%vxBb)

       !Check all field to not contain Nan
       !mapping
       if (m%nCi .ne. m%nCi) then
          call xerrab('Nan in nCi')
       endif   
       if (m%nCg .ne. m%nCg) then
          call xerrab('Nan in nCg')
       endif
       if (m%nCv .ne. m%nCv) then
          call xerrab('Nan in nCv')
       endif 
       if (m%nFc .ne. m%nFc) then
          call xerrab('Nan in nFc')
       endif   
       if (m%nVx .ne. m%nVx) then
          call xerrab('Nan in nVx')
       endif
       if (m%nFt .ne. m%nFt) then
          call xerrab('Nan in nFt')
       endif                        
       if (m%nFs .ne. m%nFs) then
          call xerrab('Nan in nFs')
       endif

       if (any(m%cvFcP .ne. m%cvFcP)) then
          call xerrab('Nan in cvFcP')
       endif   
       if (any(m%cvFc .ne. m%cvFc)) then
          call xerrab('Nan in cvFc')
       endif
       if (any(m%fcCv .ne. m%fcCv)) then
          call xerrab('Nan in fcCv')
       endif 
       if (any(m%fcVx .ne. m%fcVx)) then
          call xerrab('Nan in fcVx')
       endif   
       if (any(m%cvVxP .ne. m%cvVxP)) then
          call xerrab('Nan in cvVxP')
       endif
       if (any(m%cvVx .ne. m%cvVx)) then
          call xerrab('Nan in cvVx')
       endif                        
       if (any(m%vxFcP .ne. m%vxFcP)) then
          call xerrab('Nan in vxFcP')
       endif

       if (any(m%vxFc .ne. m%vxFc)) then
          call xerrab('Nan in vxFc')
       endif   
       if (any(m%vxCvP .ne. m%vxCvP)) then
          call xerrab('Nan in vxCvP')
       endif
       if (any(m%vxCv .ne. m%vxCv)) then
          call xerrab('Nan in vxCv')
       endif 
       if (any(m%ftCvP .ne. m%ftCvP)) then
          call xerrab('Nan in ftCvP')
       endif   
       if (any(m%ftCv .ne. m%ftCv)) then
          call xerrab('Nan in ftCv')
       endif
       if (any(m%ftFcP .ne. m%ftFcP)) then
          call xerrab('Nan in ftFcP')
       endif                        
       if (any(m%ftFc .ne. m%ftFc)) then
          call xerrab('Nan in ftFc')
       endif
       if (any(m%cvFt .ne. m%cvFt)) then
          call xerrab('Nan in cvFt')
       endif
       if (any(m%fsFcP .ne. m%fsFcP)) then
          call xerrab('Nan in fsFcP')
       endif                        
       if (any(m%fsFc .ne. m%fsFc)) then
          call xerrab('Nan in fsFc')
       endif

       if (any(m%fcReg .ne. m%fcReg)) then
          call xerrab('Nan in fcReg')
       endif
       if (any(m%cvReg .ne. m%cvReg)) then
          call xerrab('Nan in cvReg')
       endif                        
       if (any(m%ftReg .ne. m%ftReg)) then
          call xerrab('Nan in ftReg')
       endif
       if (any(m%intcellP .ne. m%intcellP)) then
          call xerrab('Nan in intcellP')
       endif
       if (any(m%intcellR .ne. m%intcellR)) then
          call xerrab('Nan in intcellR')
       endif       
       if (any(m%fcLbl .ne. m%fcLbl)) then
          call xerrab('Nan in fcLbl')
       endif                        
       if (any(m%cvLbl .ne. m%cvLbl)) then
          call xerrab('Nan in cvLbl')
       endif       
       if (any(m%ftLbl .ne. m%ftLbl)) then
          call xerrab('Nan in ftLbl')
       endif 
       !geometry
       if (any(g%cvBb .ne. g%cvBb)) then
          call xerrab('Nan in cvBb')
       endif 
       if (any(g%cvEb .ne. g%cvEb)) then
          call xerrab('Nan in cvEb')
       endif 
       if (any(g%cvX .ne. g%cvX)) then
          call xerrab('Nan in cvX')
       endif 
       if (any(g%cvY .ne. g%cvY)) then
          call xerrab('Nan in cvY')
       endif 
       if (any(g%cvSz .ne. g%cvSz)) then
          call xerrab('Nan in cvSz')
       endif 
       if (any(g%cvHz .ne. g%cvHz)) then
          call xerrab('Nan in cvHz')
       endif 
       if (any(g%cvQgam .ne. g%cvQgam)) then
          call xerrab('Nan in cvQgam')
       endif 
       if (any(g%cvVol .ne. g%cvVol)) then
          call xerrab('Nan in cvVol')
       endif
       if (any(g%fcBb .ne. g%fcBb)) then
          call xerrab('Nan in fcBb')
       endif  
       if (any(g%fcS .ne. g%fcS)) then
          call xerrab('Nan in fcS')
       endif   
       if (any(g%fcHc .ne. g%fcHc)) then
          call xerrab('Nan in fcHc')
       endif  
       if (any(g%fcHt .ne. g%fcHt)) then
          call xerrab('Nan in fcHt')
       endif  
       if (any(g%fcQgam .ne. g%fcQgam)) then
          call xerrab('Nan in fcQgam')
       endif  
       if (any(g%fcQalf .ne. g%fcQalf)) then
          call xerrab('Nan in fcQalf')
       endif   
       if (any(g%fcQbet .ne. g%fcQbet)) then
          call xerrab('Nan in fcQbet')
       endif  
       if (any(g%fcPbs .ne. g%fcPbs)) then
          call xerrab('Nan in fcPbs')
       endif
       if (any(g%vxBb .ne. g%vxBb)) then
          call xerrab('Nan in vxBb')
       endif       
       if (any(g%vxX .ne. g%vxX)) then
          call xerrab('Nan in vxX')
       endif
       if (any(g%vxY .ne. g%vxY)) then
          call xerrab('Nan in vxY')
       endif
       if (any(g%vxFfbz .ne. g%vxFfbz)) then
          call xerrab('Nan in vxFfbz')
       endif       
       if (any(g%vxFpsi .ne. g%vxFpsi)) then
          call xerrab('Nan in vxFpsi')
       endif
       if (any(g%cvConn .ne. g%cvConn)) then
          call xerrab('Nan in cvConn')
       endif  
       if (any(g%fsPsi .ne. g%fsPsi)) then
          call xerrab('Nan in fsPsi')
       endif     

* ..output - see module b2us_io and b2us.F
*   ..write headers for output files - code from b2us.F
      call cfverw (nout(1), newversion)
      idum(0) = m%nCi
      idum(1) = m%nCg
      idum(2) = m%nCv
      idum(3) = m%nFc
      idum(4) = m%nVx
      idum(5) = m%nFs
      idum(6) = m%nFt
      call cfwuin (nout(1),7,idum,'nCi,nCg,nCv,nFc,nVx,nFs,nFt')
      idum(0) = m%nCmxVx
      idum(1) = m%nCmxFc
      idum(2) = m%nVmxCv
      idum(3) = m%nVmxFc
      idum(4) = m%nCmxNv
      call cfwuin (nout(1),5,idum,'nCmxVx,nCmxFc,nVmxCv,nVmxFc,nCmxNv')
      idum(0) = m%isClassicalGrid
      call cfwuin (nout(1),1,idum,'isClassicalGrid')
      if (m%isClassicalGrid == 1) then
        idum(0) = m%nx
        idum(1) = m%ny
        idum(2) = m%nncut
      else
        ! dummy values for unstructured grid
        idum(0) = -1
        idum(1) = -1
        idum(2) = -1
      endif
      call cfwuin (nout(1), 3, idum, 'nx,ny,nncut')
      call cfwuch (nout(1), 120, lblgm, 'label')
      idum(0) = isymm
      call cfwuin (nout(1), 1, idum, 'isymm')
      
*   ..write summary of geometry
        call write_b2fgmtry (nout(1), m, g)

      call write_f30

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.error conditions

   91 call xerrab ('error trying to read heading *param')
   92 call xerrab ('error trying to read parg')
   93 call xerrab
     1 ('error trying to read non-default internal parameters')

*-----------------------------------------------------------------------
*.end b2agdr

      end subroutine b2agdr

      subroutine set_connectivity_us(g,m)
      use b2mod_types
      use b2us_geo
      use b2us_map
      implicit none

      type(geometry) :: g 
      type(mapping)  :: m
      integer ncvl, ncvr, nfpl, nfpr, nfrl, nfrr, nvl1, nvl2,
     &  nvr1, nvr2, cvl(m%nCv), cvr(m%nCv),
     & fpl(m%nFc), fpr(m%nFc), frl(m%nFc), frr(m%nFc),
     & vl1(m%nVx), vl2(m%nVx), vr1(m%nVx), vr2(m%nVx), 
     & iFc, cvs(1:2), i, iCgl, iCgr
      intrinsic pack

      !Identify cells, faces and vertices along the cut
      ncvl = 0
      ncvr = 0
      nfpl = 0
      nfpr = 0
      nfrl = 0
      nfrr = 0
      nvl1 = 0
      nvl2 = 0
      nvr1 = 0
      nvr2 = 0

      do iFc = 1,m%nFc
          cvs = m%fcCv(iFc,:)
          if (m%cvReg(cvs(1)).ne.m%cvReg(cvs(2))) then
            if ((m%cvReg(cvs(1)).eq.1) .and.
     &          (m%cvReg(cvs(2)).eq.3)) then
            !left cut
            cvl(ncvl+1) = cvs(2)
            ncvl = ncvl + 1
            fpl(nfpl+1) = iFc
            nfpl = nfpl + 1
            vl1(nvl1+1) = m%cvVx(m%cvVxP(cvs(2),1)+1) 
            vl2(nvl2+1) = m%cvVx(m%cvVxP(cvs(2),1)+2)
            nvl1 = nvl1 +1
            nvl2 = nvl2 +1
            frl(nfrl+1) = m%cvFc(m%cvFcP(cvs(2),1)) 
            nfrl = nfrl + 1
            elseif ((m%cvReg(cvs(2)).eq.1) .and.
     &            (m%cvReg(cvs(1)).eq.3)) then
            !left cut
            cvl(ncvl+1) = cvs(1)
            ncvl = ncvl + 1
            fpl(nfpl+1) = iFc
            nfpl = nfpl + 1
            vl1(nvl1+1) = m%cvVx(m%cvVxP(cvs(1),1)+1)
            vl2(nvl2+1) = m%cvVx(m%cvVxP(cvs(1),1)+2)
            nvl1 = nvl1 +1
            nvl2 = nvl2 +1
            frl(nfrl+1) = m%cvFc(m%cvFcP(cvs(1),1))
            nfrl = nfrl + 1
            elseif ((m%cvReg(cvs(1)).eq.1) .and.
     &            (m%cvReg(cvs(2)).eq.4)) then 
            !right cut
            cvr(ncvr+1) = cvs(1)
            ncvr = ncvr + 1
            fpr(nfpr+1) = iFc
            nfpr = nfpr + 1
            vr1(nvr1+1) = m%cvVx(m%cvVxP(cvs(1),1)+1)
            vr2(nvr2+1) = m%cvVx(m%cvVxP(cvs(1),1)+2)
            nvr1 = nvr1 +1
            nvr2 = nvr2 +1
            frr(nfrr+1) = m%cvFc(m%cvFcP(cvs(1),1))
            nfrr = nfrr + 1

            elseif ((m%cvReg(cvs(1)).eq.4) .and.
     &            (m%cvReg(cvs(2)).eq.1)) then
            !right cut
            cvr(ncvr+1) = cvs(2)
            fpr(nfpr+1) = iFc
            nfpr = nfpr + 1
            vr1(nvr1+1) = m%cvVx(m%cvVxP(cvs(2),1)+1)
            vr2(nvr2+1) = m%cvVx(m%cvVxP(cvs(2),1)+2)
            nvr1 = nvr1 +1
            nvr2 = nvr2 +1
            frr(nfrr+1) = m%cvFc(m%cvFcP(cvs(2),1)) 
            nfrr = nfrr + 1   
            endif
          endif
      enddo

        !Cut the arrays
        cvl = cvl(1:ncvl)
        cvr = cvr(1:ncvr)
        fpl = fpl(1:nfpl)
        fpr = fpr(1:nfpr)
        frl = frl(1:nfrl)
        frr = frr(1:nfrr)
        vl1  = vl1(1:nvl1)
        vl2  = vl2(1:nvl2)
        vr1  = vr1(1:nvr1)
        vr2  = vr2(1:nvr2)

      !Change the connectivity
      !connect the left and right information via psi values (more general of
      !psi is accurate enough) or based on the original structured ordening
      !Here based on original structured ordening
      do i = 1,ncvl  
        !Cell left of left cut
        m%cvVx(m%cvVxP(cvl(i),1)+1) = vr1(i)
        m%cvVx(m%cvVxP(cvl(i),1)+2) = vr2(i)
        m%cvFc(m%cvFcP(cvl(i),1)+1) = fpr(i)

        !Poloidal face at left cut
        !m%fcCv(fpl(i),1) = cvr(i)

        !Radial face left of left cut
        ! replace vl1 with vr1
        if (m%fcVx(frl(i),1).eq.vl1(i)) then
           m%fcVx(frl(i),1) = vr1(i)
        elseif (m%fcVx(frl(i),2).eq.vl1(i)) then
           m%fcVx(frl(i),2) = vr1(i)
        else
          call xerrab('set_connectivity_us: verts of radial face')           
        endif

        !Cell left of right cut
        m%cvVx(m%cvVxP(cvr(i),1)+1) = vl1(i)
        m%cvVx(m%cvVxP(cvr(i),1)+2) = vl2(i)
        m%cvFc(m%cvFcP(cvr(i),1)+1) = fpl(i)

        !Poloidal face at right cut
        !m%fcCv(fpr(i),1) = cvl(i)

        !Radial face left of right cut
        if (m%fcVx(frr(i),1) .eq. vr1(i)) then
           m%fcVx(frr(i),1) = vl1(i)
        elseif (m%fcVx(frr(i),2) .eq. vr1(i)) then
           m%fcVx(frr(i),2) = vl1(i)
        else
          call xerrab('set_connectivity_us: verts of radial face')
        endif
        
        
      enddo

      !Correct vertices of guard cells
      !left guard cell
       cvs = m%fcCv(frl(1),:)
       if ((cvs(1) .ne. cvl(1)) .and. (cvs(2) .eq. cvl(1))) then
         iCgl = cvs(1)
       elseif((cvs(2) .ne. cvl(1)) .and. (cvs(1) .eq. cvl(1))) then
         iCgl = cvs(2)
       else
         call xerrab('set_connectivity_us: no icgl found')
       endif
      !right cut
       cvs = m%fcCv(frr(1),:)
       if ((cvs(1) .ne. cvr(1)) .and. (cvs(2) .eq. cvr(1))) then
         iCgr = cvs(1)
       elseif ((cvs(2) .ne. cvr(1)) .and. (cvs(1) .eq. cvr(1))) then
         iCgr = cvs(2)
       else
         call xerrab('set_connectivity_us: no icgr found')
       endif  
    
       !vertex replacement for guard cells
       if (m%cvVx(m%cvVxP(iCgl,1)) .eq. vl1(1)) then
         m%cvVx(m%cvVxP(iCgl,1)) = vr1(1)
       elseif (m%cvVx(m%cvVxP(iCgl,1)+1) .eq. vl1(1)) then
         m%cvVx(m%cvVxP(iCgl,1)+1) = vr1(1)
       else
         call xerrab('set_connectivity_us: vl1(1) not in guard')
       endif
       if (m%cvVx(m%cvVxP(iCgr,1)) .eq. vr1(1)) then
         m%cvVx(m%cvVxP(iCgr,1)) = vl1(1)
       elseif (m%cvVx(m%cvVxP(iCgr,1)+1) .eq. vr1(1)) then
         m%cvVx(m%cvVxP(iCgr,1)+1) = vl1(1)
       else
         call xerrab('set_connectivity_us: vr1(1) not in guard')
       endif       


       !Correction at X-point
       do i = 1,m%cvVxP(m%nCv,1)+m%cvVxP(m%nCv,2)-1
           if (m%cvVx(i) .eq. vr2(ncvl)) then
              m%cvVx(i) = vl2(ncvl) 
           endif
       enddo
       do i = 1,m%nFc
          if (m%fcVx(i,1).eq.vr2(ncvl)) then
              m%fcVx(i,1) = vl2(ncvl)
          elseif (m%fcVx(i,2).eq.vr2(ncvl)) then
              m%fcVx(i,2) = vl2(ncvl)
          endif
       enddo
       !m%cvVx(pack(1:m%nCmVx, m%cvVx == vr2(ncvl))) = vl2(ncvl)
       !m%fcVx(pack(1:m%nFc, m%fcVx(:,1) == vr2(ncvl)),1) = vl2(ncvl)
       !m%fcVx(pack(1:m%nFc, m%fcVx(:,2) == vr2(ncvl)),2) = vl2(ncvl)

       call calc_add_connectivity(g,m)

      end subroutine set_connectivity_us
