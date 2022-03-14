*-----------------------------------------------------------------------
*.specification

      program b2grid
      use b2mod_grid_namelist
      implicit none
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2GRID (main program) is a pre-processor for b2.5. It generates a 
*     B2.5 grid for cylindrical symmetric geometries (Pilot-PSI, 
*     Magnum-PSI, Pisces and others). Next to that, it generates an input 
*     file for triangulation of the grid by Triangle.
*
*
*  2. specification
*
*     Main program.
*
*
*  3. description
*
*     this program generates a cylindrical symmetric B2.5 grid for 
*     modelling of linear plasma devices like Pilot-PSI, Magnum-PSI 
*     and Pisces. The input file is b2.grid.parameters. This gives 
*     the number of concentric coils, and for each coil:
*     r0: the inner radius of the coil
*     r1: the outer radius of the coil
*     z0: the axial start position
*     z1: the axial end position
*     current: the current through the coil
*     windings: the number of windings
*     nx: number of cells in axial direction
*     ny: number of cells in radial direction
*     nlayers: number of layers for each coil
*     stepsize: the initial spacial step size used in calculating the field lines							   																		 
*     Also the names of the 2-column files consisting of the profiles for 
*     radial and axial cell size distributions are given. Use fn_z and fn_r for that.
*
*     Optionally one could generate a triangular grid. This triangular grid 
*     generation uses Triangle (http://www.cs.cmu.edu/~quake/triangle.html). 
*     This tool has been added to the repository, and is called from within
*     b2grid. One can define a number of external region. For both the 
*     internal and external regions, Triangle will generate grids, which
*     will be automatically merged to one triangular grid.

*     Input options:
*     triangularGrid: logical, turning on or off the generation of the triangular grid
*     triangulationOption: 0 (default) splits B2 cells in 2 triangles, 1 splits the B2 cells in 4 triangles (by adding the centroid position as extra node.
*     nExternalRegions: defines the number of external regions.
*     nNodesInRegion(0:nExternalRegions-1): number of nodes used to define the external regions (minimally 3 per region)
*     area(0:nExternalRegions-1): upper limit to the surface area (in m^2) of triangular cells, to be defined per external region.
*     angle(0:nExternalRegions-1): lower limit to the angles (in degrees) of a triangular cell, to be defined per external region, typically around 30-35.
*     regionR(0:nNodesInRegion-1,0:nExternalRegions-1): radial positions of the nodes defining the external regions. 
*     regionZ(0:nNodesInRegion-1,0:nExternalRegions-1): see regionR, but now the z component.
*     startCellR(0:nExternalRegions-1): Radial (Y) index of the node on the B2.5 boundary
*     startCellZ(0:nExternalRegions-1): Axial (X) index of the node on the B2.5 boundary
*     endCellR(0:nExternalRegions-1): Radial (Y) index of the node on the B2.5 boundary
*     endCellZ(0:nExternalRegions-1): Axial (X) index of the node on the B2.5 boundary
* startCell and endCell should be nodes on the same B2.5 boundary. These nodes link the external grid to the B2.5 grid. The nodes, including the start and end cell, should be defined in counter-clockwise direction. The nodes at the boundary have index (0, :), (nx, :), (:, 0) or (0, ny) for the west, east, south and north boundary, respectively.
*
*-----------------------------------------------------------------------
*.declarations
c  Common dimensions
c
c  version : 01.12.98 21:42
c

c
c parameters that are common to Eirene and B2
c

c NOTE: DEF_NXD should not include the additional cells to handle the cuts
c*** Max. number of groups of Eirene surfaces for which the data can
c*** be transferred from B2 (DG specification "Surface special")

c new! [2002.04.22]
c new! [2002.06.14]

c
c parameters that are unique to B2
c


c
c parameters that are unique to Eirene
c


      integer ix, iy, nr, nz
      real (kind=R8) :: foo
      real (kind=R8) :: r_min
      real (kind=R8) :: r_max
      real (kind=R8), allocatable :: mr(:) !!@param starting positions in radial direction
      real (kind=R8), allocatable :: mz(:) !!@param starting positions in axial direction
      real (kind=R8), allocatable :: mz_unit(:) !!@param normalized distribution of axial positions on a magnetic field line
      real (kind=R8), allocatable :: cr(:) !!@param center between node/mesh positions in radial direction
      real (kind=R8), allocatable :: cz(:) !!@param center between node/mesh positions in axial direction
      real (kind=R8), allocatable, save :: mesh(:,:,:)
      real (kind=R8), allocatable :: r_profile(:), r_profile_pos_r(:), 
     *     r_profile_pos_z(:), z_profile(:), z_profile_pos(:)
      character*260 filename

*.computation
      filename = 'b2.grid.parameters'
      call alloc_b2mod_grid_namelist(200,100,NCOIL_MAX,
     *     NEXTERNAL_REGIONS_MAX, NNODES_IN_EXTERNAL_REGION_MAX)
      call read_b2mod_grid_namelist(filename,.true.)
      call write_b2mod_grid_namelist


      allocate(mesh(-1:nx+1,-1:ny+1,0:1))
      allocate(mr(0:ny))
      allocate(mz(0:ny))
      allocate(mz_unit(0:nx))
      allocate(cr(0:ny-1))
      allocate(cz(0:nx-1))

      if (b25grid) then
*     determine the maximum value r_max an r-coordinate of mesh point may get 
*     (equals the inner radius of the smallest coil)
        r_max = minval(r1(0:ncoil-1))

* read r and z distribution files.
        call read_first_line(fn_z, nz)
        call read_first_line(fn_r, nr)
*allocation of variables for profile
        allocate(r_profile_pos_r(0:nr-1))
        allocate(r_profile_pos_z(0:nr-1))
        allocate(r_profile(0:nr-1))
        allocate(z_profile_pos(0:nz-1))
        allocate(z_profile(0:nz-1))
        write (*, *) 'b2grid: nr, nz=', nr, nz
        call get_profile(fn_z, nz, z_profile_pos, z_profile)
        call get_profile(fn_r, nr, r_profile_pos_r, r_profile)

* determine r and z starting positions
        if (r_start.lt.r_end) then
          call cell_dist(mr, cr, ny, r_profile, r_profile_pos_r, nr,
     *         r_start, r_end)
        else
          call cell_dist(mr, cr, ny, r_profile, r_profile_pos_r, nr, 
     *         r_end, r_start)
        end if
        ! add z-positions (they are on the line rz_start - rz_end)
        mz(:) = z_start + (z_start - z_end) *
     *       (mr(:) - r_start) / (r_end - r_start)

        ! status: obtained the starting positions

        ! determine z_profile, normalized from position 0.0 to 1.0
        call cell_dist(mz_unit, cz, nx, z_profile, z_profile_pos, nz, 
     *       0.0_R8, 1.0_R8)

* start @ (r_profile_pos_r, ,r_profile_pos_z) calculating the fieldline in east and west direction, until z-coordinate passed z_east, z_west. Done in calc_mesh
        mesh(0:nx,0:ny, 0:1) = calc_mesh(mr, mz)

* grid is made, ready to write it into the different output formats.      
        call WRITE_OUT('OUT')

      else
        ! Read in OUT.ASCII
        call READ_IN('OUT')
      end if

      ! print some useful data
      ! the four corners
      WRITE (*, *) 'grid corners'
      WRITE (*, *) 'SW: ', mesh(0, 0, 1), mesh(0, 0, 0)
      WRITE (*, *) 'SE: ', mesh(nx, 0, 1), mesh(nx, 0, 0)
      WRITE (*, *) 'NW: ', mesh(0, ny, 1), mesh(0, ny, 0)
      WRITE (*, *) 'NE: ', mesh(nx, ny, 1), mesh(nx, ny, 0)
      if (triangularGrid) then
        call writeTriangleInput(mesh)
      end if

      contains
  
      function calc_mesh(mr_in, mz_in) result (meshOut)
      use b2mod_grid_namelist
      implicit none
*input variables
      real (kind=R8), intent(in) :: mr_in(0:ny), mz_in(0:ny)

*output variables
      real (kind=R8) :: meshOut(0:nx,0:ny, 0:1)

*local variables
      integer :: dir, ix, iy, z_index
      real (kind=R8) :: step, eps, 
     *     rEast, rWest, zEast, zWest, r, z,
     *     rEastOld, rWestOld, zEastOld, zWestOld
      real (kind=R8) :: zPositions(0:nx)
      external xertst

      eps = 1.0e-7_R8
      call xertst (z_nw.lt.z_ne, 
     *     'z_nw should be smaller than z_ne at this stage')

      ! determine the index of the node just before the start position z_start.
      do iy = 0, ny
        write (*, *) ' '
        write (*, *) 'calculate new magnetic field line', iy

        ! calculate initial fieldline
        step = stepsize

        ! obtain positions at wall
        ! from starting point in positive direction, towards east wall
        call calc_fieldline(mr_in(iy), mz_in(iy), -step, 
     *       r_se, z_se, r_ne, z_ne, rEast, zEast)
        ! from starting point in negative direction, towards west wall
        call calc_fieldline(mr_in(iy), mz_in(iy), step,
     *       r_sw, z_sw, r_nw, z_nw, rWest, zWest)

        do ! iterate, halving the stepsize until the changes to the fieldlines are small enough.
           ! this is measured by comparing the endpoints of two fieldlines of successive iterations.
          step = step*0.5_R8
          rWestOld = rWest
          rEastOld = rEast
          zWestOld = zWest
          zEastOld = zEast

          ! from starting point in positive direction, towards east wall
          call calc_fieldline(mr_in(iy), mz_in(iy), step, 
     *         r_se, z_se, r_ne, z_ne, rEast, zEast)
          
          ! from starting point in negative direction, towards west wall
          call calc_fieldline(mr_in(iy), mz_in(iy), -step,
     *         r_sw, z_sw, r_nw, z_nw, rWest, zWest)

          if (abs(rEast - rEastOld) < eps .and.
     *         abs(rWest - rWestOld) < eps .and.
     *         abs(zEast - zEastOld) < eps .and.
     *         abs(zWest - zWestOld) < eps) then
            ! inaccuracy is sufficiently small
            exit
          else
            write (*, *) 'stepsize is reduced to', step*0.5_R8
          endif
        end do

        ! we have correct points at the wall
        ! set z-positions
        meshOut(0:nx, iy, 0) = (mz_unit(:) * (zEast - zWest)) + zWest
        ! set first and last r position
!        meshOut(0, iy, 1) = rWest
!        meshOut(nx, iy, 1) = rEast

        ! find r positions belonging to the other z positions
        ! positive direction
        r = mr_in(iy)
        z = mz_in(iy)
        do ix = 0, nx
          if (z > meshOut(ix, iy, 0)) then
            cycle
          else
            ! determine r positions (using line perpendicular to z-axis)
            call calc_fieldline(r, z, step, 
     *            0.0_R8, meshOut(ix, iy, 0),
     *            1.0_R8, meshOut(ix, iy, 0), rEast, zEast)
            r = rEast
            z = zEast
            meshOut(ix, iy, 1) = rEast
          end if
        end do
        
        ! negative direction
        r = mr_in(iy)
        z = mz_in(iy)
        do ix = nx, 0, -1
          if (z < meshOut(ix, iy, 0)) then
            cycle
          else
            ! determine r positions (using line perpendicular to z-axis)
            call calc_fieldline(r, z, -step, 
     *            0.0_R8, meshOut(ix, iy, 0),
     *            1.0_R8, meshOut(ix, iy, 0), rWest, zWest)
            r = rWest
            z = zWest
            meshOut(ix, iy, 1) = rWest
          end if
        end do

      end do 
 
      return
      end function calc_mesh
 
      subroutine calc_fieldline(r_in, z_in, step, lr0, lz0, lr1, lz1, 
     *     rWall,zWall) 
      use b2mod_grid_namelist
      implicit none

*input variables
      integer index
      real (kind=R8), INTENT(IN) :: r_in, z_in !initial position, during calculation updated, following the fieldline
      real (kind=R8), intent(in) :: lr0, lz0 , lr1, lz1 ! line to which has to be calculated
      real (kind=R8), INTENT(IN) :: step !stepsize (positive or negative, depending on direction

*output variables
      real (kind=R8) :: rWall, zWall ! r and z position of field line at wall

*local variables
      real (kind=R8) :: b(0:1) !B-field (Br, Bz)
      real (kind=R8) :: weight, deriv, r, z
      real (kind=R8) :: rOld, zOld, crossProduct, crossProductOld
*calculation

      r = r_in
      z = z_in

      crossProduct = 
     .     (r - lr0) * (lz1 - lz0) - (z - lz0) * (lr1 - lr0)
      do ! while line did not hit wall
        ! store current position and crossProduct
        rOld = r
        zOld = z
        crossProductOld = crossProduct

        ! calc magnetic field in r,z
        b(0:1) = calc_b(r,z)
        ! update position
        deriv = b(0) / b(1)
        r = r + deriv * step
        z = z + step        

        ! calculate cross product with wall
        crossProduct = 
     .       (r - lr0) * (lz1 - lz0) - (z - lz0) * (lr1 - lr0)
        if (crossProduct * step < 0.0_R8) then
          ! interpolate to wall
          weight = crossProductOld / (crossProductOld - crossProduct)
          rWall = rOld + deriv * weight * step
          zWall = zOld + weight * step
          
          exit
        end if
      end do

      end subroutine calc_fieldline
      

      subroutine cell_dist(nodes,cc, n, profile, profile_pos, p,
     *     vmin, vmax)
* This subroutine calculates n+1 nodes and n cell centered positions (cc) 
* given a profile on p positions (profile_pos). 
* vmin and vmax are the edge positions of the profile, so finally their values should correspond with node(0) and node(n).   
*input variables
      implicit none
      integer n, p
      real (kind=R8) :: vmin, vmax
      real (kind=R8) :: profile(0:p-1), profile_pos(0:p-1)

*output variables
      real (kind=R8) :: nodes(0:n), cc(0:n-1)

*local variables
      real (kind=R8) :: delta, error, eps, tmp
      real (kind=R8) :: deriv(0:p-1), prof(0:n-1)
      integer is, ifail, cntr

*calculation
      eps=1.0e-10_R8
      error=10.0_R8*eps !at least larger than eps 
      call E01BEF(p,profile_pos(0:p-1),profile(0:p-1),
     *     deriv(0:p-1),ifail)
      delta = (vmax-vmin)/n
      nodes(0) = vmin
      nodes(n) = vmax
      do is = 1, n
        nodes(is) = vmin + delta * is
        prof(is-1) = nodes(is) - 0.5_R8 * delta
      end do
      cntr=0
      do
        if (cntr.gt.100) then
        write (*, *) 
     *          'calc_dist: 100 iterations done, break out of loop'
          exit
        end if
        do is = 0, n-1
          cc(is) = prof(is) 
        end do 
       
        call E01BFF(p,profile_pos(0:p-1),profile(0:p-1),deriv(0:p-1),
     *       n,cc(0:n-1),prof(0:n-1),ifail)
        nodes(0) = 0.0_R8
        do is = 0, n-1
          nodes(is+1) = nodes(is) + prof(is) 
        end do 
        nodes(0) = vmin
        tmp= (vmax - vmin) / nodes(n)
        error=0.0_R8
        do is = 1, n
          nodes(is) = nodes(is) * tmp + vmin
          prof(is-1) = prof(is-1)*tmp
          prof(is-1) = (nodes(is) + nodes(is-1)) * 0.5_R8
          error = error + (cc(is-1) - prof(is-1))**2.0_R8
        end do 
        error = error**0.5_R8 / n
        if(error.lt.eps) exit
        cntr = cntr + 1
      end do
      end subroutine

      function calc_b(r,z)
* this function returns the Br and Bz @ (r,z)
*input variables
      use b2mod_grid_namelist
      use b2mod_constants
      implicit none
      integer ncoils
      real (kind=R8) :: r,z

*output variables
      real (kind=R8) :: calc_b(0:1)

*local variables
      real (kind=R8) :: factor, dz1, dz2, dr, rlayer, bfield1(0:1), 
     *     bfield2(0:1), b1(0:1), b2(0:1)
      integer is, it

      calc_b = 0.0_R8
      DO is = 0, ncoil-1
        bfield1 = 0.0_R8
        bfield2 = 0.0_R8
        factor = current(is)*windings(is)/(Z2(is)-Z1(is))/nlayers(is)  !changed nlayers to array, so just changed to nlayers(is). 
        dr = (R2(is)-R1(is))/nlayers(is)
        DO it = 1, nlayers(is)
          rlayer=r1(is)+(it-0.5_R8)*dr
          b1(0:1) = calc_bb(z-z1(is),rlayer-r,rlayer+r)
          b2(0:1) = calc_bb(z-z2(is),rlayer-r,rlayer+r)
          bfield1(0) = bfield1(0) + (b2(0) - b1(0)) * factor * rlayer
          bfield1(1) = bfield1(1) + (b1(1) - b2(1)) * factor
        END DO
        factor = current(is)*windings(is)/(Z2(is)-Z1(is))/
     *   (2*nlayers(is))
        dr = (R2(is)-R1(is))/(2.0_R8*nlayers(is))
        DO it = 1, 2*nlayers(is)
          rlayer=r1(is)+(it-0.5_R8)*dr
          b1(0:1) = calc_bb(z-z1(is),rlayer-r,rlayer+r)
          b2(0:1) = calc_bb(z-z2(is),rlayer-r,rlayer+r)
          bfield2(0) = bfield2(0) + (b2(0) - b1(0)) * factor * rlayer
          bfield2(1) = bfield2(1) + (b1(1) - b2(1)) * factor
        END DO
*??? calc = calc + 2/3 * bfield2 + 1/3 * bfield1??? why assymetric? because of the 2 do loops (nlayers and 2nlayers)? but shouldn't it be the other way around then?
        calc_b(0) = calc_b(0) + bfield2(0) + 
     *       (bfield2(0) - bfield1(0)) / 3.0_R8  
        calc_b(1) = calc_b(1) + bfield2(1) + 
     *       (bfield2(1) - bfield1(1)) / 3.0_R8  
      END DO
      calc_b(0) = calc_b(0) * PI * 1.0e-7_R8
      calc_b(1) = calc_b(1) * PI * 1.0e-7_R8
      end function calc_b
     
      function calc_bb(zz,rm,rp)
      use b2mod_grid_namelist
      implicit none
* input variables
* input (rm,zz) and (rp,zz) are two coordinates...
      real (kind=R8) :: zz, rm, rp

*output variables
      real (kind=R8) :: calc_bb(0:1)

*local variables
      integer cntr, cntr_max
      real (kind=R8) :: tmp1, tmp2, eps, a, b, b0, c, d, d0, m, n, p, 
     *     xmn, xmp
*.computation
      cntr_max=50
      eps = 1.0e-12_R8
      tmp1 = sqrt(rp*rp + zz*zz)
      n = sqrt(rm*rm + zz*zz)/tmp1
      if (n.lt.eps) n = eps
      m = 1.0_R8
      a = -1.0_R8
      b = 1.0_R8
      p = rm / rp
      c = p + 1.0_R8
      if(abs(p).le.eps) p = 1.0_R8
      d = c
      if (p.lt.0.0_R8) then
        d = -d
        p = -p
      end if
      cntr=0
      do                        !recursively converge to correct value of B
        if (abs(m - n).lt.eps * m) EXIT !converged
        IF (cntr >= cntr_max) then
          WRITE (*, *) 'zz, cntr, m, n', zz, cntr, m, n
        END IF
        call xertst(cntr.lt.cntr_max, 
     *       'B calculation did not converge!!!')
        xmn = m * n
        xmp = xmn / p
        d0 = d
        d = 2.0_R8 * (d + xmp * c)
        c = c + d0 / p
        b0 = b
        b = 2.0_R8 * (b + n * a)
        a = a + b0 / m
        p = p + xmp
        m = m + n
        n = 2.0_R8 * xmn**0.5_R8
        cntr = cntr + 1
      end do
      calc_bb(0) = (a + b / m) / (m * tmp1)
      calc_bb(1) = (c + d / m) / ((m + p) * tmp1) * zz
      return
      end function calc_bb

      SUBROUTINE WRITE_OUT(FN)
      use b2mod_grid_namelist
      IMPLICIT NONE
      CHARACTER*(*) FN
      CHARACTER*1 S1,S2
      INTEGER ix, iy, is, col, cntr, row
      real (kind=R8) :: bToroidal, one, cr, cz
      real (kind=R8) :: b(0:1) !B-field (Br, Bz)

      ! mesh(0:nx,0:ny,0:1) contains nodes of computational grid
      ! to add the guard cells

* construct guard cells around the grid
* The guard cells are added around the currrent mesh. 
* The additional nodes are calculated by a linear extrapolation 
* from the last two nodes at the edge.
      write (*, *) 'gc_nsew:', gc_n, gc_s, gc_e, gc_w
      forall (ix=0:nx, is=0:1)
        mesh(ix,-1,is) = (1.0_R8 + gc_s) * mesh(ix,0,is) - 
     *     gc_s * mesh(ix,1,is)
        mesh(ix,ny+1,is) = (1.0_R8 + gc_n) * mesh(ix,ny,is) - 
     *     gc_n * mesh(ix,ny-1,is)
      end forall
      forall (iy=-1:ny+1, is=0:1) ! -1:ny+1 to add the corner cells too.
        mesh(-1,iy,is) = (1.0_R8 + gc_w) * mesh(0,iy,is) - 
     *     gc_w * mesh(1,iy,is)
        mesh(nx+1,iy,is) = (1.0_R8 + gc_e) * mesh(nx,iy,is) - 
     *     gc_e * mesh(nx-1,iy,is)
      end forall

* write the mesh to the output in 3 formats
* .ASCII has been validated, other formats might be incorrect
      S1='='
      S2='-'
      ! assuming zero toroidal field everywhere
      bToroidal = bTor
      one = 1.0_R8
      OPEN(UNIT=16,FILE=FN//'.ELEMENT',STATUS='UNKNOWN')
      OPEN(UNIT=17,FILE=FN//'.ASCII',STATUS='UNKNOWN')
      OPEN(UNIT=30,FILE='fort.30',STATUS='UNKNOWN')

      write (30,*) 'Geometry from B2'
      write (30,*) ' '
      write (30,*) nx, ny, 0
      write (30,*) ' 0 0 0 0 '
      write (30,*) ' '
      
      WRITE(17,'(I3,I3)') nx, ny      !RW + 2 guard cells
      WRITE(16,*)
      WRITE(16,*) '   Element output: '
      WRITE(16,*)
      WRITE(16,'(A3,88A1)')	'   ',(S1,is=1,88)
      col = -1
      cntr = -1 ! for index in OUT.ELEMENTS

!      DO iy = ny+1, 0, -1
      DO iy = 0, ny+1
        col = col + 1
        row=-1
!        DO ix = 0, nx+1
        DO ix = nx+1, 0, -1
          row = row + 1
          cntr = cntr + 1
          write (30,'(4e15.7)') mesh(IX-1,IY,1), mesh(IX,IY,1),
     *          mesh(IX-1,IY-1,1), mesh(IX,IY-1,1)  
          write (30,'(4e15.7)') mesh(IX-1,IY,0), mesh(IX,IY,0),
     *         mesh(IX-1,IY-1,0), mesh(IX,IY-1,0)  
          cr = (mesh(IX-1,IY-1,1) + mesh(IX,IY-1,1) + mesh(IX-1,IY,1) + 
     *         mesh(IX,IY,1)) * 0.25_R8
          cz = (mesh(IX-1,IY-1,0) + mesh(IX,IY-1,0) + mesh(IX-1,IY,0) + 
     *         mesh(IX,IY,0)) * 0.25_R8
          b(0:1) = calc_b(cr, cz)
          WRITE(17,'(I3,I3,12E28.18)') row, col, CR, CZ,
     *         mesh(IX,IY-1,1), mesh(IX,IY-1,0),
     *         mesh(IX-1,IY-1,1), mesh(IX-1,IY-1,0),
     *         mesh(IX,IY,1), mesh(IX,IY,0),
     *         mesh(IX-1,IY,1), mesh(IX-1,IY,0), b(1), bToroidal

          WRITE(16,'(A10,I5,A4,I3,A1,I3,A4,E17.10,A1,E17.10,A8,
     *         E17.10,A1,E17.10,A1)')
     *         '   Element',CNTR,' = (',ix,',',col,'): (',
     *         mesh(IX-1,IY-1,1),',',mesh(IX-1,IY-1,0),')      (',
     *         mesh(IX,IY-1,1),',', mesh(IX,IY-1,0),')'
          WRITE(16,'(A18,E17.10,A14,E17.10,A1,E17.10,A1)') 
     *         '   Field ratio  = ',one,'             (',CR,',',CZ,')'
          WRITE(16,'(A30,E17.10,A1,E17.10,A8,E17.10,A1,E17.10,A1)') 
     *         '  (',mesh(IX-1,IY,1),',',mesh(IX-1,IY,0),
     *         ')      (',mesh(IX,IY,1),',',mesh(IX,IY,0),')'
          WRITE(16,'(A3,88A1)')	'   ',(S2,is=1,88)
        end do
      end do

      CLOSE(17)
      CLOSE(16)
      CLOSE(30)

      RETURN

      END SUBROUTINE WRITE_OUT


      SUBROUTINE READ_IN(FN)
      use b2mod_grid_namelist
      IMPLICIT NONE
      CHARACTER*(*) FN
      CHARACTER*1 S1,S2
      INTEGER ix, iy, is, col, cntr, row
      real (kind=R8) :: bToroidal, one, cr, cz
      real (kind=R8) :: b(0:1) !B-field (Br, Bz)

* write the mesh to the output in 3 formats
* .ASCII has been validated, other formats might be incorrect
      S1='='
      S2='-'
      ! assuming zero toroidal field everywhere
      OPEN(UNIT=17,FILE=FN//'.ASCII',STATUS='UNKNOWN')

      
      READ(17,'(I3,I3)') nx, ny      !RW + 2 guard cells
      col = -1
      cntr = -1 ! for index in OUT.ELEMENTS

!      DO iy = ny+1, 0, -1
      DO iy = 0, ny+1
        col = col + 1
        row=-1
!        DO ix = 0, nx+1
        DO ix = nx+1, 0, -1
          ! read statement similar to write statement in WRITE_OUT to avoid errors
          ! disadvantage: reading in every node 4 times (if it is a disadvantage)
          READ(17,'(I3,I3,12E28.18)') row, col, CR, CZ,
     *         mesh(IX,IY-1,1), mesh(IX,IY-1,0),
     *         mesh(IX-1,IY-1,1), mesh(IX-1,IY-1,0),
     *         mesh(IX,IY,1), mesh(IX,IY,0),
     *         mesh(IX-1,IY,1), mesh(IX-1,IY,0), b(1), bToroidal

        end do
      end do

      CLOSE(17)

      RETURN

      END SUBROUTINE READ_IN


      !> \brief This routine writes an input file for the internal trianglular grid and when present, also input files for the external triangular grids.
      !! The input file for the internal grid consists of 
      !! - a list of all nodes
      !! - a list of all cell sides (to prevent triangles to cross cell sides)
      subroutine writeTriangleInput(mesh)
      use b2mod_grid_namelist

      implicit none
      character*120 fnNode0, fnNode1, fnPoly0, fnPoly1, regionChar, line
      integer ix, iy, is, cntr, step, additionalNodes, wallIndex
      integer, dimension(-1:nExternalRegions - 1) :: nrSides, nrNodes
      integer, dimension(0:3) :: intData
      integer firstIX, firstIY, lastIX, lastIY
      real (kind=R8) :: mesh(-1:nx+1,-1:ny+1,0:1)
      real (kind=R8), dimension(1:2) :: realData

      ! first create the file for the internal triangular grid
      fnPoly0 = 'internal.poly'
      open(unit=20, file=fnPoly0, status='UNKNOWN')

      ! write first line
      write (20, *) '# Set the nodes of the internal grid'

      select case (triangulationOption)
      case (0)
        ! number of nodes to be written is (nx+1) * (ny+1)
        write (20, *) (nx+1)*(ny+1), 2, 0, 0

      case (1)
        ! also add centroid positions
        ! number of nodes to be written is (nx+1) * (ny+1) + nx*ny
        write (20, *) (nx+1)*(ny+1) + nx*ny, 2, 0, 0

      case default
        call xerrab('Invalid triangulation option')
      end select

      ! write nodes
      cntr = 0
      do iy = 0, ny
        do ix = 0, nx
          write (20, *) cntr, mesh(ix, iy, 1), mesh(ix, iy, 0)
          cntr = cntr + 1
        end do
      end do
      
      select case (triangulationOption)
      case (0)
        ! nothing to be done

      case (1)
        ! also add centroid positions
        do iy = 0, ny - 1
          do ix = 0, nx - 1
            write (20, *) cntr,
     *           0.25_R8 * (mesh(ix, iy, 1) + mesh(ix+1, iy, 1) + 
     *           mesh(ix+1, iy+1, 1) + mesh(ix, iy+1, 1)),
     *           0.25_R8 * (mesh(ix, iy, 0) + mesh(ix+1, iy, 0) + 
     *           mesh(ix+1, iy+1, 0) + mesh(ix, iy+1, 0))
            cntr = cntr + 1
          end do
        end do

      case default
        call xerrab ('Invalid triangulation option')
      end select

      ! add blank line
      write (20, *) ''

      ! add the list of b2.5 cell sides
      ! write comments line first
      write (20, *) '# set the sides'
      write (20, *) '# format: counter, 2 nodes, wall index'
      
      ! write data
      ! nr of sides to be written is nx*(ny+1) + ny*(nx+1) (horizontally and vertically)
      nrSides(-1) = (nx)*(ny+1) + ny*(nx+1)
      write (20, *) nrSides(-1), 1
      cntr = 0
      ! horizontal sides
      do iy = 0, ny
        do ix = 0, nx - 1
          write (20, *) cntr, ix + iy*(nx+1), ix + iy*(nx+1) + 1,
     *         cntr + 1000
          cntr = cntr + 1
        end do
      end do
      ! vertical sides
      do ix = 0, nx
        do iy = 0, ny - 1
          write (20, *) cntr, ix + iy*(nx+1), ix + (iy+1)*(nx+1),
     *         cntr + 1000
          cntr = cntr + 1
        end do
      end do

      ! add blank line and final line containing a zero
      write (20, *) ''
      write (20, *) 0

      ! close file containing information for the internal triangular grid
      close(20) 

      ! set symbolic link to triangle in current directory
      call system('ln -sf ${SOLPSTOP}/src/Triangle/triangle triangle')
      ! run triangle for the internal grid
      call system('./triangle -zepn internal.poly')
     
      ! print some useful info to the screen
      write (*, *) 'Checking external region definitions'

      ! create input files for the external regions
      do is = 0, nExternalRegions - 1

        ! first and last node are obtained from input
        ! should coincide with boundary nodes of the B2.5 grid.
        firstIX = startCellZ(is)
        firstIY = startCellR(is)
        lastIX = endCellZ(is)
        lastIY = endCellR(is)

        cntr = nNodesInRegion(is) - 1 ! index of the last node

        WRITE (*, *) 'indices of first and last node on b2.5 boundary'
        WRITE (*, *) firstIX, firstIY
        WRITE (*, *) lastIX, lastIY
        WRITE (*, *) 'regionR', regionR(0:10, is)
        WRITE (*, *) 'regionZ', regionZ(0:10, is)


        ! first node should be on some boundary
        call xertst(firstIX == 0 .or. firstIX == nx .or.
     *       firstIY == 0 .or. firstIY == ny,
     *       'startCellR or startCellZ should be on a b2 boundary node')

        call xertst(firstIX == lastIX .or. firstIY == lastIY,
     *       'first and last node should be on same b2 boundary')

        ! points along boundary of b2 grid in between first and last node
        ! have to be added to triangle file.
        ! determine number of nodes between first and last node 
        ! (including first and last node)
        additionalNodes = lastIX + lastIY - firstIX - firstIY
        if (additionalNodes .lt. 0) then
          additionalNodes = -additionalNodes
          step = 1 ! index of last < first, so loop from last node to first node is in positive direction
        else
          step = -1 ! index of last > first, so loop from last node to first node is in negative direction
        end if
        additionalNodes = additionalNodes + 1 ! adding 1 to include first and last node.

        ! create the file for the external triangular grid 
        ! (1 for each external region)
        write (regionChar, *) is
        fnPoly0 = 'external' // trim(adjustl(regionChar)) // '.poly'
        open(unit=20, file=fnPoly0, status='UNKNOWN')

        ! write first line
        write (20, *) '# Set the nodes of this external grid'

        write (20, *) nNodesInRegion(is) + additionalNodes, 2, 0, 0
        cntr = 0
        ! write external nodes
        do ix = 0, nNodesInRegion(is) - 1
          write (20, *) cntr, regionR(ix, is), regionZ(ix, is)
          cntr = cntr + 1
        end do
        ! write nodes on B2 boundary
        if (firstIX == lastIX) then
          ! loop over nodes in between lastIY and firstIY
          do iy = lastIY, firstIY, step
            write (20, *) cntr, mesh(firstIX, iy, 1), 
     *           mesh(firstIX, iy, 0)
            cntr = cntr + 1
          end do
        else
          ! loop over nodes in between lastIY and firstIY
          do ix = lastIX, firstIX, step
            write (20, *) cntr, mesh(ix, firstIY, 1), 
     *           mesh(ix, firstIY, 0)
            cntr = cntr + 1
          end do
        end if

        ! add blank line
        write (20, *) ''

        ! add the list of b2.5 cell sides
        ! write comments line first
        write (20, *) '# set the sides'
        write (20, *) '# format: counter, 2 nodes, wall index'
 
        ! write sides of external grid
        ! number of sides equal number of nodes (for closed polygon)
        nrSides(is) = nNodesInRegion(is) + additionalNodes
        write (20, *) nrSides(is), 1
        cntr = 0
        do ix = 0, nNodesInRegion(is) + additionalNodes - 2
          write (20, *) cntr, ix, ix + 1, cntr + 1000
          cntr = cntr + 1
        end do
        ! close polygon (last node to first node)
        write (20, *) cntr, nrSides(is) - 1, 0,
     *       cntr + 1000

        ! add blank line and final line containing a zero
        write (20, *) ''
        write (20, *) 0

        close(20)
        call system('./triangle -zepn -a' // trim(area(is)) //
     *       ' -q' // trim(angle(is)) // ' ' // trim(fnPoly0))
      end do

      ! Merge the internal triangular grid and the external triangular grids.
      ! Hereto, we merge the *.poly and *.node files.
      ! The merged file will be called 'tria'
      
      fnNode0 = 'tria.node'
      fnPoly0 = 'tria.poly'

      open(unit=20, file=fnNode0, status='UNKNOWN')
      open(unit=30, file=fnPoly0, status='UNKNOWN')

      !!! merge node files

      ! determine total number of nodes of the different regions
      fnNode1 = 'internal.1.node'
      open(unit=120, file=fnNode1, status='UNKNOWN')
      read (120, *) nrNodes(-1), intData(1:3)
      close (120)
 
      do is = 0, nExternalRegions - 1
        ! construct filename of external grid
        write (regionChar, *) is
        fnNode1 = 'external' // trim(adjustl(regionChar)) // '.1.node'
        open(unit=120, file=fnNode1, status='UNKNOWN')

        read (120, *) nrNodes(is), intData(1:3)
        close (120)
      end do

      ! determine total number of sides of the different regions
      fnPoly1 = 'internal.1.poly'
      open(unit=130, file=fnPoly1, status='UNKNOWN')
      ! skip first line
      read (130, *) line
      ! extract nr of sides
      read (130, *) nrSides(-1), intData(1)
      close (130)
 
      do is = 0, nExternalRegions - 1
        ! construct filename of external grid
        write (regionChar, *) is
        fnPoly1 = 'external' // trim(adjustl(regionChar)) // '.1.poly'
        open(unit=130, file=fnPoly1, status='UNKNOWN')

        ! skip first line
        read (130, *) line
        ! extract nr of sides
        read (130, *) nrSides(is), intData(1)
        close (130)
      end do

      !!!! start constructing merged node file
      ! open internal grid
      fnNode1 = 'internal.1.node'
      open(unit=120, file=fnNode1, status='UNKNOWN')

      ! some comments in file
      write (20, *) '# set the nodes of the merged grid'
      ! write line containing info on the number of nodes.
      write (20, *) sum(nrNodes(-1:nExternalRegions - 1)), 2, 0, 1

      ! skip first line
      read (120, *) line

      wallIndex = 1000
      write (20, *) '# from ', trim(fnNode1)
      do ix = 0, nrNodes(-1) - 1
        read (120, *) cntr, realData(1:2), intData(3)
        write (20, *) ix, realData(1:2), wallIndex
        wallIndex = wallIndex + 1
      end do
      close(120)

      ! add the nodes of the external regions
      do is = 0, nExternalRegions - 1
        ! construct filename of external grid
        write (regionChar, *) is
        fnNode1 = 'external' // trim(adjustl(regionChar)) // '.1.node'
        write (20, *) '# from ', trim(fnNode1)        
        ! open file
        open(unit=120, file=fnNode1, status='UNKNOWN')

        ! skip first line
        read (120, *) line
        do ix = sum(nrNodes(-1:is-1)), sum(nrNodes(-1:is)) - 1
          read (120, *) cntr, realData(1:2), intData(3)
          write (20, *) ix, realData(1:2), wallIndex
          wallIndex = wallIndex + 1
        end do
        close(120)
      end do

      !!!! merge poly files
      fnPoly1 = 'internal.1.poly'
      open(unit=130, file=fnPoly1, status='UNKNOWN')
      write (30, *) '# set the sides of the merged grid'
      ! copy first line of poly file
      read (130, *) intData(0:3)
      write (30, *) intData(0:3)
      ! construct second line
      read (130, *) line ! This is the line with the number of sides of the internal grid. Skip it
      write (30, *) sum(nrSides(:)), 1
      write (30, *) '# format: counter, 2 nodes, wall index'
      write (30, *) '# from ', trim(fnPoly1)

      cntr = 0
      wallIndex = 1000
      ! copy all side info
      do ix = 1, nrSides(-1)
        read (130, *) intData(0:3)
        write (30, *) intData(0:2), wallIndex
        wallIndex = wallIndex + 1
       cntr = cntr + 1
      end do
      ! finshed copying side info of internal grid
      close(130)
      
      ! copy sides of all external grids.
      do is = 0, nExternalRegions - 1
        ! construct filename of external grid
        write (regionChar, *) is
        fnPoly1 = 'external' // trim(adjustl(regionChar)) // '.1.poly'
        
        ! open file
        open(unit=130, file=fnPoly1, status='UNKNOWN')

        ! skip first 2 lines
        read (130, *) line
        read (130, *) line

        ! add comments
        write (30, *) '# from ', trim(fnPoly1)

        ! copy data
        do ix = 1, nrSides(is)
          read (130, *) intData(0:3)
          write (30, *) cntr, intData(1:2) + sum(nrNodes(-1:is-1)),
     *         intData(3) + (is+1)*10000
          cntr = cntr + 1
        end do
        
        close(130)
      end do
      ! finished copying all external sides.
      ! write last line
      write (30, *) 0
      
      close(20)
      close(30)

      call system('./triangle -zepnj tria.poly')

      return

      end subroutine writeTriangleInput

      end program
