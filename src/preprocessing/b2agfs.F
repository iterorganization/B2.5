*-----------------------------------------------------------------------
*.specification

*      subroutine b2agfs (nrid,nCv,nFc,nVx,nx,ny,nx1,ny1,istyle,
*     & crx1,cry1,fpsi1,ffbz1,bb1,wbbc1,
*     1 cflag,lblgm,leftix1,leftiy1,rightix1,rightiy1,
*     2 topix1,topiy1,bottomix1,bottomiy1,
*     3 leftcut1,rightcut1,bottomcut1,topcut1,
*     4 periodic_bc,nncut,nncutmax)

      subroutine b2agfs (nrid,nCv,nFc,nVx,istyle,lblgm,g,m)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
!      use b2mod_connectivity
      use b2mod_indirect
      use b2ag_parameters
      use b2ag_ghostcells
      use carre_constants
      use b2mod_geo
      implicit none
#include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nrid, nCv, nFc, nVx, istyle
      type(geometry_ag), intent(inout) :: g
      type(mapping_ag), intent(inout) :: m
*   .. output arguments
!      real (kind=R8) ::
!     * crx1(-1:nx1,-1:ny1,0:3), cry1(-1:nx1,-1:ny1,0:3),
!     1 bb1(-1:nx1,-1:ny1,0:3),fpsi1(-1:nx1,-1:ny1,0:3,0:2),
!     2 ffbz1(-1:nx1,-1:ny1,0:3),wbbc1(-1:nx1+1,-1:ny1+1,0:3)
!     integer
!     *  leftix1(-1:nx1,-1:ny1),leftiy1(-1:nx1,-1:ny1),
!     1  rightix1(-1:nx1,-1:ny1),rightiy1(-1:nx1,-1:ny1),
!     2  topix1(-1:nx1,-1:ny1),topiy1(-1:nx1,-1:ny1),
!     3  bottomix1(-1:nx1,-1:ny1),bottomiy1(-1:nx1,-1:ny1),
!     4  leftcut1(nncutmax),rightcut1(nncutmax),
!     5  bottomcut1(nncutmax),topcut1(nncutmax),
!     6  periodic_bc,nncut,cflag(-1:nx1,-1:ny1,GRID_N_CELLFLAGS)
!      integer periodic_bc
      character*(*) lblgm
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGFS obtains a description of the flux surface structure.
*
*  5. parameters (see also routine b2cdcv)
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     (parameters in common:)
*
*     /b2cmfs/ - output.
*     /b2cmfs/ is reserved for numerical mhd equilibrium data.
*     (Presently empty--to be used in conjunction with grid generation
*     based on a numerical mhd equilibrium.)
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
*     periodic_bc - integer, read from file.
*     periodic_bc specifies if this is either an island or limiter geometry. 
*     If periodic_bc.eq.1 then island/limiter treatment is turned on.
*     We differentiate between the two case through nncut:
*      nncut.eq.0 = limiter case
*      nncut.ge.1 = island divertor case (there should be nncut islands then)
*     The case periodic_bc.eq.-1 is used to remove tranverse physics in
*     1-D cases.
*     If periodic_bc.eq.0 then we have the usual treatment.
*     No other values are allowed.
*
*     xoffset, yoffset - real*8, read from file.
*     xoffset and yoffset are offsets of the mesh in the x- and y-direction
*     respectively. The mesh will be translated by (xoffset,yoffset) if any
*     of the two is nonzero.
*
*     Bt_adjust - integer, read from file.
*     If the mesh is offset and Bt_adjust.eq.1, then the magnetic field is
*     recomputed, keeping the pitch constant but changing the toroidal field
*     magnitude according to Bt = Bt0/R.
*
*     Bt_rescale - real*8, read from file.
*     The magnetic field will be multiplied by Bt_rescale.
*     Default value is 1. It will hold that 0.lt.Bt_rescale.
*
*     min_pitch - real*8, read from file.
*     Minimum allowed value for the pitch angle (in degrees) at the plates.
*
*  6. periodic boundary condition treatment
*
*     added by XPB on 20.09.1999
*     This assumes only one cut and is coded with the diverted island
*     and limiter geometries in mind when periodic_bc.eq.1.
*     The program looks for non-matching sides and assigns the
*     neighbourhoods array from one cut area to the other.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ix, iy, i, ic, ixl, ixr, iyl, iyr, inseltop, inselbot
      integer Bt_adjust, Bt_reversal
      character*256 local_sonnet
      data local_sonnet /'upgrade.geometry'/
      real (kind=R8) :: wbbc_us(0,nCv-1,0:3)
!     *     crxs(-1:nx1,-1:ny1,0:4),
!     1     crys(-1:nx1,-1:ny1,0:4),
!     2     bp(-1:nx1,-1:ny1),
!     3     bt(-1:nx1,-1:ny1),
!     4     psi(-1:nx1,-1:ny1,0:4,0:2), ! Dims: 1,2: cell position, 3: center=0/corners=1..4, 4: value, ddx, ddy
!     5     ffbz(-1:nx1,-1:ny1,0:3), ! ffbz at cell corners
!     6     wbbc(-1:nx1+1,-1:ny1+1,0:3)
      real (kind=R8) ::
     *  xoffset, yoffset, xrescale, yrescale, Bt_rescale, pit_rescale,
     .  min_pitch, geom_match_dist, distval,
     .  r1,r2,r3,rf,z1,z2,z3,z4,numer,denom,aaa,bbb,
     .  nhatr,nhatz,bphatr,bphatz,bdotn,beta,cosa,cosb,cose,sinb
      integer ireg,errors
!      character*1 zone(-1:nx1)
      logical found, concave
      data xoffset /0.0_R8/, yoffset /0.0_R8/, Bt_rescale/1.0_R8/,
     .  Bt_adjust/0/, Bt_reversal/0/, xrescale/1.0_R8/, 
     .  yrescale/1.0_R8/, pit_rescale/1.0_R8/, min_pitch/1.0_R8/, 
     .  geom_match_dist/1.0e-6_R8/

*   ..procedures
      external subini, subend, ipgetc, ipgeti, ipgetr, xertst,
     * b2xvsg, b2gxdr_dummy
      intrinsic abs, sqrt, cos, sin, min, max, sign, acos
!      real (kind=R8) ::
!     *  dist
!      integer ix1,iy1,ip1,ix2,iy2,ip2
!      dist(ix1,iy1,ip1,ix2,iy2,ip2)=
!     1 sqrt((crx1(ix1,iy1,ip1)-crx1(ix2,iy2,ip2))**2+
!     2 (cry1(ix1,iy1,ip1)-cry1(ix2,iy2,ip2))**2)
!      logical match
!      match(ix1,iy1,ix2,iy2)=
!     1 (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt.
!     2 geom_match_dist

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2agfs')
*   ..test nCv, nFc, nVx, nrid
      call xertst (0.le.nCv.and.0.le.nFc.and.0.le.nVx,
     &   'faulty argument nCv, nFc, nVx')
      call xertst (0.lt.nrid, 'faulty argument nrid')

      if(istyle.ge.0) goto 9999

      call ipgetc ('b2agfs_geometry', local_sonnet)
      lblgm=local_sonnet
      periodic_bc=0
      call ipgeti ('b2agfs_periodic_bc',periodic_bc)
      call xertst (abs(periodic_bc).le.1, 'faulty argument periodic_bc')
      call ipgetr ('b2agfs_xoffset',xoffset)
      call ipgetr ('b2agfs_yoffset',yoffset)
      call ipgetr ('b2agfs_xrescale',xrescale)
      call xertst (0.0_R8.ne.xrescale,'faulty argument xrescale')
      call ipgetr ('b2agfs_yrescale',yrescale)
      call xertst (0.0_R8.ne.yrescale,'faulty argument yrescale')
      call ipgeti ('b2agfs_Bt_adjust',Bt_adjust)
      call xertst (Bt_adjust.eq.0 .or. Bt_adjust.eq.1,
     &  'faulty argument Bt_adjust')
      call ipgetr ('b2agfs_Bt_rescale',Bt_rescale)
      call xertst (0.0_R8.ne.Bt_rescale,'faulty argument Bt_rescale')
      call ipgeti ('b2agfs_Bt_reversal', Bt_reversal)
      call xertst (Bt_reversal.eq.0 .or. Bt_reversal.eq.1,
     &  'faulty argument Bt_reversal')
      call ipgetr ('b2agfs_pit_rescale', pit_rescale)
      call ipgetr ('b2agfs_min_pitch', min_pitch)
      call xertst (min_pitch.ge.0.0_R8 .and. min_pitch.lt.90.0_R8,
     &  'faulty argument min_pitch')
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')

      ! Read in grid input file (as typically produced by traduit)
      ! This routine also ensures the existence of ghost cells
      call b2gxdr_dummy(local_sonnet,nCv,nFc,nVx,istyle,
     &     wbbc_us,isymm,g,m)


!       ! all on the nx1, ny1 grid, just do all on regular grid


      ! Set up magnetic field and shift coordinates
      do i = 0,nCv-1
        if ( g%cvCflags(i,0) .eq. GRID_UNDEFINED .or. 
     &     g%cvCflags(i,0) .eq. GRID_EXTERNAL .or. 
     &     g%cvCflags(i,0) .eq. GRID_DEAD ) cycle  !isUnusedCell works only for unstructured
        g%cvBb(i,0) = g%cvBp(i)
        g%cvBb(i,1) = 0.0_R8
        g%cvBb(i,2) = g%cvBt(i)
        g%cvBb(i,3) = sqrt( g%cvBp(i)**2 + g%cvBt(i)**2 ) 
        g%vxX(i) = g%vxX(i)*xrescale+xoffset
        g%vxY(i) = g%vxY(i)*yrescale+xoffset
        !g%vxFpsi(i) = g%vxFpsi(i)
        !g%vxFfbz(i) = g%vxFfbz(i)
        !wbbc1_us(i) = wbbc_us(i)
!          if (ix.eq.nx1 .and. iy.lt.ny1) then
!           wbbc1_us(nx1+1,iy,0:3) = wbbc_us(nx1+1,iy,0:3)
!          else if (ix.lt.nx1 .and. iy.eq.ny1) then
!            wbbc1_us(ix,ny1+1,0:3) = wbbc1_us(ix,ny1+1,0:3)
!          else if (ix.eq.nx1 .and. iy.eq.ny1) then
!            wbbc1_us(nx1+1,ny1+1,0:3) = wbbc_us(nx1+1,ny1+1,0:3)
!          end if
        if (Bt_adjust.eq.1 .and. xoffset.ne.0.0_R8) then
          g%cvBb(i,2) = g%cvBt(i)*g%cvX(i)/(g%cvX(i)+xoffset)
          g%cvBb(i,0) = g%cvBp(i)/g%cvBt(i)*g%cvBb(i,2)
          g%cvBb(i,3) = sqrt(g%cvBb(i,0)**2 + g%cvBb(i,2)**2)
        endif
      enddo
      g%cvBb(:,0)=g%cvBb(:,0)*Bt_rescale*pit_rescale
      g%cvBb(:,1)=g%cvBb(:,1)*Bt_rescale*pit_rescale
      if (Bt_reversal.eq.1) then
         g%cvBb(:,2) = -g%cvBb(:,2)*Bt_rescale
      else
         g%cvBb(:,2) = g%cvBb(:,2)*Bt_rescale
      endif
      g%cvBb(:,3) = sqrt(g%cvBb(:,0)**2+g%cvBb(:,1)**2+g%cvBb(:,2)**2)
          
 9999 continue
      ! main connectivity already established in Make_traduit_b2us
      ! MAYBE HERE BUILD OTHER CONNECTIVITIES NEEDED FOR B2FGMTRY_us? OR TO EARLY
!      call init_connectivity(nx1,ny1,crx1,cry1,cflag,
!     & leftix1,leftiy1,rightix1,rightiy1, 
!     & topix1,topiy1,bottomix1,bottomiy1, 
!     & leftcut1,rightcut1,bottomcut1,topcut1, 
!     & periodic_bc,nncut,nncutmax, inseltop, inselbot, geom_match_dist, 
!     & istyle )

!      call test_connectivity(nx1,ny1,crx1,cry1,cflag,
!     & leftix1,leftiy1,rightix1,rightiy1, 
!     & topix1,topiy1,bottomix1,bottomiy1 )      

      !END OF CONNECTIVITIES CALCULATIONS


cxpb For classical grids only
cxpb Check the pitch angle at the plates and make sure it is at least min_pitch
cxpb Do not correct corner cells
cxpb ONLY AT TARGET PLATES
      if (min_pitch.gt.0.0_R8 .and.m%isClassicalGrid.eq.1) then
!        do ix = -1, nx1
!          do iy = 0, ny1-1
!!             Skip unused cells
!             if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
!
cxpb Western edge
!            if (.not.isInDomain(nx1,ny1,leftix1(ix,iy),leftiy1(ix,iy)))
!     >       then
!              r1 = crx1(rightix1(ix,iy),rightiy1(ix,iy),2)
!              z1 = cry1(rightix1(ix,iy),rightiy1(ix,iy),2)
!              r2 = crx1(rightix1(ix,iy),rightiy1(ix,iy),0)
!              z2 = cry1(rightix1(ix,iy),rightiy1(ix,iy),0)
!              r3 = crx1(rightix1(ix,iy),rightiy1(ix,iy),1)
!              z3 = cry1(rightix1(ix,iy),rightiy1(ix,iy),1)
!              rf = crx1(rightix1(ix,iy),rightiy1(ix,iy),3)
!              z4 = cry1(rightix1(ix,iy),rightiy1(ix,iy),3)
!              nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
!              nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
!              bphatr = (r1-rf+r2-r3)/
!     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!              bphatz = (z1-z4+z2-z3)/
!     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!              bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
!              cosa = nhatr*bphatr+nhatz*bphatz
!              concave = ((r2-r1)*(rf-r1)+(z2-z1)*(z4-z1)).lt.0.0_R8 
!     &            .and. ((r2-r1)*(z4-z1)-(z2-z1)*(rf-r1)).lt.0.0_R8
!              if(concave) write(*,'(a,i4,a,i4,a)') 
!     >         'Cell ',ix,',',iy,' is concave !'
!              if(abs(bdotn).lt.
!     >         bb1(ix,iy,3)*cos(pi*(90.0_R8-min_pitch)/180.0_R8)) then
!                write(*,*) 'Western guard cell',iy,': Pitch too small',
!     .           abs(bdotn)/bb1(ix,iy,3)
cxpb  We first attempt a correction of the guard cell coordinates, if this
cxpb  fails, we change the poloidal field strength
!                cose = sign(sin(min_pitch*pi/180.0_R8),cosa)*
!     .                 bb1(ix,iy,3)/abs(bb1(ix,iy,0))
!                if(abs(cose).le.1.0_R8.and..not.concave) then
!                  beta = acos(cose) - acos(cosa)
!                  cosb = cos(beta)
!                  sinb = sin(beta)
cxpb  Modify r1,z1
cxpb  We rotate point 1 about point 2 by an angle beta
cxpb  and keep the new coordinates as the intersection of the rotated line
cxpb  and the parallel to Bp line between point 4 and the original point 1
!                  numer = (r1 - r2) * cosb - (z1 - z2) * sinb
!                  denom = (r1 - r2) * sinb + (z1 - z2) * cosb
!                  if(denom.ne.0.0_R8) then 
!                    aaa = numer / denom
!                    if (rf.ne.r1) then
!                      bbb = (z4 - z1)/(rf - r1)
!                      r1 = (r2 + (z4-z2)*aaa - rf*bbb*aaa)/
!     .                 (1.0_R8-bbb*aaa)
!                      z1 =  z4 - (rf-r1)*bbb
!                    else
!                      call xertst(aaa.ne.0.0_R8,'unexpected geometry !')
!                      r1 = rf
!                      z1 = z2 + (rf-r2)/aaa
!                    endif
!                  else
!                    if (rf.ne.r1) then
!                      bbb = (z4 - z1)/(rf - r1)
!                      r1 = rf - (z4-z2)/bbb  
!                    else
!                      r1 = rf
!                    endif   
!                    z1 = z2
!                  endif
cxpb  Copy the changes onto crx1, cry1
!                  crx1(rightix1(ix,iy),rightiy1(ix,iy),2) = r1
!                  cry1(rightix1(ix,iy),rightiy1(ix,iy),2) = z1
!                  crx1(ix,iy,3) = r1
!                  cry1(ix,iy,3) = z1
cxpb  ALSO MOVE THE BACK OF THE GUARD CELLS!
!                  crx1(ix,iy,2) = r1 + (r1-rf)/1000.0_R8
!                  cry1(ix,iy,2) = z1 + (z1-z4)/1000.0_R8
!                  if (isInDomain(nx1,ny1,topix1(ix,iy),
!     .                                   topiy1(ix,iy))) then
!                    crx1(topix1(ix,iy),topiy1(ix,iy),1) = r1
!                    cry1(topix1(ix,iy),topiy1(ix,iy),1) = z1
!                    crx1(topix1(ix,iy),topiy1(ix,iy),0) = 
!     .               r1 + (r1-rf)/1000.0_R8
!                    cry1(topix1(ix,iy),topiy1(ix,iy),0) = 
!     .               z1 + (z1-z4)/1000.0_R8
!                    crx1(topix1(rightix1(ix,iy),rightiy1(ix,iy)),
!     .                   topiy1(rightix1(ix,iy),rightiy1(ix,iy)),0) 
!     .               = r1
!                    cry1(topix1(rightix1(ix,iy),rightiy1(ix,iy)),
!     .                   topiy1(rightix1(ix,iy),rightiy1(ix,iy)),0) 
!     .               = z1
!                  endif
cxpb  if just below the top corner, make the corner cell follow
!                  if (iy.eq.ny1-1) then
!                    crx1(ix,ny1,3) = r1 + (r1-r2)/1000.0_R8
!                    cry1(ix,ny1,3) = z1 + (z1-z2)/1000.0_R8
!                    crx1(ix,ny1,2) = r1 + (r1-r3)/1000.0_R8
!                    cry1(ix,ny1,2) = z1 + (z1-z3)/1000.0_R8
!                    crx1(rightix1(ix,ny1),rightiy1(ix,ny1),2) = 
!     .               r1 + (r1-r2)/1000.0_R8
!                    cry1(rightix1(ix,ny1),rightiy1(ix,ny1),2) = 
!     .               z1 + (z1-z2)/1000.0_R8
!                  end if
!                  nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
!                  nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
!                  bphatr = (r1-rf+r2-r3)/
!     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!                  bphatz = (z1-z4+z2-z3)/
!     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)     
!                else  ! rotation correction impossible
!                  bb1(ix,iy,0) = sign(
!     .             cos(pi*(90.0_R8-min_pitch)/180.0_R8)/
!     .             (nhatr*bphatr+nhatz*bphatz), 
!     1             bb1(ix,iy,0))
!                  bb1(ix,iy,3) = sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
!                endif
!                bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
!                write(*,*) 'Corrected pitch :', abs(bdotn)/bb1(ix,iy,3)
!              endif
!            endif
cxpb Eastern edge
!            if (
!     .       .not.isInDomain(nx1,ny1,rightix1(ix,iy),rightiy1(ix,iy)))
!     >       then
!              rf = crx1(leftix1(ix,iy),leftiy1(ix,iy),2)
!              z4 = cry1(leftix1(ix,iy),leftiy1(ix,iy),2)
!              r3 = crx1(leftix1(ix,iy),leftiy1(ix,iy),0)
!              z3 = cry1(leftix1(ix,iy),leftiy1(ix,iy),0)
!              r2 = crx1(leftix1(ix,iy),leftiy1(ix,iy),1)
!              z2 = cry1(leftix1(ix,iy),leftiy1(ix,iy),1)
!              r1 = crx1(leftix1(ix,iy),leftiy1(ix,iy),3)
!              z1 = cry1(leftix1(ix,iy),leftiy1(ix,iy),3)
!              nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
!              nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
!              bphatr = (r1-rf+r2-r3)/
!     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!              bphatz = (z1-z4+z2-z3)/
!     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!              bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
!              cosa = nhatr*bphatr+nhatz*bphatz
!              concave = ((rf-r1)*(r2-r1)+(z4-z1)*(z2-z1)).lt.0.0_R8 
!     &            .and. ((rf-r1)*(z2-z1)-(z4-z1)*(r2-r1)).lt.0.0_R8
!              if(concave) write(*,'(a,i4,a,i4,a)') 
!     >         'Cell ',ix,',',iy,' is concave !'
!              if(abs(bdotn).lt.
!     >         bb1(ix,iy,3)*cos(pi*(90.0_R8-min_pitch)/180.0_R8)) then
!                write(*,*) 'Eastern guard cell',iy,': Pitch too small',
!     .           abs(bdotn)/bb1(ix,iy,3)
cxpb  We first attempt a correction of the guard cell coordinates, if this
cxpb  fails, we change the poloidal field strength
!                cose = sign(sin(min_pitch*pi/180.0_R8),cosa)*
!     .                 bb1(ix,iy,3)/bb1(ix,iy,0)
!                if(abs(cose).le.1.0_R8.and..not.concave) then
!                  beta = acos(cose) - acos(cosa)
!                  cosb = cos(beta)
!                  sinb = sin(beta)
cxpb  Modify r1,z1
cxpb  We rotate point 1 about point 2 by an angle beta
cxpb  and keep the new coordinates as the intersection of the rotated line
cxpb  and the parallel to Bp line between point 4 and the original point 1
!                  numer = (r1 - r2) * cosb - (z1 - z2) * sinb
!                  denom = (r1 - r2) * sinb + (z1 - z2) * cosb
!                  if(denom.ne.0.0_R8) then 
!                    aaa = numer / denom
!                    if (rf.ne.r1) then
!                      bbb = (z4 - z1)/(rf - r1)
!                      r1 = (r2 + (z4-z2)*aaa - rf*bbb*aaa)/
!     .                 (1.0_R8-bbb*aaa)
!                      z1 =  z4 - (rf-r1)*bbb
!                    else
!                      call xertst(aaa.ne.0.0_R8,'unexpected geometry !')
!                      r1 = rf
!                      z1 = z2 + (rf-r2)/aaa
!                    endif
!                  else
!                    if (rf.ne.r1) then
!                      bbb = (z4 - z1)/(rf - r1)
!                      r1 = rf - (z4-z2)/bbb   
!                    else
!                      r1 = rf
!                    endif   
!                    z1 = z2
!                  endif
cxpb  Copy the changes onto crx1, cry1
!                  crx1(leftix1(ix,iy),leftiy1(ix,iy),3) = r1
!                  cry1(leftix1(ix,iy),leftiy1(ix,iy),3) = z1
!                  crx1(ix,iy,2) = r1
!                  cry1(ix,iy,2) = z1
cxpb ALSO MOVE THE BACK OF THE GUARD CELLS!
!                  crx1(ix,iy,3) = r1 + (r1-rf)/1000.0_R8
!                  cry1(ix,iy,3) = z1 + (z1-z4)/1000.0_R8
!                  if (isInDomain(nx1,ny1,topix1(ix,iy),
!     .                                   topiy1(ix,iy))) then
!                    crx1(topix1(ix,iy),topiy1(ix,iy),0) = r1
!                    cry1(topix1(ix,iy),topiy1(ix,iy),0) = z1
!                    crx1(topix1(ix,iy),topiy1(ix,iy),1) = 
!     .               r1 + (r1-rf)/1000.0_R8
!                    cry1(topix1(ix,iy),topiy1(ix,iy),1) = 
!    .               z1 + (z1-z4)/1000.0_R8
!                    crx1(topix1(leftix1(ix,iy),leftiy1(ix,iy)),
!     .                   topiy1(leftix1(ix,iy),leftiy1(ix,iy)),1) 
!     .               = r1
!                    cry1(topix1(leftix1(ix,iy),leftiy1(ix,iy)),
!     .                   topiy1(leftix1(ix,iy),leftiy1(ix,iy)),1) 
!     .               = z1
!                  endif
cxpb  if just below the top corner, make the corner cell follow
!                  if (iy.eq.ny1-1) then
!                    crx1(ix,ny1,2) = r1 + (r1-r2)/1000.0_R8
!                    cry1(ix,ny1,2) = z1 + (z1-z2)/1000.0_R8
!                    crx1(ix,ny1,3) = r1 + (r1-r3)/1000.0_R8
!                    cry1(ix,ny1,3) = z1 + (z1-z3)/1000.0_R8
!                    crx1(leftix1(ix,ny1),leftiy1(ix,ny1),3) = 
!     .               r1 + (r1-r2)/1000.0_R8
!                    cry1(leftix1(ix,ny1),leftiy1(ix,ny1),3) = 
!     .               z1 + (z1-z2)/1000.0_R8
!                  end if
!                  nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
!                  nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
!                  bphatr = (r1-rf+r2-r3)/
!     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
!                  bphatz = (z1-z4+z2-z3)/
!     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)     
!                else  ! rotation correction impossible
!                  bb1(ix,iy,0) = sign(
!     .             cos(pi*(90.0_R8-min_pitch)/180.0_R8)/
!     .             (nhatr*bphatr+nhatz*bphatz), 
!     1             bb1(ix,iy,0))
!                  bb1(ix,iy,3) = sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
!                endif
!                bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
!                write(*,*) 'Corrected pitch :', abs(bdotn)/bb1(ix,iy,3)
!              endif
!            endif
!          enddo
!        enddo

      endif  ! if of min_pithc and isClassicalGrid
!
*   ..test sign of bb3
!      call b2xvsg1(nx1, ny1, bb1(-1,-1,3), 1, 'bb3', '.gt.')
       call b2xvs1(m%nCv,g%cvBb(:,3), 1, 'cvBb3', '.gt.')

!      write(*,*) 'X regions (from rightix<>ix+1)'
!      do iy=ny1,-1,-1
!        zone(-1:nx1)='0'
!        ireg=1
!        do ix=-1,nx1
!          if(isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
!          if(isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
!            zone(ix)='G'
!          else
!            write(zone(ix),'(i1)') ireg
!            if(rightix1(ix,iy).ne.ix+1) ireg=ireg+1
!          endif
!        enddo
!        write(*,'(3276(10a1,1x))') zone
!      enddo
!      write(*,*) 'X regions (from connections)'
!      do iy=ny1,-1,-1
!        ireg=0
!        zone(-1:nx1)='0'
!        do ix=-1,nx1
!          if(isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
!          if(zone(ix).eq.'0') then
!            if(isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
!              zone(ix) = 'G'
!            else
!              ireg=ireg+1
!              write(zone(ix),'(i1)') ireg
!              ixr=rightix1(ix,iy)
!              iyr=rightiy1(ix,iy)
!              do while(ixr.ne.ix.and.ixr.ge.-1.and.ixr.le.nx1.and.
!     &                 iyr.eq.iy)
!                if(.not.isRealCell(cflag(ixr,iyr,CELLFLAG_TYPE))) exit
!                write(zone(ixr),'(i1)') ireg
!                ixr=rightix1(ixr,iy)
!                if (ixr.ge.-1 .and. ixr.le.nx1)
!     >           iyr=rightiy1(ixr,iy)
!              enddo
!              ixl=leftix1(ix,iy)
!              iyl=leftiy1(ix,iy)
!              do while(ixl.ne.ix.and.ixl.ge.-1.and.ixl.le.nx1.and.
!     &                 iyl.eq.iy)
!                if(.not.isRealCell(cflag(ixl,iyl,CELLFLAG_TYPE))) exit
!                write(zone(ixl),'(i1)') ireg
!                ixl=leftix1(ixl,iy)
!                if (ixl.ge.-1 .and. ixl.le.nx1)
!     >           iyl=leftiy1(ixl,iy)
!              enddo
!            endif
!          endif
!        enddo
!        write(*,'(3276(10a1,1x))') zone
!      enddo

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agfs

      end subroutine b2agfs

      subroutine b2gxdr_dummy(local_sonnet,nCv,nFc,nVx,istyle,
     &  wbbc_us,isymm,g,m)
      use b2mod_types
      use b2mod_version
      use b2mod_indirect
     . ,only: GRID_N_CELLFLAGS
      use b2ag_parameters
      use carre_constants
      use b2ag_ghostcells
      use carre_b2ag
      use b2mod_geo
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
      integer nCv, nFc, nVx, istyle, isymm
   
      character*(*) local_sonnet
      type(geometry_ag) :: g
      type(mapping_ag) :: m

      !nnCv, nnFc, nnVx: dimension of the grid read from file
      !nCv, nFc, nVx: dimension of the grid as used by B2
      integer nnCv, nnFc, nnVx, nGhostCell,
     1 iCv, iFc, iVx, i, j, lun, line_length, iostat, idum(0:2)
#ifdef USE_PXFGETENV
      integer lenval, ierror     
#else
      intrinsic get_environment_variable
#endif
      real (kind=R8) ::
     * ubr, ubt, ubp, upsi
      real(kind=R8) :: bb_tmp_Vx(0:nVx-1) ! lengte onder voorbehoud
      real(kind=R8) :: wbbc_us(0:nCv-1,0:3)
      parameter (line_length=256)
      character*256 meshpath
      character*80 format_string(4)
      integer, save :: new_format
      character*12 device
      logical open_file, uop, first_read
      character uline*(line_length)
      external xertst, subini, subend, open_file
      intrinsic index, sqrt

      data new_format/1/
      data format_string/
     1 '(19x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     2 '(19x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     3 '(21x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     4 '(21x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)'/

*   ..subprogram start-up calls
      call subini('b2gxdr_dummy')
*   ..test nCv, nFc, nVx, istyle
      call xertst (0.le.nCv.or.0.le.nFc.or.0.le.nVx, 
     & 'faulty arggument nCv, nFc, nVx')
      call xertst (istyle.eq.-1 .or. istyle.eq.-2,
     & 'faulty input parameter istyle')

      i=-1 
      j=-1
      lun=96
      meshpath = 'data.local/meshes|data/meshes'
#ifdef JET
      device = 'jet'
#else
      device = 'iter'
#endif
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('DEVICE', 0, device, lenval, ierror)
#else
      call get_environment_variable ('DEVICE', device)
#endif
#endif
      if (istyle.eq.-1) then  !xpb  DG-CARRE-format geometry input
        meshpath = trim(meshpath)
     .   //'|modules/Carre/meshes'//trim(device)
        if(open_file(lun,local_sonnet,meshpath)) then
          write(*,*)
     1     "Reading the geometry file assuming the 'new' format"
          rewind lun
          call cfverr(lun,grid_version)
          call cfruin (lun, 3, idum, 'nCv,nFc,nVx')
          nnCv = idum(0)
          nnFc = idum(1)
	    nnVx = idum(2)
          if (grid_version.eq."03.002.000") then
           !ghostcells are not provided, nnCv does not match nCv
           !but will match after guard cells are created

           !.. read geometry file
           rewind lun
           call read_traduitoutb2us(lun,nnCv,nnFc,nnVx,g,m)

           !calculate number of ghostcells
           nGhostCell = nCv - nnCv

           call create_guard_cells
     1         (nnCv,nnFc,nnVx,nCv,nFc,nVx,g,m)

!     Compute magnetic field with the same recipe as carre...
           call carre_b2agbb(nCv,nVx,g,m,wbbc_us,isymm)

!     ..and transfer it to the ghost cells with sanity checks
             do i = 0, nCv-1
                select case (g%cvCflags(i,0))
                case(GRID_GUARD)  ! 
                   g%cvBp(i) = g%cvBb(i,0)
                   g%cvBt(i) = g%cvBb(i,2)
                case(GRID_INTERNAL,GRID_BOUNDARY)   
                   if (abs(g%cvBp(i) - g%cvBb(i,0))  > 1e-7 .or.
     $                 abs(g%cvBt(i) - g%cvBb(i,2))  > 1e-7) then
                      write(*,*) "b2agfs: not maching carre bb 
     $                at", i    
                   endif
                end select
             enddo       

           endif  ! if of the grid_version

        else
           call xerrab('Error in reading geometry file')
        endif     ! if of open_file
      elseif(istyle.eq.-2) then !xpb SONNET-format geometry input
           call xerrab
     $       ('Sonnet format not supported in unstructured
     $         format. Use the structured format.')

      endif       ! if of istyle .eq. -1

      call subend ()
      return

      !last lines with labels 10, 700-720 and 990 are related to Sonnet

      end subroutine b2gxdr_dummy


*----------------------------------------------------------------------
      
      subroutine read_traduitoutb2us(nrid,nnCv,nnFc,nnVx,g,m)
      use b2mod_types
      use b2mod_version
      use b2mod_indirect
     . ,only: GRID_N_CELLFLAGS
      use b2ag_parameters
      use carre_constants
      use b2ag_ghostcells
      use carre_b2ag
      use b2mod_geo
      implicit none
      !input arguments 
      integer nrid, nnCv, nnFc, nnVx
      !output argument
      type(geometry_ag) :: g 
      type(mapping_ag)  :: m  


      !internal
      character idcod*8, id*32
      integer i, ivx(0:nnVx-1), icv(0:nnCv-1), idum(0:2)
      logical streql
      external subini, subend, xertst, xerrab, streql, cfverr, cfruin

*   ..subprogram start-up calls  
      call subini('read_traduitoutb2us')
*   ..test nrid
      call xertst (1.le.nrid, 'faulty argument nrid')
*   ..read grid version
      !rewind nrid
      call cfverr(nrid,grid_version)
      
      if (grid_version.eq."03.002.000") then
        !read nCv, nFc, nVx (without ghostcells)
        call cfruin (nrid, 3, idum, 'nCv,nFc,nVx')
        nnCv = idum(0)
        nnFc = idum(1)
	  nnVx = idum(2)
        !read nCmVx, nCmxFc
        call cfruin (nrid, 2, idum, 'nCmxVx,nCmxFc')
        m%nCmxVx = idum(0)
        m%nCmxFc = idum(1)
        !read nx, ny
        call cfruin (nrid, 2, idum, 'nx,ny')
        m%nx = idum(0)
        m%ny = idum(1)
        call cfruin (nrid, 1, idum, 'isClassicalGrid')
        m%isClassicalGrid = idum(0)
        !read vertex information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('Vx',id)) then
           do i = 0, nnVx-1
           ! vx vxX vxY vxPsi vxBx vxBy vxFfbz
             read(nrid,*,err=91) ivx(i), g%vxX(i), g%vxY(i),
     &      g%vxFpsi(i), g%vxBx(i), g%vxBy(i), g%vxFfbz(i)
           enddo
        endif

        !read cell information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('cv',id)) then
           do i = 0, nnCv-1
           ! cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(0:4) ctype
              read(nrid,*,err=91) icv(i), m%cvVxP(i,0), m%cvVxP(i,1),
     &     g%cvX(i), g%cvY(i), g%cvFpsi(i), g%cvBp(i), g%cvBt(i),
     &     g%cvCflags(i,0), g%cvCflags(i,1), g%cvCflags(i,2),
     &     g%cvCflags(i,3), g%cvCflags(i,4), g%cvCtype(i)
           enddo
        endif
        m%cvFcP = m%cvVxP

       ! look whether these are correct variable allocations

       !read cvVx
       call cfruin(nrid,m%nCmxVx,m%cvVx,'cvVx')
       !read cvFc
       call cfruin(nrid,m%nCmxFc,m%cvFc,'cvFc')
       !read fcVx(:,1)
       call cfruin(nrid,nnFc,m%fcVx(:,0),'fcVx(:,1)')
       !read fcVx(:,2)
       call cfruin(nrid,nnFc,m%fcVx(:,1),'fcVx(:,2)')


      else
        call xerrab
     1    ('error read_traduitoutb2us, version not yet implemented')

      endif

      call subend ()
      return

91    call xerrab('read_traduitoutb2us: error reading geometry file')


      end subroutine read_traduitoutb2us