*-----------------------------------------------------------------------
*.specification

      subroutine b2agfs (nrid,istyle,lblgm,g,m,artificial_slab)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2ag_parameters
      use b2ag_ghostcells
      use carre_constants
      use b2us_geo
      use b2us_map
      implicit none
#include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nrid, istyle
      type(geometry), intent(inout) :: g
      type(mapping), intent(inout) :: m
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
      integer ix, iy, i, ic, inseltop, inselbot, iCv, iVx
      integer Bt_adjust, Bt_reversal, periodic_bc, artificial_slab,
     .  intcellP_new 
      character*256 local_sonnet
      data local_sonnet /'upgrade.geometry'/
      real (kind=R8) ::
     *  xoffset, yoffset, xrescale, yrescale, Bt_rescale, pit_rescale,
     .  min_pitch, geom_match_dist, distval,
     .  r1,r2,r3,rf,z1,z2,z3,z4,numer,denom,aaa,bbb,
     .  nhatr,nhatz,bphatr,bphatz,bdotn,beta,cosa,cosb,cose,sinb
      real (kind=R8) :: cvBt(m%nCv), cvBp(m%nCv)
      integer ireg,errors
      logical found, concave
      data xoffset /0.0_R8/, yoffset /0.0_R8/, Bt_rescale/1.0_R8/,
     .  Bt_adjust/0/, Bt_reversal/0/, xrescale/1.0_R8/, 
     .  yrescale/1.0_R8/, pit_rescale/1.0_R8/, min_pitch/1.0_R8/, 
     .  geom_match_dist/1.0e-6_R8/, intcellP_new/0/

*   ..procedures
      external subini, subend, ipgetc, ipgeti, ipgetr, xertst,
     * b2xvsg, b2gxdr_dummy
      intrinsic abs, sqrt, cos, sin, min, max, sign, acos

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2agfs')
*   ..test nCv, nFc, nVx, nrid
      call xertst (0.le.m%nCv.and.0.le.m%nFc.and.0.le.m%nVx,
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
      call ipgeti ('b2agfs_intcellP_new' , intcellP_new)
      call xertst (intcellP_new.eq.0 .or. intcellP_new.eq.1,
     &  'faulty argument intcellP_new')      


      ! Read in grid input file (as typically produced by traduit)
      ! This routine also ensures the existence of ghost cells
      call b2gxdr_dummy(local_sonnet,istyle,
     &     isymm,g,m,cvBp,cvBt,artificial_slab,intcellP_new)


      ! Set up magnetic field and shift coordinates
      do iCv = 1,m%nCv
        if ( m%cvLbl(iCv) .eq. GRID_UNDEFINED .or. 
     &     m%cvLbl(iCv) .eq. GRID_EXTERNAL .or. 
     &     m%cvLbl(iCv) .eq. GRID_DEAD ) cycle  !isUnusedCell works only for unstructured
        g%cvBb(iCv,0) = abs(cvBp(iCv))
        g%cvBb(iCv,1) = 0.0_R8
        g%cvBb(iCv,2) = cvBt(iCv) !keep the sign of the toroidal field
        g%cvBb(iCv,3) = sqrt( cvBp(iCv)**2 + cvBt(iCv)**2 )
! WIP B2AG_US
        !verts = m%cvVx(m%cvVxP(i,0)-1:m%cvVxP(i,0)-1+m%cvVxP(i,1)-1)-1
        !verts = m%cvVx(m%cvVxP(i,1):m%cvVxP(i,1)+m%cvVxP(i,2))
        ! array with length 2, 3 and 4 need to be made
        ! to distinguish between ghost, tria en quad
          !g%vxX(i) = g%vxX(i)*xrescale+xoffset
          !g%vxY(i) = g%vxY(i)*yrescale+xoffset

        !g%vxFpsi(i) = g%vxFpsi(i)
        !g%vxFfbz(i) = g%vxFfbz(i)
        if (Bt_adjust.eq.1 .and. xoffset.ne.0.0_R8) then
          g%cvBb(iCv,2) = cvBt(iCv)*g%cvX(iCv)/(g%cvX(iCv)+xoffset)
          g%cvBb(iCv,0) = abs(cvBp(iCv))/cvBt(iCv)*g%cvBb(iCv,3)
          g%cvBb(iCv,3) = sqrt(g%cvBb(iCv,0)**2 + g%cvBb(iCv,2)**2)
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
!WIP B2AG_US: SHOULD BE THIS, BUT GIVE ERROR BY ABSENT BB at guards (wbbc)
       !call b2xvs1(nCv,g%cvBb(:,3), 1, 'cvBb3', '.gt.')
       

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

      subroutine b2gxdr_dummy(local_sonnet,istyle,
     &  isymm,g,m,cvBp,cvBt,artificial_slab,intcellP_new)
      use b2mod_types
      use b2mod_version
      use b2ag_parameters
      use carre_constants
      use b2ag_ghostcells
      use b2us_geo
      use b2us_map
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
      integer  istyle, isymm
   
      character*(*) local_sonnet
      type(geometry) :: g
      type(mapping) :: m

      !nnCv, nnFc, nnVx: dimension of the grid read from file
      !nCv, nFc, nVx: dimension of the grid as used by B2
      integer nnCv, nnFc, nnVx, nnFs, nnFt,
     1 iCv, iFc, iVx, i, j, lun, line_length, iostat, idum(0:5),
     2 artificial_slab, intcellP_new
#ifdef USE_PXFGETENV
      integer lenval, ierror     
#else
      intrinsic get_environment_variable
#endif

      real(kind=R8) :: cvBp(m%nCv), cvBt(m%nCv), cvFpsi(m%nCv)
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
      call xertst (0.le.m%nCv.or.0.le.m%nFc.or.0.le.m%nVx, 
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
          call cfruin (lun, 6, idum, 'nCi,nFc,nVx,nCg,nFs,nFt')
          nnCv = idum(0)
          nnFc = idum(1)
	    nnVx = idum(2)
          nnFs = idum(4)
          nnFt = idum(5)
          if (grid_version.eq."03.002.000") then
           !ghostcells are not provided, nnCv does not match nCv
           !but will match after guard cells are created

           !.. read geometry file
           rewind lun
           call read_traduitoutb2us(lun,nnCv,nnFc,nnVx,nnFs,nnFt,g,m,
     1      cvBp,cvBt,cvFpsi)

!      .. Create the guard cells
           call create_guard_cells
     1         (nnCv,nnFc,nnVx,g,m)

!     .. Construct additional connectivity data (fcCv, vxCv, vxFc)
           call calc_add_connectivity(g,m)

!     .. Order the vertices and faces as rightturning
!     .. and consequently check convexity
         call order_vertices(g,m)
         call check_convexity(g,m)


!     .. Construct intcellP and intcellR
           call calc_intcell(g,m,intcellP_new)     

          
!     Compute magnetic field with the same recipe as carre...
            if (artificial_slab.eq.0) then
              call b2agbb(g,m,isymm)

!     ..and transfer it to the ghost cells with sanity checks

             do iCv = 1, m%nCv
                select case (m%cvFcP(iCv,2))
                case(1)  ! guard
                   cvBp(iCv) = g%cvBb(iCv,0)
                   cvBt(iCv) = g%cvBb(iCv,2)
                case(3) ! internal and boundary triangels   
                   if (abs(abs(cvBp(iCv)) - abs(g%cvBb(iCv,0))) 
     $                   > 1e-7 .or.
     $                 abs(abs(cvBt(iCv)) - abs(g%cvBb(iCv,2))) 
     $                   > 1e-7) then
                      write(*,*) "b2agfs: not maching carre bb 
     $                at", iCv , "because is triangle"   
                   endif
                case(4,5) ! internal and boundary triangels 
                   if (abs(abs(cvBp(iCv)) - abs(g%cvBb(iCv,0))) 
     $                   > 1e-7 .or.
     $                 abs(abs(cvBt(iCv)) - abs(g%cvBb(iCv,2))) 
     $                   > 1e-7) then
                      write(*,*) "b2agfs: not maching carre bb 
     $                at", iCv    
                   endif                  
                end select
             enddo  

            else ! for artificial slab constant magnetic field
                g%cvBb(:,0) = cvBp(1)
                g%cvBb(:,1) = 0
                g%cvBb(:,2) = cvBt(1)
                g%cvBb(:,3) = sqrt( cvBp(1)**2 + cvBt(1)**2)
                g%fcBb(:,0) = cvBp(1)
                g%fcBb(:,1) = 0
                g%fcBb(:,2) = cvBt(1)
                g%fcBb(:,3) = sqrt( cvBp(1)**2 + cvBt(1)**2)
            endif        
           else
             call xerrab('Grid version not equal to 03.002.00')
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
      
      subroutine read_traduitoutb2us(nrid,nnCv,nnFc,nnVx,nnFs,nnFt,g,m,
     &  cvBp,cvBt,cvFpsi)
      use b2mod_types
      use b2mod_version
      use b2ag_parameters
      use carre_constants
      use b2ag_ghostcells
      use b2us_map
      use b2us_geo
      implicit none
      !input arguments 
      integer nrid, nnCv, nnFc, nnVx, nnFt, nnFs,
     &     nnCmxVx,nnCmxFc,nnFmxCv,nnVmxCv,nnVmxFc
      !output argument
      type(geometry) :: g 
      type(mapping)  :: m  


      !internal
      character idcod*8, id*32
      integer :: i, iCv, iVx, iFc, iv(nnVx), ic(nnCv), if(nnFc),
     &   idum(0:5) , iFt, iFs, ifcv(nnFt), ifsv(nnFs), n2
      real(kind=R8) :: cvBp(m%nCv) , cvBt(m%nCv), cvFpsi(m%nCv)
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
        call cfruin (nrid, 6, idum, 'nCi,nFc,nVx,nCg,nFs,nFt')
        nnCv = idum(0)
        nnFc = idum(1)
	  nnVx = idum(2)
        nnFs = idum(4)
        nnFt = idum(5)
        !read nCmVx, nCmxFc
        call cfruin (nrid, 5, idum, 
     &    'nCmxVx,nCmxFc,nFmxCv,nVmxCv,nVmxFc')
        nnCmxVx = idum(0)
        nnCmxFc = idum(1)
        nnFmxCv = idum(2)
        nnVmxCv = idum(3)
        nnVmxFc = idum(4)
        
        ! read isClassicalGrid
        call cfruin (nrid, 1, idum, 'isClassicalGrid')
        m%isClassicalGrid = idum(0)
        
        if (m%isClassicalGrid.eq.1) then
        !read nx, ny, nncut
        call cfruin (nrid, 3, idum, 'nx,ny,nncut')
        m%nx = idum(0)
        m%ny = idum(1)
        m%nncut = idum(2)
        endif

        !read vertex information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('Vx',id)) then
           do iVx = 1, nnVx
           ! vx vxX vxY vxPsi vxBx vxBy vxFfbz
             read(nrid,*,err=91) iv(iVx), g%vxX(iVx), g%vxY(iVx),
     &      g%vxFpsi(iVx), g%vxBb(iVx,0), g%vxBb(iVx,1), g%vxFfbz(iVx)
           enddo
        endif

        !read cell information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('cv',id)) then
           do iCv = 1, nnCv
           ! cv cvVxP(:,1) cvVxP(:,2) cvX cvY psi bp bt cflags(1) cvReg cvFt
              read(nrid,*,err=91) ic(iCv), m%cvVxP(iCv,1), 
     &      m%cvVxP(iCv,2), g%cvX(iCv), g%cvY(iCv), cvFpsi(iCv),
     &      cvBp(iCv), cvBt(iCv), m%cvLbl(iCv), m%cvReg(iCv), 
     &      m%cvFt(iCv)
           enddo
        endif
        ! put cvFcP equal to cvVxP
        m%cvFcP = m%cvVxP

       ! look whether these are correct variable allocations

       !read cvVx
       call cfruin(nrid,nnCmxVx,m%cvVx,'cvVx')
       !read cvFc
       call cfruin(nrid,nnCmxFc,m%cvFc,'cvFc')

       !Read face information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('fc',id)) then
           do iFc = 1, nnFc
           ! fc fcVx(:,1) fcVx(:,1) fcLbl(:) fcReg fcAligned
              read(nrid,*,err=91) if(iFc), m%fcVx(iFc,1), 
     &      m%fcVx(iFc,2), m%fcLbl(iFc), m%fcReg(iFc), 
     &      m%fcAligned(iFc)
           enddo
        endif

       !Read flux tube information
        if (nnFt .gt. 0.0_R8) then
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('ft',id)) then
          do iFt = 1, nnFt
          ! ft ftCvP(:,1) ftCvP(:,2) ftFcP(:,1) ftFcP(:,2), ftReg
            read(nrid,*,err=91) ifcv(iFt), m%ftCvP(iFt,1), 
     &      m%ftCvP(iFt,2), m%ftFcP(iFt,1), m%ftFcP(iFt,2), 
     &      m%ftReg(iFt)
          enddo   
        endif

        !read ftCv, length nnCv is an estimate
        call cfruin(nrid,sum(m%ftCvP(nnFt,:))-1,m%ftCv,'ftCv')
        !read ftFc, length nnFc is an estimate
        call cfruin(nrid,sum(m%ftFcP(nnFt,:))-1,m%ftFc,'ftFc')

       !Read flux surface information
        read(nrid,*,err=91) idcod, id
        if (streql('*cf:',idcod).and.streql('fs',id)) then
          do iFs = 1, nnFs
          !fs fsFcP(:,1) fsFcP(:,2) fsPsi
            read(nrid,*,err=91) ifsv(iFs),  m%fsFcP(iFs,1),
     &      m%fsFcP(iFs,2), g%fsPsi(iFs)
          enddo
        endif

        !read fsFc
        call cfruin(nrid,sum(m%fsFcP(nnFs,:))-1,m%fsFc,'fsFc')
        endif

        if (m%isClassicalGrid.eq.1) then
           n2 = (m%nx+2)*(m%ny+2)
           call cfruin(nrid,n2,m%imapCv,'imapCv')
           call cfruin(nrid,n2,m%imapFcx,'imapFcx')
           call cfruin(nrid,n2,m%imapFcy,'imapFcy')
           call cfruin(nrid,n2,m%imapVx,'imapVx')
        endif

      else
        call xerrab
     1    ('error read_traduitoutb2us, version not yet implemented')

      endif

      call subend ()
      return

91    call xerrab('read_traduitoutb2us: error reading geometry file')


      end subroutine read_traduitoutb2us

*-----------------------------------------------------------------------
      
      subroutine order_vertices(g,m)
      use b2mod_types
      use b2us_cellhelper
      use b2us_map
      use b2us_geo
      implicit none 
      !output argument
      type(geometry) :: g 
      type(mapping)  :: m 

      ! internal

      logical :: is_ordered(m%nCv)
      integer :: iCv, i, j, vs(1:2), nv, nf, index(1:2)
      real (kind=R8) ::  vec_start(1:2), vec_face1(1:2),
     &     vec_face2(1:2), sin1, sin2, fc1X, fc1Y, fc2X, fc2Y
      real (kind=R8),allocatable :: vec_vxX(:), vec_vxY(:),
     &     vec_fcX(:), vec_fcY(:) , fcX(:), fcY(:)
      integer, allocatable :: verts(:), faces(:),  new_verts(:),
     &     new_faces(:) , verts_of_faces(:,:), ind_vof(:), vof_h(:)

    


      !computation
      is_ordered = 1
      ! first CheckVertOrder

      do iCv = 1,m%nCv   ! or use upper boundary m%nCv-m%nCg
         nv = m%cvVxP(iCv,2)
         nf = m%cvFcP(iCv,2)
         if (nv.lt.3) cycle
         allocate(verts(nv))
         allocate(faces(nf))
         allocate(new_verts(nv))
         allocate(new_faces(nf))
         allocate(verts_of_faces(nf,2))

         verts = 
     &    m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1) + nv-1 )
         faces = 
     &    m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1) + nf-1 )

      !first check whether vertices and faces are connected according to the correct ordering    
         do i = 1,nv-1
            if (( .not. isVertInFace(verts(i),faces(i),m)).or. 
     &            ( .not. isVertInFace(verts(i+1),faces(i),m))) then
                is_ordered(iCv) = 0
            endif

         enddo
         i = nv
         if (( .not. isVertInFace(verts(i),faces(i),m)) .or. 
     &      ( .not. isVertInFace(verts(1),faces(i),m))) then
                is_ordered(iCv) = 0
         endif

         if (is_ordered(iCv).eq.1) then ! if above is still OK, check the directions
           allocate(vec_vxX(nv))
           allocate(vec_vxY(nv))

          !Calculate vectors between centroids and vertices
          vec_vxX = g%vxX(verts) - g%cvX(iCv)
          vec_vxY = g%vxY(verts) - g%cvY(iCv)

          !Calculate vectors between centoid and face centers
          !face centers
          fcX = 0.5_R8*(g%vxX(m%fcVx(faces,1)) 
     &               + g%vxX(m%fcVx(faces,2)) )
          fcY = 0.5_R8*(g%vxY(m%fcVx(faces,1))
     &               + g%vxY(m%fcVx(faces,2)) )
          !vectors
          vec_fcX = fcX - g%cvX(iCv)
          vec_fcY = fcY - g%cvY(iCv)

          do j =1,m%cvFcP(iCv,2)-1
            !check sinus between vector to vertex (a) and vector to face
            !center (b)
            !calculate angles (sin = |a x b| / norm(a)*norm(b))
            !with |a x b | = ax*by - bx*ay

               sin1 = vec_vxX(j)*vec_fcY(j) 
     &               - vec_fcX(j)*vec_vxY(j)

               if (sin1 .lt. 0.0_R8) then
                  is_ordered(iCv) = 0
               else
                !check sinus between vector to face center (b) and vector
                !to second vertex (a)
                sin2 = vec_fcX(j) *vec_vxY(j+1) 
     &             - vec_vxX(j+1)*vec_fcY(j)
                if (sin2 .lt. 0.0_R8) then
                   is_ordered(iCv)  = 0

                endif
               endif
   
          enddo
          deallocate(vec_vxX)
          deallocate(vec_vxY)

         endif

         ! second ReorderCellConn
         if (is_ordered(iCv).eq.0) then

             !Starting point 
             new_verts(1:nv) = 0
             new_faces(1:nf) = 0
             new_verts(1) = verts(1)

             verts_of_faces = m%fcVx(faces,:)

             ! find the first face
             allocate(ind_vof(2*nf))
             ind_vof(1:nf)      = (/ (i, i = 1,nf) /)
             ind_vof(nf+1:2*nf) =  (/ (i, i=1,nf) /)

             allocate(vof_h(2*nf))
             vof_h(1:nf)      = verts_of_faces(1:nf,1)
             vof_h(nf+1:2*nf) = verts_of_faces(1:nf,2)

             index(1:2) = pack(ind_vof,vof_h == new_verts(1))         
             ![index,~] = find(verts_of_faces == new_verts(1));

             !Define vector from cell center to first vertex
             vec_start(1) = g%vxX(new_verts(1)) - g%cvX(iCv)
             vec_start(2) = g%vxY(new_verts(1)) - g%cvY(iCv)
             
             !Define vector from cell center to first face center
             fc1X = 0.5_R8*( g%vxX(m%fcVx(faces(index(1)),1)) +
     &            g%vxX(m%fcVx(faces(index(1)),2))  )
             fc1Y = 0.5_R8*( g%vxY(m%fcVx(faces(index(1)),1)) +  
     &          g%vxY(m%fcVx(faces(index(1)),2))  )

             vec_face1(1) = fc1X - g%cvX(iCv)
             vec_face1(2) = fc1Y - g%cvY(iCv)

             fc2X = 0.5_R8*( g%vxX(m%fcVx(faces(index(2)),1)) +  
     &           g%vxX(m%fcVx(faces(index(2)),2))  )
             fc2Y = 0.5_R8*( g%vxY(m%fcVx(faces(index(2)),1)) +  
     &           g%vxY(m%fcVx(faces(index(2)),2))  )

             vec_face2(1) = fc2X - g%cvX(iCv)
             vec_face2(2) = fc2Y - g%cvY(iCv)

             !calculate angles
             sin1 = vec_start(1) * vec_face1(2) - 
     &               vec_face1(1) * vec_start(2)

             sin2 = vec_start(1) * vec_face2(2) - 
     &              vec_face2(1) * vec_start(2)

             !Choose the one with positive sign, which is right turning
             if (sin1 .gt. 0.0_R8) then
                 new_faces(1) =  faces(index(1))
             elseif (sin2.gt.0.0_R8) then
                 new_faces(1) =  faces(index(2))
             else
                 call xerrab('order_vertices: no start face assigned')
             endif 

             !Chain the other verts and face accordingly


             do i = 2,m%cvVxP(iCv,2)
                !find the next vertex
                vs = m%fcVx(new_faces(i-1),:)
                if (.not. isMember(vs(1), new_verts,nv)) then
                    new_verts(i) = vs(1)
                else 
                    new_verts(i) = vs(2)                
                endif

                !find the next face
                index(1:2) = pack(ind_vof,vof_h == new_verts(i))
                !index = find(verts_of_faces == new_verts(i))

                if (.not. isMember(faces(index(1)),new_faces,nf)) then
                   new_faces(i) = faces(index(1))
                else
                   new_faces(i) = faces(index(2))
                endif
             enddo

             !Plug new_verts and new_faces in cvVx an cvFc
             m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+m%cvVxP(iCv,2)-1) =
     &               new_verts
             m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1)+m%cvFcP(iCv,2)-1) =
     &                new_faces

            deallocate(ind_vof)
            deallocate(vof_h)                  
       

         endif

         deallocate(verts)
         deallocate(faces)
         deallocate(new_verts)
         deallocate(new_faces)
         deallocate(verts_of_faces)

      enddo

      end subroutine order_vertices

*-----------------------------------------------------------------------

      subroutine check_convexity(g,m)
      use b2mod_types
      use b2us_map
      use b2us_geo
      implicit none 
      !output argument 
      type(mapping)  :: m
      type(geometry)  :: g
      integer iCv, nv, i, j, k
      integer, allocatable ::  vx(:)
      real (kind=R8) :: sin, dx1, dx2, dy1, dy2

      do iCv = 1,m%nCv
         nv = m%cvVxP(iCv,2)
         if (nv.ge.4) then
           allocate(vx(nv))
           vx = m%cvVx(m%cvVxP(iCv,1):m%cvVxP(iCv,1)+nv-1) 
           do i = 1,nv
              j = i-1
              if (j .lt. 1) then
                 j = nv
              endif
              k = i+1
              if (k.gt.nv) then
                 k = 1
              endif

              dx1 = g%vxX(vx(j)) - g%vxX(vx(i))
              dy1 = g%vxY(vx(j)) - g%vxY(vx(i))
              dx2 = g%vxX(vx(k)) - g%vxX(vx(i))
              dy2 = g%vxY(vx(k)) - g%vxY(vx(i)) 

              sin = dx1*dy2 - dy1*dx2

              if (sin.gt.1.0e-7_R8) then
                write(*,*) "b2agfs: cell no convex", iCv
                !call xerrab('Check convexity: cell is not convex') 
              endif           
           enddo
           deallocate(vx)
         endif
      enddo


      end subroutine check_convexity


*-----------------------------------------------------------------------
   
      subroutine calc_intcell(g,m,intcellP_new)

      use b2mod_types
      use b2us_cellhelper
      use b2us_map
      use b2us_geo
      implicit none 
      !output argument
      type(geometry) :: g 
      type(mapping)  :: m

      !internal 
      integer :: iCv, faces4(0:3), i1, i2, i, verts3(0:2), f3(0:2),
     &    v3(0:2) , f5(0:4) , v5(0:4), counter, faces3(0:2), v4(0:3),
     &    v1, v2, intcellP_new
      real (kind=R8) :: dpsi4(0:3), dpsi_max, dpsi3(0:2), intP, 
     &   max_length, fcS5(0:4)
      intrinsic sqrt

      do iCv = 1, m%nCv

         if (m%cvFcP(iCv,2) == 4) then !quad

            faces4 = m%cvFc(m%cvFcP(iCv,1):m%cvFcP(iCv,1)+3)

            if (intcellP_new.eq.1) then
              !NEW METHOD to account for trapezoids
              v4 = m%cvVx(m%cvVxP(iCv,1): m%cvVxP(iCv,1)+3)

              i1 = minloc(g%vxFpsi(v4),1) -1
              i2 = maxloc(g%vxFpsi(v4),1) -1
              dpsi_max = g%vxFpsi(v4(i2)) - g%vxFpsi(v4(i1)) 

              do i = 0,3
               intP = 
     &         abs( g%vxFpsi(m%fcVx(faces4(i),1)) - 
     &             g%vxFpsi(m%fcVx(faces4(i),2)) ) /  dpsi_max 
               if (intP .gt. 0.95_R8) then
                m%intcellP(m%cvFcP(iCv,1) + i) = 1.0_R8
                m%intcellR(m%cvFcP(iCv,1) + i) = 0.0_R8
               elseif (intP .lt. 0.05_R8) then
                m%intcellP(m%cvFcP(iCv,1) + i) = 0.0_R8
                m%intcellR(m%cvFcP(iCv,1) + i) = 1.0_R8
               else
                m%intcellP(m%cvFcP(iCv,1)+i) = intP
                !no possible to assign intR value based on intP
               endif 
         
              enddo

              ! To find correct intR
              !find the most aligned face
               dpsi4 = abs(g%vxFpsi(m%fcVx(faces4,1)) - 
     &         g%vxFpsi(m%fcVx(faces4,2)))
               i1 = minloc(dpsi4(0:3),1)-1
               !find the opposite face index of the most aligned faces
               !as faces are right turned ordered
               i2 = i1 + 2
               if (i2 .gt. 3) then
                 i2 = i2 - 4
               endif 
  

               do i = 0,3
                 if ((i .eq. i1) .or. (i .eq. i2)) then 
                     ! aligned face or face with span in radial direction
                     m%intcellR(m%cvFcP(iCv,1)+i) = 1.0_R8
                 else 
                     ! not aligned faces => poloidal faces
                     m%intcellR(m%cvFcP(iCv,1)+i) = 0.0_R8
                 endif
               enddo
           
            elseif (intcellP_new.eq.0) then

             !OLD-METHOD - not correct for trapezoids
             !find the most aligned face
             dpsi4 = abs(g%vxFpsi(m%fcVx(faces4,1)) - 
     &        g%vxFpsi(m%fcVx(faces4,2)))
             i1 = minloc(dpsi4(0:3),1)-1
             !find the opposite face index of the most aligned faces
             ! as faces are right turned ordered
             ! find the second most aligned face
             counter = -1
             do i = 0,3
               if (i.ne.i1) then
                faces3(counter + 1) = faces4(i)
                counter = counter + 1
               endif
             enddo

             dpsi3 = abs(g%vxFpsi(m%fcVx(faces3,1)) - 
     &       g%vxFpsi(m%fcVx(faces3,2)))
             i2 = minloc(dpsi3(0:2),1)-1  

             if (i2 .ge. i1) then ! compensate for removed face out of the array
              i2 = i2 + 1
             endif           

             do i = 0,3
                 if ((i .eq. i1) .or. (i .eq. i2)) then 
                     ! aligned faces
                     !intcellP = 0 and intcellR = 1
                     m%intcellP(m%cvFcP(iCv,1)+i) = 0.0_R8
                     m%intcellR(m%cvFcP(iCv,1)+i) = 1.0_R8
                 else 
                     ! not aligned faces => poloidal faces
                     m%intcellP(m%cvFcP(iCv,1)+i) = 1.0_R8
                     m%intcellR(m%cvFcP(iCv,1)+i) = 0.0_R8
                 endif
             enddo

            endif 

         elseif (m%cvFcP(iCv,2) == 3) then !tria
             v3 = m%cvVx(m%cvVxP(iCv,1):
     &            m%cvVxP(iCv,1)+2)
             f3 = m%cvFc(m%cvFcP(iCv,1):
     &            m%cvFcP(iCv,1)+2)

             i1 = minloc(g%vxFpsi(v3),1) -1
             i2 = maxloc(g%vxFpsi(v3),1) -1
             dpsi_max = g%vxFpsi(v3(i2)) - g%vxFpsi(v3(i1)) 

             do i = 0,2
             m%intcellP(m%cvFcP(iCv,1) + i) = 
     &        abs( g%vxFpsi(m%fcVx(f3(i),1)) - 
     &             g%vxFpsi(m%fcVx(f3(i),2)) ) / (2.0_R8 * dpsi_max)         
              m%intcellR(m%cvFcP(iCv,1)+i) = 
     &         m%intcellP(m%cvFcP(iCv,1)+i)
             enddo

         elseif (m%cvFcP(iCv,2) == 5) then !pent
            write (*,*) 'intcell only for quasi-quad pents'
             v5 = m%cvVx(m%cvVxP(iCv,1):
     &            m%cvVxP(iCv,1)+4)
             f5 = m%cvFc(m%cvFcP(iCv,1):
     &            m%cvFcP(iCv,1)+4)
             
             i1 = minloc(g%vxFpsi(v5),1) -1
             i2 = maxloc(g%vxFpsi(v5),1) -1
             dpsi_max = g%vxFpsi(v5(i2)) - g%vxFpsi(v5(i1)) 
             max_length = 0.0_R8
             

             do i = 0,4
                v1 = m%fcVx(f5(i),1)
                v2 = m%fcVx(f5(i),2) 
                fcS5(i) =  sqrt( (g%vxX(v1) - g%vxX(v2) )**2  + 
     &                    (g%vxY(v1) - g%vxY(v2) )**2 )  
                if ( ( isAligned(f5(i),m)).and.
     &            (fcS5(i).gt.max_length) ) then
                   max_length = fcS5(i)
                endif
                 
             enddo
             do i = 0,4
               m%intcellP(m%cvFcP(iCv,1) + i) = 
     &        abs( g%vxFpsi(m%fcVx(f5(i),1)) - 
     &             g%vxFpsi(m%fcVx(f5(i),2)) ) / ( dpsi_max)
              if (m%intcellP(m%cvFcP(iCv,1) + i) .gt. 1.0e-1) then
                 m%intcellR(m%cvFcP(iCv,1) + i) = 0.0_R8
              else 
                 m%intcellR(m%cvFcP(iCv,1) + i) = 
     &            1.0_R8* fcS5(i)/max_length
              endif
      
             enddo          


         elseif (m%cvFcP(iCv,2) == 1) then !guard
            !For guard cell both intcellP and intcellR at one

            m%intcellP(m%cvFcP(iCv,1)) = 1.0_R8
            m%intcellR(m%cvFcP(iCv,1)) = 1.0_R8

         else
           call xerrab('calc_intcell: type of cell not implemented')
           
         endif

      enddo
      
   

      end subroutine calc_intcell


*-----------------------------------------------------------------------      
