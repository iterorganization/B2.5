*-----------------------------------------------------------------------
*.specification

      subroutine b2agfs_st (nrid,nx,ny,nx1,ny1,istyle,
     1 cflag,lblgm,periodic_bc,nncut,ccut1)

      use b2mod_geo
     . ,only: geom_match_dist
      use b2mod_geo2
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_ipmain
      use b2mod_connectivity
      use b2mod_indirect
     . ,only: GRID_N_CELLFLAGS
      use b2ag_parameters
      use b2ag_ghostcells
      use carre_constants
      use b2mod_subsys
      use b2mod_dimensions
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nx, ny, nx1, ny1, istyle, ccut1
*   ..output arguments
      integer :: periodic_bc,nncut,cflag(-1:nx1,-1:ny1,GRID_N_CELLFLAGS)
      character*(*) lblgm
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2AGFS_ST obtains a description of the flux surface structure.
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
      integer ix, iy, ixl, ixr, iyl, iyr, inseltop, inselbot
      integer Bt_adjust, Bt_reversal
      character*256 local_sonnet
      data local_sonnet /'upgrade.geometry'/
      real (kind=R8) ::
     *     crxs(-1:nx1,-1:ny1,0:4),
     1     crys(-1:nx1,-1:ny1,0:4),
     2     bp(-1:nx1,-1:ny1),
     3     bt(-1:nx1,-1:ny1),
     4     psi(-1:nx1,-1:ny1,0:4,0:2), ! Dims: 1,2: cell position, 3: center=0/corners=1..4, 4: value, ddx, ddy
     5     ffbz(-1:nx1,-1:ny1,0:3), ! ffbz at cell corners
     6     wbbc(-1:nx1+1,-1:ny1+1,0:3)
      real (kind=R8) ::
     *  xoffset, yoffset, xrescale, yrescale, Bt_rescale, pit_rescale,
     .  min_pitch,
     .  r1,r2,r3,rf,z1,z2,z3,z4,numer,denom,aaa,bbb,
     .  nhatr,nhatz,bphatr,bphatz,bdotn,beta,cosa,cosb,cose,sinb
      integer ireg
      character*1 zone(-1:nx1)
      logical concave
      data xoffset /0.0_R8/, yoffset /0.0_R8/, Bt_rescale/1.0_R8/,
     .  Bt_adjust/0/, Bt_reversal/0/, xrescale/1.0_R8/,
     .  yrescale/1.0_R8/, pit_rescale/1.0_R8/, min_pitch/1.0_R8/

*   ..procedures
      external ipgeti, ipgetr, xertst, b2xvsg, b2gxdr_dummy_st
      intrinsic abs, cos, sin, sign, acos

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2agfs_st')
*   ..test nx, ny, nrid
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
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
      call b2gxdr_dummy_st(local_sonnet,
     & nx1, ny1, istyle, isymm,
     & crxs, crys, psi, bp, bt, ffbz, wbbc, cflag)

      ! Set up magnetic field and shift coordinates
      do  iy = -1, ny1
        do  ix = -1, nx1
          if (isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
          bb1(ix,iy,0)=bp(ix,iy)
          bb1(ix,iy,1)=0.0_R8
          bb1(ix,iy,2)=bt(ix,iy)
          bb1(ix,iy,3)=sqrt(bp(ix,iy)**2 + bt(ix,iy)**2)
          crx1(ix,iy,0:3)=crxs(ix,iy,1:4)*xrescale+xoffset
          cry1(ix,iy,0:3)=crys(ix,iy,1:4)*yrescale+yoffset
          fpsi1(ix,iy,0:3,0:2) = psi(ix,iy,1:4,0:2)
          ffbz1(ix,iy,0:3) = ffbz(ix,iy,0:3)
          wbbc1(ix,iy,0:3) = wbbc(ix,iy,0:3)
          if (ix.eq.nx1 .and. iy.lt.ny1) then
            wbbc1(nx1+1,iy,0:3) = wbbc(nx1+1,iy,0:3)
          else if (ix.lt.nx1 .and. iy.eq.ny1) then
            wbbc1(ix,ny1+1,0:3) = wbbc(ix,ny1+1,0:3)
          else if (ix.eq.nx1 .and. iy.eq.ny1) then
            wbbc1(nx1+1,ny1+1,0:3) = wbbc(nx1+1,ny1+1,0:3)
          end if
          if (Bt_adjust.eq.1) then
            if (isymm.eq.1 .and. xoffset.ne.0.0_R8) then
              bb1(ix,iy,2)=bt(ix,iy)*
     &          crxs(ix,iy,0)/(crxs(ix,iy,0)+xoffset)
              wbbc1(ix,iy,2)=wbbc(ix,iy,2)*
     &          crxs(ix,iy,1)/(crxs(ix,iy,1)+xoffset)
              if (ix.eq.nx1 .and. iy.lt.ny1) then
                wbbc1(nx1+1,iy,2) = wbbc(nx1+1,iy,2)*
     &            crxs(nx1,iy,2)/(crxs(nx1,iy,2)+xoffset)
              else if (ix.lt.nx1 .and. iy.eq.ny1) then
                wbbc1(ix,ny1+1,2) = wbbc(ix,ny1+1,2)*
     &            crxs(ix,ny1,3)/(crxs(ix,ny1,3)+xoffset)
              else if (ix.eq.nx1 .and. iy.eq.ny1) then
                wbbc1(nx1+1,ny1+1,2) = wbbc(nx1+1,ny1+1,2)*
     &            crxs(nx1,ny1,3)/(crxs(nx1,ny1,4)+xoffset)
              end if
              fpsi1(ix,iy,0:3,0)=psi(ix,iy,1:4,0)/crxs(ix,iy,1:4)*
     .                               (crxs(ix,iy,1:4)+xoffset)
            else if (isymm.eq.3 .and. yoffset.ne.0.0_R8) then
              bb1(ix,iy,2)=bt(ix,iy)*
     &          crys(ix,iy,0)/(crys(ix,iy,0)+yoffset)
              wbbc1(ix,iy,2)=wbbc(ix,iy,2)*
     &          crys(ix,iy,1)/(crys(ix,iy,1)+yoffset)
              if (ix.eq.nx1 .and. iy.lt.ny1) then
                wbbc1(nx1+1,iy,2) = wbbc(nx1+1,iy,2)*
     &            crys(nx1,iy,2)/(crys(nx1,iy,2)+yoffset)
              else if (ix.lt.nx1 .and. iy.eq.ny1) then
                wbbc1(ix,ny1+1,2) = wbbc(ix,ny1+1,2)*
     &            crys(ix,ny1,3)/(crys(ix,ny1,3)+yoffset)
              else if (ix.eq.nx1 .and. iy.eq.ny1) then
                wbbc1(nx1+1,ny1+1,2) = wbbc(nx1+1,ny1+1,2)*
     &            crys(nx1,ny1,3)/(crys(nx1,ny1,4)+yoffset)
              end if
              fpsi1(ix,iy,0:3,0)=psi(ix,iy,1:4,0)/crys(ix,iy,1:4)*
     .                               (crys(ix,iy,1:4)+yoffset)
            end if
            bb1(ix,iy,0)=bp(ix,iy)/bt(ix,iy)*bb1(ix,iy,2)
            bb1(ix,iy,3)=sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
            wbbc1(ix,iy,0)=wbbc(ix,iy,0)/wbbc(ix,iy,2)*wbbc1(ix,iy,2)
            wbbc1(ix,iy,3)=sqrt(wbbc1(ix,iy,0)**2+wbbc1(ix,iy,2)**2)
          endif
        enddo
      enddo
      bb1(:,:,0)=bb1(:,:,0)*Bt_rescale*pit_rescale
      bb1(:,:,1)=bb1(:,:,1)*Bt_rescale*pit_rescale
      wbbc1(:,:,0)=wbbc1(:,:,0)*Bt_rescale*pit_rescale
      wbbc1(:,:,1)=wbbc1(:,:,1)*Bt_rescale*pit_rescale
      fpsi1(:,:,:,:)=fpsi1(:,:,:,:)*Bt_rescale*pit_rescale
      if (Bt_reversal.eq.1) then
        bb1(:,:,2)=-bb1(:,:,2)*Bt_rescale
        wbbc1(:,:,2)=-wbbc1(:,:,2)*Bt_rescale
        if (isymm.eq.0) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale
        else if (isymm.eq.1.or.isymm.eq.2) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale*xrescale
        else if (isymm.eq.3.or.isymm.eq.4) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale*yrescale
        end if
      else
        bb1(:,:,2)=bb1(:,:,2)*Bt_rescale
        wbbc1(:,:,2)=wbbc1(:,:,2)*Bt_rescale
        if (isymm.eq.0) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale
        else if (isymm.eq.1.or.isymm.eq.2) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale*xrescale
        else if (isymm.eq.3.or.isymm.eq.4) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale*yrescale
        end if
      endif
      bb1(:,:,3)=sqrt(bb1(:,:,0)**2+bb1(:,:,1)**2+bb1(:,:,2)**2)
      wbbc1(:,:,3)=sqrt(wbbc1(:,:,0)**2+wbbc1(:,:,1)**2+wbbc1(:,:,2)**2)

 9999 continue

      call init_connectivity( nx1, ny1, cflag,
     & periodic_bc, nncut, ccut1, inseltop, inselbot, istyle )

      call test_connectivity(nx1,ny1,crx1,cry1,cflag,
     & leftix1,leftiy1,rightix1,rightiy1,
     & topix1,topiy1,bottomix1,bottomiy1 )

cxpb For classical grids only
cxpb Check the pitch angle at the plates and make sure it is at least min_pitch
cxpb Do not correct corner cells
      if (min_pitch.gt.0.0_R8 .and. isClassicalGrid(cflag)) then
        do ix = -1, nx1
          do iy = 0, ny1-1
!             Skip unused cells
             if (isUnusedCell(cflag(ix, iy, CELLFLAG_TYPE))) cycle
cxpb Western edge
            if (.not.isInDomain(nx1,ny1,leftix1(ix,iy),leftiy1(ix,iy)))
     >       then
              r1 = crx1(rightix1(ix,iy),rightiy1(ix,iy),2)
              z1 = cry1(rightix1(ix,iy),rightiy1(ix,iy),2)
              r2 = crx1(rightix1(ix,iy),rightiy1(ix,iy),0)
              z2 = cry1(rightix1(ix,iy),rightiy1(ix,iy),0)
              r3 = crx1(rightix1(ix,iy),rightiy1(ix,iy),1)
              z3 = cry1(rightix1(ix,iy),rightiy1(ix,iy),1)
              rf = crx1(rightix1(ix,iy),rightiy1(ix,iy),3)
              z4 = cry1(rightix1(ix,iy),rightiy1(ix,iy),3)
              nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
              nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
              bphatr = (r1-rf+r2-r3)/
     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
              bphatz = (z1-z4+z2-z3)/
     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
              bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
              cosa = nhatr*bphatr+nhatz*bphatz
              concave = ((r2-r1)*(rf-r1)+(z2-z1)*(z4-z1)).lt.0.0_R8
     &            .and. ((r2-r1)*(z4-z1)-(z2-z1)*(rf-r1)).lt.0.0_R8
              if(concave) write(*,'(a,i4,a,i4,a)')
     >         'Cell ',ix,',',iy,' is concave !'
              if(abs(bdotn).lt.
     >         bb1(ix,iy,3)*cos(pi*(90.0_R8-min_pitch)/180.0_R8)) then
                write(*,*) 'Western guard cell',iy,': Pitch too small',
     .           abs(bdotn)/bb1(ix,iy,3)
cxpb  We first attempt a correction of the guard cell coordinates, if this
cxpb  fails, we change the poloidal field strength
                cose = sign(sin(min_pitch*pi/180.0_R8),cosa)*
     .                 bb1(ix,iy,3)/abs(bb1(ix,iy,0))
                if(abs(cose).le.1.0_R8.and..not.concave) then
                  beta = acos(cose) - acos(cosa)
                  cosb = cos(beta)
                  sinb = sin(beta)
cxpb  Modify r1,z1
cxpb  We rotate point 1 about point 2 by an angle beta
cxpb  and keep the new coordinates as the intersection of the rotated line
cxpb  and the parallel to Bp line between point 4 and the original point 1
                  numer = (r1 - r2) * cosb - (z1 - z2) * sinb
                  denom = (r1 - r2) * sinb + (z1 - z2) * cosb
                  if(denom.ne.0.0_R8) then
                    aaa = numer / denom
                    if (rf.ne.r1) then
                      bbb = (z4 - z1)/(rf - r1)
                      r1 = (r2 + (z4-z2)*aaa - rf*bbb*aaa)/
     .                 (1.0_R8-bbb*aaa)
                      z1 =  z4 - (rf-r1)*bbb
                    else
                      call xertst(aaa.ne.0.0_R8,'unexpected geometry !')
                      r1 = rf
                      z1 = z2 + (rf-r2)/aaa
                    endif
                  else
                    if (rf.ne.r1) then
                      bbb = (z4 - z1)/(rf - r1)
                      r1 = rf - (z4-z2)/bbb
                    else
                      r1 = rf
                    endif
                    z1 = z2
                  endif
cxpb  Copy the changes onto crx1, cry1
                  crx1(rightix1(ix,iy),rightiy1(ix,iy),2) = r1
                  cry1(rightix1(ix,iy),rightiy1(ix,iy),2) = z1
                  crx1(ix,iy,3) = r1
                  cry1(ix,iy,3) = z1
cxpb  ALSO MOVE THE BACK OF THE GUARD CELLS!
                  crx1(ix,iy,2) = r1 + (r1-rf)/1000.0_R8
                  cry1(ix,iy,2) = z1 + (z1-z4)/1000.0_R8
                  if (isInDomain(nx1,ny1,topix1(ix,iy),
     .                                   topiy1(ix,iy))) then
                    crx1(topix1(ix,iy),topiy1(ix,iy),1) = r1
                    cry1(topix1(ix,iy),topiy1(ix,iy),1) = z1
                    crx1(topix1(ix,iy),topiy1(ix,iy),0) =
     .               r1 + (r1-rf)/1000.0_R8
                    cry1(topix1(ix,iy),topiy1(ix,iy),0) =
     .               z1 + (z1-z4)/1000.0_R8
                    crx1(topix1(rightix1(ix,iy),rightiy1(ix,iy)),
     .                   topiy1(rightix1(ix,iy),rightiy1(ix,iy)),0)
     .               = r1
                    cry1(topix1(rightix1(ix,iy),rightiy1(ix,iy)),
     .                   topiy1(rightix1(ix,iy),rightiy1(ix,iy)),0)
     .               = z1
                  endif
cxpb  if just below the top corner, make the corner cell follow
                  if (iy.eq.ny1-1) then
                    crx1(ix,ny1,3) = r1 + (r1-r2)/1000.0_R8
                    cry1(ix,ny1,3) = z1 + (z1-z2)/1000.0_R8
                    crx1(ix,ny1,2) = r1 + (r1-r3)/1000.0_R8
                    cry1(ix,ny1,2) = z1 + (z1-z3)/1000.0_R8
                    crx1(rightix1(ix,ny1),rightiy1(ix,ny1),2) =
     .               r1 + (r1-r2)/1000.0_R8
                    cry1(rightix1(ix,ny1),rightiy1(ix,ny1),2) =
     .               z1 + (z1-z2)/1000.0_R8
                  end if
                  nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
                  nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
                  bphatr = (r1-rf+r2-r3)/
     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
                  bphatz = (z1-z4+z2-z3)/
     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
                else  ! rotation correction impossible
                  bb1(ix,iy,0) = sign(
     .             cos(pi*(90.0_R8-min_pitch)/180.0_R8)/
     .             (nhatr*bphatr+nhatz*bphatz),
     1             bb1(ix,iy,0))
                  bb1(ix,iy,3) = sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
                endif
                bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
                write(*,*) 'Corrected pitch :', abs(bdotn)/bb1(ix,iy,3)
              endif
            endif
cxpb Eastern edge
            if (
     .       .not.isInDomain(nx1,ny1,rightix1(ix,iy),rightiy1(ix,iy)))
     >       then
              rf = crx1(leftix1(ix,iy),leftiy1(ix,iy),2)
              z4 = cry1(leftix1(ix,iy),leftiy1(ix,iy),2)
              r3 = crx1(leftix1(ix,iy),leftiy1(ix,iy),0)
              z3 = cry1(leftix1(ix,iy),leftiy1(ix,iy),0)
              r2 = crx1(leftix1(ix,iy),leftiy1(ix,iy),1)
              z2 = cry1(leftix1(ix,iy),leftiy1(ix,iy),1)
              r1 = crx1(leftix1(ix,iy),leftiy1(ix,iy),3)
              z1 = cry1(leftix1(ix,iy),leftiy1(ix,iy),3)
              nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
              nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
              bphatr = (r1-rf+r2-r3)/
     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
              bphatz = (z1-z4+z2-z3)/
     .         sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
              bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
              cosa = nhatr*bphatr+nhatz*bphatz
              concave = ((rf-r1)*(r2-r1)+(z4-z1)*(z2-z1)).lt.0.0_R8
     &            .and. ((rf-r1)*(z2-z1)-(z4-z1)*(r2-r1)).lt.0.0_R8
              if(concave) write(*,'(a,i4,a,i4,a)')
     >         'Cell ',ix,',',iy,' is concave !'
              if(abs(bdotn).lt.
     >         bb1(ix,iy,3)*cos(pi*(90.0_R8-min_pitch)/180.0_R8)) then
                write(*,*) 'Eastern guard cell',iy,': Pitch too small',
     .           abs(bdotn)/bb1(ix,iy,3)
cxpb  We first attempt a correction of the guard cell coordinates, if this
cxpb  fails, we change the poloidal field strength
                cose = sign(sin(min_pitch*pi/180.0_R8),cosa)*
     .                 bb1(ix,iy,3)/bb1(ix,iy,0)
                if(abs(cose).le.1.0_R8.and..not.concave) then
                  beta = acos(cose) - acos(cosa)
                  cosb = cos(beta)
                  sinb = sin(beta)
cxpb  Modify r1,z1
cxpb  We rotate point 1 about point 2 by an angle beta
cxpb  and keep the new coordinates as the intersection of the rotated line
cxpb  and the parallel to Bp line between point 4 and the original point 1
                  numer = (r1 - r2) * cosb - (z1 - z2) * sinb
                  denom = (r1 - r2) * sinb + (z1 - z2) * cosb
                  if(denom.ne.0.0_R8) then
                    aaa = numer / denom
                    if (rf.ne.r1) then
                      bbb = (z4 - z1)/(rf - r1)
                      r1 = (r2 + (z4-z2)*aaa - rf*bbb*aaa)/
     .                 (1.0_R8-bbb*aaa)
                      z1 =  z4 - (rf-r1)*bbb
                    else
                      call xertst(aaa.ne.0.0_R8,'unexpected geometry !')
                      r1 = rf
                      z1 = z2 + (rf-r2)/aaa
                    endif
                  else
                    if (rf.ne.r1) then
                      bbb = (z4 - z1)/(rf - r1)
                      r1 = rf - (z4-z2)/bbb
                    else
                      r1 = rf
                    endif
                    z1 = z2
                  endif
cxpb  Copy the changes onto crx1, cry1
                  crx1(leftix1(ix,iy),leftiy1(ix,iy),3) = r1
                  cry1(leftix1(ix,iy),leftiy1(ix,iy),3) = z1
                  crx1(ix,iy,2) = r1
                  cry1(ix,iy,2) = z1
cxpb ALSO MOVE THE BACK OF THE GUARD CELLS!
                  crx1(ix,iy,3) = r1 + (r1-rf)/1000.0_R8
                  cry1(ix,iy,3) = z1 + (z1-z4)/1000.0_R8
                  if (isInDomain(nx1,ny1,topix1(ix,iy),
     .                                   topiy1(ix,iy))) then
                    crx1(topix1(ix,iy),topiy1(ix,iy),0) = r1
                    cry1(topix1(ix,iy),topiy1(ix,iy),0) = z1
                    crx1(topix1(ix,iy),topiy1(ix,iy),1) =
     .               r1 + (r1-rf)/1000.0_R8
                    cry1(topix1(ix,iy),topiy1(ix,iy),1) =
     .               z1 + (z1-z4)/1000.0_R8
                    crx1(topix1(leftix1(ix,iy),leftiy1(ix,iy)),
     .                   topiy1(leftix1(ix,iy),leftiy1(ix,iy)),1)
     .               = r1
                    cry1(topix1(leftix1(ix,iy),leftiy1(ix,iy)),
     .                   topiy1(leftix1(ix,iy),leftiy1(ix,iy)),1)
     .               = z1
                  endif
cxpb  if just below the top corner, make the corner cell follow
                  if (iy.eq.ny1-1) then
                    crx1(ix,ny1,2) = r1 + (r1-r2)/1000.0_R8
                    cry1(ix,ny1,2) = z1 + (z1-z2)/1000.0_R8
                    crx1(ix,ny1,3) = r1 + (r1-r3)/1000.0_R8
                    cry1(ix,ny1,3) = z1 + (z1-z3)/1000.0_R8
                    crx1(leftix1(ix,ny1),leftiy1(ix,ny1),3) =
     .               r1 + (r1-r2)/1000.0_R8
                    cry1(leftix1(ix,ny1),leftiy1(ix,ny1),3) =
     .               z1 + (z1-z2)/1000.0_R8
                  end if
                  nhatr = (z1-z2)/sqrt((z1-z2)**2+(r2-r1)**2)
                  nhatz = (r2-r1)/sqrt((z1-z2)**2+(r2-r1)**2)
                  bphatr = (r1-rf+r2-r3)/
     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
                  bphatz = (z1-z4+z2-z3)/
     .             sqrt((r1-rf+r2-r3)**2+(z1-z4+z2-z3)**2)
                else  ! rotation correction impossible
                  bb1(ix,iy,0) = sign(
     .             cos(pi*(90.0_R8-min_pitch)/180.0_R8)/
     .             (nhatr*bphatr+nhatz*bphatz),
     1             bb1(ix,iy,0))
                  bb1(ix,iy,3) = sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
                endif
                bdotn = (nhatr*bphatr+nhatz*bphatz)*bb1(ix,iy,0)
                write(*,*) 'Corrected pitch :', abs(bdotn)/bb1(ix,iy,3)
              endif
            endif
          enddo
        enddo
      endif

*   ..test sign of bb3
      call b2xvsg1(nx1, ny1, bb1(-1,-1,3), 1, 'bb3', '.gt.')

      write(*,*) 'X regions (from rightix<>ix+1)'
      do iy=ny1,-1,-1
        zone(-1:nx1)='0'
        ireg=1
        do ix=-1,nx1
          if(isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
          if(isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
            zone(ix)='G'
          else
            write(zone(ix),'(i1)') ireg
            if(rightix1(ix,iy).ne.ix+1) ireg=ireg+1
          endif
        enddo
        write(*,'(3276(10a1,1x))') zone
      enddo
      write(*,*) 'X regions (from connections)'
      do iy=ny1,-1,-1
        ireg=0
        zone(-1:nx1)='0'
        do ix=-1,nx1
          if(isUnusedCell(cflag(ix,iy,CELLFLAG_TYPE))) cycle
          if(zone(ix).eq.'0') then
            if(isGhostCell(cflag(ix,iy,CELLFLAG_TYPE))) then
              zone(ix) = 'G'
            else
              ireg=ireg+1
              write(zone(ix),'(i1)') ireg
              ixr=rightix1(ix,iy)
              iyr=rightiy1(ix,iy)
              do while(ixr.ne.ix.and.ixr.ge.-1.and.ixr.le.nx1.and.
     &                 iyr.eq.iy)
                if(.not.isRealCell(cflag(ixr,iyr,CELLFLAG_TYPE))) exit
                write(zone(ixr),'(i1)') ireg
                ixr=rightix1(ixr,iy)
                if (ixr.ge.-1 .and. ixr.le.nx1)
     >           iyr=rightiy1(ixr,iy)
              enddo
              ixl=leftix1(ix,iy)
              iyl=leftiy1(ix,iy)
              do while(ixl.ne.ix.and.ixl.ge.-1.and.ixl.le.nx1.and.
     &                 iyl.eq.iy)
                if(.not.isRealCell(cflag(ixl,iyl,CELLFLAG_TYPE))) exit
                write(zone(ixl),'(i1)') ireg
                ixl=leftix1(ixl,iy)
                if (ixl.ge.-1 .and. ixl.le.nx1)
     >           iyl=leftiy1(ixl,iy)
              enddo
            endif
          endif
        enddo
        write(*,'(3276(10a1,1x))') zone
      enddo

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agfs_st

      contains

      logical function match(ix1,iy1,ix2,iy2)
      use b2mod_geo2
      implicit none
      integer, intent(in) :: ix1, iy1, ix2, iy2
      match=
     1 (distc1(ix1,iy1,1,ix2,iy2,0)+distc1(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist
      return
      end function match

      end subroutine b2agfs_st
c
      subroutine b2gxdr_dummy_st(local_sonnet,
     & nx, ny, istyle, isymm,
     & crxs, crys, psi, bp, bt, ffbz, wbbc, cflag)
      use b2mod_types
      use b2mod_version
      use b2mod_indirect
     . ,only: GRID_N_CELLFLAGS
      use b2ag_parameters
      use carre_constants
      use b2ag_ghostcells
      use carre_b2ag
      use b2mod_subsys
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
      integer nx, ny, istyle, isymm
      integer :: cflag(0:nx+1,0:ny+1,GRID_N_CELLFLAGS)
      real (kind=R8) ::
     * crxs(0:nx+1,0:ny+1,0:4),
     1 crys(0:nx+1,0:ny+1,0:4),
     2 bp(0:nx+1,0:ny+1),bt(0:nx+1,0:ny+1),
     4 psi(0:nx+1,0:ny+1,0:4,0:2),ffbz(0:nx+1,0:ny+1,0:3),
     6 wbbc(0:nx+2,0:ny+2,0:3)
      character*(*) local_sonnet

      logical no_flux_data            ! IYS 21.02.2020
      ! nnx, nny: dimension of the grid as read from the file
      ! nx, ny: dimension of the grid as used by B2
      ! NOTE: the nx, ny here corresponds to nx1, ny1 in b2agfs_st

      integer nnx, nny,
     1 i, j, ix, iy, ixd, iyd, lun, line_length, mx, my, iostat
#ifndef NO_GETENV
      character*12 device_env
      logical streql
      external streql
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
      real (kind=R8) ::
     *  ubr, ubt, ubp, upsi
      real(kind=R8) :: bb_tmp(0:nx+1,0:ny+1,0:3)
      parameter (line_length=256)
      character*256 meshpath
      character*80 format_string(4)
      integer, save :: new_format
      character*12 device
      character*16 username
      logical open_file, uop, first_read, ex
      character uline*(line_length)
      character*16 usrnam
      external xertst, open_file, usrnam
      intrinsic index, sqrt
      real (kind=R8) :: t0

      data new_format/1/
      data format_string/
     1 '(19x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     2 '(19x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     3 '(21x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     4 '(21x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)'/

*   ..subprogram start-up calls
      call subini('b2gxdr_dummy_st')
*   ..test nx, ny, istyle
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (istyle.eq.-1 .or. istyle.eq.-2,
     & 'faulty input parameter istyle')

      i=-1
      j=-1
      lun=96
      meshpath = 'data.local/meshes'
#ifdef JET
      device = 'jet'
#else
      device = 'iter'
#endif
#ifndef NO_GETENV
      device_env = ' '
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
#else
      call get_environment_variable('DEVICE',
     . status=ierror, length=lenval)
      if (ierror.eq.0)
     > call get_environment_variable('DEVICE', value=device_env)
#endif
      if (.not.streql(device_env,' ')) device = device_env
#endif
      username = usrnam()
      if(istyle.eq.-1) then        !xpb  DG-CARRE-format geometry input
        meshpath = trim(meshpath)
     .   //'|modules/Carre2/meshes/'//trim(device)
     .   //'|modules/Carre2/meshes/'//trim(device)//'/'//trim(username)
        call find_file(local_sonnet,ex)
        if(open_file(lun,local_sonnet,meshpath)) then
          write(*,*)
     1     "Reading the geometry file assuming the Carre format"
          call cfverr(lun,grid_version)
          read(lun,*) nnx,nny
          if (grid_version.lt."01.001.028") then
             ! Standard (Pre-2012) DG/B2.5-format grid input
             ! Ghost cells are given, nnx,nny - nx,ny have to be the same

             if(nnx.ne.nx.or.nny.ne.ny) then
               write(*,'(a,4i5)')
     1       'b2gxdr_dummy_st: nx,ny code and data disagree',
     2       nx,ny,nnx,nny
               call xerrab ('Mismatch nx/ny and nnx/nny')
             endif
             do iy=0,ny+1
               do ix=0,nx+1
                 read(lun,*,err=10,iostat=iostat) ixd,iyd,
     1            (crxs(ix,iy,i),crys(ix,iy,i),i=0,4),
     2            bp(ix,iy),bt(ix,iy)
               enddo
             enddo
             ! Set fields not provided by input file
             psi = 0.0_R8
             ! We do not get any cell type information from the old style grid input file
             ! Initially, set all to internal. init_connectivity will then properly mark
             ! the ghost cells
             cflag(:,:,1) = GRID_INTERNAL
             cflag(:,:,2:CARREOUT_NCELLFLAGS) = GRID_UNDEFINED
          else
             ! 2011/12 extended grid input with object classification and no ghost cells
             ! nnx, nny do not match nx, ny here, but will after the guard cells are created
             do iy=0,nny-1
                do ix=0,nnx-1
                   read(lun,*,err=10,iostat=iostat) ixd,iyd,
     1                  crxs(ix,iy,0),crys(ix,iy,0),psi(ix,iy,0,0),
     1                  (crxs(ix,iy,i),crys(ix,iy,i),
     $                   psi(ix,iy,i,0),psi(ix,iy,i,1),
     $                   psi(ix,iy,i,2),ffbz(ix,iy,i-1),i=1,4),
     3                  bp(ix,iy),bt(ix,iy),
     4                  (cflag(ix,iy,i),i=1,GRID_N_CELLFLAGS)
                enddo
             enddo

             call create_guard_cells_st(nnx, nny, nx, ny, niso, nxiso,
     $            crxs, crys, psi, bp, bt, ffbz, cflag, .false.)

!     Compute magnetic field with the same recipe as carre...
             call carre_b2agbb(nx,ny,cflag, 
     .            psi(:,:,1:4,0),ffbz,wbbc,bb_tmp, 
     &            crxs(:,:,1:4),crys(:,:,1:4),
     $            psi(:,:,1:4,1),psi(:,:,1:4,2),isymm )

!     ...and transfer it to ghost cells, with sanity checks
             do iy=0,ny+1
                do ix=0,nx+1
                   select case( cflag(ix,iy,CELLFLAG_TYPE) )
                   case(GRID_GUARD)
                      bp(ix,iy) = bb_tmp(ix,iy,0)
                      bt(ix,iy) = bb_tmp(ix,iy,2)
                   case(GRID_INTERNAL, GRID_BOUNDARY)
                      !tmp: less strict test, until computation in carre & b2ag streamlined
                      if (abs(bp(ix,iy) - bb_tmp(ix,iy,0)) > 1e-3 .or.
     $                    abs(bt(ix,iy) - bb_tmp(ix,iy,2)) > 1e-3) then
                         write (*,'(a,2i5,4e15.8)')
     $                    "b2agfs_st: not matching carre bb at",
     $                        ix, iy, bp(ix,iy), bt(ix,iy),
     $                                bb_tmp(ix,iy,0), bb_tmp(ix,iy,2)
                      endif
                   end select
                enddo
             enddo

          end if
          close(lun)
! If the toroidal flux function is not provided, compute it
          if (maxval(ffbz).eq.0.0_R8 .and. minval(ffbz).eq.0.0_R8) then
            if (isymm.eq.1 .or. isymm.eq.3) then
              write(*,*) 'Toroidal flux function missing!'
              write(*,*) 'Trying to fix this automatically'
            end if
            t0 = 0.0_R8
            do iy=0,ny+1
              do ix=0,nx+1
                if (isymm.eq.1) then
                  t0 = t0 + bt(ix,iy) * crxs(ix,iy,4)
                else if (isymm.eq.3) then
                  t0 = t0 + bt(ix,iy) * crys(ix,iy,4)
                end if
              end do
            end do
            ffbz = 2.0_R8 * pi * t0 / (nx + 2) / (ny + 2)
          endif
        else
          call xerrab('Error in reading geometry file')
        endif
      elseif(istyle.eq.-2) then   !xpb  SONNET-format geometry input
        meshpath = trim(meshpath)
     .   //'|modules/DivGeo/device/'//trim(device)
     .   //'|modules/DivGeo/device/'//trim(device)//'/'//trim(username)
        if(open_file(lun,local_sonnet,meshpath)) then
          write(*,*)
     1     "Reading the geometry file assuming the Sonnet format"
          write(*,'(a,i3,2a)') '==== b2gxdr_dummy_st: read unit',lun,
     1     '  file  ',trim(local_sonnet)
   20     read (lun,'(a)', end=990) uline
          if(index(uline,'R*Btor').ne.0) then
            read(uline(index(uline,('R*Btor'))+8:line_length),
     ,       '(e30.0)') ubr
            write(*,*) ' R*Btor = ',ubr
            go to 20
          elseif(index(uline,'===============================').eq.0)
     >     then
            go to 20
          end if
c----------------------------------------------------------------------
          if(ubr.eq.0.0_R8.and.(isymm.eq.1.or.isymm.eq.3)) then
            write(*,*) 'b2gxdr_dummy_st: the value of R*Btor undefined.'
            write(*,*)
     1       'It must be defined in the header of the grid file '
            call xerrab 
     1       ('b2gxdr_dummy_st: the value of R*Btor undefined')
          end if
c----------------------------------------------------------------------
c*** Start reading the grid data
c
          write(*,*) 'Reading the grid data ...'
          first_read = .true.
          mx=0
          my=0
   30     uop=.false.
          if(first_read)
     &         write(*,*) 'Trying with format_string option ',
     &         new_format
          read (lun,format_string(new_format),end=40,err=31)
     ,       i,j,crxs(i,j,3),crys(i,j,3),crxs(i,j,4),crys(i,j,4)
          go to 32
   31     continue
          if (first_read) then
            backspace(lun)
            new_format = new_format+1
            if (new_format.le.4) then
              go to 30
            else
              call xerrab ('Unrecognized format in geometry file')
            end if
          else
            call xerrab ('Unrecognized format in geometry file')
          end if
   32     continue
          first_read = .false.
#ifdef DBG
          uop = i.le.5 .and. j.le.1
#endif
          if(uop) then
            print *,'*3* ',i,j,crxs(i,j,3),crys(i,j,3)
            print *,'*4* ',i,j,crxs(i,j,4),crys(i,j,4)
          end if
          if(i.gt.mx) then
            mx=i
            if(mx.gt.nx+1) then
              write(*,*) ' mx = ',mx,' nx = ',nx
              call xerrab ('b2gxdr_dummy_st: mx > nx+1')
            endif
          end if
          if(j.gt.my) then
            my=j
            if(my.gt.ny+1) then
              write(*,*) ' my = ',my,' ny = ',ny
              call xerrab ('b2gxdr_dummy_st: my > ny+1')
            endif
          end if
          read (lun,710,end=990) upsi,crxs(i,j,0),crys(i,j,0)
          read (lun,'(a)',end=990) uline
          if (uline(32:32).eq.'(') then
            read (uline,720) crxs(i,j,1),crys(i,j,1),
     ,        crxs(i,j,2),crys(i,j,2)
          else if (uline(30:30).eq.'(') then
            read (uline,721) crxs(i,j,1),crys(i,j,1),
     ,        crxs(i,j,2),crys(i,j,2)
          end if
          read (lun,'()')
          if(uop) then
            print *,'*0* ',i,j,crxs(i,j,0),crys(i,j,0)
            print *,'*1* ',i,j,crxs(i,j,1),crys(i,j,1)
            print *,'*2* ',i,j,crxs(i,j,2),crys(i,j,2)
          end if
          if(isymm.ne.2.and.isymm.ne.4) then
            ubt=ubr/crxs(i,j,0)
            ubp=ubt*upsi/sqrt(1.0_R8-upsi**2)
            bp(i,j)=ubp
            bt(i,j)=ubt
          endif
          go to 30
c----------------------------------------------------------------------
   40     continue
c      print *,'Finished reading. nx,ny,mx,my = ',nx,ny,mx,my
          if(isymm.eq.2) then
            do i = 0, mx
              do j = 0, my
                bt(i,j) = 0.0_R8
                if(i.eq.0) then
                  bp(i,j)=ubr
                else
                  bp(i,j)=ubr*(crxs(0,j,1)**2-crxs(0,j,3)**2)/
     .             (crxs(i,j,1)**2-crxs(i,j,3)**2)
                  bp(i,j)=min(ubr,bp(i,j))
                endif
                if(j.gt.0) bp(i,j)=max(bp(i,j),bp(i,j-1))
              enddo
            enddo
          else if (isymm.eq.4) then
            do i = 0, mx
              do j = 0, my
                bt(i,j) = 0.0_R8
                if(j.eq.0) then
                  bp(i,j)=ubr
                else
                  bp(i,j)=ubr*(crys(i,0,1)**2-crys(i,0,2)**2)/
     .             (crys(i,j,1)**2-crys(i,j,2)**2)
                  bp(i,j)=min(ubr,bp(i,j))
                endif
                if(i.gt.0) bp(i,j)=max(bp(i-1,j),bp(i,j))
              enddo
            enddo
          endif
          mx=mx-1
          my=my-1
          if(nx.ne.mx) then
            write(*,*) 'b2gxdr_dummy_st: mx.ne.nx:  ',mx,nx
            call xerrab
     .       ('b2gxdr_dummy_st: Sonnet file inconsistent with the grid')
          end if
          if(ny.ne.my) then
            write(*,*) 'b2gxdr_dummy_st: my.ne.ny:  ',my,ny
            call xerrab
     .       ('b2gxdr_dummy_st: Sonnet file inconsistent with the grid')
          end if
          psi=0.0_R8
          ffbz=0.0_R8

!          write(*,*) "Somebody has filled ffbz with zeroes"
!          write(*,*) " ... and is going to read it from some file again"

! IYS 21.02.2020
! I could not find a solution without 'GOTO' here
! Sorry for that ...
          read(lun,'(a)',end=210,err=210) uline
          no_flux_data = .FALSE.
          GOTO 211    ! Wenn es meine Wille gaebe...
  210     no_flux_data = .TRUE.
  211     if (no_flux_data) then
            if (isymm.eq.1.or.isymm.eq.3) then
              ffbz = 2.0_R8 * pi * ubr  ! IYS 14.05.2019
              write(*,*) "No flux data is found the 'old' format file"
              write(*,*) "Trying to recover ffbz using Rtor*Btor"
            endif
          else
            backspace (lun)
            do iy=0,nny-1
              do ix=0,nnx-1
                read(lun,*,err=10,iostat=iostat) ixd,iyd,
     1                  crxs(ix,iy,0),crys(ix,iy,0),
     1                  (psi(ix,iy,i,0),i=1,4),
     $                  (ffbz(ix,iy,i),i=0,3)
              enddo
            enddo
          endif
          close(lun)
! IYS 21.02.2020
        else
          call xerrab('Error in reading geometry file')
        endif
      endif

      call subend ()
      return

   10 continue
      write(*,*) 'Error ',iostat,' in reading ',local_sonnet
      write(*,*) 'ix=',ix,' iy=',iy,' ixd=',ixd,' iyd=',iyd
      call xerrab ('Error in reading grid file')

  710 format(18x,e17.10,14x,e17.10,1x,e17.10)
  720 format(32x,e17.10,1x,e17.10,8x,e17.10,1x,e17.10)
  721 format(30x,e17.10,1x,e17.10,8x,e17.10,1x,e17.10)

  990 continue
      write(*,*) 'Error: i,j ',i,j
      call xerrab ('b2gxdr_dummy_st: Premature end in the grid file')

      end subroutine b2gxdr_dummy_st

