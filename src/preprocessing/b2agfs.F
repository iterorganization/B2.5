*-----------------------------------------------------------------------
*.specification

      subroutine b2agfs (nrid,nx,ny,nx1,ny1,istyle,
     & crx1,cry1,fpsi1,ffbz1,bb1,wbbc1,
     1 cflag,lblgm,leftix1,leftiy1,rightix1,rightiy1,
     2 topix1,topiy1,bottomix1,bottomiy1,
     3 leftcut1,rightcut1,bottomcut1,topcut1,
     4 periodic_bc,nncut,nncutmax)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmfs
      use b2mod_connectivity
      use b2mod_indirect
     . ,only: GRID_N_CELLFLAGS
      use carre_constants
      implicit none
#include "DIMENSIONS.F"
*   ..input arguments (unchanged on exit)
      integer nrid, nx, ny, nx1, ny1, istyle, nncutmax
*   .. output arguments
      real (kind=R8) ::
     * crx1(-1:nx1,-1:ny1,0:3), cry1(-1:nx1,-1:ny1,0:3),
     1 bb1(-1:nx1,-1:ny1,0:3),fpsi1(-1:nx1,-1:ny1,0:3),
     2 ffbz1(-1:nx1,-1:ny1,0:3),wbbc1(-1:nx1+1,-1:ny1+1,0:3)
      integer
     *  leftix1(-1:nx1,-1:ny1),leftiy1(-1:nx1,-1:ny1),
     1  rightix1(-1:nx1,-1:ny1),rightiy1(-1:nx1,-1:ny1),
     2  topix1(-1:nx1,-1:ny1),topiy1(-1:nx1,-1:ny1),
     3  bottomix1(-1:nx1,-1:ny1),bottomiy1(-1:nx1,-1:ny1),
     4  leftcut1(nncutmax),rightcut1(nncutmax),
     5  bottomcut1(nncutmax),topcut1(nncutmax),
     6  periodic_bc,nncut,cflag(-1:nx1,-1:ny1,GRID_N_CELLFLAGS)
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
      integer ix, iy, i, ic, n2, ixl, ixr, iyt, iyb, inseltop, inselbot
      integer Bt_adjust, Bt_reversal
      character*256 local_sonnet
      data local_sonnet /'upgrade.geometry'/
      real (kind=R8) ::
     *  bp((nx1+2)*(ny1+2)),bt((nx1+2)*(ny1+2)),
     .  crxs((nx1+2)*(ny1+2)*5),crys((nx1+2)*(ny1+2)*5),
     .  fpsi((nx1+2)*(ny1+2)*4),ffbz((nx1+2)*(ny1+2)*4)
      real (kind=R8) ::
     *  xoffset, yoffset, xrescale, yrescale, Bt_rescale, pit_rescale,
     .  min_pitch, geom_match_dist, distval,
     .  r1,r2,r3,rf,z1,z2,z3,z4,numer,denom,aaa,bbb,
     .  nhatr,nhatz,bphatr,bphatz,bdotn,beta,cosa,cosb,cose,sinb
      integer ireg,zone(-1:nx1),errors
      logical found, concave
      data xoffset /0.0_R8/, yoffset /0.0_R8/, Bt_rescale/1.0_R8/,
     .  Bt_adjust/0/, Bt_reversal/0/, xrescale/1.0_R8/,
     .  yrescale/1.0_R8/, pit_rescale/1.0_R8/, min_pitch/1.0_R8/,
     .  geom_match_dist/1.0e-6_R8/

*   ..procedures
      external subini, subend, ipgetc, ipgeti, ipgetr, xertst,
     * b2xvsg, b2gxdr_dummy
      intrinsic abs, sqrt, cos, sin, min, max, sign, acos
      real (kind=R8) ::
     *  dist
      integer ix1,iy1,ip1,ix2,iy2,ip2
      dist(ix1,iy1,ip1,ix2,iy2,ip2)=
     1 sqrt((crx1(ix1,iy1,ip1)-crx1(ix2,iy2,ip2))**2+
     2 (cry1(ix1,iy1,ip1)-cry1(ix2,iy2,ip2))**2)
      logical match
      match(ix1,iy1,ix2,iy2)=
     1 (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist

*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2agfs')
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

      call b2gxdr_dummy(local_sonnet, nx1, ny1, istyle, isymm,
     & crxs, crys, bp, bt, fpsi, ffbz)
      n2=(nx1+2)*(ny1+2)
      ic=1
      do  iy = -1, ny1
        do  ix = -1, nx1
          bb1(ix,iy,0)=bp(ic)
          bb1(ix,iy,1)=0
          bb1(ix,iy,2)=bt(ic)
          bb1(ix,iy,3)=sqrt(bp(ic)*bp(ic)+bt(ic)*bt(ic))
          crx1(ix,iy,1)=crxs(ic+2*n2)*xrescale+xoffset   ! bottom right
          cry1(ix,iy,1)=crys(ic+2*n2)*yrescale+yoffset
          crx1(ix,iy,3)=crxs(ic+4*n2)*xrescale+xoffset   ! top right
          cry1(ix,iy,3)=crys(ic+4*n2)*yrescale+yoffset
          crx1(ix,iy,2)=crxs(ic+3*n2)*xrescale+xoffset   ! top left
          cry1(ix,iy,2)=crys(ic+3*n2)*yrescale+yoffset
          crx1(ix,iy,0)=crxs(ic+1*n2)*xrescale+xoffset   ! bottom left
          cry1(ix,iy,0)=crys(ic+1*n2)*yrescale+yoffset
          if (Bt_adjust.eq.1) then
            if (isymm.eq.1 .and. xoffset.ne.0.0_R8) then
              bb1(ix,iy,2)=bt(ic)*crxs(ic)/(crxs(ic)+xoffset)
            else if (isymm.eq.3 .and. yoffset.ne.0.0_R8) then
              bb1(ix,iy,2)=bt(ic)*crys(ic)/(crys(ic)+yoffset)
            end if
            bb1(ix,iy,0)=bp(ic)/bt(ic)*bb1(ix,iy,2)
            bb1(ix,iy,3)=sqrt(bb1(ix,iy,0)**2+bb1(ix,iy,2)**2)
          endif
          fpsi1(ix,iy,0)=fpsi(ic)
          fpsi1(ix,iy,1)=fpsi(ic+1*n2)
          fpsi1(ix,iy,2)=fpsi(ic+2*n2)
          fpsi1(ix,iy,3)=fpsi(ic+3*n2)
          ffbz1(ix,iy,0)=ffbz(ic)
          ffbz1(ix,iy,1)=ffbz(ic+1*n2)
          ffbz1(ix,iy,2)=ffbz(ic+2*n2)
          ffbz1(ix,iy,3)=ffbz(ic+3*n2)
          if (Bt_adjust.eq.1) then
            if (isymm.eq.1 .and. xoffset.ne.0.0_R8) then
              fpsi1(ix,iy,0)=fpsi(ic)/crxs(ic+1*n2)*
     .                               (crxs(ic+1*n2)+xoffset)
              fpsi1(ix,iy,1)=fpsi(ic+1*n2)/crxs(ic+2*n2)*
     .                                    (crxs(ic+2*n2)+xoffset)
              fpsi1(ix,iy,2)=fpsi(ic+2*n2)/crxs(ic+3*n2)*
     .                                    (crxs(ic+3*n2)+xoffset)
              fpsi1(ix,iy,3)=fpsi(ic+3*n2)/crxs(ic+4*n2)*
     .                                    (crxs(ic+4*n2)+xoffset)
            else if (isymm.eq.3 .and. yoffset.ne.0.0_R8) then
              fpsi1(ix,iy,0)=fpsi(ic)/crys(ic+1*n2)*
     .                               (crys(ic+1*n2)+yoffset)
              fpsi1(ix,iy,1)=fpsi(ic+1*n2)/crys(ic+2*n2)*
     .                                    (crys(ic+2*n2)+yoffset)
              fpsi1(ix,iy,2)=fpsi(ic+2*n2)/crys(ic+3*n2)*
     .                                    (crys(ic+3*n2)+yoffset)
              fpsi1(ix,iy,3)=fpsi(ic+3*n2)/crys(ic+4*n2)*
     .                                    (crys(ic+4*n2)+yoffset)
            end if
          end if
          ic=ic+1
        enddo
      enddo
      bb1(:,:,0)=bb1(:,:,0)*Bt_rescale*pit_rescale
      bb1(:,:,1)=bb1(:,:,1)*Bt_rescale*pit_rescale
      fpsi1(:,:,:)=fpsi1(:,:,:)*Bt_rescale*pit_rescale
      if (Bt_reversal.eq.1) then
        bb1(:,:,2)=-bb1(:,:,2)*Bt_rescale
        if (isymm.eq.0) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale
        else if (isymm.eq.1.or.isymm.eq.2) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale*xrescale
        else if (isymm.eq.3.or.isymm.eq.4) then
          ffbz1(:,:,:)=-ffbz1(:,:,:)*Bt_rescale*yrescale
        end if
      else
        bb1(:,:,2)=bb1(:,:,2)*Bt_rescale
        if (isymm.eq.0) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale
        else if (isymm.eq.1.or.isymm.eq.2) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale*xrescale
        else if (isymm.eq.3.or.isymm.eq.4) then
          ffbz1(:,:,:)=ffbz1(:,:,:)*Bt_rescale*yrescale
        end if
      endif
      bb1(:,:,3)=sqrt(bb1(:,:,0)**2+bb1(:,:,1)**2+bb1(:,:,2)**2)
!!! FIXME: For now, until the cflag array is provided by the mesh generator
!!!        as it will be for extended grids
      cflag=0

 9999 continue

      nncut=0
      do ic=1,nncutmax
        bottomcut1(ic)=ny1+1
        topcut1(ic)=-2
        rightcut1(ic)=-2
        leftcut1(ic)=nx1+1
      enddo
      inselbot=ny1+1
      inseltop=-2
      ic=0
      do iy=-1,ny1
        do ix=-1,nx1
          leftix1(ix,iy)=ix-1
          if(ix.gt.-1) then
            if(.not.match(leftix1(ix,iy),iy,ix,iy)) then
              leftix1(ix,iy)=-2
              do ixl=-1,nx1
                if(match(ixl,iy,ix,iy)) then
                  leftix1(ix,iy)=ixl
                  exit
                endif
              enddo
              if(leftix1(ix,iy).eq.-2.and.periodic_bc.eq.0) then
                write(*,*) 'Left neighbour not found for ',ix,iy
              elseif(periodic_bc.eq.1) then
                do ixl=nx1-1,-1,-1
                  if (.not.match(ixl,iy,ixl+1,iy).and.ixl.ne.ix-1) exit
                enddo
                leftix1(ix,iy)=ixl
                if(ic.eq.0) then
                  found=.false.
                else
                  found=.false.
                  do i = 1, ic
                    if(ixl+1.eq.rightcut1(i) .or.
     &                 ixl+1.eq.leftcut1(i)) then
                      found = .true.
                      exit
                    endif
                  end do
                endif
                if (.not.found) then
                  ic = ic + 1
                  write(*,*) 'ic',ic
                  i = ic
                endif
                write(*,'(a,3i4)') 'iy,ix,ixl',iy,ix,ixl
                bottomcut1(i)=min(iy,bottomcut1(i))
                topcut1(i)=max(iy+1,topcut1(i))
                rightcut1(i)=max(ix,rightcut1(i))
              else
                write(*,'(a,3i4)') 'iy,ix,ixl',iy,ix,ixl
                if(ic.eq.0) then
                  found=.false.
                else
                  found=.false.
                  do i = 1, ic
                    if(ixl+1.eq.rightcut1(i) .or.
     &                 ixl+1.eq.leftcut1(i)) then
                      found = .true.
                      exit
                    endif
                  end do
                endif
                if (.not.found) then
                  ic = ic + 1
                  write(*,*) 'ic',ic
                  i = ic
                endif
                bottomcut1(i)=min(iy,bottomcut1(i))
                topcut1(i)=max(iy+1,topcut1(i))
                rightcut1(i)=max(ix,rightcut1(i))
              endif
            endif
          endif
          rightix1(ix,iy)=ix+1
          if(ix.lt.nx1) then
            if(.not.match(ix,iy,rightix1(ix,iy),iy)) then
              rightix1(ix,iy)=nx1+1
              do ixr=-1,nx1
                if(match(ix,iy,ixr,iy)) then
                  rightix1(ix,iy)=ixr
                  exit
                endif
              enddo
              if(rightix1(ix,iy).eq.nx1+1.and.periodic_bc.eq.0) then
                write(*,*) 'Right neighbour not found for ',ix,iy
              elseif(periodic_bc.eq.1) then
                do ixr=0,nx1
                  if(.not.match(ixr-1,iy,ixr,iy).and.ixr.ne.ix+1) exit
                enddo
                rightix1(ix,iy)=ixr
                write(*,'(a,3i4)') 'iy,ix,ixr',iy,ix,ixr
                if(ic.eq.0) then
                  found=.false.
                else
                  found=.false.
                  do i = 1, ic
                    if(ixr.eq.rightcut1(i) .or. ixr.eq.leftcut1(i)) then
                      found = .true.
                      exit
                    endif
                  end do
                endif
                if (.not.found) then
                  ic = ic + 1
                  write(*,*) 'ic',ic
                  i = ic
                endif
                bottomcut1(i)=min(iy,bottomcut1(i))
                topcut1(i)=max(iy+1,topcut1(i))
                leftcut1(i)=min(ixr,leftcut1(i))
              else
                write(*,'(a,3i4)') 'iy,ix,ixr',iy,ix,ixr
                if(ic.eq.0) then
                  found=.false.
                else
                  found=.false.
                  do i = 1, ic
                    if(ixr.eq.rightcut1(i) .or. ixr.eq.leftcut1(i)) then
                      found = .true.
                      exit
                    endif
                  end do
                endif
                if (.not.found) then
                  ic = ic + 1
                  write(*,*) 'ic',ic
                  i = ic
                endif
                bottomcut1(i)=min(iy,bottomcut1(i))
                topcut1(i)=max(iy+1,topcut1(i))
                leftcut1(i)=min(ixr,leftcut1(i))
              endif
            endif
          endif

c Assume rest OK

          topix1(ix,iy)=ix
          bottomix1(ix,iy)=ix
          leftiy1(ix,iy)=iy
          rightiy1(ix,iy)=iy
          topiy1(ix,iy)=iy+1
          bottomiy1(ix,iy)=iy-1

        end do

cxpb Must now re-direct neighbours for periodic boundary BC case when the
cxpb two edges of the mesh coincide in real space.
cxpb The exact position of the "cut" depends on problem implementation.
cxpb The dead region is made into a topological torus so as not to
cxpb interact with the rest of the plasma domain.

        if (periodic_bc.eq.1) then
          ixl=-2
          do i=1,-1,-1
            if(match(nx1-1-i,iy,i,iy)) then
              inselbot=min(inselbot,iy)
              inseltop=max(inseltop,iy)
              leftix1(i,iy)=nx1-1-i
              write(*,'(a,3i4)') 'iy,ix,ixl',iy,i,nx1-1-i
              ixl=i
              rightix1(nx1-1-i,iy)=i
              write(*,'(a,3i4)') 'iy,ix,ixr',iy,nx1-1-i,i
              exit
            end if
          end do
          do i=-1,ixl-1
            rightix1(i,iy)=nx1-1-i
            write(*,'(a,2i4)') 'Dummy cell ix,iy',i,iy
            leftix1(nx1-1-i,iy)=i
            write(*,'(a,2i4)') 'Dummy cell ix,iy',nx1-1-i,iy
          end do
        end if

        nncut=max(nncut,ic)
        if(nncut.gt.nncutmax) then
          write(*,*) ' Increase DEF_NCUT in DIMENSIONS.F!'
          write(*,*) ' nncut = ',nncut,' DEF_NCUT = ',DEF_NCUT
          call xertst (nncut.le.DEF_NCUT,'faulty parameter DEF_NCUT')
        end if
      end do

cxpb  Once we know the vertical boundaries of the periodic region, we can
cxpb  now finish reconnecting the top and bottom rows of the dead regions
cxpb  to each other, and, in the case of the diverted island geometry, also
cxpb  reconnect the guard cells on both sides of the island.

      if (periodic_bc.eq.1) then
        ixl=-2
        do i=1,-1,-1
          if(match(nx1-1-i,(inseltop+inselbot)/2,
     &                   i,(inseltop+inselbot)/2)) then
            ixl=i
            exit
          end if
        end do
        do i=-1,ixl-1
          bottomiy1(i,inselbot)=inseltop
          write(*,*) 'ix,iy,iyb',i,inselbot,inseltop
          bottomiy1(nx1-1-i,inselbot)=inseltop
          write(*,*) 'ix,iy,iyb',nx1-1-i,inselbot,inseltop
          if(inselbot.gt.-1) then                  ! connected region on North side
            topix1(i,inselbot-1)=nx1-1-i           ! connect to the other side
            topiy1(i,inselbot-1)=inselbot-1        ! of the island
            write(*,*) 'ix,iy,ixt,iyt',i,inselbot-1,nx1-1-i,inselbot-1
            topix1(nx1-1-i,inselbot-1)=i           ! connect to the other side
            topiy1(nx1-1-i,inselbot-1)=inselbot-1  ! of the island
            write(*,*) 'ix,iy,ixt,iyt',nx1-1-i,inselbot-1,i,inselbot-1
          end if
          if(inseltop.lt.ny1) then                 ! connected region on South side
            bottomix1(i,inseltop+1)=nx1-1-i        ! connect to the other side
            bottomiy1(i,inseltop+1)=inseltop+1     ! of the limited SOL
            write(*,*) 'ix,iy,ixb,iyb',i,inseltop+1,nx1-1-i,inseltop+1
            bottomix1(nx1-1-i,inseltop+1)=i        ! connect to the other side
            bottomiy1(nx1-1-i,inseltop+1)=inseltop+1 ! of the limited SOL
            write(*,*) 'ix,iy,ixb,iyb',nx1-1-i,inseltop+1,i,inseltop+1
cxpb  Remove neighbours for guard cells next to tangency point
!            bottomiy1(i,inseltop+1)=-2
!            bottomiy1(nx1-1-i,inseltop+1)=-2
cxpb
          endif
          topiy1(i,inseltop)=inselbot
          write(*,*) 'ix,iy,iyt',i,inseltop,inselbot
          topiy1(nx1-1-i,inseltop)=inselbot
          write(*,*) 'ix,iy,iyt',nx1-1-i,inseltop,inselbot
        end do
      end if

cxpb For 1-D cases where we apply periodic boundary conditions to eliminate
cxpb tranverse physics
      if(periodic_bc.eq.-1) then
        if(nx1.eq.1) then
          do iy = -1, ny1
            rightix1(-1,iy) = nx1
            leftix1(nx1,iy) = -1
            rightix1(0,iy) = 0
            leftix1(0,iy) = 0
          enddo
        endif
        if(ny1.eq.1) then
          do ix = -1, nx1
            topiy1(ix,-1) = ny1
            bottomiy1(ix,ny1) = -1
            topiy1(ix,0) = 0
            bottomiy1(ix,0) = 0
          enddo
        endif
      endif

      write(*,*) 'istyle',istyle,'nncut',nncut
      if(nncut.eq.0) write(*,'(a)') 'No cuts found'
      if(nncut.ge.1) write(*,'(a,4i4)')
     1 'Calculated leftcut1, rightcut1, topcut1, bottomcut1 = ',
     1 leftcut1(1), rightcut1(1), topcut1(1), bottomcut1(1)
      if(nncut.ge.2) write(*,'(a,4i4)')
     1 'Calculated leftcut2, rightcut2, topcut2, bottomcut2 = ',
     1 leftcut1(2), rightcut1(2), topcut1(2), bottomcut1(2)

cxpb Check the pitch angle at the plates and make sure it is at least min_pitch
cxpb Do not correct corner cells
      if (min_pitch.gt.0.0_R8) then
        do ix = -1, nx1
          do iy = 0, ny1-1
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
cxpb
cxpb  Additional checks on basis mesh
cxpb
      write(*,*) 'Checking CRX1'
      errors=0
      do iy=-1,ny1
c check right neighbours
        do ix=-1,nx1-1
          ixr=rightix1(ix,iy)
          if(ixr.lt.-1.or.ixr.gt.nx1) then
            write(*,*) 'No right neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. ixr.eq.ix+1) then
            distval=
     1       dist(ix,iy,1,ixr,rightiy1(ix,iy),0)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1         'Error with right neighbour of ',
     2         ix,iy,' (lower point), ','dist=',distval
            endif
            distval=
     1       dist(ix,iy,3,ixr,rightiy1(ix,iy),2)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1         'Error with right neighbour of ',
     2         ix,iy,' (upper point), ','dist=',distval
            endif
          endif
        enddo
c check left neighbours
        do ix=0,nx1
          ixl=leftix1(ix,iy)
          if(ixl.lt.-1.or.ixl.gt.nx1) then
            write(*,*) 'No left neighbour for  ',ix,iy
          elseif(periodic_bc.eq.0 .or. ixl.eq.ix-1) then
            distval=
     1       dist(ix,iy,0,ixl,leftiy1(ix,iy),1)
            if(distval.ne.0) then
            if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1         'Error with left neighbour of ',ix,iy,' (lower point), ',
     2         'dist=',distval
            endif
            distval=
     1       dist(ix,iy,2,ixl,leftiy1(ix,iy),3)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1         'Error with left neighbour of ',ix,iy,' (upper point), ',
     2         'dist=',distval
            endif
          endif
        enddo
      enddo
      do ix=-1,nx1
c check top neighbours
        do iy=-1,ny1-1
          iyt=topiy1(ix,iy)
          if(iyt.lt.-1.or.iyt.gt.ny1) then
            write(*,*) 'No top neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. iyt.eq.iy+1) then
            distval=
     1       dist(ix,iy,2,topix1(ix,iy),iyt,0)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with top neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       dist(ix,iy,3,topix1(ix,iy),iyt,1)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with top neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
c check bottom neighbours
        do iy=0,ny1
          iyb=bottomiy1(ix,iy)
          if(iyb.lt.-1.or.iyb.gt.ny1) then
            write(*,*) 'No bottom neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. iyb.eq.iy-1) then
            distval=
     1       dist(ix,iy,0,bottomix1(ix,iy),iyb,2)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with bottom neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       dist(ix,iy,1,bottomix1(ix,iy),iyb,3)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with bottom neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
      write(*,*) 'Checking CRY1'
      do iy=-1,ny1
c check right neighbours
        do ix=-1,nx1-1
          ixr=rightix1(ix,iy)
          if(ixr.lt.-1.or.ixr.gt.nx1) then
            write(*,*) 'No right neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. ixr.eq.ix+1) then
            distval=
     1       dist(ix,iy,1,ixr,rightiy1(ix,iy),0)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with right neighbour of ',ix,iy,' (lower point), ',
     2       'dist=',distval
            endif
            distval=
     1     dist(ix,iy,3,ixr,rightiy1(ix,iy),2)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with right neighbour of ',ix,iy,' (upper point), ',
     2       'dist=',distval
           endif
          endif
        enddo
c check left neighbours
        do ix=0,nx1
          ixl=leftix1(ix,iy)
          if(ixl.lt.-1.or.ixl.gt.nx1) then
            write(*,*) 'No left neighbour for  ',ix,iy
          elseif(periodic_bc.eq.0 .or. ixl.eq.ix-1) then
            distval=
     1     dist(ix,iy,0,ixl,leftiy1(ix,iy),1)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with left neighbour of ',ix,iy,' (lower point), ',
     2       'dist=',distval
            endif
            distval=
     1       dist(ix,iy,2,ixl,leftiy1(ix,iy),3)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with left neighbour of ',ix,iy,' (upper point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
      do ix=-1,nx1
c check top neighbours
        do iy=-1,ny1-1
          iyt=topiy1(ix,iy)
          if(iyt.lt.-1.or.iyt.gt.ny1) then
            write(*,*) 'No top neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. iyt.eq.iy+1) then
            distval=
     1       dist(ix,iy,2,topix1(ix,iy),iyt,0)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with top neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       dist(ix,iy,3,topix1(ix,iy),iyt,1)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with top neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
c check bottom neighbours
        do iy=0,ny1
          iyb=bottomiy1(ix,iy)
          if(iyb.lt.-1.or.iyb.gt.ny1) then
            write(*,*) 'No bottom neighbour for ',ix,iy
          elseif(periodic_bc.eq.0 .or. iyb.eq.iy-1) then
            distval=
     1       dist(ix,iy,0,bottomix1(ix,iy),iyb,2)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with bottom neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       dist(ix,iy,1,bottomix1(ix,iy),iyb,3)
            if(distval.ne.0) then
              if(distval.gt.geom_match_dist) errors=errors+1
              write(*,'(a,2i4,a,a,1p,1g12.4)')
     1       'Error with bottom neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
      if(errors.eq.0) then
        write(*,*) 'No errors found in basis geometry'
      elseif(errors.eq.1) then
        write(*,*) errors,' error found in the basis geometry'
      else
        write(*,*) errors,' errors found in the basis geometry'
      endif
      call xertst (errors.eq.0, 'Errors found in basis geometry')
cxpb

*   ..test sign of bb3
      call b2xvsg (n2, bb1(-1,-1,3), 1, 'bb3', '.gt.')

      write(*,*) 'X regions (from rightix<>ix+1)'
      do iy=ny1,-1,-1
        ireg=1
        do ix=-1,nx1
          zone(ix)=ireg
          if(rightix1(ix,iy).ne.ix+1) ireg=ireg+1
        enddo
        write(*,'(3276(10i1,1x))') zone
      enddo
      write(*,*) 'X regions (from connections)'
      do iy=ny1,-1,-1
        ireg=0
        zone=0
        do ix=-1,nx1
          if(zone(ix).eq.0) then
            ireg=ireg+1
            zone(ix)=ireg
            ixr=rightix1(ix,iy)
            do while(ixr.ne.ix.and.ixr.ge.-1.and.ixr.le.nx1)
              zone(ixr)=ireg
              call xertst(rightix1(ixr,iy).ne.ixr,
     .         'Problem in connectivity to the right!') 
              ixr=rightix1(ixr,iy)
            enddo
            ixl=leftix1(ix,iy)
            do while(ixl.ne.ix.and.ixl.ge.-1.and.ixl.le.nx1)
              zone(ixl)=ireg
              call xertst(leftix1(ixl,iy).ne.ixl,
     .         'Problem in connectivity to the left!')
              ixl=leftix1(ixl,iy)
            enddo
          endif
        enddo
        write(*,'(3276(10i1,1x))') zone
      enddo

      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2agfs

      end subroutine b2agfs
c
      subroutine b2gxdr_dummy(local_sonnet, nx, ny, istyle, isymm,
     & crxs, crys, bp, bt, fpsi, ffbz)
      use b2mod_types
      use b2mod_constants
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
      integer nx, ny, istyle, isymm
      real (kind=R8) ::
     * crxs(0:nx+1,0:ny+1,0:4),
     1 crys(0:nx+1,0:ny+1,0:4),
     2 bp(0:nx+1,0:ny+1),bt(0:nx+1,0:ny+1),
     3 fpsi(0:nx+1,0:ny+1,0:3),ffbz(0:nx+1,0:ny+1,0:3)
      character*(*) local_sonnet
 
      logical no_flux_data            ! IYS 21.02.2020

      integer nnx, nny,
     1 i, j, ix, iy, ixd, iyd, lun, line_length, mx, my, iostat
#ifdef USE_PXFGETENV
      integer lenval, ierror
#else
#ifdef NAGFOR
      integer lenval, ierror
#endif
#endif
#ifndef NO_GETENV
      character*12 device_env
      logical streql
      external streql
#endif
      real (kind=R8) ::
     *  ubr, ubt, ubp, upsi
      parameter (line_length=256)
      character*256 meshpath
      character*80 format_string(4)
      integer, save :: new_format
      character*12 device
      logical open_file, uop, first_read, ex
      character uline*(line_length)
      external xertst, subini, subend, open_file
      intrinsic index, sqrt
      real (kind=R8) :: t0

      data new_format/1/
      data format_string/
     1 '(19x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     2 '(19x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     3 '(21x,i3,1x,i3,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)',
     4 '(21x,i4,1x,i4,4x,e17.0,1x,e17.0,8x,e17.0,1x,e17.0)'/

*   ..subprogram start-up calls
      call subini('b2gxdr_dummy')
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
#ifdef NAGFOR
      call get_environment_variable('DEVICE',
     . status=ierror, length=lenval)
      if (ierror.eq.0)
     > call get_environment_variable('DEVICE',value=device_env)
#else
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('DEVICE', 0, device_env, lenval, ierror)
#else
      call getenv ('DEVICE', device_env)
#endif
#endif
      if (.not.streql(device_env,' ')) device = device_env
#endif
      if(istyle.eq.-1) then        !xpb  DG-CARRE-format geometry input
        meshpath = trim(meshpath)
     .   //'|modules/Carre/meshes/'//trim(device)
        call find_file(local_sonnet,ex)
        if(open_file(lun,local_sonnet,meshpath)) then
          write(*,*)
     1     "Reading the geometry file assuming the 'new' format"
          read(lun,*) nnx,nny
          if(nnx.ne.nx.or.nny.ne.ny) then
            write(*,'(a,4i5)')
     1       'b2gxdr_dummy: nx,ny code and data disagree',nx,ny,nnx,nny
            call xerrab ('Mismatch nx/ny and nnx/nny')
          endif
          do iy=0,ny+1
            do ix=0,nx+1
              read(lun,*,err=10,iostat=iostat) ixd,iyd,
     1         (crxs(ix,iy,i),crys(ix,iy,i),i=0,4),
     2         bp(ix,iy),bt(ix,iy)
            enddo
          enddo
          fpsi=0.0_R8
          ffbz=0.0_R8
          read(lun,'(a)',end=110) uline
          backspace (lun)
          do iy=0,ny+1
            do ix=0,nx+1
              read(lun,*,err=10,iostat=iostat) ixd,iyd,
     1         (fpsi(ix,iy,i),i=0,3),(ffbz(ix,iy,i),i=0,3)
            end do
          end do
  110     continue
          close(lun)
! If the toroidal flux function is not provided, compute it
          if (maxval(ffbz).eq.0.0_R8 .and. minval(ffbz).eq.0.0_R8) then
            write(*,*) 'Toroidal flux function missing!'
            write(*,*) 'Trying to fix this automatically'
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
            ffbz =2.0_R8 * pi * t0 / (nx + 2) / (ny + 2)
          endif
        else
          call xerrab('Error in reading geometry file')
        endif
      elseif(istyle.eq.-2) then   !xpb  SONNET-format geometry input
        meshpath = trim(meshpath)
     .   //'|modules/DivGeo/device/'//trim(device)
        if(open_file(lun,local_sonnet,meshpath)) then
          write(*,*)
     1     "Reading the geometry file assuming the 'old' format"
          write(*,'(a,i3,2a)') '==== b2gxdr_dummy: read unit',lun,
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
          if(ubr.eq.0.0_R8) then
            write(*,*) 'b2gxdr_dummy: the value of R*Btor undefined.'
            write(*,*)
     1       'It must be defined in the header of the grid file '
            call xerrab ('b2gxdr_dummy: the value of R*Btor undefined')
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
              call xerrab ('b2gxdr_dummy: mx > nx+1')
            endif
          end if
          if(j.gt.my) then
            my=j
            if(my.gt.ny+1) then
              write(*,*) ' my = ',my,' ny = ',ny
              call xerrab ('b2gxdr_dummy: my > ny+1')
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
            write(*,*) 'b2gxdr_dummy: mx.ne.nx:  ',mx,nx
            call xerrab
     .       ('b2gxdr_dummy: Sonnet file inconsistent with the grid')
          end if
          if(ny.ne.my) then
            write(*,*) 'b2gxdr_dummy: my.ne.ny:  ',my,ny
            call xerrab
     .       ('b2gxdr_dummy: Sonnet file inconsistent with the grid')
          end if
          fpsi=0.0_R8
          ffbz=0.0_R8

          write(*,*) "Somebody has filled ffbz with zeros"
          write(*,*) " ... and is going to read it from some file again"

! IYS 21.02.2020
! I slept little last night, so I could not find a solution without 'GOTO' here
! Sorry for that ...
          read(lun,'(a)',end=210) uline
          no_flux_data = .FALSE.
          GOTO 211    ! Wenn es meine Wille gaebe...
  210     no_flux_data = .TRUE.
  211     if (no_flux_data) then
            ffbz = 2.0_R8 * pi * ubr  ! IYS 14.05.2019
            write(*,*) "No flux data is found the 'old' format file"
            write(*,*) "Trying to recover ffbz using RTOR*BTOR"
          else
            backspace (lun)
            do iy=0,ny+1
              do ix=0,nx+1
                read(lun,*,err=10) ixd,iyd,
     1         (fpsi(ix,iy,i),i=0,3),(ffbz(ix,iy,i),i=0,3)
              end do
            end do
          endif
          close(lun)
! IYS 21.02.2020 ... being wanted to sleep



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
      call xerrab ('b2gxdr_dummy: Premature end in the grid file')

      end subroutine b2gxdr_dummy
c
      subroutine init_region(nx,ny,nncut,nncutmax,
     1 leftcut,rightcut,topcut,bottomcut,
     1 leftix,rightix,rightiy,topix,topiy,bottomiy,
     2 region,nnreg,resignore,crx,cry,periodic_bc)
      use b2mod_types
      implicit none
      integer nx,ny,nncut,nncutmax
      integer leftcut(nncutmax),rightcut(nncutmax),
     1 topcut(nncutmax),bottomcut(nncutmax),
     1 leftix(-1:nx,-1:ny),rightix(-1:nx,-1:ny),rightiy(-1:nx,-1:ny),
     1 topix(-1:nx,-1:ny),topiy(-1:nx,-1:ny),bottomiy(-1:nx,-1:ny),
     2 region(-1:nx,-1:ny,0:2),nnreg(0:2),
     3 resignore(-1:nx,-1:ny,1:2),
     4 periodic_bc
      real (kind=R8) ::
     1 crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
c
c            Core case
c    ++-----------------------------++
c    ||                             ||
c    ||                             ||
c    ||              1              ||
c    ||                             ||
c    ||                             ||
c    ++-----------------------------++
c
c    ++-----------------------------++
c    ||                             ||
c    ||                             ||
c    ||1                            ||
c    ||                             ||
c    ||                             ||
c    ++-----------------------------++
c
c    ++--------------2--------------++
c    ||                             ||
c    ||                             ||
c    ||                             ||
c    ||                             ||
c    ||                             ||
c    ++--------------1--------------++
c
c
c Cases with no X-point:
c region 1 === core
c
c                           +++++++++++++++++++++++++
c
c            Limiter case
c    +-------------------------------+
c    |                               |
c    |               2               |
c    |                               |
c    +/-----------------------------\+
c    ||                             ||
c    ||              1              ||
c    ||                             ||
c    ++-----------------------------++
c
c    +-------------------------------+
c    |                               |
c    1                               2
c    |                               |
c    +/-----------------------------\+
c    ||                             ||
c    ||3                            ||
c    ||                             ||
c    ++-----------------------------++
c
c    +---------------3---------------+
c    |                               |
c    |                               |
c    |                               |
c    +/--------------2--------------\+
c    ||                             ||
c    ||                             ||
c    ||                             ||
c    ++--------------1--------------++
c
c
c Cases with no X-point:
c region 1 === core
c region 2 === SOL
c
c                           +++++++++++++++++++++++++
c
c                                       Periodic BC stellarator island case
c    +-------+---------------+-------+  ++-----------------------------++
c    |       :               :       |  ||              5              ||
c    |       :       2       :       |  +\-----------------------------/+
c    |       :               :       |  |               2               |
c    |   3   +---------------+   4   |  +-------+---------------+-------+
c    |       |               |       |  |       |               |       |
c    |       |       1       |       |  |   3   |       1       |   4   |
c    |       |               |       |  |       |               |       |
c    +-------+---------------+-------+  +-------+---------------+-------+
c
c
c    +-------+---------------+-------+  ++-----------------------------++
c    |       :               :       |  ||7                            ||
c    |       2               3       |  +\-----------------------------/+
c    |       :               :       |  |       2               3       |
c    1       +---------------+       4  1-------+---------------+-------4
c    |       |               |       |  |       |               |       |
c    |       |5              |6      |  |       |5              |6      |
c    |       |               |       |  |       |               |       |
c    +-------+---------------+-------+  +-------+---------------+-------+
c
c
c    +---5---+-------6-------+---7---+  ++--------------6--------------++
c    |       :               :       |  ||                             ||
c    |       :               :       |  +\--------------8--------------/+
c    |       :               :       |  |                               |
c    |       +-------4-------+       |  +---5---+-------4-------+---7---+
c    |       |               |       |  |       |               |       |
c    |       |               |       |  |       |               |       |
c    |       |               |       |  |       |               |       |
c    +---1---+-------2-------+---3---+  +---1---+-------2-------+---3---+
c
c Cases with one X-point:
c region 1 === core
c region 2 === SOL
c region 3 === inboard divertor
c region 4 === outboard divertor
c region 5 === island
c
c                           +++++++++++++++++++++++++
c
c                      Connected double-null geometry case
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       :       2       :       ||       :       6       :       |
c    |       :               :       ||       :               :       |
c    |   3   +---------------+   4   ||   7   +---------------+   8   |
c    |       |               |       ||       |               |       |
c    |       |       1       |       ||       |       5       |       |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       2               3       ||       6               7       |
c    |       :               :       ||       :               :       |
c    1       +---------------+       45       +---------------+       8
c    |       |               |       ||       |               |       |
c    |       |9              |10     ||       |11             |12     |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |       +-------4-------+       ||       +------11-------+       |
c    |       |               |       ||       |               |       |
c    |       |               |       ||       |               |       |
c    |       |               |       ||       |               |       |
c    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
c
c Cases with two X-points:
c region 1 === left core
c region 2 === left SOL
c region 3 === left inboard divertor
c region 4 === left outboard divertor
c region 5 === right core
c region 6 === right SOL
c region 7 === right inboard divertor
c region 8 === right outboard divertor
c
c                           +++++++++++++++++++++++++
c
c                      Disconnected double-null geometry case
c                           (Bottom X-point is active)
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       :       2       :       ||       :       6       :       |
c    |       :               +...4...||...7...+               :       |
c    |       :               |       ||       |               :       |
c    |...3...+---------------+       ||       +---------------+...8...|
c    |       |       1       |       ||       |       5       |       |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       2               3       ||       6               7       |
c    |       :               +.......45.......+               :       |
c    |       :               |       ||       |13             :       |
c    1.......+---------------+10     ||       +---------------+.......8
c    |       |9              |       ||       |11             |12     |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |       :               +.......||.......+               :       |
c    |       :               |       ||       |               :       |
c    |.......+-------4-------+       ||       +------11-------+.......|
c    |       |               |       ||       |               |       |
c    |       |               |       ||       |               |       |
c    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
c
c                           - - - - - - - - - - - - -
c
c                      Disconnected double-null geometry case
c                            (Top X-point is active)
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       :       2       :       ||       :       6       :       |
c    |...3...+               :       ||       :               +...8...|
c    |       |               :       ||       :               |       |
c    |       +---------------+...4...||...7...+---------------+       |
c    |       |       1       |       ||       |       5       |       |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       2               3       ||       6               7       |
c    1.......+               :       ||       :               +.......8
c    |       |13             :       ||       :               |       |
c    |       +---------------+.......45.......+---------------+12     |
c    |       |9              |10     ||       |11             |       |
c    |       |               |       ||       |               |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
c
c    +---5---+-------6-------+---7---++--12---+------13-------+--14---+
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |.......+               :       ||       :               +.......|
c    |       |               :       ||       :               |       |
c    |       +-------4-------+.......||.......+------11-------+       |
c    |       |               |       ||       |               |       |
c    |       |               |       ||       |               |       |
c    +---1---+-------2-------+---3---++---8---+-------9-------+--10---+
c
c Cases with two X-points:
c region 1 === left core
c region 2 === left SOL
c region 3 === left inboard divertor
c region 4 === left outboard divertor
c region 5 === right core
c region 6 === right SOL
c region 7 === right inboard divertor
c region 8 === right outboard divertor
c

      integer ix,iy,inseliy,inselix,iyt
      real (kind=R8) ::
     * geom_match_dist
      data geom_match_dist/1.0e-6_R8/
c
      real (kind=R8) ::
     *  dist
      integer ix1,iy1,ip1,ix2,iy2,ip2
      dist(ix1,iy1,ip1,ix2,iy2,ip2)=
     1 sqrt((crx(ix1,iy1,ip1)-crx(ix2,iy2,ip2))**2+
     2 (cry(ix1,iy1,ip1)-cry(ix2,iy2,ip2))**2)
      logical match
      match(ix1,iy1,ix2,iy2)=
     1 (dist(ix1,iy1,1,ix2,iy2,0)+dist(ix1,iy1,3,ix2,iy2,2)).lt.
     2 geom_match_dist
      intrinsic min, max
c
      call ipgetr ('b2agfs_geom_match_dist', geom_match_dist)
      call xertst (geom_match_dist.gt.0.0_R8,
     &  'faulty argument geom_match_dist')
c if stellarator island or limiter case, find periodicity domain coordinate
      if(periodic_bc.eq.1) then
        if (nncut.eq.1) then    ! stellarator island is on North side
          inseliy=ny+1
          inselix=-2
          do iy=ny,-1,-1
            do ix=2,0,-1
              if(match(nx-ix,iy,ix-1,iy)) then
                inselix=ix-1
                inseliy=iy
                exit
              endif
            enddo
          enddo
        else if (nncut.eq.0) then   ! limited region is on South side
          inseliy=-2
          inselix=-2
          do iy=-1,ny
            do ix=2,0,-1
              if(match(nx-ix,iy,ix-1,iy)) then
                inselix=ix-1
                inseliy=iy
                exit
              endif
            enddo
          enddo
        else
          call xerrab ('Unexpected geometry!')
        endif
      endif
c
      region=0
c volume component
      if(nx.eq.1) then  ! 1d radial
        do iy=-1,ny
          do ix=-1,nx
            region(ix,iy,0)=1
          enddo
        enddo
        nnreg(0)=1
      else if(ny.eq.1) then  ! 1d parallel
        do iy=-1,ny
          do ix=-1,nx
            region(ix,iy,0)=1
          enddo
        enddo
        nnreg(0)=1
      else if (nncut.eq.0 .and. periodic_bc.eq.1) then   ! limiter case
        do iy = -1, inseliy
          do ix = inselix, nx-1-inselix
            region(ix,iy,0) = 1
          enddo
        enddo
        do iy = inseliy+1, ny
          do ix = -1, nx
            region(ix,iy,0) = 2
          enddo
        enddo
        if(inseliy.eq.ny) then
           nnreg(0) = 1
        else
           nnreg(0) = 2
        endif
      else if (nncut.eq.0 .and. periodic_bc.eq.0) then   ! 2-D slab case
        do iy = -1, ny
          do ix = -1, nx
            region(ix,iy,0) = 1
          end do
        end do
        nnreg(0) = 1
      else
        if(topcut(1).ge.-1) then
          if(nncut.eq.1) then
            do ix=-1,leftcut(1)-1
              do iy=-1,ny
                if (periodic_bc.ne.1.or.iy.lt.topcut(1))
     1           region(ix,iy,0)=3
                if(periodic_bc.eq.1.and.
     .           iy.ge.topcut(1).and.iy.lt.inseliy)
     .           region(ix,iy,0)=2
                if(periodic_bc.eq.1.and.
     .           iy.ge.inseliy.and.ix.ge.inselix)
     .           region(ix,iy,0)=5
              enddo
            enddo
            do ix=leftcut(1),rightcut(1)-1
              do iy=bottomcut(1),topcut(1)-1
                region(ix,iy,0)=1
              enddo
              if(periodic_bc.eq.0) then
                do iy=topcut(1),ny
                  region(ix,iy,0)=2
                enddo
              elseif(periodic_bc.eq.1) then
                do iy=topcut(1),inseliy-1
                  region(ix,iy,0)=2
                enddo
                do iy=inseliy,ny
                  region(ix,iy,0)=5
                enddo
              endif
            enddo
            do ix=rightcut(1),nx
              do iy=-1,ny
                if (periodic_bc.eq.0.or.iy.lt.topcut(1))
     1           region(ix,iy,0)=4
                if(periodic_bc.eq.1.and.
     .           iy.ge.topcut(1).and.iy.lt.inseliy)
     .           region(ix,iy,0)=2
                if(periodic_bc.eq.1.and.
     .           iy.ge.inseliy.and.ix.le.nx-1-inselix)
     .           region(ix,iy,0)=5
              enddo
            enddo
            if(periodic_bc.eq.1) then
              nnreg(0)=5
            else
              nnreg(0)=4
            endif
          elseif(nncut.eq.2) then
            do ix = -1, leftcut(1)-1
              do iy = -1, ny
                region(ix,iy,0)=3
              enddo
            enddo
            do ix = leftcut(1), leftcut(2)-1
              do iy = max(bottomcut(1),bottomcut(2)), ny
                if (iy.lt.min(topcut(1),topcut(2))) then
                  region(ix,iy,0)=1
                else
                  region(ix,iy,0)=2
                endif
              enddo
            enddo
            do iy = -1, ny
              ix = leftcut(2)
              do while (leftix(ix,iy).ne.-2)
                region(ix,iy,0)=4
                ix=ix+1
              enddo
              ix = rightcut(2)-1
              do while (rightix(ix,iy).ne.nx+1)
                region(ix,iy,0)=7
                ix=ix-1
              enddo
            enddo
            do ix = rightcut(2), rightcut(1)-1
              do iy = max(bottomcut(1),bottomcut(2)), ny
                if (iy.lt.min(topcut(1),topcut(2))) then
                  region(ix,iy,0)=5
                else
                  region(ix,iy,0)=6
                endif
              enddo
            enddo
            do ix = rightcut(1), nx
              do iy = -1, ny
                region(ix,iy,0)=8
              enddo
            enddo
            nnreg(0)=8
          endif
        endif
      endif
      do ix=-1,nx
        do iy=-1,ny
          if(region(ix,iy,0).eq.0) then
            write(*,*) 'REGION not set for ',ix,iy,0
          endif
        enddo
      enddo
c x-flux component
      if(nncut.eq.0 .and. periodic_bc.eq.1) then   ! limiter case
        if(inseliy.eq.ny) then
           do iy = -1, inseliy
              region(inselix+1,iy,1)=1
           enddo
           nnreg(1)=1
           do ix = inselix, nx-1-inselix
              region(topix(ix,-1),topiy(ix,-1),2)=1
              region(ix,ny,2)=2
           enddo
           nnreg(2)=2
        else
           do iy = inseliy+1, ny
              region(rightix(-1,iy),rightiy(-1,iy),1)=1
              region(nx,iy,1)=2
           enddo
           do iy = -1, inseliy
              region(inselix+1,iy,1)=3
           enddo
           nnreg(1)=3
           do ix = inselix, nx-1-inselix
              region(topix(ix,-1),topiy(ix,-1),2)=1
              region(topix(ix,inseliy),topiy(ix,inseliy),2)=2
           enddo
           region(:,ny,2)=3
           nnreg(2)=3
        endif
      elseif(topcut(1).ge.-1) then
        if(periodic_bc.eq.1) then
          iyt=inseliy-1
        else
          iyt=ny
        endif
        do iy=-1,iyt
          region(rightix(-1,iy),rightiy(-1,iy),1)=1
          if(nncut.lt.2) then
            region(nx,iy,1)=4
          elseif(nncut.eq.2) then
            do ix = 0, nx-1
              if(rightix(ix,iy).eq.nx+1) then
                region(ix,iy,1)=4
              elseif(leftix(ix,iy).eq.-2) then
                region(rightix(ix,iy),rightiy(ix,iy),1)=5
              endif
            enddo
            region(nx,iy,1)=8
          endif
        enddo
        do iy=topcut(1),iyt
          region(leftcut(1),iy,1)=2
          if(nncut.eq.1) then
            region(rightcut(1),iy,1)=3
          elseif(nncut.eq.2) then
            region(rightcut(1),iy,1)=7
          endif
        enddo
        do iy=bottomcut(1),topcut(1)-1
          if(nncut.eq.1) then
            region(leftcut(1),iy,1)=5
            region(rightcut(1),iy,1)=6
          elseif(nncut.eq.2) then
            if(iy.ge.bottomcut(2).and.iy.lt.topcut(2)) then
              region(leftcut(1),iy,1)=9
            else
              region(leftcut(1),iy,1)=13
            endif
            region(rightcut(1),iy,1)=12
          endif
        enddo
        if(nncut.eq.2) then
          do iy=bottomcut(2),topcut(2)-1
            region(leftcut(2),iy,1)=10
            if(iy.ge.bottomcut(1).and.iy.lt.topcut(1)) then
              region(rightcut(2),iy,1)=11
            else
              region(rightcut(2),iy,1)=13
            endif
          enddo
          do iy=topcut(2),iyt
            region(leftcut(2),iy,1)=3
            region(rightcut(2),iy,1)=6
          enddo
        endif
        if(nncut.eq.1.and.periodic_bc.eq.1) then
          do iy=iyt+1,ny
            region(inselix+1,iy,1)=7
          enddo
        endif
        if(nncut.eq.1.and.nnreg(0).eq.4) then    ! SN case
          nnreg(1)=6
        elseif(nncut.eq.1.and.nnreg(0).eq.5) then    ! Stellarator island case
          nnreg(1)=7
        elseif(nncut.eq.2.and.topcut(1).eq.topcut(2)) then   ! DNC case
          nnreg(1)=12
        elseif(nncut.eq.2.and.topcut(1).ne.topcut(2)) then   ! DND case
          nnreg(1)=13
        endif
c y-flux component
        if(periodic_bc.eq.1) then
          iyt=topcut(1)
        else
          iyt=ny
        endif
        do ix=-1,leftcut(1)-1
          region(topix(ix,-1),topiy(ix,-1),2)=1
          region(ix,iyt,2)=5
        enddo
        if(nncut.lt.2) then
          do ix=leftcut(1),rightcut(1)-1
            region(topix(ix,-1),topiy(ix,-1),2)=2
            region(ix,topcut(1),2)=4
          enddo
          if(periodic_bc.eq.1) then
            do ix=inselix,nx-1-inselix
              region(ix,ny,2)=6
              region(ix,inseliy,2)=8
            enddo
          else
            do ix=leftcut(1),rightcut(1)-1
              region(ix,ny,2)=6
            enddo
          endif
          do ix=rightcut(1),nx
            region(topix(ix,-1),topiy(ix,-1),2)=3
            region(ix,iyt,2)=7
          enddo
        elseif(nncut.eq.2) then
          do ix=leftcut(1),leftcut(2)-1
            region(topix(ix,-1),topiy(ix,-1),2)=2
            region(ix,min(topcut(1),topcut(2)),2)=4
            region(ix,ny,2)=6
          enddo
          ix = leftcut(2)
          do while (leftix(ix,ny/2).ne.-2)
            region(topix(ix,-1),topiy(ix,-1),2)=3
            region(ix,ny,2)=7
            ix=ix+1
          enddo
          ix = rightcut(2)-1
          do while (rightix(ix,ny/2).ne.nx+1)
            region(topix(ix,-1),topiy(ix,-1),2)=8
            region(ix,ny,2)=12
            ix=ix-1
          enddo
          do ix = rightcut(2),rightcut(1)-1
            region(topix(ix,-1),topiy(ix,-1),2)=9
            region(ix,min(topcut(1),topcut(2)),2)=11
            region(ix,ny,2)=13
          enddo
          do ix = rightcut(1),nx
            region(topix(ix,-1),topiy(ix,-1),2)=10
            region(ix,ny,2)=14
          enddo
        endif
        if(periodic_bc.eq.1) then
          nnreg(2)=8
        elseif(nncut.lt.2) then
          nnreg(2)=7
        elseif(nncut.eq.2) then
          nnreg(2)=14
        endif
      else ! 1-D and 2-D slab cases
        region(0,:,1)=1
        region(nx,:,1)=2
        nnreg(1)=2
        region(:,0,2)=1
        region(:,ny,2)=2
        nnreg(2)=2
      endif
c set up region where the residuals should be ignored (none yet)
      do iy=-1,ny
        do ix=-1,nx
          if(region(ix,iy,0).eq.0 .or.
     &     leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1 .or.
     &     bottomiy(ix,iy).eq.-2 .or. topiy(ix,iy).eq.ny+1) then
            resignore(ix,iy,1)=0
            resignore(ix,iy,2)=0
          else
            resignore(ix,iy,1)=1
            resignore(ix,iy,2)=1
          endif
        enddo
      enddo

      return
      end subroutine init_region
