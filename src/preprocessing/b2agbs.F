      subroutine b2agbs (nx, ny, nncutmax,
     1 leftix, leftiy, rightix, rightiy, region, crx, cry,
     1 ffbz, bb, wbbl, gs, qz, qc, hx, hy, vol,
     1 pbs, cell_width, cell_height, cflags)
      use b2mod_types
      use b2mod_b2cmfs
      use b2mod_constants
      use carre_constants
      use tradui_constants
      use b2mod_connectivity
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nncutmax
      integer leftix(-1:nx,-1:ny), leftiy(-1:nx,-1:ny)
      integer rightix(-1:nx,-1:ny), rightiy(-1:nx,-1:ny)
      integer region(-1:nx,-1:ny,0:2)
      integer cflags(-1:nx,-1:ny,GRID_N_CELLFLAGS)
      real (kind=R8) ::
     *  crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3),
     *  ffbz(-1:nx,-1:ny,0:3),
     *  bb(-1:nx,-1:ny,0:3), wbbl(-1:nx,-1:ny,0:3),
     *  gs(-1:nx,-1:ny,0:2), qz(-1:nx,-1:ny,0:1),
     *  qc(-1:nx,-1:ny,0:1), vol(-1:nx,-1:ny)
*   ..input/output arguments (defined on entry, changed upon exit)
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),
     *  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
*   ..local variables
      integer i, ix, iy, ixl, ixr, iyl, iyr, jxi, jxa, jsep
      integer ntubes(-1:ny), ixref(nncutmax,-1:ny), iseg(-1:nx,-1:ny),
     *  ixleft(nncutmax,-1:ny),
     *  isnum(nncutmax,-1:ny)
      real (kind=R8) :: hz(-1:nx,-1:ny),
     *                  hzl(-1:nx,-1:ny),    ! hz at the left cell face
     *                  pol_flux, width, wrk
      data width /1.0_R8/
      save width
*     ------------------------------------------------------------------
*     B2AGBS computes a correction to the pbsx array to ensure the
*     constancy of magnetic flux along flux tubes.
*     ------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbs')
      call ipgetr ('b2agmt_1d_width', width)
      
*   ..preparatory quantities
      hz(:,:) = vol(:,:)/gs(:,:,2)
      if (isymm.eq.0) then
        hzl(:,:) = width
      else if (isymm.eq.1 .or. isymm.eq.2) then
        hzl(:,:) = pi * ( crx(:,:,0) + crx(:,:,2) )
      else if (isymm.eq.3 .or. isymm.eq.4) then
        hzl(:,:) = pi * ( cry(:,:,0) + cry(:,:,2) )
      end if

*   ..toroidal magnetic field
      if (isymm.eq.1 .or. isymm.eq.3) then
        if (maxval(abs(ffbz(:,:,:))).gt.0.0_R8) then
          bb(:,:,2) =
     .     (ffbz(:,:,0) + ffbz(:,:,1) + ffbz(:,:,2) + ffbz(:,:,3)) /
     .     (4.0_R8 * hz(:,:))
          wbbl(:,:,2) = bb(:,:,2)
          where (hzl(:,:).ne.0.0_R8)
            wbbl(:,:,2) = (ffbz(:,:,0) + ffbz(:,:,2)) /
     .       (2.0_R8 * hzl(:,:))
          end where
        end if
      end if
      wbbl(:,:,1) = 0.0_R8

*   ..identify flux tubes
      iseg = 0
      isnum = 0
      ntubes = 0
      ixleft = -2
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).ne.0) cycle
          if (rightiy(ix,iy).ne.iy) cycle
          if (leftiy(ix,iy).ne.iy) cycle
          if (isGhostCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ntubes(iy) = ntubes(iy) + 1
          call xertst (ntubes(iy).le.nncutmax,
     .     'Too many flux tube segments. Increase DEF_NCUT!')
          iseg(ix,iy) = ntubes(iy)
          isnum(iseg(ix,iy),iy) = 1
          ixl = ix
          ixleft(iseg(ix,iy),iy) = ix
          ixr = ix
   20     continue
          if (leftix(ixl,iy).ne.-2.and.leftiy(ixl,iy).eq.iy) then
            ixl = leftix(ixl,iy)
            iyl = leftiy(ixl,iy)
            if (iseg(ixl,iy).ne.0.and.iyl.eq.iy) then
              call xertst(iseg(ixl,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (left)!')
              goto 30
            else if (iyl.eq.iy) then
              iseg(ixl,iy) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
              if (leftix(ixl,iy).eq.-2) ixleft(iseg(ix,iy),iy) = ixl
            else if (iyl.ne.iy) then
              goto 30
            end if
            goto 20
          end if
   30     continue
          if (rightix(ixr,iy).ne.nx+1.and.rightiy(ixr,iy).eq.iy) then
            ixr = rightix(ixr,iy)
            iyr = rightiy(ixr,iy)
            if (iseg(ixr,iy).ne.0.and.iyr.eq.iy) then
              call xertst(iseg(ixr,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (right)!')
              goto 40
            else if (iyr.eq.iy) then
              iseg(ixr,iyr) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
            else
              goto 40
            end if
            goto 30
          end if
   40     continue
        end do
      end do

*   ..select reference location for each flux tube
*   ..we use the outer midplane whenever possible, then the inner midplane,
*   ..otherwise the cell to the right of the left end of the flux tube.
      ixref = -2
      call get_jsep(nx,ny,jxi,jxa,jsep)
      do iy = -1, ny
        do i = 1, ntubes(iy)
          if (i.eq.iseg(jxa,iy)) then
            ixref(i,iy) = jxa
          else if (i.eq.iseg(jxi,iy)) then
            ixref(i,iy) = jxi
          else
            ixref(i,iy) = rightix(ixleft(i,iy),iy)
          end if
        end do
      end do

*   ..correct pbsx, hx, and hy to enforce:
*   ..B_{pol}*hy1*hz constant (on western cell faces and at cell centers)
*   ..vol = hy1*hx*hz
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).eq.0) cycle
          if (ixref(iseg(ix,iy),iy).eq.-2) cycle

          if (leftix(ixref(iseg(ix,iy),iy),iy).ne.-2) then
            pol_flux = wbbl(ixref(iseg(ix,iy),iy),iy,0) *
     &                   gs(ixref(iseg(ix,iy),iy),iy,0) *
     &                   qc(ixref(iseg(ix,iy),iy),iy,0)
          else
            pol_flux = bb(ixref(iseg(ix,iy),iy),iy,0) *
     &                 hy(ixref(iseg(ix,iy),iy),iy)
          end if

          if (leftix(ix,iy).ne.-2) then ! do not adjust pbs if leftix=-2
                                        ! since there is no western face
            wbbl(ix,iy,0) = pol_flux / (gs(ix,iy,0) * qc(ix,iy,0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
            pbs(ix,iy,0) = pol_flux / wbbl(ix,iy,3)
          else
            wbbl(ix,iy,0) = pol_flux /
     &       (gs(rightix(ix,iy),rightiy(ix,iy),0) *
     &        qc(rightix(ix,iy),rightiy(ix,iy),0))
            wbbl(ix,iy,3) =
     &       sqrt(wbbl(ix,iy,0)**2+wbbl(ix,iy,1)**2+wbbl(ix,iy,2)**2)
          end if
          wrk = pol_flux /
     .               (hz(ix,iy)*qz(ix,iy,1)*bb(ix,iy,0))
          if (wrk .lt. 0.0_R8) then
            write(*,'(a,i2,a,i2,a)') "For cell (ix,iy) (", ix,",", iy, 
     &        ") redefined hy < 0!"
            cycle
          endif
          hy(ix,iy) = wrk
          hx(ix,iy) = vol(ix,iy)/
     .               (hz(ix,iy)*hy(ix,iy)*qz(ix,iy,1))
        end do
      end do

*   ..copying over modified data
      cell_width = hx
      cell_height = hy

      bb(:,:,3)=sqrt(bb(:,:,0)**2+bb(:,:,1)**2+bb(:,:,2)**2)

*   ..return
      call subend ()
      return
      end subroutine b2agbs

