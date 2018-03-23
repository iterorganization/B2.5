      subroutine b2agbs (nx, ny, nncutmax, leftix, rightix, region,
     1 bb, wbbl, gs, qz, qc, hx, hy, vol, pbs, cell_width, cell_height)
      use b2mod_types
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, nncutmax
      integer leftix(-1:nx,-1:ny), rightix(-1:nx,-1:ny)
      integer region(-1:nx,-1:ny,0:2)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), wbbl(-1:nx,-1:ny,0:3),
     *  gs(-1:nx,-1:ny,0:2), qz(-1:nx,-1:ny,0:1),
     *  qc(-1:nx,-1:ny,0:1), vol(-1:nx,-1:ny)
*   ..input/output arguments (defined on entry, changed upon exit)
      real (kind=R8) ::
     *  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),
     *  cell_width(-1:nx,-1:ny), cell_height(-1:nx,-1:ny)
*   ..local variables
      integer i, j, ix, iy, ixl, ixr, jxi, jxa, jsep
      integer ntubes(-1:ny), ixref(nncutmax,-1:ny),
     *  iseg(-1:nx,-1:ny), nseg(-1:nx,-1:ny),
     *  ixleft(nncutmax,-1:ny), ixright(nncutmax,-1:ny),
     *  isnum(nncutmax,-1:ny)
      real (kind=R8) :: hz(-1:nx,-1:ny)
*     ------------------------------------------------------------------
*     B2AGBS computes a correction to the pbsx array to ensure the
*     constancy of magnetic flux along flux tubes.
*     ------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2agbs')
      
*   ..preparatory quantities
      hz(:,:) = vol(:,:)/gs(:,:,2)

*   ..identify flux tubes
      iseg = 0
      isnum = 0
      ixleft = -2
      ixright = -2
      do iy = -1, ny
        ntubes(iy) = 0
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).ne.0) cycle
          ntubes(iy) = ntubes(iy) + 1
          call xertst (ntubes(iy).le.nncutmax,
     .     'Too many flux tube segments. Increase DEF_NCUT!')
          iseg(ix,iy) = ntubes(iy)
          isnum(iseg(ix,iy),iy) = 1
          ixl = ix
          ixleft(iseg(ix,iy),iy) = ix
          ixr = ix
          ixright(iseg(ix,iy),iy) = ix
 20       continue
          if (leftix(ixl,iy).ne.-2) then
            ixl = leftix(ixl,iy)
            if (iseg(ixl,iy).ne.0) then
              call xertst(iseg(ixl,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (left)!')
              goto 30
            else
              iseg(ixl,iy) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
              if (leftix(ixl,iy).eq.-2) ixleft(iseg(ix,iy),iy) = ixl
            end if
            goto 20
          end if
 30       continue
          if (rightix(ixr,iy).ne.nx+1) then
            ixr = rightix(ixr,iy)
            if (iseg(ixr,iy).ne.0) then
              call xertst(iseg(ixr,iy).eq.iseg(ix,iy),
     .         'Mismatched segments (right)!')
              goto 40
            else
              iseg(ixr,iy) = iseg(ix,iy)
              isnum(iseg(ix,iy),iy) = isnum(iseg(ix,iy),iy) + 1
              if (rightix(ixr,iy).eq.nx+1) ixright(iseg(ix,iy),iy) = ixr
            end if
            goto 30
          end if
 40       continue
        end do
      end do

*   ..select reference location for each flux tube
*   ..we use the outer midplane whenever possible, then the inner midplane,
*   ..otherwise the cell half-way along the flux tube in computational space
      ixref = -2
      call get_jsep(nx,ny,jxi,jxa,jsep)
      do iy = -1, ny
        do i = 1, ntubes(iy)
          if (i.eq.iseg(jxa,iy)) then
            ixref(i,iy) = jxa
          else if (i.eq.iseg(jxi,iy)) then
            ixref(i,iy) = jxi
          else
            ix = ixleft(i,iy)
            do j = 1, isnum(i,iy)/2
              ix = rightix(ix,iy)
            end do
            ixref(i,iy) = ix
          end if
        end do
      end do

*   ..correct pbsx, hx, and hy to enforce:
*   ..B_{pol}*hy1*hz constant (at western cell face and cell centers)
*   ..vol = hy1*hx*hz
      do iy = -1, ny
        do ix = -1, nx
          if (region(ix,iy,0).eq.0) cycle
          if (iseg(ix,iy).eq.0) cycle
          if (ixref(iseg(ix,iy),iy).eq.-2) cycle
          if (leftix(ix,iy).ne.-2) then ! do not adjust pbs if leftix=-2
                                        ! since there is no western face
            pbs(ix,iy,0) =
     .         gs(ixref(iseg(ix,iy),iy),iy,0)*
     .       wbbl(ixref(iseg(ix,iy),iy),iy,0)/wbbl(ix,iy,3)*
     .         qc(ixref(iseg(ix,iy),iy),iy,0)
          end if
          hy(ix,iy) = hz(ixref(iseg(ix,iy),iy),iy)*
     .                hy(ixref(iseg(ix,iy),iy),iy)*
     .                qz(ixref(iseg(ix,iy),iy),iy,1)*
     .                bb(ixref(iseg(ix,iy),iy),iy,0)/
     .               (hz(ix,iy)*qz(ix,iy,1)*bb(ix,iy,0))
          hx(ix,iy) = vol(ix,iy)/
     .               (hz(ix,iy)*hy(ix,iy)*qz(ix,iy,1))
        end do
      end do

*   ..copying over modified data
      cell_width = hx
      cell_height = hy

*   ..return
      call subend ()
      return
      end subroutine b2agbs

