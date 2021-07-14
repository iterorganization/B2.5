*-----------------------------------------------------------------------
*.specification
      subroutine b2tfvh(nx,ny,ns,hx,bb,vol,na,ua,floi_vhx)
      use b2mod_types
      use b2mod_subsys
      use b2mod_indirect
      use b2mod_zhfrtf
      use b2mod_b2zhco
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx,ny,ns
      real (kind=R8) :: 
     *  bb(-1:nx,-1:ny,0:3),
     &  vol(-1:nx,-1:ny), 
     &  hx(-1:nx,-1:ny),
     &  ua(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1)

*   ..output arguments 
      real (kind=R8) :: 
     *   floi_vhx(-1:nx,-1:ny), cfloi_vh(-1:nx,-1:ny),                  !som 09.07.21
     *   cfloi_vhx(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     b2tfvh computes the ion velocity dependent heat flux 
*
*  3. description 
*
*     This routine computes:
*     -floi_vhx: the ion velocity dependent heat flux 
*
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer ix,iy,inucl
      real (kind=R8) :: temp

      external intfaceh
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2tfvh')
*   ..main
*   ..recalculate velocities and densities
      call b2mod_zhfrtf_prep(nx,ny,ns,na,ua)

*   ..velocity dependent heat flux in the cell center
      do iy = -1, ny
        do ix = -1, nx
        cfloi_vh(ix,iy)=0.0_R8
          do inucl = 0, nspecies-1
            cfloi_vh(ix,iy)=cfloi_vh(ix,iy)+nal(ix,iy,inucl)*
     & sum(c_hwA(ix,iy,0:nspecies-1,inucl)*
     & (av_ualpha(ix,iy,inucl)-av_ualpha(ix,iy,0:nspecies-1)))
          enddo
        cfloi_vhx(ix,iy)=cfloi_vh(ix,iy)*(vol(ix,iy)/hx(ix,iy))*
     & (bb(ix,iy,0)/bb(ix,iy,3))
        enddo
      enddo

*   ..velocity dependent heat flux in the cell boundary
      call intfaceh(nx,ny,leftix,leftiy,vol,cfloi_vhx,floi_vhx)

      call subend ()
      return


      end subroutine b2tfvh


