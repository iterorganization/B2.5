      subroutine b2nxdv (nx, ny, isb, rob, smqb, smqdu, smb)
      use b2mod_types
      use b2mod_indirect
      use b2mod_sources !djm Jan2017
     & , only : b2nxdv_smo0to3
      implicit none
      integer nx, ny, isb
      real (kind=R8) ::
     *   rob(-1:nx,-1:ny), smqb(-1:nx,-1:ny,0:3), 
     &  smqdu(-1:nx,-1:ny), smb(-1:nx,-1:ny,0:3)
      real (kind=R8) :: smbtmp(-1:nx,-1:ny,0:3) !djm Jan2017
*     ------------------------------------------------------------------
*     B2NXDV computes the correction to the parallel momentum source
*     that is required to make the total friction force cancel.
*     ------------------------------------------------------------------
      integer ncall, ix, iy
      integer, save :: style = 1                                       !srv 26.10.08
      save ncall
      external subini, subend, xertst
*   ..initialisation
      data ncall/0/
*     ------------------------------------------------------------------
      call subini ('b2nxdv')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2nxdv_style', style)                            !srv 26.10.08
      endif
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')

cdjm Jan2017 Keep source at the start to calculate total change for balance routines
      smbtmp = smb

*   ..adjust smb
      if(style.eq.0) then
       do iy = -1, ny
        do ix = -1, nx
         smb(ix,iy,0) = smb(ix,iy,0)-smqb(ix,iy,1)*smqdu(ix,iy)-
     &    smqb(ix,iy,3)*rob(ix,iy)*smqdu(ix,iy)
        enddo
       enddo
      else
       do iy = 0, ny-1
        do ix = 0, nx-1
	 if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then     !sv 26.11.02  !sv 01.07.05
          smb(ix,iy,0) = smb(ix,iy,0)-smqb(ix,iy,1)*smqdu(ix,iy)-
     &                  smqb(ix,iy,3)*rob(ix,iy)*smqdu(ix,iy)
	 endif                                                         !sv 26.11.02
        enddo
       enddo
      endif

cdjm Store the total change in linearised source for balance
      b2nxdv_smo0to3(-1:nx,-1:ny,0:3,isb) = smb-smbtmp

      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
