      subroutine b2nxdp (nx, ny, nregionv, solvereg, qe, te, nep, ne, 
     &  delpo)
      use b2mod_types
      use b2mod_indirect
      use b2mod_subsys
      implicit none
      integer nx, ny, nregionv
      real (kind=R8) ::
     *  qe, te(-1:nx,-1:ny), nep(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     &  delpo(-1:nx,-1:ny)
      logical solvereg(0:nregionv)
*     ------------------------------------------------------------------
*     B2NXDP computes the correction to the electric potential following
*     an update to the electron density.
*     ------------------------------------------------------------------
      integer ncall, n2, ix, iy, ireg
      save ncall
      external xertst
      external b2xvsg
*   ..initialisation
      data ncall/0/
*     ------------------------------------------------------------------
      call subini ('b2nxdp')
      n2 = (nx+2)*(ny+2)
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (0.0_R8.lt.qe, 'faulty argument qe')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, nep, 1, 'nep', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif
*   ..compute delpo
      do iy = -1, ny
       do ix = -1, nx
        ireg = region(ix,iy,0)
        if (solvereg(ireg)) then
         delpo(ix,iy) = (te(ix,iy)/qe)*
     &    (ne(ix,iy)-nep(ix,iy))*2.0_R8/(ne(ix,iy)+nep(ix,iy))
        else
         delpo(ix,iy) = 0.0_R8
        endif
       enddo
      enddo
      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2nxdp
