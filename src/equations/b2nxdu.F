*-----------------------------------------------------------------------
*.specification

      subroutine b2nxdu (nx, ny, ns, smq, smqdu, wrk0, wrk1)
      use b2mod_types
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     &  smq(-1:nx,-1:ny,0:3,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smqdu(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2NXDU computes the artificial velocity shift that is
*     required to make the total numerical friction force cancel.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is
      save ncall
*   ..procedures
      external xertst, sfill
      external b2xvsg, b2xvfx
*   ..initialisation
      data ncall/0/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2nxdu')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then

*    ..test sign of smq
       do is = 0, ns-1
        call b2xvsg (n2, smq(-1,-1,1,is), 1, 'smq1', '.le.')
        call b2xvsg (n2, smq(-1,-1,3,is), 1, 'smq3', '.le.')
       enddo
      endif

* ..compute smqdu
*   ..compute auxiliaries wrk0 and wrk1
      call sfill (n2, 0.0_R8, wrk0, 1)
      call sfill (n2, 0.0_R8, wrk1, 1)
      do is = 0, ns-1
       do iy = -1, ny
        do ix = -1, nx
         wrk0(ix,iy) = wrk0(ix,iy)+smq(ix,iy,1,is)+
     &     smq(ix,iy,3,is)*roxa(ix,iy,is)
         wrk1(ix,iy) = wrk1(ix,iy)+smq(ix,iy,0,is)+
     &     smq(ix,iy,1,is)*ua(ix,iy,is)+
     &     smq(ix,iy,2,is)*roxa(ix,iy,is)+
     &     smq(ix,iy,3,is)*roxa(ix,iy,is)*ua(ix,iy,is)
        enddo
       enddo
      enddo
*   ..compute smqdu
      do iy = -1, ny
       do ix = -1, nx
        if (wrk0(ix,iy).eq.0.0_R8) then
         smqdu(ix,iy) = 0.0_R8
        else
         smqdu(ix,iy) = wrk1(ix,iy)/wrk0(ix,iy)
        endif
       enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2nxdu

      end subroutine b2nxdu
