*-----------------------------------------------------------------------
*.specification

      subroutine b2xehy (nx, ny, hy, qz, te, po, ne, ehy)       !srv 13.01.17
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     &   hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), te(-1:nx,-1:ny),
     &   po(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &   ehy(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2XEHY computes the modified electrical field on y - cell faces.
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, n2, ix, iy
      integer, save :: mode_ehy = 0                                     !srv 06.04.15
      real (kind=R8) ::
     *   t0, hy1_(-1:nx,-1:ny)
      save ncall
*   ..procedures
      external subini, subend, xertst
      external b2xvsg
      real (kind=R8) ::
     *   hy1
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 07.05.08 04.07.14
*   ..initialization
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2xehy')
*   ..set internal parameters on first call                             !srv
      if (ncall.eq.0) then                                              !srv
       call ipgeti ('b2tfhe_mode_ehy',mode_ehy)                         !srv 06.04.15
      endif                                                             !srv
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, te, csig
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
      endif

      do iy=-1,ny
        do ix=-1,nx
          hy1_(ix,iy)=hy1(ix,iy)
        enddo
      enddo

*   ..compute ehy
      do iy = -1, ny
        do ix = -1, nx
          if(bottomiy(ix,iy).eq.-2) then
            ehy(ix,iy) = 0.0_R8
	  else
	    if(mode_ehy.eq.0) then                                      !srv 06.04.15 {
	      ehy(ix,iy) =
     &         -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &          (ne(ix,iy)*te(ix,iy)-
     &           ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &	         te(bottomix(ix,iy),bottomiy(ix,iy)))*
     &	        (2.0_R8/qe)/
     &          (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))
	    else if(mode_ehy.eq.1) then
              t0= (ne(bottomix(ix,iy),bottomiy(ix,iy))*hy1_(ix,iy)+
     &             ne(ix,iy)*hy1_(bottomix(ix,iy),bottomiy(ix,iy)))/
     &            (hy1_(bottomix(ix,iy),bottomiy(ix,iy))+hy1_(ix,iy))
	      ehy(ix,iy) =
     &         -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &  	(ne(ix,iy)*te(ix,iy)-
     &           ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &	         te(bottomix(ix,iy),bottomiy(ix,iy)))/(qe*t0)
	    else
	      ehy(ix,iy) =
     &         -(po(ix,iy)-po(bottomix(ix,iy),bottomiy(ix,iy)))+
     &          (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy))+
     &          (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))*
     &           0.5e0_R8*
     &          (log(ne(ix,iy))-
     &           log(ne(bottomix(ix,iy),bottomiy(ix,iy)))))/qe          !srv 06.04.15 }
            endif
          endif
	enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2xehy

      end
