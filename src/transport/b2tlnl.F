*-----------------------------------------------------------------------
*.specification

      subroutine b2tlnl (nx, ny, icase, te, ti, ne, lnlam)              !srv 01.07.09 {
      use b2mod_types
      use b2mod_constants
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, icase                                             !srv 20.09.11
      real (kind=R8) ::
     &   te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny)              !srv 20.09.11
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   lnlam(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TLNL computes the Coulomb logarithm according to Braginskii or 
*     Wesson formulas.
*
*  3. description (see also routine b2cdca)
*
*  4. references
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, n2, ix, iy
      real (kind=R8), save :: lambda=-5.0e0_R8
      save ncall
*   ..procedures
      external subini, subend, xertst, sfill
      external b2xvsg
*   ..initialization
      data ncall/0/

*-----------------------------------------------------------------------
*.documentation-internal
*
*     Description of some local variables:
*
*     lnlam is the log of the plasma parameter.
*
*     Description of some local procedures:
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tlnl')
      if (ncall.eq.0) then
       call ipgetr('b2trcl_lambda', lambda)
       call xertst(0.0_R8.ne.lambda, 'faulty internal parameter lambda') !srv 02.07.13
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

* ..compute lnlam
      if (lambda.lt.0.0_R8) then
        select case(icase)
         case(0)      ! Braginsky formula
        do iy = -1, ny
         do ix = -1, nx
            if(te(ix,iy)/ev .le. 50.0_R8) then
           lnlam(ix,iy)=max(-lambda,
     &     23.4_R8 - 1.15_R8*log10(ne(ix,iy)/1.0e6_R8) +
     &                            3.45_R8*log10(te(ix,iy)/ev))
            else
           lnlam(ix,iy)=max(-lambda,
     &             25.3_R8 - 1.15_R8*log10(ne(ix,iy)/1.0e6_R8) + 
     &                            2.30_R8*log10(te(ix,iy)/ev))
            endif
            if(lnlam(ix,iy).le.0.0_R8) then
             open(70,file='info_b2tlnl.dat')
             write(70,"(a)")
     &     ' ix   iy   ne(ix,iy)      te(ix,iy)  lnlam(ix,iy)'
             write(70,"(2i5,3es13.4)")
     &       ix,iy,ne(ix,iy),te(ix,iy),lnlam(ix,iy)
             close(70)
             write(*,"(a)")
     &     ' ix   iy   ne(ix,iy)      te(ix,iy)  lnlam(ix,iy)'
             write(*,"(2i5,3es13.4)")
     &       ix,iy,ne(ix,iy),te(ix,iy),lnlam(ix,iy)
             call xertst (lnlam(ix,iy).gt.0.0_R8,
     &                  'Negative Coulomb logarithm')     !srv 12.01.16
            endif
         enddo
        enddo
         case(1)     ! Wesson formula for ion-ion collisions
        do iy = -1, ny
         do ix = -1, nx
           lnlam(ix,iy)=max(-lambda,17.3_R8-
     &       log(ne(ix,iy)/1.0e20_R8)/2.0_R8+
     &       1.5_R8*log(ti(ix,iy)/ev/1000.0_R8))

         enddo
        enddo
         case(2)    ! Wesson formula for electron-electron collisions
        do iy = -1, ny
         do ix = -1, nx
           lnlam(ix,iy)=max(-lambda,14.9_R8-
     &     log(ne(ix,iy)/1.0e20_R8)/2.0_R8+
     &     log(te(ix,iy)/ev/1000.0_R8))
         enddo
        enddo
         case(3)   !  Wesson formula for electron-ion collisions
        do iy = -1, ny
         do ix = -1, nx
           lnlam(ix,iy)=max(-lambda,15.2-
     &    log(ne(ix,iy)/1e20_R8)/2.0_R8+
     &    log(te(ix,iy)/ev/1000.0_R8))
         enddo
        enddo
         case default
          write(*,*) ' b2tlnl: incorrect value icase =',icase
          call xerrab ('Incorrect icase value in b2tlnl')
         end select
      else
        lnlam = lambda
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tlnl

      end                                                            !srv 01.07.09 }
