*-----------------------------------------------------------------------
*.specification

      subroutine b2tstch (nx, ny, vol, hy, qz, bb, po, te, ne,          !srv 07.07.14
     *           csig_stoch, fchstoch)                                  !lk 31.01.08 {
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), po(-1:nx,-1:ny), 
     1 te(-1:nx,-1:ny), ne(-1:nx,-1:ny),
     1 hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &    csig_stoch(-1:nx,-1:ny), fchstoch(-1:nx,-1:ny)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2Tstch computes the y-component of electron current 
*       associated with stochastization of magnetic field
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, ix, iy, bottomdelta1, ixout
      real (kind=R8) :: sig_stoch, delta_stoch, width
      save ncall, sig_stoch, delta_stoch,
     &  bottomdelta1, ixout
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   hy1                                                            !srv 07.07.14
      external subini, subend, xertst, xerrab
      external b2xvsg, b2xvff, b2xvfx
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)                                !srv 04.07.14
*   ..initialization
      data ncall/0/, delta_stoch/0.0e0_R8/, sig_stoch/0.0e0_R8/ 
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tstch')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tstch_delta', delta_stoch)
       call ipgetr ('b2tstch_sig', sig_stoch)
       call xertst (0.0_R8.lt.delta_stoch,
     &	 'faulty internal parameter delta_stoch')
       call xertst (0.0_R8.le.sig_stoch,
     &	 'faulty internal parameter sig_stoch')
        
       if(nncut.eq.2) then
         ixout = rightcut(nncut)
         do ix = rightcut(nncut)+1, rightcut(nncut-1)-1
           if(abs(bb(ix,ny,2)) .lt. abs(bb(ixout,ny,2))) ixout = ix  !srv 11.01.13
         enddo
       elseif(nncut.eq.1) then
         ixout = leftcut(nncut)
         do ix = leftcut(nncut)+1, rightcut(nncut)-1
           if(abs(bb(ix,ny,2)) .lt. abs(bb(ixout,ny,2))) ixout = ix  !srv 11.01.13
         enddo
       else
         write(*,*) ' b2tqna: nncut not equal 1 or 2'
         call xerrab ('Not available outside SN and DN geometries')
       endif

       bottomdelta1 = topcut(1)
       width = 0.0e0_R8
       do iy =topcut(1)-1,0,-1
         width = width + hy(ixout,iy)
         if(width.gt.delta_stoch) then
           bottomdelta1 = iy+1
           exit
         endif
       enddo
       bottomdelta1 = bottomdelta1 + 1         !srv 17.12.13   plus 1 specially
!       write(*,"(a,3i5)") ' b2tstch: ixout, bottomdelta1, topcut(1)',
!     &                               ixout, bottomdelta1, topcut(1)
      endif

      csig_stoch = 0.0e0_R8
      fchstoch   = 0.0e0_R8

      do iy = bottomdelta1, topcut(1)
        do ix = -1, nx
         if(mod(region(ix,bottomiy(ix,iy),0),4).eq.1) then        
          csig_stoch(ix,iy) = sig_stoch*
     &    qe*(ne(ix,iy) + ne(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8*  !lk 24.03.08 17.12.13
     &    (vol(ix,iy) + vol(bottomix(ix,iy),bottomiy(ix,iy)))*2.0e0_R8
     &    /(hy1(ix,iy) + hy1(bottomix(ix,iy),bottomiy(ix,iy)))**2        !srv 07.07.14
          fchstoch(ix,iy) = -csig_stoch(ix,iy)*
     &     ((po(ix,iy) - po(bottomix(ix,iy),bottomiy(ix,iy)))-
     &     0.5e0_R8/qe*(te(ix,iy) - te(bottomix(ix,iy),bottomiy(ix,iy)))
     &    -(te(ix,iy) + te(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8/qe
     &    *(log(ne(ix,iy)) - log(ne(bottomix(ix,iy),bottomiy(ix,iy)))))
         endif
        enddo
      enddo

!      call my_out(70,nx,ny,csig_stoch,'b2tstch_csigy_stoch')
!      call my_out(70,nx,ny,fchstoch,'b2tstch_fchstochy')

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tstch

      end                                                              !lk 31.01.08 }
