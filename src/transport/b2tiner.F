*-----------------------------------------------------------------------
*.specification

      subroutine b2tiner (nx, ny, ns, vol, gs, hx, hy, bb, na, ua,    !srv 20.09.06 {
     &                    facdrift, fchinert)                         !srv 22.06.08
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmfs                                                !srv 22.06.08
      use b2mod_b2cmpa                                                !srv 22.06.08
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     & vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2),
     & na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1 hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     1 bb(-1:nx,-1:ny,0:3), facdrift(-1:nx,-1:ny)                     !srv 22.06.08
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &     fchinert(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TINER computes the X,Y-COMPONENTS OF INERTIAL CURRENT
*
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..global named constants
      
*   ..local variables
      integer ncall, n2, ix, iy, ib, is
      integer, save  :: fch_inert_style = 0                              !srv 19.04.13
      real (kind=R8) ::
     &  fch_inert
      save ncall, fch_inert
*   ..procedures
      intrinsic max, abs
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*   ..initialization
      data ncall/0/, fch_inert/1.0e0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tiner')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2tiner_inert', fch_inert)                           !srv 07.07.99
       call ipgeti ('b2tiner_fch_inert_style', fch_inert_style)           !srv 07.07.99
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
      endif
c Drift terms
	fchinert =0.0e0_R8
	 if(fch_inert.ne.0.0e0_R8)then                                    !srv 07.07.99
        if(fch_inert_style.eq.0) then	
	   do is = 0, ns-1                                                 !srv 07.07.99
	    if (.not.is_neutral(is)) then                                  !srv 07.07.99
*         X-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny-1								       !srv 07.07.99
             do ix=0,nx						                   !srv 07.07.99
	        if(leftix(ix,iy).ne.-2) then                               !srv 26.11.02
               fchinert(ix,iy,0)=fchinert(ix,iy,0)					   !srv 07.07.99
     &         +0.5e0_R8*am(is)*mp*0.5e0_R8*fch_inert*     ! signmf*     !srv 23.07.05 11.01.13
     &          facdrift(ix,iy)*                                         !srv 22.06.08
     &          (ua(ix,iy,is)**2+ua(leftix(ix,iy),leftiy(ix,iy),is)**2)  !srv 07.07.99 05.02.05 23.07.05
     &         *(na(ix,iy,is)+na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8 !srv 07.07.99
     &         *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8   !srv 07.07.99
     &         *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))            !srv 07.07.99
     &         /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))          !srv 07.07.99
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 07.07.99
	        endif                                                      !srv 26.11.02
             enddo								       !srv 07.07.99
           enddo									       !srv 07.07.99 
*         Y-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny								       !srv 07.07.99
             do ix=0,nx-1						             !srv 07.07.99
	        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then  !srv 26.11.02
               fchinert(ix,iy,1)=fchinert(ix,iy,1)			       !srv 07.07.99
     &          -0.5e0_R8*am(is)*mp*fch_inert*facdrift(ix,iy) ! signmf*  !srv 07.07.99 21.04.00 22.06.08 11.01.13
     &          *(                                                       !srv 07.07.99
     &           na(ix,iy,is)+na(bottomix(ix,iy),bottomiy(ix,iy),is)     !srv 07.07.99
     &           )/2.0e0_R8                                              !srv 07.07.99
     &          *0.5e0_R8*(ua(ix,iy,is)**2 +                             !srv 07.07.99 18.02.00 23.07.05
     &                     ua(bottomix(ix,iy),bottomiy(ix,iy),is)**2)    !srv 07.07.99 23.07.05
     &          *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          /2.0e0_R8                                                !srv 07.07.99
     &          *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 07.07.99
     &          /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 07.07.99
	        endif                                                      !srv 26.11.02
             enddo								       !srv 07.07.99
           enddo								             !srv 07.07.99
          endif                                                          !srv 07.07.99
         enddo                                                           !srv 07.07.99
        elseif(fch_inert_style.eq.1) then
	   do is = 0, ns-1                                                 !srv 07.07.99
	    if (.not.is_neutral(is)) then                                  !srv 07.07.99
*         X-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny-1								       !srv 07.07.99
             do ix=0,nx						                   !srv 07.07.99
	        if(leftix(ix,iy).ne.-2) then                               !srv 26.11.02
               fchinert(ix,iy,0)=fchinert(ix,iy,0)					   !srv 07.07.99
     &         +0.5e0_R8*am(is)*mp*0.5e0_R8*fch_inert*     ! signmf*     !srv 23.07.05 11.01.13
     &          facdrift(ix,iy)*                                         !srv 22.06.08
     &          (ua(ix,iy,is)**2*na(ix,iy,is)+                         ! IYS 01.02.2013
     &           ua(leftix(ix,iy),leftiy(ix,iy),is)**2*                ! IYS 01.02.2013
     &          na(leftix(ix,iy),leftiy(ix,iy),is))/2.0_R8
     &         *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0_R8   !srv 07.07.99
     &         *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))            !srv 07.07.99
     &         /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))          !srv 07.07.99
     &         *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))    !srv 07.07.99
	        endif                                                      !srv 26.11.02
             enddo								       !srv 07.07.99
           enddo									       !srv 07.07.99 
*         Y-COMPONENT OF INERTIA CURRENT                                 !srv 07.07.99
           do iy=0,ny								       !srv 07.07.99
             do ix=0,nx-1						             !srv 07.07.99
	        if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then  !srv 26.11.02
               fchinert(ix,iy,1)=fchinert(ix,iy,1)			       !srv 07.07.99
     &          -0.5e0_R8*am(is)*mp*fch_inert*facdrift(ix,iy) ! signmf*  !srv 07.07.99 21.04.00 22.06.08 11.01.13
     &          *0.5e0_R8*(na(ix,iy,is)*ua(ix,iy,is)**2 +                ! IYS 01.02.2013
     &          na(bottomix(ix,iy),bottomiy(ix,iy),is)*                  ! IYS 01.02.2013
     &          ua(bottomix(ix,iy),bottomiy(ix,iy),is)**2)               ! IYS 01.02.2013                                                 !sv 07.07.99
     &          *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          /2.0e0_R8                                                !srv 07.07.99
     &          *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))       !srv 07.07.99
     &          /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))     !srv 07.07.99
     &          *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-
     &            OnedBsqc(ix,iy))                                       !srv 07.07.99
	        endif                                                      !srv 26.11.02
             enddo								       !srv 07.07.99
           enddo								             !srv 07.07.99
          endif                                                          !srv 07.07.99
         enddo                                                           !srv 07.07.99
        else
          CALL XERRAB ('B2TINER -- INCORRECT VALUE OF fch_inert_style')
        endif
	 endif                                                            !srv 07.07.99
*  INTERIOR BOUNDARY AT DEAD REGIONS
       DO IB=1,NNISO
	  DO IX=NXISO1(IB), NXISO2(IB)+1
	   DO IY=NYISO1(IB), NYISO2(IB)
	    fchinert(ix,iy,0)=0.0e0_R8
	   end do
	  end do
	  DO IX=NXISO1(IB), NXISO2(IB)
	   DO IY=NYISO1(IB), NYISO2(IB)+1
	    fchinert(ix,iy,1)=0.0e0_R8
	  end do
	  end do
       end do
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tiner

      end                                                              !srv 20.09.06 }
