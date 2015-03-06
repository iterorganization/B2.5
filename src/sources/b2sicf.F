*-----------------------------------------------------------------------
*.specification

      subroutine b2sicf (nx, ny, isb, mb, bb, vol, hx, hz, nb, ub,    !srv 09.07.01 28.01.02 22.07.05
     *                  smbcf, smbcb, ctcf, wrk0)                     !srv 06.04.07
      use b2mod_types
      use b2mod_constants
      use b2mod_indirect
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, isb                                             !srv 28.01.02
      real (kind=R8) ::
     *   mb, bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny), hx(-1:nx,-1:ny),  !srv 22.07.05
     *   hz(-1:nx,-1:ny), nb(-1:nx,-1:ny), ub(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   smbcf(-1:nx,-1:ny,0:3), ctcf(-1:nx,-1:ny),                   !srv 06.04.07
     *   smbcb(-1:nx,-1:ny)
*   ..workspace
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2SICF computes centrifugal force
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, style, n2, ix, iy, k
      integer, save :: iout = 0                                     !srv 17.06.02
      character chns*3, chark*1                                     !srv 17.06.02
      real (kind=R8) ::
     *    phm0, phm1, t0
      save ncall, style, phm0, phm1
*   ..procedures
      external subini, subend, xertst, ipgeti, ipgetr
      external b2xvff, b2xvsg
*   ..initialization
      data ncall/0/, style/0/, phm0/0.0e0_R8/, phm1/0.0e0_R8/         !srv 18.01.02

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2sicf')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2sicf_style', style)
       call ipgetr ('b2sicf_phm0', phm0)
       call ipgetr ('b2sicf_phm1', phm1)
       call xertst (0.le.style.and.style.lt.2,
     &   'faulty internal parameter style')
       call ipgeti ('b2sicf_iout', iout)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test mb
      call xertst (0.le.mb, 'faulty argument mb')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of nb
       call b2xvsg (n2, nb, 1, 'nb', '.gt.')
      endif

* ..compute centrifugal force
      smbcf = 0.0e0_R8
	ctcf = 0.0e0_R8                                                 !srv 06.04.07
      if (style.eq.0) then
	 do iy = 0, ny-1
        do ix = 0, nx
	   if(leftix(ix,iy).eq.-2) then
	    wrk0(ix,iy) = 0.0e0_R8
	   else
	    wrk0(ix,iy) = (hz(ix,iy)-hz(leftix(ix,iy),leftiy(ix,iy)))/
     &      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0e0_R8      !srv 05.04.07
	   endif
        enddo
       enddo

*   .. compute smbcf
       do iy = 0, ny-1
        do ix = 0, nx-1                                               !srv 26.11.02 {
	    if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 22.07.05 {
           smbcf(ix,iy,0) = phm0*vol(ix,iy)*               !srv factor hz/hz=1  
     &       (bb(ix,iy,0)/bb(ix,iy,3))*mb*ub(ix,iy)**2                !srv 05.04.07
     &       *nb(ix,iy)*0.5e0_R8*(wrk0(ix,iy)+
     &                            wrk0(rightix(ix,iy),leftiy(ix,iy))) !srv 22.07.05 {
           smbcf(ix,iy,1) = 0.0e0_R8
           smbcf(ix,iy,2) = 0.0e0_R8
           smbcf(ix,iy,3) = 0.0e0_R8
	     ctcf(ix,iy) = phm0*phm1*vol(ix,iy)*                        !srv 06.04.07 {
     &       (bb(ix,iy,0)/bb(ix,iy,3))*mb*2.0e0_R8*ub(ix,iy)
     &       *nb(ix,iy)*0.5e0_R8*(wrk0(ix,iy)+
     &                            wrk0(rightix(ix,iy),leftiy(ix,iy))) !srv 06.04.07 }

	    else                                                        !srv 26.11.02 }
           smbcf(ix,iy,0) = 0.0e0_R8
           smbcf(ix,iy,1) = 0.0e0_R8
           smbcf(ix,iy,2) = 0.0e0_R8
           smbcf(ix,iy,3) = 0.0e0_R8
	     ctcf(ix,iy) = 0.0e0_R8                                     !srv 06.04.07
	    endif                                                       !srv 26.11.02
        enddo
       enddo

      else if (style.eq.1) then

	 do iy = 0, ny-1
        do ix = 0, nx
	   if(leftix(ix,iy).eq.-2) then
	    wrk0(ix,iy) = 0.0e0_R8
	   else
	    wrk0(ix,iy) = hz(ix,iy)-hz(leftix(ix,iy),leftiy(ix,iy))
	   endif
        enddo
       enddo

*   .. compute smbcf
       do iy = 0, ny-1
        do ix = 0, nx-1                                               !srv 26.11.02 {
	 if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then   !srv 22.07.05 {
           smbcf(ix,iy,0) = phm0*vol(ix,iy)*               !srv factor hz/hz=1  
     &       (bb(ix,iy,0)/bb(ix,iy,3))/hx(ix,iy)*mb*ub(ix,iy)**2
     &       *nb(ix,iy)*0.5e0_R8*(wrk0(ix,iy)+
     &                            wrk0(rightix(ix,iy),leftiy(ix,iy))) !srv 22.07.05 {
           smbcf(ix,iy,1) = 0.0e0_R8
           smbcf(ix,iy,2) = 0.0e0_R8
           smbcf(ix,iy,3) = 0.0e0_R8
	     ctcf(ix,iy) = phm0*phm1*vol(ix,iy)*                        !srv 06.04.07 {
     &       (bb(ix,iy,0)/bb(ix,iy,3))/hx(ix,iy)*mb*2.0e0_R8*ub(ix,iy)
     &       *nb(ix,iy)*0.5e0_R8*(wrk0(ix,iy)+
     &                            wrk0(rightix(ix,iy),leftiy(ix,iy))) !srv 06.04.07 }

	 else                                                        !srv 26.11.02 }
           smbcf(ix,iy,0) = 0.0e0_R8
           smbcf(ix,iy,1) = 0.0e0_R8
           smbcf(ix,iy,2) = 0.0e0_R8
           smbcf(ix,iy,3) = 0.0e0_R8
	     ctcf(ix,iy) = 0.0e0_R8                                     !srv 06.04.07
	 endif
         smbcb(ix,iy) =  smbcf(ix,iy,0)
        enddo
       enddo
      endif
	if (iout.ne.0) then                                            !srv 17.06.02 {
	 write (chns,'(i3.3)') isb    
	 do k=0,3
	  write (chark,'(i1)') k  
	  call my_out(70,nx,ny,smbcf(-1,-1,k),'b2sicf_smbcf'//chark//chns)
	 enddo
	 call my_out(70,nx,ny,smbcf(:,:,0)+smbcf(:,:,1)*ub(:,:)+
     &        smbcf(:,:,2)*nb(:,:)*mb+smbcf(:,:,3)*ub(:,:)*nb(:,:)*mb,
     &        'b2sicf_smocf'//chns)
	 call my_out(70,nx,ny,ctcf,'b2sicf_ctcf'//chns)
	endif                                                           !srv 17.06.02 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2sicf

      end subroutine b2sicf
