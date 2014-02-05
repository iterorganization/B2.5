      subroutine write_results(nx, ny)
      use b2mod_types
      use b2mod_constants 
      use b2mod_plasma
      use b2mod_geo
      use b2mod_geo_corner
      use b2mod_indirect
      use b2mod_exp_file_names
      implicit none
      integer nx, ny
	integer ix, iy, ic
      real (kind=R8) :: power

*   ..for inner lower
      ix = -1
      fnm = 'mod_'//trim(fnm_ne_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), ne(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), te(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), ti(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_inner_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), 
     &     qe*ne(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*
     &     sqrt((te(ix,iy)+ti(ix,iy))/(mp)) 
	 enddo
	close(70)

*   ..for inner upper
      fnm = 'mod_'//trim(fnm_ne_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) then
	     write(70) crdy_c_cell(ix,iy), ne(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) then
	     write(70) crdy_c_cell(ix,iy), te(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) then
	     write(70) crdy_c_cell(ix,iy), ti(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_inner_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(rightix(ix,iy).eq.nx+1 .and. ix.lt.nx) then
	     write(70) crdy_c_cell(ix,iy), 
     &     qe*ne(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*
     &     sqrt((te(ix,iy)+ti(ix,iy))/(mp)) 
	    endif
	   enddo
	 enddo
	close(70)

*   ..for outer midpane
      fnm = 'mod_'//trim(fnm_ne_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(nmdpl,iy), ne(nmdpl,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(nmdpl,iy), te(nmdpl,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_mdpl)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(nmdpl,iy), ti(nmdpl,iy)
	 enddo
	close(70)

*   ..for outer upper
      fnm = 'mod_'//trim(fnm_ne_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(leftix(ix,iy).eq.-2 .and. ix.gt.-1) then
	     write(70) crdy_c_cell(ix,iy), ne(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(leftix(ix,iy).eq.-2 .and. ix.gt.-1) then
	     write(70) crdy_c_cell(ix,iy), te(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(leftix(ix,iy).eq.-2 .and. ix.gt.-1) then
	     write(70) crdy_c_cell(ix,iy), ti(ix,iy)
	    endif
	   enddo
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_outer_upper)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do ix = -1,nx
	   do iy = -1,ny
	    if(leftix(ix,iy).eq.-2 .and. ix.gt.-1) then
	     write(70) crdy_c_cell(ix,iy), 
     &     qe*ne(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*
     &     sqrt((te(ix,iy)+ti(ix,iy))/(mp)) 
	    endif
	   enddo
	 enddo
	close(70)

*   ..for outer lower
      ix = nx
      fnm = 'mod_'//trim(fnm_ne_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), ne(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_te_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), te(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_ti_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), ti(ix,iy)
	 enddo
	close(70)
      fnm = 'mod_'//trim(fnm_js_outer_lower)//'.dat'
      open(70,file=fnm,form='unformatted')
	 do iy=-1,ny
	  write(70) crdy_c_cell(ix,iy), 
     &     qe*ne(ix,iy)*bb(ix,iy,0)/bb(ix,iy,3)*
     &     sqrt((te(ix,iy)+ti(ix,iy))/(mp)) 
	 enddo
	close(70)

*   ..calculate total energy flux
	iy = 0
      power=0.0_R8
      do ix = -1, nx
	  if(mod(region(ix,iy,0),4).eq.1) then                ! for CORE
	    power = power + (fhe(ix,iy,1)+fhi(ix,iy,1))
	  endif
	enddo
      fnm = 'mod_'//trim(fnm_power)//'.dat'
	open(70,file=fnm,form='unformatted')
	 write(70) power
	close(70)

	end subroutine write_results
