	subroutine my_outi(nwrite, nx, ny, iarr, first_name)
	 integer nwrite, nx, ny
	 integer iarr(-1:nx,-1:ny)
	 character *(*) first_name
	 character *64 filename
	 integer ix, iy
	 filename = trim(first_name)//'.dat'
       open (nwrite,file=filename)                         
	 write (nwrite,'(6x,<nx+2>(   i4   ))') (ix,ix=-1,nx)	 
	 write (nwrite,'(6x,<nx+2>(   i4   ))')	 
	 do 	iy=ny,-1,-1									 
	  write(nwrite,'(1x,i4,1x,<nx+2>i4)') iy,(iarr(ix,iy),ix=-1,nx) 
	 enddo											 
	 close (nwrite)
	end
