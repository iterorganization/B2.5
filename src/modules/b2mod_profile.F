      module b2mod_profile
      use b2mod_types
      use b2mod_boundary_namelist
#ifndef ALLOCATE_AND_NAMELIST
#     include "DIMENSIONS.F"
#endif
      contains
      subroutine read_profile(filename,bc1,bc2,pos_in,
     .     profile,ny,max_profile,scale)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     read_profile reads a profile from a file and interpolates this
*     profile (after normalization, if necessary) onto the cell centered positions into profile. 
*     This routine is developed to use a profile from a file as boundary
*     condition, which is for instance needed for the inlet of a linear
*     plasma generator.
*
*
*  3. description 
*
*     This routine reads in a profile from a file. This file contains
*     two columns, the position and its corresponding value.
*     This profile is interpolated onto the positions pos_in using E01BEF, 
*     and will be saved in profile.
*     If scale>0, the profile will be rescaled such that the maximum equals max_profile.
*
* E01BEF(nD,D_X,D_Y,D_D,ifail) fits a monotonicity-preserving, 
* piecewise cubic Hermitian function to the input coordinates 
* (D_X,D_Y) of length nD
* output: D_D, estimation of first derivative at coordinates (D_X,D_Y).
*
* E01BFF(nD,D_X,D_Y,D_D,nP,POS,result,ifail) uses input 
* D_X,D_y,D_D of length nD (prescribing the monotonicity-preserving, 
* piecewise cubic Hermitian function) to return nP values at x=POS(1:nP) 
* output: result(1:nP)
*      
*-----------------------------------------------------------------------
*.declarations
      implicit none
*input
      integer bc1, bc2, pos_in, ny
      real (kind=R8) :: max_profile, scale
*output
      real (kind=R8) :: profile(-1:ny)

*other
      integer iy, npoints, ifail
      real (kind=R8), allocatable :: deriv(:), pos(:), value(:)
      real (kind=R8) :: max, cntr(1), tmp(1)
      logical file_ok
      character*260 filename

*set initial values
      max = 0.0_R8
      ifail = 0

*open and read profile input file
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        read(99,*,err=100,end=101) npoints
*check input
        if (npoints.LE.3) then
          call XERRAB
     1          ('read_profile, more than 3 points needed for profile')
        endif
*allocation of variables for profile
        allocate(pos(1:npoints))
        allocate(value(1:npoints))
        allocate(deriv(1:npoints))

        do iy = 1, npoints
          read(99,*,err=100,end=101) pos(iy), value(iy)
        enddo
        close(99)
        max = maxval(value(1:npoints))

*interpolate to desired position and rescale
        call E01BEF(npoints,pos(1:npoints),value(1:npoints),
     .       deriv(1:npoints),ifail)
        do iy= bc1, bc2
          cntr(1)=center_of_cell(pos_in,iy,1)
          call E01BFF(npoints,pos,value,deriv,1,cntr,tmp,
     .         ifail)
          if (scale.gt.0.0_R8) then
c     scale the profile to max_profile
            profile(iy) = tmp(1) * max_profile / max
          else
            profile(iy) = tmp(1)
          endif
        enddo
      else
        write(*,*) ' filename for profile is', filename
        call XERRAB ('READ_PROFILE -- FILE DOES NOT EXIST')
      endif
      return
 100  call XERRAB('READ_PROFILE, FILE READING ERROR')
 101  call XERRAB('read_profile, reached EOF')
        
      end subroutine read_profile
      
      function center_of_cell(ix,iy,flag)
*****************************************************************
*     description
*     this function returns the center coordinate of cell (ix,iy)
*     flag=1 indicates that crx is used
*     flag=2 indicates that cry is used
*     So, for the linear case, flag=1 gives the r coordinate,
*     flag=2 gives the z coordinate.

      use b2mod_geo
      real (kind=R8) :: center_of_cell
      integer ix, iy, flag
      if (flag.eq.1) then
        center_of_cell=(crx(ix,iy,0)+
     .                  crx(ix,iy,1)+
     .                  crx(ix,iy,2)+
     .                  crx(ix,iy,3))/4
      elseif (flag.eq.2) then
        center_of_cell=(cry(ix,iy,0)+
     .                  cry(ix,iy,1)+
     .                  cry(ix,iy,2)+
     .                  cry(ix,iy,3))/4
      else
        write (*,*) 'incorrect value of flag in center_of_cell'
        CALL XERRAB('b2mod_profile, set value of flag to 1 or 2')
      endif

      return
      end function center_of_cell
      
      end module b2mod_profile
