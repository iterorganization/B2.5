! Fortran part of the adaptor for paraview catalyst for b2.5 simulation.
! Author: Jure Bartol
! Created on: 22.07.2016
! Modified on: 03.08.2016


subroutine coprocessorstart()
  character(len=200) :: arg
  integer :: i, ilen

!initialize with python script as an argument (this probably won't be possible). Instead we will use
! function coprocessorinitializewithpython("nameOfScript.py",<length of the file name>)
  call coprocessorinitialize() 
  do i=1, iargc()
     call getarg(i, arg)
     ilen = len_trim(arg)
     arg(ilen+1:) = char(0)
     call coprocessoraddpythonscript(arg, ilen)
  enddo
end subroutine


subroutine coprocessor(crx,cry,ncrx,nx,ny,step,time)
  use iso_c_binding ! IGNORE
  implicit none
  integer :: nx, ny, step, flag
  integer, intent(in) :: ncrx
  real :: time
  real(kind=8), dimension(nx,ny,4) :: crx, cry

  call requestdatadescription(step,time,flag) 
  if (flag .ne. 0) then
     call needtocreategrid(flag)
     if (flag .ne. 0) then
        call creategrid(crx, cry, ncrx, nx, ny)
     end if
     !call adddata(crx,"crx"//char(0),(nx+2)*(ny+2)) !//char(0) is a C++ terminating character in fortran
     !call adddata(cry,"cry"//char(0),(nx+2)*(ny+2))
     call coprocess()
  end if
end subroutine

!!!Local Variables:
!!! mode: f90
!!! End:
