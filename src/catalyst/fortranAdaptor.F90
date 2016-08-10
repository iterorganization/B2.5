! Fortran part of the adaptor for paraview catalyst for b2.5 simulation.
! Author: Jure Bartol
! Created on: 22.07.2016
! Modified on: 10.08.2016


subroutine coprocessor(crx,cry,ncrx,nx,ny,step,time,vol,hx,hy,qc,te,ti,po,qz,pbs,na,ua,uadia,fhe,fhi,fch,gs,bb)
!Coprocessor is called each time step in b2mndr.F in order to output
!simulation data to ParaView Catalyst.

  implicit none
  integer :: nx, ny, step, flag, numofcells
  integer, intent(in) :: ncrx
  real(kind=8) :: time
  real(kind=8), dimension(nx,ny,1) :: vol, hx, hy, qc, te, ti, po
  real(kind=8), dimension(nx,ny,2) :: qz, pbs, na, ua, uadia, fhe, fhi, fch
  real(kind=8), dimension(nx,ny,3) :: gs
  real(kind=8), dimension(nx,ny,4) :: crx, cry, bb

  numofcells=(nx+2)*(ny+2)

  call requestdatadescription(step,time,flag) 
  if (flag .ne. 0) then
     call needtocreategrid(flag)
     if (flag .ne. 0) then
        call creategrid(crx, cry, ncrx, nx, ny)
     end if
!TODO variables:  fna(4d matrix)
!1 component in each cell
     call adddata(vol,"vol"//char(0),numofcells,1)
     call adddata(hx,"hx"//char(0),numofcells,1)
     call adddata(hy,"hy"//char(0),numofcells,1)
     call adddata(qc,"qc"//char(0),numofcells,1)
     call adddata(te,"te"//char(0),numofcells,1)
     call adddata(ti,"ti"//char(0),numofcells,1)
     call adddata(po,"po"//char(0),numofcells,1)
!2 components in each cell
     call adddata(qz,"qz"//char(0),numofcells,2)
     call adddata(pbs,"pbs"//char(0),numofcells,2)
     call adddata(na,"na"//char(0),numofcells,2)
     call adddata(ua,"ua"//char(0),numofcells,2)
     call adddata(uadia,"uadia"//char(0),numofcells,2)
     call adddata(fhe,"fhe"//char(0),numofcells,2)
     call adddata(fhi,"fhi"//char(0),numofcells,2)
     call adddata(fch,"fch"//char(0),numofcells,2)
!3 components in each cell
     call adddata(gs,"gs"//char(0),numofcells,3) 
!4 components in each cell
     call adddata(crx,"crx"//char(0),numofcells,4)
     call adddata(cry,"cry"//char(0),numofcells,4)
!bb has 4 components, one of them is magnitude
!which can be ignored, because paraview calculates it
!on its own.
     call adddata(bb,"bb"//char(0),numofcells,3)
    
     call coprocess()
  end if
end subroutine

!!!Local Variables:
!!! mode: f90
!!! End:
