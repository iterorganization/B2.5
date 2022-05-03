! Fortran part of the adaptor for paraview catalyst
! for B2.5 simulation.
! Author: Jure Bartol


subroutine coprocessor(crx,cry,ncrx,nx,ny,ns,step,time,vol,hx,hy,qc,te,ti,po,&
     &                 bzb,OnedBsq,ne,qz,pbs,fhe,fhi,fch,pbshz,fhe_mdf,&
     &                 fhi_mdf,fchvispar,fchvisq,fchinert,fchdia,fchin,fch_p,&
     &                 fchvisper,gs,bb,na,ua,kinrgy,rra,rqa,fna,fna_mdf,&
     &                 fna_fcor,uadia,vadia,vaecrb,rlsa,rlra,rlqa,rlza,rlpt,&
     &                 rlpi, rightix, topiy)
!Coprocessor is called each time step in b2mndr_1 to output
!simulation data to ParaView Catalyst.

  use b2mod_types
  implicit none
  integer :: nx, ny, ns, step, flag, numC !, ix, iy
  integer, intent(in) :: ncrx
  real(kind=R8) :: time
  integer, dimension(-1:nx,-1:ny) :: rightix, topiy
  real(kind=R8), dimension(-1:nx,-1:ny) :: vol,hx,hy,te,ti,po,bzb,OnedBsq,ne
  real(kind=R8), dimension(-1:nx,-1:ny,0:1) :: qc,qz,pbs,pbshz
  real(kind=R8), dimension(-1:nx,-1:ny,0:1,0:1) :: fhe,fhi,fch,fhe_mdf,&
      & fhi_mdf,fchvispar,fchvisq,fchinert,fchdia,fchin,fch_p,fchvisper !,&
     !& velocity
  real(kind=R8), dimension(-1:nx,-1:ny,0:2) :: gs
  real(kind=R8), dimension(-1:nx,-1:ny,0:3) :: crx,cry,bb
  real(kind=R8), dimension(-1:nx,-1:ny,0:ns-1) :: na,ua,kinrgy,rra,rqa,rsa
  real(kind=R8), dimension(-1:nx,-1:ny,0:1,0:ns-1) :: rlsa,rlra,rlqa,&
      & rlza,rlpt,rlpi
  real(kind=R8), dimension(-1:nx,-1:ny,0:1,0:1,0:ns-1) :: fna,fna_mdf,fna_fcor,&
      & uadia,vadia,vaecrb
  real :: start, finish

  call cpu_time(start)

  !velocity(10,10,1) = 0
  !do ix = -1,nx
  !   do iy = -1,ny
  !      velocity(ix,iy,1) = 0.5 * (fna(ix,iy,0,0,0) / gs(ix,iy,1) + fna(rightix(ix,iy)+1,iy,0,0,0) / gs(rightix(ix,iy)+1,iy,1))
  !      velocity(ix,iy,2) = 0.5 * (fna(ix,iy,1,1,0) / gs(ix,iy,2) + fna(ix,topiy(ix,iy)+1,1,1,0) / gs(ix,topiy(ix,iy)+1,2))
  !   enddo
  !enddo
  !velocity(:,:,1) = velocity(:,:,1) / ne
  !velocity(:,:,2) = velocity(:,:,2) / ne

! Query Catalyst to see if there is something to do this time step.
  call requestdatadescription(step,time,flag)
  if (flag .ne. 0) then
     call needtocreategrid(flag)
     if (flag .ne. 0) then
        call creategrid(crx, cry, ncrx, nx, ny)
     end if

! Add field data to cells.
numC=(nx+2)*(ny+2) !number of cells
! 1 component in each cell
     call adddata(vol,"vol"//char(0),numC,1)
     call adddata(hx,"hx"//char(0),numC,1)
     call adddata(hy,"hy"//char(0),numC,1)
     call adddata(ti,"ti"//char(0),numC,1)
     call adddata(te,"te"//char(0),numC,1)
     call adddata(po,"po"//char(0),numC,1)
     call adddata(bzb,"bzb"//char(0),numC,1)
     call adddata(OnedBsq,"OnedBsq"//char(0),numC,1)
     call adddata(ne,"ne"//char(0),numC,1)
     !call adddata(velocity,"velocity"//char(0),numC,2)
! 2 components in each cell
     call adddata(qc,"qc"//char(0),numC,2)
     call adddata(qz,"qz"//char(0),numC,2)
     call adddata(pbs,"pbs"//char(0),numC,2)
     call adddata(pbshz,"pbshz"//char(0),numC,2)
! 3 components in each cell
     call adddata(gs,"gs"//char(0),numC,3)
     call adddata(bb,"bb"//char(0),numC,3)
! 2*2 components in each cell
     call adddata(fhe,"fhe"//char(0),numC,4)
     call adddata(fhi,"fhi"//char(0),numC,4)
     call adddata(fch,"fch"//char(0),numC,4)
     call adddata(fhe_mdf,"fhe_mdf"//char(0),numC,2)
     call adddata(fhi_mdf,"fhi_mdf"//char(0),numC,2)
     call adddata(fchvispar,"fchvispar"//char(0),numC,2)
     call adddata(fchvisq,"fchvisq"//char(0),numC,2)
     call adddata(fchinert,"fchinert"//char(0),numC,2)
     call adddata(fchdia,"fchdia"//char(0),numC,2)
     call adddata(fchin,"fchin"//char(0),numC,2)
     call adddata(fch_p,"fch_p"//char(0),numC,2)
     call adddata(fchvisper,"fchvisper"//char(0),numC,2)
! ns components
     call adddata(na,"na"//char(0),numC,ns)
     call adddata(ua,"ua"//char(0),numC,ns)
     call adddata(kinrgy,"kinrgy"//char(0),numC,ns)
     call adddata(rra,"rra"//char(0),numC,ns)
     call adddata(rqa,"rqa"//char(0),numC,ns)
     call adddata(rsa,"rsa"//char(0),numC,ns)
     call adddata(rlsa(:,:,0,:),"rlsa0"//char(0),numC,ns)
     call adddata(rlsa(:,:,1,:),"rlsa1"//char(0),numC,ns)
     call adddata(rlra(:,:,0,:),"rlra0"//char(0),numC,ns)
     call adddata(rlra(:,:,1,:),"rlra1"//char(0),numC,ns)
     call adddata(rlqa(:,:,0,:),"rlqa0"//char(0),numC,ns)
     call adddata(rlqa(:,:,1,:),"rlqa1"//char(0),numC,ns)
     call adddata(rlza(:,:,0,:),"rlza0"//char(0),numC,ns)
     call adddata(rlza(:,:,1,:),"rlza1"//char(0),numC,ns)
     call adddata(rlpt(:,:,0,:),"rlpt0"//char(0),numC,ns)
     call adddata(rlpt(:,:,1,:),"rlpt1"//char(0),numC,ns)
     call adddata(rlpi(:,:,0,:),"rlpi0"//char(0),numC,ns)
     call adddata(rlpi(:,:,1,:),"rlpi1"//char(0),numC,ns)
! 2*ns components
     call adddata(fna(:,:,0,:,:),"fnax1"//char(0),numC,2*ns)
     call adddata(fna(:,:,1,:,:),"fnay1"//char(0),numC,2*ns)
     call adddata(fna_mdf(:,:,0,:,:),"fna_mdfx"//char(0),numC,2*ns)
     call adddata(fna_mdf(:,:,1,:,:),"fna_mdfy"//char(0),numC,2*ns)
     call adddata(fna_fcor(:,:,0,:,:),"fna_fcorx"//char(0),numC,2*ns)
     call adddata(fna_fcor(:,:,1,:,:),"fna_fcory"//char(0),numC,2*ns)
     call adddata(uadia(:,:,0,:,:),"uadiax"//char(0),numC,2*ns)
     call adddata(uadia(:,:,1,:,:),"uadiay"//char(0),numC,2*ns)
     call adddata(vadia(:,:,0,:,:),"vadiax"//char(0),numC,2*ns)
     call adddata(vadia(:,:,1,:,:),"vadiay"//char(0),numC,2*ns)
     call adddata(vaecrb(:,:,0,:,:),"vaecrbx"//char(0),numC,2*ns)
     call adddata(vaecrb(:,:,1,:,:),"vaecrby"//char(0),numC,2*ns)

     call coprocess()
     call cpu_time(finish)
     print*, "Coprocessing time: ",finish-start
  end if
end subroutine

!!!Local Variables:
!!! mode: f90
!!! End:
