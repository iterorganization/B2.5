! Fortran part of the adaptor for paraview catalyst
! for b2.5 simulation.
! Author: Jure Bartol
! Created on: 22.07.2016
! Modified on: 23.08.2016


subroutine coprocessor(crx,cry,ncrx,nx,ny,ns,step,time,vol,hx,hy,qc,te,ti,po,bzb,OnedBsq,qz,pbs,fhe,fhi,fch,pbshz,fhe_mdf,fhi_mdf,fchvispar,fchvisq,fchinert,fchdia,fchin,fch_p,fchvisper,gs,bb,na,ua,kinrgy,rra,rqa,fna,fna_mdf,fna_fcor,uadia,vadia,vaecrb,rlsa,rlra,rlqa,rlza,rlpt,rlpi)
!Coprocessor is called each time step in b2mndr to output
!simulation data to ParaView Catalyst.

  implicit none
  integer :: nx, ny, ns, step, flag, numC
  integer, intent(in) :: ncrx
  real(kind=8) :: time
  real(kind=8), dimension(-1:nx,-1:ny,1) :: vol,hx,hy,qc,te,ti,po,bzb,OnedBsq
  real(kind=8), dimension(-1:nx,-1:ny,2) :: qz,pbs,fhe,fhi,fch,pbshz,fhe_mdf,fhi_mdf,fchvispar,fchvisq,fchinert,fchdia,fchin,fch_p,fchvisper
  real(kind=8), dimension(-1:nx,-1:ny,3) :: gs
  real(kind=8), dimension(-1:nx,-1:ny,4) :: crx,cry,bb
  real(kind=8), dimension(-1:nx,-1:ny,0:ns-1) :: na,ua,kinrgy,rra,rqa,rsa
  real(kind=8), dimension(-1:nx,-1:ny,0:1,0:ns-1) :: fna,fna_mdf,fna_fcor,uadia,vadia,vaecrb,rlsa,rlra,rlqa,rlza,rlpt,rlpi

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
     call adddata(vol,"vol: Volume of the cell."//char(0),numC,1)
     call adddata(hx,"hx: Diameter of the cell in x-direction."//char(0),numC,1)
     call adddata(hy,"hy: Diameter of the cell in y-direction."//char(0),numC,1)
     call adddata(qc,"qc: Cos(t); t is the complementary angle between the x-direction of the cell and its left face. "//char(0),numC,1)
     call adddata(te,"te: Electron temperature on the cell."//char(0),numC,1)
     call adddata(ti,"ti: All atom temperature on the cell."//char(0),numC,1)
     call adddata(po,"po: Electric potential of the cell."//char(0),numC,1)
     call adddata(bzb,"bzb: Bz/B**2 - Bz/<B**2> on the cell centre."//char(0),numC,1)
     call adddata(OnedBsq,"OnedBsq: 1/bb(,,3)**2 on the cell centre"//char(0),numC,1)
! 2 components in each cell
     call adddata(qz,"qz: Sin(t) {x} and cos(t) {y}; t is the angle between the x-direction and the y-direction on the cell."//char(0),numC,2)
     call adddata(pbs,"pbs: Parallel contact area, on the cell left{x}/bottom{y} faces."//char(0),numC,2)
     call adddata(fhe,"fhe: Electron heat flux through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fhi,"fhi: All atom heat flux through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fch,"fch: Electric current through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(pbshz,"pbshz: Product hz*(bx/bb)*sx {x}  and hz*(by/bb)*sy {y}, or (magnetic field pitch)*(surface area)*hz."//char(0),numC,2)
     call adddata(fhe_mdf,"fhe_mdf: Modified el. heat flux through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fhi_mdf,"fhi_mdf: Modified all atom heat flux through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fchvispar,"fchvispar: El. current driven by parallel viscosity through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fchvisq,"fchvisq: El. current produced by the viscosity tensor through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fchinert,"fchinert: Contributions from inertia and gyroviscosity current through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fchdia,"fchdia: Modified diamagnetic current through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fchin,"fchin: Contribution from ion-neutral friction through the face between the cell and left{x}/bottom{y} neighbor."//char(0),numC,2)
     call adddata(fch_p,"fch_p: Product of the paral. el. current and pol. mag. field component sign through the face between the cell and left neighbor{x}. {y}=0"//char(0),numC,2)
     call adddata(fchvisper,"fchvisper: El. current connected with contribution from the perpend. viscosity through the face between the cell and bottom neighbor{y}. {x}=0"//char(0),numC,2)
! 3 components in each cell
     call adddata(gs,"gs: Area of the face between the cell and its left{x}/bottom{y}/third (ignorable) coordinate{z} neighbor."//char(0),numC,3) 
     call adddata(bb,"bb: Magnetic field (x,y,z component and magnitude) at the center of the cell."//char(0),numC,3)
! ns components
     call adddata(na,"na: Density of atomic species (each component - one species)."//char(0),numC,ns)
     call adddata(ua,"ua: Parallel velocity of atomic species (each component - one species) on the cell."//char(0),numC,ns)
     call adddata(kinrgy,"kinrgy: Kinetic energy of a particle of species (each component - one species) in the cell."//char(0),numC,ns)
     call adddata(rra,"rra: Rate coefficient for recombination from atomic species (each component - one species)."//char(0),numC,ns)
     call adddata(rqa,"rqa: Rate coefficient for el. heat loss terms (each component - one species)."//char(0),numC,ns)
     call adddata(rsa,"rsa: Rate coefficient for ionisation from atomic species (is) to (is+1) (each component - one species)."//char(0),numC,ns)
     call adddata(fna(:,:,0,:),"fnax: Flux of atoms of species through the face between the cell and its left neighbor."//char(0),numC,ns)
     call adddata(fna(:,:,1,:),"fnay: Flux of atoms of species through the face between the cell and its bottom neighbor."//char(0),numC,ns)
    call adddata(fna_mdf(:,:,0,:),"fna_mdfx: Modified flux of atoms of species through the face between the cell and its left neighbor."//char(0),numC,ns)
    call adddata(fna_mdf(:,:,1,:),"fna_mdfy: Modified flux of atoms of species through the face between the cell and its bottom neighbor."//char(0),numC,ns)
    call adddata(fna_fcor(:,:,0,:),"fna_fcorx: Flux of atoms of species for calculating momentum transport through the face between the cell and its left neighbor."//char(0),numC,ns)
     call adddata(fna_fcor(:,:,1,:),"fna_fcory: Flux of atoms of species for calculating momentum transport through the face between the cell and its left neighbor."//char(0),numC,ns)
     call adddata(uadia(:,:,0,:),"uadiax: Total effective drift velocity of species in the cell in the diamagnetic direction."//char(0),numC,ns)
     call adddata(uadia(:,:,1,:),"uadiay: Total effective drift velocity of species in the cell in the radial direction."//char(0),numC,ns)
     call adddata(vadia(:,:,0,:),"vadiax: Effective diamagnetic drift velocity of species in the cell in the diamagnetic direction."//char(0),numC,ns)
     call adddata(vadia(:,:,1,:),"vadiay: Effective diamagnetic drift velocity of species in the cell in the radial direction."//char(0),numC,ns)
     call adddata(vaecrb(:,:,0,:),"vaecrbx: ExB drift velocity of species in the cell in the diamagnetic direction."//char(0),numC,ns)
     call adddata(vaecrb(:,:,1,:),"vaecrby: ExB drift velocity of species in the cell in the radial direction."//char(0),numC,ns)
     call adddata(rlsa(:,:,0,:),"rlsa0: Coefficient in the equation rsa = exp(rlsa0+rlsa1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlsa(:,:,1,:),"rlsa1: Coefficient in the equation rsa = exp(rlsa0+rlsa1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlra(:,:,0,:),"rlra0: Coefficient in the equation rra = exp(rlra0+rlra1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlra(:,:,1,:),"rlra: Coefficient in the equation rra = exp(rlra0+rlra1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlqa(:,:,0,:),"rlqa0: Coefficient in the equation rqa = ev*exp(rlqa0+rlqa1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlqa(:,:,1,:),"rlqa: Coefficient in the equation rqa = ev*exp(rlqa0+rlqa1*log(te(ix,iy)/ev))."//char(0),numC,ns)
     call adddata(rlza(:,:,0,:),"rlza0: Coefficient in the equation rza = rlza0+rlza1*log(te(ix,iy)/ev)."//char(0),numC,ns)
     call adddata(rlza(:,:,1,:),"rlza1: Coefficient in the equation rza = rlza0+rlza1*log(te(ix,iy)/ev)."//char(0),numC,ns)
     call adddata(rlpt(:,:,0,:),"rlpt0: Coefficient in the equation rpt = rlpt0+rlpt1*log(te(ix,iy)/ev)."//char(0),numC,ns)
     call adddata(rlpt(:,:,1,:),"rlpt1: Coefficient in the equation rpt = rlpt0+rlpt1*log(te(ix,iy)/ev)."//char(0),numC,ns)
     call adddata(rlpi(:,:,0,:),"rlpi0: Coefficient in the equation rpi = rlpi0+rlpi1*log(te(ix,iy)/ev)."//char(0),numC,ns)
     call adddata(rlpi(:,:,1,:),"rlpi1: Coefficient in the equation rpi = rlpi0+rlpi1*log(te(ix,iy)/ev)."//char(0),numC,ns)

     call coprocess()
  end if
end subroutine

!!!Local Variables:
!!! mode: f90
!!! End:
