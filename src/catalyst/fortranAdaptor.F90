! Fortran part of the adaptor for paraview catalyst for b2.5 simulation.
! Author: Jure Bartol
! Created on: 22.07.2016
! Modified on: 10.08.2016


subroutine coprocessor(crx,cry,ncrx,nx,ny,ns,step,time,vol, hx, hy, qc, te, ti, po, bzb, OnedBsq, qz, pbs, fhe, fhi,fch, pbshz, fhe_mdf, fhi_mdf, fchvispar,fchvisq, fchinert,  fchdia, fchin, fch_p, fchvisper, gs, bb,na, ua, kinrgy, rra, rqa,fna, fna_mdf, fna_fcor, uadia, vadia, vaecrb, rlsa, rlra, rlqa, rlza, rlpt, rlpi)
!Coprocessor is called each time step in b2mndr.F in order to output
!simulation data to ParaView Catalyst.

  implicit none
  integer :: nx, ny, ns, step, flag, numC
  integer, intent(in) :: ncrx
  real(kind=8) :: time
  real(kind=8), dimension(nx,ny,1) :: vol, hx, hy, qc, te, ti, po, bzb, OnedBsq
  real(kind=8), dimension(nx,ny,2) :: qz, pbs, fhe, fhi, fch, pbshz, fhe_mdf, fhi_mdf, fchvispar, fchvisq, fchinert,  fchdia, fchin, fch_p, fchvisper
  real(kind=8), dimension(nx,ny,3) :: gs
  real(kind=8), dimension(nx,ny,4) :: crx, cry, bb
  real(kind=8), dimension(nx,ny,ns) :: na, ua, kinrgy, rra, rqa, rsa
  real(kind=8), dimension(nx,ny,2,ns) :: fna, fna_mdf, fna_fcor, uadia, vadia, vaecrb, rlsa, rlra, rlqa, rlza, rlpt, rlpi

!special: fna, fna_mdf, fna_fcor, uadia, vadia, vaecrb, rlsa, rlra, rlqa, rlza, rlpt, rlpi


  numC=(nx+2)*(ny+2) !number of cells

  call requestdatadescription(step,time,flag)
  if (flag .ne. 0) then
     call needtocreategrid(flag)
     if (flag .ne. 0) then
        call creategrid(crx, cry, ncrx, nx, ny)
     end if

!1 component in each cell
!bzb, OnedBsq, nb, ub, ne, ne2, pr, pe, ro, ne2m, resmb, resmt, reshe, reshi, resht, respo, corut, corte, corti, cortt, corpo, hce0, hci0, sig0, alf0, ceqp,
     call adddata(vol,"vol"//char(0),numC,1)
     call adddata(hx,"hx"//char(0),numC,1)
     call adddata(hy,"hy"//char(0),numC,1)
     call adddata(qc,"qc"//char(0),numC,1)
     call adddata(te,"te"//char(0),numC,1)
     call adddata(ti,"ti"//char(0),numC,1)
     call adddata(po,"po"//char(0),numC,1)
     call adddata(bzb,"bzb"//char(0),numC,1)
     call adddata(OnedBsq,"OnedBsq"//char(0),numC,1)
  
!2 components in each cell
!pbs, pbshz, fhe, fhe_mdf, fhi, fhi_mdf, fch, fub, fub_mdf, fub_fcor, ni, fne, fni, fchvispar, fchvisq, fchinert, fchanml, fchdia, fchin, fch_p, fchvisper, vedia, veecrb, sne, sne0,  hce, hci, siq, alf, chce, chve, chci, chvi, csig, calf, csig_cl, calf_cl, csig_an, calf_an
     call adddata(qz,"qz"//char(0),numC,2)
     call adddata(pbs,"pbs"//char(0),numC,2)
     call adddata(fhe,"fhe"//char(0),numC,2)
     call adddata(fhi,"fhi"//char(0),numC,2)
     call adddata(fch,"fch"//char(0),numC,2)
     call adddata(pbshz,"pbshz"//char(0),numC,2)
     call adddata(fhe_mdf,"fhe_mdf"//char(0),numC,2)
     call adddata(fhi_mdf,"fhi_mdf,"//char(0),numC,2)
     call adddata(fchvispar,"fchvispar"//char(0),numC,2)
     call adddata(fchvisq,"fchvisq"//char(0),numC,2)
     call adddata(fchinert,"fchinert"//char(0),numC,2)
     call adddata(fchdia,"fchdia"//char(0),numC,2)
     call adddata(fchin,"fchin"//char(0),numC,2)
     call adddata(fch_p,"fch_p"//char(0),numC,2)
     call adddata(fchvisper,"fchvisper"//char(0),numC,2)
!3 components in each cell
!gs
     call adddata(gs,"gs"//char(0),numC,3) 
!4 components in each cell
 !bb,she, shi, sch, she0, shi0, sch0
!bb has 4 components, one of them is magnitude
!which can be ignored, because paraview calculates it
!on its own.
     call adddata(bb,"bb"//char(0),numC,3)
!ns components
!na, ua, kinrgy, resco, resmo, copra, corua, rra, rqa, rcx, dna0, dpa0, vsa0, rsa
!fna, fna_mdf, fna_fcor, uadia, vadia, vaecrb, rlsa, rlra, rlqa, rlza, rlpt, rlpi
     call adddata(na,"na"//char(0),numC,ns)
     call adddata(ua,"ua"//char(0),numC,ns)
     call adddata(kinrgy,"kinrgy"//char(0),numC,ns)
     call adddata(rra,"rra"//char(0),numC,ns)
     call adddata(rqa,"rqa"//char(0),numC,ns)
     call adddata(rsa,"rsa"//char(0),numC,ns)
     call adddata(fna(:,:,0,:),"fnax"//char(0),numC,ns)
     call adddata(fna(:,:,1,:),"fnay"//char(0),numC,ns)
     call adddata(fna_mdf(:,:,0,:),"fna_mdfx"//char(0),numC,ns)
     call adddata(fna_mdf(:,:,1,:),"fna_mdfy"//char(0),numC,ns)
     call adddata(fna_fcor(:,:,0,:),"fna_fcorx"//char(0),numC,ns)
     call adddata(fna_fcor(:,:,1,:),"fna_fcory"//char(0),numC,ns)
     call adddata(uadia(:,:,0,:),"uadiax"//char(0),numC,ns)
     call adddata(uadia(:,:,1,:),"uadiay"//char(0),numC,ns)
     call adddata(vadia(:,:,0,:),"vadiax"//char(0),numC,ns)
     call adddata(vadia(:,:,1,:),"vadiay"//char(0),numC,ns)
     call adddata(vaecrb(:,:,0,:),"vaecrbx"//char(0),numC,ns)
     call adddata(vaecrb(:,:,1,:),"vaecrby"//char(0),numC,ns)
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

     call coprocess()
  end if
end subroutine

!!!Local Variables:
!!! mode: f90
!!! End:
