! Fortran part of the adaptor for paraview catalyst
! for B2.5 simulation.
! Author: Jure Bartol


subroutine coprocessor(nVx, nFc, nCv, step, time, geo, state)
!Coprocessor is called each time step in b2mndr_1 to output
!simulation data to ParaView Catalyst.

  use b2mod_types
  use b2us_geo_diff
  use b2us_plasma_diff
#ifdef _OPENMP
  use b2mod_openmp
#endif

  implicit none
  integer, intent(in) :: nVx, nFc, nCv, step
  real(kind=R8), intent(in) :: time
  type(geometry), intent(in) :: geo
  type(B2State), intent(in) :: state
  integer :: flag
#ifdef _OPENMP
  real(kind=R8) :: start, finish
#else
  real(kind=R4) :: start, finish
#endif

#ifdef _OPENMP
  start = omp_get_wtime()
#else
  call cpu_time(start)
#endif

! Query Catalyst to see if there is something to do this time step.
  call requestdatadescription(step,time,flag)
  if (flag .ne. 0) then
     call needtocreategrid(flag)
     if (flag .ne. 0) then
        call creategrid(geo%vxX, geo%vxY, nVx)
     end if

! Add field data to cells.
! 1 component in each cell
     call adddata(geo%cvVol,"vol"//char(0),nCv,1)
     call adddata(geo%cvHx,"hx"//char(0),nCv,1)
     call adddata(state%pl%ti,"ti"//char(0),nCv,1)
     call adddata(state%pl%te,"te"//char(0),nCv,1)
     call adddata(state%pl%po,"po"//char(0),nCv,1)
     call adddata(geo%cvBzb,"bzb"//char(0),nCv,1)
     call adddata(geo%cvOnedBsq,"OnedBsq"//char(0),nCv,1)
     call adddata(state%dv%ne,"ne"//char(0),nCv,1)
     call adddata(geo%fcPbs,"pbs"//char(0),nFc,1)
     call adddata(geo%fcPbshz,"pbshz"//char(0),nFc,1)
     call adddata(geo%fcS,"gs"//char(0),nFc,1)
     call adddata(geo%cvSz,"sz"//char(0),nCv,1)
! 2 components in each cell
     call adddata(geo%fcQgam,"qc"//char(0),nFc,2)
     call adddata(geo%cvQgam,"qz"//char(0),nCv,2)
     call adddata(state%dv%fhe,"fhe"//char(0),nFc,2)
     call adddata(state%dv%fhi,"fhi"//char(0),nFc,2)
     call adddata(state%dv%fch,"fch"//char(0),nFc,2)
     call adddata(state%dv%fhe_mdf,"fhe_mdf"//char(0),nFc,2)
     call adddata(state%dv%fhi_mdf,"fhi_mdf"//char(0),nFc,2)
     call adddata(state%dv%fchvispar,"fchvispar"//char(0),nFc,2)
     call adddata(state%dv%fchvisq,"fchvisq"//char(0),nFc,2)
     call adddata(state%dv%fchinert,"fchinert"//char(0),nFc,2)
     call adddata(state%dv%fchdia,"fchdia"//char(0),nFc,2)
     call adddata(state%dv%fchin,"fchin"//char(0),nFc,2)
     call adddata(state%dv%fch_p,"fch_p"//char(0),nFc,2)
     call adddata(state%dv%fchvisper,"fchvisper"//char(0),nFc,2)
! 3 components in each cell
     call adddata(geo%cvBb(:,0:2),"bb"//char(0),nCv,3)
! ns components
     call adddata(state%pl%na,"na"//char(0),nCv,state%ns)
     call adddata(state%pl%ua,"ua"//char(0),nCv,state%ns)
     call adddata(state%dv%kinrgy,"kinrgy"//char(0),nCv,state%ns)
!WG_TODO     call adddata(state%rt%rra,"rra"//char(0),nCv,state%ns)
!WG_TODO     call adddata(state%rt%rqa,"rqa"//char(0),nCv,state%ns)
!WG_TODO     call adddata(state%rt%rsa,"rsa"//char(0),nCv,state%ns)
     call adddata(state%dv%fna(:,0,:),"fnax"//char(0),nFc,state%ns)
     call adddata(state%dv%fna(:,1,:),"fnay"//char(0),nFc,state%ns)
     call adddata(state%dv%fna_mdf(:,0,:),"fna_mdfx"//char(0),nFc,state%ns)
     call adddata(state%dv%fna_mdf(:,1,:),"fna_mdfy"//char(0),nFc,state%ns)
     call adddata(state%dv%fna_fcor(:,0,:),"fna_fcorx"//char(0),nFc,state%ns)
     call adddata(state%dv%fna_fcor(:,1,:),"fna_fcory"//char(0),nFc,state%ns)
     call adddata(state%dv%uadia(:,0,:),"uadiax"//char(0),nFc,state%ns)
     call adddata(state%dv%uadia(:,1,:),"uadiay"//char(0),nFc,state%ns)
     call adddata(state%dv%vadia(:,0,:),"vadiax"//char(0),nFc,state%ns)
     call adddata(state%dv%vadia(:,1,:),"vadiay"//char(0),nFc,state%ns)
     call adddata(state%dv%vaecrb(:,0,:),"vaecrbx"//char(0),nFc,state%ns)
     call adddata(state%dv%vaecrb(:,1,:),"vaecrby"//char(0),nFc,state%ns)
     call adddata(state%rt%rlsa(:,0,:),"rlsa0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlsa(:,1,:),"rlsa1"//char(0),nCv,state%ns)
     call adddata(state%rt%rlra(:,0,:),"rlra0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlra(:,1,:),"rlra1"//char(0),nCv,state%ns)
     call adddata(state%rt%rlqa(:,0,:),"rlqa0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlqa(:,1,:),"rlqa1"//char(0),nCv,state%ns)
     call adddata(state%rt%rlza(:,0,:),"rlza0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlza(:,1,:),"rlza1"//char(0),nCv,state%ns)
     call adddata(state%rt%rlpt(:,0,:),"rlpt0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlpt(:,1,:),"rlpt1"//char(0),nCv,state%ns)
     call adddata(state%rt%rlpi(:,0,:),"rlpi0"//char(0),nCv,state%ns)
     call adddata(state%rt%rlpi(:,1,:),"rlpi1"//char(0),nCv,state%ns)

     call coprocess()
#ifdef _OPENMP
     finish = omp_get_wtime()
#else
     call cpu_time(finish)
#endif
     print*, "Coprocessing time: ",finish-start
  end if
  return

end subroutine coprocessor

!!!Local Variables:
!!! mode: f90
!!! End:
