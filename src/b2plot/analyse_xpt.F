      subroutine analyse_xpt(lun)

c  version : 20.06.2005 02:38

c=====================================================
c*** analyse_xpt: analysis of the particle balance around the x-point
c***
c*** Input:
c***
c***  b2.parameters and B2SXDR files
c***
c*** We consider here the points of control volumes with the X-point
c*** as a vertex.
c*** jx+ fluxes are those of radial faces ending or starting at the X-point (if any)
c*** jx- fluxes are those through the opposite faces radial (if any)
c*** jy+ fluxes are those across the separatrix segment (if any)
c*** jy- fluxes are those through the opposite faces (if any)
c=====================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2us_map
      use b2us_geo
      use b2us_data
      implicit none
      integer lun
      integer n_points
      integer i,j,k,is,iVx,iFc,iCv1,iCv2
      integer, allocatable :: iCx(:),iFxm(:),iFxp(:),iFym(:),iFyp(:)
      real (kind=R8) :: fl_nx,fl_ny,fl_nt,dtim
      real (kind=R8), allocatable ::
     ,  div_nx(:),div_ny(:),div_nt(:),ddt_ni(:),tot(:),rel(:),
     ,  scl(:),s_ni(:),s_io(:),s_nt(:),f_nxp(:),f_nyp(:),
     ,  f_nxm(:),f_nym(:),nt0(:),nt(:)
      character hlp(16)*4,frm*24
      data dtim /1.0_R8/
      save dtim
      data hlp /'tot','rel','scl','d/dt','Sni','div','divx',
     ,  'divy','ni ','ni0','Snn','Sio','jx+','jx-','jy+','jy-'/
c=====================================================
c
      call ipgetr('b2mndr_dtim',dtim)
c*** Set the indices of the cells in question

      allocate(ddt_ni(nspecies))
      allocate(div_nt(nspecies),div_nx(nspecies),div_ny(nspecies))
      allocate(f_nxp(nspecies),f_nxm(nspecies))
      allocate(f_nyp(nspecies),f_nym(nspecies))
      allocate(nt(nspecies),nt0(nspecies))
      allocate(s_ni(nspecies),s_io(nspecies),s_nt(nspecies))
      allocate(rel(nspecies),scl(nspecies),tot(nspecies))
      open(lun,file='x-point.dat')
      write(lun,*) 'time step= ',dtim

c*** Calculate the balances for each cell in question
      
      iVx = mpg%Xpt(1)
      n_points = mpg%vxCvP(iVx,2)
      allocate(iCx(n_points),
     .         iFxm(n_points),iFxp(n_points),
     .         iFym(n_points),iFyp(n_points))
      do i=1,n_points !{
        iCx(i) = mpg%vxCv(mpg%vxCvP(iVx,1)+i-1)
        do j=1,mpg%cvFcP(iCx(i),2)
          iFc=mpg%cvFc(mpg%cvFcP(iCx(i),1)+j-1)
          iCv1=mpg%fcCv(iFc,1)
          iCv2=mpg%fcCv(iFc,2)
          if (mpg%fcFs(iFc).eq.mpg%iFssep.or.
     ,       (mpg%fcFs(iFc).eq.mpg%iFssep2.and.mpg%iFssep2.ne.0)) then
            iFyp(i)=iFc
          else if (mpg%fcFs(iFc).ne.0) then
            iFym(i)=iFc
          else if (mpg%cvReg(iCv1).ne.mpg%cvReg(iCv2)) then
            iFxp(i)=iFc
          else if (.not.mpg%cvOnClosedSurface(iCv1).and.
     ,             .not.mpg%cvOnClosedSurface(iCv2)) then
            iFxm(i)=iFc
          else if (mpg%fcVx(iFc,1).eq.iVx.or.mpg%fcVx(iFc,2).eq.iVx)
     >     then
            iFxp(i)=iFc
          else
            iFxm(i)=iFc
          end if
        end do
        write(lun,'(/1x,a,i2,a,i5,1p,3e14.6)')
     ,   'Point ',i,',  iCx,R,Z,V=',iCx(i),
     ,    geo%cvX(iCx(i)),geo%cvY(iCx(i)),geo%cvVol(iCx(i))
        write(frm,'(a,i2,a)') '(6x,',nspecies,'(i8,6x))'
        if (i.eq.1) print *,'header format: ',frm !###
        write(lun,frm) (j,j=1,nspecies)
        write(frm,'(a,i2,a)') '(a5,1x,1p,',nspecies,'e14.6)'
        if (i.eq.1) print *,'table format: ',frm !###
        is=0
        do j=1,nspecies !{
          ddt_ni(j)=0.0_R8
          div_nx(j)=0.0_R8
          div_ny(j)=0.0_R8
          f_nxp(j)=0.0_R8
          f_nxm(j)=0.0_R8
          f_nyp(j)=0.0_R8
          f_nym(j)=0.0_R8
          fl_nx=0.0_R8
          fl_ny=0.0_R8
          fl_nt=0.0_R8
          nt(j)=0.0_R8
          nt0(j)=0.0_R8
          s_ni(j)=0.0_R8
          s_io(j)=0.0_R8
          do k=1,nfluids(j) !{
            is=is+1
            s_ni(j)=s_ni(j)+state%srw%b2stbr_sna(iCx(i),is)
            if (iFxp(i).ne.0)
     .       f_nxp(j)=f_nxp(j)+state%dv%fna(iFxp(i),0,is)
            if (iFxm(i).ne.0)
     .       f_nxm(j)=f_nxm(j)+state%dv%fna(iFxm(i),0,is)
            if (iFyp(i).ne.0)
     .       f_nyp(j)=f_nyp(j)+state%dv%fna(iFyp(i),1,is)
            if (iFym(i).ne.0)
     .       f_nym(j)=f_nym(j)+state%dv%fna(iFym(i),1,is)
            if (iFxp(i).ne.0)
     .       fl_nx=fl_nx+abs(state%dv%fna(iFxp(i),0,is))
            if (iFxm(i).ne.0)
     .       fl_nx=fl_nx+abs(state%dv%fna(iFxm(i),0,is))
            if (iFyp(i).ne.0)
     .       fl_ny=fl_ny+abs(state%dv%fna(iFyp(i),1,is))
            if (iFym(i).ne.0)
     .       fl_ny=fl_ny+abs(state%dv%fna(iFym(i),1,is))
            nt(j)=nt(j)+state%pl%na(iCx(i),is)
            nt0(j)=nt0(j)+state%psnl%na(iCx(i),is)
            s_io(j)=s_io(j)+state%srw%rsana(iCx(i),is)
          end do !}
          is=is+1
          nt(j)=nt(j)*geo%cvVol(iCx(i))
          nt0(j)=nt0(j)*geo%cvVol(iCx(i))
          ddt_ni(j)=(nt(j)-nt0(j))/dtim
          fl_nt=nt(j)
          s_nt(j)=s_ni(j)+s_io(j)
          div_nx(j)=f_nxp(j)-f_nxm(j)
          div_ny(j)=f_nyp(j)-f_nym(j)
          div_nt(j)=div_nx(j)+div_ny(j)
          tot(j)=ddt_ni(j)-s_nt(j)+div_nt(j)
          rel(j)=tot(j)/(fl_nt+fl_nx+fl_ny+abs(s_io(j))+abs(s_nt(j)))
          scl(j)=tot(j)/nt(j)
        end do !} 
        write(lun,frm) hlp(1),(tot(j),j=1,nspecies)
        write(lun,frm) hlp(2),(rel(j),j=1,nspecies)
        write(lun,frm) hlp(3),(scl(j),j=1,nspecies)
        write(lun,frm) hlp(4),(ddt_ni(j),j=1,nspecies)
        write(lun,frm) hlp(5),(s_nt(j),j=1,nspecies)
        write(lun,frm) hlp(6),(div_nt(j),j=1,nspecies)
        write(lun,frm) hlp(7),(div_nx(j),j=1,nspecies)
        write(lun,frm) hlp(8),(div_ny(j),j=1,nspecies)
        write(lun,frm) hlp(9),(nt(j),j=1,nspecies)
        write(lun,frm) hlp(10),(nt0(j),j=1,nspecies)
        write(lun,frm) hlp(11),(s_ni(j),j=1,nspecies)
        write(lun,frm) hlp(12),(s_io(j),j=1,nspecies)
        write(lun,frm) hlp(13),(f_nxp(j),j=1,nspecies)
        write(lun,frm) hlp(14),(f_nxm(j),j=1,nspecies)
        write(lun,frm) hlp(15),(f_nyp(j),j=1,nspecies)
        write(lun,frm) hlp(16),(f_nym(j),j=1,nspecies)
      end do !}
      close(lun)
      deallocate(ddt_ni,div_nx,div_ny,div_nt)
      deallocate(f_nxp,f_nxm)
      deallocate(f_nyp,f_nym)
      deallocate(nt,nt0)
      deallocate(rel,scl,tot)
      deallocate(s_ni,s_io,s_nt)
      deallocate(iCx,iFxm,iFxp,iFyp,iFym)

      return
      end subroutine analyse_xpt
