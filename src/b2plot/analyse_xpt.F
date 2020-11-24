      subroutine analyse_xpt(lun)

c  version : 20.06.2005 02:38

c=====================================================
c*** analyse_xpt: analysis of the particle balance around the x-point
c***
c*** Input:
c***
c***  b2.parameters and B2SXDR files
c***
c*** Points around the lower x-point are considered, that is
c***  [nxcut1(1)+1,nycut1(1)],    [nxcut1(2),nycut1(1)]
c***  [nxcut1(1)+1,nycut1(1)+1],  [nxcut1(1),nycut1(1)+1]
c***  [nxcut1(1),nycut1(1)],      [nxcut1(2)+1,nycut1(1)]
c***  [nxcut1(2)+1,nycut1(1)+1],  [nxcut1(2),nycut1(1)+1]
c=====================================================
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_sources
      use b2mod_indirect
      implicit none
      integer lun
      integer n_points
      parameter (n_points=8)
      integer i,j,k,is
      integer ix(n_points),iy(n_points), ixx(n_points)
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

      ix(1)=leftcut(1)-1
      iy(1)=topcut(1)
      ix(2)=leftcut(1)
      iy(2)=topcut(1)
      ix(3)=leftcut(1)-1
      iy(3)=topcut(1)+1
      ix(4)=leftcut(1)
      iy(4)=topcut(1)+1
      ix(5)=rightcut(1)-1
      iy(5)=topcut(1)
      ix(6)=rightcut(1)
      iy(6)=topcut(1)
      ix(7)=rightcut(1)-1
      iy(7)=topcut(1)+1
      ix(8)=rightcut(1)
      iy(8)=topcut(1)+1
      ixx(1)=ix(2)
      ixx(2)=ix(2)-1
      ixx(3)=ix(4)
      ixx(4)=ix(4)-1
      ixx(5)=ix(5)-1
      ixx(6)=ix(5)
      ixx(7)=ix(8)
      ixx(8)=ix(8)-1

      allocate(ddt_ni(nspecies))
      allocate(div_nt(nspecies),div_nx(nspecies),div_ny(nspecies))
      allocate(f_nxp(nspecies),f_nxm(nspecies))
      allocate(f_nyp(nspecies),f_nym(nspecies))
      allocate(nt(nspecies),nt0(nspecies))
      allocate(s_ni(nspecies),s_io(nspecies),s_nt(nspecies))
      allocate(rel(nspecies),scl(nspecies),tot(nspecies))
      open(lun,file='x-point.dat')

c*** Calculate the balances for each cell in question
      
      write(lun,*) 'time step= ',dtim
      do i=1,n_points !{
        write(lun,'(/1x,a,i2,a,3i5,1p,3e14.6)') 'Point ',i,
     ,                ',  ix,iy,ixx,R,Z,V=',ix(i),iy(i),ixx(i),
     ,                  cr(ix(i),iy(i)),cz(ix(i),iy(i)),vol(ix(i),iy(i))
        write(frm,'(a,i2,a)') '(6x,',nspecies,'(i8,6x))'
        print *,'header format: ',frm !###
        write(lun,frm) (j,j=1,nspecies)
        write(frm,'(a,i2,a)') '(a5,1x,1p,',nspecies,'e14.6)'
        print *,'table format: ',frm !###
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
            s_ni(j)=s_ni(j)+b2stbr_sna(ix(i),iy(i),is)
            f_nxp(j)=f_nxp(j)+fna(ix(i),iy(i),0,is)
            f_nxm(j)=f_nxm(j)+fna(ixx(i),iy(i),0,is)
            f_nyp(j)=f_nyp(j)+fna(ix(i),iy(i),1,is)
            f_nym(j)=f_nym(j)+fna(ix(i),iy(i)-1,1,is)
            fl_nx=fl_nx+abs(fna(ix(i),iy(i),0,is))+
     +                  abs(fna(ixx(i),iy(i),0,is))
            fl_ny=fl_ny+abs(fna(ix(i),iy(i),1,is))+
     +                  abs(fna(ix(i),iy(i)-1,1,is))
            nt(j)=nt(j)+na(ix(i),iy(i),is)
            nt0(j)=nt0(j)+na0(ix(i),iy(i),is)
            s_io(j)=s_io(j)+rsana(ix(i),iy(i),is)
          end do !}
          is=is+1
          nt(j)=nt(j)*vol(ix(i),iy(i))
          nt0(j)=nt0(j)*vol(ix(i),iy(i))
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

      return
      end subroutine analyse_xpt
