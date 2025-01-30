      subroutine b2pplentd(lun)
c
c----------------------------------------------------------------------
c***  Produces a file with radiation emissivity mapped onto a
c***  rectangular grid - for 3D calculation of the wall loading
c***
c*** Input data - passed through the real stack
c***  11 numbers: a,b,c,d,e,f,g,h,i,j,r0rad
c***
c***  Q = A*exp(-j*z**2)*
c***   a + b*x**c + d*(1-x**e)**f     0 < x < 1  inside the separatrix
c***   a0 + a1*y + a2*y**2            0 < y < 1  outside the separatrix
c***       a0 = a + b
c***       a0 + a1 + a2 = g
c***       a0 + h*a1 + h**2*a2 = i
c***    that is,
c***      a+d is the emissivity at the plasma center
c***      e,f characterise radiation peak at the plasma center
c***      a+b is the emissivity at the separatrix
c***      c   characterises radiation peak just inside the separatrix
c***      g   is the emissivity at the grid edge
c***      h,i characterise radiation profile outside the separatrix
c***
c***   The normalisation constant A ensures the correct value of
c***   the integral: Qtot = r0rad
c***
c*** Output files:
c***
c***  plenter.dat
c----------------------------------------------------------------------
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_constants
      use b2mod_b2plot_ank
      use b2mod_user_namelist
      implicit none
      integer lun
c
c*** Rectangular grid:
      integer nnr,nnz
      real (kind=R8) :: rmn,rmx,zmn,zmx
      parameter (nnr=40,nnz=60)
      parameter (rmn=5.0_R8, rmx=11.0_R8, zmn=-3.0_R8, zmx=6.0_R8)
      integer nqch
c*** Grid sub-division for calculation of sources
      parameter (nqch=5)
c
c*** Local variables
      integer i,j,k,l,iCv,iFt,jmp,jtp
      real (kind=R8) :: x
      real (kind=R8) :: a, b, cc, d, e, f, a0, a1, a2
      real (kind=R8) :: sga,bgd,sgd,sgb
      real (kind=R8) :: initparm(10), r0rad
      real (kind=R8) :: q_tmp(nnr,nnz)
      real (kind=R8), allocatable :: qpr_tmp(:), qpz_tmp(:)
      real (kind=R8) :: r,z,bq,dr,dz,qmx,qmn,rmp,ztp,ztt,pwr_tot
      logical found
c----------------------------------------------------------------------
c*** Calculate the coefficients for radiation profiles
c
      allocate(qpr_tmp(1:DEF_NYD),qpz_tmp(1:DEF_NXD))
      stack_req=(/0,11,0,0,0,0,0/)
      if (bad_args(stack_req)) return
      r0rad=rpop(rstack,irstack)
      initparm(10)=rpop(rstack,irstack)
      initparm(9)=rpop(rstack,irstack)
      initparm(8)=rpop(rstack,irstack)
      initparm(7)=rpop(rstack,irstack)
      initparm(6)=rpop(rstack,irstack)
      initparm(5)=rpop(rstack,irstack)
      initparm(4)=rpop(rstack,irstack)
      initparm(3)=rpop(rstack,irstack)
      initparm(2)=rpop(rstack,irstack)
      initparm(1)=rpop(rstack,irstack)
      a=initparm(1)
      b=initparm(2)
      cc=initparm(3)
      d=initparm(4)
      e=initparm(5)
      f=initparm(6)
      a0=a+b
      bgd=initparm(7)-a0
      sgd=initparm(8)-a0
      sga=initparm(9)+0.5_R8
      sgb=initparm(10)
      a2=(sga*bgd-sgd)/(1.0_R8-sga)
      a1=bgd-a2
      pwr_tot=r0rad
c
c*** Calculate the radial profile qpr_tmp
c*** First, find the "outer midplane" ...
      rmp=0.0_R8
      jmp=0
      iFt=mpg%cvFt(omp(icsepomp))
      do i=1,mpg%ftCvP(iFt,2)
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
        if(mod(mpg%cvReg(iCv),4).eq.3.or.
     .     mod(mpg%cvReg(iCv),4).eq.0) cycle
        if(geo%cvX(iCv).gt.rmp) then
          jmp=i
          rmp=geo%cvX(iCv)
        end if
      end do
      if(jmp.eq.0) then
        write(*,*) 'Unexpected error in b2paux: midplane not found'
        return
      end if
c
c*** and then calculate the profiles taking x=R(midplane)/rmp
c*** and y=(R(midplane)-rmp)/(Redge-rmp)
c
      open(lun,file='fpri.dat')
      rmp=(geo%cvX(omp(icsepomp-1))+geo%cvX(omp(icsepomp)))/2.0_R8
      do j=1,icsepomp-1
        x=(geo%cvX(omp(j))-geo%cvX(omp(1)))/(rmp-geo%cvX(omp(1)))
        qpr_tmp(j)=rd_in(x)
        write(lun,*) j,x,qpr_tmp(j)
      end do
      close(lun)
      open(lun,file='fpro.dat')
      write(lun,*) icsepomp,0.0_R8,a+b
      do j=icsepomp,nomp
        x=(geo%cvX(omp(j))-rmp)/(geo%cvX(omp(nomp))-rmp)
        qpr_tmp(j)=rd_out(x)
        write(lun,*) j,x,qpr_tmp(j)
      end do
      close(lun)
c
c*** Calculate the poloidal profile qpz_tmp
c*** First, find the "top" location ...
      ztp=0.0_R8
      jtp=0
      do i=1,mpg%ftCvP(iFt,2)
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
        if(mod(mpg%cvReg(iCv),4).eq.3.or.
     .     mod(mpg%cvReg(iCv),4).eq.0) cycle
        if(geo%cvY(iCv).gt.ztp) then
          jtp=i
          ztp=geo%cvY(iCv)
        end if
      end do
      if(jtp.eq.0) then
        write(*,*) 'Unexpected error in b2paux: top not found'
        return
      end if
c
c*** ... calculate the poloidal length for normalisation ...
c
      ztt=0.0_R8
      do i=1,mpg%ftCvP(iFt,2)
        iCv=mpg%ftCv(mpg%ftCvP(iFt,1)+i-1)
        if(mod(mpg%cvReg(iCv),4).eq.3.or.
     .     mod(mpg%cvReg(iCv),4).eq.0) cycle
        ztt=ztt+geo%cvHx(iCv)
      end do
      sgb=sgb/ztt
c
c*** ... and then the poloidal profiles
c
#ifdef WG_TODO
      qpz_tmp(jtp)=1.0_R8
      j=jsep+1
      ztt=geo%hx(jtp,j)*0.5_R8
      do i=jtp+1,nx-1
        x=hx(i,j)*0.5_R8
        ztt=ztt+x
        qpz_tmp(i)=exp(-(ztt*sgb)**2)
        ztt=ztt+x
      end do
      ztt=hx(jtp,j)*0.5_R8
      do i=jtp-1,0,-1
        x=hx(i,j)*0.5_R8
        ztt=ztt+x
        qpz_tmp(i)=exp(-(ztt*sgb)**2)
        ztt=ztt+x
      end do
#endif
      open(lun,file='fprp.dat')
#ifdef WG_TODO
      do i=0,nx-1
        write(lun,*) i,qpz_tmp(i)
      end do
#endif
      close(lun)
c
c*** Now calculate the emissivity. Sample points along the sub-grid
c*** and calculate the average emissivity for each cell
c
      dz=(zmx-zmn)/(nnz*nqch)
      dr=(rmx-rmn)/(nnr*nqch)
      z=zmn-0.5_R8*dz
      do k=1,nnz
        do j=1,nqch
          z=z+dz
          r=rmn-0.5_R8*dr
          do i=1,nnr
            q_tmp(i,k)=0.0_R8
            do l=1,nqch
              r=r+dr
              call chord_find(mpg,geo,r,0.0_R8,z,found,iCv)
#ifdef WG_TODO
              if(found) then
                q_tmp(i,k)=q_tmp(i,k)+qpr_tmp(iy)*qpz_tmp(ix)
              end if
#endif
            end do
          end do
        end do
      end do
      bq=0.0_R8
      dr=(rmx-rmn)/nnr
      dz=(zmx-zmn)/nnz
      r=rmn+0.5_R8*dr
      do i=1,nnr
        do k=1,nnz
          bq=bq+q_tmp(i,k)*r
        end do
        r=r+dr
      end do
      bq=bq*dz*dr*2.0_R8*pi
      print *,'Q = ',bq
      qmn= huge(1.0_R8)
      qmx=-huge(1.0_R8)
      do k=1,nnz
        do i=1,nnr
          q_tmp(i,k)=q_tmp(i,k)*pwr_tot*1.e-3_R8/bq
          qmn=min(qmn,q_tmp(i,k))
          qmx=max(qmx,q_tmp(i,k))
        end do
      end do
c
      open(lun,file='plenter.dat')
      write(lun,'(1p,8e15.6)') ((q_tmp(i,k),i=1,nnr),k=1,nnz)
      close(lun)
c
      open(lun,file='plenter.grd')
      write(lun,'(a4/2i6,3(/2g16.5))')
     -                           'DSAA',nnr,nnz,rmn,rmx,zmn,zmx,qmn,qmx
      do k=1,nnz
        write(lun,'(1p,8e15.6)') (q_tmp(i,k),i=1,nnr)
      end do
      close(lun)
      deallocate(qpr_tmp,qpz_tmp)
      return

      contains

      function rd_in(x)
      implicit none
      real(kind=R8) :: rd_in, x
      rd_in = a+b*x**cc+d*(1.0_R8-x**e)**f
      return
      end function rd_in

      function rd_out(x)
      implicit none
      real(kind=R8) :: rd_out, x
      rd_out = a0+x*(a1+x*a2)
      return
      end function rd_out

      end subroutine b2pplentd
