      subroutine b2pwrdld
c
c  version : 18.02.98 18:24
c
c======================================================================
c*** Calculates radiation load onto the specified wall sections
c*** The list of the wall segments is to be specified in lwldrad
c======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_sources
      implicit none
#     include "DIMENSIONS.F"
#     include "b2pwcom.h"
#     include "rdload.inc"
      integer i,k,l,m,is,ix,iy,luwrki(nsgmx)
      real (kind=R8) ::
     *   u
c<<<
c      integer lun_w
c      parameter (lun_w=76)
c      character tab
c>>>
c======================================================================
c<<<
c      tab=char(9)
c>>>
c*** Fill in the input arrays for rlcomp.
c*** First, translate the information on the wall segments belonging
c*** to the "shadowing structure" and "plot zones" into the format used
c*** in rlcomp. That is, we make a "structure" polygon, which is closed,
c*** and the target list referring to the segments of this polygon.
c
      do i=1,nsgmx
        luwrki(i)=0
      end do
      nplg=nwldlis(1)
      nbndry=nplg
      do i=1,nplg
        l=iabs(lwldlis(i,1))
        luwrk(l)=i
        xp(i)=pp1(1,l)
        yp(i)=pp1(2,l)
        lbnd(i)=i
        radld(i)=0.0_R8
      end do
c
      k=0
      do i=1,nwldrad
        l=lwldrad(i)
        m=luwrk(l)
        if(m .eq. 0) then
          write(*,*) 'b2pwrdld: radiation target must belong to ',
     ,             'the shadowing structure. Target ',i,',  segment ',l
          lwldrad(i)=0
        else
          k=k+1
          ltrg(k)=m
        end if
      end do
      ntrg=k
c<<<
c      print *,'Shadowing structure: ',nplg,nbndry,nwldlis(1),lwldlis(1,1)
c      print '(5(1x,2f7.2))',(xp(i),yp(i),i=1,nplg)
c      print *,'Radiation target ',ntrg,(nwldlis(i),i=2,nwldprt)
c      print '(16i5)',(ltrg(i),i=1,ntrg)

c      open(lun_w,file='layout-w.dat')
c      write(lun_w,'(6a)') '#Wall layout : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,nplg
c        write(lun_w,*) xp(i),tab,yp(i)
c      end do
c      close(lun_w)
c>>>
c
c*** Now specify the radiation sources - in W (remember, the additional
c*** sources are in MW !!!)
c
      k=0

c      print *,'!!!! DUMMY b2pwrdld: no radiation from plasma !!!!'

      do ix=0,nx-1
        do iy=0,ny-1
          u = 0.0_R8
          do is=0,ns-1
            u=u+rqrad(ix,iy,is)+rqbrm(ix,iy,is)
          end do
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=cr(ix,iy)
          xysrc(2,k)=cz(ix,iy)
        end do
      end do

c<<<
c      print *,'nadsrc,add_rad = ',nadsrc,add_rad
c>>>
      if(nadsrc .gt. 0  .and.  add_rad .gt. 0.0_R8) then
        u=1.0e6_R8*add_rad/nadsrc
        do i=1,nadsrc
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=addsrc_xy(1,i)
          xysrc(2,k)=addsrc_xy(2,i)
        end do
      end if

      nsrc=k
c<<<
c      u=0.0_R8
c      do i=1,nsrc
c        u=u+ssrcm(i)
c      end do
c      print *,'Sources: nsrc, power = ',nsrc,u
c>>>
c
c*** Calculate the radiation load
c
      call rlcomp
c      write(*,*) 'b2pwrdld: after rlcomp'
      write(*,*) '  Total radiation      : ',ttsrc
      write(*,*) '  Reaching the targets : ',tttrg
c----------------------------------------------------------------------
c*** Extract the data and put them in the output common block.
c*** After rlcomp, the radiation loads are stored in radld (rdload.inc)
c*** in the first ntrg elements, in the order given by the ltrg array.
c*** We need them put in wldrdld (b2pwcom.h) in the positions just
c*** after lppp, in the following order:
c***  inner target, cells ny-1 to 0;
c***  outer target, cells 0 to ny-1
c*** Normalisation per unit area is also done here.
c
      k=0
      do i=1,nwldrad
        l=lwldrad(i)
        if(l .gt. 0) then
          k=k+1
          wldrdld(l)=radld(k)/segm_area(l)
        end if
      end do
      return
c======================================================================
      end subroutine b2pwrdld
