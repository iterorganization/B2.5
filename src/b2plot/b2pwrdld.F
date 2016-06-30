      subroutine b2pwrdld
c
c  version : 18.02.98 18:24
c
c======================================================================
c*** Calculates radiation load onto the specified wall sections
c*** The list of the wall segments is to be specified in lwldrad
c======================================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_sources
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
      implicit none
#     include "DIMENSIONS.F"
#     include "rdload.inc"
      integer i,j,k,l,m,is,ix,iy,n_hlperr,iatm,
     ,        l_errt(nsgmx),l_errs(nsgmx)
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
      lrdlsp=.false.
      do i=1,nsgmx !{
        luwrk(i)=0
      end do !}
      n_hlperr=0
      nplgs=nwldlis(1)
      nbndry=nplgs
      do i=1,nplgs
        l=iabs(lwldlis(i,1))
        luwrk(l)=i
        xp(i)=pp1(1,l)
        yp(i)=pp1(2,l)
        lbnd(i)=i
        radld(i)=0.0_R8
      end do
	  
c      print *,'b2pwrdld: nplgs,nbndry,nwldrad=',nplgs,nbndry,nwldrad
c      print *,'lwldlis(i,1):'
c      print '(20i5)',(lwldlis(i,1),i=1,nplgs)
c      print *,'lbnd(i):'
c      print '(20i5)',(lbnd(i),i=1,nplgs)
c      print *,'luwrk:'
c      print '(20i5)',luwrk

      k=0
      do i=1,nwldrad
        l=lwldrad(i)
        m=luwrk(l)
        if(m .eq. 0) then
          write(*,'(2a,i4,a,i4)') 
     ,             'b2pwrdld: radiation target must belong to ',
     ,             'the shadowing structure. Target ',i,',  segment ',l
          lwldrad(i)=0
          n_hlperr=n_hlperr+1
          l_errt(n_hlperr)=i
          l_errs(n_hlperr)=l
        else
          k=k+1
          ltrg(k)=m
        end if
      end do
      ntrg=k
c<<<
c      print *,'Shadowing structure: ',nplgs,nbndry,nwldlis(1),lwldlis(1,1)
c      print '(5(1x,2f7.2))',(xp(i),yp(i),i=1,nplgs)
c      print *,'Radiation target ',ntrg,(nwldlis(i),i=2,nwldprt)
c      print '(16i5)',(ltrg(i),i=1,ntrg)

c      open(lun_w,file='layout-w.dat')
c      write(lun_w,'(6a)') '#Wall layout : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,nplgs
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
		  do is=1,natmi
			u=u-eneutrad(ix+1,iy+1,is,0)
		  end do
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=cr(ix,iy)
          xysrc(2,k)=cz(ix,iy)
        end do
      end do
c	  print*,ssrcm(1),(rqbrm(0,0,is),rqrad(0,0,is),is=0,ns-1),
c     .  (eneutrad(1,1,is,0),is=1,natmi)

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
c      write(*,*) '  Total radiation      : ',ttsrc
c      write(*,*) '  Reaching the targets : ',tttrg
      if(n_hlperr.gt.0) then !{
        write(*,*) '*** ',n_hlperr,
     ,                      ' segments outside the shadowing structure'
        write(*,*) '   indices in the segment list:'
        write(*,'(1x,20i5)') (l_errs(i),i=1,n_hlperr)
        write(*,*) '   indices in the target list:'
        write(*,'(1x,20i5)') (l_errt(i),i=1,n_hlperr)
        write(*,*)
      end if !}
c----------------------------------------------------------------------
c*** Extract the data and put them in the output module.
c*** After rlcomp, the radiation loads are stored in radld (rdload.inc)
c*** in the first ntrg elements, in the order given by the ltrg array.
c*** We need them put in wldrdld (b2mod_b2plot_wall_loading.F) in the
c*** positions just after lppp, in the following order:
c***  inner target, cells ny-1 to lytg;
c*** [ inner top target, cells lytg to ny-1;
c***   outer top target, cells ny-1 to lytg ];
c***  outer target, cells lytg to ny-1
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

      do i=1,nsgmx !{
        do j=1,natm !{
          wldrdlsp(i,j)=0.0_R8
        end do !}
      end do !}

      if(ildparm.eq.99) then !{

c*** Calculate the radiation load split by impurity species from B2

        do iatm=1,natmi !{
c          write(*,*) '  Radiation split: ',textan(iatm-1)
          do i=1,nplgs !{
            radld(i)=0.0_R8
          end do !}
          k=0
          do ix=0,nx-1 !{
            do iy=0,ny-1 !{
              k=k+1
              ssrcm(k)=-eneutrad(ix+1,iy+1,iatm,0)
              do is=eb2spcr(iatm),eb2spcr(iatm)+nfluids(iatm)
                ssrcm(k)=ssrcm(k)+rqrad(ix,iy,is)+rqbrm(ix,iy,is)
              end do
              xysrc(1,k)=cr(ix,iy)
              xysrc(2,k)=cz(ix,iy)
            end do !}
          end do !}
          nsrc=k
          call rlcomp
c		  write(*,*) 'b2pwrdld: after rlcomp'
c          write(*,*) '  Total radiation      : ',ttsrc
c          write(*,*) '  Reaching the targets : ',tttrg
          k=0
          do i=1,nwldrad !{
            l=lwldrad(i)
            if(l.gt.0) then !{
              k=k+1
              wldrdlsp(l,iatm)=radld(k)/segm_area(l)
            end if !}
          end do !}
        end do !}
        lrdlsp=.true.
      end if !}

      return
c======================================================================
      end subroutine b2pwrdld
