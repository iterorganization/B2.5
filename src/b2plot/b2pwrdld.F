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
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_rdload
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
      use b2us_data
      use b2mod_dimensions
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer iatm,imol,iion
      integer i,j,k,l,m,is,iCv,n_hlperr,
     ,        l_errt(nsgmx),l_errs(nsgmx)
      real (kind=R8) ::
     *   u
#ifdef DBG
      integer lun_w
      parameter (lun_w=96)
      character tab
#endif
c======================================================================
#ifdef DBG
      tab=char(9)
#endif
c*** Fill in the input arrays for rlcomp.
c*** First, translate the information on the wall segments belonging
c*** to the "shadowing structure" and "plot zones" into the format used
c*** in rlcomp. That is, we make a "structure" polygon, which is closed,
c*** and the target list referring to the segments of this polygon.
c
      lrdlsp=.false.
      lrdlcr=.false.
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

#ifdef DBG
      print *,'b2pwrdld: nplgs,nbndry,nwldrad=',nplgs,nbndry,nwldrad
      print *,'lwldlis(i,1):'
      print '(20i5)',(lwldlis(i,1),i=1,nplgs)
      print *,'lbnd(i):'
      print '(20i5)',(lbnd(i),i=1,nplgs)
      print *,'luwrk:'
      print '(20i5)',(luwrk(i),i=1,nsgmx)
      print *,'lewrk:'
      print '(20i5)',(lewrk(i),i=1,nsgmx)
#endif

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
#ifdef DBG
      print *,'Shadowing structure: ',
     ,         nplgs,nbndry,nwldlis(1),lwldlis(1,1)
      print '(20i5)',(lbnd(i),i=1,nplgs)
      print '(5(1x,2f9.4))',(xp(i),yp(i),i=1,nplgs)
      print '(a,12i5)','Radiation target ',ntrg,
     ,       (nwldlis(i),i=2,nwldprt)
      print '(20i5)',(ltrg(i),i=1,ntrg)

      open(lun_w,file='layout-w.dat')
      write(lun_w,'(2a)') '# Wall layout : ',trim(label)
      write(lun_w,'(3a)') '# R',tab,'Z'
      do i=1,nplgs
        write(lun_w,'(1p,2(1x,e15.8))') xp(i),yp(i)
      end do
      close(lun_w)
#endif
c
c*** Now specify the radiation sources - in W (remember, the additional
c*** sources are in MW !!!)
c
      do iCv=1,mpg%nCi
        u = 0.0_R8
        do is=0,state%ns-1
          u=u+state%srw%rqrad(iCv,is)+state%srw%rqbrm(iCv,is)
        end do
        do is=1,natmi
          u=u-eneutrad(iCv,is,0)
        end do
        do is=1,nmoli
          u=u-emolrad(iCv,is,0)
        end do
        do is=1,nioni
          u=u-eionrad(iCv,is,0)
        end do
        ssrcm(iCv)=u
        xysrc(1,iCv)=geo%cvX(iCv)
        xysrc(2,iCv)=geo%cvY(iCv)
      end do

#ifdef DBG
      print *,'nadsrc,add_rad = ',nadsrc,add_rad
#endif
      k=mpg%nCi
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
#ifdef DBG
      u=0.0_R8
      do i=1,nsrc
        u=u+ssrcm(i)
      end do
      print *,'Sources: nsrc, power = ',nsrc,u
#endif
c
c*** Calculate the radiation load
c
      call rlcomp
c      write(*,*) 'b2pwrdld: after rlcomp'
      write(*,*)   '  Total radiation         : ',ttsrc
      if (add_rad.gt.0.0_R8)
     >  write(*,*) '(of which add. core rad.) : ',add_rad*1.0e6
      if (irdparm.eq.1) then
        write(*,*) '  Reaching the targets    : ',tttrg
      else if (irdparm.eq.2 .and. nwldprt.gt.1) then
        write(*,*) '  Reaching the plot zones : ',tttrg
      end if
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
c*** After rlcomp, the radiation loads are stored in radld (b2mod_rdload)
c*** in the first ntrg elements, in the order given by the ltrg array.
c*** We need them put in wldrdld (b2mod_b2plot_wall_loading.F) in the
c*** positions just after lppp, in the following order:
c***  Western target, cells ny-1 to lytg;
c*** [ inner top target, cells lytg to ny-1;
c***   outer top target, cells ny-1 to lytg ];
c***  Eastern target, cells lytg to ny-1
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
        do j=1,nspecies !{
          wldrdlsp(i,j)=0.0_R8
        end do !}
      end do !}

      if(ildparm.eq.99) then !{

c*** Calculate the radiation load split by impurity species from B2

        do j=1,nspecies !{
          write(*,*) '  Radiation split: ',species_list(j)
          do i=1,nplgs !{
            radld(i)=0.0_R8
          end do !}
          do iCv=1,mpg%nCi !{
            ssrcm(iCv)=0.0_R8
            do iatm=1,natmi !{
              if (latmscl(iatm).eq.j)
     >         ssrcm(iCv)=ssrcm(iCv)-eneutrad(iCv,iatm,0)
            end do !}
            do imol=1,nmoli !{
              iatm = lmolscl(imol)
              if (iatm.eq.0) cycle
              if (latmscl(iatm).eq.j)
     >         ssrcm(iCv)=ssrcm(iCv)-emolrad(iCv,imol,0)
            end do !}
            do iion=1,nioni !{
              iatm = lionscl(iion)
              if (iatm.eq.0) cycle
              if (latmscl(iatm).eq.j)
     >         ssrcm(iCv)=ssrcm(iCv)-eionrad(iCv,iion,0)
            end do !}
            do is=eb2spcr(j),eb2spcr(j)+nfluids(j)
              ssrcm(iCv)=ssrcm(iCv)+state%srw%rqrad(iCv,is)+
     .                              state%srw%rqbrm(iCv,is)
            end do
            xysrc(1,iCv)=geo%cvX(iCv)
            xysrc(2,iCv)=geo%cvY(iCv)
          end do !}
          nsrc=mpg%nCi
          call rlcomp
          write(*,*)   '  Total radiation         : ',ttsrc
          if (irdparm.eq.1) then
            write(*,*) '  Reaching the targets    : ',tttrg
          else if (irdparm.eq.2 .and. nwldprt.gt.1) then
            write(*,*) '  Reaching the plot zones : ',tttrg
          end if
          k=0
          do i=1,nwldrad !{
            l=lwldrad(i)
            if(l.gt.0) then !{
              k=k+1
              wldrdlsp(l,j)=radld(k)/segm_area(l)
            end if !}
          end do !}
        end do !}
        lrdlsp=.true.
      end if !}

#ifdef DBG
      print *,'ildparm,nadsrc,add_rad: ',ildparm,nadsrc,add_rad
#endif
      if(ildparm.eq.98 .and. nadsrc.gt.0 .and. add_rad.gt.0.0_R8) then !{

c*** Calculate the radiation load from the core sources

        write(*,*) '  Extra radiation from the core       : ',
     .   1.0e6_R8*add_rad
        do i=1,nplgs !{
          radld(i)=0.0_R8
        end do !}
        k=0
        u=1.0e6_R8*add_rad/nadsrc
        do i=1,nadsrc !{
          k=k+1
          ssrcm(k)=u
          xysrc(1,k)=addsrc_xy(1,i)
          xysrc(2,k)=addsrc_xy(2,i)
        end do !}
        nsrc=k
        call rlcomp
#ifdef DBG
        write(*,*) '  Total additional core radiation     : ',ttsrc
#endif
        write(*,*) '  Core radiation reaching the targets : ',tttrg
        k=0
        is=1
        do i=1,nwldrad !{
          l=lwldrad(i)
          if(l.gt.0) then !{
            k=k+1
            wldrdlsp(l,is)=radld(k)/segm_area(l)
          end if !}
        end do !}
        lrdlcr=.true.
      end if !}
      return
c======================================================================
      end subroutine b2pwrdld
