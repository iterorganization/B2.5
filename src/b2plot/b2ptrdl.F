      subroutine b2ptrdl(nwe,ntg,rdl,tg_only)
c
c  version : 07.11.97 13:22
c
c======================================================================
c*** Calculates the radiation power load onto the targets
c*** rdl(i,1,:) corresponds to the West and rdl(i,2,:) to the East target
c*** Last index 0 is a sum over all species
c*** Power is in W/m**2 (previously, MW/m**2)
c======================================================================
      use b2us_data
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
      use b2mod_dimensions
#ifdef B25_EIRENE
      use eirmod_eirbra
#else
      use b2mod_eirbra
#endif
      implicit none
      integer, intent(in) :: nwe, ntg
      real (kind=R8) ::
     *   rdl(nwe,ntg,0:DEF_NATM)
      logical, intent(in) :: tg_only
      integer i,ii,j,jj,k,m,n,nRFc
#ifdef DBG
      integer istart,iend
#endif
c======================================================================
c*** If the load is already calculated, just fill in the output array
c
      if(.not.l2rdldt) then !{
c----------------------------------------------------------------------
#ifdef DBG
        print *,'b2ptrdl starting...'
#endif
c*** Input the DG data and prepare the geometry
c
        call b2pwinp
c<<<
#ifdef DBG
        print *,' !!!!!!!!!'
        print '(a/5a10/1p,l10,3i10,e10.2)','After b2pwinp',
     ,   'l2pwinp','lppp','nwldprt','nadsrc','add_rad',
     ,    l2pwinp , lppp , nwldprt , nadsrc , add_rad
        if (nwldlis(1).gt.0) then !{
          print '(a/(20i5))','Shadowing structure',
     ,                                   (lwldlis(i,1),i=1,nwldlis(1))
          iend=lppp
          do j=1,ntrgts !{
            if (itrg(j).eq.0) cycle
            istart=iend+1
            iend=iend+mpg%divFcP(itrg(j),2)
            print '(a,i2/(5(2f8.3,1x)))','Target',j,
     ,       (pp1(1,i),pp1(2,i),i=istart,iend),
     ,        pp2(1,iend),pp2(2,iend)
          end do !}
        else !}{
          print *,' No shadowing structure found!'
        end if !}
#else
        if (nwldlis(1).eq.0) write(*,*) 'No shadowing structure found!'
#endif
c>>>
        if(.not.l2pwinp) then !{
          write(*,*) '**** b2ptrdl: error detected in the input data ',
     ,                            '==> no radiation load calculated!!!'
          do k=0,DEF_NATM !{
            do j=1,ntg !{
              do i=1,maxval(mpg%divFcP(:,2)) !{
                rdl(i,j,k)=0.0_R8
              end do !}
            end do !}
          end do !}
          return
        end if !}
c----------------------------------------------------------------------
c*** Calculate the load

        call b2pwrdld

        if(lrdlsp.and.lrdlcr) then !{
          print *,'#### b2ptrdl: both lrdlsp and lrdlcr are true - ',
     ,                                                'internal error'
          lrdlsp=.false.
          lrdlcr=.false.
        end if !}
#ifdef DBG
        print *,'after b2pwrdld: lrdlsp=',lrdlsp !###
#endif
        l2rdldt=.true.
#ifdef DBG
        iend=lppp
        if (tg_only) then !{
          do j=1,ntrgts !{
            if (itrg(j).eq.0) cycle
            istart=iend+1
            iend=iend+mpg%divFcP(itrg(j),2)
            print '(a,i2,a/(1p,5e17.4))','Target',j,': load',
     ,       (wldrdld(i),i=istart,iend)
          end do !}
        else !}{
          do j=ncns+1,nsts !{
            if (eirdiag_nds_ind(j)+1.eq.eirdiag_nds_ind(j+1)) cycle
            istart=iend+1
            iend=iend+eirdiag_nds_ind(j+1)-(eirdiag_nds_ind(j)+1)
            print '(a,i3,a/(1p,5e17.4))',
     ,       'Eirene non-standard surface',j,': load',
     ,       (wldrdld(i),i=istart,iend)
          end do !}
        end if !}
#endif
      end if !}
c----------------------------------------------------------------------
c*** Arrange the output data

#ifdef B25_EIRENE
      n=natmi
#else
      n=nspecies
#endif
      k=0
      if (tg_only) then !{
        do j=1,ntrgts !{
          if (itrg(j).eq.0) cycle
          do i=1,mpg%divFcP(itrg(j),2) !{
            k=k+1
            rdl(i,itrg(j),0)=wldrdld(lwldrad(k))
            if(lrdlsp) then !{
              do m=1,n !{
                rdl(i,itrg(j),m)=wldrdlsp(lwldrad(k),m)
              end do !}
            end if !}
            if(lrdlcr) then !{
              rdl(i,itrg(j),1)=wldrdlsp(lwldrad(k),1)
            end if !}
          end do !}
        end do !}
      else if (switch%use_eirene.ne.0) then !}{
        jj=0
        do j=ncns+1,nsts !{
          if (eirdiag_nds_ind(j)+1.eq.eirdiag_nds_ind(j+1)) cycle
          jj=jj+1
          do i=eirdiag_nds_ind(j)+1,eirdiag_nds_ind(j+1)-1
            k=k+1
            if (lwldrad(k).eq.0) cycle
            ii=i-eirdiag_nds_ind(j)
            rdl(ii,jj,0)=wldrdld(lwldrad(k))
            if(lrdlsp) then !{
              do m=1,n !{
                rdl(ii,jj,m)=wldrdlsp(lwldrad(k),m)
              end do !}
            end if !}
            if(lrdlcr) then !{
              rdl(ii,jj,1)=wldrdlsp(lwldrad(k),1)
            end if !}
          end do
        end do !}
      else if (grid_mode.eq.2) then !}{
        do j=1,nstsi !{
          nRFc=count(mpg%fcLbl==j)
          do i=1,nRFc !{
            k=k+1
            if (lwldrad(k).eq.0) cycle
            rdl(i,j,0)=wldrdld(lwldrad(k))
            if(lrdlsp) then !{
              do m=1,n !{
                rdl(i,j,m)=wldrdlsp(lwldrad(k),m)
              end do !}
            end if !}
            if(lrdlcr) then !{
              rdl(i,j,1)=wldrdlsp(lwldrad(k),1)
            end if !}
          end do !}
        end do !}
      end if !}
c======================================================================
      return
      end subroutine b2ptrdl
