      subroutine b2ptrdl(rdl)
c
c  version : 07.11.97 13:22
c
c======================================================================
c*** Calculates the radiation power load onto the targets
c*** rdl(i,1,:) corresponds to the West and rdl(i,2,:) to the East target
c*** Last index 0 is a sum over all species
c*** Power is in W/m**2 (previously, MW/m**2)
c======================================================================
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
      implicit none
#     include <DIMENSIONS.F>
      real (kind=R8) ::
     *   rdl(0:ny-1,ntrgsx,0:DEF_NATM)
      integer i,j,k,l,m,n
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
c        print '(a/7a10/1p,l10,5i10,e10.2)','After b2pwinp','l2pwinp',
c     ,             'lppp','nwldprt','lytg','nytg','nadsrc','add_rad',
c     ,                   l2pwinp,lppp,nwldprt,lytg,nytg,nadsrc,add_rad
c        print '(a/(16i5))','Shadowing structure',
c     ,                                   (lwldlis(i,1),i=1,nwldlis(1))
c        print '(a/(5(2f8.3,1x)))','West target',
c     ,                  (pp1(1,lppp+i),pp1(2,lppp+i),i=1,nytg)
c        print '(a/(5(2f8.3,1x)))','East target',
c     ,                (pp1(1,lppp+nytg+i),pp1(2,lppp+nytg+i),i=1,nytg)
c
c>>>
        if(.not.l2pwinp) then !{
          write(*,*) '**** b2ptrdl: error detected in the input data ',
     ,                            '==> no radiation load calculated!!!'
          do j=1,ntrgsx !{
            do i=0,ny-1 !{
              do k=0,DEF_NATM !{
                rdl(i,j,k)=0.0_R8
              end do !}
            end do !}
          end do !}
          return
        end if !}
c----------------------------------------------------------------------
c*** Fill in the input arrays for b2pwrdld
c
        nwldrad=0
        do j=1,ntrgts !{
          do i=1,nytg !{
            lwldrad(nwldrad+i)=lppp+nwldrad+i
          end do !}
          nwldrad=nwldrad+nytg
        end do !}
#ifdef DBG
        print *,'%% b2ptrdl: nwldrad,nytg / lwldrad',nwldrad,nytg
        k=0
        do j=1,ntrgts !{
          print '(30i4)',(lwldrad(k+i),i=1,nytg)
          k=k+nytg
        end do !}
        print '(1x,a,i5/(30i4))','nwldlis / lwldlis',nwldlis(1),
     ,                             (lwldlis(i,1),i=1,nwldlis(1))
#endif

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
        print '(a/(1p,5e17.4))','West target: load',
     ,    (wldrdld(lppp+i+1),i=0,ny-1)
        print '(a/(1p,5e17.4))','East target: load',
     ,    (wldrdld(lppp+ny+i+1),i=0,ny-1)
#endif
      end if !}
c----------------------------------------------------------------------
c*** Arrange the output data

#ifdef B25_EIRENE
      n=natmi
#else
      n=nspecies
#endif
      k=nytg
      l=1
      do j=1,max(1,ntrgts/2) !{
        do i=lytg,ny-1 !{
          rdl(i,l,0)=wldrdld(lppp+k-i+lytg)
          rdl(i,l+1,0)=wldrdld(lppp+k+1+i-lytg)
          if(lrdlsp) then !{
            do m=1,n !{
              rdl(i,l,m)=wldrdlsp(lppp+k-i+lytg,m)
              rdl(i,l+1,m)=wldrdlsp(lppp+k+1+i-lytg,m)
            end do !}
          end if !}
          if(lrdlcr) then !{
            rdl(i,l,1)=wldrdlsp(lppp+k+1-i+lytg-1,1)
            rdl(i,l+1,1)=wldrdlsp(lppp+k+i-lytg+1,1)
          end if !}
        end do !}
        k=3*nytg
        l=3
      end do !}
c======================================================================
      return
      end subroutine b2ptrdl
