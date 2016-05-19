      subroutine b2ptrdl(rdl)
c
c  version : 07.11.97 13:22
c
c======================================================================
c*** Calculates the radiation power load onto the targets
c*** rdl(i,1) corresponds to the inner and rdl(i,2) to the outer target
c*** Power is in W/m**2 (previously, MW/m**2)
c======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
      implicit none
#     include <DIMENSIONS.F>
      real (kind=R8) ::
     *   rdl(0:ny-1,ntrgsx,0:DEF_NATM)
      integer i,j,k,l,m
c======================================================================
c*** If the load is already calculated, just fill in the output array
c
      if(.not.l2rdldt) then !{
c----------------------------------------------------------------------
c        print *,'b2ptrdl starting...'
c*** Input the DG data and prepare the geometry
c
        call b2pwinp
c<<<
c        print '(a/5a10/1p,l10,3i10,e10.2)','After b2pwinp',
c     ,          'l2pwinp','lppp','nwldprt','nadsrc','add_rad',
c     ,                              l2pwinp,lppp,nwldprt,nadsrc,add_rad
c        print '(a/(16i5))','Shadowing structure',
c     ,                                     (lwldlis(i,1),i=1,nwldlis(1))
c        print '(a/(5(2f8.3,1x)))','Inner target',
c     ,                             (pp1(1,lppp+i),pp1(2,lppp+i),i=0,ny-1)
c        print '(a/(5(2f8.3,1x)))','Outer target',
c     ,                       (pp1(1,lppp+ny+i),pp1(2,lppp+ny+i),i=0,ny-1)
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
          do i=0,nytg-1 !{
            lwldrad(nwldrad+i+1)=lppp+nwldrad+i
          end do !}
          nwldrad=nwldrad+nytg
        end do !}
c!###{
c        print *,'%% b2ptrdl: nwldrad,nytg / lwldrad',nwldrad,nytg
c        k=0
c        do j=1,ntrgts !{
c          print '(1x,30i4)',(lwldrad(k+i),i=1,nytg)
c          k=k+nytg
c        end do !}
c        print '(1x,a,i5/(1x,30i4))','nwldlis / lwldlis',nwldlis(1),
c     ,                                    (lwldlis(i,1),i=1,nwldlis(1))
c!###}

c*** Calculate the load

        call b2pwrdld
c        print *,'after b2pwrdld: lrdlsp=',lrdlsp !###
        l2rdldt=.true.
c<<<
c        print '(a/(1p,5e17.4))','Inner target: load',
c     ,    (wldrdld(lppp+i+1),i=0,ny-1)
c        print '(a/(1p,5e17.4))','Outer target: load',
c     ,    (wldrdld(lppp+ny+i+1),i=0,ny-1)
c>>>
      end if !}
c----------------------------------------------------------------------
c*** Arrange the output data

      k=nytg
      l=1
      do j=1,ntrgts/2 !{
        do i=lytg,ny-1 !{
          rdl(i,l,0)=wldrdld(lppp+k+1-i+lytg-1)
          rdl(i,l+1,0)=wldrdld(lppp+k+i-lytg+1)
          if(lrdlsp) then !{
            do m=1,natmi !{
              rdl(i,l,m)=wldrdlsp(lppp+k+1-i+lytg-1,m)
              rdl(i,l+1,m)=wldrdlsp(lppp+k+i-lytg+1,m)
            end do !}
          end if !}
        end do !}
        k=3*nytg
        l=3
      end do !}
c======================================================================
      return
      end subroutine b2ptrdl
