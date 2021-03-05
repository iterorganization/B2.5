      subroutine b2pwlprp
c
c  version : 18.02.98 18:08
c
c======================================================================
c*** Prepares the geometry for the wall loading plots
c***
c*** Input
c***  nwldprt : number of the "plot zones"
c***  lwldstrt: the initial segments of the chains (orientation)
c***  nwldlis : lengths of the lists of the "plot zones"
c***  lwldlis : lists of the "plot zones" as read from DG
c***  pp1,pp2 : arrays of the wall segments as read from DG
c***  p1, p2  : arrays of the wall segments as saved in Eirene
c***
c*** Output
c***  lwldlis : sorted lists of the "plot zones" (DG numbering)
c***  leirlis : sorted lists of the "plot zones" (Eirene numbering)
c***  lwlddir : orientation of the "plot zones" (true = CCW)
c***  pp1,pp2 : completed arrays of the wall segments
c======================================================================
      use b2mod_types
      use b2mod_b2plot
      use b2mod_indirect
      use b2mod_boundary_namelist
      use b2mod_neutrals_namelist
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      use b2mod_subsys
      implicit none
#include <DIMENSIONS.F>
      integer i,j,k,n,il
      logical lim_cnf,target_east,target_west
c<<<
      integer lun_w
      parameter (lun_w=76)
      character tab
c>>>
      integer use_eirene
      data use_eirene/0/
      save use_eirene
      external ipgeti, ipgetr, xerrab
c======================================================================
      call subini ('b2pwlprp')
      call get_boundary_namelist
      call ipgeti ('b2mndr_eirene', use_eirene)
c<<<
      tab=char(9)
c
c      print *,'Starting b2pwlprp...',nlim,lppp,DEF_NLIM
c      if(nnreg(0).lt.8) then
c        print '(a,4i6)','jedgi1, jedgi2, jedgo1, jedgo2 =',jedge(1:4)
c      else
c        print '(a,4i6)','jedgi1, jedgi2, jedgi3, jedgi4 =',jedge(1:4)
c        print '(a,4i6)','jedgo1, jedgo2, jedgo3, jedgo4 =',jedge(5:8)
c      end if
c      print '(a,8i6)','ltrgedge = ',ltrgedge
c      print '(a,8i6)','mtrgedge = ',mtrgedge
c
cc+++ layout-a.dat: wall as read from DG
c      open(lun_w,file='layout-a.dat')
c      write(lun_w,'(2a)') '# Wall layout (dg): ',trim(label)
c      write(lun_w,'(3a)') '# R',tab,'Z'
c      do i=1,lppp !{
c       write(lun_w,'(1p,2(1x,e15.8))') pp1(1,i),pp1(2,i)
c       write(lun_w,'(1p,2(1x,e15.8))') pp2(1,i),pp2(2,i)
c       write(lun_w,*)
c      end do !}
c      close(lun_w)
c>>>
      lim_cnf=nncut.eq.0 .and. nnreg(0).eq.2
      nytg=ny
      lytg=0
      if(lim_cnf) then !{
        lytg=topcut(1)+1
        nytg=ny-topcut(1)
      end if !}
c      print *,'lim_cnf,nytg,lytg=',lim_cnf,nytg,lytg !###

c*** Make the list of references to the Eirene "Additional surfaces" -
c*** count for "unused elements", "elements not for Eirene", and
c*** "numbering offset" in the DG model

      if(lvdg20) then !{
        k=1
        do j=1,lppp !{
          lelmind(j)=k+lnumboff
          do i=1,nskipel !{
            if(lskipel(i).eq.j) lelmind(j)=0
          end do !}
          if(lelmind(j).gt.0) k=k+1
        end do !}
      else !}{
        do j=1,lppp !{
          lelmind(j)=j
        end do !}
      end if !}
c      print *,'lelmind:'                                !###
c      print '(10(i6,2h :,i4))',(lelmind(i),i,i=1,lppp)  !###

c
c*** Check the presence of the shadowing structure
c
      n=nwldlis(1)
      if(n .le. 0) then
        write(*,*) 'b2pwlprp: no shadowing structure in the data!'
        call walloff
        l2pwinp=.false.
        call subend ()
        return
      end if
c
c*** Convert the input wall coordinates into double precision
c*** and scale them to metres.
c*** We have now two sets of the wall segments, one saved by Eirene
c*** and (another) one read from DG file. Only the segments which are
c*** identical in both sets are to be used for neutral loading, the rest
c*** will be filled with the value representing zero in the list-
c*** sorting routine b2pwsort. The reference distance along the wall
c*** will however be calculated from the geometry read from DG file.
c*** This should allow one to have the same reference frame for the
c*** radiation and neutral wall load.
c
      do i=1,lppp
        pp1(1,i)=0.001_R8*ppp1(1,i)
        pp1(2,i)=0.001_R8*ppp1(2,i)
        pp2(1,i)=0.001_R8*ppp2(1,i)
        pp2(2,i)=0.001_R8*ppp2(2,i)
      end do
c<<<
c      print *,'nwldprt =',nwldprt
c      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5,1x,a/(5x,15i5))','lwldlis :',il,pltzlbl(il),
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c      end do
c
c      open(lun_w,file='layout-d.dat')
c      write(lun_w,'(2a)') '# Wall layout (DG) : ',trim(label)
c      write(lun_w,'(3a)') '# R',tab,'Z'
c      do i=1,lppp
c        write(lun_w,'(1p,2(1x,e15.8))') pp1(1,i),pp1(2,i)
c        write(lun_w,'(1p,2(1x,e15.8))') pp2(1,i),pp2(2,i)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>

c*** Check the topology of the grid and the list of target edges

      if(ntrgsa.gt.0) then !{
        ntrgs=ntrgsa
        if (ntrgs.eq.1) then
          if (mtrgedge(2,1).gt.0) ltrgedge(1,1)=mtrgedge(2,1)
          if (mtrgedge(1,1).gt.0) ltrgedge(2,1)=mtrgedge(1,1)
        else
          do i=1,ntrgs !{
            if (mtrgedge(1,i).gt.0) ltrgedge(1,i)=mtrgedge(1,i)
            if (mtrgedge(2,i).gt.0) ltrgedge(2,i)=mtrgedge(2,i)
          end do !}
        end if
      end if !}
      ntrgts=2
      if(nncut.eq.0 .and. .not.lim_cnf) then
        if (periodic_bc.eq.1) then
          ntrgts=0
          itrg=0
        else if (boundary_namelist.ne.0) then
          ntrgts=0
          do i = 1, nbc
            if(bcchar(i).eq.'E'.and.bcpos(i).eq.-1) then
              target_east = bcene(i).eq.3 .or. bcene(i).eq.12 .or.
     &                      bcene(i).eq.15
            end if
            if(bcchar(i).eq.'W'.and.bcpos(i).eq.nx) then
               target_west = bcene(i).eq.3 .or. bcene(i).eq.12 .or.
     &                       bcene(i).eq.15
            end if
          end do
          if (target_west) then
            ntrgs=1
            itrg(1)=1
          end if
          if (target_east) then
            ntrgs=ntrgs+1
            itrg(ntrgs)=2
          end if
        else if (use_eirene.ne.0) then
          itrg=0
          if (count(ktns(:)==1).eq.0) ntrgts=ntrgts-1
          if (count(ktns(:)==2).eq.0) ntrgts=ntrgts-1
          if (count(ktns(:)==1).gt.0) then
            itrg(1)=1
            if (count(ktns(:)==2).gt.0) itrg(2)=2
          else if (count(ktns(:)==2).gt.0) then
            itrg(1)=2
          end if
        end if
      end if
      if(nncut.eq.1 .or. lim_cnf) then
        ntrgts=2
        itrg(1)=1
        itrg(2)=2
      end if
      if(nncut.eq.2) then
        ntrgts=4
        itrg(1)=1
        itrg(2)=2
        itrg(3)=3
        itrg(4)=4
      end if
      if(ntrgs.gt.ntrgts) then !{
        write(*,*) 'b2pwlprp: inconsistent topology. ',
     ,                                      'Check the "param.dg" file'
        call xerrab ('Inconsistent topology in b2pwlprp')
      end if !}
      if(ntrgs.eq.0) then !{
c*** use the obsolete jedg[io][12] parameters
        k=0
        do i=1,2*max(1,nncut) !{
          if (itrg(i).gt.0) then !{
            ltrgedge(2-k,itrg(i))=jedge(i+(1-k)*ntrgsx) ! jedgo
            ltrgedge(1+k,itrg(i))=jedge(i+k*ntrgsx)     ! jedgi
          end if !}
          k=1-k
        end do !}
      end if !}
c*** Check the target edge specifications and make them clear
      do i=1,ntrgts !{
        if(ltrgedge(1,i).lt.0) ltrgedge(1,i)=0
        if(ltrgedge(2,i).lt.0) ltrgedge(2,i)=0
      end do !}
c*** No PFR for limiter configuration
      if(lim_cnf) then !{
        ltrgedge(1,1)=0
        ltrgedge(1,2)=0
      end if !}

c*** Prepare the list of grid indices for the target flux faces

      if (itrg(1).eq.1) then
        ltrgpos(1)=-1
        lseppos(1)=topcut(1)-1
        lpfrpos(1)=leftcut(1)-1
      else if (itrg(1).eq.2) then
        ltrgpos(1)=nx-1
        lseppos(1)=jsep
        lpfrpos(1)=3*nx/4
      endif
      k=2
      if(ntrgts.eq.4) then
        ltrgpos(2)=nxtl-1
        lseppos(2)=topcut(2)-1
        lpfrpos(2)=leftcut(2)
        ltrgpos(3)=nxtr
        lseppos(3)=topcut(2)-1
        lpfrpos(3)=rightcut(2)-1
        k=4
      end if !}
      if (itrg(k).eq.ntrgts) then
        ltrgpos(k)=nx-1
        lseppos(k)=topcut(1)-1
        lpfrpos(k)=rightcut(1)
      else if (itrg(k).eq.1) then
        ltrgpos(k)=-1
        lseppos(k)=jsep
        lpfrpos(k)=nx/4
      endif
#ifdef DBG
      print '(a,8i6)','ltrgedge: ',((ltrgedge(i,j),i=1,2),j=1,ntrgts) !###
      print '(a,8i6)','mtrgedge: ',((mtrgedge(i,j),i=1,2),j=1,ntrgts) !###
      print '(a,4i6)','%% ltrgpos : ',(ltrgpos(i),i=1,ntrgts)         !###
      print '(a,4i6)','%% lseppos : ',(lseppos(i),i=1,ntrgts)         !###
      print '(a,4i6)','%% lpfrpos : ',(lpfrpos(i),i=1,ntrgts)         !###
#endif

c*** Sort the segments in the lists of the "wall zones"
c*** The shadowing structure must be closed, therefore it starts from
c*** an arbitrary point.

      lwldstrt(1)=iabs(lwldlis(1,1))
      call b2pwsort

c*** The sorting routine, b2pwsort, also checks the wall segments. If a
c*** segment input from DG does not coincide with the corresponding
c*** segment saved by Eirene, its number becomes negative, and no
c*** neutral-related data will be plotted for it.

c#ifdef DBG
c      print *,'After sorting. nwldprt =',nwldprt
c      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      print '(a,(20l5))','lwlddir :',(lwlddir(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5,1x,a/(5x,15i5))','lwldlis :',il,pltzlbl(il),
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c        print '(a,i5/(5x,15i5))','leirlis :',il,
c     ,                                  (leirlis(i,il),i=1,nwldlis(il))
c      end do
c#endif
      call b2pwltrg
c======================================================================
      call subend ()
      return
      end subroutine b2pwlprp

      subroutine b2pwltrg
      use b2mod_types
      use b2mod_geo
      use b2mod_diag
      use b2mod_b2cmfs
      use b2mod_b2plot
      use b2mod_constants
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer i,i1,i2,ii,il,is,ias,ind,iss,itmp,ix,iy,j,jj,
     *  k,k1,k2,kb,ke,kk,kn,kkb,kke,l,n,cont,flagmin,flagmin1,ksel
      integer ib(ntrgsx),ie(ntrgsx),it(ntrgsx),m_hlp(ntrgsx)
      integer :: nwldlis_old
      integer, allocatable :: lwldlis_old(:), leirlis_old(:)
      integer, allocatable :: lis_ind(:)
      real (kind=R8) :: xx1,xx2,yy1,yy2,distmin,distmin1,distance
      logical lim_cnf
      logical, allocatable :: istarget(:)
      data ib,ie,it/ntrgsx*0,ntrgsx*0,ntrgsx*0/
      real (kind=R8) ::
     *   d,x
c
      real (kind=R8) :: width
      data width /1.0_R8/
      save width
c*** wsrttol : a tolerance parameter for defining coinciding segments
      real (kind=R8) ::
     *   wsrttol
      parameter (wsrttol=1.0e-6_R8)
c<<<
      integer lun_w, lfn_hlp
      parameter (lun_w=76, lfn_hlp=mlblprtx+9)
      character tab, fn_hlp*(lfn_hlp)
c
c======================================================================
c*** Add segments for the B2 targets to the wall structures
c======================================================================
      call subini ('b2pwltrg')
      call ipgetr ('b2agmt_1d_width', width)
      lim_cnf=nncut.eq.0 .and. nnreg(0).eq.2
      tab=char(9)
c
c*** Check the shadowing structure - it must be closed
c
      n=nwldlis(1)
      l=iabs(lwldlis(1,1))
      k=iabs(lwldlis(n,1))
      if(lwlddir(1)) then
        if(abs(pp1(1,l)-pp2(1,k)) .gt. wsrttol  .or. 
     .     abs(pp1(2,l)-pp2(2,k)) .gt. wsrttol) then
          write(*,*) 'b2pwltrg: the shadowing structure must be closed!'
          write(*,*) 'lwlddir(1) = ',lwlddir(1)
          write(*,*) 'Starting segment: ',l,pp1(1,l),pp1(2,l)
          write(*,*) 'Ending segment  : ',k,pp2(1,k),pp2(2,k)
          call walloff
          l2pwinp=.false.
          call subend ()
          return
        end if
      else
        if(abs(pp1(1,k)-pp2(1,l)) .gt. wsrttol  .or. 
     .     abs(pp1(2,k)-pp2(2,l)) .gt. wsrttol) then
          write(*,*) 'b2pwltrg: the shadowing structure must be closed!'
          write(*,*) 'lwlddir(1) = ',lwlddir(1)
          write(*,*) 'Starting segment: ',l,pp2(1,l),pp2(2,l)
          write(*,*) 'Ending segment  : ',k,pp1(1,k),pp1(2,k)
          call walloff
          l2pwinp=.false.
          call subend ()
          return
        end if
      end if
c
c
c*** Add the target sections from B2 grid to the segment lists
c*** All CCW, targets in the B2 order - clockwise from left bottom
c
cc<<<
c      print *,'Grid edges'
c      print '(a4,8a12)','j','r1i','r2i','z1i','z2i',
c     ,                      'r1o','r2o','z1o','z2o'
cc>>>
      do i=1,nytg
        j=lytg+i
        k=nytg*(ntrgts-1)
        if (ntrgts.gt.1.or.ltrgpos(1).eq.-1) then
          pp1(1,i+lppp)=crx(0,ny-i,2)
          pp1(2,i+lppp)=cry(0,ny-i,2)
          pp2(1,i+lppp)=crx(0,ny-i,0)
          pp2(2,i+lppp)=cry(0,ny-i,0)
        end if
        if (ntrgts.gt.1.or.ltrgpos(1).eq.nx-1) then
          pp1(1,i+lppp+k)=crx(nx,j-1,0)
          pp1(2,i+lppp+k)=cry(nx,j-1,0)
          pp2(1,i+lppp+k)=crx(nx,j-1,2)
          pp2(2,i+lppp+k)=cry(nx,j-1,2)
        end if
      end do
      if(ntrgts.eq.4) then !{
        k=nxtl
        do i=1,nytg !{
          j=lytg+i
          pp1(1,i+lppp+nytg)=crx(k,j-1,0)
          pp1(2,i+lppp+nytg)=cry(k,j-1,0)
          pp2(1,i+lppp+nytg)=crx(k,j-1,2)
          pp2(2,i+lppp+nytg)=cry(k,j-1,2)
          pp1(1,i+lppp+2*nytg)=crx(nxtr,ny-i,3)
          pp1(2,i+lppp+2*nytg)=cry(nxtr,ny-i,3)
          pp2(1,i+lppp+2*nytg)=crx(nxtr,ny-i,1)
          pp2(2,i+lppp+2*nytg)=cry(nxtr,ny-i,1)
        end do !}
      end if !}
c<<<
c      print '(1x,i4,1p,8e12.4)',
c     , (i,(pp1(j,i+1+lppp)   ,pp2(j,i+1+lppp)   ,j=1,2),
c     ,    (pp1(j,i+1+lppp+ny),pp2(j,i+1+lppp+ny),j=1,2),i=0,ny-1)
c
c      open(lun_w,file='layout-g.dat')
c      write(lun_w,'(2a)') '# Wall layout (g) : ',trim(label)
c      write(lun_w,'(3a)') '# R',tab,'Z'
c      do i=1,lppp+2*ny
c        write(lun_w,'(1p,2(1x,e15.8))') pp1(1,i),pp1(2,i)
c        write(lun_w,'(1p,2(1x,e15.8))') pp2(1,i),pp2(2,i)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>
c
c*** Modify the target edges for continuity with B2 grid edges
c
c      if (nnreg(0).lt.8) then
c        print '(a,4i6)','Before modification: ', (ltrgedge(1:2,j),j=1,2)
c        print '((2(1x,a3,1p,4e12.4)))',
c     ,          'i1',(pp1(j,ltrgedge(1,1)),pp2(j,ltrgedge(1,1)),j=1,2),
c     ,          'o1',(pp1(j,ltrgedge(2,1)),pp2(j,ltrgedge(2,1)),j=1,2),
c     ,          'i2',(pp1(j,ltrgedge(1,2)),pp2(j,ltrgedge(1,2)),j=1,2),
c     ,          'o2',(pp1(j,ltrgedge(2,2)),pp2(j,ltrgedge(2,2)),j=1,2)
c      else
c        print '(a,8i6)','Before modification: ', (ltrgedge(1:2,j),j=1,4)
c        print '((2(1x,a3,1p,4e12.4)))',
c     ,          'i1',(pp1(j,ltrgedge(1,1)),pp2(j,ltrgedge(1,1)),j=1,2),
c     ,          'o1',(pp1(j,ltrgedge(2,1)),pp2(j,ltrgedge(2,1)),j=1,2),
c     ,          'i2',(pp1(j,ltrgedge(1,2)),pp2(j,ltrgedge(1,2)),j=1,2),
c     ,          'o2',(pp1(j,ltrgedge(2,2)),pp2(j,ltrgedge(2,2)),j=1,2),
c     ,          'i3',(pp1(j,ltrgedge(1,3)),pp2(j,ltrgedge(1,3)),j=1,2),
c     ,          'o3',(pp1(j,ltrgedge(2,3)),pp2(j,ltrgedge(2,3)),j=1,2),
c     ,          'i4',(pp1(j,ltrgedge(1,4)),pp2(j,ltrgedge(1,4)),j=1,2),
c     ,          'o4',(pp1(j,ltrgedge(2,4)),pp2(j,ltrgedge(2,4)),j=1,2)
c      end if

      do i=1,ntrgts !{
        if (mod(i,2).eq.1.and.ntrgts.ge.2) then
          m_hlp(i)=ltrgpos(i)+1
        else
          m_hlp(i)=ltrgpos(i)
        end if
      end do !}

      k=1
      iy=0
      do i=1,ntrgts !{
c        print *,'i,k,ltrgedges: ',i,k,ltrgedge(1,i),ltrgedge(2,i) !###
        if(ltrgedge(1,i).gt.0 .or. ltrgedge(2,i).gt.0) then !{
          ix=m_hlp(i)
          l=ltrgedge(k,i)
c          print *,'pp1: l,ix,iy ',l,ix,iy !###
          if(l.gt.0) then !{
            if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
              iy=ny-1-iy
              pp1(1,l)=crx(ix,iy,3)
              pp1(2,l)=cry(ix,iy,3)
            else
              pp1(1,l)=crx(ix,iy,0)
              pp1(2,l)=cry(ix,iy,0)
            end if
          end if !}
          k=3-k
          iy=ny-1-iy
          l=ltrgedge(k,i)
c          print *,'pp2: l,ix,iy ',l,ix,iy !###
          if(l.gt.0) then !{
            if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
              pp2(1,l)=crx(ix,iy,0)
              pp2(2,l)=cry(ix,iy,0)
            else
              pp2(1,l)=crx(ix,iy,3)
              pp2(2,l)=cry(ix,iy,3)
            end if
          end if !}
          k=3-k
          iy=ny-1-iy
        end if !}
        k=3-k
        iy=ny-1-iy
      end do !}

c<<<
c      open(lun_w,file='layout-h.dat')
c      write(lun_w,'(2a)') '# Wall layout (h) : ',trim(label)
c      write(lun_w,'(3a)') '# R',tab,'Z'
c      do i=1,lppp !{
c        write(lun_w,'(1p,2(1x,e15.8))') pp1(1,i),pp1(2,i)
c        write(lun_w,'(1p,2(1x,e15.8))') pp2(1,i),pp2(2,i)
c        write(lun_w,*)
c      end do !}
c      close(lun_w)
c+++ layout-t.dat: the grid edges - must ~ coincide with targets
c      open(lun_w,file='layout-t.dat')
c      write(lun_w,'(2a)') '# Target layout (h) : ',trim(label)
c      write(lun_w,'(3a)') '# R',tab,'Z'
c      k=lppp
c      do j=1,ntrgts !{
c       if(pp1(1,k+1).gt.0. .or. pp2(1,k+1).gt.0) then !{
c         do i=k+1,k+nytg !{
c           write(lun_w,'(1p,2(1x,e15.8))') pp1(1,i),pp1(2,i)
c           write(lun_w,'(1p,2(1x,e15.8))') pp2(1,i),pp2(2,i)
c           write(lun_w,*)
c         end do !}
c       end if !}
c       k=k+nytg
c      end do !}
c      close(lun_w)
c>>>

c*** Check whether the modified edges coincide with the saved ones

      do k=1,ntrgts !{
        do l=1,2 !{
          j=ltrgedge(l,k)
          if(j.gt.0) then !{
            jj=lelmind(j)
            if (jj.gt.0) then !{
c!###{
c              print *,'%% Check ',l,k,j,jj
c              print '(14x,4(a2,13x))','R1','Z1','R2','Z2'
c              print '(1x,a8,1p,4e15.6)','Eirene',
c     ,                              p1(1,jj),p1(2,jj),p2(1,jj),p2(2,jj)
c              print '(1x,a8,1p,4e15.6)','b2plot',
c     ,                              pp1(1,j),pp1(2,j),pp2(1,j),pp2(2,j)
c!###}
              if(abs(p1(1,jj)-pp1(1,j)) .le. wsrttol  .and.
     .           abs(p1(2,jj)-pp1(2,j)) .le. wsrttol  .and.
     .           abs(p2(1,jj)-pp2(1,j)) .le. wsrttol  .and.
     .           abs(p2(2,jj)-pp2(2,j)) .le. wsrttol) then !{
c                print *,'== OK'
                do il=1,nwldprt !{
                  n=nwldlis(il)
                  do i=1,n !{
                    if(iabs(lwldlis(i,il)).eq.j) lwldlis(i,il)=j
                  end do !}
                end do !}
              else if (p1(1,jj).ne.0.0_R8 .or. p1(2,jj).ne.0.0_R8 .or.
     .                 p2(1,jj).ne.0.0_R8 .or. p2(2,jj).ne.0.0_R8) then !}{
                print '(a,i3,1p,4e12.4)',
     .           'Eirene  : ',jj,p1(1,jj),p1(2,jj),p2(1,jj),p2(2,jj)
                print '(a,i3,1p,4e12.4)',
     .           'DG (mod): ', j,pp1(1,j),pp1(2,j),pp2(1,j),pp2(2,j)
              end if !}
            else !}{
              write(0,*) 'Edges check failed - not in Eirene:',l,k,j,jj
              call xerrab('==> Check the DG model')
            end if !}
          end if !}
        end do !}
      end do !}
c
c<<<
c      print *,'ltrgedge: ',((ltrgedge(i,j),i=1,2),j=1,ntrgts)
c      do j=1,ntrgts !{
c       i1=ltrgedge(1,j)
c       if(i1.gt.0) print '((2(1x,a3,i1,2x,1p,4e12.4)))',
c     ,                        'pfr',j,(pp1(i,i1),pp2(i,i1),i=1,2)
c       i2=ltrgedge(2,j)
c       if(i2.gt.0) print '((2(1x,a3,i1,2x,1p,4e12.4)))',
c     ,                        'sol',j,(pp1(i,i2),pp2(i,i2),i=1,2)
c      end do !}
cc+++ layout-e.dat:  modified target edges
c      open(lun_w,file='layout-e.dat')
c      write(lun_w,'(6a)') '#Target edges : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do j=1,ntrgts !{
c        do k=1,2 !{
c          i=ltrgedge(k,j)
c          if(i.gt.0) then !{
c            write(lun_w,*) pp1(1,i),tab,pp1(2,i)
c            write(lun_w,*) pp2(1,i),tab,pp2(2,i)
c            write(lun_w,*)
c          end if !}
c        end do !}
c      end do !}
c      close(lun_w)
c>>>
c
      allocate(istarget(nwldprt))
      istarget=.false.
c*** Correct the segment lists for the refined targets
c*** All the segments between the target edges - or between a target
c*** edge and the end of the list - are replaced with the segments
c*** of the target sections of the B2 grid.
c
c*** Below, kb marks the target; ib, ie are the locations of the inner
c*** and outer edges in the lists
c
      do il=1,nwldprt !{
        n=nwldlis(il)
        if(n .gt. 0) then !{
          do i=1,ntrgsx !{
            ib(i)=0
            ie(i)=0
            it(i)=0
          end do !}
          kb=0
          if(lwlddir(il)) then !{
            is=1
            ii=1
          else !}{
            is=n
            ii=-1
          end if !}
c          print *,'il,n,is,ii/lwldlis : ',il,n,is,ii    !###
c          print '(1x,30i4)',(lwldlis(i,il),i=1,n)       !###
          do i=1,n
            k=iabs(lwldlis(is,il))
            if(k.eq.ltrgedge(2,1) .or. k.eq.ltrgedge(1,2) .or.
     .          (ntrgts.eq.4 .and.
     .        (k.eq.ltrgedge(2,3) .or. k.eq.ltrgedge(1,4)))) then !{
c*** CW edge
c              print *,'CW edge: k=',k                   !###
              if(kb.ne.0) then !{
c*** Target still open...
                write(*,'(2a,3i4)') 'b2pwltrg: unexpected error. ',
     ,                                  'Missing CCW edge ',il,i,kb
                write(*,*) 'k,ltrgedge = ',
     ,                         k,((ltrgedge(i1,i2),i1=1,2),i2=1,ntrgts)
                write(*,*) 'lwldlis : ',pltzlbl(il)
                write(*,*) (lwldlis(i1,il),i1=1,n)
                call xerrab ('Missing CCW edge in b2pwltrg')
              end if !}
              i1=2
              do i2=1,ntrgts !{
                if(k.eq.ltrgedge(i1,i2)) kb=i2
                i1=3-i1
              end do !}
              ib(kb)=is
              it(kb)=kb
c              print *,'kb,ib(kb)=',kb,ib(kb)            !###
            else if(k.eq.ltrgedge(1,1) .or. k.eq.ltrgedge(2,2) .or.
     .        (ntrgts.eq.4 .and.
     .             (k.eq.ltrgedge(1,3) .or. k.eq.ltrgedge(2,4)))) then !}{
c*** CCW edge. Check to which target this edge belongs
              i1=1
              do i2=1,ntrgts !{
                if(k.eq.ltrgedge(i1,i2)) kkb=i2
                i1=3-i1
              end do !}
              if(kb.eq.0) kb=kkb
              if (kb.ne.kkb .and. .not. lim_cnf) then !{
                write(*,*) 'b2pwltrg: unexpected error. ',
     ,                                  'Inner edge not found ',il,i,kb
                write(*,*) 'k,ltrgedge = ',
     ,                         k,((ltrgedge(i1,i2),i1=1,2),i2=1,ntrgts)
                call xerrab ('Inner edge not found in b2pwltrg')
              end if !}
              ie(kkb)=is
              it(kkb)=kkb
              kb=0
            end if !}
            is=is+ii
          end do !}
c
c          if(kb .gt. 0) then
cc*** A target left open - consider it spanning through the edge
c            ie(kb)=is
c          end if
c
c*** Replace the target section with the B2 grid edge.
c*** First, check the sequence. We have here the indices ib and ie,
c*** marking the positions of the target edges in the lists.
c*** ib, ie mark the CW and CCW edges.
c***
c*** The list orientation is given by lwlddir value (true=CCW)

          kk=0
          do i=1,ntrgts !{
            if(ib(i)+ie(i).gt.0) then !{
              kk=kk+1
            end if !}
          end do !}
c<<<
c          print *,'After the check. il,kk,lppp = ',il,kk,lppp
c          print '(2a6)','ib','ie'
c          print '(2i6)',(ib(i),ie(i),i=1,ntrgts)
c>>>

c*** Now fill in the necessary arrays

          if(kk .gt. 0) then
            istarget(il)=.true.

c!###{
c           print *,         'Before merging the targets ',
c     ,          ':      ib       :      ie       :      it       :'
c           print '(a,24i4)','                            ',ib,ie,it
c           print *,'lytg,nytg= ',lytg,nytg
c!###}

c*** Go along the list analysing the intersections. Here
c***   kn      is the counter for the new list
c***   j       is the counter for the old list
c***   k       if nonzero, the target number k is "open"
c***   kb      if nonzero, the CW edge number kb is encountered
c***   ke      if nonzero, the CCW edge number ke is encountered
c***   kke     if nonzero, the first segment was CCW edge abs(kke)
c***           if negative, the corresponding target has been added
c***
c*** It is assumed here that CW and CCW edges differ
c*** Build only the shadowing structure with the additional wall segments
c*** from target cells
c*** Other plot zones will be defined using starting/ending segments and
c*** selecting the chain from shadowing structure itself
            if (il.eq.1) then
              kke=0
              kn=0
              k=0

              if(lwlddir(il)) then !{
c*** CCW target
                do j=1,n !{
                  ke=0
                  do i=1,ntrgts !{
                    if(j.eq.ie(i)) ke=i
                  end do !}
                  kb=0
                  do i=1,ntrgts !{
                    if(j.eq.ib(i)) kb=i
                  end do !}
c                  print '(1x,a,5i5)',':: j,ke,kb,k,kn : ',j,ke,kb,k,kn!###
                  if(kb.gt.0 .and. ke.gt.0) then !{
                    write (*,*) 'b2pwltrg: a segment can appear in ',
     ,                   'the grid edge lists no more than once.'
                    write (*,*) 'Segment ',lwldlis(j,il),',  list ',
     ,                                 il,',  targets ',it(kb),it(ke)
                    write (*,*) 'ltrgedge:'
                    write (*,'(8x,12i5)') (ltrgedge(1,i),i=1,ntrgts)
                    write (*,'(8x,12i5)') (ltrgedge(2,i),i=1,ntrgts)
                    write (*,*) '==> Check the target edge ',
     ,                                 'specification in the DG model.'
                    call xerrab ('Duplicate segment!')
                  end if !}
                  if(k.eq.0) then !{
                    if(ke.ne.0) then !{
                      kke=it(ke)
                      if(j.ne.1) then !{
c*** The zone starts from a target: reset counters and insert the target
                        kn=0
                        if(lim_cnf) then !{
c*** Limiter configuration. Add first the CW target face
                          do i=1,nytg !{
                            luwrk(i)=lppp+i
                            do iss=1,nsts !{
                              if (eirdiag_nds_typ(iss).ne.2) cycle
                              if (eirdiag_nds_srf(iss).ne.
     .                                               ltrgpos(1)+1) cycle
                              if (eirdiag_nds_start(iss).gt.i) cycle
                              if (eirdiag_nds_end(iss).lt.i) cycle
                              ind=eirdiag_nds_ind(iss)
                              ias=ind+1-(eirdiag_nds_start(iss)-1)
                              lewrk(i)=nlim+ias+(i-1)
                            end do !}
                          end do !}
                          kn=kn+nytg
                        end if !}
                        k=it(ke)
                        do i=1,nytg !{
                          luwrk(kn+i)=lppp+nytg*(k-1)+i
                          do iss=1,nsts !{
                            if (eirdiag_nds_typ(iss).ne.2) cycle
                            if (eirdiag_nds_srf(iss).ne.
     .                                           ltrgpos(k)+1) cycle
                            if (eirdiag_nds_start(iss).gt.i) cycle
                            if (eirdiag_nds_end(iss).lt.i) cycle
                            ind=eirdiag_nds_ind(iss)
                            ias=ind+1-(eirdiag_nds_start(iss)-1)
                            lewrk(kn+i)=nlim+ias+(i-1)
                          end do !}
                        end do !}
                        kn=kn+nytg
                        kke=-kke
                        k=0
                      end if !}
                      kn=kn+1
                      luwrk(kn)=lwldlis(j,il)
                      lewrk(kn)=leirlis(j,il)
                    else if(kb.ne.0) then !}{
                      kn=kn+1
                      luwrk(kn)=lwldlis(j,il)
                      lewrk(kn)=leirlis(j,il)
                      k=it(kb)
                      if(k.ne.-kke) then !{
                        if(j.ne.n .or. k.eq.kke) then !{
                          do i=1,nytg !{
                            luwrk(kn+i)=lppp+nytg*(k-1)+i
                            do iss=1,nsts !{
                              if (eirdiag_nds_typ(iss).ne.2) cycle
                              if (eirdiag_nds_srf(iss).ne.
     .                                             ltrgpos(k)+1) cycle
                              if (eirdiag_nds_start(iss).gt.i) cycle
                              if (eirdiag_nds_end(iss).lt.i) cycle
                              ind=eirdiag_nds_ind(iss)
                              ias=ind+1-(eirdiag_nds_start(iss)-1)
                              lewrk(kn+i)=nlim+ias+(i-1)
                            end do !}
                          end do !}
                          kn=kn+nytg
                          if(lim_cnf) then !{
c*** Limiter configuration. Add the CCW target face
                            do i=1,nytg !{
                              luwrk(kn+i)=lppp+nytg+i
                              do iss=1,nsts !{
                                if (eirdiag_nds_typ(iss).ne.2) cycle
                                if (eirdiag_nds_srf(iss).ne.
     .                                               ltrgpos(2)+1) cycle
                                if (eirdiag_nds_start(iss).gt.i) cycle
                                if (eirdiag_nds_end(iss).lt.i) cycle
                                ind=eirdiag_nds_ind(iss)
                                ias=ind+1-(eirdiag_nds_start(iss)-1)
                                lewrk(kn+i)=nlim+ias+(i-1)
                              end do !}
                            end do !}
                            kn=kn+nytg
                          end if !}
                        end if !}
                      end if !}
                    else !}{
                      kn=kn+1
                      luwrk(kn)=lwldlis(j,il)
                      lewrk(kn)=leirlis(j,il)
                    end if !}
                  else !}{
                    if(ke.ne.0) then !{
                      if(it(ke).ne.k .and. .not.lim_cnf) then !{
                        write(*,*) 'b2pwltrg: unexpected error. ',
     ,                        'Trying to close a target different ',
     ,                        'from the opened one.'
                        write(*,*) 'Segment ',lwldlis(j,il),
     ,                                   ',  targets: ',k,it(ke)
                        write(*,*) 'List ',il
                        write(*,'(1x,30i4)') (lwldlis(i,il),i=1,n)
                        write(*,*) 'it : ',(it(i),i=1,kk)
                        write(*,*) 'ib : ',(ib(i),i=1,kk)
                        write(*,*) 'ie : ',(ie(i),i=1,kk)
                        call xerrab ('Not closing the open target!')
                      end if !}
                      kn=kn+1
                      luwrk(kn)=lwldlis(j,il)
                      lewrk(kn)=leirlis(j,il)
                      k=0
                   else if(kb.ne.0) then !}{
                      write(*,*) 'b2pwltrg: unexpected error. ',
     ,                         'Trying to open a target having one ',
     ,                         'already opened.'
                      write(*,*) 'Segment ',lwldlis(j,il),
     ,                                    ',  targets: ',k,it(kb)
                      write(*,*) 'List ',il
                      write(*,'(1x,30i4)') (lwldlis(i,il),i=1,n)
                      write(*,*) 'it : ',(it(i),i=1,kk)
                      write(*,*) 'ib : ',(ib(i),i=1,kk)
                      write(*,*) 'ie : ',(ie(i),i=1,kk)
                      call xerrab ('Target already opened!')
                    end if !}
                  end if !}
                end do !}
                do i=1,kn !{
                  lwldlis(i,il)=luwrk(i)
                  leirlis(i,il)=lewrk(i)
                end do !}
                nwldlis(il)=kn
c                print *,'il,kn = ',il,kn                  !###
              else !}{
                print *,'b2pwltrg: not yet provided for the CW target ',
     ,                                                   'orientation.'
                print *,'==> Check specification of the "Plot zones"',
     ,                            ' in the DG model.'
                call xerrab ('CW target in plot zones not allowed!')
              end if !}
            end if !} il.eq.1
          end if !} kk.gt.0
        end if !} n.gt.0
      end do !}
c!###{
c      print *,'After correction. nwldprt =',nwldprt
c      print '(a,(20i4))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt !{
c       print '(a,i4,1x,a/(20i4))','lwldlis :',il,pltzlbl(il),
c     ,                                 (lwldlis(i,il),i=1,nwldlis(il))
c      end do !}
c!###}

      do il=2,nwldprt
c*** Store original plot zones to be used for substituting proper elements inside targets
        nwldlis_old = nwldlis(il)
        allocate(lwldlis_old(size(lwldlis,1)))
        allocate(leirlis_old(size(leirlis,1)))
        lwldlis_old(:)=lwldlis(:,il)
        leirlis_old(:)=leirlis(:,il)
        allocate(lis_ind(nwldlis_old))
ccc*** Order elements of original plot zone in CCW direction
        do i=1,nwldlis_old
          k1=abs(lwldlis_old(i))
          xx1=pp1(1,k1)
          yy1=pp1(2,k1)
          distmin=1.0e5
          flagmin=1
          do i1=1,nwldlis(1)
            k2=abs(lwldlis(i1,1))
            xx2=pp1(1,k2)
            yy2=pp1(2,k2)
            distance=sqrt((xx1-xx2)**2+(yy1-yy2)**2)
            if (distance<distmin) then !{
              distmin=distance
              flagmin=i1
            end if !}
          end do
          lis_ind(i)=flagmin
        end do
        !Sort
        do i1=1,nwldlis_old-1
          cont=minloc(lis_ind(i1:nwldlis_old),1)
          itmp=lis_ind(cont+i1-1)
          lis_ind(cont+i1-1)=lis_ind(i1)
          lis_ind(i1)=itmp
          itmp=lwldlis_old(cont+i1-1)
          lwldlis_old(cont+i1-1)=lwldlis_old(i1)
          lwldlis_old(i1)=itmp
          itmp=leirlis_old(cont+i1-1)
          leirlis_old(cont+i1-1)=leirlis_old(i1)
          leirlis_old(i1)=itmp
        end do
        deallocate(lis_ind)

ccc*** Check if a plot zone has starting and ending element which is inside a target     
        if (.not. istarget(il)) then 
          !{ check if first and last element of the list are not in lwldlis(1), if so they are inside targets
          k1=abs(lwldlis_old(1));
          k2=abs(lwldlis_old(nwldlis_old))
          flagmin=0
          do i1=1,nwldlis(1)
            ksel=abs(lwldlis(i1,1))
            if ((ksel .eq. k1) .or. (ksel .eq. k2)) then
              flagmin=1
            end if
          end do
          if (flagmin .eq. 0) then
            istarget(il)=.true.
          end if
        end if

ccc*** Substitute B2 cells only when overlapping original plot zone target elements
        if (istarget(il)) then !{ inside target
          k1=lwldlis_old(1)
          k2=lwldlis_old(nwldlis_old)
          distmin=1.0e5
          flagmin=1
          distmin1=1.0e5
          flagmin1=1
          do i=1,nwldlis(1)
            xx1=pp1(1,k1)
            yy1=pp1(2,k1)
            ksel=lwldlis(i,1)
            xx2=pp1(1,ksel)
            yy2=pp1(2,ksel)
            distance=sqrt((xx1-xx2)**2+(yy1-yy2)**2)
            if (distance<distmin) then !{
              distmin=distance
              flagmin=i
            end if !}
            xx1=pp2(1,k2)
            yy1=pp2(2,k2)
            xx2=pp2(1,ksel)
            yy2=pp2(2,ksel)
            distance=sqrt((xx1-xx2)**2+(yy1-yy2)**2)
            if (distance<distmin1) then !{
              distmin1=distance
              flagmin1=i
            end if !}
          end do !}
          cont=1
          do i=flagmin,flagmin1
            lwldlis_old(cont)=lwldlis(i,1)
            leirlis_old(cont)=leirlis(i,1)
            cont=cont+1
          end do
          nwldlis_old=cont-1
          nwldlis(il)=cont-1
c*** Now change also the new plot zone according to the CW or CCW order
          if (lwlddir(il)) then !CCW order
            do i=1,nwldlis(il)
              lwldlis(i,il)=lwldlis_old(i)
              leirlis(i,il)=leirlis_old(i)
            end do
          else !CW order
            do i=1,nwldlis(il)
              lwldlis(i,il)=lwldlis_old(nwldlis(il)-i+1)
              leirlis(i,il)=leirlis_old(nwldlis(il)-i+1)
            end do
          end if
        end if !}
      deallocate(lwldlis_old)
      deallocate(leirlis_old)
      end do !}

      deallocate(istarget) 


c*** Calculate the segment areas and related quantities

      do i=1,lppp+ntrgts*nytg !{
        segm_len(i)=sqrt((pp1(1,i)-pp2(1,i))**2+(pp1(2,i)-pp2(2,i))**2)
        if (isymm.eq.0) then
          segm_area(i) = width*segm_len(i)
        else if (isymm.eq.1 .or. isymm.eq.2) then
          segm_area(i) = pi*(pp1(1,i)+pp2(1,i))*segm_len(i)
        else if (isymm.eq.3 .or. isymm.eq.4) then
          segm_area(i) = pi*(pp1(2,i)+pp2(2,i))*segm_len(i)
        end if
      end do !}
      do il=1,nwldprt
        d=0.5_R8*segm_len(iabs(lwldlis(1,il)))
        x=d
        wldrfdst(1,il)=x
        do i=2,nwldlis(il)
          x=x+d
          d=0.5_R8*segm_len(iabs(lwldlis(i,il)))
          x=x+d
          wldrfdst(i,il)=x
        end do
      end do

c*** Produce the data files containing the geometry of the shadowing
c*** structure and of the plot zones

c*** The whole wall including the refined targets (shadowing)

      open(lun_w,file='shadow.dat')
      write(lun_w,'(2a)') '# Wall layout : ',trim(label)
      write(lun_w,'(3a)') '# R',tab,'Z'
      do i=1,nwldlis(1) !{
        k=iabs(lwldlis(i,1))
        write(lun_w,'(1p,2(1x,e15.8))') pp1(1,k),pp1(2,k)
      end do !}
      write(lun_w,'(1p,2(1x,e15.8))') pp2(1,k),pp2(2,k)

c*** The plot zones: one file for each zone. Here we also write the
c*** coordinates of the segment centers and the reference distance

      do il=2,nwldprt !{
        l=len_trim(pltzlbl(il))
c*** create the filename for the output
        if(l.eq.0) then !{
c*** use the zone index to label the plot zone
          write(fn_hlp,'(a,i2.2,a)') 'wlly_',il-1,'.dat'
        else !}{
          fn_hlp='wlly_'//pltzlbl(il)(:l)//'.dat'
        end if !}
        open(lun_w,file=fn_hlp)
        write(lun_w,'(4a)') '#Layout of plot zone ',
     ,                                      pltzlbl(il)(:l),' : ',label
        write(lun_w,'(11a)') '#R',tab,'Z',tab,'Rc',tab,'Zc',tab,'x',
     ,    tab,'DG_elm'
        do i=1,nwldlis(il) !{
          k=iabs(lwldlis(i,il))
          write(lun_w,'(1p,e12.5,4(a1,e12.5),a1,i5)') pp1(1,k),
     ,      tab,pp1(2,k),
     ,      tab,0.5_R8*(pp1(1,k)+pp2(1,k)),
     ,      tab,0.5_R8*(pp1(2,k)+pp2(2,k)),
     ,      tab,wldrfdst(i,il),tab,k
        end do !}
        write(lun_w,'(1p,e12.5,4(a1,e12.5))') pp2(1,k),tab,pp2(2,k)
      end do !}
      close(lun_w)
#ifdef DBG
      print *,'Leaving b2pwltrg. lytg,nytg= ',lytg,nytg
      print *,'nwldprt, nwldlis, lwlddir = ',nwldprt
      print '(16i5)',(nwldlis(i),i=1,nwldprt)
      print '(16l5)',(lwlddir(i),i=1,nwldprt)
      do il=1,nwldprt
        print '(a,i5,1x,a/(1x,8i11))','lwldlis :',il,pltzlbl(il),
     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
        print '(a,i5/(1x,8i11))','leirlis :',il,
     ,                                  (leirlis(i,il),i=1,nwldlis(il))
        print '(a,i5/(1x,1p,8e11.3))','length  :',il,
     ,                        (segm_len(lwldlis(i,il)),i=1,nwldlis(il))
        print '(a,i5/(1x,1p,8e11.3))','distance:',il,
     ,                                 (wldrfdst(i,il),i=1,nwldlis(il))
      end do
#endif
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
      print *,'%% b2pwltrg: nwldrad,nytg / lwldrad',nwldrad,nytg
      k=0
      do j=1,ntrgts !{
        print '(20i4)',(lwldrad(k+i),i=1,nytg)
        k=k+nytg
      end do !}
      print '(1x,a,i5/(20i4))','nwldlis / lwldlis',nwldlis(1),
     ,                           (lwldlis(i,1),i=1,nwldlis(1))
#endif

c======================================================================
      call subend ()
      return
      end subroutine b2pwltrg
