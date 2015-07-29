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
c***  lwldlis : sorted lists of the "plot zones"
c***  lwlddir : orientation of the "plot zones" (true = CCW)
c***  pp1,pp2 : completed arrays of the wall segments
c======================================================================
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_constants
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#include <DIMENSIONS.F>
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
      real (kind=R8) ::
     *   d,x
      integer i,j,k,l,n,i1,i2,ii,il,is,kb,ke,kk,ks,kkb,kke,kks,
     ,        n_atg,n_targ,
     ,        ib(2),ie(2),ibb(2),iee(2),iit(2),l_atg(2)
      logical add_trg
c
c*** wsrttol : a tolerance parameter for defining coinciding segments
      real (kind=R8) ::
     *   wsrttol
      parameter (wsrttol=2.0e-4_R8)
c<<<
c      integer lun_w
c      parameter (lun_w=76)
c      character tab
c>>>
      intrinsic mod
c======================================================================
c<<<
c      tab=char(9)
c
      print *,'Starting b2pwlprp...',nlim,lppp,DEF_NLIM
      print '(a,4i6)','jedgi1, jedgi2, jedgo1, jedgo2 =',jedge(1:4)
c>>>
c
c*** Check the presence of the shadowing structure
c
      n=nwldlis(1)
      if(n .le. 0) then
        write(*,*) 'b2pwlprp: no shadowing structure in the data!'
        call walloff
        l2pwinp=.false.
        return
      end if
c
c*** Convert the input wall coordinates into double precision
c*** and scale them to meters.
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
      print *,'nwldprt =',nwldprt
      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
      do il=1,nwldprt
        print '(a,i5/(5x,15i5))','lwldlis :',il,
     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
      end do
c
c      open(lun_w,file='layout-d.dat')
c      write(lun_w,'(6a)') '#Wall layout (DG) : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,lppp
c        write(lun_w,*) pp1(1,i),tab,pp1(2,i)
c        write(lun_w,*) pp2(1,i),tab,pp2(2,i)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>
c
c*** Sort the segments in the lists of the "wall zones"
c*** The shadowing structure must be closed, therefore it starts from
c*** an arbitrary point.
c
      lwldstrt(1)=iabs(lwldlis(1,1))
      call b2pwsort
c
c*** The sorting routine, b2pwsort, also checks the wall segments. If a
c*** segment input from DG does not coincide with the corresponding
c*** segment saved by Eirene, its number becomes negative, and no
c*** neutral-related data will be plotted for it.
c
c<<<
c      print *,'After sorting. nwldprt =',nwldprt
c      print '(a,(20i5))','lwldstrt:',(lwldstrt(i),i=1,nwldprt)
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      print '(a,(20l5))','lwlddir :',(lwlddir(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5/(5x,15i5))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c      end do
c>>>
c*** Check the shadowing structure - it must be closed
c
      n=nwldlis(1)
      l=iabs(lwldlis(1,1))
      k=iabs(lwldlis(n,1))
      if(lwlddir(1)) then
        if(abs(pp1(1,l)-pp2(1,k)) .gt. wsrttol  .or. 
     .     abs(pp1(2,l)-pp2(2,k)) .gt. wsrttol) then
          write(*,*) 'b2pwlprp: the shadowing structure must be closed!'
          write(*,*) 'lwlddir(1) = ',lwlddir(1)
          write(*,*) 'Starting segment: ',l,pp1(1,l),pp1(2,l)
          write(*,*) 'Ending segment  : ',k,pp2(1,k),pp2(2,k)
          call walloff
          l2pwinp=.false.
          return
        end if
      else
        if(abs(pp1(1,k)-pp2(1,l)) .gt. wsrttol  .or. 
     .     abs(pp1(2,k)-pp2(2,l)) .gt. wsrttol) then
          write(*,*) 'b2pwlprp: the shadowing structure must be closed!'
          write(*,*) 'lwlddir(1) = ',lwlddir(1)
          write(*,*) 'Starting segment: ',l,pp2(1,l),pp2(2,l)
          write(*,*) 'Ending segment  : ',k,pp1(1,k),pp1(2,k)
          call walloff
          l2pwinp=.false.
          return
        end if
      end if
c
c
c*** Add the target sections from B2 grid to the segment lists
c
cc<<<
c      print *,'Grid edges'
c      print '(a4,8a12)','j','r1i','r2i','z1i','z2i',
c     ,                      'r1o','r2o','z1o','z2o'
cc>>>
      do i=0,ny-1
        pp1(1,i+1+lppp)=crx(0,ny-i,2)
        pp1(2,i+1+lppp)=cry(0,ny-i,2)
        pp2(1,i+1+lppp)=crx(0,ny-i,0)
        pp2(2,i+1+lppp)=cry(0,ny-i,0)
        pp1(1,i+1+lppp+ny)=crx(nx,i,0)
        pp1(2,i+1+lppp+ny)=cry(nx,i,0)
        pp2(1,i+1+lppp+ny)=crx(nx,i,2)
        pp2(2,i+1+lppp+ny)=cry(nx,i,2)
      end do
c<<<
c      print '(1x,i4,1p,8e12.4)',(i,(pp1(j,i+1+lppp),pp2(j,i+1+lppp),j=1,2),
c     ,                (pp1(j,i+1+lppp+ny),pp2(j,i+1+lppp+ny),j=1,2),i=0,ny-1)
c
c      open(lun_w,file='layout-g.dat')
c      write(lun_w,'(6a)') '#Wall layout (g) : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,lppp+2*ny
c        write(lun_w,*) pp1(1,i),tab,pp1(2,i)
c        write(lun_w,*) pp2(1,i),tab,pp2(2,i)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>
c
c*** Modify the target edges for continuity with B2 grid edges
c
      pp2(1,jedge(1))=crx(0,0,0)
      pp2(2,jedge(1))=cry(0,0,0)
      pp2(1,jedge(2))=crx(nx-1,0,1)
      pp2(2,jedge(2))=cry(nx-1,0,1)
      pp1(1,jedge(3))=crx(0,ny-1,2)
      pp1(2,jedge(3))=cry(0,ny-1,2)
      pp1(1,jedge(4))=crx(nx-1,ny-1,3)
      pp1(2,jedge(4))=cry(nx-1,ny-1,3)

c<<<
c      open(lun_w,file='layout-h.dat')
c      write(lun_w,'(6a)') '#Wall layout (h) : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,lppp+2*ny
c        write(lun_w,*) pp1(1,i),tab,pp1(2,i)
c        write(lun_w,*) pp2(1,i),tab,pp2(2,i)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>
c
c*** Check whether the modified edges coincide with the saved ones
c
      do k=1,4
        j=jedge(k)
c        print *,'%% Check ',k,j
        if(abs(p1(1,j)-pp1(1,j)) .le. wsrttol  .and.
     .     abs(p1(2,j)-pp1(2,j)) .le. wsrttol  .and.
     .     abs(p2(1,j)-pp2(1,j)) .le. wsrttol  .and.
     .     abs(p2(2,j)-pp2(2,j)) .le. wsrttol) then
c          print *,'== OK'
          do il=1,nwldprt
            n=nwldlis(il)
            do i=1,n
              if(iabs(lwldlis(i,il)) .eq. j) lwldlis(i,il)=j
            end do
          end do
        end if
      end do
c
c<<<
c      print *,jedge(1:4)
c      print '((2(1x,a3,1p,4e12.4)))',
c     ,          'i1',(pp1(j,jedge(1)),pp2(j,jedge(1)),j=1,2),
c     ,          'o1',(pp1(j,jedge(3)),pp2(j,jedge(3)),j=1,2),
c     ,          'i2',(pp1(j,jedge(2)),pp2(j,jedge(2)),j=1,2),
c     ,          'o2',(pp1(j,jedge(4)),pp2(j,jedge(4)),j=1,2)
c>>>
c
c*** Correct the segment lists for the refined targets
c*** All the segments between the target edges - or between a target
c*** edge and the end of the list - are replaced with the segments
c*** of the target sections of the B2 grid.
c
c*** Below, kb marks the target; ib, ie are the locations of the inner
c*** and outer edges in the lists
c
      do il=1,nwldprt
        n=nwldlis(il)
        if(n .gt. 0) then
          do i=1,2
            ib(i)=0
            ie(i)=0
          end do
          kb=0
          if(lwlddir(il)) then
            is=1
            ii=1
          else
            is=n
            ii=-1
          end if
          do i=1,n
            k=iabs(lwldlis(is,il))
            if(k .eq. jedge(1)  .or.  k .eq. jedge(2)) then
              if(kb .ne. 0) then
                write(*,'(a,a,3i6)') 'b2pwlprp: unexpected error. ',
     ,           'Outer edge not found ',il,i,kb
                write(*,'(a,5i6)')
     ,           'k,jedgi1,jedgi2,jedgo1,jedgo2 = ',k,jedge(1:4)
                stop
              end if
              if(k .eq. jedge(1)) then
                kb=1
              else
                kb=2
              end if
              ib(kb)=is
            else if(k .eq. jedge(3)  .or.  k .eq. jedge(4)) then
              if(kb .eq. 0) then
                if(k .eq. jedge(3)) then
                  kb=1
                else
                  kb=2
                end if
              else
                if((k .eq. jedge(3)  .and.  kb .ne. 1)  .or. 
     ,           (k .eq. jedge(4)  .and.  kb .ne. 2)) then
                  write(*,'(a,a,3i6)')
     ,             'b2pwlprp: unexpected error. ',
     ,             'Inner edge not found ',il,i,kb
                  write(*,'(a,5i6)')
     ,             'k,jedgi1,jedgi2,jedgo1,jedgo2 = ',k,jedge(1:4)
                  stop
                end if
              end if
              ie(kb)=is
              kb=0
            end if
            is=is+ii
          end do
c
c          if(kb .gt. 0) then
cc*** A target left open - consider it spanning through the edge
c            ie(kb)=is
c          end if
c
c*** Replace the target section with the B2 grid edge.
c*** First, check the sequence. We have here the indices ib and ie,
c*** marking the positions of the target edges in the lists.
c*** ib marks the inner edges (1 for inner 2 for outer targets), and
c*** ie marks the outer edges.
c*** We need to arrange new counters, ibb, iee, and iit, which would be
c*** sorted along the list direction and mark
c***  ibb   the first edge
c***  iee   the matching edge
c***  itt   the corresponding target
c*** for each passible target, that is (depending on the orientation)
c***            CCW             CW
c*** ibb =  ib(1),ib(2)     ie(2),ie(1)
c*** iee =  ie(1),ie(2)     ib(2),ib(1)
c*** iit =      1,2             2,1
c***
c*** iit=0 : corresponding target does not belong to the list
c*** ibb=0 : target starts from the begining of the list
c*** iee=0 : target continues over the end of the list
c***
c*** The list orientation is given by lwlddir value (true=CCW)
c
          k=0
          iit(1)=0
          iit(2)=0
          ibb(2)=0
          if(lwlddir(il)) then
            do i=1,2
              if(ib(i)+ie(i) .gt. 0) then
                k=k+1
                iit(k)=i
                ibb(k)=ib(i)
                iee(k)=ie(i)
              end if
            end do
          else
            do i=2,1,-1
              if(ib(i)+ie(i) .gt. 0) then
                k=k+1
                iit(k)=i
                ibb(k)=ie(i)
                iee(k)=ib(i)
              end if
            end do
          end if
          kk=k
c<<<
c          print *,'After the check. il,kk,lppp = ',il,kk,lppp
c          print '(5a6)','ib','ie','ibb','iee','iit'
c          print '(5i6)',(ib(i),ie(i),ibb(i),iee(i),iit(i),i=1,2)
c>>>
c
c*** Now fill in the necessary arrays
c
          if(kk .gt. 0) then
c
c*** We have here the arrays ibb and iee containing the starting
c*** and ending segments for each divertor plate - as seen along the
c*** target under consideration.
c***
c*** First, sort the ibb's and iee's
c
            if(ibb(1) .le. ibb(2)) then
              ib(1)=1
              ib(2)=2
            else
              ib(1)=2
              ib(2)=1
            end if
            if(iee(1) .le. iee(2)) then
              ie(1)=1
              ie(2)=2
            else
              ie(1)=2
              ie(2)=1
            end if
            do i=1,2
              if(ibb(ib(1)) .eq. 0) then
                ib(1)=ib(2)
                ib(2)=0
              end if
              if(iee(ie(1)) .eq. 0) then
                ie(1)=ie(2)
                ie(2)=0
              end if
            end do
c
c*** ... make a temporary copy of the segment list ...
c
            do i=1,n
              luwrk(i)=lwldlis(i,il)
            end do
c
c*** ... and go along the list analysing the intersections. Here
c***   k       is the counter for the new list
c***   kk      is the counter for the old list
c***   kb      if non-zero, the target number kb is "open"
c***   ks      is the counter for the target starts
c***   ke      is the counter for the target ends
c***           very begining of the list (needed to avoid duplication)
c***   l_atg   list of target sections already added to the segment list
c***   n_atg   length of this list
c
            n_targ=2
            k=0
            kb=0
            kkb=0
            ks=1
            ke=1
            kk=1
            n_atg=0
c<<<
c            print *,         'Before merging the targets.  ',
c     ,                      ':  ibb  :  iee  :  ib   :  ie   :  iit  :'
c            print '(a,10i4)','                             ',
c     ,                      ibb,iee,ib,ie,iit
c>>>

            do ii=1,2*n_targ

c<<<
c              print *,' ii   k  kb  ks  ke  kk n_atg,l_atg'
c              print '(12i4)',ii,k,kb,ks,ke,kk,n_atg,l_atg(1),l_atg(2)
c>>>
              if(ks+ke .gt. 0) then
                add_trg=.false.
                kks=nsgmx+5
                if(ks .gt. 0) then
                  kks=ibb(ib(ks))
                end if
                if(kks .eq. 0) then
                  kks=nsgmx+5
                end if
                kke=nsgmx+5
                if(ke .gt. 0) then
                  kke=iee(ie(ke))
                end if
                if(kke .eq. 0) then
                  kke=nsgmx+5
                end if
c                print *,'kks,kke = ',kks,kke
                if(kks .lt. kke) then
c*** entering a target
                  if(kb .ne. 0) then
                    write(*,*) 'Unexpected error in b2pwlprp: ',
     ,               'the target ',kb,'  not closed. Overlap?'
                    stop
                  end if
                  kb=iit(ib(ks))
                  do i=kk,kks
                    k=k+1
                    lwldlis(k,il)=luwrk(i)
                  end do
                  kk=n+1
                  add_trg=kks .lt. n
                  kkb=kb
                  ks=mod(ks+1,n_targ+1)
                else if(kks .eq. kke  .and.  kks .le. n) then
                  write(*,*) 'Unexpected error in b2pwlprp: ',
     ,             'kks = kke = ',kks
                    stop
                else
                  if(kb .eq. 0) then
                    if(n_atg .gt. 0) then
                      write(*,*) 'Unexpected error in b2pwlprp: ',
     ,                 'trying to close unopened target ',iit(ie(ke))
                      stop
                    end if
                    kb=iit(ie(ke))
                    add_trg=kke .gt. 1
                  else
                    add_trg=.true.
                  end if
                  kkb=0
                  kk=kke
                  ke=mod(ke+1,n_targ+1)
                end if
                if(add_trg) then
                  do i=1,n_atg
                    add_trg= add_trg  .and.  kb .ne. l_atg(i)
                  end do
                  if(add_trg) then
c
c*** Add the target section to the new list
c
                    n_atg=n_atg+1
                    l_atg(n_atg)=kb
                    if(lwlddir(il)) then
                      i1=lppp
                      i2=1
                    else
                      i1=lppp+ny+1
                      i2=-1
                    end if
                    l=i1+(kb-1)*ny
c                    print *,'k,i1,i2,l = ',k,i1,i2,l
                    do i=0,ny-1
                      l=l+i2
                      k=k+1
                      lwldlis(k,il)=l
                    end do
                  end if
                end if
                kb=kkb
              end if
            end do
c
c*** Add the last section if left
c
            do i=kk,n
              k=k+1
              lwldlis(k,il)=luwrk(i)
            end do

            nwldlis(il)=k
          end if
        end if
      end do
c<<<
c      print *,'After correction. nwldprt =',nwldprt
c      print '(a,(20i5))','nwldlis :',(nwldlis(i),i=1,nwldprt)
c      do il=1,nwldprt
c        print '(a,i5/(5x,15i5))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c      end do
c
c      open(lun_w,file='layout-f.dat')
c      write(lun_w,'(6a)') '#Wall layout (f) : ',label
c      write(lun_w,'(6a)') '#R',tab,'Z'
c      do i=1,nwldlis(1)
c        k=iabs(lwldlis(i,1))
c        write(lun_w,*) pp1(1,k),tab,pp1(2,k)
c        write(lun_w,*) pp2(1,k),tab,pp2(2,k)
c        write(lun_w,*)
c      end do
c      close(lun_w)
c>>>
c
c*** Calculate the segment areas and related quantities
c
      do i=1,lppp+2*ny
        segm_len(i)=sqrt((pp1(1,i)-pp2(1,i))**2+(pp1(2,i)-pp2(2,i))**2)
        segm_area(i) = pi*(pp1(1,i)+pp2(1,i))*segm_len(i)
      end do
      do il=2,nwldprt
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
c<<<
c      print *,'Leaving b2pwlprp. '
c      print *,'nwldprt, nwldlis, lwlddir = ',nwldprt
c      print '(16i5)',(nwldlis(i),i=1,nwldprt)
c      print '(16l5)',(lwlddir(i),i=1,nwldprt)
c      do il=2,nwldprt
c        print '(a,i5/(1x,8i11))','lwldlis :',il,
c     ,                                  (lwldlis(i,il),i=1,nwldlis(il))
c        print '(a,i5/(1x,1p,8e11.3))','length  :',il,
c     ,                        (segm_len(lwldlis(i,il)),i=1,nwldlis(il))
c        print '(a,i5/(1x,1p,8e11.3))','distance:',il,
c     ,                                 (wldrfdst(i,il),i=1,nwldlis(il))
c      end do
c>>>
c======================================================================
      return
      end subroutine b2pwlprp
