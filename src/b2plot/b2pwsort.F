      subroutine b2pwsort
c
c  version : 20.08.97 18:14
c
c======================================================================
c*** Sorts the lists of target sections to simply connected chains
c======================================================================
c*** Global variables
c***
c***  nwldprt : the number of lists to be sorted
c***  lwldstrt: the initial segments of the chains
c***  nwldlis : the list lengths
c***  lwldlis : the list themselves
c***  lwlddir : the orientation of the chain (true = CCW)
c
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_wall_loading
#ifdef B25_EIRENE
      use eirmod_parmmod
      use eirmod_extrab25
#endif
      implicit none
#ifndef B25_EIRENE
#     include "COUPLE/KOPPLDIM.F"
#endif
c
c*** Local variables
c***
c*** goup  : the end of the last segment of the chain is still free
c*** placed: the segment is already in the chain
c*** lwsrtwrk: a temporary array for sorting
c*** nwsrtwrk: a counter for this array
c*** wsrttol : a tolerance parameter for defining adjacent segments
c
      logical goup,placed(nlim)
      integer lwsrtwrk(nlim),il,is,ii,jj,js,nwsrtwrk,k,ik,j
      real (kind=R8) ::
     *   wsrttol
      parameter (wsrttol=2.0e-4_R8)
c======================================================================
c     integer i
c
cc%%%
c      print *,'b2pwsrt',nwldprt
c      do il=1,nwldprt
c        print '(1x,20i4)',il,nwldlis(il),lwldstrt(il)
c        print '(1x,20i4)',(lwldlis(i,il),i=1,nwldlis(il))
c        do ii=1,nwldlis(il)
c          is=lwldlis(ii,il)
c          print '(4f12.4)',pp1(1,is),pp1(2,is),pp2(1,is),pp2(2,is)
c        end do
c      end do
cc%%%
      do il=1,nwldprt
        is=lwldstrt(il)
        if(is .gt. nlim  .or.  is .le. 0) then
          write (*,*) 'Starting segment outside the limits',is,il
          nwldlis(il)=0.0_R8
        else
          goup=.false.
          do ii=1,nwldlis(il)
            placed(ii)=.false.
            if(is .eq. lwldlis(ii,il)) then
              placed(ii)=.true.
              goup=.true.
            end if
          end do
          if(.not.goup) then
            write (*,*) 'Starting segment outside the list',is,il
            nwldlis(il)=0.0_R8
          else
c*** initialise the chain counter and the chain list
            nwsrtwrk=1
            lwsrtwrk(nwsrtwrk)=is
c*** both ends of the first element will be checked (goup=.true. here)
            k=2
            do 100 ii=2,nwldlis(il)
c              print *,'do 100:',ii,is,goup
              do ik=1,k
c                print *,'do ik:',ik,is,goup
                do jj=1,nwldlis(il)
c                  print *,'do jj:',jj,is,goup,placed(jj)
                  if(.not.placed(jj)) then
                    js=lwldlis(jj,il)
c                    print *,'js:',js,is,goup
                    if(goup) then
                      if(abs(pp2(1,is)-pp1(1,js)) .le. wsrttol  .and. 
     .                 abs(pp2(2,is)-pp1(2,js)) .le. wsrttol) then
                        goup=.true.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'2->1:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp2(1,is)-pp2(1,js)) .le. wsrttol  .and. 
     .                 abs(pp2(2,is)-pp2(2,js)) .le. wsrttol) then
                        goup=.false.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'2->2:',is,nwsrtwrk
                        go to 100
                      end if
                    else
                      if(abs(pp1(1,is)-pp2(1,js)) .le. wsrttol  .and. 
     .                 abs(pp1(2,is)-pp2(2,js)) .le. wsrttol) then
                        goup=.false.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'1->2:',is,nwsrtwrk
                        go to 100
                      end if
                      if(abs(pp1(1,is)-pp1(1,js)) .le. wsrttol  .and. 
     .                 abs(pp1(2,is)-pp1(2,js)) .le. wsrttol) then
                        goup=.true.
                        is=js
                        nwsrtwrk=nwsrtwrk+1
                        lwsrtwrk(nwsrtwrk)=is
                        placed(jj)=.true.
c                        print *,'1->1:',is,nwsrtwrk
                        go to 100
                      end if
                    end if
                  end if
                end do
c*** This only works for the first element.
c*** Look in the opposite direction
                goup=.false.
              end do
c*** No connection found
              go to 50
 100        k=1
 50         if(nwsrtwrk .ne. nwldlis(il)) then
              write (*,*) 'The chain is broken: ',il
              nwldlis(il)=nwsrtwrk
            end if
            do ii=1,nwldlis(il)
              j=lwsrtwrk(ii)
              if(abs(p1(1,j)-pp1(1,j)) .le. wsrttol  .and. 
     .         abs(p1(2,j)-pp1(2,j)) .le. wsrttol  .and. 
     .         abs(p2(1,j)-pp2(1,j)) .le. wsrttol  .and. 
     .         abs(p2(2,j)-pp2(2,j)) .le. wsrttol) then
                lwldlis(ii,il)=j
              else
                lwldlis(ii,il)=-j
              end if
            end do
c            print '(a,i2/(1x,20i4))','il =',il,
c     ,                                (lwldlis(ii,il),ii=1,nwldlis(il))
          end if
          lwlddir(il)=goup
        end if
      end do
      return
c======================================================================
      end subroutine b2pwsort
