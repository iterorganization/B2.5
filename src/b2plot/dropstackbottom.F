      subroutine dropstackbottom()
      use b2mod_types
      use b2mod_b2plot
      implicit none
cxpb
      integer i, ix, iy, is, ixx, iyy
      intrinsic mod
      ixx(ix,iy)=mod(ix+1,nx+2)-1
      iyy(ix,iy)=iy+(ix+1)/(nx+2)
cxpb

      write(*,*) 'Too many data values on work stack'
cxpb
      write(*,*) 'Dropping the bottom element!'
      do i = 2,iwrk
        ndim(1,i-1) = ndim(1,i)
        ndim(2,i-1) = ndim(2,i)
        ndim(3,i-1) = ndim(3,i)
        do is = 0,ndim(3,i)-1
          do iy = -1,ndim(2,i)
            do ix = -1,ndim(1,i)
              work(ixx(ix,iy),iyy(ix,iy),is,i-1) = 
     .         work(ixx(ix,iy),iyy(ix,iy),is,i) 
            end do
          end do
          textwk(is,i-1) = textwk(is,i)
        end do
        graphlabel(i-1) = graphlabel(i)
      end do
      iwrk = iwrk-1
cxpb
      return
      end

      subroutine dropwallstackbottom()
      use b2mod_types
      use b2mod_b2plot
      implicit none
cxpb
      integer i, iw, is

      write(*,*) 'Too many data values on wall stack'
cxpb
      write(*,*) 'Dropping the bottom element!'
      do i = 2,iwwrk
        nwdim(1,i-1) = nwdim(1,i)
        nwdim(2,i-1) = nwdim(2,i)
        do is = 0,nwdim(2,i)-1
          do iw = 1,nwdim(1,i)
            wall_stack(iw,is,i-1) = wall_stack(iw,is,i) 
          end do
          textwwk(is,i-1) = textwwk(is,i)
        end do
        wall_label(i-1) = wall_label(i)
      end do
      iwwrk = iwwrk-1
cxpb
      return
      end

      subroutine droptargetstackbottom()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      implicit none
cxpb
      integer i, iw, id

      write(*,*) 'Too many data values on target stack'
cxpb
      write(*,*) 'Dropping the bottom element!'
      do i = 2,itwrk
        do id = 1,ndepth_nml
          do iw = 1,nbnd
            target_stack(iw,id,i-1) = target_stack(iw,id,i) 
          end do
        end do
        target_label(i-1) = target_label(i)
      end do
      itwrk = itwrk-1
cxpb
      return
      end
