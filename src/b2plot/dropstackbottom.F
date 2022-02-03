      subroutine dropstackbottom()
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer i, ix, iy, is

      write(*,*) 'Too many data values on work stack'
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
      return

      contains

      integer function ixx(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      ixx = mod(ix+1,nx+2)-1
      return
      end function ixx

      integer function iyy(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      iyy = iy+(ix+1)/(nx+2)
      return
      end function iyy

      end subroutine dropstackbottom

      subroutine dropeirenestackbottom()
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer i, iw, is

      write(*,*) 'Too many data values on Eirene stack'
      write(*,*) 'Dropping the bottom element!'
      do i = 2,iewrk
        nedim(1,i-1) = nedim(1,i)
        nedim(2,i-1) = nedim(2,i)
        do is = 0,nedim(2,i)-1
          do iw = 1,nedim(1,i)
            eirene_stack(iw,is,i-1) = eirene_stack(iw,is,i) 
          end do
          textewk(is,i-1) = textewk(is,i)
        end do
        eirene_label(i-1) = eirene_label(i)
      end do
      iewrk = iewrk-1
      return
      end subroutine dropeirenestackbottom

      subroutine dropwallstackbottom()
      use b2mod_types
      use b2mod_b2plot
      implicit none
      integer i, iw, is

      write(*,*) 'Too many data values on wall stack'
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
      return
      end subroutine dropwallstackbottom

      subroutine droptargetstackbottom()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      implicit none
      integer i, iw, id

      write(*,*) 'Too many data values on target stack'
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
      return
      end subroutine droptargetstackbottom
