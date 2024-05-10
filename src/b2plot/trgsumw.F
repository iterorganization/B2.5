c
      subroutine sumtrgw()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,ibnd
      if(ndepth_nml.le.1) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do i=1,ndepth_nml
        do ibnd=2,nbnd
          target_stack(1,i,itwrk)=
     .     target_stack(1,i,itwrk)+target_stack(ibnd,i,itwrk)
        enddo
      enddo
      values(32)=target_stack(1,1,itwrk)
      values(33)=target_stack(1,1,itwrk)
      do i=1,ndepth_nml
        do ibnd=2,nbnd
          target_stack(ibnd,i,itwrk)=target_stack(1,i,itwrk)
        enddo
        values(32)=min(values(32),target_stack(1,i,itwrk))
        values(33)=max(values(33),target_stack(1,i,itwrk))
      enddo

      return
      end subroutine sumtrgw
c
      subroutine sumtrgz()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer i,ibnd
      if(nbnd.le.1) then
        write(*,*) 'To apply this operator, a two dimensional array ',
     1   'must be specified'
        return
      endif
c
c apply operator
c
      do ibnd=1,nbnd
        do i=2,ndepth_nml
          target_stack(ibnd,1,itwrk)=
     .     target_stack(ibnd,1,itwrk)+target_stack(ibnd,i,itwrk)
        enddo
      enddo
      values(32)=target_stack(1,1,itwrk)
      values(33)=target_stack(1,1,itwrk)
      do ibnd=1,nbnd
        do i=2,ndepth_nml
          target_stack(ibnd,i,itwrk)=target_stack(ibnd,1,itwrk)
        enddo
        values(32)=min(values(32),target_stack(ibnd,1,itwrk))
        values(33)=max(values(33),target_stack(ibnd,1,itwrk))
      enddo

      return
      end subroutine sumtrgz
c
      subroutine inttrgw()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ibnd,iseg,i
      intrinsic min, max
c
c apply operator
c
      do i = 1, ndepth_nml
        do iseg = 1, nsegments
          do ibnd=wall_begin(iseg)+1,wall_end(iseg)
            target_stack(ibnd,i,itwrk)=
     .       target_stack(ibnd-1,i,itwrk)+
     .       target_stack(ibnd  ,i,itwrk)
          enddo
        enddo
      enddo
      values(32)=target_stack(1,1,itwrk)
      values(33)=target_stack(1,1,itwrk)
      do i=1,ndepth_nml
        do ibnd=1,nbnd
          values(32)=min(values(32),target_stack(ibnd,i,itwrk))
          values(33)=max(values(33),target_stack(ibnd,i,itwrk))
        enddo
      enddo

      return
      end subroutine inttrgw
c
      subroutine inttrgwstar()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ibnd,i
      intrinsic min, max
c
c apply operator
c
      values(32)=target_stack(1,1,itwrk)
      values(33)=target_stack(1,1,itwrk)
      do i=1,ndepth_nml
        do ibnd=2,nbnd
          target_stack(1,i,itwrk)=
     .      target_stack(1,i,itwrk)+target_stack(ibnd,i,itwrk)
        enddo
        values(32)=min(values(32),target_stack(1,i,itwrk))
        values(33)=max(values(33),target_stack(1,i,itwrk))
      enddo

      return
      end subroutine inttrgwstar
c
      subroutine inttrgz()
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
      use b2mod_b2plot_values
      implicit none
      integer ibnd,i
      intrinsic min, max
c
c apply operator
c
      do ibnd=1,nbnd
        do i = 2, ndepth_nml
          target_stack(ibnd,i,itwrk)=
     .       target_stack(ibnd,i-1,itwrk)+
     .       target_stack(ibnd,i  ,itwrk)
        enddo
      enddo
      values(32)=target_stack(1,1,itwrk)
      values(33)=target_stack(1,1,itwrk)
      do i = 1, ndepth_nml
        do ibnd=1,nbnd
          values(32)=min(values(32),target_stack(ibnd,i,itwrk))
          values(33)=max(values(33),target_stack(ibnd,i,itwrk))
        enddo
      enddo

      return
      end subroutine inttrgz
