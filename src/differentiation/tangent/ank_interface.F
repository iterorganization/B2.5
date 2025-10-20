      subroutine ank_interface(inum)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_ank_diffv
      use b2mod_wall
      implicit none
      integer inum

      integer i,in
      real (kind=R8) :: dr, dz

      in=max(min(inum,nwall),2)
      dr=1.0_R8/in
      dz=0.0_R8/in
      r1(1)=1.0_R8
      z1(1)=0.0_R8
      r2(1)=r1(1)+dr
      z2(1)=z1(1)+dz
      do i=2,in
        r1(i)=r2(i-1)
        z1(i)=z2(i-1)
        r2(i)=r1(i)+dr
        z2(i)=z1(i)+dz
        data(i)=i
      enddo
      ndata=in
      write(wall_label,*) 'wall: ',ndata
      return
      end subroutine ank_interface
