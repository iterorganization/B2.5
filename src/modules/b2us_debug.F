      module b2us_debug

      use b2mod_types

      implicit none

      private

      public :: to_struct_cell, to_struct_face

      integer, public :: nx_st, ny_st

      contains

      subroutine to_struct_cell(nCv, fun_us, fun_st)
      implicit none
      integer :: nCv, ix, iy, ic
      real (kind=R8) :: fun_us(nCv), fun_st(-1:nx_st,-1:ny_st)

      ic = 1
      do iy = 0, ny_st-1
        do ix = 0, nx_st-1
          fun_st(ix,iy) = fun_us(ic)
          ic = ic + 1
        enddo
      enddo
      do ix = 0, nx_st-1
        fun_st(ix,-1) = fun_us(ic)
        ic = ic + 1
      enddo
      do ix = 0, nx_st-1
        fun_st(ix,ny_st) = fun_us(ic)
        ic = ic + 1
      enddo
      do iy = 0, ny_st-1
        fun_st(-1,iy) = fun_us(ic)
        ic = ic +1
      enddo
      do iy = 0, ny_st-1
        fun_st(nx_st,iy) = fun_us(ic)
        ic = ic +1
      enddo
      
!     pehaps here should be the interpolated values 
      fun_st(-1,-1)       = 0.0_R8     
      fun_st(-1,ny_st)    = 0.0_R8     
      fun_st(nx_st,-1)    = 0.0_R8     
      fun_st(nx_st,ny_st) = 0.0_R8
      
      end subroutine to_struct_cell

      subroutine to_struct_face(nFc,fun_us,funf_st)
      implicit none
      integer :: nFc, ix, iy, ifc, nfx
      real (kind=R8) :: fun_us(nFc), funf_st(-1:nx_st,-1:ny_st,0:1)

      ifc = 1                              ! "poloidal" cell faces case 
      do iy = 0, ny_st-1
        do ix = 0, nx_st
          funf_st(ix,iy,0) = fun_us(ifc)
          ifc = ifc + 1
        enddo
      enddo
      funf_st(-1,:,0)    = 0.0_R8   
      funf_st(:,-1,0)    = 0.0_R8   
      funf_st(:,ny_st,0) = 0.0_R8   
      ifc = 1                             ! "radial" cell faces case
      nfx = (nx_st+1)*ny_st
      do ix = 0, nx_st-1
        do iy = 0, ny_st
          funf_st(ix,iy,1) = fun_us(ifc+nfx)
          ifc = ifc + 1
        enddo
      enddo
      funf_st(:,-1,1)    = 0.0_R8   
      funf_st(-1,:,1)    = 0.0_R8   
      funf_st(nx_st,:,1) = 0.0_R8

      end subroutine to_struct_face

      end module b2us_debug
