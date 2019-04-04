      module b2us_debug

      use b2mod_types

      implicit none

      private

      public :: to_struct_cell, to_struct_face

      integer, public :: nx_st, ny_st, isClassicalGrid
      integer, allocatable, public :: imapCv(:,:),
     .                                imapFcx(:,:),
     .                                imapFcy(:,:),
     .                                imapVx(:,:)

      contains

      subroutine to_struct_cell(nCv, fun_us, fun_st)
      implicit none
      integer :: nCv, ix, iy
      real (kind=R8) :: fun_us(nCv), fun_st(-1:nx_st,-1:ny_st)

      do iy = -1, ny_st
        do ix = -1, nx_st
          if(imapCv(ix,iy) /= 0) then          
            fun_st(ix,iy) = fun_us(imapCv(ix,iy))
          else
            ! perhaps here should be the interpolated values 
            fun_st(ix,iy) = 0.0_R8  
          endif
        enddo
      enddo
      
      end subroutine to_struct_cell

      subroutine to_struct_face(nFc,fun_us,funf_st)
      implicit none
      integer :: nFc, ix, iy
      real (kind=R8) :: fun_us(nFc), funf_st(-1:nx_st,-1:ny_st,0:1)

      do iy = -1, ny_st
        do ix = -1, nx_st
          if(imapFcx(ix,iy) /= 0) then ! poloidal faces
            funf_st(ix,iy,0) = fun_us(imapFcx(ix,iy))
          else
            ! perhaps here should be the interpolated values 
            funf_st(ix,iy,0) = 0.0_R8  
          endif
          if(imapFcy(ix,iy) /= 0) then ! radial faces
            funf_st(ix,iy,1) = fun_us(imapFcy(ix,iy))
          else
            ! perhaps here should be the interpolated values 
            funf_st(ix,iy,1) = 0.0_R8  
          endif
        enddo
      enddo

      end subroutine to_struct_face

      end module b2us_debug
