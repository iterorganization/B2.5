      subroutine my_outi_us(nwrite, nel, flag, iarr_us, first_name)
cwdk  This routine is intended to write cell-centered fields on an
c     unstructured grid similar to the functionality of the my_out-routine
c     for "original" B2.5 grids.
c     File names will be appended with "*_st.dat".
      
      use b2mod_types
      use b2mod_ad_diff
#ifdef DBG
      use b2us_debug
#endif     
      implicit none
#ifdef DBG
      integer ix, iy
      integer :: iarr_st(-1:nx_st,-1:ny_st), 
     .  iarrf_st(-1:nx_st,-1:ny_st,0:1)
      real (kind=R8) :: rarr(-1:nx_st,-1:ny_st,0:1) 
#endif     
      integer nwrite, nel, flag     ! flag = 0 for cell centers;
      integer :: iarr_us(1:nel)
      character *(*) first_name

      character *64 filename, finame
!      character(len=*), allocatable :: filename, finame
      character(24) hlp_frm
      integer ic
      integer, save :: ncall = 0
      integer, save :: iappend = 0 !, ndigits = 6

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
!        call ipgeti('my_out_digits', ndigits)
!        call xertst(0.le.ndigits, 'invalid ndigits value')
        ncall = 1
*    ..try to create a file in optional output folder
        open(70,file='output/__1__st.dat',err=10)
*    ..folder exists: future output will be done in folder
        my_out_folder='output/'
        write(70,*) 'Test file to check on output folder'
        close(70)
        go to 20
*    ..folder does not exist: future output will be done in current folder
 10     my_out_folder=''
 20     continue
      end if
      
#ifdef DBG
      if (isClassicalGrid.eq.0) then
#endif
!        filename = trim(folder)//trim(first_name)//'_us.dat'
        finame=first_name
        filename = trim(my_out_folder)//trim(finame)//'_us.dat'
        if (iappend.eq.0) then
          open (nwrite,file=trim(filename))
        else
          open (nwrite,file=trim(filename),position='append')
        end if
        write(hlp_frm,"(a13)") '(1x,i7,1x,i7)' 
        do ic = 1, nel
          write (nwrite,hlp_frm) ic, iarr_us(ic)
        end do
        close (nwrite)
#ifdef DBG
      else

      if (flag.eq.1) then

        my_out_sf(0) ='_left_st.dat'
        my_out_sf(1) ='_bottom_st.dat'

        call to_struct_face(nel,real(iarr_us,kind=R8),rarr)
        
        iarrf_st = nint(rarr)
        
        do ic = 0,1  
!          filename=trim(my_out_folder)//trim(first_name)//my_out_sf(ic)
           finame=first_name
           filename = trim(my_out_folder)//trim(finame)//
     &                                      trim(my_out_sf(ic))
!          write(*,*) 'first_name = ',first_name)
!          write(*,*) 'trim(first_name) = ',trim(first_name))
!          write(*,*) 'filename = ',filename
!          write(*,*) 'trim(filename) = ',trim(filename)
          
          if (iappend.eq.0) then
            open (nwrite,file=trim(filename))
          else
            open (nwrite,file=trim(filename),position='append')
          end if
          write(hlp_frm,"(a4,i4,a11)") 
     &       '(6x,',nx_st+2,'(2x,i5))'
!          write(*,*) 'hlp_frm(1) = ', hlp_frm
          write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
          write(hlp_frm,"(a10,i4,a5)")
     &       '(1x,i4,1x,',nx_st+2,'(i7))'
!          write(*,*) 'hlp_frm(2) = ',hlp_frm
          do iy=ny_st,-1,-1
         write(nwrite,hlp_frm) iy, (iarrf_st(ix,iy,ic), ix = -1, nx_st)
          enddo
          close (nwrite)
        enddo
        
!        filename = trim(my_out_folder)//trim(first_name)//'_sum_st.dat'                !srv 13.12.18 {
!        if (iappend.eq.0) then
!          open (nwrite,file=trim(filename))
!        else
!            open (nwrite,file=trim(filename),position='append')
!        end if
!        write(hlp_frm,"(a4,i4,a1,i2,a9)") 
!     &         '(5x,',nx_st+2,'(',ndigits,'x,i3,5x))'
!        write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
!        write(hlp_frm,"(a10,i4,a3,i2.2,a1,i2.2,a2)")
!     &         '(1x,i3,1x,',nx_st+2,'(es',ndigits+8,'.',ndigits,'))'
!        do iy=ny_st,-1,-1
!         write(nwrite,hlp_frm)iy,(funf_st(ix,iy,0)+funf_st(ix,iy,1),
!     &         ix=-1,nx_st)
!        enddo       
!        close (nwrite)
      else

        if(flag.eq.0) then  ! cell centers case    
        
          call to_struct_cell(nel,real(iarr_us,kind=R8),rarr(:,:,0))
          
          iarr_st = nint(rarr(:,:,0))

        finame=first_name

        elseif (flag.eq.2) then

          call to_struct_vertex(nel,real(iarr_us,kind=R8),rarr(:,:,0))

          iarr_st = nint(rarr(:,:,0))

!          filename = trim(my_out_folder)//trim(first_name)//'_st.dat'
          finame=first_name//'_vx'
        
        endif

        filename = trim(my_out_folder)//trim(finame)//'_st.dat'

        if (iappend.eq.0) then
          open (nwrite,file=trim(filename))
        else
          open (nwrite,file=trim(filename),position='append')     !srv 07.09.18
        end if
        write(hlp_frm,"(a4,i4,a11)") 
     &       '(6x,',nx_st+2,'(2x,i5))'
!        write(*,*) 'hlp_frm(1) = ', hlp_frm
        write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
        write(hlp_frm,"(a10,i4,a5)")
     &       '(1x,i4,1x,',nx_st+2,'(i7))'
!        write(*,*) 'hlp_frm(2) = ',hlp_frm
        do iy=ny_st,-1,-1
          write(nwrite,hlp_frm) iy, (iarr_st(ix,iy), ix = -1, nx_st)
        enddo
        close (nwrite)

      endif                                                                     !srv 13.12.18  }

      endif
#endif
      return
      end subroutine my_outi_us
