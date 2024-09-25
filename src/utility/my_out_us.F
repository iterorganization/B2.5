      subroutine my_out_us(nwrite, nel, flag, fun_us, first_name)
cwdk  This routine is intended to write cell-centered fields on an
c     unstructured grid similar to the functionality of the my_out-routine
c     for "original" B2.5 grids.
c     File names will be appended with "*_st.dat".

      use b2mod_types
      use b2mod_ad
#ifdef WG_DEBUG
      use b2us_debug
#endif
      implicit none
#ifdef WG_DEBUG
      integer ix, iy
      real (kind=R8) :: fun_st(-1:nx_st,-1:ny_st),
     .  funf_st(-1:nx_st,-1:ny_st,0:1)
#endif
      integer nwrite, nel, flag     ! flag = 0 for cell centers;
      real (kind=R8) :: fun_us(1:nel)
      character *(*) first_name

      character *64 filename, finame
!      character(len=*), allocatable :: filename, finame
      character(26) hlp_frm
      integer ic
      integer, save :: ncall = 0
      integer, save :: iappend = 0, ndigits = 6
      external ipgeti, xertst

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        call ipgeti('my_out_digits', ndigits)
        call xertst(0.le.ndigits, 'invalid ndigits value')
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

#ifdef WG_DEBUG
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
        write(hlp_frm,"(a,i2,a,i2.2,a)")
     .        '(1x,i9,1x,es',ndigits+9,'.',ndigits,'e3)'
        do ic = 1, nel
          write (nwrite,hlp_frm) ic, fun_us(ic)
        end do
        close (nwrite)
#ifdef WG_DEBUG
      else

      if (flag.eq.1) then

        my_out_sf(0) ='_left_st.dat'
        my_out_sf(1) ='_bottom_st.dat'

        call to_struct_face(nel,fun_us,funf_st)

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
          write(hlp_frm,"(a4,i4,a1,i2,a5,i2,a3)")
     &         '(5x,',nx_st+2,
     &         '(',6+ndigits-(ndigits+6)/2,'x,i3,',(ndigits+6)/2,'x))'
          write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
          write(hlp_frm,"(a10,i4,a3,i2.2,a1,i2.2,a4)")
     &         '(1x,i3,1x,',nx_st+2,'(es',ndigits+9,'.',ndigits,'e3))'
          do iy=ny_st,-1,-1
           write(nwrite,hlp_frm) iy, (funf_st(ix,iy,ic), ix = -1, nx_st)
          enddo
          close (nwrite)
        enddo

!        filename = trim(my_out_folder)//trim(first_name)//'_sum_st.dat' !srv 13.12.18 {
!        if (iappend.eq.0) then
!          open (nwrite,file=trim(filename))
!        else
!            open (nwrite,file=trim(filename),position='append')
!        end if
!        write(hlp_frm,"(a4,i4,a1,i2,a5,i2,a3)")
!     &         '(5x,',nx_st+2,
!     &         '(',6+ndigits-(ndigits+6)/2,'x,i3,',(ndigits+6)/2,'x))'
!        write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
!        write(hlp_frm,"(a10,i4,a3,i2.2,a1,i2.2,a4)")
!     &         '(1x,i3,1x,',nx_st+2,'(es',ndigits+9,'.',ndigits,'e3))'
!        do iy=ny_st,-1,-1
!         write(nwrite,hlp_frm)iy,(funf_st(ix,iy,0)+funf_st(ix,iy,1),
!     &         ix=-1,nx_st)
!        enddo
!        close (nwrite)
      else

        if(flag.eq.0) then  ! cell centers case

          call to_struct_cell(nel,fun_us,fun_st)

        finame=first_name

        elseif (flag.eq.2) then

          call to_struct_vertex(nel,fun_us,fun_st)

!          filename = trim(my_out_folder)//trim(first_name)//'_st.dat'
          finame=first_name//'_vx'

        endif

        filename = trim(my_out_folder)//trim(finame)//'_st.dat'

        if (iappend.eq.0) then
          open (nwrite,file=trim(filename))
        else
          open (nwrite,file=trim(filename),position='append')     !srv 07.09.18
        end if
        write(hlp_frm,"(a4,i4,a1,i2,a5,i2,a3)")
     &       '(5x,',nx_st+2,
     &       '(',6+ndigits-(ndigits+6)/2,'x,i3,',(ndigits+6)/2,'x))'
        write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
        write(hlp_frm,"(a10,i4,a3,i2.2,a1,i2.2,a4)")
     &       '(1x,i3,1x,',nx_st+2,'(es',ndigits+9,'.',ndigits,'e3))'
        do iy=ny_st,-1,-1
          write(nwrite,hlp_frm) iy, (fun_st(ix,iy), ix = -1, nx_st)
        enddo
        close (nwrite)

      endif                                                       !srv 13.12.18  }

      endif
#endif
      return
      end subroutine my_out_us
