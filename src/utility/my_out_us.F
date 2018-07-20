      subroutine my_out_us(nwrite, nCv, flag, fun_us, first_name)
cwdk  This routine is intended to write cell-centered fields on an
c     unstructured grid similar to the functionality of the my_out-routine
c     for "original" B2.5 grids.
c     File names will be appended with "*_st.dat".
      
!     This subroutine sould be modified for double null case

      use b2mod_types
#ifdef WG_DEBUG
      use b2us_debug
#endif     
      implicit none
#ifdef WG_DEBUG
      real (kind=R8) :: fun_st(-1:nx_st,-1:ny_st), 
     .  funf_st(-1:nx_st,-1:ny_st,0:1)
      character(8), save :: sf(0:1)
#endif     
      integer nwrite, nCv, flag     ! flag = 0 for cell centers;
      real (kind=R8) :: fun_us(1:nCv)
      character *(*) first_name

      character *64 filename
      character(23) hlp_frm
      integer ic, ix, iy
      integer, save :: ncall = 0
      integer, save :: iappend = 0
      logical ex
      character(7), save :: folder

      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
        ncall = 1
*    ..try to create a file in optional output folder
        open(70,file='output/__1__st.dat',err=10)
*    ..folder exists: future output will be done in folder
        folder='output/'
        write(70,*) 'Test file to check on output folder'
        close(70)
        go to 20
*    ..folder does not exist: future output will be done in current folder
 10     folder=''
 20     continue
      end if
      
#ifndef WG_DEBUG         
        filename = trim(folder)//trim(first_name)//'_us.dat'
        if (iappend.eq.0) then
          open (nwrite,file=filename)
        else
          open (nwrite,file=filename,position='append')
        end if
!       write(hlp_frm,"(a4,i4,a11)") '(5x,',nCv,'(6x,i3,5x))'
!       write (nwrite,hlp_frm) (ic,ic=1,nCv)
!       write (hlp_frm,"(a10,i4,a9)") '(5x,',nCv,'(es14.6))'
!       write (nwrite,hlp_frm) (fun(ic), ic = 1, nCv)
        do ic = 1, nCv
          write (nwrite,'(1x,i9,1x,es14.6)') ic, fun_us(ic)
        end do
        close (nwrite)
#else

      if(flag.eq.0) then  ! cell centers case    
        
        call to_struct_cell(nCv,fun_us,fun_st)

        filename = trim(folder)//trim(first_name)//'_st.dat'
        if (iappend.eq.0) then
          open (nwrite,file=filename)
        else
          open (nwrite,file=filename,position='append')
        end if
        write(hlp_frm,"(a4,i4,a11)") '(5x,',nx_st+2,'(6x,i3,5x))'
        write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
        write(hlp_frm,"(a10,i4,a9)") '(1x,i3,1x,',nx_st+2,'(es14.6))'
        do iy=ny_st,-1,-1
          write(nwrite,hlp_frm) iy, (fun_st(ix,iy), ix = -1, nx_st)
        enddo
        close (nwrite)
        
      elseif(flag.eq.1) then
        sf(0) ='x_st.dat'
        sf(1) ='y_st.dat'

        call to_struct_face(nFc,fun_us,funf_st)
        
        do ic = 0,1  
          filename = trim(folder)//trim(first_name)//sf(ic)
          if (iappend.eq.0) then
            open (nwrite,file=filename)
          else
            open (nwrite,file=filename,position='append')
          end if
          write(hlp_frm,"(a4,i4,a11)") '(5x,',nx_st+2,'(6x,i3,5x))'
          write (nwrite,hlp_frm) (ix,ix=-1,nx_st)
          write(hlp_frm,"(a10,i4,a9)") '(1x,i3,1x,',nx_st+2,'(es14.6))'
          do iy=ny_st,-1,-1
            write(nwrite,hlp_frm) iy, (funf_st(ix,iy,ic), ix = -1, nx_st)
          enddo
          close (nwrite)
        enddo
      endif
#endif
      return
      end
