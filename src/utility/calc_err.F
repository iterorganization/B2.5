      subroutine set_exact_solution(nx, ny, ns)                         !srv 26.02.18 {
       use b2mod_types
	 use b2mod_plasma
	 implicit none
       integer nx, ny, ns

       integer is, ix, iy
	 integer, save :: set_na_numerical = 0, set_ua_numerical = 0,
     *                  set_te_numerical = 0, set_ti_numerical = 0,
     *                  set_po_numerical = 0
	 logical l_exist

	 call ipgeti ('b2mndr_set_na_numerical', set_na_numerical)
	 call ipgeti ('b2mndr_set_ua_numerical', set_ua_numerical)
	 call ipgeti ('b2mndr_set_te_numerical', set_te_numerical)
	 call ipgeti ('b2mndr_set_ti_numerical', set_ti_numerical)
	 call ipgeti ('b2mndr_set_po_numerical', set_po_numerical)

       if(set_na_numerical.ne.0) then 
        inquire(file='../exact_na.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_na.dat')
	 open(70,file='../exact_na.dat')
         read (70,*) (((na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
	 close(70)
        inquire(file='../exact_ne.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_ne.dat')
	 open(70,file='../exact_ne.dat')
         read (70,*) ((ne(ix,iy),ix=-1,nx),iy=-1,ny)
	 close(70)
	 endif

       if(set_ua_numerical.ne.0) then
        inquire(file='../exact_ua.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_ua.dat')
	 open(70,file='../exact_ua.dat')
         read (70,*) (((ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
	 close(70)
	 endif

       if(set_te_numerical.ne.0) then
        inquire(file='../exact_te.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_te.dat')
	 open(70,file='../exact_te.dat')
         read (70,*) ((te(ix,iy),ix=-1,nx),iy=-1,ny)
	 close(70)
	 endif

       if(set_ti_numerical.ne.0) then
        inquire(file='../exact_ti.dat',exist=l_exist)
        call xertst(l_exist,'No ../exact_ti.dat')
	 open(70,file='../exact_ti.dat')
         read (70,*) ((ti(ix,iy),ix=-1,nx),iy=-1,ny)
	 close(70)
	 endif

       if(set_po_numerical.ne.0) then
        inquire(file='../exact_po.dat',exist=l_exist)
        call xertst(l_exist,'No ../exact_po.dat')
	 open(70,file='../exact_po.dat')
         read (70,*) ((po(ix,iy),ix=-1,nx),iy=-1,ny)
	 close(70)
       endif
       
	end subroutine set_exact_solution

	subroutine calc_err(nx, ny, ns)
       use b2mod_types
       use b2mod_constants 
	 use b2mod_plasma
       use b2mod_geo
	 use b2mod_indirect
	 implicit none
       integer nx, ny, ns

       integer is, ix, iy
       real (kind=R8) ::
     &    exact_na(-1:nx,-1:ny,0:ns-1), exact_ua(-1:nx,-1:ny,0:ns-1),
     &    exact_ua_st(-1:nx,-1:ny,0:ns-1),
     &    exact_te(-1:nx,-1:ny), exact_ti(-1:nx,-1:ny),
     &    exact_po(-1:nx,-1:ny), err(-1:nx,-1:ny), era(-1:nx,-1:ny),
     &    uacf(-1:nx,-1:ny)
	 logical l_exist
       character chns*2, chnx*3, chny*3

       write (chnx,'(i3.3)') nx
       write (chny,'(i3.3)') ny
	 open(71, file='output/'//'errors_'//chnx//'x'//chny//'.dat')

       inquire(file='../exact_na.dat',exist=l_exist)
       call xertst(l_exist,'No '//'../exact_na.dat')
	 open(70,file='../exact_na.dat')
        read (70,*) (((exact_na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
	 close(70)

	 do is =0,ns-1
!	  if(.not.is_neutral(is)) then
	   write (chns,'(i2.2)') is
         call my_out(70,nx,ny,na(-1,-1,is),'na'//chns)
         call my_out(70,nx,ny,exact_na(-1,-1,is),'na'//chns//'_precise')
         era = na(:,:,is)-exact_na(:,:,is)
         call my_out(70,nx,ny,era,'na'//chns//'_aerr')
	   write(71,"(' Max|eraNa',a2,'|',8x,'=',es10.3)")
     &              chns,maxval(abs(era))
!	   write(71,"(' Max|eraNa',a2,'|',8x,'=',es10.3,
!     &              ' ix=',i4,' iy=',i4)")
!     &              chns,maxval(abs(era)), maxloc(abs(era))-(/2,2/)
	   if(all(exact_na(:,:,is).ne.0.0_R8)) then
	    err = abs(era)/exact_na(:,:,is)
	    write(71,"(' Max|eraNa',a2,'/Na|     =',es10.3)")
     &          chns,maxval(err)
          call my_out(70,nx,ny,err,'na'//chns//'_rerr')
	   else
	    write(71,"(' Max|eraNa',a2,'/Na|     =',a4)") chns, ' inf'
	   endif
	   write(71,*)
!	  endif
       enddo

       call calc_err_cc(nx, ny, 'Ti', ti, exact_ti, qe) 

       call calc_err_cc(nx, ny, 'Te', te, exact_te, qe) 

       inquire(file='../exact_ua.dat',exist=l_exist)
       call xertst(l_exist,'No '//'../exact_ua.dat')
	 open(70,file='../exact_ua.dat')
        read (70,*) (((exact_ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
	 close(70)

	 do is =0,ns-1
!	  if(.not.is_neutral(is)) then
	   write (chns,'(i2.2)') is
         call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)
         call my_out(70,nx,ny,exact_ua(-1,-1,is),'ua'//chns//'_precise')
         era = ua(:,:,is)-exact_ua(:,:,is)
         call my_out(70,nx,ny,era,'ua'//chns//'_aerr')
	   write(71,"(' Max|eraUa',a2,'|',8x,'=',es10.3)")
     &              chns,maxval(abs(era))
!	   write(71,"(' Max|eraUa',a2,'|',8x,'=',es10.3,
!     &              ' ix=',i4,' iy=',i4)")
!     &              chns,maxval(abs(era)), maxloc(abs(era))-(/2,2/)
	   if(all(exact_ua(:,:,is).ne.0.0_R8)) then
	    err = abs(era/exact_ua(:,:,is))
	    write(71,"(' Max|eraUa',a2,'/ua|     =',es10.3)")
     &          chns,maxval(err)
          call my_out(70,nx,ny,err,'ua'//chns//'_rerr')
	   else
	    write(71,"(' Max|eraUa/Ua|     =',a4)") ' inf'
	   endif
	   write(71,*)
!	  endif
       enddo

!       inquire(file='../exact_ua_st.dat',exist=l_exist)
!       call xertst(l_exist,'No '//'../exact_ua_st.dat')
!	 open(70,file='../exact_ua_st.dat')
!        read (70,*) 
!     &    (((exact_ua_st(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
!	 close(70)

!	 do is =0,ns-1
!!	  if(.not.is_neutral(is)) then
!	   write (chns,'(i2.2)') is
!         call intfaceh(nx,ny,leftix,leftiy,vol,ua(-1,-1,is),uacf)
!         call my_out(70,nx,ny,uacf,'uacf'//chns)
!         call my_out(70,nx,ny,exact_ua_st(-1,-1,is),'uacf'
!     &                                       //chns//'_precise')
!         era = uacf(:,:)-exact_ua_st(:,:,is)
!         call my_out(70,nx,ny,era,'uacf'//chns//'_aerr')
!	   write(71,"(' Max|eraUacf',a2,'|',6x,'=',es10.3)")
!     &              chns,maxval(abs(era))
!!	   write(71,"(' Max|eraUacf',a2,'|',6x,'=',es10.3,
!!     &              ' ix=',i4,' iy=',i4)")
!!     &              chns,maxval(abs(era)), maxloc(abs(era))
!	   if(all(exact_ua_st(:,:,is).ne.0.0_R8)) then
!	    err = abs((era)/exact_ua_st(:,:,is))
!	    write(71,"(' Max|eraUacf',a2,'/uacf  =',es10.3)")
!     &          chns,maxval(err)
!          call my_out(70,nx,ny,err,'uacf'//chns//'_rerr')
!	   else
!	    write(71,"(' Max|eraUacf',a2,'/uacf| =',a4)") chns,'inf'
!	   endif
!!        endif
!       enddo

       call calc_err_cc(nx, ny, 'Po', po, exact_po, 1.0_R8) 

	 close(71)

	end subroutine calc_err

	subroutine calc_err_transport(nx, ny, ns)
       use b2mod_types
       use b2mod_constants 
	 use b2mod_plasma
	 use b2mod_transport
       use b2mod_geo
	 use b2mod_indirect
	 implicit none
       integer nx, ny, ns

       integer is, ix, iy
       real (kind=R8) ::
     &    exact_cdna(-1:nx,-1:ny,0:1,0:ns-1),
     &    exact_cdpa(-1:nx,-1:ny,0:1,0:ns-1),
     &    exact_cddi(-1:nx,-1:ny,0:1,0:ns-1),
!     &    exact_cdnahz(-1:nx,-1:ny,0:1,0:ns-1),
!     &    exact_cdpahz(-1:nx,-1:ny,0:1,0:ns-1),
     &    exact_cvsa(-1:nx,-1:ny,0:1,0:ns-1),
     &    exact_chce(-1:nx,-1:ny,0:1),
     &    exact_chci(-1:nx,-1:ny,0:1),
     &    exact_csig(-1:nx,-1:ny,0:1),
     &    exact_calf(-1:nx,-1:ny,0:1),
!     &    exact_csigin(-1:nx,-1:ny,0:1,0:ns-1,0:ns-1),
     &    err(-1:nx,-1:ny), era(-1:nx,-1:ny)

      end subroutine calc_err_transport

	subroutine calc_err_cc(nx, ny, nmvar, var, exact_var, scl)
       use b2mod_types
	 implicit none
       integer nx, ny
	 character(*) :: nmvar

       real (kind=R8) ::
     &    var(-1:nx,-1:ny), exact_var(-1:nx,-1:ny), scl,
     &    err(-1:nx,-1:ny), era(-1:nx,-1:ny)

       logical l_exist
       integer ix, iy
       character(32) :: pass

	 pass='../exact_'//trim(nmvar)//'.dat'
       inquire(file=trim(pass),exist=l_exist)
       call xertst(l_exist,'No '//trim(pass))

	 open(70,file=pass)
        read (70,*) ((exact_var(ix,iy),ix=-1,nx),iy=-1,ny)
	 close(70)

       call my_out(70,nx,ny,var,trim(nmvar))
       call my_out(70,nx,ny,exact_var,trim(nmvar)//'_precise')

       era = var(:,:)-exact_var(:,:)
       call my_out(70,nx,ny,era,trim(nmvar)//'_aerr')
	 write(71,"(' Max|era',a,'|',10x,'=',es10.3)")
     &              trim(nmvar),maxval(abs(era))/scl
!	 write(71,"(' Max|eravar|',10x,'=',es10.3,' ix=',i4,' iy=',i4)")
!     &              maxval(abs(era)), maxloc(abs(era))-(/2,2/)
	 if(all(var.ne.0.0_R8)) then
	  err = abs(era)/exact_var
	  write(71,"(' Max|era',a,'/',a,'|       =',es10.3)")
     &              trim(nmvar),trim(nmvar),maxval(err)
        call my_out(70,nx,ny,err,trim(nmvar)//'_rerr')
	 else
	   write(71,"(' Max|era',a,'/',a'|     =',a4)")
     &              trim(nmvar), trim(nmvar),' inf'
	 endif
	 write(71,*)

      end subroutine calc_err_cc                                        !srv 26.02.18 }
