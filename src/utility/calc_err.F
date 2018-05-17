      subroutine set_exact_solution(nCv, ns, pl, dv)                         !srv 26.02.18 {
       use b2mod_types
       use b2us_plasma
       implicit none
       integer nCv, ns
       type (B2Plasma), intent (inout) :: pl
       type (B2Derivatives), intent (inout) :: dv

       integer is, iCv
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
        read (70,*) ((pl%na(iCv,is),iCv=1,nCv),is=0,ns-1)
        close(70)
        inquire(file='../exact_ne.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_ne.dat')
        open(70,file='../exact_ne.dat')
        read (70,*) (dv%ne(iCv),iCv=1,nCv)
        close(70)
       endif

       if(set_ua_numerical.ne.0) then
        inquire(file='../exact_ua.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_ua.dat')
        open(70,file='../exact_ua.dat')
        read (70,*) ((pl%ua(iCv,is),iCv=1,nCv),is=0,ns-1)
        close(70)
       endif

       if(set_te_numerical.ne.0) then
        inquire(file='../exact_te.dat',exist=l_exist)
        call xertst(l_exist,'No '//'../exact_te.dat')
        open(70,file='../exact_te.dat')
        read (70,*) (pl%te(iCv),iCv=1,nCv)
        close(70)
       endif

       if(set_ti_numerical.ne.0) then
        inquire(file='../exact_ti.dat',exist=l_exist)
        call xertst(l_exist,'No ../exact_ti.dat')
        open(70,file='../exact_ti.dat')
        read (70,*) (pl%ti(iCv),iCv=1,nCv)
        close(70)
       endif

       if(set_po_numerical.ne.0) then
        inquire(file='../exact_po.dat',exist=l_exist)
        call xertst(l_exist,'No ../exact_po.dat')
        open(70,file='../exact_po.dat')
        read (70,*) (pl%po(iCv),iCv=1,nCv)
        close(70)
       endif
       
      end subroutine set_exact_solution

      subroutine calc_err(nCv, ns, pl, dv)
       use b2mod_types
       use b2mod_constants 
       use b2us_plasma
       implicit none
       integer nCv, ns
       type (B2Plasma), intent (in) :: pl
       type (B2Derivatives), intent (in) :: dv

       integer is, iCv
       real (kind=R8) ::
     &    exact_na(nCv,0:ns-1), exact_ua(nCv,0:ns-1),
     &    exact_ua_st(nCv,0:ns-1),
     &    exact_te(nCv), exact_ti(nCv),
     &    exact_po(nCv), err(nCv), era(nCv),
     &    uacf(nCv)
       logical l_exist
       character chns*2, chnCv*3

       write (chnCv,'(i3.3)') nCv
       open(71, file='output/'//'errors_'//chnCv//'.dat')

       inquire(file='../exact_na.dat',exist=l_exist)
       call xertst(l_exist,'No '//'../exact_na.dat')
       open(70,file='../exact_na.dat')
       read (70,*) ((exact_na(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       do is =0,ns-1
!        if(.not.is_neutral(is)) then
         write (chns,'(i2.2)') is
         call my_out_us(70,nCv,0,pl%na(1,is),'na'//chns)
         call my_out_us(70,nCv,0,exact_na(1,is),'na'//chns//'_precise')
         era = pl%na(:,is)-exact_na(:,is)
         call my_out_us(70,nCv,0,era,'na'//chns//'_aerr')
         write(71,"(' Max|eraNa',a2,'|',8x,'=',es10.3)")
     &              chns,maxval(abs(era))
!	   write(71,"(' Max|eraNa',a2,'|',8x,'=',es10.3,
!     &              ' ix=',i4,' iy=',i4)")
!     &              chns,maxval(abs(era)), maxloc(abs(era))-(/2,2/)
         if(all(exact_na(:,is).ne.0.0_R8)) then
          err = abs(era)/exact_na(:,is)
          write(71,"(' Max|eraNa',a2,'/Na|     =',es10.3)")
     &          chns,maxval(err)
          call my_out_us(70,nCv,0,err,'na'//chns//'_rerr')
         else
          write(71,"(' Max|eraNa',a2,'/Na|     =',a4)") chns, ' inf'
         endif
         write(71,*)
!        endif
       enddo

       call calc_err_cc(nCv, 'Ti', pl%ti, exact_ti, qe) 

       call calc_err_cc(nCv, 'Te', pl%te, exact_te, qe) 

       inquire(file='../exact_ua.dat',exist=l_exist)
       call xertst(l_exist,'No '//'../exact_ua.dat')
       open(70,file='../exact_ua.dat')
       read (70,*) ((exact_ua(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       do is =0,ns-1
!        if(.not.is_neutral(is)) then
         write (chns,'(i2.2)') is
         call my_out_us(70,nCv,0,pl%ua(1,is),'ua'//chns)
         call my_out_us(70,nCv,0,exact_ua(1,is),'ua'//chns//'_precise')
         era = pl%ua(:,is)-exact_ua(:,is)
         call my_out_us(70,nCv,0,era,'ua'//chns//'_aerr')
         write(71,"(' Max|eraUa',a2,'|',8x,'=',es10.3)")
     &              chns,maxval(abs(era))
!	   write(71,"(' Max|eraUa',a2,'|',8x,'=',es10.3,
!     &              ' ix=',i4,' iy=',i4)")
!     &              chns,maxval(abs(era)), maxloc(abs(era))-(/2,2/)
         if(all(exact_ua(:,is).ne.0.0_R8)) then
          err = abs(era/exact_ua(:,is))
          write(71,"(' Max|eraUa',a2,'/ua|     =',es10.3)")
     &          chns,maxval(err)
          call my_out_us(70,nCv,0,err,'ua'//chns//'_rerr')
         else
          write(71,"(' Max|eraUa/Ua|     =',a4)") ' inf'
         endif
         write(71,*)
!        endif
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

       call calc_err_cc(nCv, 'Po', pl%po, exact_po, 1.0_R8) 

       close(71)

      end subroutine calc_err

      subroutine calc_err_transport(nFc, ns)
       use b2mod_types
       use b2mod_constants 
       use b2us_plasma
       implicit none
       integer nFc, ns

       integer is, iFc
       real (kind=R8) ::
     &    exact_cdna(nFc,0:1,0:ns-1),
     &    exact_cdpa(nFc,0:1,0:ns-1),
     &    exact_cddi(nFc,0:1,0:ns-1),
!     &    exact_cdnahz(nFc,0:1,0:ns-1),
!     &    exact_cdpahz(nFc,0:1,0:ns-1),
     &    exact_cvsa(nFc,0:1,0:ns-1),
     &    exact_chce(nFc,0:1),
     &    exact_chci(nFc,0:1),
     &    exact_csig(nFc,0:1),
     &    exact_calf(nFc,0:1),
!     &    exact_csigin(nFc,0:1,0:ns-1,0:ns-1),
     &    err(nFc), era(nFc)

      end subroutine calc_err_transport

      subroutine calc_err_cc(nCv, nmvar, var, exact_var, scl)
       use b2mod_types
       implicit none
       integer nCv
       character(*) :: nmvar

       real (kind=R8) ::
     &    var(nCv), exact_var(nCv), scl,
     &    err(nCv), era(nCv)

       logical l_exist
       integer iCv
       character(32) :: pass

       pass='../exact_'//trim(nmvar)//'.dat'
       inquire(file=trim(pass),exist=l_exist)
       call xertst(l_exist,'No '//trim(pass))

       open(70,file=pass)
       read (70,*) (exact_var(iCv),iCv=1,nCv)
       close(70)

       call my_out_us(70,nCv,0,var,trim(nmvar))
       call my_out_us(70,nCv,0,exact_var,trim(nmvar)//'_precise')

       era = var-exact_var
       call my_out_us(70,nCv,0,era,trim(nmvar)//'_aerr')
       write(71,"(' Max|era',a,'|',10x,'=',es10.3)")
     &              trim(nmvar),maxval(abs(era))/scl
!       write(71,"(' Max|eravar|',10x,'=',es10.3,' ix=',i4,' iy=',i4)")
!     &              maxval(abs(era)), maxloc(abs(era))-(/2,2/)
       if(all(var.ne.0.0_R8)) then
        err = abs(era)/exact_var
        write(71,"(' Max|era',a,'/',a,'|       =',es10.3)")
     &              trim(nmvar),trim(nmvar),maxval(err)
        call my_out_us(70,nCv,0,err,trim(nmvar)//'_rerr')
       else
        write(71,"(' Max|era',a,'/',a'|     =',a4)")
     &              trim(nmvar), trim(nmvar),' inf'
       endif
       write(71,*)

      end subroutine calc_err_cc                                        !srv 26.02.18 }
