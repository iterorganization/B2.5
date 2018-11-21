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

       integer is, iCv, ind                                             !srv 09.07.18
       real (kind=R8) ::
     &    exact_na(nCv,0:ns-1), exact_ua(nCv,0:ns-1),
     &    exact_ua_st(nCv,0:ns-1),
     &    exact_te(nCv), exact_ti(nCv),
     &    exact_po(nCv), err(nCv), era(nCv),
     &    uacf(nCv), L2norm, Cnorm                                      !srv 09.07.18
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
         write (chns,'(i2.2)') is
         call my_out_us(70,nCv,0,pl%na(1,is),'na'//chns)
         call my_out_us(70,nCv,0,exact_na(1,is),'na'//chns//'_precise')
         era = pl%na(:,is)-exact_na(:,is)
         call my_out_us(70,nCv,0,era,'na'//chns//'_aerr')
         ind = maxloc(abs(era),1)                                       !srv 09.07.18
         write(71,"(' Max|eraNa',a2,'|',8x,'=',es10.3,a,i4)")
     &              chns,abs(era(ind)),
     &        '   is in iCv = ', ind
         if(all(exact_na(:,is).ne.0.0_R8)) then
          err = abs(era)/exact_na(:,is)
          ind = maxloc(err,1)
          write(71,"(' Max|eraNa',a2,'/Na|',a2,'|   =',es10.3,a,i4)")
     &          chns, chns, err(ind),
     &      '   is in iCv = ', ind
          call my_out_us(70,nCv,0,err,'na'//chns//'_rerr')
         else
          write(71,"(' Max|eraNa',a2,'/Na',a2,'|   =',a4)") 
     &          chns, chns,' inf'                                       !srv 09.07.18
         endif

         L2norm=norm2(exact_na(:,is))                                   !srv 09.07.18 {
         if(L2norm.ne.0.0_R8) then
          write(71,"(' L2eraNa',a2,'/L2Na',a2,'    =',es10.3)")
     &       chns, chns, norm2(era)/L2norm
          write(71,"(' C_eraNa',a2,'/C_Na',a2,'    =',es10.3)")
     &       chns, chns, maxval(abs(era))/maxval(abs(exact_na(:,is)))
         else
          write(71,"(' L2eraNa',a2,'/L2Na',a2,'    =',a4)") 
     &       chns, chns,' inf'
          write(71,"(' C_eraNa',a2,'/C_Na',a2,'    =',a4)") 
     &       chns, chns,' inf'
         endif                                                           !srv 09.07.18 }
         
         write(71,*)

       enddo

       call calc_err_cc(nCv, 'Ti', pl%ti, exact_ti, qe) 

       call calc_err_cc(nCv, 'Te', pl%te, exact_te, qe) 

       inquire(file='../exact_ua.dat',exist=l_exist)
       call xertst(l_exist,'No '//'../exact_ua.dat')
       open(70,file='../exact_ua.dat')
       read (70,*) ((exact_ua(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       do is =0,ns-1
         write (chns,'(i2.2)') is
         call my_out_us(70,nCv,0,pl%ua(1,is),'ua'//chns)
         call my_out_us(70,nCv,0,exact_ua(1,is),'ua'//chns//'_precise')
         era = pl%ua(:,is)-exact_ua(:,is)
         call my_out_us(70,nCv,0,era,'ua'//chns//'_aerr')
         ind = maxloc(abs(era),1)                                       !srv 09.07.18
         write(71,"(' Max|eraUa',a2,'|',8x,'=',es10.3,a,i4)")
     &              chns,abs(era(ind)),
     &        '   is in iCv = ', ind
         if(all(exact_ua(:,is).ne.0.0_R8)) then
          err = abs(era/exact_ua(:,is))
          ind = maxloc(err,1)
          write(71,"(' Max|eraUa',a2,'/Ua',a2,'|   =',es10.3,a,i4)")
     &          chns, chns, err(ind),                                   !srv 09.07.18
     &      '   is in iCv = ', ind
          call my_out_us(70,nCv,0,err,'ua'//chns//'_rerr')
         else
          write(71,"(' Max|eraUa',a2,'/Ua',a2,'|   =',a4)") 
     &          chns, chns,' inf'                                       !srv 09.07.18
         endif
         
         L2norm=norm2(exact_ua(:,is))                                   !srv 09.07.18 {
         if(L2norm.ne.0.0_R8) then
          write(71,"(' L2eraUa',a2,'/L2Ua',a2,'    =',es10.3)")
     &          chns, chns, norm2(era)/L2norm
          write(71,"(' C_eraUa',a2,'/C_Ua',a2,'    =',es10.3)")
     &       chns, chns, maxval(abs(era))/maxval(abs(exact_ua(:,is)))
         else
          write(71,"(' L2eraUa',a2,'/L2Ua',a2,'    =',a4)") 
     &       chns, chns,' inf'
          write(71,"(' C_eraUa',a2,'/C_Ua',a2,'    =',a4)") 
     &       chns, chns,' inf'
         endif                                                          !srv 09.07.18 }

         write(71,*)

       enddo

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
     &    err(nCv), era(nCv), L2norm, Cnorm                             !srv 09.07.18

       logical l_exist
       integer iCv, ind                                                 !srv 09.07.18
       character(32) :: pass

       pass='../exact_'//trim(nmvar)//'.dat'
       inquire(file=trim(pass),exist=l_exist)
       call xertst(l_exist,'No '//trim(pass))

       open(70,file=pass)
       read (70,*) (exact_var(iCv),iCv=1,nCv)
       close(70)

       call my_out_us(70,nCv,0,var/scl,trim(nmvar))
       call my_out_us(70,nCv,0,exact_var/scl,trim(nmvar)//'_precise')

       era = var-exact_var
       call my_out_us(70,nCv,0,era/scl,trim(nmvar)//'_aerr')
       ind = maxloc(abs(era),1)                                         !srv 09.07.18
       write(71,"(' Max|era',a,'|',10x,'=',es10.3,a,i4)")
     &              trim(nmvar),abs(era(ind))/scl,
     &          '   is in iCv = ', ind
       if(all(var.ne.0.0_R8)) then
        err = abs(era)/exact_var
        ind = maxloc(abs(err),1)                                        !srv 09.07.18
        write(71,"(' Max|era',a,'/',a,'|       =',es10.3,a,i4)")
     &              trim(nmvar),trim(nmvar),err(ind),                   !srv 09.07.18
     &          '   is in iCv = ', ind
        call my_out_us(70,nCv,0,err,trim(nmvar)//'_rerr')
       else
        write(71,"(' Max|era',a,'/',a'|     =',a4)")
     &              trim(nmvar), trim(nmvar),' inf'
       endif

       L2norm=norm2(exact_var)                                          !srv 09.07.18 {
       if(L2norm.ne.0.0_R8) then
        write(71,"(' L2era',a,'/L2',a,'        =',es10.3)")
     &   trim(nmvar),trim(nmvar),norm2(era)/L2norm
        write(71,"(' C_era',a,'/C_',a,'        =',es10.3)")
     &   trim(nmvar),trim(nmvar),maxval(abs(era))/maxval(abs(exact_var))
       else
        write(71,"(' L2era',a,'/L2',a,'        =',a4)") 
     &              trim(nmvar),trim(nmvar),' inf'
        write(71,"(' C_era',a,'/C_',a,'        =',a4)") 
     &              trim(nmvar),trim(nmvar),' inf'
       endif                                                             !srv 09.07.18 }

       write(71,*)

      end subroutine calc_err_cc                                        !srv 26.02.18 }
