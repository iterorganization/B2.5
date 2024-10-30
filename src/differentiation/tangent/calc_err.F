      subroutine set_exact_solution(nCv, ns, switch, pl, dv)             !srv 26.02.18 {
       use b2mod_types
       use b2us_plasma_diffv
       use b2mod_switches_diffv
       implicit none
       integer nCv, ns
       type (switches), intent (in) :: switch
       type (B2Plasma), intent (inout) :: pl
       type (B2Derivatives), intent (inout) :: dv

       integer is, iCv
       character*27 filename
       logical l_exist

       if(switch%set_na_numerical.ne.0) then
        filename='exact_na_us.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_na_us.dat file')
        open(70,file=trim(filename))
        read (70,*) ((pl%na(iCv,is),iCv=1,nCv),is=0,ns-1)
        close(70)
        filename='exact_ne_us.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ne_us.dat file')
        open(70,file=trim(filename))
        read (70,*) (dv%ne(iCv),iCv=1,nCv)
        close(70)
       endif

       if(switch%set_ua_numerical.ne.0) then
        filename='exact_ua_us.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ua_us.dat file')
        open(70,file=trim(filename))
        read (70,*) ((pl%ua(iCv,is),iCv=1,nCv),is=0,ns-1)
        close(70)
       endif

       if(switch%set_te_numerical.ne.0) then
        filename='exact_te_us.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_te_us.dat file')
        open(70,file=trim(filename))
        read (70,*) (pl%te(iCv),iCv=1,nCv)
        close(70)
       endif

       if(switch%set_ti_numerical.ne.0) then
        filename='exact_ti_us,dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ti_us.dat file')
        open(70,file=trim(filename))
        read (70,*) (pl%ti(iCv),iCv=1,nCv)
        close(70)
       endif

       if(switch%set_po_numerical.ne.0) then
        filename='exact_po_us.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_po_us.dat file')
        open(70,file=trim(filename))
        read (70,*) (pl%po(iCv),iCv=1,nCv)
        close(70)
       endif

       return
       end subroutine set_exact_solution

       subroutine calc_err(nCv, ns, pl)
       use b2mod_types
       use b2mod_constants
       use b2us_plasma_diffv
#ifdef LEGACYCOMP
       use b2mod_math_diffv
#endif
       implicit none
       integer nCv, ns
       type (B2Plasma), intent (in) :: pl

       integer is, iCv, ind                                             !srv 09.07.18
       real (kind=R8) ::
     &    exact_na(nCv,0:ns-1), exact_ua(nCv,0:ns-1),
     &    exact_te(nCv), exact_ti(nCv),
     &    exact_po(nCv), err(nCv), era(nCv),
     &    L2norm                                                        !srv 09.07.18
       logical l_exist
       character chns*3, chnCv*6
       character*27 filename
#ifndef LEGACYCOMP
       intrinsic norm2
#endif

       write (chnCv,'(i6.6)') nCv
       open(71, file='output/'//'errors_'//chnCv//'_us.dat')

       filename='exact_na_us.dat'
       call find_file(filename,l_exist)
       call xertst(l_exist,'No exact_na_us.dat file')
       open(70,file=trim(filename))
       read (70,*) ((exact_na(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       do is = 0, ns-1
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,0,pl%na(1,is),'na'//chns)
         call my_out_us(70,nCv,0,exact_na(1,is),'na'//chns//'_precise')
         era = pl%na(:,is)-exact_na(:,is)
         call my_out_us(70,nCv,0,era,'na'//chns//'_aerr')
         ind = maxloc(abs(era),1)                                       !srv 09.07.18
         write(71,'(a,a3,a,8x,a,es10.3,a,i7)')
     &        ' Max|eraNa',chns,'|','=',abs(era(ind)),
     &        '   is in iCv = ', ind
         if(all(exact_na(:,is).ne.0.0_R8)) then
           err = abs(era)/exact_na(:,is)
           ind = maxloc(err,1)
           write(71,'(a,a3,a,a3,a,es10.3,a,i7)')
     &      ' Max|eraNa',chns,'/Na',chns,'|   =',err(ind),              !srv 09.07.18
     &      '   is in iCv = ', ind
           call my_out_us(70,nCv,0,err,'na'//chns//'_rerr')
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' Max|eraNa',chns,'/Na',chns,'|   =',' Inf'                 !srv 09.07.18
         endif

#ifndef LEGACYCOMP
         L2norm=norm2(exact_na(:,is))                                   !srv 09.07.18 {
#else
         L2norm=norm2(nCv,exact_na(:,is))
#endif
         if(L2norm.ne.0.0_R8) then
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' L2eraNa',chns,'/L2Na',chns,'    =',
#ifndef LEGACYCOMP
     &        norm2(era)/L2norm
#else
     &        norm2(nCv,era)/L2norm
#endif
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' C_eraNa',chns,'/C_Na',chns,'    =',
     &        maxval(abs(era))/maxval(abs(exact_na(:,is)))
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' L2eraNa',chns,'/L2Na',chns,'    =',' Inf'
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' C_eraNa',chns,'/C_Na',chns,'    =',' Inf'
         endif                                                          !srv 09.07.18 }

         write(71,*)
       enddo

       call calc_err_cc(nCv, 'Te', pl%te, exact_te, qe)                 !srv 29.12.18

       call calc_err_cc(nCv, 'Ti', pl%ti, exact_ti, qe) 


       filename='exact_ua_us.dat'
       call find_file(filename,l_exist)
       call xertst(l_exist,'No exact_ua_us.dat file')
       open(70,file=trim(filename))
       read (70,*) ((exact_ua(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       do is = 0, ns-1
         write (chns,'(i3.3)') is
         call my_out_us(70,nCv,0,pl%ua(1,is),'ua'//chns)
         call my_out_us(70,nCv,0,exact_ua(1,is),'ua'//chns//'_precise')
         era = pl%ua(:,is)-exact_ua(:,is)
         call my_out_us(70,nCv,0,era,'ua'//chns//'_aerr')
         ind = maxloc(abs(era),1)                                       !srv 09.07.18
         write(71,'(a,a3,a,8x,a,es10.3,a,i7)')
     &    ' Max|eraUa',chns,'|','=',abs(era(ind)),
     &    '   is in iCv = ', ind
         if(all(exact_ua(:,is).ne.0.0_R8)) then
           err = abs(era/exact_ua(:,is))
           ind = maxloc(err,1)
           write(71,'(a,a3,a,a3,a,es10.3,a,i7)')
     &      ' Max|eraUa',chns,'/Ua',chns,'|   =',err(ind),              !srv 09.07.18
     &      '   is in iCv = ', ind
           call my_out_us(70,nCv,0,err,'ua'//chns//'_rerr')
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' Max|eraUa',chns,'/Ua',chns,'|   =',' Inf'                 !srv 09.07.18
         endif

#ifndef LEGACYCOMP
         L2norm=norm2(exact_ua(:,is))                                   !srv 09.07.18 {
#else
         L2norm=norm2(nCv,exact_ua(:,is))
#endif
         if(L2norm.ne.0.0_R8) then
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' L2eraUa',chns,'/L2Ua',chns,'    =',
#ifndef LEGACYCOMP
     &        norm2(era)/L2norm
#else
     &        norm2(nCv,era)/L2norm
#endif
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' C_eraUa',chns,'/C_Ua',chns,'    =',
     &        maxval(abs(era))/maxval(abs(exact_ua(:,is)))
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' L2eraUa',chns,'/L2Ua',chns,'    =',' Inf'
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' C_eraUa',chns,'/C_Ua',chns,'    =',' Inf'
         endif                                                          !srv 09.07.18 }

         write(71,*)
       enddo

       call calc_err_cc(nCv, 'Po', pl%po, exact_po, 1.0_R8)

       close(71)
       return
       end subroutine calc_err

       subroutine calc_err_transport(nFc, ns)
       use b2mod_types
       use b2mod_constants
       use b2us_plasma_diffv
       implicit none
       integer nFc, ns

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

      return
      end subroutine calc_err_transport

      subroutine calc_err_cc(nCv, nmvar, var, exact_var, scl)
       use b2mod_types
#ifdef LEGACYCOMP
       use b2mod_math_diffv
#endif
       implicit none
       integer nCv
       character(*) :: nmvar

       real (kind=R8) ::
     &    var(nCv), exact_var(nCv), scl,
     &    err(nCv), era(nCv), L2norm                                    !srv 09.07.18
       logical l_exist
       integer iCv, ind                                                 !srv 09.07.18
       character(32) :: pass
#ifndef LEGACYCOMP
       intrinsic norm2
#endif

       pass='exact_'//trim(nmvar)//'_us.dat'
       call find_file(pass,l_exist)
       call xertst(l_exist,'Did not find '//trim(pass)//' file')

       open(70,file=pass)
       read (70,*) (exact_var(iCv),iCv=1,nCv)
       close(70)

       call my_out_us(70,nCv,0,var/scl,trim(nmvar))
       call my_out_us(70,nCv,0,exact_var/scl,trim(nmvar)//'_precise')

       era = var-exact_var
       call my_out_us(70,nCv,0,era/scl,trim(nmvar)//'_aerr')
       ind = maxloc(abs(era),1)                                         !srv 09.07.18
       write(71,'(a,a,a,10x,a,es10.3,a,i7)')
     &  ' Max|era',trim(nmvar),'|','=',abs(era(ind))/scl,
     &  '   is in iCv = ', ind
       if(all(var.ne.0.0_R8)) then
        err = abs(era)/exact_var
        ind = maxloc(abs(err),1)                                        !srv 09.07.18
        write(71,'(a,a,a,a,a,es10.3,a,i7)')
     &  ' Max|era',trim(nmvar),'/',trim(nmvar),'|       =',
     &   err(ind),
     &  '   is in iCv = ', ind
        call my_out_us(70,nCv,0,err,trim(nmvar)//'_rerr')
       else
        write(71,'(a,a,a,a,a,a4)')
     &   ' Max|era',trim(nmvar),'/',trim(nmvar),'|     =',' Inf'
       endif

#ifndef LEGACYCOMP
       L2norm=norm2(exact_var)                                          !srv 09.07.18 {
#else
       L2norm=norm2(nCv,exact_var)
#endif
       if(L2norm.ne.0.0_R8) then
        write(71,'(a,a,a,a,a,es10.3)')
     &   ' L2era',trim(nmvar),'/L2',trim(nmvar),'        =',
#ifndef LEGACYCOMP
     &   norm2(era)/L2norm
#else
     &   norm2(nCv,era)/L2norm
#endif
        write(71,'(a,a,a,a,a,es10.3)')
     &   ' C_era',trim(nmvar),'/C_',trim(nmvar),'        =',
     &   maxval(abs(era))/maxval(abs(exact_var))
       else
        write(71,'(a,a,a,a,a,a4)')
     &   ' L2era',trim(nmvar),'/L2',trim(nmvar),'        =',' Inf'
        write(71,'(a,a,a,a,a,a4)')
     &   ' C_era',trim(nmvar),'/C_',trim(nmvar),'        =',' Inf'
       endif                                                            !srv 09.07.18 }

       write(71,*)

      return
      end subroutine calc_err_cc                                        !srv 26.02.18 }
