      subroutine set_exact_solution(nx, ny, ns)                         !srv 26.02.18 {
       use b2mod_types
       use b2mod_plasma
       implicit none
       integer nx, ny, ns

       integer is, ix, iy
       integer, save :: set_na_numerical = 0, set_ua_numerical = 0,
     *                  set_te_numerical = 0, set_ti_numerical = 0,
     *                  set_po_numerical = 0
       character*24 filename
       logical l_exist

       call ipgeti ('b2mndr_set_na_numerical', set_na_numerical)
       call ipgeti ('b2mndr_set_ua_numerical', set_ua_numerical)
       call ipgeti ('b2mndr_set_te_numerical', set_te_numerical)
       call ipgeti ('b2mndr_set_ti_numerical', set_ti_numerical)
       call ipgeti ('b2mndr_set_po_numerical', set_po_numerical)

       if(set_na_numerical.ne.0) then
        filename='exact_na.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_na.dat file')
        open(70,file=trim(filename))
        read (70,*) (((na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        close(70)
        filename='exact_ne.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ne.dat file')
        open(70,file=trim(filename))
        read (70,*) ((ne(ix,iy),ix=-1,nx),iy=-1,ny)
        close(70)
       endif

       if(set_ua_numerical.ne.0) then
        filename='exact_ua.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ua.dat file')
        open(70,file=trim(filename))
        read (70,*) (((ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
        close(70)
       endif

       if(set_te_numerical.ne.0) then
        filename='exact_te.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_te.dat file')
        open(70,file=trim(filename))
        read (70,*) ((te(ix,iy),ix=-1,nx),iy=-1,ny)
        close(70)
       endif

       if(set_ti_numerical.ne.0) then
        filename='exact_ti,dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_ti.dat file')
        open(70,file=trim(filename))
        read (70,*) ((ti(ix,iy),ix=-1,nx),iy=-1,ny)
        close(70)
       endif

       if(set_po_numerical.ne.0) then
        filename='exact_po.dat'
        call find_file(filename,l_exist)
        call xertst(l_exist,'No exact_po.dat file')
        open(70,file=trim(filename))
        read (70,*) ((po(ix,iy),ix=-1,nx),iy=-1,ny)
        close(70)
       endif

       return
       end subroutine set_exact_solution

       subroutine calc_err(nx, ny, ns)
       use b2mod_types
       use b2mod_constants
       use b2mod_plasma
       use b2mod_geo
       use b2mod_indirect
       implicit none
       integer nx, ny, ns

       integer is, ix, iy, ind(2)                                        !srv 09.07.18
       real (kind=R8) ::
     &    exact_na(-1:nx,-1:ny,0:ns-1), exact_ua(-1:nx,-1:ny,0:ns-1),
     &    exact_te(-1:nx,-1:ny), exact_ti(-1:nx,-1:ny),
     &    exact_po(-1:nx,-1:ny), err(-1:nx,-1:ny), era(-1:nx,-1:ny),
     &    L2norm                                                         !srv 09.07.18
       logical l_exist
       character chns*3, chnx*4, chny*4
       character*24 filename

       write (chnx,'(i4.4)') nx
       write (chny,'(i4.4)') ny
       open(71, file='output/'//'errors_'//chnx//'x'//chny//'.dat')

       filename='exact_na.dat'
       call find_file(filename,l_exist)
       call xertst(l_exist,'No exact_na.dat file')
       open(70,file=trim(filename))
       read (70,*) (((exact_na(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
       close(70)

       do is = 0, ns-1
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,na(-1,-1,is),'na'//chns)
         call my_out(70,nx,ny,exact_na(-1,-1,is),'na'//chns//'_precise')
         era = na(:,:,is)-exact_na(:,:,is)
         call my_out(70,nx,ny,era,'na'//chns//'_aerr')
         ind = maxloc(abs(era))                                          !srv 09.07.18
         write(71,'(a,a3,a,8x,a,es10.3,a,i4,a,i4)')
     &        ' Max|eraNa',chns,'|','=',abs(era(ind(1)-2,ind(2)-2)),
     &        '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
         if(all(exact_na(:,:,is).ne.0.0_R8)) then
           err = abs(era)/exact_na(:,:,is)
           ind = maxloc(err)
           write(71,'(a,a3,a,a3,a,es10.3,a,i4,a,i4)')
     &      ' Max|eraNa',chns,'/Na',chns,'|   =',err(ind(1)-2,ind(2)-2), !srv 09.07.18
     &      '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
           call my_out(70,nx,ny,err,'na'//chns//'_rerr')
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' Max|eraNa',chns,'/Na',chns,'|   =',' Inf'                  !srv 09.07.18
         endif

         L2norm=norm2(exact_na(:,:,is))                                  !srv 09.07.18 {
         if(L2norm.ne.0.0_R8) then
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' L2eraNa',chns,'/L2Na',chns,'    =',norm2(era)/L2norm
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' C_eraNa',chns,'/C_Na',chns,'    =',
     &        maxval(abs(era))/maxval(abs(exact_na(:,:,is)))
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' L2eraNa',chns,'/L2Na',chns,'    =',' Inf'
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' C_eraNa',chns,'/C_Na',chns,'    =',' Inf'
         endif                                                          !srv 09.07.18 }

         write(71,*)
       enddo

       call calc_err_cc(nx, ny, 'Te', te, exact_te, qe)                 !srv 09.07.18

       call calc_err_cc(nx, ny, 'Ti', ti, exact_ti, qe)                 !srv 09.07.18

       filename='exact_ua.dat'
       call find_file(filename,l_exist)
       call xertst(l_exist,'No exact_ua.dat file')
       open(70,file=trim(filename))
       read (70,*) (((exact_ua(ix,iy,is),ix=-1,nx),iy=-1,ny),is=0,ns-1)
       close(70)

       do is = 0, ns-1
         write (chns,'(i3.3)') is
         call my_out(70,nx,ny,ua(-1,-1,is),'ua'//chns)
         call my_out(70,nx,ny,exact_ua(-1,-1,is),'ua'//chns//'_precise')
         era = ua(:,:,is)-exact_ua(:,:,is)
         call my_out(70,nx,ny,era,'ua'//chns//'_aerr')
         ind = maxloc(abs(era))                                           !srv 09.07.18
         write(71,'(a,a3,a,8x,a,es10.3,a,i4,a,i4)')
     &    ' Max|eraUa',chns,'|','=',abs(era(ind(1)-2,ind(2)-2)),
     &    '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
         if(all(exact_ua(:,:,is).ne.0.0_R8)) then
           err = abs(era/exact_ua(:,:,is))
           ind = maxloc(err)
           write(71,'(a,a3,a,a3,a,es10.3,a,i4,a,i4)')
     &      ' Max|eraUa',chns,'/Ua',chns,'|   =',err(ind(1)-2,ind(2)-2),  !srv 09.07.18
     &      '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
           call my_out(70,nx,ny,err,'ua'//chns//'_rerr')
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' Max|eraUa',chns,'/Ua',chns,'|   =',' Inf'                   !srv 09.07.18
         endif

         L2norm=norm2(exact_ua(:,:,is))                                   !srv 09.07.18 {
         if(L2norm.ne.0.0_R8) then
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' L2eraUa',chns,'/L2Ua',chns,'    =',norm2(era)/L2norm
           write(71,'(a,a3,a,a3,a,es10.3)')
     &      ' C_eraUa',chns,'/C_Ua',chns,'    =',
     &        maxval(abs(era))/maxval(abs(exact_ua(:,:,is)))
         else
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' L2eraUa',chns,'/L2Ua',chns,'    =',' Inf'
           write(71,'(a,a3,a,a3,a,a4)')
     &      ' C_eraUa',chns,'/C_Ua',chns,'    =',' Inf'
         endif                                                           !srv 09.07.18 }

         write(71,*)
       enddo

       call calc_err_cc(nx, ny, 'Po', po, exact_po, 1.0_R8)

       close(71)
       return

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

      end subroutine calc_err_transport

      subroutine calc_err_cc(nx, ny, nmvar, var, exact_var, scl)
       use b2mod_types
       implicit none
       integer nx, ny
       character(*) :: nmvar

       real (kind=R8) ::
     &    var(-1:nx,-1:ny), exact_var(-1:nx,-1:ny), scl,
     &    err(-1:nx,-1:ny), era(-1:nx,-1:ny), L2norm                   !srv 09.07.18

       logical l_exist
       integer ix, iy, ind(2)                                          !srv 09.07.18
       character(32) :: pass

       pass='exact_'//trim(nmvar)//'.dat'
       call find_file(pass,l_exist)
       call xertst(l_exist,'Did not find '//trim(pass)//' file')

       open(70,file=pass)
       read (70,*) ((exact_var(ix,iy),ix=-1,nx),iy=-1,ny)
       close(70)

       call my_out(70,nx,ny,var/scl,trim(nmvar))
       call my_out(70,nx,ny,exact_var/scl,trim(nmvar)//'_precise')

       era = var(:,:)-exact_var(:,:)
       call my_out(70,nx,ny,era/scl,trim(nmvar)//'_aerr')
       ind = maxloc(abs(era))                                          !srv 09.07.18
       write(71,'(a,a,a,10x,a,es10.3,a,i4,a,i4)')
     &  ' Max|era',trim(nmvar),'|','=',abs(era(ind(1)-2,ind(2)-2))/scl,
     &  '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
       if(all(var.ne.0.0_R8)) then
        err = abs(era)/exact_var
        ind = maxloc(abs(err))                                          !srv 09.07.18
        write(71,'(a,a,a,a,a,es10.3,a,i4,a,i4)')
     &  ' Max|era',trim(nmvar),'/',trim(nmvar),'|       =',
     &   abs(err(ind(1)-2,ind(2)-2)),
     &  '   is at ix = ', ind(1)-2,' iy = ', ind(2)-2
        call my_out(70,nx,ny,err,trim(nmvar)//'_rerr')
       else
        write(71,'(a,a,a,a,a,a4)')
     &   ' Max|era',trim(nmvar),'/',trim(nmvar),'|     =',' Inf'
       endif

       L2norm=norm2(exact_var)                                          !srv 09.07.18 {
       if(L2norm.ne.0.0_R8) then
        write(71,'(a,a,a,a,a,es10.3)')
     &   ' L2era',trim(nmvar),'/L2',trim(nmvar),'        =',
     &   norm2(era)/L2norm
        write(71,'(a,a,a,a,a,es10.3)')
     &   ' C_era',trim(nmvar),'/C_',trim(nmvar),'        =',
     &   maxval(abs(era))/maxval(abs(exact_var))
       else
        write(71,'(a,a,a,a,a,a4)')
     &   ' L2era',trim(nmvar),'/L2',trim(nmvar),'        =',' Inf'
        write(71,'(a,a,a,a,a,a4)')
     &   ' C_era',trim(nmvar),'/C_',trim(nmvar),'        =',' Inf'
       endif                                                             !srv 09.07.18 }

       write(71,*)

      end subroutine calc_err_cc                                        !srv 26.02.18 }
