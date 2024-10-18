c
      subroutine fy(ix)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_user_namelist
      use b2mod_subsys
      use b2us_data
      implicit none
      integer, intent(in) :: ix
      character*12 hlp_frm
      character*44 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring
      real (kind=R8) ::
     *   lmin,lmax,xmin,xmax,psign,poffset,pnorm,x1,x2
#ifdef WG_TODO
      integer ixc
      real (kind=R8) :: y_transform
#endif
      real (kind=R8), allocatable :: ds(:), wd(:), wrk(:)
      integer, allocatable :: cvlist(:)
      integer i,j,iyl,iyr,iy,is,isep,iCv,iCn,iFc,iVx,nlist
      integer iystart,iyend
      integer len_of_digits
      external xertst, len_of_digits

      call subini ('fy')
      if(mpg%isClassicalGrid.eq.0) return
      iystart=-1
      iCv=mpg%imapCv(ix,iystart)
      do while (iCv.eq.0 .and. iystart.le.mpg%ny)
        iystart=iystart+1
        iCv=mpg%imapCv(ix,iystart)
      enddo
      call xertst(iystart.le.mpg%ny, 'faulty parameter iystart')
      iyend=mpg%ny
      iCv=mpg%imapCv(ix,iyend)
      do while (iCv.eq.0 .and. iyend.gt.-1)
        iyend=iyend-1
        iCv=mpg%imapCv(ix,iyend)
      enddo
      call xertst(iystart.le.iyend, 'faulty parameter iyend')
      if (b2plot_debug.ge.6)
     > write(6,*) ' iystart, iyend : ', iystart, iyend
      write(hlp_frm,'(a,i1,a)') '(a,i',len_of_digits(ix),',a)'
      iCv = mpg%imapCv(ix,iystart)
      if (.not.map_to_midplane.or.nomp.eq.0) then
        nlist = iyend-iystart+1
        allocate(ds(nlist),wd(nlist))
        allocate(cvlist(nlist))
        isep = 0
        cvlist(1) = iCv
        do iy = iystart+1, iyend
          cvlist(iy-iystart+1) = mpg%imapCv(ix,iy)
          do i = 1, mpg%cvFcP(cvlist(iy-iystart),2)
            iFc = mpg%cvFc(mpg%cvFcP(cvlist(iy-iystart),1)+i-1)
            if (mpg%fcFs(iFc).eq.0) cycle
            if ((mpg%fcCv(iFc,1).eq.cvlist(iy-iystart).and.
     &           mpg%fcCv(iFc,2).eq.cvlist(iy-iystart+1)).or.
     &          (mpg%fcCv(iFc,2).eq.cvlist(iy-iystart).and.
     &           mpg%fcCv(iFc,1).eq.cvlist(iy-iystart+1))) then
              if (mpg%fcFs(iFc).eq.mpg%iFssep.or.
     &            mpg%fcFs(iFc).eq.mpg%iFssep2) isep = iy-iystart
              if (isep.eq.0 .and. mpg%nXpt.gt.1) then
                do j = 1, mpg%nXpt
                  if (mpg%vxFs(mpg%Xpt(j)).eq.mpg%fcFs(iFc)) then
                    if (isep.eq.0) isep = iy-iystart
                  end if
                end do
              end if
            end if
          end do
        end do
        call set_distances(mpg,geo,nlist,cvlist,isep,ds,wd)
      else
        if (mpg%nnreg(0).eq.8 .and. mpg%cvReg(iCv).le.4 .and.
     >      nimp.gt.0) then
          if (b2plot_debug.ge.5)
     >     write(*,*) 'Mapping to inner midplane'
          allocate(ds(nimp),cvlist(nimp))
          isep = 0
          nlist = 0
          iyl = iyend
          iyr = iystart
          do iy = nimp, 1, -1
            if (cvIy(imp(iy)).ge.iystart.and.
     &          cvIy(imp(iy)).le.iyend) then
              iyl = min(iyl,cvIy(imp(iy)))
              iyr = max(iyr,cvIy(imp(iy)))
              nlist = nlist+1
              ds(nlist) = dsimp(iy)
              cvlist(nlist) = imp(iy)
              if (iy.eq.icsepimp-1) isep = nlist
            end if
          end do
          iystart = iyl
          iyend = iyr
        else
          if (b2plot_debug.ge.5)
     >     write(*,*) 'Mapping to outer midplane'
          allocate(ds(nomp),cvlist(nomp))
          nlist = 0
          iyl = iyend
          iyr = iystart
          do iy = 1, nomp
            if (cvIy(omp(iy)).ge.iystart.and.
     &          cvIy(omp(iy)).le.iyend) then
              iyl = min(iyl,cvIy(omp(iy)))
              iyr = max(iyr,cvIy(omp(iy)))
              nlist = nlist+1
              ds(nlist) = dsomp(iy)
              cvlist(nlist) = omp(iy)
              if (iy.eq.icsepomp-1) isep = nlist
            end if
          end do
          iystart = iyl
          iyend = iyr
        end if
#ifdef WG_TODO
!cxpb Special mapping for profiles within PFRs
!cxpb Need to decide on a convention
        if (mpg%cvReg(iCv).eq.3) then
          ixc = leftcut(1)-1
        else if (mpg%nnreg(0).eq.8 .and. mpg%cvReg(iCv).eq.4) then
          ixc = rightcut(2)
        else if (mpg%nnreg(0).eq.7 .and. mpg%cvReg(iCv).eq.4) then
          ixc = rightcut(1)
        else if (mpg%nnreg(0).eq.4 .and. mpg%cvReg(iCv).eq.4) then
          ixc = rightcut(1)
        else if (mpg%nnreg(0).eq.7 .and. mpg%cvReg(iCv).eq.5) then
          ixc = leftcut(2)
        else if (mpg%nnreg(0).eq.7 .and. mpg%cvReg(iCv).eq.6) then
          ixc = rightcut(2)-1
        else if (mpg%nnreg(0).eq.7 .and. mpg%cvReg(iCv).eq.7) then
          ixc = rightcut(1)
        else if (mpg%nnreg(0).eq.8 .and. mpg%cvReg(iCv).eq.7) then
          ixc = leftcut(2)-1
        else if (mpg%nnreg(0).eq.8 .and. mpg%cvReg(iCv).eq.8) then
          ixc = rightcut(1)
        end if
#endif
      endif
      allocate(wrk(iystart:iyend))
      do is=0,ndim(2,iwrk)-1
        tmpstring1=textwk(is,iwrk)
        iCv=mpg%imapCv(ix,iystart)
        lmin=work(iCv,is,iwrk)
        lmax=work(iCv,is,iwrk)
        if(y_transform_l .ne. 0) then
#ifdef WG_TODO
          do iy=iystart,jsep
            iCv=mpg%imapCv(ix,iy)
            wrk(iy)=y_transform(cr(ixc,iy),cz(ixc,iy))
            write(*,*) wrk(iy),cz(ixc,iy)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          enddo
          do iy=jsep+1,iyend
            iCv=mpg%imapCv(ix,iy)
            wrk(iy)=y_transform(cr(ixx,iy),cz(ixx,iy))
            write(*,*) wrk(iy),cz(ixx,iy)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          enddo
#endif
        else if (plot_vs_norm_psi.or.plot_vs_full_psi.or.
     &           plot_vs_root_psi) then
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(cvlist(1),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(cvlist(1),1)))
          do i = 2, mpg%cvVxP(cvlist(1),2)
            iVx = mpg%cvVx(mpg%cvVxP(cvlist(1),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x1 = (xmin+xmax)/2.0_R8
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(cvlist(nlist),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(cvlist(nlist),1)))
          do i = 2, mpg%cvVxP(cvlist(nlist),2)
            iVx = mpg%cvVx(mpg%cvVxP(cvlist(nlist),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x2 = (xmin+xmax)/2.0_R8
          if (x1.eq.x2) then
            write(6,*) 'No psi variation over the range to plot'
            write(6,*) 'Graph abandoned'
            if(allocated(ds)) deallocate(ds)
            if(allocated(wd)) deallocate(wd)
            if(allocated(wrk)) deallocate(wrk)
            if(allocated(cvlist)) deallocate(cvlist)
            call subend ()
            return
          else
            psign = sign(1.0_R8,x2-x1)
            if (plot_vs_norm_psi.or.plot_vs_root_psi) then
              poffset = 1.0_R8
              pnorm = geo%fsPsi(mpg%iFssep)
            end if
            do iy = iystart, iyend
              iCv=mpg%imapCv(ix,iy)
              xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
              xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
              do i = 2, mpg%cvVxP(iCv,2)
                iVx = mpg%cvVx(mpg%cvVxP(iCv,1)+i-1)
                xmin = min(xmin,geo%vxFpsi(iVx))
                xmax = max(xmax,geo%vxFpsi(iVx))
              end do
              if (plot_vs_norm_psi.or.plot_vs_root_psi) then
                wrk(iy)=(((xmin+xmax)/2.0_R8-pnorm)*psign/pnorm)+poffset
              else
                wrk(iy)=(xmin+xmax)/2.0_R8
              end if
              if (plot_vs_root_psi) wrk(iy)=sqrt(wrk(iy))
              lmin=min(lmin,work(iCv,is,iwrk))
              lmax=max(lmax,work(iCv,is,iwrk))
            end do
          endif
        else if (plot_vs_rmaj) then
          x1 = geo%cvX(cvlist(1))
          x2 = geo%cvX(cvlist(nlist))
          if (x1.eq.x2) then
            write(6,*)
     .       'No major radius variation over the range to plot'
            write(6,*) 'Graph abandoned'
            if(allocated(ds)) deallocate(ds)
            if(allocated(wd)) deallocate(wd)
            if(allocated(wrk)) deallocate(wrk)
            if(allocated(cvlist)) deallocate(cvlist)
            call subend ()
            return
          else
            do iy=iystart,iyend
              iCv=mpg%imapCv(ix,iy)
              wrk(iy)=geo%cvX(iCv)
              lmin=min(lmin,work(iCv,is,iwrk))
              lmax=max(lmax,work(iCv,is,iwrk))
            end do
          endif
        else
          do iy=iystart,iyend
            wrk(iy)=ds(iy-iystart+1)
            iCv=mpg%imapCv(ix,iy)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          enddo
        endif
        if(values(9) .ne. values(10)) then
          lmin=values(9)
          lmax=values(10)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9_R8,lmin*1.1_R8)
          lmax=max(lmax*0.9_R8,lmax*1.1_R8)
          if(lmin.eq.lmax) then
            lmin=-1.0_R8
            lmax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(mod(typeofplot,2) .eq. 1) then
          if (plot_vs_norm_psi) then
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'> (vs \Gw\R)'
          else if (plot_vs_full_psi) then
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'> (vs \GW\R)'
          else if (plot_vs_root_psi) then
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'> (vs \Gq\R_{pol})'
          else if (plot_vs_rmaj) then
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'> (vs R)'
          else if (map_to_midplane) then
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'> '//
     &       '(mapped to midplane)'
          else
            write(tmpstring2,hlp_frm) ' <ix = ',ix,'>'
          endif
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmin=values(11)
          else
            xmin=min(wrk(iystart),wrk(iyend))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(iystart),wrk(iyend))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.true.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do iy=iystart,iyend
              iCv=mpg%imapCv(ix,iy)
              write(nwtf,*) wrk(iy),work(iCv,is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do iy=iystart,iyend-1
            iCv=mpg%imapCv(ix,iy)
            iCn=mpg%imapCv(ix,iy+1)
            call lined_d(
     1       wrk(iy),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(iy+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
            if(cvIy(iCn).ne.cvIy(iCv)+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(iy)+wrk(iy+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(iy)+wrk(iy+1))/2.0_R8,lmin,
     1                       (wrk(iy)+wrk(iy+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,lfsz)
#endif
            endif
          enddo
          if (plotting_sep.and.(isep.ne.0.or.
     &       (iyref.ge.iystart.and.iyref.le.iyend))) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if (isep.ne.0) then
              if ((wrk(isep+iystart-1)+wrk(isep+iystart)).gt.0.0_R8
     &            .or.linlog.le.2)
     >          call lined_d(
     1           (wrk(isep+iystart-1)+wrk(isep+iystart))/2.0_R8,lmin,
     1           (wrk(isep+iystart-1)+wrk(isep+iystart))/2.0_R8,lmax)
            else
              if ((wrk(iyref)+wrk(iyref+1)).gt.0.0_R8.or.linlog.le.2)
     >         call lined_d(
     1           (wrk(iyref)+wrk(iyref+1))/2.0_R8,lmin,
     1           (wrk(iyref)+wrk(iyref+1))/2.0_R8,lmax)
            end if
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        endif
        if(typeofplot/2 .eq. 1 .and. mpg%isClassicalGrid .eq. 1) then
          do iy=iystart,iyend
            iCv=mpg%imapCv(ix,iy)
            wrk(iy)=yy(iCv,2)
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          enddo
          if(values(9) .ne. values(10)) then
            lmin=values(9)
            lmax=values(10)
          endif
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9_R8,lmin*1.1_R8)
            lmax=max(lmax*0.9_R8,lmax*1.1_R8)
            if(lmin.eq.lmax) then
              lmin=-1.0_R8
              lmax=+1.0_R8
            end if
          end if
          if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(3).ne.values(4)) then
            xmin=values(3)
            xmax=values(4)
          else
            xmin=min(wrk(iystart),wrk(iyend))
            xmax=max(wrk(iystart),wrk(iyend))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          write(tmpstring2,hlp_frm) ' <ix = ',ix,'>'
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do iy=iystart,iyend
              iCv=mpg%imapCv(ix,iy)
              write(nwtf,*) wrk(iy),work(iCv,is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            if(typeofplot.eq.3.and.mpg%isClassicalGrid.eq.1)
     >       dcol=dcol+1
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do iy=iystart,iyend-1
            if(wrk(iy+1).lt.xmin) cycle
            if(wrk(iy).gt.xmax) cycle
            iCv=mpg%imapCv(ix,iy)
            iCn=mpg%imapCv(ix,iy+1)
            call lined_d(
     1       wrk(iy),
     2       min(lmax,max(lmin,work(iCv,is,iwrk))),
     3       wrk(iy+1),
     4       min(lmax,max(lmin,work(iCn,is,iwrk))))
            if(cvIy(iCn).ne.cvIy(iCv)+1) then
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
              call gslwsc_d(sepwidth)
              call dashdc('$',1,0)
#endif
              if ((wrk(iy)+wrk(iy+1)).gt.0.0_R8.or.linlog.le.2)
     >          call lined_d((wrk(iy)+wrk(iy+1))/2.0_R8,lmin,
     1                       (wrk(iy)+wrk(iy+1))/2.0_R8,lmax)
#ifndef PLPLOT
              call lastd
              if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
              call gslwsc_d(lwidth)
              call dashdc(linefmt,8,lfsz)
#endif
            endif
          enddo
          if (plotting_sep.and.(isep.ne.0.or.
     &       (iyref.ge.iystart.and.iyref.le.iyend))) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if (isep.ne.0) then
              if ((wrk(isep+iystart-1)+wrk(isep+iystart)).gt.0.0_R8
     &            .or.linlog.le.2)
     >          call lined_d(
     1           (wrk(isep+iystart-1)+wrk(isep+iystart))/2.0_R8,lmin,
     1           (wrk(isep+iystart-1)+wrk(isep+iystart))/2.0_R8,lmax)
            else
              if ((wrk(iyref)+wrk(iyref+1)).gt.0.0_R8.or.linlog.le.2)
     >         call lined_d(
     1           (wrk(iyref)+wrk(iyref+1))/2.0_R8,lmin,
     1           (wrk(iyref)+wrk(iyref+1))/2.0_R8,lmax)
            end if
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        endif
      enddo
      if(allocated(ds)) deallocate(ds)
      if(allocated(wd)) deallocate(wd)
      if(allocated(wrk)) deallocate(wrk)
      if(allocated(cvlist)) deallocate(cvlist)
      call subend ()
      return
      end subroutine fy
c
      subroutine fimp()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_user_namelist
      use b2us_data
      implicit none
      integer i,j,is,iCv,iVx
      real (kind=R8) :: x1,x2,xmin,xmax,lmin,lmax,psign,pnorm,poffset
      real (kind=R8), allocatable :: wrk(:)
      character*44 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring

      if (nimp.eq.0) return
      if (b2plot_debug.ge.5)
     > write(*,*) 'Plotting along the inner midplane'
      allocate(wrk(1:nimp))
      do is=0,ndim(2,iwrk)-1
        tmpstring1=textwk(is,iwrk)
        lmin=work(imp(1),is,iwrk)
        lmax=work(imp(1),is,iwrk)
        if (plot_vs_norm_psi.or.plot_vs_full_psi.or.
     &      plot_vs_root_psi) then
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(imp(1),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(imp(1),1)))
          do i = 2, mpg%cvVxP(imp(1),2)
            iVx = mpg%cvVx(mpg%cvVxP(imp(1),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x1 = (xmin+xmax)/2.0_R8
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(imp(nimp),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(imp(nimp),1)))
          do i = 2, mpg%cvVxP(imp(nimp),2)
            iVx = mpg%cvVx(mpg%cvVxP(imp(nimp),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x2 = (xmin+xmax)/2.0_R8
          psign = sign(1.0_R8,x2-x1)
          if (plot_vs_norm_psi.or.plot_vs_root_psi) then
            poffset = 1.0_R8
            pnorm = geo%fsPsi(mpg%iFssep)
          end if
          do j = 1, nimp
            iCv=imp(j)
            xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
            xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
            do i = 2, mpg%cvVxP(iCv,2)
              iVx = mpg%cvVx(mpg%cvVxP(iCv,1)+i-1)
              xmin = min(xmin,geo%vxFpsi(iVx))
              xmax = max(xmax,geo%vxFpsi(iVx))
            end do
            if (plot_vs_norm_psi.or.plot_vs_root_psi) then
              wrk(j)=(((xmin+xmax)/2.0_R8-pnorm)*psign/pnorm)+poffset
            else
              wrk(j)=(xmin+xmax)/2.0_R8
            end if
            if (plot_vs_root_psi) wrk(j)=sqrt(wrk(j))
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          end do
        else if (plot_vs_rmaj) then
          do i=1,nimp
            wrk(i)=geo%cvX(imp(i))
            lmin=min(lmin,work(imp(i),is,iwrk))
            lmax=max(lmax,work(imp(i),is,iwrk))
          end do
        else
          do i=1,nimp
            wrk(i)=dsimp(i)
            lmin=min(lmin,work(imp(i),is,iwrk))
            lmax=max(lmax,work(imp(i),is,iwrk))
          enddo
        end if
        if(values(9) .ne. values(10)) then
          lmin=values(9)
          lmax=values(10)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9_R8,lmin*1.1_R8)
          lmax=max(lmax*0.9_R8,lmax*1.1_R8)
          if(lmin.eq.lmax) then
            lmin=-1.0_R8
            lmax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(mod(typeofplot,2) .eq. 1) then
          if (plot_vs_norm_psi) then
            write(tmpstring2,'(a)') ' (IMP vs \Gw\R)'
          else if (plot_vs_full_psi) then
            write(tmpstring2,'(a)') ' (IMP vs \GW\R)'
          else if (plot_vs_root_psi) then
            write(tmpstring2,'(a)') ' (IMP vs \Gq\R_{pol})'
          else if (plot_vs_rmaj) then
            write(tmpstring2,'(a)') ' (IMP vs R)'
          else
            write(tmpstring2,'(a)') ' (IMP)'
          endif
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmin=values(11)
          else
            xmin=min(wrk(1),wrk(nimp))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(1),wrk(nimp))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.true.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do i=1,nimp
              write(nwtf,*) wrk(i),work(imp(i),is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do i=1,nimp-1
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(imp(i),is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(imp(i+1),is,iwrk))))
          end do
          if (plotting_sep.and.icsepimp.gt.1) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if ((wrk(icsepimp-1)+wrk(icsepimp)).gt.0.0_R8
     &          .or.linlog.le.2)
     >        call lined_d(
     1         (wrk(icsepimp-1)+wrk(icsepimp))/2.0_R8,lmin,
     1         (wrk(icsepimp-1)+wrk(icsepimp))/2.0_R8,lmax)
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          end if
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
        if(typeofplot/2 .eq. 1) then
          do i=1,nimp
            wrk(i)=real(i,R8)
            lmin=min(lmin,work(imp(i),is,iwrk))
            lmax=max(lmax,work(imp(i),is,iwrk))
          end do
          if(values(9) .ne. values(10)) then
            lmin=values(9)
            lmax=values(10)
          end if
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9_R8,lmin*1.1_R8)
            lmax=max(lmax*0.9_R8,lmax*1.1_R8)
            if(lmin.eq.lmax) then
              lmin=-1.0_R8
              lmax=+1.0_R8
            end if
          end if
          if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(3).ne.values(4)) then
            xmin=values(3)
            xmax=values(4)
          else
            xmin=min(wrk(1),wrk(nimp))
            xmax=max(wrk(1),wrk(nimp))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          write(tmpstring2,'(a)') ' (IMP)'
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do i=1,nimp
              write(nwtf,*) wrk(i),work(imp(i),is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            if(typeofplot.eq.3) dcol=dcol+1
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do i=1,nimp-1
            if(wrk(i+1).lt.xmin) cycle
            if(wrk(i).gt.xmax) cycle
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(imp(i),is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(imp(i+1),is,iwrk))))
          end do
          if (plotting_sep.and.icsepimp.gt.1) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if ((wrk(icsepimp-1)+wrk(icsepimp)).gt.0.0_R8
     &          .or.linlog.le.2)
     >        call lined_d(
     1         (wrk(icsepimp-1)+wrk(icsepimp))/2.0_R8,lmin,
     1         (wrk(icsepimp-1)+wrk(icsepimp))/2.0_R8,lmax)
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          end if
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
      end do
      deallocate(wrk)

      return
      end subroutine fimp
c
      subroutine fomp()
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_user_namelist
      use b2us_data
      implicit none
      integer i,j,is,iCv,iVx
      real (kind=R8) :: x1,x2,xmin,xmax,lmin,lmax,psign,pnorm,poffset
      real (kind=R8), allocatable :: wrk(:)
      character*44 tmpstring1,tmpstring2
      character*256 tmpstring,titlestring

      if (nomp.eq.0) return
      if (b2plot_debug.ge.5)
     > write(*,*) 'Plotting along the outer midplane'
      allocate(wrk(1:nomp))
      do is=0,ndim(2,iwrk)-1
        tmpstring1=textwk(is,iwrk)
        lmin=work(omp(1),is,iwrk)
        lmax=work(omp(1),is,iwrk)
        if (plot_vs_norm_psi.or.plot_vs_full_psi.or.
     &      plot_vs_root_psi) then
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(omp(1),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(omp(1),1)))
          do i = 2, mpg%cvVxP(omp(1),2)
            iVx = mpg%cvVx(mpg%cvVxP(omp(1),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x1 = (xmin+xmax)/2.0_R8
          xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(omp(nomp),1)))
          xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(omp(nomp),1)))
          do i = 2, mpg%cvVxP(omp(nomp),2)
            iVx = mpg%cvVx(mpg%cvVxP(omp(nomp),1)+i-1)
            xmin = min(xmin,geo%vxFpsi(iVx))
            xmax = max(xmax,geo%vxFpsi(iVx))
          end do
          x2 = (xmin+xmax)/2.0_R8
          psign = sign(1.0_R8,x2-x1)
          if (plot_vs_norm_psi.or.plot_vs_root_psi) then
            poffset = 1.0_R8
            pnorm = geo%fsPsi(mpg%iFssep)
          end if
          do j = 1, nomp
            iCv=omp(j)
            xmin = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
            xmax = geo%vxFpsi(mpg%cvVx(mpg%cvVxP(iCv,1)))
            do i = 2, mpg%cvVxP(iCv,2)
              iVx = mpg%cvVx(mpg%cvVxP(iCv,1)+i-1)
              xmin = min(xmin,geo%vxFpsi(iVx))
              xmax = max(xmax,geo%vxFpsi(iVx))
            end do
            if (plot_vs_norm_psi.or.plot_vs_root_psi) then
              wrk(j)=(((xmin+xmax)/2.0_R8-pnorm)*psign/pnorm)+poffset
            else
              wrk(j)=(xmin+xmax)/2.0_R8
            end if
            if (plot_vs_root_psi) wrk(j)=sqrt(wrk(j))
            lmin=min(lmin,work(iCv,is,iwrk))
            lmax=max(lmax,work(iCv,is,iwrk))
          end do
        else if (plot_vs_rmaj) then
          do i=1,nomp
            wrk(i)=geo%cvX(omp(i))
            lmin=min(lmin,work(omp(i),is,iwrk))
            lmax=max(lmax,work(omp(i),is,iwrk))
          end do
        else
          do i=1,nomp
            wrk(i)=dsomp(i)
            lmin=min(lmin,work(omp(i),is,iwrk))
            lmax=max(lmax,work(omp(i),is,iwrk))
          enddo
        end if
        if(values(9) .ne. values(10)) then
          lmin=values(9)
          lmax=values(10)
        endif
        if(lmin.eq.lmax) then
          lmin=min(lmin*0.9_R8,lmin*1.1_R8)
          lmax=max(lmax*0.9_R8,lmax*1.1_R8)
          if(lmin.eq.lmax) then
            lmin=-1.0_R8
            lmax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
        if(mod(typeofplot,2) .eq. 1) then
          if (plot_vs_norm_psi) then
            write(tmpstring2,'(a)') ' (OMP vs \Gw\R)'
          else if (plot_vs_full_psi) then
            write(tmpstring2,'(a)') ' (OMP vs \GW\R)'
          else if (plot_vs_root_psi) then
            write(tmpstring2,'(a)') ' (OMP vs \Gq\R_{pol})'
          else if (plot_vs_rmaj) then
            write(tmpstring2,'(a)') ' (OMP vs R)'
          else
            write(tmpstring2,'(a)') ' (OMP)'
          endif
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmin=values(11)
          else
            xmin=min(wrk(1),wrk(nomp))
          endif
          if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8)) then
            xmax=values(12)
          else
            xmax=max(wrk(1),wrk(nomp))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.true.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do i=1,nomp
              write(nwtf,*) wrk(i),work(omp(i),is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do i=1,nomp-1
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(omp(i),is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(omp(i+1),is,iwrk))))
          end do
          if (plotting_sep.and.icsepomp.gt.1) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if ((wrk(icsepomp-1)+wrk(icsepomp)).gt.0.0_R8
     &          .or.linlog.le.2)
     >        call lined_d(
     1         (wrk(icsepomp-1)+wrk(icsepomp))/2.0_R8,lmin,
     1         (wrk(icsepomp-1)+wrk(icsepomp))/2.0_R8,lmax)
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          end if
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
        if(typeofplot/2 .eq. 1) then
          do i=1,nomp
            wrk(i)=real(i,R8)
            lmin=min(lmin,work(omp(i),is,iwrk))
            lmax=max(lmax,work(omp(i),is,iwrk))
          end do
          if(values(9) .ne. values(10)) then
            lmin=values(9)
            lmax=values(10)
          end if
          if(lmin.eq.lmax) then
            lmin=min(lmin*0.9_R8,lmin*1.1_R8)
            lmax=max(lmax*0.9_R8,lmax*1.1_R8)
            if(lmin.eq.lmax) then
              lmin=-1.0_R8
              lmax=+1.0_R8
            end if
          end if
          if(mod(linlog,2).eq.0.and.lmin.le.0.0_R8) lmin=1.0e-30_R8
          if(mod(linlog,2).eq.0.and.lmax.le.0.0_R8) lmax=1.0_R8
          if(values(3).ne.values(4)) then
            xmin=values(3)
            xmax=values(4)
          else
            xmin=min(wrk(1),wrk(nomp))
            xmax=max(wrk(1),wrk(nomp))
          endif
          if(xmin.eq.xmax) then
            xmin=min(xmin*0.9_R8,xmin*1.1_R8)
            xmax=max(xmax*0.9_R8,xmax*1.1_R8)
            if(xmin.eq.xmax) then
              xmin=-1.0_R8
              xmax=+1.0_R8
            end if
          end if
          if(linlog.gt.2.and.xmin.le.0.0_R8) xmin=1.0e-30_R8
          if(linlog.gt.2.and.xmax.le.0.0_R8) xmax=1.0_R8
          write(tmpstring2,'(a)') ' (OMP)'
          tmpstring =
     1     trim(graphlabel(iwrk))//trim(tmpstring2)//trim(tmpstring1)
          call customise_label(tmpstring,titlestring)
          call init_page(titlestring,
     1     xmin,xmax,lmin,lmax,
     4     linlog,.false.,.false.,.false.)
          if(writing_to_file) then
            write(nwtf,'(''# '',a)') trim(titlestring)
            do i=1,nomp
              write(nwtf,*) wrk(i),work(omp(i),is,iwrk)
            enddo
          endif
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(lwidth)
          if (dcol.gt.1) then
            call gsplci(dcol)
          elseif(dcol.lt.0) then
            if(typeofplot.eq.3) dcol=dcol+1
            call gsplci(-dcol)
            dcol=dcol-1
            if (dcol.lt.-17) then
              dcol=-mod(abs(dcol),17)
              if (dcol.eq.0) dcol=-17
            endif
          endif
#endif
          do i=1,nomp-1
            if(wrk(i+1).lt.xmin) cycle
            if(wrk(i).gt.xmax) cycle
            call lined_d(
     1       wrk(i),
     2       min(lmax,max(lmin,work(omp(i),is,iwrk))),
     3       wrk(i+1),
     4       min(lmax,max(lmin,work(omp(i+1),is,iwrk))))
          end do
          if (plotting_sep.and.icsepomp.gt.1) then
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(scol)
            call gslwsc_d(sepwidth)
            call dashdc('$',1,0)
#endif
            if ((wrk(icsepomp-1)+wrk(icsepomp)).gt.0.0_R8
     &          .or.linlog.le.2)
     >        call lined_d(
     1         (wrk(icsepomp-1)+wrk(icsepomp))/2.0_R8,lmin,
     1         (wrk(icsepomp-1)+wrk(icsepomp))/2.0_R8,lmax)
#ifndef PLPLOT
            call lastd
            if (scol.gt.1 .or. dcol.ne.1) call gsplci(abs(dcol))
            call gslwsc_d(lwidth)
            call dashdc(linefmt,8,lfsz)
#endif
          end if
#ifndef PLPLOT
          call lastd
          if (lwidth.ne.1.0_R8) call gslwsc_d(1.0_R8)
          if (dcol.ne.1) call gsplci(1)
#endif
          call end_page
        end if
      end do
      deallocate(wrk)

      return
      end subroutine fomp
c
      function y_transform(r,z)
      use b2mod_types
      use b2mod_subsys
      real (kind=R8) ::
     *   y_transform
      real (kind=R8) ::
     *   r,z
      real (kind=R8) ::
     *   z1(100),s1(100),z2(100),s2(100)
      logical first,open_file
      integer n1,n2,i1,i2
      real (kind=R8) ::
     *   x, r0r
      external xertst
      save z1,s1,z2,s2,n1,n2,first
      data first /.true./

      call subini('y_transform')
      write(*,*) 'y_transform'
      if(first) then
        call xertst(open_file(99,'n.s.r.z.a_inner',
     1   'data.local|data'), 'Could not open n.s.r.z.a_inner')
        n1=0
        read(99,*)
   10   read(99,*,end=20) x,s1(n1+1),x,z1(n1+1)
        n1=n1+1
        goto 10
   20   continue
        write(*,*) 'Read n.s.r.z.a_inner',n1
        close(99)
        call xertst(open_file(99,'n.s.r.z.a_outer',
     1   'data.local|data'), 'Could not open n.s.r.z.a_outer')
        n2=0
        read(99,*)
   30   read(99,*,end=40) x,s2(n2+1),x,z2(n2+1)
        n2=n2+1
        goto 30
   40   continue
        write(*,*) 'Read n.s.r.z.a_outer',n2
        close(99)
      endif
      first=.false.

c     if(r0r .eq. 0.0_R8) r0r=1.5_R8
      r0r=1.5_R8
      if(r .lt. r0r) then
        write(*,*) 'r <'
        i1=int((z-z1(1))/(z1(n1)-z1(1))*n1)+1
        i1=min(n1-1,max(1,i1))
        do while(z1(i1) .lt. z .and. i1 .gt. 1)
          i1=i1-1
          write(*,*) '-: i1',i1
        enddo
        do while(z1(i1+1) .gt. z .and. i1 .lt. n1-1)
          i1=i1+1
          write(*,*) '+: i1',i1
        enddo
        y_transform=(z-z1(i1))/(z1(i1+1)-z1(i1))*
     1   (s1(i1+1)-s1(i1))+s1(i1)
      else
        write(*,*) 'r >'
        i2=int((z-z2(1))/(z2(n2)-z2(1))*n2)+1
        i2=min(n2-1,max(1,i2))
        do while(z2(i2) .gt. z .and. i2 .gt. 1)
          i2=i2-1
          write(*,*) '-: i1',i1
        enddo
        do while(z2(i2+1) .lt. z .and. i2 .lt. n2-1)
          i2=i2+1
          write(*,*) '+: i1',i1
        enddo
        y_transform=(z-z2(i2))/(z2(i2+1)-z2(i2))*
     1   (s2(i2+1)-s2(i2))+s2(i2)
      endif

      call subend ()
      return
      end function y_transform

      subroutine plot_expfile
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      character*256 line, title
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax
      logical open_file, is_comment
      external is_comment
      integer lngarr,i,inc
      parameter (lngarr=100 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr)

      write(*,*) 'Called plot_expfile'

      if(.not.open_file(99,expfile,
     1 'data.local/exp|data/exp')) then
        write(*,*) 'Could not open ',trim(expfile)
      else
        title = ''
   10   read(99,'(a)',end=40) line
        if(is_comment(line)) then
          title = trim(line)
          write(*,*) trim(line)
          goto 10
        endif
        backspace(99)
        inc=0
   20   read(99,*,end=30) xvalues(inc+1),yvalues(inc+1)
        inc=inc+1
        goto 20
   30   continue
        close(99)
        xmin=xvalues(1)
        xmax=xvalues(1)
        ymin=yvalues(1)
        ymax=yvalues(1)
        do i=2,inc
          xmin=min(xmin,xvalues(i))
          xmax=max(xmax,xvalues(i))
          ymin=min(ymin,yvalues(i))
          ymax=max(ymax,yvalues(i))
        enddo
        if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8))
     >   xmin=values(11)
        if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8))
     >   xmax=values(12)
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(values(14) .ne. values(15)) then
          ymin=values(14)
          ymax=values(15)
        endif
        if(ymin.eq.ymax) then
          ymin=min(ymin*0.9_R8,ymin*1.1_R8)
          ymax=max(ymax*0.9_R8,ymax*1.1_R8)
          if(ymin.eq.ymax) then
            ymin=-1.0_R8
            ymax=+1.0_R8
          end if
        end if
        if(mod(linlog,2).eq.0.and.ymin.le.0.0_R8) ymin=1.0e-30_R8
        if(mod(linlog,2).eq.0.and.ymax.le.0.0_R8) ymax=1.0_R8
        call init_page(
     1   trim(title),
     1   xmin,xmax,ymin,ymax,linlog,.false.,.false.,.true.)
         if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') trim(title)
          do i=1,inc
            write(nwtf,*) xvalues(i),yvalues(i)
          enddo
        endif
#ifndef PLPLOT
        call lastd
        if(ecol .gt. 1) then
          if (expline) then
            call gsplci(ecol)
          else
            call gspmci(ecol)
          endif
        elseif(ecol .lt. 0) then
          if (expline) then
            call gsplci(-ecol)
          else
            call gspmci(-ecol)
          endif
          ecol=ecol-1
        endif
        if (ewidth.ne.1.0_R8) then
          call lastd
          if (expline) then
            call gslwsc_d(ewidth)
          else
            call gsmksc_d(ewidth)
          endif
        endif
#endif
        if (expline) then
          do i=1,inc-1
            call lined_d(xvalues(i),yvalues(i),
     .                   xvalues(i+1),yvalues(i+1))
          enddo
        else
          do i=1,inc
            call grpts_d(xvalues(i),yvalues(i),expsymbol)
          enddo
        endif
#ifndef PLPLOT
        if (ewidth.ne.1.0_R8) then
          call lastd
          if (expline) then
            call gslwsc_d(1.0_R8)
          else
            call gsmksc_d(1.0_R8)
          endif
        endif
#endif
        call end_page
   40   continue
        if(ecol .ne. 1) then
#ifndef PLPLOT
          call lastd
          if (expline) then
            call gsplci(1)
          else
            call gspmci(1)
          endif
#endif
        endif
      endif
      return
      end subroutine plot_expfile

      subroutine plot_wallfile
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      implicit none
      character*256 line, title
      real (kind=R8) ::
     *   xmin,xmax,ymin,ymax
      logical open_file, is_comment
      external is_comment
      integer lngarr,i,inc
      parameter (lngarr=100 000)
      real (kind=R8) ::
     *   xvalues(lngarr),yvalues(lngarr)

      write(*,*) 'Called plot_wallfile'

      if(.not.open_file(99,wallfile,
     1 'data.local/exp|data/exp')) then
        write(*,*) 'Could not open ',trim(wallfile)
      else
        title = ''
        if(ecol .gt. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(ecol)
#endif
        elseif(ecol.lt.0) then
#ifndef PLPLOT
          call lastd
          call gsplci(-ecol)
#endif
          ecol=ecol-1
        endif
   10   read(99,'(a)',end=40) line
        if(is_comment(line)) then
          title = trim(line)
          write(*,*) trim(line)
          goto 10
        endif
        backspace(99)
        inc=0
   20   read(99,*,end=30) xvalues(inc+1),yvalues(inc+1)
        inc=inc+1
        goto 20
   30   continue
        close(99)
        xmin=xvalues(1)
        xmax=xvalues(1)
        ymin=yvalues(1)
        ymax=yvalues(1)
        do i=2,inc
          xmin=min(xmin,xvalues(i))
          xmax=max(xmax,xvalues(i))
          ymin=min(ymin,yvalues(i))
          ymax=max(ymax,yvalues(i))
        enddo
        if(values(11).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8))
     >   xmin=values(11)
        if(values(12).ne.init_value .and.
     &      .not.(values(11).eq.0.0_R8 .and. values(12).eq.0.0_R8))
     >   xmax=values(12)
        if(xmin.eq.xmax) then
          xmin=min(xmin*0.9_R8,xmin*1.1_R8)
          xmax=max(xmax*0.9_R8,xmax*1.1_R8)
          if(xmin.eq.xmax) then
            xmin=-1.0_R8
            xmax=+1.0_R8
          end if
        end if
        if(values(14) .ne. values(15)) then
          ymin=values(14)
          ymax=values(15)
        endif
        if(ymin.eq.ymax) then
          ymin=min(ymin*0.9_R8,ymin*1.1_R8)
          ymax=max(ymax*0.9_R8,ymax*1.1_R8)
          if(ymin.eq.ymax) then
            ymin=-1.0_R8
            ymax=+1.0_R8
          end if
        end if
        if(mod(wlinlog,2).eq.0.and.ymin.le.0.0_R8) ymin=1.0e-30_R8
        if(mod(wlinlog,2).eq.0.and.ymax.le.0.0_R8) ymax=1.0_R8
        call init_page(
     1   trim(title),
     1   xmin,xmax,ymin,ymax,wlinlog,.false.,.false.,.true.)
        if(b2plot_debug.ge.4)
     >   write(*,*) 'xmin,xmax,ymin,ymax ',xmin,xmax,ymin,ymax
        if(writing_to_file) then
          write(nwtf,'(''# '',a)') trim(title)
          do i=1,inc
            write(nwtf,*) xvalues(i),yvalues(i)
          enddo
        endif
#ifndef PLPLOT
        if (wwidth.ne.1.0_R8) then
          call lastd
          call gslwsc_d(wwidth)
        endif
#endif
        do i=1,inc-1
          call lined_d(xvalues(i),yvalues(i),xvalues(i+1),yvalues(i+1))
        enddo
#ifndef PLPLOT
        if (wwidth.ne.1.0_R8) then
          call lastd
          call gslwsc_d(1.0_R8)
        endif
#endif
        call end_page
   40   continue
        if(ecol .ne. 1) then
#ifndef PLPLOT
          call lastd
          call gsplci(1)
#endif
        endif
      endif
      return
      end subroutine plot_wallfile
