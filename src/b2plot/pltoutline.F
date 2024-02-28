c
      subroutine pltoutline(physical)
      use b2mod_types
      use b2mod_b2plot
      use b2us_data
      implicit none
      integer i, iCv, iFc, iVx1, iVx2
      logical physical

c       call gslwsc_d(1.0_R8)
#ifndef PLPLOT
      call lastd
      if(ocol .gt. 1) call gsplci(ocol)
      call gslwsc_d(owidth)
#endif
      if(physical) then ! mesh physical domain
        do iCv = mpg%nCi+1, mpg%nCv
          do i = 1, mpg%cvFcP(iCv,2)
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
            iVx1 = mpg%fcVx(iFc,1)
            iVx2 = mpg%fcVx(iFc,2)
            call lined_d(geo%vxX(iVx1),geo%vxY(iVx1),
     &                   geo%vxX(iVx2),geo%vxY(iVx2))
          end do
        end do
      else if (mpg%isClassicalGrid.eq.1) then
        do iCv = mpg%nCi+1, mpg%nCv
          do i = 1, mpg%cvFcP(iCv,2)
            iFc = mpg%cvFc(mpg%cvFcP(iCv,1)+i-1)
            if (fcxIx(iFc).gt.-2  .and. fcxIy(iFc).gt.-2) then
              call lined_d(real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)-0.5_R8,
     &                     real(fcxIx(iFc),R8)-0.5_R8,
     &                     real(fcxIy(iFc),R8)+0.5_R8)
            else if (fcyIx(iFc).gt.-2  .and. fcyIy(iFc).gt.-2) then
              call lined_d(real(fcyIx(iFc),R8)-0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8,
     &                     real(fcyIx(iFc),R8)+0.5_R8,
     &                     real(fcyIy(iFc),R8)-0.5_R8)
            endif
          end do
        end do
      endif
#ifndef PLPLOT
      call lastd
      if(ocol .gt. 1) call gsplci(1)
      call gslwsc_d(lwidth)
#endif
      return
      end subroutine pltoutline
c
      subroutine plot_arrows(is)
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_debug
      implicit none
      integer is,iCv
      real (kind=R8) ::
     * fmax,fmay,fmmax,rmagst,rmaguv,r,s,t,u,v,sb,cb,rl,rd,
     1 uu1,vv1,uu2,vv2,brx,bzy,dfx,dfy,uuu,vvv
      intrinsic max, sqrt

      fmax=arrowx(1,is)
      fmay=arrowy(1,is)
      fmmax=sqrt(fmax**2+fmay**2)
      do iCv=2,mpg%nCi
        fmax=arrowx(iCv,is)
        fmay=arrowy(iCv,is)
        fmmax=max(fmmax,sqrt(fmax*fmax+fmay*fmay))
      enddo
      if(fmmax.eq.0.0_R8) then
        write(*,*) 'Cannot plot arrows for a zero field'
        return
      endif
      if (b2plot_debug.ge.2) write (*,*) 'drawing arrows'
#ifndef PLPLOT
      call lastd
      call gslwsc_d(arroww) ! increase arrow line width
      if (acol .gt. 1) call gsplci(acol) ! change arrow colour
#endif
      do iCv=1,mpg%nCi
        brx=geo%cvX(iCv)
        bzy=geo%cvY(iCv)
        s= geo%cvEb(iCv,0)
        t= geo%cvEb(iCv,1)
        u= geo%cvEb(iCv,1)*signmf
        v=-geo%cvEb(iCv,0)*signmf
        rmagst=sqrt(s**2+t**2)
        rmaguv=sqrt(u**2+v**2)

        dfx=arrowx(iCv,is)/fmmax
        dfy=arrowy(iCv,is)/fmmax

        uuu=arrowl*(dfx*s/rmagst+dfy*u/rmaguv)
        vvv=arrowl*(dfx*t/rmagst+dfy*v/rmaguv)
        r=sqrt(uuu**2+vvv**2)
        sb=vvv/r
        cb=uuu/r
        rl=0.2_R8*r
        rd=0.1_R8*r
        uu1=uuu-rl*cb+rd*sb
        vv1=vvv-rl*sb-rd*cb
        uu2=uuu-rl*cb-rd*sb
        vv2=vvv-rl*sb+rd*cb
        call lined_d(brx,bzy,brx+uuu,bzy+vvv)
        call lined_d(brx+uuu,bzy+vvv,brx+uu1,bzy+vv1)
        call lined_d(brx+uuu,bzy+vvv,brx+uu2,bzy+vv2)
c       call lined_d(brx+uu1,bzy+vv1,brx+uu2,bzy+vv2)
      enddo
#ifndef PLPLOT
      call gslwsc_d(1.0_R8)
      if (acol .gt. 1) then
        call lastd
        call gsplci(1)
      endif
#endif
      if (b2plot_debug.ge.2) write (*,*) 'finished drawing arrows'
      return
      end subroutine plot_arrows
c
      subroutine plot_trig_arrows(is)
      use b2mod_types
      use b2mod_b2plot
      use b2mod_b2plot_debug
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif
      implicit none
      integer is
#ifdef B25_EIRENE
      integer i
      real (kind=R8) ::
     * fmax,fmay,fmmax,dfx,dfy,uuu,vvv,r,sb,cb,rl,rd,brx,bzy,
     * uu1,vv1,uu2,vv2
      intrinsic max, sqrt

      fmax=earrowx(1,is)
      fmay=earrowy(1,is)
      fmmax=sqrt(fmax**2+fmay**2)
      do i=1,ntrii
        fmax=earrowx(i,is)
        fmay=earrowy(i,is)
        fmmax=max(fmmax,sqrt(fmax*fmax+fmay*fmay))
      enddo
      if(fmmax.eq.0.0_R8) then
        write(*,*) 'Cannot plot arrows for a zero field'
        return
      endif
      if (b2plot_debug.ge.2) write (*,*) 'drawing arrows'
#ifndef PLPLOT
      call lastd
      call gslwsc_d(arroww) ! increase arrow line width
      if (acol .gt. 1) call gsplci(acol) ! change arrow colour
#endif
      do i=1,ntrii
        brx=
     .   (xtrian(necke(1,i))+xtrian(necke(2,i))+xtrian(necke(3,i)))/3.
        bzy=
     .   (ytrian(necke(1,i))+ytrian(necke(2,i))+ytrian(necke(3,i)))/3.
        dfx=earrowx(i,is)/fmmax
        dfy=earrowy(i,is)/fmmax
        uuu=earrowl*dfx
        vvv=earrowl*dfy
        r=sqrt(uuu**2+vvv**2)
        sb=vvv/r
        cb=uuu/r
        rl=0.2_R8*r
        rd=0.1_R8*r
        uu1=uuu-rl*cb+rd*sb
        vv1=vvv-rl*sb-rd*cb
        uu2=uuu-rl*cb-rd*sb
        vv2=vvv-rl*sb+rd*cb
        call lined_d(brx,bzy,brx+uuu,bzy+vvv)
        call lined_d(brx+uuu,bzy+vvv,brx+uu1,bzy+vv1)
        call lined_d(brx+uuu,bzy+vvv,brx+uu2,bzy+vv2)
c       call lined_d(brx+uu1,bzy+vv1,brx+uu2,bzy+vv2)
      enddo
#ifndef PLPLOT
      call gslwsc_d(1.0_R8)
      if (acol .gt. 1) then
        call lastd
        call gsplci(1)
      endif
#endif
      if (b2plot_debug.ge.2) write (*,*) 'finished drawing arrows'
#endif
      return
      end subroutine plot_trig_arrows
