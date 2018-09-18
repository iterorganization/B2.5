c
      subroutine pltoutline(physical)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_indirect
      implicit none
      integer ix,iy
      logical physical

c       call gslwsc_d(1.0_R8)
#ifndef PLPLOT
      call lastd
      if(ocol .gt. 1) call gsplci(ocol)
      call gslwsc_d(owidth)
#endif
      if(physical) then ! mesh physical domain
        do ix=-1,nx
          do iy=-1,ny
            if(leftix(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     >       call lined_d(crx(ix,iy,1),cry(ix,iy,1),
     1                    crx(ix,iy,3),cry(ix,iy,3))
            if(topiy(ix,iy).eq.ny+1 .and. region(ix,iy,0).ne.0 .and.
     &         leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >       call lined_d(crx(ix,iy,0),cry(ix,iy,0),
     1                    crx(ix,iy,1),cry(ix,iy,1))
            if(rightix(ix,iy).eq.nx+1 .and. region(ix,iy,0).ne.0 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     >       call lined_d(crx(ix,iy,0),cry(ix,iy,0),
     1                    crx(ix,iy,2),cry(ix,iy,2))
            if(bottomiy(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &         leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >       call lined_d(crx(ix,iy,2),cry(ix,iy,2),
     1                    crx(ix,iy,3),cry(ix,iy,3))
          enddo
        enddo
      else
        do ix=-1,nx
          do iy=-1,ny
            if(leftix(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     >       call lined_d(xx(ix,1)+1.0_R8,yy(iy,1),
     1                    xx(ix,1)+1.0_R8,yy(iy,1)+1.0_R8)
            if(topiy(ix,iy).eq.ny+1 .and. region(ix,iy,0).ne.0 .and.
     &         leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >       call lined_d(xx(ix,1)       ,yy(iy,1),
     1                    xx(ix,1)+1.0_R8,yy(iy,1))
            if(rightix(ix,iy).eq.nx+1 .and. region(ix,iy,0).ne.0 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1)
     >       call lined_d(xx(ix,1),yy(iy,1),
     1                    xx(ix,1),yy(iy,1)+1.0_R8)
            if(bottomiy(ix,iy).eq.-2 .and. region(ix,iy,0).ne.0 .and.
     &         leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)
     >       call lined_d(xx(ix,1)       ,yy(iy,1)+1.0_R8,
     1                    xx(ix,1)+1.0_R8,yy(iy,1)+1.0_R8)
          enddo
        enddo
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
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_indirect
      implicit none
      integer is,ix,iy
      real (kind=R8) ::
     * fmax,fmay,fmmax,rmagst,rmaguv,r,s,t,u,v,sb,cb,rl,rd,
     1 uu1,vv1,uu2,vv2,brx,bzy,dfx,dfy,uuu,vvv
      intrinsic max, sqrt

      fmax=arrowx(0,0,is)
      fmay=arrowy(0,0,is)
      fmmax=sqrt(fmax**2+fmay**2)
      do iy=-1,ny
        do ix=-1,nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2) then
            fmax=arrowx(ix,iy,is)
            fmay=arrowy(ix,iy,is)
            fmmax=max(fmmax,sqrt(fmax*fmax+fmay*fmay))
          endif
        enddo
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
      do iy=-1,ny
        do ix=-1,nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2) then
            brx=cr(ix,iy)
            bzy=cz(ix,iy)
            s=cr(rightix(ix,iy),rightiy(ix,iy))-
     .        cr(leftix(ix,iy),leftiy(ix,iy))
            t=cz(rightix(ix,iy),rightiy(ix,iy))-
     .        cz(leftix(ix,iy),leftiy(ix,iy))
            u=cr(topix(ix,iy),topiy(ix,iy))-
     .        cr(bottomix(ix,iy),bottomiy(ix,iy))
            v=cz(topix(ix,iy),topiy(ix,iy))-
     .        cz(bottomix(ix,iy),bottomiy(ix,iy))
            rmagst=sqrt(s**2+t**2)
            rmaguv=sqrt(u**2+v**2)

            dfx=arrowx(ix,iy,is)/fmmax
            dfy=arrowy(ix,iy,is)/fmmax

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
c           call lined_d(brx+uu1,bzy+vv1,brx+uu2,bzy+vv2)
          endif
        enddo
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
