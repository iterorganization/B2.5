      subroutine clvp(npplvp,npplvm)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot

c version : 04.02.2002 17:40

c=======================================================================
c*** This subroutine produces a set of level values, in ascending order,
c*** for the data array on the top of the stack.
c*** The levels are arranged so that the integral of absolute value of 
c*** the data over the region between each two neighbouring level lines 
c*** is the same.
c*** If npplvm < 0, then this integral is equal (approximately) 
c*** to 1/abs(npplvp) of the integral over the whole volume.
c*** If npplvm.ge.0, then the regions of positive and negative data
c*** values are treated separately, with abs(npplvp) data intervals for
c*** the positive and npplvm for the negative data range. In this case,
c*** if npplvp.ge.0 then 0 is added to the level list, otherwise, the 
c*** biggest negative and smallest positive values fall in one zone.

c*** The data dimensions must be nx, ny, and 1, respectively.

      implicit none
      integer npplvp,npplvm
      integer i,j,k,l,m,n,kk,ll,mm,nn,npoplv
      real (kind=R8) :: tot_hlp,neg_hlp,vol_hlp,d,d_hlp,f,t,v,w_hlp
c=======================================================================
c      print *,'=== clvp: ',npplvp,npplvm !###

c      print *,'=== clvp: iwrk=',iwrk     !###
      npoplv=iabs(npplvp)
      if(iwrk.lt.1) then !{
        print *,'clvp: too few elements on the data stack ,',
     ,                                            '- no data selected?'
        return
      end if !}
c      print *,'=== clvp: nx,ny,ns,ndim',nx,ny,ns,               !###
c     ,                                    (ndim(i,iwrk),i=1,3)  !###
      if(ndim(3,iwrk).ne.1) then !{
        print *,'clvp: ns must be 1 - is',ndim(3,iwrk)
        return
      end if !}
      if(ndim(1,iwrk).ne.nx .or. ndim(2,iwrk).ne.ny) then !{
        print *,'clvp: the data array not conforming to the grid'
        return
      end if !}
      
c*** Put the data in a 1D array and sort it

      n=0
      do j=0,ny-1 !{
        do i=0,nx-1 !{
          n=n+1
          wrk_hlp(n)=work(i,j,0,iwrk)
        end do !}
      end do !}
      call sortda(wrk_hlp,ind_hlp,n)

c*** Calculate the integral over the whole volume

      tot_hlp=0.0_R8
      neg_hlp=0.0_R8
      vol_hlp=0.0_R8
      do j=0,ny-1 !{
        do i=0,nx-1 !{
          tot_hlp=tot_hlp+abs(work(i,j,0,iwrk)*vol(i,j))
          if(work(i,j,0,iwrk).lt.0.0_R8) neg_hlp=neg_hlp-
     -                                        work(i,j,0,iwrk)*vol(i,j)
          vol_hlp=vol_hlp+vol(i,j)
        end do !}
      end do !}
c      print *,'clvp: volume, integral, negative part ',     !###
c     ,                              vol_hlp,tot_hlp,neg_hlp !###
      if(tot_hlp.le.0.0_R8) then !{
        print *,'clvp: no nonzero data found'
        return
      end if !}

c*** Calculate the necessary levels and put them on the real stack

      t=0.0_R8
      v=0.0_R8
      l=0
      mm=1
      nn=npoplv-1
      if(npplvm.gt.0 .and. neg_hlp.gt.0.0_R8) then !{
        nn=nn+npplvm-1
        d=neg_hlp/npplvm
        f=d
        do m=mm,n !{
          if(l.lt.npplvm-1) then !{
            k=ind_hlp(m)
            i=ixx(k)
            j=iyy(k)
            t=t+abs(wrk_hlp(k)*vol(i,j))
            v=v+vol(i,j)
            if(t.ge.f) then !{
              l=l+1
              if(l.eq.npplvm-1) mm=m+1
              wrk_hlp(l)=wrk_hlp(k)
c              print *,'=== clvp: l,m,f,v,level = ',                !###
c     ,                                          l,m,f,v,wrk_hlp(l) !###
              print *,'level',l,wrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/neg_hlp
              f=f+d
            end if !}
          end if !}
        end do !}
      end if !}
      if(npplvm.ge.0) then !{
        i=0
        do m=mm,n !{
          if(i.eq.0) then !{
            k=ind_hlp(m)
            if(wrk_hlp(k).gt.0.0_R8) i=m
          end if !}
        end do !}
        mm=i
        if(npplvp.ge.0) then !{
          nn=nn+1
          l=l+1
          wrk_hlp(l)=0.0_R8
        end if !}
        tot_hlp=tot_hlp-neg_hlp
      end if !}
c      print *,'=== mm,nn,l,tot_hlp,W(mm),W(mm-1)',mm,nn,l,tot_hlp, !###
c     ,                wrk_hlp(ind_hlp(mm)),wrk_hlp(ind_hlp(mm-1))  !###
      if(npoplv.gt.0 .and. tot_hlp.gt.0.0_R8) then !{
        d=tot_hlp/npoplv
        f=d
        k=0
        do m=mm,n !{
          if(l.lt.nn) then !{
            kk=k
            k=ind_hlp(m)
            i=ixx(k)
            j=iyy(k)
            t=t+abs(wrk_hlp(k)*vol(i,j))
            v=v+vol(i,j)
            if(t.ge.f) then !{
              ll=int(t/d)
              if(ll-1.gt.l+1) then !{
                w_hlp=0.0_R8
                if(kk.gt.0) w_hlp=wrk_hlp(kk)
                d_hlp=(wrk_hlp(k)-w_hlp)/(ll-l)
                w_hlp=wrk_hlp(kk)
c                print *,'l,ll,k,kk,w_hlp,d_hlp,wrk_hlp(k)',       !###
c     ,                          l,ll,k,kk,w_hlp,d_hlp,wrk_hlp(k)  !###
                do j=l+1,ll-1 !{
                  l=l+1
                  w_hlp=w_hlp+d_hlp
                  wrk_hlp(l)=w_hlp
c                  print *,'=== clvp: l,m,f,v,level = ',            !###
c     ,                                          l,m,f,v,wrk_hlp(l) !###
                end do !}
              end if !}
              if(l.lt.nn) then !{
                l=l+1
                wrk_hlp(l)=wrk_hlp(k)
c                print *,'=== clvp: l,m,f,v,level = ',              !###
c     ,                                          l,m,f,v,wrk_hlp(l) !###
                print *,'level',l,wrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/tot_hlp
                f=f+d
              end if !}
            end if !}
          end if !}
        end do !}
      end if !}

c*** Transfer the data to the plotting routines

      call clvs(l,wrk_hlp,l)
      return

      contains

      integer function ixx(i)
      implicit none
      integer, intent(in) :: i
      ixx = mod(i-1,nx)+1-1
      return
      end function ixx

      integer function iyy(i)
      implicit none
      integer, intent(in) :: i
      iyy = (i-1)/nx+1-1
      return
      end function iyy
c=======================================================================
      end subroutine clvp

      subroutine eclvp(npplvp,npplvm)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot
#ifdef B25_EIRENE
      use eirmod_ctrig
#endif

c version : 04.02.2002 17:40

c=======================================================================
c*** This subroutine produces a set of level values, in ascending order,
c*** for the data array on the top of the Eirene stack.
c*** The levels are arranged so that the integral of absolute value of 
c*** the data over the region between each two neighbouring level lines 
c*** is the same.
c*** If npplvm < 0, then this integral is equal (approximately) 
c*** to 1/abs(npplvp) of the integral over the whole wall area.
c*** If npplvm.ge.0, then the regions of positive and negative data
c*** values are treated separately, with abs(npplvp) data intervals for
c*** the positive and npplvm for the negative data range. In this case,
c*** if npplvp.ge.0 then 0 is added to the level list, otherwise, the 
c*** biggest negative and smallest positive values fall in one zone.

c*** The data dimensions must be ntrii, and 1, respectively.

      implicit none
      integer npplvp,npplvm
#ifdef B25_EIRENE
      integer i,j,k,l,m,kk,ll,mm,nn,npoplv
      real (kind=R8) :: tot_hlp,neg_hlp,vol_hlp,d,d_hlp,f,t,v,e_hlp
c=======================================================================
c      print *,'=== eclvp: ',npplvp,npplvm !###
c*** Check the data dimensions and input parameter     

c      print *,'=== eclvp: iwrk=',iwrk     !###
      npoplv=iabs(npplvp)
      if(iewrk.lt.1) then !{
        print *,'eclvp: too few elements on the data stack ,',
     ,                                            '- no data selected?'
        return
      end if !}
c      print *,'=== eclvp: ntrii,nedim',ntrii,             !###
c     ,                            (nedim(i,iwwrk),i=1,2)  !###
      if(nedim(2,iewrk).ne.1) then !{
        print *,'eclvp: ns must be 1 - is',nedim(2,iewrk)
        return
      end if !}
      if(nedim(1,iewrk).ne.ntrii) then !{
        print *,'eclvp: the data array not conforming to the triangles'
        return
      end if !}
      
c*** Sort the data

      ewrk_hlp(1:ntrii) = eirene_stack(1:ntrii,0,iewrk)
      call sortda(ewrk_hlp,ind_hlp,ntrii)

c*** Calculate the integral over the whole volume

      tot_hlp=0.0_R8
      neg_hlp=0.0_R8
      vol_hlp=0.0_R8
      do i=1,ntrii !{
        tot_hlp=tot_hlp+abs(eirene_stack(i,0,iewrk)*triangle_vol(i))
        if(eirene_stack(i,0,iewrk).lt.0.0_R8) neg_hlp=neg_hlp-
     -      eirene_stack(i,0,iewrk)*triangle_vol(i)
        vol_hlp=vol_hlp+triangle_vol(i)
      end do !}
c      print *,'eclvp: volume, integral, negative part ',     !###
c     ,                               vol_hlp,tot_hlp,neg_hlp !###
      if(tot_hlp.le.0.0_R8) then !{
        print *,'eclvp: no nonzero data found'
        return
      end if !}

c*** Calculate the necessary levels and put them on the real stack

      t=0.0_R8
      v=0.0_R8
      l=0
      mm=1
      nn=npoplv-1
      if(npplvm.gt.0 .and. neg_hlp.gt.0.0_R8) then !{
        nn=nn+npplvm-1
        d=neg_hlp/npplvm
        f=d
        do m=mm,ntrii !{
          if(l.lt.npplvm-1) then !{
            k=ind_hlp(m)
            t=t+abs(ewrk_hlp(k)*triangle_vol(k))
            v=v+triangle_vol(k)
            if(t.ge.f) then !{
              l=l+1
              if(l.eq.npplvm-1) mm=m+1
              ewrk_hlp(l)=ewrk_hlp(k)
c              print *,'=== eclvp: l,m,f,v,level = ',                !###
c     ,                                          l,m,f,v,ewrk_hlp(l) !###
              print *,'level',l,ewrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/neg_hlp
              f=f+d
            end if !}
          end if !}
        end do !}
      end if !}
      if(npplvm.ge.0) then !{
        i=0
        do m=mm,ntrii !{
          if(i.eq.0) then !{
            k=ind_hlp(m)
            if(ewrk_hlp(k).gt.0.0_R8) i=m
          end if !}
        end do !}
        mm=i
        if(npplvp.ge.0) then !{
          nn=nn+1
          l=l+1
          ewrk_hlp(l)=0.0_R8
        end if !}
        tot_hlp=tot_hlp-neg_hlp
      end if !}
c      print *,'=== mm,nn,l,tot_hlp,E(mm),E(mm-1)',mm,nn,l,tot_hlp, !###
c     ,              ewrk_hlp(ind_hlp(mm)),ewrk_hlp(ind_hlp(mm-1))  !###
      if(npoplv.gt.0 .and. tot_hlp.gt.0.) then !{
        d=tot_hlp/npoplv
        f=d
        k=0
        do m=mm,ntrii !{
          if(l.lt.nn) then !{
            kk=k
            k=ind_hlp(m)
            t=t+abs(ewrk_hlp(k)*triangle_vol(k))
            v=v+triangle_vol(k)
            if(t.ge.f) then !{
              ll=int(t/d)
              if(ll-1.gt.l+1) then !{
                e_hlp=0.0_R8
                if(kk.gt.0) e_hlp=ewrk_hlp(kk)
                d_hlp=(ewrk_hlp(k)-e_hlp)/(ll-l)
                e_hlp=ewrk_hlp(kk)
c                print *,'l,ll,k,kk,e_hlp,d_hlp,ewrk_hlp(k)',       !###
c     ,                          l,ll,k,kk,e_hlp,d_hlp,ewrk_hlp(k)  !###
                do j=l+1,ll-1 !{
                  l=l+1
                  e_hlp=e_hlp+d_hlp
                  ewrk_hlp(l)=e_hlp
c                  print *,'=== eclvp: l,m,f,v,level = ',             !###
c     ,                                          l,m,f,vw,ewrk_hlp(l) !###
                end do !}
              end if !}
              if(l.lt.nn) then !{
                l=l+1
                ewrk_hlp(l)=ewrk_hlp(k)
c                print *,'=== eclvp: l,m,f,v,level = ',               !###
c     ,                                           l,m,f,v,ewrk_hlp(l) !###
                print *,'level',l,ewrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/tot_hlp
                f=f+d
              end if !}
            end if !}
          end if !}
        end do !}
      end if !}

c*** Transfer the data to the plotting routines

      call eclvs(l,ewrk_hlp,l)
#endif
      return
c=======================================================================
      end subroutine eclvp

      subroutine wclvp(npplvp,npplvm)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot

c version : 04.02.2002 17:40

c=======================================================================
c*** This subroutine produces a set of level values, in ascending order,
c*** for the data array on the top of the wall stack.
c*** The levels are arranged so that the integral of absolute value of 
c*** the data over the region between each two neighbouring level lines 
c*** is the same.
c*** If npplvm < 0, then this integral is equal (approximately) 
c*** to 1/abs(npplvp) of the integral over the whole wall area.
c*** If npplvm.ge.0, then the regions of positive and negative data
c*** values are treated separately, with abs(npplvp) data intervals for
c*** the positive and npplvm for the negative data range. In this case,
c*** if npplvp.ge.0 then 0 is added to the level list, otherwise, the 
c*** biggest negative and smallest positive values fall in one zone.

c*** The data dimensions must be nbnd, and 1, respectively.

      implicit none
      integer npplvp,npplvm
      integer i,j,k,l,m,kk,ll,mm,nn,npoplv
      real (kind=R8) :: tot_hlp,neg_hlp,vol_hlp,d,d_hlp,f,t,v,w_hlp
c=======================================================================
c      print *,'=== wclvp: ',npplvp,npplvm !###
c*** Check the data dimensions and input parameter     

c      print *,'=== wclvp: iwrk=',iwrk     !###
      npoplv=iabs(npplvp)
      if(iwwrk.lt.1) then !{
        print *,'wclvp: too few elements on the data stack ,',
     ,                                            '- no data selected?'
        return
      end if !}
c      print *,'=== wclvp: nbnd,nwdim',nbnd,               !###
c     ,                            (nwdim(i,iwwrk),i=1,2)  !###
      if(nwdim(2,iwwrk).ne.1) then !{
        print *,'wclvp: ns must be 1 - is',nwdim(2,iwwrk)
        return
      end if !}
      if(nwdim(1,iwwrk).ne.nbnd) then !{
        print *,'wclvp: the data array not conforming to the grid'
        return
      end if !}
      
c*** Sort the data

      wwrk_hlp(1:nbnd) = wall_stack(1:nbnd,0,iwwrk)
      call sortda(wwrk_hlp,ind_hlp,nbnd)

c*** Calculate the integral over the whole volume

      tot_hlp=0.0_R8
      neg_hlp=0.0_R8
      vol_hlp=0.0_R8
      do i=1,nbnd !{
        tot_hlp=tot_hlp+abs(wall_stack(i,0,iwwrk)*wall_area(i))
        if(wall_stack(i,0,iwwrk).lt.0.0_R8) neg_hlp=neg_hlp-
     -      wall_stack(i,0,iwwrk)*wall_area(i)
        vol_hlp=vol_hlp+wall_area(i)
      end do !}
c      print *,'wclvp: volume, integral, negative part ',     !###
c     ,                               vol_hlp,tot_hlp,neg_hlp !###
      if(tot_hlp.le.0.0_R8) then !{
        print *,'wclvp: no nonzero data found'
        return
      end if !}

c*** Calculate the necessary levels and put them on the real stack

      t=0.0_R8
      v=0.0_R8
      l=0
      mm=1
      nn=npoplv-1
      if(npplvm.gt.0 .and. neg_hlp.gt.0.0_R8) then !{
        nn=nn+npplvm-1
        d=neg_hlp/npplvm
        f=d
        do m=mm,nbnd !{
          if(l.lt.npplvm-1) then !{
            k=ind_hlp(m)
            t=t+abs(wwrk_hlp(k)*wall_area(k))
            v=v+wall_area(k)
            if(t.ge.f) then !{
              l=l+1
              if(l.eq.npplvm-1) mm=m+1
              wwrk_hlp(l)=wwrk_hlp(k)
c              print *,'=== wclvp: l,m,f,v,level = ',                !###
c     ,                                          l,m,f,v,wwrk_hlp(l) !###
              print *,'level',l,wwrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/neg_hlp
              f=f+d
            end if !}
          end if !}
        end do !}
      end if !}
      if(npplvm.ge.0) then !{
        i=0
        do m=mm,nbnd !{
          if(i.eq.0) then !{
            k=ind_hlp(m)
            if(wwrk_hlp(k).gt.0.0_R8) i=m
          end if !}
        end do !}
        mm=i
        if(npplvp.ge.0) then !{
          nn=nn+1
          l=l+1
          wwrk_hlp(l)=0.0_R8
        end if !}
        tot_hlp=tot_hlp-neg_hlp
      end if !}
c      print *,'=== mm,nn,l,tot_hlp,W(mm),W(mm-1)',mm,nn,l,tot_hlp, !###
c     ,              wwrk_hlp(ind_hlp(mm)),wwrk_hlp(ind_hlp(mm-1))  !###
      if(npoplv.gt.0 .and. tot_hlp.gt.0.) then !{
        d=tot_hlp/npoplv
        f=d
        k=0
        do m=mm,nbnd !{
          if(l.lt.nn) then !{
            kk=k
            k=ind_hlp(m)
            t=t+abs(wwrk_hlp(k)*wall_area(k))
            v=v+wall_area(k)
            if(t.ge.f) then !{
              ll=int(t/d)
              if(ll-1.gt.l+1) then !{
                w_hlp=0.0_R8
                if(kk.gt.0) w_hlp=wwrk_hlp(kk)
                d_hlp=(wwrk_hlp(k)-w_hlp)/(ll-l)
                w_hlp=wwrk_hlp(kk)
c                print *,'l,ll,k,kk,w_hlp,d_hlp,wwrk_hlp(k)',       !###
c     ,                          l,ll,k,kk,w_hlp,d_hlp,wwrk_hlp(k)  !###
                do j=l+1,ll-1 !{
                  l=l+1
                  w_hlp=w_hlp+d_hlp
                  wwrk_hlp(l)=w_hlp
c                  print *,'=== wclvp: l,m,f,v,level = ',             !###
c     ,                                          l,m,f,vw,wwrk_hlp(l) !###
                end do !}
              end if !}
              if(l.lt.nn) then !{
                l=l+1
                wwrk_hlp(l)=wwrk_hlp(k)
c                print *,'=== wclvp: l,m,f,v,level = ',               !###
c     ,                                           l,m,f,v,wwrk_hlp(l) !###
                print *,'level',l,wwrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/tot_hlp
                f=f+d
              end if !}
            end if !}
          end if !}
        end do !}
      end if !}

c*** Transfer the data to the plotting routines

      call wclvs(l,wwrk_hlp,l)
      return
c=======================================================================
      end subroutine wclvp

      subroutine tclvp(npplvp,npplvm)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2plot

c version : 04.02.2002 17:40

c=======================================================================
c*** This subroutine produces a set of level values, in ascending order,
c*** for the data array on the top of the target stack.
c*** The levels are arranged so that the integral of absolute value of 
c*** the data over the region between each two neighbouring level lines 
c*** is the same.
c*** If npplvm < 0, then this integral is equal (approximately) 
c*** to 1/abs(npplvp) of the integral over the whole volume.
c*** If npplvm.ge.0, then the regions of positive and negative data
c*** values are treated separately, with abs(npplvp) data intervals for
c*** the positive and npplvm for the negative data range. In this case,
c*** if npplvp.ge.0 then 0 is added to the level list, otherwise, the 
c*** biggest negative and smallest positive values fall in one zone.

c*** The data dimensions must be nbnd and ndepth_nml, respectively.

      implicit none
      integer npplvp,npplvm
      integer i,j,k,l,m,n,kk,ll,mm,nn,npoplv
      real (kind=R8) :: tot_hlp,neg_hlp,vol_hlp,d,d_hlp,f,t,v,w_hlp
c=======================================================================
c      print *,'=== tclvp: ',npplvp,npplvm !###

c      print *,'=== tclvp: iwrk=',iwrk     !###
      npoplv=iabs(npplvp)
      if(itwrk.lt.1) then !{
        print *,'tclvp: too few elements on the data stack ,',
     ,                                            '- no data selected?'
        return
      end if !}
c      print *,'=== tclvp: nbnd,ndepth_nml',nbnd,ndepth_nml               !###
      
c*** Put the data in a 1D array and sort it

      n=0
      do j=1,ndepth_nml !{
        do i=1,nbnd !{
          n=n+1
          wrk_hlp(n)=target_stack(i,j,itwrk)
        end do !}
      end do !}
      call sortda(twrk_hlp,indt_hlp,n)

c*** Calculate the integral over the whole volume

      tot_hlp=0.0_R8
      neg_hlp=0.0_R8
      vol_hlp=0.0_R8
      do j=1,ndepth_nml !{
        do i=1,nwall !{
          tot_hlp=tot_hlp+abs(target_stack(i,j,itwrk)*target_vol(i,j))
          if(target_stack(i,j,itwrk).lt.0.0_R8) neg_hlp=neg_hlp-
     -     target_stack(i,j,itwrk)*target_vol(i,j)
          vol_hlp=vol_hlp+target_vol(i,j)
        end do !}
      end do !}
c      print *,'tclvp: volume, integral, negative part ',     !###
c     ,                               vol_hlp,tot_hlp,neg_hlp !###

      if(tot_hlp.le.0.0_R8) then !{
        print *,'tclvp: no nonzero data found'
        return
      end if !}

c*** Calculate the necessary levels and put them on the real stack

      t=0.0_R8
      v=0.0_R8
      l=0
      mm=1
      nn=npoplv-1
      if(npplvm.gt.0 .and. neg_hlp.gt.0.0_R8) then !{
        nn=nn+npplvm-1
        d=neg_hlp/npplvm
        f=d
        do m=mm,n !{
          if(l.lt.npplvm-1) then !{
            k=indt_hlp(m)
            i=ixx(k)
            j=iyy(k)
            t=t+abs(twrk_hlp(k)*target_vol(i,j))
            v=v+target_vol(i,j)
            if(t.ge.f) then !{
              l=l+1
              if(l.eq.npplvm-1) mm=m+1
              twrk_hlp(l)=twrk_hlp(k)
c              print *,'=== tclvp: l,m,f,v,level = ',                !###
c     ,                                          l,m,f,v,twrk_hlp(l) !###
              print *,'level',l,twrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/neg_hlp
              f=f+d
            end if !}
          end if !}
        end do !}
      end if !}
      if(npplvm.ge.0) then !{
        i=0
        do m=mm,n !{
          if(i.eq.0) then !{
            k=indt_hlp(m)
            if(twrk_hlp(k).gt.0.0_R8) i=m
          end if !}
        end do !}
        mm=i
        if(npplvp.ge.0) then !{
          nn=nn+1
          l=l+1
          twrk_hlp(l)=0.0_R8
        end if !}
        tot_hlp=tot_hlp-neg_hlp
      end if !}
c      print *,'=== mm,nn,l,tot_hlp,W(mm),W(mm-1)',mm,nn,l,tot_hlp, !###
c     ,            twrk_hlp(indt_hlp(mm)),twrk_hlp(indt_hlp(mm-1))  !###
      if(npoplv.gt.0 .and. tot_hlp.gt.0.0_R8) then !{
        d=tot_hlp/npoplv
        f=d
        k=0
        do m=mm,n !{
          if(l.lt.nn) then !{
            kk=k
            k=indt_hlp(m)
            i=ixx(k)
            j=iyy(k)
            t=t+abs(twrk_hlp(k)*target_vol(i,j))
            v=v+target_vol(i,j)
            if(t.ge.f) then !{
              ll=int(t/d)
              if(ll-1.gt.l+1) then !{
                w_hlp=0.0_R8
                if(kk.gt.0) w_hlp=twrk_hlp(kk)
                d_hlp=(twrk_hlp(k)-w_hlp)/(ll-l)
                w_hlp=twrk_hlp(kk)
c                print *,'l,ll,k,kk,w_hlp,d_hlp,twrk_hlp(k)',       !###
c     ,                          l,ll,k,kk,w_hlp,d_hlp,twrk_hlp(k)  !###
                do j=l+1,ll-1 !{
                  l=l+1
                  w_hlp=w_hlp+d_hlp
                  twrk_hlp(l)=w_hlp
c                  print *,'=== tclvp: l,m,f,v,level = ',            !###
c     ,                                         l,m,f,v,twrk_hlp(l) !###
                end do !}
              end if !}
              if(l.lt.nn) then !{
                l=l+1
                twrk_hlp(l)=twrk_hlp(k)
c                print *,'=== tclvp: l,m,f,v,level = ',              !###
c     ,                                          l,m,f,v,twrk_hlp(l) !###
                print *,'level',l,twrk_hlp(k),
     ,                ';  volume & value fractions ',v/vol_hlp,f/tot_hlp
                f=f+d
              end if !}
            end if !}
          end if !}
        end do !}
      end if !}

c*** Transfer the data to the plotting routines

      call tclvs(l,twrk_hlp,l)
      return

      contains

      integer function ixx(i)
      implicit none
      integer, intent(in) :: i
      ixx = mod(i-1,nbnd)+1
      return
      end function ixx

      integer function iyy(i)
      implicit none
      integer, intent(in) :: i
      iyy = (i-1)/nbnd+1
      return
      end function iyy
c=======================================================================
      end subroutine tclvp

