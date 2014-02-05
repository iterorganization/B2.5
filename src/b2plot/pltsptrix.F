c
      subroutine pltsptrix(physical)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_boundary_namelist
      implicit none
      logical physical
      real (kind=R8) ::
     *   x1,x2,y1,y2
      integer ix,iy,ic,ix1,iiso
      logical iniso

#ifndef PLPLOT
      call lastd
      if(scol .gt. 1) call gsplci(scol)
      call gslwsc_d(sepwidth)
#endif
      if(.not.physical) then
        do iy = -1, ny-1
          ix1 = -1
          ix = 0
          do while (ix.lt.nx)
            if(rightix(ix,iy).ne.nx+1 .and.
     &         topix(rightix(ix,iy),rightiy(ix,iy)).ne.
     &	       rightix(topix(ix,iy),topiy(ix,iy))) then
              ix1 = ix
              x1=ix1-0.5_R8
              y1=iy +0.5_R8
              x2=ix1+0.5_R8
              y2=iy +0.5_R8
              call lined_d(x1,y1,x2,y2)
              do while(leftix(ix1,iy).gt.-2 .and.
     &                 leftix(topix(ix1,iy),topiy(ix1,iy)).eq.ix1-1)
                ix1 = ix1-1
                iniso = .false.
                do iiso = 1, nniso
                  iniso = iniso .or. 
     1             (ix1.ge.nxiso1(iiso).and.ix1.le.nxiso2(iiso).and.
     2		    iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
                end do
                if(.not.iniso) then
                  x1=ix1-0.5_R8
                  y1=iy +0.5_R8
                  x2=ix1+0.5_R8
                  y2=iy +0.5_R8
                  call lined_d(x1,y1,x2,y2)
        	end if
              end do
              ix1=ix+1
              x1=ix1-0.5_R8
              y1=iy +0.5_R8
              x2=ix1+0.5_R8
              y2=iy +0.5_R8
              call lined_d(x1,y1,x2,y2)
              do while(rightix(ix1,iy).lt.nx+1 .and.
     &                 rightix(topix(ix1,iy),topiy(ix1,iy)).eq.ix1+1)
                ix1=ix1+1
                iniso = .false.
                do iiso = 1, nniso
                  iniso = iniso .or. 
     1             (ix1.ge.nxiso1(iiso).and.ix1.le.nxiso2(iiso).and.
     2		    iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
                end do
                if(.not.iniso) then
                  x1=ix1-0.5_R8
                  y1=iy +0.5_R8
                  x2=ix1+0.5_R8
                  y2=iy +0.5_R8
                  call lined_d(x1,y1,x2,y2)
                end if
              end do                
              ix1=-1                
            end if
            if (ix1.ne.-1 .and. rightix(ix1,iy).ne.nx+1) then
              ix=ix1+1
            else
              ix=ix+1
            end if
          end do
        end do
        do ix = -1, nx
          do iy = -1, ny
            if(rightix(ix,iy).ne.ix+1) then
              x1=ix+0.5_R8
              y1=iy-0.5_R8
              x2=ix+0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(leftix(ix,iy).ne.ix-1) then
              x1=ix-0.5_R8
              y1=iy-0.5_R8
              x2=ix-0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(topiy(ix,iy).ne.iy+1) then
              x1=ix-0.5_R8
              y1=iy+0.5_R8
              x2=ix+0.5_R8
              y2=iy+0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
            if(bottomiy(ix,iy).ne.iy-1) then
              x1=ix-0.5_R8
              y1=iy-0.5_R8
              x2=ix+0.5_R8
              y2=iy-0.5_R8
              call lined_d(x1,y1,x2,y2)
            end if
          end do
        end do
        do ic=1,nniso
          x1=nxiso1(ic)-0.5_R8
          y1=nyiso1(ic)-0.5_R8
          x2=nxiso1(ic)-0.5_R8
          y2=nyiso2(ic)+0.5_R8
          call lined_d(x1,y1,x2,y2)
          x1=nxiso2(ic)+0.5_R8
          y1=nyiso1(ic)-0.5_R8
          x2=nxiso2(ic)+0.5_R8
          y2=nyiso2(ic)+0.5_R8
          call lined_d(x1,y1,x2,y2)
          x1=nxiso1(ic)-0.5_R8
          y1=nyiso1(ic)-0.5_R8
          x2=nxiso2(ic)+0.5_R8
          y2=nyiso1(ic)-0.5_R8
          call lined_d(x1,y1,x2,y2)
          x1=nxiso1(ic)-0.5_R8
          y1=nyiso2(ic)+0.5_R8
          x2=nxiso2(ic)+0.5_R8
          y2=nyiso2(ic)+0.5_R8
          call lined_d(x1,y1,x2,y2)
        enddo
      endif
      if(physical) then
        do iy = -1, ny-1
          ix1 = -1
          ix = 0
          do while (ix.lt.nx)
            if(rightix(ix,iy).ne.nx+1 .and.
     &         topix(rightix(ix,iy),rightiy(ix,iy)).ne.
     &	       rightix(topix(ix,iy),topiy(ix,iy))) then
              ix1 = ix
              call lined_d(crx(ix1,iy,3),cry(ix1,iy,3),
     1	       crx(ix1,iy,2),cry(ix1,iy,2))
              do while(leftix(ix1,iy).gt.-2 .and.
     &           leftix(topix(ix1,iy),topiy(ix1,iy)).eq.ix1-1)
                ix1 = ix1-1
                iniso = .false.
                do iiso = 1, nniso
                  iniso = iniso .or. 
     1             (ix1.ge.nxiso1(iiso).and.ix1.le.nxiso2(iiso).and.
     2		    iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
                end do
                if(.not.iniso)
     1		 call lined_d(crx(ix1,iy,3),cry(ix1,iy,3),
     2		  crx(ix1,iy,2),cry(ix1,iy,2))
              end do
              ix1=ix+1
              call lined_d(crx(ix1,iy,3),cry(ix1,iy,3),
     1	       crx(ix1,iy,2),cry(ix1,iy,2))
              do while(rightix(ix1,iy).lt.nx+1 .and.
     &                 rightix(topix(ix1,iy),topiy(ix1,iy)).eq.ix1+1)
                ix1=ix1+1
                iniso = .false.
                do iiso = 1, nniso
                  iniso = iniso .or. 
     1             (ix1.ge.nxiso1(iiso).and.ix1.le.nxiso2(iiso).and.
     2		    iy .ge.nyiso1(iiso).and.iy .le.nyiso2(iiso))
                end do
                if(.not.iniso)
     1		 call lined_d(crx(ix1,iy,3),cry(ix1,iy,3),
     2		  crx(ix1,iy,2),cry(ix1,iy,2))
              end do                
              ix1=-1                
            end if
            if (ix1.ne.-1 .and. rightix(ix1,iy).ne.nx+1) then
              ix=ix1+1
            else
              ix=ix+1
            end if
          end do
        end do
      endif
#ifndef PLPLOT
      call gslwsc_d(lwidth)
      call lastd
      if(scol .gt. 1) call gsplci(1)
#endif
      return
      end
