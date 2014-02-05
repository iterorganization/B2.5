      module b2mod_first_flight
      use b2mod_types
      implicit none

      integer, save :: no_of_slices, size_of_table, no_of_start_points
      real (kind=R8), allocatable, save ::
     *  angle(:)
      integer, allocatable, save :: 
     *     cell(:,:),cellx(:),celly(:)
      integer, allocatable, save :: table(:,:),index(:,:,:)
      integer, save :: ptr_table, ptr_index
      real (kind=R8), save :: dl, eion

      contains

      subroutine alloc_b2mod_first_flight(in_no_of_slices,
     *     in_size_of_table,in_no_of_start_points,nx,ny,in_dl,in_eion)
      use b2mod_constants
      implicit none
      integer in_no_of_slices,in_size_of_table,
     *     in_no_of_start_points,nx,ny
      real (kind=R8) ::
     *     t1,t2, in_dl, in_eion
      integer ia,ic,ix,iy
      no_of_slices=in_no_of_slices
      size_of_table=in_size_of_table
      no_of_start_points=in_no_of_start_points
      allocate(angle(no_of_slices))
      t1=-pi/2.0_R8
      do ia=1,no_of_slices
         t2=asin(real(ia)/real(no_of_slices)*2+sin(-pi/2.0_R8))
         angle(ia)=0.5_R8*(t1+t2)
         t1=t2
         write(*,*) angle(ia)/pi*180
      enddo
      allocate(cell(-1:nx,-1:ny))
      allocate(cellx((nx+2)*(ny+2)))
      allocate(celly((nx+2)*(ny+2)))
      ic=0
      do ix=-1,nx
         do iy=-1,ny
            ic=ic+1
            cell(ix,iy)=ic
            cellx(ic)=ix
            celly(ic)=iy
         enddo
      enddo
      allocate(table(size_of_table,0:1))
      ptr_table=0
      allocate(index(no_of_start_points,no_of_slices,0:1))
      ptr_index=0
      dl=in_dl
      eion=in_eion
      end subroutine alloc_b2mod_first_flight

      subroutine setup_b2mod_first_flight
     * (ix0,iy0,nx,ny,orientation,output)
      use b2mod_geo
      use b2mod_plasma
      use b2mod_constants
      implicit none
      integer ix0,iy0,nx,ny,output
      integer ix,iy,ia,curr,cnt,new

      real (kind=R8) ::
     *    t1,x,y,dx,dy
      logical inside
      character*1 orientation

      ptr_index=ptr_index+1
      if(ptr_index.gt.no_of_start_points) 
     *     stop 'b2mod_first_flight: increase no_of_start_points'
      write(*,*) 'setup_b2mod_first_flight: ',
     *     ix0,iy0,orientation,ptr_index,ptr_table
      do ia=1,no_of_slices
         index(ptr_index,ia,0)=0
         index(ptr_index,ia,1)=-1
         cnt=0
         curr=-1
         new=-1
         ix=ix0
         iy=iy0
         if(orientation.eq.'W') then
            t1=atan2(cry(ix,iy,2)-cry(ix,iy,0),
     *           crx(ix,iy,2)-crx(ix,iy,0))-
     *           pi/2e0_R8+angle(ia)
            x=0.5_R8*(crx(ix,iy,2)+crx(ix,iy,0))
            y=0.5_R8*(cry(ix,iy,2)+cry(ix,iy,0))
         elseif(orientation.eq.'E') then
            t1=atan2(cry(ix,iy,1)-cry(ix,iy,3),
     *           crx(ix,iy,1)-crx(ix,iy,3))-
     *           pi/2e0_R8+angle(ia)
            x=0.5_R8*(crx(ix,iy,3)+crx(ix,iy,1))
            y=0.5_R8*(cry(ix,iy,3)+cry(ix,iy,1))
         elseif(orientation.eq.'N') then
            t1=atan2(cry(ix,iy,3)-cry(ix,iy,2),
     *           crx(ix,iy,3)-crx(ix,iy,2))-
     *           pi/2e0_R8+angle(ia)
            x=0.5_R8*(crx(ix,iy,3)+crx(ix,iy,2))
            y=0.5_R8*(cry(ix,iy,3)+cry(ix,iy,2))
         elseif(orientation.eq.'S') then
            t1=atan2(cry(ix,iy,0)-cry(ix,iy,1),
     *           crx(ix,iy,0)-crx(ix,iy,1))-
     *           pi/2e0_R8+angle(ia)
            x=0.5_R8*(crx(ix,iy,1)+crx(ix,iy,0))
            y=0.5_R8*(cry(ix,iy,1)+cry(ix,iy,0))
         else
            write(*,*) 'Unknown orientation ',orientation,' at ',ix0,iy0
            stop
         endif
         inside=.true.
         if (output.ge.1) write(88,*) x,y,ia,ix,iy
         dx=dl*cos(t1)
         dy=dl*sin(t1)
         do while(inside)
            x=x+dx
            y=y+dy
            call chord_find(nx,ny,crx,cry,x,0.0_R8,y,inside,ix,iy)
            if(inside) then
               new=cell(ix,iy)
               if(new.eq.curr) then
                  cnt=cnt+1
               else
                  if(curr.eq.-1) then
                     index(ptr_index,ia,0)=ptr_table+1
                  else
                     ptr_table=ptr_table+1
                     if(ptr_table.gt.size_of_table) 
     *        stop 'b2mod_first_flight: increase size_of_table'
                     table(ptr_table,0)=curr
                     table(ptr_table,1)=cnt
                  endif
                  curr=new
                  cnt=1
               endif
            else
               if(cnt.gt.0) then
                  ptr_table=ptr_table+1
                  table(ptr_table,0)=curr
                  table(ptr_table,1)=cnt
                  index(ptr_index,ia,1)=ptr_table
               endif
            endif
            if (output.ge.1) 
     >       write(88,*) x,y,inside,ix,iy,new,cnt,curr,ptr_table
         enddo
         if (output.ge.1) write(88,'()')
      enddo
      return
      end subroutine setup_b2mod_first_flight

      subroutine process_b2mod_first_flight(
     *     indptr,ix0,iy0,nx,ny,ns,nscx,iscx,isneut,src,sri,work,wrk2,
     *     divpuff,mainpuff,corepuff)
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer indptr,ix0,iy0,nx,ny,ns,nscx,iscx(0:nscxmax-1),isneut
      integer ix,iy,ia,ptr,ision,k
      real (kind=R8) ::
     *     dist_start,dist_end,dist_deposit,vcx,
     *     rcx,src,sri,work(-1:nx,-1:ny),wrk2(-1:nx,-1:ny),
     *     divpuff,mainpuff,corepuff
      real (kind=R8) ::
     *   machsfr, expu, cutlo, cutll, t0
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
*     (expu is the exponential function with underflow set to 0)
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
      if((region(ix0,iy0,0).eq.1 .and. nnreg(0).ge.4) .or.
     &   (region(ix0,iy0,0).eq.5 .and. nnreg(0).eq.8)) return   ! ignore core
      dist_start=src/no_of_slices
      k=0
      do while (iscx(k).ne.isneut .and. k.lt.nscx)
        k = k+1
      enddo
      call xertst (k.lt.nscx, 'CX species index not found!')
      do ia=1,no_of_slices
         dist_end=dist_start
         do ptr=index(indptr,ia,0),index(indptr,ia,1)
            ix=cellx(table(ptr,0))
            iy=celly(table(ptr,0))
            vcx=0
            do ision = 0, ns-1 
              rcx=expu(rlcx(ix,iy,0,ision,k)+
     *         rlcx(ix,iy,1,ision,k)*log(ti(ix,iy)/(am(isneut)*ev)))
              vcx=vcx+rcx*na(ix,iy,ision)
            enddo
            dist_deposit=vcx*dl/
     *       sqrt(eion*ev/(am(isneut)*mp))*table(ptr,1)
            dist_deposit=min(dist_end*dist_deposit,dist_end)
            work(ix,iy)=work(ix,iy)+dist_deposit
            if (src.ne.0.0_R8) 
     >       wrk2(ix,iy)=wrk2(ix,iy)+dist_deposit*sri/src
            dist_end=dist_end-dist_deposit
            if(dist_end.le.0) exit
         enddo
         if(dist_end.gt.0) then
            if(iy.gt.ny/2) then
               mainpuff=mainpuff+dist_end
            else
               if((region(ix,iy,0).eq.1 .and. nnreg(0).ge.4) .or.
     &            (region(ix,iy,0).eq.5 .and. nnreg(0).eq.8)) then
                  corepuff=corepuff+dist_end
               else
                  divpuff=divpuff+dist_end
               endif
            endif
         endif
      enddo
      end subroutine process_b2mod_first_flight

      subroutine read_b2mod_first_flight(nread)
      implicit none
      integer nread
      end subroutine read_b2mod_first_flight

      subroutine write_b2mod_first_flight(nwrite)
      implicit none
      integer nwrite
      integer is,ia,ptr
      do is=1,ptr_index
         write(nwrite,*) 'is = ',is
         do ia=1,no_of_slices
            write(nwrite,*) 'index: ',index(is,ia,0),index(is,ia,1)
            do ptr=index(is,ia,0),index(is,ia,1)
               write(nwrite,*) '        table: ',table(ptr,:),'  --  ',
     *              cellx(table(ptr,0)),celly(table(ptr,0))
            enddo
         enddo
      enddo
      end subroutine write_b2mod_first_flight

      end module b2mod_first_flight
