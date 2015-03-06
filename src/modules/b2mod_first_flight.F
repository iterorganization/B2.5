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
      real (kind=R8), save :: dl

      contains

      subroutine alloc_b2mod_first_flight(in_no_of_slices,
     *     in_size_of_table,in_no_of_start_points,nx,ny,in_dl)
      use b2mod_constants
      implicit none
      integer in_no_of_slices,in_size_of_table,
     *     in_no_of_start_points,nx,ny
      real (kind=R8) ::
     *     t1,t2, in_dl
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
      write(6,'(a,2i4,1x,a1,1x,2i6)') 'setup_b2mod_first_flight: ',
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
     *     indptr,ix0,iy0,nx,ny,ns,nscx,iscx,isneut,jsep,
     *     dtim,src,sri,work,wrk2,wrk3,
     *     pfrpuff,mainpuff,corepuff)
      use b2mod_geo
      use b2mod_plasma
      use b2mod_rates
      use b2mod_indirect
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      integer indptr,ix0,iy0,nx,ny,ns,nscx,iscx(0:nscxmax-1),isneut,jsep
      integer ix,iy,ia,ptr,ik,ision,k,ixend,iyend
      real (kind=R8) ::
     *     dtim,dist_start,dist_end,dist_deposit,vcx,vion,vnn,vexc,
     *     rcx,rion,rexc,rnn,src,sri,eion,sxb,
     *     flight_time,flown_time,thru_time,
     *     work(-1:nx,-1:ny,0:1),wrk2(-1:nx,-1:ny),wrk3(-1:nx,-1:ny),
     *     pfrpuff,mainpuff,corepuff
      real (kind=R8) ::
     *   machsfr, expu, cutlo, cutll, t0
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
*     (expu is the exponential function with underflow set to 0)
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
      if((region(ix0,iy0,0).eq.1 .and. nnreg(0).ge.2) .or.
     &   (region(ix0,iy0,0).eq.5 .and. nnreg(0).eq.8)) return   ! ignore core
      if (abs(src).lt.cutlo) return  ! jump over zero sources
      dist_start=src/no_of_slices
      eion=sri/src
      flight_time=dtim*dtco(isneut,region(ix0,iy0,0))*
     .            time_factor(ix0,iy0)
      if (eion.eq.0.0_R8) eion=1.0_R8 ! recover the old default value of 1.0 eV
      k=0
      do while (iscx(k).ne.isneut .and. k.lt.nscx)
        k = k+1
      enddo
      call xertst (k.lt.nscx, 'CX species index not found!')
cxpb
cxpb  We follow the launched neutrals over their pre-computed trajectories
cxpb  and distribute the associated sources along the way.
cxpb  We must account for:
cxpb   a. CX processes with ions
cxpb   b. electron-impact ionization
cxpb   c. electron-impact excitation (we use the S/XB factor to convert 
cxpb      electron heat loss rate into number of collision events)
cxpb   d. elastic collisions with neutrals
cxpb  We ditribute the sources into
cxpb   work(:,:,0:1) => sna(:,:,0,is0:is0+1)
cxpb   wrk2(:,:) => shi(:,:,0)
cxpb   wrk3(:,:) => she(:,:,0)
cxpb  We must also take care that we do not allow the neutrals to travel longer
cxpb  that the timestep permits. We set the average energy of the launched
cxpb  neutrals to be eion=shi/sna.
      do ia=1,no_of_slices
         dist_end=dist_start
         flown_time=0.0_R8
	 ixend=-2
	 iyend=-2
         ptr=index(indptr,ia,0)
         do while (flown_time.lt.flight_time .and. 
     1     ptr.le.index(indptr,ia,1))
            ix=cellx(table(ptr,0))
            iy=celly(table(ptr,0))
            vcx=0.0_R8
            vnn=0.0_R8
            do ision = 0, ns-1 
              rcx=expu(rlcx(ix,iy,0,ision,k)+
     *         rlcx(ix,iy,1,ision,k)*log(ti(ix,iy)/(am(isneut)*ev)))
              vcx=vcx+rcx*na(ix,iy,ision)
            enddo
            rion=expu(rlsa(ix,iy,0,isneut)+
     *         rlsa(ix,iy,1,isneut)*log(te(ix,iy)/ev))
            vion=rion*ne(ix,iy)
            rexc=expu(rlqa(ix,iy,0,isneut)+
     *         rlqa(ix,iy,1,isneut)*log(te(ix,iy)/ev))-
     *         rion*rpt(ix,iy,isneut+1)
            sxb=(expu(rlqa(ix,iy,0,isneut)+
     *         rlqa(ix,iy,1,isneut)*log(te(ix,iy)/ev)))/
     *         (rion*rpt(ix,iy,isneut+1))
            vexc=sxb*vion
            do ik = 0, nscx-1
              rnn=2.0e-19_R8*sqrt(ti(ix,iy)/(mp*am(iscx(ik))))
              vnn=vnn+rnn*na(ix,iy,iscx(ik))
            enddo
            thru_time=dl/sqrt(eion/(am(isneut)*mp))*table(ptr,1)
            if (flown_time+thru_time.lt.flight_time) then
              dist_deposit=(vcx+vion+vexc+vnn)*dl/
     *         sqrt(eion/(am(isneut)*mp))*table(ptr,1)
              dist_deposit=min(dist_end*dist_deposit,dist_end)
            else
              thru_time=flight_time-flown_time
              dist_deposit=dist_end
            endif
            flown_time=flown_time+thru_time
            work(ix,iy,0)=work(ix,iy,0)+dist_deposit*
     *       (vcx+vnn+vexc)/(vcx+vion+vexc+vnn)
            work(ix,iy,1)=work(ix,iy,1)+dist_deposit*
     *       vion/(vcx+vion+vexc+vnn)
            wrk2(ix,iy)=wrk2(ix,iy)+dist_deposit*eion
            wrk3(ix,iy)=wrk3(ix,iy)-dist_deposit*
     *       sxb*rpt(ix,iy,isneut+1)*ev*
     *       (vion+vexc)/(vcx+vion+vexc+vnn) 
            dist_end=dist_end-dist_deposit
            if(dist_end.le.0.0_R8) then
	      ixend = ix
	      iyend = iy
              exit
            endif
            ptr = ptr + 1
         enddo
         if(dist_end.gt.0.0_R8.and.ixend.ne.-2.and.iyend.ne.-2) then
            if(iyend.gt.jsep) then
               mainpuff=mainpuff+dist_end
            else
               if((region(ixend,iyend,0).eq.1 .and. nnreg(0).ge.2) .or.
     &            (region(ixend,iyend,0).eq.5 .and. nnreg(0).eq.8)) then
                  corepuff=corepuff+dist_end
               else if (nnreg(0).ge.4) then
                  pfrpuff=pfrpuff+dist_end
               endif
            endif
         endif
      enddo
      return
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
      return
      end subroutine write_b2mod_first_flight

      end module b2mod_first_flight
