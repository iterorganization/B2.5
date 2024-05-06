      module b2mod_first_flight
      use b2mod_types
      implicit none

      integer, save :: no_of_slices, size_of_table, no_of_start_points
      real (kind=R8), allocatable, save ::
     *  angle(:)
      integer, allocatable, save :: table(:,:), index(:,:,:)
      integer, save :: ptr_table, ptr_index
      real (kind=R8), save :: dl

      contains

      subroutine alloc_b2mod_first_flight(in_no_of_slices,
     *     in_size_of_table, in_no_of_start_points, in_dl)
      use b2mod_constants
      implicit none
      integer in_no_of_slices,in_size_of_table,
     *     in_no_of_start_points
      real (kind=R8) ::
     *     t1,t2, in_dl
      integer ia

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
      allocate(table(size_of_table,0:1))
      ptr_table=0
      allocate(index(no_of_start_points,no_of_slices,0:1))
      ptr_index=0
      dl=in_dl
      end subroutine alloc_b2mod_first_flight

      subroutine dealloc_b2mod_first_flight
      implicit none
      if(.not.allocated(angle)) return
      deallocate(angle)
      deallocate(table)
      deallocate(index)
      return
      end subroutine dealloc_b2mod_first_flight

      subroutine setup_b2mod_first_flight
     *  (mpg,geo,iCv0,iFc0,output)
      use b2us_map
      use b2us_geo
      use b2mod_constants
      use b2mod_subsys
      implicit none
      type(mapping) :: mpg
      type(geometry) :: geo
      integer iCv0,iFc0,output
      integer iCv,ia,curr,cnt,new

      real (kind=R8) ::
     *    t1,x,y,dx,dy
      logical inside
      external xerrab

      call subini ('setup_b2mod_first_flight')
      ptr_index=ptr_index+1
      call xertst(ptr_index.le.no_of_start_points,
     * 'b2mod_first_flight: increase no_of_start_points')
      write(*,'(a,4i6)') 'setup_b2mod_first_flight: ',
     *     iCv0,iFc0,ptr_index,ptr_table
      do ia=1,no_of_slices
         index(ptr_index,ia,0)=0
         index(ptr_index,ia,1)=-1
         cnt=0
         curr=-1
         new=-1
         iCv=iCv0
cxpb We travel the boundary face so that the real cell is to the left
cxpb Recall iCv is the guard cell
         if (mpg%fcCv(iFc0,1).eq.iCv) then
           t1=atan2(geo%cvY(mpg%fcVx(iFc0,1))-geo%cvY(mpg%fcVx(iFc0,2)),
     *              geo%cvX(mpg%fcVx(iFc0,1))-geo%cvX(mpg%fcVx(iFc0,2)))
     *       -pi/2.0_R8+angle(ia)
         else if (mpg%fcCv(iFc0,2).eq.iCv) then
           t1=atan2(geo%cvY(mpg%fcVx(iFc0,2))-geo%cvY(mpg%fcVx(iFc0,1)),
     *              geo%cvX(mpg%fcVx(iFc0,2))-geo%cvX(mpg%fcVx(iFc0,1)))
     *       -pi/2.0_R8+angle(ia)
         end if
         x=0.5_R8*(geo%cvX(mpg%fcVx(iFc0,1))+geo%cvX(mpg%fcVx(iFc0,2)))
         y=0.5_R8*(geo%cvY(mpg%fcVx(iFc0,1))+geo%cvY(mpg%fcVx(iFc0,2)))
         inside=.true.
         if (output.ge.1) write(88,*) x,y,ia,iCv
         dx=dl*cos(t1)
         dy=dl*sin(t1)
         do while(inside)
            x=x+dx
            y=y+dy
            call chord_find(mpg,geo,x,0.0_R8,y,inside,iCv)
            if(inside) then
               new=iCv
               if(new.eq.curr) then
                  cnt=cnt+1
               else
                  if(curr.eq.-1) then
                     index(ptr_index,ia,0)=ptr_table+1
                  else
                     ptr_table=ptr_table+1
                     call xertst(ptr_table.le.size_of_table,
     *                'b2mod_first_flight: increase size_of_table')
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
     >       write(88,*) x,y,inside,iCv,new,cnt,curr,ptr_table
         enddo
         if (output.ge.1) write(88,'()')
      enddo
      call subend ()
      return
      end subroutine setup_b2mod_first_flight

      subroutine process_b2mod_first_flight(
     *     mpg,st,indptr,iCv0,ns,isneut,
     *     dtim,src,sri,work,wrk2,wrk3,
     *     pfrpuff,mainpuff,corepuff)
      use b2us_map
      use b2us_plasma
      use b2mod_math
      use b2mod_numerics_namelist
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
      type(mapping) :: mpg
      type(B2State) :: st
      integer indptr,iCv0,ns,isneut
      integer iCv,ia,ptr,ik,ision,k,iCvend
      real (kind=R8) ::
     *     dtim,dist_start,dist_end,dist_deposit,vcx,vion,vnn,vexc,
     *     rcx,rion,rnn,src,sri,eion,sxb,
     *     flight_time,flown_time,thru_time,
     *     work(mpg%nCv,0:1),wrk2(mpg%nCv),wrk3(mpg%nCv),
     *     pfrpuff,mainpuff,corepuff

      if (mpg%cvOnClosedSurface(iCv0)) return   ! ignore core
      if (.not.b2mod_math_initialised) call b2mod_math_init
      if (abs(src).lt.cutlo) return  ! jump over zero sources
      dist_start=src/no_of_slices
      eion=sri/src
      flight_time=dtim*dtco(isneut,mpg%cvReg(iCv0))*time_factor(iCv0) !WG_TODO!
      if (eion.eq.0.0_R8) eion=1.0_R8 ! recover the old default value of 1.0 eV
      k=0
      do while (st%rt%iscx(k).ne.isneut .and. k.lt.st%rt%nscx)
        k = k+1
      enddo
      if (k.ge.st%rt%nscx) write(*,*) 'CX species index not found!'
cxpb
cxpb  We follow the launched neutrals over their pre-computed trajectories
cxpb  and distribute the associated sources along the way.
cxpb  We must account for:
cxpb   a. CX processes with ions
cxpb   b. electron-impact ionisation
cxpb   c. electron-impact excitation (we use the S/XB factor to convert 
cxpb      electron heat loss rate into number of collision events)
cxpb   d. elastic collisions with neutrals
cxpb  We distribute the sources into
cxpb   work(:,0:1) => sna(:,0,is0:is0+1)
cxpb   wrk2(:) => shi(:,0)
cxpb   wrk3(:) => she(:,0)
cxpb  We must also take care that we do not allow the neutrals to travel longer
cxpb  that the timestep permits. We set the average energy of the launched
cxpb  neutrals to be eion=shi/sna.
      do ia=1,no_of_slices
        dist_end=dist_start
        flown_time=0.0_R8
        iCvend=0
        ptr=index(indptr,ia,0)
        do while (flown_time.lt.flight_time .and.
     1    ptr.le.index(indptr,ia,1))
           iCv=table(ptr,0)
           vcx=0.0_R8
           vnn=0.0_R8
           if (k.lt.st%rt%nscx) then
             do ision = 0, ns-1
               rcx=expu(st%rt%rlcx(iCv,0,ision,k)+
     *                  st%rt%rlcx(iCv,1,ision,k)*
     *              log(st%pl%ti(iCv)/(am(isneut)*ev)))
               vcx=vcx+rcx*st%pl%na(iCv,ision)
             enddo
           end if
           rion=expu(st%rt%rlsa(iCv,0,isneut)+
     *               st%rt%rlsa(iCv,1,isneut)*log(st%pl%te(iCv)/ev))
           vion=rion*st%dv%ne(iCv)
           sxb=(expu(st%rt%rlqa(iCv,0,isneut)+
     *               st%rt%rlqa(iCv,1,isneut)*log(st%pl%te(iCv)/ev)))/
     *         (rion*st%rt%rpt(iCv,isneut+1))
           vexc=sxb*vion
           do ik = 0, st%rt%nscx-1
             rnn=2.0e-19_R8*sqrt(st%pl%tn(iCv)/(mp*am(st%rt%iscx(ik))))
             vnn=vnn+rnn*st%pl%na(iCv,st%rt%iscx(ik))
           enddo
           thru_time=dl/sqrt(eion/(am(isneut)*mp))*table(ptr,1)
           if (flown_time+thru_time.lt.flight_time) then
             dist_deposit=(vcx+vion+vexc+vnn)*dl/
     *        sqrt(eion/(am(isneut)*mp))*table(ptr,1)
             dist_deposit=min(dist_end*dist_deposit,dist_end)
           else
             thru_time=flight_time-flown_time
             dist_deposit=dist_end
           endif
           flown_time=flown_time+thru_time
           work(iCv,0)=work(iCv,0)+dist_deposit*
     *      (vcx+vnn+vexc)/(vcx+vion+vexc+vnn)
           work(iCv,1)=work(iCv,1)+dist_deposit*
     *      vion/(vcx+vion+vexc+vnn)
           wrk2(iCv)=wrk2(iCv)+dist_deposit*eion
           wrk3(iCv)=wrk3(iCv)-dist_deposit*
     *      sxb*st%rt%rpt(iCv,isneut+1)*ev*
     *      (vion+vexc)/(vcx+vion+vexc+vnn)
           dist_end=dist_end-dist_deposit
           if(dist_end.le.0.0_R8) then
             iCvend = iCv
             exit
           endif
           ptr = ptr + 1
         enddo
         if(dist_end.gt.0.0_R8.and.iCvend.ne.0) then
           if(mpg%cvOnClosedSurface(iCvend)) then
             corepuff=corepuff+dist_end
           elseif(mod(mpg%cvReg(iCvend),4).eq.2) then
             mainpuff=mainpuff+dist_end
           else
             pfrpuff=pfrpuff+dist_end
           endif
        endif
      enddo
      return
      end subroutine process_b2mod_first_flight

      subroutine read_b2mod_first_flight(nread)
      implicit none
      integer nread
      call xertst(1.le.nread, 'faulty argument nread')
      return
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
            write(nwrite,*) '        table: ',table(ptr,:)
          enddo
        enddo
      enddo
      return
      end subroutine write_b2mod_first_flight

      end module b2mod_first_flight
