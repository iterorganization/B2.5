      module b2mod_equation_sources                                     !srv 26.02.18 {
      use b2mod_types
      implicit none
       real (kind=R8), allocatable, save ::
     *  art_sch(:,:),
     *  art_she(:,:),
     *  art_shi(:,:),
     *  art_smo(:,:,:),
     *  art_sna(:,:,:)

      contains

      subroutine alloc_b2mod_equation_sources(nCv,nsd)
       implicit none
       integer nCv, nsd

       allocate(art_sch(nCv,0:3))
       allocate(art_she(nCv,0:3))
       allocate(art_shi(nCv,0:3))
       allocate(art_smo(nCv,0:3,0:nsd-1))
       allocate(art_sna(nCv,0:1,0:nsd-1))

       art_sch =0.0_R8
       art_she =0.0_R8
       art_shi =0.0_R8
       art_smo =0.0_R8
       art_sna =0.0_R8

       write(*,*)'Allocated space for b2mod_equation_sources'
       return

      end subroutine alloc_b2mod_equation_sources

      subroutine dealloc_b2mod_equation_sources
       implicit none

       deallocate(art_sch)
       deallocate(art_she)
       deallocate(art_shi)
       deallocate(art_smo)
       deallocate(art_sna)

       return
      end subroutine dealloc_b2mod_equation_sources

      subroutine read_equation_sources(nCv, ns, switch, itim)
       use b2mod_types
       use b2mod_switches
       implicit none
       type (switches), intent (inout) :: switch
       integer nCv, ns
       integer is, iCv, i, itim
       logical b2news_solving(4)
       character*27 filename
       logical art_exists
       character*4 chitim

       b2news_solving(1:4)=switch%no_solve.eq.0
       do i=1,4
         b2news_solving(i)=b2news_solving(i).or.
     &    (switch%no_solve.lt.0.and.
     &     mod(abs(switch%no_solve),2**i).ge.2**(i-1))
       enddo

       if (b2news_solving(2)) then
         if (switch%use_mms.eq.1) then
           filename = 'art_sna_us.dat'
           call find_file(filename,art_exists)
           if(art_exists) then
             open(70,FILE=trim(filename))
             read (70,"(1x,e23.16)")
     &        ((art_sna(iCv,0,is),iCv=1,nCv),is=0,ns-1)
             close(70)
           else
             call xerrab ('There is no art_sna_us.dat file')
           endif
         elseif (switch%use_mms.eq.2) then
           write (chitim,'(i4.4)') itim
           filename = 'art_sna_us'//chitim//'.dat'
           call find_file(filename,art_exists)
           if(art_exists) then
             open(70,FILE=trim(filename))
             read (70,"(1x,e23.16)")
     &        ((art_sna(iCv,0,is),iCv=1,nCv),is=0,ns-1)
             close(70)
           else
             call xerrab ('There is no art_sna_us.dat file')
           endif
         endif
       endif

       if (b2news_solving(1)) then
         filename='art_smo_us.dat'
         call find_file(filename,art_exists)
         if(art_exists) then
           open(70,FILE=trim(filename))
           read (70,"(1x,e23.16)")
     &      ((art_smo(iCv,0,is),iCv=1,nCv),is=0,ns-1)
           close(70)
         else
           call xerrab('There is no art_smo_us.dat file')
         endif
       endif

       if (b2news_solving(4)) then
         filename='art_she_us.dat'
         call find_file(filename,art_exists)
         if(art_exists) then
           open(70,FILE=trim(filename))
           read (70,"(1x,e23.16)")
     &      (art_she(iCv,0),iCv=1,nCv)
           close(70)
         else
           call xerrab('There is no art_she_us.dat file')
         endif
         filename='art_shi_us.dat'
         call find_file(filename,art_exists)
         if(art_exists) then
           open(70,FILE=trim(filename))
           read (70,"(1x,e23.16)")
     &      (art_shi(iCv,0),iCv=1,nCv)
           close(70)
         else
           call xerrab('There is no art_shi_us.dat file')
         endif
       endif

       if (b2news_solving(3)) then
         filename='art_sch_us.dat'
         call find_file(filename,art_exists)
         if(art_exists) then
           open(70,FILE=trim(filename))
           read (70,"(1x,e23.16)")
     &      (art_sch(iCv,0),iCv=1,nCv)
           close(70)
         else
           call xerrab('There is no art_sch_us.dat file')
         endif
       endif

       return
      end subroutine read_equation_sources


      end module b2mod_equation_sources                                 !srv 26.02.18 }
