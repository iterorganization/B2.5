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

      end subroutine alloc_b2mod_equation_sources 

      subroutine dealloc_b2mod_equation_sources
       implicit none
       
       deallocate(art_sch)
       deallocate(art_she)
       deallocate(art_shi)
       deallocate(art_smo)
       deallocate(art_sna)

      end subroutine dealloc_b2mod_equation_sources
      
      subroutine read_equation_sources(nCv, ns)
       use b2mod_types
       implicit none
       integer nCv, ns

       integer is, iCv
       logical art_exists

       INQUIRE (FILE = '../art_sna.dat', EXIST = art_exists)
       if(art_exists) then
         open(70,FILE = '../art_sna.dat')
         read (70,"(<nCv>(2x,e23.16))") 
     &       ((art_sna(iCv,0,is),iCv=1,nCv),is=0,ns-1)
         close(70)
       else
         write(*,*) ' There is no file ../art_sna.dat'
         stop
       endif 

       INQUIRE (FILE = '../art_smo.dat', EXIST = art_exists)
       if(art_exists) then
         open(70,FILE = '../art_smo.dat')
         read (70,"(<nCv>(2x,e23.16))") 
     &       ((art_smo(iCv,0,is),iCv=1,nCv),is=0,ns-1)
         close(70)
       else
         write(*,*) ' There is no file ../art_smo.dat'
         stop
       endif 

       INQUIRE (FILE = '../art_she.dat', EXIST = art_exists)
       if(art_exists) then
         open(70,FILE = '../art_she.dat')
         read (70,"(<nCv>(2x,e23.16))") 
     &        (art_she(iCv,0),iCv=1,nCv)
         close(70)
       else
         write(*,*) ' There is no file ../art_she.dat'
         stop
       endif 

       INQUIRE (FILE = '../art_shi.dat', EXIST = art_exists)
       if(art_exists) then
         open(70,FILE = '../art_shi.dat')
         read (70,"(<nCv>(2x,e23.16))") 
     &        (art_shi(iCv,0),iCv=1,nCv)
         close(70)
       else
         write(*,*) ' There is no file ../art_she.dat'
         stop
       endif 

       INQUIRE (FILE = '../art_sch.dat', EXIST = art_exists)
       if(art_exists) then
         open(70,FILE = '../art_sch.dat')
         read (70,"(<nCv>(2x,e23.16))") 
     &        (art_sch(iCv,0),iCv=1,nCv)
         close(70)
       else
         write(*,*) ' There is no file ../art_sch.dat'
         stop
       endif 

      end subroutine read_equation_sources


      end module b2mod_equation_sources                                 !srv 26.02.18 }
