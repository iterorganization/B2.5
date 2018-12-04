      module b2mod_boundary_sources                          !srv 26.02.18 { 
      use b2mod_types
      implicit none
      real (kind=R8), allocatable, save ::
     *  bv_po(:), 
     *  bv_te(:),               
     *  bv_ti(:), 
     *  bv_ua(:,:), 
     *  bv_na(:,:)

      contains

      subroutine alloc_b2mod_boundary_sources(nCv,nsd)
       implicit none
       integer nCv, nsd 
       allocate(bv_po(nCv))
       allocate(bv_te(nCv))
       allocate(bv_ti(nCv))
       allocate(bv_ua(nCv,0:nsd-1))
       allocate(bv_na(nCv,0:nsd-1))

       bv_po = 0.0_R8
       bv_te = 0.0_R8
       bv_ti = 0.0_R8
       bv_ua = 0.0_R8
       bv_na = 0.0_R8

       write(*,*)'Allocated space for b2mod_boundary_sources' 

      end subroutine alloc_b2mod_boundary_sources 

      subroutine dealloc_b2mod_boundary_sources
       implicit none
       
       deallocate(bv_po)
       deallocate(bv_te)
       deallocate(bv_ti)
       deallocate(bv_ua)
       deallocate(bv_na)

      end subroutine dealloc_b2mod_boundary_sources
      
      subroutine read_boundary_sources(nCv, ns)
       use b2mod_types
       use b2mod_boundary_namelist
       implicit none
       integer nCv, ns

       integer ib, is, iCv
       logical ltemp, lexist

       inquire(file='../bv_na_us.dat',exist=lexist)
       call xertst(lexist,'../bv_na_us.dat does not exist')
       open(70,file='../bv_na_us.dat')
       read (70,*) ((bv_na(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       inquire(file='../bv_ua_us.dat',exist=lexist)
       call xertst(lexist,'../bv_ua_us.dat  does not exist')
       open(70,file='../bv_ua_us.dat')
       read (70,*) ((bv_ua(iCv,is),iCv=1,nCv),is=0,ns-1)
       close(70)

       inquire(file='../bv_te_us.dat',exist=lexist)
       call xertst(lexist,'../bv_te_us.dat  does not exist')
       open(70,file='../bv_te_us.dat')
       read (70,*) (bv_te(iCv),iCv=1,nCv)
       close(70)

       inquire(file='../bv_ti_us.dat',exist=lexist)
       call xertst(lexist,'../bv_ti_us.dat  does not exist')
       open(70,file='../bv_ti_us.dat')
       read (70,*) (bv_ti(iCv),iCv=1,nCv)
       close(70)

       inquire(file='../bv_po_us.dat',exist=lexist)
       call xertst(lexist,'../bv_po_us.dat does not exist')
       open(70,file='../bv_po_us.dat')
       read (70,*) (bv_po(iCv),iCv=1,nCv)
       close(70)
 
      end subroutine read_boundary_sources

      end module b2mod_boundary_sources                      !srv 26.02.18 }
