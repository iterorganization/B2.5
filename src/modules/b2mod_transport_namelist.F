      module b2mod_transport_namelist
      use b2mod_types
      use b2mod_indirect
      implicit none

      integer, save ::
     *   flag_dna,
     *   flag_dpa,
     *   flag_vla,
     *   flag_vsa,
     *   flag_hci,
     *   flag_hce,
     *   flag_sig,
     *   flag_alf
      real (kind=R8), save ::
     *   parm_hce,
     *   parm_sig,
     *   parm_alf
#ifdef ALLOCATE_AND_NAMELIST
      real (kind=R8), allocatable, save ::
     *   parm_dna(:), 
     *   parm_dpa(:), 
     *   parm_vla(:), 
     *   parm_vsa(:), 
     *   parm_hci(:)
#else
#     include "DIMENSIONS.F"
      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)
      real (kind=R8), save :: parm_dna(0:NSDMAX-1)
      real (kind=R8), save :: parm_dpa(0:NSDMAX-1)
      real (kind=R8), save :: parm_vla(0:NSDMAX-1)
      real (kind=R8), save :: parm_vsa(0:NSDMAX-1)
      real (kind=R8), save :: parm_hci(0:NSDMAX-1)
#endif
      character*256,  save :: transport_filename
      real (kind=R8), save :: transport_time_mod
      real (kind=R8), save :: transport_time_switch
      real (kind=R8), save :: cflme, cflmi, cflmv
      real (kind=R8), save :: vout_cnv(0:NSDMAX-1), 
     .  pw0_cnv(0:NSDMAX-1), pw1_cnv(0:NSDMAX-1), pw2_cnv(0:NSDMAX-1)

c*** vout_cnv  :  value of the "blob" convection velocity at the outer 
c***              grid edge
c*** pw0_cnv   :  exponent in a radial profile of the "blob" velocity
c*** pw[12]_cnv:  exponents in a poloidal profile of the "blob" velocity

      logical, save :: allocated_transport_namelist = .false.
      private :: transport
      logical, save :: write_nml_transp = .true.
      namelist /transport/ 
     1   flag_dna,flag_dpa,flag_vla,flag_vsa,flag_hci,
     2   flag_hce,flag_sig,flag_alf,
     3   parm_dna,parm_dpa,parm_vla,parm_vsa,parm_hci,
     4   parm_hce,parm_sig,parm_alf,
     5   transport_filename,transport_time_mod,transport_time_switch,
     6   cflme,cflmi,cflmv,
     7   vout_cnv,pw0_cnv,pw1_cnv,pw2_cnv, write_nml_transp

      
      contains

      subroutine alloc_transport_namelist(nsd)
      implicit none
      integer nsd
      if (allocated_transport_namelist) return
#ifdef ALLOCATE_AND_NAMELIST
      allocate(parm_dna(0:nsd-1))
      allocate(parm_dpa(0:nsd-1))
      allocate(parm_vla(0:nsd-1))
      allocate(parm_vsa(0:nsd-1))
      allocate(parm_hci(0:nsd-1))
#else
      call xertst (nsd.le.nsdmax,
     1  'increase NSDMAX in b2mod_transport_namelist')
#endif
      flag_dna=0
      flag_dpa=0
      flag_vla=0
      flag_vsa=0
      flag_hci=0
      flag_hce=0
      flag_sig=0
      flag_alf=0
      transport_filename='b2.transport.parameters'
      transport_time_mod=0.0_R8
      transport_time_switch=0.0_R8
      vout_cnv=0.0_R8 
      pw0_cnv=1.0_R8 
      pw1_cnv=5.0_R8 
      pw2_cnv=2.0_R8
      allocated_transport_namelist = .true. 
      return
      end subroutine alloc_transport_namelist

      subroutine dealloc_transport_namelist
      implicit none
      if (.not.allocated_transport_namelist) return
#ifdef ALLOCATE_AND_NAMELIST
      deallocate(parm_dna)
      deallocate(parm_dpa)
      deallocate(parm_vla)
      deallocate(parm_vsa)
      deallocate(parm_hci)
#endif
      allocated_transport_namelist = .false.
      return
      end subroutine dealloc_transport_namelist

      subroutine read_b2mod_transport_namelist

      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none

      logical :: file_ok
      character*260 filename

      filename='../' // transport_filename
      inquire(file=filename,exist=file_ok)
      if(file_ok) then
        cflme=cflim(0)
        cflmi=cflim(1)
        cflmv=cflim(2)
        open(99,file=filename)
        write(*,*) 'Reading namelist transport from ',trim(filename)
#ifdef USE_MPI
        read(99,transport,err=999)
#else
        read(99,transport)
#endif
        write(*,*) 'Read namelist transport from ',trim(filename)
        cflim(0)=cflme
        cflim(1)=cflmi
        cflim(2)=cflmv
      else
        write(*,*) 'No ',trim(filename)
        cflme=-1.0_R8
        cflmi=-1.0_R8
        cflmv=-1.0_R8
      endif

      return
 999  call xerrab ('Error while reading TRANSPORT namelist!')
      end subroutine read_b2mod_transport_namelist

      subroutine write_b2mod_transport_namelist

      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none 

      cflme=cflim(0)
      cflmi=cflim(1)
      cflmv=cflim(2)
      if (write_nml_transp) write(*,transport)

      return
      end subroutine write_b2mod_transport_namelist

      subroutine transform_transport(flag,parm,cf)
      implicit none
      integer flag
      real (kind=R8) ::
     *   parm,cf(0:7)

      if(flag.eq.0) then           ! use Bas
      elseif(flag.eq.1) then       ! constant
        cf(:)=0.0_R8
        cf(0)=parm
      elseif(flag.eq.2) then       ! 1/N
        cf(:)=0.0_R8
        cf(1)=parm
      elseif(flag.eq.3) then       ! Bohm
        cf(:)=0.0_R8
        cf(3)=parm
      elseif(flag.eq.4) then       ! Flux-scaled
        cf(:)=0.0_R8
        cf(0)=parm
        cf(7)=1.0_R8
      else
        write(*,*) 'Option ',flag,' not implemented'
      endif

      return
      end subroutine transform_transport

      end module b2mod_transport_namelist
