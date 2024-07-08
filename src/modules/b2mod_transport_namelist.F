      module b2mod_transport_namelist
      use b2mod_types
      use b2mod_ad
     & , only : nsdmax
#ifndef ALLOCATE_AND_NAMELIST
      use b2mod_dimensions
#endif
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
      real (kind=R8), save :: parm_dna(0:NSDMAX-1)
      real (kind=R8), save :: parm_dpa(0:NSDMAX-1)
      real (kind=R8), save :: parm_vla(0:NSDMAX-1)
      real (kind=R8), save :: parm_vsa(0:NSDMAX-1)
      real (kind=R8), save :: parm_hci(0:NSDMAX-1)
#endif
      character*256,  save :: transport_filename
      real (kind=R8), save :: transport_time_mod
      real (kind=R8), save :: transport_time_switch
      real (kind=R8), save :: cflme, cflmi, cflmv, cflal, cflab
      real (kind=R8), save :: cflme_core, cflme_sol,
     .                        cflme_tanh_a, cflme_tanh_b
      real (kind=R8), save :: cflmi_core, cflmi_sol,
     .                        cflmi_tanh_a, cflmi_tanh_b
      real (kind=R8), save :: cflmv_core, cflmv_sol,
     .                        cflmv_tanh_a, cflmv_tanh_b
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
     6   cflme,cflmi,cflmv,cflal,cflab,
     7   vout_cnv,pw0_cnv,pw1_cnv,pw2_cnv,write_nml_transp

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

      filename=transport_filename
      call find_file(filename,file_ok)
      if(file_ok) then
        cflme=cflim(0)
        cflmi=cflim(1)
        cflmv=cflim(2)
        cflal=cflim(3)
        cflab=cflim(4)
        open(99,file=filename)
        write(*,*) 'Reading namelist transport from ',trim(filename)
#ifdef USE_MPI
        read(99,transport,err=999)
#else
        read(99,transport)
#endif
        close(99)
        write(*,*) 'Read namelist transport from ',trim(filename)
        cflim(0)=cflme
        cflim(1)=cflmi
        cflim(2)=cflmv
        cflim(3)=cflal
        cflim(4)=cflab
      else
        write(*,*) 'No ',trim(filename)
        cflme=0.0_R8
        cflmi=0.0_R8
        cflmv=0.0_R8
        cflal=0.0_R8
        cflab=0.0_R8
      endif

      return
#ifdef USE_MPI
  999 call xerrab ('Error while reading TRANSPORT namelist!')
#endif
      end subroutine read_b2mod_transport_namelist

      subroutine write_b2mod_transport_namelist

      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none 

      cflme=cflim(0)
      cflmi=cflim(1)
      cflmv=cflim(2)
      cflal=cflim(3)
      cflab=cflim(4)
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

      subroutine cfl_tanh_profile(nCv, geo, mpg, cflm, core, sol, a, b,
     &  cfl)
      use b2us_geo
      use b2us_map
      integer, intent(in) :: nCv
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8), intent(in) :: cflm, core, sol, a, b
      real (kind=R8), intent(inout) :: cfl(nCv)
cxpb  This routine applies a tanh profile, centered at 'a' and of width 'b'
cxpb  along the OMP, to the flux limit multipliers, if the limit values in
cxpb  the deep core ('core') and far-SOL ('sol') are provided
cxpb  The tanh profile applies to the core and CFR.
cpxb  In the PFR, we apply the value at the relevant separatrix.
      integer :: iCv

      cfl(nCv) = cflm
      if (core.gt.0.0_R8.and.sol.gt.0.0_R8) then
        if (core.eq.sol) then
          cfl(nCv) = core
        else
          call xerrab('cfl_tanh_profile not coded yet')
#ifdef WG_TODO
          do iCv=1,nCv
            iFt = mpg%cvFt(iCv)
            if (mpg%cvOnClosedSurface(iCv).or.
     &          mod(mpg%cvReg(iCv,0),4).eq.2) then
              cfl(iCv)=core+(sol-core)*0.5_R8*
     &         (1.0_R8+tanh((geo%ds(iFt)-a)/b))
            else if (nncut.eq.1) then
              if (iy.gt.topcut(1)) then
                cfl(iCv)=core+(sol-core)*0.5_R8*
     &           (1.0_R8+tanh((geo%ds(iFt)-a)/b))
              else
                cfl(iCv)=core+(sol-core)*0.5_R8*
     &           (1.0_R8+tanh((-a)/b))
              end if
            else if (nncut.eq.2) then
              if (mpg%cvReg(iCv,0).eq.3.or.
     &            mpg%cvReg(iCv,0).eq.8) then
                if (iy.gt.topcut(1)) then
                  cfl(iCv)=core+(sol-core)*0.5_R8*
     &             (1.0_R8+tanh((geo%ds(iFt)-a)/b))
                else
                  cfl(iCv)=core+(sol-core)*0.5_R8*
     &             (1.0_R8+tanh((0.5_R8*(geo%ds(topcut(1))+
     &                                   geo%ds(topcut(1)+1))-a)/b))
                end if
              else if (mpg%cvReg(iCv,0).eq.4.or.
     &                 mpg%cvReg(iCv,0).eq.7) then
                if (iy.gt.topcut(2)) then
                  cfl(iCv)=core+(sol-core)*0.5_R8*
     &             (1.0_R8+tanh((geo%ds(iFt)-a)/b))
                else
                  cfl(iCv)=core+(sol-core)*0.5_R8*
     &             (1.0_R8+tanh((0.5_R8*(geo%ds(topcut(2))+
     &                                   geo%ds(topcut(2)+1))-a)/b))
                end if
              end if
            end if
          end do
#endif
        end if
      end if

      return
      end subroutine cfl_tanh_profile

      end module b2mod_transport_namelist
