      module b2mod_feedback
      use b2mod_types
      use b2mod_b2cmpa

      integer, parameter :: nvac=4, nvacreg=4
#     include "DIMENSIONS.F"
      integer NSDMAX
      parameter (NSDMAX=DEF_NSD)

*   ..variables
      real (kind=R8), save ::
     * cum_volrec
      logical, save :: feedback_namelist_used
      logical, save :: feedback_read, control_read
      integer, allocatable, save :: feedback_strata(:)

*   ..namelists
      integer, save :: vacuum_communication,
     & vacuum_communication_nreg(nvac), 
     & vacuum_communication_method(nvac),
     & vacuum_communication_iy(nvacreg,nvac),
     & vacuum_communication_ix1(nvacreg,nvac),
     & vacuum_communication_ix2(nvacreg,nvac),
     & na_feedback_choice(0:DEF_NATM-1),
     & na_feedback_option(0:DEF_NATM-1),
     & na_feedback_actuator(0:DEF_NATM-1),
     & na_feedback_ix1(0:DEF_NATM-1), na_feedback_ix2(0:DEF_NATM-1), 
     & na_feedback_iy1(0:DEF_NATM-1), na_feedback_iy2(0:DEF_NATM-1), 
     & na_feedback_ib(0:DEF_NATM-1)
      real (kind=R8), save ::
     * saved_cbshe_core, saved_cbshi_core, saved_cbsna_core,
     * saved_cbsch_core, saved_cbsna_pfr1, saved_cbsna_pfr2,
     * saved_cbsna_sol,
     & vacuum_communication_alpha(0:nsdmax-1,nvac),
     & vacuum_communication_beta(0:nsdmax-1,nvac),
     & na_feedback_target(0:DEF_NATM-1), na_feedback_time(0:DEF_NATM-1),
     & saved_na_feedback_actuator(0:DEF_NATM-1),
     & na_feedback_alpha(0:DEF_NATM-1), na_feedback_beta(0:DEF_NATM-1),
     & na_feedback_const(0:DEF_NATM-1),
     & na_feedback_puff_min(0:DEF_NATM-1),
     & na_feedback_puff_max(0:DEF_NATM-1),
     & na_feedback_overshoot(0:DEF_NATM-1)
      namelist /feedback/ saved_cbshe_core, saved_cbshi_core,
     1 saved_cbsna_core, saved_cbsch_core, saved_cbsna_pfr1, 
     2 saved_cbsna_pfr2, saved_cbsna_sol, saved_na_feedback_actuator
! new namelist for specifying more complicated feedback options
      namelist /feedback_control/ vacuum_communication,
     & vacuum_communication_nreg,
     & vacuum_communication_iy, vacuum_communication_method,
     & vacuum_communication_ix1, vacuum_communication_ix2,
     & vacuum_communication_alpha, vacuum_communication_beta,
     & na_feedback_target, na_feedback_choice, na_feedback_option,
     & na_feedback_actuator, na_feedback_alpha, na_feedback_beta,
     & na_feedback_const, na_feedback_time,
     & na_feedback_ix1, na_feedback_ix2,
     & na_feedback_iy1, na_feedback_iy2, 
     & na_feedback_ib,
     & na_feedback_puff_min, na_feedback_puff_max,
     & na_feedback_overshoot
      data 
     & feedback_namelist_used/.false./, feedback_read/.false./,
     & control_read/.false./,
     1 saved_cbshe_core/0.0_R8/, saved_cbshi_core/0.0_R8/,
     2 saved_cbsna_core/0.0_R8/, saved_cbsch_core/0.0_R8/,
     3 saved_cbsna_pfr1/0.0_R8/, saved_cbsna_pfr2/0.0_R8/,
     4 cum_volrec/0.0_R8/, saved_na_feedback_actuator/DEF_NATM*0.0_R8/,
     5 na_feedback_target/DEF_NATM*0.0_R8/,
     * na_feedback_time/DEF_NATM*0.0_R8/,
     * na_feedback_choice/DEF_NATM*0/,
     6 na_feedback_option/DEF_NATM*0/, na_feedback_actuator/DEF_NATM*0/,
     7 na_feedback_ix1/DEF_NATM*-2/, na_feedback_ix2/DEF_NATM*-2/,
     & na_feedback_iy1/DEF_NATM*-2/, na_feedback_iy2/DEF_NATM*-2/, 
     & na_feedback_ib/DEF_NATM*-1/,
     8 na_feedback_alpha/DEF_NATM*0.001_R8/, 
     & na_feedback_beta/DEF_NATM*1.0_R8/,
     & na_feedback_const/DEF_NATM*0.0_R8/,
     & na_feedback_puff_min/DEF_NATM*0.0_R8/,
     & na_feedback_puff_max/DEF_NATM*0.0_R8/,
     & na_feedback_overshoot/DEF_NATM*0.0_R8/

      contains

      subroutine read_b2mod_feedback(nread,edition,found)
      implicit none
      integer nread
      character*4 edition
      logical found
      character*256 filename
      logical file_ok
      external find_file

      if (edition.eq.'    ') then
        filename = 'b2.feedback_save.parameters'
      else
        filename = 'b2.feedback_save.parameters.'//edition
      endif
      call find_file(filename,file_ok)
      if(file_ok.and.(.not.feedback_read.or.edition.ne.'    ')) then
        open(nread,file=trim(filename))
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,feedback,err=999)
#else
        read(nread,feedback)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
        if (edition.eq.'    ') feedback_read = .true.
      else
        if (.not.file_ok) write(6,*) 'No ',trim(filename)
      endif
      found = max(saved_cbshe_core, saved_cbshi_core, 
     .            saved_cbsna_core, saved_cbsch_core, 
     .            saved_cbsna_pfr1, saved_cbsna_pfr2,
     .            saved_cbsna_sol, maxval(saved_na_feedback_actuator))
     .            .gt. 0.0_R8
      return

#ifdef USE_MPI
  999 call xerrab ('Error while reading FEEDBACK namelist!')
#endif
      end subroutine read_b2mod_feedback

      subroutine write_b2mod_feedback(nwrite,edition)
      implicit none
      integer nwrite
      character*4 edition
      integer is
      character*26 hlp_frm
      character*256 filename

      if (edition.eq.'    ') then
        filename = 'b2.feedback_save.parameters'
        feedback_read = .false.
      else
        filename = 'b2.feedback_save.parameters.'//edition
      endif
cxpb We do the write explicitly to avoid array overflows when the problem
cxpb size is changed in DIMENSIONS.F
      open(nwrite,file=trim(filename))
      write(nwrite,'(a9)') '&FEEDBACK'
      if (saved_cbshe_core.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSHE_CORE=',saved_cbshe_core,','
      if (saved_cbshi_core.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSHI_CORE=',saved_cbshi_core,','
      if (saved_cbsna_core.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSNA_CORE=',saved_cbsna_core,','
      if (saved_cbsch_core.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSCH_CORE=',saved_cbsch_core,','
      if (saved_cbsna_pfr1.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSNA_PFR1=',saved_cbsna_pfr1,','
      if (saved_cbsna_pfr2.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSNA_PFR2=',saved_cbsna_pfr2,','
      if (saved_cbsna_sol.ne.0.0_R8)
     > write(nwrite,*) 'SAVED_CBSNA_SOL= ',saved_cbsna_sol,','
      write(hlp_frm,'(a,i2,a)')
     > '(a,1p,',nspecies,'(1x,1e15.8,'',''))'
      if (maxval(abs(saved_na_feedback_actuator)).gt.0.0_R8)
     > write(nwrite,hlp_frm) 'SAVED_NA_FEEDBACK_ACTUATOR=',
     > (saved_na_feedback_actuator(is),is=0,nspecies-1)
      write(nwrite,'(a1)') '/'
      close(nwrite)
      write(*,*) 'Saved ',trim(filename)
      return

      end subroutine write_b2mod_feedback

      subroutine read_b2mod_feedback_control(nread)
      implicit none
      integer nread
      character*256 filename
      logical file_ok
      external find_file

      if (control_read) return
      vacuum_communication = 0
      vacuum_communication_nreg = 0
      vacuum_communication_method = 0
      vacuum_communication_iy = -2
      vacuum_communication_ix1 = -2
      vacuum_communication_ix2 = -2
      vacuum_communication_alpha = 0.0
      vacuum_communication_beta = 1.0
      filename='b2.feedback_control.parameters'
      call find_file(filename,file_ok)
      if(file_ok) then
        open(nread,file=filename)
        write(6,*) 'Reading ',trim(filename)
#ifdef USE_MPI
        read(nread,feedback_control,err=999)
#else
        read(nread,feedback_control)
#endif
        write(6,*) 'Read ',trim(filename)
        close(nread)
        control_read = .true.
      else
        write(6,*) 'No ',trim(filename)
      endif
      return

#ifdef USE_MPI
  999 call xerrab ('Error while reading FEEDBACK_CONTROL namelist!')
#endif
      end subroutine read_b2mod_feedback_control

      subroutine write_b2mod_feedback_control(nwrite)
      implicit none
      integer nwrite
      open(nwrite,file='b2.feedback_control.parameters')
      write(nwrite,feedback_control)
      close(nwrite)
      control_read = .false.
      end subroutine write_b2mod_feedback_control

      end module b2mod_feedback
