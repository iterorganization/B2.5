      module b2mod_feedback
      use b2mod_types

*   ..variables
      real (kind=R8), save ::
     * cum_volrec

*   ..namelists
      real (kind=R8), save ::
     * saved_cbshe_core, saved_cbshi_core, saved_cbsna_core,
     * saved_cbsch_core, saved_cbsna_pfr1, saved_cbsna_pfr2,
     * saved_cbsna_sol
      namelist /feedback/ saved_cbshe_core, saved_cbshi_core,
     1 saved_cbsna_core, saved_cbsch_core, saved_cbsna_pfr1, 
     2 saved_cbsna_pfr2, saved_cbsna_sol
      
      data 
     1 saved_cbshe_core/0.0_R8/, saved_cbshi_core/0.0_R8/,
     2 saved_cbsna_core/0.0_R8/, saved_cbsch_core/0.0_R8/,
     3 saved_cbsna_pfr1/0.0_R8/, saved_cbsna_pfr2/0.0_R8/,
     4 cum_volrec/0.0_R8/
      contains

      subroutine read_b2mod_feedback(nread)
      integer nread
      logical file_ok
      inquire(file='../b2.feedback_save.parameters',exist=file_ok)
      if(file_ok) then
        open(nread,file='../b2.feedback_save.parameters')
        read(nread,feedback)
        close(nread)
      endif
      end subroutine read_b2mod_feedback

      subroutine write_b2mod_feedback(nwrite)
      integer nwrite
      open(nwrite,file='../b2.feedback_save.parameters')
      write(nwrite,feedback)
      close(nwrite)
      end subroutine write_b2mod_feedback

      end module b2mod_feedback
