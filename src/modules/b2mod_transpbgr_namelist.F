
      module b2mod_transpbgr_namelist
      use b2mod_types

      real (kind=R8), save ::
     *   fix_dna, fix_hci, fix_hce

      contains

      subroutine read_b2mod_transpbgr_namelist

      namelist /transpbgr/ fix_dna,fix_hci,fix_hce
      character*256 filename
      logical :: file_ok

*   ..default: no background transport
      fix_dna=0._R8
      fix_hci=0._R8
      fix_hce=0._R8

      filename='b2.transpbgr.parameters'
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist transpbgr'
#ifdef USE_MPI
        read(99,transpbgr,err=999)
#else
        read(99,transpbgr)
#endif
        write(*,*) 'Read namelist transpbgr'
      else
        write(*,*) 'No b2.transpbgr.parameters file found'
      endif
      return

  999 call xerrab ('Error while reading TRANSPBGR namelist!')
      end subroutine read_b2mod_transpbgr_namelist

      subroutine write_b2mod_transpbgr_namelist

      namelist /transpbgr/  fix_dna,fix_hci,fix_hce

      write(*,transpbgr)

      end subroutine write_b2mod_transpbgr_namelist

      end module b2mod_transpbgr_namelist
