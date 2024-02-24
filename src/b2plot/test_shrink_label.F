      program test_shrink_label
      use b2mod_types
      implicit none
      character*20 label

    1 read(*,'(a20)') label
      call shrink_label(label)
      write(*,'(a20)') trim(label)
      if (label.ne.'') goto 1
      stop
      end program test_shrink_label
