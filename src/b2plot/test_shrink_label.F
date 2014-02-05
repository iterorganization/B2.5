      program test_shrink_label
      use b2mod_types
      integer lnblnk
      external lnblnk
      character*20 label
    1 read(*,'(a20)') label
      call shrink_label(label)
      write(*,'(a20)') label(1:lnblnk(label))
      if (label.ne.'') goto 1
      end
