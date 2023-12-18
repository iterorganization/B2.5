      subroutine check_string_duplicates(filename)
      ! Description
      !============
      ! Subroutine to check that the strings in b2mn.dat are not
      ! repeated, to avoid unwanted overwriting during reading.
      ! To do so, read the file line by line and check whether
      ! switches are repeated.
      ! If a switch is repeated, an error is returned
      ! Initialize
      !===========
      ! Modules
      use b2mod_subsys

      implicit none
      ! Input variables
      character, intent(in) :: filename*8

      ! Allocatable variables
      character, allocatable :: strings*256(:)
      integer, allocatable :: lastchar(:)

      ! Local variables
      integer :: fstatus, file_length, lstchar_end
      character :: current_line*256, stringtmp*256
      integer :: nlines, ii, jj, charindx
      logical :: repeated

      ! External functions
c      external xerrab, subini, subend

      ! Main computation
      !=============
      ! Call subini
      call subini ('check_string_duplicates')

      open (unit=10, file=filename, status='old', action='read',
     &      iostat=fstatus)
      if (fstatus /= 0) then
         call xerrab ('Error opening file: ' // trim(filename))
      endif

      ! First reading to get file length
      file_length = 0
      do
        read(10, '(A)', end=98) current_line
        if(current_line(1:1)=="'") then
          file_length = file_length + 1
        endif
      enddo

98    continue
      close(10)

      allocate(strings(file_length+1))
      allocate(lastchar(file_length+1))

      nlines = 0
      charindx = 0
      repeated = .FALSE.

      ! Second reading to check switches
      open (unit=11, file=filename, status='old', action='read',
     &      iostat=fstatus)
      if (fstatus /= 0) then
         call xerrab ('Error opening file: ' // trim(filename))
      endif

      do
        read(11, '(A)', end=99) current_line
        if(current_line(1:1)=="'") then
          charindx = charindx + 1
          lastchar(charindx) = -1
          jj = 1
          do while(lastchar(charindx).lt.0)
            jj = jj+1
            if(current_line(jj:jj)=="'") then
              lastchar(charindx) = jj
              lstchar_end = jj
            endif
          enddo
          nlines = nlines + 1
          strings(nlines) = current_line(1:lstchar_end)
          ! Check
          do jj = 1, nlines-1
            if(lastchar(jj).eq.lstchar_end) then
              stringtmp = strings(jj)
              if (trim(current_line(1:lstchar_end)) ==
     &            trim(stringtmp(1:lstchar_end))) then
                repeated = .TRUE.
                close(11)
                call xerrab ('Duplicate switch ' //
     &                      trim(current_line(1:lstchar_end)) //
     &                      ' found in b2mn.dat')
              endif
            endif
          enddo
        endif
      enddo

99    continue

      if(repeated==.FALSE.) then
        close(11)
      endif

      ! Deallocation
      if(allocated(strings)) then
        deallocate(strings)
      endif
      if(allocated(lastchar)) then
        deallocate(lastchar)
      endif

      ! Call subend and return
      call subend ()
      return
      end subroutine check_string_duplicates
