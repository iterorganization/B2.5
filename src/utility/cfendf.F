      logical function cfendf (nget)
      use b2mod_types
      implicit none
      integer nget
*     ------------------------------------------------------------------
*     CFENDF tests whether an un*formatted file is positioned at
*     end-of-file.
*     ------------------------------------------------------------------
      logical lgendf
      character chf*12, idcod*8
      logical streql
      external xerrab, streql
*     ------------------------------------------------------------------
*     Note: It is assumed that the file is positioned either at a
*     header record, or at end-of-file. All error conditions are
*     treated as an end-of-file.
*     ------------------------------------------------------------------
*   ..get form property
      inquire (nget,form=chf,err=91)
*   ..test next record for end-of-file
      lgendf = .true.
      if (streql(chf,'formatted')) then
       read (nget,'(a8)',end=1,err=1) idcod
       lgendf = .false.
      else if (streql(chf,'unformatted')) then
       read (nget,end=1,err=1) idcod
       lgendf = .false.
      else
       call xerrab ('cfendf--invalid file form-property: '//chf)
      endif
    1 continue
*   ..undo the read
      backspace (nget)
*   ..set return value
      cfendf = lgendf
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfendf--error return from inquire statement')
*     ------------------------------------------------------------------
      end function cfendf
