      character*26 function jobnam ()
      use b2mod_types
      implicit none
*     ------------------------------------------------------------------
*     JOBNAM returns an identification of the present job.
*!!   This routine may contain system dependencies.
*     This version combines the date, time and user-name to create a
*     job name.
*     ------------------------------------------------------------------
      character jobnm*26
      save jobnm
      character daytim*17, usrnam*16
      external daytim, usrnam
      data jobnm/' '/
*     ------------------------------------------------------------------
*     (Care is taken to ensure that successive calls to jobnam produce
*     identical results; as jobnam invokes a time routine this is not
*     automatic.)
      if (jobnm.eq.' ') then
*      (this is the first call to jobnam)
       jobnm = daytim()//' '//usrnam()
*      (note that usrnam() will be truncated to eight characters.)
      endif
      jobnam = jobnm
      return
*     ------------------------------------------------------------------
*     Note: When porting this routine to another operating system,
*     use some suitable combination from among day, time, user-id,
*     job-id, and machine-id. Do not lightly change the character*26
*     length specification.
*     ------------------------------------------------------------------
      end
