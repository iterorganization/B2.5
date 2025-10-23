      subroutine cfwuin (nput, n, infun, id)
      use b2mod_types
      use b2mod_ad_diffv_diffv
     , , only: fmtwuin
      implicit none
      integer nput, n, infun(0:n-1)
      character id*(*)
*     ------------------------------------------------------------------
*     CFWUIN writes an integer array to an un*formatted file.
*     ------------------------------------------------------------------
      integer ifmt, i
      parameter (ifmt=2)
      character chf*12, idcod*8, idtyp*8, id0*32
      logical streql
      external xerrab, streql
*     ------------------------------------------------------------------
*!!   Note: set ifmt to select the desired format.
*     ------------------------------------------------------------------
*   ..preliminaries
      idcod = '*cf:'
      idtyp = 'int'
      id0 = id
      inquire (nput,form=chf,err=91)
*   ..branch on form
      if (streql(chf,'formatted')) then
*    ..formatted write
       write (nput,'(2a8,i12,4x,a32)') idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput,fmtwuin(ifmt)) (infun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted write
       write (nput) idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput) (infun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfwuin--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfwuin--error return from inquire statement')
*     ------------------------------------------------------------------
      end subroutine cfwuin
