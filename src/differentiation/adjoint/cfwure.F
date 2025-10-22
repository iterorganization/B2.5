      subroutine cfwure (nput, n, refun, id)
      use b2mod_types
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad_diff
     & , only : fmt
      implicit none
      integer nput, n
      real (kind=R8) ::
     *   refun(0:n-1)
      character id*(*)
*     ------------------------------------------------------------------
*     CFWURE writes a real*8 array to an un*formatted file.
*!!   Note: If formatted output is selected, real*8 values must be
*     representable using three digits in the exponent.
*     ------------------------------------------------------------------
      integer ifmt, i
      parameter (ifmt=1)
      character chf*12, idcod*8, idtyp*8, id0*32
      logical streql
      external xerrab, streql
*     ------------------------------------------------------------------
*!!   Note: set ifmt to select the desired format.
*     ------------------------------------------------------------------
*   ..preliminaries
      idcod = '*cf:'
      idtyp = 'real'
      id0 = id
      inquire (nput,form=chf,err=91)
*   ..branch on form
      if (streql(chf,'formatted')) then
*    ..formatted write
       write (nput,'(2a8,i12,4x,a32)') idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput,fmt(ifmt)) (refun(i),i=0,n-1)
       endif
      else if (streql(chf,'unformatted')) then
*    ..unformatted write
       write (nput) idcod, idtyp, n, id0
       if (0.lt.n) then
        write (nput) (refun(i),i=0,n-1)
       endif
      else
       call xerrab ('cfwure--invalid file form-property: '//chf)
      endif
      return
*     ------------------------------------------------------------------
   91 call xerrab ('cfwure--error return from inquire statement')
*     ------------------------------------------------------------------
      end subroutine cfwure
