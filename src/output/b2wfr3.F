      subroutine b2wfr3 (nprt)
      use b2mod_types
      implicit none
      integer nprt
*     ------------------------------------------------------------------
*     B2WFR3 writes to a formatted file a text containing credits for
*     the AMNS atomic physics rate coefficients.
*     ------------------------------------------------------------------
*$$$  call subini ('b2wfr3')
      write (nprt,101)
      write (nprt,102)
      write (nprt,'()')
      flush (nprt)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
  101 format (/
     &  ' This table of rate coefficients was obtained from the AMNS ',
     &  '(Atomic, Molecular, Nuclear and Surface) database, ',
     &  'maintained by David Coster and Martin O''Mullane. ')
  102 format ()
*     ------------------------------------------------------------------
      end subroutine b2wfr3
