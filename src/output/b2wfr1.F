      subroutine b2wfr1 (nprt)
      use b2mod_types
      implicit none
      integer nprt
*     ------------------------------------------------------------------
*     B2WFR1 writes to a formatted file a text containing credits for
*     the STRAHL atomic physics rate coefficients.
*     ------------------------------------------------------------------
*$$$  call subini ('b2wfr1')
      write (nprt,101)
      write (nprt,102)
      write (nprt,'()')
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
  101 format (/
     &  ' This table of rate coefficients was computed using the STRA',
     &  'HL package, described in Ref. /1/. Please address questions'/,
     &  ' about these rates to the author of STRAHL, K. Behringer of ',
     &  'the University at Stuttgart. K. Behringer will appreciate'/,
     &  ' to be informed about ongoing research that relies on this s',
     &  'et of rate coefficients. In publications, as a minimum,'/,
     &  ' refer to /1/.'//,
     &  ' The CX rate coefficients are provisional, not part of STRAH',
     &  'L, supplied by BJB from /2/; check before use!'/,
     &  ' Hydrogen rate coefficients are taken from code by Jon C. We',
     &  'isheit, described in Refs. /3/ and /4/. These data were'/,
     &  ' made available by D. E. Post and D. K. Stotler of Princeton',
     &  ' PPPL. Please contact D. K. Stotler for further information'/,
     &  ' about these rate coefficients. In publications, please refe',
     &  'r to Refs. /3/ and /4/.')
  102 format (/
     &  ' /1/ K. Behringer, ''Description of the Impurity Transport C',
     &  'ode STRAHL'', JET Report R(87)08, 1987.'/,
     &  ' /2/ M. E. Puiatti et al., Plasma Physics, vol. 23 (1981) pp',
     &  '. 1075-1092.'/,
     &  ' /3/ Jon C. Weisheit, ''Recombination in Dense Plasmas'', J.',
     &  ' Phys. B (Atom. Molec. Phys.) 8 (1975) 2556.'/,
     &  ' /4/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B. ',
     &  'Heifetz and J. C. Weisheit, ''Survey of Atomic Processes in',
     &  ' Edge Plasmas'', J. Nuclear Mater. 121 (1984) 10-16.')
*     ------------------------------------------------------------------
      end subroutine b2wfr1
