      subroutine b2wfr2 (nprt)
      use b2mod_types
      implicit none
      integer nprt
*     ------------------------------------------------------------------
*     B2WFR2 writes to a formatted file a text containing credits for
*     the ADAS atomic physics rate coefficients.
*     ------------------------------------------------------------------
*$$$  call subini ('b2wfr2')
      write (nprt,101)
      write (nprt,102)
      write (nprt,'()')
      flush (nprt)
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
  101 format (/
     &  ' This table of rate coefficients was obtained from the ADAS',
     &  '(Atomic Data and Analysis Structure) ',
     &  'project, led by H. P. Summers of the University of ',
     &  'Strathclyde. The collisional-radiative model used is ',
     &  'described in /1/. The data quality is characterised by a ',
     &  '"year" number, the higher numbers indicating more complete ',
     &  'physics models. ADAS also supports bundling of charge states ',
     &  'which can be applied to SOLPS runs.'//,
     &  'The database can be consulted at ',
     &  'http://adas.phys.strath.ac.uk/' //,
     &  ' Hydrogen rate coefficients are taken from code by Jon C. We',
     &  'isheit, described in Refs. /2/ and /3/. These data were'/,
     &  ' made available by D. E. Post and D. K. Stotler of Princeton',
     &  ' PPPL. Please contact D. K. Stotler for further information'/,
     &  ' about these rate coefficients. In publications, please refe',
     &  'r to Refs. /2/ and /3/.')
  102 format (/
     &  ' /1/ H. P. Summers, W. J. Dickson, M. G. O''Mullane, N. R. ',
     &  'Badnell, A. D. Whiteford, D. H. Brooks, J. Lang, S. D. Loch ',
     &  'and D. C. Griffin, Ionization state, excited populations and ',
     &  'emission of impurities in dynamic finite density plasmas: I. ',
     &  'The generalized collisional-radiative model for light ',
     &  'elements, Plasma Physics and Controlled Fusion 48 (2006) pp. ',
     &  '263-293.'/,
     &  ' /2/ Jon C. Weisheit, ''Recombination in Dense Plasmas'', J.',
     &  ' Phys. B (Atom. Molec. Phys.) 8 (1975) 2556.'/,
     &  ' /3/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B. ',
     &  'Heifetz and J. C. Weisheit, ''Survey of Atomic Processes in',
     &  ' Edge Plasmas'', J. Nuclear Mater. 121 (1984) 10-16.')
*     ------------------------------------------------------------------
      end subroutine b2wfr2
