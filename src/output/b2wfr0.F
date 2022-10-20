      subroutine b2wfr0 (nprt)
      use b2mod_types
      implicit none
      integer nprt
*     ------------------------------------------------------------------
*     B2WFR0 writes to a formatted file a text containing credits for
*     the ADPAK atomic physics rate coefficients.
*     ------------------------------------------------------------------
*$$$  call subini ('b2wfr0')
      write (nprt,101)
      write (nprt,102)
      write (nprt,'()')
*$$$  call subend ()
      return
*     ------------------------------------------------------------------
  101 format (/
     &  ' This table of rate coefficients was obtained using ADPAK: a',
     &  ' collection of routines developed by R. A. Hulse of the'/,
     &  ' Princeton Plasma Physics Laboratory. ADPAK is described in ',
     &  'Ref. /1/; it builds upon work by Post et al. /2/, Hulse'/,
     &  ' et al. /3/, and Bell et al. /4/. Russel Hulse will apprecia',
     &  'te to be informed about ongoing research that relies on'/,
     &  ' this set of rate coefficients. Contact him also in case of ',
     &  'questions about these rates. In publications, as a minimum,'/,
     &  ' refer to /1/.'//,
     &  ' Hydrogen rate coefficients are taken from code by Jon C. We',
     &  'isheit, described in Refs. /5/ and /6/. These data were'/,
     &  ' made available by D. E. Post and D. K. Stotler of Princeton',
     &  ' PPPL. Please contact D. K. Stotler for further information'/,
     &  ' about these rate coefficients. In publications, please refe',
     &  'r to Refs. /5/ and /6/.')
  102 format (/
     &  ' /1/ R. A. Hulse, Numerical Studies of Impurities in Fusion ',
     &  'Plasmas, Nucl. Technol./Fusion, vol.3 (1983) pp.259-272.'/,
     &  ' /2/ D. E. Post, R. V. Jensen, C. B. Tarter, W. H. Grasberge',
     &  'r, and W. A. Lokke, Steady-State Radiative Cooling Rates'/,
     &  ' for Low-Density, High-Temperature Plasma, At. Data Nucl. Da',
     &  'ta Tables vol.20 (1977) pp.397-.'/,
     &  ' /3/ R. A. Hulse, D. E. Post, and D. R. Mikkelsen, Charge Ex',
     &  'change as a Recombination Mechanism in High-Temperature'/,
     &  ' Plasmas, J. Phys. B, vol.13 (1980) pp.3895-.'/,
     &  ' /4/ K. L. Bell, H. B. Gilbody, J. G. Hughes, A. E. Kingston',
     &  ', and F. J. Smith, Recommended Data on the Electron Impact'/,
     &  ' Ionization of Light Atoms and Ions, J. Phys. Chem. Ref. Dat',
     &  'a, vol. 12 (1983) pp. 891-916.'/,
     &  ' /5/ Jon C. Weisheit, ''Recombination in Dense Plasmas'', J.',
     &  ' Phys. B (Atom. Molec. Phys.) 8 (1975) 2556.'/,
     &  ' /6/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B. ',
     &  'Heifetz and J. C. Weisheit, ''Survey of Atomic Processes in'/,
     &  ' Edge Plasmas'', J. Nuclear Mater. 121 (1984) 10-16.')
*     ------------------------------------------------------------------
      end subroutine b2wfr0
