      module b2mod_b2plot_values
      use b2mod_types

      integer nwtf
#ifdef WINDOWS
cxpb  Done to avoid a clash with a reserved file for Windows NCAR
      parameter (nwtf=87)
#else
      parameter (nwtf=97)
#endif
      integer nparameters
      parameter (nparameters=39)
      real (kind=R8) init_value
      parameter (init_value=-66.0_R8)

      real (kind=R8), save :: values(nparameters)

      character*8, save :: parameters(nparameters)

      data parameters /
     1 'xmin', 'xmax', 'ymin', 'ymax', 'rmin',
     2 'rmax', 'zmin', 'zmax', 'fmin', 'fmax',
     3 'omin', 'omax', 'xref', 'wmin', 'wmax',
     4 't_pt', 'rshi', 'zshi', 'r0'  , 'z0'  ,
     5 'a0'  , 'b0'  , 'tt0' , 'trgxmin', 'trgxmax', 
     6 'trgymin', 'trgymax', 'trgrmin', 'trgrmax', 'trgzmin', 
     7 'trgzmax', 'trgfmin', 'trgfmax', 'trgomin', 'trgomax',
     8 'smin', 'smax', 'emin', 'emax'/

      end module b2mod_b2plot_values
