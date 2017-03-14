*-----------------------------------------------------------------------
*.specification

      subroutine b2cdca ()
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CDCA provides some general documentation of the b2 code.
*     No computation is performed in this routine.
*
*-----------------------------------------------------------------------
*.text
*
*  0. Introduction
*
*     The b2 code /1,2/ solves a system of fluid equations that models
*     a multi-species plasma in a two-dimensional geometry. Its main
*     application is simulation of the edge plasma in tokamaks.
*
*     B2 is composed of several independent programs that communicate
*     through external files. The main sequence involves the programs
*     b2ag (prepares a geometry and magnetic field), b2ah (prepares
*     a table of default physics parameters), b2ai (prepares an
*     initial plasma state), b2mn (principal calculation), and b2yp
*     (postprocessing).
*     Additional programs may prepare tables of rate coefficients for
*     atomic processes, tables of plasma-wall interaction processes,
*     external source terms, and other physics, or may perform further
*     postprocessing.
*
*     ...
*
*  1. Mathematical model
*
*     ...
*
*  2. Geometry
*
*     B2 employs a two-dimensional, curvilinear, topologically
*     rectangular mesh. The coordinates of the mesh are labelled 'x'
*     and 'y'. There are nx+2 cells along the x-coordinate and ny+2
*     cells along the y-coordinate. The cells of the mesh are indexed
*     by (ix,iy) in (-1:nx,-1:ny); the interior cells have index in
*     (0:nx-1,0:ny-1). The mesh may contain cuts, which typically
*     originate from an x-point in the real*8 geometry.
*
*     The primary data describing the mesh are the coordinates of the
*     four vertices of all the cells of the mesh. These data are stored
*     in the arrays crx and cry. The metric parameters (vol, hx, hy, qz,
*     qc and gs; described elsewhere) are derived from crx and cry.
*
*     The primary data describing the magnetic field are the flux
*     functions fpsi and ffbz at the cell vertices. The magnetic field
*     components at cell centers are derived from fpsi and ffbz.
*
*     Singular regions in the geometry are specified by the data in the
*     commmon block /b2cmgs/.
*
*     ...
*
*  3. Numerical method
*
*     ...
*
*  4. References
*
*     (Refs. /1/ and /2/ are out of date; !! to be replaced.)
*     /1/ B. J. Braams, "Computational Studies in Tokamak Equilibrium
*     and Transport", thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, "A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks", NET Report Nr 68, January 1987.
*     /3/ R. A. Hulse, "Numerical Studies of Impurities in Fusion
*     Plasmas", Nuclear Technol./Fusion, vol.3 (1983) pp.259-272.
*     /4/ R. A. Hulse, !! (detailed description?).
*     /5/ K. Behringer, "Description of the Impurity Transport Code
*     STRAHL", JET Report R(87)08, 1987.
*     /6/ K. Behringer, !! (journal article?).
*     /7/ Jon C. Weisheit, "Recombination in Dense Plasmas", J. Phys. B
*     (Atom. Molec. Phys.) 8 (1975) 2556.
*     /8/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B.
*     Heifetz and J. C. Weisheit, "Survey of Atomic Processes in Edge
*     Plasmas", J. Nuclear Mater. 121 (1984) 10-16.
*
*-A1. Notes
*
*   ..Adherence to ANSI Fortran-77.
*     We respect the standard with the following exceptions.
*     - 'implicit none' is used to force type declarations for all
*       constants, variables and functions. We also declare, as external
*       or intrinsic, all procedures.
*     - 'do, enddo' syntax is used; a pre-processor based on Emacs
*       (courtesy of Charles Karney) may convert to labelled do loops.
*     - CPP 'include' statements are used for common blocks and for some
*       parameter statements.
*     - symbolic names may have up to eight characters.
*     These extensions to ANSI Fortran-77 are accepted in the identical
*     syntax on all popular workstations as well as on the Cray.
*     The code can be distributed using do-enddo syntax or using
*     labelled do loops, and with or without use of include statements.
*
*   ..Scientific libraries.
*     The BLAS (basic linear algebra subroutines) and LAPACK (linear
*     algebra package) are employed; these are included in the standard
*     scientific libraries on many computer systems, and are also in
*     NAG and SLATEC; moreover, they are available by electronic mail
*     or ftp from netlib@ornl.gov.  Some simple additional routines
*     from the NAG library may be employed, for constants of the machine
*     arithmetic and for time and date.
*     An implementation of these scientific libraries is required in
*     which real*8 variables occupy eight bytes.
*
*   ..Graphics libraries.
*     As distributed, B2 relies on NCAR and GKS for graphics output.
*     The graphics postprocessing codes for B2 are completely separate
*     from the main calculation.
*
*   ..Concerning 'un*formatted' files.
*     The programs that make up b2 communicate through external files,
*     which could be unformatted if all programs run on the same
*     computer system, and otherwise would be formatted. In the latter
*     case the format is very simple, and not designed for convenience
*     of a human reader. We refer to these files as 'un*formatted'.
*     Reading from and writing to an un*formatted file is done through
*     the routines cfruch and cfwuch for character variables, cfruin
*     and cfwuin for integer variables, and cfrure and cfwure for real
*     variables.
*     The un*formatted files could easily be replaced by files written
*     using 'netCDF' or 'HDF' routines if either one of those packages
*     is desired as a standard interface.
*
*   ..Concerning internal documentation.
*     The style of the internal documentation follows that of the NAG
*     library routine documents.
*     It has been the intent to specify the meaning of any non-trivial
*     variable in the outermost routine in which it is declared, and
*     also in the routine in which it is first computed.
*
*   ..Special patterns in comments.
*     The pattern !! is used to draw attention to segments of code
*     that might have to be changed. The pattern *$$$ in columns 1-4 is
*     used for code that is commented out with the text retained as a
*     reminder.
*
*   ..Contacting the author.
*     Use e-mail to braams@cims.nyu.edu, phone to 212-998-3276,
*     FAX to 212-995-4121, or regular mail to
*       B. J. Braams
*       Courant Institute of Mathematical Sciences
*       New York University
*       251 Mercer Street
*       New York, NY 10012
*
*-----------------------------------------------------------------------
*.end b2cdca

      end subroutine b2cdca
