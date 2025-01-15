*-----------------------------------------------------------------------
*.specification

      subroutine b2cdca ()
      implicit none

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2CDCA provides some general documentation of the B2.5 code.
*     No computation is performed in this routine.
*
*-----------------------------------------------------------------------
*.text
*
*  0. Introduction
*
*     The B2.5 code solves a system of fluid equations that models
*     a multi-species plasma and optionally fluid neutrals in a two-dimensional
*     geometry. Its main application is simulation of the edge plasma in
*     tokamaks.
*
*     B2.5 is composed of several independent programs that communicate
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
*     B2.5 employs a two-dimensional, curvilinear mesh. Through much of its
*     history, the mesh was structured and topologically rectangular. The
*     current version uses an unstructured mesh such that it becomes possible to
*     mesh up to the real vessel walls and to have better control over local
*     refinement. There are nCv (1:nCv) cells, comprised of nCi internal cells
*     and nCg guard cells.
*
*     The primary data describing the mesh are the coordinates of the
*     vertices and cell centers of all the cells of the mesh.
*
*     The primary data describing the magnetic field are the flux
*     functions vxFpsi and vxFfbz at the cell vertices. The magnetic field
*     components at cell centers are derived from vxFpsi and vxFfbz.
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
*     /1/ B. J. Braams, Computational Studies in Tokamak Equilibrium
*     and Transport, thesis, State University at Utrecht, June 1986.
*     /2/ B. J. Braams, A Multi-Fluid Code for Simulation of the Edge
*     Plasma in Tokamaks, NET Report Nr 68, January 1987.
*     /3/ R. A. Hulse, Numerical Studies of Impurities in Fusion
*     Plasmas, Nuclear Technol./Fusion, vol.3 (1983) pp.259-272.
*     /4/ R. A. Hulse, !! (detailed description?).
*     /5/ K. Behringer, Description of the Impurity Transport Code
*     STRAHL, JET Report R(87)08, 1987.
*     /6/ K. Behringer, !! (journal article?).
*     /7/ Jon C. Weisheit, Recombination in Dense Plasmas, J. Phys. B
*     (Atom. Molec. Phys.) 8 (1975) 2556.
*     /8/ R. K. Janev, D. E. Post, W. D. Langer, K. Evans, D. B.
*     Heifetz and J. C. Weisheit, Survey of Atomic Processes in Edge
*     Plasmas, J. Nuclear Mater. 121 (1984) 10-16.
*     /9/ W. Dekeyser, W., et al. "Plasma edge simulations including realistic
*     wall geometry with SOLPS-ITER."
*     Nuclear Materials and Energy 27 (2021): 100999.
*
*-A1. Notes
*
*   ..Originally the code adhered to ANSI Fortran-77.
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
*   ..At present, the code mainly uses Fortran 90, 95 and 2003.
*     Some users HPC systems may not yet support 2003, although this is becoming
*     an exception. The -DLEGACYCOMP preprocessor pragma deals with this.
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
*     The programs that make up B2 communicate through external files,
*     which could be unformatted if all programs run on the same
*     computer system, and otherwise would be formatted. In the latter
*     case the format is very simple, and not designed for convenience
*     of a human reader. We refer to these files as 'un*formatted'.
*     Reading from and writing to an un*formatted file is done through
*     the routines cfruch and cfwuch for character variables, cfruin
*     and cfwuin for integer variables, and cfrure and cfwure for real
*     variables.
*     The un*formatted files could easily be replaced by files written
*     using 'NetCDF' or 'HDF' routines if either one of those packages
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
*   ..Code authors
*     The original author of the code is Dr. Bastiaan J. Braams.
*     Dozens of other authors have contributed since then. See references in
*     the SOLPS-ITER manual.
*
*
*-----------------------------------------------------------------------
*.end b2cdca

      end subroutine b2cdca
