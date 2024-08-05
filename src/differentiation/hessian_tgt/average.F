
*-----------------------------------------------------------------------
*.specification

      subroutine average (e,b,w0,phi,mu)
      use b2mod_types
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      real (kind=R8) ::
     *  e, b, w0, phi
*   ..output argument (unspecified on entry)
      real (kind=R8) ::
     *  mu

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     This subroutine calculates the average value over the emission angle
*     of the ratio mu = g/g0, where g0 is the secondary electron emission
*     coefficient for a given surface (characteristic of the material) and
*     g is the effective secondary electron emission coefficient, dependent
*     on the emission angle.
*
*     A description of the used parameters follows. Units are not specified.
*     It is assumed that MKS units are used everywhere.
*     Input parameters:
*     e     - Modulus of the electric field.
*     b     - Modulus of the magnetic field.
*     w0    - Work function for the material constituting the given surface.
*     phi   - Incidence angle of the magnetic field with respect to the solid
*             surface (in degrees).
*
*     Output parameter:
*     mu    - <g/g0>, averaged value.
*
*     Internal parameters
*     nv    - Number of points in velocity space - 1.
*     np    - Number of points in angle space - 1.
*     gg0   - gg0(ip,iv) is the ratio g/g0 for the (ip,iv) point of the
*             angle-velocity space.
*     v     - Normalised velocity array.
*     v0    - Physical velocity array. It will hold that v0=v*e/b.
*     ph    - Incidence angle array (in degrees).
*     acc   - Required accuracy.
*     vwrk0 - Initial velocity work array.
*     gwrk0 - Initial gg0 work array.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nv, np, ite
      parameter (nv=11, np=27)
      real (kind=R8) :: gg0(0:np,0:nv)
      real (kind=R8) :: v(0:nv),vwrk0(0:nv),gwrk0(0:nv),ph(0:np)
      real (kind=R8) :: acc
      parameter (acc=1.e-2_R8)
*   ..procedures
      external xertst, ggfill, setwrk0, integrate
*   ..initialisation
*
*     The following table assign the values of v, ph and gg0. Note some
*     values of gg0 are set to -1. They correspond to missing data. At
*     the start of the computation, the missing values will be calculated
*     by linear interpolation.
*
      data           v(0),      v(1),      v(2),      v(3),      v(4),
     &               v(5),      v(6),      v(7),      v(8),      v(9),
     &              v(10),     v(11)
     &
     &          /0.0e0_R8, 5.0e-2_R8, 1.0e-1_R8, 2.0e-1_R8, 3.0e-1_R8,
     &          6.0e-1_R8,  1.0e0_R8,  2.0e0_R8,  4.0e0_R8,  8.0e0_R8,
     &           2.0e1_R8, 1.0e12_R8/

      data          ph(0),     ph(1),     ph(2),     ph(3),     ph(4),
     &              ph(5),     ph(6),     ph(7),     ph(8),     ph(9),
     &             ph(10),    ph(11),    ph(12),    ph(13),    ph(14),
     &             ph(15),    ph(16),    ph(17),    ph(18),    ph(19),
     &             ph(20),    ph(21),    ph(22),    ph(23),    ph(24),
     &             ph(25),    ph(26),    ph(27)
     &
     &          /1.0e0_R8,  1.5e0_R8,  2.0e0_R8,  3.0e0_R8,  4.0e0_R8,
     &           5.0e0_R8,  5.4e0_R8,  6.0e0_R8,  7.0e0_R8,  8.0e0_R8,
     &           1.0e1_R8,  1.1e1_R8,  1.4e1_R8,  1.8e1_R8,  1.9e1_R8,
     &           2.0e1_R8,  2.5e1_R8,  3.0e1_R8,  3.5e1_R8,  4.0e1_R8,
     &           4.5e1_R8,  5.0e1_R8,  6.0e1_R8,  6.5e1_R8,  7.0e1_R8,
     &           7.5e1_R8,  8.0e1_R8,  8.5e1_R8/

      data  gg0( 0, 0), gg0( 0, 1), gg0( 0, 2), gg0( 0, 3), gg0( 0, 4),
     &      gg0( 0, 5), gg0( 0, 6), gg0( 0, 7), gg0( 0, 8), gg0( 0, 9),
     &      gg0( 0,10), gg0( 0,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,  6.6e-1_R8,  3.4e-1_R8,  2.3e-1_R8,
     &       1.1e-1_R8,  7.0e-2_R8,  4.0e-2_R8,  3.0e-2_R8,  2.0e-2_R8,
     &       2.0e-2_R8,  2.0e-2_R8/

      data  gg0( 1, 0), gg0( 1, 1), gg0( 1, 2), gg0( 1, 3), gg0( 1, 4),
     &      gg0( 1, 5), gg0( 1, 6), gg0( 1, 7), gg0( 1, 8), gg0( 1, 9),
     &      gg0( 1,10), gg0( 1,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,  9.1e-1_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  4.0e-2_R8,
     &       4.0e-2_R8,  3.0e-2_R8/

      data  gg0( 2, 0), gg0( 2, 1), gg0( 2, 2), gg0( 2, 3), gg0( 2, 4),
     &      gg0( 2, 5), gg0( 2, 6), gg0( 2, 7), gg0( 2, 8), gg0( 2, 9),
     &      gg0( 2,10), gg0( 2,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  6.6e-1_R8,  4.5e-1_R8,
     &       -1.0e0_R8,  1.4e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  5.0e-2_R8/

      data  gg0( 3, 0), gg0( 3, 1), gg0( 3, 2), gg0( 3, 3), gg0( 3, 4),
     &      gg0( 3, 5), gg0( 3, 6), gg0( 3, 7), gg0( 3, 8), gg0( 3, 9),
     &      gg0( 3,10), gg0( 3,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  9.1e-1_R8,  6.6e-1_R8,
     &       3.5e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  7.0e-2_R8/

      data  gg0( 4, 0), gg0( 4, 1), gg0( 4, 2), gg0( 4, 3), gg0( 4, 4),
     &      gg0( 4, 5), gg0( 4, 6), gg0( 4, 7), gg0( 4, 8), gg0( 4, 9),
     &      gg0( 4,10), gg0( 4,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  8.4e-1_R8,
     &       -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.00_R8/

      data  gg0( 5, 0), gg0( 5, 1), gg0( 5, 2), gg0( 5, 3), gg0( 5, 4),
     &      gg0( 5, 5), gg0( 5, 6), gg0( 5, 7), gg0( 5, 8), gg0( 5, 9),
     &      gg0( 5,10), gg0( 5,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  9.7e-1_R8,
     &       5.5e-1_R8,  3.4e-1_R8,  2.0e-1_R8,  1.4e-1_R8,  1.1e-1_R8,
     &       1.0e-1_R8,  9.0e-2_R8/

      data  gg0( 6, 0), gg0( 6, 1), gg0( 6, 2), gg0( 6, 3), gg0( 6, 4),
     &      gg0( 6, 5), gg0( 6, 6), gg0( 6, 7), gg0( 6, 8), gg0( 6, 9),
     &      gg0( 6,10), gg0( 6,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0( 7, 0), gg0( 7, 1), gg0( 7, 2), gg0( 7, 3), gg0( 7, 4),
     &      gg0( 7, 5), gg0( 7, 6), gg0( 7, 7), gg0( 7, 8), gg0( 7, 9),
     &      gg0( 7,10), gg0( 7,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0( 8, 0), gg0( 8, 1), gg0( 8, 2), gg0( 8, 3), gg0( 8, 4),
     &      gg0( 8, 5), gg0( 8, 6), gg0( 8, 7), gg0( 8, 8), gg0( 8, 9),
     &      gg0( 8,10), gg0( 8,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       7.5e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  1.2e-1_R8/

      data  gg0( 9, 0), gg0( 9, 1), gg0( 9, 2), gg0( 9, 3), gg0( 9, 4),
     &      gg0( 9, 5), gg0( 9, 6), gg0( 9, 7), gg0( 9, 8), gg0( 9, 9),
     &      gg0( 9,10), gg0( 9,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       8.3e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(10, 0), gg0(10, 1), gg0(10, 2), gg0(10, 3), gg0(10, 4),
     &      gg0(10, 5), gg0(10, 6), gg0(10, 7), gg0(10, 8), gg0(10, 9),
     &      gg0(10,10), gg0(10,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       9.6e-1_R8,  6.5e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  2.2e-1_R8,
     &       -1.0e0_R8,  1.7e-1_R8/

      data  gg0(11, 0), gg0(11, 1), gg0(11, 2), gg0(11, 3), gg0(11, 4),
     &      gg0(11, 5), gg0(11, 6), gg0(11, 7), gg0(11, 8), gg0(11, 9),
     &      gg0(11,10), gg0(11,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(12, 0), gg0(12, 1), gg0(12, 2), gg0(12, 3), gg0(12, 4),
     &      gg0(12, 5), gg0(12, 6), gg0(12, 7), gg0(12, 8), gg0(12, 9),
     &      gg0(12,10), gg0(12,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,  8.4e-1_R8,  5.2e-1_R8,  3.7e-1_R8,   3.e-1_R8,
     &       2.7e-1_R8,  2.4e-1_R8/

      data  gg0(13, 0), gg0(13, 1), gg0(13, 2), gg0(13, 3), gg0(13, 4),
     &      gg0(13, 5), gg0(13, 6), gg0(13, 7), gg0(13, 8), gg0(13, 9),
     &      gg0(13,10), gg0(13,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,  9.8e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(14, 0), gg0(14, 1), gg0(14, 2), gg0(14, 3), gg0(14, 4),
     &      gg0(14, 5), gg0(14, 6), gg0(14, 7), gg0(14, 8), gg0(14, 9),
     &      gg0(14,10), gg0(14,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,  9.9e-1_R8,  -1.0e0_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(15, 0), gg0(15, 1), gg0(15, 2), gg0(15, 3), gg0(15, 4),
     &      gg0(15, 5), gg0(15, 6), gg0(15, 7), gg0(15, 8), gg0(15, 9),
     &      gg0(15,10), gg0(15,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  6.9e-1_R8,  5.1e-1_R8,  4.3e-1_R8,
     &       3.7e-1_R8,  3.4e-1_R8/

      data  gg0(16, 0), gg0(16, 1), gg0(16, 2), gg0(16, 3), gg0(16, 4),
     &      gg0(16, 5), gg0(16, 6), gg0(16, 7), gg0(16, 8), gg0(16, 9),
     &      gg0(16,10), gg0(16,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  8.1e-1_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  4.2e-1_R8/

      data  gg0(17, 0), gg0(17, 1), gg0(17, 2), gg0(17, 3), gg0(17, 4),
     &      gg0(17, 5), gg0(17, 6), gg0(17, 7), gg0(17, 8), gg0(17, 9),
     &      gg0(17,10), gg0(17,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  8.9e-1_R8,  7.1e-1_R8,  6.1e-1_R8,
     &       5.4e-1_R8,  5.0e-1_R8/

      data  gg0(18, 0), gg0(18, 1), gg0(18, 2), gg0(18, 3), gg0(18, 4),
     &      gg0(18, 5), gg0(18, 6), gg0(18, 7), gg0(18, 8), gg0(18, 9),
     &      gg0(18,10), gg0(18,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  9.5e-1_R8,  -1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(19, 0), gg0(19, 1), gg0(19, 2), gg0(19, 3), gg0(19, 4),
     &      gg0(19, 5), gg0(19, 6), gg0(19, 7), gg0(19, 8), gg0(19, 9),
     &      gg0(19,10), gg0(19,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  9.8e-1_R8,  8.6e-1_R8,  7.6e-1_R8,
     &       6.9e-1_R8,  6.4e-1_R8/

      data  gg0(20, 0), gg0(20, 1), gg0(20, 2), gg0(20, 3), gg0(20, 4),
     &      gg0(20, 5), gg0(20, 6), gg0(20, 7), gg0(20, 8), gg0(20, 9),
     &      gg0(20,10), gg0(20,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,  9.9e-1_R8,  -1.0e0_R8,  8.2e-1_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(21, 0), gg0(21, 1), gg0(21, 2), gg0(21, 3), gg0(21, 4),
     &      gg0(21, 5), gg0(21, 6), gg0(21, 7), gg0(21, 8), gg0(21, 9),
     &      gg0(21,10), gg0(21,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  9.4e-1_R8,  8.8e-1_R8,
     &       8.1e-1_R8,  7.7e-1_R8/

      data  gg0(22, 0), gg0(22, 1), gg0(22, 2), gg0(22, 3), gg0(22, 4),
     &      gg0(22, 5), gg0(22, 6), gg0(22, 7), gg0(22, 8), gg0(22, 9),
     &      gg0(22,10), gg0(22,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  9.9e-1_R8,  9.5e-1_R8,
     &       9.1e-1_R8,  8.7e-1_R8/

      data  gg0(23, 0), gg0(23, 1), gg0(23, 2), gg0(23, 3), gg0(23, 4),
     &      gg0(23, 5), gg0(23, 6), gg0(23, 7), gg0(23, 8), gg0(23, 9),
     &      gg0(23,10), gg0(23,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  -1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(24, 0), gg0(24, 1), gg0(24, 2), gg0(24, 3), gg0(24, 4),
     &      gg0(24, 5), gg0(24, 6), gg0(24, 7), gg0(24, 8), gg0(24, 9),
     &      gg0(24,10), gg0(24,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  9.9e-1_R8,
     &       9.7e-1_R8,  9.4e-1_R8/

      data  gg0(25, 0), gg0(25, 1), gg0(25, 2), gg0(25, 3), gg0(25, 4),
     &      gg0(25, 5), gg0(25, 6), gg0(25, 7), gg0(25, 8), gg0(25, 9),
     &      gg0(25,10), gg0(25,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &       -1.0e0_R8,  -1.0e0_R8/

      data  gg0(26, 0), gg0(26, 1), gg0(26, 2), gg0(26, 3), gg0(26, 4),
     &      gg0(26, 5), gg0(26, 6), gg0(26, 7), gg0(26, 8), gg0(26, 9),
     &      gg0(26,10), gg0(26,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,  -1.0e0_R8,
     &        1.0e0_R8,  9.9e-1_R8/
     &
      data  gg0(27, 0), gg0(27, 1), gg0(27, 2), gg0(27, 3), gg0(27, 4),
     &      gg0(27, 5), gg0(27, 6), gg0(27, 7), gg0(27, 8), gg0(27, 9),
     &      gg0(27,10), gg0(27,11)
     &
     &       /1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,   1.0e0_R8,
     &        1.0e0_R8,   1.0e0_R8/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('average')
*   ..test input parameters
      call xertst (0.lt.b, 'faulty input b')
      call xertst (0.lt.phi, 'faulty input phi')

*   ..computation
      call ggfill (ph,gg0,np,nv)                               !Tested ok
      call setwrk0 (phi,ph,v,gg0,np,nv,e,b,vwrk0,gwrk0)        !Tested ok
      call integrate (vwrk0,gwrk0,w0,acc,nv,mu,ite)            !Tested ok

*   ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end average

      end subroutine average

