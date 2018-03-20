      module b2mod_constants
      use b2mod_types
c     -------------------------------------------------------------------------
c     09/12/2007 xpb : source CODATA 2006 (http://www.nist.gov/)
c     08/19/2011 xpb : source CODATA 2010 (http://www.nist.gov/)
c     03/19/2018 xpb : source CODATA 2014 (http://www.nist.gov/)
c     -------------------------------------------------------------------------
      real (kind=R8) :: pi
      parameter (pi=3.14159 26535 89793 23846 26433 83280 e0_R8)
      real (kind=R8) :: 
     &  c, me, mp, ev, qe, eps0, mu0, avogr, KBolt
      parameter (c=2.99792458e8_R8, me=9.10938356e-31_R8,
     &  mp=1.672621898e-27_R8, ev=1.6021766028e-19_R8, qe=ev,
     &  mu0=4.0e-7_R8*pi, eps0=1.0_R8/(mu0*c*c),
     &  avogr=6.022140857e23_R8, KBolt=1.38064852e-23_R8)
      end module b2mod_constants
