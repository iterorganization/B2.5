*-----------------------------------------------------------------------
*.specification

      subroutine ratamjhyd_mol(nrid, nms, am1, am2, rtt, nt, rtn, nn,
     &  rttr, ntr, rt_mcoll, n_mcoll, 
     &  rt_nd, rt_d, rt_di1, rt_ec_m, rt_ec_E, rt_ic, rt_dr,  
     &  rt_de, rt_di2)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nms, nt, nn, ntr, n_mcoll
      real (kind=R8) ::
     *  am1, am2, rtt(0:nt), rtn(0:nn), rttr(0:ntr),
     *  rt_mcoll(0:n_mcoll)
*   ..input/output arguments
      real (kind=R8) ::
     *  rt_nd(0:nt,0:nn,0:nms),  rt_d(0:nt,0:nn,0:nms),
     *  rt_di1(0:nt,0:nn,0:nms),
     *  rt_ic(0:nt,0:nn,0:nms),  rt_dr(0:nt,0:nn,0:nms),
     *  rt_de(0:nt,0:nn,0:nms),  rt_di2(0:nt,0:nn,0:nms),
     *  rt_ec_m(0:nt,0:ntr,0:n_mcoll,0:nms),
     *  rt_ec_E(0:nt,0:ntr,0:nms)
      integer :: iend
*-----------------------------------------------------------------------
*.documentation
*     TODO: provide warning when outside valid range of EIRENE data
*
*  1. purpose
*
*     RATAMJHYD_MOL uses Amjuel-Hydhel data for the rate coefficients
*     for molecular processes.
*     The call to ratamjhyd_mol must be preceded by a call to one of the
*     comprehensive routines.
*
*  3. description
*
*     rtt represents temperature or temperature per amu; units are eV.
*     rtn represents density; units are m**(-3).
*     rt_nd(i,j,k)  is the rate coefficient for non-dissociative ionization; m**3/s.
*     rt_d(i,j,k)   is the rate coefficient for dissociation; m**3/s.
*     rt_di1(i,j,k) is the rate coefficient for dissociative ionization; m**3/s.
*     rt_ic(i,j,k)  is the rate coefficient for ion conversion; m**3/s.
*     rt_ec_m(i,j,l,k)  is the rate coefficient for elastic collisions; m**3/s. (momentum)
*     rt_ec_E(i,j,k)  is the rate coefficient for elastic collisions; m**3/s. (energy)
*     rt_dr(i,j,k)  is the rate coefficient for dissociative recombination; m**3/s.
*     rt_de(i,j,k)  is the rate coefficient for dissociative excitation; m**3/s.
*     rt_di2(i,j,k) is the rate coefficient for dissociative ionization (from mol. ion); m**3/s.
*
*     All att coordinates rtt(i), rtn(j), and ims(k) (index molecular species). 
*     Furthermore ec_m has dimension mcoll(l) (mach number bulk velocity difference).
*
*  4. references
*
*  5. parameters
*
*     nrid - integer, input.
*     nrid specifies a unit number for formatted input.
*
*     nz - integer, input.
*     nz specifies the nuclear charge of the species for which the rate
*     coefficients are to be evaluated.  In the present code it is
*     required that nz.eq.1.
*
*     am1, am2 - real*8, input.
*     am1 and am2 specify the nuclear mass of the tow atoms making up the molecule 
*     for which the rate coefficients are to be evaluated, in amu. am must be positive.
*
*     rtt - (0:nt) real*8 array, input.
*     rtt(0:nt) specifies the sequence of values of the temperature
*     variable, expressed in eV, for which the rate coefficients are
*     evaluated. rtt(0:nt) must be positive.
*
*     nt - integer, input.
*     nt specifies the upper dimension bound for the array rtt.
*     It will hold that 0.le.nt.
*
*     rtn - (0:nn) real*8 array, input.
*     rtn(0:nn) specifies the sequence of values of the density
*     variable, expressed in m**-3, for which the rate coefficients are
*     evaluated. rtn(0:nn) must be positive.
*
*     nn - integer, input.
*     nn specifies the upper dimension bound for the array rtn.
*     It will hold that 0.le.nn.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical found
      integer n3, i, j
      real (kind=R8) :: rates(0:nt,0:nn), 
     .                  rates_ec_m(0:nt,0:ntr,0:n_mcoll),
     .                  rates_ec_E(0:nt,0:ntr)
      character(420) :: row
      character :: FILNAM*8,H123*4,REAC*9,CRC*3
      real (kind = R8) :: tmin, tmax, nmin, nmax, bremsf, depot
      integer :: nreac, ireac, itok, ier
      parameter (bremsf=1.69e-38_R8/ev, depot=-13.5984340_R8)

*   ..procedures
      logical streql
*   ..initialisation
      data tmin/0.0_R8/, tmax/1E+30_R8/, nmin/0.0_R8/, nmax/1E+30_R8/
      data nreac/0/
      external xertst, b2xvsg, streql, xerrab

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
      call subini ('ratamjhyd_mol')
*   ..Get density and temperature bounds
      call ipgetr('b2amdr_t_min',tmin)
      call ipgetr('b2amdr_t_max',tmax)
      call ipgetr('b2amdr_nn_min',nmin)
      call ipgetr('b2amdr_nn_max',nmax)
*   ..Number of explicitly specified AMJUEL/HYDHEL reactions
      call ipgeti('b2amdr_nreac',nreac)
*   ..set n3
      n3 = (nms+1)*(nt+1)*(nn+1)
*   ..test nt, nn
      call xertst (0.le.nt.and.0.le.nn, 'faulty argument nt, nn')
      call xertst (am1.gt.0.0_R8, 'faulty argument am1')
      call xertst (am2.gt.0.0_R8, 'faulty argument am2')
*   ..test rtt, rtn
      call b2xvsg (nt+1, rtt, 1, 'rtt', '.gt.')
      call b2xvsg (nn+1, rtn, 1, 'rtn', '.gt.')

*   ..Format for reaction card reading
6666  FORMAT (A6,1X,A4,A9,A3)

*   .. Set up default reactions

!! TODO what to do with am2?

!! non-dissociative ionization
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.9    ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_nd(:,:,0) = rates
!! dissociation
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.5g   ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_d(:,:,0) = rates
!! dissociative ionization (from neutral molecule)
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.10    ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_di1(:,:,0) = rates
!! ion conversion 
      call get_amjhyd_rates('AMJUEL  ','H.2 ','3.2.3     ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_ic(:,:,0) = rates

!! elastic collisions
      call get_amjhyd_rates_h3('AMJUEL  ','H.3 ','0.3T      ','EL ',am1, 
     . nt,ntr, n_mcoll,rtt,rttr,rt_mcoll,rates_ec_m, rates_ec_E)
      rt_ec_m(:,:,:,0) = rates_ec_m
      rt_ec_E(:,:,0) = rates_ec_E

!! dissociative recombination
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.14    ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_dr(:,:,0) = rates
!! dissociative excitation
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.12    ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_de(:,:,0) = rates
!! dissociative ionization (from molecular ion)
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.2.11    ','EI ',am1, 
     . nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rt_di2(:,:,0) = rates

* ..return
      call subend ()
      return

   99 call xerrab('No Amjuel/Hydhel reaction header found in b2am.dat!')
      return

*----------------------------------------------------------------------
*.end ratamjhyd_mol

      end subroutine ratamjhyd_mol

