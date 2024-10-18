      module b2mod_create_afn_bc_tables
      use b2mod_types
      use b2mod_constants
      use b2mod_math
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
c**************************************************************************************
c*********   Module to create new databases for advanced fluid neutral BCs  ***********
c**************************************************************************************

c     This module creates the pre-evaluated integrals of TRIM data to consistently calculate
c     recycling and reflection coefficients based on the physics models used by EIRENE.
c     The routines in this module are called from the pre-processor routine b2abdr.F
c     This module is based on the theory presented in Chapter 4 of the PhD of N. Horsten [1],
c     and the original Nuclear Fusion paper [2].

c     References:
c     1.
c     Fluid and Hybrid Fluid-Kinetic Models for the Neutral Particles in the Plasma Edge of Nuclear Fusion Devices
c     PhD Thesis N. Horsten, KU Leuven, 2019
c     2.
c     Horsten, Niels, Giovanni Samaey, and M. Baelmans.
c     "Development and assessment of 2D fluid neutral models that include atomic databases and a microscopic reflection model."
c     Nuclear Fusion 57.11 (2017): 116043.

cc global variables
      real (kind=R8), allocatable, save ::
     &  M(:),T(:), Sh(:)
      real (kind=R8), save :: dM
      real (kind=R8), allocatable, save ::
     &   R1(:,:,:), R2(:,:,:), Rmx1(:,:,:), Rmx2(:,:,:),
     &  Rmy1(:,:,:), Rmy2(:,:,:), Rmz1(:,:,:), Rmz2(:,:,:), RE1(:,:,:),
     &  RE2(:,:,:), Rp11(:,:,:), Rp12(:,:,:), Rp21(:,:,:), Rp22(:,:,:),
     &  Rp31(:,:,:), Rp32(:,:,:)
      real (kind=R8), save :: Tmin, Tmax, Mmin, Mmax, Shmin, Shmax,
     &                        Eomin, mass
      real (kind=R8), save :: TRIMdata(2688,7)
      integer, save :: nv, ntheta, ntp, nT, nM, nSh, nphi
      integer, save :: nE = 12, nA = 7, cE = 5, cA = 6, cR = 7, !! these values are hard-coded because the format of the TRIM data is not expected to change
     &                 nB = 32
      real (kind=R8), allocatable, save :: Rv(:,:),
     &  Rmomparv(:,:),Rmomperpv(:,:),Rpv1(:,:),Rpv2(:,:),
     &  Rpv3(:,:),REv(:,:)
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif

cc description of global variables
c    M: discretization vector toroidal Mach number
c    T: discretization vector ion temperature (eV)
c    Sh:discretization vector sheath potential factor
c    dM: distance between discretization points of M

c    R1: macroscopic particle reflection coefficient scaling with density
c    R2: macroscopic particle reflection coefficient scaling with normal gradient of density
c    Rmx1: macroscopic toroidal momentum coefficient scaling with density
c    Rmx2: macroscopic toroidal momentum coefficient scaling with normal gradient of density
c    Rmy1: zero
c    Rmy2: macroscopic tangential momentum coefficient scaling with tangential gradient of density
c    Rmz1: macroscopic normal momentum coefficient scaling with density
c    Rmz2: macroscopic normal momentum coefficient scaling with normal gradient of density
c    RE1: macroscopic energy reflection coefficient scaling with density
c    RE2: macroscopic energy reflection coefficient scaling with normal gradient of density
c    Rp11: coefficient for macroscopic contribution in pressure of reflected particles, see PhD N. Horsten Eq. 4.62, part that scales with density
c    Rp12: coefficient for macroscopic contribution in pressure of reflected particles, see PhD N. Horsten Eq. 4.62, part that scales with normal gradient of density
c    Rp21: similar as Rp11, but for the contribution from Eq 4.67
c    Rp22: similar as Rp12, but for the contribution from Eq 4.67
c    Rp31: similar as Rp11, but for the contribution from Eq 4.70
c    Rp32: similar as Rp11, but for the contribution from Eq 4.70

c    Tmin: minimum ion temperature in discretization
c    Tmax: maximum ion temperature in discretization
c    Mmin: minimum toroidal Mach number in discretization
c    Mmax: maximum toroidal Mach number in discretization
c    Shmin: minimum sheath potential factor in discretization
c    Shmax: maximum sheath potential factor in discretization
c    Eomin: minimum energy of reflected particles (eV)
c    mass: mass of the incident particle (kg)

c    TRIMdata: local copy of a TRIM data file from EIRENE

c    nv: number of discretization points in incident velocity v
c    ntheta: number of discretization points in incident polar angle
c    ntp: number of Monte Carlo particles used in Phase 1
c    nT: number of elements in T
c    nM: number of elements in M
c    nSh: number of elements in Sh
c    nphi: number of discretization points in azimuthal angle
c

c    nE Number of incident energy levels (as defined in TRIM)
c    nA Number of incident angles (as defined in TRIM)
c    cE, cA: Column-index of energy and angle of incidence in TRIM-matrix
c    cR: Column-index of reflection coefficient
c    nB: Number of data rows in one block for given angle and energy of incidence (as defined in TRIM)

c    Rv: fast particle reflection coefficient, see PhD N. Horsten Eq. 4.13
c    Rmomparv: fast reflected parallel momentum factor, see PhD N. Horsten Eq. 4.55
c    Rmomperpv: fast reflected perpendicular momentum factor, see PhD N. Horsten Eq. 4.56
c    Rpv1: contribution in pressure of reflected particles, see PhD N. Horsten Eq. 4.62
c    Rpv2: contribution in pressure of reflected particles, see PhD N. Horsten Eq. 4.67
c    Rpv3: contribution in pressure of reflected particles, see PhD N. Horsten Eq. 4.70
c    REv: energy reflection coefficient, see PHD N. Horsten Eq. 4.73

      contains

c**************************************************************************************
c************************   Phase 0: setting up inputs      ***************************
c**************************************************************************************

      subroutine SetInput(nT_in, nM_in, nSh_in, ntp_in, nv_in,
     &                  ntheta_in, nphi_in, Tmin_in, Tmax_in, Mmin_in,
     &                  Mmax_in, Shmin_in, Shmax_in, Eomin_in)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nT_in, nM_in, nSh_in,
     &                       ntp_in, nv_in, ntheta_in, nphi_in
      real (kind=R8), intent(in) :: Tmin_in, Tmax_in, Mmin_in, Mmax_in,
     &                              Shmin_in, Shmax_in, Eomin_in
      call subini ('SetInput')

cc set the global variables of this module to the values read by b2abdr.F
      ntp = ntp_in
      nv = nv_in
      ntheta = ntheta_in
      nphi = nphi_in

      nT = nT_in
      nM = nM_in
      nSh = nSh_in

      Tmin = Tmin_in
      Tmax = Tmax_in
      Mmin = Mmin_in
      Mmax = Mmax_in
      Shmin = Shmin_in
      Shmax = Shmax_in
      Eomin = Eomin_in

      call subend ()
      return
      end subroutine SetInput


      subroutine SetMass(incident)
      use b2mod_types
      use b2mod_subsys
      use b2mod_xerset
      implicit none
      character(len=*), intent(in) :: incident
      logical streql
      external streql, xerrab

      call subini ('SetMass')
cc set the correct mass depending on the incident particle
      if (streql(incident,'H')) then
        mass = 1*mp
      elseif (streql(incident,'D')) then
        mass = 2*mp
      elseif (streql(incident,'T')) then
        mass = 3*mp
      else
        call xerrab('Only hydrogen isotopes supported for now')
      endif

      call subend ()
      return
      end subroutine SetMass

      subroutine SetV(v)
      use b2mod_types
      use b2mod_subsys
      implicit none
** output
      real (kind=R8), intent(out) :: v(nv)
** local
      real (kind=R8) :: logv(nv), mD, x
      real (kind=R8) :: Tminl, Tmaxl
      integer :: i

cc  description of variables

c   v:  the discretization vector of the incident particle velocity (m/s)
c   logv: natural logarithm of v
c   mD: mass of deuterium
c   x: dummy variable
c   Tminl: local value for Tmin used for discretizing v (J)
c   Tmaxl: local value for Tmax used for discretizing v (J)
c   i: loop variable

      call subini ('SetV')

cc set the discretization vector v according to eq. 4.75 of PhD N. Horsten (KU Leuven, 2019)

cc set local values voor backward compatibility
      Tminl = Tmin*eV
      Tmaxl = Tmax*eV
cc the discretization of v is for now always calculated using the mass of deuterium in this equation, such that the discretization space does not depend on the incident particle.
      mD = 2*mp
cc see PhD N. Horsten eq. 4.75
      do i = 1,nv
        x = (real(i)-1)/(nv-1)
        logv(i) = log(0.1*sqrt(Tminl/mD)) + x*
     &        (log(10*sqrt(Tmaxl/mD))-log(0.1*sqrt(Tminl/mD)))
        v(i) = exp(logv(i))
      enddo

      call subend ()
      return
      end subroutine SetV


c**************************************************************************************
c*********   Phase 1: calculate average response of an incident particle    ***********
c**************************************************************************************

      subroutine CalcCoefficientsTRIM()
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
#ifdef LEGACYCOMP
      use b2mod_math
#endif
      implicit none
** local variables
      real (kind=R8) :: v(nv), dtheta, theta(ntheta), Ei,
     &                  thetai,
     &                  vi(3), n(3), t1(3), t2(3), vo(ntp,3),
     &                  voavg(3), ER(ntp), ERavg, vout(3),
     &                  finish_ratio
      integer :: i, j, k, ii
#ifndef LEGACYCOMP
      intrinsic norm2
#endif

cc description of variables

c   v: discretization of incident velocity
c   dtheta: distance between discretization points theta
c   theta: discretization vector incident polar angle
c   Ei: energy of incident particle
c   thetai: polar angle of incident particle
c   vi: velocity of incident particle in vector form
c   n: surface normal unit vector (pointing away from wall, into plasma)
c   t1: toroidal unit vector
c   t2: tangential unit vector (perpendicular to n and t1)
c   vo: velocity vectors of reflected particles
c   vout: velocity vector of a reflected particle after random sampling in CalcVelTRIM
c   voavg: average velocity vector of reflected particles (average over the Monte Carlo samples)
c   ER: energy of reflected particle
c   ERavg: average energy of reflected particles (average over the Monte Carlo samples)
c   finish_ratio: to output how far we are in the calculation
c   i, j, k, ii: loop variables

** documentation

* This routine calculates the reflection, momentum reflection and energy
* reflection coefficients from the TRIM code database

      call subini ('CalcCoefficientsTRIM')

** computation

cc discretize theta
      dtheta = pi/2/ntheta
      do i = 1, ntheta
       theta(i) = dtheta/2 + (i - 1) * dtheta
      end do

cc discretize v
      call SetV(v)

cc allocate arrays
      allocate(Rv(nv,ntheta))
      allocate(Rmomparv(nv,ntheta))
      allocate(Rmomperpv(nv,ntheta))
      allocate(Rpv1(nv,ntheta))
      allocate(Rpv2(nv,ntheta))
      allocate(Rpv3(nv,ntheta))
      allocate(REv(nv,ntheta))

cc Interpolate the reflection data from TRIM to the format we need, i.e. as a function of v and theta
      call InterpolateRN(v,Rv)

      do i = 1,nv
        finish_ratio = 100.0*i/nv
        write(*,*) finish_ratio, '% of Phase 1 done'
        do j = 1,ntheta
          !! Incident energy and angle
          Ei = 0.5*mass*v(i)**2
          thetai = theta(j)
          !! Velocity vector of the incident particle
          vi(1) = sqrt(2*Ei/mass)*sin(thetai)
          vi(2) = 0
          vi(3) = -sqrt(2*Ei/mass)*cos(thetai)
          !! Define normal and tangential directions
          n(1) = 0.0_R8
          n(2) = 0.0_R8
          n(3) = 1.0_R8
          t1(1) = 1.0_R8
          t1(2) = 0.0_R8
          t1(3) = 0.0_R8
          t2(1) = 0.0_R8
          t2(2) = 1.0_R8
          t2(3) = 0.0_R8
          do k = 1,ntp
            call CalcVelTRIM(mass,n,t1,t2,vi,vout)
            !! check for NaNs
            do ii = 1,3
!! should not happen anymore after bug fixes, but can stay here just in case.
#ifndef NO_IEEE
              if (ieee_is_nan(vout(ii))) then
#else
              if (my_isnan(vout(ii))) then
#endif
                write(*,*) 'NaN error in vout!'
                write(*,*) 'k', k
                write(*,*) 'ii', ii
                write(*,*) 'putting value to zero'
                vout(ii) = 0.0_R8
              endif
            enddo
            vo(k,:) = vout
          enddo
          !! take average
          voavg(1) = SUM(vo(:,1))/SIZE(vo(:,1));
          voavg(2) = SUM(vo(:,2))/SIZE(vo(:,2));
          voavg(3) = SUM(vo(:,3))/SIZE(vo(:,3));
          !! energy of reflected particles
          do k = 1,ntp
#ifndef LEGACYCOMP
            ER(k) = 0.5*mass*(NORM2(vo(k,:)))**2
#else
            ER(k) = 0.5*mass*(NORM2(3,vo(k,:)))**2
#endif
          enddo
          !! take average
          ERavg = SUM(ER)/SIZE(ER)
          !! final contributions, see Eqs 4.62, 4.67 and 4.73 in PhD N. Horsten
          Rmomparv(i,j) = Rv(i,j)*voavg(1)/vi(1)
          Rmomperpv(i,j) = -Rv(i,j)*voavg(3)/vi(3)
          REv(i,j) = Rv(i,j)*ERavg/Ei
        enddo
      enddo

cc deterministic way to calculate the Rpv coefficients, because they were much too noisy with the Monte Carl procedure.
      call CalcPressureCoefficientsTRIM(v)

* ..return
      call subend ()
      return
      end subroutine CalcCoefficientsTRIM


      subroutine InterpolateRN(v,RN)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
** input
      real (kind=R8), intent(in) :: v(nv)
** output
      real (kind=R8), intent(out) :: RN(nv,ntheta)
** local variables
      real (kind=R8), allocatable :: A(:), E(:), Ei_mesh(:,:),
     &         Thetai_mesh(:,:), Ei_v(:), Thetai_v(:),
     &         Thetai_v_mesh(:,:), Ei_v_mesh(:,:), RN_EiThetai(:,:),
     &         RN_EiThetai_trans(:,:)
      real (kind=R8) :: dtheta
      integer :: i, j

cc description of variables

c   v: discretization vector of incident velocity
c   RN: reflection probability as a function of v and theta
c   A: discretization vector of the incident polar angle in TRIM
c   E: discretization vector of the incident energy in TRIM
c   Ei_mesh: meshgrid of E
c   Thetai_mesh: meshgrid of A
c   Ei_v: discretization vector of the incident energy corresponding to v
c   Thetai_v: discretization vector of the incident polar angle corresponding to theta
c   Thetai_v_mesh: meshgrid of Thetai_v
c   Ei_v_mesh: meshgrid for Ei_v
c   RN_EiThetai: reflection probability as a function of E and A
c   RN_EiThetai_trans: transpose of RN_EiThetai
c   dtheta: distance between discretization points of theta
c   i,j: loop variables

** documentation
cc This function is used to interpolate the reflection coefficient from the
cc TRIM code database to the discrete velocities and polar angles where the
cc reflection coefficient has to be known.

      call subini ('InterpolateRN')

** computation

cc Energies and angles from the TRIM code database
      allocate(E(nE))
      allocate(A(nA+1))

      do i = 1,nE
        E(i) = TRIMdata((i-1)*nB*nA+1,cE)
      enddo
      do i = 1,nA
        A(i) = TRIMdata((i-1)*nB+1,cA)
      enddo

      A(nA+1) = 90

      allocate(Ei_mesh(nA+1,nE))
      allocate(Thetai_mesh(nA+1,nE))
      call meshgrid(nE,nA+1,E,A,Ei_mesh,Thetai_mesh)

cc Incident energies and polar angles to where the data has to be
cc interpolated
      allocate(Ei_v(nv))
      Ei_v   = mass/2*v**2/ev

      allocate(Thetai_v(ntheta))
      dtheta = pi/2/ntheta
      do i = 1,ntheta
        Thetai_v(i) = ((i-1)*dtheta+dtheta/2)*180/pi
      enddo

      allocate(Thetai_v_mesh(ntheta,nv))
      allocate(Ei_v_mesh(ntheta,nv))
      call meshgrid(ntheta,nv,Thetai_v,Ei_v,Thetai_v_mesh,Ei_v_mesh)

cc Particle reflection coefficient as a function of the incident energy and the polar angle from the TRIM
cc code database
      allocate(RN_EiThetai(nE,nA+1))

      do i= 1,nE
        do j = 1,nA
          RN_EiThetai(i,j) = TRIMdata(nB*nA*(i-1)+nB*(j-1)+1,cR)
        enddo
        RN_EiThetai(i,nA+1) = 1
      enddo
cc take transpose to match the expected format of interp2
      allocate(RN_EiThetai_trans(nA+1,nE))
      RN_EiThetai_trans = TRANSPOSE(RN_EiThetai)

      call interp2(nA+1,nE,Ei_mesh,Thetai_mesh,RN_EiThetai_trans,
     &                       nv,ntheta, Ei_v_mesh,Thetai_v_mesh,RN)

cc If the incident energy is smaller than the minimum energy, linear
cc extrapolation is used
      call extrapolate_RN_minenergy(nv,ntheta,RN,Ei_v,
     &     nE, nA+1, E, A, RN_EiThetai,Thetai_v)

cc If the incident energy is larger than the maximum energy, linear
cc extrapolation is used
      call extrapolate_Rpv_maxenergy(nv,ntheta,RN,Ei_v,
     &     nE, nA+1, E, A, RN_EiThetai,Thetai_v)

      do i=1,nv
        do j=1,ntheta
          RN(i,j) = max(0.0_R8,RN(i,j))
          RN(i,j) = min(1.0_R8,RN(i,j))
        enddo
      enddo

cc deallocate
      deallocate(E)
      deallocate(A)
      deallocate(Ei_v)
      deallocate(Thetai_v)
      deallocate(Ei_v_mesh)
      deallocate(Thetai_v_mesh)
      deallocate(Ei_mesh)
      deallocate(Thetai_mesh)
      deallocate(RN_EiThetai)
      deallocate(RN_EiThetai_trans)

* ..return
      call subend ()
      return
      end subroutine InterpolateRN

      subroutine CalcVelTRIM(mass,n,t1,t2,v,vout)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
#ifdef LEGACYCOMP
      use b2mod_math
#endif
      implicit none

** input
      real (kind=R8), intent(in) :: mass, n(3), t1(3), t2(3),
     &                v(3)
** output
      real (kind=R8), intent(out) :: vout(3)
** local variables
      real (kind=R8), allocatable :: E(:), A(:)
      integer :: i
      real (kind=R8) :: theta, Ei, phi, dotvn, dotvt2, dotvt1
      integer :: jE, iE, jA, iA, jEo, iEo, iTo, jTo, iPo, jPo
      real (kind=R8) :: wE1, wE2, wA1, wA2, wEo1, wEo2, wTo1, wTo2,
     &                  wPo1, wPo2, temprand
      integer :: index11, index12, index21, index22
      real (kind=R8) :: P11(31,5), P12(31,5), P21(31,5), P22(31,5),
     &                  Q(5), Et(5), Tt(5,5), Ttt(5),
     &                  Phit(25,5), Phitt(5,5), Phittt(5),
     &                  Eo, To, Phio, Es, Ts, Phis, vo,
     &                  vn(3), vt1(3), vt2(3), P(31,5)
#ifndef LEGACYCOMP
      intrinsic norm2
#endif

cc description of variables

c   mass: mass of the incident particle (kg)
c   n: surface normal unit vector (pointing away from wall, into plasma)
c   t1: toroidal unit vector
c   t2: tangential unit vector (perpendicular to n and t1)
c   v: discretization vector of incident velocity

c   vout: velocity vector of a reflected particle

c   A: discretization vector of the incident polar angle in TRIM
c   E: discretization vector of the incident energy in TRIM

c   i: loop variable

c   theta: polar angle incident particle
c   Ei: energy incident particle
c   phi: azimuthal angle of reflected particle
c   dotvn: dot product of v and n
c   dotvt2: dot product of v and t2
c   dotvt1: dot product of v and t1

c   jE: index of point in E just after Ei
c   iE: index of point in E just before Ei
c   jA: index of point in A just after theta
c   iA: index of point in A just before theta
c   jEo, iEo: indices to find reflected energies corresponding to random number Es
c   iTo, jTo: indices to find reflected polar angle corresponding to random number Ts
c   iPo, jPo: indices to find reflected azimuth angle corresponding to random number Ts
c   wE1, wE2: interpolation weights related to incident energy
c   wA1, wA2:: interpolation weights related to incident polar angle
c   wEo1, wEo2: interpolation weights related to reflected energy
c   wTo1, wTo2: interpolation weights related to reflected polar angle
c   wPo1, wPo2: interpolation weights related to reflected azimuth angle
c   temprand: temporary variable to assign random number

c   index11, index12, index21, index22: indices to find the TRIM reflection data for the points surrounding (v,theta)

c   P11, P12, P21, P22 reflection data from TRIM for the points surrounding (v,theta)

c   Q: quantile matrix (see line just below Eq. 4.19 in PhD N. Horsten)
c   Et: probability distribution for reflected energy following for a given incident v and theta
c   Tt: conditional probability distribution for reflected polar angle for a given incident v and theta, and depending on the reflected energy
c   Ttt: conditional probability distribution for reflected polar angle for a given incident v and theta, knwowing the reflected energy
c   Phit: conditional probability distribution for reflected azimuth angle for a given incident v and theta, and depending on the reflected energy and reflected polar angle
c   Phittt: conditional probability distribution for reflected azimuth angle for a given incident v and theta, and depending on the reflected polar angle, knowing the reflected energy
c   Phitttt: conditional probability distribution for reflected azimuth angle for a given incident v and theta, knowing the reflected energy and reflected polar angle
c   Eo: energy of reflected particle
c   To: polar angle of reflected particle
c   Phio: azimuth angle of reflected particle
c   Es: random number for energy
c   Ts: random number for polar angle (theta)
c   Phis: random number of azimuth angle (phi)
c   vo: speed of reflected particle
c   vn: normal component of reflected particle
c   vt1: toroidal component of reflected particle
c   vt2: tangential component of reflected particle
c   P: reflection data, interpolated from TRIM to incident particle with certain incident v and theta

** Documentation
* This routine calculates the speed of the neutral if it is reflected as a
* fast neutral as described by the TRIM code database. To understand this routine,
* you must understand how the TRIM database from EIRENE is read. This is described in for example
* Section 4.1.2 of the PhD thesis of Niels Horsten (KU Leuven, 2019). More sp

      call subini ('CalcVelTRIM')
** Computation

cc Allocate
      allocate(E(nE))
      allocate(A(nA))

      do i=1,nE
        E(i) = TRIMdata((i-1)*nB*nA+1,cE)
      enddo

      do i=1,nA
        A(i) = TRIMdata((i-1)*nB+1,cA)
      enddo
cc Calculate angle and energy of incidence
      dotvn = dot_product(v,n)
#ifndef LEGACYCOMP
      theta = acos(-dotvn/NORM2(v))
#else
      theta = acos(-dotvn/NORM2(3,v))
#endif
      theta = theta/pi*180.0_R8
#ifndef LEACYCOMP
      Ei = 0.5*mass*NORM2(v)**2/eV
#else
      Ei = 0.5*mass*NORM2(3,v)**2/eV
#endif
      dotvt2 = dot_product(v,t2)
      dotvt1 = dot_product(v,t1)
      phi = atan(dotvt2/dotvt1)
      if (dotvt1.le.0) then
        phi = phi + pi
      endif
cc Determine weight coefficients and reflection coefficients and parameters
cc P contains reflection parameters
      if (Ei.lt.E(1)) then
        jE = 2
        iE = 1
        wE1 = (Ei-E(jE))/(E(jE)-E(iE))+1
        wE2 = -(Ei-E(jE))/(E(jE)-E(iE))
        if (theta.gt.A(nA)) then
          iA = SIZE(A)
          wA1 = (theta-A(iA))/(A(iA)-A(iA-1))+1
          wA2 = -(theta-A(iA))/(A(iA)-A(iA-1))
          index11 = ((jE-1)*nA+iA-1)*nB
          index12 = ((jE-1)*nA+iA-2)*nB
          index21 = ((iE-1)*nA+iA-1)*nB
          index22 = ((iE-1)*nA+iA-2)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        else
          call search(theta,nA,A,iA,jA)
          wA1 = (theta-A(iA))/(A(jA)-A(iA))
          wA2 = (A(jA)-theta)/(A(jA)-A(iA))
          index11 = ((jE-1)*nA+jA-1)*nB
          index12 = ((jE-1)*nA+iA-1)*nB
          index21 = ((iE-1)*nA+jA-1)*nB
          index22 = ((iE-1)*nA+iA-1)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        endif
      elseif (Ei.gt.E(nE)) then
        jE = SIZE(E)
        iE = SIZE(E)-1
        if (theta.gt.A(nA)) then
          iA = SIZE(A)
          wA1 = (theta-A(iA))/(A(iA)-A(iA-1))+1
          wA2 = -(theta-A(iA))/(A(iA)-A(iA-1))
          index11 = ((jE-1)*nA+iA-1)*nB
          index12 = ((jE-1)*nA+iA-2)*nB
          index21 = ((iE-1)*nA+iA-1)*nB
          index22 = ((iE-1)*nA+iA-2)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        else
          call search(theta,nA,A,iA,jA)
          wA1 = (theta-A(iA))/(A(jA)-A(iA))
          wA2 = (A(jA)-theta)/(A(jA)-A(iA))
          index11 = ((jE-1)*nA+jA-1)*nB
          index12 = ((jE-1)*nA+iA-1)*nB
          index21 = ((iE-1)*nA+jA-1)*nB
          index22 = ((iE-1)*nA+iA-1)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        endif
        wE1 = (Ei-E(jE))/(E(jE)-E(iE))+1
        wE2 = -(Ei-E(jE))/(E(jE)-E(iE))
      else
        call search(Ei,nE,E,iE,jE)
        wE1 = (Ei-E(iE))/(E(jE)-E(iE))
        wE2 = (E(jE)-Ei)/(E(jE)-E(iE))
        if (theta.gt.A(nA)) then
          iA = SIZE(A)
          wA1 = (theta-A(iA))/(A(iA)-A(iA-1))+1
          wA2 = -(theta-A(iA))/(A(iA)-A(iA-1))
          index11 = ((jE-1)*nA+iA-1)*nB
          index12 = ((jE-1)*nA+iA-2)*nB
          index21 = ((iE-1)*nA+iA-1)*nB
          index22 = ((iE-1)*nA+iA-2)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        else
          call search(theta,nA,A,iA,jA)
          wA1 = (theta-A(iA))/(A(jA)-A(iA))
          wA2 = (A(jA)-theta)/(A(jA)-A(iA))
          index11 = ((jE-1)*nA+jA-1)*nB
          index12 = ((jE-1)*nA+iA-1)*nB
          index21 = ((iE-1)*nA+jA-1)*nB
          index22 = ((iE-1)*nA+iA-1)*nB
          P11 = TRIMdata(index11+2:index11+32,1:5)
          P12 = TRIMdata(index12+2:index12+32,1:5)
          P21 = TRIMdata(index21+2:index21+32,1:5)
          P22 = TRIMdata(index22+2:index22+32,1:5)
        endif
      endif
cc Interpolate to obtain reflection coefficient and weighted parameter matrix
      P = wE1*(wA1*P11+wA2*P12)+wE2*(wA1*P21+wA2*P22)
cc define quantile matrix
      Q(1) = 0.1_R8
      Q(2) = 0.3_R8
      Q(3) = 0.5_R8
      Q(4) = 0.7_R8
      Q(5) = 0.9_R8

cc sample Eout, Thetaout and Phiout (s is samples, o is out, t is temp)
      call random_number(Es)
      call random_number(Ts)
      call random_number(Phis)

      Et = P(1,:)
      Tt = P(2:6,:)
      Phit = P(7:,:)
cc note: TRIM-code contains the cosine of the angles, which is what will result here as well

cc The energy of the reflected particle has to be larger than 0 eV
      Et = max(Eomin,Et)

cc The energy of the reflected particle has to be smaller than the incident energy
      Et = min(Et,Ei)

cc Calculate weights for Eo
      if (Es.gt.Q(5)) then !! Maximum value of the reflected energy is the incident energy
        jEo = SIZE(Q)
        wEo1 = (1-Es)/0.1_R8
        wEo2 = (Es-Q(jEo))/0.1_R8
        Eo = wEo1*Et(jEo) + WEo2*Ei
      elseif (Es.lt.Q(1)) then  !! Minimum value of the reflected energy is 0 eV
        iEo = 1
        wEo2 = Es/Q(1)
        Eo = wEo2*Et(iEo)
      else
        call search(Es,5,Q,iEo,jEo)
        wEo1 = (Es-Q(iEo))/(Q(jEo)-Q(iEo));
        wEo2 = (Q(jEo)-Es)/(Q(jEo)-Q(iEo));
        Eo = wEo1*Et(jEo)+wEo2*Et(iEo);
      endif

cc Determine the quantiles of the angles corresponding to Eo
      if (Es.gt.Q(5)) then
        jEo = SIZE(Q)
        iEo = jEo-1
        wEo1 = (Es-Q(jEo))/(Q(jEo)-Q(iEo)) + 1
        wEo2 = -(Es-Q(jEo))/(Q(jEo)-Q(iEo))
      elseif (Es.lt.Q(1)) then
        jEo = 2
        iEo = 1
        wEo1 = (Es-Q(jEo))/(Q(jEo)-Q(iEo)) + 1
        wEo2 = -(Es-Q(jEo))/(Q(jEo)-Q(iEo))
      else
        call search(Es,5,Q,iEo,jEo)
        wEo1 = (Es-Q(iEo))/(Q(jEo)-Q(iEo))
        wEo2 = (Q(jEo)-Es)/(Q(jEo)-Q(iEo))
      endif
      Ttt   = wEo1*Tt(jEo,:) + wEo2*Tt(iEo,:)
      Phitt = wEo1*Phit((jEo-1)*5+1:(jEo-1)*5+5,:) +
     &        wEo2*Phit((iEo-1)*5+1:(iEo-1)*5+5,:)

cc The minimum value of To is 0
      Ttt = max(Ttt,0.0_R8)
cc The maximum value of To is 1
      Ttt = min(Ttt,1.0_R8)

cc Calculate weights for To
      if (Ts.gt.Q(5)) then
        jTo  = SIZE(Q)
        iTo  = jTo-1
        wTo1 = (Ts-Q(jTo))/(Q(jTo)-Q(iTo))+1
        wTo2 = -(Ts-Q(jTo))/(Q(jTo)-Q(iTo))
        !! interpolate To
        To = (Ts-Q(5))/(1-Q(5))*1+(1-Ts)/(1-Q(5))*Ttt(5)
      elseif(Ts.lt.Q(1)) then
        jTo  = 2
        iTo  = 1
        wTo1 = (Ts-Q(jTo))/(Q(jTo)-Q(iTo))+1
        wTo2 = -(Ts-Q(jTo))/(Q(jTo)-Q(iTo))
        !! interpolate To
        To = Ts/Q(1)*Ttt(1);
      else
        call search(Ts,5,Q,iTo,jTo)
        wTo1      = (Ts-Q(iTo))/(Q(jTo)-Q(iTo))
        wTo2      = (Q(jTo)-Ts)/(Q(jTo)-Q(iTo))
        !! interpolate To
        To   = wTo1*Ttt(jTo) + wTo2*Ttt(iTo)
      endif

cc interpolate Phio
      Phittt = wTo1*Phitt(jTo,:)+wTo2*Phitt(iTo,:)

cc The minimum value of Phio is -1
      Phittt = max(Phittt,-0.999_R8) !! (setting this to -1.0 sometimes gives NaN errors when taking the acos later on ...)

cc The maximum value of To is 1
      Phittt = min(Phittt,1.0_R8)

cc Calculate weights for Phio
      if(Phis.gt.Q(5)) then
cc     interpolate phi:
        Phio = (Phis-Q(5))/(1-Q(5))*1+(1-Phis)/(1-Q(5))*Phittt(5)
      elseif(Phis.lt.Q(1)) then
cc     interpolate phi:
       Phio = (Q(1)-Phis)/Q(1)*(-1)+Phis/Q(1)*Phittt(1)
      else
        call search(Phis,5,Q,iPo,jPo)
        wPo1      = (Phis-Q(iPo))/(Q(jPo)-Q(iPo))
        wPo2      = (Q(jPo)-Phis)/(Q(jPo)-Q(iPo))
cc     interpolate phi:
        Phio = wPo1*Phittt(jPo)+wPo2*Phittt(iPo)
      endif

      Phio = acos(Phio)
      call random_number(temprand)
      if (temprand.lt.0.5_R8) then
        Phio = -Phio
      endif

cc For perfect normal incidence, the distribution of the azimuthal angle is uniform between 0 and 2*pi.
      if (dotvt1.eq.0.and.dotvt2.eq.0) then
        Phio = Phis*2*pi
        phi = 0.0_R8
      endif

cc Calculate outgoing velocity components
      vo  = sqrt(2*Eo*eV/mass)
      vn  = vo*To*n
      if (To.gt.1.0_R8) then
        write(*,*) 'Error: To gt 1 !!!'
        write(*,*) 'To =', To
      endif
      vt1 = vo*sqrt(1-To**2)*cos(phi + Phio)*t1
      vt2 = vo*sqrt(1-To**2)*sin(phi + Phio)*t2
      vout = vn + vt1 + vt2

cc deallocate
      deallocate(A)
      deallocate(E)
* ..return
      call subend ()
      return
      end subroutine CalcVelTRIM


      subroutine CalcPressureCoefficientsTRIM(v)
      use b2mod_constants
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      real (kind=R8), intent(in) :: v(nv)
** local
      integer :: iE, iA, i , j
      real (kind=R8), allocatable :: Rp1(:,:), Rp2(:,:), Rp3(:,:),
     &                Rp1trans(:,:), Rp2trans(:,:), Rp3trans(:,:),
     &                   E(:),A(:),Ei_mesh(:,:), Thetai_mesh(:,:)
      real (kind=R8) :: ER, vR, thetaR(5), thetaRinv(5), thetaRj,
     &                  phiR(5), Ei_v(nv), dtheta, Thetai_v(ntheta),
     &                  Thetai_vv(ntheta,nv), Ei_vv(ntheta,nv),
     &                  fac

cc description of variables

c   v: discretization of incident particle velocity

c   iE: index in E
c   iA: index in A
c   i, j: loop variables
c   Rp1: same as Rpv1, but based on E and A instead of v and theta
c   Rp2: same as Rpv2, but based on E and A instead of v and theta
c   Rp3: same as Rpv3, but based on E and A instead of v and theta
c   Rp1trans: transpose of Rp1
c   Rp2trans: transpose of Rp2
c   Rp3trans: transpose of Rp3
c   E: discretization vector of the incident energy in TRIM
c   A: discretization vector of the incident polar angle in TRIM
c   Ei_mesh: meshgrid base on E
c   Thetai_mesh: meshgrid based on A

c   ER: energy of reflected particles
c   vR: velocity of reflected particles (m/s)
c   thetaR: cosine of polar angle of reflected particles (TRIM tabulates the cosine of the angle, not the angle itself)
c   thetaRinv: inverse of thetaR
c   thetaRj: polar angle of reflected particles
c   phiR: azimuthal angle of reflected particles
c   Ei_v: energies corresponding to v
c   dtheta: distance between discretization points theta
c   Thetai_v: discretization vector for theta
c   Thetai_vv: meshgrid of Thetai_v
c   Ei_vv: meshgrid of Ei_v
c   fac: temporary variable for scaling factors

** Documentation
* This function calculates the pressure reflection coefficients from the
* TRIM code database in a deterministic way due to issues when the
* coefficients are calculated with a Monte Carlo approach.
      call subini ('CalcPressureCoefficientsTRIM')

** Computation

cc Pressure reflection coefficients for the energy and angles given in the database
      allocate(Rp1(nE,nA))
      allocate(Rp2(nE,nA))
      allocate(Rp3(nE,nA))
      allocate(Rp1trans(nA,nE))
      allocate(Rp2trans(nA,nE))
      allocate(Rp3trans(nA,nE))

      Rp1 = 0.0_R8
      Rp2 = 0.0_R8
      Rp3 = 0.0_R8
      Rp1trans = 0.0_R8
      Rp2trans = 0.0_R8
      Rp3trans = 0.0_R8

      do iE = 1,nE
        do iA = 1,nA
          do i = 1,5
            ER = TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 2,i)
            vR = sqrt(2*ER*eV/mass)
            thetaR = TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 2 + i,1:5)
            call inverseVector(thetaR,thetaRinv)
            Rp1(iE,iA) = Rp1(iE,iA) +
     &                   1.0/5.0*vR*SUM(thetaRinv)/SIZE(thetaRinv)
            Rp2(iE,iA) = Rp2(iE,iA) +
     &                   1.0/5.0*SUM(thetaRinv)/SIZE(thetaRinv)/vR
            do j = 1,5
              thetaRj = acos(thetaR(j))
              phiR  =
     &        TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 7 + (i-1)*5 + j,1:5)
              Rp3(iE,iA) = Rp3(iE,iA) +
     &                     1.0/25.0*tan(thetaRj)*SUM(phiR)/SIZE(phiR)
            enddo
          enddo
          Rp1(iE,iA) = Rp1(iE,iA)*
     &                 TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 1,cR)
          Rp2(iE,iA) = Rp2(iE,iA)*
     &                 TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 1,cR)
          Rp3(iE,iA) = Rp3(iE,iA)*
     &                 TRIMdata((iE-1)*nA*nB + (iA-1)*nB + 1,cR)
        enddo
      enddo
cc Interpolation to v and theta
cc Incident energies and polar angles to where the data has to be interpolated
      Ei_v = mass/2*v**2/eV
      dtheta = pi/2/ntheta
      do i = 1,ntheta
        Thetai_v(i) = ((i-1)*dtheta+dtheta/2.0)*180.0/pi
      enddo

      call meshgrid(ntheta,nv,Thetai_v,Ei_v,Thetai_vv,Ei_vv)

      allocate(E(nE))
      allocate(A(nA))

      do i = 1,nE
        E(i) = TRIMdata((i-1)*nB*nA+1,cE)
      enddo
      do i = 1,nA
        A(i) = TRIMdata((i-1)*nB+1,cA)
      enddo

      allocate(Ei_mesh(nA,nE))
      allocate(Thetai_mesh(nA,nE))
      call meshgrid(nE,nA,E,A,Ei_mesh,Thetai_mesh)

      Rp1trans = TRANSPOSE(Rp1)
      Rp2trans = TRANSPOSE(Rp2)
      Rp3trans = TRANSPOSE(Rp3)
      call interp2(nA,nE,Ei_mesh,Thetai_mesh,Rp1trans,
     &                        nv,ntheta,Ei_vv,Thetai_vv,Rpv1)
      call interp2(nA,nE,Ei_mesh,Thetai_mesh,Rp2trans,
     &                        nv,ntheta,Ei_vv,Thetai_vv,Rpv2)
      call interp2(nA,nE,Ei_mesh,Thetai_mesh,Rp3trans,
     &                        nv,ntheta,Ei_vv,Thetai_vv,Rpv3)

cc Extrapolation
cc If the incident angle is larger than the maximum angle, linear extrapolation is used.
      call extrapolate_Rpv_angle(nv, ntheta, Rpv1, Thetai_v,
     &                                    nA, nE, A, E, Rp1, Ei_v)
      call extrapolate_Rpv_angle(nv, ntheta, Rpv2, Thetai_v,
     &                                    nA, nE, A, E, Rp2, Ei_v)
      call extrapolate_Rpv_angle(nv, ntheta, Rpv3, Thetai_v,
     &                                    nA, nE, A, E, Rp3, Ei_v)

cc If the incident energy is smaller than the minimum energy, linear extrapolation is used
      call extrapolate_Rpv_minenergy(nv, ntheta, Rpv1, Ei_v,
     &                                    nE, nA, E, A, Rp1, Thetai_v)
      call extrapolate_Rpv_minenergy(nv, ntheta, Rpv2, Ei_v,
     &                                    nE, nA, E, A, Rp2, Thetai_v)
      call extrapolate_Rpv_minenergy(nv, ntheta, Rpv3, Ei_v,
     &                                    nE, nA, E, A, Rp3, Thetai_v)

cc If the incident energy is larger than the maximum energy, linear extrapolation is used
      call extrapolate_Rpv_maxenergy(nv, ntheta, Rpv1, Ei_v,
     &                                    nE, nA, E, A, Rp1, Thetai_v)
      call extrapolate_Rpv_maxenergy(nv, ntheta, Rpv2, Ei_v,
     &                                    nE, nA, E, A, Rp2, Thetai_v)
      call extrapolate_Rpv_maxenergy(nv, ntheta, Rpv3, Ei_v,
     &                                    nE, nA, E, A, Rp3, Thetai_v)

      Rpv1 = max(0.0_R8,Rpv1)
      Rpv2 = max(0.0_R8,Rpv2)

cc Scaling
      do i = 1,nv
        do j = 1,ntheta
          fac = cos(Thetai_v(j)*pi/180)/(2*Ei_v(i)*eV/mass)**(0.5)
          Rpv1(i,j) = Rpv1(i,j)*fac

          fac = cos(Thetai_v(j)*pi/180)*(2*Ei_v(i)*eV/mass)**(0.5)
          Rpv2(i,j) = Rpv2(i,j)*fac

          fac = 1/(tan(Thetai_v(j)*pi/180))
          Rpv3(i,j) = Rpv3(i,j)*fac
        enddo
      enddo

      deallocate(E)
      deallocate(A)
      deallocate(Ei_mesh)
      deallocate(Thetai_mesh)

* ..return
      call subend ()
      return
      end subroutine CalcPressureCoefficientsTRIM

c**************************************************************************************
c*********   Phase 2: Integrate response over incident distributions        ***********
c**************************************************************************************

      subroutine Discretize_M_T_Sh()
      use b2mod_types
      use b2mod_subsys
      implicit none
** local
      integer :: i
      real (kind=R8) :: dSh, logt, dlogt, logTmin, logTmax

cc description of variables:

c   i: loop variable
c   dSh: (linear) distance between Sh values
c   logt: log10 of T
c   dlogt: distance between logt values
c   logTmin: log10 of Tmin
c   logTmax: log10 of Tmax

      call subini ('Discretize_M_T_Sh')
cc create the discretization vectors Tv (temperature) and Mv (toroidal mach number) based on the minimum, maximum, and number of discretization points

cc allocate arrays
      allocate(M(nM))
      allocate(T(nT))
      allocate(Sh(nSh))

cc discretize the toroidal Mach number
      dM = (Mmax - Mmin)/nM
      do i = 1,nM
        M(i) = Mmin + dM/2 + (i-1)*dM
      enddo

cc discretize the temperature
      logTmin = log10(Tmin)
      logTmax = log10(Tmax)
      dlogt = (logTmax-logTmin)/nT
      logt = logTmin
      do i = 1,nT
        T(i) = 10**logt
        logt = logt + dlogt
      enddo

cc discretize the sheath potential drop
      dSh = (Shmax-Shmin)/(nSh-1)
      do i = 1,nSh
        Sh(i) = Shmin + (i-1)*dSh
      enddo

* ..return
      call subend ()
      return
      end subroutine Discretize_M_T_Sh


      subroutine CalcFluidReflCoefficients()
      use b2mod_types
      use b2mod_subsys
      implicit none
** local
      integer :: i, j
      real (kind=R8) :: v(nv), theta(ntheta), dtheta,
     &                  dv(nv), deltaP

cc description of variables

c   i,j: loop variables
c   theta: discretization vector polar angle
c   dtheta: distance between discretization points theta
c   dv: distance between discretization points v
c   deltaP: factor for sheath potential drop

      call subini ('CalcFluidReflCoefficients')
cc This function is used to interpolate the TRIM reflection coefficients to
cc the discrete velocities for a certain value of the sheath potential.

cc computation

cc create discretization for M and T
      call Discretize_M_T_Sh()
      call setV(v)
      do i = 1,nv
        if (i.gt.1.and.i.lt.nv) then
          dv(i) = (v(i+1)-v(i-1))/2.0_R8
        endif
      enddo
      dv(1) = v(2) - v(1)
      dv(nv) = v(nv) - v(nv-1)

cc discretize theta
      dtheta = pi/2/ntheta
      do i = 1, ntheta
       theta(i) = dtheta/2 + (i - 1) * dtheta
      end do

cc allocate arrays
      allocate(R1(nM,nT,nSh))
      allocate(R2(nM,nT,nSh))
      allocate(Rmx1(nM,nT,nSh))
      allocate(Rmx2(nM,nT,nSh))
      allocate(Rmy1(nM,nT,nSh))
      allocate(Rmy2(nM,nT,nSh))
      allocate(Rmz1(nM,nT,nSh))
      allocate(Rmz2(nM,nT,nSh))
      allocate(Rp11(nM,nT,nSh))
      allocate(Rp12(nM,nT,nSh))
      allocate(Rp21(nM,nT,nSh))
      allocate(Rp22(nM,nT,nSh))
      allocate(Rp31(nM,nT,nSh))
      allocate(Rp32(nM,nT,nSh))
      allocate(RE1(nM,nT,nSh))
      allocate(RE2(nM,nT,nSh))

      do j=1,nSh
        deltaP = Sh(j)
        write(*,*) 'performing sheath potential', j, 'of', nSh
        do i=1,nT
          write(*,*) 'performing temperature', i, 'of', nT
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rv,T(i),deltaP,'particle',R1(:,i,j),R2(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rmomparv,T(i),deltaP,'momx',Rmx1(:,i,j),Rmx2(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rmomparv,T(i),deltaP,'momy',Rmy1(:,i,j),Rmy2(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &     Rmomperpv,T(i),deltaP,'momz',Rmz1(:,i,j),Rmz2(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      REv,T(i),deltaP,'energy',RE1(:,i,j),RE2(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rpv1,T(i),deltaP,'press1',Rp11(:,i,j),Rp12(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rpv2,T(i),deltaP,'press2',Rp21(:,i,j),Rp22(:,i,j))
          call CalcReflectionCoefficientsFluidTRIM(v,dv,theta,
     &      Rpv3,T(i),deltaP,'press3',Rp31(:,i,j),Rp32(:,i,j))
        enddo
      enddo

* ..return
      call subend ()
      return
      end subroutine CalcFluidReflCoefficients

      subroutine CalcReflectionCoefficientsFluidTRIM(
     &           vin,dvin,theta,R,Tin,deltaP,mom,R1,R2)
      use b2mod_types
      use b2mod_subsys
      use b2mod_xerset
      implicit none
** input
cc      integer, intent(in) ::
      real (kind=R8), intent(in) :: vin(nv), dvin(nv),
     &         theta(ntheta), Tin, R(nv,ntheta), deltaP
      character(len=*), intent(in), optional :: mom
** output
      real (kind=R8), intent(out) :: R1(nM), R2(nM)
** local
      integer :: i, k, nthetai, ithetai, j, l, iv, nvl
      real (kind=R8) :: vmin, phi(nphi), dphi, cs, u,
     &              thetamax, dtheta, tref
      real (kind=R8), allocatable :: thetai(:), dthetai(:),
     &          A1(:,:), A2(:,:), int_temp1(:), int_temp2(:),
     &          dv(:), v(:), int1(:), int2(:)

cc description of variables

c   vin: input discretization vector v (can be altered)
c   dvin: distance between points of vin
c   thetain: input discretization vector theta (can be altered)
c   Tin: input ion temperature
c   R: either Rv, REv, Rmomparv, Rmomperpv, Rp1, Rp2, or Rp3, depending on 'mom'
c   deltaP: factor for sheath potential drop
c   mom: the reflection coefficient that will be calculated

c   R1: either R1, Rmx1, (Rmy1 is zero), Rmz1, RE1, Rp11, Rp21, or Rp31, depending on 'mom'
c   R2: either R2, Rmx2, Rmy2, Rmz2, RE2, Rp12, Rp22, or Rp32, depending on 'mom'

c   i, k: loop variables
c   nthetai: number of elements in thetai
c   j, l: loop variables
c   iv: helper index to find v from vin
c   nvl: lenght of v
c   vmin: minimum velocity of an incident particle (if it was only accelerated by sheath)
c   phi: discretization of azimuth angle
c   dphi: distance between points of phi
c   cs: sound speed corresponding to Ti
c   u: toroidal bulk velocity in incident particles
c   thetamax:
c   detha:
c   tref: Tin expressed in joule instead of eV
c   thetai: thatain, but with the values below thetamax discarded
c   dthetai: distance between points of thetai
c   A1: integrand to be integrated over phi and theta, related to R1
c   A2: integrand to be integrated over phi and theta, related to R2
c   int_temp1: temporary variable in integration
c   int_temp2: temporary variable in integration
c   v: vin, but with the values below vmin discarded
c   dv: distances between points of v
c   int1: integrand for integration over v, related to R1
c   int2: integrand for integration over v, related to R2

** documentation
c This function calculates the dimensionless reflection coefficients from
c the TRIM code database needed for the fluid neutral models.
c nM = number of Mach numbers for which the reflection coefficients are
c calculated

      call subini ('CaldReflectionCoefficientsFluidTRIM')

      tref = Tin*eV
      vmin = sqrt(2*deltaP*tref/mass)

      if (vmin.le.vin(size(vin)-1)) then
        !Find the index of the first element in v that is greater than or equal to vmin
        do i=1,nv
          iv = i
          if(vin(i).ge.vmin) exit
        enddo
        allocate(v(size(vin)-iv+1))
        allocate(dv(size(vin)-iv+1))
        v = vin(iv:size(vin)) ! Remove elements from v that are less than or equal to vmin
        dv = dvin(size(dvin) - size(v) + 1:) ! Trim dv to match the size of the updated v
        dv(1) = (v(1) + v(2))/2 - vmin
      else
        iv = size(vin) - 1 ! Set iv as an array with the value [length(v) - 1]
        v = vin(size(vin) - 1:size(vin)) ! Extract the last two elements of v
        dv = dvin(size(dvin) - size(v) + 2:) ! Trim dv to match the size of the updated v
        dv(1) = (v(1) + v(2))/2 - vmin
      end if
      nvl = size(v)
      allocate(int1(nvl))
      allocate(int2(nvl))

      dtheta = pi/2/ntheta
      dphi = 2*pi/nphi
      do i=1,nphi
        phi(i) = dphi/2+(i-1)*dphi
      enddo

c Speed of sound for reference temperature
      cs = sqrt(2*tref/mass)

      do i=1,nM
        u = M(i)*cs

        int1 = 0.0_R8
        int2 = 0.0_R8
        do k = 1,nvl
          if (vmin.le.v(k)) then
            thetamax = acos(vmin/v(k))
          else
            thetamax = 0.0_R8
          endif
cc make a thetai array which contains only the values smaller than thetamax
          nthetai = 0
          do j = 1, ntheta
            if (theta(j).le.thetamax) then
              nthetai = nthetai + 1
            endif
          enddo
          allocate(thetai(nthetai))
          allocate(dthetai(nthetai))
          ithetai = 0
          do j = 1, ntheta
            if (theta(j).le.thetamax) then
              ithetai = ithetai + 1
              thetai(ithetai) = theta(j)
            endif
          enddo
          if (ntheta.ne.0) then
            if (size(thetai).gt.1) then
              dthetai(:)   = dtheta
              dthetai(nthetai) = thetamax -
     &        (thetai(nthetai) + thetai(nthetai-1))/2
            else
              dthetai(1) = thetamax
            endif
            allocate(A1(nthetai,nphi))
            allocate(A2(nthetai,nphi))
            allocate(int_temp1(nthetai))
            allocate(int_temp2(nthetai))

            select case (mom)
              case ('particle')
                do j = 1,nthetai
                  do l=1,nphi
                    A1(j,l) = v(k)**3*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     cos(thetai(j))*sin(thetai(j))
                    A2(j,l) = v(k)**4*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))**2*sin(thetai(j))
                  enddo
                enddo
              case ('momx')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**4*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     cos(thetai(j))*sin(thetai(j))**2*cos(phi(l))
                    A2(j,l) = v(k)**5*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                 (cos(thetai(j))*sin(thetai(j)))**2*cos(phi(l))
                  enddo
                enddo
              case ('momy')
                do j = 1,nthetai
                  do l = 1,nphi
                    A2(j,l) = v(k)**5*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                 cos(thetai(j))*sin(thetai(j))**3*sin(phi(l))**2
                  enddo
                enddo
              case ('momz')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**4*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     cos(thetai(j))**2*sin(thetai(j))
                    A2(j,l) = v(k)**5*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))**3*sin(thetai(j))
                  enddo
                enddo
              case ('energy')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**5*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     cos(thetai(j))*sin(thetai(j))
                    A2(j,l) = v(k)**6*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))**2*sin(thetai(j))
                  enddo
                enddo
              case ('press1')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**4*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     sin(thetai(j))
                    A2(j,l) = v(k)**5*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))*sin(thetai(j))
                  enddo
                enddo
              case ('press2')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**2*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     sin(thetai(j))
                    A2(j,l) = v(k)**3*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))*sin(thetai(j))
                  enddo
                enddo
              case ('press3')
                do j = 1,nthetai
                  do l = 1,nphi
                    A1(j,l) = v(k)**3*exp(-mass*(v(k)**2-2*v(k)*u*
     &                     sin(thetai(j))*cos(phi(l)) + u**2 -
     &                     2*deltaP*tref/mass)/(2*tref))*
     &                     sin(thetai(j))**2*cos(phi(l))
                    A2(j,l) = v(k)**4*exp(-mass*(v(k)**2 - 2*v(k)*u*
     &                    sin(thetai(j))*cos(phi(l)) + u**2 -
     &                    2*deltaP*tref/mass)/(2*tref))*
     &                    cos(thetai(j))*sin(thetai(j))**2*cos(phi(l))
                  enddo
                enddo

              case default
                call xerrab('Error: momentum type not found')
            end select
            int_temp1 = R(iv+k-1,1:nthetai)*dthetai
            int_temp2 = matmul(int_temp1,A1)
            int1(k) = sum(int_temp2)*dphi
            int_temp2 = matmul(int_temp1,A2)
            int2(k) = sum(int_temp2)*dphi
          endif

        deallocate(thetai)
        deallocate(dthetai)
        deallocate(A1)
        deallocate(A2)
        deallocate(int_temp1)
        deallocate(int_temp2)

        enddo !! loop over nvl

cc remove NaN values if any
        call remove_NaNs(nvl,int1)
        call remove_NaNs(nvl,int2)
cc final integral
        select case (mom)
          case ('particle')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**4
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**5
          case ('momx')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**5
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**6
          case ('momy')
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**6
          case ('momz')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**5
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**6
          case ('energy')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**6
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**7
          case ('press1')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**5
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**6
          case ('press2')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**3
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**4
          case ('press3')
            R1(i) = (1/pi)**(1.5)*dot_product(int1,dv)/cs**4
            R2(i) = (1/pi)**(1.5)*dot_product(int2,dv)/cs**5
          case default
        end select

      enddo !! end loop over nM

      deallocate(v)
      deallocate(dv)
      deallocate(int1)
      deallocate(int2)

* ..return
      call subend ()
      return
      end subroutine CalcReflectionCoefficientsFluidTRIM

c**************************************************************************************
c*********                        Helper functions                          ***********
c**************************************************************************************

      recursive subroutine search(c,na,A,i,j)
      use b2mod_types
      use b2mod_subsys
      implicit none

** input
      integer, intent(in) :: na
      real (kind=R8), intent(in) :: c, A(na)
** output
      integer, intent(out) :: i, j
** local variables
      integer :: pivot, nnew

** documentation
cc find the coordinate points i and j such that c >= A(i) c and c < A(j),
cc using a recursive binary search

      call subini ('search')
      pivot = na+mod(na,2)
      pivot = pivot/2

      if (c.lt.A(pivot)) then
        call search(c,pivot,A(1:pivot),i,j)
      elseif (A(pivot).le.c.and.c.lt.A(pivot+1)) then
        i = pivot
        j = pivot + 1
      else
        nnew = na-pivot+1
        call search(c,nnew,A(pivot:na),i,j)
        i = pivot+i-1;
        j = pivot+j-1;
      endif

* ..return
      call subend ()
      return
      end subroutine search

      subroutine inverseVector(vec_in,vec_out)
      use b2mod_types
      use b2mod_subsys
      implicit none
      !! Subroutine to calculate the inverse of a vector
      real (kind=R8), intent(in) :: vec_in(:)
      real (kind=R8), intent(out) :: vec_out(:)
      integer :: i

      call subini ('inverseVector')

      do i = 1, size(vec_in)
        vec_out(i) = 1.0 / vec_in(i)
      enddo

* ..return
      call subend ()
      return
      end subroutine inverseVector


      subroutine meshgrid(nx, ny, x, y, xx, yy)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nx, ny
      real (kind=R8), intent(in) :: x(nx)
      real (kind=R8), intent(in) :: y(ny)
** output
      real (kind=R8), intent(out) :: xx(ny,nx)
      real (kind=R8), intent(out) :: yy(ny,nx)
** local
      integer :: i, j

cc meshgrid takes as input two vectors x and y, and produces two matrices xx and yy,
cc such that xx and yy both have dimension ny, nx, and each row of xx is equal to x,
cc and each column of yy is equal to y.

      call subini ('meshgrid')

      ! Generate the grid coordinates
      do j = 1, nx
	do i = 1, ny
          xx(i, j) = x(j)
	  yy(i, j) = y(i)
	end do
      end do

* ..return
      call subend ()
      return
      end subroutine meshgrid


      subroutine interp2(nx, ny, x, y, z, nxi, nyi, xi, yi, zi)
      use b2mod_types
      use b2mod_subsys
      implicit none
cc this is a Fortran translation of the internal Matlab function interp2
** input
      integer, intent(in) :: nx, ny, nxi, nyi
      real (kind=R8), intent(in) :: x(nx,ny)
      real (kind=R8), intent(in) :: y(nx,ny)
      real (kind=R8), intent(in) :: z(nx,ny)
      real (kind=R8), intent(in) :: xi(nxi,nyi)
      real (kind=R8), intent(in) :: yi(nxi,nyi)
** output
      real (kind=R8), intent(out) :: zi(nxi,nyi)
** local
      real (kind=R8) :: t, u
      real (kind=R8) :: dx, dy
      real (kind=R8) :: wx(2), wy(2)
      real (kind=R8) :: x0, y0, x1, y1
      real (kind=R8) :: z00, z01, z10, z11
      real (kind=R8) :: m00, m01, m10, m11
      integer :: i, j, ii, jj

cc description of variables

c   nx: number of elements in x
c   ny: number of elements in y
c   nxi: number of x coordinates of query points
c   nyi: number of y coordinates of query points
c   x: x-coordinate of input function
c   y: y-coordinate of input function
c   z: value of input function as a function of x and y
c   xi: x-coordinates of query points
c   yi: y-coordinates of query points

c   zi: values of query points

c   t, u: temporary storage of x(i,j) and y(i,j)
c   dx: dimensionless number to indicate where the query point is between the two surrounding x coordinates
c   dy: dimensionless number to indicate where the query point is between the two surrounding y coordinates
c   wx: weights for x direction
c   wy: weights for y direction
c   x0: x-point 'left' of query point
c   y0: y-point 'top' of query point
c   x1: x-point 'right' of query point
c   y1: y-point 'bottom' of query point
c   z00: value of 'top left' point
c   z01: value of 'top right' point
c   z10: value of 'bottom left' point
c   z11: value of 'bottom right' point
c   m00: contribution from 'top left' point
c   m01: contribution from 'top right' point
c   m10: contribution from 'bottom left' point
c   m11: contribution from 'bottom right' point

c   i, j: loop variables over query points
c   ii, jj: variables to find the x and y values surrounding a query point

      call subini ('interp2')

      ! Loop over output/query points
      do j = 1, nyi
        do i = 1, nxi
          t = xi(i, j)
          u = yi(i, j)

          ! Find the appropriate indices for interpolation
          if (t < x(1, 1) .or. t > x(1, ny) .or. u < y(1, 1)
     &                                   .or. u > y(nx, 1)) then
            zi(i, j) = 0.0  ! Point outside the grid, assign zero
          else
            ii = 0
            jj = 0
            do while (t >= x(1, 1+ii))
              ii = ii + 1
            end do
            do while (u >= y(1+jj, 1))
              jj = jj + 1
            end do

        ! Get the four surrounding points for interpolation
            x0 = x(1, ii)
            y0 = y(jj, 1)
            x1 = x(1, ii+1)
            y1 = y(jj+1, 1)
            z00 = z(jj, ii)
            z01 = z(jj, ii+1)
            z10 = z(jj+1, ii)
            z11 = z(jj+1, ii+1)

        ! Compute interpolation weights
            dx = (t - x0) / (x1 - x0)
            dy = (u - y0) / (y1 - y0)
        ! linear
            wx = [1.0 - dx, dx]
            wy = [1.0 - dy, dy]


        ! Perform interpolation
            m00 = z00 * wx(1) * wy(1)
            m01 = z01 * wx(2) * wy(1)
            m10 = z10 * wx(1) * wy(2)
            m11 = z11 * wx(2) * wy(2)
            zi(i, j) = m00 + m01 + m10 + m11
          end if
        end do
      end do

* ..return
      call subend ()
      return

      end subroutine interp2


      subroutine extrapolate_RN_minenergy(nv, ntheta, RNv, Ei_v,
     &                                nE, nA, E, A, RN, Thetai_v)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nv, ntheta, nE, nA
      real (kind=R8), intent(in) :: Ei_v(nv)
      real (kind=R8), intent(in) :: E(nE)
      real (kind=R8), intent(in) :: A(nA)
      real (kind=R8), intent(in) :: RN(nE,nA)
      real (kind=R8), intent(in) :: Thetai_v(ntheta)
** input/output
      real (kind=R8), intent(inout) :: RNv(nv,ntheta)
** local
      real (kind=R8) :: extrap_value1, RN_A(nA)
      integer :: i, j

      call subini ('extrapolate_RN_minenergy')

      ! Interpolate Rpv1 for Ei_v < E(1)
      RN_A = RN(1,:)
      do i = 1, nv
        if (Ei_v(i).lt.E(1)) then
          do j = 1,ntheta
            extrap_value1 = interp1(nA, A, RN_A,
     &                            Thetai_v(j))
            RNv(i,j) = extrap_value1*Ei_v(i)/E(1)
          enddo
        endif
      enddo

* ..return
      call subend ()
      return
      end subroutine extrapolate_RN_minenergy

      subroutine extrapolate_Rpv_angle(nv,ntheta,Rpv, Thetai_v,
     &                                      nA, nE, A, E, Rp, Ei_v)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer , intent(in) :: nv, ntheta, nA, nE

      real (kind=R8), intent(in) :: Thetai_v(ntheta)
      real (kind=R8), intent(in) :: A(nA)
      real (kind=R8), intent(in) :: E(nE)
      real (kind=R8), intent(in) :: Rp(nE,nA)
      real (kind=R8), intent(in) :: Ei_v(nv)
** input/output
      real (kind=R8), intent(inout) :: Rpv(nv,ntheta)
** local
      real (kind=R8) :: extrap_value1, extrap_value2
      real (kind=R8) :: slope, delta_theta
      integer :: i , j

cc description of variables
c    nv: number of elments in v
c    ntheta: number of elements in theta
c    nA: number of elements in A
c    nE: number of elements in E

c    Thetai_v: discretization vector of the incident polar angle corresponding to theta
c    A: discretization vector of the incident polar angle in TRIM
c    E: discretization vector of the incident energy in TRIM
c    Rp: pressure coefficient as a function of E and A
c    Ei_v: energy corresponding to v

c    Rpv: pressure coefficient as a function of v and theta

c    extrap_value1: value of Rp for largest angle in A, interpolated to correct energy
c    extrap_value2: value of Rp for second-to-largest angle in A, interpolated to correct enegy
c    slope: slope between extrap_value1 and extrap_value2
c    delta_theta: distance over which we need to extrapolate
c    i, j: loop variables

      call subini ('extrapolate_Rpv_angle')
      ! Interpolate Rpv for Thetai_v > A(end)
      do i = 1, nv
        do j = 1, ntheta
          if (Thetai_v(j).gt.A(size(A))) then
           ! Extrapolate using the last two columns of Rp1 and Ei_v
            extrap_value1 = interp1(nE, E, Rp(:, nA),
     &                            Ei_v(i))
            extrap_value2 = interp1(nE, E, Rp(:, nA-1),
     &                            Ei_v(i))

            ! Calculate slope and delta_theta
            slope = (extrap_value1 - extrap_value2) /
     &                                     (A(size(A)) - A(size(A)-1))
            delta_theta = Thetai_v(j) - A(size(A))

            ! Interpolate Rpv using the last column of Rp1 and Ei_v
            Rpv(i,j) = slope * delta_theta +
     &      interp1(nE, E, Rp(:, size(Rp, 2)),Ei_v(i))
          endif
        enddo
      enddo


* ..return
      call subend ()
      return
      end subroutine extrapolate_Rpv_angle


      subroutine extrapolate_Rpv_minenergy(nv, ntheta, Rpv, Ei_v,
     &                                nE, nA, E, A, Rp, Thetai_v)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nv, ntheta, nE, nA

      real (kind=R8), intent(in) :: Ei_v(nv)
      real (kind=R8), intent(in) :: E(nE)
      real (kind=R8), intent(in) :: A(nA)
      real (kind=R8), intent(in) :: Rp(nE,nA)
      real (kind=R8), intent(in) :: Thetai_v(ntheta)
** input/output
      real (kind=R8), intent(inout) :: Rpv(nv,ntheta)
** local
      real (kind=R8) :: extrap_value1, extrap_value2
      real (kind=R8) :: slope, delta_Ei, Rp_a(nA), Rp_b(nA)
      integer :: i, j

cc description of variables
c    nv: number of elements in v
c    ntheta: number of elements in theta
c    nE: number of elements in E
c    nA: number of elments in A

c    Ei_v: energy corresponding to v
c    E: discretization vector of the incident energy in TRIM
c    A: discretization vector of the incident polar angle in TRIM
c    Rp: pressure coefficient as a function of E and A
c    Thetai_v: discretization vector of the incident polar angle corresponding to theta

c    Rpv: pressure coefficient as a function of v and theta

c    extrap_value1: value at lowest energy, interpolated to the correct angle
c    extrap_value2: value at second-to-lowest energy, interpolated to the correct angle
c    slope: slope between extrap_value1 and extrap_value2
c    delta_Ei: distance over which we need to extrapolate
c    Rp_a: array of Rp values for the lowest energy
c    Rp_b: array of Rp values for the second-to-lowest energy
c    i, j: loop variables

      call subini ('extrapolate_Rpv_minenergy')

      ! Interpolate Rpv1 for Ei_v < E(1)
      Rp_a = Rp(1, :)
      Rp_b = Rp(2, :)
      do i = 1, nv
        if (Ei_v(i).lt.E(1)) then
          do j = 1,ntheta
            extrap_value1 = interp1(nA, A, Rp_a,
     &                            Thetai_v(j))
            extrap_value2 = interp1(nA, A, Rp_b,
     &                            Thetai_v(j))

            slope = (extrap_value1 - extrap_value2) / (E(1) - E(2))
            delta_Ei = Ei_v(i) - E(1)

            Rpv(i,j) = slope * delta_Ei + interp1(nA, A,
     &      Rp_a, Thetai_v(j))
          enddo
        endif
      enddo

* ..return
      call subend ()
      return
      end subroutine extrapolate_Rpv_minenergy

      subroutine extrapolate_Rpv_maxenergy(nv, ntheta, Rpv, Ei_v,
     &                             nE, nA, E, A, Rp, Thetai_v)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nv, ntheta, nE, nA

      real (kind=R8), intent(in) :: Ei_v(nv)
      real (kind=R8), intent(in) :: E(nE)
      real (kind=R8), intent(in) :: A(nA)
      real (kind=R8), intent(in) :: Rp(nE,nA)
      real (kind=R8), intent(in) :: Thetai_v(ntheta)
** input/output
      real (kind=R8), intent(inout) :: Rpv(nv,ntheta)
** local
      real (kind=R8) :: extrap_value1, extrap_value2
      real (kind=R8) :: slope, delta_Ei
      integer :: i,j
      real (kind=R8) :: Rp_a(nA), Rp_b(nA)

cc similar as subroutine extrapolate_Rpv_minenergy, but for maximum

      call subini ('extrapolate_Rpv_maxenergy')

      Rp_a = Rp(size(Rp, 1), :)
      Rp_b = Rp(size(Rp, 1)-1, :)
      ! Interpolate Rpv for Ei_v > E(end)
      do i = 1,nv
        do j = 1,ntheta
          if(Ei_v(i).gt.E(size(E))) then
            extrap_value1 = interp1(nA, A, Rp_a, Thetai_v(j))
            extrap_value2 = interp1(nA, A, Rp_b, Thetai_v(j))

            slope = (extrap_value1 - extrap_value2) /
     &                         (E(size(E)) - E(size(E) - 1))
            delta_Ei = Ei_v(i) - E(size(E))

            Rpv(i,j) = slope * delta_Ei +
     &        interp1(nA, A, Rp_a, Thetai_v(j))
          endif
        enddo
      enddo

* ..return
      call subend ()
      return
      end subroutine extrapolate_Rpv_maxenergy

      function interp1(nx, x, y, xq) result(yq)
      use b2mod_types
      use b2mod_subsys
      implicit none
** input
      integer, intent(in) :: nx
      real (kind=R8), intent(in) :: x(nx)  ! Independent variable
      real (kind=R8), intent(in) :: y(nx) ! Dependent variable
      real (kind=R8), intent(in) :: xq ! Query point
** output
      real (kind=R8) :: yq ! Interpolated values
** local
      real (kind=R8) :: xp, yp
      integer ::  j

cc description of variables

c  nx: number of points in x
c  x: coordinate
c  y: function value as a function of the x-coordinate
c  xq: query point along x

c  yq: interpolated value corresponding to query point xq

c  xp: last point along x which is smaller than the query point xq
c  yp: function value corresponding to xp

c  j: loop variable

** documentation
c simple 1D linear interpolation / extrapolation

      ! Check if query point is outside the range of x
        if (xq < x(1) .or. xq > x(size(x))) then
        ! Extrapolate using the nearest available points
          if (xq < x(1)) then
            yq = y(1) - (y(2)-y(1))*(x(2)-x(1))*(x(1)-xq)
          else !! true extrapolation
            yq = y(nx) + (y(nx)-y(nx-1))/(x(nx)-x(nx-1))*(xq-x(nx))
          end if
        else
      ! Find the indices of the nearest points in x
          do j = 2, size(x)
            if (xq <= x(j)) then
              xp = x(j-1)
              yp = y(j-1)
              exit
            end if
          end do
          ! Perform linear interpolation
          yq = yp + (y(j) - yp) * (xq - xp) / (x(j) - xp)
        end if

      return
      end function interp1

      subroutine remove_NaNs(nx,x)
      use b2mod_types
      use b2mod_subsys
#ifndef NO_IEEE
      use ieee_arithmetic    ! IGNORE
#endif
      implicit none
** input/output
      integer :: nx
      real (kind=R8), intent(inout) :: x(nx)
** local
      integer :: i
#ifdef NO_IEEE
      logical, external :: my_isnan
#endif

      call subini ('remove_NaNs')

** documentation
c replace NaN values by 0 (otherwise a single NaN value in x could cause SUM(x) to be NaN)

      do i = 1,nx
#ifndef NO_IEEE
        if (ieee_is_nan(x(i))) then
#else
        if (my_isnan(x(i))) then
#endif
          x(i) = 0.0_R8
        endif
      enddo

* ..return
      call subend ()
      return
      end subroutine remove_NaNs

c**************************************************************************************
c*************************   Read & write routines     ********************************
c**************************************************************************************

      subroutine read_TRIM(incident,material)
      use b2mod_types
      use b2mod_xerset
      use b2mod_subsys
      implicit none
** input
      character(len=*), intent(in) :: incident, material
** local
      character*256 :: filename
      logical :: open_file
      integer :: i,j,k,l
      character*80 header

      external open_file

cc description of variables
c   incident: name of incident particle (at present: H, D or T)
c   material: chemical symbol of wall material (e.g. C, Be, Cu, Fe, Mo, W, ...)

c   filename: name of the TRIM file
c   open_file: opens the file, and the logical open_file will be TRUE if opening the file was successful.
c   i, j, k, l: loop variables
c   header: file header

      call subini ('read_TRIM')
cc set correct name
      filename = incident(1:len_trim(incident))//'_on_'//
     &           material(1:len_trim(material))

cc the format of the TRIM data files is hard-coded because it is not expected to change (until some new database is created)

cc open the TRIM file stored in the EIRENE source code
      if (open_file(9,filename,
     &      'modules/Eirene/Database/Surfacedata/TRIM|' //
     &      'modules/B2.5/Database/trim.data/refl.data.local')) then
      else
        write(*,*) 'If you do not have access to EIRENE,',
     &      'request the required TRIM file from an IO member ',
     &      'and place it in',
     &      'modules/B2.5/Database/trim.data/refl.data.local'
        call xerrab ('TRIM file not found')
      endif

      TRIMdata = 0.0_R8

cc read header
      read (9, '(A)') header
      write(*,*) 'Start reading TRIM file'
      write(*,*) header

cc read data blocks
      k = 0
      do i = 1,84 !! loop over the 84 blocks with different combinations of incident energy and incident angles
        k = k + 1
        !! read first single line (incident energy, incident angle, reflection coefficient)
        read(9,*,end=101) (TRIMdata(k,j), j=1,7)
        k = k + 1
        !! read second single line (conditional distribution reflected energy)
        read(9,*,end=101) (TRIMdata(k,j), j=1,5)

        !! read block of 5 (conditional distribution of reflected polar angle for a given reflected energy)
        do l=1,5
          k = k + 1
          read(9,*,end=101) (TRIMdata(k,j), j=1,5)
        enddo

        !! read block of 25 (conditional distribution of reflected azimuthal angle for a given reflected energy and reflected polar angle)
        do l=1,25
          k = k + 1
          read(9,*,end=101) (TRIMdata(k,j), j=1,5)
        enddo

      enddo

cc close the file
101   close (9)
      write (*,*) 'Finished reading TRIM file'

* ..return
      call subend ()
      return
      end subroutine read_TRIM


      subroutine write_fluid_TRIM(incident,material)
      use b2mod_types
      use b2mod_xerset
      use b2mod_b2cmpa
      use b2mod_elements
      use b2mod_subsys
      use b2mod_version
      implicit none
** input
      character(len=*), intent(in) :: incident, material
** local
      character*256 :: filename
      integer :: nout, idum(3), io, ntot
      character*80 version
      character*120 :: label
      external cfwuin, cfwuch, cfwure
#ifdef LEGACYCOMP
      integer newunit
      external newunit
#endif

cc description of variables
c   incident: incident particle
c   material: wall material
c   filename: name of output file (e.g. hc.fnn.coarse)
c   nout: dummy numbers for the output files
c   idum: dummy variable to group output
c   io: error number
c   ntot: total number of values per field
c   version: TODO
c   label: description of the file

      call subini ('write_fluid_TRIM')
cc set correct name
      filename = incident(1:len_trim(incident))//
     &           material(1:len_trim(material))//
     &           '.fnn.coarse'
      call lower_case(filename)

cc open new file
#ifndef LEGACYCOMP
      open (newunit=nout,file=filename,iostat=io)
#else
      nout = newunit()
      open (unit=nout,file=filename,iostat=io)
#endif

      version = 'VERSION'//default_version//' Fortran'

      call cfverw (nout, version)
      idum(1) = nM
      idum(2) = nT
      idum(3) = nSh
      call cfwuin (nout, 3, idum, 'nm,nt,nshpottau')

      label = 'Coefficients for advanced fluid boundary conditions '//
     &        'based on TRIM data'
      call cfwuch (nout, 120, label, 'label')

      call cfwure (nout, nM, M, 'M')
      call cfwure (nout, nT, T, 'T')
      call cfwure (nout, nSh, Sh, 'shpottau')

      ntot = nM*nT*nSh
      call cfwure (nout, ntot, R1, 'R1')
      call cfwure (nout, ntot, R2, 'R2')
      call cfwure (nout, ntot, Rmx1, 'Rmx1')
      call cfwure (nout, ntot, Rmx2, 'Rmx2')
      call cfwure (nout, ntot, Rmy1, 'Rmy1')
      call cfwure (nout, ntot, Rmy2, 'Rmy2')
      call cfwure (nout, ntot, Rmz1, 'Rmz1')
      call cfwure (nout, ntot, Rmz2, 'Rmz2')
      call cfwure (nout, ntot, RE1, 'RE1')
      call cfwure (nout, ntot, RE2, 'RE2')
      call cfwure (nout, ntot, Rp11, 'Rp11')
      call cfwure (nout, ntot, Rp12, 'Rp12')
      call cfwure (nout, ntot, Rp21, 'Rp21')
      call cfwure (nout, ntot, Rp22, 'Rp22')
      call cfwure (nout, ntot, Rp31, 'Rp31')
      call cfwure (nout, ntot, Rp32, 'Rp32')
cc close the file
      close (nout)

* ..return
      call subend ()
      return
      end subroutine write_fluid_TRIM

      subroutine deallocate_all_TRIM()
      use b2mod_subsys
      implicit none

      call subini ('deallocate_all_TRIM')

** documentation
cc deallocate all allocatable global variables

      deallocate(M)
      deallocate(T)
      deallocate(Sh)
      deallocate(R1)
      deallocate(R2)
      deallocate(Rmx1)
      deallocate(Rmx2)
      deallocate(Rmy1)
      deallocate(Rmy2)
      deallocate(Rmz1)
      deallocate(Rmz2)
      deallocate(RE1)
      deallocate(RE2)
      deallocate(Rp11)
      deallocate(Rp12)
      deallocate(Rp21)
      deallocate(Rp22)
      deallocate(Rp31)
      deallocate(Rp32)

      deallocate(Rv)
      deallocate(Rmomparv)
      deallocate(Rmomperpv)
      deallocate(Rpv1)
      deallocate(Rpv2)
      deallocate(Rpv3)
      deallocate(REv)

      call subend ()
      return
      end subroutine deallocate_all_TRIM

      end module b2mod_create_afn_bc_tables
