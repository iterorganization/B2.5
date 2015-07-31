C.....................................................................
c Program to solve the heat diffusion equation in 2 dimensions

      subroutine prehd2D(ns,isize,list_x,list_y,dt,isign,
     1  ixf,iyf,idir,area,BoRiS)
      use b2mod_types
      use b2mod_wall
      use b2mod_rates
      use b2mod_plasma
      use b2mod_constants
      use b2mod_external
      implicit none

      integer ns,isize,isign,idir
      integer list_x(isize),list_y(isize)
      real (kind=R8) :: dt,BoRiS
      integer ixf(isize),iyf(isize)
      real (kind=R8) :: area(isize)

      integer iwall, is
      integer bndx(isize)
      real (kind=R8) :: qin, qel
      real (kind=R8) :: Q(isize), Tp(isize,ndepth_nml)

      external subini, subend, xertst
      intrinsic max

      call subini ('prehd2D')
      call xertst (isize.ge.1, 'faulty input isize')

      do iwall = 1, isize
        bndx(iwall) = xymap(list_x(iwall),list_y(iwall))
        qin = 
     .   max(0.0_R8,isign*(fhi(ixf(iwall),iyf(iwall),idir)+
     .                     fhi_ext(ixf(iwall),iyf(iwall),idir)))/
     .   area(iwall)
        do is=0,ns-1
         qin=qin+max(0.0_R8,isign*fna(ixf(iwall),iyf(iwall),idir,is))
     &        *(kinrgy(list_x(iwall),list_y(iwall),is)*(1.0_R8-BoRiS) + 
     &          rpt(ixf(iwall),iyf(iwall),is)*ev)
     &        /area(iwall)
        enddo
        do is=0,ns_ext-1
         qin=qin+max(0.0_R8,isign*fa_ext(ixf(iwall),iyf(iwall),idir,is))
     &    *(0.5_R8*am_ext(is)*mp*
     &      ua_ext(list_x(iwall),list_y(iwall),is)**2*(1.0_R8-BoRiS) + 
     &      pt_ext(ixf(iwall),iyf(iwall),is)*ev)
     &    /area(iwall)
        enddo
        qel =
     .   max(0.0_R8,isign*fhe(ixf(iwall),iyf(iwall),idir))/area(iwall)
        Q(iwall)=qin+qel
        Tp(iwall,1:ndepth_nml) = target_temp(bndx(iwall),1:ndepth_nml)    
      enddo

      call heatdiff2D(isize,ns,bndx,Q,dt,Tp)

      do iwall = 1, isize
        target_temp(bndx(iwall),1:ndepth_nml) = Tp(iwall,1:ndepth_nml)
      enddo

      call subend ()
      return
      end subroutine prehd2D

c.....................................................................

c This is the subroutine which solves for a heat diffusion equation
c of the form:
c              doe T         doe^2 T     doe D   doe T
c             ------- =  D  --------- + ------- ------- + sources/sinks
c              doe t         doe^2 x     doe x   doe x
c T is the temperature, t the time variable, x the space variable and
c D the heat diffusion coefficient.

c i.e., ....
c       DU = (1/rho/Cp) * (K*UXX + gradK*UX + S)
c       DU is the rate of change of temperature
c       UXX is "doe^2 T / doe^2 x"   and UX is "doe T / doe x"
c with Neumann bondary condition at the plasma facing surface side
c i.e.  UX is specified.
c and a Dirichlet type boundary condition at the coolant side.
c i.e. T is specified.

c Inputs:
c  Q = incident total heat flux (Watts/m^2/sec).
c  Tcbc = Temperature at the cooling end (K).
c  dt = Time step of simulation (sec)

c Outputs:
c  Tp = Target surface temperature (K).
      subroutine heatdiff2D(nw,ns,bndx,Q,dt,Tp)
      use b2mod_types
      use b2mod_wall
      implicit none

      integer nw, ns, n2
      integer :: bndx(nw)
      real (kind=R8) :: Q(nw), dt, Tp(nw, ndepth_nml)

      real (kind=R8) :: Qtot, Tcbc
      common /phyparams/ Qtot,Tcbc

      integer :: i, iwall, iters, itmax
      parameter (itmax=100)
      real (kind=R8) :: U(nw, ndepth_nml), rhs(nw, ndepth_nml)
      real (kind=R8), allocatable, save ::
     1 aaf(:),wrk0(:),wrk1(:),wrk2(:)
      integer, allocatable, save :: iwrk(:)
      real (kind=R8) ::
     *   aa5(nw,ndepth_nml,-2:2), nrmc(0:itmax-1)
      real (kind=R8) :: t0, ta, tb, tc, x(ndepth_nml), y(nw,ndepth_nml), 
     *   rho

      real (kind=R8) :: 
     *   Cp, Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT, 
     *   qevap, qrad, mass_density
      external Cp, Kappa, dKappabydT, Kappa_perp, dKappa_perpbydT,
     *   qevap, qrad, mass_density
      external subini, subend
      external slv5pt, fic

      call subini('heatdf2D')

c Checking if it is the first call and if it is, call fic(x,u), else
c use previous Tp profile as initial guess.

      if (dt.le.0.0_R8) then
c       Calling subroutine fic(ibnd,x,u) which supplies the
c       initial condition u(x) i.e.. the temperature "u" as a
c       function of the position "x" at ndepth_nml grid points.
        do iwall = 1, nw
          Qtot = Q(iwall)
          Tcbc = backplate_temp(bndx(iwall))
          call fic(bndx(iwall),x,U(iwall,1:ndepth_nml))
          do i=1,ndepth_nml
            Tp(iwall,i) = U(iwall,i)
          enddo
          call subend ()
          return
        enddo
      else
        do iwall = 1, nw
          do i = 1, ndepth_nml
            U(iwall,i) = Tp(iwall,i)
          end do
        enddo
      endif

*   ..set n2
      n2 = nw*ndepth_nml
*    ..obtain workspace for slv5pt
      allocate(aaf(1:n2*(3*min(nw,ndepth_nml)+1)))
      allocate(wrk0(1:n2))
      allocate(wrk1(1:1))
      allocate(wrk2(1:1))
      allocate(iwrk(1:n2))

*   ..initialise the aa5 matrix and rhs vector
      aa5(1:nw,1:ndepth_nml,-2:2) = 0.0_R8
      rhs(1:nw,1:ndepth_nml) = 0.0_R8

*   ..build depth neighbours (i,j-1) and (i,j+1) terms
      do iwall = 1, nw
        x(1:ndepth_nml) = target_depth(bndx(iwall),1:ndepth_nml)
*    ..treat incoming heat flux b.c.
        i = 1
        Qtot = Q(iwall)-qevap(bndx(iwall),ns,Tp(iwall,i))
     &                 - qrad(bndx(iwall),Tp(iwall,i),.true.)
        t0 = Kappa(bndx(iwall),i,Tp(iwall,i)) / (x(2)-x(1))
        aa5(iwall,i,-2) = 0.0_R8
        aa5(iwall,i, 0) = aa5(iwall,i, 0) + t0
        aa5(iwall,i, 2) = aa5(iwall,i, 2) - t0
        rhs(iwall,i)    = rhs(iwall,i) + Qtot
*    ..main terms
        do i = 2, ndepth_nml-1
          ta = dKappabydT(bndx(iwall),i,Tp(iwall,i)) * 
     1     (Tp(iwall,i+1)-Tp(iwall,i-1)) / 
     2     (x(i+1)-x(i-1))**2
          tb = 2.0_R8 / (x(i+1)-x(i)) / (x(i+1)-x(i-1))
          tc = 2.0_R8 / (x(i)-x(i-1)) / (x(i+1)-x(i-1))
          aa5(iwall,i,-2) = aa5(iwall,i,-2) + ta      - tc
          aa5(iwall,i, 0) = aa5(iwall,i, 0)      + tb + tc
          aa5(iwall,i, 2) = aa5(iwall,i, 2) - ta - tb
        enddo
*    ..treat backplate cooling b.c.
        i = ndepth_nml
        if (inertial_cooling(iwall)) then
          t0 = Kappa(bndx(iwall),i,Tp(iwall,i)) / (x(i)-x(i-1))
          aa5(iwall,i,-2) = aa5(iwall,i,-2) - t0
          aa5(iwall,i, 0) = aa5(iwall,i, 0) + t0
          aa5(iwall,i, 2) = 0.0_R8
!! FIXME : Could also include radiative and evaporative cooling term here !!
          rhs(iwall,i) = rhs(iwall,i) - 
     .     100.0_R8*(Tp(iwall,i)-backplate_temp(bndx(iwall))) ! we assume the contact heat conductivity of 100 W.m-2.K-1
        else
          aa5(iwall,i,-2) = 0.0_R8
          aa5(iwall,i, 0) = 1.0_R8
          aa5(iwall,i, 2) = 0.0_R8
          rhs(iwall,i) = backplate_temp(bndx(iwall))
        endif
      enddo

*   ..build side neighbours (i-1,j) and (i+1,j) terms
      do i = 1, ndepth_nml
        y(1:nw,i) = target_width(bndx(1:nw),i)
      enddo
*    ..treat radiative cooling b.c.
      iwall = 1
      do i = 1, ndepth_nml-1
!! FIXME: Careful about the area over which the radiative cooloing occurs, it should be here the side area, not the top area !!
!! Can also include evaporative cooling term !!
        Qtot = -qrad(bndx(iwall),Tp(iwall,i),.false.)
        t0 = Kappa_perp(bndx(iwall),i,Tp(iwall,i)) /
     2   (0.5_R8*(y(iwall,i)+y(iwall+1,i))) 
        aa5(iwall,i,-1) = 0.0_R8
        aa5(iwall,i, 0) = aa5(iwall,i, 0) + t0
        aa5(iwall,i, 1) = aa5(iwall,i, 1) - t0
        rhs(iwall,i)    = rhs(iwall,i) + Qtot
      enddo
*    ..main terms
      do iwall = 2, nw-1
        do i = 1, ndepth_nml-1
          if (y(iwall,i).gt.0.0_R8) then
            ta = dKappa_perpbydT(bndx(iwall),i,Tp(iwall,i)) * 
     1       (Tp(iwall+1,i)-Tp(iwall-1,i)) / 
     2       (y(iwall,i)+0.5_R8*(y(iwall-1,i)+y(iwall+1,i)))**2
            tb = 2.0_R8 / y(iwall,i) / (y(iwall+1,i)+y(iwall,i)) 
            tc = 2.0_R8 / y(iwall,i) / (y(iwall,i)+y(iwall-1,i))
          endif
          aa5(iwall,i,-1) = aa5(iwall,i,-1) + ta      - tc
          aa5(iwall,i, 0) = aa5(iwall,i, 0)      + tb + tc
          aa5(iwall,i, 1) = aa5(iwall,i, 1) - ta - tb
        enddo
      enddo
*    ..treat radiative cooling b.c.
      iwall = nw
      do i = 1, ndepth_nml-1
!! FIXME: Careful about the area over which the radiative cooloing occurs, it should be here the side area, not the top area !!
!! Can also include evaporative cooling term !!
        Qtot = -qrad(bndx(iwall),Tp(iwall,i),.false.)
        t0 = Kappa_perp(bndx(iwall),i,Tp(iwall,i)) /
     2   (0.5_R8*(y(iwall,i)+y(iwall-1,i)))
        aa5(iwall,i,-1) = aa5(iwall,i,-1) - t0
        aa5(iwall,i, 0) = aa5(iwall,i, 0) + t0
        aa5(iwall,i, 1) = 0.0_R8
        rhs(iwall,i)    = rhs(iwall,i) + Qtot
      enddo

*   ..add the d/dt term contribution
      do iwall = 1, nw
        do i = 1, ndepth_nml-1
          if (target_vol(bndx(iwall),i).gt.0.0_R8) then
            rho = mass_density(bndx(iwall),i) *
!! FIXME : Need to test if proper quantity here is
     .            (x(i+1)-x(i))
! or may be :
!    .             target_vol(bndx(iwall),i)
            t0 = rho * Cp(bndx(iwall),i,Tp(iwall,i)) / dt
            aa5(iwall, i, 0) = aa5(iwall,i, 0) + t0
            rhs(iwall, i) = rhs(iwall,i) + t0 * Tp(iwall,i)
          endif
        enddo
      enddo

*   ..correct equations for special case of zero or negative width
      do iwall = 2, nw-1
        do i = 1, ndepth_nml
          if (target_width(bndx(iwall),i).le.0.0_R8) then
            aa5(iwall, i, -2) = 0.0_R8
            aa5(iwall, i, -1) =-0.5_R8
            aa5(iwall, i,  0) = 1.0_R8
            aa5(iwall, i,  1) =-0.5_R8
            aa5(iwall, i,  2) = 0.0_R8
            rhs(iwall, i)     = 0.0_R8
          endif
        enddo
      enddo

*   ..solve the temperature equation
*    ..using slv5pt
      iters = 1
      call slv5pt (nw, ndepth_nml, aa5, rhs, iters, U, nrmc,
     &   min(nw,ndepth_nml), aaf, wrk0, iwrk)

*   ..copy over new solution
      do iwall = 1, nw
        do i = 1, ndepth_nml
          Tp(iwall,i) = min( U(iwall,i), 
     .                       matsurf_vapor_temp(bndx(iwall)))
        enddo
      enddo

      deallocate(aaf)
      deallocate(wrk0)
      deallocate(wrk1)
      deallocate(wrk2)
      deallocate(iwrk)
      call subend ()
      return
      end subroutine heatdiff2D
