C.....................................................................
c Program to solve the heat diffusion equation

c  Last modified by Manoj Warrier on (08-02-2004)

c Using the numerical method written by M. Kundu to
c SOLVE PARABOLIC EQUATIONS USING THE METHOD OF LINES

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
      subroutine heatdiff1D(ibnd,ns0,Q,dt,Tp,output)
      use b2mod_types
      use b2mod_wall
      implicit none
      external lines2c

      integer i
      integer ibnd, ns0, iwall, output, ns

      logical active
      real (kind=R8) :: Q, dt, Qtot, Tcbc
      real (kind=R8) :: Tp(ndepth_nml)
      real (kind=R8) :: u(7*ndepth_nml), u1(7*ndepth_nml)
      real (kind=R8) :: x(ndepth)
      real (kind=R8) :: x0, xn, t0, tn, dx(ndepth-1)
      real (kind=R8) :: Cp, Kappa, mass_density
      external Cp, Kappa, mass_density

#ifdef NO_NAG
      integer ier, iflg, iopt, neq, nmax, nstep
      real (kind=R8) :: hstep, rerr
      real (kind=R8) :: du(7*ndepth_nml), yf(ndepth_nml), ut(ndepth_nml)
      real (kind=R8) :: wk(3*ndepth_nml + ndepth_nml*2*ndepth_nml)
      integer iwk(ndepth_nml)
      external mstep
#else
      external pdedef, bndary, D03PCF
#endif
      integer ncall
      data ncall/0/
      common /zzline/ dx,x0,xn,x
      common /phyparams/ Qtot,Tcbc,active
      common /location/ iwall, ns

      integer, save :: npde, nw, niw
      integer :: itask, itrace, ind, ifail
      integer, allocatable, save :: iw(:)
      real (kind=R8), allocatable, save :: U_nag(:,:), w(:)
      real (kind=R8) :: acc

      external xertst, subini, subend

      call subini('heatdiff1D')
      if (ncall.eq.0) then
        ncall = ncall + 1
        npde = 1
        nw = (10+6*npde)*npde*ndepth_nml+
     .       (21+3*npde)*npde+7*ndepth_nml+56
        niw = npde*ndepth_nml+25
        allocate(U_nag(npde,ndepth_nml))
        allocate(iw(niw))
        allocate(w(nw))
      endif

      if (output.ge.10) then
        write(*,*)'Wall element number =',ibnd
        write(*,*)'Total energy flux =',Q,'Joules/m^2/s'
        write(*,*)'Temp at the cold end =',backplate_temp(ibnd),'K'
        write(*,*)'Temp at the warm end =',Tp(1),'K'
        write(*,*)'Thickness of the plate =',plate_thickness(ibnd),'m'
        if (inertial_cooling(ibnd)) then
          write(*,*)'Plate element is inertially cooled'
        else
          write(*,*)'Plate element is actively cooled'
        endif
        write(*,*)'Time step of simulation =',dt,'sec'
      endif

      Qtot = Q
      Tcbc = backplate_temp(ibnd)
      active = .not.inertial_cooling(ibnd)
      t0 = 0.0_R8
      tn = dt
      iwall = ibnd
      ns = ns0

c  Compute positions and intervals
      x(1:ndepth_nml) = target_depth(ibnd,1:ndepth_nml)
      x0 = target_depth(ibnd,1)
      xn = target_depth(ibnd,ndepth_nml)
      do i=1,ndepth_nml-1
        dx(i) = x(i+1)-x(i)
      enddo

c Checking if it is the first call and if it is, call fic(x,u), else
c use previous Tp profile as initial guess.
      if (dt.le.0.0_R8) then
c       Calling subroutine fic(ibnd,x,u) which supplies the
c       initial condition u(x) i.e.. the temperature "u" as a
c       function of the position "x" at ndepth_nml grid points.
        call fic(ibnd,x,u)
        do i=1,ndepth_nml
          Tp(i) = u(i)
        enddo
        call subend ()
        return
      else
        do i=1,ndepth_nml-1
          u1(i)= Tp(i)
        enddo
      endif

c No. of equations to be solved.
#ifdef NO_NAG
      neq = ndepth_nml - 1   ! since we specify T at one end and grad(T) at the other
      iflg = 0
      iopt = 1
      hstep = dt
      nmax = 100000000
      rerr = 1.0e-8_R8
      call mstep(ibnd,ns,neq,u1,ut,du,lines2c,hstep,t0,tn,yf,
     &  rerr,nstep,nmax,ier,iflg,iopt,wk,iwk)
      call xertst(ier.eq.0,'ier.ne.0 in heatdiff1D')
      if (output.ge.10) write(6,*) 'No. of successful steps = ',nstep
      Tp(ndepth_nml) = Tcbc
      do i=1,ndepth_nml-1
        Tp(i) = min( yf(i), matsurf_vapor_temp(ibnd))
      enddo
#else
cxpb   New NAG solution to replace MSTEP calls
      ind = 0
      ifail = 1
      itask = 1
      itrace = -1
      do i = 1, ndepth_nml
        U_nag(1,i)=Tp(i)
      enddo
      acc = 1.0e-8_R8

      call D03PCF(npde, 0, t0, tn, pdedef,
     &   bndary, U_nag, ndepth_nml, x, acc, w, nw, iw, niw,
     &   itask, itrace, ind, ifail)
      if(ifail.ne.0.and.ifail.ne.10) write(*,*) 'D03PCF:ifail',ifail
      call xertst (ifail.eq.0.or.ifail.eq.10,
     & 'Abnormal exit from D03PCF')
      do i=1,ndepth_nml
        Tp(i) = min( U_nag(1,i), matsurf_vapor_temp(ibnd))
      enddo
#endif

      call subend ()
      return
      end subroutine heatdiff1D
c_______________________________________________________________________

c               SUBROUTINE FIC

c Subroutine to specify the initial condition: Temperature "uu" as a
c function of the spatial coordinate "xx".

c Both xx and uu are outputs from the subroutine.

c Thermal conductivity is assumed to vary as 1/(KA*T+KB) with
c temperature T (= uu in this routine).
c The above dependence follows the reference
      subroutine fic(ibnd,xx,uu)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none

      integer i, ibnd
      logical active
      real (kind=R8) :: Qtot, Tcbc, KAb, KBb, Tc, KAc, KBc, Tl, KAl, KBl
      real (kind=R8) :: dx(ndepth-1), x0, xn
      real (kind=R8) :: x(ndepth), xx(ndepth_nml), uu(ndepth_nml)

      common /zzline/ dx,x0,xn,x
      common /phyparams/ Qtot,Tcbc,active

      KAb = matbulk_Kappa_KA(ibnd)
      KBb = matbulk_Kappa_KB(ibnd)
cxpb  Will need to update this later for a better computation
      if (coating_thickness(ibnd).gt.0.0_R8.and.
     &   (coating_thickness(ibnd)+layer_thickness(ibnd)).gt.0.0_R8)
     > then
        KAl = coating_Kappa_KA(ibnd)
        KBl = coating_Kappa_KB(ibnd)
      else
        KAl = matbulk_Kappa_KA(ibnd)
        KBl = matbulk_Kappa_KB(ibnd)
      endif
      if (coating_thickness(ibnd).gt.0.0_R8) then
        KAc = coating_Kappa_KA(ibnd)
        KBc = coating_Kappa_KB(ibnd)
        Tc = max(Tcbc, min( matsurf_vapor_temp(ibnd), (Tcbc + KBb/KAb) *
     .       dexp(-KAb*Qtot*(coating_thickness(ibnd)-xn)) - KBb/KAb ))
      endif
      do i = ndepth_nml, 2, -1
        xx(i) = x(i)
        if (target_depth(ibnd,i).ge.coating_thickness(ibnd)) then
          uu(i) = max(Tcbc, min( matsurf_vapor_temp(ibnd),
     .     (Tcbc + KBb/KAb) * dexp(-KAb*Qtot*(x(i)-xn)) - KBb/KAb ))
        else
          uu(i) = max(Tcbc, min( matsurf_vapor_temp(ibnd),
     .     (Tc + KBc/KAc) *
     .     dexp(-KAc*Qtot*(x(i)-coating_thickness(ibnd))) - KBc/KAc ))
        endif
c        print*,x(i),' ',uu(i),' ',i,' ',xn,' ',x0
      enddo
      xx(1) = x(1)
      Tl = uu(2)
      uu(1) = max(Tl, min( matsurf_vapor_temp(ibnd),
     .     (Tl + Kbl/KAl) *
     .     dexp(-KAl*Qtot*(x(1)-x(2))) - KBl/KAl ))

      return
      end subroutine fic

c_______________________________________________________________________
c Dirichlet condition(u) at xn. !Neumann condition(du/dx) at x0.
c No. of equations to be solved = ndepth_nml - 1
c_______________________________________________________________________
      subroutine lines2c(ibnd,ns,t,n,u,du)
      use b2mod_types
      use b2mod_wall

      implicit none

      integer ns,n,i,j,ibnd
      real (kind=R8) :: ux,uxx
      real (kind=R8) :: x(ndepth),u(ndepth_nml),uu(ndepth_nml+1),
     .                  du(ndepth_nml),t
      real (kind=R8) :: a0, b0, c0, an, bn, cn
      real (kind=R8) :: dx(ndepth-1), x0, xn, Ts, Td
      common /zzline/ dx,x0,xn,x

c      print*,t,' ',n,' InpVars ',u
c Maximum no. of equations that can be solved is np.

c Call to coefbc to get the b.cs
      Ts = u(1)
      Td = u(ndepth_nml)
      call coefbc(ibnd,ndepth_nml,ns,t,Ts,Td,a0,b0,c0,an,bn,cn)

c Neumann condition at x = x0.
c u0 = u(x0-dx,t)
      uu(1) = u(2) + 2.0_R8*dx(1)*(a0*u(1) - c0)/b0

c Storing intermediate values.
      do i = 1,ndepth_nml-1
        uu(i+1) = u(i)
      enddo

c Dirichlet condition(u) at xn.
      uu(ndepth_nml+1) = cn/an

c Calculate du/dt at all grid points.
      do j = 1,ndepth_nml-1
        i = j + 1
        ux=(uu(i+1) - uu(i-1))/(dx(i)+dx(i-1))
        uxx=((uu(i+1)-uu(i))/dx(i) - (uu(i)-uu(i-1))/dx(i-1))/
     .      ((dx(i)+dx(i-1))/2.0_R8)
        call fcn(ibnd,j,x(j),t,uu(i),ux,uxx,du(j))
      enddo

      return
      end subroutine lines2c

c..................................................................
c                     SUBROUTINE FCN

c  The diffusion equation in our case looks like ->
c       du = (1/rho/Cp(T)) * (K(T)*UXX + gradK(T)*UX + S)
c  where gradK = dK/dx and S is a heat flux source which is a delta
c  function at the plasma facing surface, which is implemented as a
c  boundary condition. Therefore here S = 0.
c  The second term on RHS is grad(K*UX) = dK/dx * UX = dK/dU*(UX)^2

      subroutine fcn(ibnd,i,x,t,u,ux,uxx,du)
      use b2mod_types
      use b2mod_wall

      implicit none
      integer ibnd, i
      real (kind=R8) :: x, t, u, du, ux, uxx
      real (kind=R8) :: Kappa,dKappabydT,Cp,mass_density, Ts

      Ts = u
      du = (1.0_R8/mass_density(ibnd,i)/Cp(ibnd,i,Ts))*
     &     (Kappa(ibnd,i,Ts)*uxx + dKappabydT(ibnd,i,Ts)*ux*ux)

      return
      end subroutine fcn
c______________________________________________________________________

c                 subroutine coefbc

c Subroutine to specify the boundary condition.
c general boundary condition.
c a0*u + b0*du/dx = c0  at the plasma facing surface.
c an*u + bn*du/dx = cn  at the cooling tube surface.

c t, Ts and Td are inputs to the subroutine and a0, b0, c0, an, bn, and cn
c are the outputs from the subroutine.

      subroutine coefbc(ibnd,nd,ns,t,Ts,Td,a0,b0,c0,an,bn,cn)
      use b2mod_types

      implicit none
      integer ibnd, nd, ns
      real (kind=R8) :: a0,b0,c0,an,bn,cn,t,Ts,Td
      real (kind=R8) :: Qtot, Tcbc
      logical active
      common /phyparams/ Qtot,Tcbc,active
      real (kind=R8) :: Kappa, qevap, qrad, cond_coef
      external Kappa, qevap, qrad, cond_coef

c  Coefficients at x=x0
      a0 = 0.0_R8
      b0 = 1.0_R8
      c0 = 1.0_R8/Kappa(ibnd,1,Ts)*
     . (-Qtot+qevap(ibnd,ns,Ts)+qrad(ibnd,Ts,.true.))

c  Coefficients at x=xn
      an = 0.0_R8
      bn = 1.0_R8
      cn = 1.0_R8/Kappa(ibnd,nd,Td)*
     . (-cond_coef(.not.active,Td)*(Td-Tcbc))
!! FIXME: Could use additionally the code below once we have a way of applying
!!        the correct emissivity
cxpb        if (.not.active)
cxpb     >    cn=cn-qrad(ibnd,Td,.false.)
cxpb     .        -qevap(ibnd,ns,Td)
c$$$ Old treatment where conductivity to coolant was assumed to be perfect
c$$$      if (active) then
c$$$        an = 1.0_R8
c$$$        bn = 0.0_R8
c$$$        cn = Tcbc
c$$$      else
c$$$        an = 0.0_R8
c$$$        bn = 1.0_R8
c$$$        cn = 1.0_R8/Kappa(ibnd,nd,Td)* !xpb we assume the contact heat conductivity of 100 W.m-2.K-1
c$$$     .   (-100.0_R8*(Td-Tcbc))
c$$$      endif

      return
      end subroutine coefbc

c.................................................................
c                      function Kappa
c  function Kappa(T) gives the coefficient of thermal conductivity
c  at temperature T.

      function Kappa(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Kappa, T
      real (kind=R8) :: a(5), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: layer_Kappa
      external layer_Kappa, xerrab

c Units of Kappa are in "W/m/K"
      a(1) = matbulk_Kappa_A1(ibnd)
      a(2) = matbulk_Kappa_A2(ibnd)
      a(3) = matbulk_Kappa_A3(ibnd)
      a(4) = matbulk_Kappa_A4(ibnd)
      a(5) = matbulk_Kappa_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Kappa) !')
      endif
      bulk_value = a(1) +
     . (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_A1(ibnd)
        a(2) = coating_Kappa_A2(ibnd)
        a(3) = coating_Kappa_A3(ibnd)
        a(4) = coating_Kappa_A4(ibnd)
        a(5) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (coat Kappa) !')
        endif
        coating_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa = coating_value
      else
        Kappa = layer_Kappa(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Kappa = Kappa+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Kappa = Kappa+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Kappa

c.................................................................
c                      function dKappabydT
c  function Kappa(T) gives the coefficient of thermal conductivity
c  at temperature T. We use the analytical expression resulting
c  from the first derivative of Kappa above.

      function dKappabydT(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer ibnd, i
      real (kind=R8) :: dKappabydT, T
      real (kind=R8) :: a(4), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: dlayer_KappabydT
      external dlayer_KappabydT, xerrab

      a(1) = matbulk_Kappa_A2(ibnd)
      a(2) = matbulk_Kappa_A3(ibnd)
      a(3) = matbulk_Kappa_A4(ibnd)
      a(4) = matbulk_Kappa_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk dKappadT) !')
      endif
      bulk_value = a(1) +
     . (2.0_R8*a(2) + (3.0_R8*a(3) + 4.0_R8*a(4) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_A2(ibnd)
        a(2) = coating_Kappa_A3(ibnd)
        a(3) = coating_Kappa_A4(ibnd)
        a(4) = coating_Kappa_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat dKappadT) !')
        endif
        coating_value = a(1) +
     .   (2.0_R8*a(2) + (3.0_R8*a(3) +
     .                   4.0_R8*a(4) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappabydT = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappabydT = coating_value
      else
        dKappabydT = dlayer_KappabydT(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          dKappabydT = dKappabydT+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          dKappabydT = dKappabydT+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function dKappabydT

c.................................................................
c                      function Kappa_perp
c  function Kappa_perp(T) gives the coefficient of perpendicular
c  thermal conductivity at temperature T.

      function Kappa_perp(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Kappa_perp, T
      real (kind=R8) :: a(5), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: layer_Kappa
      external layer_Kappa, xerrab

c Units are in "W/m/K"
      a(1) = matbulk_Kappa_perp_A1(ibnd)
      a(2) = matbulk_Kappa_perp_A2(ibnd)
      a(3) = matbulk_Kappa_perp_A3(ibnd)
      a(4) = matbulk_Kappa_perp_A4(ibnd)
      a(5) = matbulk_Kappa_perp_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Kappa) !')
      endif
      bulk_value = a(1) +
     . (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_perp_A1(ibnd)
        a(2) = coating_Kappa_perp_A2(ibnd)
        a(3) = coating_Kappa_perp_A3(ibnd)
        a(4) = coating_Kappa_perp_A4(ibnd)
        a(5) = coating_Kappa_perp_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat Kappa) !')
        endif
        coating_value = a(1) +
     .   (a(2) + (a(3) + (a(4) +  a(5) * Teff) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa_perp = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Kappa_perp = coating_value
      else
c-xpb  The layer is expected to be amorphous and therefore not to support
c-xpb  the distinction between perp and parallel directions
        Kappa_perp = layer_Kappa(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Kappa_perp = Kappa_perp+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Kappa_perp = Kappa_perp+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Kappa_perp

c.................................................................
c                      function dKappa_perpbydT
c  function Kappa_perp(T) gives the coefficient of thermal conductivity
c  at temperature T in the perpendicular direction. We use the analytical
c  expression resulting from the first derivative of Kappa_perp above.

      function dKappa_perpbydT(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: dKappa_perpbydT, T
      real (kind=R8) :: a(4), Tcels, Teff
      real (kind=R8) :: bulk_value, coating_value
      real (kind=R8) :: dlayer_KappabydT
      external dlayer_KappabydT, xerrab

      a(1) = matbulk_Kappa_perp_A2(ibnd)
      a(2) = matbulk_Kappa_perp_A3(ibnd)
      a(3) = matbulk_Kappa_perp_A4(ibnd)
      a(4) = matbulk_Kappa_perp_A5(ibnd)
      Tcels = T - 273.15_R8
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk dKappadT) !')
      endif
      bulk_value = a(1) +
     . (2.0_R8*a(2) + (3.0_R8*a(3) + 4.0_R8*a(4) * Teff) * Teff) * Teff
      if (coating_thickness(ibnd).gt.0.0_R8) then
        a(1) = coating_Kappa_perp_A2(ibnd)
        a(2) = coating_Kappa_perp_A3(ibnd)
        a(3) = coating_Kappa_perp_A4(ibnd)
        a(4) = coating_Kappa_perp_A5(ibnd)
        if (coating_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (coating_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab
     .     ('Unrecognized temperature units (coat dKappadT) !')
        endif
        coating_value = a(1) +
     .   (2.0_R8*a(2) + (3.0_R8*a(3) +
     .                   4.0_R8*a(4) * Teff) * Teff) * Teff
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappa_perpbydT = bulk_value
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        dKappa_perpbydT = coating_value
      else
c-xpb  The layer is expected to be amorphous and therefore not to support
c-xpb  the distinction between perp and parallel directions
        dKappa_perpbydT = dlayer_Kappabydt(ibnd,T)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          dKappa_perpbydT = dKappa_perpbydT+bulk_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          dKappa_perpbydT = dKappa_perpbydT+coating_value*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function dKappa_perpbydT

c....................................................................
c                      function Cp
c  function Cp(T) gives the specific heat capacity at const pressure
c  at temperature T.

      function Cp(ibnd,i,T)
      use b2mod_types
      use b2mod_wall
      use b2mod_layer
      implicit none
      integer ibnd, i
      real (kind=R8) :: Cp, T
      integer mf, cf
      real (kind=R8) :: b(4), Tcels, Teff
      real (kind=R8) :: Cp_bulk, Cp_coating
      real (kind=R8) :: layer_Cp
      external layer_Cp, xerrab

      Tcels = T - 273.15_R8
      b(1) = matbulk_Cp_B1(ibnd)
      b(2) = matbulk_Cp_B2(ibnd)
      b(3) = matbulk_Cp_B3(ibnd)
      b(4) = matbulk_Cp_B4(ibnd)
      mf = matbulk_format(ibnd)
      if (bulk_temp_unit(ibnd).eq.'K') then
        Teff = T
      else if (bulk_temp_unit(ibnd).eq.'C') then
        Teff = Tcels
      else
        call xerrab ('Unrecognized temperature units (bulk Cp) !')
      endif
      if (mf.eq.1) then
        Cp_bulk = b(1) + ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
      else
        Cp_bulk = b(1) + b(2) * Teff - b(3)/Teff**2
      end if

      if (coating_thickness(ibnd).gt.0.0_R8) then
        b(1) = coating_Cp_B1(ibnd)
        b(2) = coating_Cp_B2(ibnd)
        b(3) = coating_Cp_B3(ibnd)
        b(4) = coating_Cp_B4(ibnd)
        cf = coating_format(ibnd)
        if (bulk_temp_unit(ibnd).eq.'K') then
          Teff = T
        else if (bulk_temp_unit(ibnd).eq.'C') then
          Teff = Tcels
        else
          call xerrab ('Unrecognized temperature units (coat Cp) !')
        endif
        if (cf.eq.1) then
          Cp_coating =
     .     b(1) + ( b(2) + ( b(3) + b(4) * Teff) * Teff ) * Teff
        else
          Cp_coating = b(1) + b(2) * Teff - b(3)/Teff**2
        end if
      end if
      if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .    target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Cp = Cp_bulk
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        Cp = Cp_coating
      else
        Cp = (layer_Cp(ibnd,T))*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          Cp = Cp+Cp_bulk*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          Cp = Cp+Cp_coating*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function Cp

c....................................................................
c                      function mass_density
c  function mass_density(ibnd,i) gives the mass density of the target
c  element ibnd at depth i.

      function mass_density(ibnd,i)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer ibnd, i
      real (kind=R8) :: mass_density

      if (i.gt.ndepth_nml) then
        mass_density = matbulk_target_density(ibnd)
      else if (target_depth(ibnd,i).ge.coating_thickness(ibnd).and.
     .         target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        mass_density = matbulk_target_density(ibnd)
      else if (target_depth(ibnd,i).ge.layer_depth(ibnd)) then
        mass_density = coating_target_density(ibnd)
      else
        mass_density = layer_mass_density(ibnd)*
     .   (layer_depth(ibnd)-target_depth(ibnd,i))/
     .   (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        if (layer_depth(ibnd).ge.coating_thickness(ibnd)) then
          mass_density = mass_density+
     .      matbulk_target_density(ibnd)*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        else
          mass_density = mass_density+
     .      coating_target_density(ibnd)*
     .     (target_depth(ibnd,i+1)-layer_depth(ibnd))/
     .     (target_depth(ibnd,i+1)-target_depth(ibnd,i))
        end if
      end if

      return
      end function mass_density

c Getting the heat loss due to evaporation and radiation from surface.
c....................................................................
      function qrad(ibnd,Ts,top_layer)
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
      implicit none
      integer ibnd
#     include "emissivity.constants"
      real (kind=R8) :: Ts, Tw, SBconst, qrad
      logical top_layer
      parameter (SBconst =
     . pi**5 * (4.0_R8/15.0_R8) * (KBolt/const_h)**3
     . * KBolt / c**2)
      real (kind=R8) :: layer_emissivity
      external layer_emissivity

      Tw = backplate_temp(ibnd)      ! Temperature of surrounding environment
      if (layer_depth(ibnd).eq.0.0_R8 .or. .not.top_layer) then
        qrad = matsurf_blackbody_emiss(ibnd)*SBconst*(Ts**4-Tw**4)
      else
        qrad = layer_emissivity(ibnd)*SBconst*(Ts**4-Tw**4)
      endif

      return
      end function qrad

c....................................................................
c   Ref: 1) C. Garc'ia Rosales, Erosion processes in plasma-wall
c           interactions, J. Nucl. Mater. 211 (1994) Pg.202-214.
c        2) R.A. Langley, Evaporation, Data Compendium for Plasma
c           Surface Interactions, Special Issue, NUCLEAR FUSION,
c           IAEA, Vienna, (1984) 55-61.
c  Table of A and B for some first wall materials in Fusion research:
c  (Ref (2) above).
c  (note that A and B are provided such that the vappres is in Pascals.
c  element    |  A       |  B      |
c  Aluminum   | -16165   | 10.91   |
c  Beryllium  | -16720   | 11.61   |
c  Graphite   | -40181   | 14.80   |
c  Copper     | -17079   | 11.26   |
c  Molybdenum | -33069   | 12.02   |
c  Titanium   | -23340   | 9.79    |
c  Tungsten   | -44485   | 12.74   |
c  Titanium - |          |         |
c     Carbide | -39117.9 | 15.23   |
c  Stainless -|          |         |
c     Steel   | -14321.9 | 9.44    |

      function qevap(ibnd,ns,Ts)
      use b2mod_types
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_sputter
      implicit none
      integer ibnd, ns, is, itrack
      logical, save :: lfirst
      real (kind=R8) :: Ts, Tv, flxthev, qevap, Vapor_temp, qthev
      real (kind=R8), save :: therm_evap
      data lfirst /.true./, therm_evap /0.0_R8/
      external subini, subend, ipgetr
      intrinsic min

c  Ts       ->  temperature of the Graphite plate in Kelvin.
c  flxthev  ->  Flux of target atoms evaporated ( / m^2 /s )
c  qevap -> Energy lost (Joules / m^2 / s)

      call subini ('qevap')
      if (lfirst) then
        call ipgetr('b2stbr_therm_evap', therm_evap)
        call xertst(0.0_R8.le.therm_evap.and.
     &              1.0_R8.ge.therm_evap,
     .    'faulty parameter therm_evap')
        lfirst = .false.
      endif

      do is = 0, ns-1
        itrack = track_index(is)
        if (is_neutral(is)) then
c  Checking inputs
          if (layer_depth(ibnd).eq.0.0_R8 .or. itrack.eq.0) then
            Vapor_temp = matsurf_vapor_temp(ibnd)
            Tv = min(Vapor_temp, Ts)
          else if (itrack.gt.0 .and. itrack.le.ntrack) then
            Vapor_temp = species_vapor_temp(itrack)
            Tv = min(Vapor_temp, Ts)
          else
            Tv = Ts
          endif
          if (Ts.gt.Vapor_temp) write(*,'(a,1p,2e12.5)')
     >     'Input error; Correct so that Ts <= Vapor_temp',
     .                                   Ts,   Vapor_temp
          call thermev(ibnd,is,Tv,flxthev,qthev)
          qevap = qthev * therm_evap
        endif
      enddo

      call subend ()
      return
      end function qevap

      subroutine pdedef(npde, T, X, U, UX, P, Q, R, ires)
      use b2mod_types
      use b2mod_wall
      implicit none
      integer npde, ires
      real (kind=R8) :: T, X, U(NPDE), UX(NPDE),
     .  P(NPDE,NPDE), Q(NPDE), R(NPDE)
      integer i, ib, it, ih

      real (kind=R8) :: Cp, Kappa, mass_density
      external Cp, Kappa, mass_density

      logical active
      real (kind=R8) :: Qtot, Tcbc
      common /phyparams/ Qtot,Tcbc,active
      integer :: iwall, ns
      common /location/ iwall, ns

      if (X.le.target_depth(iwall,1)) then
        ih = 1
      else if (X.ge.target_depth(iwall,ndepth_nml)) then
        ih = ndepth_nml
      else
        ib = 1
        it = ndepth_nml
        do while (it-ib.gt.1)
          ih = (ib+it)/2
          if (X.lt.target_depth(iwall,ih)) then
            it = ih
          else if (X.gt.target_depth(iwall,ih)) then
            ib = ih
          else
            ib = ih
            it = ih
          endif
        end do
      end if

      do i = 1, npde
        P(i,1:npde) = mass_density(iwall, ih) * Cp(iwall, ih, U(i))
        Q(i) = 0.0_R8
        R(i) = Kappa(iwall, ih, U(i)) * UX(i)
      enddo

      return
      end subroutine pdedef

      subroutine bndary (npde, T, U, UX, iposb, beta, gamma, ires)
      use b2mod_types
      implicit none
      integer npde, iposb, ires
      real (kind=R8) :: T, U(npde), UX(npde), beta(npde), gamma(npde)
      real (kind=R8) :: Qtot, Tcbc, cond_coef
      logical active

      real (kind=R8) :: qrad, qevap
      external qrad, qevap

      common /phyparams/ Qtot,Tcbc,active
      integer :: iwall, ns
      common /location/ iwall, ns

      if (iposb.eq.0) then  ! LHS BC : x=0 : Heat flux imposed
        beta(1)=-1.0_R8
        gamma(1)=Qtot-qrad(iwall,U(1),.true.)-qevap(iwall,ns,U(1))
      else if (iposb.eq.1) then ! RHS BC : Cooling with some contact heat conductivity (low for inertial cooling, high for active cooling)
        beta(1)= 1.0_R8
        gamma(1)=-cond_coef(.not.active,U(1))*(U(1)-Tcbc)
!! FIXME: Could use additionally the code below once we have a way of applying 
!!        the correct emissivity
cxpb        if (.not.active) 
cxpb     >    gamma(1)=gamma(1)-qrad(iwall,U(1),.false.)
cxpb     .                     -qevap(iwall,ns,U(1))
c$$$ Old treatment where conductivity to coolant was assumed to be perfect
c$$$      else if (iposb.eq.1.and.active) then ! RHS BC : x=xn : Cooling temperature imposed
c$$$        beta(1)= 0.0_R8
c$$$        gamma(1)= U(1)-Tcbc
c$$$      else if (iposb.eq.1.and..not.active) then ! RHS BC : x=xn : Inertial Cooling
c$$$        beta(1)= 1.0_R8
c$$$        gamma(1)= -100.0_R8*(U(1)-Tcbc) ! we assume the contact heat conductivity of 100 W.m-2.K-1
      end if

      return
      end subroutine bndary

      function cond_coef(inertial,T)
      use b2mod_types
      implicit none
      real(kind=R8) :: cond_coef
      logical inertial
      real(kind=R8) :: T, Tc

      Tc = T - 273.15_R8

      if (inertial) then
! we assume a residual contact heat conductivity of 100 W.m-2.K-1
        cond_coef = 100.0_R8
      else
! The following are fits to the result obtained by Nukiyama for the heat conduction
! transfer coefficient between the CuCrZr pipe and the water coolant,
! applicable to the ITER divertor targets
! Full reference can be found at: https://user.iter.org/?uid=35L29Y
!      cond_coef = 1.64908e-8_R8*Tc**6 - 1.24649e-5_R8*Tc**5
!     @       + 3.4938e-3_R8*Tc**4 - 0.44568_R8*Tc**3 + 25.5532_R8*Tc**2
!     @       - 502.333_R8*Tc + 77703.2_R8
!      cond_coef = 1.357954e-3_R8*Tc**3 + 6.616505e-1_R8*Tc**2 - 24.41418_R8*Tc
!     @       + 7.725710e4_R8
        cond_coef = 8.0e4_R8
      end if

      return
      end function cond_coef
