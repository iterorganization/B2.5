*-----------------------------------------------------------------------
*.specification

      subroutine b2t21t(nx, ny, ns, lnlam, bb, vol, hx, pbs,        !srv 01.07.09
     1 ne, na, ne2, te, ti, ua, fch_p, chce, chci, cthe, cthi)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  lnlam(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny),     !srv 01.07.09
     &  hx(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1), fch_p(-1:nx,-1:ny,0:1), 
     &  na(-1:nx,-1:ny,0:ns-1), 
     &  ua(-1:nx,-1:ny,0:ns-1), te(-1:nx,-1:ny), ti(-1:nx,-1:ny), 
     &  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  chce(-1:nx,-1:ny,0:1), chci(-1:nx,-1:ny,0:1),
     &  cthe(-1:nx,-1:ny,0:ns-1), cthi(-1:nx,-1:ny,0:ns-1)

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*       parallel thermal force, heat flux, and heat conduct. coeff.
*
*       uses subroutine THERMF
*                                                       Andreas Bergmann
*
*-----------------------------------------------------------------------
*.declarations

c     Variables in COMPLASMA:
c                 hfe, hfi:	poloidal comp. of par. heat flux
c                 hce, hci:	poloidal heat conduction coeff.
c                 thfi:	        parallel thermal force

c     local variables:

      real (kind=R8) ::
     *   xscal, xpos
      parameter (xscal=100.0_R8, xpos=0.1_R8)

      real (kind=R8) ::
     *          u(ns+1), n(ns+1), r(ns+1), mm(ns+1),
     *          tef, tif, gte, gti, he, hi, ce, ci, pitface,
     *          t0, t0i, t0e, cond, fac, 
     *          z(-1:nx,-1:ny,ns+1), z2(-1:nx,-1:ny,ns+1)

      integer 	nz(ns+1)
      integer 	i, ix, iy, ifail, nspec, ncharged

      external subini, subend, xertst, b2xvff, b2xvsg, ipgeti

      integer, save :: icase = 0
      integer ncall, n2, failcount, failcount0
      save ncall, failcount
      data ncall/0/, failcount/0/

*=======================================================================
*//COMPUTATION//

*  ..start-up calls
      call subini ('b2t21t')
*  ..test nx, ny, ns
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx,ny')
      call xertst (1.le.ns, 'faulty argument ns')
*  ..set n2
      n2=(nx+2)*(ny+2)
*  ..read switches
      if(ncall.eq.0) then
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
      endif
*  ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, bb
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test sign of na, ne, te, ti, ne2 
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
*    ..test edge values of fch_p
       call b2xvff (nx, ny, fch_p, 'fch_p')
      endif

*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)              !srv 20.09.11

      failcount0=failcount
      mm(1) = me
      nz(1) = 1
      nspec=2
      nz(nspec)=0
      mm(nspec)=am(0)*mp 	 ! xpb <------- new line
      do i=1,ns-1
        if(.not.is_neutral(i-1)) then
          nz(nspec)=nz(nspec)+1
          if(zn(i).ne.zn(i-1)) then
            mm(nspec)=am(i-1)*mp
            nspec=nspec+1
            nz(nspec)=0
            mm(nspec)=am(i)*mp
          endif
        endif
      enddo
      nz(nspec)=nz(nspec)+1

      z(-1:nx,-1:ny,1) = 1.e0_R8
      z2(-1:nx,-1:ny,1) = 1.e0_R8

      ncharged=1
      do i=0,ns-1
        if(.not.is_neutral(i)) then
          ncharged=ncharged+1
          z(-1:nx,-1:ny,ncharged) = rza(-1:nx,-1:ny,i)
          z2(-1:nx,-1:ny,ncharged) = rz2(-1:nx,-1:ny,i)
        endif
      end do

      do  iy=-1,ny
        do  ix=-1,nx
          if(leftix(ix,iy).ne.-2) then
            pitface=(abs(bb(leftix(ix,iy),leftiy(ix,iy),0))/
     1       bb(leftix(ix,iy),leftiy(ix,iy),3)*hx(ix,iy)+
     2       abs(bb(ix,iy,0))/bb(ix,iy,3)*
     2       hx(leftix(ix,iy),leftiy(ix,iy)))/
     3       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))

c         parallel temperature gradients
c         ------------------------------
            gte = (te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))*
     1       pitface/(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8
            gti = (ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))*
     1       pitface/(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8

c	  temperature at cell faces
c         -------------------------
c !? seems to have a different averaging scheme for te/i and ne/i
            tef = (te(leftix(ix,iy),leftiy(ix,iy))*
     1       hx(leftix(ix,iy),leftiy(ix,iy)) + te(ix,iy)*hx(ix,iy))/
     *       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            tif = (ti(leftix(ix,iy),leftiy(ix,iy))*
     1       hx(leftix(ix,iy),leftiy(ix,iy)) + ti(ix,iy)*hx(ix,iy))/
     *       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
c         (This weighting reduces the influence of a very low
c          temperature near the target plates, where hx is very
c          small -- same as aver. of hcx in b2cene, s. below)

c         densities at cell faces and velocities
c         --------------------------------------
            n(1) = (ne(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy)+
     1       ne(ix,iy)*hx(leftix(ix,iy),leftiy(ix,iy)))/
     *       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            u(1) = 0.e0_R8
            ncharged=1
            do i=0,ns-1
              if(.not.is_neutral(i)) then
                ncharged=ncharged+1
                u(ncharged) = 
     1           (ua(leftix(ix,iy),leftiy(ix,iy),i)*hx(ix,iy)+
     1            ua(ix,iy,i)*hx(leftix(ix,iy),leftiy(ix,iy)))/
     *           (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
                n(ncharged) =
     1           (na(leftix(ix,iy),leftiy(ix,iy),i)*hx(ix,iy)+
     1            na(ix,iy,i)*hx(leftix(ix,iy),leftiy(ix,iy)))/
     *           (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
                u(1) = u(1) + u(ncharged)*z(ix,iy,ncharged)*n(ncharged)
              endif
            end do
            u(1) = u(1) + 
     .       (ue_ext(leftix(ix,iy),leftiy(ix,iy))*
     .        ne_ext(leftix(ix,iy),leftiy(ix,iy))*hx(ix,iy)+
     .        ue_ext(ix,iy)*ne_ext(ix,iy)*
     .        hx(leftix(ix,iy),leftiy(ix,iy)))/
     .       (hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))
            u(1) = u(1) - fch_p(ix,iy,0)/pbs(ix,iy,0)/qe  ! parallel current
c !? do we need to apply a 21-moment ansatz for 
c    the electrical conductivities
c
c         electron velocity
c         -----------------
            u(1) = u(1)/n(1)
        
c           -------------------------------------------------------------
c                 solve 21-moments equations
c           -------------------------------------------------------------
c       write(*,*) 'ix,iy',ix,iy
c       write(*,*) 'u',u
c       write(*,*) 'n',n
c       write(*,*) 'z',z
c       write(*,*) 'mm',mm
c       write(*,*) 'tef,tif',tef,tif
c       write(*,*) 'gte,gti',gte,gti
c       write(*,*) 'nz',nz
c       write(*,*) 'nspec',nspec
c       write(*,*) 'ncharged',ncharged
          call thermf( u, n, z(ix,iy,1), z2(ix,iy,1), 
     &     mm, tef, tif, gte, gti, lnlam(ix,iy),                     !srv 01.07.09
     *     nz, nspec, ncharged, ce, ci, he, hi, r, ifail, cond )
          fac=(atan(xscal*(cond-xpos))-atan(-xscal*xpos))/
     1     (atan(1.0e8_R8)-atan(-xscal*xpos))
c       write(*,*) 'b2t21t:ix,iy,cond,fac',ix,iy,cond,fac

          if( ifail .eq. 0 )  then
c            poloidal heat flux and heat conduction coefficient
c              --------------------------------------------------

c         write(*,*) 'ix,iy,old_chce(ix,iy,0),new_chce(ix,iy,0),ce',
c     1          ix,iy,chce(ix,iy,0),he/gte*pitface**2,ce*pitface**2
c         write(*,*) 'ix,iy,old_chci(ix,iy,0),new_chci(ix,iy,0),ci',
c     1     ix,iy,chci(ix,iy,0),hi/gti*pitface**2,ci*pitface**2

            chce(ix,iy,0) = fac*ce*pitface**2 + (1-fac)*chce(ix,iy,0)
            chci(ix,iy,0) = fac*ci*pitface**2 + (1-fac)*chci(ix,iy,0)

c         write(*,*) 'ix,iy,chce(ix,iy,0)',ix,iy,chce(ix,iy,0)
c         write(*,*) 'ix,iy,chci(ix,iy,0)',ix,iy,chci(ix,iy,0)

c            parallel thermal force with contribution from electron
c            therm. force (parallel electric field is eliminated)
c              --------------------------------------------------------

            if(.true.) then
              ncharged=1
              do i=0,ns-1
                if(.not.is_neutral(i)) then
                  ncharged=ncharged+1
                  t0 = -
     1             (na(leftix(ix,iy),leftiy(ix,iy),i)+na(ix,iy,i))*
     &             ((rza(ix,iy,i)+rza(leftix(ix,iy),leftiy(ix,iy),i))/
     &             2.0_R8/(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))-
     &             (rz2(ix,iy,i)+rz2(leftix(ix,iy),leftiy(ix,iy),i))/
     &             2.0_R8/(ne2(leftix(ix,iy),leftiy(ix,iy))+ne2(ix,iy)))
c
c             write(*,*) 'ix,iy,is,thfi,thfi_s',ix,iy,i,
c     1         r(ncharged)+n(ncharged)*z(ix,iy,ncharged)/n(1)*r(1),
c     2         t0*0.5_R8*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*
c     3         (cthe(ix,iy,i)*gte+cthi(ix,iy,i)*gti)
c
c electrons
                  t0e=t0*gte
                  if(t0e.ne.0) then
                    cthe(ix,iy,i)=
     .               fac * (n(ncharged)*z(ix,iy,ncharged)/n(1)*r(1) /
     1               (0.5_R8*
     2               (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*t0e))+
     2               (1-fac) * cthe(ix,iy,i)
c *
c     1         pitface/(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8
c ions
                  endif
                  t0i=t0*gti
                  if(t0i.ne.0) then
                    cthi(ix,iy,i) = fac * (r(ncharged) / 
     1               (0.5_R8*
     2               (ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*t0i)) +
     2               (1-fac) * cthi(ix,iy,i)
c *
c     1         pitface/(hx(leftix(ix,iy),leftiy(ix,iy))+hx(ix,iy))*2.0_R8
c
                  endif
                  if(cthi(ix,iy,i)*cthe(ix,iy,i).lt.-1e+3_R8.and.
     1             abs(gti/tif).gt.1e-6_R8) then
c               write(*,*) 'ix,iy,is,cthe_p,cthi_p',
c     1           ix,iy,i,cthe(ix,iy,i),cthi(ix,iy,i)
c-x               cthi(ix,iy,i)=cthi(ix,iy,i)+cthe(ix,iy,i)*gte/gti
c-x               cthe(ix,iy,i)=0
c               write(*,*) 'ix,iy,is,cthe_r,cthi_r',
c     1           ix,iy,i,cthe(ix,iy,i),cthi(ix,iy,i)
                  endif
c             write(*,*) 'ix,iy,is,t0,gte,gti,cthe,cthi',
c     1         ix,iy,i,t0,gte,gti,cthe(ix,iy,i),cthi(ix,iy,i)
c             write(*,*) 'ix,iy,is,thfi,thfi_c',ix,iy,i,
c     1         r(ncharged)+n(ncharged)*z(ix,iy,ncharged)/n(1)*r(1),
c     2         t0*0.5_R8*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))*
c     3         (cthe(ix,iy,i)*gte+cthi(ix,iy,i)*gti)
                endif
              end do
            endif
          else
            if(failcount.lt.1000) then
              write(6,*)  ' THFHC:  matrix inversion failed, code ',
     *         ifail,', Braams formula used, ix,iy = ',ix,iy
            elseif(failcount.eq.1000) then
              write(6,*)  ' THFHC:  matrix inversion failed, ',
     *         'Further messages suppressed'
            endif
            failcount=failcount+1
          end if
c          ( ifail )        
        endif
      end do
c         ( iy )        
      end do
c       ( ix )

      write(*,*) 'Failures per b2t21t call = ',failcount-failcount0

      ncall = ncall+1
      call subend ()
      return
      end subroutine b2t21t

*-----------------------------------------------------------------------
*.specification

      subroutine b2t21v (nx, ny, ns, lnlam, bb, vol, hx, pbs,        !srv 01.07.09
     1 ne, na, te, ti, ua, cvsa)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  lnlam(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), vol(-1:nx,-1:ny),    !srv 01.07.09
     &  hx(-1:nx,-1:ny), pbs(-1:nx,-1:ny,0:1),
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  cvsa(-1:nx,-1:ny,0:1,0:ns-1)

*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*       classical parallel (multifluid) viscosity
*
*       uses subroutine VISCALC
*                                                       Andreas Bergmann
*
*-----------------------------------------------------------------------
*.declarations

      real (kind=R8) ::
     *                  wwpp(ns), n(ns), z(-1:nx,-1:ny,ns), visc(ns), 
     *                  eta(ns), mm(ns), coef, pit, cond,
     *                  z2(-1:nx,-1:ny,ns)
      integer	        nz(ns)
      integer	        i, ix, iy, ifail, nspec, ncharged

      external subini, subend, xertst, b2xvsg

      integer, save :: icase = 0
      integer ncall,n2,failcount,failcount0
      save failcount
      data ncall/0/, failcount /0/

*=======================================================================
*//COMPUTATION//

*  ..start-up calls
      call subini ('b2t21v')
*  ..test nx, ny, ns
      call xertst (0.le.nx .and. 0.le.ny, 'faulty argument nx,ny')
      call xertst (1.le.ns, 'faulty argument ns')
*  ..set n2
      n2=(nx+2)*(ny+2)
*  ..read switches
      if(ncall.eq.0) then
       call ipgeti ('b2tlnl_ii', icase)
       call xertst (0.eq.icase .or. 1.eq.icase, 
     &               'faulty input b2tlnl_ii')
      endif
*  ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of vol, hx, bb
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test sign of na, ne, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif

*   ..compute the Coulomb logarithm
      call b2tlnl (nx, ny, icase, te, ti, ne, lnlam)              !srv 20.09.11

      failcount0=failcount
      nspec=1
      nz(nspec)=0
      mm(nspec)=am(0)*mp 	 ! xpb <------- new line
      do i=1,ns-1
        if(.not.is_neutral(i-1)) then
          nz(nspec)=nz(nspec)+1
          if(zn(i).ne.zn(i-1)) then
            mm(nspec)=am(i-1)*mp
            nspec=nspec+1
            nz(nspec)=0
            mm(nspec)=am(i)*mp
          endif
        endif
      enddo
      nz(nspec)=nz(nspec)+1

      ncharged=0
      do i=0,ns-1
        if(.not.is_neutral(i)) then
          ncharged=ncharged+1
          z(-1:nx,-1:ny,ncharged) = rza(-1:nx,-1:ny,i)
          z2(-1:nx,-1:ny,ncharged) = rz2(-1:nx,-1:ny,i)
        endif
      end do

      coef = 2.e0_R8
c       as in phythc, but this value is questionable (div u != 0)
c       (w_pp = coef * grad_|| up,  e.g. coef=2, if div u = 0;
c        coef=4/3, if div u = grad_|| up)

      do  iy=-1,ny
        do  ix=-1,nx
          if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1) then
            pit=abs(bb(ix,iy,0))/bb(ix,iy,3)

            ncharged=0
            do i=0,ns-1

              if(.not.is_neutral(i)) then
                ncharged=ncharged+1
                n(ncharged) = na(ix,iy,i)
C                       parallel strain (velocity gradient)
                wwpp(ncharged) = coef*
     1           (ua(rightix(ix,iy),rightiy(ix,iy),i)-ua(ix,iy,i))/
     2           (0.5_R8*(hx(rightix(ix,iy),rightiy(ix,iy))+hx(ix,iy)))
              endif
            end do

c            -------------------------------------------------------------
c                 solve 21-moments equations
c            -------------------------------------------------------------
c         write(*,*) 'ix,iy',ix,iy
c         write(*,*) 'n',n
c         write(*,*) 'z',z(ix,iy,1:ncharged)
c         write(*,*) 'z2',z2(ix,iy,1:ncharged)
c         write(*,*) 'mm',mm
c         write(*,*) 'ti(ix,iy)',ti(ix,iy)
c         write(*,*) 'nz',nz
c         write(*,*) 'nspec,ncharged',nspec,ncharged
c         write(*,*) 'wwpp',wwpp
            call viscalc(n, z(ix,iy,1), z2(ix,iy,1), mm, ti(ix,iy), nz, 
     *       nspec, ncharged, wwpp, visc, eta, lnlam(ix,iy), ifail,   !srv 01.07.09
     *       cond)
c         write(*,*) 'b2t21v:ix,iy,cond',ix,iy,cond

            if( ifail .eq. 0 )  then

              ncharged=0
              do i=0,ns-1
                if(.not.is_neutral(i)) then
                  ncharged=ncharged+1
                  cvsa(ix,iy,0,i)=eta(ncharged) * pit**2 * coef 
                endif
              end do

            else
c            -------------------------------------------------------------
c                 ( ifail != 0 )        
c               switch to previous simpler calculation, when matrix
c               inversion with NAG routine fails
c            -------------------------------------------------------------
              if(failcount.lt.1000) then
                write(6,*)  ' VISC:  matrix inversion failed, code ',
     *           ifail,', Braams formula used, ix,iy = ',ix,iy
              elseif(failcount.eq.1000) then
                write(6,*)  ' VISC:  matrix inversion failed, ',
     *          'Further messages suppressed'
              endif
              failcount=failcount+1
            end if
          endif
        end do
c        ( iy )        
      end do
c       ( ix )        

      write(*,*) 'Failures per b2t21v call = ',failcount-failcount0

      ncall = ncall+1
      call subend ()
      return
      end subroutine b2t21v
c
*//VISCALC//
*=======================================================================
*                      S U B R O U T I N E   V I S C A L C
*=======================================================================
*
*       VISCALC calculates the parallel viscosity for a multi-fluid
*       plasma (a fluid for each charge state) at one position in
*       21 moments approximation according to Igitkhanov/Zhdanov.
*       Required input are masses, charges, parallel strain (~ velocity
*       gradient), number densities, ion temperature.
*
*       The INPUT VARIABLES ARE REQUIRED TO BE and the output var. are
*       combined into vectors, e.g. n = (n11, n12, ..., n1nz(1), n21,
*       ..., n2nz(2), ...), where the species number is the first index,
*       and the charge state the second index.
*
*       Input:
*       --------
*       n  = ( ns  (s=1,...,nsz) ) :    density                 in m-3
*       z  = ( ns  (s=1,...,nsz) ) :  	charge                  in e
*       z2 = ( ns  (s=1,...,nsz) ) :  	square charge           in e**2
*       mm = ( ms  (s=1,...,ns) )  :  	mass                	in kg
*       ti                         :  	k*temperature           in J
*       w  = ( wps (s=1,...,nsz) ) :    velocity gradient       in s-1
*       nz = ( nzs (s=1,...,ns) )  :  	number of charge states
*       ns                         :  	number of species
*       nsz	                   :	total number of fluids
*
*       Output:
*       --------
*       pi  = ( pi_pp (s=1,...,nsz) ) :  parallel viscosity     in Pa
*       eta = ( eta_0 (s=1,...,nsz) ) :  parallel visc. coeff.  in Pa s
*
*       real: 	min. 8 Bytes precision
*
*       Subroutines used:
*       -----------------
*       F01AAF from the NAG library  (8 Byte reals)
*
*                                                       Andreas Bergmann
*=======================================================================
*//DECLARATIONS//

      subroutine viscalc(n, z, z2, mm, ti, nz, ns, nsz, w, pi_pp,        !xpb
     1 eta, lnlam, ifail, cond)                                          !srv 01.07.09
      use b2mod_types
      use b2mod_constants                                        	 !xpb

      implicit none

      integer           ns, nsz, nz(ns), ifail
      real (kind=R8) ::
     *                  n(nsz), z(nsz), z2(nsz), mm(ns),
     1 eta(nsz), pi_pp(nsz), w(nsz), ti	                                 !xpb

C       internal variables:
C       ww:     Sum  w*n    over Z
C       p:	Sum  n*T    over Z         
C       rho:	Sum  mm*n   over Z        
C       d:	rho/mm
C       z2n:	Sum  n*Z**2 over Z

      real (kind=R8) ::
     *            ww(ns), d(ns), p(ns), z2n(ns)

      real (kind=R8) ::
     *            g(ns,ns,0:6), mu(ns,ns), mu2(ns,ns)
      real (kind=R8) ::
     *            lnlam, const, cond, condval
      real (kind=R8) ::
     *            s1, s3a, s3b, s5, ss

      real (kind=R8) ::
     *             a(2*ns,2*ns), x(2*ns)
      real (kind=R8) ::
     *             inva(2*ns,2*ns), wk1(2*ns)
#ifndef NO_NAG
#ifndef NO_F01AAF
#else
      integer   IPIV(2*ns)
#endif
#else
      integer   IPIV(2*ns)
#endif

      integer          n1(ns)
      integer          i, j, k


*=======================================================================
*//COMPUTATION//

C       n1(j): 	index of vector element corresponding to first charge
C               state of species j in vectors n, z, w, pi_pp, eta
      n1(1)= 1
      do i=2,ns
        n1(i) = n1(i-1) + nz(i-1)
      end do

c       calculating  d(), z2n(), p(), ww()

      do j=1,ns
        d(j)   = n(n1(j))
        z2n(j) = n(n1(j))*z2(n1(j))
        ww(j)  = w(n1(j))*n(n1(j))
        do k=1, nz(j)-1
          d(j)   = d(j)   + n(n1(j)+k)
          z2n(j) = z2n(j) + n(n1(j)+k)*z2(n1(j)+k)
          ww(j)  = ww(j)  + n(n1(j)+k)*w(n1(j)+k)
        end do
        p(j)   = d(j)*ti
        ww(j)  = ww(j)/d(j)
      end do

      do j=1,ns
        do i=1,ns
          mu(i,j)  = mm(i)*mm(j)/(mm(i)+mm(j))
          mu2(i,j) = mu(i,j)/(mm(i)+mm(j))
        end do
      end do

c         const = sqrt(2)*e**4 / (12*pi**1.5_R8*eps0**2)
!xpb      const = 1.79e-55_R8
        const = sqrt(2.0_R8)*qe**4 / (12*pi**1.5_R8*eps0**2)		 !xpb
!xpb      lnlam = 12.e0_R8

      do j=1,ns
        do i=1,ns
          g(i,j,0) = const * z2n(i)*z2n(j) * sqrt(mu(i,j)) /
     *     ti**1.5_R8 * lnlam
          g(i,j,1) = -2.e0_R8/5.e0_R8*
     *     (5.e0_R8*mm(i)+3.e0_R8*mm(j))/(mm(i)+mm(j))*g(i,j,0)
          g(i,j,2) = -4.e0_R8/5.e0_R8 * mu(i,j)/mm(i) * g(i,j,0)
          g(i,j,3) =  6.e0_R8/35.e0_R8*
     *     (7.e0_R8*mm(i)+3.e0_R8*mm(j))/(mm(i)+mm(j))*g(i,j,0)
          g(i,j,4) = -6.e0_R8/7.e0_R8 * g(i,j,2)
          g(i,j,5) = - (51.e0_R8/35.e0_R8*mm(j)**2/mm(i) + 
     *     37.e0_R8/7.e0_R8*mm(j)
     *     + 22.e0_R8/5.e0_R8*mm(i) + 4.e0_R8*mm(i)**2/mm(j) )
     *     * mu2(i,j)/(mm(i)+mm(j)) * g(i,j,0)
          g(i,j,6) = - 24.e0_R8/7.e0_R8*mm(j)*mu2(i,j)/
     *     (mm(i)+mm(j))*g(i,j,0)
        end do
      end do

c  ---------------------------------------------------------------------
c       linear system  A * x = b with 2*ns  equations
c       x(i) = pi_pp(i)/rho(i), x(i+ns) = higher order moment, i=1...ns
c       b(i) = p(i)*ww(i), b(i+ns) = 0, i=1...ns
c  ---------------------------------------------------------------------
      do j=1,ns
        do i=1,ns
          a(i,j)       = - g(i,j,2)
          a(i+ns,j)    = - 7.e0_R8/2.e0_R8*mu(i,j)/mm(i) * g(i,j,4)
          a(i,j+ns)    = - g(i,j,4)/ti * mu(i,j)
          a(i+ns,j+ns) = - g(i,j,6)*mm(j)/ti
        end do
        do k=1,ns
          a(j,j)       = a(j,j)       + g(j,k,1)
          a(j+ns,j)    = a(j+ns,j)    + 7.0_R8/2.0_R8*mu(j,k)/
     1           mm(j)*g(j,k,3)
          a(j,j+ns)    = a(j,j+ns)    + g(j,k,3)/ti * mu(j,k)
          a(j+ns,j+ns) = a(j+ns,j+ns) + g(j,k,5)*mm(j)/ti
        end do
      end do

c  ---------------------------------------------------------------------
c       solving A*x=b by inverting A (Inverse matrix needed for coeff.)
c  ---------------------------------------------------------------------
      ifail = +1

#ifndef NO_NAG
#ifndef NO_F01AAF
      call f01aaf(a, 2*ns, 2*ns, inva, 2*ns, wk1, ifail)
#else
      call f07adf(2*ns, 2*ns, a, 2*ns, IPIV, ifail)
      if( ifail .eq. 0 )  then
        call f06qff('General', 2*ns, 2*ns, a, 2*ns, inva, 2*ns)
        call f07ajf(2*ns, inva, 2*ns, IPIV, wk1, 2*ns, ifail)
      endif
#endif
#else
      call sgetrf(2*ns, 2*ns, a, 2*ns, IPIV, ifail)
      if( ifail .eq. 0 )  then
        inva(1:2*ns,1:2*ns) = a(1:2*ns,1:2*ns)
        call sgetri(2*ns, inva, 2*ns, IPIV, wk1, 2*ns, ifail)
      endif
#endif
      cond=condval(a,2*ns,ifail)

      if( ifail .ne. 0 )  return
      
      do i=1,ns
        x(i) = inva(i,1)*p(1)*ww(1)
        do j=2,ns
          x(i) = x(i) + inva(i,j)*p(j)*ww(j)
        end do
      end do

c  ---------------------------------------------------------------------
c        viscosity and viscosity coefficient
c  ---------------------------------------------------------------------
      do i=1,ns
        s1 = 0.e0_R8
        do k=1,ns
          s1 = s1 + g(i,k,1)
        end do
        s3a = 0.e0_R8
        do k=1,ns
          s3a = s3a + mu(i,k)/ti*g(i,k,3)
        end do
        s3b = 0.e0_R8
        do k=1,ns
          s3b = s3b + 3.5_R8*mu(i,k)/mm(i)*g(i,k,3)
        end do
        s5 = 0.e0_R8
        do k=1,ns
          s5 = s5 + mm(i)*g(i,k,5)/ti
        end do
        ss = 1.e0_R8/(s1*s5 - s3a*s3b)
        do k=0,nz(i)-1
          pi_pp(n1(i)+k) = n(n1(i)+k)*mm(i) * ( x(i) + s5*ss*ti *
     *     (z2n(i)/z2(n1(i)+k)*w(n1(i)+k) - d(i)*ww(i)) )
          eta(n1(i)+k) = - n(n1(i)+k)*mm(i) * ( inva(i,i)*d(i)
     *     + s5*ss*(z2n(i)/z2(n1(i)+k)-d(i)) ) * ti
        end do
      end do

      return
      end subroutine viscalc
c
c NOTE from XPB, about the upcoming removal of F01AAF from the NAG library
c
c  Replace: CALL F01AAF(A,IA,N,X,IX,WKSPCE,IFAIL)
c  by     : CALL F07ADF(N,N,A,IA,IPIV,IFAIL)
c           CALL F06QFF('General',N,N,A,IA,X,IX)
c           CALL F07AJF(N,X,IX,IPIV,WKSPCE,LWORK,IFAIL)
c IPIV is an INTEGER array of length N, and
c LWORK is the dimension of WKSPCE (one must have LWORK .ge. MAX(1,N)).
c The routine just does a matrix inversion, and you will find all relevant
c information in section 4 of the preamble to chapter F01 in the NAG documentation.
c
c
*//THERMF//
*=======================================================================
*                     S U B R O U T I N E   T H E R M F
*=======================================================================
*
*       THERMF calculates the parallel heat flow and thermal force for
*       all the charge states of all species in 21 moments approximation
*       according to Igitkhanov/Zhdanov for given v_||, number densities,
*       temperatures, and grad_|| T.
*
*       The input VARIABLES ARE required to be and the output var. are
*       COMBINED TO VECTORS WITH THE ELECTRON VARIABLE FIRST, e.g.
*       n = (ne, ni11, ni12, ..., ni21, ni22, ... [first index: species
*       (H,D,C,O,...), second index: charge state number]).
*
*       Input:
*       --------
*       up = ( ups (s=1,...,nfl) ) :  	||-velocity 	        in m/s
*       n  = ( ns  (s=1,...,nfl) ) :  	density     	        in m-3
*       z  = ( ns  (s=1,...,nfl) ) :  	charge                  in e
*       te, ti 	                   :  	k*temperature 	        in J
*       gte, gti 	           :  	grad_|| kT	        in J/m
*       mm = ( ms  (s=1,...,ns) )  :  	mass	                in kg
*	nz = ( nzs (s=1,...,ns) )  :  	number of charge states
*                                       of species s
*       ns  	                   :  	number of species
*       nfl  	                   :  	number of fluids
*
*       Output:
*       --------
*       r = ( rth (s=1,...,nfl) )  :  	||-thermal force      in N m-3
*       he, hi   	           :  	||-heat flux	      in Jm-2s-1
*       ce, ci    	           :  	||-heat conductivity  in m-1s-1
*
*       real:   min. 8 Bytes precision
*
*       Subroutines used:
*       -----------------
*       F01AAF from the NAG library  (8 Byte reals)
*
*                                                       Andreas Bergmann
*=======================================================================
*//DECLARATIONS//

      subroutine thermf(up, n, z, z2, mm, te, ti, gte, gti, lnlam,    !srv 01.07.09
     *                  nz, ns, nfl, ce, ci, he, hi, r, ifail, cond)

      use b2mod_types
      use b2mod_constants	                                	 !xpb
      implicit none

        integer     ns, nfl, nz(ns)
        real (kind=R8) :: 
     *    up(nfl), n(nfl), z(nfl), z2(nfl), mm(ns), r(nfl)
        real (kind=R8) :: te, ti, gte, gti, lnlam, he, hi, ce, ci     !srv 01.07.09
        integer     ifail

C       internal variables:

C       heat flux h and a higher moment hh
      real (kind=R8) :: h(nfl+1), hh(nfl+1)

C       d:	sum  of n      	                number density
C       z2n:	sum  of n*Z**2 	                d * <Z**2>
C       upm:	(sum of n*up*Z**2 )/z2n         mean parallel velocity
C       p:	p(i) = d(i)*t(i) 	        pressure
C       t:	                                temperature
C       gt:	                                temperature gradient

      real (kind=R8) ::
     *     z2n(ns), d(ns), upm(ns), p(ns)
      real (kind=R8) ::
     *     t(ns), gt(ns)

      real (kind=R8) ::
     *          g(ns,ns,9), tt(ns,ns)
      real (kind=R8) ::
     *     mu(ns,ns), mu2(ns,ns)
      real (kind=R8) ::
     *          s2, s3, s4, s6a, s6b, s8, ss, cs, cond, condval

      real (kind=R8) ::
     *          a(2*ns,2*ns), x(2*ns), b(2*ns)
      real (kind=R8) ::
     *          inva(2*ns,2*ns), wk1(2*ns)
#ifndef NO_NAG
#ifndef NO_F01AAF
#else
      integer   IPIV(2*ns)
#endif
#else
      integer   IPIV(2*ns)
#endif

      integer   n1(ns)
      integer   i, j, k, l
      real (kind=R8) ::
     *          i1, i2, h1, hh1, const, cii


      if( nz(1) .ne. 1 )  stop
     *    'thermf:  nz(1) != 1 --> electron var. not first element'

*=======================================================================
*//COMPUTATION//

C     n1(j): 	index of vector element corresponding to first charge state
C               of species j in vectors n, z, up, r
      n1(1)= 1
      do i=2,ns
        n1(i) = n1(i-1) + nz(i-1)
      end do

      t(1) = te
      gt(1) = gte

      do j=2,ns
        t(j)  = ti
        gt(j) = gti
      end do

      do j=1,ns
        d(j)   = n(n1(j))
        z2n(j) = n(n1(j))*z2(n1(j))
        upm(j) = n(n1(j))*up(n1(j))*z2(n1(j))
        do k=2, nz(j)
          d(j)   = d(j)   + n(n1(j)+k-1)
          z2n(j) = z2n(j) + n(n1(j)+k-1)*z2(n1(j)+k-1)
          upm(j) = upm(j) + n(n1(j)+k-1)*up(n1(j)+k-1)*z2(n1(j)+k-1)
        end do
        upm(j) = upm(j)/z2n(j)
        p(j)  = d(j)*t(j)
      end do

C     ----  const = sqrt(2)*e**4 / (12*pi**1.5_R8*eps0**2)
!xpb      const = 1.79e-55_R8
        const = sqrt(2.0_R8)*qe**4 / (12*pi**1.5_R8*eps0**2)		 !xpb
!xpb      lnlam = 12.e0_R8

      do j=1,ns
        do i=1,ns
          tt(i,j)  = (t(i)*mm(j)+t(j)*mm(i))/(mm(i)+mm(j))
          mu(i,j)  = mm(i)*mm(j)/(mm(i)+mm(j))
          mu2(i,j) = mu(i,j)/(mm(i)+mm(j))
        end do
      end do

      do j=1,ns
        do i=1,ns
          g(i,j,1) = -const * z2n(i)*z2n(j) * sqrt(mu(i,j)) /
     *     tt(i,j)**1.5_R8 * lnlam
          g(i,j,2) = -3.e0_R8/5.e0_R8 * g(i,j,1)
          g(i,j,3) =  3.e0_R8/14.e0_R8 * g(i,j,1)
          g(i,j,5) =  2.7_R8 * mu2(i,j) * g(i,j,1)
          g(i,j,4) =  (1.3_R8*mm(j)/mm(i) + 1.6_R8 + 3.e0_R8*mm(i)/
     1     mm(j))/2.7_R8
     *                                                    * g(i,j,5)
          g(i,j,6) = -2.e0_R8/9.e0_R8*
     1     ( 23.e0_R8/28.e0_R8*mm(j)/mm(i) +8.e0_R8/7.e0_R8 +
     *     3.e0_R8*mm(i)/mm(j) )*g(i,j,5)
          g(i,j,7) = -25.e0_R8/42.e0_R8 * g(i,j,5)
          g(i,j,9) = -35.e0_R8/6.e0_R8 * mu2(i,j) * g(i,j,7)
          g(i,j,8) =  8.e0_R8/75.e0_R8 *
     1     ( 433.e0_R8/280.e0_R8*(mm(j)/mm(i))**2 +
     *     136.e0_R8/35.e0_R8*mm(j)/mm(i) + 459.e0_R8/35.e0_R8 +
     *     32.e0_R8/5.e0_R8*mm(i)/mm(j) + 5.e0_R8*(mm(i)/mm(j))**2 )
     *                  * g(i,j,9)
        end do
      end do


c  ---------------------------------------------------------------------
c       linear system  A * x = b with 2*ns  Eqs.:
c       x(1,...,ns) is q/p, x(ns+1,...,2*ns) is qq/p
c  ---------------------------------------------------------------------

      do j=1,ns
        do i=1,ns
           a(i,j)       = - g(i,j,5)*tt(i,j)
           a(i+ns,j)    = - 7.e0_R8*mu(i,j)/mm(i)*tt(i,j)**2*g(i,j,7)
           a(i,j+ns)    = - g(i,j,7)*tt(i,j)/t(j)*mu(i,j)
           a(i+ns,j+ns) = - mm(j)/t(j)*g(i,j,9)*tt(i,j)**2
        end do
        do k=1,ns
           a(j,j)       = a(j,j)       + tt(j,k)*g(j,k,4)
           a(j+ns,j)    = a(j+ns,j)    + 7.e0_R8*mu(j,k)/mm(j)*
     *                                       tt(j,k)**2*g(j,k,6)
           a(j,j+ns)    = a(j,j+ns)    + mu(j,k)*tt(j,k)/t(j)*g(j,k,6)
           a(j+ns,j+ns) = a(j+ns,j+ns) + mm(j)/t(j)*g(j,k,8)*tt(j,k)**2
        end do
      end do

c  ---------------------------------------------------------------------
c       solving A*x=b by inverting A
c       inverse of A is necessary for calc. heat conduct. coeff.
c  ---------------------------------------------------------------------
      ifail = +1

#ifndef NO_NAG
#ifndef NO_F01AAF
      call f01aaf(a, 2*ns, 2*ns, inva, 2*ns, wk1, ifail)
#else
      call f07adf(2*ns, 2*ns, a, 2*ns, IPIV, ifail)
      if( ifail .eq. 0 )  then
        call f06qff('General', 2*ns, 2*ns, a, 2*ns, inva, 2*ns)
        call f07ajf(2*ns, inva, 2*ns, IPIV, wk1, 2*ns, ifail)
      endif
#endif
#else
      call sgetrf(2*ns, 2*ns, a, 2*ns, IPIV, ifail)
      if( ifail .eq. 0 )  then
        inva(1:2*ns,1:2*ns) = a(1:2*ns,1:2*ns)
        call sgetri(2*ns, inva, 2*ns, IPIV, wk1, 2*ns, ifail)
      endif
#endif
      cond=condval(a,2*ns,ifail)

      if( ifail .ne. 0 )   return
        
      do i=1,ns
        b(i) = p(i)*gt(i)
        do k=1,ns
          b(i) = b(i) + mu(i,k)/mm(i)*tt(i,k)*g(i,k,2)*(upm(k)-upm(i))
        end do
        b(i) = 2.5_R8*b(i)
        b(i+ns) = 0.e0_R8
        do k=1,ns
          b(i+ns) = b(i+ns) +
     *      (mu(i,k)/mm(i)*tt(i,k))**2 * g(i,k,3)*(upm(k)-upm(i))
        end do
        b(i+ns) = 17.5_R8*b(i+ns)
      end do
        
      do i=1,2*ns
        x(i) = inva(i,1)*b(1)
        do j=2,2*ns
          x(i) = x(i) + inva(i,j)*b(j)
        end do
      end do

c  ---------------------------------------------------------------------
c       approx. heat conduction coeff.
c       ( gte in b(1), gti in b(2)...b(ns) )
c  ---------------------------------------------------------------------
      ce = - inva(1,1)*2.5_R8*p(1)*p(1)
      ci = 0.e0_R8
      do i=2,ns
        cii = 0.e0_R8
        do j=2,ns
          cii = cii - inva(i,j)*p(j)
        end do
        ci = ci + cii*p(i)
      end do
      ci = ci*2.5_R8

c  ---------------------------------------------------------------------
c       h(i,z), hh(i,z)
c       h steht fuer q/n (q: Waermestrom), hh fuer qq/n
c  ---------------------------------------------------------------------
      h(1)  = x(1)*t(1)
      hh(1) = x(1+ns)*t(1)

      do i=2,ns
        s2  = 0.e0_R8
        do k=1,ns
          s2 = s2 + 2.5_R8*mu(i,k)*g(i,k,2)*tt(i,k)/mm(i)
        end do
        s3  = 0.e0_R8
        do k=1,ns
          s3 = s3 + 17.5_R8*(mu(i,k)/mm(i))**2*g(i,k,3)*tt(i,k)**2
        end do
        s4  = 0.e0_R8
        do k=1,ns
          s4 = s4 + g(i,k,4)*tt(i,k)
        end do
        s6a = 0.e0_R8
        s6b = 0.e0_R8
        do k=1,ns
          s6a = s6a + mu(i,k)/mm(i)*g(i,k,6)*tt(i,k)
          s6b = s6b + mu(i,k)/mm(i)*g(i,k,6)*tt(i,k)**2
        end do
        s8  = 0.e0_R8
        do k=1,ns
          s8 = s8 + g(i,k,8)*tt(i,k)**2
        end do
        ss = 1.e0_R8/(7.e0_R8*s6a*s6b - s4*s8)
        do k=0,nz(i)-1
          cs = ss*2.5_R8*s8*t(i)*( z2n(i)/z2(n1(i)+k)-d(i) )*t(i)
          ci = ci + cs*n(n1(i)+k)
          h(n1(i)+k) = t(i) * ( x(i)
     *                 + ss*(s3*s6a-s2*s8)*(upm(i)-up(n1(i)+k)) )
     *                 - cs*gt(i)
          hh(n1(i)+k) = t(i) * ( x(i+ns)
     *       + ss*t(i)/mm(i)*((7.0_R8*s2*s6b-s3*s4)*(upm(i)-up(n1(i)+k))
     *       + 17.5_R8*s6b*t(i)*(z2n(i)/z2(n1(i)+k)-d(i))*gt(i) ) )
        end do
      end do

c  ---------------------------------------------------------------------
c       thermal forces r(i)
c  ---------------------------------------------------------------------

      do i=1,ns
        do k=0,nz(i)-1
          i1  = n(n1(i)+k)*z2(n1(i)+k)/z2n(i)
          h1  = h(n1(i)+k)/mm(i)
          hh1 = hh(n1(i)+k)/mm(i)
          r(n1(i)+k) = 0.e0_R8
          do j=1,ns
            do l=0,nz(j)-1
              i2 = n(n1(j)+l)*z2(n1(j)+l)/z2n(j)
              r(n1(i)+k) = r(n1(i)+k)
     *                   + g(i,j,2)*i1*i2 * mu(i,j)/tt(i,j) *
     *                     ( h1 - h(n1(j)+l)/mm(j) )
     *                   + g(i,j,3)*i1*i2 * (mu(i,j)/tt(i,j))**2 *
     *                     ( hh1 - hh(n1(j)+l)/mm(j) )
              if( i.eq.1 .or. j.eq.1 )  then
                r(n1(i)+k) = r(n1(i)+k) + g(i,j,1)*i1*i2*
     *                                      (up(n1(i)+k)-up(n1(j)+l))
              end if
            end do
          end do
        end do
      end do
        
c  ---------------------------------------------------------------------
c       heat fluxes he, hi
c  ---------------------------------------------------------------------
      he = h(1)*n(1)
      hi = h(2)*n(2)
      do i=3,nfl
        hi = hi + h(i)*n(i)
      end do
        
      return
      end subroutine thermf
c
c NOTE from XPB, about the upcoming removal of F01AAF from the NAG library
c
c  Replace: CALL F01AAF(A,IA,N,X,IX,WKSPCE,IFAIL)
c  by     : CALL F07ADF(N,N,A,IA,IPIV,IFAIL)
c           CALL F06QFF('General',N,N,A,IA,X,IX)
c           CALL F07AJF(N,X,IX,IPIV,WKSPCE,LWORK,IFAIL)
c IPIV is an INTEGER array of length N, and
c LWORK is the dimension of WKSPCE (one must have LWORK .ge. MAX(1,N)).
c The routine just does a matrix inversion, and you will find all relevant
c information in section 4 of the preamble to chapter F01 in the NAG documentation.
c

      function condval(a,ns,ifail)
      use b2mod_types
      implicit none
      integer ns,ifail
      real (kind=R8) ::
     *   a(ns,ns),condval
      integer i,j
      real (kind=R8) ::
     *   a1,ainv1,asum,ainvsum

      if(ifail.ne.0) then
        condval=0
      else
        a1=0
        ainv1=0
        do i=1,ns
          asum=0
          ainvsum=0
          do j=1,ns
            asum=asum+abs(a(i,j))
            ainvsum=ainvsum+abs(a(i,j))
          enddo
          a1=max(a1,asum)
          ainv1=max(ainv1,ainvsum)
        enddo
        condval=1.0_R8/(a1*ainv1)
      endif
      return
      end function condval
