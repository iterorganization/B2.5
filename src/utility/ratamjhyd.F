*-----------------------------------------------------------------------
*.specification

      subroutine ratamjhyd (nrid, nz, am, rtt, nt, rtn, nn,
     &  rtsa, rtra, rtqa, rtqr, rtcx, rtrd, rtbr,
     &  fix_recomb, zmax_recomb)
      use b2mod_types
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nrid, nz, nt, nn, fix_recomb, zmax_recomb
      real (kind=R8) ::
     *  am, rtt(0:nt), rtn(0:nn)
*   ..input/output arguments
      real (kind=R8) ::
     *  rtsa(0:nt,0:nn,0:nz), rtra(0:nt,0:nn,0:nz),
     &  rtqa(0:nt,0:nn,0:nz), rtcx(0:nt,0:nn,0:nz),
     &  rtrd(0:nt,0:nn,0:nz), rtbr(0:nt,0:nn,0:nz),
     &  rtqr(0:nt,0:nn,0:nz)
      integer :: iend
*-----------------------------------------------------------------------
*.documentation
*     TODO: provide warning when outside valid range of EIRENE data
*
*  1. purpose
*
*     RATAMJHYD uses Amjuel-Hydhel data for the rate coefficients
*     for ionisation and charge; as function of temperature and density.
*     The call to ratamjhyd must be preceded by a call to one of the
*     comprehensive routines.
*
*     Maarten Blommaert - June 2017
*
*  3. description
*
*     This routine employs ...
*
*     rtt represents temperature or temperature per amu; units are eV.
*     rtn represents density; units are m**(-3).
*     rtsa(i,j,k) is ionisation rate coefficient; m**3/s.
*     rtra(i,j,k) is recombination rate coefficient; m**3/s.
*     rtqa(i,j,k) is energy loss rate coefficient; eV*m**3/s.
*     rtcx(i,j,k) is charge exchange rate coefficient; m**3/s.
*     rtrd(i,j,k) is line radiation rate coefficient; eV*m**3/s.
*     rtbr(i,j,k) is bremsstrahlung radiation rate coefficient; eV*m**3/s.
*     All at coordinates rtt(i) and rtn(j), and for the transition
*     starting from charge state k.
*     For rtsa, rtra, rtqa, rtrd, and rtbr, the quantities rtt and rtn are
*     interpreted as the electron temperature and density. For rtcx,
*     the quantity rtt is interpreted as the neutral hydrogen
*     temperature per atomic mass, and rtn is interpreted as the
*     electron density.
*
*     Amjuel-Hydhel is used for hydrogen data only, and only for
*     the processes of ionisation and charge-exchange for neutral hydrogen.
*     Other rate coefficients must be obtained by a preceding
*     call to one of the comprehensive packages, such as STRAHL, ADPAK,
*     or ADAS.
*     The call to the present routine will then cause the rate
*     coefficients for the selected processes to be overwritten.
*
*  4. references
*
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
*     am - real*8, input.
*     am specifies the nuclear mass of the species for which the rate
*     coefficients are to be evaluated, in amu. am must be positive.
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
*     rtsa - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtsa(i,j,k) specifies the
*     ionisation rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k+1. rtsa has
*     dimension m**3/s.
*     This subroutine will override rtsa(,,0).
*
*     rtra - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtra(i,j,k) specifies the
*     recombination rate coefficient at electron temperature rtt(i) and
*     electron density rtn(j) for the process k->k-1. rtra has
*     dimension m**3/s.
*     This subroutine will override rtra(,,nz).
*
*     rtqa - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtqa(i,j,k) specifies the
*     electron heat loss rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k, excluding the contributions coming from potential energy
*     of recombination and other recombination process contributions.
*     In short, it includes the processes not accounted for in Eirene.
*     rtqa has dimension eV*m**3/s.
*     This subroutine will override rtqa(,,nz).
*
*     rtqr - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtqr(i,j,k) specifies the
*     electron heat loss rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtqr only includes the contribution of recombination
*     which replaces AMJUEL when Eirene is not used. rtqr has
*     dimension eV*m**3/s.
*
*     rtcx - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtcx(i,j,k) specifies the charge
*     exchange rate coefficient at neutral hydrogen temperature rtt(i)
*     and electron density rtn(j) for the process k->k-1. rtcx has
*     dimension m**3/s.
*     This subroutine will override rtcx(,,nz).
*
*     rtrd - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtrd(i,j,k) specifies the
*     line radiation rate coefficient at electron temperature rtt(i)
*     and electron density rtn(j) for processes starting from charge
*     state k. rtrd has dimension eV*m**3/s.
*     This subroutine will leave rtrd(,,0) unchanged.
*
*     rtbr - (0:nt,0:nn,0:nz) real*8 array, input/output.
*     For (i,j,k) in (0:nt,0:nn,0:nz), rtbr(i,j,k) specifies the
*     bremsstrahlung radiation rate coefficient at electron temperature
*     rtt(i) and electron density rtn(j) emitted from charge
*     state k. rtbr has dimension eV*m**3/s.
*     This subroutine will overwrite rtbr(,,nz).
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      logical found
      integer n3, i, j
      real (kind=R8) :: rates(0:nt,0:nn)
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
      call subini ('ratamjhyd')
*   ..Get density and temperature bounds
      call ipgetr('b2ardr_t_min',tmin)
      call ipgetr('b2ardr_t_max',tmax)
      call ipgetr('b2ardr_nn_min',nmin)
      call ipgetr('b2ardr_nn_max',nmax)
*   ..Number of explicitly specified AMJUEL/HYDHEL reactions
      call ipgeti('b2ardr_nreac',nreac)
*   ..set n3
      n3 = (nz+1)*(nt+1)*(nn+1)
*   ..test nz, nt, nn
      call xertst (nz.eq.1, 'faulty argument nz')
      call xertst (0.le.nt.and.0.le.nn, 'faulty argument nt, nn')
      call xertst (am.gt.0.0_R8, 'faulty argument am')
*   ..test rtt, rtn
      call b2xvsg (nt+1, rtt, 1, 'rtt', '.gt.')
      call b2xvsg (nn+1, rtn, 1, 'rtn', '.gt.')
*   ..test rtsa, rtra, rtqa, rtcx, rtrd, rtbr
      call b2xvsg (n3, rtsa, 1, 'rtsa', '.ge.')
      call b2xvsg (n3, rtra, 1, 'rtra', '.ge.')
      call b2xvsg (n3, rtqa, 1, 'rtqa', '.ge.')
      call b2xvsg (n3, rtqr, 1, 'rtqr', '.ge.')
      call b2xvsg (n3, rtcx, 1, 'rtcx', '.ge.')
      call b2xvsg (n3, rtrd, 1, 'rtrd', '.ge.')
      call b2xvsg (n3, rtbr, 1, 'rtbr', '.ge.')
*   ..Format for reaction card reading
6666  FORMAT (A6,1X,A4,A9,A3)

*   .. Set up default reactions
cnh ionization
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.1.5    ','EI ',am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rtsa(:,:,0) = rates
      call get_amjhyd_rates('AMJUEL  ','H.10','2.1.5    ','EI ',am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rtqa(:,:,0) = rates
cnh charge exchange
      call get_amjhyd_rates('HYDHEL  ','H.1 ','3.1.8    ','CX ',am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rtcx(:,:,nz) = rates
cnh recombination
      call get_amjhyd_rates('AMJUEL  ','H.4 ','2.1.8    ','RC ',am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      rtra(:,:,nz) = rates
      call get_amjhyd_rates('AMJUEL  ','H.10','2.1.8    ','RC ',am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
cmb   In fact, potential energy contribution should still be added like in commented code.
c     This will, however lead to negative values in rtra which will give NaNs in the
c     logarithmic representation. Instead, one should add the potential energy contribution afterwards.
      if (fix_recomb.ne.0.and.nz.le.zmax_recomb) then
        rtqr(:,:,nz) = rates
      else
        rtqr(:,:,nz) = rates +
     1    depot*(am*(mp+me)/(me+am*mp))*rtra(:,:,nz)
        call xertst(any(rtqr(:,:,nz)<0.0_R8),
     1   'Negative rates present, try turning on b2ardr_fix_recomb')
      end if

c     Add Bremsstrahlung
      do i = 0,nt
        do j = 0,nn
          rtbr(i,j,nz) = bremsf * sqrt(rtt(i))
          rtqa(i,j,nz) = rtqa(i,j,nz) + bremsf * sqrt(rtt(i))
        enddo
      enddo

*   .. Possibility to overwrite default reactions in b2ar.dat
      if (nreac.gt.0) then
        found = .false.
        do while (.not.found)
          READ (nrid,'(A5)',end=99) row
          if (INDEX(row,'*reac')>0) then
            found = .true.
          endif
        end do
      endif
      do ireac = 1, nreac
        READ (nrid,'(A72)') row
        call read_token(row,' ',filnam,itok, ier,.false.)
        iend = itok
        call read_token(row(iend:),' ',h123,itok, ier,.false.)
        iend = iend + itok
        call read_token(row(iend:),' ',reac,itok, ier,.false.)
        iend = iend + itok
        call read_token(row(iend:),' ',crc,itok, ier,.false.)
        call get_amjhyd_rates(filnam,h123,reac,crc,am,
     1    nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
        select case (H123)
          case ('H.1')
            if (streql(crc,'CX')) then
              rtcx(:,:,nz) = rates
            else
              call xerrab('Not possible yet to specify a non-CX ' //
     1          'reaction of H.1')
            endif
          case ('H.4')
            if (streql(crc,'EI')) then
              rtsa(:,:,0) = rates
            elseif (streql(crc,'RC')) then
              rtra(:,:,nz) = rates
            elseif (streql(crc,'CX')) then
              rtcx(:,:,nz) = rates
            else
              call xerrab('Wrong value for CRC')
            endif
          case ('H.10')
            if (streql(crc,'EI')) then
              rtqa(:,:,0) = rates
            elseif (streql(crc,'RC')) then
              if (fix_recomb.ne.0.and.nz.le.zmax_recomb) then
                rtqr(:,:,nz) = rates
              else
                rtqr(:,:,nz) = rates +
     1            depot*(am*(mp+me)/(me+am*mp))*rtra(:,:,nz)
                call xertst(any(rtqr(:,:,nz)<0.0_R8),
     1            'Negative rates present, try turning on ' //
     2            'b2ardr_fix_recomb')
              end if
            else
              call xerrab('Not allowed combination of H.10 and CRC')
            endif
          case default
            write(*,*) 'Do not yet know about ',h123
            call xerrab('H-number on reaction card not yet implemented')
        end select
      enddo


* ..return
      call subend ()
      return

   99 call xerrab('No Amjuel/Hydhel reaction header found in b2ar.dat!')
      return

*----------------------------------------------------------------------
*.end ratamjhyd

      end subroutine ratamjhyd

      subroutine get_amjhyd_rates(filnam,h123,reac,crc,am,
     1 nt,nn,rtt,rtn,tmin,tmax,nmin,nmax,rates)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
c      use b2mod_b2cmpa
      IMPLICIT NONE
#     include <eir_reac.F>
c     Input/output variables
      integer nt,nn
      real (kind=R8) ::
     1 rtn(0:nn),rtt(0:nt),
     2 rates(0:nt,0:nn),D11(0:8,0:8)
      real (kind = R8) :: tmin, tmax, nmin, nmax
      real(kind=R8), intent(in) :: am
#ifdef DBG
      real (kind=R8) :: product
#endif

      CHARACTER :: FILNAM*8,H123*4,REAC*9,CRC*3

c     Local variables
      integer ir,i,j,k,it,jn,iu,iv,ii
      real (kind=R8) ::
     1 tef,def,tei,dej,dnnn,duder,
     2 dlog,vmin,vmax,Emin
      integer, parameter :: nv=100
      real (kind=R8) ::
     1 v(1:nv),vf(1:nv+1),t(1:2*nv),dt(1:2*nv),du(0:1),
     2 energy(1:2*nv,1:2*nv,1:2*nv), s(1:2*nv,1:2*nv,1:2*nv),
     3 sv(1:2*nv,1:2*nv,1:2*nv),smo(0:nt,0:1),
     4 eei(1:2*nv,1:2*nv,1:2*nv),eef(1:2*nv,1:2*nv,1:2*nv),
     5 mdr(1:2*nv),mpol(1:2*nv)
      real (kind=R8) ::
     *  rttt(0:nt), rtnt(0:nn)
      data duder/2e3_R8/

*  1. purpose
*     This routine gets the amjuel/hydhel rates based on the reaction card
*     information and returns it to B2.5 as a function of density and
*     temperature for usage in fluid neutral models.
*     The energy-dependence in charge-exchange rates is integrated away
*     by assuming a Maxwellian distribution for the neutrals and a
*     momentum-linearized CX-rate coefficient is extracted. It is
*     assumed that the neutral species index is = 0.
*  3. description
*     rtt represents temperature or temperature per amu; units are eV.
*     rtn represents density; units are m**(-3).
*  5. parameters
*     rtt - (0:nt) real*8 array, input.
*     rtt(0:nt) specifies the sequence of values of the ion temperature
*     per amu or electrion temperature (not per amu), expressed in eV,
*     for which the rate coefficients are evaluated. rtt(0:nt) must be
*     positive.
*
*     rttt(0:nt) specifies the sequence of values of the temperature
*     after truncation (and not per amu), expressed in eV, for which
*     the rate coefficients are
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
*     am - real*8, input.
*     am specifies the nuclear mass of the species for which the rate
*     coefficients are to be evaluated, in amu. am must be positive.
*
*     rates - (0:nt,0:nn) real*8 array, output.
*     nn specifies the upper dimension bound for the array rtn.
*     It will hold that 0.le.nn.

*   ..set cutlo
      if (.not.b2mod_math_initialised) call b2mod_math_init
*     ..Read Amjuel data
      IR=1
      CREAC = 0.0_R8
      CALL SLREAC_b2plot(IR,FILNAM,H123,REAC,CRC)
      select case (H123)
       case ('H.1')
        write(*,*) 'Calculate momentum-linearized CX-rate coefficient'
*   ..Truncate tables to user-defined bounds
        rttt = rtt*am !Conversion of ion temperature per atomic mass unit to ion temperature
        rtnt = rtn
        do ii = 0,nn
         if (rtn(ii)<nmin) then
           rtnt(ii) = nmin
           if (ii.EQ.0) then
            write(*,*) 'densities bounded below at ', nmin, ' m^-3'
           endif
         else if (rtn(ii)>nmax) then
           rtnt(ii) = nmax
           if (ii.EQ.nn) then
            write(*,*) 'densities bounded above at ', nmax, ' m^-3'
           endif
         endif
        enddo
        do ii = 0,nt
         if (rttt(ii)<tmin) then
           rttt(ii) = tmin
           if (ii.EQ.0) then
            write(*,*) 'temperatures bounded below at ', tmin, ' eV'
           endif
         else if (rttt(ii)>tmax) then
           rttt(ii) = tmax
           if (ii.EQ.nt) then
            write(*,*) 'temperatures bounded above at ', tmax, ' eV'
           endif
         endif
        enddo

        du(0)= -duder*0.5_R8
        du(1)= duder*0.5_R8
        do J=1,9
          D11(J-1,0)=CREAC(J,0,IR)
        enddo
*        -- Set up energy distribution
*    --First sample velocity distribution
        vmin = 0.2_R8*sqrt(rttt(0)*ev/mp/am)
        vmax = 5.0_R8*sqrt(rttt(nt)*ev/mp/am)
        dlog = (log(vmax)-log(vmin))/(nv-1)
        v = (/((exp(log(vmin)+iv*dlog)), iv= 0,nv-1)/)
        vf(1) = 0.0_R8
        vf(2:nv) = (0.5_R8*(v(1:nv-1)+v(2:nv)))
        vf(nv+1) = exp(log(vmax)+0.5_R8*dlog)
        t(1:nv) = -v(nv:1:-1)
        t(nv+1:2*nv) = v
        dt(1:nv) = vf(nv+1:2:-1)-vf(nv:1:-1)
        dt(nv+1:2*nv) = vf(2:nv+1)-vf(1:nv)
*    --Then calculate 3D energy distribution
        Emin = mp*am*0.5_R8*3*sqrt(tmin*ev/mp/am)**2
        do i = 1,2*nv
         do j = 1,2*nv
          do k = 1,2*nv
           energy(i,j,k) = mp*am*0.5_R8*(t(i)**2+t(j)**2+t(k)**2)
           if (energy(i,j,k)<Emin) then
             energy(i,j,k)=Emin
           end if
          enddo
         enddo
        enddo
        write(*,*) 'Extract energy-dependent cx rates'
*    -- Extract energy dependent rates
        EEF=LOG(1.0_R8/am*energy/ev)
        s=0.0_R8
        do I=0,8
         EEI=EEF**I
         s=s+D11(I,0)*EEI
        enddo
*    --Calculate CX-rate in m**2
        s=exp(s)*1.0d-4
*    --pre-calculate product of rate and velocity (slow)
        sv = s*sqrt(2*energy/mp/am)
        smo = 0.0_R8
        do iu = 0,1
         write(*,*) 'Calculate momentum source for drift ',
     1    'velocity index', iu,'of 1'
         do it=0,nt
*    --Assemble Maxwellian distributions and integrate
c       write(*,*) 'Assemble Maxwellians for temperature index', it,
c     1    'of', nt
          mdr = maxwell1d(t,mp*am,2*rttt(it)*ev,0.0_R8,nv)
#ifdef DBG
          do i = 1,2*nv
            if (mdr(i).lt.tiny(1.0_R8)) mdr(i) = 0.0_R8
          enddo
#endif
          mdr = mdr/dot_product(mdr,dt)
          mpol = maxwell1d(t,mp*am,2*rttt(it)*ev,du(iu),nv)
          mpol = mpol/dot_product(mpol,dt)
c         write(*,*) 'Integrate momentum source'
          do i = 1,2*nv
           do j = 1,2*nv
            do k = 1,2*nv
             if (mdr(j).ne.0.0_R8.and.mdr(k).ne.0.0_R8.and.
     &          mpol(i).ne.0.0_R8) then
#ifndef DBG
               smo(it,iu) = smo(it,iu)+sv(i,j,k)*mdr(j)*dt(j)*
     1          mdr(k)*dt(k)*t(i)*mpol(i)*dt(i)
#else
               if (mdr(j).ge.tiny(1.0_R8).and.
     &             mdr(k).ge.tiny(1.0_R8).and.
     &            mpol(i).ge.tiny(1.0_R8)) then
                 if (t(i).gt.tiny(1.0_R8)) then
                   product = log(sv(i,j,k)) + log(mdr(j)) + log(dt(j))
     &               + log(mdr(k)) + log(dt(k)) + log(t(i))
     &               + log(mpol(i)) + log(dt(i))
                   if (product.ge.cutll)
     >              smo(it,iu) = smo(it,iu) + exp(product)
                 else if (t(i).lt.-tiny(1.0_R8)) then
                   product = log(sv(i,j,k)) + log(mdr(j)) + log(dt(j))
     &               + log(mdr(k)) + log(dt(k)) + log(abs(t(i)))
     &               + log(mpol(i)) + log(dt(i))
                   if (product.ge.cutll)
     >              smo(it,iu) = smo(it,iu) - exp(product)
                 end if
               endif
#endif
             endif
            enddo
           enddo
          enddo
         enddo
        enddo
        write(*,*) 'linearize momentum source'
c       write(*,*) 'smo(:,0)'
c       write (*,*) (smo(it,0),it=0,nt)
c       write(*,*) 'smo(:,1)'
c       write (*,*) (smo(it,1),it=0,nt)
        rates(:,0) =(smo(:,1)-smo(:,0))/duder
c       write(*,*) 'rates(:,0)'
c       write (*,*) (rates(it,0),it=0,nt)
        rates = spread(rates(:,0),2,nn+1)

c       write(*,*) 'CX rate table'
c       do it = 0,nt
c         write (*,*) (rates(it,jn),jn=0,nn)
c       enddo
       case ('H.4','H.10')
*   ..Truncate tables to user-defined bounds
        rttt = rtt
        rtnt = rtn
        do ii = 0,nn
         if (rtn(ii)<nmin) then
           rtnt(ii) = nmin
           if (ii.EQ.0) then
            write(*,*) 'densities bounded below at ', nmin, ' m^-3'
           endif
         else if (rtn(ii)>nmax) then
           rtnt(ii) = nmax
           if (ii.EQ.nn) then
            write(*,*) 'densities bounded above at ', nmax, ' m^-3'
           endif
         endif
        enddo
        do ii = 0,nt
         if (rttt(ii)<tmin) then
           rttt(ii) = tmin
           if (ii.EQ.0) then
            write(*,*) 'temperatures bounded below at ', tmin, ' eV'
           endif
         else if (rttt(ii)>tmax) then
           rttt(ii) = tmax
           if (ii.EQ.nt) then
            write(*,*) 'temperatures bounded above at ', tmax, ' eV'
           endif
         endif
        enddo

        do J=1,9
         do I=1,9
           D11(J-1,I-1)=CREAC(J,I,IR)
         enddo
        enddo
        do it=0,nt
          do jn=0,nn
            TEF=LOG(rttt(it))
            DEF=LOG(rtnt(jn)*1.0e-8_R8*1.0e-6_R8)
            dnnn=0.0
            do I=0,8
              TEI=TEF**I
              do J=0,8
                DEJ=DEF**J
                DNNN=DNNN+D11(I,J)*TEI*DEJ
              enddo
            enddo
c         --Convert from cm**3/s to m**3/s
            rates(it,jn)=exp(dnnn)*1.0e-6_R8
          enddo
        enddo
       case default
        write(*,*) 'Do not yet know about ',h123
        call xerrab('H-number on reaction card not yet implemented')
       end select

      return

      contains

      function maxwell1d(v,m,T,u,nv)
      IMPLICIT NONE

*     Input variables
      integer nv
      real (kind=R8) ::
     1  m,T,u
      real (kind=R8) ::
     1  v(1:2*nv)

*     Output variables
      real (kind=R8) ::
     1  maxwell1d(1:2*nv)

*     Local variables
      integer i
#ifdef DBG
      real (kind=R8) :: product
#endif
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     maxwell1d evaluates a 1D drifting Maxwellian on the velocities v
*     given a particle mass m, temperature T, and drift velocity u.
*
*  5. parameters
*     t - real*8 , input.
*     t specifies the temperature, expressed in Joule, for which the
*     drifting Maxwellian has to be evaluated.
*
*     nv - integer, input.
*     nv specifies the number of discrete velocities with the same sign
*     used to evaluate the (drifting) Maxwellian distribution.
*     It will hold that 1.le.nt.
*
*     v - real*8 array, input
*     v specifies the velocity distribution for which the Maxwellian with
*     drift velocity u is to be evaluated. v hereby contains both the
*     positive and negative parts of the velocity distribution.
*
      do i = 1, 2*nv
#ifndef DBG
        maxwell1d(i) =
     .   sqrt(m/(2.0_R8*pi*T))*expu(-m*(v(i)-u)**2/(2.0_R8*T))
#else
        product = 0.5_R8*(log(m)-log(2.0_R8*pi*T))+
     .   (-m*(v(i)-u)**2/(2.0_R8*T))
        if (product.gt.cutll) then
          maxwell1d(i) = exp(product)
        else
          maxwell1d(i) = 0.0_R8
        end if
#endif
      end do
      return
      end function maxwell1d

      end subroutine get_amjhyd_rates

