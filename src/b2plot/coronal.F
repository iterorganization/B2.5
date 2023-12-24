c     =========================================================================
c     SUBROUTINES FOR COMPARING B2 WITH CORONAL EQUILIBRIUM DATA
c     =========================================================================

c     -------------------------------------------------------------------------
c     subroutine getsa:
c     calculating the ionisation and recombination coefficents rsa and rra for
c     impurity belonging to the same isonuclear sequence as species i, for
c     all charge states, at temperature T in eV and ne=1e16/m3.
c     -------------------------------------------------------------------------

      subroutine getsa(ns,rcsa,rcra,T)
      use b2mod_types
      use b2mod_math
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_subsys
      implicit none
*   ..common blocks

      integer ns
      real (kind=R8) ::
     *   T
      real (kind=R8) ::
     *   rcra(0:ns-1), rcsa(0:ns-1), rlte, rlne,
     *   rclra(0:1,0:ns-1), rclsa(0:1,0:ns-1)
      real (kind=R8) ::
     *   t0, t1, tt, wrk0, ne0, fxte, fxne

      integer is,ite,ine
      logical ltst0
      intrinsic log

      call subini('getsa')

      wrk0 = log(T)
      ne0 = log(1.0e16_R8)
*   ..compute rate coefficients
*   ..initialise indices for interpolation
      ite = 0
      ine = 0
*   ..set ltst0
      ltst0 = rtlt(0).le.wrk0
*   ..compute abscissae
      rlte = max(rtlt(0),min(rtlt(rtnt),wrk0))
      rlne = max(rtln(0),min(rtln(rtnn),ne0))
*   ..find ite
*     (linear search starting from previous ite, to establish
*      0.le.ite.lt.rtnt;
*      ite.eq.0.or.rtlt(ite).le.rlte;
*      ite.eq.rtnt-1.or.rlte.lt.rtlt(ite+1).)
    1 if (ite.lt.rtnt-1) then
        if (rtlt(ite+1).le.rlte) then
          ite = ite+1
          goto 1
        endif
      endif
    2 if (0.lt.ite) then
        if (rlte.lt.rtlt(ite)) then
          ite = ite-1
          goto 2
        endif
      endif
*   ..compute coefficient for linear interpolation
      fxte = (rlte-rtlt(ite))/(rtlt(ite+1)-rtlt(ite))
    3 if (ine.lt.rtnn-1) then
        if (rtln(ine+1).le.rlne) then
          ine = ine+1
          goto 3
        endif
      endif
    4 if (0.lt.ine) then
        if (rlne.lt.rtln(ine)) then
          ine = ine-1
          goto 4
        endif
      endif
*   ..compute coefficient for linear interpolation
      fxne = (rlne-rtln(ine))/(rtln(ine+1)-rtln(ine))
*   ..find ine
*     (linear search starting from previous ine, to establish
*       0.le.ine.lt.rtnn;
*       ine.eq.0.or.rtln(ine).le.rlne;
*       ine.eq.rtnn-1.or.rlne.lt.rtln(ine+1).)
*   ..compute log-log rate coefficients
      do is = 0, ns-1
*    ..compute coronal ionisation rate coefficient
          t0 = (1-fxne)*rtlsa(ite,ine,is)+fxne*rtlsa(ite,ine+1,is)
          t1 = (1-fxne)*rtlsa(ite+1,ine,is)+fxne*rtlsa(ite+1,ine+1,is)
          tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
          if (ltst0) then
            rclsa(1,is) = trimg(-20.0_R8,20.0_R8,tt)
          else
            rclsa(1,is) = trimg(2.0_R8,20.0_R8,tt)
          endif
          rclsa(0,is) = (1-fxte)*t0+fxte*t1-rclsa(1,is)*rlte
*    ..compute recombination rate coefficient
          t0 = (1-fxne)*rtlra(ite,ine,is)+fxne*rtlra(ite,ine+1,is)
          t1 = (1-fxne)*rtlra(ite+1,ine,is)+fxne*rtlra(ite+1,ine+1,is)
          tt = (t1-t0)/(rtlt(ite+1)-rtlt(ite))
          if (ltst0) then
            rclra(1,is) = trimg(-20.0_R8,20.0_R8,tt)
          else
            rclra(1,is) = trimg(-20.0_R8,20.0_R8,tt)
          endif
          rclra(0,is) = (1-fxte)*t0+fxte*t1-rclra(1,is)*rlte
*      ..compute ionisation rate coefficient
          rcsa(is) = 
     &     expu(rclsa(0,is)+rclsa(1,is)*wrk0)
*      ..compute recombination rate coefficient
          rcra(is) = 
     &     expu(rclra(0,is)+rclra(1,is)*wrk0)
      enddo

      call subend()
      return
      end subroutine getsa

c     -------------------------------------------------------------------------
c     subroutine callgetsa:
c     This subroutine is called by the b2plot command 
c     lt i coro
c     lt: (Te/eV)
c     i: nuclear charge of impurity
c     Output:
c     CORONAL:n(0) n(1) ... n(i)
c     -------------------------------------------------------------------------

      subroutine callgetsa(ns,T,i)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_subsys
      implicit none
*   ..common blocks

      integer ns,i
      integer is,isstart,isend

      real (kind=R8) ::
     *   sum(0:ns-1)
      real (kind=R8) ::
     *   n(0:ns-1)
      real (kind=R8) ::
     *   rcra(0:ns-1), rcsa(0:ns-1)
      real (kind=R8) ::
     *   T,prod,divid
      integer c1,c2
      real (kind=R8) ::
     *   result
      external xertst
      intrinsic abs

      call subini('callgetsa')
      call getsa(ns,rcsa,rcra,T)
      isstart=-1
      do is = 0, ns-1
        if (isstart.eq.-1.and.nint(zn(is)).eq.nint(zn(i)).and.
     &      abs(am(is)-am(i)).lt.1.0e-3_R8) then
          isstart=is
        endif
      enddo
      isend=-1
      do is = ns-1, isstart, -1
        if (isend.eq.-1.and.nint(zn(is)).eq.nint(zn(i)).and.
     &      abs(am(is)-am(i)).lt.1.0e-3_R8) then
          isend=is
        endif
      enddo    
      call xertst(0.le.isstart.and.isstart.lt.ns,
     .  'starting species index not found !')
      call xertst(isstart.le.isend.and.isend.lt.ns,
     .  'ending species index not found !')   

      do c1=isstart,isend
         prod=1.0_R8
         do c2=isstart,c1-1
            prod=prod*rcsa(c2)/rcra(c2+1)
         enddo
         sum(c1)=prod
      enddo
      divid=0
      do c1=isstart,isend
         divid=divid+sum(c1)
      enddo
      n(isstart)=1.0_R8/divid
      do c1=isstart+1,isend
         n(c1)=sum(c1)*n(isstart)
      enddo

      write(*,*)'CORONAL:'
      do c1=isstart,isend
         result=log10(n(c1))
         write(*,*)result
      enddo
      call subend()
      return
      end subroutine callgetsa

c     -------------------------------------------------------------------------
c     subroutine setcorora
c     calculates the array for coronal equilibrium radiation and puts it 
c     on top of the stack
c     b2plot command: cora
c     -------------------------------------------------------------------------

      subroutine setcorora(ismain)
      use b2mod_types
      use b2us_data
      use b2us_plasma
      use b2mod_diag
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_b2plot
      use b2mod_subsys
      implicit none
      integer ismain
*   ..common blocks

      integer iCv, iss
      real (kind=R8) ::
     *   ssum,gsum
      real (kind=R8) ::
     *   neco(mpg%nCv)
      type(B2Rates) :: rtco
      type(B2RatesWork) :: rtwco
      integer is
      external b2spel, b2sqel

      call subini('corora')
      neco = 1.0e16_R8
      call createB2Rates (mpg%nCv, state%ns, rtco)
      call createB2RatesWork (mpg%nCv, state%ns, rtwco)
*    ..compute log-log electron rate coefficients
      call b2spel (mpg%nCv, state%ns, ev, state%pl%te, neco, rtco)
*   ..compute electron rate coefficients
      call b2sqel (mpg%nCv, state%ns, ismain, switch, ev, state%pl%te,
     .  rtco, rtwco)
c     initialise array on stack for coronal radiation
      if(iwrk .ge. wklng) call dropstackbottom()
      iwrk=iwrk+1
      ndim(1,iwrk)=mpg%nCv
      ndim(2,iwrk)=nspecies
      graphlabel(iwrk)='Coronal radiation (W)'
      do iss=0,nspecies-1
        textwk(iss,iwrk)=textsp(iss)
      enddo
      work(:,:,iwrk)=0.0_R8
      do is=0,state%ns-1
        iss = b2espcr(is)-1
        work(:,iss,iwrk)=work(:,iss,iwrk)+
     .   (rtwco%rrd(:,is)+rtwco%rbr(:,is))
     .   *geo%cvVol(:)*state%dv%ne(:)*state%pl%na(:,is)
      enddo

      gsum=0.0_R8
      do iss=0,nspecies-1
        ssum=0.0_R8
        do iCv=1,mpg%nCv
          ssum=ssum+work(iCv,iss,iwrk)
        enddo
        write(*,*) 'Coronal radiation rate for iss = ',iss,' = ',ssum
        gsum=gsum+ssum
      enddo
      write(*,*) 'Coronal radiation rate for all species = ',gsum
      call destroyB2RatesWork( rtwco )
      call destroyB2Rates( rtco )
      call subend()
      return
      end subroutine setcorora

c     -------------------------------------------------------------------------
c     subroutine setcoroni
c     calculates the array for coronal equilibrium density and puts it 
c     on top of the stack
c     b2plot command: coni
c     -------------------------------------------------------------------------

      subroutine setcoroni(ismain)
      use b2mod_types
      use b2us_data
      use b2us_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      use b2mod_b2plot
      use b2mod_subsys
      implicit none
      integer ismain
*   ..common blocks

      integer is
      type(B2Rates) :: rtco
      type(B2RatesWork) :: rtwco
      real (kind=R8) ::
     *   neco(mpg%nCv)
      real (kind=R8) ::
     *   prod(mpg%nCv),divid(mpg%nCv)
      real (kind=R8) ::
     *   ntot(mpg%nCv,0:state%ns-1)
      real (kind=R8) ::
     *   n(mpg%nCv,0:state%ns-1)
      real (kind=R8) ::
     *   sum(mpg%nCv,0:state%ns-1)
      integer is0,isstart,isend,co1,co2
      external b2spel, b2sqel

      call subini('coroni')
      neco = 1.0e16_R8
      call createB2Rates (mpg%nCv, state%ns, rtco)
      call createB2RatesWork (mpg%nCv, state%ns, rtwco)
*    ..compute log-log electron rate coefficients
      call b2spel (mpg%nCv, state%ns, ev, state%pl%te, neco, rtco)
*   ..compute electron rate coefficients
      call b2sqel (mpg%nCv, state%ns, ismain, switch, ev, state%pl%te,
     .  rtco, rtwco)
c     initialise array on stack for coronal radiation
      if(iwrk .ge. wklng) call dropstackbottom()
      iwrk=iwrk+1
      ndim(1,iwrk)=mpg%nCv
      ndim(2,iwrk)=state%ns
      call customise_label('Coronal equilibrium ion density (m^{-3})',
     . graphlabel(iwrk))
      do is=0,state%ns-1
         textwk(is,iwrk)=textpl(is)
      enddo
      work(:,is,iwrk)=0.0_R8
      is0 = 0
      do while (is0.lt.state%ns)
        is = is0
        ntot(:,is0) = 0.0_R8
        do while (nint(zn(is)).eq.nint(zn(is0)).and.
     &     abs(am(is)-am(is0)).lt.1.0e-3_R8.and.is.lt.state%ns)
          ntot(:,is0) = ntot(:,is0) + state%pl%na(:,is)
          is = is + 1
        enddo
        isstart = is0
        isend = is - 1
c     coronal equilibrium
        do co1=isstart,isend
          prod=1.0_R8
          do co2=isstart,co1-1
            prod(:)=prod(:)*rtwco%rsa(:,co2)/rtwco%rra(:,co2+1)
          enddo
          sum(:,co1)=prod(:)
        enddo
        divid=0.0_R8
        do co1=isstart,isend
          divid(:)=divid(:)+sum(:,co1)
        enddo
        n(:,isstart)=1.0_R8/divid(:)
        do co1=isstart+1,isend
          n(:,co1)=sum(:,co1)*n(:,isstart)
        enddo
        do co1=isstart,isend
          work(:,co1,iwrk)=n(:,co1)*ntot(:,isstart)
        enddo
c Move on to next isonuclear sequence
        is0 = is 
      enddo

      call destroyB2RatesWork( rtwco )
      call destroyB2Rates( rtco )
      call subend()
      return
      end subroutine setcoroni
