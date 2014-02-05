c     =========================================================================
c     SUBROUTINES FOR COMPARING B2 WITH CORONAL EQUILIBRIUM DATA
c     =========================================================================

c     -------------------------------------------------------------------------
c     subroutine getsa:
c     calculating the ionisation and recombination coefficents rsa and rra for
c     impurity belonging to the same homonuclear sequence as species i, for 
c     all charge states, at temperature T in eV and ne=1e16/m3.
c     -------------------------------------------------------------------------

      subroutine getsa(ns,rcsa,rcra,T)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      implicit none
*   ..common blocks

      integer ns
      real (kind=R8) ::
     *   T
      real (kind=R8) ::
     *   rcra(0:ns-1), rcsa(0:ns-1), rlte, rlne,
     *   rclra(0:1,0:ns-1), rclsa(0:1,0:ns-1)
      real (kind=R8) ::
     *   cutlo, cutll, t0, t1, tt, wrk0, ne0, fxte, fxne
      real (kind=R8) ::
     *   machsfr, expu, trimg

      integer is,ite,ine
      logical ltst0
      external subini, subend
      external b2sqel
      intrinsic log, dim, exp, max, min
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
      trimg(t0,t1,tt) = max(t0,min(t1,tt))
*     (trimg will be used to trim the logarithmic gradients to avoid
*     some pathetic behaviour near the underflow limit.)

      call subini('getsa')
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))

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
      external subini, subend, xertst
      intrinsic abs

      call subini('callgetsa')
      call getsa(ns,rcsa,rcra,T)
      isstart=-1
      do is = 0, ns-1
        if (isstart.eq.-1.and.zn(is).eq.zn(i).and.
     &      abs(am(is)-am(i)).lt.1.0e-3_R8) then
          isstart=is
        endif
      enddo
      isend=-1
      do is = ns-1, isstart, -1
        if (isend.eq.-1.and.zn(is).eq.zn(i).and.
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

      subroutine setcorora
      use b2mod_types
      use b2mod_geo
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      implicit none
*   ..common blocks

      integer ix,iy,is

      real (kind=R8) ::
     *   ssum,gsum
      real (kind=R8) ::
     *   neco(-1:nx,-1:ny), rcra(-1:nx,-1:ny,0:ns-1),
     *   rcsa(-1:nx,-1:ny,0:ns-1), rcqa(-1:nx,-1:ny,0:ns-1),
     *   rcrd(-1:nx,-1:ny,0:ns-1), rcbr(-1:nx,-1:ny,0:ns-1),
     *   rcza(-1:nx,-1:ny,0:ns-1), rcz2(-1:nx,-1:ny,0:ns-1),
     *   rcpt(-1:nx,-1:ny,0:ns-1), rcpi(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   rclsa(-1:nx,-1:ny,0:1,0:ns-1), rclra(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclqa(-1:nx,-1:ny,0:1,0:ns-1), rclrd(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclbr(-1:nx,-1:ny,0:1,0:ns-1), rclza(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclz2(-1:nx,-1:ny,0:1,0:ns-1), rclpt(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclpi(-1:nx,-1:ny,0:1,0:ns-1), logT(-1:nx,-1:ny)

      external subini, subend
      external b2spel, b2sqel

      call subini('corora')
      neco(-1:nx,-1:ny) = 1.0e16_R8
*    ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns, ev, te, neco, rclsa, rclra,
     &     rclqa, rclrd, rclbr, rclza, rclz2, rclpt, rclpi)
*   ..compute electron rate coefficients
      call b2sqel (nx, ny, ns, ev, te, 
     &    rclsa, rclra, rclqa, rclrd, rclbr, rclza, rclz2, rclpt, rclpi,
     &    rcsa, rcra, rcqa, rcrd, rcbr, rcza, rcz2, rcpt, rcpi, logT)

c     initialise array on stack for coronal radiation
      if(iwrk .ge. wklng) call dropstackbottom()
      iwrk=iwrk+1
      ndim(1,iwrk)=nx
      ndim(2,iwrk)=ny
      ndim(3,iwrk)=ns
      graphlabel(iwrk)='Coronal radiation (W)'
      do is=0,ns-1
         textwk(is,iwrk)=textpl(is)
      enddo
      do is=0,ns-1
         do ix=-1,nx
            do iy=-1,ny
               work(ix,iy,is,iwrk)=(rcrd(ix,iy,is)+rcbr(ix,iy,is))
     .          *vol(ix,iy)*ne(ix,iy)*na(ix,iy,is)
            enddo
         enddo
      enddo

      gsum=0.0_R8
      do is=0,ns-1
         ssum=0.0_R8
         do iy=-1,ny
            do ix=-1,nx
               ssum=ssum+work(ix,iy,is,iwrk)
            enddo
         enddo
         write(*,*) 'Coronal radiation rate for is = ',is,' = ',ssum
         gsum=gsum+ssum
      enddo
      write(*,*) 'Coronal radiation rate for all species = ',gsum
      call subend()
      return
      end subroutine setcorora

c     -------------------------------------------------------------------------
c     subroutine setcoroni
c     calculates the array for coronal equilibrium density and puts it 
c     on top of the stack
c     b2plot command: coni
c     -------------------------------------------------------------------------

      subroutine setcoroni
      use b2mod_types
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmrc
      implicit none
*   ..common blocks

      integer ix,iy,is,is0,isstart,isend
      real (kind=R8) ::
     *   sum(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   n(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   prod(-1:nx,-1:ny),divid(-1:nx,-1:ny)
      real (kind=R8) ::
     *   ntot(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   neco(-1:nx,-1:ny), rcra(-1:nx,-1:ny,0:ns-1),
     *   rcsa(-1:nx,-1:ny,0:ns-1), rcqa(-1:nx,-1:ny,0:ns-1),
     *   rcrd(-1:nx,-1:ny,0:ns-1), rcbr(-1:nx,-1:ny,0:ns-1),
     *   rcza(-1:nx,-1:ny,0:ns-1), rcz2(-1:nx,-1:ny,0:ns-1),
     *   rcpt(-1:nx,-1:ny,0:ns-1), rcpi(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   rclsa(-1:nx,-1:ny,0:1,0:ns-1), rclra(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclqa(-1:nx,-1:ny,0:1,0:ns-1), rclrd(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclbr(-1:nx,-1:ny,0:1,0:ns-1), rclza(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclz2(-1:nx,-1:ny,0:1,0:ns-1), rclpt(-1:nx,-1:ny,0:1,0:ns-1),
     &   rclpi(-1:nx,-1:ny,0:1,0:ns-1)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny)
      integer co1,co2

      external subini, subend
      external b2spel, b2sqel

      call subini('coroni')
      neco(-1:nx,-1:ny) = 1.0e16_R8
*    ..compute log-log electron rate coefficients
      call b2spel (nx, ny, ns, ev, te, neco, 
     &    rclsa, rclra, rclqa, rclrd, rclbr, rclza, rclz2, rclpt, rclpi)
*   ..compute electron rate coefficients
      call b2sqel (nx, ny, ns, ev, te, 
     &    rclsa, rclra, rclqa, rclrd, rclbr, rclza, rclz2, rclpt, rclpi,
     &    rcsa, rcra, rcqa, rcrd, rcbr, rcza, rcz2, rcpt, rcpi, wrk0)

c     initialise array on stack for coronal radiation
      if(iwrk .ge. wklng) call dropstackbottom()
      iwrk=iwrk+1
      ndim(1,iwrk)=nx
      ndim(2,iwrk)=ny
      ndim(3,iwrk)=ns
      graphlabel(iwrk)='Coronal equilibrium ion density (m^-3)'
      do is=0,ns-1
         textwk(is,iwrk)=textpl(is)
      enddo
      do is=0,ns-1
         do ix=-1,nx
            do iy=-1,ny
               work(ix,iy,is,iwrk)=0.0_R8
            enddo
         enddo
      enddo

      is0 = 0
      do while (is0.lt.ns)
        is = is0
        ntot(-1:nx,-1:ny,is0) = 0.0_R8
        do while (zn(is).eq.zn(is0).and.
     &     abs(am(is)-am(is0)).lt.1.0e-3_R8.and.is.lt.ns)
          ntot(-1:nx,-1:ny,is0) = 
     .     ntot(-1:nx,-1:ny,is0) + na(-1:nx,-1:ny,is)
          is = is + 1
        enddo
        isstart = is0
        isend = is - 1
c     coronal equilibrium
        do co1=isstart,isend
          prod(-1:nx,-1:ny)=1.0_R8
          do co2=isstart,co1-1
            prod(-1:nx,-1:ny)=prod(-1:nx,-1:ny)*
     .       rcsa(-1:nx,-1:ny,co2)/rcra(-1:nx,-1:ny,co2+1)
          enddo
          sum(-1:nx,-1:ny,co1)=prod(-1:nx,-1:ny)
        enddo
        divid(-1:nx,-1:ny)=0.0_R8
        do co1=isstart,isend
          divid(-1:nx,-1:ny)=divid(-1:nx,-1:ny)+sum(-1:nx,-1:ny,co1)
        enddo
        n(-1:nx,-1:ny,isstart)=1.0_R8/divid(-1:nx,-1:ny)
        do co1=isstart+1,isend
          n(-1:nx,-1:ny,co1)=sum(-1:nx,-1:ny,co1)*n(-1:nx,-1:ny,isstart)
        enddo
        do co1=isstart,isend
          work(-1:nx,-1:ny,co1,iwrk)=
     .     n(-1:nx,-1:ny,co1)*ntot(-1:nx,-1:ny,isstart)
        enddo
c Move on to next homonuclear sequence
        is0 = is 
      enddo

      call subend()
      return
      end subroutine setcoroni
