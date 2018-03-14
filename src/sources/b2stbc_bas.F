*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_bas (nx, ny, ns, ismain,                         !xpb
     1  facdrift, vol, hx, hy, qz, qc, gs,                               !xpb
     1  pbs, bb, na, ua, uadia, vaecrb, te, ti, po, ne, ni,
     1  fna, fhe, fhi, fch, fne, fni, chce, chve, chci, chvi, calf,
     &  sna0, smo0, she0, shi0, sch0, sne0,
     &  wrong_flow, main_call, diagno,     !xpb
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage,
     &  b2stbc_fchy_dia_coreonly, secmodel, lbndusr, 
     &  area_fix, sheath_drift_fix)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_anomalous_transport
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, ismain, diagno
      real (kind=R8) ::
     *  facdrift(-1:nx,-1:ny), 
     &  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),    !xpb
     &  qz(-1:nx,-1:ny,0:1), pbs(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2),
     &  qc(-1:nx,-1:ny,0:1), bb(-1:nx,-1:ny,0:3)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1), vaecrb(-1:nx,-1:ny,0:1,0:ns-1), 
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1),
     &  fne(-1:nx,-1:ny,0:1,0:1), fni(-1:nx,-1:ny,0:1,0:1)
      real (kind=R8) ::
     1  chce(-1:nx,-1:ny,0:1,0:1), chve(-1:nx,-1:ny,0:1,0:1),
     2  chci(-1:nx,-1:ny,0:1,0:1), chvi(-1:nx,-1:ny,0:1,0:1),
     3  calf(-1:nx,-1:ny,0:1,0:1)
      real (kind=R8) ::
     &  b2stbc_fchy_dia, b2stbc_neoclassical, praverage(2)
      logical main_call, lbndusr                                         !xpb
      integer b2stbc_fchy_dia_coreonly, secmodel, area_fix,
     &  sheath_drift_fix
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1)
      logical wrong_flow                                                 !xpb
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_BAS provides linearised expressions for the sources of
*     density, parallel momentum, electron heat, atom heat, and charge,
*     due to the boundary conditions.
*     Here the boundary conditions are expressed by use of the B2.5
*     format.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ireg
      real (kind=R8) ::
     *  cutll, cutlo, tpr, tro, csb, cst, vte, s0, s1, fche, fchi, t0,
     1  csbm, cstm, fchy_dia, seec,
     2  naf(0:ns-1), udf(0:ns-1), nef, nif, tef, tif, pof, wdiaf,
     3  she_test, shi_test, zaf(0:ns-1), zxf, nxf, taf,
     4  prt(-1:nx,-1:ny)
      save ncall
*   ..procedures
      intrinsic min, max, abs, sqrt, exp, log, dim
      real (kind=R8) ::
     *  machsfr, expu2, hy1
      external subini, subend, xertst, machsfr, eseec0
      external b2xvsg, b2xvff, b2xvps
      expu2(t0) = dim(exp(min(20.0_R8,max(t0,cutll))),cutlo)
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
*     (expu2 is the exponential function with underflow set to 0 and
*     large values limited at a somewhat arbitrary level.)
*   ..initialisation
      data ncall/0/
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_bas')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     . 'invalid main plasma species index ismain')
      call xertst (.not.is_neutral(ismain),
     . 'invalid main plasma species ismain; must not be neutral')
*   ..test facdrift (not on first call)                                  !xpb
      if (ncall.gt.0)                                                    !xpb
     > call xertst (0.0_R8.le.minval(facdrift) .and.
     &       1.0_R8.ge.maxval(facdrift), 'faulty argument facdrift')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, hx, hy
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2, hy, 1, 'hy', '.gt.')
*    ..test state
       call b2xvps (nx, ny, ns, 
     &  ne, na, ua, te, ti, po, fna, fhe, fhi, fch)
*    ..test edge values of fne, fni, chce, chve, chci, chvi, calf
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
       call b2xvff (nx, ny, chce, 'chce')
       call b2xvff (nx, ny, chve, 'chve')
       call b2xvff (nx, ny, chci, 'chci')
       call b2xvff (nx, ny, chvi, 'chvi')
       call b2xvff (nx, ny, calf, 'calf')
*    ..test sign of chce, chci
       call b2xvsg (n2*2, chce(-1,-1,0,0), 1, 'chce', '.ge.')
       call b2xvsg (n2*2, chce(-1,-1,1,1), 1, 'chce', '.ge.')
       call b2xvsg (n2*2, chci(-1,-1,0,0), 1, 'chci', '.ge.')
       call b2xvsg (n2*2, chci(-1,-1,1,1), 1, 'chci', '.ge.')
      endif
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
*   ..initialise resignore array
      if (ncall.eq.0) then
        resignore(-1,-1:ny,1:2)=0
        resignore(nx,-1:ny,1:2)=0
        resignore(-1:nx,-1,1:2)=0
        resignore(-1:nx,ny,1:2)=0
      endif
*  ..compute total pressure
      prt(:,:) = ne(:,:)*te(:,:) + ni(:,:,0)*ti(:,:)
      do is = 0, ns_ext-1
        prt(:,:) = prt(:,:) + na_ext(:,:,is)*ta_ext(:,:,is) 
      end do

      if(ncall.eq.0) write(*,*) 'Doing B2.5-style bc'

* ..compute conditions at south boundary
*   ..initialise ireg
        ireg = cbirso
*   ..loop over boundary cells
      do ix = -1, nx
       if(region(ix,-1,0).eq.0) goto 101         ! do not do dead cells
*    ..update ireg
    1  if (cbrbrk(ireg)*nx.lt.ix+0.5_R8) then
        call xertst (ireg.lt.cbirso+cbnrso-1,
     &    'error in computing region on south boundary')
        ireg = ireg+1
        goto 1
       endif
*    ..compute s0, s1
       s0 = gs(topix(ix,-1),topiy(ix,-1),1)/qz(ix,-1,1)
       s1 = abs(pbs(topix(ix,-1),topiy(ix,-1),1))
*    ..compute cell-faced quantities
       do is = 0, ns-1
         naf(is) = (na(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &              na(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &             (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
         zaf(is) = (rza(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &              rza(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &             (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       enddo
       tef = (te(ix,-1)*hy1(topix(ix,-1),topiy(ix,-1))+
     &        te(topix(ix,-1),topiy(ix,-1))*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       tif = (ti(ix,-1)*hy1(topix(ix,-1),topiy(ix,-1))+
     &        ti(topix(ix,-1),topiy(ix,-1))*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       nef = (ne(ix,-1)*hy1(topix(ix,-1),topiy(ix,-1))+
     &        ne(topix(ix,-1),topiy(ix,-1))*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       nif = (ni(ix,-1,0)*hy1(topix(ix,-1),topiy(ix,-1))+
     &        ni(topix(ix,-1),topiy(ix,-1),0)*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       pof = (po(ix,-1)*hy1(topix(ix,-1),topiy(ix,-1))+
     &        po(topix(ix,-1),topiy(ix,-1))*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
       wdiaf=(prt(ix,-1)*bb(ix,-1,2)*
     &        hy1(topix(ix,-1),topiy(ix,-1))+
     &        prt(topix(ix,-1),topiy(ix,-1))*
     &        bb(topix(ix,-1),topiy(ix,-1),2)*hy1(ix,-1))/
     &       (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
*    ..compute cst, vte
       tpr = 0.0_R8
       tro = 0.0_R8
       do is = 0, ns-1
         tpr = tpr+naf(is)*(zaf(is)*tef+tif)
         tro = tro+am(is)*mp*naf(is)
       enddo
       do is = 0, ns_ext-1
         nxf = (na_ext(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &          na_ext(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
         taf = (ta_ext(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &          ta_ext(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
         zxf = (za_ext(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &          za_ext(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
         tpr = tpr+nxf*(zxf*tef+taf)
         tro = tro+am_ext(is)*mp*nxf
       enddo
       cst = sqrt(tpr/tro)
       vte = sqrt(tef/me)
*    ..compute Bohm fchi, fche
       t0 = 0.0_R8
       do is = 0, ns-1
        t0 = t0-zaf(is)*fna(topix(ix,-1),topiy(ix,-1),1,1,is)
       enddo
       do is = 0, ns_ext-1
        zxf = (za_ext(ix,-1,is)*hy1(topix(ix,-1),topiy(ix,-1))+
     &         za_ext(topix(ix,-1),topiy(ix,-1),is)*hy1(ix,-1))/
     &        (hy1(topix(ix,-1),topiy(ix,-1))+hy1(ix,-1))
        t0 = t0-zxf*fa_ext(topix(ix,-1),topiy(ix,-1),1,1,is)
       enddo
       fchi = qe*t0
!       if(po(ix,-1).ge.0.0_R8) then                                     !xpb
         fche = 
     .    sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1*expu2(-qe*pof/tef)
!       else                                                             !xpb
!         fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1                  !xpb
!       endif                                                            !xpb
       fchy_dia = 0.0_R8                                                 !sv 01.07.99
       if(maxval(facdrift) .ne. 0.0_R8) then                             !xpb
         if (.not.on_closed_surface(ix,-1)) then                         !sv
           if(periodic_bc.ne.1 .or. rightix(ix,-1).eq.ix+1) then         !xpb
             fchy_dia = -wdiaf                                           !xpb
     &        *(vol(topix(ix,-1),topiy(ix,-1))+vol(ix,-1))
     &        /(gs(topix(ix,-1),topiy(ix,-1),2)+gs(ix,-1,2))             !sv 01.07.99
     &        *(OnedBsqc(rightix(topix(ix,-1),topiy(ix,-1)),
     &                   rightiy(topix(ix,-1),topiy(ix,-1)))
     &         -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                     !sv
     &        *facdrift(ix,-1)*b2stbc_fchy_dia                           !xpb
           else                                                          !xpb
             fchy_dia = -wdiaf                                           !xpb
     &        *(vol(topix(ix,-1),topiy(ix,-1))+vol(ix,-1))
     &        /(gs(topix(ix,-1),topiy(ix,-1),2)+gs(ix,-1,2))             !sv 01.07.99
     &        *(1.0_R8/wbbv(topix(ix,-1),topiy(ix,-1),3)**2
     &         -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                     !xpb
     &        *facdrift(ix,-1)*b2stbc_fchy_dia                           !xpb
           endif                                                         !xpb
         elseif (on_closed_surface(ix,-1)) then                          !xpb
           if(periodic_bc.ne.1 .or. rightix(ix,-1).eq.ix+1) then         !xpb
 !              fchy_dia = -praverage(1)*wbbv(topix(ix,-1),topiy(ix,-1),2)       !xpb
 !      &        *(vol(topix(ix,-1),topiy(ix,-1))+vol(ix,-1))
 !      &        /(gs(topix(ix,-1),topiy(ix,-1),2)+gs(ix,-1,2))          !sv 01.07.99
 !      &        *(OnedBsqc(rightix(topix(ix,-1),topiy(ix,-1)),
 !      &                   rightiy(topix(ix,-1),topiy(ix,-1)))
 !      &         -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                  !sv
 !      &        *facdrift(ix,-1)*b2stbc_fchy_dia                        !xpb
              fchy_dia = -praverage(1)*facdrift(ix,-1)*b2stbc_fchy_dia   !sv 01.07.99
     &         *(OnedBsqc(rightix(topix(ix,-1),topiy(ix,-1)),
     &                    rightiy(topix(ix,-1),topiy(ix,-1)))
     &          -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                    !sv
           else                                                          !xpb
 !              fchy_dia = -praverage(1)*wbbv(topix(ix,-1),topiy(ix,-1),2)       !xpb
 !      &        *(vol(topix(ix,-1),topiy(ix,-1))+vol(ix,-1))
 !      &        /(gs(topix(ix,-1),topiy(ix,-1),2)+gs(ix,-1,2))          !sv 01.07.99
 !      &        *(1.0_R8/wbbv(topix(ix,-1),topiy(ix,-1),3)**2
 !      &         -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                  !xpb
 !      &        *facdrift(ix,-1)*b2stbc_fchy_dia                        !xpb
              fchy_dia = -praverage(1)*facdrift(ix,-1)*b2stbc_fchy_dia   !sv 01.07.99
     &        *(1.0_R8/wbbv(topix(ix,-1),topiy(ix,-1),3)**2
     &         -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                     !xpb
          endif                                                          !xpb
         endif                                                           !sv
       elseif (b2stbc_neoclassical .ne. 0.0_R8 .and.                     !xpb
     &         on_closed_surface(ix,-1)) then                            !xpb
         if(periodic_bc.ne.1 .or. rightix(ix,-1).eq.ix+1) then           !xpb
           fchy_dia = -praverage(1)                                      !xpb
     &      *(OnedBsqc(rightix(topix(ix,-1),topiy(ix,-1)),
     &                 rightiy(topix(ix,-1),topiy(ix,-1)))
     &       -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                       !xpb
     &      *b2stbc_neoclassical                                         !xpb
         else                                                            !xpb
           fchy_dia = -praverage(1)
     &      *(1.0_R8/wbbv(topix(ix,-1),topiy(ix,-1),3)**2
     &       -OnedBsqc(topix(ix,-1),topiy(ix,-1)))                       !xpb
     &      *b2stbc_neoclassical                                         !xpb
         endif                                                           !xpb
       else                                                              !xpb
         fchy_dia = 0.0_R8                                               !xpb
       endif                                                             !sv
       if(lbndusr.or.(.not.on_closed_surface(ix,-1).and.
     &                b2stbc_fchy_dia_coreonly.eq.1)) fchy_dia = 0.0_R8
*    ..loop over species
       do is = 0, ns-1
*     ..compute csb
        csb = sqrt((zaf(is)*tef+tif)/(am(is)*mp))
*     ..compute particle source
*!!     cbsna(6,,) and cbsna(7,,) to be introduced with ua*pbs(,,1)
        sna0(ix,-1,0,is) = sna0(ix,-1,0,is)+
     &    cbsna(0,is,ireg)*s0
        sna0(ix,-1,1,is) = sna0(ix,-1,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csb*s0+
     &    cbsna(3,is,ireg)*csb*s1+
     &    cbsna(4,is,ireg)*cst*s0+
     &    cbsna(5,is,ireg)*cst*s1
*     ..compute parallel momentum source
        smo0(ix,-1,0,is) = smo0(ix,-1,0,is)+cbsmo(0,is,ireg)*s0
        smo0(ix,-1,1,is) = smo0(ix,-1,1,is)+cbsmo(7,is,ireg)*s0          !rfs
        smo0(ix,-1,2,is) = smo0(ix,-1,2,is)+
     &   cbsmo(1,is,ireg)*(-csb)*pbs(topix(ix,-1),topiy(ix,-1),1)+
     &   cbsmo(2,is,ireg)*(-cst)*pbs(topix(ix,-1),topiy(ix,-1),1)+
     &   cbsmo(5,is,ireg)*
     &   max(0.0_R8,fna(topix(ix,-1),topiy(ix,-1),1,1,is))/na(ix,-1,is)
        smo0(ix,-1,3,is) = smo0(ix,-1,3,is)+
     &   cbsmo(3,is,ireg)*s0+
     &   cbsmo(4,is,ireg)*s1+
     &   cbsmo(6,is,ireg)*
     &   max(0.0_R8,-fna(topix(ix,-1),topiy(ix,-1),1,1,is))/na(ix,-1,is)
*     ..compute atom heat source
        t0 = na(ix,-1,is)/ni(ix,-1,0)
        shi0(ix,-1,0) = shi0(ix,-1,0)+
     &    cbshi(0,is,ireg)*t0*s0+
     &    cbshi(6,is,ireg)*
     &    max(0.0_R8,fna(topix(ix,-1),topiy(ix,-1),1,1,is))
        shi0(ix,-1,1) = shi0(ix,-1,1)+
     &    cbshi(1,is,ireg)*t0*s0
        shi0(ix,-1,3) = shi0(ix,-1,3)+
     &    cbshi(2,is,ireg)*t0*csb*s0+
     &    cbshi(3,is,ireg)*t0*csb*s1+
     &    cbshi(4,is,ireg)*t0*cst*s0+
     &    cbshi(5,is,ireg)*t0*cst*s1+
     &    cbshi(7,is,ireg)*
     &    max(0.0_R8,-fna(topix(ix,-1),topiy(ix,-1),1,1,is))/ni(ix,-1,0)
       enddo
*    ..compute electron heat source
       she0(ix,-1,0) = she0(ix,-1,0)+
     &   cbshe(0,ireg)*s0+
     &   cbshe(6,ireg)*max(0.0_R8,fne(topix(ix,-1),topiy(ix,-1),1,1))
       she0(ix,-1,1) = she0(ix,-1,1)+
     &   cbshe(1,ireg)*s0
       she0(ix,-1,3) = she0(ix,-1,3)+
     &   cbshe(2,ireg)*vte*s0+
     &   cbshe(3,ireg)*vte*s1+
     &   cbshe(4,ireg)*cst*s0+
     &   cbshe(5,ireg)*cst*s1+
     &   cbshe(7,ireg)*
     &    max(0.0_R8,-fne(topix(ix,-1),topiy(ix,-1),1,1))/ne(ix,-1)
*    ..compute charge source
       sch0(ix,-1,2) = sch0(ix,-1,2)+
     &   cbsch(0,ireg)*qe*vte*s0+
     &   cbsch(2,ireg)*qe*vte*s1+
     &   cbsch(4,ireg)*qe*cst*s1
       sch0(ix,-1,3) = sch0(ix,-1,3)+
     &   cbsch(1,ireg)*qe**2*vte*s0/te(ix,-1)+
     &   cbsch(3,ireg)*qe**2*vte*s1/te(ix,-1)+
     &   cbsch(5,ireg)*qe**2*cst*s1/te(ix,-1)
       if (secmodel.eq.1 .and. ncall.lt.3) write (*,*)
     . 'Kin. sec. emiss. coeff. model not available for south boundary!'
       seec = 1.0_R8 - cbsch(7,ireg)
!       if(po(ix,-1).ge.0.0_R8) then                                     !xpb
         t0 = max(cbsch(6,ireg)*fchi,(1.0_R8-seec)*fche)*qe/tef
!       else                                                             !xpb
!         t0 = 0.0_R8                                                    !xpb
!       endif                                                            !xpb
       sch0(ix,-1,0) = sch0(ix,-1,0)+fchy_dia                            !sv 01.07.99
       sch0(ix,-1,0) = sch0(ix,-1,0)+
     &  (1.0_R8-seec)*fche-cbsch(6,ireg)*fchi+t0*po(ix,-1)
       sch0(ix,-1,1) = sch0(ix,-1,1)-t0
 101   continue
      enddo

* ..compute conditions at north boundary
*   ..initialise ireg
      ireg = cbirno
*   ..loop over boundary cells
      do ix = -1, nx
       if(region(ix,ny,0).eq.0) goto 202     ! do not do dead cells
*    ..update ireg
    2  if (cbrbrk(ireg)*nx.lt.ix+0.5_R8) then
        call xertst (ireg.lt.cbirno+cbnrno-1,
     &    'error in computing region on north boundary')
        ireg = ireg+1
        goto 2
       endif
*    ..compute s0, s1
       s0 = gs(ix,ny,1)/qz(ix,ny,1)
       s1 = abs(pbs(ix,ny,1))
*    ..compute cell-faced quantities
       do is = 0, ns-1
         naf(is) = (na(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &              na(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     &             (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
         zaf(is) = (rza(ix,ny,is)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &              rza(bottomix(ix,ny),bottomiy(ix,ny),is)*hy1(ix,ny))/
     &             (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       enddo
       tef = (te(ix,ny)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &        te(bottomix(ix,ny),bottomiy(ix,ny))*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       tif = (ti(ix,ny)*hy(bottomix(ix,ny),bottomiy(ix,ny))+
     &        ti(bottomix(ix,ny),bottomiy(ix,ny))*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       nef = (ne(ix,ny)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &        ne(bottomix(ix,ny),bottomiy(ix,ny))*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       nif = (ni(ix,ny,0)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &        ni(bottomix(ix,ny),bottomiy(ix,ny),0)*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       pof = (po(ix,ny)*hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &        po(bottomix(ix,ny),bottomiy(ix,ny))*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
       wdiaf=(prt(ix,ny)*bb(ix,ny,2)*
     &        hy1(bottomix(ix,ny),bottomiy(ix,ny))+
     &        prt(bottomix(ix,ny),bottomiy(ix,ny))*
     &        bb(bottomix(ix,ny),bottomiy(ix,ny),2)*hy1(ix,ny))/
     &       (hy1(bottomix(ix,ny),bottomiy(ix,ny))+hy1(ix,ny))
*    ..compute cst, vte
       tpr = 0.0_R8
       tro = 0.0_R8
       do is = 0, ns-1
         tpr = tpr+naf(is)*(zaf(is)*tef+tif)
         tro = tro+am(is)*mp*naf(is)
       enddo
       do is = 0, ns_ext-1
         nxf = (na_ext(ix,-1,is)*hy1(bottomix(ix,-1),bottomiy(ix,-1))+
     &          na_ext(bottomix(ix,-1),bottomiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(bottomix(ix,-1),bottomiy(ix,-1))+hy1(ix,-1))
         taf = (ta_ext(ix,-1,is)*hy1(bottomix(ix,-1),bottomiy(ix,-1))+
     &          ta_ext(bottomix(ix,-1),bottomiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(bottomix(ix,-1),bottomiy(ix,-1))+hy1(ix,-1))
         zxf = (za_ext(ix,-1,is)*hy1(bottomix(ix,-1),bottomiy(ix,-1))+
     &          za_ext(bottomix(ix,-1),bottomiy(ix,-1),is)*hy1(ix,-1))/
     &         (hy1(bottomix(ix,-1),bottomiy(ix,-1))+hy1(ix,-1))
         tpr = tpr+nxf*(zxf*tef+taf)
         tro = tro+am_ext(is)*mp*nxf
       enddo
       cst = sqrt(tpr/tro)
       vte = sqrt(tef/me)
*    ..compute Bohm fchi, fche
       t0 = 0.0_R8
       do is = 0, ns-1
        t0 = t0+zaf(is)*fna(ix,ny,1,1,is)
       enddo
       do is = 0, ns_ext-1
        zxf = (za_ext(ix,-1,is)*hy1(bottomix(ix,-1),bottomiy(ix,-1))+
     &         za_ext(bottomix(ix,-1),bottomiy(ix,-1),is)*hy1(ix,-1))/
     &        (hy1(bottomix(ix,-1),bottomiy(ix,-1))+hy1(ix,-1))
        t0 = t0+zxf*fa_ext(ix,ny,1,1,is)
       enddo
       fchi = qe*t0
!       if(po(ix,ny).ge.0.0_R8) then                                     !xpb
         fche = 
     .    sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1*expu2(-qe*pof/tef)
!       else                                                             !xpb
!         fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1                  !xpb
!       endif                                                            !xpb
       fchy_dia = 0.0_R8                                                 !sv 01.07.99
       if(maxval(facdrift) .ne. 0.0_R8) then                             !xpb
         if(region(ix,ny,0).ne.5) then                                   !xpb
           fchy_dia = -wdiaf                                             !sv
     &      *(vol(ix,ny)+vol(bottomix(ix,ny),bottomiy(ix,ny)))
     &      /(gs(ix,ny,2)+gs(bottomix(ix,ny),bottomiy(ix,ny),2))         !sv 01.07.99
     &      *(OnedBsqc(rightix(ix,ny),rightiy(ix,ny))-OnedBsqc(ix,ny))   !sv 01.07.99
     &      *facdrift(ix,ny)*b2stbc_fchy_dia                             !xpb
         else                                                            !xpb
 !            fchy_dia = -praverage(2)*wbbv(ix,ny,2)                     !xpb
 !      .      *(vol(bottomix(ix,ny),bottomiy(ix,ny))+vol(ix,ny))
 !      .      /(gs(bottomix(ix,ny),bottomiy(ix,ny),2)+gs(ix,ny,2))      !xpb
 !      .      *(OnedBsqc(rightix(bottomix(ix,ny),bottomiy(ix,ny)),
 !      .                 rightiy(bottomix(ix,ny),bottomiy(ix,ny)))      !xpb
 !      .       -OnedBsqc(bottomix(ix,ny),bottomiy(ix,ny)))
 !      .      *facdrift(ix,ny)*b2stbc_fchy_dia                          !xpb
            fchy_dia = -praverage(2)*facdrift(ix,ny)*b2stbc_fchy_dia     !sv 01.07.99
     &       *(OnedBsqc(rightix(ix,ny),rightiy(ix,ny))-OnedBsqc(ix,ny))  !sv
         endif                                                           !xpb
       elseif (b2stbc_neoclassical .ne. 0.0_R8 .and.                     !xpb
     &         on_closed_surface(ix,ny)) then                            !xpb
         fchy_dia = -praverage(2)
     .    *(OnedBsqc(rightix(ix,ny),rightiy(ix,ny))-OnedBsqc(ix,ny))
     .    *b2stbc_neoclassical                                           !xpb
       else                                                              !xpb
         fchy_dia = 0.0_R8                                               !sv
       endif                                                             !sv
       if(lbndusr.or.b2stbc_fchy_dia_coreonly.eq.1) fchy_dia = 0.0_R8
*    ..loop over species
       do is = 0, ns-1
*     ..compute csb
        csb = sqrt((zaf(is)*tef+tif)/(am(is)*mp))
*     ..compute particle source
*!!     cbsna(6,,) and cbsna(7,,) to be introduced with ua*pbs(,,1)
        sna0(ix,ny,0,is) = sna0(ix,ny,0,is)+
     &    cbsna(0,is,ireg)*s0
        sna0(ix,ny,1,is) = sna0(ix,ny,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csb*s0+
     &    cbsna(3,is,ireg)*csb*s1+
     &    cbsna(4,is,ireg)*cst*s0+
     &    cbsna(5,is,ireg)*cst*s1
*     ..compute parallel momentum source
        smo0(ix,ny,0,is) = smo0(ix,ny,0,is)+cbsmo(0,is,ireg)*s0
        smo0(ix,ny,1,is) = smo0(ix,ny,1,is)+cbsmo(7,is,ireg)*s0          !rfs
        smo0(ix,ny,2,is) = smo0(ix,ny,2,is)+
     &    cbsmo(1,is,ireg)*csb*pbs(ix,ny,1)+
     &    cbsmo(2,is,ireg)*cst*pbs(ix,ny,1)+
     &    cbsmo(5,is,ireg)*max(0.0_R8,-fna(ix,ny,1,1,is))/na(ix,ny,is)
        smo0(ix,ny,3,is) = smo0(ix,ny,3,is)+
     &    cbsmo(3,is,ireg)*s0+
     &    cbsmo(4,is,ireg)*s1+
     &    cbsmo(6,is,ireg)*max(0.0_R8,fna(ix,ny,1,1,is))/na(ix,ny,is)
*     ..compute atom heat source
        t0 = na(ix,ny,is)/ni(ix,ny,0)
        shi0(ix,ny,0) = shi0(ix,ny,0)+
     &    cbshi(0,is,ireg)*t0*s0+
     &    cbshi(6,is,ireg)*max(0.0_R8,-fna(ix,ny,1,1,is))
        shi0(ix,ny,1) = shi0(ix,ny,1)+
     &    cbshi(1,is,ireg)*t0*s0
        shi0(ix,ny,3) = shi0(ix,ny,3)+
     &    cbshi(2,is,ireg)*t0*csb*s0+
     &    cbshi(3,is,ireg)*t0*csb*s1+
     &    cbshi(4,is,ireg)*t0*cst*s0+
     &    cbshi(5,is,ireg)*t0*cst*s1+
     &    cbshi(7,is,ireg)*max(0.0_R8,fna(ix,ny,1,1,is))/ni(ix,ny,0)
       enddo
*    ..compute electron heat source
       she0(ix,ny,0) = she0(ix,ny,0)+
     &   cbshe(0,ireg)*s0+
     &   cbshe(6,ireg)*max(0.0_R8,-fne(ix,ny,1,1))
       she0(ix,ny,1) = she0(ix,ny,1)+
     &   cbshe(1,ireg)*s0
       she0(ix,ny,3) = she0(ix,ny,3)+
     &   cbshe(2,ireg)*vte*s0+
     &   cbshe(3,ireg)*vte*s1+
     &   cbshe(4,ireg)*cst*s0+
     &   cbshe(5,ireg)*cst*s1+
     &   cbshe(7,ireg)*max(0.0_R8,fne(ix,ny,1,1))/ne(ix,ny)
*    ..compute charge source
       sch0(ix,ny,2) = sch0(ix,ny,2)+
     &   cbsch(0,ireg)*qe*vte*s0+
     &   cbsch(2,ireg)*qe*vte*s1+
     &   cbsch(4,ireg)*qe*cst*s1
       sch0(ix,ny,3) = sch0(ix,ny,3)+
     &   cbsch(1,ireg)*qe**2*vte*s0/te(ix,ny)+
     &   cbsch(3,ireg)*qe**2*vte*s1/te(ix,ny)+
     &   cbsch(5,ireg)*qe**2*cst*s1/te(ix,ny)
       if (secmodel.eq.1 .and. ncall.lt.3) write (*,*)
     . 'Kin. sec. emiss. coeff. model not available for north boundary!'
       seec = 1.0_R8 - cbsch(7,ireg)
!       if(po(ix,ny).ge.0.0_R8) then                                     !xpb
         t0 = max(cbsch(6,ireg)*fchi,(1.0_R8-seec)*fche)*qe/tef
!       else                                                             !xpb
!         t0 = 0.0_R8                                                    !xpb
!       endif                                                            !xpb
       sch0(ix,ny,0) = sch0(ix,ny,0)-fchy_dia             !sv  minus 01.07.99
       sch0(ix,ny,0) = sch0(ix,ny,0)+
     &  (1.0_R8-seec)*fche-cbsch(6,ireg)*fchi+t0*po(ix,ny)
       sch0(ix,ny,1) = sch0(ix,ny,1)-t0
 202   continue
      enddo

* ..compute conditions at west boundary
*   ..initialise ireg
      ireg = cbirwe
*   ..loop over boundary cells
      do iy = -1, ny
       if(region(-1,iy,0).eq.0) goto 303    ! do not do dead cells
*    ..update ireg
    3  if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
        call xertst (ireg.lt.cbirwe+cbnrwe-1,
     &    'error in computing region on west boundary')
        ireg = ireg+1
        goto 3
       endif
*    ..compute s0, s1: hb
       s0 = gs(rightix(-1,iy),rightiy(-1,iy),0)
       if (area_fix.ge.3) s0 = s0*qc(rightix(-1,iy),rightiy(-1,iy),0)       !srv 09.03.09
       s1 = abs(pbs(rightix(-1,iy),rightiy(-1,iy),0))
*    ..compute cell-faced quantities
       do is = 0, ns-1
         naf(is) = (na(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &              na(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         zaf(is) = (rza(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &              rza(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &             (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
cxpb  Following suggestion from Russians, the drift velocity should
cxpb  only include ExB terms and not diamagnetic contributions
         if (sheath_drift_fix.eq.0) then
           udf(is) =
     &      (uadia(-1,iy,0,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &       uadia(rightix(-1,iy),rightiy(-1,iy),0,is)*hx(-1,iy))/
     &      (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         else
           udf(is) =
     &      (vaecrb(-1,iy,0,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &       vaecrb(rightix(-1,iy),rightiy(-1,iy),0,is)*hx(-1,iy))/
     &      (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         endif
       enddo
       tef = (te(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &        te(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &       (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
       tif = (ti(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &        ti(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &       (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
       nef = (ne(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &        ne(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &       (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
       nif = (ni(-1,iy,0)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &        ni(rightix(-1,iy),rightiy(-1,iy),0)*hx(-1,iy))/
     &       (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
       pof = (po(-1,iy)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &        po(rightix(-1,iy),rightiy(-1,iy))*hx(-1,iy))/
     &       (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
*    ..compute cst, vte
       tpr = 0.0_R8
       tro = 0.0_R8
       do is = 0, ns-1
         tpr = tpr+naf(is)*(zaf(is)*tef+tif)
         tro = tro+am(is)*mp*naf(is)
       enddo
       do is = 0, ns_ext-1
         nxf = (na_ext(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &          na_ext(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &         (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         zxf = (za_ext(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &          za_ext(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &         (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         taf = (ta_ext(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &          ta_ext(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &         (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
         tpr = tpr+nxf*(zxf*tef+taf)
         tro = tro+am_ext(is)*mp*nxf
       enddo
       cst = sqrt(tpr/tro)
       vte = sqrt(tef/me)
*    ..compute Bohm fchi, fche
       t0 = 0.0_R8
       do is = 0, ns-1
        t0 = t0-zaf(is)*fna(rightix(-1,iy),rightiy(-1,iy),0,0,is)
       enddo
       do is = 0, ns_ext-1
        zxf = (za_ext(-1,iy,is)*hx(rightix(-1,iy),rightiy(-1,iy))+
     &         za_ext(rightix(-1,iy),rightiy(-1,iy),is)*hx(-1,iy))/
     &        (hx(rightix(-1,iy),rightiy(-1,iy))+hx(-1,iy))
        t0 = t0-zxf*fa_ext(rightix(-1,iy),rightiy(-1,iy),0,0,is)
       enddo
       fchi = qe*t0
       if(t0.lt.0.0_R8 .and. region(-1,iy,0).ne.0 .and. 
     &  bottomiy(-1,iy).ne.-2 .and. topiy(-1,iy).ne.ny+1) then
         wrong_flow = .true.                                             !xpb
         if(diagno.ge.0) then
           write(*,*) 'b2stbc_west:iy,fchi',iy,fchi
           write(*,*) zaf(ismain),
     &      fna(rightix(-1,iy),rightiy(-1,iy),0,0,ismain)
           write(*,'(a,4g14.7)') 'Bx,Bz,Ud,-Ua =',
     &      wbbr(-1,iy,0),wbbr(-1,iy,2),udf(ismain),-ua(-1,iy,ismain)
         endif
c         fchi=0.0_R8
       endif
!       if (po(-1,iy).ge.0.0_R8) then                                     !sv
         fche = 
     .    sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1*expu2(-qe*pof/tef)       !sv
!       else                                                              !sv
!         fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1                   !sv
!       endif
*    ..loop over species
       do is = 0, ns-1
*     ..compute csb
        csb = sqrt((zaf(is)*tef+tif)/(am(is)*mp))
*     .. account for diamagnetic drifts
        csbm = csb+udf(is)*abs(wbbr(-1,iy,2)/wbbr(-1,iy,0))      !xpb
        cstm = cst+udf(is)*abs(wbbr(-1,iy,2)/wbbr(-1,iy,0))      !xpb
        if(.not.is_neutral(is)) then                                      !xpb
          csbm = max(0.0_R8,csbm)                                         !xpb
          cstm = max(0.0_R8,cstm)                                         !xpb
        endif                                                             !xpb
        if(csbm.gt.csb .and. diagno.ge.2)
     >   write(*,*) 'west:iy,is,csb,csbm',iy,is,csb,csbm
        if(cstm.gt.cst .and. diagno.ge.2)
     >   write(*,*) 'west:iy,is,cst,cstm',iy,is,cst,cstm
*     ..compute particle source hb
        sna0(-1,iy,0,is) = sna0(-1,iy,0,is)+
     &    cbsna(0,is,ireg)*s0+
     &    cbsna(6,is,ireg)*
     &    max(0.0_R8,ua(-1,iy,is)*pbs(rightix(-1,iy),rightiy(-1,iy),0))
        sna0(-1,iy,1,is) = sna0(-1,iy,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csbm*s0+
     &    cbsna(3,is,ireg)*csbm*s1+
     &    cbsna(4,is,ireg)*cstm*s0+
     &    cbsna(5,is,ireg)*cstm*s1+
     &    cbsna(7,is,ireg)*
     &    max(0.0_R8,-ua(-1,iy,is)*pbs(rightix(-1,iy),rightiy(-1,iy),0))
*     ..compute parallel momentum source
        smo0(-1,iy,0,is) = smo0(-1,iy,0,is)+cbsmo(0,is,ireg)*s0
        smo0(-1,iy,1,is) = smo0(-1,iy,1,is)+cbsmo(7,is,ireg)*s0           !rfs
        smo0(-1,iy,2,is) = smo0(-1,iy,2,is)+
     &    cbsmo(1,is,ireg)*(-csbm)*
     &    pbs(rightix(-1,iy),rightiy(-1,iy),0)+                           !xpb
     &    cbsmo(2,is,ireg)*(-cstm)*
     &    pbs(rightix(-1,iy),rightiy(-1,iy),0)+                           !xpb
     &    cbsmo(5,is,ireg)*max(0.0_R8,
     &     -fna(rightix(-1,iy),rightiy(-1,iy),0,0,is))/na(-1,iy,is)         !sv
        smo0(-1,iy,3,is) = smo0(-1,iy,3,is)+
     &    cbsmo(3,is,ireg)*s0+
     &    cbsmo(4,is,ireg)*s1+
     &    cbsmo(6,is,ireg)*max(0.0_R8,
     &     -fna(rightix(-1,iy),rightiy(-1,iy),0,0,is))/na(-1,iy,is)
*     ..compute atom heat source
        t0 = na(-1,iy,is)/ni(-1,iy,0)
        shi0(-1,iy,0) = shi0(-1,iy,0)+
     &    cbshi(0,is,ireg)*t0*s0+
     &    cbshi(6,is,ireg)*
     &     max(0.0_R8,fna(rightix(-1,iy),rightiy(-1,iy),0,0,is))
        shi0(-1,iy,1) = shi0(-1,iy,1)+
     &    cbshi(1,is,ireg)*t0*s0
        shi0(-1,iy,3) = shi0(-1,iy,3)+
     &    cbshi(2,is,ireg)*t0*csbm*s0+
     &    cbshi(3,is,ireg)*t0*csbm*s1+
     &    cbshi(4,is,ireg)*t0*cstm*s0+
     &    cbshi(5,is,ireg)*t0*cstm*s1+
     &    cbshi(7,is,ireg)*max(0.0_R8,
     &     -fna(rightix(-1,iy),rightiy(-1,iy),0,0,is))/ni(-1,iy,0)
       enddo
       shi_test=
     1  shi0(-1,iy,0)+
     2  shi0(-1,iy,1)*ti(-1,iy)+
     3  shi0(-1,iy,2)*ni(-1,iy,0)+
     4  shi0(-1,iy,3)*ni(-1,iy,0)*ti(-1,iy)
       if(shi_test.gt.0.0_R8.and.diagno.ge.2) then
         write(*,*) 'shi_test>0: W ',
     1    iy,shi0(-1,iy,:),ti(-1,iy)/qe,ni(-1,iy,0),s0,s1,csbm,cstm,
     2    ireg,cbshi(:,:,ireg)
       endif
*    ..compute electron heat source
       she0(-1,iy,0) = she0(-1,iy,0)+
     &  cbshe(0,ireg)*s0+
     &  cbshe(6,ireg)*max(0.0_R8,fne(rightix(-1,iy),rightiy(-1,iy),0,0))
       she0(-1,iy,1) = she0(-1,iy,1)+
     &   cbshe(1,ireg)*s0
       she0(-1,iy,3) = she0(-1,iy,3)+
     &   cbshe(2,ireg)*vte*s0+
     &   cbshe(3,ireg)*vte*s1+
     &   cbshe(4,ireg)*cstm*s0+
     &   cbshe(5,ireg)*cstm*s1+
     &   cbshe(7,ireg)*
     &    max(0.0_R8,-fne(rightix(-1,iy),rightiy(-1,iy),0,0))/ne(-1,iy)
       she_test=she0(-1,iy,0)+
     2  she0(-1,iy,1)*te(-1,iy)+
     3  she0(-1,iy,2)*ne(-1,iy)+
     4  she0(-1,iy,3)*ne(-1,iy)*te(-1,iy)
       if(she_test.gt.0.0_R8.and.diagno.ge.2) then
         write(*,*) 'she_test>0: W ',
     1    iy,she0(-1,iy,:),te(-1,iy)/qe,ne(-1,iy),s0,s1,vte,cstm,
     2    ireg,cbshe(:,ireg)
       endif
*    ..compute charge source
       sch0(-1,iy,2) = sch0(-1,iy,2)+
     &   cbsch(0,ireg)*qe*vte*s0+
     &   cbsch(2,ireg)*qe*vte*s1+
     &   cbsch(4,ireg)*qe*cstm*s1
       sch0(-1,iy,3) = sch0(-1,iy,3)+
     &   cbsch(1,ireg)*qe**2*vte*s0/te(-1,iy)+
     &   cbsch(3,ireg)*qe**2*vte*s1/te(-1,iy)+
     &   cbsch(5,ireg)*qe**2*cstm*s1/te(-1,iy)
       if (secmodel.eq.1) then
         if (ne(-1,iy)/(rza(-1,iy,ismain)*na(-1,iy,ismain)).gt.1.1_R8) 
     >    then
           write(*,*) 'Large ne contribution from impurities!'
           write(*,*) 'Kinetic s.e.c. model may not be valid!'
         endif
         call eseec0 (tef,tif,nef,pof,wbbr(-1,iy,0),wbbr(-1,iy,3),
     &    0,-1,iy,nx,ny,ns,ismain,seec,diagno)                          !cfs
       else
         seec = 1.0_R8 - cbsch(7,ireg)
       endif
!       if (po(-1,iy).ge.0.0_R8) then                                   !sv
         t0 = max(cbsch(6,ireg)*fchi,(1.0_R8-seec)*fche)*qe/tef
!       else                                                            !sv
!         t0 = 0.0_R8                                                   !sv
!       endif                                                           !sv
       sch0(-1,iy,0) = sch0(-1,iy,0)+
     &  (1.0_R8-seec)*fche-cbsch(6,ireg)*fchi+t0*po(-1,iy)
       sch0(-1,iy,1) = sch0(-1,iy,1)-t0
 303   continue
      enddo

* ..compute conditions at east boundary
*   ..initialise ireg
      ireg = cbirea
*   ..loop over boundary cells
      do iy = -1, ny
       if(region(nx,iy,0).eq.0) goto 404   ! do not do dead cells
*    ..update ireg
    4  if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
        call xertst (ireg.lt.cbirea+cbnrea-1,
     &    'error in computing region on east boundary')
        ireg = ireg+1
        goto 4
       endif
*    ..compute s0, s1
       s0 = gs(nx,iy,0)
       if (area_fix.ge.3) s0 = s0*qc(nx,iy,0)                       !srv 09.03.09
       s1 = abs(pbs(nx,iy,0))
*    ..compute cell-faced quantities
       do is = 0, ns-1
         naf(is) = (na(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &              na(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         zaf(is) = (rza(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &              rza(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &             (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
cxpb  Following suggestion from Russians, the drift velocity should
cxpb  only include ExB terms and not diamagnetic contributions
         if (sheath_drift_fix.eq.0) then
           udf(is) =
     &      (uadia(nx,iy,0,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &       uadia(leftix(nx,iy),leftiy(nx,iy),0,is)*hx(nx,iy))/
     &      (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         else
           udf(is) =
     &      (vaecrb(nx,iy,0,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &       vaecrb(leftix(nx,iy),leftiy(nx,iy),0,is)*hx(nx,iy))/
     &      (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         endif
       enddo
       tef = (te(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &        te(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &       (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
       tif = (ti(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &        ti(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &       (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
       nef = (ne(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &        ne(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &       (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
       nif = (ni(nx,iy,0)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &        ni(leftix(nx,iy),leftiy(nx,iy),0)*hx(nx,iy))/
     &       (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
       pof = (po(nx,iy)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &        po(leftix(nx,iy),leftiy(nx,iy))*hx(nx,iy))/
     &       (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
*    ..compute cst, vte
       tpr = 0.0_R8
       tro = 0.0_R8
       do is = 0, ns-1
         tpr = tpr+naf(is)*(zaf(is)*tef+tif)
         tro = tro+am(is)*mp*naf(is)
       enddo
       do is = 0, ns_ext-1
         nxf = (na_ext(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &          na_ext(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &         (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         zxf = (za_ext(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &          za_ext(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &         (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         taf = (ta_ext(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &          ta_ext(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &         (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
         tpr = tpr+nxf*(zxf*tef+taf)
         tro = tro+am_ext(is)*mp*nxf
       enddo
       cst = sqrt(tpr/tro)
       vte = sqrt(tef/me)
*    ..compute Bohm fchi, fche
       t0 = 0.0_R8
       do is = 0, ns-1
        t0 = t0+zaf(is)*fna(nx,iy,0,0,is)
       enddo
       do is = 0, ns_ext-1
        zxf = (za_ext(nx,iy,is)*hx(leftix(nx,iy),leftiy(nx,iy))+
     &         za_ext(leftix(nx,iy),leftiy(nx,iy),is)*hx(nx,iy))/
     &        (hx(leftix(nx,iy),leftiy(nx,iy))+hx(nx,iy))
        t0 = t0-zxf*fa_ext(nx,iy,0,0,is)
       enddo
       fchi = qe*t0
       if(t0.lt.0.0_R8 .and. region(nx,iy,0).ne.0 .and.
     &  bottomiy(nx,iy).ne.-2 .and. topiy(nx,iy).ne.ny+1) then
         wrong_flow = .true.                                             !xpb
         if(diagno.ge.0) then
           write(*,*) 'b2stbc_east:iy,fchi',iy,fchi
           write(*,*) zaf(ismain),fna(nx,iy,0,0,ismain)
           write(*,'(a,4g14.7)') 'Bx,Bz,Ud,Ua =',
     &      wbbl(nx,iy,0),wbbl(nx,iy,2),-udf(ismain),ua(nx,iy,ismain)
         endif
c         fchi=0.0_R8
       endif
!       if (po(nx,iy).ge.0.0_R8) then                                    !sv
         fche = 
     .    sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1*expu2(-qe*pof/tef)
!       else                                                             !sv
!         fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*nef*vte*s1                  !sv
!       endif                                                            !sv
*    ..loop over species
       do is = 0, ns-1
*     ..compute csb
        csb = sqrt((zaf(is)*tef+tif)/(am(is)*mp))
*     .. account for diamagnetic drifts
        csbm = csb-udf(is)*abs(wbbl(nx,iy,2)/wbbl(nx,iy,0))              !xpb
        cstm = cst-udf(is)*abs(wbbl(nx,iy,2)/wbbl(nx,iy,0))              !xpb
        if(.not.is_neutral(is)) then                                     !xpb
          csbm = max(0.0_R8,csbm)                                        !xpb
          cstm = max(0.0_R8,cstm)                                        !xpb
        endif                                                            !xpb
        if(csbm.gt.csb .and. diagno.ge.2)
     >   write(*,*) 'east:iy,is,csb,csbm',iy,is,csb,csbm
        if(cstm.gt.cst .and. diagno.ge.2)
     >   write(*,*) 'east:iy,is,cst,cstm',iy,is,cst,cstm
*     ..compute particle source
        sna0(nx,iy,0,is) = sna0(nx,iy,0,is)+
     &    cbsna(0,is,ireg)*s0+
     &    cbsna(6,is,ireg)*max(0.0_R8,-ua(nx,iy,is)*pbs(nx,iy,0))
        sna0(nx,iy,1,is) = sna0(nx,iy,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csbm*s0+
     &    cbsna(3,is,ireg)*csbm*s1+
     &    cbsna(4,is,ireg)*cstm*s0+
     &    cbsna(5,is,ireg)*cstm*s1+
     &    cbsna(7,is,ireg)*max(0.0_R8,ua(nx,iy,is)*pbs(nx,iy,0))
*     ..compute parallel momentum source
        smo0(nx,iy,0,is) = smo0(nx,iy,0,is)+cbsmo(0,is,ireg)*s0
        smo0(nx,iy,1,is) = smo0(nx,iy,1,is)+cbsmo(7,is,ireg)*s0          !rfs
        smo0(nx,iy,2,is) = smo0(nx,iy,2,is)+
     &    cbsmo(1,is,ireg)*csbm*pbs(nx,iy,0)+                            !xpb
     &    cbsmo(2,is,ireg)*cstm*pbs(nx,iy,0)+                            !xpb
     &    cbsmo(5,is,ireg)*max(0.0_R8,+fna(nx,iy,0,0,is))/na(nx,iy,is)     !sv
        smo0(nx,iy,3,is) = smo0(nx,iy,3,is)+
     &    cbsmo(3,is,ireg)*s0+
     &    cbsmo(4,is,ireg)*s1+
     &    cbsmo(6,is,ireg)*max(0.0_R8,fna(nx,iy,0,0,is))/na(nx,iy,is)
*     ..compute atom heat source
        t0 = na(nx,iy,is)/ni(nx,iy,0)
        shi0(nx,iy,0) = shi0(nx,iy,0)+
     &    cbshi(0,is,ireg)*t0*s0+
     &    cbshi(6,is,ireg)*max(0.0_R8,-fna(nx,iy,0,0,is))
        shi0(nx,iy,1) = shi0(nx,iy,1)+
     &    cbshi(1,is,ireg)*t0*s0
        shi0(nx,iy,3) = shi0(nx,iy,3)+
     &    cbshi(2,is,ireg)*t0*csbm*s0+
     &    cbshi(3,is,ireg)*t0*csbm*s1+
     &    cbshi(4,is,ireg)*t0*cstm*s0+
     &    cbshi(5,is,ireg)*t0*cstm*s1+
     &    cbshi(7,is,ireg)*max(0.0_R8,fna(nx,iy,0,0,is))/ni(nx,iy,0)
       enddo
       shi_test=
     1  shi0(nx,iy,0)+
     2  shi0(nx,iy,1)*ti(nx,iy)+
     3  shi0(nx,iy,2)*ni(nx,iy,0)+
     4  shi0(nx,iy,3)*ni(nx,iy,0)*ti(nx,iy)
       if(shi_test.gt.0.0_R8.and.diagno.ge.2) then
         write(*,*) 'shi_test>0: E ',
     1    iy,shi0(nx,iy,:),ti(nx,iy)/qe,ni(nx,iy,0),s0,s1,csbm,cstm,
     2    ireg,cbshi(:,:,ireg)
       endif
*    ..compute electron heat source
       she0(nx,iy,0) = she0(nx,iy,0)+
     &   cbshe(0,ireg)*s0+
     &   cbshe(6,ireg)*max(0.0_R8,-fne(nx,iy,0,0))
       she0(nx,iy,1) = she0(nx,iy,1)+
     &   cbshe(1,ireg)*s0
       she0(nx,iy,3) = she0(nx,iy,3)+
     &   cbshe(2,ireg)*vte*s0+
     &   cbshe(3,ireg)*vte*s1+
     &   cbshe(4,ireg)*cstm*s0+
     &   cbshe(5,ireg)*cstm*s1+
     &   cbshe(7,ireg)*max(0.0_R8,fne(nx,iy,0,0))/ne(nx,iy)
       she_test=she0(nx,iy,0)+
     2  she0(nx,iy,1)*te(nx,iy)+
     3  she0(nx,iy,2)*ne(nx,iy)+
     4  she0(nx,iy,3)*ne(nx,iy)*te(nx,iy)
       if(she_test.gt.0.0_R8.and.diagno.ge.2) then
         write(*,*) 'she_test>0: E ',
     1    iy,she0(nx,iy,:),te(nx,iy)/qe,ne(nx,iy),s0,s1,vte,cstm,
     2    ireg,cbshe(:,ireg)
       endif
*    ..compute charge source
       sch0(nx,iy,2) = sch0(nx,iy,2)+
     &   cbsch(0,ireg)*qe*vte*s0+
     &   cbsch(2,ireg)*qe*vte*s1+
     &   cbsch(4,ireg)*qe*cstm*s1
       sch0(nx,iy,3) = sch0(nx,iy,3)+
     &   cbsch(1,ireg)*qe**2*vte*s0/te(nx,iy)+
     &   cbsch(3,ireg)*qe**2*vte*s1/te(nx,iy)+
     &   cbsch(5,ireg)*qe**2*cstm*s1/te(nx,iy)
       if (secmodel.eq.1) then
         if (ne(nx,iy)/(rza(nx,iy,ismain)*na(nx,iy,ismain)).gt.1.1_R8) 
     >    then
           write(*,*) 'Large ne contribution from impurities!'
           write(*,*) 'Kinetic s.e.c. model may not be valid!'
         endif
         call eseec0 (tef,tif,nef,pof,wbbl(nx,iy,0),wbbl(nx,iy,3),
     &    0,nx,iy,nx,ny,ns,ismain,seec,diagno)                           !cfs
       else
         seec = 1.0_R8 - cbsch(7,ireg)
       endif
!       if (po(nx,iy).ge.0.0_R8) then                                    !sv
         t0 = max(cbsch(6,ireg)*fchi,(1.0_R8-seec)*fche)*qe/tef
!       else                                                             !sv
!         t0 = 0.0_R8                                                    !sv
!       endif                                                            !sv
       sch0(nx,iy,0) = sch0(nx,iy,0)+
     &  (1.0_R8-seec)*fche-cbsch(6,ireg)*fchi+t0*po(nx,iy)
       sch0(nx,iy,1) = sch0(nx,iy,1)-t0
 404   continue
      enddo

* ..return
      if(ncall.eq.0 .or. main_call) ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_bas

      end subroutine b2stbc_bas

