*-----------------------------------------------------------------------
*.specification

      subroutine b2tfrn (nx, ny, ns, guard_flows,                        !srv 14.07.10 {
     &  facdrift, fac_ExB, rza, rz2,                                     !srv 13.01.17
     &  vol, hx, hy, hz, qz, gs, qc, bb, pbs, ua, na, ne, te, ti, po,   
     &  fch_p, fchvispar, fchvisper, fchinert, fchanml, fchvisq,         !srv 20.09.06
     &  fne, cvla, cdna, cdpa, chce, chci, vaecrb, wadia, ve, wedia,
     &  uadia, fna_eir, fne_eir, fhe_eir, fhi_eir)
      use b2mod_types
      use b2mod_boundary_namelist
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmfs
      use b2mod_transport_fun                                            !srv 13.01.17
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, guard_flows
      real (kind=R8) ::
     &  rza(-1:nx,-1:ny,0:ns-1), rz2(-1:nx,-1:ny,0:ns-1),                !srv 13.01.17
     &  pbs(-1:nx,-1:ny,0:1), po(-1:nx,-1:ny), fch_p(-1:nx,-1:ny,0:1),
     &  fchvispar(-1:nx,-1:ny,0:1), fchinert(-1:nx,-1:ny,0:1),           !srv 20.09.06
     &  fchvisper(-1:nx,-1:ny,0:1), fchvisq(-1:nx,-1:ny,0:1),
     &  fchanml(-1:nx,-1:ny,0:1),
     &  ua(-1:nx,-1:ny,0:ns-1), na(-1:nx,-1:ny,0:ns-1), 
     &  fne(-1:nx,-1:ny,0:1), 
     &  cvla(-1:nx,-1:ny,0:1,0:ns-1), cdna(-1:nx,-1:ny,0:1,0:ns-1),
     &  cdpa(-1:nx,-1:ny,0:1,0:ns-1),
     &  chce(-1:nx,-1:ny,0:1), chci(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny), vaecrb(-1:nx,-1:ny,0:1,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1), uadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  veecrb(-1:nx,-1:ny,0:1),
     &  vol(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), 
     2  hz(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1), gs(-1:nx,-1:ny,0:2), 
     &  qc(-1:nx,-1:ny), bb(-1:nx,-1:ny,0:3), ve(-1:nx,-1:ny), 
     2  facdrift(-1:nx,-1:ny), fac_ExB(-1:nx,-1:ny), 
     3  te(-1:nx,-1:ny), ti(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  wedia(-1:nx,-1:ny,0:1), fna_eir(-1:nx,-1:ny,0:1, 0:ns-1),
     &  fne_eir(-1:nx,-1:ny,0:1),
     &  fhe_eir(-1:nx,-1:ny,0:1), fhi_eir(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TFRN computes the particle and energy fluxes between
*     neighbouring cells which are to be provided to Eirene
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the particle fluxes
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, isb, ib
      integer, save :: no_Ptncr_x_he = 0, no_Ptncr_y_he = 0              !srv 06.07.06
      integer, save :: no_Ptncr_x_hi = 0, no_Ptncr_y_hi = 0              !srv 06.07.06
      integer, save :: no_Ptncr_x_co = 0, no_Ptncr_y_co = 0              !srv 06.07.06
      integer, save :: no_hybr = 0
      integer, save :: drift_style = 1
      integer, save :: iout = 0,
     &                 iout_b2npc9 = 0, iout_b2npco = 0, iout_b2npmo =0
      integer, save :: style = 2                                         !srv 13.01.17
      real (kind=R8), save :: xvecrb = 1.0_R8, xwdia = 1.0_R8            !srv 27.01.12
      character chns*3                                                   !srv 16.09.02
      real (kind=R8), save :: fch_pTe = 1.0_R8                           !srv 30.12.09
      real (kind=R8), save :: xfrhie = 1.0_R8                            !srv 05.11.09 30.01.10
      real (kind=R8), save :: xcur = 0.0_R8                              !xpb 27.09.17
      real (kind=R8), save :: ycur = 1.0_R8                              !srv 07.07.08
      real (kind=R8), save :: vis_per = 1.0_R8, vis_q = 1.0_R8           !srv 02.10.17
      real (kind=R8) ::
     *   dif, flo, cutlo, fle_eir, cone_eir, floe_eir, chce_eir, scur,   !srv 14.07.10
     &   floi_eir, chci_eir, flo_eir, con_eir, term(-1:nx,-1:ny),
     &   vsprcur                                                         !srv 02.04.18
      real (kind=R8) :: 
     &   wedia_(-1:nx,-1:ny,0:1), nel(-1:nx,-1:ny), pb(-1:nx,-1:ny)      !srv 14.07.10 }
      parameter (cutlo=1.0e-60_R8)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny), wrk1(-1:nx,-1:ny), wrk2(-1:nx,-1:ny),
     *   zeff(-1:nx,-1:ny), c071(-1:nx,-1:ny)                            !srv 13.01.17
      real (kind=R8) :: t0,t1,t2,cvmgp
      save ncall
*   ..procedures
      real (kind=R8) ::
     *   hybr, hy1
      external subini, subend, xertst, ipgeti, ipgetr
      external b2xvsg, b2xvff, intfaceh, intfacev, cvmgp
      intrinsic abs, max, min, mod, sqrt
      hy1(ix,iy) = hy(ix,iy)*qz(ix,iy,1)
      hybr(flo,dif) = max(dif,abs(flo/2.),sqrt(dif**2+flo**2/4.)) !!!
*   ..initialisation
      data ncall/0/ 

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tfrn')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tfrn_iout', iout)
        call ipgeti ('b2npc9_iout', iout_b2npc9)
        call ipgeti ('b2npco_iout', iout_b2npco)
        call ipgeti ('b2npmo_iout', iout_b2npmo)
        call ipgetr ('b2tfhe_fch_pTe', fch_pTe)                          !srv 30.12.09
        call ipgetr ('b2tfnb_xfrhie', xfrhie)
        call ipgetr ('b2tfrn_xvecrb', xvecrb)                            !srv 27.01.12
        call ipgetr ('b2tfrn_xwdia', xwdia)                              !srv 27.01.12
        call ipgeti ('b2tfnb_no_Ptncr_x', no_Ptncr_x_co)                 !srv 06.07.06
        call ipgeti ('b2tfnb_no_Ptncr_y', no_Ptncr_y_co)                 !srv 06.07.06
        call ipgeti ('b2tfhi_no_Ptncr_x', no_Ptncr_x_hi)                 !srv 06.07.06
        call ipgeti ('b2tfhi_no_Ptncr_y', no_Ptncr_y_hi)                 !srv 06.07.06
        call ipgeti ('b2tfhe_no_Ptncr_x', no_Ptncr_x_he)                 !srv 06.07.06
        call ipgeti ('b2tfhe_no_Ptncr_y', no_Ptncr_y_he)                 !srv 06.07.06
        call ipgetr ('b2tfnb_xcur', xcur)                                !xpb 27.09.17
        call ipgetr ('b2tfnb_ycur', ycur)                                !srv 17.07.08
        call ipgetr ('b2tfnb_vis_per', vis_per)                          !srv 02.10.17
        call ipgetr ('b2tfnb_vis_q', vis_q)                              !srv 02.10.17
        call ipgeti ('b2tfnb_no_hybr', no_hybr)
        call ipgeti ('b2tfnb_drift_style', drift_style)
        call ipgeti ('b2sigp_style', style)                              !srv 13.01.17
      endif
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of ne, cdnb
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

      do is = 0, ns-1
*    ..compute partial pressure
       call b2xppb (nx, ny, rza(-1,-1,is), na(-1,-1,is), te, ti, pb)     !srv 01.10.99
*   ..compute wrk0
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
           wrk0(ix,iy) = 0.0_R8
         else
           wrk0(ix,iy) = !pbs(ix,iy,0)*                                  !srv 20.06.08
     &      (vol(ix,iy)*ua(leftix(ix,iy),leftiy(ix,iy),is)+
     &       vol(leftix(ix,iy),leftiy(ix,iy))*ua(ix,iy,is))/
     &      (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
         endif
        enddo
       enddo
       call b2xxgs (nx, ny, 0.0_R8, wrk0, 1)
*   ..compute wrk1
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk1(ix,iy) = 0.0_R8
         else
          wrk1(ix,iy) = (pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))/
     &     (2.0_R8*(pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy)))
         endif
        enddo
       enddo
       call b2xxgs (nx, ny, 0.0_R8, wrk1, 1)
*   ..compute wrk2
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          wrk2(ix,iy) = 0.0_R8
         elseif(leftix(leftix(ix,iy),leftiy(ix,iy)).eq.-2) then
          wrk2(ix,iy) = 0.0_R8
         elseif(rightix(ix,iy).eq.nx+1) then
          wrk2(ix,iy) = 0.0_R8
         else
          wrk2(ix,iy) = wrk1(leftix(ix,iy),leftiy(ix,iy))-
     &     2.0_R8*wrk1(ix,iy)+wrk1(rightix(ix,iy),rightiy(ix,iy))
         endif
        enddo
       enddo
       call b2xxgs (nx, ny, 0.0_R8, wrk2, 1)
*   ..compute fna_eirx
       do iy = -1, ny
        do ix = -1, nx
         if(leftix(ix,iy).eq.-2) then
          fna_eir(ix,iy,0,is)     = 0.0_R8                               !srv 14.07.10
         else

csw 30sep2013, added check for no_hybr
          if (nint(rza(ix,iy,is)).eq.1.and.nint(zamax(is)).eq.1) then
           scur =
     &      (fchinert(ix,iy,0)+fchvispar(ix,iy,0)+
     &       vis_per*fchvisper(ix,iy,0)                                  !srv 02.10.17
     &      +fchanml(ix,iy,0)+vis_q*fchvisq(ix,iy,0))*xcur/qe/           !srv 02.10.17
     &     ((na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))/2.0_R8)
          else
            scur = 0.0_R8
          end if
          if(no_hybr .ne. 0) then
           term(ix,iy) = xfrhie*wrk2(ix,iy)*0.5_R8*abs(pbs(ix,iy,0))*
     &     sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/
     &     (am(is)*mp*
     &     (na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))))*
     &     min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &     (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
           con_eir = cdna(ix,iy,0,is)
           flo_eir = cvla(ix,iy,0,is)+pbs(ix,iy,0)*wrk0(ix,iy)+          !srv 01.07.05
     &      term(ix,iy)                                                  !iyv 19.04.07
           if (drift_style.eq.0) then
            flo_eir = flo_eir + 0.5_R8*
     1       ((vaecrb(ix,iy,0,is)*xvecrb+wadia(ix,iy,0,is)*xwdia)*
     2         vol(ix,iy)/hx(ix,iy)+
     3        (vaecrb(leftix(ix,iy),leftiy(ix,iy),0,is)*xvecrb+
     4         wadia(leftix(ix,iy),leftiy(ix,iy),0,is)*xwdia)*
     5         vol(leftix(ix,iy),leftiy(ix,iy))/
     6         hx(leftix(ix,iy),leftiy(ix,iy)))
           else
            flo_eir = flo_eir +
     4      (vaecrb(ix,iy,0,is)*xvecrb+wadia(ix,iy,0,is)*xwdia)*
     5      (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/               !srv 01.03.00
     6      (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))
           end if
           flo_eir = flo_eir + scur
           t0 = (na(ix,iy,is) 
     .        + na(leftix(ix,iy),leftiy(ix,iy),is)) / 2.0_R8
           t1 = (2.0_R8 * na(ix,iy,is) 
     .                  * na(leftix(ix,iy),leftiy(ix,iy),is)) /
     .      (cutlo + na(ix,iy,is) + na(leftix(ix,iy),leftiy(ix,iy),is))
           t2 = cvmgp(t0, t1, flo_eir*
     1      (na(ix,iy,is)-na(leftix(ix,iy),leftiy(ix,iy),is)))
           fna_eir(ix,iy,0,is) = flo_eir * t2-
     &      con_eir*(na(ix,iy,is)-na(leftix(ix,iy),leftiy(ix,iy),is))-
     &      cdpa(ix,iy,0,is)*(pb(ix,iy)-pb(leftix(ix,iy),leftiy(ix,iy)))

          else
           term(ix,iy) = xfrhie*wrk2(ix,iy)*0.5_R8*abs(pbs(ix,iy,0))*    !iyv 19.04.07 11.01.13
     &      sqrt((pb(leftix(ix,iy),leftiy(ix,iy))+pb(ix,iy))/(am(is)*mp*
     &       (na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))))*
     &      min(vol(leftix(ix,iy),leftiy(ix,iy)),vol(ix,iy))/
     &      (vol(leftix(ix,iy),leftiy(ix,iy))+vol(ix,iy))
           flo_eir = cvla(ix,iy,0,is)+pbs(ix,iy,0)*wrk0(ix,iy)+          !srv 14.07.10 {
     &      term(ix,iy)
           if (drift_style.eq.0) then
            flo_eir = flo_eir + 0.5_R8*
     1       ((vaecrb(ix,iy,0,is)*xvecrb+wadia(ix,iy,0,is)*xwdia)*
     2         vol(ix,iy)/hx(ix,iy)+
     3        (vaecrb(leftix(ix,iy),leftiy(ix,iy),0,is)*xvecrb+
     4         wadia(leftix(ix,iy),leftiy(ix,iy),0,is)*xwdia)*
     5         vol(leftix(ix,iy),leftiy(ix,iy))/
     6         hx(leftix(ix,iy),leftiy(ix,iy)))
           else
            flo_eir = flo_eir +
     4       (vaecrb(ix,iy,0,is)*xvecrb+wadia(ix,iy,0,is)*xwdia)*
     5       (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     6       (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))                 !srv 14.07.10 }
           endif
           flo_eir = flo_eir + scur
           if(no_Ptncr_x_co.eq.0) then                                   !srv 13.10.06 {
            con_eir = hybr(flo_eir,cdna(ix,iy,0,is))
           else
            con_eir = cdna(ix,iy,0,is)
           endif
           fna_eir(ix,iy,0,is) = flo_eir*                                !srv 14.07.10
     1      (na(leftix(ix,iy),leftiy(ix,iy),is)+na(ix,iy,is))/2.0_R8-
     &      con_eir*(na(ix,iy,is)-na(leftix(ix,iy),leftiy(ix,iy),is))
          endif ! no_hybr = 0
         endif
        enddo
       enddo
*   ..compute fna_eiry
       do iy = -1, ny
        do ix = -1, nx
         if(bottomiy(ix,iy).eq.-2) then
           fna_eir(ix,iy,1,is)     = 0.0_R8                              !srv 14.07.10
         else

csw 30sep2013, added check for no_hybr
          if(no_hybr .ne. 0) then
           con_eir = cdna(ix,iy,1,is)
           flo_eir = cvla(ix,iy,1,is)
           if(drift_style.eq.0) then
            flo_eir = flo_eir + 0.5_R8*
     1       ((vaecrb(ix,iy,1,is)*xvecrb+wadia(ix,iy,1,is)*xwdia)*
     2         vol(ix,iy)/hy1(ix,iy)+
     3        (vaecrb(bottomix(ix,iy),bottomiy(ix,iy),1,is)*xvecrb+
     4         wadia(bottomix(ix,iy),bottomiy(ix,iy),1,is)*xwdia)*
     5         vol(bottomix(ix,iy),bottomiy(ix,iy))/
     6         hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           else
            flo_eir = flo_eir +
     1       (vaecrb(ix,iy,1,is)*xvecrb+wadia(ix,iy,1,is)*xwdia)*
     2       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     3       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           end if
           if(nint(rza(ix,iy,is)).eq.1.and.nint(zamax(is)).eq.1) then    !srv 20.09.06
            scur =
     &        (fchvispar(ix,iy,1)+fchinert(ix,iy,1)+
     &         vis_q*fchvisq(ix,iy,1)+fchanml(ix,iy,1))/                 !srv 02.10.17
     &        (qe*(na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &        na(ix,iy,is))/2.0_R8)
            vsprcur = vis_per*fchvisper(ix,iy,1)/qe                      !srv 02.04.18
           else
            scur = 0.0_R8
            vsprcur = 0.0_R8                                             !srv 02.04.18
          endif
           flo_eir = flo_eir + scur
           t0 = (na(ix,iy,is) +na(bottomix(ix,iy),bottomiy(ix,iy),is))/
     1      2.0_R8
           t1 = 
     .      (2.0_R8 * na(ix,iy,is) 
     .              * na(bottomix(ix,iy),bottomiy(ix,iy),is)) /
     .      (cutlo + na(ix,iy,is) 
     .             + na(bottomix(ix,iy),bottomiy(ix,iy),is))
           t2 = cvmgp(t0, t1, 
     .      flo_eir*(na(ix,iy,is)
     .              -na(bottomix(ix,iy),bottomiy(ix,iy),is)))
           fna_eir(ix,iy,1,is) = flo_eir * t2-
     &      con_eir*(na(ix,iy,is)
     &              -na(bottomix(ix,iy),bottomiy(ix,iy),is))-
     &      cdpa(ix,iy,1,is)*(pb(ix,iy)
     &              -pb(bottomix(ix,iy),bottomiy(ix,iy))) + vsprcur       !srv 02.04.18
          else
           flo_eir = cvla(ix,iy,1,is)
           if(drift_style.eq.0) then
            flo_eir = flo_eir + 0.5_R8*
     1       ((vaecrb(ix,iy,1,is)*xvecrb+wadia(ix,iy,1,is)*xwdia)*
     2         vol(ix,iy)/hy1(ix,iy)+
     3        (vaecrb(bottomix(ix,iy),bottomiy(ix,iy),1,is)*xvecrb+
     4         wadia(bottomix(ix,iy),bottomiy(ix,iy),1,is)*xwdia)*
     5         vol(bottomix(ix,iy),bottomiy(ix,iy))/
     6         hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           else
            flo_eir = flo_eir +
     1       (vaecrb(ix,iy,1,is)*xvecrb+wadia(ix,iy,1,is)*xwdia)*
     2       (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     3       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
           end if
           if(nint(rza(ix,iy,is)).eq.1.and.nint(zamax(is)).eq.1) then    !srv 20.09.06
             scur =
     &        (fchinert(ix,iy,1)+fchvispar(ix,iy,1)+
     &          fchanml(ix,iy,1)+vis_q*fchvisq(ix,iy,1))*ycur/          !srv 02.10.17
     &        (qe*(na(bottomix(ix,iy),bottomiy(ix,iy),is)+
     &             na(ix,iy,is))/2.0_R8)
             vsprcur = vis_per*fchvisper(ix,iy,1)/qe                    !srv 02.04.18
           else
             scur = 0.0_R8
             vsprcur = 0.0_R8                                            !srv 02.04.18
          endif
           flo_eir = flo_eir + scur
           if(no_Ptncr_y_co.eq.0) then                                   !srv 06.07.06 {
             con_eir = hybr(flo_eir,cdna(ix,iy,1,is))
           else
             con_eir = cdna(ix,iy,1,is)
           endif
           fna_eir(ix,iy,1,is) = flo_eir*                                !srv 14.07.10
     &      (na(bottomix(ix,iy),bottomiy(ix,iy),is)+na(ix,iy,is))
     &      /2.0_R8-con_eir*
     &      (na(ix,iy,is)-na(bottomix(ix,iy),bottomiy(ix,iy),is))
     &      + vsprcur                                                    !srv 02.04.18
         endif !no_hybr =0
         endif
        enddo
       enddo
      enddo !is
c Drift terms
      wedia =0.0_R8
      wedia_=0.0_R8

      if((maxval(facdrift).ne.0.0_R8 .or. 
     &    maxval(fac_ExB).ne.0.0_R8)) then
        wrk0=0.0_R8  ! compute .... on the cell faces                    !srv 14.07.10 {
        do iy=0,ny
          do ix=-1,nx
            wrk0(ix,iy)=(ne(ix,iy)*te(ix,iy)-
     &       ne(bottomix(ix,iy),bottomiy(ix,iy))*
     &       te(bottomix(ix,iy),bottomiy(ix,iy)))/
     &       (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))*2.0_R8
          enddo
        enddo
        do iy=-1,ny  ! compute ... on the cell centers
          do ix=-1,nx
            if(bottomiy(ix,iy).eq.-2) then
              wedia_(ix,iy,0) = wrk0(ix,topiy(ix,iy))
            elseif(topiy(ix,iy).eq.ny+1) then
              wedia_(ix,iy,0) = wrk0(ix,iy)
            else
              wedia_(ix,iy,0) = 0.5_R8*
     &         (wrk0(topix(ix,iy),topiy(ix,iy))+wrk0(ix,iy))
            endif
            wedia_(ix,iy,0) = wedia_(ix,iy,0)/
     &       (qe*ne(ix,iy))*(bb(ix,iy,2)/bb(ix,iy,3)**2)                 !srv 15.02.12
          enddo
        enddo
        call intfaceh(nx,ny,leftix,leftiy,vol,wedia_(-1,-1,0),
     &                                         wedia(-1,-1,0))

        wrk0=0.0_R8  ! compute ... on the cell faces
        do iy=-1,ny
          do ix=-1,nx
            if(leftix(ix,iy).ne.-2) then
              wrk0(ix,iy)=(ne(ix,iy)*te(ix,iy)-
     &         ne(leftix(ix,iy),leftiy(ix,iy))*
     &         te(leftix(ix,iy),leftiy(ix,iy)))/
     &         (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))*2.0_R8
            endif
          enddo
        enddo
        do iy=-1,ny  ! compute .... on the cell centers
          do ix=-1,nx
            if(leftix(ix,iy).eq.-2) then
              wedia_(ix,iy,1) = wrk0(rightix(ix,iy),iy)
            elseif(rightix(ix,iy).eq.nx+1) then
              wedia_(ix,iy,1) = wrk0(ix,iy)
            else
              wedia_(ix,iy,1) = 0.5_R8*
     &         (wrk0(rightix(ix,iy),rightiy(ix,iy))+wrk0(ix,iy))
            endif
              wedia_(ix,iy,1) = - wedia_(ix,iy,1)/
     &         (qe*ne(ix,iy))*(bb(ix,iy,2)/bb(ix,iy,3)**2)               !srv 15.02.12
          enddo
        enddo
        call intfacev(nx,ny,bottomix,bottomiy,vol,wedia_(-1,-1,1),
     &                                            wedia(-1,-1,1))        !srv 14.07.10 }
      endif

* ..compute electron velocity
      call intfaceh(nx,ny,leftix,leftiy,vol,ne,nel)                      !srv 14.07.10 {
      do iy = -1, ny
        do ix = -1, nx
          if(leftix(ix,iy).eq.-2) then
            ve(ix,iy) = 0.0e0_R8
          else     ! ve has been myltiplied by vol*bx/hx
            ve(ix,iy)=fne(ix,iy,0)/nel(ix,iy)
          endif
        enddo
      enddo

      if(style.eq.2) then
        call b2xzef (nx, ny, ns, rz2, na, ne, zeff)
        do iy = -1, ny
          do ix = -1, nx
            c071(ix,iy) = fce2n(zeff(ix,iy))
          enddo
        enddo
      else
        c071 = 0.71e0_R8                                                 !srv 13.01.17
      endif

* ..compute fluxes
      do is = 0,ns-1
        if(nint(zamax(is)).eq.1) then
*   ..compute fne_eirx, fhe_eirx, fhi_eirx
          do iy = -1, ny
            do ix = -1, nx
              if(leftix(ix,iy).eq.-2) then
                fne_eir(ix,iy,0)     = 0.0_R8                            !srv 14.07.10
                fhe_eir(ix,iy,0)     = 0.0_R8                            !srv 14.07.10
                fhi_eir(ix,iy,0)     = 0.0_R8                            !srv 14.07.10
              else     ! ve has been myltiplied by vol*bx/hx
                fle_eir = ve(ix,iy)                                      !srv 14.07.10 {
                if(drift_style.eq.0) then
                 fle_eir = fle_eir + 0.5_R8*
     1            ((veecrb(ix,iy,0)*xvecrb+wedia(ix,iy,0)*xwdia)*
     2              vol(ix,iy)/hx(ix,iy)+
     3             (veecrb(leftix(ix,iy),leftiy(ix,iy),0)*xvecrb+
     4              wedia(leftix(ix,iy),leftiy(ix,iy),0)*xwdia)*
     5              vol(leftix(ix,iy),leftiy(ix,iy))/
     6              hx(leftix(ix,iy),leftiy(ix,iy)))
                else
                 fle_eir = fle_eir +
     &            (veecrb(ix,iy,0)*xvecrb+wedia(ix,iy,0)*xwdia)*
     &            (vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))/
     &            (hx(ix,iy)+hx(leftix(ix,iy),leftiy(ix,iy)))            !srv 14.07.10 }
                end if
                if(no_Ptncr_x_co.eq.0) then                              !srv 13.10.06 {
                  cone_eir = hybr(fle_eir,cdna(ix,iy,0,is))
                else
                  cone_eir = cdna(ix,iy,0,is)
                endif
                fne_eir(ix,iy,0) =                                       !srv 14.07.10
     &           fle_eir*(ne(leftix(ix,iy),leftiy(ix,iy))+ne(ix,iy))
     &            /2.0_R8-
     &            cone_eir*(ne(ix,iy)-ne(leftix(ix,iy),leftiy(ix,iy))) 
                floe_eir = 2.5_R8*fne_eir(ix,iy,0)
     &           -c071(ix,iy)*fch_p(ix,iy,0)/qe*fch_pTe                  !srv 30.12.09 13.01.17
                if(no_Ptncr_x_he.eq.0) then                              !srv 06.07.06 {
                  chce_eir = hybr(floe_eir,chce(ix,iy,0))
                else
                  chce_eir = chce(ix,iy,0)
                endif
                fhe_eir(ix,iy,0) = floe_eir*                             !srv 14.07.10
     &           (te(leftix(ix,iy),leftiy(ix,iy))+te(ix,iy))
     &           /2.0_R8-
     &           chce_eir*(te(ix,iy)-te(leftix(ix,iy),leftiy(ix,iy)))
                floi_eir = 0.0_R8
                do isb = 0, ns-1
                  if(.not.is_neutral(isb)) then
                    floi_eir = floi_eir + fna_eir(ix,iy,0,isb)
                  endif
                enddo
                floi_eir = 2.5_R8*floi_eir
                if(no_Ptncr_x_hi.eq.0) then                              !srv 06.07.06 {
                  chci_eir = hybr(floi_eir,chci(ix,iy,0))
                else
                  chci_eir = chci(ix,iy,0)
                endif
                fhi_eir(ix,iy,0) = floi_eir*                             !srv 14.07.10
     &           (ti(leftix(ix,iy),leftiy(ix,iy))+ti(ix,iy))
     &           /2.0_R8-
     &           chci_eir*(ti(ix,iy)-ti(leftix(ix,iy),leftiy(ix,iy)))
              endif
            enddo
          enddo

*   ..  compute fne_eiry, fhe_eiry
          do iy = -1, ny
            do ix = -1, nx
              if(bottomiy(ix,iy).eq.-2) then
                fne_eir(ix,iy,1)     = 0.0_R8
                fhe_eir(ix,iy,1)     = 0.0_R8
                fhi_eir(ix,iy,1)     = 0.0_R8
              else
!             flo_eir = cvlb(ix,iy,1)+(wedia(ix,iy,1)+vbecrb(ix,iy,1))*  !srv 14.07.10
!     &     (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
!     &      (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
!     &        (fchinert(ix,iy,1)+fchvispar(ix,iy,1)+fchvisper(ix,iy,1)
!     &         +fchanml(ix,iy,1)+fchvisq(ix,iy,1))*ycur/
!     &        (qe*(nb(bottomix(ix,iy),bottomiy(ix,iy))+nb(ix,iy))/2)
                if(drift_style.eq.0) then
                  fle_eir = 0.5_R8*
     1             ((veecrb(ix,iy,1)*xvecrb+wedia(ix,iy,1)*xwdia)*
     2               vol(ix,iy)/hy1(ix,iy)+
     3              (veecrb(bottomix(ix,iy),bottomiy(ix,iy),1)*xvecrb+
     4               wedia(bottomix(ix,iy),bottomiy(ix,iy),1)*xwdia)*
     5               vol(bottomix(ix,iy),bottomiy(ix,iy))/
     6               hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                else
                  fle_eir =
     1             (veecrb(ix,iy,1)*xvecrb+wedia(ix,iy,1)*xwdia)*
     2             (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/
     3             (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))
                end if
                if(no_Ptncr_y_co.eq.0) then                              !srv 13.10.06 {
                  cone_eir = hybr(fle_eir,cdna(ix,iy,1,is))
                else
                  cone_eir = cdna(ix,iy,1,is)
                endif
                fne_eir(ix,iy,1) = fle_eir*                              !srv 14.07.10
     &           (ne(bottomix(ix,iy),bottomiy(ix,iy))+ne(ix,iy))/
     &           2.0_R8-cone_eir*
     &           (ne(ix,iy)-ne(bottomix(ix,iy),bottomiy(ix,iy)))
                floe_eir = 2.5_R8*fne_eir(ix,iy,1)
                if(no_Ptncr_y_he.eq.0) then                              !srv 13.10.06 {
                  chce_eir = hybr(floe_eir,chce(ix,iy,1))
                else
                  chce_eir = chce(ix,iy,1)
                endif
                fhe_eir(ix,iy,1) = floe_eir*                             !srv 14.07.10
     &           (te(bottomix(ix,iy),bottomiy(ix,iy))+te(ix,iy))/
     &           2.0_R8-chce_eir*
     &           (te(ix,iy)-te(bottomix(ix,iy),bottomiy(ix,iy)))
                floi_eir = 0.0_R8
                do isb = 0, ns-1
                  if(.not.is_neutral(isb)) then
                    floi_eir = floi_eir + fna_eir(ix,iy,1,isb)
                  endif
                enddo
                floi_eir = 2.5_R8*floi_eir
                if(no_Ptncr_y_hi.eq.0) then                              !srv 13.10.06 {
                  chci_eir = hybr(floi_eir,chci(ix,iy,1))
                else
                  chci_eir = chci(ix,iy,1)
                endif
                fhi_eir(ix,iy,1) = floi_eir*                             !srv 14.07.10
     &           (ti(bottomix(ix,iy),bottomiy(ix,iy))+ti(ix,iy))/
     &           2.0_R8-chci_eir*
     &           (ti(ix,iy)-ti(bottomix(ix,iy),bottomiy(ix,iy)))
              endif
            enddo
          enddo
        else

        endif
      enddo

      if(guard_flows.eq.0) then                                          !xpb
cxpb turn off flows between boundary cells
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fna_eir(ix,iy,1,:)=0.0_R8
              fne_eir(ix,iy,1)=0.0_R8
              fhe_eir(ix,iy,1)=0.0_R8
              fhi_eir(ix,iy,1)=0.0_R8
              wedia(ix,iy,1)=0.0_R8
            endif
          enddo
        enddo

cxpb turn off poloidal flows in core boundary
        do ix = -1, nx
          if(on_closed_surface(ix,-1)) then
            fna_eir(ix,-1,0,:) = 0.0_R8
            fne_eir(ix,-1,0) = 0.0_R8
            fhe_eir(ix,-1,0) = 0.0_R8
            fhi_eir(ix,-1,0) = 0.0_R8
            wedia(ix,-1,0)   = 0.0_R8
          endif
        enddo
      elseif(guard_flows.eq.1) then                                      !srv 20.06.08 {
*   .. do nothing
      elseif(guard_flows.eq.2) then
*     turn off flows between boundary cells in west and east
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).eq.-2 .or. rightix(ix,iy).eq.nx+1) then
              fna_eir(ix,iy,1,:) = 0.0_R8
!             fne_eir(ix,iy,1)   = 0.0_R8
!             fhe_eir(ix,iy,1)   = 0.0_R8
!             fhi_eir(ix,iy,1)   = 0.0_R8
!             wedia(ix,iy,1)=0.0_R8
!             ubdia(ix,iy,1)=0.0_R8
            endif
          enddo
        enddo
*     turn off flows between boundary cells in south and in north
        do iy = -1, ny, ny+1
          fna_eir(:,iy,0,:)       = 0.0_R8
!         fne_eir(:,iy,0)         = 0.0_R8
!         fhe_eir(:,iy,0)         = 0.0_R8
!         fhi_eir(:,iy,0)         = 0.0_R8
!         wedia(:,iy,0) = 0.0_R8
!         ubdia(:,iy,0) = 0.0_R8
        enddo                                                            !srv 20.06.08 }
      endif                                                              !xpb

*  INTERIOR BOUNDARY AT DEAD REGIONS
      DO IB=1,NNISO
        DO IX=NXISO1(IB), NXISO2(IB)
          DO IY=NYISO1(IB), NYISO2(IB)
            wedia(ix,iy,:)  =0.0_R8
            fna_eir(ix,iy,:,:)=0.0_R8
            fne_eir(ix,iy,:)=0.0_R8
            fhe_eir(ix,iy,:)=0.0_R8
            fhi_eir(ix,iy,:)=0.0_R8
          enddo
        enddo
      enddo

csw 18feb2014 remove fluxes between guard cells 
      do ix=-1,nx
        do iy=-1,ny
          if(iffboundx(ix,iy).ne.0) then
            fna_eir(ix,iy,0,:)=0.0_R8
            fne_eir(ix,iy,0)=0.0_R8
            fhe_eir(ix,iy,0)=0.0_R8
            fhi_eir(ix,iy,0)=0.0_R8
          endif
          if(iffboundy(ix,iy).ne.0) then
            fna_eir(ix,iy,1,:)=0.0_R8
            fne_eir(ix,iy,1)=0.0_R8
            fhe_eir(ix,iy,1)=0.0_R8
            fhi_eir(ix,iy,1)=0.0_R8
          endif
        enddo
      enddo

      if(iout.ne.0) then                                                 !srv 18.06.08 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out(70,nx,ny,fna_eir(-1,-1,0,is),
     .     'b2tfrn_fna_eirx'//chns)
          call my_out(70,nx,ny,fna_eir(-1,-1,1,is),
     .     'b2tfrn_fna_eiry'//chns)
          call my_out(70,nx,ny,vaecrb(-1,-1,0,is),'b2tfrn_vecrbx'//chns)
          call my_out(70,nx,ny,vaecrb(-1,-1,1,is),'b2tfrn_vecrby'//chns)
          call my_out(70,nx,ny,wadia(-1,-1,0,is),'b2tfrn_wadiax'//chns)
          call my_out(70,nx,ny,wadia(-1,-1,1,is),'b2tfrn_wadiay'//chns)
          call my_out(70,nx,ny,cdna(-1,-1,0,is),'b2tfrn_cdnax'//chns)
          call my_out(70,nx,ny,cdna(-1,-1,1,is),'b2tfrn_cdnay'//chns)
          call my_out(70,nx,ny,ua(-1,-1,is),'b2tfrn_ua'//chns)
          call my_out(70,nx,ny,na(-1,-1,is),'b2tfrn_na'//chns)
        enddo
        call my_out(70,nx,ny,fne_eir(-1,-1,0),'b2tfrn_fne_eirx')
        call my_out(70,nx,ny,fne_eir(-1,-1,1),'b2tfrn_fne_eiry')
        call my_out(70,nx,ny,fhe_eir(-1,-1,0),'b2tfrn_fhe_eirx')
        call my_out(70,nx,ny,fhe_eir(-1,-1,1),'b2tfrn_fhe_eiry')
        call my_out(70,nx,ny,fhi_eir(-1,-1,0),'b2tfrn_fhi_eirx')
        call my_out(70,nx,ny,fhi_eir(-1,-1,1),'b2tfrn_fhi_eiry')
        call my_out(70,nx,ny,fch_p(-1,-1,0),'b2tfrn_fch_px')
        call my_out(70,nx,ny,fch_p(-1,-1,1),'b2tfrn_fch_py')
        call my_out(70,nx,ny,fchvispar(-1,-1,0),'b2tfrn_fchvisparx')
        call my_out(70,nx,ny,fchvispar(-1,-1,1),'b2tfrn_fchvispary')
        call my_out(70,nx,ny,fchinert(-1,-1,0),'b2tfrn_fchinertx')
        call my_out(70,nx,ny,fchinert(-1,-1,1),'b2tfrn_fchinerty')
        call my_out(70,nx,ny,fchanml(-1,-1,0),'b2tfrn_fchanmlx')
        call my_out(70,nx,ny,fchanml(-1,-1,1),'b2tfrn_fchanmly')
        call my_out(70,nx,ny,wedia(-1,-1,0),'b2tfrn_wediax')
        call my_out(70,nx,ny,wedia(-1,-1,1),'b2tfrn_wediay')
        call my_out(70,nx,ny,ve,'b2tfrn_ve')
        call my_out(70,nx,ny,ne,'b2tfrn_ne')
        call my_out(70,nx,ny,te,'b2tfrn_te')
        call my_out(70,nx,ny,ti,'b2tfrn_ti')
        call my_out(70,nx,ny,po,'b2tfrn_po')
      endif                                                              !srv 18.06.08 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tfrn

      end subroutine b2tfrn
