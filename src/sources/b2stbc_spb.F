*-----------------------------------------------------------------------
*.specification
      subroutine b2stbc_spb (nx, ny, ns, facdrift,                      !sv 09.01.01
     1  na, ua, uadia, vaecrb, vadia, te, ti, po,                       !sv 14.10.99 12.05.00 09.01.01
     &  ne, ni, floe_noc, floi_noc, fna, fna_fcor, fnaPSch, fhe, fhi,   !sv 13.10.06
     &  fhePSch, fhiPSch, fne, fni, fch,                                !sv 25.02.00 21.10.02 28.10.02 13.10.06
     &  fchdia, fchin, fch_p, fchvispar, fchvisper, fchvisq, fchinert,  !sv 21.10.02
     &  sna0, smo0, she0, shi0, sch0, wrong_flow)                       !sv 09.07.99 01.07.05
      use b2mod_types
      use b2mod_geo
      use b2mod_math
      use b2mod_indirect
      use b2mod_anomalous_transport
      use b2mod_boundary_namelist
      use b2mod_constants
      use b2mod_geo_corner
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmfs
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  facdrift(-1:nx,-1:ny)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     1  uadia(-1:nx,-1:ny,0:1,0:ns-1),
     1  vaecrb(-1:nx,-1:ny,0:1,0:ns-1),                                 !sv 14.10.99
     1  vadia(-1:nx,-1:ny,0:1,0:ns-1),                                  !sv 12.05.00
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *   fna(-1:nx,-1:ny,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1),
     &   fhi(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &   fhePSch(-1:nx,-1:ny,0:1), fhiPSch(-1:nx,-1:ny,0:1),            !sv 13.10.06
     &   fch_p(-1:nx,-1:ny,0:1), fchvispar(-1:nx,-1:ny,0:1),            !sv 11.10.02 21.10.02 {
     &   fchvisper(-1:nx,-1:ny,0:1), fchin(-1:nx,-1:ny,0:1),
     &   fchvisq(-1:nx,-1:ny,0:1), fchdia(-1:nx,-1:ny,0:1),
     &   fchinert(-1:nx,-1:ny,0:1),                                     !sv 21.10.02 }
     &   fne(-1:nx,-1:ny,0:1), fni(-1:nx,-1:ny,0:1),
     &   fna_fcor(-1:nx,-1:ny,0:1,0:ns-1),                              !sv 25.02.00
     &   fnaPSch(-1:nx,-1:ny,0:1,0:ns-1),                               !sv 13.10.06
     &   floe_noc(-1:nx,-1:ny,0:1), floi_noc(-1:nx,-1:ny,0:1)           !sv 28.10.02
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBC_SPB provides linearized expressions for the sources of
*     density, momentum, electron heat, atom heat, and charge, due to
*     the boundary conditions.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ireg
      logical wrong_flow
      real (kind=R8) ::
     *  sna0ep, she0ep, shi0ep,
     1  tpz, trz, csb, cst, vte, s0, s1, fche, fchi, t0,
     4  csbm, cstm,   ! cstm1                                         !sv 14.10.99
     6  s0hz, s1hz, vwest, veast, vminus, vplus  !, v0                !sv 09.01.01 18.07.03 28.01.05
      real (kind=R8) ::
     *  pz(-1:nx,-1:ny)                                               !sv 18.01.02
      integer, save :: istyle_average_pr = 0, istyle_fchi = 0         !sv 08.08.02 14.02.05
      integer, save :: idrifts_to_the_wall = 0                        !sv 02.03.05
      integer, save :: ivminus = 1, ivplus = 1                        !sv 24.03.05
      integer, save :: iform_corr_dia = 0, istyle_cs = 0              !sv 25.07.05 25.04.08
      integer, save :: mdf_fnb = 0, mdf_fhe = 0, mdf_fhi = 0          !sv 03.10.06
      integer       :: num_inner                                      !sv 26.11.02
      real (kind=R8) ::
     1   facdriftm
     &  ,pzaverage, pzsum, volsum                                     !sv 24.03.99 18.01.02
     &  ,ua2sum, ua2average(0:ns-1)                                   !sv 07.04.00 15.02.02
     &  ,uasum, uaaverage(0:ns-1)                                     !sv 15.02.02
     &  ,nesum, tisum, tiaverage, tesum                               !sv 15.12.00 25.01.02
     &  ,nasum, naaverage(0:ns-1), bbsum, bbaverage                   !sv 22.01.02 18.10.04
      real (kind=R8) :: bzhzsum, bzhzaverage                          !lk 20.11.07
      real (kind=R8) :: fchy_inertia                                  !sv 07.04.00
      real (kind=R8) :: fchy_biasing_cur                              !sv 31.08.00
      real (kind=R8) :: fchy_dia                                      !sv
      real (kind=R8) , save :: y_dia = 1.0e0_R8, y_inertia = 1.0e0_R8 !sv 25.04.00
      real (kind=R8) , save :: biasing_cur = 0.0e0_R8                 !sv 31.08.00
      real (kind=R8) , save :: fchy_dia_at_PR = 0.0e0_R8,             !sv 09.07.99
     &    fchy_in_at_PR = 0.0e0_R8, fchy_vispar_at_PR = 0.0e0_R8,     !sv 21.10.02
     &    fchy_visq_at_PR = 0.0e0_R8, fchy_visper_at_PR = 0.0e0_R8,   !sv 21.10.02
     &    fchy_inert_at_PR = 0.0e0_R8                                 !sv 21.10.02
      real (kind=R8) , save :: fchy_dia_at_north = 0.0e0_R8,          !sv 09.07.99
     &                         fchy_in_at_north = 0.0e0_R8,           !sv 21.10.02
     &                         fchy_vispar_at_north = 0.0e0_R8,       !sv 21.10.02
     &                         fchy_visq_at_north = 0.0e0_R8,         !sv 21.10.02
     &                         fchy_visper_at_north = 0.0e0_R8,       !sv 21.10.02
     &                         fchy_inert_at_north = 0.0e0_R8         !sv 21.10.02
      real (kind=R8) , save :: cbc = 1.0e0_R8                         !sv 14.10.99
      real (kind=R8) , save :: po_in_lower = 0.0e0_R8,
     &                         po_in_upper = 0.0e0_R8,                !sv 16.11.04
     &                         corr_flux   = 0.0e0_R8                 !sv 25.07.05
      real (kind=R8) , save :: phm0 = 1.0e0_R8, phm1 = 1.0e0_R8,
     &                         phm2 = 1.0e0_R8, ubig0 = 0.0e0_R8      !iyv 22.04.08
      character chns*3, chk*1                                         !sv 11.07.05
      integer, save :: iout = 0, k                                    !sv 11.07.05
      save ncall, sna0ep, she0ep, shi0ep
*   ..procedures
      intrinsic abs, sqrt
      real (kind=R8) ::
     *   b2ssum
      external xertst, sfill, b2ssum, b2saxpy
      external b2xvsg, b2xvff, b2xvps, b2xxgs, b2xppz                  !sv 18.01.02
*   ..initialization
      data ncall/0/, sna0ep/1.0e-36_R8/, she0ep/1.0e-36_R8/,
     1 shi0ep/1.0e-36_R8/

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbc_')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stbc_phm0', phm0)                                  !sv 17.01.07
       call ipgetr ('b2stbc_phm1', phm1)                                  !sv 17.01.07
       call ipgetr ('b2stbc_phm2', phm2)                                  !sv 17.01.07
       call ipgetr ('b2stbc_cbc', cbc)                                    !sv 14.10.99
       call ipgetr ('b2stbc_y_dia', y_dia)                                !sv 25.04.00
       call ipgetr ('b2stbc_y_inertia', y_inertia)                        !sv 25.04.00
       call ipgetr ('b2stbc_biasing_cur', biasing_cur)                    !sv 31.08.00
       call ipgeti ('b2stbc_average_pr', istyle_average_pr)               !sv 08.08.02
       call ipgeti ('b2stbc_istyle_fchi', istyle_fchi)                    !sv 14.02.05
       call ipgeti ('b2stbc_idrifts_to_the_wall', idrifts_to_the_wall)    !sv 02.03.05
       call ipgeti ('b2stbc_ivplus', ivplus)                              !sv 24.03.05
       call ipgeti ('b2stbc_ivminus', ivminus)                            !sv 24.03.05
       call ipgeti ('b2stbc_iform_corr_dia', iform_corr_dia)              !sv 25.07.05
       call ipgetr ('b2stbc_corr_flux_in_core', corr_flux)                !sv 25.07.05
       call ipgeti ('b2tfnb_mdf', mdf_fnb)                                !sv 13.10.06 {
       call ipgeti ('b2tfhe_mdf', mdf_fhe)
       call ipgeti ('b2tfhi_mdf', mdf_fhi)                                !sv 13.10.06 }
       call get_boundary_namelist
       call ipgetr ('b2stbc_sna0ep', sna0ep)
       call xertst (0.le.sna0ep, 'faulty internal parameter sna0ep')
       call ipgetr ('b2stbc_she0ep', she0ep)
       call xertst (0.le.she0ep, 'faulty internal parameter she0ep')
       call ipgetr ('b2stbc_shi0ep', shi0ep)
       call xertst (0.le.shi0ep, 'faulty internal parameter shi0ep')
       call ipgetr ('b2stbc_fchy_dia_at_PR', fchy_dia_at_PR)             !sv 09.07.99
       call ipgetr ('b2stbc_fchy_in_at_PR', fchy_in_at_PR)               !sv 21.10.02
       call ipgetr ('b2stbc_fchy_vispar_at_PR', fchy_vispar_at_PR)       !sv 21.10.02
       call ipgetr ('b2stbc_fchy_visq_at_PR', fchy_visq_at_PR)           !sv 21.10.02
       call ipgetr ('b2stbc_fchy_visper_at_PR', fchy_visper_at_PR)       !sv 21.10.02
       call ipgetr ('b2stbc_fchy_inert_at_PR', fchy_inert_at_PR)         !sv 21.10.02
       call ipgetr ('b2stbc_fchy_dia_at_north', fchy_dia_at_north)       !sv 09.07.99
       call ipgetr ('b2stbc_fchy_in_at_north', fchy_in_at_north)         !sv 21.10.02
       call ipgetr ('b2stbc_fchy_vispar_at_north', fchy_vispar_at_north) !sv 21.10.02
       call ipgetr ('b2stbc_fchy_visq_at_north', fchy_visq_at_north)     !sv 21.10.02
       call ipgetr ('b2stbc_fchy_visper_at_north', fchy_visper_at_north) !sv 21.10.02
       call ipgetr ('b2stbc_fchy_inert_at_north', fchy_inert_at_north)   !sv 21.10.02
       call ipgetr ('b2stbc_po_in_lower', po_in_lower)                   !sv 16.11.04
       call ipgetr ('b2stbc_po_in_upper', po_in_upper)                   !sv 16.11.04
       call ipgeti ('b2stbc_iout', iout)                                 !sv 11.07.05
       call ipgetr ('b2stbc_ubig0', ubig0)                               !iyv 22.04.08
       call ipgeti ('b2stbc_istyle_cs', istyle_cs)                       !sv 25.04.08
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test state
       call b2xvps (nx, ny, ns, ne, na, ua, te, ti, po,
     &              fna, fhe, fhi, fch)
*    ..test edge values of fne, fni
       call b2xvff (nx, ny, fne, 'fne')
       call b2xvff (nx, ny, fni, 'fni')
      endif
*   ..initialize sources to 0
      call sfill (n2*2*ns, 0.0e0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0e0_R8, smo0, 1)
      call sfill (n2*4, 0.0e0_R8, she0, 1)
      call sfill (n2*4, 0.0e0_R8, shi0, 1)
      call sfill (n2*4, 0.0e0_R8, sch0, 1)
      facdriftm=maxval(facdrift)                                      !srv 29.07.08

* ..compute standard form volume sources
*   (This code placed here for want of a better location)
*   ..low-level particle source
      do is = 0, ns-1
       call b2saxpy (n2, sna0ep, vol, 1, sna0(-1,-1,0,is), 1)
      enddo
      call b2saxpy (n2, she0ep, vol, 1, she0(-1,-1,0), 1)
      call b2saxpy (n2, shi0ep, vol, 1, shi0(-1,-1,0), 1)

      if(boundary_namelist.eq.0) then
       write(*,*) 'Doing spb B2.5-style bc'

* ..compute conditions at south boundary
*   ..initialize ireg
       ireg = cbirso
       call b2xppz (nx, ny, ns, ne, na, te, ti, pz)               !sv 18.01.02
*   ..find the average pressure for inner flux surface            !sv 24.03.99
       if (istyle_average_pr.eq.0) then                           !sv 08.08.02
        iy = -1                                                   !sv 24.03.99
        pzsum =0.0e0_R8                                           !sv 24.03.99 18.01.02
        volsum=0.0e0_R8                                           !sv 24.03.99
        do ix = -1, nx                                            !sv 26.11.02
         if (on_closed_surface(ix,iy)) then                       !sv 26.11.02
          pzsum=pzsum + pz(ix,iy)*vol(ix,iy)                      !sv 24.03.99 18.01.02
          volsum=volsum + vol(ix,iy)                              !sv 24.03.99
         endif                                                    !sv 26.11.02
        enddo                                                     !sv 24.03.99
        pzaverage=pzsum/volsum                                    !sv 24.03.99 18.01.02
       else                                                       !sv 08.08.02
        iy = 0                                                    !sv 24.03.99  07.06.02 {
        pzsum =0.0e0_R8                                           !sv 24.03.99 18.01.02
        volsum=0.0e0_R8                                           !sv 24.03.99
        num_inner=0                                               !sv 26.11.02
        do ix = -1, nx                                            !sv 26.11.02
         if (on_closed_surface(ix,iy)) then                       !sv 26.11.02
          pzsum=pzsum
     *    + (pz(ix,iy)+pz(bottomix(ix,iy),bottomiy(ix,iy)))       !sv 24.03.99 18.01.02
          volsum=volsum + vol(ix,iy)                              !sv 24.03.99
          num_inner=num_inner+1                                   !sv 26.11.02
         endif                                                    !sv 26.11.02
        enddo                                                     !sv 24.03.99
        pzaverage=pzsum*0.5e0_R8/num_inner                        !sv 24.03.99 18.01.02 07.06.02 } 26.11.02
       endif                                                      !sv 08.08.02
*   ..find the average ne, te, ti, bb for inner flux surface      !sv 15.12.00 18.10.04
       bbsum = 0.0e0_R8                                           !sv 18.10.04
       nesum = 0.0e0_R8                                           !sv 15.12.00 25.01.02
       tesum = 0.0e0_R8                                           !sv 15.12.00
       tisum = 0.0e0_R8                                           !sv 15.12.00
       bzhzsum  = 0.0e0_R8                                        !lk 20.11.07
       bbaverage  = 0.0e0_R8                                      !lk 20.11.07
       do ix = -1, nx                                             !sv 26.11.02
        if (on_closed_surface(ix,iy)) then                        !sv 26.11.02
         bbsum = bbsum + bb(ix,iy,3)*vol(ix,iy)                   !sv 18.10.04
         nesum = nesum + ne(ix,iy)*vol(ix,iy)                     !sv 15.12.00 25.01.02
         tisum = tisum + ti(ix,iy)*vol(ix,iy)                     !sv 15.12.00
         tesum = tesum + te(ix,iy)*vol(ix,iy)                     !sv 15.12.00
         bzhzsum = bzhzsum+bb(ix,iy,2)*vol(ix,iy)**2/gs(ix,iy,2)  !lk 20.11.07
        endif                                                     !sv 26.11.02
       enddo                                                      !sv 15.12.00
       bbaverage=bbsum/volsum                                     !sv 18.10.04
       tiaverage=tisum/volsum                                     !sv 15.12.00
       bzhzaverage=bzhzsum/volsum                                 !lk 20.11.07
*   ..find the average ua**2, na for inner flux surface           !sv 22.01.02 {
       do is = 0, ns-1
        uasum  = 0.0e0_R8                                         !sv 15.02.02
        ua2sum = 0.0e0_R8                                         !sv 07.04.00 15.02.02
        nasum  = 0.0e0_R8
        do ix = -1, nx                                            !sv 26.11.02
         if (on_closed_surface(ix,iy)) then                       !sv 26.11.02
          ua2sum = ua2sum + ua(ix,iy,is)**2*vol(ix,iy)            !sv 07.04.00 15.02.02
          uasum  = uasum  + ua(ix,iy,is)*vol(ix,iy)               !sv 15.02.02
          nasum  = nasum  + na(ix,iy,is)*vol(ix,iy)
         endif                                                    !sv 26.11.02
        enddo
        ua2average(is) = ua2sum/volsum                            !sv 07.04.00 15.02.02
        uaaverage(is) = uasum/volsum                              !sv 07.04.00 15.02.02 16.04.02
        naaverage(is) = nasum/volsum
       enddo                                                      !sv 22.01.02 }
*   ..loop over boundary cells
       do ix = -1, nx
*    ..update ireg
    1   if (cbrbrk(ireg)*nx.lt.ix+0.5e0_R8) then
         call xertst (ireg.lt.cbirso+cbnrso-1,
     &    'error in computing region on south boundary')
         ireg = ireg+1
         goto 1
        endif
*    ..compute s0, s1
        s0 = vol(ix,-1)/hy1(ix,-1)
        s1 = abs(pbs(ix,0,1))
*    ..compute s0hz, s1hz                                         !sv 09.01.01
        s0hz = hz(ix,-1)*vol(ix,-1)/hy1(ix,-1)                    !sv 09.01.01
        s1hz = abs(pbshz(ix,0,1))                                 !sv 09.01.01
*    ..compute cst, vte
        tpz = 0
        trz = 0
        do is = 0, ns-1
         if (is_neutral(is)) cycle
         tpz = tpz+na(ix,-1,is)*(rza(ix,-1,is)*te(ix,-1)+ti(ix,-1))
         trz = trz+am(is)*mp*na(ix,-1,is)
        enddo
        cst = sqrt(tpz/trz)
        vte = sqrt(te(ix,-1)/me)
*    ..compute Bohm fchi, fche
        t0 = 0
        do is = 0, ns-1
         t0 = t0-rza(ix,0,is)*fna(ix,0,1,is)
        enddo
        fchi = qe*t0
        fche = sqrt(1/(2*pi))*qe*ne(ix,-1)*vte*s1*
     &   expu2(-qe*po(ix,-1)/te(ix,-1))
        fchy_dia = 0.0e0_R8                                       !sv 01.07.99
        fchy_inertia = 0.0e0_R8                                   !sv 07.04.00
        fchy_biasing_cur = 0.0e0_R8                               !sv 31.08.00
        if(facdriftm.ne.0.0e0_R8) then                            !sv 09.07.99
         if (.not.on_closed_surface(ix,-1)) then                  !sv 26.11.02
*           for PR
           if(fchy_dia_at_PR .ne. 0.0e0_R8) then                  !sv 09.07.99
             iy = 0                                               !sv 21.10.02
             fchy_dia = fchdia(ix,iy,1)                           !sv 21.10.02
           else
             fchy_dia = 0.0e0_R8                                  !sv 09.07.99
           endif
         else                                                     !sv
*        for inner surfaces
           iy=0                                                   !sv 01.07.99
           fchy_dia =                                             !sv 01.07.99
     &        -pzaverage*facdrift(ix,iy)*phm0      ! *signmf    !sv 24.03.99 01.07.99 21.04.00 18.01.02 2 11.01.13
!     &      *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2)) !sv 01.07.99 20.11.07 {
!     ^       /2.0e0_R8                                           !sv 01.07.99
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))   !sv 01.07.99
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2)) !sv 01.07.99
     &      *bzhzaverage                                          !lk 20.11.07 }
     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !sv 01.07.99
           fchy_inertia = 0.0e0_R8                                !sv 25.06.03
           do is = 0, ns-1                                        !sv 25.06.03
             if (.not.is_neutral(is)) then                        !sv 25.06.03
               fchy_inertia = fchy_inertia                        !sv 07.04.00 25.06.03
     &        -0.5e0_R8*am(is)*mp*facdrift(ix,iy)*phm1  ! *signmf !sv 07.04.00 21.04.00 21.10.02 11.01.13
!     &        *(                                                 !sv 07.04.00
!     &          na(ix,iy,is)+na(bottomix(ix,iy),bottomiy(ix,iy),is)   !sv 07.04.00
!     &         )/2.0e0_R8                                        !sv 07.04.00
     &        *naaverage(is)                                      !lk 20.11.07                       !sv 07.04.00
     &        *ua2average(is)                                     !sv 07.04.00 15.02.02
!     &        *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2)) !sv 07.04.00
!     &         /2.0e0_R8                                         !sv 07.04.00
!     &        *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy))) !sv 07.04.00
!     &        /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2)) !sv 07.04.00
     &        *bzhzaverage                                        !lk 20.11.07
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !sv 07.04.00
             endif
           enddo
           fchy_biasing_cur = biasing_cur                         !sv 31.08.00
     &        *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))  !sv 31.08.00
     &        /(hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))  !sv 31.08.00 07.05.08
         endif                                                    !sv
        endif                                                     !sv 09.07.99
*    ..loop over species
        do is = 0, ns-1
*     ..compute csb
         csb = sqrt((rza(ix,-1,is)*te(ix,-1)+ti(ix,-1))/(am(is)*mp))
*     ..compute particle source
*!!     cbsna(6,,) and cbsna(7,,) to be introduced with ua*pbs(,,1)
         iy = 0                                                   !sv 23.04.01
         sna0(ix,-1,0,is) = sna0(ix,-1,0,is)+
     &    cbsna(0,is,ireg)*s0
         if(mdf_fnb.ne.0 .and. .not.is_neutral(is)) then          !sv 09.11.06 {
           if(.not.on_closed_surface(ix,-1))      ! for PR region only
     &      sna0(ix,-1,0,is) = sna0(ix,-1,0,is) - fnaPSch(ix,0,1,is)
           endif                                                  !sv 09.11.06 }
           if(is.eq.1) then                                       !sv 08.02.02
          sna0(ix,-1,0,is) = sna0(ix,-1,0,is)                     !sv 08.02.02
     &    +cbsna(6,is,ireg)*(                                     !sv 08.02.02
!     &     cvla(ix,iy,1)+   ! remember if you need
     &     uadia(ix,iy,1,is)*                                     !sv 08.02.02
     &    (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))/      !sv 01.03.00
     &     (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))      !sv 01.03.00 07.05.08
     &                      )                                     !sv 08.02.02
     &    *(na(bottomix(ix,iy),bottomiy(ix,iy),is) + na(ix,iy,is))/2  !sv 23.04.01
         endif                                                    !sv 08.02.02
         sna0(ix,-1,1,is) = sna0(ix,-1,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csb*s0+
     &    cbsna(3,is,ireg)*csb*s1+
     &    cbsna(4,is,ireg)*cst*s0+
     &    cbsna(5,is,ireg)*cst*s1
*    ..drifts to the wall                                         !sv 02.03.05 {
         if(idrifts_to_the_wall .ne. 0) then
           if(.not.is_neutral(is) .and.
     &        cbsna(2,is,ireg).ne.0.0e0_R8) then
             if(mdf_fnb.eq.0) then
               sna0(ix,-1,1,is) = sna0(ix,-1,1,is)+
     *         min(0.0e0_R8,s0*(vadia(ix,0,1,is)+vaecrb(ix,0,1,is)))
             else
               sna0(ix,-1,1,is) = sna0(ix,-1,1,is)+
     *         min(0.0e0_R8,s0*vaecrb(ix,0,1,is))
             endif
           endif
         endif                                                    !sv 02.03.05 }
*     ..compute momentum source
         if (on_closed_surface(ix,iy)) then ! only core region and drifts !sv 27.07.05
           if (.not.is_neutral(is)) then                          !sv 15.02.02 {
            select case (iform_corr_dia)                          !sv 25.07.05
            case (0)
             iy = -1                                              !sv 15.02.02
             smo0(ix,-1,0,is) = smo0(ix,-1,0,is)+                 !sv 09.01.01 01.07.05
     &       s0hz*(am(is)*mp)*naaverage(is)                       !sv 15.02.02
!     &     *(2.0e0_R8*vadia(ix,iy,1,is))
     &        *(-2.0e0_R8*facdrift(ix,iy)*phm2*       ! *signmf   !sv 07.06.00 17.05.99 25.04.00 15.02.02 11.01.13
     &           tiaverage*bb(ix,iy,2)/(qe*rza(ix,iy,is)*hx(ix,iy)) !sv 15.02.02
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !sv 05.03.99 15.02.02
     &        )*uaaverage(is)                                     !sv 15.02.02
     &        +corr_flux*s0hz                                     !sv 08.09.03 25.07.05
            case (1)
             iy = -1                                              !sv 15.02.02
             smo0(ix,-1,0,is) = smo0(ix,-1,0,is)+                 !sv 09.01.01 01.07.05
     &       s0hz*(am(is)*mp)*naaverage(is)                       !sv 15.02.02
!     &     *(2.0e0_R8*vadia(ix,iy,1,is))
     &       *(-2.0e0_R8*facdrift(ix,iy)*phm2*        ! *signmf   !sv 07.06.00 17.05.99 25.04.00 15.02.02 11.01.13
     &          tiaverage*bb(ix,iy,2)/(qe*rza(ix,iy,is)*hx(ix,iy))  !sv 15.02.02
     &        *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy)) !sv 05.03.99 15.02.02
     &       )*ua(ix,iy,is)                                       !sv 25.07.05
     &       +corr_flux*s0hz                                      !sv 08.09.03 25.07.05
            end select                                            !sv 25.07.05
           endif                                                  !sv 15.02.02
         endif                                                    !sv 27.07.05
         smo0(ix,-1,0,is) = smo0(ix,-1,0,is)+cbsmo(0,is,ireg)*s0hz  !sv 09.01.01 01.07.05
         smo0(ix,-1,1,is) = smo0(ix,-1,1,is)                      !sv 01.07.05
         smo0(ix,-1,2,is) = smo0(ix,-1,2,is)+                     !sv 01.07.05
     &    cbsmo(1,is,ireg)*(-csb)*pbshz(ix,0,1)+                  !sv 09.01.01
     &    cbsmo(2,is,ireg)*(-cst)*pbshz(ix,0,1)+                  !sv 09.01.01
     &    cbsmo(5,is,ireg)*max(0.0e0_R8,fna_fcor(ix,0,1,is))
     &    /na(ix,-1,is)                                           !sv 25.02.00
         smo0(ix,-1,2,is) = smo0(ix,-1,2,is)                      !sv 18.10.04 { 01.07.05
     &    +cbsmo(7,is,ireg)*s0hz
     &    *(bbaverage/bb(ix,iy,3))                                !sv 18.10.04 }
         smo0(ix,-1,3,is) = smo0(ix,-1,3,is) +                    !sv 01.07.05
     &    cbsmo(3,is,ireg)*s0hz+                                  !sv 09.01.01
     &    cbsmo(4,is,ireg)*s1hz+                                  !sv 09.01.01
     &    cbsmo(6,is,ireg)                                        !sv 14.10.02
     &             *max(0.0e0_R8,-fna_fcor(ix,0,1,is))/na(ix,-1,is) !sv 25.02.00
*     ..compute atom heat source
         t0 = na(ix,-1,is)/ni(ix,-1,0)
         shi0(ix,-1,0) = shi0(ix,-1,0)+
     &    cbshi(0,is,ireg)*t0*s0+
     &    cbshi(6,is,ireg)*max(0.0e0_R8,fna(ix,0,1,is))
         shi0(ix,-1,1) = shi0(ix,-1,1)+
     &    cbshi(1,is,ireg)*t0*s0
         shi0(ix,-1,3) = shi0(ix,-1,3)+
     &    cbshi(2,is,ireg)*t0*csb*s0+
     &    cbshi(3,is,ireg)*t0*csb*s1+
     &    cbshi(4,is,ireg)*t0*cst*s0+
     &    cbshi(5,is,ireg)*t0*cst*s1
!     &   +cbshi(7,is,ireg)*max(0.0e0_R8,-fna(ix,0,1,is))/ni(ix,-1,0) !sv 28.10.02
        enddo
        shi0(ix,-1,1) = shi0(ix,-1,1)+                            !sv 28.10.02
     &    cbshi(7,0,ireg)*max(0.0e0_R8,-floi_noc(ix,0,1))         !sv 28.10.02
*    ..compute electron heat source
        she0(ix,-1,0) = she0(ix,-1,0)+
     &    cbshe(0,ireg)*s0+
     &    cbshe(6,ireg)*max(0.0e0_R8,fne(ix,0,1))
        she0(ix,-1,1) = she0(ix,-1,1)+
     &    cbshe(1,ireg)*s0
     &   +cbshe(7,ireg)*max(0.0e0_R8,-floe_noc(ix,0,1))           !sv 28.10.02
        she0(ix,-1,3) = she0(ix,-1,3)+
     &    cbshe(2,ireg)*vte*s0+
     &    cbshe(3,ireg)*vte*s1+
     &    cbshe(4,ireg)*cst*s0+
     &    cbshe(5,ireg)*cst*s1
!     &  +cbshe(7,ireg)*max(0.0e0_R8,-fne(ix,0,1))/ne(ix,-1)      !sv 28.10.02
*    ..compute charge source
        sch0(ix,-1,2) = sch0(ix,-1,2)+
!     &   cbsch(0,ireg)*qe*vte*s0+                                !sv 11.03.01
     &   cbsch(2,ireg)*qe*vte*s1+
     &   cbsch(4,ireg)*qe*cst*s1
        sch0(ix,-1,3) = sch0(ix,-1,3)+
!     &   cbsch(1,ireg)*qe**2*vte*s0/te(ix,-1)+                   !sv 11.03.01
     &   cbsch(3,ireg)*qe**2*vte*s1/te(ix,-1)+
     &   cbsch(5,ireg)*qe**2*cst*s1/te(ix,-1)
        t0 = max(cbsch(6,ireg)*fchi,cbsch(7,ireg)*fche)*qe/te(ix,-1)
        sch0(ix,-1,0) = sch0(ix,-1,0) + fchy_dia*y_dia            !sv 01.07.99 07.04.00 25.04.00
     &                  + fchy_inertia*y_inertia                  !sv 25.04.00
     &                  + fchy_biasing_cur                        !sv 31.08.00
     &                  + cbsch(0,ireg)*s0                        !sv 11.03.01
!       if (ix.lt.leftcut(1) .or. ix.ge.rightcut(1)) then         !sv 21.10.02
        if (.not.on_closed_surface(ix,-1)) then                   !sv 26.11.02
          ! for PR
          if(fchy_in_at_PR .ne. 0.0e0_R8) then                    !sv 21.10.02
           sch0(ix,-1,0) = sch0(ix,-1,0) +fchin(ix,0,1)
          endif                                                   !sv 21.10.02
          if(fchy_vispar_at_PR .ne. 0.0e0_R8) then                !sv 21.10.02
           sch0(ix,-1,0) = sch0(ix,-1,0) +fchvispar(ix,0,1)
          endif                                                   !sv 21.10.02
          if(fchy_visq_at_PR .ne. 0.0e0_R8) then                  !sv 21.10.02
           sch0(ix,-1,0) = sch0(ix,-1,0) +fchvisq(ix,0,1)
          endif                                                   !sv 21.10.02
          if(fchy_visper_at_PR .ne. 0.0e0_R8) then                !sv 21.10.02
           sch0(ix,-1,0) = sch0(ix,-1,0) +fchvisper(ix,0,1)
          endif                                                   !sv 21.10.02
          if(fchy_inert_at_PR .ne. 0.0e0_R8) then                 !sv 21.10.02
           sch0(ix,-1,0) = sch0(ix,-1,0) +fchinert(ix,0,1)
          endif                                                   !sv 21.10.02
        endif
        sch0(ix,-1,0) = sch0(ix,-1,0)+
     &   cbsch(7,ireg)*fche-cbsch(6,ireg)*fchi+t0*po(ix,-1)
        sch0(ix,-1,1) = sch0(ix,-1,1)-t0
        sch0(ix,-1,1) = sch0(ix,-1,1)+cbsch(1,ireg)*s0            !sv 11.03.01
       enddo

* ..compute conditions at north boundary
*   ..initialize ireg
       ireg = cbirno
*   ..loop over boundary cells
       do ix = -1, nx
*    ..update ireg
    2   if (cbrbrk(ireg)*nx.lt.ix+0.5e0_R8) then
         call xertst (ireg.lt.cbirno+cbnrno-1,
     &    'error in computing region on north boundary')
         ireg = ireg+1
         goto 2
        endif
*    ..compute s0, s1
        s0 = vol(ix,ny)/hy1(ix,ny)
        s1 = abs(pbs(ix,ny,1))
*    ..compute s0hz, s1hz                                         !sv 09.01.01
        s0hz = hz(ix,ny)*vol(ix,ny)/hy1(ix,ny)                    !sv 09.01.01
        s1hz = abs(pbshz(ix,ny,1))                                !sv 09.01.01
*    ..compute cst, vte
        tpz = 0
        trz = 0
        do is = 0, ns-1
         if (is_neutral(is)) cycle
         tpz = tpz+na(ix,ny,is)*(rza(ix,ny,is)*te(ix,ny)+ti(ix,ny))
         trz = trz+am(is)*mp*na(ix,ny,is)
        enddo
        cst = sqrt(tpz/trz)
        vte = sqrt(te(ix,ny)/me)
*    ..compute Bohm fchi, fche
        t0 = 0
        do is = 0, ns-1
         t0 = t0+rza(ix,ny,is)*fna(ix,ny,1,is)
        enddo
        fchi = qe*t0
        fche = sqrt(1/(2*pi))*qe*ne(ix,ny)*vte*s1*
     &   expu2(-qe*po(ix,ny)/te(ix,ny))
!        fchy_dia = 0.0e0_R8                                           !sv 01.07.99
!        if(facdriftm.ne.0.0e0_R8) then                                !sv 09.07.99
        if(fchy_dia_at_north .ne. 0.0e0_R8) then                       !sv 09.07.99
!         iy=ny                                                        !sv 01.07.99 { 18.05.02
!         fchy_dia =-signmf*pz(ix,iy)*bb(ix,iy,2)                      !sv 31.03.00 18.01.02
!     &      *vol(ix,iy)/gs(ix,iy,2)                                   !sv 01.07.99
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))!sv 01.07.99 } 18.05.02
          iy=ny                                                        !sv 21.10.02
          fchy_dia = fchdia(ix,iy,1)                                   !sv 21.10.02
!         iy=ny                                                        !sv 01.07.99
!         fchy_dia =-signmf*                                           !sv 09.07.99 21.04.00
!     &       (pz(ix,iy)*bb(ix,iy,2)+                                  !sv 31.03.00 18.01.02 08.06.02
!     &        pz(bottomix(ix,iy),bottomiy(ix,iy))                     !sv 31.03.00 18.01.02 08.06.02
!     &       *bb(bottomix(ix,iy),bottomiy(ix,iy),2))*0.5e0_R8         !sv 31.03.00
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))        !sv 01.07.99
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))      !sv 01.07.99
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))!sv 01.07.99
         else                                                          !sv
          fchy_dia=0.0e0_R8                                            !sv
         endif                                                         !sv
!        fchy_dia = fchy_dia*facdrift(ix,iy)                           !sv 09.07.99
!       endif                                                          !sv 09.07.99
*    ..loop over species
        do is = 0, ns-1
*     ..compute csb
         csb = sqrt((rza(ix,ny,is)*te(ix,ny)+ti(ix,ny))/(am(is)*mp))
*     ..compute particle source
*!!     cbsna(6,,) and cbsna(7,,) to be introduced with ua*pbs(,,1)
         sna0(ix,ny,0,is) = sna0(ix,ny,0,is)+
     &    cbsna(0,is,ireg)*s0
         if(mdf_fnb.ne.0 .and. .not.is_neutral(is)) then               !sv 09.11.06 {
          sna0(ix,ny,0,is) = sna0(ix,ny,0,is) + fnaPSch(ix,ny,1,is)
         endif                                                         !sv 09.11.06 }
         sna0(ix,ny,1,is) = sna0(ix,ny,1,is)+
     &    cbsna(1,is,ireg)*s0+
     &    cbsna(2,is,ireg)*csb*s0+
     &    cbsna(3,is,ireg)*csb*s1+
     &    cbsna(4,is,ireg)*cst*s0+
     &    cbsna(5,is,ireg)*cst*s1
*    ..drifts to the wall                                              !sv 02.03.05 {
         if(idrifts_to_the_wall .ne. 0) then
          if(.not.is_neutral(is) .and.
     &        cbsna(2,is,ireg).ne.0.0e0_R8) then
           if(mdf_fnb.eq.0) then                                       !sv 09.11.06 {
            sna0(ix,ny,1,is) = sna0(ix,ny,1,is)-
     *       max(0.0e0_R8,s0*(vadia(ix,ny,1,is)+vaecrb(ix,ny,1,is)))
           else
            sna0(ix,ny,1,is) = sna0(ix,ny,1,is)-                       !lk 12.12.07 {
     *      max(0.0e0_R8,
     8         (vol(ix,ny)+vol(bottomix(ix,ny),bottomiy(ix,ny)))/
     9         (hy1(ix,ny)+hy1(bottomix(ix,ny),bottomiy(ix,ny)))       !sv 07.05.08
     &         *vaecrb(ix,ny,1,is))                                    !lk 12.12.07 }
           endif                                                       !sv 09.11.06 }
          endif
         endif                                                         !sv 02.03.05 {
*     ..compute momentum source
*       (in temporary array, see below for assignment to smo)
         smo0(ix,ny,0,is) = smo0(ix,ny,0,is)+cbsmo(0,is,ireg)*s0hz     !sv 09.01.01 01.07.05
         smo0(ix,ny,1,is) = smo0(ix,ny,1,is)
         smo0(ix,ny,2,is) = smo0(ix,ny,2,is)+                          !sv 01.07.05
     &    cbsmo(1,is,ireg)*csb*pbshz(ix,ny,1)+                         !sv 09.01.01
     &    cbsmo(2,is,ireg)*cst*pbshz(ix,ny,1)+                         !sv 09.01.01
     &    cbsmo(5,is,ireg)*max(0.0e0_R8,-fna_fcor(ix,ny,1,is))/
     &                                       na(ix,ny,is)              !sv 25.02.00
         smo0(ix,ny,3,is) = smo0(ix,ny,3,is)+                          !sv 01.07.05
     &    cbsmo(3,is,ireg)*s0hz+                                       !sv 09.01.01
     &    cbsmo(4,is,ireg)*s1hz+                                       !sv 09.01.01
     &    cbsmo(6,is,ireg)                                             !sv 14.10.02
     &               *max(0.0e0_R8,fna_fcor(ix,ny,1,is))/na(ix,ny,is)  !sv 25.02.00
*     ..compute atom heat source
         t0 = na(ix,ny,is)/ni(ix,ny,0)
         shi0(ix,ny,0) = shi0(ix,ny,0)+
     &     cbshi(0,is,ireg)*t0*s0
     &    +cbshi(6,is,ireg)*max(0.0e0_R8,-fna(ix,ny,1,is))
         shi0(ix,ny,1) = shi0(ix,ny,1)+
     &     cbshi(1,is,ireg)*t0*s0
         shi0(ix,ny,3) = shi0(ix,ny,3)+
     &     cbshi(2,is,ireg)*t0*csb*s0+
     &     cbshi(3,is,ireg)*t0*csb*s1+
     &     cbshi(4,is,ireg)*t0*cst*s0+
     &     cbshi(5,is,ireg)*t0*cst*s1
!     &   +cbshi(7,is,ireg)*max(0.0e0_R8,fna(ix,ny,1,is))/ni(ix,ny,0)  !sv 28.10.02
        enddo
        shi0(ix,ny,1) = shi0(ix,ny,1)+                                 !sv 28.10.02 {
     &    cbshi(7,0,ireg)*max(0.0e0_R8,floi_noc(ix,ny,1))              !sv 28.10.02 }
*    ..compute electron heat source
        she0(ix,ny,0) = she0(ix,ny,0)+
     &    cbshe(0,ireg)*s0+
     &    cbshe(6,ireg)*max(0.0e0_R8,-fne(ix,ny,1))
        she0(ix,ny,1) = she0(ix,ny,1)+
     &    cbshe(1,ireg)*s0
        she0(ix,ny,3) = she0(ix,ny,3)+
     &    cbshe(2,ireg)*vte*s0+
     &    cbshe(3,ireg)*vte*s1+
     &    cbshe(4,ireg)*cst*s0+
     &    cbshe(5,ireg)*cst*s1
!     &  +cbshe(7,ireg)*max(0.0e0_R8,fne(ix,ny,1))/ne(ix,ny)           !sv 28.10.02
        she0(ix,ny,1) = she0(ix,ny,1)+                                 !sv 28.10.02 {
     &    cbshe(7,ireg)*max(0.0e0_R8,floe_noc(ix,ny,1))                !sv 28.10.02 }
*    ..compute charge source
        sch0(ix,ny,2) = sch0(ix,ny,2)+
     &    cbsch(0,ireg)*qe*vte*s0+
     &    cbsch(2,ireg)*qe*vte*s1+
     &    cbsch(4,ireg)*qe*cst*s1
        sch0(ix,ny,3) = sch0(ix,ny,3)+
     &    cbsch(1,ireg)*qe**2*vte*s0/te(ix,ny)+
     &    cbsch(3,ireg)*qe**2*vte*s1/te(ix,ny)+
     &    cbsch(5,ireg)*qe**2*cst*s1/te(ix,ny)
        t0 = max(cbsch(6,ireg)*fchi,cbsch(7,ireg)*fche)*qe/te(ix,ny)
!       sch0(ix,ny,0) = sch0(ix,ny,0)-fchy_dia*y_dia                   !sv  minus 01.07.99 25.04.00
        if(fchy_dia_at_north .ne. 0.0e0_R8) then                       !sv 21.10.02
            sch0(ix,ny,0) = sch0(ix,ny,0) - fchdia(ix,ny,1)*y_dia      !sv 21.10.02
        endif                                                          !sv 21.10.02
        if(fchy_in_at_north .ne. 0.0e0_R8) then                        !sv 21.10.02
            sch0(ix,ny,0) = sch0(ix,ny,0) - fchin(ix,ny,1)             !sv 21.10.02
        endif                                                          !sv 21.10.02
        if(fchy_vispar_at_north .ne. 0.0e0_R8) then                    !sv 21.10.02
          sch0(ix,ny,0) = sch0(ix,ny,0) - fchvispar(ix,ny,1)           !sv 21.10.02
        endif                                                          !sv 21.10.02
        if(fchy_visq_at_north .ne. 0.0e0_R8) then                      !sv 21.10.02
          sch0(ix,ny,0) = sch0(ix,ny,0) - fchvisq(ix,ny,1)             !sv 21.10.02
        endif                                                          !sv 21.10.02
        if(fchy_visper_at_north .ne. 0.0e0_R8) then                    !sv 21.10.02
          sch0(ix,ny,0) = sch0(ix,ny,0) - fchvisper(ix,ny,1)           !sv 21.10.02
        endif                                                          !sv 21.10.02
        if(fchy_inert_at_north .ne. 0.0e0_R8) then                     !sv 21.10.02
          sch0(ix,ny,0) = sch0(ix,ny,0) - fchinert(ix,ny,1)            !sv 21.10.02
        endif                                                          !sv 21.10.02
        sch0(ix,ny,0) = sch0(ix,ny,0)+
     &   cbsch(7,ireg)*fche-cbsch(6,ireg)*fchi+t0*po(ix,ny)
        sch0(ix,ny,1) = sch0(ix,ny,1)-t0
       enddo
*   ..compute momentum source

* ..compute conditions at west boundary
       do ix = -1,nx                                       !sv 26.11.02
        if(leftix(ix,-1).eq.-2) then                       !sv 26.11.02
*   ..initialize ireg
         ireg = cbirwe
*   ..loop over boundary cells
         do iy = -1, ny
*    ..update ireg
    3     if (cbrbrk(ireg)*ny.lt.iy+0.5e0_R8) then
           call xertst (ireg.lt.cbirwe+cbnrwe-1,
     &     'error in computing region on west boundary')
           ireg = ireg+1
           goto 3
          endif
*    ..compute s0, s1
          s0 = vol(ix,iy)/hx(ix,iy)                                      !sv 26.11.02
          s1 = abs(pbs(rightix(ix,iy),rightiy(ix,iy),0))                 !sv 26.11.02
*    ..compute s0hz, s1hz                                                !sv 09.01.01
          s0hz = hz(ix,iy)*vol(ix,iy)/hx(ix,iy)                          !sv 09.01.01 26.11.02
          s1hz = abs(pbshz(rightix(ix,iy),rightiy(ix,iy),0))             !sv 09.01.01 26.11.02
*    ..compute cst, vte
          tpz = 0
          trz = 0
          do is = 0, ns-1
           if (is_neutral(is)) cycle
           tpz = tpz+na(ix,iy,is)*(rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))    !sv 26.11.02
           trz = trz+am(is)*mp*na(ix,iy,is)                              !sv 26.11.02
          enddo
          cst = sqrt(tpz/trz)
          vte = sqrt(te(ix,iy)/me)                                       !sv 26.11.02
*    ..compute Bohm fchi, fche
          wrong_flow=.false.
          if(istyle_fchi.eq.0) then                                      !sv 14.02.05
           t0 = 0.0e0_R8
           do is = 0, ns-1
            t0=t0-rza(rightix(ix,iy),rightiy(ix,iy),is)*
     .            fna(rightix(ix,iy),rightiy(ix,iy),0,is)                !sv 26.11.02
           enddo
          else                                                           !sv 14.02.05 {
*    ..set boundary condition in explicit form (the same form as fna)
           t0 = 0.0e0_R8
           do is = 0, ns-1
            t0 = t0-rza(ix,iy,is)*cbsna(7,is,ireg)*cst*s1*na(ix,iy,is)
           enddo
          endif                                                          !sv 14.02.05 }
crfs eventually include here diamagnetic drifts
crfs  (b_phi/b^2 dpi/dr)
          fchi = qe*t0
!       if(t0.lt.0) then
!        wrong_flow=.true.
!        write(*,'(a,i3,a,i3,a,es9.2)')
!     &    ' West. Wrong flow at ix=',ix,' iy=',iy,' fchi=',fchi
!        do is = 0, ns-1
!         write(*,"(4x,' is=',i2,' fnax(ix,iy,is)=',es9.2)")
!     &      is, fna(rightix(ix,iy),rightiy(ix,iy),0,is)                 !sv 26.11.02
!        enddo                                                           !sv 20.02.03 }
c         fchi=0.0e0_R8                            ! was commented
!       endif
          if(ix.eq.-1) then                                              !sv 16.11.04 {
*  .. set the potential at inner lower plate for DND case with respect to outer plate
           fche = sqrt(1/(2*pi))*qe*ne(ix,iy)*vte*s1*
     &       expu2(-qe*(po(ix,iy)-po_in_lower)/te(ix,iy))
          else
           fche = sqrt(1/(2*pi))*qe*ne(ix,iy)*vte*s1*
     &       expu2(-qe*po(ix,iy)/te(ix,iy))
          endif                                                          !sv 16.11.04 }
*    ..loop over species
          do is = 0, ns-1
*     ..compute csb
           csb = sqrt((rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))/(am(is)*mp))   !sv 26.11.02
*     .. account for diamagnetic drifts
           csbm = max(csb,0.1e0_R8*csb+uadia(0,iy,0,is)*    !signmf*     !sv 17.02.04 {  11.01.13
     &                  bb(ix,iy,2)/bb(ix,iy,0))
           cstm = max(cst,0.1e0_R8*cst+uadia(0,iy,0,is)*    !signmf*                 11.01.13
     &                 bb(ix,iy,2)/bb(ix,iy,0))                          !sv 17.02.04 }
           if(istyle_cs.eq.0) then                                       !sv 25.04.08 {
            vwest=(-cst)*pbshz(rightix(ix,iy),rightiy(ix,iy),0)+         !sv 17.02.04 19.02.08
     &      hz(ix,iy)*vol(rightix(ix,iy),rightiy(ix,iy))/
     &                 hx(rightix(ix,iy),rightiy(ix,iy))      !*signmf*bb(ix,iy,2)/bb(ix,iy,3)   !sv 05.02.05
     &       *(-vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is))*cbc          !sv 17.02.04 19.02.08
!        vwest=(-cst)*pbshz(rightix(ix,iy),rightiy(ix,iy),0)+            !sv 17.02.04 {
!     &    hz(ix,iy)*vol(ix,iy)/hx(ix,iy)  !*signmf*bb(ix,iy,2)/bb(ix,iy,3)   !sv 05.02.05
 !    &    *(-vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is))*cbc            !sv 17.02.04 }
           else
            vwest=(-csb)*pbshz(rightix(ix,iy),rightiy(ix,iy),0)+         !sv 17.02.04 19.02.08
     &       hz(ix,iy)*vol(rightix(ix,iy),rightiy(ix,iy))/
     &                  hx(rightix(ix,iy),rightiy(ix,iy))
     &       *(-vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is))*cbc          !sv 17.02.04 19.02.08
           endif                                                         !sv 25.04.08 }
           if(ivminus.eq.0) then                                         !sv 24.03.05
             vminus= vwest                                               !sv 28.01.05
           else
             vminus=(vwest-abs(vwest))*0.5e0_R8                          !sv 28.01.05
           endif                                                         !sv 24.03.05
*     ..compute particle source
           sna0(ix,iy,0,is) = sna0(ix,iy,0,is)+                          !sv 26.11.02
     &      cbsna(0,is,ireg)*s0+
     &      cbsna(6,is,ireg)*max(0.0e0_R8,ua(ix,iy,is)                   !sv 26.11.02 01.07.05
     &          *pbs(rightix(ix,iy),rightiy(ix,iy),0))                   !sv 26.11.02
           if(mdf_fnb .ne. 0) then                                       !sv 13.10.06 {
              sna0(ix,iy,0,is) = sna0(ix,iy,0,is)-
     &                     fnaPSch(rightix(ix,iy),rightiy(ix,iy),0,is)
           endif                                                         !sv 13.10.06 }
           if(ubig0 .gt. 0.0e0_R8)                                       !iyv 22.04.08 {
     &       sna0(ix,iy,0,is) = sna0(ix,iy,0,is) - cbsna(7,is,ireg)*
     &       ubig0*s0*na(rightix(ix,iy),rightiy(ix,iy),is)
           sna0(ix,iy,1,is) = sna0(ix,iy,1,is)+                          !sv 26.11.02
     &      cbsna(1,is,ireg)*s0+
     &      cbsna(2,is,ireg)*csbm*s0+
     &      cbsna(3,is,ireg)*csbm*s1+
     &      cbsna(4,is,ireg)*cstm*s0+
     &      cbsna(5,is,ireg)*cstm*s1
!    &      cbsna(7,is,ireg)*max(0.0e0_R8,-ua(0,iy,is)*pbs(0,iy,0))
           if( ubig0 .gt. 0.0e0_R8) then                                 !iyv 22.04.08 {
             sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &           cbsna(7,is,ireg) * ubig0 * s0
           else
             if(istyle_cs.eq.0) then                                     !sv 25.04.08 {
               sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &           cbsna(7,is,ireg)*cst*s1 ! pbs(rightix(ix,iy),rightiy(ix,iy),0)   !sv 02.11.99 26.11.02 14.02.05
             else
               sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &           cbsna(7,is,ireg)*csb*s1 ! pbs(rightix(ix,iy),rightiy(ix,iy),0)   !sv 02.11.99 26.11.02 14.02.05
             endif                                                       !sv  25.04.08
           endif                                                         !iyv 22.04.08 }
!       if(wrong_flow.gt.0) then
!         write(*,*) 'b2stbc_west density set to zero:iy,fchi',iy,fchi
crfs force density to 0.0e0_R8 if flux has wrong direction
c         sna0(-1,iy,0,1) = 0.0e0_R8
c         sna0(-1,iy,1,1) = -1.0d60
!       endif
*     ..compute momentum source
!       if (0.0e0_R8 .le. vwest) then                                    !sv 02.11.99
!         write(*,'(a,i3,a,i3,a,i3)') ' west: ix =', rightix(ix,iy),     !sv 25.11.04
!     *       '  iy =',rightiy(ix,iy), '  Switch in bc for Upar, is =',is!sv 02.11.99 20.02.03
!         write(*,*) '     cst*pbshz(ix,iy,0) = ',
!     &                     cst*pbshz(rightix(ix,iy),rightiy(ix,iy),0)
!          write(*,*) '       hz*vol/hx*Bz/B*(-VaExBx) = ',
!     &           hz(ix,iy)*vol(ix,iy)/hx(ix,iy) ! *signmf*bb(ix,iy,2)/bb(ix,iy,3) !sv 05.02.05
!     &           *(-vaecrb(rightix(ix,iy),rightiy(ix,iy),0,is))*cbc
!!         vwest = 0.0e0_R8                                              !sv 02.11.99
!        endif
           smo0(ix,iy,0,is) = smo0(ix,iy,0,is)+  !sv 26.11.02
     &      cbsmo(0,is,ireg)*s0hz                                        !sv 09.01.01
           smo0(ix,iy,2,is) = smo0(ix,iy,2,is)+ cbsmo(1,is,ireg)*        !sv 26.11.02
     &      (-csbm)*pbshz(rightix(ix,iy),rightiy(ix,iy),0)+              !sv 09.01.01 26.11.02
     &      cbsmo(2,is,ireg)*vminus+                                     !sv 02.11.99 28.01.05
     &      cbsmo(5,is,ireg)* max(0.0e0_R8,
     &       fna_fcor(rightix(ix,iy),rightiy(ix,iy),0,is))/
     &       na(ix,iy,is)                                                !sv 25.02.00 26.11.02 01.07.05
           smo0(ix,iy,3,is) = smo0(ix,iy,3,is)+                          !sv 26.11.02
     &      cbsmo(3,is,ireg)*s0hz+                                       !sv 09.01.01
     &      cbsmo(4,is,ireg)*s1hz+                                       !sv 09.01.01
     &      cbsmo(6,is,ireg)*max(0.0e0_R8,
     &     -fna_fcor(rightix(ix,iy),rightiy(ix,iy),0,is))/na(ix,iy,is)   !sv 25.02.00 26.11.02 01.07.05
*     ..compute atom heat source
           t0 = na(ix,iy,is)/ni(ix,iy,0)                                 !sv 26.11.02
           shi0(ix,iy,0) = shi0(ix,iy,0)+
     &      cbshi(0,is,ireg)*t0*s0+cbshi(6,is,ireg)*
     &      max(0.0e0_R8,fna(rightix(ix,iy),rightiy(ix,iy),0,is))        !sv 26.11.02
           shi0(ix,iy,1) = shi0(ix,iy,1)+                                !sv 26.11.02
     &      cbshi(1,is,ireg)*t0*s0
           shi0(ix,iy,3) = shi0(ix,iy,3)+                                !sv 26.11.02
     &      cbshi(2,is,ireg)*t0*csbm*s0+
     &      cbshi(3,is,ireg)*t0*csbm*s1+
     &      cbshi(4,is,ireg)*t0*cstm*s0+
     &      cbshi(5,is,ireg)*t0*cstm*s1+
     &      cbshi(7,is,ireg)*t0*cst*pbs(rightix(ix,iy),rightiy(ix,iy),0) !sv 26.11.02 28.01.03
!     &    cbshi(7,is,ireg)*max(0.0e0_R8,-fna(rightix(ix,iy),rightiy(ix,iy),0,is)) !sv 26.11.02 28.01.03
!     &                                  /ni(ix,iy,0)                    !sv 26.11.02 28.01.03
          enddo
          if(mdf_fhi.ne.0) then                                          !sv 13.10.06 {
           shi0(ix,iy,0) = shi0(ix,iy,0) -
     .      fhiPSch(rightix(ix,iy),rightiy(ix,iy),0)
          endif                                                          !sv 13.10.06 }
*    ..compute electron heat source
          she0(ix,iy,0) = she0(ix,iy,0)+                                 !sv 26.11.02
     &     cbshe(0,ireg)*s0+cbshe(6,ireg)*
     &     max(0.0e0_R8,fne(rightix(ix,iy),rightiy(ix,iy),0))            !sv 26.11.02
          she0(ix,iy,1) = she0(ix,iy,1)+                                 !sv 26.11.02
     &     cbshe(1,ireg)*s0
          she0(ix,iy,3) = she0(ix,iy,3)+                                 !sv 26.11.02
     &     cbshe(2,ireg)*vte*s0+
     &     cbshe(3,ireg)*vte*s1+
     &     cbshe(4,ireg)*cstm*s0+
     &     cbshe(5,ireg)*cstm*s1+                                        !sv 11.02.00
!     &   cbshe(7,ireg)*((2.0e0_R8-cbsch(7,ireg))/cbsch(7,ireg)+3.0e0_R8)!sv 11.02.00 !iyv 19.05.08 !sv 03.10.08
     &     cbshe(7,ireg)*((2.0e0_R8-1.0e0_R8)/1.0e0_R8+3.0e0_R8)         !sv 11.02.00 !iyv 19.05.08 !sv 03.10.08
!     &    *max(0.0e0_R8,-fne(rightix(ix,iy),rightiy(ix,iy),0))/ne(ix,iy)!sv 11.02.00 26.11.02 25.12.06
     &    *max(0.0e0_R8, (pbs(rightix(ix,iy),rightiy(ix,iy),0)*cst       !sv 25.12.06
     &     +fch_p(rightix(ix,iy),rightiy(ix,iy),0)/qe/ne(ix,iy)))        !sv 25.12.06
!         write(*,"('pbs*cst=',es10.3,' fch_p/qe/ne=',es10.3)")
!     &             pbs(rightix(ix,iy),rightiy(ix,iy),0)*cst,
!     &           fch_p(rightix(ix,iy),rightiy(ix,iy),0)/qe/ne(ix,iy)
          if(mdf_fhe.ne.0) then                                          !sv 01.06.06 {
           she0(ix,iy,0) = she0(ix,iy,0) -
     .      fhePSch(rightix(ix,iy),rightiy(ix,iy),0)
          endif                                                          !sv 01.06.06 }
*    ..compute charge source
          sch0(ix,iy,2) = sch0(ix,iy,2)+                                 !sv 26.11.02
     &     cbsch(0,ireg)*qe*vte*s0+
     &     cbsch(2,ireg)*qe*vte*s1+
     &     cbsch(4,ireg)*qe*cstm*s1
          sch0(ix,iy,3) = sch0(ix,iy,3)+                                 !sv 26.11.02 {
     &     cbsch(1,ireg)*qe**2*vte*s0/te(ix,iy)+
     &     cbsch(3,ireg)*qe**2*vte*s1/te(ix,iy)+
     &     cbsch(5,ireg)*qe**2*cstm*s1/te(ix,iy)                         !sv 26.11.02 }
          t0 = max(cbsch(6,ireg)*fchi,cbsch(7,ireg)*fche)*qe/te(ix,iy)   !sv 26.11.02
          sch0(ix,iy,0) = sch0(ix,iy,0)+                                 !sv 26.11.02
     &     cbsch(7,ireg)*fche-cbsch(6,ireg)*fchi+t0*po(ix,iy)            !sv 26.11.02
          sch0(ix,iy,1) = sch0(ix,iy,1)-t0                               !sv 26.11.02
         enddo
        endif                                                            !sv 26.11.02
       enddo                                                             !sv 26.11.02

* ..compute conditions at east boundary
       do ix=-1,nx                                                       !sv 26.11.02
        if(rightix(ix,-1).eq.nx+1) then                                  !sv 26.11.02
*   ..initialize ireg
         ireg = cbirea
*   ..loop over boundary cells
         do iy = -1, ny
*    ..update ireg
    4     if (cbrbrk(ireg)*ny.lt.iy+0.5e0_R8) then
           call xertst (ireg.lt.cbirea+cbnrea-1,
     &      'error in computing region on east boundary')
           ireg = ireg+1
           goto 4
          endif
*    ..compute s0, s1
          s0 = vol(ix,iy)/hx(ix,iy)                                      !sv 26.11.02
          s1 = abs(pbs(ix,iy,0))                                         !sv 26.11.02
*    ..compute s0hz, s1hz                                                !sv 09.01.01
          s0hz = hz(ix,iy)*vol(ix,iy)/hx(ix,iy)                          !sv 09.01.01 26.11.02
          s1hz = abs(pbshz(ix,iy,0))                                     !sv 09.01.01 26.11.02
*    ..compute cst, vte
          tpz = 0
          trz = 0
          do is = 0, ns-1
           if (is_neutral(is)) cycle
           tpz = tpz+na(ix,iy,is)*(rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))    !sv 26.11.02
           trz = trz+am(is)*mp*na(ix,iy,is)                              !sv 26.11.02
          enddo
          cst = sqrt(tpz/trz)
          vte = sqrt(te(ix,iy)/me)                                       !sv 26.11.02
*    ..compute Bohm fchi, fche
          wrong_flow=.false.
          if(istyle_fchi.eq.0) then                                      !sv 14.02.05
           t0 = 0
           do is = 0, ns-1
            t0 = t0+rza(ix,iy,is)*fna(ix,iy,0,is)                        !sv 26.11.02
           enddo
          else                                                           !sv 14.02.05 {
*    ..set boundary condition in explicit form (the same form as fna)
           t0 = 0.0e0_R8
           do is = 0, ns-1
            t0 = t0-rza(ix,iy,is)*cbsna(7,is,ireg)*cst*s1*na(ix,iy,is)
           enddo
          endif                                                          !sv 14.02.05 }
crfs eventually include here diamagnetic drifts
crfs  (b_phi/b^2 dpi/dr)
          fchi = qe*t0
!       if(t0.lt.0) then
!         wrong_flow=.true.
!         write(*,'(a,i3,a,i3,a,es9.2)')
!     &    ' East. Wrong flow at ix=',ix,' iy=',iy,' fchi=',fchi
!         do is = 0, ns-1
!           write(*,"(4x,' is=',i2,' fnax(ix,iy,is)=',es9.2)")
!     &           is, fna(ix,iy,0,is)                                    !sv 26.11.02
!         enddo                                                          !sv 20.02.03 }
c         fchi=0.0e0_R8                      !sv was commented
!       endif
         if(ix.ne.nx) then                                              !sv 16.11.04 {
*  .. set the potential at inner upper plate for DND case with respect to outer plate
           fche = sqrt(1/(2*pi))*qe*ne(ix,iy)*vte*s1*
     &       expu2(-qe*(po(ix,iy)-po_in_upper)/te(ix,iy))
          else
           fche = sqrt(1/(2*pi))*qe*ne(ix,iy)*vte*s1*
     &      expu2(-qe*po(ix,iy)/te(ix,iy))
          endif                                                          !sv 16.11.04 }
*    ..loop over species
          do is = 0, ns-1
*     ..compute csb
           csb = sqrt((rza(ix,iy,is)*te(ix,iy)+ti(ix,iy))/(am(is)*mp))   !sv 26.11.02
*     .. account for diamagnetic drifts
           csbm = max(csb,0.1e0_R8*csb-uadia(ix,iy,0,is)*                !sv 17.02.04 {
     1           bb(ix,iy,2)/bb(ix,iy,0))
           cstm = max(cst,0.1e0_R8*cst-uadia(ix,iy,0,is)*
     1           bb(ix,iy,2)/bb(ix,iy,0))                                !sv 17.02.04 }
           if(istyle_cs.eq.0) then                                       !sv 25.04.08 {
            veast = cst*pbshz(ix,iy,0)+                                  !sv 17.02.04 19.02.08 25.04.08 {
     &       hz(ix,iy)*vol(leftix(ix,iy),iy)/hx(leftix(ix,iy),iy)     !*signmf*bb(ix,iy,2)/bb(ix,iy,3) !sv 05.02.05
     &           *(-vaecrb(ix,iy,0,is))*cbc                              !sv 17.02.04 19.02.08
!        veast = cst*pbshz(ix,iy,0)+                                     !sv 17.02.04 {
!     &       hz(ix,iy)*vol(ix,iy)/hx(ix,iy) !*signmf*bb(ix,iy,2)/bb(ix,iy,3) !sv 05.02.05
!     &           *(-vaecrb(ix,iy,0,is))*cbc                             !sv 17.02.04 }
           else
            veast = csb*pbshz(ix,iy,0)+                                  !sv 17.02.04 19.02.08
     &       hz(ix,iy)*vol(leftix(ix,iy),iy)/hx(leftix(ix,iy),iy)
     &           *(-vaecrb(ix,iy,0,is))*cbc                              !sv 17.02.04 19.02.08
           endif                                                         !sv 25.04.08 }
           if(ivplus.eq.0) then                                          !sv 24.03.05
             vplus=veast                                                 !sv 28.01.05
           else
             vplus=(veast+abs(veast))*0.5e0_R8                           !sv 28.01.05
           endif                                                         !sv 24.03.05
!        v0 =vol(ix,iy)/hx(ix,iy)*bb(ix,iy,2)/bb(ix,iy,3)                !sv 18.07.03 {
!     &      *vaecrb(ix,iy,0,is)*cbc*signmf
*     ..compute particle source
           sna0(ix,iy,0,is) = sna0(ix,iy,0,is)+                          !sv 26.11.02
     &      cbsna(0,is,ireg)*s0+
     &      cbsna(6,is,ireg)*max(0.0e0_R8,-ua(ix,iy,is)*pbs(ix,iy,0))    !sv 26.11.02
           if(mdf_fnb .ne. 0) then                                       !sv 13.10.06 {
            sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + fnaPSch(ix,iy,0,is)
           endif                                                         !sv 13.10.06 }
           if(ubig0 .gt. 0.0e0_R8)                                       !iyv 22.04.08
     &      sna0(ix,iy,0,is) = sna0(ix,iy,0,is) -
     &       cbsna(7,is,ireg)*ubig0*s0*na(leftix(ix,iy),iy,is)
           sna0(ix,iy,1,is) = sna0(ix,iy,1,is)+                          !sv 26.11.02
     &      cbsna(1,is,ireg)*s0+
     &      cbsna(2,is,ireg)*csbm*s0+
     &      cbsna(3,is,ireg)*csbm*s1+
     &      cbsna(4,is,ireg)*cstm*s0+
     &      cbsna(5,is,ireg)*cstm*s1
!     &   +cbsna(7,is,ireg)*max(0.0e0_R8,ua(nx,iy,is)*pbs(nx,iy,0))
           if(ubig0 .gt. 0.0e0_R8) then                                  !iyv 22.04.08 {
             sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &               cbsna(7,is,ireg)*ubig0*s0
           else
             if(istyle_cs.eq.0) then                                     !sv 25.04.08 {
               sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &          cbsna(7,is,ireg)*cst*s1     ! pbs(ix,iy,0)               !sv 14.02.05
            else
              sna0(ix,iy,1,is) = sna0(ix,iy,1,is) +
     &         cbsna(7,is,ireg)*csb*s1     ! pbs(ix,iy,0)                !sv 14.02.05
            endif                                                        !sv 25.04.08
          endif                                                          !iyv 22.04.08 }
!       if(wrong_flow) then
!         write(*,*) 'b2stbc_east density set to zero:iy,fchi',iy,fchi
crfs force density to 0.0e0_R8 if flux has wrong direction
c         sna0(nx,iy,0,1) = 0.0e0_R8
c         sna0(nx,iy,1,1) = -1.0d160
!       endif
*     ..compute momentum source
!       if (0.0e0_R8 .ge. veast) then                                    !sv 02.11.99
!         write(*,'(a,i3,a,i3,a,i3)') ' east: ix =', ix,'  iy =',iy,     !sv 25.11.04
!     *    '  Switch in bc for Upar, is =',is                            !sv 02.11.99 20.02.03
!         write(*,*) '     cst*pbshz(ix,iy,0) = ',cst*pbshz(ix,iy,0)
!         write(*,*) '       hz*vol/hx*Bz/B*(-VaExBx) = ',
!     &           hz(ix,iy)*vol(ix,iy)/hx(ix,iy) ! *signmf*bb(ix,iy,2)/bb(ix,iy,3) !sv 05.02.05
!     &           *(-vaecrb(ix,iy,0,is))*cbc
!!        veast = 0.0e0_R8                                               !sv 02.11.99
!        endif
           smo0(ix,iy,0,is) = smo0(ix,iy,0,is)+                          !sv 26.11.02
     &      cbsmo(0,is,ireg)*s0hz                                        !sv 09.01.01
           smo0(ix,iy,2,is) = smo0(ix,iy,2,is)+                          !sv 26.11.02
     &      cbsmo(1,is,ireg)*csbm*pbshz(ix,iy,0)+                        !sv 09.01.01 26.11.02
     &      cbsmo(2,is,ireg)*vplus+                                      !sv 02.11.99 28.01.05
     &      cbsmo(5,is,ireg)*max(0.0e0_R8,-fna_fcor(ix,iy,0,is))*        !sv 25.02.00 26.11.02
     &       2/(na(leftix(ix,iy),iy,is)+na(ix,iy,is))                    !sv 26.11.02
           smo0(ix,iy,3,is) = smo0(ix,iy,3,is)+                          !sv 26.11.02
     &      cbsmo(3,is,ireg)*s0hz+                                       !sv 09.01.01
     &      cbsmo(4,is,ireg)*s1hz+                                       !sv 09.01.01
     &      cbsmo(6,is,ireg)*max(0.0e0_R8,fna_fcor(ix,iy,0,is))*         !sv 25.02.00 26.11.02
     &       2/(na(leftix(ix,iy),iy,is)+na(ix,iy,is))                    !sv 14.10.02 26.11.02
*     ..compute atom heat source
           t0 = na(ix,iy,is)/ni(ix,iy,0)                                 !sv 26.11.02
           shi0(ix,iy,0) = shi0(ix,iy,0)+                                !sv 26.11.02
     &      cbshi(0,is,ireg)*t0*s0+
     &      cbshi(6,is,ireg)*max(0.0e0_R8,-fna(ix,iy,0,is))              !sv 26.11.02
           shi0(ix,iy,1) = shi0(ix,iy,1)+                                !sv 26.11.02
     &      cbshi(1,is,ireg)*t0*s0
           shi0(ix,iy,3) = shi0(ix,iy,3)+                                !sv 26.11.02
     &      cbshi(2,is,ireg)*t0*csbm*s0+
     &      cbshi(3,is,ireg)*t0*csbm*s1+
     &      cbshi(4,is,ireg)*t0*cstm*s0+
     &      cbshi(5,is,ireg)*t0*cstm*s1+
     &      cbshi(7,is,ireg)*t0*cst*pbs(ix,iy,0)                         !sv 26.11.02 28.01.03
!     &    cbshi(7,is,ireg)*max(0.0e0_R8,fna(ix,iy,0,is))/ni(ix,iy,0)    !sv 26.11.02
          enddo
          if(mdf_fhi.ne.0) then                                          !sv 13.10.06 {
            shi0(ix,iy,0) = shi0(ix,iy,0)+fhiPSch(ix,iy,0)
          endif                                                          !sv 13.10.06 }
*    ..compute electron heat source
          she0(ix,iy,0) = she0(ix,iy,0)+                                 !sv 26.11.02
     &     cbshe(0,ireg)*s0+
     &     cbshe(6,ireg)*max(0.0e0_R8,-fne(ix,iy,0))                     !sv 26.11.02
          she0(ix,iy,1) = she0(ix,iy,1)+                                 !sv 26.11.02
     &     cbshe(1,ireg)*s0
          she0(ix,iy,3) = she0(ix,iy,3)+                                 !sv 26.11.02
     &     cbshe(2,ireg)*vte*s0+
     &     cbshe(3,ireg)*vte*s1+
     &     cbshe(4,ireg)*cstm*s0+
     &     cbshe(5,ireg)*cstm*s1+                                        !sv 11.02.00
!     &   cbshe(7,ireg)*((2.0e0_R8-cbsch(7,ireg))/cbsch(7,ireg)+3.0e0_R8)!sv 11.02.00 !iyv 19.05.08 !sv 03.10.08
     &     cbshe(7,ireg)*((2.0e0_R8-1.0e0_R8)/1.0e0_R8+3.0e0_R8)         !sv 11.02.00 !iyv 19.05.08 !sv 03.10.08
!     &                *max(0.0e0_R8,fne(ix,iy,0))/ne(ix,iy)             !sv 11.02.00 26.11.02 25.12.06
     &    *max(0.0e0_R8,
     &                (pbs(ix,iy,0)*cst-fch_p(ix,iy,0)/qe/ne(ix,iy)))    !sv 25.12.06
          if(mdf_fhe.ne.0) then                                          !sv 13.10.06 {
            she0(ix,iy,0) = she0(ix,iy,0)+fhePSch(ix,iy,0)               !sv 09.11.08
          endif                                                          !sv 13.10.06 }
*    ..compute charge source
          sch0(ix,iy,2) = sch0(ix,iy,2)+                                 !sv 26.11.02
     &     cbsch(0,ireg)*qe*vte*s0+
     &     cbsch(2,ireg)*qe*vte*s1+
     &     cbsch(4,ireg)*qe*cstm*s1
          sch0(ix,iy,3) = sch0(ix,iy,3)+                                 !sv 26.11.02 {
     &     cbsch(1,ireg)*qe**2*vte*s0/te(ix,iy)+
     &     cbsch(3,ireg)*qe**2*vte*s1/te(ix,iy)+
     &     cbsch(5,ireg)*qe**2*cstm*s1/te(ix,iy)                         !sv 26.11.02 }
!         if (po(ix,iy).ge.0.0e0_R8) then                                !sv 15.09.98 26.11.02
          t0 = max(cbsch(6,ireg)*fchi,cbsch(7,ireg)*fche)*qe/te(ix,iy)   !sv 26.11.02
!         else                                                           !sv 15.09.98
!           t0 = 0.0e0_R8                                                !sv 15.09.98
!         endif                                                          !sv 15.09.98
          sch0(ix,iy,0) = sch0(ix,iy,0)+                                 !sv 26.11.02
     &     cbsch(7,ireg)*fche-cbsch(6,ireg)*fchi+t0*po(ix,iy)            !sv 26.11.02
          sch0(ix,iy,1) = sch0(ix,iy,1)-t0                               !sv 26.11.02
         enddo
        endif                                                            !sv 26.11.02
       enddo                                                             !sv 26.11.02
      endif
* ..compute conditions on the cuts
*     (no flow)
*   ..compute suitable large coefficient
      t0 = 1.0e20_R8*b2ssum(n2,vol,1)
*     (the dimension of t0 is m**3/s)
*   ..loop over species
      do is = 0, ns-1
*    ..set momentum source
       call b2xxgs (nx, ny, 0.0e0_R8, smo0(-1,-1,0,is), 1)
       call b2xxgs (nx, ny, 0.0e0_R8, smo0(-1,-1,1,is), 1)
       call b2xxgs (nx, ny, 0.0e0_R8, smo0(-1,-1,2,is), 1)
       call b2xxgs (nx, ny, -t0, smo0(-1,-1,3,is), 1)
      enddo

      if (iout.eq.1) then                                             !sv 11.07.05 {
       do is = 0, ns-1
         write (chns,'(i3.3)') is
         do k=0,3
           write (chk,'(i1)') k
           call my_out(70,nx,ny,smo0(-1,-1,k,is),
     &               'b2stbc_spb_smo'//chk//chns)
         enddo
         do k=0,1
           write (chk,'(i1)') k
           call my_out(70,nx,ny,sna0(-1,-1,k,is),
     &               'b2stbc_spb_sna'//chk//chns)
         enddo
       enddo
       do k=0,3
         write (chk,'(i1)') k
         call my_out(70,nx,ny,she0(-1,-1,k),'b2stbc_spb_she'//chk)
         call my_out(70,nx,ny,shi0(-1,-1,k),'b2stbc_spb_shi'//chk)
         call my_out(70,nx,ny,sch0(-1,-1,k),'b2stbc_spb_sch'//chk)
       enddo
      endif                                                           !sv 11.07.05 }

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stbc_spb

      end subroutine b2stbc_spb
