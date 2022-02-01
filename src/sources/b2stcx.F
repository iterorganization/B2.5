*-----------------------------------------------------------------------
*.specification

      subroutine b2stcx (nx, ny, ns, is0, ismain, hz, BoRiS, vol, na,   !srv 09.01.01
     &  ua, wadia, vadia, vaecrb, ti, ni, rcx, sna0, smq0, shi0)        !srv 07.07.09
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_tallies
      use b2mod_sources
      use b2mod_b2cmpa
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif
      use b2mod_plasma
     & , only : roxa
      use b2mod_balance !djm Jan2017
     & , only : b2stcx_sna0to1, b2stcx_smq0to3, b2stcx_shi0to3,
     &          balance_netcdf
      use b2mod_subsys
#ifdef WITH_FTIMINGS
      use timings
#endif     
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, is0, ismain
      real (kind=R8) ::
     *  hz(-1:nx,-1:ny), BoRiS, vol(-1:nx,-1:ny),                       !srv 09.01.01
     &  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1), rcx(-1:nx,-1:ny,0:ns-1),
     &  wadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  vadia(-1:nx,-1:ny,0:1,0:ns-1),
     &  vaecrb(-1:nx,-1:ny,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smq0(-1:nx,-1:ny,0:3,0:ns-1),
     &  shi0(-1:nx,-1:ny,0:3)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STCX provides linearised expressions for the sources of atom
*     density and atom parallel momentum, due to charge exchange.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, styl0, n2, ix, iy, is, k
      integer, save :: iout = 0, !iout_b2wdat = 0,                      !srv 30.06.08 11.09.09
     &                 iout_b2npc9 = 0, iout_b2npco = 0,
     &                 iout_b2npc7 = 0, iout_b2npc11 = 0,               !srv 02.04.18
     &                 iout_b2nph9 = 0, iout_b2npht = 0,
     &                 iout_b2npmo =0
      character chns*3, chcx*3, chk*1                                   !srv 11.09.09
      real (kind=R8) ::
     *  rg0, t0, t1, t2, du0, du1
#ifdef ASTRA
      real (kind=R8) :: wrk0(-1:nx,-1:ny)
#endif
      save ncall, styl0, rg0
*   ..procedures
      real (kind=R8) ::
     *  samax
      external xertst, ipgeti, ipgetr, samax, sfill
      external b2xvsg
*   ..initialisation
      data ncall/0/, styl0/0/, rg0/1.0_R8/
*!!   (rg0 for numerical stabilisation; needs experiments.)

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stcx')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgetr ('b2stcx_rg0', rg0)
       call xertst (0.0_R8.le.rg0, 'faulty internal parameter rg0')
       call ipgeti ('b2stcx_styl0', styl0)
       call xertst (0.le.styl0.and.styl0.lt.4,
     &   'faulty internal parameter styl0')
       call ipgeti ('b2stcx_iout', iout)                              !srv 11.09.09 {
       call ipgeti ('b2npc9_iout', iout_b2npc9)
       call ipgeti ('b2npco_iout', iout_b2npco)
       call ipgeti ('b2npc7_iout', iout_b2npc7)                       !srv 02.04.18
       call ipgeti ('b2npc11_iout', iout_b2npc11)                     !srv 02.04.18
       call ipgeti ('b2npmo_iout', iout_b2npmo)
       call ipgeti ('b2nph9_iout', iout_b2nph9)
       call ipgeti ('b2npht_iout', iout_b2npht)                       !srv 11.09.09 }
!       call ipgeti ('b2wdat_iout', iout_b2wdat)                      !iys 17.04.17
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns, ismain
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.ismain.and.ismain.lt.ns,
     * 'faulty argument ismain')
*   ..test is0
      call xertst (0.le.is0.and.is0.lt.ns-1, 'faulty argument is0')
      call xertst (lnext(is0,is0+1),
     &  'faulty argument is0: not related to species is0+1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, na, ti, ni
       call b2xvsg (n2, vol, 1, 'vol', '.gt.')
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
*    ..test rate coefficients
       do is = 1, ns-1
        if (.not.lnext(is-1,is)) then
         call xertst (samax(n2,rcx(-1,-1,is),1).eq.0.0_R8,
     &     'faulty argument rcx: nonzero for unrelated species')
        endif
       enddo
      endif

ctmp.dpc
c       write(*,*) 'DPC: b2stcx: is0 ',is0
c       do is = 1, ns-1
c        write(*,*) 'is,lnext,samax(rcx) ',
c     1   is,lnext(is-1,is),samax(n2,rcx(-1,-1,is),1)
c       enddo
ctmp.dpc

* ..compute source terms
*   ..initialise sources to 0
      call sfill (n2*2*ns, 0.0_R8, sna0, 1)
      call sfill (n2*4*ns, 0.0_R8, smq0, 1)
      call sfill (n2*4, 0.0_R8, shi0, 1)
*   ..compute sources due to charge exchange
#ifndef NO_OPENMP_B2STCX
#ifdef WITH_FTIMINGS
      call timer%start("parallel")
#endif
!$OMP PARALLEL DO default(none) private(is, ix, iy)
!$OMP& private(t0, t1, t2, du0, du1)
!$OMP& shared(ns, nx, ny, is0, BoRiS, rg0, styl0)
!$OMP& shared(rightix, leftix, topiy, bottomiy, region, hz)
!$OMP& shared(rcx, vol, sna0, na, rcxna, am, shi0, rcxhi, smq0, ua)
!$OMP& reduction(+: rcxnareg, rcxhireg, rcxmoreg)
!$OMP& shared(rcxmo, zamin, zamax, zn)
#endif
      do iy = -1, ny
      do is = 1, ns-1
       if (lnext(is-1,is)) then
*     ..particle source
        if (is-1.eq.is0) then
*        (This is a special case: charge exchange between H0 and H1,
*        for which there is no particle source.)
        else
          do ix = -1, nx
           if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &        topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &        region(ix,iy,0).ne.0) then
*        ..compute rate for is->is-1 and is0->is0+1
            t0 = rcx(ix,iy,is)*vol(ix,iy)
*        ..compute particle source
            sna0(ix,iy,0,is-1) = sna0(ix,iy,0,is-1)+
     &       t0*na(ix,iy,is0)*na(ix,iy,is)
            sna0(ix,iy,0,is0+1) = sna0(ix,iy,0,is0+1)+
     &       t0*na(ix,iy,is0)*na(ix,iy,is)
            sna0(ix,iy,1,is) = sna0(ix,iy,1,is)-
     &       t0*na(ix,iy,is0)
            sna0(ix,iy,1,is0) = sna0(ix,iy,1,is0)-
     &       t0*na(ix,iy,is)
            rcxna(ix,iy,is)=rcxna(ix,iy,is)+
     1       t0*na(ix,iy,is0)*na(ix,iy,is)
            rcxna(ix,iy,is0)=rcxna(ix,iy,is0)-
     1       t0*na(ix,iy,is0)*na(ix,iy,is)
            rcxnareg(region(ix,iy,0),is)=
     1       rcxnareg(region(ix,iy,0),is)+t0*na(ix,iy,is0)*na(ix,iy,is)
            rcxnareg(region(ix,iy,0),is0)=
     1       rcxnareg(region(ix,iy,0),is0)-t0*na(ix,iy,is0)*na(ix,iy,is)
           endif
         enddo
        endif
*     ..heat source
         do ix = -1, nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &       topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &       region(ix,iy,0).ne.0) then
*       ..compute rate for is->is-1 and is0->is0+1
           t0 = rcx(ix,iy,is)*vol(ix,iy)
*       ..compute atom heat source
           du0 = ua(ix,iy,is)-ua(ix,iy,is-1)
           du1 = ua(ix,iy,is0)-ua(ix,iy,is0+1)
           t1 = t0*na(ix,iy,is0)*na(ix,iy,is)*(mp/2.0_R8)*
     &      (((am(is)+am(is-1))/2.0_R8)*(du0**2)+
     &       ((am(is0)+am(is0+1))/2.0_R8)*(du1**2))*
     &      (1.0_R8-BoRiS)
           shi0(ix,iy,0) = shi0(ix,iy,0)+t1
           rcxhi(ix,iy,is)=rcxhi(ix,iy,is)+t1
           rcxhireg(region(ix,iy,0),is)=rcxhireg(region(ix,iy,0),is)+t1
          endif
        enddo
*     ..parallel momentum source
         do ix = -1, nx
          if(rightix(ix,iy).ne.nx+1 .and. leftix(ix,iy).ne.-2 .and.
     &        topiy(ix,iy).ne.ny+1 .and. bottomiy(ix,iy).ne.-2 .and.
     &        region(ix,iy,0).ne.0) then
*       ..compute rate for is->is-1 and is0->is0+1
            t0 = rcx(ix,iy,is)*vol(ix,iy)*na(ix,iy,is0)*na(ix,iy,is)
     &           *hz(ix,iy)                                          !srv 09.01.01
            t1 = ((am(is)+am(is-1))/2.0_R8)*mp*t0
            t2 = ((am(is0)+am(is0+1))/2.0_R8)*mp*t0
            if (styl0.eq.0) then
*       ..compute parallel momentum source
              smq0(ix,iy,0,is) = smq0(ix,iy,0,is)+
     &         rg0*t1*ua(ix,iy,is)
              smq0(ix,iy,1,is) = smq0(ix,iy,1,is)-
     &         (1.0_R8+rg0)*t1
              smq0(ix,iy,0,is-1) = smq0(ix,iy,0,is-1)+
     &         t1*(ua(ix,iy,is)+rg0*ua(ix,iy,is-1))
              smq0(ix,iy,1,is-1) = smq0(ix,iy,1,is-1)-
     &         rg0*t1
              smq0(ix,iy,0,is0+1) = smq0(ix,iy,0,is0+1)+
     &         t2*(ua(ix,iy,is0)+rg0*ua(ix,iy,is0+1))
              smq0(ix,iy,1,is0+1) = smq0(ix,iy,1,is0+1)-
     &         rg0*t2
              smq0(ix,iy,0,is0) = smq0(ix,iy,0,is0)+
     &         rg0*t2*ua(ix,iy,is0)
              smq0(ix,iy,1,is0) = smq0(ix,iy,1,is0)-
     &         (1.0_R8+rg0)*t2
            else if (styl0.eq.1) then
*       ..compute parallel momentum source
              smq0(ix,iy,2,is) = smq0(ix,iy,2,is)+
     &         rg0*t1*ua(ix,iy,is)/roxa(ix,iy,is)
              smq0(ix,iy,3,is) = smq0(ix,iy,3,is)-
     &         (1.0_R8+rg0)*t1/roxa(ix,iy,is)
              smq0(ix,iy,2,is-1) = smq0(ix,iy,2,is-1)+
     &         t1*(ua(ix,iy,is)+rg0*ua(ix,iy,is-1))/roxa(ix,iy,is-1)
              smq0(ix,iy,3,is-1) = smq0(ix,iy,3,is-1)-
     &         rg0*t1/roxa(ix,iy,is-1)
              smq0(ix,iy,2,is0+1) = smq0(ix,iy,2,is0+1)+
     &         t2*(ua(ix,iy,is0)+rg0*ua(ix,iy,is0+1))/roxa(ix,iy,is0+1)
              smq0(ix,iy,3,is0+1) = smq0(ix,iy,3,is0+1)-
     &         rg0*t2/roxa(ix,iy,is0+1)
              smq0(ix,iy,2,is0) = smq0(ix,iy,2,is0)+
     &         rg0*t2*ua(ix,iy,is0)/roxa(ix,iy,is0)
              smq0(ix,iy,3,is0) = smq0(ix,iy,3,is0)-
     &         (1.0_R8+rg0)*t2/roxa(ix,iy,is0)
            else if (styl0.eq.2) then
*       ..compute parallel momentum source
              smq0(ix,iy,0,is) = smq0(ix,iy,0,is)+
     &         rg0*t1*ua(ix,iy,is)
              smq0(ix,iy,3,is) = smq0(ix,iy,3,is)-
     &         (1.0_R8+rg0)*t1/roxa(ix,iy,is)
              smq0(ix,iy,0,is-1) = smq0(ix,iy,0,is-1)+
     &         t1*(ua(ix,iy,is)+rg0*ua(ix,iy,is-1))
              smq0(ix,iy,3,is-1) = smq0(ix,iy,3,is-1)-
     &         rg0*t1/roxa(ix,iy,is-1)
              smq0(ix,iy,0,is0+1) = smq0(ix,iy,0,is0+1)+
     &         t2*(ua(ix,iy,is0)+rg0*ua(ix,iy,is0+1))
              smq0(ix,iy,3,is0+1) = smq0(ix,iy,3,is0+1)-
     &         rg0*t2/roxa(ix,iy,is0+1)
              smq0(ix,iy,0,is0) = smq0(ix,iy,0,is0)+
     &         rg0*t2*ua(ix,iy,is0)
              smq0(ix,iy,3,is0) = smq0(ix,iy,3,is0)-
     &         (1.0_R8+rg0)*t2/roxa(ix,iy,is0)
            else if (styl0.eq.3) then
*       ..compute parallel momentum source
              smq0(ix,iy,2,is) = smq0(ix,iy,2,is)+
     &         rg0*t1*ua(ix,iy,is)/roxa(ix,iy,is)
              smq0(ix,iy,1,is) = smq0(ix,iy,1,is)-
     &         (1.0_R8+rg0)*t1
              smq0(ix,iy,2,is-1) = smq0(ix,iy,2,is-1)+
     &         t1*(ua(ix,iy,is)+rg0*ua(ix,iy,is-1))/roxa(ix,iy,is-1)
              smq0(ix,iy,1,is-1) = smq0(ix,iy,1,is-1)-
     &         rg0*t1
              smq0(ix,iy,2,is0+1) = smq0(ix,iy,2,is0+1)+
     &         t2*(ua(ix,iy,is0)+rg0*ua(ix,iy,is0+1))/roxa(ix,iy,is0+1)
              smq0(ix,iy,1,is0+1) = smq0(ix,iy,1,is0+1)-
     &         rg0*t2
              smq0(ix,iy,2,is0) = smq0(ix,iy,2,is0)+
     &         rg0*t2*ua(ix,iy,is0)/roxa(ix,iy,is0)
              smq0(ix,iy,1,is0) = smq0(ix,iy,1,is0)-
     &         (1.0_R8+rg0)*t2
            endif
            rcxmo(ix,iy,is)=rcxmo(ix,iy,is)+t1*ua(ix,iy,is)
            rcxmo(ix,iy,is0)=rcxmo(ix,iy,is0)-t2*ua(ix,iy,is)
            rcxmoreg(region(ix,iy,0),is)=
     1       rcxmoreg(region(ix,iy,0),is)+t1*ua(ix,iy,is)
            rcxmoreg(region(ix,iy,0),is0)=
     1       rcxmoreg(region(ix,iy,0),is0)-t2*ua(ix,iy,is0)
          endif
         enddo
        endif
       enddo
      enddo
#ifndef NO_OPENMP_B2STCX
!$OMP END PARALLEL DO
#ifdef WITH_FTIMINGS
      call timer%stop("parallel")
#endif
#endif
*   ..zero out contributions in dead cells and core boundary
      do is = 0, ns-1
        call b2xzdd (nx, ny, 1, sna0(-1,-1,0,is))
        call b2xzdd (nx, ny, 3, smq0(-1,-1,0,is))
      enddo
      call b2xzdd (nx, ny, 3, shi0)

#ifdef ASTRA
      do is = 0, ns-1                     !srv 28.07.11 {
        if(is.eq.ismain) then
          wrk0=sna0(:,:,0,is)+sna0(:,:,1,is)*na(:,:,is)
          call b2_to_astra_src(nx, ny, vol, wrk0, svol, se)
        endif
      enddo
      wrk0=shi0(:,:,0)+
     &     shi0(:,:,1)*ti+shi0(:,:,2)*ni(:,:,0)+
     &     shi0(:,:,3)*ni(:,:,0)*ti
      call b2_to_astra_src(nx, ny, vol, wrk0, svol, sih)                !srv 11.06.13
#endif

      write (chcx,'(i3.3)') is0
!      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1
!     &              .or.  iout_b2npc7.eq.1 .or. iout_b2npc11.eq.1 )    !srv 02.04.18
!     &        .or. iout_b2wdat.eq.4 ) then
!        do is = 0, ns-1
!          write (chns,'(i3.3)') is
!!         do k=0,1
!!          write (chk,'(i1)') k
!!          call my_out(70,nx,ny,sna0(-1,-1,k,is),
!!     &     'b2stcx_sna'//chk//'_'//chcx//'_'//chns)
!!         enddo
!          wrk0=sna0(:,:,0,is)+sna0(:,:,1,is)*na(:,:,is)
!!          call my_out(70,nx,ny,wrk0,'b2stcx_sna_'//chcx//'_'//chns)
!!          call my_out(70,nx,ny,wrk0,'b2stcx_sna_'//chns)
!        enddo
!      endif
!      if (iout.eq.1 .or. iout_b2npmo.eq.1 .or. iout_b2wdat .eq.4) then
!        do is = 0, ns-1
!          write (chns,'(i3.3)') is
!!         do k=0,1
!!          write (chk,'(i1)') k
!!          call my_out(70,nx,ny,smq0(-1,-1,k,is),
!!     &     'b2stcx_smq'//chk//'_'//chcx//'_'//chns)
!!         enddo
!          wrk0=smq0(:,:,0,is)+                                          !srv 28.02.17
!     &         smq0(:,:,1,is)*ua(:,:,is)+
!     &         smq0(:,:,2,is)*na(:,:,is)*mp*am(is)+
!     &         smq0(:,:,3,is)*ua(:,:,is)*na(:,:,is)*mp*am(is)
!!	  call my_out(70,nx,ny,wrk0,'b2stcx_smq'//chns)
!!          call my_out(70,nx,ny,wrk0,'b2stcx_smq_'//chcx//'_'//chns)
!        enddo
!      endif
!      if (iout.eq.1 .or. (iout_b2nph9.eq.1 .or. iout_b2npht.eq.1)
!     & .or. iout_b2wdat .eq.4) then
!!        do k=0,1
!!         write (chk,'(i1)') k
!!         call my_out(70,nx,ny,shi0(-1,-1,k),'b2stcx_shi'//chk//'_'//chcx)
!!        enddo
!        wrk0=shi0(:,:,0)+shi0(:,:,1)*ti+
!     &       shi0(:,:,2)*ni(:,:,0)+shi0(:,:,3)*ti*ni(:,:,0)
!!        call my_out(70,nx,ny,wrk0,'b2stcx_shi_'//chcx)
!!        call my_out(70,nx,ny,wrk0,'b2stcx_shi')
!      endif
      if (iout.gt.0) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          do k=0,1
            write (chk,'(i1)') k
            call my_out(70,nx,ny,smq0(-1,-1,k,is),
     &        'b2stcx_smq'//chk//'_'//chcx//'_'//chns)
           enddo
        enddo
      endif

cdjm Jan2017 Keep linearised sources for balance
      if (balance_netcdf.ne.0) then
        b2stcx_sna0to1 = sna0
        b2stcx_smq0to3 = smq0
        b2stcx_shi0to3 = shi0
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stcx

      end subroutine b2stcx
