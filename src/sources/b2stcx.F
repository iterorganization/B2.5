*-----------------------------------------------------------------------
*.specification

      subroutine b2stcx (nCv, nFc, ns, is0, ismain, BoRiS, switch, geo,   !srv 09.01.01
     &  mpg, na, ua, ti, tn, ni, rcx, sna0, smq0, shi0, shn0)        !srv 07.07.09
      use b2mod_types
      use b2mod_constants
!      use b2mod_tallies
!      use b2mod_sources
      use b2mod_b2cmpa
#ifdef ASTRA
      use b2mod_b2_to_astra
#endif
!      use b2mod_balance !djm Jan2017
!     & , only : b2stcx_sna0to1, b2stcx_smq0to3, b2stcx_shi0to3,
!     &          balance_netcdf
      use b2mod_switches
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, ns, is0, ismain
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  BoRiS, na(nCv,0:ns-1), ua(nCv,0:ns-1), ti(nCv), tn(nCv),
     *   ni(nCv,0:1), rcx(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(nCv,0:1,0:ns-1), smq0(nCv,0:3,0:ns-1), shi0(nCv,0:3),
     *  shn0(nCv,0:3)
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
      integer ncall, iCv, is, js, k
      character chns*3, chcx*3, chk*1                                   !srv 11.09.09
      real (kind=R8) ::
     *  rg0, t0, t1, t2, du0, du1, tkin, t1i, t1n
      real (kind=R8) :: wrk0(nCv)
      save ncall
*   ..procedures
      intrinsic abs
      logical lnext
      real (kind=R8) ::
     *  samax, roxa
      external subini, subend, xertst, samax, sfill
      external b2xvsg
      lnext(is,js) = zamax(is).lt.zamin(js).and.zn(is).eq.zn(js).and.
     &  abs(am(is)-am(js)).lt.1.0e-3_R8
      roxa(iCv,is) = am(is)*mp*na(iCv,is)
*   ..initialisation
      data ncall/0/


*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stcx')
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test is0
      call xertst (0.le.is0.and.is0.lt.ns-1, 'faulty argument is0')
      call xertst (lnext(is0,is0+1),
     &  'faulty argument is0: not related to species is0+1')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test sign of vol, na, ti, ni
       call b2xvsg (nCv, geo%cvVol, 1, 'vol', '.gt.')
       call b2xvsg (nCv*ns, na, 1, 'na', '.gt.')
       call b2xvsg (nCv, ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv, tn, 1, 'tn', '.gt.')
       call b2xvsg (nCv*2, ni, 1, 'ni', '.gt.')
*    ..test rate coefficients
       do is = 1, ns-1
        if (.not.lnext(is-1,is)) then
         call xertst (samax(nCv,rcx(1,is),1).eq.0.0_R8,
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
      call sfill (nCv*2*ns, 0.0_R8, sna0, 1)
      call sfill (nCv*4*ns, 0.0_R8, smq0, 1)
      call sfill (nCv*4, 0.0_R8, shi0, 1)
      call sfill (nCv*4, 0.0_R8, shn0, 1)
*   ..compute sources due to charge exchange
      do is = 1, ns-1
       if (lnext(is-1,is)) then
*     ..particle source
        if (is-1.eq.is0) then
*        (This is a special case: charge exchange between H0 and H1,
*        for which there is no particle source.)
        else
         do iCv = 1, mpg%nCi
           if(mpg%cvReg(iCv).ne.0) then
*        ..compute rate for is->is-1 and is0->is0+1
            t0 = rcx(iCv,is)*geo%cvVol(iCv)
*        ..compute particle source
            sna0(iCv,0,is-1) = sna0(iCv,0,is-1)+
     &       t0*na(iCv,is0)*na(iCv,is)
            sna0(iCv,0,is0+1) = sna0(iCv,0,is0+1)+
     &       t0*na(iCv,is0)*na(iCv,is)
            sna0(iCv,1,is) = sna0(iCv,1,is)-
     &       t0*na(iCv,is0)
            sna0(iCv,1,is0) = sna0(iCv,1,is0)-
     &       t0*na(iCv,is)
!WG_TODO            rcxna(iCv,is)=rcxna(iCv,is)+
!WG_TODO     1       t0*na(iCv,is0)*na(iCv,is)
!WG_TODO            rcxna(iCv,is0)=rcxna(iCv,is0)-
!WG_TODO     1       t0*na(iCv,is0)*na(iCv,is)
!WG_TODO            rcxnareg(mpg%cvReg(iCv),is)=
!WG_TODO     1       rcxnareg(mpg%cvReg(iCv),is)+t0*na(iCv,is0)*na(iCv,is)
!WG_TODO            rcxnareg(mpg%cvReg(iCv),is0)=
!WG_TODO     1       rcxnareg(mpg%cvReg(iCv),is0)-t0*na(iCv,is0)*na(iCv,is)
           endif
         enddo
        endif
*     ..heat source
        do iCv = 1, mpg%nCi
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute rate for is->is-1 and is0->is0+1
           t0 = rcx(iCv,is)*geo%cvVol(iCv)
*       ..compute atom heat source
           du0 = ua(iCv,is)-ua(iCv,is-1)
           du1 = ua(iCv,is0)-ua(iCv,is0+1)
          if (switch%tn_style.eq.0) then
            t1 = t0*na(iCv,is0)*na(iCv,is)*(mp/2.0_R8)*
     &        (((am(is)+am(is-1))/2.0_R8)*(du0**2)+
     &        ((am(is0)+am(is0+1))/2.0_R8)*(du1**2))*
     &        (1.0_R8-BoRiS)
             shi0(iCv,0) = shi0(iCv,0)+t1
          elseif (switch%tn_style.eq.1) then
            tkin = (mp/4.0_R8)*
     &        (((am(is)+am(is-1))/2.0_R8)*(du0**2)+
     &        ((am(is0)+am(is0+1))/2.0_R8)*(du1**2))
            t1i = t0*na(iCv,is0)*na(iCv,is)*tkin
            shi0(iCv,0) = shi0(iCv,0)+t1i
          else
            tkin = (mp/4.0_R8)*
     &        (((am(is)+am(is-1))/2.0_R8)*(du0**2)+
     &        ((am(is0)+am(is0+1))/2.0_R8)*(du1**2))
            t1i = t0*na(iCv,is0)*na(iCv,is)*
     &            (1.5_R8*(tn(iCv)-ti(iCv)) + tkin)
            t1n = t0*na(iCv,is0)*na(iCv,is)*
     &            (1.5_R8*(ti(iCv)-tn(iCv)) + tkin)
            shi0(iCv,0) = shi0(iCv,0)+t1i
            shn0(iCv,0) = shn0(iCv,0)+t1n
          endif

!WG_TODO           rcxhi(iCv,is)=rcxhi(iCv,is)+t1
!WG_TODO           rcxhireg(mpg%cvReg(iCv),is)=rcxhireg(mpg%cvReg(iCv),is)+t1


          endif
        enddo
*     ..parallel momentum source
        rg0 = switch%b2stcx_rg0
!!   (rg0 for numerical stabilisation; needs experiments.) (default = 1.0)
        do iCv = 1,mpg%nCi
          if(mpg%cvReg(iCv).ne.0) then
*       ..compute rate for is->is-1 and is0->is0+1
            t0 = rcx(iCv,is)*geo%cvVol(iCv)*na(iCv,is0)*na(iCv,is)
     &           *geo%cvHz(iCv)                                          !srv 09.01.01
            t1 = ((am(is)+am(is-1))/2.0_R8)*mp*t0
            t2 = ((am(is0)+am(is0+1))/2.0_R8)*mp*t0
            if (switch%b2stcx_styl0.eq.0) then
*       ..compute parallel momentum source
              smq0(iCv,0,is) = smq0(iCv,0,is)+
     &         rg0*t1*ua(iCv,is)
              smq0(iCv,1,is) = smq0(iCv,1,is)-
     &         (1.0_R8+rg0)*t1
              smq0(iCv,0,is-1) = smq0(iCv,0,is-1)+
     &         t1*(ua(iCv,is)+rg0*ua(iCv,is-1))
              smq0(iCv,1,is-1) = smq0(iCv,1,is-1)-
     &         rg0*t1
              smq0(iCv,0,is0+1) = smq0(iCv,0,is0+1)+
     &         t2*(ua(iCv,is0)+rg0*ua(iCv,is0+1))
              smq0(iCv,1,is0+1) = smq0(iCv,1,is0+1)-
     &         rg0*t2
              smq0(iCv,0,is0) = smq0(iCv,0,is0)+
     &         rg0*t2*ua(iCv,is0)
              smq0(iCv,1,is0) = smq0(iCv,1,is0)-
     &         (1.0_R8+rg0)*t2
            else if (switch%b2stcx_styl0.eq.1) then
*       ..compute parallel momentum source
              smq0(iCv,2,is) = smq0(iCv,2,is)+
     &         rg0*t1*ua(iCv,is)/roxa(iCv,is)
              smq0(iCv,3,is) = smq0(iCv,3,is)-
     &         (1.0_R8+rg0)*t1/roxa(iCv,is)
              smq0(iCv,2,is-1) = smq0(iCv,2,is-1)+
     &         t1*(ua(iCv,is)+rg0*ua(iCv,is-1))/roxa(iCv,is-1)
              smq0(iCv,3,is-1) = smq0(iCv,3,is-1)-
     &         rg0*t1/roxa(iCv,is-1)
              smq0(iCv,2,is0+1) = smq0(iCv,2,is0+1)+
     &         t2*(ua(iCv,is0)+rg0*ua(iCv,is0+1))/roxa(iCv,is0+1)
              smq0(iCv,3,is0+1) = smq0(iCv,3,is0+1)-
     &         rg0*t2/roxa(iCv,is0+1)
              smq0(iCv,2,is0) = smq0(iCv,2,is0)+
     &         rg0*t2*ua(iCv,is0)/roxa(iCv,is0)
              smq0(iCv,3,is0) = smq0(iCv,3,is0)-
     &         (1.0_R8+rg0)*t2/roxa(iCv,is0)
            else if (switch%b2stcx_styl0.eq.2) then
*       ..compute parallel momentum source
              smq0(iCv,0,is) = smq0(iCv,0,is)+
     &         rg0*t1*ua(iCv,is)
              smq0(iCv,3,is) = smq0(iCv,3,is)-
     &         (1.0_R8+rg0)*t1/roxa(iCv,is)
              smq0(iCv,0,is-1) = smq0(iCv,0,is-1)+
     &         t1*(ua(iCv,is)+rg0*ua(iCv,is-1))
              smq0(iCv,3,is-1) = smq0(iCv,3,is-1)-
     &         rg0*t1/roxa(iCv,is-1)
              smq0(iCv,0,is0+1) = smq0(iCv,0,is0+1)+
     &         t2*(ua(iCv,is0)+rg0*ua(iCv,is0+1))
              smq0(iCv,3,is0+1) = smq0(iCv,3,is0+1)-
     &         rg0*t2/roxa(iCv,is0+1)
              smq0(iCv,0,is0) = smq0(iCv,0,is0)+
     &         rg0*t2*ua(iCv,is0)
              smq0(iCv,3,is0) = smq0(iCv,3,is0)-
     &         (1.0_R8+rg0)*t2/roxa(iCv,is0)
            else if (switch%b2stcx_styl0.eq.3) then
*       ..compute parallel momentum source
              smq0(iCv,2,is) = smq0(iCv,2,is)+
     &         rg0*t1*ua(iCv,is)/roxa(iCv,is)
              smq0(iCv,1,is) = smq0(iCv,1,is)-
     &         (1.0_R8+rg0)*t1
              smq0(iCv,2,is-1) = smq0(iCv,2,is-1)+
     &         t1*(ua(iCv,is)+rg0*ua(iCv,is-1))/roxa(iCv,is-1)
              smq0(iCv,1,is-1) = smq0(iCv,1,is-1)-
     &         rg0*t1
              smq0(iCv,2,is0+1) = smq0(iCv,2,is0+1)+
     &         t2*(ua(iCv,is0)+rg0*ua(iCv,is0+1))/roxa(iCv,is0+1)
              smq0(iCv,1,is0+1) = smq0(iCv,1,is0+1)-
     &         rg0*t2
              smq0(iCv,2,is0) = smq0(iCv,2,is0)+
     &         rg0*t2*ua(iCv,is0)/roxa(iCv,is0)
              smq0(iCv,1,is0) = smq0(iCv,1,is0)-
     &         (1.0_R8+rg0)*t2
            endif
!WG_TODO            rcxmo(iCv,is)=rcxmo(iCv,is)+t1*ua(iCv,is)
!WG_TODO            rcxmo(iCv,is0)=rcxmo(iCv,is0)-t2*ua(iCv,is)
!WG_TODO            rcxmoreg(mpg%cvReg(iCv),is)=
!WG_TODO     1       rcxmoreg(mpg%cvReg(iCv),is)+t1*ua(iCv,is)
!WG_TODO            rcxmoreg(mpg%cvReg(iCv),is0)=
!WG_TODO     1       rcxmoreg(mpg%cvReg(iCv),is0)-t2*ua(iCv,is0)
          endif
        enddo
       endif
      enddo

#ifdef ASTRA
      do is = 0, ns-1                     !srv 28.07.11 {
        if(is.eq.ismain) then
          wrk0=sna0(:,0,is)+sna0(:,1,is)*na(:,is)
          call b2_to_astra_src(nx, ny, vol, wrk0, svol, se)
        endif
      enddo
      wrk0=shi0(:,:,0)+
     &     shi0(:,:,1)*ti+shi0(:,:,2)*ni(:,:,0)+
     &     shi0(:,:,3)*ni(:,:,0)*ti
      call b2_to_astra_src(nx, ny, vol, wrk0, svol, sih)                !srv 11.06.13
#endif

      write (chcx,'(i3.3)') is0
!      if (iout.eq.1 .or. (iout_b2npc9.eq.1 .or. iout_b2npco.eq.1)
!     &              .or. iout_b2wdat.eq.4) then
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
      if (switch%b2stcx_iout.gt.0) then
        do is = 0, ns-1
          write (chns,'(i3.3)') is
          do k=0,1
            write (chk,'(i1)') k
            call my_out_us(70,nCv,0,smq0(1,k,is),
     &        'b2stcx_smq'//chk//'_'//chcx//'_'//chns)
           enddo
        enddo
      endif


cdjm Jan2017 Keep linearised sources for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2stcx_sna0to1 = sna0
!WG_TODO        b2stcx_smq0to3 = smq0
!WG_TODO        b2stcx_shi0to3 = shi0
!WG_TODO      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2stcx

      end subroutine b2stcx
