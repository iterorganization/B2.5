*-----------------------------------------------------------------------
*.specification

      subroutine b2tqna (nCv, nFc, nVx, ns, nscx, nscxmax, iscx, ismain,
     &   switch, geo, mpg, pl, dv, rt, st_ext, dna0, dpa0, vla0, vma0,
     &   vsa0, hci0, hcn0, hcib, hce0, sig0, alf0, hvi0, hve0, dkt0,
     &   dzt0, dna_ExB, hce_ExB, hci_ExB)
      use b2mod_types
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
!      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_b2cmpb
!      use b2mod_ranges
!      use b2mod_neoclassical
!      use b2mod_geo_corner
!     &                    , only : lengthy, nmdpl                 !srv 26.11.04
      use b2mod_user_namelist
     , , only : omp
      use b2mod_switches
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_par_opt
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2tqna, ncall_transp_keps
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, nscxmax, iscx(0:nscxmax-1),
     *  ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dna0(nCv,0:ns-1), dpa0(nCv,0:ns-1),
     &  vla0(nCv,0:1,0:ns-1), vsa0(nCv,0:ns-1),
     &  hci0(nCv), hcn0(nCv), hcib(nCv,0:ns-1),
     &  vma0(nCv,0:1,0:ns-1)                              !srv 15.12.05
      real (kind=R8) ::
     *  hce0(nCv), sig0(nCv), alf0(nCv),
     &  hvi0(nCv,0:1), hve0(nCv,0:1), dkt0(nCv), dzt0(nCv),
     &  dna_ExB(nCv), hce_ExB(nCv), hci_ExB(nCv)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(nCv)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TQNA computes the anomalous contributions to the transport
*     coefficients.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine computes the anomalous contributions to the following
*     transport coefficients.
*     - dna: particle transport coefficients for grad(na).
*     - dpa: particle transport coefficients for grad(pa).
*     - vla: anomalous velocity.
*     - vma: anomalous velocity in parallel momentum balance equation. !srv 15.12.05
*     - vsa: viscosity.
*     - hci: atom thermal conductivity.
*     - hce: electron thermal conductivity.
*     - sig: electrical conductivity.
*     - alf: thermo-electric coefficient.
*     - hvi0: all atom thermal strange velocity.
*     - hve0: electron thermal strange velocity.
*     These coefficients occur in the equations as indicated below.
*     Particle transport.
*       fnax(is) = -dna0(is)*(d/dx)na(is)-dpa0(is)*(d/dx)pa(is)+
*         vla0*na(is)
*       fnay(is) = -dna0(is)*(d/dy)na(is)-dpa0(is)*(d/dy)pa(is)+
*         vla0*na(is)
*     Parallel momentum transport.
*       momx(is) = na(is)*up(is)*vma0 - vsa0(is)*(d/dx)up(is)     !srv 15.12.05
*       momy(is) = na(is)*up(is)*vma0 - vsa0(is)*(d/dy)up(is)     !srv 15.12.05
*     Atom heat flux.
*       fhix = -hci0*(d/dx)ti+hvi0*na*ti
*       fhiy = -hci0*(d/dy)ti+hvi0*na*ti
*     Electron heat flux.
*       fhex = alf0*te*ehx-hce0*(d/dx)te+hve0*ne*te
*       fhey = alf0*te*ehy-hce0*(d/dy)te+hve0*ne*te
*     Electric current.
*       fchx = sig0*ehx-alf0*(d/dx)te
*       fchy = sig0*ehy-alf0*(d/dy)te
*     Here, ehx and ehy are the components of the modified electric
*     field: ehx = -(d/dx)po+(1/(qe*ne))*(d/dx)(ne*te),
*     ehy = -(d/dy)po+(1/(ne*te))*(d/dy)(ne*te). See Balescu /1/.
*     The particle and parallel momentum transport is a property of
*     each atomic species, the other transport laws involve summation
*     over all species.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iCv, iFc, is, ic, k,
     2  jsep, jxi, jxa
      character chns*3                                            !srv 10.11.09
      real (kind=R8) ::
     *  df0, dfvsb, dfhcb, dfhce, dfsig, dfalf,
     1  dcx, dion, vcx, lcx, lion, din, dnn,
     1  poloidal_scale_factor
      real (kind=R8) ::
     *  t0, t1,
     *  min_df0_tmp, max_df0_tmp, tav
      integer elm_count, iv
      save jxi, jxa, jsep
      logical catch_up
*   ..procedures
      intrinsic sqrt, max, exp
      external subini, subend, xertst, sfill, ipgeti, ipgetr, expu
*   ..initialisation
*   ..statement functions
      real (kind=R8) ::
     *  expu, addup
      addup(t0,t1)=t0+t1

*-----------------------------------------------------------------------
*.computation
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqna')
*   ..set internal parameters on first call
      if (ncall_b2tqna.eq.0) then
!WG_TODO        call get_jsep(nx,ny,jxi,jxa,jsep)
!WG_TODO        switch%b2tqna_ixref=jxa
!WG_TODO        call xertst (-1.le.switch%b2tqna_ixref.and.switch%b2tqna_ixref.le.nx,
!WG_TODO     1   'faulty argument b2tqna_ixref')
        if(switch%b2tqna_transport_namelist.eq.1) then
          call write_b2mod_transport_namelist
          if(transport_time_mod.gt.0.0_R8) then
            elm_count=tim/transport_time_mod
          else
            elm_count=0
          endif
          catch_up=transport_time_mod.gt.0.0_R8
          do while(catch_up)
            if(mod(tim,transport_time_mod).gt.transport_time_switch)
     1       then
              write(*,*) 'b2tqna: catching up ...'
              call read_b2mod_transport_namelist
              call write_b2mod_transport_namelist
              catch_up=transport_time_mod.gt.0.0_R8
            else
              catch_up=.false.
            endif
          enddo
        endif
        if(switch%b2tqna_ballooning.ne.0.0_R8.and.
     &     switch%b2tqna_bb_ref.eq.0.0_R8) then
          switch%b2tqna_bb_ref=sum(geo%cvBb(:,3))/(nCv)
          write(*,*) 'b2tqna: bb_ref set to ',switch%b2tqna_bb_ref
        endif
      else
        if(switch%b2tqna_transport_namelist.eq.1.and.
     1   transport_time_mod.gt.0.0_R8) then
          write(*,*) 'ELM: ',tim,mod(tim,transport_time_mod),
     1     transport_time_switch,int(tim/transport_time_mod),elm_count
          if(mod(tim,transport_time_mod).gt.transport_time_switch.or.
     1     int(tim/transport_time_mod).gt.elm_count) then
            elm_count=tim/transport_time_mod
            call read_b2mod_transport_namelist
            call write_b2mod_transport_namelist
          endif
        endif
      endif
*   ..test nCv, nFc, ns
      call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall_b2tqna.lt.3) then
*    ..test bb
       call b2xvsg (nCv, geo%fcBb(1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (nCv*ns, pl%na, 1, 'na', '.gt.')
       call b2xvsg (nCv, pl%te, 1, 'te', '.gt.')
       call b2xvsg (nCv, pl%ti, 1, 'ti', '.gt.')
       call b2xvsg (nCv*2, dv%ni, 1, 'ni', '.gt.')
       call b2xvsg (nCv, dv%ne, 1, 'ne', '.gt.')
       call b2xvsg (nCv, dv%ne2, 1, 'ne2', '.gt.')
      endif

*   ..compute auxiliary wrk0
      wrk0 = sqrt(pl%ti/mp)
      if (flag_optim) then
        do iv = 1, nnvar !loop over all design variables
          select case (xtype(iv))
          case (1) ! dna
            parm_dna(1) = par_opt_phys(iv) ! FIXME to improve for multispecies
          case (2) ! hce
            parm_hce = par_opt_phys(iv)
          case (3) ! hci
            parm_hci(1) = par_opt_phys(iv) ! FIXME to improve for multispecies
          case default
            write(*,*) xtype(iv)
            call xerrab ('xtype out of bounds')
          end select
        end do
      endif
      if(switch%b2tqna_transport_namelist.eq.1) then
        do is=0,ns-1
          if(.not.is_neutral(is)) then
* dna
            call transform_transport(flag_dna,parm_dna(is),cfdna(0,is))
* dpa
            call transform_transport(flag_dpa,parm_dpa(is),cfdpa(0,is))
* vla
            call transform_transport(flag_vla,parm_vla(is),cfvla(0,is))
* vsa
            call transform_transport(flag_vsa,parm_vsa(is),cfvsa(0,is))
* hci
            call transform_transport(flag_hci,parm_hci(is),cfhci(0,is))
          endif
        enddo
* hce
        call transform_transport(flag_hce,parm_hce,cfhce)
* sig
        call transform_transport(flag_sig,parm_sig,cfsig)
* alf
        call transform_transport(flag_alf,parm_alf,cfalf)
      endif
#ifdef PAR_OPT
      cfhce(0)   = par_opt_phys(1)                               !srv 05.04.09
      cfhci(0,1) = par_opt_phys(2)                               !srv 05.04.09
      cfdna(0,1) = par_opt_phys(3)
      cfvsa(0,1) = cfdna(0,1)                                    !srv 22.07.10
#endif
* ..compute anomalous transport coefficients
*   ..compute dna0, dpa0, vla0, vma0, vsa0, hci0
*     (initialise hci0 to 0)
      call sfill (nCv, 0.0_R8, hci0, 1)
      call sfill (nCv, 0.0_R8, hcn0, 1)
      call sfill (nCv*ns, 0.0_R8, dna0, 1)
      call sfill (nCv*ns, 0.0_R8, dpa0, 1)
      call sfill (nCv*ns*2, 0.0_R8, vla0, 1)
      call sfill (nCv*ns, 0.0_R8, vsa0, 1)
      call sfill (nCv*ns, 0.0_R8, hcib, 1)
      call sfill (nCv*ns*2, 0.0_R8, vma0, 1)
      call sfill (nCv, 0.0_R8, hce0, 1)
      call sfill (nCv, 0.0_R8, sig0, 1)
      call sfill (nCv, 0.0_R8, alf0, 1)
      call sfill (nCv*2, 0.0_R8, hvi0, 1)
      call sfill (nCv*2, 0.0_R8, hve0, 1)
      call sfill (nCv, 0.0_R8, dkt0, 1)
      call sfill (nCv, 0.0_R8, dzt0, 1)
      call sfill (nCv, 0.0_R8, dna_ExB, 1)
      call sfill (nCv, 0.0_R8, hce_ExB, 1)
      call sfill (nCv, 0.0_R8, hci_ExB, 1)
      do is = 0, ns-1
        max_df0_tmp=0.0_R8
        min_df0_tmp=1.0e30_R8
        if (is_neutral(is).and.zn(is).eq.1.0_R8) then
          k=0
          do while (iscx(k).ne.is .and. k.lt.nscx)
            k = k+1
          enddo
          call xertst (k.lt.nscx, 'CX species index not found!')
        endif
        do iCv = 1, nCv
*      ..compute auxiliary quantity df0
            if (cfdf0(2,is).eq.0.0_R8) then
              df0 = cfdf0(0,is)+cfdf0(1,is)/
     &         ((dv%ni(iCv,0)+st_ext%ni(iCv,0))/1.0e20_R8)
            else
              df0 = cfdf0(0,is)+cfdf0(1,is)/
     &         ((dv%ni(iCv,0)+st_ext%ni(iCv,0))/1.0e20_R8)+
     &         (cfdf0(2,is)/sqrt(am(is)))*wrk0(iCv)/
     &         (cfdf0(3,is)*(dv%ni(iCv,0)+st_ext%ni(iCv,0))+
     &          cfdf0(4,is)*dv%ne(iCv))
            endif
            if(switch%b2tqna_new_df0.eq.1) then
              if(is_neutral(is)) then
                dcx=0.0_R8
                vcx=0.0_R8
                tav = (pl%ti(iCv)+pl%tn(iCv))/2.0_R8
                do ic = 0, ns-1
                  dcx = dcx + expu(rt%rlcx(iCv,0,ic,k) +
     &             rt%rlcx(iCv,1,ic,k)*log(tav/(am(is)*ev)))
                  vcx = vcx + pl%na(iCv,ic) * expu(rt%rlcx(iCv,0,ic,k) +
     &             rt%rlcx(iCv,1,ic,k)*log(tav/(am(is)*ev)))
                enddo
                lcx=wrk0(iCv)/sqrt(am(is))/vcx
                dion=expu(
     2           rt%rlsa(iCv,0,is)+
     3           rt%rlsa(iCv,1,is)*log(pl%te(iCv)/ev))
                lion=wrk0(iCv)/sqrt(am(is))/dv%ne(iCv)/dion
                din=(wrk0(iCv)/sqrt(am(is)))**2/
     1           (vcx+dion*dv%ne(iCv))
c               if(lcx.gt.lion) then
c                 din=lion*wrk0(iCv)/sqrt(am(is))
c               else
c                 din=sqrt(lcx*lion)*wrk0(iCv)/sqrt(am(is))
c               endif
c-wrong(Bonnin,Owen)                dnn=1.0_R8/na(iCv,is)/2.0e-19_R8
c$$$                dnn=(wrk0(iCv)/sqrt(am(is)))/
c$$$     1           (na(iCv,is)/switch%neutral_sources_rescale)/2.0e-19_R8
cxpb Corrected for dimensionality 
                dnn=1.0_R8/(pl%na(iCv,is)/
     &            switch%neutral_sources_rescale)/2.0e-19_R8
                if (switch%b2tqna_diagno.gt.1) then
                  write(*,*) 'new_df0',
     1             iCv,is,lcx,lion,din,dnn,addup(din,dnn),df0
                else if (addup(din,dnn)/df0.gt.100.0_R8
     &           .and. switch%b2tqna_diagno.eq.1) then
                  write(*,*) 'new_df0:problem',
     1             iCv,is,lcx,lion,din,dnn,addup(din,dnn),df0
                endif
                df0=addup(din,dnn)
              endif
            endif
            min_df0_tmp=min(min_df0_tmp,df0)
            max_df0_tmp=max(max_df0_tmp,df0)
            df0=max(switch%b2tqna_min_df0,
     &          min(switch%b2tqna_max_df0,df0))
*      ..compute dna0, dpa0, vla0, vma0
            dna0(iCv,is) = cfdna(0,is)+
     &       cfdna(1,is)/(max(1.0e-50_R8,dv%ne(iCv))/1.0e20_R8)+    ! dpc changed na(iCv,is) to ne(iCv)
     &       cfdna(2,is)*df0+
     &       cfdna(3,is)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
            if(cfdna(7,is).ne.0.0_R8)
     &       dna0(iCv,is) = dna0(iCv,is)!WG_TODO*
!WG_TODO     &        cfdna(7,is)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
            dpa0(iCv,is) =
     &       (1.0_R8/(rt%rza(iCv,is)*pl%te(iCv)+pl%ti(iCv)))*
     &       (cfdpa(0,is)+cfdpa(1,is)/
     &       (max(1.0e-50_R8,dv%ne(iCv))/1.0e20_R8)+    ! dpc changed na(iCv,is) to ne(iCv)
     &        cfdpa(2,is)*df0+
     &        cfdpa(3,is)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3)))
            if(cfdpa(7,is).ne.0.0_R8)
     &       dpa0(iCv,is) = dpa0(iCv,is)!WG_TODO*
!WG_TODO     &        cfdpa(7,is)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
            vla0(iCv,0,is) = 0.0_R8
            vla0(iCv,1,is) = cfvla(0,is)+
     &       cfvla(1,is)/(max(1.0e-50_R8,dv%ne(iCv))/1.0e20_R8)+    ! dpc changed na(iCv,is) to ne(iCv)
     &       cfvla(2,is)*df0+
     &       cfvla(3,is)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
            if(cfvla(7,is).ne.0.0_R8)
     &       vla0(iCv,1,is) = vla0(iCv,1,is)!WG_TODO*
!WG_TODO     &        cfvla(7,is)*hy1(iCv)/hy1(switch%b2tqna_ixref,iy)
            if(.not.is_neutral(is)) then                             !srv 15.12.05 {
              vma0(iCv,0,is) = 0.0_R8
              vma0(iCv,1,is) = 0.0_R8
            else
              vma0(iCv,0,is) = 0.0_R8
!             vma0(iCv,1,is) = cfvma(0,is)
              vma0(iCv,1,is) = switch%b2tqna_cfvma                   !srv 15.12.05 }
            endif
*      ..compute dfvsb and dfhcb
*        (dfvsb and dfhcb specify effective diffusivities, of
*         dimension m**2/s, for the anomalous transport coefficients
*         vsa0 and hci0)
            dfvsb = cfvsa(0,is)+
     &       cfvsa(1,is)/(max(1.0e-50_R8,dv%ne(iCv))/1.0e20_R8)+    ! dpc changed na(iCv,is) to ne(iCv)
     &       cfvsa(2,is)*df0+
     &       cfvsa(3,is)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
            if(cfvsa(7,is).ne.0.0_R8)
     &       dfvsb = dfvsb!WG_TODO*cfvsa(7,is)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
            dfhcb = cfhci(0,is)+
     &       cfhci(1,is)/(max(1.0e-50_R8,dv%ne(iCv))/1.0e20_R8)+    ! dpc changed na(iCv,is) to ne(iCv)
     &       cfhci(2,is)*df0+
     &       cfhci(3,is)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
            if(cfhci(7,is).ne.0.0_R8)
     &       dfhcb = dfhcb!WG_TODO*cfhci(7,is)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
*      ..compute vsa0, hci0
            vsa0(iCv,is) = dfvsb*mp*am(is)*pl%na(iCv,is)
            hcib(iCv,is)=dfhcb*pl%na(iCv,is)
            if (.not.is_neutral(is)) then
              hci0(iCv) = hci0(iCv)+hcib(iCv,is)
            else
              hcn0(iCv) = hcn0(iCv)+hcib(iCv,is)
            endif
*      ..compute hvi0
            hvi0(iCv,0) = 0.0_R8
            hvi0(iCv,1) = 0.0_R8
        enddo
        write(*,'(a,i3,1p,2g14.7)') 'is,min_df0,max_df0 ',
     .   is, min_df0_tmp, max_df0_tmp
      enddo
      if (switch%tn_style.eq.0) then     !! default combined ion-neutral, add neutral contribution to "total" hci0
        hci0 = hci0+hcn0
      endif
*   ..compute hce0, sig0, alf0
      do iCv = 1, nCv
*     ..compute dfhce, dfsig, dfalf
*       (dfhce, dfsig, and dfalf specify effective diffusivity,
*        dimension m**2/s, for the anomalous transport coefficients
*        hce0, sig0, and alf0)
         dfhce = cfhce(0)+cfhce(1)/(dv%ne(iCv)/1.0e20_R8)+
     1    cfhce(3)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
         if(cfhce(7).ne.0.0_R8)
     &    dfhce=dfhce!WG_TODO*cfhce(7)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
         dfsig = cfsig(0)+cfsig(1)/(dv%ne(iCv)/1.0e20_R8)+
     1    cfsig(3)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))!WG_TODO+
!WG_TODO     2    cfsig(4)*exp(-cfsig(5)*(lengthy(iCv)/lengthy(ix,ny))**2) !srv 26.11.04
         if(cfsig(7).ne.0.0_R8)
     &    dfsig=dfsig!WG_TODO*cfsig(7)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
         dfalf = cfalf(0)+cfalf(1)/(dv%ne(iCv)/1.0e20_R8)+
     1    cfalf(3)*6.25e-2_R8*pl%te(iCv)/ev/abs(geo%cvBb(iCv,3))
         if(cfalf(7).ne.0.0_R8)
     &    dfalf=dfalf!WG_TODO*cfalf(7)*(hy1(iCv)/hy1(switch%b2tqna_ixref,iy))**2
*     ..compute hce0, sig0, alf0
         hce0(iCv) = dfhce*dv%ne(iCv)
         if(switch%model_sig.eq.0) then                !srv 02.12.04 {
           sig0(iCv) = dfsig*qe*dv%ne(iCv)
         else
!          sig0(iCv) = dfsig*qe*(cbsna(0,1,1)/abs(cbsna(1,1,1))) ! use constant density as core region
           sig0(iCv) = dfsig*qe*dv%ne(omp(1))  !srv 10.06.09 ! use constant density as core region
         endif                                  !srv 02.12.04 }
         alf0(iCv) = dfalf*dv%ne(iCv)*sqrt(qe/pl%te(iCv))
      enddo

      if(switch%b2tqna_transport_inputfile.eq.1) then
        call transport_input(nCv, nFc, ns, geo, mpg, pl, dv, rt, hcib,
     &   dna0, dpa0, vla0, vsa0, hci0, hce0, sig0, alf0, ncall_b2tqna)
        do iCv = 1, nCv
          hci0(iCv)=0.0_R8
          do is=0,ns-1
            hci0(iCv)=hci0(iCv)+hcib(iCv,is)
          enddo
        enddo
      endif

      if(switch%b2tqna_iout.ne.0) then                               !srv 10.11.09 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,dna0(1,is),'b2tqna_in_dna0'//chns)
          call my_out_us(70,nCv,0,dpa0(1,is),'b2tqna_in_dpa0'//chns)
!          wrk0=vsa0(:,is)/(mp*am(is)*pl%na(:,is))
          wrk0=vsa0(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2tqna_in_vsa0'//chns)
          call my_out_us(70,nCv,0,vla0(1,0,is),'b2tqna_in_vla0x'//chns)
          call my_out_us(70,nCv,0,vla0(1,1,is),'b2tqna_in_vla0y'//chns)
!          wrk0=hcib(:,is)/pl%na(:,is)
          wrk0=hcib(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2tqna_in_hcib'//chns)
        enddo
       call my_out_us(70,nCv,0,hci0,'b2tqna_in_hci0')
!        wrk0=hce0/dv%ne
        wrk0=hce0
        call my_out_us(70,nCv,0,wrk0,'b2tqna_in_hce0')
        call my_out_us(70,nCv,0,sig0,'b2tqna_in_sig0')
        call my_out_us(70,nCv,0,alf0,'b2tqna_in_alf0')
      endif                                                           !srv 10.11.09 }

c lvn1

      if (switch%user_transport.eq.0) then
        !! the default transport model calculated above is used. Nothing to do.
      elseif (switch%user_transport.eq.10) then                              !mb 24.04.2017
        call set_transport_kul(nCv, nFc,
     &   nVx, ns, nscx, iscx, switch%tn_style, geo, 
     &   mpg, pl, dv, rt, st_ext, dna0, dpa0, vla0, vma0, vsa0, hci0,
     &   hcn0,hcib)
      else
        call xerrab('transport model undefined or not yet converted')
      endif



#ifdef WG_TODO
      if(user_transport.eq.0) then
c do nothing
      elseif(user_transport.eq.1) then
        call set_transport(nx, ny, ns, bb, na, te, ti, ni, ne, ne2,
     &  hy, qz, dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0,
     &  wrk0)
      elseif(user_transport.eq.2) then
        call set_transport_swy(nx, ny, ns,
     &   dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0)
      elseif(user_transport.eq.3) then
        call set_transport_eta(nx, ny, ns, bb, na, te, ti, ni, ne, ne2,
     &   hy, dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0)
      elseif(user_transport.eq.4) then
        call set_transport_xpb(nx, ny, ns, bb, na, ne, te, ti, ni,
     &   hx, hy, qz, dna0, vla0)
c mw added transport
      elseif(user_transport.eq.5) then
        call set_transport_mw(nx, ny, ns, ne, jsep, jxa,
     &   dna0, vla0, hci0, hcib, hce0)
c mw end
      elseif(user_transport.eq.6) then                                    !srv 13.06.08
!! FIXME: Need to fix clash with ank blobby transport model
        call set_transport_srv(nx, ny, ns, bb, na, te, ti, ni, ne, ne2,
     &  hy, qz, dna0, dpa0, vla0, vma0, vsa0, hci0, hcib, hce0, sig0,
     &  alf0, wrk0)
!! FIXME: in routines below, remove hvi0, hve0 ???
      elseif(user_transport.eq.7) then
        call set_transport_disruption(nx, ny, ns,
     1   dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0,
     2   na, ne, te, ti, crx, hy1)
      elseif(user_transport.eq.8) then
        if(ncall_b2tqna.eq.0) then
          call alloc_b2mod_neoclassical(nx,ny,nsdecl,mmx)
        endif
        call set_transport_neo(nx, ny, ns, hx, hy, vol, crx, cry,
     &   conn, qz, bb, na, te, ti, ni, ne, dna0, dpa0, vla0, vsa0,
     &   hci0, hcib, hce0, sig0, alf0, hvi0, hve0)
c ank added blobby transport
!! FIXME: Was user_transport.eq.6 in 4.x, 5.0/5.1
      elseif(user_transport.eq.9) then
        call set_transport_blob(nx, ny, ns, jsep, jxa, vla0)
      else
        write(*,*) 'user_transport =',switch%user_transport,' not coded'
        call xerrab ('user_transport value not coded')
      endif
#endif
c sc  implementation of turbulence closure models 29.11.2018-23/9/2020
      if(switch%transport_keps.ne.0) then
        call set_transport_keps( 
     &    nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &    pl, dv, rt, dna0, vsa0, hce0, hci0, hcib, sig0, alf0,
     &    dkt0, dzt0, dna_ExB,
     &    hce_ExB, hci_ExB)
      endif

cdpc
      if(switch%b2tqna_ballooning.ne.0.0_R8) then
        do iCv = 1, nCv
            poloidal_scale_factor=
     1       abs(switch%b2tqna_bb_ref/geo%cvBb(iCv,3))**
     1       switch%b2tqna_ballooning*
     1       switch%b2tqna_ballooning_rescale
            do is=0,ns-1
              if(.not.is_neutral(is)) then
                dna0(iCv,is)=dna0(iCv,is)*poloidal_scale_factor
                dpa0(iCv,is)=dpa0(iCv,is)*poloidal_scale_factor
                vla0(iCv,:,is)=vla0(iCv,:,is)*poloidal_scale_factor
                vsa0(iCv,is)=vsa0(iCv,is)*poloidal_scale_factor
                hci0(iCv)=hci0(iCv)+
     1           hcib(iCv,is)*(poloidal_scale_factor-1.0_R8)
                hcib(iCv,is)=hcib(iCv,is)*poloidal_scale_factor
              endif
            enddo
            hce0(iCv)=hce0(iCv)*poloidal_scale_factor
            if (switch%b2tqna_ballooning_sig.ne.0)  then
              sig0(iCv)=sig0(iCv)*poloidal_scale_factor
              alf0(iCv)=alf0(iCv)*poloidal_scale_factor
            endif
        enddo
      endif
#ifdef WG_TODO
      if(switch%b2tqna_pfr_rescale.ne.1.0_R8) then
        do iCv = 1, nCv
            if(mpg%cvReg(iCv,0).ne.0.and.
     1       mod(mpg%cvReg(iCv,0),4).ne.1.and.
     1       mod(mpg%cvReg(iCv,0),4).ne.2.and.
     2       (((ix.lt.leftcut(1).or.ix.ge.rightcut(1)).and.
     3       iy.lt.topcut(1)).or.(nnreg(0).eq.8.and.
     4       ix.ge.leftcut(2).and.ix.lt.rightcut(2).and.
     5       iy.lt.topcut(2)))) then
              do is=0,ns-1
                if(.not.is_neutral(is)) then
                  dna0(iCv,is)=dna0(iCv,is)*switch%b2tqna_pfr_rescale
                  dpa0(iCv,is)=dpa0(iCv,is)*switch%b2tqna_pfr_rescale
                  vla0(iCv,:,is)=vla0(iCv,:,is)*
     1             switch%b2tqna_pfr_rescale
                  vsa0(iCv,is)=vsa0(iCv,is)*switch%b2tqna_pfr_rescale
                  hci0(iCv)=hci0(iCv)+
     1             hcib(iCv,is)*(switch%b2tqna_pfr_rescale-1.0_R8)
                  hcib(iCv,is)=hcib(iCv,is)*switch%b2tqna_pfr_rescale
                endif
              enddo
              hce0(iCv)=hce0(iCv)*switch%b2tqna_pfr_rescale
              sig0(iCv)=sig0(iCv)*switch%b2tqna_pfr_rescale
              alf0(iCv)=alf0(iCv)*switch%b2tqna_pfr_rescale
            endif
        enddo
      endif
      if(switch%b2tqna_divsol_rescale.ne.1.0_R8) then
        do iy=jsep+1,ny
          do ix=-1,nx
            if(mpg%cvReg(iCv,0).ne.0.and.
     1       mod(mpg%cvReg(iCv,0),4).ne.1.and.
     1       mod(mpg%cvReg(iCv,0),4).ne.2.and.
     2       (((ix.lt.leftcut(1).or.ix.ge.rightcut(1)).and.
     3       iy.ge.topcut(1)).or.(nnreg(0).eq.8.and.
     4       ix.ge.leftcut(2).and.ix.lt.rightcut(2).and.
     5       iy.ge.topcut(2)))) then
              do is=0,ns-1
                if(.not.is_neutral(is)) then
                  dna0(iCv,is)=dna0(iCv,is)*switch%b2tqna_divsol_rescale
                  dpa0(iCv,is)=dpa0(iCv,is)*switch%b2tqna_divsol_rescale
                  vla0(iCv,:,is)=vla0(iCv,:,is)*
     1             switch%b2tqna_divsol_rescale
                  vsa0(iCv,is)=vsa0(iCv,is)*switch%b2tqna_divsol_rescale
                  hci0(iCv)=hci0(iCv)+
     1             hcib(iCv,is)*(switch%b2tqna_divsol_rescale-1.0_R8)
                  hcib(iCv,is)=hcib(iCv,is)*switch%b2tqna_divsol_rescale
                endif
              enddo
              hce0(iCv)=hce0(iCv)*switch%b2tqna_divsol_rescale
              sig0(iCv)=sig0(iCv)*switch%b2tqna_divsol_rescale
              alf0(iCv)=alf0(iCv)*switch%b2tqna_divsol_rescale
            endif
          enddo
        enddo
      endif
#endif
      if(switch%iout_b2wdat.eq.4) then                                 !srv 10.11.09 {
        do is = 0,ns-1
          write (chns,'(i3.3)') is
          call my_out_us(70,nCv,0,dna0(1,is),'b2tqna_dna0'//chns)
          call my_out_us(70,nCv,0,dpa0(1,is),'b2tqna_dpa0'//chns)
          wrk0 = vsa0(:,is)/(mp*am(is)*pl%na(:,is))
          call my_out_us(70,nCv,0,wrk0,'b2tqna_vsa0'//chns)
!          call my_out_us(70,nCv,0,vla0(1,0,is),'b2tqna_vla0x'//chns)
          call my_out_us(70,nCv,0,vla0(1,1,is),'b2tqna_vla0y'//chns)
          wrk0 = hcib(:,is)/pl%na(:,is)
          call my_out_us(70,nCv,0,wrk0,'b2tqna_hcib'//chns)
        enddo
!        call my_out_us(70,nCv,0,hci0,'b2tqna_in_hci0')
        wrk0=hce0/dv%ne
        call my_out_us(70,nCv,0,wrk0,'b2tqna_hce0')
        call my_out_us(70,nCv,0,sig0,'b2tqna_sig0')
        call my_out_us(70,nCv,0,alf0,'b2tqna_alf0')
      endif                                                           !srv 10.11.09 }
cdpc
* ..return
      ncall_b2tqna = ncall_b2tqna+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqna

      end subroutine b2tqna

#ifdef WG_TODO
      subroutine set_transport(nx, ny, ns, bb, na, te, ti, ni, ne, ne2,
     &  hy, qz, dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0,
     &  wrk0)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &  vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &  hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
*   ..common blocks
*   ..local variables
      integer ncall, n2, ix, iy, is, ixref, iyref
      real (kind=R8) ::
     *  local_te_gradient, required_te_gradient,
     1  d_old, d_new, scaling_factor, ratio
      save ncall, ixref, iyref, required_te_gradient,
     1  d_old, ratio
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external ipgeti, ipgetr, xertst, b2xvsg
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, required_te_gradient/5.0e4_R8/

*   ..subprogram start-up calls
      call subini ('set_transport')
      if (ncall.eq.0) then
        if(nnreg(0).gt.1 .and. nnreg(0).le.5) then
          ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
          iyref=2*topcut(1)/3
        elseif(nnreg(0).eq.8) then
          ixref=(rightcut(1)+rightcut(2))/2
          iyref=2*min(topcut(1),topcut(2))/3
        else
          ixref=3*nx/4
          iyref=ny/2
        endif
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
        d_old=hce0(ixref,iyref)/ne(ixref,iyref)
        ratio=dna0(ixref,iyref,1)/d_old
        write(*,*) 'ratio ',ratio
        call ipgetr ('set_transport_required_te_gradient',
     1   required_te_gradient)
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif

*   ..main
      local_te_gradient=(te(ixref,iyref-1)-te(ixref,iyref+1))/ev/
     1 (0.5_R8*hy1(ixref,iyref-1)+hy1(ixref,iyref)+
     2  0.5_R8*hy1(ixref,iyref+1))
      scaling_factor=local_te_gradient/required_te_gradient
      scaling_factor=min(2.0_R8,max(0.5_R8,scaling_factor))
      d_new=d_old*scaling_factor
      d_new=min(100.0_R8,max(0.01_R8,d_new))
      write(*,*) 'd_old,d_new ',d_old,d_new
      do iy=-1,ny
        do ix=-1,nx
          hce0(ix,iy)=d_new*ne(ix,iy)
cdpc     hci0(ix,iy)=d_new  ! cxpb *ni(ix,iy,0) ?
cxpb     hcib(ix,iy,is)=d_new*na(in,is)
        enddo
      enddo
      do is=0,ns-1
        if(.not.is_neutral(is)) then
          do iy=-1,ny
            do ix=-1,nx
              dna0(ix,iy,is)=d_new*ratio
            enddo
          enddo
        endif
      enddo
      d_old=d_new
*   ..return
      ncall = ncall+1
      call subend ()
      return
      end subroutine set_transport

      subroutine set_transport_srv(nx, ny, ns, bb, na, te, ti, ni, ne,     !srv 12.06.08 {
     &  ne2, hy, qz, dna0, dpa0, vla0, vma0, vsa0, hci0, hcib, hce0,
     &  sig0, alf0, wrk0)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo
     &             , only: cry, vol
      use b2mod_geo_corner
      use b2mod_b2cmpa
      use b2mod_b2cmpt
#ifdef PAR_OPT
      use b2mod_par_opt
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &  qz(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &  vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &  hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1),
     &  vma0(-1:nx,-1:ny,0:1,0:ns-1)
      real (kind=R8) ::
     *  hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
*   ..common blocks
*   ..local variables
      integer ncall, n2, ix, iy, is
      integer, save :: itype = 0, iout = 0                            !srv 24.04.02 05.11.03 15.04.11
      real (kind=R8) :: fdnc(-1:nx,-1:ny,0:ns-1)                      !srv 13.11.02 26.11.04 01.03.07
      real (kind=R8) :: fvsc(-1:nx,-1:ny,0:ns-1),
     &                  ftec(-1:nx,-1:ny), ftic(-1:nx,-1:ny,0:ns-1),  !srv 17.04.06 01.03.07
     &                  fvlc(-1:nx,-1:ny,0:1,0:ns-1)                  !srv 19.04.10
      real (kind=R8), save :: cdf = 0.0e0_R8, alf = 0.0e0_R8,         !srv 24.04.02 05.11.03
     &      alfa1 = 0.0e0_R8, alfa2 = 0.0e0_R8,
     &      delta1 = 0.0e0_R8, delta2 = 0.0e0_R8                      !srv 13.11.02 24.06.03
      real (kind=R8), save :: vs_alfa1 = 1.0e0_R8, vs_alfa2 = 1.0e0_R8,
     &                       vs_delta1 = 1.0e0_R8, vs_delta2= 1.0e0_R8 !srv 17.04.06 {
      real (kind=R8), save :: te_alfa1 = 1.0e0_R8, te_alfa2 = 1.0e0_R8,
     &                       te_delta1 = 1.0e0_R8, te_delta2= 1.0e0_R8
      real (kind=R8), save :: ti_alfa1 = 1.0e0_R8, ti_alfa2 = 1.0e0_R8,
     &                       ti_delta1 = 1.0e0_R8, ti_delta2= 1.0e0_R8 !srv 17.04.06 }
      real (kind=R8), save :: vl_alfa1 = 1.0e0_R8, vl_alfa2 = 1.0e0_R8,
     &                       vl_delta1 = 1.0e0_R8, vl_delta2= 1.0e0_R8 !srv 19.04.10
      real (kind=R8), save ::    alfa3 = 1.0e0_R8                      !srv 13.03.07 {
      real (kind=R8), save :: vs_alfa3 = 1.0e0_R8
      real (kind=R8), save :: te_alfa3 = 1.0e0_R8
      real (kind=R8), save :: ti_alfa3 = 1.0e0_R8                      !srv 13.03.07 }
      real (kind=R8), save :: vl_alfa3 = 1.0e0_R8                      !srv 19.04.10
! { IYS 29.11.2012
      real (kind=R8), save ::    alfa4 = 1.0e0_R8                      !srv 19.04.13
      real (kind=R8), save :: vs_alfa4 = 1.0e0_R8
      real (kind=R8), save :: te_alfa4 = 1.0e0_R8
      real (kind=R8), save :: ti_alfa4 = 1.0e0_R8
      real (kind=R8), save :: vl_alfa4 = 1.0e0_R8
! IYS 29.12.2012}
      character chns*3                                                !srv 05.11.03
      save ncall
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external ipgeti, ipgetr, xertst, b2xvsg
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/

*   ..subprogram start-up calls
      call subini ('set_transport')
      if (ncall.eq.0) then
       call ipgeti ('b2tqna_iout', iout)                          !srv 05.11.03
       call ipgetr ('b2tqna_cdf', cdf)                            !srv 24.04.02
       call ipgetr ('b2tqna_alf', alf)                            !srv 05.11.03
       call ipgeti ('b2tqna_itype', itype)                        !srv 24.04.02
       call ipgetr ('b2tqna_alfa1', alfa1)                        !srv 13.11.02 {
       call xertst (0.le.alfa1,'faulty argument b2tqna_alfa1')
       call ipgetr ('b2tqna_alfa2', alfa2)
       call xertst (0.le.alfa2,'faulty argument b2tqna_alfa2')
       call ipgetr ('b2tqna_alfa3', alfa3)                        !srv 13.03.07
       call xertst (0.le.alfa3,'faulty argument b2tqna_alfa3')    !srv 13.03.07
       alfa4 = alfa3                                              !IYS 29.11.2012
       call ipgetr ('b2tqna_alfa4', alfa4)                        !IYS 29.11.2012
       call xertst (0.le.alfa4,'faulty argument b2tqna_alfa4')    !IYS 29.11.2012
       call ipgetr ('b2tqna_delta1', delta1)
       call ipgetr ('b2tqna_delta2', delta2)                      !srv 24.06.03
       call xertst (0.le.delta1,'faulty argument b2tqna_delta1')  !srv 13.11.02 }
       call xertst (0.le.delta2,'faulty argument b2tqna_delta2')  !srv 24.06.03 }

       vs_alfa1 = alfa1; te_alfa1 = alfa1; ti_alfa1 = alfa1;      !srv 17.04.06 {
       vs_alfa2 = alfa2; te_alfa2 = alfa2; ti_alfa2 = alfa2;
       vs_alfa3 = alfa3; te_alfa3 = alfa3; ti_alfa3 = alfa3;
       vs_alfa4 = alfa4; te_alfa4 = alfa4; ti_alfa4 = alfa4;      !IYS 29.11.2012
       vs_delta1=delta1; te_delta1=delta1; ti_delta1=delta1;
       vs_delta2=delta2; te_delta2=delta2; ti_delta2=delta2;

       vl_alfa1 = alfa1;                                          !srv 19.04.10 {
       vl_alfa2 = alfa2;
       vl_alfa3 = alfa3;
       vl_alfa4 = alfa4;                                          !IYS 29.11.2012
       vl_delta1=delta1;
       vl_delta2=delta2;                                          !srv 19.04.10 }
c                      Connected double-null geometry case
c    Figure shows values of coefficients which are applied to anomalous
c    transport coefficients in regions.
c    alfa1 is applied on cells which have radial coordinates less than delta1 [m] in outmidplane
c    alfa2 is applied on cells which have radial coordinates greater than -delta2 [m] in outmidplane
c
c    +-------+---------------+-------++-------+---------------+-------+
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |       :     alfa3     :       ||       :     alfa3     :       |
c    |       :               :       ||       :               :       |
c    |       :               :       ||       :               :       |
c    |       :...............:       ||       :...............:       |
c    |       :               :       ||       :               :       |
c    |       :     alfa2     :       ||       :     alfa2     :       |
c    |       :               :       ||       :               :       |
c    | alfa4 +---------------+ alfa4 || alfa4 +---------------+ alfa4 |
c    |       |               |       ||       |               |       |
c    |       |     alfa1     |       ||       |     alfa1     |       |
c    |       |               |       ||       |               |       |
c    |       |...............|       ||       |...............|       |
c    |       |               |       ||       |               |       |
c    |       |       1       |       ||       |       1       |       |
c    +-------+---------------+-------++-------+---------------+-------+
c
       call ipgetr ('b2tqna_vs_alfa1', vs_alfa1)
       call xertst (0.le.vs_alfa1,'faulty argument b2tqna_vs_alfa1')
       call ipgetr ('b2tqna_vs_alfa2', vs_alfa2)
       call xertst (0.le.vs_alfa2,'faulty argument b2tqna_vs_alfa2')
       call ipgetr ('b2tqna_vs_alfa3', vs_alfa3)                       !srv 13.03.07
       call xertst (0.le.vs_alfa3,'faulty argument b2tqna_vs_alfa3')   !srv 13.03.07
       vs_alfa4 = vs_alfa3                                             !IYS 29.11.2012
       call ipgetr ('b2tqna_vs_alfa4', vs_alfa4)                       !IYS 29.11.2012
       call xertst (0.le.vs_alfa4,'faulty argument b2tqna_vs_alfa4')   !IYS 29.11.2012
       call ipgetr ('b2tqna_vs_delta1', vs_delta1)
       call ipgetr ('b2tqna_vs_delta2', vs_delta2)
       call xertst (0.le.vs_delta1,'faulty argument b2tqna_vs_delta1')
       call xertst (0.le.vs_delta2,'faulty argument b2tqna_vs_delta2')

       call ipgetr ('b2tqna_te_alfa1', te_alfa1)
       call xertst (0.le.te_alfa1,'faulty argument b2tqna_te_alfa1')
       call ipgetr ('b2tqna_te_alfa2', te_alfa2)
       call xertst (0.le.te_alfa2,'faulty argument b2tqna_te_alfa2')
       call ipgetr ('b2tqna_te_alfa3', te_alfa3)                       !srv 13.03.07
       call xertst (0.le.te_alfa3,'faulty argument b2tqna_te_alfa3')   !srv 13.03.07
       te_alfa4 = te_alfa3                                             !IYS 29.11.2012
       call ipgetr ('b2tqna_te_alfa4', te_alfa4)                       !IYS 29.11.2012
       call xertst (0.le.te_alfa4,'faulty argument b2tqna_te_alfa4')   !IYS 29.11.2012
       call ipgetr ('b2tqna_te_delta1', te_delta1)
       call ipgetr ('b2tqna_te_delta2', te_delta2)
       call xertst (0.le.te_delta1,'faulty argument b2tqna_te_delta1')
       call xertst (0.le.te_delta2,'faulty argument b2tqna_te_delta2')

       call ipgetr ('b2tqna_ti_alfa1', ti_alfa1)
       call xertst (0.le.ti_alfa1,'faulty argument b2tqna_ti_alfa1')
       call ipgetr ('b2tqna_ti_alfa2', ti_alfa2)
       call xertst (0.le.ti_alfa2,'faulty argument b2tqna_ti_alfa2')
       call ipgetr ('b2tqna_ti_alfa3', ti_alfa3)                       !srv 13.03.07
       call xertst (0.le.ti_alfa3,'faulty argument b2tqna_ti_alfa3')   !srv 13.03.07
       ti_alfa4 = ti_alfa3                                             !IYS 29.11.2012
       call ipgetr ('b2tqna_ti_alfa4', ti_alfa4)                       !IYS 29.11.2012
       call xertst (0.le.ti_alfa4,'faulty argument b2tqna_ti_alfa4')   !IYS 29.11.2012
       call ipgetr ('b2tqna_ti_delta1', ti_delta1)
       call ipgetr ('b2tqna_ti_delta2', ti_delta2)
       call xertst (0.le.ti_delta1,'faulty argument b2tqna_ti_delta1')
       call xertst (0.le.ti_delta2,'faulty argument b2tqna_ti_delta2') !srv 17.04.06 }

       call ipgetr ('b2tqna_vl_alfa1', vl_alfa1)  !sign any            !srv 19.04.10 {
       call ipgetr ('b2tqna_vl_alfa2', vl_alfa2)
       call ipgetr ('b2tqna_vl_alfa3', vl_alfa3)
       vl_alfa4 = vl_alfa3                         !IYS 29.11.2012
       call ipgetr ('b2tqna_vl_alfa4', vl_alfa4)   !IYS 29.11.2012
       call ipgetr ('b2tqna_vl_delta1', vl_delta1)
       call ipgetr ('b2tqna_vl_delta2', vl_delta2)
       call xertst (0.le.vl_delta1,'faulty argument b2tqna_vl_delta1')
       call xertst (0.le.vl_delta2,'faulty argument b2tqna_vl_delta2') !srv 19.04.10 }
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif

*   ..main
#ifdef PAR_OPT
      if(npar_opt.eq.16) then
        delta1 = par_opt_phys(9)
        delta2 = par_opt_phys(10)
        alfa1  = par_opt_phys(11)
        alfa2  = par_opt_phys(12)
        alfa3  = par_opt_phys(13)

        vs_alfa1 = alfa1
        vs_alfa2 = alfa2
        vs_alfa3 = alfa3
        vs_delta1=delta1
        vs_delta2=delta2

        te_alfa1 = par_opt_phys(14)
        te_alfa2 = par_opt_phys(15)
        te_alfa3 = par_opt_phys(16)
        te_delta1=delta1
        te_delta2=delta2

        ti_alfa1 =te_alfa1
        ti_alfa2 =te_alfa2
        ti_alfa3 =te_alfa3
        ti_delta1=te_delta1
        ti_delta2=te_delta2
      endif
#endif
*   ..compute special factor for dna0, vsa0, hci0, hce0               !srv 13.11.02
      call defdc(nx, ny, itype, ns, fdnc, fvlc, fvsc, ftec, ftic)     !srv 13.11.02 17.04.06 01.03.07 09.03.17
      if(iout.ne.0 .and. ncall.eq.0) then                             !srv 09.10.01
       do is = 0,ns-1
        write (chns,'(i3.3)') is                                      !srv 09.10.01
        call my_out(70,nx,ny,fdnc(-1,-1,is),'Fdnc'//chns)             !srv 28.06.02
        call my_out(70,nx,ny,fvsc(-1,-1,is),'Fvsc'//chns)             !srv 17.04.06
        call my_out(70,nx,ny,ftic(-1,-1,is),'Ftic'//chns)             !srv 17.04.06
        call my_out(70,nx,ny,fvlc(-1,-1,0,is),'Fvlcx'//chns)          !srv 19.04.10
        call my_out(70,nx,ny,fvlc(-1,-1,1,is),'Fvlcy'//chns)          !srv 19.04.10
       enddo
       call my_out(70,nx,ny,ftec,'Ftec')                              !srv 17.04.06
       call my_out(70,nx,ny,sig0,'sig0')
      endif
      do iy=-1,ny
        do ix=-1,nx
          hce0(ix,iy) = hce0(ix,iy)*ftec(ix,iy)
        enddo
      enddo
      do iy=-1,ny
        do ix=-1,nx
          hci0(ix,iy) = 0.0_R8
          do is=0,ns-1
            dna0(ix,iy,is) = dna0(ix,iy,is)*fdnc(ix,iy,is)
            vla0(ix,iy,:,is) = vla0(ix,iy,:,is)*fvlc(ix,iy,:,is)      !srv 19.04.10
            vsa0(ix,iy,is) = vsa0(ix,iy,is)*fvsc(ix,iy,is)
            hcib(ix,iy,is) = hcib(ix,iy,is)*ftic(ix,iy,is)            !srv 19.05.03
            hci0(ix,iy)    = hci0(ix,iy) + hcib(ix,iy,is)
          enddo
        enddo
      enddo

      if(iout.ne.0) then                             !srv 15.06.08 {
       do is = 0,ns-1
        write (chns,'(i3.3)') is
        call my_out(70,nx,ny,dna0(-1,-1,is),'dna0'//chns)
        call my_out(70,nx,ny,dpa0(-1,-1,is),'dpa0'//chns)
        call my_out(70,nx,ny,vsa0(-1,-1,is),'vsa0'//chns)
        call my_out(70,nx,ny,vla0(-1,-1,0,is),'vla0x'//chns)
        call my_out(70,nx,ny,vla0(-1,-1,1,is),'vla0y'//chns)
        call my_out(70,nx,ny,vma0(-1,-1,0,is),'vma0x'//chns)
        call my_out(70,nx,ny,vma0(-1,-1,1,is),'vma0y'//chns)
        call my_out(70,nx,ny,hcib(-1,-1,is),'hcib'//chns)
       enddo
       call my_out(70,nx,ny,hci0,'hci0')
       call my_out(70,nx,ny,hce0,'hce0')
       call my_out(70,nx,ny,sig0,'sig0')
       call my_out(70,nx,ny,alf0,'alf0')
      endif                                                           !srv 15.06.08

*   ..return
      ncall = ncall+1
      call subend ()
      return

      contains                          !srv 24.04.02 {

       function barier1(x, x1, level1, x2, level2, x3, level3)
        implicit none
        real (kind=R8) :: x, x1, level1, x2, level2, x3, level3
        real (kind=R8) :: barier1
*   ..   x1 < x2 < x3
         if(x .le. x1) then
          barier1 = level1
         elseif(x .le. x2) then
          barier1 = level1+(x - x1)*(level2 - level1)/(x2 - x1)
         elseif(x .le. x3) then
          barier1 = level2
         else
          barier1 = level3
         endif
       end function barier1

       function barier2(x, x1, level1, x2, level2, x3, level3)
        implicit none
        real (kind=R8) :: x, x1, level1, x2, level2, x3, level3
        real (kind=R8) :: barier2
*   ..   x1 > x2 > x3
         if(x .le. x1) then
          barier2 = level1
         elseif(x .le. x2) then
          barier2 = level2
         elseif(x .le. x3) then
          barier2 = level2+(x - x2)*(level3 - level2)/(x3 - x2)
         else
          barier2 = level3
         endif
       end function barier2

       subroutine defdc(nx, ny, itype, ns, fdnc, fvlc, fvsc, ftec, ftic) !srv 13.11.02 01.03.07 09.03.17
        implicit none
        integer :: nx, ny, itype, ns
        real(kind=R8) :: fdnc(-1:nx,-1:ny,0:ns-1)                    !srv 17.04.06 {
        real(kind=R8) :: fvlc(-1:nx,-1:ny,0:1,0:ns-1)                !srv 09.03.17
        real(kind=R8) :: fvsc(-1:nx,-1:ny,0:ns-1)
        real(kind=R8) :: ftec(-1:nx,-1:ny), ftic(-1:nx,-1:ny,0:ns-1) !srv 17.04.06 }

        logical :: exist_file(0:2)
        character*256 :: filename(0:2)
        integer ix, iy, is, nc_c_m_f(1)                           !srv 05.07.04
        integer, save :: ncall = 0, bottomdelta1, topdelta2, mintopcut !srv 28.07.10
        real (kind=R8) :: width, bbaverage, t0, t1, step
        real (kind=R8), allocatable, save :: coord(:), fun(:,:)

        fdnc = 1.0e0_R8
        fvsc = 1.0e0_R8                                           !srv 17.04.06 {
        ftec = 1.0e0_R8
        ftic = 1.0e0_R8                                           !srv 17.04.06 }
        fvlc = 1.0e0_R8                                           !srv 19.04.10 }

       if(nncut > 0) then                                         !srv 25.02.16
        select case (itype)
         case(9)                                                !srv 22.05.09 {
          mintopcut=minval(topcut(1:nncut))
*   ..modification of case (5). It is used barrier function
          do iy = -1, ny                  ! for SOL and PR region
           do ix = -1, nx
            if(.not.on_closed_surface(ix,iy)) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is) =    alfa3
               fvsc(ix,iy,is) = vs_alfa3
               ftic(ix,iy,is) = ti_alfa3
              endif
             enddo
             ftec(ix,iy) = te_alfa3
            endif
           enddo
          enddo
          step = maxval(hy(nmdpl,mintopcut:ny))
          do iy = mintopcut,ny          ! for SOL region
           do ix=-1,nx
            if(mod(region(ix,iy,0),4).eq.2) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then

               fdnc(ix,iy,is)= barier2(crdy_c_cell(nmdpl,iy),
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8, delta2-step,
     &          alfa2, delta2, alfa3)
               fvsc(ix,iy,is)= barier2(crdy_c_cell(nmdpl,iy),
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8, vs_delta2-step,
     &          vs_alfa2, vs_delta2, vs_alfa3)
               ftic(ix,iy,is)= barier2(crdy_c_cell(nmdpl,iy),
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8, ti_delta2-step,
     &          ti_alfa2, ti_delta2, ti_alfa3)
              endif
             enddo
             ftec(ix,iy)= barier2(crdy_c_cell(nmdpl,iy),
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8, te_delta2-step,
     &          te_alfa2, te_delta2, te_alfa3)
            endif
           enddo
          enddo

          step = maxval(hy(nmdpl,-1:mintopcut-1))
          do iy = -1, mintopcut-1   ! for inner surfaces
           do ix=-1,nx
            if(on_closed_surface(ix,iy)) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is)= barier1(crdy_c_cell(nmdpl,iy),
     &          -delta1, 1.0_R8, -delta1+step, alfa1,
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8)
               fvsc(ix,iy,is)= barier1(crdy_c_cell(nmdpl,iy),
     &          -vs_delta1, 1.0_R8, -vs_delta1+step, vs_alfa1,
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8)
               ftic(ix,iy,is)= barier1(crdy_c_cell(nmdpl,iy),
     &          -ti_delta1, 1.0_R8, -ti_delta1+step, ti_alfa1,
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8)
              endif
             enddo
             ftec(ix,iy)= barier1(crdy_c_cell(nmdpl,iy),
     &          -te_delta1, 1.0_R8, -te_delta1+step, te_alfa1,
     &          crdy_surf(nmdpl,mintopcut), 1.0_R8)
            endif
           enddo
          enddo                                                   !srv 2.05.09 }

         case (8)                                                 !srv 01.03.07 {
*   ..setting of transport coefficients from files
          if(ncall.eq.0) then
           if(nncut.eq.2) then
            nmdpl = rightcut(nncut)
            do ix = rightcut(nncut)+1, rightcut(nncut-1)-1
             if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix   !srv 11.01.13
            enddo
           elseif(nncut.eq.1) then
            nmdpl = leftcut(nncut)
            do ix = leftcut(nncut)+1, rightcut(nncut)-1
             if(abs(bb(ix,ny,2)) .lt. abs(bb(nmdpl,ny,2))) nmdpl = ix   !srv 11.01.13
            enddo
           else
            write(*,*) ' b2tqna: nncut not equal 1 or 2'
            call xerrab ('Not available outside SN and DN geometries')
           endif
           write(*,"(/' Parameters of H-regime')")
           write(*,"(' Number of midplane cell =',i4)") nmdpl
           call xertst(0.0e0_R8.le.delta1,
     &                 'faulty argument b2tqna_delta1') !srv 15.05.09
           call xertst(0.0e0_R8.le.delta2,
     &                 'faulty argument b2tqna_delta2') !srv 15.05.09

           allocate(coord(-1:ny), fun(-1:ny,0:2))

*          shift coordinates in order to separatrix is located in zero
           coord(-1:ny) = lengthy(nmdpl,-1:ny)-lengthy(nmdpl,topcut(1))

           filename(0) = 'b2tdnc.dat'
           filename(1) = 'b2thic.dat'
           filename(2) = 'b2thec.dat'
           call find_file(filename(0),exist_file(0))
           call find_file(filename(1),exist_file(1))
           call find_file(filename(2),exist_file(2))
           call xertst (all(exist_file),
     1      'There are no b2tdnc.dat or b2thic.dat or b2thec.dat files')

           call set_profile(filename(0), ny, coord, fun(-1,0))
           do iy = -1, ny
             call xertst (fun(iy,0).gt. 0.0e0_R8,
     1        'Approximation of data from file '//trim(filename(0))//
     2        ' has given '//'negative coefficient')
           enddo

           call set_profile(filename(1), ny, coord, fun(-1,1))
           do iy = -1, ny
             call xertst (fun(iy,1).gt. 0.0e0_R8,
     1        'Approximation of data from file '//trim(filename(1))//
     2        ' has given '//'negative coefficient')
           enddo

           call set_profile(filename(2), ny, coord, fun(-1,2))
           do iy = -1, ny-1
             call xertst (fun(iy,2).gt. 0.0e0_R8,
     1        'Approximation of data from file '//trim(filename(2))//
     2        ' has given '//'negative coefficient')
           enddo

           write(*,"(/' Exchanges of transport coefficients',
     &               ' have been done'/)")
           do is = 0, ns-1
             if(.not.is_neutral(is)) then
               cfdna(2,is)   = 1.0e0_R8
               cfhci(2,is)   = 1.0e0_R8
               write(*,"(' Now cfdna(2,',i3,') =',es10.2)")
     &          is, cfdna(2,is)
               write(*,"(' Now cfhci(2,',i3,') =',es10.2)")
     &          is, cfhci(2,is)
             endif
           enddo

           cfhce(0) = 1.0e0_R8
           write(*,"( ' Now cfhce(0)    =',es10.2/)") cfhce(0)

          endif

          do is = 0, ns-1
           if(.not.is_neutral(is)) then
            do ix = -1,nx
             fdnc(ix,:,is)  = fun(:,0)
             ftic(ix,:,is) = fun(:,1)    ! ????
            enddo
           endif
          enddo

          do ix = -1,nx
           ftec(ix,:) = fun(:,2)
          enddo

          nc_c_m_f(:)=maxloc(cry(-1:nx,-1,0))-2 ! number of cell for the centre of magnetic axis
          write(*,*)' number of cell for the centre of magnetic axis =',
     *              nc_c_m_f(1)
          do iy = -1, ny
           do ix = -1, nx
            t0 = bb(nc_c_m_f(1),iy,2)/bb(ix,iy,2)
            do is = 0, ns-1
             if(.not.is_neutral(is)) then
              fdnc(ix,iy,is)  = fdnc(ix,iy,is)*t0
              ftic(ix,iy,is) = ftic(ix,iy,is)*t0
             endif
            enddo
            ftec(ix,iy) = ftec(ix,iy)*t0
           enddo
          enddo

         case (7)                                                  !srv 05.07.04 {
*   ..calculate
          do is = 0, ns-1
           if(.not.is_neutral(is)) then
            t0 = 0.0e0_R8
            t1 = 0.0e0_R8
            mintopcut=min(topcut(1)-1,topcut(2)-1)
            do iy = 0,mintopcut    ! for inner surfaces
             do ix=0,nx
              if(on_closed_surface(ix,iy)) then
               t0 = t0 +bb(ix,iy,3)*vol(ix,iy)/hy1(ix,iy)
               t1 = t1 +vol(ix,iy)/hy1(ix,iy)
              endif
             enddo
            enddo
            bbaverage=t0/t1
            do iy = -1, ny                   ! for all surfaces
             do ix = -1, nx
              if(on_closed_surface(ix,iy) .and. iy.le.mintopcut) then
               fdnc(ix,iy,is)=1.0e0_R8
              else
               fdnc(ix,iy,is)=1.0e0_R8+cdf*(bbaverage/bb(ix,iy,3))**alf
              endif
             enddo
            enddo
           endif                                                     !srv 05.07.04 }
          enddo
         case (6)                                                  !srv 05.11.03 {

*   ..calculate poloidal variation of transport coefficients
          do is = 0, ns-1
           if(.not.is_neutral(is)) then
            if(alf .ne. 0.0e0_R8) then
             t0 = 0.0e0_R8
             t1 = 0.0e0_R8
             do iy = 0, topcut(1)-1   ! for inner surfaces
!           do ix = leftcut(1), rightcut(1)-1
              do ix=0,nx
               if(on_closed_surface(ix,iy)) then
                t0 = t0 +bb(ix,iy,3)**(-alf)*vol(ix,iy)/hy1(ix,iy)**2
                t1 = t1 +vol(ix,iy)/hy1(ix,iy)**2
               endif
              enddo
             enddo
             bbaverage = (t0/t1)**(-1.0e0_R8/alf)
             do iy = -1, ny                   ! for all surfaces
              do ix = -1, nx
               fdnc(ix,iy,is) = (bbaverage/bb(ix,iy,3))**alf
              enddo
             enddo
            else
             fdnc(:,:,is) = 1.0e0_R8
            endif
           endif                                     !srv 05.11.03 }
          enddo
         case (5)                                               !srv 13.11.02 12.07.10 {

          if(ncall.eq.0) then
           mintopcut=minval(topcut(1:nncut))
           call xertst(0.0e0_R8.le.delta1,
     *         'faulty argument b2tqna_delta1') !srv 15.05.09
           call xertst(0.0e0_R8.le.delta2,
     *         'faulty argument b2tqna_delta2') !srv 15.05.09
           if(delta1.eq.0) then
            bottomdelta1 = mintopcut
           elseif(delta1.lt.lengthy(nmdpl,mintopcut)) then              !srv 15.05.09 {
            bottomdelta1 = mintopcut
            width = 0.0e0_R8
            do iy =mintopcut-1,-1,-1               ! srv 15.05.09 was topcut(1)-1,0,-1
             width = width + hy(nmdpl,iy)
             if(width.ge.delta1) then
              bottomdelta1 = iy
              exit
             endif
            enddo
           else
            bottomdelta1 = -1
           endif

           if(delta2.eq.0) then
            topdelta2 = mintopcut-1
           elseif(delta2.lt.
     &           lengthy(nmdpl,ny)-lengthy(nmdpl,mintopcut)) then
            topdelta2 = mintopcut                 !srv 24.06.03 {
            width = 0.0e0_R8
            do iy =mintopcut+1,ny
             width = width + hy(nmdpl,iy)
             if(width.ge.delta2) then
              topdelta2 = iy
              exit
             endif
            enddo                                !srv 24.06.03 {
           else
            topdelta2 = ny
           endif                                                       !srv 15.05.09 }
           if(iout.ne.0) then
            open(70,file= 'Parameters_of_H-regime.dat')
            write(70,"(/' Parameters of H-regime')")
            write(70,"(' Number of midplane cell =',i4)") nmdpl
            write(70,"(' Number of first cell in CORE where transport ',
     &         ' coefficients are changed',2i5)") bottomdelta1
            write(70,"(' Number of last  cell in CORE where transport ',
     &         ' coefficients are changed',2i5)") topdelta2
            close(70)
           endif
          endif

          do iy = -1, ny                  ! for SOL and PR region  !srv 13.03.07
           do ix = -1, nx
            if(.not.on_closed_surface(ix,iy)) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is)  =    alfa3
               fvsc(ix,iy,is) = vs_alfa3
               ftic(ix,iy,is) = ti_alfa3
               fvlc(ix,iy,1,is) = vl_alfa3                      !srv 19.04.10
              endif
             enddo
             ftec(ix,iy) = te_alfa3
            endif
           enddo
          enddo

          do iy = mintopcut,topdelta2          ! for SOL region !srv 24.06.03
           do ix=-1,nx                                          !srv 10.03.06
            if(mod(region(ix,iy,0),4).eq.2) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is)  = alfa2
               fvsc(ix,iy,is) = vs_alfa2
               ftic(ix,iy,is) = ti_alfa2
               fvlc(ix,iy,1,is) = vl_alfa2                        !srv 19.04.10
              endif
             enddo
             ftec(ix,iy)    = te_alfa2
            endif
           enddo
          enddo

          do iy = bottomdelta1, mintopcut-1   ! for inner surfaces
           do ix=-1,nx                                          !srv 10.03.06
            if(on_closed_surface(ix,iy)) then
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is) = alfa1
               fvsc(ix,iy,is) = vs_alfa1
               ftic(ix,iy,is) = ti_alfa1
               fvlc(ix,iy,1,is) = vl_alfa1                        !srv 19.04.10
              endif
             enddo
             ftec(ix,iy) = te_alfa1
            endif
           enddo
          enddo

          do iy = -1, ny                  ! correction for   !IYS 29.11.12
           do ix = -1, nx
            if(mod(region(ix,iy,0),4).ne.1 .and.    ! NOT Core
     &         mod(region(ix,iy,0),4).ne.2) then    ! NOT SOL
             do is = 0, ns-1
              if(.not.is_neutral(is)) then
               fdnc(ix,iy,is)   =    alfa4
               fvsc(ix,iy,is)   = vs_alfa4
               ftic(ix,iy,is)   = ti_alfa4
               fvlc(ix,iy,1,is) = vl_alfa4
              endif
             enddo
             ftec(ix,iy) = te_alfa4
            endif
           enddo
          enddo

*   ..take into account poloidal variation of transport coefficients
          if(alf .ne. 0.0e0_R8) then                                     !srv 04.07.06 {
           t0 = 0.0e0_R8
           t1 = 0.0e0_R8
           do iy = 0, mintopcut-1   ! for inner surfaces
            do ix=0,nx
             if(on_closed_surface(ix,iy)) then
              t0 = t0 +bb(ix,iy,3)**(-alf)*vol(ix,iy)/hy1(ix,iy)**2
              t1 = t1 +vol(ix,iy)/hy1(ix,iy)**2
             endif
            enddo
           enddo
           bbaverage = (t0/t1)**(-1.0e0_R8/alf)

           do iy = topdelta2+1, ny                    ! for SOL region
            do ix = -1,nx
             if(mod(region(ix,iy,0),4).eq.2) then
              do is = 0, ns-1
               if(.not.is_neutral(is)) then
                fdnc(ix,iy,is)  = fdnc(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                fvsc(ix,iy,is) = fvsc(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                ftic(ix,iy,is) = ftic(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                fvlc(ix,iy,1,is) = fvlc(ix,iy,1,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf         !srv 19.04.10
               endif
              enddo
              ftec(ix,iy) = ftec(ix,iy)*(bbaverage/bb(ix,iy,3))**alf
             endif
            enddo
           enddo

           do iy = -1, bottomdelta1-1                 ! for inner surfaces
            do ix = -1,nx
             if(on_closed_surface(ix,iy)) then
              do is = 0, ns-1
               if(.not.is_neutral(is)) then
                fdnc(ix,iy,is)  = fdnc(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                fvsc(ix,iy,is) = fvsc(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                ftic(ix,iy,is) = ftic(ix,iy,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf
                fvlc(ix,iy,1,is) = fvlc(ix,iy,1,is)*
     &                           (bbaverage/bb(ix,iy,3))**alf         !srv 19.04.10
               endif
              enddo
             ftec(ix,iy) = ftec(ix,iy)*(bbaverage/bb(ix,iy,3))**alf
             endif
            enddo
           enddo
          endif                                                          !srv 04.07.06 12.07.10 }

         case (4)

          do iy = topcut(1),ny              ! for SOL region
           do ix = -1,nx
            fdnc(ix,iy,:) = fdnc(ix,iy,:) + cdf            !srv 01.03.07 ???
           enddo
          enddo

         case (3)

          do iy = topcut(1),ny              ! for SOL region
           do ix = -1,nx
            fdnc(ix,iy,:) = fdnc(ix,iy,:) + cdf*                !srv 01.03.07 ???
     *                 sqrt((lengthy(nx/2,iy)-lengthy(nx/2,topcut(1)))
     *                     /(lengthy(nx/2,ny)-lengthy(nx/2,topcut(1))))

           enddo
          enddo

         case (2)

          do iy = topcut(1),ny              ! for SOL region
           do ix = -1,nx
            fdnc(ix,iy,:) = fdnc(ix,iy,:) +               !srv 01.03.07 ???
     *                 cdf*((lengthy(nx/2,iy)-lengthy(nx/2,topcut(1)))
     *                 /(lengthy(nx/2,ny)-lengthy(nx/2,topcut(1))))**2
           enddo
          enddo

         case (1)

          do iy = topcut(1),ny              ! for SOL region
           do ix = -1,nx
            fdnc(ix,iy,:) = fdnc(ix,iy,:) +               !srv 01.03.07 ???
     *                  cdf*((lengthy(nx/2,iy)-lengthy(nx/2,topcut(1)))
     *                     /(lengthy(nx/2,ny)-lengthy(nx/2,topcut(1))))
           enddo
          enddo
         case default

          fdnc = 1.0e0_R8
        end select
       endif                                                !srv 25.02.16

        ncall = ncall + 1

      end subroutine defdc                  !srv 24.02.04 13.11.02}

      end subroutine set_transport_srv                            !srv 12.06.08 }

      subroutine set_transport_swy(nx, ny, ns,
     & dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0)
      use b2mod_types
      use b2mod_rates
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
*   ..in/out arguments
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &   hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..local variables
      integer ncall, ix, iy, is
      integer iy1, iy2, ix1, ix2, elm_on, elm_bal, elm_ix1, elm_ix2    !swy
      real (kind=R8) ::                                                !swy
     *     fac0_dna0, fac1_dna0, fac2_dna0,                            !swy
     1     fac0_hce0, fac1_hce0, fac2_hce0,                            !swy
     2     fac0_vla0, fac1_vla0, fac2_vla0,                            !swy
     3     elm_t1, elm_t2, elm_fac0, elm_fac1, elm_fac2                !swy
      save ncall,
     *     fac0_dna0, fac1_dna0, fac2_dna0,
     1     fac0_hce0, fac1_hce0, fac2_hce0,
     2     fac0_vla0, fac1_vla0, fac2_vla0,
     3     elm_t1, elm_t2, elm_fac0, elm_fac1, elm_fac2
*   ..initialisation
      data ncall/0/,
     2 fac0_dna0/1.0_R8/, fac1_dna0/1.0_R8/, fac2_dna0/1.0_R8/,        !swy
     4 fac0_hce0/1.0_R8/, fac1_hce0/1.0_R8/, fac2_hce0/1.0_R8/,        !swy
     5 fac0_vla0/1.0_R8/, fac1_vla0/1.0_R8/, fac2_vla0/1.0_R8/ ,       !swy
     6 elm_on/0/, elm_t1/-1.0_R8/, elm_t2/-1.0_R8/,                    !swy
     7 elm_fac0/1.0_R8/, elm_fac1/1.0_R8/, elm_fac2/1.0_R8/,           !swy
     8 elm_bal/0/, elm_ix1/-1/, elm_ix2/-1/                            !swy

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('set_transport_swy')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tqna_iy1', iy1)
        call ipgeti ('b2tqna_iy2', iy2)
        call ipgetr ('b2tqna_fac0_dna0', fac0_dna0)
        call ipgetr ('b2tqna_fac1_dna0', fac1_dna0)
        call ipgetr ('b2tqna_fac2_dna0', fac2_dna0)
        call ipgetr ('b2tqna_fac0_hce0', fac0_hce0)
        call ipgetr ('b2tqna_fac1_hce0', fac1_hce0)
        call ipgetr ('b2tqna_fac2_hce0', fac2_hce0)
        call ipgetr ('b2tqna_fac0_vla0', fac0_vla0)
        call ipgetr ('b2tqna_fac1_vla0', fac1_vla0)
        call ipgetr ('b2tqna_fac2_vla0', fac2_vla0)
        call ipgeti ('b2tqna_ix1', ix1)
        call ipgeti ('b2tqna_ix2', ix2)

        call ipgeti ('b2tqna_elm_on', elm_on)
        call ipgetr ('b2tqna_elm_t1', elm_t1)
        call ipgetr ('b2tqna_elm_t2', elm_t2)
        call ipgetr ('b2tqna_elm_fac0', elm_fac0)
        call ipgetr ('b2tqna_elm_fac1', elm_fac1)
        call ipgetr ('b2tqna_elm_fac2', elm_fac2)

        call ipgeti ('b2tqna_elm_bal', elm_bal)
        call ipgeti ('b2tqna_elm_ix1', elm_ix1)
        call ipgeti ('b2tqna_elm_ix2', elm_ix2)

      endif

*--------setting transport coefficients factor
      do iy=-1, ny
        do ix=ix1, ix2
          if(iy.ge.-1.and.iy.lt.iy1) then
            dna0(ix,iy,:)=fac0_dna0 * dna0(ix,iy,:)
            hce0(ix,iy)=fac0_hce0 * hce0(ix,iy)
            hcib(ix,iy,:)=fac0_hce0 * hcib(ix,iy,:)
            vla0(ix,iy,:,:)=fac0_vla0 * vla0(ix,iy,:,:)
          endif
          if(iy.ge.iy1.and.iy.le.iy2) then
            dna0(ix,iy,:)=fac1_dna0 * dna0(ix,iy,:)
            hce0(ix,iy)=fac1_hce0 * hce0(ix,iy)
            hcib(ix,iy,:)=fac1_hce0 * hcib(ix,iy,:)
            vla0(ix,iy,:,:)=fac1_vla0 * vla0(ix,iy,:,:)
          endif
          if(iy.gt.iy2.and.iy.le.ny) then
            dna0(ix,iy,:)=fac2_dna0 * dna0(ix,iy,:)
            hce0(ix,iy)=fac2_hce0 * hce0(ix,iy)
            hcib(ix,iy,:)=fac2_hce0 * hcib(ix,iy,:)
            vla0(ix,iy,:,:)=fac2_vla0 * vla0(ix,iy,:,:)
          endif
        enddo !end ix loop
      enddo !end iy loop
      do ix=-1,nx
        do iy=-1,ny
          hci0(ix,iy)=0.0_R8
          do is=0,ns-1
            hci0(ix,iy)=hci0(ix,iy)+hcib(ix,iy,is)
          enddo
        enddo
      enddo

*-------applying ELM characteristics
      if(elm_on.eq.1) then
        write(*,*) 'DURING ELM', 'time', tim
        if(tim.ge.elm_t1.and.tim.le.elm_t2) then
          if(elm_bal.ne.1) then
            elm_ix1=ix1
            elm_ix2=ix2
          endif
          do iy=-1, ny
            do ix=elm_ix1, elm_ix2
              if(iy.ge.-1.and.iy.lt.iy1) then
                dna0(ix,iy,:)=elm_fac0 * dna0(ix,iy,:)
                hce0(ix,iy)=elm_fac0 * hce0(ix,iy)
                hcib(ix,iy,:)=elm_fac0 * hcib(ix,iy,:)
              endif
              if(iy.ge.iy1.and.iy.le.iy2) then
                dna0(ix,iy,:)=elm_fac1 * dna0(ix,iy,:)
                hce0(ix,iy)=elm_fac1 * hce0(ix,iy)
                hcib(ix,iy,:)=elm_fac1 * hcib(ix,iy,:)
              endif
              if(iy.gt.iy2.and.iy.le.ny) then
                dna0(ix,iy,:)=elm_fac2 * dna0(ix,iy,:)
                hce0(ix,iy)=elm_fac2 * hce0(ix,iy)
                hcib(ix,iy,:)=elm_fac2 * hcib(ix,iy,:)
              endif
            enddo !end ix loop
          enddo !end iy loop
          do ix=-1,nx
            do iy=-1,ny
              hci0(ix,iy)=0.0_R8
              do is=0,ns-1
                hci0(ix,iy)=hci0(ix,iy)+hcib(ix,iy,is)
              enddo
            enddo
          enddo
        endif
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

      end subroutine set_transport_swy

      subroutine set_transport_eta(nx, ny, ns, bb, na,
     1  te, ti, ni, ne, ne2, hy,
     &  dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_transport_models
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), hy(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &  vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &  hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..common blocks
*   ..local variables
      integer ncall, n2, ix, iy, ixref, iyref
      real (kind=R8) ::
     *   eta, current_eta, eta_floor, eta_ceiling, eta_alpha
      logical models_found
      save ncall, ixref, iyref, eta, eta_floor, eta_ceiling, eta_alpha
*   ..procedures
      external ipgeti, ipgetr, xertst, b2xvsg
*   ..initialisation
      data ncall/0/, eta /2.0_R8/,
     1 eta_floor/0.1_R8/, eta_ceiling/10.0_R8/, eta_alpha/0.5_R8/

*   ..subprogram start-up calls
      call subini ('set_transport_eta')
      if (ncall.eq.0) then
        call alloc_b2mod_transport_models(nx, ny, ns)
        call read_b2mod_transport_models(99,'    ',models_found)
        if(nnreg(0).gt.1 .and. nnreg(0).le.5) then
          ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
          iyref=2*topcut(1)/3
        elseif(nnreg(0).eq.8) then
          ixref=(rightcut(1)+rightcut(2))/2
          iyref=2*min(topcut(1),topcut(2))/3
        else
          ixref=3*nx/4
          iyref=ny/2
        endif
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
        write(*,*) 'set_transport_eta: ixref, iyref = ',ixref, iyref
        call ipgetr ('set_transport_eta',
     1   eta)
        call ipgetr ('set_transport_eta_floor',
     1   eta_floor)
        call ipgetr ('set_transport_eta_ceiling',
     1   eta_ceiling)
        call ipgetr ('set_transport_eta_alpha',
     1   eta_alpha)
        call xertst (0.0_R8.le.eta_alpha.and.eta_alpha.le.1.0_R8,
     1   'faulty argument set_transport_eta_alpha')
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, ne2, 1, 'ne2', '.gt.')
      endif

*   ..main
      do iy=-0,ny-1
        current_eta=
     1   ((te(ixref,iy-1)-te(ixref,iy+1))/
     2    (te(ixref,iy-1)+te(ixref,iy+1)))/
     3   ((ne(ixref,iy-1)-ne(ixref,iy+1))/
     4    (ne(ixref,iy-1)+ne(ixref,iy+1)))
        eta_hce_mult(iy)=eta_hce_mult(iy)**(1.0_R8-eta_alpha)*
     1   (min(eta_ceiling,
     2   max(eta_floor,eta_hce_mult(iy)*current_eta/eta)))**eta_alpha
        write(*,*) 'set_transport_eta ',
     1   iy,eta,current_eta,eta_hce_mult(iy)
      enddo

      do iy=-1,ny
        do ix=-1,nx
          if(on_closed_surface(ix,iy) .or.
     .       mod(region(ix,iy,0),4).eq.2) then
            hce0(ix,iy)=hce0(ix,iy)*eta_hce_mult(iy)
          endif
        enddo
      enddo
*   ..return
      ncall = ncall+1
      call subend ()
      return
      end subroutine set_transport_eta

      subroutine set_transport_xpb(nx, ny, ns,
     & bb, na, ne, te, ti, ni, hx, hy, qz, dna0, vla0)
      use b2mod_types
      use b2mod_rates
      use b2mod_indirect
      use b2mod_constants
      use b2mod_geo
     & , only : crx
      use b2mod_b2cmpa
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1), ne(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1)
*   ..input/output arguments (changed on exit)
      real (kind=R8) ::
     *  dna0(-1:nx,-1:ny,0:ns-1),
     &  vla0(-1:nx,-1:ny,0:1,0:ns-1)
*   ..common blocks
*   ..local variables
      integer, save ::  ncall, n2, ix, iy, is, ir, jxi, jxa, jsep
      real (kind=R8), save ::
     *  alpha,
     *  wb, seglen, perimeter, L_perp
      real (kind=R8) ::
     *  pe(-1:nx,-1:ny), pri(-1:nx,-1:ny), pp(-1:nx,-1:ny),
     *  coss(-1:nx,-1:ny), R, width(-1:nx), rzb(-1:nx,-1:ny)
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external subini, subend, ipgetr, xertst, xerrab
      external b2xvsg, b2xppe, b2xppb
      intrinsic sin, log, mod, max, min
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, alpha/1.0_R8/

cxpb
cxpb  Apply a space-varying pinch velocity according to
cxpb  V. Naulin, Phys. Rev. E 71, 015402(R), (2005).
cxpb  The poloidal coordinate is made to run from X-point
cxpb  to X-point in order to avoid strong discontinuities
cxpb  at the divertor throats, so that the cos(s) factor
cxpb  vanishes there.
cxpb  The model is applied inside the separatrix and an
cxpb  effective diffusion coefficient is used outside.
cxpb

*   ..subprogram start-up calls
      call subini ('set_transport_xpb')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*ns, dna0, 1, 'dna0', '.ge.')
      endif

*    ..initial set-up
      if (ncall.eq.0) then
        write(*,*) 'Using set_transport_xpb'
        call ipgetr ('set_transport_alpha', alpha)
        call get_jsep(nx,ny,jxi,jxa,jsep)
      endif
*   ..compute plasma pressure from charged species only
      call b2xppe (nx, ny, ne, te, pe)
      rzb = 0.0_R8
      call b2xppb (nx, ny, rzb, ni(-1,-1,1), te, ti, pri)
      pp = pe + pri      ! plasma pressure
*   ..compute geometrical quantities
      coss = 0.0_R8
      L_perp = 0.0_R8
      width = 0.0_R8
      perimeter = 0.0_R8
      R = 0.0_R8
      ir = 0
      do ix = -1, nx
        if (on_closed_surface(ix,0)) then
          ir = ir + 1
          perimeter = perimeter + hx(ix,0)
          R = R + 0.5_R8*(crx(ix,0,0)+crx(ix,0,1))*hx(ix,0)
          width(ix) = 0.5_R8 * hy1(ix,0)
          do iy = 0, jsep
            width(ix) = width(ix) + 0.5_R8 * (hy1(ix,iy-1)+hy1(ix,iy))
          end do
          width(ix) = width(ix) + 0.5_R8 * hy1(ix,jsep)
          if (pp(ix,0).gt.pp(ix,jsep)) then
            L_perp = L_perp + width(ix) / log(pp(ix,0)/pp(ix,jsep))
          else
            L_perp = L_perp + width(ix)
          endif
        end if
      end do
      R = R / perimeter
      L_perp = L_perp / ir
cxpb  Compute cos(s) factor
      if (nnreg(0).eq.8) then
         do iy = -1, jsep
cxpb  Compute inner SOL part (coss > 0 by convention)
            coss(leftcut(1),iy) = 0.5_R8 * hx(leftcut(1),iy)/
     1         (sqrt(bb(leftcut(1),iy,0)**2+bb(leftcut(1),iy,1)**2)/
     2          bb(leftcut(1),iy,3))
            do ix = leftcut(1)+1, leftcut(2)-1
              coss(ix,iy) = coss(ix-1,iy) + 0.5_R8 * hx(ix-1,iy)/
     1         (sqrt(bb(ix-1,iy,0)**2+bb(ix-1,iy,1)**2)/bb(ix-1,iy,3))
     2         + 0.5_R8 * hx(ix,iy)/
     3         (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
            end do
            seglen = coss(leftcut(2)-1,iy)
     .       + 0.5_R8 * hx(leftcut(2)-1,iy)/
     1       (sqrt(bb(leftcut(2)-1,iy,0)**2+bb(leftcut(2)-1,iy,1)**2)/
     2       bb(leftcut(2)-1,iy,3))
            coss(leftcut(1):leftcut(2)-1,iy) =
     .       sin(pi*coss(leftcut(1):leftcut(2)-1,iy)/seglen)
cxpb  Compute outer SOL part (coss < 0 by convention)
            coss(rightcut(2),iy) = 0.5_R8 * hx(rightcut(2),iy)/
     1         (sqrt(bb(rightcut(2),iy,0)**2+bb(rightcut(2),iy,1)**2)/
     2          bb(rightcut(2),iy,3))
            do ix = rightcut(2)+1, rightcut(1)-1
              coss(ix,iy) = coss(ix-1,iy) + 0.5_R8 * hx(ix-1,iy)/
     1         (sqrt(bb(ix-1,iy,0)**2+bb(ix-1,iy,1)**2)/bb(ix-1,iy,3))
     2         + 0.5_R8 * hx(ix,iy)/
     3         (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
            end do
            seglen = coss(rightcut(1)-1,iy) +
     1       0.5_R8 * hx(rightcut(1)-1,iy)/
     1       (sqrt(bb(rightcut(1)-1,iy,0)**2+bb(rightcut(1)-1,iy,1)**2)/
     2       bb(rightcut(1)-1,iy,3))
            coss(rightcut(2):rightcut(1)-1,iy) =
     .       -sin(pi*coss(rightcut(2):rightcut(1)-1,iy)/seglen)
         end do
      else if (nnreg(0).eq.4) then
cxpb  Compute cos(s) factor
         do iy = -1, jsep
            coss(leftcut(1),iy) = 0.5_R8 * hx(leftcut(1),iy)/
     1         (sqrt(bb(leftcut(1),iy,0)**2+bb(leftcut(1),iy,1)**2)/
     2          bb(leftcut(1),iy,3))
            do ix=leftcut(1)+1,rightcut(1)-1
              coss(ix,iy) = coss(ix-1,iy) + 0.5_R8 * hx(ix-1,iy)/
     1         (sqrt(bb(ix-1,iy,0)**2+bb(ix-1,iy,1)**2)/bb(ix-1,iy,3))
     2         + 0.5_R8 * hx(ix,iy)/
     3         (sqrt(bb(ix,iy,0)**2+bb(ix,iy,1)**2)/bb(ix,iy,3))
            end do
            seglen = coss(rightcut(1)-1,iy) +
     1       0.5_R8 * hx(rightcut(1)-1,iy)/
     1       (sqrt(bb(rightcut(1)-1,iy,0)**2+bb(rightcut(1)-1,iy,1)**2)/
     2       bb(rightcut(1)-1,iy,3))
            coss(leftcut(1):rightcut(1)-1,iy) =
     .       sin(2.0_R8*pi*coss(leftcut(1):rightcut(1)-1,iy)/seglen)
         end do
      else
         call xerrab
     .     ('Invalid number of regions : Case not yet coded!')
      end if
*   ..computation
      wb = 2.0_R8 * L_perp / R
      do is = 0, ns-1
        if (.not.is_neutral(is)) then
          do ix = -1, nx
            do iy = -1, jsep
              vla0(ix,iy,1,is) = -alpha*wb*coss(ix,iy)*dna0(ix,iy,is)/
     .         (sqrt(am(is)*mp*te(ix,iy))/qe/bb(ix,iy,3))  ! conversion factor rho_s
            enddo
            do iy = jsep+1, ny
              vla0(ix,iy,1,is) = 0.0_R8
              dna0(ix,iy,is) = max(switch%b2tqna_min_df0,
     .         min(switch%b2tqna_max_df0,
     .         dna0(ix,jsep,is)+vla0(ix,jsep,1,is)*hy1(ix,jsep)))
            enddo
          enddo
        endif
      enddo
*   ..return
      ncall = ncall+1
      call subend ()
      return
      end subroutine set_transport_xpb

c mw Added routine set_transport_mw
c idea: include radial velocity profile for ions between the top
c of the machine and the X-point along the low field side
c under development

      subroutine set_transport_mw(nx, ny, ns, ne, jsep, jxa,
     & dna0, vla0, hci0, hcib, hce0)
      use b2mod_types
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, jsep, jxa
*   ..inout arguments
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1),
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1),
     &   hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1),
     &   ne(-1:nx,-1:ny)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny)
*   ..common blocks
      integer ncall, ix, iy, is, xtop, xbottom
      integer mw_interpol, mw_jxa, mw_jsep                           !mw
      real (kind=R8) ::                                              !mw
     *     vsep, vwall,                                              !mw
c     *     vtop, vox,
     &     mw_flux, mw_alpha, mw_flux_wall, vold, mw_fluxalpha       !mw
       save ncall,
     *     vsep, vwall,                                              !mw
c     *     vtop, vox,
     &     mw_interpol, mw_flux, mw_alpha, mw_flux_wall              !mw
      data ncall/0/,
     * vsep/0.0_R8/, vwall/0.0_R8/,                                  !mw
c     * vtop/0.0_R8/, vox/0.0_R8/,                                    !mw
     & mw_flux/0.0_R8/, mw_alpha/1.0_R8/,mw_fluxalpha/0.0_R8/,       !mw
     & mw_interpol/0/, mw_flux_wall/0.0_R8/,                         !mw
     & vold/0.0_R8/                                                  !mw
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('set_transport_mw')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call ipgeti ('b2tqna_mw_ix1', xtop)
        call ipgeti ('b2tqna_mw_ix2', xbottom)
        call ipgetr ('b2tqna_mw_vsep', vsep)
        call ipgetr ('b2tqna_mw_vwall', vwall)
        call ipgetr ('b2tqna_mw_flux', mw_flux)
        mw_flux_wall = mw_flux
        call ipgetr ('b2tqna_mw_flux_wall', mw_flux_wall)
        call ipgetr ('b2tqna_mw_alpha', mw_alpha)
c        call ipgetr ('b2tqna_mw_vtop', vtop)
c        call ipgetr ('b2tqna_mw_voxpoint', vox)
        call ipgeti ('b2tqna_mw_interpol', mw_interpol)
        mw_jxa = jxa !pre set value for mw_jxa = standard ix for midplane
        call ipgeti ('b2tqna_mw_jxa', mw_jxa)
        mw_jsep = jsep
        call ipgeti ('b2tqna_mw_jsep', mw_jsep)
        call ipgetr ('b2tqna_mw_fluxalpha', mw_fluxalpha)
      endif

*--------setting diffusive and convective transport coefficients

      do iy=-1, ny
        do ix=-1, nx
c         inside separatrix
          if(ix.ge.xtop.and.ix.le.xbottom) then
c         in the transition zone
          if(iy.ge.-1.and.iy.le.mw_jsep) then
            dna0(ix,iy,:)=dna0(ix,iy,:)
            hce0(ix,iy)=hce0(ix,iy)
            hcib(ix,iy,:)=hcib(ix,iy,:)
            vla0(ix,iy,0,:) = vla0(ix,iy,0,:)
            vla0(ix,iy,1,1) = vsep
            if (mw_interpol.ge.3) then
              vla0(ix,iy,1,1) = mw_flux * mw_alpha / ne(ix,iy)
            endif
          endif
c         in the SOL
          if(iy.gt.mw_jsep.and.iy.le.ny) then
            dna0(ix,iy,:)=dna0(ix,iy,:)
            hce0(ix,iy)=hce0(ix,iy)
            hcib(ix,iy,:)=hcib(ix,iy,:)
            vla0(ix,iy,0,:) = vla0(ix,iy,0,:)
            if(mw_interpol.eq.0) then !default value
              vla0(ix,iy,1,1) = (vwall-vsep)/(ny-(mw_jsep+1))*
     &         (iy-(mw_jsep+1))+vsep
            elseif(mw_interpol.eq.1.and.ix.le.mw_jxa) then !sin poloidal
              vla0(ix,iy,1,1)=(((vwall-vsep)/(ny-(mw_jsep+1))*
     &         (iy-(mw_jsep+1)))+vsep)*
     &         sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
            elseif(mw_interpol.eq.1.and.mw_jxa.lt.ix) then !sin poloidal
              vla0(ix,iy,1,1)=(((vwall-vsep)/(ny-(mw_jsep+1))*
     &         (iy-(mw_jsep+1)))+vsep)*sin(pi/2.0_R8*(ix-xbottom)/
     &         (mw_jxa-xbottom))
            elseif(mw_interpol.eq.2) then !sin radial
              vla0(ix,iy,1,1) = (vwall-vsep)*sin(pi/2.0_R8*(iy-mw_jsep)/
     &         (ny-mw_jsep))+vsep
            elseif(mw_interpol.eq.3) then
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
            elseif(mw_interpol.eq.4.and.ix.le.mw_jxa) then !sin poloidal
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &         sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
            elseif(mw_interpol.eq.4.and.mw_jxa.lt.ix) then !sin poloidal
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &         sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
            elseif(mw_interpol.eq.5.and.ix.le.mw_jxa) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy) < ne(topix(ix,iy),topiy(ix,iy))) then
                  vla0(ix,iy,1,1) =
     &             0.9_R8*vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &             sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
                endif
              else
                if(ne(ix,iy).gt.ne(bottomix(ix,iy),bottomiy(ix,iy)))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)*
     &             sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
                endif
              endif
            elseif(mw_interpol.eq.5.and.mw_jxa.lt.ix) then
              if(topiy(ix,iy).le.ny) then
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &             sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
                endif
              else
                if(ne(ix,iy).gt.ne(bottomix(ix,iy),bottomiy(ix,iy)))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)*
     &             sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
                endif
              endif
            elseif(mw_interpol.eq.6.and.ix.le.mw_jxa) then
              if(topiy(ix,iy).le.ny) then
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
              endif
            elseif(mw_interpol.eq.6.and.mw_jxa.lt.ix) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy).lt.ne(topix(ix,iy),topiy(ix,iy))) then
                  vla0(ix,iy,1,1) =
     &             0.9_R8*vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
              endif
c MW working
            elseif(mw_interpol.eq.7.and.ix.le.mw_jxa) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy).lt.ne(topix(ix,iy),topiy(ix,iy))) then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)*
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))/ne(ix,iy)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)*
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))/ne(ix,iy)
                else
                  vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
                endif
              endif
            elseif(mw_interpol.eq.7.and.mw_jxa.lt.ix) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy).le.ne(topix(ix,iy),topiy(ix,iy))) then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)*
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))/ne(ix,iy)
                else
                  vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &             0.9_R8*vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) =
     &             vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)*
     &             ne(bottomix(ix,iy),bottomiy(ix,iy))/ne(ix,iy)
                else
                  vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
                endif
              endif
c MW working end
c MW working 2
            elseif(mw_interpol.eq.8) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy).le.ne(topix(ix,iy),topiy(ix,iy))) then
                   vold = vla0(ix,iy,1,1)
                   vla0(ix,iy,1,1) = vold - vold * 0.5_R8 * mw_fluxalpha
                 if(vla0(ix,iy,1,1).lt.0.0_R8) then
                   vla0(ix,iy,1,1) = 0.0_R8
                 endif
                else
                  vla0(ix,iy,1,1) = (mw_flux / ne(ix,iy) -
     &            vla0(ix,iy,1,1)) * mw_fluxalpha + vla0(ix,iy,1,1)
                endif
              else
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vold = vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) = vold - vold * 0.5_R8 * mw_fluxalpha
                  if(vla0(ix,iy,1,1).lt.0.0_R8) then
                    vla0(ix,iy,1,1) = 0.0_R8
                  endif
                else
                  vla0(ix,iy,1,1) = (mw_flux_wall / ne(ix,iy) -
     &            vla0(ix,iy,1,1)) * mw_fluxalpha + vla0(ix,iy,1,1)
                endif
              endif
c MW working end 2
c MW working 3
c MW same as int8 + mw_flux also for cells at mw_jsep - transport
c through separatrix
C This is the one really used...
            elseif(mw_interpol.eq.9) then
              if(topiy(ix,iy).le.ny) then
                if(ne(ix,iy).le.ne(topix(ix,iy),topiy(ix,iy))) then
                   vold = vla0(ix,iy,1,1)
                   vla0(ix,iy,1,1) = vold - vold * 0.5_R8 * mw_fluxalpha
                  if(vla0(ix,iy,1,1).lt.0.0_R8) then
                    vla0(ix,iy,1,1) = 0.0_R8
                  endif
                else
                  vla0(ix,iy,1,1) = (mw_flux / ne(ix,iy) -
     &            vla0(ix,iy,1,1)) * mw_fluxalpha + vla0(ix,iy,1,1)
                endif
                if(iy.eq.mw_jsep+1) then
                  if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     &             then
                    vold = vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                    vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &               vold - vold * 0.5_R8 * mw_fluxalpha
                    if(vla0(bottomix(ix,iy),
     &                      bottomiy(ix,iy),1,1).lt.0.0_R8) then
                      vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) = 0.0_R8
                    endif
                  else
                    vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1) =
     &              (mw_flux / ne(bottomix(ix,iy),bottomiy(ix,iy)) -
     &               vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)) *
     &               mw_fluxalpha +
     &               vla0(bottomix(ix,iy),bottomiy(ix,iy),1,1)
                  endif
                endif
              else
                if(ne(bottomix(ix,iy),bottomiy(ix,iy)).le.ne(ix,iy))
     >           then
                  vold = vla0(ix,iy,1,1)
                  vla0(ix,iy,1,1) = vold - vold * mw_fluxalpha
                  if(vla0(ix,iy,1,1).lt.0.0_R8) then
                    vla0(ix,iy,1,1) = 0.0_R8
                  endif
                else
                  vla0(ix,iy,1,1) = (mw_flux_wall / ne(ix,iy) -
     &            vla0(ix,iy,1,1)) * mw_fluxalpha + vla0(ix,iy,1,1)
                endif
              endif
c MW working end 3
            else
              write(*,*) 'mw_interpol = ',mw_interpol,' not coded'
              call xerrab ('mw_interpol value not coded')
            endif
          endif
          else
            dna0(ix,iy,:)=dna0(ix,iy,:)
            hce0(ix,iy)=hce0(ix,iy)
            hcib(ix,iy,:)=hcib(ix,iy,:)
            vla0(ix,iy,:,:) = vla0(ix,iy,:,:)
          endif
        enddo !end ix loop
      enddo !end iy loop
      do ix=-1,nx
        do iy=-1,ny
          hci0(ix,iy)=0.0_R8
          do is=0,ns-1
            hci0(ix,iy)=hci0(ix,iy)+hcib(ix,iy,is)
          enddo
        enddo
      enddo

* ..return
      ncall = ncall+1
      call subend ()
      return

      end subroutine set_transport_mw

      subroutine set_transport_blob(nx, ny, ns, jsep, jxa, vla0)
      use b2mod_types
      use b2mod_geo
      use b2mod_b2cmpa
      use b2mod_indirect
      use b2mod_transport_namelist
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, jsep, jxa
*   ..output arguments
      real (kind=R8) ::
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1)
*   ..local variables
      integer ncall
      integer ix, iy, is, ibegin, iend
      real (kind=R8) :: dx(-1:nx), dy(-1:ny), p_p(-1:nx), p_r(-1:ny)
      real (kind=R8) :: dbegin, dend, fac
*   ..procedures
      real (kind=R8) :: hy1
      intrinsic sqrt, max
      external subini, subend, xerrab
*   ..initialisation
      save ncall, ibegin, iend
      data ncall /0/
*   ..statement functions
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)

c*** "Blobby" convection following S. Krasheninnikov
c*** Here we assume the power radial profile of the outward velocity,
c*** starting from zero at the separatrix and reaching the specified
c*** value at the grid edge. The poloidal profile is of a bell shape
c*** determined by two powers: (1-x**n1)**n2.

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('set_transport_blob')
      if (ncall.eq.0) then
        if (nnreg(0).eq.2.or.nnreg(0).eq.5) then
          ibegin=-1
          iend=nx
        else if (nnreg(0).eq.4) then
          ibegin=leftcut(1)
          iend=rightcut(1)-1
        else if (nnreg(0).eq.8) then
          ibegin=leftcut(2)
          iend=rightcut(1)-1
        else
          call xerrab('Unexpected geometry for blobby transport!')
        endif
      endif

c*** We compute the radial distance away from the separatrix
      dy(jsep+1)=0.5_R8*hy1(jxa,jsep+1)
      do iy=jsep+2,ny
        dy(iy)=dy(iy-1)+0.5_R8*(hy1(jxa,iy-1)+hy1(jxa,iy))
      enddo
c*** Calculate the radial profile as a function of the cell index
      do iy=jsep+1,ny
        p_r(iy)=dy(iy)/dy(ny)
      enddo

c*** We compute the poloidal distance away from the midplane
      dx(jxa)=0.0_R8
      do ix=jxa-1,ibegin,-1
        dx(ix)=dx(ix+1)+0.5_R8*(hx(ix+1,jsep+1)+hx(ix,jsep+1))
      enddo
      dbegin=dx(ibegin)+0.5_R8*hx(ibegin,jsep+1)
      do ix=jxa+1,iend
        dx(ix)=dx(ix-1)+0.5_R8*(hx(ix-1,jsep+1)+hx(ix,jsep+1))
      enddo
      dend=dx(iend)+0.5_R8*hx(iend,jsep+1)
c*** Calculate the poloidal profile as a function of the cell index
      do ix=ibegin,jxa-1
        p_p(ix)=dx(ix)/dbegin
      enddo
      p_p(jxa)=0.0_R8
      do ix=jxa+1,iend
        p_p(ix)=dx(ix)/dend
      enddo

c*** Combine the profiles
c*** We only apply this in the outer SOL region
      do is=0, ns-1
        if (.not.is_neutral(is)) then
          do ix=ibegin, iend
            fac=max(0.0_R8,(1.0_R8-p_p(ix)**pw1_cnv(is))**pw2_cnv(is))*
     .       vout_cnv(is)
            do iy=jsep+1, ny
              vla0(ix,iy,1,is)=fac*p_r(iy)**pw0_cnv(is)
            end do
          end do
        end if
      end do

* ..return
      ncall = ncall+1
      call subend ()
      return

      end subroutine set_transport_blob

      subroutine set_transport_disruption(nx, ny, ns,
     & dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0,
     & na, ne, te, ti, crx, hy1)
      use b2mod_types
      use b2mod_rates
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
      use b2mod_b2cmpt
      use b2mod_transport_disruption
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   te(-1:nx,-1:ny), crx(-1:nx,-1:ny,0:3),
     *   ti(-1:nx,-1:ny), hy1(-1:nx,-1:ny)
*   ..in/out arguments
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &   hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) :: dfhce, dfhcb
      real (kind=R8) :: ne(-1:nx,-1:ny), na(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..local variables
      integer ncall, ix, iy, is, ixref, iyref
      integer get_quenching
      real (kind=R8):: get_time_quenching
      save ncall, ixref, iyref
*   ..initialisation
      data ncall/0/

*   ..procedures
      external ipgeti, ipgetr, xertst
* ..preliminaries
*   ..subprogram start-up calls
      call subini ('set_transport_disruption')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call alloc_b2mod_transport_disruption(nx,ny)
        if(nnreg(0).gt.1 .and. nnreg(0).le.5) then
          ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
          iyref=2*topcut(1)/3
        elseif(nnreg(0).eq.8) then
          ixref=(rightcut(1)+rightcut(2))/2
          iyref=2*min(topcut(1),topcut(2))/3
        else
          ixref=3*nx/4
          iyref=ny/2
        endif
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')

        quenching = 0
        time_quenching = 0._R8
        call ipgeti ('b2tqna_quenching', get_quenching)
        call ipgetr ('b2tqna_time_quenching', get_time_quenching)
        quenching = get_quenching
        time_quenching = get_time_quenching

        write(*,*) 'Quenching: ', quenching

*   ..calculate profile f_quench(ix,iy) for chi_e, chi_i, and Dn_i   !cjk
*         f_quench=1
*         i.e. constant profile
        if (quenching .eq. 1) then
          do iy=-1,ny
*         ..flux surface label
            rho_(iy) = (crx(ixref,iy,0)+crx(ixref,iy,2))/2._R8
            do ix=-1,nx
              f_quench(ix,iy)=1._R8
            enddo
          enddo
        endif

      endif

*   ..strongly enhanced perpendicular heat and particle transport for  !cjk
*       plasma disruption (quenching=1, time > time_quenching)
*
*       chi_e=1000*f_quench*ne_ref/ne*Te^5/2*geometry_factor
*       chi_i=1000/sqrt(mi/me)*f_quench*ni_ref/ni*Ti^5/2*geometry_factor
*       Dn_i=1000/sqrt(mi/me)*f_quench*Ti^5/2*geometry_factor
*
*       i.e. kappa_perp,e and kappa_perp,i independent of ne and ni

      if (quenching .eq. 1) then
         if (tim .ge. time_quenching) then
           do ix=-1,nx
             do iy=-1,ny
               dfhce=1.e3_R8*f_quench(ix,iy)*(ne(ixref,iyref)/
     *               ne(ix,iy))*(te(ix,iy)/(1.e3_R8*ev))**2.5_R8*
     *               (hy1(ix,iy)/hy1(ixref,iy))**2
               do is=0, ns-1
*          charged ions only
                 if (.not. is_neutral(is)) then
                   dfhcb=1.e3_R8*f_quench(ix,iy)/sqrt(am(is)*mp/me)*
     *                   (na(ixref,iyref,is)/na(ix,iy,is))*
     *                   (ti(ix,iy)/(1.e3_R8*ev))**2.5_R8*
     *                   (hy1(ix,iy)/hy1(ixref,iy))**2
                   hcib(ix,iy,is)=hcib(ix,iy,is)+dfhcb*na(ix,iy,is)
                   hci0(ix,iy)=hci0(ix,iy)+dfhcb*na(ix,iy,is)
                   dna0(ix,iy,is)=dna0(ix,iy,is)+1.e3_R8*
     *                   f_quench(ix,iy)/sqrt(am(is)*mp/me)*
     *                   (ti(ix,iy)/(1.e3_R8*ev))**2.5_R8*
     *                   (hy1(ix,iy)/hy1(ixref,iy))**2
                 endif
               enddo
               hce0(ix,iy)=hce0(ix,iy)+dfhce*ne(ix,iy)
             enddo
           enddo
         endif
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return

      end subroutine set_transport_disruption

      subroutine set_transport_neo(nx, ny, ns, hx, hy, vol,
     &  crx, cry, conn, qz, bb, na, te, ti, ni, ne,
     &  dna0, dpa0, vla0, vsa0, hci0, hcib, hce0, sig0, alf0,
     &  hvi0, hve0)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_ranges
      use b2mod_rates
      use b2mod_neoclassical
      use b2mod_time
      use b2mod_transpbgr_namelist
      use b2mod_external
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns
      real (kind=R8) ::
     *   bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     &   te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &   ne(-1:nx,-1:ny), hx(-1:nx,-1:ny), hy(-1:nx,-1:ny),
     &   vol(-1:nx,-1:ny), conn(-1:ny), qz(-1:nx,-1:ny,0:1),
     &   crx(-1:nx,-1:ny,0:3), cry(-1:nx,-1:ny,0:3)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &   hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny),
     &   hvi0(-1:nx,-1:ny,0:1), hve0(-1:nx,-1:ny,0:1)

      integer ::
     *  nar,nzm
*     nints=number of intervals for integration
      integer nints
      parameter (nints=10)
      real (kind=R8) ::
     *   bav, b2av, bi2a, rbt, bgradp, dpsidr, rnq, fc, gclass,
     &   fm(0:nmaxgr-1), r2i, cosbl(0:nmaxgr-1),
     &   cosb2(0:nmaxgr-1), sinbl(0:nmaxgr-1), sinb2(0:nmaxgr-1)
      common /b2ngeom/ bav,b2av,bi2a,rbt,bgradp,dpsidr,rnq,fc,
     *                 gclass,fm,r2i,cosbl,cosb2,sinbl,sinb2

*   ..local variables
      integer ncall, n2, ix, iy, is, is2, ixref, iyref,
     *        transpbgr_namelist
      save ncall, ixref, iyref, transpbgr_namelist
      data transpbgr_namelist/1/
*   ..local variables for neoclassical transport
      integer, save :: ic
      integer :: istate,ispecies,ishift
      integer :: isorth(0:nsdecl,0:nzdecl-1), nclh(0:nsdecl), nslh
      real (kind=R8) ::
     *   zah(0:nsdecl-1), amh(0:nsdecl-1), t(0:nsdecl),
     &   den(0:nsdecl,0:nzdecl-1), ds(0:nsdecl,0:nzdecl-1,2),
     &   denl(0:nsdecl,0:nzdecl-1), tl(0:nsdecl),
     &   denr(0:nsdecl,0:nzdecl-1), tr(0:nsdecl),
     &   vtn(0:nsdecl,0:nzdecl-1,-1:ny),
     &   kappa(0:nsdecl,0:nzdecl-1,-1:ny),
     &   vw(0:nsdecl,0:nzdecl-1,-1:ny),
     &   dfhcb, dfhce,
     &   hthet(-1:nx,-1:ny), quant(-1:nx,-1:ny),
     &   cellface(0:nsdecl,0:nzdecl-1,-1:ny,1:5)
      real (kind=R8) ::
     *   zacomp, amcomp, epsn
*   ..set the accuracy
      parameter (epsn=1.e-3)
      real (kind=R8) ::
     *   thet, ilpb, hxface, ilp(-1:ny), bmax, nlambda, nfactor,
     &   hxpface, hxnface, eparr, rho, prefactor, postfactor
      save eparr
      real (kind=R8) ::
     *   fluxav, sfluxav, tfluxav
      logical closedf
*   ..for consistency check
      real (kind=R8) ::
     *   densav(0:nsdecl,0:nzdecl-1,-1:ny),
     &   gradnsav(0:nsdecl,0:nzdecl-1,-1:ny),
     &   tforce(0:nsdecl,0:nzdecl-1,1:2,-1:ny)

*   ..procedures
      real (kind=R8) :: hy1
      external ipgeti, ipgetr, xertst, b2xvsg, fluxav, sfluxav,
     &         tfluxav, neodv, neochi
*   ..statement functions
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)

*   ..set the contribution ic
*     0 --- classical particle flux
*     1 --- banana plateau contribution
*     2 --- Pfirsch Schlueter contribution
*     3 --- both banana and PS
*     4 --- all contributions
*     mind that B2 already calculates the classical transport !
*     avoid double transport, 0+4 for cross checks only !
      data ic /3/

*   ..subprogram start-up calls
      call subini ('set_transport_neo')
      if (ncall.eq.0) then
*   ..neoclassical transport only valid for single charge states
*     bundled charge states yet to be implemented
        if(nnreg(0).gt.1 .and. nnreg(0).le.5) then
          ixref=rightcut(1)-(rightcut(1)-leftcut(1))/4
          iyref=2*topcut(1)/3
        elseif(nnreg(0).eq.8) then
          ixref=(rightcut(1)+rightcut(2))/2
          iyref=2*min(topcut(1),topcut(2))/3
        else
          ixref=3*nx/4
          iyref=ny/2
        endif
        call ipgeti ('set_transport_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument set_transport_ixref')
        call ipgeti ('set_transport_iyref', iyref)
        call xertst (-1.le.iyref.and.iyref.le.ny,
     1   'faulty argument set_transport_iyref')
	call ipgeti ('neoclassical_ic', ic)
        call xertst(0.le.ic .and. ic.le.4,
     .   'invalid parameter neoclassical_ic')
	write(*,*) 'NEOART option ', ic, ' has been used'
	if(ic.eq.0)
     1	 write(*,*) '  NEOART: classical'
	if(ic.eq.1)
     1	 write(*,*) '  NEOART: banana plateau'
	if(ic.eq.2)
     1	 write(*,*) '  NEOART: Pfirsch-Schlueter'
	if(ic.eq.3)
     1	 write(*,*) '  NEOART: banana plateau + Pfirsch-Schlueter'
	if(ic.eq.4)
     1	 write(*,*) '  NEOART: all'
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*    ..test bb
       call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*    ..test state
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
      endif

*   ..main

*   ..initialize densities and temperatures
      do ispecies=0,nsdecl
        t(ispecies)=0.0_R8
        do istate=0,nzdecl-1
          den(ispecies,istate)=0.0_R8
        enddo
      enddo

*   ..species reordering and calculation of geometric parameters
*       on first call only

      if (ncall.eq.0) then

*  ..regions
c      do ix=-1,nx
c        write(*,*) '0: ix = ',ix,' in region ', region(ix,0,0)
c        write(*,*) 'ny: ix = ',ix,' in region ', region(ix,ny,0)
c      enddo


*   ..obtain input
      call ipgeti ('b2tqna_transpbgr_namelist', transpbgr_namelist)
      if(transpbgr_namelist.ge.1) then
        call read_b2mod_transpbgr_namelist
        call write_b2mod_transpbgr_namelist
      endif

*    ..loop voltage eparr (to be included in a separate routine!)
*       default value:
      eparr=0.0_R8
*       read scalar from b2mn.dat
      call ipgetr ('b2tqna_eparr', eparr)
      write(*,*) 'eparr =', eparr

*   ..initialize geometric parameters
      do iy=-1,ny
        bav_(iy)=0.0_R8
        b2av_(iy)=0.0_R8
        bi2a_(iy)=0.0_R8
        rbt_(iy)=0.0_R8
        bgradp_(iy)=0.0_R8
        dpsidr_(iy)=0.0_R8
        rnq_(iy)=0.0_R8
        fc_(iy)=0.0_R8
        gclass_(iy)=0.0_R8
        do is=1,mmx
          fm_(is,iy)=0.0_R8
        enddo
        r2i_(iy)=0.0_R8
        rho_(iy)=0.0_R8
      enddo

*   ..split species into species and charge states
      nar=nsdecl
      nzm=nzdecl-1
      do ispecies=0,nsdecl
        ncl(ispecies)=0
        do istate=0,nzdecl-1
          isort(ispecies,istate)=0
        enddo
      enddo
      do is=0,ns-1
        zah(is)=(zamin(is)+zamax(is))/2.0_R8
        if (zamin(is).ne.zamax(is)) then
          write(6,*) 'Warning: '//
     .     'neoclassical transport model may not be accurate '//
     .     'for bundled charge states!'
        end if
        amh(is)=am(is)
      enddo
      do is=0,ns-2
        zacomp=zah(is)
        amcomp=amh(is)
        ispecies=is
        do is2=is+1,ns-1
          if((abs(amh(is2)-amcomp).lt.0.5_R8).and.
     &       (zah(is2).lt.zacomp)) then
            zacomp=zah(is2)
            ispecies=is2
          elseif(amh(is2).lt.amcomp) then
            amcomp=amh(is2)
            zacomp=zah(is2)
            ispecies=is2
          endif
        enddo
        zah(ispecies)=zah(is)
        amh(ispecies)=amh(is)
        zah(is)=zacomp
        amh(is)=amcomp
      enddo
      ispecies=1
      istate=1
      amcomp=amh(0)
      do is=1,ns-1
        if((amh(is)-amcomp).gt.0.5_R8) then
          ncl(ispecies-1)=istate
          ispecies=ispecies+1
          amcomp=amh(is)
          istate=1
        elseif(abs(amh(is)-amcomp).lt.0.5_R8) then
          istate=istate+1
        endif
      enddo

      ncl(ispecies-1)=istate

      nsl=ispecies

      is=0
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          do is2=0,ns-1
            if((abs(am(is2)-amh(is)).lt.0.5_R8).and.
     &         (abs((zamin(is2)+zamax(is2))/2.0_R8-zah(is)).lt.
     .            0.5_R8)) then
              isort(ispecies,istate)=is2
            endif
          enddo
          is=is+1
        enddo
      enddo

*   ..remove neutrals from list
      do ispecies=0,nsl-1
        nclh(ispecies)=ncl(ispecies)
      enddo
      nslh=nsl
      isorth=isort
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          if (is_neutral(isort(ispecies,istate))) then
            nclh(ispecies)=nclh(ispecies)-1
            if (nclh(ispecies).eq.0) nslh=nslh-1
          endif
        enddo
      enddo
      is=0
      do ispecies=0,nsl-1
        if (nclh(ispecies).ne.0) then
           nclh(is)=nclh(ispecies)
           is=is+1
        endif
      enddo
      is=0
      is2=0
      do ispecies=0,nsl-1
        do istate=0,ncl(ispecies)-1
          if (.not. is_neutral(isorth(ispecies,istate))) then
            isort(is,is2)=isorth(ispecies,istate)
            is2=is2+1
            if (is2.gt.(nclh(is)-1)) then
                is=is+1
                is2=0
              endif
            endif
          enddo
        enddo

        do ispecies=0,nsl-1
          if (ispecies.lt.nslh) then
            ncl(ispecies)=nclh(ispecies)
          else
            ncl(ispecies)=0
          endif
        enddo
        nsl=nslh

*   ..mass
        do ispecies=0,nsl-1
          m(ispecies)=am(isort(ispecies,0))*mp
        enddo

*   ..add electrons as a separate species
        nsl=nsl+1
        ncl(nsl-1)=1
        m(nsl-1)=me
        zsp(nsl-1,0)=-1.0_R8

*   ..calculate generalized poloidal angle theta
*      default value 0 for grid points outside the core and on
*      the inner boundary
        do iy=-1,ny
          do ix=-1,nx
            theta(ix,iy)=0._R8
          enddo
        enddo
        do iy=0,ny
          thet=0._R8
          ilp(iy)=0._R8
          do ix=ixref,-1,-1
           if (on_closed_surface(ix,iy)) then
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
          enddo
          do ix=nx,ixref+1,-1
           if (on_closed_surface(ix,iy)) then
             theta(ix,iy)=thet
             hxface=(vol(ix,iy)*hx(bottomix(ix,iy),bottomiy(ix,iy))+
     &             vol(bottomix(ix,iy),bottomiy(ix,iy))*hx(ix,iy))/
     &            (vol(bottomix(ix,iy),bottomiy(ix,iy))+vol(ix,iy))
             hxpface=(vol(leftix(ix,iy),leftiy(ix,iy))*
     &              hx(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))*
     &              hx(leftix(ix,iy),leftiy(ix,iy)))/
     &              (vol(bottomix(leftix(ix,iy),leftiy(ix,iy)),
     &              bottomiy(leftix(ix,iy),leftiy(ix,iy)))+
     &              vol(leftix(ix,iy),leftiy(ix,iy)))
             hxnface=(vol(rightix(ix,iy),rightiy(ix,iy))*
     &              hx(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))*
     &              hx(rightix(ix,iy),rightiy(ix,iy)))/
     &              (vol(bottomix(rightix(ix,iy),rightiy(ix,iy)),
     &              bottomiy(rightix(ix,iy),rightiy(ix,iy)))+
     &              vol(rightix(ix,iy),rightiy(ix,iy)))
            thet=thet+0.5_R8*bb(ix,iy,3)/bb(ix,iy,0)*
     &           ((hxpface+hxnface)/2._R8+hxface)
            ilp(iy)=ilp(iy)+0.5_R8*((hxpface+hxnface)/2._R8+hxface)/
     &              bb(ix,iy,0)
           endif
          enddo
          theta(ixref,iy)=thet
          ilpb=thet
          do ix=-1,nx
           if (on_closed_surface(ix,iy)) then
            theta(ix,iy)=2._R8*pi*theta(ix,iy)/ilpb
           endif
          enddo
        enddo

*   ..calculate geometric parameters

*   ..loop over flux surfaces
        do iy=0,ny

*   ..calculate neoclassical transport only for closed flux surfaces
         closedf=.false.
         do ix=-1,nx
          if (on_closed_surface(ix,iy)) then
            closedf=.true.
          endif
         enddo
         if (closedf) then
           write(*,*) 'closed flux surface', iy
           closedf=.false.

*   ..flux surface averages

           quant(:,:)=bb(:,:,3)
           bav_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)
           quant(:,:)=bb(:,:,3)**2
           b2av_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)
           quant(:,:)=1.0_R8/bb(:,:,3)**2
           bi2a_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)
           quant(:,:)=bb(:,:,2)*(crx(:,:,0)+crx(:,:,1)+
     &                           crx(:,:,2)+crx(:,:,3))/4.0_R8
           rbt_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)
           bgradp_(iy)=2._R8*pi/(bav_(iy)*ilp(iy))
           dpsidr_(iy)=bb(ixref,iy,0)*(crx(ixref,iy,0)+
     &            crx(ixref,iy,2))/2._R8
           rnq_(iy)=conn(iy)/(2._R8*pi)

           bmax=0._R8
           do ix=-1,nx
             if (on_closed_surface(ix,iy)) then
               if (bb(ix,iy,3).gt.bmax) then
                 bmax=bb(ix,iy,3)
               endif
             endif
           enddo
           nlambda=1.0_R8/nints

           fc_(iy)=0
           do is2=0,nints
             if ((is2.eq.0).or.(is2.eq.nints)) then
               nfactor=0.5_R8
             else
               nfactor=1.0_R8
             endif
             quant(:,:)=bb(:,:,3)
             fc_(iy)=fc_(iy)+is2*nlambda/sfluxav(quant,bb,hx,vol,
     &            nx,ny,iy,is2*nlambda,bmax)*nlambda*nfactor
           enddo
           fc_(iy)=0.75_R8*b2av_(iy)/bmax**2*fc_(iy)

           quant=((crx(:,:,0)+crx(:,:,1)+crx(:,:,2)+crx(:,:,3))
     &             /4.0_R8*bb(:,:,0)/bb(:,:,3))**2
           gclass_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)

           do is2=1,mmx
             do is=-1,ny
               do ix=-1,nx
                 hthet(ix,is)=is2*theta(ix,iy)
               enddo
             enddo
             quant(:,:)=cos(hthet(:,:))*bb(:,:,3)*log(bb(:,:,3))
             cosbl(is2)=fluxav(quant,bb,hx,vol,nx,ny,iy)
             quant(:,:)=cos(hthet(:,:))*bb(:,:,3)**2
             cosb2(is2)=fluxav(quant,bb,hx,vol,nx,ny,iy)
             quant(:,:)=sin(hthet(:,:))*bb(:,:,3)*log(bb(:,:,3))
             sinbl(is2)=fluxav(quant,bb,hx,vol,nx,ny,iy)
             quant(:,:)=sin(hthet(:,:))*bb(:,:,3)**2
             sinb2(is2)=fluxav(quant,bb,hx,vol,nx,ny,iy)
           enddo
           do is2=1,mmx
             fm_(is2,iy)=2.0_R8*(is2*bgradp_(iy))**2/
     &               (b2av_(iy)*bav_(iy))*
     &               (cosb2(is2)*cosbl(is2)+sinb2(is2)*sinbl(is2))
           enddo
           quant(:,:)=(4.0_R8/(crx(:,:,0)+crx(:,:,1)+crx(:,:,2)+
     &         crx(:,:,3)))**2
           r2i_(iy)=fluxav(quant,bb,hx,vol,nx,ny,iy)

*    ..flux surface label
           rho_(iy)=(crx(ixref,iy,0)+crx(ixref,iy,2))/2.0_R8

         else
           write(*,*) 'open flux surface ',iy
         endif

        enddo
       else
        write(*,*) 'call nr. ',ncall
      endif


*   ..now calculate the neoclassical transport (if newneo set)

      if (newneo) then

*   ..loop over flux surfaces
       do iy=0,ny

*   ..calculate neoclassical transport only for closed flux surfaces
        closedf=.false.
        do ix=-1,nx
         if (on_closed_surface(ix,iy)) then
           closedf=.true.
         endif
        enddo
        if (closedf) then
          write(*,*) 'closed flux surface', iy
          closedf=.false.

*   ..flux surfaces averages
*     temperatures in keV, densities in 10^19 m^-3
          do ispecies=0,nsl-2
            t(ispecies)=fluxav(ti,bb,hx,vol,nx,ny,iy)/
     &                  (1.0e3_R8*ev)
            tl(ispecies)=fluxav(ti,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &                  (1.0e3_R8*ev)
            tr(ispecies)=tfluxav(ti,bb,hx,vol,nx,ny,iy)/
     &                  (1.0e3_R8*ev)
            do istate=0,ncl(ispecies)-1
              quant(:,:)=na(:,:,isort(ispecies,istate))
              den(ispecies,istate)=
     &          fluxav(quant,bb,hx,vol,nx,ny,iy)/1.0e19_R8
              denl(ispecies,istate)=
     &          fluxav(quant,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &          1.0e19_R8
              denr(ispecies,istate)=
     &          tfluxav(quant,bb,hx,vol,nx,ny,iy)/1.0e19_R8
              quant(:,:)=rza(:,:,isort(ispecies,istate))
              zsp(ispecies,istate)=fluxav(quant,bb,hx,vol,nx,ny,iy)
            enddo
          enddo
          t(nsl-1)=fluxav(te,bb,hx,vol,nx,ny,iy)/(1.0e3_R8*ev)
          den(nsl-1,0)=fluxav(ne,bb,hx,vol,nx,ny,iy)/1.0e19_R8
          tl(nsl-1)=fluxav(te,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &               (1.0e3_R8*ev)
          denl(nsl-1,0)=fluxav(ne,bb,hx,vol,nx,ny,bottomiy(ixref,iy))/
     &                1.0e19_R8
          tr(nsl-1)=tfluxav(te,bb,hx,vol,nx,ny,iy)/(1.0e3_R8*ev)
          denr(nsl-1,0)=tfluxav(ne,bb,hx,vol,nx,ny,iy)/1.0e19_R8

*   ..thermodynamic forces
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              ds(ispecies,istate,1)=(log(denr(ispecies,istate)*
     &           tr(ispecies))-log(denl(ispecies,istate)*
     &           tl(ispecies)))/(hy1(bottomix(ixref,iy),
     &           bottomiy(ixref,iy))+hy1(ixref,iy))
              ds(ispecies,istate,2)=(log(tr(ispecies))-
     &           log(tl(ispecies)))/(hy1(bottomix(ixref,iy),
     &           bottomiy(ixref,iy))+hy1(ixref,iy))


*   ..for consistency check
              densav(ispecies,istate,iy)=den(ispecies,istate)
              gradnsav(ispecies,istate,iy)=den(ispecies,istate)*
     &            (ds(ispecies,istate,1)-ds(ispecies,istate,2))
              tforce(ispecies,istate,1,iy)=ds(ispecies,istate,1)
              tforce(ispecies,istate,2,iy)=ds(ispecies,istate,2)
            enddo
          enddo

*   ..recover geometric parameters for the flux surface iy
          bav=bav_(iy)
          b2av=b2av_(iy)
          bi2a=bi2a_(iy)
          rbt=rbt_(iy)
          bgradp=bgradp_(iy)
          dpsidr=dpsidr_(iy)
          rnq=rnq_(iy)
          fc=fc_(iy)
          gclass=gclass_(iy)
          do is2=1, mmx
            fm(is2)=fm_(is2,iy)
          enddo
          r2i=r2i_(iy)
          rho=rho_(iy)

*    ..calculate the diffusion coefficients
          call neodv(nsdecl,nzdecl,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &               epsn,dd(:,:,iy),vvn(:,:,iy),vtn(:,:,iy))

*    ..calculate the heat conduction coefficients
          call neochi(nsdecl,nzdecl,nsl,ncl,zsp,ic,m,t,den,ds,rho,eparr,
     &                epsn,kappa(:,:,iy),chi(:,:,iy),vw(:,:,iy))

*    ..calculate thermal conductivity and heat pinch from
*        q(i,j) = - kappa(i,j)*grad.T(i) + vw(i,j)*T(i)
*               = - kappa(i,j)*grad.T(i) + den(i,j)*vwt(i,j)*T(i)
          do ispecies=0, nsl-1
            do istate=0, ncl(ispecies)-1
              vwt(ispecies,istate,iy)=vw(ispecies,istate,iy)/
     &                                (den(ispecies,istate)*1e19)
           enddo
          enddo

         else
          write(*,*) 'open flux surface ',iy
          do ispecies=0,nsl-1
            do istate=0,ncl(ispecies)-1
              dd(ispecies,istate,iy)=0.0_R8
              vvn(ispecies,istate,iy)=0.0_R8
              chi(ispecies,istate,iy)=0.0_R8
              kappa(ispecies,istate,iy)=0.0_R8
              vwt(ispecies,istate,iy)=0.0_R8
              densav(ispecies,istate,iy)=0.0_R8
              gradnsav(ispecies,istate,iy)=0.0_R8
              tforce(ispecies,istate,1,iy)=0.0_R8
              tforce(ispecies,istate,2,iy)=0.0_R8
            enddo
          enddo
        endif

       enddo

*    ..set transport coefficients on boundary by symmetric condition
       do ispecies=0,nsl-1
         do istate=0,ncl(ispecies)-1
           densav(ispecies,istate,-1)=densav(ispecies,istate,0)
           gradnsav(ispecies,istate,-1)=gradnsav(ispecies,istate,0)
           tforce(ispecies,istate,1,-1)=tforce(ispecies,istate,1,0)
           tforce(ispecies,istate,2,-1)=tforce(ispecies,istate,2,0)
           dd(ispecies,istate,-1)=dd(ispecies,istate,0)
           vvn(ispecies,istate,-1)=vvn(ispecies,istate,0)
           chi(ispecies,istate,-1)=chi(ispecies,istate,0)
           kappa(ispecies,istate,-1)=kappa(ispecies,istate,0)
           vwt(ispecies,istate,-1)=vwt(ispecies,istate,0)
         enddo
       enddo

c      open(14,file='transp.dat',form='unformatted')
cc       write(*,*) 'transport output nr.', ncall
cc       write(14) tim
cc       write(14) ny,nsl
cc       write(14) (ncl(is),is=0,nsl-1)
cc       do ispecies=0,nsl-1
cc         do istate=0,ncl(ispecies)-1
cc           write(14) (dd(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (vvn(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (chi(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (vwt(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (densav(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (gradnsav(ispecies,istate,iy),iy=-1,ny)
cc           write(14) (tforce(ispecies,istate,1,iy),iy=-1,ny)
cc           write(14) (tforce(ispecies,istate,2,iy),iy=-1,ny)
cc         enddo
cc       enddo
c      close(14)

*    ..calculate the transport coefficients on cell centers
       do iy=-1,ny
         do ispecies=0,nsl-1
           do istate=0,ncl(ispecies)-1
             cellface(ispecies,istate,iy,1)=dd(ispecies,istate,iy)
             cellface(ispecies,istate,iy,2)=vvn(ispecies,istate,iy)
             cellface(ispecies,istate,iy,3)=chi(ispecies,istate,iy)
             cellface(ispecies,istate,iy,4)=kappa(ispecies,istate,iy)
             cellface(ispecies,istate,iy,5)=vwt(ispecies,istate,iy)
           enddo
         enddo
       enddo

       do iy=0,ny
         if (on_closed_surface(ixref,iy)) then
           if (iy.eq.ny) then
             prefactor=-hy1(ixref,iy)/(2.0_R8*
     &                  hy1(bottomix(ixref,iy),bottomiy(ixref,iy)))
             postfactor=1.0_R8+hy1(ixref,iy)/(2.0_R8*
     &                  hy1(bottomix(ixref,iy),bottomiy(ixref,iy)))
             ishift=bottomiy(ixref,iy)
            else
             if (on_closed_surface(topix(ixref,iy),topiy(ixref,iy)))
     >        then
               prefactor=0.5_R8
               postfactor=0.5_R8
               ishift=iy
              else
               prefactor=-hy1(ixref,iy)/(2.0_R8*
     &                    hy1(bottomix(ixref,iy),bottomiy(ixref,iy)))
               postfactor=1.0_R8+hy1(ixref,iy)/(2.0_R8*
     &                    hy1(bottomix(ixref,iy),bottomiy(ixref,iy)))
               ishift=bottomiy(ixref,iy)
             endif
           endif
           do ispecies=0, nsl-1
             do istate=0, ncl(ispecies)-1
               dd(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,1)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),1)*postfactor
               vvn(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,2)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),2)*postfactor
               chi(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,3)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),3)*postfactor
               kappa(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,4)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),4)*postfactor
               vwt(ispecies,istate,iy)=cellface(ispecies,istate,
     &           ishift,5)*prefactor+cellface(ispecies,istate,
     &           topiy(ixref,ishift),5)*postfactor
             enddo
           enddo
         endif
       enddo

*     ..calculate transport coefficients at inner boundary by
*       a mixture of extrapolation and symmetry
       prefactor=0.5_R8*
     1  (hy1(ixref,-1)+hy1(topix(ixref,-1),topiy(ixref,-1)))/
     &  (hy1(topix(ixref,-1),topiy(ixref,-1))+
     &   hy1(topix(topix(ixref,-1),topiy(ixref,-1)),
     &       topiy(topix(ixref,-1),topiy(ixref,-1))))
       do ispecies=0, nsl-1
         do istate=0, ncl(ispecies)-1
           if (on_closed_surface(ixref,-1)) then
             dd(ispecies,istate,-1)=dd(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(dd(ispecies,istate,
     &          topiy(ixref,-1))-dd(ispecies,istate,
     &          topiy(topix(ixref,-1),topiy(ixref,-1))))
             vvn(ispecies,istate,-1)=vvn(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(vvn(ispecies,istate,
     &          topiy(ixref,-1))-vvn(ispecies,istate,
     &          topiy(topix(ixref,-1),topiy(ixref,-1))))
             chi(ispecies,istate,-1)=chi(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(chi(ispecies,istate,
     &          topiy(ixref,-1))-chi(ispecies,istate,
     &          topiy(topix(ixref,-1),topiy(ixref,-1))))
             kappa(ispecies,istate,-1)=kappa(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(kappa(ispecies,istate,
     &          topiy(ixref,-1))-kappa(ispecies,istate,
     &          topiy(topix(ixref,-1),topiy(ixref,-1))))
             vwt(ispecies,istate,-1)=vwt(ispecies,istate,
     &          topiy(ixref,-1))+prefactor*(vwt(ispecies,istate,
     &          topiy(ixref,-1))-vwt(ispecies,istate,
     &          topiy(topix(ixref,-1),topiy(ixref,-1))))
            else
             dd(ispecies,istate,-1)=0.0_R8
             vvn(ispecies,istate,-1)=0.0_R8
             chi(ispecies,istate,-1)=0.0_R8
             kappa(ispecies,istate,-1)=0.0_R8
             vwt(ispecies,istate,-1)=0.0_R8
           endif
         enddo
       enddo

c      open(14,file='transp.dat',form='unformatted')
c       write(*,*) 'transport output nr.', ncall
c       write(14) tim
c       write(14) ny,nsl
c       write(14) (ncl(is),is=0,nsl-1)
c       do ispecies=0,nsl-1
c         do istate=0,ncl(ispecies)-1
c           write(14) (dd(ispecies,istate,iy),iy=-1,ny)
c           write(14) (vvn(ispecies,istate,iy),iy=-1,ny)
c           write(14) (chi(ispecies,istate,iy),iy=-1,ny)
c           write(14) (vwt(ispecies,istate,iy),iy=-1,ny)
c           write(14) (densav(ispecies,istate,iy),iy=-1,ny)
c           write(14) (gradnsav(ispecies,istate,iy),iy=-1,ny)
c           write(14) (tforce(ispecies,istate,1,iy),iy=-1,ny)
c           write(14) (tforce(ispecies,istate,2,iy),iy=-1,ny)
c         enddo
c       enddo
c      close(14)

*    ..set newneo to false (no further calculation of neoclassical
*          transport); will be reset in b2mndt to .true.
       newneo=.false.
       write(*,*) 'recalculated neoclassical transport coefficients'
      else
       write(*,*) 'use previous neoclassical transport coefficients'
      endif


 !       if (mod(ncall,200).eq.0) then
 !         write(*,*) 'temperature output', ncall
 ! c       open(15,file='temp.dat',form='unformatted')
 !         write(15) nx,ny
 !         write(15) nsl
 !         write(15) (ncl(is),is=0,nsl-1)
 !         if (ncall.lt.1) then
 !           write(15) ((0.25_R8*(crx(ix,iy,0)+crx(ix,iy,1)+crx(ix,iy,2)
 !      &               +crx(ix,iy,3)),ix=-1,nx),iy=-1,ny)
 !           write(15) ((0.25_R8*(cry(ix,iy,0)+cry(ix,iy,1)+cry(ix,iy,2)
 !      &               +cry(ix,iy,3)),ix=-1,nx),iy=-1,ny)
 !           write(15) ((bb(ix,iy,0),ix=-1,nx),iy=-1,ny)
 !           write(15) ((bb(ix,iy,1),ix=-1,nx),iy=-1,ny)
 !           write(15) ((bb(ix,iy,2),ix=-1,nx),iy=-1,ny)
 !           write(15) ((bb(ix,iy,3),ix=-1,nx),iy=-1,ny)
 !           write(15) ((hy(ix,iy),ix=-1,nx),iy=-1,ny)
 !           do ispecies=0, nsl-1
 !             write(15) m(ispecies)
 !             if (ispecies.ne.nsl-1) then
 !               do istate=0, ncl(ispecies)-1
 !                 write(15) ((rza(ix,iy,isort(ispecies,istate)),ix=-1,nx),
 !      &                    iy=-1,ny)
 !               enddo
 !             else
 !               write(15) -1.0_R8
 !             endif
 !           enddo
 !         endif
 !         write(15) tim
 !         write(15) ((te(ix,iy)/ev,ix=-1,nx),iy=-1,ny)
 !         write(15) ((ti(ix,iy)/ev,ix=-1,nx),iy=-1,ny)
 !         do ispecies=0, nsl-2
 !           do istate=0, ncl(ispecies)-1
 !             write(15) ((na(ix,iy,isort(ispecies,istate)),ix=-1,nx),
 !      &                 iy=-1,ny)
 !           enddo
 !         enddo
 !         write(15) ((ne(ix,iy),ix=-1,nx),iy=-1,ny)
 !         write(15) (fc_(iy),iy=-1,ny)
 !       endif
 ! c       close(15)


*    ..add the neoclassical contribution to the transport coefficients
*       assuming the following structure for the particle and heat flux
*         gamma(i,j) = - dd(i,j)*grad.n(i,j) + vvn(i,j)*den(i,j)
*           already including the ware pinch, and
*         q(i,j) = - kappa(i,j)*grad.T(i) + den(i,j)*vwt(i,j)*T(i)
*       also including the ware pinch
*       therefore: no contributions to dpa0, vsa0, sig0, and alf0

      do ix=-1,nx
       do iy=-1,ny
        if (on_closed_surface(ix,iy)) then
         do ispecies=0,nsl-2
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dna0(ix,iy,is)=dna0(ix,iy,is)+dd(ispecies,istate,iy)*
     &                    (hy1(ix,iy)/hy1(ixref,iy))**2
           vla0(ix,iy,1,is)=vla0(ix,iy,1,is)+
     &        vvn(ispecies,istate,iy)*(hy1(ix,iy)/hy1(ixref,iy))
           dfhcb=chi(ispecies,istate,iy)*
     &           (hy1(ix,iy)/hy1(ixref,iy))**2
           hci0(ix,iy)=hci0(ix,iy)+dfhcb*na(ix,iy,is)
           hcib(ix,iy,is)=hcib(ix,iy,is)+dfhcb*na(ix,iy,is)
           hvi0(ix,iy,1)=hvi0(ix,iy,1)+vwt(ispecies,istate,iy)*
     &           (hy1(ix,iy)/hy1(ixref,iy))*na(ix,iy,is)/
     &           (ni(ix,iy,0)+ni_ext(ix,iy,0))
          enddo
         enddo
         dfhce=chi(nsl-1,0,iy)*(hy1(ix,iy)/hy1(ixref,iy))**2
         if (dfhce.lt.0) then
           write(*,*) 'dfhce(',ix,',',iy,')=',dfhce
         endif
         hce0(ix,iy)=hce0(ix,iy)+dfhce*ne(ix,iy)
         hve0(ix,iy,1)=hve0(ix,iy,1)+
     &      vwt(nsl-1,0,iy)*(hy1(ix,iy)/hy1(ixref,iy))

*   ..set constant anomalous transport in the SOL and divertor regions
        elseif (.not.on_closed_surface(ix,iy)) then
         do ispecies=0,nsl-2
          do istate=0,ncl(ispecies)-1
           is=isort(ispecies,istate)
           dna0(ix,iy,is)=dna0(ix,iy,is)+fix_dna
           dfhcb=fix_hci
           hci0(ix,iy)=hci0(ix,iy)+dfhcb*na(ix,iy,is)
           hcib(ix,iy,is)=hcib(ix,iy,is)+dfhcb*na(ix,iy,is)
          enddo
         enddo
         dfhce=fix_hce
         hce0(ix,iy)=hce0(ix,iy)+dfhce*ne(ix,iy)
        endif

       enddo
      enddo
*   ..return
      ncall = ncall+1
      call subend ()
      return

      end subroutine set_transport_neo






      subroutine set_profile(fnm, n, coord, fun)                        !srv 01.03.07 {
      use b2mod_types
      implicit none
      character*(*) fnm
      integer :: n, iy
      real (kind=R8) :: coord(-1:n), fun(-1:n)

      real (kind=R8), allocatable, save :: r(:), f(:), b(:), c(:), d(:)
      integer :: ndat, i
      real (kind=R8) :: t0, t1
      real (kind=R8), external :: dseval

*   ..read data in order to determine ndat
      open(70,file=trim(fnm))
      ndat = 0
      do
        read(70,*,end=10) t0, t1
        ndat = ndat + 1
      enddo
   10 close(70)

      write(*,"(/' The length of file ',a,' is ',i4/)") trim(fnm), ndat

      allocate (r(1:ndat), f(1:ndat), b(1:ndat), c(1:ndat), d(1:ndat))

      open(70,file=trim(fnm))
      do iy = 1, ndat
        read(70,*) r(iy), f(iy)
      enddo
      close(70)
      do iy = 1, ndat-1
        call xertst (r(iy).lt. r(iy+1),
     1  'Sequence of coordinates in file '//fnm//' is not monotonic')
      enddo

!            call dcscon (ndat, r, f, ibreak, break, cscoef)
!            write(*,*) ' ibreak = ', ibreak
!            ideriv = 0
!            nintv  = ndat-1
!            call dcs1gd (ideriv, ndat, r, nintv, break, cscoef, pf)
!            do iy = 1,ndat
!             pf(iy-2) = dcsval(r(iy), nintv, break, cscoef)
!             write(*,*) r(iy), pf(iy-2)
!            enddo

!        CALL DSPLINE(NDAT, R, F, B, C, D)
*   ..the simplest linear interpolation
      do i = 1, ndat-1
        b(i) = (f(i+1)-f(i))/(r(i+1)-r(i))
      enddo
      b(ndat) = b(ndat-1)
      c = 0.0e0_R8
      d = 0.0e0_R8
      do i = -1, n
        fun(i) = dseval(ndat, coord(i), r, f, b, c, d)
        write(*,*) coord(i), fun(i)
      enddo

      deallocate(r, f, b, c, d)

      end subroutine set_profile                                      !srv 01.03.07 }
#endif

c**************************************************************************************
c*****************     New KU Leuven transport model for neutrals *********************
c**************************************************************************************
      subroutine set_transport_kul(nCv, nFc,
     &            nVx, ns, nscx, iscx, tn_style, geo, 
     &   mpg, pl, dv, rt, st_ext, dna0, dpa0, vla0, vma0, vsa0, hci0,
     &   hcn0,hcib)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_rates
      use b2mod_transport_namelist
      use b2us_geo
      use b2us_map
      use b2us_plasma
      implicit none

*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, nscx, iscx(0:nscxmax-1)
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (inout) :: dv
      type (B2Rates), intent (in) :: rt
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) :: t0, t_av
*   ..input/output arguments
      real (kind=R8), intent(inout) ::
     &  hci0(nCv), hcn0(nCv), hcib(nCv,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dna0(nCv,0:ns-1), dpa0(nCv,0:ns-1),
     &  vla0(nCv,0:1,0:ns-1), vsa0(nCv,0:ns-1),
     &  vma0(nCv,0:1,0:ns-1)                              !srv 15.12.05
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(nCv), dcx, vcx, dion
      integer is, k, ic, iCv, tn_style
      intrinsic max, exp
      external expu
*   ..statement functions
      real (kind=R8) :: expu

*   ..subprogram start-up calls
      call subini ('set_transport_kul')
      	do is = 0, ns-1
	  if (is_neutral(is).and.zn(is).eq.1.0_R8) then
            k=0
            do while (iscx(k).ne.is .and. k.lt.nscx)
              k = k+1
            enddo

            call xertst (k.lt.nscx, 'CX species index not found!')
          endif

          do iCv = 1, nCv
            if (is_neutral(is).and.zn(is).eq.1.0_R8) then !! only for H D T
              dcx=0.0_R8
              vcx=0.0_R8
              do ic = 0, ns-1
                t_av = 0.5_R8*(pl%ti(iCv)+pl%tn(iCv))
                dcx = dcx + expu(rt%rlcx(iCv,0,ic,k) +
     &            rt%rlcx(iCv,1,ic,k)*log(t_av/(am(is)*ev)))
                vcx = vcx + pl%na(iCv,ic) * 
     &            expu(rt%rlcx(iCv,0,ic,k) +
     &            rt%rlcx(iCv,1,ic,k)*log(t_av/(am(is)*ev)))
              enddo
              dion=expu(
     2         rt%rlsa(iCv,0,is)+
     3         rt%rlsa(iCv,1,is)*log(pl%te(iCv)/ev))
              dpa0(iCv,is)=1.0_R8/mp/am(is)/(vcx+dion*dv%ne(iCv))

*    ..For now, no flux limiters imposed.
*    ..For now, no neutral-neutral collisions inlcuded
              vcx = vcx + dion*dv%ne(iCv)
              vsa0(iCv,is) = pl%na(iCv,is)*pl%tn(iCv)/
     &         (vcx) 
              hcib(iCv,is) = 2.5_R8*pl%na(iCv,is)*pl%tn(iCv)/
     &         (mp*am(is)*(vcx))
              
            endif
          enddo !! end loop over cells



*   ..Check whether dna0, vla0 and vma0 are zero for neutrals or print warning
        if (is_neutral(is).and.zn(is).eq.1.0_R8) then
          if (ANY( dna0(:,is)/=0 )) then
            write(*,*) 'Warning: '//
     .       'dna0 is not zero for neutrals'//
     .       'Recommended choice is zero!'
          endif
          if (ANY( vla0(:,:,is)/=0 )) then
            write(*,*) 'Warning: '//
     .       'vla0 is not zero for neutrals'//
     .       'Recommended choice is zero!'
          endif
          if (ANY( vma0(:,:,is)/=0 )) then
            write(*,*) 'Warning: '//
     .       'vma0 is not zero for neutrals'//
     .       'Recommended choice is zero!'
          endif
        endif 
      enddo !! end loop over species

c recompute hci0 and hcn0
      hci0 = 0.0_R8
      hcn0 = 0.0_R8
      do is =0, ns-1
        if (.not.is_neutral(is)) then
          hci0 = hci0 + hcib(:,is)
        else
          hcn0 = hcn0 + hcib(:,is)
        endif
      enddo
      if (tn_style.eq.0) then     !! default combined ion-neutral, add neutral contribution to "total" hci0
        hci0 = hci0+hcn0
      endif

*   ..return
      call subend ()
      return
      end subroutine set_transport_kul

      subroutine set_transport_keps(
     &     nCv, nFc, nVx, ns, ismain, switch, geo, mpg,
     &     pl, dv, rt, dna0, vsa0, hce0, hci0, hcib, sig0, alf0,
     &     dkt0, dzt0, dna_ExB, hce_ExB, hci_ExB)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_switches
      use b2mod_user_namelist
     , , only : nomp, omp, icsep
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : b2tqna_keps_eps, ncall_transp_keps
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns, ismain
      type (switches), intent (inout) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Rates), intent (in) :: rt
*   ..input/output arguments
      real (kind=R8) :: dna0(nCv,0:ns-1),
     &   vsa0(nCv,0:ns-1),
     &   hci0(nCv), hcib(nCv,0:ns-1),
     &   sig0(nCv), alf0(nCv),
     &   hce0(nCv), dkt0(nCv), dzt0(nCv),
     &   dna_ExB(nCv), hce_ExB(nCv),
     &   hci_ExB(nCv), rhol(nCv)
*   ..local variables
      integer :: is
      real(kind=R8) ::
     1   wrkf(nFc), wrkc(nCv), shear(nCv)

*   ..subprogram start-up calls
      call subini ('set_transport_keps')
      if(ncall_transp_keps.eq.0) then
        if (switch%keps_local.ne.1) then
         call xertst(nomp.gt.0,
     &    'No CVs in omp list, check rzomp in b2.user.parameters')
         call xertst(icsep.gt.0,
     &    'Invalid icsep value, check rzomp in b2.user.parameters')
        endif
      endif
      hci_ExB = 0.0_R8
      do is=0,ns-1
        if(.not.is_neutral(is)) then
          if(switch%keps_local.eq.1) then
*   ..compute local larmor radius
            rhol = am(ismain)*mp*sqrt(2.0_R8*pl%ti/(am(ismain)*mp))/
     &       abs(rt%rza(:,ismain)*qe*geo%cvBb(:,3))
          else
            rhol = am(ismain)*mp*sqrt(2.0_R8*pl%ti(omp(icsep))/
     &            (am(ismain)*mp))/abs(rt%rza(omp(icsep),ismain)*qe*
     &             geo%cvBb(omp(icsep),3))
          endif
*       ..compute radial shear of diamagnetic ExB velocity
          wrkf = dv%vaecrb(:,0,ismain)*geo%fcBb(:,3)/geo%fcBb(:,2)
          call intcell(nFc,nCv,mpg,mpg%intcellR,wrkf,wrkc)
          call gradc_div_r(nCv, nFc, nVx, 1, geo, mpg, wrkc, wrkf,
     &      shear)
          if(switch%transport_keps.eq.1) then
*   ..compute D according to KUL, using kt only
cwdk at the moment: assumption that kt is related to main ion species;
cwdk same dna_ExB for all ion species
            dna_ExB = switch%keps_cd*(pl%kt/(am(ismain)*mp))/
     &      (sqrt(pl%kt/(am(ismain)*mp))/rhol + switch%keps_shear*
     &       abs(shear) + b2tqna_keps_eps)
            dna0(:,is) = (dna_ExB + switch%dna_min)*switch%keps_fac
     &        +(1.0_R8-switch%keps_fac)*dna0(:,is)
          elseif(switch%transport_keps.eq.2) then
*   ..compute D according to KUL, using kt and zt
            dna_ExB = switch%keps_cd*(pl%kt/(am(ismain)*mp))/
     &      (sqrt(pl%zt/(am(ismain)*mp)) + switch%keps_shear*
     &       abs(shear) + b2tqna_keps_eps)
            dna0(:,is) = (dna_ExB + switch%dna_min)*switch%keps_fac
     &        +(1.0_R8-switch%keps_fac)*dna0(:,is)
          else
            call xerrab ('wrong value transport_keps')
          endif
          if (switch%keps_visc.gt.0.0_R8) then
            vsa0(:,is) = switch%keps_visc*am(is)*mp*pl%na(:,is)*
     &        (dna_ExB + switch%dna_min)*switch%keps_fac + 
     &        (1.0_R8-switch%keps_fac)*vsa0(:,is)
          endif
          if (switch%keps_heat.gt.0.0_R8) then
            hce_ExB = switch%keps_heat*dna_ExB*dv%ne
            hce0 = switch%keps_heat*(dna_ExB + switch%hce_min)*dv%ne*
     &         switch%keps_fac + (1.0_R8-switch%keps_fac)*hce0
          endif
          if (switch%keps_heat_i.gt.0.0_R8) then
            hci_ExB = hci_ExB + switch%keps_heat_i*dna_ExB*pl%na(:,is)
            hcib(:,is) = switch%keps_heat_i*switch%keps_fac*
     &        (dna_ExB + switch%hci_min)*
     &        pl%na(:,is) + (1.0_R8-switch%keps_fac)*hcib(:,is)
          endif
          if (switch%keps_sig.gt.0.0_R8) then
            sig0 = switch%keps_sig*(dna_ExB + switch%dna_min)*qe*
     &         dv%ne(omp(icsep))*
     &         switch%keps_fac + (1.0_R8-switch%keps_fac)*sig0
          endif
          if (switch%keps_alf.gt.0.0_R8) then
            alf0 = switch%keps_alf*(dna_ExB + switch%dna_min)*
     &         dv%ne(omp(icsep))*sqrt(qe/pl%te)*
     &         switch%keps_fac + (1.0_R8-switch%keps_fac)*sig0
          endif
          dkt0 = switch%keps_dkt*dna0(:,ismain)*pl%na(:,ismain)
          dzt0 = switch%keps_dzt*dna0(:,ismain)*pl%na(:,ismain)
        endif
      enddo
      hci0 = 0.0_R8
      do is = 0, ns-1
        hci0 = hci0 + hcib(:,is)
      enddo
 
*   ..increase percentage of coefficients from k-eps
      switch%keps_fac = min(switch%keps_fac*switch%keps_inc,1.0_R8)
      if (switch%keps_inc.gt.1.0_R8) 
     > write(*,*) 'b2tqna_keps_fac = ',switch%keps_fac

      if (switch%keps_iout.eq.1) then
        call my_out_us(70,nCv,0,rhol,'b2tqna_keps_rhol')
        call my_out_us(70,nCv,0,shear,'b2tqna_keps_shear')
      endif

*   ..return
      ncall_transp_keps = ncall_transp_keps +1
      call subend ()
      return

      end subroutine set_transport_keps

