*-----------------------------------------------------------------------
*.specification

      subroutine b2tqna (nx, ny, ns, nscx, iscx, bb, na,
     &   te, ti, ni, ne, ne2, hx, hy, qz, dna0, dpa0, vla0, vma0,    !srv 15.12.05
     &   vsa0, hci0, hcib, hce0, sig0, alf0, wrk0)
      use b2mod_types
      use b2mod_rates
      use b2mod_transport_namelist
      use b2mod_input_profile
      use b2mod_constants
      use b2mod_time
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_b2cmpb                                               !srv 02.12.04
      use b2mod_geo_corner
     &                    , only : lengthy, nmdpl                    !srv 26.11.04
#ifdef PAR_OPT
	use b2mod_par_opt
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, iscx(0:nscxmax-1)
      real (kind=R8) ::
     *  bb(-1:nx,-1:ny,0:3), na(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1),
     &  ne(-1:nx,-1:ny), ne2(-1:nx,-1:ny), 
     &  hx(-1:nx,-1:ny), hy(-1:nx,-1:ny), qz(-1:nx,-1:ny,0:1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &  vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &  hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1),
     &  vma0(-1:nx,-1:ny,0:1,0:ns-1)                                  !srv 15.12.05
      real (kind=R8) ::
     *  hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..workspace arguments (unspecified on entry and on exit)
      real (kind=R8) ::
     *  wrk0(-1:nx,-1:ny)
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
*     - vma: anomalous velocity in momentum balance equation.        !srv 15.12.05
*     - vsa: viscosity.
*     - hci: atom thermal conductivity.
*     - hce: electron thermal conductivity.
*     - sig: electrical conductivity.
*     - alf: thermo-electric coefficient.
*     These coefficients occur in the equations as indicated below.
*     Particle transport.
*       fnax(is) = -dna0(is)*(d/dx)na(is)-dpa0(is)*(d/dx)pa(is)+
*         vla0*na(is)
*       fnay(is) = -dna0(is)*(d/dy)na(is)-dpa0(is)*(d/dy)pa(is)+
*         vla0*na(is)
*     Momentum transport.
*       momx(is) = na(is)*up(is)*vma0 - vsa0(is)*(d/dx)up(is)        !srv 15.12.05
*       momy(is) = na(is)*up(is)*vma0 - vsa0(is)*(d/dy)up(is)        !srv 15.12.05
*     Atom heat flux.
*       fhix = -hci0*(d/dx)ti
*       fhiy = -hci0*(d/dy)ti
*     Electron heat flux.
*       fhex = alf0*te*ehx-hce0*(d/dx)te
*       fhey = alf0*te*ehy-hce0*(d/dy)te
*     Electric current.
*       fchx = sig0*ehx-alf0*(d/dx)te
*       fchy = sig0*ehy-alf0*(d/dy)te
*     Here, ehx and ehy are the components of the modified electric
*     field: ehx = -(d/dx)po+(1/(qe*ne))*(d/dx)(ne*te),
*     ehy = -(d/dy)po+(1/(ne*te))*(d/dy)(ne*te). See Balescu /1/.
*     The particle and momentum transport is a property of each atomic
*     species, the other transport laws involve summation over all
*     species.
*
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer ncall, n2, ix, iy, is, ic, k, new_df0, ixref, diagno,
     1  transport_namelist, user_transport, transport_inputfile,
     2  jsep, jxi, jxa
	integer, save :: model_sig = 0, iout = 0                        !srv 02.12.04
	character chns*2                                                !srv 10.11.09
      real (kind=R8) ::
     *  ballooning, bb_ref, ballooning_rescale, pfr_rescale,
     1  divsol_rescale
      real (kind=R8) ::
     *  df0, dfvsb, dfhcb, dfhce, dfsig, dfalf,
     1  dcx, dion, vcx, lcx, lion, din, dnn,
     1  cutlo, cutll, poloidal_scale_factor
      real (kind=R8) ::
     *  t0, t1, min_df0, max_df0,
     *  min_df0_tmp, max_df0_tmp
      real (kind=R8), save :: cfvma = 0.0e0_R8                       !srv 15.12.05      
      integer elm_count
      save ncall, new_df0, ixref, transport_namelist, user_transport,
     1  ballooning, bb_ref, ballooning_rescale, elm_count, pfr_rescale,
     2  divsol_rescale, jxi, jxa, jsep, min_df0, max_df0, diagno
      logical catch_up
*   ..procedures
      intrinsic sqrt
      external subini, subend, xertst, sfill, machsfr, ipgeti, ipgetr
      real (kind=R8) ::
     *  expu, addup, machsfr, hy1
*   ..initialisation
      data ncall/0/, new_df0/0/, transport_namelist/0/, diagno/0/,
     1  user_transport/0/, transport_inputfile/0/,
     2  ballooning/0.0_R8/, bb_ref/0.0_R8/, ballooning_rescale/1.0_R8/,
     3  pfr_rescale/1.0_R8/, divsol_rescale/1.0_R8/,
     4  min_df0/0.0_R8/, max_df0/1.0e30_R8/
*   ..statement functions
      addup(t0,t1)=t0+t1
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
      expu(t0) = dim(exp(max(t0,cutll)),cutlo)
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tqna')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
        call get_jsep(nx,ny,jxi,jxa,jsep)
        ixref=jxa
        call ipgeti ('b2tqna_diagno', diagno)
        call ipgeti ('b2tqna_new_df0', new_df0)
        call ipgetr ('b2tqna_min_df0', min_df0)
        call ipgetr ('b2tqna_max_df0', max_df0)
        call xertst (0.0_R8.le.min_df0.and.min_df0.le.max_df0.and.
     &   max_df0.gt.0.0_R8,
     *   'Invalid range for b2tqna_min_df0 and/or b2tqna_max_df0')
        call ipgeti ('b2tqna_user_transport', user_transport)
        call ipgeti ('b2tqna_inputfile', transport_inputfile)
        call ipgeti ('b2tqna_ixref', ixref)
        call xertst (-1.le.ixref.and.ixref.le.nx,
     1   'faulty argument b2tqna_ixref')
        call ipgeti ('b2tqna_transport_namelist', transport_namelist)
        if(transport_namelist.eq.1) then
          call alloc_b2mod_transport(ns)
          call read_b2mod_transport_namelist
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
        call ipgetr ('b2tqna_ballooning', ballooning)
        call ipgetr ('b2tqna_bb_ref', bb_ref)
        call ipgetr ('b2tqna_ballooning_rescale', ballooning_rescale)
        if(ballooning.ne.0.0_R8.and.bb_ref.eq.0.0_R8) then
          bb_ref=sum(bb(:,:,3))/((nx+2)*(ny+2))
          write(*,*) 'b2tqna: bb_ref set to ',bb_ref
        endif
        call ipgetr ('b2tqna_pfr_rescale', pfr_rescale)
        call xertst (pfr_rescale.gt.0.0_R8, 'faulty input pfr_rescale')
        call ipgetr ('b2tqna_divsol_rescale', divsol_rescale)
        call xertst (divsol_rescale.gt.0.0_R8,
     1   'faulty input divsol_rescale')
        call ipgetr ('b2tqna_cfvma', cfvma)                             !srv 15.12.05
        call ipgeti ('b2tqna_model_sig', model_sig)                     !srv 02.12.04
        call ipgeti ('b2tqna_iout', iout)                               !srv 10.11.09
      else
        if(transport_namelist.eq.1.and.
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
*   ..set cutll, cutlo
      cutll = log(machsfr())
      cutlo = max(machsfr(),exp(cutll))
*   ..compute auxiliary wrk0
      do iy = -1, ny
       do ix = -1, nx
        wrk0(ix,iy) = sqrt(ti(ix,iy)/mp)
       enddo
      enddo

      if(transport_namelist.eq.1) then
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
	cfhce(0)   = par_opt_phys(1)                                 !srv 05.04.09
	cfhci(0,1) = par_opt_phys(2)                                 !srv 05.04.09
	cfdna(0,1) = par_opt_phys(3)
	cfvsa(0,1) = cfdna(0,1)                                      !srv 22.07.10
#endif
* ..compute anomalous transport coefficients
*   ..compute dna0, dpa0, vla0, vma0, vsa0, hci0
*     (initialise hci0 to 0)
      call sfill (n2, 0.0_R8, hci0, 1)
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
        do iy = -1, ny
          do ix = -1, nx
*      ..compute auxiliary quantity df0
            if (cfdf0(2,is).eq.0.0_R8) then
              df0 = cfdf0(0,is)+cfdf0(1,is)/(ni(ix,iy,0)/1.0e20_R8)
            else
              df0 = cfdf0(0,is)+cfdf0(1,is)/(ni(ix,iy,0)/1.0e20_R8)+
     &         (cfdf0(2,is)/sqrt(am(is)))*wrk0(ix,iy)/
     &         (cfdf0(3,is)*ni(ix,iy,0)+cfdf0(4,is)*ne(ix,iy))
            endif
            if(new_df0.eq.1) then
              if(is_neutral(is)) then
                dcx=0.0_R8
                vcx=0.0_R8
                do ic = 0, ns-1
                  dcx = dcx + expu(rlcx(ix,iy,0,ic,k) +
     &             rlcx(ix,iy,1,ic,k)*log(ti(ix,iy)/(am(is)*ev)))
                  vcx = vcx + na(ix,iy,ic) * expu(rlcx(ix,iy,0,ic,k) +
     &             rlcx(ix,iy,1,ic,k)*log(ti(ix,iy)/(am(is)*ev)))
                enddo
                lcx=wrk0(ix,iy)/sqrt(am(is))/vcx
                dion=expu(
     2           rlsa(ix,iy,0,is)+
     3           rlsa(ix,iy,1,is)*log(te(ix,iy)/ev))
                lion=wrk0(ix,iy)/sqrt(am(is))/ne(ix,iy)/dion
                din=(wrk0(ix,iy)/sqrt(am(is)))**2/
     1           (vcx+dion*ne(ix,iy))
c               if(lcx.gt.lion) then
c                 din=lion*wrk0(ix,iy)/sqrt(am(is))
c               else
c                 din=sqrt(lcx*lion)*wrk0(ix,iy)/sqrt(am(is))
c               endif
c-wrong(Bonnin,Owen)                dnn=1.0_R8/na(ix,iy,is)/2.0e-19_R8
		dnn=(wrk0(ix,iy)/sqrt(am(is)))/na(ix,iy,is)/2.0e-19_R8
                if (diagno.gt.1) then
                  write(*,*) 'new_df0',
     1             ix,iy,is,lcx,lion,din,dnn,addup(din,dnn),df0
                else if (addup(din,dnn)/df0.gt.100.0_R8 
     &           .and. diagno.eq.1) then
                  write(*,*) 'new_df0:problem',
     1             ix,iy,is,lcx,lion,din,dnn,addup(din,dnn),df0
                endif
                df0=addup(din,dnn)
              endif
            endif
            min_df0_tmp=min(min_df0_tmp,df0)
            max_df0_tmp=max(max_df0_tmp,df0)
            df0=max(min_df0,min(max_df0,df0))
*      ..compute dna0, dpa0, vla0, vma0
            dna0(ix,iy,is) = cfdna(0,is)+
     &       cfdna(1,is)/(max(1.0e-50_R8,ne(ix,iy))/1.0e20_R8)+    ! dpc changed na(ix,iy,is) to ne(ix,iy)
     &       cfdna(2,is)*df0+
     &       cfdna(3,is)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
            if(cfdna(7,is).ne.0.0_R8)
     &       dna0(ix,iy,is) = dna0(ix,iy,is)*
     &        cfdna(7,is)*(hy1(ix,iy)/hy1(ixref,iy))**2
            dpa0(ix,iy,is) = 
     &       (1.0_R8/(rza(ix,iy,is)*te(ix,iy)+ti(ix,iy)))*
     &       (cfdpa(0,is)+cfdpa(1,is)/
     &       (max(1.0e-50_R8,ne(ix,iy))/1.0e20_R8)+    ! dpc changed na(ix,iy,is) to ne(ix,iy)
     &        cfdpa(2,is)*df0+
     &        cfdpa(3,is)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3)))
            if(cfdpa(7,is).ne.0.0_R8)
     &       dpa0(ix,iy,is) = dpa0(ix,iy,is)*
     &        cfdpa(7,is)*(hy1(ix,iy)/hy1(ixref,iy))**2
            vla0(ix,iy,0,is) = 0.0_R8
            vla0(ix,iy,1,is) = cfvla(0,is)+
     &       cfvla(1,is)/(max(1.0e-50_R8,ne(ix,iy))/1.0e20_R8)+    ! dpc changed na(ix,iy,is) to ne(ix,iy)
     &       cfvla(2,is)*df0+
     &       cfvla(3,is)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
            if(cfvla(7,is).ne.0.0_R8)
     &       vla0(ix,iy,1,is) = vla0(ix,iy,1,is)*
     &        cfvla(7,is)*hy1(ix,iy)/hy1(ixref,iy)
	      if(.not.is_neutral(is)) then                             !srv 15.12.05 {
	        vma0(ix,iy,0,is) = 0.0_R8
	        vma0(ix,iy,1,is) = 0.0_R8
	      else
	        vma0(ix,iy,0,is) = 0.0_R8
!             vma0(ix,iy,1,is) = cfvma(0,is)
              vma0(ix,iy,1,is) = cfvma                               !srv 15.12.05 }
	      endif
*      ..compute dfvsb and dfhcb
*        (dfvsb and dfhcb specify effective diffusivities, of
*         dimension m**2/s, for the anomalous transport coefficients
*         vsa0 and hci0)
            dfvsb = cfvsa(0,is)+
     &       cfvsa(1,is)/(max(1.0e-50_R8,ne(ix,iy))/1.0e20_R8)+    ! dpc changed na(ix,iy,is) to ne(ix,iy)
     &       cfvsa(2,is)*df0+
     &       cfvsa(3,is)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
            if(cfvsa(7,is).ne.0.0_R8)
     &       dfvsb = dfvsb*cfvsa(7,is)*(hy1(ix,iy)/hy1(ixref,iy))**2
            dfhcb = cfhci(0,is)+
     &       cfhci(1,is)/(max(1.0e-50_R8,ne(ix,iy))/1.0e20_R8)+    ! dpc changed na(ix,iy,is) to ne(ix,iy)
     &       cfhci(2,is)*df0+
     &       cfhci(3,is)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
            if(cfhci(7,is).ne.0.0_R8)
     &       dfhcb = dfhcb*cfhci(7,is)*(hy1(ix,iy)/hy1(ixref,iy))**2
*      ..compute vsa0, hci0
            vsa0(ix,iy,is) = dfvsb*mp*am(is)*na(ix,iy,is)
            hci0(ix,iy) = hci0(ix,iy)+dfhcb*na(ix,iy,is)
            hcib(ix,iy,is)=dfhcb*na(ix,iy,is)
          enddo
        enddo
!        write(*,*) 'is,min_df0,max_df0 ',is, min_df0_tmp, max_df0_tmp  !srv 18.05.09
      enddo
*   ..compute hce0, sig0, alf0
      do iy = -1, ny
       do ix = -1, nx
*     ..compute dfhce, dfsig, dfalf
*       (dfhce, dfsig, and dfalf specify effective diffusivity,
*        dimension m**2/s, for the anomalous transport coefficients
*        hce0, sig0, and alf0)
         dfhce = cfhce(0)+cfhce(1)/(ne(ix,iy)/1.0e20_R8)+
     1    cfhce(3)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
         if(cfhce(7).ne.0.0_R8)
     &    dfhce=dfhce*cfhce(7)*(hy1(ix,iy)/hy1(ixref,iy))**2
         dfsig = cfsig(0)+cfsig(1)/(ne(ix,iy)/1.0e20_R8)+
     1    cfsig(3)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
     2   +cfsig(4)*exp(-cfsig(5)*(lengthy(ix,iy)/lengthy(ix,ny))**2) !srv 26.11.04
         if(cfsig(7).ne.0.0_R8)
     &    dfsig=dfsig*cfsig(7)*(hy1(ix,iy)/hy1(ixref,iy))**2
         dfalf = cfalf(0)+cfalf(1)/(ne(ix,iy)/1.0e20_R8)+
     1    cfalf(3)*6.25e-2_R8*te(ix,iy)/ev/abs(bb(ix,iy,3))
         if(cfalf(7).ne.0.0_R8)
     &    dfalf=dfalf*cfalf(7)*(hy1(ix,iy)/hy1(ixref,iy))**2
*     ..compute hce0, sig0, alf0
         hce0(ix,iy) = dfhce*ne(ix,iy)
	   if(model_sig.eq.0) then                !srv 02.12.04 {
           sig0(ix,iy) = dfsig*qe*ne(ix,iy)
	   else
!	     sig0(ix,iy) = dfsig*qe*(cbsna(0,1,1)/abs(cbsna(1,1,1))) ! use constant density as core region
	     sig0(ix,iy) = dfsig*qe*ne(nmdpl,-1)  !srv 10.06.09 ! use constant density as core region 
	   endif                                  !srv 02.12.04 } 
         alf0(ix,iy) = dfalf*ne(ix,iy)*sqrt(qe/te(ix,iy))
       enddo
      enddo
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
        call set_transport_srv(nx, ny, ns, bb, na, te, ti, ni, ne, ne2,
     &  hy, qz, dna0, dpa0, vla0, vma0, vsa0, hci0, hcib, hce0, sig0,  
     &  alf0, wrk0)
      else
        write(*,*) 'user_transport = ',user_transport,' not coded'
        stop
      endif
      if(transport_inputfile.eq.1) then
        call transport_input(nx, ny, ns, rza, na, ne, te, ti, hcib,
     &   dna0, dpa0, vla0, vsa0, hci0, hce0, sig0, alf0, ncall)
        do ix=-1,nx
           do iy=-1,ny
              hci0(ix,iy)=0.0_R8
              do is=0,ns-1
                hci0(ix,iy)=hci0(ix,iy)+hcib(ix,iy,is)
              enddo
           enddo
        enddo
        if(iout.ne.0) then                                            !srv 10.11.09 {
	   do is = 0,ns-1
	    write (chns,'(i2.2)') is
	    call my_out(70,nx,ny,dna0(-1,-1,is),'b2tqna_in_dna0'//chns)
	    call my_out(70,nx,ny,dpa0(-1,-1,is),'b2tqna_in_dpa0'//chns)
	    call my_out(70,nx,ny,vsa0(-1,-1,is),'b2tqna_in_vsa0'//chns)
	    call my_out(70,nx,ny,vla0(-1,-1,0,is),'b2tqna_in_vla0x'//chns)
	    call my_out(70,nx,ny,vla0(-1,-1,1,is),'b2tqna_in_vla0y'//chns)
	    call my_out(70,nx,ny,hcib(-1,-1,is),'b2tqna_in_hcib'//chns)
	   enddo
	   call my_out(70,nx,ny,hci0,'b2tqna_in_hci0')
	   call my_out(70,nx,ny,hce0,'b2tqna_in_hce0')
	   call my_out(70,nx,ny,sig0,'b2tqna_in_sig0')
	   call my_out(70,nx,ny,alf0,'b2tqna_in_alf0')
	  endif                                                           !srv 10.11.09 }
      endif
ccc-dpc
      if(ballooning.ne.0.0_R8) then
        do iy=-1,ny
          do ix=-1,nx
            poloidal_scale_factor=
     1       abs(bb_ref/bb(ix,iy,3))**ballooning*
     1       ballooning_rescale
            do is=0,ns-1
              if(.not.is_neutral(is)) then
                dna0(ix,iy,is)=dna0(ix,iy,is)*poloidal_scale_factor
                dpa0(ix,iy,is)=dpa0(ix,iy,is)*poloidal_scale_factor
                vla0(ix,iy,:,is)=vla0(ix,iy,:,is)*poloidal_scale_factor
                vsa0(ix,iy,is)=vsa0(ix,iy,is)*poloidal_scale_factor
                hci0(ix,iy)=hci0(ix,iy)+
     1           hcib(ix,iy,is)*(poloidal_scale_factor-1.0_R8)
                hcib(ix,iy,is)=hcib(ix,iy,is)*poloidal_scale_factor
              endif
            enddo
            hce0(ix,iy)=hce0(ix,iy)*poloidal_scale_factor
            sig0(ix,iy)=sig0(ix,iy)*poloidal_scale_factor
            alf0(ix,iy)=alf0(ix,iy)*poloidal_scale_factor
          enddo
        enddo
      endif
      if(pfr_rescale.ne.1.0_R8) then
        do iy=-1,maxval(topcut)
          do ix=-1,nx
            if(region(ix,iy,0).ne.0.and.
     1       mod(region(ix,iy,0),4).ne.1.and.
     1       mod(region(ix,iy,0),4).ne.2.and.
     2       (((ix.lt.leftcut(1).or.ix.ge.rightcut(1)).and.
     3       iy.lt.topcut(1)).or.(nnreg(0).eq.8.and.
     4       ix.ge.leftcut(2).and.ix.lt.rightcut(2).and.
     5       iy.lt.topcut(2)))) then
              do is=0,ns-1
                if(.not.is_neutral(is)) then
                  dna0(ix,iy,is)=dna0(ix,iy,is)*pfr_rescale
                  dpa0(ix,iy,is)=dpa0(ix,iy,is)*pfr_rescale
                  vla0(ix,iy,:,is)=vla0(ix,iy,:,is)*pfr_rescale
                  vsa0(ix,iy,is)=vsa0(ix,iy,is)*pfr_rescale
                  hci0(ix,iy)=hci0(ix,iy)+
     1             hcib(ix,iy,is)*(pfr_rescale-1.0_R8)
                  hcib(ix,iy,is)=hcib(ix,iy,is)*pfr_rescale
                endif
              enddo
              hce0(ix,iy)=hce0(ix,iy)*pfr_rescale
              sig0(ix,iy)=sig0(ix,iy)*pfr_rescale
              alf0(ix,iy)=alf0(ix,iy)*pfr_rescale
            endif
          enddo
        enddo
      endif
      if(divsol_rescale.ne.1.0_R8) then
        do iy=jsep+1,ny
          do ix=-1,nx
            if(region(ix,iy,0).ne.0.and.
     1       mod(region(ix,iy,0),4).ne.1.and.
     1       mod(region(ix,iy,0),4).ne.2.and.
     2       (((ix.lt.leftcut(1).or.ix.ge.rightcut(1)).and.
     3       iy.ge.topcut(1)).or.(nnreg(0).eq.8.and.
     4       ix.ge.leftcut(2).and.ix.lt.rightcut(2).and.
     5       iy.ge.topcut(2)))) then
              do is=0,ns-1
                if(.not.is_neutral(is)) then
                  dna0(ix,iy,is)=dna0(ix,iy,is)*divsol_rescale
                  dpa0(ix,iy,is)=dpa0(ix,iy,is)*divsol_rescale
                  vla0(ix,iy,:,is)=vla0(ix,iy,:,is)*divsol_rescale
                  vsa0(ix,iy,is)=vsa0(ix,iy,is)*divsol_rescale
                  hci0(ix,iy)=hci0(ix,iy)+
     1             hcib(ix,iy,is)*(divsol_rescale-1.0_R8)
                  hcib(ix,iy,is)=hcib(ix,iy,is)*divsol_rescale
                endif
              enddo
              hce0(ix,iy)=hce0(ix,iy)*divsol_rescale
              sig0(ix,iy)=sig0(ix,iy)*divsol_rescale
              alf0(ix,iy)=alf0(ix,iy)*divsol_rescale
            endif
          enddo
        enddo
      endif
ccc-dpc
* ..return
      ncall = ncall+1
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2tqna

      end subroutine b2tqna

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
c-dpc     hci0(ix,iy)=d_new  ! cxpb *ni(ix,iy,0) ?
c-xpb     hcib(ix,iy,is)=d_new*na(in,is)
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
      use b2mod_geo_corner                                                 !srv 24.04.02
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
      real (kind=R8) :: fdnc(-1:nx,-1:ny,0:ns-1), wrk1(-1:nx,-1:ny)    !srv 13.11.02 26.11.04 01.03.07
      real (kind=R8) :: fvsc(-1:nx,-1:ny,0:ns-1), 
     &                  ftec(-1:nx,-1:ny), ftic(-1:nx,-1:ny,0:ns-1),  !srv 17.04.06 01.03.07
     &                  fvlc(-1:nx,-1:ny,0:1,0:ns-1)                  !srv 19.04.10
      real (kind=R8), save :: cdf = 0.0e0_R8, alf = 0.0e0_R8,         !srv 24.04.02 05.11.03
     &      alfa1 = 0.0e0_R8, alfa2 = 0.0e0_R8,
     &      delta1 = 0.0e0_R8, delta2 = 0.0e0_R8                      !srv 13.11.02 24.06.03
      real (kind=R8), save :: vs_alfa1 = 1.0e0_R8, vs_alfa2 = 1.0e0_R8,
     &                       vs_delta1 = 1.0e0_R8, vs_delta2 = 1.0e0_R8 !srv 17.04.06 {
      real (kind=R8), save :: te_alfa1 = 1.0e0_R8, te_alfa2 = 1.0e0_R8,
     &                       te_delta1 = 1.0e0_R8, te_delta2 = 1.0e0_R8
      real (kind=R8), save :: ti_alfa1 = 1.0e0_R8, ti_alfa2 = 1.0e0_R8,
     &                       ti_delta1 = 1.0e0_R8, ti_delta2 = 1.0e0_R8 !srv 17.04.06 }
      real (kind=R8), save :: vl_alfa1 = 1.0e0_R8, vl_alfa2 = 1.0e0_R8,
     &                       vl_delta1 = 1.0e0_R8, vl_delta2 = 1.0e0_R8 !srv 19.04.10
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
	character chns*2                                               !srv 05.11.03
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
c    Figure show values of coefficiets which are applyed to transport anomalous
c    coefficients in regions.
c    alfa1 is applied on cells which have radial coordinates less than delta1 [m] in outmidplane
c    alfa2 is applied jn cells which have radial coordinates great than -delta2 [m] in outmidplane
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
      call defdc(nx, ny, itype, ns, fdnc, fvsc, ftec, ftic)           !srv 13.11.02 17.04.06 01.03.07
	if(iout.ne.0 .and. ncall.eq.0) then                             !srv 09.10.01
	 do is = 0,ns-1
	  write (chns,'(i2.2)') is                                      !srv 09.10.01
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
	  write (chns,'(i2.2)') is
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

	 subroutine defdc(nx, ny, itype, ns, fdnc, fvsc, ftec, ftic)    !srv 13.11.02 01.03.07
	  implicit none
	  integer :: nx, ny, itype, ns
        real(kind=R8) :: fdnc(-1:nx,-1:ny,0:ns-1)                     !srv 17.04.06 {
        real(kind=R8) :: fvsc(-1:nx,-1:ny,0:ns-1)
        real(kind=R8) :: ftec(-1:nx,-1:ny), ftic(-1:nx,-1:ny,0:ns-1) !srv 17.04.06 }

        logical :: exist_file(0:2)
	  integer ix, iy, is, nc_c_m_f(1)                           !srv 05.07.04
	  integer, save :: ncall = 0, bottomdelta1, topdelta2, mintopcut !srv 28.07.10
        real (kind=R8) :: width, bbaverage, t0, t1, step
        real (kind=R8), allocatable, save :: coord(:), fun(:,:)

	  fdnc = 1.0e0_R8
	  fvsc = 1.0e0_R8                                           !srv 17.04.06 {
	  ftec = 1.0e0_R8
	  ftic = 1.0e0_R8                                           !srv 17.04.06 }
	  fvlc = 1.0e0_R8                                           !srv 19.04.10 }

	  select case (itype)
	   case(9)                                                !srv 22.05.09 {
	    mintopcut=minval(topcut(1:nncut))
*   ..modification of case (5). It is used barer function
	    do iy = -1, ny                  ! for SOL and PR region
	     do ix = -1, nx
	      if(mod(region(ix,iy,0),4).ne.1) then 
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
     &		  crdy_surf(nmdpl,mintopcut), 1.0_R8, delta2-step,
     &		  alfa2, delta2, alfa3)
	         fvsc(ix,iy,is)= barier2(crdy_c_cell(nmdpl,iy),
     &		  crdy_surf(nmdpl,mintopcut), 1.0_R8, vs_delta2-step,
     &		  vs_alfa2, vs_delta2, vs_alfa3)
	         ftic(ix,iy,is)= barier2(crdy_c_cell(nmdpl,iy),
     &		  crdy_surf(nmdpl,mintopcut), 1.0_R8, ti_delta2-step,
     &		  ti_alfa2, ti_delta2, ti_alfa3)
	        endif
	       enddo
	       ftec(ix,iy)= barier2(crdy_c_cell(nmdpl,iy),
     &		  crdy_surf(nmdpl,mintopcut), 1.0_R8, te_delta2-step,
     &		  te_alfa2, te_delta2, te_alfa3)
	      endif
	     enddo
	    enddo

          step = maxval(hy(nmdpl,-1:mintopcut-1))
	    do iy = -1, mintopcut-1   ! for inner surfaces 
	     do ix=-1,nx
	      if(mod(region(ix,iy,0),4).eq.1) then
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
	      stop
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

           inquire (file = 'b2tdnc.dat', exist = exist_file(0))
           inquire (file = 'b2thic.dat', exist = exist_file(1))
           inquire (file = 'b2thec.dat', exist = exist_file(2))
	     call xertst (all(exist_file),
     1	 'There are no b2tdnc.dat or b2thic.dat or b2thec.dat files')

           call set_profile('b2tdnc.dat', ny, coord, fun(-1,0))
 	     do iy = -1, ny
	      call xertst (fun(iy,0).gt. 0.0e0_R8,
     1	 'Approximation of data from file '//'b2tdnc.dat'//
     2      ' has gave '//'negative coefficient')
 	     enddo

           call set_profile('b2thic.dat', ny, coord, fun(-1,1))
 	     do iy = -1, ny
	      call xertst (fun(iy,1).gt. 0.0e0_R8,
     1	 'Approximation of data from file '//'b2thic.dat'//
     2      ' has gave '//'negative coefficient')
 	     enddo

           call set_profile('b2thec.dat', ny, coord, fun(-1,2))
 	     do iy = -1, ny-1
	      call xertst (fun(iy,2).gt. 0.0e0_R8,
     1	 'Approximation of data from file '//'b2thec.dat'//
     2      ' has gave '//'negative coefficient')
 	     enddo

           write(*,"(/' Exchanges of transport coefficients',
     &               ' have been done'/)")
	     do is = 0, ns-1
	      if(.not.is_neutral(is)) then
	       cfdna(2,is)   = 1.0e0_R8
	       cfhci(2,is)   = 1.0e0_R8
             write(*,"(' Now cfdna(2,',i2,') =',es10.2)")is, cfdna(2,is) 
             write(*,"(' Now cfhci(2,',i2,') =',es10.2)")is, cfhci(2,is) 
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
	        if(mod(region(ix,iy,0),4).eq.1) then
	         t0 = t0 +bb(ix,iy,3)*vol(ix,iy)/hy1(ix,iy)
	         t1 = t1 +vol(ix,iy)/hy1(ix,iy)
	        endif
	       enddo
	      enddo
            bbaverage=t0/t1
	      do iy = -1, ny                   ! for all surfaces 
	       do ix = -1, nx
	        if(mod(region(ix,iy,0),4).eq.1 .and. iy.le.mintopcut) then
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
!	     do ix = leftcut(1), rightcut(1)-1
	        do ix=0,nx
	         if(mod(region(ix,iy,0),4).eq.1) then
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
	      if(mod(region(ix,iy,0),4).ne.1) then 
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
	      if(mod(region(ix,iy,0),4).eq.1) then
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
     &                mod(region(ix,iy,0),4).ne.2) then   ! NOT SOL
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
	       if(mod(region(ix,iy,0),4).eq.1) then
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
	       if(mod(region(ix,iy,0),4).eq.1) then
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
*   ..inout arguments
      real (kind=R8) ::
     *   dna0(-1:nx,-1:ny,0:ns-1), dpa0(-1:nx,-1:ny,0:ns-1),
     &   vla0(-1:nx,-1:ny,0:1,0:ns-1), vsa0(-1:nx,-1:ny,0:ns-1),
     &   hci0(-1:nx,-1:ny), hcib(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *   hce0(-1:nx,-1:ny), sig0(-1:nx,-1:ny), alf0(-1:nx,-1:ny)
*   ..common blocks
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
      integer ncall, n2, ix, iy, is, ixref, iyref
      real (kind=R8) ::
     *   eta, current_eta, eta_floor, eta_ceiling, eta_alpha
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
        call read_b2mod_transport_models
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
          if(mod(region(ix,iy,0),4).le.2) then
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
      integer ncall, n2, ix, iy, is, ir, jxi, jxa, jsep
      real (kind=R8) ::
     *  alpha, min_df0, max_df0,
     *  wb, seglen, perimeter, L_perp, 
     *  pe(-1:nx,-1:ny), pri(-1:nx,-1:ny), pp(-1:nx,-1:ny), 
     *  coss(-1:nx,-1:ny), R, width(-1:nx), rzb(-1:nx,-1:ny)
      save 
*   ..procedures
      real (kind=R8) ::
     *  hy1
      external subini, subend, ipgetr, xertst, xerrab
      external b2xvsg, b2xppe, b2xppb
      intrinsic sin, log, mod, max, min
      hy1(ix,iy)=hy(ix,iy)*qz(ix,iy,1)
*   ..initialisation
      data ncall/0/, alpha/1.0_R8/, min_df0/0.0_R8/, max_df0/1.0e30_R8/

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
        call ipgetr ('b2tqna_min_df0', min_df0)
        call ipgetr ('b2tqna_max_df0', max_df0)
        call xertst (0.0_R8.le.min_df0.and.min_df0.le.max_df0.and.
     &   max_df0.gt.0.0_R8,
     *   'Invalid range for b2tqna_min_df0 and/or b2tqna_max_df0')
        call get_jsep(nx,ny,jxi,jxa,jsep)
      endif
*   ..compute plasma pressure from charged species only
      call b2xppe (nx, ny, ns, ne, te, pe)
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
        if (mod(region(ix,0,0),4).eq.1) then
          ir = ir + 1
          perimeter = perimeter + hx(ix,0)
          R = R + 0.5_R8*(crx(ix,0,0)+crx(ix,0,1))*hx(ix,0)
          iy = 0
          width(ix) = 0.5_R8 * hy1(ix,0)
          do iy = 0, jsep-1
            width(ix) = width(ix) + 0.5_R8 * (hy1(ix,iy-1)+hy1(ix,iy))
          end do 
          width(ix) = width(ix) + 0.5_R8 * hy1(ix,iy)
          if (pp(ix,0).gt.pp(ix,iy)) then
            L_perp = L_perp + width(ix) / log(pp(ix,0)/pp(ix,iy))
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
              dna0(ix,iy,is) = max(min_df0,min(max_df0,
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
     &        (iy-(mw_jsep+1))+vsep
            elseif(mw_interpol.eq.1.and.ix.le.mw_jxa) then !sin poloidal
              vla0(ix,iy,1,1)=(((vwall-vsep)/(ny-(mw_jsep+1))*
     &        (iy-(mw_jsep+1)))+vsep)*
     &        sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
            elseif(mw_interpol.eq.1.and.mw_jxa.lt.ix) then !sin poloidal
              vla0(ix,iy,1,1)=(((vwall-vsep)/(ny-(mw_jsep+1))*
     &        (iy-(mw_jsep+1)))+vsep)*sin(pi/2.0_R8*(ix-xbottom)/
     &        (mw_jxa-xbottom))
            elseif(mw_interpol.eq.2) then !sin radial
              vla0(ix,iy,1,1) = (vwall-vsep)*sin(pi/2.0_R8*(iy-mw_jsep)/
     &        (ny-mw_jsep))+vsep
            elseif(mw_interpol.eq.3) then
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
            elseif(mw_interpol.eq.4.and.ix.le.mw_jxa) then !sin poloidal
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &        sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
            elseif(mw_interpol.eq.4.and.mw_jxa.lt.ix) then !sin poloidal
              vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &        sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
            elseif(mw_interpol.eq.5.and.ix.le.mw_jxa) then
              if(iy.lt.ny) then
                if(ne(ix,iy) < ne(ix,iy+1)) then
                 vla0(ix,iy,1,1) = 0.9_R8*vla0(ix,iy-1,1,1)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &           sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
                endif
              else
                if(ne(ix,ny).gt.ne(ix,ny-1)) then
                 vla0(ix,ny-1,1,1) = 0.9_R8*vla0(ix,ny,1,1)
                 vla0(ix,ny,1,1) = vla0(ix,ny-1,1,1)
                else
                 vla0(ix,ny,1,1) = mw_flux_wall / ne(ix,ny)*
     &           sin(pi/2.0_R8*(ix-xtop)/(mw_jxa-xtop))
                endif
              endif
            elseif(mw_interpol.eq.5.and.mw_jxa.lt.ix) then
              if(iy.lt.ny) then
                if(ne(ix,iy-1).le.ne(ix,iy)) then
                 vla0(ix,iy-1,1,1) = 0.9_R8*vla0(ix,iy,1,1)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)*
     &           sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
                endif
              else
                if(ne(ix,ny).gt.ne(ix,ny-1)) then
                 vla0(ix,iy-1,1,1) = 0.9_R8*vla0(ix,iy,1,1)
                 vla0(ix,iy,1,1) = vla0(ix,iy-1,1,1)
                else
                 vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)*
     &           sin(pi/2.0_R8*(ix-xbottom)/(mw_jxa-xbottom))
                endif
              endif
            elseif(mw_interpol.eq.6.and.ix.le.mw_jxa) then
              if(iy.le.ny) then
                if(ne(ix,iy-1).le.ne(ix,iy)) then
                 vla0(ix,iy-1,1,1) = 0.9_R8*vla0(ix,iy,1,1)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
              endif
            elseif(mw_interpol.eq.6.and.mw_jxa.lt.ix) then
              if(iy.lt.ny) then
                if(ne(ix,iy).lt.ne(ix,iy+1)) then
                 vla0(ix,iy,1,1) = 0.9_R8*vla0(ix,iy-1,1,1)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
              endif 
c MW working                                 
            elseif(mw_interpol.eq.7.and.ix.le.mw_jxa) then
              if(iy.lt.ny) then
                if(ne(ix,iy).lt.ne(ix,iy+1)) then
                 vla0(ix,iy-1,1,1) = 0.9_R8*vla0(ix,iy,1,1)
                 vla0(ix,iy,1,1) = vla0(ix,iy-1,1,1)*ne(ix,iy-1)/
     &           ne(ix,iy)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                if(ne(ix,ny-1).le.ne(ix,ny)) then
                 vla0(ix,ny-1,1,1) = 0.9_R8*vla0(ix,ny,1,1)
                 vla0(ix,ny,1,1) = vla0(ix,ny-1,1,1)*ne(ix,ny-1)/
     &           ne(ix,ny)
                else
                 vla0(ix,ny,1,1) = mw_flux_wall / ne(ix,ny)
                endif
              endif
            elseif(mw_interpol.eq.7.and.mw_jxa.lt.ix) then
              if(iy.lt.ny) then
                if(ne(ix,iy).le.ne(ix,iy+1)) then
                 vla0(ix,iy-1,1,1) = 0.9_R8*vla0(ix,iy,1,1)
                 vla0(ix,ny,1,1) = vla0(ix,ny-1,1,1)*ne(ix,iy-1)/
     &           ne(ix,iy)
                else
                 vla0(ix,iy,1,1) = mw_flux / ne(ix,iy)
                endif
              else
                if(ne(ix,ny-1).le.ne(ix,ny)) then
                 vla0(ix,ny-1,1,1) = 0.9_R8*vla0(ix,ny,1,1)
                 vla0(ix,ny,1,1) = vla0(ix,iy-1,1,1)*ne(ix,ny-1)/
     &           ne(ix,ny)
                else
                 vla0(ix,iy,1,1) = mw_flux_wall / ne(ix,iy)
                endif
              endif 
c MW working end
c MW working 2                                
            elseif(mw_interpol.eq.8) then
              if(iy.lt.ny) then
                if(ne(ix,iy).le.ne(ix,iy+1)) then
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
                if(ne(ix,ny-1).le.ne(ix,ny)) then
                 vold = vla0(ix,ny,1,1)
                 vla0(ix,ny,1,1) = vold - vold * 0.5_R8 * mw_fluxalpha
                 if(vla0(ix,ny,1,1).lt.0.0_R8) then
                    vla0(ix,ny,1,1) = 0.0_R8
                 endif                 
                else
                  vla0(ix,ny,1,1) = (mw_flux_wall / ne(ix,ny) - 
     &            vla0(ix,ny,1,1)) * mw_fluxalpha + vla0(ix,ny,1,1)                  
                endif
              endif 
c MW working end 2
c MW working 3
c MW same as int8 + mw_flux also for cells at mw_jsep - transport
c through separatrix
C This is the one really used...                                
            elseif(mw_interpol.eq.9) then
              if(iy.lt.ny) then
                if(ne(ix,iy).le.ne(ix,iy+1)) then
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
                 if(ne(ix,iy-1).le.ne(ix,iy)) then
                   vold = vla0(ix,iy-1,1,1)
                   vla0(ix,iy-1,1,1) = 
     &              vold - vold * 0.5_R8 * mw_fluxalpha
                 if(vla0(ix,iy-1,1,1).lt.0.0_R8) then
                    vla0(ix,iy-1,1,1) = 0.0_R8
                 endif
                else
                  vla0(ix,iy-1,1,1) = (mw_flux / ne(ix,iy-1) - 
     &            vla0(ix,iy-1,1,1)) * mw_fluxalpha + vla0(ix,iy-1,1,1)
                endif               
                endif
              else
                if(ne(ix,ny-1).le.ne(ix,ny)) then
                 vold = vla0(ix,ny,1,1)
                 vla0(ix,ny,1,1) = vold - vold * mw_fluxalpha
                 if(vla0(ix,ny,1,1).lt.0.0_R8) then
                    vla0(ix,ny,1,1) = 0.0_R8
                 endif
                else
                  vla0(ix,ny,1,1) = (mw_flux_wall / ne(ix,ny) - 
     &            vla0(ix,ny,1,1)) * mw_fluxalpha + vla0(ix,ny,1,1)                  
                endif
              endif 
c MW working end 3
            else
             write(*,*) 'mw_interpol = ',mw_interpol,' not coded'
             stop
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
   10	 close(70)

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
!	      write(*,*) ' ibreak = ', ibreak
!            ideriv = 0
!	      nintv  = ndat-1
!            call dcs1gd (ideriv, ndat, r, nintv, break, cscoef, pf)
!	      do iy = 1,ndat 
!	       pf(iy-2) = dcsval(r(iy), nintv, break, cscoef)
!	       write(*,*) r(iy), pf(iy-2)
!	      enddo

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


