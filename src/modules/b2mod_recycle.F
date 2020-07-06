      module b2mod_recycle
      use b2mod_types
      use b2mod_constants
      use b2mod_wall
c**************************************************************************************
c*************   Module for new KU Leuven recycling model for neutrals  ***************
c**************************************************************************************
!mb   06/07/2017
c   This module is only loaded for the advanced recycling boundary conditions 
c     that make use of pre-evaluated integrals of TRIM data to consistently calculate 
c     recycling coefficients. To turn on the model, the flag b2stbr_kul_recycle is set to '1'   
c     Mixed materials are not considered in this recycling model

      logical, save :: recycle_initialised
      data recycle_initialised /.false./


      integer, save :: nm, nt, nshpottau,nx1,ny1,ns1

c     Recycle_index and target_index refer for each species and boundary in which integer
c     the corresponding recycling coefficients are stored. 
      integer, allocatable, save ::
     *  recycle_index(:),target_index(:)
c     Boolean to indicate whether trim data is present for species/material combination
      logical, allocatable, save :: 
     *  trimrecyc_present(:,:)
      real (kind=R8), allocatable, save ::
     &  M(:),T(:),dM(:),dT(:),shpottau(:),rc1(:,:,:,:,:),rc2(:,:,:,:,:),
     &  rmx1(:,:,:,:,:),rmx2(:,:,:,:,:),rmy1(:,:,:,:,:),rmy2(:,:,:,:,:),
     &  rmz1(:,:,:,:,:),rmz2(:,:,:,:,:),re1(:,:,:,:,:),re2(:,:,:,:,:),
     &  rp11(:,:,:,:,:),rp12(:,:,:,:,:),rp21(:,:,:,:,:),
     &  rp22(:,:,:,:,:),rp31(:,:,:,:,:),rp32(:,:,:,:,:),
     &  M_dummy(:), T_dummy(:), shpottau_dummy(:)        !nh 21.02.2020

      contains

      subroutine read_recycle_fnn(nsd,nbnd)
      use b2mod_b2cmpa
      use b2mod_elements
      implicit none
      integer,intent(in):: nsd, nbnd
      character*256 filename
      integer ib, isp, is, il,it,itm,ndata
c      integer nmax, itargsp, ic, i
      integer idum(0:9)
      logical open_file, found

      character*10::fnn_version
      character lblfnn*120
      character lblfnn_dummy*120              !nh 21.02.2020
      character*2 target_code(0:nbnd-1)
cmb size cannot become bigger than nsd-1?
c#ifndef PSR
c      character*2 recycle_code(0:nsd-1+6)
      character*2 recycle_code(0:nsd-1)
c#else
c      character*2 recycle_code(0:100)
c#endif
      external open_file, ipgeti
      external subini, subend, cfrure, cfruin, machsfr
      real (kind=R8) ::
     *  machsfr
      if(recycle_initialised) return
      call subini ('read_recycle_fnn')
!!     nx1 = nxd
!!     ny1 = nyd
      ns1 = nsd
      allocate(recycle_index(0:nsd-1))
      allocate(target_index(0:nbnd-1))
cmb   Check how many different recycled species there are
      recycle_index=-1
      recycle_code='  '
      isp=-1
      do is=0,nsd-1
        do il=0,isp
          if(is_codes(is).eq.recycle_code(il)) then
            recycle_index(is)=il
            goto 5
          endif
        enddo
        isp=isp+1
c#ifdef PSR
c        call xertst(isp.le.100, 
c     &   'read_recycle_fnn (1): increase size of recycle_code')
c#endif
        recycle_index(is)=isp
        recycle_code(isp)=is_codes(is)
    5   continue
        write(*,'(1x,i3,1x,i3,1x,a)') is,recycle_index(is),
     1   recycle_code(recycle_index(is))
      enddo

cmb   Check how many different surface materials
cmb   Mixed materials not taken into account.
      target_index=-1
      target_code='  '
      itm=-1
      do ib=0,nbnd-1
        do it=0,itm
          if(surface_material_name(ib+1).eq.target_code(it)) then !mb ibnd = iwall = ib+1 (arrays shifted)
            target_index(ib)=it
            goto 6
          endif
        enddo
        itm=itm+1
c#ifdef PSR
c        call xertst(itm.le.20, 
c     &   'read_recycle_fnn (1): increase size of recycle_code')
c#endif
        target_index(ib)=itm
        target_code(itm)=surface_material_name(ib+1)
    6   continue
        write(*,'(1x,i3,1x,i3,1x,a)') ib,target_index(ib),
     1   target_code(target_index(ib))
      enddo

      nM = -1
      nT = -1
      nshpottau = -1
      allocate(trimrecyc_present(0:isp,0:itm))
      trimrecyc_present = .false.
      do is=0,isp
        do it=0,itm
          filename = 
     1     recycle_code(is)(1:len_trim(recycle_code(is)))//
     2     target_code(it)(1:len_trim(target_code(it)))//
     3     '.fnn'
!!          filename = 'DC.fnn'
          call lower_case(filename)
          if(open_file(96,filename,
     1      'data.local/trim.data/refl.data|' // 
     2      'modules/B2.5/Database/trim.data/refl.data')) then
            trimrecyc_present(is,it) = .true.
            call cfverr (96,fnn_version)
            call cfruin (96, 3, idum, 'nm,nt,nshpottau')
c           Check whether array has same dimensions as previous, otherwise give error
            if (nM.eq.-1.and.nt.eq.-1.and.nshpottau.eq.-1) then
c           Allocate arrays with sizes from first read data file
              nM = idum(0)-1
              nT = idum(1)-1
              nshpottau = idum(2)-1
              ndata = (nM+1)*(nT+1)*(nshpottau+1)
c             ALLOCATE ARRAYS
              allocate(M(0:nM))
              allocate(T(0:nT))
              allocate(dM(0:nM-1))
              allocate(dT(0:nT-1))
              allocate(shpottau(0:nshpottau))
              allocate(rc1(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rc2(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmx1(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmx2(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmy1(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmy2(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmz1(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rmz2(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(re1(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(re2(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp11(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp12(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp21(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp22(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp31(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
              allocate(rp32(0:nM,0:nT,0:nshpottau,0:isp,0:itm))
c             Read label
              call cfruch (96,120,lblfnn,'label')
              call cfrure (96, nM+1, M, 'M')
              call cfrure (96, nT+1, T, 'T')
              call cfrure (96, nshpottau+1, shpottau, 'shpottau')
              dM(:) = 1.0_R8/(M(1:nM)-M(0:nM-1))
              dT(:) = 1.0_R8/(T(1:nT)-T(0:nT-1))
               
            else
              call xertst((nM.eq.idum(0)-1).and.(nT.eq.idum(1)-1).and.
     &          (nshpottau.eq.idum(2)-1), 
     &          'read_recycle_fnn: all .fnn files should have'//
     &          ' same matrix dimensions')
    
cnh  After checking the dimensions of M, T and shpottau, we should skip the 
c    first lines in the .fnn file. They are allocated to dummy variables.
              allocate(M_dummy(0:nM))
              allocate(T_dummy(0:nT))
              allocate(shpottau_dummy(0:nshpottau))
              
              call cfruch (96,120,lblfnn_dummy,'label')
              call cfrure (96, nM+1, M_dummy, 'M')
              call cfrure (96, nT+1, T_dummy, 'T')
              call cfrure (96, nshpottau+1, shpottau_dummy, 'shpottau')

            endif

c         Read preintegrated recycling coefficient contributions          
          call cfrure (96,ndata,rc1(:,:,:,is,it),'r1')
          call cfrure (96,ndata,rc2(:,:,:,is,it),'r2')
          call cfrure (96,ndata,rmx1(:,:,:,is,it),'rmx1')
          call cfrure (96,ndata,rmx2(:,:,:,is,it),'rmx2')
          call cfrure (96,ndata,rmy1(:,:,:,is,it),'rmy1')
          call cfrure (96,ndata,rmy2(:,:,:,is,it),'rmy2')
          call cfrure (96,ndata,rmz1(:,:,:,is,it),'rmz1')
          call cfrure (96,ndata,rmz2(:,:,:,is,it),'rmz2')
          call cfrure (96,ndata,re1(:,:,:,is,it),'re1')
          call cfrure (96,ndata,re2(:,:,:,is,it),'re2')
          call cfrure (96,ndata,rp11(:,:,:,is,it),'rp11')
          call cfrure (96,ndata,rp12(:,:,:,is,it),'rp12')
          call cfrure (96,ndata,rp21(:,:,:,is,it),'rp21')
          call cfrure (96,ndata,rp22(:,:,:,is,it),'rp22')
          call cfrure (96,ndata,rp31(:,:,:,is,it),'rp31')
          call cfrure (96,ndata,rp32(:,:,:,is,it),'rp32')
          close(96)
          else
c         Use standard recycling model for this species/wall combination
            write(*,*) 'No TRIM recycling data for ',
     1       recycle_code(is),' onto ',target_code(it)
          endif
        enddo 
      enddo
c     Other initializations

      recycle_initialised=.true.
      call subend ()
      return
      end subroutine read_recycle_fnn


      subroutine CalcIncidentFluxes(iCv,iCn,iFc,isign,
     &            isn,iscx0,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,                           
     &            rt,dnndy2, dnndz2,cosa,sina,
     &            dnndz,dnndzc,dnndzv,use_uy_uz_0,use_uy_uz_1,
     &            print_output,split_fluxes,
     &            fnni, fnnic,fnniv,fmomni,nni,nnwwni,feneni)
c     This routine calculates the incident neutral particle,momentum and energy fluxes based 
c     on the diffusion approximation (see papers of Niels Horsten)
c main reference: 
c   Horsten, Niels, Giovanni Samaey, and M. Baelmans.
c  "Development and assessment of 2D fluid neutral models that include atomic databases and a microscopic reflection model."
c   Nuclear Fusion 57.11 (2017): 116043.
      use b2mod_rates
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2rates), intent (in) :: rt
      integer,intent(in):: iCv,iCn, iFc,
     &    isn,iscx0,istra, use_uy_uz_0, use_uy_uz_1
      integer,intent(in)::print_output, split_fluxes  
      real (kind=R8),intent(in):: Phi_app, area, tnf, tif, tef,
     2    pof,nef,isign,dnndy2,dnndz2, dnndzc, dnndzv, cosa, sina
      real (kind=R8) :: dnndz
c     Output variables: Incident neutral flux has been broken up in two parts:
c     A constant contribution fnnic and a contribution varying linearly with the 
c     neutral density fnniv. fmomni: momentum flux of 
c     incident neutrals.
      real (kind=R8),intent(out):: fnni, fnnic,fnniv,fmomni,nni,nnwwni,
     &                             feneni
c     Internal variables
      real (kind=R8):: dnndyc,dnndyv,vcx,vion,
     2    vtot,fnni_nodrifts,pit,pit2,mn,nnf,ux,uy,uz,
     3    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     4    fmox,fmoy,fmoz,fmomnip,fmomnit,pn,expu,t1,
     5    dnndy,t_av,mfp,c_th,upf, v_theta, v_r, v_phi
      integer :: ic,ix2
*   ..statement functions
      external subini, subend, expu

      call subini ('CalcIncidentFluxes')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCn,isn+1) !!AFN_TODO correct species index
      ux = upf*pit2
      uy = -sina*isign*upf*pit
      uz = -cosa*isign*upf*pit 

      if (print_output.ne.0) then
        write(*,*) '-------------------------------------------'
        write(*,*) 'iCv = ', iCv
        write(*,*) 'iCn = ', iCn
        write(*,*) 'pit = ', pit
        write(*,*) 'ux = ', ux
        write(*,*) 'uy = ', uy
        write(*,*) 'uz = ', uz
        write(*,*) '-------------------------------------------'
      endif

      dnndyc = isign*dnndy2
      dnndyv = 0
      dnndy = dnndyc
c     Calculate rate coefficients at boundary
      vcx=0.0_R8
      do ic = 0, pl%ns-1
        t_av = (tif+tnf)/2.0_R8
        vcx = vcx + pl%na(iCv,ic) * expu(rt%rlcx(iCv,0,ic,iscx0) +
     &   rt%rlcx(iCv,1,ic,iscx0)*log(t_av/(am(isn)*ev)))
      enddo
      vion=expu(
     2 rt%rlsa(iCv,0,isn)+
     3 rt%rlsa(iCv,1,isn)*log(tef/ev))*nef
      vtot =vion + vcx
      c_th = sqrt(tnf/mn)
      mfp = c_th/vtot

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Incident particle flux density
      if (use_uy_uz_0.eq.0) then
      uy = 0.0_R8
      uz = 0.0_R8
      endif

      if (split_fluxes.eq.1) then
        fnnic = vcx/vtot**2*(dnndyc*uy*i1l+
     &     dnndzc*i2l)
        fnniv = vcx/vtot*(-i1l+
     &     1.0_R8/vtot*(dnndyv*uy*i1l+dnndzv*i2l))
        fnni = fnnic + fnniv * nnf
      else
        dnndz = -isign*dnndz2
        fnni = vcx/vtot**2*(dnndy*uy*i1l + dnndz*i2l) - 
     &         nnf*vcx/vtot*i1l
      endif

c     Assume poloidal and radial drifts are zero for momentum and energy
c     equations to reduce the dimensionality of the pre-evaluates integrals
c     of TRIM tables. Only toroidal drift remains!
      if (use_uy_uz_1.eq.0) then
      uy = 0.0_R8
      uz = 0.0_R8
      endif

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Incident particle flux density without poloidal and radial drift velocities
c     if split_fluxes = 1 the particle flux term is split into a constant part and a part that scales with nn in the guard cell but then
c     the gradient components dnndy and dnndz are simply the '5-point' version which may be erronous for steep targets.
c     if split_fluxes = 0 the particle flux term is not split-up, but dnndy and dnndz are then the true tangential and normal components.
      if (split_fluxes.eq.1) then
        fnni_nodrifts = vcx/vtot**2*(dnndyc*uy*i1l+
     2     dnndzc*i2l) + vcx/vtot*(-i1l+
     3     1.0_R8/vtot*(dnndyv*uy*i1l+dnndzv*i2l)) * nnf
      else
        fnni_nodrifts = vcx/vtot**2*(dnndy*uy*i1l + dnndz*i2l) - 
     &                  nnf*vcx/vtot*i1l
      endif

c     Incident momentum flux densities
      fmox= mn*vcx/vtot*(-nnf*ux*I1l + 
     2     1.0_R8/vtot*(dnndy*ux*uy*I1l + 
     3     dnndz*ux*I2l))


      fmoy= mn*vcx/vtot*(-nnf*uy*I1l + 
     2     1.0_R8/vtot*(dnndy*(tif/mn + uy**2)*I1l + 
     3     dnndz*uy*I2l))


      fmoz= mn*vcx/vtot*(-nnf*I2l + 
     2     1.0_R8/vtot*(dnndy*uy*I2l + 
     3     dnndz*I3l))


c       Pressure of incident neutrals (total including contribution from the drift velocity)
      pn= mn/3.0_R8*vcx/vtot*(nnf*((2.0_R8*tif/mn + ux**2 + 
     2     uy**2)*I0l + I2l) - 
     3     1.0_R8/vtot*(dnndy*((4.0_R8*tif/mn*uy + 
     4     ux**2*uy + uy**3)*I0l + uy*I2l) + 
     5     dnndz*((2.0_R8*tif/mn + ux**2 + uy**2)*I1l + 
     6     I3l)))


c       Density of incident neutrals
      nni = vcx/vtot*(nnf*I0l - 
     2     1.0_R8/vtot*(dnndy*uy*I0l + dnndz*I1l))

      nnwwni = nni*ux

c     Incident total energy flux density of neutrals
      feneni = vcx/vtot*(-nnf*((tif + mn/2.0_R8*(ux**2 + uy**2))*
     2     I1l + mn/2.0_R8*I3l) + 
     3     1.0_R8/vtot*(dnndy*((2.0_R8*tif + 
     4     mn/2.0_R8*(ux**2 + uy**2))*uy*I1l + mn/2.0_R8*uy*I3l) + 
     5     dnndz*((tif + mn/2.0_R8*(ux**2 + uy**2))*I2l + 
     6     mn/2.0_R8*I4l)))


c     Convert particle flux densities in fluxes
      fnnic = fnnic*area
      fnniv = fnniv*area
      fnni  = fnni*area

c     Calculate incident poloidal momentum flux of neutrals
      fmomnip = isign*(fmoz*cosa + fmoy*sina)* 
     2     area +isign*pn*area*cosa
      fmomnit = -fmox*area
c     calculate parallel momentum boundary flux
      fmomni = pit*fmomnip + pit2*fmomnit        
c     Convert energy flux densities in net energy flux contributions
      feneni = -feneni*area

c     If flux has incorrect sign, resort to Maxwellian distribution for BC's
      if (fnni.le.0.0_R8) then
!!        call CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
!!     &            area,pl,tif,geo,1,
!!     &            fnnic,fnniv,fmomni,
!!     &            nni,nnwwni,feneni)
        fnnic  = 0.0_R8
        fnniv  = 0.0_R8
        fnni   = 0.0_R8 
        fmomni = 0.0_R8
        nni    = 0.0_R8
        nnwwni = 0.0_R8
        feneni = 0.0_R8
      endif

      call subend ()
      return
      end subroutine CalcIncidentFluxes

      subroutine CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tf,geo,isrc, cosa, sina,
     &            fnni,fnnic,fnniv,fmomni,
     &            nni,nnwwni,feneni)
!nh   08.05.2018
c     This routine calculates the incident neutral particle fluxes based on 
c     a truncated Maxwellian for the incident neutrals. These particle fluxes
c     are necessary for the hybrid delta f model. It is assumed that only the
c     neutral parallel drift velocity is non-negligible (unf). This has to be
c     adapted. 
      use b2mod_b2cmpa
      use b2us_plasma
      use b2us_geo
      implicit none
c     Input variables
c     isn is the neutral species index
c     isrc is the sign of the source needed for the hybrid delta f method
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      integer,intent(in) :: iCv, iFc, isn, istra,
     &    isrc
      real (kind=R8),intent(in) :: area, tf, isign, cosa, sina
c     Output variables: Incident neutral flux has been broken up in two parts:
c     A constant contribution fnnic and a contribution varying linearly with the 
c     neutral density fnniv.
      real (kind=R8),intent(out):: fnni, fnnic,fnniv,fmomni,nni,nnwwni,
     &                             feneni
c     Internal variables
      real (kind=R8):: pit,pit2,mn,nnf,ux,uy,uz,
     1    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     2    fmox,fmoy,fmoz,pn,fmomnip,fmomnit,unf
      integer :: ic 
*   ..statement functions
      external subini, subend

      call subini ('CalcIncidentFluxesMaxwellian')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cnh   neutral velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and only neutral parallel
c     velocity component
      unf = pl%ua(iCv,isn) 
      ux = unf*pit2
      uy = -sina*isign*unf*pit
      uz = -cosa*isign*unf*pit 

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tf,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Incident particle flux density
      fnnic = 0.0_R8
      fnniv = 0.0_R8
      if (isrc.eq.1) then
      	fnniv = -i1l
      elseif (isrc.eq.-1) then
	fnniv = i1r
      else    
	call xerrab('Wrong value for isrc.')
      endif
      fnni = fnnic + nnf*fnniv

c     Incident momentum flux densities
      fmox = -mn*nnf*ux*I1l 
      fmoy = -mn*nnf*uy*I1l
      fmoz = -mn*nnf*I2l

c     Pressure of incident neutrals
      pn = mn/3.0_R8*(nnf*((2.0_R8*tf/mn + ux**2 + 
     2     uy**2)*I0l + I2l))

c     Density of incident neutrals
      nni = nnf*I0l 

      nnwwni = nni*ux

c     Incident total energy flux density of neutrals
      feneni = -nnf*((tf + mn/2.0_R8*(ux**2 + uy**2))*
     2     I1l + mn/2.0_R8*I3l)

c     Convert particle flux densities in fluxes
      fnnic = fnnic*area
      fnniv = fnniv*area
      fnni  = fnni*area

c     Calculate incident poloidal momentum flux of neutrals 
      fmomnip = isign*(fmoz*cosa + fmoy*sina)*
     2   area + isign*pn*area*cosa
      fmomnit = -fmox*area
      fmomni = pit*fmomnip + pit2*fmomnit   
     
c     Convert energy flux densities in net energy flux contributions
      feneni = -feneni*area

      call subend ()
      return
      end subroutine CalcIncidentFluxesMaxwellian






      subroutine CalcRecycledFluxes(iCv,nCi,iFc,
     &             isign,
     &            isn,iscx0,istra,iwall,
     &            Phi_app,area,recyc0,pl,tif,tef,pof,nif,nef,
     &            geo,t0, cosa, sina, use_uy_uz_0,
     &            fmomrec,nnrec,nnwwnrec,fenerec)
!mb   17.07.17
c     This routine calculates the recycled neutral particle,momentum and energy fluxes based 
c     on the diffusion approximation (see papers of Niels Horsten).
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist , only:E_fc, accel_ion
      use b2us_geo
      use b2us_plasma
      implicit none
c     nif = plasma density on face, isn neutral species index
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
!!      type (B2Rates), intent (in) :: rt
      integer, intent(in):: iCv, nCi, iFc,
     &     isn, iscx0, iwall,istra,
     &     use_uy_uz_0
      real (kind=R8),intent(in):: Phi_app, area,recyc0,
     2    pof,nif,nef,t0,isign, tif, tef, cosa, sina
c     Output variables: Incident neutral flux has been broken up in two parts:
c     A constant contribution fnnic and a contribution varying linearly with the 
c     neutral density fnniv. fmomni: momentum flux of 
c     incident neutrals.
      real (kind=R8),intent(out):: fmomrec,nnrec,nnwwnrec,fenerec
c     Internal variables
      real (kind=R8):: cs,Mj,
     &    shpottauj,ishpottau,pf,pcorf,pt,vt,
     &     pit,pit2,mn,nnf,ux,uz,
     &    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &    rc1i,rmx1i,rmy1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &    fmox,fmoy,fmoz,fmomrecp,fmomrect,pnc,fenec,expu,t1,fdni,
     &    intp_2dtable,intp_2dtable0,tifev,upf
      real (kind = R8):: tab(0:nM,0:nT)
      integer :: ic,is,it,ispt,ispt2(1),ilT,ilM,ip1,ip2
      external subini, subend,intp_2dtable,intp_2dtable0, expu

      call subini ('CalcRecycledFluxes')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCv,isn+1) !! AFN_TODO species index
      ux = upf*pit2
      uz = -cosa*isign*upf*pit
c       ion particle flux density
c     AFN_TODO compare both options
      fdni = t0/area

c     Assume poloidal and radial drifts are zero for momentum and energy
c     equations to reduce the dimensionality of the pre-evaluates integrals
c     of TRIM tables. Only toroidal drift remains!
      if (use_uy_uz_0.eq.0) then
      uz = 0.0_R8
      endif

c AFN_TODO uz should never be put to zero here ? 
c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)


c     Interpolate tables to find recycling values
c       Isothermal speed of sound and Mach number
      cs = sqrt(2.0_R8*tif/mn)
      Mj = abs(ux)/cs
      tifev = tif/eV
      is = recycle_index(isn) 
      it = target_index(iCv-nCi-1)
c       Shpottau value: sheath_potential devided by tau (Ti/Te)
      shpottauj = accel_ion(istra)*qe*(pof-phi_app)/tif
c       For now, neirest neigbour interpolation for shpottau
      ispt2 = minloc(abs(shpottau-shpottauj)) !mb: adjust for index starting from 0
c       Interpolate
      ispt = ispt2(1) -1 !mb: adjust for index starting from 0
      tab = rc1(:,:,ispt,is,it)
      rc1i = intp_2dtable(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      tab = rmx1(:,:,ispt,is,it)
      rc1i = max(rc1i,0.0_R8)! Make sure thy don't become negative when extrapolating
      rmx1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmx1i = max(rmx1i,0.0_R8)
      tab = rmz1(:,:,ispt,is,it)
      rmz1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmz1i = max(rmz1i,0.0_R8)
      tab = rp11(:,:,ispt,is,it)
      rp11i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp11i = max(rp11i,0.0_R8)
      tab = rp21(:,:,ispt,is,it)
      rp21i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp21i = max(rp21i,0.0_R8)
      tab = rp31(:,:,ispt,is,it)
      rp31i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp31i = max(rp31i,0.0_R8)
      tab = re1(:,:,ispt,is,it)
      re1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      re1i = max(re1i,0.0_R8)
c     Calculate fast and thermal recycled parts
      pf = min(rc1i*cs/(-i1l),recyc0)+1.0E-30_R8
      pcorf = min(pf,recyc0)/pf 
      pt = recyc0-pf+1.0E-30_R8
      vt = sqrt(2.0_R8*E_fc*eV/mn)
c     The reflection coefficients below already account for the probability of
c     fast particle reflection. Like in Eirene, we will only pump fast particles
c     once all thermal particles are pumped.


c     Recycled momentum flux densities
      fmox = pcorf*mn*rmx1i*sign(1.0_R8,ux)*cs**2*fdni/(-i1l)

      fmoy = 0.0_R8

      fmoz = pcorf*mn*rmz1i*cs**2*fdni/(-i1l)
    
      fmoz = fmoz + 2.0_R8/3.0_R8*mn*pt*fdni*vt
c       Pressure of recycled neutrals 
      pnc = pcorf*mn/3.0_R8*rp11i*cs**2*fdni/(-i1l)

      pnc = pnc + 2.0_R8/3.0_R8*mn*pt*fdni*vt
c       Density of recycled neutrals
      nnrec = pcorf*rp21i*fdni/(-i1l)
      nnrec = nnrec+pt*fdni*2.0_R8/vt
c       product of density and toroidal velocity
      nnwwnrec = pcorf*rp31i*cs*sign(1.0_R8,ux)*fdni/(-i1l)

c     recycled total energy flux density of neutrals
      fenec = 0.5_R8*pcorf*mn*re1i*cs**3*fdni/(-i1l)
      fenec = fenec + 0.5_R8*mn*pt*fdni*vt**2


c     Calculate recycled poloidal momentum flux of neutrals
      fmomrecp = -isign*(fmoz*cosa + fmoy*sina)*
     2              area +isign*pnc*area*cosa
      fmomrect = fmox*area
      fmomrec = pit*fmomrecp + pit2*fmomrect         
c     Convert energy flux densities in energy fluxes
      fenerec = fenec*area
      call subend ()
      return
      end subroutine CalcRecycledFluxes

      subroutine CalcReflectedFluxes(iCv,iCn,nCi,iFc,isign,
     &            isn,iscx0,istra,iwall,
     &            Phi_app,area,recyc0,pl,tif,tnf,tef,pof,                           
     &            nef,geo,rt,dnndy2,
     &            dnndz,
     &            cosa, sina, use_uy_uz_0,
     &            fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
!mb   17.07.17
c     This routine calculates the reflected neutral particle,momentum and energy fluxes with incident neutral flux based on the diffusion approximation (see papers of Niels Horsten).
      use b2mod_b2cmpa
      use b2mod_rates
      use b2mod_neutrals_namelist , only:E_fc
      use b2us_geo
      use b2us_plasma
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2Rates), intent (in) :: rt

      integer,intent(in):: iCv, iCn, nCi, iFc,
     &    isn, iscx0, iwall,istra, use_uy_uz_0
      real (kind=R8),intent(in):: Phi_app, area,recyc0,
     2    pof,nef,isign,tif,tnf,tef,
     &    sina,cosa, dnndy2
      real (kind=R8) :: dnndz
c     Output variables: Incident neutral flux has been broken up in two parts:
c     A constant contribution fnnic and a contribution varying linearly with the 
c     neutral density fnniv. fmomni: momentum flux of 
c     incident neutrals.
      real (kind=R8),intent(out):: fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &                             fene_el
c     Internal variables
      real (kind=R8):: cs,Mj,
     &    shpottauj,ishpottau,pf,pcorf,pt,vt,vion,vcx,vtot,
     &     pit,pit2,mn,nnf,ux,uy,uz,
     &    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &    rc1i,rmx1i,rmy1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &    rc2i,rmx2i,rmy2i,rmz2i,rp12i,rp22i,rp32i,re2i,
     &    fnni,fnni_nodrifts,fmox,fmoy,fmoz,fmomreflp,fmomreflt,pn,fene,
     &   dnndy,expu,t1,intp_2dtable,intp_2dtable0,tifev,t_av,
     &   upf, nif
      real (kind = R8):: tab(0:nM,0:nT)
      integer :: ic,is,it,ispt,ispt2(1),ilT,ilM,ip1,ip2,ix2
      external subini, subend,intp_2dtable,intp_2dtable0, expu


      call subini ('CalcReflectedFluxes')
c     Initialize/calculate some variables
      nif = pl%na(iCv,isn+1) !! AFN_TODO species index
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCn,isn+1) !! AFN_TODO species index
                             !! We use iCn instead of iCv because EIRENE doesn't see the value in the guard cell and the difference can be large.
      ux = upf*pit2
      uy = -sina*isign*upf*pit
      uz = -cosa*isign*upf*pit 

      dnndz = -isign*dnndz

c     Calculate tangental and perpendicular gradients of neutral density at boundary
      dnndy = isign*dnndy2

c     Calculate rate coefficients at boundary
      vcx=0.0_R8
      do ic = 0, pl%ns-1
        t_av = (tif+tnf)/2.0_R8
        vcx = vcx + pl%na(iCv,ic) * expu(rt%rlcx(iCv,0,ic,iscx0) +
     &   rt%rlcx(iCv,1,ic,iscx0)*log(t_av/(am(isn)*ev)))
      enddo
      vion=expu(
     2 rt%rlsa(iCv,0,isn)+
     3 rt%rlsa(iCv,1,isn)*log(tef/ev))*nef
      vtot =vion + vcx

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Calculate incident particle flux density
      fnni = vcx/vtot*(-i1l*nnf+(dnndz*i2l+dnndy*uy*i1l)/vtot)

c     Assume poloidal and radial drifts are zero for momentum and energy
c     equations to reduce the dimensionality of the pre-evaluates integrals
c     of TRIM tables. Only toroidal drift remains!
      if (use_uy_uz_0.eq.0) then
        uy = 0.0_R8
        uz = 0.0_R8
      endif

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)


c     Interpolate tables to find recycling values
c       Isothermal speed of sound and Mach number
      cs = sqrt(2.0_R8*tif/mn)
      Mj = abs(ux)/cs
      tifev = tif/eV
      is = recycle_index(isn)
      it = target_index(iCv-nCi-1)
      ispt = 0 ! Neutrals are not accelerated to sheath
c       Interpolate
      tab = rc1(:,:,ispt,is,it)
      rc1i = intp_2dtable(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rc1i = max(rc1i,0.0_R8) ! Make sure thy don't become negative when extrapolating
      tab = rc2(:,:,ispt,is,it)
      rc2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rc2i = max(rc2i,0.0_R8)
      tab = rmx1(:,:,ispt,is,it)
      rmx1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmx1i = max(rmx1i,0.0_R8)
      tab = rmx2(:,:,ispt,is,it)
      rmx2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmx2i = max(rmx2i,0.0_R8)
c      rmy1i = intp_2dtable0(rmy1(:,:,ispt,is,it),M,T,dM,dT,nM,nT,
c     & Mj,tifev,ilM,ilT,ip1,ip2)
      tab = rmy2(:,:,ispt,is,it)
      rmy2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmy2i = max(rmy2i,0.0_R8)
      tab = rmz1(:,:,ispt,is,it)
      rmz1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmz1i = max(rmz1i,0.0_R8)
      tab = rmz2(:,:,ispt,is,it)
      rmz2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rmz2i = max(rmz2i,0.0_R8)
      tab = rp11(:,:,ispt,is,it)
      rp11i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp11i = max(rp11i,0.0_R8)
      tab = rp21(:,:,ispt,is,it)
      rp21i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp21i = max(rp21i,0.0_R8)
      tab = rp31(:,:,ispt,is,it)
      rp31i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp31i = max(rp31i,0.0_R8)
      tab = rp12(:,:,ispt,is,it)
      rp12i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp12i = max(rp12i,0.0_R8)
      tab = rp22(:,:,ispt,is,it)
      rp22i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp22i = max(rp22i,0.0_R8)
      tab = rp32(:,:,ispt,is,it)
      rp32i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rp32i = max(rp32i,0.0_R8)
      tab = re1(:,:,ispt,is,it)
      re1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      re1i = max(re1i,0.0_R8)
      tab = re2(:,:,ispt,is,it)
      re2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      re2i = max(re2i,0.0_R8)
c     Calculate incident particle flux density
      fnni_nodrifts = vcx/vtot*(-i1l*nnf+dnndz*i2l/vtot)
c     Calculate fast and thermal reflected parts
      pf = min(vcx/vtot*(nnf*rc1i*cs+dnndz*rc2i*cs**2/vtot)
     &   /fnni_nodrifts,recyc0)+1.0E-30_R8
      pcorf = min(pf,recyc0)/pf 
      
      pt = recyc0-pf+1.0E-30_R8
      vt = sqrt(2.0_R8*E_fc*eV/mn)
c     The reflection coefficients below already account for the probability of
c     fast particle reflection. Like in Eirene, we will only pump fast particles
c     once all thermal particles are pumped.

c     reflected momentum flux densities
      fmox = pcorf*mn*vcx/vtot*(nnf*rmx1i*cs**2 + 
     2     1.0_R8/vtot*dnndz*rmx2i*cs**3)*sign(1.0_R8,ux)
      fmoy = -pcorf*mn*vcx/vtot**2*dnndy*rmy2i*cs**3

      fmoz = pcorf*mn*vcx/vtot*(nnf*rmz1i*cs**2 + 
     2     1.0_R8/vtot*dnndz*rmz2i*cs**3)
    
      fmoz = fmoz + pt*2.0_R8/3.0_R8*mn*fnni_nodrifts*vt
c       Pressure of reflected neutrals 
      pn = pcorf*mn/3.0_R8*vcx/vtot*(nnf*rp11i*cs**2+
     2     1.0_R8/vtot*dnndz*rp12i*cs**3)

      pn = pn + pt*2.0_R8/3.0_R8*mn*fnni_nodrifts*vt
c       Density of reflected neutrals
      nnrefl = pcorf*vcx/vtot*(nnf*rp21i+1.0_R8/vtot*dnndz*rp22i*cs)
      nnrefl = nnrefl+pt*fnni_nodrifts*2.0_R8/vt
c       product of density and toroidal velocity
      nnwwnrefl = pcorf*vcx/vtot*(nnf*rp31i*cs+
     2     1.0_R8/vtot*dnndz*rp32i*cs**2)*sign(1.0_R8,ux)
c     reflected total energy flux density of neutrals
      fene = 0.5_R8*pcorf*mn*vcx/vtot*(nnf*re1i*cs**3+
     2     1.0_R8/vtot*dnndz*re2i*cs**4)
      fene = fene + pt*fnni_nodrifts*E_fc*ev
      fene_el =   - pt*fnni_nodrifts*E_fc*ev    !! the Franck-Condon energy must at least be provided by the electrons.



c     Calculate reflected poloidal momentum flux of neutrals
      fmomreflp = -isign*(fmoz*cosa + fmoy*sina)* 
     2   area +isign*pn*area*cosa
      fmomreflt = fmox*area
      fmomrefl = pit*fmomreflp + pit2*fmomreflt         
c     Convert energy flux densities in energy fluxes
      fenerefl = fene*area
      fene_el =  fene_el*area

!!    resort to maxwellian distribution of incident flux has wrong sign
      if (fnni.le.0.0_R8) then
!!        call CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
!!     &            isn,istra,iwall,area,recyc0,pl,tif,                          
!!     &            geo,sign_cosa, sign_sina, use_uy_uz,
!!     &            fmomrefl,nnrefl,nnwwnrefl,fenerefl)
             fmomrefl = 0.0_R8
             nnrefl = 0.0_R8
             nnwwnrefl = 0.0_R8
             fenerefl = 0.0_R8
      endif
      call subend ()
      return
      end subroutine CalcReflectedFluxes

      subroutine CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,pl,tf,                         
     &            geo, cosa, sina,
     &            fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
!nh   11.06.18
c     This routine calculates the reflected neutral particle, momentum and energy fluxes with incident neutral flux based on a truncated Maxwellian
c     for the incident neutrals.
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist , only:E_fc
      use b2us_plasma
      use b2us_geo
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      integer,intent(in):: iCv, nCi, iFc, isn, iwall, istra
      real (kind=R8),intent(in):: area, recyc0, isign, tf, sina, cosa
c     Output variables
      real (kind=R8),intent(out):: fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &                             fene_el
c     Internal variables
      real (kind=R8):: cs,Mj,pf,pcorf,pt,vt,
     &     pit,pit2,mn,nnf,ux,uz,
     &     i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &     rc1i,rmx1i,rmy1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &     rc2i,rmx2i,rmy2i,rmz2i,rp12i,rp22i,rp32i,re2i,
     &     fnni,fmox,fmoy,fmoz,fmomreflp,fmomreflt,pn,fene,
     &     intp_2dtable,intp_2dtable0,tfev,unf
      real (kind = R8):: tab(0:nM,0:nT)
      integer :: ic,is,it,ispt,ispt2(1),ilT,ilM,ip1,ip2
*   ..statement functions
      external subini, subend, intp_2dtable, intp_2dtable0

      call subini ('CalcReflectedFluxesMaxwellian')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
c     Neutral velocity components: toroidal (x),tangential (y), and normal (z)
      unf = pl%ua(iCv,isn)
      ux =  unf*pit2
      uz = -cosa*isign*unf*pit
      
c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tf,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Calculate incident particle flux density
      fnni = -i1l*nnf

c     Interpolate tables to find recycling values
c     Isothermal speed of sound and Mach number
      cs = sqrt(2.0_R8*tf/mn)
      Mj = abs(ux)/cs
      tfev = tf/eV
      is = recycle_index(isn)
      it = target_index(iCv-nCi-1)
      ispt = 0 ! Neutrals are not accelerated to sheath
c       Interpolate
      tab = rc1(:,:,ispt,is,it)
      rc1i = intp_2dtable(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rc1i = max(rc1i,0.0_R8) ! Make sure they don't become negative when extrapolating
      tab = rc2(:,:,ispt,is,it)
      rc2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rc2i = max(rc2i,0.0_R8)
      tab = rmx1(:,:,ispt,is,it)
      rmx1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rmx1i = max(rmx1i,0.0_R8)
      tab = rmx2(:,:,ispt,is,it)
      rmx2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rmx2i = max(rmx2i,0.0_R8)
c      rmy1i = intp_2dtable0(rmy1(:,:,ispt,is,it),M,T,dM,dT,nM,nT,
c     & Mj,tfev,ilM,ilT,ip1,ip2)
      tab = rmy2(:,:,ispt,is,it)
      rmy2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rmy2i = max(rmy2i,0.0_R8)
      tab = rmz1(:,:,ispt,is,it)
      rmz1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rmz1i = max(rmz1i,0.0_R8)
      tab = rmz2(:,:,ispt,is,it)
      rmz2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rmz2i = max(rmz2i,0.0_R8)
      tab = rp11(:,:,ispt,is,it)
      rp11i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp11i = max(rp11i,0.0_R8)
      tab = rp21(:,:,ispt,is,it)
      rp21i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp21i = max(rp21i,0.0_R8)
      tab = rp31(:,:,ispt,is,it)
      rp31i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp31i = max(rp31i,0.0_R8)
      tab = rp12(:,:,ispt,is,it)
      rp12i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp12i = max(rp12i,0.0_R8)
      tab = rp22(:,:,ispt,is,it)
      rp22i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp22i = max(rp22i,0.0_R8)
      tab = rp32(:,:,ispt,is,it)
      rp32i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rp32i = max(rp32i,0.0_R8)
      tab = re1(:,:,ispt,is,it)
      re1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      re1i = max(re1i,0.0_R8)
      tab = re2(:,:,ispt,is,it)
      re2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      re2i = max(re2i,0.0_R8)


c     Calculate fast and thermal reflected parts
      pf = min(nnf*rc1i/fnni,recyc0)+1.0E-30_R8
      pcorf = min(pf,recyc0)/pf 


      pt = recyc0-pf+1.0E-30_R8
      vt = sqrt(2.0_R8*E_fc*eV/mn)
c     The reflection coefficients below already account for the probability of
c     fast particle reflection. Like in Eirene, we will only pump fast particles
c     once all thermal particles are pumped.

c     reflected momentum flux densities
      fmox = pcorf*mn*nnf*rmx1i*cs**2*sign(1.0_R8,ux)
      fmoy = 0._R8
      fmoz = pcorf*mn*nnf*rmz1i*cs**2
      fmoz = fmoz + pt*2.0_R8/3.0_R8*mn*fnni*vt
c     Pressure of reflected neutrals 
      pn = pcorf*mn/3.0_R8*nnf*rp11i*cs**2
      pn = pn + pt*2.0_R8/3.0_R8*mn*fnni*vt
c     Density of reflected neutrals
      nnrefl = pcorf*nnf*rp21i
      nnrefl = nnrefl+pt*fnni*2.0_R8/vt
c     Product of density and toroidal velocity
      nnwwnrefl = pcorf*nnf*rp31i*cs*sign(1.0_R8,ux)
c     reflected total energy flux density of neutrals
      fene = 0.5_R8*pcorf*mn*nnf*re1i*cs**3
      fene = fene + pt*fnni*E_fc*ev
      fene_el =   - pt*fnni*E_fc*ev  !! the Franck-Condon energy is provided by the electrons



c     Calculate reflected poloidal momentum flux of neutrals
      fmomreflp = -isign*(fmoz*cosa + fmoy*sina)* 
     2   area +isign*pn*area*cosa
      fmomreflt = fmox*area
      fmomrefl = pit*fmomreflp + pit2*fmomreflt         
c     Convert energy flux densities in energy fluxes
      fenerefl = fene*area
      fene_el = fene_el*area
      call subend ()
      return
      end subroutine CalcReflectedFluxesMaxwellian







      subroutine CalcMomentsMaxwellian(u,T,m,i0l,i0r,i1l,i1r,
     2     i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)
c mb This subroutine calculates the moments of a 1D Maxwellian distribution (the
c zeroth until the sixth moment)
c u = drifting velocity
c T = temperature
c m = mass of a particle

      implicit none
      !input variables
      real (kind=R8),intent(in)::u,T,m
      !output variables
      real (kind=R8),intent(out)::i0l,i0r,
     &    i1l,i1r,i2l,i2r
      real (kind=R8),intent(out)::
     &    i3l,i3r,i4l,i4r,
     &    i5l,i5r,i6l,i6r  
      !internal variables
      real (kind=R8),save::uold=0.0_R8,told=0.0_R8,mold=0.0_R8
      real (kind=R8),save::int0l,int0r,
     &    int1l,int1r,int2l,int2r,
     &    int3l,int3r,int4l,int4r,
     &    int5l,int5r,int6l,int6r 
      real (kind=R8):: vth,arg
c      real (kind=R8):: S15AEC
      intrinsic erf
      external subini, subend
c      external S15AEC !mb NAG erf function
      call subini('CalcMomentsMaxwellian')
      
      if (uold.eq.u.and.Told.eq.t.and.mold.eq.m) then
        !Don't do the same thing twice
        i0l = int0l
        i0r = int0r
        i1l = int1l
        i1r = int1r
        i2l = int2l
        i2r = int2r
        i3l = int3l
        i3r = int3r
        i4l = int4l
        i4r = int4r
        i5l = int5l
        i5r = int5r
        i6l = int6l
        i6r = int6r
        call subend ()
        return
      endif
      vth = sqrt(2.0_R8*T/m)
      arg = -(u/vth)
c     int0l = int_{-inf}^{0} M(v) dv
      int0l = 0.5_R8*(1.0_R8+erf(arg)) 

c     int0r = int_{0}^{inf} M(v) dv
      int0r = 0.5_R8*(1.0_R8-erf(arg))

c     int1l = int_{-inf}^{0} v M(v) dv
      int1l = 0.5_R8*u*(1.0_R8+erf(arg))- 
     2   0.5_R8/sqrt(pi)*vth*exp(arg**2)

c     int1r = int_{0}^{inf} v M(v) dv
      int1r = 0.5_R8*u*(1.0_R8-erf(arg))+ 
     2   0.5_R8/sqrt(pi)*vth*exp(arg**2)

c     int2l = int_{-inf}^{0} v^{2} M(v) dv
      int2l = (0.5_R8*T/m+0.5_R8*u**2)*(1.0_R8+erf(arg))-
     2   0.5_R8/sqrt(pi)*vth*u*exp(arg**2)

c     int2r = int_{0}^{inf} v^{2} M(v) dv
      int2r = (0.5_R8*T/m+0.5_R8*u**2)*(1.0_R8-erf(arg))+
     2   0.5_R8/sqrt(pi)*vth*u*exp(arg**2) 

c Implementation as optional arguments gave segmentation errors...
c     int3l = int_{-inf}^{0} v^{3} M(v) dv
       int3l = (1.5_R8*T/m*u+0.5_R8*u**3)*(1.0_R8+erf(arg)) - 
     2   1.0_R8/sqrt(pi)*(0.5_R8*vth**3+0.5_R8*vth*u**2)*exp(arg**2) 

c     int3r = int_{0}^{inf} v^{3} M(v) dv
       int3r = (1.5_R8*T/m*u+0.5_R8*u**3)*(1-erf(arg)) + 
     2   1.0_R8/sqrt(pi)*(0.5_R8*vth**3+0.5_R8*vth*u**2)*exp(arg**2) 

c     int4l = int_{-inf}^{0} v^{4} M(v) dv
       int4l = (1.5_R8*(T/m)**2+3.0_R8*T/m*u**2+u**4/2.0_R8)*(1.0_R8+
     2    erf(arg)) - 1.0_R8/sqrt(pi)*
     3    (1.25_R8*vth**3*u+0.5_R8*vth*u**3)*exp(arg**2) 

c     int4r = int_{0}^{inf} v^{4} M(v) dv
       int4r = (1.5_R8*(T/m)**2+3.0_R8*T/m*u**2+u**4/2.0_R8)*(1.0_R8-
     2    erf(arg)) + 1.0_R8/sqrt(pi)*
     3    (1.25_R8*vth**3*u+0.5_R8*vth*u**3)*exp(arg**2)

c     int5l = int_{-inf}^{0} v^{5} M(v) dv
       int5l = (7.5_R8*(T/m)**2*u+5.0_R8*T/m*u**3+u**5/2.0_R8)*(1.0_R8+
     2    erf(arg)) -1.0_R8/sqrt(pi)*
     3    (vth**5+2.25_R8*vth**3*u**2+0.5_R8*vth*u**4)*exp(arg**2) 

c     int5r = int_{0}^{inf} v^{5} M(v) dv
       int5r = (7.5_R8*(T/m)**2*u+5.0_R8*T/m*u**3+u**5/2.0_R8)*(1.0_R8-
     2    erf(arg)) +1.0_R8/sqrt(pi)*
     3    (vth**5+2.25_R8*vth**3*u**2+0.5_R8*vth*u**4)*exp(arg**2)

c     int6l = int_{-inf}^{0} v^{6} M(v) dv
       int6l = (7.5_R8*(T/m)**3+22.5_R8*(T/m)**2*u**2+7.5_R8*T/m*u**4+
     2    0.5_R8*u**6)*(1.0_R8+erf(arg))-
     3    1.0_R8/sqrt(pi)*(4.125_R8*vth**5*u+3.5_R8*vth**3*u**3 + 
     4    0.5_R8*vth*u**5)*exp(arg**2) 

c     int6r = int_{0}^{inf} v^{6} M(v) dv
       int6r = (7.5_R8*(T/m)**3+22.5_R8*(T/m)**2*u**2+7.5_R8*T/m*u**4+
     2    0.5_R8*u**6)*(1.0_R8-erf(arg))+
     3    1.0_R8/sqrt(pi)*(4.125_R8*vth**5*u+3.5_R8*vth**3*u**3 + 
     4    0.5_R8*vth*u**5)*exp(arg**2) 
        i0l = int0l
        i0r = int0r
        i1l = int1l
        i1r = int1r
        i2l = int2l
        i2r = int2r
        i3l = int3l
        i3r = int3r
        i4l = int4l
        i4r = int4r
        i5l = int5l
        i5r = int5r
        i6l = int6l
        i6r = int6r
      call subend ()
      return
      end subroutine CalcMomentsMaxwellian

      subroutine dealloc_recycle_fnn
      implicit none
      if(.not.recycle_initialised) return
      deallocate(trimrecyc_present)
      deallocate(target_index)
      deallocate(recycle_index)
      deallocate(M,T,shpottau,rc1,rc2,rmx1,rmx2,rmy1,rmy2,rmz1,rmz2)
      deallocate(M_dummy,T_dummy,shpottau_dummy)  !nh 21.02.2020
      deallocate(re1,re2,rp11,rp12,rp21,rp22,rp31,rp32)
      recycle_initialised=.false.
      return
      end subroutine dealloc_recycle_fnn




      end module b2mod_recycle

