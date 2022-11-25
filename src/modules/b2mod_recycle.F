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
c     recycling coefficients. To turn on the model, the flag b2stbr_recycle_afn is set to '1'
c     Mixed materials are not considered in this recycling model

      logical, save :: recycle_initialised
      data recycle_initialised /.false./

      real (kind=R8),save::uold=0.0_R8,told=0.0_R8,mold=0.0_R8
      integer, save :: nm, nt, nshpottau,nx1,ny1,ns1

c     Recycle_index and target_index refer for each species and boundary in which integer
c     the corresponding recycling coefficients are stored.
      integer, allocatable, save ::
     *  recycle_index(:),target_index(:)
c     Boolean to indicate whether trim data is present for species/material combination
      logical, allocatable, save ::
     *  trimrecyc_present(:,:)
      real (kind=R8), allocatable, save ::
     &  M(:),T(:),dM(:),dT(:),shpottau(:),dSh(:),
     &  rc1(:,:,:,:,:),rc2(:,:,:,:,:),
     &  rmx1(:,:,:,:,:),rmx2(:,:,:,:,:),rmy1(:,:,:,:,:),rmy2(:,:,:,:,:),
     &  rmz1(:,:,:,:,:),rmz2(:,:,:,:,:),re1(:,:,:,:,:),re2(:,:,:,:,:),
     &  rp11(:,:,:,:,:),rp12(:,:,:,:,:),rp21(:,:,:,:,:),
     &  rp22(:,:,:,:,:),rp31(:,:,:,:,:),rp32(:,:,:,:,:),
     &  M_dummy(:), T_dummy(:), shpottau_dummy(:),        !nh 21.02.2020
     &  fna_mol(:,:)                                    !nh 30.03.2021

      real (kind=R8),save::int0l,int0r,
     &    int1l,int1r,int2l,int2r,
     &    int3l,int3r,int4l,int4r,
     &    int5l,int5r,int6l,int6r

      contains

      subroutine alloc_b2mod_recycle(nFcd,nsd)
      implicit none
      integer,intent(in) :: nFcd,nsd

      allocate(fna_mol(nFcd,0:nsd-1))

      fna_mol = 0.0_R8

      return
      end subroutine alloc_b2mod_recycle

      subroutine read_recycle_fnn(nsd,nbnd_local)
      use b2mod_b2cmpa
      use b2mod_elements
      use b2mod_neutrals_namelist
      use b2mod_subsys
      implicit none
      integer,intent(in):: nsd, nbnd_local
      character*256 filename
      integer istra, isp, is, il,it,itm,ndata
      integer idum(0:9)
      logical open_file

      character*10::fnn_version
      character lblfnn*120
      character lblfnn_dummy*120              !nh 21.02.2020
      character*2 target_code(0:nbnd_local-1)

      character*2 recycle_code(0:nsd-1)
      external open_file
      external cfrure, cfruin, machsfr
      real (kind=R8) ::
     *  machsfr
      if(recycle_initialised) return
      call subini ('read_recycle_fnn')
      ns1 = nsd
      allocate(recycle_index(0:nsd-1))
      allocate(target_index(0:nbnd_local-1))
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
      do istra=1,nstrai
        do it=0,itm
          if(surf_mat(istra).eq.target_code(it)) then
            target_index(istra)=it
            goto 6
          endif
        enddo
        itm=itm+1
        target_index(istra)=itm
        target_code(itm)=surf_mat(istra)
    6   continue
        write(*,'(1x,i3,1x,i3,1x,a)') istra,target_index(istra),
     1   target_code(target_index(istra))
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
              allocate(dSh(0:nshpottau-1))
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
              dSh(:) = 1.0_R8/
     &               (shpottau(1:nshpottau)-shpottau(0:nshpottau-1))

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

              deallocate(M_dummy,T_dummy,shpottau_dummy)

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


      subroutine CalcIncidentFluxesDiffusion(iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,safeguard,
     &            fnni,fmomni,nni,nnwwni,feneni)
c     This routine calculates the incident neutral particle,momentum and energy fluxes based
c     on the diffusion approximation (see papers of Niels Horsten)
c main reference:
c   Horsten, Niels, Giovanni Samaey, and M. Baelmans.
c  "Development and assessment of 2D fluid neutral models that include atomic databases and a microscopic reflection model."
c   Nuclear Fusion 57.11 (2017): 116043.
      use b2mod_math
      use b2mod_rates
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2rates), intent (in) :: rt
      integer,intent(in):: iCv,iCn, iFc,
     &    isn,iscx0,isi,istra, use_uy_uz_0, use_uy_uz_1, safeguard
      real (kind=R8),intent(in):: Phi_app, area, tnf, tif, tef,
     2    pof,nef,isign,dnndy,dnndz, cosa, sina
c     Output variables:
c     fnni: incident neutral flux
c     fmomni: momentum flux of incident neutrals.
      real (kind=R8),intent(out):: fnni,fmomni,nni,nnwwni,
     &                             feneni
c     Internal variables
      real (kind=R8):: vcx,vion,
     2    vtot,fnni_nodrifts,pit,pit2,mn,nnf,ux,uy,uz,
     3    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     4    fmox,fmoy,fmoz,fmomnip,fmomnit,pn,
     5    t_av,mfp,c_th,upf
      integer :: ic

      call subini ('CalcIncidentFluxesDiffusion')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCn,isi)
      ux = upf*pit2
      uy = -sina*isign*upf*pit
      uz = -cosa*isign*upf*pit

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

      if (use_uy_uz_0.eq.0) then
        uy = 0.0_R8
        uz = 0.0_R8
      endif

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Incident particle flux density

      fnni = vcx/vtot**2*(dnndy*uy*i1l + dnndz*i2l) -
     &         nnf*vcx/vtot*i1l

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
      fnni_nodrifts = vcx/vtot**2*(dnndy*uy*i1l + dnndz*i2l) -
     &                  nnf*vcx/vtot*i1l

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
      fnni  = fnni*area

c     Calculate incident poloidal momentum flux of neutrals
      fmomnip = isign*(fmoz*cosa + fmoy*sina)*
     2     area +isign*pn*area*cosa
      fmomnit = -fmox*area
c     calculate parallel momentum boundary flux
      fmomni = pit*fmomnip + pit2*fmomnit
c     Convert energy flux densities in net energy flux contributions
      feneni = -feneni*area

c if incident flux is negative (this can happen if the diffusion approximation is used where it is not valid), use one-sided maxwellian anyway
      if (fnni.lt.0.0_R8 .and. safeguard.eq.1) then
        call CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tnf,geo, cosa, sina,
     &            fnni,fmomni,
     &            nni,nnwwni,feneni)
      endif

      call subend ()
      return
      end subroutine CalcIncidentFluxesDiffusion

      subroutine CalcIncidentFluxesKn(Kn_b1,Kn_b2,
     &            iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,L,
     &            fnni,fmomni,nni,nnwwni,feneni)
c     This routine is an interface to CalcIncidentFluxesMaxwellian and CalcIncidentFluxesDiffusion
c     It will decide which one to use based on the Knudsen number (Kn) of the neutral atoms.
      use b2mod_math
      use b2mod_rates
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2rates), intent (in) :: rt
      integer,intent(in):: iCv,iCn, iFc,
     &    isn,iscx0,isi,istra, use_uy_uz_0, use_uy_uz_1
      real (kind=R8),intent(in):: Phi_app, area, tnf, tif, tef,
     &    pof,nef,isign,dnndy,dnndz,cosa, sina,
     &    Kn_b1, Kn_b2,L
      real (kind=R8) :: t_av
c     Output variables
      real (kind=R8),intent(out):: fnni,fmomni,nni,nnwwni,
     &                             feneni
c     Internal variables
      real (kind=R8):: vcx,vion,vtot,mn,mfp,c_th,Kn
      integer :: ic
      real (kind=R8):: fnniM,nnwwniM,nniM,feneniM,fmomniM
      real (kind=R8):: fnniD,nnwwniD,nniD,feneniD,fmomniD

      call subini ('CalcIncidentFluxesKn')

c     Calculate rate coefficients at boundary
      mn = am(isn)*mp
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
      Kn = mfp/L !! L = macroscopic length scale

c     if Kn is low, we use the diffusion approximation
c     if Kn is larger, we use the drifting maxwellian approximation
c     if the Kn is in between, we use a linear interpolation
      if (Kn.lt.Kn_b1) then

        call CalcIncidentFluxesDiffusion(iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,1,!! <- '1': = safeguard for diffusion BCs
     &            fnni,fmomni,nni,nnwwni,feneni)

      elseif (Kn.gt.Kn_b2) then

        call CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tnf,geo, cosa, sina,
     &            fnni,fmomni,
     &            nni,nnwwni,feneni)

      else
        call CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tnf,geo, cosa, sina,
     &            fnniM,fmomniM,
     &            nniM,nnwwniM,feneniM)
        call CalcIncidentFluxesDiffusion(iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,1,!! <- '1': = safeguard for diffusion BCs
     &            fnniD,fmomniD,nniD,nnwwniD,feneniD)
        call CalcLinComb(Kn,Kn_b1,Kn_b2,fnniD,  fnniM,   fnni)
        call CalcLinComb(Kn,Kn_b1,Kn_b2,fmomniD,fmomniM, fmomni)
        call CalcLinComb(Kn,Kn_b1,Kn_b2,nniD,   nniM,    nni)
        call CalcLinComb(Kn,Kn_b1,Kn_b2,nnwwniD,nnwwniM, nnwwni)
        call CalcLinComb(Kn,Kn_b1,Kn_b2,feneniD,feneniM, feneni)
      endif

      call subend ()
      return
      end subroutine CalcIncidentFluxesKn

      subroutine CalcIncidentFluxesMinKn(iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,
     &            fnni,fmomni,nni,nnwwni,feneni)
c     This routine is an interface to either CalcIncidentFluxesMaxwellian or CalcIncidentFluxesDiffusion
c     It will decide which one to use based on the minimum Knudsen of the neutral atoms.
c     If the choice is based on Kn at each boundary face (as is done in CalcIncidentFluxesKn), additional convergence problems may arise.
c     In this function, the choice is the same for each boundary face of the same group (istra) and the choice is based on the lowest Kn of that "stratum".
c     Here, the choice is made based on maxw_eff (effective), which has been set in b2stbr_pĥys.F
      use b2mod_rates
      use b2mod_b2cmpa
      use b2us_geo
      use b2us_plasma
      use b2mod_neutrals_namelist , only: maxw_eff
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2rates), intent (in) :: rt
      integer,intent(in):: iCv,iCn, iFc,
     &    isn,iscx0,isi,istra, use_uy_uz_0, use_uy_uz_1
      real (kind=R8),intent(in):: Phi_app, area, tnf, tif, tef,
     &    pof,nef,isign,dnndy,dnndz,cosa, sina
c     Output variables
      real (kind=R8),intent(out):: fnni,fmomni,nni,nnwwni,
     &                             feneni

      call subini ('CalcIncidentFluxesMinKn')

      if (maxw_eff(istra).eq.0) then
        call CalcIncidentFluxesDiffusion(iCv,iCn,iFc,isign,
     &            isn,iscx0,isi,istra,
     &            Phi_app,area,pl,tif,tnf,tef,pof,nef,geo,
     &            rt,dnndy, dnndz,cosa,sina,
     &            use_uy_uz_0,use_uy_uz_1,1, !! <- '1': = safeguard for diffusion BCs
     &            fnni,fmomni,nni,nnwwni,feneni)
      else
        call CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tnf,geo, cosa, sina,
     &            fnni,fmomni,
     &            nni,nnwwni,feneni)
      endif

      call subend ()
      return
      end subroutine CalcIncidentFluxesMinKn

      subroutine CalcIncidentFluxesMaxwellian(iCv,iFc,isn,isign,istra,
     &            area,pl,tf,geo, cosa, sina,
     &            fnni,fmomni,
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
      use b2mod_subsys
      implicit none
c     Input variables
c     isn is the neutral species index
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      integer,intent(in) :: iCv, iFc, isn, istra
      real (kind=R8),intent(in) :: area, tf, isign, cosa, sina
c     Output variables
      real (kind=R8),intent(out):: fnni,fmomni,nni,nnwwni,
     &                             feneni
c     Internal variables
      real (kind=R8):: pit,pit2,mn,nnf,ux,uy,uz,
     1    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     2    fmox,fmoy,fmoz,pn,fmomnip,fmomnit,unf

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
      fnni = -nnf*i1l

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
     &            isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tef,pof,nif,nef,
     &            geo,mpg,t0, cosa, sina, use_uy_uz_0,
     &            fnnrec,fmomrec,nnrec,nnwwnrec,fenerec)
!mb   17.07.17
c     This routine calculates the recycled neutral particle, momentum and energy fluxes based
c     on the diffusion approximation (see papers of Niels Horsten).
      use b2mod_rates
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist , only:E_fc, accel_ion
      use b2us_geo
      use b2us_plasma
      use b2us_map
      use b2mod_subsys
      implicit none
c     nif = plasma density on face, isn neutral species index
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
!!      type (B2Rates), intent (in) :: rt
      integer, intent(in):: iCv, nCi, iFc,
     &     isn, iscx0, isi,iwall,istra,
     &     use_uy_uz_0
      real (kind=R8),intent(in):: Phi_app, area,recyc0,recycm,
     2    pof,nif,nef,t0,isign, tif, tef, cosa, sina,
     &    fluid_frac_hyb
!! fluid_frac_hyb: in fluid neutral cases or standard spatially hybrid cases this multiplier is always equal to 1.
!! only in the generalized form of the spatially hybrid approach, where the fluid and kinetic treatment are mixed
!! on a single target stratum, this value can deviate from 1.
c     Output variables:
c     fmomni: momentum flux of incident neutrals.
      real (kind=R8),intent(out):: fnnrec,fmomrec,nnrec,nnwwnrec,fenerec
c     Internal variables
      real (kind=R8):: cs,Mj,
     &    shpottauj,pf,pcorf,pt,vt,
     &    pit,pit2,mn,nnf,ux,uz,
     &    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &    rc1i,rmx1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &    fmox,fmoy,fmoz,fmomrecp,fmomrect,pnc,fenec,fdni,
     &    intp_3dtable,intp_3dtable0,intp_3dtable00,
     &    tifev,upf,
     &    fna_mol_rec
      real (kind = R8):: tab(0:nM,0:nT,0:nshpottau),ff(0:1,0:1,0:1)
      integer :: is,it,ilT,ilM,ilSh,ip1,ip2,ip3
      external intp_3dtable, intp_3dtable0, intp_3dtable00

      call subini ('CalcRecycledFluxes')
c     Initialize/calculate some variables
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCv,isi)
      ux = upf*pit2
      uz = -cosa*isign*upf*pit
c       ion particle flux density
      fdni = t0/area

c     Assume poloidal and radial drifts are zero for momentum and energy
c     equations to reduce the dimensionality of the pre-evaluates integrals
c     of TRIM tables. Only toroidal drift remains!
      if (use_uy_uz_0.eq.0) then
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
      it = target_index(istra)
c       Shpottau value: sheath_potential devided by tau (Ti/Te)
      shpottauj = accel_ion(istra)*qe*(pof-phi_app)/tif
      tab = rc1(:,:,:,is,it)
      rc1i = intp_3dtable(tab,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh,ip1,ip2,ip3)
      rc1i = max(rc1i,0.0_R8)! Make sure thy don't become negative when extrapolating

      ff = rmx1(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      rmx1i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      rmx1i = max(rmx1i,0.0_R8)

      ff = rmz1(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      rmz1i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      rmz1i = max(rmz1i,0.0_R8)

      ff = rp11(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      rp11i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      rp11i = max(rp11i,0.0_R8)

      ff = rp21(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      rp21i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      rp21i = max(rp21i,0.0_R8)

      ff = rp31(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      rp31i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      rp31i = max(rp31i,0.0_R8)

      ff = re1(ilM:ilM+1,ilT:ilT+1,ilSh:ilSh+1,is,it)
      re1i = intp_3dtable00(ff,M,T,shpottau,dM,dT,dSh,nM,nT,nshpottau,
     & Mj,tifev,shpottauj,ilM,ilT,ilSh)
      re1i = max(re1i,0.0_R8)

c     Calculate fast and thermal recycled parts
      pf = min(rc1i*cs/(-i1l),recyc0)+1.0E-30_R8
      pcorf = min(pf,recyc0)/pf
      pt = recyc0-pf+1.0E-30_R8

      vt = sqrt(2.0_R8*E_fc*eV/mn)
c     The reflection coefficients below already account for the probability of
c     fast particle reflection. Like in Eirene, we will only pump fast particles
c     once all thermal particles are pumped.

cnh 22.02.2020
c     Calculate the flux emitted as molecules. At the moment, molecules will be
c     treated kinetically.
      fna_mol_rec = recycm*pt*fdni*area*0.5_R8
      fna_mol(iFc,isn) = fna_mol(iFc,isn) + fna_mol_rec

      fnnrec = max(0.0_R8,(recyc0*t0 - 2.0_R8*fna_mol_rec)*
     1  fluid_frac_hyb)


c     Recycled momentum flux densities
      fmox = pcorf*mn*rmx1i*sign(1.0_R8,ux)*cs**2*fdni/(-i1l)*
     1  fluid_frac_hyb

      fmoy = 0.0_R8

      fmoz = pcorf*mn*rmz1i*cs**2*fdni/(-i1l)*fluid_frac_hyb

      fmoz = fmoz + 2.0_R8/3.0_R8*mn*(1-recycm)*pt*fdni*vt*
     1  fluid_frac_hyb
c       Pressure of recycled neutrals
      pnc = pcorf*mn/3.0_R8*rp11i*cs**2*fdni/(-i1l)*fluid_frac_hyb

      pnc = pnc + 2.0_R8/3.0_R8*mn*(1-recycm)*pt*fdni*vt*fluid_frac_hyb
c       Density of recycled neutrals
      nnrec = pcorf*rp21i*fdni/(-i1l)*fluid_frac_hyb
      nnrec = nnrec+(1-recycm)*pt*fdni*2.0_R8/vt*fluid_frac_hyb
c       product of density and toroidal velocity
      nnwwnrec = pcorf*rp31i*cs*sign(1.0_R8,ux)*fdni/(-i1l)*
     1  fluid_frac_hyb

c     recycled total energy flux density of neutrals
      fenec = 0.5_R8*pcorf*mn*re1i*cs**3*fdni/(-i1l)*fluid_frac_hyb
      fenec = fenec + 0.5_R8*mn*(1-recycm)*pt*fdni*vt**2*fluid_frac_hyb


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

      subroutine CalcReflectedFluxesDiffusion(iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tnf,tef,pof,
     &            nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0,use_uy_uz_1,safeguard,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
!mb   17.07.17
c     This routine calculates the reflected neutral particle,momentum and energy fluxes with incident neutral flux based on the diffusion approximation (see papers of Niels Horsten).
      use b2mod_b2cmpa
      use b2mod_math
      use b2mod_rates
      use b2mod_neutrals_namelist , only:E_fc
      use b2us_geo
      use b2us_plasma
      use b2us_map
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2Rates), intent (in) :: rt
      type (mapping), intent (in) :: mpg
      integer,intent(in):: iCv, iCn, nCi, iFc,
     &    isn, iscx0, isi, iwall,istra, use_uy_uz_0, use_uy_uz_1,
     &    safeguard
      real (kind=R8),intent(in):: Phi_app, area,recyc0,recycm,
     2    pof,nef,isign,tif,tnf,tef,
     &    sina,cosa, dnndy, dnndz,
     &    fluid_frac_hyb
!! fluid_frac_hyb: in fluid neutral cases or standard spatially hybrid cases this multiplier is always equal to 1.
!! only in the generalized form of the spatially hybrid approach, where the fluid and kinetic treatment are mixed
!! on a single target stratum, this value can deviate from 1.
c     Output variables:
c     fmomni: momentum flux of
c     incident neutrals.
      real (kind=R8),intent(out):: fnnrefl,fmomrefl,nnrefl,
     &                             nnwwnrefl,fenerefl,fene_el
c     Internal variables
      real (kind=R8):: cs,Mj,
     &    pf,pcorf,pt,vt,vion,vcx,vtot,
     &     pit,pit2,mn,nnf,ux,uy,uz,
     &    i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &    rc1i,rmx1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &    rc2i,rmx2i,rmy2i,rmz2i,rp12i,rp22i,rp32i,re2i,
     &    fnni,fnni_nodrifts,fmox,fmoy,fmoz,fmomreflp,fmomreflt,pn,fene,
     &    intp_2dtable,intp_2dtable0,tifev,t_av,
     &    upf, nif,mfp,c_th,
     &    fna_mol_refl
      real (kind = R8):: tab(0:nM,0:nT)
      integer :: ic,is,it,ispt,ilT,ilM,ip1,ip2
      external intp_2dtable,intp_2dtable0

      call subini ('CalcReflectedFluxesDiffusion')
c     Initialize/calculate some variables
      nif = pl%na(iCv,isi)
      nnf = pl%na(iCv,isn)
      mn = am(isn)*mp
      pit = abs(geo%cvbb(iCv,0)/geo%cvbb(iCv,3))
      pit2 = abs(geo%cvbb(iCv,2)/geo%cvbb(iCv,3))
cmb   plasma velocities toroidal (x),tangential (y), and normal (z) to boundary assuming absence of 9-point stencil and drifts
      upf = pl%ua(iCn,isi)!! We use iCn instead of iCv because EIRENE doesn't see the value in the guard cell and the difference can be large.
      ux = upf*pit2
      uy = -sina*isign*upf*pit
      uz = -cosa*isign*upf*pit

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

      if (use_uy_uz_0.eq.0) then
        uy = 0.0_R8
        uz = 0.0_R8
      endif

c     Calculate the moments of a truncated drifting Maxwellian distribution
      call CalcMomentsMaxwellian(uz,tif,mn,i0l,i0r,i1l,i1r,i2l,i2r,
     2     i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r)

c     Calculate incident particle flux density
      fnni = vcx/vtot*(-i1l*nnf+(dnndz*i2l+dnndy*uy*i1l)/vtot)

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


c     Interpolate tables to find recycling values
c       Isothermal speed of sound and Mach number
      cs = sqrt(2.0_R8*tif/mn)
      Mj = abs(ux)/cs
      tifev = tif/eV
      is = recycle_index(isn)
      it = target_index(istra)
      ispt = 0 ! Neutrals are not accelerated to sheath
c       Interpolate
      tab = rc1(:,:,ispt,is,it)
      rc1i = intp_2dtable(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT,ip1,ip2)
      rc1i = max(rc1i,0.0_R8) ! Make sure thy don't become negative when extrapolating
      tab = rc2(:,:,ispt,is,it)
      rc2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rc2i = max(rc2i,0.0_R8)
      tab = rmx1(:,:,ispt,is,it)
      rmx1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rmx1i = max(rmx1i,0.0_R8)
      tab = rmx2(:,:,ispt,is,it)
      rmx2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rmx2i = max(rmx2i,0.0_R8)
      tab = rmy2(:,:,ispt,is,it)
      rmy2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rmy2i = max(rmy2i,0.0_R8)
      tab = rmz1(:,:,ispt,is,it)
      rmz1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rmz1i = max(rmz1i,0.0_R8)
      tab = rmz2(:,:,ispt,is,it)
      rmz2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rmz2i = max(rmz2i,0.0_R8)
      tab = rp11(:,:,ispt,is,it)
      rp11i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp11i = max(rp11i,0.0_R8)
      tab = rp21(:,:,ispt,is,it)
      rp21i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp21i = max(rp21i,0.0_R8)
      tab = rp31(:,:,ispt,is,it)
      rp31i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp31i = max(rp31i,0.0_R8)
      tab = rp12(:,:,ispt,is,it)
      rp12i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp12i = max(rp12i,0.0_R8)
      tab = rp22(:,:,ispt,is,it)
      rp22i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp22i = max(rp22i,0.0_R8)
      tab = rp32(:,:,ispt,is,it)
      rp32i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      rp32i = max(rp32i,0.0_R8)
      tab = re1(:,:,ispt,is,it)
      re1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
      re1i = max(re1i,0.0_R8)
      tab = re2(:,:,ispt,is,it)
      re2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tifev,ilM,ilT)
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

cnh 22.02.2020
c     Calculate the part of the reflected flux emitted as molecules. At the moment,
c     the molecules are treated kinetically.
      if (fnni.ge.0.0_R8) then
        fna_mol_refl = recycm*pt*fnni_nodrifts*area*0.5_R8
        fna_mol(iFc,isn) = fna_mol(iFc,isn) + fna_mol_refl

        fnnrefl = max(0.0_R8,(recyc0*fnni*area - 2.0_R8*fna_mol_refl)*
     1    fluid_frac_hyb)
      else
        fnnrefl = 0.0_R8
      endif


c     reflected momentum flux densities
      fmox = pcorf*mn*vcx/vtot*(nnf*rmx1i*cs**2 +
     2     1.0_R8/vtot*dnndz*rmx2i*cs**3)*sign(1.0_R8,ux)*fluid_frac_hyb
      fmoy = -pcorf*mn*vcx/vtot**2*dnndy*rmy2i*cs**3*fluid_frac_hyb

      fmoz = pcorf*mn*vcx/vtot*(nnf*rmz1i*cs**2 +
     2     1.0_R8/vtot*dnndz*rmz2i*cs**3)*fluid_frac_hyb

      fmoz = fmoz + (1-recycm)*pt*2.0_R8/3.0_R8*mn*fnni_nodrifts*vt*
     1  fluid_frac_hyb
c       Pressure of reflected neutrals
      pn = pcorf*mn/3.0_R8*vcx/vtot*(nnf*rp11i*cs**2+
     2     1.0_R8/vtot*dnndz*rp12i*cs**3)*fluid_frac_hyb

      pn = pn + (1-recycm)*pt*2.0_R8/3.0_R8*mn*fnni_nodrifts*vt*
     1  fluid_frac_hyb
c       Density of reflected neutrals
      nnrefl = pcorf*vcx/vtot*(nnf*rp21i+1.0_R8/vtot*dnndz*rp22i*cs)*
     1  fluid_frac_hyb
      nnrefl = nnrefl+(1-recycm)*pt*fnni_nodrifts*2.0_R8/vt*
     1  fluid_frac_hyb
c       product of density and toroidal velocity
      nnwwnrefl = pcorf*vcx/vtot*(nnf*rp31i*cs+
     2     1.0_R8/vtot*dnndz*rp32i*cs**2)*sign(1.0_R8,ux)*fluid_frac_hyb
c     reflected total energy flux density of neutrals
      fene = 0.5_R8*pcorf*mn*vcx/vtot*(nnf*re1i*cs**3+
     2     1.0_R8/vtot*dnndz*re2i*cs**4)*fluid_frac_hyb
      fene = fene + (1-recycm)*pt*fnni_nodrifts*E_fc*ev*fluid_frac_hyb
      fene_el =   - (1-recycm)*pt*fnni_nodrifts*E_fc*ev*fluid_frac_hyb    !! the Franck-Condon energy must at least be provided by the electrons.



c     Calculate reflected poloidal momentum flux of neutrals
      fmomreflp = -isign*(fmoz*cosa + fmoy*sina)*
     2   area +isign*pn*area*cosa
      fmomreflt = fmox*area
      fmomrefl = pit*fmomreflp + pit2*fmomreflt
c     Convert energy flux densities in energy fluxes
      fenerefl = fene*area
      fene_el =  fene_el*area

c if incident flux is negative (this can happen if the diffusion approximation is used where it is not valid), use one-sided maxwellian anyway
      if (fnni.lt.0.0_R8 .and. safeguard.eq.1) then
        call CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tnf, geo, mpg,cosa, sina,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
      endif



      call subend ()
      return
      end subroutine CalcReflectedFluxesDiffusion

      subroutine CalcReflectedFluxesKn(Kn_b1,Kn_b2,
     &             iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tnf,tef,pof,
     &            nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0, use_uy_uz_1, L,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
c     This routine is an interface to CalcReflectedFluxesMaxwellian and CalcReflectedFluxesDiffusion
c     It will decide which one to use based on the Knudsen number (Kn) of the neutral atoms.
      use b2mod_b2cmpa
      use b2mod_math
      use b2mod_rates
      use b2us_geo
      use b2us_plasma
      use b2us_map
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2Rates), intent (in) :: rt
      type (mapping), intent (in) :: mpg
      integer,intent(in):: iCv, iCn, nCi, iFc,
     &    isn, iscx0, isi, iwall,istra, use_uy_uz_0, use_uy_uz_1
      real (kind=R8),intent(in):: Phi_app, area,recyc0,recycm,
     2    pof,nef,isign,tif,tnf,tef,L,
     &    sina,cosa, dnndy, dnndz, Kn_b1 , Kn_b2, fluid_frac_hyb
c     Output variables:
      real (kind=R8),intent(out):: fnnrefl,fmomrefl,nnrefl,
     &                             nnwwnrefl,fenerefl,fene_el
c     Internal variables
      real (kind=R8)::
     &    vion,vcx,vtot,
     &    mn,t_av,mfp,c_th,Kn
      integer :: ic
      real (kind = R8):: fnnreflM,fmomreflM,nnreflM,nnwwnreflM,
     &    fenereflM,fene_elM
      real (kind = R8):: fnnreflD,fmomreflD,nnreflD,nnwwnreflD,
     &    fenereflD,fene_elD

      call subini ('CalcReflectedFluxesKn')

c     Calculate rate coefficients at boundary
      mn = am(isn)*mp
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
      Kn = mfp/L !! L = macroscopic length scale

!!    zero out if incident flux has wrong sign
      if (Kn.lt.Kn_b1) then
        call CalcReflectedFluxesDiffusion(iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tnf,tef,pof,
     &            nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0, use_uy_uz_1, 1,!! <- '1': = safeguard for diffusion BCs
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
      elseif (Kn.gt.Kn_b2) then
        call CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tnf,geo, mpg,cosa, sina,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,fene_el)
      else
        call CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,recycm,fluid_frac_hyb,pl,
     &            tnf,geo, mpg,cosa, sina,
     &            fnnreflM,fmomreflM,nnreflM,nnwwnreflM,fenereflM,
     &            fene_elM)

        call CalcReflectedFluxesDiffusion(iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,pl,tif,tnf,
     &            tef,pof,nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0, use_uy_uz_1, 1,!! <- '1': = safeguard for diffusion BCs
     &            fnnreflD,fmomreflD,nnreflD,nnwwnreflD,fenereflD,
     &            fene_elD)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,fnnreflD,fnnreflM,fnnrefl)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,fmomreflD,fmomreflM,
     &                     fmomrefl)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,nnreflD,nnreflM,nnrefl)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,nnwwnreflD,  nnwwnreflM,
     &                     nnwwnrefl)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,fenereflD,fenereflM,
     &                     fenerefl)
          call CalcLinComb(Kn,Kn_b1,Kn_b2,fene_elD,fene_elM,fene_el)
      endif

      call subend ()
      return
      end subroutine CalcReflectedFluxesKn

      subroutine CalcReflectedFluxesMinKn(iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tnf,tef,pof,
     &            nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0, use_uy_uz_1,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
c     This routine is an interface to either CalcReflectedFluxesMaxwellian or CalcReflectedFluxesDiffusion
c     It will decide which one to use based on the minimum Kn of the neutral atoms.
c     If the choice is based on the Kn at each boundary face (as is done in CalcIncidentFluxesKn), additional convergence problems may arise.
c     In this function, the choice is the same for each boundary face of the same group (istra) and the choice is based on the lowest Kn of that "stratum".
c     Here, the choice is made based on maxw_eff (effective), which has been set in b2stbr_pĥys.F
      use b2mod_b2cmpa
      use b2mod_rates
      use b2us_geo
      use b2us_plasma
      use b2us_map
      use b2mod_neutrals_namelist , only: maxw_eff
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (B2Rates), intent (in) :: rt
      type (mapping), intent (in) :: mpg
      integer,intent(in):: iCv, iCn, nCi, iFc,
     &    isn, iscx0, isi, iwall,istra, use_uy_uz_0, use_uy_uz_1
      real (kind=R8),intent(in):: Phi_app, area,recyc0,recycm,
     2    pof,nef,isign,tif,tnf,tef,
     &    sina,cosa, dnndy, dnndz, fluid_frac_hyb
c     Output variables:
      real (kind=R8),intent(out):: fnnrefl,fmomrefl,nnrefl,
     &                             nnwwnrefl,fenerefl,fene_el

      call subini ('CalcReflectedFluxesMinKn')

!!    zero out if incident flux has wrong sign
      if (maxw_eff(istra).eq.0) then
        call CalcReflectedFluxesDiffusion(iCv,iCn,nCi,iFc,
     &            isign,isn,iscx0,isi,istra,iwall,
     &            Phi_app,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tif,tnf,tef,pof,
     &            nef,geo,mpg,rt,dnndy,dnndz,
     &            cosa, sina, use_uy_uz_0, use_uy_uz_1, 1,!! <- '1': = safeguard for diffusion BCs
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
      else
        call CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,recycm,fluid_frac_hyb,
     &            pl,tnf,geo, mpg,cosa, sina,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
      endif

      call subend ()
      return
      end subroutine CalcReflectedFluxesMinKn


      subroutine CalcReflectedFluxesMaxwellian(iCv,nCi,iFc,isign,
     &            isn,istra,iwall,area,recyc0,recycm,fluid_frac_hyb,pl,
     &            tf,geo,mpg,cosa,sina,
     &            fnnrefl,fmomrefl,nnrefl,nnwwnrefl,fenerefl,
     &            fene_el)
!nh   11.06.18
c     This routine calculates the reflected neutral particle, momentum and energy fluxes with incident neutral flux based on a truncated Maxwellian
c     for the incident neutrals.
      use b2mod_b2cmpa
      use b2mod_neutrals_namelist , only:E_fc
      use b2us_plasma
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      integer,intent(in):: iCv, nCi, iFc, isn, iwall, istra
      real (kind=R8),intent(in):: area, recyc0, recycm, isign, tf, sina,
     &                            cosa, fluid_frac_hyb
c     Output variables
      real (kind=R8),intent(out):: fnnrefl,fmomrefl,nnrefl,
     &                             nnwwnrefl,fenerefl,fene_el
c     Internal variables
      real (kind=R8):: cs,Mj,pf,pcorf,pt,vt,
     &     pit,pit2,mn,nnf,ux,uz,
     &     i0l,i0r,i1l,i1r,i2l,i2r,i3l,i3r,i4l,i4r,i5l,i5r,i6l,i6r,
     &     rc1i,rmx1i,rmz1i,rp11i,rp21i,rp31i,re1i,
     &     rc2i,rmx2i,rmy2i,rmz2i,rp12i,rp22i,rp32i,re2i,
     &     fnni,fmox,fmoy,fmoz,fmomreflp,fmomreflt,pn,fene,
     &     intp_2dtable,intp_2dtable0,tfev,unf,
     &     fna_mol_refl
      real (kind = R8):: tab(0:nM,0:nT)
      integer :: is,it,ispt,ilT,ilM,ip1,ip2
*   ..statement functions
      external intp_2dtable, intp_2dtable0

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
      it = target_index(istra)
      ispt = 0 ! Neutrals are not accelerated to sheath
c       Interpolate
      tab = rc1(:,:,ispt,is,it)
      rc1i = intp_2dtable(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT,ip1,ip2)
      rc1i = max(rc1i,0.0_R8) ! Make sure they do not become negative when extrapolating
      tab = rc2(:,:,ispt,is,it)
      rc2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rc2i = max(rc2i,0.0_R8)
      tab = rmx1(:,:,ispt,is,it)
      rmx1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rmx1i = max(rmx1i,0.0_R8)
      tab = rmx2(:,:,ispt,is,it)
      rmx2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rmx2i = max(rmx2i,0.0_R8)
      tab = rmy2(:,:,ispt,is,it)
      rmy2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rmy2i = max(rmy2i,0.0_R8)
      tab = rmz1(:,:,ispt,is,it)
      rmz1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rmz1i = max(rmz1i,0.0_R8)
      tab = rmz2(:,:,ispt,is,it)
      rmz2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rmz2i = max(rmz2i,0.0_R8)
      tab = rp11(:,:,ispt,is,it)
      rp11i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp11i = max(rp11i,0.0_R8)
      tab = rp21(:,:,ispt,is,it)
      rp21i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp21i = max(rp21i,0.0_R8)
      tab = rp31(:,:,ispt,is,it)
      rp31i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp31i = max(rp31i,0.0_R8)
      tab = rp12(:,:,ispt,is,it)
      rp12i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp12i = max(rp12i,0.0_R8)
      tab = rp22(:,:,ispt,is,it)
      rp22i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp22i = max(rp22i,0.0_R8)
      tab = rp32(:,:,ispt,is,it)
      rp32i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      rp32i = max(rp32i,0.0_R8)
      tab = re1(:,:,ispt,is,it)
      re1i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      re1i = max(re1i,0.0_R8)
      tab = re2(:,:,ispt,is,it)
      re2i = intp_2dtable0(tab,M,T,dM,dT,nM,nT,
     & Mj,tfev,ilM,ilT)
      re2i = max(re2i,0.0_R8)


c     Calculate fast and thermal reflected parts
      pf = min(nnf*rc1i*cs/fnni,recyc0)+1.0E-30_R8
      pcorf = min(pf,recyc0)/pf
      pt = recyc0-pf+1.0E-30_R8


      vt = sqrt(2.0_R8*E_fc*eV/mn)
c     The reflection coefficients below already account for the probability of
c     fast particle reflection. Like in Eirene, we will only pump fast particles
c     once all thermal particles are pumped.

cnh 22.02.2020
c     Calculate the part of the reflected flux emitted as molecules. At the moment,
c     the molecules are treated kinetically.
      fna_mol_refl = recycm*pt*fnni*area*0.5_R8
      fna_mol(iFc,isn) = fna_mol(iFc,isn) + fna_mol_refl

      fnnrefl = max(0.0_R8,(recyc0*fnni*area - 2.0_R8*fna_mol_refl)*
     1  fluid_frac_hyb)


c     reflected momentum flux densities
      fmox = pcorf*mn*nnf*rmx1i*cs**2*sign(1.0_R8,ux)*fluid_frac_hyb
      fmoy = 0._R8
      fmoz = pcorf*mn*nnf*rmz1i*cs**2*fluid_frac_hyb
      fmoz = fmoz + (1-recycm)*pt*2.0_R8/3.0_R8*mn*fnni*vt*
     1   fluid_frac_hyb
c     Pressure of reflected neutrals
      pn = pcorf*mn/3.0_R8*nnf*rp11i*cs**2*fluid_frac_hyb
      pn = pn + (1-recycm)*pt*2.0_R8/3.0_R8*mn*fnni*vt*fluid_frac_hyb
c     Density of reflected neutrals
      nnrefl = pcorf*nnf*rp21i*fluid_frac_hyb
      nnrefl = nnrefl+(1-recycm)*pt*fnni*2.0_R8/vt*fluid_frac_hyb
c     Product of density and toroidal velocity
      nnwwnrefl = pcorf*nnf*rp31i*cs*sign(1.0_R8,ux)*fluid_frac_hyb
c     reflected total energy flux density of neutrals
      fene = 0.5_R8*pcorf*mn*nnf*re1i*cs**3*fluid_frac_hyb
      fene = fene + (1-recycm)*pt*fnni*E_fc*ev*fluid_frac_hyb
      fene_el =   - (1-recycm)*pt*fnni*E_fc*ev*fluid_frac_hyb  !! the Franck-Condon energy is provided by the electrons



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
      use b2mod_subsys
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
      real (kind=R8):: vth,arg
c      real (kind=R8):: S15AEC
      intrinsic erf
c      external S15AEC !mb NAG erf function
      call subini('CalcMomentsMaxwellian')

      if (uold.eq.u.and.Told.eq.t.and.mold.eq.m) then
        !Do not do the same thing twice
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

      subroutine CalcLinComb(Kn,Kn1,Kn2,F1,F2,F)
      implicit none
      !input variables
      real (kind=R8),intent(in)::Kn,Kn1,Kn2,F1,F2
      !output variables
      real (kind=R8),intent(out)::F
c                      F2
c              F       |
c       F1     |       |
c       |      |       |
c      Kn1----Kn------Kn2
      if (Kn.lt.Kn1.or.Kn.gt.Kn2) then
        call xerrab('Kn should be between Kn1 and Kn2')
      endif

      F = ((Kn-Kn1)*F2+(Kn2-Kn)*F1)/(Kn2-Kn1)

      return
      end subroutine CalcLinComb

      subroutine CalcMinKn(iCv,isn,istra,
     &            pl,dv,rt,L,min_Kn,Kn,mfp)
c     This routine calculates the minimum mean free path for a certain stratum
      use b2mod_rates
      use b2mod_math
      use b2mod_b2cmpa
      use b2us_plasma
      use b2mod_neutrals_namelist
      use b2mod_subsys
      implicit none
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2rates), intent (in) :: rt
      integer,intent(in):: iCv, istra, isn
      real (kind=R8),intent(in):: L
c     InOut variables
      real (kind=R8),intent(inout) :: min_Kn(nstrat)
c     Internal variables
      real (kind=R8):: vcx,vion,vtot,mn,mfp,c_th,t_av
      real (kind=R8) :: tn, ti, ne, te, Kn
      integer :: ic

      call subini ('CalcMinKn')

c     Calculate rate coefficients at boundary
      tn = pl%tn(iCv)
      ti = pl%ti(iCv)
      ne = dv%ne(iCv)
      te = pl%te(iCv)
      mn = am(isn)*mp
      vcx=0.0_R8
      do ic = 0, pl%ns-1
        t_av = (ti+tn)/2.0_R8
        vcx = vcx + pl%na(iCv,ic) * expu(rt%rlcx(iCv,0,ic,isn) +
     &   rt%rlcx(iCv,1,ic,isn)*log(t_av/(am(isn)*ev)))
      enddo
      vion=expu(
     2 rt%rlsa(iCv,0,isn)+
     3 rt%rlsa(iCv,1,isn)*log(te/ev))*ne
      vtot =vion + vcx
      c_th = sqrt(tn/mn)
      mfp = c_th/vtot
      Kn = mfp/L

      min_Kn(istra) = min(Kn,min_Kn(istra))


      call subend ()
      return
      end subroutine CalcMinKn

      subroutine dealloc_b2mod_recycle
      implicit none
      deallocate(fna_mol)
      return
      end subroutine dealloc_b2mod_recycle


      subroutine dealloc_recycle_fnn
      implicit none
      if(.not.recycle_initialised) return
      deallocate(trimrecyc_present)
      deallocate(target_index)
      deallocate(recycle_index)
      deallocate(M,T,shpottau,rc1,rc2,rmx1,rmx2,rmy1,rmy2,rmz1,rmz2)
      if (allocated(M_dummy)) then
        deallocate(M_dummy,T_dummy,shpottau_dummy)  !nh 21.02.2020
      endif
      deallocate(re1,re2,rp11,rp12,rp21,rp22,rp31,rp32)
      recycle_initialised=.false.
      return
      end subroutine dealloc_recycle_fnn

      subroutine precompute_kul_quant(nCv, nFc, nVx, switch, geo, mpg,
     &  pl, dv, co, rt, dnn, fluid_frac_hyb, kin_frac_hyb)
      use b2mod_types
      use b2mod_neutrals_namelist
      use b2us_geo
      use b2us_map
      use b2us_plasma
      use b2mod_switches
      implicit none
      integer, intent(in) :: nCv, nFc, nVx
      type (switches), intent (in) :: switch
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      type (B2Plasma), intent (in) :: pl
      type (B2Derivatives), intent (in) :: dv
      type (B2Coeff), intent (inout) :: co
      type (B2Rates), intent (in) :: rt
      real (kind=R8), intent(inout) :: dnn(nFc,0:1)
      real (kind=R8), intent (inout) :: fluid_frac_hyb(nFc),
     &                                  kin_frac_hyb(nFc)
c     local variables
      integer :: istra, i, iRc1, isn, iFc
      integer, save :: ncall = 0, style_visc = 0
      real(kind=R8) :: nnVx(nVx),
     &                 neVx(nVx),tiVx(nVx),
     &                 min_Kn(1:nstrat), mfp, mfp_out(nCv),
     &                 Kn,Kn_ti,Kn_ne,
     &                 Kn_out(nCv),Kn_ti_out(nCv),Kn_ne_out(nCv),
     &                 L_ne_out(nCv),L_ti_out(nCv),
     &                 dti(nFc,0:1),dne(nFc,0:1),L_ti,L_ne

      if (ncall.eq.0) then
        call ipgeti ('b2sihs_style_visc', style_visc)
      endif
      isn = 0 !! AFN_TODO generalize species index for multi-speciescases?
              !! For now KUL AFN models have only been developed for single species cases.
      call intvertex(nCv,nVx,mpg,geo%vxVol,pl%na(:,isn),nnVx)
      call intvertex(nCv,nVx,mpg,geo%vxVol,dv%ne,neVx)
      call intvertex(nCv,nVx,mpg,geo%vxVol,pl%ti,tiVx)
      call grad_nt(nCv, nFc, nVx, 1, geo, mpg, pl%na(:,isn), nnVx, dnn)
      call grad(nCv, nFc, nVx, 1, geo, mpg, pl%ti, tiVx, dti)
      call grad(nCv, nFc, nVx, 1, geo, mpg, dv%ne, neVx, dne)

      mfp_out = 0.0_R8
c     pre-computation of the minimal mean-free path along a boundary stratum
      min_Kn(:) = 1.0_R8
      do istra = 1, nstrai
*   ..loop over boundary cells
        do i = 1, mpg%rcCvP(istra,2)
          iRc1 = mpg%rcCv(mpg%rcCvP(istra,1)+i-1,1) ! the guard cell
          iFc  = mpg%cvFc(mpg%cvFcP(iRc1,1))        ! the guard face
          call CalcMinKn(iRc1,isn,istra,
     &        pl,dv,rt,switch%L_macro_AFN,min_Kn,Kn,mfp)
          mfp_out(iRc1) = mfp

cc compute Knudsen number
          L_ti = pl%ti(iRc1)/sqrt(dti(iFc,0)**2+dti(iFc,1)**2)  !! gradient length of Ti
          L_ne = dv%ne(iRc1)/sqrt(dne(iFc,0)**2+dne(iFc,1)**2)  !! gradient length of ne
c          L = min(L_ti,L_ne) !! shortest gradient length
          Kn = mfp/switch%L_macro_AFN
          Kn_ti = mfp/L_ti
          Kn_ne = mfp/L_ne
          Kn_out(iRc1) = Kn
          Kn_ti_out(iRc1) = Kn_ti
          Kn_ne_out(iRc1) = Kn_ne
          L_ti_out(iRc1) = L_ti
          L_ne_out(iRc1) = L_ne


cc compute fluid/kinetic fraction
          if (switch%use_auto_spatial_hyb.eq.1.and.
     &                                HYB_TYPE(istra).eq.'A') then !! 'A' = automatic
            if (Kn.le.switch%auto_spatial_hyb_Kn_1) then !! fluid treatment
              fluid_frac_hyb(iFc) = 1.0_R8
              kin_frac_hyb(iFc) = 0.0_R8
            elseif (Kn.ge.switch%auto_spatial_hyb_Kn_2) then  !! kinetic treatment
              fluid_frac_hyb(iFc) = 0.0_R8
              kin_frac_hyb(iFc) = 1.0_R8
            else
!! avoid discontinuities if a target face switches from fluid to kinetic. Linear interpolation in betweeen auto_Kn_1 and auto_Kn_2
              kin_frac_hyb(iFc) = (Kn-switch%auto_spatial_hyb_Kn_1)/
     &     (switch%auto_spatial_hyb_Kn_2-switch%auto_spatial_hyb_Kn_1)
              fluid_frac_hyb(iFc) = 1.0_R8 - kin_frac_hyb(iFc)
            endif
          else !! keep scaling factors both equal to 1.
            fluid_frac_hyb(iFc) = 1.0_R8
            kin_frac_hyb(iFc) = 1.0_R8
          endif
        enddo

cc  choose the type of neutral reflection boundary condition, based on the minimum Kn along that boundary (only used if maxw(istra) = 3).
        if (min_Kn(istra).lt.switch%Kn_b1) then
          maxw_eff(istra) = 0
        elseif (min_Kn(istra).gt.switch%Kn_b2) then
          maxw_eff(istra) = 1
        else
c         do nothing to prevent oscillations if minimum Kn is just around Kn_b1
        endif
      enddo
      if (switch%AFN_out.ne.0) then
        call my_out_us(70,nCv,0,mfp_out,'AFN_mfp')
        call my_out_us(70,nCv,0,Kn_out,'AFN_Kn')
        call my_out_us(70,nCv,0,Kn_ti_out,'AFN_Kn_ti')
        call my_out_us(70,nCv,0,Kn_ne_out,'AFN_Kn_ne')
        call my_out_us(70,nCv,0,L_ti_out,'L_ti')
        call my_out_us(70,nCv,0,L_ne_out,'L_ne')
      endif

      ncall = ncall + 1
      return
      end subroutine precompute_kul_quant


      end module b2mod_recycle

