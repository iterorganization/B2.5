      subroutine b2paux(mpg,gm,st,st_ext,ig)
c
c  version : 20.06.2005 02:25
c
c=======================================================================
c*** b2paux: post-processor for B2-Eirene
c***
c*** ig = -20 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  densities, temperatures for all species
c***
c*** ig = -19 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, drift velocity)
c***  in these points - for J. Brooks
c***
c*** ig = -18 :
c***
c***  Prepares output for DIVIMP run.
c***
c*** ig = -17 :
c***
c***  Prints detailed particle balance for the cells around the lower
c***  x-point
c***
c*** ig = -16 :
c***
c***  Produces a file with the plasma densities, temperatures,
c***  fluid velocities in the cell centers - for K. Krieger
c***
c*** ig = -15 :
c***
c***  Produces data files for the neutral fluxes on the specified
c***  surfaces. File names are "ndt_*.dat"
c***
c*** ig = -14 :
c***
c***  Produces data files for the cell-face-integrated fluxes of ions
c***  on every target. File names are "fi_tg_*.dat"
c***
c*** ig = -13 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the radiation intensity (total and components)
c***  in these points - for bolometry optimisation
c***
c*** ig = -12 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (densities, temperatures, fluid velocities)
c***  in the Eastern divertor - for H. Wuerz
c***
c*** ig = -11 :
c***
c***  Creates tables with coordinates of the grid cell corners along
c***  the separatrices and along the magnetic surfaces, together with
c***  the wall geometry in specified zones - for separate plotting
c***
c*** ig = -10 :
c***
c***  Creates table of fluxes and temperatures at the grid edges
c***
c*** ig = -9 :
c***
c***  Prints out the total number of nuclei of each kind in the system
c***  in file nuclei.tot
c***
c*** ig = -8 :
c***
c***  Creates tables of plasma pressure and neutral momentum source
c***  in files pr.dat and sm.dat
c***
c*** ig = -7 :
c***
c***  Prints out the neutral data read from fort.44 for debugging
c***
c*** ig = -6 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density)
c***  in these points - for V. Zaveryaev
c***
c*** ig = -5 :
c***
c***  Writes the CII sources in the plasma for DIVIMP
c***
c*** ig = -4 :
c***
c***  Produces a map of model radiation emissivity on a rectangular grid
c***  for 3D calculation of the wall loading by Romano
c***
c*** ig = -3 :
c***
c***  Prints out the particle flux integrals onto the targets
c***
c*** ig = -2 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density,
c***  and so on) in these points - for L. Horton
c***
c*** ig = -1 :
c***
c***  Produces a file with coordinates of the grid cell centers and
c***  the plasma parameters (density, temperatures, neutral density)
c***  in these points - for V. Abramov
c***
c*** ig > 999 :
c***
c***  Adds a line to the file ../sum_tmp/spec_pwr.dat with the power
c***  entering both divertors, power radiated along the divertor legs,
c***  power radiated near the targets, areas of the separatrix branches
c***  and mean specific power radiated along the divertor legs per
c***  square meter. ig-1000 cells are treated as "close to the targets"
c***
c*** ig >= DEF_NYD :
c***
c***  Produces two files with estimated helium concentration and
c***  neutral pressure in the PFR for summary plots
c***
c*** ig >= 0 :
c***
c***  Produces a file with plasma parameters (density, temperature,
c***  parallel flow velocity, parallel heat flux) against the
c***  distance along the magnetic field along flux surface index
c***  ig - for Oleg
c***
c*** Is to be called from b2plot!
c***
c*** Output files:
c***
c***  brooks.dat with the data for J. Brooks;
c***  fraction.he and pressure.ntr for summary plots;
c***  oleg.dat with the data for Oleg
c=======================================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2plot
      use b2mod_b2cmpa
      use b2mod_constants
      use b2mod_b2plot_ank
      use b2mod_user_namelist
      use b2mod_neutrals_namelist
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      use b2us_map
      use b2us_geo
      use b2us_plasma
      use b2mod_dimensions
      implicit none
      type(mapping), intent(in) :: mpg
      type(geometry), intent(in) :: gm
      type(B2State), intent(in) :: st
      type(B2StateExt), intent(in) :: st_ext
      integer, intent(in) :: ig
c
      integer ih(3)
      real (kind=R8), allocatable ::
     . hlp_pi(:), hlp_po(:), hlp_at(:,:), hlp_ml(:,:),
     . hlp_na(:,:), hlp_nm(:,:), mach(:,:), pz(:)
      real (kind=R8) :: hlp_a, hlp_m, hlp_n, hlp_p, hlp_s(3), hlp_ms(3),
     . mass, vel, vth, hlp_sp(2), hlp_mh(3)
      character*2, allocatable :: hlp_al(:)
      character*5 :: trg_lbl
      character*11, allocatable :: hlp_sg(:)
      character hlp_frmr*600
      integer, allocatable :: mz_hlp(:)
      integer lun
      parameter (lun=96)
      integer indsrf(stacklen)
      integer i, j, im, in, is, iCv, iFc, k, n, nindsrf
      logical found
      character frm*128, filename*256, extralabel*256
c=======================================================================
c
      if(ig.eq.-1) then
        call b2pabrmv(lun)
      else if(ig.eq.-2) then
        call b2phortn(lun,mpg,gm,st)
      else if(ig.eq.-3) then
c-----------------------------------------------------------------------
        allocate(hlp_pi(nspecies),hlp_po(nspecies))
        allocate(hlp_at(natmi,2),hlp_na(natmi,2))
        if (nmoli.gt.0) allocate(hlp_ml(nmoli,2),hlp_nm(nmoli,2))
        hlp_sp(1:2)=0.0_R8
        do k=1,nspecies
          hlp_pi(k)=0.0_R8
          hlp_po(k)=0.0_R8
        end do
        do k=1,natmi
          hlp_at(k,1:2)=0.0_R8
          hlp_na(k,1:2)=0.0_R8
        end do
        do k=1,nmoli
          hlp_ml(k,1:2)=0.0_R8
          hlp_nm(k,1:2)=0.0_R8
        end do
        do is=0,st%ns-1
          k=b2espcr(is)
          do i=1,mpg%divFcP(1,2)
            iFc=mpg%divFc(mpg%divFcP(1,1)+i-1)
            hlp_pi(k)=hlp_pi(k)+mpg%divFcOr(mpg%divFcP(1,1)+i-1)*
     .       (st%dv%fna(iFc,0,is)+st%dv%fna(iFc,1,is))
          end do
          do i=1,mpg%divFcP(maxval(mpg%strDiv),2)
            iFc=mpg%divFc(mpg%divFcP(maxval(mpg%strDiv),1)+i-1)
            hlp_po(k)=hlp_po(k)+
     .        mpg%divFcOr(mpg%divFcP(maxval(mpg%strDiv),1)+i-1)*
     .       (st%dv%fna(iFc,0,is)+st%dv%fna(iFc,1,is))
          end do
        end do
        n=nspecies+2*(natmi+nmoli)+1
        open(lun,file='particle.flx')
        write(lun,'(2a)')
     ,    '# Total fluxes of neutrals and plasma onto, ',
     ,    'and of neutrals and sputtered material from the targets'
        write(lun,'(2a)') '# ',trim(label)
        write(frm,720) n
  720   format('(a2,4x,',i3,'(1x,a16))')
        if (nmoli.gt.0) then
          write(lun, frm) '# ',
     ,     (trim(textsp(i))//'_pls_in ',i=0,nspecies-1),
     ,     (trim(textan(i))//'_atm_in ',i=0,natmi-1),
     ,     (trim(textmn(i))//'_mol_in ',i=0,nmoli-1),
     ,     (trim(textan(i))//'_atm_out',i=0,natmi-1),
     ,     (trim(textmn(i))//'_mol_out',i=0,nmoli-1),
     ,      'sput.'
        else
          write(lun, frm) '# ',
     ,     (trim(textsp(i))//'_pls_in ',i=0,nspecies-1),
     ,     (trim(textan(i))//'_atm_in ',i=0,natmi-1),
     ,     (trim(textan(i))//'_atm_out',i=0,natmi-1),
     ,      'sput.'
        end if
        write(frm,730) n
  730   format('(1p,1x,a5,',i3,'(5x,e12.5))')
        trg_lbl = 'inner'
        if (.not.gm%LSN) trg_lbl = 'outer'
        do k=1,ntns
          if (ktns(k).ne.1) cycle
          do i=1,natmi
            hlp_na(i,1)=hlp_na(i,1)+wldna(nlim+ltns(k),i,0)
            hlp_at(i,1)=hlp_at(i,1)+wldra(nlim+ltns(k),i,0)+
     .                              wldpa(nlim+ltns(k),i,0)
          end do
          do i=1,nmoli
            hlp_nm(i,1)=hlp_nm(i,1)+wldnm(nlim+ltns(k),i,0)
            hlp_ml(i,1)=hlp_ml(i,1)+wldrm(nlim+ltns(k),i,0)+
     .                              wldpm(nlim+ltns(k),i,0)
          end do
          hlp_sp(1)=hlp_sp(1)+wldspt(nlim+ltns(k),0)
        end do
        if (nmoli.gt.0) then
          write(lun, frm) trg_lbl,
     ,      (hlp_pi(i),i=1,nspecies),
     ,      (hlp_na(i,1),i=1,natmi),
     ,      (hlp_nm(i,1),i=1,nmoli),
     ,      (hlp_at(i,1),i=1,natmi),
     ,      (hlp_ml(i,1),i=1,nmoli),
     ,       hlp_sp(1)
        else
          write(lun, frm) trg_lbl,
     ,      (hlp_pi(i),i=1,nspecies),
     ,      (hlp_na(i,1),i=1,natmi),
     ,      (hlp_at(i,1),i=1,natmi),
     ,       hlp_sp(1)
        end if
        trg_lbl = 'outer'
        if (.not.gm%LSN) trg_lbl = 'inner'
        do k=1,ntns
          if (ktns(k).ne.2*max(1,mpg%nncut)) cycle
          do i=1,natmi
            hlp_na(i,2)=hlp_na(i,2)+wldna(nlim+ltns(k),i,0)
            hlp_at(i,2)=hlp_at(i,2)+wldra(nlim+ltns(k),i,0)+
     .                              wldpa(nlim+ltns(k),i,0)
          end do
          do i=1,nmoli
            hlp_nm(i,2)=hlp_nm(i,2)+wldnm(nlim+ltns(k),i,0)
            hlp_ml(i,2)=hlp_ml(i,2)+wldrm(nlim+ltns(k),i,0)+
     .                              wldpm(nlim+ltns(k),i,0)
          end do
          hlp_sp(2)=hlp_sp(2)+wldspt(nlim+ltns(k),0)
        end do
        if (nmoli.gt.0) then
          write(lun, frm) trg_lbl,
     ,      (hlp_po(i),i=1,nspecies),
     ,      (hlp_na(i,2),i=1,natmi),
     ,      (hlp_nm(i,2),i=1,nmoli),
     ,      (hlp_at(i,2),i=1,natmi),
     ,      (hlp_ml(i,2),i=1,nmoli),
     ,       hlp_sp(2)
        else
          write(lun, frm) trg_lbl,
     ,      (hlp_po(i),i=1,nspecies),
     ,      (hlp_na(i,2),i=1,natmi),
     ,      (hlp_at(i,2),i=1,natmi),
     ,       hlp_sp(2)
        end if
        deallocate(hlp_pi,hlp_po)
        deallocate(hlp_at,hlp_na)
        if (nmoli.gt.0) deallocate(hlp_ml,hlp_nm)
c-----------------------------------------------------------------------
      else if(ig.eq.-4) then
        call b2pplentd(lun)
      else if(ig.eq.-5) then
c-----------------------------------------------------------------------
c*** Write the CII sources in the plasma for DIVIMP
        found=.false.
        do is=0,st%ns-1
          if (found) cycle
          found = nint(zn(is)).eq.6 .and. nint(zamin(is)).eq.1
          if (found) k=is
        end do
        if(.not.found) then !{
          write(*,*) 'b2paux (-5): no CII species found'
          return
        end if !}
        open(lun,file='CII_source.dvi')
        write(lun,*) 'CII source from B2-Eirene [1/m**3/s]'
        write(lun,*) 'nCv,nncut: ',mpg%nCv,mpg%nncut
        write(lun,*) (st%srw%b2stbr_sna(iCv,k)/gm%cvVol(iCv),
     .                                                  iCv=1,mpg%nCv)
        close(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-6) then
        call b2pzavrv(lun,mpg,gm,st)
      else if(ig.eq.-7) then
c-----------------------------------------------------------------------
        print *
c*** print hydrogen isotopes fluxes across the surfaces

        hlp_s=0.0_R8
        in = 0
        ih = 0
        do is = 1, 3
          if (j_h_at(is).ne.0) then
            in = in+1
            ih(in) = is
          end if
        end do
        allocate(hlp_sg(5*in))
        hlp_mh=0.0_R8
        hlp_ms=0.0_R8
        i = 0
        do is = 1, 3
          if (j_h_at(is).eq.0) cycle
          i = i+1
          hlp_sg(5*(i-1)+1)=is_codes(eb2atcr(j_h_at(is)))//'atoms on'
          hlp_sg(5*(i-1)+2)=is_codes(eb2atcr(j_h_at(is)))//'atoms off'
          hlp_sg(5*(i-1)+3)=is_codes(eb2atcr(j_h_at(is)))//'mols on'
          hlp_sg(5*(i-1)+4)=is_codes(eb2atcr(j_h_at(is)))//'mols off'
          hlp_sg(5*(i-1)+5)=is_codes(eb2atcr(j_h_at(is)))//'saldo'
        end do
        write(hlp_frmr,'(a,i2,a)') '(1x,a7,2x,',5*in,'(1x,a11,2x))'
        write(6,hlp_frmr) 'surface',(hlp_sg(is),is=1,5*in)
        write(6,*) '  additional:'
        do i=1,nlimi+nstsi !{
          j=i
          if(i.gt.nlimi) then !{
            j=i-nlimi+nlim
            if(i.eq.nlimi+1) write(6,*) '  non-default:'
          end if !}
          do is = 1, 3 !{
            if (j_h_at(is).eq.0) cycle
            hlp_s(is)=wldna(j,j_h_at(is),0)-wldra(j,j_h_at(is),0)
     ,               -wldpa(j,j_h_at(is),0)
            do im = 1, nmoli !{
              if (l_h_mol(im,is).gt.0) then !{
                hlp_s(is)=hlp_s(is)+l_h_mol(im,is)*
     ,           (wldnm(j,im,0)-wldrm(j,im,0)-wldpm(j,im,0))
                hlp_mh(is)=l_h_mol(im,is)*wldnm(j,im,0)
                hlp_ms(is)=l_h_mol(im,is)*(wldrm(j,im,0)+wldpm(j,im,0))
              end if !}
            end do !}
          end do !}
          write(hlp_frmr,'(a,i2,a)') '(2x,i4,2x,1p,',5*in,'e14.6)'
          write(6,hlp_frmr)
     ,     j,(wldna(j,j_h_at(ih(is)),0),
     ,        wldra(j,j_h_at(ih(is)),0)+wldpa(j,j_h_at(ih(is)),0),
     ,        hlp_mh(ih(is)),hlp_ms(ih(is)),hlp_s(ih(is)),is=1,in)
        end do !}
        print *
        deallocate(hlp_sg)
c-----------------------------------------------------------------------
      else if(ig.eq.-8) then
c-----------------------------------------------------------------------
        allocate(mach(1:mpg%nCv,2))
        allocate(pz(1:mpg%nCv))
        call b2xppz(mpg%nCv,st%ns,st%dv%ne,st%pl%na,st%pl%te,st%pl%ti,
     .              pz,st_ext)
        do iCv=1,mpg%nCi
          mass=0.0_R8
          vel=0.0_R8
          do is=0,st%ns-1
            if (is_neutral(is)) cycle
            mass=mass+am(is)*mp*st%pl%na(iCv,is)
            vel=vel+st%rt%rza(iCv,is)*st%pl%na(iCv,is)*st%pl%ua(iCv,is)
          end do
          vth=sqrt(pz(iCv)/mass)
          mach(iCv,1)=vel/(vth*st%dv%ne(iCv))
          mach(iCv,1)=pz(iCv)*(1.0_R8+mach(iCv,1)**2)
        end do
        open(lun,file='pr.dat')
        write(lun,'(2a)') '#  Plasma pressure:  ',trim(label)
        write(frm,710) mpg%nCi
  710   format('((1p,2x,i6,1x,',i6,'(e11.4,2x)))')
        write(lun,frm) (iCv,mach(iCv,1),iCv=1,mpg%nCi)

        do iCv=1,mpg%nCi
          mach(iCv,2)=0.0_R8
        end do
        do iCv=1,mpg%nCi ! Only hydrogen atoms appear in momentum balance
          do is = 1, 3
            if (j_h_at(is).ne.0) then
              mach(iCv,2)=mach(iCv,2)+tab2(iCv,j_h_at(is),1)*
     .                                dab2(iCv,j_h_at(is),1)
            end if
          end do
        end do
        open(lun,file='prn.dat')
        write(lun,'(2a)') '#  Neutral pressure:  ',trim(label)
        write(frm,710) mpg%nCi
        write(lun,frm) (iCv,mach(iCv,2),iCv=1,mpg%nCi)
        open(lun,file='prt.dat')
        write(lun,'(2a)') '#  Total pressure:  ',trim(label)
        write(frm,710) mpg%nCi
        write(lun,frm) (iCv,mach(iCv,1)+mach(iCv,2),iCv=1,mpg%nCi)

        do iCv=1,mpg%nCv
          mach(iCv,2)=0.0_R8
        end do
        do iCv=1,mpg%nCv
          do k=0,st%ns-1
            mach(iCv,2)=mach(iCv,2)+st%srw%b2stbr_smo(iCv,k)
          end do
          mach(iCv,2)=mach(iCv,2)/gm%cvVol(iCv)
        end do
        open(lun,file='sm.dat')
        write(lun,'(2a)') '#  Neutral momentum source:  ',trim(label)
        write(frm,710) mpg%nCi
        write(lun,frm) (iCv,mach(iCv,2),j=1,mpg%nCi)

        close(lun)
        deallocate(mach)
        deallocate(pz)
c-----------------------------------------------------------------------
      else if(ig.eq.-9) then
c-----------------------------------------------------------------------
        open(lun,file='nuclei.tot')
        write(lun,'(2a)') '# Case: ',trim(label)
        write(lun,'(2a)') '# Total number of nuclei, including those ',
     ,                          'in plasma ions, atoms, and molecules:'
        do i = 1, nspecies
          hlp_p = 0.0_R8
          hlp_a = 0.0_R8
          hlp_m = 0.0_R8
          do is = 0, st%ns-1
            if (b2espcr(is).eq.i) then
              do iCv = 1, mpg%nCi
                hlp_p = hlp_p + gm%cvVol(iCv)*st%pl%na(iCv,is)
              end do
            end if
          end do
c
c*** Add the neutrals and molecules
c***  (for now, the test ions are neglected)
c
          do j=1,natmi
            if (latmscl(j).eq.i) then
              do iCv=1,mpg%nCi
                hlp_a = hlp_a + gm%cvVol(iCv)*dab2(iCv,j,1)
              end do
            end if
          end do
          do j=1,nmoli
            do k=1,natmi
              if (latmscl(k).eq.i) then
                if (mlcmp(k,j).gt.0) then
                  do iCv=1,mpg%nCi
                    hlp_m =
     .               hlp_m + gm%cvVol(iCv)*mlcmp(k,j)*dmb2(iCv,j,1)
                  end do
                end if
              end if
            end do
          end do
          hlp_n = hlp_a + hlp_m + hlp_p
          write(lun,'(a,1p,4e12.4)')
     .     species_list(i),hlp_n,hlp_p,hlp_a,hlp_m
        end do
        close(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-10) then
        call b2pedges(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-11) then
        call b2pgeom(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-12) then
        call b2pwuerz(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-13) then
        call b2pblmtr(lun,mpg,gm,st)
c-----------------------------------------------------------------------
      else if(ig.eq.-14) then
        call b2pgeier(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-15) then
        stack_req=(/2,0,0,0,0,0,0/)
        if(bad_args(stack_req)) then !{
          write(*,'(1x,a,i4,a)')
     ,      'b2paux: option',ig,' requires surface specification'
          write(*,'(1x,a,i4,a)')
     ,      'usage: ind1 ind2 ... indn nind',ig,' user'
          write(*,'(1x,2a)')
     ,      'where ind1,ind2,...,indn are the surface indices',
     ,      ' (like in EIRENE), and nind is the number of them'
          write(*,*)
          stack_req=(/1,0,0,0,0,0,0/)
          return
        end if !}
        k=ipop(istack,iistack)
        nindsrf=ipop(istack,iistack)
        stack_req=(/nindsrf,0,0,0,0,0,0/)
        if(bad_args(stack_req)) then !{
          write(*,*) 'b2paux: insufficient integer data on the stack'
          write(*,*) 'usage: ind1 ind2 ... indn nind',ig,' user'
          write(*,*) 'where ind1,ind2,...,indn are the surface indices',
     ,      ' (like in EIRENE), and nind is the number of them'
          call ipush(istack,iistack,nindsrf)
          call ipush(istack,iistack,k)
          stack_req=(/2,0,0,0,0,0,0/)
          return
        end if !}
        do i=1,nindsrf !{
          indsrf(i)=ipop(istack,iistack)
        end do !}
        call ipush(istack,iistack,k)
        call b2pntrdt(lun,indsrf,nindsrf)
c-----------------------------------------------------------------------
      else if(ig.eq.-16) then
        call b2krieger(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-17) then
        call analyse_xpt(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-18) then
        stack_req=(/0,0,1,1,0,0,0/)
        if (bad_args(stack_req)) return
        filename=spop(sstack,isstack)
        filename='../'//trim(filename)
        call get_extralabel(extralabel)
        open(unit=lun,file=filename)
        write(lun,'(2a)') '# ', trim(extralabel)
        write(lun,'(2a)') '# ', trim(graphlabel(iwrk))
        write(lun,'(2x,a3,6x,a1,11x,a1,11x,a3)')
     .   'iCv','x','y','z01'
        do iCv = 1, mpg%nCv
          write(lun,'(i5,1p,3e12.4)')
     .       iCv, gm%cvX(iCv), gm%cvY(iCv), work(iCv,0,iwrk)
        end do
        close(lun)
        stack_req=(/1,0,0,0,0,0,0/)
c-----------------------------------------------------------------------
      else if(ig.eq.-19) then
        call b2pbrooks(lun)
c-----------------------------------------------------------------------
      else if(ig.eq.-20) then
        allocate(mz_hlp(0:st%ns-1),hlp_al(0:st%ns-1))
        do is=0,st%ns-1
          mz_hlp(is)=nint((zamin(is)+zamax(is))/2.0_R8)
          hlp_al(is)=is_codes(is)
          if (hlp_al(is)(2:2).eq.' ') hlp_al(is)(2:2) = '_'
        end do
        open(lun,file='diagno.dat')
        write (lun,'(a)') '# '//trim(label)
        write (lun,'(a)')
     ,  '# iCv      : cell indices in the grid;'
        write (lun,'(a)')
     ,  '# R , Z    : R and Z coordinates of the cell centres [m];'
        write (lun,'(a)')
     ,  '# ne       : electron density [m^-3];'
        write (lun,'(a)')
     ,  '# Te, Ti   : electron and ion temperatures [eV];'
        write (lun,'(a)')
     ,  '# n...     : densities of the ion species [m^-3];'
        write (lun,*)
        write (hlp_frmr,'(a,i3,a)')
     ,   '(a6,5x,a1,11x,a1,10x,a2,10x,a2,10x,a2,4x,',
     ,     st%ns,'(4x,a3,i2.2,3x))'
        write (lun,hlp_frmr)
     ,    '#  iCv','R','Z','ne','Te','Ti',
     ,     ('n'//hlp_al(i),mz_hlp(i),i=0,st%ns-1)
        write (hlp_frmr,'(a,i3,a)') '(1x,1p,i6,',st%ns+5,'e12.4)'
        do iCv = 1, mpg%nCi
          write (lun,hlp_frmr) iCv,gm%cvX(iCv),gm%cvY(iCv),
     ,                         st%dv%ne(iCv),
     ,                         st%pl%te(iCv)/ev,st%pl%ti(iCv)/ev,
     ,                        (st%pl%na(iCv,is),is=0,st%ns-1)
        end do
        close(lun)
c-----------------------------------------------------------------------
      else if(ig.ge.1000) then
        call b2paspwr(ig-1000)
      else if(ig.gt.DEF_NYD) then
        call b2phel(lun)
      else if(ig.ge.0) then
        call b2poleg(lun,ig)
      else
        print *,'*** b2paux: parameter must be > -20 ***'
      end if
c=======================================================================
      return
      end subroutine b2paux
