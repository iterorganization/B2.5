*-----------------------------------------------------------------------
*.specification

      subroutine b2stbr_bas (nx, ny, ns, nscx, iscx, dtim, BoRiS,
     &  hx, gs,
     &  na, ua, te, ti, po, ne, ni, fna, fhe, fhi, fht, fch, fne, fni,
     &  sna0, smo0, she0, shi0, sch0, sne0, tchem, tchee, tphys, tphye,
     &  thevp, thvpe, trese, tresn, trfln, trfle,
     &  coreregno, coreregn2, output,
     &  neutral_sources_rescale, core_sources_rescale,
     &  sput_src, sput_frc, sput_phys,
     &  sput_chem_alpha, sput_phys_alpha, alpha,
     &  sput_chem_model, therm_evap, sput_res,
     &  reflection_on, sputter_energy_on, main_call,
     &  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     &  shi0_ff, f_redep, potential_at_guard_cell,
     &  recycled_neutrals_contr)

      use b2mod_types
      use b2mod_rates
      use b2mod_sputter
      use b2mod_tallies
      use b2mod_elements
      use b2mod_indirect
      use b2mod_constants
      use b2mod_boundary_namelist
      use b2mod_wall
      use b2mod_b2cmpa
      use b2mod_b2cmpb
      use b2mod_b2cmfs
      use b2mod_external
      use b2mod_geo
     , , only : hy1
      use b2mod_math
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny, ns, nscx, iscx(0:nscxmax-1)
      integer coreregno, coreregn2, sput_src,
     *   sput_chem_model, reflection_on,
     *   sputter_energy_on, output, potential_at_guard_cell
      real (kind=R8) :: dtim, BoRiS,
     *  hx(-1:nx,-1:ny),
     &  gs(-1:nx,-1:ny,0:2)
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  ne(-1:nx,-1:ny), ni(-1:nx,-1:ny,0:1)
      real (kind=R8) ::
     *  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1),
     &  fne(-1:nx,-1:ny,0:1,0:1), fni(-1:nx,-1:ny,0:1,0:1),
     *  fht(-1:nx,-1:ny,0:1,0:1), f_redep(-1:nx,-1:ny,0:ns-1)
      real (kind=R8) ::
     *  neutral_sources_rescale, core_sources_rescale,
     &  sput_frc, sput_phys, sput_chem_alpha, sput_phys_alpha,
     &  therm_evap, sput_res, alpha,
     &  sput_chem_cutoff_alpha, sput_chem_cutoff_beta,
     &  recycled_neutrals_contr
      logical main_call
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     *  sna0(-1:nx,-1:ny,0:1,0:ns-1), smo0(-1:nx,-1:ny,0:3,0:ns-1),
     &  she0(-1:nx,-1:ny,0:3), shi0(-1:nx,-1:ny,0:3),
     &  sch0(-1:nx,-1:ny,0:3), sne0(-1:nx,-1:ny,0:1),
     &  shi0_ff(-1:nx,-1:ny,0:nscx-1),
     &  tchem(0:ns-1), tchee(0:ns-1), tphys(0:ns-1), tphye(0:ns-1),
     &  thevp(0:ns-1), thvpe(0:ns-1), trese(0:ns-1), tresn(0:ns-1),
     &  trfln(0:ns-1), trfle(0:ns-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2STBR_BAS provides linearized expressions for the sources of
*     density, parallel momentum and heat, due to special boundary
*     conditions and the neutral gas recycling model, using the default
*     Bas computational geometry.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, i, ix, iy, is, ireg, ncall, ibegin, iend, isize, k
      integer, allocatable :: list_x(:), list_y(:), rix(:), riy(:)
      integer, save :: iout = 0                          !srv 09.07.04
      real (kind=R8), allocatable :: area(:)
      character charns*3, chark*1
      character*15 hlp_frm
      real (kind=R8) ::
     *   t0, t1, t2, t3, qpl, qex, surface_temp,
     *   Ei, Pflx, vbar, yldchm, yldres, GammaEvap, Phi_app,
     *   redep_max, yself_max
      logical lcore
*   ..procedures
      intrinsic abs, max, exp
      external xertst
      external b2xvsg, b2xvff
*   ..initialisation
      save ncall
      data ncall/0/
*-----------------------------------------------------------------------
*.documentation-internal
*
*     Note: This version of b2stbr is for the case when recycling is
*     calculated via the fluid model.
*
*     cbrec - (0:7,0:nsdecl-1,0:cbnrmx-1) real*8 array.
*     For 0.le.is.lt.ns-1 and 0.le.ireg.lt.nreg, the eight coefficients
*     cbrec(0:7,is,ireg) specify the recycling process for species (is)
*     on the boundary region ireg.
*     cbrec(0,is,ireg) gives the recycling coefficient for particles
*         of charged species (is) to neutrals of the same species
*     cbrec(1,is,ireg) gives the associated energy recycling coefficient
*     cbrec(3,is,ireg) gives an additional contribution to physical
*         sputtering from species (is)
*     cbrec(4,is,ireg) gives an additional contribution to chemical
*         sputtering from species (is)
*     cbrec(6,is,ireg) gives the recycling coefficient for particles
*         of charged species (is) to (is+1) (intended mainly for the
*         core boundary)
*
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2stbr_bas')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       write(*,*) 'First call of b2stbr_bas'
       call ipgeti ('b2stbr_bas_iout', iout)                         !srv 11.07.05
      endif
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..extensive tests on first few calls
      if(ncall.lt.3) then
*    ..test sign of hx
       call b2xvsg (n2, hx, 1, 'hx', '.gt.')
*    ..test sign of na, ni, ne, te, ti
       call b2xvsg (n2*ns, na, 1, 'na', '.gt.')
       call b2xvsg (n2*2, ni, 1, 'ni', '.gt.')
       call b2xvsg (n2, ne, 1, 'ne', '.gt.')
       call b2xvsg (n2, te, 1, 'te', '.gt.')
       call b2xvsg (n2, ti, 1, 'ti', '.gt.')
      endif
*   ..test edge values of fluxes
      do is = 0, ns-1
       call b2xvff (nx, ny, fna(-1,-1,0,0,is), 'fna')
      enddo
      call b2xvff (nx, ny, fhe, 'fhe')
      call b2xvff (nx, ny, fhi, 'fhi')
      call b2xvff (nx, ny, fch, 'fch')
      call b2xvff (nx, ny, fne, 'fne')
      call b2xvff (nx, ny, fni, 'fni')

      redep_max=0.0_R8
      yself_max=0.0_R8

* ..recycling model at south boundary
*   ..initialize ireg
      ireg = cbirso
*   ..loop over boundary cells
      do ix = -1, nx
        if(region(ix,-1,0).ne.0.and.   ! do not do dead cells or corners
     &     leftix(ix,-1).ne.-2.and.rightix(ix,-1).ne.nx+1) then
*    ..update ireg
    1     if (cbrbrk(ireg)*nx.lt.ix+0.5_R8) then
            call xertst (ireg.lt.cbirso+cbnrso-1,
     &       'error in computing region on south boundary')
            ireg = ireg+1
            goto 1
          endif
          if (ireg.eq.cbirso) then
            ibegin = 0
          else
            ibegin = max(0, nint(cbrbrk(ireg-1)*nx))
          endif
          iend = min(nx, nint(cbrbrk(ireg)*nx))-1
          Phi_app = get_phi_applied(boundary_namelist,ix,-1)
          lcore=((ireg.eq.coreregno .and. isymm.ne.0) .or.
     &           (ireg.eq.coreregn2 .and. nncut.ge.2))   ! xpb: no sputtering in core!
          if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &     .and. ix.eq.ibegin) then
            isize = iend - ibegin + 1
            allocate(list_x(isize),list_y(isize),area(isize))
            do iy = 1, isize
              list_x(iy) = ibegin + iy - 1
              area(iy) =
     .         gs(topix(ibegin+iy-1,-1),topiy(ibegin+iy-1,-1),1)
            end do
            list_y(1:isize) = -1
            call prehd2D(ns,isize,list_x,list_y,dtim,-1,
     &       topix(ibegin:iend,-1),topiy(ibegin:iend,-1),1,area,BoRiS)
            deallocate(list_x,list_y,area)
          endif

          call recycle(ix,-1,topix(ix,-1),topiy(ix,-1),
     &         -1,topix(ix,-1),topiy(ix,-1),1,
     &         hy1(ix,-1),hy1(topix(ix,-1),topiy(ix,-1)),
     &         Phi_app,gs(topix(ix,-1),topiy(ix,-1),1),lcore)
        endif  ! region.ne.0 and corners
      enddo

* ..recycling model at north boundary
*   ..initialize ireg
      ireg = cbirno
*   ..loop over boundary cells
      do ix = -1, nx
        if(region(ix,ny,0).ne.0.and.    ! do not do dead cells or corners
     &     leftix(ix,ny).ne.-2.and.rightix(ix,ny).ne.nx+1) then
*    ..update ireg
    2     if (cbrbrk(ireg)*nx.lt.ix+0.5_R8) then
           call xertst (ireg.lt.cbirno+cbnrno-1,
     &      'error in computing region on north boundary')
           ireg = ireg+1
           goto 2
          endif
          if (ireg.eq.cbirno) then
            ibegin = 0
          else
            ibegin = max(0, nint(cbrbrk(ireg-1)*nx))
          endif
          iend = min(nx, nint(cbrbrk(ireg)*nx))-1
          Phi_app = get_phi_applied(boundary_namelist,ix,ny)
          lcore=on_closed_surface(ix,ny)
          if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &     .and. ix.eq.ibegin) then
            isize = iend - ibegin + 1
            allocate(list_x(isize),list_y(isize))
            do iy = 1, isize
              list_x(iy) = ibegin + iy - 1
            end do
            list_y(1:isize) = ny
            call prehd2D(ns,isize,list_x,list_y,dtim,+1,
     &        list_x,list_y,1,
     &        gs(ibegin:iend,ny,1),BoRiS)
            deallocate(list_x,list_y)
          endif
          call recycle(ix,ny,bottomix(ix,ny),bottomiy(ix,ny),
     &         +1,ix,ny,1,
     &         hy1(ix,ny),hy1(bottomix(ix,ny),bottomiy(ix,ny)),
     &         Phi_app,gs(ix,ny,1),lcore)
        endif   ! region.ne.0 and corners
      enddo

* ..recycling model at west boundary
      do ix = -1,nx                 !for double null case             !srv 20.07.04
       if(leftix(ix,-1).eq.-2) then                                   !srv 20.07.04
*   ..initialize ireg
        ireg = cbirwe
*   ..loop over boundary cells
        do iy = -1, ny
          if(region(ix,iy,0).ne.0.and.  ! do not do dead cells or corners
     &       bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1) then
*    ..update ireg
    3       if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
              call xertst (ireg.lt.cbirwe+cbnrwe-1,
     &        'error in computing region on west boundary')
              ireg = ireg+1
              goto 3
            endif
            if (ireg.eq.cbirwe) then
              ibegin = 0
            else
              ibegin = max(0, nint(cbrbrk(ireg-1)*ny))
            endif
            iend = min(ny, nint(cbrbrk(ireg)*ny))-1
            lcore=.false.
            Phi_app = get_phi_applied(boundary_namelist,-1,iy)
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. iy.eq.ibegin) then
              isize = iend - ibegin + 1
              allocate(list_x(isize),list_y(isize),
     &                 rix(isize),riy(isize))
              allocate(area(isize))
              list_x(1:isize) = ix
              do i = 1, isize
                list_y(i) = ibegin + i - 1
                rix(i) = rightix(ix,ibegin + i - 1)
                riy(i) = rightiy(ix,ibegin + i - 1)
                area(i) = gs(rix(i),riy(i),0)
              end do
              call prehd2D(ns,isize,list_x,list_y,dtim,-1,
     &          rix,riy,0,area,BoRiS)
              deallocate(list_x,list_y,rix,riy,area)
            endif
            call recycle(ix,iy,rightix(ix,iy),rightiy(ix,iy),
     &         -1,rightix(ix,iy),rightiy(ix,iy),0,
     &         hx(ix,iy),hx(rightix(ix,iy),rightiy(ix,iy)),
     &         Phi_app,gs(rightix(ix,iy),rightiy(ix,iy),0),lcore)
          endif   ! region.ne.0 and corners
        enddo
       endif
      enddo

* ..recycling model at east boundary
      do ix = -1,nx                 !for double null case             !srv 20.07.04
       if(rightix(ix,-1).eq.nx+1) then                                !srv 20.07.04
*   ..initialize ireg
        ireg = cbirea
*   ..loop over boundary cells
        do iy = -1, ny
          if(region(ix,iy,0).ne.0.and.   ! do not do dead cells or corners
     &       bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1) then
*    ..update ireg
    4       if (cbrbrk(ireg)*ny.lt.iy+0.5_R8) then
              call xertst (ireg.lt.cbirea+cbnrea-1,
     &         'error in computing region on east boundary')
              ireg = ireg+1
              goto 4
            endif
            if (ireg.eq.cbirea) then
              ibegin = 0
            else
              ibegin = max(0, nint(cbrbrk(ireg-1)*ny))
            endif
            iend = min(ny, nint(cbrbrk(ireg)*ny))-1
            Phi_app = get_phi_applied(boundary_namelist,nx,iy)
            lcore=.false.
            if(plate_model.eq.2 .and. .not.lcore .and. main_call
     &       .and. iy.eq.ibegin) then
              isize = iend - ibegin + 1
              allocate(list_x(isize),list_y(isize),area(isize))
              list_x(1:isize) = ix
              do i = 1, isize
                list_y(i) = ibegin + i - 1
                area(i) = gs(ix,ibegin + i - 1,0)
              end do
              call prehd2D(ns,isize,list_x,list_y,dtim,+1,
     &          list_x,list_y,0,area,BoRiS)
              deallocate(list_x,list_y,area)
            endif
            call recycle(ix,iy,leftix(ix,iy),leftiy(ix,iy),
     &         +1,ix,iy,0,
     &         hx(ix,iy),hx(leftix(ix,iy),leftiy(ix,iy)),
     &         Phi_app,gs(ix,iy,0),lcore)
         endif   ! region.ne.0 and corners
        enddo
       endif
      enddo

      if (iout.eq.1) then                                          !srv 11.07.05 {
        do is = 0, ns-1
          write (charns,'(i3.3)') is
          do k=0,3
            write (chark,'(i1)') k
            call my_out(70,nx,ny,smo0(-1,-1,k,is),
     &                       'smo_br_bas'//chark//charns)
          enddo
          do k=0,1
            write (chark,'(i1)') k
            call my_out(70,nx,ny,sna0(-1,-1,k,is),
     &                       'sna_br_bas'//chark//charns)
          enddo
        enddo
        do k=0,3
          write (chark,'(i1)') k
          call my_out(70,nx,ny,she0(-1,-1,k),'she_br_bas'//chark)
          call my_out(70,nx,ny,shi0(-1,-1,k),'shi_br_bas'//chark)
          call my_out(70,nx,ny,sch0(-1,-1,k),'sch_br_bas'//chark)
        enddo
        open(70,file='t_files.dat')
        write (hlp_frm,'(a4,i4,a7)') '(a7,',ns,'es13.5)'            !srv 11.03.13 {
        write(70,hlp_frm) ' tchem=',tchem
        write(70,hlp_frm) ' tchee=',tchee
        write(70,hlp_frm) ' tphys=',tphys
        write(70,hlp_frm) ' tphye=',tphye
        write(70,hlp_frm) ' thevp=',thevp
        write(70,hlp_frm) ' thvpe=',thvpe
        write(70,hlp_frm) ' trese=',trese
        write(70,hlp_frm) ' tresn=',tresn
        write(70,hlp_frm) ' trfln=',trfln
        write(70,hlp_frm) ' trfle=',trfle                           !srv 11.03.13 }
        close(70)
      endif                                                         !srv 11.07.05 }

      if (redep_max.gt.0.0_R8 .and. (main_call.or.output.gt.1))
     .  write(*,*) 'Maximum redeposition coefficient ', redep_max
      if(redep_max.gt.0.0_R8 .and. (yself_max.gt.1.0_R8 .or.
     .  ((main_call.and.output.ge.1).or.output.gt.1)))
     .   write(*,*) 'Maximum self-sputtering coefficient ', yself_max
      write(*,*) 'Maximum redeposition coefficient ', redep_max
      write(*,*) 'Maximum self-sputtering coefficient ', yself_max
* ..return
      ncall = ncall + 1
      call subend ()
      return

*-----------------------------------------------------------------------

      contains

      subroutine recycle(ix,iy,ixn,iyn,
     &     isign,ixf,iyf,idir,
     &     hh,hhn,Phi_app,area,
     &     lcore)
c ix,iy           position of the cell
c ixn,iyn         position of the neighbour cell
c isign           multiplies fluxes to get the right sign
c ixf,iyf         position where fluxes are evaluated
c idir            index for direction of flux (x or y)
c hh              either hx or hy1 @ ix,iy
c hhn             either hx or hy1 @ ixn,iyn
c Phi_app         applied potential
c area            area element
c lcore           is this a core element
      implicit none
      integer ix, iy, ixn, iyn, isign, ixf, iyf, idir
      real (kind=R8) hh, hhn, Phi_app, area, pof,
     *  redep_self_sputter, redep_self_energy,
     *  wrk(ndepth_nml)
#ifdef DBG
      real (kind=R8) product
#endif
      logical lcore, found

      integer is0, is1, is2, is, k, iwall, itrack
      real (kind=R8) :: tef, tif, taf, zaf(0:ns-1), ptf

      iwall=xymap(ix,iy)
      tef = (te(ix,iy)*hhn+te(ixn,iyn)*hh)/(hhn+hh)
      tif = (ti(ix,iy)*hhn+ti(ixn,iyn)*hh)/(hhn+hh)
      if (potential_at_guard_cell.eq.1) then
!        write(6,*) 'Potential taken at guard cell (bas)'
        pof = po(ix,iy)
      else
!        write(6,*) 'Potential taken at real cell (bas)'
        pof = po(ixn,iyn)
      endif
      zaf(0:ns-1) =
     &  (rza(ix,iy,0:ns-1)*hhn+rza(ixn,iyn,0:ns-1)*hh)/(hhn+hh)
      qpl = max(0.0_R8,isign*(fht(ixf,iyf,idir,0)+
     &                        fht(ixf,iyf,idir,1)))/area
      qex = max(0.0_R8,isign*(fhi_ext(ixf,iyf,idir,0)+
     &                        fhi_ext(ixf,iyf,idir,1)))/area
      do is=0,ns_ext-1
         ptf =
     &    (pt_ext(ix,iy,is)*hhn+pt_ext(ixn,iyn,is)*hh)/(hhn+hh)
         taf =
     &    (ta_ext(ix,iy,is)*hhn+ta_ext(ixn,iyn,is)*hh)/(hhn+hh)
         qex = qex + max(0.0_R8,isign*(fa_ext(ixf,iyf,idir,0,is)+
     &                                 fa_ext(ixf,iyf,idir,1,is)))
     &        *((0.5_R8*am_ext(is)*mp*ua_ext(ix,iy,is)**2+
     &           taf)*(1.0_R8-BoRiS)
     &          + ptf*ev)/area
      enddo
      if(lcore) then
         surface_temp = 0.0_R8
      else
         if(iwall.eq.0) then
            write(*,*) 'xymap returned zero for ix,iy = ',ix,iy
            call xerrab ('Invalid iwall setting!')
         endif
         if (plate_model.eq.0) then
            call heatdiff0D(iwall,qpl+qex,surface_temp)
            target_temp(iwall,1) = surface_temp
         else if (plate_model.eq.1) then
            if (main_call) then
              wrk(1:ndepth_nml) = target_temp(iwall,1:ndepth_nml)
              call heatdiff1D(iwall,ns,qpl+qex,
     .           dtim*plate_time_factor(iwall),wrk,output)
              target_temp(iwall,1:ndepth_nml) = wrk(1:ndepth_nml)
            endif
            surface_temp = target_temp(iwall,1)
         else if (plate_model.eq.2) then
           surface_temp = target_temp(iwall,1)
         endif
      endif

*    ..loop over species
      is0 = 0
      is1 = 0
      do is = 0, ns-1
         if (zn(is0).ne.zn(is).or.amtol.lt.abs(am(is0)-am(is))) then
            is0 = is
         endif
cdpc
         if(is_neutral(is)) is1=is
         if(zn(is1).ne.zn(is).or.amtol.lt.abs(am(is1)-am(is))) then
            write(*,*) 'b2stbr problem: is, is0, is1 ',is,is0,is1
         endif
         if(is0.ne.is1) then
            if(ncall.eq.0)
     1           write(*,*) 'b2stbr: setting is0 to is1 ',is0,is1
            is0=is1
         endif
         is2=-1
         k=0
         found=.false.
         do while (.not.found .and. k.lt.nscx)
           if (is0.eq.iscx(k)) then
             found = .true.
             is2 = k
           end if
           k = k + 1
        enddo
cdpc
*       (species (is) recycles into is0)
*       (is0 corresponds to CX species is2)
         t0 = max(0.0_R8,isign*(fna(ixf,iyf,idir,0,is)+
     &                          fna(ixf,iyf,idir,1,is))) *
     1        neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         Ei=max(0.0_R8,tif/ev+(pof-Phi_app)*zaf(is))
                                 ! xpb: Incident energy (eV)
         Pflx = t0/area          ! xpb: Incident particle flux (part/m^2/s)
         if (is_neutral(is)) then
            vbar = sqrt(8.0_R8 * ti(ix,iy) / (pi * am(is) * mp))
            Pflx = Pflx + na(ix,iy,is)*vbar*alpha *
     1           neutral_sources_rescale ! dpc: Eirene coupling 1999.06.28
         endif
         sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0)+cbrec(0,is,ireg)*t0
         itrack=track_index(is0)
         if(itrack.gt.0.and.itrack.le.ntrack
     1        .and.iwall.gt.0.and.main_call) then
            deposition(iwall,itrack)=deposition(iwall,itrack)+
     1           dtim*plate_time_factor(iwall)*
     2           (1.0_R8-cbrec(0,is,ireg))*t0
         endif
         shi0(ix,iy,0) = shi0(ix,iy,0)+recycled_neutrals_contr*max(
     &        cbrec(1,is,ireg)*(1.5_R8+BoRiS)*Ei*ev,
     &        (1.5_R8+BoRiS)*surface_temp*KBolt)*cbrec(0,is,ireg)*t0
         if (is2.ne.-1) shi0_ff(ix,iy,is2) = shi0_ff(ix,iy,is2)+
     &        recycled_neutrals_contr*max(
     &        cbrec(1,is,ireg)*(1.5_R8+BoRiS)*Ei*ev,
     &        (1.5_R8+BoRiS)*surface_temp*KBolt)*cbrec(0,is,ireg)*t0
         ! xpb: recycled neutrals come back with at least thermal energy
         if(.not.lcore) then     ! xpb: no sputtering in core!
c chemical sputtering
            if(sput_frc.gt.0.0_R8.and.nint(zn(is)).eq.sput_src.and.
     &         sput_dst.ge.0.and.sput_dst.lt.ns) then
               if(sput_chem_model.eq.1) then
                  if (Ei.gt.0.0_R8) then
#ifdef DBG
                     product =-((sput_chem_cutoff_alpha/Ei)
     &                         **sput_chem_cutoff_beta)
                     if(product.gt.cutll) then
                        t1=t0 * (sput_frc+cbrec(4,is,ireg)) *
     &                    exp(-(sput_chem_cutoff_alpha/Ei)**
     &                    sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
                        chemical_sputtering(iwall,is) =
     .                    (sput_frc+cbrec(4,is,ireg)) *
     &                    exp(-(sput_chem_cutoff_alpha/Ei)**
     &                    sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
                     else
                        t1=0.0_R8
                        chemical_sputtering(iwall,is) = 0.0_R8
                     end if
#else
                     t1=t0 * (sput_frc+cbrec(4,is,ireg)) *
     &                    exp(-(sput_chem_cutoff_alpha/Ei)**
     &                    sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
                     chemical_sputtering(iwall,is) =
     .                    (sput_frc+cbrec(4,is,ireg)) *
     &                    exp(-(sput_chem_cutoff_alpha/Ei)**
     &                    sput_chem_cutoff_beta) ! dpc: cut off chem sputter at 1eV
#endif
                  else
                     t1=0.0_R8
                     chemical_sputtering(iwall,is) = 0.0_R8
                  endif
               else
                  call chmsput(iwall,Ei,zaf(is),zn(is),
     &                 zn(sput_dst),am(is),am(sput_dst),surface_temp,
     &                 Pflx/neutral_sources_rescale,yldchm)
                  t1=Pflx*neutral_sources_rescale*
     &                 yldchm*(sput_frc+cbrec(4,is,ireg))*
     &                 area
                  chemical_sputtering(iwall,is) = yldchm
               endif
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               if(ncall .ne. 0 .or. sputter_namelist_used) then
                  t1=(sputter_yield(ix,iy,is,1)+sput_chem_alpha*t1)/
     *                 (1.0_R8+sput_chem_alpha)
               endif
               sputter_yield(ix,iy,is,1)=t1
               tchem(is)=tchem(is)+t1
               sna0(ix,iy,0,sput_dst) = sna0(ix,iy,0,sput_dst) + t1
               itrack=track_index(sput_dst)
               if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call)
     1              erosion(iwall,itrack) = erosion(iwall,itrack) +
     1              dtim*plate_time_factor(iwall)*t1
               if(sputter_energy_on.ne.0) then
                  if(ncall .ne. 0 .or. sputter_namelist_used) then
                     t2=(sputter_yield2(ix,iy,is,1)+sput_chem_alpha*t2)/
     *                    (1.0_R8+sput_chem_alpha)
                  endif
                  sputter_yield2(ix,iy,is,1)=t2
                  tchee(is)=tchee(is)+t2
                  shi0(ix,iy,0) = shi0(ix,iy,0)+t2
                  if (is2.ne.-1)
     >             shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               else
                  sputter_yield2(ix,iy,is,1)=0.0_R8
               endif
            else
               sputter_yield(ix,iy,is,1)=0.0_R8
               sputter_yield2(ix,iy,is,1)=0.0_R8
               chemical_sputtering(iwall,is) = 0.0_R8
            endif
c physical sputtering
            if(sput_phys.gt.0.0_R8.and.
     &         sput_dst.ge.0.and.sput_dst.lt.ns) then
               redep_max=max(redep_max,f_redep(ixn,iyn,sput_dst))
               redep_self_sputter =
     .              sputter(sput_dst+1,zaf(sput_dst+1),sput_dst,
     .                      (pof-Phi_app)*zaf(sput_dst+1)*ev)*
     1              (sput_phys+cbrec(3,sput_dst+1,ireg))
               if(f_redep(ixn,iyn,sput_dst).gt.0.0_R8 .and.
     1            redep_self_sputter.gt.0.0_R8 .and.
     2            (output.gt.2 .or. (main_call.and.output.ge.2) .or.
     3            redep_self_sputter.ge.1.0_R8))
     >           write(*,'(a,1x,a2,2i4,1p,2g12.5)')
     .            'Self-sputtering coefficient ', is_codes(sput_dst),
     .            ix, iy, pof-Phi_app, redep_self_sputter
               yself_max=max(yself_max,redep_self_sputter)
               t1=t0 * sputter(is,zaf(is),sput_dst,Ei*ev) *
     1           (sput_phys+cbrec(3,is,ireg)) *
     2           (1.0_R8-f_redep(ixn,iyn,sput_dst))/
     3           (1.0_R8-f_redep(ixn,iyn,sput_dst)*redep_self_sputter)
! the two lines above are the reduction of sputtering due to prompt redeposition
! and then the multiplication due to self-sputtering of the redeposited particles
               physical_sputtering(iwall,is) =
     .              sputter(is,zaf(is),sput_dst,Ei*ev) *
     1              (sput_phys+cbrec(3,is,ireg)) *
     2           (1.0_R8-f_redep(ixn,iyn,sput_dst))/
     3           (1.0_R8-f_redep(ixn,iyn,sput_dst)*redep_self_sputter)
               if(ncall .ne. 0 .or. sputter_namelist_used) then
                  t1=(sputter_yield(ix,iy,is,2)+sput_phys_alpha*t1)/
     *                 (1.0_R8+sput_phys_alpha)
               endif
               sputter_yield(ix,iy,is,2)=t1
               tphys(is)=tphys(is)+t1
               sna0(ix,iy,0,sput_dst) = sna0(ix,iy,0,sput_dst) + t1
               itrack=track_index(sput_dst)
               if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call) then
cxpb  We count here as eroded all the sputtered material, including the promptly redeposited fraction
cxpb  We add this promptly redeposited fraction in the deposited material
                 erosion(iwall,itrack) = erosion(iwall,itrack) +
     1            dtim*plate_time_factor(iwall)*t1/
     2            (1.0_R8-f_redep(ixn,iyn,sput_dst))*
     3            (1.0_R8-f_redep(ixn,iyn,sput_dst)*redep_self_sputter)
                 deposition(iwall,itrack) = deposition(iwall,itrack) +
     1            dtim*plate_time_factor(iwall)*t1*(1.0_R8-1.0_R8/
     2            (1.0_R8-f_redep(ixn,iyn,sput_dst))*
     3            (1.0_R8-f_redep(ixn,iyn,sput_dst)*redep_self_sputter))
               endif
               if(sputter_energy_on.ne.0) then
                  redep_self_energy =
     .              sputter2(sput_dst+1,sput_dst,
     .                       (pof-Phi_app)*zaf(sput_dst+1)*ev)*
     1              (sput_phys+cbrec(3,sput_dst+1,ireg))
                  t2=t0 * sputter2(is,sput_dst,Ei*ev) *
     1              (sput_phys+cbrec(3,is,ireg))*Ei*ev *
     2              (1.0_R8-f_redep(ixn,iyn,sput_dst))
cxpb  We now must remove the energy lost in ionizing the promptly redeposited material
cxpb  but from the electrons!
                  t3=-t0*ev*sputter(is,zaf(is),sput_dst,Ei*ev)*
     1              (sput_phys+cbrec(3,is,ireg))*
     1              rpi(ixn,iyn,sput_dst+1)*f_redep(ixn,iyn,sput_dst)/
     2              (1.0_R8-
     3                 f_redep(ixn,iyn,sput_dst)*redep_self_sputter)
cxpb  We must now add the energy brought back in by the self-sputtered material
                  t2=t2+t0*sputter(is,zaf(is),sput_dst,Ei*ev)*
     1              (sput_phys+cbrec(3,is,ireg))*
     .              f_redep(ixn,iyn,sput_dst)*
     1              redep_self_energy*(pof-Phi_app)*ev/
     2              (1.0_R8-
     3                 f_redep(ixn,iyn,sput_dst)*redep_self_sputter)
                  physical_sputtering_energy(iwall,is) =
     .                 sputter2(is,sput_dst,Ei*ev) *
     1                 (sput_phys+cbrec(3,is,ireg))
                  if(ncall .ne. 0 .or. sputter_namelist_used) then
                     t2=(sputter_yield2(ix,iy,is,2)+sput_phys_alpha*t2)/
     *                    (1+sput_phys_alpha)
                  endif
                  sputter_yield2(ix,iy,is,2)=t2
                  tphye(is)=tphye(is)+t2+t3
                  shi0(ix,iy,0) = shi0(ix,iy,0) + t2
                  she0(ix,iy,0) = she0(ix,iy,0) + t3
                  if (is2.ne.-1)
     >             shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               else
                  physical_sputtering_energy(iwall,is) = 0.0_R8
                  sputter_yield2(ix,iy,is,2)=0.0_R8
               endif
            else
               physical_sputtering(iwall,is) = 0.0_R8
               physical_sputtering_energy(iwall,is) = 0.0_R8
               sputter_yield(ix,iy,is,2)=0.0_R8
               sputter_yield2(ix,iy,is,2)=0.0_R8
            endif
c thermal evaporation of species is if bulk or tracked
            itrack=track_index(is)
            if(therm_evap.gt.0.0_R8.and.iwall.gt.0.and.
     .        (is.eq.sput_dst .or.
     1        (is_neutral(is).and.(itrack.gt.0.and.itrack.le.ntrack)))
     >        ) then
               call thermev(iwall, is, surface_temp, GammaEvap, t2)
               thermal_evaporation(iwall,is) = GammaEvap
               t1=GammaEvap*area*therm_evap
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
               if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call)
     1              erosion(iwall,itrack) = erosion(iwall,itrack) +
     1              dtim*plate_time_factor(iwall)*t1
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               thevp(is)=thevp(is)+t1
               thvpe(is)=thvpe(is)+t2
            else
               thermal_evaporation(iwall,is) = 0.0_R8
            endif
c radiation enhanced sublimation  ! the plate is assumed to be carbon !!!
            if(sput_res.gt.0.0_R8.and.
     &         sput_dst.ge.0.and.sput_dst.lt.ns) then
               call ressput(Ei,zaf(is),zn(is),zn(sput_dst),
     &              am(is),am(sput_dst),surface_temp,
     &              Pflx/neutral_sources_rescale,sput_dst,yldres)
               RES_sputtering(iwall,is) = yldres
               t1=yldres*neutral_sources_rescale*
     &              area*sput_res
               t2=t1*(1.5_R8+BoRiS)*surface_temp*KBolt
               sna0(ix,iy,0,sput_dst) = sna0(ix,iy,0,sput_dst) + t1
               itrack=track_index(sput_dst)
               if(itrack.gt.0.and.itrack.le.ntrack.and.
     1              iwall.gt.0.and.main_call)
     1              erosion(iwall,itrack) = erosion(iwall,itrack) +
     1              dtim*plate_time_factor(iwall)*t1
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               tresn(is)=tresn(is)+t1
               trese(is)=trese(is)+t2
            else
               RES_sputtering(iwall,is) = 0.0_R8
            endif
c back scattering
            if (cbrec(0,is,ireg).eq.0.0_R8 .and.
     &          cbrec(1,is,ireg).eq.0.0_R8 .and. .not.is_neutral(is)
     &    .and. reflection_on.ne.0) then ! use reflection model
               t1=t0 * reflect(is,Ei*ev,iwall)
               t2=t0 * reflect2(is,Ei*ev,iwall)*Ei*ev
               itrack=track_index(is0)
               if(itrack.gt.0.and.itrack.le.ntrack
     1              .and.iwall.gt.0.and.main_call) then
                  deposition(iwall,itrack)=deposition(iwall,itrack)-
     1                 dtim*plate_time_factor(iwall)*t1
               endif
               backscattering(iwall,is) = reflect(is,Ei*ev,iwall)
               backscattering_energy(iwall,is) =
     &              reflect2(is,Ei*ev,iwall)
               if(reflection_on.eq.1) then
                 sna0(ix,iy,0,is) = sna0(ix,iy,0,is) + t1
               else
                 sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0) + t1    ! reflection to the neutral
               endif
               shi0(ix,iy,0) = shi0(ix,iy,0) + t2
               if (is2.ne.-1)
     >           shi0_ff(ix,iy,is2)=shi0_ff(ix,iy,is2)+t2
               trfln(is)=trfln(is)+t1
               trfle(is)=trfle(is)+t2
            else
               backscattering(iwall,is) = 0.0_R8
               backscattering_energy(iwall,is) = 0.0_R8
            endif
         endif
*    ..recycling into next ionized species
         if(is.lt.ns-1) then
           if(lnext(is,is+1))
     1        sna0(ix,iy,0,is+1) = sna0(ix,iy,0,is+1)+
     2        cbrec(6,is,ireg)*t0/neutral_sources_rescale
         endif
      enddo
*    ..recycling at core into highest ionized stage
      if (lcore) then
*    ..loop over species
         is0 = ns-1
         is1 = ns-1
         do is = ns-1, 0, -1
            if (zn(is0).ne.zn(is).or.
     &           amtol.lt.abs(am(is0)-am(is))) then
               is0 = is
            endif
cxpb
            if(zamax(is).eq.zn(is)) is1=is
            if(zn(is1).ne.zn(is).or.amtol.lt.abs(am(is1)-am(is))) then
               write(*,*) 'b2stbr problem: is, is0, is1 ',is,is0,is1
            endif
            if(is0.ne.is1) then
               if(ncall.eq.0)
     1              write(*,*) 'b2stbr: setting is0 to is1 ',is0,is1
               is0=is1
            endif
cxpb
*          (species (is) recycles into is0)
            t0 = max(0.0_R8,isign*(fna(ixf,iyf,idir,0,is)+
     &                             fna(ixf,iyf,idir,1,is))) *
     .           core_sources_rescale ! xpb: 2002.02.01
            if (is.ne.is0) then
               sna0(ix,iy,0,is0) = sna0(ix,iy,0,is0)+
     &          cbrec(0,is,ireg)*t0
               shi0(ix,iy,0) = shi0(ix,iy,0)+recycled_neutrals_contr*
     &          cbrec(1,is,ireg)*(1.5_R8+BoRiS)*(zaf(is)*tef+tif)*t0
            endif
         enddo
      endif

      end subroutine recycle

*-----------------------------------------------------------------------
*.end b2stbr_bas

      end subroutine b2stbr_bas
