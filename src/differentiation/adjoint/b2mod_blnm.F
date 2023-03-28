      module b2mod_blnm
      use b2mod_types
      use b2mod_diag_diff
      use b2mod_time
      implicit none

c======================================================================
c*** Computes the parallel momentum balance (forces) in the plasma in
c*** the inner and outer divertors (only lower divertors for
c*** double-null case)
c----------------------------------------------------------------------

c
c*** lv    :  version number to appear in the file
c*** lvm   :  version number corresponding to the last major change
c
      integer lv, lvm
      parameter (lv=970911)
      parameter (lvm=970909)

c
c----------------------------------------------------------------------
c*** Parallel momentum balance in the plasma:
c***
c***      header  var       Legend
c***                                                        (divertors)
c***              -- parallel forces on the plasma in the divertor (N):
c***     prs_xpnt prsup   pressure force upstream
c***     prs_trgt prstg   pressure force at the targets
c***     mfl_xpnt flmup   parallel momentum flux integral upstream
c***     mfl_trgt flmtg   parallel momentum flux integral at the targets
c***     fric_neu prsnu   friction with neutrals
c***     prs_magn prsmf   pressure taken by magnetic field
c***     res_frce rsfrc   resulting force onto divertor leg
c***                            -- normalised to the upstream pressure:
c***     rel_prtg rlptg   pressure force at the targets (relative)
c***     rel_flxp rlmup   parallel momentum flux integral upstream (relative)
c***     rel_fltg rlmtg   parallel momentum flux integral at the targets (rel.)
c***     rel_fneu rlfnu   friction with neutrals (relative)
c***     rel_fmag rlfmf   pressure taken by magnetic field (relative)
c***     rel_resf rlres   resulting force onto divertor leg (relative)
c----------------------------------------------------------------------

      integer lun, ltrc
      parameter (ltrc=ntrc_blm)
c*** Values by divertors. Here nqudiv is the number of values,
c*** ndivx is the maximum, and ndiv the actual, number
c*** of divertor legs (at present, the same)
      integer ndivx, ndiv, nqudiv, nd_tmp, ndv_tmp
      parameter (nqudiv=13, ndivx=2)

      character hlp_frm*80, hlp_frf*160
      character hlp_div(ndivx)*2, hlp_hdv(nqudiv)*8
c
      data hlp_div / 'il', 'ol' /
c----------------------------------------------------------------------
c*** Column headers:
      data hlp_hdv /
     1      'prs_xpnt','prs_trgt','mfl_xpnt','mfl_trgt','fric_neu',
     2      'prs_magn','res_frce','rel_prtg','rel_flxp','rel_fltg',
     3      'rel_fneu','rel_fmag','rel_resf'/
c----------------------------------------------------------------------

      save lun, ndiv, nd_tmp, ndv_tmp,
     & hlp_frm, hlp_frf, hlp_div, hlp_hdv

      contains

      subroutine b2blnm(nx,ny,ns)
      use b2mod_geo_diff
      use b2mod_rates
      use b2mod_plasma_diff
      use b2mod_sources_diff
      use b2mod_indirect
      use b2mod_constants
      use b2mod_external_diff
      use b2mod_b2cmpa_diff
      use b2mod_numerics_namelist_diff
#ifdef B25_EIRENE
      use eirmod_extrab25
#endif
      implicit none
c
c
      integer nx, ny, ns
c
      real (kind=R8) ::
     *     equdiv(ndivx,nqudiv)
      real (kind=R8) ::
     *     prsup(ndivx), prstg(ndivx), flmup(ndivx),
     ,     flmtg(ndivx), prsnu(ndivx), prsmf(ndivx),
     ,     rsfrc(ndivx), rlptg(ndivx), rlmup(ndivx),
     ,     rlmtg(ndivx), rlfnu(ndivx), rlfmf(ndivx),
     ,     rlres(ndivx)
      equivalence (equdiv(1,1), prsup), (equdiv(1,2), prstg),
     ,            (equdiv(1,3), flmup), (equdiv(1,4), flmtg),
     ,            (equdiv(1,5), prsnu), (equdiv(1,6), prsmf),
     ,            (equdiv(1,7), rsfrc), (equdiv(1,8), rlptg),
     ,            (equdiv(1,9), rlmup), (equdiv(1,10), rlmtg),
     ,            (equdiv(1,11), rlfnu), (equdiv(1,12), rlfmf),
     ,            (equdiv(1,13), rlres)
c----------------------------------------------------------------------
c
      integer i, j, ic, is, ix, iy
      real (kind=R8) ::
     *     hlp_pi, hlp_po, hlp_fi, hlp_fo, hlp_nu, hlp_mf,
     *     pr(-1:nx, -1:ny)
c----------------------------------------------------------------------
c*** To be taken from the common blocks (not used yet):
c***  lcbs: boundary segments representing the core boundary in B2
c
      save
     +  equdiv,
     +  prsup, prstg, flmup,
     ,  flmtg, prsnu, prsmf,
     ,  rsfrc, rlptg, rlmup,
     ,  rlmtg, rlfnu, rlfmf,
     ,  rlres,
     +  i, j, ic, is, ix, iy,
     +  hlp_pi, hlp_po, hlp_fi, hlp_fo, hlp_nu, hlp_mf
c======================================================================

      do i=1,ndivx
        do j=1,nqudiv
          equdiv(i,j)=0.0_R8
        end do
      end do
      if(nncut.ge.1) then !{
c*** We only need the pressure balance for divertors
c
c*** Pressure upstream
c
      hlp_pi=0.0_R8
      hlp_po=0.0_R8
      hlp_fi=0.0_R8
      hlp_fo=0.0_R8
      pr=ne*te
      do is = 0, ns-1
        if (is_neutral(is)) cycle
        pr(-1:nx,-1:ny)=pr(-1:nx,-1:ny)+
     .   na(-1:nx,-1:ny,is)*ti(-1:nx,-1:ny)
      end do
      do is = 0, ns_ext-1
        if (is_ext_neutral(is)) cycle
        pr(-1:nx,-1:ny)=pr(-1:nx,-1:ny)+
     .   na_ext(-1:nx,-1:ny,is)*ta_ext(-1:nx,-1:ny,is)
      end do
      do iy = 0, ny-1
        ix=leftcut(1)-1
        do is=0,ns-1
          if (.not.is_neutral(is)) then
            hlp_fi=hlp_fi+
     .       fmo(rightix(ix,iy),rightiy(ix,iy),0,0,is)
          end if
        end do
        hlp_pi=hlp_pi+pr(ix,iy)*sx(ix,iy)*pit(ix,iy)
        ix=rightcut(1)
        do is=0,ns-1
          if (.not.is_neutral(is)) then
            hlp_fo=hlp_fo+fmo(ix,iy,0,0,is)
           end if
        end do
        hlp_po=hlp_po+
     .   pr(ix,iy)*sx(leftix(ix,iy),leftiy(ix,iy))*pit(ix,iy)
      end do
      prsup(1)=hlp_pi
      prsup(2)=hlp_po
      flmup(1)=hlp_fi
      flmup(2)=hlp_fo
      hlp_pi=0.0_R8
      hlp_po=0.0_R8
      hlp_fi=0.0_R8
      hlp_fo=0.0_R8
      do iy = 0, ny-1
        hlp_pi=hlp_pi+pr(0,iy)*sx(0,iy)*pit(0,iy)
        hlp_po=hlp_po+pr(nx-1,iy)*sx(nx-1,iy)*pit(nx-1,iy)
        do is=0,ns-1
          if (.not.is_neutral(is)) then
            hlp_fi=hlp_fi+fmo(0,iy,0,0,is)
            hlp_fo=hlp_fo+fmo(nx,iy,0,0,is)
          end if
       end do
      end do
      prstg(1)=hlp_pi
      prstg(2)=hlp_po
      flmtg(1)=hlp_fi
      flmtg(2)=hlp_fo
c
c*** Volumetric sources
c
      hlp_nu=0.0_R8
      hlp_mf=0.0_R8
      do iy=0,ny-1
        do ix=0,leftcut(1)-1
          if(resignore(ix,iy,1).ne.0) then
            hlp_mf=hlp_mf+pit(ix,iy)*pr(ix,iy)*vol(ix,iy)
          end if
          if(resignore(ix,iy,2).ne.0) then
            do is=0,ns-1
              if (.not.is_neutral(is))
     >         hlp_nu=hlp_nu+
     1               smo(ix,iy,0,is)+
     2               smo(ix,iy,1,is)*ua(ix,iy,is)+
     3               smo(ix,iy,2,is)*am(is)*mp*na(ix,iy,is)+
     4               smo(ix,iy,3,is)*am(is)*mp*na(ix,iy,is)*
     5               ua(ix,iy,is)
            end do
          end if
        end do
      end do
      prsnu(1)=hlp_nu
      prsmf(1)=hlp_mf
      hlp_nu=0.0_R8
      hlp_mf=0.0_R8
      do iy=0,ny-1
        do ix=rightcut(1),nx-1
          if(resignore(ix,iy,1).ne.0) then
            hlp_mf=hlp_mf+pit(ix,iy)*pr(ix,iy)*vol(ix,iy)
          end if
        end do
        do ix=rightcut(1),nx-1
          if(resignore(ix,iy,2).ne.0) then
            do is=0,ns-1
              if (.not.is_neutral(is))
     >         hlp_nu=hlp_nu+
     1               smo(ix,iy,0,is)+
     2               smo(ix,iy,1,is)*ua(ix,iy,is)+
     3               smo(ix,iy,2,is)*am(is)*mp*na(ix,iy,is)+
     4               smo(ix,iy,3,is)*am(is)*mp*na(ix,iy,is)*
     5               ua(ix,iy,is)
            end do
          end if
        end do
      end do
      prsnu(2)=hlp_nu
      prsmf(2)=hlp_mf
c
c*** Totals and relative quantities
c
      rsfrc(1)=prstg(1)+flmtg(1)+prsnu(1)+prsmf(1)-prsup(1)-flmup(1)
      rsfrc(2)=prsup(2)+flmup(2)+prsnu(2)+prsmf(2)-prstg(2)-flmtg(2)
      do ic=1,2
        rlptg(ic)=prstg(ic)/prsup(ic)
        rlmup(ic)=flmup(ic)/prsup(ic)
        rlmtg(ic)=flmtg(ic)/prsup(ic)
        rlfnu(ic)=prsnu(ic)/prsup(ic)
        rlfmf(ic)=prsmf(ic)/prsup(ic)
        rlres(ic)=rsfrc(ic)/prsup(ic)
      end do
      end if !}
      if(ank_tracing.eq.0) return
      if(filedata(ltrc)) then
c
c*** Write the data in a separate file for plotting...
c
        write(lun,hlp_frm) tim,((equdiv(i,j),j=1,nd_tmp),i=1,ndv_tmp)
      else if(nncut.ge.1) then !{
c
c*** ... or just print out the data - if not a limiter configuration
c
        write(*,*) 'Parallel momentum balance in divertors:'
        write(*,hlp_frf) 'div',(hlp_hdv(i),i=1,nqudiv)
        write(*,hlp_frm) (hlp_div(i),(equdiv(i,j),j=1,nqudiv),i=1,ndiv)
      end if !}
      return

      contains
c----------------------------------------------------------------------
      function sx(ix,iy)
      implicit none
      integer, intent(in) :: ix, iy
      real (kind=R8) ::
     *     sx
csw 29sep2014: xpb fix on sx calculation
      sx=0.5_R8*(gs(ix,iy,0)*qc(ix,iy,0)+
     .  gs(rightix(ix,iy),rightiy(ix,iy),0)*
     .  qc(rightix(ix,iy),rightiy(ix,iy),0))
      return
      end function sx

      end subroutine b2blnm
c======================================================================
c
      subroutine b2blnmi(switch)
      use b2mod_switches_diff
      implicit none
      type (switches), intent (in) :: switch
      integer i, j, nv
      character f*(lenfnm), hlp_c
      logical :: ex, warn_tmp
c----------------------------------------------------------------------
c*** Initialise the file for tracing data
c*** Only called with filedata(ntrc_blm)=.true.
c----------------------------------------------------------------------
c
      ndiv=ndivx
c
c*** lun   :  logical unit number for the output file (=ntrc(ltrc))
c*** f     :  the associated file name
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
c
      call find_file(f,ex)
      open(lun,file=f,err=90)
      rewind(lun)
      if(ex) read(lun,*,end=91) hlp_c
      goto 93
   91 ex = .false.
   93 continue
      rewind(lun)
      if(ex.and.switch%b2mndr_stim.lt.0.0_R8) then
        read(lun,*,err=92) hlp_c,nv,nd_tmp,ndv_tmp
  700   format('#   ',3i8/'#  ',a)
        lv_trc(ltrc)=nv
        warn_tmp=.false.
        if(nv.ne.lv) then
          warn_tmp=.true.
          write(*,'(a,i8,3a,i8)')
     ,     '!!! b2blnm: current version ',lv,
     ,     ' differs from that in the file ',trim(f),' : ',nv
          if(nv .lt. lvm) then
c*** It was a major change of the file format - impossible to continue
            write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
            filedata(ltrc)=.false.
            go to 98
          end if
        end if
        if(ndv_tmp.ne.ndiv) then
          warn_tmp=.true.
          write(*,'(a,i2,3a,i2)')
     ,     '!!! b2blnm: current number of divertor legs ',ndiv,
     ,     ' differs from that in the file ',trim(f),' : ',ndv_tmp
          filedata(ltrc)=.false.
        end if
        if(nd_tmp.ne.nqudiv) then
          warn_tmp=.true.
          write(*,'(2a,i4,3a,i4)')
     ,     '!!! b2blnm: current number of columns',
     ,     ' for divertor data ',nqudiv,
     ,     ' differs from that in the file ',trim(f),' : ',nd_tmp
          if(nd_tmp.gt.nqudiv) then
            filedata(ltrc)=.false.
          end if
        end if
c
        if(.not.filedata(ltrc)) then
          write(*,*) 'No ',trim(f),' tracing file - ',
     ,                                'remove the existing file first!'
          go to 98
        end if
        if(warn_tmp) write(*,*) ' You should remove the file!'
c*** Read to the end of the file
   50   read(lun,'(a1)',end=55) hlp_c
        go to 50
   55   backspace lun
      else
c*** New file
        nd_tmp=nqudiv
        ndv_tmp=ndiv
        write(lun,700) lv,nd_tmp,ndv_tmp,label
        write(hlp_frf,'(a,i3,a)')
     ,   '(a6,12x,',nd_tmp*ndv_tmp,'(2x,a11,1x))'
        write(lun,hlp_frf) '# time',
     -   ((hlp_hdv(i)//'_'//hlp_div(j),i=1,nd_tmp),j=1,ndv_tmp)
      end if
      write(hlp_frm,'(a,i3,a)')
     , '(1p,e16.7,',nd_tmp*ndv_tmp,'(2x,e12.3))'
      return
c======================================================================
   90 write(*,*) 'b2blnm: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      go to 98
c
   92 write(*,*) 'b2blnm: wrong format of the file ',trim(f)
      filedata(ltrc)=.false.
      go to 98
c
c*** Prepare the formats for the printout
   98 write(hlp_frf,710) nqudiv
      write(hlp_frm,720) nqudiv
      return
  710 format('(1x,a3,',i3,'(2x,a8,1x))')
  720 format('(1p,1x,a2,1x,',i3,'e11.3)')
c======================================================================
      end subroutine b2blnmi
c
      subroutine b2blnmc(iline)
      implicit none
      integer, intent(out) :: iline
      real (kind=R8) :: tdummy
      character f*(lenfnm)
      character hlp_c
      logical ex
c----------------------------------------------------------------------
c*** - to be called by b2co for locating truncation for savefile
c----------------------------------------------------------------------
c
      if(.not.filedata(ntrc_blm)) return
      lun=ntrc(ltrc)
      f=fn_ntrc(ltrc)
      call find_file(f,ex)
      open(lun,file=f,err=190)
      rewind(lun)
      iline = 0
      read(lun,'(a1)',err=155,end=155) hlp_c
      do while (hlp_c.eq.'#')
        iline = iline+1
        read(lun,'(a1)',err=155,end=155) hlp_c
      end do
      backspace lun
      read(lun,*,err=155,end=155) tdummy
      do while (tdummy.le.tim)
        iline = iline+1
        read(lun,*,err=155,end=155) tdummy
      end do
  155 continue
!      close(lun)
      return
c======================================================================
  190 write(*,*) 'b2blnm: cannot open file ',trim(f)
      filedata(ltrc)=.false.
      iline=-1
      return
      end subroutine b2blnmc
c======================================================================
c
      subroutine b2blnms
c----------------------------------------------------------------------
c*** - to be called BEFORE each timestep
c----------------------------------------------------------------------
c
      return
      end subroutine b2blnms
c======================================================================
c
      subroutine b2blnmf
c----------------------------------------------------------------------
c*** - Finalization step
c----------------------------------------------------------------------
c
      logical file_open

      lun=ntrc(ltrc)
      inquire(unit=lun,opened=file_open)
      if (file_open) close(lun)

      return
      end subroutine b2blnmf
c======================================================================
      end module b2mod_blnm
