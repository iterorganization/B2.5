      subroutine b2pgeier(lun)

c  version : 15.02.2002 00:38

c=======================================================================
c*** Calculate and write down the profiles of ion fluxes on the targets
c=======================================================================
      use b2us_data
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot_wall_loading
      use b2us_data
      implicit none
      integer lun

      integer i,j,is,iy,jj,ibf,iFc
      real (kind=R8) :: u
      real (kind=R8), allocatable :: hlp_flx(:,:,:)
      character ftrg(ntrgsx)*12,trg_label(ntrgsx)*9,fmt_peak*80
      data ftrg /
     ,       'fi_tg_i.dat','fi_tg_iu.dat','fi_tg_ou.dat','fi_tg_o.dat'/
      data trg_label /'inner','inner_top','outer_top','outer'/
c=======================================================================
c      print *,'Starting b2pgeier...'      !###

c*** Determine topology and calculate some indices

      call b2pwinp
c      print *,'After b2pwinp: ntrgts,lrdlsp = ',  !###
c     ,                        ntrgts,lrdlsp       !###
      allocate(hlp_flx(mpg%nFc,0:state%ns-1,ntrgsx))

      do j=1,ntrgts !{
        if (itrg(j).eq.0) cycle
        do iy=1,mpg%divFcP(itrg(j),2) !{
          do is=0,state%ns-1 !{
            iFc=mpg%divFc(mpg%divFcP(itrg(j),1)+iy-1)
            u=mpg%divFcOr(mpg%divFcP(itrg(j),1)+iy-1)
            hlp_flx(iy,is,j)=(state%dv%fna(iFc,0,is)+
     .                        state%dv%fna(iFc,1,is))*u
          end do !}
        end do !}
      end do !}

      do j=1,ntrgts !{
        if (itrg(j).eq.0) cycle
        jj=j
        if(j.eq.ntrgts) jj=ntrgsx
        open(lun,file=ftrg(jj))
        write(fmt_peak,'(a,i3,a)') '(a,5x,',state%ns,'(a2,i2.2,8x))'
c        print *,trim(fmt_peak)  !###
        write(lun,740) 'Ion fluxes on the target - ',
     ,      trim(trg_label(jj)),trim(label)
        write(lun,fmt_peak) '# face',
     ,    ((trim(species_list(is)),i,i=0,nfluids(is)),
     ,                                               is=1,nspecies)
        write(fmt_peak,'(a,i3,a)') '(2x,i4,1p,',state%ns,'(e12.4))'
c        print *,trim(fmt_peak)  !###
        do iy=1,mpg%divFcP(itrg(j),2) !{
          write(lun,fmt_peak) iy,(hlp_flx(iy,is,j),is=0,state%ns-1)
        end do !}
        close (lun)
      end do !}
      deallocate(hlp_flx)
      return
c=======================================================================
  740 format('# ',2a,' : ',a)
c=======================================================================
      end subroutine b2pgeier
