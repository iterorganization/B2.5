      subroutine b2pgeier(lun)

c  version : 15.02.2002 00:38

c=======================================================================
c*** Calculate and write down the profiles of ion fluxes on the targets
c=======================================================================
      use b2mod_types
      use b2mod_diag
      use b2mod_b2cmpa
      use b2mod_b2plot
      use b2mod_plasma
      use b2mod_indirect
      use b2mod_b2plot_wall_loading
      implicit none
      integer lun

      integer i,j,k,l,is,iy,jj
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
c      print *,'After b2pwinp: ntrgts,lytg,nytg,lrdlsp = ',  !###
c     ,                        ntrgts,lytg,nytg,lrdlsp       !###
      allocate(hlp_flx(0:ny-1,0:ns-1,ntrgsx))

      l=1
      u=-1
      do j=1,ntrgts !{
        if (ntrgts.eq.1 .and. ltrgpos(1).eq.nx-1) then
          l=0
          u=1
        end if
        k=ltrgpos(j)
        do iy=lytg,ny-1 !{
          do is=0,ns-1 !{
            hlp_flx(iy,is,j)=fna(k+1-l,iy,0,is)*u
          end do !}
        end do !}
        l=1-l
        u=-u
      end do !}

      do j=1,ntrgts !{
        jj=j
        if(j.eq.ntrgts) jj=ntrgsx
        if(.not.LSN) jj=ntrgsx+1-j
        open(lun,file=ftrg(jj))
        write(fmt_peak,'(a,i3,a)') '(a,5x,',ns,'(a,i2.2,9x))'
c        print *,trim(fmt_peak)  !###
        write(lun,740) 'Ion fluxes on the target - ',
     ,      trim(trg_label(jj)),trim(label)
        write(lun,fmt_peak) '# cell',
     ,    ((trim(adjustl(species_list(is))),i,i=0,nfluids(is)),
     ,                                               is=1,nspecies)
        write(fmt_peak,'(a,i3,a)') '(2x,i4,1p,',ns,'(e12.4))'
c        print *,trim(fmt_peak)  !###
        if(ntrgts.eq.4) then !{
          k=j
        else !}{
          k=ntrgts
        end if !}
        do iy=lytg,ny-1 !{
          write(lun,fmt_peak) iy,(hlp_flx(iy,is,j),is=0,ns-1)
        end do !}
        close (lun)
      end do !}
      deallocate(hlp_flx)
      return
c=======================================================================
  740 format('# ',2a,' : ',a)
c=======================================================================
      end subroutine b2pgeier
