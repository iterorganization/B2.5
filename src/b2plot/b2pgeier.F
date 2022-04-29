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
#ifdef B25_EIRENE
      use eirmod_wneutrals
#endif
      implicit none
      integer lun

      integer i,j,k,l,is,iy,jj
      integer iatm,ind,ias,iss
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
      do j=1,2*max(1,nncut) !{ Do not use ntrgts here to allow for sign flips
        if(itrg(j).gt.0) then !{
          k=ltrgpos(itrg(j))
          do iy=lytg,ny-1 !{
            do is=0,ns-1 !{
              hlp_flx(iy,is,itrg(j))=(fna(k+1-l,iy,0,0,is)+
     .                                fna(k+1-l,iy,0,1,is))*u
              if (is_neutral(is)) then !{
                iatm=b2eatcr(is)
                do iss=1,nsts !{
                  if (eirdiag_nds_typ(iss).ne.2) cycle
                  if (eirdiag_nds_srf(iss).ne.ltrgpos(itrg(j))+1) cycle
                  if (eirdiag_nds_start(iss).gt.iy+1) cycle
                  if (eirdiag_nds_end(iss).lt.iy+1) cycle
                  ind=eirdiag_nds_ind(iss)
                  ias=ind+1-(eirdiag_nds_start(iss)-1)
                  hlp_flx(iy,is,itrg(j))=hlp_flx(iy,is,itrg(j))+
     .             wldna_res(iatm,iy+ias)
                end do !}
              end if !}
            end do !}
          end do !}
        end if !}
        l=1-l
        u=-u
      end do !}

      do j=1,ntrgts !{
        jj=itrg(j)
        if(j.eq.ntrgts.and..not.(ntrgts.eq.1.and.itrg(1).eq.1))
     >   jj=ntrgsx
        if(.not.LSN) jj=ntrgsx+1-j
        open(lun,file=ftrg(jj))
        write(fmt_peak,'(a,i3,a)') '(a,5x,',ns,'(a2,i2.2,8x))'
c        print *,trim(fmt_peak)  !###
        write(lun,740) 'Ion fluxes on the target - ',
     ,      trim(trg_label(jj)),trim(label)
        write(lun,fmt_peak) '# cell',
     ,    ((trim(species_list(is)),i,i=0,nfluids(is)),
     ,                                               is=1,nspecies)
        write(fmt_peak,'(a,i3,a)') '(2x,i4,1p,',ns,'(e12.4))'
c        print *,trim(fmt_peak)  !###
        do iy=lytg,ny-1 !{
          write(lun,fmt_peak) iy,(hlp_flx(iy,is,itrg(j)),is=0,ns-1)
        end do !}
        close (lun)
      end do !}
      deallocate(hlp_flx)
      return
c=======================================================================
  740 format('# ',2a,' : ',a)
c=======================================================================
      end subroutine b2pgeier
